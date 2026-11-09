import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { sigmaInterfaces, SIGMA_GRIDS } from '../js/dynamics/sigmaCore.module.js';

const LEVELS = ['bl60', ...(SIGMA_GRIDS.includes('bl48') ? ['bl48'] : [])];
import { deckProfile } from './helpers/deckProfile.mjs';
import { decisionTracker, cpuDecisions, gpuDecisions, neighbourhood, radiationTraces } from './helpers/decisions.mjs';
import { LATENT_HEAT } from '../js/physics/moist.module.js';
import { PREVIOUS_RADIATION, PREVIOUS_MOIST, PREVIOUS_SURFACE, PREVIOUS_BOUNDARY_LAYER } from './helpers/previousDefaults.mjs';
import { F32_RMS } from './helpers/distribution.mjs';
import { groupLayout, bandLayers, spreadGroups, groupChecks } from './helpers/spreadGroups.mjs';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuCore } = gpuAvailable ? await import('../js/gpu/core.gpu.js') : {};
const N = 16;
const MOIST = { ...PREVIOUS_MOIST, cloudLifetime: 3 * 3600, plumeCape: 70 };

function meanTheta(model) {
  const { K } = model.core, C = model.mesh.nCells, theta = model.state[1];
  return Float64Array.from({ length: K }, (_, k) => { let s = 0, a = 0; for (let i = 0; i < C; i++) { s += model.mesh.areaCell[i] * theta[k * C + i]; a += model.mesh.areaCell[i]; } return s / a; });
}
function stats(cpu, gpu) {
  let maxDiff = 0, sumSq = 0, sumRef = 0, at = -1;
  for (let x = 0; x < cpu.length; x++) { const d = Math.abs(cpu[x] - gpu[x]); if (d > maxDiff || (Number.isNaN(d) && !Number.isNaN(maxDiff))) { maxDiff = d; at = x; } sumSq += d * d; sumRef += cpu[x] * cpu[x]; }
  return { maxDiff, at, rms: Math.sqrt(sumSq / cpu.length), rmsRel: Math.sqrt(sumSq / Math.max(sumRef, 1e-300)) };
}

function columnEnergy(after, C, K, held) {
  let worst = 0;
  for (let i = 0; i < C; i++) {
    let sum = 0;
    if (held) {
      for (let k = 0; k < K; k++) sum += after.RADSW[k * C + i];
      worst = Math.max(worst, Math.abs(sum - after.RADATM[i]));
    } else {
      for (let k = 0; k < K; k++) sum += after.LWH[k * C + i];
      worst = Math.max(worst, Math.abs(sum + after.OLR[i] + after.SFLUX[i] - after.ABS[i] + after.ATMSW[i] + after.SH[i] + LATENT_HEAT * after.EVAP[i]));
    }
  }
  return worst;
}

async function deckPair(steps, { levels = 'bl60', turbulence = 'moist', surface = {}, boundaryLayer = {}, ...options } = {}) {
  const physics = { ...PREVIOUS_RADIATION, mixedLayerDeck: true, deckRest: 'depth', minimumInversion: 2, ...options }, sigma = sigmaInterfaces(levels);
  const model = createModel(new Grid(N), { ocean: false, levels: sigma, radiation: physics, moist: MOIST, boundaryLayer: { turbulence, velocityScale: 'capped', ...boundaryLayer }, surface: { exchange: 'roughness', ...PREVIOUS_SURFACE, ...surface } });
  const { ambiguous } = deckProfile(model), traces = radiationTraces(model);
  const gpu = await createGpuCore(model.mesh, { levels: sigma, nu4: model.core.nu4, nu4Theta: model.core.nu4Theta, divergenceDamping: model.core.divergenceDamping, referenceTheta: meanTheta(model), physics: { ...physics, ...MOIST, ...surface, turbulence, velocityScale: 'capped', ...boundaryLayer, surfaceExchange: 'roughness', exchangeOptions: PREVIOUS_SURFACE } });
  gpu.upload(model.state);
  gpu.uploadPhysics({ mlmSubsidence: model.radiation.mlmSubsidence, mlmHeight: model.radiation.mlmHeight });
  const C = model.mesh.nCells, K = model.core.K, decisions = decisionTracker(K, C), runs = new Set();
  for (let n = 0; n < steps; n++) {
    const time = model.time;
    model.step(900);
    await gpu.stepModel(900, time);
    const now = await gpu.downloadPhysics();
    decisions.check(cpuDecisions(model, traces), gpuDecisions(await gpu.download(), now));
    for (let i = 0; i < C; i++) {
      const a = model.boundaryLayer.cloudTopCooling[i], b = now.CTCOOL[i];
      if (Math.abs(a - b) > 1e-3 * Math.max(Math.abs(a), Math.abs(b))) runs.add(i);
    }
  }
  const parted = new Set([...decisions.parted, ...runs]);
  const after = await gpu.downloadPhysics(), state = await gpu.download(), r = model.radiation, bl = model.boundaryLayer;
  const near = neighbourhood(model.mesh, parted, 2), kept = Array.from({ length: C }, (_, i) => i).filter((i) => !near.has(i));
  const cells = (cpu, name) => stats(Float64Array.from(kept, (i) => cpu[i]), Float64Array.from(kept, (i) => after[name][i]));
  const keptColumns = (cpu, gpu) => stats(Float64Array.from({ length: K * kept.length }, (_, x) => cpu[Math.floor(x / kept.length) * C + kept[x % kept.length]]), Float64Array.from({ length: K * kept.length }, (_, x) => gpu[Math.floor(x / kept.length) * C + kept[x % kept.length]]));
  const columns = (cpu, name) => keptColumns(cpu, after[name]);
  const meanGap = (cpu, gpu, layers = 1) => { let a = 0, b = 0; for (let k = 0; k < layers; k++) for (const i of kept) { a += Math.abs(cpu[k * C + i]); b += gpu[k * C + i] - cpu[k * C + i]; } return Math.abs(b) / a; };
  const nearMean = (cpu, gpu) => Math.max(...[true, false].map((north) => { let sum = 0, n = 0; for (const i of near) if ((model.mesh.latCell[i] >= 0) === north) { sum += Math.abs(gpu[i] - cpu[i]); n++; } return n ? sum / n : 0; }));
  const bounds = DECK_S9[`${levels} ${turbulence} ${JSON.stringify({ surface, boundaryLayer, ...options })}`], layout = bounds && groupLayout(model.mesh, { land: model.geography?.land ?? null });
  const fields = { q: [model.state[4], state[4], K], theta: [model.state[1], state[1], K], mixing: [bl.mixing, after.MIX, K], longwave: [r.longwave, after.LWH, K], entrainment: [r.mlmEntrainment, after.MLMENT, 1], blEntrainment: [bl.entrainment, after.ENTRAIN, 1], sw: [r.surfaceShortwave, after.SWDN, 1] };
  const s9 = Object.fromEntries(Object.entries(bounds ?? {}).map(([field, groups]) => { const [cpu, gpu, layers] = fields[field]; return [field, groupChecks(spreadGroups(cpu, gpu.subarray(0, layers * C), layout, bandLayers(model.core.sigmaMid, layers)), groups)]; }));
  let decked = 0, gpuDecked = 0, layerMiss = 0, homeMiss = 0, fromMiss = 0, sealed = 0, moved = 0;
  for (let i = 0; i < C; i++) {
    if (r.mlmTop[i] > 0) decked++;
    if (after.MLMTOP[i] > 0) gpuDecked++;
    if (r.mlmLayer[i] !== after.MLMK[i]) layerMiss++;
    if (r.mlmTop[i] > 0 && r.deckLayer[i] !== after.DECKK[i]) homeMiss++;
    if (r.mlmTop[i] > 0 && r.deckFrom[i] !== after.DECKFROM[i]) fromMiss++;
    if (r.mlmLayer[i] > 0 && bl.mixing[(r.mlmLayer[i] - 1) * C + i] === 0 && after.MIX[(r.mlmLayer[i] - 1) * C + i] === 0 && bl.entrainment[i] === r.mlmEntrainment[i]) sealed++;
    if (r.mlmTop[i] > 0 && r.deckLayer[i] !== r.stratusLayer) moved++;
  }
  return {
    C, K, ambiguous, decked, gpuDecked, layerMiss, homeMiss, fromMiss, sealed, moved, parted: parted.size, runs: runs.size, near: near.size,
    share: cells(r.mlmShare, 'MLMM'), top: cells(r.mlmTop, 'MLMTOP'), height: cells(r.mlmHeight, 'MLMH'), water: cells(r.mlmWater, 'MLMWATER'), entrainment: cells(r.mlmEntrainment, 'MLMENT'),
    deck: cells(r.stratus, 'DECK'), fraction: cells(r.stratusFraction, 'DECKF'), olr: cells(r.outgoing, 'OLR'), sw: cells(r.surfaceShortwave, 'SWDN'), blEntrainment: cells(bl.entrainment, 'ENTRAIN'),
    gpuSealed: Array.from({ length: C }, (_, i) => after.MLMK[i]).filter((k) => k > 0.5).length, gpuEnergy: columnEnergy(after, C, K, (physics.radiationEvery ?? 1) > 1),
    means: { q: meanGap(model.state[4], state[4], K), mixing: meanGap(bl.mixing, after.MIX, K), entrainment: meanGap(r.mlmEntrainment, after.MLMENT), blEntrainment: meanGap(bl.entrainment, after.ENTRAIN), sw: meanGap(r.surfaceShortwave, after.SWDN), longwave: turbulence === 'moist' ? meanGap(r.longwave, after.LWH, K) : null },
    kept, cpuWater: Float64Array.from(r.mlmWater), gpuWater: Float64Array.from(after.MLMWATER),
    mixing: columns(bl.mixing, 'MIX'), longwave: turbulence === 'moist' ? columns(r.longwave, 'LWH') : null, theta: keptColumns(model.state[1], state[1]), q: keptColumns(model.state[4], state[4]),
    s9,
    nearMeans: { share: nearMean(r.mlmShare, after.MLMM), top: nearMean(r.mlmTop, after.MLMTOP), height: nearMean(r.mlmHeight, after.MLMH), water: nearMean(r.mlmWater, after.MLMWATER), entrainment: nearMean(r.mlmEntrainment, after.MLMENT), deck: nearMean(r.stratus, after.DECK), fraction: nearMean(r.stratusFraction, after.DECKF), olr: nearMean(r.outgoing, after.OLR), sw: nearMean(r.surfaceShortwave, after.SWDN), blEntrainment: nearMean(bl.entrainment, after.ENTRAIN) },
  };
}

const show = (label, run) => console.log(`${label}: deck on ${run.decked} of ${run.C} cells (GPU ${run.gpuDecked}), ${run.parted} cells with moist decisions or the boundary layer's cloud-top run apart (${run.runs} the run; ${run.near} cells left out with their two rings); inversion layer differs on ${run.layerMiss}, deck layer on ${run.homeMiss}; ${run.sealed} sealed alike, ${run.moved} decks off σ 0.92; rms relative differences: share ${run.share.rmsRel.toExponential(1)} (max ${run.share.maxDiff.toExponential(1)}), top ${run.top.rmsRel.toExponential(1)}, height ${run.height.rmsRel.toExponential(1)}, water ${run.water.rmsRel.toExponential(1)}, w_e ${run.entrainment.rmsRel.toExponential(1)}, mixing ${run.mixing.rmsRel.toExponential(1)}, boundary-layer w_e ${run.blEntrainment.rmsRel.toExponential(1)}, θ ${run.theta.rmsRel.toExponential(1)}, q ${run.q.rmsRel.toExponential(1)}, OLR ${run.olr.rmsRel.toExponential(1)}, surface shortwave ${run.sw.rmsRel.toExponential(1)}, longwave heating ${run.longwave ? run.longwave.rmsRel.toExponential(1) : "not stored by the GPU under the dry boundary layer"}`);

const DECK_S9 = JSON.parse(readFileSync(new URL('./helpers/subgridDeckS9.json', import.meta.url), 'utf8'));
const SEAL_BACKSTOP = {
  bl60: { share: 3.35e-05, top: 0.0081, height: 0.00739, water: 3.91e-05, entrainment: 3.51e-07, deck: 1.19e-08, fraction: F32_RMS, olr: 0.000363, sw: 0.00945, blEntrainment: 3.52e-07 },
  bl48: { share: 3.2e-05, top: 0.0076, height: 0.0074, water: 1.91e-05, entrainment: 4.11e-07, deck: 1.2e-08, fraction: F32_RMS, olr: 0.00181, sw: 0.00272, blEntrainment: 4.11e-07 },
};

function agree(label, run, onMeans = [], backstop = null) {
  const pick = (name) => (run.s9[name] ? 0 : onMeans.includes(name) ? run.means[name] : run[name].rmsRel);
  assert.ok(run.decked > 0.5 * run.C && run.gpuDecked === run.decked, `${label}: deck on ${run.decked} cells, ${run.gpuDecked} on the GPU`);
  assert.ok(run.layerMiss === 0 && run.homeMiss === 0, `${label}: inversion layer apart on ${run.layerMiss}, deck layer on ${run.homeMiss}`);
  assert.ok(run.share.maxDiff < 1e-3 && run.top.rmsRel < 1e-5 && run.height.rmsRel < 1e-5, `${label}: share ${run.share.maxDiff}, top ${run.top.rmsRel}, height ${run.height.rmsRel}`);
  assert.ok(run.water.rmsRel < 1e-4 && pick('entrainment') < 1e-4 && run.deck.maxDiff < 5e-5 && run.fraction.maxDiff < 1e-3, `${label}: water ${run.water.rmsRel}, w_e ${pick('entrainment')}, deck ${run.deck.maxDiff}, fraction ${run.fraction.maxDiff}`);
  assert.ok(pick('mixing') < 1e-4 && pick('blEntrainment') < 1e-4, `${label}: mixing ${pick('mixing')}, boundary-layer w_e ${pick('blEntrainment')}`);
  assert.ok((run.s9.theta ? 0 : run.theta.rmsRel) < 1e-6 && pick('q') < 1e-5, `${label}: θ ${run.theta.rmsRel}, q ${pick('q')}`);
  assert.ok(run.olr.rmsRel < 3e-5 && pick('sw') < 3e-5 && (!run.longwave || pick('longwave') < 1e-4), `${label}: OLR ${run.olr.rmsRel}, surface shortwave ${pick('sw')}, longwave heating ${run.longwave && pick('longwave')}`);
  const groups = Object.entries(run.s9), held = Object.entries(backstop ?? {}).map(([field, bound]) => [field, run.nearMeans[field], bound]);
  if (groups.length) console.log(`${label}: each column's largest gap over its band's layers, by group: ${groups.map(([field, check]) => `${field} ${check.rows.join('; ')}`).join('; ')}`);
  if (held.length) console.log(`${label}: over the ${run.near} cells left out, the larger hemisphere's mean gap: ${held.map(([field, mean, bound]) => `${field} ${mean.toExponential(1)} of ${bound}`).join(', ')}`);
  assert.deepEqual([...groups.flatMap(([field, check]) => check.failed.map((text) => `${field} ${text}`)), ...held.filter(([, mean, bound]) => !(mean <= bound)).map(([field, mean, bound]) => `${field}: mean ${mean} over the cells left out, past ${bound}`)], [], `${label}: the engines' gaps past the bounds the CPU twins set`);
}

test('on bl60, and bl48 where it is built, the sub-grid inversion, the deck on it and the boundary layer it seals match between the engines, under the moist and the dry boundary layer, the dry one diagnosed after the deposit; at the default order humidity, θ, mixing, longwave heating, the deck\'s and the boundary layer\'s entrainment and the surface shortwave part the engines no further than CPU twins part themselves, in each band of reference pressure each column\'s largest gap in the largest, a ladder of counts and each hemisphere\'s mean, trimmed mean and trimmed signed mean, and over the cells left out beside the parted decisions each other field\'s mean gap stays within twice the twins\'', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  for (const [levels, turbulence] of [...LEVELS.map((levels) => [levels, 'moist']), ['bl60', 'dry']]) {
    const held = turbulence === 'dry', name = `${levels} subgrid, ${turbulence}${held ? ', with the boundary layer diagnosed after the deposit' : ''}`, run = await deckPair(3, { levels, deckInversion: 'subgrid', turbulence, boundaryLayer: held ? PREVIOUS_BOUNDARY_LAYER : {} });
    show(name, run);
    agree(name, run, [], held ? null : SEAL_BACKSTOP[levels]);
    assert.ok(run.sealed > 0.5 * run.decked, `${turbulence}: ${run.sealed} of ${run.decked} sealed alike`);
  }
});

test('on bl60, and bl48 where it is built, the deck\'s cloud at its own top, on the grid ceiling and on the sub-grid inversion, with the radiation every step and held, matches between the engines, bl48\'s on the sub-grid inversion with the boundary layer diagnosed after the deposit; outside that, humidity, θ, mixing, longwave heating, entrainment and surface shortwave part the engines no further than CPU twins part themselves, group by group', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  for (const options of [{ deckCloud: 'top' }, ...LEVELS.map((levels) => ({ levels, deckCloud: 'top', deckInversion: 'subgrid' })), { deckCloud: 'top', deckInversion: 'subgrid', radiationEvery: 2 }]) {
    const plain = options.levels === 'bl48' && !options.radiationEvery, label = Object.entries(options).map(([key, value]) => `${key} ${value}`).join(', ') + (plain ? ', with the boundary layer diagnosed after the deposit' : ''), run = await deckPair(3, plain ? { ...options, boundaryLayer: PREVIOUS_BOUNDARY_LAYER } : options);
    show(label, run);
    agree(label, run);
    assert.ok(run.moved > 0.5 * run.decked, `${label}: ${run.moved} of ${run.decked} decks off σ 0.92`);
  }
});

test('on bl48 the deck that replaces the resolved cloud at and below its top inside its fraction matches between the engines, on the grid ceiling and the sub-grid inversion, with the radiation every step and held, the grid ceiling\'s with the boundary layer diagnosed after the deposit; outside that, humidity, θ, mixing, longwave heating, entrainment and surface shortwave part the engines no further than CPU twins part themselves, group by group', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const subgrid = { deckCloud: 'top', deckInversion: 'subgrid' };
  for (const options of [{ deckResolved: 'above' }, { deckResolved: 'above', ...subgrid }, { deckResolved: 'above', ...subgrid, radiationEvery: 2 }]) {
    const held = Object.keys(options).length === 1, label = Object.entries(options).map(([key, value]) => `${key} ${value}`).join(', ') + (held ? ', with the boundary layer diagnosed after the deposit' : '');
    const run = await deckPair(3, { levels: 'bl48', cloudColumns: true, cloudSpread: 0.75, boundaryLayer: held ? PREVIOUS_BOUNDARY_LAYER : {}, ...options });
    show(label, run);
    agree(label, run);
    console.log(`${label}: the replaced layers start apart on ${run.fromMiss} columns; worst GPU column off its boundary fluxes by ${run.gpuEnergy.toExponential(2)} W/m²`);
    assert.ok(run.fromMiss === 0 && run.gpuEnergy < 1e-2, `${label}: replaced layers apart on ${run.fromMiss}, ${run.gpuEnergy} W/m² off`);
  }
});

test('on the GPU too the inversion\'s layer hands its heating down inside its column: each column\'s longwave heating sums to its outgoing longwave less the surface\'s net, and with the radiation held its shortwave heating to the sunlight its air absorbs, to single precision', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  for (const options of [...LEVELS.map((levels) => ({ levels })), { radiationEvery: 2 }]) {
    const label = Object.entries(options).map(([key, value]) => `${key} ${value}`).join(', '), run = await deckPair(3, { ...options, deckCloud: 'top', deckInversion: 'subgrid' });
    console.log(`${label}: ${run.gpuSealed} of ${run.C} columns sealed on the GPU, worst column off its boundary fluxes by ${run.gpuEnergy.toExponential(2)} W/m²`);
    assert.ok(run.gpuSealed > 0.5 * run.C && run.gpuEnergy < 1e-2, `${label}: ${run.gpuSealed} sealed, ${run.gpuEnergy} W/m² off`);
  }
});

test('with surface.virtualDensity the deck\'s surface fluxes take the moist density on both engines: the engines agree as closely as without it, and the option moves the deck\'s water alike, with the boundary layer diagnosed after the deposit', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const on = await deckPair(3, { surface: { virtualDensity: true }, boundaryLayer: PREVIOUS_BOUNDARY_LAYER }), off = await deckPair(3, { boundaryLayer: PREVIOUS_BOUNDARY_LAYER });
  show('virtualDensity', on);
  agree('virtualDensity', on);
  const both = new Set(off.kept);
  let gap = 0, move = 0, cells = 0;
  for (const i of on.kept) {
    if (!both.has(i)) continue;
    const cpuMove = on.cpuWater[i] - off.cpuWater[i], gpuMove = on.gpuWater[i] - off.gpuWater[i];
    gap += (gpuMove - cpuMove) ** 2; move += cpuMove ** 2; cells++;
  }
  const rms = (x) => Math.sqrt(x / Math.max(cells, 1)), scale = Math.sqrt(on.cpuWater.reduce((sum, w) => sum + w * w, 0) / on.C);
  console.log(`over ${cells} cells the option moves the deck's water by ${rms(move).toExponential(2)} rms (the water ${scale.toExponential(2)}); the engines' moves differ by ${rms(gap).toExponential(2)}`);
  assert.ok(rms(move) > 1e-4 * scale && gap < 0.01 * move, `the deck water's move ${rms(move)}, its gap between the engines ${rms(gap)}`);
});

test('at the default order the deck matches between the engines on bl60 under the dry boundary layer, on bl48 with its cloud at the top on the sub-grid inversion and on the grid ceiling with the resolved cloud replaced, and with surface.virtualDensity: the deck cell by cell outside two rings of the parted decisions, and humidity, θ, mixing, longwave heating, entrainment and surface shortwave group by group no further than CPU twins part themselves, but for the dry boundary layer\'s mixing, held cell by cell outside the rings', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const means = (label, run) => console.log(`${label} at the default order, means over the ${run.kept.length} kept cells apart by: q ${run.means.q.toExponential(1)} (rms ${run.q.rmsRel.toExponential(1)}), mixing ${run.means.mixing.toExponential(1)} (rms ${run.mixing.rmsRel.toExponential(1)}), w_e ${run.means.entrainment.toExponential(1)} (rms ${run.entrainment.rmsRel.toExponential(1)}), boundary-layer w_e ${run.means.blEntrainment.toExponential(1)} (rms ${run.blEntrainment.rmsRel.toExponential(1)}), surface shortwave ${run.means.sw.toExponential(1)} (rms ${run.sw.rmsRel.toExponential(1)})${run.longwave ? `, longwave heating ${run.means.longwave.toExponential(1)} (rms ${run.longwave.rmsRel.toExponential(1)})` : ''}`);
  for (const [label, options, onMeans, extra] of [
    ['bl60 subgrid, dry', { levels: 'bl60', deckInversion: 'subgrid', turbulence: 'dry' }, ['q'], (run) => assert.ok(run.sealed > 0.5 * run.decked, `${run.sealed} of ${run.decked} sealed alike`)],
    ['bl48 deckCloud top, deckInversion subgrid', { levels: 'bl48', deckCloud: 'top', deckInversion: 'subgrid' }, ['q'], (run) => assert.ok(run.moved > 0.5 * run.decked, `${run.moved} of ${run.decked} decks off σ 0.92`)],
    ['bl48 deckResolved above', { levels: 'bl48', cloudColumns: true, cloudSpread: 0.75, deckResolved: 'above' }, ['q', 'mixing', 'entrainment', 'blEntrainment'], (run) => assert.ok(run.fromMiss === 0 && run.gpuEnergy < 1e-2, `replaced layers apart on ${run.fromMiss}, ${run.gpuEnergy} W/m² off`)],
  ]) {
    const run = await deckPair(3, options);
    show(label, run);
    means(label, run);
    agree(label, run, onMeans);
    extra(run);
  }
  const on = await deckPair(3, { surface: { virtualDensity: true } }), off = await deckPair(3);
  show('virtualDensity', on);
  means('virtualDensity', on);
  agree('virtualDensity', on, ['q', 'blEntrainment', 'sw', 'longwave']);
  const both = new Set(off.kept);
  let gap = 0, move = 0, cells = 0;
  for (const i of on.kept) {
    if (!both.has(i)) continue;
    const cpuMove = on.cpuWater[i] - off.cpuWater[i], gpuMove = on.gpuWater[i] - off.gpuWater[i];
    gap += (gpuMove - cpuMove) ** 2; move += cpuMove ** 2; cells++;
  }
  const rms = (x) => Math.sqrt(x / Math.max(cells, 1)), scale = Math.sqrt(on.cpuWater.reduce((sum, w) => sum + w * w, 0) / on.C);
  console.log(`at the default order over ${cells} cells the option moves the deck's water by ${rms(move).toExponential(2)} rms (the water ${scale.toExponential(2)}); the engines' moves differ by ${rms(gap).toExponential(2)}`);
  assert.ok(rms(move) > 1e-4 * scale && gap < 0.01 * move, `the deck water's move ${rms(move)}, its gap between the engines ${rms(gap)}`);
});
