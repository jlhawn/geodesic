import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { sigmaInterfaces, SIGMA_GRIDS } from '../js/dynamics/sigmaCore.module.js';

const LEVELS = ['bl60', ...(SIGMA_GRIDS.includes('bl48') ? ['bl48'] : [])];
import { deckProfile } from './helpers/deckProfile.mjs';
import { decisionTracker, cpuDecisions, gpuDecisions, neighbourhood, radiationTraces } from './helpers/decisions.mjs';
import { LATENT_HEAT } from '../js/physics/moist.module.js';
import { PREVIOUS_RADIATION, PREVIOUS_MOIST } from './helpers/previousDefaults.mjs';

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
  for (let x = 0; x < cpu.length; x++) { const d = Math.abs(cpu[x] - gpu[x]); if (d > maxDiff) { maxDiff = d; at = x; } sumSq += d * d; sumRef += cpu[x] * cpu[x]; }
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

async function deckPair(steps, { levels = 'bl60', turbulence = 'moist', ...options } = {}) {
  const physics = { ...PREVIOUS_RADIATION, mixedLayerDeck: true, deckRest: 'depth', minimumInversion: 2, ...options }, sigma = sigmaInterfaces(levels);
  const model = createModel(new Grid(N), { ocean: false, levels: sigma, radiation: physics, moist: MOIST, boundaryLayer: { turbulence }, surface: { exchange: 'roughness' } });
  const { ambiguous } = deckProfile(model), traces = radiationTraces(model);
  const gpu = await createGpuCore(model.mesh, { levels: sigma, nu4: model.core.nu4, nu4Theta: model.core.nu4Theta, divergenceDamping: model.core.divergenceDamping, referenceTheta: meanTheta(model), physics: { ...physics, ...MOIST, turbulence, surfaceExchange: 'roughness', exchangeOptions: {} } });
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
    mixing: columns(bl.mixing, 'MIX'), longwave: turbulence === 'moist' ? columns(r.longwave, 'LWH') : null, theta: keptColumns(model.state[1], state[1]), q: keptColumns(model.state[4], state[4]),
  };
}

const show = (label, run) => console.log(`${label}: deck on ${run.decked} of ${run.C} cells (GPU ${run.gpuDecked}), ${run.parted} cells with moist decisions or the boundary layer's cloud-top run apart (${run.runs} the run; ${run.near} cells left out with their two rings); inversion layer differs on ${run.layerMiss}, deck layer on ${run.homeMiss}; ${run.sealed} sealed alike, ${run.moved} decks off σ 0.92; rms relative differences: share ${run.share.rmsRel.toExponential(1)} (max ${run.share.maxDiff.toExponential(1)}), top ${run.top.rmsRel.toExponential(1)}, height ${run.height.rmsRel.toExponential(1)}, water ${run.water.rmsRel.toExponential(1)}, w_e ${run.entrainment.rmsRel.toExponential(1)}, mixing ${run.mixing.rmsRel.toExponential(1)}, boundary-layer w_e ${run.blEntrainment.rmsRel.toExponential(1)}, θ ${run.theta.rmsRel.toExponential(1)}, q ${run.q.rmsRel.toExponential(1)}, OLR ${run.olr.rmsRel.toExponential(1)}, surface shortwave ${run.sw.rmsRel.toExponential(1)}, longwave heating ${run.longwave ? run.longwave.rmsRel.toExponential(1) : "not stored by the GPU under the dry boundary layer"}`);

function agree(label, run) {
  assert.ok(run.decked > 0.5 * run.C && run.gpuDecked === run.decked, `${label}: deck on ${run.decked} cells, ${run.gpuDecked} on the GPU`);
  assert.ok(run.layerMiss === 0 && run.homeMiss === 0, `${label}: inversion layer apart on ${run.layerMiss}, deck layer on ${run.homeMiss}`);
  assert.ok(run.share.maxDiff < 1e-3 && run.top.rmsRel < 1e-5 && run.height.rmsRel < 1e-5, `${label}: share ${run.share.maxDiff}, top ${run.top.rmsRel}, height ${run.height.rmsRel}`);
  assert.ok(run.water.rmsRel < 1e-4 && run.entrainment.rmsRel < 1e-4 && run.deck.maxDiff < 5e-5 && run.fraction.maxDiff < 1e-3, `${label}: water ${run.water.rmsRel}, w_e ${run.entrainment.rmsRel}, deck ${run.deck.maxDiff}, fraction ${run.fraction.maxDiff}`);
  assert.ok(run.mixing.rmsRel < 1e-4 && run.blEntrainment.rmsRel < 1e-4, `${label}: mixing ${run.mixing.rmsRel}, boundary-layer w_e ${run.blEntrainment.rmsRel}`);
  assert.ok(run.theta.rmsRel < 1e-6 && run.q.rmsRel < 1e-5, `${label}: θ ${run.theta.rmsRel}, q ${run.q.rmsRel}`);
  assert.ok(run.olr.rmsRel < 3e-5 && run.sw.rmsRel < 3e-5 && (!run.longwave || run.longwave.rmsRel < 1e-4), `${label}: OLR ${run.olr.rmsRel}, surface shortwave ${run.sw.rmsRel}, longwave heating ${run.longwave?.rmsRel}`);
}

test('on bl60, and bl48 where it is built, the sub-grid inversion, the deck on it and the boundary layer it seals match between the engines, under the moist and the dry boundary layer', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  for (const [levels, turbulence] of [...LEVELS.map((levels) => [levels, 'moist']), ['bl60', 'dry']]) {
    const run = await deckPair(3, { levels, deckInversion: 'subgrid', turbulence });
    show(`${levels} subgrid, ${turbulence}`, run);
    agree(`${levels} subgrid, ${turbulence}`, run);
    assert.ok(run.sealed > 0.5 * run.decked, `${turbulence}: ${run.sealed} of ${run.decked} sealed alike`);
  }
});

test('on bl60, and bl48 where it is built, the deck\'s cloud at its own top, on the grid ceiling and on the sub-grid inversion, with the radiation every step and held, matches between the engines', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  for (const options of [{ deckCloud: 'top' }, ...LEVELS.map((levels) => ({ levels, deckCloud: 'top', deckInversion: 'subgrid' })), { deckCloud: 'top', deckInversion: 'subgrid', radiationEvery: 2 }]) {
    const label = Object.entries(options).map(([key, value]) => `${key} ${value}`).join(', '), run = await deckPair(3, options);
    show(label, run);
    agree(label, run);
    assert.ok(run.moved > 0.5 * run.decked, `${label}: ${run.moved} of ${run.decked} decks off σ 0.92`);
  }
});

test('on bl48 the deck that replaces the resolved cloud at and below its top inside its fraction matches between the engines, on the grid ceiling and the sub-grid inversion, with the radiation every step and held', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const subgrid = { deckCloud: 'top', deckInversion: 'subgrid' };
  for (const options of [{ deckResolved: 'above' }, { deckResolved: 'above', ...subgrid }, { deckResolved: 'above', ...subgrid, radiationEvery: 2 }]) {
    const label = Object.entries(options).map(([key, value]) => `${key} ${value}`).join(', ');
    const run = await deckPair(3, { levels: 'bl48', cloudColumns: true, cloudSpread: 0.75, ...options });
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
