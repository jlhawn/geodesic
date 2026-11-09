import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { saturationHumidity } from '../js/physics/moist.module.js';
import { decisionTracker, cpuDecisions, gpuDecisions } from './helpers/decisions.mjs';
import { PREVIOUS_LEVELS, PREVIOUS_BOUNDARY_LAYER } from './helpers/previousDefaults.mjs';
import { distributionCheck, magnitude, hemispheres, S9_TIERS } from './helpers/distribution.mjs';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const topography = syntheticTopography(90, 180, (lat, lon) => ((Math.cos(lon) > 0 && Math.abs(lat) < 1.2) || lat < -1.15 ? 300 : -4000));

function stats(cpu, gpu, skip = null) {
  let maxDiff = 0, at = -1, sumSq = 0, sumRef = 0;
  for (let x = 0; x < cpu.length; x++) { if (skip && skip.has(x)) continue; const d = Math.abs(cpu[x] - gpu[x]); if (!Number.isNaN(maxDiff) && !(d <= maxDiff)) { maxDiff = d; at = x; } sumSq += d * d; sumRef += cpu[x] * cpu[x]; }
  return { maxDiff, at, rmsRel: Math.sqrt(sumSq / Math.max(sumRef, 1e-300)) };
}

function seedFog(model) {
  const { mesh } = model, C = mesh.nCells, land = model.geography.land, { K, sigmaMid } = model.core, [pi, theta, , , q, qc] = model.state;
  const coast = new Uint8Array(C);
  for (let e = 0; e < mesh.nEdges; e++) { const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1]; if (land[a] !== land[b]) coast[land[a] ? b : a] = 1; }
  let seeded = 0, cold = 0;
  for (let i = 0; i < C; i++) {
    if (!(land[i] || coast[i]) || i % 3 !== 0) continue;
    const x = (K - 1) * C + i, p = pi[i] * sigmaMid[K - 1], T = theta[x] * Math.pow(p / 1e5, 287.04 / 1004.64);
    q[x] = saturationHumidity(T, p);
    qc[x] = 0.5e-3 * Math.pow(10, (i % 7) / 6);
    seeded++;
    if (mesh.latCell[i] > 1.0) cold++;
  }
  return { seeded, cold };
}

function prepare(model, fog = false) {
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i]) model.state[6][i] = 0;
  if (fog) model.seeded = seedFog(model);
  if (model.load) model.load();
  model.ocean.initialize(model.state[3], model.state[6]);
  model.land.initialize();
  for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i]) { model.land.soil[i] = 40; model.land.snow[i] = model.mesh.latCell[i] > 1.0 ? 5 : 0; }
  model.land.load({ soil: Float64Array.from(model.land.soil), snow: Float64Array.from(model.land.snow), vegetation: new Float64Array(model.mesh.nCells).fill(0.5) });
  return model;
}

/*
 * Each group's gaps against 30 CPU twins of the CPU run per order under S9
 * (test/helpers/distribution.mjs): 20 with U4 and U8 nudges at the start or
 * after the first step and every carried array rounded to f32 each step,
 * and 10 with every carried array nudged by up to 4 f32 ulps each step.
 * The groups are Ts over land and over the sea, and soil water, snow and the
 * surface layer over land; over the sea, where the CPU and every twin hold
 * exactly 0, soil water, snow and the surface layer must be exactly 0. One
 * tier sits above twice the twins: the held surface layer's northern
 * trimmed signed mean, at 2.3.
 */
const LAND_S9 = {
  default: {
    tsLand: {largest: 0.541, counts: [[0.27, 4], [0.135, 6], [0.0676, 6], [0.0338, 6], [0.0169, 14], [0.00845, 18], [0.00422, 28], [0.00211, 50]], hemisphere: {north: 0.00106, south: 0.00164}, bulk: {north: {trimmed: 0.000267, signed: 0.000267}, south: {trimmed: 0.000267, signed: 0.000267}}},
    tsSea: {largest: 0.0189, counts: [[0.00944, 4], [0.00472, 10], [0.00236, 18], [0.00118, 22]], hemisphere: {north: 0.000332, south: 0.000327}, bulk: {north: {trimmed: 0.000278, signed: 0.000278}, south: {trimmed: 0.000278, signed: 0.000278}}},
    soil: {largest: 0.0548, counts: [[0.0273, 4], [0.0136, 4], [0.00684, 4], [0.00342, 8], [0.00171, 12], [0.000855, 14], [0.000427, 28], [0.000213, 44]], hemisphere: {north: 0.000121, south: 0.000147}, bulk: {north: {trimmed: 5.61e-5, signed: 3.82e-5}, south: {trimmed: 5.66e-5, signed: 3.82e-5}}},
    snow: {largest: 9.87e-5, counts: [[4.93e-5, 4], [2.46e-5, 4]], hemisphere: {north: 1.19e-6, south: 9.53e-7}, bulk: {north: {trimmed: 9.53e-7, signed: 9.53e-7}, south: {trimmed: 9.53e-7, signed: 9.53e-7}}},
    surface: {largest: 0.076, counts: [[0.0379, 4], [0.0189, 6], [0.00949, 6], [0.00474, 6], [0.00237, 8], [0.00118, 10], [0.000593, 20], [0.000296, 28], [0.000148, 42], [7.42e-5, 56], [3.71e-5, 82]], hemisphere: {north: 1.99e-5, south: 0.000224}, bulk: {north: {trimmed: 2.93e-7, signed: 4.73e-8}, south: {trimmed: 1.95e-7, signed: 2.45e-8}}},
  },
  held: {
    tsLand: {largest: 0.697, counts: [[0.348, 4], [0.174, 4], [0.087, 4], [0.0435, 8], [0.0217, 10], [0.0108, 10], [0.00544, 16], [0.00272, 20], [0.00136, 40]], hemisphere: {north: 0.00173, south: 0.000489}, bulk: {north: {trimmed: 0.000267, signed: 0.000267}, south: {trimmed: 0.000267, signed: 0.000267}}},
    tsSea: {largest: 0.016, counts: [[0.00799, 4], [0.00399, 6], [0.00199, 6]], hemisphere: {north: 0.000313, south: 0.000301}, bulk: {north: {trimmed: 0.000277, signed: 0.000277}, south: {trimmed: 0.000277, signed: 0.000277}}},
    soil: {largest: 0.0563, counts: [[0.0281, 4], [0.014, 4], [0.00703, 4], [0.00351, 6], [0.00175, 6], [0.000879, 8], [0.000439, 18], [0.000219, 48]], hemisphere: {north: 0.000108, south: 0.000103}, bulk: {north: {trimmed: 5.77e-5, signed: 3.82e-5}, south: {trimmed: 5.62e-5, signed: 3.82e-5}}},
    snow: {largest: 9.29e-5, counts: [[4.64e-5, 4], [2.32e-5, 6]], hemisphere: {north: 1.1e-6, south: 9.52e-7}, bulk: {north: {trimmed: 9.52e-7, signed: 9.52e-7}, south: {trimmed: 9.52e-7, signed: 9.52e-7}}},
    surface: {largest: 0.0289, counts: [[0.0144, 4], [0.0072, 4], [0.0036, 6], [0.0018, 6], [0.0009, 8], [0.00045, 14], [0.000225, 20], [0.000112, 34], [5.62e-5, 52], [2.81e-5, 74], [1.4e-5, 90]], hemisphere: {north: 5.8e-5, south: 3.06e-5}, bulk: {north: {trimmed: 2.62e-7, signed: 5.15e-8}, south: {trimmed: 1.81e-7, signed: 2.49e-8}}},
  },
};

for (const held of [true, false]) test(`eight GPU steps over a continent${held ? ', with the boundary layer diagnosed after the deposit,' : ' at the default order'} track the CPU model within the spread of its f32 twins under S9: the surface temperature over land and sea, soil water, snow and surface layer by their largest gap, a ladder of counts, and full and trimmed means per hemisphere, the vegetation and the coast-bound ocean`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const land = { growthTime: 3 * 3600, declineTime: 2 * 3600, snowDeclineTime: 4 * 3600, percolationTime: 1800, stomatalResistance: 140, growthColdest: 283.15, growthWarmest: 303.15 };
  const cpu = prepare(createModel(new Grid(16), { topography, land, levels: PREVIOUS_LEVELS, boundaryLayer: { ...(held ? PREVIOUS_BOUNDARY_LAYER : {}), velocityScale: 'capped' } }), true);
  const gpu = prepare(await createGpuModel(new Grid(16), { topography, land, levels: PREVIOUS_LEVELS, boundaryLayer: { ...(held ? PREVIOUS_BOUNDARY_LAYER : {}), velocityScale: 'capped' } }), true);
  const decisions = decisionTracker(cpu.core.K, cpu.mesh.nCells, { mix: true });
  for (let n = 0; n < 8; n++) {
    cpu.step(900); await gpu.step(900);
    decisions.check(cpuDecisions(cpu), gpuDecisions(await gpu.gpu.download(), await gpu.gpu.downloadPhysics()));
  }
  await gpu.sync();
  const fogGpu = Float64Array.from((await gpu.gpu.downloadPhysics()).FOG.subarray(0, cpu.mesh.nCells)), fogCpu = Float64Array.from(cpu.moist.fogPrecipitation);
  const fog = stats(fogCpu, fogGpu);
  let fogged = 0, coldFog = 0;
  for (let i = 0; i < cpu.mesh.nCells; i++) if (fogCpu[i] > 0) { fogged++; if (cpu.geography.land[i] && cpu.mesh.latCell[i] > 1.0) coldFog++; }
  console.log(`lowest-layer cloud of 0.5–5 g/kg seeded in ${cpu.seeded.seeded} land and coastal cells (${cpu.seeded.cold} poleward of 1 rad): over eight steps fog reached the ground in ${fogged} cells (${coldFog} cold land cells), at most ${Math.max(...fogCpu).toExponential(2)} kg/m², the engines differ by at most ${fog.maxDiff.toExponential(1)} kg/m² (rms ${fog.rmsRel.toExponential(1)} of it)`);
  assert.ok(fogged > cpu.seeded.seeded / 2 && coldFog > 0, `fog from ${fogged} cells, ${coldFog} cold`);
  assert.ok(fog.maxDiff < 1e-3 * Math.max(...fogCpu) + 1e-6 && fog.rmsRel < 1e-4, `fog differs by ${fog.maxDiff} at ${fog.at}, rms ${fog.rmsRel}`);
  const saved = await gpu.land.serialize();
  const d = await gpu.diagnostics(), dc = cpu.diagnostics();
  const theta = stats(cpu.state[1], gpu.state[1]), vegetation = stats(cpu.land.vegetation, saved.vegetation);
  let moved = 0;
  for (let i = 0; i < cpu.mesh.nCells; i++) if (cpu.geography.land[i]) moved = Math.max(moved, Math.abs(cpu.land.vegetation[i] - 0.5));
  assert.ok(moved > 0.2, `the vegetation moved at most ${moved} from its start`);
  assert.ok(vegetation.maxDiff < 1e-3, `vegetation max ${vegetation.maxDiff} at ${vegetation.at}`);
  if (held) assert.equal(gpu.boundaryLayer.diagnoseFirst, false);
  assert.ok(theta.rmsRel < 1e-4, `θ rms ${theta.rmsRel}`);
  const landCells = [...cpu.geography.land.keys()].filter((i) => cpu.geography.land[i]);
  const seaCells = [...cpu.geography.land.keys()].filter((i) => !cpu.geography.land[i]);
  const groups = { tsLand: [cpu.state[3], gpu.state[3], landCells], tsSea: [cpu.state[3], gpu.state[3], seaCells], soil: [cpu.land.soil, saved.soil, landCells], snow: [cpu.land.snow, saved.snow, landCells], surface: [cpu.land.surface, saved.surface, landCells] };
  const checks = Object.entries(groups).map(([name, [a, b, cells]]) => [name, distributionCheck(cells.map((i) => b[i] - a[i]), hemispheres(cpu.mesh, cells), LAND_S9[held ? 'held' : 'default'][name], magnitude(cells.map((i) => a[i])), { cells, require: S9_TIERS })]);
  console.log(`${held ? 'held' : 'default'} order:\n${checks.map(([name, check]) => `  ${name}: ${check.text}`).join('\n')}\n  complete mixes parted at ${[...decisions.kinds.mix].join(', ') || 'none'}; land T ${dc.landMeanT.toFixed(2)} vs ${d.landMeanT.toFixed(2)}, soil ${dc.soilWater.toFixed(2)} vs ${d.soilWater.toFixed(2)}, ocean h1 ${dc.oceanUpperDepth.toFixed(2)} vs ${d.oceanUpperDepth.toFixed(2)}`);
  const seaWater = ['soil', 'snow', 'surface'].flatMap((name) => seaCells.filter((i) => !(cpu.land[name][i] === 0 && saved[name][i] === 0)).map((i) => `${name} at sea cell ${i}: ${cpu.land[name][i]}/${saved[name][i]}`));
  assert.deepEqual([...checks.flatMap(([name, check]) => check.failed.map((failure) => `${name}: ${failure}`)), ...seaWater], [], 'groups past the twins\x27 spread, or land water at sea');
  assert.ok(Math.abs(dc.landMeanT - d.landMeanT) < 0.01 && Math.abs(dc.soilWater - d.soilWater) < 0.01 && Math.abs(dc.oceanUpperDepth - d.oceanUpperDepth) < 0.2);
  for (let i = 0; i < cpu.mesh.nCells; i++) if (cpu.geography.land[i]) assert.equal(gpu.state[6][i], 0);
});

test('the stratiform share that tapers the boundary layer\'s entrainment is zero over land in both engines, and the engines agree on it and on the entrainment over the sea', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const inverted = (model) => {
    const { K, sigmaMid } = model.core, C = model.mesh.nCells;
    const init = initializeState(model, {});
    for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
    for (let k = 0; k < K; k++) if (sigmaMid[k] < 0.75) for (let i = 0; i < C; i++) model.state[1][k * C + i] += 9;
    if (model.geography) for (let i = 0; i < C; i++) if (model.geography.land[i]) model.state[6][i] = 0;
    if (model.load) model.load();
    return model;
  };
  const cpu = inverted(createModel(new Grid(16), { topography })), gpu = inverted(await createGpuModel(new Grid(16), { topography }));
  const allSea = inverted(createModel(new Grid(16), {}));
  cpu.step(900); await gpu.step(900); allSea.step(900);
  await gpu.sync();
  const ph = await gpu.gpu.downloadPhysics(), C = cpu.mesh.nCells, land = cpu.geography.land;
  const share = stats(cpu.radiation.stratiform, ph.STRAT.subarray(0, C)), entrain = stats(cpu.boundaryLayer.entrainment, ph.ENTRAIN.subarray(0, C));
  let landRamp = 0, seaRamp = 0, landEntraining = 0, lands = 0;
  for (let i = 0; i < C; i++) {
    if (land[i]) {
      lands++;
      assert.equal(cpu.radiation.stratiform[i], 0, `cell ${i}: CPU share over land`);
      assert.equal(ph.STRAT[i], 0, `cell ${i}: GPU share over land`);
      if (allSea.radiation.stratiform[i] > 0) landRamp++;
      if (cpu.boundaryLayer.entrainment[i] > 0) landEntraining++;
    } else if (cpu.radiation.stratiform[i] > 0) seaRamp++;
  }
  console.log(`one step at N=16 over a continent under a 9 K inversion: the share would be positive on ${landRamp} land cells and is on ${seaRamp} sea cells, engines to ${share.maxDiff.toExponential(1)}; ${landEntraining} land cells entrain, w_e engines to ${(1000 * entrain.maxDiff).toExponential(1)} mm/s`);
  assert.ok(landRamp > 0.1 * lands && seaRamp > 0.1 * (C - lands), `share on ${landRamp} of ${lands} land and ${seaRamp} of ${C - lands} sea cells`);
  assert.ok(share.maxDiff < 2e-3, `share differs by ${share.maxDiff} at ${share.at}`);
  assert.ok(entrain.maxDiff < 2e-5, `w_e differs by ${entrain.maxDiff} at ${entrain.at}`);
});

test('both engines give each land cell the same albedo from its surface layer, soil water, vegetation and snow, with the wet-soil darkening by the surface layer, by the root zone and without it', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const run = async (land) => {
    const cpu = createModel(new Grid(16), { topography, land }), gpu = await createGpuModel(new Grid(16), { topography, land });
    const C = cpu.mesh.nCells;
    let s = 7;
    const rnd = () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; };
    const init = initializeState(cpu, {});
    for (let i = 0; i < C; i++) if (cpu.geography.land[i]) init[6][i] = 0;
    cpu.land.initialize();
    const soil = Float64Array.from({ length: C }, (_, i) => (cpu.geography.land[i] ? 300 * rnd() : 0));
    const snow = Float64Array.from({ length: C }, (_, i) => (cpu.geography.land[i] && rnd() < 0.2 ? 30 * rnd() : 0));
    const vegetation = Float64Array.from({ length: C }, () => rnd()), surface = Float64Array.from({ length: C }, () => 15 * rnd());
    for (const m of [cpu, gpu]) {
      for (let a = 0; a < init.length; a++) m.state[a].set(init[a]);
      m.land.load({ soil, snow, vegetation, surface });
    }
    gpu.load();
    const expected = Float64Array.from({ length: C }, (_, i) => (cpu.geography.land[i] ? cpu.land.albedo(i) : NaN));
    await gpu.step(900);
    const ph = await gpu.gpu.downloadPhysics();
    let worst = 0, cells = 0;
    for (let i = 0; i < C; i++) if (cpu.geography.land[i]) { worst = Math.max(worst, Math.abs(ph.ADIF[i] - expected[i])); cells++; }
    gpu.destroy();
    return { worst, cells, expected, land: cpu.geography.land };
  };
  const dark = await run({}), roots = await run({ soilDarkening: 'rootZone' }), plain = await run({ soilDarkening: false }), sheets = await run({ iceSheetAlbedo: 0.7 });
  assert.ok(sheets.worst < 1e-6 && sheets.expected.some((a) => a === 0.7), `a chosen ice-sheet albedo: engines apart by ${sheets.worst}`);
  for (const [name, ran] of [['surface layer', dark], ['root zone', roots]]) {
    let darkened = 0, most = 0;
    for (let i = 0; i < ran.expected.length; i++) if (ran.land[i]) { const d = plain.expected[i] - ran.expected[i]; if (d > 1e-6) darkened++; most = Math.max(most, d); }
    console.log(`${ran.cells} land cells with a random surface layer, soil water, vegetation and snow: albedo engines apart by ${ran.worst.toExponential(1)} darkened by the ${name} and ${plain.worst.toExponential(1)} without; it darkens ${darkened} cells, by up to ${most.toFixed(3)}`);
    assert.ok(ran.worst < 1e-6 && plain.worst < 1e-6, `${name}: engines apart by ${ran.worst} and ${plain.worst}`);
    assert.ok(darkened > ran.cells / 3 && most > 0.05, `${name}: ${darkened} cells darkened by up to ${most}`);
  }
});

test('a trace of snow on land leaves the bucket its own wetness in both engines, a snow cover evaporates at the potential rate', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const evaporation = async (engine, snow) => {
    const model = engine === 'cpu' ? createModel(new Grid(16), { topography }) : await createGpuModel(new Grid(16), { topography });
    const C = model.mesh.nCells, init = initializeState(model, {});
    for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
    for (let i = 0; i < C; i++) if (model.geography.land[i]) { model.state[6][i] = 0; model.state[3][i] = 268; }
    if (model.load) model.load();
    model.ocean.initialize(model.state[3], model.state[6]);
    model.land.load({ soil: new Float64Array(C).fill(40), snow: new Float64Array(C).fill(snow), vegetation: new Float64Array(C).fill(0.5) });
    if (engine === 'cpu') { model.step(900); return Float64Array.from(model.radiation.evaporation); }
    await model.step(900);
    const ph = await model.gpu.downloadPhysics();
    model.destroy();
    return Float64Array.from(ph.EVAP);
  };
  const land = createModel(new Grid(16), { topography }).geography.land, cells = [...land.keys()].filter((i) => land[i]);
  for (const engine of ['cpu', 'gpu']) {
    const [bare, trace, covered] = [await evaporation(engine, 0), await evaporation(engine, 1e-20), await evaporation(engine, 5)];
    let apart = 0, wetter = 0;
    for (const i of cells) { apart = Math.max(apart, Math.abs(trace[i] - bare[i])); if (covered[i] > bare[i]) wetter++; }
    console.log(`${engine}: ${cells.length} land cells; 1e-20 kg/m² of snow moves their evaporation by at most ${apart.toExponential(1)} kg/m²/s, 5 kg/m² raises it in ${wetter}`);
    assert.equal(apart, 0, `${engine}: a trace of snow moves the evaporation by ${apart}`);
    assert.ok(wetter > cells.length / 4, `${engine}: a snow cover raises the evaporation in ${wetter} cells`);
  }
});
