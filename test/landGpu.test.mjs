import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const topography = syntheticTopography(90, 180, (lat, lon) => ((Math.cos(lon) > 0 && Math.abs(lat) < 1.2) || lat < -1.15 ? 300 : -4000));

function stats(cpu, gpu) {
  let maxDiff = 0, at = -1, sumSq = 0, sumRef = 0;
  for (let x = 0; x < cpu.length; x++) { const d = Math.abs(cpu[x] - gpu[x]); if (d > maxDiff) { maxDiff = d; at = x; } sumSq += d * d; sumRef += cpu[x] * cpu[x]; }
  return { maxDiff, at, rmsRel: Math.sqrt(sumSq / Math.max(sumRef, 1e-300)) };
}

function prepare(model) {
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i]) model.state[6][i] = 0;
  if (model.load) model.load();
  model.ocean.initialize(model.state[3], model.state[6]);
  model.land.initialize();
  for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i]) { model.land.soil[i] = 40; model.land.snow[i] = model.mesh.latCell[i] > 1.0 ? 5 : 0; }
  model.land.load({ soil: Float64Array.from(model.land.soil), snow: Float64Array.from(model.land.snow), vegetation: new Float64Array(model.mesh.nCells).fill(0.5) });
  return model;
}

test('eight GPU steps over a continent track the CPU model: surface, soil, snow, vegetation, and the coast-bound ocean', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const land = { growthTime: 3 * 3600, declineTime: 2 * 3600, snowDeclineTime: 4 * 3600 };
  const cpu = prepare(createModel(new Grid(6), { topography, land }));
  const gpu = prepare(await createGpuModel(new Grid(6), { topography, land }));
  for (let n = 0; n < 8; n++) { cpu.step(900); await gpu.step(900); }
  await gpu.sync();
  const saved = await gpu.land.serialize();
  const d = await gpu.diagnostics(), dc = cpu.diagnostics();
  const ts = stats(cpu.state[3], gpu.state[3]), theta = stats(cpu.state[1], gpu.state[1]);
  const soil = stats(cpu.land.soil, saved.soil), snow = stats(cpu.land.snow, saved.snow), vegetation = stats(cpu.land.vegetation, saved.vegetation), surface = stats(cpu.land.surface, saved.surface);
  assert.ok(surface.maxDiff < 1e-2, `surface layer max ${surface.maxDiff} at ${surface.at}`);
  let moved = 0;
  for (let i = 0; i < cpu.mesh.nCells; i++) if (cpu.geography.land[i]) moved = Math.max(moved, Math.abs(cpu.land.vegetation[i] - 0.5));
  assert.ok(moved > 0.2, `the vegetation moved at most ${moved} from its start`);
  assert.ok(vegetation.maxDiff < 1e-3, `vegetation max ${vegetation.maxDiff} at ${vegetation.at}`);
  console.log(`eight steps over a continent at N=6: Ts max ${ts.maxDiff.toExponential(1)} K, θ rms ${theta.rmsRel.toExponential(1)}, soil max ${soil.maxDiff.toExponential(1)} kg/m², snow max ${snow.maxDiff.toExponential(1)} kg/m²; land T ${dc.landMeanT.toFixed(2)} vs ${d.landMeanT.toFixed(2)}, soil ${dc.soilWater.toFixed(2)} vs ${d.soilWater.toFixed(2)}, ocean h1 ${dc.oceanUpperDepth.toFixed(2)} vs ${d.oceanUpperDepth.toFixed(2)}`);
  assert.ok(ts.maxDiff < 0.02, `Ts max ${ts.maxDiff} at ${ts.at}`);
  assert.ok(theta.rmsRel < 1e-4, `θ rms ${theta.rmsRel}`);
  assert.ok(soil.maxDiff < 1e-2, `soil max ${soil.maxDiff} at ${soil.at}`);
  assert.ok(snow.maxDiff < 1e-2, `snow max ${snow.maxDiff} at ${snow.at}`);
  assert.ok(Math.abs(dc.landMeanT - d.landMeanT) < 0.01 && Math.abs(dc.soilWater - d.soilWater) < 0.01 && Math.abs(dc.oceanUpperDepth - d.oceanUpperDepth) < 1e-3);
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
  const cpu = inverted(createModel(new Grid(6), { topography })), gpu = inverted(await createGpuModel(new Grid(6), { topography }));
  const allSea = inverted(createModel(new Grid(6), {}));
  cpu.step(900); await gpu.step(900); allSea.step(900);
  await gpu.sync();
  const ph = await gpu.gpu.downloadPhysics(), C = cpu.mesh.nCells, land = cpu.geography.land;
  const share = stats(cpu.radiation.stratiform, ph.STRAT.subarray(0, C)), entrain = stats(cpu.boundaryLayer.entrainment, ph.ENTRAIN.subarray(0, C));
  let landRamp = 0, seaRamp = 0, landEntraining = 0;
  for (let i = 0; i < C; i++) {
    if (land[i]) {
      assert.equal(cpu.radiation.stratiform[i], 0, `cell ${i}: CPU share over land`);
      assert.equal(ph.STRAT[i], 0, `cell ${i}: GPU share over land`);
      if (allSea.radiation.stratiform[i] > 0) landRamp++;
      if (cpu.boundaryLayer.entrainment[i] > 0) landEntraining++;
    } else if (cpu.radiation.stratiform[i] > 0) seaRamp++;
  }
  console.log(`one step at N=6 over a continent under a 9 K inversion: the share would be positive on ${landRamp} land cells and is on ${seaRamp} sea cells, engines to ${share.maxDiff.toExponential(1)}; ${landEntraining} land cells entrain, w_e engines to ${(1000 * entrain.maxDiff).toExponential(1)} mm/s`);
  assert.ok(landRamp > 20 && seaRamp > 20, `share on ${landRamp} land and ${seaRamp} sea cells`);
  assert.ok(share.maxDiff < 2e-3, `share differs by ${share.maxDiff} at ${share.at}`);
  assert.ok(entrain.maxDiff < 2e-5, `w_e differs by ${entrain.maxDiff} at ${entrain.at}`);
});

test('both engines give each land cell the same albedo from its soil water, vegetation and snow, with the wet-soil darkening and without it', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const run = async (land) => {
    const cpu = createModel(new Grid(6), { topography, land }), gpu = await createGpuModel(new Grid(6), { topography, land });
    const C = cpu.mesh.nCells;
    let s = 7;
    const rnd = () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; };
    const init = initializeState(cpu, {});
    for (let i = 0; i < C; i++) if (cpu.geography.land[i]) init[6][i] = 0;
    cpu.land.initialize();
    const soil = Float64Array.from({ length: C }, (_, i) => (cpu.geography.land[i] ? 300 * rnd() : 0));
    const snow = Float64Array.from({ length: C }, (_, i) => (cpu.geography.land[i] && rnd() < 0.2 ? 30 * rnd() : 0));
    const vegetation = Float64Array.from({ length: C }, () => rnd());
    for (const m of [cpu, gpu]) {
      for (let a = 0; a < init.length; a++) m.state[a].set(init[a]);
      m.land.load({ soil, snow, vegetation });
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
  const dark = await run({}), plain = await run({ soilDarkening: false });
  let darkened = 0, most = 0;
  for (let i = 0; i < dark.expected.length; i++) if (dark.land[i]) { const d = plain.expected[i] - dark.expected[i]; if (d > 1e-6) darkened++; most = Math.max(most, d); }
  console.log(`${dark.cells} land cells with random soil water, vegetation and snow: albedo engines apart by ${dark.worst.toExponential(1)} with the darkening and ${plain.worst.toExponential(1)} without; it darkens ${darkened} cells, by up to ${most.toFixed(3)}`);
  assert.ok(dark.worst < 1e-6 && plain.worst < 1e-6, `engines apart by ${dark.worst} and ${plain.worst}`);
  assert.ok(darkened > dark.cells / 3 && most > 0.05, `${darkened} cells darkened by up to ${most}`);
});
