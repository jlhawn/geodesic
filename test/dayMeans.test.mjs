import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { SUMMED } from '../js/physics/radiation.module.js';
import { encodeState, decodeState } from '../js/stateFile.module.js';
import { RADIATION_FIELDS, savedRadiationField } from '../js/physics/regrid.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const SUM_SLOTS = { absorbedSolar: 'ABSSUM', atmosphereSolar: 'ATMSUM', outgoingLongwave: 'OLRSUM', insolation: 'INSSUM', reflectedSolar: 'REFLSUM' };
const STEP_SLOTS = { absorbedSolar: 'ABS', atmosphereSolar: 'ATMSW', outgoingLongwave: 'OLR', insolation: 'INS', reflectedSolar: 'REFL' };
const MEAN_SLOTS = { meanAbsorbedSolar: 'ASRMEAN', meanOutgoingLongwave: 'OLRMEAN', meanPlanetaryAlbedo: 'ALBMEAN' };
const DT = 900;

function stats(cpu, gpu) {
  let maxDiff = 0, sumSq = 0, sumRef = 0;
  for (let x = 0; x < cpu.length; x++) { const d = Math.abs(cpu[x] - gpu[x]); maxDiff = Math.max(maxDiff, d); sumSq += d * d; sumRef += cpu[x] * cpu[x]; }
  return { maxDiff, rmsRel: Math.sqrt(sumSq / Math.max(sumRef, 1e-300)) };
}
const relative = (got, expected) => Math.abs(got - expected) / Math.abs(expected);

async function engines() {
  const cpu = createModel(new Grid(6), { ocean: false }), gpu = await createGpuModel(new Grid(6), { ocean: false });
  const init = initializeState(cpu, {});
  for (let a = 0; a < init.length; a++) { cpu.state[a].set(init[a]); gpu.state[a].set(init[a]); }
  gpu.load();
  let area = 0;
  for (let i = 0; i < cpu.mesh.nCells; i++) area += cpu.mesh.areaCell[i];
  return { cpu, gpu, C: cpu.mesh.nCells, area };
}

async function stepBoth({ cpu, gpu, area }, n) {
  const cpuSteps = [], gpuSteps = [];
  for (let s = 0; s < n; s++) {
    if (cpu) { cpu.step(DT); cpuSteps.push(Object.fromEntries(SUMMED.map((name) => [name, cpu.totals[name] / area]))); }
    if (gpu) {
      await gpu.step(DT);
      const ph = await gpu.gpu.downloadPhysics();
      gpuSteps.push(Object.fromEntries(SUMMED.map((name) => {
        let sum = 0;
        for (let i = 0; i < gpu.mesh.nCells; i++) sum += gpu.mesh.areaCell[i] * ph[STEP_SLOTS[name]][i];
        return [name, sum / area];
      })));
    }
  }
  return { cpuSteps, gpuSteps };
}

function meansOf(steps) {
  const mean = (name) => steps.reduce((s, x) => s + x[name], 0) / steps.length;
  return { absorbedSolar: mean('absorbedSolar'), atmosphereSolar: mean('atmosphereSolar'), outgoingLongwave: mean('outgoingLongwave'), planetaryAlbedo: mean('reflectedSolar') / mean('insolation') };
}

function assertReadout(label, d, steps, tolerance) {
  const expected = meansOf(steps), last = steps[steps.length - 1];
  for (const name of ['absorbedSolar', 'atmosphereSolar', 'outgoingLongwave', 'planetaryAlbedo']) assert.ok(relative(d[name], expected[name]) <= tolerance, `${label}: ${name} ${d[name]} against the mean of the steps ${expected[name]}`);
  for (const name of ['absorbedSolar', 'atmosphereSolar', 'outgoingLongwave']) assert.ok(relative(d.instantaneous[name], last[name]) <= tolerance, `${label}: the last step's ${name} ${d.instantaneous[name]} against ${last[name]}`);
  assert.ok(relative(d.instantaneous.planetaryAlbedo, last.reflectedSolar / last.insolation) <= tolerance, `${label}: the last step's albedo`);
  return expected;
}

test('both engines sum each cell\'s radiation over the steps alike and read out day means: the mean of the per-step global values, the albedo the ratio of the summed reflected to the summed incoming sunlight, the last step\'s kept apart, the sums starting again at each read-out', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const run = await engines(), { cpu, gpu, C } = run;
  const first = await stepBoth(run, 24);
  const device = await gpu.gpu.downloadPhysics();
  const cpuSums = Object.fromEntries(SUMMED.map((name) => [name, Float64Array.from(cpu.radiation.summed[name])]));
  const gpuSums = Object.fromEntries(SUMMED.map((name) => [name, Float64Array.from(device[SUM_SLOTS[name]].subarray(0, C))]));
  const parity = Object.fromEntries(SUMMED.map((name) => [name, stats(cpuSums[name], gpuSums[name])]));
  console.log(`24 steps at N=6, per-cell sums CPU against GPU: ${SUMMED.map((name) => `${name} rms ${parity[name].rmsRel.toExponential(1)} (max ${parity[name].maxDiff.toExponential(1)} W/m²)`).join(', ')}`);
  for (const name of SUMMED) assert.ok(parity[name].rmsRel < 1e-4, `${name}: per-cell rms ${parity[name].rmsRel}`);
  assert.ok(parity.insolation.rmsRel < 1e-6, `the incoming sunlight differs by ${parity.insolation.rmsRel}`);

  const dc = cpu.diagnostics(), dg = await gpu.diagnostics();
  await gpu.sync();
  assertReadout('CPU, 24 steps', dc, first.cpuSteps, 1e-12);
  const expectedGpu = assertReadout('GPU, 24 steps', dg, first.gpuSteps, 2e-6);
  const last = first.cpuSteps[23];
  console.log(`day means over the 24 steps (CPU / GPU): ASR ${dc.absorbedSolar.toFixed(3)} / ${dg.absorbedSolar.toFixed(3)} (atmosphere ${dc.atmosphereSolar.toFixed(3)} / ${dg.atmosphereSolar.toFixed(3)}), OLR ${dc.outgoingLongwave.toFixed(3)} / ${dg.outgoingLongwave.toFixed(3)} W/m², albedo ${dc.planetaryAlbedo.toFixed(4)} / ${dg.planetaryAlbedo.toFixed(4)}; the last step's ASR ${last.absorbedSolar.toFixed(1)}, albedo ${(last.reflectedSolar / last.insolation).toFixed(4)}; GPU against the mean of its own steps to ${relative(dg.absorbedSolar, expectedGpu.absorbedSolar).toExponential(1)}`);
  assert.ok(relative(dg.absorbedSolar, dc.absorbedSolar) < 1e-4 && relative(dg.outgoingLongwave, dc.outgoingLongwave) < 1e-4 && Math.abs(dg.planetaryAlbedo - dc.planetaryAlbedo) < 1e-4, 'the engines agree on the day means');
  assert.ok(Math.abs(dc.absorbedSolar - last.absorbedSolar) > 1, `the day mean ${dc.absorbedSolar} is not the last step's ${last.absorbedSolar}`);

  const after = await gpu.gpu.downloadPhysics();
  for (let i = 0; i < C; i++) {
    for (const name of SUMMED) assert.ok(cpu.radiation.summed[name][i] === 0 && after[SUM_SLOTS[name]][i] === 0, `cell ${i}: ${name} starts again`);
    assert.ok(Math.abs(cpu.radiation.meanAbsorbedSolar[i] - cpuSums.absorbedSolar[i] / 24) <= 1e-12 * cpuSums.absorbedSolar[i] && Math.abs(cpu.radiation.meanOutgoingLongwave[i] - cpuSums.outgoingLongwave[i] / 24) <= 1e-12 * cpuSums.outgoingLongwave[i], `cell ${i}: the CPU's per-cell means`);
    assert.equal(cpu.radiation.meanPlanetaryAlbedo[i], cpuSums.insolation[i] > 0 ? cpuSums.reflectedSolar[i] / cpuSums.insolation[i] : 0, `cell ${i}: the CPU's per-cell albedo`);
    assert.ok(Math.abs(gpu.radiation.meanAbsorbedSolar[i] - gpuSums.absorbedSolar[i] / 24) <= 1e-6 * gpuSums.absorbedSolar[i] && Math.abs(gpu.radiation.meanOutgoingLongwave[i] - gpuSums.outgoingLongwave[i] / 24) <= 1e-6 * gpuSums.outgoingLongwave[i], `cell ${i}: the GPU's per-cell means`);
    assert.ok(Math.abs(gpu.radiation.meanPlanetaryAlbedo[i] - (gpuSums.insolation[i] > 0 ? gpuSums.reflectedSolar[i] / gpuSums.insolation[i] : 0)) <= 1e-6, `cell ${i}: the GPU's per-cell albedo`);
    for (const [name, slot] of Object.entries(MEAN_SLOTS)) assert.equal(gpu.radiation[name][i], after[slot][i], `cell ${i}: ${name} mirrored`);
  }
  const albedo = stats(cpu.radiation.meanPlanetaryAlbedo, gpu.radiation.meanPlanetaryAlbedo);
  assert.ok(albedo.maxDiff < 1e-3, `per-cell albedo apart by ${albedo.maxDiff}`);

  const second = await stepBoth(run, 8);
  const dc2 = cpu.diagnostics(), dg2 = await gpu.diagnostics();
  assertReadout('CPU, the next 8 steps', dc2, second.cpuSteps, 1e-12);
  assertReadout('GPU, the next 8 steps', dg2, second.gpuSteps, 2e-6);
  const all = meansOf([...first.cpuSteps, ...second.cpuSteps]);
  assert.ok(Math.abs(dc2.absorbedSolar - all.absorbedSolar) > 1, `the second read-out ${dc2.absorbedSolar} covers its own 8 steps, not all 32 (${all.absorbedSolar})`);
  console.log(`the next 8 steps read out alone: ASR ${dc2.absorbedSolar.toFixed(3)} / ${dg2.absorbedSolar.toFixed(3)} W/m² against ${all.absorbedSolar.toFixed(3)} over all 32, albedo ${dc2.planetaryAlbedo.toFixed(4)} / ${dg2.planetaryAlbedo.toFixed(4)}`);

  await gpu.sync();
  const kept = Float64Array.from(gpu.radiation.meanAbsorbedSolar), keptCpu = Float64Array.from(cpu.radiation.meanAbsorbedSolar);
  const dc3 = cpu.diagnostics(), dg3 = await gpu.diagnostics();
  await gpu.sync();
  for (const d of [dc3, dg3]) for (const name of ['absorbedSolar', 'atmosphereSolar', 'outgoingLongwave', 'planetaryAlbedo']) assert.equal(d[name], d.instantaneous[name], `with no step since the read-out ${d === dc3 ? 'the CPU' : 'the GPU'} gives the last step's ${name}`);
  assert.deepEqual(cpu.radiation.meanAbsorbedSolar, keptCpu, 'the CPU keeps its per-cell means');
  assert.deepEqual(gpu.radiation.meanAbsorbedSolar, kept, 'the GPU keeps its per-cell means');
});

test('the GPU model sends the per-cell means of a saved state to the device on load, an older state loads them and the sums as zeros, and the first read-out after a load covers only the steps since', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const run = await engines(), { cpu, gpu, C } = run;
  await stepBoth(run, 12);
  cpu.diagnostics(); await gpu.diagnostics();
  await gpu.sync();
  const header = { N: 6, K: gpu.core.K, day: 0, time: gpu.time };
  const state = Object.fromEntries(['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'].map((name, a) => [name, gpu.state[a]]));
  const withMeans = await decodeState(encodeState({ ...header, ...state, ...Object.fromEntries(Object.keys(RADIATION_FIELDS).map((name) => [name, gpu.radiation[name]])) }));
  const older = await decodeState(encodeState({ ...header, ...state }));
  assert.ok(gpu.radiation.meanAbsorbedSolar.some((x) => x > 0) && gpu.radiation.meanPlanetaryAlbedo.some((x) => x > 0), 'the run has per-cell means to save');

  const fresh = await createGpuModel(new Grid(6), { ocean: false });
  ['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'].forEach((name, a) => fresh.state[a].set(withMeans[name]));
  for (const name of Object.keys(RADIATION_FIELDS)) fresh.radiation[name].set(savedRadiationField(withMeans, name, fresh));
  fresh.load();
  let device = await fresh.gpu.downloadPhysics();
  for (let i = 0; i < C; i++) for (const [name, slot] of Object.entries(MEAN_SLOTS)) assert.equal(device[slot][i], Math.fround(gpu.radiation[name][i]), `cell ${i}: the saved ${name} on the device`);

  await stepBoth({ gpu, area: run.area }, 5);
  for (const name of Object.keys(RADIATION_FIELDS)) gpu.radiation[name].set(savedRadiationField(older, name, gpu));
  gpu.load();
  device = await gpu.gpu.downloadPhysics();
  for (let i = 0; i < C; i++) {
    for (const slot of Object.values(MEAN_SLOTS)) assert.equal(device[slot][i], 0, `cell ${i}: ${slot} of an older state`);
    for (const slot of Object.values(SUM_SLOTS)) assert.equal(device[slot][i], 0, `cell ${i}: ${slot} after the load`);
  }
  const idle = await gpu.diagnostics();
  assert.equal(idle.absorbedSolar, idle.instantaneous.absorbedSolar, 'a read-out straight after the load has no steps to average');

  await stepBoth({ cpu, area: run.area }, 5);
  ['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'].forEach((name, a) => cpu.state[a].set(older[name]));
  for (const name of Object.keys(RADIATION_FIELDS)) cpu.radiation[name].set(savedRadiationField(older, name, cpu));
  cpu.restartPrecipitation();
  for (const name of SUMMED) assert.ok(cpu.radiation.summed[name].every((x) => x === 0), `the CPU's ${name} sums after the load`);
  assert.ok(Object.keys(RADIATION_FIELDS).every((name) => cpu.radiation[name].every((x) => x === 0)), 'the CPU\'s per-cell means of an older state');

  const since = await stepBoth(run, 6);
  const dc = cpu.diagnostics(), dg = await gpu.diagnostics();
  assertReadout('CPU, six steps after the load', dc, since.cpuSteps, 1e-12);
  assertReadout('GPU, six steps after the load', dg, since.gpuSteps, 2e-6);
  await gpu.sync();
  assert.ok(gpu.radiation.meanOutgoingLongwave.every((x) => x > 0) && cpu.radiation.meanOutgoingLongwave.every((x) => x > 0), 'the per-cell means of the six steps');
  console.log(`six steps after loading a state saved without the means: ASR ${dc.absorbedSolar.toFixed(3)} / ${dg.absorbedSolar.toFixed(3)}, OLR ${dc.outgoingLongwave.toFixed(3)} / ${dg.outgoingLongwave.toFixed(3)} W/m², albedo ${dc.planetaryAlbedo.toFixed(4)} / ${dg.planetaryAlbedo.toFixed(4)} (CPU / GPU)`);
  fresh.destroy(); gpu.destroy();
});

test('a state carries the per-cell radiation means, regridded at another resolution, and one saved without them gives zeros', async () => {
  const model = createModel(new Grid(6), { ocean: false, physics: false }), target = createModel(new Grid(8), { ocean: false, physics: false });
  const C = model.mesh.nCells, albedo = Float64Array.from({ length: C }, (_, i) => 0.2 + 0.01 * (i % 7));
  const saved = await decodeState(encodeState({ N: 6, K: model.core.K, day: 0, time: 0, pi: model.state[0], meanAbsorbedSolar: albedo.map((a) => 400 * (1 - a)), meanOutgoingLongwave: new Float64Array(C).fill(240), meanPlanetaryAlbedo: albedo }));
  assert.deepEqual(Object.keys(RADIATION_FIELDS), ['meanAbsorbedSolar', 'meanOutgoingLongwave', 'meanPlanetaryAlbedo']);
  for (const name of Object.keys(RADIATION_FIELDS)) model.radiation[name].set(savedRadiationField(saved, name, model));
  for (let i = 0; i < C; i++) assert.equal(model.radiation.meanPlanetaryAlbedo[i], Math.fround(albedo[i]), `albedo of cell ${i}`);
  const moved = savedRadiationField(saved, 'meanOutgoingLongwave', target, model);
  assert.equal(moved.length, target.mesh.nCells);
  assert.ok(moved.every((x) => Math.abs(x - 240) < 1e-9), 'a uniform field regrids to itself');
  const legacy = await decodeState(encodeState({ N: 6, K: model.core.K, day: 0, time: 0, pi: model.state[0] }));
  assert.ok(Object.keys(RADIATION_FIELDS).every((name) => savedRadiationField(legacy, name, model).every((x) => x === 0) && savedRadiationField(null, name, target).every((x) => x === 0)), 'a state saved without them, or a fresh start, starts from zero');
});
