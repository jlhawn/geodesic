import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { sunDistance, DAY } from '../js/physics/radiation.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuCore } = gpuAvailable ? await import('../js/gpu/core.gpu.js') : {};
const { readBuffer } = gpuAvailable ? await import('../js/gpu/device.module.js') : {};
const N = 16, DT = 1350;

function meanTheta(model) {
  const { K } = model.core, C = model.mesh.nCells, theta = model.state[1];
  return Float64Array.from({ length: K }, (_, k) => { let s = 0, a = 0; for (let i = 0; i < C; i++) { s += model.mesh.areaCell[i] * theta[k * C + i]; a += model.mesh.areaCell[i]; } return s / a; });
}
function stats(cpu, gpu) {
  let maxDiff = 0, sumSq = 0, sumRef = 0, at = -1;
  const errors = new Float64Array(cpu.length);
  for (let x = 0; x < cpu.length; x++) { const d = Math.abs(cpu[x] - gpu[x]); errors[x] = d; if (d > maxDiff) { maxDiff = d; at = x; } sumSq += d * d; sumRef += cpu[x] * cpu[x]; }
  const scale = Math.sqrt(sumRef / cpu.length);
  return { maxDiff, at, rmsRel: Math.sqrt(sumSq / Math.max(sumRef, 1e-300)), p99Rel: errors.sort()[Math.floor(0.99 * (cpu.length - 1))] / Math.max(scale, 1e-300) };
}

async function pair(radiation, { steps = 4, day = 100, offset = 0, batched = false } = {}) {
  const model = createModel(new Grid(N), { ocean: false, radiation: { stratus: true, clearSkyPass: true, ...radiation }, surface: { exchange: 'roughness' } });
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  model.time = day * DAY + offset * DT;
  const gpu = await createGpuCore(model.mesh, { nu4: model.core.nu4, nu4Theta: model.core.nu4Theta, divergenceDamping: model.core.divergenceDamping, referenceTheta: meanTheta(model), physics: { stratus: true, clearSkyPass: true, ...radiation } });
  gpu.upload(model.state);
  gpu.uploadPhysics({});
  const times = Array.from({ length: steps }, (_, n) => model.time + n * DT);
  for (const time of times) { model.step(DT); if (!batched) await gpu.stepModel(DT, time); }
  if (batched) await gpu.batched(async () => { for (const time of times) await gpu.stepModel(DT, time); });
  return { model, gpu, state: await gpu.download(), physics: await gpu.downloadPhysics() };
}

async function bits(gpu) {
  const { device, buffers } = gpu;
  return { S: await readBuffer(device, buffers.S, buffers.S.size, Uint32Array), PH: await readBuffer(device, buffers.PH, buffers.PH.size, Uint32Array) };
}

function measure(label, { model, state, physics }) {
  const C = model.mesh.nCells, r = model.radiation;
  const m = {
    insolation: stats(r.summed.insolation, physics.INSSUM.subarray(0, C)), absorbed: stats(r.summed.absorbedSolar, physics.ABSSUM.subarray(0, C)), olr: stats(r.outgoing, physics.OLR.subarray(0, C)),
    theta: stats(model.state[1], state[1]), ts: stats(model.state[3], state[3]),
  };
  console.log(`${label}: insolation sums rms ${m.insolation.rmsRel.toExponential(1)}; against the CPU at the 99th percentile of the cells: absorbed sums ${m.absorbed.p99Rel.toExponential(1)}, OLR ${m.olr.p99Rel.toExponential(1)}, θ ${m.theta.p99Rel.toExponential(1)}; rms absorbed ${m.absorbed.rmsRel.toExponential(1)}, OLR ${m.olr.rmsRel.toExponential(1)}, θ ${m.theta.rmsRel.toExponential(1)}; Ts max ${m.ts.maxDiff.toExponential(1)} K`);
  return m;
}

function agrees(kepler, circular) {
  assert.ok(kepler.insolation.rmsRel < 1e-6 && circular.insolation.rmsRel < 1e-6, `insolation sums rms ${kepler.insolation.rmsRel} and ${circular.insolation.rmsRel}`);
  for (const key of ['absorbed', 'olr', 'theta']) assert.ok(kepler[key].p99Rel < 3 * circular[key].p99Rel + 1e-5, `${key} 99th percentile ${kepler[key].p99Rel} against ${circular[key].p99Rel} on the circular orbit`);
}

test('the kepler orbit on the GPU matches the CPU model at N = 16 as the circular orbit does, the sunlight to single precision', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const kepler = await pair({ orbit: 'kepler' }), circular = await pair({ orbit: 'circular' });
  const k = measure('kepler, radiation every step', kepler), c = measure('circular, radiation every step', circular);
  const C = kepler.model.mesh.nCells, sums = { kepler: kepler.model.radiation.summed.insolation, circular: circular.model.radiation.summed.insolation };
  let largest = 0;
  for (let i = 0; i < C; i++) if (sums.circular[i] > 0) largest = Math.max(largest, Math.abs(sums.kepler[i] / sums.circular[i] - 1));
  console.log(`day 100: (a/r)² ${sunDistance(kepler.model.time, 'kepler').toFixed(4)}; the kepler cells' sunlight differs from the circular orbit's by up to ${(100 * largest).toFixed(1)} %`);
  assert.ok(largest > 0.02, `the two orbits differ by ${largest}`);
  agrees(k, c);
});

test('the kepler orbit on the GPU holds its distance factor in the held radiation too: a call mid-interval, then the scaled steps', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const options = { steps: 9, day: 200, offset: 1 };
  const held = await pair({ orbit: 'kepler', radiationEvery: 4 }, options), circular = await pair({ orbit: 'circular', radiationEvery: 4 }, options), C = held.model.mesh.nCells;
  const k = measure('kepler, radiationEvery 4', held), c = measure('circular, radiationEvery 4', circular);
  const mu = stats(held.model.radiation.held.mu, held.physics.RADMU.subarray(0, C)), muCircular = stats(circular.model.radiation.held.mu, circular.physics.RADMU.subarray(0, C));
  console.log(`held cosine times distance factor: rms ${mu.rmsRel.toExponential(1)}, max ${mu.maxDiff.toExponential(1)}; circular rms ${muCircular.rmsRel.toExponential(1)}`);
  assert.ok(mu.rmsRel < 1e-5, `RADMU rms ${mu.rmsRel}`);
  agrees(k, c);
});

test('a batched run of the kepler orbit leaves the device where the same steps one at a time do, and the circular orbit leaves the factor at 1', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const stepped = await pair({ orbit: 'kepler', radiationEvery: 2 }, { steps: 6, day: 250 }), batched = await pair({ orbit: 'kepler', radiationEvery: 2 }, { steps: 6, day: 250, batched: true });
  const a = await bits(stepped.gpu), b = await bits(batched.gpu);
  for (const name of ['S', 'PH']) {
    let differ = 0;
    for (let x = 0; x < a[name].length; x++) if (a[name][x] !== b[name][x]) differ++;
    assert.equal(differ, 0, `${name} differs in ${differ} words between the stepped and the batched run`);
  }
  const { buffers, device } = batched.gpu, params = await readBuffer(device, buffers.P, buffers.P.size, Float32Array);
  assert.equal(params.length, 16);
  assert.equal(params[8], Math.fround(sunDistance(250 * DAY + 5 * DT, 'kepler')), 'the last step\'s distance factor in P[8]');
  const circular = await pair({ orbit: 'circular' }, { steps: 2 }), untouched = await readBuffer(circular.gpu.device, circular.gpu.buffers.P, 64, Float32Array);
  assert.equal(untouched[8], 1, 'P[8] stays 1 on the circular orbit');
});
