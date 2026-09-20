import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuCore } = gpuAvailable ? await import('../js/gpu/core.gpu.js') : {};

function meanTheta(model) {
  const { K } = model.core, C = model.mesh.nCells, theta = model.state[1];
  return Float64Array.from({ length: K }, (_, k) => { let s = 0, a = 0; for (let i = 0; i < C; i++) { s += model.mesh.areaCell[i] * theta[k * C + i]; a += model.mesh.areaCell[i]; } return s / a; });
}

function compare(label, cpu, gpu, tolerance) {
  let maxDiff = 0, scale = 0;
  for (let x = 0; x < cpu.length; x++) { maxDiff = Math.max(maxDiff, Math.abs(cpu[x] - gpu[x])); scale = Math.max(scale, Math.abs(cpu[x])); }
  const relative = scale > 0 ? maxDiff / scale : maxDiff;
  assert.ok(Number.isFinite(relative) && relative < tolerance, `${label}: max difference ${maxDiff} against scale ${scale} (${relative.toExponential(2)} relative)`);
  return relative;
}

test('the GPU tendency matches the CPU core to single precision', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const model = createModel(new Grid(8), { physics: false });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  const out = model.state.map((s) => new Float64Array(s.length));
  model.core.tendency(model.state, out);
  const gpu = await createGpuCore(model.mesh, { dragCoefficient: 0, topDragDays: 0, referenceTheta: meanTheta(model) });
  gpu.upload(model.state);
  await gpu.tendency();
  const got = await gpu.download(gpu.buffers.K1);
  const r = ['dπ', 'dθ', 'du', '', 'dq', 'dqc'].map((name, a) => (name ? compare(name, out[a], got[a], name === 'du' ? 2e-3 : 1e-4) : 0));
  console.log(`GPU tendency at N=8: relative differences dπ ${r[0].toExponential(1)} dθ ${r[1].toExponential(1)} du ${r[2].toExponential(1)} dq ${r[4].toExponential(1)}`);
});

test('a rest state stays at rest on the GPU', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const model = createModel(new Grid(6), { physics: false });
  const [pi, theta, u] = model.state;
  pi.fill(101325); u.fill(0);
  for (let k = 0; k < model.core.K; k++) for (let i = 0; i < model.mesh.nCells; i++) theta[k * model.mesh.nCells + i] = 300 + 2 * (model.core.K - 1 - k);
  const gpu = await createGpuCore(model.mesh, { dragCoefficient: 0, topDragDays: 0 });
  gpu.upload(model.state);
  for (let n = 0; n < 20; n++) await gpu.step(450);
  const [piAfter, , uAfter] = await gpu.download();
  let maxU = 0, maxDpi = 0;
  for (const x of uAfter) maxU = Math.max(maxU, Math.abs(x));
  for (const x of piAfter) maxDpi = Math.max(maxDpi, Math.abs(x - 101325));
  assert.ok(maxU < 1e-6, `wind ${maxU}`);
  assert.ok(maxDpi < 1e-2, `surface pressure drift ${maxDpi} Pa`);
});

test('twenty GPU steps track twenty CPU steps', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const model = createModel(new Grid(8), { physics: false });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  const gpu = await createGpuCore(model.mesh, { dragCoefficient: 0, topDragDays: 0, nu4: model.core.nu4, nu4Theta: model.core.nu4Theta, referenceTheta: meanTheta(model) });
  gpu.upload(model.state);
  const dt = 900;
  for (let n = 0; n < 20; n++) { model.step(dt); await gpu.step(dt); }
  const got = await gpu.download();
  const piRel = compare('π after 20 steps', model.state[0], got[0], 1e-4);
  const thetaRel = compare('θ after 20 steps', model.state[1], got[1], 1e-4);
  let rms = 0; for (let i = 0; i < model.mesh.nCells; i++) rms += (model.state[0][i] - got[0][i]) ** 2;
  rms = Math.sqrt(rms / model.mesh.nCells);
  let uMax = 0, uDiff = 0; for (let x = 0; x < got[2].length; x++) { uMax = Math.max(uMax, Math.abs(model.state[2][x])); uDiff = Math.max(uDiff, Math.abs(model.state[2][x] - got[2][x])); }
  assert.ok(uDiff < 1e-2 * Math.max(uMax, 1), `wind difference ${uDiff} m/s against ${uMax}`);
  console.log(`GPU vs CPU after 20 steps at N=8: ps RMS ${rms.toFixed(3)} Pa, π ${piRel.toExponential(1)}, θ ${thetaRel.toExponential(1)}, wind ${uDiff.toExponential(1)} m/s of ${uMax.toFixed(1)}`);
});
