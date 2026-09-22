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
function stats(cpu, gpu) {
  let maxDiff = 0, sumSq = 0, sumRef = 0, at = -1;
  for (let x = 0; x < cpu.length; x++) { const d = Math.abs(cpu[x] - gpu[x]); if (d > maxDiff) { maxDiff = d; at = x; } sumSq += d * d; sumRef += cpu[x] * cpu[x]; }
  return { maxDiff, at, rms: Math.sqrt(sumSq / cpu.length), rmsRel: Math.sqrt(sumSq / Math.max(sumRef, 1e-300)) };
}
async function pair(N, steps, dt) {
  const model = createModel(new Grid(N), { ocean: false, ice: { oceanDiffusivity: 0, oceanHeatFlux: 0 } });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  const gpu = await createGpuCore(model.mesh, { nu4: model.core.nu4, nu4Theta: model.core.nu4Theta, referenceTheta: meanTheta(model) });
  gpu.upload(model.state);
  gpu.uploadPhysics();
  for (let n = 0; n < steps; n++) { const time = model.time; model.step(dt); await gpu.stepModel(dt, time); }
  return { model, gpu, state: await gpu.download(), physics: await gpu.downloadPhysics() };
}

test('one full GPU step with physics matches the CPU model', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { model, state } = await pair(6, 1, 900);
  const theta = stats(model.state[1], state[1]), q = stats(model.state[4], state[4]), qc = stats(model.state[5], state[5]);
  const ts = stats(model.state[3], state[3]), ice = stats(model.state[6], state[6]), u = stats(model.state[2], state[2]);
  console.log(`one step at N=6: θ rms ${theta.rmsRel.toExponential(1)} max ${theta.maxDiff.toExponential(1)} K; q rms ${q.rmsRel.toExponential(1)} max ${q.maxDiff.toExponential(1)}; qc max ${qc.maxDiff.toExponential(1)}; Ts max ${ts.maxDiff.toExponential(1)} K; ice max ${ice.maxDiff.toExponential(1)} m; wind max ${u.maxDiff.toExponential(1)} m/s`);
  assert.ok(theta.rmsRel < 1e-5, `θ rms ${theta.rmsRel}`);
  assert.ok(theta.maxDiff < 0.05, `θ max ${theta.maxDiff} K at ${theta.at}`);
  assert.ok(q.rmsRel < 5e-4, `q rms ${q.rmsRel}`);
  assert.ok(ts.maxDiff < 0.02, `Ts max ${ts.maxDiff} K at ${ts.at}`);
  assert.ok(ice.maxDiff < 1e-3, `ice max ${ice.maxDiff} m`);
  assert.ok(u.maxDiff < 1e-2, `wind max ${u.maxDiff} m/s`);
  const { gpu } = await pair(6, 1, 900);
  const physics = await gpu.downloadPhysics();
  const olr = stats(model.radiation.outgoing, physics.OLR.subarray(0, model.mesh.nCells)), sw = stats(model.radiation.surfaceShortwave, physics.SWDN.subarray(0, model.mesh.nCells));
  assert.ok(olr.rmsRel < 1e-5 && sw.rmsRel < 1e-5, `per-cell OLR rms ${olr.rmsRel}, surface shortwave rms ${sw.rmsRel}`);
});

test('twelve full GPU steps track the CPU model and its energy budget', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { model, state, physics } = await pair(6, 12, 900);
  const d = model.diagnostics();
  const C = model.mesh.nCells;
  let area = 0, ts = 0, abs = 0, olr = 0, rain = 0;
  for (let i = 0; i < C; i++) { const a = model.mesh.areaCell[i]; area += a; ts += a * state[3][i]; abs += a * physics.ABS[i]; olr += a * physics.OLR[i]; rain += a * physics.RAIN[i]; }
  ts /= area; abs /= area; olr /= area; rain /= area;
  const theta = stats(model.state[1], state[1]), tsStat = stats(model.state[3], state[3]);
  console.log(`twelve steps at N=6: mean Ts ${d.meanSurfaceT.toFixed(3)} vs ${ts.toFixed(3)} K; solar ${d.absorbedSolar.toFixed(2)} vs ${abs.toFixed(2)}; OLR ${d.outgoingLongwave.toFixed(2)} vs ${olr.toFixed(2)} W/m²; θ rms ${theta.rmsRel.toExponential(1)}, Ts max ${tsStat.maxDiff.toExponential(1)} K`);
  assert.ok(Math.abs(d.meanSurfaceT - ts) < 0.02, `mean Ts ${d.meanSurfaceT} vs ${ts}`);
  assert.ok(Math.abs(d.absorbedSolar - abs) < 0.5, `absorbed solar ${d.absorbedSolar} vs ${abs}`);
  assert.ok(Math.abs(d.outgoingLongwave - olr) < 0.5, `OLR ${d.outgoingLongwave} vs ${olr}`);
  assert.ok(theta.rmsRel < 1e-4, `θ rms ${theta.rmsRel}`);
});
