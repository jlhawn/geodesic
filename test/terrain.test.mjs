import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { syntheticTopography, surfaceGeopotential } from '../js/geography.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const T0 = 280, P0 = 101325;
const mountain = syntheticTopography(180, 360, (lat, lon) => {
  const d = Math.hypot(lat - 0.6, (lon - 0.5) * Math.cos(0.6));
  return d < 0.5 ? 1 + 3000 * Math.exp(-((d / 0.18) ** 2)) : -4000;
});

function isothermalAtRest(model) {
  const { mesh, core, state } = model;
  const { K, C, R } = core.diagnostics;
  const phis = model.surfaceGeopotential;
  for (let i = 0; i < C; i++) state[0][i] = P0 * Math.exp(-phis[i] / (R * T0));
  core.diagnose(state[0], state[1], null, null);
  const exner = core.arrays.exnerLayer;
  for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) state[1][k * C + i] = T0 / exner[k * C + i];
  for (let a = 2; a < state.length; a++) state[a].fill(0);
  state[3].fill(T0);
  if (model.load) model.load();
  return model;
}

function maxWind(u) { let m = 0; for (const x of u) m = Math.max(m, Math.abs(x)); return m; }

test('the smoothed surface geopotential is zero over the sea, positive on the mountain, and gentler than the raw elevation', () => {
  const model = createModel(new Grid(8), { physics: false, topography: mountain });
  const { geography, mesh } = model;
  const phis = surfaceGeopotential(mesh, geography);
  let peak = 0, rawPeak = 0;
  for (let i = 0; i < mesh.nCells; i++) {
    if (!geography.land[i]) assert.equal(phis[i], 0);
    peak = Math.max(peak, phis[i] / 9.80616);
    rawPeak = Math.max(rawPeak, geography.elevation[i]);
  }
  assert.ok(peak > 1500 && peak < rawPeak, `peak ${peak.toFixed(0)} m of raw ${rawPeak.toFixed(0)} m`);
});

test('an isothermal atmosphere at rest over a 3 km mountain stays nearly at rest for five days', () => {
  const model = isothermalAtRest(createModel(new Grid(16), { physics: false, topography: mountain }));
  const dt = 1350;
  const mass0 = model.core.mass(model.state[0]);
  let worst = 0;
  for (let n = 0; n < Math.round(5 * 86400 / dt); n++) { model.step(dt); worst = Math.max(worst, maxWind(model.state[2])); }
  const mass = model.core.mass(model.state[0]);
  console.log(`rest over a mountain at N=16: max |u| ${worst.toFixed(3)} m/s over five days, mass drift ${((mass - mass0) / mass0).toExponential(1)}`);
  assert.ok(worst < 0.5, `spurious wind ${worst} m/s`);
  assert.ok(Math.abs(mass - mass0) / mass0 < 1e-12);
});

test('the GPU engine reproduces the CPU engine over terrain', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const cpu = isothermalAtRest(createModel(new Grid(8), { physics: false, topography: mountain }));
  const gpu = isothermalAtRest(await createGpuModel(new Grid(8), { topography: mountain }));
  gpu.gpu.hooks.beforePhysics = null;
  for (let n = 0; n < 6; n++) { cpu.step(1800); await gpu.gpu.step(1800); }
  const state = await gpu.gpu.download();
  let du = 0, dpi = 0;
  const E = cpu.mesh.nEdges, sigmaMid = cpu.core.sigmaMid;
  for (let x = 0; x < cpu.state[2].length; x++) if (sigmaMid[Math.floor(x / E)] > 0.05) du = Math.max(du, Math.abs(cpu.state[2][x] - state[2][x]));
  for (let i = 0; i < cpu.mesh.nCells; i++) dpi = Math.max(dpi, Math.abs(cpu.state[0][i] - state[0][i]));
  console.log(`terrain dynamics at N=8: GPU vs CPU max |Δu| ${du.toExponential(1)} m/s, max |Δπ| ${dpi.toExponential(1)} Pa, CPU max |u| ${maxWind(cpu.state[2]).toFixed(3)}`);
  assert.ok(du < 2e-3 && dpi < 2, `Δu ${du} Δπ ${dpi}`);
});
