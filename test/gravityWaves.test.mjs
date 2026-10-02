import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { gravityWaveSpectrum, GRAVITY_WAVES } from '../js/physics/gravityWaves.module.js';
import { P0 } from '../js/dynamics/sigmaCore.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuCore } = gpuAvailable ? await import('../js/gpu/core.gpu.js') : {};

// A state with a westerly jet that strengthens upward through the stratosphere and a wave on it.
function windyModel(N, options = {}) {
  const model = createModel(new Grid(N), { ocean: false, gravityWaves: {}, ...options });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  const { mesh, core } = model, E = mesh.nEdges;
  for (let k = 0; k < core.K; k++) {
    const boost = core.sigmaMid[k] < 0.2 ? 60 * Math.log(0.2 / core.sigmaMid[k]) / Math.log(0.2 / core.sigmaMid[0]) : 0;
    for (let e = 0; e < E; e++) {
      const x = mesh.xEdge[3 * e], y = mesh.xEdge[3 * e + 1], z = mesh.xEdge[3 * e + 2], lat = Math.asin(z / Math.hypot(x, y, z)), lon = Math.atan2(y, x), n = mesh.nEdge.subarray(3 * e, 3 * e + 3);
      const eastward = boost * Math.cos(lat) ** 2 * (1 + 0.3 * Math.cos(2 * lon));
      model.state[2][k * E + e] += eastward * (-Math.sin(lon) * n[0] + Math.cos(lon) * n[1]);
    }
  }
  return model;
}

test('the source spectrum carries the flux and no net momentum', () => {
  const amplitude = gravityWaveSpectrum(GRAVITY_WAVES);
  assert.equal(amplitude.length, Math.floor(GRAVITY_WAVES.maxSpeed / GRAVITY_WAVES.speedStep));
  assert.ok(Math.abs(2 * amplitude.reduce((s, x) => s + x, 0) - GRAVITY_WAVES.flux) < 1e-15);
  for (let j = 1; j < amplitude.length; j++) assert.ok(amplitude[j] < amplitude[j - 1]);
});

test('each column keeps its momentum: what the drag takes from one layer it gives to another, and the westward waves that rise through a strengthening westerly reach the top layer', () => {
  const model = windyModel(8), { mesh, core, gravityWaves, state } = model;
  const { K, C, dSigma, g } = core.diagnostics;
  gravityWaves.compute(state);
  let worst = 0, scale = 0, topWestward = 0, columns = 0;
  for (let i = 0; i < C; i++) {
    let east = 0, north = 0, absolute = 0;
    for (let k = 0; k < K; k++) {
      const mass = state[0][i] * dSigma[k] / g;
      east += mass * gravityWaves.east[k * C + i]; north += mass * gravityWaves.north[k * C + i];
      absolute += mass * (Math.abs(gravityWaves.east[k * C + i]) + Math.abs(gravityWaves.north[k * C + i]));
    }
    worst = Math.max(worst, Math.abs(east), Math.abs(north)); scale = Math.max(scale, absolute);
    if (Math.abs(mesh.latCell[i]) < 1.2) { columns++; if (gravityWaves.east[i] < 0) topWestward++; }
  }
  console.log(`N=8: the largest column momentum change ${worst.toExponential(1)} Pa against ${scale.toExponential(1)} Pa of absolute deposit; the top layer is pushed westward in ${topWestward} of ${columns} columns within 69 degrees of the equator`);
  assert.ok(worst < 1e-12 * scale, `a column gains ${worst} Pa`);
  assert.ok(topWestward > 0.9 * columns, `westward in ${topWestward} of ${columns}`);
  for (let k = gravityWaves.source; k < K; k++) for (let i = 0; i < C; i++) assert.equal(gravityWaves.east[k * C + i], 0);
});

test('the model returns the kinetic energy the gravity-wave drag removes as heat', () => {
  const model = windyModel(4, { nu4Hours: Infinity, divergenceDamping: 0, boundaryLayer: false, moist: false });
  const { state, mesh, core } = model, { K, C, E, dSigma, g, cp } = core.diagnostics;
  state[0].fill(P0);
  core.diagnose(state[0], state[1]);
  model.gravityWaves.compute(state);
  const kinetic = () => { let s = 0; for (let k = 0; k < K; k++) for (let e = 0; e < E; e++) s += 0.25 * mesh.dcEdge[e] * mesh.dvEdge[e] * (state[0][mesh.cellsOnEdge[2 * e]] + state[0][mesh.cellsOnEdge[2 * e + 1]]) * dSigma[k] / g * state[2][k * E + e] ** 2; return s; };
  const energy = () => { let s = kinetic(); for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) s += mesh.areaCell[i] * state[0][i] * dSigma[k] / g * cp * state[1][k * C + i] * core.diagnostics.exnerLayer[k * C + i]; return s; };
  core.arrays.dissipation.fill(0);
  const before = energy(), kineticBefore = kinetic();
  model.phases.mixMomentum(0, E, 86400);
  const changed = Math.abs(kineticBefore - kinetic());
  model.phases.dissipate(0, C);
  assert.ok(changed > 0);
  assert.ok(Math.abs(energy() - before) < 1e-9 * changed, `energy changed by ${energy() - before} J against ${changed} J of kinetic energy moved`);
});

async function engineAgreement(waves) {
  const model = windyModel(8, { gravityWaves: waves });
  const { K } = model.core, C = model.mesh.nCells;
  const meanTheta = Float64Array.from({ length: K }, (_, k) => { let s = 0, a = 0; for (let i = 0; i < C; i++) { s += model.mesh.areaCell[i] * model.state[1][k * C + i]; a += model.mesh.areaCell[i]; } return s / a; });
  const gpu = await createGpuCore(model.mesh, { nu4: model.core.nu4, nu4Theta: model.core.nu4Theta, divergenceDamping: model.core.divergenceDamping, referenceTheta: meanTheta, gravityWaves: waves });
  gpu.upload(model.state);
  gpu.uploadPhysics();
  const dt = 900, time = model.time;
  model.step(dt);
  await gpu.stepModel(dt, time);
  const [, , u] = await gpu.download(), physics = await gpu.downloadPhysics();
  let worst = 0, scale = 0, uDiff = 0;
  for (let x = 0; x < K * C; x++) {
    for (const [cpu, gpuValue] of [[model.gravityWaves.east[x], physics.GWE[x]], [model.gravityWaves.north[x], physics.GWN[x]]]) { worst = Math.max(worst, Math.abs(cpu - gpuValue)); scale = Math.max(scale, Math.abs(cpu)); }
  }
  let rms = 0;
  for (let x = 0; x < K * C; x++) rms += (model.gravityWaves.east[x] - physics.GWE[x]) ** 2;
  rms = Math.sqrt(rms / (K * C));
  for (let x = 0; x < u.length; x++) uDiff = Math.max(uDiff, Math.abs(model.state[2][x] - u[x]));
  console.log(`one step at N=8 with ${JSON.stringify(waves)}: the accelerations differ by ${(86400 * worst).toExponential(1)} m/s/day at most (rms ${(86400 * rms).toExponential(1)}) of up to ${(86400 * scale).toFixed(2)}; the wind by ${uDiff.toExponential(1)} m/s`);
  assert.ok(rms < 1e-3 * scale, `rms difference ${rms} against ${scale}`);
  assert.ok(uDiff < 1e-2, `wind differs by ${uDiff} m/s`);
  return model;
}

test('the GPU gravity-wave drag matches the CPU model after a step', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  await engineAgreement({ breakingAmplitude: null });
  await engineAgreement({ breakingAmplitude: 0.4 });
});

test('the intermittent breaking amplitude makes waves break below the top layer that otherwise reach it', () => {
  const absoluteAt = (waves) => {
    const model = windyModel(8, { gravityWaves: { ...waves, diagnose: true } }), { gravityWaves, state, mesh } = model;
    gravityWaves.compute(state);
    const C = mesh.nCells;
    let top = 0, source = 0;
    for (let i = 0; i < C; i++) { top += mesh.areaCell[i] * gravityWaves.absoluteFlux[C + i]; source += mesh.areaCell[i] * gravityWaves.absoluteFlux[(gravityWaves.source - 1) * C + i]; }
    return { top, source };
  };
  const mean = absoluteAt({ breakingAmplitude: null }), present = absoluteAt({ breakingAmplitude: 0.4 });
  console.log(`N=8: the share of the flux that rises from the source layer into the top layer, ${(present.top / present.source).toFixed(3)} with B_w 0.4 m²/s², ${(mean.top / mean.source).toFixed(3)} with the grid-box mean flux tested`);
  assert.ok(present.top < 0.9 * mean.top, `top-layer flux ${present.top} against ${mean.top}`);
});
