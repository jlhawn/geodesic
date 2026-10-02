import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { gravityWaveSpectrum, gravityWaveFlux, gravityWaveSource, GRAVITY_WAVES } from '../js/physics/gravityWaves.module.js';
import { P0, sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';

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

test('the source flux follows Garfinkel et al. (2022) eq. A3, uniform by default, and the source layer descends with latitude', () => {
  const deg = Math.PI / 180, banded = { ...GRAVITY_WAVES, flux: 4.3e-3, equatorialFlux: 2e-3, northFlux: 3.5e-3, southFlux: 3.5e-3, edge: 15, width: 10 };
  for (const sign of [1, -1]) {
    assert.ok(Math.abs(gravityWaveFlux(sign * 5 * deg, banded) - 2e-3) < 1e-15);
    assert.ok(Math.abs(gravityWaveFlux(sign * 12.5 * deg, banded) - 3.15e-3) < 1e-15);
    assert.ok(Math.abs(gravityWaveFlux(sign * 60 * deg, banded) - (4.3e-3 + 0.5 * 3.5e-3 * (1 + Math.tanh(4.5)) + 0.5 * 3.5e-3 * (1 + Math.tanh(-7.5)))) < 1e-15);
  }
  for (let lat = -89; lat <= 89; lat += 7) assert.equal(gravityWaveFlux(lat * deg, GRAVITY_WAVES), GRAVITY_WAVES.flux);
  const levels = sigmaInterfaces('bl34'), mid = Float64Array.from({ length: levels.length - 1 }, (_, k) => 0.5 * (levels[k] + levels[k + 1]));
  const nearest = (sigma) => mid.reduce((best, s, k) => (Math.abs(s - sigma) < Math.abs(mid[best] - sigma) ? k : best), 0);
  assert.equal(gravityWaveSource(mid, 31500, P0, 0, true), nearest(0.315));
  assert.equal(gravityWaveSource(mid, 31500, P0, 60 * deg, true), nearest(Math.sqrt(0.315)));
  assert.equal(gravityWaveSource(mid, 31500, P0, 60 * deg, false), nearest(0.315));
  assert.equal(gravityWaveSource(mid, 31500, P0, 89.9 * deg, true), mid.length - 2);
});

test('on bl36 the flux that rises into the layers above 0.85 hPa is spread over them at one acceleration, whatever their winds, and each column keeps its momentum, with the lid layers tested as well', () => {
  const model = createModel(new Grid(8), { ocean: false, levels: sigmaInterfaces('bl36'), gravityWaves: { breakingAmplitude: null, flux: 1e-9, equatorialFlux: 1e-9 } });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  model.state[2].fill(0);
  const { gravityWaves, core } = model, { C } = core.diagnostics;
  assert.equal(gravityWaves.lid, 2);
  gravityWaves.compute(model.state);
  let launched = 0;
  for (let i = 0; i < C; i++) {
    if (gravityWaves.columns.scale[i] === 0) continue;
    launched++;
    const top = gravityWaves.east[i], next = gravityWaves.east[C + i];
    assert.ok(Math.abs(top - next) < 1e-9 * Math.abs(top) + 1e-30, `column ${i}: ${top} against ${next} m/s²`);
    for (let k = 2; k < core.K; k++) assert.equal(gravityWaves.east[k * C + i], 0);
  }
  assert.ok(launched > 0);
  const bl34 = createModel(new Grid(8), { ocean: false, levels: sigmaInterfaces('bl34'), gravityWaves: {} });
  assert.equal(bl34.gravityWaves.lid, 1);
  const windy = windyModel(8, { levels: sigmaInterfaces('bl36') }), { dSigma, g } = windy.core.diagnostics;
  windy.gravityWaves.compute(windy.state);
  let worst = 0, scale = 0;
  for (let i = 0; i < C; i++) {
    const top = windy.gravityWaves.east[i], next = windy.gravityWaves.east[C + i];
    assert.ok(Math.abs(top - next) <= 1e-9 * Math.abs(top) + 1e-30, `windy column ${i}: the lid layers are pushed at ${top} and ${next} m/s²`);
    let east = 0, absolute = 0;
    for (let k = 0; k < windy.core.K; k++) { const mass = windy.state[0][i] * dSigma[k] / g; east += mass * windy.gravityWaves.east[k * C + i]; absolute += mass * Math.abs(windy.gravityWaves.east[k * C + i]); }
    worst = Math.max(worst, Math.abs(east)); scale = Math.max(scale, absolute);
  }
  assert.ok(worst < 1e-12 * scale, `a column gains ${worst} Pa`);
  const tested = windyModel(8, { levels: sigmaInterfaces('bl36'), gravityWaves: { lidTests: true } });
  tested.gravityWaves.compute(tested.state);
  let apart = 0, testedWorst = 0;
  for (let i = 0; i < C; i++) {
    if (Math.abs(tested.gravityWaves.east[i] - tested.gravityWaves.east[C + i]) > 1e-6 * Math.abs(tested.gravityWaves.east[i])) apart++;
    let east = 0, absolute = 0;
    for (let k = 0; k < tested.core.K; k++) { const mass = tested.state[0][i] * dSigma[k] / g; east += mass * tested.gravityWaves.east[k * C + i]; absolute += mass * Math.abs(tested.gravityWaves.east[k * C + i]); }
    testedWorst = Math.max(testedWorst, Math.abs(east) / Math.max(absolute, 1e-300));
  }
  assert.ok(apart > 0, 'with lidTests the waves that break in a lid layer push it apart from the other');
  assert.ok(testedWorst < 1e-12, `with lidTests a column gains ${testedWorst} of its absolute deposit`);
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
  const kinetic = (u) => { let s = 0; for (let k = 0; k < K; k++) for (let e = 0; e < E; e++) s += 0.25 * mesh.dcEdge[e] * mesh.dvEdge[e] * (state[0][mesh.cellsOnEdge[2 * e]] + state[0][mesh.cellsOnEdge[2 * e + 1]]) * dSigma[k] / g * u[k * E + e] ** 2; return s; };
  const enthalpyChange = (theta) => { let s = 0; for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) s += mesh.areaCell[i] * state[0][i] * dSigma[k] / g * cp * (state[1][k * C + i] - theta[k * C + i]) * core.diagnostics.exnerLayer[k * C + i]; return s; };
  core.arrays.dissipation.fill(0);
  const uBefore = Float64Array.from(state[2]), thetaBefore = Float64Array.from(state[1]);
  model.phases.mixMomentum(0, E, 86400);
  const kineticChange = kinetic(state[2]) - kinetic(uBefore);
  model.phases.dissipate(0, C);
  const heat = enthalpyChange(thetaBefore);
  assert.ok(Math.abs(kineticChange) > 0);
  assert.ok(Math.abs(kineticChange + heat) < 1e-9 * Math.abs(kineticChange), `energy changed by ${kineticChange + heat} J against ${-kineticChange} J of kinetic energy moved`);
});

async function engineAgreement(waves, levels = undefined) {
  const model = windyModel(8, { gravityWaves: waves, ...(levels ? { levels } : {}) });
  const { K } = model.core, C = model.mesh.nCells;
  const meanTheta = Float64Array.from({ length: K }, (_, k) => { let s = 0, a = 0; for (let i = 0; i < C; i++) { s += model.mesh.areaCell[i] * model.state[1][k * C + i]; a += model.mesh.areaCell[i]; } return s / a; });
  const gpu = await createGpuCore(model.mesh, { levels: model.core.levels, nu4: model.core.nu4, nu4Theta: model.core.nu4Theta, divergenceDamping: model.core.divergenceDamping, referenceTheta: meanTheta, gravityWaves: waves });
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
  console.log(`one step at N=8 on ${K} layers with ${JSON.stringify(waves)}: the accelerations differ by ${(86400 * worst).toExponential(1)} m/s/day at most (rms ${(86400 * rms).toExponential(1)}) of up to ${(86400 * scale).toFixed(2)}; the wind by ${uDiff.toExponential(1)} m/s`);
  assert.ok(rms < 1e-3 * scale, `rms difference ${rms} against ${scale}`);
  assert.ok(uDiff < 1e-2, `wind differs by ${uDiff} m/s`);
  return model;
}

test('the GPU gravity-wave drag matches the CPU model after a step', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  await engineAgreement({ breakingAmplitude: null });
  await engineAgreement({ breakingAmplitude: 0.4 });
  await engineAgreement({}, sigmaInterfaces('bl36'));
  await engineAgreement({ lidTests: true }, sigmaInterfaces('bl36'));
});

test('the intermittent breaking amplitude makes waves break below the top layer that otherwise reach it', () => {
  const absoluteAt = (waves) => {
    const model = windyModel(8, { gravityWaves: { ...waves, diagnose: true } }), { gravityWaves, state, mesh } = model;
    gravityWaves.compute(state);
    const C = mesh.nCells;
    let top = 0, source = 0;
    for (let i = 0; i < C; i++) { top += mesh.areaCell[i] * gravityWaves.absoluteFlux[C + i]; source += mesh.areaCell[i] * gravityWaves.absoluteFlux[(gravityWaves.columns.source[i] - 1) * C + i]; }
    return { top, source };
  };
  const mean = absoluteAt({ breakingAmplitude: null }), present = absoluteAt({ breakingAmplitude: 0.4 });
  console.log(`N=8: the share of the flux that rises from the source layer into the top layer, ${(present.top / present.source).toFixed(3)} with B_w 0.4 m²/s², ${(mean.top / mean.source).toFixed(3)} with the grid-box mean flux tested`);
  assert.ok(present.top < 0.9 * mean.top, `top-layer flux ${present.top} against ${mean.top}`);
});
