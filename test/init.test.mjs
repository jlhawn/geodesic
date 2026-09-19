import { SOLAR_CONSTANT } from '../js/physics/radiation.module.js';
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState, equilibriumProfile, surfaceTemperature } from '../js/physics/init.module.js';
import { P0 } from '../js/dynamics/sigmaCore.module.js';

const N = +(process.env.INIT_TEST_N ?? 8);
const DAYS = +(process.env.INIT_TEST_DAYS ?? 3);

test('the radiative-convective equilibrium column is stable, warm at the ground, and in balance', () => {
  const model = createModel(new Grid(2));
  const { core, radiation } = model;
  const { K, C, dSigma, cp, g, exnerLayer, sigmaMid } = core.diagnostics;
  const surfaceT = 305.086;
  const profile = equilibriumProfile(model, { surfaceT });
  for (let k = 0; k < K - 1; k++) assert.ok(profile[k] >= profile[k + 1] * (1 - 1e-9));
  const pi = new Float64Array(C).fill(P0);
  const theta = new Float64Array(K * C);
  for (let k = 0; k < K; k++) theta[k * C] = profile[k];
  core.diagnoseColumn(0, pi, theta);
  const airT = profile[K - 1] * exnerLayer[(K - 1) * C];
  assert.ok(surfaceT - airT > 0 && surfaceT - airT < 20);
  radiation.column(0, P0, theta, surfaceT, 3, radiation.opticalDepth(Math.asin(Math.sqrt(1 / 3))), SOLAR_CONSTANT / 4);
  let column = 0, scale = 0;
  for (let k = 0; k < K; k++) { column += radiation.layerFlux[k]; scale += Math.abs(radiation.layerFlux[k]); }
  const longer = equilibriumProfile(model, { surfaceT, days: 2400 });
  let drift = 0;
  for (let k = 0; k < K; k++) drift = Math.max(drift, Math.abs(longer[k] - profile[k]) / profile[k]);
  console.log(`equilibrium column: surface air ${airT.toFixed(1)} K under a ${surfaceT} K surface; column net flux ${column.toFixed(3)} W/m² (scale ${scale.toFixed(0)}); θ at 850/500/200/50 hPa ≈ ${[0.85, 0.5, 0.2, 0.05].map((s) => profile[sigmaMid.findIndex((m) => m > s)].toFixed(0)).join('/')} K; 1200→2400 day drift ${drift.toExponential(1)}`);
  assert.ok(Math.abs(column) / scale < 1e-3);
  assert.ok(drift < 1e-3);
  assert.ok(Math.abs(surfaceTemperature(0) - 303) < 1e-12 && Math.abs(surfaceTemperature(Math.PI / 2) - 258) < 1e-12);
});

test(`N=${N}: the initial state is balanced, stable, and has the intended mass and pressure structure`, () => {
  const model = createModel(new Grid(N));
  const { mesh, core } = model;
  const [pi, theta, u, surfaceT] = initializeState(model, { seedAmplitude: 0 });
  const C = mesh.nCells, K = core.K;
  let area = 0, mass = 0;
  for (let i = 0; i < C; i++) { area += mesh.areaCell[i]; mass += mesh.areaCell[i] * pi[i]; }
  assert.ok(Math.abs(mass / area - P0) / P0 < 1e-9);
  let equator = 0, pole = 0, nEq = 0, nPole = 0, maxWind = 0;
  for (let i = 0; i < C; i++) {
    const lat = Math.abs(mesh.latCell[i]) * 180 / Math.PI;
    if (lat < 10) { equator += pi[i]; nEq++; }
    if (lat > 80) { pole += pi[i]; nPole++; }
    for (let k = 0; k < K - 1; k++) assert.ok(theta[(k + 1) * C + i] < theta[k * C + i]);
  }
  for (const x of u) maxWind = Math.max(maxWind, Math.abs(x));
  console.log(`init N=${N}: ps equator ${(equator / nEq / 100).toFixed(1)} hPa, pole ${(pole / nPole / 100).toFixed(1)} hPa, max geostrophic wind ${maxWind.toFixed(1)} m/s`);
  assert.ok(pole / nPole > equator / nEq);
  assert.ok(maxWind > 5 && maxWind < 80);
});

test(`N=${N}: the balanced initial state rings quietly for ${DAYS} days with the physics off`, () => {
  const model = createModel(new Grid(N), { physics: false });
  const { mesh, core, state } = model;
  const init = initializeState(model, { seedAmplitude: 0 });
  for (let a = 0; a < 4; a++) state[a].set(init[a]);
  const pi0 = Float64Array.from(state[0]);
  const dt = 450 * 16 / N;
  const perDay = Math.round(86400 / dt);
  const rows = [];
  for (let d = 1; d <= DAYS; d++) {
    for (let n = 0; n < perDay; n++) model.step(dt);
    let worst = 0;
    for (let i = 0; i < mesh.nCells; i++) worst = Math.max(worst, Math.abs(state[0][i] - pi0[i]));
    rows.push(`d${d} max|Δps| ${(worst / 100).toFixed(2)} hPa, max|u| ${model.diagnostics().maxWind.toFixed(1)}`);
  }
  console.log(`ringing N=${N}: ${rows.join(' | ')}`);
  const d = model.diagnostics();
  assert.ok(Number.isFinite(d.maxWind) && d.maxWind < 100);
  assert.ok(Math.abs(d.mass - P0) / P0 < 1e-12);
  let worst = 0;
  for (let i = 0; i < mesh.nCells; i++) worst = Math.max(worst, Math.abs(state[0][i] - pi0[i]));
  assert.ok(worst < 300);
});
