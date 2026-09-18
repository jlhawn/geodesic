import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState, referenceTheta, surfaceTemperature } from '../js/physics/init.module.js';
import { P0 } from '../js/dynamics/sigmaCore.module.js';

const N = +(process.env.INIT_TEST_N ?? 8);
const DAYS = +(process.env.INIT_TEST_DAYS ?? 3);

test('the reference profile and surface temperature reproduce the A-grid values', () => {
  assert.equal(referenceTheta(0.9537), 283.107);
  assert.equal(referenceTheta(0.0005), 2054.310);
  assert.ok(referenceTheta(0.5) > 300 && referenceTheta(0.5) < 330);
  assert.ok(Math.abs(surfaceTemperature(0) - 303) < 1e-12);
  assert.ok(Math.abs(surfaceTemperature(Math.PI / 2) - 258) < 1e-12);
});

test(`N=${N}: the initial state is balanced, stable, and has the intended mass and pressure structure`, () => {
  const model = createModel(new Grid(N));
  const { mesh, core } = model;
  const [pi, theta, u, surfaceT] = initializeState(mesh, core, { seedAmplitude: 0 });
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
  const init = initializeState(mesh, core, { seedAmplitude: 0 });
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
