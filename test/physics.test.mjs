import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createSigmaCore, P0, CP_DRY } from '../js/dynamics/sigmaCore.module.js';
import { createRadiation, sunDirection, AXIAL_TILT, DAY, YEAR } from '../js/physics/radiation.module.js';
import { createSurface } from '../js/physics/surface.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';

const N = +(process.env.PHYSICS_TEST_N ?? 6);
const grid = new Grid(N);
const mesh = buildMesh(grid);
const core = createSigmaCore(mesh);
const { K, C, E } = core.diagnostics;
const EPS = 1e-12;

function random(seed) {
  let s = seed >>> 0;
  return () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; };
}

function sampleState(seed) {
  const rnd = random(seed);
  const pi = Float64Array.from({ length: C }, () => P0 * (0.95 + 0.1 * rnd()));
  const theta = new Float64Array(K * C);
  for (let i = 0; i < C; i++) for (let k = 0; k < K; k++) theta[k * C + i] = 280 + 250 * (1 - core.sigmaMid[k]) ** 2 + 5 * (rnd() - 0.5);
  const surfaceT = Float64Array.from({ length: C }, () => 270 + 30 * rnd());
  const u = Float64Array.from({ length: K * E }, () => 20 * (rnd() - 0.5));
  return [pi, theta, u, surfaceT];
}

test('the sun follows the equinox, the solstice and the daily rotation', () => {
  const s0 = sunDirection(0);
  assert.ok(Math.abs(s0[0] - 1) < EPS && Math.abs(s0[1]) < EPS && Math.abs(s0[2]) < EPS);
  const solstice = sunDirection(YEAR / 4);
  assert.ok(Math.abs(solstice[2] - Math.sin(AXIAL_TILT)) < 1e-9);
  const noonLater = sunDirection(DAY / 2);
  assert.ok(noonLater[0] < -0.9999 && Math.abs(noonLater[1]) < 1e-9);
});

test('radiation column: layer and surface fluxes sum to absorbed solar minus outgoing longwave', () => {
  const radiation = createRadiation(mesh, core);
  const [pi, theta, u, surfaceT] = sampleState(3);
  core.diagnose(pi, theta);
  radiation.setTime(0.3 * DAY);
  let worst = 0;
  for (let i = 0; i < C; i++) {
    const surfaceFlux = radiation.column(i, pi[i], theta, surfaceT[i], 5);
    let layers = 0, scale = 0;
    for (let k = 0; k < K; k++) { layers += radiation.layerFlux[k]; scale += Math.abs(radiation.layerFlux[k]); }
    const residual = layers + surfaceFlux - (radiation.budget.absorbedSolar - radiation.budget.outgoingLongwave);
    worst = Math.max(worst, Math.abs(residual) / (scale + Math.abs(surfaceFlux)));
    assert.ok(radiation.budget.outgoingLongwave > 0 && radiation.budget.outgoingLongwave < 600);
  }
  assert.ok(worst < EPS);
  const pi0 = new Float64Array(C).fill(P0);
  core.diagnose(pi0, theta);
  radiation.column(0, P0, theta, 288, 5);
  let transmitted = 1;
  for (let k = 0; k < K; k++) transmitted *= 1 - radiation.emissivity[k];
  assert.ok(Math.abs((1 - transmitted) - 0.78) < 1e-12);
});

test('radiation warms the column where it absorbs and cools the slab where it emits', () => {
  const radiation = createRadiation(mesh, core);
  const [pi, theta, u, surfaceT] = sampleState(5);
  surfaceT.fill(320);
  core.diagnose(pi, theta);
  radiation.setTime(DAY / 4);
  const out = [new Float64Array(C), new Float64Array(K * C), new Float64Array(K * E), new Float64Array(C)];
  const wind = new Float64Array(C).fill(5);
  radiation.apply([pi, theta, u, surfaceT], out, wind);
  let night = -1;
  for (let i = 0; i < C; i++) if (radiation.insolation(i) === 0) { night = i; break; }
  assert.ok(night >= 0);
  assert.ok(out[3][night] < 0);
  assert.ok(out[1].every(Number.isFinite) && out[3].every(Number.isFinite));
});

test('convective adjustment leaves a stable column and conserves enthalpy', () => {
  const surface = createSurface(mesh, core);
  const [pi, theta] = sampleState(7);
  for (let i = 0; i < C; i++) for (let k = 0; k < K; k++) theta[k * C + i] = 300 + 40 * core.sigmaMid[k] + 20 * Math.sin(3 * k + i);
  core.diagnose(pi, theta);
  const { exnerLayer, dSigma } = core.diagnostics;
  const enthalpy = (i) => { let h = 0; for (let k = 0; k < K; k++) h += CP_DRY * theta[k * C + i] * exnerLayer[k * C + i] * dSigma[k]; return h; };
  const before = Float64Array.from({ length: C }, (_, i) => enthalpy(i));
  const mixes = surface.convectiveAdjustment(pi, theta);
  assert.ok(mixes > 0);
  for (let i = 0; i < C; i++) {
    assert.ok(Math.abs(enthalpy(i) - before[i]) / before[i] < EPS);
    for (let k = 0; k < K - 1; k++) assert.ok(theta[(k + 1) * C + i] <= theta[k * C + i] * (1 + 1e-9));
  }
});

test('surface and boundary-layer drag only remove kinetic energy', () => {
  const surface = createSurface(mesh, core);
  const state = sampleState(11);
  const [pi, theta, u] = state;
  core.diagnose(pi, theta);
  surface.lowestWindSpeed(u);
  const out = [new Float64Array(C), new Float64Array(K * C), new Float64Array(K * E), new Float64Array(C)];
  surface.apply(state, out);
  for (let k = 0; k < K; k++) {
    let work = 0;
    for (let e = 0; e < E; e++) work += mesh.dcEdge[e] * mesh.dvEdge[e] * u[k * E + e] * out[2][k * E + e];
    if (core.sigmaMid[k] > 0.7 || k === K - 1) assert.ok(work < 0); else assert.equal(work, 0);
  }
});

test('the assembled model steps a uniform atmosphere without blowing up', () => {
  const model = createModel(new Grid(4));
  const init = initializeState(model, { seedAmplitude: 0, geostrophic: false });
  for (let a = 0; a < 4; a++) model.state[a].set(init[a]);
  model.state[0].fill(P0);
  for (let n = 0; n < 20; n++) model.step(600);
  const d = model.diagnostics();
  assert.ok(Number.isFinite(d.maxWind) && d.maxWind < 30);
  assert.ok(Math.abs(d.mass - P0) / P0 < 1e-12);
  assert.ok(d.absorbedSolar > 150 && d.absorbedSolar < 300);
  assert.ok(d.outgoingLongwave > 100 && d.outgoingLongwave < 400);
});
