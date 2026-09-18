import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createSigmaCore, sigmaInterfaces, stretchedSigmaInterfaces, CP_DRY, P0 } from '../js/dynamics/sigmaCore.module.js';
import { createRK4Arrays } from '../js/dynamics/integrators.module.js';

const N = +(process.env.SIGMA_TEST_N ?? 8);
const grid = new Grid(N);
const mesh = buildMesh(grid);
const C = mesh.nCells, E = mesh.nEdges;

test(`N=${N}: an atmosphere at rest with horizontally uniform theta stays at rest to roundoff`, () => {
  const core = createSigmaCore(mesh);
  const K = core.K;
  const pi = new Float64Array(C).fill(P0);
  const theta = new Float64Array(K * C);
  for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) theta[k * C + i] = 300 + 200 * (1 - core.sigmaMid[k]);
  const u = new Float64Array(K * E);
  const state = [pi, theta, u];
  const step = createRK4Arrays([C, K * C, K * E]);
  const out = [new Float64Array(C), new Float64Array(K * C), new Float64Array(K * E)];
  core.tendency(state, out);
  assert.equal(Math.max(...out[0].map(Math.abs)), 0);
  assert.equal(Math.max(...out[1].map(Math.abs)), 0);
  assert.equal(Math.max(...out[2].map(Math.abs)), 0);
  for (let n = 0; n < 50; n++) step(core.tendency, state, 300);
  assert.equal(Math.max(...pi.map((p) => Math.abs(p - P0))), 0);
  assert.equal(Math.max(...u.map(Math.abs)), 0);
});

test(`N=${N}: an isentropic column reproduces the analytic hydrostatic geopotential`, () => {
  const core = createSigmaCore(mesh);
  const K = core.K;
  const thetaConst = 320;
  const pi = new Float64Array(C).fill(P0);
  const theta = new Float64Array(K * C).fill(thetaConst);
  core.diagnose(pi, theta);
  const { exnerLayer, geopotential } = core.arrays;
  const kappa = core.diagnostics.kappa;
  let worst = 0;
  for (let k = 0; k < K; k++) {
    const exact = CP_DRY * thetaConst * (1 - exnerLayer[k * C]);
    worst = Math.max(worst, Math.abs(geopotential[k * C] - exact) / Math.max(1, Math.abs(exact)));
    const upper = core.diagnostics.sigmaUpper[k], lower = core.diagnostics.sigmaLower[k];
    const exactMean = (lower ** (1 + kappa) - upper ** (1 + kappa)) / ((1 + kappa) * (lower - upper)) * (P0 / P0) ** kappa;
    assert.ok(Math.abs(exnerLayer[k * C] - exactMean) < 1e-12);
  }
  assert.ok(worst < 1e-12);
});

test('sigma interfaces are monotone from 0 to 1', () => {
  for (const levels of [sigmaInterfaces(), stretchedSigmaInterfaces()]) {
    const K = levels.length - 1;
    assert.equal(levels[0], 0);
    assert.equal(levels[K], 1);
    for (let k = 1; k <= K; k++) assert.ok(levels[k] > levels[k - 1]);
  }
  assert.equal(sigmaInterfaces().length, 28);
  assert.equal(stretchedSigmaInterfaces().length, 23);
});
