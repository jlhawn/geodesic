import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { createSigmaCore, sigmaInterfaces, dynamicsOptions, standardTemperature, standardSigma, DIVERGENCE_DAMPING, DIVERGENCE_SURFACES, P0 } from '../js/dynamics/sigmaCore.module.js';

const N = 16, dt = 1350 * 16 / N, levels = sigmaInterfaces('bl48'), K = levels.length - 1;
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const model = createModel(new Grid(N), { topography, levels, physics: false });
const { mesh } = model, { sigmaMid } = model.core, C = mesh.nCells, E = mesh.nEdges, phis = model.surfaceGeopotential, g = model.core.diagnostics.g;
const start = initializeState(model, {});
const nu4 = (model.core.spacing / Math.PI) ** 4 / (3 * 3600);
const core = (options) => createSigmaCore(mesh, { levels, splitClosure: true, surfaceGeopotential: phis, ...options });
const largest = (a) => a.reduce((m, x) => Math.max(m, Math.abs(x)), 0);

function closed(options, state, step = dt) {
  const c = core(options), s = state.map((a) => Float64Array.from(a));
  c.diagnose(s[0], s[1], s[4], s[5]);
  c.phaseClosure(s, 0, K, step);
  return s;
}

function divergentWind(pi) {
  const u = Float64Array.from(start[2]);
  for (let k = 0; k < K; k++) for (let e = 0; e < E; e++) u[k * E + e] += 3 * Math.sin(5 * mesh.latEdge[e] + 0.3 * k) * Math.cos(7 * mesh.xEdge[3 * e]) * (pi[mesh.cellsOnEdge[2 * e]] > 9e4 ? 1 : 0.5);
  return u;
}

test('the dynamics options take divergenceViscosity and temperatureClosure and return them only when given', () => {
  assert.deepEqual(dynamicsOptions(), { divergenceSurfaces: 'sigma' });
  assert.deepEqual(dynamicsOptions({ divergenceViscosity: 6e5, temperatureClosure: 'standard' }), { divergenceSurfaces: 'sigma', divergenceViscosity: 6e5, temperatureClosure: 'standard' });
  assert.throws(() => dynamicsOptions({ divergenceViscosity: -1 }), /divergenceViscosity/);
  assert.throws(() => dynamicsOptions({ temperatureClosure: 'pressure' }), /pressure/);
  assert.throws(() => createSigmaCore(mesh, { levels, divergenceViscosity: 1e5 }), /splitClosure/);
});

test(`divergenceViscosity ν_d = c d²/dt damps as the coefficient c does at that step on both surfaces, and its step grows with Δt at N=${N}`, () => {
  const state = start.map((a) => Float64Array.from(a));
  state[2] = divergentWind(state[0]);
  const viscosity = DIVERGENCE_DAMPING * model.core.spacing ** 2 / dt;
  for (const divergenceSurfaces of DIVERGENCE_SURFACES) {
    const byStep = closed({ divergenceDamping: DIVERGENCE_DAMPING, divergenceSurfaces }, state), byViscosity = closed({ divergenceViscosity: viscosity, divergenceDamping: DIVERGENCE_DAMPING, divergenceSurfaces }, state);
    const change = byStep[2].map((x, n) => x - state[2][n]), gap = byViscosity[2].map((x, n) => x - byStep[2][n]);
    assert.ok(largest(change) > 1e-3, `${divergenceSurfaces}: the damping moves the wind ${largest(change)} m/s`);
    assert.ok(largest(gap) <= 1e-12 * largest(change), `${divergenceSurfaces}: ν_d against c differ by ${largest(gap)} m/s`);
  }
  const once = closed({ divergenceViscosity: viscosity }, state, dt), half = closed({ divergenceViscosity: viscosity }, state, dt / 2);
  const ratio = largest(once[2].map((x, n) => x - state[2][n])) / largest(half[2].map((x, n) => x - state[2][n]));
  assert.ok(Math.abs(ratio - 2) < 1e-9, `the step at Δt is ${ratio} times the step at Δt/2`);
});

test('the US Standard Atmosphere\'s temperature is continuous across its layers and matches its heights', () => {
  for (const z of [0, 5000, 10999, 11001, 25000, 40000, 60000]) {
    const T = standardTemperature(standardSigma(z)), byHeight = z < 11000 ? 288.15 - 0.0065 * z : z < 20000 ? 216.65 : z < 32000 ? 216.65 + 0.001 * (z - 20000) : z < 47000 ? 228.65 + 0.0028 * (z - 32000) : z < 51000 ? 270.65 : 270.65 - 0.0028 * (z - 51000);
    assert.ok(Math.abs(T - byHeight) < 1e-6, `${z} m: ${T} K against ${byHeight}`);
  }
  assert.ok(Math.abs(standardTemperature(1.05) - 288.15 * 1.05 ** 0.190263) < 1e-3, 'below sea level the lowest layer goes on');
});

test(`temperatureClosure 'standard' leaves air in the standard atmosphere at rest over the terrain at N=${N} where 'sigma' heats the high ground, and both keep c_p T`, () => {
  const pi = Float64Array.from(phis, (phi) => P0 * standardSigma(phi / g));
  const look = core({}), theta = new Float64Array(K * C);
  look.diagnose(pi, theta);
  const exner = look.diagnostics.exnerLayer;
  for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) theta[k * C + i] = standardTemperature(pi[i] * sigmaMid[k] / P0) / exner[k * C + i];
  const state = [pi, theta, new Float64Array(K * E), start[3], new Float64Array(K * C), new Float64Array(K * C), start[6]];
  const heat = (s, i) => { let H = 0; for (let k = 0; k < K; k++) H += model.core.diagnostics.cp * s[1][k * C + i] * exner[k * C + i] * s[0][i] * (levels[k + 1] - levels[k]) / g; return H; };
  const heating = {};
  for (const temperatureClosure of ['sigma', 'standard']) {
    const after = closed({ nu4, nu4Theta: nu4, temperatureClosure }, state);
    let high = 0, area = 0, total = 0, before = 0;
    for (let i = 0; i < C; i++) {
      const gain = (heat(after, i) - heat(state, i)) / dt;
      total += mesh.areaCell[i] * gain; before += mesh.areaCell[i] * heat(state, i);
      if (phis[i] > 1500 * g) { high += mesh.areaCell[i] * gain; area += mesh.areaCell[i]; }
    }
    heating[temperatureClosure] = { high: high / area, total, before, change: largest(after[1].map((x, n) => x - theta[n])) };
  }
  console.log(`a standard atmosphere over the terrain, one closure step: cells above 1500 m gain ${heating.sigma.high.toFixed(3)} W/m² along σ, ${heating.standard.high.toExponential(2)} from the standard reference`);
  assert.ok(heating.sigma.high > 1, `along σ the high ground gains ${heating.sigma.high} W/m²`);
  assert.ok(heating.standard.change < 1e-9, `the standard reference moves θ ${heating.standard.change} K`);
  for (const [name, h] of Object.entries(heating)) assert.ok(Math.abs(h.total * dt) <= 1e-12 * h.before, `${name}: c_p T changes by ${h.total * dt} J of ${h.before}`);
});

test(`temperatureClosure 'standard' keeps c_p T on a disturbed state at N=${N} and is the σ closure of T less the reference`, () => {
  const state = start.map((a) => Float64Array.from(a));
  for (let n = 0; n < K * C; n++) state[1][n] *= 1 + 0.01 * Math.sin(1.7 * n);
  const look = core({});
  const total = (s) => { look.diagnose(s[0], s[1]); const exner = look.diagnostics.exnerLayer; let H = 0; for (let i = 0; i < C; i++) for (let k = 0; k < K; k++) H += mesh.areaCell[i] * s[0][i] * (levels[k + 1] - levels[k]) * s[1][k * C + i] * exner[k * C + i]; return H; };
  const sigma = closed({ nu4, nu4Theta: nu4 }, state), standard = closed({ nu4, nu4Theta: nu4, temperatureClosure: 'standard' }, state);
  assert.ok(Math.abs(total(standard) - total(state)) <= 1e-13 * total(state), 'c_p T is kept');
  const gap = largest(standard[1].map((x, n) => x - sigma[1][n])), moved = largest(sigma[1].map((x, n) => x - state[1][n]));
  assert.ok(gap > 0 && gap < moved, `the reference changes the step by ${gap} K of ${moved}`);
  assert.deepEqual(Array.from(closed({ nu4, nu4Theta: nu4, temperatureClosure: 'standard' }, state)[4]), Array.from(sigma[4]), 'q takes the same closure');
});
