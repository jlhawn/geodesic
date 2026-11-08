import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { createSigmaCore, sigmaInterfaces, dynamicsOptions, standardTemperature, standardSigma, vanLeerFace, DIVERGENCE_DAMPING, DIVERGENCE_SURFACES, VERTICAL_ADVECTION, THETA_FACES, P0 } from '../js/dynamics/sigmaCore.module.js';

const N = 16, dt = 1350 * 16 / N, levels = sigmaInterfaces('bl48'), K = levels.length - 1;
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const model = createModel(new Grid(N), { topography, levels, physics: false });
const { mesh } = model, { sigmaMid } = model.core, C = mesh.nCells, E = mesh.nEdges, phis = model.surfaceGeopotential, g = model.core.diagnostics.g;
const start = initializeState(model, {});
const nu4 = (model.core.spacing / Math.PI) ** 4 / (3 * 3600);
const core = (options) => createSigmaCore(mesh, { levels, splitClosure: true, surfaceGeopotential: phis, ...options });
const largest = (a) => a.reduce((m, x) => Math.max(m, Math.abs(x)), 0);
const difference = (a, b) => largest(a.map((x, n) => x - b[n]));

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
  assert.equal(difference(closed({ nu4, nu4Theta: nu4, temperatureClosure: 'standard' }, state)[4], sigma[4]), 0, 'q takes the same closure');
});

test('verticalAdvection takes \'centered\' or \'vanLeerWater\'; van Leer\'s face is the centered value on a profile straight in its spacing, even or not, the upwind value at an extremum or a step, finite on subnormal differences, and never past the downwind value', () => {
  assert.deepEqual(dynamicsOptions({ verticalAdvection: 'vanLeerWater' }), { divergenceSurfaces: 'sigma', verticalAdvection: 'vanLeerWater' });
  assert.throws(() => dynamicsOptions({ verticalAdvection: 'vanLeer' }), /not vanLeer$/);
  assert.throws(() => core({ verticalAdvection: 'ppm' }), /ppm/);
  assert.equal(vanLeerFace(1, 2, 3, 2.5), 2.5);
  assert.equal(vanLeerFace(3, 2, 1, 1.5), 1.5);
  assert.equal(vanLeerFace(3, 2, 3, 2.5), 2, 'an extremum takes the upwind value');
  assert.equal(vanLeerFace(2, 2, 3, 2.5), 2, 'a step takes the upwind value');
  assert.equal(vanLeerFace(0, 2, 3, 2.5), 2 + (4 / 3) * 0.5, 'r = 2 gives φ = 4/3');
  assert.equal(vanLeerFace(-10, 2, 3, 2.9), 3, 'never past the downwind value');
  const line = (x) => 3 + 2 * x;
  assert.ok(Math.abs(vanLeerFace(line(0.5), line(0.6), line(0.8), line(0.7), 0.1, 0.2) - line(0.7)) < 1e-12, 'uneven spacing keeps a straight profile centered');
  assert.ok(Math.abs(vanLeerFace(line(0.5), line(0.6), line(0.8), line(0.7)) - line(0.7)) > 0.05, 'plain differences would not');
  assert.equal(vanLeerFace(1e-9, 1e-320, 0, 5e-321, 0.05, 0.05), 1e-320, 'a subnormal difference across gives the upwind value');
  const noise = (n) => { const x = Math.sin(12.9898 * n + 78.233) * 43758.5453; return x - Math.floor(x) - 0.5; };
  for (let n = 0; n < 6000; n += 6) {
    const upwind = noise(n), downwind = noise(n + 1), face = vanLeerFace(noise(n + 2), upwind, downwind, upwind + (noise(n + 3) + 0.5) * (downwind - upwind), noise(n + 4) + 0.51, noise(n + 5) + 0.51);
    assert.ok((face - upwind) * (face - downwind) <= 0, `${face} outside ${upwind} to ${downwind}`);
  }
});

test(`on a profile of water straight in each column's Exner function at N=${N} on bl48 'vanLeerWater' transports as 'centered' does away from the top and bottom interfaces`, () => {
  const state = start.map((a) => Float64Array.from(a)), look = core({});
  look.diagnose(state[0], state[1]);
  const exner = look.diagnostics.exnerLayer;
  for (let n = 0; n < K * C; n++) { state[4][n] = 0.02 * exner[n] - 0.012; state[5][n] = 1e-4 * exner[n]; }
  const [centered, limited] = ['centered', 'vanLeerWater'].map((verticalAdvection) => {
    const out = state.map((a) => new Float64Array(a.length));
    core({ verticalAdvection }).tendency(state, out);
    return out;
  });
  for (const s of [4, 5]) {
    const inner = (out) => Array.from(out[s].subarray(2 * C, (K - 2) * C)), scale = largest(inner(centered));
    const apart = largest(inner(limited).map((x, n) => x - inner(centered)[n]));
    console.log(`field ${s}: 'vanLeerWater' against 'centered' ${apart.toExponential(2)} of ${scale.toExponential(2)}`);
    assert.ok(scale > 0 && apart < 1e-9 * scale, `field ${s}: ${apart} against ${scale}`);
  }
});

test(`over a sharp step in water at N=${N} over the terrain the centered faces dry the dry layer above it and moisten the moist layer below it past its value; 'vanLeerWater' does neither, and both keep the water`, () => {
  const step = sigmaMid.findIndex((s) => s > 0.85), q0 = 0.01, dSigma = Float64Array.from({ length: K }, (_, k) => levels[k + 1] - levels[k]);
  const state = start.map((a) => Float64Array.from(a));
  for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) { state[4][k * C + i] = k >= step ? q0 : 0; state[5][k * C + i] = k >= step ? 1e-4 : 0; }
  const rates = Object.fromEntries(VERTICAL_ADVECTION.map((verticalAdvection) => {
    const out = state.map((a) => new Float64Array(a.length));
    core({ verticalAdvection }).tendency(state, out);
    return [verticalAdvection, out];
  }));
  const layer = (out, k) => Array.from(out[4].subarray(k * C, (k + 1) * C));
  const scale = largest(layer(rates.centered, step)), count = (a, keep) => a.filter(keep).length;
  const dried = (v) => count(layer(rates[v], step - 1), (x) => x < 0), moistened = (v) => count(layer(rates[v], step), (x) => x > 1e-9 * scale);
  console.log(`σ ${sigmaMid[step - 1].toFixed(3)} over ${sigmaMid[step].toFixed(3)}: the dry layer dried in ${dried('centered')} of ${C} columns centered, ${dried('vanLeerWater')} van Leer; the moist layer moistened past q0 in ${moistened('centered')} centered, ${moistened('vanLeerWater')} van Leer`);
  assert.ok(dried('centered') > 0 && moistened('centered') > 0, 'the centered faces overshoot on both sides of the step');
  assert.equal(dried('vanLeerWater'), 0);
  assert.equal(moistened('vanLeerWater'), 0);
  for (const v of VERTICAL_ADVECTION) for (const s of [4, 5]) {
    let total = 0, size = 0;
    for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) { const w = mesh.areaCell[i] * state[0][i] * dSigma[k] * rates[v][s][k * C + i]; total += w; size += Math.abs(w); }
    assert.ok(Math.abs(total) < 1e-13 * size, `${v} field ${s}: ${total} against ${size}`);
  }
  assert.equal(difference(rates.vanLeerWater[2], rates.centered[2]), 0, 'the wind does not see it');
  assert.equal(difference(rates.vanLeerWater[1], rates.centered[1]), 0, 'vanLeerWater leaves θ centered');
});

test(`thetaFace takes 'exner' or 'arakawaSuarez'; at N=${N} on bl48 over the terrain c_p θ ΔΠ at the 'arakawaSuarez' face is the hydrostatic step between two dry layers, which the 'exner' face misses, and in moist air only θ's tendency sees the face`, () => {
  assert.deepEqual(dynamicsOptions({ thetaFace: 'arakawaSuarez' }), { divergenceSurfaces: 'sigma', thetaFace: 'arakawaSuarez' });
  assert.throws(() => dynamicsOptions({ thetaFace: 'linear' }), /not linear$/);
  assert.throws(() => core({ thetaFace: 'upwind' }), /upwind/);
  assert.equal(createModel(new Grid(N), { levels, physics: false, dynamics: { thetaFace: 'arakawaSuarez' } }).core.thetaFace, 'arakawaSuarez');
  const state = start.map((a) => Float64Array.from(a));
  state[4].fill(0);
  state[5].fill(0);
  const wet = start.map((a) => Float64Array.from(a)), misses = {}, rates = {};
  for (let n = 0; n < K * C; n++) wet[5][n] = 0.01 * wet[4][n];
  for (const thetaFace of THETA_FACES) {
    const c = core({ thetaFace }), out = state.map((a) => new Float64Array(a.length));
    c.tendency(state, out);
    const { exnerLayer, geopotential, cp } = c.diagnostics, face = c.arrays.thetaLower;
    let miss = 0;
    for (let k = 0; k < K - 1; k++) for (let i = 0; i < C; i++) {
      const n = k * C + i, rise = geopotential[n] - geopotential[n + C];
      miss = Math.max(miss, Math.abs(cp * face[n] * (exnerLayer[n + C] - exnerLayer[n]) - rise) / Math.abs(rise));
    }
    misses[thetaFace] = miss;
    rates[thetaFace] = wet.map((a) => new Float64Array(a.length));
    c.tendency(wet, rates[thetaFace]);
  }
  console.log(`c_p θ ΔΠ against the hydrostatic step: 'exner' ${misses.exner.toExponential(2)} at most, 'arakawaSuarez' ${misses.arakawaSuarez.toExponential(2)}`);
  assert.ok(misses.arakawaSuarez < 1e-12, `'arakawaSuarez' ${misses.arakawaSuarez}`);
  assert.ok(misses.exner > 1e-2, `'exner' ${misses.exner}`);
  for (const [s, name] of [[4, 'q'], [5, 'q_c']]) assert.ok(largest(rates.exner[s]) > 0, `${name} does not move`);
  for (const [s, name] of [[0, 'π'], [2, 'the wind'], [4, 'q'], [5, 'q_c']]) assert.equal(difference(rates.arakawaSuarez[s], rates.exner[s]), 0, `${name} sees the face`);
  assert.ok(difference(rates.arakawaSuarez[1], rates.exner[1]) > 0, 'θ does not see the face');
});
