import { test } from 'node:test';
import assert from 'node:assert/strict';
import { classPoints, water, nearestClass, mixture, densityMaximum, startSurface, classesByLatitude, columnAt, layerAt, thickClasses, nearestLatitude, FOUR_CLASSES, CELSIUS } from '../explain/seaClasses.module.js';
import { SEAWATER, seawaterDensity, salinityForDensity } from '../js/ocean/seawater.module.js';
import { createOcean, LAYER_DENSITIES, THERMOCLINE_DENSITY } from '../js/ocean/layered.module.js';
import { FREEZING_POINT } from '../js/physics/ice.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';

const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} against ${b}`);
const C = (celsius) => celsius + CELSIUS;

test('the class points sit on their label densities, none floored at freezing', () => {
  const points = classPoints();
  assert.equal(points.length, 46);
  for (const p of points) close(seawaterDensity(p.t, p.s), p.rho, 1e-9, `class ${p.rho}`);
  assert.ok(points.every((p) => p.t > FREEZING_POINT));
  const coldest = points.reduce((a, b) => (b.t < a.t ? b : a));
  assert.equal(coldest.rho, 1027.3);
  close(coldest.t - CELSIUS, -1.664, 5e-4, 'coldest label');
  close(points[LAYER_DENSITIES.indexOf(THERMOCLINE_DENSITY)].t - CELSIUS, 19.435, 5e-4, 'the 1024.0 label');
  const four = classPoints(FOUR_CLASSES);
  assert.deepEqual(four.map((p) => p.rho), [1023, 1025, 1026.5, 1026.95]);
  four.forEach((p, n) => close(p.t - CELSIUS, [23.055, 15.229, 5.797, 0.775][n], 5e-4, `four ${p.rho}`));
  assert.equal(points[nearestClass(water(C(25), 35).rho, points)].rho, 1022.5);
});

test('density, expansion and contraction at the prose waters', () => {
  const table = [[-1.8, 34, 1026.473, 0.456], [10, 35, 1026.0, 1.613], [25, 35, 1022.409, 3.053], [28, 34, 1020.703, 3.318]];
  for (const [t, s, rho, alpha] of table) {
    const w = water(C(t), s);
    close(w.rho, rho, 5e-4, `ρ at ${t} °C`);
    close(w.alpha * 1e4, alpha, 5e-4, `α at ${t} °C`);
    close(w.perKelvin, (seawaterDensity(C(t) + 1e-3, s) - seawaterDensity(C(t) - 1e-3, s)) / 2e-3, 1e-9, `∂ρ/∂T at ${t} °C`);
    close(w.perPsu, (seawaterDensity(C(t), s + 1e-3) - seawaterDensity(C(t), s - 1e-3)) / 2e-3, 1e-9, `∂ρ/∂S at ${t} °C`);
  }
  close(water(C(25), 35).alpha / water(C(-1.8), 34).alpha, 6.69, 0.005, 'expansion near 25 °C against freezing');
  close(water(C(25), 35).perKelvin, -0.313, 5e-4, 'a kelvin at 25 °C');
  close(water(C(0), 35).perKelvin, -0.067, 5e-4, 'a kelvin at 0 °C');
  close(water(C(10), 35).perPsu, 0.766, 5e-4, 'a psu at 10 °C');
  let low = Infinity, high = -Infinity;
  for (let t = -1.8; t <= 30; t += 0.1) for (let s = 33; s <= 37; s += 0.05) { const b = water(C(t), s).beta; low = Math.min(low, b); high = Math.max(high, b); }
  assert.ok(low > 6.97e-4 && high < 7.76e-4, `β from ${low} to ${high}`);
});

test('two waters of equal density mix into denser water', () => {
  const a = { t: C(0), s: salinityForDensity(1026.8, C(0)) }, b = { t: C(15), s: salinityForDensity(1026.8, C(15)) };
  close(a.s, 34.541, 5e-4, 'cold water');
  close(b.s, 37.326, 5e-4, 'warm water');
  const m = mixture(a, b), { a0, b0, lambda1, lambda2, nu } = SEAWATER, dT = b.t - a.t, dS = b.s - a.s;
  close(m.rho, 1027.103, 5e-4, 'the mix');
  close(m.excess, (a0 * lambda1 * dT * dT + b0 * lambda2 * dS * dS + 2 * nu * dT * dS) / 8, 1e-9, 'the curvature of the equation of state');
  close(m.overDenser, m.excess, 1e-9, 'equal densities');
  for (const share of [0.1, 0.3, 0.7, 0.9]) assert.ok(mixture(a, b, share).excess > 0 && mixture(a, b, share).excess < m.excess);
});

test('the density maximum lies below freezing in seawater', () => {
  close(densityMaximum(35) - CELSIUS, -6.80, 5e-3, 'at 35 psu');
  close(densityMaximum(0) - CELSIUS, 1.85, 5e-3, 'in fresh water');
  for (const s of [0, 30, 35, 40]) { const t = densityMaximum(s); close(seawaterDensity(t + 1e-3, s), seawaterDensity(t - 1e-3, s), 1e-9, `flat at ${s} psu`); }
  close(salinityForDensity(1026.95, FREEZING_POINT), 34.601, 5e-4, 'freezing water in the 1026.95 class');
});

test('startSurface is initializeState’s sea', () => {
  const model = createModel(buildMesh(new Grid(4)), { ocean: false });
  const state = initializeState(model, { profile: new Float64Array(model.core.diagnostics.K).fill(300) });
  const { surfaceT, ice } = startSurface(model.mesh);
  assert.deepEqual(Array.from(surfaceT), Array.from(state[3]));
  assert.deepEqual(Array.from(ice), Array.from(state[6]));
});

test('the start depends on latitude alone', () => {
  const mesh = buildMesh(new Grid(8)), { surfaceT, ice } = startSurface(mesh), ocean = createOcean(mesh), n = mesh.nCells;
  ocean.initialize(surfaceT, ice);
  const first = new Map();
  let worst = 0;
  for (let i = 0; i < n; i++) {
    const key = mesh.latCell[i].toFixed(9);
    if (!first.has(key)) { first.set(key, i); continue; }
    const j = first.get(key);
    for (let k = 0; k < ocean.layers; k++) worst = Math.max(worst, Math.abs(ocean.h[k * n + i] - ocean.h[k * n + j]));
  }
  assert.ok(worst < 1e-6, `${worst} m`);
});

test('the classes by latitude at the start', () => {
  const [all, four] = classesByLatitude(16);
  assert.equal(all.latitude.length, 289);
  assert.equal(all.layers, 47);
  const L = all.layers, k1024 = Array.from(all.densities).indexOf(1024), k1026 = Array.from(all.densities).indexOf(1026.2);
  const equator = nearestLatitude(all, 0), bases = columnAt(all, equator);
  assert.equal(thickClasses(all, equator), 34);
  close(bases[k1024 - 1], 351.0, 0.05, 'top of 1024.0 on the equator');
  close(all.eta[equator], 1.118, 5e-4, 'sea surface on the equator');
  close(all.temperature[equator * L] - CELSIUS, 28, 1e-9, 'SST on the equator');
  close(bases[k1026 - 1], 954.2, 0.05, 'top of 1026.2 on the equator');
  close(bases[k1026], 1041.3, 0.05, 'base of 1026.2 on the equator');
  [[15, 30], [30, 23], [45, 20], [60, 16], [75, 13], [-75, 13]].forEach(([lat, count]) => assert.equal(thickClasses(all, nearestLatitude(all, lat)), count, `${lat}°`));
  const sixty = nearestLatitude(all, 60), b60 = columnAt(all, sixty);
  close(all.latitude[sixty], 60.09, 0.005, 'the cell near 60°N');
  close(b60[k1026 - 1], 60.2, 0.05, 'top of 1026.2 at 60°N');
  close(b60[k1026], 246.1, 0.05, 'base of 1026.2 at 60°N');
  close(all.temperature[sixty * L + k1026] - CELSIUS, -1, 1e-9, '1026.2 at 60°N');
  close(all.salinity[sixty * L + k1026], 33.708, 5e-4, '1026.2 at 60°N');
  let edge = 0, low = Infinity, high = -Infinity;
  for (let j = 0; j < all.latitude.length; j++) {
    if (all.thickness[j * L + k1024] > 5) edge = Math.max(edge, Math.abs(all.latitude[j]));
    low = Math.min(low, all.eta[j]); high = Math.max(high, all.eta[j]);
    assert.equal(all.thickness[j * L], 60);
    close(columnAt(all, j)[L - 1], 4000 + all.eta[j], 1e-6, 'column');
  }
  close(edge, 25.3, 0.05, 'the 1024.0 class is thick within');
  close(low, -1.543, 5e-4, 'lowest sea surface');
  close(high, 1.118, 5e-4, 'highest sea surface');
  assert.equal(four.layers, 5);
  const fourBases = columnAt(four, equator);
  [60, 220.3, 650.0, 1305.2].forEach((z, k) => close(fourBases[k], z, 0.05, `four-class base ${k}`));
  [1023, 1025, 1026.5].forEach((rho, k) => close(fourBases[k + 1], bases[Array.from(all.densities).indexOf(rho, 1)], 0.7, `the four nearly share the 46’s base of ${rho}`));
  assert.equal(layerAt(all, equator, 30).k, 0);
  assert.equal(layerAt(all, equator, 1000).k, k1026);
  assert.equal(layerAt(all, equator, 4002), null);
});
