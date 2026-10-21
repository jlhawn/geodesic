import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  HYBRID_LAYERS, SURFACE_LAYERS, SURFACE_DEPTH, H_MIN, water, labelOf, labelRange, isopycnalLine, mixture, densityMaximum, startSurface, layersByLatitude,
  thickness, modeAt, layerAt, modeCounts, emptyEverywhere, nearestLatitude, depthFraction, depthAt, CELSIUS,
} from '../explain/seaClasses.module.js';
import { SEAWATER, seawaterDensity, salinityForDensity } from '../js/ocean/seawater.module.js';
import { createOcean } from '../js/ocean/layered.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';

const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} against ${b}`);
const C = (celsius) => celsius + CELSIUS;
const { densities: LABELS, minimum: MINIMUM, depths: ZGRID } = HYBRID_LAYERS;

const unesco = {
  surface: (S, T) => {
    const rw = 999.842594 + 6.793952e-2 * T - 9.095290e-3 * T ** 2 + 1.001685e-4 * T ** 3 - 1.120083e-6 * T ** 4 + 6.536332e-9 * T ** 5;
    const A = 8.24493e-1 - 4.0899e-3 * T + 7.6438e-5 * T ** 2 - 8.2467e-7 * T ** 3 + 5.3875e-9 * T ** 4;
    const B = -5.72466e-3 + 1.0227e-4 * T - 1.6546e-6 * T ** 2;
    return rw + A * S + B * S ** 1.5 + 4.8314e-4 * S * S;
  },
  bulk: (S, T, P) => {
    const p = P / 10;
    const K0 = 19652.21 + 148.4206 * T - 2.327105 * T ** 2 + 1.360477e-2 * T ** 3 - 5.155288e-5 * T ** 4 + (54.6746 - 0.603459 * T + 1.09987e-2 * T ** 2 - 6.1670e-5 * T ** 3) * S + (7.944e-2 + 1.6483e-2 * T - 5.3009e-4 * T ** 2) * S ** 1.5;
    const A = 3.239908 + 1.43713e-3 * T + 1.16092e-4 * T ** 2 - 5.77905e-7 * T ** 3 + (2.2838e-3 - 1.0981e-5 * T - 1.6078e-6 * T ** 2) * S + 1.91075e-4 * S ** 1.5;
    const B = 8.50935e-5 - 6.12293e-6 * T + 5.2787e-8 * T ** 2 + (-9.9348e-7 + 2.0816e-8 * T + 9.1697e-10 * T ** 2) * S;
    return K0 + A * p + B * p * p;
  },
  density(S, T, P = 0) { return this.surface(S, T) / (1 - P / 10 / this.bulk(S, T, P)); },
  lapse: (S, T, P) => (3.5803e-5 + (8.5258e-6 + (-6.8360e-8 + 6.6228e-10 * T) * T) * T) + (1.8932e-6 - 4.2393e-8 * T) * (S - 35)
    + ((1.8741e-8 + (-6.7795e-10 + (8.7330e-12 - 5.4481e-14 * T) * T) * T) + (-1.1351e-10 + 2.7759e-12 * T) * (S - 35)) * P
    + (-4.6206e-13 + (1.8676e-14 - 2.1687e-16 * T) * T) * P * P,
  potential(S, T0, P0, P1, n = 100) {
    let T = T0, P = P0;
    const h = (P1 - P0) / n;
    for (let i = 0; i < n; i++) {
      const k1 = this.lapse(S, T, P), k2 = this.lapse(S, T + 0.5 * h * k1, P + 0.5 * h), k3 = this.lapse(S, T + 0.5 * h * k2, P + 0.5 * h), k4 = this.lapse(S, T + h * k3, P + h);
      T += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
      P += h;
    }
    return T;
  },
  atPressure(S, theta, P) {
    let low = theta - 1, high = theta + 3;
    for (let i = 0; i < 50; i++) { const mid = 0.5 * (low + high); if (this.potential(S, mid, P, 0, 40) < theta) low = mid; else high = mid; }
    return this.density(S, 0.5 * (low + high), P);
  },
};

test('the UNESCO reference the real-seawater numbers are measured against reproduces its published values', () => {
  close(unesco.surface(35, 5), 1027.67547, 5e-5, 'ρ(35, 5, 0)');
  close(unesco.density(35, 5, 10000), 1069.48914, 5e-5, 'ρ(35, 5, 1000 bar)');
  close(unesco.lapse(40, 40, 10000), 3.255976e-4, 1e-9, 'adiabatic gradient');
  close(unesco.potential(40, 40, 10000, 0), 36.89073, 1e-4, 'potential temperature');
});

test('the model’s formula against real seawater', () => {
  const model = (t, s) => seawaterDensity(C(t), s);
  close(unesco.surface(35, 25), 1023.343, 5e-4, 'real water at 25 °C and 35 psu');
  close(unesco.surface(35, 25) - model(25, 35), 0.93, 5e-3, 'the formula’s deficit at 25 °C');
  close(unesco.surface(35, 10) - 1026, 0.95, 5e-3, 'the 1026 reference against real water at 10 °C and 35 psu');
  let low = Infinity, high = -Infinity;
  for (let t = -1.8; t <= 30 + 1e-9; t += 0.2) for (let s = 33; s <= 37 + 1e-9; s += 0.1) { const d = unesco.surface(s, t) - model(t, s); low = Math.min(low, d); high = Math.max(high, d); }
  close(low, 0.87, 5e-3, 'least deficit');
  close(high, 1.08, 5e-3, 'greatest deficit');
  const slope = (f, t, s) => -(f(t + 1e-3, s) - f(t - 1e-3, s)) / 2e-3 / f(t, s);
  close(slope((t, s) => unesco.surface(s, t), -1.8, 34) * 1e4, 0.249, 5e-4, 'real expansion at freezing');
  close(slope((t, s) => model(t, s), -1.8, 34) * 1e4, 0.456, 5e-4, 'the formula’s expansion at freezing');
  let best = 0;
  for (let t = 0; t < 8; t += 0.002) if (unesco.surface(0, t) > unesco.surface(0, best)) best = t;
  close(best, 3.98, 5e-3, 'real fresh water’s density maximum');
  close(densityMaximum(0) - CELSIUS, 1.85, 5e-3, 'the formula’s fresh-water maximum');
});

test('deep water: real seawater needs pressure to order the two, the model’s formula does not', () => {
  const nadw = { t: 2.5, s: 34.92 }, aabw = { t: -0.5, s: 34.66 };
  close(unesco.surface(nadw.s, nadw.t), 1027.866, 5e-4, 'real NADW at the surface');
  close(unesco.surface(nadw.s, nadw.t) - unesco.surface(aabw.s, aabw.t), 0.01, 5e-3, 'the warmer heavier at the surface');
  close(unesco.atPressure(aabw.s, aabw.t, 4000) - unesco.atPressure(nadw.s, nadw.t, 4000), 0.33, 5e-3, 'the colder heavier at 4000 dbar');
  const wn = water(C(nadw.t), nadw.s), wa = water(C(aabw.t), aabw.s);
  close(wn.rho, 1026.901, 5e-4, 'NADW in the formula');
  close(wa.rho, 1026.926, 5e-4, 'AABW in the formula');
  close(wa.rho - wn.rho, 0.024, 5e-4, 'the colder heavier in the formula');
  assert.equal(LABELS[labelOf(wn.rho)], 1026.89);
  assert.equal(LABELS[labelOf(wa.rho)], 1026.92);
  assert.equal(labelOf(wa.rho) - labelOf(wn.rho), 1);
  const squeeze = (unesco.atPressure(nadw.s, nadw.t, 4000) / unesco.surface(nadw.s, nadw.t) - 1) * 100;
  assert.ok(squeeze > 1.6 && squeeze < 2, `${squeeze} percent denser 4 km down`);
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
  const warm = water(C(25), 35), cold = water(C(-1.8), 34);
  close(warm.alpha / cold.alpha, 6.69, 0.005, 'expansion near 25 °C against freezing');
  close(warm.perKelvin, -0.313, 5e-4, 'a kelvin at 25 °C');
  close(warm.perPsu, 0.729, 5e-4, 'a psu at 25 °C');
  close(warm.perPsu / -warm.perKelvin, 2.3, 0.05, 'a psu in kelvins at 25 °C');
  close(cold.perKelvin, -0.047, 5e-4, 'a kelvin at freezing');
  close(cold.perPsu, 0.795, 5e-4, 'a psu at freezing');
  close(cold.perPsu / -cold.perKelvin, 17, 0.2, 'a psu in kelvins at freezing');
  close(water(C(0), 35).perKelvin, -0.067, 5e-4, 'a kelvin at 0 °C');
  close(water(C(10), 35).perPsu, 0.766, 5e-4, 'a psu at 10 °C');
  let low = Infinity, high = -Infinity;
  for (let t = -1.8; t <= 30; t += 0.1) for (let s = 33; s <= 37; s += 0.05) { const b = water(C(t), s).beta; low = Math.min(low, b); high = Math.max(high, b); }
  assert.ok(low > 6.97e-4 && high < 7.76e-4, `β from ${low} to ${high}`);
});

test('two waters of equal density mix into denser water, and into a layer seven below theirs', () => {
  const a = { t: C(0), s: salinityForDensity(1026.8, C(0)) }, b = { t: C(15), s: salinityForDensity(1026.8, C(15)) };
  close(a.s, 34.541, 5e-4, 'cold water');
  close(b.s, 37.326, 5e-4, 'warm water');
  const m = mixture(a, b), { a0, b0, lambda1, lambda2, nu } = SEAWATER, dT = b.t - a.t, dS = b.s - a.s;
  close(m.t - CELSIUS, 7.5, 1e-9, 'the mix’s temperature');
  close(m.s, 35.933, 5e-4, 'the mix’s salinity');
  close(m.rho, 1027.103, 5e-4, 'the mix');
  close(m.excess, (a0 * lambda1 * dT * dT + b0 * lambda2 * dS * dS + 2 * nu * dT * dS) / 8, 1e-9, 'the curvature of the equation of state');
  close(m.excess, 0.303, 5e-4, 'denser than either water');
  close(m.overDenser, m.excess, 1e-9, 'equal densities');
  for (const share of [0.1, 0.3, 0.7, 0.9]) assert.ok(mixture(a, b, share).excess > 0 && mixture(a, b, share).excess < m.excess);
  const ka = labelOf(seawaterDensity(a.t, a.s)), kb = labelOf(seawaterDensity(b.t, b.s)), km = labelOf(m.rho);
  assert.equal(ka, kb);
  assert.equal(LABELS[ka], 1026.8);
  assert.equal(LABELS[km], 1027.1);
  assert.equal(km - ka, 7);
});

test('the density maximum lies below freezing in seawater', () => {
  close(densityMaximum(35) - CELSIUS, -6.80, 5e-3, 'at 35 psu');
  close(densityMaximum(0) - CELSIUS, 1.85, 5e-3, 'in fresh water');
  for (const s of [0, 30, 35, 40]) { const t = densityMaximum(s); close(seawaterDensity(t + 1e-3, s), seawaterDensity(t - 1e-3, s), 1e-9, `flat at ${s} psu`); }
});

test('the 58 labels and their minimum thicknesses', () => {
  assert.equal(HYBRID_LAYERS.K, 58);
  assert.equal(SURFACE_LAYERS, 12);
  assert.deepEqual(Array.from(LABELS.slice(0, 12)), Array.from({ length: 12 }, (_, k) => 1016.5 + 0.25 * k));
  assert.equal(LABELS[11], 1019.25);
  assert.equal(LABELS[57], 1027.5);
  const gaps = {};
  for (let k = 1; k < 58; k++) { const g = (LABELS[k] - LABELS[k - 1]).toFixed(3); gaps[g] = (gaps[g] ?? 0) + 1; }
  assert.deepEqual(gaps, { '0.250': 30, '0.500': 3, '0.125': 4, '0.100': 6, '0.050': 1, '0.030': 11, '0.120': 1, '0.400': 1 });
  assert.deepEqual(Array.from(LABELS.slice(12, 21)), [1019.5, 1020, 1020.5, 1021, 1021.25, 1021.5, 1021.625, 1021.75, 1021.875]);
  assert.deepEqual(Array.from(LABELS.slice(38, 44)), [1026.1, 1026.2, 1026.3, 1026.4, 1026.5, 1026.6]);
  assert.deepEqual(Array.from(LABELS.slice(44, 56)), [1026.65, 1026.68, 1026.71, 1026.74, 1026.77, 1026.8, 1026.83, 1026.86, 1026.89, 1026.92, 1026.95, 1026.98]);
  close(MINIMUM[0], 3, 1e-12, 'the top layer');
  close(MINIMUM[11], 6.99, 5e-3, 'the twelfth');
  for (let k = 1; k < 12; k++) close(MINIMUM[k] / MINIMUM[k - 1], 1.08, 1e-12, `growth at ${k}`);
  close(SURFACE_DEPTH, 56.93, 5e-3, 'the surface layers together');
  assert.deepEqual(Array.from(MINIMUM.slice(12, 21)), [8, 9, 10, 10, 10, 10, 10, 10, 10]);
  assert.ok(MINIMUM.slice(21, 44).every((m) => m === 10));
  close(ZGRID[21], 143.93, 5e-3, 'the warm pool’s base');
  close(ZGRID[44], 373.93, 5e-3, 'the thermocline’s base');
  close(MINIMUM[44], 40, 1e-9, 'the first deep layer');
  close(MINIMUM[57], 500, 1e-9, 'the last');
  close(ZGRID[58], 3019.0, 0.05, 'all the minima');
});

test('each label’s line is its isopycnal, and its water lies between the half-way densities to its neighbors', () => {
  for (let k = 0; k < 58; k++) {
    const line = isopycnalLine(LABELS[k], -2, 32);
    assert.equal(line.length, 137);
    for (const [s, t] of line) { assert.ok(Number.isFinite(s)); close(seawaterDensity(t + CELSIUS, s), LABELS[k], 1e-9, `label ${LABELS[k]}`); }
    const { below, above } = labelRange(k);
    if (k > 0) close(below, 0.5 * (LABELS[k - 1] + LABELS[k]), 1e-9, 'lower boundary');
    else assert.equal(below, -Infinity);
    if (k < 57) close(above, 0.5 * (LABELS[k] + LABELS[k + 1]), 1e-9, 'upper boundary');
    else assert.equal(above, Infinity);
  }
  for (let n = 0; n < 2000; n++) {
    const rho = 1015 + 14 * ((n * 0.6180339887) % 1);
    let best = 0;
    for (let k = 1; k < 58; k++) if (Math.abs(LABELS[k] - rho) < Math.abs(LABELS[best] - rho)) best = k;
    assert.equal(labelOf(rho), best, `ρ ${rho}`);
  }
  assert.equal(labelOf(1016), 0);
  assert.equal(labelOf(1030), 57);
  const salinityAt = (rho, t) => salinityForDensity(rho, t + CELSIUS);
  const k = labelOf(water(C(2.5), 34.92).rho), { below, above } = labelRange(k);
  const width = salinityAt(above, 2.5) - salinityAt(below, 2.5), spacing = salinityAt(LABELS[k + 1], 2.5) - salinityAt(LABELS[k], 2.5);
  close(above - below, 0.03, 1e-9, 'the deep layers’ water spans 0.03 kg/m³');
  assert.ok(width > 0.03 && width < 0.045, `${width} psu wide at 2.5 °C`);
  close(spacing, width, 2e-4, 'a label’s neighbor is a band’s width away');
});

test('the probe’s opening water and the plane’s corner', () => {
  const probe = water(C(25), 35);
  assert.equal(LABELS[labelOf(probe.rho)], 1022.5);
  assert.equal(labelOf(probe.rho), 23);
  const { below, above } = labelRange(23);
  close(below, 1022.375, 1e-9, 'below');
  close(above, 1022.625, 1e-9, 'above');
  assert.equal(LABELS[labelOf(water(C(-1.8), 34).rho)], 1026.5);
  let low = Infinity, high = -Infinity;
  for (let t = -2; t <= 32; t += 0.25) for (let s = 30; s <= 40; s += 0.1) { const r = seawaterDensity(t + CELSIUS, s); low = Math.min(low, r); high = Math.max(high, r); }
  close(low, 1016.41, 5e-3, 'the lightest water on the plane');
  assert.ok(high > 1031 && low < LABELS[0]);
  for (let k = 0; k < 58; k++) assert.ok(isopycnalLine(LABELS[k], -2, 32).some(([s]) => s >= 30 && s <= 40), `label ${LABELS[k]} crosses the plane`);
});

test('startSurface is initializeState’s sea', () => {
  const model = createModel(buildMesh(new Grid(16)), { ocean: false });
  const state = initializeState(model, { profile: new Float64Array(model.core.diagnostics.K).fill(300) });
  const { surfaceT, ice } = startSurface(model.mesh);
  let worst = 0;
  for (let i = 0; i < surfaceT.length; i++) worst = Math.max(worst, Math.abs(surfaceT[i] - state[3][i]));
  assert.ok(worst < 1e-9, `${worst} K`);
  assert.deepEqual(Array.from(ice), Array.from(state[6]));
});

test('the start depends on latitude alone', () => {
  const mesh = buildMesh(new Grid(16)), { surfaceT, ice } = startSurface(mesh), ocean = createOcean(mesh), n = mesh.nCells;
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

test('the start’s first regrid leaves every layer within half a meter of its target equatorward of 66 degrees', () => {
  const mesh = buildMesh(new Grid(16)), { surfaceT, ice } = startSurface(mesh), ocean = createOcean(mesh), n = mesh.nCells;
  ocean.initialize(surfaceT, ice);
  const before = Float64Array.from(ocean.h);
  ocean.regrid(Infinity);
  let tropics = 0, polar = 0;
  for (let i = 0; i < n; i++) for (let k = 0; k < ocean.layers; k++) {
    const change = Math.abs(ocean.h[k * n + i] - before[k * n + i]);
    if (Math.abs(mesh.latCell[i]) * 180 / Math.PI < 66) tropics = Math.max(tropics, change); else polar = Math.max(polar, change);
  }
  assert.ok(tropics < 0.5, `${tropics} m`);
  assert.ok(polar > 1 && polar < 40, `${polar} m`);
});

const section = layersByLatitude(16), K = section.layers;
const at = (lat) => nearestLatitude(section, lat);
const cell = (j, k) => ({ top: section.interface[j * (K + 1) + k], base: section.interface[j * (K + 1) + k + 1], t: section.temperature[j * K + k] - CELSIUS, s: section.salinity[j * K + k], rho: section.density[j * K + k], mode: modeAt(section, j, k), weight: section.weight[j * K + k] });

test('the layers by latitude are a stack from the surface to the floor in every column', () => {
  assert.equal(K, 58);
  assert.equal(section.latitude.length, 289);
  close(section.latitude[0], -90, 1e-9, 'south pole');
  close(section.latitude.at(-1), 90, 1e-9, 'north pole');
  for (let j = 1; j < section.latitude.length; j++) assert.ok(section.latitude[j] > section.latitude[j - 1]);
  for (let j = 0; j < section.latitude.length; j++) {
    assert.equal(section.interface[j * (K + 1)], 0);
    close(section.floor[j], 4000, 2, `floor at ${section.latitude[j]}`);
    for (let k = 0; k < K; k++) {
      const h = thickness(section, j, k);
      assert.ok(h >= H_MIN - 1e-9, `layer ${k} at ${section.latitude[j]} is ${h} m`);
      if (modeAt(section, j, k) !== 'empty') assert.ok(h >= MINIMUM[k] - 1e-6, `layer ${k} thinner than its minimum at ${section.latitude[j]}`);
      if (k > 0) assert.ok(section.density[j * K + k] >= section.density[j * K + k - 1] - 1e-9, `density falls at layer ${k}`);
      assert.ok(section.interface[j * (K + 1) + k + 1] >= section.interface[j * (K + 1) + k]);
    }
  }
});

test('the twelve surface layers sit at their minimum thicknesses in every column', () => {
  for (let j = 0; j < section.latitude.length; j++) {
    for (let k = 0; k < SURFACE_LAYERS; k++) { close(thickness(section, j, k), MINIMUM[k], 1e-6, `layer ${k}`); assert.equal(modeAt(section, j, k), 'fixed'); close(section.weight[j * K + k], 0, 1e-9, 'weight'); }
    close(section.interface[j * (K + 1) + SURFACE_LAYERS], SURFACE_DEPTH, 1e-6, 'the base of the surface layers');
  }
  const eq = at(0);
  close(cell(eq, 0).rho, 1021.032, 5e-4, 'the lightest water of the start');
  close(cell(eq, 11).rho, 1021.861, 5e-4, 'the twelfth layer at the equator');
  close(cell(eq, 0).t, 27.93, 5e-3, 'the equator’s surface layer');
  close(cell(eq, 11).t, 25.46, 5e-3, 'the equator’s twelfth layer');
  let lightest = Infinity;
  for (const rho of section.density) lightest = Math.min(lightest, rho);
  close(lightest, 1021.032, 5e-4, 'nothing lighter anywhere');
  assert.ok(LABELS[SURFACE_LAYERS - 1] < lightest - 1.5);
  assert.equal(LABELS.filter((r) => r < lightest).length, 16);
});

test('the weight and the mode of every layer agree', () => {
  for (let j = 0; j < section.latitude.length; j++) for (let k = 0; k < K; k++) {
    const { weight, mode } = cell(j, k), h = thickness(section, j, k);
    close(weight, Math.min(1, Math.max(0, (h - MINIMUM[k]) / 5)), 1e-12, 'weight');
    if (mode === 'fixed') close(weight, 0, 1e-6, 'fixed');
    if (mode === 'isopycnal') close(weight, 1, 1e-6, 'isopycnal');
    if (mode === 'between') assert.ok(weight > 0 && weight < 1);
    if (mode === 'empty') assert.ok(h <= H_MIN + 1e-6 && weight === 0);
  }
});

test('the nine densest layers hold no water in this start and lie on the floor', () => {
  assert.equal(emptyEverywhere(section), 9);
  for (let j = 0; j < section.latitude.length; j++) {
    for (let k = 49; k < 58; k++) assert.equal(modeAt(section, j, k), 'empty');
    assert.notEqual(modeAt(section, j, 48), 'empty');
    close(section.floor[j] - section.interface[j * (K + 1) + 49], 4.5, 1e-6, 'nine sheets of 0.5 m');
  }
  const eq = at(0), deep = cell(eq, 48);
  assert.equal(LABELS[48], 1026.77);
  close(deep.top, 2519.0, 0.05, 'where the deep water begins');
  close(deep.base - deep.top, 1477.7, 0.05, 'the deep water’s thickness');
  close(deep.rho, 1026.767, 5e-4, 'the deep water’s density');
  close(deep.t, 2.07, 5e-3, 'the deep water’s temperature');
  close(deep.s, 34.70, 5e-3, 'the deep water’s salinity');
  assert.equal(LABELS[labelOf(deep.rho)], 1026.77);
  assert.ok(deep.rho < labelRange(49).below);
});

test('the layer counts of the prose columns', () => {
  const walk = [[0, 28, 21], [30.3, 37, 12], [60.1, 43, 6], [79.9, 47, 2]];
  for (const [lat, fixed, isopycnal] of walk) {
    const j = at(lat), counts = modeCounts(section, j);
    close(section.latitude[j], lat, 0.05, 'latitude');
    assert.deepEqual(counts, { fixed, between: 0, isopycnal, empty: 9 }, `${lat}°`);
  }
  const eq = at(0);
  assert.equal(modeAt(section, eq, 27), 'fixed');
  assert.equal(modeAt(section, eq, 28), 'isopycnal');
  close(cell(eq, 28).top, 213.9, 0.05, 'the top of the equator’s first isopycnal layer');
  close(thickness(section, eq, 28), 21.8, 0.05, 'its thickness');
  assert.equal(LABELS[28], 1023.75);
  close(ZGRID[28], 213.93, 5e-3, 'the z grid at that layer');
  close(cell(eq, 28).top, ZGRID[28], 1e-6, 'its top is the z grid’s depth, as the layers above are at their minima');
});

test('the probe opens on the equator in the 1026.2 layer, isopycnal', () => {
  const eq = at(0), layer = layerAt(section, eq, 870);
  assert.equal(layer.k, 39);
  assert.equal(LABELS[39], 1026.2);
  assert.equal(layer.mode, 'isopycnal');
  close(layer.weight, 1, 1e-9, 'weight');
  close(layer.top, 836.8, 0.05, 'top');
  close(layer.base, 912.1, 0.05, 'base');
  close(layer.base - layer.top, 75.2, 0.05, 'thickness');
  const c = cell(eq, 39);
  close(c.t, 6.84, 5e-3, 'temperature');
  close(c.s, 34.65, 5e-3, 'salinity');
  close(c.rho, 1026.2, 5e-3, 'density');
  assert.equal(layerAt(section, eq, 4100), null);
  assert.equal(layerAt(section, eq, 0).k, 0);
  assert.equal(layerAt(section, eq, 3997).mode, 'empty');
});

test('following the 1026.2 layer poleward: it rises, reaches the top of its z grid and turns fixed-depth', () => {
  const k = 39, track = (lat) => cell(at(lat), k), shallowest = ZGRID[k];
  close(shallowest, 323.93, 5e-3, 'the shallowest the layer can lie');
  const north = track(30.3);
  close(north.top, 377.0, 0.05, 'top at 30.3°N');
  close(north.base, 429.8, 0.05, 'base at 30.3°N');
  close(cell(at(0), k).top - north.top, 460, 0.5, 'how far it has risen by 30.3°N');
  close(north.t - cell(at(0), k).t, 2.8, 0.05, 'warmer than at the equator');
  close(north.s - cell(at(0), k).s, 0.54, 5e-3, 'saltier than at the equator');
  close(north.t, 9.63, 5e-3, 'temperature at 30.3°N');
  close(north.s, 35.18, 5e-3, 'salinity at 30.3°N');
  assert.equal(north.mode, 'isopycnal');
  close(north.rho, 1026.2, 5e-3, 'density at 30.3°N');
  let previous = Infinity, first = null, between = [];
  for (let j = at(0); j < section.latitude.length; j++) {
    const c = cell(j, k);
    assert.ok(c.top <= previous + 1e-6, `the layer sinks at ${section.latitude[j]}`);
    assert.ok(c.top >= shallowest - 1e-6);
    previous = c.top;
    if (c.mode === 'between') between.push(section.latitude[j]);
    if (c.mode === 'fixed' && first === null) first = section.latitude[j];
    if (first !== null) { assert.equal(c.mode, 'fixed'); close(c.base - c.top, 10, 1e-6, 'held at its minimum'); close(c.top, shallowest, 1e-6, 'held at the top of its z grid'); }
  }
  close(first, 43.9, 0.05, 'where it turns fixed-depth');
  let reach = null;
  for (let j = at(0); j < section.latitude.length && reach === null; j++) if (cell(j, k).top <= shallowest + 1e-6) reach = j;
  close(section.latitude[reach], 37.0, 0.05, 'where it reaches the top of its z grid');
  close(thickness(section, reach, k), 50, 0.1, 'and how thick it is there');
  assert.deepEqual(between.map((lat) => lat.toFixed(1)), ['43.4', '43.6']);
  const top = cell(at(39.8), k);
  close(top.top, shallowest, 1e-6, 'it is at the top of its z grid by 39.8°N');
  close(top.base - top.top, 33.8, 0.05, 'and 33.8 m thick');
  const polar = track(60.1);
  assert.equal(polar.mode, 'fixed');
  close(polar.t, 2.55, 5e-3, 'temperature at 60.1°N');
  close(polar.s, 34.47, 5e-3, 'salinity at 60.1°N');
  close(polar.rho, 1026.544, 5e-4, 'density at 60.1°N');
  assert.equal(LABELS[labelOf(polar.rho)], 1026.5);
  close(polar.rho - LABELS[k], 0.34, 5e-3, 'denser than its label');
  const surface = cell(at(60.1), 0);
  close(surface.rho, 1026.167, 5e-4, 'the surface water at 60.1°N');
  assert.equal(LABELS[labelOf(surface.rho)], 1026.2);
  assert.equal(modeAt(section, at(60.1), 0), 'fixed');
});

test('at 79.9°N the surface water is denser than the labels of the 42 lightest layers', () => {
  const j = at(79.9), surface = cell(j, 0);
  close(surface.rho, 1026.476, 5e-4, 'surface density');
  close(surface.t, -1.8, 1e-9, 'under the ice');
  assert.equal(labelOf(surface.rho), 42);
  assert.equal(LABELS[42], 1026.5);
  assert.ok(labelRange(42).below < surface.rho);
  close(cell(j, 46).rho, 1026.731, 5e-4, 'density of the last fixed-depth layer');
  assert.equal(modeAt(section, j, 46), 'fixed');
  assert.equal(modeAt(section, j, 47), 'isopycnal');
  close(cell(j, 47).top, ZGRID[47], 1e-6, 'the first isopycnal layer starts at the end of the z grid above it');
  close(ZGRID[47], 521.5, 0.05, 'the z grid at that layer');
  close(cell(j, 46).rho - surface.rho, 0.255, 5e-3, 'the density rises by this much through 521 m');
});

test('the equator’s upper layers: twelve surface, nine warm-pool and seven thermocline layers are fixed-depth', () => {
  const eq = at(0);
  for (let k = 0; k < 28; k++) assert.equal(modeAt(section, eq, k), 'fixed');
  assert.equal(modeAt(section, eq, 28), 'isopycnal');
  assert.equal(SURFACE_LAYERS + 9 + 7, 28);
});

test('the section’s depth axis shows the twelve surface layers', () => {
  for (const z of [0, 3, 57, 100, 1000, 4000]) close(depthAt(depthFraction(z)), z, 1e-9, `round trip at ${z}`);
  assert.equal(depthFraction(0), 0);
  assert.equal(depthFraction(4000), 1);
  assert.equal(depthFraction(5000), 1);
  const plot = 370;
  close(depthFraction(SURFACE_DEPTH), 0.254, 1e-3, 'the surface layers take a quarter of the axis');
  for (let k = 0; k < SURFACE_LAYERS; k++) {
    const px = (depthFraction(ZGRID[k + 1]) - depthFraction(ZGRID[k])) * plot;
    assert.ok(px > 6.5, `surface layer ${k} is ${px} px`);
  }
  close((depthFraction(3) - depthFraction(0)) * plot, 9.8, 0.1, 'the top layer');
  close((depthFraction(ZGRID[12]) - depthFraction(ZGRID[11])) * plot, 6.6, 0.1, 'the twelfth');
});
