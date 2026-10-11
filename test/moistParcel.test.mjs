import { test } from 'node:test';
import assert from 'node:assert/strict';
import { saturationHumidity } from '../js/physics/moist.module.js';
import { heightOf, KAPPA, CP, P0 } from '../explain/physics.module.js';
import { moistLapse, dewPoint, liftParcel, parcelAt, moistAdjustColumn, DRY_LAPSE } from '../explain/moistParcel.module.js';

const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const C0 = 273.15;
const at = (parcel, hPa) => {
  const j = parcel.p.findIndex((p) => Math.abs(p - hPa * 100) < 1e-6);
  assert.ok(j >= 0, `a sample at ${hPa} hPa`);
  return Object.fromEntries(['p', 'z', 'T', 'q', 'condensate', 'theta', 'dry', 'dew'].map((key) => [key, parcel[key][j]]));
};

test('a parcel at 15 °C and 1000 hPa: its vapor, dew point and cloud base at 50, 70 and 90 percent', () => {
  for (const [rh, q, dew, lcl, height, temperature] of [[0.5, 5.33, 4.7, 856, 1289, 2.5], [0.7, 7.47, 9.6, 922, 678, 8.4], [0.9, 9.60, 13.4, 976, 204, 13.0]]) {
    const parcel = liftParcel({ T0: 288.15, rh }, 200e2);
    close(parcel.q0 * 1000, q, 0.005, `q at ${rh}`);
    close(parcel.dew0 - C0, dew, 0.05, `dew point at ${rh}`);
    close(parcel.lcl.pressure / 100, lcl, 0.5, `cloud base pressure at ${rh}`);
    close(parcel.lcl.height, height, 1, `cloud base height at ${rh}`);
    close(parcel.lcl.temperature - C0, temperature, 0.05, `cloud base temperature at ${rh}`);
    close(parcel.lcl.height, heightOf(parcel.lcl.pressure), 0, 'heights are Part 1\'s standard atmosphere');
  }
});

test('at 70 percent the cloudy parcel arrives at 500 hPa 14.0 K warmer than a dry one, with 6.25 g/kg condensed', () => {
  const parcel = liftParcel({ T0: 288.15, rh: 0.7 }, 200e2);
  const mid = at(parcel, 700), half = at(parcel, 500), high = at(parcel, 200);
  close(mid.T - C0, -4.5, 0.05, '700 hPa');
  close(mid.dry - C0, -12.95, 0.005, '700 hPa, dry');
  close(half.T - C0, -22.9, 0.05, '500 hPa');
  close(half.dry - C0, -36.8, 0.05, '500 hPa, dry');
  close(half.T - half.dry, 13.97, 0.005, 'warmer than a dry parcel');
  close(half.condensate * 1000, 6.25, 0.005, 'condensed by 500 hPa, g/kg');
  close(half.condensate * 2.5e6 / CP, 15.6, 0.05, 'its latent heat, K');
  close(half.theta, 305.2, 0.05, 'θ at 500 hPa');
  close(high.theta, 308.9, 0.05, 'θ at 200 hPa');
  close(half.z, 5478, 1, '500 hPa height');
  close((parcel.lcl.temperature - mid.T) / (mid.z - parcel.lcl.height) * 1000, 5.8, 0.05, 'mean lapse rate from cloud base to 700 hPa, K/km');
  close((parcel.lcl.temperature - half.T) / (half.z - parcel.lcl.height) * 1000, 6.5, 0.05, 'mean lapse rate from cloud base to 500 hPa, K/km');
});

test('a tropical parcel at 300 K and 80 percent: 17.8 g/kg, cloud base 948 hPa, −1.7 °C at 500 hPa, θ 352.2 K at 200 hPa', () => {
  const parcel = liftParcel({ T0: 300, rh: 0.8 }, 200e2);
  close(parcel.q0 * 1000, 17.8, 0.05, 'q');
  close(parcel.lcl.pressure / 100, 948, 0.5, 'cloud base');
  close(parcel.lcl.height, 452, 1, 'cloud base height');
  close(at(parcel, 500).T - C0, -1.7, 0.05, '500 hPa');
  close(at(parcel, 500).dry - C0, -27.1, 0.05, '500 hPa, dry');
  close(at(parcel, 200).theta, 352.2, 0.05, 'θ at 200 hPa');
});

test('the moist adiabatic lapse rate: 3.71 K/km in warm air near the ground, 9.21 in cold air aloft, against 9.77 dry', () => {
  for (const [T, p, rate] of [[300, 1000e2, 3.71], [288, 1000e2, 4.74], [270, 700e2, 6.19], [250, 500e2, 8.10], [230, 300e2, 9.21]]) close(moistLapse(T, p) * 1000, rate, 0.005, `${T} K, ${p / 100} hPa`);
  close(DRY_LAPSE * 1000, 9.77, 0.005, 'dry');
  for (let T = 200; T <= 310; T += 10) assert.ok(moistLapse(T, 800e2) < DRY_LAPSE, `slower than dry at ${T} K`);
});

test('the dew point is Bolton\'s saturation inverted', () => {
  for (const T of [240, 263.15, 273.15, 288.15, 303.15]) for (const p of [300e2, 700e2, 1000e2]) close(dewPoint(saturationHumidity(T, p), p), T, 1e-9, `${T} K at ${p / 100} hPa`);
  assert.ok(Number.isNaN(dewPoint(0, 1e5)), 'perfectly dry air has none');
});

test('below cloud base the parcel is dry adiabatic and keeps its vapor; above it it is saturated, and its water and θ only grow', () => {
  for (const [T0, rh] of [[288.15, 0.5], [288.15, 0.9], [300, 0.8], [273.15, 0.6]]) {
    const parcel = liftParcel({ T0, rh });
    const { lcl } = parcel;
    for (let j = 0; j < parcel.p.length; j++) {
      close(parcel.q[j] + parcel.condensate[j], parcel.q0, 1e-15, 'water is kept or condensed');
      if (parcel.p[j] > lcl.pressure + 10) {
        close(parcel.T[j], parcel.dry[j], 1e-6, `dry below cloud base at ${parcel.p[j]} Pa`);
        assert.equal(parcel.condensate[j], 0);
        close(parcel.theta[j], T0 * (P0 / 1e5) ** KAPPA, 1e-6, 'θ is kept');
      } else if (parcel.p[j] < lcl.pressure - 10) {
        close(parcel.q[j], saturationHumidity(parcel.T[j], parcel.p[j]), 1e-15, `saturated at ${parcel.p[j]} Pa`);
        close(parcel.dew[j], parcel.T[j], 0, 'cloud is at its dew point');
        assert.ok(parcel.theta[j] > parcel.theta[j - 1] - 1e-6 && parcel.condensate[j] >= parcel.condensate[j - 1], `θ and condensate grow at ${parcel.p[j]} Pa`);
      }
    }
    close(dewPoint(parcel.q0, lcl.pressure), lcl.temperature, 0.1, `cloud base temperature is its dew point, ${T0} K and ${rh}`);
  }
  const dry = liftParcel({ rh: 0 });
  assert.equal(dry.lcl, null);
  assert.ok(dry.condensate.every((c) => c === 0) && dry.T.every((T, j) => Math.abs(T - dry.dry[j]) < 1e-5), 'a parcel without vapor stays dry');
});

test('parcelAt reads the path at a height between samples', () => {
  const parcel = liftParcel({ T0: 288.15, rh: 0.7 });
  const j = 300, mid = 0.5 * (parcel.z[j] + parcel.z[j + 1]);
  close(parcelAt(parcel, parcel.z[j]).T, parcel.T[j], 1e-12, 'at a sample');
  close(parcelAt(parcel, mid).T, 0.5 * (parcel.T[j] + parcel.T[j + 1]), 1e-9, 'between samples');
  assert.equal(parcelAt(parcel, 0).cloudy, false);
  assert.equal(parcelAt(parcel, 5000).cloudy, true);
  close(parcelAt(parcel, -100).p, parcel.p0, 0, 'below the start');
});

function testColumn() {
  const K = 30, pressure = Float64Array.from({ length: K }, (_, k) => 1e5 * (k + 0.5) / K), mass = new Float64Array(K).fill(1e5 / K / 9.806);
  const T = Float64Array.from(pressure, (p) => Math.max(205, 300 * (p / 1e5) ** (KAPPA * 1.05)));
  return { K, pressure, mass, T };
}

test('moistAdjustColumn conserves c_p T plus the heat it is given and puts the unstable layers on one adiabat', () => {
  for (const [moist, heat] of [[true, 0], [true, 60 * 10800], [false, 0], [false, 2e6]]) {
    const { K, pressure, mass, T } = testColumn(), start = Float64Array.from(T), out = new Float64Array(K);
    const top = moistAdjustColumn(T, pressure, mass, { rh: 0.8, heat, moist, out });
    let energy = 0;
    for (let k = 0; k < K; k++) energy += CP * mass[k] * (T[k] - start[k]);
    close(energy, heat, 1e-3, `energy, moist ${moist}, heat ${heat}`);
    assert.ok(top < K - 3, `the adjustment reaches up, moist ${moist}`);
    for (let k = 0; k < top; k++) assert.equal(T[k], start[k], 'layers above are untouched');
    assert.ok(T[top - 1] >= out[top - 1], 'the layer above the top is no colder than the parcel would be there');
    if (!moist) for (let k = top; k < K; k++) close(T[k] * (1e5 / pressure[k]) ** KAPPA, T[K - 1] * (1e5 / pressure[K - 1]) ** KAPPA, 2e-3, `dry adiabat at layer ${k}`);
  }
});

test('moistAdjustColumn: a moist adiabat cools more slowly than a dry one, and a stable column only takes the heat in its lowest layer', () => {
  const dry = testColumn(), moist = testColumn();
  moistAdjustColumn(dry.T, dry.pressure, dry.mass, { moist: false });
  moistAdjustColumn(moist.T, moist.pressure, moist.mass, { rh: 0.8 });
  const k = 15;
  assert.ok(moist.T[k] - moist.T[moist.K - 1] > dry.T[k] - dry.T[dry.K - 1], 'less cooling with height');
  const { K, pressure, mass } = testColumn(), T = Float64Array.from(pressure, (p) => 280 + 20 * (1 - p / 1e5)), start = Float64Array.from(T);
  const top = moistAdjustColumn(T, pressure, mass, { rh: 0.8, heat: 1e5 });
  assert.equal(top, K - 1);
  close(T[K - 1] - start[K - 1], 1e5 / (CP * mass[K - 1]), 1e-6, 'the heat warms the lowest layer');
  for (let j = 0; j < K - 1; j++) assert.equal(T[j], start[j]);
});

test('moistAdjustColumn starts one layer below the top it is given', () => {
  const { K, pressure, mass, T } = testColumn(), stable = Float64Array.from(pressure, (p) => 280 + 20 * (1 - p / 1e5));
  const top = moistAdjustColumn(Float64Array.from(stable), pressure, mass, { rh: 0.8, top: 10 });
  assert.equal(top, 11, 'a stable column keeps the layers below the old top on the adiabat');
  assert.ok(moistAdjustColumn(T, pressure, mass, { rh: 0.8, top: K - 1 }) < K - 3);
});
