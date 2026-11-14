import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createGeography, readTopography, encodeSoilWaterMonthly, decodeSoilWaterMonthly, decodeSoilWater, soilWaterUrl, soilWaterMonthlyUrl, soilWaterMonthlyOnMesh, SOIL_SCALE, SOIL_NONE } from '../js/geography.module.js';
import { monthWeights, MID_MONTH_DAYS, YEAR } from '../js/monthCalendar.module.js';
import { fillMonths } from '../scripts/soilWater.mjs';

const monthly = decodeSoilWaterMonthly(readFileSync(soilWaterMonthlyUrl));
const march = decodeSoilWater(readFileSync(soilWaterUrl));
const size = monthly.rows * monthly.cols;
const at = (m, lat, lon) => Math.round((lat + 89.75) / 0.5) * monthly.cols + Math.round((lon + 179.75) / 0.5) + m * size;
const fillAt = (m, lat, lon) => monthly.fill[at(m, lat, lon)] / SOIL_SCALE;

test('a monthly soil-water file round-trips its grid and bytes and refuses a wrong tag, month count or size', () => {
  const fill = Uint8Array.from({ length: 12 * 6 * 3 }, (_, n) => (n % 7 === 0 ? SOIL_NONE : (n * 37) % (SOIL_SCALE + 1)));
  const bytes = new Uint8Array(encodeSoilWaterMonthly({ cols: 6, rows: 3, months: 12, fill })), back = decodeSoilWaterMonthly(bytes);
  assert.deepEqual([back.cols, back.rows, back.months, [...back.fill]], [6, 3, 12, [...fill]]);
  assert.throws(() => decodeSoilWaterMonthly(bytes.slice(0, -1)), /wrong size/);
  assert.throws(() => encodeSoilWaterMonthly({ cols: 6, rows: 3, months: 12, fill: fill.slice(1) }), /values for/);
  const eleven = new Uint8Array(encodeSoilWaterMonthly({ cols: 6, rows: 3, months: 11, fill: fill.slice(0, 11 * 18) }));
  assert.throws(() => decodeSoilWaterMonthly(eleven), /11 months/);
  bytes[0] = 0;
  assert.throws(() => decodeSoilWaterMonthly(bytes), /not a monthly soil-water file/);
  const march = new Uint8Array(readFileSync(soilWaterUrl));
  assert.throws(() => decodeSoilWaterMonthly(march), /not a monthly soil-water file/);
});

test('data/soilwater_monthly.bin holds twelve fills with a dry Sahara, a wet Amazon, the Indian monsoon cycle and no Antarctica', () => {
  assert.deepEqual([monthly.cols, monthly.rows, monthly.months, readFileSync(soilWaterMonthlyUrl).length], [720, 360, 12, 16 + 12 * 720 * 360]);
  assert.ok(monthly.fill.every((v) => v <= SOIL_SCALE || v === SOIL_NONE), 'every byte is a fill or none');
  for (let m = 0; m < 12; m++) {
    assert.ok(fillAt(m, 23, 10) < 0.1, `Sahara month ${m + 1}: ${fillAt(m, 23, 10)}`);
    assert.ok(fillAt(m, -3, -60) >= 0.4, `Amazon month ${m + 1}: ${fillAt(m, -3, -60)}`);
    assert.equal(monthly.fill[at(m, -80, 0)], SOIL_NONE, `Antarctica month ${m + 1}`);
    assert.equal(monthly.fill[at(m, 0, -140)], SOIL_NONE, `equatorial Pacific month ${m + 1}`);
  }
  assert.ok(fillAt(4, 22, 78) < fillAt(7, 22, 78) - 0.3, `India May ${fillAt(4, 22, 78)}, August ${fillAt(7, 22, 78)}`);
});

test('every point has all twelve months or none, so the fill leaves no month missing where another holds a value', () => {
  for (let n = 0; n < size; n++) {
    let held = 0;
    for (let m = 0; m < 12; m++) if (monthly.fill[m * size + n] !== SOIL_NONE) held++;
    assert.ok(held === 0 || held === 12, `point ${n} holds ${held} months`);
  }
});

test('interpolating the months to day 0 reproduces data/soilwater_march.bin within 1/250 where both are valid', () => {
  const [m0, m1, f] = monthWeights(0);
  let compared = 0;
  for (let n = 0; n < size; n++) {
    const a = march.fill[n], b = monthly.fill[m0 * size + n], c = monthly.fill[m1 * size + n];
    if (a === SOIL_NONE || b === SOIL_NONE || c === SOIL_NONE) continue;
    compared++;
    assert.ok(Math.abs(((1 - f) * b + f * c) / SOIL_SCALE - a / SOIL_SCALE) <= 1 / 250 + 1e-12, `point ${n}: ${a} against ${b}, ${c}`);
  }
  assert.ok(compared > 60000, `${compared} points compared`);
});

test('fillMonths takes a missing month from the nearest valid months, cyclically, and leaves a cell with none alone', () => {
  const v = (...x) => x.map((e) => (e === 0 ? null : e));
  assert.deepEqual(fillMonths(v(1, 0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)), [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]);
  const wrapped = fillMonths(v(0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0));
  assert.ok(Math.abs(wrapped[0] - 5) < 1e-12 && Math.abs(wrapped[11] - 8) < 1e-12, `${wrapped[0]}, ${wrapped[11]}`);
  assert.deepEqual(fillMonths(v(0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0)), Array(12).fill(6));
  assert.deepEqual(fillMonths(Array(12).fill(null)), Array(12).fill(null));
});

test('the month weights bracket a model time by its mid-month days, cyclic across December and January', () => {
  const near = (day, [a, b, f]) => {
    const got = monthWeights(day);
    assert.deepEqual(got.slice(0, 2), [a, b], `day ${day}`);
    assert.ok(Math.abs(got[2] - f) < 1e-12, `day ${day}: ${got[2]} against ${f}`);
  };
  near(0, [2, 3, 5 / 31]);
  near(15, [2, 3, 20 / 31]);
  near(45, [3, 4, 19 / 30]);
  near(79, [4, 5, 23 / 31]);
  near(300, [11, 0, 30 / 31]);
  near(270, [11, 0, 0]);
  near(285, [11, 0, 15 / 31]);
  near(286, [11, 0, 16 / 31]);
  near(301, [0, 1, 0]);
  near(-5, [2, 3, 0]);
  near(-365, [2, 3, 5 / 31]);
  near(365 * 3 + 45, [3, 4, 19 / 30]);
  assert.deepEqual(MID_MONTH_DAYS, [-64, -34, -5, 26, 56, 87, 117, 148, 179, 209, 240, 270]);
  MID_MONTH_DAYS.forEach((day, m) => near(day, [m, (m + 1) % 12, 0]));
  assert.equal(YEAR, 365);
});

test('the monthly fields on the mesh hold a value on every open land cell at N = 16 and N = 64, NaN on the sea, and createGeography reads the file only when asked', () => {
  const topography = readTopography('etopo1');
  for (const N of [16, 64]) {
    const mesh = buildMesh(new Grid(N)), plain = createGeography(mesh, topography);
    assert.equal(plain.wetnessMonthly, undefined, 'a plain geography holds no monthly field');
    const geography = createGeography(mesh, topography, { soilWaterMonthly: true });
    assert.equal(geography.wetnessMonthly.length, 12);
    assert.deepEqual(geography.wetness, plain.wetness);
    for (let m = 0; m < 12; m++) {
      const w = geography.wetnessMonthly[m];
      assert.ok(w instanceof Float64Array && w.length === mesh.nCells);
      for (let i = 0; i < mesh.nCells; i++) {
        if (geography.land[i] && !geography.iceSheet[i]) assert.ok(w[i] >= 0 && w[i] <= 1, `N = ${N}, month ${m + 1}, cell ${i}: ${w[i]}`);
        if (!geography.land[i]) assert.ok(Number.isNaN(w[i]), `N = ${N}, month ${m + 1}, sea cell ${i}`);
      }
    }
    const again = soilWaterMonthlyOnMesh(monthly, mesh, geography.land);
    assert.deepEqual(again, geography.wetnessMonthly);
  }
});
