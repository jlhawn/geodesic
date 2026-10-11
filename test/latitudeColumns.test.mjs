import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { meanSunlight, K } from '../explain/radiativeColumn.module.js';
import { readLatitudes, createLatitudes, areaWeights, areaMean, heldSuarezSurface, effectiveTemperature, temperatureAt, CALM_DAYS, STEPS_PER_DAY } from '../explain/latitudeColumns.module.js';

const DATA = readLatitudes(JSON.parse(readFileSync(new URL('../explain/data/latitudes.json', import.meta.url))));
const D = Math.PI / 180;
const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);

test('the ten columns sit at the N = 16 grid’s cells nearest 0°, 10°, … 90°, each under its annual-mean sunlight', () => {
  const mesh = buildMesh(new Grid(16)), cells = new Set(Array.from(mesh.latCell, (lat) => Math.abs(lat / D).toFixed(4)));
  const nearest = (target) => Math.min(...Array.from(mesh.latCell, (lat) => Math.abs(Math.abs(lat / D) - target)));
  assert.equal(DATA.latitudes.length, 10);
  close(DATA.latitudes[0], 0, 1e-9, 'first column');
  close(DATA.latitudes[9], 90, 1e-9, 'last column');
  DATA.latitudes.forEach((lat, n) => {
    assert.ok(cells.has(lat.toFixed(4)), `${lat}° is a cell latitude`);
    close(Math.abs(lat - 10 * n), nearest(10 * n), 1e-6, `${lat}° is the cell latitude nearest ${10 * n}°`);
    const sun = meanSunlight(lat);
    close(DATA.insolation[n], sun.insolation, 1e-5, `insolation at ${lat}°`);
    close(DATA.mu[n], sun.mu, 1e-7, `μ at ${lat}°`);
  });
  close(DATA.insolation[0], 415.6, 0.05, 'equator');
  close(DATA.insolation[9], 174.1, 0.05, 'pole');
});

test('the area weights integrate a profile linear in latitude exactly against cos φ', () => {
  const weights = areaWeights(DATA.latitudes);
  close(weights.reduce((a, b) => a + b, 0), 1, 1e-12, 'weights sum to one');
  close(areaMean(weights, DATA.latitudes.map(() => 7)), 7, 1e-12, 'a constant');
  close(areaMean(weights, DATA.latitudes.map((lat) => lat * D)), Math.PI / 2 - 1, 1e-12, '∫ φ cos φ dφ');
  close(areaMean(weights, DATA.insolation), 340.5, 0.5, 'the columns\' sunlight averages to S₀/4');
});

test('temperatureAt interpolates linearly in log pressure', () => {
  const pressure = Float64Array.from({ length: K }, (_, k) => 30 * 1.2 ** k * 100), T = Float64Array.from(pressure, (p) => 200 + 10 * Math.log(p));
  close(temperatureAt({ pressure, T }, 50000), 200 + 10 * Math.log(50000), 1e-9, '500 hPa');
});

test('every latitude balances its own sunlight with its own infrared', () => {
  for (const column of createLatitudes(DATA).columns) {
    close(column.budget.imbalance, 0, 0.01, 'absorbed sunlight − OLR');
    close(column.budget.groundNet, 0, 0.01, 'net flux into the ground');
  }
});

test('the dry planet, every latitude for itself: 298.4 K at the equator, 234.7 K at the pole', () => {
  const { rows, mean, contrast, targetContrast } = createLatitudes(DATA).summary();
  for (const [n, Ts] of [[0, 298.4], [3, 288.2], [6, 256.9], [9, 234.7]]) close(rows[n].Ts, Ts, 0.1, `Ts at ${rows[n].latitude}°`);
  close(rows[0].T500, 252.0, 0.1, '500 hPa at the equator');
  close(rows[0].absorbed, 365.8, 0.05, 'absorbed at the equator');
  close(rows[9].absorbed, 145.3, 0.1, 'absorbed at the pole');
  close(heldSuarezSurface(0), 315, 1e-12, 'Held–Suarez at the equator');
  close(heldSuarezSurface(90), 255, 1e-12, 'Held–Suarez at the pole');
  for (const r of rows) assert.ok(r.Ts < r.target - 9, `${r.latitude}° is more than 9 K colder than the Held–Suarez target`);
  close(contrast, 63.7, 0.1, 'equator to pole');
  close(targetContrast, 60, 1e-12, 'Held–Suarez equator to pole');
  close(mean.Ts, 281.4, 0.1, 'area-mean ground temperature');
  close(mean.absorbed, 297.3, 0.1, 'area-mean absorbed sunlight');
  close(mean.Te, effectiveTemperature(mean.absorbed), 1e-12, 'effective temperature');
  close(mean.Te, 269.1, 0.1, 'effective temperature');
});

test('the columns run until each holds its balance for ten days, and warm when CO₂ doubles', () => {
  const band = createLatitudes(DATA);
  assert.equal(band.settled, false);
  assert.equal(band.advance(CALM_DAYS * STEPS_PER_DAY), 0);
  assert.ok(band.settled);
  close(band.days, CALM_DAYS, 0, 'days run');
  const before = band.columns.map((column) => column.Ts);
  band.set({ carbonDioxide: 2 });
  assert.equal(band.settled, false);
  assert.equal(band.advance(50 * STEPS_PER_DAY), DATA.latitudes.length);
  band.columns.forEach((column, n) => assert.ok(column.Ts > before[n] + 0.1, `${DATA.latitudes[n]}° warms`));
});
