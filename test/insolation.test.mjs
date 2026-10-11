import { test } from 'node:test';
import assert from 'node:assert/strict';
import { insolationMap, dailyMean, daylightHours, subsolarLatitude, noonElevation, seasonPhrase, DAYS } from '../explain/insolation.module.js';
import { SOLAR_CONSTANT, AXIAL_TILT } from '../js/physics/radiation.module.js';
import { insolationCycle } from '../js/physics/land.module.js';

const D = Math.PI / 180;
const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const map = insolationMap();
const at = (lat, day) => map.daily[(lat + 90) * DAYS + day];
const tilt = AXIAL_TILT / D;

function analytic(lat, day) {
  const dec = subsolarLatitude(day) * D, phi = lat * D, c = -Math.tan(phi) * Math.tan(dec);
  const h = c >= 1 ? 0 : c <= -1 ? Math.PI : Math.acos(c);
  return SOLAR_CONSTANT / Math.PI * (h * Math.sin(phi) * Math.sin(dec) + Math.cos(phi) * Math.cos(dec) * Math.sin(h));
}

test('the map has a row for every degree of latitude and a column for every day', () => {
  assert.equal(map.rows, 181);
  assert.equal(map.daily.length, 181 * DAYS);
  close(map.latitude[0], -90, 0, 'first row'); close(map.latitude[90], 0, 0, 'equator row'); close(map.latitude[180], 90, 0, 'last row');
  close(at(37, 150), dailyMean(37, 150), 1e-4, 'a cell against dailyMean');
});

test('the sun stands over the equator at the equinoxes and over 23.44° at the solstices', () => {
  close(subsolarLatitude(0), 0, 1e-12, 'March equinox');
  close(subsolarLatitude(182.5), 0, 1e-9, 'September equinox');
  close(subsolarLatitude(91.25), tilt, 1e-9, 'June solstice');
  close(subsolarLatitude(273.75), -tilt, 1e-9, 'December solstice');
});

test('the equator gets S₀/π at the equinoxes and S₀ cos ε/π at the solstices', () => {
  close(at(0, 0), SOLAR_CONSTANT / Math.PI, 0.02, 'March equinox');
  close(at(0, 183), SOLAR_CONSTANT / Math.PI, 0.02, 'September equinox');
  close(at(0, 91), SOLAR_CONSTANT / Math.PI * Math.cos(AXIAL_TILT), 0.02, 'June solstice');
  close(at(0, 274), SOLAR_CONSTANT / Math.PI * Math.cos(AXIAL_TILT), 0.02, 'December solstice');
  for (let d = 0; d < DAYS; d++) assert.ok(at(0, d) <= SOLAR_CONSTANT / Math.PI + 0.02, `equator on day ${d}`);
});

test('45°N gets 499.7 W/m² on the June solstice and 116.6 on the December one', () => {
  close(at(45, 91), 499.7, 0.05, 'June');
  close(at(45, 274), 116.6, 0.05, 'December');
  close(at(-45, 274), at(45, 91), 0.01, 'the southern mirror');
});

test('the summer pole gets S₀ sin ε, more than any other place on any day', () => {
  close(at(90, 91), SOLAR_CONSTANT * Math.sin(AXIAL_TILT), 0.02, 'North Pole, June');
  close(at(-90, 274), SOLAR_CONSTANT * Math.sin(AXIAL_TILT), 0.02, 'South Pole, December');
  const most = map.daily.reduce((a, b) => Math.max(a, b), 0);
  close(most, at(90, 91), 0.01, 'the largest daily mean');
  for (let lat = 9; lat <= 90; lat++) assert.ok(at(lat, 91) > at(0, 0), `${lat}°N on the June solstice beats the equator's best`);
});

test('every daily mean matches the analytic form at the model’s declination', () => {
  for (let i = 0; i < map.rows; i++) if (Math.abs(map.latitude[i]) <= 80) for (let d = 0; d < DAYS; d++) close(map.daily[i * DAYS + d], analytic(map.latitude[i], d), 0.1, `${map.latitude[i]}° on day ${d}`);
});

test('the annual means match insolationCycle and average to S₀/4 over the sphere', () => {
  for (let i = 0; i < map.rows; i++) close(map.annual[i], insolationCycle(map.latitude[i] * D, 365).mean, 0.01, `${map.latitude[i]}°`);
  close(map.annual[90], 415.6, 0.05, 'equator'); close(map.annual[135], 307.1, 0.05, '45°N'); close(map.annual[180], 174.1, 0.05, 'North Pole');
  let sum = 0, weight = 0;
  for (let i = 0; i < map.rows; i++) { const w = Math.cos(map.latitude[i] * D) * (i % (map.rows - 1) ? 1 : 0.5); sum += map.annual[i] * w; weight += w; }
  close(sum / weight, SOLAR_CONSTANT / 4, 0.02, 'area mean');
});

test('the winter polar night starts at 90° − 23.44°', () => {
  for (let lat = 67; lat <= 90; lat++) { assert.equal(at(lat, 274), 0, `${lat}°N in December`); assert.equal(at(-lat, 91), 0, `${lat}°S in June`); }
  assert.ok(at(66, 274) > 0, '66°N in December');
  let dark = 0;
  for (let d = 0; d < DAYS; d++) if (at(90, d) === 0) dark++;
  close(dark, DAYS / 2, 1, 'days of night at the North Pole');
});

test('daylight is 12 hours at the equator and follows the hour angle elsewhere', () => {
  for (const d of [0, 50, 91, 183, 274]) close(daylightHours(0, d), 12, 1 / 60, `equator on day ${d}`);
  close(daylightHours(90, 91), 24, 0, 'pole in June'); close(daylightHours(90, 274), 0, 0, 'pole in December');
  for (const lat of [-60, -30, 10, 45, 60, 66]) for (const d of [30, 91, 150, 274, 330]) {
    const c = -Math.tan(lat * D) * Math.tan(subsolarLatitude(d) * D), hours = 24 / Math.PI * (c >= 1 ? 0 : c <= -1 ? Math.PI : Math.acos(c));
    close(daylightHours(lat, d), hours, 2 / 60, `${lat}° on day ${d}`);
  }
  close(noonElevation(0, 91), 90 - tilt, 0.001, 'equator at noon in June');
  close(noonElevation(90, 91), tilt, 0.001, 'pole in June');
});

test('days are named from the nearest equinox or solstice', () => {
  assert.equal(seasonPhrase(0), 'the March equinox');
  assert.equal(seasonPhrase(91), 'the June solstice');
  assert.equal(seasonPhrase(120), '29 days after the June solstice');
  assert.equal(seasonPhrase(182), '1 day before the September equinox');
  assert.equal(seasonPhrase(274), 'the December solstice');
  assert.equal(seasonPhrase(364), '1 day before the March equinox');
});
