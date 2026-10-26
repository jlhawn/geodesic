import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createSigmaCore } from '../js/dynamics/sigmaCore.module.js';
import { createRadiation, sunDirection, AXIAL_TILT, DAY, YEAR } from '../js/physics/radiation.module.js';
import { SOLAR_CONSTANT, timeOf, subsolarPoint, sunlightAt, sunlightOnCells, areaMean } from '../explain/subsolar.module.js';

const D = Math.PI / 180;
const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const wrap = (lon) => ((lon / D + 540) % 360) - 180;
const unit = (lat, lon) => [Math.cos(lat) * Math.cos(lon), Math.cos(lat) * Math.sin(lon), Math.sin(lat)];

test('day 0 at noon at Greenwich is t = 0, the March equinox with the sun over 0° N, 0° E', () => {
  assert.equal(timeOf(0, 12), 0);
  for (const orbit of ['circular', 'kepler']) {
    const sun = sunDirection(0, undefined, orbit), { lat, lon } = subsolarPoint(sun);
    close(lat, 0, 1e-12, `${orbit} latitude`);
    close(lon, 0, 0, `${orbit} longitude`);
    close(sunlightAt(...unit(0, 0), sun), 1362, 1e-9, `${orbit} sunlight under the sun`);
    close(sunlightAt(...unit(45 * D, 0), sun), 1362 * Math.SQRT1_2, 1e-9, `${orbit} sunlight at 45° N on the Greenwich meridian`);
    close(sunlightAt(...unit(0, 90 * D), sun), 0, 1e-12, `${orbit} sunlight on the terminator`);
  }
});

test('the sun moves west 15° an hour and is over Greenwich at noon every day', () => {
  for (const day of [0, 45, 91, 200, 364]) for (let hour = 0; hour <= 24; hour += 0.25) {
    const { lon } = subsolarPoint(sunDirection(timeOf(day, hour)));
    close(wrap(lon + 15 * (hour - 12) * D), 0, 1e-9, `day ${day} hour ${hour}`);
  }
});

test('on the circular orbit the subsolar latitude is 23.44°·sin(2πt/YEAR): 23.44° N at the June solstice, 23.44° S at the December one', () => {
  close(timeOf(91, 18), YEAR / 4, 0, 'June solstice time');
  close(subsolarPoint(sunDirection(YEAR / 4, undefined, 'circular')).lat / D, 23.44, 1e-12, 'June solstice');
  close(subsolarPoint(sunDirection(3 * YEAR / 4, undefined, 'circular')).lat / D, -23.44, 1e-12, 'December solstice');
  for (let n = 0; n < 200; n++) {
    const t = n * 1.83 * DAY + n * 977;
    close(subsolarPoint(sunDirection(t, undefined, 'circular')).lat, AXIAL_TILT * Math.sin(2 * Math.PI * (t % YEAR) / YEAR), 1e-12, `t = ${t}`);
  }
});

test('every cell gets exactly the sunlight the model\'s own radiation gives it', () => {
  const mesh = buildMesh(new Grid(16)), radiation = createRadiation(mesh, createSigmaCore(mesh), { orbit: 'circular' }), out = new Float64Array(mesh.nCells);
  for (const t of [0, 0.3 * DAY, 91.25 * DAY + 5000, 200.7 * DAY, 300 * DAY + 77]) {
    radiation.setTime(t);
    sunlightOnCells(mesh.xCell, sunDirection(t, undefined, 'circular'), out);
    for (let i = 0; i < mesh.nCells; i++) assert.equal(out[i], radiation.insolation(i), `cell ${i} at t = ${t}`);
  }
});

test('on the N = 32 grid the mean over the globe is a quarter of 1,362 W/m² at every hour of every day, with half the area lit', () => {
  const mesh = buildMesh(new Grid(32)), values = new Float64Array(mesh.nCells), lit = new Float64Array(mesh.nCells), sun = new Float64Array(3);
  let lowest = Infinity, highest = -Infinity;
  for (let day = 0; day < 365; day += 3) for (const hour of [0, 5.25, 9.5, 12, 17.75]) {
    sunlightOnCells(mesh.xCell, sunDirection(timeOf(day, hour), sun), values);
    const mean = areaMean(values, mesh.areaCell);
    lowest = Math.min(lowest, mean);
    highest = Math.max(highest, mean);
    for (let i = 0; i < values.length; i++) lit[i] = values[i] > 0 ? 1 : 0;
    close(areaMean(lit, mesh.areaCell), 0.5, 5e-4, `lit share, day ${day} hour ${hour}`);
  }
  close(lowest, SOLAR_CONSTANT / 4, 0.045, 'lowest mean');
  close(highest, SOLAR_CONSTANT / 4, 0.045, 'highest mean');
});
