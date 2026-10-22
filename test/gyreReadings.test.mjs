import { test } from 'node:test';
import assert from 'node:assert/strict';
import { cellVector } from '../js/dynamics/operators.module.js';
import { basinCase, caseMesh, meanSpacing, RHO } from '../explain/oceanCases.module.js';
import { closureCoefficient } from '../js/ocean/layered.module.js';
import { gyreReadings, westernCurrent, zonalExtremes, rowWidths, DEPTH } from '../explain/gyreReadings.module.js';
import { EARTH_RADIUS } from '../js/mesh.module.js';
import { runThreaded } from './helpers/oceanThreads.mjs';

const DEG = Math.PI / 180, PER_DEGREE = EARTH_RADIUS * DEG * Math.cos(30 * DEG);
const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);

test('a western current is read from the coast to where the row first turns south', () => {
  const longitude = [-36, -32, -28, -24, -20], north = [40, 20, -10, -5, -5];
  assert.deepEqual(Array.from(rowWidths(longitude)), [4, 4, 4, 4, 4]);
  const current = westernCurrent({ longitude, north });
  close(current.coast, -38, 1e-12, 'coast');
  close(current.end, -32 + 4 * 20 / 30, 1e-12, 'end');
  close(current.width, (current.end + 38) * PER_DEGREE, 1e-6, 'width');
  close(current.peak, 40 / DEPTH, 1e-15, 'peak');
  close(current.transport, (40 + 20) * 4 * PER_DEGREE / 1e6, 1e-9, 'transport');
  assert.equal(westernCurrent({ longitude, north: [-1, 2, 3, 4, 5] }), null);
  const allNorth = westernCurrent({ longitude, north: [1, 1, 1, 1, 1] });
  close(allNorth.end, -18, 1e-12, 'a row that never turns south ends at its eastern edge');
});

test('the zonal extremes are the strongest westward and eastward bands between 10 and 50 degrees', () => {
  const latitude = [2.5, 12.5, 17.5, 37.5, 47.5, 52.5], east = [-900, -100, -269, 231, 100, 900];
  const { westward, eastward } = zonalExtremes({ latitude, east });
  assert.deepEqual(westward, { speed: -269 / DEPTH, latitude: 17.5 });
  assert.deepEqual(eastward, { speed: 231 / DEPTH, latitude: 37.5 });
});

test('after 20 days the readings agree with the ocean\'s own layers, a direct row integral and Sverdrup\'s formula', async () => {
  const sim = basinCase({ days: 20 });
  const frame = await runThreaded(sim), { mesh, ocean } = sim, reading = gyreReadings(frame);
  const C = mesh.nCells, E = mesh.nEdges, L = ocean.layers;
  const row = sim.layout().parallel;

  for (const i of row) {
    let column = 0;
    for (let k = 0; k < L; k++) column += ocean.h[k * C + i];
    close(column, DEPTH + ocean.eta[i], 1e-6, `water column in cell ${i}`);
  }

  const north = new Float64Array(row.length), vector = new Float64Array(3 * C);
  for (let k = 0; k < L; k++) {
    cellVector(mesh, ocean.u.subarray(k * E, (k + 1) * E), vector);
    row.forEach((i, s) => {
      const x = mesh.xCell[3 * i], y = mesh.xCell[3 * i + 1], z = mesh.xCell[3 * i + 2], r = Math.hypot(x, y);
      north[s] += ocean.h[k * C + i] * (vector[3 * i] * -z * x / r + vector[3 * i + 1] * -z * y / r + vector[3 * i + 2] * r);
    });
  }
  const gap = (a, b) => Math.acos(Math.min(1, (mesh.xCell[3 * a] * mesh.xCell[3 * b] + mesh.xCell[3 * a + 1] * mesh.xCell[3 * b + 1] + mesh.xCell[3 * a + 2] * mesh.xCell[3 * b + 2]) / Math.hypot(mesh.xCell[3 * a], mesh.xCell[3 * a + 1], mesh.xCell[3 * a + 2]) / Math.hypot(mesh.xCell[3 * b], mesh.xCell[3 * b + 1], mesh.xCell[3 * b + 2]))) * mesh.radius;
  let k = 0, transport = 0;
  const widthOf = (s) => (s === 0 ? gap(row[0], row[1]) : s === row.length - 1 ? gap(row[s - 1], row[s]) : 0.5 * (gap(row[s - 1], row[s]) + gap(row[s], row[s + 1])));
  transport += north[0] * widthOf(0);
  while (north[k + 1] > 0) { k++; transport += north[k] * widthOf(k); }
  const end = gap(row[0], row[k]) + gap(row[k], row[k + 1]) * north[k] / (north[k] - north[k + 1]), width = 0.5 * gap(row[0], row[1]) + end;

  const current = reading.current;
  console.log(`day ${reading.day}: western current ${(current.peak * 100).toFixed(3)} cm/s at its fastest, ${(current.width / 1e3).toFixed(0)} km wide, ${current.transport.toFixed(1)} Sv; directly from the layers ${(Math.max(...north.slice(0, k + 1)) / DEPTH * 100).toFixed(3)} cm/s, ${(width / 1e3).toFixed(0)} km, ${(transport / 1e6).toFixed(1)} Sv; band ${(reading.west * 100).toFixed(3)} cm/s north, interior ${(reading.interior * 100).toFixed(4)} cm/s, ${reading.ratio.toFixed(1)} times; Sverdrup ${(reading.sverdrup.speed * 100).toFixed(4)} cm/s, ${reading.sverdrup.transport.toFixed(1)} Sv`);
  close(current.width, width, 0.1 * width, 'width of the western current');
  close(current.transport, transport / 1e6, 0.1 * transport / 1e6, 'transport of the western current');
  close(current.peak, Math.max(...north.slice(0, k + 1)) / DEPTH, 0.1 * current.peak, 'fastest depth-mean speed');

  const beta = 2 * mesh.omega * Math.cos(30 * DEG) / mesh.radius, curl = -0.1 * Math.PI / (30 * DEG * mesh.radius);
  close(reading.sverdrup.speed, curl / (RHO * beta) / DEPTH, 1e-12, 'Sverdrup speed');
  close(reading.sverdrup.transport, curl / (RHO * beta) * 80 * PER_DEGREE / 1e6, 1e-6, 'Sverdrup transport');
  assert.ok(reading.interior < 0 && reading.ratio >= 12 && reading.ratio <= 16, `the western current runs ${reading.ratio} times as fast as the interior's return`);
  assert.ok(current.width > 600e3 && current.width < 1300e3, `the western current is about two cells wide: ${current.width}`);
});

test('with no walls the water runs around the planet with the wind and no gyre forms, while the smoothing sets the width of the western current', async () => {
  const open = basinCase({ walls: false }), f = await runThreaded(open, 100), r = gyreReadings(f);
  console.log(`no walls, day 100: ${f.northward.toFixed(2)} Sv north and ${f.southward.toFixed(2)} Sv south across 30°N; westward ${(100 * r.zonal.westward.speed).toFixed(2)} cm/s at ${r.zonal.westward.latitude}°, eastward ${(100 * r.zonal.eastward.speed).toFixed(2)} cm/s at ${r.zonal.eastward.latitude}°`);
  assert.equal(r.current, null);
  assert.ok(f.northward < 3 && -f.southward < 3, `${f.northward} Sv north and ${f.southward} Sv south`);
  assert.ok(r.zonal.westward.latitude === 17.5 && r.zonal.eastward.latitude === 37.5);
  assert.ok(r.zonal.westward.speed < -0.055 && r.zonal.westward.speed > -0.08 && r.zonal.eastward.speed > 0.05 && r.zonal.eastward.speed < 0.065);
  const early = gyreReadings(await runThreaded(basinCase({ walls: false }), 50));
  assert.ok(-r.zonal.westward.speed > -early.zonal.westward.speed && r.zonal.eastward.speed > early.zonal.eastward.speed, 'still speeding up on day 100');

  const spacing = meanSpacing(caseMesh(16)), beta = 2 * caseMesh(16).omega * Math.cos(30 * DEG) / caseMesh(16).radius;
  const nu4 = closureCoefficient(spacing, 12), natural = (nu4 / beta) ** 0.2;
  console.log(`ν₄ ${nu4.toExponential(3)} m⁴/s at N = 16 (cells ${(spacing / 1e3).toFixed(0)} km apart), natural width ${(natural / 1e3).toFixed(0)} km; at N = 128 ${(closureCoefficient(meanSpacing(caseMesh(128)), 12) / 1e0).toExponential(3)} m⁴/s, natural width ${((closureCoefficient(meanSpacing(caseMesh(128)), 12) / beta) ** 0.2 / 1e3).toFixed(0)} km`);
  assert.ok(Math.abs(nu4 / 1.27e16 - 1) < 0.01);
  assert.ok(Math.abs(natural / 1e3 - 230) < 5);
  const weak = gyreReadings(await runThreaded(basinCase({ days: 20, ocean: { closureHours: 48 } })));
  console.log(`smoothing four times weaker, day 20: western current ${(weak.current.width / 1e3).toFixed(0)} km wide, ${weak.current.transport.toFixed(1)} Sv`);
  assert.ok(Math.abs(closureCoefficient(spacing, 48) * 4 / nu4 - 1) < 1e-9);
  assert.ok(weak.current.width > 600e3 && weak.current.width < 780e3, `${weak.current.width}`);
});
