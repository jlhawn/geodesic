import { test } from 'node:test';
import assert from 'node:assert/strict';
import { cellVector } from '../js/dynamics/operators.module.js';
import { basinCase, runCase, RHO } from '../explain/oceanCases.module.js';
import { gyreReadings, westernCurrent, zonalExtremes, rowWidths, DEPTH } from '../explain/gyreReadings.module.js';
import { EARTH_RADIUS } from '../js/mesh.module.js';

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

test('after 20 days the readings agree with the ocean\'s own layers, a direct row integral and Sverdrup\'s formula', () => {
  const sim = basinCase({ days: 20 });
  const frame = runCase(sim), { mesh, ocean } = sim, reading = gyreReadings(frame);
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
