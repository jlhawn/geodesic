import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { readFrames } from '../explain/frames.module.js';
import { classOverturning, classSection } from '../explain/oceanRun.module.js';
import { contours, overturningGrid, overturningCells, northwardByClass, THICK } from '../explain/seaZones.module.js';

const bytes = readFileSync(new URL('../explain/data/ocean_N16.bin', import.meta.url));
const { header } = readFrames('ocean_N16.bin', bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength));
const L = header.layers, M = L + 1;
const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const grid = (ni, nj, f) => Float64Array.from({ length: ni * nj }, (_, n) => f(Math.floor(n / nj), n % nj));
const onEdge = (values, nj, [x, y]) => {
  const at = (a, b) => values[a * nj + b];
  if (Number.isInteger(x)) { const j = Math.floor(y), t = y - j; return (1 - t) * at(x, j) + t * at(x, Math.min(nj - 1, j + 1)); }
  const i = Math.floor(x), s = x - i;
  return (1 - s) * at(i, y) + s * at(i + 1, y);
};
const highOnRight = (values, nj, line, level) => {
  const above = ([i, j]) => values[i * nj + j] > level;
  for (let n = 0; n + 1 < line.points.length; n++) {
    const [a, b] = [line.points[n], line.points[n + 1]], d = [b[0] - a[0], b[1] - a[1]];
    const side = ([i, j]) => (i - a[0]) * -d[1] + (j - a[1]) * d[0];
    for (const [x, y] of [a, b]) {
      const ends = Number.isInteger(x) ? [[x, Math.floor(y)], [x, Math.floor(y) + 1]] : [[Math.floor(x), y], [Math.floor(x) + 1, y]];
      for (const corner of ends) if (side(corner) !== 0 && (side(corner) > 0) !== above(corner)) return false;
    }
  }
  return true;
};

test('a contour around a peak is closed, lies on the level and runs with the larger values on its right', () => {
  const ni = 21, nj = 17, peak = grid(ni, nj, (i, j) => Math.exp(-((i - 10) ** 2 + (j - 8) ** 2) / 20));
  const lines = contours(peak, ni, nj, 0.5);
  assert.equal(lines.length, 1);
  const [line] = lines;
  assert.equal(line.closed, true);
  for (const p of line.points) close(onEdge(peak, nj, p), 0.5, 1e-12, `value at ${p}`);
  assert.ok(highOnRight(peak, nj, line, 0.5));
  const top = line.points.reduce((best, p) => (p[1] < best[1] ? p : best)), next = line.points[(line.points.indexOf(top) + 1) % line.points.length];
  assert.ok(next[0] > top[0], 'over a peak drawn with j down the page, the line runs clockwise');
  const trough = grid(ni, nj, (i, j) => -peak[i * nj + j]), [around] = contours(trough, ni, nj, -0.5);
  assert.ok(highOnRight(trough, nj, around, -0.5));
  const low = around.points.reduce((best, p) => (p[1] < best[1] ? p : best)), after = around.points[(around.points.indexOf(low) + 1) % around.points.length];
  assert.ok(after[0] < low[0], 'around a trough it runs counterclockwise');
});

test('an open contour runs edge to edge, and a value growing down the page makes it run to the right', () => {
  const ni = 6, nj = 5, ramp = grid(ni, nj, (i, j) => j);
  const lines = contours(ramp, ni, nj, 2.5);
  assert.equal(lines.length, 1);
  assert.equal(lines[0].closed, false);
  assert.equal(lines[0].points.length, ni);
  lines[0].points.forEach((p, n) => { close(p[0], n, 1e-12, 'i'); close(p[1], 2.5, 1e-12, 'j'); });
  const across = contours(grid(ni, nj, (i) => i), ni, nj, 2.5);
  assert.ok(across[0].points.at(-1)[1] < across[0].points[0][1], 'a value growing to the right makes the line run up the page');
});

test('a saddle splits into two lines, each with the larger values on its right', () => {
  const ni = 4, nj = 4, saddle = grid(ni, nj, (i, j) => (i - 1.5) * (j - 1.5) + 0.01 * i);
  for (const level of [0.3, -0.3]) {
    const lines = contours(saddle, ni, nj, level);
    assert.equal(lines.length, 2, `level ${level}`);
    for (const line of lines) assert.ok(highOnRight(saddle, nj, line, level), `level ${level}`);
  }
});

test('the overturning grid is zero at the poles and the surface, carries only the free surface’s own flow at the floor, and keeps the class streamfunction between them', () => {
  for (const month of [null, 'March', 'June', 'September', 'December']) {
    const g = overturningGrid(header, month), flow = classOverturning(header, month);
    assert.equal(g.ni, 37);
    assert.equal(g.nj, M);
    assert.deepEqual([g.lats[0], g.lats[1], g.lats[18], g.lats[36]], [-90, -85, 0, 90]);
    for (let j = 0; j < M; j++) { assert.equal(g.psi[j], 0); assert.equal(g.psi[36 * M + j], 0); }
    for (let i = 0; i < g.ni; i++) { assert.equal(g.psi[i * M], 0); close(g.psi[i * M + L], i && i < 36 ? flow.etaNet[i - 1] : 0, 0.05, `floor at ${g.lats[i]}`); }
    if (!month) for (let i = 0; i < g.ni; i++) assert.ok(Math.abs(g.psi[i * M + L]) < 0.2);
    for (let c = 0; c < 35; c++) for (let k = 0; k <= L; k++) { assert.equal(g.psi[(c + 1) * M + k], flow.psi[c * M + k]); assert.equal(g.depth[(c + 1) * M + k], flow.depth[c * M + k]); }
    const base = classSection(header, month).base;
    for (let k = 1; k <= L; k++) { assert.equal(g.depth[k], base[k - 1]); assert.equal(g.depth[36 * M + k], base[35 * L + k - 1]); }
    for (let i = 0; i < g.ni; i++) for (let k = 0; k < L; k++) assert.ok(g.depth[i * M + k + 1] >= g.depth[i * M + k] - 1e-9);
  }
});

test('the overturning cells the figure and the page read: the shallow tropical cells, the cells under the westerlies and the deep polar cells', () => {
  const year = overturningCells(header);
  close(year.south, -58.4, 0.05, 'south across 10°S');
  close(year.north, 52.4, 0.05, 'north across 10°N');
  close(year.mixed(-35), 35.1, 0.05, 'north across 35°S');
  close(year.mixed(35), -34.6, 0.05, 'south across 35°N');
  assert.deepEqual([year.deepNorth.latitude, Math.round(year.deepNorth.depth)], [55, 3190]);
  assert.deepEqual([year.deepSouth.latitude, Math.round(year.deepSouth.depth)], [-55, 3187]);
  close(year.deepNorth.psi, 32.7, 0.05, 'the northern deep cell');
  close(year.deepSouth.psi, -26.2, 0.05, 'the southern deep cell');
  const months = Object.fromEntries(['March', 'June', 'September', 'December'].map((m) => [m, overturningCells(header, m)]));
  close(months.September.south, -90.1, 0.05, 'September, 10°S');
  close(months.September.north, 24.1, 0.05, 'September, 10°N');
  close(months.March.south, -26.9, 0.05, 'March, 10°S');
  close(months.March.north, 81.2, 0.05, 'March, 10°N');
  close(months.December.deepNorth.psi, 44.9, 0.05, 'December, the northern deep cell');
  close(months.June.deepNorth.psi, 21.3, 0.05, 'June, the northern deep cell');
  const g = overturningGrid(header), inner = contours(g.psi, g.ni, g.nj, -55);
  assert.equal(inner.length, 1);
  assert.ok(inner[0].closed);
  for (const [i] of inner[0].points) assert.ok(g.lats[0] + 5 * i > -11.5 && g.lats[0] + 5 * i < -8.5);
});

test('the northward current by class is each class’s transport over its cross-section, and matches the section where both exist', () => {
  const R = header.constants.radius, v = northwardByClass(header), base = header.section.base;
  assert.equal(v.values.length, 35 * L);
  const c = header.overturning.caps.indexOf(-10), h = 0.5 * (base[c * L] + base[(c + 1) * L]);
  close(v.values[c * L], 100 * header.overturning.mass[c * L] * 1e6 / (2 * Math.PI * R * Math.cos(10 * Math.PI / 180) * h), 1e-12, 'mixed layer at 10°S');
  close(v.values[c * L], -2.60, 0.005, 'mixed layer at 10°S, cm/s');
  close(v.values[header.overturning.caps.indexOf(10) * L], 2.41, 0.005, 'mixed layer at 10°N, cm/s');
  for (const lat of [-35, -20, -10, 10, 20, 35]) {
    const cap = header.overturning.caps.indexOf(lat), sectionMean = 0.5 * (header.section.v[cap * L] + header.section.v[(cap + 1) * L]);
    close(v.values[cap * L], sectionMean, 0.5, `mixed layer at ${lat}°`);
  }
  for (let n = 0; n < v.values.length; n++) {
    const cap = Math.floor(n / L), k = n % L, thick = (b) => base[b * L + k] - (k ? base[b * L + k - 1] : 0);
    assert.equal(v.values[n] === null, 0.5 * (thick(cap) + thick(cap + 1)) < THICK);
  }
  const september = northwardByClass(header, 'September'), flow = header.overturning.months.find((m) => m.name === 'September'), sb = header.section.months.find((m) => m.name === 'September').base;
  close(september.values[c * L], 100 * flow.mass[c * L] * 1e6 / (2 * Math.PI * R * Math.cos(10 * Math.PI / 180) * 0.5 * (sb[c * L] + sb[(c + 1) * L])), 1e-12, 'September, mixed layer at 10°S');
  assert.ok(september.values[c * L] < -3.5, 'in September the mixed layer runs south across 10°S faster than over the year');
});
