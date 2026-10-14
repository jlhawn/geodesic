import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { topographyFromInt16, createGeography } from '../js/geography.module.js';
import { createOcean, EPS, RESTORE_TOLERANCE } from '../js/ocean/layered.module.js';
import { loadClimatology, CLIMATOLOGY_FILE } from '../js/ocean/climatology.module.js';
import { seawaterDensity } from '../js/ocean/seawater.module.js';
import { refitBand, bandWeight, BAND } from '../scripts/bandRefit.mjs';

const DEG = 180 / Math.PI;
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const atlas = await loadClimatology(new URL(`../${CLIMATOLOGY_FILE}`, import.meta.url).pathname);
const mesh = buildMesh(new Grid(16)), C = mesh.nCells;
const geography = createGeography(mesh, topography);
const analytic = createOcean(mesh, { geography });
const ice = new Float64Array(C), surfaceT = Float64Array.from(mesh.latCell, (lat) => 273.15 + 28 * Math.cos(lat) ** 2);
analytic.initialize(surfaceT, ice);
const start = analytic.serialize(), L = start.h.length / C;
const clone = (o) => ({ ...o, h: Float64Array.from(o.h), u: Float64Array.from(o.u), T: Float64Array.from(o.T), S: Float64Array.from(o.S), eta: Float64Array.from(o.eta) });
const lat = (i) => mesh.latCell[i] * DEG;
const column = (o, i) => { let s = 0; for (let k = 0; k < L; k++) s += o.h[k * C + i]; return s; };
const sea = Array.from({ length: C }, (_, i) => i).filter((i) => column(start, i) > 1);
const run = (options = {}) => { const o = clone(start); const report = refitBand(o, mesh, atlas, { ice, ...options }); return { o, report }; };

test('the weight is 1 within 10 degrees, 0 from 15 and linear between', () => {
  assert.equal(bandWeight(0), 1);
  assert.equal(bandWeight(-10), 1);
  assert.equal(bandWeight(15), 0);
  assert.equal(bandWeight(-20), 0);
  assert.ok(Math.abs(bandWeight(12.5) - 0.5) < 1e-12);
  assert.deepEqual(BAND, { inner: 10, outer: 15 });
});

test('the refit keeps each column\'s water, its mixed layer, sea level and velocities, and only the band changes', () => {
  const { o, report } = run();
  assert.ok(report.stats.refit > 20, `columns refit: ${report.stats.refit}`);
  let changed = 0;
  for (const i of sea) {
    const w = bandWeight(lat(i));
    assert.ok(Math.abs(column(o, i) - column(start, i)) <= 1e-9 * column(start, i), `column ${i} water`);
    for (const key of ['h', 'T', 'S']) assert.equal(o[key][i], start[key][i], `mixed layer ${key} of cell ${i}`);
    for (let k = 1; k < L; k++) {
      const n = k * C + i;
      assert.ok(o.h[n] >= EPS - 1e-12 && Number.isFinite(o.T[n]) && Number.isFinite(o.S[n]), `class ${k} of cell ${i}`);
      if (w === 0) assert.ok(o.h[n] === start.h[n] && o.T[n] === start.T[n] && o.S[n] === start.S[n], `cell ${i} at ${lat(i).toFixed(1)} left alone`);
      else if (o.h[n] !== start.h[n]) changed++;
    }
  }
  assert.ok(changed > 100, `${changed} class values moved in the band`);
  assert.deepEqual(Array.from(o.eta), Array.from(start.eta));
  assert.deepEqual(Array.from(o.u), Array.from(start.u));
  for (let i = 0; i < C; i++) if (!sea.includes(i)) for (let k = 0; k < L; k++) assert.equal(o.h[k * C + i], start.h[k * C + i]);
});

test('within 10 degrees the classes are the atlas\'s, near their labels, and a second pass changes nothing', () => {
  const { o } = run(), again = clone(o);
  refitBand(again, mesh, atlas, { ice });
  const inner = sea.filter((i) => Math.abs(lat(i)) <= BAND.inner);
  assert.ok(inner.length > 10);
  const classes = start.densities;
  let checked = 0;
  for (const i of inner) for (let k = 1; k < L; k++) {
    const n = k * C + i;
    assert.ok(Math.abs(again.h[n] - o.h[n]) <= 1e-9 * Math.max(1, o.h[n]), `idempotent h, cell ${i} class ${k}`);
    assert.ok(Math.abs(again.T[n] - o.T[n]) <= 1e-9 && Math.abs(again.S[n] - o.S[n]) <= 1e-9, `idempotent T and S, cell ${i} class ${k}`);
    if (o.h[n] > 5) { checked++; assert.ok(Math.abs(seawaterDensity(o.T[n], o.S[n]) - classes[k - 1]) <= RESTORE_TOLERANCE + 1e-9, `class ${k} of cell ${i} keeps its density`); }
  }
  assert.ok(checked > 50, `${checked} classes over 5 m checked`);
});

test('between 10 and 15 degrees the classes are the weighted blend of the old and the atlas\'s', () => {
  const { o } = run(), { o: pure } = run({ inner: 90, outer: 91 });
  const middle = sea.filter((i) => Math.abs(lat(i)) > BAND.inner + 0.2 && Math.abs(lat(i)) < BAND.outer - 0.2);
  assert.ok(middle.length > 5, `${middle.length} cells in the ramp`);
  for (const i of middle) {
    const w = bandWeight(lat(i));
    for (let k = 1; k < L; k++) {
      const n = k * C + i, oh = start.h[n], ph = pure.h[n];
      const expectH = w * ph + (1 - w) * oh, expectQ = w * ph * pure.T[n] + (1 - w) * oh * start.T[n];
      assert.ok(Math.abs(o.h[n] - expectH) <= 1e-6 * Math.max(1, expectH) + 1e-9, `h, cell ${i} class ${k}: ${o.h[n]} against ${expectH}`);
      assert.ok(Math.abs(o.h[n] * o.T[n] - expectQ) <= 1e-6 * Math.max(1, Math.abs(expectQ)) + 1e-6, `heat content, cell ${i} class ${k}`);
    }
  }
});

test('the heat the report gives is the change of the columns\' content, and heat and salt are free to change', () => {
  const { o, report } = run();
  let heat = 0, moved = 0;
  for (const i of sea) {
    let before = 0, after = 0;
    for (let k = 1; k < L; k++) { before += start.h[k * C + i] * start.T[k * C + i]; after += o.h[k * C + i] * o.T[k * C + i]; }
    heat += mesh.areaCell[i] * (after - before) * 1025 * 3985;
    if (Math.abs(after - before) > 1e-6) moved++;
  }
  assert.ok(moved > 0 && Math.abs(report.totalJ - heat) <= 1e-9 * Math.abs(heat) + 1e6, `report ${report.totalJ} against ${heat}`);
  assert.equal(report.bands.length, 6);
  assert.ok(report.warmClasses.length > 3);
});

test('the refit loads into the layered ocean without a refit of its own and without moving the classes', () => {
  const { o } = run();
  const target = createOcean(mesh, { geography });
  const surface = Float64Array.from(surfaceT);
  target.load(o, surface, ice);
  const back = target.serialize();
  for (const i of sea) {
    assert.ok(Math.abs(column(back, i) - column(o, i)) <= 1e-6, `column ${i}`);
    for (let k = 0; k < L; k++) assert.ok(Math.abs(back.h[k * C + i] - o.h[k * C + i]) <= 1e-6, `h, cell ${i} layer ${k}`);
  }
});
