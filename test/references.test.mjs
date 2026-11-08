import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { gridMesh } from '../js/stateGrid.module.js';
import { gpcp, era5Fields, woa23Sst, globalTs, seaIce, tao, ceres, checkSha256, sampleGrid, boxMean, gridMean, TAO_STATIONS } from '../scripts/references.mjs';

test('each reference loads with its recorded sha256 and shape, and a changed byte fails the check', () => {
  const p = gpcp(), e = era5Fields();
  assert.deepEqual([p.grid.months, p.grid.rows, p.grid.columns, p.grid.values.length], [12, 72, 144, 12 * 72 * 144]);
  assert.deepEqual(Object.fromEntries(Object.entries(e.fields).map(([k, f]) => [k, [f.months, f.rows, f.columns, f.south]])), { t2m: [12, 41, 360, 30], msl: [12, 21, 360, -10], metss: [12, 21, 360, -10], lsm: [1, 81, 360, -10] });
  assert.ok(woa23Sst().sst);
  const bytes = Buffer.from(p.grid.values.buffer.slice(0));
  assert.equal(checkSha256('gpcp', p.meta, bytes), bytes);
  bytes[1000] ^= 1;
  assert.throws(() => checkSha256('gpcp', p.meta, bytes), /sha256/);
});

test('known values: GPCP\x27s global mean, ERA5\x27s easterly stress at 0N 140W, the month-weighted Ts, the March Arctic, the undercurrent cores', () => {
  const p = gpcp(), mean = gridMean(p.grid);
  assert.ok(Math.abs(mean - 2.691) < 5e-4 && Math.abs(mean - p.meta.globalMean) < 5e-4, `GPCP global mean ${mean}`);
  const { metss, msl, t2m, lsm } = era5Fields().fields;
  for (let m = 0; m < 12; m++) {
    const tau = sampleGrid(metss, m, 0, -140);
    assert.ok(tau < 0 && tau > -0.15, `metss ${tau} N/m² at 0N 140W in month ${m + 1}`);
    assert.ok(Math.abs(sampleGrid(msl, m, 0, 120) - 101000) < 1000);
  }
  assert.ok(sampleGrid(t2m, 0, 55, 100) < sampleGrid(t2m, 6, 55, 100) - 20, 'Siberian winter');
  assert.ok(sampleGrid(lsm, 0, 50, 10, { lsm: null }) > 0.99 && sampleGrid(lsm, 0, 0, -140) < 0.01);
  const ts = globalTs();
  assert.equal(ts.annual, 15.259);
  assert.equal(seaIce().extent.N[2], 15.027);
  const t = tao();
  assert.deepEqual(TAO_STATIONS.map((s) => t.stations[s].core.speed), [0.498, 0.697, 1.076, 0.952]);
  assert.ok(TAO_STATIONS.every((s) => t.stations[s].profile.every((d) => d.years >= 10 && d.days >= 30)));
  assert.ok(Object.keys(ceres().fields).length > 0);
});

test('with lsm a land cell takes the land points around it and a sea cell the sea points, and a box averages its cells by area', () => {
  const g = { south: 0, west: 0, step: 1, rows: 2, columns: 360, months: 1, values: new Float32Array(720) };
  const lsm = { south: 0, west: 0, step: 1, rows: 2, columns: 360, months: 1, values: new Float32Array(720) };
  g.values[0] = 10; g.values[1] = 20; g.values[360] = 30; g.values[361] = 40;
  lsm.values[0] = 1; lsm.values[1] = 0; lsm.values[360] = 0.5; lsm.values[361] = 0;
  assert.equal(sampleGrid(g, 0, 0.5, 0.5), 25);
  assert.ok(Math.abs(sampleGrid(g, 0, 0.5, 0.5, { lsm, land: true }) - (10 + 0.5 * 30) / 1.5) < 1e-6);
  assert.ok(Math.abs(sampleGrid(g, 0, 0.5, 0.5, { lsm, land: false }) - (20 + 0.5 * 30 + 40) / 2.5) < 1e-6);
  assert.equal(sampleGrid(g, 0, 0.5, 359.5), sampleGrid(g, 0, 0.5, -0.5));
  assert.throws(() => sampleGrid(g, 0, 0.5, 0.5, { lsm }), /needs land/);
  const mesh = gridMesh(new Grid(16)), cells = [0, 1, 2];
  const flat = { south: -90, west: 0, step: 1, rows: 181, columns: 360, months: 1, values: new Float32Array(181 * 360).fill(7) };
  assert.ok(Math.abs(boxMean(mesh, flat, 0, cells) - 7) < 1e-6);
});
