import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { ROWS, BOXES, SIGMA_FILE, logDays, undercurrentCore, dayWeighted, monthAt, windowDays, windowMeans, calibrate, twoSidedTail, compareScores } from '../scripts/scorecard64.mjs';

test('the log parser reads the day-mean TOA line and the daily sea-ice cut', () => {
  const text = [
    'day 410 (1.4 min): Ts 16.14 °C, ASR 231.5 (atmosphere 77.0) OLR 234.5 W/m², ps 539–1035 hPa, albedo 0.309, SWCRE -50.7 LWCRE 26.2, clear-sky reflectance 0.1573',
    'sea ice on day 410\'s mean concentration, sea cells, 15 %: N extent 12.14 area 10.73, S extent 1.59 area 0.89 Mkm²',
    'regions after 2 days',
  ].join('\n');
  const { days, ice } = logDays(text);
  assert.deepEqual(days.get(410), { asr: 231.5, olr: 234.5, swcre: -50.7, lwcre: 26.2 });
  assert.deepEqual(ice.get(410), { nExtent: 12.14, nArea: 10.73, sExtent: 1.59, sArea: 0.89 });
  assert.equal(days.size + ice.size, 2);
});

test('the undercurrent core is the profile\x27s eastward maximum refined by a parabola', () => {
  const profile = Array.from({ length: 60 }, (_, j) => { const z = 20 + 5 * j; return { z, u: 0.9 - ((z - 107.5) / 40) ** 2 }; });
  const { speed, depth } = undercurrentCore(profile);
  assert.ok(Math.abs(speed - 0.9) < 1e-12 && Math.abs(depth - 107.5) < 1e-9, `${speed} m/s at ${depth} m`);
  assert.deepEqual(undercurrentCore([{ z: 20, u: 0.5 }, { z: 25, u: 0.2 }, { z: 30, u: NaN }]), { speed: 0.5, depth: 20 });
  assert.deepEqual(undercurrentCore([{ z: 20, u: NaN }]), { speed: NaN, depth: NaN });
});

test('window weights: day 0 is the March equinox and a window\x27s reference weighs each month by its hours', () => {
  assert.deepEqual([monthAt(0), monthAt(12), monthAt(365), monthAt(-1)], [3, 4, 3, 3]);
  assert.deepEqual(windowDays(425).slice(0, 2), [396, 397]);
  assert.equal(dayWeighted(395, () => 7), 7);
  const april = dayWeighted(395, (m) => (m === 4 ? 1 : 0));
  assert.ok(april > 0.6 && april < 0.75, `${april} of days 366–395 fall in April`);
});

test('sigma is the saturated spread of the two-window mean, with its serial correlation, degrees of freedom and the nulls at the gate\x27s leads', () => {
  const rows = (model, proxy = NaN) => Object.fromEntries(ROWS.map(([key]) => [key, { model, ref: 0, proxy }]));
  const year = (start, late, early = 0.1) => Array.from({ length: 12 }, (_, k) => { const lead = k === 11 ? 365 : 30 * (k + 1); return { day: start + lead, rows: rows(lead >= 180 ? late(k) : early) }; });
  const zero = (start) => Array.from({ length: 12 }, (_, k) => ({ day: start + (k === 11 ? 365 : 30 * (k + 1)), rows: rows(NaN, 0) }));
  const twins = [{ label: 'S1', start: 365, a: year(365, (k) => (k % 2 ? 1 : -1)), b: zero(365) }, { label: 'S2', start: 548, a: year(548, () => 2), b: zero(548) }];
  const { rows: sigma, nulls } = calibrate(twins);
  const r = sigma.sst30, pairMeans = [...Array(6).fill(0), ...Array(6).fill(2)];
  assert.equal(r.windows, 14);
  assert.ok(Math.abs(r.sigma - Math.sqrt(pairMeans.reduce((a, x) => a + x * x, 0) / 12)) < 1e-12, `sigma ${r.sigma}`);
  assert.ok(Math.abs(r.perWindow - Math.sqrt((7 + 28) / 14)) < 1e-12);
  assert.ok(Math.abs(r.lag1 - (-6 + 24) / 12 / (35 / 14)) < 1e-12 && Math.abs(r.dof - 14 * (1 - r.lag1) / (1 + r.lag1)) < 1e-9);
  assert.deepEqual(r.early.map((x) => +x.toFixed(12)), [0.1, 0.1]);
  assert.ok(Math.abs(r.saturation - 0.1 / r.sigma) < 1e-12 && r.proxy === false);
  assert.deepEqual(r.nullAtLeads, [false, false]);
  assert.deepEqual(nulls.atLeads.map((n) => [n.moved, n.rows]), [[0, ROWS.length], [0, ROWS.length]]);
  assert.ok(Math.abs(twoSidedTail(2) - 0.0455) < 2e-4 && Math.abs(twoSidedTail(0) - 1) < 1e-6);
  const base = windowMeans([{ day: 395, rows: rows(1) }, { day: 425, rows: rows(3) }]), moved = windowMeans([{ day: 395, rows: rows(2 + 2.01 * r.sigma) }]), inside = windowMeans([{ day: 395, rows: rows(2 + 1.99 * r.sigma) }]);
  assert.equal(base.sst30.model, 2);
  assert.ok(compareScores(base, moved, sigma).every((c) => c.verdict === 'moved'));
  assert.ok(compareScores(base, inside, sigma).every((c) => c.verdict === 'unresolved'));
  assert.equal(compareScores(base, windowMeans([{ day: 395, rows: rows(NaN) }]), sigma)[0].verdict, 'no value');
});

test('the committed sigma covers every row with its design, twins, degrees of freedom and expiry, and the rows and boxes are named once', () => {
  const record = JSON.parse(readFileSync(SIGMA_FILE, 'utf8'));
  assert.ok(ROWS.every(([key]) => Number.isFinite(record.rows[key]?.sigma) && record.rows[key].sigma > 0 && record.rows[key].dof >= 1 && record.rows[key].windows > 0), ROWS.filter(([key]) => !(record.rows[key]?.sigma > 0)).map(([key]) => key).join(', '));
  assert.ok(record.twins.length >= 1 && record.twins.every((t) => Number.isFinite(t.start) && t.days.length > 0) && record.design && record.expires && record.measured && record.definition && record.nulls);
  assert.equal(new Set(ROWS.map(([key]) => key)).size, ROWS.length);
  assert.ok(Object.values(BOXES).every(([kind, s, n, w, e]) => ['land', 'sea', 'all'].includes(kind) && s < n && w < e));
});

test('reduceState on an N=16 fixture: boxes, masks, units, the undercurrent\x27s core, the ice cut and the log rows come out as built', async () => {
  const { Grid } = await import('../js/grid.module.js');
  const { gridHeader, stateGrid } = await import('../js/stateGrid.module.js');
  const { stateGeography, createGeography, FRESH_GEOGRAPHY } = await import('../js/geography.module.js');
  const { indexCells, seaIceIndexCut } = await import('../js/seaIceIndex.module.js');
  const { reduceState } = await import('../scripts/scorecard64.mjs');
  const R = await import('../scripts/references.mjs');
  const header = { N: 16, ...gridHeader(new Grid(16)), ...FRESH_GEOGRAPHY }, { topography, geography } = stateGeography(header);
  const { mesh } = stateGrid(header, topography), made = createGeography(mesh, topography, geography), land = made.land;
  const C = mesh.nCells, E = mesh.nEdges, L = 58, deg = 180 / Math.PI, D = 425;
  const fill = (f) => Float32Array.from({ length: C }, (_, i) => f(i));
  const lon = (i) => ((mesh.lonCell[i] * deg + 540) % 360) - 180, lat = (i) => mesh.latCell[i] * deg;
  const concentration = fill((i) => (land[i] ? 0 : lat(i) > 75 ? 0.5 : lat(i) < -70 ? 1 : 0));
  const means = {
    fields: ['sst', 'surfaceT', 'gridBoxSkin', 'temperature', 'precipitation', 'zonalStress', 'mslp', 'concentration'],
    sst: fill((i) => (land[i] ? NaN : 300)), surfaceT: fill(() => 290), gridBoxSkin: fill(() => 288.15), temperature: fill(() => 280),
    precipitation: fill(() => 3), zonalStress: fill(() => -0.05), mslp: fill((i) => (lon(i) < -140 ? 100100 : 100000)), concentration,
  };
  const core = (z) => 1 - ((z - 120) / 60) ** 2, u = new Float64Array(L * E), h = new Float64Array(L * C);
  for (let i = 0; i < C; i++) if (!land[i]) for (let k = 0; k < L; k++) h[k * C + i] = 10;
  for (let k = 0; k < L; k++) for (let e = 0; e < E; e++) {
    const [x, y] = [mesh.xEdge[3 * e], mesh.xEdge[3 * e + 1]], r = Math.hypot(x, y), east = [-y / r, x / r, 0];
    u[k * E + e] = core(10 * k + 5) * (east[0] * mesh.nEdge[3 * e] + east[1] * mesh.nEdge[3 * e + 1]);
  }
  const days = Array.from({ length: 30 }, (_, j) => D - 29 + j);
  const log = days.flatMap((d) => [
    `day ${d} (1.0 min): Ts 15.00 °C, ASR 240.0 (atmosphere 70.0) OLR 238.0 W/m², ps 500–1040 hPa, albedo 0.300, SWCRE -45.0 LWCRE 25.0, clear-sky reflectance 0.15`,
    `sea ice on day ${d}'s mean concentration, sea cells on the index's grids, 15 %: N extent 10.00 area 8.00, S extent 3.00 area 2.00 Mkm²`,
  ]).join('\n');
  const { rows } = reduceState({ ...header, day: D, means, ocean: { h, u } }, log);
  const near = (key, value, tol = 1e-4) => assert.ok(Math.abs(rows[key].model - value) <= tol, `${key}: ${rows[key].model}, built ${value}`);
  for (const key of ['sst30', 'sstCold', 'sstWarm']) near(key, 26.85);
  near('ts', 15); near('landT', 280);
  for (const box of ['amazon', 'india', 'tropLand', 'tropSea', 'io', 'wp', 'global']) near(`rain_${box}`, 3);
  near('stress', -50); near('slp', 1);
  near('asr', 240); near('net', 2); near('swcre', -45); near('lwcre', 25);
  near('iceNExtent', 10); near('iceNArea', 8); near('iceSExtent', 3); near('iceSArea', 2);
  const cut = seaIceIndexCut(mesh, indexCells(mesh, made), concentration);
  assert.ok(cut.N.extent > 0 && Math.abs(rows.iceNExtent.proxy - cut.N.extent) < 1e-9 && Math.abs(rows.iceNArea.proxy - 0.5 * cut.N.extent) < 1e-9 && Math.abs(rows.iceSArea.proxy - cut.S.extent) < 1e-9);
  for (const s of R.TAO_STATIONS) {
    assert.ok(rows[`euc_${s}`].cells > 0 && Math.abs(rows[`euc_${s}`].model - 1) < 0.05, `${s}: core ${rows[`euc_${s}`].model} m/s`);
    near(`eucDepth_${s}`, 120, 3);
  }
  assert.ok(Math.abs(rows.rain_global.ref - R.gridMean(R.gpcp().grid, 3)) < 0.2 && Math.abs(rows.ts.ref - 15.5) < 0.6, 'references for late April');
});
