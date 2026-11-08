import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { ROWS, BOXES, SIGMA_FILE, logDays, undercurrentCore, dayWeighted, monthAt, windowDays, windowMeans, calibrate, twoSidedTail, studentTail, studentCritical, sigmaDof, designProblems, compareScores, identityCheck, inMooringWindow, mooringProfile, PROFILE_DEPTHS } from '../scripts/scorecard64.mjs';

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

test('Student\x27s t and the degrees of freedom of an RMS over overlapping pairs from independent starts', () => {
  assert.deepEqual([1, 2, 10, 1e7].map((nu) => +studentCritical(nu).toFixed(3)), [12.706, 4.303, 2.228, 1.96]);
  assert.ok(Math.abs(studentTail(2, 1e7) - twoSidedTail(2)) < 2e-6 && Math.abs(studentTail(2.228, 10) - 0.05) < 1e-4);
  assert.ok(Math.abs(sigmaDof([7, 7], 0) - 36 / 4.25) < 1e-12, 'six pairs a start at r = 0: tr C = 3, tr C² = 2.125 each');
  assert.ok(Math.abs(sigmaDof([2, 2, 2], 0.5) - 3) < 1e-12, 'one pair a start: a start a degree of freedom');
  assert.ok(sigmaDof([7, 7], 0.9) > 2 && sigmaDof([7, 7], 0.9) < sigmaDof([7, 7], 0.5) && sigmaDof([7, 7], 0.5) < sigmaDof([7, 7], 0));
});

test('sigma is the saturated spread of the two-window mean, with its serial correlation, degrees of freedom, threshold and the nulls at the gate\x27s leads', () => {
  const rows = (model, proxy = NaN) => Object.fromEntries(ROWS.map(([key]) => [key, { model, ref: 0, proxy }]));
  const year = (start, late, early = 0.1) => Array.from({ length: 12 }, (_, k) => { const lead = k === 11 ? 365 : 30 * (k + 1); return { day: start + lead, rows: rows(lead >= 180 ? late(k) : early) }; });
  const zero = (start) => Array.from({ length: 12 }, (_, k) => ({ day: start + (k === 11 ? 365 : 30 * (k + 1)), rows: rows(NaN, 0) }));
  const twins = [{ label: 'S1', start: 365, a: year(365, (k) => (k % 2 ? 1 : -1)), b: zero(365) }, { label: 'S2', start: 548, a: year(548, () => 2), b: zero(548) }];
  const { rows: sigma, nulls } = calibrate(twins);
  const r = sigma.sst30, pairMeans = [...Array(6).fill(0), ...Array(6).fill(2)];
  assert.equal(r.windows, 14);
  assert.ok(Math.abs(r.sigma - Math.sqrt(pairMeans.reduce((a, x) => a + x * x, 0) / 12)) < 1e-12, `sigma ${r.sigma}`);
  assert.ok(Math.abs(r.perWindow - Math.sqrt((7 + 28) / 14)) < 1e-12);
  assert.ok(Math.abs(r.lag1 - (-6 + 24) / 12 / (35 / 14)) < 1e-12 && Math.abs(r.dof - sigmaDof([7, 7], r.lag1)) < 1e-12);
  assert.ok(Math.abs(r.threshold - studentCritical(r.dof) * r.sigma) < 1e-12 && r.critical > 2);
  assert.deepEqual(r.early.map((x) => +x.toFixed(12)), [0.1, 0.1]);
  assert.ok(Math.abs(r.saturation - 0.1 / r.sigma) < 1e-12 && r.proxy === false);
  assert.deepEqual(r.nullAtLeads, [false, false]);
  assert.deepEqual(nulls.atLeads.map((n) => [n.moved, n.rows]), [[0, ROWS.length], [0, ROWS.length]]);
  assert.ok(Math.abs(twoSidedTail(2) - 0.0455) < 2e-4 && Math.abs(twoSidedTail(0) - 1) < 1e-6);
  const design = { start: { sha256: 'a'.repeat(64), day: 365 }, days: [395, 425] }, score = (reductions) => ({ ...design, means: windowMeans(reductions) });
  const base = score([{ day: 395, rows: rows(1) }, { day: 425, rows: rows(3) }]), moved = score([{ day: 395, rows: rows(2 + 1.01 * r.threshold) }]), inside = score([{ day: 395, rows: rows(2 + 0.99 * r.threshold) }]);
  assert.equal(base.means.sst30.model, 2);
  assert.ok(compareScores(base, moved, sigma).every((c) => c.verdict === 'moved'));
  assert.ok(compareScores(base, inside, sigma).every((c) => c.verdict === 'unresolved'));
});

test('identity: a base read back from JSON, bit-exact rows, rows with no value on either side, the starts and the saves', () => {
  const sigma = Object.fromEntries(ROWS.map(([key]) => [key, { sigma: 1, dof: 8 }])), viaJson = (x) => JSON.parse(JSON.stringify(x));
  const record = (value, extra = {}) => ({ start: { sha256: 'a'.repeat(64), day: 365 }, days: [395, 425], saves: ['s1', 's2'], means: Object.fromEntries(ROWS.map(([key]) => [key, { model: value, ref: 0 }])), ...extra });
  const verdicts = (base, now) => { const comparison = compareScores(viaJson(base), now, sigma); return { comparison, check: identityCheck(viaJson(base), now, comparison) }; };
  let { comparison, check } = verdicts(record(2), record(2));
  assert.ok(comparison.every((c) => c.identical) && check.ok && check.sameSaves);
  ({ check } = verdicts(record(2), record(2, { saves: ['t1', 't2'] })));
  assert.ok(check.ok && !check.sameSaves, 'a fresh bit-identical run');
  const ulp = record(2);
  ulp.means.sst30.model = 2 + 2 ** -51;
  ({ comparison, check } = verdicts(record(2), ulp));
  assert.ok(!comparison[0].identical && comparison[0].verdict === 'unresolved' && !check.ok && check.changed.join() === 'sst30', 'one ulp is a change');
  ({ comparison, check } = verdicts(record(NaN), record(NaN)));
  assert.ok(comparison.every((c) => c.verdict === 'no value' && !c.identical && c.missing.join() === 'base,candidate') && !check.ok && check.empty.length === ROWS.length, 'empty on both sides is no value, not identical');
  const oneSided = record(2);
  oneSided.means.asr.model = NaN;
  ({ comparison, check } = verdicts(record(2), oneSided));
  const asr = comparison.find((c) => c.key === 'asr');
  assert.ok(asr.verdict === 'no value' && !asr.identical && asr.missing.join() === 'candidate' && !check.ok && check.empty.length === 1);
  ({ comparison } = verdicts(oneSided, record(2)));
  assert.deepEqual(comparison.find((c) => c.key === 'asr').missing, ['base']);
  const missingKey = record(2);
  delete missingKey.means.slp;
  ({ comparison } = verdicts(missingKey, record(2)));
  assert.equal(comparison.find((c) => c.key === 'slp').verdict, 'no value');
  const other = record(2, { start: { sha256: 'b'.repeat(64), day: 365 } }), withAny = compareScores(record(2), other, sigma, { anyStart: true });
  assert.ok(withAny.every((c) => c.identical) && !identityCheck(record(2), other, withAny).ok, '--identical refuses another start even under --any-start');
});

test('a base from another start, other leads or no recorded start is refused unless --any-start', () => {
  const means = {}, base = { start: { sha256: 'a'.repeat(64), day: 365 }, days: [395, 425], means };
  assert.deepEqual(designProblems(base, base), []);
  assert.match(designProblems(base, { ...base, start: { sha256: 'b'.repeat(64), day: 365 } }).join(), /started from a{12}…, the candidate from b{12}…/);
  assert.match(designProblems(base, { ...base, start: { sha256: 'a'.repeat(64), day: 548 }, days: [578, 638] }).join(), /leads are \[30,60\], the candidate's \[30,90\]/);
  assert.match(designProblems({ means }, base).join(), /base records no start/);
  assert.throws(() => compareScores({ means }, base, {}), /not the same gate design/);
});

test('the mooring window: the N=64 floors of 1° and 1.5° bind below them, the mesh spacing above; the profile weighs columns by area', () => {
  const at = (lat, lon, spacing) => inMooringWindow(lat, lon, -140, spacing);
  assert.deepEqual([at(0.95, -140, 0.9), at(1.05, -140, 0.9), at(0, -138.55, 0.9), at(0, -138.45, 0.9), at(0, -141.45, 0.9)], [true, false, true, false, true]);
  assert.deepEqual([at(3.5, -140, 3.6), at(3.7, -140, 3.6), at(0, -136.5, 3.6), at(0, -136.3, 3.6)], [true, false, true, false]);
  assert.ok(at(0, 220, 0.9) && !at(0, 140, 3.6), 'longitudes wrap');
  const Z = [20, 30], column = (area, u) => ({ area, z: [10, 40], u: [u, u] });
  assert.deepEqual(mooringProfile([column(3, 1), column(1, 5)], Z).map((p) => p.u), [2, 2]);
  assert.deepEqual(mooringProfile([column(1, 1), { area: 1, z: [25, 40], u: [3, 3] }], Z).map((p) => p.u), [1, 2], 'a column counts only within its own depths');
  assert.equal(PROFILE_DEPTHS[0], 20);
});

test('the committed sigma covers every row with its design, twins, degrees of freedom and expiry, and the rows and boxes are named once', () => {
  const record = JSON.parse(readFileSync(SIGMA_FILE, 'utf8'));
  assert.ok(ROWS.every(([key]) => Number.isFinite(record.rows[key]?.sigma) && record.rows[key].sigma > 0 && record.rows[key].dof >= 2 && record.rows[key].threshold > record.rows[key].sigma && record.rows[key].windows > 0), ROWS.filter(([key]) => !(record.rows[key]?.sigma > 0)).map(([key]) => key).join(', '));
  assert.ok(record.twins.length >= 1 && record.twins.every((t) => Number.isFinite(t.start) && t.days.length > 0) && record.design && record.expires && record.measured && record.definition && record.nulls);
  assert.equal(new Set(ROWS.map(([key]) => key)).size, ROWS.length);
  assert.ok(Object.values(BOXES).every(([kind, s, n, w, e]) => ['land', 'sea', 'all'].includes(kind) && s < n && w < e));
});

test('reduceState on an N=16 fixture whose fields vary by place and mask: every box, mask, unit, the undercurrent\x27s core, the ice cut and the log rows come out as built', async () => {
  const { Grid } = await import('../js/grid.module.js');
  const { gridHeader, stateGrid } = await import('../js/stateGrid.module.js');
  const { stateGeography, createGeography, FRESH_GEOGRAPHY } = await import('../js/geography.module.js');
  const { indexCells, seaIceIndexCut } = await import('../js/seaIceIndex.module.js');
  const { climatologyAt } = await import('../js/ocean/climatology.module.js');
  const { windowMonths } = await import('../scripts/woaSstBias.mjs');
  const { reduceState } = await import('../scripts/scorecard64.mjs');
  const { cellVector } = await import('../js/dynamics/operators.module.js');
  const R = await import('../scripts/references.mjs');
  const header = { N: 16, ...gridHeader(new Grid(16)), ...FRESH_GEOGRAPHY }, { topography, geography } = stateGeography(header);
  const { mesh } = stateGrid(header, topography), made = createGeography(mesh, topography, geography), land = made.land;
  const C = mesh.nCells, E = mesh.nEdges, L = 58, rad = Math.PI / 180, D = 425;
  const fill = (f) => Float32Array.from({ length: C }, (_, i) => f(i));
  const lon = (i) => ((mesh.lonCell[i] / rad + 540) % 360) - 180, lat = (i) => mesh.latCell[i] / rad;
  const wave = (i) => Math.sin(lon(i) * rad) * Math.cos(lat(i) * rad) + lat(i) / 90;
  const concentration = fill((i) => (land[i] ? 0 : lat(i) > 75 ? 0.5 : lat(i) < -70 ? 1 : 0));
  const means = {
    fields: ['sst', 'surfaceT', 'gridBoxSkin', 'temperature', 'precipitation', 'zonalStress', 'mslp', 'concentration'],
    sst: fill((i) => (land[i] ? 250 : 300 + 3 * wave(i))), surfaceT: fill(() => 290), gridBoxSkin: fill((i) => 288 + 5 * wave(i)),
    temperature: fill((i) => (land[i] ? 270 + 4 * wave(i) : 300)), precipitation: fill((i) => 3 + 2 * wave(i) + (land[i] ? 1 : 0)),
    zonalStress: fill((i) => -0.05 + 0.03 * wave(i)), mslp: fill((i) => 100000 + 400 * wave(i)), concentration,
  };
  const tao = R.tao(), wrap = (x) => ((x + 540) % 360) - 180, core = (z) => 1 - ((z - 120) / 60) ** 2;
  const offStation = (lonDeg) => Math.min(...R.TAO_STATIONS.map((st) => Math.abs(wrap(lonDeg - tao.stations[st].longitude))));
  const amplitude = (latDeg, lonDeg) => (1 + 0.3 * Math.sin(lonDeg * rad)) * (1 + 0.004 * offStation(lonDeg) ** 2 + 0.004 * latDeg ** 2);
  const u = new Float64Array(L * E), h = new Float64Array(L * C);
  for (let i = 0; i < C; i++) if (!land[i]) for (let k = 0; k < L; k++) h[k * C + i] = 10;
  for (let e = 0; e < E; e++) {
    const [x, y, z] = [mesh.xEdge[3 * e], mesh.xEdge[3 * e + 1], mesh.xEdge[3 * e + 2]], r = Math.hypot(x, y), across = (-y * mesh.nEdge[3 * e] + x * mesh.nEdge[3 * e + 1]) / r;
    for (let k = 0; k < L; k++) u[k * E + e] = amplitude(Math.atan2(z, r) / rad, Math.atan2(y, x) / rad) * core(10 * k + 5) * across;
  }
  const days = Array.from({ length: 30 }, (_, j) => D - 29 + j);
  const log = days.flatMap((d, j) => [
    `day ${d} (1.0 min): Ts 15.00 °C, ASR ${(230 + j).toFixed(1)} (atmosphere 70.0) OLR ${(238 - j / 2).toFixed(1)} W/m², ps 500–1040 hPa, albedo 0.300, SWCRE ${(-40 - j / 2).toFixed(1)} LWCRE ${(20 + j / 4).toFixed(2)}, clear-sky reflectance 0.15`,
    `sea ice on day ${d}'s mean concentration, sea cells on the index's grids, 15 %: N extent ${(10 + j / 10).toFixed(2)} area ${(8 + j / 20).toFixed(2)}, S extent ${(3 - j / 20).toFixed(2)} area ${(2 - j / 50).toFixed(2)} Mkm²`,
  ]).join('\n');
  const { rows } = reduceState({ ...header, day: D, means, ocean: { h, u } }, log);
  const gap = reduceState({ ...header, day: D, means, ocean: { h, u } }, log.split('\n').filter((line) => !line.includes(`day ${D - 10}`)).join('\n')).rows;
  assert.ok(['asr', 'net', 'swcre', 'lwcre', 'iceNExtent', 'iceNArea', 'iceSExtent', 'iceSArea'].every((key) => Number.isNaN(gap[key].model)), 'a window missing a day has no log rows');

  const inBox = (s, n, w, e) => (i) => lat(i) >= s && lat(i) <= n && lon(i) >= w && lon(i) <= e;
  const mean = (field, keep) => { let sum = 0, area = 0; for (let i = 0; i < C; i++) if (keep(i) && Number.isFinite(field[i])) { sum += mesh.areaCell[i] * field[i]; area += mesh.areaCell[i]; } return sum / area; };
  const woa = R.woa23Sst().sst, months = windowMonths(woa, days), observed = (i) => Number.isFinite(climatologyAt(woa, months, mesh.latCell[i], mesh.lonCell[i]));
  const seaBox = (s, n, w = -180, e = 180) => (i) => !land[i] && observed(i) && inBox(s, n, w, e)(i);
  const near = (key, value, tol = 1e-5) => assert.ok(Math.abs(rows[key].model - value) <= tol, `${key}: ${rows[key].model}, built ${value}`);
  near('sst30', mean(means.sst, seaBox(-30, 30)) - 273.15);
  near('sstCold', mean(means.sst, seaBox(-2, 2, -150, -90)) - 273.15);
  near('sstWarm', mean(means.sst, seaBox(-10, 10, 130, 160)) - 273.15);
  near('ts', mean(means.gridBoxSkin, () => true) - 273.15);
  near('landT', mean(means.temperature, (i) => land[i] && lat(i) >= 40 && lat(i) <= 60));
  const boxes = { amazon: ['land', -10, 3, -70, -50], india: ['land', 15, 28, 72, 88], tropLand: ['land', -15, 15, -180, 180], tropSea: ['sea', -15, 15, -180, 180], io: ['sea', -10, 10, 50, 100], wp: ['sea', -10, 10, 120, 160], global: ['all', -90, 90, -180, 180] };
  for (const [box, [kind, s, n, w, e]] of Object.entries(boxes)) near(`rain_${box}`, mean(means.precipitation, (i) => inBox(s, n, w, e)(i) && (kind === 'all' || (kind === 'land') === !!land[i])));
  const equator = (w, e) => (i) => !land[i] && Math.abs(lat(i)) <= 3 && lon(i) >= w && lon(i) <= e;
  near('stress', 1e3 * mean(means.zonalStress, equator(-160, -100)));
  near('slp', (mean(means.mslp, equator(-180, -160)) - mean(means.mslp, equator(-120, -100))) / 100);
  const meanOver = (g) => days.reduce((a, _, j) => a + g(j), 0) / days.length;
  near('asr', meanOver((j) => 230 + j)); near('net', meanOver((j) => 230 + j - (238 - j / 2))); near('swcre', meanOver((j) => -40 - j / 2)); near('lwcre', meanOver((j) => 20 + j / 4));
  near('iceNExtent', meanOver((j) => 10 + j / 10)); near('iceNArea', meanOver((j) => 8 + j / 20)); near('iceSExtent', meanOver((j) => 3 - j / 20)); near('iceSArea', meanOver((j) => 2 - j / 50));
  const cut = seaIceIndexCut(mesh, indexCells(mesh, made), concentration);
  assert.ok(cut.N.extent > 0 && Math.abs(rows.iceNExtent.proxy - cut.N.extent) < 1e-9 && Math.abs(rows.iceNArea.proxy - 0.5 * cut.N.extent) < 1e-9 && Math.abs(rows.iceSArea.proxy - cut.S.extent) < 1e-9);
  const spacing = 0.9 * Math.sqrt(4 * Math.PI / C) / rad;
  assert.ok(spacing > 1.5 && spacing < 5, `N=16 spacing ${spacing}°: the window is the mesh's, not the N=64 floors`);
  const unit = new Float64Array(E);
  for (let e = 0; e < E; e++) {
    const [x, y, z] = [mesh.xEdge[3 * e], mesh.xEdge[3 * e + 1], mesh.xEdge[3 * e + 2]], r = Math.hypot(x, y);
    unit[e] = amplitude(Math.atan2(z, r) / rad, Math.atan2(y, x) / rad) * (-y * mesh.nEdge[3 * e] + x * mesh.nEdge[3 * e + 1]) / r;
  }
  const v = cellVector(mesh, unit), rebuilt = Float64Array.from({ length: C }, (_, i) => -Math.sin(mesh.lonCell[i]) * v[3 * i] + Math.cos(mesh.lonCell[i]) * v[3 * i + 1]);
  const centers = Array.from({ length: L }, (_, k) => 10 * k + 5), layered = (z) => { const k = Math.min(Math.floor((z - 5) / 10), L - 2), t = (z - centers[k]) / 10; return core(centers[k]) + t * (core(centers[k + 1]) - core(centers[k])); };
  for (const s of R.TAO_STATIONS) {
    const window = (i) => !land[i] && Math.abs(lat(i)) <= Math.max(1, spacing) && Math.abs(wrap(lon(i) - tao.stations[s].longitude)) <= Math.max(1.5, spacing);
    const cs = Array.from({ length: C }, (_, i) => i).filter(window), built = mean(fill((i) => amplitude(lat(i), lon(i))), window);
    assert.equal(rows[`euc_${s}`].cells, cs.length, `${s}: the window holds ${rows[`euc_${s}`].cells} cells, built ${cs.length}`);
    const reconstructed = mean(rebuilt, window), expected = undercurrentCore(Array.from({ length: 67 }, (_, j) => 20 + 5 * j).map((z) => ({ z, u: reconstructed * layered(z) })));
    assert.ok(cs.length > 0 && Math.abs(reconstructed - built) < 0.02 * built, `${s}: the reconstruction's mean amplitude ${reconstructed}, built ${built}`);
    near(`euc_${s}`, expected.speed, 1e-9 * expected.speed); near(`eucDepth_${s}`, expected.depth, 1e-6);
  }
  assert.ok(Math.abs(rows.rain_global.ref - R.gridMean(R.gpcp().grid, 3)) < 0.2 && Math.abs(rows.ts.ref - 15.5) < 0.6, 'references for late April');
});
