import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, writeFileSync, mkdtempSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { spawnSync } from 'node:child_process';
import { SECTIONS, sectionName, sectionOf, monthOfDay, utcHourOfDay, era5Sections, MONTHS, sectionTable, readSections, pairTable } from '../scripts/equatorialSections.mjs';

const era5 = JSON.parse(readFileSync(new URL('../data/era5_equatorial_stress.json', import.meta.url), 'utf8'));

test('the sections run in 10 degrees from 120E to 80W, whatever the longitude\'s turn', () => {
  assert.deepEqual(SECTIONS.map(sectionName), ['120E-130E', '130E-140E', '140E-150E', '150E-160E', '160E-170E', '170E-180', '180-170W', '170W-160W', '160W-150W', '150W-140W', '140W-130W', '130W-120W', '120W-110W', '110W-100W', '100W-90W', '90W-80W']);
  for (const [lon, s] of [[119.99, -1], [120, 0], [129.99, 0], [130, 1], [179.99, 5], [180, 6], [-180, 6], [-170, 7], [-80.01, 15], [-80, -1], [280, -1], [480, 0], [-240, 0], [0, -1]]) assert.equal(sectionOf(lon), s, `longitude ${lon}`);
});

test('a model day falls in its calendar month, day 0 the March equinox on March 20, and begins at 12 UTC', () => {
  for (const [day, month] of [[0, 'March'], [11, 'March'], [12, 'April'], [-1, 'March'], [-20, 'February'], [882, 'August'], [892.9, 'August'], [894, 'August'], [895, 'September'], [286, 'December'], [287, 'January'], [652, 'January']]) assert.equal(MONTHS[monthOfDay(day)], month, `day ${day}`);
  for (const [day, hour] of [[882, 12], [882.25, 18], [882.5, 0], [882.75, 6], [-0.25, 6], [882.999, 12]]) assert.equal(utcHourOfDay(day), hour, `day ${day}`);
  assert.throws(() => monthOfDay(undefined), /model day/);
  assert.throws(() => utcHourOfDay(NaN), /model day/);
});

test('ERA5\'s climatology holds every month and section, its means the means of its hours, with the trades westward east of the date line', () => {
  assert.deepEqual(era5.sections, SECTIONS);
  assert.equal(era5.latitude, 3);
  for (const key of ['tauX', 'slp']) {
    assert.equal(era5[key].length, 12);
    era5[key].forEach((month, m) => month.forEach((v, s) => {
      const hours = era5[`${key}ByHour`][m].map((h) => h[s]);
      assert.equal(hours.length, 24);
      assert.ok(Math.abs(v - hours.reduce((a, x) => a + x, 0) / 24) < 0.006, `${key} ${MONTHS[m]} ${sectionName(SECTIONS[s])}`);
    }));
  }
  for (const month of era5.tauX) for (let s = 6; s <= 13; s++) assert.ok(month[s] < -10, `tau_x ${month[s]} at ${sectionName(SECTIONS[s])}`);
  for (const month of era5.slp) for (const v of month) assert.ok(v > 1005 && v < 1016, `SLP ${v}`);
  assert.ok(era5.seaPoints.every((n) => n > 500 && n <= 1000), era5.seaPoints.join(' '));
  const august = era5Sections(era5, 882.25);
  assert.deepEqual([august.era5Tau, august.era5Slp, august.era5SlpAt], [era5.tauX[7], era5.slp[7], era5.slpByHour[7][18]]);
  assert.deepEqual(august.era5TauAt, era5.tauXByHour[7][18].map((v, s) => (v + era5.tauXByHour[7][19][s]) / 2));
  assert.equal(era5Sections(era5, 882.99).era5TauAt[0], (era5.tauXByHour[7][12][0] + era5.tauXByHour[7][13][0]) / 2);
});

test('a pair\'s table gives each side\'s mean over its runs ± its standard error and the difference ± its own', () => {
  const table = (tau0, slp0) => sectionTable({ tau: SECTIONS.map((_, k) => (k ? -50 : tau0)), slp: SECTIONS.map((_, k) => (k ? 1010 : slp0)), era5Tau: SECTIONS.map((_, k) => (k ? -40 : -7.5)), era5Slp: SECTIONS.map(() => 1011), era5TauAt: SECTIONS.map(() => -6), era5SlpAt: SECTIONS.map(() => 1012) }, '12 UTC');
  assert.deepEqual(readSections(table(-12, 1009.5))[0], { tau: -12, era5Tau: -7.5, era5TauAt: -6, slp: 1009.5, era5Slp: 1011, era5SlpAt: 1012 });
  const lines = pairTable([table(-10, 1009), table(-12, 1010), table(-14, 1011)], [table(-20, 1010), table(-20, 1010), table(-23, 1013)]).split('\n');
  assert.equal(lines.find((l) => l.startsWith('| section |')), '| section | ERA5 | ERA5 at 12 UTC | base | base/ERA5 | feature | feature/ERA5 | feature - base |');
  assert.equal(lines.find((l) => l.startsWith('| 120E-130E |')), '| 120E-130E | -7.5 | -6.0 | -12.0 ± 1.2 | 1.60 | -21.0 ± 1.0 | 2.80 | -9.0 ± 1.5 |');
  assert.equal(lines.find((l) => l.startsWith('| 130E-140E |')), '| 130E-140E | -40.0 | -6.0 | -50.0 ± 0.0 | 1.25 | -50.0 ± 0.0 | 1.25 | 0.0 ± 0.0 |');
  assert.equal(lines.filter((l) => l.startsWith('| 120E-130E |'))[1], '| 120E-130E | 1011.00 | 1012.00 | 1010.00 ± 0.58 | 1011.00 ± 1.00 | 1.00 ± 1.15 |');
});

test('armMetrics.mjs ends with the sections against ERA5 in the saves\' month and hour and leaves its bins\' table as it was', async () => {
  const { createModel } = await import('../js/model.module.js');
  const { Grid } = await import('../js/grid.module.js');
  const { stateGeography, geographyHeader } = await import('../js/geography.module.js');
  const { gridHeader } = await import('../js/stateGrid.module.js');
  const { initializeState } = await import('../js/physics/init.module.js');
  const { encodeState } = await import('../js/stateFile.module.js');
  const { record, topography, geography, land } = stateGeography(null);
  const model = createModel(new Grid(16), { topography, geography, land, ocean: false });
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  model.land.initialize();
  for (let n = 0; n < 4; n++) model.step(1350);
  const [pi, theta, u, surfaceT, q, qc, ice] = model.state, { mlmSubsidence, mlmHeight, mlmGate } = model.radiation;
  const dir = mkdtempSync(join(tmpdir(), 'armMetrics-')), file = join(dir, 'arm16_day0882.bin');
  writeFileSync(file, encodeState({ N: 16, ...gridHeader(model.mesh), K: model.core.K, day: 882, time: 882 * 86400, terrain: true, ...geographyHeader(record), levels: model.core.levels, pi, theta, u, surfaceT, q, qc, ice, concentration: model.seaIce.concentration, mlmSubsidence, mlmHeight, mlmGate, land: model.land.serialize() }));
  const run = spawnSync(process.execPath, [new URL('../scripts/armMetrics.mjs', import.meta.url).pathname, 'arm16', file], { encoding: 'utf8' });
  rmSync(dir, { recursive: true, force: true });
  assert.equal(run.status, 0, run.stderr);
  const lines = run.stdout.split('\n');
  assert.equal(lines.find((l) => l.startsWith('| | ')), '| | 140E-160E | 160E-180 | 180-160W | 160W-140W | 140W-120W | 120W-100W | 100W-80W |');
  assert.equal(lines.find((l) => l.startsWith('| stress tau_x')).split('|').length, 10);
  assert.match(run.stdout, /^summary arm16: SLP 180-160W less 120W-100W -?[\d.]+ hPa; stress 160W-140W/m);
  assert.match(run.stdout, /^sections, sea cells within 3 deg, mean of 1 states; ERA5 2007-2019 in August, its sea points within 3 deg$/m);
  assert.ok(lines.includes('| section | tau_x (mN/m2) | ERA5 tau_x | ratio | ERA5 tau_x at 12 UTC | SLP (hPa) | ERA5 SLP | ERA5 SLP at 12 UTC |'));
  const rows = readSections(run.stdout), at = era5Sections(era5, 882);
  rows.forEach((row, s) => {
    for (const key of ['era5Tau', 'era5Slp', 'era5TauAt', 'era5SlpAt']) assert.equal(row[key], Number(at[key][s].toFixed(2)), `${key} at ${sectionName(SECTIONS[s])}`);
    const ratio = Number(lines.find((l) => l.startsWith(`| ${sectionName(SECTIONS[s])} |`)).split('|')[4]);
    if (Number.isFinite(row.tau)) assert.ok(Math.abs(ratio - row.tau / row.era5Tau) < 0.01, `${sectionName(SECTIONS[s])}: ratio ${ratio}, ${row.tau} over ${row.era5Tau}`);
  });
  assert.ok(rows.filter((row) => Number.isFinite(row.tau) && row.tau !== 0 && row.slp > 950 && row.slp < 1070).length >= 12, rows.map((row) => `${row.tau} ${row.slp}`).join(', '));
});
