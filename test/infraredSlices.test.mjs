import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { clearLongwave, longwaveTableFor, gasPaths, opticalDepth } from '../js/physics/longwave.module.js';
import { STEFAN_BOLTZMANN } from '../js/physics/radiation.module.js';
import { makeColumn, readColumn } from '../explain/radiativeColumn.module.js';
import { slices, dryColumn, airColumn, vaporTable, GASES, VAPOR, ABSORBERS, CLEAR } from '../explain/infraredSlices.module.js';

const MLS = readColumn(JSON.parse(readFileSync(new URL('../explain/data/mls_bl36.json', import.meta.url))));
const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const sum = (values) => values.reduce((a, b) => a + b, 0);
const only = (gas) => Object.fromEntries(GASES.map((g) => [g, g === gas]));
const CASES = [
  ['2010 gases', {}],
  ['carbon dioxide alone', { gases: only('carbonDioxide') }],
  ['no carbon dioxide', { carbonDioxide: 0 }],
  ['half the carbon dioxide', { carbonDioxide: 0.5 }],
  ['doubled carbon dioxide', { carbonDioxide: 2 }],
  ['four times the carbon dioxide', { carbonDioxide: 4 }],
  ['no ozone', { gases: { ozone: false } }],
  ['no methane', { gases: { methane: false } }],
  ['no nitrous oxide', { gases: { nitrousOxide: false } }],
  ['no gases', { gases: only(null) }],
];

test('the slices sum to clearLongwave and to the model\'s own column on the same dry air', () => {
  const table = longwaveTableFor(MLS.levels);
  for (const [name, options] of CASES) {
    const s = slices(MLS, options), flux = clearLongwave(dryColumn(MLS, options), { table });
    close(sum(s.points.map((p) => p.up)), flux.up[0], 1e-9, `${name}: OLR from the slices`);
    close(sum(s.points.map((p) => p.down)), flux.down[MLS.temperature.length], 1e-9, `${name}: DLR from the slices`);
    const column = makeColumn({ profile: MLS, beam: 0, carbonDioxide: options.carbonDioxide ?? 1, gases: { ...options.gases, ...(options.carbonDioxide === 0 ? { carbonDioxide: false } : {}) } });
    close(s.olr, column.budget.olr, 1e-9, `${name}: OLR against radiation.column`);
    close(s.dlr, column.budget.dlr, 1e-9, `${name}: DLR against radiation.column`);
    close(s.escaping, column.budget.escaping, 1e-12, `${name}: escaping share against radiation.column`);
    close(s.groundEmission, column.budget.groundEmission, 1e-9, `${name}: σTs⁴`);
  }
});

test('the slices share out σTs⁴ exactly, and no slice sends more to space than the ground gives it', () => {
  for (const [name, options] of CASES) {
    const s = slices(MLS, options);
    assert.equal(s.points.length, 34);
    close(sum(s.points.map((p) => p.share)), 1, 1e-12, `${name}: shares`);
    close(sum(s.points.map((p) => p.ground)), STEFAN_BOLTZMANN * MLS.groundTemperature ** 4, 1e-9, `${name}: ground emission`);
    for (const p of s.points) {
      assert.ok(p.up <= p.ground + 1e-12 && p.up >= 0, `${name}: slice ${p.index} sends ${p.up} of ${p.ground}`);
      close(p.escaping, p.share * Math.exp(-p.depth), 1e-15, `${name}: slice ${p.index} escaping`);
    }
  }
});

test('the dry midlatitude summer: 356.2 W/m² leaves the top of 424.8, 127.6 comes back down and 65.3 percent escapes directly', () => {
  const s = slices(MLS);
  close(s.groundEmission, 424.80, 0.005, 'σTs⁴');
  close(s.olr, 356.18, 0.005, 'OLR');
  close(s.greenhouse, 68.62, 0.005, 'G');
  close(s.dlr, 127.63, 0.005, 'DLR');
  close(s.escaping, 0.6534, 5e-5, 'escaping');
  const alone = slices(MLS, { gases: only('carbonDioxide') }), doubled = slices(MLS, { carbonDioxide: 2 }), none = slices(MLS, { gases: only(null) });
  close(alone.olr, 379.60, 0.005, 'OLR, carbon dioxide alone');
  close(alone.greenhouse, 45.19, 0.005, 'G, carbon dioxide alone');
  close(doubled.olr, 352.86, 0.005, 'OLR, doubled carbon dioxide');
  close(none.olr, none.groundEmission, 1e-9, 'OLR with no gases');
  close(none.dlr, 0, 1e-12, 'DLR with no gases');
  assert.ok(none.points.every((p) => p.level === null && p.dominant === null), 'with no gases every slice is clear');
});

test('the optical depth from the top to each slice\'s level is 1, and clear slices never reach 1', () => {
  const table = longwaveTableFor(MLS.levels);
  for (const [name, options] of CASES) {
    const s = slices(MLS, options), column = dryColumn(MLS, options), paths = gasPaths(column, table), K = MLS.temperature.length;
    for (const p of s.points) {
      const row = table.points[p.index], tau = Array.from({ length: K }, (_, k) => opticalDepth(row, paths, k));
      close(p.depth, sum(tau), 1e-9 * Math.max(1, p.depth), `${name}: slice ${p.index} column depth`);
      close(sum(Object.values(p.gas)), p.depth, 1e-9 * Math.max(1, p.depth), `${name}: slice ${p.index} split by gas`);
      if (p.level === null) { assert.ok(p.depth < 1, `${name}: slice ${p.index} is clear with depth ${p.depth}`); continue; }
      const k = p.layer, top = MLS.levels[k] * MLS.surfacePressure, bottom = MLS.levels[k + 1] * MLS.surfacePressure;
      assert.ok(p.level >= top && p.level <= bottom, `${name}: slice ${p.index} level inside layer ${k}`);
      close(sum(tau.slice(0, k)) + tau[k] * (p.level - top) / (bottom - top), 1, 1e-12, `${name}: slice ${p.index} depth above its level`);
      if (p.dominant) assert.ok(p.gas[p.dominant] >= Math.max(...Object.values(p.gas)), `${name}: slice ${p.index} dominant gas`);
    }
  }
});

test('more carbon dioxide lifts every slice it dominates to lower pressure', () => {
  const runs = [0.5, 1, 2, 4].map((c) => slices(MLS, { carbonDioxide: c }));
  for (let j = 0; j < 34; j++) {
    const levels = runs.map((s) => (s.points[j].dominant === 'carbonDioxide' ? s.points[j].level : undefined));
    for (let n = 1; n < levels.length; n++) if (levels[n - 1] != null && levels[n] != null) assert.ok(levels[n] < levels[n - 1], `slice ${j}: ${levels[n]} after ${levels[n - 1]}`);
  }
});

test('the order runs from the clearest slice to the one that leaves from highest', () => {
  for (const [name, options] of CASES) {
    const { points, order } = slices(MLS, options), seen = order.map((j) => points[j]);
    assert.deepEqual([...order].sort((a, b) => a - b), points.map((p) => p.index), `${name}: a permutation`);
    for (let n = 1; n < seen.length; n++) {
      const a = seen[n - 1], b = seen[n];
      if (a.level === null && b.level === null) assert.ok(a.depth <= b.depth, `${name}: clear slices by depth`);
      else assert.ok(a.level === null || (b.level !== null && b.level <= a.level), `${name}: slices by level`);
    }
    for (const p of points) assert.equal(p.dominant === null, p.depth < CLEAR, `${name}: slice ${p.index} clear threshold`);
  }
});

test('the prose numbers: carbon dioxide\'s opaque slices, the clear window slice and the clear share', () => {
  const s = slices(MLS), opaque = s.points.filter((p) => p.dominant === 'carbonDioxide' && p.level !== null && p.level < 0.7 * MLS.surfacePressure);
  close(sum(opaque.map((p) => p.share)), 0.1672, 5e-4, 'carbon dioxide opaque share');
  const wide = opaque.filter((p) => p.share > 0.01).map((p) => p.level / 100);
  close(Math.max(...wide), 467.2, 0.1, 'lowest wide carbon dioxide level, hPa');
  close(Math.min(...wide), 11.7, 0.1, 'highest wide carbon dioxide level, hPa');
  const window = s.points.reduce((best, p) => (p.share > best.share ? p : best));
  close(window.share, 0.2136, 5e-4, 'the widest slice');
  assert.equal(window.dominant, 'ozone');
  close(window.depth, 0.0616, 5e-4, 'its optical depth');
  close(sum(s.points.filter((p) => p.level === null).map((p) => p.share)), 0.7642, 5e-4, 'share in slices that see the ground');
});

const DRY_GASES = Object.fromEntries(GASES.map((g) => [g, false]));
const MOIST = [
  ['the midlatitude summer\'s vapor', { humidity: 1 }],
  ['half its vapor', { humidity: 0.5 }],
  ['twice its vapor', { humidity: 2 }],
  ['its vapor and doubled carbon dioxide', { humidity: 1, carbonDioxide: 2 }],
  ['its vapor without ozone', { humidity: 1, gases: { ozone: false } }],
  ['vapor alone', { humidity: 1, gases: DRY_GASES }],
  ['vapor lines alone', { humidity: 1, gases: { ...DRY_GASES, continuum: false } }],
  ['the continuum alone', { humidity: 1, gases: { ...DRY_GASES, vaporLines: false } }],
  ['its vapor without the continuum', { humidity: 1, gases: { continuum: false } }],
  ['vapor given per layer', { humidity: MLS.humidity }],
];

test('with vapor the slices still sum to clearLongwave on the same table, and to the model\'s own column', () => {
  for (const [name, options] of MOIST) {
    const s = slices(MLS, options), table = vaporTable(MLS.levels, options.gases), column = airColumn(MLS, options), flux = clearLongwave(column, { table }), K = MLS.temperature.length;
    close(sum(s.points.map((p) => p.up)), flux.up[0], 1e-9, `${name}: OLR from the slices`);
    close(sum(s.points.map((p) => p.down)), flux.down[K], 1e-9, `${name}: DLR from the slices`);
    close(sum(s.points.map((p) => p.share)), 1, 1e-12, `${name}: shares`);
    const paths = gasPaths(column, table);
    for (const p of s.points) {
      close(sum(Object.values(p.gas)), p.depth, 1e-9 * Math.max(1, p.depth), `${name}: slice ${p.index} split by absorber`);
      close(p.depth, sum(Array.from({ length: K }, (_, k) => opticalDepth(table.points[p.index], paths, k))), 1e-9 * Math.max(1, p.depth), `${name}: slice ${p.index} depth`);
      if (options.gases?.vaporLines === false) assert.equal(p.gas.vaporLines, 0);
      if (options.gases?.continuum === false) assert.equal(p.gas.continuum, 0);
      if (p.dominant) assert.ok(ABSORBERS.includes(p.dominant) && p.gas[p.dominant] >= Math.max(...Object.values(p.gas)), `${name}: slice ${p.index} dominant absorber`);
    }
    if (options.gases?.vaporLines === false || options.gases?.continuum === false) continue;
    const humidity = typeof options.humidity === 'number' ? MLS.humidity.map((q) => q * options.humidity) : options.humidity;
    const own = makeColumn({ profile: MLS, beam: 0, humidity, carbonDioxide: options.carbonDioxide ?? 1, gases: options.gases }).budget;
    close(s.olr, own.olr, 1e-9, `${name}: OLR against radiation.column`);
    close(s.dlr, own.dlr, 1e-9, `${name}: DLR against radiation.column`);
    close(s.escaping, own.escaping, 1e-12, `${name}: escaping share against radiation.column`);
    close(s.precipitableWater, own.precipitableWater, 1e-9, `${name}: precipitable water against the column's`);
  }
  assert.deepEqual(VAPOR, ['vaporLines', 'continuum']);
  assert.equal(vaporTable(MLS.levels), longwaveTableFor(MLS.levels), 'with both vapor absorbers on the table is the model\'s');
  assert.equal(vaporTable(MLS.levels, { continuum: false }), vaporTable(MLS.levels, { continuum: false }), 'variants are made once');
  for (const p of slices(MLS).points) assert.ok(p.gas.vaporLines === 0 && p.gas.continuum === 0, 'dry air has no vapor depth');
});

test('vapor\'s greenhouse in the midlatitude summer: G 143.6 W/m² against 68.6 dry, 95.6 from vapor alone and 15.4 percent escaping', () => {
  const moist = slices(MLS, { humidity: 1 });
  close(moist.olr, 281.2, 0.05, 'OLR');
  close(moist.dlr, 348.1, 0.05, 'DLR');
  close(moist.greenhouse, 143.6, 0.05, 'G');
  close(moist.escaping, 0.154, 5e-4, 'share escaping directly');
  close(moist.precipitableWater, 28.8, 0.05, 'precipitable water, kg/m²');
  close(slices(MLS, { humidity: 1, gases: DRY_GASES }).greenhouse, 95.6, 0.05, 'G, vapor alone');
  close(slices(MLS, { humidity: 1, gases: { ...DRY_GASES, vaporLines: false } }).greenhouse, 41.1, 0.05, 'G, the continuum alone');
  close(slices(MLS, { humidity: 1, gases: { ...DRY_GASES, continuum: false } }).greenhouse, 86.5, 0.05, 'G, vapor lines alone');
  const without = slices(MLS, { humidity: 1, gases: { continuum: false } });
  close(moist.dlr - without.dlr, 44.9, 0.05, 'the continuum\'s back radiation');
  close(without.olr - moist.olr, 6.2, 0.05, 'the continuum\'s OLR');
  const window = moist.points[15];
  assert.equal(window.index, moist.points.reduce((best, p) => (p.share > best.share ? p : best)).index, 'slice 15 is the widest');
  close(window.depth, 0.61, 0.005, 'the window slice\'s optical depth');
  close(window.gas.continuum, 0.48, 0.005, 'of which the continuum');
  close(window.gas.vaporLines, 0.064, 0.001, 'of which the vapor lines');
  assert.equal(window.dominant, 'continuum');
  const doubling = [0.25, 0.5, 1, 2].map((f) => slices(MLS, { humidity: f }).greenhouse);
  [121.36, 131.15, 143.57, 160.67].forEach((G, n) => close(doubling[n], G, 0.01, `G with ×${[0.25, 0.5, 1, 2][n]} vapor`));
});
