import { test, after } from 'node:test';
import assert from 'node:assert/strict';
import { mkdtempSync, rmSync, writeFileSync, symlinkSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { spawnSync } from 'node:child_process';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { encodeState as encode } from '../js/stateFile.module.js';
import { gridHeader } from '../js/stateGrid.module.js';
import { stateGeography, createGeography, encodeSoilWater, FRESH_GEOGRAPHY, SOIL_SCALE, SOIL_NONE } from '../js/geography.module.js';
import { referenceScale, scaleFill, startGeography } from '../scripts/soilWater.mjs';
import { reseedLand } from '../scripts/reseedSoil.mjs';
import { decodeState } from '../js/stateFile.module.js';

const dir = mkdtempSync(join(tmpdir(), 'soil-seed-'));
after(() => rmSync(dir, { recursive: true, force: true }));

const COLS = 720, ROWS = 360, ROOT_ZONE = 300, DEG = 180 / Math.PI;
const grid = new Grid(16), header = { N: 16, ...gridHeader(grid), ...FRESH_GEOGRAPHY };
const { topography, geography, land } = stateGeography(header);
const model = createModel(grid, { topography, geography, land }), mesh = model.mesh, C = mesh.nCells;
const raw = { cols: COLS, rows: ROWS, fill: new Uint8Array(COLS * ROWS).fill(Math.round(0.3 * SOIL_SCALE)) };
const goal = (lat) => 0.3 + 0.5 * Math.sin(lat) ** 2;
const save = (name, soil, extra = {}) => {
  const path = join(dir, name);
  writeFileSync(path, encode({ ...header, K: model.core.K, day: 0, time: 0, land: { soil }, ...extra }));
  return path;
};
const references = [0.9, 1.1].map((f, n) => save(`ref${n}.bin`, Float64Array.from({ length: C }, (_, i) => f * goal(mesh.latCell[i]) * ROOT_ZONE)));

test('a seed scaled to reference states takes each 10° band of open land to the states\x27 mean fill there', async () => {
  const factor = await referenceScale(raw, references), seed = createGeography(mesh, topography, { ...geography, landCover: null, soilWater: { ...raw, fill: scaleFill(Uint8Array.from(raw.fill), factor) } });
  const sums = new Map();
  for (let i = 0; i < C; i++) {
    if (!seed.land[i] || (seed.iceSheet && seed.iceSheet[i]) || !Number.isFinite(seed.wetness[i])) continue;
    const band = Math.floor((mesh.latCell[i] * DEG + 90) / 10), s = sums.get(band) ?? { a: 0, seed: 0, goal: 0 };
    s.a += mesh.areaCell[i]; s.seed += mesh.areaCell[i] * seed.wetness[i]; s.goal += mesh.areaCell[i] * goal(mesh.latCell[i]);
    sums.set(band, s);
  }
  const all = [...sums.values()], area = all.reduce((x, s) => x + s.a, 0);
  let worst = 0;
  for (const s of all) if (s.a > 0.02 * area) worst = Math.max(worst, Math.abs(s.seed - s.goal) / s.a);
  const mean = (key) => all.reduce((x, s) => x + s[key], 0) / area;
  console.log(`${sums.size} bands: seed ${mean('seed').toFixed(3)} against the states\x27 ${mean('goal').toFixed(3)}, the worst band off by ${worst.toFixed(4)}; factors ${factor(-60).toFixed(2)} at 60S, ${factor(0).toFixed(2)} at the equator, ${factor(60).toFixed(2)} at 60N`);
  assert.ok(Math.abs(mean('seed') - mean('goal')) < 0.005 && worst < 0.02, `seed ${mean('seed')} against ${mean('goal')}, worst band ${worst}`);
  assert.ok(factor(0) < factor(60) && factor(60) < 5, 'the factor follows the states\x27 fill poleward within its range');
});

test('the reference states must share one grid and ground and hold the model\x27s own water, and a scaled seed must go to its own file', async () => {
  for (const [name, extra] of [['relaxed.bin', { grid: { relax: 3 } }], ['n32.bin', { N: 32 }], ['raster.bin', { mask: 'raster' }]]) {
    await assert.rejects(referenceScale(raw, [references[0], save(name, new Float64Array(C), extra)]), /the first reference is N = 16/, name);
  }
  await assert.rejects(referenceScale(raw, [references[0], save('observed.bin', new Float64Array(C).fill(0.3 * ROOT_ZONE))]), /holds the observed seed/);
  const script = new URL('../scripts/soilWater.mjs', import.meta.url).pathname, link = join(dir, 'linked.mjs');
  symlinkSync(script, link);
  for (const path of [script, link]) {
    const run = spawnSync(process.execPath, [path, join(dir, 'none.dods'), `--reference=${references[0]}`], { encoding: 'utf8' });
    assert.notEqual(run.status, 0, path);
    assert.match(run.stderr, /goes to its own file \(--out\)/, path);
  }
});

test('a scaled fill is held at 1 where the factor would carry it past', async () => {
  const wet = { cols: COLS, rows: ROWS, fill: Uint8Array.from({ length: COLS * ROWS }, (_, n) => (n % 997 === 0 ? SOIL_NONE : Math.round(0.6 * SOIL_SCALE))) };
  const factor = await referenceScale(wet, [save('full.bin', new Float64Array(C).fill(0.98 * ROOT_ZONE))]);
  const fitted = scaleFill(Uint8Array.from(wet.fill), factor), pushed = scaleFill(Uint8Array.from(wet.fill), () => 5);
  wet.fill.forEach((v, n) => {
    if (v === SOIL_NONE) { assert.ok(fitted[n] === SOIL_NONE && pushed[n] === SOIL_NONE, `point ${n}: no value`); return; }
    assert.ok(fitted[n] >= 0.95 * SOIL_SCALE && fitted[n] <= SOIL_SCALE, `point ${n}: fitted ${fitted[n]}`);
    assert.equal(pushed[n], SOIL_SCALE, `point ${n}: 0.6 times 5`);
  });
});

test('SOIL_SEED\x27s file replaces the observed fill of a fresh start cell for cell, is recorded, and is refused where the land comes from elsewhere', async () => {
  const fill = Uint8Array.from({ length: COLS * ROWS }, (_, n) => (n % 997 === 0 ? SOIL_NONE : Math.round((0.2 + 0.6 * Math.floor(n / COLS) / ROWS) * SOIL_SCALE)));
  const path = join(dir, 'seed.bin');
  writeFileSync(path, new Uint8Array(encodeSoilWater({ cols: COLS, rows: ROWS, fill })));
  const start = startGeography(geography, { SOIL_SEED: path });
  assert.equal(startGeography(geography, {}).geography, geography);
  assert.equal(startGeography(geography, {}).soilSeed, null);
  assert.match(start.soilSeed, /^seed\.bin sha256:[0-9a-f]{64}$/);
  assert.equal(start.soilNote, `soil seed ${start.soilSeed}`);
  assert.equal(startGeography(geography, {}).soilNote, 'soil observed');
  assert.throws(() => startGeography(geography, { SOIL_SEED: path, FROM: 'x.bin' }, { from: {} }), /FROM x\.bin gives this start/);
  assert.throws(() => startGeography(geography, { SOIL_SEED: path, LAND_FROM: 'y.bin' }, { landFrom: {} }), /LAND_FROM y\.bin gives this start/);
  for (const placed of ['neutral', 'bare', 'green']) {
    assert.throws(() => startGeography(geography, { SOIL_SEED: path, LAND: JSON.stringify({ start: placed }) }), new RegExp(`the '${placed}' start`));
    assert.equal(startGeography(geography, { LAND: JSON.stringify({ start: placed }) }).soilNote, `soil at the '${placed}' start`);
  }
  const inherited = startGeography(geography, { FROM: 'a.bin' }, { from: { soilSeed: start.soilSeed } }), landInherited = startGeography(geography, { FROM: 'a.bin', LAND_FROM: 'b.bin' }, { from: {}, landFrom: { soilSeed: start.soilSeed } });
  assert.equal(inherited.geography, geography);
  assert.ok(inherited.soilSeed === start.soilSeed && landInherited.soilSeed === start.soilSeed, 'a start carries the soil seed of the state its land comes from');
  assert.match(inherited.soilNote, /^soil from FROM a\.bin, seeded from seed\.bin sha256:/);
  assert.match(landInherited.soilNote, /^soil from LAND_FROM b\.bin, seeded from/);
  assert.equal(startGeography(geography, { FROM: 'a.bin' }, { from: {} }).soilNote, 'soil from FROM a.bin');
  assert.throws(() => startGeography(geography, { SOIL_SEED: path }, { saved: {} }), /not seeded from a file/);
  const recorded = await decodeState(encode({ ...header, K: model.core.K, day: 0, time: 0, soilSeed: start.soilSeed, land: { soil: new Float64Array(C) } }));
  assert.equal(recorded.soilSeed, start.soilSeed);
  assert.equal(startGeography(geography, { SOIL_SEED: path }, { saved: recorded }).soilSeed, start.soilSeed);
  assert.equal(startGeography(geography, {}, { saved: recorded }).soilSeed, start.soilSeed);
  const soilOf = (geo) => { const m = createModel(grid, { topography, geography: geo, land }); m.land.initialize(); return { wetness: m.geography.wetness, soil: m.land.serialize().soil, land: m.geography.land }; };
  const seeded = soilOf(start.geography), twin = soilOf({ ...geography, soilWater: { cols: COLS, rows: ROWS, fill } }), observed = soilOf(geography);
  let checked = 0, moved = 0;
  for (let i = 0; i < C; i++) {
    if (!seeded.land[i]) continue;
    assert.ok(Object.is(seeded.wetness[i], twin.wetness[i]) && seeded.soil[i] === twin.soil[i], `cell ${i}: fill ${seeded.wetness[i]} and soil ${seeded.soil[i]} against the seed's ${twin.wetness[i]} and ${twin.soil[i]}`);
    checked++;
    if (seeded.soil[i] !== observed.soil[i]) moved++;
  }
  console.log(`${checked} land cells take the seed, ${moved} of them off the observed fill`);
  assert.ok(checked > 200 && moved > 0.5 * checked, `${checked} land cells, ${moved} moved`);
});

test('a saved state re-seeded from a soil-water file holds the land water an observed fresh start on that file would', async () => {
  const fill = Uint8Array.from({ length: COLS * ROWS }, (_, n) => (n % 991 === 0 ? SOIL_NONE : Math.round((0.15 + 0.7 * ((n * 7) % ROWS) / ROWS) * SOIL_SCALE)));
  const soilWater = { cols: COLS, rows: ROWS, fill };
  const running = createModel(grid, { topography, geography, land });
  running.land.initialize();
  const serialized = running.land.serialize();
  for (let i = 0; i < C; i++) { serialized.soil[i] *= 0.6; if (serialized.surface) serialized.surface[i] *= 0.3; if (serialized.deep) serialized.deep[i] *= 0.9; }
  const saved = await decodeState(encode({ ...header, K: running.core.K, day: 365, time: 365 * 86400, soilSeed: 'earlier.bin sha256:0', land: serialized }));
  const kept = Object.fromEntries(['soil', 'surface', 'deep'].map((name) => [name, Float64Array.from(saved.land[name])])), { cells } = reseedLand(saved, soilWater, 'later.bin sha256:1');
  assert.equal(saved.soilSeed, 'later.bin sha256:1');
  const smaller = await decodeState(encode({ ...header, K: running.core.K, day: 365, time: 365 * 86400, soilSeed: 'earlier.bin sha256:0', options: { land: { rootZoneCapacity: 200, surfaceCapacity: 10 } }, land: running.land.serialize() }));
  reseedLand(smaller, soilWater);
  assert.equal(smaller.soilSeed, undefined);
  const fresh = createModel(grid, { topography, geography: { ...geography, soilWater }, land });
  fresh.land.initialize();
  const start = fresh.land.serialize();
  let open = 0, untouched = 0;
  for (let i = 0; i < C; i++) {
    const seeded = fresh.geography.land[i] && !(fresh.geography.iceSheet && fresh.geography.iceSheet[i]) && Number.isFinite(fresh.geography.wetness[i]);
    if (!seeded) { for (const name of ['soil', 'surface', 'deep']) assert.equal(saved.land[name][i], kept[name][i], `cell ${i}: ${name} kept`); untouched++; continue; }
    for (const name of ['soil', 'surface', 'deep']) if (start[name] && saved.land[name]) assert.ok(Math.abs(saved.land[name][i] - start[name][i]) < 1e-4, `cell ${i}: ${name} ${saved.land[name][i]} against the fresh start's ${start[name][i]}`);
    assert.ok(Math.abs(smaller.land.soil[i] - start.soil[i] * 2 / 3) < 1e-4 && Math.abs(smaller.land.surface[i] - start.surface[i] * 2 / 3) < 1e-4, `cell ${i}: the state's own capacities`);
    open++;
  }
  console.log(`${cells} cells re-seeded (${open} checked against a fresh start, ${untouched} kept); arrays ${['soil', 'surface', 'deep'].filter((n) => saved.land[n]).join(', ')}`);
  assert.ok(cells === open && open > 600 && saved.land.surface, `${cells} re-seeded, ${open} open`);
});
