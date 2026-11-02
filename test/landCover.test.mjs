import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createGeography, readTopography, geographyOptions, FRESH_GEOGRAPHY, decodeLandCover, encodeLandCover, landCoverUrl, landCoverOnMesh, COVER_NONE, COVER_SCALE } from '../js/geography.module.js';
import { createLandSurface, startPlaceholders, freshJumpDue, START_CODES } from '../js/physics/land.module.js';

const DAY = 86400, DEG = 180 / Math.PI;
const bytes = readFileSync(landCoverUrl), cover = decodeLandCover(bytes);
const cellAt = (lat, lon) => Math.round((lat + 89.875) / 0.25) * cover.cols + Math.round((lon + 179.875) / 0.25);
const mesh64 = buildMesh(new Grid(64)), geo64 = createGeography(mesh64, readTopography('etopo2022'), geographyOptions(FRESH_GEOGRAPHY).geography);
const mesh16 = buildMesh(new Grid(16)), geo16 = createGeography(mesh16, readTopography('etopo1'));
const areaMean = (mesh, field, keep) => { let a = 0, s = 0; for (let i = 0; i < mesh.nCells; i++) if (keep(i)) { a += mesh.areaCell[i]; s += mesh.areaCell[i] * field[i]; } return s / a; };

test('data/landCover.bin holds the tree cover and the green maximum on the 0.25° grid, 255 where a cell has no land point', () => {
  assert.equal(String.fromCharCode(...bytes.subarray(0, 4)), 'LCV1');
  assert.equal(bytes.length, 12 + 2 * 1440 * 720);
  assert.deepEqual([cover.cols, cover.rows], [1440, 720]);
  for (const field of [cover.treeCover, cover.greenMax]) {
    let none = 0;
    for (const v of field) { if (v === COVER_NONE) none++; else assert.ok(v <= COVER_SCALE, `value ${v}`); }
    assert.ok(none > 0.5 * field.length && none < 0.75 * field.length, `${none} cells without land`);
  }
  assert.equal(cover.treeCover[cellAt(0, -150)], COVER_NONE, 'open Pacific');
  assert.equal(cover.greenMax[cellAt(0, -150)], COVER_NONE);
  assert.ok(cover.treeCover[cellAt(-3.1, -60.1)] >= 0.5 * COVER_SCALE, 'central Amazon');
  assert.ok(cover.greenMax[cellAt(-3.1, -60.1)] >= 0.7 * COVER_SCALE);
  assert.equal(cover.treeCover[cellAt(23.1, 10.1)], 0, 'Libyan desert');
  assert.ok(cover.greenMax[cellAt(23.1, 10.1)] <= 0.1 * COVER_SCALE);
  assert.ok(cover.treeCover[cellAt(62.1, 100.1)] > 0.1 * COVER_SCALE && cover.treeCover[cellAt(62.1, 100.1)] < 0.7 * COVER_SCALE, 'central Siberia');
  const copy = decodeLandCover(encodeLandCover(cover));
  assert.deepEqual([copy.cols, copy.rows], [1440, 720]);
  assert.ok(Buffer.from(copy.treeCover).equals(Buffer.from(cover.treeCover)) && Buffer.from(copy.greenMax).equals(Buffer.from(cover.greenMax)));
  assert.throws(() => decodeLandCover(new Uint8Array(40)), /not a land-cover file/);
  assert.throws(() => decodeLandCover(bytes.subarray(0, bytes.length - 1)), /wrong size/);
});

test('on the N = 64 mesh the cover has the observed means over the land, 0 over the sea, and an Amazon of forest', () => {
  assert.ok(geo64.treeCover && geo64.greenMax);
  let sea = 0;
  for (let i = 0; i < mesh64.nCells; i++) {
    for (const f of [geo64.treeCover, geo64.greenMax]) assert.ok(f[i] >= 0 && f[i] <= 1, `cell ${i}: ${f[i]}`);
    if (!geo64.land[i] && (geo64.treeCover[i] !== 0 || geo64.greenMax[i] !== 0)) sea++;
  }
  assert.equal(sea, 0);
  const open = (i) => geo64.land[i] && !geo64.iceSheet[i], amazon = (i) => { const lat = mesh64.latCell[i] * DEG, lon = mesh64.lonCell[i] * DEG; return geo64.land[i] && lat > -10 && lat < 5 && lon > -75 && lon < -50; };
  const tree = areaMean(mesh64, geo64.treeCover, open), green = areaMean(mesh64, geo64.greenMax, open);
  assert.ok(Math.abs(tree - 0.26) < 0.015 && Math.abs(green - 0.51) < 0.015, `land means ${tree.toFixed(3)} and ${green.toFixed(3)}`);
  const forest = areaMean(mesh64, geo64.treeCover, amazon), leaves = areaMean(mesh64, geo64.greenMax, amazon);
  assert.ok(Math.abs(forest - 0.75) < 0.02 && Math.abs(leaves - 0.84) < 0.02, `Amazon ${forest.toFixed(3)} and ${leaves.toFixed(3)}`);
});

test('a land cell without a point takes the mean of its neighbors that have one, ring by ring, and the sea is 0', () => {
  const none = () => new Uint8Array(1440 * 720).fill(COVER_NONE), sparse = { cols: 1440, rows: 720, treeCover: none(), greenMax: none() };
  sparse.treeCover[cellAt(10.125, 20.125)] = 100; sparse.greenMax[cellAt(10.125, 20.125)] = 200;
  const lat = 10.125 / DEG, lon = 20.125 / DEG, x = Math.cos(lat) * Math.cos(lon), y = Math.cos(lat) * Math.sin(lon), z = Math.sin(lat);
  let seed = 0, best = -2;
  for (let i = 0; i < mesh16.nCells; i++) { const d = x * Math.cos(mesh16.latCell[i]) * Math.cos(mesh16.lonCell[i]) + y * Math.cos(mesh16.latCell[i]) * Math.sin(mesh16.lonCell[i]) + z * Math.sin(mesh16.latCell[i]); if (d > best) { best = d; seed = i; } }
  const ring = (cells) => new Set(cells.flatMap((i) => Array.from({ length: mesh16.nEdgesOnCell[i] }, (_, m) => mesh16.cellsOnCell[mesh16.maxEdges * i + m])));
  const first = ring([seed]), second = new Set([...ring([...first])].filter((i) => i !== seed && !first.has(i)));
  const everywhere = landCoverOnMesh(sparse, mesh16, new Uint8Array(mesh16.nCells).fill(1));
  assert.ok(everywhere.treeCover.every((v) => Math.abs(v - 0.5) < 1e-12) && everywhere.greenMax.every((v) => Math.abs(v - 1) < 1e-12), 'every land cell is reached');
  const land = new Uint8Array(mesh16.nCells);
  for (const i of [...first, ...second]) land[i] = 1;
  land[seed] = 0;
  const some = landCoverOnMesh(sparse, mesh16, land);
  for (let i = 0; i < mesh16.nCells; i++) assert.equal(some.treeCover[i], land[i] ? 0.5 : 0, `cell ${i}${first.has(i) ? ' (first ring)' : second.has(i) ? ' (second ring)' : ''}`);
  const apart = new Uint8Array(mesh16.nCells), far = mesh16.nCells - 1 - seed;
  apart[far] = 1;
  assert.equal(landCoverOnMesh(sparse, mesh16, apart).greenMax[far], first.has(far) ? 1 : 0, 'a land cell with only sea around it stays at 0');
});

test('the observed start sets the cover to the green maximum, the trees to the lower of the tree cover and that cover, the carbon to the neutral placeholder and holds the carbon only', () => {
  const land = createLandSurface(mesh64, geo64), p = startPlaceholders();
  land.initialize();
  assert.deepEqual([...land.record], [0, START_CODES.observed, START_CODES.observed]);
  assert.equal(START_CODES.observed, 4);
  assert.deepEqual(p.observed, { cover: null, share: null, carbon: p.neutral.carbon });
  let cells = 0;
  for (let i = 0; i < mesh64.nCells; i++) {
    if (!geo64.land[i] || geo64.iceSheet[i]) { assert.deepEqual([land.vegetation[i], land.canopy[i], land.soilCarbon[i], land.snowFreeCover[i]], [0, 0, 0, 0], `cell ${i}`); continue; }
    cells++;
    assert.equal(land.vegetation[i], Math.min(1, geo64.greenMax[i]));
    assert.equal(land.snowFreeCover[i], land.vegetation[i]);
    assert.equal(land.canopy[i], Math.min(geo64.treeCover[i], land.vegetation[i]));
    assert.equal(land.soilCarbon[i], p.neutral.carbon);
    assert.equal(land.soil[i], geo64.wetness[i] * 300);
  }
  assert.ok(cells > 0.2 * mesh64.nCells, `${cells} land cells`);
  const plain = createLandSurface(mesh64, geo64, { treeline: false });
  plain.initialize();
  for (let i = 0; i < mesh64.nCells; i++) if (geo64.land[i] && !geo64.iceSheet[i]) assert.equal(plain.canopy[i], plain.vegetation[i], 'without the treeline the trees are the cover');
  const bare = createLandSurface(mesh64, geo64, { vegetation: false });
  bare.initialize();
  assert.ok(bare.vegetation.every((v) => v === 0) && bare.record[1] === 0);
});

test('a geography without the cover throws an error naming data/landCover.bin for the observed start and the observed goal, not for a placeholder start under the water goal, and a land without vegetation needs neither', () => {
  const without = { ...geo16, treeCover: undefined, greenMax: undefined };
  assert.throws(() => createLandSurface(mesh16, without).initialize(), /data\/landCover\.bin/);
  assert.throws(() => createLandSurface(mesh16, without, { start: 'observed', coverGoal: 'water' }).initialize(), /data\/landCover\.bin/);
  for (const start of ['neutral', 'bare', 'green']) createLandSurface(mesh16, without, { start }).initialize();
  createLandSurface(mesh16, without, { vegetation: false }).initialize();
  assert.throws(() => createLandSurface(mesh16, without, { start: 'neutral', coverGoal: 'observed' }), /data\/landCover\.bin.*coverGoal 'water'/);
  assert.throws(() => createLandSurface(mesh16, geo16, { coverGoal: 'dry' }), /coverGoal/);
  assert.equal(createGeography(mesh16, readTopography('etopo1'), { landCover: null }).treeCover, undefined, 'null turns the cover off');
});

const FAST = { treeMoisture: false, growthTime: 3600, declineTime: 3600, treeGrowthTime: 3600, treeDeclineTime: 3600, snowDeclineTime: 3600 };
function run(land, steps, { dt = 900, air = 288 } = {}) {
  const C = mesh16.nCells, flux = new Float64Array(C), surfaceT = new Float64Array(C);
  for (let n = 0; n < steps; n++) {
    surfaceT.fill(air);
    for (let i = 0; i < C; i++) if (geo16.land[i]) land.update(i, surfaceT, flux, 0, dt, air, 1e-5);
    land.advance(dt);
  }
}
const vegetated = (i) => geo16.land[i] && !geo16.iceSheet[i];

test('under hold 4 the trees run free from day 0 and the carbon stays at the placeholder until the jump, which sets the carbon and leaves the trees', () => {
  const land = createLandSurface(mesh16, geo16, { ...FAST, coverGoal: 'water' }), neutral = startPlaceholders().neutral.carbon;
  land.initialize();
  const trees0 = Float64Array.from(land.canopy), open = Array.from({ length: mesh16.nCells }, (_, i) => i).filter(vegetated);
  run(land, 24);
  let freed = 0, held = 0;
  for (let i = 0; i < mesh16.nCells; i++) if (vegetated(i)) { if (land.canopy[i] !== trees0[i]) freed++; if (land.soilCarbon[i] === neutral) held++; }
  assert.ok(freed > 0.9 * open.length, `${freed} of ${open.length} cells with moving trees`);
  assert.equal(held, open.length, 'the carbon is held');
  assert.equal(land.record[1], START_CODES.observed);
  const trees = Float64Array.from(land.canopy), report = land.jump();
  let changed = 0;
  for (let i = 0; i < mesh16.nCells; i++) { assert.equal(land.canopy[i], trees[i], `trees of cell ${i}`); if (vegetated(i) && land.soilCarbon[i] !== neutral) changed++; }
  assert.ok(changed > 0.9 * open.length, `${changed} cells with new carbon`);
  assert.deepEqual([land.record[1], land.record[2]], [0, START_CODES.observed], 'the jump ends the hold');
  assert.ok(report.global.carbon[1] !== report.global.carbon[0] && report.global.trees[1] === report.global.trees[0]);
  assert.equal(freshJumpDue([365 * DAY, 4, 4], 365 * DAY - 900, 365 * DAY), true, 'an observed start jumps with the others');
  const other = createLandSurface(mesh16, geo16, { ...FAST, start: 'neutral' });
  other.initialize(); run(other, 24);
  const before = Float64Array.from(other.canopy);
  other.jump();
  assert.ok(Array.from(other.canopy).some((v, i) => v !== before[i]), 'a neutral start still sets its trees');
});

test('a land saved under hold 4 loads with its hold and keeps the carbon at the placeholder', () => {
  const land = createLandSurface(mesh16, geo16, FAST);
  land.initialize(); run(land, 6);
  const saved = land.serialize(), back = createLandSurface(mesh16, geo16, FAST);
  assert.deepEqual([...saved.record].slice(1), [4, 4]);
  back.load(saved);
  assert.deepEqual([...back.record], [...land.record]);
  run(back, 3);
  const carbon = startPlaceholders().neutral.carbon;
  assert.ok(Array.from(back.soilCarbon).every((c, i) => !vegetated(i) || c === carbon));
});

test("coverGoal 'observed' relaxes the cover toward the green maximum and the trees toward the lower of the tree cover and the cover, and under snow the trees hold", () => {
  const land = createLandSurface(mesh16, geo16, { ...FAST, start: 'neutral', coverGoal: 'observed' });
  land.initialize(); land.jump();
  const open = Array.from({ length: mesh16.nCells }, (_, i) => i).filter(vegetated), wet = createLandSurface(mesh16, geo16, { ...FAST, start: 'neutral', coverGoal: 'water' });
  wet.initialize(); wet.jump();
  const gap = (l) => open.reduce((s, i) => s + Math.abs(l.vegetation[i] - geo16.greenMax[i]), 0) / open.length;
  const before = gap(land);
  run(land, 120); run(wet, 120);
  assert.ok(before > 0.1, `the start is ${before.toFixed(3)} from the map`);
  assert.ok(gap(land) < 1e-6, `the cover is ${gap(land)} from the map`);
  assert.ok(gap(wet) > 0.1, 'the bucket-fill goal does not follow the map');
  for (const i of open) assert.ok(Math.abs(land.canopy[i] - Math.min(geo16.treeCover[i], geo16.greenMax[i])) < 1e-6, `trees of cell ${i}: ${land.canopy[i]}`);
  const cell = open.find((i) => geo16.treeCover[i] > 0.3), trees = land.canopy[cell];
  land.snow[cell] = 5;
  const flux = new Float64Array(mesh16.nCells), surfaceT = new Float64Array(mesh16.nCells).fill(260);
  for (let n = 0; n < 8; n++) land.update(cell, surfaceT, flux, 0, 900, 260, 1e-5);
  assert.equal(land.canopy[cell], trees, 'the trees hold under snow');
  assert.ok(land.vegetation[cell] < 0.9 * geo16.greenMax[cell], 'the cover decays under snow');
});
