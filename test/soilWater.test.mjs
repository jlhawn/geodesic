import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createGeography, readTopography, encodeSoilWater, decodeSoilWater, soilWaterUrl, soilWaterOnMesh, SOIL_SCALE, SOIL_NONE } from '../js/geography.module.js';
import { createLandSurface } from '../js/physics/land.module.js';

const DEG = 180 / Math.PI;
const soil = decodeSoilWater(readFileSync(soilWaterUrl));
const cellAt = (lat, lon) => Math.round((lat + 89.75) / 0.5) * soil.cols + Math.round((lon + 179.75) / 0.5);
const mesh16 = buildMesh(new Grid(16)), geo16 = createGeography(mesh16, readTopography('etopo1')), C = mesh16.nCells;
const free = (i) => geo16.land[i] && !geo16.iceSheet[i];

test('a soil-water file round-trips its grid and bytes and refuses a wrong tag or size', () => {
  const fill = Uint8Array.from({ length: 6 * 3 }, (_, n) => (n % 5 === 0 ? SOIL_NONE : (n * 37) % (SOIL_SCALE + 1)));
  const bytes = new Uint8Array(encodeSoilWater({ cols: 6, rows: 3, fill })), back = decodeSoilWater(bytes);
  assert.deepEqual([back.cols, back.rows, [...back.fill]], [6, 3, [...fill]]);
  assert.throws(() => decodeSoilWater(bytes.slice(0, -1)), /wrong size/);
  bytes[0] = 0;
  assert.throws(() => decodeSoilWater(bytes), /not a soil-water file/);
});

test('data/soilwater_march.bin holds the March 20 fill on the 0.5° grid: a dry Sahara, a wet Amazon, a middling Europe and nothing on the sea', () => {
  assert.deepEqual([soil.cols, soil.rows, readFileSync(soilWaterUrl).length], [720, 360, 12 + 720 * 360]);
  assert.ok(soil.fill.every((v) => v <= SOIL_SCALE || v === SOIL_NONE), 'every byte is a fill or none');
  const fill = (lat, lon) => soil.fill[cellAt(lat, lon)] / SOIL_SCALE;
  assert.ok(fill(23, 10) < 0.1, `Sahara ${fill(23, 10)}`);
  assert.ok(fill(-3, -60) >= 0.5, `Amazon ${fill(-3, -60)}`);
  assert.ok(fill(50, 10) > 0.3 && fill(50, 10) < 0.6, `Europe ${fill(50, 10)}`);
  assert.equal(soil.fill[cellAt(0, -140)], SOIL_NONE, 'the equatorial Pacific');
});

test('at N = 16 the geography carries wetness, the fill on every land cell and NaN on the sea, and null turns it off', () => {
  assert.ok(geo16.wetness instanceof Float64Array && geo16.wetness.length === C);
  for (let i = 0; i < C; i++) {
    if (geo16.land[i]) assert.ok(geo16.wetness[i] >= 0 && geo16.wetness[i] <= 1, `land cell ${i}: ${geo16.wetness[i]}`);
    else assert.ok(Number.isNaN(geo16.wetness[i]), `sea cell ${i}: ${geo16.wetness[i]}`);
  }
  const box = (S, N, W, E) => {
    let a = 0, s = 0;
    for (let i = 0; i < C; i++) {
      const lat = mesh16.latCell[i] * DEG, lon = mesh16.lonCell[i] * DEG;
      if (free(i) && lat >= S && lat <= N && lon >= W && lon <= E) { a += mesh16.areaCell[i]; s += mesh16.areaCell[i] * geo16.wetness[i]; }
    }
    return s / a;
  };
  const amazon = box(-10, 5, -75, -50), sahara = box(18, 30, -10, 30);
  assert.ok(amazon > 0.6 && sahara < 0.05, `Amazon ${amazon}, Sahara ${sahara}`);
  assert.equal(createGeography(mesh16, readTopography('etopo1'), { soilWater: null }).wetness, undefined);
});

test('a land cell without a point takes the mean of its neighbors that have one, ring by ring, and one out of reach and the sea stay NaN', () => {
  const sparse = { cols: 720, rows: 360, fill: new Uint8Array(720 * 360).fill(SOIL_NONE) };
  sparse.fill[cellAt(10.25, 20.25)] = 200;
  const lat = 10.25 / DEG, lon = 20.25 / DEG, x = Math.cos(lat) * Math.cos(lon), y = Math.cos(lat) * Math.sin(lon), z = Math.sin(lat);
  let seed = 0, best = -2;
  for (let i = 0; i < C; i++) { const d = x * Math.cos(mesh16.latCell[i]) * Math.cos(mesh16.lonCell[i]) + y * Math.cos(mesh16.latCell[i]) * Math.sin(mesh16.lonCell[i]) + z * Math.sin(mesh16.latCell[i]); if (d > best) { best = d; seed = i; } }
  const ring = (cells) => new Set(cells.flatMap((i) => Array.from({ length: mesh16.nEdgesOnCell[i] }, (_, m) => mesh16.cellsOnCell[mesh16.maxEdges * i + m])));
  const first = ring([seed]), second = new Set([...ring([...first])].filter((i) => i !== seed && !first.has(i)));
  assert.ok(soilWaterOnMesh(sparse, mesh16, new Uint8Array(C).fill(1)).wetness.every((v) => Math.abs(v - 0.8) < 1e-12), 'every land cell is reached');
  const land = new Uint8Array(C);
  for (const i of [...first, ...second]) land[i] = 1;
  land[seed] = 0;
  const some = soilWaterOnMesh(sparse, mesh16, land).wetness;
  for (let i = 0; i < C; i++) assert.ok(land[i] ? Math.abs(some[i] - 0.8) < 1e-12 : Number.isNaN(some[i]), `cell ${i}: ${some[i]}`);
  const apart = new Uint8Array(C), far = C - 1 - seed;
  apart[far] = 1;
  assert.ok(!first.has(far) && Number.isNaN(soilWaterOnMesh(sparse, mesh16, apart).wetness[far]), 'a land cell with only sea around it has no value');
});

test('the observed start fills the root zone and the surface layer to the fill and the deep store to the root zone\'s fill, and keeps the half fill where a cell has none', () => {
  const wetness = Float64Array.from(geo16.wetness), holes = [...Array(C).keys()].filter(free).filter((_, n) => n % 7 === 0);
  for (const i of holes) wetness[i] = NaN;
  const land = createLandSurface(mesh16, { ...geo16, wetness });
  land.initialize();
  let seeded = 0;
  for (let i = 0; i < C; i++) {
    const fill = free(i) && Number.isFinite(wetness[i]) ? wetness[i] : null, cap = geo16.land[i] ? 300 : 0;
    if (fill !== null) seeded++;
    assert.equal(land.soil[i], fill === null ? 0.5 * cap : fill * 300, `soil ${i}`);
    assert.equal(land.surface[i], fill === null ? 0 : fill * 15, `surface ${i}`);
    assert.equal(land.deepWater[i], free(i) ? Math.min(1, Math.max(0, land.soil[i] / 300)) * 300 : 0, `deep ${i}`);
    if (fill !== null) assert.ok(Math.abs(land.deepWater[i] - fill * 300) < 1e-12, `deep ${i} at the fill`);
  }
  assert.ok(seeded > 400 && holes.length > 50, `${seeded} seeded cells, ${holes.length} without a fill`);
  const bucket = createLandSurface(mesh16, geo16, { vegetation: false });
  bucket.initialize();
  for (let i = 0; i < C; i++) assert.equal(bucket.soil[i], free(i) ? geo16.wetness[i] * 150 : geo16.land[i] ? 75 : 0, `bucket ${i}`);
});

test('without the soil water the observed start throws an error naming data/soilwater_march.bin, the placeholder starts and a land without vegetation do not, and the cover is named first', () => {
  const without = { ...geo16, wetness: undefined };
  assert.throws(() => createLandSurface(mesh16, without).initialize(), /data\/soilwater_march\.bin.*start 'neutral'/);
  for (const start of ['neutral', 'bare', 'green']) createLandSurface(mesh16, without, { start }).initialize();
  const bucket = createLandSurface(mesh16, without, { vegetation: false });
  bucket.initialize();
  for (let i = 0; i < C; i++) assert.equal(bucket.soil[i], geo16.land[i] ? 75 : 0);
  assert.throws(() => createLandSurface(mesh16, { ...without, treeCover: undefined, greenMax: undefined }).initialize(), /data\/landCover\.bin/);
});

test('the placeholder starts keep the half fill, an empty surface layer and the deep store at half whatever the geography holds', () => {
  for (const start of ['neutral', 'bare', 'green']) {
    const land = createLandSurface(mesh16, geo16, { start });
    land.initialize();
    for (let i = 0; i < C; i++) {
      assert.equal(land.soil[i], geo16.land[i] ? 0.5 * 300 : 0, `${start} soil ${i}`);
      assert.equal(land.surface[i], 0, `${start} surface ${i}`);
      assert.equal(land.deepWater[i], free(i) ? 150 : 0, `${start} deep ${i}`);
    }
  }
});
