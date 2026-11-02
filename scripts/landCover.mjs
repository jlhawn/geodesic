#!/usr/bin/env node
/*
 * Writes data/landCover.bin (encodeLandCover in js/geography.module.js), the tree cover and the annual maximum of the green
 * fraction of the land on a 0.25° grid: 1440 columns from −180° and 720 rows from −90°, cell centers at −179.875 + 0.25 i and
 * −89.875 + 0.25 j, each field a uint8 of the fraction times 200 and 255 where the cell holds no land point. The sources are two
 * WPS geogrid datasets under WPS_DIR (default ~/.cache/webgcm-scratch/obsland/data/wps), each a folder of tiles named
 * XXXXX-XXXXX.YYYYY-YYYYY (1-based point ranges, rows from the south) and an index file:
 *
 *   modis_landuse_20class_30s_with_lakes: the MODIS IGBP land-use class of each 30″ point, 648 tiles of 1200 × 1200 unsigned
 *   bytes, no border. A point is land unless its class is 0 (missing), 17 (water) or 21 (lake). Its tree cover is by class: the
 *   forests 1–5 0.8, woody savanna 8 0.45, savanna 9 0.2, cropland/natural mosaic 14 0.15, wooded tundra 18 0.2, permanent
 *   wetland 11 and urban 13 0.1, every other class 0. A cell's treeCover is the mean over its land points.
 *
 *   greenfrac: the AVHRR monthly green fraction (scale 0.01) of 0.144° points, 2 tiles of 1250 × 1250 × 12 unsigned bytes with a
 *   3-point border. A point's maximum is the largest of its 12 months; a cell's greenMax is the mean of those maxima over its
 *   points that are land by the 30″ classes, read every sixth point (a 0.05° mask).
 *
 *   node scripts/landCover.mjs
 */
import { readFileSync, readdirSync, writeFileSync } from 'node:fs';
import { homedir } from 'node:os';
import { encodeLandCover, landCoverUrl, COVER_SCALE, COVER_NONE } from '../js/geography.module.js';

const WPS = (process.env.WPS_DIR ?? `${homedir()}/.cache/webgcm-scratch/obsland/data/wps`).replace(/\/?$/, '/');
const COLS = 1440, ROWS = 720, SIZE = 0.25, STRIDE = 6;
const TREES = { 1: 0.8, 2: 0.8, 3: 0.8, 4: 0.8, 5: 0.8, 8: 0.45, 9: 0.2, 11: 0.1, 13: 0.1, 14: 0.15, 18: 0.2 };
const isLand = new Uint8Array(256), treeOf = new Float64Array(256);
for (let c = 1; c <= 21; c++) if (c !== 17 && c !== 21) { isLand[c] = 1; treeOf[c] = TREES[c] ?? 0; }

function dataset(name) {
  const dir = `${WPS}${name}/${name}/`;
  const index = Object.fromEntries(readFileSync(dir + 'index', 'utf8').split(/\s+/).filter(Boolean).map((kv) => { const [k, v] = kv.split('='); return [k, v?.replace(/"/g, '')]; }));
  const w = { dir, dx: +index.dx, dy: +index.dy, lat0: +index.known_lat, lon0: +index.known_lon, tx: +index.tile_x, ty: +index.tile_y, tz: +(index.tile_z ?? 1), border: +(index.tile_bdr ?? 0), scale: +(index.scale_factor ?? 1) };
  w.tiles = readdirSync(dir).filter((f) => /^\d{5}-\d{5}\.\d{5}-\d{5}$/.test(f)).sort();
  w.width = w.tx + 2 * w.border; w.plane = w.width * (w.ty + 2 * w.border);
  return w;
}

function* tiles(w) {
  for (const t of w.tiles) {
    const [xr, yr] = t.split('.'), x0 = +xr.split('-')[0] - 1, y0 = +yr.split('-')[0] - 1, buf = readFileSync(w.dir + t);
    if (buf.length !== w.plane * w.tz) throw new Error(`${t}: ${buf.length} bytes, expected ${w.plane * w.tz}`);
    yield { buf, x0, y0 };
  }
}

const t0 = Date.now(), elapsed = () => `${((Date.now() - t0) / 1000).toFixed(0)} s`;
const cellOf = (lat, lon) => {
  const col = Math.min(COLS - 1, Math.max(0, Math.floor(((lon > 180 ? lon - 360 : lon) + 180) / SIZE))), row = Math.min(ROWS - 1, Math.max(0, Math.floor((lat + 90) / SIZE)));
  return row * COLS + col;
};

const lu = dataset('modis_landuse_20class_30s_with_lakes'), per = Math.round(SIZE / lu.dx);
if (Math.abs(per * lu.dx - SIZE) > 1e-6 || Math.abs(per * lu.dy - SIZE) > 1e-6 || lu.border !== 0 || lu.tz !== 1) throw new Error(`the land-use points do not tile the ${SIZE}° cells: dx ${lu.dx}, dy ${lu.dy}, border ${lu.border}, ${lu.tz} levels`);
if (Math.abs(lu.lon0 - (-180 + lu.dx / 2)) > 1e-4 || Math.abs(lu.lat0 - (-90 + lu.dy / 2)) > 1e-4) throw new Error(`the land-use points start at ${lu.lat0}, ${lu.lon0}, not half a point inside −90°, −180°`);
const MW = Math.round(360 / (lu.dx * STRIDE)), MH = Math.round(180 / (lu.dy * STRIDE)), mask = new Uint8Array(MW * MH);
const landN = new Float64Array(COLS * ROWS), treeSum = new Float64Array(COLS * ROWS);
for (const { buf, x0, y0 } of tiles(lu)) {
  for (let jy = 0; jy < lu.ty; jy++) {
    const row = Math.floor((y0 + jy) / per) * COLS, base = jy * lu.tx, sampled = (y0 + jy) % STRIDE === 0;
    for (let ix = 0; ix < lu.tx; ix++) {
      const c = buf[base + ix];
      if (!isLand[c]) continue;
      const cell = row + Math.floor((x0 + ix) / per);
      landN[cell]++; treeSum[cell] += treeOf[c];
      if (sampled && (x0 + ix) % STRIDE === 0) mask[((y0 + jy) / STRIDE) * MW + (x0 + ix) / STRIDE] = 1;
    }
  }
}
console.log(`land use: ${lu.tiles.length} tiles, ${elapsed()}`);
const landAt = (lat, lon) => {
  const cx = ((Math.round((lon - lu.lon0) / (lu.dx * STRIDE)) % MW) + MW) % MW, cy = Math.min(MH - 1, Math.max(0, Math.round((lat - lu.lat0) / (lu.dy * STRIDE))));
  return mask[cy * MW + cx] === 1;
};

const gf = dataset('greenfrac'), greenN = new Float64Array(COLS * ROWS), greenSum = new Float64Array(COLS * ROWS);
for (const { buf, x0, y0 } of tiles(gf)) {
  for (let jy = 0; jy < gf.ty; jy++) {
    const lat = gf.lat0 + (y0 + jy) * gf.dy;
    for (let ix = 0; ix < gf.tx; ix++) {
      const lon = gf.lon0 + (x0 + ix) * gf.dx;
      if (!landAt(lat, lon > 180 ? lon - 360 : lon)) continue;
      const off = (jy + gf.border) * gf.width + ix + gf.border;
      let top = 0;
      for (let z = 0; z < 12; z++) { const v = buf[off + z * gf.plane]; if (v > top) top = v; }
      const cell = cellOf(lat, lon);
      greenN[cell]++; greenSum[cell] += top * gf.scale;
    }
  }
}
console.log(`green fraction: ${gf.tiles.length} tiles, ${elapsed()}`);

const stored = (sum, n) => Uint8Array.from(sum, (s, k) => (n[k] > 0 ? Math.min(COVER_SCALE, Math.round(COVER_SCALE * s / n[k])) : COVER_NONE));
const treeCover = stored(treeSum, landN), greenMax = stored(greenSum, greenN);
writeFileSync(landCoverUrl, new Uint8Array(encodeLandCover({ cols: COLS, rows: ROWS, treeCover, greenMax })));

const mean = (field, south, north, west, east) => {
  let a = 0, s = 0;
  for (let r = 0; r < ROWS; r++) {
    const lat = -89.875 + SIZE * r;
    if (lat < south || lat > north) continue;
    for (let c = 0; c < COLS; c++) {
      const lon = -179.875 + SIZE * c;
      if (lon < west || lon > east || field[r * COLS + c] === COVER_NONE) continue;
      const w = Math.cos(lat * Math.PI / 180); a += w; s += w * field[r * COLS + c] / COVER_SCALE;
    }
  }
  return (s / a).toFixed(3);
};
const count = (field) => field.reduce((n, v) => n + (v === COVER_NONE ? 0 : 1), 0);
console.log(`data/landCover.bin: ${count(treeCover)} cells with a tree cover, ${count(greenMax)} with a green maximum; area-weighted means over those cells: tree cover ${mean(treeCover, -90, 90, -180, 180)}, green maximum ${mean(greenMax, -90, 90, -180, 180)}; Amazon 10S-5N 75-50W: ${mean(treeCover, -10, 5, -75, -50)} and ${mean(greenMax, -10, 5, -75, -50)}; ${elapsed()}`);
