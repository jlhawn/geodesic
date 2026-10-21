#!/usr/bin/env node
/*
 * Writes data/topography2022_0p25.bin and data/landice_0p25.bin, the
 * elevation and the land ice on the 0.25° grid of data/topography_0p25.bin
 * (1440 columns from −180°, 720 rows from the north), from ETOPO 2022's
 * 60″ ice-surface and bedrock elevations (NOAA
 * NCEI, ETOPO_2022_v1_60s_N90W180_{surface,bed}.tif in CACHE; the bed
 * under the Greenland and Antarctic ice sheets is BedMachine's). Each
 * 0.25° pixel is a block of 15 × 15 points of 60″. A point's ice
 * thickness is s − b, the surface less the bed, where the ice rests on
 * the bed and s/(1 − ρi/ρw), the thickness that floats with freeboard s,
 * where the bed lies deeper than that ice would reach (an ice shelf,
 * whose s − b would count the sea beneath it). A point holds ice where
 * that thickness is at least MINIMUM_THICKNESS: on ice-free land the
 * surface (the Copernicus DEM) and the bed (BedMachine's) disagree by up
 * to tens of meters either way, and 10 m keeps that noise from flipping
 * coastal cells while Antarctica's ice-free cells stay its rock outcrop.
 * The elevation is each pixel's mean surface as int16 meters, the format
 * of data/topography_0p25.bin: the bed with the ice on top where the ice
 * is grounded, and the floating ice's top where it floats. The land ice
 * (encodeLandIce in js/geography.module.js) holds the pixels with ice:
 * their counts of ice points and of those afloat, the mean thickness of
 * each, and the mean bed under the floating points. The pixel means of
 * the bed go to CACHE as an int16 raster on the same grid.
 *
 *   node scripts/landIce.mjs [CACHE]
 */
import { openSync, readSync, writeFileSync, existsSync } from 'node:fs';
import { inflateSync } from 'node:zlib';
import { encodeLandIce, landIceUrl, topographyUrl, ICE_DENSITY, SEAWATER_DENSITY } from '../js/geography.module.js';

const ROWS = 720, COLS = 1440, BLOCK = 15, SOURCE_ROWS = ROWS * BLOCK, SOURCE_COLS = COLS * BLOCK;
const FILES = { surface: 'ETOPO_2022_v1_60s_N90W180_surface.tif', bed: 'ETOPO_2022_v1_60s_N90W180_bed.tif' };
const SOURCE = 'https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/data/60s';
const FLOATING = 1 / (1 - ICE_DENSITY / SEAWATER_DENSITY), MINIMUM_THICKNESS = 10;

export function openTiff(path) {
  const fd = openSync(path, 'r');
  const read = (position, length) => { const b = Buffer.alloc(length); readSync(fd, b, 0, length, position); return b; };
  const head = read(0, 8);
  if (head.toString('latin1', 0, 2) !== 'II' || head.readUInt16LE(2) !== 42) throw new Error(`${path}: not a little-endian classic TIFF`);
  const at = head.readUInt32LE(4), count = read(at, 2).readUInt16LE(0), body = read(at + 2, 12 * count);
  const SIZE = { 1: 1, 2: 1, 3: 2, 4: 4, 12: 8 }, tags = new Map();
  for (let e = 0; e < count; e++) {
    const o = 12 * e, tag = body.readUInt16LE(o), type = body.readUInt16LE(o + 2), n = body.readUInt32LE(o + 4);
    if (!SIZE[type]) continue;
    const size = SIZE[type] * n, data = size <= 4 ? body.subarray(o + 8, o + 8 + size) : read(body.readUInt32LE(o + 8), size);
    tags.set(tag, type === 3 ? Array.from({ length: n }, (_, k) => data.readUInt16LE(2 * k))
      : type === 4 ? Array.from({ length: n }, (_, k) => data.readUInt32LE(4 * k))
        : type === 12 ? Array.from({ length: n }, (_, k) => data.readDoubleLE(8 * k))
          : data.toString('latin1').replace(/\0+$/, ''));
  }
  const one = (tag, fallback) => (tags.has(tag) ? tags.get(tag)[0] : fallback);
  const tiff = {
    width: one(256), height: one(257), bits: one(258), compression: one(259, 1), predictor: one(317, 1), sampleFormat: one(339, 1),
    tileWidth: one(322), tileLength: one(323), offsets: tags.get(324), counts: tags.get(325), scale: tags.get(33550), tiepoint: tags.get(33922),
    keys: tags.get(34735) ?? [], nodata: tags.has(42113) ? Number(tags.get(42113)) : null,
  };
  const rasterType = (() => { for (let k = 4; k + 3 < tiff.keys.length; k += 4) if (tiff.keys[k] === 1025) return tiff.keys[k + 3]; return null; })();
  const expect = [
    ['60″ global grid', tiff.width === SOURCE_COLS && tiff.height === SOURCE_ROWS],
    ['float32 samples', tiff.bits === 32 && tiff.sampleFormat === 3],
    ['DEFLATE', tiff.compression === 8 || tiff.compression === 32946],
    ['no or floating-point predictor', tiff.predictor === 1 || tiff.predictor === 3],
    ['tiles', tiff.tileWidth > 0 && tiff.tileLength > 0 && tiff.offsets && tiff.counts],
    ['pixel scale 1/60°', tiff.scale && Math.abs(tiff.scale[0] - 1 / 60) < 1e-12 && Math.abs(tiff.scale[1] - 1 / 60) < 1e-12],
    ['corner at −180°, 90°', tiff.tiepoint && tiff.tiepoint[0] === 0 && tiff.tiepoint[1] === 0 && tiff.tiepoint[3] === -180 && tiff.tiepoint[4] === 90],
    ['pixels are areas', rasterType === 1],
  ];
  for (const [what, ok] of expect) if (!ok) throw new Error(`${path}: expected ${what}`);
  const across = Math.ceil(tiff.width / tiff.tileWidth), shuffled = new Uint8Array(4 * tiff.tileWidth);
  tiff.across = across;
  tiff.readTile = (tx, ty) => {
    const n = ty * across + tx, raw = inflateSync(read(tiff.offsets[n], tiff.counts[n])), tw = tiff.tileWidth, th = tiff.tileLength;
    if (raw.length !== 4 * tw * th) throw new Error(`${path}: tile ${n} inflates to ${raw.length} bytes`);
    const out = new Float32Array(tw * th), bytes = new Uint8Array(out.buffer);
    if (tiff.predictor === 1) { bytes.set(raw); return out; }
    for (let r = 0; r < th; r++) {
      const base = 4 * tw * r;
      let sum = 0;
      for (let k = 0; k < 4 * tw; k++) { sum = (sum + raw[base + k]) & 255; shuffled[k] = sum; }
      for (let k = 0; k < tw; k++) {
        const o = base + 4 * k;
        bytes[o] = shuffled[3 * tw + k]; bytes[o + 1] = shuffled[2 * tw + k]; bytes[o + 2] = shuffled[tw + k]; bytes[o + 3] = shuffled[k];
      }
    }
    return out;
  };
  return tiff;
}

/*
 * The 0.25° pixel sums of the two rasters: per pixel the points read,
 * the points above sea level, the ice points, those afloat, and the sums
 * of the surface, the bed, the ice points' thickness, and the thickness
 * and the bed of those afloat.
 */
export function accumulate(surface, bed, minimum = MINIMUM_THICKNESS) {
  const P = ROWS * COLS;
  const sums = { points: new Int32Array(P), land: new Int32Array(P), ice: new Int32Array(P), floating: new Int32Array(P), surface: new Float64Array(P), bed: new Float64Array(P), thickness: new Float64Array(P), floatThickness: new Float64Array(P), floatBed: new Float64Array(P), missing: 0 };
  const { tileWidth: tw, tileLength: th, across } = surface;
  for (let ty = 0; ty * th < SOURCE_ROWS; ty++) {
    const top = Array.from({ length: across }, (_, tx) => surface.readTile(tx, ty)), under = Array.from({ length: across }, (_, tx) => bed.readTile(tx, ty));
    for (let r = ty * th; r < Math.min(SOURCE_ROWS, (ty + 1) * th); r++) {
      const row = Math.floor(r / BLOCK) * COLS, y = (r - ty * th) * tw;
      for (let c = 0; c < SOURCE_COLS; c++) {
        const tx = Math.floor(c / tw), x = c - tx * tw, s = top[tx][y + x], b = under[tx][y + x];
        if (s === surface.nodata || b === bed.nodata || !Number.isFinite(s) || !Number.isFinite(b)) { sums.missing++; continue; }
        const p = row + Math.floor(c / BLOCK);
        sums.points[p]++; sums.surface[p] += s; sums.bed[p] += b;
        if (s > 0) sums.land[p]++;
        if (!(s > 0 && s > b)) continue;
        const afloat = s * FLOATING, thickness = Math.min(s - b, afloat);
        if (thickness < minimum) continue;
        sums.ice[p]++; sums.thickness[p] += thickness;
        if (s - b > afloat) { sums.floating[p]++; sums.floatThickness[p] += thickness; sums.floatBed[p] += b; }
      }
    }
  }
  return sums;
}

function report(sums) {
  const R = 6371.0088, dLon = Math.PI / 720, region = (lat, lon) => (lat < -50 ? 'Antarctica' : lat > 58 && lon > -75 && lon < -10 ? 'Greenland box' : 'elsewhere');
  const out = {};
  for (let r = 0; r < ROWS; r++) {
    const north = Math.PI / 2 - r * Math.PI / ROWS, area = R * R * dLon * (Math.sin(north) - Math.sin(north - Math.PI / ROWS)), lat = 90 - (r + 0.5) * 0.25;
    for (let c = 0; c < COLS; c++) {
      const p = r * COLS + c;
      if (!sums.ice[p]) continue;
      const o = (out[region(lat, -180 + (c + 0.5) * 0.25)] ??= { area: 0, volume: 0, shelfArea: 0, shelfVolume: 0 }), share = area / sums.points[p];
      o.area += sums.ice[p] * share; o.volume += sums.thickness[p] * share / 1000;
      o.shelfArea += sums.floating[p] * share;
    }
  }
  for (const [name, o] of Object.entries(out)) console.log(`  ${name}: ice ${(o.area / 1e6).toFixed(3)} million km², ${(o.volume / 1e6).toFixed(2)} million km³; afloat ${(o.shelfArea / 1e6).toFixed(3)} million km²`);
}

if (import.meta.url === `file://${process.argv[1]}`) {
  const cache = process.argv[2] ?? new URL('../../geodesic-terrain-cache/etopo2022', import.meta.url).pathname;
  for (const name of Object.values(FILES)) if (!existsSync(`${cache}/${name}`)) throw new Error(`${cache}/${name} is missing: fetch it from ${SOURCE}/60s_${name.includes('surface') ? 'surface' : 'bed'}_elev_gtif/${name}`);
  const start = performance.now();
  const sums = accumulate(openTiff(`${cache}/${FILES.surface}`), openTiff(`${cache}/${FILES.bed}`));
  const P = ROWS * COLS, held = [];
  for (let p = 0; p < P; p++) {
    if (sums.points[p] !== BLOCK * BLOCK) throw new Error(`pixel ${p} holds ${sums.points[p]} points`);
    if (sums.ice[p]) held.push(p);
  }
  const mean = (sum, n) => (n ? Math.round(sum / n) : 0);
  const bytes = encodeLandIce({
    rows: ROWS, cols: COLS, points: BLOCK * BLOCK, index: Uint32Array.from(held), count: Uint8Array.from(held, (p) => sums.ice[p]), floating: Uint8Array.from(held, (p) => sums.floating[p]),
    thickness: Uint16Array.from(held, (p) => Math.min(65535, mean(sums.thickness[p], sums.ice[p]))), floatThickness: Uint16Array.from(held, (p) => Math.min(65535, mean(sums.floatThickness[p], sums.floating[p]))),
    floatBed: Int16Array.from(held, (p) => mean(sums.floatBed[p], sums.floating[p])),
  });
  writeFileSync(landIceUrl, new Uint8Array(bytes));
  const raster = (name) => new Uint8Array(Int16Array.from({ length: P }, (_, p) => mean(sums[name][p], sums.points[p])).buffer);
  writeFileSync(topographyUrl('etopo2022'), raster('surface'));
  writeFileSync(`${cache}/etopo2022_bed_0p25.bin`, raster('bed'));
  const name = (url) => url.pathname.split('/').slice(-2).join('/');
  console.log(`${name(topographyUrl('etopo2022'))} and ${name(landIceUrl)} (${bytes.byteLength} bytes, ${held.length} pixels with ice): ${((performance.now() - start) / 1000).toFixed(1)} s, ${sums.missing} points without data`);
  report(sums);
}
