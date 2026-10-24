#!/usr/bin/env node
/*
 * Writes data/geography1024.bin (js/fineGeography.module.js): the relaxed
 * N mesh (1024 and 30 Lloyd passes by default) with ETOPO 2022's 60″
 * surface and bed (CACHE, as for scripts/landIce.mjs) on it. Every 60″
 * point goes to its nearest cell, weighted by its area; its ice follows
 * scripts/landIce.mjs. A fine edge's sill is found from the points of its
 * two cells as createGeography finds a model edge's from the 0.25° raster,
 * but stepping only to a point's four neighbors: a diagonal step slips past
 * a crest one point wide (Gibraltar's Camarinal Sill reads 341 m that way,
 * 283 m four-connected). A point under ice afloat is water to its bed for
 * the sea fields. With LAKES, a GeoJSON file of lakes (scripts/lakes.mjs),
 * each lake's points carry its level and their depth below it.
 * The file is generated and not committed.
 *
 *   [LAKES=ne_10m_lakes.geojson] node --max-old-space-size=20000 scripts/fineGeography.mjs [CACHE] [N] [OUT]
 */
import { openSync, writeSync, closeSync, readFileSync, writeFileSync } from 'node:fs';
import { lakeFeatures, rasterizeLakes, lakeLevels } from './lakes.mjs';
import { openTiff } from './landIce.mjs';
import { relaxedMesh, nearestCell, fineEdges, encodeFineGeography, FINE_GEOGRAPHY_URL } from '../js/fineGeography.module.js';
import { ICE_DENSITY, SEAWATER_DENSITY } from '../js/geography.module.js';
import { DEFAULT_RELAX } from '../js/grid.module.js';

const ROWS = 10800, COLS = 21600, FLOATING = 1 / (1 - ICE_DENSITY / SEAWATER_DENSITY), MINIMUM_THICKNESS = 10, FIRN_AIR = 15;
const cache = process.argv[2] ?? '/Users/jlhawn/git_repos/jlhawn/geodesic-terrain-cache/etopo2022';
const N = Number(process.argv[3] ?? 1024), out = process.argv[4] ?? FINE_GEOGRAPHY_URL;
const clock = performance.now(), stamp = (what) => console.log(`${((performance.now() - clock) / 1000).toFixed(1)} s  ${what}  rss ${(process.memoryUsage().rss / 2 ** 30).toFixed(1)} GB`);

const mesh = relaxedMesh(N, DEFAULT_RELAX), C = mesh.nCells;
stamp(`relaxed mesh N=${N}: ${C} cells`);

const surface = openTiff(`${cache}/ETOPO_2022_v1_60s_N90W180_surface.tif`), bedTiff = openTiff(`${cache}/ETOPO_2022_v1_60s_N90W180_bed.tif`);
const sums = Object.fromEntries(['weight', 'land', 'surface', 'landSurface', 'seaSurface', 'bed', 'ice', 'iceThickness', 'floating', 'floatThickness', 'floatSurface', 'floatBed', 'wet', 'wetSea'].map((k) => [k, new Float64Array(C)]));
const points = new Uint32Array(C), maxDepth = new Float64Array(C), deepPoint = new Int32Array(C).fill(-1), maxDepthSea = new Float64Array(C), deepPointSea = new Int32Array(C).fill(-1);
const depth = new Uint16Array(ROWS * COLS), depthSea = new Uint16Array(ROWS * COLS), height = new Int16Array(ROWS * COLS), owner = new Int32Array(ROWS * COLS);
let missing = 0, cell = 0, rowStart = 0;
const { tileWidth: tw, tileLength: th, across } = surface;
for (let ty = 0; ty * th < ROWS; ty++) {
  const top = Array.from({ length: across }, (_, tx) => surface.readTile(tx, ty)), under = Array.from({ length: across }, (_, tx) => bedTiff.readTile(tx, ty));
  for (let r = ty * th; r < Math.min(ROWS, (ty + 1) * th); r++) {
    const north = Math.PI / 2 - r * Math.PI / ROWS, weight = Math.sin(north) - Math.sin(north - Math.PI / ROWS);
    const lat = Math.PI / 2 - (r + 0.5) * Math.PI / ROWS, cosLat = Math.cos(lat), sinLat = Math.sin(lat), y = (r - ty * th) * tw;
    cell = rowStart;
    for (let c = 0; c < COLS; c++) {
      const lon = -Math.PI + (c + 0.5) * 2 * Math.PI / COLS, n = r * COLS + c;
      cell = nearestCell(cosLat * Math.cos(lon), cosLat * Math.sin(lon), sinLat, cell, mesh.xCell, mesh.nEdgesOnCell, mesh.cellsOnCell, 6);
      if (c === 0) rowStart = cell;
      owner[n] = cell;
      const tx = Math.floor(c / tw), x = c - tx * tw, s = top[tx][y + x], b = under[tx][y + x];
      if (s === surface.nodata || b === bedTiff.nodata || !Number.isFinite(s) || !Number.isFinite(b)) { missing++; owner[n] = -1; continue; }
      const d = Math.max(0, -s);
      depth[n] = depthSea[n] = Math.min(65535, Math.round(d));
      height[n] = Math.max(-32768, Math.min(32767, Math.round(s)));
      points[cell]++;
      sums.weight[cell] += weight; sums.surface[cell] += weight * s; sums.bed[cell] += weight * b; sums.wet[cell] += weight * d;
      if (s > 0) { sums.land[cell] += weight; sums.landSurface[cell] += weight * s; } else sums.seaSurface[cell] += weight * s;
      if (deepPoint[cell] < 0 || d > maxDepth[cell]) { maxDepth[cell] = d; deepPoint[cell] = n; }
      let wet = d;
      if (s > 0 && s > b) {
        const afloat = (s - Math.min(FIRN_AIR, s / 2)) * FLOATING, thickness = Math.min(s - b, afloat);
        if (thickness >= MINIMUM_THICKNESS) {
          sums.ice[cell] += weight; sums.iceThickness[cell] += weight * thickness;
          if (s - b > afloat) {
            sums.floating[cell] += weight; sums.floatThickness[cell] += weight * thickness; sums.floatSurface[cell] += weight * s; sums.floatBed[cell] += weight * b;
            wet = Math.max(0, -b); depthSea[n] = Math.min(65535, Math.round(wet));
          }
        }
      }
      sums.wetSea[cell] += weight * wet;
      if (deepPointSea[cell] < 0 || wet > maxDepthSea[cell]) { maxDepthSea[cell] = wet; deepPointSea[cell] = n; }
    }
  }
  if (ty % 20 === 0) stamp(`row ${Math.min(ROWS, (ty + 1) * th)} of ${ROWS}`);
}
stamp(`points assigned (${missing} without data)`);

let empty = 0;
const fields = { points, maxDepth, maxDepthSea };
const ratio = (num, den) => Float64Array.from(num, (v, i) => (den[i] > 0 ? v / den[i] : 0));
fields.landFraction = ratio(sums.land, sums.weight);
fields.elevation = ratio(sums.surface, sums.weight);
fields.landElevation = ratio(sums.landSurface, sums.land);
fields.seaElevation = Float64Array.from(sums.seaSurface, (v, i) => (sums.weight[i] > sums.land[i] ? v / (sums.weight[i] - sums.land[i]) : 0));
fields.bed = ratio(sums.bed, sums.weight);
fields.iceFraction = ratio(sums.ice, sums.weight);
fields.iceThickness = ratio(sums.iceThickness, sums.ice);
fields.floatFraction = ratio(sums.floating, sums.weight);
fields.floatThickness = ratio(sums.floatThickness, sums.floating);
fields.floatSurface = ratio(sums.floatSurface, sums.floating);
fields.floatBed = ratio(sums.floatBed, sums.floating);
fields.wetDepth = ratio(sums.wet, sums.weight);
fields.wetDepthSea = ratio(sums.wetSea, sums.weight);
for (let i = 0; i < C; i++) if (!points[i]) empty++;
if (empty) throw new Error(`${empty} fine cells hold no 60″ point`);

const start = new Int32Array(C + 1), list = new Int32Array(ROWS * COLS - missing);
for (let n = 0; n < ROWS * COLS; n++) if (owner[n] >= 0) start[owner[n] + 1]++;
for (let i = 0; i < C; i++) start[i + 1] += start[i];
{ const fill = Int32Array.from(start.subarray(0, C)); for (let n = 0; n < ROWS * COLS; n++) if (owner[n] >= 0) list[fill[owner[n]]++] = n; }
stamp('point lists');

const lakeFraction = new Float64Array(C), lakeLevel = new Float64Array(C), lakeDepth = new Float64Array(C), dryMaxDepth = new Float64Array(C), lakeId = new Uint32Array(C);
if (process.env.LAKES) {
  const { lakes, left } = lakeFeatures(JSON.parse(readFileSync(process.env.LAKES, 'utf8')));
  const id = rasterizeLakes(lakes, ROWS, COLS), levels = lakeLevels(id, height, ROWS, COLS, lakes.length), shares = new Map();
  stamp(`${lakes.length} lakes rasterized (${left.length} intermittent left out: ${[...new Set(left.map((l) => l.featurecla))].join(', ')})`);
  for (let r = 0; r < ROWS; r++) {
    const north = Math.PI / 2 - r * Math.PI / ROWS, weight = Math.sin(north) - Math.sin(north - Math.PI / ROWS);
    for (let n = r * COLS; n < (r + 1) * COLS; n++) {
      const i = owner[n], k = id[n];
      if (i < 0) continue;
      if (!k || !Number.isFinite(levels[k - 1].level)) { dryMaxDepth[i] = Math.max(dryMaxDepth[i], depth[n]); continue; }
      const level = levels[k - 1].level;
      lakeFraction[i] += weight; lakeLevel[i] += weight * level; lakeDepth[i] += weight * Math.max(0, level - height[n]);
      if (!shares.has(i)) shares.set(i, new Map());
      shares.get(i).set(k, (shares.get(i).get(k) ?? 0) + weight);
    }
  }
  for (let i = 0; i < C; i++) if (lakeFraction[i] > 0) { lakeLevel[i] /= lakeFraction[i]; lakeDepth[i] /= lakeFraction[i]; lakeFraction[i] /= sums.weight[i]; }
  for (const [i, m] of shares) { let most = 0; for (const [k, w] of m) if (w > most) { most = w; lakeId[i] = k; } }
  const listed = lakes.map((lake, k) => ({ id: k + 1, name: lake.name, featurecla: lake.featurecla, ...levels[k] })).sort((a, b) => b.points - a.points);
  console.log(JSON.stringify(listed.slice(0, 30).map(({ name, featurecla, level, points: p, flat }) => [name, featurecla, Math.round(level), p, flat])));
  if (process.env.LAKE_REPORT) writeFileSync(process.env.LAKE_REPORT, JSON.stringify({ lakes: listed, left }, null, 1));
  stamp('lake fields');
} else dryMaxDepth.set(maxDepth);
Object.assign(fields, { dryMaxDepth, lakeFraction, lakeLevel, lakeDepth, lakeId });
const local = owner.fill(-1);

const { nEdges: E, cellsOnEdge } = fineEdges(mesh);
let keys = new Float64Array(4096), parent = new Int32Array(4096), members = new Int32Array(4096);
const find = (p) => { while (parent[p] !== p) { parent[p] = parent[parent[p]]; p = parent[p]; } return p; };
const SHIFT = 2 ** 28, FOUR = [[-1, 0], [1, 0], [0, -1], [0, 1]];
const sills = (depth, deepPoint, maxDepth, only = null, given = null) => {
  const edgeSill = given ? Uint16Array.from(given) : new Uint16Array(E);
  for (let e = 0; e < E; e++) {
    const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
    if (only && !only(a) && !only(b)) continue;
    edgeSill[e] = 0;
    if (maxDepth[a] <= 0 || maxDepth[b] <= 0) continue;
    const n = start[a + 1] - start[a] + start[b + 1] - start[b];
    if (n > keys.length) { keys = new Float64Array(2 * n); parent = new Int32Array(2 * n); members = new Int32Array(2 * n); }
    let m = 0;
    for (const i of [a, b]) for (let p = start[i]; p < start[i + 1]; p++) { members[m] = list[p]; keys[m] = depth[list[p]] * SHIFT + m; m++; }
    const order = keys.subarray(0, n).sort();
    const da = deepPoint[a], db = deepPoint[b];
    let sill = 0;
    for (let q = n - 1, added = 0; q >= 0; q--, added++) {
      const point = members[order[q] % SHIFT], r = Math.floor(point / COLS), c = point - r * COLS;
      local[point] = added; parent[added] = added;
      for (const [dr, dc] of FOUR) {
        const rr = r + dr;
        if (rr < 0 || rr >= ROWS) continue;
        const k = rr * COLS + ((c + dc + COLS) % COLS);
        if (local[k] < 0) continue;
        const x = find(added), y = find(local[k]);
        if (x !== y) parent[x] = y;
      }
      if (local[da] >= 0 && local[db] >= 0 && find(local[da]) === find(local[db])) { sill = depth[point]; break; }
    }
    for (let q = 0; q < n; q++) local[members[q]] = -1;
    edgeSill[e] = sill;
    if (e % 5e6 === 0) stamp(`edge ${e} of ${E}`);
  }
  return edgeSill;
};
const edgeSill = sills(depth, deepPoint, maxDepth);
stamp(`fine edge sills: ${E} edges`);
const edgeSillSea = sills(depthSea, deepPointSea, maxDepthSea, (i) => fields.floatFraction[i] > 0, edgeSill);
stamp('fine edge sills with the ice afloat as water');

const fd = openSync(out, 'w');
let bytes = 0;
for (const part of encodeFineGeography(mesh, fields, edgeSill, edgeSillSea)) { writeSync(fd, part); bytes += part.byteLength; }
closeSync(fd);
stamp(`${out.pathname ?? out}: ${(bytes / 2 ** 20).toFixed(0)} MiB`);
