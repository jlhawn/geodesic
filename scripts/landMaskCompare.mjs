#!/usr/bin/env node
/*
 * Compares the geography of each model mesh built from data/geography1024.bin
 * (js/fineGeography.module.js: the area rule alone, and with the
 * connectivity rule) with createGeography's from the 0.25° ETOPO 2022
 * raster: the cells that change class, the Panama, Cuba and San Francisco
 * boxes, and the ocean's widest path through named straits on each mesh
 * and on the fine mesh. Writes the details as JSON to OUT.
 *
 *   node --max-old-space-size=12000 scripts/landMaskCompare.mjs OUT [N ...]
 */
import { writeFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { gridMesh } from '../js/stateGrid.module.js';
import { createGeography, readTopography, geographyOptions, FRESH_GEOGRAPHY } from '../js/geography.module.js';
import { readFineGeography, coarseGeography, connectedGeography } from '../js/fineGeography.module.js';

const STRAITS = {
  gibraltar: { a: [36.0, -9.0], b: [37.0, -1.0], box: [34.5, 38.5, -10, 0] },
  babElMandeb: { a: [16.5, 40.5], b: [12.0, 45.0], box: [11, 17, 40, 46] },
  hormuz: { a: [26.6, 53.5], b: [25.0, 57.8], box: [23.8, 27.8, 52.5, 59] },
  bering: { a: [62.5, -172], b: [67.5, -168.5], box: [62, 69, -177, -163] },
  lombok: { a: [-8.2, 115.85], b: [-9.0, 115.8], box: [-9.6, -7.8, 115.5, 116.1] },
  panama: { a: [8.0, -79.5], b: [9.7, -79.5], box: [6.5, 11, -84, -76] },
  danish: { a: [57.3, 11.5], b: [55.3, 15.5], box: [54.2, 58.0, 9.5, 16.5] },
  bosporus: { a: [42.2, 29.5], b: [40.8, 28.3], box: [40.4, 42.6, 27.5, 30.5] },
  hudson: { a: [60.0, -84.0], b: [61.0, -62.5], box: [57.5, 65.5, -86, -61] },
  tsushima: { a: [37.0, 131.0], b: [32.5, 127.5], box: [31.5, 37.5, 126.5, 131.3] },
  tatar: { a: [47.0, 139.8], b: [54.2, 141.2], box: [46.5, 55, 137, 142.3] },
  florida: { a: [25.0, -85.0], b: [27.5, -79.8], box: [22.6, 28.5, -86, -79.0] },
  denmarkStrait: { a: [64.5, -31.0], b: [68.0, -22.0], box: [63.5, 69.0, -36, -20.5] },
  faroeBank: { a: [61.8, -3.0], b: [61.0, -12.0], box: [59.5, 62.4, -13, -2.5] },
};
const BOXES = { darien: [7, 10, -80, -76.5], cuba: [19.8, 23.3, -85, -74.1], sanFrancisco: [37.3, 38.3, -122.8, -121.8] };
const deg = 180 / Math.PI, R = 6371.0088;
const out = process.argv[2], Ns = process.argv.slice(3).map(Number);
const fine = readFineGeography(), F = fine.nCells;
const fineLat = Float64Array.from({ length: F }, (_, f) => Math.asin(Math.max(-1, Math.min(1, fine.xCell[3 * f + 2]))) * deg);
const fineLon = Float64Array.from({ length: F }, (_, f) => Math.atan2(fine.xCell[3 * f + 1], fine.xCell[3 * f]) * deg);
const topography = readTopography('etopo2022'), { geography: options } = geographyOptions(FRESH_GEOGRAPHY);

const inBox = (lat, lon, [s, n, w, e], margin = 0) => lat >= s - margin && lat <= n + margin && lon >= w - margin / Math.cos((s + n) / 2 / deg) && lon <= e + margin / Math.cos((s + n) / 2 / deg);
let last = null;
function widest({ count, lat, lon, valid, links, deep }, { a, b, box }, margin, ends = null) {
  const nodes = [];
  for (let i = 0; i < count; i++) if (valid(i) && inBox(lat[i], lon[i], box, margin)) nodes.push(i);
  if (!deep) {
    for (const [la, lo] of [a, b]) {
      let best = -1, d = Infinity;
      for (let i = 0; i < count; i++) { if (!valid(i)) continue; const x = (lat[i] - la) ** 2 + ((lon[i] - lo) * Math.cos(la / deg)) ** 2; if (x < d) { d = x; best = i; } }
      if (best >= 0 && !nodes.includes(best)) nodes.push(best);
    }
  }
  const near = ([la, lo]) => {
    let best = -1, d = Infinity, deepest = -1, most = -1;
    for (const i of nodes) {
      const x = (lat[i] - la) ** 2 + ((lon[i] - lo) * Math.cos(la / deg)) ** 2;
      if (x < d) { d = x; best = i; }
      if (deep && x < 0.3 ** 2 && deep(i) > most) { most = deep(i); deepest = i; }
    }
    return deepest >= 0 ? deepest : best;
  };
  const [s, t] = ends ?? [near(a), near(b)];
  if (s < 0 || t < 0) return null;
  if (ends) for (const i of ends) if (!nodes.includes(i)) nodes.push(i);
  last = [s, t];
  const member = new Set(nodes), width = new Map([[s, Infinity]]), heap = [[Infinity, s]];
  const push = (x) => { heap.push(x); let k = heap.length - 1; while (k > 0) { const p = (k - 1) >> 1; if (heap[p][0] >= heap[k][0]) break; [heap[p], heap[k]] = [heap[k], heap[p]]; k = p; } };
  const pop = () => { const top = heap[0], last = heap.pop(); if (heap.length) { heap[0] = last; let k = 0; for (;;) { const l = 2 * k + 1, r = l + 1; let m = k; if (l < heap.length && heap[l][0] > heap[m][0]) m = l; if (r < heap.length && heap[r][0] > heap[m][0]) m = r; if (m === k) break; [heap[m], heap[k]] = [heap[k], heap[m]]; k = m; } } return top; };
  const done = new Set();
  while (heap.length) {
    const [w, i] = pop();
    if (done.has(i)) continue;
    done.add(i);
    if (i === t) return Number.isFinite(w) ? Math.round(w) : 'one cell';
    for (const [j, s2] of links(i)) {
      if (!member.has(j) || !(s2 > 0)) continue;
      const v = Math.min(w, s2);
      if (v > (width.get(j) ?? 0)) { width.set(j, v); push([v, j]); }
    }
  }
  return 0;
}
const fineGraph = (resolved) => ({
  count: F, lat: fineLat, lon: fineLon,
  valid: (f) => fine.maxDepth[f] > 0 && (!resolved || fine.landFraction[f] <= 0.5), deep: (f) => fine.maxDepth[f],
  links: function* (f) { for (let k = 0; k < fine.nEdgesOnCell[f]; k++) yield [fine.cellsOnCell[6 * f + k], fine.edgeSill[fine.edgeOfSlot[6 * f + k]]]; },
});
const coarseGraph = (mesh, lat, lon, land, open, sill) => ({
  count: mesh.nCells, lat, lon, valid: (i) => !land[i],
  links: function* (i) { for (let m = 0; m < mesh.nEdgesOnCell[i]; m++) { const e = mesh.edgesOnCell[mesh.maxEdges * i + m]; if (open[e]) yield [mesh.cellsOnCell[mesh.maxEdges * i + m], sill[e]]; } },
});

const report = { fine: {}, meshes: {} };
for (const [name, strait] of Object.entries(STRAITS)) report.fine[name] = { sill: widest(fineGraph(false), strait, 0), resolved: widest(fineGraph(true), strait, 0) };
let fineLand = 0, fineArea = 0;
for (let f = 0; f < F; f++) { fineArea += fine.area[f]; fineLand += fine.area[f] * fine.landFraction[f]; }
report.fine.landMillionKm2 = +(fineLand / fineArea * 4 * Math.PI * R * R / 1e6).toFixed(2);
console.log('fine N=1024 widest sills (sill graph / resolved):', JSON.stringify(report.fine));

for (const N of Ns) {
  const clock = performance.now(), mesh = gridMesh(new Grid(N, { relax: 30 })), C = mesh.nCells;
  const lat = Float64Array.from(mesh.latCell, (v) => v * deg), lon = Float64Array.from(mesh.lonCell, (v) => v * deg);
  const today = createGeography(mesh, topography, options), raw = coarseGeography(mesh, fine), rule = connectedGeography(mesh, raw, fine);
  const rawOpen = Uint8Array.from({ length: mesh.nEdges }, (_, e) => (!raw.land[mesh.cellsOnEdge[2 * e]] && !raw.land[mesh.cellsOnEdge[2 * e + 1]] && raw.edgeSill[e] > 0 ? 1 : 0));
  const area = (i) => mesh.areaCell[i] / mesh.radius ** 2 * R * R;
  const r = { cells: C, land: {}, flips: {}, positiveSea: {}, boxes: {}, straits: {}, ruleFlips: [], conflicts: rule.conflicts.map((i) => [+lat[i].toFixed(2), +lon[i].toFixed(2)]) };
  for (const [key, g] of Object.entries({ today, raw, rule })) {
    let n = 0, a = 0, pos = 0;
    for (let i = 0; i < C; i++) { if (g.land[i]) { n++; a += area(i); } else if (g.elevation[i] > 0) pos++; }
    r.land[key] = { cells: n, millionKm2: +(a / 1e6).toFixed(2) };
    r.positiveSea[key] = pos;
  }
  for (const [key, g] of Object.entries({ raw, rule })) {
    let toSea = 0, toLand = 0, km2 = 0;
    for (let i = 0; i < C; i++) if (!today.land[i] !== !g.land[i]) { if (today.land[i]) toSea++; else toLand++; km2 += area(i); }
    r.flips[key] = { toSea, toLand, km2: Math.round(km2) };
  }
  let d2 = 0, dl = 0, nl = 0;
  for (let i = 0; i < C; i++) { d2 += (raw.elevation[i] - today.elevation[i]) ** 2; if (today.land[i] && raw.land[i]) { dl += (raw.elevation[i] - today.elevation[i]) ** 2; nl++; } }
  r.elevationRms = { all: +Math.sqrt(d2 / C).toFixed(1), landBoth: +Math.sqrt(dl / nl).toFixed(1) };
  for (const [name, box] of Object.entries(BOXES)) {
    const cells = [];
    for (let i = 0; i < C; i++) if (inBox(lat[i], lon[i], box)) cells.push(i);
    r.boxes[name] = {
      cells: cells.length, landToday: cells.filter((i) => today.land[i]).length, landRaw: cells.filter((i) => raw.land[i]).length, landRule: cells.filter((i) => rule.land[i]).length,
      list: cells.length <= 400 ? cells.map((i) => ({ lat: +lat[i].toFixed(2), lon: +lon[i].toFixed(2), todayLand: today.land[i], todayFraction: +today.landFraction[i].toFixed(2), todayElevation: Math.round(today.elevation[i]), fraction: +raw.landFraction[i].toFixed(3), elevation: Math.round(raw.elevation[i]), rule: rule.land[i], flipped: rule.flipped[i] })) : undefined,
    };
  }
  const margin = 0.5 * Math.sqrt(4 * Math.PI / (10 * N * N)) * deg;
  for (const [name, strait] of Object.entries(STRAITS)) {
    r.straits[name] = {
      today: widest(coarseGraph(mesh, lat, lon, today.land, today.edgeOcean, today.edgeSill), strait, margin),
      raw: widest(coarseGraph(mesh, lat, lon, raw.land, rawOpen, raw.edgeSill), strait, margin),
      rule: widest(coarseGraph(mesh, lat, lon, rule.land, rule.edgeOcean, rule.edgeSill), strait, margin),
    };
    const [s, t] = last ?? [-1, -1];
    r.straits[name].fineBetweenCells = s >= 0 && raw.rep[s] >= 0 && raw.rep[t] >= 0 ? widest(fineGraph(false), strait, margin, [raw.rep[s], raw.rep[t]]) : null;
  }
  for (let i = 0; i < C; i++) if (rule.flipped[i]) r.ruleFlips.push({ lat: +lat[i].toFixed(2), lon: +lon[i].toFixed(2), to: rule.flipped[i] < 0 ? 'sea' : 'land', fraction: +raw.landFraction[i].toFixed(2), elevation: Math.round(rule.elevation[i]), todayLand: today.land[i] });
  r.components = { sea: rule.seaComponents, land: rule.landComponents };
  const describe = (cells) => ({ cells: cells.length, lat: +lat[cells[0]].toFixed(2), lon: +lon[cells[0]].toFixed(2), seaToday: cells.filter((i) => !today.land[i]).length, maxDepth: Math.round(Math.max(...cells.map((i) => raw.maxDepth[i]))) });
  r.unresolved = rule.unresolved.map(describe);
  r.isolated = rule.isolated.map(describe);
  report.meshes[N] = r;
  const { ruleFlips, boxes, ...brief } = r;
  console.log(`N=${N} (${((performance.now() - clock) / 1000).toFixed(1)} s)`, JSON.stringify(brief), 'rule flips', ruleFlips.length, 'boxes', JSON.stringify(Object.fromEntries(Object.entries(boxes).map(([k, v]) => [k, [v.cells, v.landToday, v.landRaw, v.landRule]]))));
}
writeFileSync(out, JSON.stringify(report, null, 1));
