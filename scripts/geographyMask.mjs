#!/usr/bin/env node
/*
 * Writes data/geography_n1024_N<N>.bin (js/geography.module.js's mask
 * file) for each N from data/geography1024.bin: fineGeography.module.js's
 * maskGeography with the ice shelves as land, and with their water as sea
 * for the cavities. With REPORT, the per-N table as JSON there; with
 * LAKE_REPORT, scripts/fineGeography.mjs's lake list, for the lakes' names.
 *
 *   node --max-old-space-size=16000 scripts/geographyMask.mjs [N ...]
 */
import { writeFileSync, readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { gridMesh } from '../js/stateGrid.module.js';
import { readFineGeography, maskGeography, seaFine } from '../js/fineGeography.module.js';
import { encodeMask, decodeMask, maskUrl, createGeography } from '../js/geography.module.js';

const Ns = process.argv.length > 2 ? process.argv.slice(2).map(Number) : [16, 32, 64, 128, 192];
const deg = 180 / Math.PI, R = 6371.0088, clock = performance.now();
const names = process.env.LAKE_REPORT ? new Map(JSON.parse(readFileSync(process.env.LAKE_REPORT, 'utf8')).lakes.map((lake) => [lake.id, lake])) : null;
const fine = readFineGeography(), sea = seaFine(fine), report = {};
console.log(`fine geography read (${((performance.now() - clock) / 1000).toFixed(1)} s)`);

for (const N of Ns) {
  const start = performance.now(), mesh = gridMesh(new Grid(N)), C = mesh.nCells;
  const shelf = maskGeography(mesh, fine), water = maskGeography(mesh, sea);
  for (const g of [shelf, water]) g.saline = Uint8Array.from(g.lake, (l, i) => (l && (g.seawater[i] || names?.get(g.coarse.lakeId[i])?.featurecla === 'Alkaline Lake') ? 1 : 0));
  const parts = encodeMask(mesh, shelf, water), bytes = new Uint8Array(parts.reduce((n, p) => n + p.byteLength, 0));
  parts.reduce((at, p) => { bytes.set(p, at); return at + p.byteLength; }, 0);
  writeFileSync(maskUrl(N), bytes);
  const maskData = decodeMask(bytes), land = createGeography(mesh, null, { mask: 'n1024', maskData }), cavity = createGeography(mesh, null, { mask: 'n1024', maskData, iceShelves: 'cavity' });
  const area = (i) => mesh.areaCell[i] / mesh.radius ** 2 * R * R, where = (i) => [+(mesh.latCell[i] * deg).toFixed(2), +(mesh.lonCell[i] * deg).toFixed(2)];
  const describe = (cells) => ({ cells: cells.length, at: where(cells[0]), km2: Math.round(cells.reduce((s, i) => s + area(i), 0)), depth: Math.round(Math.max(...cells.map((i) => shelf.coarse.dryMaxDepth[i]))), elevation: Math.round(Math.min(...cells.map((i) => shelf.coarse.elevation[i]))) });
  const lakes = new Map();
  for (let i = 0; i < C; i++) {
    if (!shelf.lake[i] || shelf.seawater[i]) continue;
    const k = shelf.coarse.lakeId[i], entry = lakes.get(k) ?? { name: names?.get(k)?.name ?? `lake ${k}`, featurecla: names?.get(k)?.featurecla ?? null, cells: 0, km2: 0, level: 0, depth: 0 };
    entry.cells++; entry.km2 += area(i); entry.level += shelf.lakeLevel[i]; entry.depth += -shelf.elevation[i];
    lakes.set(k, entry);
  }
  let landCells = 0, landKm2 = 0, cavities = 0;
  for (let i = 0; i < C; i++) { if (shelf.land[i]) { landCells++; landKm2 += area(i); } if (cavity.cavity[i]) cavities++; }
  let mismatched = 0;
  for (let i = 0; i < C; i++) if (land.land[i] !== shelf.land[i] || Math.abs(land.elevation[i] - shelf.elevation[i]) > 0.5) mismatched++;
  report[N] = {
    cells: C, land: { cells: landCells, millionKm2: +(landKm2 / 1e6).toFixed(2) }, bytes: bytes.byteLength, seconds: +((performance.now() - start) / 1000).toFixed(1),
    lakes: [...lakes.values()].sort((a, b) => b.km2 - a.km2).map((l) => ({ ...l, km2: Math.round(l.km2), level: Math.round(l.level / l.cells), depth: Math.round(l.depth / l.cells) })),
    depressions: shelf.depressions.map(describe), inlets: shelf.inlets.map(describe), closed: shelf.closed.map(describe),
    unresolved: shelf.rule.unresolved.length, conflicts: shelf.rule.conflicts.length, cavities, seaModeCells: maskData.sea.index.length, seaModeEdges: maskData.sea.edges.length, roundTripMismatches: mismatched,
  };
  const { lakes: l, depressions: d, inlets: s, closed: c, ...brief } = report[N];
  console.log(`N=${N}`, JSON.stringify(brief), `lakes ${l.length}, depressions ${JSON.stringify(d.map((x) => [x.at, x.cells, x.depth]))}, seawater ${JSON.stringify(s.map((x) => [x.at, x.cells, x.depth]))}, closed ${JSON.stringify(c.map((x) => [x.at, x.cells, x.depth]))}`);
}
if (process.env.REPORT) writeFileSync(process.env.REPORT, JSON.stringify(report, null, 1));
