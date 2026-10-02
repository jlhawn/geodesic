#!/usr/bin/env node
/*
 * The mesh around a few cells for scripts/subgridTerrainHand.py, as JSON:
 * per cell its centre, the centres of its two rings of neighbours, every
 * triangle of cell centres around them with the resolved orography at its
 * corners (the surface geopotential over g), and the cell's values in
 * data/subgrid_N<N>.bin.
 *
 *   node scripts/subgridTerrainCells.mjs N lat,lon [lat,lon ...] > cells.json
 */
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { topographyFromInt16, meshSubgrid } from '../js/geography.module.js';
import { createModel } from '../js/model.module.js';

const N = Number(process.argv[2]);
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const model = createModel(new Grid(N), { physics: false, topography });
const { mesh, core } = model, g = core.diagnostics.g;
const { xCell, latCell, lonCell, cellsOnCell, nEdgesOnCell, maxEdges, verticesOnCell, cellsOnVertex, nCells } = mesh;
const fields = meshSubgrid(mesh);
const height = (i) => model.surfaceGeopotential[i] / g;
const point = (i) => [xCell[3 * i], xCell[3 * i + 1], xCell[3 * i + 2]];
const out = [];
for (const spec of process.argv.slice(3)) {
  const [lat, lon] = spec.split(',').map((v) => Number(v) * Math.PI / 180);
  const x = [Math.cos(lat) * Math.cos(lon), Math.cos(lat) * Math.sin(lon), Math.sin(lat)];
  let cell = 0, best = -2;
  for (let i = 0; i < nCells; i++) { const d = x[0] * xCell[3 * i] + x[1] * xCell[3 * i + 1] + x[2] * xCell[3 * i + 2]; if (d > best) { best = d; cell = i; } }
  const ring = new Set([cell]);
  for (let pass = 0; pass < 2; pass++) for (const i of [...ring]) for (let m = 0; m < nEdgesOnCell[i]; m++) ring.add(cellsOnCell[maxEdges * i + m]);
  const triangles = new Map();
  for (const i of ring) for (let m = 0; m < nEdgesOnCell[i]; m++) {
    const v = verticesOnCell[maxEdges * i + m];
    if (v >= 0) triangles.set(v, [0, 1, 2].map((j) => cellsOnVertex[3 * v + j]));
  }
  out.push({
    cell, lat: latCell[cell] * 180 / Math.PI, lon: lonCell[cell] * 180 / Math.PI, land: model.geography.land[cell],
    centre: point(cell), ring: [...ring].map((i) => ({ i, x: point(i) })),
    triangles: [...triangles.values()].map((t) => t.map((i) => ({ x: point(i), h: height(i) }))),
    file: { deviation: fields.deviation[cell], anisotropy: fields.anisotropy[cell], orientation: fields.orientation[cell], slope: fields.slope[cell], filtered: fields.filtered[cell] },
  });
}
console.log(JSON.stringify(out));
