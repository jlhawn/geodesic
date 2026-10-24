import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { gridMesh } from '../js/stateGrid.module.js';
import { relaxedMesh, fineEdges, encodeFineGeography, decodeFineGeography, coarseGeography, connectedGeography, FIELDS } from '../js/fineGeography.module.js';

test('relaxedMesh reproduces Grid bit for bit at N=16', () => {
  const grid = new Grid(16, { relax: 30 }), mesh = relaxedMesh(16, 30);
  for (const cell of grid) {
    const i = cell.index, c = cell.centerVertex;
    assert.deepEqual([mesh.xCell[3 * i], mesh.xCell[3 * i + 1], mesh.xCell[3 * i + 2], mesh.area[i]], [c.x, c.y, c.z, cell.area]);
    assert.deepEqual(Array.from(mesh.cellsOnCell.subarray(6 * i, 6 * i + cell.neighbors.length)), cell.neighbors.map((n) => n.index));
  }
});

test('a fine geography survives its file, and a coarse mesh takes its land by area and its sills from the fine edges', () => {
  const fine = relaxedMesh(64, 30), C = fine.nCells, { nEdges, cellsOnEdge } = fineEdges(fine), z = (f) => fine.xCell[3 * f + 2];
  const fields = Object.fromEntries(Object.keys(FIELDS).map((name) => [name, new Float64Array(C)]));
  for (let f = 0; f < C; f++) {
    const land = z(f) > 0.3 ? 1 : 0;
    fields.points[f] = 14; fields.landFraction[f] = land; fields.elevation[f] = land ? 500 : -4000;
    fields.landElevation[f] = land ? 500 : 0; fields.seaElevation[f] = land ? 0 : -4000; fields.bed[f] = fields.elevation[f]; fields.maxDepth[f] = land ? 0 : 4000;
  }
  const edgeSill = Float64Array.from({ length: nEdges }, (_, e) => (fields.maxDepth[cellsOnEdge[2 * e]] > 0 && fields.maxDepth[cellsOnEdge[2 * e + 1]] > 0 ? 3000 : 0));
  const parts = encodeFineGeography(fine, fields, edgeSill), bytes = new Uint8Array(parts.reduce((n, p) => n + p.byteLength, 0));
  parts.reduce((at, p) => { bytes.set(p, at); return at + p.byteLength; }, 0);
  const read = decodeFineGeography(bytes);
  assert.equal(read.nCells, C);
  assert.equal(read.nEdges, nEdges);
  assert.deepEqual(Array.from(read.landFraction), Array.from(fields.landFraction));
  const mesh = gridMesh(new Grid(16)), coarse = coarseGeography(mesh, read);
  let fineLand = 0, coarseLand = 0, total = 0;
  for (let f = 0; f < C; f++) { fineLand += read.area[f] * read.landFraction[f]; total += read.area[f]; }
  for (let i = 0; i < mesh.nCells; i++) coarseLand += coarse.area[i] * coarse.landFraction[i];
  assert.ok(Math.abs(coarseLand - fineLand) < 1e-9 * total);
  for (let i = 0; i < mesh.nCells; i++) {
    const zi = mesh.xCell[3 * i + 2];
    if (zi > 0.45) assert.equal(coarse.land[i], 1);
    if (zi < 0.15) assert.equal(coarse.land[i], 0);
  }
  for (let e = 0; e < mesh.nEdges; e++) {
    const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1];
    if (coarse.landFraction[a] === 0 && coarse.landFraction[b] === 0) assert.equal(coarse.edgeSill[e], 3000);
    if (coarse.landFraction[a] === 1 || coarse.landFraction[b] === 1) assert.equal(coarse.edgeSill[e], 0);
  }
  const rule = connectedGeography(mesh, coarse, read);
  assert.equal(rule.flipped.filter(Boolean).length, 0);
});
