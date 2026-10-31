import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createGeography, syntheticTopography } from '../js/geography.module.js';
import { oceanPacking } from '../js/gpu/oceanPacking.gpu.js';
import { meshSections, fillMeshTables } from '../js/gpu/meshTables.gpu.js';
import { oceanCells } from '../js/ocean/layered.module.js';

/*
 * The set of cells, edges and vertices the GPU ocean keeps rows for
 * (js/gpu/oceanPacking.gpu.js): it holds every element a sea kernel
 * reads a nonzero row of, its tables point everything else at a ghost, and
 * they are the mesh's tables where the set is the mesh.
 */
const mesh = buildMesh(new Grid(16), { radius: 6.371e6, omega: 7.29e-5 });
const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 300 : -4000));
const geography = createGeography(mesh, topography, {});
const { water, cellOcean } = oceanCells(mesh.nCells, geography), edgeOcean = water.edgeOcean;

test('the packed set holds the sea\'s cells, every edge of one, and the vertices of the sea\'s edges, each kind followed by one ghost', () => {
  const packing = oceanPacking(mesh, geography);
  assert.ok(packing.packed);
  assert.deepEqual(Array.from(packing.cells), Array.from({ length: mesh.nCells }, (_, i) => i).filter((i) => cellOcean[i]));
  const edges = new Set(packing.edges), vertices = new Set(packing.vertices);
  for (let e = 0; e < mesh.nEdges; e++) {
    const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1];
    assert.equal(edges.has(e), !!(cellOcean[a] || cellOcean[b]), `edge ${e}`);
    if (edgeOcean[e]) for (const v of [mesh.verticesOnEdge[2 * e], mesh.verticesOnEdge[2 * e + 1]]) {
      assert.ok(vertices.has(v), `vertex ${v} of the sea edge ${e}`);
      for (let m = 0; m < 3; m++) assert.ok(edges.has(mesh.edgesOnVertex[3 * v + m]), `edge ${mesh.edgesOnVertex[3 * v + m]} of the vertex ${v}`);
    }
  }
  assert.equal(packing.C, packing.cells.length + 1);
  assert.equal(packing.E, packing.edges.length + 1);
  assert.equal(packing.V, packing.vertices.length + 1);
  assert.ok(packing.C < mesh.nCells && packing.E < mesh.nEdges && packing.V < mesh.nVertices);
  const ghost = { cell: packing.C - 1, edge: packing.E - 1, vertex: packing.V - 1 }, sub = packing.mesh;
  assert.ok(cellOcean[packing.cells.at(-1)] && !cellOcean[packing.cellMap[ghost.cell]], 'the ghost cell reads a cell that is not the sea\'s');
  assert.ok(sub.areaCell[ghost.cell] > 0 && sub.dcEdge[ghost.edge] > 0 && sub.dvEdge[ghost.edge] > 0 && sub.areaTriangle[ghost.vertex] > 0, 'a ghost has a size');
  assert.equal(sub.fVertex[ghost.vertex], 0);
  for (let m = 0; m < 6; m++) assert.equal(sub.edgesOnCell[6 * ghost.cell + m], ghost.edge);
  for (let m = 0; m < 3; m++) assert.equal(sub.edgesOnVertex[3 * ghost.vertex + m], ghost.edge);
  assert.deepEqual(Array.from(sub.cellsOnEdge.subarray(2 * ghost.edge)), [ghost.cell, ghost.cell]);
  assert.deepEqual(Array.from(sub.verticesOnEdge.subarray(2 * ghost.edge)), [ghost.vertex, ghost.vertex]);
});

test('the packed set\'s tables hold the mesh\'s neighbors of each element and the ghost for the others', () => {
  const packing = oceanPacking(mesh, geography), { cells, edges, vertices, cellAt, edgeAt, vertexAt } = packing, ghost = { cell: packing.C - 1, edge: packing.E - 1, vertex: packing.V - 1 };
  const sections = meshSections(packing.C, packing.E, packing.V);
  const names = { MI: ['COE', 'VOE', 'EOC', 'ESC', 'COC', 'NEC', 'EOV', 'ESV', 'EOE'], MF: ['AREA', 'ATRI', 'DC', 'DV', 'FV', 'PVW', 'NEDGE'] };
  const seq = (list) => { const out = {}; let off = 0; for (const [name, n] of list) { out[name] = off; off += n; } out.total = off; return out; };
  const MI = seq(sections.MI.filter(([name]) => names.MI.includes(name))), MF = seq(sections.MF.filter(([name]) => names.MF.includes(name)));
  const mi = new Int32Array(MI.total), mf = new Float32Array(MF.total);
  fillMeshTables(packing.mesh, MI, MF, mi, mf);
  const full = meshSections(mesh.nCells, mesh.nEdges, mesh.nVertices), reference = { MI: seq(full.MI), MF: seq(full.MF) };
  const wholeInt = new Int32Array(reference.MI.total), wholeReal = new Float32Array(reference.MF.total);
  fillMeshTables(mesh, reference.MI, reference.MF, wholeInt, wholeReal);
  cells.forEach((g, s) => {
    assert.equal(mi[MI.NEC + s], wholeInt[reference.MI.NEC + g]);
    assert.equal(mf[MF.AREA + s], wholeReal[reference.MF.AREA + g]);
    for (let m = 0; m < 6; m++) {
      assert.equal(mi[MI.EOC + 6 * s + m], edgeAt[wholeInt[reference.MI.EOC + 6 * g + m]], `edge ${m} of cell ${g}`);
      assert.equal(mi[MI.ESC + 6 * s + m], wholeInt[reference.MI.ESC + 6 * g + m]);
      assert.equal(mi[MI.COC + 6 * s + m], cellAt[wholeInt[reference.MI.COC + 6 * g + m]]);
    }
  });
  edges.forEach((g, e) => {
    for (let side = 0; side < 2; side++) {
      assert.equal(mi[MI.COE + 2 * e + side], cellAt[wholeInt[reference.MI.COE + 2 * g + side]]);
      assert.equal(mi[MI.VOE + 2 * e + side], vertexAt[wholeInt[reference.MI.VOE + 2 * g + side]]);
    }
    assert.equal(mf[MF.DC + e], wholeReal[reference.MF.DC + g]);
    assert.equal(mf[MF.DV + e], wholeReal[reference.MF.DV + g]);
    for (let s = 0; s < 10; s++) assert.equal(mi[MI.EOE + 10 * e + s], s < mesh.nEdgesOnEdge[g] ? edgeAt[wholeInt[reference.MI.EOE + 10 * g + s]] : e);
    if (edgeOcean[g]) for (let s = 0; s < mesh.nEdgesOnEdge[g]; s++) assert.equal(mf[MF.PVW + 10 * e + s], wholeReal[reference.MF.PVW + 10 * g + s], `weight ${s} of the sea edge ${g}`);
  });
  vertices.forEach((g, v) => {
    for (let m = 0; m < 3; m++) {
      assert.equal(mi[MI.EOV + 3 * v + m], edgeAt[wholeInt[reference.MI.EOV + 3 * g + m]]);
      assert.notEqual(mi[MI.EOV + 3 * v + m], ghost.edge, 'a vertex of the set has its three edges in it');
    }
    assert.equal(mf[MF.FV + v], wholeReal[reference.MF.FV + g]);
  });
  assert.equal(vertexAt.filter((v) => v === ghost.vertex).length, mesh.nVertices - vertices.length);
  assert.equal(cellAt.filter((c) => c === ghost.cell).length, mesh.nCells - cells.length);
});

test('with the packing off, or where every element is kept, the set is the mesh and has no ghost', () => {
  const off = oceanPacking(mesh, geography, { packed: false });
  assert.ok(!off.packed);
  assert.equal(off.C, mesh.nCells); assert.equal(off.E, mesh.nEdges); assert.equal(off.V, mesh.nVertices);
  assert.deepEqual(Array.from(off.cellAt), Array.from({ length: mesh.nCells }, (_, i) => i));
  assert.deepEqual(Array.from(off.edgeMap), Array.from({ length: mesh.nEdges }, (_, e) => e));
  const aqua = oceanPacking(mesh, null);
  assert.ok(!aqua.packed);
  assert.equal(aqua.C, mesh.nCells); assert.equal(aqua.E, mesh.nEdges); assert.equal(aqua.V, mesh.nVertices);
});
