import { oceanCells } from '../ocean/layered.module.js';
import { MAX_EDGES } from './meshTables.gpu.js';

/*
 * The cells, edges and vertices the GPU ocean keeps rows for. A layered
 * row costs a float for every layer, so the ocean's buffers hold the cells
 * of the sea, every edge of one of them, and the vertices of the edges the
 * sea flows through, in the mesh's order, each of the three followed by a
 * ghost. Everything the ocean's kernels read beyond that set is a row the
 * kernels hold at zero, since the velocity of an edge that is not the
 * sea's is zero and so are the divergence and the curl of the velocity at
 * a cell or vertex without such an edge; the tables point every such
 * element at its kind's ghost, which has a unit size and no rotation.
 * The ∇⁴ closure's second pass leaves a nonzero divergence on the land
 * beside a coast, and no sea edge reads it. When every element is kept
 * there are no ghosts.
 *
 * cells, edges and vertices list the mesh index of each packed element,
 * before the ghosts. cellAt, edgeAt and vertexAt give the packed index of
 * each mesh element, the ghost's where it is left out. cellMap and edgeMap
 * give for every packed cell and edge, the ghosts included, the mesh index
 * the atmosphere's PH, S and D are read at; a ghost reads a cell or edge
 * outside the set. mesh is the set as a mesh for the table builder.
 */
export function oceanPacking(mesh, geography, { packed = true } = {}) {
  const C = mesh.nCells, E = mesh.nEdges, V = mesh.nVertices, M = mesh.maxEdges, ME = mesh.maxEdgesOnEdge;
  const { water, cellOcean } = oceanCells(C, geography || null);
  const edgeOcean = water ? water.edgeOcean : new Uint8Array(E).fill(1);
  const keepCell = Uint8Array.from(cellOcean, (wet) => (packed ? wet : 1)), keepEdge = new Uint8Array(E), keepVertex = new Uint8Array(V).fill(packed ? 0 : 1);
  for (let e = 0; e < E; e++) {
    keepEdge[e] = !packed || keepCell[mesh.cellsOnEdge[2 * e]] || keepCell[mesh.cellsOnEdge[2 * e + 1]] ? 1 : 0;
    if (edgeOcean[e]) keepVertex[mesh.verticesOnEdge[2 * e]] = keepVertex[mesh.verticesOnEdge[2 * e + 1]] = 1;
  }
  const kept = (flags) => Int32Array.from(flags.flatMap((keep, n) => (keep ? [n] : [])));
  const cells = kept(Array.from(keepCell)), edges = kept(Array.from(keepEdge)), vertices = kept(Array.from(keepVertex));
  const whole = cells.length === C && edges.length === E && vertices.length === V, ghosts = whole ? 0 : 1;
  const nC = cells.length + ghosts, nE = edges.length + ghosts, nV = vertices.length + ghosts;
  const ghostCell = cells.length, ghostEdge = edges.length, ghostVertex = vertices.length;
  const at = (list, count, ghost) => { const out = new Int32Array(count).fill(ghost); list.forEach((g, s) => { out[g] = s; }); return out; };
  const cellAt = at(cells, C, ghostCell), edgeAt = at(edges, E, ghostEdge), vertexAt = at(vertices, V, ghostVertex);

  const sub = {
    nCells: nC, nEdges: nE, nVertices: nV, maxEdges: M, maxEdgesOnEdge: ME,
    nEdgesOnCell: new Int32Array(nC), edgesOnCell: new Int32Array(M * nC), edgeSignOnCell: new Int32Array(M * nC), cellsOnCell: new Int32Array(M * nC), areaCell: new Float64Array(nC),
    cellsOnEdge: new Int32Array(2 * nE), verticesOnEdge: new Int32Array(2 * nE), dcEdge: new Float64Array(nE), dvEdge: new Float64Array(nE), nEdge: new Float64Array(3 * nE),
    nEdgesOnEdge: new Int32Array(nE), edgesOnEdge: new Int32Array(ME * nE), weightsOnEdge: new Float64Array(ME * nE),
    cellsOnVertex: new Int32Array(3 * nV), edgesOnVertex: new Int32Array(3 * nV), edgeSignOnVertex: new Int32Array(3 * nV), areaTriangle: new Float64Array(nV), fVertex: new Float64Array(nV),
  };
  cells.forEach((g, s) => {
    sub.nEdgesOnCell[s] = mesh.nEdgesOnCell[g]; sub.areaCell[s] = mesh.areaCell[g];
    for (let m = 0; m < mesh.nEdgesOnCell[g]; m++) {
      sub.edgesOnCell[M * s + m] = edgeAt[mesh.edgesOnCell[M * g + m]]; sub.edgeSignOnCell[M * s + m] = mesh.edgeSignOnCell[M * g + m]; sub.cellsOnCell[M * s + m] = cellAt[mesh.cellsOnCell[M * g + m]];
    }
  });
  edges.forEach((g, s) => {
    for (let side = 0; side < 2; side++) { sub.cellsOnEdge[2 * s + side] = cellAt[mesh.cellsOnEdge[2 * g + side]]; sub.verticesOnEdge[2 * s + side] = vertexAt[mesh.verticesOnEdge[2 * g + side]]; }
    sub.dcEdge[s] = mesh.dcEdge[g]; sub.dvEdge[s] = mesh.dvEdge[g];
    for (let c = 0; c < 3; c++) sub.nEdge[3 * s + c] = mesh.nEdge[3 * g + c];
    sub.nEdgesOnEdge[s] = mesh.nEdgesOnEdge[g];
    for (let m = 0; m < mesh.nEdgesOnEdge[g]; m++) { sub.edgesOnEdge[ME * s + m] = edgeAt[mesh.edgesOnEdge[ME * g + m]]; sub.weightsOnEdge[ME * s + m] = mesh.weightsOnEdge[ME * g + m]; }
  });
  vertices.forEach((g, s) => {
    for (let m = 0; m < 3; m++) { sub.cellsOnVertex[3 * s + m] = cellAt[mesh.cellsOnVertex[3 * g + m]]; sub.edgesOnVertex[3 * s + m] = edgeAt[mesh.edgesOnVertex[3 * g + m]]; sub.edgeSignOnVertex[3 * s + m] = mesh.edgeSignOnVertex[3 * g + m]; }
    sub.areaTriangle[s] = mesh.areaTriangle[g]; sub.fVertex[s] = mesh.fVertex[g];
  });
  const cellMap = new Int32Array(nC), edgeMap = new Int32Array(nE);
  cellMap.set(cells); edgeMap.set(edges);
  if (ghosts) {
    const outside = (flags) => Math.max(0, flags.indexOf(0));
    cellMap[ghostCell] = outside(Array.from(keepCell)); edgeMap[ghostEdge] = outside(Array.from(keepEdge));
    sub.nEdgesOnCell[ghostCell] = MAX_EDGES; sub.areaCell[ghostCell] = 1;
    for (let m = 0; m < MAX_EDGES; m++) { sub.edgesOnCell[M * ghostCell + m] = ghostEdge; sub.edgeSignOnCell[M * ghostCell + m] = 1; sub.cellsOnCell[M * ghostCell + m] = ghostCell; }
    for (let side = 0; side < 2; side++) { sub.cellsOnEdge[2 * ghostEdge + side] = ghostCell; sub.verticesOnEdge[2 * ghostEdge + side] = ghostVertex; }
    sub.dcEdge[ghostEdge] = 1; sub.dvEdge[ghostEdge] = 1; sub.nEdge[3 * ghostEdge] = 1;
    for (let m = 0; m < 3; m++) { sub.cellsOnVertex[3 * ghostVertex + m] = ghostCell; sub.edgesOnVertex[3 * ghostVertex + m] = ghostEdge; sub.edgeSignOnVertex[3 * ghostVertex + m] = 1; }
    sub.areaTriangle[ghostVertex] = 1;
  }
  return { packed: !whole, cells, edges, vertices, C: nC, E: nE, V: nV, cellAt, edgeAt, vertexAt, cellMap, edgeMap, mesh: sub, cellCount: C, edgeCount: E, vertexCount: V };
}
