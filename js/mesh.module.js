import * as THREE from "./three.module.js";

export const EARTH_RADIUS = 6371e3;
export const SIDEREAL_OMEGA = 7.2921159e-5;
export const MAX_EDGES = 6;
export const MAX_EDGES_ON_EDGE = 2 * (MAX_EDGES - 1);

const _p = new THREE.Vector3();
const _q = new THREE.Vector3();
const _r = new THREE.Vector3();

function arcLength(a, b) {
  return Math.atan2(_p.crossVectors(a, b).length(), a.dot(b));
}

function sphericalExcess(a, b, c) {
  return 2 * Math.atan2(_p.crossVectors(b, c).dot(a), 1 + a.dot(b) + b.dot(c) + c.dot(a));
}

function latitude(p) {
  return Math.atan2(p.z, Math.hypot(p.x, p.y));
}

function longitude(p) {
  return Math.atan2(p.y, p.x);
}

function setVector(array, index, v) {
  array[3 * index] = v.x;
  array[3 * index + 1] = v.y;
  array[3 * index + 2] = v.z;
}

function getVector(array, index, out) {
  return out.set(array[3 * index], array[3 * index + 1], array[3 * index + 2]);
}

/*
 * Structure-of-arrays mesh for the C-grid, following the MPAS mesh
 * specification's names. Positions are unit vectors; lengths and areas are
 * in metres from `radius`. Cell and triangle areas are sums of their
 * spherical kites, so the kite partition identities hold to roundoff.
 * Orientation conventions:
 *  - nEdge points from cellsOnEdge[2e] to cellsOnEdge[2e+1]; tEdge = xEdge × nEdge.
 *  - verticesOnEdge[2e+1] is the vertex on the tEdge side.
 *  - edgeSignOnCell is +1 where nEdge points out of the cell.
 *  - edgeSignOnVertex is +1 where travelling along nEdge circulates
 *    counter-clockwise (viewed from outside) around the vertex.
 *  - All *OnCell and *OnVertex rings run counter-clockwise viewed from
 *    outside; verticesOnCell[6i+k] lies between edgesOnCell[6i+k] and
 *    [6i+k+1]; edgesOnVertex[3v+m] joins cellsOnVertex[3v+m] and [3v+m+1].
 *  - weightsOnEdge[10e+s] applies to edgesOnEdge[10e+s]:
 *      w(e,e') = n(e,i) n(e',i) (1/2 − Σ R(i,v)/A(i)) over the vertices v of
 *      cell i passed walking counter-clockwise from e to e'.
 */
export function buildMesh(grid, { radius = EARTH_RADIUS, omega = SIDEREAL_OMEGA } = {}) {
  const cells = [...grid];
  const nCells = cells.length;
  const nVertices = grid.vertices.length;
  const nEdges = 3 * nCells - 6;
  const R2 = radius * radius;
  const coriolis = (lat) => 2 * omega * Math.sin(lat);

  const xCell = new Float64Array(3 * nCells);
  const latCell = new Float64Array(nCells);
  const lonCell = new Float64Array(nCells);
  const areaCell = new Float64Array(nCells);
  const fCell = new Float64Array(nCells);
  const nEdgesOnCell = new Int32Array(nCells);
  const cellsOnCell = new Int32Array(MAX_EDGES * nCells).fill(-1);
  const verticesOnCell = new Int32Array(MAX_EDGES * nCells).fill(-1);
  const edgesOnCell = new Int32Array(MAX_EDGES * nCells).fill(-1);
  const edgeSignOnCell = new Int32Array(MAX_EDGES * nCells);
  const kiteAreasOnCell = new Float64Array(MAX_EDGES * nCells);

  const xVertex = new Float64Array(3 * nVertices);
  const latVertex = new Float64Array(nVertices);
  const lonVertex = new Float64Array(nVertices);
  const areaTriangle = new Float64Array(nVertices);
  const fVertex = new Float64Array(nVertices);
  const cellsOnVertex = new Int32Array(3 * nVertices).fill(-1);
  const edgesOnVertex = new Int32Array(3 * nVertices).fill(-1);
  const edgeSignOnVertex = new Int32Array(3 * nVertices);
  const kiteAreasOnVertex = new Float64Array(3 * nVertices);

  const cellsOnEdge = new Int32Array(2 * nEdges);
  const verticesOnEdge = new Int32Array(2 * nEdges);
  const dcEdge = new Float64Array(nEdges);
  const dvEdge = new Float64Array(nEdges);
  const xEdge = new Float64Array(3 * nEdges);
  const nEdge = new Float64Array(3 * nEdges);
  const tEdge = new Float64Array(3 * nEdges);
  const latEdge = new Float64Array(nEdges);
  const fEdge = new Float64Array(nEdges);
  const nEdgesOnEdge = new Int32Array(nEdges);
  const edgesOnEdge = new Int32Array(MAX_EDGES_ON_EDGE * nEdges).fill(-1);
  const weightsOnEdge = new Float64Array(MAX_EDGES_ON_EDGE * nEdges);

  for (const cell of cells) {
    const i = cell.index;
    const n = cell.neighbors.length;
    setVector(xCell, i, cell.centerVertex);
    latCell[i] = latitude(cell.centerVertex);
    lonCell[i] = longitude(cell.centerVertex);
    fCell[i] = coriolis(latCell[i]);
    nEdgesOnCell[i] = n;
    for (let k = 0; k < n; k++) {
      const vertex = cell.vertices[k];
      const v = vertex.index;
      cellsOnCell[MAX_EDGES * i + k] = cell.neighbors[k].index;
      verticesOnCell[MAX_EDGES * i + k] = v;
      if (cellsOnVertex[3 * v] !== -1) continue;
      const a = cell.neighbors[k];
      const b = cell.neighbors[(k + 1) % n];
      cellsOnVertex[3 * v] = i;
      cellsOnVertex[3 * v + 1] = a.index;
      cellsOnVertex[3 * v + 2] = b.index;
      setVector(xVertex, v, vertex);
      latVertex[v] = latitude(vertex);
      lonVertex[v] = longitude(vertex);
      fVertex[v] = coriolis(latVertex[v]);
    }
  }

  let edgeCount = 0;
  for (const cell of cells) {
    const i = cell.index;
    const n = cell.neighbors.length;
    for (let k = 0; k < n; k++) {
      const other = cell.neighbors[k];
      const j = other.index;
      if (j < i) {
        edgesOnCell[MAX_EDGES * i + k] = edgesOnCell[MAX_EDGES * j + other.neighbors.indexOf(cell)];
        edgeSignOnCell[MAX_EDGES * i + k] = -1;
        continue;
      }
      const e = edgeCount++;
      const v1 = cell.vertices[(k - 1 + n) % n];
      const v2 = cell.vertices[k];
      cellsOnEdge[2 * e] = i;
      cellsOnEdge[2 * e + 1] = j;
      verticesOnEdge[2 * e] = v1.index;
      verticesOnEdge[2 * e + 1] = v2.index;
      dcEdge[e] = arcLength(cell.centerVertex, other.centerVertex) * radius;
      dvEdge[e] = arcLength(v1, v2) * radius;
      const midpoint = _q.addVectors(cell.centerVertex, other.centerVertex).normalize();
      const normal = _r.subVectors(other.centerVertex, cell.centerVertex).normalize();
      setVector(xEdge, e, midpoint);
      setVector(nEdge, e, normal);
      setVector(tEdge, e, _p.crossVectors(midpoint, normal));
      latEdge[e] = latitude(midpoint);
      fEdge[e] = coriolis(latEdge[e]);
      edgesOnCell[MAX_EDGES * i + k] = e;
      edgeSignOnCell[MAX_EDGES * i + k] = 1;
    }
  }
  if (edgeCount !== nEdges) {
    throw new Error(`built ${edgeCount} edges, expected ${nEdges}`);
  }

  for (let v = 0; v < nVertices; v++) {
    for (let m = 0; m < 3; m++) {
      const c = cellsOnVertex[3 * v + m];
      const next = cellsOnVertex[3 * v + (m + 1) % 3];
      let e = -1;
      for (let k = 0; k < nEdgesOnCell[c]; k++) {
        if (cellsOnCell[MAX_EDGES * c + k] === next) e = edgesOnCell[MAX_EDGES * c + k];
      }
      edgesOnVertex[3 * v + m] = e;
      edgeSignOnVertex[3 * v + m] = verticesOnEdge[2 * e + 1] === v ? 1 : -1;
    }
  }

  const xi = new THREE.Vector3();
  const xv = new THREE.Vector3();
  const before = new THREE.Vector3();
  const after = new THREE.Vector3();
  for (let i = 0; i < nCells; i++) {
    const n = nEdgesOnCell[i];
    getVector(xCell, i, xi);
    for (let k = 0; k < n; k++) {
      const v = verticesOnCell[MAX_EDGES * i + k];
      getVector(xVertex, v, xv);
      getVector(xEdge, edgesOnCell[MAX_EDGES * i + k], before);
      getVector(xEdge, edgesOnCell[MAX_EDGES * i + (k + 1) % n], after);
      const kite = (sphericalExcess(xi, before, xv) + sphericalExcess(xi, xv, after)) * R2;
      kiteAreasOnCell[MAX_EDGES * i + k] = kite;
      areaCell[i] += kite;
      for (let m = 0; m < 3; m++) {
        if (cellsOnVertex[3 * v + m] === i) {
          kiteAreasOnVertex[3 * v + m] = kite;
          areaTriangle[v] += kite;
        }
      }
    }
  }

  for (let e = 0; e < nEdges; e++) {
    for (let side = 0; side < 2; side++) {
      const c = cellsOnEdge[2 * e + side];
      const n = nEdgesOnCell[c];
      let start = -1;
      for (let k = 0; k < n; k++) {
        if (edgesOnCell[MAX_EDGES * c + k] === e) start = k;
      }
      const sign = edgeSignOnCell[MAX_EDGES * c + start];
      let passed = 0;
      for (let m = 1; m < n; m++) {
        passed += kiteAreasOnCell[MAX_EDGES * c + (start + m - 1) % n] / areaCell[c];
        const k = (start + m) % n;
        const slot = MAX_EDGES_ON_EDGE * e + nEdgesOnEdge[e]++;
        edgesOnEdge[slot] = edgesOnCell[MAX_EDGES * c + k];
        weightsOnEdge[slot] = sign * edgeSignOnCell[MAX_EDGES * c + k] * (0.5 - passed);
      }
    }
  }

  return {
    nCells, nEdges, nVertices, radius, omega,
    maxEdges: MAX_EDGES, maxEdgesOnEdge: MAX_EDGES_ON_EDGE,
    xCell, latCell, lonCell, areaCell, fCell, nEdgesOnCell,
    cellsOnCell, verticesOnCell, edgesOnCell, edgeSignOnCell, kiteAreasOnCell,
    xVertex, latVertex, lonVertex, areaTriangle, fVertex,
    cellsOnVertex, edgesOnVertex, edgeSignOnVertex, kiteAreasOnVertex,
    cellsOnEdge, verticesOnEdge, dcEdge, dvEdge, xEdge, nEdge, tEdge, latEdge, fEdge,
    nEdgesOnEdge, edgesOnEdge, weightsOnEdge,
  };
}
