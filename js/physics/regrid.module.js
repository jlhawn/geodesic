import { cellVector } from '../dynamics/operators.module.js';

/*
 * Barycentric weights of p in the plane through unit vectors a, b, c:
 * solve [a b c] w = p and normalize so the weights sum to one.
 */
function barycentric(a, b, c, p) {
  const det = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
  if (Math.abs(det) < 1e-18) return null;
  const wa = (p[0] * (b[1] * c[2] - b[2] * c[1]) - p[1] * (b[0] * c[2] - b[2] * c[0]) + p[2] * (b[0] * c[1] - b[1] * c[0])) / det;
  const wb = (a[0] * (p[1] * c[2] - p[2] * c[1]) - a[1] * (p[0] * c[2] - p[2] * c[0]) + a[2] * (p[0] * c[1] - p[1] * c[0])) / det;
  const wc = (a[0] * (b[1] * p[2] - b[2] * p[1]) - a[1] * (b[0] * p[2] - b[2] * p[0]) + a[2] * (b[0] * p[1] - b[1] * p[0])) / det;
  const sum = wa + wb + wc;
  return [wa / sum, wb / sum, wc / sum];
}

/*
 * Linear interpolation weights from a source mesh's cell centers onto a
 * set of unit-vector points: each point is located in the Delaunay
 * triangle (three cells around a source vertex) that contains it, found
 * among the triangles touching its nearest source cell.
 */
export function interpolationWeights(source, points) {
  const { nCells, xCell, nEdgesOnCell, verticesOnCell, cellsOnVertex, maxEdges } = source;
  const count = points.length / 3;
  const cells = new Int32Array(3 * count), weights = new Float64Array(3 * count);
  const vector = (i) => [xCell[3 * i], xCell[3 * i + 1], xCell[3 * i + 2]];
  for (let n = 0; n < count; n++) {
    const p = [points[3 * n], points[3 * n + 1], points[3 * n + 2]];
    let nearest = 0, best = -Infinity;
    for (let i = 0; i < nCells; i++) {
      const dot = p[0] * xCell[3 * i] + p[1] * xCell[3 * i + 1] + p[2] * xCell[3 * i + 2];
      if (dot > best) { best = dot; nearest = i; }
    }
    let chosen = null, chosenScore = -Infinity;
    for (let k = 0; k < nEdgesOnCell[nearest]; k++) {
      const v = verticesOnCell[maxEdges * nearest + k];
      const tri = [cellsOnVertex[3 * v], cellsOnVertex[3 * v + 1], cellsOnVertex[3 * v + 2]];
      const w = barycentric(vector(tri[0]), vector(tri[1]), vector(tri[2]), p);
      if (!w) continue;
      const score = Math.min(...w);
      if (score > chosenScore) { chosenScore = score; chosen = { tri, w }; }
    }
    const w = chosen.w.map((x) => Math.max(0, x));
    const sum = w[0] + w[1] + w[2];
    for (let m = 0; m < 3; m++) { cells[3 * n + m] = chosen.tri[m]; weights[3 * n + m] = w[m] / sum; }
  }
  return { cells, weights };
}

function apply(field, offset, { cells, weights }, out, outOffset, count) {
  for (let n = 0; n < count; n++) {
    out[outOffset + n] = weights[3 * n] * field[offset + cells[3 * n]] + weights[3 * n + 1] * field[offset + cells[3 * n + 1]] + weights[3 * n + 2] * field[offset + cells[3 * n + 2]];
  }
}

/*
 * Carries a state [pi, theta, u, surfaceT] from one model to another with
 * the same sigma levels: scalars are interpolated to the target cell
 * centers, and the cell-center wind vectors of each layer are
 * interpolated to the target edge midpoints and projected onto the edge
 * normals.
 */
export function regridState(source, target, state) {
  const [pi, theta, u, surfaceT, q = null, qc = null] = state;
  const K = source.core.K;
  if (K !== target.core.K) throw new Error(`layer counts differ: ${K} vs ${target.core.K}`);
  const sm = source.mesh, tm = target.mesh;
  const atCells = interpolationWeights(sm, tm.xCell);
  const atEdges = interpolationWeights(sm, tm.xEdge);
  const outPi = new Float64Array(tm.nCells), outTheta = new Float64Array(K * tm.nCells), outU = new Float64Array(K * tm.nEdges), outSurfaceT = new Float64Array(tm.nCells);
  const outQ = q ? new Float64Array(K * tm.nCells) : null;
  const outQc = qc ? new Float64Array(K * tm.nCells) : null;
  apply(pi, 0, atCells, outPi, 0, tm.nCells);
  apply(surfaceT, 0, atCells, outSurfaceT, 0, tm.nCells);
  const vector = new Float64Array(3 * sm.nCells);
  const component = new Float64Array(sm.nCells);
  const edgeComponent = new Float64Array(tm.nEdges);
  for (let k = 0; k < K; k++) {
    apply(theta, k * sm.nCells, atCells, outTheta, k * tm.nCells, tm.nCells);
    if (q) apply(q, k * sm.nCells, atCells, outQ, k * tm.nCells, tm.nCells);
    if (qc) apply(qc, k * sm.nCells, atCells, outQc, k * tm.nCells, tm.nCells);
    cellVector(sm, u.subarray(k * sm.nEdges, (k + 1) * sm.nEdges), vector);
    for (let axis = 0; axis < 3; axis++) {
      for (let i = 0; i < sm.nCells; i++) component[i] = vector[3 * i + axis];
      apply(component, 0, atEdges, edgeComponent, 0, tm.nEdges);
      for (let e = 0; e < tm.nEdges; e++) outU[k * tm.nEdges + e] += edgeComponent[e] * tm.nEdge[3 * e + axis];
    }
  }
  const out = [outPi, outTheta, outU, outSurfaceT];
  for (const tracer of [outQ, outQc]) {
    if (!tracer) continue;
    for (let x = 0; x < tracer.length; x++) if (tracer[x] < 0) tracer[x] = 0;
    out.push(tracer);
  }
  return out;
}
