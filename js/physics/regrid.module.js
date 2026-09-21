import { cellVector } from '../dynamics/operators.module.js';
import { EARTH_RADIUS } from '../mesh.module.js';

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
export function interpolationWeights(source, points, progress = null) {
  const { nCells, xCell, nEdgesOnCell, verticesOnCell, cellsOnVertex, cellsOnCell, maxEdges } = source;
  const count = points.length / 3;
  const cells = new Int32Array(3 * count), weights = new Float64Array(3 * count), tiles = new Int32Array(count);
  const vector = (i) => [xCell[3 * i], xCell[3 * i + 1], xCell[3 * i + 2]];
  const dotWith = (p, i) => p[0] * xCell[3 * i] + p[1] * xCell[3 * i + 1] + p[2] * xCell[3 * i + 2];
  const bruteForce = (p) => { let nearest = 0, best = -Infinity; for (let i = 0; i < nCells; i++) { const d = dotWith(p, i); if (d > best) { best = d; nearest = i; } } return nearest; };
  const walk = (p, start) => {
    let here = start, best = dotWith(p, here);
    for (;;) {
      let next = here;
      for (let k = 0; k < nEdgesOnCell[here]; k++) { const j = cellsOnCell[maxEdges * here + k]; const d = dotWith(p, j); if (d > best) { best = d; next = j; } }
      if (next === here) return here;
      here = next;
    }
  };
  const every = Math.max(1, Math.floor(count / 50));
  let previous = 0;
  for (let n = 0; n < count; n++) {
    if (progress && n % every === 0) progress(n / count);
    const p = [points[3 * n], points[3 * n + 1], points[3 * n + 2]];
    let nearest = walk(p, previous);
    let chosen = null, chosenScore = -Infinity;
    for (let attempt = 0; attempt < 2 && !chosen; attempt++) {
      if (attempt === 1) nearest = bruteForce(p);
    for (let k = 0; k < nEdgesOnCell[nearest]; k++) {
      const v = verticesOnCell[maxEdges * nearest + k];
      const tri = [cellsOnVertex[3 * v], cellsOnVertex[3 * v + 1], cellsOnVertex[3 * v + 2]];
      const w = barycentric(vector(tri[0]), vector(tri[1]), vector(tri[2]), p);
      if (!w) continue;
      const score = Math.min(...w);
      if (score > chosenScore) { chosenScore = score; chosen = { tri, w }; }
    }
    }
    previous = nearest;
    tiles[n] = nearest;
    const w = chosen.w.map((x) => Math.max(0, x));
    const sum = w[0] + w[1] + w[2];
    for (let m = 0; m < 3; m++) { cells[3 * n + m] = chosen.tri[m]; weights[3 * n + m] = w[m] / sum; }
  }
  return { cells, weights, tiles };
}

/*
 * The source cell nearest to `start` that the mask admits, searched ring
 * by ring through the neighbours; -1 when none lies within `rings`.
 */
function nearestAdmitted(mesh, start, mask, rings) {
  if (mask[start]) return start;
  const { cellsOnCell, nEdgesOnCell, maxEdges } = mesh;
  let frontier = [start];
  const seen = new Set(frontier);
  for (let ring = 0; ring < rings && frontier.length; ring++) {
    const next = [];
    for (const i of frontier) {
      for (let k = 0; k < nEdgesOnCell[i]; k++) {
        const j = cellsOnCell[maxEdges * i + k];
        if (seen.has(j)) continue;
        if (mask[j]) return j;
        seen.add(j);
        next.push(j);
      }
    }
    frontier = next;
  }
  return -1;
}

const ringsWithin = (mesh, reach) => Math.ceil(reach / (EARTH_RADIUS * Math.sqrt(4 * Math.PI / mesh.nCells)));

function apply(field, offset, { cells, weights }, out, outOffset, count) {
  for (let n = 0; n < count; n++) {
    out[outOffset + n] = weights[3 * n] * field[offset + cells[3 * n]] + weights[3 * n + 1] * field[offset + cells[3 * n + 1]] + weights[3 * n + 2] * field[offset + cells[3 * n + 2]];
  }
}

/*
 * Interpolation that only draws on source cells the mask admits: the
 * weights of the others are dropped and the rest renormalized, and a
 * point whose whole triangle is excluded takes the nearest admitted
 * cell's value.
 */
function applyMasked(mesh, field, { cells, weights, tiles }, mask, out, count, reach = 1e6) {
  const rings = ringsWithin(mesh, reach);
  for (let n = 0; n < count; n++) {
    let sum = 0, value = 0, plain = 0;
    for (let m = 0; m < 3; m++) {
      const j = cells[3 * n + m];
      plain += weights[3 * n + m] * field[j];
      if (!mask[j]) continue;
      sum += weights[3 * n + m];
      value += weights[3 * n + m] * field[j];
    }
    if (sum > 0) { out[n] = value / sum; continue; }
    const near = nearestAdmitted(mesh, tiles[n], mask, rings);
    out[n] = near >= 0 ? field[near] : plain;
  }
}

/*
 * Every target point takes the value of the source tile it lies in, or
 * of the nearest admitted tile within `reach` when a mask is given, and
 * `fill` when there is none that close. Ice, snow and soil are carried
 * this way so that an ice edge or a snow line stays where the coarser
 * run had it instead of being smeared across the neighbours, so that
 * sea values never leak onto land or land values into the sea, and so
 * that land the finer mesh resolves for the first time starts with the
 * state of the coast next to it.
 */
export function sampleTiles(source, target, field, mask = null, weights = interpolationWeights(source.mesh, target.mesh.xCell), { reach = 1e6, fill = 0 } = {}) {
  const { tiles } = weights;
  const rings = ringsWithin(source.mesh, reach);
  const out = new Float64Array(target.mesh.nCells);
  for (let n = 0; n < out.length; n++) {
    const tile = mask ? nearestAdmitted(source.mesh, tiles[n], mask, rings) : tiles[n];
    out[n] = tile >= 0 ? field[tile] : fill;
  }
  return out;
}

const seaMask = (model) => (model.geography ? Uint8Array.from(model.geography.land, (l) => 1 - l) : null);
const landMask = (model) => (model.geography ? model.geography.land : null);

/*
 * Carries a state [pi, theta, u, surfaceT] from one model to another with
 * the same sigma levels: scalars are interpolated to the target cell
 * centers, and the cell-center wind vectors of each layer are
 * interpolated to the target edge midpoints and projected onto the edge
 * normals.
 */
export function regridCellField(source, target, field, weights = interpolationWeights(source.mesh, target.mesh.xCell), mask = null) {
  const out = new Float64Array(target.mesh.nCells);
  if (mask) applyMasked(source.mesh, field, weights, mask, out, target.mesh.nCells);
  else apply(field, 0, weights, out, 0, target.mesh.nCells);
  return out;
}

export function regridEdgeField(source, target, u, weights = interpolationWeights(source.mesh, target.mesh.xEdge)) {
  const sm = source.mesh, tm = target.mesh;
  const out = new Float64Array(tm.nEdges);
  const vector = cellVector(sm, u), component = new Float64Array(sm.nCells), edgeComponent = new Float64Array(tm.nEdges);
  for (let axis = 0; axis < 3; axis++) {
    for (let i = 0; i < sm.nCells; i++) component[i] = vector[3 * i + axis];
    apply(component, 0, weights, edgeComponent, 0, tm.nEdges);
    for (let e = 0; e < tm.nEdges; e++) out[e] += edgeComponent[e] * tm.nEdge[3 * e + axis];
  }
  return out;
}

export function regridOcean(source, target, ocean, progress = null) {
  if (source.mesh.nCells === target.mesh.nCells) return Object.fromEntries(Object.entries(ocean).map(([k, v]) => [k, Float64Array.from(v)]));
  if (progress) progress(0, 'the ocean');
  const atCells = interpolationWeights(source.mesh, target.mesh.xCell), atEdges = interpolationWeights(source.mesh, target.mesh.xEdge);
  const sea = seaMask(source);
  const cell = (v) => regridCellField(source, target, Float64Array.from(v), atCells, sea), edge = (v) => regridEdgeField(source, target, Float64Array.from(v), atEdges);
  return { h1: cell(ocean.h1), h2: cell(ocean.h2), u1: edge(ocean.u1), u2: edge(ocean.u2), T2: cell(ocean.T2) };
}

export function regridLand(source, target, land, progress = null) {
  if (source.mesh.nCells === target.mesh.nCells) return { soil: Float64Array.from(land.soil), snow: Float64Array.from(land.snow) };
  if (progress) progress(0, 'the land');
  const atCells = interpolationWeights(source.mesh, target.mesh.xCell), onLand = landMask(source);
  const halfBucket = 0.5 * (target.land && target.land.bucketCapacity ? target.land.bucketCapacity : 150);
  return { soil: sampleTiles(source, target, Float64Array.from(land.soil), onLand, atCells, { fill: halfBucket }), snow: sampleTiles(source, target, Float64Array.from(land.snow), onLand, atCells) };
}

export function regridState(source, target, state, progress = null) {
  const [pi, theta, u, surfaceT, q = null, qc = null, ice = null] = state;
  const K = source.core.K;
  if (K !== target.core.K) throw new Error(`layer counts differ: ${K} vs ${target.core.K}`);
  const sm = source.mesh, tm = target.mesh;
  const report = (fraction, text) => { if (progress) progress(fraction, text); };
  const atCells = interpolationWeights(sm, tm.xCell, (f) => report(0.1 * f, `locating the cells, ${Math.round(100 * f)}%`));
  const atEdges = interpolationWeights(sm, tm.xEdge, (f) => report(0.1 + 0.2 * f, `locating the edges, ${Math.round(100 * f)}%`));
  const outPi = new Float64Array(tm.nCells), outTheta = new Float64Array(K * tm.nCells), outU = new Float64Array(K * tm.nEdges), outSurfaceT = new Float64Array(tm.nCells);
  const outQ = q ? new Float64Array(K * tm.nCells) : null;
  const outQc = qc ? new Float64Array(K * tm.nCells) : null;
  const outIce = ice ? sampleTiles(source, target, ice, seaMask(source), atCells) : null;
  apply(pi, 0, atCells, outPi, 0, tm.nCells);
  apply(surfaceT, 0, atCells, outSurfaceT, 0, tm.nCells);
  const vector = new Float64Array(3 * sm.nCells);
  const component = new Float64Array(sm.nCells);
  const edgeComponent = new Float64Array(tm.nEdges);
  for (let k = 0; k < K; k++) {
    report(0.3 + 0.7 * k / K, `layer ${k + 1} of ${K}`);
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
  for (const tracer of [outQ, outQc, outIce]) {
    if (!tracer) continue;
    for (let x = 0; x < tracer.length; x++) if (tracer[x] < 0) tracer[x] = 0;
    out.push(tracer);
  }
  return out;
}
