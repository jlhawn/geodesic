import { gradient, divergence, tangential } from '../js/dynamics/operators.module.js';

const SQRT3 = Math.sqrt(3);
const RING = [[[1, 0], [0, -1], [-1, -1], [-1, 0], [-1, 1], [0, 1]], [[1, 0], [1, -1], [0, -1], [-1, 0], [0, 1], [1, 1]]];

export function hexTorus(cols, rows, spacing) {
  if (cols % 3 || rows % 2) throw new Error(`a ${cols} × ${rows} patch cannot wrap the three colors`);
  const nCells = cols * rows, nEdges = 3 * nCells, maxEdgesOnEdge = 10;
  const place = (i) => { const c = i % cols; return [c, (i - c) / cols]; };
  const neighbor = (i, k) => { const [c, r] = place(i), [dc, dr] = RING[r & 1][k]; return ((r + dr + rows) % rows) * cols + (c + dc + cols) % cols; };
  const xy = Array.from({ length: nCells }, (_, i) => { const [c, r] = place(i); return [spacing * (c + (r & 1) / 2), -spacing * SQRT3 / 2 * r]; });
  const color = Int8Array.from({ length: nCells }, (_, i) => { const [c, r] = place(i); return (((c - (r >> 1) - r) % 3) + 3) % 3; });
  const nEdgesOnCell = new Int32Array(nCells).fill(6), edgesOnCell = new Int32Array(6 * nCells), edgeSignOnCell = new Int32Array(6 * nCells);
  const cellsOnEdge = new Int32Array(2 * nEdges), nEdge = new Float64Array(3 * nEdges);
  for (let i = 0; i < nCells; i++) for (let k = 0; k < 6; k++) {
    const own = k < 3, e = own ? 3 * i + k : 3 * neighbor(i, k) + k - 3;
    edgesOnCell[6 * i + k] = e;
    edgeSignOnCell[6 * i + k] = own ? 1 : -1;
    if (!own) continue;
    cellsOnEdge[2 * e] = i; cellsOnEdge[2 * e + 1] = neighbor(i, k);
    nEdge[3 * e] = Math.cos(k * Math.PI / 3); nEdge[3 * e + 1] = Math.sin(k * Math.PI / 3);
  }
  const nEdgesOnEdge = new Int32Array(nEdges), edgesOnEdge = new Int32Array(maxEdgesOnEdge * nEdges), weightsOnEdge = new Float64Array(maxEdgesOnEdge * nEdges);
  for (let e = 0; e < nEdges; e++) for (let side = 0; side < 2; side++) {
    const i = cellsOnEdge[2 * e + side], start = [0, 1, 2, 3, 4, 5].find((k) => edgesOnCell[6 * i + k] === e), sign = edgeSignOnCell[6 * i + start];
    for (let m = 1; m < 6; m++) {
      const k = (start + m) % 6, slot = maxEdgesOnEdge * e + nEdgesOnEdge[e]++;
      edgesOnEdge[slot] = edgesOnCell[6 * i + k];
      weightsOnEdge[slot] = sign * edgeSignOnCell[6 * i + k] * (0.5 - m / 6);
    }
  }
  const mesh = {
    nCells, nEdges, maxEdges: 6, maxEdgesOnEdge, nEdgesOnCell, edgesOnCell, edgeSignOnCell, cellsOnEdge, nEdge, nEdgesOnEdge, edgesOnEdge, weightsOnEdge,
    dcEdge: new Float64Array(nEdges).fill(spacing), dvEdge: new Float64Array(nEdges).fill(spacing / SQRT3), areaCell: new Float64Array(nCells).fill(SQRT3 / 2 * spacing ** 2),
  };
  return { cols, rows, spacing, mesh, xy, color, neighbor };
}

export function threeColors({ color }, amplitude = 1) {
  return Float64Array.from(color, (k) => amplitude * [1, 0, -1][k]);
}

export function broadHigh({ cols, rows, spacing, xy }, radius, amplitude = 1) {
  const width = cols * spacing, height = rows * spacing * SQRT3 / 2, [cx, cy] = xy[(rows >> 1) * cols + (cols >> 1)];
  const wrap = (d, span) => d - span * Math.round(d / span);
  const p = Float64Array.from(xy, ([x, y]) => amplitude * Math.exp(-(wrap(x - cx, width) ** 2 + wrap(y - cy, height) ** 2) / radius ** 2));
  const mean = p.reduce((a, b) => a + b, 0) / p.length;
  return p.map((v) => v - mean);
}

export function linearWater(mesh, { c, f, rho, cd = false }) {
  const { nCells, nEdges } = mesh;
  const push = new Float64Array(nEdges), along = new Float64Array(nEdges), pushAlong = new Float64Array(nEdges), spread = new Float64Array(nCells);
  const blank = () => ({ p: new Float64Array(nCells), u: new Float64Array(nEdges), v: new Float64Array(nEdges) });
  const stages = Array.from({ length: 4 }, blank), trial = blank();
  function rate(p, u, v, out) {
    gradient(mesh, p, push);
    divergence(mesh, u, spread);
    if (cd) tangential(mesh, push, pushAlong);
    else tangential(mesh, u, along);
    for (let e = 0; e < nEdges; e++) {
      out.u[e] = -push[e] / rho + f * (cd ? v[e] : along[e]);
      if (cd) out.v[e] = -f * u[e] - pushAlong[e] / rho;
    }
    for (let i = 0; i < nCells; i++) out.p[i] = -rho * c * c * spread[i];
  }
  return function step(state, dt) {
    const v = state.v ??= new Float64Array(nEdges);
    rate(state.p, state.u, v, stages[0]);
    for (let n = 1; n < 4; n++) {
      const h = n < 3 ? dt / 2 : dt, last = stages[n - 1];
      for (let i = 0; i < nCells; i++) trial.p[i] = state.p[i] + h * last.p[i];
      for (let e = 0; e < nEdges; e++) { trial.u[e] = state.u[e] + h * last.u[e]; trial.v[e] = v[e] + h * last.v[e]; }
      rate(trial.p, trial.u, trial.v, stages[n]);
    }
    const [a, b, d, g] = stages;
    for (let i = 0; i < nCells; i++) state.p[i] += dt / 6 * (a.p[i] + 2 * b.p[i] + 2 * d.p[i] + g.p[i]);
    for (let e = 0; e < nEdges; e++) {
      state.u[e] += dt / 6 * (a.u[e] + 2 * b.u[e] + 2 * d.u[e] + g.u[e]);
      v[e] += dt / 6 * (a.v[e] + 2 * b.v[e] + 2 * d.v[e] + g.v[e]);
    }
  };
}
