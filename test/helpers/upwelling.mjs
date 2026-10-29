import { cellUpwelling } from '../../js/ocean/layered.module.js';

// The upward velocity at the base of every layer of every sea cell, L × C, from the model's mass fluxes at the state {h, u, eta, sillEdge, edgeOcean}.
export function layerVelocities(mesh, L, state, cellOcean) {
  const C = mesh.nCells, out = new Float64Array(L * C), column = new Float64Array(L);
  for (let i = 0; i < C; i++) {
    if (!cellOcean[i]) continue;
    cellUpwelling(mesh, L, state, i, column);
    for (let k = 0; k < L; k++) out[k * C + i] = column[k];
  }
  return out;
}

// The velocities V (L × C, at the layers' bases) at `upwellingDepth`, linear in depth between the top and the base of the layer that holds it (zero at the surface); NaN on land and where no layer holds `depth`, the value at the floor where only `upwellingDepth` lies below it.
export function velocityAt(C, L, h, cellOcean, V, depth, upwellingDepth = depth) {
  const out = new Float64Array(C);
  for (let i = 0; i < C; i++) {
    if (!cellOcean[i]) { out[i] = NaN; continue; }
    let top = 0, previous = 0, held = false, value = null;
    for (let k = 0; k < L; k++) {
      const hk = h[k * C + i], here = V[k * C + i];
      if (depth < top + hk) held = true;
      if (value === null && upwellingDepth < top + hk) value = previous + (upwellingDepth - top) / hk * (here - previous);
      top += hk; previous = here;
    }
    out[i] = !held ? NaN : value === null ? previous : value;
  }
  return out;
}

// The largest difference between a and b over the cells where a is a number, with where it falls; the two must agree on which cells are NaN.
export function worstDifference(a, b) {
  let max = 0, at = -1, scale = 0, sumSquares = 0, sumReference = 0, n = 0;
  for (let i = 0; i < a.length; i++) {
    if (Number.isNaN(a[i]) || Number.isNaN(b[i])) {
      if (Number.isNaN(a[i]) !== Number.isNaN(b[i])) return { max: Infinity, at: i, scale, rms: Infinity };
      continue;
    }
    const d = Math.abs(a[i] - b[i]);
    if (d > max) { max = d; at = i; }
    scale = Math.max(scale, Math.abs(a[i]));
    sumSquares += d * d; sumReference += a[i] * a[i]; n++;
  }
  return { max, at, scale, rms: Math.sqrt(sumSquares / Math.max(n, 1)), rmsScale: Math.sqrt(sumReference / Math.max(n, 1)) };
}

// Makes the top layer of every third sea cell 40% thicker at the expense of its thickest other layer, column sums kept, so that neighboring top layers differ and the donor cell's thickness limits some of the edge fluxes.
export function unevenTopLayers(h, C, L, cellOcean) {
  for (let i = 0; i < C; i += 3) {
    if (!cellOcean[i]) continue;
    let thick = 1;
    for (let k = 2; k < L; k++) if (h[k * C + i] > h[thick * C + i]) thick = k;
    const extra = 0.4 * h[i];
    h[thick * C + i] -= extra;
    h[i] += extra;
  }
}
