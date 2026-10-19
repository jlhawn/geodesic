/*
 * The water and heat the hybrid ocean moves through one edge, layer by layer,
 * as its tendency computes them from a state (h, u, T, eta): the edge
 * thickness of each layer (layer 0 the mean of its two cells, every other
 * layer the smaller, the column clipped from the top at the sill depth plus
 * the mean free surface), the thickness flux donor-limited, the heat it
 * carries at the face value of the ocean's tracer scheme (the limited
 * third-order value from the upwind cell, the cell beyond it and the
 * downwind cell, or the upwind value where the scheme falls back to it), and
 * the lateral diffusion of heat, κ = lateralDiffusionSpeed·d plus, in the
 * fixed-depth layers, the shear-scaled term min(cap, coefficient·|ζ|·d²)
 * with ζ the mean of the edge's two vertex vorticities. Volumes are m³/s
 * and heats K·m³/s from cell a = cellsOnEdge[2e] toward cell b, T in kelvin
 * as given.
 */
import { curl } from '../js/dynamics/operators.module.js';
import { beyondCells, limitedEdgeValue, EPS } from '../js/ocean/layered.module.js';

export function edgeLayerFluxes(mesh, { layers: L, fixed, options, cellOcean, depth, sill = null }) {
  const { nCells: C, nEdges: E, nVertices: V, cellsOnEdge, dcEdge, dvEdge, verticesOnEdge } = mesh;
  const beyond = options.tracerAdvection === 'third' ? beyondCells(mesh) : null;
  const lateral = options.lateralDiffusionSpeed ?? 0, shear = options.shearDiffusivity ?? 0, cap = options.shearDiffusivityCap ?? Infinity;
  const zeta = Array.from({ length: shear > 0 ? fixed : 0 }, () => new Float64Array(V)), column = new Float64Array(L);
  let state = null;
  return {
    load(next) {
      state = next;
      for (let k = 0; k < zeta.length; k++) curl(mesh, next.u.subarray(k * E, (k + 1) * E), zeta[k]);
    },
    edge(e, out) {
      const { h, u, T, eta } = state, a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1], dc = dcEdge[e], dv = dvEdge[e];
      column[0] = 0.5 * (h[a] + h[b]);
      for (let k = 1; k < L; k++) column[k] = Math.min(h[k * C + a], h[k * C + b]);
      const floor = Math.max(EPS, (sill ? sill[e] : Math.min(depth[a], depth[b])) + 0.5 * (eta[a] + eta[b]));
      let above = 0;
      for (let k = 0; k < L; k++) { const hk = column[k]; column[k] = Math.min(hk, Math.max(0, floor - above)); above += hk; }
      for (let k = 0; k < L; k++) {
        const o = k * C, ue = u[k * E + e], side = ue > 0 ? 0 : 1, up = side ? b : a, down = side ? a : b;
        const flux = Math.min(column[k], Math.max(0, h[o + up])) * ue;
        let face = T[o + up];
        if (beyond && flux !== 0) {
          const uu = beyond.cell[2 * e + side];
          if (uu >= 0 && cellOcean[uu] && h[o + uu] > 1e-6 && h[o + down] > 1e-6) face = limitedEdgeValue(T.subarray(o, o + C), up, down, uu, dcEdge[beyond.edge[2 * e + side]], dc);
        }
        const kappa = lateral * dc + (k < zeta.length ? Math.min(cap, shear * Math.abs(0.5 * (zeta[k][verticesOnEdge[2 * e]] + zeta[k][verticesOnEdge[2 * e + 1]])) * dc * dc) : 0);
        out.volume[k] = dv * flux;
        out.net[k] = dv * column[k] * ue;
        out.heat[k] = dv * flux * face;
        out.diffusion[k] = -dv * kappa * column[k] / dc * (T[o + b] - T[o + a]);
        out.donor[k] = up;
      }
      return out;
    },
  };
}
