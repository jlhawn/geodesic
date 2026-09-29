import { Grid } from '../../js/grid.module.js';
import { buildMesh } from '../../js/mesh.module.js';
import { cellVector } from '../../js/dynamics/operators.module.js';
import { createOcean } from '../../js/ocean/layered.module.js';
import { labelTemperature } from '../../js/ocean/seawater.module.js';

export const RHO_AIR = 1.2, DRAG = 1.5e-3, RHO = 1025, DEG = Math.PI / 180;
export const mesh = buildMesh(new Grid(8));
export const { nCells: C, nEdges: E } = mesh;

export function zonalWindOnEdges(m, speedAt) {
  const u = new Float64Array(m.nEdges);
  for (let e = 0; e < m.nEdges; e++) {
    const x = m.xEdge[3 * e], y = m.xEdge[3 * e + 1], r = Math.hypot(x, y);
    const east = r > 0 ? [-y / r, x / r, 0] : [0, 0, 0];
    u[e] = speedAt(m.latEdge[e]) * (east[0] * m.nEdge[3 * e] + east[1] * m.nEdge[3 * e + 1] + east[2] * m.nEdge[3 * e + 2]);
  }
  return u;
}

// Area-weighted sums of h·T and h·S over every layer: what the module's own
// flux-divergence advection and mixing conserve exactly, independent of the
// rhoCp-scaled diagnostics().oceanHeat proxy.
export function totalHeatSalt(ocean, m) {
  let heat = 0, salt = 0;
  for (let i = 0; i < m.nCells; i++) {
    if (!ocean.cellOcean[i]) continue;
    for (let k = 0; k < ocean.layers; k++) { heat += m.areaCell[i] * ocean.Q[k * m.nCells + i]; salt += m.areaCell[i] * ocean.W[k * m.nCells + i]; }
  }
  return { heat, salt };
}

export function northwardMixedTransport(ocean, m) {
  const vector = cellVector(m, ocean.u.subarray(0, m.nEdges));
  const out = new Float64Array(m.nCells);
  for (let i = 0; i < m.nCells; i++) {
    const x = m.xCell[3 * i], y = m.xCell[3 * i + 1], z = m.xCell[3 * i + 2], r = Math.hypot(x, y);
    const north = [-z * x / r, -z * y / r, r];
    out[i] = ocean.h[i] * (vector[3 * i] * north[0] + vector[3 * i + 1] * north[1] + vector[3 * i + 2] * north[2]);
  }
  return out;
}

// A flat ocean the same in every column, 34.5 psu water at 278 K over the
// interior classes denser than it, with a 50 m mixed layer: no pressure
// gradient anywhere until a column is changed.
export function uniformOcean(options = {}, m = mesh) {
  const ocean = createOcean(m, { everySteps: 1, thermoclineTilt: 0, salinityProfile: () => 34.5, mixedDepth: 50, ...options });
  const surfaceT = new Float64Array(m.nCells).fill(278), ice = new Float64Array(m.nCells), flux = new Float64Array(m.nCells);
  ocean.initialize(surfaceT, ice);
  return { ocean, surfaceT, ice, flux, calm: new Float64Array(m.nEdges) };
}

// Mixes each cell's interior into its mixed layer from the top down until the
// layer is `depth` deep, gives it the density `offset` from the first interior
// layer left beneath it, and reloads, so the surface temperature and the
// mixed layer's last temperature are that water's.
export function deepen({ ocean, surfaceT, ice }, cells, depth, offset) {
  const { h, Q, W, densities, layers: L } = ocean, n = h.length / L;
  for (const i of cells) {
    for (let k = 1; k < L && h[i] < depth; k++) {
      const a = k * n + i, take = Math.min(h[a] - 0.01, depth - h[i]);
      if (take <= 0) continue;
      const f = take / h[a];
      h[i] += take; Q[i] += Q[a] * f; W[i] += W[a] * f;
      h[a] -= take; Q[a] -= Q[a] * f; W[a] -= W[a] * f;
    }
    let below = 1;
    while (below < L - 1 && h[below * n + i] <= 5) below++;
    const t = labelTemperature(densities[below] + offset, W[i] / h[i]);
    Q[i] = h[i] * t; surfaceT[i] = t;
  }
  ocean.load(ocean.serialize(), surfaceT, ice);
}

export function buriedBump(options) {
  const ocean = createOcean(mesh, { bathymetry: new Float64Array(C).fill(4000), thermoclineTilt: 0, salinityProfile: () => 35, ...options });
  ocean.initialize(new Float64Array(C).fill(290), new Float64Array(C));
  const centre = [Math.cos(-45 * DEG), 0, Math.sin(-45 * DEG)];
  const reshape = (k, amount) => {
    for (let i = 0; i < C; i++) {
      const x = mesh.xCell[3 * i], y = mesh.xCell[3 * i + 1], z = mesh.xCell[3 * i + 2], r = Math.hypot(x, y, z);
      const distance = Math.acos(Math.min(1, (x * centre[0] + y * centre[1] + z * centre[2]) / r)) * mesh.radius;
      const n = k * C + i, t = ocean.Q[n] / ocean.h[n], s = ocean.W[n] / ocean.h[n];
      ocean.h[n] += amount * Math.exp(-((distance / 2e6) ** 2));
      ocean.Q[n] = ocean.h[n] * t; ocean.W[n] = ocean.h[n] * s;
    }
  };
  return { ocean, reshape };
}

export function interfaceRange(ocean, k) {
  let low = Infinity, high = -Infinity;
  for (let i = 0; i < C; i++) {
    let z = 0;
    for (let j = 0; j <= k; j++) z += ocean.h[j * C + i];
    low = Math.min(low, z); high = Math.max(high, z);
  }
  return high - low;
}

export function classContents(ocean) {
  return Array.from({ length: ocean.layers }, (_, k) => {
    let volume = 0, heat = 0, salt = 0;
    for (let i = 0; i < C; i++) { const n = k * C + i; volume += mesh.areaCell[i] * ocean.h[n]; heat += mesh.areaCell[i] * ocean.Q[n]; salt += mesh.areaCell[i] * ocean.W[n]; }
    return { volume, heat, salt };
  });
}
