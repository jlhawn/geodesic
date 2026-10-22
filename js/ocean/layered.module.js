import { divergence, gradient, curl, kineticEnergy, laplacianVelocity, cellVector } from '../dynamics/operators.module.js';
import { FREEZING_POINT } from '../physics/ice.module.js';
import { SEAWATER, seawaterDensity, thermalExpansion as expansionOf, halineContraction as contractionOf, labelTemperature, salinityForDensity } from './seawater.module.js';
import { profileAt } from './climatology.module.js';
import { FRICTION_PER_WIND, CALM_WIND, RETREAT_EKMAN, ROTATION_LIMITS, LANGMUIR_ICE, ICE_WAVE_DAMPING, erfc, stokesSurfaceLayer, langmuirFactor, JERLOV, PENETRATION_DEPTH, transmitted, deposited, HYBRID_LAYERS, H_MIN, MODE_DEPTH, labelBoundaries, columnScratch, targetInterfaces, remapColumn, solveColumn, mixColumn, mixedLayerDepth, spaceInterfaces, layerDensities } from './column.module.js';
import { basalMelt, baseFreezing, SHELF_MELT } from './iceShelfMelt.module.js';

/*
 * The ocean on the C-grid: the ALE hybrid vertical coordinate of
 * docs/ocean-hybrid-coordinate.md, HYBRID_LAYERS fixed-depth surface
 * layers over isopycnal classes that hold their minimum thickness where
 * their water is absent, on the real bathymetry, with the density of
 * seawater/seawater.module.js. Each layer is a TRiSK shallow-water layer
 * in the vector-invariant form, carrying thickness, edge velocity, heat
 * h·T and salt h·S; the layers move with the flow through one RK4 step,
 * the free surface sub-stepped, are regridded to the generator's targets
 * and remapped conservatively, and then take the column physics
 * (column.module.js). Every column's thicknesses sum to its depth plus
 * the free surface; land columns are zero.
 *
 * The interior classes of the saved slab oceans, placed where the World Ocean Atlas holds its water
 * in this equation of state: 1019.5–1021.5 for the warm pool, the two
 * lightest at WOA23's tropical surface salinities, 0.25 apart
 * through the thermocline, 0.1 apart for the intermediate water to 1026.6,
 * 0.03 apart from 1026.65 to 1026.98 so that the Southern Ocean's winter
 * water (1026.79), Circumpolar Deep Water (1026.85), North Atlantic Deep
 * Water (1026.89) and Antarctic Bottom Water (1026.93) lie in separate
 * classes, and a few steps for the Arctic's deep and shelf water and the
 * Red Sea's and Mediterranean's. LAYER_SALINITIES label each class; a saved
 * slab ocean is converted onto the hybrid layers on load (convertSlab).
 */
export const LAYER_DENSITIES = [
  1019.5, 1020.0, 1020.5, 1021.0, 1021.5,
  1022.0, 1022.25, 1022.5, 1022.75, 1023.0, 1023.25, 1023.5, 1023.75, 1024.0, 1024.25, 1024.5, 1024.75, 1025.0, 1025.25, 1025.5, 1025.75, 1026.0,
  1026.1, 1026.2, 1026.3, 1026.4, 1026.5, 1026.6,
  1026.65, 1026.68, 1026.71, 1026.74, 1026.77, 1026.8, 1026.83, 1026.86, 1026.89, 1026.92, 1026.95, 1026.98,
  1027.05, 1027.15, 1027.3, 1027.5, 1027.75, 1028.0,
];
export const LAYER_SALINITIES = [
  33.0, 33.5, 34.0, 34.5, 34.85,
  35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 34.98, 34.95, 34.92, 34.9,
  34.89, 34.88, 34.875, 34.87, 34.86, 34.85,
  34.843, 34.839, 34.835, 34.831, 34.826, 34.82, 34.814, 34.809, 34.806, 34.803, 34.8, 34.82,
  34.91, 34.94, 35.05, 39.9, 39.0, 38.65,
];
/*
 * The interior classes of saved oceans that carry no list of their own,
 * told apart by their layer count: the seven of the page's saved runs,
 * the twenty-three of later ones and the forty-four from 1020.5.
 */
export const UNLISTED_LAYER_DENSITIES = [
  [1022.0, 1023.0, 1024.0, 1025.0, 1026.0, 1026.6, 1026.95],
  [1022.0, 1022.25, 1022.5, 1022.75, 1023.0, 1023.25, 1023.5, 1023.75, 1024.0, 1024.25, 1024.5, 1024.75, 1025.0, 1025.25, 1025.5, 1025.75, 1026.0, 1026.2, 1026.4, 1026.6, 1026.75, 1026.85, 1026.95],
  [
    1020.5, 1021.0, 1021.5,
    1022.0, 1022.25, 1022.5, 1022.75, 1023.0, 1023.25, 1023.5, 1023.75, 1024.0, 1024.25, 1024.5, 1024.75, 1025.0, 1025.25, 1025.5, 1025.75, 1026.0,
    1026.1, 1026.2, 1026.3, 1026.4, 1026.5, 1026.6,
    1026.65, 1026.68, 1026.71, 1026.74, 1026.77, 1026.8, 1026.83, 1026.86, 1026.89, 1026.92, 1026.95, 1026.98,
    1027.05, 1027.15, 1027.3, 1027.5, 1027.75, 1028.0,
  ],
];
export const THERMOCLINE_DENSITY = 1024.0;

/*
 * The fields at a depth below the surface, per cell: the layer holding
 * that depth (the first whose base lies below it; interior layers no
 * thicker than `thin`, 0 by default, are passed over), its temperature, its
 * current as a cell vector and as a speed, and the vertical velocity
 * there, upward positive, as the divergence of the transport above the
 * depth with the free surface held: h and u are the L layers'
 * thicknesses and edge velocities, temperature(k, i) the layer
 * temperature. Cells that are land or whose water column ends above the
 * depth are NaN, their current vector zero.
 */
export function depthFields(mesh, L, { h, u, temperature, cellOcean }, depth, { thin = 0 } = {}) {
  const { nCells: C, nEdges: E, maxEdges, nEdgesOnCell, edgesOnCell, edgeSignOnCell, cellsOnEdge, dcEdge, dvEdge, nEdge, areaCell } = mesh;
  const out = { temperature: new Float32Array(C), current: new Float32Array(3 * C), speed: new Float32Array(C), upwelling: new Float32Array(C) };
  for (let i = 0; i < C; i++) {
    let layer = -1, top = 0;
    if (cellOcean[i]) for (let k = 0; k < L; k++) { const hk = h[k * C + i]; if (k > 0 && hk <= thin) continue; if (depth < top + hk) { layer = k; break; } top += hk; }
    if (layer < 0) { out.temperature[i] = NaN; out.speed[i] = NaN; out.upwelling[i] = NaN; continue; }
    out.temperature[i] = temperature(layer, i);
    let x = 0, y = 0, z = 0, w = 0;
    for (let m = 0; m < nEdgesOnCell[i]; m++) {
      const e = edgesOnCell[maxEdges * i + m], s = edgeSignOnCell[maxEdges * i + m];
      const f = 0.5 * dcEdge[e] * dvEdge[e] * u[layer * E + e];
      x += f * nEdge[3 * e]; y += f * nEdge[3 * e + 1]; z += f * nEdge[3 * e + 2];
      const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
      let above = 0, transport = 0;
      for (let k = 0; k < L && above < depth; k++) {
        const he = 0.5 * (h[k * C + a] + h[k * C + b]);
        if (k > 0 && he <= thin) continue;
        transport += u[k * E + e] * Math.min(he, depth - above);
        above += he;
      }
      w += s * dvEdge[e] * transport;
    }
    out.current[3 * i] = x / areaCell[i]; out.current[3 * i + 1] = y / areaCell[i]; out.current[3 * i + 2] = z / areaCell[i];
    out.speed[i] = Math.hypot(x, y, z) / areaCell[i];
    out.upwelling[i] = w / areaCell[i];
  }
  return out;
}

export const EPS = 0.01, THIN = 5, PV_FLOOR = 20, SPEED_LIMIT = 5, RESTORE_TOLERANCE = 0.01;

export { FRICTION_PER_WIND, CALM_WIND, RETREAT_EKMAN, erfc, stokesSurfaceLayer, langmuirFactor, JERLOV, PENETRATION_DEPTH, transmitted, deposited };
/*
 * The ∇⁴ closure coefficient. At closureSpacing and coarser, the
 * grid-scale wave decays in closureHours; on finer meshes the
 * coefficient falls only in proportion to the spacing, so the
 * grid-scale decay time shortens as its cube (1.5 h at N=128). Jets a
 * few cells across over ridges and seamounts outgrow anything weaker at
 * N=128 (docs/c-grid-dynamical-core.md, M18).
 */
export const CLOSURE_SPACING = 120e3;
export function closureCoefficient(spacing, closureHours, closureSpacing = CLOSURE_SPACING) {
  return closureHours > 0 ? Math.pow(spacing / Math.PI, 4) / (closureHours * 3600) * Math.max(1, closureSpacing / spacing) ** 3 : 0;
}

/*
 * The longest step of the dynamics: one radian of the fastest inertial
 * rotation on the mesh. The layers carry the Coriolis force explicitly,
 * and the part of its depth sum that the barotropic solver's f·U⊥ misses
 * (a large part where the depth changes between neighboring edges) is
 * held over the step, so longer steps grow inertial oscillations. The
 * closure's grid-scale damping, about 2Ω on meshes of closureSpacing and
 * coarser, stays inside RK4's stability region at this step.
 */
export function inertialStepLimit(mesh) {
  let f = 0;
  for (let e = 0; e < mesh.nEdges; e++) f = Math.max(f, Math.abs(mesh.fEdge[e]));
  for (let v = 0; v < mesh.nVertices; v++) f = Math.max(f, Math.abs(mesh.fVertex[v]));
  return 1 / f;
}

export const EDDY_WAVE_SPEED = 2, EDDY_BOTTOM_TAPER = 100, EDDY_SLACK = 1e-6;
export function eddyDiffusivities(mesh, kappa, waveSpeed = EDDY_WAVE_SPEED) {
  const { nEdges: E, latEdge, dcEdge, radius, omega } = mesh;
  const out = new Float64Array(E);
  if (!(kappa > 0)) return out;
  for (let e = 0; e < E; e++) {
    const f = 2 * omega * Math.sin(latEdge[e]), beta = 2 * omega * Math.cos(latEdge[e]) / radius;
    const deformation = waveSpeed / Math.sqrt(f * f + 2 * beta * waveSpeed);
    out[e] = kappa / (1 + (deformation / dcEdge[e]) ** 2);
  }
  return out;
}
export function eddyDiffusionLimit(mesh) {
  const { nCells: C, maxEdges, nEdgesOnCell, edgesOnCell, dvEdge, dcEdge, areaCell } = mesh;
  let rate = 0;
  for (let i = 0; i < C; i++) {
    let sum = 0;
    for (let m = 0; m < nEdgesOnCell[i]; m++) { const e = edgesOnCell[maxEdges * i + m]; sum += dvEdge[e] / dcEdge[e]; }
    rate = Math.max(rate, sum / areaCell[i]);
  }
  return 0.25 / rate;
}

/*
 * For each edge, the cell beyond each of its two cells along the edge's
 * own direction: beyond cell a, the neighbor of a whose direction from a
 * is most nearly opposite to b's, and beyond b the neighbor of b most
 * nearly continuing from a. With the edge between them these make the
 * four-point line the upwind-biased tracer reconstruction reads; `edge`
 * holds the edge from the cell to its beyond cell, whose length is the
 * upwind spacing.
 */
export function beyondCells(mesh) {
  const { nEdges: E, cellsOnEdge, cellsOnCell, edgesOnCell, nEdgesOnCell, maxEdges, xCell } = mesh;
  const cell = new Int32Array(2 * E).fill(-1), edge = new Int32Array(2 * E).fill(-1);
  const pick = (from, exclude, dx, dy, dz, sign) => {
    let best = -Infinity, cellAt = -1, edgeAt = -1;
    for (let m = 0; m < nEdgesOnCell[from]; m++) {
      const j = cellsOnCell[maxEdges * from + m];
      if (j === exclude || j < 0) continue;
      const jx = xCell[3 * j] - xCell[3 * from], jy = xCell[3 * j + 1] - xCell[3 * from + 1], jz = xCell[3 * j + 2] - xCell[3 * from + 2];
      const score = sign * (jx * dx + jy * dy + jz * dz) / Math.hypot(jx, jy, jz);
      if (score > best) { best = score; cellAt = j; edgeAt = edgesOnCell[maxEdges * from + m]; }
    }
    return [cellAt, edgeAt];
  };
  for (let e = 0; e < E; e++) {
    const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
    const dx = xCell[3 * b] - xCell[3 * a], dy = xCell[3 * b + 1] - xCell[3 * a + 1], dz = xCell[3 * b + 2] - xCell[3 * a + 2];
    [cell[2 * e], edge[2 * e]] = pick(a, b, dx, dy, dz, -1);
    [cell[2 * e + 1], edge[2 * e + 1]] = pick(b, a, dx, dy, dz, 1);
  }
  return { cell, edge };
}

/*
 * The edge value of a tracer carried from the cell `up` to the cell
 * `down`: the upwind cell's value plus half the central difference scaled
 * by Koren's limiter on the ratio of the upwind gradient (over the spacing
 * `dUp` to the cell `uu` beyond the upwind one) to the central gradient
 * (over `dCentral`), which is the third-order upwind-biased (κ = 1/3)
 * scheme where the field is smooth and monotone everywhere. The caller
 * falls back to the upwind value where `uu` is land or either cell beyond
 * the upwind one holds no water in the layer.
 */
export const TRACER_ADVECTIONS = ['upwind', 'third'], EDGE_TRANSPORTS = ['remapped', 'kept'], SHEAR_DIFFUSIVITY = 0.2, HYCOM_DIFFUSION_SPEED = 0.005;
export function limitedEdgeValue(field, up, down, uu, dUp, dCentral) {
  const phiU = field[up], central = field[down] - phiU;
  if (central === 0) return phiU;
  const r = (phiU - field[uu]) * dCentral / dUp / central;
  return phiU + 0.5 * Math.max(0, Math.min(2 * r, (2 + r) / 3, 2)) * central;
}

/*
 * The open ocean's abyss: the sea cells at least `depth` deep that connect
 * to one another through such cells over the largest area. A deep basin
 * behind a shallower sill, as the Mediterranean behind Gibraltar or the
 * Arctic behind Fram Strait, holds water of its own density below the
 * sill and is left out.
 */
export function abyssalCells(mesh, D, cellOcean, depth) {
  const { nCells: C, maxEdges, nEdgesOnCell, cellsOnCell, areaCell } = mesh;
  const component = new Int32Array(C).fill(-1), areas = [];
  for (let start = 0; start < C; start++) {
    if (!cellOcean[start] || D[start] < depth || component[start] >= 0) continue;
    const label = areas.length, stack = [start];
    component[start] = label;
    let area = 0;
    while (stack.length) {
      const i = stack.pop();
      area += areaCell[i];
      for (let m = 0; m < nEdgesOnCell[i]; m++) {
        const j = cellsOnCell[maxEdges * i + m];
        if (cellOcean[j] && D[j] >= depth && component[j] < 0) { component[j] = label; stack.push(j); }
      }
    }
    areas.push(area);
  }
  const largest = areas.indexOf(Math.max(...areas));
  return Uint8Array.from(component, (label) => (label >= 0 && label === largest ? 1 : 0));
}

/*
 * The interior classes a saved ocean is on: the list it carries, else the
 * one of UNLISTED_LAYER_DENSITIES with its layer count, else null.
 */
export function savedDensities(saved, C) {
  if (saved && saved.densities && saved.densities.length) return Array.from(saved.densities);
  if (saved && saved.h) return UNLISTED_LAYER_DENSITIES.find((list) => saved.h.length === (list.length + 1) * C) ?? null;
  return null;
}
export function sameDensities(a, b) {
  return a.length === b.length && a.every((r, k) => Math.abs(r - b[k]) < 1e-3);
}

/*
 * The sea cell each cell's runoff reaches, following the terrain: water
 * leaves a land cell for the lowest of its neighbours that is lower
 * than it; when none is, for the lowest lower cell two steps away, then
 * three, and so on, so a pit does not hold it; the path ends at the sea.
 * Sea cells are their own outlet; -1 without a geography.
 */
export function runoffOutlets(mesh, geography) {
  const { nCells: C, maxEdges, nEdgesOnCell, cellsOnCell } = mesh;
  const outlet = new Int32Array(C).fill(-1);
  if (!geography) return outlet;
  const { land, elevation } = geography;
  const next = new Int32Array(C).fill(-1);
  const seen = new Int32Array(C).fill(-1);
  let ring = [], following = [];
  for (let i = 0; i < C; i++) {
    if (!land[i]) { next[i] = i; continue; }
    ring = [i]; seen[i] = i;
    while (ring.length && next[i] < 0) {
      following = [];
      let lowest = elevation[i], choice = -1;
      for (const c of ring) {
        for (let m = 0; m < nEdgesOnCell[c]; m++) {
          const j = cellsOnCell[maxEdges * c + m];
          if (seen[j] === i) continue;
          seen[j] = i;
          following.push(j);
          if (elevation[j] < lowest) { lowest = elevation[j]; choice = j; }
        }
      }
      if (choice >= 0) next[i] = choice;
      ring = following;
    }
  }
  for (let i = 0; i < C; i++) {
    let c = i;
    while (land[c] && next[c] >= 0 && next[c] !== c) c = next[c];
    outlet[i] = land[c] ? -1 : c;
  }
  return outlet;
}

/*
 * The model's bathymetry: the cell-mean ETOPO depth of every sea cell,
 * at least `minimumDepth`, and never shallower than `neighbourRatio`
 * times its deepest sea neighbour, so that no shelf break or trench
 * wall drops by more than that ratio across one edge. A cell 110 km
 * across cannot hold a real shelf, and a ten-to-one step in depth
 * between two cells makes the free surface swing wildly at the coast.
 */
export function bathymetryFrom(mesh, geography, { minimumDepth = 50, neighbourRatio = 0.5, flatDepth = 4000 } = {}) {
  const { nCells: C, maxEdges, nEdgesOnCell, cellsOnCell } = mesh;
  const D = new Float64Array(C);
  for (let i = 0; i < C; i++) D[i] = geography ? (geography.land[i] ? 0 : Math.max(minimumDepth, -geography.elevation[i])) : flatDepth;
  if (!geography) return D;
  for (let sweep = 0; sweep < 50; sweep++) {
    let changed = false;
    for (let i = 0; i < C; i++) {
      if (geography.land[i]) continue;
      let deepest = 0;
      for (let m = 0; m < nEdgesOnCell[i]; m++) { const j = cellsOnCell[maxEdges * i + m]; if (!geography.land[j]) deepest = Math.max(deepest, D[j]); }
      const floor = neighbourRatio * deepest;
      if (D[i] < floor) { D[i] = floor; changed = true; }
    }
    if (!changed) break;
  }
  return D;
}

/*
 * The depth each edge's column is clipped to: the shallower cell, and
 * the geography's edge sill where it is shallower still,
 * no shallower than `floor`, so raster land parting two sea cells narrows
 * the edge rather than closing it.
 */
export function edgeSills(mesh, D, edgeOcean, edgeSill, floor) {
  const { nEdges: E, cellsOnEdge } = mesh, out = new Float64Array(E);
  for (let e = 0; e < E; e++) {
    const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
    out[e] = Math.min(D[a], D[b]);
    if (edgeSill && edgeOcean[e]) out[e] = Math.min(out[e], Math.max(floor, edgeSill[e]));
  }
  return out;
}

/*
 * The ocean's cells: wet where the water's geography (`ocean` of a
 * geography with ice-shelf cavities, ./geography.module.js) is sea, and
 * open to the air only where the atmosphere's is; `shelf` marks the
 * cavities between, and `draft` the depth of their shelves' base.
 */
export function oceanCells(C, geography) {
  const water = geography ? geography.ocean ?? geography : null;
  const cellOcean = water ? Uint8Array.from(water.land, (l) => (l ? 0 : 1)) : new Uint8Array(C).fill(1);
  const shelf = geography && geography.cavity ? Uint8Array.from(geography.cavity) : new Uint8Array(C);
  const cellOpen = Uint8Array.from(cellOcean, (wet, i) => (wet && !shelf[i] ? 1 : 0));
  const draft = Float64Array.from(shelf, (s, i) => (s ? geography.iceDraft[i] : 0));
  return { water, cellOcean, shelf, cellOpen, draft };
}

/*
 * A column under a shelf rests between the shelf's draft and the bed: its
 * resting thickness D is the bed's depth less the draft, and an edge
 * passes water from the deeper of its two drafts down to its sill, or
 * through the sill floor where that is less: across a cell the shelf's
 * base slopes, and the water beneath it follows.
 */
export function restingColumns(mesh, bed, draft, edgeOcean, edgeSill, floor) {
  const D = Float64Array.from(bed, (d, i) => d - draft[i]), sillEdge = edgeSills(mesh, bed, edgeOcean, edgeSill, floor);
  for (let e = 0; e < mesh.nEdges; e++) { const top = Math.max(draft[mesh.cellsOnEdge[2 * e]], draft[mesh.cellsOnEdge[2 * e + 1]]); if (top > 0) sillEdge[e] = Math.max(edgeOcean[e] ? floor : 0, sillEdge[e] - top); }
  return { D, sillEdge };
}

/*
 * The weight per unit area (kg/m²) of the water from the surface down to a
 * shelf's draft in the profile sample(z) → [T, S] under the equation of
 * state density(T, S): the shelf's load when it floats in that water.
 */
export function displacedWeight(draft, sample, density = seawaterDensity, n = 64) {
  let sum = 0;
  for (let m = 0; m < n; m++) { const [t, s] = sample((m + 0.5) * draft / n); sum += density(t, s); }
  return sum * draft / n;
}

export function createOcean(mesh, {
  minimumDepth = 50, flatDepth = 4000, thermoclineTilt = 0.3,
  salinityProfile = (lat) => 34 + 2 * Math.exp(-(((Math.abs(lat) * 180 / Math.PI - 25) / 20) ** 2)),
  density = 1025, specificHeat = 3985, referenceS = 35, gravity = 9.81,
  vorticityCentring = 0.5, iceStressTransmission = 0.8, iceSalinity = 5, iceDensity = 917, latentHeatFusion = 3.34e5, langmuir = 2, langmuirIce = 'none',
  momentumDepth = 50, shearViscosity = 1e-2, backgroundViscosity = 1e-4, bottomDrag = 3e-3, topDrag = 2.5e-3, closureHours = 12, closureSpacing = CLOSURE_SPACING, everySteps = 4,
  eddyDiffusivity = 1000, eddyTaperDepth = 200, penetration = 'IB',
  geography = null, bathymetry = null, edgeSill = undefined, sillFloor = minimumDepth, buffers = null, climatology = null,
  vertical = 'hybrid', generator = 'hybrid', layerTable = HYBRID_LAYERS, regridTime = 0, minimumLayer = H_MIN, modeDepth = MODE_DEPTH, mstar = 1.25, ekmanFactor = 0.5, rotationLimit = 'algebraic', nstar = 0.2, backgroundDiffusivity = 1e-5, convectiveDiffusivity = 0.1,
  tracerAdvection = 'third', shearDiffusivity = SHEAR_DIFFUSIVITY, shearDiffusivityCap = 4000, lateralDiffusionSpeed = HYCOM_DIFFUSION_SPEED,
  edgeTransport = 'kept', shelfMelt = true, tidalSpeed = SHELF_MELT.tidalSpeed, shelfHeatTransfer = SHELF_MELT.heatTransfer, shelfSaltTransfer = SHELF_MELT.saltTransfer, shelfMeltLayer = 0,
} = {}) {
  const {
    nCells: C, nEdges: E, nVertices: V, maxEdgesOnEdge, nEdgesOnEdge, edgesOnEdge, weightsOnEdge,
    cellsOnEdge, verticesOnEdge, areaCell, dcEdge, dvEdge, fVertex, fEdge, fCell, radius, latCell,
    maxEdges, nEdgesOnCell, edgesOnCell, edgeSignOnCell, cellsOnCell, nEdge,
  } = mesh;
  if (climatology && typeof climatology.columnAt !== 'function') throw new Error('the ocean climatology must be decoded (loadClimatology in ./climatology.module.js) before the ocean is built');
  if (penetration !== null && !JERLOV[penetration]) throw new Error(`penetration is one of ${Object.keys(JERLOV).join(', ')} or null, not ${penetration}`);
  if (vertical !== 'hybrid') throw new Error(`the ocean's vertical coordinate is 'hybrid'; a saved slab ocean is converted on load, not run (vertical ${vertical})`);
  if (langmuir !== 0 && !(langmuir >= 1)) throw new Error(`langmuir caps a factor of at least 1, so it is 0 or at least 1, not ${langmuir}`);
  if (!ROTATION_LIMITS.includes(rotationLimit)) throw new Error(`rotationLimit is one of ${ROTATION_LIMITS.join(', ')}, not ${rotationLimit}`);
  if (!LANGMUIR_ICE.includes(langmuirIce)) throw new Error(`langmuirIce is one of ${LANGMUIR_ICE.join(', ')}, not ${langmuirIce}`);
  const advection = tracerAdvection;
  if (!TRACER_ADVECTIONS.includes(advection)) throw new Error(`tracerAdvection is one of ${TRACER_ADVECTIONS.join(', ')}, not ${advection}`);
  if (!EDGE_TRANSPORTS.includes(edgeTransport)) throw new Error(`edgeTransport is one of ${EDGE_TRANSPORTS.join(', ')}, not ${edgeTransport}`);
  const beyond = advection === 'third' ? beyondCells(mesh) : null;
  const shearDiffusion = shearDiffusivity;
  if (!(shearDiffusion >= 0) || !(shearDiffusivityCap > 0)) throw new Error(`shearDiffusivity is a non-negative coefficient and shearDiffusivityCap a positive diffusivity, not ${shearDiffusion} and ${shearDiffusivityCap}`);
  if (advection === 'third' && !(shearDiffusion > 0)) throw new Error('the limited third-order tracer flux admits grid-scale noise that the shear-scaled diffusion of the surface layers damps, so tracerAdvection third needs shearDiffusivity above zero');
  const lateralSpeed = lateralDiffusionSpeed;
  if (!(lateralSpeed >= 0)) throw new Error(`lateralDiffusionSpeed is a non-negative velocity scale, not ${lateralSpeed}`);
  const table = layerTable;
  const L = table.K, g = gravity, rho0 = density, rhoCp = density * specificHeat;
  const rho = Array.from(table.densities);
  const labels = Array.from(table.densities);
  const labelS = rho.map((r) => { const j = LAYER_DENSITIES.indexOf(r); return j >= 0 ? LAYER_SALINITIES[j] : referenceS; });
  const labelT = rho.map((r, k) => Math.max(FREEZING_POINT, labelTemperature(r, labelS[k])));
  let spacing = 0, minSpacing = Infinity;
  for (let e = 0; e < E; e++) { spacing += dcEdge[e]; minSpacing = Math.min(minSpacing, dcEdge[e]); }
  spacing /= E;
  const nu4 = closureCoefficient(spacing, closureHours, closureSpacing);
  const eddyKappa = eddyDiffusivities(mesh, eddyDiffusivity), eddyLimit = eddyDiffusivity > 0 ? eddyDiffusionLimit(mesh) : 0;
  let eddy = null;

  const { water, cellOcean, shelf, cellOpen, draft } = oceanCells(C, geography);
  const edgeOcean = water ? water.edgeOcean : new Uint8Array(E).fill(1);
  const outlet = runoffOutlets(mesh, geography);
  const bed = bathymetry ? Float64Array.from(bathymetry, (d, i) => (cellOcean[i] ? d : 0)) : bathymetryFrom(mesh, water, { minimumDepth, flatDepth });
  const { D, sillEdge } = restingColumns(mesh, bed, draft, edgeOcean, edgeSill === undefined ? (water ? water.edgeSill ?? null : null) : edgeSill, sillFloor);
  const shelfEdge = Float64Array.from({ length: E }, (_, e) => 0.5 * (shelf[mesh.cellsOnEdge[2 * e]] + shelf[mesh.cellsOnEdge[2 * e + 1]])), hasShelves = shelf.some((s) => s);
  let deepest = 0;
  for (let i = 0; i < C; i++) deepest = Math.max(deepest, D[i]);
  const substepLimit = 0.35 * minSpacing / Math.sqrt(g * Math.max(deepest, 1));
  const inertialLimit = inertialStepLimit(mesh);

  const adopting = !!(buffers && buffers.h);
  const shared = (name, n) => new Float64Array(buffers && buffers[name] ? buffers[name] : new SharedArrayBuffer(8 * n));
  const h = shared('h', L * C), u = shared('u', L * E), Q = shared('Q', L * C), W = shared('W', L * C);
  const state = [h, u, Q, W];
  const eta = shared('eta', C);
  const T0 = new Float64Array(C), S0 = new Float64Array(C), rhoMl = shared('rhoMl', C), previousT0 = new Float64Array(C), surfaceIn = new Float64Array(C), iceVolume = new Float64Array(C), previousVolume = new Float64Array(C);
  let seaCover = null;
  const capacity = shared('capacity', C);
  if (!adopting) capacity.fill(rhoCp * table.depths[1]);
  const transmission = shared('transmission', C), penetrated = shared('penetrated', C), shares = new Float64Array(L);
  const waves = shared('waves', C);
  if (!(buffers && buffers.waves)) waves.fill(1);
  const stress = shared('stress', E), fresh = new Float64Array(C), calved = new Float64Array(C);
  const iced = new Uint8Array(C);
  const hEdge = shared('hEdge', L * E), pressure = shared('pressure', L * C), flux = new Float64Array(E), fluxPV = new Float64Array(E), tracerFlux = new Float64Array(E), edgeS = new Float64Array(E);
  const T = new Float64Array(C), S = new Float64Array(C);
  const zeta = new Float64Array(V), qEdge = new Float64Array(E);
  const K = new Float64Array(C), gradPhi = new Float64Array(E), gradEta = shared('gradEta', E);
  const lap = new Float64Array(E), lap2 = new Float64Array(E), divScratch = new Float64Array(C), curlScratch = new Float64Array(V);
  const slow = shared('baroSlow', E), U = shared('baroU', E), depthEdge = shared('baroDepth', E), etaB = shared('baroEta', C), avgU = shared('baroAvgU', E), avgEta = shared('baroAvgEta', C);
  const countBuffer = buffers && buffers.counts ? buffers.counts : new SharedArrayBuffer(16), counts = new Int32Array(countBuffer);
  const stages = [0, 1, 2, 3].map((s) => [shared(`stage${s}h`, L * C), shared(`stage${s}u`, L * E), shared(`stage${s}Q`, L * C), shared(`stage${s}W`, L * C)]);
  const trial = [shared('trialh', L * C), shared('trialu', L * E), shared('trialQ', L * C), shared('trialW', L * C)];
  const params = shared('params', 8);
  if (!adopting) { params[0] = 1 / 3600; params[1] = 3600; }
  const buoyancyLoss = new Float64Array(C);
  let counter = 0, limited = 0, initialised = false, layerRunner = null, columnRunner = null, dried = 0;
  const relaxRate = () => params[0];
  const shelfLoad = shared('shelfLoad', C);
  const melted = shared('melted', C), meltBase = shared('meltBase', C), shelfMelted = shared('shelfMelted', C), shelfDraw = shared('shelfDraw', C);
  const meltOn = shelfMelt && shelf.some((s) => s), meltPool = (i) => (fCell[i] < 0 ? 0 : 1), meltConstants = { rho0, rhoCp, latentHeat: latentHeatFusion, heatTransfer: shelfHeatTransfer, saltTransfer: shelfSaltTransfer };
  const shelfBudget = { discharge: [0, 0], melt: [0, 0], moved: [0, 0], drawn: [0, 0], keep: [1, 1], share: [0, 0] };
  const layerRho = shared('layerRho', L * C), columnAbove = shared('columnAbove', L * C), targetZ = shared('targetZ', (L + 1) * C), viscosity = shared('viscosity', L * C), edgeBottom = shared('edgeBottom', E);
  const mixedDepths = shared('mixedDepths', C), boundaryDepths = shared('boundaryDepths', C);
  const edgeAbove = shared('edgeAbove', L * E), eddyFlux = shared('eddyFlux', L * E), eddyVolume = shared('eddyVolume', L * C), eddyHeat = shared('eddyHeat', L * C), eddySalt = shared('eddySalt', L * C);
  const pvFloor = Float64Array.from(table.minimum, (m) => Math.min(PV_FLOOR, m)), tokenEdge = Float64Array.from(table.minimum, (m) => 2 * Math.min(m, table.spacing ?? m)), pinned = new Uint8Array(E);
  const boundaries = labelBoundaries(rho), cs = columnScratch(L);
  const cols = Object.fromEntries(['h', 'T', 'S', 'rho', 'shear', 'kappa', 'nu', 'hNew', 'amount', 'hOld', 'hEdgeNew', 'uCol', 'nuEdge'].map((name) => [name, new Float64Array(L)]));
  cols.z = new Float64Array(L + 1); cols.w = new Float64Array(L + 1); cols.vector = new Float64Array(3 * L);
  const interfaceWeight = (k, i) => (generator === 'isopycnal' ? 1 : generator === 'zstar' ? 0 : Math.min(1, Math.max(0, (h[at(k - 1, i)] - table.minimum[k - 1]) / modeDepth)));

  const eos = seawaterDensity;
  const at = (k, i) => k * C + i;
  const ae = (k, e) => k * E + e;

  /*
   * The top layer, present everywhere, takes the centred thickness at an
   * edge; a layer below takes the smaller of the two, so it never flows
   * into a cell where it has no water, whether the layer has vanished
   * there or the bottom lies above it; and the column is cut at the
   * edge's sill, the layers below it carrying nothing across.
   */
  function edgeColumn(hIn, e, out, base, stride) {
    const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
    out[base] = 0.5 * (hIn[a] + hIn[b]);
    for (let k = 1; k < L; k++) out[base + k * stride] = Math.min(hIn[at(k, a)], hIn[at(k, b)]);
    const sill = Math.max(EPS, sillEdge[e] + 0.5 * (eta[a] + eta[b]));
    let above = 0;
    for (let k = 0; k < L; k++) { const hk = out[base + k * stride]; out[base + k * stride] = Math.min(hk, Math.max(0, sill - above)); above += hk; }
  }
  let flowDepth = 0;
  function edgeFlow(thickness, e) {
    const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
    const sill = Math.max(EPS, sillEdge[e] + 0.5 * (eta[a] + eta[b]));
    let above = 0, transport = 0, depth = 0;
    for (let k = 0; k < L; k++) {
      const hk = k === 0 ? 0.5 * (thickness(0, a) + thickness(0, b)) : Math.min(thickness(k, a), thickness(k, b)), he = Math.min(hk, Math.max(0, sill - above));
      above += hk; transport += he * u[ae(k, e)]; depth += he;
    }
    flowDepth = depth;
    return transport;
  }
  function keepFlow(thickness, e, before) {
    const after = edgeFlow(thickness, e);
    if (!(flowDepth > EPS)) return;
    const shift = (before - after) / flowDepth;
    for (let k = 0; k < L; k++) u[ae(k, e)] += shift;
  }
  const stateThickness = (k, i) => h[at(k, i)], targetThickness = (k, i) => targetZ[(k + 1) * C + i] - targetZ[k * C + i];
  function freeSurface() {
    for (let i = 0; i < C; i++) {
      let sum = 0;
      for (let k = 0; k < L; k++) sum += h[at(k, i)];
      eta[i] = cellOcean[i] ? sum - D[i] : 0;
    }
  }

  /*
   * The fields every layer's tendency reads: each cell-layer's density,
   * the running sums of ρh and h above it, the free-surface gradient and
   * the edge thicknesses with their running sums down the edge column.
   */
  function mixedDepthsAll() {
    const { h: hCol, rho: rhoCol } = cols;
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) { mixedDepths[i] = 0; continue; }
      for (let k = 0; k < L; k++) { hCol[k] = h[at(k, i)]; rhoCol[k] = layerRho[at(k, i)]; }
      mixedDepths[i] = mixedLayerDepth(hCol, rhoCol, L);
    }
  }
  function layerDensitiesAll(hIn = h, QIn = Q, WIn = W) {
    for (let k = 0; k < L; k++) for (let i = 0; i < C; i++) { const n = at(k, i), hv = hIn[n]; layerRho[n] = hv > 1e-6 ? eos(QIn[n] / hv, WIn[n] / hv) : rho[k]; }
    for (let i = 0; i < C; i++) rhoMl[i] = layerRho[i];
  }
  function prepareCells(from, to, input) {
    const [hIn, , QIn, WIn] = input;
    for (let i = from; i < to; i++) {
      let weighted = shelfLoad[i], total = draft[i];
      for (let k = 0; k < L; k++) {
        const n = at(k, i), hv = hIn[n];
        layerRho[n] = hv > 1e-6 ? eos(QIn[n] / hv, WIn[n] / hv) : rho[k];
        pressure[n] = weighted; columnAbove[n] = total;
        weighted += layerRho[n] * hIn[n]; total += hIn[n];
      }
      rhoMl[i] = layerRho[i];
    }
  }
  function prepareEdges(from, to, input) {
    const [hIn] = input;
    for (let e = from; e < to; e++) {
      gradEta[e] = (eta[cellsOnEdge[2 * e + 1]] - eta[cellsOnEdge[2 * e]]) / dcEdge[e];
      edgeColumn(hIn, e, hEdge, e, E);
      let sum = 0;
      for (let k = 0; k < L; k++) { edgeAbove[ae(k, e)] = sum; sum += hEdge[ae(k, e)]; }
      edgeBottom[e] = sum;
    }
  }
  function prepare(input) {
    if (columnRunner) { const useTrial = input === trial ? 1 : 0; columnRunner('prepareCells', useTrial); columnRunner('prepareEdges', useTrial); return; }
    prepareCells(0, C, input);
    prepareEdges(0, E, input);
  }

  function tendency(input, out) {
    prepare(input);
    if (layerRunner) layerRunner(input === trial ? 1 : 0, stages.indexOf(out));
    else tendencyLayers(input, out, 0, L);
  }

  /*
   * The layers' tendencies, kFrom to kTo, every layer alike, reading the
   * prepared fields and the input state alone so that layers can be
   * computed by any thread in any order; `part` 'tracers' writes only
   * dh, dQ and dW and 'momentum' only du. The flux through an
   * edge is the edge's mean thickness, no more than the donor cell's, times
   * the normal velocity; the potential vorticity's thickness is no less than
   * vorticityCentring times the cells' mean nor than the layer's minimum
   * thickness (PV_FLOOR at most), so the fixed-depth surface layers feel
   * their full Coriolis force; the pressure force is the hydrostatic
   * gradient at fixed depth, −g∇η − (g/ρ₀)Δ[P_k + ρ_k (z − H_k)]/dc with
   * P_k = Σ_{j<k} ρ_j h_j and H_k = Σ_{j<k} h_j the sums above the layer,
   * each cell's layer line taken at the middle z of the depths the layer
   * occupies in both cells, or of the gap between them where it occupies
   * none in common, which is the Montgomery gradient for a class at its
   * label and the layer-mean term ½ h Δρ for a layer at the same depths
   * on both sides. The edge thickness of a layer is the smaller of its
   * two cells', so a layer that has vanished on either side carries next
   * to nothing, and a layer no thicker than twice minimumLayer at an edge,
   * or one that occupies no depth in common between its two cells, each
   * measured down from its column's top, while
   * no thicker than tokenEdge (twice its minimum, at most twice the
   * thermocline spacing: a class held at its floor on one side against
   * water of that class elsewhere in the other column), follows the
   * layer above, relaxing to its velocity here and copied from it after
   * the barotropic shift; the stress is a flux into the top layer; the bottom drag acts on the part
   * of each layer within momentumDepth of the floor, both measured down
   * the edge column's own thicknesses (edgeAbove, edgeBottom); the ∇⁴
   * closure on every layer.
   */
  function tendencyLayers(input, out, kFrom, kTo, part = 'all') {
    const [hIn, uIn, QIn, WIn] = input;
    const [dh, du, dQ, dW] = out;
    const relax = relaxRate(), tracers = part !== 'momentum', momentum = part !== 'tracers', thinEdge = 2 * minimumLayer;
    for (let k = kFrom; k < kTo; k++) {
      const oc = k * C, oe = k * E;
      for (let e = 0; e < E; e++) {
        const he = hEdge[oe + e], ue = uIn[oe + e];
        fluxPV[e] = edgeOcean[e] ? he * ue : 0;
        flux[e] = edgeOcean[e] ? Math.min(he, Math.max(0, hIn[oc + cellsOnEdge[2 * e + (ue > 0 ? 0 : 1)]])) * ue : 0;
      }
      if (tracers) {
        divergence(mesh, flux, divScratch);
        for (let i = 0; i < C; i++) {
          dh[oc + i] = -divScratch[i];
          const hv = hIn[oc + i];
          if (hv > 1e-6) { T[i] = QIn[oc + i] / hv; S[i] = WIn[oc + i] / hv; } else { T[i] = labelT[k]; S[i] = labelS[k]; }
        }
        if (beyond) {
          for (let e = 0; e < E; e++) {
            const f = flux[e], side = f > 0 ? 0 : 1, up = cellsOnEdge[2 * e + side], down = cellsOnEdge[2 * e + 1 - side], uu = beyond.cell[2 * e + side];
            if (f === 0 || !cellOcean[uu] || hIn[oc + uu] <= 1e-6 || hIn[oc + down] <= 1e-6) { tracerFlux[e] = f * T[up]; edgeS[e] = S[up]; continue; }
            const dUp = dcEdge[beyond.edge[2 * e + side]], dCentral = dcEdge[e];
            tracerFlux[e] = f * limitedEdgeValue(T, up, down, uu, dUp, dCentral);
            edgeS[e] = limitedEdgeValue(S, up, down, uu, dUp, dCentral);
          }
          divergence(mesh, tracerFlux, divScratch);
          for (let i = 0; i < C; i++) dQ[oc + i] = -divScratch[i];
          for (let e = 0; e < E; e++) tracerFlux[e] = flux[e] * edgeS[e];
        } else {
          for (let e = 0; e < E; e++) tracerFlux[e] = flux[e] * T[cellsOnEdge[2 * e + (flux[e] > 0 ? 0 : 1)]];
          divergence(mesh, tracerFlux, divScratch);
          for (let i = 0; i < C; i++) dQ[oc + i] = -divScratch[i];
          for (let e = 0; e < E; e++) tracerFlux[e] = flux[e] * S[cellsOnEdge[2 * e + (flux[e] > 0 ? 0 : 1)]];
        }
        divergence(mesh, tracerFlux, divScratch);
        for (let i = 0; i < C; i++) dW[oc + i] = -divScratch[i];
        if ((shearDiffusion > 0 && k < table.fixed) || lateralSpeed > 0) {
          const sheared = shearDiffusion > 0 && k < table.fixed;
          if (sheared) curl(mesh, uIn.subarray(oe, oe + E), zeta);
          for (let e = 0; e < E; e++) {
            const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1], dc = dcEdge[e];
            const kappa = lateralSpeed * dc + (sheared ? Math.min(shearDiffusivityCap, shearDiffusion * Math.abs(0.5 * (zeta[verticesOnEdge[2 * e]] + zeta[verticesOnEdge[2 * e + 1]])) * dc * dc) : 0);
            const weight = edgeOcean[e] ? kappa * hEdge[oe + e] / dc : 0;
            tracerFlux[e] = -weight * (T[b] - T[a]); edgeS[e] = -weight * (S[b] - S[a]);
          }
          divergence(mesh, tracerFlux, divScratch);
          for (let i = 0; i < C; i++) dQ[oc + i] -= divScratch[i];
          divergence(mesh, edgeS, divScratch);
          for (let i = 0; i < C; i++) dW[oc + i] -= divScratch[i];
        }
        for (let i = 0; i < C; i++) if (!cellOcean[i]) { dh[oc + i] = 0; dQ[oc + i] = 0; dW[oc + i] = 0; }
      }
      if (!momentum) continue;
      curl(mesh, uIn.subarray(oe, oe + E), zeta);
      for (let e = 0; e < E; e++) {
        const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
        qEdge[e] = 0.5 * (zeta[verticesOnEdge[2 * e]] + fVertex[verticesOnEdge[2 * e]] + zeta[verticesOnEdge[2 * e + 1]] + fVertex[verticesOnEdge[2 * e + 1]]) / Math.max(hEdge[oe + e], vorticityCentring * 0.5 * (hIn[oc + a] + hIn[oc + b]), pvFloor[k]);
      }
      kineticEnergy(mesh, uIn.subarray(oe, oe + E), K);
      gradient(mesh, K, gradPhi);
      for (let e = 0; e < E; e++) {
        let sum = 0;
        for (let s = 0; s < nEdgesOnEdge[e]; s++) {
          const other = edgesOnEdge[maxEdgesOnEdge * e + s];
          sum += weightsOnEdge[maxEdgesOnEdge * e + s] * dvEdge[other] * fluxPV[other] * 0.5 * (qEdge[e] + qEdge[other]);
        }
        const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1], na = oc + a, nb = oc + b;
        const Ha = columnAbove[na], Hb = columnAbove[nb], top = Math.max(Ha, Hb), bottom = Math.min(Ha + hIn[na], Hb + hIn[nb]);
        const zR = bottom > top ? 0.5 * (top + bottom) : 0.5 * (Ha + 0.5 * hIn[na] + Hb + 0.5 * hIn[nb]);
        const layerForce = pressure[nb] + layerRho[nb] * (zR - Hb) - pressure[na] - layerRho[na] * (zR - Ha);
        du[oe + e] = sum / dcEdge[e] - gradPhi[e] - g * gradEta[e] - g / rho0 * layerForce / dcEdge[e];
        const he = hEdge[oe + e];
        const ra = Ha - draft[a], rb = Hb - draft[b];
        pinned[e] = he <= thinEdge || (he <= tokenEdge[k] && !(Math.min(ra + hIn[na], rb + hIn[nb]) > Math.max(ra, rb))) ? 1 : 0;
        if (k === 0) du[e] += stress[e] / rho0 / Math.max(he, minimumLayer) - topDrag * shelfEdge[e] * Math.abs(uIn[e]) * uIn[e] / Math.max(he, minimumLayer);
        const layerTop = edgeAbove[oe + e], floor = edgeBottom[e];
        const overlap = Math.min(layerTop + he, floor) - Math.max(layerTop, floor - momentumDepth);
        if (overlap > 0) du[oe + e] -= bottomDrag * Math.abs(uIn[oe + e]) * uIn[oe + e] * (overlap / momentumDepth) / Math.max(he, minimumLayer);
      }
      if (nu4 > 0) {
        laplacianVelocity(mesh, uIn.subarray(oe, oe + E), lap, divScratch, curlScratch);
        laplacianVelocity(mesh, lap, lap2, divScratch, curlScratch);
        for (let e = 0; e < E; e++) du[oe + e] -= nu4 * lap2[e];
      }
      if (k > 0) for (let e = 0; e < E; e++) if (pinned[e]) du[oe + e] = (uIn[ae(k - 1, e)] - uIn[oe + e]) * relax;
      for (let e = 0; e < E; e++) if (!edgeOcean[e]) du[oe + e] = 0;
    }
  }

  function combine(target, base, stage, dt) {
    for (let a = 0; a < state.length; a++) {
      const t = target[a], b = base[a], s = stage[a];
      for (let i = 0; i < t.length; i++) t[i] = b[i] + dt * s[i];
    }
  }

  function baroPrep(from, to) {
    const [, du] = stages[0];
    for (let e = from; e < to; e++) {
      let sumH = 0, transport = 0, forcing = 0;
      for (let k = 0; k < L; k++) { const he = hEdge[ae(k, e)]; sumH += he; transport += he * u[ae(k, e)]; forcing += he * (du[ae(k, e)] + g * gradEta[e]); }
      U[e] = edgeOcean[e] ? transport : 0;
      slow[e] = edgeOcean[e] ? forcing : 0;
      depthEdge[e] = sumH;
    }
  }
  function edgeCoriolis(Ue, e) {
    let sum = 0;
    for (let s = 0; s < nEdgesOnEdge[e]; s++) {
      const other = edgesOnEdge[maxEdgesOnEdge * e + s];
      sum += weightsOnEdge[maxEdgesOnEdge * e + s] * dvEdge[other] * Ue[other] * 0.5 * (fEdge[e] + fEdge[other]);
    }
    return sum / dcEdge[e];
  }
  function baroCoriolis(from, to) {
    for (let e = from; e < to; e++) slow[e] -= edgeCoriolis(U, e);
  }
  const STAGE_FACTOR = [0.5, 0.5, 1, 0];
  const baroInput = (s) => (s === 0 ? [etaB, U] : s === 2 ? bTrialB : bTrial);
  const baroOutput = (s) => (s === 1 ? bTrialB : bTrial);
  function baroEdges(from, to) {
    const s = params[4], dtb = params[5], weight = params[6], [etaIn, uIn] = baroInput(s), dU = bStages[s][1], out = baroOutput(s)[1];
    for (let e = from; e < to; e++) {
      const d = edgeOcean[e] ? -g * depthEdge[e] * (etaIn[cellsOnEdge[2 * e + 1]] - etaIn[cellsOnEdge[2 * e]]) / dcEdge[e] + edgeCoriolis(uIn, e) + slow[e] : 0;
      dU[e] = d;
      if (s < 3) out[e] = U[e] + STAGE_FACTOR[s] * dtb * d;
      else { U[e] += dtb / 6 * (bStages[0][1][e] + 2 * bStages[1][1][e] + 2 * bStages[2][1][e] + d); avgU[e] += U[e] * weight; }
    }
  }
  function baroCells(from, to) {
    const s = params[4], dtb = params[5], weight = params[6], [, uIn] = baroInput(s), dEta = bStages[s][0], out = baroOutput(s)[0];
    for (let i = from; i < to; i++) {
      let sum = 0;
      for (let m = 0; m < nEdgesOnCell[i]; m++) { const e = edgesOnCell[maxEdges * i + m]; sum += edgeSignOnCell[maxEdges * i + m] * uIn[e] * dvEdge[e]; }
      const d = cellOcean[i] ? -(sum / areaCell[i]) : 0;
      dEta[i] = d;
      if (s < 3) out[i] = etaB[i] + STAGE_FACTOR[s] * dtb * d;
      else { etaB[i] += dtb / 6 * (bStages[0][0][i] + 2 * bStages[1][0][i] + 2 * bStages[2][0][i] + d); avgEta[i] += etaB[i] * weight; }
    }
  }
  function barotropic(dt) {
    const M = Math.max(1, Math.ceil(dt / substepLimit)), dtb = dt / M;
    params[5] = dtb; params[6] = 1 / M;
    const pass = (name, n, fn) => (columnRunner ? columnRunner(name) : fn(0, n));
    pass('baroPrep', E, baroPrep);
    pass('baroCoriolis', E, baroCoriolis);
    etaB.set(eta); avgU.fill(0); avgEta.fill(0);
    for (let m = 0; m < M; m++) for (let stage = 0; stage < 4; stage++) {
      params[4] = stage;
      pass('baroEdges', E, baroEdges);
      pass('baroCells', C, baroCells);
    }
  }
  const bStages = [0, 1, 2, 3].map((s) => [shared(`baroStage${s}Eta`, C), shared(`baroStage${s}U`, E)]);
  const bTrial = [shared('baroTrialEtaA', C), shared('baroTrialUA', E)], bTrialB = [shared('baroTrialEtaB', C), shared('baroTrialUB', E)];
  function dynamics(dt) {
    const [k1, k2, k3, k4] = stages;
    params[0] = Math.min(1 / 3600, 1 / dt);
    params[1] = dt;
    freeSurface();
    tendency(state, k1);
    barotropic(dt);
    combine(trial, state, k1, dt / 2);
    tendency(trial, k2);
    combine(trial, state, k2, dt / 2);
    tendency(trial, k3);
    combine(trial, state, k3, dt);
    tendency(trial, k4);
    const w = dt / 6;
    for (let a = 0; a < state.length; a++) {
      const s = state[a], a1 = k1[a], a2 = k2[a], a3 = k3[a], a4 = k4[a];
      for (let i = 0; i < s.length; i++) s[i] += w * (a1[i] + 2 * a2[i] + 2 * a3[i] + a4[i]);
    }
  }
  function rescaleCells(from, to) {
    for (let i = from; i < to; i++) {
      if (!cellOcean[i]) continue;
      let sum = 0;
      for (let k = 0; k < L; k++) { const n = at(k, i); if (!(h[n] > 0)) { if (h[n] < 0) Atomics.add(counts, 0, 1); h[n] = 0; Q[n] = 0; W[n] = 0; } sum += h[n]; }
      const scale = (D[i] + avgEta[i]) / sum;
      for (let k = 0; k < L; k++) { const n = at(k, i); h[n] *= scale; Q[n] *= scale; W[n] *= scale; }
      eta[i] = avgEta[i];
    }
  }
  function rescale() {
    if (columnRunner) columnRunner('rescale'); else rescaleCells(0, C);
    dried = Atomics.load(counts, 0);
  }
  function shiftEdges(from, to) {
    for (let e = from; e < to; e++) {
      edgeColumn(h, e, hEdge, e, E);
      if (!edgeOcean[e]) continue;
      let sumH = 0, transport = 0, clamped = false;
      const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
      let Ha = h[a], Hb = h[b];
      for (let k = 1; k < L; k++) {
        const ha = h[at(k, a)], hb = h[at(k, b)], he = hEdge[ae(k, e)];
        if (he < 2 * minimumLayer || (he <= tokenEdge[k] && Math.min(Ha + ha, Hb + hb) <= Math.max(Ha, Hb))) u[ae(k, e)] = u[ae(k - 1, e)];
        Ha += ha; Hb += hb;
      }
      for (let k = 0; k < L; k++) { sumH += hEdge[ae(k, e)]; transport += hEdge[ae(k, e)] * u[ae(k, e)]; }
      const shift = (avgU[e] - transport) / Math.max(sumH, EPS);
      for (let k = 0; k < L; k++) { const n = ae(k, e); u[n] += shift; if (Math.abs(u[n]) > SPEED_LIMIT) { u[n] = Math.sign(u[n]) * SPEED_LIMIT; clamped = true; } }
      if (clamped) Atomics.add(counts, 1, 1);
    }
  }
  function shiftVelocities() {
    Atomics.store(counts, 1, 0);
    if (columnRunner) columnRunner('shift'); else shiftEdges(0, E);
    limited = Atomics.load(counts, 1);
  }
  function step(dt) {
    const parts = Math.max(1, Math.ceil(dt / inertialLimit)), part = dt / parts;
    for (let m = 0; m < parts; m++) {
      dynamics(part);
      rescale();
      if (eddyLimit > 0) eddyTransport(part);
      shiftVelocities();
    }
  }
  const phases = { dynamics, rescale, eddy: eddyTransport, shift: shiftVelocities, regrid, column: columnPhysics };

  function eddyEdges(from, to) {
    const dt = params[1];
    if (!eddy) eddy = { flux: new Float64Array(L), share: new Float64Array(L) };
    const { flux: classFlux, share } = eddy;
    const kappaLimit = eddyLimit / dt, ramp = (x, width) => Math.max(0, Math.min(1, x / width));
    const allowance = (k, i) => (Math.max(0, h[at(k, i)] - EPS) + EDDY_SLACK) * areaCell[i] / (nEdgesOnCell[i] * dt);
    for (let e = from; e < to; e++) {
      for (let k = 0; k < L; k++) eddyFlux[k * E + e] = 0;
      const kappa = Math.min(eddyKappa[e], kappaLimit);
      if (!edgeOcean[e] || !(kappa > 0)) continue;
      const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
      let shared = 0;
      for (let k = 1; k < L; k++) shared += Math.max(0, Math.min(h[at(k, a)], h[at(k, b)]) - EPS);
      const sill = sillEdge[e] + 0.5 * (eta[a] + eta[b]), coefficient = kappa * dvEdge[e] / dcEdge[e];
      let za = h[a], zb = h[b], aboveA = 0, aboveB = 0, sharedAbove = 0, upper = 0;
      for (let k = 1; k < L; k++) {
        const ha = h[at(k, a)], hb = h[at(k, b)], common = Math.max(0, Math.min(ha, hb) - EPS);
        za += ha; zb += hb;
        aboveA += Math.max(0, ha - EPS); aboveB += Math.max(0, hb - EPS); sharedAbove += common;
        let lower = 0;
        if (k < L - 1) {
          const top = ramp(Math.min(za, zb), eddyTaperDepth);
          const bottom = ramp(Math.min(sill - Math.max(za, zb), shared - sharedAbove), EDDY_BOTTOM_TAPER);
          const ceiling = ramp(Math.max(aboveA, aboveB), THIN);
          lower = -coefficient * top * bottom * ceiling * (zb - za);
          lower *= Math.min(interfaceWeight(k + 1, a), interfaceWeight(k + 1, b));
        }
        const weight = ramp(Math.max(ha, hb) - EPS, THIN);
        classFlux[k] = weight * (lower - upper);
        share[k] = weight * common;
        upper = lower;
      }
      let residual = 0, carriers = 0;
      for (let round = 0; ; round++) {
        residual = 0; carriers = 0;
        for (let k = 1; k < L; k++) { residual += classFlux[k]; carriers += share[k]; }
        if (round === 3 || !(carriers > THIN)) break;
        let held = false;
        for (let k = 1; k < L; k++) {
          const f = classFlux[k] - residual * share[k] / carriers, allowed = allowance(k, f > 0 ? a : b);
          if (Math.abs(f) > allowed) { classFlux[k] = f > 0 ? allowed : -allowed; share[k] = 0; held = true; }
        }
        if (!held) break;
      }
      if (!(carriers > THIN)) continue;
      let scale = 1;
      for (let k = 1; k < L; k++) {
        const f = classFlux[k] - residual * share[k] / carriers, allowed = allowance(k, f > 0 ? a : b);
        classFlux[k] = f;
        if (Math.abs(f) > allowed) scale = Math.min(scale, allowed / Math.abs(f));
      }
      for (let k = 1; k < L; k++) eddyFlux[k * E + e] = scale * classFlux[k];
    }
  }
  function eddyCells(from, to) {
    for (let i = from; i < to; i++) {
      for (let k = 0; k < L; k++) { const n = at(k, i); eddyVolume[n] = 0; eddyHeat[n] = 0; eddySalt[n] = 0; }
      if (!cellOcean[i]) continue;
      for (let m = 0; m < nEdgesOnCell[i]; m++) {
        const e = edgesOnCell[maxEdges * i + m], a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1], sign = i === b ? 1 : -1;
        for (let k = 1; k < L; k++) {
          const f = eddyFlux[k * E + e];
          if (f === 0) continue;
          const d = at(k, f > 0 ? a : b), n = at(k, i);
          eddyVolume[n] += sign * f; eddyHeat[n] += sign * (f * Q[d] / h[d]); eddySalt[n] += sign * (f * W[d] / h[d]);
        }
      }
    }
  }
  function eddyApply(from, to) {
    const dt = params[1];
    for (let i = from; i < to; i++) {
      const factor = dt / areaCell[i];
      for (let k = 1; k < L; k++) { const n = at(k, i); h[n] += factor * eddyVolume[n]; Q[n] += factor * eddyHeat[n]; W[n] += factor * eddySalt[n]; }
    }
  }
  function eddyTransport(dt) {
    params[1] = dt;
    if (columnRunner) { columnRunner('eddyEdges'); columnRunner('eddyCells'); columnRunner('eddyApply'); return; }
    eddyEdges(0, E); eddyCells(0, C); eddyApply(0, C);
  }

  /*
   * The ALE step of the hybrid layers (docs/ocean-hybrid-coordinate.md): each
   * sea column's target interfaces from its generator, then the edges'
   * velocities remapped piecewise-linearly onto the mean of their cells'
   * new thicknesses and the cells' heat and salt remapped by the parabolas,
   * which sets the thicknesses; shiftVelocities then restores the
   * barotropic transport. The column physics follow on the new layers:
   * mixColumn on each sea column with the cell's stress and the shear of
   * its reconstructed velocities, publishing the viscosity profile that the
   * edges' implicit solve of their normal velocity averages from its two
   * cells, and the mixed-layer and boundary-layer depths. Each pass runs
   * over a range of cells or edges (regridTargets, regridEdges,
   * regridCells, columnCells, columnEdges; prepareCells and prepareEdges
   * for the tendency's layer densities, pressure sums, free-surface
   * gradient and edge thicknesses; rescaleCells and shiftEdges; eddyEdges,
   * eddyCells and eddyApply, the eddy transport as each edge's class
   * fluxes, each cell's gather of them with its edges' signs, and their
   * application; baroPrep, baroCoriolis, baroEdges and baroCells, the
   * barotropic sub-stepping with each RK4 stage's trial state written to
   * alternating buffers and the stage, sub-step and averaging weight in
   * params[4..6]; the columns are then rescaled to the sub-step mean of
   * the free surface and shifted to the mean transport, which lags the
   * barotropic mode by half a step but is the only
   * filter the split has, the end-of-window state alone running to the
   * speed limit within months) reading the step's relax from params[2]
   * and dt from
   * params[1], so a column runner (setColumnRunner) can hand the ranges
   * to worker threads that adopt the shared arrays; counts holds the dried
   * and clamped tallies.
   */
  function cellColumn(i, hOut, TOut, SOut) {
    for (let k = 0; k < L; k++) {
      const n = at(k, i), hv = h[n];
      hOut[k] = hv;
      if (hv > 1e-9) { TOut[k] = Q[n] / hv; SOut[k] = W[n] / hv; } else if (k > 0) { TOut[k] = TOut[k - 1]; SOut[k] = SOut[k - 1]; } else { TOut[k] = labelT[k]; SOut[k] = labelS[k]; }
    }
  }
  function regridTargets(from, to) {
    const relax = params[2], { h: hCol, T: TCol, S: SCol, rho: rhoCol, z: zCol, w: wCol } = cols;
    for (let i = from; i < to; i++) {
      if (!cellOcean[i]) continue;
      cellColumn(i, hCol, TCol, SCol);
      layerDensities(TCol, SCol, L, rhoCol);
      targetInterfaces(generator, { h: hCol, rho: rhoCol, K: L, depth: D[i], table, boundaries, hMin: minimumLayer, relax, modeDepth }, zCol, wCol, cs);
      for (let k = 0; k <= L; k++) targetZ[k * C + i] = zCol[k];
    }
  }
  function regridEdges(from, to) {
    const { amount, hOld, hEdgeNew, uCol } = cols, kept = edgeTransport === 'kept';
    for (let e = from; e < to; e++) {
      if (!edgeOcean[e]) continue;
      const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
      for (let k = 0; k < L; k++) {
        hOld[k] = 0.5 * (h[at(k, a)] + h[at(k, b)]);
        hEdgeNew[k] = 0.5 * (targetZ[(k + 1) * C + a] - targetZ[k * C + a] + targetZ[(k + 1) * C + b] - targetZ[k * C + b]);
        uCol[k] = u[ae(k, e)];
      }
      const before = kept ? edgeFlow(stateThickness, e) : 0;
      remapColumn(hOld, uCol, L, hEdgeNew, amount, cs, { order: 2 });
      for (let k = 0; k < L; k++) u[ae(k, e)] = hEdgeNew[k] > 0 ? amount[k] / hEdgeNew[k] : uCol[k];
      if (kept) keepFlow(targetThickness, e, before);
    }
  }
  function regridCells(from, to) {
    const { h: hCol, T: TCol, S: SCol, hNew, amount } = cols;
    for (let i = from; i < to; i++) {
      if (!cellOcean[i]) continue;
      cellColumn(i, hCol, TCol, SCol);
      for (let k = 0; k < L; k++) hNew[k] = targetZ[(k + 1) * C + i] - targetZ[k * C + i];
      remapColumn(hCol, TCol, L, hNew, amount, cs);
      for (let k = 0; k < L; k++) Q[at(k, i)] = amount[k];
      remapColumn(hCol, SCol, L, hNew, amount, cs);
      for (let k = 0; k < L; k++) { W[at(k, i)] = amount[k]; h[at(k, i)] = hNew[k]; }
    }
  }
  function regrid(dt) {
    params[2] = regridTime > 0 && Number.isFinite(dt) ? Math.min(1, dt / regridTime) : 1;
    if (columnRunner) { columnRunner('targets'); columnRunner('uremap'); columnRunner('remap'); return; }
    regridTargets(0, C);
    regridEdges(0, E);
    regridCells(0, C);
  }
  function cellVectorAt(k, i, out, o) {
    let x = 0, y = 0, z = 0;
    for (let m = 0; m < nEdgesOnCell[i]; m++) {
      const e = edgesOnCell[maxEdges * i + m], w = 0.5 * dcEdge[e] * dvEdge[e] * u[ae(k, e)];
      x += w * nEdge[3 * e]; y += w * nEdge[3 * e + 1]; z += w * nEdge[3 * e + 2];
    }
    out[o] = x / areaCell[i]; out[o + 1] = y / areaCell[i]; out[o + 2] = z / areaCell[i];
  }
  function stressAt(i) {
    let x = 0, y = 0, z = 0;
    for (let m = 0; m < nEdgesOnCell[i]; m++) {
      const e = edgesOnCell[maxEdges * i + m], w = 0.5 * dcEdge[e] * dvEdge[e] * stress[e];
      x += w * nEdge[3 * e]; y += w * nEdge[3 * e + 1]; z += w * nEdge[3 * e + 2];
    }
    return Math.hypot(x, y, z) / areaCell[i];
  }
  function columnCells(from, to) {
    const dt = params[1], { h: hCol, T: TCol, S: SCol, rho: rhoCol, shear, kappa, nu, vector } = cols, out = { kappa, nu };
    const options = { rho0, gravity: g, mstar, ekman: ekmanFactor, rotationLimit, nstar, langmuir, shearViscosity, backgroundViscosity, backgroundDiffusivity, convectiveDiffusivity };
    for (let i = from; i < to; i++) {
      if (!cellOcean[i]) { mixedDepths[i] = 0; boundaryDepths[i] = 0; continue; }
      cellColumn(i, hCol, TCol, SCol);
      for (let k = 0; k < L; k++) cellVectorAt(k, i, vector, 3 * k);
      for (let k = 0; k < L - 1; k++) shear[k] = (vector[3 * k] - vector[3 * k + 3]) ** 2 + (vector[3 * k + 1] - vector[3 * k + 4]) ** 2 + (vector[3 * k + 2] - vector[3 * k + 5]) ** 2;
      options.guess = boundaryDepths[i]; options.langmuir = shelf[i] ? 0 : langmuir; options.waves = waves[i];
      let meltHeat = 0, meltSalt = 0;
      if (shelf[i] && meltOn) {
        const speed2 = vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2];
        let Tw = TCol[0], Sw = SCol[0], slab = 0;
        if (shelfMeltLayer > 0) {
          let sumT = 0, sumS = 0;
          for (let k = 0; k < L && slab < shelfMeltLayer; k++) { const part = Math.min(hCol[k], shelfMeltLayer - slab); sumT += part * TCol[k]; sumS += part * SCol[k]; slab += part; }
          Tw = sumT / slab; Sw = sumS / slab;
        }
        const { melt, Tb } = basalMelt(Tw, Sw, baseFreezing(g * shelfLoad[i]), Math.sqrt(topDrag * (speed2 + tidalSpeed * tidalSpeed)), meltConstants);
        const mass = melt * dt, water = mass / SHELF_MELT.freshDensity;
        meltHeat = -(latentHeatFusion * mass / rhoCp + water * (Tw - Tb)); meltSalt = -Sw * water;
        if (shelfMeltLayer > 0) {
          let given = 0;
          for (let k = 0; k < L && given < slab; k++) {
            const part = Math.min(hCol[k], slab - given);
            if (part > 0) { TCol[k] += meltHeat * part / slab / hCol[k]; SCol[k] += meltSalt * part / slab / hCol[k]; }
            given += part;
          }
          meltHeat = 0; meltSalt = 0;
        }
        melted[i] = mass; meltBase[i] = Tb; shelfMelted[i] += mass;
      }
      const result = mixColumn({ h: hCol, T: TCol, S: SCol, shear, K: L, ustar: shelf[i] ? Math.sqrt(topDrag) * Math.hypot(vector[0], vector[1], vector[2]) : Math.sqrt(stressAt(i) / rho0), f: fCell[i], dt, heat: meltHeat, salt: meltSalt, capped: shelf[i] === 1 && !(meltOn && melted[i] < 0), options }, out, cs);
      for (let k = 0; k < L; k++) { const n = at(k, i); Q[n] = hCol[k] * TCol[k]; W[n] = hCol[k] * SCol[k]; viscosity[n] = nu[k]; }
      layerDensities(TCol, SCol, L, rhoCol);
      mixedDepths[i] = mixedLayerDepth(hCol, rhoCol, L);
      boundaryDepths[i] = result.depth;
    }
  }
  function columnEdges(from, to) {
    const dt = params[1], { hOld, uCol, nuEdge } = cols, kept = edgeTransport === 'kept';
    for (let e = from; e < to; e++) {
      if (!edgeOcean[e]) continue;
      const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
      for (let k = 0; k < L; k++) { hOld[k] = 0.5 * (h[at(k, a)] + h[at(k, b)]); nuEdge[k] = 0.5 * (viscosity[at(k, a)] + viscosity[at(k, b)]); uCol[k] = u[ae(k, e)]; }
      if (!solveColumn(hOld, uCol, L, nuEdge, dt, cs)) continue;
      const before = kept ? edgeFlow(stateThickness, e) : 0;
      for (let k = 0; k < L; k++) u[ae(k, e)] = uCol[k];
      if (kept) keepFlow(stateThickness, e, before);
    }
  }
  function columnPhysics(dt) {
    params[1] = dt;
    if (columnRunner) { columnRunner('column'); columnRunner('usolve'); return; }
    columnCells(0, C);
    columnEdges(0, E);
  }

  const coverOf = (ice, concentration, i) => (ice[i] > 0 ? (concentration && concentration[i] > 0 ? concentration[i] : 1) : 0);
  function setStress(total, ice, concentration = null) {
    const cover = (i) => (shelf[i] ? 0 : coverOf(ice, concentration, i));
    if (langmuirIce === 'tanh') for (let i = 0; i < C; i++) waves[i] = 1 - Math.tanh(ICE_WAVE_DAMPING * cover(i));
    for (let e = 0; e < E; e++) stress[e] = edgeOcean[e] ? total[e] * (1 - 0.5 * (cover(cellsOnEdge[2 * e]) + cover(cellsOnEdge[2 * e + 1])) * (1 - iceStressTransmission)) * (1 - 0.5 * (shelf[cellsOnEdge[2 * e]] + shelf[cellsOnEdge[2 * e + 1]])) : 0;
  }
  function measureIce(ice, concentration, out) {
    for (let i = 0; i < C; i++) out[i] = coverOf(ice, concentration, i) * ice[i];
  }
  function readSurface(surfaceT, ice, concentration = seaCover) {
    for (let i = 0; i < C; i++) {
      if (shelf[i]) { iced[i] = 0; continue; }
      iced[i] = ice[i] > 0 ? 1 : 0;
      T0[i] = iced[i] ? FREEZING_POINT : surfaceT[i];
      surfaceIn[i] = T0[i];
      Q[i] = h[i] * T0[i];
      if (penetrated[i] !== 0) depositLight(i);
    }
  }
  function depositLight(i) {
    const light = penetrated[i];
    penetrated[i] = 0;
    if (!(transmission[i] > 0)) { Q[i] += light / rhoCp; return; }
    let z = h[i], top = transmitted(z, penetration), given = 0, last = 1;
    for (let k = 1; k < L; k++) { z += h[at(k, i)]; const bottom = transmitted(z, penetration), part = (top - bottom) / transmission[i]; shares[k] = part; given += part; top = bottom; last = k; }
    shares[last] += 1 - given;
    for (let k = 1; k < L; k++) Q[at(k, i)] += light * shares[k] / rhoCp;
  }
  function transmissionAt(i) {
    return penetration ? transmitted(h[i], penetration) : 0;
  }
  function publishTransmission() {
    transmission.fill(0);
    if (penetration) for (let i = 0; i < C; i++) if (cellOpen[i]) transmission[i] = transmissionAt(i);
  }
  /*
   * `discharge`, the part of `runoff` that leaves the land as ice, is
   * melted at the next ocean step by the top layer of its outlet.
   */
  function accumulate(evaporation, rain, dt, runoff = null, discharge = null) {
    for (let i = 0; i < C; i++) if (cellOpen[i]) fresh[i] += (evaporation ? evaporation[i] * dt : 0) - (rain ? rain[i] : 0);
    if (runoff) for (let j = 0; j < C; j++) if (!cellOpen[j] && runoff[j] > 0 && outlet[j] >= 0) fresh[outlet[j]] -= runoff[j] * areaCell[j] / areaCell[outlet[j]];
    if (discharge) for (let j = 0; j < C; j++) if (!cellOpen[j] && discharge[j] > 0 && outlet[j] >= 0) calved[outlet[j]] += discharge[j] * areaCell[j] / areaCell[outlet[j]];
  }
  function meltCalved() {
    for (let i = 0; i < C; i++) {
      if (cellOcean[i] && calved[i] !== 0) Q[i] -= latentHeatFusion * calved[i] / rhoCp;
      calved[i] = 0;
    }
  }
  /*
   * The shelves stay as they are, so what melts from them is ice the
   * sheets discharge: each sheet's calving at its outlets is its discharge
   * less its shelves' melt, and the outlets give up that much of their
   * fresh water and latent heat. Melt beyond the discharge comes out of
   * the sheet's store, kept per cell in shelfDraw; freezing beyond it
   * goes into the store. Antarctica's shelves draw on the southern
   * discharge, the rest on the northern.
   */
  function budgetShelves() {
    const { discharge, melt, moved, drawn, keep, share } = shelfBudget;
    discharge.fill(0); melt.fill(0); keep.fill(1); share.fill(0);
    for (let i = 0; i < C; i++) {
      if (shelf[i]) melt[meltPool(i)] += melted[i] * areaCell[i];
      else if (cellOpen[i] && calved[i] !== 0) discharge[meltPool(i)] += calved[i] * areaCell[i];
    }
    for (let p = 0; p < 2; p++) {
      moved[p] = discharge[p] > 0 ? Math.min(melt[p], discharge[p]) : 0;
      drawn[p] = melt[p] - moved[p];
      if (discharge[p] > 0) keep[p] = (discharge[p] - moved[p]) / discharge[p];
      if (melt[p] !== 0) share[p] = drawn[p] / melt[p];
    }
    for (let i = 0; i < C; i++) {
      if (shelf[i]) shelfDraw[i] += share[meltPool(i)] * melted[i];
      else if (cellOpen[i] && calved[i] !== 0) { const p = meltPool(i); fresh[i] += (1 - keep[p]) * calved[i]; calved[i] *= keep[p]; }
    }
  }
  function salt() {
    if (meltOn) budgetShelves();
    for (let i = 0; i < C; i++) {
      if (!cellOpen[i]) continue;
      const s = W[i] / h[i];
      W[i] += s * fresh[i] / 1000;
      const grown = (iceVolume[i] - previousVolume[i]) * iceDensity / 1000;
      W[i] += (s - iceSalinity) * grown;
      W[i] = Math.max(0, W[i]);
      fresh[i] = 0;
      previousVolume[i] = iceVolume[i];
    }
  }
  function writeSurface(surfaceT, oceanFlux, dt) {
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) { oceanFlux[i] = 0; continue; }
      if (shelf[i]) { oceanFlux[i] = 0; S0[i] = W[i] / h[i]; T0[i] = Q[i] / h[i]; previousT0[i] = T0[i]; continue; }
      capacity[i] = rhoCp * Math.max(h[i], 1);
      if (penetration) transmission[i] = transmissionAt(i);
      S0[i] = W[i] / h[i];
      if (iced[i]) {
        oceanFlux[i] = rhoCp * (Q[i] - h[i] * FREEZING_POINT) / dt;
        Q[i] = h[i] * FREEZING_POINT;
        T0[i] = FREEZING_POINT;
      } else {
        T0[i] = Q[i] / h[i];
        surfaceT[i] = T0[i];
        oceanFlux[i] = 0;
      }
      previousT0[i] = T0[i];
    }
  }

  function advance(surfaceT, ice, oceanFlux, totalStress, dt, concentration = seaCover, time = null) {
    if (++counter % everySteps !== 0) return false;
    if (!initialised) initialize(surfaceT, ice, { concentration });
    const dtOcean = everySteps * dt;
    readSurface(surfaceT, ice, concentration);
    measureIce(ice, concentration, iceVolume);
    setStress(typeof totalStress === 'function' ? totalStress() : totalStress, ice, concentration);
    step(dtOcean);
    regrid(dtOcean); columnPhysics(dtOcean);
    salt(); meltCalved();
    writeSurface(surfaceT, oceanFlux, dtOcean);
    return true;
  }

  /*
   * The start from rest: every sea column on the depth table clipped to
   * its floor (the z grid of a column held at its minima), with the
   * atlas's layer means of temperature and salinity where the atlas
   * (`climatology`, decoded by ./climatology.module.js) covers the cell,
   * the nearest covered column's where it does not, under a shelf the
   * nearest open column's at the same depth (its bottom water below its
   * floor, no colder than the in-situ freezing point), and without an atlas
   * the analytic profile T = T_d + (T₀ − T_d) e^(−z/λ), T_d 2 °C, λ 400 m
   * stretched by the thermocline tilt, the salinity from salinityProfile
   * to 34.7 the same way; the steric free surface; then one regrid, so
   * the first state is already on the generator's targets. `initialize`
   * writes the top layer's temperature into surfaceT over open sea and
   * returns how many sea cells took the atlas and how many did not.
   */
  function initialize(surfaceT, ice, { climatology: atlas = climatology, concentration = seaCover } = {}) {
    const filled = build(surfaceT, ice, atlas, concentration);
    if (!filled) return null;
    let covered = 0, other = 0;
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) continue;
      if (filled[i]) covered++; else other++;
      if (filled[i] && !iced[i]) surfaceT[i] = T0[i];
    }
    return { atlas: covered, analytic: other };
  }
  function hybridGrid(i, zCol) {
    const H = D[i] + eta[i];
    for (let k = 1; k < L; k++) zCol[k] = Math.min(table.depths[k], H);
    spaceInterfaces(zCol, L, H, minimumLayer);
    return H;
  }
  function build(surfaceT, ice, atlas, concentration = seaCover) {
    u.fill(0); eta.fill(0); fresh.fill(0); calved.fill(0); buoyancyLoss.fill(0);
    measureIce(ice, concentration, previousVolume);
    const filled = new Uint8Array(C), { z: zCol } = cols, zs = [], ts = [], ss = [];
    const integral = (v, a, b) => {
      let sum = 0;
      for (let j = 1; j < zs.length; j++) {
        const lo = Math.max(a, zs[j - 1]), hi = Math.min(b, zs[j]);
        if (hi <= lo) continue;
        sum += (v[j - 1] + (v[j] - v[j - 1]) / (zs[j] - zs[j - 1]) * (0.5 * (lo + hi) - zs[j - 1])) * (hi - lo);
      }
      return sum;
    };
    const pending = [];
    for (let i = 0; i < C; i++) {
      iced[i] = ice[i] > 0 || shelf[i] ? 1 : 0;
      if (!cellOcean[i]) { for (let k = 0; k < L; k++) { h[at(k, i)] = 0; Q[at(k, i)] = 0; W[at(k, i)] = 0; } T0[i] = iced[i] ? FREEZING_POINT : surfaceT[i]; S0[i] = referenceS; continue; }
      hybridGrid(i, zCol);
      const lid = draft[i], bed = lid + D[i];
      let column = atlas && !shelf[i] ? atlas.columnAt(latCell[i], mesh.lonCell[i]) : null;
      if (column) {
        const reached = column.depths.length;
        if (atlas.columnReaching && atlas.depths && reached < atlas.depths.length && atlas.depths[reached] <= bed) {
          const deep = atlas.columnReaching(latCell[i], mesh.lonCell[i], bed);
          if (deep && deep.depths.length > reached) {
            const Td = Float64Array.from(deep.T), Sd = Float64Array.from(deep.S);
            Td.set(column.T); Sd.set(column.S);
            column = { depths: deep.depths, T: Td, S: Sd };
          }
        }
        filled[i] = 1;
        zs.length = ts.length = ss.length = 0;
        for (let j = 0; j < column.depths.length && column.depths[j] < bed; j++) { zs.push(column.depths[j]); ts.push(column.T[j]); ss.push(column.S[j]); }
        const [tb, sb] = profileAt(column, bed);
        zs.push(bed); ts.push(tb); ss.push(sb);
        for (let k = 0; k < L; k++) { const n = at(k, i), hk = zCol[k + 1] - zCol[k]; h[n] = hk; Q[n] = integral(ts, lid + zCol[k], lid + zCol[k + 1]); W[n] = integral(ss, lid + zCol[k], lid + zCol[k + 1]); }
        if (lid > 0) shelfLoad[i] = displacedWeight(lid, (z) => profileAt(column, z), eos);
      } else if (atlas) {
        pending.push(i);
        for (let k = 0; k < L; k++) h[at(k, i)] = zCol[k + 1] - zCol[k];
      } else {
        const lat = latCell[i], s0 = salinityProfile(lat), t0 = iced[i] ? FREEZING_POINT : surfaceT[i], scale = 400 * (1 - thermoclineTilt + 2 * thermoclineTilt * Math.cos(lat) ** 2);
        const tDeep = 275.15, sDeep = 34.7, decay = (a, b) => scale * (Math.exp(-a / scale) - Math.exp(-b / scale)) / (b - a);
        if (lid > 0) shelfLoad[i] = displacedWeight(lid, (z) => [tDeep + (t0 - tDeep) * Math.exp(-z / scale), sDeep + (s0 - sDeep) * Math.exp(-z / scale)], eos);
        for (let k = 0; k < L; k++) {
          const n = at(k, i), hk = zCol[k + 1] - zCol[k], f = decay(lid + zCol[k], lid + zCol[k + 1]);
          h[n] = hk; Q[n] = hk * (tDeep + (t0 - tDeep) * f); W[n] = hk * (sDeep + (s0 - sDeep) * f);
        }
      }
    }
    for (const i of pending) {
      let source = -1, ring = [i];
      const seen = new Set(ring);
      for (let r = 0; r < 40 && source < 0 && ring.length; r++) {
        const next = [];
        for (const j of ring) for (let m = 0; m < nEdgesOnCell[j]; m++) {
          const c = cellsOnCell[maxEdges * j + m];
          if (seen.has(c)) continue;
          seen.add(c);
          if (filled[c]) { source = c; break; }
          if (cellOcean[c]) next.push(c);
        }
        ring = next;
      }
      const sourceAt = (y) => {
        let top = draft[source], j = 0, last = at(0, source);
        while (j < L - 1 && top + h[at(j, source)] < y) { if (h[at(j, source)] > 1e-9) last = at(j, source); top += h[at(j, source)]; j++; }
        return shelf[i] && !(h[at(j, source)] > 1e-9) ? last : at(j, source);
      };
      if (draft[i] > 0) shelfLoad[i] = displacedWeight(draft[i], (y) => { const m = source >= 0 ? sourceAt(y) : -1; return m >= 0 && h[m] > 1e-9 ? [Q[m] / h[m], W[m] / h[m]] : [FREEZING_POINT, referenceS]; }, eos);
      let z = draft[i];
      for (let k = 0; k < L; k++) {
        const n = at(k, i), hk = h[n], center = z + 0.5 * hk;
        z += hk;
        let t = iced[i] ? FREEZING_POINT : surfaceT[i], s = referenceS;
        if (source >= 0) {
          const m = sourceAt(center);
          if (h[m] > 1e-9) { t = Q[m] / h[m]; s = W[m] / h[m]; }
        }
        if (shelf[i]) t = Math.max(t, SHELF_MELT.liquidusSlope * s + baseFreezing(g * (shelfLoad[i] + rho0 * (center - draft[i]))));
        Q[n] = hk * t; W[n] = hk * s;
      }
    }
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) continue;
      T0[i] = shelf[i] ? Q[i] / h[i] : iced[i] ? FREEZING_POINT : Math.max(FREEZING_POINT, Q[i] / h[i]);
      if (!shelf[i]) Q[i] = h[i] * T0[i];
      S0[i] = W[i] / h[i];
      previousT0[i] = T0[i];
      capacity[i] = rhoCp * Math.max(h[i], 1);
    }
    layerDensitiesAll();
    mixedDepthsAll();
    stericSurface();
    regrid(Infinity);
    publishTransmission();
    penetrated.fill(0);
    counter = 0;
    initialised = true;
    return atlas ? filled : null;
  }
  /*
   * The free surface that levels the pressure at `referenceDepth` in every
   * column of the open ocean's abyss (abyssalCells), so the deep ocean
   * starts without a barotropic pressure gradient; every other column takes
   * the free surface of the water around it, found by relaxation from the
   * abyss's mean, as a shelf's sea level follows the ocean beside it and a
   * basin behind a sill follows the ocean outside. The ocean-mean height is
   * zero, and the height goes into the deepest layer that holds water.
   */
  function stericSurface(referenceDepth = 3500) {
    const deep = abyssalCells(mesh, D, cellOcean, referenceDepth);
    let deepArea = 0, deepMean = 0;
    for (let i = 0; i < C; i++) {
      eta[i] = 0;
      if (!deep[i]) continue;
      let budget = referenceDepth, anomaly = 0;
      for (let k = 0; k < L && budget > 0; k++) {
        const part = Math.min(h[at(k, i)], budget);
        anomaly += (layerRho[at(k, i)] - rho0) * part;
        budget -= part;
      }
      eta[i] = -anomaly / rho0;
      deepArea += areaCell[i]; deepMean += areaCell[i] * eta[i];
    }
    if (deepArea > 0) for (let i = 0; i < C; i++) if (deep[i]) eta[i] -= deepMean / deepArea;
    const next = new Float64Array(C);
    for (let sweep = 0; sweep < 2000; sweep++) {
      let moved = 0;
      for (let i = 0; i < C; i++) {
        if (!cellOcean[i] || deep[i]) { next[i] = eta[i]; continue; }
        let sum = 0, weight = 0;
        for (let m = 0; m < nEdgesOnCell[i]; m++) {
          const j = cellsOnCell[maxEdges * i + m];
          if (!cellOcean[j]) continue;
          sum += eta[j]; weight++;
        }
        next[i] = weight ? sum / weight : eta[i];
        moved = Math.max(moved, Math.abs(next[i] - eta[i]));
      }
      eta.set(next);
      if (moved < 1e-7) break;
    }
    let area = 0, mean = 0;
    for (let i = 0; i < C; i++) if (cellOcean[i]) { area += areaCell[i]; mean += areaCell[i] * eta[i]; }
    mean = area > 0 ? mean / area : 0;
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) { eta[i] = 0; continue; }
      eta[i] -= mean;
      let deepest = 0;
      for (let k = L - 1; k >= 1; k--) if (h[at(k, i)] > THIN - eta[i]) { deepest = k; break; }
      const n = at(deepest, i), t = Q[n] / h[n], sal = W[n] / h[n];
      h[n] += eta[i]; Q[n] = h[n] * t; W[n] = h[n] * sal;
    }
  }

  /*
   * As on the GPU (layeredOcean.gpu.js, restartable), a state that
   * carries every array of the GPU's restartArrays, its columns within a
   * meter of this bathymetry and dry on land, brings back the surface's
   * last temperature, the ice volume it last saw, the freshwater not yet
   * taken and the heat capacity, the sunlight not yet deposited when it
   * carries `penetrated`, the ice not yet melted when it carries
   * `calved`, the boundary-layer depths the column physics
   * search from when it carries `boundaryDepths`, and the heat and salt
   * stores as they were when it carries `Q` and `W`.
   */
  const RESTART = ['previousT0', 'previousIceVolume', 'fresh', 'buoyancy', 'flux', 'capacity'];
  function restartable(saved) {
    if (![saved.Qref ?? saved.Q, saved.Wref ?? saved.W].every((store) => store && store.length === L * C) || !RESTART.every((name) => saved[name] && saved[name].length === C)) return false;
    for (let i = 0; i < C; i++) {
      let sum = 0;
      for (let k = 0; k < L; k++) sum += saved.h[at(k, i)];
      if (cellOcean[i] ? !(Math.abs(sum - D[i] - saved.eta[i]) < 1) : sum !== 0) return false;
    }
    return true;
  }
  /*
   * A saved ocean on these labels loads as it is, its columns
   * stretched to this bathymetry's depth; a saved slab ocean is converted
   * (convertSlab): each column's temperature and salinity remapped by the
   * parabolas from its mixed layer and classes onto the depth table at its
   * floor, each edge's velocity piecewise-linearly onto the mean of its
   * cells' new thicknesses, then one regrid to the generator's targets.
   */
  function convertSlab(saved, classes) {
    const Ls = classes.length + 1;
    if (saved.h.length !== Ls * C) throw new Error(`a slab state on ${classes.length} classes holds ${Ls * C} thicknesses, not ${saved.h.length}`);
    const scratch = columnScratch(Math.max(Ls, L)), hs = new Float64Array(Ls), qs = new Float64Array(Ls), { z: zCol, hNew, amount } = cols;
    eta.set(saved.eta);
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) { for (let k = 0; k < L; k++) { h[at(k, i)] = 0; Q[at(k, i)] = 0; W[at(k, i)] = 0; } continue; }
      const H = hybridGrid(i, zCol);
      let sum = 0;
      for (let j = 0; j < Ls; j++) { hs[j] = Math.max(0, saved.h[j * C + i]); sum += hs[j]; }
      if (!(sum > 0)) { for (let k = 0; k < L; k++) { const n = at(k, i), hk = zCol[k + 1] - zCol[k]; h[n] = hk; Q[n] = hk * surfaceTAt(i); W[n] = hk * referenceS; } continue; }
      for (let j = 0; j < Ls; j++) hs[j] *= H / sum;
      for (let k = 0; k < L; k++) hNew[k] = zCol[k + 1] - zCol[k];
      for (let j = 0; j < Ls; j++) qs[j] = saved.T[j * C + i];
      remapColumn(hs, qs, Ls, hNew, amount, scratch, { layers: L });
      for (let k = 0; k < L; k++) Q[at(k, i)] = amount[k];
      for (let j = 0; j < Ls; j++) qs[j] = saved.S[j * C + i];
      remapColumn(hs, qs, Ls, hNew, amount, scratch, { layers: L });
      for (let k = 0; k < L; k++) { W[at(k, i)] = amount[k]; h[at(k, i)] = hNew[k]; }
    }
    const us = new Float64Array(Ls), hEdgeNew = cols.hEdgeNew;
    for (let e = 0; e < E; e++) {
      if (!edgeOcean[e]) { for (let k = 0; k < L; k++) u[ae(k, e)] = 0; continue; }
      const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
      for (let j = 0; j < Ls; j++) { hs[j] = 0.5 * (Math.max(0, saved.h[j * C + a]) + Math.max(0, saved.h[j * C + b])); us[j] = saved.u[j * E + e]; }
      for (let k = 0; k < L; k++) hEdgeNew[k] = 0.5 * (h[at(k, a)] + h[at(k, b)]);
      remapColumn(hs, us, Ls, hEdgeNew, amount, scratch, { order: 2, layers: L });
      for (let k = 0; k < L; k++) u[ae(k, e)] = hEdgeNew[k] > 0 ? amount[k] / hEdgeNew[k] : 0;
    }
  }
  let surfaceTAt = () => 288;
  function load(saved, surfaceT, ice, concentration = seaCover) {
    if (!saved || !saved.h) { build(surfaceT, ice, climatology, concentration); return; }
    surfaceTAt = (i) => (ice[i] > 0 ? FREEZING_POINT : surfaceT[i]);
    const classes = savedDensities(saved, C), own = classes && sameDensities(classes, labels) && saved.h.length === L * C;
    if (own) {
      h.set(saved.h); u.set(saved.u); eta.set(saved.eta);
      for (let n = 0; n < L * C; n++) { Q[n] = h[n] * saved.T[n]; W[n] = h[n] * saved.S[n]; }
      if (saved.Q && saved.Q.length === L * C && saved.W && saved.W.length === L * C) { Q.set(saved.Q); W.set(saved.W); }
    } else convertSlab(saved, classes);
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) { for (let k = 0; k < L; k++) { h[at(k, i)] = 0; Q[at(k, i)] = 0; W[at(k, i)] = 0; } eta[i] = 0; continue; }
      let sum = 0;
      for (let k = 0; k < L; k++) sum += h[at(k, i)];
      const H = D[i] + eta[i];
      if (!(sum > 0)) { hybridGrid(i, cols.z); for (let k = 0; k < L; k++) { const n = at(k, i), hk = cols.z[k + 1] - cols.z[k]; h[n] = hk; Q[n] = hk * surfaceTAt(i); W[n] = hk * referenceS; } continue; }
      if (Math.abs(sum - H) > 1e-9) { const scale = H / sum; for (let k = 0; k < L; k++) { const n = at(k, i); h[n] *= scale; Q[n] *= scale; W[n] *= scale; } }
    }
    for (let k = 0; k < L; k++) for (let e = 0; e < E; e++) if (!edgeOcean[e]) u[ae(k, e)] = 0;
    if (hasShelves) for (let i = 0; i < C; i++) shelfLoad[i] = !shelf[i] ? 0 : saved.shelfLoad && saved.shelfLoad.length === C ? saved.shelfLoad[i] : draft[i] * eos(saved.T[i], saved.S[i]);
    for (const [name, field] of [['shelfMelted', shelfMelted], ['shelfDraw', shelfDraw]]) field.set(saved[name] && saved[name].length === C ? saved[name] : new Float64Array(C));
    if (!own) regrid(Infinity);
    measureIce(ice, concentration, previousVolume);
    for (let i = 0; i < C; i++) { iced[i] = ice[i] > 0 ? 1 : 0; T0[i] = cellOcean[i] ? Q[i] / h[i] : surfaceTAt(i); S0[i] = cellOcean[i] ? W[i] / h[i] : referenceS; previousT0[i] = T0[i]; capacity[i] = rhoCp * Math.max(h[i], 1); }
    fresh.fill(0); calved.fill(0); buoyancyLoss.fill(0); penetrated.fill(0);
    layerDensitiesAll();
    mixedDepthsAll();
    publishTransmission();
    if (restartable(saved)) {
      previousT0.set(saved.previousT0); previousVolume.set(saved.previousIceVolume); fresh.set(saved.fresh); capacity.set(saved.capacity);
      if (saved.penetrated && saved.penetrated.length === C) penetrated.set(saved.penetrated);
      if (saved.calved && saved.calved.length === C) calved.set(saved.calved);
      if (saved.boundaryDepths && saved.boundaryDepths.length === C) boundaryDepths.set(saved.boundaryDepths);
    }
    counter = 0;
    initialised = true;
  }
  function restartArrays() {
    const copy = (values) => Float64Array.from(values);
    return {
      Q: copy(Q), W: copy(W), previousT0: copy(previousT0), previousIceVolume: copy(previousVolume), fresh: copy(fresh), calved: copy(calved), buoyancy: copy(buoyancyLoss), capacity: copy(capacity), penetrated: copy(penetrated), boundaryDepths: copy(boundaryDepths),
      ...(penetration ? { transmission: copy(transmission) } : {}),
    };
  }

  function serialize() {
    const Tall = new Array(L * C), Sall = new Array(L * C);
    for (let n = 0; n < L * C; n++) { const hh = Math.max(EPS, h[n]); Tall[n] = Q[n] / hh; Sall[n] = W[n] / hh; }
    return { h: Array.from(h), u: Array.from(u), T: Tall, S: Sall, eta: Array.from(eta), densities: Array.from(labels), vertical, generator, ...(hasShelves ? { shelfLoad: Array.from(shelfLoad), shelfMelted: Array.from(shelfMelted), shelfDraw: Array.from(shelfDraw) } : {}) };
  }

  const thermoclineDepth = new Float64Array(C), sst = T0, sss = S0;
  function thermocline(i) {
    let z = 0, previous = -Infinity, center = 0;
    for (let k = 0; k < L; k++) {
      const n = at(k, i), hv = h[n], rk = Math.max(previous, layerRho[n]), c = z + 0.5 * hv;
      if (rk >= THERMOCLINE_DENSITY) return k === 0 ? 0 : center + (THERMOCLINE_DENSITY - previous) / Math.max(rk - previous, 1e-12) * (c - center);
      previous = rk; center = c; z += hv;
    }
    return z;
  }
  function fields(depth = 0) {
    layerDensitiesAll();
    for (let i = 0; i < C; i++) thermoclineDepth[i] = cellOcean[i] ? thermocline(i) : NaN;
    const layerT = (k, i) => (k === 0 ? T0[i] : Q[at(k, i)] / h[at(k, i)]);
    return { h1: mixedDepths, T1: T0, S1: S0, u1: u.subarray(0, E), eta, thermoclineDepth, ...depthFields(mesh, L, { h, u, temperature: layerT, cellOcean }, depth, { thin: 0 }) };
  }

  function diagnostics() {
    layerDensitiesAll();
    let area = 0, depth = 0, heat = 0, thermo = 0, salinity = 0, ssh = 0, speed = 0, transport = 0, interiorT = 0, interiorH = 0;
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) continue;
      const a = areaCell[i];
      area += a; depth += a * mixedDepths[i]; salinity += a * (W[i] / h[i]); ssh = Math.max(ssh, Math.abs(eta[i]));
      thermo += a * thermocline(i);
      for (let k = 0; k < L; k++) heat += a * rhoCp * Q[at(k, i)];
      for (let k = 1; k < L; k++) { interiorT += Q[at(k, i)] * a; interiorH += h[at(k, i)] * a; }
    }
    for (let e = 0; e < E; e++) {
      speed = Math.max(speed, Math.abs(u[e]));
      let t = 0;
      for (let k = 0; k < L; k++) t += hEdge[ae(k, e)] * u[ae(k, e)];
      transport = Math.max(transport, Math.abs(t) * dvEdge[e]);
    }
    return { oceanUpperDepth: depth / area, oceanHeat: heat / area, oceanInteriorT: interiorH > 0 ? interiorT / interiorH : 0, oceanSpeed: speed, oceanThermoclineDepth: thermo / area, oceanSalinity: salinity / area, oceanSSH: ssh, oceanTransport: transport / 1e6, oceanLimited: limited, oceanDried: dried };
  }

  if (!adopting) build(new Float64Array(C).fill(288), new Float64Array(C), null);
  initialised = false;
  const sharedBuffers = { h: h.buffer, u: u.buffer, Q: Q.buffer, W: W.buffer, eta: eta.buffer, rhoMl: rhoMl.buffer, capacity: capacity.buffer, transmission: transmission.buffer, penetrated: penetrated.buffer, stress: stress.buffer, waves: waves.buffer, hEdge: hEdge.buffer, pressure: pressure.buffer, gradEta: gradEta.buffer, params: params.buffer, trialh: trial[0].buffer, trialu: trial[1].buffer, trialQ: trial[2].buffer, trialW: trial[3].buffer };
  Object.assign(sharedBuffers, { baroSlow: slow.buffer, baroU: U.buffer, baroDepth: depthEdge.buffer, baroEta: etaB.buffer, baroAvgU: avgU.buffer, baroAvgEta: avgEta.buffer, baroTrialEtaA: bTrial[0].buffer, baroTrialUA: bTrial[1].buffer, baroTrialEtaB: bTrialB[0].buffer, baroTrialUB: bTrialB[1].buffer, counts: countBuffer });
  bStages.forEach((stage, s) => { sharedBuffers[`baroStage${s}Eta`] = stage[0].buffer; sharedBuffers[`baroStage${s}U`] = stage[1].buffer; });
  Object.assign(sharedBuffers, { layerRho: layerRho.buffer, columnAbove: columnAbove.buffer, targetZ: targetZ.buffer, viscosity: viscosity.buffer, edgeBottom: edgeBottom.buffer, edgeAbove: edgeAbove.buffer, mixedDepths: mixedDepths.buffer, boundaryDepths: boundaryDepths.buffer, shelfLoad: shelfLoad.buffer, eddyFlux: eddyFlux.buffer, eddyVolume: eddyVolume.buffer, eddyHeat: eddyHeat.buffer, eddySalt: eddySalt.buffer });
  stages.forEach((stage, s) => { sharedBuffers[`stage${s}h`] = stage[0].buffer; sharedBuffers[`stage${s}u`] = stage[1].buffer; sharedBuffers[`stage${s}Q`] = stage[2].buffer; sharedBuffers[`stage${s}W`] = stage[3].buffer; });
  return { state, trial, stages, layers: L, everySteps, h, u, Q, W, eta, D, sillEdge, hEdge, T0, S0, capacity, transmission, penetrated, penetration, stress, fresh, calved, avgU, tendency, prepare, tendencyLayers, setLayerRunner(fn) { layerRunner = fn; }, setColumnRunner(fn) { columnRunner = fn; }, useSeaIce(concentration) { seaCover = concentration; }, advance, accumulate, meltCalved, melted, meltBase, shelfMelted, shelfDraw, shelfBudget, initialize, load, serialize, restartArrays, diagnostics, fields, setStress, readSurface, salt, writeSurface, phases, eddyTransport, rhoCp, edgeOcean, cellOcean, densities: rho, shared: sharedBuffers, vertical, generator, table, mixedDepths, boundaryDepths, viscosity, layerRho, regrid, columnPhysics, regridTargets, regridEdges, regridCells, columnCells, columnEdges, prepareCells, prepareEdges, rescaleCells, shiftEdges, eddyEdges, eddyCells, eddyApply, baroPrep, baroCoriolis, baroEdges, baroCells };
}
