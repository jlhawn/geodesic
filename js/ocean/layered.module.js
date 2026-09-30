import { divergence, gradient, curl, kineticEnergy, laplacianVelocity, cellVector } from '../dynamics/operators.module.js';
import { FREEZING_POINT } from '../physics/ice.module.js';
import { SEAWATER, seawaterDensity, thermalExpansion as expansionOf, halineContraction as contractionOf, labelTemperature, salinityForDensity } from './seawater.module.js';
import { profileAt } from './climatology.module.js';

/*
 * A layered ocean on the C-grid: a bulk mixed layer with its own
 * temperature and salinity over interior layers of fixed density, the
 * hybrid isopycnal design of MICOM, on the real bathymetry, with the
 * density of seawater/seawater.module.js. Each layer
 * is a TRiSK shallow-water layer in the vector-invariant form, carrying
 * thickness, edge velocity, heat h·T and salt h·S. The pressure force
 * in an interior layer k is the gradient of a cell potential,
 *   Φ_k = g η + (g/ρ₀)[ρ_ml h₀ − ρ_k h₀ + Σ_{j<k} (ρ_j − ρ_k) h_j],
 * which is exact on the mesh whatever the bathymetry; the mixed layer's
 * force is g∇η plus (g/ρ₀)(h₀/2)∇ρ_ml, the exact layer mean of the
 * hydrostatic gradient whatever its thickness. A layer
 * that has outcropped, or lies below the bottom, keeps a token thickness
 * at its label, made up from the mixed layer, which takes the difference
 * in water, heat and salt, and follows the velocity of the layer above.
 *
 * A mixed layer far deeper than its neighbours is therefore no trouble
 * for the pressure force, nor for the free surface, which follows the
 * transport of the centred edge thicknesses while the mixed layer's own
 * mass moves with the donor cell's (the rescaling below carries the
 * difference through the whole column). What keeps such a column
 * balanced is the Coriolis term. The mixed layer's takes the centred
 * flux h_e u: with the donor-limited flux, where a thin mixed layer
 * feeds a deep one, it feels a fraction of the Coriolis force its
 * potential vorticity implies and runs down the sea-level gradient. An
 * interior layer's potential vorticity is built on the smaller of its
 * two cells' thicknesses, as its flux is, but on no less than
 * `vorticityCentring` times their mean: the layer under a deep mixed
 * layer is often a remnant a few metres thick beside the full layer next
 * door, and on the smaller thickness alone the Coriolis force it
 * assembles from its thick neighbours' fluxes is many times what its
 * pressure gradient balances, a jet at metres a second.
 *
 * The free surface η = Σh − D moves at √(gD), too fast for the ocean's
 * step, so the depth-integrated flow is sub-stepped: the baroclinic
 * layers take one RK4 step with η frozen at its starting value, the
 * barotropic transport and η take RK4 sub-steps within the same
 * interval forced by the depth integral of the slow tendencies, and the
 * layers are then rescaled to the sub-stepped η and shifted to the
 * averaged transport.
 *
 * The mixed layer exchanges mass with the interior after each step,
 * driven by its surface buoyancy loss B, from the step's heat, freshwater
 * and ice growth, remembered over `buoyancyMemory` so that a day's
 * sunshine does not undo a winter's convection. Each interior class
 * stands for water spanning half-way to its neighbours' labels: a mixed
 * layer denser than all of it swallows it at up to `convectiveRate`; one
 * within that span erodes it while losing buoyancy at B/(h₀N²), N² the
 * class's remaining span of density over its thickness, the deepening of
 * convection into stratified water, so that its depth follows the
 * season's buoyancy loss, not a whole class at a time where the surface
 * water crosses a label; and it entrains the layer below at the Kraus–Turner
 * wind-stirring rate, the stirring fading with depth over
 * `stirringDepth`. It deepens within `maximumMixedDepth` and, when
 * `mixedNeighbourRatio` is set, within that multiple of the mean depth of
 * the neighbouring mixed layers. Otherwise it holds its depth,
 * convectively neutral or not: it detrains at once what lies beyond the
 * maximum, and, over `detrainmentTime` and never above
 * `shallowestMixedDepth`, what lies beyond the Monin–Obukhov depth
 * 2 m u*³/(−B) while it gains buoyancy and what lies beyond its
 * neighbours' reach. `neutralSnap` returns everything below
 * `shallowestMixedDepth` whenever the layer is as dense as the water
 * beneath it, and `convectiveErosion: false` swallows every class no
 * denser than the mixed layer at `convectiveRate`.
 * Detrained water goes to the interior layer whose density is nearest
 * its own, so water swallowed from a layer returns to that layer; a
 * mass-conserving split between the two bracketing layers instead
 * ratcheted a fraction of every swallow-and-return cycle into the
 * denser class. After the exchanges each interior layer more than
 * RESTORE_TOLERANCE from its label density mixes in water from the
 * nearest layer that lies clearly on the other side of the label, a
 * fraction dt/restoreTime of the full correction per step, so its
 * temperature and salinity stay those of its class; the curvature of the
 * equation of state makes a mixture denser than the linear estimate, and
 * the gradual, dead-banded correction keeps that from overshooting.
 * Tracers are carried by the flux with the donor cell's value. Surface
 * freshwater (evaporation minus rain, and runoff where it flows down to the sea) and ice
 * growth or melt act on its salinity as virtual salt fluxes. The mixed
 * layer's temperature is the sea surface temperature the atmosphere sees, its
 * heat capacity is published per cell, and the heat converged under ice
 * is handed to the ice base.
 *
 * Mesoscale eddies, 10–30 km across and unresolved on these meshes, flatten
 * the interfaces between the interior classes (Gent and McWilliams 1990,
 * in the interface-height form of MICOM and HYCOM). After the dynamics of
 * each step, the water above interior interface k+½, at depth
 * z = Σ_{j≤k} h_j below the free surface, crosses edge e from cell a to b at
 *   G_{k+½} = −κ τ (z_b − z_a) dv_e/dc_e,
 * and class k carries G_{k+½} − G_{k−½}, so every interface diffuses and the
 * column sum is unchanged; the base of the mixed layer and the sea floor
 * carry none, so the mixed layer is left to its own exchanges. Heat and salt
 * go with the donor cell's water and the velocities are left as they are.
 * κ is `eddyDiffusivity` wherever the mesh is much coarser than the first
 * baroclinic deformation radius L_d = c/√(f² + 2βc), c = EDDY_WAVE_SPEED,
 * and falls as 1/(1 + (L_d/dc)²) where it begins to resolve it (Hallberg
 * 2013), which is only near the equator at 112 or 56 km. The taper τ is
 * the product of linear ramps: over the top `eddyTaperDepth` metres, on the
 * shallower cell's depth, so the scheme does not act in or just under the
 * mixed layer; over EDDY_BOTTOM_TAPER metres above the edge's sill and of
 * the interior water below the interface that both columns hold (the
 * classes that can carry the return flow across the edge), so an interface
 * that meets the bottom on either side carries nothing; and over THIN
 * metres of interior water above the interface in the fuller column, so
 * interfaces that are really the mixed-layer base carry nothing, while a
 * class that has outcropped on one side can still spread under the mixed
 * layer from the other. A class with less than THIN metres on both sides
 * carries that fraction of its flux, and the rest is spread over the
 * classes present on both sides in proportion to their thickness there, so
 * a token class keeps exactly its token and the interfaces on either side
 * of it move together; an edge whose two columns share less than THIN
 * metres of interior water carries nothing. A class that would send more
 * than 1/nEdges of its water above the token thickness out of a cell
 * through one edge in the step is held to that and the difference spread
 * in the same way, for up to three rounds, after which any excess scales
 * the edge's fluxes together; κ is held below a quarter of the explicit
 * diffusion limit 1/(dt·max Σ dv/(dc·A)).
 */
/*
 * The interior classes, placed where the World Ocean Atlas holds its water
 * in this equation of state: 1020.5–1021.5 for the warm pool, 0.25 apart
 * through the thermocline, 0.1 apart for the intermediate water to 1026.6,
 * 0.03 apart from 1026.65 to 1026.98 so that the Southern Ocean's winter
 * water (1026.79), Circumpolar Deep Water (1026.85), North Atlantic Deep
 * Water (1026.89) and Antarctic Bottom Water (1026.93) lie in separate
 * classes, and a few steps for the Arctic's deep and shelf water and the
 * Red Sea's and Mediterranean's. LAYER_SALINITIES label each class and
 * LAYER_BOTTOMS give its base in the subtropics for the analytic start;
 * the classes from 1026.95 lie below any sea floor there, so the 1026.95
 * class fills the deep columns and the denser ones start as tokens.
 */
export const LAYER_DENSITIES = [
  1020.5, 1021.0, 1021.5,
  1022.0, 1022.25, 1022.5, 1022.75, 1023.0, 1023.25, 1023.5, 1023.75, 1024.0, 1024.25, 1024.5, 1024.75, 1025.0, 1025.25, 1025.5, 1025.75, 1026.0,
  1026.1, 1026.2, 1026.3, 1026.4, 1026.5, 1026.6,
  1026.65, 1026.68, 1026.71, 1026.74, 1026.77, 1026.8, 1026.83, 1026.86, 1026.89, 1026.92, 1026.95, 1026.98,
  1027.05, 1027.15, 1027.3, 1027.5, 1027.75, 1028.0,
];
export const LAYER_BOTTOMS = [
  65, 74, 83,
  90, 110, 130, 150, 170, 205, 235, 270, 300, 350, 400, 450, 500, 550, 600, 650, 689,
  734, 801, 867, 933, 1004, 1061,
  1092, 1180, 1300, 1420, 1540, 1735, 2005, 2275, 2650, 3550, 12000, 12000,
  12000, 12000, 12000, 12000, 12000,
];
export const LAYER_SALINITIES = [
  34.0, 34.5, 34.85,
  35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 34.98, 34.95, 34.92, 34.9,
  34.89, 34.88, 34.875, 34.87, 34.86, 34.85,
  34.843, 34.839, 34.835, 34.831, 34.826, 34.82, 34.814, 34.809, 34.806, 34.803, 34.8, 34.82,
  34.91, 34.94, 35.05, 39.9, 39.0, 38.65,
];
/*
 * The interior classes of saved oceans that carry no list of their own,
 * told apart by their layer count: the seven of the page's saved runs and
 * the twenty-three of later ones.
 */
export const UNLISTED_LAYER_DENSITIES = [
  [1022.0, 1023.0, 1024.0, 1025.0, 1026.0, 1026.6, 1026.95],
  [1022.0, 1022.25, 1022.5, 1022.75, 1023.0, 1023.25, 1023.5, 1023.75, 1024.0, 1024.25, 1024.5, 1024.75, 1025.0, 1025.25, 1025.5, 1025.75, 1026.0, 1026.2, 1026.4, 1026.6, 1026.75, 1026.85, 1026.95],
];
export const THERMOCLINE_DENSITY = 1024.0;
export const POLAR_INTERIOR_T = 272.15;

/*
 * The fields at a depth below the surface, per cell: the layer holding
 * that depth (the first whose base lies below it, interior layers no
 * thicker than THIN being tokens and passed over), its temperature, its
 * current as a cell vector and as a speed, and the vertical velocity
 * there, upward positive, as the divergence of the transport above the
 * depth with the free surface held: h and u are the L layers'
 * thicknesses and edge velocities, temperature(k, i) the layer
 * temperature. Cells that are land or whose water column ends above the
 * depth are NaN, their current vector zero.
 */
export function depthFields(mesh, L, { h, u, temperature, cellOcean }, depth) {
  const { nCells: C, nEdges: E, maxEdges, nEdgesOnCell, edgesOnCell, edgeSignOnCell, cellsOnEdge, dcEdge, dvEdge, nEdge, areaCell } = mesh;
  const out = { temperature: new Float32Array(C), current: new Float32Array(3 * C), speed: new Float32Array(C), upwelling: new Float32Array(C) };
  for (let i = 0; i < C; i++) {
    let layer = -1, top = 0;
    if (cellOcean[i]) for (let k = 0; k < L; k++) { const hk = h[k * C + i]; if (k > 0 && hk <= THIN) continue; if (depth < top + hk) { layer = k; break; } top += hk; }
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
        if (k > 0 && he <= THIN) continue;
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

/*
 * The starting water of an interior class at a latitude: its class
 * temperature and salinity equatorward of 45°, blending over 15° of
 * latitude toward POLAR_INTERIOR_T at the salinity that keeps its
 * density, as the polar oceans hold cold, fresh water on the same
 * density surfaces as the warm, salty subtropical thermocline; the
 * blend is complete at 60°, where the winter ice edges lie.
 */
export function interiorWater(rho, t, s, lat) {
  const x = Math.max(0, Math.min(1, (Math.abs(lat) * 180 / Math.PI - 45) / 15));
  const weight = x * x * (3 - 2 * x);
  if (weight === 0 || t <= POLAR_INTERIOR_T) return [t, s];
  const tp = t - (t - POLAR_INTERIOR_T) * weight;
  return [tp, salinityForDensity(rho, tp)];
}
export const EPS = 0.01, THIN = 5, PV_FLOOR = 20, SPEED_LIMIT = 5, DENSITY_TOLERANCE = 0.005, RESTORE_TOLERANCE = 0.01;

/*
 * Columns from an ocean climatology (./climatology.module.js), sampled at
 * each sea cell's centre down to its bottom D, the atlas's deepest values
 * carried on down where the model is deeper. The profile's potential
 * density, held from decreasing with depth, places the interfaces: class
 * k holds the water denser than half-way to the label above it and no
 * denser than half-way to the label below, the lightest class taking all
 * lighter water and the densest all denser. The mixed layer reaches to
 * where the density first exceeds the surface's by `mixedExcess`, within
 * shallowestMixedDepth, maximumMixedDepth and D, and the classes begin
 * below it. A class with less than EPS of water keeps the EPS token at its
 * label, the thickest layer giving up the tokens' metres so the column
 * still sums to D. Every other layer takes the atlas's mean temperature
 * and salinity over its depths, and a class whose mean lies more than
 * `tolerance` from its label density keeps that temperature at the
 * salinity that gives the label's density. The mixed layer is at the
 * freezing point under ice and no colder in open water. Fills h, Q and W
 * (layer-major, L = rho.length) and T0 of the cells the atlas covers and
 * returns which cells those are.
 */
export function atlasColumns(mesh, atlas, { D, cellOcean, ice, rho, labelT, labelS, h, Q, W, T0, shallowestMixedDepth = 50, maximumMixedDepth = 600, mixedExcess = 0.03, tolerance = RESTORE_TOLERANCE }) {
  const C = mesh.nCells, L = rho.length, at = (k, i) => k * C + i;
  const filled = new Uint8Array(C), top = new Float64Array(L + 1);
  const z = [], t = [], s = [], sigma = [];
  const integral = (v, a, b) => {
    let sum = 0;
    for (let j = 1; j < z.length; j++) {
      const lo = Math.max(a, z[j - 1]), hi = Math.min(b, z[j]);
      if (hi <= lo) continue;
      const slope = (v[j] - v[j - 1]) / (z[j] - z[j - 1]);
      sum += (v[j - 1] + slope * (0.5 * (lo + hi) - z[j - 1])) * (hi - lo);
    }
    return sum;
  };
  const depthOf = (target) => {
    if (sigma[0] >= target) return 0;
    for (let j = 1; j < z.length; j++) if (sigma[j] >= target) return z[j - 1] + (target - sigma[j - 1]) / (sigma[j] - sigma[j - 1]) * (z[j] - z[j - 1]);
    return z[z.length - 1];
  };
  for (let i = 0; i < C; i++) {
    if (!cellOcean[i]) continue;
    const column = atlas.columnAt(mesh.latCell[i], mesh.lonCell[i]);
    if (!column) continue;
    filled[i] = 1;
    const bottom = D[i];
    z.length = t.length = s.length = sigma.length = 0;
    for (let j = 0; j < column.depths.length && column.depths[j] < bottom; j++) { z.push(column.depths[j]); t.push(column.T[j]); s.push(column.S[j]); }
    const [tb, sb] = profileAt(column, bottom);
    z.push(bottom); t.push(tb); s.push(sb);
    for (let j = 0; j < z.length; j++) sigma.push(Math.max(j > 0 ? sigma[j - 1] : -Infinity, seawaterDensity(t[j], s[j])));
    const mixed = Math.min(bottom, maximumMixedDepth, Math.max(shallowestMixedDepth, depthOf(sigma[0] + mixedExcess)));
    top[1] = mixed; top[L] = bottom;
    for (let k = 2; k < L; k++) top[k] = Math.min(bottom, Math.max(mixed, depthOf(0.5 * (rho[k - 1] + rho[k]))));
    const iced = ice[i] > 0;
    T0[i] = iced ? FREEZING_POINT : Math.max(FREEZING_POINT, integral(t, 0, mixed) / mixed);
    h[i] = mixed; Q[i] = mixed * T0[i]; W[i] = integral(s, 0, mixed);
    let excess = 0, thickest = 0;
    for (let k = 1; k < L; k++) {
      const n = at(k, i), hk = top[k + 1] - top[k];
      if (hk < EPS) { h[n] = EPS; Q[n] = EPS * labelT[k]; W[n] = EPS * labelS[k]; excess += EPS - hk; continue; }
      const tk = integral(t, top[k], top[k + 1]) / hk;
      let sk = integral(s, top[k], top[k + 1]) / hk;
      if (Math.abs(seawaterDensity(tk, sk) - rho[k]) > tolerance) sk = salinityForDensity(rho[k], tk);
      h[n] = hk; Q[n] = hk * tk; W[n] = hk * sk;
      if (hk > h[at(thickest, i)]) thickest = k;
    }
    const n = at(thickest, i), f = (h[n] - excess) / h[n];
    h[n] -= excess; Q[n] *= f; W[n] *= f;
  }
  return filled;
}

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
 * The model's bathymetry: the cell-mean ETOPO depth of every sea cell,
 * at least `minimumDepth`, and never shallower than `neighbourRatio`
 * times its deepest sea neighbour, so that no shelf break or trench
 * wall drops by more than that ratio across one edge. A cell 110 km
 * across cannot hold a real shelf, and a ten-to-one step in depth
 * between two cells makes the free surface swing wildly at the coast.
 */
/*
 * Fits loaded columns to this mesh's bathymetry: a column carried over
 * from another resolution keeps its interface depths from the top down
 * and is cut or extended at the bottom, a sea cell with no usable water
 * (new coast, or a stencil with no sea source) takes the climatology
 * column (built on first need), the mixed layer keeps its floor, and a
 * column that already fits is left exactly as it is.
 */
export function fitColumns({ h, Q, W, eta }, climatology, { D, cellOcean, L, C, labelT, labelS, minimumThickness }) {
  const at = (k, i) => k * C + i;
  for (let i = 0; i < C; i++) {
    if (!cellOcean[i]) { for (let k = 0; k < L; k++) { h[at(k, i)] = 0; Q[at(k, i)] = 0; W[at(k, i)] = 0; } eta[i] = 0; continue; }
    let sum = 0, valid = Number.isFinite(eta[i]);
    for (let k = 0; k < L && valid; k++) { const n = at(k, i); if (!(h[n] >= 0) || !Number.isFinite(Q[n]) || !Number.isFinite(W[n])) valid = false; else sum += h[n]; }
    if (valid && Math.abs(sum - D[i] - eta[i]) < 1e-6) continue;
    if (!valid || sum < 1) {
      const clim = typeof climatology === 'function' ? (climatology = climatology()) : climatology;
      let s = 0;
      for (let k = 0; k < L; k++) { const n = at(k, i); h[n] = clim.h[n]; Q[n] = clim.Q[n]; W[n] = clim.W[n]; s += h[n]; }
      eta[i] = s - D[i];
      continue;
    }
    eta[i] = Math.max(-5, Math.min(5, eta[i]));
    let remaining = D[i] + eta[i];
    for (let k = 0; k < L; k++) {
      const n = at(k, i), take = Math.min(h[n], Math.max(remaining, 0));
      if (take < EPS - 1e-6) { h[n] = EPS; Q[n] = EPS * labelT[k]; W[n] = EPS * labelS[k]; }
      else { const f = take / h[n]; Q[n] *= f; W[n] *= f; h[n] = take; }
      remaining -= h[n];
    }
    if (Math.abs(remaining) < 1e-9) remaining = 0;
    if (remaining > 0) {
      let deepest = 0;
      for (let k = L - 1; k > 0; k--) if (h[at(k, i)] > THIN) { deepest = k; break; }
      const n = at(deepest, i), f = (h[n] + remaining) / h[n];
      Q[n] *= f; W[n] *= f; h[n] += remaining;
    } else {
      for (let k = L - 1; k >= 0 && remaining < 0; k--) {
        const n = at(k, i), cut = Math.min(h[n] - EPS, -remaining);
        if (cut > 0) { const f = (h[n] - cut) / h[n]; Q[n] *= f; W[n] *= f; h[n] -= cut; remaining += cut; }
      }
    }
    for (let k = 1; k < L && h[i] < minimumThickness; k++) {
      const n = at(k, i), available = h[n] - EPS;
      if (available <= 0) continue;
      const amount = Math.min(available, minimumThickness - h[i]), f = amount / h[n];
      Q[i] += Q[n] * f; W[i] += W[n] * f; h[i] += amount;
      Q[n] -= Q[n] * f; W[n] -= W[n] * f; h[n] -= amount;
    }
  }
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
 * A saved ocean (serialize's layer-major h, u, T and S, and eta) carried
 * from the interior classes `from` onto `to`, column by column, the mixed
 * layer copied through. An old class stands for water spread over its
 * span, half-way to its neighbours' labels; the lightest and densest,
 * which hold all lighter and all denser water, span the same width
 * centred on their own water's density. Each new class, the lightest and
 * densest open-ended, takes the share of that span within its own, the
 * shares tilted linearly across the span (and none below zero) so that
 * their labels average to the density of the water. Each share keeps the
 * old class's temperature and salinity moved by the least change, counted
 * in 4 K and 0.5 psu (their spread along a density surface), that changes
 * the density by ρ_new − ρ̄, ρ̄ the shares' mean label: it lies near its
 * new label, warmer or colder in the thermocline and saltier or fresher
 * in cold water, and the column keeps its water, heat and salt. A class
 * left with less than the EPS token is made up to it at its label from the
 * column's thickest class, or from the mixed layer when no class can spare
 * the water. A new class moves with the old classes whose label spans
 * overlap its own, weighted by the overlap, or else with the old class of
 * the nearest label.
 */
const SPREAD_T = 4, SPREAD_S = 0.5;
export function rebinOcean(saved, from, to, C, { labelT, labelS }) {
  const Lo = from.length + 1, Ln = to.length + 1, E = saved.u.length / Lo;
  const spans = (labels) => labels.map((r, k) => {
    const below = k > 0 ? 0.5 * (labels[k - 1] + r) : null, above = k < labels.length - 1 ? 0.5 * (r + labels[k + 1]) : null;
    return [below ?? 2 * r - above, above ?? 2 * r - below];
  });
  const oldSpans = spans(from), newSpans = spans(to), receiving = newSpans.map(([a, b], n) => [n === 0 ? -Infinity : a, n === to.length - 1 ? Infinity : b]);
  const h = new Float64Array(Ln * C), Q = new Float64Array(Ln * C), W = new Float64Array(Ln * C);
  const shares = [];
  for (let i = 0; i < C; i++) {
    h[i] = saved.h[i]; Q[i] = h[i] > 0 ? h[i] * saved.T[i] : 0; W[i] = h[i] > 0 ? h[i] * saved.S[i] : 0;
    let column = h[i];
    for (let j = 1; j < Lo; j++) {
      const n = j * C + i, hj = saved.h[n], t = saved.T[n], s = saved.S[n];
      if (!(hj >= 0) || (hj > 0 && !(Number.isFinite(t) && Number.isFinite(s)))) { h[i] = NaN; continue; }
      if (hj === 0) continue;
      column += hj;
      let [a, b] = oldSpans[j - 1];
      if (j === 1 || j === Lo - 1) { const r = seawaterDensity(t, s), half = 0.5 * (b - a); a = r - half; b = r + half; }
      shares.length = 0;
      let total = 0;
      for (let m = 0; m < to.length; m++) {
        const f = Math.min(b, receiving[m][1]) - Math.max(a, receiving[m][0]);
        if (f > 0) { shares.push([m, f]); total += f; }
      }
      let centre = 0, spread = 0, mean = 0, kept = 0;
      for (const share of shares) { share[1] /= total; centre += share[1] * to[share[0]]; }
      for (const [m, f] of shares) spread += f * (to[m] - centre) ** 2;
      const tilt = spread > 0 ? (seawaterDensity(t, s) - centre) / spread : 0;
      for (const share of shares) { share[1] = Math.max(0, share[1] * (1 + tilt * (to[share[0]] - centre))); kept += share[1]; }
      for (const share of shares) { share[1] /= kept; mean += share[1] * to[share[0]]; }
      const byT = -expansionOf(t, s) * SEAWATER.rho0, byS = contractionOf(t, s) * SEAWATER.rho0, norm = (byT * SPREAD_T) ** 2 + (byS * SPREAD_S) ** 2;
      for (const [m, f] of shares) {
        const k = (m + 1) * C + i, part = f * hj, change = (to[m] - mean) / norm;
        h[k] += part; Q[k] += part * (t + byT * SPREAD_T ** 2 * change); W[k] += part * (s + byS * SPREAD_S ** 2 * change);
      }
    }
    if (!(column > 0)) continue;
    let thickest = 1, needed = 0;
    for (let k = 1; k < Ln; k++) {
      if (h[k * C + i] > h[thickest * C + i]) thickest = k;
      needed += Math.max(0, EPS - h[k * C + i]);
    }
    const giver = h[thickest * C + i] > needed + EPS ? thickest * C + i : i;
    for (let k = 1; k < Ln; k++) {
      const n = k * C + i, need = EPS - h[n];
      if (!(need > 0)) continue;
      h[n] += need; Q[n] += need * labelT[k]; W[n] += need * labelS[k];
      h[giver] -= need; Q[giver] -= need * labelT[k]; W[giver] -= need * labelS[k];
    }
  }
  const u = new Float64Array(Ln * E);
  for (let e = 0; e < E; e++) u[e] = saved.u[e];
  for (let m = 0; m < to.length; m++) {
    const [A, B] = newSpans[m], weights = [];
    let total = 0, nearest = 0;
    for (let j = 0; j < from.length; j++) {
      const w = Math.min(B, oldSpans[j][1]) - Math.max(A, oldSpans[j][0]);
      if (w > 0) { weights.push([j, w]); total += w; }
      if (Math.abs(from[j] - to[m]) < Math.abs(from[nearest] - to[m])) nearest = j;
    }
    if (!weights.length) weights.push([nearest, total = 1]);
    const out = (m + 1) * E;
    for (const [j, w] of weights) for (let e = 0; e < E; e++) u[out + e] += (w / total) * saved.u[(j + 1) * E + e];
  }
  const T = Float64Array.from(Q, (q, n) => (h[n] > 0 ? q / h[n] : 0)), S = Float64Array.from(W, (w, n) => (h[n] > 0 ? w / h[n] : 0));
  return { h, u, T, S, eta: Float64Array.from(saved.eta), densities: Array.from(to) };
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

export function createOcean(mesh, {
  densities = LAYER_DENSITIES, salinities = LAYER_SALINITIES, bottoms = LAYER_BOTTOMS, mixedDepth = 60, minimumDepth = 50, flatDepth = 4000, thermoclineTilt = 0.3,
  salinityProfile = (lat) => 34 + 2 * Math.exp(-(((Math.abs(lat) * 180 / Math.PI - 25) / 20) ** 2)),
  density = 1025, specificHeat = 3985, referenceS = 35, gravity = 9.81,
  minimumThickness = 50, shallowestMixedDepth = 50, maximumMixedDepth = 600, convectiveRate = 100 / 86400, neutralSnap = false, convectiveErosion = true, buoyancyMemory = 86400, mixedNeighbourRatio = 0, vorticityCentring = 0.5, stirring = 0.8, stirringDepth = 100, detrainmentTime = 86400, restoreTime = 2 * 86400, iceStressTransmission = 0.8, iceSalinity = 5, iceDensity = 917,
  interfacialDrag = 2e-4, bottomDrag = 3e-3, closureHours = 12, closureSpacing = CLOSURE_SPACING, diffusivity = 0.01, everySteps = 4,
  eddyDiffusivity = 1000, eddyTaperDepth = 200,
  geography = null, bathymetry = null, buffers = null, climatology = null,
} = {}) {
  const {
    nCells: C, nEdges: E, nVertices: V, maxEdgesOnEdge, nEdgesOnEdge, edgesOnEdge, weightsOnEdge,
    cellsOnEdge, verticesOnEdge, areaCell, dcEdge, dvEdge, fVertex, fEdge, radius, latCell,
    maxEdges, nEdgesOnCell, edgesOnCell, cellsOnCell,
  } = mesh;
  if (climatology && typeof climatology.columnAt !== 'function') throw new Error('the ocean climatology must be decoded (loadClimatology in ./climatology.module.js) before the ocean is built');
  const L = densities.length + 1, g = gravity, rho0 = density, rhoCp = density * specificHeat;
  const rho = [rho0, ...densities];
  const labelS = [referenceS, ...salinities];
  const labelT = rho.map((r, k) => Math.max(FREEZING_POINT, labelTemperature(r, labelS[k])));
  const thermoclineLayers = rho.filter((r, k) => k > 0 && r < THERMOCLINE_DENSITY).length;
  const diffusion = diffusivity * radius * radius / rhoCp;
  let spacing = 0, minSpacing = Infinity;
  for (let e = 0; e < E; e++) { spacing += dcEdge[e]; minSpacing = Math.min(minSpacing, dcEdge[e]); }
  spacing /= E;
  const nu4 = closureCoefficient(spacing, closureHours, closureSpacing);
  const eddyKappa = eddyDiffusivities(mesh, eddyDiffusivity), eddyLimit = eddyDiffusivity > 0 ? eddyDiffusionLimit(mesh) : 0;
  let eddy = null;

  const edgeOcean = geography ? geography.edgeOcean : new Uint8Array(E).fill(1);
  const cellOcean = geography ? Uint8Array.from(geography.land, (l) => (l ? 0 : 1)) : new Uint8Array(C).fill(1);
  const outlet = runoffOutlets(mesh, geography);
  const D = bathymetry ? Float64Array.from(bathymetry, (d, i) => (cellOcean[i] ? d : 0)) : bathymetryFrom(mesh, geography, { minimumDepth, flatDepth });
  let deepest = 0;
  for (let i = 0; i < C; i++) deepest = Math.max(deepest, D[i]);
  const substepLimit = 0.35 * minSpacing / Math.sqrt(g * Math.max(deepest, 1));

  const adopting = !!(buffers && buffers.h);
  const shared = (name, n) => new Float64Array(buffers && buffers[name] ? buffers[name] : new SharedArrayBuffer(8 * n));
  const h = shared('h', L * C), u = shared('u', L * E), Q = shared('Q', L * C), W = shared('W', L * C);
  const state = [h, u, Q, W];
  const eta = shared('eta', C);
  const T0 = new Float64Array(C), S0 = new Float64Array(C), rhoMl = shared('rhoMl', C), previousT0 = new Float64Array(C), surfaceIn = new Float64Array(C), previousIce = new Float64Array(C);
  const capacity = shared('capacity', C);
  if (!adopting) capacity.fill(rhoCp * mixedDepth);
  const stress = shared('stress', E), fresh = new Float64Array(C);
  const iced = new Uint8Array(C);
  const hEdge = shared('hEdge', L * E), pressure = shared('pressure', L * C), flux = new Float64Array(E), fluxPV = new Float64Array(E), tracerFlux = new Float64Array(E);
  const T = new Float64Array(C), S = new Float64Array(C), lapT = new Float64Array(C);
  const zeta = new Float64Array(V), qEdge = new Float64Array(E);
  const K = new Float64Array(C), phi = new Float64Array(C), gradPhi = new Float64Array(E), gradEta = shared('gradEta', E), gradRho = shared('gradRho', E);
  const lap = new Float64Array(E), lap2 = new Float64Array(E), divScratch = new Float64Array(C), curlScratch = new Float64Array(V);
  const slow = new Float64Array(E), U = new Float64Array(E), depthEdge = new Float64Array(E), etaB = new Float64Array(C), avgU = new Float64Array(E), avgEta = new Float64Array(C), divU = new Float64Array(C);
  const stages = [0, 1, 2, 3].map((s) => [shared(`stage${s}h`, L * C), shared(`stage${s}u`, L * E), shared(`stage${s}Q`, L * C), shared(`stage${s}W`, L * C)]);
  const trial = [shared('trialh', L * C), shared('trialu', L * E), shared('trialQ', L * C), shared('trialW', L * C)];
  const params = shared('params', 4);
  if (!adopting) params[0] = 1 / 3600;
  const tauCell = new Float64Array(3 * C), reach = new Float64Array(C).fill(Infinity), buoyancyLoss = new Float64Array(C);
  const lightest = rho.map((r, k) => r - 0.5 * (k > 1 ? r - rho[k - 1] : rho[k + 1] - r));
  const densest = rho.map((r, k) => r + 0.5 * (k < L - 1 ? rho[k + 1] - r : r - rho[k - 1]));
  let counter = 0, limited = 0, initialised = false, layerRunner = null;
  const relaxRate = () => params[0];

  const eos = seawaterDensity;
  const at = (k, i) => k * C + i;
  const ae = (k, e) => k * E + e;

  function maskedLaplacian(field, out) {
    for (let i = 0; i < C; i++) {
      let sum = 0;
      for (let m = 0; m < nEdgesOnCell[i]; m++) {
        const e = edgesOnCell[maxEdges * i + m];
        if (edgeOcean[e]) sum += dvEdge[e] * (field[cellsOnCell[maxEdges * i + m]] - field[i]) / dcEdge[e];
      }
      out[i] = sum / areaCell[i];
    }
  }
  /*
   * The mixed layer, present everywhere, takes the centred thickness at
   * an edge; an interior layer takes the smaller of the two, so it never
   * flows into a cell where it has no water, whether the layer has
   * outcropped there or the bottom lies above it. A column can only
   * flow through the water that exists on both sides of an edge, so at a
   * shelf break or a coast the thicknesses are scaled to the shallower
   * side's depth.
   */
  function edgeThicknesses(hIn) {
    for (let e = 0; e < E; e++) {
      const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
      hEdge[e] = 0.5 * (hIn[a] + hIn[b]);
      for (let k = 1; k < L; k++) hEdge[ae(k, e)] = Math.min(hIn[at(k, a)], hIn[at(k, b)]);
    }
    for (let e = 0; e < E; e++) {
      const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
      const sill = Math.max(EPS, Math.min(D[a], D[b]) + 0.5 * (eta[a] + eta[b]));
      let sum = 0;
      for (let k = 0; k < L; k++) sum += hEdge[ae(k, e)];
      if (sum > sill) { const f = sill / sum; for (let k = 0; k < L; k++) hEdge[ae(k, e)] *= f; }
    }
  }
  function coriolisTransport(Ue, out) {
    for (let e = 0; e < E; e++) {
      let sum = 0;
      for (let s = 0; s < nEdgesOnEdge[e]; s++) {
        const other = edgesOnEdge[maxEdgesOnEdge * e + s];
        sum += weightsOnEdge[maxEdgesOnEdge * e + s] * dvEdge[other] * Ue[other] * 0.5 * (fEdge[e] + fEdge[other]);
      }
      out[e] = sum / dcEdge[e];
    }
  }

  function surfaceDensity(hIn, QIn, WIn) {
    for (let i = 0; i < C; i++) {
      const h0 = Math.max(EPS, hIn[i]);
      rhoMl[i] = eos(QIn[i] / h0, WIn[i] / h0);
    }
  }

  /*
   * One tendency evaluation. Interior layer potentials telescope through
   * the layers above, so a token-thickness layer contributes nothing and
   * the layers it separates feel the density step across it.
   */
  function freeSurface() {
    for (let i = 0; i < C; i++) {
      let sum = 0;
      for (let k = 0; k < L; k++) sum += h[at(k, i)];
      eta[i] = cellOcean[i] ? sum - D[i] : 0;
    }
  }

  /*
   * The 2D fields every layer's tendency reads: the surface density and
   * the gradients of the free surface and of that density, the edge
   * thicknesses, and each interior layer's pressure potential
   * (ρ_ml − ρ_k) h₀ + Σ_{j<k} (ρ_j − ρ_k) h_j as two running sums down
   * the column.
   */
  function prepare(input) {
    const [hIn, , QIn, WIn] = input;
    surfaceDensity(hIn, QIn, WIn);
    gradient(mesh, eta, gradEta);
    gradient(mesh, rhoMl, gradRho);
    edgeThicknesses(hIn);
    for (let i = 0; i < C; i++) {
      let weighted = 0, total = 0;
      for (let k = 1; k < L; k++) {
        pressure[at(k, i)] = (rhoMl[i] - rho[k]) * hIn[i] + weighted - rho[k] * total;
        weighted += rho[k] * hIn[at(k, i)]; total += hIn[at(k, i)];
      }
    }
  }

  function tendency(input, out) {
    prepare(input);
    if (layerRunner) layerRunner(input === trial ? 1 : 0, stages.indexOf(out));
    else tendencyLayers(input, out, 0, L);
  }

  /*
   * The tendencies of layers kFrom to kTo from the prepared 2D fields —
   * the free surface and its gradient, the surface density gradient and
   * the edge thicknesses of every layer — reading the input state alone,
   * so that layers can be computed by any thread in any order; `part`
   * 'tracers' writes only dh, dQ and dW and 'momentum' only du.
   */
  function tendencyLayers(input, out, kFrom, kTo, part = 'all') {
    const [hIn, uIn, QIn, WIn] = input;
    const [dh, du, dQ, dW] = out;
    const relax = relaxRate(), tracers = part !== 'momentum', momentum = part !== 'tracers';
    for (let k = kFrom; k < kTo; k++) {
      const oc = k * C, oe = k * E;
      for (let e = 0; e < E; e++) {
        let he = hEdge[oe + e];
        fluxPV[e] = edgeOcean[e] ? he * uIn[oe + e] : 0;
        if (k === 0) he = Math.min(he, Math.max(0, hIn[cellsOnEdge[2 * e + (uIn[e] > 0 ? 0 : 1)]]));
        flux[e] = edgeOcean[e] ? he * uIn[oe + e] : 0;
      }
      if (tracers) {
      divergence(mesh, flux, divScratch);
      for (let i = 0; i < C; i++) {
        dh[oc + i] = -divScratch[i];
        const hv = hIn[oc + i];
        if (hv > 1e-6) { T[i] = QIn[oc + i] / hv; S[i] = WIn[oc + i] / hv; } else { T[i] = labelT[k]; S[i] = labelS[k]; }
      }
      for (let e = 0; e < E; e++) tracerFlux[e] = flux[e] * T[cellsOnEdge[2 * e + (flux[e] > 0 ? 0 : 1)]];
      divergence(mesh, tracerFlux, divScratch);
      for (let i = 0; i < C; i++) dQ[oc + i] = -divScratch[i];
      for (let e = 0; e < E; e++) tracerFlux[e] = flux[e] * S[cellsOnEdge[2 * e + (flux[e] > 0 ? 0 : 1)]];
      divergence(mesh, tracerFlux, divScratch);
      for (let i = 0; i < C; i++) dW[oc + i] = -divScratch[i];
      if (k === 0 && diffusion > 0) {
        maskedLaplacian(T, lapT);
        for (let i = 0; i < C; i++) dQ[i] += diffusion * lapT[i];
        maskedLaplacian(S, lapT);
        for (let i = 0; i < C; i++) dW[i] += diffusion * lapT[i];
      }
      for (let i = 0; i < C; i++) if (!cellOcean[i]) { dh[oc + i] = 0; dQ[oc + i] = 0; dW[oc + i] = 0; }
      }
      if (!momentum) continue;
      curl(mesh, uIn.subarray(oe, oe + E), zeta);
      for (let e = 0; e < E; e++) qEdge[e] = 0.5 * (zeta[verticesOnEdge[2 * e]] + fVertex[verticesOnEdge[2 * e]] + zeta[verticesOnEdge[2 * e + 1]] + fVertex[verticesOnEdge[2 * e + 1]]) / Math.max(hEdge[oe + e], k > 0 ? vorticityCentring * 0.5 * (hIn[oc + cellsOnEdge[2 * e]] + hIn[oc + cellsOnEdge[2 * e + 1]]) : 0, PV_FLOOR);
      kineticEnergy(mesh, uIn.subarray(oe, oe + E), K);
      if (k === 0) {
        for (let i = 0; i < C; i++) phi[i] = K[i] + g * eta[i];
      } else {
        for (let i = 0; i < C; i++) phi[i] = K[i] + g * eta[i] + g * pressure[oc + i] / rho0;
      }
      gradient(mesh, phi, gradPhi);
      for (let e = 0; e < E; e++) {
        let sum = 0;
        for (let s = 0; s < nEdgesOnEdge[e]; s++) {
          const other = edgesOnEdge[maxEdgesOnEdge * e + s];
          sum += weightsOnEdge[maxEdgesOnEdge * e + s] * dvEdge[other] * fluxPV[other] * 0.5 * (qEdge[e] + qEdge[other]);
        }
        du[oe + e] = sum / dcEdge[e] - gradPhi[e];
      }
      if (k === 0) for (let e = 0; e < E; e++) du[e] -= g / rho0 * 0.5 * hEdge[e] * gradRho[e];
      for (let e = 0; e < E; e++) {
        const he = Math.max(hEdge[oe + e], minimumThickness);
        let force = 0, drag = 0;
        if (k === 0) force += stress[e] / rho0;
        if (k > 0) {
          let j = k - 1;
          while (j > 0 && hEdge[ae(j, e)] < THIN) j--;
          drag += interfacialDrag * (uIn[ae(j, e)] - uIn[oe + e]);
        }
        if (k < L - 1) {
          let j = k + 1;
          while (j < L - 1 && hEdge[ae(j, e)] < THIN) j++;
          if (hEdge[ae(j, e)] >= THIN) drag -= interfacialDrag * (uIn[oe + e] - uIn[ae(j, e)]);
        }
        let bottom = k === L - 1;
        if (!bottom) { bottom = true; for (let j = k + 1; j < L; j++) if (hEdge[ae(j, e)] >= THIN) { bottom = false; break; } }
        if (bottom) force -= bottomDrag * Math.abs(uIn[oe + e]) * uIn[oe + e];
        du[oe + e] += force / he + drag / (k === 0 ? he : Math.max(hEdge[oe + e], THIN));
      }
      if (nu4 > 0) {
        laplacianVelocity(mesh, uIn.subarray(oe, oe + E), lap, divScratch, curlScratch);
        laplacianVelocity(mesh, lap, lap2, divScratch, curlScratch);
        for (let e = 0; e < E; e++) du[oe + e] -= nu4 * lap2[e];
      }
      if (k > 0) for (let e = 0; e < E; e++) if (hEdge[oe + e] < THIN) du[oe + e] = (uIn[ae(k - 1, e)] - uIn[oe + e]) * relax;
      for (let e = 0; e < E; e++) if (!edgeOcean[e]) du[oe + e] = 0;
    }
  }

  function combine(target, base, stage, dt) {
    for (let a = 0; a < state.length; a++) {
      const t = target[a], b = base[a], s = stage[a];
      for (let i = 0; i < t.length; i++) t[i] = b[i] + dt * s[i];
    }
  }

  function barotropic(dt) {
    const [k1] = stages;
    const [, du] = k1;
    for (let e = 0; e < E; e++) {
      let sumH = 0, transport = 0, forcing = 0;
      for (let k = 0; k < L; k++) { const he = hEdge[ae(k, e)]; sumH += he; transport += he * u[ae(k, e)]; forcing += he * (du[ae(k, e)] + g * gradEta[e]); }
      U[e] = edgeOcean[e] ? transport : 0;
      slow[e] = edgeOcean[e] ? forcing : 0;
      depthEdge[e] = sumH;
    }
    coriolisTransport(U, lap);
    for (let e = 0; e < E; e++) slow[e] -= lap[e];
    const M = Math.max(1, Math.ceil(dt / substepLimit)), dtb = dt / M;
    etaB.set(eta); avgU.fill(0); avgEta.fill(0);
    for (let m = 0; m < M; m++) {
      rk4Barotropic(dtb);
      for (let i = 0; i < C; i++) avgEta[i] += etaB[i] / M;
      for (let e = 0; e < E; e++) avgU[e] += U[e] / M;
    }
  }

  const bStages = [0, 1, 2, 3].map(() => [new Float64Array(C), new Float64Array(E)]);
  const bTrial = [new Float64Array(C), new Float64Array(E)];
  function barotropicTendency(etaIn, uIn, dEta, dU) {
    gradient(mesh, etaIn, gradPhi);
    coriolisTransport(uIn, lap);
    for (let e = 0; e < E; e++) {
      if (!edgeOcean[e]) { dU[e] = 0; continue; }
      dU[e] = -g * depthEdge[e] * gradPhi[e] + lap[e] + slow[e];
    }
    divergence(mesh, uIn, divU);
    for (let i = 0; i < C; i++) dEta[i] = cellOcean[i] ? -divU[i] : 0;
  }
  function rk4Barotropic(dtb) {
    const [s1, s2, s3, s4] = bStages, [tEta, tU] = bTrial;
    barotropicTendency(etaB, U, s1[0], s1[1]);
    for (let i = 0; i < C; i++) tEta[i] = etaB[i] + 0.5 * dtb * s1[0][i];
    for (let e = 0; e < E; e++) tU[e] = U[e] + 0.5 * dtb * s1[1][e];
    barotropicTendency(tEta, tU, s2[0], s2[1]);
    for (let i = 0; i < C; i++) tEta[i] = etaB[i] + 0.5 * dtb * s2[0][i];
    for (let e = 0; e < E; e++) tU[e] = U[e] + 0.5 * dtb * s2[1][e];
    barotropicTendency(tEta, tU, s3[0], s3[1]);
    for (let i = 0; i < C; i++) tEta[i] = etaB[i] + dtb * s3[0][i];
    for (let e = 0; e < E; e++) tU[e] = U[e] + dtb * s3[1][e];
    barotropicTendency(tEta, tU, s4[0], s4[1]);
    for (let i = 0; i < C; i++) etaB[i] += dtb / 6 * (s1[0][i] + 2 * s2[0][i] + 2 * s3[0][i] + s4[0][i]);
    for (let e = 0; e < E; e++) U[e] += dtb / 6 * (s1[1][e] + 2 * s2[1][e] + 2 * s3[1][e] + s4[1][e]);
  }

  function step(dt) {
    const [k1, k2, k3, k4] = stages;
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
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) continue;
      let sum = 0, dh = 0, dQ = 0, dW = 0;
      for (let k = 1; k < L; k++) {
        const n = at(k, i);
        if (h[n] < EPS) {
          const held = h[n] > 1e-9, t = labelT[k], s = labelS[k];
          const tHeld = held ? Math.max(t - 30, Math.min(t + 30, Q[n] / h[n])) : t, sHeld = held ? Math.max(s - 5, Math.min(s + 5, W[n] / h[n])) : s;
          dh += EPS - h[n]; dQ += EPS * t - h[n] * tHeld; dW += EPS * s - h[n] * sHeld;
          h[n] = EPS; Q[n] = EPS * t; W[n] = EPS * s;
        }
        sum += h[n];
      }
      h[i] -= dh; Q[i] -= dQ; W[i] -= dW;
      sum += h[i];
      const scale = (D[i] + avgEta[i]) / sum;
      for (let k = 0; k < L; k++) { const n = at(k, i); h[n] *= scale; Q[n] *= scale; W[n] *= scale; }
      eta[i] = avgEta[i];
    }
    if (eddyLimit > 0) eddyTransport(dt);
    edgeThicknesses(h);
    limited = 0;
    for (let e = 0; e < E; e++) {
      if (!edgeOcean[e]) continue;
      let sumH = 0, transport = 0, clamped = false;
      for (let k = 0; k < L; k++) { sumH += hEdge[ae(k, e)]; transport += hEdge[ae(k, e)] * u[ae(k, e)]; }
      const shift = (avgU[e] - transport) / Math.max(sumH, EPS);
      for (let k = 0; k < L; k++) { const n = ae(k, e); u[n] += shift; if (Math.abs(u[n]) > SPEED_LIMIT) { u[n] = Math.sign(u[n]) * SPEED_LIMIT; clamped = true; } }
      if (clamped) limited++;
      for (let k = 1; k < L; k++) if (hEdge[ae(k, e)] < THIN) u[ae(k, e)] = u[ae(k - 1, e)];
    }
  }

  function eddyTransport(dt) {
    if (!eddy) eddy = { volume: new Float64Array(L * C), heat: new Float64Array(L * C), salt: new Float64Array(L * C), flux: new Float64Array(L), share: new Float64Array(L) };
    const { volume, heat, flux: classFlux, share } = eddy, salinity = eddy.salt;
    const kappaLimit = eddyLimit / dt, ramp = (x, width) => Math.max(0, Math.min(1, x / width));
    const allowance = (k, i) => (Math.max(0, h[at(k, i)] - EPS) + EDDY_SLACK) * areaCell[i] / (nEdgesOnCell[i] * dt);
    volume.fill(0); heat.fill(0); salinity.fill(0);
    for (let e = 0; e < E; e++) {
      const kappa = Math.min(eddyKappa[e], kappaLimit);
      if (!edgeOcean[e] || !(kappa > 0)) continue;
      const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
      let shared = 0;
      for (let k = 1; k < L; k++) shared += Math.max(0, Math.min(h[at(k, a)], h[at(k, b)]) - EPS);
      const sill = Math.min(D[a], D[b]) + 0.5 * (eta[a] + eta[b]), coefficient = kappa * dvEdge[e] / dcEdge[e];
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
      for (let k = 1; k < L; k++) {
        const f = scale * classFlux[k];
        if (f === 0) continue;
        const d = at(k, f > 0 ? a : b), na = at(k, a), nb = at(k, b), fq = f * Q[d] / h[d], fw = f * W[d] / h[d];
        volume[na] -= f; volume[nb] += f;
        heat[na] -= fq; heat[nb] += fq;
        salinity[na] -= fw; salinity[nb] += fw;
      }
    }
    for (let n = C; n < L * C; n++) {
      const factor = dt / areaCell[n % C];
      h[n] += factor * volume[n]; Q[n] += factor * heat[n]; W[n] += factor * salinity[n];
    }
  }

  function move(i, from, to, amount) {
    const a = at(from, i), b = at(to, i);
    const f = amount / h[a];
    const dQ = Q[a] * f, dW = W[a] * f;
    h[a] -= amount; Q[a] -= dQ; W[a] -= dW;
    h[b] += amount; Q[b] += dQ; W[b] += dW;
  }

  function detrain(i, amount, rm) {
    if (amount <= 0) return;
    let k = 1;
    for (let j = 2; j < L; j++) if (Math.abs(rho[j] - rm) < Math.abs(rho[k] - rm)) k = j;
    move(i, 0, k, amount);
  }

  function neighbourReach() {
    for (let i = 0; i < C; i++) {
      let sum = 0, count = 0;
      for (let m = 0; m < nEdgesOnCell[i]; m++) {
        const j = cellsOnCell[maxEdges * i + m];
        if (cellOcean[j]) { sum += h[j]; count++; }
      }
      reach[i] = count ? mixedNeighbourRatio * sum / count : Infinity;
    }
  }

  function mixedLayer(dt, ice) {
    cellVector(mesh, stress, tauCell);
    if (mixedNeighbourRatio > 0) neighbourReach();
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) continue;
      const tau = Math.hypot(tauCell[3 * i], tauCell[3 * i + 1], tauCell[3 * i + 2]);
      const ustar3 = Math.pow(tau / rho0, 1.5), stir = stirring * Math.exp(-h[i] / stirringDepth);
      const deepest = Math.min(maximumMixedDepth, reach[i]), s0 = W[i] / h[i];
      const salted = s0 * fresh[i] / 1000 + (s0 - iceSalinity) * (ice[i] - previousIce[i]) * iceDensity / 1000;
      const buoyancy = g * (expansionOf(surfaceIn[i], s0) * (previousT0[i] - surfaceIn[i]) * capacity[i] / rhoCp + contractionOf(surfaceIn[i], s0) * salted) / dt;
      buoyancyLoss[i] += (buoyancy - buoyancyLoss[i]) * (buoyancyMemory > 0 ? Math.min(1, dt / buoyancyMemory) : 1);
      let rm = eos(Q[i] / h[i], W[i] / h[i]), budget = convectiveRate * dt;
      for (let k = 1; k < L && h[i] < deepest && budget > 0; k++) {
        const n = at(k, i);
        if (h[n] <= EPS) continue;
        let take = Math.min(h[n] - EPS, deepest - h[i], budget);
        const eroding = convectiveErosion && rm < densest[k];
        if (eroding) take = rm >= lightest[k] && buoyancyLoss[i] > 0 ? Math.min(take, buoyancyLoss[i] * dt * rho0 * h[n] / (g * h[i] * (densest[k] - rm))) : 0;
        else if (!convectiveErosion && rho[k] > rm) continue;
        if (take > 0) { move(i, k, 0, take); budget -= take; rm = eos(Q[i] / h[i], W[i] / h[i]); }
        if (eroding) break;
      }
      let below = -1;
      for (let k = 1; k < L; k++) if (h[at(k, i)] > THIN) { below = k; break; }
      if (below > 0 && h[i] < deepest) {
        const db = Math.max(1e-3, g * (rho[below] - rm) / rho0);
        const entrain = Math.min(2 * stir * ustar3 / (h[i] * db) * dt, h[at(below, i)] - EPS, reach[i] - h[i]);
        if (entrain > 0) { move(i, below, 0, entrain); rm = eos(Q[i] / h[i], W[i] / h[i]); }
      }
      below = -1;
      for (let k = 1; k < L; k++) if (h[at(k, i)] > THIN) { below = k; break; }
      let excess = Math.max(0, h[i] - maximumMixedDepth);
      if (neutralSnap && below > 0 && rm >= rho[below] - DENSITY_TOLERANCE) excess = Math.max(excess, h[i] - shallowestMixedDepth);
      detrain(i, excess, rm);
      if (h[i] > reach[i]) detrain(i, (h[i] - Math.max(reach[i], shallowestMixedDepth)) * Math.min(1, dt / detrainmentTime), rm);
      if (buoyancyLoss[i] < -1e-9) {
        const monin = Math.max(shallowestMixedDepth, 2 * stir * ustar3 / -buoyancyLoss[i]);
        if (h[i] > monin) detrain(i, (h[i] - monin) * Math.min(1, dt / detrainmentTime), rm);
      }
      if (h[i] < minimumThickness) {
        for (let k = 1; k < L && h[i] < minimumThickness; k++) {
          const available = h[at(k, i)] - EPS;
          if (available > 0) move(i, k, 0, Math.min(available, minimumThickness - h[i]));
        }
      }
      restoreDensities(i, Math.min(1, dt / restoreTime));
    }
  }

  function restoreDensities(i, fraction) {
    for (let k = 1; k < L; k++) {
      const n = at(k, i), target = rho[k];
      if (h[n] <= THIN) continue;
      const r = eos(Q[n] / h[n], W[n] / h[n]);
      if (Math.abs(r - target) <= RESTORE_TOLERANCE) continue;
      const step = r < target ? 1 : -1;
      for (let j = k + step; j >= 1 && j < L; j += step) {
        const d = at(j, i);
        if (h[d] <= THIN) continue;
        const rd = eos(Q[d] / h[d], W[d] / h[d]);
        if (step * (rd - target) <= RESTORE_TOLERANCE) continue;
        const amount = Math.min(fraction * h[n] * (target - r) / (rd - target), h[d] - EPS);
        if (amount > 0) move(i, j, k, amount);
        break;
      }
    }
  }

  function setStress(total, ice, concentration = null) {
    const cover = (i) => (ice[i] > 0 ? (concentration && concentration[i] > 0 ? concentration[i] : 1) : 0);
    for (let e = 0; e < E; e++) stress[e] = edgeOcean[e] ? total[e] * (1 - 0.5 * (cover(cellsOnEdge[2 * e]) + cover(cellsOnEdge[2 * e + 1])) * (1 - iceStressTransmission)) : 0;
  }
  function readSurface(surfaceT, ice) {
    for (let i = 0; i < C; i++) {
      iced[i] = ice[i] > 0 ? 1 : 0;
      T0[i] = iced[i] ? FREEZING_POINT : surfaceT[i];
      surfaceIn[i] = T0[i];
      Q[i] = h[i] * T0[i];
    }
  }
  function accumulate(evaporation, rain, dt, runoff = null) {
    for (let i = 0; i < C; i++) if (cellOcean[i]) fresh[i] += (evaporation ? evaporation[i] * dt : 0) - (rain ? rain[i] : 0);
    if (runoff) for (let j = 0; j < C; j++) if (!cellOcean[j] && runoff[j] > 0 && outlet[j] >= 0) fresh[outlet[j]] -= runoff[j] * areaCell[j] / areaCell[outlet[j]];
  }
  function salt(dt, ice) {
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) continue;
      const s = W[i] / h[i];
      W[i] += s * fresh[i] / 1000;
      const grown = (ice[i] - previousIce[i]) * iceDensity / 1000;
      W[i] += (s - iceSalinity) * grown;
      W[i] = Math.max(0, W[i]);
      fresh[i] = 0;
      previousIce[i] = ice[i];
    }
  }
  function writeSurface(surfaceT, oceanFlux, dt) {
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) { oceanFlux[i] = 0; continue; }
      capacity[i] = rhoCp * Math.max(h[i], 1);
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

  function advance(surfaceT, ice, oceanFlux, totalStress, dt, concentration = null) {
    if (++counter % everySteps !== 0) return false;
    if (!initialised) initialize(surfaceT, ice);
    const dtOcean = everySteps * dt;
    params[0] = Math.min(1 / 3600, 1 / dtOcean);
    readSurface(surfaceT, ice);
    setStress(typeof totalStress === 'function' ? totalStress() : totalStress, ice, concentration);
    step(dtOcean);
    mixedLayer(dtOcean, ice);
    salt(dtOcean, ice);
    writeSurface(surfaceT, oceanFlux, dtOcean);
    return true;
  }

  /*
   * The start from rest. With a climatology (`climatology`, decoded by
   * ./climatology.module.js), every sea cell it covers takes its column
   * (atlasColumns), and `initialize` writes the mixed layer's
   * temperature into surfaceT over open sea and returns how many sea
   * cells took the atlas and how many did not. The sea cells it does not
   * cover, and all of them without one, start from the analytic
   * climatology: a mixed layer over interior layers whose bases sit at
   * their subtropical depths, shallower toward the poles and rising to
   * the surface over one layer's density step as the surface water
   * approaches each layer's density, the first whose base lies below the
   * sea floor filling to it and the deepest layer holding water giving up
   * the metres of the tokens beneath. The free surface starts at its
   * steric height, so the
   * pressure below the thermocline is level and the ocean does not begin
   * with a barotropic shock.
   */
  function initialize(surfaceT, ice, { climatology: atlas = climatology } = {}) {
    const filled = build(surfaceT, ice, atlas);
    if (!filled) return null;
    let covered = 0, other = 0;
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) continue;
      if (filled[i]) covered++; else other++;
      if (filled[i] && !iced[i]) surfaceT[i] = T0[i];
    }
    return { atlas: covered, analytic: other };
  }
  function build(surfaceT, ice, atlas) {
    u.fill(0); eta.fill(0); fresh.fill(0); buoyancyLoss.fill(0);
    const filled = atlas ? atlasColumns(mesh, atlas, { D, cellOcean, ice, rho, labelT, labelS, h, Q, W, T0, shallowestMixedDepth, maximumMixedDepth }) : null;
    for (let i = 0; i < C; i++) {
      previousIce[i] = ice[i];
      iced[i] = ice[i] > 0 ? 1 : 0;
      const lat = latCell[i];
      if (filled && filled[i]) {
        S0[i] = W[i] / h[i];
        previousT0[i] = T0[i];
        capacity[i] = rhoCp * Math.max(h[i], 1);
        continue;
      }
      const s0 = salinityProfile(lat);
      T0[i] = iced[i] ? FREEZING_POINT : surfaceT[i];
      S0[i] = s0;
      if (!cellOcean[i]) { for (let k = 0; k < L; k++) { h[at(k, i)] = 0; Q[at(k, i)] = 0; W[at(k, i)] = 0; } continue; }
      const rm = eos(T0[i], s0);
      const stretch = 1 - thermoclineTilt + 2 * thermoclineTilt * Math.cos(lat) ** 2;
      let cumulative = Math.min(mixedDepth, D[i]);
      h[i] = cumulative; Q[i] = cumulative * T0[i]; W[i] = cumulative * s0;
      for (let k = 1; k < L; k++) {
        let hk;
        if (k === L - 1) hk = Math.max(0, D[i] - cumulative);
        else {
          const taper = Math.max(0, Math.min(1, (rho[k] - rm) / (rho[k + 1] - rho[k])));
          hk = Math.max(0, Math.min(D[i], bottoms[k - 1] * stretch * taper) - cumulative);
        }
        hk = Math.max(EPS, hk);
        cumulative += hk;
        const [tk, sk] = hk > EPS ? interiorWater(rho[k], labelT[k], labelS[k], lat) : [labelT[k], labelS[k]];
        h[at(k, i)] = hk; Q[at(k, i)] = hk * tk; W[at(k, i)] = hk * sk;
      }
      const excess = cumulative - D[i];
      let giver = L - 1;
      while (giver > 0 && !(h[at(giver, i)] > excess + EPS)) giver--;
      const n = at(giver, i), f = (h[n] - excess) / h[n];
      h[n] -= excess; Q[n] *= f; W[n] *= f;
      previousT0[i] = T0[i];
      capacity[i] = rhoCp * Math.max(h[i], 1);
    }
    surfaceDensity(h, Q, W);
    stericSurface();
    counter = 0;
    initialised = true;
    return filled;
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
        anomaly += ((k === 0 ? rhoMl[i] : rho[k]) - rho0) * part;
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
   * A saved ocean on other classes is carried onto these (rebinOcean)
   * before its columns are fitted to the bathymetry.
   */
  function load(saved, surfaceT, ice) {
    const classes = savedDensities(saved, C);
    if (classes && !sameDensities(classes, densities)) saved = rebinOcean(saved, classes, densities, C, { labelT, labelS });
    if (!saved.h || saved.h.length !== L * C) {
      build(surfaceT, ice, climatology);
      return;
    }
    h.set(saved.h); u.set(saved.u); eta.set(saved.eta);
    for (let n = 0; n < L * C; n++) { Q[n] = h[n] * saved.T[n]; W[n] = h[n] * saved.S[n]; }
    const start = () => {
      const kept = { h: Float64Array.from(h), u: Float64Array.from(u), Q: Float64Array.from(Q), W: Float64Array.from(W), eta: Float64Array.from(eta) };
      build(surfaceT, ice, climatology);
      const built = { h: Float64Array.from(h), Q: Float64Array.from(Q), W: Float64Array.from(W) };
      h.set(kept.h); u.set(kept.u); Q.set(kept.Q); W.set(kept.W); eta.set(kept.eta);
      return built;
    };
    fitColumns({ h, Q, W, eta }, start, { D, cellOcean, L, C, labelT, labelS, minimumThickness });
    for (let k = 0; k < L; k++) for (let e = 0; e < E; e++) if (!edgeOcean[e]) u[ae(k, e)] = 0;
    for (let i = 0; i < C; i++) { previousIce[i] = ice[i]; iced[i] = ice[i] > 0 ? 1 : 0; T0[i] = Q[i] / Math.max(EPS, h[i]); S0[i] = W[i] / Math.max(EPS, h[i]); previousT0[i] = T0[i]; capacity[i] = rhoCp * Math.max(h[i], 1); }
    fresh.fill(0); buoyancyLoss.fill(0);
    counter = 0;
    initialised = true;
  }

  function serialize() {
    const Tall = new Array(L * C), Sall = new Array(L * C);
    for (let n = 0; n < L * C; n++) { const hh = Math.max(EPS, h[n]); Tall[n] = Q[n] / hh; Sall[n] = W[n] / hh; }
    return { h: Array.from(h), u: Array.from(u), T: Tall, S: Sall, eta: Array.from(eta), densities: Array.from(densities) };
  }

  const thermoclineDepth = new Float64Array(C), sst = T0, sss = S0;
  function thermocline(i) {
    let depth = 0;
    for (let k = 0; k <= thermoclineLayers; k++) depth += h[at(k, i)];
    return depth;
  }
  function fields(depth = 0) {
    for (let i = 0; i < C; i++) thermoclineDepth[i] = cellOcean[i] ? thermocline(i) : NaN;
    const layerT = (k, i) => (k === 0 ? T0[i] : Q[at(k, i)] / h[at(k, i)]);
    return { h1: h.subarray(0, C), T1: T0, S1: S0, u1: u.subarray(0, E), eta, thermoclineDepth, ...depthFields(mesh, L, { h, u, temperature: layerT, cellOcean }, depth) };
  }

  function diagnostics() {
    let area = 0, depth = 0, heat = 0, thermo = 0, salinity = 0, ssh = 0, speed = 0, transport = 0, interiorT = 0, interiorH = 0;
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) continue;
      const a = areaCell[i];
      area += a; depth += a * h[i]; salinity += a * (W[i] / h[i]); ssh = Math.max(ssh, Math.abs(eta[i]));
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
    return { oceanUpperDepth: depth / area, oceanHeat: heat / area, oceanInteriorT: interiorH > 0 ? interiorT / interiorH : 0, oceanSpeed: speed, oceanThermoclineDepth: thermo / area, oceanSalinity: salinity / area, oceanSSH: ssh, oceanTransport: transport / 1e6, oceanLimited: limited };
  }

  if (!adopting) build(new Float64Array(C).fill(288), new Float64Array(C), null);
  initialised = false;
  const sharedBuffers = { h: h.buffer, u: u.buffer, Q: Q.buffer, W: W.buffer, eta: eta.buffer, rhoMl: rhoMl.buffer, capacity: capacity.buffer, stress: stress.buffer, hEdge: hEdge.buffer, pressure: pressure.buffer, gradEta: gradEta.buffer, gradRho: gradRho.buffer, params: params.buffer, trialh: trial[0].buffer, trialu: trial[1].buffer, trialQ: trial[2].buffer, trialW: trial[3].buffer };
  stages.forEach((stage, s) => { sharedBuffers[`stage${s}h`] = stage[0].buffer; sharedBuffers[`stage${s}u`] = stage[1].buffer; sharedBuffers[`stage${s}Q`] = stage[2].buffer; sharedBuffers[`stage${s}W`] = stage[3].buffer; });
  return { state, trial, stages, layers: L, h, u, Q, W, eta, D, T0, S0, capacity, stress, fresh, tendency, tendencyLayers, setLayerRunner(fn) { layerRunner = fn; }, advance, accumulate, initialize, load, serialize, diagnostics, fields, setStress, readSurface, eddyTransport, rhoCp, edgeOcean, cellOcean, densities: rho, shared: sharedBuffers };
}
