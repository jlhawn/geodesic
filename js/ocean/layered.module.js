import { divergence, gradient, curl, kineticEnergy, laplacianVelocity, cellVector } from '../dynamics/operators.module.js';
import { FREEZING_POINT } from '../physics/ice.module.js';
import { SEAWATER, seawaterDensity, thermalExpansion as expansionOf, halineContraction as contractionOf, labelTemperature, salinityForDensity } from './seawater.module.js';
import { profileAt } from './climatology.module.js';
import { FRICTION_PER_WIND, CALM_WIND, RETREAT_EKMAN, erfc, stokesSurfaceLayer, langmuirFactor, JERLOV, PENETRATION_DEPTH, transmitted, deposited, HYBRID_LAYERS, H_MIN, MODE_DEPTH, labelBoundaries, columnScratch, targetInterfaces, remapColumn, solveColumn, mixColumn, mixedLayerDepth, spaceInterfaces, layerDensities } from './column.module.js';

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
 * pressure gradient balances, a jet at meters a second. The potential
 * vorticity's thickness is no less than PV_FLOOR for a class and THIN for
 * the mixed layer.
 *
 * The wind stress τ is a body force over the top H = max(H_m, h₀) of each
 * edge's column, tokens counted for depth, on the mixed layer and the
 * classes at least `stressThickness` thick: a layer with o_k of its edge
 * thickness h_k within H gains (τ/ρ₀)(o_k/h_k)/Σo, so the column takes
 * τ/ρ₀ whatever the mixed layer's depth and a mixed layer thinner than H
 * runs at the speed of one H thick. H_m is momentumDepth or, with
 * `ekmanDepth` c_E, the wind's Ekman depth c_E u* / |f| from the edge's
 * normal and reconstructed tangential stress, held between
 * `momentumFloor` and momentumDepth. A thinner class, which feels only a
 * fraction of its Coriolis force, takes its share through the drag. Each
 * layer divides its interfacial drag by its own edge thickness, no less
 * than THIN, so the drag moves momentum between layers without making
 * any, and the bottom drag by no less than momentumDepth.
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
 * `stirringDepth`. With `shearEntrainment` a critical bulk Richardson
 * number Ri_c, it first takes in the class beneath to h₀√(Ri_c/Ri_b), at
 * most that class in a step, where Ri_b = Δb h₀/|Δu|² and Δu is the
 * cell's reconstruction of each edge's velocity less that of its first
 * class at least THIN thick (Price, Weller and Pinkel 1986); the exchange
 * of `carryMomentum`, which it needs, then brings Ri_b to Ri_c, and a
 * retreat in the same step stops at that depth. It deepens within
 * `maximumMixedDepth` and, when
 * `mixedNeighbourRatio` is set, within that multiple of the mean depth of
 * the neighboring mixed layers. Otherwise it holds its depth,
 * convectively neutral or not: it detrains at once what lies beyond the
 * maximum, and, over `detrainmentTime` and never above
 * `shallowestMixedDepth`, what lies beyond the Monin–Obukhov depth
 * C_r ε u*³/(−B) while it gains buoyancy, and what lies beyond its
 * neighbors' reach. C_r is the entrainment's 2m, or `retreatStirring`
 * without the fading with depth; ε is 1 or, with `langmuir` as its cap,
 * the Langmuir factor at the depth the step began with, its excess over 1
 * scaled by the cell's open water, as pack ice damps the waves. With
 * `retreatEkman` c that depth is at most c u* / |f|, the Ekman bound of
 * Large, McWilliams and Doney (1994), and unbounded where f is 0.
 * With `retreatDaily` the retreat alone takes its B from the previous
 * UTC day instead: the day's mean absorbed sunlight S̄ (the light handed
 * below the mixed layer over the share τ_k that passed it) and mean dark
 * loss D̄ (B with the sunlight the mixed layer kept, b S (1 − τ_k), put
 * back, b = gα/ρc_p), means over a whole day of ocean calls ending at a
 * UTC midnight (closeDay; a partial day after a start or a load is
 * dropped), so that B = D̄ − b S̄ (1 − τ(h₀)) at the depth h₀ the layer
 * has now. With `retreatLight` λ it retreats only while
 * D̄ − b S̄ (1 − τ(h₀) − λ h₀ q(h₀)) < 0, q the light absorbed per meter
 * at its base (deposited): it holds the depth to which the night's
 * convection carries the day's light, the zero-wind trapping depth.
 * `neutralSnap` returns everything below
 * `shallowestMixedDepth` whenever the layer is as dense as the water
 * beneath it, and `convectiveErosion: false` swallows every class no
 * denser than the mixed layer at `convectiveRate`.
 * With `carryMomentum`, momentum moves with the water: after the
 * exchanges each edge's velocities are made uniform down to the deeper of
 * its old and new mixed-layer base, and its column's Σ h_e u is then
 * restored by a uniform shift, so the water a retreat leaves keeps the
 * mixed layer's momentum and a refill or convection gives the mixed layer
 * the mean of what it now holds.
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
 * Tracers are carried by the flux with the donor cell's value (the hybrid's
 * with the limited third-order edge value of limitedEdgeValue); the mixed
 * layer's heat and salt also diffuse, as ∇·(h κ ∇T) with κ
 * `mixedDiffusivity` and the harmonic mean of the two cells' thickness on
 * each edge, or, with `diffusivity` alone, as a flux of heat per unit
 * temperature gradient that does not see the thickness. Surface
 * freshwater (evaporation minus rain, and runoff where it flows down to the sea) and ice
 * growth or melt, counted in the ice's volume per unit area of the cell
 * (its cover times its thickness), act on its salinity as virtual salt
 * fluxes. The mixed
 * layer's temperature is the sea surface temperature the atmosphere sees, its
 * heat capacity is published per cell, and the heat converged under ice
 * is handed to the ice base. With `penetration` a Jerlov water type, the
 * ocean also publishes `transmission`, the share of the absorbed sunlight
 * that passes the mixed layer (penetrationShares); the surface update
 * keeps that part of open water's sunlight out of the sea surface
 * temperature and adds it to `penetrated` (J/m²), which the next
 * readSurface deposits in the classes below by depth.
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
 * With `ekmanBlock`, after the eddy transport each edge carries the steady
 * Ekman transport of its own stress,
 *   U_E = (r τ_n + f τ_t) / (ρ₀ (r² + f²)),
 * r = interfacialDrag/H with H the reach of the stress (stressReach), as
 * the balance of the edge's Coriolis term against that drag. It crosses
 * with the composition of the donor cell's top B = max(wind depth, h₀)
 * meters, h₀ its mixed layer's thickness, each layer's part o_k in
 * proportion to its thickness there above the token thickness, the shares
 * ω_k = o_k/Σo, where the dynamics gave the edge's own layers the stress
 * shares σ_k. Layer k's flux changes by
 * U_E dv (ω_k − σ_k), which sums to zero over the column, so an edge's
 * column transport and the free surface are untouched, and the change
 * takes the heat and salt of its upwind cell, as the eddy transport's do.
 * One factor scales an edge's changes together so that no layer gives up
 * more than 1/nEdges of its water above the token thickness in the step.
 * Where each donor's mixed layer is as deep as the wind depth, ω = σ and
 * the pass moves nothing.
 *
 * With `upwellingSource` (which needs `ekmanBlock`), in a sea cell with a
 * land neighbor the mixed layer takes the water it entrains and refills
 * from the class that holds the base of the Ekman block, B = max(mean
 * wind depth of the cell's sea edges, h₀), the first class from there
 * down more than THIN thick, for as much as the Ekman divergence of the
 * cell, w = Σ s U_E dv/A for U_E the pass's transport and s the edge's
 * outward sign, times `upwellingSource` and the step, upwells; the rest
 * comes from where it would otherwise. The upwelling that a coast's
 * narrow strip carries is smaller than a cell, so the water rising there
 * is not the cell's top class but the one the block's base is in.
 */

/*
 * With `intakeLedger` (true for every sea cell, or a mask of the cells to
 * record) the ocean's `intake` holds, for every recorded cell, the water the
 * mixed layer takes from each class and the heat content it carries, by the
 * process that took it: the convective loop's erosion and swallowing, the
 * shear entrainment, the stirring entrainment, the refill to the minimum
 * thickness, and the draw from the upwelling source (`upwellingSource`),
 * which meets a share of the shear, stirring and refill demands. The water,
 * in meters, and its temperature times meters, in K m, for process p in
 * INTAKE_PROCESSES, class k and recorded cell i are at
 * (p L + k) cells.length + slot[i] of `water` and `degrees`, summed since the
 * arrays were last filled with zeros. The ledger changes no state.
 */
export const INTAKE_PROCESSES = ['convective', 'shear', 'stirring', 'refill', 'upwelled'];
const [INTAKE_CONVECTIVE, INTAKE_SHEAR, INTAKE_STIRRING, INTAKE_REFILL, INTAKE_UPWELLED] = INTAKE_PROCESSES.map((name, p) => p);

function intakeArrays(C, L, mask, cellOcean) {
  const slot = new Int32Array(C).fill(-1), cells = [];
  for (let i = 0; i < C; i++) if (cellOcean[i] && (mask === true || mask[i])) { slot[i] = cells.length; cells.push(i); }
  const size = INTAKE_PROCESSES.length * L * cells.length;
  return { processes: INTAKE_PROCESSES, cells: Int32Array.from(cells), slot, water: new Float64Array(size), degrees: new Float64Array(size) };
}
/*
 * The interior classes, placed where the World Ocean Atlas holds its water
 * in this equation of state: 1019.5–1021.5 for the warm pool, the two
 * lightest at WOA23's tropical surface salinities, 0.25 apart
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
  1019.5, 1020.0, 1020.5, 1021.0, 1021.5,
  1022.0, 1022.25, 1022.5, 1022.75, 1023.0, 1023.25, 1023.5, 1023.75, 1024.0, 1024.25, 1024.5, 1024.75, 1025.0, 1025.25, 1025.5, 1025.75, 1026.0,
  1026.1, 1026.2, 1026.3, 1026.4, 1026.5, 1026.6,
  1026.65, 1026.68, 1026.71, 1026.74, 1026.77, 1026.8, 1026.83, 1026.86, 1026.89, 1026.92, 1026.95, 1026.98,
  1027.05, 1027.15, 1027.3, 1027.5, 1027.75, 1028.0,
];
export const LAYER_BOTTOMS = [
  60, 62, 65, 74, 83,
  90, 110, 130, 150, 170, 205, 235, 270, 300, 350, 400, 450, 500, 550, 600, 650, 689,
  734, 801, 867, 933, 1004, 1061,
  1092, 1180, 1300, 1420, 1540, 1735, 2005, 2275, 2650, 3550, 12000, 12000,
  12000, 12000, 12000, 12000, 12000,
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
export const POLAR_INTERIOR_T = 272.15;

/*
 * The fields at a depth below the surface, per cell: the layer holding
 * that depth (the first whose base lies below it, interior layers no
 * thicker than `thin`, THIN for the slab's tokens, passed over), its temperature, its
 * current as a cell vector and as a speed, and the vertical velocity
 * there, upward positive, as the divergence of the transport above the
 * depth with the free surface held: h and u are the L layers'
 * thicknesses and edge velocities, temperature(k, i) the layer
 * temperature. Cells that are land or whose water column ends above the
 * depth are NaN, their current vector zero.
 */
export function depthFields(mesh, L, { h, u, temperature, cellOcean }, depth, { thin = THIN } = {}) {
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
export const EPS = 0.01, THIN = 5, PV_FLOOR = 20, SPEED_LIMIT = 5, DENSITY_TOLERANCE = 0.005, RESTORE_TOLERANCE = 0.01, THINNED_MARGIN = 0.01;

export { FRICTION_PER_WIND, CALM_WIND, RETREAT_EKMAN, erfc, stokesSurfaceLayer, langmuirFactor, JERLOV, PENETRATION_DEPTH, transmitted, deposited };
export function mixingOptionCheck({ stressThickness, langmuir, retreatEkman = 0, retreatDaily = false, retreatLight = 0, penetration = null, blockDrag = false, shearMixing = false, interiorShearMixing = false }) {
  if (blockDrag && (shearMixing || interiorShearMixing)) throw new Error('blockDrag sets the viscosity below the Ekman block itself, so it excludes shearMixing and interiorShearMixing');
  if (!(retreatEkman >= 0 && retreatEkman < Infinity)) throw new Error(`retreatEkman bounds the retreat at that multiple of u*/|f|, so it is 0 (no bound) or positive, not ${retreatEkman}`);
  if (!(retreatLight >= 0 && retreatLight <= 1)) throw new Error(`retreatLight weighs the light the night's convection carries down to the mixed layer's base, from 0 (no gate) to 1, not ${retreatLight}`);
  if (retreatLight > 0 && !retreatDaily) throw new Error('retreatLight gates the retreat on the previous day\'s sunlight and dark loss, so it needs retreatDaily');
  if (retreatLight > 0 && !penetration) throw new Error('retreatLight weighs the sunlight absorbed at the mixed layer\'s base, so it needs penetration');
  if (!(stressThickness >= THIN)) throw new Error(`stressThickness is at least THIN (${THIN} m), the thinnest class that keeps its own velocity, not ${stressThickness}`);
  if (langmuir !== 0 && !(langmuir >= 1)) throw new Error(`langmuir caps a factor of at least 1, so it is 0 or at least 1, not ${langmuir}`);
}
/*
 * The day's end for retreatDaily, from the model time at the end of an
 * ocean call `dt` long (day 0 starts at 00 UTC): closeDay adds the call to
 * `day.span` and returns 0 while the day goes on, and at the call that
 * reaches a UTC midnight 1/span when the span is a whole day, else −1,
 * starting the next span. Both engines take the flag from here.
 */
export const DAY_SECONDS = 86400;
export const dayEnds = (time, dt) => Math.floor((time + 1) / DAY_SECONDS) > Math.floor((time + 1 - dt) / DAY_SECONDS);
export function closeDay(day, time, dt) {
  if (!Number.isFinite(time)) throw new Error('retreatDaily needs the model time at the end of each ocean call, advance\'s last argument');
  day.span += dt;
  if (!dayEnds(time, dt)) return 0;
  const closing = Math.abs(day.span - DAY_SECONDS) < 0.5 * dt ? 1 / day.span : -1;
  day.span = 0;
  return closing;
}

/*
 * The fraction of the light a column takes below its mixed layer,
 * transmitted(h₀) − transmitted(depth), 0 unless an interior class is
 * thicker than THIN (a receiver), and in `share` each class's part of
 * it: a receiver takes what is absorbed over its own depths and over the
 * thin classes above it, and the last receiver also what passes it, so
 * the shares sum to 1. layeredOcean.gpu.js's oReadSurface follows the
 * same loop.
 */
export function penetrationShares(h, L, depth, type = 'IB', share = null) {
  let receiver = false;
  for (let k = 1; k < L && !receiver; k++) receiver = h[k] > THIN;
  const taken = receiver ? transmitted(h[0], type) - transmitted(depth, type) : 0;
  if (!share) return taken;
  share.fill(0);
  if (!(taken > 0)) return taken;
  let z = h[0], top = transmitted(z, type), carry = 0, given = 0, last = 0;
  for (let k = 1; k < L; k++) {
    z += h[k];
    const bottom = transmitted(z, type), part = carry + (top - bottom) / taken;
    if (h[k] > THIN) { share[k] = part; given += part; carry = 0; last = k; } else carry = part;
    top = bottom;
    if (z >= PENETRATION_DEPTH && last > 0) break;
  }
  share[last] += 1 - given;
  return taken;
}

/*
 * Columns from an ocean climatology (./climatology.module.js), sampled at
 * each sea cell's centre down to its bottom D. Where the atlas column
 * ends a standard level or more above D, as off coasts that the cell's
 * mean depth runs past, the levels below it come from the nearest atlas
 * column that reaches them (`columnReaching`); the atlas's deepest values
 * are carried on down where the model is deeper still. The profile's potential
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
    let column = atlas.columnAt(mesh.latCell[i], mesh.lonCell[i]);
    if (!column) continue;
    const reached = column.depths.length;
    if (atlas.columnReaching && atlas.depths && reached < atlas.depths.length && atlas.depths[reached] <= D[i]) {
      const deep = atlas.columnReaching(mesh.latCell[i], mesh.lonCell[i], D[i]);
      if (deep && deep.depths.length > reached) {
        const T = Float64Array.from(deep.T), S = Float64Array.from(deep.S);
        T.set(column.T); S.set(column.S);
        column = { depths: deep.depths, T, S };
      }
    }
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

/*
 * The velocity an interior class's ∇⁴ closure acts on: its own u, except
 * on token sea edges, which take `weight` of the normal component of the
 * uniform flow that best fits, in least squares with the ridge
 * CLOSURE_RIDGE, the normal velocities of their neighbours (edgesOnEdge)
 * where the class holds at least THIN metres, and the rest of their own.
 * A token edge's own velocity is the layer's above, the mixed layer's
 * where the class has outcropped, toward which the closure would pull the
 * class. With `rings` ({ deepest, k, valid, second }, `deepest` per edge
 * the deepest class holding more than THIN metres in both of its cells),
 * only the token edges with such a class denser than k beneath them are
 * fitted, the class pinched out between thicker ones and not cut off by
 * the sea floor, and then the token edges next to those from the fitted
 * and the thick edges; without it, every token edge beside the class.
 * `valid` is left 1 on the thick edges, 2 on the first ring and 3 on the
 * second, for closureAdjoint.
 */
export const CLOSURE_RIDGE = 0.1;
function closureWeights(mesh, e, usable, visit) {
  const { nEdgesOnEdge, edgesOnEdge, maxEdgesOnEdge, nEdge, tEdge } = mesh;
  const nx = nEdge[3 * e], ny = nEdge[3 * e + 1], nz = nEdge[3 * e + 2], tx = tEdge[3 * e], ty = tEdge[3 * e + 1], tz = tEdge[3 * e + 2];
  const along = (o) => [nEdge[3 * o] * nx + nEdge[3 * o + 1] * ny + nEdge[3 * o + 2] * nz, nEdge[3 * o] * tx + nEdge[3 * o + 1] * ty + nEdge[3 * o + 2] * tz];
  let saa = 0, sab = 0, sbb = 0, found = false;
  for (let s = 0; s < nEdgesOnEdge[e]; s++) {
    const o = edgesOnEdge[maxEdgesOnEdge * e + s];
    if (!usable(o)) continue;
    const [a, b] = along(o);
    saa += a * a; sab += a * b; sbb += b * b; found = true;
  }
  if (!found) return false;
  const p = saa + CLOSURE_RIDGE, q = sbb + CLOSURE_RIDGE, det = p * q - sab * sab;
  for (let s = 0; s < nEdgesOnEdge[e]; s++) {
    const o = edgesOnEdge[maxEdgesOnEdge * e + s];
    if (!usable(o)) continue;
    const [a, b] = along(o);
    visit(o, (q * a - sab * b) / det);
  }
  return true;
}
function closureFit(mesh, e, usable, values) {
  let fit = 0;
  return closureWeights(mesh, e, usable, (o, c) => { fit += c * values[o]; }) ? fit : null;
}
export function closureVelocity(mesh, u, hEdge, edgeOcean, weight, out, rings = null) {
  const E = mesh.nEdges, thick = (o) => edgeOcean[o] && hEdge[o] >= THIN;
  const inside = rings ? (e) => rings.deepest[e] > rings.k : () => true;
  for (let e = 0; e < E; e++) {
    out[e] = u[e];
    if (rings) rings.valid[e] = thick(e) ? 1 : 0;
    if (!edgeOcean[e] || hEdge[e] >= THIN || !inside(e)) continue;
    const fit = closureFit(mesh, e, thick, u);
    if (fit === null) continue;
    out[e] = weight * fit + (1 - weight) * u[e];
    if (rings) rings.valid[e] = 2;
  }
  if (!rings) return out;
  const { valid, second } = rings, usable = (o) => valid[o] === 1 || valid[o] === 2;
  second.set(out);
  for (let e = 0; e < E; e++) {
    if (valid[e] || !edgeOcean[e] || !inside(e)) continue;
    const fit = closureFit(mesh, e, usable, out);
    if (fit !== null) { second[e] = weight * fit + (1 - weight) * u[e]; valid[e] = 3; }
  }
  out.set(second);
  return out;
}
/*
 * The closure's ∇⁴ of the filled velocity, `closure`, made dissipative:
 * what it gives the fitted token edges goes back to the edges they were
 * fitted from, through the transpose of the fit under the edges' inner
 * product dc·dv, the second ring's to the first and the first's to the
 * thick edges. Without it the extrapolated tokens give ∇⁴ growing
 * grid-scale modes along the class's edge. Reads the rings that
 * closureVelocity left and adds to `closure` on the thick edges.
 */
export function closureAdjoint(mesh, closure, weight, rings) {
  const { nEdges: E, dcEdge, dvEdge } = mesh, { valid, second: carried } = rings;
  carried.fill(0);
  for (const ring of [3, 2]) {
    for (let e = 0; e < E; e++) {
      if (valid[e] !== ring) continue;
      const r = weight * (closure[e] + carried[e]) * dcEdge[e] * dvEdge[e];
      closureWeights(mesh, e, (o) => valid[o] === 1 || valid[o] === ring - 1, (o, c) => { carried[o] += c * r / (dcEdge[o] * dvEdge[o]); });
    }
  }
  for (let e = 0; e < E; e++) if (valid[e] === 1) closure[e] += carried[e];
  return closure;
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
export const TRACER_ADVECTIONS = ['upwind', 'third'], SHEAR_DIFFUSIVITY = 0.2, HYCOM_DIFFUSION_SPEED = 0.005;
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
 * and is cut or extended at the bottom, the mixed layer keeps its floor,
 * and a column that already fits is left exactly as it is. A sea cell
 * with no usable water (new coast, or a stencil with no sea source) takes
 * the area-weighted mean column and sea level of its sea neighbors that
 * have water, sweep by sweep inward, and is then fitted like any other;
 * one no sweep reaches takes the climatology column (built on first
 * need). Returns the cells filled from their neighbors and from the
 * climatology.
 */
export function fitColumns({ h, Q, W, eta }, climatology, { D, cellOcean, L, C, labelT, labelS, minimumThickness, mesh = null }) {
  const at = (k, i) => k * C + i;
  const usable = (i) => {
    if (!Number.isFinite(eta[i])) return false;
    let sum = 0;
    for (let k = 0; k < L; k++) { const n = at(k, i); if (!(h[n] >= 0) || !Number.isFinite(Q[n]) || !Number.isFinite(W[n])) return false; sum += h[n]; }
    return sum >= 1;
  };
  const fit = (i) => {
    let sum = 0;
    for (let k = 0; k < L; k++) sum += h[at(k, i)];
    if (Math.abs(sum - D[i] - eta[i]) < 1e-6) return;
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
  };
  const wet = new Uint8Array(C);
  let waiting = [];
  for (let i = 0; i < C; i++) {
    if (!cellOcean[i]) { for (let k = 0; k < L; k++) { h[at(k, i)] = 0; Q[at(k, i)] = 0; W[at(k, i)] = 0; } eta[i] = 0; continue; }
    if (usable(i)) { fit(i); wet[i] = 1; } else waiting.push(i);
  }
  const fromNeighbors = [];
  if (mesh) {
    const { cellsOnCell, nEdgesOnCell, maxEdges, areaCell } = mesh;
    while (waiting.length) {
      const reached = [], left = [];
      for (const i of waiting) {
        let area = 0;
        for (let m = 0; m < nEdgesOnCell[i]; m++) { const j = cellsOnCell[maxEdges * i + m]; if (wet[j]) area += areaCell[j]; }
        if (area > 0) reached.push(i); else left.push(i);
      }
      if (!reached.length) break;
      for (const i of reached) {
        let area = 0, level = 0;
        for (let k = 0; k < L; k++) { h[at(k, i)] = 0; Q[at(k, i)] = 0; W[at(k, i)] = 0; }
        for (let m = 0; m < nEdgesOnCell[i]; m++) {
          const j = cellsOnCell[maxEdges * i + m];
          if (!wet[j]) continue;
          const a = areaCell[j];
          area += a; level += a * eta[j];
          for (let k = 0; k < L; k++) { h[at(k, i)] += a * h[at(k, j)]; Q[at(k, i)] += a * Q[at(k, j)]; W[at(k, i)] += a * W[at(k, j)]; }
        }
        for (let k = 0; k < L; k++) { h[at(k, i)] /= area; Q[at(k, i)] /= area; W[at(k, i)] /= area; }
        eta[i] = level / area;
      }
      for (const i of reached) { fit(i); wet[i] = 1; fromNeighbors.push(i); }
      waiting = left;
    }
  }
  for (const i of waiting) {
    const clim = typeof climatology === 'function' ? (climatology = climatology()) : climatology;
    let s = 0;
    for (let k = 0; k < L; k++) { const n = at(k, i); h[n] = clim.h[n]; Q[n] = clim.Q[n]; W[n] = clim.W[n]; s += h[n]; }
    eta[i] = s - D[i];
  }
  return { fromNeighbors, fromClimatology: waiting };
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
 * span, half-way to its neighbors' labels. The lightest and densest,
 * which hold all lighter and all denser water, keep that inner boundary;
 * when their water lies beyond the label by more than RESTORE_TOLERANCE,
 * the span reaches as far beyond the water's density as the boundary lies
 * within it, else as far beyond the label (a lone class spans the same
 * width centered on its water). Each new class, the lightest and densest
 * open-ended, takes the share of that span within its own, the shares
 * tilted linearly across the span (and none below zero) so that their
 * labels average to the density of the water. Each share keeps the old
 * class's temperature and salinity moved by the least change, counted in
 * 4 K and 0.5 psu (their spread along a density surface), that changes
 * the density by ρ_new − ρ̄, ρ̄ the shares' mean label: it lies near its
 * new label, warmer or colder in the thermocline and saltier or fresher
 * in cold water, and the column keeps its water, heat and salt. A class
 * left with less than the EPS token is made up to it at its label from
 * the nearest class holding more than the token and two EPS besides, or
 * from the mixed layer when none does. A new class moves with the old
 * classes whose label spans overlap its own, weighted by the overlap, or
 * else with the old class of the nearest label. The per-cell arrays a
 * restart takes (REBIN_CARRIED) are carried as they are and Q and W given
 * as h·T and h·S, so a rebinned ocean restarts on either engine.
 */
const SPREAD_T = 4, SPREAD_S = 0.5;
export const DAILY = ['sunSum', 'darkSum', 'sunDay', 'darkDay'];
export const REBIN_CARRIED = ['previousT0', 'previousIceVolume', 'fresh', 'buoyancy', 'remainder', 'flux', 'capacity', 'penetrated', ...DAILY, 'daySpan'];
export function rebinOcean(saved, from, to, C, { labelT, labelS }) {
  const Lo = from.length + 1, Ln = to.length + 1, E = saved.u.length / Lo;
  const spans = (labels) => labels.map((r, k) => {
    const below = k > 0 ? 0.5 * (labels[k - 1] + r) : null, above = k < labels.length - 1 ? 0.5 * (r + labels[k + 1]) : null;
    return [below ?? 2 * r - above, above ?? 2 * r - below];
  });
  const oldSpans = spans(from), newSpans = spans(to), receiving = newSpans.map(([a, b], n) => [n === 0 ? -Infinity : a, n === to.length - 1 ? Infinity : b]);
  const h = new Float64Array(Ln * C), Q = new Float64Array(Ln * C), W = new Float64Array(Ln * C);
  const overlaps = (a, b) => {
    const out = [];
    let total = 0;
    for (let m = 0; m < to.length; m++) {
      const f = Math.min(b, receiving[m][1]) - Math.max(a, receiving[m][0]);
      if (f > 0) { out.push([m, f]); total += f; }
    }
    for (const share of out) share[1] /= total;
    return out;
  };
  const fixed = oldSpans.map(([a, b]) => overlaps(a, b));
  for (let i = 0; i < C; i++) {
    h[i] = saved.h[i]; Q[i] = h[i] > 0 ? h[i] * saved.T[i] : 0; W[i] = h[i] > 0 ? h[i] * saved.S[i] : 0;
    let column = h[i];
    for (let j = 1; j < Lo; j++) {
      const n = j * C + i, hj = saved.h[n], t = saved.T[n], s = saved.S[n];
      if (!(hj >= 0) || (hj > 0 && !(Number.isFinite(t) && Number.isFinite(s)))) { h[i] = NaN; continue; }
      if (hj === 0) continue;
      column += hj;
      let shares = fixed[j - 1], r = NaN;
      if (j === 1 || j === Lo - 1) {
        const label = from[j - 1];
        let [a, b] = oldSpans[j - 1];
        r = seawaterDensity(t, s);
        if (Lo === 2) { const half = 0.5 * (b - a); a = r - half; b = r + half; } else if (j === 1 && r < label - RESTORE_TOLERANCE) a = 2 * r - b; else if (j === Lo - 1 && r > label + RESTORE_TOLERANCE) b = 2 * r - a;
        if (a !== oldSpans[j - 1][0] || b !== oldSpans[j - 1][1]) shares = overlaps(a, b);
      }
      if (shares.length === 1) {
        const k = (shares[0][0] + 1) * C + i;
        h[k] += hj; Q[k] += hj * t; W[k] += hj * s;
        continue;
      }
      shares = shares.map(([m, f]) => [m, f]);
      let centre = 0, spread = 0, mean = 0, kept = 0;
      for (const [m, f] of shares) centre += f * to[m];
      for (const [m, f] of shares) spread += f * (to[m] - centre) ** 2;
      const tilt = spread > 0 ? ((Number.isNaN(r) ? seawaterDensity(t, s) : r) - centre) / spread : 0;
      for (const share of shares) { share[1] = Math.max(0, share[1] * (1 + tilt * (to[share[0]] - centre))); kept += share[1]; }
      for (const share of shares) { share[1] /= kept; mean += share[1] * to[share[0]]; }
      const byT = -expansionOf(t, s) * SEAWATER.rho0, byS = contractionOf(t, s) * SEAWATER.rho0, norm = (byT * SPREAD_T) ** 2 + (byS * SPREAD_S) ** 2;
      for (const [m, f] of shares) {
        const k = (m + 1) * C + i, part = f * hj, change = (to[m] - mean) / norm;
        h[k] += part; Q[k] += part * (t + byT * SPREAD_T ** 2 * change); W[k] += part * (s + byS * SPREAD_S ** 2 * change);
      }
    }
    if (!(column > 0)) continue;
    for (let k = 1; k < Ln; k++) {
      const n = k * C + i, need = EPS - h[n];
      if (!(need > 0)) continue;
      let giver = i;
      for (let d = 1; d < Ln && giver === i; d++) {
        if (k - d >= 1 && h[(k - d) * C + i] > need + 2 * EPS) giver = (k - d) * C + i;
        else if (k + d < Ln && h[(k + d) * C + i] > need + 2 * EPS) giver = (k + d) * C + i;
      }
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
  const T = new Float64Array(Ln * C), S = new Float64Array(Ln * C);
  for (let n = 0; n < Ln * C; n++) if (h[n] > 0) { T[n] = Q[n] / h[n]; S[n] = W[n] / h[n]; }
  const carried = {};
  for (const name of REBIN_CARRIED) if (saved[name]) carried[name] = Float64Array.from(saved[name]);
  return { h, u, T, S, Q, W, eta: Float64Array.from(saved.eta), densities: Array.from(to), ...carried };
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

/*
 * The depth each edge's column is clipped to: the shallower cell, and for
 * the hybrid ocean the geography's edge sill where it is shallower still,
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

export function createOcean(mesh, {
  densities = LAYER_DENSITIES, salinities = LAYER_SALINITIES, bottoms = LAYER_BOTTOMS, mixedDepth = 60, minimumDepth = 50, flatDepth = 4000, thermoclineTilt = 0.3,
  salinityProfile = (lat) => 34 + 2 * Math.exp(-(((Math.abs(lat) * 180 / Math.PI - 25) / 20) ** 2)),
  density = 1025, specificHeat = 3985, referenceS = 35, gravity = 9.81,
  minimumThickness = 50, shallowestMixedDepth = 50, maximumMixedDepth = 600, convectiveRate = 100 / 86400, neutralSnap = false, convectiveErosion = true, buoyancyMemory = 86400, mixedNeighbourRatio = 0, vorticityCentring = 0.5, stirring = 0.8, stirringDepth = 100, detrainmentTime = 86400, restoreTime = 2 * 86400, iceStressTransmission = 0.8, iceSalinity = 5, iceDensity = 917,
  momentumDepth = 50, carryMomentum = true, refillRate = 0.00027778, refillFloor = 5, shearEntrainment = 0.65,
  ekmanDepth = 0.4, momentumFloor = 10, stressThickness = PV_FLOOR, retreatStirring = 2.5, langmuir = 2, retreatEkman = RETREAT_EKMAN, retreatDaily = false, retreatLight = 0,
  interfacialDrag = 2e-4, shearMixing = false, interiorShearMixing = false, blockDrag = false, shearViscosity = 1e-2, backgroundViscosity = 1e-4, bottomDrag = 3e-3, closureHours = 12, closureSpacing = CLOSURE_SPACING, closureFill = 0, closureTokens = 'interior', diffusivity = 0.01, mixedDiffusivity = 0, everySteps = 4,
  eddyDiffusivity = 1000, eddyTaperDepth = 200, ekmanBlock = false, upwellingSource = 0, intakeLedger = false, penetration = 'IB',
  geography = null, bathymetry = null, edgeSill = undefined, sillFloor = minimumDepth, buffers = null, climatology = null,
  vertical = 'slab', generator = 'hybrid', layerTable = HYBRID_LAYERS, regridTime = 0, minimumLayer = H_MIN, modeDepth = MODE_DEPTH, mstar = 1.25, ekmanFactor = 0.5, nstar = 0.2, backgroundDiffusivity = 1e-5, convectiveDiffusivity = 0.1,
  tracerAdvection = null, shearDiffusivity = null, shearDiffusivityCap = 4000, lateralDiffusionSpeed = null,
} = {}) {
  const {
    nCells: C, nEdges: E, nVertices: V, maxEdgesOnEdge, nEdgesOnEdge, edgesOnEdge, weightsOnEdge,
    cellsOnEdge, verticesOnEdge, areaCell, dcEdge, dvEdge, fVertex, fEdge, fCell, radius, latCell,
    maxEdges, nEdgesOnCell, edgesOnCell, edgeSignOnCell, cellsOnCell, nEdge,
  } = mesh;
  if (climatology && typeof climatology.columnAt !== 'function') throw new Error('the ocean climatology must be decoded (loadClimatology in ./climatology.module.js) before the ocean is built');
  if (penetration !== null && !JERLOV[penetration]) throw new Error(`penetration is one of ${Object.keys(JERLOV).join(', ')} or null, not ${penetration}`);
  if (shearEntrainment > 0 && !carryMomentum) throw new Error('shearEntrainment needs carryMomentum: without the entrained water\'s momentum mixed in, the base\'s Richardson number never rises and the mixed layer swallows a class a step');
  if (!(upwellingSource >= 0) || (upwellingSource > 0 && !ekmanBlock)) throw new Error(`upwellingSource is a share of the Ekman divergence, not negative, and needs ekmanBlock for the transport it divides; not ${upwellingSource} with ekmanBlock ${ekmanBlock}`);
  mixingOptionCheck({ stressThickness, langmuir, retreatEkman, retreatDaily, retreatLight, penetration, blockDrag, shearMixing, interiorShearMixing });
  const hybrid = vertical === 'hybrid';
  if (vertical !== 'slab' && !hybrid) throw new Error(`vertical is 'slab' or 'hybrid', not ${vertical}`);
  if (hybrid && (ekmanBlock || upwellingSource > 0 || intakeLedger || retreatDaily || closureFill > 0 || shearMixing || interiorShearMixing || blockDrag)) throw new Error('ekmanBlock, upwellingSource, intakeLedger, retreatDaily, closureFill, shearMixing, interiorShearMixing and blockDrag belong to the slab ocean');
  const advection = tracerAdvection ?? (hybrid ? 'third' : 'upwind');
  if (!TRACER_ADVECTIONS.includes(advection)) throw new Error(`tracerAdvection is one of ${TRACER_ADVECTIONS.join(', ')}, not ${advection}`);
  if (!hybrid && advection !== 'upwind') throw new Error("tracerAdvection 'third' belongs to the hybrid ocean");
  const beyond = hybrid && advection === 'third' ? beyondCells(mesh) : null;
  const shearDiffusion = shearDiffusivity ?? (hybrid ? SHEAR_DIFFUSIVITY : 0);
  if (!hybrid && shearDiffusion > 0) throw new Error('shearDiffusivity belongs to the hybrid ocean');
  if (!(shearDiffusion >= 0) || !(shearDiffusivityCap > 0)) throw new Error(`shearDiffusivity is a non-negative coefficient and shearDiffusivityCap a positive diffusivity, not ${shearDiffusion} and ${shearDiffusivityCap}`);
  if (advection === 'third' && !(shearDiffusion > 0)) throw new Error('the limited third-order tracer flux admits grid-scale noise that the shear-scaled diffusion of the surface layers damps, so tracerAdvection third needs shearDiffusivity above zero');
  const lateralSpeed = lateralDiffusionSpeed ?? (hybrid ? HYCOM_DIFFUSION_SPEED : 0);
  if (!hybrid && lateralSpeed > 0) throw new Error('lateralDiffusionSpeed belongs to the hybrid ocean');
  if (!(lateralSpeed >= 0)) throw new Error(`lateralDiffusionSpeed is a non-negative velocity scale, not ${lateralSpeed}`);
  const table = hybrid ? layerTable : null;
  const L = hybrid ? table.K : densities.length + 1, g = gravity, rho0 = density, rhoCp = density * specificHeat;
  const rho = hybrid ? Array.from(table.densities) : [rho0, ...densities];
  const labels = hybrid ? Array.from(table.densities) : densities;
  const labelS = hybrid ? rho.map((r) => { const j = LAYER_DENSITIES.indexOf(r); return j >= 0 ? LAYER_SALINITIES[j] : referenceS; }) : [referenceS, ...salinities];
  const labelT = rho.map((r, k) => Math.max(FREEZING_POINT, labelTemperature(r, labelS[k])));
  const thermoclineLayers = rho.filter((r, k) => k > 0 && r < THERMOCLINE_DENSITY).length;
  const diffusion = hybrid ? 0 : diffusivity * radius * radius / rhoCp;
  let spacing = 0, minSpacing = Infinity;
  for (let e = 0; e < E; e++) { spacing += dcEdge[e]; minSpacing = Math.min(minSpacing, dcEdge[e]); }
  spacing /= E;
  const nu4 = closureCoefficient(spacing, closureHours, closureSpacing);
  const eddyKappa = eddyDiffusivities(mesh, eddyDiffusivity), eddyLimit = eddyDiffusivity > 0 ? eddyDiffusionLimit(mesh) : 0;
  let eddy = null, ekman = null, ekmanLimited = 0;
  const ekmanFlux = ekmanBlock ? new Float64Array(E) : null, ekmanMoved = ekmanBlock ? new Float64Array(L) : null;

  const edgeOcean = geography ? geography.edgeOcean : new Uint8Array(E).fill(1);
  const cellOcean = geography ? Uint8Array.from(geography.land, (l) => (l ? 0 : 1)) : new Uint8Array(C).fill(1);
  const outlet = runoffOutlets(mesh, geography);
  const coastal = upwellingSource > 0 ? Uint8Array.from({ length: C }, (_, i) => {
    if (!cellOcean[i]) return 0;
    for (let m = 0; m < nEdgesOnCell[i]; m++) if (!cellOcean[cellsOnCell[maxEdges * i + m]]) return 1;
    return 0;
  }) : null;
  const upwelled = upwellingSource > 0 ? new Float64Array(C) : null;
  const intake = intakeLedger ? intakeArrays(C, L, intakeLedger, cellOcean) : null;
  const D = bathymetry ? Float64Array.from(bathymetry, (d, i) => (cellOcean[i] ? d : 0)) : bathymetryFrom(mesh, geography, { minimumDepth, flatDepth });
  const sillEdge = edgeSills(mesh, D, edgeOcean, hybrid ? (edgeSill === undefined ? (geography ? geography.edgeSill ?? null : null) : edgeSill) : null, sillFloor);
  let deepest = 0;
  for (let i = 0; i < C; i++) deepest = Math.max(deepest, D[i]);
  const substepLimit = 0.35 * minSpacing / Math.sqrt(g * Math.max(deepest, 1));

  const adopting = !!(buffers && buffers.h);
  const shared = (name, n) => new Float64Array(buffers && buffers[name] ? buffers[name] : new SharedArrayBuffer(8 * n));
  const h = shared('h', L * C), u = shared('u', L * E), Q = shared('Q', L * C), W = shared('W', L * C);
  const state = [h, u, Q, W];
  const eta = shared('eta', C);
  const T0 = new Float64Array(C), S0 = new Float64Array(C), rhoMl = shared('rhoMl', C), previousT0 = new Float64Array(C), surfaceIn = new Float64Array(C), iceVolume = new Float64Array(C), previousVolume = new Float64Array(C);
  let seaCover = null;
  const capacity = shared('capacity', C);
  if (!adopting) capacity.fill(rhoCp * mixedDepth);
  const transmission = shared('transmission', C), penetrated = shared('penetrated', C), column = new Float64Array(L), shares = new Float64Array(L);
  const stress = shared('stress', E), fresh = new Float64Array(C);
  const daily = retreatDaily ? { sunStep: new Float64Array(C), sunSum: new Float64Array(C), darkSum: new Float64Array(C), sunDay: new Float64Array(C), darkDay: new Float64Array(C) } : null, day = { span: 0 };
  const iced = new Uint8Array(C), iceCover = langmuir > 0 ? new Float64Array(C) : null;
  const hEdge = shared('hEdge', L * E), pressure = shared('pressure', L * C), flux = new Float64Array(E), fluxPV = new Float64Array(E), tracerFlux = new Float64Array(E), edgeS = new Float64Array(E);
  const stressDepth = shared('stressDepth', E), stressLayer = shared('stressLayer', E), stressFraction = shared('stressFraction', E);
  const windDepth = ekmanDepth > 0 ? shared('windDepth', E) : null;
  if (windDepth && !adopting) windDepth.fill(momentumDepth);
  const T = new Float64Array(C), S = new Float64Array(C), lapT = new Float64Array(C);
  const zeta = new Float64Array(V), qEdge = new Float64Array(E);
  const K = new Float64Array(C), phi = new Float64Array(C), gradPhi = new Float64Array(E), gradEta = shared('gradEta', E), gradRho = shared('gradRho', E);
  const lap = new Float64Array(E), lap2 = new Float64Array(E), divScratch = new Float64Array(C), curlScratch = new Float64Array(V), closureU = new Float64Array(E);
  if (closureTokens !== 'interior' && closureTokens !== 'beside') throw new Error(`closureTokens is 'interior' or 'beside', not ${closureTokens}`);
  const deepestEdge = shared('deepestEdge', E), closureRings = closureTokens === 'interior' && closureFill > 0 ? { deepest: deepestEdge, k: 0, valid: new Uint8Array(E), second: new Float64Array(E) } : null;
  const slow = shared('baroSlow', E), U = shared('baroU', E), depthEdge = shared('baroDepth', E), etaB = shared('baroEta', C), avgU = shared('baroAvgU', E), avgEta = shared('baroAvgEta', C), divU = new Float64Array(C);
  const countBuffer = buffers && buffers.counts ? buffers.counts : new SharedArrayBuffer(16), counts = new Int32Array(countBuffer);
  const stages = [0, 1, 2, 3].map((s) => [shared(`stage${s}h`, L * C), shared(`stage${s}u`, L * E), shared(`stage${s}Q`, L * C), shared(`stage${s}W`, L * C)]);
  const trial = [shared('trialh', L * C), shared('trialu', L * E), shared('trialQ', L * C), shared('trialW', L * C)];
  const params = shared('params', 8);
  if (!adopting) { params[0] = 1 / 3600; params[1] = 3600; }
  const tauCell = new Float64Array(3 * C), reach = new Float64Array(C).fill(Infinity), buoyancyLoss = new Float64Array(C);
  const refilled = { water: new Float64Array(C), heat: new Float64Array(C) };
  const sheared = new Float64Array(C), retreatDepth = new Float64Array(C), baseShear = shearEntrainment > 0 ? new Float64Array(E) : null, shearCell = shearEntrainment > 0 ? new Float64Array(3 * C) : null;
  const lightest = rho.map((r, k) => r - 0.5 * (k > 1 ? r - rho[k - 1] : rho[k + 1] - r));
  const densest = rho.map((r, k) => r + 0.5 * (k < L - 1 ? rho[k + 1] - r : r - rho[k - 1]));
  let counter = 0, limited = 0, initialised = false, layerRunner = null, columnRunner = null, dried = 0;
  const relaxRate = () => params[0];
  const layerRho = hybrid ? shared('layerRho', L * C) : null, columnAbove = hybrid ? shared('columnAbove', L * C) : null, targetZ = hybrid ? shared('targetZ', (L + 1) * C) : null, viscosity = hybrid ? shared('viscosity', L * C) : null, edgeBottom = hybrid ? shared('edgeBottom', E) : null;
  const mixedDepths = hybrid ? shared('mixedDepths', C) : null, boundaryDepths = hybrid ? shared('boundaryDepths', C) : null;
  const edgeAbove = hybrid ? shared('edgeAbove', L * E) : null, eddyFlux = hybrid ? shared('eddyFlux', L * E) : null, eddyVolume = hybrid ? shared('eddyVolume', L * C) : null, eddyHeat = hybrid ? shared('eddyHeat', L * C) : null, eddySalt = hybrid ? shared('eddySalt', L * C) : null;
  const pvFloor = hybrid ? Float64Array.from(table.minimum, (m) => Math.min(PV_FLOOR, m)) : null, tokenEdge = hybrid ? Float64Array.from(table.minimum, (m) => 2 * Math.min(m, table.spacing ?? m)) : null, pinned = hybrid ? new Uint8Array(E) : null;
  const boundaries = hybrid ? labelBoundaries(rho) : null, cs = hybrid ? columnScratch(L) : null;
  const cols = hybrid ? Object.fromEntries(['h', 'T', 'S', 'rho', 'shear', 'kappa', 'nu', 'hNew', 'amount', 'hOld', 'hEdgeNew', 'uCol', 'nuEdge'].map((name) => [name, new Float64Array(L)])) : null;
  if (hybrid) { cols.z = new Float64Array(L + 1); cols.w = new Float64Array(L + 1); cols.vector = new Float64Array(3 * L); }
  const interfaceWeight = (k, i) => (generator === 'isopycnal' ? 1 : generator === 'zstar' ? 0 : Math.min(1, Math.max(0, (h[at(k - 1, i)] - table.minimum[k - 1]) / modeDepth)));

  const eos = seawaterDensity;
  const at = (k, i) => k * C + i;
  const ae = (k, e) => k * E + e;

  function maskedLaplacian(field, out, thickness = null) {
    for (let i = 0; i < C; i++) {
      let sum = 0;
      for (let m = 0; m < nEdgesOnCell[i]; m++) {
        const e = edgesOnCell[maxEdges * i + m], j = cellsOnCell[maxEdges * i + m];
        if (!edgeOcean[e]) continue;
        const weight = thickness ? 2 * thickness[i] * thickness[j] / Math.max(thickness[i] + thickness[j], 1e-6) : 1;
        sum += weight * dvEdge[e] * (field[j] - field[i]) / dcEdge[e];
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
   * side's depth; the hybrid ocean instead cuts the column at the edge's
   * sill, the layers below it carrying nothing across.
   */
  function edgeColumn(hIn, e, out, base, stride) {
    const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
    out[base] = 0.5 * (hIn[a] + hIn[b]);
    for (let k = 1; k < L; k++) out[base + k * stride] = Math.min(hIn[at(k, a)], hIn[at(k, b)]);
    const sill = Math.max(EPS, sillEdge[e] + 0.5 * (eta[a] + eta[b]));
    if (hybrid) {
      let above = 0;
      for (let k = 0; k < L; k++) { const hk = out[base + k * stride]; out[base + k * stride] = Math.min(hk, Math.max(0, sill - above)); above += hk; }
      return;
    }
    let sum = 0;
    for (let k = 0; k < L; k++) sum += out[base + k * stride];
    if (sum > sill) { const f = sill / sum; for (let k = 0; k < L; k++) out[base + k * stride] *= f; }
  }
  function edgeThicknesses(hIn) {
    for (let e = 0; e < E; e++) edgeColumn(hIn, e, hEdge, e, E);
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

  function tangential(uIn, k, e) {
    let sum = 0;
    for (let s = 0; s < nEdgesOnEdge[e]; s++) { const o = edgesOnEdge[maxEdgesOnEdge * e + s]; sum += weightsOnEdge[maxEdgesOnEdge * e + s] * dvEdge[o] * uIn[ae(k, o)]; }
    return sum / dcEdge[e];
  }
  const dragThickness = (k, e) => Math.max(hEdge[ae(k, e)], THIN);
  function stressReach() {
    for (let e = 0; e < E; e++) {
      stressLayer[e] = -1; stressDepth[e] = 0; stressFraction[e] = 0;
      if (!edgeOcean[e]) continue;
      const reach = Math.max(windDepth ? windDepth[e] : momentumDepth, hEdge[e]);
      let top = 0, depth = 0;
      for (let k = 0; k < L; k++) {
        const hk = hEdge[ae(k, e)], part = Math.min(hk, Math.max(0, reach - top));
        if (k === 0 || hk >= stressThickness) depth += part;
        top += hk;
        if (top >= reach || k === L - 1) { stressLayer[e] = depth > 0 ? k : -1; stressFraction[e] = hk > 0 ? part / hk : 1; break; }
      }
      stressDepth[e] = depth;
    }
  }
  function windDepths() {
    for (let e = 0; e < E; e++) {
      const f = Math.abs(fEdge[e]);
      if (!edgeOcean[e] || !(f > 0)) { windDepth[e] = momentumDepth; continue; }
      const ustar = Math.sqrt(Math.hypot(stress[e], tangential(stress, 0, e)) / rho0);
      windDepth[e] = Math.min(momentumDepth, Math.max(momentumFloor, ekmanDepth * ustar / f));
    }
  }
  function windAcceleration(k, e) {
    const last = stressLayer[e];
    if (k > last || (k > 0 && hEdge[ae(k, e)] < stressThickness)) return 0;
    return stress[e] / rho0 * (k < last ? 1 : stressFraction[e]) / stressDepth[e];
  }
  /*
   * The drag coefficient r (m/s) between layer `up` and the layer `down`
   * beneath it at edge e, from the edge thicknesses and mixed-layer
   * density the tendency has prepared: interfacialDrag at the mixed
   * layer's base; between two interior classes under
   * `interiorShearMixing`, and at every interface under shearMixing, the
   * Pacanowski and Philander (1981) viscosity
   * ν = shearViscosity/(1 + 5 Ri)² + backgroundViscosity over the distance
   * Δz between the layers' middles, Ri = Δb Δz/|Δu|² from their buoyancy
   * step and the difference of their full velocities, the tangential part
   * reconstructed as the Coriolis term's, no less than interfacialDrag at
   * the mixed layer's base; r is held to half of what either layer's drag
   * thickness can take in one step (params[1] seconds). Under `blockDrag`,
   * interfacialDrag holds wherever the interface, the summed edge
   * thickness down to the lower layer, lies the distance Δz between the
   * layers' middles or more above the wind depth, the viscosity with no
   * floor from the wind depth down, the mixed layer's base included, and
   * the two blend linearly in the Δz above it.
   */
  function interfaceRate(uIn, e, up, down) {
    const interior = up > 0 && interiorShearMixing;
    if (!shearMixing && !interior && !blockDrag) return interfacialDrag;
    const dz = Math.max(THIN, 0.5 * (hEdge[ae(up, e)] + hEdge[ae(down, e)]));
    let held = 0;
    if (blockDrag) {
      let depth = 0;
      for (let j = 0; j < down; j++) depth += hEdge[ae(j, e)];
      held = Math.min(1, ((windDepth ? windDepth[e] : momentumDepth) - depth) / dz);
      if (held >= 1) return interfacialDrag;
      held = Math.max(0, held);
    }
    const upper = up === 0 ? 0.5 * (rhoMl[cellsOnEdge[2 * e]] + rhoMl[cellsOnEdge[2 * e + 1]]) : rho[up];
    const buoyancy = Math.max(0, g * (rho[down] - upper) / rho0);
    const du = uIn[ae(up, e)] - uIn[ae(down, e)], dv = tangential(uIn, up, e) - tangential(uIn, down, e);
    const richardson = buoyancy * dz / (du * du + dv * dv + 1e-12);
    const nu = shearViscosity / (1 + 5 * richardson) ** 2 + backgroundViscosity;
    const rate = Math.min(Math.max(interior || blockDrag ? 0 : interfacialDrag, nu / dz), 0.5 * Math.min(dragThickness(up, e), dragThickness(down, e)) / params[1]);
    return blockDrag ? rate + (interfacialDrag - rate) * held : rate;
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
      let weighted = 0, total = 0;
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
    const [hIn, , QIn, WIn] = input;
    if (hybrid) {
      if (columnRunner) { const useTrial = input === trial ? 1 : 0; columnRunner('prepareCells', useTrial); columnRunner('prepareEdges', useTrial); return; }
      prepareCells(0, C, input);
      prepareEdges(0, E, input);
      return;
    }
    surfaceDensity(hIn, QIn, WIn);
    gradient(mesh, eta, gradEta);
    gradient(mesh, rhoMl, gradRho);
    edgeThicknesses(hIn);
    stressReach();
    if (closureRings) {
      for (let e = 0; e < E; e++) {
        const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
        let ka = 0, kb = 0;
        for (let k = L - 1; k > 0 && !ka; k--) if (hIn[at(k, a)] > THIN) ka = k;
        for (let k = L - 1; k > 0 && !kb; k--) if (hIn[at(k, b)] > THIN) kb = k;
        deepestEdge[e] = Math.min(ka, kb);
      }
    }
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
    else if (hybrid) hybridTendencyLayers(input, out, 0, L);
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
      if (k === 0 && (mixedDiffusivity > 0 || diffusion > 0)) {
        const coefficient = mixedDiffusivity > 0 ? mixedDiffusivity : diffusion, thickness = mixedDiffusivity > 0 ? hIn : null;
        maskedLaplacian(T, lapT, thickness);
        for (let i = 0; i < C; i++) dQ[i] += coefficient * lapT[i];
        maskedLaplacian(S, lapT, thickness);
        for (let i = 0; i < C; i++) dW[i] += coefficient * lapT[i];
      }
      for (let i = 0; i < C; i++) if (!cellOcean[i]) { dh[oc + i] = 0; dQ[oc + i] = 0; dW[oc + i] = 0; }
      }
      if (!momentum) continue;
      curl(mesh, uIn.subarray(oe, oe + E), zeta);
      for (let e = 0; e < E; e++) qEdge[e] = 0.5 * (zeta[verticesOnEdge[2 * e]] + fVertex[verticesOnEdge[2 * e]] + zeta[verticesOnEdge[2 * e + 1]] + fVertex[verticesOnEdge[2 * e + 1]]) / Math.max(hEdge[oe + e], k > 0 ? vorticityCentring * 0.5 * (hIn[oc + cellsOnEdge[2 * e]] + hIn[oc + cellsOnEdge[2 * e + 1]]) : 0, k > 0 ? PV_FLOOR : THIN);
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
        const he = Math.max(hEdge[oe + e], momentumDepth);
        let force = 0, drag = 0, push = 0;
        if (windDepth) {
          if (k === 0 && hEdge[e] >= windDepth[e]) push = stress[e] / rho0 / Math.max(hEdge[e], windDepth[e]);
          else push = windAcceleration(k, e);
        } else if (k === 0 && hEdge[e] >= momentumDepth) force += stress[e] / rho0;
        else push = windAcceleration(k, e);
        if (k > 0) {
          let j = k - 1;
          while (j > 0 && hEdge[ae(j, e)] < THIN) j--;
          drag += interfaceRate(uIn, e, j, k) * (uIn[ae(j, e)] - uIn[oe + e]);
        }
        if (k < L - 1) {
          let j = k + 1;
          while (j < L - 1 && hEdge[ae(j, e)] < THIN) j++;
          if (hEdge[ae(j, e)] >= THIN) drag -= interfaceRate(uIn, e, k, j) * (uIn[oe + e] - uIn[ae(j, e)]);
        }
        let bottom = k === L - 1;
        if (!bottom) { bottom = true; for (let j = k + 1; j < L; j++) if (hEdge[ae(j, e)] >= THIN) { bottom = false; break; } }
        if (bottom) force -= bottomDrag * Math.abs(uIn[oe + e]) * uIn[oe + e];
        du[oe + e] += force / he + drag / dragThickness(k, e) + push;
      }
      if (nu4 > 0) {
        const uk = uIn.subarray(oe, oe + E);
        if (closureRings) closureRings.k = k;
        laplacianVelocity(mesh, k === 0 || !(closureFill > 0) ? uk : closureVelocity(mesh, uk, hEdge.subarray(oe, oe + E), edgeOcean, closureFill, closureU, closureRings), lap, divScratch, curlScratch);
        laplacianVelocity(mesh, lap, lap2, divScratch, curlScratch);
        if (closureRings && k > 0) closureAdjoint(mesh, lap2, closureFill, closureRings);
        for (let e = 0; e < E; e++) du[oe + e] -= nu4 * lap2[e];
      }
      if (k > 0) for (let e = 0; e < E; e++) if (hEdge[oe + e] < THIN) du[oe + e] = (uIn[ae(k - 1, e)] - uIn[oe + e]) * relax;
      for (let e = 0; e < E; e++) if (!edgeOcean[e]) du[oe + e] = 0;
    }
  }

  /*
   * The hybrid layers' tendencies: every layer alike. The flux through an
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
   * or one that occupies no depth in common between its two cells while
   * no thicker than tokenEdge (twice its minimum, at most twice the
   * thermocline spacing: a class held at its floor on one side against
   * water of that class elsewhere in the other column), follows the
   * layer above, relaxing to its velocity here and copied from it after
   * the barotropic shift, as the slab's tokens do; the
   * stress is a flux into the top layer; the bottom drag acts on the part
   * of each layer within momentumDepth of the floor, both measured down
   * the edge column's own thicknesses (edgeAbove, edgeBottom); the ∇⁴
   * closure on every layer.
   */
  function hybridTendencyLayers(input, out, kFrom, kTo, part = 'all') {
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
        if (k === 0 && mixedDiffusivity > 0) {
          maskedLaplacian(T, lapT, hIn);
          for (let i = 0; i < C; i++) dQ[i] += mixedDiffusivity * lapT[i];
          maskedLaplacian(S, lapT, hIn);
          for (let i = 0; i < C; i++) dW[i] += mixedDiffusivity * lapT[i];
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
        pinned[e] = he <= thinEdge || (he <= tokenEdge[k] && !(bottom > top)) ? 1 : 0;
        if (k === 0) du[e] += stress[e] / rho0 / Math.max(he, minimumLayer);
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
  function hybridBarotropic(dt) {
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
  function barotropic(dt) {
    if (hybrid) { hybridBarotropic(dt); return; }
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

  const bStages = [0, 1, 2, 3].map((s) => [shared(`baroStage${s}Eta`, C), shared(`baroStage${s}U`, E)]);
  const bTrial = [shared('baroTrialEtaA', C), shared('baroTrialUA', E)], bTrialB = [shared('baroTrialEtaB', C), shared('baroTrialUB', E)];
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
    if (hybrid) {
      if (columnRunner) columnRunner('rescale'); else rescaleCells(0, C);
      dried = Atomics.load(counts, 0);
      return;
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
    if (hybrid) {
      Atomics.store(counts, 1, 0);
      if (columnRunner) columnRunner('shift'); else shiftEdges(0, E);
      limited = Atomics.load(counts, 1);
      return;
    }
    edgeThicknesses(h);
    limited = 0;
    for (let e = 0; e < E; e++) {
      if (!edgeOcean[e]) continue;
      let sumH = 0, transport = 0, clamped = false;
      for (let k = 1; k < L; k++) if (hEdge[ae(k, e)] < THIN) u[ae(k, e)] = u[ae(k - 1, e)];
      for (let k = 0; k < L; k++) { sumH += hEdge[ae(k, e)]; transport += hEdge[ae(k, e)] * u[ae(k, e)]; }
      const shift = (avgU[e] - transport) / Math.max(sumH, EPS);
      for (let k = 0; k < L; k++) { const n = ae(k, e); u[n] += shift; if (Math.abs(u[n]) > SPEED_LIMIT) { u[n] = Math.sign(u[n]) * SPEED_LIMIT; clamped = true; } }
      if (clamped) limited++;
    }
  }
  function step(dt) {
    dynamics(dt);
    rescale();
    if (eddyLimit > 0) eddyTransport(dt);
    if (ekmanBlock) ekmanTransport(dt);
    shiftVelocities();
  }
  const phases = { dynamics, rescale, eddy: eddyTransport, ekman: ekmanTransport, shift: shiftVelocities, ...(hybrid ? { regrid, column: columnPhysics } : {}) };

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
    if (hybrid) {
      params[1] = dt;
      if (columnRunner) { columnRunner('eddyEdges'); columnRunner('eddyCells'); columnRunner('eddyApply'); return; }
      eddyEdges(0, E); eddyCells(0, C); eddyApply(0, C);
      return;
    }
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
   * the free surface and shifted to the mean transport, as the slab's
   * are, which lags the barotropic mode by half a step but is the only
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
    const { amount, hOld, hEdgeNew, uCol } = cols;
    for (let e = from; e < to; e++) {
      if (!edgeOcean[e]) continue;
      const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
      for (let k = 0; k < L; k++) {
        hOld[k] = 0.5 * (h[at(k, a)] + h[at(k, b)]);
        hEdgeNew[k] = 0.5 * (targetZ[(k + 1) * C + a] - targetZ[k * C + a] + targetZ[(k + 1) * C + b] - targetZ[k * C + b]);
        uCol[k] = u[ae(k, e)];
      }
      remapColumn(hOld, uCol, L, hEdgeNew, amount, cs, { order: 2 });
      for (let k = 0; k < L; k++) u[ae(k, e)] = hEdgeNew[k] > 0 ? amount[k] / hEdgeNew[k] : uCol[k];
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
    const options = { rho0, gravity: g, mstar, ekman: ekmanFactor, nstar, langmuir, shearViscosity, backgroundViscosity, backgroundDiffusivity, convectiveDiffusivity };
    for (let i = from; i < to; i++) {
      if (!cellOcean[i]) { mixedDepths[i] = 0; boundaryDepths[i] = 0; continue; }
      cellColumn(i, hCol, TCol, SCol);
      for (let k = 0; k < L; k++) cellVectorAt(k, i, vector, 3 * k);
      for (let k = 0; k < L - 1; k++) shear[k] = (vector[3 * k] - vector[3 * k + 3]) ** 2 + (vector[3 * k + 1] - vector[3 * k + 4]) ** 2 + (vector[3 * k + 2] - vector[3 * k + 5]) ** 2;
      options.guess = boundaryDepths[i];
      const result = mixColumn({ h: hCol, T: TCol, S: SCol, shear, K: L, ustar: Math.sqrt(stressAt(i) / rho0), f: fCell[i], dt, options }, out, cs);
      for (let k = 0; k < L; k++) { const n = at(k, i); Q[n] = hCol[k] * TCol[k]; W[n] = hCol[k] * SCol[k]; viscosity[n] = nu[k]; }
      layerDensities(TCol, SCol, L, rhoCol);
      mixedDepths[i] = mixedLayerDepth(hCol, rhoCol, L);
      boundaryDepths[i] = result.depth;
    }
  }
  function columnEdges(from, to) {
    const dt = params[1], { hOld, uCol, nuEdge } = cols;
    for (let e = from; e < to; e++) {
      if (!edgeOcean[e]) continue;
      const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
      for (let k = 0; k < L; k++) { hOld[k] = 0.5 * (h[at(k, a)] + h[at(k, b)]); nuEdge[k] = 0.5 * (viscosity[at(k, a)] + viscosity[at(k, b)]); uCol[k] = u[ae(k, e)]; }
      if (solveColumn(hOld, uCol, L, nuEdge, dt, cs)) for (let k = 0; k < L; k++) u[ae(k, e)] = uCol[k];
    }
  }
  function columnPhysics(dt) {
    params[1] = dt;
    if (columnRunner) { columnRunner('column'); columnRunner('usolve'); return; }
    columnCells(0, C);
    columnEdges(0, E);
  }

  function ekmanTransport(dt) {
    if (!ekmanBlock) return;
    if (!ekman) ekman = { volume: new Float64Array(L * C), heat: new Float64Array(L * C), salt: new Float64Array(L * C), column: new Float64Array(L), shares: new Float64Array(L), block: new Float64Array(L), flux: new Float64Array(L) };
    const { volume, heat, salt: salinity, column, shares, block, flux: classFlux } = ekman;
    const allowance = (k, i) => (Math.max(0, h[at(k, i)] - EPS) + EDDY_SLACK) * areaCell[i] / (nEdgesOnCell[i] * dt);
    volume.fill(0); heat.fill(0); salinity.fill(0);
    ekmanFlux.fill(0); ekmanMoved.fill(0); ekmanLimited = 0;
    for (let e = 0; e < E; e++) {
      if (!edgeOcean[e]) continue;
      const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
      edgeColumn(h, e, column, 0, 1);
      const wind = windDepth ? windDepth[e] : momentumDepth, reach = Math.max(wind, column[0]);
      const rate = interfacialDrag / reach, f = fEdge[e], denominator = rho0 * (rate * rate + f * f);
      if (!(denominator > 0)) continue;
      const transport = (rate * stress[e] + f * tangential(stress, 0, e)) / denominator;
      ekmanFlux[e] = transport;
      if (transport === 0) continue;
      let top = 0, depth = 0;
      shares.fill(0);
      for (let k = 0; k < L; k++) {
        const hk = column[k], part = Math.min(hk, Math.max(0, reach - top));
        if (k === 0 || hk >= stressThickness) { shares[k] = part; depth += part; }
        top += hk;
        if (top >= reach) break;
      }
      if (!(depth > 0)) continue;
      const donor = transport > 0 ? a : b, blockDepth = Math.max(wind, h[donor]);
      let z = 0, total = 0;
      for (let k = 0; k < L; k++) {
        const hk = h[at(k, donor)], water = k === 0 ? hk : Math.max(0, hk - EPS);
        block[k] = Math.min(water, Math.max(0, blockDepth - z));
        total += block[k]; z += hk;
      }
      if (!(total > 0)) continue;
      const volumeFlux = transport * dvEdge[e];
      let scale = 1;
      for (let k = 0; k < L; k++) {
        const change = volumeFlux * (block[k] / total - shares[k] / depth);
        classFlux[k] = change;
        if (change === 0) continue;
        const allowed = allowance(k, change > 0 ? a : b);
        if (Math.abs(change) > allowed) scale = Math.min(scale, allowed / Math.abs(change));
      }
      if (scale < 1) ekmanLimited++;
      for (let k = 0; k < L; k++) {
        const change = scale * classFlux[k];
        if (change === 0) continue;
        const d = at(k, change > 0 ? a : b), na = at(k, a), nb = at(k, b), fq = change * Q[d] / h[d], fw = change * W[d] / h[d];
        volume[na] -= change; volume[nb] += change;
        heat[na] -= fq; heat[nb] += fq;
        salinity[na] -= fw; salinity[nb] += fw;
        ekmanMoved[k] += Math.abs(change) * 86400 / 1e9;
      }
    }
    for (let n = 0; n < L * C; n++) {
      const factor = dt / areaCell[n % C];
      h[n] += factor * volume[n]; Q[n] += factor * heat[n]; W[n] += factor * salinity[n];
    }
  }

  function move(i, from, to, amount, process = 0) {
    const a = at(from, i), b = at(to, i);
    const f = amount / h[a];
    const dQ = Q[a] * f, dW = W[a] * f;
    h[a] -= amount; Q[a] -= dQ; W[a] -= dW;
    h[b] += amount; Q[b] += dQ; W[b] += dW;
    if (intake && to === 0 && intake.slot[i] >= 0) {
      const n = (process * L + from) * intake.cells.length + intake.slot[i];
      intake.water[n] += amount; intake.degrees[n] += dQ;
    }
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

  function baseShears() {
    for (let e = 0; e < E; e++) {
      let j = 0;
      if (edgeOcean[e]) for (let k = 1; k < L; k++) if (hEdge[ae(k, e)] >= THIN) { j = k; break; }
      baseShear[e] = j > 0 ? u[e] - u[ae(j, e)] : 0;
    }
    cellVector(mesh, baseShear, shearCell);
  }

  function mixedLayer(dt, closing = 0) {
    if (carryMomentum) edgeThicknesses(h);
    cellVector(mesh, stress, tauCell);
    if (shearEntrainment > 0) baseShears();
    if (mixedNeighbourRatio > 0) neighbourReach();
    let left = 0, blockWind = 0;
    const fromSource = (i, want) => {
      if (!(left > 0 && want > 0)) return 0;
      const depth = Math.max(blockWind, h[i]);
      let top = h[i], source = -1;
      for (let k = 1; k < L; k++) {
        const hk = h[at(k, i)];
        if (hk > THIN && top + hk > depth) { source = k; break; }
        top += hk;
      }
      if (source < 0) return 0;
      const part = Math.min(want, left, h[at(source, i)] - EPS);
      if (!(part > 0)) return 0;
      move(i, source, 0, part, INTAKE_UPWELLED);
      left -= part; upwelled[i] += part;
      return part;
    };
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) continue;
      left = 0;
      if (coastal && coastal[i]) {
        let divergence = 0, wind = 0, sea = 0;
        for (let m = 0; m < nEdgesOnCell[i]; m++) {
          const e = edgesOnCell[maxEdges * i + m];
          divergence += edgeSignOnCell[maxEdges * i + m] * ekmanFlux[e] * dvEdge[e];
          if (edgeOcean[e]) { wind += windDepth ? windDepth[e] : momentumDepth; sea++; }
        }
        if (sea > 0) { left = upwellingSource * Math.max(0, divergence / areaCell[i]) * dt; blockWind = wind / sea; }
      }
      const tau = Math.hypot(tauCell[3 * i], tauCell[3 * i + 1], tauCell[3 * i + 2]);
      const ustar3 = Math.pow(tau / rho0, 1.5), stir = stirring * Math.exp(-h[i] / stirringDepth), began = h[i];
      const deepest = Math.min(maximumMixedDepth, reach[i]), s0 = W[i] / h[i];
      const salted = s0 * fresh[i] / 1000 + (s0 - iceSalinity) * (iceVolume[i] - previousVolume[i]) * iceDensity / 1000;
      const buoyancy = g * (expansionOf(surfaceIn[i], s0) * (previousT0[i] - surfaceIn[i]) * capacity[i] / rhoCp + contractionOf(surfaceIn[i], s0) * salted) / dt;
      buoyancyLoss[i] += (buoyancy - buoyancyLoss[i]) * (buoyancyMemory > 0 ? Math.min(1, dt / buoyancyMemory) : 1);
      const lightBuoyancy = daily ? g * expansionOf(surfaceIn[i], s0) / rhoCp : 0;
      if (daily) {
        daily.sunSum[i] += daily.sunStep[i];
        daily.darkSum[i] += buoyancy * dt + lightBuoyancy * daily.sunStep[i] * (1 - transmission[i]);
        if (closing !== 0) {
          if (closing > 0) { daily.sunDay[i] = daily.sunSum[i] * closing; daily.darkDay[i] = daily.darkSum[i] * closing; }
          daily.sunSum[i] = 0; daily.darkSum[i] = 0;
        }
      }
      let rm = eos(Q[i] / h[i], W[i] / h[i]), budget = convectiveRate * dt;
      for (let k = 1; k < L && h[i] < deepest && budget > 0; k++) {
        const n = at(k, i);
        if (h[n] <= EPS) continue;
        let take = Math.min(h[n] - EPS, deepest - h[i], budget);
        const eroding = convectiveErosion && rm < densest[k];
        if (eroding) take = rm >= lightest[k] && buoyancyLoss[i] > 0 ? Math.min(take, buoyancyLoss[i] * dt * rho0 * h[n] / (g * h[i] * (densest[k] - rm))) : 0;
        else if (!convectiveErosion && rho[k] > rm) continue;
        if (take > 0) { move(i, k, 0, take, INTAKE_CONVECTIVE); budget -= take; rm = eos(Q[i] / h[i], W[i] / h[i]); }
        if (eroding) break;
      }
      let below = -1, hold = 0;
      for (let k = 1; k < L; k++) if (h[at(k, i)] > THIN) { below = k; break; }
      if (shearEntrainment > 0 && below > 0) {
        const h0 = h[i], shear = shearCell[3 * i] ** 2 + shearCell[3 * i + 1] ** 2 + shearCell[3 * i + 2] ** 2;
        const db = Math.max(g * (rho[below] - rm) / rho0, g * DENSITY_TOLERANCE / rho0), ri = db * h0 / (shear + 1e-12);
        if (ri < shearEntrainment) {
          const reached = h0 * Math.sqrt(shearEntrainment / ri), want = Math.min(reached - h0, deepest - h0);
          const drawn = fromSource(i, want), take = Math.min(want - drawn, h[at(below, i)] - EPS);
          if (take > 0) move(i, below, 0, take, INTAKE_SHEAR);
          if (drawn + take > 0) { sheared[i] += drawn + take; rm = eos(Q[i] / h[i], W[i] / h[i]); }
          hold = Math.min(reached, h[i]);
          below = -1;
          for (let k = 1; k < L; k++) if (h[at(k, i)] > THIN) { below = k; break; }
        }
      }
      if (below > 0 && h[i] < deepest) {
        const db = Math.max(1e-3, g * (rho[below] - rm) / rho0);
        const want = Math.min(2 * stir * ustar3 / (h[i] * db) * dt, reach[i] - h[i]);
        const drawn = fromSource(i, want), entrain = Math.min(want - drawn, h[at(below, i)] - EPS);
        if (entrain > 0) move(i, below, 0, entrain, INTAKE_STIRRING);
        if (drawn + entrain > 0) rm = eos(Q[i] / h[i], W[i] / h[i]);
      }
      below = -1;
      for (let k = 1; k < L; k++) if (h[at(k, i)] > THIN) { below = k; break; }
      let excess = Math.max(0, h[i] - maximumMixedDepth);
      if (neutralSnap && below > 0 && rm >= rho[below] - DENSITY_TOLERANCE) excess = Math.max(excess, h[i] - shallowestMixedDepth);
      detrain(i, excess, rm);
      if (h[i] > reach[i]) detrain(i, (h[i] - Math.max(reach[i], shallowestMixedDepth)) * Math.min(1, dt / detrainmentTime), rm);
      let drive = buoyancyLoss[i], gate = drive;
      if (daily) {
        const passing = penetration ? transmissionAt(i) : 0, sun = lightBuoyancy * daily.sunDay[i];
        drive = daily.darkDay[i] - sun * (1 - passing);
        gate = retreatLight > 0 && passing > 0 ? drive + retreatLight * sun * h[i] * deposited(h[i], D[i] + eta[i], penetration) : drive;
      }
      retreatDepth[i] = 0;
      if (drive < -1e-9 && gate < -1e-9) {
        retreatDepth[i] = retreatStirring > 0 || langmuir > 0
          ? (retreatStirring > 0 ? retreatStirring : 2 * stir) * (langmuir > 0 ? 1 + (langmuirFactor(Math.sqrt(tau / rho0), began, langmuir) - 1) * (1 - iceCover[i]) : 1) * ustar3 / -drive
          : 2 * stir * ustar3 / -drive;
        if (retreatEkman > 0 && fCell[i] !== 0) retreatDepth[i] = Math.min(retreatDepth[i], retreatEkman * Math.sqrt(tau / rho0) / Math.abs(fCell[i]));
        const monin = Math.max(shallowestMixedDepth, retreatDepth[i], hold);
        if (h[i] > monin) detrain(i, (h[i] - monin) * Math.min(1, dt / detrainmentTime), rm);
      }
      if (h[i] < minimumThickness) {
        const h0 = h[i], q0 = Q[i];
        const refillTo = refillRate > 0 ? Math.max(Math.min(minimumThickness, h0 + refillRate * dt), Math.min(minimumThickness, refillFloor)) : minimumThickness;
        fromSource(i, refillTo - h[i]);
        for (let k = 1; k < L && h[i] < refillTo; k++) {
          const available = h[at(k, i)] - EPS;
          if (available > 0) move(i, k, 0, Math.min(available, refillTo - h[i]), INTAKE_REFILL);
        }
        refilled.water[i] += h[i] - h0;
        refilled.heat[i] += rhoCp * (Q[i] - q0 - (h[i] - h0) * q0 / h0);
      }
      restoreDensities(i, Math.min(1, dt / restoreTime));
    }
    if (carryMomentum) exchangeMomentum();
  }

  function exchangeMomentum() {
    for (let e = 0; e < E; e++) {
      if (!edgeOcean[e]) continue;
      edgeColumn(h, e, column, 0, 1);
      let carried = 0;
      for (let k = 0; k < L; k++) carried += hEdge[ae(k, e)] * u[ae(k, e)];
      const reach = Math.max(hEdge[e], column[0]);
      let top = 0, within = 0, mixed = 0;
      for (let k = 0; k < L && top < reach; k++) {
        const old = hEdge[ae(k, e)], part = Math.min(old, reach - top);
        mixed += part * u[ae(k, e)]; within += part; top += old;
      }
      const uMixed = within > 0 ? mixed / within : u[e];
      let total = 0, held = 0;
      top = 0;
      for (let k = 0; k < L; k++) {
        const n = ae(k, e), hk = column[k], inside = Math.max(0, Math.min(hk, reach - top));
        if (k > 0 && hk < THIN) u[n] = u[n - E];
        else if (hk > 0) u[n] = (inside * uMixed + (hk - inside) * u[n]) / hk;
        hEdge[n] = hk; total += hk; held += hk * u[n]; top += hk;
      }
      const shift = (carried - held) / Math.max(total, EPS);
      let clamped = false;
      for (let k = 0; k < L; k++) { const n = ae(k, e); u[n] += shift; if (Math.abs(u[n]) > SPEED_LIMIT) { u[n] = Math.sign(u[n]) * SPEED_LIMIT; clamped = true; } }
      if (clamped) limited++;
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

  const coverOf = (ice, concentration, i) => (ice[i] > 0 ? (concentration && concentration[i] > 0 ? concentration[i] : 1) : 0);
  function setStress(total, ice, concentration = null) {
    const cover = (i) => coverOf(ice, concentration, i);
    for (let e = 0; e < E; e++) stress[e] = edgeOcean[e] ? total[e] * (1 - 0.5 * (cover(cellsOnEdge[2 * e]) + cover(cellsOnEdge[2 * e + 1])) * (1 - iceStressTransmission)) : 0;
    if (windDepth) windDepths();
  }
  function measureIce(ice, concentration, out) {
    for (let i = 0; i < C; i++) out[i] = coverOf(ice, concentration, i) * ice[i];
  }
  function readSurface(surfaceT, ice, concentration = seaCover) {
    for (let i = 0; i < C; i++) {
      iced[i] = ice[i] > 0 ? 1 : 0;
      if (iceCover) iceCover[i] = coverOf(ice, concentration, i);
      T0[i] = iced[i] ? FREEZING_POINT : surfaceT[i];
      surfaceIn[i] = T0[i];
      Q[i] = h[i] * T0[i];
      if (daily) daily.sunStep[i] = transmission[i] > 1e-6 ? penetrated[i] / transmission[i] : 0;
      if (penetrated[i] !== 0) depositLight(i);
    }
  }
  function depositLight(i) {
    const light = penetrated[i];
    penetrated[i] = 0;
    if (!(transmission[i] > 0)) { Q[i] += light / rhoCp; return; }
    if (hybrid) {
      let z = h[i], top = transmitted(z, penetration), given = 0, last = 1;
      for (let k = 1; k < L; k++) { z += h[at(k, i)]; const bottom = transmitted(z, penetration), part = (top - bottom) / transmission[i]; shares[k] = part; given += part; top = bottom; last = k; }
      shares[last] += 1 - given;
      for (let k = 1; k < L; k++) Q[at(k, i)] += light * shares[k] / rhoCp;
      return;
    }
    for (let k = 0; k < L; k++) column[k] = h[at(k, i)];
    penetrationShares(column, L, D[i] + eta[i], penetration, shares);
    for (let k = 1; k < L; k++) if (shares[k] !== 0) Q[at(k, i)] += light * shares[k] / rhoCp;
  }
  function transmissionAt(i) {
    if (hybrid) return penetration ? transmitted(h[i], penetration) : 0;
    for (let k = 0; k < L; k++) column[k] = h[at(k, i)];
    return penetrationShares(column, L, D[i] + eta[i], penetration);
  }
  function publishTransmission() {
    transmission.fill(0);
    if (penetration) for (let i = 0; i < C; i++) if (cellOcean[i]) transmission[i] = transmissionAt(i);
  }
  function accumulate(evaporation, rain, dt, runoff = null) {
    for (let i = 0; i < C; i++) if (cellOcean[i]) fresh[i] += (evaporation ? evaporation[i] * dt : 0) - (rain ? rain[i] : 0);
    if (runoff) for (let j = 0; j < C; j++) if (!cellOcean[j] && runoff[j] > 0 && outlet[j] >= 0) fresh[outlet[j]] -= runoff[j] * areaCell[j] / areaCell[outlet[j]];
  }
  function salt() {
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) continue;
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
    const dtOcean = everySteps * dt, closing = daily ? closeDay(day, time, dtOcean) : 0;
    readSurface(surfaceT, ice, concentration);
    measureIce(ice, concentration, iceVolume);
    setStress(typeof totalStress === 'function' ? totalStress() : totalStress, ice, concentration);
    step(dtOcean);
    if (hybrid) { regrid(dtOcean); columnPhysics(dtOcean); } else mixedLayer(dtOcean, closing);
    salt();
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
  /*
   * The hybrid ocean's start: every sea column on the depth table clipped
   * to its floor (the z grid of a column held at its minima), with the
   * atlas's layer means of temperature and salinity where the atlas
   * covers the cell, the nearest covered column's where it does not, and
   * without an atlas the analytic profile T = T_d + (T₀ − T_d) e^(−z/λ),
   * T_d 2 °C, λ 400 m stretched by the thermocline tilt as the slab's
   * classes are, the salinity from salinityProfile to 34.7 the same way;
   * the steric free surface; then one regrid, so the first state is
   * already on the generator's targets.
   */
  function hybridGrid(i, zCol) {
    const H = D[i] + eta[i];
    for (let k = 1; k < L; k++) zCol[k] = Math.min(table.depths[k], H);
    spaceInterfaces(zCol, L, H, minimumLayer);
    return H;
  }
  function buildHybrid(surfaceT, ice, atlas, concentration) {
    u.fill(0); eta.fill(0); fresh.fill(0); buoyancyLoss.fill(0);
    forgetDays();
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
      iced[i] = ice[i] > 0 ? 1 : 0;
      if (!cellOcean[i]) { for (let k = 0; k < L; k++) { h[at(k, i)] = 0; Q[at(k, i)] = 0; W[at(k, i)] = 0; } T0[i] = iced[i] ? FREEZING_POINT : surfaceT[i]; S0[i] = referenceS; continue; }
      hybridGrid(i, zCol);
      let column = atlas ? atlas.columnAt(latCell[i], mesh.lonCell[i]) : null;
      if (column) {
        const reached = column.depths.length;
        if (atlas.columnReaching && atlas.depths && reached < atlas.depths.length && atlas.depths[reached] <= D[i]) {
          const deep = atlas.columnReaching(latCell[i], mesh.lonCell[i], D[i]);
          if (deep && deep.depths.length > reached) {
            const Td = Float64Array.from(deep.T), Sd = Float64Array.from(deep.S);
            Td.set(column.T); Sd.set(column.S);
            column = { depths: deep.depths, T: Td, S: Sd };
          }
        }
        filled[i] = 1;
        zs.length = ts.length = ss.length = 0;
        for (let j = 0; j < column.depths.length && column.depths[j] < D[i]; j++) { zs.push(column.depths[j]); ts.push(column.T[j]); ss.push(column.S[j]); }
        const [tb, sb] = profileAt(column, D[i]);
        zs.push(D[i]); ts.push(tb); ss.push(sb);
        for (let k = 0; k < L; k++) { const n = at(k, i), hk = zCol[k + 1] - zCol[k]; h[n] = hk; Q[n] = integral(ts, zCol[k], zCol[k + 1]); W[n] = integral(ss, zCol[k], zCol[k + 1]); }
      } else if (atlas) {
        pending.push(i);
        for (let k = 0; k < L; k++) h[at(k, i)] = zCol[k + 1] - zCol[k];
      } else {
        const lat = latCell[i], s0 = salinityProfile(lat), t0 = iced[i] ? FREEZING_POINT : surfaceT[i], scale = 400 * (1 - thermoclineTilt + 2 * thermoclineTilt * Math.cos(lat) ** 2);
        const tDeep = 275.15, sDeep = 34.7, decay = (a, b) => scale * (Math.exp(-a / scale) - Math.exp(-b / scale)) / (b - a);
        for (let k = 0; k < L; k++) {
          const n = at(k, i), hk = zCol[k + 1] - zCol[k], f = decay(zCol[k], zCol[k + 1]);
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
      let z = 0;
      for (let k = 0; k < L; k++) {
        const n = at(k, i), hk = h[n], center = z + 0.5 * hk;
        z += hk;
        let t = iced[i] ? FREEZING_POINT : surfaceT[i], s = referenceS;
        if (source >= 0) {
          let top = 0, j = 0;
          while (j < L - 1 && top + h[at(j, source)] < center) { top += h[at(j, source)]; j++; }
          const m = at(j, source);
          if (h[m] > 1e-9) { t = Q[m] / h[m]; s = W[m] / h[m]; }
        }
        Q[n] = hk * t; W[n] = hk * s;
      }
    }
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) continue;
      T0[i] = iced[i] ? FREEZING_POINT : Math.max(FREEZING_POINT, Q[i] / h[i]);
      Q[i] = h[i] * T0[i];
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
  function build(surfaceT, ice, atlas, concentration = seaCover) {
    if (hybrid) return buildHybrid(surfaceT, ice, atlas, concentration);
    u.fill(0); eta.fill(0); fresh.fill(0); buoyancyLoss.fill(0);
    forgetDays();
    measureIce(ice, concentration, previousVolume);
    const filled = atlas ? atlasColumns(mesh, atlas, { D, cellOcean, ice, rho, labelT, labelS, h, Q, W, T0, shallowestMixedDepth, maximumMixedDepth }) : null;
    for (let i = 0; i < C; i++) {
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
    publishTransmission();
    penetrated.fill(0);
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
        anomaly += ((hybrid ? layerRho[at(k, i)] : k === 0 ? rhoMl[i] : rho[k]) - rho0) * part;
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
   * before its columns are fitted to the bathymetry. An open-water cell
   * filled from its neighbors puts its mixed layer's temperature into
   * `surfaceT`. As on the GPU (layeredOcean.gpu.js, restartable), a state
   * that carries every array of the GPU's restartArrays, its columns within
   * a meter of this bathymetry and dry on land, brings back the mixed
   * layer's last temperature, the ice volume it last saw, the freshwater
   * not yet taken, the remembered buoyancy loss and the heat capacity,
   * and the sunlight not yet deposited when it carries `penetrated`.
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
   * A saved hybrid ocean on these labels loads as it is, its columns
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
  function loadHybrid(saved, surfaceT, ice, concentration) {
    if (!saved.h) { build(surfaceT, ice, climatology, concentration); return; }
    surfaceTAt = (i) => (ice[i] > 0 ? FREEZING_POINT : surfaceT[i]);
    const classes = savedDensities(saved, C), own = classes && sameDensities(classes, labels) && saved.h.length === L * C;
    if (own) {
      h.set(saved.h); u.set(saved.u); eta.set(saved.eta);
      for (let n = 0; n < L * C; n++) { Q[n] = h[n] * saved.T[n]; W[n] = h[n] * saved.S[n]; }
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
    if (!own) regrid(Infinity);
    measureIce(ice, concentration, previousVolume);
    for (let i = 0; i < C; i++) { iced[i] = ice[i] > 0 ? 1 : 0; T0[i] = cellOcean[i] ? Q[i] / h[i] : surfaceTAt(i); S0[i] = cellOcean[i] ? W[i] / h[i] : referenceS; previousT0[i] = T0[i]; capacity[i] = rhoCp * Math.max(h[i], 1); }
    fresh.fill(0); buoyancyLoss.fill(0); penetrated.fill(0);
    forgetDays();
    layerDensitiesAll();
    mixedDepthsAll();
    publishTransmission();
    if (restartable(saved)) {
      previousT0.set(saved.previousT0); previousVolume.set(saved.previousIceVolume); fresh.set(saved.fresh); capacity.set(saved.capacity);
      if (saved.penetrated && saved.penetrated.length === C) penetrated.set(saved.penetrated);
    }
    counter = 0;
    initialised = true;
  }
  function load(saved, surfaceT, ice, concentration = seaCover) {
    if (hybrid) return loadHybrid(saved, surfaceT, ice, concentration);
    if (saved && saved.vertical === 'hybrid') throw new Error("a hybrid ocean state loads into the ocean with vertical: 'hybrid'");
    const classes = savedDensities(saved, C);
    if (classes && !sameDensities(classes, densities)) saved = rebinOcean(saved, classes, densities, C, { labelT, labelS });
    if (!saved.h || saved.h.length !== L * C) {
      build(surfaceT, ice, climatology, concentration);
      return;
    }
    h.set(saved.h); u.set(saved.u); eta.set(saved.eta);
    for (let n = 0; n < L * C; n++) { Q[n] = h[n] * saved.T[n]; W[n] = h[n] * saved.S[n]; }
    const start = () => {
      const kept = { h: Float64Array.from(h), u: Float64Array.from(u), Q: Float64Array.from(Q), W: Float64Array.from(W), eta: Float64Array.from(eta) };
      build(surfaceT, ice, climatology, concentration);
      const built = { h: Float64Array.from(h), Q: Float64Array.from(Q), W: Float64Array.from(W) };
      h.set(kept.h); u.set(kept.u); Q.set(kept.Q); W.set(kept.W); eta.set(kept.eta);
      return built;
    };
    const { fromNeighbors } = fitColumns({ h, Q, W, eta }, start, { D, cellOcean, L, C, labelT, labelS, minimumThickness, mesh });
    for (let k = 0; k < L; k++) for (let e = 0; e < E; e++) if (!edgeOcean[e]) u[ae(k, e)] = 0;
    measureIce(ice, concentration, previousVolume);
    for (let i = 0; i < C; i++) { iced[i] = ice[i] > 0 ? 1 : 0; T0[i] = Q[i] / Math.max(EPS, h[i]); S0[i] = W[i] / Math.max(EPS, h[i]); previousT0[i] = T0[i]; capacity[i] = rhoCp * Math.max(h[i], 1); }
    for (const i of fromNeighbors) if (!iced[i]) surfaceT[i] = T0[i];
    fresh.fill(0); buoyancyLoss.fill(0); penetrated.fill(0);
    forgetDays();
    publishTransmission();
    if (restartable(saved)) {
      previousT0.set(saved.previousT0); previousVolume.set(saved.previousIceVolume); fresh.set(saved.fresh); buoyancyLoss.set(saved.buoyancy); capacity.set(saved.capacity);
      if (saved.penetrated && saved.penetrated.length === C) penetrated.set(saved.penetrated);
      if (daily && DAILY.every((name) => saved[name] && saved[name].length === C)) {
        for (const name of DAILY) daily[name].set(saved[name]);
        day.span = saved.daySpan ? saved.daySpan[0] : 0;
      }
    }
    counter = 0;
    initialised = true;
  }

  function forgetDays() {
    if (daily) for (const values of Object.values(daily)) values.fill(0);
    day.span = 0;
  }
  function restartArrays() {
    const copy = (values) => Float64Array.from(values);
    return {
      previousT0: copy(previousT0), previousIceVolume: copy(previousVolume), fresh: copy(fresh), buoyancy: copy(buoyancyLoss), capacity: copy(capacity), penetrated: copy(penetrated),
      ...(penetration ? { transmission: copy(transmission) } : {}),
      ...(daily ? { ...Object.fromEntries(DAILY.map((name) => [name, copy(daily[name])])), daySpan: Float64Array.of(day.span) } : {}),
    };
  }

  function serialize() {
    const Tall = new Array(L * C), Sall = new Array(L * C);
    for (let n = 0; n < L * C; n++) { const hh = Math.max(EPS, h[n]); Tall[n] = Q[n] / hh; Sall[n] = W[n] / hh; }
    return { h: Array.from(h), u: Array.from(u), T: Tall, S: Sall, eta: Array.from(eta), densities: Array.from(labels), ...(hybrid ? { vertical, generator } : {}) };
  }

  const thermoclineDepth = new Float64Array(C), sst = T0, sss = S0;
  function thermocline(i) {
    if (hybrid) {
      let z = 0, previous = -Infinity, center = 0;
      for (let k = 0; k < L; k++) {
        const n = at(k, i), hv = h[n], rk = Math.max(previous, layerRho[n]), c = z + 0.5 * hv;
        if (rk >= THERMOCLINE_DENSITY) return k === 0 ? 0 : center + (THERMOCLINE_DENSITY - previous) / Math.max(rk - previous, 1e-12) * (c - center);
        previous = rk; center = c; z += hv;
      }
      return z;
    }
    let depth = 0;
    for (let k = 0; k <= thermoclineLayers; k++) depth += h[at(k, i)];
    return depth;
  }
  function fields(depth = 0) {
    if (hybrid) layerDensitiesAll();
    for (let i = 0; i < C; i++) thermoclineDepth[i] = cellOcean[i] ? thermocline(i) : NaN;
    const layerT = (k, i) => (k === 0 ? T0[i] : Q[at(k, i)] / h[at(k, i)]);
    return { h1: hybrid ? mixedDepths : h.subarray(0, C), T1: T0, S1: S0, u1: u.subarray(0, E), eta, thermoclineDepth, ...depthFields(mesh, L, { h, u, temperature: layerT, cellOcean }, depth, { thin: hybrid ? 0 : THIN }) };
  }

  function diagnostics() {
    if (hybrid) layerDensitiesAll();
    let area = 0, depth = 0, heat = 0, thermo = 0, salinity = 0, ssh = 0, speed = 0, transport = 0, interiorT = 0, interiorH = 0, thinned = 0;
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) continue;
      const a = areaCell[i];
      area += a; depth += a * (hybrid ? mixedDepths[i] : h[i]); salinity += a * (W[i] / h[i]); ssh = Math.max(ssh, Math.abs(eta[i]));
      let below = 0;
      if (!hybrid) for (let k = 1; k < L; k++) below += h[at(k, i)] - EPS;
      if (!hybrid && h[i] < Math.min(minimumThickness, h[i] + below) - THINNED_MARGIN) thinned += a;
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
    return { oceanUpperDepth: depth / area, oceanHeat: heat / area, oceanInteriorT: interiorH > 0 ? interiorT / interiorH : 0, oceanSpeed: speed, oceanThermoclineDepth: thermo / area, oceanSalinity: salinity / area, oceanSSH: ssh, oceanTransport: transport / 1e6, oceanLimited: limited, oceanThinned: thinned / area, ...(ekmanBlock ? { oceanEkmanLimited: ekmanLimited } : {}), ...(hybrid ? { oceanDried: dried } : {}) };
  }

  if (!adopting) build(new Float64Array(C).fill(288), new Float64Array(C), null);
  initialised = false;
  const sharedBuffers = { h: h.buffer, u: u.buffer, Q: Q.buffer, W: W.buffer, eta: eta.buffer, rhoMl: rhoMl.buffer, capacity: capacity.buffer, transmission: transmission.buffer, penetrated: penetrated.buffer, stress: stress.buffer, hEdge: hEdge.buffer, stressDepth: stressDepth.buffer, stressLayer: stressLayer.buffer, stressFraction: stressFraction.buffer, ...(windDepth ? { windDepth: windDepth.buffer } : {}), deepestEdge: deepestEdge.buffer, pressure: pressure.buffer, gradEta: gradEta.buffer, gradRho: gradRho.buffer, params: params.buffer, trialh: trial[0].buffer, trialu: trial[1].buffer, trialQ: trial[2].buffer, trialW: trial[3].buffer };
  Object.assign(sharedBuffers, { baroSlow: slow.buffer, baroU: U.buffer, baroDepth: depthEdge.buffer, baroEta: etaB.buffer, baroAvgU: avgU.buffer, baroAvgEta: avgEta.buffer, baroTrialEtaA: bTrial[0].buffer, baroTrialUA: bTrial[1].buffer, baroTrialEtaB: bTrialB[0].buffer, baroTrialUB: bTrialB[1].buffer, counts: countBuffer });
  bStages.forEach((stage, s) => { sharedBuffers[`baroStage${s}Eta`] = stage[0].buffer; sharedBuffers[`baroStage${s}U`] = stage[1].buffer; });
  if (hybrid) Object.assign(sharedBuffers, { layerRho: layerRho.buffer, columnAbove: columnAbove.buffer, targetZ: targetZ.buffer, viscosity: viscosity.buffer, edgeBottom: edgeBottom.buffer, edgeAbove: edgeAbove.buffer, mixedDepths: mixedDepths.buffer, boundaryDepths: boundaryDepths.buffer, eddyFlux: eddyFlux.buffer, eddyVolume: eddyVolume.buffer, eddyHeat: eddyHeat.buffer, eddySalt: eddySalt.buffer });
  stages.forEach((stage, s) => { sharedBuffers[`stage${s}h`] = stage[0].buffer; sharedBuffers[`stage${s}u`] = stage[1].buffer; sharedBuffers[`stage${s}Q`] = stage[2].buffer; sharedBuffers[`stage${s}W`] = stage[3].buffer; });
  return { state, trial, stages, layers: L, h, u, Q, W, eta, D, sillEdge, hEdge, T0, S0, capacity, transmission, penetrated, penetration, stress, fresh, refilled, sheared, upwelled, intake, retreatDepth, windDepth, avgU, tendency, prepare, tendencyLayers: hybrid ? hybridTendencyLayers : tendencyLayers, interfaceRate: (e, up, down, input = state) => interfaceRate(input[1], e, up, down), windAcceleration, setLayerRunner(fn) { layerRunner = fn; }, setColumnRunner(fn) { columnRunner = fn; }, useSeaIce(concentration) { seaCover = concentration; }, advance, accumulate, initialize, load, serialize, restartArrays, daily, diagnostics, fields, setStress, readSurface, mixedLayer, salt, writeSurface, phases, eddyTransport, ekmanTransport, ekmanFlux, ekmanMoved, rhoCp, edgeOcean, cellOcean, densities: rho, shared: sharedBuffers, vertical, ...(hybrid ? { generator, table, mixedDepths, boundaryDepths, viscosity, layerRho, regrid, columnPhysics, hybridTendencyLayers, regridTargets, regridEdges, regridCells, columnCells, columnEdges, prepareCells, prepareEdges, rescaleCells, shiftEdges, eddyEdges, eddyCells, eddyApply, baroPrep, baroCoriolis, baroEdges, baroCells } : {}) };
}
