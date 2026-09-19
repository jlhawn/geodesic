import { buildMesh } from './mesh.module.js';
import { createSigmaCore } from './dynamics/sigmaCore.module.js';
import { createRK4Arrays } from './dynamics/integrators.module.js';
import { createRadiation } from './physics/radiation.module.js';
import { createSurface } from './physics/surface.module.js';
import { createMoistPhysics } from './physics/moist.module.js';

export const SIDEREAL_DAY = 86164.0905;
export const STATE_NAMES = ['pi', 'theta', 'u', 'surfaceT', 'q', 'qc'];
export const stateLengths = ({ K, C, E }) => ({ pi: C, theta: K * C, u: K * E, surfaceT: C, q: K * C, qc: K * C });

/*
 * The climate model: the sigma-coordinate core with gray radiation, a
 * slab-ocean surface, bulk surface fluxes of heat and moisture,
 * boundary-layer drag, large-scale condensation, Betts–Miller and dry
 * convective adjustment, and a 5-day Rayleigh sponge above σ = 0.02 in
 * the role of gravity-wave drag on the polar-night jet. State is
 * [pi, theta, u, surfaceT, q, qc] (vapour and cloud condensate), each
 * on a SharedArrayBuffer; with `moist: false` q and qc are carried but
 * never sourced, so they stay zero.
 *
 * One time-step tendency is four phases; each takes an index range so
 * the same code runs whole on one thread or sliced across workers:
 *   flux(layers)   mass flux and divergence
 *   column(cells)  dπ/dt, σ̇, Exner and geopotential, lowest-layer wind
 *   vertex(verts)  kite-weighted π
 *   layer(layers)  θ and momentum tendencies, closures, drag
 *   cell(cells)    radiation and the slab surface
 * Arrays read across phases live in `shared`; `buffers` adopts another
 * instance's so a worker computes on the same memory.
 */
export function createModel(gridOrMesh, {
  radius, core: coreOptions = {}, radiation: radiationOptions = {}, surface: surfaceOptions = {}, moist: moistOptions = {},
  physics = true, moist = true, nu4Hours = 3, buffers = null,
} = {}) {
  const mesh = gridOrMesh.nCells ? gridOrMesh : buildMesh(gridOrMesh, { radius, omega: 2 * Math.PI / SIDEREAL_DAY });
  let spacing = 0;
  for (let e = 0; e < mesh.nEdges; e++) spacing += mesh.dcEdge[e];
  spacing /= mesh.nEdges;
  const nu4 = Math.pow(spacing / Math.PI, 4) / (nu4Hours * 3600);
  const core = createSigmaCore(mesh, { nu4, nu4Theta: nu4, buffers: buffers ? buffers.core : null, ...coreOptions });
  const { K, C, E, V } = core.diagnostics;
  const radiation = createRadiation(mesh, core, radiationOptions);
  const surface = createSurface(mesh, core, { topSigma: 0.02, topDragDays: 5, buffers: buffers ? buffers.surface : null, ...surfaceOptions });
  const moistPhysics = createMoistPhysics(mesh, core, { buffers: buffers ? buffers.moist : null, ...moistOptions });
  const totals = { absorbedSolar: 0, outgoingLongwave: 0, sensibleHeat: 0, evaporation: 0 };

  const lengths = stateLengths({ K, C, E });
  const stateArray = (name) => new Float64Array(buffers && buffers.state && buffers.state[name] ? buffers.state[name] : new SharedArrayBuffer(8 * lengths[name]));
  const state = STATE_NAMES.map(stateArray);

  const phases = {
    flux(input, kFrom, kTo) { core.phaseFlux(input, kFrom, kTo); },
    column(input, out, iFrom, iTo) {
      core.phaseColumn(input, out, iFrom, iTo);
      if (physics) surface.lowestWindSpeed(input[2], iFrom, iTo);
    },
    vertex(input, vFrom, vTo) { core.phaseVertex(input, vFrom, vTo); },
    layer(input, out, kFrom, kTo) {
      core.phaseLayer(input, out, kFrom, kTo);
      if (physics) surface.applyLayers(input, out, kFrom, kTo);
    },
    cell(input, out, iFrom, iTo, sums) {
      out[3].fill(0, iFrom, iTo);
      if (physics) radiation.apply(moist ? input : input.slice(0, 4), out, surface.windSpeed, sums, iFrom, iTo);
    },
    adjust(iFrom, iTo, dt) {
      if (!physics) return;
      if (moist) moistPhysics.adjust(state, iFrom, iTo, dt);
      surface.convectiveAdjustment(state[0], state[1], iFrom, iTo, moist ? state[4] : null, moist ? state[5] : null);
    },
  };

  function tendency(input, out) {
    phases.flux(input, 0, K);
    phases.column(input, out, 0, C);
    phases.vertex(input, 0, V);
    phases.layer(input, out, 0, K);
    phases.cell(input, out, 0, C, totals);
  }

  let rk4 = null;
  const model = {
    mesh, core, radiation, surface, moist: moistPhysics, state, totals, phases, tendency, physics, moistOn: physics && moist, time: 0,
    shared: { core: core.shared, surface: surface.shared, moist: moistPhysics.shared, state: Object.fromEntries(STATE_NAMES.map((name, a) => [name, state[a].buffer])) },
  };

  model.step = function step(dt) {
    rk4 ??= createRK4Arrays(STATE_NAMES.map((name) => lengths[name]));
    radiation.setTime(model.time);
    rk4(tendency, state, dt);
    phases.adjust(0, C, dt);
    model.time += dt;
  };

  let lastPrecipTime = 0;
  model.diagnostics = function diagnostics(sums = totals) {
    const [pi, theta, u, surfaceT, q, qc] = state;
    const precipitation = moistPhysics.precipitation;
    let area = 0, mass = 0, meanSurfaceT = 0, piMin = Infinity, piMax = -Infinity, maxWind = 0, water = 0, cloud = 0, rain = 0;
    for (let i = 0; i < C; i++) {
      const a = mesh.areaCell[i];
      area += a;
      mass += a * pi[i];
      meanSurfaceT += a * surfaceT[i];
      piMin = Math.min(piMin, pi[i]);
      piMax = Math.max(piMax, pi[i]);
      water += a * moistPhysics.columnWater(pi, q, i);
      cloud += a * moistPhysics.columnWater(pi, qc, i);
      rain += a * precipitation[i];
    }
    for (let x = 0; x < u.length; x++) maxWind = Math.max(maxWind, Math.abs(u[x]));
    const interval = model.time - lastPrecipTime;
    const result = {
      mass: mass / area, meanSurfaceT: meanSurfaceT / area, piMin, piMax, maxWind,
      absorbedSolar: sums.absorbedSolar / area, outgoingLongwave: sums.outgoingLongwave / area, sensibleHeat: sums.sensibleHeat / area,
      evaporation: sums.evaporation / area, latentHeat: moistPhysics.latentHeat * sums.evaporation / area,
      columnWater: water / area, columnCloud: cloud / area, precipitation: interval > 0 ? rain / area / interval : 0,
    };
    precipitation.fill(0);
    lastPrecipTime = model.time;
    return result;
  };

  return model;
}
