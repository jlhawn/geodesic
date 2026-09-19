import { buildMesh } from './mesh.module.js';
import { createSigmaCore } from './dynamics/sigmaCore.module.js';
import { createRK4Arrays } from './dynamics/integrators.module.js';
import { createRadiation } from './physics/radiation.module.js';
import { createSurface } from './physics/surface.module.js';

export const SIDEREAL_DAY = 86164.0905;

/*
 * The climate model: the sigma-coordinate core with gray radiation, a
 * slab-ocean surface, bulk surface fluxes, boundary-layer drag, dry
 * convective adjustment, and a 10-day Rayleigh drag in the cap layer
 * above CAM's lid (σ < 0.005), where nothing else bounds the winter
 * jet. State is [pi, theta, u, surfaceT], each on a SharedArrayBuffer.
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
  radius, core: coreOptions = {}, radiation: radiationOptions = {}, surface: surfaceOptions = {},
  physics = true, nu4Hours = 3, buffers = null,
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
  const totals = { absorbedSolar: 0, outgoingLongwave: 0, sensibleHeat: 0 };

  const stateArray = (name, n) => new Float64Array(buffers && buffers.state && buffers.state[name] ? buffers.state[name] : new SharedArrayBuffer(8 * n));
  const state = [stateArray('pi', C), stateArray('theta', K * C), stateArray('u', K * E), stateArray('surfaceT', C)];

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
      if (physics) radiation.apply(input, out, surface.windSpeed, sums, iFrom, iTo);
    },
    adjust(iFrom, iTo) { if (physics) surface.convectiveAdjustment(state[0], state[1], iFrom, iTo); },
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
    mesh, core, radiation, surface, state, totals, phases, tendency, physics, time: 0,
    shared: { core: core.shared, surface: surface.shared, state: { pi: state[0].buffer, theta: state[1].buffer, u: state[2].buffer, surfaceT: state[3].buffer } },
  };

  model.step = function step(dt) {
    rk4 ??= createRK4Arrays([C, K * C, K * E, C]);
    radiation.setTime(model.time);
    rk4(tendency, state, dt);
    phases.adjust(0, C);
    model.time += dt;
  };

  model.diagnostics = function diagnostics(sums = totals) {
    const [pi, theta, u, surfaceT] = state;
    let area = 0, mass = 0, meanSurfaceT = 0, piMin = Infinity, piMax = -Infinity, maxWind = 0;
    for (let i = 0; i < C; i++) {
      area += mesh.areaCell[i];
      mass += mesh.areaCell[i] * pi[i];
      meanSurfaceT += mesh.areaCell[i] * surfaceT[i];
      piMin = Math.min(piMin, pi[i]);
      piMax = Math.max(piMax, pi[i]);
    }
    for (let x = 0; x < u.length; x++) maxWind = Math.max(maxWind, Math.abs(u[x]));
    return { mass: mass / area, meanSurfaceT: meanSurfaceT / area, piMin, piMax, maxWind, absorbedSolar: sums.absorbedSolar / area, outgoingLongwave: sums.outgoingLongwave / area };
  };

  return model;
}
