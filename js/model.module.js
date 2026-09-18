import { buildMesh } from './mesh.module.js';
import { createSigmaCore } from './dynamics/sigmaCore.module.js';
import { createRK4Arrays } from './dynamics/integrators.module.js';
import { createRadiation } from './physics/radiation.module.js';
import { createSurface } from './physics/surface.module.js';

export const SIDEREAL_DAY = 86164.0905;

/*
 * The climate model: the sigma-coordinate core with gray radiation, a
 * slab-ocean surface, bulk surface fluxes, boundary-layer drag, and dry
 * convective adjustment. State is [pi, theta, u, surfaceT].
 */
export function createModel(grid, {
  radius, core: coreOptions = {}, radiation: radiationOptions = {}, surface: surfaceOptions = {},
  physics = true, nu4Hours = 3,
} = {}) {
  const mesh = buildMesh(grid, { radius, omega: 2 * Math.PI / SIDEREAL_DAY });
  let spacing = 0;
  for (let e = 0; e < mesh.nEdges; e++) spacing += mesh.dcEdge[e];
  spacing /= mesh.nEdges;
  const nu4 = Math.pow(spacing / Math.PI, 4) / (nu4Hours * 3600);
  const core = createSigmaCore(mesh, { nu4, nu4Theta: nu4, ...coreOptions });
  const { K, C, E } = core.diagnostics;
  const radiation = createRadiation(mesh, core, radiationOptions);
  const surface = createSurface(mesh, core, surfaceOptions);
  const totals = { absorbedSolar: 0, outgoingLongwave: 0, sensibleHeat: 0 };

  core.setForcing((state, out) => {
    const [, , u] = state;
    const dSurfaceT = out[3];
    dSurfaceT.fill(0);
    if (!physics) return;
    const windSpeed = surface.lowestWindSpeed(u);
    radiation.apply(state, out, windSpeed, totals);
    surface.apply(state, out);
  });

  const state = [new Float64Array(C), new Float64Array(K * C), new Float64Array(K * E), new Float64Array(C)];
  const rk4 = createRK4Arrays([C, K * C, K * E, C]);
  const model = { mesh, core, radiation, surface, state, totals, time: 0, physics };

  model.step = function step(dt) {
    radiation.setTime(model.time);
    rk4(core.tendency, state, dt);
    if (physics) surface.convectiveAdjustment(state[0], state[1]);
    model.time += dt;
  };

  model.diagnostics = function diagnostics() {
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
    return { mass: mass / area, meanSurfaceT: meanSurfaceT / area, piMin, piMax, maxWind, absorbedSolar: totals.absorbedSolar / area, outgoingLongwave: totals.outgoingLongwave / area };
  };

  return model;
}
