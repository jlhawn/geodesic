import { buildMesh } from './mesh.module.js';
import { createSigmaCore, sigmaInterfaces, DIVERGENCE_DAMPING } from './dynamics/sigmaCore.module.js';
import { createRK4Arrays } from './dynamics/integrators.module.js';
import { createRadiation } from './physics/radiation.module.js';
import { createSurface } from './physics/surface.module.js';
import { createMoistPhysics } from './physics/moist.module.js';
import { createSeaIce, MELTING_POINT, FREEZING_POINT } from './physics/ice.module.js';
import { createOcean } from './ocean/layered.module.js';
import { createBoundaryLayer } from './physics/boundaryLayer.module.js';
import { createGeography, surfaceGeopotential } from './geography.module.js';
import { createLandSurface } from './physics/land.module.js';

export const SIDEREAL_DAY = 86164.0905;
export const STATE_NAMES = ['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'];
export const stateLengths = ({ K, C, E }) => ({ pi: C, theta: K * C, u: K * E, surfaceT: C, q: K * C, qc: K * C, ice: C });

/*
 * The climate model: the sigma-coordinate core with gray radiation, a
 * slab-ocean surface, bulk surface fluxes of heat and moisture,
 * boundary-layer drag, large-scale condensation, the convective plume
 * and dry convective adjustment, and a 5-day Rayleigh sponge above σ = 0.02 in
 * the role of gravity-wave drag on the polar-night jet. State is
 * [pi, theta, u, surfaceT, q, qc, ice] (vapour, cloud condensate, sea
 * ice thickness), each on a SharedArrayBuffer; with `moist: false` q
 * and qc are carried but never sourced, so they stay zero. surfaceT is
 * the skin temperature: the mixed layer's over open water, the ice's
 * over ice.
 *
 * One time-step tendency is four phases; each takes an index range so
 * the same code runs whole on one thread or sliced across workers:
 *   flux(layers)   mass flux and divergence
 *   column(cells)  dπ/dt, σ̇, Exner and geopotential
 *   vertex(verts)  kite-weighted π
 *   layer(layers)  θ and momentum tendencies, drag
 * After the RK4 step, once per step and applied to the state directly:
 *   ocean()        the dynamic ocean step (main thread only, before
 *                  physics)
 *   physics(cells) radiation, surface fluxes, evaporation, sea ice
 *   closure(layers) the ∇⁴ closures and the divergence damping
 *                  (`divergenceDamping`, DIVERGENCE_DAMPING)
 *   adjust(cells)  boundary-layer mixing, condensation, convection, filler
 *   mixMomentum(edges) boundary-layer mixing of the normal velocity
 *   dissipate(cells) the kinetic energy the closure and the mixing
 *                  removed, returned as heat
 * Arrays read across phases live in `shared`; `buffers` adopts another
 * instance's so a worker computes on the same memory.
 */
export function createModel(gridOrMesh, {
  radius, core: coreOptions = {}, radiation: radiationOptions = {}, surface: surfaceOptions = {}, moist: moistOptions = {}, ice: iceOptions = {}, ocean: oceanOptions = {}, boundaryLayer: boundaryLayerOptions = {},
  topography = null, geography: geographyOptions = {}, land: landOptions = {}, terrain = true,
  physics = true, moist = true, nu4Hours = 3, divergenceDamping = DIVERGENCE_DAMPING, buffers = null, levels = sigmaInterfaces(),
} = {}) {
  const mesh = gridOrMesh.nCells ? gridOrMesh : buildMesh(gridOrMesh, { radius, omega: 2 * Math.PI / SIDEREAL_DAY });
  const geography = topography ? createGeography(mesh, topography, geographyOptions) : null;
  const phis = geography && terrain ? surfaceGeopotential(mesh, geography) : null;
  const dragCoefficients = geography ? Float64Array.from(geography.land, (l) => (l ? landOptions.dragCoefficient ?? 1.5e-3 : surfaceOptions.dragCoefficient ?? 1.5e-3)) : null;
  let spacing = 0;
  for (let e = 0; e < mesh.nEdges; e++) spacing += mesh.dcEdge[e];
  spacing /= mesh.nEdges;
  const nu4 = Math.pow(spacing / Math.PI, 4) / (nu4Hours * 3600);
  const core = createSigmaCore(mesh, { levels, nu4, nu4Theta: nu4, divergenceDamping, splitClosure: true, buffers: buffers ? buffers.core : null, surfaceGeopotential: phis, ...coreOptions });
  const { K, C, E, V } = core.diagnostics;
  const radiation = createRadiation(mesh, core, { buffers: buffers ? buffers.radiation : null, exchangeCoefficients: dragCoefficients, ...radiationOptions });
  const boundaryLayer = physics && boundaryLayerOptions !== false ? createBoundaryLayer(mesh, core, { buffers: buffers ? buffers.boundaryLayer : null, dragCoefficients, land: geography ? geography.land : null, deckTop: radiation.mlmTop, deckGate: radiation.mlmGate, stratiform: radiation.stratiform, longwave: radiation.longwave, ...boundaryLayerOptions }) : null;
  if (boundaryLayer && boundaryLayer.turbulence === 'moist') radiation.useBoundaryLayer(boundaryLayer.regime, boundaryLayer.mixingTop, boundaryLayer.buoyancyFlux);
  const surface = createSurface(mesh, core, { topSigma: 0.02, topDragDays: 5, buffers: buffers ? buffers.surface : null, dragCoefficients, ...surfaceOptions });
  const gustiness = surfaceOptions.gustiness ?? 3;
  const moistPhysics = createMoistPhysics(mesh, core, {
    buffers: buffers ? buffers.moist : null, boundaryDepth: boundaryLayer ? boundaryLayer.depth : null, boundaryRegime: boundaryLayer ? boundaryLayer.regime : null, deckGate: radiation.mlmGate,
    boundaryTop: boundaryLayer && boundaryLayer.turbulence === 'moist' ? boundaryLayer.mixingTop : null, stratiform: radiation.stratiform,
    surfaceBuoyancy: boundaryLayer ? boundaryLayer.buoyancyFlux : null, frictionVelocity: boundaryLayer ? boundaryLayer.friction : null, ...moistOptions,
  });
  if (moist) radiation.useCumulus(moistPhysics.cumulusCover, moistPhysics.cumulusWater);
  const ocean = physics && oceanOptions !== false ? createOcean(mesh, { buffers: buffers ? buffers.ocean : null, geography, ...oceanOptions }) : null;
  const land = physics && geography ? createLandSurface(mesh, geography, { buffers: buffers ? buffers.land : null, ...landOptions }) : null;
  const landMask = geography ? geography.land : null;
  const sharedCapacity = !ocean && buffers && buffers.ocean ? new Float64Array(buffers.ocean.capacity) : null;
  const seaIce = createSeaIce(mesh, {
    buffers: { ...(buffers && buffers.ice ? buffers.ice : {}), ...(land ? { snow: land.shared.snow } : {}) },
    ...(ocean || sharedCapacity ? { heatCapacity: ocean ? ocean.capacity : sharedCapacity } : {}),
    ...iceOptions,
  });
  moistPhysics.useSeaIce(seaIce.concentration);
  const totals = { absorbedSolar: 0, atmosphereSolar: 0, outgoingLongwave: 0, sensibleHeat: 0, evaporation: 0, insolation: 0, reflectedSolar: 0 };
  const surfaceAlbedo = new Float64Array(C), diffuseAlbedo = new Float64Array(C), wetness = new Float64Array(C).fill(1), openSea = new Float64Array(C), stressScratch = new Float64Array(E);
  const fluxT = new Float64Array(C), directContrast = new Float64Array(C), diffuseContrast = new Float64Array(C);
  const runoffSeen = land ? new Float64Array(C) : null, runoffStep = land ? new Float64Array(C) : null;

  const lengths = stateLengths({ K, C, E });
  const stateArray = (name) => new Float64Array(buffers && buffers.state && buffers.state[name] ? buffers.state[name] : new SharedArrayBuffer(8 * lengths[name]));
  const state = STATE_NAMES.map(stateArray);
  const forcing = STATE_NAMES.map((name) => new Float64Array(lengths[name]));
  const fluxState = state.map((array, a) => (a === 3 ? fluxT : array)), dryFluxState = fluxState.slice(0, 4);

  const phases = {
    flux(input, kFrom, kTo) { core.phaseFlux(input, kFrom, kTo); },
    column(input, out, iFrom, iTo) { core.phaseColumn(input, out, iFrom, iTo); },
    vertex(input, vFrom, vTo) { core.phaseVertex(input, vFrom, vTo); },
    layer(input, out, kFrom, kTo, part = 'all') {
      core.phaseLayer(input, out, kFrom, kTo, part);
      if (physics && part !== 'tracers') surface.applyLayers(input, out, kFrom, kTo);
      if (physics && part !== 'momentum') surface.heatLayers(input, out, kFrom, kTo);
    },
    ocean(dt) {
      if (!physics || !ocean) return;
      if (land) for (let i = 0; i < C; i++) { const total = land.runoff[i]; runoffStep[i] = total >= runoffSeen[i] ? total - runoffSeen[i] : total; runoffSeen[i] = total; }
      ocean.accumulate(radiation.evaporation, moist ? moistPhysics.rain : null, dt, runoffStep);
      ocean.advance(state[3], state[6], seaIce.oceanFlux, () => surface.stress(state, stressScratch), dt, seaIce.concentration);
    },
    physics(iFrom, iTo, dt, sums) {
      if (!physics) return;
      surface.lowestWindSpeed(state[2], iFrom, iTo);
      const bottom = (K - 1) * C;
      forcing[3].fill(0, iFrom, iTo);
      for (let k = 0; k < K; k++) { forcing[1].fill(0, k * C + iFrom, k * C + iTo); forcing[4].fill(0, k * C + iFrom, k * C + iTo); }
      for (let i = iFrom; i < iTo; i++) {
        fluxT[i] = state[3][i];
        if (land && landMask[i]) { surfaceAlbedo[i] = diffuseAlbedo[i] = land.albedo(i); wetness[i] = land.wetness(i, dragCoefficients[i] * Math.max(surface.windSpeed[i], gustiness), state[3][i]); continue; }
        const h = state[6][i], area = seaIce.cover(i, h), mu = radiation.cosZenith(i);
        surfaceAlbedo[i] = seaIce.albedo(h, mu, seaIce.snow[i], area); diffuseAlbedo[i] = seaIce.albedo(h, null, seaIce.snow[i], area);
        openSea[i] = 1 - area;
        if (h > 0 && area < 1) {
          fluxT[i] = area * state[3][i] + (1 - area) * FREEZING_POINT;
          directContrast[i] = seaIce.albedoContrast(h, mu, seaIce.snow[i]); diffuseContrast[i] = seaIce.albedoContrast(h, null, seaIce.snow[i]);
        }
      }
      radiation.apply(moist ? fluxState : dryFluxState, forcing, surface.windSpeed, sums, iFrom, iTo, surfaceAlbedo, diffuseAlbedo, land ? wetness : null, moist ? openSea : null, boundaryLayer ? boundaryLayer.depth : null, dt);
      for (let k = 0; k < K; k++) for (let i = k * C + iFrom; i < k * C + iTo; i++) state[1][i] += dt * forcing[1][i];
      const { surfaceShortwave, surfaceDirect } = radiation;
      for (let i = iFrom; i < iTo; i++) {
        if (land && landMask[i]) { land.update(i, state[3], radiation.surfaceFlux, radiation.evaporation[i], dt); continue; }
        const h = state[6][i], partial = h > 0 && seaIce.cover(i, h) < 1;
        seaIce.update(state[3], state[6], radiation.surfaceFlux, i, dt, partial ? surfaceDirect[i] * directContrast[i] + (surfaceShortwave[i] - surfaceDirect[i]) * diffuseContrast[i] : 0);
      }
      if (moist) for (let i = bottom + iFrom; i < bottom + iTo; i++) state[4][i] += dt * forcing[4][i];
      if (boundaryLayer) boundaryLayer.diagnose(state, iFrom, iTo);
    },
    closure(kFrom, kTo, dt, part = 'all') { core.phaseClosure(state, kFrom, kTo, dt, part); },
    adjust(iFrom, iTo, dt) {
      if (!physics) return;
      if (boundaryLayer) for (let i = iFrom; i < iTo; i++) boundaryLayer.mixColumn(i, state[0], state[1], moist ? state[4] : null, moist ? state[5] : null, dt);
      if (moist) moistPhysics.adjust(state, iFrom, iTo, dt);
      if (moist) {
        const bottom = (K - 1) * C, { exnerLayer: exner, cp, g, dSigma } = core.diagnostics;
        for (let i = iFrom; i < iTo; i++) {
          const airTemperature = state[1][bottom + i] * exner[bottom + i], amount = moistPhysics.rain[i];
          let frozen;
          if (land && landMask[i]) { land.deposit(i, amount, airTemperature); frozen = airTemperature < MELTING_POINT; }
          else frozen = seaIce.deposit(i, amount, airTemperature, state[6], state[3]);
          if (frozen) state[1][bottom + i] += seaIce.latentHeatFusion * amount * g / (cp * state[0][i] * dSigma[K - 1] * exner[bottom + i]);
        }
      }
      surface.convectiveAdjustment(state[0], state[1], iFrom, iTo, moist ? state[4] : null, moist ? state[5] : null);
    },
    mixMomentum(eFrom, eTo, dt) {
      if (physics && boundaryLayer) boundaryLayer.mixEdges(state[0], state[2], eFrom, eTo, dt, core.arrays.dissipation);
      if (physics && moist) moistPhysics.transportMomentum(state[0], state[2], eFrom, eTo, dt, core.arrays.dissipation);
    },
    dissipate(iFrom, iTo) {
      if (!physics) return;
      const theta = state[1], lost = core.arrays.dissipation, { exnerLayer, cp } = core.diagnostics;
      const { nEdgesOnCell, edgesOnCell, maxEdges, dcEdge, dvEdge, areaCell } = mesh;
      for (let k = 0; k < K; k++) {
        for (let i = iFrom; i < iTo; i++) {
          let sum = 0;
          for (let m = 0; m < nEdgesOnCell[i]; m++) { const e = edgesOnCell[maxEdges * i + m]; sum += dcEdge[e] * dvEdge[e] * lost[k * E + e]; }
          theta[k * C + i] += 0.25 * sum / areaCell[i] / (cp * exnerLayer[k * C + i]);
        }
      }
    },
  };

  function tendency(input, out) {
    phases.flux(input, 0, K);
    phases.column(input, out, 0, C);
    phases.vertex(input, 0, V);
    phases.layer(input, out, 0, K);
    out[3].fill(0);
    out[6].fill(0);
  }

  let rk4 = null;
  const model = {
    mesh, core, radiation, surface, moist: moistPhysics, seaIce, ocean, boundaryLayer, geography, land, surfaceGeopotential: phis, surfaceAlbedo, state, totals, phases, tendency, physics, moistOn: physics && moist, time: 0,
    shared: { core: core.shared, surface: surface.shared, moist: moistPhysics.shared, ice: seaIce.shared, radiation: radiation.shared, ocean: ocean ? ocean.shared : (buffers && buffers.ocean ? buffers.ocean : null), boundaryLayer: boundaryLayer ? boundaryLayer.shared : null, land: land ? land.shared : null, state: Object.fromEntries(STATE_NAMES.map((name, a) => [name, state[a].buffer])) },
  };

  model.oceanFields = (depth = 0) => (ocean ? ocean.fields(depth) : null);
  model.cloudWater = (i) => moistPhysics.columnWater(state[0], state[5], i) + radiation.stratusFraction[i] * radiation.stratus[i];

  model.step = function step(dt) {
    rk4 ??= createRK4Arrays(STATE_NAMES.map((name) => lengths[name]));
    radiation.setTime(model.time);
    rk4(tendency, state, dt);
    phases.ocean(dt);
    phases.physics(0, C, dt, totals);
    phases.closure(0, K, dt);
    phases.adjust(0, C, dt);
    phases.mixMomentum(0, E, dt);
    phases.dissipate(0, C);
    model.time += dt;
  };

  /*
   * Global means, the precipitation over the interval since the last
   * restartPrecipitation, which each call ends; over the same interval
   * the per-cell convective and large-scale rain go to
   * moist.convectiveRain and moist.largeScaleRain in mm/d.
   */
  let lastPrecipTime = 0;
  model.diagnostics = function diagnostics(sums = totals) {
    const [pi, theta, u, surfaceT, q, qc, ice] = state;
    const precipitation = moistPhysics.precipitation;
    let area = 0, mass = 0, meanSurfaceT = 0, piMin = Infinity, piMax = -Infinity, maxWind = 0, water = 0, cloud = 0, rain = 0, iceArea = 0, iceVolume = 0, albedoSum = 0;
    let landArea = 0, landT = 0, snowArea = 0, soilSum = 0;
    for (let i = 0; i < C; i++) {
      const a = mesh.areaCell[i];
      area += a;
      if (land && landMask[i]) { landArea += a; landT += a * surfaceT[i]; soilSum += a * land.soil[i]; if (land.snow[i] > 1) snowArea += a; }
      mass += a * pi[i];
      meanSurfaceT += a * surfaceT[i];
      piMin = Math.min(piMin, pi[i]);
      piMax = Math.max(piMax, pi[i]);
      water += a * moistPhysics.columnWater(pi, q, i);
      cloud += a * moistPhysics.columnWater(pi, qc, i);
      rain += a * precipitation[i];
      const cover = seaIce.cover(i, ice[i]);
      if (ice[i] > 0) { iceArea += a * cover; iceVolume += a * cover * ice[i]; }
      albedoSum += a * (land && landMask[i] ? land.albedo(i) : seaIce.albedo(ice[i], null, seaIce.snow[i], cover));
    }
    for (let x = 0; x < u.length; x++) maxWind = Math.max(maxWind, Math.abs(u[x]));
    const interval = model.time - lastPrecipTime;
    if (interval > 0) moistPhysics.readRain(interval);
    const result = {
      mass: mass / area, meanSurfaceT: meanSurfaceT / area, piMin, piMax, maxWind,
      absorbedSolar: sums.absorbedSolar / area, atmosphereSolar: sums.atmosphereSolar / area, outgoingLongwave: sums.outgoingLongwave / area, sensibleHeat: sums.sensibleHeat / area,
      evaporation: sums.evaporation / area, latentHeat: moistPhysics.latentHeat * sums.evaporation / area,
      columnWater: water / area, columnCloud: cloud / area, precipitation: interval > 0 ? rain / area / interval : 0,
      iceFraction: iceArea / area, iceThickness: iceArea > 0 ? iceVolume / iceArea : 0, surfaceAlbedo: albedoSum / area,
      planetaryAlbedo: sums.insolation > 0 ? sums.reflectedSolar / sums.insolation : 0,
      ...(ocean ? ocean.diagnostics() : {}),
      ...(land ? { landFraction: landArea / area, landMeanT: landArea > 0 ? landT / landArea : 0, snowFraction: landArea > 0 ? snowArea / landArea : 0, soilWater: landArea > 0 ? soilSum / landArea : 0, runoff: land.budget.runoff / area } : {}),
    };
    model.restartPrecipitation();
    return result;
  };
  model.restartPrecipitation = function restartPrecipitation() {
    moistPhysics.precipitation.fill(0);
    moistPhysics.convectivePrecipitation.fill(0);
    moistPhysics.largeScalePrecipitation.fill(0);
    lastPrecipTime = model.time;
  };

  return model;
}
