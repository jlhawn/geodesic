import { buildMesh } from '../mesh.module.js';
import { createSigmaCore, sigmaInterfaces } from '../dynamics/sigmaCore.module.js';
import { createSeaIce } from '../physics/ice.module.js';
import { createSurface } from '../physics/surface.module.js';
import { createRadiation } from '../physics/radiation.module.js';
import { createMoistPhysics } from '../physics/moist.module.js';
import { LATENT_HEAT } from '../physics/moist.module.js';
import { SIDEREAL_DAY } from '../model.module.js';
import { createGpuCore } from './core.gpu.js';
import { createLayeredOcean } from './layeredOcean.gpu.js';
import { createGeography, surfaceGeopotential } from '../geography.module.js';
import { createLandSurface } from '../physics/land.module.js';

/*
 * The whole model on the GPU behind the CPU model's interface: `state`
 * holds double-precision mirrors that `sync` refreshes from the device,
 * `step` is asynchronous, `diagnostics` reads the per-cell energy terms
 * back and sums them on the CPU, and `ocean` carries the same
 * initialize/load/serialize contract. A CPU sigma core on the same mesh
 * serves the page's level fields from the mirrored state.
 */
export async function createGpuModel(gridOrMesh, {
  radius, nu4Hours = 3, radiation = {}, ice = {}, moist = {}, boundaryLayer = {}, ocean: oceanOptions = {}, surface = {},
  topography = null, geography: geographyOptions = {}, land: landOptions = {}, terrain = true,
} = {}) {
  const mesh = gridOrMesh.nCells ? gridOrMesh : buildMesh(gridOrMesh, { radius, omega: 2 * Math.PI / SIDEREAL_DAY });
  const geography = topography ? createGeography(mesh, topography, geographyOptions) : null;
  const phis = geography && terrain ? surfaceGeopotential(mesh, geography) : null;
  const dragCoefficients = geography ? Float64Array.from(geography.land, (l) => (l ? landOptions.dragCoefficient ?? 3e-3 : surface.dragCoefficient ?? 1.5e-3)) : null;
  let spacing = 0;
  for (let e = 0; e < mesh.nEdges; e++) spacing += mesh.dcEdge[e];
  spacing /= mesh.nEdges;
  const nu4 = Math.pow(spacing / Math.PI, 4) / (nu4Hours * 3600);
  const core = createSigmaCore(mesh, { surfaceGeopotential: phis });
  const { K, C, E } = core.diagnostics;
  const physics = {
    ...radiation, ...ice, ...moist, ...boundaryLayer,
    landed: !!geography, landHeatCapacity: landOptions.heatCapacity ?? 1e6, bucketCapacity: landOptions.bucketCapacity ?? 150, wetnessThreshold: landOptions.wetnessThreshold ?? 0.75,
    landAlbedo: landOptions.albedo ?? 0.2, snowAlbedo: landOptions.snowAlbedo ?? 0.55, fullSnow: landOptions.fullSnow ?? 20,
  };
  const gpu = await createGpuCore(mesh, { nu4, nu4Theta: nu4, physics, topSigma: surface.topSigma ?? 0.02, topDragDays: surface.topDragDays ?? 5, surfaceGeopotential: phis });
  const seaIce = createSeaIce(mesh, ice);
  const radiationCpu = createRadiation(mesh, core, radiation);
  const surfaceCpu = createSurface(mesh, core, { topSigma: 0.02, topDragDays: 5, ...surface });
  const moistCpu = createMoistPhysics(mesh, core, moist);
  const gpuOcean = oceanOptions === false ? null : createLayeredOcean(gpu, { ...oceanOptions, geography });
  const landCpu = geography ? createLandSurface(mesh, geography, landOptions) : null;
  let oceanCounter = 0;
  if (gpuOcean) gpu.hooks.beforePhysics = async (dt) => {
    gpuOcean.accumulateFreshwater(dt);
    if (++oceanCounter % gpuOcean.everySteps !== 0) return;
    await gpuOcean.advanceCoupled(gpuOcean.everySteps * dt);
  };
  const lengths = [C, K * C, K * E, C, K * C, K * C, C];
  const state = lengths.map((n) => new Float64Array(n));
  const precipitation = new Float64Array(C);
  let lastOcean = null;
  let dirty = true, lastPrecipTime = 0;

  const model = { mesh, core, seaIce, radiation: radiationCpu, surface: surfaceCpu, geography, surfaceGeopotential: phis, state, time: 0, physics: true, moistOn: true, gpu, engine: 'gpu' };
  model.moist = { precipitation, columnWater: moistCpu.columnWater, latentHeat: LATENT_HEAT, budget: moistCpu.budget };
  model.oceanFields = () => lastOcean;
  model.oceanEngine = gpuOcean;

  function pushState() {
    gpu.upload(state);
    gpu.uploadPhysics({ land: geography ? geography.land : null, drag: dragCoefficients, soil: landCpu ? landCpu.soil : null, snow: landCpu ? landCpu.snow : null });
    if (gpuOcean) gpuOcean.initialize(state[3], state[6]);
    dirty = false;
  }
  model.load = function load() { pushState(); };

  async function sync() {
    if (!dirty) return;
    const arrays = await gpu.download();
    arrays.forEach((a, i) => state[i].set(a));
    core.diagnose(state[0], state[1], state[4], state[5]);
    dirty = false;
  }
  model.sync = sync;

  model.step = async function step(dt) {
    await gpu.stepModel(dt, model.time);
    model.time += dt;
    dirty = true;
  };

  function refreshLand(ph) {
    if (!landCpu) return;
    landCpu.soil.set(ph.SOIL.subarray(0, C));
    landCpu.snow.set(ph.SNOW.subarray(0, C));
    for (let i = 0; i < C; i++) { landCpu.runoff[i] += ph.RUNOFF[i]; landCpu.budget.runoff += mesh.areaCell[i] * ph.RUNOFF[i]; }
    gpu.device.queue.writeBuffer(gpu.buffers.PH, 4 * gpu.layout.PH.RUNOFF, new Float32Array(C));
    if (gpuOcean) gpuOcean.forgetAccumulated('runoff');
  }

  model.diagnostics = async function diagnostics() {
    await sync();
    if (gpuOcean) gpuOcean.accumulateFreshwater(0);
    const ph = await gpu.downloadPhysics();
    refreshLand(ph);
    const [pi, , u, surfaceT, q, qc, iceField] = state;
    let area = 0, mass = 0, meanSurfaceT = 0, piMin = Infinity, piMax = -Infinity, maxWind = 0, water = 0, cloud = 0, rain = 0, iceArea = 0, iceVolume = 0, albedoSum = 0;
    let landArea = 0, landT = 0, snowArea = 0, soilSum = 0;
    const sums = { absorbedSolar: 0, outgoingLongwave: 0, sensibleHeat: 0, evaporation: 0, insolation: 0, reflectedSolar: 0 };
    for (let i = 0; i < C; i++) {
      const a = mesh.areaCell[i];
      area += a; mass += a * pi[i]; meanSurfaceT += a * surfaceT[i];
      if (landCpu && geography.land[i]) { landArea += a; landT += a * surfaceT[i]; soilSum += a * landCpu.soil[i]; if (landCpu.snow[i] > 1) snowArea += a; }
      piMin = Math.min(piMin, pi[i]); piMax = Math.max(piMax, pi[i]);
      water += a * model.moist.columnWater(pi, q, i); cloud += a * model.moist.columnWater(pi, qc, i);
      precipitation[i] = ph.RAIN[i];
      rain += a * ph.RAIN[i];
      if (iceField[i] > 0) { iceArea += a; iceVolume += a * iceField[i]; }
      albedoSum += a * ph.ADIF[i];
      sums.absorbedSolar += a * ph.ABS[i]; sums.outgoingLongwave += a * ph.OLR[i]; sums.sensibleHeat += a * ph.SH[i];
      sums.evaporation += a * ph.EVAP[i]; sums.insolation += a * ph.INS[i]; sums.reflectedSolar += a * ph.REFL[i];
    }
    radiationCpu.outgoing.set(ph.OLR.subarray(0, C)); radiationCpu.surfaceShortwave.set(ph.SWDN.subarray(0, C));
    for (let x = 0; x < u.length; x++) maxWind = Math.max(maxWind, Math.abs(u[x]));
    const interval = model.time - lastPrecipTime;
    gpu.device.queue.writeBuffer(gpu.buffers.PH, 4 * gpu.layout.PH.RAIN, new Float32Array(C));
    if (gpuOcean) gpuOcean.forgetAccumulated('rain');
    lastPrecipTime = model.time;
    const result = {
      mass: mass / area, meanSurfaceT: meanSurfaceT / area, piMin, piMax, maxWind,
      absorbedSolar: sums.absorbedSolar / area, outgoingLongwave: sums.outgoingLongwave / area, sensibleHeat: sums.sensibleHeat / area,
      evaporation: sums.evaporation / area, latentHeat: LATENT_HEAT * sums.evaporation / area,
      columnWater: water / area, columnCloud: cloud / area, precipitation: interval > 0 ? rain / area / interval : 0,
      iceFraction: iceArea / area, iceThickness: iceArea > 0 ? iceVolume / iceArea : 0, surfaceAlbedo: albedoSum / area,
      planetaryAlbedo: sums.insolation > 0 ? sums.reflectedSolar / sums.insolation : 0,
      ...(landCpu ? { landFraction: landArea / area, landMeanT: landArea > 0 ? landT / landArea : 0, snowFraction: landArea > 0 ? snowArea / landArea : 0, soilWater: landArea > 0 ? soilSum / landArea : 0, runoff: landCpu.budget.runoff / area } : {}),
    };
    if (gpuOcean) {
      const o = await gpuOcean.download();
      lastOcean = o;
      Object.assign(result, gpuOcean.diagnosticsFrom(o));
    }
    return result;
  };

  model.ocean = gpuOcean ? {
    initialize(surfaceT, iceField) { gpuOcean.initialize(surfaceT, iceField); },
    load(saved, surfaceT, iceField) { gpuOcean.upload(saved, surfaceT, iceField); },
    async serialize() { return gpuOcean.serialize(); },
  } : null;

  model.land = landCpu ? {
    soil: landCpu.soil, snow: landCpu.snow, runoff: landCpu.runoff, land: geography.land, budget: landCpu.budget, albedo: landCpu.albedo, wetness: landCpu.wetness, water: landCpu.water, bucketCapacity: landCpu.bucketCapacity,
    initialize() { landCpu.initialize(); gpu.uploadLand({ soil: landCpu.soil, snow: landCpu.snow }); },
    load(saved) { landCpu.load(saved); gpu.uploadLand({ soil: landCpu.soil, snow: landCpu.snow }); },
    async serialize() { if (gpuOcean) gpuOcean.accumulateFreshwater(0); refreshLand(await gpu.downloadPhysics()); return landCpu.serialize(); },
  } : null;

  return model;
}
