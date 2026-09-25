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
import { readRanges } from './device.module.js';
import { RAIN_MEMORY } from '../frames.module.js';
import { createLandSurface } from '../physics/land.module.js';

/*
 * The whole model on the GPU behind the CPU model's interface: `state`
 * holds double-precision mirrors that only `sync` refreshes from the
 * device (the land's soil and snow only `land.serialize`, its runoff
 * when the diagnostics are taken), `step` only
 * queues work, `beginFrame` computes the page's fields and the
 * diagnostics on the device, and `ocean` carries the same
 * initialize/load/serialize contract.
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
  let dirty = true, lastFrameTime = 0;

  const model = { mesh, core, seaIce, radiation: radiationCpu, surface: surfaceCpu, geography, surfaceGeopotential: phis, state, time: 0, physics: true, moistOn: true, gpu, engine: 'gpu' };
  model.moist = { columnWater: moistCpu.columnWater, latentHeat: LATENT_HEAT, budget: moistCpu.budget };
  model.oceanEngine = gpuOcean;

  function pushState() {
    gpu.upload(state);
    gpu.uploadPhysics({ land: geography ? geography.land : null, drag: dragCoefficients, soil: landCpu ? landCpu.soil : null, snow: landCpu ? landCpu.snow : null });
    gpu.clearFrame();
    if (gpuOcean) gpuOcean.initialize(state[3], state[6]);
    lastFrameTime = model.time;
    dirty = false;
  }
  model.load = function load() { pushState(); };

  async function sync() {
    if (!dirty) return;
    const arrays = await gpu.download();
    arrays.forEach((a, i) => state[i].set(a));
    dirty = false;
  }
  model.sync = sync;

  model.step = async function step(dt) {
    await gpu.stepModel(dt, model.time);
    model.time += dt;
    dirty = true;
  };
  model.settle = () => gpu.device.queue.onSubmittedWorkDone();
  model.destroy = () => { for (const buffer of [...Object.values(gpu.buffers), ...(gpuOcean ? Object.values(gpuOcean.buffers) : [])]) buffer.destroy(); };

  /*
   * One frame of output: the named fields (see frames.module.js) at the
   * pressure level, and the diagnostics when asked, all computed on the
   * device and read back alone. The work is queued before this returns,
   * so a caller can queue steps behind it; the promise resolves once the
   * copies land. Every frame also advances the three-hour rain and the
   * runoff tally.
   */
  model.beginFrame = function beginFrame({ level = 'surface', fields = [], diagnostics: summarize = false } = {}) {
    if (gpuOcean) gpuOcean.accumulateFreshwater(0);
    const time = model.time, interval = time - lastFrameTime;
    lastFrameTime = time;
    const atmosphere = gpu.frame({ pressure: level === 'surface' ? 0 : 100 * level, keep: Math.exp(-Math.max(0, interval) / RAIN_MEMORY), fields, diagnostics: summarize });
    if (gpuOcean) { gpuOcean.forgetAccumulated('rain'); gpuOcean.forgetAccumulated('runoff'); }
    const ocean = gpuOcean ? gpuOcean.frame({ fields, diagnostics: summarize }) : null;
    return Promise.all([atmosphere, ocean]).then(([a, o]) => {
      const out = { time, level, fields: { ...a.fields, ...(o ? o.fields : {}) }, diagnostics: null };
      if (!a.sums) return out;
      const s = a.sums, area = s.area;
      if (landCpu) {
        let total = 0;
        for (let i = 0; i < C; i++) { landCpu.runoff[i] += a.runoff[i]; total += mesh.areaCell[i] * a.runoff[i]; }
        landCpu.budget.runoff += total;
      }
      out.diagnostics = {
        mass: s.mass / area, meanSurfaceT: s.surfaceT / area, piMin: s.piMin, piMax: s.piMax, maxWind: s.maxWind,
        absorbedSolar: s.absorbedSolar / area, outgoingLongwave: s.outgoingLongwave / area, sensibleHeat: s.sensibleHeat / area,
        evaporation: s.evaporation / area, latentHeat: LATENT_HEAT * s.evaporation / area,
        columnWater: s.water / area, columnCloud: s.cloud / area, precipitation: interval > 0 ? s.rain / area / interval : s.recentRain / area / RAIN_MEMORY,
        iceFraction: s.iceArea / area, iceThickness: s.iceArea > 0 ? s.iceVolume / s.iceArea : 0, surfaceAlbedo: s.albedo / area,
        planetaryAlbedo: s.insolation > 0 ? s.reflectedSolar / s.insolation : 0,
        ...(landCpu ? { landFraction: s.landArea / area, landMeanT: s.landArea > 0 ? s.landT / s.landArea : 0, snowFraction: s.landArea > 0 ? s.snowArea / s.landArea : 0, soilWater: s.landArea > 0 ? s.soil / s.landArea : 0, runoff: landCpu.budget.runoff / area } : {}),
        ...(o && o.diagnostics ? o.diagnostics : {}),
      };
      return out;
    });
  };
  model.diagnostics = async function diagnostics() { return (await model.beginFrame({ diagnostics: true })).diagnostics; };

  model.ocean = gpuOcean ? {
    initialize(surfaceT, iceField) { gpuOcean.initialize(surfaceT, iceField); },
    load(saved, surfaceT, iceField) { gpuOcean.upload(saved, surfaceT, iceField); },
    async serialize() { return gpuOcean.serialize(); },
  } : null;

  model.land = landCpu ? {
    soil: landCpu.soil, snow: landCpu.snow, runoff: landCpu.runoff, land: geography.land, budget: landCpu.budget, albedo: landCpu.albedo, wetness: landCpu.wetness, water: landCpu.water, bucketCapacity: landCpu.bucketCapacity,
    initialize() { landCpu.initialize(); gpu.uploadLand({ soil: landCpu.soil, snow: landCpu.snow }); },
    load(saved) { landCpu.load(saved); gpu.uploadLand({ soil: landCpu.soil, snow: landCpu.snow }); },
    async serialize() {
      const [soil, snow] = await readRanges(gpu.device, gpu.buffers.PH, [{ offset: gpu.layout.PH.SOIL, length: C }, { offset: gpu.layout.PH.SNOW, length: C }]);
      landCpu.soil.set(soil); landCpu.snow.set(snow);
      return landCpu.serialize();
    },
  } : null;

  return model;
}
