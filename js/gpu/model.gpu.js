import { buildMesh } from '../mesh.module.js';
import { createSigmaCore, sigmaInterfaces, sigmaGridName, DIVERGENCE_DAMPING } from '../dynamics/sigmaCore.module.js';
import { createSeaIce } from '../physics/ice.module.js';
import { createSurface, SEA_DRAG, LAND_DRAG, TOP_DRAG } from '../physics/surface.module.js';
import { spongeRates, spongeSigmaFor, SPONGE } from '../dynamics/sponge.module.js';
import { createRadiation } from '../physics/radiation.module.js';
import { createMoistPhysics } from '../physics/moist.module.js';
import { LATENT_HEAT } from '../physics/moist.module.js';
import { SIDEREAL_DAY } from '../model.module.js';
import { createGpuCore } from './core.gpu.js';
import { createLayeredOcean } from './layeredOcean.gpu.js';
import { createGeography, surfaceGeopotential } from '../geography.module.js';
import { readRanges } from './device.module.js';
import { RAIN_MEMORY, VERTICAL_MEMORY } from '../frames.module.js';
import { createLandSurface } from '../physics/land.module.js';
import { loadClimatology } from '../ocean/climatology.module.js';

const VEGETATION_OPTIONS = ['vegetation', 'bareAlbedo', 'vegetatedAlbedo', 'soilDarkening', 'wetSoilAlbedo', 'darkeningWetness', 'rootZoneCapacity', 'dryWetness', 'wetWetness', 'growthTime', 'declineTime', 'snowDeclineTime'];

/*
 * The whole model on the GPU behind the CPU model's interface: `state`,
 * `seaIce.concentration` and the deck's carried state — the running-mean
 * subsidence, inversion height and gate `radiation.mlmSubsidence`,
 * `mlmHeight` and `mlmGate` — the boundary layer's depth
 * `boundaryLayer.depth` (the height of its top, which the next step's
 * physics reads before it diagnoses its own) with the moist scheme's
 * mixing top, regime and surface buoyancy flux `boundaryLayer.mixingTop`,
 * `regime` and `buoyancyFlux` (the radiation's variance cover reads the
 * flux of the step before), the per-cell convective and
 * large-scale rain `moist.convectiveRain` and `moist.largeScaleRain`,
 * the means in mm/d that each diagnostics frame takes over the interval
 * since the frame before, and the per-cell absorbed sunlight, outgoing
 * longwave and planetary albedo `radiation.meanAbsorbedSolar`,
 * `meanOutgoingLongwave` and `meanPlanetaryAlbedo` (with clearSkyPass
 * also the cloud effects `meanShortwaveCloudEffect` and
 * `meanLongwaveCloudEffect`), the means each
 * diagnostics frame takes over the steps since the frame before, hold
 * double-precision mirrors that only
 * `sync` refreshes from the device and `load` sends to it (the land's
 * soil, snow and vegetation only `land.serialize`, its runoff when the
 * diagnostics are taken), `step` only queues work, `beginFrame`
 * computes the page's fields and the diagnostics on the device, and
 * `ocean` carries the same initialize/load/serialize contract (starting
 * from a climatology, initialize also sends the device the sea-surface
 * temperature it sets).
 */
export async function createGpuModel(gridOrMesh, {
  radius, nu4Hours = 3, divergenceDamping = DIVERGENCE_DAMPING, radiation = {}, ice = {}, moist = {}, boundaryLayer = {}, ocean: oceanOptions = {}, surface = {},
  topography = null, geography: geographyOptions = {}, land: landOptions = {}, terrain = true, levels = sigmaInterfaces(), gravityWaves = {},
} = {}) {
  const mesh = gridOrMesh.nCells ? gridOrMesh : buildMesh(gridOrMesh, { radius, omega: 2 * Math.PI / SIDEREAL_DAY });
  const geography = topography ? createGeography(mesh, topography, geographyOptions) : null;
  const phis = geography && terrain ? surfaceGeopotential(mesh, geography) : null;
  const dragCoefficients = geography ? Float64Array.from(geography.land, (l) => (l ? landOptions.dragCoefficient ?? LAND_DRAG : surface.dragCoefficient ?? SEA_DRAG)) : null;
  let spacing = 0;
  for (let e = 0; e < mesh.nEdges; e++) spacing += mesh.dcEdge[e];
  spacing /= mesh.nEdges;
  const nu4 = Math.pow(spacing / Math.PI, 4) / (nu4Hours * 3600);
  const core = createSigmaCore(mesh, { levels, surfaceGeopotential: phis });
  const { K, C, E } = core.diagnostics;
  const physics = {
    ...radiation, ...ice, ...moist, ...boundaryLayer,
    landed: !!geography, landHeatCapacity: landOptions.heatCapacity ?? 1e6, bucketCapacity: landOptions.bucketCapacity ?? 150, wetnessThreshold: landOptions.wetnessThreshold ?? 0.75,
    landAlbedo: landOptions.albedo ?? 0.2, snowAlbedo: landOptions.snowAlbedo ?? 0.55, fullSnow: landOptions.fullSnow ?? 20,
    ...Object.fromEntries(VEGETATION_OPTIONS.filter((key) => landOptions[key] !== undefined).map((key) => [key, landOptions[key]])),
  };
  const sponge = spongeRates(core.sigmaMid, surface.spongeSigma ?? spongeSigmaFor(sigmaGridName(levels)), surface.spongeDays ?? SPONGE.days);
  const gpu = await createGpuCore(mesh, { levels, nu4, nu4Theta: nu4, divergenceDamping, physics, topSigma: surface.topSigma ?? TOP_DRAG.sigma, topDragDays: surface.topDragDays ?? TOP_DRAG.days, spongeRates: sponge, gravityWaves, surfaceGeopotential: phis });
  const seaIce = createSeaIce(mesh, ice);
  const radiationCpu = createRadiation(mesh, core, radiation);
  const surfaceCpu = createSurface(mesh, core, { topSigma: TOP_DRAG.sigma, topDragDays: TOP_DRAG.days, ...surface });
  const moistCpu = createMoistPhysics(mesh, core, moist);
  const gpuOcean = oceanOptions === false ? null : createLayeredOcean(gpu, { ...oceanOptions, climatology: await loadClimatology(oceanOptions.climatology ?? null), geography });
  const landCpu = geography ? createLandSurface(mesh, geography, landOptions) : null;
  let oceanCounter = 0;
  if (gpuOcean) gpu.hooks.beforePhysics = async (dt) => {
    gpuOcean.accumulateFreshwater(dt);
    if (++oceanCounter % gpuOcean.everySteps !== 0) return;
    await gpuOcean.advanceCoupled(gpuOcean.everySteps * dt);
  };
  const lengths = [C, K * C, K * E, C, K * C, K * C, C];
  const state = lengths.map((n) => new Float64Array(n));
  let dirty = true, lastFrameTime = 0, lastFrameStep = 0;

  const model = { mesh, core, seaIce, radiation: radiationCpu, surface: surfaceCpu, geography, surfaceGeopotential: phis, state, time: 0, physics: true, moistOn: true, gpu, engine: 'gpu', get oceanCounter() { return oceanCounter; } };
  model.moist = { columnWater: moistCpu.columnWater, latentHeat: LATENT_HEAT, budget: moistCpu.budget, convectiveRain: moistCpu.convectiveRain, largeScaleRain: moistCpu.largeScaleRain };
  model.boundaryLayer = { depth: new Float64Array(C), mixingTop: new Float64Array(C), regime: new Float64Array(C), buoyancyFlux: new Float64Array(C) };
  model.oceanEngine = gpuOcean;

  function pushState() {
    gpu.upload(state);
    gpu.uploadPhysics({ land: geography ? Float32Array.from(geography.land, (l, i) => (l ? (geography.iceSheet && geography.iceSheet[i] ? 2 : 1) : 0)) : null, drag: dragCoefficients, soil: landCpu ? landCpu.soil : null, snow: landCpu ? landCpu.snow : null, vegetation: landCpu ? landCpu.vegetation : null, surface: landCpu ? landCpu.surface : null, concentration: seaIce.concentration, mlmSubsidence: radiationCpu.mlmSubsidence, mlmHeight: radiationCpu.mlmHeight, mlmGate: radiationCpu.mlmGate, convectiveRain: model.moist.convectiveRain, largeScaleRain: model.moist.largeScaleRain, meanAbsorbedSolar: radiationCpu.meanAbsorbedSolar, meanOutgoingLongwave: radiationCpu.meanOutgoingLongwave, meanPlanetaryAlbedo: radiationCpu.meanPlanetaryAlbedo, meanShortwaveCloudEffect: radiationCpu.meanShortwaveCloudEffect, meanLongwaveCloudEffect: radiationCpu.meanLongwaveCloudEffect, boundaryDepth: model.boundaryLayer.depth, mixingTop: model.boundaryLayer.mixingTop, regime: model.boundaryLayer.regime, buoyancyFlux: model.boundaryLayer.buoyancyFlux });
    gpu.clearFrame();
    if (gpuOcean) gpuOcean.initialize(state[3], state[6], { climatology: null });
    lastFrameTime = model.time;
    lastFrameStep = gpu.stepCount;
    dirty = false;
  }
  model.load = function load() { pushState(); };

  async function sync() {
    if (!dirty) return;
    const [arrays, [concentration, mean, height, gate, convective, largeScale, absorbed, outgoing, albedo, shortwaveEffect, longwaveEffect, depth, mixingTop, regime, buoyancy]] = await Promise.all([gpu.download(), readRanges(gpu.device, gpu.buffers.PH, ['CONC', 'MLMSUB', 'MLMH', 'MLMGATE', 'CONVMEAN', 'CONDMEAN', 'ASRMEAN', 'OLRMEAN', 'ALBMEAN', 'SWCREMEAN', 'LWCREMEAN', 'DEPTH', 'MIXTOP', 'REGIME', 'BUOY'].map((name) => ({ offset: gpu.layout.PH[name], length: C })))]);
    arrays.forEach((a, i) => state[i].set(a));
    seaIce.concentration.set(concentration);
    radiationCpu.mlmSubsidence.set(mean);
    radiationCpu.mlmHeight.set(height);
    radiationCpu.mlmGate.set(gate);
    model.moist.convectiveRain.set(convective);
    model.moist.largeScaleRain.set(largeScale);
    radiationCpu.meanAbsorbedSolar.set(absorbed);
    radiationCpu.meanOutgoingLongwave.set(outgoing);
    radiationCpu.meanPlanetaryAlbedo.set(albedo);
    radiationCpu.meanShortwaveCloudEffect.set(shortwaveEffect);
    radiationCpu.meanLongwaveCloudEffect.set(longwaveEffect);
    model.boundaryLayer.depth.set(depth);
    model.boundaryLayer.mixingTop.set(mixingTop);
    model.boundaryLayer.regime.set(regime);
    model.boundaryLayer.buoyancyFlux.set(buoyancy);
    dirty = false;
  }
  model.sync = sync;

  model.step = async function step(dt) {
    await gpu.stepModel(dt, model.time);
    model.time += dt;
    dirty = true;
  };
  /*
   * `count` steps recorded exactly as `step` queues them, submitted
   * together (see `batched` in core.gpu.js) and waited for once; resolves
   * to the number of submissions. `afterStep` runs after each step's
   * recording, and whatever it records through the core lands after that
   * step in the same submission. The mirrors stay lazy, as after `step`.
   */
  model.stepBatch = async function stepBatch(count, dt, afterStep = null) {
    const submissions = await gpu.batched(async () => {
      for (let n = 0; n < count; n++) {
        await gpu.stepModel(dt, model.time);
        model.time += dt;
        if (afterStep) afterStep();
      }
    });
    dirty = true;
    await gpu.device.queue.onSubmittedWorkDone();
    return submissions;
  };
  model.settle = () => gpu.device.queue.onSubmittedWorkDone();
  model.destroy = () => { for (const buffer of [...Object.values(gpu.buffers), ...(gpuOcean ? Object.values(gpuOcean.buffers) : [])]) buffer.destroy(); };

  /*
   * One frame of output: the named fields (see frames.module.js) at the
   * pressure level, and the diagnostics when asked, all computed on the
   * device and read back alone. The work is queued before this returns,
   * so a caller can queue steps behind it; the promise resolves once the
   * copies land. Every frame also advances the three-hour rain and the
   * runoff tally. The diagnostics' absorbed sunlight, its part absorbed
   * in the atmosphere, outgoing longwave and planetary albedo are means
   * over the steps since the last frame (the albedo their reflected over
   * their incoming sunlight), the last step's under `instantaneous`, and
   * the last step's alone when no step came between; with clearSkyPass
   * and a step between, the clear-sky absorbed sunlight and outgoing
   * longwave and the two cloud effects are means over the same steps.
   */
  let verticalLevel = null;
  model.beginFrame = function beginFrame({ level = 'surface', depth = 'surface', fields = [], diagnostics: summarize = false } = {}) {
    if (gpuOcean) gpuOcean.accumulateFreshwater(0);
    const time = model.time, interval = time - lastFrameTime, steps = gpu.stepCount - lastFrameStep;
    lastFrameTime = time;
    lastFrameStep = gpu.stepCount;
    const vertical = fields.includes('vertical');
    const keepVertical = vertical && verticalLevel === level ? Math.exp(-Math.max(0, interval) / VERTICAL_MEMORY) : 0;
    verticalLevel = vertical ? level : null;
    const rainScale = summarize && interval > 0 ? 86400 / interval : 0, meanScale = summarize && steps > 0 ? 1 / steps : 0;
    if (rainScale || meanScale) dirty = true;
    const atmosphere = gpu.frame({ pressure: level === 'surface' ? 0 : 100 * level, keep: Math.exp(-Math.max(0, interval) / RAIN_MEMORY), keepVertical, fields, diagnostics: summarize, rainScale, meanScale });
    if (gpuOcean) { gpuOcean.forgetAccumulated('rain'); gpuOcean.forgetAccumulated('runoff'); }
    const ocean = gpuOcean ? gpuOcean.frame({ fields, diagnostics: summarize, depth: depth === 'surface' ? 0 : Number(depth) }) : null;
    return Promise.all([atmosphere, ocean]).then(([a, o]) => {
      const out = { time, level, depth, fields: { ...a.fields, ...(o ? o.fields : {}) }, diagnostics: null };
      if (!a.sums) return out;
      const s = a.sums, area = s.area;
      if (landCpu) {
        let total = 0;
        for (let i = 0; i < C; i++) { landCpu.runoff[i] += a.runoff[i]; total += mesh.areaCell[i] * a.runoff[i]; }
        landCpu.budget.runoff += total;
      }
      const instantaneous = {
        absorbedSolar: s.absorbedSolar / area, atmosphereSolar: s.atmosphereSolar / area, outgoingLongwave: s.outgoingLongwave / area,
        planetaryAlbedo: s.insolation > 0 ? s.reflectedSolar / s.insolation : 0,
      };
      out.diagnostics = {
        mass: s.mass / area, meanSurfaceT: s.surfaceT / area, piMin: s.piMin, piMax: s.piMax, maxWind: s.maxWind,
        ...(steps > 0 ? {
          absorbedSolar: s.absorbedSum / area / steps, atmosphereSolar: s.atmosphereSum / area / steps, outgoingLongwave: s.outgoingSum / area / steps,
          planetaryAlbedo: s.insolationSum > 0 ? s.reflectedSum / s.insolationSum : 0,
        } : instantaneous),
        ...(radiationCpu.clearSkyPass && steps > 0 ? {
          clearAbsorbedSolar: s.clearAbsorbedSum / area / steps, clearOutgoingLongwave: s.clearOutgoingSum / area / steps,
          shortwaveCloudEffect: (s.absorbedSum - s.clearAbsorbedSum) / area / steps, longwaveCloudEffect: (s.clearOutgoingSum - s.outgoingSum) / area / steps,
        } : {}),
        instantaneous, sensibleHeat: s.sensibleHeat / area,
        evaporation: s.evaporation / area, latentHeat: LATENT_HEAT * s.evaporation / area,
        columnWater: s.water / area, columnCloud: s.cloud / area, precipitation: interval > 0 ? s.rain / area / interval : s.recentRain / area / RAIN_MEMORY,
        iceFraction: s.iceArea / area, iceThickness: s.iceArea > 0 ? s.iceVolume / s.iceArea : 0, surfaceAlbedo: s.albedo / area,
        ...(landCpu ? { landFraction: s.landArea / area, landMeanT: s.landArea > 0 ? s.landT / s.landArea : 0, snowFraction: s.landArea > 0 ? s.snowArea / s.landArea : 0, soilWater: s.landArea > 0 ? s.soil / s.landArea : 0, runoff: landCpu.budget.runoff / area } : {}),
        ...(o && o.diagnostics ? o.diagnostics : {}),
      };
      return out;
    });
  };
  model.diagnostics = async function diagnostics() { return (await model.beginFrame({ diagnostics: true })).diagnostics; };

  model.ocean = gpuOcean ? {
    initialize(surfaceT, iceField, options) {
      const started = gpuOcean.initialize(surfaceT, iceField, options);
      if (started) gpu.uploadSurfaceTemperature(surfaceT);
      return started;
    },
    load(saved, surfaceT, iceField) { gpuOcean.upload(saved, surfaceT, iceField); },
    async serialize({ restart = false } = {}) { return restart ? { ...(await gpuOcean.serialize()), ...(await gpuOcean.restartArrays()) } : gpuOcean.serialize(); },
  } : null;

  model.land = landCpu ? {
    soil: landCpu.soil, surface: landCpu.surface, snow: landCpu.snow, runoff: landCpu.runoff, vegetation: landCpu.vegetation, capacity: landCpu.capacity, land: geography.land, budget: landCpu.budget, albedo: landCpu.albedo, wetness: landCpu.wetness, water: landCpu.water, bucketCapacity: landCpu.bucketCapacity,
    initialize() { landCpu.initialize(); gpu.uploadLand({ soil: landCpu.soil, snow: landCpu.snow, vegetation: landCpu.vegetation, surface: landCpu.surface }); },
    load(saved, ice = state[6]) { landCpu.load(saved, ice); gpu.uploadLand({ soil: landCpu.soil, snow: landCpu.snow, vegetation: landCpu.vegetation, surface: landCpu.surface }); },
    async serialize() {
      const [soil, snow, vegetation, surface] = await readRanges(gpu.device, gpu.buffers.PH, [{ offset: gpu.layout.PH.SOIL, length: C }, { offset: gpu.layout.PH.SNOW, length: C }, { offset: gpu.layout.PH.VEG, length: C }, { offset: gpu.layout.PH.SURF, length: C }]);
      landCpu.soil.set(soil); landCpu.snow.set(snow); landCpu.vegetation.set(vegetation); landCpu.surface.set(surface);
      return landCpu.serialize();
    },
  } : null;

  return model;
}
