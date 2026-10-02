import { MELTING_POINT } from './ice.module.js';

/*
 * The land surface: a skin of small heat capacity, a bucket of soil
 * water in the manner of Manabe (1969) that limits evaporation as it
 * dries and spills what it cannot hold into runoff, and a snow cover in
 * water equivalent that precipitation builds when the lowest air is
 * below freezing and the surface energy melts. Evaporation draws on the
 * snow while there is any, else on the bucket; snow raises the albedo
 * toward snowAlbedo over fullSnow kg/m². update() applies the surface
 * flux and the evaporation of the physics phase; deposit() adds the
 * step's precipitation in the adjustment phase, when it is known. The
 * ocean cells of `snow` belong to the sea ice, which keeps its snow there.
 *
 * With `vegetation` each land cell also carries a vegetation cover v
 * between 0 (bare ground) and 1 (dense forest) that the climate grows:
 * snow-free, it relaxes toward a goal set by how full the bucket is (0
 * below dryWetness of the capacity, 1 above wetWetness), over
 * growthTime when rising and declineTime when falling, and under snow it
 * decays toward 0 over snowDeclineTime. The snow-free albedo runs from
 * the bare soil's to vegetatedAlbedo with v; with soilDarkening (the
 * default) the bare soil darkens linearly with the bucket's fill, as in
 * BATS and CLM, from bareAlbedo dry to wetSoilAlbedo at darkeningWetness
 * of the capacity and beyond, else it is bareAlbedo. The bucket holds a fixed
 * rootZoneCapacity whatever the cover, because a soil keeps its water
 * capacity when its plants die and a browned region can therefore
 * regreen when the rain returns. Without it the bucket is bucketCapacity and the albedo
 * `albedo` everywhere. A cell of the geography's `iceSheet` grows no
 * vegetation and keeps iceSheetAlbedo whatever lies on it.
 *
 * With vegetation the soil has two stores. Rain fills a surface layer
 * of surfaceCapacity kg/m² first; what it cannot hold infiltrates the
 * root zone (the bucket), a share (soil/capacity)⁴ of it running off,
 * and the surface layer seeps into the root zone over percolationTime.
 * Bare ground evaporates from the surface layer alone, so a desert
 * dries within a day of rain; the cover transpires from the root zone
 * through stomata, at the aerodynamic rate times 1/(1 + r_s g_a) with
 * r_s = stomatalResistance divided by a warmth that rises from 0 at
 * growthColdest to 1 at growthWarmest, so cold or dry roots close
 * them. wetness(i, g_a, T) gives the radiation column that factor and
 * remembers how much of it was bare ground, which update() draws from
 * the surface layer. Growth toward the goal needs the same warmth;
 * decline does not.
 */
export function createLandSurface(mesh, geography, {
  heatCapacity = 1e6, bucketCapacity = 150, wetnessThreshold = 0.75, albedo = 0.2, snowAlbedo = 0.55, fullSnow = 20,
  latentHeatFusion = 3.34e5, vegetation: vegetated = true, bareAlbedo = 0.30, vegetatedAlbedo = 0.13, rootZoneCapacity = 300, dryWetness = 0.1, wetWetness = 0.6, growthTime = 180 * 86400, declineTime = 365 * 86400,
  snowDeclineTime = 720 * 86400, iceSheetAlbedo = 0.8, surfaceCapacity = 15, percolationTime = 86400, stomatalResistance = 70,
  growthColdest = 278.15, growthWarmest = 288.15, soilDarkening = true, wetSoilAlbedo = 0.15, darkeningWetness = 0.5, buffers = null,
} = {}) {
  const C = mesh.nCells;
  const shared = (name) => new Float64Array(buffers && buffers[name] ? buffers[name] : new SharedArrayBuffer(8 * C));
  const soil = shared('soil'), snow = shared('snow'), runoff = shared('runoff'), vegetation = shared('vegetation'), surface = shared('surface');
  const bareShare = new Float64Array(C);
  const warmth = (t) => Math.min(1, Math.max(0, (t - growthColdest) / (growthWarmest - growthColdest)));
  const { land, iceSheet = null } = geography;
  const onIceSheet = (i) => iceSheet !== null && iceSheet[i] > 0;
  const budget = { runoff: 0, melt: 0 };

  const capacity = () => (vegetated ? rootZoneCapacity : bucketCapacity);
  const soilAlbedo = (i) => (soilDarkening ? bareAlbedo - (bareAlbedo - wetSoilAlbedo) * Math.min(1, soil[i] / (darkeningWetness * capacity(i))) : bareAlbedo);
  const bareGround = (i) => { if (!vegetated) return albedo; const s = soilAlbedo(i); return s + (vegetatedAlbedo - s) * vegetation[i]; };

  function overflow(i) {
    const cap = capacity(i);
    if (soil[i] > cap) {
      const excess = soil[i] - cap;
      runoff[i] += excess;
      budget.runoff += mesh.areaCell[i] * excess;
      soil[i] = cap;
    }
  }

  function wetness(i, aero = 0, temperature = growthWarmest) {
    if (snow[i] > 0) { bareShare[i] = 0; return 1; }
    const roots = Math.min(1, soil[i] / (wetnessThreshold * capacity(i)));
    if (!vegetated) { bareShare[i] = 0; return roots; }
    const bare = (1 - vegetation[i]) * Math.min(1, surface[i] / surfaceCapacity);
    const canopy = vegetation[i] * roots / (1 + stomatalResistance * aero / Math.max(0.05, warmth(temperature)));
    const total = bare + canopy;
    bareShare[i] = total > 0 ? bare / total : 0;
    return total;
  }

  function surfaceAlbedo(i) {
    if (onIceSheet(i)) return iceSheetAlbedo;
    const bare = bareGround(i);
    return bare + Math.min(1, snow[i] / fullSnow) * (snowAlbedo - bare);
  }

  function grow(i, dt, temperature) {
    if (onIceSheet(i)) { vegetation[i] = 0; return; }
    if (snow[i] > 0) { vegetation[i] *= Math.exp(-dt / snowDeclineTime); return; }
    const cap = capacity(i);
    const goal = Math.min(1, Math.max(0, (Math.min(soil[i], cap) / cap - dryWetness) / (wetWetness - dryWetness)));
    if (goal > vegetation[i]) vegetation[i] += (goal - vegetation[i]) * (1 - Math.exp(-dt * warmth(temperature) / growthTime));
    else vegetation[i] += (goal - vegetation[i]) * (1 - Math.exp(-dt / declineTime));
  }

  function update(i, surfaceT, flux, evaporation, dt) {
    const temperature = surfaceT[i];
    surfaceT[i] += dt * flux[i] / heatCapacity;
    let left = evaporation * dt;
    const fromSnow = Math.min(snow[i], left);
    snow[i] -= fromSnow; left -= fromSnow;
    surfaceT[i] -= latentHeatFusion * fromSnow / heatCapacity;
    if (vegetated) {
      const fromSurface = Math.min(surface[i], left * bareShare[i]);
      surface[i] -= fromSurface; left -= fromSurface;
      const seep = surface[i] * (1 - Math.exp(-dt / percolationTime));
      surface[i] -= seep; soil[i] += seep;
    }
    soil[i] = Math.max(0, soil[i] - left);
    if (snow[i] > 0 && surfaceT[i] > MELTING_POINT) {
      const energy = (surfaceT[i] - MELTING_POINT) * heatCapacity;
      const melt = Math.min(snow[i], energy / latentHeatFusion);
      snow[i] -= melt;
      soil[i] += melt;
      budget.melt += mesh.areaCell[i] * melt;
      surfaceT[i] = MELTING_POINT + (energy - melt * latentHeatFusion) / heatCapacity;
    }
    if (vegetated) grow(i, dt, temperature);
    overflow(i);
  }

  function deposit(i, rain, airTemperature) {
    if (airTemperature < MELTING_POINT) { snow[i] += rain; return; }
    if (!vegetated) { soil[i] += rain; overflow(i); return; }
    surface[i] += rain;
    if (surface[i] > surfaceCapacity) {
      const infiltration = surface[i] - surfaceCapacity;
      surface[i] = surfaceCapacity;
      const shed = infiltration * Math.pow(Math.min(1, soil[i] / capacity(i)), 4);
      runoff[i] += shed;
      budget.runoff += mesh.areaCell[i] * shed;
      soil[i] += infiltration - shed;
      overflow(i);
    }
  }

  /*
   * A fresh land surface is fully vegetated with full buckets, so that
   * deserts emerge where the climate cannot keep the ground wet rather
   * than forests having to emerge from bare ground.
   */
  function initialize() {
    for (let i = 0; i < C; i++) {
      vegetation[i] = land[i] && vegetated && !onIceSheet(i) ? 0.5 : 0;
      soil[i] = land[i] ? 0.5 * capacity(i) : 0;
      snow[i] = 0; runoff[i] = 0; surface[i] = 0;
    }
  }

  function water() {
    let total = 0;
    for (let i = 0; i < C; i++) if (land[i]) total += mesh.areaCell[i] * (soil[i] + surface[i] + snow[i] + runoff[i]);
    return total;
  }

  /*
   * A saved land state without vegetation starts as initialize() would
   * where it is free of snow (vegetated, bucket full) and bare under snow.
   * A sea cell keeps the saved snow only where `ice` has ice to hold it.
   */
  function load(saved, ice = null) {
    for (let i = 0; i < C; i++) {
      soil[i] = land[i] ? saved.soil[i] : 0;
      surface[i] = land[i] && saved.surface ? Math.min(surfaceCapacity, saved.surface[i]) : 0;
      snow[i] = land[i] || (ice && ice[i] > 0) ? saved.snow[i] : 0;
      if (!land[i] || !vegetated || onIceSheet(i)) vegetation[i] = 0;
      else if (saved.vegetation) vegetation[i] = Math.min(1, Math.max(0, saved.vegetation[i]));
      else { vegetation[i] = snow[i] > 0 ? 0 : 1; if (snow[i] <= 0) soil[i] = capacity(i); }
    }
    runoff.fill(0);
  }

  return {
    soil, surface, snow, runoff, vegetation, land, budget, heatCapacity, latentHeatFusion, bucketCapacity, capacity, wetness, albedo: surfaceAlbedo, update, deposit, initialize, water, load,
    serialize() { return { soil: Float64Array.from(soil), snow: Float64Array.from(snow), ...(vegetated ? { vegetation: Float64Array.from(vegetation), surface: Float64Array.from(surface) } : {}) }; },
    shared: { soil: soil.buffer, surface: surface.buffer, snow: snow.buffer, runoff: runoff.buffer, vegetation: vegetation.buffer },
  };
}
