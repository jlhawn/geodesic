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
 * step's precipitation in the adjustment phase, when it is known.
 *
 * With `vegetation` each land cell also carries a vegetation cover v
 * between 0 (bare ground) and 1 (dense forest) that the climate grows:
 * snow-free, it relaxes toward a goal set by how full the bucket is (0
 * below dryWetness of the capacity, 1 above wetWetness), over
 * growthTime when rising and declineTime when falling, and under snow it
 * decays toward 0 over snowDeclineTime. The bare-ground albedo runs from
 * bareAlbedo to vegetatedAlbedo and the bucket from minimumCapacity to
 * maximumCapacity (deeper roots) with v; water above a shrinking bucket
 * runs off. Without it the bucket is bucketCapacity and the albedo
 * `albedo` everywhere.
 */
export function createLandSurface(mesh, geography, {
  heatCapacity = 1e6, bucketCapacity = 150, wetnessThreshold = 0.75, albedo = 0.2, snowAlbedo = 0.55, fullSnow = 20,
  latentHeatFusion = 3.34e5, vegetation: vegetated = true, bareAlbedo = 0.35, vegetatedAlbedo = 0.13, minimumCapacity = 50,
  maximumCapacity = 500, dryWetness = 0.1, wetWetness = 0.6, growthTime = 180 * 86400, declineTime = 180 * 86400,
  snowDeclineTime = 720 * 86400, buffers = null,
} = {}) {
  const C = mesh.nCells;
  const shared = (name) => new Float64Array(buffers && buffers[name] ? buffers[name] : new SharedArrayBuffer(8 * C));
  const soil = shared('soil'), snow = shared('snow'), runoff = shared('runoff'), vegetation = shared('vegetation');
  const { land } = geography;
  const budget = { runoff: 0, melt: 0 };

  const capacity = (i) => (vegetated ? minimumCapacity + (maximumCapacity - minimumCapacity) * vegetation[i] : bucketCapacity);
  const bareGround = (i) => (vegetated ? bareAlbedo + (vegetatedAlbedo - bareAlbedo) * vegetation[i] : albedo);

  function overflow(i) {
    const cap = capacity(i);
    if (soil[i] > cap) {
      const excess = soil[i] - cap;
      runoff[i] += excess;
      budget.runoff += mesh.areaCell[i] * excess;
      soil[i] = cap;
    }
  }

  function wetness(i) {
    return snow[i] > 0 ? 1 : Math.min(1, soil[i] / (wetnessThreshold * capacity(i)));
  }

  function surfaceAlbedo(i) {
    const bare = bareGround(i);
    return bare + Math.min(1, snow[i] / fullSnow) * (snowAlbedo - bare);
  }

  function grow(i, dt) {
    if (snow[i] > 0) { vegetation[i] *= Math.exp(-dt / snowDeclineTime); return; }
    const cap = capacity(i);
    const goal = Math.min(1, Math.max(0, (Math.min(soil[i], cap) / cap - dryWetness) / (wetWetness - dryWetness)));
    vegetation[i] += (goal - vegetation[i]) * (1 - Math.exp(-dt / (goal > vegetation[i] ? growthTime : declineTime)));
  }

  function update(i, surfaceT, flux, evaporation, dt) {
    surfaceT[i] += dt * flux[i] / heatCapacity;
    const fromSnow = Math.min(snow[i], evaporation * dt);
    snow[i] -= fromSnow;
    surfaceT[i] -= latentHeatFusion * fromSnow / heatCapacity;
    soil[i] = Math.max(0, soil[i] - (evaporation * dt - fromSnow));
    if (snow[i] > 0 && surfaceT[i] > MELTING_POINT) {
      const energy = (surfaceT[i] - MELTING_POINT) * heatCapacity;
      const melt = Math.min(snow[i], energy / latentHeatFusion);
      snow[i] -= melt;
      soil[i] += melt;
      budget.melt += mesh.areaCell[i] * melt;
      surfaceT[i] = MELTING_POINT + (energy - melt * latentHeatFusion) / heatCapacity;
    }
    if (vegetated) grow(i, dt);
    overflow(i);
  }

  function deposit(i, rain, airTemperature) {
    if (airTemperature < MELTING_POINT) snow[i] += rain;
    else { soil[i] += rain; overflow(i); }
  }

  /*
   * A fresh land surface is fully vegetated with full buckets, so that
   * deserts emerge where the climate cannot keep the ground wet rather
   * than forests having to emerge from bare ground.
   */
  function initialize() {
    for (let i = 0; i < C; i++) {
      vegetation[i] = land[i] && vegetated ? 1 : 0;
      soil[i] = land[i] ? (vegetated ? capacity(i) : 0.5 * bucketCapacity) : 0;
      snow[i] = 0; runoff[i] = 0;
    }
  }

  function water() {
    let total = 0;
    for (let i = 0; i < C; i++) if (land[i]) total += mesh.areaCell[i] * (soil[i] + snow[i] + runoff[i]);
    return total;
  }

  /*
   * A saved land state without vegetation starts as initialize() would
   * where it is free of snow (vegetated, bucket full) and bare under snow.
   */
  function load(saved) {
    for (let i = 0; i < C; i++) {
      soil[i] = land[i] ? saved.soil[i] : 0;
      snow[i] = land[i] ? saved.snow[i] : 0;
      if (!land[i] || !vegetated) vegetation[i] = 0;
      else if (saved.vegetation) vegetation[i] = Math.min(1, Math.max(0, saved.vegetation[i]));
      else { vegetation[i] = snow[i] > 0 ? 0 : 1; if (snow[i] <= 0) soil[i] = capacity(i); }
    }
    runoff.fill(0);
  }

  return {
    soil, snow, runoff, vegetation, land, budget, heatCapacity, latentHeatFusion, bucketCapacity, capacity, wetness, albedo: surfaceAlbedo, update, deposit, initialize, water, load,
    serialize() { return { soil: Float64Array.from(soil), snow: Float64Array.from(snow), ...(vegetated ? { vegetation: Float64Array.from(vegetation) } : {}) }; },
    shared: { soil: soil.buffer, snow: snow.buffer, runoff: runoff.buffer, vegetation: vegetation.buffer },
  };
}
