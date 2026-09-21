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
 */
export function createLandSurface(mesh, geography, {
  heatCapacity = 1e6, bucketCapacity = 150, wetnessThreshold = 0.75, albedo = 0.2, snowAlbedo = 0.55, fullSnow = 20,
  latentHeatFusion = 3.34e5, buffers = null,
} = {}) {
  const C = mesh.nCells;
  const shared = (name) => new Float64Array(buffers && buffers[name] ? buffers[name] : new SharedArrayBuffer(8 * C));
  const soil = shared('soil'), snow = shared('snow'), runoff = shared('runoff');
  const { land } = geography;
  const budget = { runoff: 0, melt: 0 };

  function overflow(i) {
    if (soil[i] > bucketCapacity) {
      const excess = soil[i] - bucketCapacity;
      runoff[i] += excess;
      budget.runoff += mesh.areaCell[i] * excess;
      soil[i] = bucketCapacity;
    }
  }

  function wetness(i) {
    return snow[i] > 0 ? 1 : Math.min(1, soil[i] / (wetnessThreshold * bucketCapacity));
  }

  function surfaceAlbedo(i) {
    return albedo + Math.min(1, snow[i] / fullSnow) * (snowAlbedo - albedo);
  }

  function update(i, surfaceT, flux, evaporation, dt) {
    surfaceT[i] += dt * flux[i] / heatCapacity;
    const fromSnow = Math.min(snow[i], evaporation * dt);
    snow[i] -= fromSnow;
    soil[i] = Math.max(0, soil[i] - (evaporation * dt - fromSnow));
    if (snow[i] > 0 && surfaceT[i] > MELTING_POINT) {
      const energy = (surfaceT[i] - MELTING_POINT) * heatCapacity;
      const melt = Math.min(snow[i], energy / latentHeatFusion);
      snow[i] -= melt;
      soil[i] += melt;
      budget.melt += mesh.areaCell[i] * melt;
      surfaceT[i] = MELTING_POINT + (energy - melt * latentHeatFusion) / heatCapacity;
    }
    overflow(i);
  }

  function deposit(i, rain, airTemperature) {
    if (airTemperature < MELTING_POINT) snow[i] += rain;
    else { soil[i] += rain; overflow(i); }
  }

  function initialize() {
    for (let i = 0; i < C; i++) { soil[i] = land[i] ? 0.5 * bucketCapacity : 0; snow[i] = 0; runoff[i] = 0; }
  }

  function water() {
    let total = 0;
    for (let i = 0; i < C; i++) if (land[i]) total += mesh.areaCell[i] * (soil[i] + snow[i] + runoff[i]);
    return total;
  }

  return {
    soil, snow, runoff, land, budget, heatCapacity, bucketCapacity, wetness, albedo: surfaceAlbedo, update, deposit, initialize, water,
    serialize() { return { soil: Float64Array.from(soil), snow: Float64Array.from(snow) }; },
    load(saved) { soil.set(saved.soil); snow.set(saved.snow); runoff.fill(0); },
    shared: { soil: soil.buffer, snow: snow.buffer, runoff: runoff.buffer },
  };
}
