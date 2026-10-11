/*
 * One cloud in Part 6's column of the model's own radiation, measured
 * against the same column with the cloud taken out. A cloud is one layer's
 * in-cloud water path at a cover; radiativeColumn blends its overcast call
 * and its clear call by the cover, so every effect is linear in the cover
 * and a sweep at cover 1 holds them all. uniformCloud is the cloud the
 * model makes from a layer's humidity: moist.module.js's uniform
 * total-water distribution (condenseColumn), which warms the layer by
 * L q_c/c_p, and the radiation's cover of its condensate, √(q_c/b) floored
 * at 1%.
 */
import { makeColumn, layerPressures, interfacePressures, K, GLOBAL_BEAM, GLOBAL_MU } from './radiativeColumn.module.js';
import { SOLAR_CONSTANT } from '../js/physics/radiation.module.js';
import { cloudSaturation, criticalHumidityAt, uniformCondensate, uniformCover, LATENT_HEAT, MOIST_DEFAULTS } from '../js/physics/moist.module.js';
import { CP_DRY, GRAVITY } from '../js/dynamics/sigmaCore.module.js';

export const HIGHEST = 13;
export const LAYERS = Array.from({ length: K - HIGHEST }, (_, n) => K - 1 - n);
export const PATHS = [5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 250, 300, 350, 400].map((g) => g / 1000);
export const COVER_FLOOR = 0.01;
const LOW_SUN = Math.sin(15 * Math.PI / 180);
export const SUNS = {
  average: { beam: GLOBAL_BEAM, mu: GLOBAL_MU },
  overhead: { beam: SOLAR_CONSTANT, mu: 1 },
  low: { beam: SOLAR_CONSTANT * LOW_SUN, mu: LOW_SUN },
  night: { beam: 0, mu: GLOBAL_MU },
};

export function createSweep(profile) {
  const column = makeColumn({ profile, humidity: profile.humidity }), rows = LAYERS.length, cols = PATHS.length, total = rows * cols;
  const shortwave = new Float64Array(total), longwave = new Float64Array(total);
  const order = [total - 1, cols - 1, ...Array.from({ length: total }, (_, n) => n).filter((n) => n !== total - 1 && n !== cols - 1)];
  let next = 0;
  function start({ humidity = profile.humidity, sun = SUNS.average } = {}) {
    column.set({ humidity, ...sun, cloud: null });
    shortwave.fill(NaN);
    longwave.fill(NaN);
    next = 0;
  }
  function advance(budget = Infinity) {
    const begun = performance.now();
    while (next < total && performance.now() - begun < budget) {
      const n = order[next++], r = Math.floor(n / cols), c = n % cols;
      column.set({ cloud: { pressure: column.pressure[LAYERS[r]], path: PATHS[c], cover: 1 } });
      shortwave[n] = column.budget.shortwaveCloudEffect;
      longwave[n] = column.budget.longwaveCloudEffect;
    }
    return next === total;
  }
  function largest() {
    let most = 0;
    for (const n of [total - 1, cols - 1]) if (Number.isFinite(shortwave[n])) most = Math.max(most, Math.abs(shortwave[n]), Math.abs(longwave[n]), Math.abs(shortwave[n] + longwave[n]));
    return most;
  }
  start();
  return { rows, cols, shortwave, longwave, start, advance, largest, get done() { return next === total; }, get count() { return next; } };
}

export function zeroCrossings(values, rows, cols, paths = PATHS) {
  return Array.from({ length: rows }, (_, r) => {
    const found = [];
    for (let c = 1; c < cols; c++) {
      const a = values[r * cols + c - 1], b = values[r * cols + c];
      if (!(Number.isFinite(a) && Number.isFinite(b)) || (a > 0) === (b > 0) || a === b) continue;
      const t = a / (a - b);
      found.push(Math.exp(Math.log(paths[c - 1]) + t * (Math.log(paths[c]) - Math.log(paths[c - 1]))));
    }
    return found;
  });
}

export function uniformCloud({ temperature, pressure, surfacePressure, humidity, mass = 0 }, options = MOIST_DEFAULTS) {
  const saturated = cloudSaturation(temperature, pressure, options.iceSaturation, options.liquidTemperature, options.iceTemperature, {});
  const a = 1 / (1 + LATENT_HEAT * saturated.slope / CP_DRY);
  const critical = criticalHumidityAt(pressure, surfacePressure, options.surfaceCriticalHumidity, options.topCriticalHumidity, options.criticalExponent);
  const halfWidth = a * (1 - critical) * saturated.qs, total = humidity * saturated.qs;
  const condensate = uniformCondensate(a * (total - saturated.qs), halfWidth);
  const cover = condensate > 0 ? Math.min(1, Math.max(COVER_FLOOR, uniformCover(condensate, halfWidth))) : 0;
  return {
    critical, spread: 1 - critical, saturation: saturated.qs, liquid: saturated.liquid, a, halfWidth, total, condensate, cover,
    vapor: total - condensate, warming: LATENT_HEAT * condensate / CP_DRY, gridPath: condensate * mass, path: cover > 0 ? condensate * mass / cover : 0,
  };
}

export function cloudFromHumidity(profile, k, humidity, vapor = profile.humidity) {
  const ps = profile.surfacePressure, pressure = layerPressures(ps)[k], edges = interfacePressures(ps);
  const made = uniformCloud({ temperature: profile.temperature[k], pressure, surfacePressure: ps, humidity, mass: (edges[k + 1] - edges[k]) / GRAVITY });
  const q = Float64Array.from(vapor ?? new Float64Array(K)), temperature = Float64Array.from(profile.temperature);
  q[k] = made.vapor;
  temperature[k] += made.warming;
  return { made, humidity: q, temperature, cloud: made.cover > 0 ? { pressure, path: made.path, cover: made.cover } : null };
}
