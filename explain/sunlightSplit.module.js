/*
 * Where a beam of sunlight goes in a standard column: one column call with
 * the chosen ground and one over a black ground, whose reflection is the
 * air's own; the ground's reflection is the difference. The air is dry
 * unless `humidity` scales the profile's vapor (or gives q per layer).
 * Reflected by the air and by the ground, absorbed by ozone, by vapor, by
 * O₂ + CO₂ (with the near infrared the gases take on its way back up), by
 * aerosol and by the ground sum to the sunlight arriving at the top, which
 * is the full beam at the sun's elevation unless `beam` gives it.
 */
import { makeColumn, K, SEA_DIFFUSE_ALBEDO } from './radiativeColumn.module.js';
import { SOLAR_CONSTANT } from '../js/physics/radiation.module.js';
import { openWaterAlbedo } from '../js/physics/ice.module.js';

export const GROUNDS = { sea: 'sea', grassland: 0.2, forest: 0.13, soil: 0.3, snow: 0.85, ice: 0.62 };
export const LOWER = 10000;

export function lapseTropopause(T, height, pressure, lapse = 2e-3) {
  let below = null;
  for (let k = K - 2; k >= 0; k--) {
    const rate = (T[k + 1] - T[k]) / (height[k] - height[k + 1]), at = 0.5 * (height[k] + height[k + 1]);
    if (rate < lapse) {
      const z = below ? below.at + (below.rate - lapse) / (below.rate - rate) * (at - below.at) : at;
      const j = height[k + 1] > z ? k + 1 : k, i = j + 1;
      const p = Math.exp(Math.log(pressure[i]) + (z - height[i]) / (height[j] - height[i]) * (Math.log(pressure[j]) - Math.log(pressure[i])));
      return { height: z, pressure: p };
    }
    below = { rate, at };
  }
  return null;
}

export function createSplit(profile) {
  const lit = makeColumn({ profile, beam: SOLAR_CONSTANT, mu: 1, albedo: 'sea' });
  const black = makeColumn({ profile, beam: SOLAR_CONSTANT, mu: 1, albedo: 0 });
  const heating = new Float64Array(K), top = lapseTropopause(lit.T, lit.height, lit.pressure);
  return function split({ elevation = 90, ground = 'sea', ozone = true, humidity = 0, beam: arriving = null } = {}) {
    const mu = Math.sin(elevation * Math.PI / 180), beam = arriving ?? SOLAR_CONSTANT * mu, albedo = typeof ground === 'string' ? GROUNDS[ground] : ground;
    const vapor = typeof humidity === 'number' ? (humidity > 0 ? Float64Array.from(profile.humidity, (q) => humidity * q) : null) : humidity;
    lit.set({ beam, mu, albedo, gases: { ozone }, humidity: vapor });
    black.set({ beam, mu, albedo: 0, gases: { ozone }, humidity: vapor });
    const b = lit.budget, byAir = black.budget.reflected, { pressure, interfacePressure, height } = lit;
    heating.set(lit.heating.sunlight);
    let peak = -1, low = Infinity, high = -Infinity;
    for (let k = 0; k < K; k++) {
      if (interfacePressure[k] > 0 && (peak < 0 || heating[k] > heating[peak])) peak = k;
      if (height[k] < LOWER) { low = Math.min(low, heating[k]); high = Math.max(high, heating[k]); }
    }
    return {
      elevation, mu, arriving: beam, albedo: albedo === 'sea' ? [openWaterAlbedo(mu), SEA_DIFFUSE_ALBEDO] : [albedo, albedo],
      reflected: b.reflected, byAir, byGround: b.reflected - byAir,
      air: b.air, ozone: b.ozone, vapor: b.vapor, gases: b.oxygen + b.carbonDioxide + b.upwardNearInfrared, aerosol: b.aerosol, ground: b.ground, precipitableWater: b.precipitableWater,
      heating, pressure, interfacePressure, height,
      peak: { k: peak, rate: heating[peak], pressure: pressure[peak], height: height[peak] },
      lower: { low, high }, tropopause: top,
    };
  };
}
