/*
 * What a figure reads from one plumeColumn result. The plume's excess of
 * virtual temperature over the air in each layer is T_v B/g, with T_v the
 * sounding's T (1 + 0.608 q) before the step, so R Σ excess Δp/p over the
 * buoyant cloud layers is the model's CAPE wherever the step's first
 * condensation left the sounding as it was. The source air is the mass mean
 * of c_p T + g z and q over the layers within 50 hPa of the ground, taken to
 * the lowest layer and lifted dry to its Bolton LCL as moist.module.js lifts
 * it; from there the path follows the plume's own temperature at its
 * interfaces. Drying is the vapor removed counted as the latent heat it
 * carries, L/c_p per kg/kg. Layers run top first.
 */
import { liftingCondensationLevel, LATENT_HEAT, MOIST_DEFAULTS } from '../js/physics/moist.module.js';
import { VIRTUAL_FACTOR, GRAVITY, CP_DRY, R_DRY } from '../js/dynamics/sigmaCore.module.js';

const KAPPA = R_DRY / CP_DRY, DRY_STEPS = 12;

export function heightAt(r, p) {
  const { pressure, height } = r, K = pressure.length;
  if (p >= pressure[K - 1]) return height[K - 1] * Math.log(r.surfacePressure / p) / Math.log(r.surfacePressure / pressure[K - 1]);
  for (let k = K - 1; k > 0; k--) {
    if (p >= pressure[k - 1]) return height[k] + (height[k - 1] - height[k]) * Math.log(pressure[k] / p) / Math.log(pressure[k] / pressure[k - 1]);
  }
  return height[0];
}

export function readPlume(r) {
  const K = r.pressure.length, deep = r.type === 'deep';
  const excess = new Float64Array(K).fill(NaN), heating = new Float64Array(K), drying = new Float64Array(K);
  let peak = 0, heat = 0, dried = 0;
  for (let k = 0; k < K; k++) {
    const B = r.plume.buoyancy[k];
    if (Number.isFinite(B)) excess[k] = r.temperature[k] * (1 + VIRTUAL_FACTOR * r.vapor[k]) * B / GRAVITY;
    heating[k] = r.heating.convection[k] + r.heating.largeScale[k];
    drying[k] = -(r.moistening.convection[k] + r.moistening.largeScale[k]) * 1e-3 * LATENT_HEAT / CP_DRY;
    if (heating[k] > heating[peak]) peak = k;
    heat += r.mass[k] * CP_DRY * heating[k] / 86400;
    dried += r.mass[k] * CP_DRY * drying[k] / 86400;
  }

  const bottom = K - 1, ps = r.surfacePressure;
  let mass = 0, energy = 0, water = 0;
  for (let k = bottom; k >= 0; k--) {
    if (k < bottom && !(r.pressure[k] >= ps - MOIST_DEFAULTS.cumulusSourceDepth && r.pressure[k] > MOIST_DEFAULTS.shallowTop)) break;
    mass += r.mass[k];
    energy += r.mass[k] * (CP_DRY * r.temperature[k] + GRAVITY * r.height[k]);
    water += r.mass[k] * r.vapor[k];
  }
  const start = { pressure: r.pressure[bottom], temperature: (energy / mass - GRAVITY * r.height[bottom]) / CP_DRY };
  const vapor = water / mass, lcl = liftingCondensationLevel(start.temperature, vapor, start.pressure, KAPPA);
  const path = [];
  if (deep && lcl) {
    for (let n = 0; n <= DRY_STEPS; n++) {
      const p = start.pressure * (lcl.pressure / start.pressure) ** (n / DRY_STEPS);
      path.push({ pressure: p, temperature: start.temperature * (p / start.pressure) ** KAPPA });
    }
    for (let k = K - 1; k >= 0; k--) {
      const T = r.plume.interfaceTemperature[k];
      if (Number.isFinite(T) && r.interfaces[k] < lcl.pressure) path.push({ pressure: r.interfaces[k], temperature: T });
    }
  }

  let cape = 0;
  for (let k = 0; k < K; k++) if (excess[k] > 0 && r.pressure[k] <= r.base) cape += R_DRY * excess[k] * r.mass[k] * GRAVITY / r.pressure[k];

  return {
    excess, heating, drying, path, cape,
    source: { ...start, vapor, mass, lcl },
    cap: MOIST_DEFAULTS.cumulusBoundaryLoss * mass / r.dt,
    peak: { k: peak, rate: heating[peak], pressure: r.pressure[peak] },
    topHeight: r.top === null ? null : heightAt(r, r.top),
    column: { heating: heat, drying: dried, rain: LATENT_HEAT * r.rain / 86400 },
  };
}
