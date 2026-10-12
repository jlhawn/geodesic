/*
 * Ekman's steady theory beside the model's slab, for an eastward wind
 * stress τ in N/m² (negative for a westward wind) and a Coriolis
 * parameter f in 1/s (negative in the south). Velocities are [east, north]
 * in m/s at a depth z ≤ 0 in meters, transports are in m²/s. The spiral is
 * not in the model, whose mixed layer moves as one slab.
 */
import { charnockRoughness, transfer, airViscosity, EXCHANGE_DEFAULTS } from '../js/physics/exchange.module.js';
import { sigmaInterfaces, R_DRY } from '../js/dynamics/sigmaCore.module.js';
import { SIDEREAL_OMEGA } from '../js/mesh.module.js';

export const RHO = 1025, AIR_DENSITY = 1.2, DAY = 86400;
const levels = sigmaInterfaces('bl36');
export const WIND_HEIGHT = R_DRY * 288 / 9.81 * Math.log(2 / (levels[levels.length - 2] + 1));

export const coriolis = (latitude) => 2 * SIDEREAL_OMEGA * Math.sin(latitude * Math.PI / 180);
export const inertialPeriod = (latitude) => 2 * Math.PI / Math.abs(coriolis(latitude));

export function ekmanSpiral({ tau, f, viscosity, rho = RHO }) {
  const depth = Math.sqrt(2 * viscosity / Math.abs(f)), surface = tau / (rho * Math.sqrt(viscosity * Math.abs(f))), turn = Math.sign(f);
  const at = (z) => {
    const decay = Math.exp(z / depth), angle = turn * (z / depth - Math.PI / 4);
    return [surface * decay * Math.cos(angle), surface * decay * Math.sin(angle)];
  };
  return { depth, surface, reversal: Math.PI * depth, at, transport: [0, -tau / (rho * f)] };
}

const roughness = { momentum: 0, heat: 0 }, coefficients = { drag: 0, heat: 0, zeta: 0 };
export function stressOfWind(speed, { height = WIND_HEIGHT, celsius = 15 } = {}) {
  charnockRoughness(speed, height, airViscosity(celsius), EXCHANGE_DEFAULTS, 4, roughness);
  transfer(0, height, roughness.momentum, roughness.heat, EXCHANGE_DEFAULTS.iterations, coefficients);
  return AIR_DENSITY * coefficients.drag * speed * speed;
}

export function windForStress(tau, options) {
  let low = 0.1, high = 60;
  for (let n = 0; n < 60; n++) { const mid = 0.5 * (low + high); if (stressOfWind(mid, options) < Math.abs(tau)) low = mid; else high = mid; }
  return 0.5 * (low + high);
}

export function tracks({ time, east, north, relativeEast, relativeNorth }) {
  const n = time.length, slab = new Float64Array(2 * n), beneath = new Float64Array(2 * n);
  for (let s = 1; s < n; s++) {
    const half = 0.5 * (time[s] - time[s - 1]);
    slab[2 * s] = slab[2 * s - 2] + half * (east[s] + east[s - 1]);
    slab[2 * s + 1] = slab[2 * s - 1] + half * (north[s] + north[s - 1]);
    beneath[2 * s] = beneath[2 * s - 2] + half * (east[s] - relativeEast[s] + east[s - 1] - relativeEast[s - 1]);
    beneath[2 * s + 1] = beneath[2 * s - 1] + half * (north[s] - relativeNorth[s] + north[s - 1] - relativeNorth[s - 1]);
  }
  return { slab, beneath };
}
