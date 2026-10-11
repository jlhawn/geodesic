/*
 * The model's bulk exchange between an open sea and the lowest layer, as
 * createSurfaceExchange's cell and radiation.column compute it for a sea
 * cell with no ice: COARE 3.5 roughness and Monin–Obukhov stability from
 * the bulk Richardson number of the virtual potential temperatures, then
 * E = ρ C_H U (q_s(T_s, p_s) − q), never below 0, and the sensible heat of
 * the dry static energy, ρ C_H U (c_p (T_s − T) − g z). The air fills a
 * layer shaped like the lowest bl36 layer (its midpoint near 20 m) whose
 * layer-mean Exner function Π makes the core's z = c_p θ_v (Π_s − Π)/g,
 * and holds rh times saturation at its midpoint pressure. U is the wind
 * alone unless `gust` gives a boundary layer depth (m above the layer),
 * which adds the model's free-convection gust over the sea.
 */
import { transfer, charnockRoughness, airViscosity, createBlend, EXCHANGE_DEFAULTS } from '../js/physics/exchange.module.js';
import { saturationHumidity, LATENT_HEAT } from '../js/physics/moist.module.js';
import { R_DRY, CP_DRY, GRAVITY, P0, VIRTUAL_FACTOR, sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';

export const LOWEST_HEIGHT = 20;
export const WATTS_PER_MM_DAY = LATENT_HEAT / 86400;
const KAPPA = R_DRY / CP_DRY;
const [UPPER, LOWER] = sigmaInterfaces('bl36').slice(-2);
const MIDPOINT_EXNER = Math.pow(0.5 * (UPPER + LOWER), KAPPA) * (1 + KAPPA) * (LOWER - UPPER) / (Math.pow(LOWER, 1 + KAPPA) - Math.pow(UPPER, 1 + KAPPA));
const pressureOf = (exner) => P0 * Math.pow(MIDPOINT_EXNER * exner, 1 / KAPPA);

export function bulkFluxes({ seaT, airT, rh, wind, ps = P0, z = LOWEST_HEIGHT, gust = null }) {
  if (!(wind > 0) && gust === null) throw new Error(`wind must be above 0 m/s without a gust, not ${wind}`);
  const surfaceExner = Math.pow(ps / P0, KAPPA);
  let exner = surfaceExner, qAir = 0;
  for (let n = 0; n < 4; n++) {
    qAir = rh * saturationHumidity(airT, pressureOf(exner));
    exner = surfaceExner / (1 + GRAVITY * z / (CP_DRY * airT * (1 + VIRTUAL_FACTOR * qAir)));
  }
  const pAir = pressureOf(exner);
  qAir = rh * saturationHumidity(airT, pAir);
  const qSea = saturationHumidity(seaT, ps);
  const airV = airT / exner * (1 + VIRTUAL_FACTOR * qAir), seaV = seaT / surfaceExner * (1 + VIRTUAL_FACTOR * qSea);
  const viscosity = airViscosity(airT - 273.15);
  const [seaGust, , floor] = EXCHANGE_DEFAULTS.convectiveGust;
  const blend = createBlend(EXCHANGE_DEFAULTS.blendingHeight), rough = { momentum: 0, heat: 0 }, c = { drag: 0, heat: 0, zeta: 0 };
  let speed = gust === null ? wind : Math.hypot(wind, floor), ri = 0;
  for (let pass = 0, passes = gust === null ? 1 : EXCHANGE_DEFAULTS.gustIterations + 1; pass < passes; pass++) {
    if (pass > 0) {
      const buoyancy = -ri * speed * speed * speed * c.heat / z;
      speed = Math.hypot(wind, buoyancy > 0 ? Math.max(floor, seaGust * Math.cbrt(buoyancy * Math.max(gust, z))) : floor);
    }
    charnockRoughness(speed, z, viscosity, EXCHANGE_DEFAULTS, 4, rough);
    blend.clear().add(1, rough.momentum, rough.heat).finish();
    ri = GRAVITY * z * (airV - seaV) / (0.5 * (airV + seaV) * speed * speed);
    transfer(ri, z, blend.momentum, blend.heat, EXCHANGE_DEFAULTS.iterations, c);
  }
  const neutral = transfer(0, z, blend.momentum, blend.heat);
  const density = pAir / (R_DRY * airT), conductance = density * c.heat * speed;
  const E = Math.max(0, conductance * (qSea - qAir));
  const SH = conductance * (CP_DRY * (seaT - airT) - GRAVITY * z);
  const LE = LATENT_HEAT * E;
  return {
    E, mmPerDay: 86400 * E, LE, SH, bowen: LE > 0 ? SH / LE : null,
    CH: c.heat, CD: c.drag, CHneutral: neutral.heat, Ri: ri, zeta: c.zeta, z0m: blend.momentum, z0h: blend.heat,
    speed, density, pAir, qAir, qSea, z,
  };
}
