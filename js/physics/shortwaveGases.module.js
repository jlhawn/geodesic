/*
 * Gas absorption of sunlight after CLIRAD-SW (Chou & Suarez 1999, NASA
 * TM-1999-104606 vol. 15, section 3). Every path is along the beam, in the
 * units of their tables: water vapour in g/cm2 and the gases in cm-atm (STP),
 * the water vapour, O2 and CO2 amounts scaled by (p/300 hPa)^0.8 and the
 * water vapour's also by 1 + 0.00135 (T - 240 K) (their eq. 3.5).
 *
 * Ozone: the eight ultraviolet and visible bands of their Table 3, each the
 * share s_b of the beam absorbed as 1 - exp(-k_b O) over the ozone path O.
 * Water vapour: in the visible band (0.4-0.7 um, s = 0.39081) the single
 * coefficient 0.00075 g-1 cm2; in the near infrared (0.7-10 um) the
 * ten-term k-distribution of their Table 2, its three bands' weights
 * summed, which absorbs sum_i h_i (1 - exp(-k_i w)) of the beam. Every water
 * vapour coefficient is multiplied by VAPOR_STRENGTH (1.48), fitted to RRTMG's
 * clear-sky atmospheric absorption over the tropical, midlatitude and
 * subarctic atmospheres (scripts/radiationBenchmark.mjs, VAPOR_FIT=1): the
 * tables come from HITRAN-96 lines cut off 10 cm-1 from their centres with
 * no continuum, and absorb 8 % less than RRTMG.
 * O2: 0.0633 of the beam (its 7600-8050, 12850-13190, 14310-14590 and
 * 15730-15930 cm-1 bands) absorbed as 1 - exp(-0.000145 sqrt(w)) (eq. 3.16).
 * CO2: the strong-line square root a sqrt(u) of the scaled path, with
 * a = CO2_COEFFICIENT giving the midlatitude summer column at 350 ppmv
 * under a 60 degree sun the 3.30 W/m2 their Table 7 gives (with the
 * overlap of water vapour), as scripts/radiationBenchmark.mjs checks.
 */
export const OZONE_SHARES = [0.00057, 0.00367, 0.00083, 0.00417, 0.00600, 0.00556, 0.05913, 0.39081];
export const OZONE_COEFFICIENTS = [30.47, 187.24, 301.92, 42.83, 7.09, 1.25, 0.0345, 0.0572];
export const VISIBLE_VAPOR = { share: 0.39081, coefficient: 0.00075 };
export const VAPOR_COEFFICIENTS = [0.0010, 0.0133, 0.0422, 0.1334, 0.4217, 1.3340, 5.6230, 31.620, 177.80, 1000.0];
export const VAPOR_WEIGHTS = [
  0.20673 + 0.08236 + 0.01074, 0.03497 + 0.01157 + 0.00360, 0.03011 + 0.01133 + 0.00411, 0.02260 + 0.01143 + 0.00421, 0.01336 + 0.01240 + 0.00389,
  0.00696 + 0.01258 + 0.00326, 0.00441 + 0.01381 + 0.00499, 0.00115 + 0.00650 + 0.00465, 0.00026 + 0.00244 + 0.00245, 0.00000 + 0.00094 + 0.00145,
];
export const OXYGEN = { share: 0.0633, coefficient: 0.000145, mixingRatio: 0.2095 };
export const CO2_COEFFICIENT = 1.7e-4;
export const VAPOR_STRENGTH = 1.48;
// STP_DEPTH: the depth (cm) at 273.15 K and 1 atm of 1 kg/m2 of air; OZONE_CM_ATM: kg/m2 of ozone per cm-atm
export const SCALING_PRESSURE = 30000, SCALING_EXPONENT = 0.8;
export const STP_DEPTH = 100 * 287.06 * 273.15 / 101325;
export const OZONE_CM_ATM = 2.1415e-2;

export function pressureScaling(p) {
  return Math.pow(p / SCALING_PRESSURE, SCALING_EXPONENT);
}

export function vaporScaling(p, T) {
  return pressureScaling(p) * (1 + 0.00135 * (T - 240));
}

export function ozoneAbsorptivity(path) {
  let a = 0;
  for (let b = 0; b < OZONE_SHARES.length; b++) a -= OZONE_SHARES[b] * Math.expm1(-OZONE_COEFFICIENTS[b] * path);
  return a;
}

export function visibleVaporAbsorptivity(path, strength = VAPOR_STRENGTH) {
  return -VISIBLE_VAPOR.share * Math.expm1(-VISIBLE_VAPOR.coefficient * strength * path);
}

export function nearInfraredVaporAbsorptivity(path, strength = VAPOR_STRENGTH) {
  let a = 0;
  for (let i = 0; i < VAPOR_COEFFICIENTS.length; i++) a -= VAPOR_WEIGHTS[i] * Math.expm1(-VAPOR_COEFFICIENTS[i] * strength * path);
  return a;
}

export function oxygenAbsorptivity(path) {
  return -OXYGEN.share * Math.expm1(-OXYGEN.coefficient * Math.sqrt(path));
}

export function carbonDioxideAbsorptivity(path) {
  return CO2_COEFFICIENT * Math.sqrt(path);
}
