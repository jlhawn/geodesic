import { NEAR_INFRARED_TABLE } from './nearInfraredTable.module.js';

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

/*
 * The near infrared's gases in bands (radiation option nearInfrared 'bands'):
 * the k-distributions of nearInfraredTable.module.js, fitted by
 * scripts/nearInfraredFit.mjs to RRTMG_SW in its bands 16-24 and 29
 * (820-16000 cm-1). Band b carries the share s_b of the beam; each of its
 * gases (water vapor, CO2, CH4, O2) transmits
 * w0 exp(-m S) + sum_j w_j exp(-m (tau_j + S)) of it along a path of
 * magnification m, tau_j summing over the layers above
 * k_j u (p_V/p_ref)^n_j exp(t_j (T - 250)) for each layer's mass u (kg/m2)
 * of the gas, p_V = sqrt(p^2 + p_D^2), and S, water vapor's self continuum,
 * summing self u e/p_ref (296/T)^selfT with e the vapor pressure; the
 * band's gases' transmissions multiply. The visible band's water vapor
 * takes the strength VISIBLE_VAPOR_BANDS in place of VAPOR_STRENGTH.
 */
export const NEAR_INFRARED_BANDS = NEAR_INFRARED_TABLE.bands;
export const NEAR_INFRARED_P_REF = NEAR_INFRARED_TABLE.pRef;
export const VISIBLE_VAPOR_BANDS = NEAR_INFRARED_TABLE.visibleVaporStrength;
export const NEAR_INFRARED_MOLAR = { co2: 44.01 / 28.964, ch4: 16.04 / 28.964, o2: 31.999 / 28.964 };

let bandScratch = null;

/*
 * air: per layer from the top, the masses (kg/m2) of water vapor (`vapour`)
 * and dry air (`dry`), the mid-pressure p (Pa), temperature T and vapor
 * pressure e (Pa); mixing: the mass of CO2, CH4 and O2 per mass of dry air;
 * vaporFactor scales water vapor's absorptivity in every band. Fills
 * taken[k], the share of the beam absorbed above the bottom of layer k, and
 * upward[k], the share absorbed in layer k by the beam that, reflected whole
 * at the surface, climbs back at the path D per unit of the vertical; split,
 * when given, gets each gas's part of the beam's absorption (a band's in
 * proportion to its gases' absorptivities alone). Returns the upward total.
 */
export function nearInfraredBands(air, m, D, mixing, vaporFactor, taken, upward, split = null) {
  const { K, vapour, dry, p, T, e } = air;
  if (!bandScratch || bandScratch.down.length < K + 1) bandScratch = { down: new Float64Array(K + 1), up: new Float64Array(K + 1), increments: new Float64Array(8 * K), total: new Float64Array(8), sum: new Float64Array(8) };
  const { down, up, increments, total, sum } = bandScratch;
  taken.fill(0, 0, K);
  upward.fill(0, 0, K);
  if (split) { split.h2o = 0; split.co2 = 0; split.ch4 = 0; split.o2 = 0; }
  let loss = 0;
  for (const band of NEAR_INFRARED_BANDS) {
    down.fill(1, 0, K + 1);
    up.fill(1, 0, K + 1);
    let absorptivities = 0;
    const alone = {};
    for (const gas in band.gases) {
      const par = band.gases[gas], J = par.w.length, isVapor = gas === 'h2o', self = isVapor ? par.self : 0;
      let w0 = 1;
      for (let j = 0; j < J; j++) w0 -= par.w[j];
      total.fill(0);
      for (let k = 0; k < K; k++) {
        const u = isVapor ? vapour[k] : mixing[gas] * dry[k], lnV = Math.log(Math.sqrt(p[k] * p[k] + par.pD * par.pD) / NEAR_INFRARED_P_REF), dT = T[k] - 250;
        for (let j = 0; j < J; j++) { increments[j * K + k] = par.k[j] * u * Math.exp(par.n[j] * lnV + par.t[j] * dT); total[j] += increments[j * K + k]; }
        increments[J * K + k] = self > 0 ? self * u * e[k] / NEAR_INFRARED_P_REF * Math.exp(par.selfT * Math.log(296 / T[k])) : 0;
        total[J] += increments[J * K + k];
      }
      sum.fill(0);
      for (let l = 0; l <= K; l++) {
        if (l > 0) for (let j = 0; j <= J; j++) sum[j] += increments[j * K + l - 1];
        const S = sum[J], ST = total[J];
        let through = w0 * Math.exp(-m * S), back = w0 * Math.exp(-(m * ST + D * (ST - S)));
        for (let j = 0; j < J; j++) {
          through += par.w[j] * Math.exp(-m * (sum[j] + S));
          back += par.w[j] * Math.exp(-(m * (total[j] + ST) + D * (total[j] - sum[j] + ST - S)));
        }
        if (isVapor) { through = 1 - vaporFactor * (1 - through); back = 1 - vaporFactor * (1 - back); }
        down[l] *= through;
        up[l] *= back;
        if (l === K) { alone[gas] = 1 - through; absorptivities += 1 - through; }
      }
    }
    for (let k = 0; k < K; k++) {
      taken[k] += band.share * (1 - down[k + 1]);
      upward[k] += band.share * (up[k + 1] - up[k]);
    }
    loss += band.share * (up[K] - up[0]);
    if (split && absorptivities > 0) for (const gas in alone) split[gas] += band.share * (1 - down[K]) * alone[gas] / absorptivities;
  }
  return loss;
}

let shieldedScratch = null;

/*
 * The near infrared's gases around a cloud top at level `top` (the interface
 * above layer `top`, 0 to K), by the k-distributions of nearInfraredBands, for
 * the radiation option vaporShielding 'cloudTop': light crosses the layers above
 * the top along the beam's path m and the layers below along the diffuse path D.
 * Fills, as shares of the beam, for each layer k at or below the top, down[k], the
 * share absorbed in k by light that has passed the cloud, and up[k], the share
 * absorbed in k by that light once the surface has reflected it whole and it
 * climbs back to the cloud's base, and for each layer k above the top, climb[k],
 * the share absorbed in k by light the cloud has reflected whole, climbing along D.
 */
export function nearInfraredShielded(air, m, D, mixing, vaporFactor, top, shielded) {
  const { K, vapour, dry, p, T, e } = air;
  if (!shieldedScratch || shieldedScratch.passed.length < K + 1) shieldedScratch = { passed: new Float64Array(K + 1), round: new Float64Array(K + 1), reflected: new Float64Array(K + 1), increments: new Float64Array(8 * K), total: new Float64Array(8), above: new Float64Array(8), sum: new Float64Array(8) };
  const { passed, round, reflected, increments, total, above, sum } = shieldedScratch;
  const { down, up, climb } = shielded;
  down.fill(0, 0, K);
  up.fill(0, 0, K);
  climb.fill(0, 0, K);
  for (const band of NEAR_INFRARED_BANDS) {
    passed.fill(1, 0, K + 1);
    round.fill(1, 0, K + 1);
    reflected.fill(1, 0, K + 1);
    for (const gas in band.gases) {
      const par = band.gases[gas], J = par.w.length, isVapor = gas === 'h2o', self = isVapor ? par.self : 0;
      let w0 = 1;
      for (let j = 0; j < J; j++) w0 -= par.w[j];
      total.fill(0);
      above.fill(0);
      for (let k = 0; k < K; k++) {
        const u = isVapor ? vapour[k] : mixing[gas] * dry[k], lnV = Math.log(Math.sqrt(p[k] * p[k] + par.pD * par.pD) / NEAR_INFRARED_P_REF), dT = T[k] - 250;
        for (let j = 0; j < J; j++) increments[j * K + k] = par.k[j] * u * Math.exp(par.n[j] * lnV + par.t[j] * dT);
        increments[J * K + k] = self > 0 ? self * u * e[k] / NEAR_INFRARED_P_REF * Math.exp(par.selfT * Math.log(296 / T[k])) : 0;
        for (let j = 0; j <= J; j++) {
          total[j] += increments[j * K + k];
          if (k < top) above[j] += increments[j * K + k];
        }
      }
      const SA = above[J], ST = total[J];
      sum.fill(0);
      for (let l = 0; l <= K; l++) {
        if (l > 0) for (let j = 0; j <= J; j++) sum[j] += increments[j * K + l - 1];
        const S = sum[J];
        if (l >= top) {
          let through = w0 * Math.exp(-(m * SA + D * (S - SA))), back = w0 * Math.exp(-(m * SA + D * (2 * ST - SA - S)));
          for (let j = 0; j < J; j++) {
            through += par.w[j] * Math.exp(-(m * (above[j] + SA) + D * (sum[j] - above[j] + S - SA)));
            back += par.w[j] * Math.exp(-(m * (above[j] + SA) + D * (2 * total[j] - above[j] - sum[j] + 2 * ST - SA - S)));
          }
          if (isVapor) { through = 1 - vaporFactor * (1 - through); back = 1 - vaporFactor * (1 - back); }
          passed[l] *= through;
          round[l] *= back;
        }
        if (l <= top) {
          let climbing = w0 * Math.exp(-(m * SA + D * (SA - S)));
          for (let j = 0; j < J; j++) climbing += par.w[j] * Math.exp(-(m * (above[j] + SA) + D * (above[j] - sum[j] + SA - S)));
          if (isVapor) climbing = 1 - vaporFactor * (1 - climbing);
          reflected[l] *= climbing;
        }
      }
    }
    for (let k = 0; k < K; k++) {
      if (k >= top) {
        down[k] += band.share * (passed[k] - passed[k + 1]);
        up[k] += band.share * (round[k + 1] - round[k]);
      } else climb[k] += band.share * (reflected[k + 1] - reflected[k]);
    }
  }
}
