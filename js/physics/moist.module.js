import { R_DRY } from '../dynamics/sigmaCore.module.js';

export const LATENT_HEAT = 2.5e6;
export const EPSILON = 0.622;
export const R_VAPOR = R_DRY / EPSILON;

export function saturationVaporPressure(T) {
  return 611.2 * Math.exp(17.67 * (T - 273.15) / (T - 29.65));
}

export function saturationHumidity(T, p) {
  const es = saturationVaporPressure(T);
  const dry = p - (1 - EPSILON) * es;
  return dry > 0 ? EPSILON * es / dry : 1;
}

/*
 * Moist physics for the sigma core, applied to the state after each
 * step: large-scale condensation of supersaturation (rained out at
 * once, latent heat to the layer), the simplified Betts–Miller
 * convection of Frierson (2007) — a conditionally unstable column
 * relaxes over relaxationTime toward the moist adiabat of its lowest
 * layer and a fixed relative humidity, with the reference temperature
 * shifted so the column's enthalpy change equals the latent heat of
 * the rain it produces, or with no rain when the column would have to
 * moisten — and a filler that removes negative humidity by borrowing
 * from the layer below. Precipitation accumulates per cell (kg/m²);
 * the budget sums are area-weighted masses (kg).
 */
export function createMoistPhysics(mesh, core, { latentHeat = LATENT_HEAT, relaxationTime = 7200, referenceHumidity = 0.7, buffers = null } = {}) {
  const { K, C, dSigma, sigmaMid, cp, R, g, kappa, exnerLayer } = core.diagnostics;
  const precipBuffer = buffers && buffers.precipitation ? buffers.precipitation : new SharedArrayBuffer(8 * C);
  const precipitation = new Float64Array(precipBuffer);
  const T = new Float64Array(K), p = new Float64Array(K), dp = new Float64Array(K), Tref = new Float64Array(K), qref = new Float64Array(K);
  const budget = { condensation: 0, convection: 0, lost: 0 };

  function moistLapse(temperature, pressure) {
    const qs = saturationHumidity(temperature, pressure);
    return (R * temperature + latentHeat * qs) / (cp + latentHeat * latentHeat * EPSILON * qs / (R * temperature * temperature));
  }

  function condenseColumn(i, pi, theta, q) {
    let rain = 0;
    for (let k = 0; k < K; k++) {
      const idx = k * C + i;
      const ex = exnerLayer[idx];
      const temperature = theta[idx] * ex;
      const pressure = pi[i] * sigmaMid[k];
      const qs = saturationHumidity(temperature, pressure);
      if (q[idx] <= qs) continue;
      const slope = qs * latentHeat / (R_VAPOR * temperature * temperature);
      const removed = (q[idx] - qs) / (1 + latentHeat * slope / cp);
      q[idx] -= removed;
      theta[idx] += latentHeat * removed / (cp * ex);
      rain += pi[i] * dSigma[k] / g * removed;
    }
    return rain;
  }

  /*
   * The reference profile: dry adiabat from the lowest layer up to its
   * lifting condensation level, the moist adiabat above. The convecting
   * column reaches the highest level above the LCL where the parcel is
   * still warmer than the air (its level of zero buoyancy; stable layers
   * in between are ignored, as in Frierson 2007). Returns the index of
   * that top layer, or -1 when the parcel never saturates or is never
   * buoyant above its LCL.
   */
  function referenceProfile(i, pi, theta, q) {
    const bottom = K - 1;
    for (let k = 0; k < K; k++) {
      const idx = k * C + i;
      T[k] = theta[idx] * exnerLayer[idx];
      p[k] = pi[i] * sigmaMid[k];
      dp[k] = pi[i] * dSigma[k];
    }
    const Tp = T[bottom], pb = p[bottom];
    const qp = Math.min(q[bottom * C + i], saturationHumidity(Tp, pb));
    const dryT = (pressure) => Tp * Math.pow(pressure / pb, kappa);
    if (saturationHumidity(dryT(p[0]), p[0]) > qp) return -1;
    let lo = Math.log(p[0]), hi = Math.log(pb);
    for (let n = 0; n < 40; n++) {
      const mid = 0.5 * (lo + hi), pressure = Math.exp(mid);
      if (saturationHumidity(dryT(pressure), pressure) > qp) hi = mid; else lo = mid;
    }
    const pLcl = Math.exp(0.5 * (lo + hi));
    let temperature = dryT(pLcl), pressure = pLcl;
    let top = -1;
    for (let k = bottom; k >= 0; k--) {
      if (p[k] >= pLcl) {
        Tref[k] = dryT(p[k]);
      } else {
        const steps = 4, dlnp = (Math.log(p[k]) - Math.log(pressure)) / steps;
        for (let n = 0; n < steps; n++) {
          const k1 = moistLapse(temperature, pressure);
          const k2 = moistLapse(temperature + 0.5 * dlnp * k1, pressure * Math.exp(0.5 * dlnp));
          temperature += dlnp * k2;
          pressure *= Math.exp(dlnp);
        }
        Tref[k] = temperature;
        if (Tref[k] > T[k]) top = k;
      }
    }
    if (top < 0) return -1;
    for (let k = top; k <= bottom; k++) qref[k] = referenceHumidity * saturationHumidity(Tref[k], p[k]);
    return top;
  }

  function convectColumn(i, pi, theta, q, dt) {
    const bottom = K - 1;
    const top = referenceProfile(i, pi, theta, q);
    if (top < 0 || top === bottom) return 0;
    let heating = 0, drying = 0, depth = 0;
    for (let k = top; k <= bottom; k++) {
      heating += cp * (Tref[k] - T[k]) * dp[k];
      drying -= (qref[k] - q[k * C + i]) * dp[k];
      depth += dp[k];
    }
    if (heating <= 0) return 0;
    const rate = dt / relaxationTime;
    let rain = 0;
    if (drying > 0) {
      const shift = (latentHeat * drying - heating) / (cp * depth);
      for (let k = top; k <= bottom; k++) Tref[k] += shift;
      rain = drying / g * rate;
    } else {
      const shiftQ = drying / depth, shiftT = -heating / (cp * depth);
      for (let k = top; k <= bottom; k++) { qref[k] += shiftQ; Tref[k] += shiftT; }
    }
    for (let k = top; k <= bottom; k++) {
      const idx = k * C + i;
      theta[idx] += (Tref[k] - T[k]) * rate / exnerLayer[idx];
      q[idx] += (qref[k] - q[idx]) * rate;
    }
    return rain;
  }

  function fillColumn(i, pi, q) {
    for (let k = 0; k < K - 1; k++) {
      const idx = k * C + i;
      if (q[idx] < 0) {
        q[idx + C] += q[idx] * dSigma[k] / dSigma[k + 1];
        q[idx] = 0;
      }
    }
    const bottom = (K - 1) * C + i;
    if (q[bottom] < 0) {
      budget.lost -= mesh.areaCell[i] * pi[i] * dSigma[K - 1] / g * q[bottom];
      q[bottom] = 0;
    }
  }

  function adjust(state, iFrom, iTo, dt) {
    const [pi, theta, , , q] = state;
    for (let i = iFrom; i < iTo; i++) {
      core.diagnoseColumn(i, pi, theta, q);
      const condensed = condenseColumn(i, pi, theta, q);
      const convected = convectColumn(i, pi, theta, q, dt);
      fillColumn(i, pi, q);
      precipitation[i] += condensed + convected;
      budget.condensation += mesh.areaCell[i] * condensed;
      budget.convection += mesh.areaCell[i] * convected;
    }
  }

  function columnWater(pi, q, i) {
    let water = 0;
    for (let k = 0; k < K; k++) water += pi[i] * dSigma[k] / g * q[k * C + i];
    return water;
  }

  return { adjust, condenseColumn, convectColumn, fillColumn, referenceProfile, columnWater, precipitation, budget, latentHeat, shared: { precipitation: precipBuffer }, reference: { T: Tref, q: qref } };
}
