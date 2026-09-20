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
 * Lifting condensation level of a parcel (T, q, p) by Bolton (1980):
 * the dew point from the vapour pressure, the LCL temperature from his
 * eq. 15, and the pressure along the dry adiabat. Returns null when the
 * parcel is already saturated or has no vapour.
 */
export function liftingCondensationLevel(T, q, p, kappa) {
  if (q <= 0) return null;
  const e = q * p / (EPSILON + (1 - EPSILON) * q);
  const y = Math.log(e / 611.2);
  const dewPoint = (273.15 * 17.67 - 29.65 * y) / (17.67 - y);
  if (dewPoint >= T) return { temperature: T, pressure: p };
  const temperature = 1 / (1 / (dewPoint - 56) + Math.log(T / dewPoint) / 800) + 56;
  return { temperature, pressure: p * Math.pow(temperature / T, 1 / kappa) };
}

/*
 * Moist physics for the sigma core, applied to the state after each
 * step: saturation adjustment between vapour and cloud condensate
 * (supersaturation condenses into cloud water, cloud water evaporates
 * into subsaturated air, latent heat to the layer), Kessler
 * autoconversion of cloud water above a threshold into rain that falls
 * out at once, the simplified Betts–Miller
 * convection of Frierson (2007) — a conditionally unstable column
 * relaxes over relaxationTime toward the moist adiabat of its lowest
 * layer and a fixed relative humidity, with the reference temperature
 * shifted so the column's enthalpy change equals the latent heat of
 * the condensate it produces, or with no rain when the column would
 * have to moisten — of which the fraction `detrainment` stays in the
 * column as cloud water spread through the anvil, the top `anvilDepth`
 * of pressure below the level of zero buoyancy, and the rest falls as
 * rain — and a filler that removes negative humidity by borrowing from
 * the layer below. Precipitation accumulates per cell (kg/m²); the
 * budget sums are area-weighted masses (kg).
 */
export function createMoistPhysics(mesh, core, {
  latentHeat = LATENT_HEAT, relaxationTime = 7200, referenceHumidity = 0.7,
  autoconversionThreshold = 2e-4, autoconversionRate = 1e-3, cloudLifetime = 3 * 3600,
  detrainment = 0.25, anvilDepth = 150e2, buffers = null,
} = {}) {
  const { K, C, dSigma, sigmaMid, cp, R, g, kappa, exnerLayer } = core.diagnostics;
  const precipBuffer = buffers && buffers.precipitation ? buffers.precipitation : new SharedArrayBuffer(8 * C);
  const precipitation = new Float64Array(precipBuffer);
  const T = new Float64Array(K), p = new Float64Array(K), dp = new Float64Array(K), Tref = new Float64Array(K), qref = new Float64Array(K);
  const budget = { condensation: 0, convection: 0, lost: 0 };

  function moistLapse(temperature, pressure) {
    const qs = saturationHumidity(temperature, pressure);
    return (R * temperature + latentHeat * qs) / (cp + latentHeat * latentHeat * EPSILON * qs / (R * temperature * temperature));
  }

  /*
   * Saturation adjustment: supersaturated vapour condenses into cloud
   * water and cloud water evaporates into subsaturated air, each with
   * one implicit step, so afterwards a layer is either saturated or
   * cloud-free. Returns the condensate formed (kg/m², negative when
   * cloud evaporated); nothing rains here.
   */
  function condenseColumn(i, pi, theta, q, qc) {
    let formed = 0;
    for (let k = 0; k < K; k++) {
      const idx = k * C + i;
      const ex = exnerLayer[idx];
      const temperature = theta[idx] * ex;
      const pressure = pi[i] * sigmaMid[k];
      const qs = saturationHumidity(temperature, pressure);
      const slope = qs * latentHeat / (R_VAPOR * temperature * temperature);
      let change = (q[idx] - qs) / (1 + latentHeat * slope / cp);
      if (change < 0) change = Math.max(change, -qc[idx]);
      if (change === 0) continue;
      q[idx] -= change;
      qc[idx] += change;
      theta[idx] += latentHeat * change / (cp * ex);
      formed += pi[i] * dSigma[k] / g * change;
    }
    return formed;
  }

  /*
   * Kessler autoconversion: cloud water above the threshold turns into
   * rain at autoconversionRate, and all cloud water decays over
   * cloudLifetime; the rain leaves the column at once.
   */
  function autoconvertColumn(i, pi, qc, dt) {
    let rain = 0;
    for (let k = 0; k < K; k++) {
      const idx = k * C + i;
      if (qc[idx] <= 0) continue;
      const excess = Math.max(0, qc[idx] - autoconversionThreshold);
      const converted = Math.min(qc[idx], excess * (1 - Math.exp(-autoconversionRate * dt)) + qc[idx] * (1 - Math.exp(-dt / cloudLifetime)));
      qc[idx] -= converted;
      rain += pi[i] * dSigma[k] / g * converted;
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
    const lcl = liftingCondensationLevel(Tp, qp, pb, kappa);
    if (!lcl || lcl.pressure < p[0]) return -1;
    let temperature = lcl.temperature, pressure = lcl.pressure;
    let top = -1;
    for (let k = bottom; k >= 0; k--) {
      if (p[k] >= lcl.pressure) {
        Tref[k] = dryT(p[k]);
      } else {
        const steps = 2, dlnp = (Math.log(p[k]) - Math.log(pressure)) / steps;
        for (let n = 0; n < steps; n++) {
          const k1 = moistLapse(temperature, pressure);
          const k2 = moistLapse(temperature + 0.5 * dlnp * k1, pressure * Math.exp(0.5 * dlnp));
          temperature += dlnp * k2;
          pressure *= Math.exp(dlnp);
        }
        Tref[k] = temperature;
        if (Tref[k] > T[k]) top = k;
        else if (T[k] - Tref[k] > 10) break;
      }
    }
    if (top < 0) return -1;
    for (let k = top; k <= bottom; k++) qref[k] = referenceHumidity * saturationHumidity(Tref[k], p[k]);
    return top;
  }

  function convectColumn(i, pi, theta, q, dt, qc = null) {
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
    if (qc === null || rain <= 0 || detrainment <= 0) return rain;
    let anvilMass = 0, anvilBottom = top;
    for (let k = top; k <= bottom && (k === top || anvilMass < anvilDepth); k++) { anvilMass += dp[k]; anvilBottom = k; }
    const detrained = detrainment * rain;
    for (let k = top; k <= anvilBottom; k++) qc[k * C + i] += detrained * g / anvilMass;
    return rain - detrained;
  }

  function fillColumn(i, pi, q) {
    if (!q) return;
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
    const [pi, theta, , , q, qc] = state;
    for (let i = iFrom; i < iTo; i++) {
      core.diagnoseColumn(i, pi, theta, q, qc);
      condenseColumn(i, pi, theta, q, qc);
      const convected = convectColumn(i, pi, theta, q, dt, qc);
      const rained = autoconvertColumn(i, pi, qc, dt);
      fillColumn(i, pi, q);
      fillColumn(i, pi, qc);
      precipitation[i] += rained + convected;
      budget.condensation += mesh.areaCell[i] * rained;
      budget.convection += mesh.areaCell[i] * convected;
    }
  }

  function columnWater(pi, q, i) {
    let water = 0;
    for (let k = 0; k < K; k++) water += pi[i] * dSigma[k] / g * q[k * C + i];
    return water;
  }

  return { adjust, condenseColumn, autoconvertColumn, convectColumn, fillColumn, referenceProfile, columnWater, precipitation, budget, latentHeat, shared: { precipitation: precipBuffer }, reference: { T: Tref, q: qref } };
}
