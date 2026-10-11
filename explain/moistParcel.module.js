/*
 * A parcel of the model's moist air: Bolton's saturation humidity and
 * lifting condensation level from js/physics/moist.module.js, a fixed latent
 * heat L and the dry R, c_p and g everywhere, as in the model. The ascent is
 * a liquid-only pseudo-adiabat: water condenses at saturation and leaves at
 * once. Heights are explain/physics.module.js's standard atmosphere, Part 1's.
 * moistAdjustColumn is explainer code, not the model's convection: it puts
 * the lowest layer and the unstable layers above it on the adiabat of the
 * lowest layer's air, conserving Σ c_p T Δm plus the heat it is given.
 * Pressures are Pa, temperatures K, humidities and condensate kg/kg, masses
 * kg/m²; a column's layers run top first, as in radiativeColumn.module.js.
 */
import { saturationHumidity, liftingCondensationLevel, LATENT_HEAT, R_VAPOR, EPSILON } from '../js/physics/moist.module.js';
import { heightOf, R, CP, G, KAPPA, P0 } from './physics.module.js';

export const DRY_LAPSE = G / CP;

export function moistLapse(T, p) {
  const qs = saturationHumidity(T, p);
  return G / CP * (1 + LATENT_HEAT * qs / (R * T)) / (1 + LATENT_HEAT * LATENT_HEAT * qs / (CP * R_VAPOR * T * T));
}

export function adiabatSlope(T, p, saturated) {
  return saturated ? moistLapse(T, p) * R * T / (G * p) : KAPPA * T / p;
}

export function adiabatThrough(T0, q0, p0, pressures, out = new Float64Array(pressures.length), moist = q0 > 0) {
  let T = T0, p = p0, saturated = moist && saturationHumidity(T, p) <= q0;
  for (let n = 0; n < pressures.length; n++) {
    const dp = (pressures[n] - p) / 8;
    for (let s = 0; s < 8; s++) {
      if (moist && !saturated && saturationHumidity(T, p) <= q0) saturated = true;
      const a = adiabatSlope(T, p, saturated), c = adiabatSlope(T + a * dp, p + dp, saturated);
      T += 0.5 * (a + c) * dp;
      p += dp;
    }
    p = pressures[n];
    out[n] = T;
  }
  return out;
}

export function dewPoint(q, p) {
  if (!(q > 0)) return NaN;
  const e = q * p / (EPSILON + (1 - EPSILON) * q), y = Math.log(e / 611.2);
  return (273.15 * 17.67 - 29.65 * y) / (17.67 - y);
}

export function liftParcel({ T0 = 288.15, rh = 0.7, p0 = P0 } = {}, pTop = 100e2, { step = 10, every = 100 } = {}) {
  const q0 = rh * saturationHumidity(T0, p0), found = liftingCondensationLevel(T0, q0, p0, KAPPA);
  const lcl = found && { pressure: found.pressure, temperature: found.temperature, height: heightOf(found.pressure) };
  const n = Math.floor((p0 - pTop) / every + 1e-9) + 1, stride = Math.round(every / step);
  const path = { p: new Float64Array(n), z: new Float64Array(n), T: new Float64Array(n), q: new Float64Array(n), condensate: new Float64Array(n), theta: new Float64Array(n), dry: new Float64Array(n), dew: new Float64Array(n) };
  let T = T0, p = p0, saturated = false;
  const record = (j) => {
    const q = saturated ? Math.min(q0, saturationHumidity(T, p)) : q0;
    path.p[j] = p;
    path.z[j] = heightOf(p);
    path.T[j] = T;
    path.q[j] = q;
    path.condensate[j] = q0 - q;
    path.theta[j] = T * (P0 / p) ** KAPPA;
    path.dry[j] = T0 * (p / p0) ** KAPPA;
    path.dew[j] = saturated ? T : dewPoint(q0, p);
  };
  record(0);
  for (let j = 1; j < n; j++) for (let s = 0; s < stride; s++) {
    if (!saturated && lcl && p <= lcl.pressure) saturated = true;
    const a = adiabatSlope(T, p, saturated), b = adiabatSlope(T - a * step, p - step, saturated);
    T -= 0.5 * (a + b) * step;
    p -= step;
    if (s === stride - 1) record(j);
  }
  return { T0, rh, p0, q0, dew0: dewPoint(q0, p0), lcl, ...path };
}

export function parcelAt(parcel, height) {
  const { z } = parcel, n = z.length;
  if (height <= z[0]) return pick(parcel, 0, 0);
  if (height >= z[n - 1]) return pick(parcel, n - 1, 0);
  let j = 1;
  while (z[j] < height) j++;
  return pick(parcel, j - 1, (height - z[j - 1]) / (z[j] - z[j - 1]));
}

function pick(parcel, j, f) {
  const out = {};
  for (const key of ['p', 'z', 'T', 'q', 'condensate', 'theta', 'dry', 'dew']) out[key] = f > 0 ? parcel[key][j] + f * (parcel[key][j + 1] - parcel[key][j]) : parcel[key][j];
  out.cloudy = parcel.lcl !== null && out.p <= parcel.lcl.pressure;
  return out;
}

export function moistAdjustColumn(T, pressure, mass, { heat = 0, rh = null, q = null, moist = true, top = T.length - 1, out = new Float64Array(T.length) } = {}) {
  const b = T.length - 1, up = Float64Array.from({ length: b }, (_, n) => pressure[b - 1 - n]), path = new Float64Array(b);
  const humidity = (Tb) => (rh !== null ? rh * saturationHumidity(Tb, pressure[b]) : q ?? 0);
  const excess = (Tb, t) => {
    const n = b - Math.max(0, t - 1);
    adiabatThrough(Tb, humidity(Tb), pressure[b], up.subarray(0, n), path, moist);
    out[b] = Tb;
    for (let j = 0; j < n; j++) out[b - 1 - j] = path[j];
    let sum = -heat;
    for (let k = t; k <= b; k++) sum += CP * mass[k] * (out[k] - T[k]);
    return sum;
  };
  const solve = (t) => {
    let x0 = T[b], f0 = excess(x0, t), x1 = T[b] + 1, f1 = excess(x1, t);
    for (let n = 0; n < 20 && Math.abs(f1) > 1e-3 && f1 !== f0; n++) {
      const x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
      x0 = x1; f0 = f1; x1 = x2; f1 = excess(x1, t);
    }
  };
  let t = Math.min(b, Math.max(0, top + 1));
  solve(t);
  while (t > 0 && T[t - 1] < out[t - 1]) solve(--t);
  for (let k = t; k <= b; k++) T[k] = out[k];
  return t;
}
