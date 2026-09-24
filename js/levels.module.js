import { DAY, YEAR } from './physics/radiation.module.js';
import { saturationHumidity } from './physics/moist.module.js';

export const LEVELS = ['surface', 1000, 850, 700, 500, 250, 70, 10];

/*
 * Wind, temperature and geopotential height on a pressure surface (hPa),
 * interpolated in ln p between layer midpoints; 'surface' is the lowest
 * layer. Below the lowest midpoint the wind and temperature are held and
 * the height is extrapolated hydrostatically. layerWind(k) returns the
 * cell-center wind vectors of layer k so that only the layers a level
 * touches are reconstructed. With q, relative humidity comes too.
 */
export function levelFields(core, pi, theta, layerWind, level, q = null) {
  const { K, C, sigmaMid, g, R, exnerLayer, geopotential } = core.diagnostics;
  const speed = new Float32Array(C), vector = new Float32Array(3 * C), temperature = new Float32Array(C), height = new Float32Array(C);
  const humidity = q ? new Float32Array(C) : null;
  const pressure = level === 'surface' ? Infinity : 100 * level;
  for (let i = 0; i < C; i++) {
    let k = 0;
    while (k < K - 2 && pi[i] * sigmaMid[k + 1] < pressure) k++;
    const pk = pi[i] * sigmaMid[k], pk1 = pi[i] * sigmaMid[k + 1];
    const t = pressure === Infinity ? 1 : (Math.log(pressure) - Math.log(pk)) / (Math.log(pk1) - Math.log(pk));
    const tw = Math.min(1, Math.max(0, t));
    const wk = layerWind(k), wk1 = layerWind(k + 1);
    const vx = wk[3 * i] + tw * (wk1[3 * i] - wk[3 * i]);
    const vy = wk[3 * i + 1] + tw * (wk1[3 * i + 1] - wk[3 * i + 1]);
    const vz = wk[3 * i + 2] + tw * (wk1[3 * i + 2] - wk[3 * i + 2]);
    vector[3 * i] = vx; vector[3 * i + 1] = vy; vector[3 * i + 2] = vz;
    speed[i] = Math.hypot(vx, vy, vz);
    const tk = theta[k * C + i] * exnerLayer[k * C + i], tk1 = theta[(k + 1) * C + i] * exnerLayer[(k + 1) * C + i];
    temperature[i] = tk + tw * (tk1 - tk);
    if (q) {
      const qk = q[k * C + i] + tw * (q[(k + 1) * C + i] - q[k * C + i]);
      const pressureHere = pressure === Infinity ? pi[i] * sigmaMid[K - 1] : Math.min(pressure, pi[i] * sigmaMid[K - 1]);
      humidity[i] = Math.min(1.5, qk / saturationHumidity(temperature[i], pressureHere));
    }
    if (t > 1) {
      height[i] = (geopotential[(K - 1) * C + i] - R * tk1 * Math.log(pressure === Infinity ? 1 : pressure / (pi[i] * sigmaMid[K - 1]))) / g;
    } else {
      height[i] = (geopotential[k * C + i] + t * (geopotential[(k + 1) * C + i] - geopotential[k * C + i])) / g;
    }
  }
  return { speed, vector, temperature, height, humidity };
}

/*
 * Comfort measures from temperature (°C), relative humidity (a fraction)
 * and wind (m/s). Dew point by the Magnus formula; wet-bulb by Stull's
 * fit; the misery index is the NWS heat index above 26.7 °C, the wind
 * chill below 10 °C, and the air temperature between.
 */
export function dewPoint(t, rh) {
  const r = Math.max(1e-3, Math.min(1, rh)), a = 17.625, b = 243.04, g = Math.log(r) + a * t / (b + t);
  return b * g / (a - g);
}
export function wetBulb(t, rh) {
  const p = 100 * Math.max(0.05, Math.min(0.99, rh));
  return t * Math.atan(0.151977 * Math.sqrt(p + 8.313659)) + Math.atan(t + p) - Math.atan(p - 1.676331) + 0.00391838 * Math.pow(p, 1.5) * Math.atan(0.023101 * p) - 4.686035;
}
export function heatIndex(t, rh) {
  const T = t * 9 / 5 + 32, R = 100 * Math.max(0, Math.min(1, rh));
  let hi = -42.379 + 2.04901523 * T + 10.14333127 * R - 0.22475541 * T * R - 6.83783e-3 * T * T - 5.481717e-2 * R * R + 1.22874e-3 * T * T * R + 8.5282e-4 * T * R * R - 1.99e-6 * T * T * R * R;
  if (R < 13 && T <= 112) hi -= ((13 - R) / 4) * Math.sqrt((17 - Math.abs(T - 95)) / 17);
  else if (R > 85 && T <= 87) hi += ((R - 85) / 10) * ((87 - T) / 5);
  return (hi - 32) * 5 / 9;
}
export function windChill(t, v) {
  const k = 3.6 * v;
  if (k < 4.8) return t;
  const p = Math.pow(k, 0.16);
  return 13.12 + 0.6215 * t - 11.37 * p + 0.3965 * t * p;
}
export function miseryIndex(t, rh, v) {
  if (t >= 26.7) return Math.max(t, heatIndex(t, rh));
  if (t <= 10) return Math.min(t, windChill(t, v));
  return t;
}
const SEASONS = [[0, 'spring equinox'], [0.25, 'summer solstice'], [0.5, 'autumn equinox'], [0.75, 'winter solstice'], [1, 'spring equinox']];

/*
 * The model's season as a phrase: the number of days to the nearest
 * equinox or solstice (northern names; t = 0 is the spring equinox).
 */
export function seasonPhrase(time) {
  const day = (time % YEAR) / DAY, year = YEAR / DAY;
  let best = null;
  for (const [fraction, name] of SEASONS) {
    const offset = day - fraction * year;
    if (!best || Math.abs(offset) < Math.abs(best.offset)) best = { offset, name };
  }
  const n = Math.round(Math.abs(best.offset));
  if (n === 0) return `the Northern ${best.name}`;
  return `${n} day${n === 1 ? '' : 's'} ${best.offset > 0 ? 'past' : 'before'} the Northern ${best.name}`;
}
