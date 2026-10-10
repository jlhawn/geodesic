import { rampRGB, sequentialRGB, DARK_NEUTRAL } from './runtime.module.js';

const R = 287, g = 9.80616, KAPPA = 287 / 1004.5, a = 6.37122e6, omega = 7.292e-5;
const DARK = 'rgb(52, 55, 62)';

export function jw06Point(lat, p, u0 = 35) {
  const eta0 = 0.252, T0 = 288, gamma = 0.005, deltaT = 4.8e5, etaT = 0.2;
  const eta = p / 1000, ev = (eta - eta0) * Math.PI / 2, sn = Math.sin(lat), cs = Math.cos(lat);
  const A = -2 * sn ** 6 * (cs * cs + 1 / 3) + 10 / 63, B = 8 / 5 * cs ** 3 * (sn * sn + 2 / 3) - Math.PI / 4;
  let mean = T0 * eta ** (R * gamma / g);
  if (eta < etaT) mean += deltaT * (etaT - eta) ** 5;
  const T = mean + 0.75 * (eta * Math.PI * u0 / R) * Math.sin(ev) * Math.sqrt(Math.cos(ev)) * (A * 2 * u0 * Math.cos(ev) ** 1.5 + B * a * omega);
  return { T, theta: T * (1 / eta) ** KAPPA, u: u0 * Math.cos(ev) ** 1.5 * Math.sin(2 * lat) ** 2 };
}

export function jw06Section(u0 = 35) {
  const lats = Array.from({ length: 73 }, (_, k) => -90 + 2.5 * k), pressures = Array.from({ length: 48 }, (_, k) => 50 + 20 * k);
  const temperature = [], theta = [], wind = [];
  for (const latDeg of lats) {
    const T = [], th = [], u = [];
    for (const p of pressures) { const point = jw06Point(latDeg * Math.PI / 180, p, u0); T.push(point.T); th.push(point.theta); u.push(point.u); }
    temperature.push(T); theta.push(th); wind.push(u);
  }
  return { lats, pressures, temperature, theta, wind };
}

export function heldSuarezEquilibrium() {
  const lats = Array.from({ length: 73 }, (_, k) => -90 + 2.5 * k), pressures = Array.from({ length: 48 }, (_, k) => 50 + 20 * k);
  const temperature = [], theta = [];
  for (const latDeg of lats) {
    const s = Math.sin(latDeg * Math.PI / 180), c2 = 1 - s * s, T = [], th = [];
    for (const p of pressures) { const t = Math.max(200, (315 - 60 * s * s - 10 * Math.log(p / 1000) * c2) * (p / 1000) ** KAPPA); T.push(t); th.push(t * (1000 / p) ** KAPPA); }
    temperature.push(T); theta.push(th);
  }
  return { lats, pressures, temperature, theta };
}

export function zonalMeanSection(zonalMean) {
  const { bins, sigma, u, T } = zonalMean, K = sigma.length;
  const lats = Array.from({ length: bins }, (_, b) => -90 + (b + 0.5) * 180 / bins), pressures = sigma.map((s) => s * 1000);
  const temperature = [], theta = [], wind = [];
  for (let b = 0; b < bins; b++) {
    const Tb = [], th = [], ub = [];
    for (let k = 0; k < K; k++) { const t = T[b * K + k]; Tb.push(t); th.push(t * (1 / sigma[k]) ** KAPPA); ub.push(u[b * K + k]); }
    temperature.push(Tb); theta.push(th); wind.push(ub);
  }
  return { lats, pressures, temperature, theta, wind };
}

export const LAYER_FIELDS = [
  { name: 'surfacePressure', label: 'surface pressure', unit: 'hPa', digits: 0, color: (v, rgb) => rampRGB(0.5 + (v - 1000) / 60, rgb, DARK_NEUTRAL), legend: [['ramp', 'surface pressure from 970 to 1030 hPa', 'cool', 'warm', DARK]] },
  { name: 'temperature850', label: 'temperature 1.5 km up', unit: 'K', digits: 0, color: (v, rgb) => rampRGB((v - 240) / 60, rgb, DARK_NEUTRAL), legend: [['ramp', 'temperature at 850 hPa, about 1.5 km up, from 240 to 300 K', 'cool', 'warm', DARK]] },
  { name: 'wind250', label: 'wind speed 10 km up', unit: 'm/s', digits: 0, color: (v, rgb) => sequentialRGB(v / 60, rgb), legend: [['ramp', 'wind speed at 250 hPa, about 10 km up, from 0 to 60 m/s', DARK, 'warm']] },
];
