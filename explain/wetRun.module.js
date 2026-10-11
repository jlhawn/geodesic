import { heatTransport } from './sunlitRun.module.js';
export { calendarDate, season } from './sunlitRun.module.js';

export const DATA_URL = new URL('./data/water_N16.bin', import.meta.url).href;
export const PART6_TRANSPORT = [0, -0.021, -0.087, -0.197, -0.349, -0.532, -0.711, -0.868, -0.985, -1.053, -1.08, -1.059, -0.991, -0.896, -0.792, -0.679, -0.498, -0.254, 0.02, 0.294, 0.539, 0.721, 0.839, 0.942, 1.03, 1.088, 1.101, 1.062, 0.986, 0.86, 0.694, 0.51, 0.326, 0.178, 0.076, 0.018, 0];
const DAY = 86400, D = Math.PI / 180, YEAR_DAYS = 365;

export function storage({ storage: { start, end }, constants: { latentHeat } }) {
  const seconds = end.time - start.time, rate = (name, scale = 1) => end[name].map((v, b) => scale * (v - start[name][b]) / seconds);
  const surface = rate('surfaceEnergy'), air = rate('airEnergy'), vapor = rate('precipitableWater', latentHeat), cloud = rate('cloudCondensate', latentHeat);
  return { seconds, surface, air, vapor, cloud, total: surface.map((v, b) => v + air[b] + vapor[b]) };
}

export function recordedYear({ daily }) {
  const from = Math.round(-daily.from / DAY), mean = (name) => daily[name].slice(from, from + YEAR_DAYS).reduce((sum, v) => sum + v, 0) / YEAR_DAYS;
  const year = Object.fromEntries(Object.keys(daily).filter((name) => Array.isArray(daily[name])).map((name) => [name, mean(name)]));
  return {
    ...year,
    balance: year.absorbedSolar - year.outgoingLongwave,
    cloudSunlight: year.absorbedSolar - year.clearAbsorbedSolar,
    cloudInfrared: year.clearOutgoingLongwave - year.outgoingLongwave,
    residenceDays: year.precipitableWater / year.precipitation,
  };
}

export function peakLatitude(values, from = 0, to = values.length, sign = 1) {
  const width = 180 / values.length;
  let b = from;
  for (let c = from; c < to; c++) if (sign * values[c] > sign * values[b]) b = c;
  const left = values[b - 1] ?? values[b], mid = values[b], right = values[b + 1] ?? values[b], curve = left - 2 * mid + right, shift = curve ? 0.5 * (left - right) / curve : 0;
  return { latitude: -90 + (b + 0.5 + shift) * width, value: mid - 0.25 * (left - right) * shift };
}

const signChanges = (values, at) => {
  const out = [];
  for (let n = 0; n + 1 < values.length; n++) if (values[n] * values[n + 1] < 0) out.push(at(n) + (at(n + 1) - at(n)) * values[n] / (values[n] - values[n + 1]));
  return out;
};

export function waterTransport(header) {
  const { zonalBudget: budget, constants: { latentHeat, radius } } = header, kept = storage(header), bins = budget.bins;
  const total = heatTransport(budget, kept.total, radius);
  const latent = heatTransport({ bins, absorbedSolar: budget.evaporation.map((v) => latentHeat * v / DAY), outgoingLongwave: budget.precipitation.map((v) => latentHeat * v / DAY) }, kept.vapor.map((v, b) => v + kept.cloud[b]), radius);
  const dry = total.transport.map((v, b) => v - latent.transport[b]);
  const edges = total.edges, gap = budget.evaporation.map((v, b) => v - budget.precipitation[b]);
  const along = (values, lat) => { const x = Math.min(bins, Math.max(0, (lat + 90) / (180 / bins))), b = Math.min(bins - 1, Math.floor(x)), t = x - b; return values[b] + t * (values[b + 1] - values[b]); };
  const extreme = (values, sign, from = -90, to = 90) => {
    let best = null;
    values.forEach((v, b) => { if (edges[b] >= from && edges[b] <= to && (best === null || sign * v > sign * values[best])) best = b; });
    return { latitude: edges[best], transport: values[best] };
  };
  return {
    storage: kept, total, latent, dry, edges, centers: total.centers, evaporation: budget.evaporation, precipitation: budget.precipitation, gap,
    imbalance: total.imbalance, waterImbalance: latent.imbalance * DAY / latentHeat,
    gapCrossings: signChanges(gap, (b) => total.centers[b]),
    zeros: { total: signChanges(total.transport, (b) => edges[b]), latent: signChanges(latent.transport, (b) => edges[b]), dry: signChanges(dry, (b) => edges[b]) },
    peak: (which, sign, from, to) => extreme(which === 'total' ? total.transport : which === 'latent' ? latent.transport : which === 'dry' ? dry : PART6_TRANSPORT, sign, from, to),
    at(lat) {
      const water = latent.at(lat), toDays = DAY / latentHeat;
      return { evaporation: water.absorbed * toDays, precipitation: water.outgoing * toDays, total: total.at(lat).transport, latent: water.transport, dry: along(dry, lat), part6: along(PART6_TRANSPORT, lat) };
    },
  };
}

export function overturning({ zonalMean, levels, constants: { radius, g } }, month = null) {
  const zm = month ? zonalMean.months.find((m) => m.name === month) : zonalMean, bins = zonalMean.bins, K = zonalMean.sigma.length;
  const dSigma = levels.slice(1).map((v, k) => v - levels[k]), width = 180 / bins;
  const lats = Array.from({ length: bins }, (_, b) => -90 + (b + 0.5) * width);
  const v = new Float64Array(bins * K), psi = new Float64Array(bins * (K + 1)), omega = new Float64Array(bins * (K + 1));
  for (let b = 0; b < bins; b++) {
    let mean = 0;
    for (let k = 0; k < K; k++) mean += zm.v[b * K + k] * dSigma[k];
    const ps = zm.surfacePressure[b] * 100, scale = 2 * Math.PI * radius * Math.cos(lats[b] * D) / g;
    for (let k = 0; k < K; k++) {
      v[b * K + k] = zm.v[b * K + k] - mean;
      psi[b * (K + 1) + k + 1] = psi[b * (K + 1) + k] + scale * v[b * K + k] * ps * dSigma[k];
    }
  }
  for (let b = 0; b < bins; b++) {
    const lo = Math.max(0, b - 1), hi = Math.min(bins - 1, b + 1), dPhi = (lats[hi] - lats[lo]) * D;
    for (let k = 0; k <= K; k++) omega[b * (K + 1) + k] = -g / (2 * Math.PI * radius * radius * Math.cos(lats[b] * D)) * (psi[hi * (K + 1) + k] - psi[lo * (K + 1) + k]) / dPhi;
  }
  const cell = (sign, from, to) => {
    let best = { psi: 0 };
    for (let b = 0; b < bins; b++) if (lats[b] >= from && lats[b] <= to) for (let k = 0; k <= K; k++) { const value = psi[b * (K + 1) + k]; if (sign * value > sign * best.psi) best = { psi: value, latitude: lats[b], sigma: levels[k] }; }
    return best;
  };
  return { lats, sigma: zonalMean.sigma, levels, K, v, psi, omega, surfacePressure: zm.surfacePressure, cell };
}
