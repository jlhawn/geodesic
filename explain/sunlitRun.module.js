import { MARKERS } from './insolation.module.js';

export const DATA_URL = new URL('./data/radiation_N16.bin', import.meta.url).href;
export const RADIUS = 6.371e6;
export const AIR_HEAT_CAPACITY = 1004.5 * 1e5 / 9.80616;
const D = Math.PI / 180, YEAR_DAYS = 365, EQUINOX_DATE = 78;
const MONTHS = [['January', 31], ['February', 28], ['March', 31], ['April', 30], ['May', 31], ['June', 30], ['July', 31], ['August', 31], ['September', 30], ['October', 31], ['November', 30], ['December', 31]];
const SEASONS = [['spring', 'fall'], ['summer', 'winter'], ['fall', 'spring'], ['winter', 'summer']];

export function calendarDate(day) {
  let date = ((Math.floor(day + 0.5) + EQUINOX_DATE) % YEAR_DAYS + YEAR_DAYS) % YEAR_DAYS, month = 0;
  while (date >= MONTHS[month][1]) date -= MONTHS[month++][1];
  return { month: MONTHS[month][0], date: date + 1, hour: ((day * 24 + 12) % 24 + 24) % 24 };
}

export function season(day) {
  const d = ((day % YEAR_DAYS) + YEAR_DAYS) % YEAR_DAYS;
  let n = 0;
  while (n < 3 && d >= MARKERS[n + 1][0]) n++;
  return { north: SEASONS[n][0], south: SEASONS[n][1] };
}

export function storedHeat(run, cells) {
  const { header } = run, bins = header.zonalBudget.bins, times = run.times('year'), last = times.length - 1, seconds = times[last] - times[0];
  const change = (name) => {
    const k = header.fields.findIndex((f) => f.name === name), start = run.values('year', 0, k), end = run.values('year', last, k);
    const sum = new Float64Array(bins), area = new Float64Array(bins);
    cells.forEach((c, i) => { const b = Math.min(bins - 1, Math.floor((c.lat + 90) / 180 * bins)); sum[b] += (end[i] - start[i]) * c.area; area[b] += c.area; });
    return Array.from(sum, (v, b) => v / area[b]);
  };
  const ground = change('surfaceTemperature').map((dT) => header.ground.heatCapacity * dT / seconds);
  const air = change('temperature850').map((dT) => AIR_HEAT_CAPACITY * dT / seconds);
  return { ground, air, total: ground.map((v, b) => v + air[b]) };
}

export function heatTransport({ bins, absorbedSolar, outgoingLongwave }, stored = null, radius = RADIUS) {
  const width = 180 / bins, sphere = 4 * Math.PI * radius * radius, kept = stored ?? absorbedSolar.map(() => 0);
  const edges = Array.from({ length: bins + 1 }, (_, b) => -90 + b * width);
  const centers = Array.from({ length: bins }, (_, b) => -90 + (b + 0.5) * width);
  const area = centers.map((_, b) => 2 * Math.PI * radius * radius * (Math.sin(edges[b + 1] * D) - Math.sin(edges[b] * D)));
  const net = absorbedSolar.map((v, b) => v - outgoingLongwave[b] - kept[b]);
  const mean = (values) => values.reduce((sum, v, b) => sum + v * area[b], 0) / sphere;
  const imbalance = mean(net), transport = [0];
  for (let b = 0; b < bins; b++) transport.push(transport[b] + (net[b] - imbalance) * area[b] / 1e15);
  const crossings = [];
  for (let b = 0; b + 1 < bins; b++) if (net[b] * net[b + 1] < 0) crossings.push(centers[b] + width * net[b] / (net[b] - net[b + 1]));
  const binned = (values, lat) => {
    const x = Math.min(bins - 1, Math.max(0, (lat + 90) / width - 0.5)), b = Math.min(bins - 2, Math.floor(x)), t = x - b;
    return values[b] + t * (values[b + 1] - values[b]);
  };
  return {
    edges, centers, absorbed: absorbedSolar, outgoing: outgoingLongwave, stored: kept, net, transport, imbalance, crossings,
    absorbedMean: mean(absorbedSolar), outgoingMean: mean(outgoingLongwave), storedMean: mean(kept),
    at(lat) {
      const x = Math.min(bins, Math.max(0, (lat + 90) / width)), b = Math.min(bins - 1, Math.floor(x)), t = x - b;
      const absorbed = binned(absorbedSolar, lat), outgoing = binned(outgoingLongwave, lat), held = binned(kept, lat);
      return { absorbed, outgoing, stored: held, net: absorbed - outgoing - held, transport: transport[b] + t * (transport[b + 1] - transport[b]) };
    },
    peak(sign) {
      let best = 0;
      transport.forEach((v, b) => { if (sign * v > sign * transport[best]) best = b; });
      return { latitude: edges[best], transport: transport[best] };
    },
  };
}
