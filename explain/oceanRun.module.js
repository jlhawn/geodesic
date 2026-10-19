/* Band values in ocean_N16.bin are area means over the cells whose centers fall in each 5° bin, so every
   integral here weights them by the file's binArea, never by sin-band areas. Each transport is integrated
   northward from the South Pole after its own global mean is spread evenly over the planet, as Part 6 did.
   Sections and overturning are indexed by the hybrid ocean's layers, bin * layers + layer: layers 0 to
   fixedLayers - 1 are always fixed-depth, and a layer's mean thickness over its minimum gives its mode weight. */
import { calendarDate, season, peakLatitude, recordedYear, PART6_TRANSPORT } from './wetRun.module.js';
export { calendarDate, season, peakLatitude, recordedYear, PART6_TRANSPORT };

export const DATA_URL = new URL('./data/ocean_N16.bin', import.meta.url).href;
export const PART7_TRANSPORT = {
  air: [0, -0.09, -0.361, -0.809, -1.407, -2.12, -2.878, -3.637, -4.31, -4.791, -5.038, -5.027, -4.824, -4.474, -3.974, -3.306, -2.487, -1.556, -0.554, 0.441, 1.309, 2.083, 2.861, 3.561, 4.006, 4.221, 4.249, 4.115, 3.842, 3.381, 2.763, 2.072, 1.383, 0.801, 0.366, 0.094, 0],
  latent: [0, -0.008, -0.031, -0.076, -0.153, -0.274, -0.477, -0.796, -1.211, -1.643, -2.025, -2.117, -1.633, -0.7, 0.218, 0.567, 0.633, 0.787, 1.306, 1.634, 1.479, 1.234, 0.95, 0.979, 1.569, 2.064, 2.052, 1.721, 1.334, 0.983, 0.684, 0.432, 0.252, 0.137, 0.057, 0.014, 0],
};
export const GROUPS = ['the twelve fixed-depth surface layers', 'layers lighter than 1026', 'layers 1026 to 1026.8', 'layers from 1026.8'];
export const MODE_DEPTH = 5;
const DAY = 86400, YEAR_DAYS = 365;

const clamp = (v, lo, hi) => Math.min(hi, Math.max(lo, v));
export const groupOf = (k, label, fixed) => (k < fixed ? 0 : label < 1026 ? 1 : label < 1026.8 ? 2 : 3);
export const modeWeight = (thickness, minimum) => clamp((thickness - minimum) / MODE_DEPTH, 0, 1);
export const modeName = (weight) => (weight < 0.02 ? 'fixed-depth' : weight > 0.98 ? 'isopycnal' : 'between');
export const globalMean = (values, area) => values.reduce((sum, v, b) => sum + v * area[b], 0) / area.reduce((sum, a) => sum + a, 0);
const edgesOf = (bins) => Array.from({ length: bins + 1 }, (_, b) => -90 + b * 180 / bins);
const centersOf = (bins) => Array.from({ length: bins }, (_, b) => -90 + (b + 0.5) * 180 / bins);
const pick = (block, month) => (month ? block.months.find((m) => m.name === month) : block);

export function integrate(net, area) {
  const imbalance = globalMean(net, area), transport = [0];
  net.forEach((v, b) => transport.push(transport[b] + (v - imbalance) * area[b] / 1e15));
  return { net, imbalance, transport };
}

export const latentNet = ({ evaporation, precipitation }, vaporRate, latentHeat) => evaporation.map((e, b) => latentHeat * ((e - precipitation[b]) / DAY - vaporRate[b]));

const along = (values, lat) => {
  const n = values.length - 1, x = Math.min(n, Math.max(0, (lat + 90) / (180 / n))), b = Math.min(n - 1, Math.floor(x)), t = x - b;
  return values[b] + t * (values[b + 1] - values[b]);
};

const signChanges = (values, at) => {
  const out = [];
  for (let n = 0; n + 1 < values.length; n++) if (values[n] * values[n + 1] < 0) out.push(at(n) + (at(n + 1) - at(n)) * values[n] / (values[n] - values[n + 1]));
  return out;
};

export function heatBooks(header) {
  const { heatBooks: books, zonalBudget: budget, books: { record } } = header, area = books.binArea;
  const first = books.boundaries[0], last = books.boundaries.at(-1), seconds = last.time - first.time;
  const stores = Object.fromEntries(Object.keys(first.bands).map((name) => [name, last.bands[name].map((v, b) => (v - first.bands[name][b]) / seconds)]));
  const sea = stores.ocean.map((v, b) => v + stores.ice[b]), air = stores.airEnthalpy.map((v, b) => v + stores.latent[b]), total = sea.map((v, b) => v + air[b]);
  const vapor = stores.vapor.map((v, b) => v + stores.condensate[b]), latent = vapor.map((v) => header.constants.latentHeat * v);
  const toaNet = budget.absorbedSolar.map((v, b) => v - budget.outgoingLongwave[b]), seaInput = budget.netSurfaceHeat.map((v, b) => v - budget.snowFusion[b]);
  const mean = (values) => globalMean(values, area), g = Object.fromEntries(Object.entries(stores).map(([name, v]) => [name, mean(v)]));
  const global = {
    ...g, sea: mean(sea), air: mean(air), total: mean(total), toaNet: mean(toaNet), surfaceNet: mean(budget.netSurfaceHeat), snowFusion: mean(budget.snowFusion), seaInput: mean(seaInput),
    residual: mean(total) - mean(toaNet), surfaceResidual: mean(sea) - mean(seaInput), atmosphereResidual: mean(air) - mean(toaNet) + mean(seaInput),
  };
  const negative = (values) => integrate(values.map((v) => -v), area).transport;
  return {
    seconds, area, stores, sea, air, total, vapor, latent, toaNet, seaInput, global, record,
    corrections: { total: negative(total), sea: negative(sea), air: negative(air), latent: negative(latent) },
  };
}

export function layerGroups(header, month = null) {
  const flow = pick(header.overturning, month), L = header.layers, caps = header.overturning.caps;
  const groups = GROUPS.map((name) => ({ name, mass: caps.map(() => 0), eddyMass: caps.map(() => 0), heat: caps.map(() => 0), eddyHeat: caps.map(() => 0) }));
  caps.forEach((_, c) => {
    for (let k = 0; k < L; k++) {
      const g = groups[groupOf(k, header.classes[k], header.fixedLayers)], i = c * L + k;
      g.mass[c] += flow.mass[i]; g.eddyMass[c] += flow.eddyMass[i]; g.heat[c] += flow.heat[i]; g.eddyHeat[c] += flow.eddyHeat[i];
    }
  });
  return { caps, groups };
}

export function cellHeat(header) {
  const ov = header.overturning, sec = header.section, L = header.layers, rc = header.constants.rhoCp;
  return ov.caps.map((_, c) => {
    let sum = 0;
    for (let k = 0; k < L; k++) {
      const thick = (b) => sec.base[b * L + k] - (k ? sec.base[b * L + k - 1] : 0), T0 = sec.T[c * L + k], T1 = sec.T[(c + 1) * L + k];
      const T = T0 === null ? T1 ?? 0 : T1 === null ? T0 : (T0 * thick(c) + T1 * thick(c + 1)) / (thick(c) + thick(c + 1));
      sum += rc * (ov.mass[c * L + k] + ov.eddyMass[c * L + k]) * 1e6 * T / 1e15;
    }
    return sum;
  });
}

export function transportSplit(header) {
  const budget = header.zonalBudget, bins = budget.bins, area = budget.binArea, kept = heatBooks(header), edges = edgesOf(bins), centers = centersOf(bins);
  const total = integrate(kept.toaNet.map((v, b) => v - kept.total[b]), area);
  const sea = integrate(kept.seaInput.map((v, b) => v - kept.sea[b]), area);
  const air = { net: total.net.map((v, b) => v - sea.net[b]), imbalance: total.imbalance - sea.imbalance, transport: total.transport.map((v, b) => v - sea.transport[b]) };
  const latent = integrate(latentNet(budget, kept.vapor, header.constants.latentHeat), area);
  const dry = { net: air.net.map((v, b) => v - latent.net[b]), imbalance: air.imbalance - latent.imbalance, transport: air.transport.map((v, b) => v - latent.transport[b]) };
  const bare = { total: integrate(kept.toaNet, area).transport, sea: integrate(kept.seaInput, area).transport };
  const ov = header.overturning, pole = (values) => [0, ...values, 0], { groups } = layerGroups(header), perKelvin = header.constants.rhoCp * 1e6 / 1e15;
  const direct = {
    transport: pole(ov.caps.map((_, c) => groups.reduce((sum, g) => sum + g.heat[c] + g.eddyHeat[c], 0))),
    mass: pole(ov.caps.map((_, c) => groups.reduce((sum, g) => sum + g.mass[c] + g.eddyMass[c], 0))),
    cells: pole(cellHeat(header)), imbalance: pole(ov.imbalance), etaNet: pole(ov.etaNet), perKelvin: pole(ov.imbalance.map((v) => v * perKelvin)),
    groups: groups.map((g) => ({ name: g.name, mass: pole(g.mass.map((v, c) => v + g.eddyMass[c])), heat: pole(g.heat.map((v, c) => v + g.eddyHeat[c])) })),
  };
  const part7 = { air: PART7_TRANSPORT.air, latent: PART7_TRANSPORT.latent, dry: PART7_TRANSPORT.air.map((v, b) => v - PART7_TRANSPORT.latent[b]) };
  const series = { total: total.transport, sea: sea.transport, air: air.transport, latent: latent.transport, dry: dry.transport, direct: direct.transport, cells: direct.cells, bareTotal: bare.total, bareSea: bare.sea, part7: part7.air, part7Latent: part7.latent, part7Dry: part7.dry, part6: PART6_TRANSPORT };
  const zero = (values) => signChanges(values.slice(1, -1), (n) => edges[n + 1]);
  return {
    edges, centers, area, books: kept, total, sea, air, latent, dry, bare, direct, part7, part6: PART6_TRANSPORT,
    zeros: { total: zero(total.transport), sea: zero(sea.transport), air: zero(air.transport), latent: zero(latent.transport), dry: zero(dry.transport) },
    peak(which, sign, from = -90, to = 90) {
      const values = series[which];
      let best = null;
      values.forEach((v, e) => { if (edges[e] >= from && edges[e] <= to && (best === null || sign * v > sign * values[best])) best = e; });
      return { latitude: edges[best], transport: values[best] };
    },
    at(lat) {
      const out = Object.fromEntries(Object.entries(series).map(([name, values]) => [name, along(values, lat)]));
      return { ...out, groups: direct.groups.map((g) => ({ name: g.name, mass: along(g.mass, lat), heat: along(g.heat, lat) })) };
    },
  };
}

export function waterCycle(header) {
  const budget = header.zonalBudget, area = budget.binArea, kept = heatBooks(header), mean = (values) => globalMean(values, area);
  const gap = budget.evaporation.map((v, b) => v - budget.precipitation[b]), centers = centersOf(budget.bins);
  const vapor = integrate(budget.evaporation.map((v, b) => gap[b] / DAY - kept.vapor[b]), area);
  return {
    centers, area, gap, evaporation: budget.evaporation, precipitation: budget.precipitation,
    global: Object.fromEntries(['evaporation', 'precipitation', 'convectiveRain', 'largeScaleRain', 'snowfall', 'latentHeat', 'sensibleHeat', 'netSurfaceHeat', 'surfaceShortwave'].map((name) => [name, mean(budget[name])])),
    gapCrossings: signChanges(gap, (b) => centers[b]),
    vaporSverdrups: vapor.transport.map((v) => v * 1e15 / 1000 / 1e6), waterImbalance: vapor.imbalance * DAY,
    salt: kept.stores.salt,
  };
}

export function yearlyMeans({ daily }, name) {
  const values = daily[name], out = [];
  for (let n = 0; n + YEAR_DAYS <= values.length; n += YEAR_DAYS) {
    let sum = 0, count = 0;
    for (let d = n; d < n + YEAR_DAYS; d++) if (values[d] !== null) { sum += values[d]; count++; }
    out.push(count ? sum / count : null);
  }
  return out;
}

export function zonalMeans(header, month = null) {
  const block = pick(header.zonalMean, month), bins = header.zonalMean.bins;
  return {
    ...block, name: month ?? 'year', centers: centersOf(bins), area: header.zonalMean.binArea };
}

export function layerSection(header, month = null) {
  const block = pick(header.section, month), L = header.layers, bins = header.section.bins, centers = centersOf(bins);
  const columns = centers.map((_, b) => {
    const out = [];
    for (let k = 0; k < L; k++) {
      const i = b * L + k, top = k ? block.base[i - 1] : 0, thickness = block.base[i] - top, weight = k < header.fixedLayers ? 0 : modeWeight(thickness, header.minimumThickness[k]);
      if (block.T[i] !== null) out.push({ layer: k, label: header.classes[k], top, base: block.base[i], thickness, weight, mode: modeName(weight), T: block.T[i], S: block.S[i], u: block.u[i], v: block.v ? block.v[i] : null });
    }
    return out;
  });
  return { name: month ?? 'year', centers, labels: header.classes, layers: L, base: block.base, T: block.T, S: block.S, u: block.u, v: block.v ?? null, columns, samples: block.samples };
}

export function layerOverturning(header, month = null, { closed = true } = {}) {
  const flow = pick(header.overturning, month), base = pick(header.section, month).base, L = header.layers, caps = header.overturning.caps, M = L + 1;
  const psi = new Float64Array(caps.length * M), depth = new Float64Array(caps.length * M), bottom = new Float64Array(caps.length);
  caps.forEach((_, c) => {
    const thick = [], at = (b, k) => (k < 0 ? 0 : base[b * L + k]);
    for (let k = 0; k < L; k++) thick.push(0.5 * (at(c, k) - at(c, k - 1) + at(c + 1, k) - at(c + 1, k - 1)));
    const column = thick.reduce((sum, v) => sum + v, 0), fix = closed ? flow.imbalance[c] / column : 0;
    for (let k = 0; k < L; k++) {
      const i = c * L + k;
      psi[c * M + k + 1] = psi[c * M + k] + flow.mass[i] + flow.eddyMass[i] - fix * thick[k];
      depth[c * M + k + 1] = depth[c * M + k] + thick[k];
    }
    bottom[c] = psi[c * M + L];
  });
  return {
    name: month ?? 'year', caps, labels: header.classes, layers: L, psi, depth, bottom, imbalance: flow.imbalance, etaNet: flow.etaNet, net: flow.net,
    cell(sign, from = -90, to = 90) {
      let best = { psi: 0 };
      caps.forEach((lat, c) => { if (lat >= from && lat <= to) for (let k = 0; k <= L; k++) if (sign * psi[c * M + k] > sign * best.psi) best = { psi: psi[c * M + k], latitude: lat, boundary: k, depth: depth[c * M + k] }; });
      return best;
    },
  };
}
