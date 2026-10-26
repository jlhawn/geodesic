// The near infrared's k-distributions (js/physics/nearInfraredTable.module.js) folded to POINTS (32)
// g-points for the layered shortwave (radiation option shortwave 'layered'):
//   node scripts/layeredFit.mjs            the fold and its clear-sky misses against the full table
//   WRITE=1 node scripts/layeredFit.mjs    also rewrites js/physics/layeredTable.module.js
//
// A band's gases multiply, so its transmission is a sum over every combination of one term per gas
// (the window term w0 included) of the combination's weight times exp(-m (depth of the terms
// chosen + the self continuum)); the table's bands have 5 to 25 combinations, 230 together. A g-point
// is a group of a band's combinations taken as one absorber whose depth in each layer is the group's
// weight-mean of its members' depths, times a scale below 1 that takes back the absorption that the
// mean of unlike depths adds (Jensen's inequality). The groups come from merging, one pair at a time,
// the two groups of a band whose merger adds least to the misses of the direct beam's transmission
// down to each interface and of the surface-reflected light's climb up from it, summed over the
// standard atmospheres of data/radiationBenchmark.json (TROP, MLS, MLW, SAW) with their humidity
// times 0.25, 1 and 2, at magnifications 1, 2 and 4; the 32 g-points go to the bands in the order of
// the misses (times each band's share of the beam) that one more would take away, and each group's
// scale is the one that makes its misses least. The report is the fold's miss of the direct beam's
// absorption and of the layer absorptions against nearInfraredBands, at those atmospheres.
//
// Each band also gets the Rayleigh scattering depth of its centre wavelength (Hansen and Travis
// 1974: 0.008569 l^-4 (1 + 0.0113 l^-2 + 0.00013 l^-4), l in um, at 1013.25 hPa) scaled so that the
// bands' shares weigh to the model's NEAR_INFRARED_RAYLEIGH, and the number of the band of Slingo's
// (1989) and Ebert and Curry's (1992) four (0.25-0.69, 0.69-1.19, 1.19-2.38, 2.38-4 um) that
// holds its centre.
import { writeFileSync } from 'node:fs';
import { NEAR_INFRARED_TABLE } from '../js/physics/nearInfraredTable.module.js';
import { nearInfraredBands, NEAR_INFRARED_MOLAR } from '../js/physics/shortwaveGases.module.js';
import { GREENHOUSE_GASES, NEAR_INFRARED_RAYLEIGH } from '../js/physics/radiation.module.js';
import { sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { BENCHMARK, modelColumn } from './standardAtmospheres.mjs';

const POINTS = Number(process.env.POINTS ?? 32), GRAVITY = 9.80665, D = 5 / 3, P_REF = NEAR_INFRARED_TABLE.pRef;
const FIT_ATMOSPHERES = ['TROP', 'MLS', 'MLW', 'SAW'], FIT_HUMIDITY = [0.25, 1, 2], TEST_HUMIDITY = [0.3, 1, 2], MAGNIFICATIONS = [1, 2, 4];
const SLINGO_EDGES = [0.25, 0.69, 1.19, 2.38, 4];
const MIXING = { co2: GREENHOUSE_GASES.carbonDioxide * NEAR_INFRARED_MOLAR.co2, ch4: GREENHOUSE_GASES.methane * NEAR_INFRARED_MOLAR.ch4, o2: 0.2095 * NEAR_INFRARED_MOLAR.o2 };
const levels = sigmaInterfaces('bl48'), K = levels.length - 1;

function airOf(atmosphere, humidity) {
  const column = modelColumn(atmosphere, levels), air = { K, vapour: new Float64Array(K), dry: new Float64Array(K), p: new Float64Array(K), T: new Float64Array(K), e: new Float64Array(K) };
  for (let k = 0; k < K; k++) {
    const mass = column.ps * (levels[k + 1] - levels[k]) / GRAVITY, p = column.ps * 0.5 * (levels[k] + levels[k + 1]), wet = Math.max(0, column.q[k] * humidity);
    air.vapour[k] = wet * mass; air.dry[k] = mass * (1 - wet); air.p[k] = p; air.T[k] = column.T[k]; air.e[k] = wet * p / (0.622 + 0.378 * wet);
  }
  return air;
}

function terms(band, air) {
  const out = {};
  for (const gas in band.gases) {
    const par = band.gases[gas], J = par.w.length, isVapor = gas === 'h2o';
    const depth = Array.from({ length: J }, () => new Float64Array(K)), self = new Float64Array(K);
    for (let k = 0; k < K; k++) {
      const u = isVapor ? air.vapour[k] : MIXING[gas] * air.dry[k], lnV = Math.log(Math.sqrt(air.p[k] ** 2 + par.pD ** 2) / P_REF), dT = air.T[k] - 250;
      for (let j = 0; j < J; j++) depth[j][k] = par.k[j] * u * Math.exp(par.n[j] * lnV + par.t[j] * dT);
      if (isVapor && par.self > 0) self[k] = par.self * u * air.e[k] / P_REF * Math.exp(par.selfT * Math.log(296 / air.T[k]));
    }
    out[gas] = { depth, self, weights: [1 - par.w.reduce((a, b) => a + b, 0), ...par.w] };
  }
  return out;
}

function combinations(band, air) {
  const base = terms(band, air), gases = Object.keys(band.gases), list = [];
  const walk = (g, picks, weight) => {
    if (g === gases.length) {
      const depth = new Float64Array(K);
      gases.forEach((gas, i) => { if (picks[i] > 0) for (let k = 0; k < K; k++) depth[k] += base[gas].depth[picks[i] - 1][k]; });
      for (const gas of gases) for (let k = 0; k < K; k++) depth[k] += base[gas].self[k];
      list.push({ picks: picks.slice(), weight, depth });
      return;
    }
    const w = base[gases[g]].weights;
    for (let j = 0; j < w.length; j++) { if (w[j] <= 1e-9) continue; picks[g] = j; walk(g + 1, picks, weight * w[j]); }
  };
  walk(0, [], 1);
  return { gases, list };
}

function transmissions(depth, m, out) {
  let total = 0;
  for (let k = 0; k < K; k++) total += depth[k];
  let s = 0;
  for (let l = 0; l <= K; l++) {
    if (l > 0) s += depth[l - 1];
    out.down[l] = Math.exp(-m * s);
    out.up[l] = Math.exp(-(m * total + D * (total - s)));
  }
}

const fitAirs = FIT_ATMOSPHERES.flatMap((name) => FIT_HUMIDITY.map((humidity) => airOf(BENCHMARK.atmospheres[name], humidity)));
const scratch = { down: new Float64Array(K + 1), up: new Float64Array(K + 1) };

function bandContext(band) {
  const per = fitAirs.map((air) => combinations(band, air));
  return { per, n: per[0].list.length, picks: per[0].list.map((c) => c.picks), weights: per[0].list.map((c) => c.weight), gases: per[0].gases };
}

function groupMisses(context, members, scale) {
  let cost = 0;
  for (const { list } of context.per) {
    let W = 0;
    for (const c of members) W += list[c].weight;
    const mean = new Float64Array(K);
    for (const c of members) for (let k = 0; k < K; k++) mean[k] += list[c].weight * list[c].depth[k] / W * scale;
    for (const m of MAGNIFICATIONS) {
      const down = new Float64Array(K + 1), up = new Float64Array(K + 1);
      for (const c of members) { transmissions(list[c].depth, m, scratch); for (let l = 0; l <= K; l++) { down[l] += list[c].weight * scratch.down[l]; up[l] += list[c].weight * scratch.up[l]; } }
      transmissions(mean, m, scratch);
      for (let l = 0; l <= K; l++) cost += Math.abs(down[l] - W * scratch.down[l]) + Math.abs(up[l] - W * scratch.up[l]);
    }
  }
  return cost / (fitAirs.length * MAGNIFICATIONS.length);
}

function merging(context) {
  let groups = Array.from({ length: context.n }, (_, i) => [i]);
  const memo = new Map(), cost = (g) => { const key = g.join(); if (!memo.has(key)) memo.set(key, g.length === 1 ? 0 : groupMisses(context, g, 1)); return memo.get(key); };
  const history = { [groups.length]: { cost: 0, groups: groups.map((g) => g.slice()) } };
  while (groups.length > 1) {
    let best = null;
    for (let a = 0; a < groups.length; a++) for (let b = a + 1; b < groups.length; b++) {
      const merged = groups[a].concat(groups[b]).sort((x, y) => x - y), delta = cost(merged) - cost(groups[a]) - cost(groups[b]);
      if (!best || delta < best.delta) best = { a, b, merged, delta };
    }
    groups = groups.filter((_, i) => i !== best.a && i !== best.b).concat([best.merged]);
    history[groups.length] = { cost: groups.reduce((s, g) => s + cost(g), 0), groups: groups.map((g) => g.slice()) };
  }
  return history;
}

function bestScale(context, members) {
  if (members.length === 1) return 1;
  let lo = 0.05, hi = 1.05;
  const phi = (Math.sqrt(5) - 1) / 2;
  let x1 = hi - phi * (hi - lo), x2 = lo + phi * (hi - lo), f1 = groupMisses(context, members, x1), f2 = groupMisses(context, members, x2);
  for (let it = 0; it < 30; it++) {
    if (f1 < f2) { hi = x2; x2 = x1; f2 = f1; x1 = hi - phi * (hi - lo); f1 = groupMisses(context, members, x1); } else { lo = x1; x1 = x2; f1 = f2; x2 = lo + phi * (hi - lo); f2 = groupMisses(context, members, x2); }
  }
  return 0.5 * (lo + hi);
}

const bands = NEAR_INFRARED_TABLE.bands.map((band) => { const context = bandContext(band); return { band, context, history: merging(context) }; });
const count = bands.map(() => 1);
const missesAt = (b, n) => bands[b].history[Math.min(n, bands[b].context.n)].cost * bands[b].band.share;
while (count.reduce((a, b) => a + b, 0) < POINTS) {
  let best = -1, gain = 0;
  for (let b = 0; b < bands.length; b++) { if (count[b] >= bands[b].context.n) continue; const g = missesAt(b, count[b]) - missesAt(b, count[b] + 1); if (g > gain) { gain = g; best = b; } }
  if (best < 0) break;
  count[best]++;
}

const centre = (band) => 1e4 / (0.5 * (band.range[0] + band.range[1]));
const rayleighAt = (um) => 0.008569 * um ** -4 * (1 + 0.0113 * um ** -2 + 0.00013 * um ** -4);
const shareSum = NEAR_INFRARED_TABLE.bands.reduce((s, b) => s + b.share, 0);
const rayleighNorm = NEAR_INFRARED_RAYLEIGH * shareSum / NEAR_INFRARED_TABLE.bands.reduce((s, b) => s + b.share * rayleighAt(centre(b)), 0);
const round = (x, digits = 5) => Number(x.toPrecision(digits));

const table = {
  pRef: P_REF,
  bands: bands.map(({ band, context, history }, b) => {
    const members = history[Math.min(count[b], context.n)].groups.map((g) => g.slice().sort((x, y) => x - y));
    const total = members.reduce((s, g) => s + g.reduce((t, c) => t + context.weights[c], 0), 0);
    const points = members.map((g) => {
      const weight = g.reduce((t, c) => t + context.weights[c], 0);
      const mix = context.gases.map((gas, i) => Array.from({ length: band.gases[gas].w.length }, (_, j) => round(g.reduce((t, c) => t + (context.picks[c][i] === j + 1 ? context.weights[c] : 0), 0) / weight)));
      return { weight: round(weight / total), scale: round(bestScale(context, g)), mix };
    });
    points.sort((x, y) => x.mix[0].reduce((a, v, j) => a + v * (j + 1), 0) - y.mix[0].reduce((a, v, j) => a + v * (j + 1), 0));
    return { band: band.band, share: band.share, gases: context.gases, slingo: Math.min(3, SLINGO_EDGES.findLastIndex((edge) => centre(band) >= edge)), rayleigh: round(rayleighNorm * rayleighAt(centre(band))), points };
  }),
};

function fold(air, m, taken, upward) {
  taken.fill(0); upward.fill(0);
  for (const band of table.bands) {
    const source = NEAR_INFRARED_TABLE.bands.find((b) => b.band === band.band), base = terms(source, air);
    const down = new Float64Array(K + 1), up = new Float64Array(K + 1);
    for (const point of band.points) {
      const depth = new Float64Array(K);
      band.gases.forEach((gas, i) => { for (let j = 0; j < point.mix[i].length; j++) for (let k = 0; k < K; k++) depth[k] += point.mix[i][j] * base[gas].depth[j][k]; });
      for (let k = 0; k < K; k++) depth[k] = point.scale * (depth[k] + (base.h2o.self[k]));
      transmissions(depth, m, scratch);
      for (let l = 0; l <= K; l++) { down[l] += point.weight * scratch.down[l]; up[l] += point.weight * scratch.up[l]; }
    }
    for (let k = 0; k < K; k++) { taken[k] += band.share * (1 - down[k + 1]); upward[k] += band.share * (up[k + 1] - up[k]); }
  }
}

const rows = [], taken = new Float64Array(K), upward = new Float64Array(K), takenFolded = new Float64Array(K), upwardFolded = new Float64Array(K);
for (const name of FIT_ATMOSPHERES) for (const humidity of TEST_HUMIDITY) {
  const air = airOf(BENCHMARK.atmospheres[name], humidity);
  for (const m of MAGNIFICATIONS) {
    nearInfraredBands(air, m, D, MIXING, 1, taken, upward);
    fold(air, m, takenFolded, upwardFolded);
    let layer = 0, up = 0, upTotal = 0;
    for (let k = 0; k < K; k++) {
      layer = Math.max(layer, Math.abs((taken[k] - (k ? taken[k - 1] : 0)) - (takenFolded[k] - (k ? takenFolded[k - 1] : 0))));
      up += upward[k] - upwardFolded[k]; upTotal += upward[k];
    }
    rows.push({ name, humidity, m, total: taken[K - 1], miss: takenFolded[K - 1] - taken[K - 1], layer, up, upTotal });
  }
}
const W = 1360;
console.log(`allocation ${table.bands.map((b) => `${b.band}:${b.points.length}`).join(' ')} (${table.bands.reduce((s, b) => s + b.points.length, 0)} g-points)`);
for (const r of rows) if (r.m === 1) console.log(`${r.name.padEnd(4)} humidity x${r.humidity}: direct absorption ${(r.total * W).toFixed(1)} W/m2 of 1360, fold ${(r.miss * W).toFixed(2)}, worst layer ${(r.layer * W).toFixed(2)}, reflected light's ${(r.upTotal * W).toFixed(1)}, fold ${(r.up * W).toFixed(2)}`);
console.log(`direct absorption miss: mean |${(rows.reduce((s, r) => s + Math.abs(r.miss), 0) / rows.length * W).toFixed(2)}| worst |${(Math.max(...rows.map((r) => Math.abs(r.miss))) * W).toFixed(2)}| W/m2 per 1360; worst layer ${(Math.max(...rows.map((r) => r.layer)) * W).toFixed(2)}`);

if (process.env.WRITE) {
  const header = '// Generated by scripts/layeredFit.mjs: for each of RRTMG_SW\'s bands 16-24 and 29, its share of the solar beam, its gases in\n// the order of the mix lists, the Rayleigh depth of its centre wavelength at 1013.25 hPa, the number (0-3) of the band of Slingo (1989)\n// and Ebert and Curry (1992) that holds its centre and its g-points, each with its weight in the band, the scale of its depth and, for each\n// gas, the shares of its weight that the gas\'s k-distribution terms (nearInfraredTable.module.js) take.\n';
  writeFileSync(new URL('../js/physics/layeredTable.module.js', import.meta.url), `${header}export const LAYERED_TABLE = ${JSON.stringify(table)};\n`);
  console.log('wrote js/physics/layeredTable.module.js');
}
