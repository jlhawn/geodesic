// The cloud-radiative effects of a saved state by cloud class and region,
// on the single-thread CPU engine with the ocean off:
//   node scripts/cloudClasses.mjs <state.bin>
// One model step is taken first (the plumes' cumulus is not saved), then
// the state is held fixed and lit at TIMES (24) instants spread over day
// DAY (the day ending DAY days after the equinox; the state's day by
// default). At each instant every column's radiation runs once with all
// its cloud (and its clear-sky pass) and once with each class taken away
// (the radiation's cloud mask), the deck held at the cover and water the
// pass at the day's start diagnosed. A class's effect is the full pass's flux less the
// pass without it: shortwave, absorbed sunlight with it less without it;
// longwave, outgoing longwave without it less with it. The classes:
// the deck; the plumes' cumulus; and the resolved cloud (q_c) in runs of
// adjacent cloudy layers, each run classed by the pressure of its top
// layer's midpoint: low below 680 hPa, middle 680-440 hPa, high above
// 440 hPa. Effects do not add: the sum of the classes is printed beside
// the total.
// Cover: the radiation's own, each layer's cover times its visibility
// overlapped maximum-random as the radiation does, of the class's layers
// alone; the deck's is its fraction; "all" the column's resolved and
// cumulus cover combined at random with the deck's, as the radiation's two
// columns are. Water paths are grid means and in-cloud (over the class's
// cover), split by the layer's temperature: warmer than 273 K, 273-235 K,
// colder than 235 K. tau is the in-cloud mid-visible optical depth of the
// physical optics of cloudOptics (droplet and crystal radii by the
// radiation's options), its distribution over the class's cover in the
// ISCCP bins below 3.6, 3.6-23 and above 23; "two-stream" is the depth
// the radiation's shortwave uses (cloudScattering times the path when that
// is set, else the sum of each layer's (1 - g) tau), in-cloud.
// Regions: global, 30S-30N, 30-60N, 30-60S, 60-90N, 60-90S, land, sea
// (sea ice included). RADIATION, MOIST, BOUNDARY_LAYER and SURFACE (JSON)
// pass options as to scripts/verticalAudit.mjs.
//
// Earth (annual means): total cloud 0.65-0.68 (ISCCP, MODIS, CALIPSO);
// high cloud 0.2-0.3 by passive sensors (ISCCP 0.20, MODIS ~0.27); liquid
// water path over the oceans 50-90 g/m2 (O'Dell et al. 2008, MAC-LWP);
// ice water path of order 20-70 g/m2 (CloudSat 2C-ICE, Waliser et al.
// 2009); cloud effects -47 +- 4 and +26 +- 3 W/m2 (CERES EBAF).
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { createModel, STATE_NAMES } from '../js/model.module.js';
import { decodeState, savedLevels } from '../js/stateFile.module.js';
import { savedDeckField, DECK_FIELDS } from '../js/physics/regrid.module.js';
import { DAY, VISIBLE_PATH, LOW_CLOUD_PRESSURE, CLOUD_OPTICS, cloudOptics } from '../js/physics/radiation.module.js';
import { FREEZING_POINT } from '../js/physics/ice.module.js';

const FILE = process.argv[2];
if (!FILE) { console.error('usage: node scripts/cloudClasses.mjs <state.bin>'); process.exit(1); }
const RADIATION = JSON.parse(process.env.RADIATION ?? '{}'), MOIST = JSON.parse(process.env.MOIST ?? '{}'), BOUNDARY_LAYER = JSON.parse(process.env.BOUNDARY_LAYER ?? '{}'), SURFACE = JSON.parse(process.env.SURFACE ?? '{}');
const TIMES = Number(process.env.TIMES ?? 24);
const HIGH_CLOUD_PRESSURE = 440e2;
const t0 = performance.now();

const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
const AT = Number(process.env.DAY ?? saved.day);
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const model = createModel(new Grid(saved.N), { topography, levels: savedLevels(saved), ocean: false, radiation: { clearSkyPass: true, ...RADIATION }, moist: MOIST, boundaryLayer: BOUNDARY_LAYER, surface: SURFACE });
const { mesh, core, state, radiation, boundaryLayer: bl, moist, seaIce, land, geography } = model;
const { K, C, sigmaMid, dSigma, g, exnerLayer } = core.diagnostics;
STATE_NAMES.forEach((name, a) => state[a].set(saved[name]));
seaIce.load(state[6], saved.concentration ?? null);
for (const field of Object.keys(DECK_FIELDS)) radiation[field].set(savedDeckField(saved, field, model));
land.load({ soil: Float64Array.from(saved.land.soil), snow: Float64Array.from(saved.land.snow), ...(saved.land.vegetation ? { vegetation: Float64Array.from(saved.land.vegetation) } : {}), ...(saved.land.surface ? { surface: Float64Array.from(saved.land.surface) } : {}) }, state[6]);
model.time = saved.time;
const zs = Float64Array.from({ length: C }, (_, i) => (model.surfaceGeopotential ? model.surfaceGeopotential[i] / g : 0));
if (saved.boundaryDepth) bl.depth.set(Float64Array.from(saved.boundaryDepth, (z, i) => z + zs[i]));
if (saved.mixingTop) bl.mixingTop.set(Float64Array.from(saved.mixingTop, (z, i) => z + zs[i]));
if (saved.boundaryRegime) bl.regime.set(saved.boundaryRegime);
if (saved.boundaryBuoyancy) bl.buoyancyFlux.set(saved.boundaryBuoyancy);
model.step(1350 * 16 / saved.N);

const [pi, theta, , surfaceT, q, qc, ice] = state;
const landMask = geography.land, deg = 180 / Math.PI, area = mesh.areaCell;
const fluxT = new Float64Array(C), adir = new Float64Array(C), adif = new Float64Array(C), openSea = new Float64Array(C);
const fluxState = [pi, theta, state[2], fluxT, q, qc];
const out = [null, new Float64Array(K * C), null, null, new Float64Array(K * C)];
const temperature = Float64Array.from({ length: K * C }, (_, idx) => theta[idx] * exnerLayer[idx]);
const mass = (k, i) => pi[i] * dSigma[k] / g;
const pressure = (k, i) => pi[i] * sigmaMid[k];
const cuCover = moist.cumulusCover, cuWater = moist.cumulusWater;

function surfaces() {
  for (let i = 0; i < C; i++) {
    fluxT[i] = surfaceT[i];
    if (landMask[i]) { adir[i] = adif[i] = land.albedo(i); openSea[i] = 0; continue; }
    const h = ice[i], cover = seaIce.cover(i, h), mu = radiation.cosZenith(i);
    adir[i] = seaIce.albedo(h, mu, seaIce.snow[i], cover); adif[i] = seaIce.albedo(h, null, seaIce.snow[i], cover);
    openSea[i] = 1 - cover;
    if (h > 0 && cover < 1) fluxT[i] = cover * surfaceT[i] + (1 - cover) * FREEZING_POINT;
  }
}

const CLASSES = ['deck', 'cumulus', 'low', 'middle', 'high'];
const runClass = new Int8Array(K * C).fill(-1);
for (let i = 0; i < C; i++) {
  for (let k = 0; k < K;) {
    if (!(qc[k * C + i] > 0)) { k++; continue; }
    const p = pressure(k, i), c = p > LOW_CLOUD_PRESSURE ? 2 : p > HIGH_CLOUD_PRESSURE ? 3 : 4;
    for (; k < K && qc[k * C + i] > 0; k++) runClass[k * C + i] = c;
  }
}

const ones = new Float64Array(K * C).fill(1), zeros = new Float64Array(K * C);
const fullDeck = { fraction: null, deck: null };
function mask(without) {
  const resolved = Float64Array.from({ length: K * C }, (_, idx) => (runClass[idx] === without ? 0 : 1));
  return { resolved: without >= 2 ? resolved : ones, cumulus: without === 1 ? zeros : ones, fraction: without === 0 ? new Float64Array(C) : fullDeck.fraction, deck: without === 0 ? new Float64Array(C) : fullDeck.deck };
}
function perCell(m) {
  if (!m) return null;
  const column = { resolved: new Float64Array(K), cumulus: new Float64Array(K), fraction: 0, deck: 0 };
  return (i) => { for (let k = 0; k < K; k++) { column.resolved[k] = m.resolved[k * C + i]; column.cumulus[k] = m.cumulus[k * C + i]; } column.fraction = m.fraction[i]; column.deck = m.deck[i]; return column; };
}
function pass(m) {
  const columnMask = perCell(m);
  radiation.restartSums();
  for (let i = 0; i < C; i++) {
    radiation.useCloudMask(columnMask ? columnMask(i) : null);
    radiation.apply(fluxState, out, model.surface.windSpeed, null, i, i + 1, adir, adif, null, openSea, bl.depth, 0);
  }
  radiation.useCloudMask(null);
  return { absorbed: Float64Array.from(radiation.summed.absorbedSolar), outgoing: Float64Array.from(radiation.summed.outgoingLongwave), clearAbsorbed: Float64Array.from(radiation.summed.clearAbsorbedSolar), clearOutgoing: Float64Array.from(radiation.summed.clearOutgoingLongwave), insolation: Float64Array.from(radiation.summed.insolation) };
}

const surfaceOf = (i) => (landMask[i] && !(geography.iceSheet && geography.iceSheet[i]) ? 1 : 0);
const optics = { ...CLOUD_OPTICS, ...Object.fromEntries(Object.keys(CLOUD_OPTICS).filter((key) => key in RADIATION).map((key) => [key, RADIATION[key]])) };
const layerOptics = { liquid: 0, visible: 0, solar: 0, infrared: 0 };
const BANDS = 3, band = (T) => (T > 273.15 ? 0 : T > 235.15 ? 1 : 2);
const cover = CLASSES.map(() => new Float64Array(C)), allCover = new Float64Array(C);
const water = CLASSES.map(() => Array.from({ length: BANDS }, () => new Float64Array(C)));
const liquidPath = CLASSES.map(() => new Float64Array(C)), visible = CLASSES.map(() => new Float64Array(C)), twoStream = CLASSES.map(() => new Float64Array(C));
surfaces();
radiation.setTime((AT - 1) * DAY);
let coverCheck = 0;
for (let i = 0; i < C; i++) {
  radiation.apply(fluxState, out, model.surface.windSpeed, null, i, i + 1, adir, adif, null, openSea, bl.depth, 0);
  const layerCover = Float64Array.from(radiation.layerCover), resolvedCover = new Float64Array(K);
  const full = radiation.budget.cloudCover, fraction = radiation.stratusFraction[i], deck = radiation.stratus[i];
  if (!fullDeck.fraction) { fullDeck.fraction = new Float64Array(C); fullDeck.deck = new Float64Array(C); }
  fullDeck.fraction[i] = fraction; fullDeck.deck[i] = deck;
  radiation.useCloudMask({ resolved: Float64Array.from({ length: K }, () => 1), cumulus: new Float64Array(K), fraction, deck });
  radiation.apply(fluxState, out, model.surface.windSpeed, null, i, i + 1, adir, adif, null, openSea, bl.depth, 0);
  resolvedCover.set(radiation.layerCover);
  radiation.useCloudMask(null);
  const blocks = CLASSES.map(() => ({ block: 0, clear: 1 })), together = { block: 0, clear: 1 };
  const close = (b, seen, last) => { if (seen > 0) b.block = Math.max(b.block, seen); if (b.block > 0 && (!(seen > 0) || last)) { b.clear *= 1 - b.block; b.block = 0; } };
  for (let k = 0; k < K; k++) {
    const idx = k * C + i, m = mass(k, i), T = temperature[idx];
    cloudOptics(T, surfaceOf(i) === 1, optics, layerOptics);
    const resolved = Math.max(0, qc[idx]) * m, cumulus = cuCover[idx] * cuWater[idx] * m, c = runClass[idx];
    const seenAll = resolved + cumulus > 0 ? layerCover[k] * -Math.expm1(-(resolved + cumulus) / VISIBLE_PATH) : 0;
    close(together, seenAll, k === K - 1);
    for (let n = 1; n < CLASSES.length; n++) {
      const w = n === 1 ? cumulus : c === n ? resolved : 0, f = n === 1 ? cuCover[idx] : resolvedCover[k];
      close(blocks[n], w > 0 ? f * -Math.expm1(-w / VISIBLE_PATH) : 0, k === K - 1);
      if (w > 0) { water[n][band(T)][i] += w; liquidPath[n][i] += layerOptics.liquid * w; visible[n][i] += layerOptics.visible * w; twoStream[n][i] += radiation.solarDepth[k] * w; }
    }
    if (k === radiation.stratusLayer && fraction > 0) { water[0][band(T)][i] += fraction * deck; liquidPath[0][i] += layerOptics.liquid * fraction * deck; visible[0][i] += layerOptics.visible * fraction * deck; twoStream[0][i] += radiation.solarDepth[k] * fraction * deck; }
  }
  for (let n = 1; n < CLASSES.length; n++) cover[n][i] = 1 - blocks[n].clear;
  cover[0][i] = fraction;
  const columnCover = 1 - together.clear;
  if (full > 0) coverCheck = Math.max(coverCheck, Math.abs(columnCover - full));
  allCover[i] = 1 - (1 - columnCover) * (1 - fraction);
}

const sums = { full: null, without: CLASSES.map(() => null) };
const add = (into, from) => { if (!into) return Object.fromEntries(Object.entries(from).map(([key, x]) => [key, Float64Array.from(x)])); for (const key of Object.keys(from)) for (let i = 0; i < C; i++) into[key][i] += from[key][i]; return into; };
const masks = CLASSES.map((_, n) => mask(n));
let checkFull = 0, checkMean = 0;
for (let n = 0; n < TIMES; n++) {
  radiation.setTime((AT - 1) * DAY + (n + 0.5) * DAY / TIMES);
  surfaces();
  const full = pass(mask(-1)), free = n === 0 ? pass(null) : null;
  if (free) {
    let sum = 0, total = 0;
    for (let i = 0; i < C; i++) { checkFull = Math.max(checkFull, Math.abs(full.absorbed[i] - free.absorbed[i]), Math.abs(full.outgoing[i] - free.outgoing[i])); sum += area[i] * (full.absorbed[i] - free.absorbed[i] - full.outgoing[i] + free.outgoing[i]); total += area[i]; }
    checkMean = sum / total;
  }
  sums.full = add(sums.full, full);
  masks.forEach((m, c) => { sums.without[c] = add(sums.without[c], pass(m)); });
}

const REGIONS = [
  ['global', () => true], ['30S-30N', (lat) => lat >= -30 && lat <= 30], ['30-60N', (lat) => lat > 30 && lat <= 60], ['30-60S', (lat) => lat < -30 && lat >= -60],
  ['60-90N', (lat) => lat > 60], ['60-90S', (lat) => lat < -60], ['land', (lat, i) => landMask[i] === 1 || landMask[i] === true], ['sea', (lat, i) => !landMask[i]],
];
const f = (x, d = 1) => (Number.isFinite(x) ? x.toFixed(d) : 'n/a');
const mean = (x, inside) => { let s = 0, a = 0; for (let i = 0; i < C; i++) if (inside[i]) { s += area[i] * x(i); a += area[i]; } return s / a; };
console.log(`cloud classes of ${FILE.split('/').pop()} (N=${saved.N}, K=${K}) after one step, lit over day ${AT} at ${TIMES} instants; RADIATION ${JSON.stringify(RADIATION)}`);
console.log(`checks: column cover as the radiation's to ${coverCheck.toExponential(1)}; the deck held against the deck re-diagnosed at the first instant: ${checkFull.toExponential(1)} W/m2 in the column most changed, ${checkMean.toExponential(1)} W/m2 in the global ASR - OLR; ${((performance.now() - t0) / 1000).toFixed(0)} s`);
for (const [name, test] of REGIONS) {
  const inside = Uint8Array.from({ length: C }, (_, i) => (test(mesh.latCell[i] * deg, i) ? 1 : 0));
  const S = sums.full, sw = mean((i) => (S.absorbed[i] - S.clearAbsorbed[i]) / TIMES, inside), lw = mean((i) => (S.clearOutgoing[i] - S.outgoing[i]) / TIMES, inside);
  console.log(`\n${name}: SWCRE ${f(sw)} LWCRE ${f(lw)} W/m2, cloud cover ${f(mean((i) => allCover[i], inside), 3)}, ASR ${f(mean((i) => S.absorbed[i] / TIMES, inside))} OLR ${f(mean((i) => S.outgoing[i] / TIMES, inside))}, clear-sky albedo ${f(1 - mean((i) => S.clearAbsorbed[i], inside) / mean((i) => S.insolation[i], inside), 3)}`);
  console.log('class    cover  SW     LW    | grid path g/m2 >273/273-235/<235 (total) | in-cloud g/m2 | tau in-cloud  share <3.6 3.6-23 >23 | two-stream | liquid / ice g/m2');
  let sumSw = 0, sumLw = 0, sumLiquid = 0;
  const sumPaths = [0, 0, 0];
  CLASSES.forEach((cls, c) => {
    const W = sums.without[c];
    const csw = mean((i) => (S.absorbed[i] - W.absorbed[i]) / TIMES, inside), clw = mean((i) => (W.outgoing[i] - S.outgoing[i]) / TIMES, inside);
    sumSw += csw; sumLw += clw;
    const cv = mean((i) => cover[c][i], inside), paths = water[c].map((w) => 1000 * mean((i) => w[i], inside)), total = paths.reduce((a, b) => a + b, 0);
    paths.forEach((x, b) => { sumPaths[b] += x; });
    sumLiquid += 1000 * mean((i) => liquidPath[c][i], inside);
    const bins = [0, 0, 0];
    let weight = 0;
    for (let i = 0; i < C; i++) { if (!inside[i] || !(cover[c][i] > 0)) continue; const tau = visible[c][i] / cover[c][i], w = area[i] * cover[c][i]; weight += w; bins[tau < 3.6 ? 0 : tau < 23 ? 1 : 2] += w; }
    const tauMean = mean((i) => visible[c][i], inside) / cv, depth = mean((i) => twoStream[c][i], inside) / cv;
    console.log(`${cls.padEnd(8)} ${f(cv, 3)} ${f(csw).padStart(6)} ${f(clw).padStart(5)} | ${paths.map((x) => f(x)).join(' / ')} (${f(total)}) | ${f(total / cv, 0).padStart(5)} | ${f(tauMean).padStart(6)}  ${bins.map((b) => f(b / weight, 2)).join(' ')} | ${f(depth)} | ${f(1000 * mean((i) => liquidPath[c][i], inside))} / ${f(total - 1000 * mean((i) => liquidPath[c][i], inside))}`);
  });
  const allPath = sumPaths.reduce((a, b) => a + b, 0);
  console.log(`sum of classes   SW ${f(sumSw)} LW ${f(sumLw)}; all cloud: grid path ${sumPaths.map((x) => f(x)).join(' / ')} (${f(allPath)}) g/m2, liquid ${f(sumLiquid)} ice ${f(allPath - sumLiquid)} g/m2`);
}
