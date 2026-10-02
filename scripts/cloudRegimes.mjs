// The cloud of a saved state by regime and class, beside observed values,
// on the single-thread CPU engine with the ocean off:
//   node scripts/cloudRegimes.mjs <state.bin>
// One model step is taken first (the plumes' cumulus is not saved), then
// every column's radiation runs once with all its cloud and once without
// the cumulus, at the state's time, to read the radiation's own covers.
// The classes are those of scripts/cloudClasses.mjs: the deck, the plumes'
// cumulus, and the resolved cloud (q_c) in runs of adjacent cloudy layers
// classed by their top layer's pressure (low below 680 hPa, middle
// 680-440 hPa, high above 440 hPa); covers are overlapped as the radiation
// overlaps them (cloudOverlap), "total" combines the resolved and cumulus
// cover at random with the deck's, and is given beside under
// maximum-random and exponential-random overlap. Per class: cover; grid-mean and in-cloud (over the
// class's cover) path; the share of its cover whose in-cloud mid-visible
// optical depth (cloudOptics) lies below 3.6, 3.6-23 and above 23 (ISCCP's
// thin, medium and thick bins). The humidity lines: relative humidity over
// liquid water (Bolton, as the model's) and over ice (the IFS form
// e_i = 611.21 exp(22.587 (T - 273.16)/(T + 0.7)) Pa) by layer mass and
// area in the upper troposphere (150-350 hPa) and below (350-700, 700 hPa
// to the surface), the share of the upper troposphere's layer area above
// ice saturation, and the shares of humid layers that hold no condensate:
// upper troposphere at RHi above 1, lower at RH above 0.9. "cumulus
// updraught" is the largest layer cumulus fraction of each column, the
// plumes' active area M/(rho w_u); the resolved condensate is split into
// liquid and ice by the optics' phase ramp, and the high class's
// condensate colder than 235 K is given apart. The saved day means of
// the cloud effects are printed where the state carries them.
// RADIATION, MOIST, BOUNDARY_LAYER and SURFACE (JSON) pass options as to
// scripts/verticalAudit.mjs. Earth's values are printed beside each regime
// with their sources; "from memory" marks a value not looked up for this
// table.
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { createModel, STATE_NAMES } from '../js/model.module.js';
import { decodeState, savedLevels } from '../js/stateFile.module.js';
import { savedDeckField, DECK_FIELDS } from '../js/physics/regrid.module.js';
import { LOW_CLOUD_PRESSURE, CLOUD_OPTICS, cloudOptics, overlapped, DECORRELATION_LENGTH, DECORRELATION_SLOPE } from '../js/physics/radiation.module.js';
import { saturationHumidity } from '../js/physics/moist.module.js';
import { BOXES, inLongitudes } from '../js/audit.module.js';

const FILE = process.argv[2];
if (!FILE) { console.error('usage: node scripts/cloudRegimes.mjs <state.bin>'); process.exit(1); }
const RADIATION = JSON.parse(process.env.RADIATION ?? '{}'), MOIST = JSON.parse(process.env.MOIST ?? '{}'), BOUNDARY_LAYER = JSON.parse(process.env.BOUNDARY_LAYER ?? '{}'), SURFACE = JSON.parse(process.env.SURFACE ?? '{}');
const HIGH_CLOUD_PRESSURE = 440e2, VISIBLE = 1e-3;
const t0 = performance.now();

const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
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
const dt = 1350 * 16 / saved.N;
const [pi, theta, , surfaceT, q, qc, ice] = state;
const mass = (k, i) => pi[i] * dSigma[k] / g;
model.step(dt);

const landMask = geography.land, deg = 180 / Math.PI, area = mesh.areaCell;
const lat = Float64Array.from(mesh.latCell, (x) => x * deg), lon = Float64Array.from(mesh.lonCell, (x) => x * deg);
const fluxT = new Float64Array(C), adir = new Float64Array(C), adif = new Float64Array(C), openSea = new Float64Array(C);
const fluxState = [pi, theta, state[2], fluxT, q, qc];
const out = [null, new Float64Array(K * C), null, null, new Float64Array(K * C)];
const temperature = Float64Array.from({ length: K * C }, (_, idx) => theta[idx] * exnerLayer[idx]);
const pressure = (k, i) => pi[i] * sigmaMid[k];
const iceSaturation = (T, p) => { const e = 611.21 * Math.exp(22.587 * (T - 273.16) / (T + 0.7)); return 0.622 * e / Math.max(1, p - 0.378 * e); };
const cuCover = moist.cumulusCover, cuWater = moist.cumulusWater;
for (let i = 0; i < C; i++) {
  fluxT[i] = surfaceT[i];
  if (landMask[i]) { adir[i] = adif[i] = land.albedo(i); openSea[i] = 0; continue; }
  const h = ice[i], cover = seaIce.cover(i, h), mu = radiation.cosZenith(i);
  adir[i] = seaIce.albedo(h, mu, seaIce.snow[i], cover); adif[i] = seaIce.albedo(h, null, seaIce.snow[i], cover);
  openSea[i] = 1 - cover;
  if (h > 0 && cover < 1) fluxT[i] = cover * surfaceT[i] + (1 - cover) * 271.35;
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
const optics = { ...CLOUD_OPTICS, ...Object.fromEntries(Object.keys(CLOUD_OPTICS).filter((key) => key in RADIATION).map((key) => [key, RADIATION[key]])) };
const layerOptics = { liquid: 0, visible: 0, solar: 0, infrared: 0 };
const cover = CLASSES.map(() => new Float64Array(C)), allCover = new Float64Array(C), path = CLASSES.map(() => new Float64Array(C)), visible = CLASSES.map(() => new Float64Array(C));
const coldHigh = new Float64Array(C), updraught = new Float64Array(C), icePath = new Float64Array(C), liquidPath = new Float64Array(C), exponentialCover = new Float64Array(C), randomCover = new Float64Array(C);
const exponential = (RADIATION.cloudOverlap ?? 'exponentialRandom') === 'exponentialRandom';
const humidity = { upper: [new Float64Array(C), new Float64Array(C), new Float64Array(C)], middle: [new Float64Array(C), new Float64Array(C)], lower: [new Float64Array(C), new Float64Array(C)] };
const issr = new Float64Array(C), upperArea = new Float64Array(C), humidUpper = new Float64Array(C), humidUpperDry = new Float64Array(C), humidLower = new Float64Array(C), humidLowerDry = new Float64Array(C);
radiation.setTime(model.time);
for (let i = 0; i < C; i++) {
  radiation.apply(fluxState, out, model.surface.windSpeed, null, i, i + 1, adir, adif, null, openSea, bl.depth, 0);
  const layerCover = Float64Array.from(radiation.layerCover), fraction = radiation.stratusFraction[i], deck = radiation.stratus[i];
  radiation.useCloudMask({ resolved: new Float64Array(K).fill(1), cumulus: new Float64Array(K), fraction, deck });
  radiation.apply(fluxState, out, model.surface.windSpeed, null, i, i + 1, adir, adif, null, openSea, bl.depth, 0);
  const resolvedCover = Float64Array.from(radiation.layerCover);
  radiation.useCloudMask(null);
  const blocks = CLASSES.map(() => ({ block: 0, clear: 1, above: 0, cumulative: 0 })), together = { block: 0, clear: 1 };
  const decorrelation = (RADIATION.decorrelationLength ?? DECORRELATION_LENGTH) - (RADIATION.decorrelationSlope ?? DECORRELATION_SLOPE) * Math.abs(lat[i]);
  let above = 0, cumulative = 0, alpha = 0;
  const close = (b, seen, last) => {
    if (exponential && b !== together) { b.cumulative = overlapped(b.cumulative, b.above, seen, alpha); b.above = seen; b.clear = 1 - b.cumulative; return; }
    if (seen > 0) b.block = Math.max(b.block, seen);
    if (b.block > 0 && (!(seen > 0) || last)) { b.clear *= 1 - b.block; b.block = 0; }
  };
  for (let k = 0; k < K; k++) {
    const idx = k * C + i, m = mass(k, i), T = temperature[idx], p = pressure(k, i);
    cloudOptics(T, landMask[i] && !(geography.iceSheet && geography.iceSheet[i]), optics, layerOptics);
    const resolved = Math.max(0, qc[idx]) * m, cumulus = cuCover[idx] * cuWater[idx] * m, c = runClass[idx];
    liquidPath[i] += layerOptics.liquid * resolved; icePath[i] += (1 - layerOptics.liquid) * resolved;
    const seen = resolved + cumulus > 0 ? layerCover[k] * -Math.expm1(-(resolved + cumulus) / VISIBLE) : 0;
    alpha = k > 0 ? Math.exp(-(core.diagnostics.geopotential[(k - 1) * C + i] - core.diagnostics.geopotential[idx]) / g / decorrelation) : 0;
    close(together, seen, k === K - 1);
    cumulative = overlapped(cumulative, above, seen, alpha);
    above = seen;
    updraught[i] = Math.max(updraught[i], cuCover[idx]);
    for (let n = 1; n < CLASSES.length; n++) {
      const w = n === 1 ? cumulus : c === n ? resolved : 0, f = n === 1 ? cuCover[idx] : resolvedCover[k];
      close(blocks[n], w > 0 ? f * -Math.expm1(-w / VISIBLE) : 0, k === K - 1);
      if (w > 0) { path[n][i] += w; visible[n][i] += layerOptics.visible * w; if (n === 4 && T < 235.15) coldHigh[i] += w; }
    }
    if (k === radiation.stratusLayer && fraction > 0) { path[0][i] += fraction * deck; visible[0][i] += layerOptics.visible * fraction * deck; }
    const rh = Math.max(0, q[idx]) / saturationHumidity(T, p), rhi = Math.max(0, q[idx]) / iceSaturation(T, p), cloudy = qc[idx] > 0 || cuCover[idx] > 0;
    if (p >= 150e2 && p < 350e2) {
      humidity.upper[0][i] += m * rh; humidity.upper[1][i] += m * rhi; humidity.upper[2][i] += m;
      upperArea[i] += 1; if (rhi > 1) issr[i] += 1;
      if (rhi > 1) { humidUpper[i] += 1; if (!cloudy) humidUpperDry[i] += 1; }
    } else if (p >= 350e2 && p < 700e2) { humidity.middle[0][i] += m * rh; humidity.middle[1][i] += m; } else if (p >= 700e2) {
      humidity.lower[0][i] += m * rh; humidity.lower[1][i] += m;
      if (rh > 0.9) { humidLower[i] += 1; if (!cloudy) humidLowerDry[i] += 1; }
    }
  }
  for (let n = 1; n < CLASSES.length; n++) cover[n][i] = 1 - blocks[n].clear;
  cover[0][i] = fraction;
  randomCover[i] = 1 - together.clear * (1 - fraction);
  exponentialCover[i] = 1 - (1 - cumulative) * (1 - fraction);
  allCover[i] = exponential ? exponentialCover[i] : randomCover[i];
}

const sea = (i) => !landMask[i];
const box = ([south, north, west, east], keep = () => true) => (i) => keep(i) && lat[i] >= south && lat[i] <= north && inLongitudes(lon[i], west, east);
const REGIMES = [
  ['global', () => true, 'total 0.65-0.68 (ISCCP, MODIS, CALIPSO); high 0.2-0.3; thin share of high ~0.6'],
  ['warm pool 10S-10N 120-170E sea', box([-10, 10, 120, 170], sea), 'total 0.80-0.90, high 0.55-0.70, thin share of high ~0.5 (ISCCP, CALIPSO; from memory)'],
  ['Pacific ITCZ 5-12N 160E-100W', box(BOXES.itcz), 'total 0.70-0.85, high 0.45-0.60 (ISCCP, CALIPSO; from memory)'],
  ['N Pacific trades 15-25N 170-130W', box([15, 25, -170, -130], sea), 'total 0.35-0.55, low 0.2-0.4 (MODIS, CALIPSO-GOCCP), updraught area 0.02-0.05 (LES of BOMEX; from memory)'],
  ['S Pacific trades 10-20S 160-120W', box([-20, -10, -160, -120], sea), 'as the northern trades'],
  ['Atlantic trades 10-20N 50-25W', box([10, 20, -50, -25], sea), 'as the northern trades'],
  ['SE Pacific 10-30S 110-80W', box(BOXES.sePacific, sea), 'low 0.6-0.8 (Klein and Hartmann 1993, SON)'],
  ['Peru 5-20S 90-75W', box(BOXES.peru, sea), 'low 0.6-0.8 (Klein and Hartmann 1993)'],
  ['Namibia 10-20S 0-10E', box(BOXES.namibia, sea), 'low 0.6-0.8 (Klein and Hartmann 1993)'],
  ['California 20-30N 130-120W', box(BOXES.california, sea), 'low 0.5-0.7 (Klein and Hartmann 1993, JJA higher)'],
  ['Southern Ocean 40-60S', box([-60, -40, -180, 180], sea), 'total 0.80-0.90, low 0.5-0.7 (CloudSat-CALIPSO; from memory)'],
  ['N Atlantic 40-60N 50-10W', box([40, 60, -50, -10], sea), 'total 0.75-0.85 (ISCCP; from memory)'],
  ['N Pacific 40-60N 150E-140W', box([40, 60, 150, -140], sea), 'total 0.80-0.90 (ISCCP; from memory)'],
  ['60-90N', (i) => lat[i] > 60, 'total 0.80-0.90 in September, 0.75-0.85 in June (CALIPSO, Kay and Gettelman 2009; from memory)'],
  ['60-90S', (i) => lat[i] < -60, 'total 0.65-0.80 (CALIPSO; from memory)'],
  ['land 60S-60N', (i) => landMask[i] && Math.abs(lat[i]) <= 60, 'total 0.50-0.60 (MODIS, King et al. 2013; from memory)'],
  ['sea 60S-60N', (i) => !landMask[i] && Math.abs(lat[i]) <= 60, 'total 0.68-0.75 (MODIS, King et al. 2013; from memory)'],
];
const f = (x, d = 2) => (Number.isFinite(x) ? x.toFixed(d) : 'n/a');
const mean = (x, inside) => { let s = 0, a = 0; for (let i = 0; i < C; i++) if (inside[i]) { s += area[i] * x(i); a += area[i]; } return s / a; };
const ratio = (num, den, inside) => mean(num, inside) / mean(den, inside);
console.log(`cloud regimes of ${FILE.split('/').pop()} (N=${saved.N}, K=${K}, day ${saved.day}) after one step; RADIATION ${JSON.stringify(RADIATION)} MOIST ${JSON.stringify(MOIST)}; ${((performance.now() - t0) / 1000).toFixed(0)} s`);
for (const [name, test, earth] of REGIMES) {
  const inside = Uint8Array.from({ length: C }, (_, i) => (test(i) ? 1 : 0));
  const effects = saved.meanShortwaveCloudEffect ? `, day-mean SWCRE ${f(mean((i) => saved.meanShortwaveCloudEffect[i], inside), 1)} LWCRE ${f(mean((i) => saved.meanLongwaveCloudEffect[i], inside), 1)} W/m2` : '';
  console.log(`\n${name}: total cover ${f(mean((i) => allCover[i], inside))} (maximum-random ${f(mean((i) => randomCover[i], inside))}, exponential-random ${f(mean((i) => exponentialCover[i], inside))})${effects}\n  Earth: ${earth}`);
  for (let n = 0; n < CLASSES.length; n++) {
    const cv = mean((i) => cover[n][i], inside), grid = 1000 * mean((i) => path[n][i], inside);
    const bins = [0, 0, 0];
    let weight = 0;
    for (let i = 0; i < C; i++) { if (!inside[i] || !(cover[n][i] > 0)) continue; const tau = visible[n][i] / cover[n][i], w = area[i] * cover[n][i]; weight += w; bins[tau < 3.6 ? 0 : tau < 23 ? 1 : 2] += w; }
    console.log(`  ${CLASSES[n].padEnd(8)} cover ${f(cv, 3)}  grid ${f(grid, 1).padStart(6)} g/m2  in-cloud ${f(cv > 0 ? grid / cv : NaN, 0).padStart(5)} g/m2  tau ${f(cv > 0 ? mean((i) => visible[n][i], inside) / cv : NaN, 1).padStart(5)}  shares <3.6 / 3.6-23 / >23 ${bins.map((b) => f(b / weight)).join(' / ')}`);
  }
  console.log(`  humidity: RH 150-350 hPa ${f(ratio((i) => humidity.upper[0][i], (i) => humidity.upper[2][i], inside))} (over ice ${f(ratio((i) => humidity.upper[1][i], (i) => humidity.upper[2][i], inside))}, layer area above ice saturation ${f(ratio((i) => issr[i], (i) => upperArea[i], inside))}, of it without cloud ${f(ratio((i) => humidUpperDry[i], (i) => humidUpper[i], inside))}); RH 350-700 hPa ${f(ratio((i) => humidity.middle[0][i], (i) => humidity.middle[1][i], inside))}; below 700 hPa ${f(ratio((i) => humidity.lower[0][i], (i) => humidity.lower[1][i], inside))} (layers above 0.9 without cloud ${f(ratio((i) => humidLowerDry[i], (i) => humidLower[i], inside))})`);
  console.log(`  cumulus updraught ${f(mean((i) => updraught[i], inside), 3)}; resolved liquid ${f(1000 * mean((i) => liquidPath[i], inside), 1)} and ice ${f(1000 * mean((i) => icePath[i], inside), 1)} g/m2 by the phase ramp; high cloud colder than 235 K ${f(1000 * mean((i) => coldHigh[i], inside), 1)} g/m2`);
}
