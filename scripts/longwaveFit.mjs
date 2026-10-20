// The longwave gas optics of js/physics/longwave.module.js: a simple spectral
// model fitted to RRTMG over the standard atmospheres, then reduced to the
// model's g-points.
//   node scripts/longwaveFit.mjs            report the spectral model and the g-points
//   FIT=3000 node scripts/longwaveFit.mjs   refit the spectral model (Nelder-Mead iterations)
//   WRITE=1 node scripts/longwaveFit.mjs    rewrite js/physics/longwaveTable.module.js
//   LEVELS=bl36 UPPER_WEIGHT=100 FIT=4000 WRITE=1 node scripts/longwaveFit.mjs   the bl36 table
//   LEVELS=bl36 UPPER_WEIGHT=100 PARTITION_WEIGHT=30 FIT=4000 WRITE=bl36ByGas node scripts/longwaveFit.mjs
//                                           the bl36 table refitted gas by gas (longwaveLines 'scaledByGas')
// LEVELS (one of SIGMA_GRIDS, bl34 by default) names the levels the columns are laid on, the
// spectral model the report starts from (longwaveTableFor's) and the table WRITE rewrites; WRITE
// set to a key of LONGWAVE_TABLE_FILES writes that table instead.
//
// The spectral model is that of Jeevanjee & Fueglistaler (2020, JAS 77, 479)
// and Williams et al. (2025, arXiv 2508.09353): the water vapour lines' mass
// absorption coefficient at p_ref = 500 hPa is kRot below 200 cm-1, falls as
// exp(-(nu - 200)/lRot) above it, rises toward the vibration-rotation band as
// exp(-(1450 - nu)/lVr1), is kVr over 1450-1700 cm-1 and falls as
// exp(-(nu - 1700)/lVr2) beyond; CO2's is kCo2 exp(-|nu - 667.5|/l) over
// 500-850 cm-1, l = lCo2 below the band centre and lCo2Hi above it, with
// kLaser over 900-1100 cm-1 for the 9.4 and 10.4 um bands and kCo2 over
// 2200-2400 cm-1; every line coefficient scales as p/p_ref. To these are
// added, as named below: the self continuum of water vapour with the
// spectral shape of Roberts et al. (1976), C(nu) = 1.25e-22 + 1.67e-19
// exp(-7.77e-3 nu) cm2 molec-1 atm-1 at 296 K, times `roberts`, scaled by
// the vapour pressure and exp(tSelf (1/T - 1/296)); CO2's temperature
// dependence exp(tCo2 (T - 250)); ozone's 9.6 um band kO3 exp(-|nu - 1042|/lO3)
// over 980-1100 cm-1 scaling as (p/p_ref)^nO3; methane at 1306 and nitrous
// oxide at 1285 and 589 cm-1. The lines' structure within each spectral
// interval is two equal halves at e^(+-spread) (water vapour) and
// e^(+-spreadC) (CO2) times the envelope, and, for each absorber with a
// tail share, its line cores (TAILS below; the vapour's below tailH2oTo
// cm-1, CO2's on an envelope of e-folding lTailCo2). The vapour lines', CO2's
// and ozone's paths scale with the Voigt pressure of
// js/physics/longwave.module.js. Fluxes are two-stream with the
// diffusivity 1.66 from 10 to 3250 cm-1, at STEP (5) cm-1 where the
// spectral model is scored by band.
//
// The fit (Nelder-Mead in the logarithms of the coefficients) minimises the
// misses of the g-points against RRTMG for the tropical, midlatitude summer
// and winter and subarctic winter atmospheres (OLR, surface downward flux,
// net flux at 200 hPa, the cooling-rate profile, with UPPER_WEIGHT on the
// 3-30 hPa layers and the relative miss of every layer above 3 hPa whose mass
// the reference covers to 90 %: the top layer on bl34, three on bl36) and against LBLRTM's
// forcing for doubled CO2 (DOUBLING_WEIGHT) and for methane and nitrous
// oxide from none to their 1860 amounts (MINOR_WEIGHT) on the midlatitude
// summer profile (Iacono et al. 2008), with PARTITION_WEIGHT on each of CO2's,
// ozone's and the vapor's cooling there (the change without it) against
// RRTMG's bands where it absorbs, 0.3-200 hPa (above 3 hPa as a share of the
// layer's cooling), with the spectral model's own misses
// in the midlatitude summer's sixteen RRTMG bands at the top and the surface
// and in each band's cooling above 200 hPa (BAND_WEIGHT), in LBLRTM's
// doubled-CO2 forcing by band (Mlawer et al. 1997, Table 6), and the change
// of OLR and surface downward flux without methane or nitrous oxide (Chou et
// al. 2001, Table 16).
//
// The reduction: the 1 cm-1 sub-intervals are binned into g-points by the
// midlatitude summer column's optical depth (binKey; above THICK in wider
// bins of THICK_BIN, so that the cores that reach the top layer keep their
// own g-points up to OPAQUE); each g-point keeps the
// mean coefficients of its sub-intervals, scaled to keep their mean
// transmission where the bin's mean optical depth is 1, and the share of the
// Planck emission they hold, fitted as a quartic in temperature.
import { writeFileSync } from 'node:fs';
import { GRAVITY, sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { BENCHMARK, modelColumn, referenceAt, referenceHeating, layerHeating, interfaceAt, MOLAR } from './standardAtmospheres.mjs';
import { LONGWAVE_CONSTANTS, gasPaths, clearLongwave, normalizedPoints, longwaveTableFor, LONGWAVE_TABLE_FILES, LONGWAVE_ABSORBERS, LONGWAVE_CORES, CORE_EDGES, CORE_GASES, longwaveCores } from '../js/physics/longwave.module.js';

const H = 6.62607015e-34, CL = 2.99792458e8, KB = 1.380649e-23, SIGMA = 5.670374419e-8;
const { pRef, diffusivity } = LONGWAVE_CONSTANTS;
export const SPECTRAL_MODEL = {
  kRot: 28.1946, lRot: 54.8622, kVr: 4.503, lVr1: 34.1124, lVr2: 95.9406, kCo2: 165.667, lCo2: 8.86948, lCo2Hi: 9.98274, kLaser: 0.00118713,
  spread: 1.76148, spreadC: 0.511931, roberts: 0.961601, tSelf: 2417.52, tCo2: 0.00462344, kO3: 479.29, lO3: 5.44702, nO3: 0.936731, kCh4: 587.647,
  kN2o: 735.849, tailCo2: 0.15201, depthCo2: 2.36796, lTailCo2: 16.4622, tailH2o: 0.0122694, depthH2o: 8.18674, dopplerO3: 1012, tailO3: 0.504104,
  depthO3: 16.8623, dopplerCo2: 727, dopplerH2o: 400, tailH2oTo: 500,
};
// The g-points' membership is fixed by this parameter set (a fit's starting
// point), so that refitting SPECTRAL_MODEL moves the coefficients of fixed
// g-points.
export const BINNING_MODEL = {
  kRot: 96.53, lRot: 48.03, kVr: 7.937, lVr1: 37.02, lVr2: 65.21, kCo2: 472.6, lCo2: 10.54, lCo2Hi: 10.99, kLaser: 6.691e-4,
  spread: 1.585, spreadC: 1.707, roberts: 0.9778, tSelf: 3219, tCo2: 0.0211, kO3: 2337, lO3: 6.281, nO3: 0.5772, kCh4: 20, kN2o: 20,
};
// The line-core tails the g-points' membership is fixed by, for a model that has them.
export const BINNING_TAILS = { tailH2o: 0.05, depthH2o: 6, tailCo2: 0.05, depthCo2: 6, tailO3: 0.05, depthO3: 6 };
// The Doppler cores the g-points' membership is fixed by, for a model that has them.
export const BINNING_CORES = { coreH2o: 1, coreCo2: 3, coreO3: 2, spacingH2o: 5, spacingCo2: 1.5, spacingO3: 0.1, dopplerCoreH2o: 400, dopplerCoreCo2: 727, dopplerCoreO3: 1012 };
export const binningModel = (P) => ({ ...BINNING_MODEL, ...Object.fromEntries(Object.entries(BINNING_TAILS).filter(([key]) => P[key.startsWith('tail') ? key : `tail${key.slice(5)}`] > 0)), ...(P.tailH2oTo ? { tailH2oTo: P.tailH2oTo } : {}), ...(longwaveCores(P) ? BINNING_CORES : {}) });
const LINEAR = new Set(['tCo2', 'nO3', 'spread', 'spreadC']);

// Each absorber's line cores: the share g0 of every interval where it absorbs
// holds the high-k tail of a Lorentz line's k-distribution, k(g) = k0 (g0/g)^2
// from the upper half's k0 at g = g0 down to g0 exp(-depth), in TAIL_NODES
// nodes of equal width in ln g, each with its mean k.
const TAILS = [['line', 'tailH2o', 'depthH2o'], ['co2', 'tailCo2', 'depthCo2'], ['o3', 'tailO3', 'depthO3']];
const TAIL_NODES = Number(process.env.TAIL_NODES ?? 4);

// A model with Doppler cores (longwaveCores) gives each absorber's lines a
// Voigt shape. Their far wings, the two halves and the tails above, scale
// with the Lorentz width p/p_ref at every pressure (its dopplerH2o,
// dopplerCo2 and dopplerO3 are 0). Each line's center holds its strength S in
// a Gaussian of the Doppler 1/e half-width alpha_D = nu sqrt(2 k T/m)/c at
// 250 K: lines spaced spacingH2o, spacingCo2 and spacingO3 (cm-1) apart, of
// mean strength S/delta = coreH2o, coreCo2 and coreO3 times the absorber's
// envelope at p_ref (the vapor's lines, CO2's 500-850 cm-1 band, ozone's
// 9.6 um band), in nodes of |nu - nu_line|/alpha_D between CORE_EDGES, each
// with the Gaussian's mean. Each node's share of the interval is
// 2 alpha_D dx/delta (all of an absorber's together at most CORE_FILL, its
// strength kept), and its coefficient that mean, which coreFactor carries to
// the Voigt function's at every pressure.
const CORES = [['line', 'lineCore', 'coreH2o', 'spacingH2o', 1.6024e-6], ['co2', 'co2Core', 'coreCo2', 'spacingCo2', 1.0251e-6], ['o3', 'o3Core', 'coreO3', 'spacingO3', 9.817e-7]];
const CORE_FILL = 0.25;
const erf = (x) => { const t = 1 / (1 + 0.3275911 * x); return 1 - ((((1.061405429 * t - 1.453152027) * t + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x); };
function dopplerCores(P, nu, envelope) {
  if (!longwaveCores(P)) return [];
  const out = [];
  for (const [gas, name, strength, spacing, doppler] of CORES) {
    if (!(envelope[gas] > 0 && P[strength] > 0)) continue;
    const alpha = doppler * nu, fill = Math.min(1, CORE_FILL * P[spacing] / (2 * CORE_EDGES[CORE_EDGES.length - 1] * alpha));
    for (let j = 1; j < CORE_EDGES.length; j++) {
      const share = 2 * alpha * (CORE_EDGES[j] - CORE_EDGES[j - 1]) / P[spacing] * fill;
      out.push({ name, node: `${name}${j}`, share, k: P[strength] * envelope[gas] * (erf(CORE_EDGES[j]) - erf(CORE_EDGES[j - 1])) / share });
    }
  }
  return out;
}

export function planckDensity(nu, T) {
  const m = nu * 100;
  return Math.PI * 2 * H * CL * CL * m ** 3 / Math.expm1(H * CL * m / (KB * T)) * 100;
}

// Absorption coefficients at p_ref of one wavenumber's two halves.
export function subIntervals(P, nu, width) {
  let line;
  if (nu < 200) line = P.kRot;
  else if (nu < 1450) line = Math.max(P.kRot * Math.exp(-(nu - 200) / P.lRot), P.kVr * Math.exp(-(1450 - nu) / P.lVr1));
  else if (nu < 1700) line = P.kVr;
  else line = P.kVr * Math.exp(-(nu - 1700) / P.lVr2);
  const continuum = (1.25e-22 + 1.67e-19 * Math.exp(-7.77e-3 * nu)) * 1e-4 * 6.02214076e23 / (MOLAR.h2o * 1e-3) / 101325 * P.roberts;
  const co2 = nu >= 500 && nu <= 850 ? P.kCo2 * Math.exp(-Math.abs(nu - 667.5) / (nu < 667.5 ? P.lCo2 : P.lCo2Hi)) : nu >= 900 && nu <= 1100 ? P.kLaser : nu >= 2200 && nu <= 2400 ? P.kCo2 : 0;
  const o3 = nu >= 980 && nu <= 1100 ? P.kO3 * Math.exp(-Math.abs(nu - 1042) / P.lO3) : 0;
  const ch4 = nu >= 1200 && nu <= 1400 ? P.kCh4 * Math.exp(-Math.abs(nu - 1306) / 25) : 0;
  const n2o = (nu >= 1200 && nu <= 1350 ? P.kN2o * Math.exp(-Math.abs(nu - 1285) / 20) : 0) + (nu >= 550 && nu <= 630 ? 0.1 * P.kN2o * Math.exp(-Math.abs(nu - 589) / 15) : 0);
  const bulk = { line: line * Math.cosh(P.spread), co2: co2 * Math.cosh(P.spreadC), o3 };
  const tails = TAILS.filter(([gas, share]) => P[share] > 0 && bulk[gas] > 0 && !(gas === 'line' && nu > (P.tailH2oTo ?? Infinity)));
  const cores = dopplerCores(P, nu, { line, co2: nu >= 500 && nu <= 850 ? co2 : 0, o3 });
  const kept = 1 - tails.reduce((s, [, share, depth]) => s + P[share] * -Math.expm1(-P[depth]), 0) - cores.reduce((s, c) => s + c.share, 0);
  const out = [1, -1].map((sign) => ({ nu, width: width * kept / 2, line: line * Math.exp(sign * P.spread), continuum, co2: co2 * Math.exp(sign * P.spreadC), o3, ch4, n2o }));
  const co2Core = P.lTailCo2 > 0 && nu >= 500 && nu <= 850 ? P.kCo2 * Math.exp(-Math.abs(nu - 667.5) / P.lTailCo2) : co2;
  for (const [gas, share, depth] of tails) {
    const g0 = P[share], top = gas === 'line' ? line * Math.exp(P.spread) : gas === 'co2' ? co2Core * Math.exp(P.spreadC) : o3;
    for (let j = 1; j <= TAIL_NODES; j++) {
      const b = g0 * Math.exp(-P[depth] * (j - 1) / TAIL_NODES), a = g0 * Math.exp(-P[depth] * j / TAIL_NODES);
      out.push({ nu, width: width * (b - a), ...bulk, continuum, ch4, n2o, [gas]: top * g0 * g0 / (a * b) });
    }
  }
  for (const c of cores) out.push({ nu, width: width * c.share, ...bulk, continuum, ch4, n2o, [c.node]: c.k, core: c.name, node: c.node });
  return out;
}

export function spectrum(P, step = 1, from = 10, to = 3250) {
  const out = [];
  for (let nu = from + step / 2; nu < to; nu += step) out.push(...subIntervals(P, nu, step));
  return out;
}

// Two-stream fluxes of a set of sub-intervals (each with its own Planck density) on a column.
export function spectralFluxes(P, column, intervals) {
  const paths = gasPaths(column, P);
  const K = column.T.length, up = new Float64Array(K + 1), down = new Float64Array(K + 1);
  const tr = new Float64Array(K), B = new Float64Array(K);
  for (const s of intervals) {
    const node = s.node ? paths[s.node] : null, kNode = s.node ? s[s.node] : 0;
    for (let k = 0; k < K; k++) {
      let tau = 0;
      for (const n of LONGWAVE_ABSORBERS) tau += s[n] * paths[n][k];
      if (node) tau += kNode * node[k];
      tau *= diffusivity;
      tr[k] = Math.exp(-tau);
      B[k] = planckDensity(s.nu, column.T[k]) * s.width;
    }
    let d = 0;
    for (let k = 0; k < K; k++) { d = d * tr[k] + B[k] * (1 - tr[k]); down[k + 1] += d; }
    let u = planckDensity(s.nu, column.Ts) * s.width;
    up[K] += u;
    for (let k = K - 1; k >= 0; k--) { u = u * tr[k] + B[k] * (1 - tr[k]); up[k] += u; }
  }
  return { up, down, net: Array.from(up, (x, k) => x - down[k]) };
}

const RRTMG = ['TROP', 'MLS', 'MLW', 'SAW'];
const UPPER_WEIGHT = Number(process.env.UPPER_WEIGHT ?? 30), OLR_WEIGHT = Number(process.env.OLR_WEIGHT ?? 1), BAND_WEIGHT = Number(process.env.BAND_WEIGHT ?? 3), MINOR_WEIGHT = Number(process.env.MINOR_WEIGHT ?? 20), DOUBLING_WEIGHT = Number(process.env.DOUBLING_WEIGHT ?? 150), DLR_WEIGHT = Number(process.env.DLR_WEIGHT ?? 3), TROPOSPHERE_WEIGHT = Number(process.env.TROPOSPHERE_WEIGHT ?? 1000);
const LEVELS = process.env.LEVELS ?? 'bl34', LINES = process.env.LINE_SHAPE ?? 'scaled', TOP_WEIGHT = Number(process.env.TOP_WEIGHT ?? 0);
// With FAIR the RRTMG atmospheres' columns end where the reference's does
// (6.7 Pa in the midlatitude summer, 10 Pa in the subarctic winter), not at 0.
const FAIR = process.env.FAIR ? process.env.FAIR === '1' : LINES === 'voigt';
const columns = Object.fromEntries(Object.keys(BENCHMARK.atmospheres).map((a) => {
  const column = modelColumn(BENCHMARK.atmospheres[a], sigmaInterfaces(LEVELS)), ref = BENCHMARK.rrtmgLongwave[a]?.levels;
  if (!FAIR || !ref) return [a, column];
  const levels = Float64Array.from(column.levels);
  levels[0] = ref[ref.length - 1].p / column.ps;
  return [a, { ...column, levels }];
}));
const mls = columns.MLS;
const scaled = (column, gas, factor) => ({ ...column, [gas]: column[gas].map((x) => x * factor) });
const withGas = (column, gas, vmr) => ({ ...column, [gas]: column[gas].map((_, k) => vmr * MOLAR[gas] / MOLAR.air * (1 - column.q[k])) });
const rtmip = (co2) => withGas(withGas(withGas(mls, 'co2', co2), 'ch4', 806e-9), 'n2o', 275e-9);
const PARTITION_WEIGHT = Number(process.env.PARTITION_WEIGHT ?? 0), PARTITION_BANDS = { co2: [500, 630, 700, 2250, 2380], o3: [980, 1080], q: [10, 350, 1390, 1480, 1800] };
const partition = Object.entries(PARTITION_BANDS).map(([gas, from]) => [gas, BENCHMARK.rrtmgLongwave.MLS.bands.filter((b) => from.includes(b.from)).map((b) => referenceHeating(mls, b.levels)).reduce((a, h) => a.map((v, k) => v + h[k]))]);

const linearAt = (ref, p) => {
  const j = ref.findIndex((l, i) => i < ref.length - 1 && p <= l.p && p >= ref[i + 1].p);
  if (j < 0) return p > ref[0].p ? ref[0].net : ref[ref.length - 1].net;
  return ref[j].net + (ref[j + 1].net - ref[j].net) * (ref[j].p - p) / (ref[j].p - ref[j + 1].p);
};

// Misses of a flux function (column -> {up, down, net}) against the references.
// A layer above 3 hPa counts only where RRTMG's levels resolve it: its heating
// from the reference fluxes interpolated in p and in ln p agree to 0.25 K/day
// (not so at 0.65 and 1.6 hPa in the midlatitude summer, which has no level
// between 1.3 and 0.59 hPa).
export function misses(fluxes) {
  const out = { atmospheres: {}, bands: [], doubling: {} };
  for (const a of RRTMG) {
    const column = columns[a], ref = BENCHMARK.rrtmgLongwave[a].levels, r = fluxes(column);
    const K = column.T.length, top = ref[ref.length - 1], sfc = ref[0];
    const heat = layerHeating(column, r.net), refHeat = referenceHeating(column, ref);
    const linear = layerHeating(column, Array.from(column.levels, (s) => linearAt(ref, Math.max(s * column.ps, top.p))));
    let trop = 0, nt = 0, strat = 0, ns = 0, upper = 0, nu = 0, tops = 0, topMiss = 0;
    for (let k = 0; k < K; k++) {
      const p = 0.5 * (column.levels[k] + column.levels[k + 1]) * column.ps, d = (heat[k] - refHeat[k]) ** 2;
      if (p > 20000) { trop += d; nt++; } else if (p > 300) { strat += d; ns++; }
      if (p > 300 && p < 3000) { upper += d; nu++; }
      const covered = (column.levels[k + 1] * column.ps - Math.max(column.levels[k] * column.ps, top.p)) / ((column.levels[k + 1] - column.levels[k]) * column.ps);
      if (p < 300 && covered >= 0.9 && Math.abs(refHeat[k] - linear[k]) < 0.25) { tops += (heat[k] / refHeat[k] - 1) ** 2; topMiss += (heat[k] - refHeat[k]) ** 2; }
    }
    out.atmospheres[a] = { olr: r.up[0] - top.up, dlr: r.down[K] - sfc.down, net200: interfaceAt(column, r.net, 20000) - referenceAt(ref, 20000), troposphere: Math.sqrt(trop / nt), stratosphere: Math.sqrt(strat / ns), upper: Math.sqrt(upper / nu), top: heat[0] / refHeat[0] - 1, tops, topMiss, heat, refHeat };
  }
  return out;
}

export function score(P, verbose = false, keys = null) {
  const intervals = spectrum(P, Number(process.env.STEP ?? 5));
  const table = keys ? tableOf(P, reduce(P, keys)) : null;
  const total = (c) => (table ? clearLongwave(c, { table }) : spectralFluxes(P, c, intervals));
  const m = misses(total);
  let s = 0;
  for (const a of RRTMG) { const e = m.atmospheres[a]; s += OLR_WEIGHT * e.olr ** 2 + DLR_WEIGHT * e.dlr ** 2 + e.net200 ** 2 + TROPOSPHERE_WEIGHT * e.troposphere ** 2 + 2 * e.stratosphere ** 2 * 10 + UPPER_WEIGHT * (4 * e.upper ** 2 + e.tops) + TOP_WEIGHT * e.topMiss; }
  const K = mls.T.length;
  let bandLine = '', bandCooling = '';
  for (const b of BENCHMARK.rrtmgLongwave.MLS.bands) {
    const sub = intervals.filter((x) => x.nu >= b.from && x.nu < b.to);
    const r = spectralFluxes(P, mls, sub);
    const eO = r.up[0] - b.levels[b.levels.length - 1].up, eD = r.down[K] - b.levels[0].down;
    s += 0.5 * (eO ** 2 + eD ** 2);
    bandLine += ` ${b.from}:${eO.toFixed(1)}/${eD.toFixed(1)}`;
    const heat = layerHeating(mls, r.net), refHeat = referenceHeating(mls, b.levels);
    let miss = 0, n = 0;
    for (let k = 0; k < K && 0.5 * (mls.levels[k] + mls.levels[k + 1]) * mls.ps < 20000; k++) { miss += (heat[k] - refHeat[k]) ** 2; n++; }
    s += BAND_WEIGHT * miss;
    bandCooling += ` ${b.from}:${Math.sqrt(miss / n).toFixed(2)}`;
  }
  if (verbose) console.log(`MLS by RRTMG band, cooling rms above 200 hPa, K/day:${bandCooling}`);
  if (table && (PARTITION_WEIGHT > 0 || verbose)) {
    const all = layerHeating(mls, total(mls).net), refAll = referenceHeating(mls, BENCHMARK.rrtmgLongwave.MLS.levels);
    let line = '';
    for (const [gas, want] of partition) {
      const without = layerHeating(mls, total(gas === 'q' ? { ...mls, q: mls.q.map(() => 1e-9) } : scaled(mls, gas, 0)).net);
      let miss = 0, n = 0, top = '';
      for (let k = 0; k < K; k++) {
        const p = 0.5 * (mls.levels[k] + mls.levels[k + 1]) * mls.ps;
        if (p < 30 || p > 20000) continue;
        const d = all[k] - without[k] - (p > 300 ? want[k] : want[k] * all[k] / refAll[k]);
        miss += d * d; n++;
        if (p < 300) top += ` ${(all[k] - without[k]).toFixed(2)} (${(want[k] * all[k] / refAll[k]).toFixed(2)})`;
      }
      s += PARTITION_WEIGHT * miss;
      line += ` ${gas === 'q' ? 'h2o' : gas} ${Math.sqrt(miss / n).toFixed(3)}${top}`;
    }
    if (verbose) console.log(`MLS cooling by gas against its RRTMG bands, rms 0.3-200 hPa and the layers above 3 hPa, K/day:${line}`);
  }
  let doublingLine = '';
  const twice = scaled(mls, 'co2', 2);
  for (const [a, b, tT, tP, tS] of BENCHMARK.mlawerMls.doubling.bands) {
    const sub = intervals.filter((x) => x.nu >= a && x.nu < b);
    const g1 = spectralFluxes(P, mls, sub), g2 = spectralFluxes(P, twice, sub);
    const eT = g1.net[0] - g2.net[0], eP = interfaceAt(mls, g1.net, 17900) - interfaceAt(mls, g2.net, 17900), eS = g2.down[K] - g1.down[K];
    s += 20 * ((eT - tT) ** 2 + (eP - tP) ** 2 + (eS - tS) ** 2);
    doublingLine += ` ${a}-${b}: ${eT.toFixed(2)} (${tT}) / ${eP.toFixed(2)} (${tP}) / ${eS.toFixed(2)} (${tS})`;
  }
  let minorLine = '';
  for (const [a, gas, tO, tD] of [['MLS', 'ch4', 2.22, -0.86], ['MLS', 'n2o', 1.83, -0.58], ['SAW', 'ch4', 1.00, -1.21], ['SAW', 'n2o', 1.10, -1.23]]) {
    const c = columns[a], g1 = total(c), g2 = total(scaled(c, gas, 0)), K2 = c.T.length;
    const eO = g2.up[0] - g1.up[0], eD = g2.down[K2] - g1.down[K2];
    s += 5 * ((eO - tO) ** 2 + (eD - tD) ** 2);
    minorLine += ` ${a} ${gas} ${eO.toFixed(2)} (${tO}) / ${eD.toFixed(2)} (${tD})`;
  }
  if (verbose) console.log(`removing a minor gas, OLR / surface down change (Chou et al. 2001, Table 16):${minorLine}`);
  const f1 = total(rtmip(287e-6)), f2 = total(rtmip(574e-6));
  const ref = BENCHMARK.iacono.co2Doubling.longwave;
  const dT = f1.net[0] - f2.net[0], d2 = interfaceAt(mls, f1.net, 20000) - interfaceAt(mls, f2.net, 20000), dS = f2.down[K] - f1.down[K];
  s += DOUBLING_WEIGHT * ((dT - ref.toa) ** 2 + (d2 - ref.p20000) ** 2 + (dS - ref.surface) ** 2);
  for (const [, share] of TAILS) if (P[share] > 0.5) s += 1e3 * (P[share] - 0.5) ** 2;
  const overlap = 1 - (P.tailCo2 ?? 0) * -Math.expm1(-(P.depthCo2 ?? 0)) - (P.tailO3 ?? 0) * -Math.expm1(-(P.depthO3 ?? 0));
  if (overlap < 0.1) s += 1e3 * (0.1 - overlap) ** 2;
  const n1 = total(withGas(withGas(rtmip(287e-6), 'ch4', 0), 'n2o', 0)), refMinor = BENCHMARK.iacono.minorGases.longwave;
  const mT = n1.net[0] - f1.net[0], m2 = interfaceAt(mls, n1.net, 20000) - interfaceAt(mls, f1.net, 20000), mS = f1.down[K] - n1.down[K];
  s += MINOR_WEIGHT * ((mT - refMinor.toa) ** 2 + (m2 - refMinor.p20000) ** 2 + (mS - refMinor.surface) ** 2);
  if (verbose) {
    console.log(`CH4 0 -> 806 ppbv and N2O 0 -> 275 ppbv: TOA ${mT.toFixed(2)} (${refMinor.toa}) 200 hPa ${m2.toFixed(2)} (${refMinor.p20000}) surface ${mS.toFixed(2)} (${refMinor.surface})`);
    for (const a of RRTMG) {
      const e = m.atmospheres[a], c = columns[a], above = [];
      for (let k = 0; 0.5 * (c.levels[k] + c.levels[k + 1]) * c.ps < 300; k++) above.push(`${e.heat[k].toFixed(2)} (${e.refHeat[k].toFixed(2)})`);
      console.log(`${a.padEnd(4)} OLR ${e.olr.toFixed(2)}  DLR ${e.dlr.toFixed(2)}  net 200 hPa ${e.net200.toFixed(2)}  cooling rms ${e.troposphere.toFixed(3)} (p > 200 hPa) ${e.stratosphere.toFixed(3)} (3-200 hPa) ${e.upper.toFixed(3)} (3-30 hPa) K/day; layers above 3 hPa ${above.join(' / ')}`);
    }
    console.log(`MLS by RRTMG band, OLR/DLR misses:${bandLine}`);
    console.log(`doubled CO2 by band, TOA / tropopause / surface (LBLRTM):${doublingLine}`);
    console.log(`doubled CO2 287 -> 574 ppmv: TOA ${dT.toFixed(2)} (${ref.toa}) 200 hPa ${d2.toFixed(2)} (${ref.p20000}) surface ${dS.toFixed(2)} (${ref.surface})`);
  }
  return s;
}

function nelderMead(f, x0, steps, iterations) {
  const n = x0.length;
  let simplex = [x0.slice(), ...x0.map((_, i) => x0.map((x, j) => x + (i === j ? steps[i] : 0)))];
  let values = simplex.map(f);
  for (let it = 0; it < iterations; it++) {
    const order = values.map((_, i) => i).sort((a, b) => values[a] - values[b]);
    simplex = order.map((i) => simplex[i]); values = order.map((i) => values[i]);
    const centre = Array.from({ length: n }, (_, j) => simplex.slice(0, n).reduce((s, x) => s + x[j], 0) / n);
    const at = (t) => centre.map((c, j) => c + t * (simplex[n][j] - c));
    const xr = at(-1), fr = f(xr);
    if (fr < values[0]) { const xe = at(-2), fe = f(xe); [simplex[n], values[n]] = fe < fr ? [xe, fe] : [xr, fr]; }
    else if (fr < values[n - 1]) { simplex[n] = xr; values[n] = fr; }
    else {
      const xc = at(0.5), fc = f(xc);
      if (fc < values[n]) { simplex[n] = xc; values[n] = fc; }
      else for (let i = 1; i <= n; i++) { simplex[i] = simplex[i].map((x, j) => simplex[0][j] + 0.5 * (x - simplex[0][j])); values[i] = f(simplex[i]); }
    }
  }
  return simplex[0];
}

// The midlatitude summer column's optical depth in each absorber decides a
// sub-interval's g-point: the absorber that dominates it, and the logarithm
// of its total in bins of BIN_WIDTH, with the opaque (above OPAQUE) and the
// transparent (below CLEAR) in one g-point per absorber.
const BIN_WIDTH = Number(process.env.BIN ?? 1.5), THICK = Number(process.env.THICK ?? 1e4), THICK_WIDTH = Number(process.env.THICK_BIN ?? 4.5), OPAQUE = Number(process.env.OPAQUE ?? 1e9), CLEAR = Number(process.env.CLEAR ?? 0.02);
// A Doppler core's g-point is decided by its own optical depth above 10 hPa
// in bins of CORE_BIN, those over CORE_OPAQUE together; one under CORE_CLEAR,
// thin where it differs from the far wings, joins its interval's other nodes.
const CORE_BIN = Number(process.env.CORE_BIN ?? 1.6), CORE_CLEAR = Number(process.env.CORE_CLEAR ?? 0.1), CORE_OPAQUE = Number(process.env.CORE_OPAQUE ?? 300), CORE_TOP = 1000;
function binKey(s, columnPaths, topPaths) {
  const coreTau = s.core ? s[s.node] * topPaths[s.node] : 0;
  if (coreTau >= CORE_CLEAR) return `${s.core}:${coreTau > CORE_OPAQUE ? 'opaque' : Math.floor(Math.log(coreTau / CORE_CLEAR) / CORE_BIN)}`;
  const parts = LONGWAVE_ABSORBERS.map((n) => s[n] * columnPaths[n]);
  const total = parts.reduce((a, b) => a + b, 0), kind = ['h2o', 'cont', 'co2', 'o3', 'minor', 'minor'][parts.indexOf(Math.max(...parts))];
  if (total > OPAQUE) return `${kind}:opaque`;
  if (total < CLEAR) return `${kind}:clear`;
  if (total > THICK) return `${kind}:thick${Math.floor(Math.log(total / THICK) / THICK_WIDTH)}`;
  return `${kind}:${Math.floor(Math.log(total / CLEAR) / BIN_WIDTH)}`;
}

function pathSums(P) {
  const paths = gasPaths(mls, P), top = Array.from(mls.levels.slice(1), (s) => s * mls.ps <= CORE_TOP);
  const sum = (v, only) => diffusivity * v.reduce((a, b, k) => a + (!only || top[k] ? b : 0), 0);
  return [Object.fromEntries(Object.entries(paths).map(([n, v]) => [n, sum(v)])), Object.fromEntries(Object.entries(paths).map(([n, v]) => [n, sum(v, true)]))];
}
export function binning(P, step = 1) {
  const [columnPaths, topPaths] = pathSums(P);
  return spectrum(P, step).map((s) => binKey(s, columnPaths, topPaths));
}

const shareCache = new Map();
export function reduce(P, keys = binning(P), step = 1) {
  const bins = new Map(), [columnPaths, topPaths] = pathSums(P), names = Object.keys(columnPaths);
  spectrum(P, step).forEach((s, j) => {
    if (!bins.has(keys[j])) bins.set(keys[j], []);
    bins.get(keys[j]).push(s);
  });
  const temperatures = Array.from({ length: 39 }, (_, j) => 150 + 5 * j);
  const points = [];
  for (const [key, members] of bins) {
    const weight = members.reduce((a, s) => a + s.width, 0);
    const mean = (f) => members.reduce((a, s) => a + s.width * f(s), 0) / weight;
    const cacheKey = `${key}:${members.length}:${members[0].nu}`;
    if (!shareCache.has(cacheKey)) shareCache.set(cacheKey, members.map((s) => temperatures.map((T) => planckDensity(s.nu, T) / (SIGMA * T ** 4))));
    const density = shareCache.get(cacheKey), planck = fitQuartic(temperatures, temperatures.map((_, t) => members.reduce((a, s, j) => a + density[j][t] * s.width, 0)));
    const core = CORE_GASES.find((n) => key.startsWith(`${n}:`)), tau = core ? (s) => s[s.node] * topPaths[s.node] : (s) => LONGWAVE_ABSORBERS.reduce((a, n) => a + s[n] * columnPaths[n], 0), meanTau = mean(tau);
    const scale = meanTau > 0 ? -Math.log(mean((s) => Math.exp(-tau(s) / meanTau))) : 1;
    const g = { key, nu: mean((s) => s.nu), planck };
    for (const n of LONGWAVE_ABSORBERS) g[n] = (!core ? scale : 1) * mean((s) => s[n]);
    if (names.length > LONGWAVE_ABSORBERS.length) {
      for (const n of LONGWAVE_CORES) g[n] = 0;
      for (const s of members) if (s.node) g[s.node] += s.width * s[s.node];
      for (const n of LONGWAVE_CORES) g[n] *= (!core || n.startsWith(core) ? scale : 1) / weight;
    }
    points.push(g);
  }
  return points.sort((a, b) => a.nu - b.nu);
}

const rowOf = (g) => [g.line, g.continuum, g.co2, g.o3, g.ch4, g.n2o, ...g.planck, ...(g.co2Core1 !== undefined ? LONGWAVE_CORES.map((n) => g[n]) : [])];
export const tableOf = (P, points) => ({ ...P, points: normalizedPoints(points.map(rowOf)) });

function fitQuartic(x, y) {
  const t = x.map((T) => (T - 250) / 100), n = 5;
  const A = Array.from({ length: n }, () => new Float64Array(n + 1));
  for (let i = 0; i < t.length; i++) for (let r = 0; r < n; r++) { for (let c = 0; c < n; c++) A[r][c] += t[i] ** (r + c); A[r][n] += t[i] ** r * y[i]; }
  for (let c = 0; c < n; c++) {
    for (let r = c + 1; r < n; r++) { const f = A[r][c] / A[c][c]; for (let j = c; j <= n; j++) A[r][j] -= f * A[c][j]; }
  }
  const out = new Array(n).fill(0);
  for (let r = n - 1; r >= 0; r--) { let s = A[r][n]; for (let c = r + 1; c < n; c++) s -= A[r][c] * out[c]; out[r] = s / A[r][r]; }
  return out;
}

if (import.meta.url === `file://${process.argv[1]}`) {
  const { points: _points, ...levelModel } = longwaveTableFor(sigmaInterfaces(LEVELS), LINES);
  let P = { ...SPECTRAL_MODEL, ...levelModel, ...JSON.parse(process.env.P ?? '{}') };
  const keys = binning(process.env.BIN_BY === 'model' ? P : binningModel(P));
  if (process.env.FIT) {
    const names = (process.env.KEYS ?? `kRot,lRot,kVr,kCo2,lCo2,lCo2Hi,spread,spreadC,kO3,nO3,kCh4,kN2o,tCo2,kLaser,roberts,tSelf,tailCo2,depthCo2,lTailCo2,tailH2o,depthH2o,tailO3,depthO3${longwaveCores(P) ? ',coreH2o,coreCo2,coreO3,spacingH2o,spacingCo2,spacingO3' : ''}`).split(',');
    const toX = (p) => names.map((k) => (LINEAR.has(k) ? p[k] : Math.log(p[k])));
    const fromX = (x) => ({ ...P, ...Object.fromEntries(names.map((k, i) => [k, LINEAR.has(k) ? x[i] : Math.exp(x[i])])) });
    P = fromX(nelderMead((x) => { const v = score(fromX(x), false, keys); return Number.isFinite(v) ? v : 1e12; }, toX(P), names.map((k) => (LINEAR.has(k) ? 0.05 : 0.2)), Number(process.env.FIT)));
    console.log('fitted', JSON.stringify(Object.fromEntries(Object.entries(P).map(([k, v]) => [k, +v.toPrecision(6)]))));
  }
  console.log(`spectral model at ${process.env.STEP ?? 5} cm-1:`);
  score(P, true);
  const points = reduce(P, keys);
  console.log(`${points.length} g-points:`);
  score(P, true, keys);
  const table = { ...P, points: points.map((g) => rowOf(g).map((v) => +v.toPrecision(6))) };
  for (const a of Object.keys(BENCHMARK.iae)) {
    if (a === 'note') continue;
    const c = columns[a], g = clearLongwave(c, { table: tableOf(P, points) }), [, down, olr] = BENCHMARK.iae[a];
    console.log(`${a.padEnd(8)} g-points against ICRCCM line-by-line (Feigelson et al. 1991; CO2 300 ppmv, no CH4 or N2O): OLR ${(g.up[0] - olr).toFixed(2)} DLR ${(g.down[c.T.length] - down).toFixed(2)}`);
  }
  if (process.env.WRITE) {
    const lines = table.points.map((p) => `  [${p.join(', ')}],`).join('\n');
    const file = LONGWAVE_TABLE_FILES[process.env.WRITE] ?? (LINES === 'voigt' ? LONGWAVE_TABLE_FILES[`${LEVELS}Voigt`] : LONGWAVE_TABLE_FILES[LEVELS] ?? LONGWAVE_TABLE_FILES.bl34);
    const cores = longwaveCores(P) ? ',\n// then the means (m2/kg) in the Doppler limit of the vapor\'s, CO2\'s and ozone\'s core nodes' : '';
    writeFileSync(new URL(`../js/physics/${file}`, import.meta.url), `// Generated by scripts/longwaveFit.mjs; each row a g-point: mass absorption coefficients (m2/kg) at\n// p_ref of the water vapour lines, the self continuum (per Pa of vapour pressure), CO2, ozone, methane\n// and nitrous oxide, then the quartic in (T - 250)/100 of its share of sigma T^4${cores}.\nexport const LONGWAVE_SPECTRAL_MODEL = ${JSON.stringify(Object.fromEntries(Object.entries(P).map(([k, v]) => [k, +v.toPrecision(6)])))};\nexport const LONGWAVE_POINTS = [\n${lines}\n];\n`);
  }
}
