// The longwave gas optics of js/physics/longwave.module.js: a simple spectral
// model fitted to RRTMG over the standard atmospheres, then reduced to the
// model's g-points.
//   node scripts/longwaveFit.mjs            report the spectral model and the g-points
//   FIT=3000 node scripts/longwaveFit.mjs   refit the spectral model (Nelder-Mead iterations)
//   WRITE=1 node scripts/longwaveFit.mjs    rewrite js/physics/longwaveTable.module.js
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
// e^(+-spreadC) (CO2) times the envelope. Fluxes are two-stream with the
// diffusivity 1.66 from 10 to 3250 cm-1, at STEP (5) cm-1 where the
// spectral model is scored by band.
//
// The fit (Nelder-Mead in the logarithms of the coefficients) minimises the
// misses of the g-points against RRTMG for the tropical, midlatitude summer
// and winter and subarctic winter atmospheres (OLR, surface downward flux,
// net flux at 200 hPa, the cooling-rate profile) and against LBLRTM's forcing
// for doubled CO2 on the midlatitude summer profile (Iacono et al. 2008),
// with the spectral model's own misses in the midlatitude summer's sixteen
// RRTMG bands at the top and the surface and in LBLRTM's doubled-CO2 forcing
// by band (Mlawer et al. 1997, Table 6), and the change of OLR and surface
// downward flux without methane or nitrous oxide (Chou et al. 2001, Table 16).
//
// The reduction: the 1 cm-1 sub-intervals are binned into g-points by the
// midlatitude summer column's optical depth (binKey); each g-point keeps the
// mean coefficients of its sub-intervals, scaled to keep their mean
// transmission where the bin's mean optical depth is 1, and the share of the
// Planck emission they hold, fitted as a quartic in temperature.
import { writeFileSync } from 'node:fs';
import { GRAVITY } from '../js/dynamics/sigmaCore.module.js';
import { BENCHMARK, modelColumn, referenceAt, referenceHeating, layerHeating, interfaceAt, MOLAR } from './standardAtmospheres.mjs';
import { LONGWAVE_CONSTANTS, gasPaths, clearLongwave, normalizedPoints } from '../js/physics/longwave.module.js';

const H = 6.62607015e-34, CL = 2.99792458e8, KB = 1.380649e-23, SIGMA = 5.670374419e-8;
const { pRef, diffusivity } = LONGWAVE_CONSTANTS;
export const SPECTRAL_MODEL = {
  kRot: 122.103, lRot: 47.2827, kVr: 4.87534, lVr1: 34.1124, lVr2: 95.9406, kCo2: 310.468, lCo2: 10.2552, lCo2Hi: 10.5568, kLaser: 8.97986e-4,
  spread: 1.75362, spreadC: 0.96586, roberts: 0.972618, tSelf: 2948.68, tCo2: 0.0307253, kO3: 3394.08, lO3: 5.44702, nO3: 0.168415, kCh4: 410.363, kN2o: 450.778,
};
// The g-points' membership is fixed by this parameter set (a fit's starting
// point), so that refitting SPECTRAL_MODEL moves the coefficients of fixed
// g-points.
export const BINNING_MODEL = {
  kRot: 96.53, lRot: 48.03, kVr: 7.937, lVr1: 37.02, lVr2: 65.21, kCo2: 472.6, lCo2: 10.54, lCo2Hi: 10.99, kLaser: 6.691e-4,
  spread: 1.585, spreadC: 1.707, roberts: 0.9778, tSelf: 3219, tCo2: 0.0211, kO3: 2337, lO3: 6.281, nO3: 0.5772, kCh4: 20, kN2o: 20,
};
const LINEAR = new Set(['tCo2', 'nO3', 'spread', 'spreadC']);

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
  return [1, -1].map((sign) => ({ nu, width: width / 2, line: line * Math.exp(sign * P.spread), continuum, co2: co2 * Math.exp(sign * P.spreadC), o3, ch4, n2o }));
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
    for (let k = 0; k < K; k++) {
      const tau = diffusivity * (s.line * paths.line[k] + s.continuum * paths.continuum[k] + s.co2 * paths.co2[k] + s.o3 * paths.o3[k] + s.ch4 * paths.ch4[k] + s.n2o * paths.n2o[k]);
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
const columns = Object.fromEntries(Object.keys(BENCHMARK.atmospheres).map((a) => [a, modelColumn(BENCHMARK.atmospheres[a])]));
const mls = columns.MLS;
const scaled = (column, gas, factor) => ({ ...column, [gas]: column[gas].map((x) => x * factor) });
const withGas = (column, gas, vmr) => ({ ...column, [gas]: column[gas].map((_, k) => vmr * MOLAR[gas] / MOLAR.air * (1 - column.q[k])) });
const rtmip = (co2) => withGas(withGas(withGas(mls, 'co2', co2), 'ch4', 806e-9), 'n2o', 275e-9);

// Misses of a flux function (column -> {up, down, net}) against the references.
export function misses(fluxes) {
  const out = { atmospheres: {}, bands: [], doubling: {} };
  for (const a of RRTMG) {
    const column = columns[a], ref = BENCHMARK.rrtmgLongwave[a].levels, r = fluxes(column);
    const K = column.T.length, top = ref[ref.length - 1], sfc = ref[0];
    const heat = layerHeating(column, r.net), refHeat = referenceHeating(column, ref);
    let trop = 0, nt = 0, strat = 0, ns = 0;
    for (let k = 0; k < K; k++) {
      const p = 0.5 * (column.levels[k] + column.levels[k + 1]) * column.ps, d = (heat[k] - refHeat[k]) ** 2;
      if (p > 20000) { trop += d; nt++; } else if (p > 300) { strat += d; ns++; }
    }
    out.atmospheres[a] = { olr: r.up[0] - top.up, dlr: r.down[K] - sfc.down, net200: interfaceAt(column, r.net, 20000) - referenceAt(ref, 20000), troposphere: Math.sqrt(trop / nt), stratosphere: Math.sqrt(strat / ns), heat, refHeat };
  }
  return out;
}

export function score(P, verbose = false, keys = null) {
  const intervals = spectrum(P, Number(process.env.STEP ?? 5));
  const table = keys ? tableOf(P, reduce(P, keys)) : null;
  const total = (c) => (table ? clearLongwave(c, { table }) : spectralFluxes(P, c, intervals));
  const m = misses(total);
  let s = 0;
  for (const a of RRTMG) { const e = m.atmospheres[a]; s += e.olr ** 2 + e.dlr ** 2 + e.net200 ** 2 + 20 * e.troposphere ** 2 * 22 + 2 * e.stratosphere ** 2 * 10; }
  const K = mls.T.length;
  let bandLine = '';
  for (const b of BENCHMARK.rrtmgLongwave.MLS.bands) {
    const sub = intervals.filter((x) => x.nu >= b.from && x.nu < b.to);
    const r = spectralFluxes(P, mls, sub);
    const eO = r.up[0] - b.levels[b.levels.length - 1].up, eD = r.down[K] - b.levels[0].down;
    s += 0.5 * (eO ** 2 + eD ** 2);
    bandLine += ` ${b.from}:${eO.toFixed(1)}/${eD.toFixed(1)}`;
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
  s += 20 * ((dT - ref.toa) ** 2 + (d2 - ref.p20000) ** 2 + (dS - ref.surface) ** 2);
  if (verbose) {
    for (const a of RRTMG) { const e = m.atmospheres[a]; console.log(`${a.padEnd(4)} OLR ${e.olr.toFixed(2)}  DLR ${e.dlr.toFixed(2)}  net 200 hPa ${e.net200.toFixed(2)}  cooling rms ${e.troposphere.toFixed(3)} (p > 200 hPa) ${e.stratosphere.toFixed(3)} (3-200 hPa) K/day`); }
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
const BIN_WIDTH = Number(process.env.BIN ?? 1.5), OPAQUE = Number(process.env.OPAQUE ?? 1e4), CLEAR = Number(process.env.CLEAR ?? 0.02);
function binKey(s, columnPaths) {
  const parts = ['line', 'continuum', 'co2', 'o3', 'ch4', 'n2o'].map((n) => s[n] * columnPaths[n]);
  const total = parts.reduce((a, b) => a + b, 0), kind = ['h2o', 'cont', 'co2', 'o3', 'minor', 'minor'][parts.indexOf(Math.max(...parts))];
  if (total > OPAQUE) return `${kind}:opaque`;
  if (total < CLEAR) return `${kind}:clear`;
  return `${kind}:${Math.floor(Math.log(total / CLEAR) / BIN_WIDTH)}`;
}

export function binning(P, step = 1) {
  const paths = gasPaths(mls, P);
  const columnPaths = Object.fromEntries(Object.entries(paths).map(([n, v]) => [n, diffusivity * v.reduce((a, b) => a + b, 0)]));
  return spectrum(P, step).map((s) => binKey(s, columnPaths));
}

const shareCache = new Map();
export function reduce(P, keys = binning(P), step = 1) {
  const bins = new Map(), paths = gasPaths(mls, P);
  const columnPaths = Object.fromEntries(Object.entries(paths).map(([n, v]) => [n, diffusivity * v.reduce((a, b) => a + b, 0)]));
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
    if (!shareCache.has(cacheKey)) shareCache.set(cacheKey, fitQuartic(temperatures, temperatures.map((T) => members.reduce((a, s) => a + planckDensity(s.nu, T) * s.width, 0) / (SIGMA * T ** 4))));
    const tau = (s) => ['line', 'continuum', 'co2', 'o3', 'ch4', 'n2o'].reduce((a, n) => a + s[n] * columnPaths[n], 0), meanTau = mean(tau);
    const scale = meanTau > 0 ? -Math.log(mean((s) => Math.exp(-tau(s) / meanTau))) : 1;
    const g = { key, nu: mean((s) => s.nu), planck: shareCache.get(cacheKey) };
    for (const n of ['line', 'continuum', 'co2', 'o3', 'ch4', 'n2o']) g[n] = scale * mean((s) => s[n]);
    points.push(g);
  }
  return points.sort((a, b) => a.nu - b.nu);
}

export const tableOf = (P, points) => ({ ...P, points: normalizedPoints(points.map((g) => [g.line, g.continuum, g.co2, g.o3, g.ch4, g.n2o, ...g.planck])) });

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
  let P = { ...SPECTRAL_MODEL, ...JSON.parse(process.env.P ?? '{}') };
  const keys = binning(BINNING_MODEL);
  if (process.env.FIT) {
    const names = (process.env.KEYS ?? 'kRot,lRot,kVr,kCo2,lCo2,lCo2Hi,spread,spreadC,kO3,nO3,kCh4,kN2o,tCo2,kLaser,roberts,tSelf').split(',');
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
  const table = { ...P, points: points.map((g) => [g.line, g.continuum, g.co2, g.o3, g.ch4, g.n2o, ...g.planck].map((v) => +v.toPrecision(6))) };
  for (const a of Object.keys(BENCHMARK.iae)) {
    if (a === 'note') continue;
    const c = columns[a], g = clearLongwave(c, { table: tableOf(P, points) }), [, down, olr] = BENCHMARK.iae[a];
    console.log(`${a.padEnd(8)} g-points against ICRCCM line-by-line (Feigelson et al. 1991; CO2 300 ppmv, no CH4 or N2O): OLR ${(g.up[0] - olr).toFixed(2)} DLR ${(g.down[c.T.length] - down).toFixed(2)}`);
  }
  if (process.env.WRITE) {
    const lines = table.points.map((p) => `  [${p.join(', ')}],`).join('\n');
    writeFileSync(new URL('../js/physics/longwaveTable.module.js', import.meta.url), `// Generated by scripts/longwaveFit.mjs; each row a g-point: mass absorption coefficients (m2/kg) at\n// p_ref of the water vapour lines, the self continuum (per Pa of vapour pressure), CO2, ozone, methane\n// and nitrous oxide, then the quartic in (T - 250)/100 of its share of sigma T^4.\nexport const LONGWAVE_SPECTRAL_MODEL = ${JSON.stringify(Object.fromEntries(Object.entries(P).map(([k, v]) => [k, +v.toPrecision(6)])))};\nexport const LONGWAVE_POINTS = [\n${lines}\n];\n`);
  }
}
