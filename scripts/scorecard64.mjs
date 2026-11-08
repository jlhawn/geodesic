#!/usr/bin/env node
/*
 * The N=64 observational scorecard, the per-merge gate's second half (the first is scripts/physicsDigest.mjs). Each
 * row is a 30-day day-mean window ending at a save, scored against its definition-matched reference in
 * data/references/ (README.md there gives each one's definition, period and box; scripts/references.mjs reads them).
 *
 * The gate run starts from the standard start, seventeen64_day0365.bin (the March equinox; sha256 STANDARD_START, kept at
 * ~/.cache/webgcm-scratch/seventeen/out64/ on this Mac: the page's runs/seventeen64_day0365 parts drop its means and SST
 * ring). Under gpu_run (about 6 GPU minutes, level 2):
 *   node scripts/scorecard64.mjs --run START OUT [--base BASE.json] [--json OUT.json] [--any-start] [--identical]
 * refuses any other start without --any-start, copies START into the empty directory OUT as p64_day0365.bin, runs
 * LEVELS=bl48 N=64 TAG=p64 spinup.mjs on this tree to days 395 and 425, and scores the two saves; without the copy
 * spinup.mjs would start fresh and run 395 days. Saves already made are scored with
 *   node scripts/scorecard64.mjs [--json OUT.json] [--base BASE.json] [--log FILE] [--start START] [--commit HASH] [--any-start] [--identical] DIR/p64_day0395.bin DIR/p64_day0425.bin
 * which takes the start from --start or from the save 30 days before the first one beside it. Each state's daily lines
 * (TOA fluxes and the sea-ice cuts) are read from TAG.log beside it, or from --log FILE; a missing log is an error.
 * An older log's sea-ice line cut over every sea cell ("sea cells, 15 %") is not the index's cut and is not scored:
 * its rows read no value, and only calibration falls back on the proxy, the cut of the 30-day mean on the index's
 * grids.
 * The JSON records the start's
 * sha256 and day, the save days, each save's sha256, the saves' options, the scoring tree's HEAD (scoredOn) and the
 * commit that made the saves (commit: with --run the tree's HEAD before the first segment, else --commit, else null,
 * since a save does not record it), and with --run the spin-up's environment and treeChangedDuringRun when the fingerprint
 * of HEAD, the tracked changes and the untracked files under js/, scripts/ and data/ (a symlink by its target, a
 * directory by name) differs before a segment or at
 * scoring. A --base whose start or leads differ is refused without --any-start, which prints what differs; --json may
 * not name the --base file (the same device and inode), and --run takes no other arguments.
 *
 * With --base, a row moves when its two-window mean changed by more than t·sigma, t the two-sided 5 % quantile of
 * Student's t at sigma's degrees of freedom, and is unresolved otherwise; a re-pin cites the moved rows. Unresolved is not
 * unchanged: the count of rows whose two-window mean is bit-identical to the base's is printed beside the moved rows.
 * A row with no number on either side (a log without the window's every day, a base that JSON wrote as null) is named
 * as no value, never identical and never moved. --identical, the check for a merge whose changes are off by default,
 * exits 0 only when every gated row has a value on both sides, the two are bit-identical and the recorded starts
 * match, and says whether the saves are byte-identical to the base's (its own, or a run that reproduced them).
 * The null is a
 * merge that changes a default and has no effect on the row: it perturbs from the first step, so the faster rows
 * (rain, TOA, clouds) have decorrelated within one window, and the change over sigma is t-distributed at sigma's
 * degrees of freedom; each such row then moves 5 % of the time where its March–May spread matches the year's that sigma
 * pools, about 1.5 of the 30 rows a merge if all had decorrelated. Where spring is noisier pooling narrows the
 * threshold. The March–June pairs carry 68 % of the 160–100W stress's sigma², about 3 % likely were spring no noisier
 * than the year (a season-free AR(1) null), so that guard likely moves more often than 5 %; they carry 52 % of the SLP
 * contrast's, about 10 % likely under that null. S2's saturated pairs at the gate's own windows, 1.67 sigma for the
 * stress and 1.44 for the SLP contrast, are single draws (about 9 % and 16 % likely with no season). All of it rests on
 * one start's overlapping spring pairs.
 * A fixed 2 sigma would move 7.5–14 % at these degrees of freedom.
 *
 * Sigma (data/scorecard64-sigma.json) is the pair sigma of that two-window mean: the RMS, over adjacent windows at
 * leads of 180 days or more, of the difference between twin controls (the second perturbed by 1e-6 at its start), from
 * the two starts of the pairs (day 365 in March and day 548 in September), a year each: 12 overlapping pairs, 6 a
 * start. By 180 days the twins have decorrelated (the RMS at shorter leads is a tenth of it at 30 days and half at
 * 60). Its degrees of freedom are those of an uncentered sum of squares of the 12 pair means under a Gaussian AR(1) at
 * the row's lag-1 correlation, (tr C)² / tr(C²), at least 2 for two independent starts: under 5 for the three SSTs,
 * tropical-sea rain, net TOA, the Arctic area, the Antarctic extent, the 0n140w core speed and the 0n170w core depth,
 * where sigma's 95 % interval spans a factor of 0.56 to 3.9 at 2.9 degrees of freedom. Saturated and pooled over the
 * year, sigma costs power: a row reads unresolved up to its threshold (in the current record the cold tongue 2.3 K,
 * the 160–100W stress 23 mN/m², the SLP contrast 0.71 hPa, the 30S–30N SST 0.49 K, ASR 1.4 W/m²), and the seasons
 * that set a pooled sigma need not be the gate's: the June–July pairs give 80 % of India rain's sigma², the melt and
 * freeze pairs (Jun/Jul, Sep/Oct and Oct/Nov) 63 % of the Arctic extent's, and with the second July pair and Nov/Dec
 * 90 %. Twins at the gate's own leads move 0 of 30 rows.
 *   node scripts/scorecard64.mjs --sigma [--from LEAD] [--out PATH] [--cache FILE] A_DIR B_DIR START ...
 * rewrites it from twin pairs; rows the twins cannot give (no gridBoxSkin or daily ice line) take the named proxy.
 *
 * The 40–60N land rows: landT2m, the saves' t2m day means (the 2 m air) against ERA5's 2 m air, is definition-matched;
 * landT, the lowest layer's air against the same, is kept for continuity and is not. A save without t2m means reads
 * no value on landT2m. The twin pairs carry no t2m, so landT2m has no sigma and is ungated: reported with its bias,
 * never moved, and outside --identical, until a calibration from twins with t2m exists.
 *
 * Not definition-matched: landT, and the undercurrent rows, which compare the window-end ocean's velocity, a snapshot,
 * with TAO's annual-mean climatology, not day means of the window's calendar months.
 */
import { readFileSync, writeFileSync, existsSync, readdirSync, mkdirSync, copyFileSync, statSync, lstatSync, readlinkSync } from 'node:fs';
import { spawnSync, execSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { basename, dirname, join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { decodeState } from '../js/stateFile.module.js';
import { stateGeography, createGeography } from '../js/geography.module.js';
import { stateGrid } from '../js/stateGrid.module.js';
import { climatologyAt } from '../js/ocean/climatology.module.js';
import { cellVector } from '../js/dynamics/operators.module.js';
import { indexCells, seaIceIndexCut } from '../js/seaIceIndex.module.js';
import { windowMonths } from './woaSstBias.mjs';
import * as R from './references.mjs';

export const STANDARD_START = '1a3454be3d9c263e846e4cb061f5436b3e2ca8da0ce00f11b65a77f64dca5c44';
export const SIGMA_FILE = fileURLToPath(new URL('../data/scorecard64-sigma.json', import.meta.url));
const KELVIN = 273.15, DEG = 180 / Math.PI;
const wrap = (lon) => ((lon + 540) % 360) - 180;
const avg = (xs) => xs.reduce((a, x) => a + x, 0) / xs.length;

export const BOXES = {
  amazon: ['land', -10, 3, -70, -50], india: ['land', 15, 28, 72, 88], tropLand: ['land', -15, 15, -180, 180],
  tropSea: ['sea', -15, 15, -180, 180], io: ['sea', -10, 10, 50, 100], wp: ['sea', -10, 10, 120, 160], global: ['all', -90, 90, -180, 180],
};

export const ROWS = [
  ['sst30', 'SST 30S–30N sea', '°C'], ['sstCold', 'SST cold tongue 2S–2N 150–90W', '°C'], ['sstWarm', 'SST warm pool 10S–10N 130–160E', '°C'],
  ['ts', 'global Ts, grid-box skin', '°C'], ['landT', '40–60N land air: lowest layer against 2 m (not matched; continuity)', 'K'], ['landT2m', '40–60N land 2 m air, day means', 'K'],
  ...Object.keys(BOXES).map((b) => [`rain_${b}`, `rain ${b}`, 'mm/d']),
  ['stress', 'eastward stress 160–100W within 3°', 'mN/m²'], ['slp', 'SLP 180–160W less 120–100W within 3°', 'hPa'],
  ['asr', 'ASR', 'W/m²'], ['net', 'net TOA', 'W/m²'], ['swcre', 'SWCRE', 'W/m²'], ['lwcre', 'LWCRE', 'W/m²'],
  ['iceNExtent', 'Arctic extent, daily cut', 'Mkm²'], ['iceNArea', 'Arctic area, daily cut', 'Mkm²'],
  ['iceSExtent', 'Antarctic extent, daily cut', 'Mkm²'], ['iceSArea', 'Antarctic area, daily cut', 'Mkm²'],
  ...R.TAO_STATIONS.flatMap((s) => [[`euc_${s}`, `undercurrent core speed ${s}, snapshot (not matched)`, 'm/s'], [`eucDepth_${s}`, `undercurrent core depth ${s}, snapshot (not matched)`, 'm']]),
];

export const monthAt = (t) => new Date(Date.UTC(2001, 2, 20, 12) + (((t % 365) + 365) % 365) * 86400e3).getUTCMonth() + 1;
export const windowDays = (D) => Array.from({ length: 30 }, (_, j) => D - 29 + j);
export const dayWeighted = (D, f) => windowDays(D).reduce((s, d) => { let h = 0; for (let k = 0; k < 24; k++) h += f(monthAt(d - 1 + (k + 0.5) / 24)); return s + h / 24; }, 0) / 30;

export function logDays(text) {
  const days = new Map(), ice = new Map(), seaCellIce = new Map();
  for (const line of text.split('\n')) {
    let m = line.match(/^day (\d+) \([^)]*\): Ts [-\d.]+ °C, ASR ([-\d.]+) \(atmosphere [-\d.]+\) OLR ([-\d.]+) W\/m²/);
    if (m) {
      const g = (re) => { const x = line.match(re); return x ? Number(x[1]) : NaN; };
      days.set(+m[1], { asr: +m[2], olr: +m[3], swcre: g(/SWCRE ([-\d.]+)/), lwcre: g(/LWCRE ([-\d.]+)/) });
    }
    m = line.match(/sea ice on day (\d+)'s mean concentration, sea cells( on the index's grids)?, 15 %: N extent ([\d.]+) area ([\d.]+), S extent ([\d.]+) area ([\d.]+)/);
    if (m) (m[2] ? ice : seaCellIce).set(+m[1], { nExtent: +m[3], nArea: +m[4], sExtent: +m[5], sArea: +m[6] });
  }
  return { days, ice, seaCellIce };
}

export function undercurrentCore(profile) {
  const p = profile.filter((x) => Number.isFinite(x.u));
  if (!p.length) return { speed: NaN, depth: NaN };
  let k = 0;
  for (let j = 1; j < p.length; j++) if (p[j].u > p[k].u) k = j;
  let depth = p[k].z, speed = p[k].u;
  if (k > 0 && k < p.length - 1) {
    const [a, b, c] = [p[k - 1].u, p[k].u, p[k + 1].u], h = p[k + 1].z - p[k].z, curvature = a - 2 * b + c;
    if (curvature < 0 && Math.abs(p[k].z - p[k - 1].z - h) < 1e-9) { const shift = 0.5 * (a - c) / curvature; depth += shift * h; speed = b - 0.25 * (a - c) * shift; }
  }
  return { speed, depth };
}

export const PROFILE_DEPTHS = Array.from({ length: 67 }, (_, j) => 20 + 5 * j);
export const meshSpacing = (cells) => 0.9 * Math.sqrt(4 * Math.PI / cells) * DEG;
export const inMooringWindow = (latDeg, lonDeg, stationLon, spacing) => Math.abs(latDeg) <= Math.max(1, spacing) && Math.abs(wrap(lonDeg - stationLon)) <= Math.max(1.5, spacing);

export function mooringProfile(columns, Z = PROFILE_DEPTHS) {
  const sum = new Float64Array(Z.length), weight = new Float64Array(Z.length);
  for (const { area, z: zc, u: uc } of columns) {
    for (let j = 0; j < Z.length; j++) {
      const z = Z[j];
      if (zc.length < 2 || z < zc[0] || z > zc[zc.length - 1]) continue;
      let n = 0;
      while (zc[n + 1] < z) n++;
      const f = (z - zc[n]) / (zc[n + 1] - zc[n]);
      sum[j] += area * (uc[n] + f * (uc[n + 1] - uc[n])); weight[j] += area;
    }
  }
  return Z.map((z, j) => ({ z, u: weight[j] ? sum[j] / weight[j] : NaN }));
}

const meshes = new Map();
function meshOf(saved) {
  const { topography, geography } = stateGeography(saved), key = `${saved.N} ${JSON.stringify(saved.grid ?? null)} ${JSON.stringify(geography)}`;
  if (!meshes.has(key)) { const { mesh } = stateGrid(saved, topography), made = createGeography(mesh, topography, geography); meshes.set(key, { mesh, land: made.land, geography: made }); }
  return meshes.get(key);
}

let refs = null;
const references = () => (refs ??= { gpcp: R.gpcp().grid, era5: R.era5Fields().fields, woa: R.woa23Sst().sst, ts: R.globalTs(), ice: R.seaIce(), tao: R.tao(), ceres: R.ceres().fields });

export function reduceState(saved, logText, { spacing: mooringSpacing = null } = {}) {
  const { gpcp, era5: { t2m, msl, metss, lsm }, woa, ts, ice: iceRef, tao, ceres } = references();
  const { mesh, land } = meshOf(saved), C = mesh.nCells, E = mesh.nEdges, M = saved.means, area = mesh.areaCell, D = saved.day;
  if (!M?.fields) throw new Error(`the day-${D} state carries no means`);
  const lat = (i) => mesh.latCell[i] * DEG, lon = (i) => wrap(mesh.lonCell[i] * DEG);
  const L = saved.ocean.h.length / C, sea = new Uint8Array(C);
  for (let i = 0; i < C; i++) { let column = 0; for (let k = 0; k < L; k++) column += saved.ocean.h[k * C + i]; sea[i] = column > 1 ? 1 : 0; }
  const cells = (keep) => { const out = []; for (let i = 0; i < C; i++) if (keep(i)) out.push(i); return out; };
  const over = (f, cs) => { let s = 0, a = 0; for (const i of cs) if (Number.isFinite(f[i])) { s += area[i] * f[i]; a += area[i]; } return s / a; };
  const rows = {}, row = (key, model, ref, extra = {}) => { rows[key] = { model, ref, ...extra }; };

  const months = windowMonths(woa, windowDays(D));
  for (const [key, b] of Object.entries({ sst30: { lat: [-30, 30] }, sstCold: { lat: [-2, 2], lon: [-150, -90] }, sstWarm: { lat: [-10, 10], lon: [130, 160] } })) {
    let a = 0, m = 0, w = 0, missing = 0;
    for (let i = 0; i < C; i++) {
      if (!sea[i] || !Number.isFinite(M.sst[i]) || lat(i) < b.lat[0] || lat(i) > b.lat[1] || (b.lon && (lon(i) < b.lon[0] || lon(i) > b.lon[1]))) continue;
      const c = climatologyAt(woa, months, mesh.latCell[i], mesh.lonCell[i]);
      if (!Number.isFinite(c)) { missing++; continue; }
      a += area[i]; m += area[i] * (M.sst[i] > 100 ? M.sst[i] - KELVIN : M.sst[i]); w += area[i] * (c > 100 ? c - KELVIN : c);
    }
    row(key, m / a, w / a, { missing });
  }
  const all = cells(() => true);
  row('ts', M.gridBoxSkin ? over(M.gridBoxSkin, all) - KELVIN : NaN, dayWeighted(D, (m) => ts.monthly[m - 1]), { proxy: over(M.surfaceT, all) - KELVIN });
  const band = cells((i) => land[i] && lat(i) >= 40 && lat(i) <= 60);
  const landRef = dayWeighted(D, (m) => R.boxMean(mesh, t2m, m - 1, band, { lsm, land }));
  row('landT', over(M.temperature, band), landRef);
  row('landT2m', M.t2m ? over(M.t2m, band) : NaN, landRef);
  for (const [box, [kind, s, n, w, e]] of Object.entries(BOXES)) {
    const cs = cells((i) => lat(i) >= s && lat(i) <= n && lon(i) >= w && lon(i) <= e && (kind === 'all' || (kind === 'land') === !!land[i]));
    row(`rain_${box}`, over(M.precipitation, cs), box === 'global' ? dayWeighted(D, (m) => R.gridMean(gpcp, m - 1)) : dayWeighted(D, (m) => R.boxMean(mesh, gpcp, m - 1, cs)), { cells: cs.length });
  }
  const equator = (w, e) => cells((i) => sea[i] && Math.abs(lat(i)) <= 3 && lon(i) >= w && lon(i) <= e);
  const stressCells = equator(-160, -100), west = equator(-180, -160), east = equator(-120, -100);
  row('stress', 1e3 * over(M.zonalStress, stressCells), 1e3 * dayWeighted(D, (m) => R.boxMean(mesh, metss, m - 1, stressCells, { lsm, land })));
  row('slp', (over(M.mslp, west) - over(M.mslp, east)) / 100, dayWeighted(D, (m) => R.boxMean(mesh, msl, m - 1, west, { lsm, land }) - R.boxMean(mesh, msl, m - 1, east, { lsm, land })) / 100);

  const { days, ice } = logDays(logText ?? ''), win = windowDays(D);
  const fromLog = (map, f) => { const xs = win.map((d) => map.get(d)).filter(Boolean); return xs.length === win.length ? avg(xs.map(f)) : NaN; };
  const cer = (f) => dayWeighted(D, (m) => f(m - 1));
  row('asr', fromLog(days, (x) => x.asr), cer((m) => ceres.rsdt[m] - ceres.rsut[m]));
  row('net', fromLog(days, (x) => x.asr - x.olr), cer((m) => ceres.rsdt[m] - ceres.rsut[m] - ceres.rlut[m]));
  row('swcre', fromLog(days, (x) => x.swcre), cer((m) => ceres.rstcre[m]));
  row('lwcre', fromLog(days, (x) => x.lwcre), cer((m) => ceres.rltcre[m]));
  const { N: cutN, S: cutS } = seaIceIndexCut(mesh, indexCells(mesh, meshOf(saved).geography), M.concentration), nsidc = (k, h) => dayWeighted(D, (m) => iceRef[k][h][m - 1]);
  row('iceNExtent', fromLog(ice, (x) => x.nExtent), nsidc('extent', 'N'), { proxy: cutN.extent });
  row('iceNArea', fromLog(ice, (x) => x.nArea), nsidc('area', 'N'), { proxy: cutN.area });
  row('iceSExtent', fromLog(ice, (x) => x.sExtent), nsidc('extent', 'S'), { proxy: cutS.extent });
  row('iceSArea', fromLog(ice, (x) => x.sArea), nsidc('area', 'S'), { proxy: cutS.area });

  const eastward = Array.from({ length: L }, (_, k) => cellVector(mesh, Float64Array.from(saved.ocean.u.subarray(k * E, (k + 1) * E))));
  const spacing = mooringSpacing ?? meshSpacing(C);
  for (const s of R.TAO_STATIONS) {
    const st = tao.stations[s], cs = cells((i) => sea[i] && inMooringWindow(lat(i), lon(i), st.longitude, spacing));
    const columns = cs.map((i) => {
      const z = [], u = [], sinL = Math.sin(mesh.lonCell[i]), cosL = Math.cos(mesh.lonCell[i]);
      let top = 0;
      for (let k = 0; k < L; k++) {
        const h = saved.ocean.h[k * C + i];
        if (h > 1e-3) { z.push(top + h / 2); u.push(-sinL * eastward[k][3 * i] + cosL * eastward[k][3 * i + 1]); }
        top += h;
      }
      return { area: area[i], z, u };
    });
    const core = undercurrentCore(mooringProfile(columns));
    row(`euc_${s}`, core.speed, st.core.speed, { cells: cs.length });
    row(`eucDepth_${s}`, core.depth, st.core.depth, { cells: cs.length });
  }
  return { day: D, N: saved.N, rows };
}

export function windowMeans(reductions, { proxy = false } = {}) {
  const out = {};
  for (const [key] of ROWS) {
    const pick = (r) => { const x = r.rows[key] ?? {}; return Number.isFinite(x.model) ? x.model : proxy && Number.isFinite(x.proxy) ? x.proxy : NaN; };
    out[key] = { model: avg(reductions.map(pick)), ref: avg(reductions.map((r) => r.rows[key]?.ref ?? NaN)) };
  }
  return out;
}

const pickOf = (key, proxy = true) => (r) => { const x = r.rows[key] ?? {}; return Number.isFinite(x.model) ? x.model : proxy && Number.isFinite(x.proxy) ? x.proxy : NaN; };
const rms = (xs) => Math.sqrt(avg(xs.map((x) => x * x)));
export const GATE_LEADS = [30, 60];

export function twoSidedTail(z) {
  const x = Math.abs(z) / Math.SQRT2, t = 1 / (1 + 0.3275911 * x);
  return t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429)))) * Math.exp(-x * x);
}

const lgamma = (z) => {
  const c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
  if (z < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * z)) - lgamma(1 - z);
  let x = c[0];
  for (let i = 1; i < 9; i++) x += c[i] / (z - 1 + i);
  const t = z + 6.5;
  return 0.5 * Math.log(2 * Math.PI) + (z - 0.5) * Math.log(t) - t + Math.log(x);
};

function betaFraction(a, b, x) {
  const tiny = 1e-300, clamp = (v) => (Math.abs(v) < tiny ? tiny : v);
  let c = 1, d = 1 / clamp(1 - (a + b) * x / (a + 1)), h = d;
  for (let m = 1; m <= 500; m++) {
    const even = m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m));
    d = 1 / clamp(1 + even * d); c = clamp(1 + even / c); h *= d * c;
    const odd = -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1));
    d = 1 / clamp(1 + odd * d); c = clamp(1 + odd / c);
    const step = d * c;
    h *= step;
    if (Math.abs(step - 1) < 1e-14) break;
  }
  return h;
}

function incompleteBeta(a, b, x) {
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  const front = Math.exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * Math.log(x) + b * Math.log(1 - x));
  return x < (a + 1) / (a + b + 2) ? front * betaFraction(a, b, x) / a : 1 - front * betaFraction(b, a, 1 - x) / b;
}

export const studentTail = (t, nu) => incompleteBeta(nu / 2, 0.5, nu / (nu + t * t));

export function studentCritical(nu, p = 0.05) {
  if (!(nu > 0)) return 2;
  let lo = 0, hi = 1e4;
  for (let i = 0; i < 200; i++) { const mid = (lo + hi) / 2; if (studentTail(mid, nu) > p) lo = mid; else hi = mid; }
  return (lo + hi) / 2;
}

export function sigmaDof(windowsPerStart, r) {
  const rho = (k) => (k === 0 ? 1 : r ** Math.abs(k));
  let trace = 0, square = 0;
  for (const n of windowsPerStart) for (let i = 0; i < n - 1; i++) for (let j = 0; j < n - 1; j++) {
    const c = (2 * rho(j - i) + rho(j - i - 1) + rho(j - i + 1)) / 4;
    if (i === j) trace += c;
    square += c * c;
  }
  return square > 0 ? trace * trace / square : NaN;
}

export function calibrate(twins, { from = 180, leads = GATE_LEADS } = {}) {
  const rows = {};
  for (const [key] of ROWS) {
    const pick = pickOf(key), series = twins.map((t) => t.a.map((r, k) => ({ lead: r.day - t.start, d: pick(r) - pick(t.b[k]) })).filter((x) => Number.isFinite(x.d)));
    const late = series.map((sr) => sr.filter((x) => x.lead >= from).map((x) => x.d)), all = late.flat();
    if (!all.length) { rows[key] = { sigma: null, windows: 0, ungated: 'no twin pair carries this row' }; continue; }
    const pairs = late.flatMap((ds) => ds.slice(1).map((d, k) => (ds[k] + d) / 2));
    const lag = late.flatMap((ds) => ds.slice(1).map((d, k) => d * ds[k]));
    const perWindow = rms(all), lag1 = Math.max(-0.9, Math.min(0.9, lag.length ? avg(lag) / (perWindow * perWindow) : 0));
    const sigma = rms(pairs), dof = sigmaDof(late.map((ds) => ds.length), lag1), critical = studentCritical(dof), threshold = critical * sigma;
    const at = (sr, lead) => sr.find((x) => x.lead === lead)?.d;
    const early = leads.map((lead) => rms(series.map((sr) => at(sr, lead)).filter(Number.isFinite)));
    const gate = series.map((sr) => avg(leads.map((lead) => at(sr, lead))));
    const spread = Math.sqrt(early[0] ** 2 + early[1] ** 2 + 2 * lag1 * early[0] * early[1]) / 2;
    rows[key] = {
      sigma, perWindow, lag1, windows: all.length, pairs: pairs.length, dof, critical, threshold, early, earlyGate: rms(gate.filter(Number.isFinite)), saturation: rms(gate.filter(Number.isFinite)) / sigma,
      nullAtLeads: gate.map((g) => (Number.isFinite(g) ? Math.abs(g) > threshold : null)), normalTailRate: spread > 0 ? twoSidedTail(threshold / spread) : 0, inSampleRate: pairs.filter((m) => Math.abs(m) > threshold).length / pairs.length,
      proxy: twins.some((t) => t.a.some((r) => !Number.isFinite(r.rows[key]?.model))),
    };
  }
  const values = Object.values(rows).filter((r) => !r.ungated), mean = (f) => avg(values.map(f));
  const nulls = {
    atLeads: twins.map((t, k) => ({ start: t.label, moved: values.filter((r) => r.nullAtLeads[k] === true).length, rows: values.filter((r) => r.nullAtLeads[k] !== null).length })),
    normalTailRate: mean((r) => r.normalTailRate), inSampleRate: mean((r) => r.inSampleRate), saturatedRateAtTwoSigma: mean((r) => studentTail(2, r.dof)),
  };
  return { rows, nulls };
}

export function designProblems(base, now) {
  const problems = [], leads = (x) => JSON.stringify((x.days ?? []).map((d) => d - (x.start?.day ?? NaN)));
  if (!base.start?.sha256 || !now.start?.sha256) problems.push(`the ${base.start?.sha256 ? 'candidate' : 'base'} records no start`);
  else if (base.start.sha256 !== now.start.sha256) problems.push(`the base started from ${base.start.sha256.slice(0, 12)}…, the candidate from ${now.start.sha256.slice(0, 12)}…`);
  if (leads(base) !== leads(now)) problems.push(`the base's leads are ${leads(base)}, the candidate's ${leads(now)}`);
  return problems;
}

const valueOf = (x) => (typeof x === 'number' && Number.isFinite(x) ? x : null);

export function compareScores(base, now, sigma, { anyStart = false } = {}) {
  const problems = designProblems(base, now);
  if (problems.length && !anyStart) throw new Error(`not the same gate design: ${problems.join('; ')} (--any-start compares them anyway)`);
  if (problems.length) console.warn(`--any-start: comparing anyway, though ${problems.join('; ')}`);
  return ROWS.map(([key]) => {
    const a = valueOf(base.means?.[key]?.model), b = valueOf(now.means?.[key]?.model), s = valueOf(sigma[key]?.sigma), threshold = s === null ? NaN : s * studentCritical(sigma[key]?.dof);
    const gated = Number.isFinite(threshold);
    if (a === null || b === null) return { key, delta: NaN, sigma: s, threshold, gated, identical: false, missing: [...(a === null ? ['base'] : []), ...(b === null ? ['candidate'] : [])], verdict: 'no value' };
    const delta = b - a;
    return { key, delta, sigma: s, threshold, gated, identical: a === b, verdict: !gated ? 'ungated' : Math.abs(delta) > threshold ? 'moved' : 'unresolved' };
  });
}

export function identityCheck(base, now, comparison) {
  const gated = comparison.filter((c) => c.gated), describe = (c) => (c.verdict === 'no value' ? `${c.key} (no value in the ${c.missing.join(' and ')})` : `${c.key} (${c.identical ? 'identical' : 'changed'})`);
  const empty = gated.filter((c) => c.verdict === 'no value').map(describe);
  const changed = gated.filter((c) => c.verdict !== 'no value' && !c.identical).map((c) => c.key);
  const ungated = comparison.filter((c) => !c.gated).map(describe);
  const sameStart = !!base.start?.sha256 && base.start.sha256 === now.start?.sha256;
  const sameSaves = Array.isArray(base.saves) && Array.isArray(now.saves) && base.saves.length === now.saves.length && base.saves.every((h, k) => h === now.saves[k]);
  return { ok: !empty.length && !changed.length && sameStart, empty, changed, ungated, sameStart, sameSaves, baseHasSaves: Array.isArray(base.saves) };
}

export function formatTable(reductions, sigma, comparison = null) {
  const f = (x, n = 2) => (Number.isFinite(x) ? x.toFixed(n) : '—');
  const head = `| row | ${reductions.map((r) => `day ${r.day} model / ref / bias`).join(' | ')} | pair σ / moves past |${comparison ? ' Δ vs base | verdict |' : ''}`;
  const lines = [head, `|---|${reductions.map(() => '---|').join('')}---|${comparison ? '---|---|' : ''}`];
  for (const [key, name, units] of ROWS) {
    const cells = reductions.map((r) => { const x = r.rows[key]; return `${f(x.model)} / ${f(x.ref)} / **${f(x.model - x.ref)}**`; });
    const c = comparison?.find((x) => x.key === key);
    lines.push(`| ${name} (${units}) | ${cells.join(' | ')} | ${f(sigma[key]?.sigma, 3)}${sigma[key]?.proxy ? '*' : ''} / ${f(valueOf(sigma[key]?.sigma) === null ? NaN : sigma[key].sigma * studentCritical(sigma[key].dof), 3)} |${c ? ` ${f(c.delta, 3)} | ${c.verdict} |` : ''}`);
  }
  return lines.join('\n');
}

const sha256 = (bytes) => createHash('sha256').update(bytes).digest('hex');
const readLog = (file) => { if (!existsSync(file)) throw new Error(`no log at ${file}; give the run's log with --log FILE`); return readFileSync(file, 'utf8'); };
const logBeside = (path) => join(dirname(path), `${basename(path).replace(/_day\d+\.bin$/, '')}.log`);
async function reduceFile(path, logFile) {
  const bytes = readFileSync(path), saved = await decodeState(new Uint8Array(bytes));
  return { ...reduceState(saved, readLog(logFile ?? logBeside(path))), sha256: sha256(bytes), options: saved.options ?? null, given: saved.given ?? null };
}

async function startOf(path) {
  if (!path || !existsSync(path)) return null;
  const bytes = readFileSync(path);
  return { path, sha256: sha256(bytes), day: (await decodeState(new Uint8Array(bytes))).day };
}

const ROOT = fileURLToPath(new URL('..', import.meta.url));
const GIT = { cwd: ROOT, maxBuffer: 1 << 30, stdio: ['ignore', 'pipe', 'ignore'] };
const git = (cmd, cwd = ROOT) => execSync(cmd, { ...GIT, cwd, encoding: 'utf8' }).trim();

function untrackedEntry(root, path) {
  const full = join(root, path);
  try {
    const st = lstatSync(full);
    if (st.isSymbolicLink()) return `link ${readlinkSync(full)}`;
    if (st.isDirectory()) return 'directory';
    return readFileSync(full);
  } catch (error) { return `unreadable ${error.code ?? error.message}`; }
}

export function fingerprint(root = ROOT) {
  try {
    const head = git('git rev-parse --short=12 HEAD', root), diff = execSync('git diff HEAD -- js scripts data', { ...GIT, cwd: root });
    const untracked = execSync('git ls-files -z --others --exclude-standard -- js scripts data', { ...GIT, cwd: root }).toString().split('\0').filter(Boolean);
    const hash = createHash('sha256').update(head).update(diff);
    for (const path of untracked) hash.update(`\0${path}\0`).update(untrackedEntry(root, path));
    return { scoredOn: `${head}${diff.length || untracked.length ? ' with changes' : ''}`, tree: hash.digest('hex').slice(0, 16) };
  } catch (error) {
    console.warn(`warning: no fingerprint of the tree (${String(error.message).split('\n')[0]}); scoredOn and tree are null`);
    return { scoredOn: null, tree: null };
  }
}

function spinupEnv(env) {
  const source = readFileSync(join(ROOT, 'scripts/spinup.mjs'), 'utf8');
  return Object.fromEntries(Object.keys(env).filter((k) => /^[A-Z][A-Z0-9_]+$/.test(k) && new RegExp(`\\b${k}\\b`).test(source)).sort().map((k) => [k, env[k]]));
}

function commitArg(value) {
  if (value === null) return null;
  if (!/^[0-9a-f]{7,40}$/.test(value)) throw new Error(`--commit ${value} is not a commit hash`);
  try { return git(`git rev-parse --verify --short=12 ${value}^{commit}`); } catch {}
  if (value.length < 12) throw new Error(`--commit ${value} is not a commit in this repository and too short to record as 12 characters`);
  return value.slice(0, 12);
}

export function sameFile(x, y) {
  if (existsSync(x) && existsSync(y)) { const [a, b] = [statSync(x), statSync(y)]; return a.dev === b.dev && a.ino === b.ino; }
  return resolve(x) === resolve(y);
}

export function runSegments(days, run, { first, take = () => fingerprint().tree } = {}) {
  const segments = [];
  days.forEach((day, k) => { segments.push({ day, tree: k === 0 ? first : take() }); run(day); });
  return segments;
}

function loadBase(baseFile, jsonOut) {
  if (!baseFile) return null;
  if (jsonOut && sameFile(baseFile, jsonOut)) throw new Error(`--json ${jsonOut} is the --base file; it would overwrite the base`);
  return JSON.parse(readFileSync(baseFile, 'utf8'));
}

function score({ reductions, start, made = null, commit = null, base = null, jsonOut, anyStart, identical = false, extra = {} }) {
  const sigma = JSON.parse(readFileSync(SIGMA_FILE, 'utf8')).rows, means = windowMeans(reductions), scored = fingerprint();
  const trees = made ? [made.tree, ...made.segments.map((x) => x.tree), scored.tree] : [];
  const now = {
    start: start && { sha256: start.sha256, day: start.day, path: start.path }, days: reductions.map((r) => r.day), saves: reductions.map((r) => r.sha256),
    commit: made ? made.scoredOn : commit, scoredOn: scored.scoredOn, tree: scored.tree, env: made?.env ?? null, ...(made ? { segments: made.segments } : {}),
    ...(made && new Set(trees).size > 1 ? { treeChangedDuringRun: { before: made.scoredOn, trees, after: scored.scoredOn } } : {}),
    options: reductions[0]?.options ?? null, given: reductions[0]?.given ?? null, ...extra, means,
  };
  if (now.treeChangedDuringRun) console.warn(`warning: the tree changed during the run (fingerprints ${trees.join(' → ')}); the commit recorded is the one at the start`);
  const write = (comparison) => jsonOut && writeFileSync(jsonOut, JSON.stringify({ ...now, reductions, comparison }, null, 1) + '\n');
  write(null);
  const comparison = base ? compareScores(base, now, sigma, { anyStart }) : null;
  if (comparison) write(comparison);
  console.log(formatTable(reductions, sigma, comparison));
  if (!comparison) return;
  const check = identityCheck(base, now, comparison);
  console.log(`moved: ${comparison.filter((c) => c.verdict === 'moved').map((c) => c.key).join(', ') || 'none'}; identical (Δ exactly 0): ${comparison.filter((c) => c.gated && c.identical).length} of ${comparison.filter((c) => c.gated).length} gated rows${check.empty.length ? `; no value: ${check.empty.join(', ')}` : ''}${check.ungated.length ? `; ungated, not checked: ${check.ungated.join(', ')}` : ''}`);
  if (!identical) return;
  const saves = !check.baseHasSaves ? 'the base records no save hashes' : check.sameSaves ? 'byte-identical saves (the base\'s own, or a run that reproduced them)' : 'different saves';
  if (check.ok) console.log(`--identical: all ${comparison.filter((c) => c.gated).length} gated rows bit-identical from the same start; ${saves}`);
  else {
    console.error(`--identical: ${[!check.sameStart && 'the recorded starts differ', check.changed.length && `${check.changed.join(', ')} changed`, check.empty.length && `no value: ${check.empty.join(', ')}`].filter(Boolean).join('; ')}`);
    process.exitCode = 1;
  }
}

if (process.argv[1] && fileURLToPath(import.meta.url) === (await import('node:fs')).realpathSync(process.argv[1])) {
  const args = process.argv.slice(2), opt = (name) => { const k = args.indexOf(name); return k < 0 ? null : args.splice(k, 2)[1]; };
  const flag = (name) => { const k = args.indexOf(name); if (k >= 0) args.splice(k, 1); return k >= 0; };
  const anyStart = flag('--any-start'), identical = flag('--identical');
  if (args[0] === '--run') {
    const jsonOut = opt('--json'), baseFile = opt('--base'), [, start, outDir] = args, usage = 'usage: --run START OUT [--base BASE.json] [--json OUT.json] [--any-start] [--identical]';
    if (!outDir || args.length !== 3 || [start, outDir].some((a) => a.startsWith('--'))) throw new Error(`${usage}; not ${JSON.stringify(args.slice(1))}`);
    if (identical && !baseFile) throw new Error('--identical needs --base');
    const base = loadBase(baseFile, jsonOut);
    const bytes = readFileSync(start), hash = sha256(bytes), first = await decodeState(new Uint8Array(bytes));
    if (first.N !== 64) throw new Error(`${start} is N=${first.N}, not 64`);
    if (hash !== STANDARD_START && !anyStart) throw new Error(`${start} is not the standard start (sha256 ${hash.slice(0, 12)}…, the standard ${STANDARD_START.slice(0, 12)}…); --any-start runs it anyway`);
    mkdirSync(outDir, { recursive: true });
    if (readdirSync(outDir).some((name) => /_day\d+\.bin$/.test(name))) throw new Error(`${outDir} already holds saves`);
    const D0 = first.day, days = [D0 + 30, D0 + 60], spinup = fileURLToPath(new URL('./spinup.mjs', import.meta.url));
    copyFileSync(start, join(outDir, `p64_day${String(D0).padStart(4, '0')}.bin`));
    const env = { ...process.env, LEVELS: 'bl48', N: '64', TAG: 'p64', OUT: outDir, MINUTES: '100000', KEEP: '20' }, made = { ...fingerprint(), env: spinupEnv(env) };
    made.segments = runSegments(days, (d) => {
      const run = spawnSync(process.execPath, ['--max-old-space-size=6000', spinup], { stdio: 'inherit', env: { ...env, DAYS: String(d) } });
      if (run.status !== 0) throw new Error(`spinup.mjs to day ${d} exited ${run.status}`);
    }, { first: made.tree });
    const reductions = [];
    for (const d of days) reductions.push(await reduceFile(join(outDir, `p64_day${String(d).padStart(4, '0')}.bin`)));
    score({ reductions, start: { path: start, sha256: hash, day: D0 }, made, base, jsonOut, anyStart, identical });
  } else if (args[0] === '--sigma') {
    const out = opt('--out') ?? SIGMA_FILE, from = Number(opt('--from') ?? 180), cacheFile = opt('--cache');
    const groups = args.slice(1), twins = [], cache = cacheFile && existsSync(cacheFile) ? JSON.parse(readFileSync(cacheFile, 'utf8')) : {};
    if (!groups.length || groups.length % 3) throw new Error('usage: --sigma [--from LEAD] [--out PATH] [--cache FILE] A_DIR B_DIR START_DAY ...');
    const tagOf = (dir) => (existsSync(join(dir, 'p64.log')) ? 'p64' : basename(dir));
    const savesAfter = (dir, start) => readdirSync(dir).map((f) => f.match(/_day(\d+)\.bin$/)).filter(Boolean).map((m) => Number(m[1])).filter((d) => d > start).sort((x, y) => x - y);
    const reduceDir = async (dir, days) => {
      const out = [];
      for (const d of days) {
        const key = `${dir}:${d}`;
        cache[key] ??= await reduceFile(join(dir, `${tagOf(dir)}_day${String(d).padStart(4, '0')}.bin`), join(dir, `${tagOf(dir)}.log`));
        out.push(cache[key]);
      }
      return out;
    };
    for (let k = 0; k < groups.length; k += 3) {
      const [a, b, start] = [groups[k], groups[k + 1], Number(groups[k + 2])], days = savesAfter(a, start).filter((d) => savesAfter(b, start).includes(d));
      twins.push({ label: `${basename(dirname(a))}/${basename(a)} − ${basename(dirname(b))}/${basename(b)} from day ${start}`, start, a: await reduceDir(a, days), b: await reduceDir(b, days) });
      console.log(`${twins.at(-1).label}: ${days.length} windows, days ${days[0]}–${days.at(-1)}`);
    }
    if (cacheFile) writeFileSync(cacheFile, JSON.stringify(cache));
    const { rows, nulls } = calibrate(twins, { from });
    const record = {
      name: 'Pair sigma of the N=64 scorecard\'s gate statistic',
      design: `twin controls (the second perturbed by 1e-6 at its start) from ${twins.length} starts, each a year of 30-day windows; sigma is the RMS of the two-window mean difference (the gate's statistic) over the windows at leads of ${from} days or more, where the twins have decorrelated, pooled over the starts`,
      definition: 'per row: sigma, the pair sigma of the mean of two adjacent windows; perWindow, the RMS of one window\'s difference; lag1, the correlation of adjacent windows\' differences; dof, sigma\'s degrees of freedom, (tr C)² / tr(C²) for the pair means\' covariance C under a Gaussian AR(1) at lag1 with the starts independent; critical, the two-sided 5 % quantile of Student\'s t at dof; threshold, critical × sigma, past which a row moves; early, the RMS at leads 30 and 60 days (the gate\'s), earlyGate their two-window mean, saturation earlyGate / sigma; nullAtLeads, whether each start\'s twins at the gate\'s leads would read moved; normalTailRate, P(|mean| > threshold) for normal differences at the early RMS with lag1, a closed form; inSampleRate, the share of the late pairs past threshold. nulls: atLeads, the rows each start\'s twins move at the gate\'s leads; normalTailRate and inSampleRate, their means over rows; saturatedRateAtTwoSigma, the mean over rows of P(|t| > 2) at dof, what a fixed 2 sigma would move under a decorrelated null. ungated: a row no twin pair carries, with no sigma, reported but never moved or checked by --identical. proxy: the twins carry no gridBoxSkin or daily sea-ice line, so Ts takes surfaceT and sea ice the 30-day-mean cut on the index\'s grids',
      twins: twins.map((t) => ({ label: t.label, start: t.start, days: t.a.map((r) => r.day) })), from, measured: new Date().toISOString().slice(0, 10),
      expires: 'remeasure whenever a row changes its instrument or the control physics changes, and with more starts when they exist',
      script: 'scripts/scorecard64.mjs --sigma', nulls, rows,
    };
    writeFileSync(out, JSON.stringify(record, null, 1) + '\n');
    const f = (x, n = 3) => (Number.isFinite(x) ? x.toFixed(n) : '—');
    console.log('| row | σ (gate) | per window | lag1 | dof | early 30 / 60 | saturation | null at leads | normal tail |');
    console.log('|---|---|---|---|---|---|---|---|---|');
    for (const [key, name, units] of ROWS) { const r = rows[key]; if (r.ungated) { console.log(`| ${name} (${units}) | ungated: ${r.ungated} |`); continue; } console.log(`| ${name} (${units}) | ${f(r.sigma)}${r.proxy ? '*' : ''} | ${f(r.perWindow)} | ${f(r.lag1, 2)} | ${f(r.dof, 1)} | ${f(r.early[0])} / ${f(r.early[1])} | ${f(r.saturation, 2)} | ${r.nullAtLeads.map((x) => (x ? 'moved' : x === false ? '·' : '—')).join(' ')} | ${f(100 * r.normalTailRate, 1)} % |`); }
    console.log(`nulls: ${nulls.atLeads.map((n) => `${n.start}: ${n.moved} of ${n.rows} moved`).join('; ')}; normal tail at the twins' spread ${(100 * nulls.normalTailRate).toFixed(1)} %, in-sample ${(100 * nulls.inSampleRate).toFixed(1)} %; a saturated null at 2σ would move ${(100 * nulls.saturatedRateAtTwoSigma).toFixed(1)} %`);
  } else {
    const jsonOut = opt('--json'), baseFile = opt('--base'), logFile = opt('--log'), startFile = opt('--start'), commit = commitArg(opt('--commit'));
    if (identical && !baseFile) throw new Error('--identical needs --base');
    const base = loadBase(baseFile, jsonOut);
    if (args.length !== 2 || !args.every((a) => /_day\d+\.bin$/.test(a))) throw new Error(`give the gate's two saves, not ${JSON.stringify(args)}: [--json OUT.json] [--base BASE.json] [--log FILE] [--start START] [--commit HASH] [--any-start] [--identical] DIR/p64_day0395.bin DIR/p64_day0425.bin`);
    const reductions = [];
    for (const path of args) reductions.push(await reduceFile(path, logFile));
    const first = args[0] && reductions[0], beside = first && join(dirname(args[0]), basename(args[0]).replace(/_day\d+\.bin$/, `_day${String(first.day - GATE_LEADS[0]).padStart(4, '0')}.bin`));
    score({ reductions, start: await startOf(startFile ?? beside), commit, base, jsonOut, anyStart, identical, extra: { states: args } });
  }
}
