#!/usr/bin/env node
/*
 * The N=64 observational scorecard, the per-merge gate's second half (the first is scripts/physicsDigest.mjs). Each
 * row is a 30-day day-mean window ending at a save, scored against its definition-matched reference in
 * data/references/ (README.md there gives each one's definition, period and box; scripts/references.mjs reads them).
 *
 * The gate run starts from the standard start, seventeen64_day0365.bin (the March equinox; sha256 STANDARD_START, kept at
 * ~/.cache/webgcm-scratch/seventeen/out64/ on this Mac: the page's runs/seventeen64_day0365 parts drop its means and SST
 * ring). Under gpu_run (about 6 GPU minutes, level 2):
 *   node scripts/scorecard64.mjs --run START OUT [--base BASE.json] [--json OUT.json]
 * copies START into the empty directory OUT as p64_day0365.bin, runs LEVELS=bl48 N=64 TAG=p64 spinup.mjs on this tree
 * to days 395 and 425, and scores the two saves; without the copy spinup.mjs would start fresh and run 395 days. Saves
 * already made are scored with
 *   node scripts/scorecard64.mjs [--json OUT.json] [--base BASE.json] DIR/p64_day0395.bin DIR/p64_day0425.bin
 * Each state's daily lines (TOA fluxes and the sea-ice cuts) are read from TAG.log beside it, or from --log FILE.
 * With --base, the rows whose two-window mean moved by more than twice the pair sigma are listed as moved and the rest
 * as unresolved; a re-pin cites the moved rows.
 *
 * Sigma (data/scorecard64-sigma.json) is the pair sigma of that two-window mean: the RMS, over adjacent windows at
 * leads of 180 days or more, of the difference between twin controls (the second perturbed by 1e-6 at its start), from
 * the two starts of the pairs (day 365 in March and day 548 in September), a year each. By 180 days the twins have
 * decorrelated (the RMS at shorter leads is a tenth of it at 30 days and half at 60), so sigma is the saturated
 * spread, the conservative choice for a merge whose changes perturb more than 1e-6: twins at the gate's own leads move
 * 0 of 30 rows at 2 sigma, and normal differences at their measured spread 0.2 %. Each row records its degrees of
 * freedom, its 14 windows less their serial correlation: under 5 for the SSTs, the warm-pool and tropical-sea rain,
 * net TOA, SWCRE, the Arctic area, the Antarctic extent and four undercurrent rows, so their sigma is itself uncertain
 * by tens of percent. --sigma [--from LEAD] [--cache FILE] A_DIR B_DIR START ... rewrites it from twin pairs; rows the
 * twins cannot give (no gridBoxSkin or daily ice line) take the named proxy.
 *
 * Not definition-matched: the 40–60N land row compares the lowest layer's air with ERA5's 2 m air until the model has a
 * 2 m diagnostic, and the undercurrent rows compare the window-end ocean's velocity, a snapshot, with TAO's annual-mean
 * climatology, not day means of the window's calendar months.
 */
import { readFileSync, writeFileSync, existsSync, readdirSync, mkdirSync, copyFileSync } from 'node:fs';
import { spawnSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { basename, dirname, join } from 'node:path';
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
  ['ts', 'global Ts, grid-box skin', '°C'], ['landT', '40–60N land air: lowest layer against 2 m (not matched)', 'K'],
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
  const days = new Map(), ice = new Map();
  for (const line of text.split('\n')) {
    let m = line.match(/^day (\d+) \([^)]*\): Ts [-\d.]+ °C, ASR ([-\d.]+) \(atmosphere [-\d.]+\) OLR ([-\d.]+) W\/m²/);
    if (m) {
      const g = (re) => { const x = line.match(re); return x ? Number(x[1]) : NaN; };
      days.set(+m[1], { asr: +m[2], olr: +m[3], swcre: g(/SWCRE ([-\d.]+)/), lwcre: g(/LWCRE ([-\d.]+)/) });
    }
    m = line.match(/sea ice on day (\d+)'s mean concentration, sea cells(?: on the index's grids)?, 15 %: N extent ([\d.]+) area ([\d.]+), S extent ([\d.]+) area ([\d.]+)/);
    if (m) ice.set(+m[1], { nExtent: +m[2], nArea: +m[3], sExtent: +m[4], sArea: +m[5] });
  }
  return { days, ice };
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

const meshes = new Map();
function meshOf(saved) {
  const { topography, geography } = stateGeography(saved), key = `${saved.N} ${JSON.stringify(saved.grid ?? null)} ${JSON.stringify(geography)}`;
  if (!meshes.has(key)) { const { mesh } = stateGrid(saved, topography), made = createGeography(mesh, topography, geography); meshes.set(key, { mesh, land: made.land, geography: made }); }
  return meshes.get(key);
}

let refs = null;
const references = () => (refs ??= { gpcp: R.gpcp().grid, era5: R.era5Fields().fields, woa: R.woa23Sst().sst, ts: R.globalTs(), ice: R.seaIce(), tao: R.tao(), ceres: R.ceres().fields });

export function reduceState(saved, logText) {
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
  row('landT', over(M.temperature, band), dayWeighted(D, (m) => R.boxMean(mesh, t2m, m - 1, band, { lsm, land })));
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
  const Z = Array.from({ length: 67 }, (_, j) => 20 + 5 * j), spacing = 0.9 * Math.sqrt(4 * Math.PI / C) * DEG;
  for (const s of R.TAO_STATIONS) {
    const st = tao.stations[s], cs = cells((i) => sea[i] && Math.abs(lat(i)) <= Math.max(1, spacing) && Math.abs(wrap(lon(i) - st.longitude)) <= Math.max(1.5, spacing));
    const sum = new Float64Array(Z.length), weight = new Float64Array(Z.length);
    for (const i of cs) {
      const zc = [], uc = [], sinL = Math.sin(mesh.lonCell[i]), cosL = Math.cos(mesh.lonCell[i]);
      let top = 0;
      for (let k = 0; k < L; k++) {
        const h = saved.ocean.h[k * C + i];
        if (h > 1e-3) { zc.push(top + h / 2); uc.push(-sinL * eastward[k][3 * i] + cosL * eastward[k][3 * i + 1]); }
        top += h;
      }
      for (let j = 0; j < Z.length; j++) {
        const z = Z[j];
        if (zc.length < 2 || z < zc[0] || z > zc[zc.length - 1]) continue;
        let n = 0;
        while (zc[n + 1] < z) n++;
        const f = (z - zc[n]) / (zc[n + 1] - zc[n]);
        sum[j] += area[i] * (uc[n] + f * (uc[n + 1] - uc[n])); weight[j] += area[i];
      }
    }
    const core = undercurrentCore(Z.map((z, j) => ({ z, u: weight[j] ? sum[j] / weight[j] : NaN })));
    row(`euc_${s}`, core.speed, st.core.speed, { cells: cs.length });
    row(`eucDepth_${s}`, core.depth, st.core.depth, { cells: cs.length });
  }
  return { day: D, N: saved.N, rows };
}

export function windowMeans(reductions, { proxy = false } = {}) {
  const out = {};
  for (const [key] of ROWS) {
    const pick = (r) => { const x = r.rows[key]; return Number.isFinite(x.model) ? x.model : proxy && Number.isFinite(x.proxy) ? x.proxy : NaN; };
    out[key] = { model: avg(reductions.map(pick)), ref: avg(reductions.map((r) => r.rows[key].ref)) };
  }
  return out;
}

const pickOf = (key, proxy = true) => (r) => { const x = r.rows[key]; return Number.isFinite(x.model) ? x.model : proxy && Number.isFinite(x.proxy) ? x.proxy : NaN; };
const rms = (xs) => Math.sqrt(avg(xs.map((x) => x * x)));
export const GATE_LEADS = [30, 60];

export function twoSidedTail(z) {
  const x = Math.abs(z) / Math.SQRT2, t = 1 / (1 + 0.3275911 * x);
  return t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429)))) * Math.exp(-x * x);
}

export function calibrate(twins, { from = 180, leads = GATE_LEADS } = {}) {
  const rows = {};
  for (const [key] of ROWS) {
    const pick = pickOf(key), series = twins.map((t) => t.a.map((r, k) => ({ lead: r.day - t.start, d: pick(r) - pick(t.b[k]) })).filter((x) => Number.isFinite(x.d)));
    const late = series.map((sr) => sr.filter((x) => x.lead >= from).map((x) => x.d)), all = late.flat();
    const pairs = late.flatMap((ds) => ds.slice(1).map((d, k) => (ds[k] + d) / 2));
    const lag = late.flatMap((ds) => ds.slice(1).map((d, k) => d * ds[k]));
    const perWindow = rms(all), lag1 = Math.max(-0.9, Math.min(0.9, lag.length ? avg(lag) / (perWindow * perWindow) : 0));
    const sigma = rms(pairs), dof = Math.max(1, Math.min(all.length, all.length * (1 - lag1) / (1 + lag1)));
    const at = (sr, lead) => sr.find((x) => x.lead === lead)?.d;
    const early = leads.map((lead) => rms(series.map((sr) => at(sr, lead)).filter(Number.isFinite)));
    const gate = series.map((sr) => avg(leads.map((lead) => at(sr, lead))));
    const spread = Math.sqrt(early[0] ** 2 + early[1] ** 2 + 2 * lag1 * early[0] * early[1]) / 2;
    rows[key] = {
      sigma, perWindow, lag1, windows: all.length, pairs: pairs.length, dof, early, earlyGate: rms(gate.filter(Number.isFinite)), saturation: rms(gate.filter(Number.isFinite)) / sigma,
      nullAtLeads: gate.map((g) => (Number.isFinite(g) ? Math.abs(g) > 2 * sigma : null)), monteCarloRate: spread > 0 ? twoSidedTail(2 * sigma / spread) : 0, inSampleRate: pairs.filter((m) => Math.abs(m) > 2 * sigma).length / pairs.length,
      proxy: twins.some((t) => t.a.some((r) => !Number.isFinite(r.rows[key].model))),
    };
  }
  const values = Object.values(rows), mean = (f) => avg(values.map(f));
  const nulls = {
    atLeads: twins.map((t, k) => ({ start: t.label, moved: values.filter((r) => r.nullAtLeads[k] === true).length, rows: values.filter((r) => r.nullAtLeads[k] !== null).length })),
    monteCarloRate: mean((r) => r.monteCarloRate), inSampleRate: mean((r) => r.inSampleRate),
  };
  return { rows, nulls };
}

export function compareScores(base, now, sigma) {
  return ROWS.map(([key]) => {
    const delta = now[key].model - base[key].model, s = sigma[key]?.sigma;
    return { key, delta, sigma: s, verdict: !Number.isFinite(delta) || !Number.isFinite(s) ? 'no value' : Math.abs(delta) > 2 * s ? 'moved' : 'unresolved' };
  });
}

export function formatTable(reductions, sigma, comparison = null) {
  const f = (x, n = 2) => (Number.isFinite(x) ? x.toFixed(n) : '—');
  const head = `| row | ${reductions.map((r) => `day ${r.day} model / ref / bias`).join(' | ')} | pair σ |${comparison ? ' Δ vs base | verdict |' : ''}`;
  const lines = [head, `|---|${reductions.map(() => '---|').join('')}---|${comparison ? '---|---|' : ''}`];
  for (const [key, name, units] of ROWS) {
    const cells = reductions.map((r) => { const x = r.rows[key]; return `${f(x.model)} / ${f(x.ref)} / **${f(x.model - x.ref)}**`; });
    const c = comparison?.find((x) => x.key === key);
    lines.push(`| ${name} (${units}) | ${cells.join(' | ')} | ${f(sigma[key]?.sigma, 3)}${sigma[key]?.proxy ? '*' : ''} |${c ? ` ${f(c.delta, 3)} | ${c.verdict} |` : ''}`);
  }
  return lines.join('\n');
}

const readLog = (file) => (existsSync(file) ? readFileSync(file, 'utf8') : '');
const logBeside = (path) => join(dirname(path), `${basename(path).replace(/_day\d+\.bin$/, '')}.log`);
async function reduceFile(path, logFile) {
  const saved = await decodeState(new Uint8Array(readFileSync(path)));
  return reduceState(saved, readLog(logFile ?? logBeside(path)));
}

if (process.argv[1] && fileURLToPath(import.meta.url) === (await import('node:fs')).realpathSync(process.argv[1])) {
  const args = process.argv.slice(2), opt = (name) => { const k = args.indexOf(name); return k < 0 ? null : args.splice(k, 2)[1]; };
  if (args[0] === '--run') {
    const [, start, outDir] = args, jsonOut = opt('--json'), baseFile = opt('--base');
    if (!outDir) throw new Error('usage: --run START OUT [--base BASE.json] [--json OUT.json]');
    const bytes = readFileSync(start), hash = createHash('sha256').update(bytes).digest('hex'), first = await decodeState(new Uint8Array(bytes));
    if (first.N !== 64) throw new Error(`${start} is N=${first.N}, not 64`);
    if (hash !== STANDARD_START) console.warn(`warning: ${start} is not the standard start (sha256 ${hash.slice(0, 12)}…)`);
    mkdirSync(outDir, { recursive: true });
    if (readdirSync(outDir).some((name) => /_day\d+\.bin$/.test(name))) throw new Error(`${outDir} already holds saves`);
    const D0 = first.day, days = [D0 + 30, D0 + 60], spinup = fileURLToPath(new URL('./spinup.mjs', import.meta.url));
    copyFileSync(start, join(outDir, `p64_day${String(D0).padStart(4, '0')}.bin`));
    for (const d of days) {
      const run = spawnSync(process.execPath, ['--max-old-space-size=6000', spinup], { stdio: 'inherit', env: { ...process.env, LEVELS: 'bl48', N: '64', TAG: 'p64', OUT: outDir, DAYS: String(d), MINUTES: '100000', KEEP: '20' } });
      if (run.status !== 0) throw new Error(`spinup.mjs to day ${d} exited ${run.status}`);
    }
    const reductions = [];
    for (const d of days) reductions.push(await reduceFile(join(outDir, `p64_day${String(d).padStart(4, '0')}.bin`)));
    const sigma = JSON.parse(readFileSync(SIGMA_FILE, 'utf8')).rows, means = windowMeans(reductions);
    const comparison = baseFile ? compareScores(JSON.parse(readFileSync(baseFile, 'utf8')).means, means, sigma) : null;
    if (jsonOut) writeFileSync(jsonOut, JSON.stringify({ start: { path: start, sha256: hash }, states: days, reductions, means, comparison }, null, 1) + '\n');
    console.log(formatTable(reductions, sigma, comparison));
    if (comparison) console.log(`moved: ${comparison.filter((c) => c.verdict === 'moved').map((c) => c.key).join(', ') || 'none'}`);
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
      definition: 'per row: sigma, the pair sigma of the mean of two adjacent windows (a merge moves a row past 2 sigma); perWindow, the RMS of one window\'s difference; lag1, the correlation of adjacent windows\' differences; dof, the windows less their serial correlation, windows (1 − lag1) / (1 + lag1); early, the RMS at leads 30 and 60 days (the gate\'s), earlyGate their two-window mean, saturation earlyGate / sigma; nullAtLeads, whether each start\'s twins at the gate\'s leads would read moved; monteCarloRate, P(|mean| > 2 sigma) for normal differences at the early RMS with lag1; inSampleRate, the share of the late pairs past 2 sigma. proxy: the twins carry no gridBoxSkin or daily sea-ice line, so Ts takes surfaceT and sea ice the 30-day-mean cut on the index\'s grids',
      twins: twins.map((t) => ({ label: t.label, start: t.start, days: t.a.map((r) => r.day) })), from, measured: new Date().toISOString().slice(0, 10),
      expires: 'remeasure whenever a row changes its instrument or the control physics changes, and with more starts when they exist',
      script: 'scripts/scorecard64.mjs --sigma', nulls, rows,
    };
    writeFileSync(out, JSON.stringify(record, null, 1) + '\n');
    const f = (x, n = 3) => (Number.isFinite(x) ? x.toFixed(n) : '—');
    console.log('| row | σ (gate) | per window | lag1 | dof | early 30 / 60 | saturation | null at leads | MC rate |');
    console.log('|---|---|---|---|---|---|---|---|---|');
    for (const [key, name, units] of ROWS) { const r = rows[key]; console.log(`| ${name} (${units}) | ${f(r.sigma)}${r.proxy ? '*' : ''} | ${f(r.perWindow)} | ${f(r.lag1, 2)} | ${f(r.dof, 1)} | ${f(r.early[0])} / ${f(r.early[1])} | ${f(r.saturation, 2)} | ${r.nullAtLeads.map((x) => (x ? 'moved' : x === false ? '·' : '—')).join(' ')} | ${f(100 * r.monteCarloRate, 1)} % |`); }
    console.log(`nulls: ${nulls.atLeads.map((n) => `${n.start}: ${n.moved} of ${n.rows} moved`).join('; ')}; Monte Carlo ${(100 * nulls.monteCarloRate).toFixed(1)} %, in-sample ${(100 * nulls.inSampleRate).toFixed(1)} %`);
  } else {
    const jsonOut = opt('--json'), baseFile = opt('--base'), logFile = opt('--log');
    const reductions = [];
    for (const path of args) reductions.push(await reduceFile(path, logFile));
    const sigma = JSON.parse(readFileSync(SIGMA_FILE, 'utf8')).rows, means = windowMeans(reductions);
    const comparison = baseFile ? compareScores(JSON.parse(readFileSync(baseFile, 'utf8')).means, means, sigma) : null;
    if (jsonOut) writeFileSync(jsonOut, JSON.stringify({ states: args, reductions, means, comparison }, null, 1) + '\n');
    console.log(formatTable(reductions, sigma, comparison));
    if (comparison) console.log(`moved: ${comparison.filter((c) => c.verdict === 'moved').map((c) => c.key).join(', ') || 'none'}`);
  }
}
