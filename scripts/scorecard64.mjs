#!/usr/bin/env node
/*
 * The N=64 observational scorecard, the per-merge gate's second half (the first is scripts/physicsDigest.mjs). Each
 * row is a 30-day day-mean window ending at a save, scored against its definition-matched reference in
 * data/references/ (README.md there gives each one's definition, period and box; scripts/references.mjs reads them).
 *
 * The gate run: the standard start, seventeen64 day 365 (the March equinox, the page's N=64 default state), run on the merge's tree with
 *   LEVELS=bl48 N=64 TAG=p64 OUT=DIR DAYS=395 node scripts/spinup.mjs, then DAYS=425 (about 6 GPU minutes)
 * and scored with
 *   node scripts/scorecard64.mjs [--json OUT.json] [--base BASE.json] DIR/p64_day0395.bin DIR/p64_day0425.bin
 * Each state's daily lines (TOA fluxes and the sea-ice cuts) are read from TAG.log beside it, or from --log FILE.
 * With --base, the rows whose window-mean model value moved by more than twice the pair sigma are listed as moved and
 * the rest as unresolved; a re-pin cites the moved rows.
 *
 * Sigma (data/scorecard64-sigma.json) is the pair sigma: the RMS over the two windows of the difference between twin
 * controls from the same start, the spread of one run against another. --sigma A_DIR B_DIR DAY... writes it from two
 * twins' saves, a row the twins cannot give (no gridBoxSkin or daily ice line) taking the named proxy.
 *
 * The 40–60N land row compares the lowest layer's air with ERA5's 2 m air until the model has a 2 m diagnostic; the
 * undercurrent rows read the window-end ocean, not a mean.
 */
import { readFileSync, writeFileSync, existsSync } from 'node:fs';
import { basename, dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { decodeState } from '../js/stateFile.module.js';
import { stateGeography, createGeography } from '../js/geography.module.js';
import { stateGrid } from '../js/stateGrid.module.js';
import { climatologyAt } from '../js/ocean/climatology.module.js';
import { cellVector } from '../js/dynamics/operators.module.js';
import { windowMonths } from './woaSstBias.mjs';
import * as R from './references.mjs';

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
  ['ts', 'global Ts, grid-box skin', '°C'], ['landT', '40–60N land air: lowest layer against 2 m', 'K'],
  ...Object.keys(BOXES).map((b) => [`rain_${b}`, `rain ${b}`, 'mm/d']),
  ['stress', 'eastward stress 160–100W within 3°', 'mN/m²'], ['slp', 'SLP 180–160W less 120–100W within 3°', 'hPa'],
  ['asr', 'ASR', 'W/m²'], ['net', 'net TOA', 'W/m²'], ['swcre', 'SWCRE', 'W/m²'], ['lwcre', 'LWCRE', 'W/m²'],
  ['iceNExtent', 'Arctic extent, daily cut', 'Mkm²'], ['iceNArea', 'Arctic area, daily cut', 'Mkm²'],
  ['iceSExtent', 'Antarctic extent, daily cut', 'Mkm²'], ['iceSArea', 'Antarctic area, daily cut', 'Mkm²'],
  ...R.TAO_STATIONS.flatMap((s) => [[`euc_${s}`, `undercurrent core speed ${s}`, 'm/s'], [`eucDepth_${s}`, `undercurrent core depth ${s}`, 'm']]),
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
    m = line.match(/sea ice on day (\d+)'s mean concentration, sea cells, 15 %: N extent ([\d.]+) area ([\d.]+), S extent ([\d.]+) area ([\d.]+)/);
    if (m) ice.set(+m[1], { nExtent: +m[2], nArea: +m[3], sExtent: +m[4], sArea: +m[5] });
  }
  return { days, ice };
}

export function undercurrentCore(profile) {
  const p = profile.filter((x) => Number.isFinite(x.u));
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
  if (!meshes.has(key)) { const { mesh } = stateGrid(saved, topography); meshes.set(key, { mesh, land: createGeography(mesh, topography, geography).land }); }
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
  const meanCut = (north) => {
    let e = 0, a = 0;
    for (let i = 0; i < C; i++) { const c = M.concentration[i]; if (land[i] || !sea[i] || !(c >= 0.15) || (lat(i) >= 0) !== north) continue; e += area[i]; a += c * area[i]; }
    return { extent: e / 1e12, area: a / 1e12 };
  };
  const [cutN, cutS] = [meanCut(true), meanCut(false)], nsidc = (k, h) => dayWeighted(D, (m) => iceRef[k][h][m - 1]);
  row('iceNExtent', fromLog(ice, (x) => x.nExtent), nsidc('extent', 'N'), { proxy: cutN.extent });
  row('iceNArea', fromLog(ice, (x) => x.nArea), nsidc('area', 'N'), { proxy: cutN.area });
  row('iceSExtent', fromLog(ice, (x) => x.sExtent), nsidc('extent', 'S'), { proxy: cutS.extent });
  row('iceSArea', fromLog(ice, (x) => x.sArea), nsidc('area', 'S'), { proxy: cutS.area });

  const eastward = Array.from({ length: L }, (_, k) => cellVector(mesh, Float64Array.from(saved.ocean.u.subarray(k * E, (k + 1) * E))));
  const Z = Array.from({ length: 67 }, (_, j) => 20 + 5 * j);
  for (const s of R.TAO_STATIONS) {
    const st = tao.stations[s], cs = cells((i) => sea[i] && Math.abs(lat(i)) <= 1 && Math.abs(wrap(lon(i) - st.longitude)) <= 1.5);
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

export function pairSigma(a, b, { proxy = true } = {}) {
  const out = {};
  for (const [key] of ROWS) {
    const pick = (r) => { const x = r.rows[key]; return Number.isFinite(x.model) ? x.model : proxy && Number.isFinite(x.proxy) ? x.proxy : NaN; };
    const d = a.map((r, k) => pick(r) - pick(b[k]));
    out[key] = { sigma: Math.sqrt(avg(d.map((x) => x * x))), proxy: a.some((r) => !Number.isFinite(r.rows[key].model)) };
  }
  return out;
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
  if (args[0] === '--sigma') {
    const [, a, b, ...days] = args, tagOf = (dir) => (existsSync(join(dir, 'p64.log')) ? 'p64' : basename(dir));
    const reduce = (dir) => Promise.all(days.map((d) => reduceFile(join(dir, `${tagOf(dir)}_day${String(d).padStart(4, '0')}.bin`), join(dir, `${tagOf(dir)}.log`))));
    const [ra, rb] = [await reduce(a), await reduce(b)], sigma = pairSigma(ra, rb);
    const record = {
      name: 'Pair sigma of the N=64 scorecard rows',
      definition: 'per row, the RMS over the windows of the difference between two twin controls from one start, so a change between two runs moves a row only past 2 sigma; proxy: rows the twins could not give (no gridBoxSkin means, no daily sea-ice line) take surfaceT for Ts and the 30-day-mean concentration cut for sea ice',
      twins: [a, b].map((dir) => dir.replace(process.env.HOME, '~')), days: days.map(Number), n: days.length, measured: new Date().toISOString().slice(0, 10),
      expires: 'remeasure with new twins whenever a row changes its instrument or the pairs are rerun',
      script: 'scripts/scorecard64.mjs --sigma', rows: sigma,
    };
    writeFileSync(opt('--out') ?? SIGMA_FILE, JSON.stringify(record, null, 1) + '\n');
    console.log(formatTable(ra, sigma));
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
