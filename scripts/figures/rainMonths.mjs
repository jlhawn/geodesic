// The per-cell rain of a quarter's three saves against GPCP, for rainMonths.py: each save's 30-day
// mean precipitation as a total in mm, GPCP v2.3 (1991-2020 monthly climatology) over the same window
// weighted by the days of each calendar month, and the area-weighted global, land and sea means.
//   node scripts/figures/rainMonths.mjs <out.json> <save1.bin> <save2.bin> <save3.bin>
//   node scripts/figures/rainMonths.mjs --saves <report number> [reports a year]
import { writeFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { calendarDate, dayWeighted } from '../scorecard64.mjs';

export const SAVES_PER_YEAR = 12;
export const saveDay = (m) => Math.floor((m * 365 + 6) / SAVES_PER_YEAR);

export function quarterSaves(k, perYear = 4) {
  const n = SAVES_PER_YEAR / perYear;
  if (!Number.isInteger(n) || !Number.isInteger(k) || k < 1) throw new Error(`no saves for report ${k} at ${perYear} a year`);
  return Array.from({ length: n }, (_, j) => saveDay((k - 1) * n + j + 1));
}

export const monthWeights = (D) => Array.from({ length: 12 }, (_, m) => dayWeighted(D, (month) => (month === m + 1 ? 1 : 0)));

const NAMES = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];
const dateAt = (t) => { const d = calendarDate(t); return `${NAMES[d.getUTCMonth()]} ${d.getUTCDate()}`; };
export const windowDates = (D, days = 30) => `${dateAt(D - days)} – ${dateAt(D)}`;
export const weightLabel = (weights) => weights.map((w, m) => [w, m]).filter(([w]) => w > 1e-9).sort((a, b) => a[1] - b[1]).map(([w, m]) => `${NAMES[m]} ${w.toFixed(2)}`).join(', ');

export function areaStats(area, land, model, obs) {
  const parts = { global: () => true, land: (i) => land[i], sea: (i) => !land[i] };
  const out = {};
  for (const [name, take] of Object.entries(parts)) {
    let a = 0, m = 0, o = 0, d = 0;
    for (let i = 0; i < area.length; i++) if (take(i)) { a += area[i]; m += area[i] * model[i]; o += area[i] * obs[i]; d += area[i] * Math.abs(model[i] - obs[i]); }
    out[name] = { model: m / a, obs: o / a, diff: (m - o) / a, absDiff: d / a };
  }
  return out;
}

async function main(args) {
  if (args[0] === '--saves') { console.log(quarterSaves(Number(args[1]), args[2] ? Number(args[2]) : 4).join(' ')); return; }
  const [out, ...files] = args;
  if (!out || files.length === 0) { console.error('usage: node scripts/figures/rainMonths.mjs <out.json> <save.bin>...'); process.exit(2); }
  const { refuseUnreadOptions } = await import('../optionEnvironment.mjs');
  refuseUnreadOptions(import.meta.url);
  const { createGeography, stateGeography } = await import('../../js/geography.module.js');
  const { stateGrid } = await import('../../js/stateGrid.module.js');
  const { decodeState, stateName } = await import('../../js/stateFile.module.js');
  const { readFileSync } = await import('node:fs');
  const { basename } = await import('node:path');
  const R = await import('../references.mjs');
  const grid = R.gpcp().grid, DEG = 180 / Math.PI, rounded = (x, n) => Math.round(x * 10 ** n) / 10 ** n;
  let cells = null, tag = null, N = null, saves = [];
  for (const file of files) {
    const saved = await decodeState(new Uint8Array(readFileSync(file)));
    const days = saved.means?.days, rate = saved.means?.precipitation;
    if (!rate) throw new Error(`${file} holds no means.precipitation`);
    if (!cells) {
      const { topography, geography } = stateGeography(saved), mesh = stateGrid(saved, topography).mesh, geo = createGeography(mesh, topography, geography);
      cells = { mesh, land: Uint8Array.from(geo.land, (v) => (v ? 1 : 0)), lat: Float64Array.from(mesh.latCell, (v) => v * DEG), lon: Float64Array.from(mesh.lonCell, (v) => v * DEG) };
      tag = basename(stateName(file)).replace(/_day\d+(?:_step\d+)?$/, ''); N = saved.N;
    }
    const D = Math.round(saved.day);
    const named = /_day(\d+)/.exec(basename(file));
    if (named && Number(named[1]) !== D) throw new Error(`${file} holds day ${D}`);
    const weights = monthWeights(D), C = cells.mesh.nCells;
    const model = new Float64Array(C), obs = new Float64Array(C);
    for (let m = 0; m < 12; m++) {
      if (weights[m] < 1e-9) continue;
      for (let i = 0; i < C; i++) obs[i] += weights[m] * R.sampleGrid(grid, m, cells.lat[i], cells.lon[i]);
    }
    for (let i = 0; i < C; i++) { model[i] = rate[i] * days; obs[i] *= days; }
    saves.push({ day: D, days, window: windowDates(D, days), weights: weightLabel(weights), stats: areaStats(cells.mesh.areaCell, cells.land, model, obs), model: Array.from(model, (v) => rounded(v, 1)), obs: Array.from(obs, (v) => rounded(v, 1)) });
  }
  writeFileSync(out, JSON.stringify({ tag, N, day: saves[saves.length - 1].day, lat: Array.from(cells.lat, (v) => rounded(v, 2)), lon: Array.from(cells.lon, (v) => rounded(v, 2)), land: Array.from(cells.land), saves }));
  console.log(`wrote ${out} (${cells.mesh.nCells} cells, ${saves.length} saves)`);
}

if (process.argv[1] && fileURLToPath(import.meta.url) === (await import('node:fs')).realpathSync(process.argv[1])) await main(process.argv.slice(2));
