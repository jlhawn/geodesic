// The sea surface temperature of saved states against WOA23's monthly
// climatology (data/woa_sst_monthly_1deg.bin), in the boxes of the thin
// mixed-layer judging. A state that carries a whole window of daily means
// (sstWindowOf in js/sstMeans.module.js) is judged by each cell's mean
// over the window's days against the climatology averaged over the same
// days, each day's mean dated to its middle; any other state, and every
// state under --instant, by the top layer's temperature at the instant it
// was saved against the climatology on its date.
//
//   node scripts/woaSstBias.mjs STATE.bin[@DAY] [STATE.bin[@DAY] ...] [--instant] [--woa FILE] [--json OUT]
//
// DAY is the model day (day 0 the March equinox, March 20) of an instant
// comparison, the state's own by default; a state given a DAY is compared
// at its instant. The climatology is interpolated linearly between the
// month centers to each date. With several states each cell's model and
// climatology values are averaged over them first. A sea cell holds more
// than 1 m of water; "coastal" is within two cells of land; the open-water
// box leaves out cells iced in any of the states.
import { readFileSync, writeFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { decodeState } from '../js/stateFile.module.js';
import { stateGrid } from '../js/stateGrid.module.js';
import { stateGeography } from '../js/geography.module.js';
import { decodeMonthlySst, monthWeights, climatologyAt, MONTHLY_SST_FILE } from '../js/ocean/climatology.module.js';
import { sstWindowOf } from '../js/sstMeans.module.js';

const EQUINOX_DAY_OF_YEAR = 78, KELVIN = 273.15, DEG = 180 / Math.PI;
export const dayOfYear = (modelDay) => (((modelDay + EQUINOX_DAY_OF_YEAR) % 365) + 365) % 365;

export function windowMonths(sst, days) {
  const total = new Map();
  for (const d of days) for (const [m, w] of monthWeights(sst, dayOfYear(d - 0.5))) total.set(m, (total.get(m) ?? 0) + w / days.length);
  return [...total];
}

/*
 * A state's SST as this table judges it: { kind: 'mean', days, model,
 * months } over its window, or { kind: 'instant', day, model, months },
 * `model` holding the cells' temperatures in K (for the instant the
 * state's ocean.T, the top layer first) and `months` the climatology's
 * month weights.
 */
export function stateSst(st, sst, { day, instant = false } = {}) {
  const window = instant || day !== undefined ? null : sstWindowOf(st);
  if (window) return { kind: 'mean', days: window.days, model: window.mean, months: windowMonths(sst, window.days) };
  const d = day ?? st.day;
  return { kind: 'instant', day: d, model: st.ocean.T, months: monthWeights(sst, dayOfYear(d)) };
}

const inLon = (lon, a, b) => (a <= b ? lon >= a && lon <= b : lon >= a || lon <= b);
export const BOXES = [
  { name: '60S–60N', lat: [-60, 60] },
  { name: '60–40S', lat: [-60, -40] },
  { name: '40–20S', lat: [-40, -20] },
  { name: '20S–20N', lat: [-20, 20] },
  { name: '20–40N', lat: [20, 40] },
  { name: '40–60N', lat: [40, 60] },
  { name: 'Cold tongue 5S–5N 130–90W', lat: [-5, 5], lon: [-130, -90] },
  { name: 'Cold tongue 2S–2N 110–90W', lat: [-2, 2], lon: [-110, -90] },
  { name: 'Warm pool 10S–10N 120–160E', lat: [-10, 10], lon: [120, 160] },
  { name: 'Equatorial W 5S–5N 140–170E', lat: [-5, 5], lon: [140, 170], west: true },
  { name: 'Equatorial E 5S–5N 120–90W', lat: [-5, 5], lon: [-120, -90], east: true },
  { name: 'Benguela 30–15S 8–20E, coastal', lat: [-30, -15], lon: [8, 20], coastal: true },
  { name: 'Peru 20–5S 85–70W, coastal', lat: [-20, -5], lon: [-85, -70], coastal: true },
  { name: 'California 30–42N 130–115W, coastal', lat: [30, 42], lon: [-130, -115], coastal: true },
  { name: 'Canary 15–30N 25–10W, coastal', lat: [15, 30], lon: [-25, -10], coastal: true },
  { name: 'Open water north of 66N', lat: [66, 90], open: true },
];

function coastalCells(mesh, sea, rings = 2) {
  const { nCells: C, nEdgesOnCell, cellsOnCell, maxEdges } = mesh;
  let front = [];
  const near = new Uint8Array(C);
  for (let i = 0; i < C; i++) if (!sea[i]) { near[i] = 1; front.push(i); }
  for (let r = 0; r < rings; r++) {
    const next = [];
    for (const i of front) for (let m = 0; m < nEdgesOnCell[i]; m++) { const j = cellsOnCell[maxEdges * i + m]; if (!near[j]) { near[j] = 1; next.push(j); } }
    front = next;
  }
  return near;
}

export async function biasTable(entries, sst, { instant = false } = {}) {
  let mesh = null, C = 0, model = null, clim = null, iced = null, sea = null;
  const dates = [];
  for (const { path, day } of entries) {
    const st = await decodeState(new Uint8Array(readFileSync(path)));
    const grid = stateGrid(st, stateGeography(st).topography).mesh;
    if (!mesh) { mesh = grid; C = mesh.nCells; model = new Float64Array(C); clim = new Float64Array(C); iced = new Float64Array(C); sea = new Uint8Array(C).fill(1); }
    else if (grid.nCells !== C) throw new Error(`${path} is on another grid`);
    const ocean = st.ocean, L = ocean.h.length / C, judged = stateSst(st, sst, { day, instant }), { months } = judged;
    dates.push(judged.kind === 'mean' ? { path, kind: 'mean', days: [judged.days[0], judged.days[judged.days.length - 1]], length: judged.days.length, months } : { path, kind: 'instant', day: judged.day, doy: dayOfYear(judged.day), months });
    for (let i = 0; i < C; i++) {
      let column = 0;
      for (let k = 0; k < L; k++) column += ocean.h[k * C + i];
      if (!(column > 1)) { sea[i] = 0; continue; }
      model[i] += judged.model[i] - KELVIN;
      iced[i] += st.ice[i] > 0 ? 1 : 0;
      clim[i] += climatologyAt(sst, months, mesh.latCell[i], mesh.lonCell[i]) - KELVIN;
    }
  }
  for (let i = 0; i < C; i++) { model[i] /= entries.length; clim[i] /= entries.length; iced[i] /= entries.length; }
  const coast = coastalCells(mesh, sea);
  const rows = BOXES.map((box) => {
    let area = 0, m = 0, w = 0, squared = 0, cells = 0, missing = 0, icedArea = 0;
    for (let i = 0; i < C; i++) {
      if (!sea[i]) continue;
      const lat = mesh.latCell[i] * DEG, lon = ((mesh.lonCell[i] * DEG + 540) % 360) - 180;
      if (lat < box.lat[0] || lat > box.lat[1] || (box.lon && !inLon(lon, box.lon[0], box.lon[1])) || (box.coastal && !coast[i]) || (box.open && iced[i] > 0)) continue;
      if (!Number.isFinite(clim[i])) { missing++; continue; }
      const a = mesh.areaCell[i];
      area += a; m += a * model[i]; w += a * clim[i]; squared += a * (model[i] - clim[i]) ** 2; cells++; icedArea += a * iced[i];
    }
    return { box: box.name, cells, missing, model: m / area, woa: w / area, bias: (m - w) / area, rms: Math.sqrt(squared / area), iced: icedArea / area };
  });
  const west = rows[BOXES.findIndex((b) => b.west)], east = rows[BOXES.findIndex((b) => b.east)];
  rows.splice(BOXES.findIndex((b) => b.east) + 1, 0, { box: 'W−E', cells: null, missing: 0, model: west.model - east.model, woa: west.woa - east.woa, bias: west.model - east.model - (west.woa - east.woa), rms: NaN, iced: NaN });
  let spacing = 0;
  for (let e = 0; e < mesh.nEdges; e++) spacing += mesh.dcEdge[e];
  return { dates, rows, cells: C, spacing: spacing / mesh.nEdges };
}

const MONTH_NAMES = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];
function dateLabel(doy) {
  const lengths = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
  let m = 0, d = Math.floor(doy);
  while (d >= lengths[m]) d -= lengths[m++];
  return `${MONTH_NAMES[m]} ${d + 1}`;
}

export function markdown({ dates, rows, cells, spacing }) {
  const out = [], weights = (months) => months.map(([m, w]) => `${MONTH_NAMES[m]} ×${w.toFixed(3)}`).join(' + ');
  for (const date of dates) {
    const name = date.path.split('/').pop(), [first, last] = date.days ?? [];
    out.push(date.kind === 'mean'
      ? `- ${name}: model days ${first}–${last} (${dateLabel(dayOfYear(first - 1))} – ${dateLabel(dayOfYear(last))}), the ${date.length}-day mean SST against WOA23 over the same days, ${weights(date.months)}`
      : `- ${name}: model day ${date.day}, ${dateLabel(date.doy)} (day of year ${date.doy.toFixed(1)}), the state's instant SST against WOA23 ${weights(date.months)}`);
  }
  out.push(`- ${cells} cells, mean spacing ${(spacing / 1000).toFixed(0)} km; "coastal" is within two cells of land.`);
  const kinds = new Set(dates.map((date) => date.kind)), notes = [];
  if (kinds.has('mean')) notes.push('Each cell\'s mean SST averages its daily means, each over every ocean step of its day, which removes the day\'s cycle and averages down the passing eddies and waves.');
  if (kinds.has('instant')) notes.push('An instant state is sampled at one UTC, so each longitude at one local time: read differences under 0.3 K between runs as noise.');
  out.push(`- ${notes.join(' ')}`);
  out.push('');
  out.push('| Box | Cells | Model (°C) | WOA23 (°C) | Bias (K) | RMS (K) | Iced |');
  out.push('|---|---|---|---|---|---|---|');
  const signed = (x) => `${x >= 0 ? '+' : ''}${x.toFixed(2)}`;
  for (const r of rows) {
    if (r.cells === null) { out.push(`| ${r.box} | | ${r.model.toFixed(2)} | ${r.woa.toFixed(2)} | ${signed(r.bias)} | | |`); continue; }
    out.push(`| ${r.box} | ${r.cells}${r.missing ? ` (+${r.missing} without WOA23)` : ''} | ${r.model.toFixed(2)} | ${r.woa.toFixed(2)} | ${signed(r.bias)} | ${r.rms.toFixed(2)} | ${(100 * r.iced).toFixed(0)} % |`);
  }
  return out.join('\n');
}

if (process.argv[1] && fileURLToPath(import.meta.url) === process.argv[1]) {
  const args = process.argv.slice(2), entries = [];
  let woaPath = fileURLToPath(new URL(`../${MONTHLY_SST_FILE}`, import.meta.url)), jsonOut = null, instant = false;
  for (let a = 0; a < args.length; a++) {
    if (args[a] === '--instant') instant = true;
    else if (args[a] === '--woa') woaPath = args[++a];
    else if (args[a] === '--json') jsonOut = args[++a];
    else { const [path, day] = args[a].split('@'); entries.push({ path, day: day === undefined ? undefined : Number(day) }); }
  }
  if (!entries.length) { console.error('usage: node scripts/woaSstBias.mjs STATE.bin[@DAY] [STATE.bin[@DAY] ...] [--instant] [--woa FILE] [--json OUT]'); process.exit(1); }
  const sst = decodeMonthlySst(readFileSync(woaPath));
  const table = await biasTable(entries, sst, { instant });
  console.log(markdown(table));
  if (jsonOut) writeFileSync(jsonOut, JSON.stringify(table, null, 1));
}
