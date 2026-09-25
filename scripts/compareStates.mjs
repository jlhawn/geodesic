// Saved states side by side, as a markdown table on stdout: for each
// state the means over its spin-up log's preceding window (Ts, the TOA
// imbalance, precipitation, ice, and the ocean's fastest current, largest
// transport and clamped edges), then fields measured from the state
// itself (SSTs, thermocline depths, surface winds, and the thermocline's
// cell-to-cell roughness in western boundary currents against the open
// ocean). With two states a third column gives the second minus the first.
//
// node scripts/compareStates.mjs runs/twin64_day0090.bin runs/twin128_day0090.bin [--window 90]
import { readFileSync, existsSync } from 'node:fs';
import { basename, dirname, join } from 'node:path';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { decodeState } from '../js/stateFile.module.js';
import { cellVector } from '../js/dynamics/operators.module.js';
import { THERMOCLINE_DENSITY, LAYER_DENSITIES } from '../js/ocean/layered.module.js';

const args = process.argv.slice(2);
const at = args.indexOf('--window');
const WINDOW = at >= 0 ? Number(args.splice(at, 2)[1]) : 90;
const files = args;
if (!files.length) { console.error('usage: node scripts/compareStates.mjs <state> [<state> ...] [--window days]'); process.exit(1); }

const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const meshes = new Map();
const modelFor = (N) => meshes.get(N) ?? meshes.set(N, createModel(new Grid(N), { physics: false, topography })).get(N);

/*
 * Means (maxima for the ocean's extremes) of the daily lines that
 * scripts/spinup.mjs wrote for days after day − window up to day.
 */
function logWindow(file, day) {
  const log = join(dirname(file), `${basename(file).replace(/_day\d+\.bin$/, '')}.log`);
  if (!existsSync(log)) return null;
  const pick = (line, pattern) => Number(line.match(pattern)?.[1]);
  const days = new Map();
  for (const line of readFileSync(log, 'utf8').split('\n')) {
    const d = pick(line, /^day (\d+) /);
    if (!(d > day - WINDOW && d <= day)) continue;
    days.set(d, {
      ts: pick(line, /Ts (-?[\d.]+) °C/), imbalance: pick(line, /ASR ([\d.]+)/) - pick(line, /OLR ([\d.]+)/),
      precip: pick(line, /precip ([\d.]+)/), ice: pick(line, /ice ([\d.]+)%/),
      speed: pick(line, /currents ≤ ([\d.]+)/), transport: pick(line, /transport ([\d.]+) Sv/), clamped: pick(line, /clamped (\d+)/),
    });
  }
  const rows = [...days.values()];
  if (!rows.length) return null;
  const mean = (key) => rows.reduce((s, r) => s + r[key], 0) / rows.length, max = (key) => Math.max(...rows.map((r) => r[key]));
  return { days: rows.length, ts: mean('ts'), imbalance: mean('imbalance'), precip: mean('precip'), ice: mean('ice'), speed: max('speed'), transport: max('transport'), clamped: max('clamped') };
}

async function measure(file) {
  const saved = await decodeState(new Uint8Array(readFileSync(file)));
  const model = modelFor(saved.N);
  const { mesh, core } = model, { K, C, E } = core.diagnostics, land = model.geography.land;
  const deg = 180 / Math.PI, lat = Float64Array.from(mesh.latCell, (x) => x * deg), lon = Float64Array.from(mesh.lonCell, (x) => x * deg);
  const inLon = (i, west, east) => (west <= east ? lon[i] >= west && lon[i] <= east : lon[i] >= west || lon[i] <= east);
  const box = (south, north, west, east, sea = true) => (i) => (!sea || !land[i]) && lat[i] >= south && lat[i] <= north && inLon(i, west, east);
  const mean = (values, where) => { let s = 0, w = 0; for (let i = 0; i < C; i++) if (where(i)) { s += mesh.areaCell[i] * values(i); w += mesh.areaCell[i]; } return w ? s / w : NaN; };

  const wind = cellVector(mesh, Float64Array.from(saved.u.subarray((K - 1) * E, K * E)), new Float64Array(3 * C));
  const east = (i) => -Math.sin(mesh.lonCell[i]) * wind[3 * i] + Math.cos(mesh.lonCell[i]) * wind[3 * i + 1];
  const band = (south, north) => mean(east, box(south, north, -180, 180, false));
  const strongest = (from, to) => {
    let best = -Infinity, where = 0;
    for (let s = from; s < to; s += 5) { const u = band(s, s + 5); if (u > best) { best = u; where = s + 2.5; } }
    return `${best.toFixed(1)} at ${Math.abs(where)}°${where < 0 ? 'S' : 'N'}`;
  };

  const h = saved.ocean.h, T = saved.ocean.T;
  const layers = LAYER_DENSITIES.filter((rho) => rho < THERMOCLINE_DENSITY).length;
  const thermocline = (i) => { let d = 0; for (let k = 0; k <= layers; k++) d += h[k * C + i]; return d; };
  const sst = (i) => T[i] - 273.15;
  const roughness = (where) => mean((i) => {
    let s = 0, n = 0;
    for (let k = 0; k < mesh.nEdgesOnCell[i]; k++) { const j = mesh.cellsOnCell[6 * i + k]; if (j >= 0 && !land[j]) { s += thermocline(j); n++; } }
    return n ? Math.abs(thermocline(i) - s / n) : 0;
  }, where);

  const warmPool = box(-10, 10, 130, 160);
  return {
    day: saved.day, N: saved.N, log: logWindow(file, saved.day),
    ts: mean((i) => saved.surfaceT[i] - 273.15, () => true),
    tropicalSst: mean(sst, box(-20, 20, -180, 180)), warmPoolSst: mean(sst, warmPool), coldTongueSst: mean(sst, box(-5, 5, -130, -90)),
    mixedLayer: mean((i) => h[i], box(-90, 90, -180, 180)),
    warmPoolD: mean(thermocline, warmPool), westD: mean(thermocline, box(-5, 5, 130, 150)), centralD: mean(thermocline, box(-5, 5, 170, -170)), eastD: mean(thermocline, box(-5, 5, -110, -90)),
    ridgeD: mean(thermocline, box(-10, -5, 50, 80)),
    nhTrades: band(10, 20), shTrades: band(-20, -10), nhWesterlies: strongest(30, 65), shWesterlies: strongest(-65, -30),
    pacificWind: mean(east, box(-5, 5, 150, -90, false)), indianWind: mean(east, box(-15, -5, 50, 90, false)),
    boundaryRough: [box(25, 40, 125, 150), box(-45, -30, 15, 40), box(30, 45, -80, -50)].reduce((sum, where) => sum + roughness(where), 0) / 3,
    openRough: roughness(box(20, 40, 180, -140)),
  };
}

const states = [];
for (const file of files) states.push({ file, ...(await measure(file)) });
const f = (digits) => (x) => (typeof x === 'number' && Number.isFinite(x) ? x.toFixed(digits) : x ?? '—');
const ROWS = [
  [`**Spin-up log, the ${WINDOW} days before**`, null],
  ['Days logged', (s) => s.log?.days, f(0)],
  ['Surface temperature, mean (°C)', (s) => s.log?.ts, f(2)],
  ['Absorbed solar − outgoing longwave (W/m²)', (s) => s.log?.imbalance, f(1)],
  ['Precipitation (mm/day)', (s) => s.log?.precip, f(2)],
  ['Sea ice, mean (% of the globe)', (s) => s.log?.ice, f(1)],
  ['Fastest current, max (m/s)', (s) => s.log?.speed, f(2)],
  ['Largest edge transport, max (Sv)', (s) => s.log?.transport, f(0)],
  ['Clamped edges, max', (s) => s.log?.clamped, f(0)],
  ['**The state on the day**', null],
  ['Surface temperature (°C)', (s) => s.ts, f(2)],
  ['SST 20°S–20°N (°C)', (s) => s.tropicalSst, f(2)],
  ['SST, west Pacific warm pool (°C)', (s) => s.warmPoolSst, f(2)],
  ['SST, east Pacific cold tongue (°C)', (s) => s.coldTongueSst, f(2)],
  ['Mixed layer, ocean mean (m)', (s) => s.mixedLayer, f(0)],
  ['Thermocline, warm pool (m)', (s) => s.warmPoolD, f(0)],
  ['Thermocline, equator 130–150°E / dateline / 110–90°W (m)', (s) => `${f(0)(s.westD)} / ${f(0)(s.centralD)} / ${f(0)(s.eastD)}`],
  ['Thermocline, Seychelles–Chagos 5–10°S 50–80°E (m)', (s) => s.ridgeD, f(0)],
  ['Thermocline roughness, boundary currents / open ocean (m)', (s) => `${f(1)(s.boundaryRough)} / ${f(1)(s.openRough)}`],
  ['Surface wind, trades 10–20°N / 10–20°S (m/s)', (s) => `${f(1)(s.nhTrades)} / ${f(1)(s.shTrades)}`],
  ['Surface wind, strongest westerlies N / S (m/s)', (s) => `${s.nhWesterlies} / ${s.shWesterlies}`],
  ['Surface wind, equatorial Pacific / Indian Ocean 5–15°S (m/s)', (s) => `${f(1)(s.pacificWind)} / ${f(1)(s.indianWind)}`],
];
const header = states.map((s) => `N=${s.N}`), pair = states.length === 2;
console.log(`\n### Day ${states.map((s) => s.day).join(' / ')}\n`);
console.log(`| | ${header.join(' | ')}${pair ? ` | N=${states[1].N} − N=${states[0].N}` : ''} |`);
console.log(`|---|${states.map(() => '---').join('|')}${pair ? '|---' : ''}|`);
for (const [name, value, format] of ROWS) {
  if (!value) { console.log(`| ${name} |${states.map(() => ' ').join('|')}${pair ? '| ' : ''}|`); continue; }
  const values = states.map(value);
  const cells = values.map((v) => (format ? format(v) : v ?? '—'));
  const difference = pair && typeof values[0] === 'number' && typeof values[1] === 'number' ? (values[1] - values[0] >= 0 ? '+' : '') + format(values[1] - values[0]) : pair ? '' : null;
  console.log(`| ${name} | ${cells.join(' | ')}${difference !== null ? ` | ${difference}` : ''} |`);
}
