// What the soil goal added and removed over a quarter, from the saves of a
// spin-up run with LAND '{"soilGoal":"observed"}' (scripts/spinup.mjs), as
// markdown tables of the open land by 20° band and by region: the root
// zone's and surface layer's water (kg/m², and Gt for the net), the cover
// and the trees' share of the observed cover goal's change beyond the
// model's own. Each save holds the sums of its window (its header's
// goalSince to its day); the quarter is the chain of closed windows ending
// at the given state, and the same quarter of each earlier year (the state
// 365 days before, 730, …) goes beside it: when the model's own water and
// vegetation approach the observed ones the increments shrink year on year.
//
//   node scripts/landIncrements.mjs runs/eighteen/spin128_day0365.bin [--days 91] [--years 3]
//
// Builds no model. Exits 3, printing nothing, when the state carries no
// increments or its window is open (a halted save); earlier years without
// a save are left out.
import { readdirSync, existsSync, realpathSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { basename, dirname, join } from 'node:path';
import { createGeography } from '../js/geography.module.js';
import { stateGrid } from '../js/stateGrid.module.js';
import { readState, readGeography } from './figures/figureState.mjs';

const DEG = 180 / Math.PI, YEAR = 365, FIELDS = ['goalAdded', 'goalRemoved', 'goalCover', 'goalTrees'];
const BANDS = [[-60, -40], [-40, -20], [-20, 0], [0, 20], [20, 40], [40, 60], [60, 80], [80, 90]];
export const REGIONS = {
  'Amazon': [-10, 3, -70, -50], 'Congo': [-5, 5, 12, 30], 'Sahel': [8, 16, -15, 35], 'India': [15, 28, 72, 88], 'Southeast Asia': [10, 25, 95, 110],
  'Europe': [45, 55, 0, 30], 'Eastern US': [32, 45, -95, -75], 'Siberia': [55, 65, 60, 120], 'Australia': [-37, -20, 115, 150],
};

export function parseArguments(argv) {
  const rest = [], options = { days: 91, years: 3 };
  for (let a = 0; a < argv.length; a++) {
    if (argv[a] === '--days' || argv[a] === '--years') options[argv[a].slice(2)] = Number(argv[++a]);
    else rest.push(argv[a]);
  }
  if (rest.length !== 1 || !(options.days > 0) || !Number.isInteger(options.years) || options.years < 0) throw new Error('usage: node scripts/landIncrements.mjs <state.bin> [--days 91] [--years 3]');
  return { file: rest[0], ...options };
}

const named = (file) => { const m = basename(file).match(/^(.*)_day(\d+)\.bin$/); return m ? { prefix: m[1], day: Number(m[2]) } : null; };

async function window(file) {
  const saved = await readState(file);
  if (!saved.land?.goalAdded || saved.goalOpen !== false) return null;
  return { saved, since: saved.goalSince, day: saved.day, sums: Object.fromEntries(FIELDS.map((name) => [name, Float64Array.from(saved.land[name])])) };
}

/*
 * The closed windows that end at `file` and, walking back from each one's
 * goalSince to the save ending there, cover `days` days (one short of them
 * is enough), summed by field; null with no increments there.
 */
async function quarter(file, days, C) {
  const dir = dirname(file), { prefix } = named(file), saves = new Map();
  for (const name of readdirSync(dir)) { const n = named(name); if (n && n.prefix === prefix) saves.set(n.day, join(dir, name)); }
  const sums = Object.fromEntries(FIELDS.map((name) => [name, new Float64Array(C)]));
  let to = null, from = null, count = 0;
  for (let at = file; at;) {
    const w = await window(at);
    if (!w || w.saved.land.goalAdded.length !== C) break;
    to ??= w.day;
    from = w.since;
    count++;
    for (const name of FIELDS) for (let i = 0; i < C; i++) sums[name][i] += w.sums[name][i];
    if (to - from >= days - 1) break;
    at = saves.get(from);
  }
  return count ? { sums, from, to, saves: count } : null;
}

function summary(cells, area, sums) {
  let a = 0, net = 0, abs = 0, gt = 0, cover = 0, coverAbs = 0, trees = 0, treesAbs = 0, gross = 0;
  for (const i of cells) {
    const n = sums.goalAdded[i] - sums.goalRemoved[i];
    a += area[i]; net += area[i] * n; abs += area[i] * Math.abs(n); gt += area[i] * n; gross += area[i] * (sums.goalAdded[i] + sums.goalRemoved[i]);
    cover += area[i] * sums.goalCover[i]; coverAbs += area[i] * Math.abs(sums.goalCover[i]);
    trees += area[i] * sums.goalTrees[i]; treesAbs += area[i] * Math.abs(sums.goalTrees[i]);
  }
  const mean = (x) => (a > 0 ? x / a : NaN);
  return { cells: cells.length, water: mean(net), waterAbs: mean(abs), gross: mean(gross), gt: gt / 1e12, cover: mean(cover), coverAbs: mean(coverAbs), trees: mean(trees), treesAbs: mean(treesAbs) };
}

const sign = (x, d) => (Number.isNaN(x) ? '' : `${x >= 0 ? '+' : '−'}${Math.abs(x).toFixed(d)}`);
const plain = (x, d) => (Number.isNaN(x) ? '' : x.toFixed(d));

export function tables(rows, earlier, label) {
  const out = [];
  const quantity = (title, key, absKey, digits, extra) => {
    const head = ['Region', 'Cells', `Net${label}`, 'Mean abs', ...(extra ? [extra.title] : []), ...earlier.map((e) => `Net / abs, ${e.label}`)];
    out.push(`### ${title}`, '', `| ${head.join(' | ')} |`, `|${head.map((h, k) => (k ? '---:' : '---')).join('|')}|`);
    rows.forEach((row, r) => {
      out.push(`| ${[row.name, row.s.cells, sign(row.s[key], digits), plain(row.s[absKey], digits), ...(extra ? [sign(row.s[extra.key], 2)] : []), ...earlier.map((e) => `${sign(e.rows[r].s[key], digits)} / ${plain(e.rows[r].s[absKey], digits)}`)].join(' | ')} |`);
    });
    out.push('');
  };
  quantity('Water of the root zone and surface layer (kg/m²)', 'water', 'waterAbs', 1, { title: 'Net (Gt)', key: 'gt' });
  quantity('Vegetation cover (fraction)', 'cover', 'coverAbs', 4);
  quantity('Tree fraction', 'trees', 'treesAbs', 4);
  return out.join('\n');
}

async function main() {
  const { file, days, years } = parseArguments(process.argv.slice(2));
  if (!named(file)) throw new Error(`${file} is not a <tag>_dayDDDD.bin state`);
  if (!existsSync(file)) throw new Error(`${file} does not exist`);
  const saved = await readState(file);
  if (!saved.land?.goalAdded || saved.goalOpen !== false) process.exit(3);
  const ground = readGeography(saved), mesh = stateGrid(saved, ground.topography).mesh, C = mesh.nCells;
  const geo = createGeography(mesh, ground.topography, ground.geography);
  const open = [...Array(C).keys()].filter((i) => geo.land[i] && !(geo.iceSheet && geo.iceSheet[i] > 0));
  const at = (box) => open.filter((i) => { const lat = mesh.latCell[i] * DEG, lon = mesh.lonCell[i] * DEG, [a, b, c, e] = box; return lat >= a && lat <= b && (c <= e ? lon >= c && lon <= e : lon >= c || lon <= e); });
  const groups = [
    ...BANDS.map(([a, b]) => ({ name: `${Math.abs(a)}${a < 0 ? 'S' : a > 0 ? 'N' : ''}–${Math.abs(b)}${b < 0 ? 'S' : b > 0 ? 'N' : ''}`, cells: open.filter((i) => { const lat = mesh.latCell[i] * DEG; return lat >= a && (lat < b || b === 90); }) })),
    ...Object.entries(REGIONS).map(([name, box]) => ({ name, cells: at(box) })),
    { name: 'Open land', cells: open },
  ].filter((g) => g.cells.length);
  const rowsOf = (sums) => groups.map((g) => ({ name: g.name, s: summary(g.cells, mesh.areaCell, sums) }));
  const now = await quarter(file, days, C);
  const { day } = named(file), earlier = [];
  for (let k = 1; k <= years; k++) {
    const target = day - k * YEAR;
    if (target < 0) break;
    const dir = dirname(file), { prefix } = named(file);
    const near = readdirSync(dir).map((name) => ({ name, n: named(name) })).filter((x) => x.n && x.n.prefix === prefix && Math.abs(x.n.day - target) <= 3).sort((p, q) => Math.abs(p.n.day - target) - Math.abs(q.n.day - target))[0];
    const before = near ? await quarter(join(dir, near.name), days, C) : null;
    if (before) earlier.push({ label: `days ${before.from}–${before.to}`, rows: rowsOf(before.sums) });
  }
  console.log(`Increments of the soil goal over days ${now.from}–${now.to} (${now.saves} save${now.saves === 1 ? '' : 's'}), by 20° band and region of the open land (${open.length} cells). Net is the area mean of added minus removed per cell, mean abs the area mean of its absolute value, both over the whole window; a positive net means the goal added water or cover. Earlier years give the same quarter's net / mean abs.`);
  console.log();
  console.log(tables(rowsOf(now.sums), earlier, ''));
}

if (process.argv[1] && realpathSync(process.argv[1]) === fileURLToPath(import.meta.url)) await main();
