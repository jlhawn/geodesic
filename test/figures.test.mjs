import { test, before, after } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, readdirSync, existsSync, mkdtempSync, rmSync, copyFileSync, statSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join, dirname } from 'node:path';
import { execFileSync, spawnSync } from 'node:child_process';
import { season } from '../scripts/figures/figureState.mjs';

const root = new URL('..', import.meta.url).pathname;
let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const plotting = spawnSync('python3', ['-c', 'import matplotlib, numpy'], { encoding: 'utf8' }).status === 0;

function pulledState() {
  if (process.env.FIGURES_STATE) return existsSync(process.env.FIGURES_STATE) ? process.env.FIGURES_STATE : null;
  const dirs = [join(root, 'runs/verda-eleven')];
  try {
    const common = execFileSync('git', ['rev-parse', '--path-format=absolute', '--git-common-dir'], { cwd: root, encoding: 'utf8' }).trim();
    dirs.push(join(dirname(common), 'runs/verda-eleven'));
  } catch {}
  for (const dir of dirs) {
    if (!existsSync(dir)) continue;
    const states = readdirSync(dir).filter((f) => /^eleven64_day\d+\.bin$/.test(f)).sort();
    if (states.length) return join(dir, states[states.length - 1]);
  }
  return null;
}
const source = pulledState();
const noState = !source && 'no eleven64_day*.bin in runs/verda-eleven (or FIGURES_STATE) to draw from';
let dir, state;
before(() => {
  dir = mkdtempSync(join(tmpdir(), 'figures-'));
  if (source) { state = join(dir, 'eleven64_day0000.bin'); copyFileSync(source, state); }
});
after(() => rmSync(dir, { recursive: true, force: true }));

const dump = (figure, env = {}) => {
  const out = join(dir, `${figure}.json`);
  execFileSync('node', [join(root, 'scripts/figures', `${figure}.mjs`), state, out], { cwd: root, encoding: 'utf8', env: { ...process.env, OCEAN: '{"everySteps":8}', ...env } });
  return { out, d: JSON.parse(readFileSync(out, 'utf8')) };
};
const plot = (figure, json) => {
  if (!plotting) return '';
  const png = join(dir, `${figure}.png`);
  const printed = execFileSync('python3', [join(root, 'scripts/figures', `${figure}.py`), json, png], { cwd: root, encoding: 'utf8' });
  const bytes = readFileSync(png);
  assert.equal(bytes.subarray(1, 4).toString(), 'PNG', `${figure}.png is a PNG`);
  assert.ok(statSync(png).size > 20000, `${figure}.png holds a picture`);
  return printed;
};
const finite = (values) => values.filter((v) => v !== null && Number.isFinite(v));
const within = (values, lo, hi, what) => { const v = finite(values); assert.ok(v.length > 0, `${what} has values`); assert.ok(Math.min(...v) >= lo && Math.max(...v) <= hi, `${what} within [${lo}, ${hi}]: ${Math.min(...v)} to ${Math.max(...v)}`); };

test('the season of a model day counts from the March equinox at day 0', () => {
  assert.equal(season(0), 'northern spring, March equinox');
  assert.equal(season(274), 'northern winter, December solstice');
  assert.equal(season(101), 'northern summer, June solstice + 10 d');
  assert.equal(season(81), 'northern spring, June solstice − 10 d');
  assert.equal(season(365 + 183), 'northern autumn, September equinox, year 2');
});

test('the state maps dump every cell\'s sea, ice, land and surface fields and plot four panels', { skip: noState }, () => {
  const { out, d } = dump('stateMaps');
  assert.equal(d.N, 64); assert.equal(d.tag, 'eleven64'); assert.match(d.season, /northern|March/);
  for (const key of ['lat', 'lon', 'land', 'iceSheet', 'ice', 'concentration', 'snow', 'sst', 'ts', 'vegetation', 'trees']) assert.equal(d[key].length, 40962, key);
  const sea = d.sst.filter((v, i) => !d.land[i]);
  within(sea, -2.5, 36, 'mixed-layer temperature °C');
  assert.ok(d.sst.every((v, i) => !d.land[i] || v === null), 'no mixed-layer temperature on land');
  within(d.ts, -90, 60, 'surface temperature °C');
  within(d.ice, 0, 15, 'sea-ice thickness m');
  within(d.vegetation, 0, 1, 'vegetation'); within(d.trees, 0, 1, 'trees');
  assert.ok(d.trees.every((t, i) => t <= d.vegetation[i] + 1e-3), 'the trees are part of the cover');
  plot('stateMaps', out);
});

test('the equatorial section averages 2S–2N into 2° bins down to 300 m', { skip: noState }, () => {
  const { out, d } = dump('eqsection');
  assert.equal(d.layers, 45);
  assert.equal(d.lons.length, 81); assert.equal(d.depths.length, 61); assert.equal(d.T.length, 81);
  for (const col of d.T) assert.equal(col.length, 61);
  const surface = finite(d.T.map((col) => col[0]));
  assert.ok(surface.length > 70, 'nearly every bin is sea at the surface');
  within(d.T.flat(), 0, 33, 'section temperature °C');
  const west = d.T[d.lons.indexOf(160)], east = d.T[d.lons.indexOf(260)];
  assert.ok(west[0] > east[0], 'the west Pacific is warmer at the surface than the east');
  assert.ok(west[0] - west[60] > 5, 'the water cools with depth');
  plot('eqsection', out);
});

test('the equatorial panels read the surface and mixed-layer frames after four GPU steps', { skip: noState || (!gpuAvailable && 'webgpu not installed') }, () => {
  const { out, d } = dump('eqpanels');
  assert.deepEqual(d.columns, ['lon', 'lat', 'land', 'u10', 'v10', 'tair', 'mslp', 'cu', 'cv', 'sst', 'thermocline']);
  assert.equal(d.rows.length, d.polys.length);
  assert.ok(d.rows.length > 5000);
  const col = (name) => d.rows.map((r) => r[d.columns.indexOf(name)]);
  within(col('lon'), 89, 291, 'longitude'); within(col('lat'), -16, 16, 'latitude');
  within(col('u10'), -40, 40, 'eastward wind'); within(col('tair'), -10, 40, 'air temperature');
  within(col('mslp'), 960, 1060, 'sea-level pressure hPa'); within(col('sst'), 10, 35, 'mixed-layer temperature');
  within(col('cu'), -3, 3, 'eastward current'); within(col('thermocline'), 0, 6000, 'class top m');
  plot('eqpanels', out);
});

test('the deck dump matches the host port of the column on the night side and plots its maps and box statistics', { skip: noState || (!gpuAvailable && 'webgpu not installed') }, () => {
  const { out, d } = dump('mlmdeck');
  assert.equal(d.rows.length, 40962);
  for (const row of d.rows) assert.equal(row.length, d.columns.length);
  const col = (name) => d.rows.map((r) => r[d.columns.indexOf(name)]);
  within(col('cover'), 0, 1, 'cover'); within(col('lwp'), 0, 5000, 'LWP g/m²'); within(col('h'), 0, 3000, 'inversion height m');
  within(col('mlmGate'), 0, 1, 'gate memory'); within(col('gate'), 0, d.gates.length - 1, 'gate code');
  assert.ok(col('active').filter((a) => a === 1).length > 500, 'the deck runs somewhere');
  const [count, lwp, , all] = d.validation;
  const both = Number(count.match(/both (\d+)/)[1]);
  assert.ok(both > 500 && /only GPU 0, only host 0/.test(count), count);
  assert.ok(Number(lwp.match(/LWP ≥ 1 g\/m² relative max ([\d.e+-]+)/)[1]) < 1e-2, lwp);
  assert.ok(/decisions \(G > 0\.5\) differing 0/.test(all), all);
  const printed = plot('mlmdeck', out);
  if (plotting) assert.match(printed, /^deck boxes: SE Pacific/m);
});

test('the driver carries on past a failing figure and exits 1', () => {
  const run = spawnSync('bash', [join(root, 'scripts/figures/snapshot.sh'), join(dir, 'missing_day0007.bin'), join(dir, 'driver')], { cwd: root, encoding: 'utf8', env: { ...process.env, GPULOCK: '/nonexistent' } });
  assert.equal(run.status, 1);
  for (const figure of ['stateMaps', 'eqsection', 'eqpanels', 'mlmdeck']) assert.match(run.stdout, new RegExp(`^${figure}: the dump failed .*missing_${figure}_day0007\\.log`, 'm'));
});
