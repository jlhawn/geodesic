// The perturbed-parameter sweep's screens:
//   node scripts/sweep/sweep.mjs
// writes a maximin Latin hypercube of POINTS (40) points over PARAMETERS
// (scripts/sweep/runs.mjs) to runs/sweep/design.csv (seeded by SEED, kept
// once written), with the defaults as point 0, and for every point not yet
// in runs/sweep/results.csv runs the two three-day N=64 screens one after
// the other on the GPU: TAG sx<NN>m from eight64_day0183 to day 186 and
// sx<NN>a from nine64_day0091 to day 94. The day-186 state is audited on
// the CPU beside the next point's screens, and one row per point goes to
// results.csv with each parameter, each score term's value, error and
// part, and the score (scripts/sweep/score.mjs). ONLY (comma-separated
// point numbers) limits the screens to those points.
import { existsSync, readFileSync, writeFileSync, appendFileSync, mkdirSync } from 'node:fs';
import { readDesign } from './design.mjs';
import { PARAMETERS, SWEEP, STATES, spinup, screenValues } from './runs.mjs';
import { TERMS, score } from './score.mjs';

const POINTS = Number(process.env.POINTS ?? 40), SEED = Number(process.env.SEED ?? 20261001);
const ONLY = process.env.ONLY ? process.env.ONLY.split(',').map(Number) : null;
mkdirSync(SWEEP, { recursive: true });

const designFile = `${SWEEP}/design.csv`;
const design = readDesign(designFile, PARAMETERS, POINTS, SEED);

const resultsFile = `${SWEEP}/results.csv`;
const EXTRA = ['dayRain', 'instantAlbedo', 'instantBalance', 'evaporation', 'sepRuns', 'peruRuns', 'namibiaLow', 'sepInversion', 'zonalPeakRain', 'iceStart', 'iceEnd', 'clamped', 'nan'];
const header = ['point', ...PARAMETERS.map((p) => p.key), ...TERMS.map((t) => t.key), ...TERMS.map((t) => `e_${t.key}`), ...TERMS.map((t) => `w_${t.key}`), 'score', ...EXTRA];
if (!existsSync(resultsFile)) writeFileSync(resultsFile, header.join(',') + '\n');
const done = new Set(readFileSync(resultsFile, 'utf8').trim().split('\n').slice(1).map((l) => Number(l.split(',')[0])));

function row(point, values) {
  const s = score(values);
  const cells = [point.point, ...PARAMETERS.map((p) => point[p.key]), ...TERMS.map((t) => values[t.key] ?? ''), ...TERMS.map((t) => s.errors[t.key]?.toFixed(4) ?? ''), ...TERMS.map((t) => s.parts[t.key]?.toFixed(4) ?? ''), s.total.toFixed(4), ...EXTRA.map((k) => values[k] ?? '')];
  return cells.join(',');
}

const pending = [];
for (const point of design) {
  if (done.has(point.point) || (ONLY && !ONLY.includes(point.point))) continue;
  const id = String(point.point).padStart(2, '0'), tagM = `sx${id}m`, tagA = `sx${id}a`, t0 = Date.now();
  const eight = await spinup({ tag: tagM, n: 64, from: `${STATES}/eight64_day0183.bin`, fromDay: 183, days: 3, options: point });
  const arctic = await spinup({ tag: tagA, n: 64, from: `${STATES}/nine64_day0091.bin`, fromDay: 91, days: 3, options: point });
  console.log(`point ${id}: screens done in ${((Date.now() - t0) / 60000).toFixed(1)} min${eight.code || arctic.code ? ' (NaN)' : ''}`);
  pending.push(screenValues(tagM, point, eight.log, eight.state, arctic.state).then((values) => {
    appendFileSync(resultsFile, row(point, values) + '\n');
    console.log(`point ${id}: score ${score(values).total.toFixed(2)}`);
  }));
}
await Promise.all(pending);
