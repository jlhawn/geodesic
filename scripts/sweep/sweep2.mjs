// The second sweep's screens:
//   node scripts/sweep/sweep2.mjs
// writes a maximin Latin hypercube of POINTS (40) points over PARAMETERS2
// (scripts/sweep/runs.mjs) to runs/sweep2/design.csv (seeded by SEED, kept
// once written), with the defaults as point 0, and for every point not yet
// in runs/sweep2/results.csv runs its three N=64 screens on the GPU: TAG
// s2<NN>f, ten days from a fresh atlas start on bl34 (the means of days
// 6-10), beside s2<NN>m, three days from eight64_day0183 (day 186's
// balance and clear-sky albedo, the equatorial stress and the audit of the
// day-186 state over three windows two hours apart), followed by s2<NN>a, three days from
// nine64_day0091 (the 60-90N ice loss). PARALLEL (1) points run at once.
// The audits run on the CPU beside the next point's screens, and one row
// per point goes to results.csv with each parameter, each TERMS2 value,
// error and part, the score and the extra readings. ONLY (comma-separated
// point numbers) limits the screens to those points.
import { existsSync, readFileSync, writeFileSync, appendFileSync, mkdirSync } from 'node:fs';
import { readDesign } from './design.mjs';
import { PARAMETERS2, SWEEP2, STATES, spinup, auditWindows, arcticStart } from './runs.mjs';
import { TERMS2, score, readLog, dayMeans, arcticVolume } from './score.mjs';

const POINTS = Number(process.env.POINTS ?? 40), SEED = Number(process.env.SEED ?? 20261002), PARALLEL = Number(process.env.PARALLEL ?? 1);
const ONLY = process.env.ONLY ? process.env.ONLY.split(',').map(Number) : null;

export const EXTRA2 = ['fAsr', 'fNan', 'fClamped', 'eAlbedo', 'eOlr', 'eSwcre', 'eLwcre', 'eDayRain', 'rain', 'evaporation', 'peruRain', 'sepLwp', 'peruLwp', 'sepRuns', 'peruRuns', 'namibiaLow', 'sepInversion', 'zonalPeak', 'iceStart', 'iceEnd', 'clamped', 'nan'];

export async function screens2(point, prefix) {
  const fresh = spinup({ tag: `${prefix}f`, n: 64, days: 10, options: point, parameters: PARAMETERS2, env: { LEVELS: 'bl34' } });
  const eight = await spinup({ tag: `${prefix}m`, n: 64, from: `${STATES}/eight64_day0183.bin`, fromDay: 183, days: 3, options: point, parameters: PARAMETERS2 });
  const arctic = await spinup({ tag: `${prefix}a`, n: 64, from: `${STATES}/nine64_day0091.bin`, fromDay: 91, days: 3, options: point, parameters: PARAMETERS2 });
  return { fresh: await fresh, eight, arctic };
}

export async function values2(point, prefix, { fresh, eight, arctic }) {
  const flog = readLog(fresh.log), f = dayMeans(flog, 6, 10);
  const log = readLog(eight.log), day = log.days.find((d) => d.day === 186);
  const values = {
    fBalance: f.balance, fAlbedo: f.albedo, fOlr: f.olr, fSwcre: f.swcre, fLwcre: f.lwcre, fRain: f.rain, fAsr: f.asr,
    fNan: flog.nan ? 1 : 0, fClamped: flog.days.reduce((s, d) => s + d.clamped, 0),
    nan: log.nan || !day ? 1 : 0, clamped: log.days.reduce((s, d) => s + d.clamped, 0), stress: log.stress ?? NaN,
  };
  if (day && !log.nan) Object.assign(values, { balance: day.asr - day.olr, clearAlbedo: day.clearAlbedo, eAlbedo: day.albedo, eOlr: day.olr, eSwcre: day.swcre, eLwcre: day.lwcre, eDayRain: day.precip });
  if (existsSync(eight.state) && !log.nan) Object.assign(values, await auditWindows(eight.state, point, `${SWEEP2}/${prefix}m`, 64));
  const iceStart = await arcticStart();
  if (existsSync(arctic.state) && !readLog(arctic.log).nan) { const iceEnd = await arcticVolume(arctic.state); Object.assign(values, { iceStart, iceEnd, arctic: (iceStart - iceEnd) / 3 }); }
  return values;
}

const resultsFile = `${SWEEP2}/results.csv`;
const header = ['point', ...PARAMETERS2.map((p) => p.key), ...TERMS2.map((t) => t.key), ...TERMS2.map((t) => `e_${t.key}`), ...TERMS2.map((t) => `w_${t.key}`), 'score', ...EXTRA2];
const cell = (x) => (Number.isFinite(x) ? Number(x.toPrecision(6)) : '');

function row(point, values) {
  const s = score(values, TERMS2);
  return [point.point, ...PARAMETERS2.map((p) => point[p.key]), ...TERMS2.map((t) => cell(values[t.key])), ...TERMS2.map((t) => s.errors[t.key].toFixed(4)), ...TERMS2.map((t) => s.parts[t.key].toFixed(4)), s.total.toFixed(4), ...EXTRA2.map((k) => cell(values[k]))].join(',');
}

if (import.meta.url === `file://${process.argv[1]}`) {
  mkdirSync(SWEEP2, { recursive: true });
  const design = readDesign(`${SWEEP2}/design.csv`, PARAMETERS2, POINTS, SEED);
  if (!existsSync(resultsFile)) writeFileSync(resultsFile, header.join(',') + '\n');
  const done = new Set(readFileSync(resultsFile, 'utf8').trim().split('\n').slice(1).map((l) => Number(l.split(',')[0])));
  const queue = design.filter((p) => !done.has(p.point) && (!ONLY || ONLY.includes(p.point)));
  const pending = [];
  const worker = async () => {
    for (let point = queue.shift(); point; point = queue.shift()) {
      const prefix = `s2${String(point.point).padStart(2, '0')}`, t0 = Date.now();
      const runs = await screens2(point, prefix);
      console.log(`point ${point.point}: screens done in ${((Date.now() - t0) / 60000).toFixed(1)} min${runs.fresh.code || runs.eight.code || runs.arctic.code ? ' (NaN)' : ''}`);
      pending.push(values2(point, prefix, runs).then((values) => {
        appendFileSync(resultsFile, row(point, values) + '\n');
        console.log(`point ${point.point}: score ${score(values, TERMS2).total.toFixed(2)}`);
      }));
    }
  };
  await Promise.all(Array.from({ length: PARALLEL }, worker));
  await Promise.all(pending);
}
