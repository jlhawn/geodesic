// The sweep's parameters and the runs it makes: short GPU spin-ups from
// copies of saved states under the shared GPU lock (GPULOCK), each copy
// removed once the run has saved its last day, and the CPU audit of a state.
// A run whose last day is not yet saved starts over from its copy, any
// other snapshot of its tag removed first so that spinup.mjs cannot
// continue from it. PARAMETERS are the first sweep's, PARAMETERS2 the
// second's; an option whose default is null takes its base from `unset`,
// and a point's dragScale multiplies the sea's drag (the replicates).
import { spawn } from 'node:child_process';
import { copyFileSync, existsSync, unlinkSync, openSync, closeSync, readdirSync, renameSync } from 'node:fs';
import { readLog, readAudit, meanAudits, arcticVolume } from './score.mjs';
import { PHYSICS_DEFAULTS } from '../../js/gpu/core.gpu.js';
import { SEA_DRAG } from '../../js/physics/surface.module.js';

export const ROOT = new URL('../../', import.meta.url).pathname;
export const OUT = process.env.OUT ?? `${ROOT}runs`;
export const SWEEP = `${OUT}/sweep`;
export const SWEEP2 = `${OUT}/sweep2`;
export const STATES = process.env.STATES ?? '/Users/jlhawn/git_repos/jlhawn/geodesic/runs';
export const GPULOCK = process.env.GPULOCK ?? '/private/tmp/claude-501/-Users-jlhawn-git-repos-jlhawn-geodesic/4ab28226-4f31-4147-a831-f542e02bf0fc/scratchpad/gpulock.sh';

const RANGES = [
  { key: 'varianceScale', module: 'radiation', option: 'varianceScale', low: 2, high: 10 },
  { key: 'mixingLength', module: 'radiation', option: 'mixingLength', low: 150, high: 600 },
  { key: 'stratiformHours', module: 'moist', option: 'stratiformLifetime', low: 1, high: 6, scale: 3600 },
  { key: 'cloudHours', module: 'moist', option: 'cloudLifetime', low: 0.5, high: 2, scale: 3600 },
  { key: 'plumeEntrainment', module: 'moist', option: 'plumeEntrainment', low: 0.05, high: 0.2 },
  { key: 'plumeCape', module: 'moist', option: 'plumeCape', low: 40, high: 200 },
  { key: 'minimumInversion', module: 'radiation', option: 'minimumInversion', low: 2, high: 6 },
  { key: 'criticalHumidity', module: 'radiation', option: 'criticalHumidity', low: 0.7, high: 0.9 },
  { key: 'seaDrag', module: 'surface', option: 'dragCoefficient', low: 1.0e-3, high: 1.5e-3 },
  { key: 'stableMixingLength', module: 'radiation', option: 'stableMixingLength', low: 10, high: 60 },
  { key: 'cumulusCeiling', module: 'radiation', option: 'cumulusCeiling', low: 1500, high: 2500 },
];
const RANGES2 = [
  { key: 'upperHours', module: 'moist', option: 'upperCloudLifetime', low: 1, high: 8, scale: 3600, unset: 'cloudLifetime' },
  { key: 'stratiformHours', module: 'moist', option: 'stratiformLifetime', low: 2, high: 8, scale: 3600 },
  { key: 'cloudHours', module: 'moist', option: 'cloudLifetime', low: 0.5, high: 2, scale: 3600 },
  { key: 'varianceScale', module: 'radiation', option: 'varianceScale', low: 2, high: 10 },
  { key: 'plumeRainRate', module: 'moist', option: 'plumeRainRate', low: 1e-3, high: 6e-3 },
  { key: 'criticalHumidity', module: 'radiation', option: 'criticalHumidity', low: 0.7, high: 0.9 },
  { key: 'stratusWaterMax', module: 'radiation', option: 'stratusWaterMax', low: 0.1, high: 0.3 },
  { key: 'cumulusCeiling', module: 'radiation', option: 'cumulusCeiling', low: 1500, high: 2500 },
];
const baseOf = (p) => (p.module === 'surface' ? SEA_DRAG : PHYSICS_DEFAULTS[p.option] ?? PHYSICS_DEFAULTS[p.unset]) / (p.scale ?? 1);
export const PARAMETERS = RANGES.map((p) => ({ ...p, base: baseOf(p) }));
export const PARAMETERS2 = RANGES2.map((p) => ({ ...p, base: baseOf(p) }));

export function optionsOf(point, parameters = PARAMETERS) {
  const o = { radiation: {}, moist: {}, surface: {} };
  for (const p of parameters) o[p.module][p.option] = Number((point[p.key] * (p.scale ?? 1)).toPrecision(6));
  if (point.dragScale) o.surface.dragCoefficient = Number((SEA_DRAG * point.dragScale).toPrecision(8));
  return { RADIATION: JSON.stringify(o.radiation), MOIST: JSON.stringify(o.moist), SURFACE: JSON.stringify(o.surface) };
}

function run(command, args, env, outFile) {
  return new Promise((resolve, reject) => {
    const fd = outFile ? openSync(outFile, 'w') : 'ignore';
    const child = spawn(command, args, { cwd: ROOT, env: { ...process.env, ...env }, stdio: ['ignore', fd, fd] });
    child.on('exit', (code) => { if (outFile) closeSync(fd); code === 0 || code === 2 ? resolve(code) : reject(new Error(`${command} ${args.join(' ')} exited ${code}`)); });
  });
}

const dayName = (tag, day) => `${OUT}/${tag}_day${String(day).padStart(4, '0')}.bin`;

export async function spinup({ tag, n, from = null, fromDay = 0, days, options, parameters = PARAMETERS, env = {} }) {
  const end = fromDay + days, last = dayName(tag, end);
  if (existsSync(last)) return { log: `${OUT}/${tag}.log`, state: last, code: 0 };
  const start = from ? dayName(tag, fromDay) : null;
  if (existsSync(`${OUT}/${tag}.log`)) unlinkSync(`${OUT}/${tag}.log`);
  for (const f of readdirSync(OUT)) if (f.startsWith(`${tag}_day`) && f.endsWith('.bin')) unlinkSync(`${OUT}/${f}`);
  if (from) copyFileSync(from, start);
  const code = await run('bash', [GPULOCK, 'shared', 'node', 'scripts/spinup.mjs'], {
    N: String(n), TAG: tag, OUT, DAYS: String(end), MINUTES: '100000', KEEP: '2', OCEAN: '{"everySteps":8}', ...optionsOf(options, parameters), ...env,
  }, `${OUT}/${tag}.out`);
  if (start && existsSync(start) && start !== last) unlinkSync(start);
  return { log: `${OUT}/${tag}.log`, state: last, code };
}

export async function audit(state, options, file, parameters = PARAMETERS, env = {}) {
  if (!existsSync(file)) {
    await run('node', ['scripts/verticalAudit.mjs', state], { ...optionsOf(options, parameters), ...env }, `${file}.partial`);
    renameSync(`${file}.partial`, file);
  }
  return file;
}

// The mean of the audits of a state over windows started 0, 2 and 4 hours
// after it (SKIP steps of the audit), run side by side.
export async function auditWindows(state, options, stem, n, parameters = PARAMETERS2) {
  const steps = [0, 1, 2].map((k) => Math.round((k * 7200 * n) / (1350 * 16)));
  const files = await Promise.all(steps.map((skip) => audit(state, options, `${stem}.w${skip}.audit`, parameters, { SKIP: String(skip) })));
  return meanAudits(files.map(readAudit));
}

let startVolume = null;
export const arcticStart = async () => (startVolume ??= await arcticVolume(`${STATES}/nine64_day0091.bin`));
export async function screenValues(tagM, point, eightLog, eightState, arcticState) {
  const log = readLog(eightLog), day = log.days.find((d) => d.day === 186);
  const values = { nan: log.nan || !day ? 1 : 0, clamped: log.days.reduce((s, d) => s + d.clamped, 0) };
  if (existsSync(eightState)) Object.assign(values, readAudit(await audit(eightState, point, `${OUT}/${tagM}.audit`)));
  if (day) Object.assign(values, { balance: day.meanAsr - day.meanOlr, albedo: day.meanAlbedo, dayRain: day.meanPrecip, instantAlbedo: day.albedo, instantBalance: day.asr - day.olr, stress: log.stress });
  const iceStart = await arcticStart();
  if (existsSync(arcticState)) { const iceEnd = await arcticVolume(arcticState); Object.assign(values, { iceStart, iceEnd, arctic: (iceStart - iceEnd) / 3 }); }
  return values;
}
