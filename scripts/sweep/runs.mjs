// The sweep's parameters and the runs it makes: short GPU spin-ups from
// copies of saved states under the shared GPU lock (GPULOCK), each copy
// removed once the run has saved its last day, and the CPU audit of a state.
import { spawn } from 'node:child_process';
import { copyFileSync, existsSync, unlinkSync, openSync, closeSync } from 'node:fs';
import { readLog, readAudit, arcticVolume } from './score.mjs';

export const ROOT = new URL('../../', import.meta.url).pathname;
export const OUT = process.env.OUT ?? `${ROOT}runs`;
export const SWEEP = `${OUT}/sweep`;
export const STATES = process.env.STATES ?? '/Users/jlhawn/git_repos/jlhawn/geodesic/runs';
export const GPULOCK = process.env.GPULOCK ?? '/private/tmp/claude-501/-Users-jlhawn-git-repos-jlhawn-geodesic/4ab28226-4f31-4147-a831-f542e02bf0fc/scratchpad/gpulock.sh';

export const PARAMETERS = [
  { key: 'varianceScale', module: 'radiation', option: 'varianceScale', base: 5, low: 2, high: 10 },
  { key: 'mixingLength', module: 'radiation', option: 'mixingLength', base: 300, low: 150, high: 600 },
  { key: 'stratiformHours', module: 'moist', option: 'stratiformLifetime', base: 3, low: 1, high: 6, scale: 3600 },
  { key: 'cloudHours', module: 'moist', option: 'cloudLifetime', base: 1, low: 0.5, high: 2, scale: 3600 },
  { key: 'plumeEntrainment', module: 'moist', option: 'plumeEntrainment', base: 0.1, low: 0.05, high: 0.2 },
  { key: 'plumeCape', module: 'moist', option: 'plumeCape', base: 120, low: 40, high: 200 },
  { key: 'minimumInversion', module: 'radiation', option: 'minimumInversion', base: 4, low: 2, high: 6 },
  { key: 'criticalHumidity', module: 'radiation', option: 'criticalHumidity', base: 0.8, low: 0.7, high: 0.9 },
  { key: 'seaDrag', module: 'surface', option: 'dragCoefficient', base: 1.2e-3, low: 1.0e-3, high: 1.5e-3 },
];

export function optionsOf(point) {
  const o = { radiation: {}, moist: {}, surface: {} };
  for (const p of PARAMETERS) o[p.module][p.option] = Number((point[p.key] * (p.scale ?? 1)).toPrecision(6));
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

export async function spinup({ tag, n, from = null, fromDay = 0, days, options, env = {} }) {
  const end = fromDay + days, last = dayName(tag, end);
  if (existsSync(last)) return { log: `${OUT}/${tag}.log`, state: last, code: 0 };
  const start = from ? dayName(tag, fromDay) : null;
  if (existsSync(`${OUT}/${tag}.log`)) unlinkSync(`${OUT}/${tag}.log`);
  if (from) copyFileSync(from, start);
  const code = await run('bash', [GPULOCK, 'shared', 'node', 'scripts/spinup.mjs'], {
    N: String(n), TAG: tag, OUT, DAYS: String(end), MINUTES: '100000', KEEP: '2', BATCH: '1', DAY_MEAN: '8', OCEAN: '{"everySteps":8}', ...optionsOf(options), ...env,
  }, `${OUT}/${tag}.out`);
  if (start && existsSync(start) && start !== last) unlinkSync(start);
  return { log: `${OUT}/${tag}.log`, state: last, code };
}

export async function audit(state, options, file) {
  if (!existsSync(file)) await run('node', ['scripts/verticalAudit.mjs', state], optionsOf(options), file);
  return file;
}

let startVolume = null;
export async function screenValues(tagM, point, eightLog, eightState, arcticState) {
  const log = readLog(eightLog), day = log.days.find((d) => d.day === 186);
  const values = { nan: log.nan || !day ? 1 : 0, clamped: log.days.reduce((s, d) => s + d.clamped, 0) };
  if (day) Object.assign(values, { balance: day.meanAsr - day.meanOlr, albedo: day.meanAlbedo, instantAlbedo: day.albedo, instantBalance: day.asr - day.olr, stress: log.stress });
  if (existsSync(eightState)) Object.assign(values, readAudit(await audit(eightState, point, `${OUT}/${tagM}.audit`)));
  const iceStart = startVolume ??= await arcticVolume(`${STATES}/nine64_day0091.bin`);
  if (existsSync(arcticState)) { const iceEnd = await arcticVolume(arcticState); Object.assign(values, { iceStart, iceEnd, arctic: (iceStart - iceEnd) / 3 }); }
  return values;
}
