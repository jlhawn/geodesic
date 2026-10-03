// What a paired spin-up will cost on this machine, measured: for each
// N:DAYS of CASES ('128:2 64:5') a fresh atlas start of scripts/spinup.mjs
// on LEVELS (bl36) with the run's OCEAN ('{"everySteps":8}') and
// STRATOSPHERE (1), as TAG bench<N> in OUT (its earlier bench<N> states and
// log are removed first, so every case starts fresh, and MINUTES (20)
// bounds it), timed from the arrival of the segment's own log lines: setup
// up to the '---' line, the first day up to the first 'day' line, the
// steady rate between the first and the last 'day' lines, and the finish
// (end-of-segment lines, saving and exit) after the last. One
// scripts/compareStates.mjs over the saved states times a round's
// comparison. The projection for the paired run (PER_YEAR snapshot days,
// pairedSpinup's segments, through UNTIL, 365·YEARS by default) is, per N,
// segments × (setup + first day − steady + finish) + days × steady, plus a
// comparison per round, and the cost that many hours at PRICE $/h (1.85).
// The report goes to stdout and REPORT, and its numbers as JSON to
// REPORT with .json in place of .txt; the bench states are removed at
// the end unless KEEP_STATES=1.
//   OUT=runs/bench REPORT=runs/benchmark.txt node scripts/verdaBenchmark.mjs
import { spawn, execFileSync } from 'node:child_process';
import { readdirSync, unlinkSync, statSync, writeFileSync, mkdirSync, existsSync } from 'node:fs';
import { hostname } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

export function segmentTargets(perYear, until) {
  const targets = [];
  for (let k = 1; ; k++) {
    const t = Math.floor(k * 365 / perYear + 0.5);
    targets.push(t);
    if (t >= until) return targets;
  }
}

export function timings(lines, exitAt) {
  const start = lines.find((l) => /^--- \S+ (fresh start|continuing from)/.test(l.line));
  const days = lines.filter((l) => /^day \d+ \(/.test(l.line));
  if (!start || !days.length) return null;
  const first = days[0].t, last = days[days.length - 1].t;
  const steady = days.length > 1 ? (last - first) / (days.length - 1) : first - start.t;
  return { days: days.length, setup: start.t, firstDay: first - start.t, steady, finish: exitAt - last };
}

export function projection({ cases, perYear, until, price, compare = 0 }) {
  const targets = segmentTargets(perYear, until), segments = targets.length, days = targets[segments - 1];
  const rows = cases.map(({ N, t, bytes }) => {
    const overhead = t.setup + t.firstDay - t.steady + t.finish, seconds = segments * overhead + days * t.steady;
    return { N, steady: t.steady, overhead, perDay: seconds / days, hours: seconds / 3600, gigabytes: (segments * (bytes ?? 0)) / 1e9 };
  });
  const compareHours = (segments * compare) / 3600, hours = rows.reduce((s, r) => s + r.hours, 0) + compareHours;
  return { segments, days, rows, compareHours, hours, cost: hours * price, price, gigabytes: rows.reduce((s, r) => s + r.gigabytes, 0) };
}

export function reportText({ header, config, measured, plan }) {
  const f = (x, d = 2) => x.toFixed(d), out = [header, config, ''];
  for (const { N, t, bytes } of measured) out.push(`N=${N}: ${t.days} days; setup ${f(t.setup, 1)} s, first day ${f(t.firstDay)} s, steady ${f(t.steady)} s a model day, finish ${f(t.finish, 1)} s, state ${(bytes / 1e6).toFixed(0)} MB`);
  if (plan.compareSeconds !== undefined) out.push(`comparison of a round's states: ${f(plan.compareSeconds, 1)} s`);
  const p = plan.projection;
  out.push('', `projection for the paired run: ${p.segments} segments a resolution to day ${p.days}, one N after the other`);
  for (const r of p.rows) out.push(`  N=${r.N}: ${p.segments} × ${f(r.overhead, 1)} s + ${p.days} × ${f(r.steady)} s = ${f(r.hours)} h (${f(r.perDay)} s a model day with the segments' overhead), ${f(r.gigabytes, 1)} GB of states`);
  out.push(`  comparisons: ${p.segments} × ${f(plan.compareSeconds ?? 0, 1)} s = ${f(p.compareHours)} h`);
  out.push(`  total ${f(p.hours)} h, at ${p.price} $/h $${f(p.cost)}; the volume holds ${f(p.gigabytes, 1)} GB of states at the end`);
  out.push('', `seconds per model day: ${p.rows.map((r) => `N=${r.N} ${f(r.steady)} steady, ${f(r.perDay)} with overhead`).join('; ')}`);
  out.push(`RESULT hours=${f(p.hours)} cost=${f(p.cost)} price=${p.price} gigabytes=${f(p.gigabytes, 1)} ${p.rows.map((r) => `n${r.N}_s_per_day=${f(r.steady)} n${r.N}_hours=${f(r.hours)}`).join(' ')}`);
  return out.join('\n') + '\n';
}

function runCase(root, env, log) {
  return new Promise((resolve) => {
    const t0 = performance.now(), lines = [];
    const child = spawn(process.execPath, [join(root, 'scripts/spinup.mjs')], { cwd: root, env, stdio: ['ignore', 'pipe', 'pipe'] });
    let pending = '', errors = '';
    child.stdout.on('data', (chunk) => {
      const t = (performance.now() - t0) / 1000, parts = (pending + chunk).split('\n');
      pending = parts.pop();
      for (const line of parts) { lines.push({ t, line }); log(line); }
    });
    child.stderr.on('data', (chunk) => { errors += chunk; });
    child.on('close', (code) => resolve({ code, lines, errors, exitAt: (performance.now() - t0) / 1000 }));
  });
}

if (process.argv[1] && import.meta.url === pathToFileURL(process.argv[1]).href) {
  const root = new URL('..', import.meta.url).pathname;
  const OUT = process.env.OUT ?? join(root, 'runs/bench'), REPORT = process.env.REPORT ?? join(OUT, 'benchmark.txt');
  const CASES = (process.env.CASES ?? '128:2 64:5').trim().split(/\s+/).map((c) => c.split(':').map(Number));
  const LEVELS = process.env.LEVELS ?? 'bl36', OCEAN = process.env.OCEAN ?? '{"everySteps":8}', STRATOSPHERE = process.env.STRATOSPHERE ?? '1';
  const PRICE = Number(process.env.PRICE ?? 1.85), PER_YEAR = Number(process.env.PER_YEAR ?? 36), YEARS = Number(process.env.YEARS ?? 3);
  const UNTIL = Number(process.env.UNTIL ?? 365 * YEARS), MINUTES = process.env.MINUTES ?? '20';
  if (CASES.some(([n, d]) => !(n > 0 && d >= 1))) throw new Error(`CASES is N:DAYS pairs such as '128:2 64:5', not ${process.env.CASES}`);
  mkdirSync(OUT, { recursive: true });
  const commit = (() => { try { return execFileSync('git', ['rev-parse', '--short', 'HEAD'], { cwd: root, encoding: 'utf8' }).trim(); } catch { return 'unknown'; } })();
  const states = (tag) => readdirSync(OUT).filter((f) => f.startsWith(`${tag}_day`) && f.endsWith('.bin'));
  const measured = [];
  let adapter = '';
  for (const [N, DAYS] of CASES) {
    const TAG = `bench${N}`;
    for (const f of states(TAG)) unlinkSync(join(OUT, f));
    if (existsSync(join(OUT, `${TAG}.log`))) unlinkSync(join(OUT, `${TAG}.log`));
    const env = { ...process.env, N: String(N), TAG, DAYS: String(DAYS), MINUTES, KEEP: '1', OUT, LEVELS, OCEAN, STRATOSPHERE };
    for (const name of ['FROM', 'OCEAN_FROM', 'LAND_FROM', 'ICE_FROM', 'RECORD', 'SYNC_CMD', 'STOP_AFTER_STEPS', 'TOP_BUDGET']) delete env[name];
    console.log(`== N=${N}, ${DAYS} days from a fresh start on ${LEVELS}`);
    const run = await runCase(root, env, (line) => console.log(line));
    const t = timings(run.lines, run.exitAt), saved = states(TAG);
    if (run.code !== 0 || !t || t.days < DAYS || !saved.length) {
      const why = `N=${N} failed (exit ${run.code}, ${t ? t.days : 0} of ${DAYS} days): ${(run.errors || run.lines.slice(-5).map((l) => l.line).join('\n')).trim().split('\n').slice(-8).join(' | ')}`;
      writeFileSync(REPORT, `WebGCM benchmark at ${commit} on ${hostname()}, ${new Date().toISOString()}\n${why}\n`);
      console.error(why);
      process.exit(1);
    }
    measured.push({ N, t, bytes: statSync(join(OUT, saved[0])).size, file: join(OUT, saved[0]) });
  }
  try {
    adapter = execFileSync(process.execPath, ['--input-type=module', '-e', "const { create } = await import('webgpu'); const a = await create([]).requestAdapter(); const i = a?.info ?? {}; console.log([i.vendor, i.architecture, i.device, i.description].filter(Boolean).join(' '));"], { cwd: root, encoding: 'utf8' }).trim();
  } catch { adapter = 'adapter not read'; }
  const c0 = performance.now();
  execFileSync(process.execPath, [join(root, 'scripts/compareStates.mjs'), ...measured.map((m) => m.file), '--window', '1'], { cwd: root, stdio: ['ignore', 'ignore', 'inherit'] });
  const compareSeconds = (performance.now() - c0) / 1000;
  const p = projection({ cases: measured, perYear: PER_YEAR, until: UNTIL, price: PRICE, compare: compareSeconds });
  const text = reportText({
    header: `WebGCM benchmark at ${commit} on ${hostname()} (${adapter}), ${new Date().toISOString()}`,
    config: `fresh atlas starts on ${LEVELS}, OCEAN=${OCEAN} STRATOSPHERE=${STRATOSPHERE}; the run: PER_YEAR=${PER_YEAR} to day ${UNTIL}; PRICE=${PRICE} $/h`,
    measured, plan: { compareSeconds, projection: p },
  });
  writeFileSync(REPORT, text);
  writeFileSync(REPORT.replace(/\.txt$/, '') + '.json', JSON.stringify({ commit, adapter, measured: measured.map(({ N, t, bytes }) => ({ N, ...t, bytes })), compareSeconds, ...p }, null, 1) + '\n');
  if (process.env.KEEP_STATES !== '1') for (const m of measured) unlinkSync(m.file);
  process.stdout.write('\n' + text);
}
