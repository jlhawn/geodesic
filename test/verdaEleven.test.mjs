import { test, before, after } from 'node:test';
import assert from 'node:assert/strict';
import { writeFileSync, readFileSync, mkdtempSync, mkdirSync, rmSync, chmodSync, existsSync, readdirSync, utimesSync } from 'node:fs';
import { tmpdir, hostname } from 'node:os';
import { join } from 'node:path';
import { spawn, spawnSync, execFileSync } from 'node:child_process';
import { segmentTargets, timings, projection, reportText } from '../scripts/verdaBenchmark.mjs';
import { wanted, listFiles } from '../scripts/verdaFiles.mjs';

const root = new URL('..', import.meta.url).pathname;
let dir;
before(() => { dir = mkdtempSync(join(tmpdir(), 'verdaEleven-')); });
after(() => rmSync(dir, { recursive: true, force: true }));

const script = (name, files) => {
  const at = join(dir, name);
  mkdirSync(at, { recursive: true });
  for (const [file, text] of Object.entries(files)) { writeFileSync(join(at, file), text); chmodSync(join(at, file), 0o755); }
  return at;
};
const run = (file, args, env, options = {}) => spawnSync('bash', [join(root, file), ...args], { cwd: root, encoding: 'utf8', env: { ...process.env, ...env }, ...options });

test('PER_YEAR 36 cuts three years into 108 segments of 10 or 11 days that end on days 365, 730 and 1095', () => {
  const t = segmentTargets(36, 1095);
  assert.equal(t.length, 108);
  for (const day of [10, 20, 30, 41, 365, 730, 1095]) assert.ok(t.includes(day), `${day}`);
  const lengths = t.map((d, i) => d - (t[i - 1] ?? 0));
  assert.deepEqual([Math.min(...lengths), Math.max(...lengths)], [10, 11]);
  assert.deepEqual(segmentTargets(4, 365), [91, 183, 274, 365]);
  assert.deepEqual(segmentTargets(4, 400), [91, 183, 274, 365, 456]);
});

test('the benchmark times setup, the first day, the steady day and the finish from its lines and projects hours and cost', () => {
  const lines = [
    { t: 0.5, line: 'land: record 0 days old' },
    { t: 30, line: '--- 2026-10-02T00:00:00.000Z fresh start at N=128 on bl36 (36 layers) (163842 cells, dt 168.75 s, 512 steps a day) after 29 s of setup' },
    { t: 40, line: 'day 1 (0.2 min): Ts 14 °C' },
    { t: 40.1, line: 'stratosphere day 1, layer-mean temperature' },
    { t: 48, line: 'day 2 (0.3 min): Ts 14 °C' },
    { t: 49, line: 'saved bench128_day0002.bin after 0.8 min' },
  ];
  const t = timings(lines, 52);
  assert.deepEqual(t, { days: 2, setup: 30, firstDay: 10, steady: 8, finish: 4 });
  assert.equal(timings([{ t: 1, line: 'no start' }], 2), null);
  assert.equal(timings([lines[1], lines[2]], 41).steady, 10, 'one day: the steady rate is the first day');

  const p = projection({ cases: [{ N: 128, t, bytes: 310e6 }, { N: 64, t: { days: 5, setup: 10, firstDay: 2, steady: 1, finish: 1 }, bytes: 100e6 }], perYear: 36, until: 1095, price: 2, compare: 9 });
  assert.equal(p.segments, 108);
  assert.equal(p.days, 1095);
  const h128 = (108 * (30 + 10 - 8 + 4) + 1095 * 8) / 3600, h64 = (108 * (10 + 2 - 1 + 1) + 1095 * 1) / 3600, hc = 108 * 9 / 3600;
  assert.ok(Math.abs(p.rows[0].hours - h128) < 1e-12);
  assert.ok(Math.abs(p.rows[1].hours - h64) < 1e-12);
  assert.ok(Math.abs(p.hours - (h128 + h64 + hc)) < 1e-12);
  assert.ok(Math.abs(p.cost - 2 * (h128 + h64 + hc)) < 1e-12);
  assert.ok(Math.abs(p.gigabytes - 108 * 0.41) < 1e-9);

  const text = reportText({ header: 'H', config: 'C', measured: [{ N: 128, t, bytes: 310e6 }], plan: { compareSeconds: 9, projection: p } });
  assert.match(text, /N=128: 2 days; setup 30\.0 s, first day 10\.00 s, steady 8\.00 s a model day, finish 4\.0 s, state 310 MB/);
  assert.match(text, /N=128: 108 × 36\.0 s \+ 1095 × 8\.00 s = 3\.51 h/);
  assert.match(text, /N=64: 108 × 12\.0 s \+ 1095 × 1\.00 s = 0\.66 h \(2\.18 s a model day with the segments' overhead\), 10\.8 GB of states/);
  assert.match(text, /comparisons: 108 × 9\.0 s = 0\.27 h\n  total 4\.45 h, at 2 \$\/h \$8\.89; the volume holds 44\.3 GB of states at the end/);
  assert.match(text, /^RESULT hours=4\.45 cost=8\.89 price=2 gigabytes=44\.3 n128_s_per_day=8\.00 n128_hours=3\.51 n64_s_per_day=1\.00 n64_hours=0\.66$/m);
});

test('the pulled files are the whole-day states, logs and reports, without partial files, in-day checkpoints, bench states or locks', () => {
  for (const [path, yes] of [
    ['eleven128_day0010.bin', true], ['eleven64_day1095.bin', true], ['eleven.log', true], ['eleven_compare.md', true], ['ENDED_eleven', true],
    ['eleven_suite/grid.tap', true], ['eleven_bench/bench64.log', true],
    ['eleven128_day0020.bin.partial', false], ['eleven64_day0012_step0040.bin', false], ['eleven_bench/bench64_day0002.bin', false],
    ['twin64_day0810.bin', false], ['eleven.lock/owner', false],
  ]) assert.equal(wanted(path, 'eleven'), yes, path);
  const out = join(dir, 'files');
  mkdirSync(join(out, 'eleven_suite'), { recursive: true });
  mkdirSync(join(out, 'eleven.lock'));
  for (const f of ['eleven64_day0010.bin', 'eleven64_day0011_step0008.bin', 'eleven.log', 'eleven_suite/a.tap', 'eleven.lock/owner']) writeFileSync(join(out, f), 'x');
  assert.deepEqual(listFiles(out, 'eleven'), ['eleven.log', 'eleven64_day0010.bin', 'eleven_suite/a.tap']);
});

test('the pulled files include the asynchronous spin-up\'s ocean-only year files, cycle ledgers and logs, not its forcing records or in-year checkpoints', () => {
  for (const [path, yes] of [
    ['twelve64_c01_year0020.bin', true], ['twelve128_c03_year0100.bin', true], ['twelve64_cycles.txt', true], ['twelve64_cycles.log', true], ['twelve64_c01.log', true], ['twelve64_c01.out', true],
    ['twelve64_c01_year0003_day120.bin', false], ['twelve64_c01_year0003_day120_step0040.bin', false], ['twelve64_c01_year0020.bin.partial', false],
    ['twelve64_c01_forcing/forcing-0731.bin', false], ['twelve64_c01_forcing/forcing-0732.bin.partial', false], ['twelve64_c01_forcing/notes.txt', false],
    ['twelve64_cycles.lock/pid', false], ['eleven64_c01_year0020.bin', false], ['twelve64_c01_year0020_copy.bin', false],
  ]) assert.equal(wanted(path, 'twelve'), yes, path);
  const out = join(dir, 'asyncFiles');
  mkdirSync(join(out, 'twelve64_c01_forcing'), { recursive: true });
  mkdirSync(join(out, 'twelve64_cycles.lock'));
  for (const f of ['twelve64_day1095.bin', 'twelve64_c01_year0007.bin', 'twelve64_c01_year0007_day030.bin', 'twelve64_cycles.txt', 'twelve64_cycles.log', 'twelve64_c01_forcing/forcing-0731.bin', 'twelve64_c01_forcing/index.txt', 'twelve64_cycles.lock/pid']) writeFileSync(join(out, f), 'x');
  assert.deepEqual(listFiles(out, 'twelve'), ['twelve64_c01_year0007.bin', 'twelve64_cycles.log', 'twelve64_cycles.txt', 'twelve64_day1095.bin']);
  const listed = spawnSync(process.execPath, [join(root, 'scripts/verdaFiles.mjs'), out, 'twelve'], { encoding: 'utf8' }).stdout;
  assert.equal(listed, ['twelve64_c01_year0007.bin', 'twelve64_cycles.log', 'twelve64_cycles.txt', 'twelve64_day1095.bin'].map((f) => `1\t-\t${f}\n`).join(''));
});

test('the suite report names every failing test with its error and output, and the file that failed to load', () => {
  const at = script('suite', {
    'bad.test.mjs': "import { test, describe } from 'node:test';\nimport assert from 'node:assert/strict';\ndescribe('outer', () => { test('fine', () => {}); test('two numbers', () => { console.log('printed by the test'); assert.equal(2.68e-4, 2.66e-4); }); });\n",
    'crash.test.mjs': 'throw new Error("cannot load");\n',
    'good.test.mjs': "import { test } from 'node:test';\ntest('fine', () => {});\n",
  });
  const report = join(at, 'report.txt');
  const result = run('scripts/suiteReport.sh', [report], { FILES: ['bad', 'crash', 'good'].map((f) => join(at, `${f}.test.mjs`)).join(' '), JOBS: '3', DIR: join(at, 'report') });
  assert.equal(result.status, 1);
  const text = readFileSync(report, 'utf8');
  assert.match(text, /3 files, 1 passed, 2 failed; tests: 2 passed, 2 failed/);
  assert.match(text, /not ok 2 - two numbers[\s\S]*actual: 0\.000268[\s\S]*printed by the test/);
  assert.match(text, /crash\.test\.mjs \(exit 1[\s\S]*Error: cannot load/);
  assert.match(text, /seconds per file, slowest first:[\s\S]*good\.test\.mjs\n/);
  assert.ok(existsSync(join(at, 'report', 'bad.tap')));
});

const PAIRED = `#!/bin/bash
env | grep -E '^(NS|PREFIX|LEVELS|PER_YEAR|KEEP|OCEAN|STRATOSPHERE|UNTIL|OUT|FROM|SYNC_CMD|RADIATION|MOIST)=' | sort > "$OUT/paired.env"
echo call >> "$OUT/paired.calls"
[ -f "$OUT/killed" ] && exit 143
echo "2026-10-02 12:00 stopped (day $UNTIL reached)" >> "$OUT/$PREFIX.log"
`;
const eleven = (out, args, env = {}) => run('scripts/verdaEleven.sh', args, { OUT: out, SKIP_GPU_CHECK: '1', PAIRED: join(dir, 'paired', 'paired.sh'), MARGIN_GB: '0', FROM: 'should-not-pass', ...env });

test('verdaEleven runs the paired spin-up with the run\'s settings, marks it, and resumes only a run that was started and has not ended or been stopped', () => {
  script('paired', { 'paired.sh': PAIRED });
  const out = join(dir, 'eleven');
  const calls = () => (existsSync(join(out, 'paired.calls')) ? readFileSync(join(out, 'paired.calls'), 'utf8').trim().split('\n').length : 0);

  let r = eleven(out, ['resume']);
  assert.equal(r.status, 0);
  assert.match(r.stdout, /resume: the run was never started/);
  assert.equal(calls(), 0);

  r = eleven(out, ['run']);
  assert.equal(r.status, 0, r.stdout + r.stderr);
  assert.equal(calls(), 1);
  assert.equal(readFileSync(join(out, 'paired.env'), 'utf8'), [`KEEP=1000`, `LEVELS=bl60`, `NS=64 128`, `OCEAN={"everySteps":8}`, `OUT=${out}`, `PER_YEAR=36`, `PREFIX=eleven`, `STRATOSPHERE=1`, `SYNC_CMD=sync "$1"`, `UNTIL=1095`, ''].join('\n'));
  assert.match(readFileSync(join(out, 'STARTED_eleven'), 'utf8'), /commit [0-9a-f]+/);
  assert.match(readFileSync(join(out, 'ENDED_eleven'), 'utf8'), /exit 0: .*stopped \(day 1095 reached\)/);
  assert.match(r.stdout, /disk: N=64 day 0, 108 segments left × 106 MB; N=128 day 0, 108 segments left × 425 MB; the states take 57\.4 GB at day 1095/);

  r = eleven(out, ['resume']);
  assert.match(r.stdout, /resume: the run ended/);
  assert.equal(calls(), 1);

  rmSync(join(out, 'ENDED_eleven'));
  writeFileSync(join(out, 'eleven128_day0010.bin'), Buffer.alloc(1000));
  writeFileSync(join(out, 'killed'), '');
  r = eleven(out, ['resume']);
  assert.equal(calls(), 2);
  assert.match(r.stdout, /N=128 day 10, 107 segments left × 0 MB/);
  assert.equal(r.status, 1);
  assert.match(r.stdout, /run: interrupted \(exit 143\) before .* logged a stop; not marked ended, so resume continues it/);
  assert.ok(!existsSync(join(out, 'ENDED_eleven')), 'an earlier stop line in the log does not mark a killed run ended');
  rmSync(join(out, 'killed'));

  r = eleven(out, ['resume']);
  assert.equal(calls(), 3);
  assert.match(readFileSync(join(out, 'ENDED_eleven'), 'utf8'), /exit 0: .*stopped \(day 1095 reached\)/);

  rmSync(join(out, 'ENDED_eleven'));
  writeFileSync(join(out, 'STOP_eleven'), '');
  assert.match(eleven(out, ['resume']).stdout, /resume: STOP_eleven is in/);
  r = eleven(out, ['run']);
  assert.equal(r.status, 1);
  assert.match(r.stdout, /remove it to run/);
  assert.equal(calls(), 3);
  rmSync(join(out, 'STOP_eleven'));

  r = eleven(out, ['run'], { MARGIN_GB: '1000000000' });
  assert.equal(r.status, 1);
  assert.match(r.stdout, /not started, the disk lacks room/);
  assert.equal(eleven(out, ['run'], { MARGIN_GB: '1000000000', FORCE_DISK: '1' }).status, 0);
  assert.equal(calls(), 4);

  r = eleven(out, ['status']);
  assert.match(r.stdout, /N=128: eleven128_day0010\.bin \(day 10 of 1095, 107 segments left\)/);
  assert.match(r.stdout, /STARTED_eleven: /);
  assert.match(readFileSync(join(out, 'eleven.stages.log'), 'utf8'), /\d{4}-\d\d-\d\d \d\d:\d\d:\d\d run: NS="64 128" PREFIX=eleven LEVELS=bl60 PER_YEAR=36 KEEP=1000/);
  assert.equal(eleven(out, ['nonsense']).status, 2);

  rmSync(join(out, 'ENDED_eleven'));
  r = eleven(out, ['run'], { UNTIL: '1825' });
  assert.equal(calls(), 5);
  assert.equal(readFileSync(join(out, 'UNTIL_eleven'), 'utf8').trim(), '1825');
  assert.match(readFileSync(join(out, 'paired.env'), 'utf8'), /UNTIL=1825/);
  assert.match(readFileSync(join(out, 'ENDED_eleven'), 'utf8'), /stopped \(day 1825 reached\)/);
  rmSync(join(out, 'ENDED_eleven'));
  r = eleven(out, ['resume']);
  assert.equal(calls(), 6);
  assert.match(readFileSync(join(out, 'paired.env'), 'utf8'), /UNTIL=1825/, 'a resume without UNTIL keeps the extended end');
  assert.match(r.stdout, /N=128 day 10, 179 segments left/);
  assert.match(eleven(out, ['status']).stdout, /day 10 of 1825, 179 segments left/);
});

test('verdaEleven keeps the run\'s RADIATION and MOIST in its settings, and the resume runs with them', () => {
  script('paired', { 'paired.sh': PAIRED });
  const out = join(dir, 'options'), RADIATION = '{"cloudColumns":true,"cloudSpread":0.75}', MOIST = '{"stratocumulusInversion":7}';
  mkdirSync(out, { recursive: true });
  writeFileSync(join(out, 'killed'), '');
  eleven(out, ['run'], { RADIATION, MOIST });
  const settings = readFileSync(join(out, 'SETTINGS_eleven'), 'utf8');
  assert.ok(settings.includes(`RADIATION=${RADIATION}\n`) && settings.includes(`MOIST=${MOIST}\n`), settings);
  rmSync(join(out, 'killed'));
  writeFileSync(join(out, 'eleven128_day0010.bin'), Buffer.alloc(1000));
  const r = eleven(out, ['resume']);
  assert.equal(r.status, 0, r.stdout + r.stderr);
  const env = readFileSync(join(out, 'paired.env'), 'utf8');
  assert.ok(env.includes(`RADIATION=${RADIATION}\n`) && env.includes(`MOIST=${MOIST}\n`), env);
});

const ASYNC = `#!/bin/bash
echo "N=$N $(env | grep -E '^(PREFIX|DRIVER|CYCLES|COUPLED_YEARS|RECORD_YEARS|OCEAN_YEARS|FINAL_YEARS|OCEAN_ONLY|RESTORE|OCEAN_KEEP|SNAPSHOT_DAYS|KEEP_FORCING|KEEP|PER_YEAR|LEVELS|OCEAN|STRATOSPHERE|SYNC_CMD|FROM|UNTIL)=' | sort | tr '\\n' ' ')" >> "$OUT/async.calls"
tag=$PREFIX$N
note() { echo "$(date '+%Y-%m-%d %H:%M') $*" >> "$OUT/\${tag}_cycles.log"; }
[ -f "$OUT/\${tag}_cycles.txt" ] || echo start=0 > "$OUT/\${tag}_cycles.txt"
note "asynchronous spin-up of $tag from day 0"
case "$(cat "$OUT/behave.$N" 2>/dev/null)" in
  signal) : > "$OUT/\${tag}_day0400.bin"; note "stopped by a signal at coupled day 400"; exit 0 ;;
  stopsignal) : > "$OUT/\${tag}_day0400.bin"; : > "$OUT/STOP_$PREFIX"; note "stopped by a signal at coupled day 400"; exit 0 ;;
  nan) note "NaN: Ts NaN on day 3"; exit 2 ;;
  stop) : > "$OUT/STOP_$PREFIX"; note "stopped at STOP_$PREFIX or STOP_$tag at coupled day 0"; exit 0 ;;
  locked) note "another driver (pid 4242) holds $OUT/\${tag}_cycles.lock; not starting"; exit 0 ;;
  fail) note "five failures in a row; giving up"; exit 1 ;;
  stopn) : > "$OUT/STOP_$tag"; note "stopped at STOP_$PREFIX or STOP_$tag at coupled day 0"; exit 0 ;;
  term) trap ': > "$OUT/\${tag}_day0400.bin"; note "stopped by a signal at coupled day 400"; exit 0' TERM; : > "$OUT/\${tag}.waiting"; sleep 60 & wait $!; exit 1 ;;
esac
: > "$OUT/\${tag}_day1825.bin"
echo done=1 >> "$OUT/\${tag}_cycles.txt"
note "done: 1 cycles, coupled day 1825"
`;
const SCHEDULE = { DRIVER: 'async', CYCLES: '1', COUPLED_YEARS: '3', RECORD_YEARS: '1', OCEAN_YEARS: '20', FINAL_YEARS: '2', OCEAN_ONLY: '{"everySteps":8}' };
const twelve = (out, args, env = {}) => run('scripts/verdaEleven.sh', args, { OUT: out, PREFIX: 'twelve', SKIP_GPU_CHECK: '1', PAIRED: '', ASYNC: join(dir, 'async', 'async.sh'), MARGIN_GB: '0', FROM: 'should-not-pass', ...env });
const calls = (out) => (existsSync(join(out, 'async.calls')) ? readFileSync(join(out, 'async.calls'), 'utf8').trim().split('\n') : []);

test('verdaEleven DRIVER=async runs the asynchronous schedule for each N in turn, keeps its settings for the resume, and the resume carries on with the N in progress', () => {
  script('async', { 'async.sh': ASYNC });
  const out = join(dir, 'twelve');
  mkdirSync(out, { recursive: true });
  writeFileSync(join(out, 'behave.128'), 'signal');

  let r = twelve(out, ['run'], SCHEDULE);
  assert.equal(r.status, 1, r.stdout + r.stderr);
  assert.equal(r.stderr, '');
  assert.match(r.stdout, /disk: N=64 from day 0 to 1825, 180 snapshots and 3 ocean-only year files × 106 MB and 0\.8 GB of forcing; N=128 from day 0 to 1825, 180 snapshots and 3 ocean-only year files × 425 MB and 3\.4 GB of forcing; the files take 101\.4 GB at day 1825, 101\.4 GB of it still to come/);
  assert.match(r.stdout, /run: NS="64 128" PREFIX=twelve .* UNTIL=1825 OUT=\S+ DRIVER=async CYCLES=1 COUPLED_YEARS=3 RECORD_YEARS=1 OCEAN_YEARS=20 FINAL_YEARS=2 OCEAN_ONLY=\{"everySteps":8\} RESTORE=30 OCEAN_KEEP=2 SNAPSHOT_DAYS=30 KEEP_FORCING=0 scripts\/asyncPair\.sh/);
  assert.match(r.stdout, /run: interrupted \(exit 143\) before scripts\/asyncPair\.sh logged a stop; not marked ended, so resume continues it/);
  assert.deepEqual(calls(out).map((c) => c.split(' ')[0]), ['N=64', 'N=128']);
  for (const word of ['PREFIX=twelve', 'CYCLES=1', 'COUPLED_YEARS=3', 'RECORD_YEARS=1', 'OCEAN_YEARS=20', 'FINAL_YEARS=2', 'OCEAN_ONLY={"everySteps":8}', 'OCEAN={"everySteps":8}', 'RESTORE=30', 'OCEAN_KEEP=2', 'SNAPSHOT_DAYS=30', 'KEEP_FORCING=0', 'KEEP=1000', 'PER_YEAR=36', 'LEVELS=bl60', 'STRATOSPHERE=1', 'SYNC_CMD=sync "$1"']) assert.ok(calls(out)[0].includes(` ${word} `), `${word} in ${calls(out)[0]}`);
  assert.doesNotMatch(calls(out)[0], /FROM=/);
  assert.ok(existsSync(join(out, 'STARTED_twelve')));
  assert.ok(!existsSync(join(out, 'ENDED_twelve')), 'an interrupted schedule is not marked ended');
  assert.equal(readFileSync(join(out, 'UNTIL_twelve'), 'utf8').trim(), '1825');
  const settings = readFileSync(join(out, 'SETTINGS_twelve'), 'utf8');
  for (const line of ['DRIVER=async', 'NS=64 128', 'PER_YEAR=36', 'KEEP=1000', 'LEVELS=bl60', 'OCEAN={"everySteps":8}', 'STRATOSPHERE=1', 'CYCLES=1', 'COUPLED_YEARS=3', 'OCEAN_YEARS=20', 'FINAL_YEARS=2', 'OCEAN_ONLY={"everySteps":8}']) assert.match(settings, new RegExp(`^${line.replace(/[{}"]/g, '\\$&')}$`, 'm'), line);
  const log = readFileSync(join(out, 'twelve.log'), 'utf8');
  assert.match(log, /N=64 ended: done at coupled day 1825\n/);
  assert.match(log, /N=128 interrupted by a signal at coupled day 400; a new start carries on from there\n/);
  assert.doesNotMatch(log, /stopped \(/);
  assert.match(twelve(out, ['status']).stdout, /N=128: twelve128_day0400\.bin; cycle 1 of 1: coupled, day 400 of 1095, recording the forcing from day 730; the schedule ends at coupled day 1825\n  .* stopped by a signal at coupled day 400\n/);

  rmSync(join(out, 'behave.128'));
  r = twelve(out, ['resume']);
  assert.equal(r.status, 0, r.stdout + r.stderr);
  assert.deepEqual(calls(out).map((c) => c.split(' ')[0]), ['N=64', 'N=128', 'N=128'], 'the resume runs the N in progress, not the first again');
  assert.ok(calls(out)[2].includes(' CYCLES=1 ') && calls(out)[2].includes(' OCEAN_YEARS=20 ') && calls(out)[2].includes(' FINAL_YEARS=2 '), 'the resume takes the schedule from SETTINGS_twelve');
  assert.match(r.stdout, /disk: N=64 from day 0 to 1825, 180 snapshots .* × 0 MB .*; N=128 from day 0 to 1825/);
  assert.match(readFileSync(join(out, 'twelve.log'), 'utf8'), /N=64 already ended \(done at coupled day 1825\); next\n/);
  assert.match(readFileSync(join(out, 'ENDED_twelve'), 'utf8'), /exit 0: .*stopped \(schedule ended: N=64 done at coupled day 1825; N=128 done at coupled day 1825\)/);
  assert.match(twelve(out, ['resume']).stdout, /resume: the run ended/);
  const status = twelve(out, ['status']).stdout;
  assert.match(status, /N=64: twelve64_day1825\.bin; done at coupled day 1825; the schedule ends at coupled day 1825\n  .* done: 1 cycles, coupled day 1825\n/);
  assert.match(status, /ENDED_twelve: .*stopped \(schedule ended/);

  rmSync(join(out, 'ENDED_twelve'));
  r = twelve(out, ['run']);
  assert.equal(r.status, 0, r.stdout + r.stderr);
  assert.equal(calls(out).length, 5, 'once the schedule has stopped, a new run gives every N its schedule again');
});

test('the asynchronous schedule names a NaN, a failure or a resolution\'s own STOP_<PREFIX><N> in its stop, stops at STOP_<PREFIX>, leaves no stop while another driver holds a resolution, and status shows each phase', () => {
  script('async', { 'async.sh': ASYNC });
  const fresh = (name, behave = {}) => {
    const out = join(dir, name);
    mkdirSync(out, { recursive: true });
    for (const [n, what] of Object.entries(behave)) writeFileSync(join(out, `behave.${n}`), what);
    return out;
  };
  let out = fresh('twelveNaN', { 64: 'nan' });
  let r = twelve(out, ['run'], SCHEDULE);
  assert.equal(r.status, 0, r.stdout + r.stderr);
  assert.match(readFileSync(join(out, 'ENDED_twelve'), 'utf8'), /exit 1: .*stopped \(schedule ended: N=64 NaN at coupled day 0; N=128 done at coupled day 1825\)/);
  assert.equal(calls(out).length, 2);

  out = fresh('twelveFail', { 64: 'fail' });
  r = twelve(out, ['run'], SCHEDULE);
  assert.equal(r.status, 0, r.stdout + r.stderr);
  assert.match(readFileSync(join(out, 'twelve.log'), 'utf8'), /N=64 ended: failed \(exit 1\) at coupled day 0, see twelve64_cycles\.log\n/);
  assert.match(readFileSync(join(out, 'ENDED_twelve'), 'utf8'), /exit 1: .*stopped \(schedule ended: N=64 failed \(exit 1\) at coupled day 0, see twelve64_cycles\.log; N=128 done at coupled day 1825\)/);
  assert.deepEqual(calls(out).map((c) => c.split(' ')[0]), ['N=64', 'N=128'], 'a failed resolution ends its own schedule only');

  out = fresh('twelveStopN', { 64: 'stopn' });
  r = twelve(out, ['run'], SCHEDULE);
  assert.equal(r.status, 0, r.stdout + r.stderr);
  assert.match(readFileSync(join(out, 'twelve.log'), 'utf8'), /N=64 ended: STOP_twelve64 at coupled day 0\n/);
  assert.match(readFileSync(join(out, 'ENDED_twelve'), 'utf8'), /exit 1: .*stopped \(schedule ended: N=64 STOP_twelve64 at coupled day 0; N=128 done at coupled day 1825\)/);
  assert.deepEqual(calls(out).map((c) => c.split(' ')[0]), ['N=64', 'N=128'], 'STOP_twelve64 ends N=64 and the next N runs');
  assert.ok(existsSync(join(out, 'STOP_twelve64')), 'the stop file stays until it is removed by hand');

  out = fresh('twelveStop', { 64: 'stop' });
  r = twelve(out, ['run'], SCHEDULE);
  assert.match(readFileSync(join(out, 'ENDED_twelve'), 'utf8'), /exit 0: .*stopped \(STOP_twelve: N=64 not started; N=128 not started\)/);
  assert.deepEqual(calls(out).map((c) => c.split(' ')[0]), ['N=64'], 'STOP_twelve stops the schedule before the next N');
  assert.match(twelve(out, ['resume']).stdout, /resume: the run ended/);

  out = fresh('twelveStopNow', { 128: 'stopsignal' });
  r = twelve(out, ['run'], SCHEDULE);
  assert.equal(r.status, 0, r.stdout + r.stderr);
  assert.match(readFileSync(join(out, 'ENDED_twelve'), 'utf8'), /exit 0: .*stopped \(STOP_twelve: N=64 done at coupled day 1825; N=128 at coupled day 400\)/, 'STOP_twelve with a SIGTERM ends the run at once');

  out = fresh('twelveLocked', { 128: 'locked' });
  r = twelve(out, ['run'], SCHEDULE);
  assert.equal(r.status, 1);
  assert.match(r.stdout, /run: interrupted \(exit 1\)/);
  assert.match(readFileSync(join(out, 'twelve.log'), 'utf8'), /N=128 not run: another driver \(pid 4242\) holds/);
  assert.ok(!existsSync(join(out, 'ENDED_twelve')));

  r = twelve(fresh('twelveFull'), ['run'], { ...SCHEDULE, MARGIN_GB: '1000000000' });
  assert.equal(r.status, 1);
  assert.match(r.stdout, /not started, the disk lacks room/);

  out = fresh('twelvePhases');
  writeFileSync(join(out, 'SETTINGS_twelve'), Object.entries({ ...SCHEDULE, NS: '64 128 192' }).map(([k, v]) => `${k}=${v}\n`).join(''));
  writeFileSync(join(out, 'twelve64_cycles.txt'), 'start=0\ndone=1\n');
  writeFileSync(join(out, 'twelve128_cycles.txt'), 'start=0\n');
  for (const f of ['twelve64_day1460.bin', 'twelve128_day1095.bin', 'twelve128_c01_year0007.bin', 'twelve128_c01_year0007_day120.bin']) writeFileSync(join(out, f), '');
  const status = twelve(out, ['status']).stdout;
  assert.match(status, /N=64: twelve64_day1460\.bin; after the last cycle: coupled, day 1460 of 1825; the schedule ends at coupled day 1825\n/);
  assert.match(status, /N=128: twelve128_day1095\.bin; cycle 1 of 1: ocean-only, 7 of 20 years done \(twelve128_c01_year0007_day120\.bin\); the schedule ends at coupled day 1825\n/);
  assert.match(status, /N=192: no state yet; not started; the schedule ends at coupled day 1825\n/);
  writeFileSync(join(out, 'twelve128_day0740.bin'), '');
  for (const f of ['twelve128_day1095.bin', 'twelve128_c01_year0007.bin', 'twelve128_c01_year0007_day120.bin']) rmSync(join(out, f));
  assert.match(twelve(out, ['status']).stdout, /N=128: twelve128_day0740\.bin; cycle 1 of 1: coupled, day 740 of 1095, recording the forcing since day 730;/);
});

test('an eviction\'s SIGTERM to the whole process group leaves the asynchronous run unmarked, and the resume runs only the interrupted N', async () => {
  script('async', { 'async.sh': ASYNC });
  const out = join(dir, 'twelveEvicted');
  mkdirSync(out, { recursive: true });
  writeFileSync(join(out, 'behave.128'), 'term');
  const until = async (check, what) => {
    const deadline = Date.now() + 20000;
    while (!check()) {
      if (Date.now() > deadline) throw new Error(`timed out waiting for ${what}`);
      await new Promise((resolve) => setTimeout(resolve, 50));
    }
  };
  const log = () => (existsSync(join(out, 'twelve.log')) ? readFileSync(join(out, 'twelve.log'), 'utf8') : '');
  const env = { ...process.env, OUT: out, PREFIX: 'twelve', SKIP_GPU_CHECK: '1', PAIRED: '', ASYNC: join(dir, 'async', 'async.sh'), MARGIN_GB: '0', ...SCHEDULE };
  const child = spawn('bash', [join(root, 'scripts/verdaEleven.sh'), 'run'], { cwd: root, env, detached: true, stdio: 'ignore' });
  const exited = new Promise((resolve) => child.on('exit', (code, signal) => resolve(code ?? signal)));
  await until(() => existsSync(join(out, 'twelve128.waiting')), 'N=128 to start');
  process.kill(-child.pid, 'SIGTERM');
  await exited;
  await until(() => /N=128 interrupted/.test(log()), 'asyncPair to log the interruption');

  assert.match(log(), /N=64 ended: done at coupled day 1825\n/);
  assert.match(log(), /N=128 interrupted by a signal at coupled day 400; a new start carries on from there\n/);
  assert.doesNotMatch(log(), /stopped \(/);
  assert.ok(!existsSync(join(out, 'ENDED_twelve')), 'an evicted run is not marked ended');
  assert.doesNotMatch(readFileSync(join(out, 'twelve.stages.log'), 'utf8'), /run: ended/);
  assert.match(readFileSync(join(out, 'twelve128_cycles.log'), 'utf8'), /stopped by a signal at coupled day 400\n/);
  assert.deepEqual(calls(out).map((c) => c.split(' ')[0]), ['N=64', 'N=128']);

  rmSync(join(out, 'behave.128'));
  const r = twelve(out, ['resume']);
  assert.equal(r.status, 0, r.stdout + r.stderr);
  assert.deepEqual(calls(out).map((c) => c.split(' ')[0]), ['N=64', 'N=128', 'N=128'], 'the resume runs the interrupted N only');
  assert.match(log(), /N=64 already ended \(done at coupled day 1825\); next\n/);
  assert.match(readFileSync(join(out, 'ENDED_twelve'), 'utf8'), /exit 0: .*stopped \(schedule ended: N=64 done at coupled day 1825; N=128 done at coupled day 1825\)/);
});

test('verdaEleven keeps one GPU mode at a time, takes over a lock whose process is gone, and detaches on request', () => {
  script('paired', { 'paired.sh': PAIRED });
  const out = join(dir, 'locked');
  mkdirSync(join(out, 'eleven.lock'), { recursive: true });
  const boot = execFileSync('bash', ['-c', "cat /proc/sys/kernel/random/boot_id 2>/dev/null || sysctl -n kern.boottime 2>/dev/null | tr -dc '0-9' | cut -c1-12"], { encoding: 'utf8' }).trim();
  writeFileSync(join(out, 'eleven.lock', 'owner'), `${process.pid} ${boot} bench\n`);
  let r = eleven(out, ['run']);
  assert.equal(r.status, 1);
  assert.match(r.stdout, /run not started: bench is running \(pid \d+\)/);

  writeFileSync(join(out, 'eleven.lock', 'owner'), `999999 ${boot} bench\n`);
  r = eleven(out, ['run']);
  assert.equal(r.status, 0, r.stdout);
  assert.match(r.stdout, /taking over the lock of bench \(pid 999999, gone\)/);
  assert.ok(!existsSync(join(out, 'eleven.lock')), 'the lock goes with the mode');

  rmSync(join(out, 'ENDED_eleven'));
  r = eleven(out, ['resume', '--detach']);
  assert.equal(r.status, 0);
  assert.match(r.stdout, /resume running in the background \(pid \d+\), output in .*eleven\.resume\.out/);
  const deadline = Date.now() + 20000;
  while (!existsSync(join(out, 'ENDED_eleven')) && Date.now() < deadline) execFileSync('sleep', ['0.2']);
  assert.ok(existsSync(join(out, 'ENDED_eleven')));
});

test('verdaEleven autostart installs and enables a unit that resumes the run at every boot', () => {
  const at = script('systemd', { systemctl: '#!/bin/bash\necho "$*" >> "$(dirname "$0")/calls"\n' });
  mkdirSync(join(at, 'units'), { recursive: true });
  const out = join(dir, 'autostart');
  const r = eleven(out, ['autostart'], { SYSTEMD_DIR: join(at, 'units'), SYSTEMCTL: join(at, 'systemctl'), HOME: '/root' });
  assert.equal(r.status, 0, r.stdout + r.stderr);
  assert.match(r.stdout, /autostart: gcm-eleven-resume\.service enabled/);
  const unit = readFileSync(join(at, 'units', 'gcm-eleven-resume.service'), 'utf8');
  assert.match(unit, new RegExp(`\\nEnvironment=HOME=/root OUT=${out} PREFIX=eleven\\n`));
  assert.match(unit, new RegExp(`\\nExecStart=/bin/bash ${join(root, 'scripts', 'verdaEleven.sh').replace(/\//g, '\\/')} resume\\n`));
  assert.match(unit, /\nType=simple\n[\s\S]*\nTimeoutStopSec=300\n[\s\S]*\nWantedBy=multi-user\.target\n$/);
  assert.equal(readFileSync(join(at, 'calls'), 'utf8'), 'daemon-reload\nenable gcm-eleven-resume.service\n');
});

const FAKE_SSH = `#!/bin/bash
while [ $# -gt 0 ]; do case "$1" in -o|-i|-p|-l|-F) shift 2 ;; -*) shift ;; *) break ;; esac; done
echo "$1" >> "$FAKE/ssh-hosts"; shift
exec bash -c "$*"
`;
const FAKE_VERDA = `#!/bin/bash
echo "$*" >> "$FAKE/verda-calls"
n=$(cat "$FAKE/round" 2>/dev/null || echo 0); n=$((n + 1)); echo $n > "$FAKE/round"
if [ -f "$FAKE/list.$n.json" ]; then cat "$FAKE/list.$n.json"; else cat "$FAKE/list.json"; fi
`;
const instance = (id, ip, status = 'running') => JSON.stringify([{ id, hostname: 'gcm-eleven', status, ...(ip ? { ip } : {}) }]);

test('the push stages HEAD without runs/, copies it over ssh to the address verda gives, and checks the commit and data there', () => {
  const fake = script('push', { ssh: FAKE_SSH, verda: FAKE_VERDA, 'list.1.json': '[]', 'list.json': instance('i-1', '10.0.0.1') });
  const env = { FAKE: fake, HOME: join(fake, 'home'), VERDA: join(fake, 'verda'), SSH: join(fake, 'ssh'), REMOTE_REPO: join(fake, 'remote', 'geodesic'), REMOTE_OUT: join(fake, 'remote', 'out'), ALLOW_DIRTY: '1' };
  let r = run('scripts/verdaPush.sh', [], env);
  assert.equal(r.status, 1);
  assert.match(r.stderr, /gcm-eleven is absent/);
  mkdirSync(env.REMOTE_OUT, { recursive: true });
  writeFileSync(join(env.REMOTE_OUT, 'STARTED_eleven'), '');
  r = run('scripts/verdaPush.sh', [], env);
  assert.equal(r.status, 1);
  assert.match(r.stderr, /was started and has not ended; set ALLOW_RUNNING=1/);
  assert.ok(!existsSync(env.REMOTE_REPO), 'nothing is pushed into a running run\'s checkout');
  writeFileSync(join(env.REMOTE_OUT, 'ENDED_eleven'), '');
  r = run('scripts/verdaPush.sh', [], env);
  assert.equal(r.status, 0, r.stderr);
  const head = execFileSync('git', ['rev-parse', 'HEAD'], { cwd: root, encoding: 'utf8' }).trim();
  assert.match(r.stdout, /staged [0-9a-f]{7} .*: \d+ files, \d+ MB \(\d+ MB of data\/, \d+ MB of \.git; runs\/ left out\)/);
  assert.match(r.stdout, new RegExp(`pushed to root@10\\.0\\.0\\.1:.*: HEAD ${head.slice(0, 7)} and \\d+ data files of the same size there`));
  assert.equal(execFileSync('git', ['rev-parse', 'HEAD'], { cwd: env.REMOTE_REPO, encoding: 'utf8' }).trim(), head);
  assert.ok(existsSync(join(env.REMOTE_REPO, 'data', 'woa_annual_1deg.bin')));
  assert.ok(existsSync(join(env.REMOTE_REPO, 'scripts', 'pairedSpinup.sh')));
  assert.ok(!existsSync(join(env.REMOTE_REPO, 'runs')));
  assert.ok(readFileSync(join(fake, 'ssh-hosts'), 'utf8').trim().split('\n').every((h) => h === 'root@10.0.0.1'));
});

test('the puller survives an absent instance and a new address, checks the states\' sizes, verifies every file by SHA-256, and with KEEP_DAYS takes the newest state of each resolution and the listed days', () => {
  const fake = script('pull', { ssh: FAKE_SSH, verda: FAKE_VERDA, 'list.1.json': '[]', 'list.2.json': instance('i-1', '10.0.0.1'), 'list.3.json': instance('i-1', null, 'discontinued'), 'list.json': instance('i-2', '10.0.0.2') });
  const out = join(fake, 'remote', 'out'), dest = join(fake, 'dest');
  mkdirSync(join(out, 'eleven_bench'), { recursive: true });
  const files = { 'eleven64_day0010.bin': 1000, 'eleven128_day0010.bin': 3000, 'eleven128_day0020.bin.partial': 50, 'eleven64_day0012_step0040.bin': 40, 'eleven_bench/bench64_day0002.bin': 30, 'eleven_bench/bench64.log': 3, 'eleven.log': 4 };
  for (const [f, n] of Object.entries(files)) writeFileSync(join(out, f), Buffer.alloc(n, f.length));
  const env = { FAKE: fake, HOME: join(fake, 'home'), VERDA: join(fake, 'verda'), SSH: join(fake, 'ssh'), REMOTE_OUT: out, REMOTE_REPO: root, DEST: dest, INTERVAL: '0' };

  let r = run('scripts/verdaPull.sh', [], { ...env, ROUNDS: '4' });
  assert.equal(r.status, 0, r.stderr);
  assert.match(r.stdout, /gcm-eleven is absent; nothing pulled/);
  assert.match(r.stdout, /instance at root@10\.0\.0\.1\n.*2 states there \(0\.0 GB\): 2 here of the same size, 0 missing, 0 of another size; 2 other files; newest: /);
  assert.match(r.stdout, /gcm-eleven is discontinued \(i-1\); nothing pulled/);
  assert.match(r.stdout, /instance at root@10\.0\.0\.2/);
  assert.deepEqual(readdirSync(dest).sort(), ['eleven.log', 'eleven128_day0010.bin', 'eleven64_day0010.bin', 'eleven_bench', 'pull.log']);
  assert.deepEqual(readdirSync(join(dest, 'eleven_bench')), ['bench64.log']);
  assert.ok(existsSync(join(out, 'eleven128_day0020.bin.partial')), 'nothing is deleted there');

  writeFileSync(join(dest, 'eleven64_day0010.bin'), Buffer.alloc(1000, 7));
  const old = new Date('2020-01-01');
  utimesSync(join(dest, 'eleven64_day0010.bin'), old, old);
  utimesSync(join(out, 'eleven64_day0010.bin'), old, old);
  writeFileSync(join(out, 'eleven128_day0010.bin'), Buffer.alloc(2000, 9));
  r = run('scripts/verdaPull.sh', [], { ...env, ROUNDS: '1' });
  assert.match(r.stdout, /2 states there \(0\.0 GB\): 1 here of the same size, 0 missing, 1 of another size; 2 other files; eleven128_day0010\.bin \(3000 of 2000 bytes\)/);
  assert.equal(readFileSync(join(dest, 'eleven128_day0010.bin')).length, 3000, 'a round never replaces a state already here');

  writeFileSync(join(out, 'ENDED_eleven'), 'exit 0: stopped (day 1095 reached)\n');
  r = run('scripts/verdaPull.sh', ['--verify'], env);
  assert.equal(r.status, 0, r.stdout + r.stderr);
  const report = readFileSync(join(dest, 'verify.txt'), 'utf8');
  assert.match(report, /^verified \S+: 5 files on root@10\.0\.0\.2:.* \(0\.00 GB\), 5 identical in .*, 0 problems\nthe run had ended: exit 0: stopped \(day 1095 reached\)\n2 files that differed were pulled again by checksum/);
  assert.match(report, /\n1000\t[0-9a-f]{64}\televen64_day0010\.bin\n/);
  assert.deepEqual(readFileSync(join(dest, 'eleven64_day0010.bin')), readFileSync(join(out, 'eleven64_day0010.bin')));
  const kept = readdirSync(dest).filter((f) => /\.bin\.replaced-\d{8}T\d{6}Z$/.test(f)).sort();
  assert.equal(kept.length, 2, readdirSync(dest).join(' '));
  assert.deepEqual(readFileSync(join(dest, kept[1])), Buffer.alloc(1000, 7), 'verify keeps the copy it replaces');
  assert.equal(readFileSync(join(dest, kept[0])).length, 3000);

  const latest = join(fake, 'latest');
  for (const [f, n] of Object.entries({ 'eleven64_day0020.bin': 1000, 'eleven128_day0030.bin': 3000, 'eleven128_day0365.bin': 3000, 'eleven64_day0365.bin': 1000, 'eleven128_day0730.bin': 3000 })) writeFileSync(join(out, f), Buffer.alloc(n, f.length));
  r = run('scripts/verdaPull.sh', [], { ...env, DEST: latest, ROUNDS: '1', KEEP_DAYS: '365,1095' });
  assert.equal(r.status, 0, r.stderr);
  assert.match(r.stdout, /3 states wanted \(0\.0 GB\): 3 here of the same size, 0 missing, 0 of another size; 3 other files; 4 states left there; newest: eleven128_day0730\.bin eleven64_day0365\.bin/);
  assert.deepEqual(readdirSync(latest).filter((f) => f.endsWith('.bin')).sort(), ['eleven128_day0365.bin', 'eleven128_day0730.bin', 'eleven64_day0365.bin']);
  r = run('scripts/verdaPull.sh', ['--verify'], { ...env, DEST: latest, KEEP_DAYS: '365,1095' });
  assert.equal(r.status, 0, r.stdout + r.stderr);
  assert.match(readFileSync(join(latest, 'verify.txt'), 'utf8'), /^verified \S+: 6 files on root@10\.0\.0\.2:.*, 6 identical in .*, 0 problems\n/);

  const yearly = join(fake, 'yearly');
  mkdirSync(join(out, 'eleven64_c01_forcing'), { recursive: true });
  for (const [f, n] of Object.entries({ 'eleven64_c01_year0020.bin': 1000, 'eleven64_c01_year0003_day120.bin': 1000, 'eleven64_c01_forcing/forcing-0001.bin': 10, 'eleven64_cycles.txt': 8 })) writeFileSync(join(out, f), Buffer.alloc(n, 1));
  r = run('scripts/verdaPull.sh', [], { ...env, DEST: yearly, ROUNDS: '1', KEEP_DAYS: '365' });
  assert.equal(r.status, 0, r.stderr);
  assert.ok(existsSync(join(yearly, 'eleven64_c01_year0020.bin')), 'a KEEP_DAYS list still pulls the ocean-only year files');
  assert.ok(existsSync(join(yearly, 'eleven64_cycles.txt')));
  assert.ok(!existsSync(join(yearly, 'eleven64_c01_year0003_day120.bin')) && !existsSync(join(yearly, 'eleven64_c01_forcing')), 'in-year checkpoints and forcing records stay there');

  writeFileSync(join(fake, 'list.json'), '[]');
  r = run('scripts/verdaPull.sh', ['--verify'], env);
  assert.equal(r.status, 1);
  assert.match(r.stdout, /cannot verify: gcm-eleven is absent/);
});

test('the startup script does nothing without a checkout and resumes from the kept volume\'s checkout otherwise, for the run its PREFIX line names', () => {
  const text = readFileSync(join(root, 'scripts/verdaElevenStartup.sh'), 'utf8');
  assert.match(text, /^PREFIX=eleven$/m);
  assert.match(text, /\[ -x "\$REPO\/scripts\/verdaEleven\.sh" \] \|\| exit 0/);
  assert.equal(spawnSync('bash', ['-n', join(root, 'scripts/verdaElevenStartup.sh')]).status, 0);
  const at = script('startup', { 'gcm-twelve.sh': text.replace(/^PREFIX=eleven$/m, 'PREFIX=twelve') });
  const repo = join(at, 'repo'), runs = join(at, 'runs'), env = { ...process.env, REPO: repo, RUNS: runs, HOME: at };
  let r = spawnSync('bash', [join(at, 'gcm-twelve.sh')], { env, encoding: 'utf8' });
  assert.equal(r.status, 0);
  assert.ok(!existsSync(runs), 'nothing happens without a checkout');
  mkdirSync(join(repo, 'scripts'), { recursive: true });
  writeFileSync(join(repo, 'scripts', 'verdaEleven.sh'), '#!/bin/bash\necho "OUT=$OUT PREFIX=$PREFIX $* in $PWD"\n');
  chmodSync(join(repo, 'scripts', 'verdaEleven.sh'), 0o755);
  r = spawnSync('bash', [join(at, 'gcm-twelve.sh')], { env, encoding: 'utf8' });
  assert.equal(r.status, 0, r.stderr);
  const boot = readFileSync(join(runs, 'twelve', 'twelve.boot.log'), 'utf8');
  assert.match(boot, new RegExp(`startup script on ${hostname()}\\n`));
  assert.ok(boot.includes(`OUT=${runs}/twelve PREFIX=twelve resume --detach in ${repo}\n`), boot);
});
