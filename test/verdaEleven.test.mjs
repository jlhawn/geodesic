import { test, before, after } from 'node:test';
import assert from 'node:assert/strict';
import { writeFileSync, readFileSync, mkdtempSync, mkdirSync, rmSync, chmodSync, existsSync, readdirSync, utimesSync } from 'node:fs';
import { tmpdir, hostname } from 'node:os';
import { join } from 'node:path';
import { spawnSync, execFileSync } from 'node:child_process';
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

test('the suite report names every failing test with its error and output, and the file that failed to load', () => {
  const at = script('suite', {
    'bad.test.mjs': "import { test, describe } from 'node:test';\nimport assert from 'node:assert/strict';\ndescribe('outer', () => { test('fine', () => {}); test('two numbers', () => { console.log('printed by the test'); assert.equal(2.68e-4, 2.66e-4); }); });\n",
    'crash.test.mjs': 'throw new Error("cannot load");\n',
    'good.test.mjs': "import { test } from 'node:test';\ntest('fine', () => {});\n",
  });
  const report = join(at, 'report.txt');
  const result = run('scripts/suiteReport.sh', [report], { FILES: ['bad', 'crash', 'good'].map((f) => join(at, `${f}.test.mjs`)).join(' '), JOBS: '3' });
  assert.equal(result.status, 1);
  const text = readFileSync(report, 'utf8');
  assert.match(text, /3 files, 1 passed, 2 failed; tests: 2 passed, 2 failed/);
  assert.match(text, /not ok 2 - two numbers[\s\S]*actual: 0\.000268[\s\S]*printed by the test/);
  assert.match(text, /crash\.test\.mjs \(exit 1[\s\S]*Error: cannot load/);
  assert.match(text, /seconds per file, slowest first:[\s\S]*good\.test\.mjs\n/);
  assert.ok(existsSync(join(at, 'report', 'bad.tap')));
});

const PAIRED = `#!/bin/bash
env | grep -E '^(NS|PREFIX|LEVELS|PER_YEAR|KEEP|OCEAN|STRATOSPHERE|UNTIL|OUT|FROM)=' | sort > "$OUT/paired.env"
echo call >> "$OUT/paired.calls"
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
  assert.equal(readFileSync(join(out, 'paired.env'), 'utf8'), [`KEEP=1000`, `LEVELS=bl36`, `NS=64 128`, `OCEAN={"everySteps":8}`, `OUT=${out}`, `PER_YEAR=36`, `PREFIX=eleven`, `STRATOSPHERE=1`, `UNTIL=1095`, ''].join('\n'));
  assert.match(readFileSync(join(out, 'STARTED_eleven'), 'utf8'), /commit [0-9a-f]+/);
  assert.match(readFileSync(join(out, 'ENDED_eleven'), 'utf8'), /exit 0: .*stopped \(day 1095 reached\)/);
  assert.match(r.stdout, /disk: N=64 day 0, 108 segments left × 106 MB; N=128 day 0, 108 segments left × 425 MB; the states take 57\.4 GB at day 1095/);

  r = eleven(out, ['resume']);
  assert.match(r.stdout, /resume: the run ended/);
  assert.equal(calls(), 1);

  rmSync(join(out, 'ENDED_eleven'));
  writeFileSync(join(out, 'eleven128_day0010.bin'), Buffer.alloc(1000));
  r = eleven(out, ['resume']);
  assert.equal(calls(), 2);
  assert.match(r.stdout, /N=128 day 10, 107 segments left × 0 MB/);

  rmSync(join(out, 'ENDED_eleven'));
  writeFileSync(join(out, 'STOP_eleven'), '');
  assert.match(eleven(out, ['resume']).stdout, /resume: STOP_eleven is in/);
  r = eleven(out, ['run']);
  assert.equal(r.status, 1);
  assert.match(r.stdout, /remove it to run/);
  assert.equal(calls(), 2);
  rmSync(join(out, 'STOP_eleven'));

  r = eleven(out, ['run'], { MARGIN_GB: '1000000000' });
  assert.equal(r.status, 1);
  assert.match(r.stdout, /not started, the disk lacks room/);
  assert.equal(eleven(out, ['run'], { MARGIN_GB: '1000000000', FORCE_DISK: '1' }).status, 0);
  assert.equal(calls(), 3);

  r = eleven(out, ['status']);
  assert.match(r.stdout, /N=128: eleven128_day0010\.bin \(day 10 of 1095, 107 segments left\)/);
  assert.match(r.stdout, /STARTED_eleven: /);
  assert.match(readFileSync(join(out, 'eleven.stages.log'), 'utf8'), /\d{4}-\d\d-\d\d \d\d:\d\d:\d\d run: NS="64 128" PREFIX=eleven LEVELS=bl36 PER_YEAR=36 KEEP=1000/);
  assert.equal(eleven(out, ['nonsense']).status, 2);
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
  const env = { FAKE: fake, HOME: join(fake, 'home'), VERDA: join(fake, 'verda'), SSH: join(fake, 'ssh'), REMOTE_REPO: join(fake, 'remote', 'geodesic'), ALLOW_DIRTY: '1' };
  let r = run('scripts/verdaPush.sh', [], env);
  assert.equal(r.status, 1);
  assert.match(r.stderr, /gcm-eleven is absent/);
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

test('the puller survives an absent instance and a new address, checks the states\' sizes, and verifies every file by SHA-256', () => {
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
  writeFileSync(join(out, 'ENDED_eleven'), 'exit 0: stopped (day 1095 reached)\n');
  r = run('scripts/verdaPull.sh', ['--verify'], env);
  assert.equal(r.status, 0, r.stdout + r.stderr);
  const report = readFileSync(join(dest, 'verify.txt'), 'utf8');
  assert.match(report, /^verified \S+: 5 files on root@10\.0\.0\.2:.* \(0\.00 GB\), 5 identical in .*, 0 problems\nthe run had ended: exit 0: stopped \(day 1095 reached\)\n1 files that differed were pulled again by checksum/);
  assert.match(report, /\n1000\t[0-9a-f]{64}\televen64_day0010\.bin\n/);
  assert.deepEqual(readFileSync(join(dest, 'eleven64_day0010.bin')), readFileSync(join(out, 'eleven64_day0010.bin')));

  writeFileSync(join(fake, 'list.json'), '[]');
  r = run('scripts/verdaPull.sh', ['--verify'], env);
  assert.equal(r.status, 1);
  assert.match(r.stdout, /cannot verify: gcm-eleven is absent/);
});

test('the startup script does nothing without a checkout and resumes from the kept volume\'s checkout otherwise', () => {
  const text = readFileSync(join(root, 'scripts/verdaElevenStartup.sh'), 'utf8');
  assert.match(text, /\[ -x "\$REPO\/scripts\/verdaEleven\.sh" \] \|\| exit 0/);
  assert.match(text, /OUT=\/root\/runs\/eleven scripts\/verdaEleven\.sh resume --detach/);
  assert.equal(spawnSync('bash', ['-n', join(root, 'scripts/verdaElevenStartup.sh')]).status, 0);
  assert.ok(hostname());
});
