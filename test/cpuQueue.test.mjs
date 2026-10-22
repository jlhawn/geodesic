import { test } from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { spawn, spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { testConcurrency } from '../scripts/cpuQueue.mjs';

const script = fileURLToPath(new URL('../scripts/cpuQueue.mjs', import.meta.url));
const skip = process.platform !== 'darwin' && 'the CPU slots are macOS-only';

const JOB = `import fs from 'node:fs';
const [name, gate] = process.argv.slice(2);
fs.appendFileSync(process.env.TRACE, \`start \${name}\\n\`);
while (gate && !fs.existsSync(gate)) await new Promise((resolve) => setTimeout(resolve, 10));
fs.appendFileSync(process.env.TRACE, \`end \${name}\\n\`);
`;

const pause = (ms) => new Promise((resolve) => setTimeout(resolve, ms));

function queue(slots) {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'cpu-queue-'));
  const job = path.join(dir, 'job.mjs'), trace = path.join(dir, 'trace');
  fs.writeFileSync(job, JOB);
  fs.writeFileSync(trace, '');
  const env = { ...process.env, CPU_QUEUE_DIR: path.join(dir, 'queue'), CPU_QUEUE_SLOTS: String(slots), CPU_QUEUE_TICK_MS: '20', TRACE: trace };
  for (const name of ['CPU_QUEUE_HELD', 'CPU_SLOTS', 'CPU_JOB', 'GPU_WHO', 'NODE_TEST_CONTEXT']) delete env[name];
  const gate = (name) => path.join(dir, `${name}.go`);
  const start = (name, slots, { held = false, labels = {} } = {}) => {
    const argv = [script, 'run', path.join(dir, `${name}.log`), process.execPath, job, name, ...(held ? [gate(name)] : [])];
    const child = spawn(process.execPath, argv, { env: { ...env, ...labels, CPU_SLOTS: String(slots) }, stdio: 'ignore' });
    return new Promise((resolve) => child.on('exit', resolve));
  };
  const lines = () => fs.readFileSync(trace, 'utf8').split('\n').filter(Boolean);
  const until = async (done) => { while (!done()) await pause(20); };
  return { dir, env, job, start, lines, until, release: (name) => fs.writeFileSync(gate(name), '') };
}

test('a run waits until its slots are free, and a smaller run queued after it does not overtake it', { skip, timeout: 20000 }, async (t) => {
  const q = queue(2);
  t.after(() => q.release('a'));
  const a = q.start('a', 1, { held: true });
  await q.until(() => q.lines().includes('start a'));
  const b = q.start('b', 2);
  await pause(300);
  const c = q.start('c', 1);
  await pause(300);
  assert.deepEqual(q.lines(), ['start a']);
  q.release('a');
  assert.deepEqual(await Promise.all([a, b, c]), [0, 0, 0]);
  assert.deepEqual(q.lines(), ['start a', 'end a', 'start b', 'end b', 'start c', 'end c']);
  assert.match(fs.readFileSync(path.join(q.dir, 'b.log'), 'utf8'), /^start \d\d:\d\d:\d\d, 2 of 2 CPU slots after \d+ s waiting\n.*exit 0 /s);
});

test('a node --test run must pass --test-concurrency no larger than its slots', { skip, timeout: 20000 }, () => {
  assert.equal(testConcurrency(['node', 'x.mjs']), null);
  assert.equal(testConcurrency(['node', '--test', 'a.test.mjs']), 0);
  assert.equal(testConcurrency(['node', '--test', '--test-concurrency=3', 'a.test.mjs']), 3);
  assert.equal(testConcurrency(['node', '--test-concurrency', '2', '--test', 'a.test.mjs']), 2);
  const q = queue(4), file = path.join(q.dir, 'ok.test.mjs');
  fs.writeFileSync(file, `import { test } from 'node:test';\ntest('ok', () => {});\n`);
  const run = (name, ...flags) => spawnSync(process.execPath, [script, 'run', path.join(q.dir, `${name}.log`), process.execPath, '--test', ...flags, file], {
    env: { ...q.env, CPU_SLOTS: '2' }, encoding: 'utf8',
  });
  const unbounded = run('unbounded'), wide = run('wide', '--test-concurrency=3'), fits = run('fits', '--test-concurrency=2');
  assert.equal(unbounded.status, 2);
  assert.match(unbounded.stderr, /--test-concurrency=N with N at most CPU_SLOTS \(2\)/);
  assert.equal(wide.status, 2);
  assert.equal(fs.existsSync(path.join(q.dir, 'wide.log')), false);
  assert.equal(fits.status, 0, fs.readFileSync(path.join(q.dir, 'fits.log'), 'utf8'));
});

test('a run inside another run\'s command goes ahead at once on the outer run\'s slots', { skip, timeout: 20000 }, async () => {
  const q = queue(1);
  const inner = [script, 'run', path.join(q.dir, 'inner.log'), process.execPath, q.job, 'inner'];
  const outer = spawn(process.execPath, [script, 'run', path.join(q.dir, 'outer.log'), process.execPath, ...inner], { env: { ...q.env, CPU_SLOTS: '1' }, stdio: 'ignore' });
  assert.equal(await new Promise((resolve) => outer.on('exit', resolve)), 0);
  assert.deepEqual(q.lines(), ['start inner', 'end inner']);
  assert.match(fs.readFileSync(path.join(q.dir, 'inner.log'), 'utf8'), /on the outer run's 1 CPU slot\n/);
});

test('status lists the runs holding slots and those waiting by their labels, and forgets a dead run', { skip, timeout: 20000 }, async (t) => {
  const q = queue(1);
  t.after(() => q.release('a'));
  const status = () => spawnSync(process.execPath, [script, 'status'], { env: q.env, encoding: 'utf8' }).stdout;
  fs.mkdirSync(path.join(q.dir, 'queue', 'runs'), { recursive: true });
  const dead = path.join(q.dir, 'queue', 'runs', '99999999.json');
  fs.writeFileSync(dead, JSON.stringify({ pid: 99999999, want: 1, slots: [0], queued: 0, started: 0, cwd: '/', command: 'gone' }));
  const a = q.start('a', 1, { held: true, labels: { GPU_WHO: 'tester', CPU_JOB: 'first job' } });
  await q.until(() => q.lines().includes('start a'));
  const b = q.start('b', 1);
  await pause(300);
  const text = status();
  assert.match(text, /: 1 of 1 in use\n/);
  assert.match(text, /\n- tester: first job; 1 slot, pid \d+, since \d\d:\d\d:\d\d\n/);
  assert.match(text, /\nWaiting:\n- [^\n]*job\.mjs b; 1 slot, pid \d+, queued at /);
  assert.doesNotMatch(text, /gone/);
  assert.equal(fs.existsSync(dead), false);
  q.release('a');
  assert.deepEqual(await Promise.all([a, b]), [0, 0]);
  assert.match(status(), /: 0 of 1 in use\n\nNo run holds a slot\.$/m);
});
