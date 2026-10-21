import { test } from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import net from 'node:net';
import http from 'node:http';
import { spawn, spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { config, gpuUsers, parseLsof, whoFrom, shortCommand, serve, appendHistory, recentRuns } from '../scripts/gpuQueue.mjs';

const script = fileURLToPath(new URL('../scripts/gpuQueue.mjs', import.meta.url));
const wrapper = fileURLToPath(new URL('../scripts/gpuQueue.zsh', import.meta.url));
const skip = process.platform !== 'darwin' && 'the GPU queue is macOS-only';

const JOB = `import fs from 'node:fs';
const [name, ms, gpu] = process.argv.slice(2);
const mark = \`\${process.env.FAKE_GPU}/\${process.pid}\`;
if (gpu === '1') fs.writeFileSync(mark, '');
fs.appendFileSync(process.env.TRACE, \`start \${name} \${Date.now()} \${process.pid}\\n\`);
await new Promise((resolve) => setTimeout(resolve, Number(ms)));
fs.appendFileSync(process.env.TRACE, \`end \${name} \${Date.now()}\\n\`);
if (gpu === '1') fs.rmSync(mark);
`;
const DEVICE = `import fs from 'node:fs';
export async function getDevice() {
  await new Promise((resolve) => setTimeout(resolve, 50));
  fs.writeFileSync(\`\${process.env.FAKE_GPU}/\${process.pid}\`, '');
  fs.appendFileSync(process.env.TRACE, \`device \${process.pid} \${Date.now()}\\n\`);
  return { device: { destroy() {
    fs.rmSync(\`\${process.env.FAKE_GPU}/\${process.pid}\`, { force: true });
    fs.appendFileSync(process.env.TRACE, \`destroy \${process.pid} \${Date.now()}\\n\`);
  } } };
}
`;
const UNLABELED = `import fs from 'node:fs';
import path from 'node:path';
import { spawnSync } from 'node:child_process';
const [role, dir, cwd, command] = process.argv.slice(2);
const { O_RDONLY, O_WRONLY, O_CREAT, O_EXCL, O_NONBLOCK } = fs.constants;
const start = spawnSync('/bin/ps', ['-o', 'lstart=', '-p', String(process.pid)], { encoding: 'utf8', env: { ...process.env, LC_ALL: 'C' } }).stdout.trim();
const self = { pid: process.pid, start, log: path.join(path.dirname(dir), \`\${role}.log\`), command, cwd, since: Date.now() };
fs.mkdirSync(path.join(dir, 'queue'), { recursive: true });
if (role === 'holder') {
  fs.openSync(path.join(dir, 'lock'), O_RDONLY | O_CREAT | 0x20 | O_NONBLOCK, 0o644);
  const record = { ...self, acquired: Date.now(), state: 'waiting', deviceReady: true, blockers: 'the page server on port 8000 (pid 4242)' };
  fs.writeFileSync(path.join(dir, 'holder.json.tmp'), JSON.stringify(record));
  fs.renameSync(path.join(dir, 'holder.json.tmp'), path.join(dir, 'holder.json'));
} else {
  const name = \`1-\${Date.now()}-\${process.pid}.json\`;
  const temp = path.join(dir, 'queue', \`.\${name}\`);
  const fd = fs.openSync(temp, O_WRONLY | O_CREAT | O_EXCL | 0x20 | O_NONBLOCK, 0o644);
  fs.writeSync(fd, JSON.stringify(self));
  fs.renameSync(temp, path.join(dir, 'queue', name));
}
setInterval(() => {}, 1 << 30);
`;

const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));
async function waitFor(check, seconds = 30) {
  for (const until = Date.now() + seconds * 1000; !check(); await sleep(10)) {
    if (Date.now() > until) throw new Error(`timed out waiting for ${check}`);
  }
}
async function until(check, seconds = 30) {
  for (const deadline = Date.now() + seconds * 1000; ; await sleep(50)) {
    const value = await check();
    if (value) return value;
    if (Date.now() > deadline) throw new Error(`timed out waiting for ${check}`);
  }
}
const alive = (pid) => { try { process.kill(pid, 0); return true; } catch { return false; } };
function get(url, headers = {}, agent = false) {
  return new Promise((resolve, reject) => {
    http.get(url, { headers, agent }, (res) => {
      let body = '';
      res.setEncoding('utf8');
      res.on('data', (chunk) => { body += chunk; });
      res.on('end', () => resolve({ status: res.statusCode, headers: res.headers, body }));
    }).on('error', reject);
  });
}
const json = async (url) => JSON.parse((await get(url)).body);

function bench(t, extra = {}) {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'gpu-queue-'));
  const gpu = path.join(dir, 'gpu');
  fs.mkdirSync(gpu);
  fs.writeFileSync(path.join(dir, 'job.mjs'), JOB);
  fs.writeFileSync(path.join(dir, 'device.mjs'), DEVICE);
  const env = {
    ...process.env, GPU_QUEUE_DIR: path.join(dir, 'queue'), GPU_QUEUE_TICK_MS: '20', GPU_QUEUE_PORT: '0',
    GPU_QUEUE_PROBE: `ls '${gpu}'`, GPU_QUEUE_DEVICE: path.join(dir, 'device.mjs'), GPU_QUEUE_ABORTS: 'off',
    GPU_QUEUE_NODE: process.execPath, GPU_WHO: '', GPU_JOB: '', FAKE_GPU: gpu, TRACE: path.join(dir, 'trace'), ...extra,
  };
  const spawned = [];
  t.after(() => {
    for (const child of spawned) { try { process.kill(-child.pid, 'SIGKILL'); } catch {} }
    fs.rmSync(dir, { recursive: true, force: true });
  });
  const start = (file, args, opts = {}) => {
    const child = spawn(file, args, { env: { ...env, ...opts.env }, cwd: opts.cwd, stdio: 'ignore', detached: true });
    child.exited = new Promise((resolve) => child.on('exit', (code, signal) => resolve(code ?? signal)));
    spawned.push(child);
    return child;
  };
  const read = (file) => { try { return fs.readFileSync(path.join(dir, file), 'utf8'); } catch { return ''; } };
  return {
    dir, env, gpu, read,
    run: (log, argv, opts) => start(process.execPath, [script, 'run', path.join(dir, log), ...argv], opts),
    foreign: (argv) => start(argv[0], argv.slice(1)),
    job: (name, ms, gpuUser = true) => [process.execPath, path.join(dir, 'job.mjs'), name, String(ms), gpuUser ? '1' : '0'],
    status: () => spawnSync('/bin/zsh', ['-c', `source '${wrapper}'; gpu_status`], { env, encoding: 'utf8' }).stdout,
    async monitor() {
      const child = spawn('/bin/zsh', ['-c', `source '${wrapper}'; gpu_serve 0`], { env, stdio: ['ignore', 'pipe', 'inherit'], detached: true });
      spawned.push(child);
      let out = '';
      child.stdout.setEncoding('utf8').on('data', (chunk) => { out += chunk; });
      await waitFor(() => /http:\/\/127\.0\.0\.1:\d+\//.test(out), 10);
      return /http:\/\/127\.0\.0\.1:\d+\//.exec(out)[0];
    },
    file: () => read('queue/status.md'),
    trace() {
      const events = { start: {}, end: {}, pid: {}, device: {}, destroy: {} };
      for (const [kind, key, ms, pid] of read('trace').split('\n').filter(Boolean).map((l) => l.split(' '))) {
        events[kind][key] = Number(ms);
        if (kind === 'start') events.pid[key] = Number(pid);
      }
      return events;
    },
  };
}

test('jobs run one at a time in order of arrival', { skip }, async (t) => {
  const b = bench(t);
  const names = ['a', 'b', 'c', 'd', 'e'];
  const agents = [];
  for (const name of names) {
    agents.push(b.run(`${name}.log`, b.job(name, name === 'a' ? 3000 : 300)));
    await waitFor(() => new RegExp(`^${name === 'a' ? 'start' : 'queued'} `, 'm').test(b.read(`${name}.log`)));
  }
  assert.deepEqual(await Promise.all(agents.map((a) => a.exited)), [0, 0, 0, 0, 0]);
  const { start, end, device, destroy } = b.trace();
  for (let i = 1; i < names.length; i++) {
    assert.ok(start[names[i]] >= end[names[i - 1]], `${names[i]} started before ${names[i - 1]} ended`);
    assert.ok(device[agents[i].pid] <= start[names[i]], `${names[i]} started before its idle device opened`);
    assert.ok(destroy[agents[i - 1].pid] >= device[agents[i].pid], `${names[i - 1]} closed its idle device before ${names[i]} opened one`);
  }
  assert.match(b.read('a.log'), /^start \d\d:\d\d:\d\d\nexit 0 \d\d:\d\d:\d\d\n$/);
  for (const name of names.slice(1)) {
    assert.match(b.read(`${name}.log`), /^queued \d\d:\d\d:\d\d behind (pid \d+ \(.*job\.mjs.*\)( and \d waiting)?|\d waiting)\nstart \d\d:\d\d:\d\d\nexit 0 \d\d:\d\d:\d\d\n$/);
  }
});

test('a driver keeps the GPU through the gaps between its segments', { skip }, async (t) => {
  const b = bench(t);
  const segments = `for i in 1 2 3; do '${process.execPath}' '${path.join(b.dir, 'job.mjs')}' seg$i 150 1; sleep 0.25; done`;
  const driver = b.run('driver.log', ['/bin/zsh', '-c', segments]);
  await waitFor(() => b.trace().start.seg1);
  const after = b.run('after.log', b.job('after', 100));
  const idle = [];
  for (let first = b.trace().device[driver.pid]; !b.trace().end.after; await sleep(5)) {
    const users = fs.readdirSync(b.gpu).map(Number).filter(alive);
    if (!users.length) idle.push(Date.now() - first);
  }
  assert.deepEqual(await Promise.all([driver.exited, after.exited]), [0, 0]);
  const { start, end } = b.trace();
  assert.ok(start.after > end.seg3, 'the queued job started inside the driver');
  assert.deepEqual(idle, [], 'a polling gate would have seen the GPU idle');
});

test('a killed holder keeps the lock until its job ends, then the queue moves on', { skip }, async (t) => {
  const b = bench(t);
  const holder = b.run('holder.log', b.job('orphan', 4000, false));
  await waitFor(() => b.trace().start.orphan);
  const next = b.run('next.log', b.job('next', 100));
  await waitFor(() => /^queued /m.test(b.read('next.log')));
  process.kill(holder.pid, 'SIGKILL');
  await holder.exited;
  const status = b.status();
  assert.match(status, new RegExp(`^Holder: no live agent; the lock is held by (.*, )?pid ${b.trace().pid.orphan} \\(`, 'm'));
  assert.match(status, new RegExp(`^- the recorded holder, pid ${holder.pid} \\(.*job\\.mjs orphan 4000 0\\), has exited$`, 'm'));
  assert.equal(await next.exited, 0);
  const { start, end } = b.trace();
  assert.ok(start.next >= end.orphan, 'the next job started beside the orphaned one');

  const second = b.run('second.log', b.job('doomed', 30000));
  await waitFor(() => b.trace().start.doomed);
  const third = b.run('third.log', b.job('third', 100));
  await waitFor(() => /^queued /m.test(b.read('third.log')));
  process.kill(second.pid, 'SIGKILL');
  await second.exited;
  process.kill(b.trace().pid.doomed, 'SIGKILL');
  const killed = Date.now();
  assert.equal(await third.exited, 0);
  assert.ok(b.trace().start.third - killed < 5000, 'the queue did not recover from a dead holder');
});

test('a stale record naming a live pid is ignored', { skip }, async (t) => {
  const b = bench(t);
  const queue = path.join(b.dir, 'queue');
  fs.mkdirSync(path.join(queue, 'queue'), { recursive: true });
  const fake = { pid: process.pid, start: 'Thu Jan  1 00:00:00 1970', log: '/dev/null', command: 'reused', acquired: 0, state: 'running', started: 0 };
  fs.writeFileSync(path.join(queue, 'holder.json'), JSON.stringify(fake));
  fs.writeFileSync(path.join(queue, 'queue', `1-1-${process.pid}.json`), JSON.stringify(fake));
  assert.match(b.status(), /^Holder: none$/m);
  assert.match(b.status(), /^Queue: empty$/m);
  const agent = b.run('reuse.log', b.job('reuse', 300));
  await waitFor(() => b.trace().start.reuse);
  assert.equal(JSON.parse(fs.readFileSync(path.join(queue, 'holder.json'), 'utf8')).pid, agent.pid);
  assert.equal(await agent.exited, 0);
  assert.doesNotMatch(b.read('reuse.log'), /nested|queued/);
});

test('CPU-only processes neither count as GPU users nor hold up the queue', { skip }, async (t) => {
  const cpu = spawn(process.execPath, ['-e', 'setTimeout(() => {}, 5000)'], { stdio: 'ignore' });
  t.after(() => cpu.kill('SIGKILL'));
  await sleep(200);
  assert.ok(!gpuUsers(config({ ...process.env, GPU_QUEUE_PROBE: '' })).has(cpu.pid));
  const b = bench(t);
  for (let i = 0; i < 4; i++) b.foreign(b.job(`cpu${i}`, 8000, false));
  await waitFor(() => Object.keys(b.trace().start).length === 4);
  const agent = b.run('gpu.log', b.job('gpu', 100));
  assert.equal(await agent.exited, 0);
  assert.equal(b.trace().end.cpu0, undefined, 'the GPU job waited for CPU jobs');
  assert.doesNotMatch(b.read('gpu.log'), /waiting/);
});

test('a run inside the holder goes ahead at once', { skip }, async (t) => {
  const b = bench(t);
  const inner = (name) => `source '${wrapper}'; gpu_run '${path.join(b.dir, `${name}.log`)}' '${process.execPath}' '${path.join(b.dir, 'job.mjs')}' ${name} 100 1`;
  const shell = b.run('outer.log', ['/bin/zsh', '-c', `${inner('inner')} && ${inner('inner2')}`]);
  await waitFor(() => b.trace().start.inner);
  const waiter = b.run('waiter.log', b.job('waiter', 100));
  const viaNode = `require('child_process').spawnSync('/bin/zsh', ['-c', ${JSON.stringify(inner('deep'))}], { stdio: 'inherit' })`;
  assert.equal(await shell.exited, 0);
  const node = b.run('outer2.log', [process.execPath, '-e', viaNode]);
  assert.deepEqual(await Promise.all([waiter.exited, node.exited]), [0, 0]);
  for (const name of ['inner', 'inner2', 'deep']) {
    assert.match(b.read(`${name}.log`), /^start \d\d:\d\d:\d\d nested in pid \d+\n(.*\n)*exit 0 /);
  }
  const { start, end } = b.trace();
  assert.ok(start.waiter >= end.inner2, 'a queued job started inside the holder');
  assert.ok(start.deep >= end.waiter, 'the second outer job jumped the queue');
});

test('the status file names the holder and each waiter by who and what, in order', { skip }, async (t) => {
  const b = bench(t);
  fs.mkdirSync(path.join(b.dir, 'queue'), { recursive: true });
  fs.writeFileSync(path.join(b.dir, 'queue', 'status.md'), 'junk\n'.repeat(3));
  const holder = b.run('hold.log', b.job('hold', 4000), { env: { GPU_WHO: 'alpha', GPU_JOB: 'N=128 page check' } });
  await waitFor(() => b.trace().start.hold);
  const worktree = path.join(b.dir, '.claude', 'worktrees', 'beta-123abc', 'scripts');
  fs.mkdirSync(worktree, { recursive: true });
  const first = b.run('first.log', ['env', 'FOO=1', ...b.job('first', 50)], { cwd: worktree });
  await waitFor(() => /^queued /m.test(b.read('first.log')));
  const second = b.run('second.log', b.job('second', 50), { env: { GPU_WHO: ' gamma ', GPU_JOB: 'two\nlines' } });
  await waitFor(() => /\*\*gamma\*\*/.test(b.file()));
  const file = b.file();
  assert.equal(file.match(/^# GPU queue at \d\d:\d\d:\d\d, \w{3} \d\d$/gm)?.length, 1);
  assert.doesNotMatch(file, /junk/);
  assert.match(file, new RegExp([
    '^Holder: \\*\\*alpha\\*\\*, N=128 page check',
    '- running since \\d\\d:\\d\\d:\\d\\d',
    `- pid ${holder.pid}, took the lock at \\d\\d:\\d\\d:\\d\\d, idle device open`,
    '- log \\S+/hold\\.log',
    '',
    'Queue:',
    '1\\. \\*\\*beta-123abc\\*\\*, node …/gpu-queue-[^/]+/job\\.mjs first 50 1',
    `   - queued since \\d\\d:\\d\\d:\\d\\d, pid ${first.pid}`,
    '   - log \\S+/first\\.log',
    '2\\. \\*\\*gamma\\*\\*, two lines',
    `   - queued since \\d\\d:\\d\\d:\\d\\d, pid ${second.pid}`,
    '   - log \\S+/second\\.log',
    '',
    'Outside the queue: none\n$',
  ].join('\n'), 'm'));
  assert.deepEqual(b.status().split('\n').slice(1), b.file().split('\n').slice(1));
  assert.match(b.read('first.log'), new RegExp(`^queued \\d\\d:\\d\\d:\\d\\d behind pid ${holder.pid} \\(alpha, N=128 page check\\)$`, 'm'));
  assert.match(b.read('second.log'), new RegExp(`^queued \\d\\d:\\d\\d:\\d\\d behind pid ${holder.pid} \\(alpha, N=128 page check\\) and 1 waiting$`, 'm'));
  assert.deepEqual(await Promise.all([holder.exited, first.exited, second.exited]), [0, 0, 0]);
  assert.match(b.file(), /^# GPU queue at .*\n\nHolder: none\n\nQueue: empty\n\nOutside the queue: none\n$/);
});

test('records without labels show their worktree and command, and each run rewrites the file as it moves', { skip }, async (t) => {
  const b = bench(t, { GPU_QUEUE_TICK_MS: '60' });
  const queue = path.join(b.dir, 'queue');
  fs.writeFileSync(path.join(b.dir, 'unlabeled.mjs'), UNLABELED);
  const worktree = path.join(b.dir, '.claude', 'worktrees', 'desktop-multi-profile-c73e68');
  const plain = path.join(b.dir, 'plain');
  fs.mkdirSync(worktree, { recursive: true });
  fs.mkdirSync(plain);
  const page = `${process.execPath} /private/tmp/claude-501/p/s/scratchpad/pageCheck.mjs ${worktree}/runs/fourteen128_day1825.parts.json gpu 8`;
  const figures = `env GPULOCK=/nonexistent 'OCEAN={"everySteps":8}' bash /private/tmp/claude-501/p/s/scratchpad/fifteen/scripts/figures/snapshot.sh 'N=6 ETOPO 2022 check'`;
  const holder = b.foreign([process.execPath, path.join(b.dir, 'unlabeled.mjs'), 'holder', queue, worktree, page]);
  await waitFor(() => fs.existsSync(path.join(queue, 'holder.json')));
  const waiter = b.foreign([process.execPath, path.join(b.dir, 'unlabeled.mjs'), 'waiter', queue, plain, figures]);
  await waitFor(() => fs.readdirSync(path.join(queue, 'queue')).some((n) => n.endsWith(`-${waiter.pid}.json`)));
  const agent = b.run('new.log', b.job('new', 300), { env: { GPU_WHO: 'newcomer', GPU_JOB: 'N=16 parity test' } });
  await waitFor(() => b.file().includes('newcomer'), 6);
  assert.match(b.file(), new RegExp([
    '^# GPU queue at \\d\\d:\\d\\d:\\d\\d, \\w{3} \\d\\d',
    '',
    'Holder: \\*\\*desktop-multi-profile-c73e68\\*\\*, node …/scratchpad/pageCheck\\.mjs …/runs/fourteen128_day1825\\.parts\\.json gpu 8',
    '- waiting for the page server on port 8000 \\(pid 4242\\)',
    `- pid ${holder.pid}, took the lock at \\d\\d:\\d\\d:\\d\\d, idle device open`,
    '- log \\S+/holder\\.log',
    '',
    'Queue:',
    "1\\. \\*\\*plain\\*\\*, bash …/figures/snapshot\\.sh 'N=6 ETOPO 2022 check'",
    `   - queued since \\d\\d:\\d\\d:\\d\\d, pid ${waiter.pid}`,
    '   - log \\S+/waiter\\.log',
    '2\\. \\*\\*newcomer\\*\\*, N=16 parity test',
    `   - queued since \\d\\d:\\d\\d:\\d\\d, pid ${agent.pid}`,
    '   - log \\S+/new\\.log',
    '',
    'Outside the queue: none\n$',
  ].join('\n')));
  const versions = [b.file()];
  const watch = setInterval(() => {
    const text = b.file();
    if (text && text !== versions.at(-1)) versions.push(text);
  }, 5);
  t.after(() => clearInterval(watch));
  process.kill(waiter.pid, 'SIGKILL');
  process.kill(holder.pid, 'SIGKILL');
  assert.equal(await agent.exited, 0);
  clearInterval(watch);
  versions.push(b.file());
  const at = (pattern) => versions.findIndex((v) => pattern.test(v));
  const taken = at(/^Holder: \*\*newcomer\*\*, N=16 parity test\n- checking for GPU users outside the queue since \d\d:\d\d:\d\d\n[^]*\nQueue: empty\n/m);
  const running = at(/^Holder: \*\*newcomer\*\*, N=16 parity test\n- running since \d\d:\d\d:\d\d\n/m);
  assert.ok(taken > 0, 'the file did not show the new holder taking the lock');
  assert.ok(running > taken, 'the file did not show the job start');
  assert.match(versions.at(-1), /\n\nHolder: none\n\nQueue: empty\n/);
  for (const version of versions) assert.equal(version.match(/^# GPU queue at /gm)?.length, 1);
});

test('the holder waits for GPU users outside the queue and for a page server', { skip }, async (t) => {
  const server = net.createServer();
  await new Promise((resolve) => server.listen(0, '127.0.0.1', resolve));
  t.after(() => server.close());
  const b = bench(t, { GPU_QUEUE_PORT: String(server.address().port) });
  const outsider = b.foreign(b.job('outsider', 3000));
  await waitFor(() => b.trace().start.outsider);
  const agent = b.run('blocked.log', b.job('blocked', 100));
  await waitFor(() => /^waiting /m.test(b.read('blocked.log')));
  const status = b.status();
  assert.match(b.read('blocked.log'), new RegExp(`^waiting \\d\\d:\\d\\d:\\d\\d for pid ${outsider.pid} \\(.*job\\.mjs outsider.* and the page server on port \\d+ \\(pid ${process.pid}\\)$`, 'm'));
  assert.match(status, new RegExp(`^- waiting since \\d\\d:\\d\\d:\\d\\d for pid ${outsider.pid} \\(node \\S+job\\.mjs outsider 3000 1\\) and the page server on port \\d+ \\(pid ${process.pid}\\)$`, 'm'));
  assert.match(status, new RegExp(`^- pid ${agent.pid}, took the lock at \\d\\d:\\d\\d:\\d\\d`, 'm'));
  assert.match(status, new RegExp(`^Outside the queue:\\n- pid ${outsider.pid} since \\d\\d:\\d\\d:\\d\\d: node \\S+job\\.mjs outsider 3000 1$`, 'm'));
  assert.match(status, new RegExp(`^Port \\d+: listening, pid ${process.pid}\\b`, 'm'));
  await outsider.exited;
  await sleep(500);
  assert.equal(b.trace().start.blocked, undefined, 'the job started beside the page server');
  const closed = Date.now();
  await new Promise((resolve) => server.close(resolve));
  assert.equal(await agent.exited, 0);
  const { start, end } = b.trace();
  assert.ok(start.blocked >= end.outsider && start.blocked >= closed);
});

test('signals end a waiter cleanly and reach a running job', { skip }, async (t) => {
  const b = bench(t);
  const holder = b.run('holder.log', b.job('hold', 30000));
  await waitFor(() => b.trace().start.hold);
  const waiter = b.run('waiter.log', b.job('never', 50));
  await waitFor(() => /^queued /m.test(b.read('waiter.log')));
  const last = b.run('last.log', b.job('last', 50));
  await waitFor(() => new RegExp(`pid ${last.pid}$`, 'm').test(b.file()));
  assert.match(b.file(), new RegExp(`pid ${waiter.pid}$`, 'm'));
  waiter.kill('SIGTERM');
  assert.equal(await waiter.exited, 143);
  assert.match(b.read('waiter.log'), /^exit 143 \d\d:\d\d:\d\d \(SIGTERM before the job started\)$/m);
  assert.doesNotMatch(b.file(), new RegExp(`pid ${waiter.pid}$`, 'm'));
  assert.match(b.status(), new RegExp(`^Queue:\\n1\\. .*\\n {3}- queued since \\d\\d:\\d\\d:\\d\\d, pid ${last.pid}\\n`, 'm'));
  holder.kill('SIGTERM');
  assert.deepEqual(await Promise.all([holder.exited, last.exited]), [143, 0]);
  assert.match(b.read('holder.log'), /^exit 143 /m);
  const { start, end } = b.trace();
  assert.equal(end.hold, undefined);
  assert.equal(start.never, undefined);
  assert.ok(start.last > start.hold);
});

test('a log that already holds exit 0 is skipped', { skip }, async (t) => {
  const b = bench(t);
  fs.writeFileSync(path.join(b.dir, 'done.log'), 'start 10:00:00\nexit 0 10:05:00\n');
  assert.equal(await b.run('done.log', b.job('again', 50)).exited, 0);
  assert.equal(b.read('trace'), '');
  assert.equal(b.read('done.log'), 'start 10:00:00\nexit 0 10:05:00\n');
});

test('lsof output, worktrees and commands parse', () => {
  const lsof = 'p100\nntxt\nn/System/Library/Extensions/AGXMetalG16X.bundle/Contents/MacOS/AGXMetalG16X\np200\nn/usr/lib/dyld\np300\nn/System/Library/Extensions/AGXMetal13_3.bundle/x\n';
  assert.deepEqual([...parseLsof(lsof)], [100, 300]);
  assert.equal(whoFrom('/Users/j/geodesic/.claude/worktrees/divergence-hybrid-2ec1bd'), 'divergence-hybrid-2ec1bd');
  assert.equal(whoFrom('/Users/j/geodesic/.claude/worktrees/divergence-hybrid-2ec1bd/scripts'), 'divergence-hybrid-2ec1bd');
  assert.equal(whoFrom('/private/tmp/claude-501/-Users-j-geodesic--claude-worktrees-gpu-queue-692ebb/0c9f/scratchpad'), 'gpu-queue-692ebb');
  assert.equal(whoFrom('/Users/j/geodesic'), 'geodesic');
  const node = '/Users/j/.nvm/versions/node/v22.22.2/bin/node';
  assert.equal(shortCommand(`${node} --test --test-reporter=tap '--test-name-pattern=on the GPU matches the CPU' test/divergenceSurfaces.test.mjs`),
    "node --test '--test-name-pattern=on the GPU matches the CPU' test/divergenceSurfaces.test.mjs");
  assert.equal(shortCommand(`env GPULOCK=/nonexistent PATH=/usr/bin:/bin 'OCEAN={"everySteps":8}' bash /private/tmp/s/scratchpad/fifteen/scripts/figures/snapshot.sh /private/tmp/s/n6/n6_day0003.bin 'N=6 ETOPO 2022 check'`),
    "bash …/figures/snapshot.sh …/n6/n6_day0003.bin 'N=6 ETOPO 2022 check'");
  assert.equal(shortCommand(`${node} x.mjs --state=/a/b/c/d.bin --out=/a/b 'it'\\''s'`), "node x.mjs --state=…/c/d.bin --out=/a/b 'it'\\''s'");
});

test('the monitor serves its page and the status as JSON', { skip }, async (t) => {
  const b = bench(t);
  const url = await b.monitor();
  const page = await get(url);
  assert.equal(page.status, 200);
  assert.match(page.headers['content-type'], /^text\/html/);
  assert.match(page.headers['content-security-policy'], /script-src 'sha256-[\w+/]+=*'; style-src 'sha256-[\w+/]+=*'/);
  assert.match(page.body, /<title>GPU queue<\/title>/);
  assert.match(page.body, /fetch\('status\.json'/);
  assert.doesNotMatch(page.body, /\b(src|href)=["']?(https?:)?\/\/|@import|url\(/);
  const status = await get(`${url}status.json`);
  assert.equal(status.status, 200);
  assert.match(status.headers['content-type'], /^application\/json/);
  const s = JSON.parse(status.body);
  assert.ok(Math.abs(s.now - Date.now()) < 5000);
  assert.deepEqual({ ...s, now: 0 }, { now: 0, holder: null, exited: null, lock: [], queue: [], outside: [], port: null, dir: path.join(b.dir, 'queue'), history: [] });
  assert.equal((await get(url, { host: 'gpu.example.com' })).status, 403);
  assert.equal((await get(`${url}queue/holder.json`)).status, 404);
  assert.equal(fs.existsSync(path.join(b.dir, 'queue')), false, 'the monitor wrote to the queue directory');
});

test('the monitor shows a labeled job as the holder, then in the history with its exit code', { skip }, async (t) => {
  const b = bench(t);
  const url = `${await b.monitor()}status.json`;
  const alpha = b.run('alpha.log', [process.execPath, '-e', 'setTimeout(() => process.exit(3), 2000)'], { env: { GPU_WHO: 'alpha', GPU_JOB: 'labeled check' } });
  const running = await until(async () => { const s = await json(url); return s.holder?.state === 'running' && s.holder; });
  assert.deepEqual([running.who, running.what, running.pid, running.log], ['alpha', 'labeled check', alpha.pid, path.join(b.dir, 'alpha.log')]);
  assert.ok(running.started >= running.acquired && running.acquired >= running.since);
  const beta = b.run('beta.log', b.job('beta', 50), { env: { GPU_WHO: 'beta', GPU_JOB: 'second in line' } });
  const [waiter] = await until(async () => { const { queue } = await json(url); return queue.length && queue; });
  assert.deepEqual([waiter.who, waiter.what, waiter.pid, waiter.log], ['beta', 'second in line', beta.pid, path.join(b.dir, 'beta.log')]);
  assert.equal((await json(url)).holder.pid, alpha.pid);
  assert.deepEqual(await Promise.all([alpha.exited, beta.exited]), [3, 0]);
  const done = await until(async () => { const s = await json(url); return s.history.length === 2 && s; });
  assert.equal(done.holder, null);
  const [second, first] = done.history;
  assert.deepEqual({ ...first, queued: 0, started: 0, ended: 0 },
    { who: 'alpha', what: 'labeled check', log: path.join(b.dir, 'alpha.log'), pid: alpha.pid, queued: 0, started: 0, ended: 0, rc: 3, aborts: null });
  assert.ok(first.queued <= first.started && first.ended - first.started >= 2000);
  assert.deepEqual([second.who, second.what, second.pid, second.rc, second.aborts], ['beta', 'second in line', beta.pid, 0, null]);
  assert.ok(second.started >= first.ended);
  assert.deepEqual(b.read('queue/history.jsonl').trim().split('\n').map((line) => JSON.parse(line)).sort((x, y) => x.ended - y.ended), [first, second]);
});

test('polling the monitor never opens the lock or a ticket, so no waiter waits on it', { skip }, async (t) => {
  const b = bench(t);
  const cfg = config(b.env);
  const touched = [];
  for (const name of ['openSync', 'writeFileSync', 'appendFileSync', 'renameSync', 'rmSync', 'unlinkSync', 'mkdirSync']) {
    const original = fs[name];
    fs[name] = function (file, ...rest) {
      if (typeof file === 'string' && path.resolve(file).startsWith(cfg.dir)) touched.push(`${name} ${file}`);
      return original.call(this, file, ...rest);
    };
    t.after(() => { fs[name] = original; });
  }
  const server = await serve(cfg, 0);
  const agent = new http.Agent({ keepAlive: true });
  const url = `http://127.0.0.1:${server.address().port}/status.json`;
  let polling = true;
  const codes = [];
  const pollers = Array.from({ length: 4 }, async () => {
    for (; polling; await sleep(5)) codes.push(await get(url, {}, agent).then((res) => res.status, (error) => error.code));
  });
  t.after(async () => {
    polling = false;
    await Promise.all(pollers);
    agent.destroy();
    server.close();
  });

  const holder = b.run('holder.log', b.job('orphan', 2500, false));
  await waitFor(() => b.trace().start.orphan);
  const first = b.run('first.log', b.job('first', 200));
  await waitFor(() => /^queued /m.test(b.read('first.log')));
  process.kill(holder.pid, 'SIGKILL');
  await holder.exited;
  const stuck = await until(async () => { const s = await json(url); return s.lock.length && s; });
  assert.equal(stuck.holder, null);
  assert.deepEqual(stuck.lock.map((p) => p.pid), [b.trace().pid.orphan]);
  assert.equal(stuck.exited.pid, holder.pid);
  assert.deepEqual(stuck.queue.map((w) => w.pid), [first.pid]);
  const second = b.run('second.log', b.job('second', 200));
  assert.deepEqual(await Promise.all([first.exited, second.exited]), [0, 0]);
  polling = false;
  await Promise.all(pollers);
  assert.deepEqual(touched, [], 'the monitor opened or wrote a file in the queue directory');
  assert.ok(codes.length > 20 && codes.every((code) => code === 200), `${codes.length} polls, codes ${[...new Set(codes)]}`);
  const { start, end } = b.trace();
  assert.ok(start.first >= end.orphan && start.second >= end.first);
  assert.ok(start.first - end.orphan < 3000 && start.second - end.first < 3000, 'a waiter was slow to take its turn');
});

test('the monitor refuses the ports the queue and the page use, and a busy one', { skip }, async (t) => {
  const b = bench(t, { GPU_QUEUE_PORT: '9123' });
  const serveOn = (port) => spawnSync(process.execPath, [script, 'serve', port], { env: b.env, encoding: 'utf8', timeout: 10000 });
  for (const port of ['8000', '8001', '9123']) {
    const refused = serveOn(port);
    assert.equal(refused.status, 2);
    assert.match(refused.stderr, new RegExp(`^gpu_serve: port ${port} is reserved`));
  }
  assert.match(serveOn('http').stderr, /^gpu_serve: not a port: http$/m);
  const busy = net.createServer();
  await new Promise((resolve) => busy.listen(0, '127.0.0.1', resolve));
  t.after(() => busy.close());
  const taken = serveOn(String(busy.address().port));
  assert.equal(taken.status, 1);
  assert.match(taken.stderr, new RegExp(`^gpu_serve: port ${busy.address().port} is in use by pid ${process.pid}$`, 'm'));
});

test('the history keeps its last 500 runs once it passes 1000', { skip }, (t) => {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'gpu-queue-history-'));
  t.after(() => fs.rmSync(dir, { recursive: true, force: true }));
  const cfg = config({ GPU_QUEUE_DIR: dir });
  fs.writeFileSync(cfg.history, Array.from({ length: 999 }, (_, n) => `${JSON.stringify({ n })}\n`).join(''));
  appendHistory(cfg, { n: 999 });
  assert.equal(fs.readFileSync(cfg.history, 'utf8').split('\n').filter(Boolean).length, 1000);
  appendHistory(cfg, { n: 1000 });
  const kept = fs.readFileSync(cfg.history, 'utf8').split('\n');
  assert.equal(kept.pop(), '');
  assert.deepEqual(kept.map((line) => JSON.parse(line).n), Array.from({ length: 500 }, (_, i) => 501 + i));
  assert.deepEqual(fs.readdirSync(dir).sort(), ['history.jsonl', 'status.lock']);
  assert.deepEqual(recentRuns(cfg, 3).map((run) => run.n), [1000, 999, 998]);
  appendHistory(cfg, { n: 'long scan', ended: 2000 });
  appendHistory(cfg, { n: 'short scan', ended: 3000 });
  fs.appendFileSync(cfg.history, '{"n":1001,"cut sh');
  assert.deepEqual(recentRuns(cfg, 3).map((run) => run.n), ['short scan', 'long scan', 1000]);
});
