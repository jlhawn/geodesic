import { test } from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import net from 'node:net';
import { spawn, spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { config, gpuUsers, parseLsof, duration } from '../scripts/gpuQueue.mjs';

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
  return { device: { destroy() { fs.appendFileSync(process.env.TRACE, \`destroy \${process.pid} \${Date.now()}\\n\`); } } };
}
`;

const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));
async function waitFor(check, seconds = 30) {
  for (const until = Date.now() + seconds * 1000; !check(); await sleep(10)) {
    if (Date.now() > until) throw new Error(`timed out waiting for ${check}`);
  }
}
const alive = (pid) => { try { process.kill(pid, 0); return true; } catch { return false; } };

function bench(t, extra = {}) {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'gpu-queue-'));
  const gpu = path.join(dir, 'gpu');
  fs.mkdirSync(gpu);
  fs.writeFileSync(path.join(dir, 'job.mjs'), JOB);
  fs.writeFileSync(path.join(dir, 'device.mjs'), DEVICE);
  const env = {
    ...process.env, GPU_QUEUE_DIR: path.join(dir, 'queue'), GPU_QUEUE_TICK_MS: '20', GPU_QUEUE_PORT: '0',
    GPU_QUEUE_PROBE: `ls '${gpu}'`, GPU_QUEUE_DEVICE: path.join(dir, 'device.mjs'), GPU_QUEUE_ABORTS: 'off',
    GPU_QUEUE_NODE: process.execPath, FAKE_GPU: gpu, TRACE: path.join(dir, 'trace'), ...extra,
  };
  const spawned = [];
  t.after(() => {
    for (const child of spawned) { try { process.kill(-child.pid, 'SIGKILL'); } catch {} }
    fs.rmSync(dir, { recursive: true, force: true });
  });
  const start = (file, args) => {
    const child = spawn(file, args, { env, stdio: 'ignore', detached: true });
    child.exited = new Promise((resolve) => child.on('exit', (code, signal) => resolve(code ?? signal)));
    spawned.push(child);
    return child;
  };
  const read = (file) => { try { return fs.readFileSync(path.join(dir, file), 'utf8'); } catch { return ''; } };
  return {
    dir, env, gpu, read,
    run: (log, argv) => start(process.execPath, [script, 'run', path.join(dir, log), ...argv]),
    foreign: (argv) => start(argv[0], argv.slice(1)),
    job: (name, ms, gpuUser = true) => [process.execPath, path.join(dir, 'job.mjs'), name, String(ms), gpuUser ? '1' : '0'],
    status: () => spawnSync('/bin/zsh', ['-c', `source '${wrapper}'; gpu_status`], { env, encoding: 'utf8' }).stdout,
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
  assert.match(status, new RegExp(`lock is held by pid ${b.trace().pid.orphan} `));
  assert.match(status, new RegExp(`holder pid ${holder.pid} .* has exited`));
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
  assert.match(b.status(), /^holder {3}none$/m);
  assert.match(b.status(), /^queue {4}empty$/m);
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
  for (let i = 0; i < 4; i++) b.foreign(b.job(`cpu${i}`, 3000, false));
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

test('status shows the holder and the queue in order', { skip }, async (t) => {
  const b = bench(t);
  const holder = b.run('holder.log', b.job('hold', 4000));
  await waitFor(() => b.trace().start.hold);
  const first = b.run('first.log', b.job('first', 50));
  await waitFor(() => /^queued /m.test(b.read('first.log')));
  const second = b.run('second.log', b.job('second', 50));
  await waitFor(() => /^queued /m.test(b.read('second.log')));
  const status = b.status();
  assert.match(status, new RegExp(`^holder {3}pid ${holder.pid}, running \\d+s; holds the lock since \\d\\d:\\d\\d:\\d\\d \\(\\d+s\\), idle device open$`, 'm'));
  assert.match(status, /^ {9}\S+node \S+job\.mjs hold 4000 1$/m);
  assert.match(status, new RegExp(`^queue {4}1\\. pid ${first.pid}, waiting \\d+s \\(since \\d\\d:\\d\\d:\\d\\d\\): \\S+node \\S+job\\.mjs first 50 1\\n {12}log \\S+first\\.log\\n {9}2\\. pid ${second.pid}, waiting \\d+s .*second 50 1$`, 'm'));
  assert.match(status, /^outside  no GPU users outside the queue$/m);
  assert.deepEqual(await Promise.all([holder.exited, first.exited, second.exited]), [0, 0, 0]);
  assert.match(b.status(), /^holder {3}none\nqueue {4}empty$/m);
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
  assert.match(status, new RegExp(`^holder {3}pid ${agent.pid}, waiting for pid ${outsider.pid} `, 'm'));
  assert.match(status, new RegExp(`^outside  pid ${outsider.pid}, up \\d+s: .*job\\.mjs outsider`, 'm'));
  assert.match(status, new RegExp(`^port {5}\\d+ has a listener: pid ${process.pid}$`, 'm'));
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
  await waitFor(() => /^queued /m.test(b.read('last.log')));
  waiter.kill('SIGTERM');
  assert.equal(await waiter.exited, 143);
  assert.match(b.read('waiter.log'), /^exit 143 \d\d:\d\d:\d\d \(SIGTERM before the job started\)$/m);
  assert.match(b.status(), new RegExp(`^queue {4}1\\. pid ${last.pid},`, 'm'));
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

test('lsof output and durations parse', () => {
  const lsof = 'p100\nntxt\nn/System/Library/Extensions/AGXMetalG16X.bundle/Contents/MacOS/AGXMetalG16X\np200\nn/usr/lib/dyld\np300\nn/System/Library/Extensions/AGXMetal13_3.bundle/x\n';
  assert.deepEqual([...parseLsof(lsof)], [100, 300]);
  assert.deepEqual([duration(4000), duration(125000), duration(3 * 3600e3 + 61e3)], ['4s', '2m05s', '3h01m']);
});
