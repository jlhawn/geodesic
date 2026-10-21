// One GPU job at a time on this Mac, in order of arrival, for every session:
// macOS aborts the command buffers of node WebGPU models that share the GPU,
// and Dawn never reports it.
//   node scripts/gpuQueue.mjs run LOG CMD [ARG...]
//   node scripts/gpuQueue.mjs status
// scripts/gpuQueue.zsh wraps both as gpu_run and gpu_status. `run` waits its
// turn, holds the queue's lock for the whole of CMD (so a driver's gaps
// between segments stay covered), appends CMD's output to LOG between a
// "start HH:MM:SS" and an "exit RC HH:MM:SS" line, then counts the GPU
// command-buffer aborts the system log recorded meanwhile, and exits with
// CMD's status. A LOG that already holds "exit 0" is skipped. A run from
// inside the holder's own process tree goes ahead at once. Before CMD
// starts, the holder waits while a node process outside the queue holds the
// Metal driver or a page server listens on port 8000. From before CMD starts
// until the next holder has one of its own, it holds an idle WebGPU device,
// so the polling gates of older scripts see the GPU busy throughout, gaps
// and handoffs included. CPU-only jobs do not use the queue and never wait
// for it.
// The lock is an O_EXLOCK flock on <dir>/lock, which the kernel drops when
// its holders exit; CMD inherits it as fd 19, so a job outliving a killed
// holder keeps it. Waiters hold flocks on their tickets in <dir>/queue.
// GPU_QUEUE_DIR (~/.cache/webgcm-gpu-queue), GPU_QUEUE_TICK_MS (200; every
// wait is a multiple of it), GPU_QUEUE_PORT (8000; 0 for none),
// GPU_QUEUE_PROBE (a shell command that prints the pids using the GPU, in
// place of lsof), GPU_QUEUE_DEVICE (off, or a module whose getDevice stands
// in for js/gpu/device.module.js) and GPU_QUEUE_ABORTS (off) serve tests.
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { spawn, spawnSync } from 'node:child_process';
import { fileURLToPath, pathToFileURL } from 'node:url';

const O_EXLOCK = 0x20;
const { O_RDONLY, O_WRONLY, O_CREAT, O_EXCL, O_NONBLOCK } = fs.constants;
const INHERITED_FD = 19;
const AGENT = /(^|[\s/])gpuQueue\.mjs run\b/;
const TICKET = /^(\d+)-(\d+)-(\d+)\.json$/;
const LSTART = /^\s*(\d+)\s+(\d+)\s+(\w{3} \w{3} [ \d]\d \d\d:\d\d:\d\d \d{4})\s+(.*)$/;
const USAGE = 'usage: gpuQueue.mjs run LOG CMD [ARG...] | gpuQueue.mjs status';

export function config(env = process.env) {
  const tick = Number(env.GPU_QUEUE_TICK_MS) > 0 ? Number(env.GPU_QUEUE_TICK_MS) : 200;
  const dir = path.resolve(env.GPU_QUEUE_DIR || path.join(os.homedir(), '.cache', 'webgcm-gpu-queue'));
  const device = env.GPU_QUEUE_DEVICE === 'off' ? null
    : env.GPU_QUEUE_DEVICE ? pathToFileURL(path.resolve(env.GPU_QUEUE_DEVICE)).href
    : new URL('../js/gpu/device.module.js', import.meta.url).href;
  return {
    dir, lock: path.join(dir, 'lock'), holder: path.join(dir, 'holder.json'), queue: path.join(dir, 'queue'),
    tick, recheck: 10 * tick, settle: 15 * tick, grace: 100 * tick, linger: 100 * tick,
    monitor: 150 * tick, nestedEvery: 25, stale: 300 * tick,
    port: env.GPU_QUEUE_PORT ? Number(env.GPU_QUEUE_PORT) : env.GPU_QUEUE_PORT === '' ? 0 : 8000,
    probe: env.GPU_QUEUE_PROBE || '', device, aborts: env.GPU_QUEUE_ABORTS !== 'off',
  };
}

const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));
const pad = (n) => String(n).padStart(2, '0');
const clock = (ms = Date.now()) => new Date(ms).toTimeString().slice(0, 8);
const cut = (text = '', n = 100) => (text.length > n ? `${text.slice(0, n - 1)}…` : text);
const brief = (command = '') => cut(command.replace(/^\S*\//, ''), 120);
const quote = (argv) => argv.map((a) => (/^[\w@%+=:,./-]+$/.test(a) ? a : `'${a.replace(/'/g, `'\\''`)}'`)).join(' ');
const pidList = (text) => [...new Set(text.split(/\s+/).filter((s) => /^\d+$/.test(s)).map(Number))].filter((pid) => pid > 1);

export function duration(ms) {
  const s = Math.max(0, Math.round(ms / 1000));
  if (s < 60) return `${s}s`;
  if (s < 3600) return `${Math.floor(s / 60)}m${pad(s % 60)}s`;
  return `${Math.floor(s / 3600)}h${pad(Math.floor(s / 60) % 60)}m`;
}

function output(file, args) {
  return spawnSync(file, args, { encoding: 'utf8', env: { ...process.env, LC_ALL: 'C' }, maxBuffer: 1 << 26 }).stdout || '';
}

export function processes() {
  const table = new Map();
  for (const line of output('/bin/ps', ['-axww', '-o', 'pid=,ppid=,lstart=,command=']).split('\n')) {
    const m = LSTART.exec(line);
    if (m) table.set(Number(m[1]), { pid: Number(m[1]), ppid: Number(m[2]), start: m[3], command: m[4] });
  }
  return table;
}

const startOf = (pid) => output('/bin/ps', ['-o', 'lstart=', '-p', String(pid)]).trim() || null;
const age = (proc) => Date.now() - Date.parse(proc.start);

export function lineage(table, pid) {
  const line = new Map();
  for (let p = table.get(pid); p && !line.has(p.pid); p = table.get(p.ppid)) {
    line.set(p.pid, p.start);
    if (p.ppid <= 1) break;
  }
  return line;
}

export function descendants(table, pid) {
  const children = new Map();
  for (const p of table.values()) children.set(p.ppid, [...(children.get(p.ppid) || []), p.pid]);
  const found = new Set();
  for (const stack = [pid]; stack.length;) {
    for (const child of children.get(stack.pop()) || []) if (!found.has(child)) { found.add(child); stack.push(child); }
  }
  return found;
}

export function parseLsof(text) {
  const pids = new Set();
  let pid = 0;
  for (const line of text.split('\n')) {
    if (line.startsWith('p')) pid = Number(line.slice(1));
    else if (line.startsWith('n') && line.includes('AGXMetal') && pid > 0) pids.add(pid);
  }
  return pids;
}

export function gpuUsers(cfg = config()) {
  if (cfg.probe) return new Set(pidList(output('/bin/sh', ['-c', cfg.probe])));
  return parseLsof(output('/usr/sbin/lsof', ['-w', '-c', 'node', '-a', '-d', 'txt', '-Fpn']));
}

const listeners = (port) => (port ? pidList(output('/usr/sbin/lsof', ['-w', '-nP', `-iTCP:${port}`, '-sTCP:LISTEN', '-t'])) : []);
const lockHolders = (cfg) => pidList(output('/usr/sbin/lsof', ['-w', '-t', '--', cfg.lock]));

function tryLock(cfg) {
  try {
    return fs.openSync(cfg.lock, O_RDONLY | O_CREAT | O_EXLOCK | O_NONBLOCK, 0o644);
  } catch (error) {
    if (error.code === 'EAGAIN') return null;
    throw error;
  }
}

function inheritsLock(cfg) {
  try {
    const held = fs.fstatSync(INHERITED_FD), lock = fs.statSync(cfg.lock);
    return held.ino === lock.ino && held.dev === lock.dev;
  } catch {
    return false;
  }
}

const order = (a, b) => a.seq - b.seq || a.ms - b.ms || a.pid - b.pid;

function tickets(cfg) {
  let names = [];
  try { names = fs.readdirSync(cfg.queue); } catch { return []; }
  return names.map((name) => TICKET.exec(name)).filter(Boolean)
    .map((m) => ({ name: m[0], file: path.join(cfg.queue, m[0]), seq: Number(m[1]), ms: Number(m[2]), pid: Number(m[3]) }))
    .sort(order);
}

function waiting(cfg, file) {
  let fd;
  try {
    fd = fs.openSync(file, O_RDONLY | O_EXLOCK | O_NONBLOCK);
  } catch (error) {
    return error.code === 'EAGAIN';
  }
  try { if (Date.now() - fs.fstatSync(fd).mtimeMs > cfg.stale) fs.unlinkSync(file); } catch {}
  fs.closeSync(fd);
  return false;
}

function sweep(cfg) {
  const leftovers = (dir, pattern) => { try { return fs.readdirSync(dir).filter((n) => pattern.test(n)).map((n) => path.join(dir, n)); } catch { return []; } };
  for (const file of [...leftovers(cfg.queue, /^\./), ...leftovers(cfg.dir, /^holder\.json\.\d+$/)]) {
    try { if (Date.now() - fs.statSync(file).mtimeMs > cfg.stale) fs.unlinkSync(file); } catch {}
  }
}

function enqueue(cfg, info) {
  fs.mkdirSync(cfg.queue, { recursive: true });
  sweep(cfg);
  const seq = 1 + Math.max(0, ...tickets(cfg).map((t) => t.seq));
  const ms = Date.now();
  const name = `${seq}-${ms}-${process.pid}.json`;
  const temp = path.join(cfg.queue, `.${name}`);
  const fd = fs.openSync(temp, O_WRONLY | O_CREAT | O_EXCL | O_EXLOCK | O_NONBLOCK, 0o644);
  fs.writeSync(fd, JSON.stringify(info));
  fs.renameSync(temp, path.join(cfg.queue, name));
  return { fd, name, file: path.join(cfg.queue, name), seq, ms, pid: process.pid };
}

function dequeue(ticket) {
  if (!ticket || ticket.fd === null) return;
  fs.rmSync(ticket.file, { force: true });
  fs.closeSync(ticket.fd);
  ticket.fd = null;
}

const readJson = (file) => { try { return JSON.parse(fs.readFileSync(file, 'utf8')); } catch { return null; } };
const readRecord = (cfg) => readJson(cfg.holder);
const liveRecord = (record) => (record && record.start && startOf(record.pid) === record.start ? record : null);

function writeRecord(cfg, record) {
  const temp = `${cfg.holder}.${process.pid}`;
  fs.writeFileSync(temp, JSON.stringify(record));
  fs.renameSync(temp, cfg.holder);
}

function holderAbove(cfg, line, probe) {
  const record = readRecord(cfg);
  if (record && line.get(record.pid) === record.start) return record.pid;
  const inherits = inheritsLock(cfg);
  if (!inherits) {
    if (liveRecord(record)) return 0;
    const fd = probe ? tryLock(cfg) : null;
    if (fd !== null) { fs.closeSync(fd); return 0; }
  }
  const holders = new Set(lockHolders(cfg));
  return [...line.keys()].reverse().find((pid) => holders.has(pid)) || (inherits ? process.pid : 0);
}

function survey(cfg, self) {
  const users = gpuUsers(cfg);
  const table = processes();
  const family = new Set([...lineage(table, self).keys(), ...descendants(table, self)]);
  const outside = [], queue = [];
  for (const pid of users) {
    const proc = table.get(pid);
    if (!proc || family.has(pid)) continue;
    if (AGENT.test(proc.command)) queue.push(pid);
    else outside.push({ pid, command: proc.command, age: age(proc) });
  }
  return { outside, port: listeners(cfg.port), covered: queue.length > 0 };
}

const pageServer = (cfg, pids) => `the page server on port ${cfg.port} (pid ${pids.join(', ')})`;

function describeBlock(cfg, net) {
  const parts = net.outside.map((u) => `pid ${u.pid} (${brief(u.command)})`);
  if (net.port.length) parts.push(pageServer(cfg, net.port));
  return parts.join(' and ');
}

export function finished(log) {
  try { return /^exit 0\b/m.test(fs.readFileSync(log, 'utf8')); } catch { return false; }
}

function launch(argv, log, extra) {
  const out = fs.openSync(log, 'a');
  const stdio = ['ignore', out, out];
  if (extra !== null) stdio.push(...Array(INHERITED_FD - 3).fill('ignore'), extra);
  const child = spawn(argv[0], argv.slice(1), { stdio });
  fs.closeSync(out);
  const done = new Promise((resolve) => {
    child.on('error', (error) => { fs.appendFileSync(log, `${error.message}\n`); resolve(127); });
    child.on('exit', (code, signal) => resolve(code ?? 128 + (os.constants.signals[signal] || 0)));
  });
  return { child, done };
}

function scanAborts(since) {
  const from = new Date(since - 5000);
  const stamp = `${from.getFullYear()}-${pad(from.getMonth() + 1)}-${pad(from.getDate())} ${clock(from.getTime())}`;
  return new Promise((resolve) => {
    const child = spawn('/usr/bin/log', ['show', '--style', 'compact', '--start', stamp, '--predicate', 'eventMessage CONTAINS "command buffer was aborted" AND process != "log"'], { stdio: ['ignore', 'pipe', 'ignore'] });
    let text = '';
    child.stdout.on('data', (chunk) => { text += chunk; });
    const timer = setTimeout(() => child.kill(), 15 * 60 * 1000);
    const finish = (code) => {
      clearTimeout(timer);
      const lines = text.split('\n').filter((l) => l.includes('command buffer was aborted'));
      resolve(code === 0 ? [`aborted ${lines.length}`, ...lines.slice(-3).map((l) => cut(l, 300))] : [`aborted ? (log show exited ${code})`]);
    };
    child.on('error', () => finish(-1));
    child.on('close', finish);
  });
}

async function openDevice(cfg) {
  const { getDevice } = await import(cfg.device);
  const { device } = await getDevice();
  return device;
}

async function linger(cfg, device) {
  if (!device) return;
  for (const until = Date.now() + cfg.linger; Date.now() < until; await sleep(cfg.tick)) {
    const next = readRecord(cfg);
    if (next && next.pid !== process.pid ? next.deviceReady : !tickets(cfg).some((t) => waiting(cfg, t.file))) break;
  }
  try { device.destroy(); } catch {}
}

export async function run(cfg, logArg, argv) {
  if (!logArg || !argv.length) { console.error(USAGE); return 2; }
  const log = path.resolve(logArg);
  if (finished(log)) { console.error(`gpu_run: skipped, ${log} already holds "exit 0"`); return 0; }
  const say = (text) => fs.appendFileSync(log, `${text}\n`);
  fs.mkdirSync(path.dirname(log), { recursive: true });
  fs.mkdirSync(cfg.queue, { recursive: true });
  const line = lineage(processes(), process.pid);
  const self = { pid: process.pid, start: line.get(process.pid) || startOf(process.pid), log, command: quote(argv), cwd: process.cwd(), since: Date.now() };

  let phase = 'queued', ticket = null, lockFd = null, device = null, job = null, rc = 0;
  const release = () => {
    if (lockFd === null) return;
    fs.rmSync(cfg.holder, { force: true });
    fs.closeSync(lockFd);
    lockFd = null;
  };
  for (const signal of ['SIGINT', 'SIGTERM', 'SIGHUP']) {
    process.on(signal, () => {
      if (phase === 'running') { job.child.kill(signal); return; }
      const code = 128 + os.constants.signals[signal];
      if (phase !== 'done') {
        dequeue(ticket);
        release();
        say(`exit ${code} ${clock()} (${signal} before the job started)`);
      }
      try { device?.destroy(); } catch {}
      process.exit(phase === 'done' ? rc : code);
    });
  }

  let above = holderAbove(cfg, line, true);
  if (!above) {
    ticket = enqueue(cfg, self);
    for (let n = 1; ; n++) {
      const ahead = tickets(cfg).filter((t) => order(t, ticket) < 0 && waiting(cfg, t.file));
      if (!ahead.length && (lockFd = tryLock(cfg)) !== null) break;
      if (n === 1) {
        const holder = liveRecord(readRecord(cfg));
        const parts = [holder && `pid ${holder.pid} (${brief(holder.command)})`, ahead.length && `${ahead.length} waiting`].filter(Boolean);
        say(`queued ${clock()} behind ${parts.join(' and ') || "the lock's holder"}`);
      }
      if (n % cfg.nestedEvery === 0 && (above = holderAbove(cfg, line, false))) break;
      await sleep(cfg.tick);
    }
  }
  if (above) {
    dequeue(ticket);
    phase = 'running';
    say(`start ${clock()} nested in pid ${above}`);
    job = launch(argv, log, inheritsLock(cfg) ? INHERITED_FD : null);
    rc = await job.done;
    phase = 'done';
    say(`exit ${rc} ${clock()}`);
    return rc;
  }

  phase = 'holding';
  const record = { ...self, acquired: Date.now(), state: 'checking', deviceReady: false };
  writeRecord(cfg, record);
  dequeue(ticket);

  let opening = null, clearSince = 0, said = '';
  for (;;) {
    const net = survey(cfg, process.pid);
    if (cfg.device && !opening && net.outside.every((u) => u.age >= cfg.grace)) {
      opening = openDevice(cfg).then((d) => {
        device = d;
        record.deviceReady = true;
        writeRecord(cfg, record);
      }, (error) => say(`device ${clock()} unavailable: ${error.message}`));
    }
    if (!net.outside.length && !net.port.length) {
      await opening;
      if (net.covered || (clearSince && Date.now() - clearSince >= cfg.settle)) break;
      clearSince ||= Date.now();
      await sleep(cfg.settle);
      continue;
    }
    clearSince = 0;
    const text = describeBlock(cfg, net);
    if (text !== said) {
      say(`waiting ${clock()} for ${text}`);
      said = text;
      Object.assign(record, { state: 'waiting', blockers: text });
      writeRecord(cfg, record);
    }
    await sleep(cfg.recheck);
  }

  const started = Date.now();
  Object.assign(record, { state: 'running', started, blockers: undefined });
  writeRecord(cfg, record);
  say(`start ${clock(started)}`);
  phase = 'running';
  job = launch(argv, log, lockFd);
  const seen = new Set();
  const monitor = setInterval(() => {
    const net = survey(cfg, process.pid);
    for (const u of net.outside) {
      if (!seen.has(u.pid)) { seen.add(u.pid); say(`overlap ${clock()} pid ${u.pid}: ${brief(u.command)}`); }
    }
    if (net.port.length && !seen.has('port')) { seen.add('port'); say(`overlap ${clock()} ${pageServer(cfg, net.port)}`); }
  }, cfg.monitor);
  rc = await job.done;
  clearInterval(monitor);
  phase = 'done';
  release();
  say(`exit ${rc} ${clock()}`);
  const [, aborts] = await Promise.all([linger(cfg, device), cfg.aborts ? scanAborts(started) : []]);
  for (const text of aborts) say(text);
  return rc;
}

export function status(cfg = config()) {
  const now = Date.now();
  const table = processes();
  const out = [`GPU queue ${cfg.dir.replace(os.homedir(), '~')}, ${clock(now)}`];
  const record = readRecord(cfg);
  const holder = record && table.get(record.pid)?.start === record.start ? record : null;
  if (holder) {
    const state = holder.state === 'running' ? `running ${duration(now - holder.started)}`
      : holder.state === 'waiting' ? `waiting for ${holder.blockers}` : 'checking for other GPU users';
    out.push(`holder   pid ${holder.pid}, ${state}; holds the lock since ${clock(holder.acquired)} (${duration(now - holder.acquired)})${holder.deviceReady ? ', idle device open' : ''}`);
    out.push(`         ${holder.command}`, `         log ${holder.log}`);
  } else {
    const holders = lockHolders(cfg).filter((pid) => table.has(pid));
    const left = record ? `; holder pid ${record.pid} (${brief(record.command)}) has exited` : '';
    out.push(holders.length ? `holder   none recorded, but the lock is held by ${holders.map((pid) => `pid ${pid} (${brief(table.get(pid).command)})`).join(', ')}${left}` : 'holder   none');
  }
  const queue = tickets(cfg).filter((t) => waiting(cfg, t.file));
  if (!queue.length) out.push('queue    empty');
  queue.forEach((t, i) => {
    const info = readJson(t.file) || {};
    out.push(`${i ? '         ' : 'queue    '}${i + 1}. pid ${t.pid}, waiting ${duration(now - t.ms)} (since ${clock(t.ms)}): ${info.command}`);
    out.push(`            log ${info.log}`);
  });
  const family = holder ? new Set([...lineage(table, holder.pid).keys(), ...descendants(table, holder.pid)]) : new Set();
  const outside = [...gpuUsers(cfg)].map((pid) => table.get(pid)).filter((p) => p && !family.has(p.pid) && !AGENT.test(p.command));
  if (!outside.length) out.push('outside  no GPU users outside the queue');
  outside.forEach((p, i) => out.push(`${i ? '         ' : 'outside  '}pid ${p.pid}, up ${duration(age(p))}: ${p.command}`));
  if (cfg.port) {
    const port = listeners(cfg.port);
    out.push(port.length ? `port     ${cfg.port} has a listener: pid ${port.join(', ')}` : `port     ${cfg.port} free`);
  }
  return out.join('\n');
}

if (process.argv[1] && fs.realpathSync(process.argv[1]) === fs.realpathSync(fileURLToPath(import.meta.url))) {
  const [command, ...rest] = process.argv.slice(2);
  if (process.platform !== 'darwin') {
    console.error('gpuQueue: the lock is a macOS O_EXLOCK flock; there is no queue on this platform');
    process.exit(2);
  }
  if (command === 'run') process.exit(await run(config(), rest[0], rest.slice(1)));
  else if (command === 'status') console.log(status());
  else { console.error(USAGE); process.exit(2); }
}
