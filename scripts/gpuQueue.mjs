// The GPU queue of this Mac, for every session, by level and then in order
// of arrival: macOS aborts the command buffers of node WebGPU models that
// share the GPU, and Dawn never reports it, so heavy jobs run alone. Levels:
//   1  a short job (under 5 minutes) that tolerates contention, such as a GPU
//      test file: up to 5 level-1 jobs share the GPU at once;
//   2  a measurement of the GPU, such as pace, of at most 10 minutes of heavy
//      use: alone, and ahead of level 1 once first in level 2 for 10 minutes;
//   3  a long job with heavy GPU use, the default: alone.
//   node scripts/gpuQueue.mjs run LOG CMD [ARG...]
//   node scripts/gpuQueue.mjs share LOG CMD [ARG...]
//   node scripts/gpuQueue.mjs lane [reset]
//   node scripts/gpuQueue.mjs status
//   node scripts/gpuQueue.mjs serve [PORT]
// gpuQueue.zsh wraps them as gpu_run, gpu_share, gpu_lane, gpu_status,
// gpu_serve. `run` waits its turn, holds the queue's lock for the whole of
// CMD (so a driver's gaps between segments stay covered), appends CMD's
// output to LOG between a "start HH:MM:SS" and an "exit RC HH:MM:SS" line,
// then counts the GPU command-buffer aborts the system log recorded
// meanwhile, and exits with CMD's status. A LOG that already holds "exit 0"
// is skipped. A run from inside the holder's own process tree goes ahead at
// once. Before CMD starts, the holder waits while a node process outside the
// queue holds the Metal driver or a page server listens on port 8000. From
// before CMD starts until the next holder has one of its own, it holds an
// idle WebGPU device, so the polling gates of older scripts see the GPU busy
// throughout, gaps and handoffs included. The holder and CMD run with
// GPU_QUEUE_HELD=1 in their environment, without which
// js/gpu/device.module.js opens no GPU on this Mac. CPU-only jobs do not use
// the queue and never wait for it.
// GPU_PRIORITY (1-3, default 3) sets a run's level; `share` runs at level 1,
// and 4 and 5, the levels of older scripts, run as 3. A waiter tries the
// lock only when no live ticket is ahead of it in (level, seq, ms, pid)
// order: a level waits until every lower-numbered one is empty, runs in
// order of arrival, and never preempts a running job. A level-2 waiter that
// has stood first in level 2 for 10 minutes is promoted: it renames its
// ticket p1-0-<ms>-<pid>.json, which older copies of this script also sort
// ahead of every level-1 ticket, so no level-1 job joins the GPU after it,
// and takes the lock alone once the jobs holding it finish. A level-3 waiter
// takes a free lock only once it has stayed free for 15 s with no ticket
// arriving ahead of it meanwhile, so a long job does not take the GPU in the
// gap between another session's runs; a finished holder keeps its idle
// device for a waiter of levels 1 and 2 only. A level-1 job that runs past 5
// minutes, or a level-2 job past 10, says so in its log and the status.
// Level-1 runs hold the GPU together, up to 5 at once, never beside a run of
// level 2 or 3; one joining others takes the GPU at once and treats their
// processes as its own. Each also takes CPU_SLOTS slots (1 by default) of
// scripts/cpuQueue.mjs, unless it runs inside a cpu_run, and gives the GPU
// back to the queue while it waits for them. A level-1 run writes its "exit"
// line only after its abort scan: a command-buffer abort in a node process
// during its job fails it with status 75 and drops level 1 to one job at a
// time until `lane reset`. `lane` prints how many level-1 jobs may share.
// GPU_WHO (a short name for the session) and GPU_JOB (one line on what CMD
// is) label a run. <dir>/status.md lists the holder and the queue in order,
// each with its level, by these labels, or by worktree and command for a run
// without them. A run rewrites it when it queues, takes the lock, starts and
// exits, the holder at each of its checks and the waiters every 30 s;
// `status` rewrites it and prints the same text.
// A run that held the lock appends a line to <dir>/history.jsonl once its
// abort scan ends: who, what, priority, log, pid, queued, started, ended, rc,
// aborts and, for a promoted run, promoted.
// The file keeps its last 500 lines once it passes 1000. `serve` answers on
// 127.0.0.1:PORT (8008) with scripts/gpuQueue.html, a page that polls
// /status.json: the status as JSON, reused for 10 ticks once computed, the
// last 20 runs, and the GPU's use by process: every 25 ticks it samples the
// driver's per-process counters (the AppUsage of ioreg's AGXDeviceUserClient
// entries) and keeps an hour of each process's share of the GPU's time, the
// holder's job as one line whichever job it is. It never opens the lock, a
// ticket or a GPU device and writes nothing, so no run ever waits on it.
// The lock is an O_EXLOCK flock on <dir>/lock, which the kernel drops when
// its holders exit; CMD inherits it as fd 19, so a job outliving a killed
// holder keeps it. A shared run holds an O_SHLOCK on it instead and an
// O_EXLOCK on one of <dir>/share-<i>, i below the limit, which its job
// inherits too, with its CPU slots. Its record is <dir>/shared/<pid>.json,
// and holder.json mirrors one shared record for runs of older copies of this
// script, which take the lock exclusively and so wait for every shared run.
// The lane is <dir>/lane.json. Waiters hold flocks on their tickets in
// <dir>/queue, named p<level>-<seq>-<ms>-<pid>.json; a ticket named
// <seq>-<ms>-<pid>.json is level 3, and a promoted ticket's JSON keeps
// priority 2 and adds promoted, the time of its promotion.
// GPU_QUEUE_DIR (~/.cache/webgcm-gpu-queue), GPU_QUEUE_TICK_MS (200; every
// wait is a multiple of it), GPU_QUEUE_PORT (8000; 0 for none),
// GPU_QUEUE_PROBE (a shell command that prints the pids using the GPU, in
// place of lsof), GPU_QUEUE_DEVICE (off, or a module whose getDevice stands
// in for js/gpu/device.module.js), GPU_QUEUE_COURTESY_MS (15000: how long
// level 3 waits for the lock to stay free), GPU_QUEUE_OVERTIME_MS
// (300000,600000: how long a job of level 1 and of level 2 runs before it
// says it is over), GPU_QUEUE_PROMOTE_MS (600000: how long a level-2 waiter
// stands first in its level before it is promoted), GPU_QUEUE_IOREG (a shell
// command whose output stands in for ioreg's) and GPU_QUEUE_ABORTS
// (off, or a shell command whose output stands in for the system log's)
// serve tests.
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import http from 'node:http';
import { createHash } from 'node:crypto';
import { spawn, spawnSync } from 'node:child_process';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { config as cpuConfig, holdSlots } from './cpuQueue.mjs';

const O_SHLOCK = 0x10, O_EXLOCK = 0x20;
const { O_RDONLY, O_WRONLY, O_CREAT, O_EXCL, O_NONBLOCK } = fs.constants;
const INHERITED_FD = 19;
const AGENT = /(^|[\s/])gpuQueue\.mjs (run|share)\b/;
const TICKET = /^(?:p([1-5])-)?(\d+)-(\d+)-(\d+)\.json$/;
const LEVELS = [1, 2, 3, 4, 5], DEFAULT_LEVEL = 3, SHARED_LEVEL = 1, COURTESY_MS = 15000;
const OVERTIME_MS = [5 * 60e3, 10 * 60e3], PROMOTE_MS = 10 * 60e3, PROMOTED_FROM = 2;
const LSTART = /^\s*(\d+)\s+(\d+)\s+(\w{3} \w{3} [ \d]\d \d\d:\d\d:\d\d \d{4})\s+(.*)$/;
const USAGE = 'usage: gpuQueue.mjs run LOG CMD [ARG...] | gpuQueue.mjs share LOG CMD [ARG...] | gpuQueue.mjs lane [reset] | gpuQueue.mjs status | gpuQueue.mjs serve [PORT]';
const LANE_MOST = 5, ABORTED = 75;
const ABORT_LINE = /^(\d{4})-(\d\d)-(\d\d) (\d\d):(\d\d):(\d\d)(?:\.(\d+))?\s+\S+\s+(.+?)\[(\d+):/;
const PAGE = new URL('./gpuQueue.html', import.meta.url);
const LOCAL_HOST = /^(127\.0\.0\.1|localhost|\[::1\])(:\d+)?$/i;
const SERVE_PORT = 8008, RESERVED_PORTS = [8000, 8001];
const HISTORY_SHOWN = 20, HISTORY_MAX = 1000, HISTORY_KEEP = 500;
const IOREG = ['-r', '-c', 'AGXDeviceUserClient', '-l'], MIN_SHARE = 0.005, LINE_PEAK = 0.01, USE_LINES = 7, USE_WINDOW_MS = 60000;

const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));
const pause = (ms) => Atomics.wait(new Int32Array(new SharedArrayBuffer(4)), 0, 0, ms);
const pad = (n) => String(n).padStart(2, '0');
const clock = (ms = Date.now()) => new Date(ms).toTimeString().slice(0, 8);
const day = (ms = Date.now()) => new Date(ms).toDateString().slice(4, 10);
const when = (ms) => (day(ms) === day() ? clock(ms) : `${day(ms)} ${clock(ms)}`);
const minutes = (ms) => (ms % 60000 ? `${ms / 1000} s` : `${ms / 60000} minutes`);
const cut = (text = '', n = 100) => (text.length > n ? `${text.slice(0, n - 1)}…` : text);
const oneLine = (text) => String(text ?? '').replace(/\s+/g, ' ').trim();
const quote = (argv) => argv.map((a) => (/^[\w@%+=:,./-]+$/.test(a) ? a : `'${a.replace(/'/g, `'\\''`)}'`)).join(' ');
const pidList = (text) => [...new Set(text.split(/\s+/).filter((s) => /^\d+$/.test(s)).map(Number))].filter((pid) => pid > 1);

export function config(env = process.env) {
  const tick = Number(env.GPU_QUEUE_TICK_MS) > 0 ? Number(env.GPU_QUEUE_TICK_MS) : 200;
  const dir = path.resolve(env.GPU_QUEUE_DIR || path.join(os.homedir(), '.cache', 'webgcm-gpu-queue'));
  const device = env.GPU_QUEUE_DEVICE === 'off' ? null
    : env.GPU_QUEUE_DEVICE ? pathToFileURL(path.resolve(env.GPU_QUEUE_DEVICE)).href
    : new URL('../js/gpu/device.module.js', import.meta.url).href;
  const courtesy = String(env.GPU_QUEUE_COURTESY_MS ?? '').trim() ? Number(env.GPU_QUEUE_COURTESY_MS) : NaN;
  const overtime = String(env.GPU_QUEUE_OVERTIME_MS ?? '').split(',').map((ms) => (ms.trim() ? Number(ms) : NaN));
  const promote = Number(env.GPU_QUEUE_PROMOTE_MS);
  return {
    dir, lock: path.join(dir, 'lock'), holder: path.join(dir, 'holder.json'), queue: path.join(dir, 'queue'),
    shared: path.join(dir, 'shared'), lane: path.join(dir, 'lane.json'),
    status: path.join(dir, 'status.md'), statusLock: path.join(dir, 'status.lock'), history: path.join(dir, 'history.jsonl'),
    tick, recheck: 10 * tick, cache: 10 * tick, settle: 15 * tick, grace: 100 * tick, linger: 100 * tick, useEvery: 25 * tick, useKeep: 18000 * tick,
    monitor: 150 * tick, nestedEvery: 25, refreshEvery: 150, stale: 300 * tick,
    port: env.GPU_QUEUE_PORT ? Number(env.GPU_QUEUE_PORT) : env.GPU_QUEUE_PORT === '' ? 0 : 8000,
    probe: env.GPU_QUEUE_PROBE || '', ioreg: env.GPU_QUEUE_IOREG || '', device, aborts: env.GPU_QUEUE_ABORTS !== 'off',
    abortLog: env.GPU_QUEUE_ABORTS && env.GPU_QUEUE_ABORTS !== 'off' ? env.GPU_QUEUE_ABORTS : '',
    who: oneLine(env.GPU_WHO), what: oneLine(env.GPU_JOB), priority: oneLine(env.GPU_PRIORITY) || String(DEFAULT_LEVEL),
    courtesy: courtesy >= 0 ? courtesy : COURTESY_MS,
    overtime: overtime.length === 2 && overtime.every((ms) => ms > 0) ? overtime : OVERTIME_MS,
    promote: promote > 0 ? promote : PROMOTE_MS,
  };
}

const delayOf = (cfg, level) => (level >= DEFAULT_LEVEL ? cfg.courtesy : 0);
const levelOf = (info) => (LEVELS.includes(info?.priority) ? info.priority : DEFAULT_LEVEL);

export function whoFrom(cwd = '') {
  const worktree = /\/\.claude\/worktrees\/([^/]+)/.exec(cwd) || /--claude-worktrees-([^/]+)/.exec(cwd);
  return worktree ? worktree[1] : path.basename(cwd) || '?';
}

function words(command) {
  const out = [];
  let word = null;
  for (let i = 0; i < command.length; i++) {
    const c = command[i];
    if (/\s/.test(c)) {
      if (word !== null) out.push(word);
      word = null;
      continue;
    }
    word ??= '';
    if (c === "'") {
      const end = command.indexOf("'", i + 1);
      word += command.slice(i + 1, end < 0 ? undefined : end);
      i = end < 0 ? command.length : end;
    } else {
      word += c === '\\' && i + 1 < command.length ? command[++i] : c;
    }
  }
  return word === null ? out : [...out, word];
}

export function shortCommand(command = '') {
  let argv = words(command);
  if (path.basename(argv[0] || '') === 'env') {
    argv = argv.slice(1);
    while (argv.length > 1 && /^[A-Za-z_]\w*=/.test(argv[0])) argv.shift();
  }
  if (argv[0]?.startsWith('/')) argv[0] = path.basename(argv[0]);
  return cut(argv.filter((a) => !a.startsWith('--test-reporter'))
    .map((a) => a.replace(/^([^/]*=)?\/[^=]*?\/([^/=]+\/[^/=]+)$/, '$1…/$2'))
    .map((a) => (!a || /[\s'"\\]/.test(a) ? `'${a.replace(/'/g, `'\\''`)}'` : a))
    .join(' '), 120);
}

const labels = (info) => ({ who: oneLine(info.who) || whoFrom(info.cwd), what: oneLine(info.what) || shortCommand(info.command) });
const named = (info) => {
  const { who, what } = labels(info);
  return `${who}, ${what}`;
};

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

const familyOf = (table, holder) => (holder ? new Set([...lineage(table, holder.pid).keys(), ...descendants(table, holder.pid)]) : new Set());

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

export function parseIoreg(text) {
  const times = new Map();
  for (const entry of text.split('+-o ').slice(1)) {
    const creator = /"IOUserClientCreator" = "pid (\d+), ([^"]*)"/.exec(entry);
    if (!creator) continue;
    const pid = Number(creator[1]);
    const ns = [...entry.matchAll(/"accumulatedGPUTime"=(\d+)/g)].reduce((sum, m) => sum + Number(m[1]), 0);
    times.set(pid, { pid, name: creator[2], ns: (times.get(pid)?.ns ?? 0) + ns });
  }
  return times;
}

export const gpuTimes = (cfg = config()) => parseIoreg(cfg.ioreg ? output('/bin/sh', ['-c', cfg.ioreg]) : output('/usr/sbin/ioreg', IOREG));

export function gpuShares(before, after, table, family = new Set()) {
  const window = (after.at - before.at) * 1e6;
  const users = [];
  for (const t of after.times.values()) {
    const ns = Math.max(0, t.ns - (before.times.get(t.pid)?.ns ?? 0));
    if (!(window > 0) || ns < MIN_SHARE * window) continue;
    const p = table.get(t.pid);
    users.push({
      pid: t.pid, name: t.name, what: p && t.name === 'node' ? shortCommand(p.command) : t.name,
      command: p?.command ?? null, share: ns / window, queued: family.has(t.pid),
    });
  }
  return users.sort((a, b) => b.share - a.share || a.pid - b.pid);
}

export function gpuRecorder(cfg, { read = () => gpuTimes(cfg), table = processes, records = (procs) => agents(cfg, procs), now = Date.now } = {}) {
  const samples = [];
  let before = null;
  const sample = () => {
    const at = now(), procs = table();
    const holders = records(procs).filter((r) => r?.start && procs.get(r.pid)?.start === r.start).map((r) => ({ pid: r.pid, ...labels(r) }));
    const after = { at, times: read() };
    if (before && at > before.at && at - before.at <= USE_WINDOW_MS) {
      const family = new Set(holders.flatMap((h) => [...familyOf(procs, h)]));
      samples.push({ since: before.at, until: at, holders, users: gpuShares(before, after, procs, family) });
      while (samples[0].until < at - cfg.useKeep) samples.shift();
    }
    before = after;
  };
  return { samples, sample };
}

const thousandths = (x) => Math.round(x * 1000) / 1000;

export function gpuHistory(samples) {
  const n = samples.length;
  const holder = { shares: samples.map((s) => (s.holders.length ? thousandths(s.users.reduce((sum, u) => sum + (u.queued ? u.share : 0), 0)) : null)), jobs: [] };
  samples.forEach((s, i) => s.holders.forEach((h) => {
    const job = holder.jobs.findLast((j) => j.pid === h.pid);
    if (job && job.to === i - 1) job.to = i;
    else holder.jobs.push({ from: i, to: i, ...h });
  }));
  const found = new Map();
  samples.forEach((s, i) => s.users.forEach((u) => {
    if (u.queued) return;
    const key = `${u.pid} ${u.name}`;
    const line = found.get(key) ?? found.set(key, { key, pid: u.pid, name: u.name, what: u.what, outside: u.name === 'node', first: i, last: i, peak: 0, at: [] }).get(key);
    Object.assign(line, { last: i, peak: Math.max(line.peak, u.share) });
    line.at[i] = u.share;
  }));
  const dense = (line) => Array.from({ length: n }, (_, i) => (i < line.first || i > line.last ? null : thousandths(line.at[i] ?? 0)));
  const ranked = [...found.values()].filter((line) => line.peak >= LINE_PEAK).sort((a, b) => b.peak - a.peak || a.first - b.first);
  const kept = ranked.length > USE_LINES ? ranked.slice(0, USE_LINES - 1) : ranked;
  const lines = kept.sort((a, b) => a.first - b.first || a.pid - b.pid)
    .map((line) => ({ key: line.key, pid: line.pid, name: line.name, what: line.what, outside: line.outside, shares: dense(line) }));
  const rest = ranked.slice(kept.length).map(dense);
  if (rest.length) {
    const shares = Array.from({ length: n }, (_, i) => (rest.some((v) => v[i] !== null) ? thousandths(rest.reduce((sum, v) => sum + (v[i] ?? 0), 0)) : null));
    lines.push({ key: 'other', pid: null, name: `${rest.length} other processes`, what: null, outside: false, shares });
  }
  return { times: samples.map((s) => s.until), holder, lines };
}

const listeners = (port) => (port ? pidList(output('/usr/sbin/lsof', ['-w', '-nP', `-iTCP:${port}`, '-sTCP:LISTEN', '-t'])) : []);
const lockHolders = (cfg) => pidList(output('/usr/sbin/lsof', ['-w', '-t', '--', cfg.lock]));

function tryLock(file, mode = O_EXLOCK) {
  try {
    return fs.openSync(file, O_RDONLY | O_CREAT | mode | O_NONBLOCK, 0o644);
  } catch (error) {
    if (error.code === 'EAGAIN') return null;
    throw error;
  }
}

function takeShare(cfg) {
  const { limit } = laneOf(cfg);
  for (let i = 0; i < limit; i++) {
    const slot = tryLock(path.join(cfg.dir, `share-${i}`));
    if (slot === null) continue;
    const lock = tryLock(cfg.lock, O_SHLOCK);
    if (lock !== null) return { lockFd: lock, slot: { fd: slot, index: i } };
    fs.closeSync(slot);
    break;
  }
  return { lockFd: null, slot: null };
}

function inheritsLock(cfg) {
  try {
    const held = fs.fstatSync(INHERITED_FD), lock = fs.statSync(cfg.lock);
    return held.ino === lock.ino && held.dev === lock.dev;
  } catch {
    return false;
  }
}

const order = (a, b) => a.priority - b.priority || a.seq - b.seq || a.ms - b.ms || a.pid - b.pid;

function tickets(cfg) {
  let names = [];
  try { names = fs.readdirSync(cfg.queue); } catch { return []; }
  return names.map((name) => TICKET.exec(name)).filter(Boolean)
    .map((m) => ({
      name: m[0], file: path.join(cfg.queue, m[0]),
      priority: m[1] ? Number(m[1]) : DEFAULT_LEVEL, seq: Number(m[2]), ms: Number(m[3]), pid: Number(m[4]),
    }))
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
  const dead = readShared(cfg).filter(({ record }) => !liveRecord(record)).map((r) => r.file);
  for (const file of [...leftovers(cfg.queue, /^\./), ...leftovers(cfg.dir, /^(holder\.json|lane\.json|\.status\.md|\.history\.jsonl)\.\d+$/), ...leftovers(cfg.shared, /\.json\.\d+$/), ...dead]) {
    try { if (Date.now() - fs.statSync(file).mtimeMs > cfg.stale) fs.unlinkSync(file); } catch {}
  }
}

function enqueue(cfg, info) {
  fs.mkdirSync(cfg.queue, { recursive: true });
  sweep(cfg);
  const seq = 1 + Math.max(0, ...tickets(cfg).map((t) => t.seq));
  const ms = Date.now();
  const name = `p${info.priority}-${seq}-${ms}-${process.pid}.json`;
  const temp = path.join(cfg.queue, `.${name}`);
  const fd = fs.openSync(temp, O_WRONLY | O_CREAT | O_EXCL | O_EXLOCK | O_NONBLOCK, 0o644);
  fs.writeSync(fd, JSON.stringify(info));
  fs.renameSync(temp, path.join(cfg.queue, name));
  return { fd, name, file: path.join(cfg.queue, name), priority: info.priority, seq, ms, pid: process.pid };
}

function dequeue(ticket) {
  if (!ticket || ticket.fd === null) return;
  fs.rmSync(ticket.file, { force: true });
  fs.closeSync(ticket.fd);
  ticket.fd = null;
}

function promote(ticket, info) {
  const name = `p${SHARED_LEVEL}-0-${ticket.ms}-${ticket.pid}.json`, file = path.join(path.dirname(ticket.file), name);
  const text = JSON.stringify(info);
  // Rewritten in place: a temp file renamed over it would drop the flock that marks the ticket live.
  fs.writeSync(ticket.fd, text, 0);
  fs.ftruncateSync(ticket.fd, Buffer.byteLength(text));
  fs.renameSync(ticket.file, file);
  Object.assign(ticket, { name, file, priority: SHARED_LEVEL, seq: 0 });
}

const readJson = (file) => { try { return JSON.parse(fs.readFileSync(file, 'utf8')); } catch { return null; } };
const readRecord = (cfg) => readJson(cfg.holder);
const liveRecord = (record) => (record && record.start && startOf(record.pid) === record.start ? record : null);

function writeRecord(cfg, record) {
  const temp = `${cfg.holder}.${process.pid}`;
  fs.writeFileSync(temp, JSON.stringify(record));
  fs.renameSync(temp, cfg.holder);
}

function readShared(cfg) {
  let names = [];
  try { names = fs.readdirSync(cfg.shared).filter((name) => /^\d+\.json$/.test(name)); } catch { return []; }
  return names.map((name) => ({ file: path.join(cfg.shared, name), record: readJson(path.join(cfg.shared, name)) })).filter((r) => r.record);
}

function sharedRecords(cfg, table = null) {
  return readShared(cfg).map((r) => r.record)
    .filter((record) => record.start && (table?.get(record.pid)?.start ?? startOf(record.pid)) === record.start)
    .sort((a, b) => a.acquired - b.acquired || a.pid - b.pid);
}

function writeShared(cfg, record) {
  serialized(cfg, () => {
    fs.mkdirSync(cfg.shared, { recursive: true });
    const file = path.join(cfg.shared, `${record.pid}.json`), temp = `${file}.${process.pid}`;
    fs.writeFileSync(temp, JSON.stringify(record));
    fs.renameSync(temp, file);
    const mirror = readRecord(cfg);
    if (!mirror || mirror.pid === record.pid || !liveRecord(mirror)) writeRecord(cfg, record);
  });
}

function dropShared(cfg, pid) {
  serialized(cfg, () => {
    fs.rmSync(path.join(cfg.shared, `${pid}.json`), { force: true });
    const mirror = readRecord(cfg);
    if (mirror && mirror.pid !== pid && liveRecord(mirror)) return;
    const next = sharedRecords(cfg)[0];
    if (next) writeRecord(cfg, next);
    else fs.rmSync(cfg.holder, { force: true });
  });
}

function agents(cfg, table) {
  const record = readRecord(cfg);
  if (record && !record.shared && record.start && (table.get(record.pid)?.start ?? startOf(record.pid)) === record.start) return [record];
  return sharedRecords(cfg, table);
}

function laneOf(cfg) {
  const lane = readJson(cfg.lane);
  return lane?.fallback ? { ...lane, limit: 1 } : { limit: LANE_MOST };
}

function editLane(cfg, edit) {
  return serialized(cfg, () => {
    const lane = edit(laneOf(cfg));
    if (!lane) return null;
    const temp = `${cfg.lane}.${process.pid}`;
    fs.writeFileSync(temp, JSON.stringify(lane));
    fs.renameSync(temp, cfg.lane);
    return lane;
  });
}

export function laneText(cfg = config()) {
  const f = laneOf(cfg).fallback;
  if (f) return `Level 1: one job at a time since ${when(f.at)}, after ${f.aborts} GPU command-buffer abort${f.aborts === 1 ? '' : 's'} in ${f.who}, ${f.what} (${f.log}); gpu_lane reset lets ${LANE_MOST} share again`;
  return `Level 1: up to ${LANE_MOST} jobs share the GPU; levels 2 and 3 run alone`;
}

export function lane(cfg, args = []) {
  if (args.length === 1 && args[0] === 'reset') editLane(cfg, () => ({ limit: LANE_MOST, reset: Date.now() }));
  else if (args.length) throw refused(`not a lane command: ${args.join(' ')}; use gpu_lane, or gpu_lane reset`);
  return laneText(cfg);
}

function holderAbove(cfg, line, probe) {
  const record = readRecord(cfg);
  if (record && line.get(record.pid) === record.start) return record.pid;
  const shared = readShared(cfg).find(({ record: r }) => line.get(r.pid) === r.start);
  if (shared) return shared.record.pid;
  const inherits = inheritsLock(cfg);
  if (!inherits) {
    if (liveRecord(record)) return 0;
    const fd = probe ? tryLock(cfg.lock) : null;
    if (fd !== null) { fs.closeSync(fd); return 0; }
  }
  const holders = new Set(lockHolders(cfg));
  return [...line.keys()].reverse().find((pid) => holders.has(pid)) || (inherits ? process.pid : 0);
}

function survey(cfg, self, share = false) {
  const users = gpuUsers(cfg);
  const table = processes();
  const others = share ? sharedRecords(cfg, table).map((r) => r.pid).filter((pid) => pid !== self) : [];
  const family = new Set([...lineage(table, self).keys(), ...descendants(table, self), ...others.flatMap((pid) => [pid, ...descendants(table, pid)])]);
  const outside = [], queue = [];
  for (const pid of users) {
    const proc = table.get(pid);
    if (!proc || family.has(pid)) continue;
    if (AGENT.test(proc.command)) queue.push(pid);
    else outside.push({ pid, command: proc.command, age: age(proc) });
  }
  return { outside, port: listeners(cfg.port), covered: queue.length > 0, table, users };
}

const pageServer = (cfg, pids) => `the page server on port ${cfg.port} (pid ${pids.join(', ')})`;

function describeBlock(cfg, net) {
  const parts = net.outside.map((u) => `pid ${u.pid} (${shortCommand(u.command)})`);
  if (net.port.length) parts.push(pageServer(cfg, net.port));
  return parts.join(' and ');
}

function holderState(record) {
  const since = record.stateSince ? ` since ${when(record.stateSince)}` : '';
  const over = record.over ? `, past level ${levelOf(record)}'s ${minutes(record.over)}` : '';
  if (record.state === 'running') return `running since ${when(record.started)}${over}`;
  if (record.state === 'waiting') return `waiting${since} for ${record.blockers}`;
  return `checking for GPU users outside the queue${since}`;
}

const described = (p) => ({ pid: p.pid, command: p.command, short: shortCommand(p.command) });

const probing = {
  waiting: (cfg, t) => waiting(cfg, t.file),
  holders(cfg) {
    const free = tryLock(cfg.lock);
    if (free === null) return lockHolders(cfg);
    fs.closeSync(free);
    return [];
  },
};

function watching() {
  let key = null, pids = [];
  return {
    waiting: (cfg, t, table) => (t.info?.start ? table.get(t.pid)?.start === t.info.start : table.has(t.pid)),
    holders(cfg, record, queue) {
      if (!record) return queue.length ? lockHolders(cfg) : [];
      const id = `${record.pid} ${record.start}`;
      if (id !== key) [key, pids] = [id, lockHolders(cfg)];
      return pids;
    },
  };
}

function snapshot(cfg, net = {}, how = probing) {
  const now = Date.now();
  const record = readRecord(cfg);
  const listed = tickets(cfg).map((t) => ({ ...t, info: readJson(t.file) }));
  const table = net.table || processes();
  const live = agents(cfg, table);
  const holder = live[0] && !live[0].shared ? live[0] : null, shared = holder ? [] : live;
  const queue = listed.filter((t) => how.waiting(cfg, t, table) && t.info)
    .map((t) => ({
      ...labels(t.info), priority: t.info.promoted ? levelOf(t.info) : t.priority, promoted: t.info.promoted || undefined,
      pid: t.pid, since: t.ms, log: t.info.log, command: t.info.command, cwd: t.info.cwd, shared: t.info.shared || undefined,
    }));
  const lock = live.length ? [] : how.holders(cfg, record, queue).filter((pid) => table.has(pid)).map((pid) => described(table.get(pid)));
  const family = new Set([...familyOf(table, holder), ...shared.flatMap((r) => [r.pid, ...descendants(table, r.pid)])]);
  const outside = [...(net.users || gpuUsers(cfg))].map((pid) => table.get(pid))
    .filter((p) => p && !family.has(p.pid) && !AGENT.test(p.command))
    .map((p) => ({ ...described(p), start: Date.parse(p.start) }));
  const listening = cfg.port ? net.port || listeners(cfg.port) : [];
  return {
    now,
    holder: holder && { ...holder, ...labels(holder), priority: levelOf(holder) },
    shared: shared.map((r) => ({ ...r, ...labels(r), priority: levelOf(r) })),
    lane: { ...laneOf(cfg), text: laneText(cfg) },
    exited: !live.length && record ? { pid: record.pid, ...labels(record) } : null,
    lock, queue, outside,
    port: cfg.port ? { port: cfg.port, listening: listening.map((pid) => (table.has(pid) ? described(table.get(pid)) : { pid, command: null, short: null })) } : null,
  };
}

const promotedAt = (r) => (r.promoted ? `, promoted at ${when(r.promoted)}` : '');

function statusText(cfg, net) {
  const s = snapshot(cfg, net);
  const out = [`# GPU queue at ${clock(s.now)}, ${day(s.now)}`, ''];
  const h = s.holder;
  if (h) {
    out.push(`Holder: **${h.who}**, ${h.what}`, `- level ${h.priority}${promotedAt(h)}, ${holderState(h)}`,
      `- pid ${h.pid}, took the lock at ${when(h.acquired)}${h.deviceReady ? ', idle device open' : ''}`, `- log ${h.log}`);
  } else if (s.shared.length) {
    out.push(`Holders: ${s.shared.length} level-1 job${s.shared.length === 1 ? '' : 's'} sharing the GPU`);
    for (const r of s.shared) {
      out.push(`- **${r.who}**, ${r.what}`, `  - level ${r.priority}, ${holderState(r)}, pid ${r.pid}${r.deviceReady ? ', idle device open' : ''}`, `  - log ${r.log}`);
    }
  } else if (!s.lock.length) out.push('Holder: none');
  else {
    out.push(`Holder: no live agent; the lock is held by ${s.lock.map((p) => `pid ${p.pid} (${p.short})`).join(', ')}`);
    if (s.exited) out.push(`- the recorded holder, pid ${s.exited.pid} (${s.exited.who}, ${s.exited.what}), has exited`);
  }
  out.push('', s.lane.text, '', s.queue.length ? 'Queue:' : 'Queue: empty');
  s.queue.forEach((t, i) => {
    const indent = ' '.repeat(`${i + 1}. `.length);
    out.push(`${i + 1}. **${t.who}**, ${t.what}`, `${indent}- level ${t.priority}${promotedAt(t)}, queued since ${when(t.since)}, pid ${t.pid}`, `${indent}- log ${t.log}`);
  });
  out.push('', s.outside.length ? 'Outside the queue:' : 'Outside the queue: none');
  for (const p of s.outside) out.push(`- pid ${p.pid} since ${when(p.start)}: ${p.short}`);
  if (s.port) {
    const listener = (p) => `pid ${p.pid}${p.short === null ? '' : ` (${p.short})`}`;
    out.push('', s.port.listening.length ? `Port ${s.port.port}: listening, ${s.port.listening.map(listener).join(', ')}` : `Port ${s.port.port}: free`);
  }
  return out.join('\n');
}

function serialized(cfg, write) {
  fs.mkdirSync(cfg.dir, { recursive: true });
  let fd = tryLock(cfg.statusLock);
  for (let tries = 0; fd === null && tries < 100; tries++) {
    pause(50);
    fd = tryLock(cfg.statusLock);
  }
  try {
    return write();
  } finally {
    if (fd !== null) fs.closeSync(fd);
  }
}

function writeStatus(cfg, net) {
  return serialized(cfg, () => {
    const text = statusText(cfg, net);
    const temp = path.join(cfg.dir, `.status.md.${process.pid}`);
    fs.writeFileSync(temp, `${text}\n`);
    fs.renameSync(temp, cfg.status);
    return text;
  });
}

const refresh = (cfg, net) => { try { writeStatus(cfg, net); } catch {} };

export function appendHistory(cfg, entry) {
  serialized(cfg, () => {
    fs.appendFileSync(cfg.history, `${JSON.stringify(entry)}\n`);
    const lines = fs.readFileSync(cfg.history, 'utf8').split('\n').filter(Boolean);
    if (lines.length <= HISTORY_MAX) return;
    const temp = path.join(cfg.dir, `.history.jsonl.${process.pid}`);
    fs.writeFileSync(temp, `${lines.slice(-HISTORY_KEEP).join('\n')}\n`);
    fs.renameSync(temp, cfg.history);
  });
}

function allRuns(cfg) {
  let text;
  try { text = fs.readFileSync(cfg.history, 'utf8'); } catch { return []; }
  const runs = [];
  for (const line of text.split('\n')) {
    try { if (line) runs.push(JSON.parse(line)); } catch {}
  }
  return runs;
}

export function recentRuns(cfg, count = HISTORY_SHOWN) {
  return allRuns(cfg).reverse().sort((a, b) => (b.ended || 0) - (a.ended || 0)).slice(0, count);
}

export function finished(log) {
  try { return /^exit 0\b/m.test(fs.readFileSync(log, 'utf8')); } catch { return false; }
}

function launch(argv, log, fds) {
  const out = fs.openSync(log, 'a');
  const stdio = ['ignore', out, out];
  if (fds.length) stdio.push(...Array(INHERITED_FD - 3).fill('ignore'), ...fds);
  const child = spawn(argv[0], argv.slice(1), { stdio });
  fs.closeSync(out);
  const done = new Promise((resolve) => {
    child.on('error', (error) => { fs.appendFileSync(log, `${error.message}\n`); resolve(127); });
    child.on('exit', (code, signal) => resolve(code ?? 128 + (os.constants.signals[signal] || 0)));
  });
  return { child, done };
}

function abortEntry(line) {
  const m = ABORT_LINE.exec(line);
  if (!m) return { at: null, process: null, line };
  const [, y, mo, d, h, mi, s, frac = '0'] = m;
  return { at: new Date(+y, mo - 1, +d, +h, +mi, +s, Number(frac.slice(0, 3).padEnd(3, '0'))).getTime(), process: m[8], line };
}

const sharedAborts = (scan, started, ended) => scan.entries
  .filter((e) => (e.process ?? 'node') === 'node' && (e.at === null || (e.at >= started - 1000 && e.at <= ended + 2000))).length;

function scanAborts(cfg, since) {
  const from = new Date(since - 5000);
  const stamp = `${from.getFullYear()}-${pad(from.getMonth() + 1)}-${pad(from.getDate())} ${clock(from.getTime())}`;
  return new Promise((resolve) => {
    const child = cfg.abortLog ? spawn('/bin/sh', ['-c', cfg.abortLog], { stdio: ['ignore', 'pipe', 'ignore'] })
      : spawn('/usr/bin/log', ['show', '--style', 'compact', '--start', stamp, '--predicate', 'eventMessage CONTAINS "command buffer was aborted" AND process != "log"'], { stdio: ['ignore', 'pipe', 'ignore'] });
    let text = '';
    child.stdout.on('data', (chunk) => { text += chunk; });
    const timer = setTimeout(() => child.kill(), 15 * 60 * 1000);
    const finish = (code) => {
      clearTimeout(timer);
      const lines = text.split('\n').filter((l) => l.includes('command buffer was aborted'));
      resolve(code === 0 ? { count: lines.length, lines: [`aborted ${lines.length}`, ...lines.slice(-3).map((l) => cut(l, 300))], entries: lines.map(abortEntry) }
        : { count: null, lines: [`aborted ? (log show exited ${code})`], entries: [] });
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
    if (next && next.pid !== process.pid ? next.deviceReady : !tickets(cfg).some((t) => t.priority < DEFAULT_LEVEL && waiting(cfg, t.file))) break;
  }
  try { device.destroy(); } catch {}
}

export async function run(cfg, logArg, argv, { share: shareCommand = false } = {}) {
  if (!logArg || !argv.length) { console.error(USAGE); return 2; }
  if (!shareCommand && !/^[1-5]$/.test(cfg.priority)) {
    console.error(`gpu_run: GPU_PRIORITY is '${cfg.priority}'; use 1 (under 5 minutes, shares the GPU), 2 (a measurement of at most 10 minutes, alone) or 3 (long, alone; the default)`);
    return 2;
  }
  const priority = shareCommand ? SHARED_LEVEL : Math.min(Number(cfg.priority), DEFAULT_LEVEL);
  const share = priority === SHARED_LEVEL;
  const cpuCfg = cpuConfig();
  const want = share && !/^\d+$/.test(cpuCfg.held) ? Number(cpuCfg.slots) : 0;
  if (share && !/^\d+$/.test(cpuCfg.held) && !(Number.isInteger(want) && want >= 1 && want <= cpuCfg.total)) {
    console.error(`gpu_run: CPU_SLOTS is '${cpuCfg.slots}'; a level-1 run takes a whole number of CPU slots from 1 to ${cpuCfg.total}, 1 by default`);
    return 2;
  }
  const log = path.resolve(logArg);
  if (finished(log)) { console.error(`gpu_run: skipped, ${log} already holds "exit 0"`); return 0; }
  const say = (text) => fs.appendFileSync(log, `${text}\n`);
  fs.mkdirSync(path.dirname(log), { recursive: true });
  if (!shareCommand && Number(cfg.priority) > DEFAULT_LEVEL) say(`level ${cfg.priority} runs as level ${DEFAULT_LEVEL}: GPU_PRIORITY takes 1, 2 or 3`);
  fs.mkdirSync(cfg.queue, { recursive: true });
  const line = lineage(processes(), process.pid);
  const self = {
    pid: process.pid, start: line.get(process.pid) || startOf(process.pid), log, command: quote(argv), cwd: process.cwd(),
    since: Date.now(), who: cfg.who || undefined, what: cfg.what || undefined, priority, shared: share || undefined,
  };

  let phase = 'queued', ticket = null, lockFd = null, slot = null, cpu = null, device = null, job = null, rc = 0, remember = null;
  const letGo = () => {
    if (lockFd !== null) fs.closeSync(lockFd);
    if (slot) fs.closeSync(slot.fd);
    lockFd = null;
    slot = null;
  };
  const release = () => {
    if (lockFd !== null) {
      if (share) dropShared(cfg, process.pid);
      else fs.rmSync(cfg.holder, { force: true });
      letGo();
    }
    cpu?.release();
    cpu = null;
  };
  for (const signal of ['SIGINT', 'SIGTERM', 'SIGHUP']) {
    process.on(signal, () => {
      if (phase === 'running') { job.child.kill(signal); return; }
      const code = 128 + os.constants.signals[signal];
      if (phase !== 'done') {
        dequeue(ticket);
        release();
        say(`exit ${code} ${clock()} (${signal} before the job started)`);
        refresh(cfg);
      }
      remember?.(null);
      try { device?.destroy(); } catch {}
      process.exit(phase === 'done' ? rc : code);
    });
  }

  const waitTurn = async () => {
    ticket = enqueue(cfg, self);
    const delay = delayOf(cfg, priority);
    const at = priority === DEFAULT_LEVEL ? '' : ` at level ${priority}`;
    let freeSince = 0, frontSince = 0;
    for (let n = 1; ; n++) {
      const ahead = tickets(cfg).filter((t) => order(t, ticket) < 0 && waiting(cfg, t.file));
      if (!ahead.length) ({ lockFd, slot } = share ? takeShare(cfg) : { lockFd: tryLock(cfg.lock), slot: null });
      const free = lockFd !== null;
      if (!free) freeSince = 0;
      else if ((share && sharedRecords(cfg).length) || Date.now() - (freeSince ||= Date.now()) >= delay) return 0;
      else letGo();
      if (n === 1) {
        const holder = liveRecord(readRecord(cfg));
        const parts = [holder && `pid ${holder.pid} (${named(holder)})`, ahead.length && `${ahead.length} waiting`].filter(Boolean);
        say(free ? `queued ${clock()}${at} until the GPU has stayed free for ${delay / 1000} s`
          : `queued ${clock()}${at} behind ${parts.join(' and ') || "the lock's holder"}`);
      }
      if (priority === PROMOTED_FROM && !self.promoted) {
        frontSince = ahead.every((t) => t.priority < priority) ? frontSince || Date.now() : 0;
        if (frontSince && Date.now() - frontSince >= cfg.promote) {
          self.promoted = Date.now();
          promote(ticket, self);
          say(`promoted ${clock(self.promoted)} after ${minutes(cfg.promote)} first in level ${priority}: goes ahead of level 1 and runs alone once the jobs holding the GPU now finish`);
          refresh(cfg);
        }
      }
      if (n % cfg.refreshEvery === 1) refresh(cfg);
      const nested = n % cfg.nestedEvery === 0 && holderAbove(cfg, line, false);
      if (nested) return nested;
      await sleep(cfg.tick);
    }
  };
  let above = holderAbove(cfg, line, true) || await waitTurn();
  if (!above && want) {
    const entry = {
      pid: process.pid, want, log, command: self.command, cwd: self.cwd, who: self.who, queued: Date.now(),
      what: `level-1 GPU job${self.what ? `: ${self.what}` : ''}`,
    };
    cpu = await holdSlots(cpuCfg, entry, { wait: false });
    if (!cpu) {
      letGo();
      dequeue(ticket);
      refresh(cfg);
      say(`waiting ${clock()} for a CPU slot, with the GPU left to the queue`);
      cpu = await holdSlots(cpuCfg, entry);
      above = await waitTurn();
    }
  }
  process.env.GPU_QUEUE_HELD = '1';
  if (above) {
    if (ticket) { dequeue(ticket); refresh(cfg); }
    cpu?.release();
    cpu = null;
    phase = 'running';
    say(`start ${clock()} nested in pid ${above}`);
    job = launch(argv, log, inheritsLock(cfg) ? [INHERITED_FD] : []);
    rc = await job.done;
    phase = 'done';
    say(`exit ${rc} ${clock()}`);
    return rc;
  }

  phase = 'holding';
  const record = { ...self, acquired: Date.now(), deviceReady: false, slot: slot?.index };
  const write = () => (share ? writeShared(cfg, record) : writeRecord(cfg, record));
  const enter = (state, extra = {}) => {
    if (record.state !== state) record.stateSince = Date.now();
    Object.assign(record, { state, blockers: undefined }, extra);
    write();
  };
  enter('checking');
  dequeue(ticket);

  let opening = null, clearSince = 0, said = '';
  for (;;) {
    const net = survey(cfg, process.pid, share);
    if (cfg.device && !opening && net.outside.every((u) => u.age >= cfg.grace)) {
      opening = openDevice(cfg).then((d) => {
        device = d;
        record.deviceReady = true;
        write();
      }, (error) => say(`device ${clock()} unavailable: ${error.message}`));
    }
    if (!net.outside.length && !net.port.length) {
      if (said) { said = ''; enter('checking'); }
      refresh(cfg, net);
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
      enter('waiting', { blockers: text });
    }
    refresh(cfg, net);
    await sleep(cfg.recheck);
  }

  const started = Date.now();
  enter('running', { started });
  say(`start ${clock(started)}${share ? `, sharing the GPU in slot ${slot.index + 1} of ${laneOf(cfg).limit}` : ''}`);
  phase = 'running';
  if (cpu) process.env.CPU_QUEUE_HELD = '1';
  job = launch(argv, log, share ? [lockFd, ...(cpu?.fds ?? []), slot.fd] : [lockFd]);
  if (cpu && job.child.pid) cpu.child(job.child.pid);
  refresh(cfg);
  const seen = new Set();
  const overtime = cfg.overtime[priority - 1];
  const monitor = setInterval(() => {
    if (overtime && !record.over && Date.now() - started > overtime) {
      record.over = overtime;
      write();
      say(`overtime ${clock()}: past level ${priority}'s ${minutes(overtime)}; queue jobs this long at level ${priority + 1 === DEFAULT_LEVEL ? DEFAULT_LEVEL : `${priority + 1} or ${DEFAULT_LEVEL}`}`);
    }
    const net = survey(cfg, process.pid, share);
    for (const u of net.outside) {
      if (!seen.has(u.pid)) { seen.add(u.pid); say(`overlap ${clock()} pid ${u.pid}: ${shortCommand(u.command)}`); }
    }
    if (net.port.length && !seen.has('port')) { seen.add('port'); say(`overlap ${clock()} ${pageServer(cfg, net.port)}`); }
    refresh(cfg, net);
  }, cfg.monitor);
  rc = await job.done;
  const ended = Date.now();
  clearInterval(monitor);
  phase = 'done';
  remember = (aborts) => {
    remember = null;
    try { appendHistory(cfg, { ...labels(self), priority, log, pid: process.pid, queued: self.since, started, ended, rc, aborts, shared: share || undefined, promoted: self.promoted }); } catch {}
  };
  release();
  if (!share) say(`exit ${rc} ${clock()}`);
  refresh(cfg);
  const [, scan] = await Promise.all([linger(cfg, device), cfg.aborts ? scanAborts(cfg, started) : { count: null, lines: [], entries: [] }]);
  for (const text of scan.lines) say(text);
  if (share) {
    const caught = sharedAborts(scan, started, ended);
    if (caught) {
      const failed = rc === 0 ? ` with status ${ABORTED}` : '';
      if (rc === 0) rc = ABORTED;
      const fell = editLane(cfg, (lane) => (lane.fallback ? null : { ...lane, limit: 1, fallback: { at: Date.now(), aborts: caught, ...labels(self), log } }));
      say(`${caught} GPU command-buffer abort${caught === 1 ? '' : 's'} in node during the job: the run fails${failed}${fell ? ', and level 1 falls back to one job at a time' : ''}`);
    }
    say(`exit ${rc} ${clock()}`);
  }
  remember(scan.count);
  return rc;
}

export function status(cfg = config()) {
  return writeStatus(cfg);
}

const refused = (message) => Object.assign(new Error(message), { usage: true });

export function serve(cfg, port = SERVE_PORT) {
  if (!Number.isInteger(port) || port < 0 || port > 65535) throw refused(`not a port: ${port}`);
  if (port && (RESERVED_PORTS.includes(port) || port === cfg.port)) {
    throw refused(`port ${port} is reserved: a listener on ${cfg.port || RESERVED_PORTS[0]} holds up every GPU job, and ${RESERVED_PORTS[1]} is the repo's static server`);
  }
  const busy = port ? listeners(port) : [];
  if (busy.length) throw new Error(`port ${port} is in use by pid ${busy.join(', ')}`);
  const page = fs.readFileSync(PAGE, 'utf8');
  const hash = (pattern) => `'sha256-${createHash('sha256').update(pattern.exec(page)?.[1] ?? '').digest('base64')}'`;
  const policy = `default-src 'none'; script-src ${hash(/<script>([^]*?)<\/script>/)}; style-src ${hash(/<style>([^]*?)<\/style>/)}; `
    + "img-src data:; connect-src 'self'; base-uri 'none'; form-action 'none'; frame-ancestors 'none'";
  const how = watching();
  const recorder = gpuRecorder(cfg);
  const sample = () => { try { recorder.sample(); } catch {} };
  sample();
  const sampler = setInterval(sample, cfg.useEvery).unref();
  let cached = null;
  const current = () => {
    if (!cached || Date.now() - cached.at >= cfg.cache) {
      const latest = recorder.samples.at(-1);
      const gpu = latest ? { since: latest.since, until: latest.until, users: latest.users } : { since: null, until: null, users: null };
      const body = JSON.stringify({ ...snapshot(cfg, {}, how), gpu, gpuHistory: gpuHistory(recorder.samples), dir: cfg.dir, history: recentRuns(cfg) });
      cached = { at: Date.now(), body };
    }
    return cached.body;
  };
  const server = http.createServer((req, res) => {
    const send = (code, type, body, headers = {}) => {
      res.writeHead(code, { 'Content-Type': `${type}; charset=utf-8`, 'Cache-Control': 'no-store', 'X-Content-Type-Options': 'nosniff', ...headers });
      res.end(req.method === 'HEAD' ? undefined : body);
    };
    const route = req.url.split('?')[0];
    if (!LOCAL_HOST.test(req.headers.host || '')) send(403, 'text/plain', 'forbidden\n');
    else if (req.method !== 'GET' && req.method !== 'HEAD') send(405, 'text/plain', 'method not allowed\n', { Allow: 'GET, HEAD' });
    else if (route === '/') send(200, 'text/html', page, { 'Content-Security-Policy': policy, 'Referrer-Policy': 'no-referrer' });
    else if (route === '/status.json') {
      try { send(200, 'application/json', current()); } catch (error) { send(500, 'application/json', JSON.stringify({ error: error.message })); }
    } else send(404, 'text/plain', 'not found\n');
  });
  server.on('close', () => clearInterval(sampler));
  return new Promise((resolve, reject) => {
    const fail = (error) => {
      clearInterval(sampler);
      reject(error);
    };
    server.once('error', fail);
    server.listen(port, '127.0.0.1', () => {
      server.off('error', fail);
      resolve(server);
    });
  });
}

if (process.argv[1] && fs.realpathSync(process.argv[1]) === fs.realpathSync(fileURLToPath(import.meta.url))) {
  const [command, ...rest] = process.argv.slice(2);
  if (process.platform !== 'darwin') {
    console.error('gpuQueue: the lock is a macOS O_EXLOCK flock; there is no queue on this platform');
    process.exit(2);
  }
  if (command === 'run') process.exit(await run(config(), rest[0], rest.slice(1)));
  else if (command === 'share') process.exit(await run(config(), rest[0], rest.slice(1), { share: true }));
  else if (command === 'lane') {
    try { console.log(lane(config(), rest)); } catch (error) { console.error(`gpu_lane: ${error.message}`); process.exit(2); }
  } else if (command === 'status') console.log(status());
  else if (command === 'serve') {
    const cfg = config();
    try {
      const server = await serve(cfg, rest[0] === undefined ? SERVE_PORT : /^\d+$/.test(rest[0]) ? Number(rest[0]) : rest[0]);
      console.log(`GPU queue monitor: http://127.0.0.1:${server.address().port}/ (${cfg.dir})`);
    } catch (error) {
      console.error(`gpu_serve: ${error.message}`);
      process.exit(error.usage ? 2 : 1);
    }
  } else { console.error(USAGE); process.exit(2); }
}
