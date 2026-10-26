// A budget of the Mac's memory for the node jobs that run on its CPU, shared
// by every session: test files and CPU models are memory-bound, and enough
// of them at once exhaust the memory and hang the machine.
//   node scripts/cpuQueue.mjs run LOG CMD [ARG...]
//   node scripts/cpuQueue.mjs status
// gpuQueue.zsh wraps them as cpu_run and cpu_status. The budget is
// CPU_QUEUE_SLOTS slots of 1 GB, by default the machine's memory less 16 GB
// for the system, the apps and the GPU queue's job: 20 on a 36 GB Mac. `run`
// waits until CPU_SLOTS of them (1 by default) are free, holds them for the
// whole of CMD, appends CMD's output to LOG between a "start HH:MM:SS" and an
// "exit RC HH:MM:SS" line and exits with CMD's status. A LOG that already
// holds "exit 0" is skipped. A CMD with node's --test must pass
// --test-concurrency at most its slots. A run inside another run's CMD goes
// ahead at once on the outer run's slots. GPU_WHO and CPU_JOB (one line on
// what CMD is) label a run; `status` lists the runs holding slots and those
// waiting. A shared GPU test file of scripts/gpuQueue.mjs takes one slot
// through holdSlots, without waiting when it holds the GPU.
// Each slot is an O_EXLOCK flock on <dir>/slot-<i>. A waiter takes slots
// only while it holds the flock on <dir>/door, so a run waiting for many is
// never overtaken by smaller ones. CMD inherits the slots from fd 20 on, so a
// job outliving a killed run keeps them; the kernel drops them when their
// last holder exits. <dir>/runs/<pid>.json describes each run.
// CPU_QUEUE_DIR (~/.cache/webgcm-cpu-queue) and CPU_QUEUE_TICK_MS (200)
// serve tests.
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { spawn } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { finished, shortCommand, whoFrom } from './gpuQueue.mjs';

const O_EXLOCK = 0x20;
const { O_RDONLY, O_CREAT, O_NONBLOCK } = fs.constants;
const FIRST_FD = 20, GB = 2 ** 30, RESERVED_GB = 16;
const SIGNALS = ['SIGINT', 'SIGTERM', 'SIGHUP'];
const USAGE = 'usage: cpuQueue.mjs run LOG CMD [ARG...] | cpuQueue.mjs status';

const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));
const clock = (ms = Date.now()) => new Date(ms).toTimeString().slice(0, 8);
const oneLine = (text) => String(text ?? '').replace(/\s+/g, ' ').trim();
const quote = (argv) => argv.map((a) => (/^[\w@%+=:,./-]+$/.test(a) ? a : `'${a.replace(/'/g, `'\\''`)}'`)).join(' ');

export function config(env = process.env) {
  const dir = path.resolve(env.CPU_QUEUE_DIR || path.join(os.homedir(), '.cache', 'webgcm-cpu-queue'));
  return {
    dir, door: path.join(dir, 'door'), runs: path.join(dir, 'runs'),
    total: Number(env.CPU_QUEUE_SLOTS) || Math.max(2, Math.floor(os.totalmem() / GB) - RESERVED_GB),
    tick: Number(env.CPU_QUEUE_TICK_MS) || 200,
    slots: oneLine(env.CPU_SLOTS) || '1', held: oneLine(env.CPU_QUEUE_HELD),
    who: oneLine(env.GPU_WHO), what: oneLine(env.CPU_JOB),
  };
}

export function testConcurrency(argv) {
  if (!argv.includes('--test')) return null;
  const at = argv.findIndex((a) => /^--test-concurrency(=|$)/.test(a));
  if (at < 0) return 0;
  const value = argv[at].includes('=') ? argv[at].slice(argv[at].indexOf('=') + 1) : argv[at + 1];
  return /^\d+$/.test(value ?? '') ? Number(value) : 0;
}

function tryLock(file) {
  try {
    return fs.openSync(file, O_RDONLY | O_CREAT | O_EXLOCK | O_NONBLOCK, 0o644);
  } catch (error) {
    if (error.code === 'EAGAIN') return null;
    throw error;
  }
}

function alive(pid) {
  try {
    process.kill(pid, 0);
    return true;
  } catch (error) {
    return error.code === 'EPERM';
  }
}

function save(cfg, entry) {
  const file = path.join(cfg.runs, `${entry.pid}.json`), temp = `${file}.tmp`;
  fs.writeFileSync(temp, JSON.stringify(entry));
  fs.renameSync(temp, file);
}

export function runs(cfg = config()) {
  let names;
  try { names = fs.readdirSync(cfg.runs).filter((name) => /^\d+\.json$/.test(name)); } catch { return []; }
  const live = [];
  for (const name of names) {
    const file = path.join(cfg.runs, name);
    let entry;
    try { entry = JSON.parse(fs.readFileSync(file, 'utf8')); } catch { continue; }
    if (alive(entry.pid) || (entry.child && alive(entry.child))) live.push(entry);
    else fs.rmSync(file, { force: true });
  }
  return live.sort((a, b) => a.queued - b.queued);
}

export function status(cfg = config()) {
  const all = runs(cfg), holding = all.filter((r) => r.slots), waiting = all.filter((r) => !r.slots);
  const used = holding.reduce((sum, r) => sum + r.slots.length, 0);
  const line = (r, at) => `- ${r.who || whoFrom(r.cwd)}: ${r.what || shortCommand(r.command)}; ${r.want} slot${r.want === 1 ? '' : 's'}, pid ${r.child || r.pid}, ${at} ${clock(r.started || r.queued)}`;
  return [
    `# CPU slots at ${clock()}: ${used} of ${cfg.total} in use`, '',
    ...(holding.length ? holding.map((r) => line(r, 'since')) : ['No run holds a slot.']),
    ...(waiting.length ? ['', 'Waiting:', ...waiting.map((r) => line(r, 'queued at'))] : []),
  ].join('\n');
}

function launch(argv, log, fds, held, started = () => {}) {
  const out = fs.openSync(log, 'a');
  const stdio = ['ignore', out, out];
  if (fds.length) stdio.push(...Array(FIRST_FD - 3).fill('ignore'), ...fds);
  const child = spawn(argv[0], argv.slice(1), { stdio, env: { ...process.env, CPU_QUEUE_HELD: String(held) } });
  fs.closeSync(out);
  if (child.pid) started(child);
  const forward = (signal) => { if (child.exitCode === null && child.signalCode === null) child.kill(signal); };
  SIGNALS.forEach((signal) => process.on(signal, forward));
  return new Promise((resolve) => {
    child.on('error', (error) => { fs.appendFileSync(log, `${error.message}\n`); resolve(127); });
    child.on('exit', (code, signal) => resolve(code ?? 128 + (os.constants.signals[signal] || 0)));
  }).finally(() => SIGNALS.forEach((signal) => process.off(signal, forward)));
}

export async function run(cfg, logArg, argv) {
  if (!logArg || !argv.length) { console.error(USAGE); return 2; }
  const nested = /^\d+$/.test(cfg.held), want = Number(nested ? cfg.held : cfg.slots);
  if (!nested && !(Number.isInteger(want) && want >= 1 && want <= cfg.total)) {
    console.error(`cpu_run: CPU_SLOTS is '${cfg.slots}'; use a whole number from 1 to ${cfg.total}, the GB the job takes at its peak`);
    return 2;
  }
  const concurrency = testConcurrency(argv);
  if (concurrency !== null && !(concurrency >= 1 && concurrency <= want)) {
    console.error(`cpu_run: node --test runs a file per core at once by default; pass --test-concurrency=N with N at most CPU_SLOTS (${want})`);
    return 2;
  }
  const log = path.resolve(logArg);
  if (finished(log)) { console.error(`cpu_run: skipped, ${log} already holds "exit 0"`); return 0; }
  fs.mkdirSync(path.dirname(log), { recursive: true });
  const say = (text) => fs.appendFileSync(log, `${text}\n`);
  if (nested) {
    say(`start ${clock()}, on the outer run's ${want} CPU slot${want === 1 ? '' : 's'}`);
    const rc = await launch(argv, log, [], want);
    say(`exit ${rc} ${clock()}`);
    return rc;
  }

  const entry = {
    pid: process.pid, want, log, command: quote(argv), cwd: process.cwd(),
    who: cfg.who || undefined, what: cfg.what || undefined, queued: Date.now(),
  };
  const slots = await holdSlots(cfg, entry);
  const waited = Math.round((entry.started - entry.queued) / 1000);
  say(`start ${clock(entry.started)}, ${want} of ${cfg.total} CPU slots${waited ? ` after ${waited} s waiting` : ''}`);
  const rc = await launch(argv, log, slots.fds, want, (child) => slots.child(child.pid));
  say(`exit ${rc} ${clock()}`);
  slots.release();
  return rc;
}

export async function holdSlots(cfg, entry, { wait = true } = {}) {
  fs.mkdirSync(cfg.runs, { recursive: true });
  if (wait) save(cfg, entry);
  const held = new Map();
  let door;
  while ((door = tryLock(cfg.door)) === null) {
    if (!wait) return null;
    await sleep(cfg.tick);
  }
  for (;;) {
    for (let i = 0; i < cfg.total && held.size < entry.want; i++) {
      if (held.has(i)) continue;
      const fd = tryLock(path.join(cfg.dir, `slot-${i}`));
      if (fd !== null) held.set(i, fd);
    }
    if (held.size === entry.want) break;
    if (!wait) {
      for (const fd of [...held.values(), door]) fs.closeSync(fd);
      return null;
    }
    await sleep(cfg.tick);
  }
  fs.closeSync(door);
  entry.slots = [...held.keys()];
  entry.started = Date.now();
  save(cfg, entry);
  return {
    fds: [...held.values()],
    child(pid) { entry.child = pid; save(cfg, entry); },
    release() {
      for (const fd of held.values()) fs.closeSync(fd);
      held.clear();
      fs.rmSync(path.join(cfg.runs, `${entry.pid}.json`), { force: true });
    },
  };
}

if (process.argv[1] && fs.realpathSync(process.argv[1]) === fs.realpathSync(fileURLToPath(import.meta.url))) {
  const [command, ...rest] = process.argv.slice(2);
  if (process.platform !== 'darwin') {
    console.error('cpuQueue: the slots are macOS O_EXLOCK flocks; there is no budget on this platform');
    process.exit(2);
  }
  if (command === 'run') process.exit(await run(config(), rest[0], rest.slice(1)));
  else if (command === 'status') console.log(status());
  else { console.error(USAGE); process.exit(2); }
}
