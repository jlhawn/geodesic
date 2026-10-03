import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, mkdtempSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { spawnSync } from 'node:child_process';
import { createEnergyRecord, balanceLine, RECORD_DAYS, FIELDS } from '../js/physics/energyRecord.module.js';
import { decodeState, encodeState } from '../js/stateFile.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const root = new URL('..', import.meta.url).pathname;

const entry = (day) => ({ asr: 240 + 10 * Math.sin((2 * Math.PI * day) / 365), olr: 236, ts: 288 + day / 1000 });

test('the ring holds the last 365 consecutive days, the windows mean over them, and a whole year cancels the seasonal cycle', () => {
  const r = createEnergyRecord();
  assert.equal(r.values.length, 2 + RECORD_DAYS * FIELDS.length);
  assert.equal(r.windows().day, null);
  assert.equal(balanceLine(r), 'day —, week —, month —, year —');
  for (let d = 1; d <= 400; d++) r.add(d, entry(d));
  assert.equal(r.count, 365);
  assert.equal(r.newest, 400);
  const w = r.windows();
  assert.deepEqual([w.day.n, w.week.n, w.month.n, w.year.n], [1, 7, 30, 365]);
  assert.deepEqual([w.year.from, w.year.to, w.week.from], [36, 400, 394]);
  assert.ok(Math.abs(w.year.net - 4) < 1e-9, `the year's net is the constant part, got ${w.year.net}`);
  assert.ok(Math.abs(w.day.net - (entry(400).asr - 236)) < 1e-12);
  assert.ok(Math.abs(w.month.ts - (288 + 385.5 / 1000)) < 1e-12);
  assert.match(balanceLine(r), /^day [+-]\d+\.\d, week [+-]\d+\.\d, month [+-]\d+\.\d, year \+4\.0$/);
});

test('a day short of a window is reported, a day re-added replaces its entry, and a day out of sequence restarts the record', () => {
  const r = createEnergyRecord();
  for (let d = 1; d <= 10; d++) r.add(d, entry(d));
  assert.equal(r.count, 10);
  assert.match(balanceLine(r), /month [+-]\d+\.\d \(10 of 30 d\), year [+-]\d+\.\d \(10 of 365 d\)/);
  r.add(10, { asr: 300, olr: 236, ts: 290 });
  assert.equal(r.count, 10);
  assert.equal(r.windows().day.asr, 300);
  r.add(8, entry(8));
  assert.deepEqual([r.count, r.newest], [8, 8]);
  r.add(9, entry(9));
  assert.deepEqual([r.count, r.newest], [9, 9]);
  r.add(20, entry(20));
  assert.deepEqual([r.count, r.newest], [1, 20]);
  r.add(21, entry(21));
  assert.equal(r.count, 2);
  r.add(22, null);
  r.add(23, entry(23));
  assert.deepEqual([r.count, r.newest, r.windows().week.n, r.windows().week.from], [4, 23, 3, 20]);
  assert.match(balanceLine(r), /week [+-]\d+\.\d \(3 of 7 d\)/);
  assert.throws(() => r.add(0, entry(1)), /a model day/);
  assert.throws(() => r.add(2.5, entry(1)), /a model day/);
});

test('the record goes through the state file as float64 and loads back whole', async () => {
  const r = createEnergyRecord();
  for (let d = 1; d <= 3; d++) r.add(d, entry(d));
  const saved = await decodeState(encodeState({ N: 4, day: 3, energyRecord: r.values, pi: new Float32Array([1, 2]) }));
  assert.ok(saved.energyRecord instanceof Float64Array);
  const back = createEnergyRecord();
  back.load(saved.energyRecord);
  assert.deepEqual([...back.values], [...r.values]);
  assert.throws(() => back.load(new Float64Array(5)), /an energy record of/);
});

test('the spin-up adds each whole day to the record, carries it across a split, and prints the balance', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const out = mkdtempSync(join(tmpdir(), 'energy-'));
  const run = (env) => spawnSync(process.execPath, [join(root, 'scripts/spinup.mjs')], { cwd: root, env: { ...process.env, N: '16', MINUTES: '100', OUT: out, ...env }, encoding: 'utf8' });
  let r = run({ TAG: 'whole', DAYS: '3' });
  assert.equal(r.status, 0, r.stderr.slice(-1500));
  const whole = await decodeState(new Uint8Array(readFileSync(join(out, 'whole_day0003.bin'))));
  assert.ok(whole.energyRecord instanceof Float64Array);
  assert.deepEqual([whole.energyRecord[0], whole.energyRecord[1]], [3, 3]);
  const log = readFileSync(join(out, 'whole.log'), 'utf8');
  const lines = log.split('\n').filter((l) => /^day \d/.test(l));
  assert.equal(lines.length, 3);
  assert.match(lines[0], /; balance day [+-]\d+\.\d, week [+-]\d+\.\d \(1 of 7 d\), month [+-]\d+\.\d \(1 of 30 d\), year [+-]\d+\.\d \(1 of 365 d\) W\/m²$/);
  const net = Number(lines[2].match(/; balance day ([+-]\d+\.\d)/)[1]);
  const asr = Number(lines[2].match(/ASR ([\d.]+)/)[1]), olr = Number(lines[2].match(/OLR ([\d.]+)/)[1]);
  assert.ok(Math.abs(net - (asr - olr)) < 0.11, `the day's balance ${net} is the line's ASR − OLR ${asr - olr}`);
  r = run({ TAG: 'split', DAYS: '2' });
  assert.equal(r.status, 0, r.stderr.slice(-1500));
  r = run({ TAG: 'split', DAYS: '3' });
  assert.equal(r.status, 0, r.stderr.slice(-1500));
  const split = await decodeState(new Uint8Array(readFileSync(join(out, 'split_day0003.bin'))));
  assert.deepEqual([...split.energyRecord], [...whole.energyRecord]);
});
