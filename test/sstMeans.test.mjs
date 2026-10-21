import { test, before, after } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, writeFileSync, mkdtempSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import { execFileSync, spawnSync } from 'node:child_process';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { gridHeader, stateGrid } from '../js/stateGrid.module.js';
import { encodeState, decodeState, groupHeader, groupArrays, f64Names } from '../js/stateFile.module.js';
import { decodeMonthlySst, climatologyAt, MONTHLY_SST_FILE } from '../js/ocean/climatology.module.js';
import { createSstRing, createSstRecorder, sstWindowOf, SST_WINDOW, SST_FIELDS } from '../js/sstMeans.module.js';
import { biasTable, markdown, windowMonths } from '../scripts/woaSstBias.mjs';

const root = new URL('..', import.meta.url).pathname;
const topography = topographyFromInt16(readFileSync(join(root, 'data/topography_0p25.bin')).buffer);
const woa = decodeMonthlySst(readFileSync(join(root, MONTHLY_SST_FILE)));
const plotting = spawnSync('python3', ['-c', 'import matplotlib, numpy'], { encoding: 'utf8' }).status === 0;
let dir;
before(() => { dir = mkdtempSync(join(tmpdir(), 'sstMeans-')); });
after(() => rmSync(dir, { recursive: true, force: true }));

function cpuModel(N = 8) {
  const model = createModel(new Grid(N), { topography });
  const { state, mesh } = model;
  initializeState(model, { geostrophic: !model.surfaceGeopotential }).forEach((values, a) => state[a].set(values));
  for (let i = 0; i < mesh.nCells; i++) if (model.geography.land[i]) state[6][i] = 0;
  model.ocean.initialize(state[3], state[6]);
  model.land.initialize();
  return model;
}
const same = (a, b) => a.length === b.length && Array.prototype.every.call(a, (v, n) => Object.is(v, b[n]));

test('the ring keeps its last days in order and starts again after a gap, and its window\'s mean is the mean of the daily means', () => {
  const C = 5, ring = createSstRing(C), daily = [];
  let seed = 7;
  const next = () => { seed = (seed * 16807) % 2147483647; return seed / 2147483647; };
  for (let day = 1; day <= 40; day++) {
    daily.push(Float32Array.from({ length: C }, () => 271 + 30 * next()));
    ring.push(day, daily[day - 1]);
  }
  assert.deepEqual(ring.days, Array.from({ length: SST_WINDOW }, (_, d) => 11 + d));
  const saved = { day: 40, ...ring.serialize() };
  assert.equal(saved.sstMeans.length, SST_WINDOW * C);
  const window = sstWindowOf(saved);
  assert.deepEqual(window.days, ring.days);
  for (let i = 0; i < C; i++) {
    let sum = 0;
    for (let d = 10; d < 40; d++) sum += daily[d][i];
    assert.ok(Math.abs(window.mean[i] - sum / SST_WINDOW) < 1e-9, `cell ${i}: ${window.mean[i]} against ${sum / SST_WINDOW}`);
  }
  assert.equal(sstWindowOf({ ...saved, day: 41 }), null, 'a ring that ends before the state\'s day');
  assert.equal(sstWindowOf({ ...saved, sstMeanDays: Float64Array.from(saved.sstMeanDays, (d, n) => (n < 5 ? d - 1 : d)) }), null, 'a ring with a gap');
  const short = createSstRing(C);
  for (let d = 1; d < SST_WINDOW; d++) short.push(d, daily[d]);
  assert.equal(sstWindowOf({ day: SST_WINDOW - 1, ...short.serialize() }), null, 'a ring short of the window');
  ring.push(42, daily[0]);
  assert.deepEqual(ring.days, [42]);
  const two = createSstRing(C, saved, 2);
  assert.deepEqual(two.days, [39, 40], 'a shorter ring keeps the saved one\'s last days');
});

test('a state keeps the ring and the part of a day bit for bit through its file, and a state without them loads without them', async () => {
  const C = 6, ring = createSstRing(C);
  for (let day = 698; day <= 730; day++) ring.push(day, Float32Array.from({ length: C }, (_, i) => 280 + 0.37 * i + 0.011 * day));
  const saved = ring.serialize(), sum = Float64Array.from({ length: C }, (_, i) => Math.PI * (i + 1) / 7);
  const base = { N: 1, K: 1, day: 730, time: 730 * 86400, surfaceT: new Float32Array(C) };
  const back = await decodeState(encodeState({ ...base, ...saved, sstSum: sum, sstSamples: 17 }));
  assert.ok(same(back.sstMeans, saved.sstMeans), 'the daily means');
  assert.ok(same(back.sstMeanDays, saved.sstMeanDays), 'their days');
  assert.ok(same(back.sstSum, sum), 'the part of a day keeps its float64 sum');
  assert.equal(back.sstSamples, 17);
  const again = createSstRing(C, back).serialize();
  assert.ok(same(again.sstMeans, saved.sstMeans) && same(again.sstMeanDays, saved.sstMeanDays), 'a ring loaded from the file saves the same');
  assert.deepEqual(sstWindowOf(back), sstWindowOf({ day: 730, ...saved }));
  const old = await decodeState(encodeState(base));
  assert.ok(!('sstMeans' in old) && !('sstMeanDays' in old) && !('sstSum' in old) && !('sstSamples' in old));
  assert.deepEqual(createSstRing(C, old).days, []);
  assert.equal(sstWindowOf(old), null);
});

test('the CPU recorder\'s day is the mean of the top layer after every ocean step, a checkpoint carries the day on, and the model steps bit for bit as without it', async () => {
  const a = cpuModel(), b = cpuModel(), C = a.mesh.nCells, dt = 2700, perDay = 32, { cellOcean } = a.ocean;
  const recorder = createSstRecorder(a), sum = new Float64Array(C);
  let taken = 0, second = null;
  for (let n = 1; n <= perDay; n++) {
    a.step(dt); recorder.step(); second?.step(); b.step(dt);
    if (n % a.ocean.everySteps === 0) { taken++; for (let i = 0; i < C; i++) sum[i] += a.ocean.T0[i]; }
    if (n === perDay / 2) second = createSstRecorder(a, await recorder.checkpoint());
  }
  const [whole, carried] = await Promise.all([recorder.day(), second.day()]);
  assert.equal(whole.samples, perDay / a.ocean.everySteps);
  assert.equal(carried.samples, whole.samples);
  assert.ok(same(whole.mean, carried.mean), 'a recorder carried on from a checkpoint ends the day with the same mean');
  let worst = 0;
  for (let i = 0; i < C; i++) {
    if (!cellOcean[i]) { assert.equal(whole.mean[i], 0); continue; }
    worst = Math.max(worst, Math.abs(whole.mean[i] - sum[i] / taken));
  }
  assert.ok(worst < 2e-5, `the recorder's day is ${worst} K off the mean of the samples`);
  for (let s = 0; s < a.state.length; s++) assert.ok(same(a.state[s], b.state[s]), `state ${s}`);
  const [oa, ob] = [a.ocean.serialize(), b.ocean.serialize()];
  for (const key of ['h', 'u', 'T', 'S', 'eta']) assert.ok(same(oa[key], ob[key]), `ocean.${key}`);
  const next = await recorder.day();
  assert.equal(next.samples, 0, 'the day ends the sum');
});

async function windowState() {
  const model = cpuModel(), C = model.mesh.nCells, { state } = model, ocean = model.ocean.serialize();
  const header = { N: 8, ...gridHeader(model.mesh), K: model.core.K, day: 730, time: 730 * 86400, terrain: true, levels: model.core.levels };
  const plain = join(dir, 'w8_day0730.bin'), windowed = join(dir, 'w8m_day0730.bin');
  writeFileSync(plain, encodeState({ ...header, surfaceT: state[3], ice: state[6], ocean }));
  const ring = createSstRing(C), days = Array.from({ length: SST_WINDOW }, (_, d) => 701 + d), offset = (day) => 0.1 * (day - 715);
  for (const day of days) ring.push(day, Float32Array.from({ length: C }, (_, i) => ocean.T[i] + offset(day)));
  writeFileSync(windowed, encodeState({ ...header, surfaceT: state[3], ice: state[6], ocean, ...ring.serialize() }));
  return { plain, windowed, days, shift: days.reduce((s, d) => s + offset(d), 0) / days.length, C };
}
let states = null;
const statesOnce = async () => (states ??= await windowState());

test('woaSstBias judges a state that carries a whole window by its mean against WOA23 averaged over the window\'s days, and one without, or under --instant, by its instant', async () => {
  const { plain, windowed, days, shift } = await statesOnce();
  const [mean, instant, forced] = await Promise.all([biasTable([{ path: windowed }], woa, topography), biasTable([{ path: plain }], woa, topography), biasTable([{ path: windowed }], woa, topography, { instant: true })]);
  const dated = await biasTable(days.map((day) => ({ path: plain, day: day - 0.5 })), woa, topography);
  assert.deepEqual(mean.dates.map((d) => [d.kind, d.days, d.length]), [['mean', [701, 730], SST_WINDOW]]);
  assert.deepEqual([instant.dates[0].kind, instant.dates[0].day, forced.dates[0].kind], ['instant', 730, 'instant']);
  let boxes = 0;
  for (const [n, row] of mean.rows.entries()) {
    if (!Number.isFinite(row.woa)) continue;
    boxes++;
    const step = row.cells === null ? 0 : shift;
    assert.ok(Math.abs(row.model - instant.rows[n].model - step) < 1e-4, `${row.box}: the mean is ${row.model - instant.rows[n].model} K over the instant, not ${step}`);
    assert.ok(Math.abs(row.woa - dated.rows[n].woa) < 1e-9, `${row.box}: WOA23 over the window ${row.woa} against the mean over its days' middles ${dated.rows[n].woa}`);
    assert.equal(forced.rows[n].model, instant.rows[n].model);
  }
  assert.ok(boxes >= 6, `${boxes} boxes with water at N=8`);
  const [line, , note] = markdown(mean).split('\n');
  assert.match(line, /^- w8m_day0730\.bin: model days 701–730 \(Feb 18 – Mar 20\), the 30-day mean SST against WOA23 over the same days, /);
  assert.match(note, /mean SST averages its daily means/);
  assert.match(markdown(instant).split('\n')[0], /^- w8_day0730\.bin: model day 730, Mar 20 \(day of year 78\.0\), the state's instant SST against WOA23 /);
});

test('the ocean figure plots the window\'s mean SST against WOA23 over the same days and names the window, or the instant when the state carries none', async () => {
  const { plain, windowed, days, C } = await statesOnce();
  const dump = (state, name) => {
    const out = join(dir, name);
    execFileSync(process.execPath, [join(root, 'scripts/figures/oceanPanels.mjs'), state, out], { cwd: root, encoding: 'utf8' });
    return { out, d: JSON.parse(readFileSync(out, 'utf8')) };
  };
  const mean = dump(windowed, 'mean.json'), instant = dump(plain, 'instant.json');
  assert.deepEqual(mean.d.sstWindow, { kind: 'mean', days: [701, 730], length: SST_WINDOW });
  assert.deepEqual(instant.d.sstWindow, { kind: 'instant', day: 730 });
  assert.ok(!('sstMean' in instant.d));
  const months = windowMonths(woa, days), saved = await decodeState(new Uint8Array(readFileSync(windowed)));
  const window = sstWindowOf(saved), { mesh } = stateGrid(saved, topography);
  let sea = 0;
  for (let i = 0; i < C; i++) {
    if (mean.d.land[i]) continue;
    sea++;
    assert.equal(mean.d.sst[i], instant.d.sst[i]);
    assert.ok(Math.abs(mean.d.sstMean[i] - (window.mean[i] - 273.15)) <= 0.0051, `cell ${i}: ${mean.d.sstMean[i]}`);
    assert.ok(Math.abs(mean.d.woa[i] - (climatologyAt(woa, months, mesh.latCell[i], mesh.lonCell[i]) - 273.15)) <= 0.0051, `cell ${i}: WOA23 ${mean.d.woa[i]}`);
  }
  assert.ok(sea > 100, `${sea} sea cells`);
  if (!plotting) return;
  const printed = (json, png) => execFileSync('python3', [join(root, 'scripts/figures/oceanPanels.py'), json, png, 'w8'], { cwd: join(root, 'scripts/figures'), encoding: 'utf8' });
  assert.match(printed(mean.out, join(dir, 'mean.png')), /^ocean: 30-day mean SST of days 701–730 − WOA 60S–60N/m);
  assert.match(printed(instant.out, join(dir, 'instant.png')), /^ocean: instant SST − WOA 60S–60N/m);
});

test('the page\'s worker starts from a state that carries an SST ring as from the same state without it, and its snapshots carry no ring', async () => {
  const messages = [], fetched = globalThis.fetch;
  globalThis.self = { postMessage: (message) => messages.push(message) };
  globalThis.fetch = async (url, init) => {
    if (!String(url).startsWith('file:')) return fetched(url, init);
    try { return new Response(readFileSync(new URL(String(url)))); } catch { return new Response(null, { status: 404 }); }
  };
  try {
    await import('../js/model.worker.js');
    const send = async (data, type) => {
      messages.length = 0;
      await self.onmessage({ data });
      for (let wait = 0; wait < 6000 && !messages.some((m) => m.type === 'frame'); wait++) await new Promise((resolve) => setTimeout(resolve, 10));
      const failed = messages.find((m) => m.type === 'status' && String(m.text).startsWith('error'));
      assert.equal(failed, undefined, failed && failed.text);
      const found = messages.find((m) => m.type === type);
      assert.ok(found, `no ${type} message`);
      return found;
    };
    const start = (from = null) => send({ type: 'start', N: 8, engine: 'cpu', paused: true, ...(from ? { from: pathToFileURL(from).href } : {}) }, 'ready');
    const snapshot = () => send({ type: 'snapshot' }, 'snapshotData');
    await start();
    const fresh = await snapshot();
    const saved = { N: fresh.N, K: fresh.K, day: 5, time: 5 * 86400, terrain: fresh.terrain, ...(fresh.grid ? { grid: fresh.grid } : {}) };
    for (const [key, buffer] of Object.entries(fresh.arrays)) saved[key] = new Float64Array(buffer);
    saved.ocean = { ...groupHeader(fresh.ocean), ...Object.fromEntries(Object.entries(groupArrays(fresh.ocean)).map(([key, buffer]) => [key, new Float64Array(buffer)])) };
    saved.land = Object.fromEntries(Object.entries(fresh.land).map(([key, buffer]) => [key, new Float64Array(buffer)]));
    const C = saved.surfaceT.length, plain = join(dir, 'plain_day0005.bin'), ringed = join(dir, 'ringed_day0005.bin');
    writeFileSync(plain, encodeState(saved, { f64: f64Names(saved) }));
    writeFileSync(ringed, encodeState({ ...saved, sstMeans: new Float32Array(2 * C).fill(290), sstMeanDays: Float64Array.of(4, 5), sstSum: new Float64Array(C).fill(0.5), sstSamples: 3 }, { f64: f64Names(saved) }));
    const fromPlain = { ready: await start(plain), snapshot: await snapshot() };
    const fromRinged = { ready: await start(ringed), snapshot: await snapshot() };
    const names = ({ arrays, ocean, land, ...rest }) => [...Object.keys(rest), ...Object.keys(arrays), ...Object.keys(ocean ?? {}), ...Object.keys(land ?? {})];
    for (const name of SST_FIELDS) assert.ok(!names(fromRinged.snapshot).includes(name), name);
    assert.deepEqual(names(fromRinged.snapshot), names(fromPlain.snapshot));
    assert.equal(fromPlain.ready.day, 5);
    assert.equal(fromRinged.ready.day, 5);
    for (const group of ['arrays', 'ocean', 'land']) {
      for (const [key, buffer] of Object.entries(groupArrays(fromPlain.snapshot[group]))) assert.ok(Buffer.from(buffer).equals(Buffer.from(fromRinged.snapshot[group][key])), `${group}.${key}`);
    }
  } finally {
    globalThis.fetch = fetched;
    delete globalThis.self;
  }
});
