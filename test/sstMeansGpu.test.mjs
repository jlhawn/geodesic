import { test, before, after } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, writeFileSync, mkdtempSync, mkdirSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { spawnSync } from 'node:child_process';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { decodeState, encodeState } from '../js/stateFile.module.js';
import { createSstRecorder } from '../js/sstMeans.module.js';

/*
 * The SST recorders of js/sstMeans.module.js and js/gpu/sstMeans.gpu.js
 * against each other and the engines' own top layers at N=8, the GPU
 * model with and without one, and the daily means scripts/spinup.mjs
 * keeps in its states at N=6. One GPU model at a time (CLAUDE.md).
 */
let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};
const { createGpuSstRecorder } = gpuAvailable ? await import('../js/gpu/sstMeans.gpu.js') : {};
const { readRanges } = gpuAvailable ? await import('../js/gpu/device.module.js') : {};
const skip = !gpuAvailable && 'webgpu not installed';

const root = new URL('..', import.meta.url).pathname;
const topography = topographyFromInt16(readFileSync(join(root, 'data/topography_0p25.bin')).buffer);
const same = (a, b) => !!a && !!b && a.length === b.length && Array.prototype.every.call(a, (v, n) => Object.is(v, b[n]));
const DT = 2700, PER_DAY = 32;
let dir;
before(() => { dir = mkdtempSync(join(tmpdir(), 'sstMeansGpu-')); });
after(() => rmSync(dir, { recursive: true, force: true }));

async function fresh(engine) {
  const model = engine === 'gpu' ? await createGpuModel(new Grid(8), { topography }) : createModel(new Grid(8), { topography });
  const { state, mesh } = model;
  initializeState(model, { geostrophic: !model.surfaceGeopotential }).forEach((values, a) => state[a].set(values));
  for (let i = 0; i < mesh.nCells; i++) if (model.geography.land[i]) state[6][i] = 0;
  if (engine === 'gpu') model.load();
  model.ocean.initialize(state[3], state[6]);
  model.land.initialize();
  return model;
}

test('the CPU and GPU recorders agree at N=8: each day\'s mean is its engine\'s top layer averaged over the day\'s ocean steps, and the engines\' means part by no more than their top layers do', { skip }, async () => {
  const days = 2;
  const run = async (engine) => {
    const model = await fresh(engine), C = model.mesh.nCells, gpu = engine === 'gpu', ocean = gpu ? model.oceanEngine : model.ocean;
    const recorder = gpu ? createGpuSstRecorder(model) : createSstRecorder(model), samples = [], means = [];
    for (let n = 1; n <= days * PER_DAY; n++) {
      await model.step(DT);
      recorder.step();
      if (n % ocean.everySteps === 0) samples.push(gpu ? Float64Array.from((await readRanges(model.gpu.device, ocean.buffers.OD, [{ offset: ocean.layout.OD.T0, length: C }]))[0]) : Float64Array.from(ocean.T0));
      if (n % PER_DAY === 0) means.push(await recorder.day());
    }
    const out = { C, samples, means, cells: Uint8Array.from(ocean.cellOcean), perDay: PER_DAY / ocean.everySteps };
    if (gpu) { recorder.destroy(); model.destroy(); }
    return out;
  };
  const cpu = await run('cpu'), gpu = await run('gpu'), { C, cells, perDay } = cpu;
  assert.ok(same(cells, gpu.cells) && perDay === gpu.perDay);
  const lines = [];
  for (let d = 0; d < days; d++) {
    const envelope = new Float64Array(C), own = { cpu: 0, gpu: 0 };
    let largest = 0, squares = 0, sea = 0;
    for (const [name, r] of Object.entries({ cpu, gpu })) {
      assert.equal(r.means[d].samples, perDay, `${name} day ${d + 1}`);
      for (let i = 0; i < C; i++) {
        if (!cells[i]) { assert.equal(r.means[d].mean[i], 0); continue; }
        let sum = 0;
        for (let s = d * perDay; s < (d + 1) * perDay; s++) sum += r.samples[s][i];
        own[name] = Math.max(own[name], Math.abs(r.means[d].mean[i] - sum / perDay));
      }
    }
    for (let i = 0; i < C; i++) {
      if (!cells[i]) continue;
      for (let s = d * perDay; s < (d + 1) * perDay; s++) envelope[i] = Math.max(envelope[i], Math.abs(cpu.samples[s][i] - gpu.samples[s][i]));
      const apart = Math.abs(cpu.means[d].mean[i] - gpu.means[d].mean[i]);
      assert.ok(apart <= envelope[i] + 1.2e-4, `day ${d + 1} cell ${i}: the means part by ${apart} K, the top layers by at most ${envelope[i]}`);
      largest = Math.max(largest, apart); squares += apart * apart; sea++;
    }
    assert.ok(own.cpu < 2e-5 && own.gpu < 1e-4, `day ${d + 1}: the recorders' means lie ${own.cpu} K (CPU) and ${own.gpu} K (GPU) off their samples' means`);
    lines.push(`day ${d + 1}: CPU − GPU mean rms ${Math.sqrt(squares / sea).toExponential(2)} K, largest ${largest.toExponential(2)} K over ${sea} sea cells; each recorder within ${own.cpu.toExponential(1)} / ${own.gpu.toExponential(1)} K of its samples' mean`);
    assert.ok(Math.sqrt(squares / sea) < 0.02, lines[d]);
  }
  console.log(lines.join('\n'));
});

test('a GPU recorder leaves the model bit for bit as it steps without one, a step at a time and in batches, across the day\'s end', { skip }, async () => {
  const run = async (record, batch) => {
    const model = await fresh('gpu'), recorder = record ? createGpuSstRecorder(model) : null;
    const steps = async (count) => {
      if (batch) for (let n = 0; n < count; n += 8) await model.stepBatch(8, DT, recorder ? () => recorder.step() : null);
      else for (let n = 0; n < count; n++) { await model.step(DT); recorder?.step(); }
    };
    await steps(PER_DAY);
    const day = recorder ? recorder.day() : null;
    await steps(16);
    await model.sync();
    const ocean = await model.ocean.serialize({ restart: true }), land = await model.land.serialize();
    const fields = {
      ...Object.fromEntries(model.state.map((values, a) => [`state ${a}`, Float64Array.from(values)])),
      ...Object.fromEntries(Object.entries(ocean).filter(([, v]) => ArrayBuffer.isView(v) || Array.isArray(v)).map(([k, v]) => [`ocean.${k}`, Float64Array.from(v)])),
      ...Object.fromEntries(Object.entries(land).filter(([, v]) => ArrayBuffer.isView(v)).map(([k, v]) => [`land.${k}`, Float64Array.from(v)])),
    };
    const samples = day ? (await day).samples : null;
    recorder?.destroy();
    model.destroy();
    return { fields, samples };
  };
  for (const batch of [false, true]) {
    const on = await run(true, batch), off = await run(false, batch);
    const apart = Object.keys(off.fields).filter((name) => !same(on.fields[name], off.fields[name]));
    assert.deepEqual(apart, [], `${batch ? 'in batches' : 'a step at a time'}: fields apart`);
    assert.ok(Object.keys(off.fields).length > 20, Object.keys(off.fields).join(', '));
    assert.equal(on.samples, PER_DAY / 4);
  }
});

const N = 6, STEPS = 24, OCEAN_STEPS = STEPS / 4;
const runDir = (name) => { const out = join(dir, name); mkdirSync(out); return out; };
function spin(out, env) {
  const result = spawnSync(process.execPath, [join(root, 'scripts/spinup.mjs')], { cwd: root, env: { ...process.env, N: String(N), TAG: 's', MINUTES: '100', OUT: out, ...env }, encoding: 'utf8' });
  if (result.status !== 0) console.log(`spinup.mjs exited ${result.status}: ${result.stderr.slice(-2000)}`);
  return result.status;
}
const load = async (file) => decodeState(new Uint8Array(readFileSync(file)));
const logOf = (out) => readFileSync(join(out, 's.log'), 'utf8');
const SST_KEYS = ['sstMeans', 'sstMeanDays', 'sstSum', 'sstSamples'];
async function withoutMeans(from, to) {
  const saved = await load(from);
  for (const key of SST_KEYS) delete saved[key];
  writeFileSync(to, encodeState(saved));
}
const arraysOf = (state) => [...Object.entries(state), ...Object.entries(state.ocean).map(([k, v]) => [`ocean.${k}`, v]), ...Object.entries(state.land).map(([k, v]) => [`land.${k}`, v])].filter(([, v]) => ArrayBuffer.isView(v));
const at = (state, name) => (name.startsWith('ocean.') ? state.ocean[name.slice(6)] : name.startsWith('land.') ? state.land[name.slice(5)] : state[name]);

let whole = null;
test('a spin-up keeps the daily means of its last SST_WINDOW days in its states, the same across a split at a day\'s end and inside a day, and saves the same model without them', { skip }, async () => {
  const out = runDir('whole'), parts = runDir('parts'), off = runDir('off');
  assert.equal(spin(out, { DAYS: '3', SST_WINDOW: '2' }), 0);
  whole = await load(join(out, 's_day0003.bin'));
  const C = whole.surfaceT.length;
  assert.deepEqual(Array.from(whole.sstMeanDays), [2, 3]);
  assert.match(logOf(out), new RegExp(`SST means: the daily mean of the ocean's top layer over its ${OCEAN_STEPS} steps a day, the last 2 days kept \\(none yet\\)`));
  assert.doesNotMatch(logOf(out), /left out/);
  let sea = 0, moved = 0;
  for (let i = 0; i < C; i++) {
    const [day2, day3] = [whole.sstMeans[i], whole.sstMeans[C + i]];
    if (day3 === 0) continue;
    sea++;
    assert.ok(day2 > 270 && day2 < 310 && day3 > 270 && day3 < 310, `cell ${i}: ${day2} and ${day3} K`);
    if (day3 !== day2) moved++;
  }
  assert.ok(sea > 150 && moved > sea / 2, `${sea} sea cells, ${moved} changed from day 2 to 3`);

  assert.equal(spin(parts, { DAYS: '1', SST_WINDOW: '2' }), 0);
  assert.deepEqual(Array.from((await load(join(parts, 's_day0001.bin'))).sstMeanDays), [1]);
  await withoutMeans(join(parts, 's_day0001.bin'), join(runDir('oldDay'), 's_day0001.bin'));
  assert.equal(spin(parts, { DAYS: '3', SST_WINDOW: '2', STOP_AFTER_STEPS: String(STEPS / 2) }), 0);
  const inside = await load(join(parts, 's_day0001_step0012.bin'));
  assert.deepEqual([inside.sstSamples, inside.sstSum.length, Array.from(inside.sstMeanDays)], [OCEAN_STEPS / 2, C, [1]]);
  await withoutMeans(join(parts, 's_day0001_step0012.bin'), join(runDir('oldInside'), 's_day0001_step0012.bin'));
  assert.equal(spin(parts, { DAYS: '3', SST_WINDOW: '2' }), 0);
  const split = await load(join(parts, 's_day0003.bin')), dayMeans = ['convectiveRain', 'largeScaleRain', 'meanAbsorbedSolar', 'meanOutgoingLongwave', 'meanPlanetaryAlbedo', 'meanShortwaveCloudEffect', 'meanLongwaveCloudEffect', 'energyRecord'];
  assert.deepEqual(arraysOf(whole).filter(([name, values]) => !dayMeans.includes(name) && !same(values, at(split, name))).map(([name]) => name), [], 'arrays the split run saved otherwise');
  assert.ok(same(split.sstMeans, whole.sstMeans) && same(split.sstMeanDays, whole.sstMeanDays));

  assert.equal(spin(off, { DAYS: '3', SST_WINDOW: '0' }), 0);
  const plain = await load(join(off, 's_day0003.bin'));
  assert.ok(SST_KEYS.every((key) => !(key in plain)));
  assert.doesNotMatch(logOf(off), /SST means/);
  const apart = arraysOf(whole).filter(([name]) => !SST_KEYS.includes(name)).filter(([name, values]) => !same(values, at(plain, name))).map(([name]) => name);
  assert.deepEqual(apart, [], 'arrays the run without the means saved otherwise');
  assert.deepEqual(arraysOf(plain).map(([name]) => name).sort(), arraysOf(whole).map(([name]) => name).filter((name) => !SST_KEYS.includes(name)).sort());
});

test('a run continued from a state saved without the means starts them afresh, and leaves out a day it stepped only in part', { skip }, async () => {
  assert.ok(whole, 'the run above saved its day 3');
  const C = whole.surfaceT.length, oldDay = join(dir, 'oldDay'), oldInside = join(dir, 'oldInside');
  assert.equal(spin(oldDay, { DAYS: '3', SST_WINDOW: '2' }), 0);
  const fromDay = await load(join(oldDay, 's_day0003.bin'));
  assert.match(logOf(oldDay), /the last 2 days kept \(none yet\)/);
  assert.deepEqual(Array.from(fromDay.sstMeanDays), [2, 3]);
  assert.ok(same(fromDay.sstMeans, whole.sstMeans), 'the means of days 2 and 3 are the whole run\'s');
  assert.equal(spin(oldInside, { DAYS: '3', SST_WINDOW: '2' }), 0);
  const fromInside = await load(join(oldInside, 's_day0003.bin')), text = logOf(oldInside);
  assert.match(text, new RegExp(`day 2's covers only its last ${STEPS / 2} steps and stays out`));
  assert.match(text, new RegExp(`SST mean of day 2 over ${OCEAN_STEPS / 2} of its ${OCEAN_STEPS} ocean steps left out; the window starts again`));
  assert.deepEqual(Array.from(fromInside.sstMeanDays), [3]);
  assert.ok(same(fromInside.sstMeans, whole.sstMeans.subarray(C)), 'day 3\'s mean is the whole run\'s');
});
