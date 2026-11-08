import { test, before, after } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, mkdtempSync, mkdirSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { spawnSync } from 'node:child_process';
import { Grid } from '../js/grid.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { decodeState } from '../js/stateFile.module.js';
import { createIntervalMeans, meansGroup, meansEvery, MEANS } from '../js/intervalMeans.module.js';
import { withoutOptionGroups } from '../scripts/optionEnvironment.mjs';

/*
 * The 30-day means of js/intervalMeans.module.js on the GPU engine at
 * N=16: a sample reads what a frame reads and changes nothing a frame
 * would, the recorder's means are the means of the samples, and
 * scripts/spinup.mjs saves them, carries the ring across a split
 * inside a day and saves the same model without them. One GPU model at a
 * time (CLAUDE.md).
 */
let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};
const skip = !gpuAvailable && 'webgpu not installed';

const root = new URL('..', import.meta.url).pathname;
const topography = topographyFromInt16(readFileSync(join(root, 'data/topography_0p25.bin')).buffer);
const same = (a, b) => !!a && !!b && a.length === b.length && Array.prototype.every.call(a, (v, n) => Object.is(v, b[n]));
const N = 16, DT = 1350 * 16 / N, PER_DAY = 86400 / DT;
let dir;
before(() => { dir = mkdtempSync(join(tmpdir(), 'intervalMeansGpu-')); });
after(() => rmSync(dir, { recursive: true, force: true }));

async function fresh() {
  const model = await createGpuModel(new Grid(N), { topography });
  const { state, mesh } = model;
  initializeState(model, { geostrophic: !model.surfaceGeopotential }).forEach((values, a) => state[a].set(values));
  for (let i = 0; i < mesh.nCells; i++) if (model.geography.land[i]) state[6][i] = 0;
  model.load();
  model.ocean.initialize(state[3], state[6]);
  model.land.initialize();
  return model;
}

test('a sample reads every default field as the frame does, and a day\'s recorder means are the means of the day\'s samples', { skip }, async () => {
  const model = await fresh(), C = model.mesh.nCells, every = meansEvery(DT), recorder = createIntervalMeans(model, { fields: MEANS, every }), sum = new Float64Array(C);
  for (let n = 1; n <= PER_DAY; n++) {
    await model.step(DT);
    recorder.step();
    if (n % every === 0) (await model.sampleFrame({ fields: ['mslp'] })).fields.mslp.forEach((p, i) => { sum[i] += p; });
  }
  recorder.day();
  const { means } = meansGroup(await recorder.checkpoint());
  assert.deepEqual([means.samples, means.days], [PER_DAY / every, 1]);
  let worst = 0;
  for (let i = 0; i < C; i++) worst = Math.max(worst, Math.abs(means.mslp[i] - sum[i] / means.samples) / (sum[i] / means.samples));
  assert.ok(worst < 1e-6, `means.mslp lies up to ${worst} off the mean of the samples`);
  const sampled = await model.sampleFrame({ fields: MEANS }), framed = await model.beginFrame({ fields: MEANS });
  assert.deepEqual(MEANS.filter((name) => !same(sampled.fields[name], framed.fields[name])), []);
  const land = model.geography.land;
  for (const name of MEANS) assert.equal(Array.from(means[name]).filter((v, i) => !Number.isFinite(v) && !(name === 'sst' && land[i])).length, 0, name);
  assert.ok(means.precipitation.some((r) => r > 1) && means.zonalStress.some((t) => Math.abs(t) > 0.01) && Math.max(...means.shortwave) > 100, 'plausible rain, stress and sunlight');
  model.destroy();
});

test('a recorder leaves the GPU model and its day-mean diagnostics bit for bit as they step without one, a step at a time and in batches', { skip }, async () => {
  const run = async (record, batch) => {
    const model = await fresh(), recorder = record ? createIntervalMeans(model, { fields: MEANS, every: 1 }) : null;
    const steps = async (count) => {
      if (batch) for (let n = 0; n < count; n += 8) { await model.stepBatch(8, DT); recorder?.step(8); }
      else for (let n = 0; n < count; n++) { await model.step(DT); recorder?.step(); }
    };
    await steps(PER_DAY);
    recorder?.day();
    const d = await model.diagnostics();
    await steps(16);
    await model.sync();
    const ocean = await model.ocean.serialize({ restart: true }), land = await model.land.serialize();
    const fields = {
      ...Object.fromEntries(model.state.map((values, a) => [`state ${a}`, Float64Array.from(values)])),
      ...Object.fromEntries(Object.entries(ocean).filter(([, v]) => ArrayBuffer.isView(v) || Array.isArray(v)).map(([k, v]) => [`ocean.${k}`, Float64Array.from(v)])),
      ...Object.fromEntries(Object.entries(land).filter(([, v]) => ArrayBuffer.isView(v)).map(([k, v]) => [`land.${k}`, Float64Array.from(v)])),
      diagnostics: Float64Array.from(['absorbedSolar', 'outgoingLongwave', 'planetaryAlbedo', 'precipitation', 'meanSurfaceT', 'sensibleHeat', 'latentHeat', 'meanSteps'], (key) => d[key]),
    };
    const part = recorder ? await recorder.checkpoint() : null, samples = part ? part.days.reduce((sum, day) => sum + day.samples, 0) + part.daySamples : 0;
    model.destroy();
    return { fields, samples };
  };
  for (const batch of [false, true]) {
    const on = await run(true, batch), off = await run(false, batch);
    assert.deepEqual(Object.keys(off.fields).filter((name) => !same(on.fields[name], off.fields[name])), [], `${batch ? 'in batches' : 'a step at a time'}: fields apart`);
    assert.ok(Object.keys(off.fields).length > 20, Object.keys(off.fields).join(', '));
    assert.equal(on.samples, batch ? (PER_DAY + 16) / 8 : PER_DAY + 16);
  }
});

const runDir = (name) => { const out = join(dir, name); mkdirSync(out); return out; };
function spin(out, env) {
  const result = spawnSync(process.execPath, [join(root, 'scripts/spinup.mjs')], { cwd: root, env: { ...withoutOptionGroups(), N: String(N), TAG: 's', MINUTES: '100', OUT: out, ...env }, encoding: 'utf8' });
  if (result.status !== 0) console.log(`spinup.mjs exited ${result.status}: ${result.stderr.slice(-2000)}`);
  return result.status;
}
const load = async (file) => decodeState(new Uint8Array(readFileSync(file)));
const logOf = (out) => readFileSync(join(out, 's.log'), 'utf8');
const dayLine = (out) => logOf(out).split('\n').find((line) => line.startsWith('day 1 (')).replace(/^day 1 \([\d.]+ min\)/, 'day 1');
const arraysOf = (state) => [...Object.entries(state), ...Object.entries(state.ocean).map(([k, v]) => [`ocean.${k}`, v]), ...Object.entries(state.land).map(([k, v]) => [`land.${k}`, v])].filter(([, v]) => ArrayBuffer.isView(v));
const at = (state, name) => (name.startsWith('ocean.') ? state.ocean[name.slice(6)] : name.startsWith('land.') ? state.land[name.slice(5)] : state[name]);

test('a spin-up saves its 30-day means, the same across a split inside the day, and the same model and day line without them', { skip }, async () => {
  const whole = runDir('whole'), parts = runDir('parts'), off = runDir('off');
  assert.equal(spin(whole, { DAYS: '1' }), 0);
  const saved = await load(join(whole, 's_day0001.bin'));
  assert.deepEqual([saved.means.fields, saved.means.every, saved.means.samples, saved.means.days], [MEANS, 1, PER_DAY, 1]);
  assert.ok(!('daySamples' in saved.means) && !('mslpSum' in saved.means) && !('mslpRing' in saved.means));
  assert.match(logOf(whole), new RegExp(`means over the 30 days to each save: ${MEANS.join(', ')}, every 1 steps \\(64\\.0 a day\\)\\n`));
  assert.match(logOf(whole), /saved s_day0001\.bin after [\d.]+ min, with 1-day means ending at the save \(64 samples\)/);
  const table = spawnSync(process.execPath, [join(root, 'scripts/stateMeans.mjs'), join(whole, 's_day0001.bin'), 'mslp', 'sst'], { cwd: root, encoding: 'utf8' }).stdout;
  const row = (name) => (table.split('\n').find((line) => line.startsWith(`| ${name} |`)) ?? '').split('|').slice(2, 6).map((cell) => Number(cell.trim()));
  assert.match(table, /1-day means ending at the save \(model days 0\.000 to 1\.000\), 64 samples every 1 steps/);
  assert.ok(row('mslp')[0] > 99000 && row('mslp')[0] < 103000 && row('sst')[0] > 271 && row('sst')[0] < 305, table);

  assert.equal(spin(parts, { DAYS: '1', STOP_AFTER_STEPS: String(PER_DAY / 2) }), 0);
  const inside = await load(join(parts, `s_day0000_step${String(PER_DAY / 2).padStart(4, '0')}.bin`));
  assert.deepEqual([inside.means.days, inside.means.samples, inside.means.daySamples], [0, 0, PER_DAY / 2]);
  assert.ok(inside.means.mslpSum instanceof Float64Array && !('mslp' in inside.means));
  assert.equal(spin(parts, { DAYS: '1' }), 0);
  assert.match(logOf(parts), /carrying on the ring of s_day0000_step0032\.bin \(0 days, and 32 samples of the day under way\)/);
  const split = await load(join(parts, 's_day0001.bin'));
  assert.deepEqual([split.means.samples, split.means.days], [PER_DAY, 1]);
  assert.deepEqual(MEANS.filter((name) => !same(split.means[name], saved.means[name])), [], 'means apart across the split');

  assert.equal(spin(off, { DAYS: '1', MEANS: 'none' }), 0);
  const plain = await load(join(off, 's_day0001.bin'));
  assert.ok(!('means' in plain));
  assert.doesNotMatch(logOf(off), /means over the 30 days/);
  assert.deepEqual(arraysOf(saved).filter(([name, values]) => !same(values, at(plain, name))).map(([name]) => name), [], 'arrays the run without the means saved otherwise');
  assert.equal(dayLine(off), dayLine(whole));
});
