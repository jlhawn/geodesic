import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { join } from 'node:path';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { decodeState, encodeState, f64Names } from '../js/stateFile.module.js';
import { FIELDS } from '../js/frames.module.js';
import { createIntervalMeans, meansGroup, savedMeansPart, meansFields, meansEvery, availableMeans, withoutMeans, MEANS, MEANS_WINDOW, MEMORY_FIELDS, DERIVED_MEANS } from '../js/intervalMeans.module.js';

const root = new URL('..', import.meta.url).pathname;
const same = (a, b) => !!a && !!b && a.length === b.length && Array.prototype.every.call(a, (v, n) => Object.is(v, b[n]));
const N = 16, DT = 1350 * 16 / N, PER_DAY = 86400 / DT;

function counting(C) {
  const model = { time: 0, mesh: { nCells: C }, reads: 0 };
  model.sampleFrame = ({ fields }) => {
    const k = ++model.reads;
    return new Promise((resolve) => setTimeout(() => resolve({ fields: Object.fromEntries(fields.map((name, f) => [name, Float32Array.from({ length: C }, (_, i) => (i === 2 && f === 1 && k === 2 ? NaN : (f + 1) * k + i))])) }), 3 - (k % 3)));
  };
  return model;
}

test('a recorder samples every `every` steps, sums in order however its reads land, keeps the last `depth` daily means and carries them and the day under way on from its checkpoint', async () => {
  const model = counting(4), recorder = createIntervalMeans(model, { fields: ['a', 'b'], every: 3, depth: 2 });
  for (let n = 1; n <= 27; n++) { recorder.step(); if (n % 6 === 0) recorder.day(); }
  const part = await recorder.checkpoint();
  assert.deepEqual([part.days.map((d) => d.samples), part.daySamples, part.every, part.fields, Array.from(part.sums.a)], [[2, 2], 1, 3, ['a', 'b'], [9, 10, 11, 12]]);
  const { means } = meansGroup(part);
  assert.deepEqual([means.fields, means.every, means.days, means.samples, Array.from(means.a), Array.from(means.b)], [['a', 'b'], 3, 2, 4, [6.5, 7.5, 8.5, 9.5], [13, 14, 15, 16]]);
  assert.ok(!('daySamples' in means) && !('aSum' in means) && !('aRing' in means));

  const carried = createIntervalMeans(model, { fields: ['a', 'b'], every: 3, depth: 2, offset: 3, from: part });
  for (let n = 0; n < 3; n++) carried.step();
  carried.day();
  const on = await carried.checkpoint();
  assert.deepEqual([on.days.length, on.daySamples, Array.from(meansGroup(on).means.a)], [2, 0, [8.5, 9.5, 10.5, 11.5]]);

  const deep = createIntervalMeans(counting(4), { fields: ['a', 'b'], every: 3 });
  for (let n = 1; n <= 12; n++) { deep.step(); if (n % 6 === 0) deep.day(); }
  const wide = meansGroup(await deep.checkpoint()).means;
  assert.ok(wide.days === 2 && Number.isNaN(wide.b[2]) && wide.b[3] === 8, 'a day with a sample that is not a number is not one in the window');

  const late = createIntervalMeans(counting(4), { fields: ['a'], every: 1, offset: 2 });
  for (let n = 0; n < 4; n++) late.step();
  late.day();
  assert.deepEqual(meansGroup(await late.checkpoint()), {}, 'a day begun partway through without a checkpoint stays out');

  const batched = createIntervalMeans(model, { fields: ['a'], every: 3 });
  batched.step(8); batched.step(1); batched.step(8);
  assert.equal((await batched.checkpoint()).daySamples, 3);
  assert.throws(() => createIntervalMeans(model, { fields: ['a'], every: 3, from: part }), /checkpoint samples a, b every 3/);
});

test('the means group round-trips through a binary state: float32 means and, saved for a run to carry on, the float32 ring and the float64 sums of the day under way', async () => {
  const model = counting(5), recorder = createIntervalMeans(model, { fields: ['mslp', 'sst'], every: 2 });
  for (let n = 1; n <= 10; n++) { recorder.step(); if (n % 4 === 0) recorder.day(); }
  const part = await recorder.checkpoint(), state = { N: 1, day: 0, time: 10 * 1350, pi: Float64Array.of(1, 2, 3, 4, 5) };
  const whole = await decodeState(encodeState({ ...state, ...meansGroup(part) }));
  assert.deepEqual([whole.means.fields, whole.means.every, whole.means.days, whole.means.samples], [['mslp', 'sst'], 2, 2, 4]);
  assert.ok(whole.means.mslp instanceof Float32Array && same(whole.means.mslp, meansGroup(part).means.mslp) && Number.isNaN(whole.means.sst[2]));
  assert.equal(savedMeansPart(whole), null);
  const inside = await decodeState(encodeState({ ...state, ...meansGroup(part, { inside: true }) }));
  assert.ok(inside.means.mslpSum instanceof Float64Array && inside.means.mslpRing instanceof Float32Array && f64Names(inside).includes('means.sstSum') && !f64Names(inside).includes('means.sstRing'));
  assert.deepEqual(savedMeansPart(inside), part);
  assert.deepEqual(withoutMeans(inside), await decodeState(encodeState(state)));
});

test('MEANS names the default list, none, or frame fields other than the memories, the default samples 64 times a day over 30 days, and a model offers the fields its frames give', () => {
  assert.deepEqual(meansFields(undefined), MEANS);
  assert.deepEqual(meansFields('none'), []);
  assert.deepEqual(meansFields('["mslp","sst","mslp"]'), ['mslp', 'sst']);
  for (const bad of ['["rain"]', '["nothing"]', '{"mslp":1}']) assert.throws(() => meansFields(bad), /MEANS is 'none' or a JSON list/);
  assert.ok(MEANS.every((name) => (name in FIELDS || name in DERIVED_MEANS) && !MEMORY_FIELDS.includes(name)));
  assert.deepEqual([16, 64, 128].map((n) => meansEvery(1350 * 16 / n)), [1, 4, 8]);
  assert.equal(MEANS_WINDOW, 30);
  const offering = { time: 0, frameFields: new Set(['mslp', 'sst']), sampleFrame() {} };
  assert.deepEqual([availableMeans(offering, ['mslp', 'zonalStress', 'sst']), availableMeans({}, ['zonalStress'])], [['mslp', 'sst'], ['zonalStress']]);
  assert.throws(() => createIntervalMeans(offering, { fields: ['mslp', 'zonalStress'] }), /frames give no zonalStress/);
  assert.deepEqual(availableMeans({ frameFields: new Set(['surfaceT', 'sst']) }, ['gridBoxSkin', 'sst']), ['sst']);
});

test('gridBoxSkin is the mean of each sample\'s c·surfaceT + (1 − c)·sst over sea cells with ice and of surfaceT elsewhere, and a recorder hands back the last day it closed', async () => {
  const frames = [
    { surfaceT: [250, 260, 300, 290], sst: [271, 272, NaN, 290], concentration: [0.5, 0, 1, 0] },
    { surfaceT: [254, 262, 302, 292], sst: [271, 273, NaN, 291], concentration: [1, 0.25, 0, 0] },
  ];
  let k = 0;
  const model = { time: 0, mesh: { nCells: 4 }, frameFields: new Set(['surfaceT', 'sst', 'concentration']), asked: null };
  model.sampleFrame = ({ fields }) => { model.asked = fields; const f = frames[k++]; return Promise.resolve({ fields: Object.fromEntries(fields.map((name) => [name, Float32Array.from(f[name])])) }); };
  const recorder = createIntervalMeans(model, { fields: ['gridBoxSkin', 'surfaceT'] });
  assert.equal(await recorder.latest(), null);
  recorder.step(); recorder.step(); recorder.day();
  assert.deepEqual(new Set(model.asked), new Set(['surfaceT', 'sst', 'concentration']));
  const day = await recorder.latest(), { means } = meansGroup(await recorder.checkpoint());
  const expected = [(0.5 * 250 + 0.5 * 271 + 254) / 2, (260 + 0.25 * 262 + 0.75 * 273) / 2, (300 + 302) / 2, (290 + 292) / 2];
  for (let i = 0; i < 4; i++) assert.ok(Math.abs(means.gridBoxSkin[i] - expected[i]) < 1e-4 && Math.abs(day.values.gridBoxSkin[i] - expected[i]) < 1e-4, `cell ${i}: ${means.gridBoxSkin[i]} against ${expected[i]}`);
  assert.deepEqual(Array.from(means.surfaceT), [252, 261, 301, 291]);
  assert.deepEqual(meansFields('["gridBoxSkin"]'), ['gridBoxSkin']);
});

const topography = topographyFromInt16(readFileSync(join(root, 'data/topography_0p25.bin')).buffer);
function fresh() {
  const model = createModel(new Grid(N), { topography });
  const { state, mesh } = model;
  initializeState(model, { geostrophic: !model.surfaceGeopotential }).forEach((values, a) => state[a].set(values));
  for (let i = 0; i < mesh.nCells; i++) if (model.geography.land[i]) state[6][i] = 0;
  model.ocean.initialize(state[3], state[6]);
  model.land.initialize();
  return model;
}
const mslpOf = (model) => {
  const { core, state: [pi, theta] } = model, C = model.mesh.nCells, phis = model.surfaceGeopotential, { R, g, exnerLayer } = core.diagnostics, bottom = (core.K - 1) * C;
  return Float64Array.from(pi, (p, i) => p * Math.exp(phis[i] / (R * (theta[bottom + i] * exnerLayer[bottom + i] + 0.00325 * phis[i] / g))));
};

test('a day of the CPU model at N=16: the day\'s mean of mslp is the mean of its per-step sea-level pressure, every default field is a number where it should be, and the means reload from a saved state', async () => {
  const model = fresh(), C = model.mesh.nCells, every = meansEvery(DT), recorder = createIntervalMeans(model, { fields: MEANS, every, step: DT }), sum = new Float64Array(C);
  for (let n = 1; n <= PER_DAY; n++) {
    model.step(DT);
    recorder.step();
    if (n % every === 0) mslpOf(model).forEach((p, i) => { sum[i] += p; });
  }
  recorder.day();
  const part = await recorder.checkpoint(), { means } = meansGroup(part);
  assert.deepEqual([means.samples, means.every, means.days], [PER_DAY / every, every, 1]);
  let worst = 0;
  for (let i = 0; i < C; i++) worst = Math.max(worst, Math.abs(means.mslp[i] - sum[i] / means.samples) / (sum[i] / means.samples));
  assert.ok(worst < 1e-6, `means.mslp lies up to ${worst} off the mean of the per-step sea-level pressure`);
  const land = model.geography.land;
  for (const name of MEANS) {
    const bad = Array.from(means[name]).filter((v, i) => !Number.isFinite(v) && !(name === 'sst' && land[i])).length;
    assert.equal(bad, 0, `${name}: ${bad} cells not a number`);
  }
  assert.ok(Array.from(means.sst).every((v, i) => Number.isNaN(v) === !!land[i]), 'sst is NaN over land alone');
  const rain = Array.from(means.precipitation), seaT = means.sst.filter(Number.isFinite);
  assert.ok(rain.every((r) => r >= 0) && rain.some((r) => r > 1), `precipitation up to ${Math.max(...rain)} mm/d`);
  assert.ok(Math.min(...seaT) > 265 && Math.max(...seaT) < 310 && Math.max(...means.shortwave) > 100 && means.zonalStress.some((t) => Math.abs(t) > 0.01), 'plausible sst, sunlight and stress');
  const state = { N, day: 1, time: model.time, pi: model.state[0] }, reloaded = await decodeState(encodeState({ ...state, ...meansGroup(part) }));
  assert.deepEqual(reloaded.means.fields, MEANS);
  for (const name of MEANS) assert.ok(same(reloaded.means[name], means[name]), name);
});

test('the CPU model steps bit for bit the same with a recorder sampling every default field every step as without one', async () => {
  const run = async (record) => {
    const model = fresh(), recorder = record ? createIntervalMeans(model, { fields: MEANS, every: 1, step: DT }) : null;
    for (let n = 1; n <= PER_DAY + 16; n++) { model.step(DT); recorder?.step(); if (n % PER_DAY === 0) recorder?.day(); }
    const part = recorder ? await recorder.checkpoint() : null, samples = part ? part.days.reduce((sum, day) => sum + day.samples, 0) + part.daySamples : 0;
    const ocean = await model.ocean.serialize({ restart: true }), land = await model.land.serialize();
    return {
      samples,
      fields: {
        ...Object.fromEntries(model.state.map((values, a) => [`state ${a}`, Float64Array.from(values)])),
        ...Object.fromEntries(Object.entries(ocean).filter(([, v]) => ArrayBuffer.isView(v)).map(([k, v]) => [`ocean.${k}`, Float64Array.from(v)])),
        ...Object.fromEntries(Object.entries(land).filter(([, v]) => ArrayBuffer.isView(v)).map(([k, v]) => [`land.${k}`, Float64Array.from(v)])),
      },
    };
  };
  const on = await run(true), off = await run(false);
  assert.equal(on.samples, PER_DAY + 16);
  assert.ok(Object.keys(off.fields).length > 20, Object.keys(off.fields).join(', '));
  assert.deepEqual(Object.keys(off.fields).filter((name) => !same(on.fields[name], off.fields[name])), []);
});
