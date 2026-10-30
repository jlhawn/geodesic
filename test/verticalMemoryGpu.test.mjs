import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { layerVerticalVelocity, layerFieldAtLevel, smoothCells } from '../js/levels.module.js';
import { memoryWeight, VERTICAL_MEMORY } from '../js/frames.module.js';
import { encodeState, decodeState } from '../js/stateFile.module.js';
import { meanVertical } from '../js/means.module.js';

/*
 * The GPU model's memory of its layers' vertical velocity against the CPU
 * model's at N=16: the fold of the device's own state, the weights, the
 * cadence, the frame that reads the memory at a level, a save and a load
 * that continue the mean, and the physics unseen by any of it. One GPU model
 * at a time (CLAUDE.md).
 */
let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const N = 16, DT = 900, OPTIONS = { ocean: false };
const skip = !gpuAvailable && 'webgpu not installed';

async function build(options = {}) {
  const model = await createGpuModel(new Grid(N), { ...OPTIONS, ...options });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  model.load();
  return model;
}
const velocities = (model) => layerVerticalVelocity(model.mesh, model.core, model.state[0], model.state[2], model.state[1], new Float64Array(model.core.K * model.mesh.nCells));
const largest = (values) => values.reduce((m, v) => Math.max(m, Math.abs(v)), 0);
const worst = (a, b) => a.reduce((m, v, n) => Math.max(m, Math.abs(v - b[n])), 0);
const same = (a, b, label) => {
  assert.equal(a.length, b.length, `${label}: lengths`);
  for (let n = 0; n < a.length; n++) if (!Object.is(a[n], b[n])) assert.fail(`${label}: ${a[n]} against ${b[n]} at ${n}`);
};

// The reference mean of the states a model passes through, folded as the engines do.
function follower(model, every = 1) {
  const reference = new Float64Array(model.core.K * model.mesh.nCells);
  let folds = 0, last = NaN, calls = 0;
  return {
    reference,
    get folds() { return folds; },
    get last() { return last; },
    fold() {
      calls++;
      if (folds > 0 && calls % every !== 0) return;
      const w = velocities(model), alpha = memoryWeight(folds > 0 ? model.time - last : 0, folds, VERTICAL_MEMORY);
      for (let m = 0; m < reference.length; m++) reference[m] = alpha === 1 ? w[m] : reference[m] + alpha * (w[m] - reference[m]);
      folds++; last = model.time;
    },
  };
}

test('stepped, the device\'s memory is the weighted fold of the velocities the CPU derives from its own state, at every cadence', { skip }, async () => {
  for (const every of [1, 3]) {
    const model = await build({ verticalEvery: every }), mean = follower(model, every);
    try {
      assert.equal(model.verticalMemory.every, every);
      for (let n = 1; n <= 8; n++) {
        await model.step(DT);
        await model.sync();
        mean.fold();
        const { values } = model.verticalMemory, scale = largest(mean.reference), gap = worst(mean.reference, values);
        if (n === 1 || n === 8) console.log(`fold every ${every}, step ${n}: largest |w| ${scale.toExponential(2)} m/s, GPU parts from the CPU fold by ${gap.toExponential(2)}`);
        assert.equal(model.verticalMemory.folds, mean.folds, `folds after step ${n}`);
        assert.equal(model.verticalMemory.time, mean.last, `time of the last fold after step ${n}`);
        assert.ok(scale > 1e-4, `the layers move at up to ${scale} m/s`);
        assert.ok(gap <= 2e-4 * scale, `fold every ${every}, step ${n}: the GPU memory parts from the CPU's fold by ${gap} against ${scale}`);
      }
    } finally { model.destroy(); }
  }
});

test('the frame reads the memory at a level, interpolated in ln p between the layers\' middles and averaged with the neighbours', { skip }, async () => {
  const model = await build();
  try {
    for (let n = 0; n < 6; n++) await model.step(DT);
    await model.sync();
    const { mesh, core, state } = model, C = mesh.nCells, { K, sigmaMid } = core.diagnostics, memory = Float64Array.from(model.verticalMemory.values);
    for (const level of ['surface', 850, 500, 250, 50]) {
      const frame = await model.beginFrame({ level, fields: ['vertical'] });
      const expected = smoothCells(mesh, layerFieldAtLevel(sigmaMid, K, C, state[0], memory, level, new Float64Array(C)), new Float64Array(C));
      const gap = worst(expected, frame.fields.vertical), scale = largest(expected);
      assert.ok(scale > 1e-4 && gap <= 1e-4 * scale, `at ${level}: ${gap} against ${scale}`);
    }
    const held = await model.beginFrame({ level: 500, fields: ['vertical'] });
    same((await model.beginFrame({ level: 500, fields: ['vertical'] })).fields.vertical, held.fields.vertical, 'with no step between, the frame shows the same memory');
    await model.step(DT);
    const later = await model.beginFrame({ level: 500, fields: ['vertical'] });
    assert.ok(worst(held.fields.vertical, later.fields.vertical) > 0, 'a step folds the memory again');
    await model.beginFrame({ level: 250, fields: ['vertical'] });
    await model.sync();
    assert.equal(model.verticalMemory.folds, 7, 'the frames fold nothing themselves once the model has stepped');
  } finally { model.destroy(); }
});

test('a frame that finds no memory folds the state once, so a field shows at once', { skip }, async () => {
  const model = await build();
  try {
    assert.equal(model.verticalMemory.folds, 0);
    const frame = await model.beginFrame({ level: 500, fields: ['vertical'] });
    await model.sync();
    assert.equal(model.verticalMemory.folds, 1);
    const { mesh, core, state } = model, C = mesh.nCells;
    const expected = smoothCells(mesh, layerFieldAtLevel(core.diagnostics.sigmaMid, core.diagnostics.K, C, state[0], velocities(model), 500, new Float64Array(C)), new Float64Array(C));
    const gap = worst(expected, frame.fields.vertical), scale = largest(expected);
    assert.ok(scale > 1e-4 && gap <= 2e-4 * scale, `the state\'s own field: ${gap} against ${scale}`);
  } finally { model.destroy(); }
});

test('the prognostic fields are the same bit for bit whatever the memory\'s cadence', { skip }, async () => {
  const states = [];
  for (const every of [1, 4]) {
    const model = await build({ verticalEvery: every });
    try {
      for (let n = 0; n < 6; n++) await model.step(DT);
      await model.sync();
      states.push(model.state.map((a) => Float64Array.from(a)));
    } finally { model.destroy(); }
  }
  states[0].forEach((a, s) => same(a, states[1][s], `state array ${s}`));
});

test('a state saved from the device carries the memory and its clock, and a model loaded from it continues the mean', { skip }, async () => {
  const source = await build();
  let saved;
  try {
    for (let n = 0; n < 5; n++) await source.step(DT);
    await source.sync();
    const carried = source.verticalMemory.serialize();
    assert.deepEqual(Array.from(carried.verticalClock), [5, 5 * DT]);
    const header = { N, K: source.core.K, pi: source.state[0], levels: source.core.levels, ...carried };
    saved = { arrays: source.state.map((a) => Float64Array.from(a)), time: source.time, decoded: await decodeState(encodeState(header)), memory: Float64Array.from(source.verticalMemory.values) };
    same(saved.decoded.verticalMemory, Float32Array.from(carried.verticalMemory), 'a float32 state holds the device\'s own single-precision memory');
    assert.equal(saved.decoded.verticalClock[1], 5 * DT);
    const view = meanVertical(saved.decoded, 500), expected = layerFieldAtLevel(source.core.diagnostics.sigmaMid, source.core.K, source.mesh.nCells, source.state[0], saved.decoded.verticalMemory, 500);
    assert.ok(worst(expected, view) <= 1e-6 * largest(expected), 'the figures\' helper reads the device\'s memory');
  } finally { source.destroy(); }

  const next = [];
  for (const withMemory of [true, false]) {
    const model = await createGpuModel(new Grid(N), OPTIONS);
    try {
      saved.arrays.forEach((a, s) => model.state[s].set(a));
      model.time = saved.time;
      assert.equal(model.verticalMemory.load(withMemory ? saved.decoded : {}), withMemory);
      model.load();
      await model.sync();
      if (withMemory) {
        same(await model.gpu.vertical.read(), Float32Array.from(saved.decoded.verticalMemory), 'the load restores the memory on the device');
        assert.equal(model.gpu.vertical.folds, 5);
        assert.equal(model.gpu.vertical.time, 5 * DT);
        const frame = await model.beginFrame({ level: 500, fields: ['vertical'] });
        const shown = smoothCells(model.mesh, layerFieldAtLevel(model.core.diagnostics.sigmaMid, model.core.K, model.mesh.nCells, model.state[0], saved.memory, 500, new Float64Array(model.mesh.nCells)), new Float64Array(model.mesh.nCells));
        assert.ok(worst(shown, frame.fields.vertical) <= 1e-5 * largest(shown), 'a loaded state shows its settled field at once');
      } else {
        assert.equal(model.gpu.vertical.folds, 0);
        assert.ok((await model.gpu.vertical.read()).every((v) => v === 0), 'a load without a memory leaves none on the device');
      }
      await model.step(DT);
      await model.sync();
      const w = velocities(model), folds = withMemory ? 5 : 0, previous = withMemory ? saved.memory : new Float64Array(w.length), alpha = memoryWeight(withMemory ? DT : 0, folds, VERTICAL_MEMORY);
      const expected = Float64Array.from(w, (x, m) => (alpha === 1 ? x : previous[m] + alpha * (x - previous[m])));
      const gap = worst(expected, model.verticalMemory.values), scale = largest(expected);
      console.log(`after a ${withMemory ? 'load with the memory' : 'load without it'} and a step: GPU parts from the CPU fold by ${gap.toExponential(2)} against ${scale.toExponential(2)}`);
      assert.ok(gap <= 2e-4 * scale, `${withMemory ? 'with' : 'without'} the memory: ${gap} against ${scale}`);
      assert.equal(model.verticalMemory.folds, folds + 1);
      next.push(model.state.map((a) => Float64Array.from(a)));
    } finally { model.destroy(); }
  }
  next[0].forEach((a, s) => same(a, next[1][s], `after a load with the memory and without it, state array ${s}`));
});

test('the two engines\' memories agree after steps from the same state', { skip }, async () => {
  const cpu = createModel(new Grid(N), OPTIONS), gpu = await build();
  try {
    const init = initializeState(cpu, {});
    for (let a = 0; a < init.length; a++) cpu.state[a].set(init[a]);
    for (let n = 0; n < 4; n++) { cpu.step(DT); await gpu.step(DT); }
    await gpu.sync();
    const a = cpu.verticalMemory.values, b = gpu.verticalMemory.values, scale = largest(a), gap = worst(a, b);
    console.log(`after 4 steps the CPU memory reaches ${scale.toExponential(2)} m/s and the GPU's parts from it by ${gap.toExponential(2)}`);
    assert.equal(gpu.verticalMemory.folds, cpu.verticalMemory.folds);
    assert.equal(gpu.verticalMemory.time, cpu.verticalMemory.time);
    assert.ok(scale > 1e-4 && gap <= 0.05 * scale, `the engines' memories part by ${gap} against ${scale}`);
  } finally { gpu.destroy(); }
});
