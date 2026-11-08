import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { encodeState, decodeState, groupHeader, groupArrays, f64Names } from '../js/stateFile.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

test("on the GPU, a page save waits for the step that ends an ocean step, a held radiation interval and a held drag interval, and a model that has run 3 steps, restored from it through encodeState, decodeState and the page's restore, steps on frame by frame as the run that did not stop over two ocean cadences and a step, with runoff and discharge pending at the save; and so again from the step-8 save of a run that loaded that save relabeled to step 6, off the ocean's cadence", { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const stand = !globalThis.self;
  if (stand) globalThis.self = { postMessage() {} };
  try {
    const { snapshotArrays, placeSaved, stepsToCadence } = await import('../js/model.worker.js');
    const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 ? 1500 + 1000 * Math.cos(lat) : -4000));
    const options = { topography, radiation: { radiationEvery: 4 }, dragEvery: 2, ocean: { everySteps: 4 } };
    const dt = 1350, names = ['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'], subscription = { level: 'surface', depth: 'surface', fields: [], diagnostics: false };
    const frame = async (model, steps) => { const capturing = model.beginFrame(subscription); await model.stepBatch(steps, dt, null, false); await capturing; };
    const run = await createGpuModel(new Grid(16), options);
    const start = { N: 16, K: run.core.K, day: 0, time: 0, terrain: true }, init = initializeState(run, {});
    names.forEach((name, a) => { start[name] = Float64Array.from(init[a]); });
    await placeSaved(run, start, 16, dt);
    assert.deepEqual(run.cadence, { radiationEvery: 4, dragEvery: 2, oceanEvery: 4 });
    const water = async (model, change) => { const held = await model.readFreshwater(); for (const name of ['runoff', 'discharge']) change(held[`${name}Total`], held[`${name}Seen`], name); model.writeFreshwater(held); };
    await frame(run, 3);
    const wait = stepsToCadence(run, dt);
    assert.equal(wait, 1, 'a save asked for after 3 steps waits for the fourth');
    await frame(run, wait);
    await water(run, (total, seen, name) => { total.fill(name === 'runoff' ? 3 : 4); seen.fill(1); });
    const save = async (model) => {
      const { arrays, ocean, land } = await snapshotArrays(model);
      const saved = { N: 16, K: model.core.K, day: model.time / 86400, time: model.time, terrain: true };
      for (const [key, buffer] of Object.entries(arrays)) saved[key] = new Float64Array(buffer.slice(0));
      saved.ocean = { ...groupHeader(ocean), ...Object.fromEntries(Object.entries(groupArrays(ocean)).map(([key, buffer]) => [key, new Float64Array(buffer.slice(0))])) };
      if (land) saved.land = Object.fromEntries(Object.entries(land).map(([key, buffer]) => [key, new Float64Array(buffer.slice(0))]));
      return decodeState(encodeState(saved, { f64: f64Names(saved) }));
    };
    const back = await save(run);
    for (const key of ['windSpeed', 'evaporation', 'exchangeWind', 'rainTotal', 'runoffTotal', 'rainSeen', 'runoffSeen']) assert.ok(back[key] && back[key].length > 0, `the snapshot carries ${key}`);
    for (const key of ['runoffTotal', 'dischargeTotal', 'runoffSeen', 'dischargeSeen']) assert.ok(Math.max(...back[key]) > 0, `${key} is not 0 at the save`);
    for (const key of ['previousT0', 'fresh', 'flux', 'boundaryDepths']) assert.ok(back.ocean[key] && back.ocean[key].length > 0, `the snapshot carries the ocean's restart array ${key}`);
    const twin = await createGpuModel(new Grid(16), options);
    await placeSaved(twin, start, 16, dt);
    await frame(twin, 3);
    await placeSaved(twin, back, 16, dt);
    const largest = (mine, theirs) => { let d = 0, s = 0; for (let i = 0; i < theirs.length; i++) { d = Math.max(d, Math.abs(mine[i] - theirs[i])); s = Math.max(s, Math.abs(theirs[i])); } return s > 0 ? d / s : d; };
    const follow = async (continued, restored, step) => {
      const steps = 2 * continued.cadence.oceanEvery + 1, worst = names.map(() => 0), parted = [], oceanWorst = { h: 0, u: 0, T: 0, S: 0, eta: 0 };
      for (let s = 0; s < steps; s++) {
        await frame(continued, 1);
        await frame(restored, 1);
        await continued.sync();
        await restored.sync();
        names.forEach((name, a) => { const gap = largest(restored.state[a], continued.state[a]); worst[a] = Math.max(worst[a], gap); if (gap > 0 && !parted.some(([field]) => field === name)) parted.push([name, s + 1]); });
        const [mine, theirs] = await Promise.all([restored.ocean.serialize(), continued.ocean.serialize()]);
        for (const key of Object.keys(oceanWorst)) { const gap = largest(mine[key], theirs[key]); oceanWorst[key] = Math.max(oceanWorst[key], gap); if (gap > 0 && !parted.some(([field]) => field === `ocean ${key}`)) parted.push([`ocean ${key}`, s + 1]); }
      }
      console.log(`over ${steps} steps after a restore at step ${step}, the largest gap from the continued run over the field's largest value: ${names.map((name, a) => `${name} ${worst[a]}`).join(', ')}; first parting steps: ${parted.map(([name, s]) => `${name} ${s}`).join(', ') || 'none'}; the ocean over the steps: ${Object.entries(oceanWorst).map(([key, gap]) => `${key} ${gap}`).join(', ')}`);
      names.forEach((name, a) => assert.equal(worst[a], 0, `step ${step}, ${name}: the restored run parts from the continued run by ${worst[a]}`));
      for (const [key, gap] of Object.entries(oceanWorst)) assert.equal(gap, 0, `step ${step}: the ocean's ${key} parts from the continued run's by ${gap}`);
    };
    await follow(run, twin, 4);
    const off = await createGpuModel(new Grid(16), options);
    await placeSaved(off, { ...back, time: 6 * dt, day: (6 * dt) / 86400 }, 16, dt);
    const ahead = stepsToCadence(off, dt);
    assert.equal(ahead, 2, 'loaded at step 6, a save waits for step 8');
    await frame(off, ahead);
    await water(off, (total, seen, name) => { total.fill(name === 'runoff' ? 5 : 6); seen.fill(2); });
    await placeSaved(twin, await save(off), 16, dt);
    await follow(off, twin, 8);
  } finally {
    if (stand) delete globalThis.self;
  }
});
