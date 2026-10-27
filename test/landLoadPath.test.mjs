import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createLandSurface } from '../js/physics/land.module.js';
import { loadSavedLand } from '../scripts/savedLand.mjs';
import { savedFixture, cpuStepped, digest, landArrays, options } from './helpers/savedLand.mjs';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

// pinned: the fixture's land after its load, and a CPU model one 900 s step after it; the snow on the sea ice is in both
const LAND_AFTER_LOAD = 'edf240830ae56066578fd689da9ecc94759cf7b9d6da46018ffa60ae5989b33d';
const CPU_ONE_STEP = 'c930fecc36dc27f81d621a6f49029dc61e84ef83b0ffcd67ceefd5aa6e5b88ff';

async function gpuLoaded(fixture) {
  const model = await createGpuModel(new Grid(16), options);
  fixture.state.forEach((a, k) => model.state[k].set(a));
  model.load({ ocean: false });
  return model;
}

test('a slab land loaded the way the scripts load it, from a save with snow on sea ice, is bit-identical to the load of 5738db9, and a CPU model steps from it to the same digest', () => {
  const f = savedFixture();
  assert.ok(f.iced.length > 10);
  const before = createLandSurface(f.mesh, f.geography, {}), now = createLandSurface(f.mesh, f.geography, { soil: 'slab' });
  before.load(f.saved.land, f.state[6]);
  now.load(f.saved.land, f.state[6], f.saved.time, f.state[3]);
  assert.equal(digest(landArrays(before.serialize())), LAND_AFTER_LOAD);
  assert.equal(digest(landArrays(now.serialize())), LAND_AFTER_LOAD);
  for (const i of f.iced) assert.equal(now.snow[i], 12, `sea ice cell ${i}`);
  assert.equal(cpuStepped(f, (model, saved) => model.land.load(saved.land, model.state[6])), CPU_ONE_STEP);
  assert.equal(cpuStepped(f, (model, saved) => model.land.load(saved.land, model.state[6], saved.time, model.state[3])), CPU_ONE_STEP);
  const dropped = createLandSurface(f.mesh, f.geography, {});
  dropped.load(f.saved.land, null, f.saved.time, f.state[3]);
  for (const i of f.iced) assert.equal(dropped.snow[i], 0);
  assert.notEqual(digest(landArrays(dropped.serialize())), LAND_AFTER_LOAD);
});

test('the scripts\' load keeps the snow a save carries on sea ice on the GPU model, where a null ice argument drops it', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const f = savedFixture(), model = await gpuLoaded(f);
  loadSavedLand(model, f.saved);
  for (const i of f.iced) assert.equal(model.land.snow[i], 12, `sea ice cell ${i}`);
  const open = Array.from({ length: f.mesh.nCells }, (_, i) => i).filter((i) => !f.geography.land[i] && !(f.state[6][i] > 0));
  assert.ok(open.length > 10);
  for (const i of open) assert.equal(model.land.snow[i], 0, `open sea cell ${i}`);
  await model.sync();
  const kept = await model.land.serialize();
  for (const i of f.iced) assert.equal(kept.snow[i], 12, `sea ice cell ${i} on the device`);
  model.land.load(f.saved.land, null, f.saved.time, model.state[3]);
  for (const i of f.iced) assert.equal(model.land.snow[i], 0);
});

test('a slab GPU model loaded the way fixedSst loads it steps once to the same state, bit for bit, as one loaded with the call 5738db9 made', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const f = savedFixture(), results = [];
  for (const load of [(model) => model.land.load(f.saved.land), (model) => loadSavedLand(model, f.saved)]) {
    const model = await gpuLoaded(f);
    load(model);
    const land = ['soil', 'surface', 'snow', 'vegetation', 'snowAlbedo', 'canopy', 'soilCarbon'].map((name) => Float64Array.from(model.land[name]));
    await model.step(900);
    await model.sync();
    results.push({ loaded: digest(land), stepped: digest([...model.state, ...landArrays(await model.land.serialize())]) });
  }
  assert.equal(results[1].loaded, results[0].loaded);
  assert.equal(results[1].stepped, results[0].stepped);
});
