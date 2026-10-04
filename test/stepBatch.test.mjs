import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};
const { readBuffer } = gpuAvailable ? await import('../js/gpu/device.module.js') : {};
const { createForcingRecorder } = gpuAvailable ? await import('../js/gpu/forcing.gpu.js') : {};

const topography = syntheticTopography(90, 180, (lat, lon) => ((Math.cos(lon) > 0 && Math.abs(lat) < 1.2) || lat < -1.15 ? 300 + 2500 * Math.exp(-(((lat - 0.3) / 0.3) ** 2)) : -4000));
const DT = 900;

async function prepared(radiation = {}) {
  const model = await createGpuModel(new Grid(6), { topography, radiation, land: { growthTime: 3 * 3600, declineTime: 2 * 3600, snowDeclineTime: 4 * 3600 } });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i]) model.state[6][i] = 0;
  model.load();
  model.ocean.initialize(model.state[3], model.state[6]);
  model.land.initialize();
  const C = model.mesh.nCells, soil = new Float64Array(C), snow = new Float64Array(C);
  for (let i = 0; i < C; i++) if (model.geography.land[i]) { soil[i] = 40; snow[i] = model.mesh.latCell[i] > 1.0 ? 5 : 0; }
  model.land.load({ soil, snow, vegetation: new Float64Array(C).fill(0.5) });
  return model;
}

async function snapshot(model) {
  const { device, buffers } = model.gpu, ocean = model.oceanEngine.buffers;
  const named = { S: buffers.S, T: buffers.T, K1: buffers.K1, K2: buffers.K2, K3: buffers.K3, K4: buffers.K4, D: buffers.D, P: buffers.P, PH: buffers.PH, OS: ocean.S, OT: ocean.T, OK1: ocean.K1, OK2: ocean.K2, OK3: ocean.K3, OK4: ocean.K4, OD: ocean.OD };
  const out = { time: model.time, steps: model.gpu.stepCount, oceanCounter: model.oceanCounter, buffers: {} };
  for (const [name, buffer] of Object.entries(named)) out.buffers[name] = await readBuffer(device, buffer, buffer.size, Uint32Array);
  return out;
}

function assertIdentical(label, reference, got) {
  assert.equal(got.time, reference.time, `${label}: model time`);
  assert.equal(got.steps, reference.steps, `${label}: steps counted by the core`);
  assert.equal(got.oceanCounter, reference.oceanCounter, `${label}: the ocean's step counter`);
  for (const [name, bits] of Object.entries(reference.buffers)) {
    const other = got.buffers[name];
    assert.equal(other.length, bits.length, `${label}: ${name} length`);
    let differ = 0, first = -1;
    for (let x = 0; x < bits.length; x++) if (bits[x] !== other[x]) { differ++; if (first < 0) first = x; }
    const f = (a, x) => new Float32Array(a.buffer, a.byteOffset, a.length)[x];
    assert.equal(differ, 0, `${label}: ${name} differs in ${differ} words, first at ${first} (${first >= 0 ? f(bits, first) : ''} against ${first >= 0 ? f(other, first) : ''})`);
  }
}

test('stepBatch leaves the device exactly where the same steps taken one at a time do, across the ocean cadence', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const single = await prepared(), whole = await prepared(), split = await prepared();
  assert.equal(single.oceanEngine.everySteps, 4);
  for (let n = 0; n < 32; n++) await single.step(DT);
  await single.settle();
  const submissions = await whole.stepBatch(32, DT);
  assertIdentical('32 steps', await snapshot(single), await snapshot(whole));

  await single.step(DT);
  await single.settle();
  await whole.stepBatch(1, DT);
  await split.stepBatch(3, DT);
  await split.stepBatch(30, DT);
  const reference = await snapshot(single);
  assert.equal(reference.oceanCounter, 33);
  assertIdentical('32 + 1 steps', reference, await snapshot(whole));
  assertIdentical('3 + 30 steps', reference, await snapshot(split));

  const [a, b, c] = [await single.diagnostics(), await whole.diagnostics(), await split.diagnostics()];
  assert.deepEqual(b, a);
  assert.deepEqual(c, a);
  assert.ok(Number.isFinite(a.evaporation) && a.evaporation > 0 && Number.isFinite(a.oceanUpperDepth));
  console.log(`N=6 with land and ocean: 32 steps in ${submissions} submission(s), and 32 + 1 and 3 + 30 batched steps, match 33 single steps bit for bit in the atmosphere, physics and ocean buffers and in the diagnostics (Ts ${(a.meanSurfaceT - 273.15).toFixed(3)} °C, OLR ${a.outgoingLongwave.toFixed(3)} W/m², evaporation ${(86400 * a.evaporation).toFixed(4)} mm/d)`);
  for (const model of [single, whole, split]) model.destroy();
});

test('the forcing recorder writes the same day inside batches as after single steps', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const single = await prepared(), batched = await prepared();
  const recorders = [];
  for (const model of [single, batched]) { await model.diagnostics(); recorders.push(await createForcingRecorder(model)); }
  for (let n = 0; n < 33; n++) { await single.step(DT); recorders[0].step(); }
  await single.settle();
  const submissions = (await batched.stepBatch(5, DT, () => recorders[1].step())) + (await batched.stepBatch(28, DT, () => recorders[1].step()));
  assertIdentical('33 recorded steps', await snapshot(single), await snapshot(batched));
  const days = [];
  for (const [n, model] of [single, batched].entries()) { await model.diagnostics(); days.push(await recorders[n].day(1)); }
  assert.equal(days[1].length, days[0].length);
  let differ = 0;
  for (let x = 0; x < days[0].length; x++) if (days[0][x] !== days[1][x]) differ++;
  assert.equal(differ, 0, `the recorded day differs in ${differ} of ${days[0].length} bytes`);
  console.log(`33 recorded steps at N=6, batched as 5 + 28 in ${submissions} submissions: the ${days[0].length}-byte forcing day matches the single steps' byte for byte`);
  for (const model of [single, batched]) model.destroy();
});

test('with the radiation held between calls every third step, steps batched without waiting for the device leave it where single steps do', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const single = await prepared({ radiationEvery: 3 }), batched = await prepared({ radiationEvery: 3 });
  for (let n = 0; n < 14; n++) await single.step(DT);
  await single.settle();
  await batched.stepBatch(5, DT, null, false);
  await batched.stepBatch(9, DT, null, false);
  await batched.settle();
  assertIdentical('5 + 9 steps', await snapshot(single), await snapshot(batched));
  const [a, b] = [await single.diagnostics(), await batched.diagnostics()];
  assert.deepEqual(b, a);
  console.log(`N=6, radiation every 3 steps: 5 + 9 batched steps queued without waiting match 14 single steps bit for bit (OLR ${a.outgoingLongwave.toFixed(3)} W/m²)`);
  for (const model of [single, batched]) model.destroy();
});
