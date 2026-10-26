import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, writeFileSync, mkdtempSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { CUMULUS_FIELDS, placeCumulus } from '../js/physics/moist.module.js';
import { encodeState, groupHeader, groupArrays, f64Names } from '../js/stateFile.module.js';

const N = 16;

const random = (length, seed) => {
  let s = seed;
  return Float64Array.from({ length }, () => { s = (Math.imul(s, 1103515245) + 12345) & 0x7fffffff; return 0.5 + s / 2147483648; });
};
const same = (a, b) => a.length === b.length && a.every((x, i) => x === b[i]);
const zero = (a) => a.every((x) => x === 0);

let shared = null;
function cpuModel() {
  if (!shared) {
    const model = createModel(new Grid(N), { ocean: false });
    shared = { model, C: model.mesh.nCells, K: model.core.K, K0: model.moist.cumulusK0 };
  }
  for (const name of CUMULUS_FIELDS) shared.model.moist[name].fill(9);
  return shared;
}

test(`placeCumulus puts a GPU-layout cumulus cloud in layers K0 to K − 1 of a CPU model's K × C arrays with zeros above, and copies a K × C cloud unchanged, at N=${N}`, () => {
  const { model, C, K, K0 } = cpuModel();
  assert.ok(K0 > 0 && K0 < K, `the plume layers are a tail of the column: K0 = ${K0} of K = ${K}`);
  const tail = { cumulusCover: random((K - K0) * C, 1), cumulusWater: random((K - K0) * C, 2) };
  placeCumulus(model, tail);
  for (const name of CUMULUS_FIELDS) {
    const target = model.moist[name];
    assert.equal(target.length, K * C);
    assert.ok(zero(target.subarray(0, K0 * C)), `${name}: zero above layer ${K0}`);
    assert.ok(same(target.subarray(K0 * C), tail[name]), `${name}: layers ${K0} to ${K - 1} are the saved plume layers`);
  }
  const whole = { cumulusCover: random(K * C, 3), cumulusWater: random(K * C, 4) };
  placeCumulus(model, whole);
  for (const name of CUMULUS_FIELDS) assert.ok(same(model.moist[name], whole[name]), `${name}: a K × C cloud is copied unchanged`);
});

test(`placeCumulus leaves a CPU model's cloud zero for no state, a state without it and a field saved in neither layout, and places each field on its own, at N=${N}`, () => {
  const { model, C, K, K0 } = cpuModel();
  const plume = (K - K0) * C;
  const cases = [['no state', null], ['a state without the fields', { pi: random(C, 5) }], ['one layer', { cumulusCover: random(C, 6), cumulusWater: random(C, 7) }], ['one value short', { cumulusCover: random(plume - 1, 8), cumulusWater: random(plume - 1, 9) }], ['one value long', { cumulusCover: random(plume + 1, 10), cumulusWater: random(K * C + 1, 11) }]];
  for (const [what, saved] of cases) {
    for (const name of CUMULUS_FIELDS) model.moist[name].fill(9);
    placeCumulus(model, saved);
    for (const name of CUMULUS_FIELDS) assert.ok(zero(model.moist[name]), `${name} is zero for ${what}`);
  }
  const cover = random(plume, 12);
  for (const name of CUMULUS_FIELDS) model.moist[name].fill(9);
  placeCumulus(model, { cumulusCover: cover });
  assert.ok(same(model.moist.cumulusCover.subarray(K0 * C), cover), 'the cover saved alone is placed');
  assert.ok(zero(model.moist.cumulusWater), 'the water not saved is zero');
});

test(`placeCumulus fills the GPU model's (K − K0) × C arrays from either layout, and no layout fills them zero, at N=${N}`, () => {
  const { C, K, K0 } = cpuModel();
  const plume = (K - K0) * C;
  const gpu = { mesh: { nCells: C }, core: { K }, moist: { cumulusK0: K0, cumulusCover: new Float64Array(plume).fill(9), cumulusWater: new Float64Array(plume).fill(9) } };
  const tail = { cumulusCover: random(plume, 13), cumulusWater: random(plume, 14) };
  placeCumulus(gpu, tail);
  for (const name of CUMULUS_FIELDS) assert.ok(same(gpu.moist[name], tail[name]), `${name}: a GPU-layout cloud is copied unchanged`);
  const whole = { cumulusCover: random(K * C, 15), cumulusWater: random(K * C, 16) };
  placeCumulus(gpu, whole);
  for (const name of CUMULUS_FIELDS) assert.ok(same(gpu.moist[name], whole[name].subarray(K0 * C)), `${name}: a K × C cloud gives its plume layers`);
  placeCumulus(gpu, { cumulusCover: random(plume - 1, 17) });
  for (const name of CUMULUS_FIELDS) assert.ok(zero(gpu.moist[name]), `${name}: zero for a field saved in neither layout`);
  placeCumulus(gpu, { cumulusCover: Float32Array.from(tail.cumulusCover) });
  assert.ok(same(gpu.moist.cumulusCover, Float64Array.from(Float32Array.from(tail.cumulusCover))), 'a float32 field is widened');
});

test("the page's worker places a state's cumulus cloud saved in either layout and saves it in the GPU's", async () => {
  const dir = mkdtempSync(join(tmpdir(), 'cumulus-loader-'));
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
      for (let wait = 0; wait < 12000 && !messages.some((m) => m.type === 'frame'); wait++) await new Promise((resolve) => setTimeout(resolve, 10));
      const failed = messages.find((m) => m.type === 'status' && String(m.text).startsWith('error'));
      assert.equal(failed, undefined, failed && failed.text);
      const found = messages.find((m) => m.type === type);
      assert.ok(found, `no ${type} message`);
      return found;
    };
    const start = (from = null) => send({ type: 'start', N, engine: 'cpu', paused: true, ...(from ? { from: pathToFileURL(from).href } : {}) }, 'ready');
    const snapshot = () => send({ type: 'snapshot' }, 'snapshotData');
    await start();
    const fresh = await snapshot();
    const C = fresh.arrays.surfaceT.byteLength / 8, K = fresh.K, K0 = K - fresh.arrays.cumulusCover.byteLength / 8 / C;
    assert.ok(K0 > 0 && K0 < K, `the saved plume layers are a tail of the column: K0 = ${K0} of K = ${K}`);
    for (const name of CUMULUS_FIELDS) assert.ok(zero(new Float64Array(fresh.arrays[name])), `${name} starts zero`);
    const saved = { N: fresh.N, K, day: fresh.day, time: fresh.time, terrain: fresh.terrain, ...(fresh.grid ? { grid: fresh.grid } : {}) };
    for (const [key, buffer] of Object.entries(fresh.arrays)) saved[key] = new Float64Array(buffer.slice(0));
    saved.ocean = { ...groupHeader(fresh.ocean), ...Object.fromEntries(Object.entries(groupArrays(fresh.ocean)).map(([key, buffer]) => [key, new Float64Array(buffer.slice(0))])) };
    saved.land = Object.fromEntries(Object.entries(fresh.land).map(([key, buffer]) => [key, new Float64Array(buffer.slice(0))]));
    const tail = { cumulusCover: random((K - K0) * C, 18), cumulusWater: random((K - K0) * C, 19) }, whole = { cumulusCover: random(K * C, 20), cumulusWater: random(K * C, 21) };
    const carried = async (cloud, tag) => {
      const file = join(dir, `${tag}_day0000.bin`), withCloud = { ...saved, ...cloud };
      writeFileSync(file, encodeState(withCloud, { f64: f64Names(withCloud) }));
      await start(file);
      return (await snapshot()).arrays;
    };
    const fromTail = await carried(tail, 'tail');
    for (const name of CUMULUS_FIELDS) assert.ok(same(new Float64Array(fromTail[name]), tail[name]), `${name}: a GPU-layout cloud comes back as saved`);
    const fromWhole = await carried(whole, 'whole');
    for (const name of CUMULUS_FIELDS) assert.ok(same(new Float64Array(fromWhole[name]), whole[name].subarray(K0 * C)), `${name}: a K × C cloud comes back as its plume layers`);
    await start();
    const again = await snapshot();
    for (const name of CUMULUS_FIELDS) assert.ok(zero(new Float64Array(again.arrays[name])), `${name}: a start without a state is zero again`);
  } finally {
    globalThis.fetch = fetched;
    delete globalThis.self;
    rmSync(dir, { recursive: true, force: true });
  }
});
