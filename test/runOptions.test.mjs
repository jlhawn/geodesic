import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, writeFileSync, mkdtempSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import { encodeState, decodeState, groupHeader, groupArrays, f64Names } from '../js/stateFile.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { syntheticTopography } from '../js/geography.module.js';

const SEED = 'soilWater_seed.bin sha256:0123abcd';
const OPTIONS = {
  radiation: { cloudColumns: true, cloudSpread: 0.75, varianceCover: 'water', buoyancyLength: 0.53 },
  moist: { stratocumulusInversion: 7, cumulusEntrainmentHeight: 1, cumulusAscent: 'kinetic' },
  surface: { relativeWind: true },
};

test("the page's worker runs a saved state with the options it records, and its snapshots and restores carry them and the soil seed the state descends from", async () => {
  const dir = mkdtempSync(join(tmpdir(), 'run-options-'));
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
    const start = (from = null) => send({ type: 'start', N: 16, engine: 'cpu', paused: true, ...(from ? { from: pathToFileURL(from).href } : {}) }, 'ready');
    const snapshot = () => send({ type: 'snapshot' }, 'snapshotData');
    const restore = ({ type, ...taken }) => send({ type: 'restore', snapshot: taken }, 'ready');
    assert.deepEqual((await start()).options, {});
    const fresh = await snapshot();
    assert.deepEqual(fresh.options, {});
    assert.equal(fresh.soilSeed, null, 'a fresh start descends from no soil seed');
    assert.equal(fresh.arrays.seaFlow, undefined);
    const saved = { N: fresh.N, K: fresh.K, day: fresh.day, time: fresh.time, terrain: fresh.terrain, ...(fresh.grid ? { grid: fresh.grid } : {}), options: OPTIONS, soilSeed: SEED };
    for (const [key, buffer] of Object.entries(fresh.arrays)) saved[key] = new Float64Array(buffer.slice(0));
    saved.ocean = { ...groupHeader(fresh.ocean), ...Object.fromEntries(Object.entries(groupArrays(fresh.ocean)).map(([key, buffer]) => [key, new Float64Array(buffer.slice(0))])) };
    saved.land = Object.fromEntries(Object.entries(fresh.land).map(([key, buffer]) => [key, new Float64Array(buffer.slice(0))]));
    saved.seaFlow = Float64Array.from({ length: saved.u.length / fresh.K }, (_, e) => Math.fround(0.25 * Math.sin(0.37 * e)));
    const file = join(dir, 'optioned_day0000.bin');
    writeFileSync(file, encodeState(saved, { f64: f64Names(saved) }));
    assert.deepEqual((await start(file)).options, OPTIONS);
    const taken = await snapshot();
    assert.deepEqual(taken.options, OPTIONS);
    assert.equal(taken.soilSeed, SEED, 'the snapshot records the soil seed of the state the run started from');
    assert.deepEqual(new Float64Array(taken.arrays.seaFlow), saved.seaFlow, 'the start places the sea\'s velocity under the air and the snapshot carries it');
    assert.deepEqual((await restore(fresh)).options, {});
    const unseeded = await snapshot();
    assert.equal(unseeded.arrays.seaFlow, undefined);
    assert.equal(unseeded.soilSeed, null, "a restore that starts over takes the snapshot's record");
    assert.deepEqual((await restore(taken)).options, OPTIONS);
    assert.deepEqual(new Float64Array((await snapshot()).arrays.seaFlow), saved.seaFlow, 'the restore places it again');
    const { soilSeed: recorded, ...unrecorded } = taken;
    assert.equal(recorded, SEED);
    await restore(unrecorded);
    assert.equal((await snapshot()).soilSeed, null, "a restore into the running model takes the snapshot's record");
    await restore(taken);
    assert.equal((await snapshot()).soilSeed, SEED);
    const { placeCarried } = await import('../js/model.worker.js');
    const placed = createModel(new Grid(16), { topography: syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 ? 1500 + 1000 * Math.cos(lat) : -4000)) }), { nCells: C } = placed.mesh, g = placed.core.diagnostics.g;
    const kept = { boundaryDepth: Float64Array.from({ length: C }, (_, i) => 300 + (i % 700)), mixingTop: Float64Array.from({ length: C }, (_, i) => 200 + (i % 500)), boundaryRegime: Float64Array.from({ length: C }, (_, i) => i % 4) };
    placeCarried(placed, kept);
    let highest = 0;
    for (let i = 0; i < C; i++) {
      const ground = placed.surfaceGeopotential[i] / g;
      highest = Math.max(highest, ground);
      assert.equal(placed.boundaryLayer.depth[i], kept.boundaryDepth[i] + ground, `cell ${i}: the CPU holds the depth above sea level, the save above the ground`);
      assert.equal(placed.boundaryLayer.mixingTop[i], kept.mixingTop[i] + ground, `cell ${i}: the mixing top likewise`);
      assert.equal(placed.boundaryLayer.regime[i], kept.boundaryRegime[i]);
    }
    assert.ok(highest > 500, `the highest ground stands ${highest} m`);
    const flux = Float64Array.from({ length: C }, (_, i) => (i % 9) - 4);
    placeCarried(placed, { iceOceanFlux: flux, ocean: {} });
    assert.deepEqual(placed.seaIce.oceanFlux, flux, 'a snapshot whose ocean has no restart flux places the top-level iceOceanFlux it carried before');
    placeCarried(placed, { iceOceanFlux: new Float64Array(C), ocean: { flux } });
    assert.deepEqual(placed.seaIce.oceanFlux, flux, 'the ocean\'s restart flux takes precedence');
  } finally {
    globalThis.fetch = fetched;
    delete globalThis.self;
    rmSync(dir, { recursive: true, force: true });
  }
});

test("a page save waits for the step that ends an ocean step, a held radiation interval and a held drag interval, and a model that has run 3 steps, restored from it through encodeState, decodeState and the page's restore, steps on as the run that did not stop over two ocean cadences and a step, with runoff and discharge in the ocean's totals at the save and more added after it", async () => {
  const stand = !globalThis.self;
  if (stand) globalThis.self = { postMessage() {} };
  try {
    const { snapshotArrays, placeSaved, stepsToCadence } = await import('../js/model.worker.js');
    const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 ? 1500 + 1000 * Math.cos(lat) : -4000));
    const options = { topography, radiation: { radiationEvery: 4 }, dragEvery: 2, ocean: { everySteps: 4 } };
    const build = () => { const model = createModel(new Grid(16), options), init = initializeState(model, {}); for (let a = 0; a < init.length; a++) model.state[a].set(init[a]); return model; };
    const run = build(), dt = 1350, names = ['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'];
    assert.deepEqual(run.cadence, { radiationEvery: 4, dragEvery: 2, oceanEvery: 4 });
    const water = async (model, change) => { const held = await model.readFreshwater(); for (const name of ['runoff', 'discharge']) change(held[`${name}Total`], held[`${name}Seen`], name); model.writeFreshwater(held); };
    await water(run, (total, seen, name) => { total.fill(name === 'runoff' ? 2 : 3); seen.fill(1); });
    for (let s = 0; s < 3; s++) run.step(dt);
    const wait = stepsToCadence(run, dt);
    assert.equal(wait, 1, 'a save asked for after 3 steps waits for the fourth');
    for (let s = 0; s < wait; s++) run.step(dt);
    assert.equal(stepsToCadence(run, dt), 0);
    const { arrays, ocean, land } = await snapshotArrays(run);
    const saved = { N: 16, K: run.core.K, day: run.time / 86400, time: run.time, terrain: true };
    for (const [key, buffer] of Object.entries(arrays)) saved[key] = new Float64Array(buffer.slice(0));
    saved.ocean = { ...groupHeader(ocean), ...Object.fromEntries(Object.entries(groupArrays(ocean)).map(([key, buffer]) => [key, new Float64Array(buffer.slice(0))])) };
    saved.land = Object.fromEntries(Object.entries(land).map(([key, buffer]) => [key, new Float64Array(buffer.slice(0))]));
    const back = await decodeState(encodeState(saved, { f64: f64Names(saved) }));
    for (const key of ['runoffTotal', 'dischargeTotal', 'runoffSeen', 'dischargeSeen']) assert.ok(Math.min(...back[key]) > 0, `${key} is not 0 at the save`);
    const carried = ['boundaryDepth', 'mixingTop', 'boundaryRegime', 'boundaryBuoyancy', 'subcloudVirtual', 'windSpeed', 'evaporation', 'exchangeHeat', 'exchangeWind', 'runoffTotal', 'dischargeTotal', 'runoffSeen', 'dischargeSeen', 'rainStep'];
    const restart = ['Q', 'W', 'previousT0', 'previousIceVolume', 'fresh', 'calved', 'buoyancy', 'capacity', 'penetrated', 'flux'];
    for (const key of carried) assert.ok(back[key] && back[key].length > 0, `the snapshot carries ${key}`);
    for (const key of restart) assert.ok(back.ocean[key] && back.ocean[key].length > 0, `the snapshot carries the ocean's restart array ${key}`);
    const restored = async (drop) => {
      const model = build(), kept = { ...back, ocean: { ...back.ocean } };
      for (let s = 0; s < 3; s++) model.step(dt);
      for (const key of drop) { delete kept[key]; delete kept.ocean[key]; }
      await placeSaved(model, kept, 16, dt);
      return model;
    };
    const twin = await restored([]), bare = await restored([...carried, ...restart]);
    const gap = (model, a) => { let d = 0, s = 0; for (let i = 0; i < run.state[a].length; i++) { d = Math.max(d, Math.abs(model.state[a][i] - run.state[a][i])); s = Math.max(s, Math.abs(run.state[a][i])); } return s > 0 ? d / s : d; };
    const oceanGap = (model) => Math.max(...['h', 'u', 'Q', 'W'].map((key) => run.ocean[key].reduce((d, v, i) => Math.max(d, Math.abs(model.ocean[key][i] - v)), 0)));
    const steps = 2 * run.cadence.oceanEvery + 1, worst = names.map(() => 0);
    let bareWorst = 0;
    for (const model of [run, twin, bare]) await water(model, (total) => { for (let i = 0; i < total.length; i++) total[i] += 3; });
    for (let s = 0; s < steps; s++) {
      for (const model of [run, twin, bare]) model.step(dt);
      names.forEach((name, a) => { worst[a] = Math.max(worst[a], gap(twin, a)); bareWorst = Math.max(bareWorst, gap(bare, a)); });
      assert.equal(oceanGap(twin), 0, `step ${s + 1}: the restored ocean parts from the continued run`);
    }
    console.log(`over ${steps} steps after a restore at step 4, the largest gap from the continued run over the field's largest value: ${names.map((name, a) => `${name} ${worst[a]}`).join(', ')}; restored without the carried fields and the ocean's restart arrays: ${bareWorst.toExponential(1)}`);
    names.forEach((name, a) => assert.equal(worst[a], 0, `${name}: the restored run parts from the continued run by ${worst[a]}`));
    assert.ok(bareWorst > 1e-6, `without the carried fields and the ocean's restart arrays the run parts by only ${bareWorst}`);
  } finally {
    if (stand) delete globalThis.self;
  }
});
