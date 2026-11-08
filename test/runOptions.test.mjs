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
import { syntheticTopography, geographyOptions, readTopography, decodeMask, maskUrl, FRESH_GEOGRAPHY } from '../js/geography.module.js';

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

test("a page save waits for the step that ends an ocean step, a held radiation interval and a held drag interval, and a model that has run 3 steps, restored from it through encodeState, decodeState and the page's restore, steps on as the run that did not stop over two ocean cadences and a step, with runoff and discharge pending at the save as in a spun-up run, which takes them before its land adds the step's", async () => {
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
    await water(run, (total, seen, name) => { for (let i = 0; i < total.length; i++) total[i] += name === 'runoff' ? 2 : 3; });
    const { arrays, ocean, land } = await snapshotArrays(run);
    const saved = { N: 16, K: run.core.K, day: run.time / 86400, time: run.time, terrain: true };
    for (const [key, buffer] of Object.entries(arrays)) saved[key] = new Float64Array(buffer.slice(0));
    saved.ocean = { ...groupHeader(ocean), ...Object.fromEntries(Object.entries(groupArrays(ocean)).map(([key, buffer]) => [key, new Float64Array(buffer.slice(0))])) };
    saved.land = Object.fromEntries(Object.entries(land).map(([key, buffer]) => [key, new Float64Array(buffer.slice(0))]));
    const back = await decodeState(encodeState(saved, { f64: f64Names(saved) }));
    for (const name of ['runoff', 'discharge']) assert.ok(Math.min(...back[`${name}Total`].map((total, i) => total - back[`${name}Seen`][i])) > 0, `${name} is pending at every cell at the save`);
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

test("at N=16 with two lakes, a restore writes the saved lakeSeen only beside the lake it places: from a save that holds lakeSeen but not the lake's arrays, a new model and a running one each keep the lake's own passed as seen, a save with the lake's outflow pending on the GPU's slot takes it off the saved passed, and a whole save carries both", async () => {
  const stand = !globalThis.self;
  if (stand) globalThis.self = { postMessage() {} };
  try {
    const { snapshotArrays, placeSaved } = await import('../js/model.worker.js');
    const lakes = [303, 933], mask = decodeMask(readFileSync(maskUrl(16)));
    mask.lakes = { index: Uint32Array.from([...mask.lakes.index, ...lakes]), level: Int16Array.from([...mask.lakes.level, ...lakes.map((i) => mask.elevation[i])]) };
    for (const i of lakes) { mask.flags[i] = 4; mask.landFraction[i] = 16384; mask.elevation[i] = 0; mask.wetDepth[i] = 0; }
    const sheets = geographyOptions({ ...FRESH_GEOGRAPHY, mask: 'n1024' }), dt = 1350;
    const options = { topography: readTopography('etopo2022'), geography: { ...sheets.geography, maskData: mask }, land: sheets.land, lake: { model: 'flake' } };
    const build = () => {
      const model = createModel(new Grid(16), options);
      initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
      for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i] || model.geography.lake[i]) model.state[6][i] = 0;
      model.ocean.initialize(model.state[3], model.state[6]);
      model.land.initialize();
      return model;
    };
    const run = build();
    assert.deepEqual(run.lake.cells, lakes);
    for (let s = 0; s < 4; s++) run.step(dt);
    const { arrays, ocean, land } = await snapshotArrays(run);
    const saved = { N: 16, K: run.core.K, day: run.time / 86400, time: run.time, terrain: true };
    for (const [key, buffer] of Object.entries(arrays)) saved[key] = new Float64Array(buffer.slice(0));
    saved.ocean = { ...groupHeader(ocean), ...Object.fromEntries(Object.entries(groupArrays(ocean)).map(([key, buffer]) => [key, new Float64Array(buffer.slice(0))])) };
    saved.land = Object.fromEntries(Object.entries(land).map(([key, buffer]) => [key, new Float64Array(buffer.slice(0))]));
    const back = await decodeState(encodeState(saved, { f64: f64Names(saved) }));
    for (const i of lakes) assert.ok(back.lakeSeen[i] !== 0 && back.lakePassed[i] !== back.lakeSeen[i], `lake cell ${i}: passed ${back.lakePassed[i]}, seen ${back.lakeSeen[i]}`);
    const lakeless = { ...back };
    for (const key of ['lakeMixed', 'lakeMean', 'lakeBottom', 'lakeShape', 'lakeMixedDepth', 'lakeIce', 'lakeIceT', 'lakePassed']) delete lakeless[key];
    const fresh = build(), running = build();
    for (let s = 0; s < 3; s++) running.step(dt);
    const own = Float64Array.from(running.lake.passed);
    await placeSaved(fresh, lakeless, 16, dt);
    await placeSaved(running, lakeless, 16, dt);
    const pending = Float64Array.from(back.lakePassed, (value, i) => (lakes.includes(i) ? 0.25 : 0)), gpuShaped = { ...back, lakeSeen: undefined, lakePending: pending };
    const other = build(), whole = build();
    for (let s = 0; s < 3; s++) whole.step(dt);
    await placeSaved(other, gpuShaped, 16, dt);
    await placeSaved(whole, back, 16, dt);
    const seen = (model) => model.readFreshwater().lakeSeen;
    for (const i of lakes) {
      assert.equal(seen(fresh)[i], fresh.lake.passed[i], `lake cell ${i}: a new model keeps its own lake's passed as seen`);
      assert.equal(seen(running)[i], own[i], `lake cell ${i}: a running model keeps its own lake's passed as seen`);
      assert.equal(other.lake.passed[i], back.lakePassed[i]);
      assert.equal(seen(other)[i], back.lakePassed[i] - 0.25, `lake cell ${i}: the GPU's pending outflow is taken off the saved passed`);
      assert.equal(whole.lake.passed[i], back.lakePassed[i], `lake cell ${i}: a whole save carries the lake's passed`);
      assert.equal(seen(whole)[i], back.lakeSeen[i], `lake cell ${i}: a whole save carries lakeSeen`);
    }
  } finally {
    if (stand) delete globalThis.self;
  }
});
