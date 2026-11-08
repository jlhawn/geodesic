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

const OPTIONS = {
  radiation: { cloudColumns: true, cloudSpread: 0.75, varianceCover: 'water', buoyancyLength: 0.53 },
  moist: { stratocumulusInversion: 7, cumulusEntrainmentHeight: 1, cumulusAscent: 'kinetic' },
  surface: { relativeWind: true },
};

test("the page's worker runs a saved state with the options it records, and its snapshots and restores carry them", async () => {
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
    assert.equal(fresh.arrays.seaFlow, undefined);
    const saved = { N: fresh.N, K: fresh.K, day: fresh.day, time: fresh.time, terrain: fresh.terrain, ...(fresh.grid ? { grid: fresh.grid } : {}), options: OPTIONS };
    for (const [key, buffer] of Object.entries(fresh.arrays)) saved[key] = new Float64Array(buffer.slice(0));
    saved.ocean = { ...groupHeader(fresh.ocean), ...Object.fromEntries(Object.entries(groupArrays(fresh.ocean)).map(([key, buffer]) => [key, new Float64Array(buffer.slice(0))])) };
    saved.land = Object.fromEntries(Object.entries(fresh.land).map(([key, buffer]) => [key, new Float64Array(buffer.slice(0))]));
    saved.seaFlow = Float64Array.from({ length: saved.u.length / fresh.K }, (_, e) => Math.fround(0.25 * Math.sin(0.37 * e)));
    const file = join(dir, 'optioned_day0000.bin');
    writeFileSync(file, encodeState(saved, { f64: f64Names(saved) }));
    assert.deepEqual((await start(file)).options, OPTIONS);
    const taken = await snapshot();
    assert.deepEqual(taken.options, OPTIONS);
    assert.deepEqual(new Float64Array(taken.arrays.seaFlow), saved.seaFlow, 'the start places the sea\'s velocity under the air and the snapshot carries it');
    assert.deepEqual((await restore(fresh)).options, {});
    assert.equal((await snapshot()).arrays.seaFlow, undefined);
    assert.deepEqual((await restore(taken)).options, OPTIONS);
    assert.deepEqual(new Float64Array((await snapshot()).arrays.seaFlow), saved.seaFlow, 'the restore places it again');
    const { placeBoundaryLayer } = await import('../js/model.worker.js');
    const placed = createModel(new Grid(16), { topography: syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 ? 1500 + 1000 * Math.cos(lat) : -4000)) }), { nCells: C } = placed.mesh, g = placed.core.diagnostics.g;
    const kept = { boundaryDepth: Float64Array.from({ length: C }, (_, i) => 300 + (i % 700)), mixingTop: Float64Array.from({ length: C }, (_, i) => 200 + (i % 500)), boundaryRegime: Float64Array.from({ length: C }, (_, i) => i % 4) };
    placeBoundaryLayer(placed, kept);
    let highest = 0;
    for (let i = 0; i < C; i++) {
      const ground = placed.surfaceGeopotential[i] / g;
      highest = Math.max(highest, ground);
      assert.equal(placed.boundaryLayer.depth[i], kept.boundaryDepth[i] + ground, `cell ${i}: the CPU holds the depth above sea level, the save above the ground`);
      assert.equal(placed.boundaryLayer.mixingTop[i], kept.mixingTop[i] + ground, `cell ${i}: the mixing top likewise`);
      assert.equal(placed.boundaryLayer.regime[i], kept.boundaryRegime[i]);
    }
    assert.ok(highest > 500, `the highest ground stands ${highest} m`);
  } finally {
    globalThis.fetch = fetched;
    delete globalThis.self;
    rmSync(dir, { recursive: true, force: true });
  }
});

test("a page snapshot carries what the first step reads before it computes its own through encodeState and decodeState to the page's restore, so that the first step after a restore matches the continued run", async () => {
  const stand = !globalThis.self;
  if (stand) globalThis.self = { postMessage() {} };
  try {
    const { snapshotArrays, placeSaved } = await import('../js/model.worker.js');
    const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 ? 1500 + 1000 * Math.cos(lat) : -4000));
    const build = () => { const model = createModel(new Grid(16), { topography }), init = initializeState(model, {}); for (let a = 0; a < init.length; a++) model.state[a].set(init[a]); return model; };
    const run = build(), dt = 1350, C = run.mesh.nCells, carried = ['boundaryDepth', 'mixingTop', 'boundaryRegime', 'boundaryBuoyancy', 'subcloudVirtual', 'iceOceanFlux'];
    for (let s = 0; s < 8; s++) run.step(dt);
    const { arrays, ocean, land } = await snapshotArrays(run);
    const saved = { N: 16, K: run.core.K, day: run.time / 86400, time: run.time, terrain: true };
    for (const [key, buffer] of Object.entries(arrays)) saved[key] = new Float64Array(buffer.slice(0));
    saved.ocean = { ...groupHeader(ocean), ...Object.fromEntries(Object.entries(groupArrays(ocean)).map(([key, buffer]) => [key, new Float64Array(buffer.slice(0))])) };
    saved.land = Object.fromEntries(Object.entries(land).map(([key, buffer]) => [key, new Float64Array(buffer.slice(0))]));
    const back = await decodeState(encodeState(saved, { f64: f64Names(saved) }));
    for (const key of carried) assert.equal(back[key]?.length, key === 'subcloudVirtual' ? run.moist.subcloudVirtual.length : C, `the snapshot carries ${key}`);
    const restored = async (drop) => { const model = build(), kept = { ...back }; for (const key of drop) delete kept[key]; await placeSaved(model, kept, 16); return model; };
    const twin = await restored([]), bare = await restored(carried);
    const g = run.core.diagnostics.g;
    for (let i = 0; i < C; i++) {
      const ground = run.surfaceGeopotential[i] / g;
      assert.ok(Math.abs(back.boundaryDepth[i] + ground - run.boundaryLayer.depth[i]) <= 1e-12 * Math.abs(run.boundaryLayer.depth[i]), `cell ${i}: the save holds the depth above the ground`);
      for (const [mine, theirs] of [[twin.boundaryLayer.depth, run.boundaryLayer.depth], [twin.boundaryLayer.mixingTop, run.boundaryLayer.mixingTop]]) assert.ok(Math.abs(mine[i] - theirs[i]) <= 1e-12 * Math.abs(theirs[i]), `cell ${i}: a height comes back to round-off`);
      for (const [mine, theirs] of [[twin.boundaryLayer.regime, run.boundaryLayer.regime], [twin.boundaryLayer.buoyancyFlux, run.boundaryLayer.buoyancyFlux], [twin.seaIce.oceanFlux, run.seaIce.oceanFlux]]) assert.equal(mine[i], theirs[i], `cell ${i}`);
    }
    assert.deepEqual(twin.moist.subcloudVirtual, run.moist.subcloudVirtual);
    const names = ['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'];
    for (const model of [run, twin, bare]) model.step(dt);
    const gap = (model, a) => { let d = 0, s = 0; for (let i = 0; i < run.state[a].length; i++) { d = Math.max(d, Math.abs(model.state[a][i] - run.state[a][i])); s = Math.max(s, Math.abs(run.state[a][i])); } return s > 0 ? d / s : d; };
    const restoredGaps = names.map((_, a) => gap(twin, a)), bareGaps = names.map((_, a) => gap(bare, a));
    console.log(`the first step after a restore against the continued run (largest gap over the field's largest value): ${names.map((name, a) => `${name} ${restoredGaps[a].toExponential(1)}`).join(', ')}; restored without the carried fields: ${names.map((name, a) => `${name} ${bareGaps[a].toExponential(1)}`).join(', ')}`);
    names.forEach((name, a) => assert.ok(restoredGaps[a] <= 1e-12, `${name}: the first step after the restore parts from the continued run by ${restoredGaps[a]}`));
    assert.ok(Math.max(...bareGaps) > 1e-3, `without the carried fields the first step parts by only ${Math.max(...bareGaps)}`);
  } finally {
    if (stand) delete globalThis.self;
  }
});
