import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, writeFileSync, mkdtempSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import { encodeState, groupHeader, groupArrays, f64Names } from '../js/stateFile.module.js';
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
