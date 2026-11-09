import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, writeFileSync, mkdtempSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import { encodeState, groupHeader, groupArrays, f64Names } from '../js/stateFile.module.js';
import { DIAGNOSE_FIRST } from '../js/physics/boundaryLayer.module.js';

test('the page worker starts a state whose options record diagnoseFirst false in that order and its snapshot records it, restores the unchanged run in place, rebuilds for a snapshot in the other order, restores an old stored snapshot without a boundaryLayer group in place on a default run, and takes a start message\'s boundaryLayer over the record, each ready message carrying the order the model it built runs', async () => {
  const dir = mkdtempSync(join(tmpdir(), 'worker-order-'));
  const messages = [], fetched = globalThis.fetch;
  globalThis.self = { postMessage: (message) => messages.push(message) };
  globalThis.fetch = async (url, init) => {
    if (!String(url).startsWith('file:')) return fetched(url, init);
    try { return new Response(readFileSync(new URL(String(url)))); } catch { return new Response(null, { status: 404 }); }
  };
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
  const inPlace = () => messages.some((m) => m.type === 'status' && String(m.text).includes('restoring the snapshot'));
  const start = (from = null, options) => send({ type: 'start', N: 16, engine: 'cpu', paused: true, ...(from ? { from: pathToFileURL(from).href } : {}), ...(options ? { options } : {}) }, 'ready');
  const snapshot = () => send({ type: 'snapshot' }, 'snapshotData');
  const restore = ({ type, ...taken }) => send({ type: 'restore', snapshot: taken }, 'ready');

  const first = await start();
  assert.deepEqual([first.options.boundaryLayer, first.diagnoseFirst], [{ diagnoseFirst: DIAGNOSE_FIRST }, DIAGNOSE_FIRST]);
  const fresh = await snapshot();
  const saved = { N: fresh.N, K: fresh.K, day: fresh.day, time: fresh.time, terrain: fresh.terrain, ...(fresh.grid ? { grid: fresh.grid } : {}), options: { boundaryLayer: { diagnoseFirst: false } } };
  for (const [key, buffer] of Object.entries(fresh.arrays)) saved[key] = new Float64Array(buffer.slice(0));
  saved.ocean = { ...groupHeader(fresh.ocean), ...Object.fromEntries(Object.entries(groupArrays(fresh.ocean)).map(([key, buffer]) => [key, new Float64Array(buffer.slice(0))])) };
  if (fresh.land) saved.land = Object.fromEntries(Object.entries(fresh.land).map(([key, buffer]) => [key, new Float64Array(buffer.slice(0))]));
  const file = join(dir, 'plain_day0000.bin');
  writeFileSync(file, encodeState(saved, { f64: f64Names(saved) }));

  const plain = await start(file);
  assert.deepEqual([plain.options.boundaryLayer, plain.diagnoseFirst], [{ diagnoseFirst: false }, false]);
  const taken = await snapshot();
  assert.deepEqual(taken.options.boundaryLayer, { diagnoseFirst: false });
  const again = await restore(taken);
  assert.deepEqual([inPlace(), again.diagnoseFirst], [true, false], 'the unchanged run restores in place');
  const back = await restore(fresh);
  assert.deepEqual([inPlace(), back.options.boundaryLayer, back.diagnoseFirst], [false, { diagnoseFirst: DIAGNOSE_FIRST }, DIAGNOSE_FIRST]);
  const old = await restore({ ...fresh, options: {} });
  assert.deepEqual([inPlace(), old.diagnoseFirst], [true, DIAGNOSE_FIRST], 'an old stored snapshot whose options hold no boundaryLayer group restores in place on a run at the default');
  const off = await start(file, { boundaryLayer: false }), on = await start(file, { boundaryLayer: { diagnoseFirst: true } });
  assert.deepEqual([off.options.boundaryLayer, off.diagnoseFirst, on.options.boundaryLayer, on.diagnoseFirst], [false, null, { diagnoseFirst: true }, true]);
});
