import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { groupArrays } from '../js/stateFile.module.js';
import { decodeMask } from '../js/geography.module.js';
import { messageState, toSnapshot, restoreMessage } from '../js/snapshots.module.js';

const LAKES = [303, 933];

function lakeMask(bytes) {
  const mask = decodeMask(bytes), { flags, landFraction, elevation, wetDepth, edgeSill, ice, sea, lakes } = mask;
  const index = Uint32Array.from([...lakes.index, ...LAKES]), level = Int16Array.from([...lakes.level, ...LAKES.map((i) => elevation[i])]);
  for (const i of LAKES) { flags[i] = 4; landFraction[i] = 16384; elevation[i] = 0; wetDepth[i] = 0; }
  const head = new DataView(new ArrayBuffer(48));
  [71, 77, 75, 49].forEach((b, j) => head.setUint8(j, b));
  [mask.N, mask.relax, mask.nCells, mask.nEdges, ice.index.length, sea.index.length, sea.edges.length, index.length].forEach((v, j) => head.setUint32(4 + 4 * j, v, true));
  const parts = [head, flags, landFraction, elevation, wetDepth, edgeSill, ice.index, ice.iceFraction, ice.shelfFraction, ice.iceThickness, ice.shelfThickness, sea.index, sea.flags, sea.elevation, sea.wetDepth, sea.edges, sea.edgeSill, index, level];
  return Buffer.concat(parts.map((part) => new Uint8Array(part.buffer, part.byteOffset, part.byteLength)));
}

test("the page's store keeps every array the worker's snapshot holds, and runs restored through the worker's own snapshot and restore, the page's store and its restore message, into the running model and into a new one, step on as the run that saved, from a state loaded off the ocean's cadence with runoff and discharge pending, and with two lakes' outflow pending at every save", async () => {
  const messages = [], fetched = globalThis.fetch;
  globalThis.self = { postMessage: (message) => messages.push(message) };
  globalThis.fetch = async (url, init) => {
    if (!String(url).startsWith('file:')) return fetched(url, init);
    let bytes;
    try { bytes = readFileSync(new URL(String(url))); } catch { return new Response(null, { status: 404 }); }
    return new Response(String(url).endsWith('/geography_n1024_N16.bin') ? lakeMask(bytes) : bytes);
  };
  try {
    const workers = {};
    for (const name of ['A', 'B', 'C', 'D']) { await import(`../js/model.worker.js?page${name}`); workers[name] = self.onmessage; }
    const arrived = async (type) => { for (let wait = 0; wait < 12000 && !messages.some((m) => m.type === type); wait++) await new Promise((resolve) => setTimeout(resolve, 10)); };
    const send = async (name, data, type) => {
      messages.length = 0;
      await workers[name]({ data });
      await arrived(type);
      const failed = messages.find((m) => m.type === 'status' && String(m.text).startsWith('error'));
      assert.equal(failed, undefined, failed && failed.text);
      const found = messages.find((m) => m.type === type);
      assert.ok(found, `${name}: no ${type} message`);
      return found;
    };
    const frame = async (name) => {
      messages.length = 0;
      const running = workers[name]({ data: { type: 'resume' } });
      workers[name]({ data: { type: 'pause' } });
      await running;
      await arrived('frame');
    };
    const parted = (mine, theirs) => {
      const out = [], groups = [['', mine.arrays, theirs.arrays], ['ocean.', groupArrays(mine.ocean), groupArrays(theirs.ocean)], ['land.', mine.land ?? {}, theirs.land ?? {}]];
      for (const [prefix, a, b] of groups) {
        for (const [key, buffer] of Object.entries(b)) {
          const x = a[key] ? new Float64Array(a[key]) : null, y = new Float64Array(buffer);
          if (!x || x.length !== y.length) { out.push(`${prefix}${key} missing`); continue; }
          let gap = 0, scale = 0;
          for (let i = 0; i < y.length; i++) { gap = Math.max(gap, Math.abs(x[i] - y[i])); scale = Math.max(scale, Math.abs(y[i])); }
          if (gap > 0) out.push(`${prefix}${key} ${(scale > 0 ? gap / scale : gap).toExponential(2)}`);
        }
      }
      return out;
    };
    const lakeOutflow = (taken) => {
      const [mean, passed, seen] = ['lakeMean', 'lakePassed', 'lakeSeen'].map((key) => new Float64Array(taken.arrays[key]));
      assert.deepEqual([...mean.keys()].filter((i) => mean[i] > 0), LAKES, 'the lakes hold water at the two cells the mask gives them and nowhere else');
      return LAKES.map((i) => passed[i] - seen[i]);
    };
    const dt = 1350, options = { radiation: { radiationEvery: 4 }, ocean: { everySteps: 4 } }, start = { type: 'start', N: 16, engine: 'cpu', paused: true, stepsPerFrame: 3, options };
    const stored = (message) => structuredClone(toSnapshot(messageState(message)));
    const stepOf = (taken) => Math.round(taken.time / dt);
    await send('A', start, 'ready');
    await frame('A');
    const first = await send('A', { type: 'snapshot' }, 'snapshotData');
    assert.equal(stepOf(first), 4, 'a save asked for after 3 steps waits for the fourth');
    const record = stored(first);
    assert.deepEqual(Object.keys(record.data.arrays).sort(), Object.keys(first.arrays).sort(), 'the store keeps every array the snapshot holds');
    assert.deepEqual(Object.keys(groupArrays(record.data.ocean)).sort(), Object.keys(groupArrays(first.ocean)).sort());
    assert.deepEqual(Object.keys(record.data.land).sort(), Object.keys(first.land).sort());
    for (const key of ['lakeMixed', 'lakeMean', 'lakeBottom', 'lakeShape', 'lakeMixedDepth', 'lakeIce', 'lakeIceT', 'lakePassed', 'lakeSeen']) assert.ok(record.data.arrays[key], `the store keeps ${key}`);
    for (const pending of lakeOutflow(first)) assert.notEqual(pending, 0, "a lake's outflow is pending at the step-4 save");
    const loaded = { meta: { ...record.meta, time: 6 * dt, day: (6 * dt) / 86400 }, data: record.data };
    for (const [key, value] of [['runoffTotal', 2], ['dischargeTotal', 3], ['runoffSeen', 1], ['dischargeSeen', 1]]) new Float64Array(loaded.data.arrays[key]).fill(value);
    const load = () => ({ type: 'restore', snapshot: restoreMessage(loaded.meta, structuredClone(loaded.data)) });
    await send('B', start, 'ready');
    await frame('B');
    await send('B', load(), 'ready');
    await frame('B');
    const taken = await send('B', { type: 'snapshot' }, 'snapshotData');
    assert.equal(stepOf(taken), 12, 'loaded at step 6, the run steps to 9 and its save waits for step 12');
    for (const pending of lakeOutflow(taken)) assert.notEqual(pending, 0, "a lake's outflow is pending at the step-12 save");
    const kept = stored(taken), restore = () => ({ type: 'restore', snapshot: restoreMessage(kept.meta, structuredClone(kept.data)) });
    await send('C', { ...start, options: { ocean: options.ocean } }, 'ready');
    const anew = await send('C', load(), 'ready');
    assert.deepEqual(anew.options, options, 'the restore starts a new model with the options the snapshot records');
    await frame('C');
    const twelve = await send('C', { type: 'snapshot' }, 'snapshotData');
    assert.equal(stepOf(twelve), 12);
    assert.deepEqual(parted(twelve, taken), [], 'loaded at step 6 into a new model, the run parts at step 12 from the run that loaded it into the running model');
    await send('D', start, 'ready');
    await frame('D');
    await send('D', restore(), 'ready');
    const later = {};
    for (const name of ['B', 'C', 'D']) {
      await frame(name);
      await frame(name);
      later[name] = await send(name, { type: 'snapshot' }, 'snapshotData');
      assert.equal(stepOf(later[name]), 20);
    }
    assert.deepEqual(parted(later.C, later.B), [], 'loaded at step 6 into a new model, the run parts at step 20 from the run that loaded it into the running model');
    assert.deepEqual(parted(later.D, later.B), [], 'restored from the step-12 save into the running model, the run parts at step 20 from the run that saved');
  } finally {
    globalThis.fetch = fetched;
    delete globalThis.self;
  }
});
