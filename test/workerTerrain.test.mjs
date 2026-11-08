import { test, after } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, writeFileSync, mkdtempSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import { createModel } from '../js/model.module.js';
import { stateGeography, rebalanceSurfacePressure } from '../js/geography.module.js';
import { encodeState, decodeState, groupHeader, groupArrays, f64Names, savedLevels } from '../js/stateFile.module.js';
import { stateGrid } from '../js/stateGrid.module.js';
import { regridState, regridCellField } from '../js/physics/regrid.module.js';
import { moveOntoTerrain } from '../scripts/stateOptions.mjs';

const TERRAIN = { conserveVolume: true };
const dir = mkdtempSync(join(tmpdir(), 'worker-terrain-'));
const messages = [], fetched = globalThis.fetch;
globalThis.self = { postMessage: (message) => messages.push(message) };
globalThis.fetch = async (url, init) => {
  if (!String(url).startsWith('file:')) return fetched(url, init);
  try { return new Response(readFileSync(new URL(String(url)))); } catch { return new Response(null, { status: 404 }); }
};
await import('../js/model.worker.js');
after(() => {
  globalThis.fetch = fetched;
  delete globalThis.self;
  rmSync(dir, { recursive: true, force: true });
});

async function send(data, type) {
  messages.length = 0;
  await self.onmessage({ data });
  for (let wait = 0; wait < 12000 && !messages.some((m) => m.type === 'frame'); wait++) await new Promise((resolve) => setTimeout(resolve, 10));
  const failed = messages.find((m) => m.type === 'status' && String(m.text).startsWith('error'));
  assert.equal(failed, undefined, failed && failed.text);
  const found = messages.find((m) => m.type === type);
  assert.ok(found, `no ${type} message`);
  return found;
}
const start = (N, file = null, extra = {}) => send({ type: 'start', N, engine: 'cpu', paused: true, ...(file ? { from: pathToFileURL(file).href } : {}), ...extra }, 'ready');
const snapshot = () => send({ type: 'snapshot' }, 'snapshotData');
const restore = ({ type, ...taken }) => send({ type: 'restore', snapshot: taken }, 'ready');
const pressure = (taken) => new Float64Array(taken.arrays.pi);
const rebuilt = () => messages.some((m) => m.type === 'status' && m.text === 'building the model…');
const copy = (buffers) => Object.fromEntries(Object.entries(buffers).map(([key, buffer]) => [key, new Float64Array(buffer.slice(0))]));

function saveFile(name, { type, arrays, ocean, land, ...header }, records = {}) {
  const saved = { ...header, ...records, ...copy(arrays), ocean: { ...groupHeader(ocean), ...copy(groupArrays(ocean)) }, land: copy(land) };
  const file = join(dir, `${name}_day0000.bin`);
  writeFileSync(file, encodeState(saved, { f64: f64Names(saved) }));
  return file;
}
const readSave = async (file) => decodeState(new Uint8Array(readFileSync(file)));

test('a page without terrain keeps its model on the restore of a snapshot over terrain and one over the default terrain on the restore of one without, each moving the surface pressure as a start from that state does, and as moveOntoTerrain does for spinup.mjs and fixedSst.mjs', async () => {
  await start(16);
  const fresh = await snapshot();
  assert.equal(fresh.terrain, true);
  await start(16, saveFile('terrain', fresh), { terrain: false });
  const lowered = await snapshot();
  assert.equal(lowered.terrain, false);
  assert.ok(pressure(lowered).some((p, i) => p !== pressure(fresh)[i]), 'a start without terrain moves the pressure down to sea level');
  await restore(fresh);
  assert.equal(rebuilt(), false);
  assert.deepEqual(pressure(await snapshot()), pressure(lowered));
  const flatFile = saveFile('flat', lowered);
  await start(16, flatFile);
  const raised = await snapshot();
  assert.equal(raised.terrain, true);
  await restore(lowered);
  assert.equal(rebuilt(), false);
  assert.deepEqual(pressure(await snapshot()), pressure(raised));
  const flat = await readSave(flatFile), { topography, geography } = stateGeography(flat);
  const model = createModel(stateGrid(flat, topography).mesh, { topography, geography, physics: false, levels: savedLevels(flat) });
  const state = [Float64Array.from(flat.pi), Float64Array.from(flat.theta)];
  assert.equal(moveOntoTerrain(flat, true, model, state, model.surfaceGeopotential), false);
  assert.deepEqual(state[0], pressure(raised));
});

test('a state saved over terrain of other options and started at another resolution has its surface pressure moved from that terrain\'s surface on its own geography, regridded, onto the page\'s', async () => {
  await start(16);
  const fresh = await snapshot();
  for (const [name, terrain] of [['default', true], ['conserveVolume', TERRAIN]]) {
    const file = saveFile(`regrid-${name}`, fresh, terrain === true ? {} : { terrainOptions: terrain });
    await start(32, file);
    const taken = await snapshot(), saved = await readSave(file), levels = savedLevels(saved), { topography, geography } = stateGeography(saved);
    const source = createModel(stateGrid(saved, topography).mesh, { physics: false, levels, terrain, topography, geography });
    const target = createModel(stateGrid(taken, topography).mesh, { physics: false, levels, terrain, topography, geography });
    const carried = regridState(source, target, ['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'].map((key) => Float64Array.from(saved[key])), null, { land: saved.land });
    target.core.diagnose(carried[0], carried[1], null, null);
    rebalanceSurfacePressure(target.core, carried[0], carried[1], regridCellField(source, target, source.surfaceGeopotential), target.surfaceGeopotential);
    assert.deepEqual(pressure(taken), carried[0], name);
  }
});
