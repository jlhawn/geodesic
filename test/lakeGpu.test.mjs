import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { geographyOptions, readTopography, FRESH_GEOGRAPHY } from '../js/geography.module.js';
import { LAKE_FIELDS } from '../js/physics/lake.module.js';
import { encodeState, decodeState, groupHeader, groupArrays, f64Names } from '../js/stateFile.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const N = 32, STEPS = 16, DT = 1350 * 16 / N, topography = readTopography('etopo2022');
const sheets = geographyOptions({ ...FRESH_GEOGRAPHY, mask: 'n1024' });
const setup = (lake, surface = null) => ({ topography, geography: sheets.geography, land: sheets.land, ...(lake ? { lake } : {}), ...(surface ? { surface } : {}) });

function prepare(model) {
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i] || model.geography.lake?.[i]) model.state[6][i] = 0;
  if (model.load) model.load();
  model.ocean.initialize(model.state[3], model.state[6]);
  model.land.initialize();
  return model;
}

const spread = (a, b, cells) => { let s = 0, r = 0; for (const i of cells) { s += (a[i] - b[i]) ** 2; r += a[i] ** 2; } return Math.sqrt(s / Math.max(r, 1e-300)); };

for (const [lake, surface] of [[{ model: 'flake' }, null], [{ model: 'flake' }, { virtualDensity: true }], [{ model: 'flake', missingBed: 'hydrolakes' }, null]]) test(`FLake at N=${N} on the n1024 mask${lake.missingBed ? ' with the HydroLAKES depths' : ''}${surface ? ' under surface.virtualDensity' : ''}: both engines carry the lakes' columns, surface temperature and ice in step over ${STEPS} steps, as closely as the sea's surface`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const cpu = prepare(createModel(new Grid(N), setup(lake, surface))), gpu = prepare(await createGpuModel(new Grid(N), setup(lake, surface)));
  assert.ok(cpu.lake.cells.length > 0, `lake cells at N=${N}`);
  assert.deepEqual(gpu.lake.cells, cpu.lake.cells);
  assert.deepEqual(Array.from(gpu.lake.depth), Array.from(cpu.lake.depth));
  for (let n = 0; n < STEPS; n++) { cpu.step(DT); await gpu.step(DT); }
  await gpu.sync();
  await gpu.lake.read();
  const lakes = cpu.lake.cells, sea = [];
  for (let i = 0; i < cpu.mesh.nCells; i++) if (!cpu.geography.land[i] && !cpu.geography.lake[i]) sea.push(i);
  const lakeT = spread(cpu.state[3], gpu.state[3], lakes), seaT = spread(cpu.state[3], gpu.state[3], sea);
  const fields = Object.fromEntries(LAKE_FIELDS.map((name) => [name, spread(cpu.lake[name], gpu.lake[name], lakes)]));
  const ice = Math.max(...lakes.map((i) => Math.abs(cpu.state[6][i] - gpu.state[6][i])));
  console.log(`${lake.missingBed ? 'HydroLAKES, ' : ''}${surface ? 'virtualDensity, ' : ''}N=${N}, ${lakes.length} lake cells, ${STEPS} steps, CPU − GPU relative rms: lake surface ${lakeT.toExponential(1)} (sea ${seaT.toExponential(1)}), ${Object.entries(fields).map(([k, v]) => `${k} ${v.toExponential(1)}`).join(', ')}; largest ice difference ${ice.toExponential(1)} m`);
  for (const name of [...LAKE_FIELDS, 'passed']) for (const i of lakes) assert.ok(Number.isFinite(gpu.lake[name][i]), `${name} at ${i}: ${gpu.lake[name][i]}`);
  assert.ok(lakeT < 3 * seaT + 1e-6, `lake surface ${lakeT} against the sea's ${seaT}`);
  for (const name of ['mixed', 'mean', 'bottom']) assert.ok(fields[name] < 3 * seaT + 1e-6, `${name}: ${fields[name]}`);
  assert.ok(ice < 1e-3, `ice ${ice} m`);
  gpu.destroy();
});

test(`at N=${N} on the n1024 mask the GPU model runs FLake unless told the lake is off, and off it holds no lake slots`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const plain = await createGpuModel(new Grid(N), setup(null)), off = await createGpuModel(new Grid(N), setup({ model: 'off' }));
  assert.ok(plain.lake && plain.lake.cells.length >= 2 && plain.gpu.layout.PH.LAKE !== undefined);
  assert.equal(off.lake, null);
  assert.equal(off.gpu.layout.PH.LAKE, undefined);
  assert.ok(plain.gpu.layout.PH.total > off.gpu.layout.PH.total);
});

test(`at N=${N} on the n1024 mask, a GPU model saved through the page's snapshot with its lakes' outflow pending, and a model that has run 3 steps restored from it through encodeState, decodeState and the page's restore, step on frame by frame as the run that did not stop, in the lakes and the ocean they drain to`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const stand = !globalThis.self;
  if (stand) globalThis.self = { postMessage() {} };
  try {
    const { snapshotArrays, placeSaved, stepsToCadence } = await import('../js/model.worker.js');
    const options = { ...setup({ model: 'flake' }), radiation: { radiationEvery: 4 }, ocean: { everySteps: 4 } };
    const names = ['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'], subscription = { level: 'surface', depth: 'surface', fields: [], diagnostics: false };
    const frame = async (model, steps) => { const capturing = model.beginFrame(subscription); await model.stepBatch(steps, DT, null, false); await capturing; };
    const run = prepare(await createGpuModel(new Grid(N), options)), twin = prepare(await createGpuModel(new Grid(N), options));
    await frame(run, 3);
    await frame(run, stepsToCadence(run, DT));
    const { arrays, ocean, land } = await snapshotArrays(run);
    const saved = { N, K: run.core.K, day: run.time / 86400, time: run.time, terrain: true };
    for (const [key, buffer] of Object.entries(arrays)) saved[key] = new Float64Array(buffer.slice(0));
    saved.ocean = { ...groupHeader(ocean), ...Object.fromEntries(Object.entries(groupArrays(ocean)).map(([key, buffer]) => [key, new Float64Array(buffer.slice(0))])) };
    if (land) saved.land = Object.fromEntries(Object.entries(land).map(([key, buffer]) => [key, new Float64Array(buffer.slice(0))]));
    const back = await decodeState(encodeState(saved, { f64: f64Names(saved) }));
    const pending = run.lake.cells.filter((i) => back.lakePending[i] !== 0).length;
    assert.ok(pending > 0, 'the save holds lake outflow not yet drained to the sea');
    await frame(twin, 3);
    await placeSaved(twin, back, N, DT);
    const largest = (mine, theirs) => { let d = 0, s = 0; for (let i = 0; i < theirs.length; i++) { d = Math.max(d, Math.abs(mine[i] - theirs[i])); s = Math.max(s, Math.abs(theirs[i])); } return s > 0 ? d / s : d; };
    const worst = {}, parted = [], note = (field, gap, s) => { worst[field] = Math.max(worst[field] ?? 0, gap); if (gap > 0 && !parted.some(([name]) => name === field)) parted.push([field, s]); };
    const steps = 2 * run.cadence.oceanEvery + 1;
    for (let s = 1; s <= steps; s++) {
      await frame(run, 1);
      await frame(twin, 1);
      await Promise.all([run.sync(), twin.sync()]);
      names.forEach((name, a) => note(name, largest(twin.state[a], run.state[a]), s));
      const [mine, theirs] = await Promise.all([twin.ocean.serialize(), run.ocean.serialize()]);
      for (const key of ['h', 'u', 'T', 'S', 'eta']) note(`ocean ${key}`, largest(mine[key], theirs[key]), s);
      const [ours, others] = await Promise.all([twin.lake.serialize(), run.lake.serialize()]);
      for (const key of [...LAKE_FIELDS, 'passed']) note(`lake ${key}`, largest(ours[key], others[key]), s);
    }
    console.log(`N=${N}, ${pending} of ${run.lake.cells.length} lake cells with outflow pending at the save; over ${steps} steps after the restore, first parting steps: ${parted.map(([name, s]) => `${name} ${s}`).join(', ') || 'none'}`);
    assert.deepEqual(Object.entries(worst).filter(([, gap]) => gap > 0), [], 'the restored run parts from the run that did not stop');
    run.destroy();
    twin.destroy();
  } finally {
    if (stand) delete globalThis.self;
  }
});
