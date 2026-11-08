import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { geographyOptions, readTopography, FRESH_GEOGRAPHY } from '../js/geography.module.js';
import { LAKE_FIELDS } from '../js/physics/lake.module.js';
import { oceanCells } from '../js/ocean/layered.module.js';
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

async function pageSave(snapshotArrays, model) {
  const { arrays, ocean, land } = await snapshotArrays(model);
  const saved = { N, K: model.core.K, day: model.time / 86400, time: model.time, terrain: true };
  for (const [key, buffer] of Object.entries(arrays)) saved[key] = new Float64Array(buffer.slice(0));
  saved.ocean = { ...groupHeader(ocean), ...Object.fromEntries(Object.entries(groupArrays(ocean)).map(([key, buffer]) => [key, new Float64Array(buffer.slice(0))])) };
  if (land) saved.land = Object.fromEntries(Object.entries(land).map(([key, buffer]) => [key, new Float64Array(buffer.slice(0))]));
  return decodeState(encodeState(saved, { f64: f64Names(saved) }));
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
    const back = await pageSave(snapshotArrays, run);
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

test(`at N=${N} on the n1024 mask, a page save with its lakes' outflow pending, restored on the other engine through encodeState, decodeState and the page's restore, holds the outflow a restore on its own engine holds: a CPU save's lakePassed less lakeSeen on the GPU, and a GPU save's lakePending on the CPU; and it holds the rain the open sea has yet to take that a restore on its own engine holds: a CPU save's rainStep on the GPU, and a GPU save's rainTotal less rainSeen on the CPU`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const stand = !globalThis.self;
  if (stand) globalThis.self = { postMessage() {} };
  try {
    const { snapshotArrays, placeSaved } = await import('../js/model.worker.js');
    const options = { ...setup({ model: 'flake' }), radiation: { radiationEvery: 4 }, ocean: { everySteps: 4 } };
    const cpu = prepare(createModel(new Grid(N), options)), gpu = prepare(await createGpuModel(new Grid(N), options));
    for (let n = 0; n < 3; n++) { cpu.step(DT); await gpu.step(DT); }
    const fromCpu = await pageSave(snapshotArrays, cpu), fromGpu = await pageSave(snapshotArrays, gpu), lakes = cpu.lake.cells;
    assert.equal(fromCpu.lakePending, undefined, 'a CPU save holds the outflow as lakePassed less lakeSeen');
    assert.equal(fromGpu.lakeSeen, undefined, 'a GPU save holds the outflow as lakePending');
    assert.ok(lakes.some((i) => fromCpu.lakePassed[i] !== fromCpu.lakeSeen[i]) && lakes.some((i) => fromGpu.lakePending[i] !== 0), 'both saves hold lake outflow not yet drained to the sea');
    const held = async (model, saved) => {
      await placeSaved(model, saved, N, DT);
      if (model === cpu) return { lake: Array.from(cpu.lake.passed, (passed, i) => passed - cpu.lake.seen[i]), rain: Array.from(cpu.readFreshwater().rainStep) };
      const { lakePending, rainTotal, rainSeen } = await gpu.readFreshwater();
      return { lake: Array.from(lakePending), rain: Array.from(rainTotal, (total, i) => Math.fround(total >= rainSeen[i] ? total - rainSeen[i] : total)) };
    };
    const cpuOwn = await held(cpu, fromCpu), cpuOnGpu = await held(gpu, fromCpu), gpuOwn = await held(gpu, fromGpu), gpuOnCpu = await held(cpu, fromGpu);
    console.log(`N=${N}, lake cells ${lakes.join(', ')}: the CPU save's outflow ${lakes.map((i) => cpuOwn.lake[i].toExponential(3)).join(', ')} kg/m², on the GPU ${lakes.map((i) => cpuOnGpu.lake[i].toExponential(3)).join(', ')}; the GPU save's ${lakes.map((i) => gpuOwn.lake[i].toExponential(3)).join(', ')}, on the CPU ${lakes.map((i) => gpuOnCpu.lake[i].toExponential(3)).join(', ')}`);
    assert.deepEqual(cpuOnGpu.lake, cpuOwn.lake.map(Math.fround), "the GPU holds the CPU save's outflow as the CPU restored from it does, to float32");
    assert.deepEqual(gpuOnCpu.lake, gpuOwn.lake, "the CPU holds the GPU save's outflow as the GPU restored from it does");
    const { cellOpen } = oceanCells(cpu.mesh.nCells, cpu.geography), sea = (rain) => rain.filter((_, i) => cellOpen[i] && cpu.ocean.cellOcean[i]);
    assert.ok(sea(cpuOwn.rain).some((rain) => rain > 0) && sea(gpuOwn.rain).some((rain) => rain > 0), 'both saves hold rain the open sea has yet to take');
    assert.deepEqual(sea(cpuOnGpu.rain), sea(cpuOwn.rain).map(Math.fround), "on the open sea the GPU holds the CPU save's rain as the CPU restored from it does, to float32");
    assert.deepEqual(sea(gpuOnCpu.rain).map(Math.fround), sea(gpuOwn.rain), "on the open sea the CPU holds the GPU save's rain as the GPU restored from it does, to float32");
    gpu.destroy();
  } finally {
    if (stand) delete globalThis.self;
  }
});
