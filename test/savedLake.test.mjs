import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readdirSync, readFileSync, writeFileSync, mkdtempSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import { encodeState, decodeState, groupHeader, groupArrays, f64Names } from '../js/stateFile.module.js';
import { decodeMask, pickGeography } from '../js/geography.module.js';
import { toSnapshot, restoreMessage, exportState } from '../js/snapshots.module.js';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { geographyOptions, readTopography, FRESH_GEOGRAPHY } from '../js/geography.module.js';
import { LAKE_FIELDS, lakeArrays, savedLakeOptions, holdsLake, savedLake } from '../js/physics/lake.module.js';
import { loadSavedLake, savedLakeOff, lakeOptionsFor } from '../scripts/savedLand.mjs';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const N = 32, DT = 1350 * 16 / N, topography = readTopography('etopo2022');
const sheets = geographyOptions({ ...FRESH_GEOGRAPHY, mask: 'n1024' });
const setup = (lake = { model: 'flake' }) => ({ topography, geography: sheets.geography, land: sheets.land, lake });

function savedColumns(lake, frozenAt = (n) => n % 2 === 0) {
  const C = lake.depth.length, columns = Object.fromEntries([...LAKE_FIELDS, 'passed'].map((name) => [name, new Float64Array(C)]));
  lake.cells.forEach((i, n) => {
    const D = lake.depth[i], frozen = frozenAt(n);
    const column = frozen
      ? { mixed: 273.15, mean: 276.8, bottom: 277.1, shape: 0.6, mixedDepth: 0, ice: 0.4 + 0.01 * (n % 30), iceT: 270.5 }
      : { mixed: 281.2, mean: 279.4, bottom: 277.6, shape: 0.6, mixedDepth: Math.min(4, 0.5 * D), ice: 0, iceT: 273.15 };
    for (const name of LAKE_FIELDS) columns[name][i] = Math.fround(column[name]);
  });
  return lakeArrays({ ...columns, depth: Float64Array.from(lake.depth) });
}

function prepare(model, saved) {
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i] || model.geography.lake?.[i]) model.state[6][i] = 0;
  if (saved) for (const i of model.lake.cells) if (saved.lakeIce[i] > 0) { model.state[6][i] = saved.lakeIce[i]; model.state[3][i] = saved.lakeIceT[i]; } else model.state[3][i] = saved.lakeMixed[i];
  model.seaIce.load(model.state[6], null);
  if (model.load) model.load();
  model.ocean.initialize(model.state[3], model.state[6]);
  model.land.initialize();
  return model;
}

test(`at N=${N} loadSavedLake puts a save's lake columns on a CPU model and refuses a load that leaves them otherwise`, () => {
  const probe = createModel(new Grid(N), setup()), saved = savedColumns(probe.lake);
  assert.ok(probe.lake.cells.length > 0, `lake cells at N=${N}`);
  const cpu = prepare(createModel(new Grid(N), setup()), saved);
  assert.equal(loadSavedLake(cpu, saved), cpu.lake.cells.length);
  for (const name of LAKE_FIELDS) for (const i of cpu.lake.cells) assert.equal(cpu.lake[name][i], saved[`lake${name[0].toUpperCase()}${name.slice(1)}`][i], `${name} at ${i}`);
  const skipped = prepare(createModel(new Grid(N), setup()), saved);
  skipped.lake.load = () => {};
  assert.throws(() => loadSavedLake(skipped, saved), /lake columns differ from the save after the load \(cells off\): mixed \d+, mean \d+/);
  assert.equal(loadSavedLake(prepare(createModel(new Grid(N), setup()), null), {}), 0);
  const off = savedLakeOff((name) => (name === 'ice' ? new Float64Array(cpu.lake.ice.length) : cpu.lake[name]), { lakeIce: saved.lakeIce, ice: saved.lakeIce, mean: saved.lakeMean }, cpu.lake.cells);
  assert.deepEqual(off, [['ice', cpu.lake.cells.filter((i) => saved.lakeIce[i] > 0).length]]);
});

test(`at N=${N} both engines take a save's frozen lakes whole: the columns read back as saved, and the frozen lakes stay covered through two steps`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const probe = createModel(new Grid(N), setup()), saved = savedColumns(probe.lake);
  const cpu = prepare(createModel(new Grid(N), setup()), saved), gpu = prepare(await createGpuModel(new Grid(N), setup()), saved);
  assert.deepEqual(gpu.lake.cells, cpu.lake.cells);
  assert.equal(loadSavedLake(cpu, saved), cpu.lake.cells.length);
  assert.equal(await loadSavedLake(gpu, saved), gpu.lake.cells.length);
  const frozen = cpu.lake.cells.filter((i) => saved.lakeIce[i] > 0);
  for (let n = 1; n <= 2; n++) {
    cpu.step(DT); await gpu.step(DT); await gpu.sync(); await gpu.lake.read();
    const concentration = Float32Array.from(await gpu.seaIce.concentration);
    const split = cpu.lake.cells.filter((i) => cpu.seaIce.concentration[i] !== concentration[i]);
    const ice = Math.max(...frozen.map((i) => Math.abs(cpu.lake.ice[i] - gpu.lake.ice[i]) / saved.lakeIce[i]));
    console.log(`N=${N} step ${n}: ${cpu.lake.cells.length} lake cells, ${frozen.length} frozen, concentration split on ${split.length}, frozen ice rel max ${ice.toExponential(1)}`);
    assert.deepEqual(split, [], `step ${n}: concentration differs at ${split.slice(0, 5)}`);
    assert.ok(frozen.every((i) => cpu.seaIce.concentration[i] === 1 && concentration[i] === 1), `step ${n}: a frozen lake lost its cover`);
    assert.ok(ice < 1e-4, `step ${n}: frozen-lake ice differs by ${ice} of its thickness`);
  }
  gpu.destroy();
});

test(`at N=${N} a save's lake options come from its given record over its options, and a model built with them takes a save of other depths that a default model refuses`, () => {
  assert.deepEqual(savedLakeOptions({ options: { lake: { freshDepth: 20, salineDepth: 4 } }, given: { lake: { freshDepth: 30 } } }), { freshDepth: 30, salineDepth: 4 });
  assert.deepEqual(savedLakeOptions({ options: { radiation: {} } }), {});
  const deeper = createModel(new Grid(N), setup({ model: 'flake', freshDepth: 20 }));
  const saved = { ...savedColumns(deeper.lake), given: { lake: { freshDepth: 20 } } };
  assert.ok(deeper.lake.cells.some((i) => deeper.lake.depth[i] === 20), 'a lake at the fresh depth');
  assert.throws(() => loadSavedLake(prepare(createModel(new Grid(N), setup()), saved), saved), /give the lake options they were made with/);
  const matched = prepare(createModel(new Grid(N), setup(savedLakeOptions(saved))), saved);
  assert.equal(loadSavedLake(matched, saved), matched.lake.cells.length);
});

test(`at N=${N} a save whose lake arrays hold no lake loads as none, and a save that holds lakes refuses a model without them`, () => {
  const model = prepare(createModel(new Grid(N), setup()), null), C = model.mesh.nCells;
  const empty = lakeArrays(Object.fromEntries([...LAKE_FIELDS, 'passed', 'depth'].map((name) => [name, new Float64Array(C)])));
  assert.ok(savedLake(empty) && !holdsLake(savedLake(empty)) && !holdsLake(null));
  assert.equal(loadSavedLake(model, empty), 0);
  assert.ok(model.lake.mean.every((v) => v === 0));
  const open = savedColumns(model.lake, () => false);
  assert.ok(holdsLake(savedLake(open)), 'a save whose lakes are all open holds lakes');
  assert.equal(loadSavedLake(model, open), model.lake.cells.length);
  assert.ok(model.lake.cells.every((i) => model.lake.ice[i] === 0 && model.lake.mean[i] === Math.fround(279.4)));
  const saved = savedColumns(model.lake);
  assert.throws(() => loadSavedLake(prepare(createModel(new Grid(N), setup({ model: 'off' })), null), saved), /the save holds lakes and the model has none/);
});

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

test("the page's worker runs a spin-up save made with HydroLAKES depths with its lake options, loads its lakes at those depths, and its store, restore and export keep the options", async () => {
  const dir = mkdtempSync(join(tmpdir(), 'lake-options-')), messages = [], fetched = globalThis.fetch;
  let depthReads = 0;
  globalThis.self = { postMessage: (message) => messages.push(message) };
  globalThis.fetch = async (url, init) => {
    if (!String(url).startsWith('file:')) return fetched(url, init);
    if (String(url).endsWith('/data/lakedepth.bin')) depthReads++;
    let bytes;
    try { bytes = readFileSync(new URL(String(url))); } catch { return new Response(null, { status: 404 }); }
    return new Response(String(url).endsWith('/geography_n1024_N16.bin') ? lakeMask(bytes) : bytes);
  };
  try {
    await import('../js/model.worker.js?lakeOptions');
    const send = async (data, type) => {
      messages.length = 0;
      await self.onmessage({ data });
      for (let wait = 0; wait < 12000 && !messages.some((m) => m.type === type); wait++) await new Promise((resolve) => setTimeout(resolve, 10));
      const failed = messages.find((m) => m.type === 'status' && String(m.text).startsWith('error'));
      assert.equal(failed, undefined, failed && failed.text);
      const found = messages.find((m) => m.type === type);
      assert.ok(found, `no ${type} message`);
      return found;
    };
    const LAKE = { missingBed: 'hydrolakes', freshDepth: 20 };
    const start = (file) => send({ type: 'start', N: 16, engine: 'cpu', paused: true, ...(file ? { from: pathToFileURL(file).href } : {}) }, 'ready');
    const snapshot = () => send({ type: 'snapshot' }, 'snapshotData');
    const write = (taken, name, edit) => {
      const saved = { N: taken.N, K: taken.K, day: taken.day, time: taken.time, terrain: taken.terrain, ...(taken.grid ? { grid: taken.grid } : {}), ...pickGeography(taken), given: { lake: LAKE } };
      for (const [key, buffer] of Object.entries(taken.arrays)) saved[key] = new Float64Array(buffer.slice(0));
      saved.ocean = { ...groupHeader(taken.ocean), ...Object.fromEntries(Object.entries(groupArrays(taken.ocean)).map(([key, buffer]) => [key, new Float64Array(buffer.slice(0))])) };
      saved.land = Object.fromEntries(Object.entries(taken.land).map(([key, buffer]) => [key, new Float64Array(buffer.slice(0))]));
      edit(saved);
      const file = join(dir, `${name}_day0000.bin`);
      writeFileSync(file, encodeState(saved, { f64: f64Names(saved) }));
      return file;
    };
    await start();
    const plain = await snapshot(), plainDepth = new Float64Array(plain.arrays.lakeDepth);
    assert.deepEqual(LAKES.map((i) => plainDepth[i]), [10, 10], 'the default depth of a lake without a bed');
    const unlaked = write(plain, 'unlaked', (saved) => { for (const key of Object.keys(saved)) if (/^lake[A-Z]/.test(key)) delete saved[key]; });
    assert.deepEqual((await start(unlaked)).options.lake, LAKE, 'the run takes the lake options the save records under given');
    assert.equal(depthReads, 1, 'HydroLAKES read once');
    const depth = new Float64Array((await snapshot()).arrays.lakeDepth);
    assert.deepEqual(LAKES.map((i) => depth[i]), [20, 20], 'the lakes take the depths of the save\x27s options');
    const frozen = write(plain, 'frozen', (saved) => {
      for (const i of LAKES) Object.assign(saved, Object.fromEntries(Object.entries({ lakeMixed: 273.15, lakeMean: 276.8, lakeBottom: 277.1, lakeShape: 0.6, lakeMixedDepth: 0, lakeIce: 0.5, lakeIceT: 270.5, lakeDepth: depth[i] }).map(([key, v]) => { saved[key][i] = Math.fround(v); return [key, saved[key]]; })));
    });
    assert.deepEqual((await start(frozen)).options.lake, LAKE);
    const loaded = await snapshot(), ice = new Float64Array(loaded.arrays.lakeIce), mean = new Float64Array(loaded.arrays.lakeMean);
    assert.deepEqual(LAKES.map((i) => [ice[i], mean[i]]), LAKES.map(() => [0.5, Math.fround(276.8)]), 'the saved columns load with no depth refusal');
    assert.deepEqual(loaded.options.lake, LAKE);
    const { meta, data } = toSnapshot(await decodeState(readFileSync(frozen)));
    assert.deepEqual(meta.options.lake, LAKE, 'the store keeps the lake options');
    assert.deepEqual(restoreMessage(meta, data).options.lake, LAKE);
    assert.deepEqual(savedLakeOptions(exportState(meta, data).saved), LAKE, 'an export keeps them');
    const restored = await send({ type: 'restore', snapshot: restoreMessage(meta, structuredClone(data)) }, 'ready');
    assert.deepEqual(restored.options.lake, LAKE, 'a restore from the store runs with them');
    assert.equal(depthReads, 1, 'the decoded depths are kept across restarts');
  } finally {
    globalThis.fetch = fetched;
    rmSync(dir, { recursive: true, force: true });
  }
});

test('the tools build their lakes with the save\x27s lake options, LAKE over them', () => {
  const was = process.env.LAKE;
  try {
    delete process.env.LAKE;
    assert.deepEqual(lakeOptionsFor({ options: { lake: { freshDepth: 20, salineDepth: 4 } }, given: { lake: { freshDepth: 30 } } }), { freshDepth: 30, salineDepth: 4 });
    process.env.LAKE = '{"model":"off"}';
    assert.deepEqual(lakeOptionsFor({ given: { lake: { freshDepth: 30 } } }), { freshDepth: 30, model: 'off' });
  } finally { if (was === undefined) delete process.env.LAKE; else process.env.LAKE = was; }
  const loaders = [...readdirSync('scripts').filter((f) => f.endsWith('.mjs')).map((f) => `scripts/${f}`), 'scripts/figures/figureState.mjs'].filter((f) => /^\s*(await )?loadSavedLake\(|= loadSavedLake\(/m.test(readFileSync(f, 'utf8')));
  const own = { 'scripts/spinup.mjs': 'records and takes its own LAKE', 'scripts/explainOcean.mjs': 'resumes only its own runs', 'scripts/fixedSst.mjs': 'rebuilds every given group, LAKE among its env groups' };
  assert.deepEqual(loaders.filter((f) => !(f in own) && !/\bcreate(Gpu)?Model\([^;{]*\{\s*lake: lakeOptionsFor\(saved\), /.test(readFileSync(f, 'utf8'))), []);
  assert.match(readFileSync('scripts/fixedSst.mjs', 'utf8'), /\['lake', 'LAKE'\]/);
  assert.ok(loaders.length >= 16, `${loaders.length} loaders`);
});

const READ_ONLY = ['adjustStages', 'armMetrics', 'cellColumn', 'clearSkyBudget', 'compareStates', 'edgeBudget', 'equatorialMomentum', 'shortwaveBudget', 'snowIceAlbedo'];
const OWN_LOAD = { compareEngines: 'loads the columns itself and reports both engines against the save at step 0', hybridOceanYear: 'steps only the ocean engine' };

test('every script that builds a model from a save and steps it loads the lake columns through loadSavedLake or gpuModelFrom', () => {
  const scripts = [...readdirSync('scripts').filter((f) => f.endsWith('.mjs')).map((f) => `scripts/${f}`), ...readdirSync('scripts/figures').filter((f) => f.endsWith('.mjs')).map((f) => `scripts/figures/${f}`)];
  const fromSaves = scripts.filter((f) => { const s = readFileSync(f, 'utf8'); return /\b(decodeState|readState)\(/.test(s) && /\b(createModel|createGpuModel|gpuModelFrom)\(/.test(s); });
  const name = (f) => f.replace(/^.*\//, '').replace(/\.mjs$/, '');
  const missing = fromSaves.filter((f) => !READ_ONLY.includes(name(f)) && !(name(f) in OWN_LOAD) && !/\b(loadSavedLake|gpuModelFrom)\(/.test(readFileSync(f, 'utf8')));
  assert.deepEqual(missing, []);
  const stepping = READ_ONLY.filter((t) => /\.step\(/.test(readFileSync(`scripts/${t}.mjs`, 'utf8')));
  assert.deepEqual(stepping, [], 'a read-only script steps a model');
  assert.ok(fromSaves.length >= 30, `${fromSaves.length} scripts build a model from a save`);
  assert.match(readFileSync('scripts/figures/figureState.mjs', 'utf8'), /await loadSavedLake\(model, saved\);\n {2}return model;/);
});
