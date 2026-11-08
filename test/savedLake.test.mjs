import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readdirSync, readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { geographyOptions, readTopography, FRESH_GEOGRAPHY } from '../js/geography.module.js';
import { LAKE_FIELDS, lakeArrays, savedLakeOptions, holdsLake, savedLake } from '../js/physics/lake.module.js';
import { loadSavedLake, savedLakeOff } from '../scripts/savedLand.mjs';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const N = 32, DT = 1350 * 16 / N, topography = readTopography('etopo2022');
const sheets = geographyOptions({ ...FRESH_GEOGRAPHY, mask: 'n1024' });
const setup = (lake = { model: 'flake' }) => ({ topography, geography: sheets.geography, land: sheets.land, lake });

function savedColumns(lake) {
  const C = lake.depth.length, columns = Object.fromEntries([...LAKE_FIELDS, 'passed'].map((name) => [name, new Float64Array(C)]));
  lake.cells.forEach((i, n) => {
    const D = lake.depth[i], frozen = n % 2 === 0;
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
  const saved = savedColumns(model.lake);
  assert.throws(() => loadSavedLake(prepare(createModel(new Grid(N), setup({ model: 'off' })), null), saved), /the save holds lakes and the model has none/);
});

test('the page and compareEngines build their lakes with the save\x27s lake options', () => {
  const worker = readFileSync('js/model.worker.js', 'utf8'), engines = readFileSync('scripts/compareEngines.mjs', 'utf8');
  assert.match(worker, /const lakeGiven = saved \? savedLakeOptions\(saved\) : \{\};\n\s*const given = runOptions\(saved \? \{ \.\.\.saved\.options, \.\.\.\(Object\.keys\(lakeGiven\)\.length \? \{ lake: lakeGiven \} : \{\}\) \}/);
  assert.match(worker, /missingBed === 'hydrolakes'\) options\.lake = \{ \.\.\.options\.lake, lakeDepths: decodeLakeDepths\(await \(await fetch\(lakeDepthsUrl\)\)\.arrayBuffer\(\)\) \}/);
  assert.match(engines, /const lakeOptions = \{ \.\.\.savedLakeOptions\(saved\), \.\.\.JSON\.parse\(process\.env\.LAKE \?\? '\{\}'\) \};\n\s*if \(Object\.keys\(lakeOptions\)\.length\) options\.lake = lakeOptions;/);
  assert.match(engines, /lake options: /);
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
