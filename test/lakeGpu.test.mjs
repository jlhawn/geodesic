import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { geographyOptions, readTopography, FRESH_GEOGRAPHY } from '../js/geography.module.js';
import { LAKE_FIELDS } from '../js/physics/lake.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const N = 32, STEPS = 16, DT = 900, topography = readTopography('etopo2022');
const sheets = geographyOptions({ ...FRESH_GEOGRAPHY, mask: 'n1024' });
const setup = (lake) => ({ topography, geography: sheets.geography, land: sheets.land, ...(lake ? { lake } : {}) });

function prepare(model) {
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i] || model.geography.lake?.[i]) model.state[6][i] = 0;
  if (model.load) model.load();
  model.ocean.initialize(model.state[3], model.state[6]);
  model.land.initialize();
  return model;
}

const spread = (a, b, cells) => { let s = 0, r = 0; for (const i of cells) { s += (a[i] - b[i]) ** 2; r += a[i] ** 2; } return Math.sqrt(s / Math.max(r, 1e-300)); };

test(`FLake at N=${N} on the n1024 mask: both engines carry the lakes' columns, surface temperature and ice in step over ${STEPS} steps, as closely as the sea's surface`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const cpu = prepare(createModel(new Grid(N), setup({ model: 'flake' }))), gpu = prepare(await createGpuModel(new Grid(N), setup({ model: 'flake' })));
  assert.ok(cpu.lake.cells.length > 0, `lake cells at N=${N}`);
  assert.deepEqual(gpu.lake.cells, cpu.lake.cells);
  for (let n = 0; n < STEPS; n++) { cpu.step(DT); await gpu.step(DT); }
  await gpu.sync();
  await gpu.lake.read();
  const lakes = cpu.lake.cells, sea = [];
  for (let i = 0; i < cpu.mesh.nCells; i++) if (!cpu.geography.land[i] && !cpu.geography.lake[i]) sea.push(i);
  const lakeT = spread(cpu.state[3], gpu.state[3], lakes), seaT = spread(cpu.state[3], gpu.state[3], sea);
  const fields = Object.fromEntries(LAKE_FIELDS.map((name) => [name, spread(cpu.lake[name], gpu.lake[name], lakes)]));
  const ice = Math.max(...lakes.map((i) => Math.abs(cpu.state[6][i] - gpu.state[6][i])));
  console.log(`N=${N}, ${lakes.length} lake cells, ${STEPS} steps, CPU − GPU relative rms: lake surface ${lakeT.toExponential(1)} (sea ${seaT.toExponential(1)}), ${Object.entries(fields).map(([k, v]) => `${k} ${v.toExponential(1)}`).join(', ')}; largest ice difference ${ice.toExponential(1)} m`);
  for (const name of [...LAKE_FIELDS, 'passed']) for (const i of lakes) assert.ok(Number.isFinite(gpu.lake[name][i]), `${name} at ${i}: ${gpu.lake[name][i]}`);
  assert.ok(lakeT < 3 * seaT + 1e-6, `lake surface ${lakeT} against the sea's ${seaT}`);
  for (const name of ['mixed', 'mean', 'bottom']) assert.ok(fields[name] < 3 * seaT + 1e-6, `${name}: ${fields[name]}`);
  assert.ok(ice < 1e-3, `ice ${ice} m`);
});

test(`at N=${N} on the n1024 mask the GPU model runs FLake unless told the lake is off, and off it holds no lake slots`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const plain = await createGpuModel(new Grid(N), setup(null)), off = await createGpuModel(new Grid(N), setup({ model: 'off' }));
  assert.ok(plain.lake && plain.lake.cells.length >= 2 && plain.gpu.layout.PH.LAKE !== undefined);
  assert.equal(off.lake, null);
  assert.equal(off.gpu.layout.PH.LAKE, undefined);
  assert.ok(plain.gpu.layout.PH.total > off.gpu.layout.PH.total);
});
