import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { geographyOptions, readTopography, FRESH_GEOGRAPHY } from '../js/geography.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const N = 16, topography = readTopography('etopo2022');
const rms = (a, b) => { let s = 0, r = 0; for (let x = 0; x < a.length; x++) { s += (a[x] - b[x]) ** 2; r += a[x] ** 2; } return Math.sqrt(s / Math.max(r, 1e-300)); };

function prepare(model) {
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i]) model.state[6][i] = 0;
  if (model.load) model.load();
  model.ocean.initialize(model.state[3], model.state[6]);
  model.land.initialize();
  return model;
}

test(`on the n1024 mask at N=${N}, shelves as land and as cavities, both engines stand on the same terrain and track each other over eight steps as closely as on the raster's mask`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const pair = async (record) => {
    const sheets = geographyOptions(record), setup = { topography, geography: sheets.geography, land: sheets.land };
    const cpu = prepare(createModel(new Grid(N), setup)), gpu = prepare(await createGpuModel(new Grid(N), setup));
    assert.deepEqual(Array.from(gpu.geography.land), Array.from(cpu.geography.land));
    assert.deepEqual(Array.from(gpu.surfaceGeopotential ?? []), Array.from(cpu.surfaceGeopotential ?? []));
    for (let n = 0; n < 8; n++) { cpu.step(900); await gpu.step(900); }
    await gpu.sync();
    return { ts: rms(cpu.state[3], gpu.state[3]), theta: rms(cpu.state[1], gpu.state[1]), pi: rms(cpu.state[0], gpu.state[0]) };
  };
  const raster = await pair(FRESH_GEOGRAPHY), mask = await pair({ ...FRESH_GEOGRAPHY, mask: 'n1024' }), cavity = await pair({ ...FRESH_GEOGRAPHY, mask: 'n1024', iceShelves: 'cavity' });
  const show = (o) => Object.entries(o).map(([k, v]) => `${k} ${v.toExponential(1)}`).join(', ');
  console.log(`N=${N}, eight steps, CPU − GPU relative rms: raster ${show(raster)}; n1024 ${show(mask)}; n1024 with cavities ${show(cavity)}`);
  for (const [name, run] of Object.entries({ mask, cavity })) for (const field of ['ts', 'theta', 'pi']) assert.ok(run[field] < 3 * raster[field] + 1e-7, `${name} ${field}: ${run[field]} against ${raster[field]} on the raster's mask`);
});
