import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { readTopography } from '../js/geography.module.js';
import { startPlaceholders } from '../js/physics/land.module.js';
import { PREVIOUS_LEVELS } from './helpers/previousDefaults.mjs';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};
const { layoutFor } = gpuAvailable ? await import('../js/gpu/core.gpu.js') : {};

const topography = readTopography('etopo1');
const FAST = { growthTime: 3 * 3600, declineTime: 2 * 3600, snowDeclineTime: 4 * 3600, treeGrowthTime: 2 * 3600, treeDeclineTime: 3600, treeMoisture: false };
const LAND = { ...FAST, deepCapacity: 300 };
const neutral = startPlaceholders().neutral.carbon;

function prepare(model) {
  const C = model.mesh.nCells, land = model.geography.land, init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  for (let i = 0; i < C; i++) if (land[i]) model.state[6][i] = 0;
  if (model.load) model.load();
  model.ocean.initialize(model.state[3], model.state[6]);
  model.land.initialize();
  const each = (value) => Float64Array.from({ length: C }, (_, i) => (land[i] ? value(i) : 0));
  model.land.load({
    soil: each((i) => (i % 2 ? 300 : 40)), surface: each((i) => (i % 2 ? 15 : 0)), deep: each((i) => [10, 120, 200, 300][(i >> 1) % 4]),
    snow: each((i) => (model.mesh.latCell[i] > 1.0 ? 5 : 0)), vegetation: each(() => 0.5), canopy: each(() => 0.2), soilCarbon: each(() => neutral),
    seasonLength: each(() => 0.8), seasonWarmth: each(() => 4), rainMean: each(() => 3), demandMean: each(() => 4), litterMean: each(() => 0.3), decayMean: each(() => 0.2),
    record: Float64Array.of(0, 4, 4),
  });
  return model;
}

const apart = (a, b, keep) => { let worst = 0; for (let i = 0; i < a.length; i++) if (keep(i)) worst = Math.max(worst, Math.abs(a[i] - b[i])); return worst; };

test('the deep store adds one slot of C floats before the 2 m temperature, the last slot of the physics buffer, and moves no slot but that one', { skip: !gpuAvailable && 'webgpu not installed' }, () => {
  const mesh = buildMesh(new Grid(16)), C = mesh.nCells, args = [8, 0, 0, 0, null, false, 0, 0, false, false, false, false, false, 3, false, true];
  const without = layoutFor(mesh, ...args), plain = layoutFor(mesh, ...args, false), deep = layoutFor(mesh, ...args, true);
  assert.equal(without.PH.DEEP, undefined);
  assert.deepEqual(plain.PH, without.PH);
  assert.equal(deep.PH.DEEP, without.PH.T2M);
  assert.equal(deep.PH.T2M, deep.PH.total - C);
  assert.equal(deep.PH.total, without.PH.total + C);
  for (const name of Object.keys(without.PH)) if (name !== 'total') assert.equal(deep.PH[name], without.PH[name] + (name === 'T2M' ? C : 0), name);
});

test('from an observed start both engines recharge the deep store from the root zone and let the trees draw on it alike, to rounding after a step and over 16 steps in all but a handful of cells, and the GPU saves and reloads it', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const cpu = prepare(createModel(new Grid(16), { topography, land: LAND, levels: PREVIOUS_LEVELS })), gpu = prepare(await createGpuModel(new Grid(16), { topography, land: LAND, levels: PREVIOUS_LEVELS }));
  const C = cpu.mesh.nCells, grown = (i) => cpu.geography.land[i] && !cpu.geography.iceSheet[i], count = Array.from({ length: C }, (_, i) => i).filter(grown).length;
  const largest = (values) => { let worst = 0; for (let i = 0; i < C; i++) if (grown(i)) worst = Math.max(worst, Math.abs(values[i])); return worst; };
  assert.ok(gpu.gpu.layout.PH.DEEP !== undefined, 'the GPU allocates the store');
  const start = Float64Array.from(cpu.land.deepWater), loaded = await gpu.land.serialize();
  assert.ok(loaded.deep && apart(loaded.deep, start, () => true) === 0, 'the GPU takes the loaded store');
  cpu.step(900); await gpu.step(900); await gpu.sync();
  const first = await gpu.land.serialize(), firstDeep = apart(cpu.land.deepWater, first.deep, grown), firstEvaporation = apart(cpu.radiation.evaporation, gpu.radiation.evaporation, grown) / largest(cpu.radiation.evaporation);
  for (let n = 1; n < 16; n++) { cpu.step(900); await gpu.step(900); }
  await gpu.sync();
  const saved = await gpu.land.serialize(), scale = largest(cpu.radiation.evaporation), parted = [];
  let drew = 0, filled = 0, gpuDrew = 0, gpuFilled = 0;
  for (let i = 0; i < C; i++) {
    if (!grown(i)) { assert.equal(saved.deep[i], 0, `cell ${i} off the grown land`); continue; }
    if (Math.abs(cpu.land.deepWater[i] - saved.deep[i]) >= 1e-2 || Math.abs(cpu.land.soil[i] - saved.soil[i]) >= 1e-2 || Math.abs(cpu.land.canopy[i] - saved.canopy[i]) >= 1e-4 || Math.abs(cpu.radiation.evaporation[i] - gpu.radiation.evaporation[i]) >= 0.05 * scale) parted.push(i);
    if (cpu.land.deepWater[i] < start[i] - 1e-3) drew++;
    if (cpu.land.deepWater[i] > start[i] + 1e-3) filled++;
    if (saved.deep[i] < start[i] - 1e-3) gpuDrew++;
    if (saved.deep[i] > start[i] + 1e-3) gpuFilled++;
  }
  const rest = (i) => grown(i) && !parted.includes(i), worst = (i) => `store ${cpu.land.deepWater[i].toFixed(3)}/${saved.deep[i].toFixed(3)}, evaporation ${cpu.radiation.evaporation[i].toExponential(2)}/${gpu.radiation.evaporation[i].toExponential(2)}`;
  console.log(`${count} vegetated cells at N=16: after a step the store differs by ${firstDeep.toExponential(1)} kg/m² and the evaporation by ${firstEvaporation.toExponential(1)} of the largest; after 16 the engines are apart beyond the tolerances in ${parted.length} of the cells (${parted.map((i) => `${i}: ${worst(i)}`).join('; ')}), and in the rest the store differs by ${apart(cpu.land.deepWater, saved.deep, rest).toExponential(1)} kg/m², the root zone ${apart(cpu.land.soil, saved.soil, rest).toExponential(1)} kg/m², the trees ${apart(cpu.land.canopy, saved.canopy, rest).toExponential(1)} and the evaporation ${(apart(cpu.radiation.evaporation, gpu.radiation.evaporation, rest) / scale).toExponential(1)} of the largest (${scale.toExponential(1)} kg/m²/s); the store fell in ${drew} cells (GPU ${gpuDrew}) and rose in ${filled} (GPU ${gpuFilled})`);
  assert.ok(firstDeep < 1e-3 && firstEvaporation < 1e-4, `after a step: store ${firstDeep}, evaporation ${firstEvaporation} of the largest`);
  assert.ok(drew > 0.1 * count && filled > 0.2 * count && gpuDrew > 0.1 * count && gpuFilled > 0.2 * count, `${drew} cells drew on the store and ${filled} recharged it (GPU ${gpuDrew} and ${gpuFilled}), of ${count}`);
  assert.ok(parted.length <= 0.005 * count, `${parted.length} of ${count} cells apart by 1e-2 kg/m² in the store or the root zone, 1e-4 in the trees or 0.05 of the largest evaporation: ${parted.join(', ')}`);
  gpu.land.load(saved);
  const again = await gpu.land.serialize();
  assert.equal(apart(again.deep, saved.deep, () => true), 0, 'a reload keeps the store bit for bit');
  const { deep: kept, ...without } = saved;
  gpu.land.load(without);
  const started = await gpu.land.serialize();
  for (let i = 0; i < C; i++) assert.equal(started.deep[i], grown(i) ? Math.fround(Math.min(1, Math.max(0, without.soil[i] / 300)) * 300) : 0, `cell ${i}: a save without the store starts it at the root zone's fill`);
  gpu.destroy();
});
