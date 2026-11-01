import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';

/*
 * The GPU model with the ocean's iceExchange 'turbulent' (physics/ice.module.js): the layout of the
 * physics buffer and the kernels carry the exchange only when it is asked for, the model steps with
 * ice, the water under the ice holds heat above freezing, and the ocean's save and load carry the
 * water's heat and the conductance through the device bit for bit.
 */
let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};
const { layoutFor } = gpuAvailable ? await import('../js/gpu/core.gpu.js') : {};
const { physicsSurfaceKernel, seaSurfaceWgsl, SEA_SURFACE_WGSL } = gpuAvailable ? await import('../js/gpu/physics.gpu.js') : {};
const { readRanges } = gpuAvailable ? await import('../js/gpu/device.module.js') : {};
const N = 16, DT = 1350;

test('the physics buffer, the surface kernel and the constants carry the turbulent exchange only when it is asked for', { skip: !gpuAvailable && 'webgpu not installed' }, () => {
  const mesh = buildMesh(new Grid(N));
  const plain = layoutFor(mesh, 8), asked = layoutFor(mesh, 8, 0, 0, 0, null, false, 0, 0, false, false, false, false, false, 0, false, true);
  assert.ok(plain.PH.WHEAT === undefined && plain.PH.CONDUCT === undefined);
  assert.equal(asked.PH.WHEAT, plain.PH.total);
  assert.equal(asked.PH.CONDUCT, plain.PH.total + mesh.nCells);
  assert.equal(asked.PH.total, plain.PH.total + 2 * mesh.nCells);
  assert.equal(seaSurfaceWgsl(false), SEA_SURFACE_WGSL);
  assert.equal(physicsSurfaceKernel(false), physicsSurfaceKernel(false, false, null, false));
  assert.ok(!/PH_WHEAT|PH_CONDUCT|XFER_/.test(physicsSurfaceKernel(false)) && /PH_WHEAT/.test(physicsSurfaceKernel(false, false, null, true)));
});

test('the GPU model steps with ice over a turbulent exchange: the water under the ice holds heat above freezing, nothing goes wrong, and a save and load of the ocean give the device its heat and conductance back', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const cpu = createModel(new Grid(N), { ocean: false }), init = initializeState(cpu, {});
  const model = await createGpuModel(new Grid(N), { ocean: { iceExchange: 'turbulent' } });
  const C = model.mesh.nCells, { device, buffers, layout } = model.gpu;
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  model.load();
  const read = async () => {
    const [heat, conductance, flux, cap] = await readRanges(device, buffers.PH, ['WHEAT', 'CONDUCT', 'OFLUX', 'CAP'].map((name) => ({ offset: layout.PH[name], length: C })));
    const [ice] = await readRanges(device, buffers.S, [{ offset: layout.S.ICE, length: C }]);
    return { heat: Float32Array.from(heat), conductance: Float32Array.from(conductance), flux: Float32Array.from(flux), cap: Float32Array.from(cap), ice: Float32Array.from(ice) };
  };
  const start = await read();
  assert.ok(start.ice.some((h) => h > 0), 'the start has ice');
  for (let n = 0; n < 16; n++) await model.step(DT);
  const after = await read();
  let iced = 0, held = 0, top = 0, floor = Infinity;
  for (let i = 0; i < C; i++) {
    assert.ok(Number.isFinite(after.heat[i]) && Number.isFinite(after.conductance[i]) && Number.isFinite(after.ice[i]), `cell ${i} is finite`);
    if (!(after.ice[i] > 0)) continue;
    iced++;
    assert.ok(after.heat[i] >= 0, `cell ${i}: the water under the ice holds ${after.heat[i]} J/m²`);
    if (after.heat[i] > 0) held++;
    top = Math.max(top, after.heat[i] / after.cap[i]); floor = Math.min(floor, after.conductance[i]);
  }
  console.log(`after 16 steps at N=${N}: ${iced} iced cells, ${held} with water above freezing under them (up to ${top.toFixed(3)} K), the conductance at least ${floor.toFixed(1)} W/m²/K`);
  assert.ok(iced > 0 && floor >= 0.99 * 4.08e6 * 0.006 * 5e-4, `${iced} iced cells, conductance floor ${floor}`);
  const saved = await model.ocean.serialize({ restart: true });
  assert.deepEqual(Array.from(saved.iceWater), Array.from(after.heat));
  assert.deepEqual(Array.from(saved.iceConductance), Array.from(after.conductance));
  const marked = Float64Array.from(saved.iceWater, (x, i) => (x > 0 ? 0.5 * x : x)), changed = { ...saved, iceWater: marked };
  model.ocean.load(changed, model.state[3], model.state[6]);
  const loaded = await read();
  for (let i = 0; i < C; i++) assert.equal(loaded.heat[i], Math.fround(marked[i]), `cell ${i}: the heat loaded`);
  assert.deepEqual(Array.from(loaded.conductance), Array.from(after.conductance), 'the conductance loaded');
  const instant = await createGpuModel(new Grid(N), { ocean: {} });
  assert.equal(instant.gpu.layout.PH.WHEAT, undefined, 'the default model has no water heat on the device');
  model.destroy(); instant.destroy();
});
