import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { SKIN_DEFAULTS, surfaceShare } from '../js/physics/skin.module.js';
import { syntheticTopography } from '../js/geography.module.js';

/*
 * The sea's skin on both engines at N = 8. A GPU file: run it alone,
 * with the GPU free.
 */
let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};
const { readRanges, readBuffer } = gpuAvailable ? await import('../js/gpu/device.module.js') : {};

const N = 8, DT = 1350, STEPS = 20, ON = { warmLayer: true, coolSkin: true };
const OPTIONS = { ocean: false, levels: sigmaInterfaces('bl36') };
const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 300 : -4000));
const LAYERED = { levels: OPTIONS.levels, topography };
const CASES = [
  ['over the slab from no skin', OPTIONS, 0],
  ['over the slab from a 1 K warm layer', OPTIONS, 1],
  ['over a layered ocean with 2 m mixed layers under IB light, from a 1 K warm layer', { ...LAYERED, ocean: { mixedDepth: 2, minimumThickness: 2, shallowestMixedDepth: 2 } }, 1],
  ['over a layered ocean under type III light, from a 1 K warm layer', { ...LAYERED, ocean: { penetration: 'III' } }, 1],
];

function start(model) {
  initializeState(model).forEach((values, s) => model.state[s].set(values));
  model.time = 0.25 * 86400;
  if (model.load) model.load();
  return model;
}

for (const [label, options, warmStart] of CASES) test(`${label} over the open sea, both engines step the warm layer and the cool skin alike, and the air sees the same skin${warmStart ? ' on the first step' : ''}`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const cpu = start(createModel(new Grid(N), { ...options, surface: ON })), gpu = start(await createGpuModel(new Grid(N), { ...options, surface: ON }));
  if (cpu.ocean) cpu.ocean.initialize(cpu.state[3], cpu.state[6]);
  const warm = Float64Array.from(cpu.state[6], (ice, i) => ((cpu.geography && cpu.geography.land[i]) || ice > 0 ? 0 : warmStart));
  cpu.skin.warm.set(warm);
  gpu.skin.load({ warmLayer: warm, coolSkin: new Float64Array(warm.length) });
  gpu.load();
  if (gpu.ocean) gpu.ocean.initialize(gpu.state[3], gpu.state[6]);
  const C = cpu.mesh.nCells, PH = gpu.gpu.layout.PH, first = { warm: 0, cool: 0, ts: 0, net: 0, sensible: 0, bulk: 0 }, worst = { ...first };
  const capacityOf = (i) => (cpu.ocean ? cpu.ocean.capacity[i] : cpu.seaIce.slabHeatCapacity);
  let thin = 0;
  for (let n = 0; n < STEPS; n++) {
    cpu.step(DT); await gpu.step(DT); await gpu.settle();
    const [warm, cool, net, sensible, capacity] = await readRanges(gpu.gpu.device, gpu.gpu.buffers.PH, ['WARM', 'COOL', 'SFLUX', 'SH', 'CAP'].map((name) => ({ offset: PH[name], length: C })));
    await gpu.sync();
    for (let i = 0; i < C; i++) {
      if (cpu.state[6][i] > 0 || gpu.state[6][i] > 0) continue;
      if (n === 0 && capacityOf(i) < 3 * 1025 * 3985) thin++;
      const cpuSurface = cpu.state[3][i] + cpu.skin.offset(i, capacityOf(i)), gpuSurface = gpu.state[3][i] + surfaceShare(capacity[i], SKIN_DEFAULTS) * warm[i] - cool[i];
      worst.warm = Math.max(worst.warm, Math.abs(warm[i] - cpu.skin.warm[i]));
      worst.cool = Math.max(worst.cool, Math.abs(cool[i] - cpu.skin.cool[i]));
      worst.ts = Math.max(worst.ts, Math.abs(gpuSurface - cpuSurface));
      worst.bulk = Math.max(worst.bulk, Math.abs(gpu.state[3][i] - cpu.state[3][i]));
      worst.net = Math.max(worst.net, Math.abs(net[i] - cpu.radiation.surfaceFlux[i]) / Math.max(1, Math.abs(cpu.radiation.surfaceFlux[i])));
      worst.sensible = Math.max(worst.sensible, Math.abs(sensible[i] - cpu.radiation.sensibleHeat[i]) / Math.max(1, Math.abs(cpu.radiation.sensibleHeat[i])));
      if (n === 0) first.net = Math.max(first.net, Math.abs(net[i] - cpu.radiation.surfaceFlux[i]) - 1e-5 * Math.abs(cpu.radiation.surfaceFlux[i]));
    }
    if (n === 0) Object.assign(first, { warm: worst.warm, cool: worst.cool, ts: worst.ts, sensible: worst.sensible });
  }
  console.log(`${label}: skin parity after one step: ΔT_w ${first.warm.toExponential(1)} K, ΔT_c ${first.cool.toExponential(1)} K, T_s ${first.ts.toExponential(1)} K, net ${first.net.toExponential(1)} W/m² beyond 1e-5 relative, sensible ${first.sensible.toExponential(1)} relative`);
  console.log(`${label}: ${thin} open cells under 3 m; skin parity over ${STEPS} steps, where the bulk temperatures drift ${worst.bulk.toExponential(1)} K apart: ΔT_w ${worst.warm.toExponential(1)} K, ΔT_c ${worst.cool.toExponential(1)} K, T_s ${worst.ts.toExponential(1)} K, net ${worst.net.toExponential(1)} and sensible ${worst.sensible.toExponential(1)} relative`);
  assert.ok(cpu.skin.warm.some((x) => x > 0.01), 'a warm layer grows');
  assert.ok(first.warm < 1e-4 && first.cool < 1e-4 && first.ts < 1e-4, 'the first step\'s skin to 1e-4 K');
  assert.ok(first.net < 0.1 && first.sensible < 1e-3, 'the first step\'s fluxes to 0.1 W/m² and 1e-3');
  assert.ok(worst.warm < 2e-3, 'the warm layer to 2e-3 K over every step');
  if (warmStart) return;
  assert.ok(worst.cool < 1e-2 && worst.ts < 1e-2, 'the surface to 1e-2 K over every step');
  assert.ok(worst.net < 0.1 && worst.sensible < 0.1, 'the fluxes to 0.1 relative over every step');
});

test('on the GPU, options absent and present but off give the same layout and the same buffers after 64 steps', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const absent = start(await createGpuModel(new Grid(N), OPTIONS)), off = start(await createGpuModel(new Grid(N), { ...OPTIONS, surface: { ...SKIN_DEFAULTS } }));
  assert.deepEqual(off.gpu.layout.PH, absent.gpu.layout.PH);
  assert.equal(absent.gpu.layout.PH.WARM, undefined);
  for (let n = 0; n < 64; n++) { await absent.step(DT); await off.step(DT); }
  await absent.settle(); await off.settle();
  for (const name of ['S', 'PH']) {
    const [a, b] = await Promise.all([readBuffer(absent.gpu.device, absent.gpu.buffers[name], 4 * absent.gpu.layout[name].total), readBuffer(off.gpu.device, off.gpu.buffers[name], 4 * off.gpu.layout[name].total)]);
    assert.ok(Buffer.from(a.buffer, a.byteOffset, a.byteLength).equals(Buffer.from(b.buffer, b.byteOffset, b.byteLength)), `buffer ${name}`);
  }
});

test('the GPU saves and loads its skin', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const gpu = start(await createGpuModel(new Grid(N), { ...OPTIONS, surface: ON }));
  for (let n = 0; n < 8; n++) await gpu.step(DT);
  await gpu.settle(); await gpu.sync();
  const saved = gpu.skin.serialize(), other = start(await createGpuModel(new Grid(N), { ...OPTIONS, surface: ON }));
  other.skin.load(saved);
  other.load();
  const [warm, cool] = await readRanges(other.gpu.device, other.gpu.buffers.PH, ['WARM', 'COOL'].map((name) => ({ offset: other.gpu.layout.PH[name], length: other.mesh.nCells })));
  assert.deepEqual(Array.from(warm), Array.from(Float32Array.from(saved.warmLayer)));
  assert.deepEqual(Array.from(cool), Array.from(Float32Array.from(saved.coolSkin)));
});
