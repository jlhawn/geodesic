import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { SKIN_DEFAULTS, surfaceShare } from '../js/physics/skin.module.js';

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

function start(model) {
  initializeState(model).forEach((values, s) => model.state[s].set(values));
  model.time = 0.25 * 86400;
  if (model.load) model.load();
  return model;
}

test('both engines step the warm layer and the cool skin alike, and the air sees the same skin', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const cpu = start(createModel(new Grid(N), { ...OPTIONS, surface: ON })), gpu = start(await createGpuModel(new Grid(N), { ...OPTIONS, surface: ON }));
  const C = cpu.mesh.nCells, PH = gpu.gpu.layout.PH, first = { warm: 0, cool: 0, ts: 0, net: 0, sensible: 0 }, worst = { ...first };
  for (let n = 0; n < STEPS; n++) {
    cpu.step(DT); await gpu.step(DT); await gpu.settle();
    const [warm, cool, net, sensible, capacity] = await readRanges(gpu.gpu.device, gpu.gpu.buffers.PH, ['WARM', 'COOL', 'SFLUX', 'SH', 'CAP'].map((name) => ({ offset: PH[name], length: C })));
    await gpu.sync();
    for (let i = 0; i < C; i++) {
      if (cpu.state[6][i] > 0 || gpu.state[6][i] > 0) continue;
      const cpuSurface = cpu.state[3][i] + cpu.skin.offset(i, cpu.seaIce.slabHeatCapacity), gpuSurface = gpu.state[3][i] + surfaceShare(capacity[i], SKIN_DEFAULTS) * warm[i] - cool[i];
      worst.warm = Math.max(worst.warm, Math.abs(warm[i] - cpu.skin.warm[i]));
      worst.cool = Math.max(worst.cool, Math.abs(cool[i] - cpu.skin.cool[i]));
      worst.ts = Math.max(worst.ts, Math.abs(gpuSurface - cpuSurface));
      worst.net = Math.max(worst.net, Math.abs(net[i] - cpu.radiation.surfaceFlux[i]) / Math.max(1, Math.abs(cpu.radiation.surfaceFlux[i])));
      worst.sensible = Math.max(worst.sensible, Math.abs(sensible[i] - cpu.radiation.sensibleHeat[i]) / Math.max(1, Math.abs(cpu.radiation.sensibleHeat[i])));
      if (n === 0) first.net = Math.max(first.net, Math.abs(net[i] - cpu.radiation.surfaceFlux[i]) - 1e-5 * Math.abs(cpu.radiation.surfaceFlux[i]));
    }
    if (n === 0) Object.assign(first, { warm: worst.warm, cool: worst.cool, ts: worst.ts, sensible: worst.sensible });
  }
  console.log(`skin parity after one step: ΔT_w ${first.warm.toExponential(1)} K, ΔT_c ${first.cool.toExponential(1)} K, T_s ${first.ts.toExponential(1)} K, net ${first.net.toExponential(1)} W/m² beyond 1e-5 relative, sensible ${first.sensible.toExponential(1)} relative`);
  console.log(`skin parity over ${STEPS} steps: ΔT_w ${worst.warm.toExponential(1)} K, ΔT_c ${worst.cool.toExponential(1)} K, T_s ${worst.ts.toExponential(1)} K, net ${worst.net.toExponential(1)} and sensible ${worst.sensible.toExponential(1)} relative`);
  assert.ok(cpu.skin.warm.some((x) => x > 0.01), 'a warm layer grows');
  assert.ok(first.warm < 1e-4 && first.cool < 1e-4 && first.ts < 1e-4, 'the first step\'s skin to 1e-4 K');
  assert.ok(first.net < 0.1 && first.sensible < 1e-3, 'the first step\'s fluxes to 0.1 W/m² and 1e-3');
  assert.ok(worst.warm < 1e-4 && worst.cool < 1e-2 && worst.ts < 1e-2, 'the warm layer to 1e-4 K and the surface to 1e-2 K over every step');
  assert.ok(worst.net < 5e-2 && worst.sensible < 1e-2, 'the fluxes to 5e-2 and 1e-2 over every step');
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
