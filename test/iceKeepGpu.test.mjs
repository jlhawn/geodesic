import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { FREEZING_POINT, MINIMUM_CONCENTRATION, NEGLIGIBLE_CONCENTRATION } from '../js/physics/ice.module.js';

/*
 * The thermodynamics' removal floors on the device against the CPU model's,
 * at N = 16: sea ice of 0.004 cover (under MINIMUM_CONCENTRATION) on a dozen
 * polar cells goes through one step of the whole model. With the sea ice's
 * `keepSmallIce` both engines keep it; without, both remove it. One GPU
 * model at a time (CLAUDE.md).
 */
let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuCore } = gpuAvailable ? await import('../js/gpu/core.gpu.js') : {};

const N = 16, DT = 900, DEG = Math.PI / 180;

function meanTheta(model) {
  const { K } = model.core, C = model.mesh.nCells, theta = model.state[1];
  return Float64Array.from({ length: K }, (_, k) => { let s = 0, a = 0; for (let i = 0; i < C; i++) { s += model.mesh.areaCell[i] * theta[k * C + i]; a += model.mesh.areaCell[i]; } return s / a; });
}

async function oneStep(keepSmallIce, cover, extra = {}) {
  const model = createModel(new Grid(N), { ocean: false, radiation: { stratus: false }, surface: { exchange: 'roughness' }, ice: { keepSmallIce, ...extra } });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  const C = model.mesh.nCells, cells = [];
  for (let i = 0; i < C && cells.length < 12; i++) if (Math.abs(model.mesh.latCell[i]) > 62 * DEG && i % 7 === 0) cells.push(i);
  for (const i of cells) { model.state[6][i] = 0.8; model.seaIce.concentration[i] = cover; model.seaIce.snow[i] = 6; model.state[3][i] = FREEZING_POINT - 3; model.seaIce.snowAlbedo[i] = 0.8; }
  const gpu = await createGpuCore(model.mesh, { nu4: model.core.nu4, nu4Theta: model.core.nu4Theta, divergenceDamping: model.core.divergenceDamping, referenceTheta: meanTheta(model), physics: { stratus: false, dragEvery: 1, keepSmallIce, ...extra } });
  gpu.upload(model.state);
  gpu.uploadPhysics({ concentration: model.seaIce.concentration, snow: model.seaIce.snow, snowAlbedo: model.seaIce.snowAlbedo });
  const time = model.time;
  model.step(DT);
  await gpu.stepModel(DT, time);
  const state = await gpu.download(), physics = await gpu.downloadPhysics();
  return { model, cells, C, thickness: state[6], concentration: physics.CONC.subarray(0, C), snow: physics.SNOW.subarray(0, C) };
}

for (const keepSmallIce of [true, false]) {
  test(`one step of the whole model with keepSmallIce ${keepSmallIce}: ice of ${keepSmallIce ? 'trace' : 'under 1 percent'} cover ${keepSmallIce ? 'stays' : 'is removed'} on the device as on the CPU`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
    const { model, cells, thickness, concentration } = await oneStep(keepSmallIce, 0.004);
    const cpuA = cells.map((i) => model.seaIce.concentration[i]), gpuA = cells.map((i) => concentration[i]);
    const cpuH = cells.map((i) => model.state[6][i]), gpuH = cells.map((i) => thickness[i]);
    console.log(`keepSmallIce ${keepSmallIce}: cover on ${cells.length} cells after a step, CPU ${cpuA.map((a) => a.toFixed(5)).join(' ')} | device ${gpuA.map((a) => a.toFixed(5)).join(' ')}`);
    assert.equal(cells.length, 12);
    for (let n = 0; n < cells.length; n++) {
      if (keepSmallIce) {
        assert.ok(cpuA[n] >= NEGLIGIBLE_CONCENTRATION && cpuH[n] > 0, `CPU cell ${cells[n]}: ${cpuA[n]} ${cpuH[n]}`);
        assert.ok(gpuA[n] >= NEGLIGIBLE_CONCENTRATION && gpuH[n] > 0, `device cell ${cells[n]}: ${gpuA[n]} ${gpuH[n]}`);
        assert.ok(Math.abs(cpuA[n] - gpuA[n]) < 1e-3 * Math.max(cpuA[n], 0.004) && Math.abs(cpuH[n] - gpuH[n]) < 1e-3 * Math.max(cpuH[n], 0.1), `cell ${cells[n]}: cover ${cpuA[n]} vs ${gpuA[n]}, thickness ${cpuH[n]} vs ${gpuH[n]}`);
      } else if (cpuA[n] < MINIMUM_CONCENTRATION) {
        assert.equal(cpuH[n], 0); assert.equal(gpuH[n], 0); assert.equal(gpuA[n], 0);
      }
    }
  });
}

test("one step of the whole model with leadClosure 'hibler': cover 0.6 under a cold sky closes as on the CPU, faster than as its open area squared", { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const squared = await oneStep(false, 0.6), hibler = await oneStep(false, 0.6, { leadClosure: 'hibler' });
  const gain = (run, n, a) => a - 0.6;
  let faster = 0;
  for (let n = 0; n < hibler.cells.length; n++) {
    const i = hibler.cells[n], cpu = hibler.model.seaIce.concentration[i], gpu = hibler.concentration[i];
    assert.ok(Math.abs(cpu - gpu) < 1e-4, `cell ${i}: cover ${cpu} vs ${gpu}`);
    assert.ok(Math.abs(hibler.model.state[6][i] - hibler.thickness[i]) < 1e-3 * hibler.model.state[6][i], `cell ${i}: thickness ${hibler.model.state[6][i]} vs ${hibler.thickness[i]}`);
    const sq = squared.model.seaIce.concentration[i];
    if (gain(hibler, n, cpu) > 0) { assert.ok(Math.abs(gain(hibler, n, cpu) - gain(squared, n, sq) / 0.4) < 1e-6 + 1e-3 * gain(hibler, n, cpu), `cell ${i}: ${cpu - 0.6} against ${(sq - 0.6) / 0.4}`); faster++; }
  }
  console.log(`leadClosure 'hibler' on ${hibler.cells.length} cells: ${faster} closing, CPU cover ${hibler.cells.map((i) => hibler.model.seaIce.concentration[i].toFixed(6)).join(' ')}`);
  assert.ok(faster > 0, 'some leads close');
});
