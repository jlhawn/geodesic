import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createOcean as createCpuOcean, PENETRATION_DEPTH, transmitted } from '../js/ocean/layered.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { readRanges } from '../js/gpu/device.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0.5 && Math.abs(lat) < 1.2 ? 300 : Math.cos(lon) > 0 ? -40 - 200 * (1 + Math.sin(3 * lat)) : -4000));
const HYBRID = { vertical: 'hybrid', everySteps: 1, iceDrift: false };
const N = 16, DT = 1350, CARRY_GRAINS = 1.25;
const grain = (x) => 2 ** (x ? Math.max(Math.floor(Math.log2(Math.abs(x))) - 23, -126) : -24);

test('a known pending light goes into the GPU hybrid ocean\'s stores below the top layer as on the CPU, over the deep sea and a shelf 40 to 440 m deep whose floor takes what light reaches it, each layer\'s change in h·(T − HEAT_REFERENCE) within the f32 grains of its store and the one above, the column\'s whole light kept, and what PH_PEN holds back, under 1.25 grains of the store the deposit fills last, waits there through the rest of the ocean step', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const cpuModel = createModel(new Grid(N), { topography });
  const init = initializeState(cpuModel, {});
  for (let a = 0; a < init.length; a++) cpuModel.state[a].set(init[a]);
  const { mesh, geography } = cpuModel, C = mesh.nCells;
  for (let i = 0; i < C; i++) if (geography.land[i]) cpuModel.state[6][i] = 0;
  const surfaceT = Float64Array.from(cpuModel.state[3]), ice = Float64Array.from(cpuModel.state[6]);
  const cpu = createCpuOcean(mesh, { geography, ...HYBRID });
  const gpuModel = await createGpuModel(new Grid(N), { topography, ocean: HYBRID });
  const gpu = gpuModel.oceanEngine, L = cpu.layers, rhoCp = cpu.rhoCp, { device, buffers, layout } = gpuModel.gpu;
  assert.ok(layout.PH.PEN !== undefined && layout.PH.TRANS !== undefined && rhoCp > 0, 'the GPU hybrid ocean carries PH_PEN and PH_TRANS');
  cpu.initialize(Float64Array.from(surfaceT), ice);
  gpu.initialize(Float64Array.from(surfaceT), ice);
  await gpu.upload({ ...cpu.serialize(), Q: cpu.Q, W: cpu.W, previousT0: cpu.T0, previousIceVolume: ice, fresh: new Float64Array(C), buoyancy: new Float64Array(C), flux: new Float64Array(C), capacity: cpu.capacity }, Float64Array.from(surfaceT), ice);
  const lit = [], light = new Float32Array(C), share = new Float32Array(C);
  for (let i = 0; i < C; i++) if (!geography.land[i] && ice[i] === 0 && cpu.h[i] > 0) {
    lit.push(i);
    light[i] = 1.5e5 * (1 + 0.7 * Math.sin(7.3 * i));
    share[i] = transmitted(cpu.h[i], 'IB');
  }
  assert.ok(lit.length > C / 4, `${lit.length} open sea cells take the light`);
  cpu.penetrated.set(light); cpu.transmission.set(share);
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.PEN, light);
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.TRANS, share);
  const before = await gpu.restartArrays(), cpuBefore = Float64Array.from(cpu.Q);
  cpu.readSurface(Float64Array.from(surfaceT), ice);
  gpu.readSurface(Float64Array.from(surfaceT), ice);
  const after = await gpu.restartArrays(), { h } = await gpu.download();
  let layerWorst = 0, columnWorst = 0, cpuColumnWorst = 0, carryWorst = 0, carried = 0, deep = 0, shallow = 0, shallowWorst = 0;
  for (const i of lit) {
    let gpuSum = 0, cpuSum = 0, cpuRounding = 0, above = 0;
    for (let k = 1; k < L; k++) {
      const n = k * C + i, dGpu = after.Qref[n] - before.Qref[n], dCpu = cpu.Q[n] - cpuBefore[n];
      cpuRounding += 2 * rhoCp * 2 ** -52 * (Math.abs(cpuBefore[n]) + Math.abs(cpu.Q[n]));
      const own = grain(Math.max(Math.abs(before.Qref[n]), Math.abs(after.Qref[n])));
      const arithmetic = 8 * 2 ** -24 * light[i] / (rhoCp * share[i]);
      layerWorst = Math.max(layerWorst, Math.abs(dGpu - dCpu) / (above + own + arithmetic));
      gpuSum += rhoCp * dGpu; cpuSum += rhoCp * dCpu;
      if (h[n] > 0) above = own;
    }
    columnWorst = Math.max(columnWorst, Math.abs(gpuSum + after.penetrated[i] - light[i]) / (2 * L * 2 ** -24 * light[i]));
    cpuColumnWorst = Math.max(cpuColumnWorst, Math.abs(cpuSum - light[i]) / cpuRounding);
    let z = h[i], last = 0;
    for (let k = 1; k < L; k++) {
      z += h[k * C + i];
      if (h[k * C + i] > 0) last = k;
      if (z >= PENETRATION_DEPTH && last > 0) break;
    }
    if (z >= PENETRATION_DEPTH) deep++;
    const ratio = Math.abs(after.penetrated[i]) / (rhoCp * grain(after.Qref[last * C + i]));
    carryWorst = Math.max(carryWorst, ratio);
    if (z < PENETRATION_DEPTH) { shallow++; shallowWorst = Math.max(shallowWorst, ratio); }
    if (after.penetrated[i] !== 0) carried++;
  }
  const remainder = Float64Array.from(after.penetrated);
  gpu.setStress(new Float64Array(mesh.nEdges), ice);
  for (const phase of ['dynamics', 'rescale', 'eddy', 'shift', 'regrid', 'column']) gpu.phases[phase](DT);
  gpu.salt(DT); gpu.writeSurface(DT);
  const [waiting] = await readRanges(device, buffers.PH, [{ offset: layout.PH.PEN, length: C }]);
  let moved = 0;
  for (const i of lit) if (waiting[i] !== remainder[i]) moved++;
  console.log(`${lit.length} lit cells, ${deep} reaching ${PENETRATION_DEPTH} m and ${shallow} shallower (PH_PEN there at most ${shallowWorst.toFixed(3)} grains): worst layer gap ${layerWorst.toFixed(3)} of its grains, column ${columnWorst.toFixed(3)} of its f32 bound (CPU ${cpuColumnWorst.toFixed(3)} of its stores' f64 rounding), PH_PEN at most ${carryWorst.toFixed(3)} grains of the last store, nonzero on ${carried}, changed by the rest of the step on ${moved}`);
  assert.ok(deep > lit.length / 2 && shallow >= 100, `${deep} of ${lit.length} columns reach ${PENETRATION_DEPTH} m and ${shallow} stop short of it`);
  assert.ok(layerWorst <= 1, `a layer's deposit parts by ${layerWorst} of its and the store above's grains`);
  assert.ok(columnWorst <= 1 && cpuColumnWorst <= 1, `the column keeps its light to ${columnWorst} of the f32 bound on the GPU and ${cpuColumnWorst} of its stores' f64 rounding on the CPU`);
  assert.ok(carryWorst <= CARRY_GRAINS, `PH_PEN holds back ${carryWorst} grains of the store the deposit fills last`);
  assert.ok(carried > 0, 'some light waits in PH_PEN');
  assert.equal(moved, 0, `the rest of the ocean step changes PH_PEN on ${moved} cells`);
  gpuModel.destroy();
});
