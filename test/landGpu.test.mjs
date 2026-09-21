import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 300 : -4000));

function stats(cpu, gpu) {
  let maxDiff = 0, at = -1, sumSq = 0, sumRef = 0;
  for (let x = 0; x < cpu.length; x++) { const d = Math.abs(cpu[x] - gpu[x]); if (d > maxDiff) { maxDiff = d; at = x; } sumSq += d * d; sumRef += cpu[x] * cpu[x]; }
  return { maxDiff, at, rmsRel: Math.sqrt(sumSq / Math.max(sumRef, 1e-300)) };
}

function prepare(model) {
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i]) model.state[6][i] = 0;
  if (model.load) model.load();
  model.ocean.initialize(model.state[3], model.state[6]);
  model.land.initialize();
  for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i]) { model.land.soil[i] = 40; model.land.snow[i] = model.mesh.latCell[i] > 1.0 ? 5 : 0; }
  model.land.load({ soil: Float64Array.from(model.land.soil), snow: Float64Array.from(model.land.snow) });
  return model;
}

test('eight GPU steps over a continent track the CPU model: surface, soil, snow, and the coast-bound ocean', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const cpu = prepare(createModel(new Grid(6), { topography }));
  const gpu = prepare(await createGpuModel(new Grid(6), { topography }));
  for (let n = 0; n < 8; n++) { cpu.step(900); await gpu.step(900); }
  await gpu.sync();
  const land = await gpu.land.serialize();
  const d = await gpu.diagnostics(), dc = cpu.diagnostics();
  const ts = stats(cpu.state[3], gpu.state[3]), theta = stats(cpu.state[1], gpu.state[1]);
  const soil = stats(cpu.land.soil, land.soil), snow = stats(cpu.land.snow, land.snow);
  console.log(`eight steps over a continent at N=6: Ts max ${ts.maxDiff.toExponential(1)} K, θ rms ${theta.rmsRel.toExponential(1)}, soil max ${soil.maxDiff.toExponential(1)} kg/m², snow max ${snow.maxDiff.toExponential(1)} kg/m²; land T ${dc.landMeanT.toFixed(2)} vs ${d.landMeanT.toFixed(2)}, soil ${dc.soilWater.toFixed(2)} vs ${d.soilWater.toFixed(2)}, ocean h1 ${dc.oceanUpperDepth.toFixed(2)} vs ${d.oceanUpperDepth.toFixed(2)}`);
  assert.ok(ts.maxDiff < 0.02, `Ts max ${ts.maxDiff} at ${ts.at}`);
  assert.ok(theta.rmsRel < 1e-4, `θ rms ${theta.rmsRel}`);
  assert.ok(soil.maxDiff < 1e-2, `soil max ${soil.maxDiff} at ${soil.at}`);
  assert.ok(snow.maxDiff < 1e-2, `snow max ${snow.maxDiff} at ${snow.at}`);
  assert.ok(Math.abs(dc.landMeanT - d.landMeanT) < 0.01 && Math.abs(dc.soilWater - d.soilWater) < 0.01 && Math.abs(dc.oceanUpperDepth - d.oceanUpperDepth) < 1e-3);
  for (let i = 0; i < cpu.mesh.nCells; i++) if (cpu.geography.land[i]) assert.equal(gpu.state[6][i], 0);
});
