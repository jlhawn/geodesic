import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { PREVIOUS_SHORTWAVE } from './helpers/previousDefaults.mjs';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { saturationHumidity } from '../js/physics/moist.module.js';
import { DAY } from '../js/physics/radiation.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuCore } = gpuAvailable ? await import('../js/gpu/core.gpu.js') : {};
const N = 16, DT = 1350;

function meanTheta(model) {
  const { K } = model.core, C = model.mesh.nCells, theta = model.state[1];
  return Float64Array.from({ length: K }, (_, k) => { let s = 0, a = 0; for (let i = 0; i < C; i++) { s += model.mesh.areaCell[i] * theta[k * C + i]; a += model.mesh.areaCell[i]; } return s / a; });
}

function cloudier(model, thickness = 1) {
  const { K, sigmaMid, exnerLayer } = model.core.diagnostics, C = model.mesh.nCells, [pi, theta, , , q, qc] = model.state;
  const layers = [[0.8, 0.95, 4e-4, (i) => i % 4 === 0], [0.45, 0.6, 2e-4, (i) => i % 7 < 3], [0.2, 0.3, 6e-5, (i) => i % 5 === 1], [0.55, 0.9, 1e-4, (i) => i % 11 === 2]];
  for (const [from, to, water, test] of layers) {
    for (let k = 0; k < K; k++) {
      if (sigmaMid[k] < from || sigmaMid[k] > to) continue;
      for (let i = 0; i < C; i++) {
        if (!test(i)) continue;
        const idx = k * C + i;
        q[idx] = saturationHumidity(theta[idx] * exnerLayer[idx], pi[i] * sigmaMid[k]);
        qc[idx] = thickness * water;
      }
    }
  }
}

function stats(cpu, gpu) {
  let maxDiff = 0, sumSq = 0, sumRef = 0;
  const errors = new Float64Array(cpu.length);
  for (let x = 0; x < cpu.length; x++) { const d = Math.abs(cpu[x] - gpu[x]); errors[x] = d; if (d > maxDiff) maxDiff = d; sumSq += d * d; sumRef += cpu[x] * cpu[x]; }
  const scale = Math.sqrt(sumRef / cpu.length);
  return { maxDiff, rmsRel: Math.sqrt(sumSq / Math.max(sumRef, 1e-300)), p99Rel: errors.sort()[Math.floor(0.99 * (cpu.length - 1))] / Math.max(scale, 1e-300) };
}

async function pair(radiation, { steps = 4, day = 100, thickness = 1 } = {}) {
  const model = createModel(new Grid(N), { ocean: false, radiation: { ...PREVIOUS_SHORTWAVE, stratus: true, clearSkyPass: true, ...radiation }, surface: { exchange: 'roughness' } });
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  cloudier(model, thickness);
  model.time = day * DAY;
  const gpu = await createGpuCore(model.mesh, { nu4: model.core.nu4, nu4Theta: model.core.nu4Theta, divergenceDamping: model.core.divergenceDamping, referenceTheta: meanTheta(model), physics: { ...PREVIOUS_SHORTWAVE, stratus: true, clearSkyPass: true, ...radiation } });
  gpu.upload(model.state);
  gpu.uploadPhysics({});
  const times = Array.from({ length: steps }, (_, n) => model.time + n * DT);
  for (const time of times) { model.step(DT); await gpu.stepModel(DT, time); }
  const physics = await gpu.downloadPhysics(), state = await gpu.download();
  return { model, gpu, state, physics };
}

function measure(label, { model, state, physics }) {
  const C = model.mesh.nCells, r = model.radiation;
  const m = {
    absorbed: stats(r.summed.absorbedSolar, physics.ABSSUM.subarray(0, C)), atmosphere: stats(r.summed.atmosphereSolar, physics.ATMSUM.subarray(0, C)), reflected: stats(r.summed.reflectedSolar, physics.REFLSUM.subarray(0, C)),
    clear: stats(r.summed.clearAbsorbedSolar, physics.ABSCLRSUM.subarray(0, C)), olr: stats(r.outgoing, physics.OLR.subarray(0, C)), theta: stats(model.state[1], state[1]),
  };
  console.log(`${label}: against the CPU, 99th percentile of the cells (rms): absorbed ${m.absorbed.p99Rel.toExponential(1)} (${m.absorbed.rmsRel.toExponential(1)}), atmosphere ${m.atmosphere.p99Rel.toExponential(1)} (${m.atmosphere.rmsRel.toExponential(1)}), reflected ${m.reflected.p99Rel.toExponential(1)}, clear-sky absorbed ${m.clear.p99Rel.toExponential(1)}, OLR ${m.olr.p99Rel.toExponential(1)}, θ ${m.theta.p99Rel.toExponential(1)}`);
  return m;
}

function means(model) {
  const r = model.radiation, C = model.mesh.nCells, a = model.mesh.areaCell;
  let area = 0, absorbed = 0, atmosphere = 0, reflected = 0;
  for (let i = 0; i < C; i++) { area += a[i]; absorbed += a[i] * r.summed.absorbedSolar[i]; atmosphere += a[i] * r.summed.atmosphereSolar[i]; reflected += a[i] * r.summed.reflectedSolar[i]; }
  return { absorbed: absorbed / area, atmosphere: atmosphere / area, reflected: reflected / area };
}

function agrees(layers, sequential) {
  for (const key of ['absorbed', 'atmosphere', 'reflected', 'clear', 'olr', 'theta']) assert.ok(layers[key].p99Rel < 3 * sequential[key].p99Rel + 2e-4, `${key} 99th percentile ${layers[key].p99Rel} against ${sequential[key].p99Rel} for the sequential scheme`);
}

test('layered on the GPU matches the CPU model at N = 16 as the sequential scheme does, and moves the same sunlight', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const sequential = await pair({}), layers = await pair({ shortwave: 'layered' });
  const s = measure('sequential', sequential), t = measure('layered', layers);
  const cpu = [means(sequential.model), means(layers.model)], gpu = [sequential, layers].map(({ model, physics }) => {
    const C = model.mesh.nCells, a = model.mesh.areaCell;
    let area = 0, absorbed = 0, atmosphere = 0;
    for (let i = 0; i < C; i++) { area += a[i]; absorbed += a[i] * physics.ABSSUM[i]; atmosphere += a[i] * physics.ATMSUM[i]; }
    return { absorbed: absorbed / area, atmosphere: atmosphere / area };
  });
  console.log(`area-mean absorbed sunlight, W/m2 per step sum: CPU ${cpu[0].absorbed.toFixed(2)} to ${cpu[1].absorbed.toFixed(2)}, GPU ${gpu[0].absorbed.toFixed(2)} to ${gpu[1].absorbed.toFixed(2)}; atmosphere CPU ${cpu[0].atmosphere.toFixed(2)} to ${cpu[1].atmosphere.toFixed(2)}, GPU ${gpu[0].atmosphere.toFixed(2)} to ${gpu[1].atmosphere.toFixed(2)}`);
  assert.ok(Math.abs(cpu[1].atmosphere - cpu[0].atmosphere) > 0.5 && Math.abs(cpu[1].absorbed - cpu[0].absorbed) > 0.5, 'the layers move absorbed sunlight in the CPU model');
  assert.ok(Math.abs((gpu[1].absorbed - gpu[0].absorbed) - (cpu[1].absorbed - cpu[0].absorbed)) < 0.05 * Math.abs(cpu[1].absorbed - cpu[0].absorbed) + 0.05, 'the GPU takes the same absorbed sunlight');
  agrees(t, s);
});

test('layered on the GPU without the lognormal nodes of cloudColumns matches the CPU model, and the nodes change the sunlight alike on both', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const sequential = await pair({ cloudColumns: false }), single = await pair({ cloudColumns: false, shortwave: 'layered' }), nodes = await pair({ shortwave: 'layered' });
  const s = measure('sequential, one node', sequential), t = measure('layered, one node', single);
  const sums = ({ model, physics }, gpu) => {
    const C = model.mesh.nCells, a = model.mesh.areaCell, r = model.radiation.summed;
    let area = 0, absorbed = 0;
    for (let i = 0; i < C; i++) { area += a[i]; absorbed += a[i] * (gpu ? physics.ABSSUM[i] : r.absorbedSolar[i]); }
    return absorbed / area;
  };
  const cpu = [nodes, single].map((run) => sums(run, false)), gpu = [nodes, single].map((run) => sums(run, true));
  console.log(`area-mean absorbed sunlight with the nodes and with one node, W/m2 per step sum: CPU ${cpu[0].toFixed(2)} and ${cpu[1].toFixed(2)}, GPU ${gpu[0].toFixed(2)} and ${gpu[1].toFixed(2)}`);
  assert.ok(Math.abs(cpu[0] - cpu[1]) > 1, 'the nodes change the absorbed sunlight in the CPU model');
  assert.ok(Math.abs((gpu[0] - gpu[1]) - (cpu[0] - cpu[1])) < 0.05 * Math.abs(cpu[0] - cpu[1]) + 0.05, 'the GPU takes the same change');
  agrees(t, s);
});

test('layered on the GPU under cloud ten times as thick, where the surface may take no less than nothing at a low sun, matches the CPU model', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const options = { thickness: 10 };
  const sequential = await pair({}, options), layers = await pair({ shortwave: 'layered' }, options);
  const s = measure('sequential, thick', sequential), t = measure('layered, thick', layers);
  const C = layers.model.mesh.nCells;
  let least = Infinity;
  for (let i = 0; i < C; i++) least = Math.min(least, layers.physics.SWDN[i]);
  console.log(`least downward shortwave at the surface in the layered GPU model: ${least.toFixed(3)} W/m2`);
  assert.ok(least >= 0, 'the surface receives sunlight, not gives it');
  agrees(t, s);
});

test('layered on the GPU holds its layers in the held radiation too: a call, then the scaled steps', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const options = { steps: 6, day: 200 };
  const sequential = await pair({ radiationEvery: 4 }, options), layers = await pair({ radiationEvery: 4, shortwave: 'layered' }, options);
  const s = measure('sequential, radiationEvery 4', sequential), t = measure('layered, radiationEvery 4', layers);
  const C = layers.model.mesh.nCells, held = stats(layers.model.radiation.held.reflected, layers.physics.RADREFL.subarray(0, C)), heldAtmosphere = stats(layers.model.radiation.held.atmosphere, layers.physics.RADATM.subarray(0, C));
  console.log(`held reflected sunlight of the call: 99th percentile ${held.p99Rel.toExponential(1)}; held atmosphere ${heldAtmosphere.p99Rel.toExponential(1)}`);
  assert.ok(held.p99Rel < 5e-4 && heldAtmosphere.p99Rel < 5e-4);
  agrees(t, s);
});
