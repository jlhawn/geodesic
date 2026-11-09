import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { sigmaInterfaces, SIGMA_GRIDS } from '../js/dynamics/sigmaCore.module.js';
import { PHYSICS_KERNELS } from '../js/gpu/physics.gpu.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { getDevice, createKernel, storageBuffer, emptyBuffer, readBuffer } = gpuAvailable ? await import('../js/gpu/device.module.js') : {};
const N = 16, LEVELS = ['bl48', 'bl60'].filter((name) => SIGMA_GRIDS.includes(name));
const ULP = 2 ** -23, THRESHOLD = 1 + 1e-6, F32_THRESHOLD = Math.fround(THRESHOLD), SLACK = 6.5 * ULP;
const BAND = [Math.min(F32_THRESHOLD * (1 - SLACK), THRESHOLD), Math.max(F32_THRESHOLD * (1 + SLACK), THRESHOLD)];

function adjustmentKernel(K, C) {
  const source = PHYSICS_KERNELS.adjustRain, from = source.indexOf('// dry convective adjustment'), to = source.indexOf('saturateOutOfLine(i, pi, dt);', from);
  assert.ok(from > 0 && to > from && source.indexOf('// dry convective adjustment', from + 1) < 0, 'adjustRain holds one dry adjustment, before its out-of-line saturation');
  return `const K = ${K}; const C = ${C}; const S_TH = 0; const S_Q = ${K * C}; const S_QC = ${2 * K * C}; const D_EXM = 0; const L_DS = 0;
@group(0) @binding(0) var<storage, read_write> IN: array<f32>;
@group(0) @binding(1) var<storage, read> D: array<f32>;
@group(0) @binding(2) var<storage, read> LV: array<f32>;
@group(0) @binding(3) var<storage, read_write> OUT: array<f32>;
@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= C) { return; }
  ${source.slice(from, to)}
  OUT[i] = f32(merges);
}`;
}

async function bothEngines(model, columns) {
  const { K, C, dSigma } = model.core.diagnostics, pi = Float64Array.from({ length: C }, (_, i) => 9.6e4 + 6e3 * i / C);
  const theta32 = new Float32Array(K * C), q32 = new Float32Array(K * C), qc32 = new Float32Array(K * C);
  for (let i = 0; i < C; i++) { const column = columns(i); for (let k = 0; k < K; k++) { theta32[k * C + i] = column.theta[k]; q32[k * C + i] = 5e-3 * (k + 1) / K; qc32[k * C + i] = 1e-5 * ((7 * k + i) % 5); } }
  const theta = Float64Array.from(theta32), q = Float64Array.from(q32), qc = Float64Array.from(qc32), exner = new Float32Array(K * C), cpuMixes = new Int32Array(C);
  for (let i = 0; i < C; i++) {
    cpuMixes[i] = model.surface.convectiveAdjustColumn(i, pi, theta, q, qc);
    for (let k = 0; k < K; k++) exner[k * C + i] = model.core.diagnostics.exnerLayer[k * C + i];
  }
  const { device } = await getDevice(), kernel = createKernel(device, adjustmentKernel(K, C), 'dry adjustment');
  const state = storageBuffer(device, Float32Array.from([...theta32, ...q32, ...qc32])), out = emptyBuffer(device, 4 * C);
  const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
  kernel.run(pass, [state, storageBuffer(device, exner), storageBuffer(device, Float32Array.from(dSigma)), out], C);
  pass.end();
  device.queue.submit([encoder.finish()]);
  const after = await readBuffer(device, state, 12 * K * C), gpuMixes = await readBuffer(device, out, 4 * C);
  return { K, C, input: theta32, cpu: { theta, q, qc, mixes: cpuMixes }, gpu: { theta: after.subarray(0, K * C), q: after.subarray(K * C, 2 * K * C), qc: after.subarray(2 * K * C), mixes: gpuMixes } };
}

const near = (a, b, terms, scale) => Math.abs(a - b) <= (terms + 6) * ULP * scale;

for (const levels of LEVELS) test(`on ${levels} the dry adjustment's merge of two layers decides alike on both engines for every ratio outside the f32 rounding of its 1 + 1e-6 threshold, mixes the pair alike, and runs a cascade of up to seven merges in one push alike`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const model = createModel(new Grid(N), { ocean: false, levels: sigmaInterfaces(levels) }), { K, C } = model.core.diagnostics;
  const stable = (i, k, top = K - 1) => Math.fround((240 + 60 * i / C) * 1.02 ** (top - k));
  const pairAt = (i) => i % (K - 1), ratioOf = (i) => 1 - 5e-7 + 3.5e-6 * i / (C - 1);
  const pairs = await bothEngines(model, (i) => { const theta = Array.from({ length: K }, (_, k) => stable(i, k)), k = pairAt(i); theta[k + 1] = Math.fround(theta[k] * ratioOf(i)); return { theta }; });
  let parted = 0, gpuOnly = 0, lowest = Infinity, highest = 0;
  const outside = [], unmixed = [], moved = [];
  for (let i = 0; i < C; i++) {
    const k = pairAt(i), up = k * C + i, low = up + C, ratio = pairs.input[low] / pairs.input[up], cpu = pairs.cpu.mixes[i] > 0, gpu = pairs.gpu.mixes[i] > 0;
    if (cpu !== ratio > THRESHOLD) outside.push(`column ${i}, ratio 1 + ${(ratio - 1).toExponential(3)}: the CPU ${cpu ? 'merges' : 'keeps'} it`);
    if (cpu !== gpu) {
      parted++; if (gpu) gpuOnly++; lowest = Math.min(lowest, ratio); highest = Math.max(highest, ratio);
      if (!(ratio > BAND[0] && ratio <= BAND[1])) outside.push(`column ${i}, ratio 1 + ${(ratio - 1).toExponential(3)}, CPU ${cpu}, GPU ${gpu}`);
    } else if (cpu) {
      for (const x of [up, low]) if (!near(pairs.gpu.theta[x], pairs.cpu.theta[x], 2, pairs.cpu.theta[x]) || !near(pairs.gpu.q[x], pairs.cpu.q[x], 2, 5e-3) || !near(pairs.gpu.qc[x], pairs.cpu.qc[x], 2, 4e-5)) unmixed.push(`column ${i} layer ${Math.floor(x / C)}: θ ${pairs.gpu.theta[x]} against ${pairs.cpu.theta[x]}, q ${pairs.gpu.q[x]} against ${pairs.cpu.q[x]}`);
      if ((pairs.cpu.mixes[i] !== 1) || (pairs.gpu.mixes[i] !== 1)) unmixed.push(`column ${i}: ${pairs.cpu.mixes[i]} CPU and ${pairs.gpu.mixes[i]} GPU merges`);
    } else if (pairs.gpu.theta[up] !== pairs.input[up] || pairs.gpu.theta[low] !== pairs.input[low] || pairs.cpu.theta[up] !== pairs.input[up] || pairs.cpu.theta[low] !== pairs.input[low]) moved.push(`column ${i}`);
  }
  console.log(`${levels}: of ${C} pairs at ratios from 1 - 5e-7 to 1 + 3e-6, ${parted} decided apart (${gpuOnly} merged by the GPU alone), at ratios from 1 + ${parted ? (lowest - 1).toExponential(3) : '-'} to 1 + ${parted ? (highest - 1).toExponential(3) : '-'}, inside the band from 1 + ${(BAND[0] - 1).toExponential(3)} to 1 + ${(BAND[1] - 1).toExponential(3)}`);
  assert.deepEqual(outside, [], `${levels}: pairs decided apart outside the f32 band`);
  assert.deepEqual(unmixed, [], `${levels}: pairs both engines merged but mixed apart`);
  assert.deepEqual(moved, [], `${levels}: pairs neither engine merged but one moved`);
  const depthOf = (i) => 2 + (i % 7), cascades = await bothEngines(model, (i) => {
    const m = depthOf(i), top = K - m, theta = Array.from({ length: K }, (_, k) => (k < top ? stable(i, k, top) * 1.02 : Math.fround(stable(i, top, top) * (1 + 1e-4 * (K - 1 - k)))));
    theta[top] = Math.fround(stable(i, top, top) * (1 - 1e-3 * m));
    return { theta };
  });
  const apart = [];
  for (let i = 0; i < C; i++) {
    const m = depthOf(i);
    if (cascades.cpu.mixes[i] !== m - 1 || cascades.gpu.mixes[i] !== m - 1) { apart.push(`column ${i}: ${cascades.cpu.mixes[i]} CPU and ${cascades.gpu.mixes[i]} GPU merges of ${m - 1}`); continue; }
    for (let k = K - m; k < K; k++) { const x = k * C + i; if (!near(cascades.gpu.theta[x], cascades.cpu.theta[x], m, cascades.cpu.theta[x]) || !near(cascades.gpu.q[x], cascades.cpu.q[x], m, 5e-3) || !near(cascades.gpu.qc[x], cascades.cpu.qc[x], m, 4e-5)) apart.push(`column ${i} layer ${k}: θ ${cascades.gpu.theta[x]} against ${cascades.cpu.theta[x]}, q ${cascades.gpu.q[x]} against ${cascades.cpu.q[x]}, qc ${cascades.gpu.qc[x]} against ${cascades.cpu.qc[x]}`); }
  }
  assert.deepEqual(apart.slice(0, 10), [], `${levels}: ${apart.length} cascades merged or mixed apart`);
});
