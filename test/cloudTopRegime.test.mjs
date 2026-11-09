import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { sigmaInterfaces, SIGMA_GRIDS } from '../js/dynamics/sigmaCore.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { saturationHumidity } from '../js/physics/moist.module.js';
import { CLOUD_TOP_DEFAULTS, REGIME } from '../js/physics/boundaryLayer.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const N = 16, LEVELS = ['bl48', 'bl60'].filter((name) => SIGMA_GRIDS.includes(name));
const THRESHOLD = CLOUD_TOP_DEFAULTS.threshold, F32_THRESHOLD = Math.fround(THRESHOLD), U = 2 ** -24, STEP = 2 ** -10;
const bits = new Float32Array(1), word = new Uint32Array(bits.buffer);
const ulps = (x, j) => { bits[0] = x; word[0] += j; return bits[0]; };
const WATERS = [0, 5e-7, ...[-3, -2, -1, 0, 1, 2, 3].map((j) => ulps(F32_THRESHOLD, j)), 2e-6, 1e-4].map(Math.fround);
const inWaterBand = (v) => v > Math.min(THRESHOLD, F32_THRESHOLD) && v <= Math.max(THRESHOLD, F32_THRESHOLD);

function build(levels) {
  const pair = createModel(new Grid(N), { ocean: false, levels: sigmaInterfaces(levels), boundaryLayer: { diagnoseFirst: false }, surface: { exchange: 'fixed' } });
  const { core, state, radiation } = pair, { K, C, exnerLayer, sigmaMid, kappa, geopotential, g } = core.diagnostics;
  const [pi, theta, u, surfaceT, q, qc] = state, longwave = radiation.longwave;
  const init = initializeState(pair, {});
  for (let a = 0; a < init.length; a++) state[a].set(init[a]);
  core.diagnose(pi, theta, q, qc);
  let seed = 11;
  const random = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; };
  const plan = [];
  for (let i = 0; i < C; i++) {
    const zb = geopotential[(K - 1) * C + i] / g, height = (k) => geopotential[k * C + i] / g - zb, family = ['water', 'top', 'exact', 'rounded'][i % 4], j = Math.floor(i / 4);
    const n = family === 'water' ? 1 : family === 'top' ? 2 : family === 'exact' ? 2 + (j % 3) : 3 + (j % 2);
    let low = K - 1;
    const wanted = 500 + 700 * random();
    while (low > 0 && height(low) < wanted) low--;
    const top = low - n + 1, surface = 285 + 12 * random();
    for (let k = 0; k < K; k++) {
      const x = k * C + i, p = pi[i] * sigmaMid[k];
      theta[x] = k > low ? surface : k >= top ? surface + 3 : Math.max(surface + 10 + 4e-3 * (height(k) - height(top)), init[1][x]);
      q[x] = (k > low ? 0.5 : k >= top ? 1 : 0.3) * saturationHumidity(theta[x] * exnerLayer[x], p);
      qc[x] = k >= top && k <= low ? 1e-4 : 0;
      longwave[x] = 0;
    }
    const terms = [], entry = { family, top, n, terms };
    if (family === 'water') { entry.water = WATERS[j % WATERS.length]; qc[top * C + i] = entry.water; terms.push(-60); }
    else if (family === 'top') { entry.water = WATERS[j % WATERS.length]; qc[top * C + i] = entry.water; terms.push(-40, -60); }
    else if (family === 'exact') {
      for (let k = 1; k < n; k++) terms.push((random() < 0.5 ? -1 : 1) * Math.floor(1000 + 59000 * random()));
      entry.offset = (j % 5) - 2;
      terms.push(entry.offset - terms.reduce((a, b) => a + b, 0));
      for (let k = 0; k < n; k++) terms[k] *= STEP;
    } else {
      for (let k = 1; k < n; k++) terms.push(Math.fround((random() < 0.5 ? -1 : 1) * (100 + 900 * random())));
      const before = terms.reduce((a, b) => a + b, 0), scale = terms.reduce((a, b) => a + Math.abs(b), 0) + Math.abs(before);
      terms.push(Math.fround((2 * random() - 1) * 3 * (n - 2) * U * scale - before));
    }
    terms.forEach((t, k) => { longwave[(top + k) * C + i] = t; });
    entry.sum = terms.reduce((a, b) => a + b, 0);
    entry.band = Math.max(0, n - 2) * U * terms.reduce((a, b) => a + Math.abs(b), 0);
    surfaceT[i] = theta[(K - 1) * C + i] * exnerLayer[(K - 1) * C + i] * Math.pow(sigmaMid[K - 1], -kappa) - 5 - 2 * random();
    plan.push(entry);
  }
  u.fill(0);
  for (const a of state) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  for (let x = 0; x < K * C; x++) longwave[x] = Math.fround(longwave[x]);
  core.diagnose(pi, theta, q, qc);
  return { pair, plan, K, C };
}

const expectedCloud = (entry) => (entry.family === 'water' ? entry.water > THRESHOLD : entry.family === 'top' ? true : entry.sum < 0);
const expectedCooling = (entry) => (entry.family === 'top' ? (entry.water > THRESHOLD ? 100 : 60) : -entry.sum);

for (const levels of LEVELS) {
  test(`on ${levels} the CPU puts a cloud-top layer where the cloud water passes 1e-6 and the run's longwave cools it, and the built columns are stable or decoupled as built`, () => {
    const { pair, plan, C } = build(levels), layer = pair.boundaryLayer;
    layer.diagnose(pair.state);
    const wrong = [];
    for (let i = 0; i < C; i++) {
      const entry = plan[i], cloud = expectedCloud(entry), regime = layer.regime[i];
      if (regime !== (cloud ? REGIME.DECOUPLED : REGIME.STABLE)) wrong.push(`column ${i} (${entry.family}, sum ${entry.sum}, water ${entry.water}): regime ${regime}`);
      else if (cloud && layer.cloudTopCooling[i] !== expectedCooling(entry)) wrong.push(`column ${i} (${entry.family}): cooling ${layer.cloudTopCooling[i]} against ${expectedCooling(entry)}`);
    }
    assert.deepEqual(wrong.slice(0, 10), [], `${levels}: ${wrong.length} columns`);
  });

  test(`on ${levels} the cloud-top decision, and with it the regime, is alike on both engines: cloud water at each f32 neighbor of the 1e-6 threshold and longwave runs summing exactly to zero or a step either side decide alike, and a rounded sum parts the engines only inside its f32 band`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
    const { createGpuCore } = await import('../js/gpu/core.gpu.js');
    const { pair, plan, K, C } = build(levels), { mesh, state, radiation, boundaryLayer: layer } = pair;
    const gpu = await createGpuCore(mesh, { levels: sigmaInterfaces(levels), physics: { diagnoseFirst: false, surfaceExchange: 'fixed' } });
    const { device, buffers, kernels } = gpu;
    gpu.upload(state);
    gpu.uploadPhysics();
    device.queue.writeBuffer(buffers.PH, 4 * gpu.layout.PH.LWH, Float32Array.from(radiation.longwave.subarray(0, K * C)));
    device.queue.writeBuffer(buffers.P, 0, Float32Array.from([900, 0, 1, 0, 0, 0, 0, 0]));
    const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
    pass.setPipeline(kernels.pblDiagnose);
    pass.setBindGroup(0, device.createBindGroup({ layout: kernels.pblDiagnose.getBindGroupLayout(0), entries: [buffers.MI, buffers.MF, buffers.LV, buffers.S, buffers.K1, buffers.D, buffers.P, buffers.PH].map((buffer, binding) => ({ binding, resource: { buffer } })) }));
    pass.dispatchWorkgroups(Math.ceil(C / 64));
    pass.end();
    device.queue.submit([encoder.finish()]);
    const ph = await gpu.downloadPhysics();
    layer.diagnose(state);
    const apart = [], counts = {};
    let parted = 0, gpuOnly = 0, inBand = 0;
    for (let i = 0; i < C; i++) {
      const entry = plan[i], cpu = layer.regime[i], gpuRegime = ph.REGIME[i], cpuCloud = cpu >= REGIME.DECOUPLED, gpuCloud = gpuRegime >= REGIME.DECOUPLED;
      const band = entry.family === 'rounded' ? Math.abs(entry.sum) <= entry.band : entry.family === 'exact' ? false : inWaterBand(entry.water);
      const count = (counts[entry.family] ??= { columns: 0, clouds: 0 });
      count.columns++; if (cpuCloud) count.clouds++; if (band) inBand++;
      if (cpuCloud !== gpuCloud) {
        parted++; if (gpuCloud) gpuOnly++;
        if (!band) apart.push(`column ${i} (${entry.family}, ${entry.n} layers, sum ${entry.sum}, band ${entry.band}, water ${entry.water}): CPU regime ${cpu}, GPU ${gpuRegime}`);
      } else if (cpu !== gpuRegime) apart.push(`column ${i} (${entry.family}): CPU regime ${cpu}, GPU ${gpuRegime}, both ${cpuCloud ? 'with' : 'without'} a cloud top`);
      else if (cpuCloud && !(Math.abs(ph.CTCOOL[i] - layer.cloudTopCooling[i]) <= (1 + 1e-6) * Math.max(0, entry.n - 1) * U * entry.terms.reduce((a, b) => a + Math.abs(b), 0))) apart.push(`column ${i} (${entry.family}): cooling ${ph.CTCOOL[i]} on the GPU, ${layer.cloudTopCooling[i]} on the CPU`);
    }
    console.log(`${levels}: ${Object.entries(counts).map(([family, { columns, clouds }]) => `${family} ${clouds} of ${columns} topped`).join(', ')}; ${inBand} columns inside their band, ${parted} decided apart (${gpuOnly} topped by the GPU alone)`);
    assert.deepEqual(apart.slice(0, 10), [], `${levels}: ${apart.length} columns decided or cooled apart outside the band`);
  });
}
