import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createParallelModel } from '../js/parallel.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { createSigmaCore, sigmaInterfaces, hybridWeights, dynamicsOptions, DIVERGENCE_DAMPING, DIVERGENCE_SURFACES } from '../js/dynamics/sigmaCore.module.js';
import { spongeRates, spongeSigmaFor, SPONGE, lidFrictionRates, lidFrictionFor } from '../js/dynamics/sponge.module.js';
import { divergence } from '../js/dynamics/operators.module.js';
import { pressureColumns, TEM_PRESSURES } from '../scripts/residualCirculation.mjs';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuCore } = gpuAvailable ? await import('../js/gpu/core.gpu.js') : {};

const N = 16, dt = 1350 * 16 / N, perDay = 86400 / dt, levels = sigmaInterfaces('bl60'), K = levels.length - 1;
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const model = createModel(new Grid(N), { topography, levels, physics: false });
const { mesh } = model, { sigmaMid } = model.core, C = mesh.nCells, E = mesh.nEdges, land = model.geography.land;
const start = initializeState(model, {});
const [startPi, startTheta] = start;
const dampingCore = (divergenceSurfaces, options = {}) => createSigmaCore(mesh, { levels, divergenceDamping: DIVERGENCE_DAMPING, divergenceSurfaces, splitClosure: true, ...options });
const rms = (a) => Math.sqrt(a.reduce((s, x) => s + x * x, 0) / a.length);
const largest = (a) => a.reduce((m, x) => Math.max(m, Math.abs(x)), 0);

// Easterly above 12 hPa and westerly below, the shear between about 5 and 30 hPa.
const SHEAR = { pressure: 1200, width: 0.45, mean: -8, half: 18 };
const qboWind = (lnp) => SHEAR.mean + SHEAR.half * Math.tanh((lnp - Math.log(SHEAR.pressure)) / SHEAR.width);
const qboCurvature = (lnp) => { const t = Math.tanh((lnp - Math.log(SHEAR.pressure)) / SHEAR.width); return -2 * SHEAR.half / SHEAR.width ** 2 * (1 - t * t) * t; };

function pressureWind(pi) {
  const u = new Float64Array(K * E);
  for (let e = 0; e < E; e++) {
    const piEdge = 0.5 * (pi[mesh.cellsOnEdge[2 * e]] + pi[mesh.cellsOnEdge[2 * e + 1]]);
    const lon = Math.atan2(mesh.xEdge[3 * e + 1], mesh.xEdge[3 * e]), east = -Math.sin(lon) * mesh.nEdge[3 * e] + Math.cos(lon) * mesh.nEdge[3 * e + 1];
    for (let k = 0; k < K; k++) u[k * E + e] = qboWind(Math.log(sigmaMid[k] * piEdge)) * east;
  }
  return u;
}

test('the hybrid weight is 1 where CAM L26 is pure pressure, falls through its hybrid levels and is 0 at the ground; the dynamics options take only divergenceSurfaces', () => {
  const w = hybridWeights(Float64Array.from([0.0005, 0.05, 0.0785, 0.1, 0.5, 0.9, 1]));
  assert.deepEqual(Array.from(w.subarray(0, 3)), [1, 1, 1]);
  assert.ok(w[3] > 0.7 && w[3] < 0.8, `w at σ 0.1: ${w[3]}`);
  assert.ok(w[4] > 0.05 && w[4] < 0.15 && w[5] < 0.02, `w at σ 0.5 and 0.9: ${w[4]}, ${w[5]}`);
  assert.ok(Math.abs(w[6]) < 1e-12);
  for (let j = 1; j < w.length; j++) assert.ok(w[j] <= w[j - 1]);
  assert.deepEqual(dynamicsOptions(), { divergenceSurfaces: 'sigma' });
  assert.deepEqual(dynamicsOptions({ divergenceSurfaces: 'hybrid' }), { divergenceSurfaces: 'hybrid' });
  assert.throws(() => dynamicsOptions({ divergenceSurface: 'hybrid' }), /divergenceSurface\b/);
  assert.throws(() => dynamicsOptions({ divergenceSurfaces: 'eta' }), /eta/);
  assert.throws(() => dampingCore('hybrid').phaseClosure([startPi, startTheta, pressureWind(startPi)], 0, 1, dt), /all of them at once/);
});

test(`an east wind on pressure surfaces over the terrain at N=${N} on bl60: damping the divergence along σ drives it at 5S–5N and 1–50 hPa as a vertical diffusion, ∝ ∂²U/∂(ln p)², and along the hybrid surfaces by a tenth of that or less`, () => {
  const pressures = TEM_PRESSURES.filter((p) => p >= 1 && p <= 50).map((p) => 100 * p);
  const tropics = [];
  for (let i = 0; i < C; i++) if (Math.abs(mesh.latCell[i]) <= 5 * Math.PI / 180) tropics.push(i);
  const area = tropics.reduce((s, i) => s + mesh.areaCell[i], 0);
  const tendency = {};
  for (const surfaces of DIVERGENCE_SURFACES) {
    const core = dampingCore(surfaces), state = [startPi, startTheta, pressureWind(startPi)];
    const before = pressureColumns(mesh, core, state, pressures).east;
    core.phaseClosure(state, 0, K, dt);
    const after = pressureColumns(mesh, core, state, pressures).east;
    tendency[surfaces] = pressures.map((p, j) => perDay * tropics.reduce((s, i) => s + mesh.areaCell[i] * (after[j * C + i] - before[j * C + i]), 0) / area);
  }
  const { sigma, hybrid } = tendency, curvature = pressures.map((p) => qboCurvature(Math.log(p)));
  const dot = (a, b) => a.reduce((s, x, j) => s + x * b[j], 0);
  const kappa = dot(sigma, curvature) / dot(curvature, curvature), correlation = dot(sigma, curvature) / Math.sqrt(dot(sigma, sigma) * dot(curvature, curvature));
  console.log(`zonal-mean tendency at 5S–5N (m/s/day), N=${N} bl60 over the terrain: ${pressures.map((p, j) => `${p / 100} hPa σ ${sigma[j].toExponential(2)} hybrid ${hybrid[j].toExponential(2)}`).join('; ')}; σ = ${kappa.toExponential(2)}/day × ∂²U/∂(ln p)² (correlation ${correlation.toFixed(4)}, K_z ${(7000 ** 2 * kappa / 86400).toFixed(3)} m²/s at H = 7 km), the hybrid's largest ${(largest(hybrid) / largest(sigma)).toFixed(3)} of the σ's`);
  assert.ok(kappa > 0 && correlation > 0.99, `σ against ∂²U/∂(ln p)²: κ ${kappa}, correlation ${correlation}`);
  assert.ok(largest(hybrid) <= 0.1 * largest(sigma), `hybrid ${largest(hybrid)} m/s/day against σ ${largest(sigma)}`);
});

test('a divergent wave packet over the open Pacific is damped along the hybrid surfaces as along σ', () => {
  const center = [0, -150 * Math.PI / 180], width = 1500e3;
  const distance = (i) => mesh.radius * Math.acos(Math.min(1, Math.sin(center[0]) * Math.sin(mesh.latCell[i]) + Math.cos(center[0]) * Math.cos(mesh.latCell[i]) * Math.cos(mesh.lonCell[i] - center[1])));
  const shape = Float64Array.from({ length: C }, (_, i) => Math.exp(-((distance(i) / width) ** 2)));
  const region = [];
  for (let i = 0; i < C; i++) if (distance(i) < 2 * width) region.push(i);
  assert.ok(region.every((i) => !land[i]), 'the packet is over the sea');
  const u0 = new Float64Array(K * E);
  for (let k = 0; k < K; k++) {
    const potential = 6e6 * (1 + 0.5 * Math.sin(2 * Math.log(sigmaMid[k])));
    for (let e = 0; e < E; e++) u0[k * E + e] = potential * (shape[mesh.cellsOnEdge[2 * e + 1]] - shape[mesh.cellsOnEdge[2 * e]]) / mesh.dcEdge[e];
  }
  const regionDivergence = (u) => {
    const d = new Float64Array(C);
    let sum = 0;
    for (let k = 0; k < K; k++) { divergence(mesh, u.subarray(k * E, (k + 1) * E), d); for (const i of region) sum += d[i] * d[i]; }
    return Math.sqrt(sum / (K * region.length));
  };
  const damped = {};
  for (const surfaces of DIVERGENCE_SURFACES) {
    const core = dampingCore(surfaces), u = Float64Array.from(u0);
    for (let n = 0; n < 20; n++) core.phaseClosure([startPi, startTheta, u], 0, K, dt);
    damped[surfaces] = u;
  }
  const before = regionDivergence(u0), decay = Object.fromEntries(DIVERGENCE_SURFACES.map((s) => [s, 1 - regionDivergence(damped[s]) / before]));
  const change = rms(damped.sigma.map((x, n) => x - u0[n])), apart = rms(damped.hybrid.map((x, n) => x - damped.sigma[n]));
  console.log(`a ${width / 1e3} km packet at 0°, 150°W, ${region.length} sea cells within ${2 * width / 1e3} km: 20 damping steps take ${(100 * decay.sigma).toFixed(2)} % of its rms divergence along σ and ${(100 * decay.hybrid).toFixed(2)} % along the hybrid surfaces; the two winds differ by ${(apart / change).toExponential(2)} of the damping's rms change`);
  assert.ok(decay.sigma > 0.02, `the σ damping takes ${decay.sigma} of the divergence`);
  assert.ok(Math.abs(decay.hybrid - decay.sigma) <= 0.01 * decay.sigma, `decay along σ ${decay.sigma}, along the hybrid surfaces ${decay.hybrid}`);
  assert.ok(apart <= 0.01 * change, `the winds differ by ${apart} m/s rms against the damping's ${change}`);
});

test(`repeated damping steps on noise over the terrain at N=${N} on bl60 converge along the hybrid surfaces as along σ, the step's largest eigenvalue far inside the explicit limit of 2`, () => {
  let seed = 12345;
  const random = () => { seed = (seed * 1664525 + 1013904223) >>> 0; return seed / 4294967296 - 0.5; };
  const noise = Float64Array.from({ length: K * E }, random), before = new Float64Array(K * E), after = new Float64Array(K * E);
  const norm = (a) => { let s = 0; for (let m = 0; m < a.length; m++) s += a[m] * a[m]; return Math.sqrt(s); };
  const found = {};
  for (const surfaces of DIVERGENCE_SURFACES) {
    const core = dampingCore(surfaces), step = (u) => core.phaseClosure([startPi, startTheta, u], 0, K, dt);
    const u = Float64Array.from(noise), changes = [];
    for (let n = 1; n <= 200; n++) {
      before.set(u);
      step(u);
      if (n !== 1 && n % 50 !== 0) continue;
      for (let m = 0; m < u.length; m++) before[m] = u[m] - before[m];
      changes.push(norm(before));
    }
    before.set(noise);
    let radius = 0;
    for (let n = 0; n < 40; n++) {
      after.set(before);
      step(after);
      for (let m = 0; m < after.length; m++) after[m] -= before[m];
      radius = norm(after) / norm(before);
      before.set(after);
    }
    found[surfaces] = { changes, radius };
  }
  const { sigma, hybrid } = found;
  console.log(`damping noise over the terrain, ‖u(n+1) − u(n)‖ at steps 1, 50, 100, 150, 200: σ ${sigma.changes.map((x) => x.toExponential(2)).join(', ')}; hybrid ${hybrid.changes.map((x) => x.toExponential(2)).join(', ')}; the step's largest eigenvalue σ ${sigma.radius.toFixed(4)}, hybrid ${hybrid.radius.toFixed(4)}`);
  for (let j = 1; j < hybrid.changes.length; j++) {
    assert.ok(hybrid.changes[j] < hybrid.changes[j - 1], `the hybrid steps grow: ${hybrid.changes}`);
    assert.ok(Math.abs(hybrid.changes[j] / sigma.changes[j] - 1) < 0.1, `hybrid ${hybrid.changes[j]} against σ ${sigma.changes[j]}`);
  }
  assert.ok(hybrid.radius < 0.5 && hybrid.radius < 1.1 * sigma.radius, `largest eigenvalue: hybrid ${hybrid.radius}, σ ${sigma.radius}`);
});

test('the hybrid closure takes every layer\'s ∇⁴ step, then the damping, then the sponge, and books before² − after² per edge and layer', () => {
  const nu4 = (model.core.spacing / Math.PI) ** 4 / (3 * 3600);
  const sponge = spongeRates(sigmaMid, spongeSigmaFor('bl60'), SPONGE.days), friction = lidFrictionRates(levels, lidFrictionFor('bl60'));
  const options = { nu4, spongeRates: sponge, spongeMeanRates: friction, surfaceGeopotential: model.surfaceGeopotential };
  const qbo = pressureWind(startPi), wind = Float64Array.from(start[2], (x, n) => x + qbo[n]);
  const run = (surfaces, extra = {}) => {
    const core = dampingCore(surfaces, { ...options, ...extra }), u = Float64Array.from(wind);
    core.phaseClosure([startPi, startTheta, u], 0, K, dt);
    return { u, dissipation: Float64Array.from(core.arrays.dissipation) };
  };
  const hybrid = run('hybrid');
  for (let n = 0; n < K * E; n++) if (hybrid.dissipation[n] !== wind[n] * wind[n] - hybrid.u[n] * hybrid.u[n]) assert.fail(`dissipation at ${n}: ${hybrid.dissipation[n]}`);
  const sigma = run('sigma'), undamped = run('sigma', { divergenceDamping: 0 });
  const correction = rms(hybrid.u.map((x, n) => x - sigma.u[n])), damping = rms(sigma.u.map((x, n) => x - undamped.u[n]));
  console.log(`one closure step on the geostrophic start plus the QBO wind: the σ damping moves the wind by ${damping.toExponential(2)} m/s rms, the hybrid surfaces by ${correction.toExponential(2)} m/s rms more`);
  assert.ok(correction > 1e-3 * damping && correction < damping, `the hybrid correction ${correction} m/s rms against the damping's ${damping}`);
});

test(`the worker-thread step with the hybrid damping reproduces the single-thread step bit for bit at N=${N} on bl60`, async () => {
  const options = { topography, levels, ocean: false, dynamics: { divergenceSurfaces: 'hybrid' } };
  const serial = createModel(new Grid(N), options);
  const parallel = await createParallelModel(new Grid(N), options, 3);
  for (let a = 0; a < 4; a++) { serial.state[a].set(start[a]); parallel.state[a].set(start[a]); }
  for (const m of [serial, parallel]) { m.seaIce.load(start[6]); m.land.initialize(); }
  try {
    assert.equal(serial.core.coupledClosure, true);
    for (let n = 0; n < 2; n++) { serial.step(dt); parallel.step(dt); }
    for (let a = 0; a < serial.state.length; a++) {
      const s = serial.state[a], p = parallel.state[a];
      for (let i = 0; i < s.length; i++) if (s[i] !== p[i]) assert.fail(`array ${a} differs at ${i}: ${s[i]} vs ${p[i]}`);
    }
    assert.deepEqual(parallel.core.arrays.dissipation, serial.core.arrays.dissipation);
  } finally {
    await parallel.close();
  }
});

test(`one closure step with the hybrid damping on the GPU matches the CPU at N=${N} on bl60 over the real terrain`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const nu4 = (model.core.spacing / Math.PI) ** 4 / (3 * 3600);
  const options = { levels, nu4, nu4Theta: nu4, divergenceDamping: DIVERGENCE_DAMPING, spongeRates: spongeRates(sigmaMid, spongeSigmaFor('bl60'), SPONGE.days), spongeMeanRates: lidFrictionRates(levels, lidFrictionFor('bl60')), surfaceGeopotential: model.surfaceGeopotential };
  const initial = start.map((a) => Float64Array.from(a));
  const qbo = pressureWind(startPi);
  for (let n = 0; n < K * E; n++) initial[2][n] += qbo[n];
  const runs = {};
  for (const surfaces of DIVERGENCE_SURFACES) {
    const cpu = createSigmaCore(mesh, { ...options, divergenceSurfaces: surfaces, splitClosure: true });
    const state = initial.map((a) => Float64Array.from(a));
    cpu.diagnose(state[0], state[1], state[4], state[5]);
    cpu.phaseClosure(state, 0, K, dt);
    const gpu = await createGpuCore(mesh, { ...options, divergenceSurfaces: surfaces, dragCoefficient: 0, topDragDays: 0, gravityWaves: false });
    gpu.upload(initial);
    await gpu.tendency();
    gpu.closure(dt);
    const got = await gpu.download(), { DISS } = await gpu.downloadDiagnostics();
    runs[surfaces] = { cpu: state, gpu: got, cpuDissipation: Float64Array.from(cpu.arrays.dissipation), gpuDissipation: Float64Array.from(DISS.subarray(0, K * E)), kernels: Object.keys(gpu.kernels), constants: gpu.preludeConstants };
  }
  assert.ok(!runs.sigma.constants.includes('L_HW') && !runs.sigma.kernels.some((name) => /Hybrid/.test(name)), 'σ compiles the kernels and the level constants it always has');
  assert.ok(runs.hybrid.constants.includes('L_HW') && runs.hybrid.kernels.includes('divergenceHybrid'));
  const start32 = Float64Array.from(initial[2], Math.fround);
  const worst = (cpu, gpu) => { let diff = 0, scale = 0; for (let n = 0; n < cpu.length; n++) { diff = Math.max(diff, Math.abs(gpu[n] - cpu[n])); scale = Math.max(scale, Math.abs(cpu[n])); } return { diff, scale, relative: diff / scale }; };
  const { cpu, gpu, cpuDissipation, gpuDissipation } = runs.hybrid;
  const step = worst(cpu[2].map((x, n) => x - initial[2][n]), gpu[2].map((x, n) => x - start32[n]));
  const correction = worst(cpu[2].map((x, n) => x - runs.sigma.cpu[2][n]), gpu[2].map((x, n) => x - runs.sigma.gpu[2][n]));
  const booked = worst(cpuDissipation, gpuDissipation);
  const theta = worst(cpu[1].map((x, n) => x - initial[1][n]), gpu[1].map((x, n) => x - Math.fround(initial[1][n])));
  console.log(`one hybrid closure step, GPU against CPU: the wind's change differs by ${step.diff.toExponential(2)} m/s of ${step.scale.toFixed(3)} (${step.relative.toExponential(2)}), the hybrid correction by ${correction.diff.toExponential(2)} of ${correction.scale.toFixed(4)} (${correction.relative.toExponential(2)}), the dissipation by ${booked.relative.toExponential(2)}, θ's change by ${theta.relative.toExponential(2)}`);
  assert.ok(step.relative < 1e-3, `wind change: ${step.diff} m/s against ${step.scale}`);
  assert.ok(correction.relative < 1e-2, `hybrid correction: ${correction.diff} m/s against ${correction.scale}`);
  assert.ok(booked.relative < 1e-2, `dissipation: ${booked.diff} against ${booked.scale}`);
  assert.ok(theta.relative < 1e-2, `θ change: ${theta.diff} K against ${theta.scale}`);
});
