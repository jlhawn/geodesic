import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createOcean as createCpuOcean } from '../js/ocean/layered.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { HEAT_REFERENCE, SALT_REFERENCE } from '../js/gpu/layeredOcean.gpu.js';
import { HYBRID_LAYERS, labelBoundaries, H_MIN } from '../js/ocean/column.module.js';
import { SEAWATER, seawaterDensity, salinityForDensity } from '../js/ocean/seawater.module.js';

/*
 * The GPU hybrid ocean (vertical: 'hybrid' in js/gpu/layeredOcean.gpu.js)
 * against the serial CPU hybrid ocean at N=16 on a sea with a continent,
 * phase by phase through a step and over twenty-one steps, and the float32
 * water, heat and salt it keeps through each phase of its step, the regrid
 * and the column physics above all. One GPU model at a time (CLAUDE.md).
 */
let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 300 : -4000));
const HYBRID = { vertical: 'hybrid', everySteps: 1, iceDrift: false };
const RESTRAT = { ...HYBRID, mixedLayerRestratification: { coefficient: 50 } };
const FIXED_ENERGY = { ...HYBRID, mstarScheme: 'fixed', langmuirScheme: 'factor' };
const RH18_ENERGY = { ...HYBRID, mstarScheme: 'rh18' };
const RH18_ENERGY_ALL = { ...RH18_ENERGY, langmuirScheme: 'rescale', mixedLayerRestratification: { coefficient: 50 } };
const N = 16, DT = 1350;

function stats(cpu, gpu, mask = null) {
  let maxDiff = 0, at = -1, sumSq = 0, sumRef = 0, n = 0;
  for (let x = 0; x < cpu.length; x++) {
    if (mask && !mask[x]) continue;
    const d = Math.abs(cpu[x] - gpu[x]);
    if (d > maxDiff) { maxDiff = d; at = x; }
    sumSq += d * d; sumRef += cpu[x] * cpu[x]; n++;
  }
  return { maxDiff, at, rmsRel: Math.sqrt(sumSq / Math.max(sumRef, 1e-300)), rms: Math.sqrt(sumSq / Math.max(n, 1)) };
}
function buildStress(mesh) {
  const stress = new Float64Array(mesh.nEdges);
  for (let e = 0; e < mesh.nEdges; e++) {
    const lon = Math.atan2(mesh.xEdge[3 * e + 1], mesh.xEdge[3 * e]);
    stress[e] = 0.08 * Math.sin(2 * mesh.latEdge[e]) * Math.cos(lon) + 0.01 * Math.sin(5 * lon);
  }
  return stress;
}
function scenario(n, options = HYBRID) {
  const cpuModel = createModel(new Grid(n), { topography });
  const init = initializeState(cpuModel, {});
  for (let a = 0; a < init.length; a++) cpuModel.state[a].set(init[a]);
  for (let i = 0; i < cpuModel.mesh.nCells; i++) if (cpuModel.geography.land[i]) cpuModel.state[6][i] = 0;
  const surfaceT = Float64Array.from(cpuModel.state[3]), ice = Float64Array.from(cpuModel.state[6]), stress = buildStress(cpuModel.mesh);
  const cpu = createCpuOcean(cpuModel.mesh, { geography: cpuModel.geography, ...options });
  return { cpuModel, mesh: cpuModel.mesh, cpu, surfaceT, ice, stress };
}
const cpuFields = (ocean, C) => {
  const T = Float64Array.from(ocean.Q, (q, n) => q / Math.max(1e-9, ocean.h[n])), S = Float64Array.from(ocean.W, (w, n) => w / Math.max(1e-9, ocean.h[n]));
  return { h: ocean.h, u: ocean.u, T, S, eta: ocean.eta, C };
};
const compare = (label, a, b, wet, limit) => {
  const results = {};
  for (const f of ['h', 'u', 'T', 'S', 'eta']) results[f] = stats(a[f], b[f], f === 'T' || f === 'S' ? wet : null);
  console.log(`${label}: ${Object.entries(results).map(([f, r]) => `${f} rms ${r.rmsRel.toExponential(2)} max ${r.maxDiff.toExponential(2)}`).join(', ')}`);
  for (const f of ['h', 'u', 'T', 'S', 'eta']) assert.ok(results[f].rmsRel < limit[f] ?? limit.any, `${label}: ${f} rms relative diff ${results[f].rmsRel} (largest ${results[f].maxDiff} at ${results[f].at})`);
  return results;
};

async function tracking(options) {
  const { mesh, cpu, surfaceT, ice, stress } = scenario(N, options), C = mesh.nCells, E = mesh.nEdges;
  const gpuModel = await createGpuModel(new Grid(N), { topography, ocean: options });
  const gpu = gpuModel.oceanEngine;
  assert.equal(gpu.layers, cpu.layers);
  cpu.initialize(Float64Array.from(surfaceT), ice);
  gpu.initialize(Float64Array.from(surfaceT), ice);
  const wet = Array.from(cpu.h, (v) => v > 1);
  const loose = { h: 2e-2, u: 2e-3, T: 2e-4, S: 2e-5, eta: 2e-3 };
  // The deep interfaces of a nearly uniform column are ill-conditioned, so the start's thicknesses part between f32 and f64 where the water does not.
  compare('the start', cpuFields(cpu, C), await gpu.download(), wet, { h: 2e-2, u: 1, T: 1e-4, S: 1e-5, eta: 1e-3 });
  const flux = new Float64Array(C);
  // From the CPU's state on both engines, one step, then a step phase by phase, both reading the same surface and stress.
  const sync = async () => { await gpu.upload({ ...cpu.serialize(), Q: cpu.Q, W: cpu.W, previousT0: cpu.T0, previousIceVolume: ice, fresh: new Float64Array(C), buoyancy: new Float64Array(C), flux: new Float64Array(C), capacity: cpu.capacity }, Float64Array.from(surfaceT), ice); };
  await sync();
  cpu.advance(Float64Array.from(surfaceT), ice, flux, stress, DT);
  await gpu.advance(Float64Array.from(surfaceT), ice, stress, DT);
  {
    const a = cpuFields(cpu, C), b = await gpu.download();
    for (const k of [0, 1, 2, 3, 4, 5, 8, 12, 20, 30, 40, 50, 57]) { const r = stats(a.u.subarray(k * E, (k + 1) * E), b.u.subarray(k * E, (k + 1) * E)); let m = 0; for (let e = 0; e < E; e++) m = Math.max(m, Math.abs(a.u[k * E + e])); console.log(`  u layer ${k}: rms diff ${r.rms.toExponential(2)} m/s, largest ${r.maxDiff.toExponential(2)} at edge ${r.at} (cpu ${a.u[k * E + r.at].toExponential(3)} gpu ${b.u[k * E + r.at].toExponential(3)}), largest |u| ${m.toExponential(2)}`); }
    compare('one step', a, b, wet, loose);
  }
  await sync();
  cpu.readSurface(Float64Array.from(surfaceT), ice); cpu.setStress(stress, ice);
  gpu.readSurface(Float64Array.from(surfaceT), ice); gpu.setStress(stress, ice);
  const phases = [
    ['dynamics', () => cpu.phases.dynamics(DT), () => gpu.phases.dynamics(DT)],
    ['rescale', () => cpu.phases.rescale(), () => gpu.phases.rescale()],
    ['eddy', () => cpu.phases.eddy(DT), () => gpu.phases.eddy(DT)],
    ['shift', () => cpu.phases.shift(), () => gpu.phases.shift()],
    ['regrid', () => cpu.phases.regrid(DT), () => gpu.phases.regrid(DT)],
    ['column', () => cpu.phases.column(DT), () => gpu.phases.column(DT)],
    ['salt and surface', () => { cpu.salt(); cpu.writeSurface(Float64Array.from(surfaceT), flux, DT); }, () => { gpu.salt(DT); gpu.writeSurface(DT); }],
  ];
  for (const [name, onCpu, onGpu] of phases) {
    onCpu(); onGpu();
    if (name === 'column' && options.mstarScheme === 'rh18') {
      const depths = stats(cpu.boundaryDepths, (await gpu.restartArrays()).boundaryDepths, wet);
      console.log(`boundary depths after the column phase: rms ${depths.rmsRel.toExponential(2)} max ${depths.maxDiff.toExponential(2)} m`);
      assert.ok(depths.rmsRel < 1e-2, `the boundary layer depths part by ${depths.rmsRel} rms relative (largest ${depths.maxDiff} m)`);
    }
    compare(`after ${name}`, cpuFields(cpu, C), await gpu.download(), wet, loose);
  }
  const state = await gpu.serialize();
  assert.equal(state.vertical, 'hybrid');
  await gpu.upload(state, Float64Array.from(surfaceT), ice);
  const again = await gpu.download();
  for (const f of ['h', 'u', 'T', 'S', 'eta']) assert.ok(stats(Float64Array.from(state[f]), again[f], f === 'T' || f === 'S' ? wet : null).rmsRel < 1e-5, `round trip ${f}`);
  for (let n = 0; n < 19; n++) {
    cpu.advance(Float64Array.from(surfaceT), ice, flux, stress, DT);
    await gpu.advance(Float64Array.from(surfaceT), ice, stress, DT);
  }
  const final = await gpu.download();
  compare('twenty-one steps', cpuFields(cpu, C), final, wet, { h: 1e-2, u: 2e-2, T: 1e-3, S: 1e-4, eta: 1e-2 });
  for (const f of ['h', 'u', 'T', 'S', 'eta']) for (const v of final[f]) assert.ok(Number.isFinite(v), `${f} has a non-finite value after 21 steps`);
  const mixed = stats(cpu.mixedDepths, final.h1), cpuDiag = cpu.diagnostics(), gpuDiag = await gpu.diagnostics();
  console.log(`mixed-layer depth rms relative diff ${mixed.rmsRel.toExponential(2)} after 21 steps; diagnostics cpu=`, cpuDiag, 'gpu=', gpuDiag);
  assert.ok(mixed.rmsRel < 2e-2, `mixed-layer depth rms relative diff ${mixed.rmsRel}`);
  for (const key of ['oceanUpperDepth', 'oceanHeat', 'oceanInteriorT', 'oceanThermoclineDepth', 'oceanSalinity']) {
    const c = cpuDiag[key], g = gpuDiag[key], rel = Math.abs(c - g) / Math.max(Math.abs(c), 1e-6);
    assert.ok(rel < 0.05, `${key} relative diff ${rel}: cpu=${c} gpu=${g}`);
  }
  for (const key of ['oceanSpeed', 'oceanSSH']) assert.ok(Math.abs(cpuDiag[key] - gpuDiag[key]) < 0.05 * Math.max(Math.abs(cpuDiag[key]), 1) + 1e-3, `${key}: cpu=${cpuDiag[key]} gpu=${gpuDiag[key]}`);
  assert.equal(gpuDiag.oceanDried, 0);
  gpuModel.destroy();
}
/*
 * An option's change to the state over a few steps on the GPU against its
 * change on the CPU, both engines from the CPU's start with a sheared
 * current in the upper layers, under the same surface and stress: it sees
 * an option that moves the state by less than the engines part. The GPU
 * misses 10 to 15 % of the change from halving shearDiffusivity or raising
 * lateralDiffusionSpeed tenfold here. The sheared start parts the engines by
 * a few 1e-3 K rms in the top four layers over six steps and by about
 * 5e-5 K below them, so fromLayer 4 measures an option whose change in the
 * top layers is of that size below them.
 */
async function optionChange(base, extra, steps = 6, fromLayer = 0) {
  const runs = [];
  for (const options of [base, { ...base, ...extra }]) {
    const { mesh, cpu, surfaceT, ice, stress } = scenario(N, options), C = mesh.nCells, E = mesh.nEdges, flux = new Float64Array(C);
    const gpuModel = await createGpuModel(new Grid(N), { topography, ocean: options }), gpu = gpuModel.oceanEngine;
    cpu.initialize(Float64Array.from(surfaceT), ice);
    for (let k = 0; k < cpu.layers; k++) for (let e = 0; e < E; e++) cpu.u[k * E + e] = 0.3 * Math.exp(-k / 4) * Math.cos(mesh.latEdge[e]) * Math.sin(3 * e);
    await gpu.upload({ ...cpu.serialize(), Q: cpu.Q, W: cpu.W, previousT0: cpu.T0, previousIceVolume: ice, fresh: new Float64Array(C), buoyancy: new Float64Array(C), flux: new Float64Array(C), capacity: cpu.capacity }, Float64Array.from(surfaceT), ice);
    for (let n = 0; n < steps; n++) { cpu.advance(Float64Array.from(surfaceT), ice, flux, stress, DT); await gpu.advance(Float64Array.from(surfaceT), ice, stress, DT); }
    runs.push({ cpu: cpuFields(cpu, C), gpu: await gpu.download(), wet: Array.from(cpu.h, (v) => v > 1), layers: cpu.layers });
    gpuModel.destroy();
  }
  const [a, b] = runs, out = {};
  for (const f of ['T', 'S', 'u']) {
    let change = 0, miss = 0;
    const length = a.cpu[f].length, first = fromLayer * length / a.layers;
    for (let n = first; n < length; n++) if (f === 'u' || a.wet[n]) { const dc = b.cpu[f][n] - a.cpu[f][n], dg = b.gpu[f][n] - a.gpu[f][n]; change += dc * dc; miss += (dg - dc) ** 2; }
    out[f] = { change: Math.sqrt(change / (length - first)), ratio: Math.sqrt(miss / Math.max(change, 1e-300)) };
  }
  console.log(`the option's change on the CPU${fromLayer ? ` from layer ${fromLayer} down` : ''}, rms, and the GPU's miss of it, relative: ${Object.entries(out).map(([f, r]) => `${f} ${r.change.toExponential(2)}, ${r.ratio.toExponential(2)}`).join('; ')}`);
  return out;
}
const OPTION_CHANGE = 'the GPU hybrid ocean changes the state under an option as the CPU does, from one start over six steps';
test(`${OPTION_CHANGE}: shearDiffusivityResolution, under a strengthened shear-scaled diffusion`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const change = await optionChange({ ...HYBRID, shearDiffusivity: 10, shearDiffusivityCap: 1e6 }, { shearDiffusivityResolution: { scale: 4, power: 2 } });
  for (const f of ['T', 'S']) assert.ok(change[f].ratio < 0.1, `${f}: the GPU misses ${change[f].ratio} of the option's change`);
});
/*
 * optionChange's measure a step at a time: before each step the GPU starts
 * from the CPU's state, so it misses the option's step from a common state
 * and not the option's change along each engine's own path. The six-step
 * change from isoneutralDiffusion moves by 20 to 40 % on the CPU itself when
 * its start's top four layers are perturbed by 3e-3 to 6e-3 K rms, what the
 * engines part there over six steps, so that option is measured this way.
 */
async function optionSteps(base, extra, steps = 6) {
  const runs = [];
  for (const options of [base, { ...base, ...extra }]) {
    const { mesh, cpu, surfaceT, ice, stress } = scenario(N, options), C = mesh.nCells, E = mesh.nEdges, flux = new Float64Array(C), run = [];
    const gpuModel = await createGpuModel(new Grid(N), { topography, ocean: options }), gpu = gpuModel.oceanEngine;
    cpu.initialize(Float64Array.from(surfaceT), ice);
    for (let k = 0; k < cpu.layers; k++) for (let e = 0; e < E; e++) cpu.u[k * E + e] = 0.3 * Math.exp(-k / 4) * Math.cos(mesh.latEdge[e]) * Math.sin(3 * e);
    for (let n = 0; n < steps; n++) {
      const start = { ...cpuFields(cpu, C), u: Float64Array.from(cpu.u) };
      await gpu.upload({ ...cpu.serialize(), Q: cpu.Q, W: cpu.W, previousT0: cpu.T0, previousIceVolume: ice, fresh: new Float64Array(C), buoyancy: new Float64Array(C), flux: new Float64Array(C), capacity: cpu.capacity, boundaryDepths: cpu.boundaryDepths, mixedDepths: cpu.mixedDepths }, Float64Array.from(surfaceT), ice);
      await gpu.advance(Float64Array.from(surfaceT), ice, stress, DT);
      const moved = await gpu.download();
      cpu.advance(Float64Array.from(surfaceT), ice, flux, stress, DT);
      const after = cpuFields(cpu, C), step = { wet: Array.from(cpu.h, (v) => v > 1) };
      for (const f of ['T', 'S', 'u']) step[f] = [Float64Array.from(after[f], (v, m) => v - start[f][m]), Float64Array.from(moved[f], (v, m) => v - start[f][m])];
      run.push(step);
    }
    runs.push(run);
    gpuModel.destroy();
  }
  const out = {};
  for (const f of ['T', 'S', 'u']) {
    let change = 0, miss = 0, count = 0;
    for (let n = 0; n < steps; n++) {
      const [cpuA, gpuA] = runs[0][n][f], b = runs[1][n], [cpuB, gpuB] = b[f];
      for (let m = 0; m < cpuA.length; m++) if (f === 'u' || b.wet[m]) { const dc = cpuB[m] - cpuA[m], dg = gpuB[m] - gpuA[m]; change += dc * dc; miss += (dg - dc) ** 2; count++; }
    }
    out[f] = { change: Math.sqrt(change / count), ratio: Math.sqrt(miss / Math.max(change, 1e-300)) };
  }
  console.log(`the option's step on the CPU, rms over ${steps} steps, and the GPU's miss of it from the CPU's state, relative: ${Object.entries(out).map(([f, r]) => `${f} ${r.change.toExponential(2)}, ${r.ratio.toExponential(2)}`).join('; ')}`);
  return out;
}
test('the GPU hybrid ocean takes each of six steps under an option as the CPU does, from the CPU state: isoneutralDiffusion', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const change = await optionSteps(HYBRID, { isoneutralDiffusion: true });
  for (const [f, bound] of [['T', 0.2], ['S', 0.25]]) assert.ok(change[f].ratio < bound, `${f}: the GPU misses ${change[f].ratio} of the CPU's step`);
});
test('under isoneutralDiffusion and 30 m boundary layers the κS² diffusivity the GPU hybrid ocean\'s column takes from a dynamics phase is the CPU\'s 1:2:2:1 mean of the four stages\', not the fourth stage\'s, and its triads move the tracers above 30 m as the CPU\'s do', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const options = { ...HYBRID, isoneutralDiffusion: true }, { mesh, cpu, surfaceT, ice, stress } = scenario(N, options), C = mesh.nCells, E = mesh.nEdges, zeros = new Float64Array(C), dt = 4 * DT;
  const off = scenario(N, HYBRID).cpu, sheared = (ocean) => {
    ocean.initialize(Float64Array.from(surfaceT), ice);
    for (let k = 0; k < ocean.layers; k++) for (let e = 0; e < E; e++) ocean.u[k * E + e] = 0.6 * Math.exp(-k / 4) * Math.cos(mesh.latEdge[e]) * Math.sin(3 * e);
    ocean.boundaryDepths.fill(30);
  };
  sheared(cpu); sheared(off);
  const start = { ...cpu.serialize(), Q: cpu.Q, W: cpu.W, previousT0: cpu.T0, previousIceVolume: ice, fresh: zeros, buoyancy: zeros, flux: zeros, capacity: cpu.capacity, boundaryDepths: cpu.boundaryDepths, mixedDepths: cpu.mixedDepths };
  const upper = Array.from(cpu.h, (hv, n) => { let z = 0; for (let k = 0; k < Math.floor(n / C); k++) z += cpu.h[k * C + n % C]; return z < 30 && hv > 1; }), moved = [];
  let kappa;
  for (const ocean of [options, HYBRID]) {
    const gpuModel = await createGpuModel(new Grid(N), { topography, ocean }), gpu = gpuModel.oceanEngine;
    await gpu.upload(start, Float64Array.from(surfaceT), ice);
    gpu.setStress(stress.map(() => 0), ice);
    gpu.phases.dynamics(dt);
    if (ocean === options) kappa = Float64Array.from(await gpu.readSection('ISOK'));
    moved.push(await gpu.download());
    gpuModel.destroy();
  }
  cpu.phases.dynamics(dt);
  off.phases.dynamics(dt);
  const mean = Float64Array.from(cpu.isoKappa), on = cpuFields(cpu, C), without = cpuFields(off, C), turned = {};
  cpu.tendency(cpu.trial, cpu.stages[0].map((a) => new Float64Array(a.length)));
  const fourth = Float64Array.from(cpu.isoKappa), toMean = stats(mean, kappa), apart = stats(mean, fourth);
  for (const f of ['T', 'S']) {
    let miss = 0, change = 0;
    for (let n = 0; n < upper.length; n++) if (upper[n]) { const c = on[f][n] - without[f][n], g = moved[0][f][n] - moved[1][f][n]; miss += (g - c) ** 2; change += c * c; }
    turned[f] = Math.sqrt(miss / change);
  }
  console.log(`κS² after a dynamics phase, rms relative to the CPU's: the GPU's from the CPU's 1:2:2:1 mean ${toMean.rmsRel.toExponential(2)}, the CPU's fourth stage's from it ${apart.rmsRel.toExponential(2)}; the GPU's miss of the triads' change above 30 m, T ${turned.T.toExponential(2)} and S ${turned.S.toExponential(2)}`);
  assert.ok(apart.rmsRel > 0 && toMean.rmsRel < 2e-3 && toMean.rmsRel < 0.5 * apart.rmsRel, `the GPU's κS² parts from the CPU's stage mean by ${toMean.rmsRel}, the fourth stage's from it by ${apart.rmsRel}`);
  assert.ok(turned.T < 2e-2 && turned.S < 2e-2, `the GPU misses the triads' change above 30 m by ${turned.T} in T and ${turned.S} in S`);
});
test('under isoneutralDiffusion the GPU hybrid ocean takes the CPU\'s κS² diffusivity on a shelf sea, where no column takes a slope from a shelf\'s thin layer at its boundary layer\'s base, nor from a base its edge\'s sill closes', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const shelf = syntheticTopography(90, 180, (lat, lon) => { const d = Math.acos(Math.max(-1, Math.min(1, Math.cos(lat) * Math.cos(lon)))) * 180 / Math.PI; return d < 7 ? 100 : -Math.min(4000, 40 * 1.6 ** ((d - 7) / 4)); });
  const cpuModel = createModel(new Grid(N), { topography: shelf }), init = initializeState(cpuModel, {});
  for (let a = 0; a < init.length; a++) cpuModel.state[a].set(init[a]);
  const { mesh, geography } = cpuModel, { latCell, lonCell, cellsOnEdge } = mesh, C = mesh.nCells, surfaceT = Float64Array.from(cpuModel.state[3]), ice = new Float64Array(C), zeros = new Float64Array(C);
  const edgeSill = Float64Array.from(geography.edgeSill), near = new Set();
  const closeness = (i) => Math.sin(latCell[i]) * Math.sin(-Math.PI / 3) + Math.cos(latCell[i]) * Math.cos(-Math.PI / 3) * Math.cos(lonCell[i] - Math.PI);
  let pocket = 0;
  for (let i = 0; i < C; i++) if (closeness(i) > closeness(pocket)) pocket = i;
  for (let e = 0; e < mesh.nEdges; e++) if (cellsOnEdge[2 * e] === pocket || cellsOnEdge[2 * e + 1] === pocket) { edgeSill[e] = 50; near.add(cellsOnEdge[2 * e]); near.add(cellsOnEdge[2 * e + 1]); }
  const options = { ...HYBRID, isoneutralDiffusion: true, edgeSill };
  const cpu = createCpuOcean(mesh, { geography, ...options }), F = cpu.table.fixed, thin = 2 * H_MIN;
  cpu.initialize(Float64Array.from(surfaceT), ice);
  cpu.boundaryDepths.fill(30);
  for (const i of near) cpu.boundaryDepths[i] = 64;
  const h = Float64Array.from(cpu.h), start = { ...cpu.serialize(), Q: cpu.Q, W: cpu.W, previousT0: cpu.T0, previousIceVolume: ice, fresh: zeros, buoyancy: zeros, flux: zeros, capacity: cpu.capacity, boundaryDepths: cpu.boundaryDepths, mixedDepths: cpu.mixedDepths };
  const gpuModel = await createGpuModel(new Grid(N), { topography: shelf, ocean: options }), gpu = gpuModel.oceanEngine;
  await gpu.upload(start, Float64Array.from(surfaceT), ice);
  gpu.setStress(new Float64Array(mesh.nEdges), ice);
  gpu.phases.dynamics(DT);
  const kappa = Float64Array.from(await gpu.readSection('ISOK'));
  gpuModel.destroy();
  cpu.phases.dynamics(DT);
  const all = [], beside = [];
  for (let n = 0; n < F * C; n++) if (cpu.isoKappa[n] > 0) all.push(n);
  for (let e = 0; e < mesh.nEdges; e++) {
    const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1];
    if (geography.land[a] || geography.land[b]) continue;
    for (const [i, j] of [[a, b], [b, a]]) {
      let base = 0, below = 0;
      while (base < F && below < 30) below += h[base++ * C + i];
      if (h[(base - 1) * C + i] > thin && h[(base - 1) * C + j] <= thin) for (let m = 1; m < base; m++) beside.push((m - 1) * C + i);
    }
  }
  const miss = (entries) => { let d = 0, s = 0; for (const n of entries) { d += (kappa[n] - cpu.isoKappa[n]) ** 2; s += cpu.isoKappa[n] ** 2; } return Math.sqrt(d / s); };
  const sill = [];
  let pocketKappa = 0;
  for (let m = 0; m < F; m++) pocketKappa = Math.max(pocketKappa, Math.abs(kappa[m * C + pocket]), cpu.isoKappa[m * C + pocket]);
  for (const i of near) if (i !== pocket) for (let m = 0; m < F; m++) if (cpu.isoKappa[m * C + i] > 0) sill.push(m * C + i);
  console.log(`κS² on a shelf sea under 30 m boundary layers, the GPU's rms miss of the CPU's: ${miss(all).toExponential(2)} over ${all.length} interfaces, ${miss(beside).toExponential(2)} over the ${beside.length} above the base of a column beside a shelf's thin layer; a column whose edges 50 m sills close under 64 m boundary layers takes κS² up to ${pocketKappa} on either engine, and the GPU misses its neighbors' by ${miss(sill).toExponential(2)} over ${sill.length}`);
  assert.ok(beside.length > 0 && sill.length > 0, 'no column beside a shelf or a sill');
  assert.equal(pocketKappa, 0);
  assert.ok(miss(all) < 1e-2 && miss(beside) < 1e-2 && miss(sill) < 1e-2, `the GPU misses the CPU's κS² by ${miss(all)}, beside the shelves by ${miss(beside)} and beside the sills by ${miss(sill)}`);
});
test(`${OPTION_CHANGE}: backgroundLatitude equatorial, under a strengthened background diffusivity`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const change = await optionChange({ ...HYBRID, backgroundDiffusivity: 1e-3 }, { backgroundLatitude: 'equatorial' });
  for (const f of ['T', 'S']) assert.ok(change[f].ratio < 0.2, `${f}: the GPU misses ${change[f].ratio} of the CPU's change`);
});
test(`${OPTION_CHANGE}: closureScale`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const change = await optionChange(HYBRID, { closureScale: 0.25 });
  assert.ok(change.u.ratio < 0.1, `u: the GPU misses ${change.u.ratio} of the option's change`);
});
test(`${OPTION_CHANGE}: interiorShear lmd94`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const change = await optionChange({ ...HYBRID, shearViscosity: 5e-3 }, { interiorShear: 'lmd94' }, 6, 4);
  for (const f of ['T', 'S', 'u']) assert.ok(change[f].ratio < 0.1, `${f}: the GPU misses ${change[f].ratio} of the CPU's change`);
});
test(`${OPTION_CHANGE}: interiorShear lmd94 with no shear term, whose change is PP81's ν_b/(1 + 5Ri) share of κ`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const change = await optionChange({ ...HYBRID, shearViscosity: 0, backgroundViscosity: 1e-3 }, { interiorShear: 'lmd94', lmd94Viscosity: 0 }, 6, 4);
  for (const f of ['T', 'S']) assert.ok(change[f].ratio < 0.1, `${f}: the GPU misses ${change[f].ratio} of the CPU's change`);
});
for (const eddyTracers of ['third', 'centered']) test(`${OPTION_CHANGE}: eddyTracers ${eddyTracers}, under a strengthened restratification`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const change = await optionChange({ ...HYBRID, mixedLayerRestratification: { coefficient: 1000 } }, { eddyTracers });
  for (const f of ['T', 'S']) assert.ok(change[f].ratio < 0.1, `${f}: the GPU misses ${change[f].ratio} of the CPU's change`);
});

const TRACKING = 'the GPU hybrid ocean tracks the serial CPU hybrid ocean at N=16: the start, every phase of a step, and twenty-one steps';
test(TRACKING, { skip: !gpuAvailable && 'webgpu not installed' }, () => tracking(HYBRID));
test(`${TRACKING}, with mixedLayerRestratification`, { skip: !gpuAvailable && 'webgpu not installed' }, () => tracking(RESTRAT));
test(`${TRACKING}, with mstarScheme rh18`, { skip: !gpuAvailable && 'webgpu not installed' }, () => tracking(RH18_ENERGY));
test(`${TRACKING}, with mstarScheme rh18, langmuirScheme rescale and mixedLayerRestratification`, { skip: !gpuAvailable && 'webgpu not installed' }, () => tracking(RH18_ENERGY_ALL));
test(`${TRACKING}, with backgroundLatitude equatorial`, { skip: !gpuAvailable && 'webgpu not installed' }, () => tracking({ ...HYBRID, backgroundLatitude: 'equatorial' }));
test(`${TRACKING}, with interiorShear lmd94`, { skip: !gpuAvailable && 'webgpu not installed' }, () => tracking({ ...HYBRID, interiorShear: 'lmd94' }));
test(`${TRACKING}, with eddyTracers third and mixedLayerRestratification taperLatitude 5`, { skip: !gpuAvailable && 'webgpu not installed' }, () => tracking({ ...HYBRID, eddyTracers: 'third', mixedLayerRestratification: { coefficient: 50, taperLatitude: 5 } }));
test(`${TRACKING}, with eddyTracers centered and mixedLayerRestratification`, { skip: !gpuAvailable && 'webgpu not installed' }, () => tracking({ ...RESTRAT, eddyTracers: 'centered' }));
test(`${TRACKING}, with shearDiffusivityResolution`, { skip: !gpuAvailable && 'webgpu not installed' }, () => tracking({ ...HYBRID, shearDiffusivityResolution: { scale: 4, power: 2 } }));
test(`${TRACKING}, with closureScale 0.5`, { skip: !gpuAvailable && 'webgpu not installed' }, () => tracking({ ...HYBRID, closureScale: 0.5 }));
test(`${TRACKING}, with isoneutralDiffusion`, { skip: !gpuAvailable && 'webgpu not installed' }, () => tracking({ ...HYBRID, isoneutralDiffusion: true }));

/*
 * What each phase of the GPU hybrid step changes of the stores' heat
 * (W/m² as a global mean over the step), salt (psu·m) and water (m),
 * summed over every layer and cell in double precision from float32
 * snapshots; the rescale's water changes by design and carries its heat and
 * salt at the columns' means, the surface read and write change the
 * surface layer's heat by design, so those rows measure only salt and water.
 * The salt is h·S, so a row's salt carries its water's rounding at 35 psu:
 * the regrid's thicknesses are differences of float32 target depths.
 */
const FLOORS = {
  dynamics: [0.05, 0.2, 3e-4, 1e-3, 3e-6, 2e-5],
  eddy: [0.02, 0.05, 1e-6, 5e-6, 1e-7, 2e-7],
  shift: [1e-9, 1e-9, 1e-12, 1e-12, 1e-15, 1e-15],
  regrid: [0.05, 0.2, 1e-3, 3e-3, 3e-5, 6e-5],
  column: [0.05, 0.2, 3e-4, 1e-3, 1e-15, 1e-15],
};
async function conservation(options) {
  const { mesh, surfaceT, ice, stress } = scenario(N, options), C = mesh.nCells, E = mesh.nEdges;
  const gpuModel = await createGpuModel(new Grid(N), { topography, ocean: options });
  const gpu = gpuModel.oceanEngine, { device } = gpuModel.gpu, L = gpu.layers, rhoCp = gpu.options.density * gpu.options.specificHeat;
  gpu.initialize(Float64Array.from(surfaceT), ice);
  for (let n = 0; n < 4; n++) await gpu.advance(Float64Array.from(surfaceT), ice, stress, DT);
  let area = 0;
  for (let i = 0; i < C; i++) if (gpu.cellOcean[i]) area += mesh.areaCell[i];
  const snapshot = async () => {
    await device.queue.onSubmittedWorkDone();
    const { h, q, w } = await gpu.readStores();
    return { h: Float64Array.from(h), q: Float64Array.from(q), w: Float64Array.from(w) };
  };
  const books = (before, after) => {
    let heat = 0, salt = 0, volume = 0, offHeat = 0, offSalt = 0;
    for (let i = 0; i < C; i++) {
      if (!gpu.cellOcean[i]) continue;
      let dh = 0, dq = 0, dw = 0, hs = 0, qs = 0, ws = 0;
      for (let k = 0; k < L; k++) { const n = k * C + i; dh += after.h[n] - before.h[n]; dq += after.q[n] - before.q[n]; dw += after.w[n] - before.w[n]; hs += before.h[n]; qs += before.q[n]; ws += before.w[n]; }
      const a = mesh.areaCell[i];
      heat += a * dq; salt += a * (dw + SALT_REFERENCE * dh); volume += a * dh;
      offHeat += a * (dq - qs / hs * dh); offSalt += a * (dw - ws / hs * dh);
    }
    return { heat: rhoCp * heat / (area * DT), salt: salt / area, volume: volume / area, offHeat: rhoCp * offHeat / (area * DT), offSalt: offSalt / area };
  };
  const phases = [
    ['readSurface', () => { gpu.readSurface(Float64Array.from(surfaceT), ice); gpu.setStress(stress, ice); }],
    ['dynamics', () => gpu.phases.dynamics(DT)],
    ['rescale', () => gpu.phases.rescale()],
    ['eddy', () => gpu.phases.eddy(DT)],
    ['shift', () => gpu.phases.shift()],
    ['regrid', () => gpu.phases.regrid(DT)],
    ['column', () => gpu.phases.column(DT)],
    ['salt and surface', () => { gpu.salt(DT); gpu.writeSurface(DT); }],
  ];
  const record = Object.fromEntries(phases.map(([name]) => [name, []]));
  let before = await snapshot();
  for (let n = 0; n < 8; n++) {
    for (const [name, run] of phases) {
      run();
      const after = await snapshot();
      record[name].push(books(before, after));
      before = after;
    }
  }
  const largest = (rows, pick) => Math.max(...rows.map((r) => Math.abs(pick(r)))), mean = (rows, pick) => rows.reduce((s, r) => s + pick(r), 0) / rows.length, show = (x) => x.toExponential(1);
  for (const [name, rows] of Object.entries(record)) {
    const heat = name === 'rescale' ? (r) => r.offHeat : (r) => r.heat, salt = name === 'rescale' ? (r) => r.offSalt : (r) => r.salt;
    console.log(`GPU hybrid ${name}: heat ${name === 'rescale' ? 'off the columns\' means' : 'changed'} by ${show(mean(rows, heat))} W/m² on average, ${show(largest(rows, heat))} at most; salt ${show(mean(rows, salt))}, ${show(largest(rows, salt))} psu·m; water ${show(mean(rows, (r) => r.volume))}, ${show(largest(rows, (r) => r.volume))} m`);
  }
  for (const [name, [meanHeat, stepHeat, meanSalt, stepSalt, meanVolume, stepVolume]] of Object.entries(FLOORS)) {
    const rows = record[name], heat = (r) => r.heat, salt = (r) => r.salt, volume = (r) => r.volume;
    assert.ok(Math.abs(mean(rows, heat)) < meanHeat && largest(rows, heat) < stepHeat, `${name}: heat ${mean(rows, heat)} W/m² on average, ${largest(rows, heat)} at most`);
    assert.ok(Math.abs(mean(rows, salt)) < meanSalt && largest(rows, salt) < stepSalt, `${name}: salt ${mean(rows, salt)} psu·m on average, ${largest(rows, salt)} at most`);
    assert.ok(Math.abs(mean(rows, volume)) < meanVolume && largest(rows, volume) < stepVolume, `${name}: water ${mean(rows, volume)} m on average, ${largest(rows, volume)} at most`);
  }
  const rescale = record.rescale;
  assert.ok(largest(rescale, (r) => r.offHeat) < 0.5 && largest(rescale, (r) => r.offSalt) < 1e-4, 'the rescale carries heat and salt at the columns\' means');
  const final = await gpu.download();
  for (const f of ['h', 'u', 'T', 'S', 'eta']) for (const v of final[f]) assert.ok(Number.isFinite(v), `${f} finite`);
  assert.equal(final.u.length, L * E);
  gpuModel.destroy();
}
const CONSERVATION = 'in single precision the GPU hybrid ocean keeps the heat, salt and water it stores through the regrid and the column physics to the stores\' rounding, and through the other phases of its step';
test(CONSERVATION, { skip: !gpuAvailable && 'webgpu not installed' }, () => conservation(HYBRID));
test(`${CONSERVATION}, with mixedLayerRestratification`, { skip: !gpuAvailable && 'webgpu not installed' }, () => conservation(RESTRAT));
test(`${CONSERVATION}, with mixedLayerRestratification and eddyTracers third`, { skip: !gpuAvailable && 'webgpu not installed' }, () => conservation({ ...RESTRAT, eddyTracers: 'third' }));
test(`${CONSERVATION}, with shearDiffusivityResolution`, { skip: !gpuAvailable && 'webgpu not installed' }, () => conservation({ ...HYBRID, shearDiffusivityResolution: { scale: 4, power: 2 } }));
test(`${CONSERVATION}, with isoneutralDiffusion`, { skip: !gpuAvailable && 'webgpu not installed' }, () => conservation({ ...HYBRID, isoneutralDiffusion: true }));

test('the GPU hybrid ocean keeps a brackish band\'s fixed-depth surface layers as the CPU does and tracks it over five steps', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const brackish = (lat) => (Math.abs(lat) < 0.35 ? 12 : 34.5);
  const { mesh, cpuModel, surfaceT, ice, stress } = scenario(N), C = mesh.nCells;
  const cpu = createCpuOcean(mesh, { geography: cpuModel.geography, ...HYBRID, salinityProfile: brackish });
  const gpuModel = await createGpuModel(new Grid(N), { topography, ocean: { ...HYBRID, salinityProfile: brackish } });
  const gpu = gpuModel.oceanEngine;
  const cpuT = Float64Array.from(surfaceT);
  cpu.initialize(cpuT, ice);
  gpu.initialize(Float64Array.from(surfaceT), ice);
  const start = await gpu.download();
  let band = 0;
  for (let i = 0; i < C; i++) {
    if (!cpu.cellOcean[i] || Math.abs(mesh.latCell[i]) >= 0.35) continue;
    band++;
    for (let k = 0; k < HYBRID_LAYERS.fixed; k++) {
      assert.ok(Math.abs(cpu.h[k * C + i] - HYBRID_LAYERS.minimum[k]) < 1e-9, `CPU surface layer ${k} of brackish column ${i}: ${cpu.h[k * C + i]} m`);
      assert.ok(Math.abs(start.h[k * C + i] - HYBRID_LAYERS.minimum[k]) < 1e-3, `GPU surface layer ${k} of brackish column ${i}: ${start.h[k * C + i]} m`);
    }
  }
  assert.ok(band > 0, 'the band holds brackish columns');
  const wet = Array.from(cpu.h, (v) => v > 1), flux = new Float64Array(C);
  for (let n = 0; n < 5; n++) {
    cpu.advance(cpuT, ice, flux, stress, DT);
    await gpu.advance(Float64Array.from(surfaceT), ice, stress, DT);
  }
  const after = await gpu.download(), fixed = HYBRID_LAYERS.fixed * C;
  // The deep classes of a uniform column swap their water between f32 and f64, so the thicknesses are held to parity in the fixed-depth layers only.
  compare('five steps with a brackish band', cpuFields(cpu, C), after, wet, { h: 1, u: 2e-2, T: 1e-3, S: 1e-4, eta: 1e-2 });
  const top = stats(cpu.h.subarray(0, fixed), after.h.subarray(0, fixed));
  assert.ok(top.rmsRel < 1e-3, `fixed-depth thicknesses rms relative diff ${top.rmsRel} (largest ${top.maxDiff} at ${top.at})`);
  const sst = stats(cpuT, after.T.subarray(0, C), cpu.cellOcean);
  assert.ok(sst.maxDiff < 0.2, `the surface temperature the engines hand back differs by ${sst.maxDiff} K at ${sst.at}`);
});

/*
 * Each edge's transport as the continuity, the barotropic solver and the
 * shift count it: the top layer's centered thickness, the smaller of the
 * two below, the column cut at the edge's sill.
 */
function edgeTransports(fields, ocean, mesh) {
  const C = mesh.nCells, E = mesh.nEdges, L = ocean.layers, out = new Float64Array(E);
  for (let e = 0; e < E; e++) {
    if (!ocean.edgeOcean[e]) continue;
    const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1], sill = Math.max(1e-9, ocean.sillEdge[e] + 0.5 * (fields.eta[a] + fields.eta[b]));
    let above = 0;
    for (let k = 0; k < L; k++) {
      const hk = k === 0 ? 0.5 * (fields.h[a] + fields.h[b]) : Math.min(fields.h[k * C + a], fields.h[k * C + b]);
      out[e] += Math.min(hk, Math.max(0, sill - above)) * fields.u[k * E + e]; above += hk;
    }
  }
  return out;
}

test('with edgeTransport kept the GPU hybrid ocean holds every edge\'s transport through its regrid and column physics and tracks the CPU at N=16', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { cpuModel, mesh, surfaceT, ice, stress } = scenario(N), C = mesh.nCells, E = mesh.nEdges, KEPT = { ...HYBRID, edgeTransport: 'kept' };
  const cpu = createCpuOcean(mesh, { geography: cpuModel.geography, ...KEPT });
  const gpuModel = await createGpuModel(new Grid(N), { topography, ocean: KEPT });
  const gpu = gpuModel.oceanEngine, flux = new Float64Array(C);
  cpu.initialize(Float64Array.from(surfaceT), ice);
  gpu.initialize(Float64Array.from(surfaceT), ice);
  const wet = Array.from(cpu.h, (v) => v > 1), loose = { h: 2e-2, u: 2e-3, T: 2e-4, S: 2e-5, eta: 2e-3 };
  for (let n = 0; n < 4; n++) cpu.advance(Float64Array.from(surfaceT), ice, flux, stress, DT);
  await gpu.upload({ ...cpu.serialize(), Q: cpu.Q, W: cpu.W, previousT0: cpu.T0, previousIceVolume: ice, fresh: new Float64Array(C), buoyancy: new Float64Array(C), flux: new Float64Array(C), capacity: cpu.capacity }, Float64Array.from(surfaceT), ice);
  cpu.readSurface(Float64Array.from(surfaceT), ice); cpu.setStress(stress, ice);
  gpu.readSurface(Float64Array.from(surfaceT), ice); gpu.setStress(stress, ice);
  cpu.phases.dynamics(DT); gpu.phases.dynamics(DT);
  cpu.phases.rescale(); gpu.phases.rescale();
  cpu.phases.eddy(DT); gpu.phases.eddy(DT);
  cpu.phases.shift(); gpu.phases.shift();
  let fields = await gpu.download(), before = edgeTransports(fields, cpu, mesh);
  let largest = 0;
  for (let e = 0; e < E; e++) largest = Math.max(largest, Math.abs(before[e]));
  for (const [name, onCpu, onGpu] of [['regrid', () => cpu.phases.regrid(DT), () => gpu.phases.regrid(DT)], ['column', () => cpu.phases.column(DT), () => gpu.phases.column(DT)]]) {
    onCpu(); onGpu();
    fields = await gpu.download();
    compare(`kept, after ${name}`, cpuFields(cpu, C), fields, wet, loose);
    const after = edgeTransports(fields, cpu, mesh);
    let worst = 0;
    for (let e = 0; e < E; e++) worst = Math.max(worst, Math.abs(after[e] - before[e]));
    console.log(`kept, ${name}: largest change of an edge's transport ${(worst / largest).toExponential(2)} of the largest transport ${largest.toFixed(2)} m²/s`);
    assert.ok(worst < 1e-5 * largest, `the GPU's ${name} moved an edge's transport by ${worst} m²/s`);
    before = after;
  }
  for (let n = 0; n < 16; n++) {
    cpu.advance(Float64Array.from(surfaceT), ice, flux, stress, DT);
    await gpu.advance(Float64Array.from(surfaceT), ice, stress, DT);
  }
  const final = await gpu.download();
  compare('kept, sixteen more steps', cpuFields(cpu, C), final, wet, { h: 1e-2, u: 2e-2, T: 1e-3, S: 1e-4, eta: 1e-2 });
  for (const f of ['h', 'u', 'T', 'S', 'eta']) for (const v of final[f]) assert.ok(Number.isFinite(v), `${f} has a non-finite value`);
  gpuModel.destroy();
});

test('at a six-hour ocean step, four times the coupled ocean\'s at N=16, the GPU hybrid ocean takes the CPU\'s parts of its dynamics, tracks it over two days at N=16 and keeps its currents slow', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const LONG = 6 * 3600;
  const { mesh, cpu, surfaceT, ice, stress } = scenario(N), C = mesh.nCells;
  const gpuModel = await createGpuModel(new Grid(N), { topography, ocean: HYBRID });
  const gpu = gpuModel.oceanEngine;
  cpu.initialize(Float64Array.from(surfaceT), ice);
  gpu.initialize(Float64Array.from(surfaceT), ice);
  const wet = Array.from(cpu.h, (v) => v > 1), flux = new Float64Array(C);
  await gpu.upload({ ...cpu.serialize(), Q: cpu.Q, W: cpu.W, previousT0: cpu.T0, previousIceVolume: ice, fresh: new Float64Array(C), buoyancy: new Float64Array(C), flux: new Float64Array(C), capacity: cpu.capacity }, Float64Array.from(surfaceT), ice);
  for (let n = 0; n < 8; n++) {
    cpu.advance(Float64Array.from(surfaceT), ice, flux, stress, LONG);
    await gpu.advance(Float64Array.from(surfaceT), ice, stress, LONG);
  }
  const final = await gpu.download();
  compare('eight six-hour steps', cpuFields(cpu, C), final, wet, { h: 1e-2, u: 2e-2, T: 1e-3, S: 1e-4, eta: 1e-2 });
  const fastest = (u) => u.reduce((m, v) => Math.max(m, Math.abs(v)), 0);
  assert.ok(fastest(cpu.u) < 1 && fastest(final.u) < 1, `currents reach ${fastest(cpu.u)} m/s on the CPU and ${fastest(final.u)} on the GPU`);
  gpuModel.destroy();
});

test('with rotationLimit rh18 the GPU hybrid ocean\'s boundary layer takes Reichl and Hallberg\'s rotational decay as the CPU\'s does and tracks it over sixteen steps at N=16', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { cpuModel, mesh, surfaceT, ice, stress } = scenario(N), C = mesh.nCells, RH18 = { ...FIXED_ENERGY, rotationLimit: 'rh18' };
  const cpu = createCpuOcean(mesh, { geography: cpuModel.geography, ...RH18 }), algebraic = createCpuOcean(mesh, { geography: cpuModel.geography, ...FIXED_ENERGY, rotationLimit: 'algebraic' });
  const gpuModel = await createGpuModel(new Grid(N), { topography, ocean: RH18 });
  const gpu = gpuModel.oceanEngine, flux = new Float64Array(C);
  for (const ocean of [cpu, algebraic, gpu]) ocean.initialize(Float64Array.from(surfaceT), ice);
  for (let n = 0; n < 16; n++) {
    cpu.advance(Float64Array.from(surfaceT), ice, flux, stress, DT);
    algebraic.advance(Float64Array.from(surfaceT), ice, flux, stress, DT);
    await gpu.advance(Float64Array.from(surfaceT), ice, stress, DT);
  }
  const final = await gpu.download(), { boundaryDepths } = await gpu.restartArrays(), wet = Array.from(cpu.h, (v) => v > 1);
  compare('rh18, sixteen steps', cpuFields(cpu, C), final, wet, { h: 1e-2, u: 2e-2, T: 1e-3, S: 1e-4, eta: 1e-2 });
  let engines = 0, limits = 0;
  for (let i = 0; i < C; i++) if (wet[i]) { engines += Math.abs(boundaryDepths[i] - cpu.boundaryDepths[i]); limits += Math.abs(algebraic.boundaryDepths[i] - cpu.boundaryDepths[i]); }
  console.log(`boundary depths: GPU against CPU ${engines.toFixed(2)} m summed, rh18 against algebraic on the CPU ${limits.toFixed(2)} m`);
  assert.ok(engines < 0.1 * limits, `the engines part by ${engines} m, the two limits by ${limits} m`);
  gpuModel.destroy();
});

test('with langmuirIce tanh the GPU hybrid ocean weights the Langmuir factor by the open water as the CPU does and tracks it over sixteen steps at N=16', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { cpuModel, mesh, surfaceT, ice, stress } = scenario(N), C = mesh.nCells, DAMPED = { ...HYBRID, langmuirIce: 'tanh' };
  const cpu = createCpuOcean(mesh, { geography: cpuModel.geography, ...DAMPED }), kept = createCpuOcean(mesh, { geography: cpuModel.geography, ...HYBRID, langmuirIce: 'none' });
  const gpuModel = await createGpuModel(new Grid(N), { topography, ocean: DAMPED });
  const gpu = gpuModel.oceanEngine, flux = new Float64Array(C);
  for (const ocean of [cpu, kept, gpu]) ocean.initialize(Float64Array.from(surfaceT), ice);
  for (let n = 0; n < 16; n++) {
    cpu.advance(Float64Array.from(surfaceT), ice, flux, stress, DT);
    kept.advance(Float64Array.from(surfaceT), ice, flux, stress, DT);
    await gpu.advance(Float64Array.from(surfaceT), ice, stress, DT);
  }
  const final = await gpu.download(), { boundaryDepths } = await gpu.restartArrays(), wet = Array.from(cpu.h, (v) => v > 1);
  compare('langmuirIce tanh, sixteen steps', cpuFields(cpu, C), final, wet, { h: 1e-2, u: 2e-2, T: 1e-3, S: 1e-4, eta: 1e-2 });
  let engines = 0, iced = 0, cells = 0;
  for (let i = 0; i < C; i++) if (wet[i]) {
    engines += Math.abs(boundaryDepths[i] - cpu.boundaryDepths[i]);
    if (ice[i] > 0) { iced += Math.abs(kept.boundaryDepths[i] - cpu.boundaryDepths[i]); cells++; }
  }
  console.log(`boundary depths: GPU against CPU ${engines.toFixed(2)} m summed; under ice (${cells} cells) tanh against none on the CPU ${iced.toFixed(2)} m`);
  assert.ok(cells > 0 && engines < 0.1 * iced, `the engines part by ${engines} m, the two options under ice by ${iced} m`);
  gpuModel.destroy();
});

test('with mixedLayerRestratification the GPU eddy phase takes the CPU\'s fluxes across a strong front at N=16 and keeps the water, heat and salt it stores', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { mesh, cpu, surfaceT, ice, stress } = scenario(N, RESTRAT), C = mesh.nCells;
  const gpuModel = await createGpuModel(new Grid(N), { topography, ocean: RESTRAT });
  const gpu = gpuModel.oceanEngine, { device } = gpuModel.gpu, L = gpu.layers, rhoCp = gpu.options.density * gpu.options.specificHeat;
  cpu.initialize(Float64Array.from(surfaceT), ice);
  const front = (i) => 288 + 4 * Math.tanh((mesh.latCell[i] - 0.8) / 0.12);
  for (let i = 0; i < C; i++) {
    let top = 0;
    for (let k = 0; k < L && top < 150; k++) {
      const n = k * C + i;
      if (cpu.h[n] > 0) { cpu.Q[n] = cpu.h[n] * front(i); cpu.W[n] = cpu.h[n] * 35; }
      top += cpu.h[n];
    }
  }
  cpu.boundaryDepths.fill(100);
  cpu.phases.regrid(Infinity);
  const zeros = new Float64Array(C);
  await gpu.upload({ ...cpu.serialize(), Q: cpu.Q, W: cpu.W, previousT0: cpu.T0, previousIceVolume: ice, fresh: zeros, buoyancy: zeros, flux: zeros, capacity: cpu.capacity, boundaryDepths: cpu.boundaryDepths, mixedDepths: cpu.mixedDepths }, Float64Array.from(surfaceT), ice);
  cpu.setStress(stress, ice); gpu.setStress(stress, ice);
  const wet = Array.from(cpu.h, (v) => v > 1), start = cpuFields(cpu, C), first = { h: Float64Array.from(start.h), T: Float64Array.from(start.T), S: Float64Array.from(start.S) };
  const g0 = await gpu.download(), begin = { h: stats(first.h, g0.h).rms, T: stats(first.T, g0.T, wet).rms, S: stats(first.S, g0.S, wet).rms };
  assert.ok(begin.h < 1e-3 && begin.T < 1e-4 && begin.S < 1e-5, `the upload parts from the CPU state: ${JSON.stringify(begin)}`);
  const area = mesh.areaCell, snapshot = async () => {
    await device.queue.onSubmittedWorkDone();
    const { h, q, w } = await gpu.readStores();
    return { h: Float64Array.from(h), q: Float64Array.from(q), w: Float64Array.from(w) };
  };
  let before = await snapshot(), worst = { volume: 0, heat: 0, salt: 0, column: 0 };
  for (let n = 0; n < 30; n++) {
    cpu.phases.eddy(DT); gpu.phases.eddy(DT);
    const after = await snapshot();
    let dh = 0, dq = 0, dw = 0, ocean = 0;
    for (let i = 0; i < C; i++) {
      if (!gpu.cellOcean[i]) continue;
      let column = 0;
      for (let k = 0; k < L; k++) { const m = k * C + i; column += after.h[m] - before.h[m]; dq += area[i] * (after.q[m] - before.q[m]); dw += area[i] * (after.w[m] + SALT_REFERENCE * after.h[m] - before.w[m] - SALT_REFERENCE * before.h[m]); }
      dh += area[i] * column; ocean += area[i]; worst.column = Math.max(worst.column, Math.abs(column));
    }
    worst.volume = Math.max(worst.volume, Math.abs(dh) / ocean); worst.heat = Math.max(worst.heat, Math.abs(rhoCp * dq / (ocean * DT))); worst.salt = Math.max(worst.salt, Math.abs(dw) / ocean);
    before = after;
  }
  const a = cpuFields(cpu, C), b = await gpu.download();
  const change = {}, diff = {};
  for (const f of ['h', 'T', 'S']) {
    change[f] = stats(first[f], a[f], f === 'h' ? null : wet).rms;
    diff[f] = stats(a[f], b[f], f === 'h' ? null : wet).rms;
  }
  console.log(`restratification start: GPU against CPU h ${begin.h.toExponential(1)} m, T ${begin.T.toExponential(1)} K, S ${begin.S.toExponential(1)} rms`);
  console.log(`restratification over 30 eddy steps: CPU changes h by ${change.h.toExponential(2)} m, T ${change.T.toExponential(2)} K rms; GPU differs by ${diff.h.toExponential(2)} m, ${diff.T.toExponential(2)} K; per application at most ${worst.volume.toExponential(1)} m water, ${worst.heat.toExponential(1)} W/m², ${worst.salt.toExponential(1)} psu·m, a column ${worst.column.toExponential(1)} m`);
  assert.ok(change.h > 1e-3 && change.T > 1e-3, 'the front moves water');
  for (const f of ['h', 'T']) assert.ok(diff[f] < 0.03 * change[f], `${f}: the engines part by ${diff[f]} against a change of ${change[f]}`);
  assert.ok(diff.S < 1e-5, `S: the engines part by ${diff.S} psu rms, the front's salinity change being ${change.S}`);
  const gpuMle = (await gpu.restartArrays()).mleDepths;
  let mle = 0;
  for (let i = 0; i < C; i++) if (wet[i]) mle = Math.max(mle, Math.abs(gpuMle[i] - cpu.mleDepths[i]));
  assert.ok(mle < 1e-3, `the filtered depths part by ${mle} m`);
  const [, floorHeat, , floorSalt, , floorVolume] = FLOORS.eddy;
  assert.ok(worst.volume < floorVolume && worst.heat < floorHeat && worst.salt < floorSalt && worst.column < 1e-4, `the books: ${JSON.stringify(worst)}`);
  gpuModel.destroy();
});

test('under recordEddyFluxes the GPU hybrid ocean keeps the volume and heat fluxes its eddy pass applies through the passes after it, and steps bit for bit as without it', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { mesh, cpu, surfaceT, ice, stress } = scenario(N), C = mesh.nCells, E = mesh.nEdges, L = cpu.layers, zeros = new Float64Array(C);
  cpu.initialize(Float64Array.from(surfaceT), ice);
  for (let k = 0; k < L; k++) for (let e = 0; e < E; e++) cpu.u[k * E + e] = 0.3 * Math.exp(-k / 4) * Math.cos(mesh.latEdge[e]) * Math.sin(3 * e);
  const start = { ...cpu.serialize(), Q: cpu.Q, W: cpu.W, previousT0: cpu.T0, previousIceVolume: ice, fresh: zeros, buoyancy: zeros, flux: zeros, capacity: cpu.capacity, boundaryDepths: cpu.boundaryDepths, mixedDepths: cpu.mixedDepths };
  const SECTIONS = ['EDDYVREC', 'EDDYQREC', 'EDDYV', 'EDDYQ'], runs = [];
  for (const recordEddyFluxes of [false, true]) {
    const gpuModel = await createGpuModel(new Grid(N), { topography, ocean: { ...HYBRID, recordEddyFluxes } }), gpu = gpuModel.oceanEngine;
    await gpu.upload(start, Float64Array.from(surfaceT), ice);
    await gpu.advance(Float64Array.from(surfaceT), ice, stress, DT);
    const run = { ...(await gpu.readStores()), eta: (await gpu.download()).eta, restart: await gpu.restartArrays() };
    if (recordEddyFluxes) {
      const read = () => Promise.all(SECTIONS.map((name) => gpu.readSection(name)));
      await gpu.upload(start, Float64Array.from(surfaceT), ice);
      gpu.setStress(stress, ice);
      run.before = await gpu.readStores();
      gpu.phases.eddy(DT);
      run.after = await gpu.readStores();
      run.eddy = await read();
      gpu.readSurface(Float64Array.from(surfaceT), ice);
      gpu.phases.shift(); gpu.phases.regrid(DT); gpu.phases.column(DT); gpu.phases.dynamics(DT);
      run.later = await read();
    }
    runs.push(run);
    gpuModel.destroy();
  }
  const [off, on] = runs, differing = {};
  const count = (a, b) => { let n = 0; for (let m = 0; m < a.length; m++) if (!Object.is(a[m], b[m])) n++; return n; };
  for (const f of ['h', 'u', 'q', 'w', 'eta']) differing[f] = count(off[f], on[f]);
  for (const [name, values] of Object.entries(off.restart)) if (ArrayBuffer.isView(values)) differing[name] = count(values, on.restart[name]);
  assert.ok(Object.values(differing).every((n) => n === 0), `recordEddyFluxes changes the step: ${JSON.stringify(differing)}`);
  const kept = Object.fromEntries(SECTIONS.map((name, s) => [name, count(on.eddy[s], on.later[s])]));
  console.log(`values the shift, regrid, column and dynamics passes after the eddy pass change, of ${L * E}: ${JSON.stringify(kept)}`);
  assert.ok(kept.EDDYVREC === 0 && kept.EDDYQREC === 0, `the passes after the eddy pass write over the record: ${JSON.stringify(kept)}`);
  const inflow = (record) => {
    const out = new Float64Array(L * C);
    for (let e = 0; e < E; e++) {
      const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1];
      for (let k = 0; k < L; k++) { out[k * C + b] += record[k * E + e]; out[k * C + a] -= record[k * E + e]; }
    }
    return out;
  };
  const applied = (f) => Float64Array.from(on.after[f], (v, m) => (v - on.before[f][m]) * mesh.areaCell[m % C] / DT);
  const volume = stats(applied('h'), inflow(on.later[0])), heat = stats(applied('q'), inflow(on.later[1]));
  console.log(`the record's convergence against the eddy pass's change of the stores, per layer and cell, rms relative: volume ${volume.rmsRel.toExponential(2)}, heat ${heat.rmsRel.toExponential(2)}`);
  assert.ok(stats(applied('h'), new Float64Array(L * C)).rms > 0, 'the eddy pass moves water');
  assert.ok(volume.rmsRel < 0.05, `the record's volume convergence parts from the water the eddy pass moved by ${volume.rmsRel} rms relative (largest ${volume.maxDiff} m³/s at ${volume.at})`);
  assert.ok(heat.rmsRel < 0.05, `the record's heat convergence parts from the heat the eddy pass moved by ${heat.rmsRel} rms relative (largest ${heat.maxDiff} K m³/s at ${heat.at})`);
});
test(`at N=${N} a thick deep layer lighter than its class boundary by less than float32 resolves at 1027 kg/m³ goes to the lighter class in the regrid on both engines`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { mesh, cpu, surfaceT, ice } = scenario(N, HYBRID), C = mesh.nCells;
  const gpuModel = await createGpuModel(new Grid(N), { topography, ocean: HYBRID });
  const gpu = gpuModel.oceanEngine;
  cpu.initialize(Float64Array.from(surfaceT), ice);
  const L = cpu.layers, boundaries = labelBoundaries(cpu.serialize().densities);
  let i = -1;
  for (let c = 0; c < C; c++) if (cpu.D[c] > 3000 && (i < 0 || cpu.D[c] > cpu.D[i])) i = c;
  let k = -1;
  for (let m = L - 1; m > L / 2; m--) if (cpu.h[m * C + i] > 200) { k = m; break; }
  assert.ok(i >= 0 && k > 0, `a deep column with a thick layer below mid-column (cell ${i}, layer ${k})`);
  const tied = (rho) => Math.fround(SEAWATER.rho0 + Math.fround(rho - SEAWATER.rho0)) === Math.fround(boundaries[k]);
  let deficit = 1e-5;
  while (deficit < 1e-4 && !tied(boundaries[k] - deficit)) deficit += 1e-6;
  assert.ok(deficit < 1e-4, 'a deficit that float32 rounds onto the boundary');
  const n = k * C + i, hk = cpu.h[n], T = cpu.Q[n] / hk, S = salinityForDensity(boundaries[k] - deficit, T);
  cpu.W[n] = hk * S;
  assert.ok(Math.abs(seawaterDensity(T, S) - (boundaries[k] - deficit)) < 1e-9);
  await gpu.upload({ ...cpu.serialize(), Q: cpu.Q, W: cpu.W, previousT0: cpu.T0, previousIceVolume: ice, fresh: new Float64Array(C), buoyancy: new Float64Array(C), flux: new Float64Array(C), capacity: cpu.capacity }, Float64Array.from(surfaceT), ice);
  const before = Array.from({ length: L }, (_, m) => cpu.h[m * C + i]);
  cpu.phases.regrid(DT); gpu.phases.regrid(DT);
  const after = await gpu.download(), column = (h) => Array.from({ length: L }, (_, m) => h[m * C + i]);
  const onCpu = column(cpu.h), onGpu = column(after.h);
  console.log(`cell ${i} layer ${k}, ${deficit.toExponential(1)} kg/m³ under the boundary: before ${before[k - 1].toFixed(1)} / ${before[k].toFixed(1)} m, CPU ${onCpu[k - 1].toFixed(1)} / ${onCpu[k].toFixed(1)} m, GPU ${onGpu[k - 1].toFixed(1)} / ${onGpu[k].toFixed(1)} m`);
  assert.ok(onCpu[k - 1] - before[k - 1] > 0.5 * before[k], `the CPU files the water in class ${k - 1}`);
  for (let m = 0; m < L; m++) assert.ok(Math.abs(onCpu[m] - onGpu[m]) < 1, `layer ${m}: CPU ${onCpu[m]} m, GPU ${onGpu[m]} m`);
  gpuModel.destroy();
});
