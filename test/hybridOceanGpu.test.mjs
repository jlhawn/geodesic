import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createOcean as createCpuOcean } from '../js/ocean/layered.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { HEAT_REFERENCE, SALT_REFERENCE } from '../js/gpu/layeredOcean.gpu.js';
import { HYBRID_LAYERS, labelBoundaries } from '../js/ocean/column.module.js';
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
const TRACKING = 'the GPU hybrid ocean tracks the serial CPU hybrid ocean at N=16: the start, every phase of a step, and twenty-one steps';
test(TRACKING, { skip: !gpuAvailable && 'webgpu not installed' }, () => tracking(HYBRID));
test(`${TRACKING}, with mixedLayerRestratification`, { skip: !gpuAvailable && 'webgpu not installed' }, () => tracking(RESTRAT));
test(`${TRACKING}, with mstarScheme rh18`, { skip: !gpuAvailable && 'webgpu not installed' }, () => tracking(RH18_ENERGY));
test(`${TRACKING}, with mstarScheme rh18, langmuirScheme rescale and mixedLayerRestratification`, { skip: !gpuAvailable && 'webgpu not installed' }, () => tracking(RH18_ENERGY_ALL));

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
