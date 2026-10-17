import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createOcean as createCpuOcean } from '../js/ocean/layered.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { HEAT_REFERENCE, SALT_REFERENCE } from '../js/gpu/layeredOcean.gpu.js';

/*
 * The GPU hybrid ocean (vertical: 'hybrid' in js/gpu/layeredOcean.gpu.js)
 * against the serial CPU hybrid ocean at N=8 on a sea with a continent,
 * phase by phase through a step and over twenty-one steps, and the float32
 * water, heat and salt it keeps through each phase of its step, the regrid
 * and the column physics above all. One GPU model at a time (CLAUDE.md).
 */
let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};
const { readRanges } = gpuAvailable ? await import('../js/gpu/device.module.js') : {};

const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 300 : -4000));
const HYBRID = { vertical: 'hybrid', everySteps: 1 };
const DT = 1350;

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
function scenario(n) {
  const cpuModel = createModel(new Grid(n), { topography });
  const init = initializeState(cpuModel, {});
  for (let a = 0; a < init.length; a++) cpuModel.state[a].set(init[a]);
  for (let i = 0; i < cpuModel.mesh.nCells; i++) if (cpuModel.geography.land[i]) cpuModel.state[6][i] = 0;
  const surfaceT = Float64Array.from(cpuModel.state[3]), ice = Float64Array.from(cpuModel.state[6]), stress = buildStress(cpuModel.mesh);
  const cpu = createCpuOcean(cpuModel.mesh, { geography: cpuModel.geography, ...HYBRID });
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

test('the GPU hybrid ocean tracks the serial CPU hybrid ocean at N=8: the start, every phase of a step, and twenty-one steps', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { mesh, cpu, surfaceT, ice, stress } = scenario(8), C = mesh.nCells, E = mesh.nEdges;
  const gpuModel = await createGpuModel(new Grid(8), { topography, ocean: HYBRID });
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
});

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
test('in single precision the GPU hybrid ocean keeps the heat, salt and water it stores through the regrid and the column physics to the stores\' rounding, and through the other phases of its step', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { mesh, surfaceT, ice, stress } = scenario(8), C = mesh.nCells, E = mesh.nEdges;
  const gpuModel = await createGpuModel(new Grid(8), { topography, ocean: HYBRID });
  const gpu = gpuModel.oceanEngine, { device } = gpuModel.gpu, { OS } = gpu.layout, L = gpu.layers, rhoCp = gpu.options.density * gpu.options.specificHeat;
  gpu.initialize(Float64Array.from(surfaceT), ice);
  for (let n = 0; n < 4; n++) await gpu.advance(Float64Array.from(surfaceT), ice, stress, DT);
  let area = 0;
  for (let i = 0; i < C; i++) if (gpu.cellOcean[i]) area += mesh.areaCell[i];
  const snapshot = async () => {
    await device.queue.onSubmittedWorkDone();
    const [h, q, w] = await readRanges(device, gpu.buffers.S, [{ offset: OS.OH, length: L * C }, { offset: OS.OQ, length: L * C }, { offset: OS.OW, length: L * C }]);
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
});
