import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createOcean as createCpuLayeredOcean } from '../js/ocean/layered.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 300 : -4000));
const OCEAN_OPTIONS = { everySteps: 1 };

function stats(cpu, gpu) {
  let maxDiff = 0, at = -1, sumSq = 0, sumRef = 0;
  for (let x = 0; x < cpu.length; x++) {
    const d = Math.abs(cpu[x] - gpu[x]);
    if (d > maxDiff) { maxDiff = d; at = x; }
    sumSq += d * d; sumRef += cpu[x] * cpu[x];
  }
  return { maxDiff, at, rmsRel: Math.sqrt(sumSq / Math.max(sumRef, 1e-300)) };
}

function buildStress(mesh) {
  const stress = new Float64Array(mesh.nEdges);
  for (let e = 0; e < mesh.nEdges; e++) {
    const lon = Math.atan2(mesh.xEdge[3 * e + 1], mesh.xEdge[3 * e]);
    stress[e] = 0.08 * Math.sin(2 * mesh.latEdge[e]) * Math.cos(lon) + 0.01 * Math.sin(5 * lon);
  }
  return stress;
}

function buildScenario(n) {
  const cpuModel = createModel(new Grid(n), { topography });
  const init = initializeState(cpuModel, {});
  for (let a = 0; a < init.length; a++) cpuModel.state[a].set(init[a]);
  for (let i = 0; i < cpuModel.mesh.nCells; i++) if (cpuModel.geography.land[i]) cpuModel.state[6][i] = 0;
  const surfaceT0 = Float64Array.from(cpuModel.state[3]);
  const ice = Float64Array.from(cpuModel.state[6]);
  const stress = buildStress(cpuModel.mesh);
  const cpuOcean = createCpuLayeredOcean(cpuModel.mesh, { geography: cpuModel.geography, ...OCEAN_OPTIONS });
  return { cpuModel, cpuOcean, surfaceT0, ice, stress };
}

test('one and twenty GPU ocean steps track the CPU layered ocean at N=8', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { cpuModel, cpuOcean, surfaceT0, ice, stress } = buildScenario(8);
  const gpuModel = await createGpuModel(new Grid(8), { topography, ocean: OCEAN_OPTIONS });
  const gpuOcean = gpuModel.oceanEngine;

  cpuOcean.initialize(surfaceT0, ice);
  gpuOcean.initialize(surfaceT0, ice);

  const dt = 1350;
  const oceanFluxScratch = new Float64Array(cpuModel.mesh.nCells);
  cpuOcean.advance(Float64Array.from(surfaceT0), ice, oceanFluxScratch, stress, dt);
  await gpuOcean.advance(Float64Array.from(surfaceT0), ice, stress, dt);

  const cpuState = cpuOcean.serialize();
  const gpuState = await gpuOcean.serialize();

  const fields = ['h', 'u', 'T', 'S', 'eta'];
  const results = {};
  // Tracers are compared only where the layer holds water: a token layer's
  // temperature is its label on one engine and the water it last held on the
  // other whenever their thicknesses straddle the token threshold.
  const wet = Array.from(cpuState.h, (v) => v > 1);
  for (const f of fields) results[f] = f === 'T' || f === 'S' ? stats(cpuState[f].filter((_, x) => wet[x]), gpuState[f].filter((_, x) => wet[x])) : stats(cpuState[f], gpuState[f]);
  console.log('one ocean step at N=8:', Object.entries(results).map(([f, r]) => `${f} rms ${r.rmsRel.toExponential(2)} max ${r.maxDiff.toExponential(2)}`).join(', '));
  for (const f of fields) assert.ok(results[f].rmsRel < 2e-3, `${f} rms relative diff ${results[f].rmsRel} at index ${results[f].at}`);

  // Round trip: upload(serialize()) reproduces the state.
  const roundTripSurface = Float64Array.from(surfaceT0);
  await gpuOcean.upload(gpuState, roundTripSurface, ice);
  const afterRoundTrip = await gpuOcean.serialize();
  for (const f of fields) {
    const r = stats(gpuState[f], afterRoundTrip[f]);
    assert.ok(r.rmsRel < 1e-4, `round trip ${f} relative diff ${r.rmsRel}`);
  }
  // Restore the pre-round-trip state so the 20-step run below continues from the real trajectory.
  await gpuOcean.upload(gpuState, roundTripSurface, ice);

  // 20 further steps: no NaN, and diagnostics stay close.
  for (let n = 0; n < 20; n++) {
    const surfaceTStep = Float64Array.from(surfaceT0);
    cpuOcean.advance(surfaceTStep, ice, oceanFluxScratch, stress, dt);
  }
  for (let n = 0; n < 20; n++) {
    const surfaceTStep = Float64Array.from(surfaceT0);
    await gpuOcean.advance(surfaceTStep, ice, stress, dt);
  }

  const gpuFinal = await gpuOcean.download();
  for (const f of ['h', 'u', 'T', 'S', 'eta']) {
    for (const v of gpuFinal[f]) assert.ok(Number.isFinite(v), `${f} has a non-finite value after 20 steps`);
  }

  const cpuDiag = cpuOcean.diagnostics();
  const gpuDiag = await gpuOcean.diagnostics();
  console.log('diagnostics after 21 total steps:', 'cpu=', cpuDiag, 'gpu=', gpuDiag);
  for (const key of ['oceanUpperDepth', 'oceanHeat', 'oceanThermoclineT', 'oceanThermoclineDepth', 'oceanSalinity']) {
    const c = cpuDiag[key], g = gpuDiag[key];
    const rel = Math.abs(c - g) / Math.max(Math.abs(c), 1e-6);
    assert.ok(rel < 0.05, `${key} relative diff ${rel}: cpu=${c} gpu=${g}`);
  }
  // Speed and SSH can be near zero this early; use a looser combined tolerance.
  for (const key of ['oceanSpeed', 'oceanSSH']) {
    const c = cpuDiag[key], g = gpuDiag[key];
    assert.ok(Math.abs(c - g) < 0.05 * Math.max(Math.abs(c), 1) + 1e-3, `${key} diff too large: cpu=${c} gpu=${g}`);
  }
});
