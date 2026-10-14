import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createOcean as createCpuLayeredOcean, LAYER_DENSITIES, LAYER_SALINITIES, LAYER_BOTTOMS, langmuirFactor, transmitted, deposited, DAILY, DAY_SECONDS } from '../js/ocean/layered.module.js';
import { UNLISTED_OCEANS, columnMomentum } from './helpers/layered.mjs';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { labelTemperature, seawaterDensity, thermalExpansion } from '../js/ocean/seawater.module.js';
import { cellVector } from '../js/dynamics/operators.module.js';
import { FREEZING_POINT } from '../js/physics/ice.module.js';
import { SALT_REFERENCE } from '../js/gpu/layeredOcean.gpu.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};
const { readRanges } = gpuAvailable ? await import('../js/gpu/device.module.js') : {};

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
    const r = f === 'T' || f === 'S' ? stats(gpuState[f].filter((_, x) => wet[x]), afterRoundTrip[f].filter((_, x) => wet[x])) : stats(gpuState[f], afterRoundTrip[f]);
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
  const mixedDepth = stats(cpuOcean.serialize().h.slice(0, cpuModel.mesh.nCells), gpuFinal.h.slice(0, cpuModel.mesh.nCells));
  assert.ok(mixedDepth.rmsRel < 1e-2, `mixed-layer depth rms relative diff ${mixedDepth.rmsRel} after 21 steps`);
  for (const f of ['h', 'u', 'T', 'S', 'eta']) {
    for (const v of gpuFinal[f]) assert.ok(Number.isFinite(v), `${f} has a non-finite value after 20 steps`);
  }

  const cpuDiag = cpuOcean.diagnostics();
  const gpuDiag = await gpuOcean.diagnostics();
  console.log('diagnostics after 21 total steps:', 'cpu=', cpuDiag, 'gpu=', gpuDiag);
  for (const key of ['oceanUpperDepth', 'oceanHeat', 'oceanInteriorT', 'oceanThermoclineDepth', 'oceanSalinity']) {
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

test('under 10 m floors, with the wind over a 50 m momentumDepth and the mixed layer carrying its momentum, one and twenty-one GPU ocean steps track the CPU layered ocean at N=8', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const thin = { ...OCEAN_OPTIONS, mixedDepth: 10, minimumThickness: 10, shallowestMixedDepth: 10, momentumDepth: 50, carryMomentum: true };
  const { cpuModel, surfaceT0, ice, stress } = buildScenario(8);
  const mesh = cpuModel.mesh, C = mesh.nCells, dt = 1350;
  const gpuModel = await createGpuModel(new Grid(8), { topography, ocean: thin });
  const gpuOcean = gpuModel.oceanEngine;
  const cpuOcean = (options) => { const ocean = createCpuLayeredOcean(mesh, { geography: cpuModel.geography, ...options }); ocean.initialize(surfaceT0, ice); return ocean; };
  const rms = (a, b, mask = null) => { let d = 0, r = 0; for (let x = 0; x < a.length; x++) if (!mask || mask[x]) { d += (a[x] - b[x]) ** 2; r += a[x] * a[x]; } return Math.sqrt(d / Math.max(r, 1e-300)); };
  const cpu = cpuOcean(thin), today = cpuOcean({ ...thin, momentumDepth: 10, carryMomentum: false, shearEntrainment: 0 }), flux = new Float64Array(C);
  gpuOcean.initialize(surfaceT0, ice);
  for (const steps of [1, 20]) {
    for (let n = 0; n < steps; n++) {
      cpu.advance(Float64Array.from(surfaceT0), ice, flux, stress, dt);
      today.advance(Float64Array.from(surfaceT0), ice, flux, stress, dt);
      await gpuOcean.advance(Float64Array.from(surfaceT0), ice, stress, dt);
    }
    const cpuState = cpu.serialize(), todayState = today.serialize(), gpuState = await gpuOcean.serialize();
    const wet = Array.from(cpuState.h, (v) => v > 1);
    let shallow = 0;
    for (let i = 0; i < C; i++) if (cpu.cellOcean[i] && cpuState.h[i] < 20) shallow++;
    const lines = [];
    for (const f of ['h', 'u', 'T', 'S', 'eta']) {
      const mask = f === 'T' || f === 'S' ? wet : null, gap = rms(cpuState[f], gpuState[f], mask), effect = rms(cpuState[f], todayState[f], mask);
      lines.push(`${f} ${gap.toExponential(2)} against ${effect.toExponential(2)}`);
      assert.ok(gap < 2e-3, `${f} rms relative diff ${gap} after ${steps === 1 ? 1 : 21} steps`);
      if (f === 'u') assert.ok(effect > 5 * gap, `u: the momentum depth and its carriage changed the CPU ocean by ${effect}, no more than the engines differ (${gap})`);
    }
    console.log(`${steps === 1 ? 1 : 21} ocean steps under 10 m floors at N=8 (${shallow} sea cells with the mixed layer under 20 m), the engines' rms relative difference against the effect of the stress over momentumDepth 50 and the carriage: ${lines.join(', ')}`);
    assert.ok(shallow > 100, `${shallow} shallow mixed layers`);
  }
  for (const v of (await gpuOcean.download()).u) assert.ok(Number.isFinite(v));
  assert.equal((await gpuOcean.diagnostics()).oceanLimited, 0);
  gpuModel.destroy();
});

test('the GPU ocean keeps its water volume over 400 wind-driven steps at N=8', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { cpuModel, surfaceT0, ice, stress } = buildScenario(8);
  const gpuModel = await createGpuModel(new Grid(8), { topography, ocean: OCEAN_OPTIONS });
  const gpuOcean = gpuModel.oceanEngine, mesh = cpuModel.mesh, C = mesh.nCells;
  gpuOcean.initialize(surfaceT0, ice);
  async function meanColumn() {
    const { h } = await gpuOcean.download();
    let volume = 0, area = 0;
    for (let i = 0; i < C; i++) {
      if (!gpuOcean.cellOcean[i]) continue;
      let column = 0;
      for (let k = 0; k < gpuOcean.layers; k++) column += h[k * C + i];
      volume += mesh.areaCell[i] * column; area += mesh.areaCell[i];
    }
    return volume / area;
  }
  const before = await meanColumn();
  for (let n = 0; n < 400; n++) await gpuOcean.advance(Float64Array.from(surfaceT0), ice, stress, 1350);
  const drift = (await meanColumn()) - before;
  assert.ok(Math.abs(drift) < 2e-3, `mean water column changed by ${drift} m over 400 steps`);
});

test('the GPU eddy transport tracks the CPU\'s at N=8, alone and through twenty-one ocean steps', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const strong = { ...OCEAN_OPTIONS, eddyDiffusivity: 1e6 };
  const { cpuModel, surfaceT0, ice, stress } = buildScenario(8);
  const mesh = cpuModel.mesh, C = mesh.nCells, dt = 1350;
  const gpuModel = await createGpuModel(new Grid(8), { topography, ocean: strong });
  const gpuOcean = gpuModel.oceanEngine;
  const cpuOcean = (options) => { const ocean = createCpuLayeredOcean(mesh, { geography: cpuModel.geography, ...options }); ocean.initialize(surfaceT0, ice); return ocean; };
  const rms = (a, b, mask = null) => { let d = 0, r = 0; for (let x = 0; x < a.length; x++) if (!mask || mask[x]) { d += (a[x] - b[x]) ** 2; r += a[x] * a[x]; } return Math.sqrt(d / Math.max(r, 1e-300)); };

  const alone = cpuOcean(strong);
  await gpuOcean.upload(alone.serialize(), Float64Array.from(surfaceT0), ice);
  const start = Float64Array.from(alone.h), gpuStart = (await gpuOcean.download()).h;
  for (let n = 0; n < 20; n++) { alone.eddyTransport(dt); gpuOcean.eddyTransport(dt); }
  const cpuAlone = alone.serialize(), gpuAlone = await gpuOcean.download();
  const cpuChange = Float64Array.from(cpuAlone.h, (v, x) => v - start[x]), gpuChange = Float64Array.from(gpuAlone.h, (v, x) => v - gpuStart[x]);
  const wetAlone = Array.from(cpuAlone.h, (v) => v > 1);
  const change = rms(cpuChange, gpuChange), largest = cpuChange.reduce((m, v) => Math.max(m, Math.abs(v)), 0);
  console.log(`twenty eddy transports alone at N=8: layers move by up to ${largest.toFixed(1)} m, the GPU's change differing by ${change.toExponential(2)} rms relative, T ${rms(cpuAlone.T, gpuAlone.T, wetAlone).toExponential(2)}, S ${rms(cpuAlone.S, gpuAlone.S, wetAlone).toExponential(2)}`);
  assert.ok(largest > 5 && change < 1e-3, `the eddy transport moved layers by up to ${largest} m, the engines' changes differing by ${change}`);
  for (const f of ['T', 'S']) assert.ok(rms(cpuAlone[f], gpuAlone[f], wetAlone) < 2e-3, f);

  const cpu = cpuOcean(strong), still = cpuOcean({ ...OCEAN_OPTIONS, eddyDiffusivity: 0 }), flux = new Float64Array(C);
  gpuOcean.initialize(surfaceT0, ice);
  for (let n = 0; n < 21; n++) {
    cpu.advance(Float64Array.from(surfaceT0), ice, flux, stress, dt);
    still.advance(Float64Array.from(surfaceT0), ice, flux, stress, dt);
    await gpuOcean.advance(Float64Array.from(surfaceT0), ice, stress, dt);
  }
  const cpuState = cpu.serialize(), stillState = still.serialize(), gpuState = await gpuOcean.serialize();
  const wet = Array.from(cpuState.h, (v) => v > 1);
  const lines = [];
  for (const f of ['h', 'u', 'T', 'S', 'eta']) {
    const mask = f === 'T' || f === 'S' ? wet : null, gap = rms(cpuState[f], gpuState[f], mask), effect = rms(cpuState[f], stillState[f], mask);
    lines.push(`${f} ${gap.toExponential(2)} against ${effect.toExponential(2)}`);
    assert.ok(gap < 2e-3, `${f} rms relative diff ${gap} after 21 steps`);
    assert.ok(effect > 5 * gap, `${f}: the eddy transport changed the CPU ocean by ${effect}, no more than the engines differ (${gap})`);
  }
  console.log(`21 ocean steps with eddyDiffusivity 1e6 at N=8, the engines' rms relative difference against the eddy transport's effect: ${lines.join(', ')}`);
});

const THIN_FLOORS = { ...OCEAN_OPTIONS, mixedDepth: 10, minimumThickness: 10, shallowestMixedDepth: 10, momentumDepth: 50, carryMomentum: true };

test('under 10 m floors the GPU Ekman block carries the CPU\'s transport in the donor\'s top block, alone and through twenty-one ocean steps at N=8', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const ekman = { ...THIN_FLOORS, ekmanBlock: true };
  const { cpuModel, surfaceT0, ice, stress: weak } = buildScenario(8);
  const stress = Float64Array.from(weak, (x) => 10 * x);
  const mesh = cpuModel.mesh, C = mesh.nCells, E = mesh.nEdges, dt = 1350;
  const gpuModel = await createGpuModel(new Grid(8), { topography, ocean: ekman });
  const gpuOcean = gpuModel.oceanEngine;
  const cpuOcean = (options) => { const ocean = createCpuLayeredOcean(mesh, { geography: cpuModel.geography, ...options }); ocean.initialize(surfaceT0, ice); return ocean; };
  const rms = (a, b, mask = null) => { let d = 0, r = 0; for (let x = 0; x < a.length; x++) if (!mask || mask[x]) { d += (a[x] - b[x]) ** 2; r += a[x] * a[x]; } return Math.sqrt(d / Math.max(r, 1e-300)); };
  const { OD } = gpuOcean.layout;
  const edgeTransport = async () => (await readRanges(gpuModel.gpu.device, gpuOcean.buffers.OD, [{ offset: OD.UEK, length: E }]))[0];

  const alone = cpuOcean(ekman), flux = new Float64Array(C);
  gpuOcean.initialize(surfaceT0, ice);
  alone.advance(Float64Array.from(surfaceT0), ice, flux, stress, dt);
  await gpuOcean.advance(Float64Array.from(surfaceT0), ice, stress, dt);
  const cpuStart = alone.serialize(), gpuStart = await gpuOcean.download();
  // One long step: on the GPU a class's flux under half the coarser store's spacing moves nothing.
  alone.ekmanTransport(20 * dt); gpuOcean.phases.ekman(20 * dt);
  const cpuAlone = alone.serialize(), gpuAlone = await gpuOcean.download();
  const cpuChange = Float64Array.from(cpuAlone.h, (v, x) => v - cpuStart.h[x]), gpuChange = Float64Array.from(gpuAlone.h, (v, x) => v - gpuStart.h[x]);
  const wetAlone = Array.from(cpuAlone.h, (v) => v > 1);
  const change = rms(cpuChange, gpuChange), largest = cpuChange.reduce((m, v) => Math.max(m, Math.abs(v)), 0);
  const transport = await edgeTransport();
  let strongest = 0, transportGap = 0, transportNorm = 0;
  for (let e = 0; e < E; e++) { strongest = Math.max(strongest, Math.abs(alone.ekmanFlux[e])); transportGap += (alone.ekmanFlux[e] - transport[e]) ** 2; transportNorm += alone.ekmanFlux[e] ** 2; }
  console.log(`an Ekman block of twenty steps alone at N=8 under 10 m floors: layers move by up to ${largest.toFixed(3)} m, the GPU's change differing by ${change.toExponential(2)} rms relative, U_E up to ${strongest.toFixed(2)} m²/s and the GPU's differing by ${Math.sqrt(transportGap / transportNorm).toExponential(2)}, T ${rms(cpuAlone.T, gpuAlone.T, wetAlone).toExponential(2)}, S ${rms(cpuAlone.S, gpuAlone.S, wetAlone).toExponential(2)}`);
  assert.ok(largest > 0.02 && change < 5e-3, `the Ekman block moved layers by up to ${largest} m, the engines' changes differing by ${change}`);
  assert.ok(Math.sqrt(transportGap / transportNorm) < 2e-3, `U_E differs by ${Math.sqrt(transportGap / transportNorm)}`);
  for (const f of ['T', 'S']) assert.ok(rms(cpuAlone[f], gpuAlone[f], wetAlone) < 2e-3, f);

  const cpu = cpuOcean(ekman), off = cpuOcean(THIN_FLOORS);
  gpuOcean.initialize(surfaceT0, ice);
  for (const steps of [1, 20]) {
    for (let n = 0; n < steps; n++) {
      cpu.advance(Float64Array.from(surfaceT0), ice, flux, stress, dt);
      off.advance(Float64Array.from(surfaceT0), ice, flux, stress, dt);
      await gpuOcean.advance(Float64Array.from(surfaceT0), ice, stress, dt);
    }
    const cpuState = cpu.serialize(), offState = off.serialize(), gpuState = await gpuOcean.serialize();
    const wet = Array.from(cpuState.h, (v) => v > 1), upper = Array.from(cpuState.h, (_, x) => x < 6 * C);
    const lines = [];
    for (const f of ['h', 'u', 'T', 'S', 'eta']) {
      const mask = f === 'T' || f === 'S' ? wet : null, gap = rms(cpuState[f], gpuState[f], mask);
      lines.push(`${f} ${gap.toExponential(2)}`);
      assert.ok(gap < 2e-3, `${f} rms relative diff ${gap} after ${steps === 1 ? 1 : 21} steps`);
    }
    const gap = rms(cpuState.h, gpuState.h, upper), effect = rms(cpuState.h, offState.h, upper);
    console.log(`${steps === 1 ? 1 : 21} ocean steps with ekmanBlock under 10 m floors at N=8, the engines' rms relative difference: ${lines.join(', ')}; in the mixed layer and the five classes under it ${gap.toExponential(2)} against the Ekman block's effect ${effect.toExponential(2)}`);
    assert.ok(effect > 5 * gap, `the Ekman block changed the upper layers by ${effect}, no more than the engines differ (${gap})`);
  }
  for (const v of (await gpuOcean.download()).u) assert.ok(Number.isFinite(v));
  assert.equal((await gpuOcean.diagnostics()).oceanLimited, 0);
  gpuModel.destroy();
});

const UPWELLING = 6, UPWELLING_STEPPED = 100;

test('under upwellingSource the GPU mixed layer draws the Ekman divergence of a coastal cell from the class under the block\'s base as the CPU\'s does, alone and through twenty-one ocean steps at N=8', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { cpuModel, surfaceT0, ice, stress: weak } = buildScenario(8);
  const stress = Float64Array.from(weak, (x) => 10 * x);
  const mesh = cpuModel.mesh, C = mesh.nCells, E = mesh.nEdges, dt = 1350, rate = 3 / 86400;
  const rms = (a, b, mask = null) => { let d = 0, r = 0; for (let x = 0; x < a.length; x++) if (!mask || mask[x]) { d += (a[x] - b[x]) ** 2; r += a[x] * a[x]; } return Math.sqrt(d / Math.max(r, 1e-300)); };
  const coastal = new Uint8Array(C);
  for (let i = 0; i < C; i++) {
    if (!cpuModel.geography.land[i]) for (let m = 0; m < mesh.nEdgesOnCell[i]; m++) if (cpuModel.geography.land[mesh.cellsOnCell[mesh.maxEdges * i + m]]) coastal[i] = 1;
  }

  const alone = { ...OCEAN_OPTIONS, ekmanDepth: 0, momentumDepth: 100, refillRate: rate, ekmanBlock: true };
  const cpuOcean = createCpuLayeredOcean(mesh, { geography: cpuModel.geography, ...alone });
  const surfaceT = Float64Array.from(surfaceT0);
  cpuOcean.initialize(surfaceT, ice);
  for (let i = 0; i < C; i++) {
    if (!cpuOcean.cellOcean[i]) continue;
    const amount = cpuOcean.h[i] - 30, f = amount / cpuOcean.h[i];
    cpuOcean.h[C + i] += amount; cpuOcean.Q[C + i] += cpuOcean.Q[i] * f; cpuOcean.W[C + i] += cpuOcean.W[i] * f;
    cpuOcean.h[i] -= amount; cpuOcean.Q[i] -= cpuOcean.Q[i] * f; cpuOcean.W[i] -= cpuOcean.W[i] * f;
  }
  const saved = cpuOcean.serialize();
  const onCpu = (source) => {
    const ocean = createCpuLayeredOcean(mesh, { geography: cpuModel.geography, ...alone, upwellingSource: source });
    ocean.load(saved, Float64Array.from(surfaceT), ice);
    ocean.setStress(stress, ice);
    ocean.readSurface(Float64Array.from(surfaceT), ice);
    ocean.ekmanTransport(dt);
    ocean.mixedLayer(dt);
    return ocean;
  };
  const onGpu = async (source) => {
    const gpuModel = await createGpuModel(new Grid(8), { topography, ocean: { ...alone, upwellingSource: source } });
    const gpuOcean = gpuModel.oceanEngine;
    await gpuOcean.upload(saved, Float64Array.from(surfaceT), ice);
    gpuOcean.setStress(stress, ice);
    gpuOcean.readSurface(Float64Array.from(surfaceT), ice);
    gpuOcean.phases.ekman(dt);
    gpuOcean.mixedLayer(dt);
    const out = await gpuOcean.download();
    gpuModel.destroy();
    return out;
  };
  const cpuOn = onCpu(UPWELLING), cpuOff = onCpu(0), gpuOn = await onGpu(UPWELLING), gpuOff = await onGpu(0);
  const cpuMove = Float64Array.from(cpuOn.h, (v, x) => v - cpuOff.h[x]), gpuMove = Float64Array.from(gpuOn.h, (v, x) => v - gpuOff.h[x]);
  let moved = 0, drawn = 0, inland = 0;
  for (let x = 0; x < cpuMove.length; x++) {
    if (Math.abs(cpuMove[x]) > 1e-9) { moved++; drawn = Math.max(drawn, Math.abs(cpuMove[x])); }
    if (Math.abs(cpuMove[x]) > 1e-9 && !coastal[x % C]) inland++;
  }
  const gap = rms(cpuMove, gpuMove), whole = rms(cpuOn.h, gpuOn.h), outside = rms(Float64Array.from(cpuMove, (v, x) => (coastal[x % C] ? 0 : v)), Float64Array.from(gpuMove, (v, x) => (coastal[x % C] ? 0 : v)));
  console.log(`upwellingSource ${UPWELLING} at N=8 on cells cut to 30 m under 100 m wind depth: the CPU moved ${moved} layer values by up to ${drawn.toFixed(3)} m, none away from a coast (${inland}); the GPU's redirect differs from the CPU's by ${gap.toExponential(2)} rms relative, the whole h by ${whole.toExponential(2)}`);
  assert.ok(moved >= 20 && drawn > 0.01 && inland === 0, `${moved} values moved by up to ${drawn} m, ${inland} away from a coast`);
  assert.ok(gap < 2e-2 && whole < 2e-3, `the redirect differs by ${gap}, h by ${whole}`);
  assert.ok(Math.max(...gpuMove.filter((_, x) => !coastal[x % C]).map(Math.abs)) < 1e-4, 'the GPU moves water only next to land');

  const thin = { ...THIN_FLOORS, ekmanBlock: true };
  const cpuStep = (source) => { const ocean = createCpuLayeredOcean(mesh, { geography: cpuModel.geography, ...thin, upwellingSource: source }); ocean.initialize(surfaceT0, ice); return ocean; };
  const gpuModel = await createGpuModel(new Grid(8), { topography, ocean: { ...thin, upwellingSource: UPWELLING_STEPPED } });
  const gpuStep = gpuModel.oceanEngine;
  const cpu = cpuStep(UPWELLING_STEPPED), base = cpuStep(0), flux = new Float64Array(C);
  gpuStep.initialize(surfaceT0, ice);
  for (const steps of [1, 20]) {
    for (let n = 0; n < steps; n++) {
      cpu.advance(Float64Array.from(surfaceT0), ice, flux, stress, dt);
      base.advance(Float64Array.from(surfaceT0), ice, flux, stress, dt);
      await gpuStep.advance(Float64Array.from(surfaceT0), ice, stress, dt);
    }
    const cpuState = cpu.serialize(), baseState = base.serialize(), gpuState = await gpuStep.serialize();
    const wet = Array.from(cpuState.h, (v) => v > 1), upper = Array.from(cpuState.h, (_, x) => x < 6 * C);
    const lines = [];
    for (const f of ['h', 'u', 'T', 'S', 'eta']) {
      const mask = f === 'T' || f === 'S' ? wet : null, difference = rms(cpuState[f], gpuState[f], mask);
      lines.push(`${f} ${difference.toExponential(2)}`);
      assert.ok(difference < 2e-3, `${f} rms relative diff ${difference} after ${steps === 1 ? 1 : 21} steps`);
    }
    const closeness = rms(cpuState.h, gpuState.h, upper), effect = rms(cpuState.h, baseState.h, upper);
    console.log(`${steps === 1 ? 1 : 21} ocean steps with upwellingSource ${UPWELLING_STEPPED} under 10 m floors at N=8, the engines' rms relative difference: ${lines.join(', ')}; in the mixed layer and the five classes under it ${closeness.toExponential(2)} against the source's effect ${effect.toExponential(2)}`);
    if (steps > 1) assert.ok(effect > 5 * closeness, `upwellingSource changed the upper layers by ${effect}, no more than the engines differ (${closeness})`);
  }
  for (const v of (await gpuStep.download()).u) assert.ok(Number.isFinite(v));
  assert.equal((await gpuStep.diagnostics()).oceanLimited, 0);
  gpuModel.destroy();
});

test('under closureFill the GPU closure acts on the class flow carried onto token edges as the CPU\'s does, through twenty-one ocean steps at N=8', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const filling = { ...OCEAN_OPTIONS, closureFill: 0.5, closureTokens: 'beside', closureHours: 1 };
  const { cpuModel, surfaceT0, ice, stress } = buildScenario(8);
  const mesh = cpuModel.mesh, C = mesh.nCells, dt = 1350;
  const gpuModel = await createGpuModel(new Grid(8), { topography, ocean: filling });
  const gpuOcean = gpuModel.oceanEngine;
  const cpuOcean = (options) => { const ocean = createCpuLayeredOcean(mesh, { geography: cpuModel.geography, ...options }); ocean.initialize(surfaceT0, ice); return ocean; };
  const rms = (a, b, mask = null) => { let d = 0, r = 0; for (let x = 0; x < a.length; x++) if (!mask || mask[x]) { d += (a[x] - b[x]) ** 2; r += a[x] * a[x]; } return Math.sqrt(d / Math.max(r, 1e-300)); };
  const cpu = cpuOcean(filling), pulled = cpuOcean({ ...filling, closureFill: 0 }), flux = new Float64Array(C);
  gpuOcean.initialize(surfaceT0, ice);
  for (let n = 0; n < 21; n++) {
    cpu.advance(Float64Array.from(surfaceT0), ice, flux, stress, dt);
    pulled.advance(Float64Array.from(surfaceT0), ice, flux, stress, dt);
    await gpuOcean.advance(Float64Array.from(surfaceT0), ice, stress, dt);
  }
  const cpuState = cpu.serialize(), pulledState = pulled.serialize(), gpuState = await gpuOcean.serialize();
  const wet = Array.from(cpuState.h, (v) => v > 1);
  const lines = [];
  for (const f of ['h', 'u', 'T', 'S', 'eta']) {
    const mask = f === 'T' || f === 'S' ? wet : null, gap = rms(cpuState[f], gpuState[f], mask), effect = rms(cpuState[f], pulledState[f], mask);
    lines.push(`${f} ${gap.toExponential(2)} against ${effect.toExponential(2)}`);
    assert.ok(gap < 2e-3, `${f} rms relative diff ${gap} after 21 steps`);
    if (f === 'u') assert.ok(effect > 5 * gap, `u: the fill changed the CPU ocean by ${effect}, no more than the engines differ (${gap})`);
  }
  console.log(`21 ocean steps with closureFill 0.5 and closureHours 1 at N=8, the engines' rms relative difference against the fill's effect: ${lines.join(', ')}`);
  gpuModel.destroy();
});

test('the GPU closure fills a class\'s token edges between thicker classes, in two rings, as the CPU\'s does, through twenty-one ocean steps at N=8', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const filling = { ...OCEAN_OPTIONS, closureFill: 1, closureTokens: 'interior', closureHours: 1 };
  const { cpuModel, surfaceT0, ice, stress } = buildScenario(8);
  const mesh = cpuModel.mesh, C = mesh.nCells, dt = 1350;
  const gpuModel = await createGpuModel(new Grid(8), { topography, ocean: filling });
  const gpuOcean = gpuModel.oceanEngine;
  const cpuOcean = (options) => { const ocean = createCpuLayeredOcean(mesh, { geography: cpuModel.geography, ...options }); ocean.initialize(surfaceT0, ice); return ocean; };
  const rms = (a, b, mask = null) => { let d = 0, r = 0; for (let x = 0; x < a.length; x++) if (!mask || mask[x]) { d += (a[x] - b[x]) ** 2; r += a[x] * a[x]; } return Math.sqrt(d / Math.max(r, 1e-300)); };
  const cpu = cpuOcean(filling), pulled = cpuOcean({ ...filling, closureFill: 0 }), flux = new Float64Array(C);
  gpuOcean.initialize(surfaceT0, ice);
  for (let n = 0; n < 21; n++) {
    cpu.advance(Float64Array.from(surfaceT0), ice, flux, stress, dt);
    pulled.advance(Float64Array.from(surfaceT0), ice, flux, stress, dt);
    await gpuOcean.advance(Float64Array.from(surfaceT0), ice, stress, dt);
  }
  const cpuState = cpu.serialize(), pulledState = pulled.serialize(), gpuState = await gpuOcean.serialize();
  const wet = Array.from(cpuState.h, (v) => v > 1);
  const lines = [];
  for (const f of ['h', 'u', 'T', 'S', 'eta']) {
    const mask = f === 'T' || f === 'S' ? wet : null, gap = rms(cpuState[f], gpuState[f], mask), effect = rms(cpuState[f], pulledState[f], mask);
    lines.push(`${f} ${gap.toExponential(2)} against ${effect.toExponential(2)}`);
    assert.ok(gap < 2e-3, `${f} rms relative diff ${gap} after 21 steps`);
    if (f === 'u') assert.ok(effect > 5 * gap, `u: the fill changed the CPU ocean by ${effect}, no more than the engines differ (${gap})`);
  }
  console.log(`21 ocean steps with the interior fill and closureHours 1 at N=8, the engines' rms relative difference against the fill's effect: ${lines.join(', ')}`);
  gpuModel.destroy();
});

test('under shearMixing the GPU interfacial drag follows the Richardson number as the CPU\'s does, through twenty-one ocean steps at N=8', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const mixing = { ...OCEAN_OPTIONS, shearMixing: true, interfacialDrag: 0 };
  const { cpuModel, surfaceT0, ice, stress } = buildScenario(8);
  const mesh = cpuModel.mesh, C = mesh.nCells, dt = 1350;
  const gpuModel = await createGpuModel(new Grid(8), { topography, ocean: mixing });
  const gpuOcean = gpuModel.oceanEngine;
  const cpuOcean = (options) => { const ocean = createCpuLayeredOcean(mesh, { geography: cpuModel.geography, ...options }); ocean.initialize(surfaceT0, ice); return ocean; };
  const rms = (a, b, mask = null) => { let d = 0, r = 0; for (let x = 0; x < a.length; x++) if (!mask || mask[x]) { d += (a[x] - b[x]) ** 2; r += a[x] * a[x]; } return Math.sqrt(d / Math.max(r, 1e-300)); };
  const cpu = cpuOcean(mixing), constant = cpuOcean({ ...mixing, shearMixing: false, interfacialDrag: 5e-5 }), flux = new Float64Array(C);
  gpuOcean.initialize(surfaceT0, ice);
  for (let n = 0; n < 21; n++) {
    cpu.advance(Float64Array.from(surfaceT0), ice, flux, stress, dt);
    constant.advance(Float64Array.from(surfaceT0), ice, flux, stress, dt);
    await gpuOcean.advance(Float64Array.from(surfaceT0), ice, stress, dt);
  }
  const cpuState = cpu.serialize(), constantState = constant.serialize(), gpuState = await gpuOcean.serialize();
  const wet = Array.from(cpuState.h, (v) => v > 1);
  const lines = [];
  for (const f of ['h', 'u', 'T', 'S', 'eta']) {
    const mask = f === 'T' || f === 'S' ? wet : null, gap = rms(cpuState[f], gpuState[f], mask), effect = rms(cpuState[f], constantState[f], mask);
    lines.push(`${f} ${gap.toExponential(2)} against ${effect.toExponential(2)}`);
    assert.ok(gap < 2e-3, `${f} rms relative diff ${gap} after 21 steps`);
    if (f === 'u') assert.ok(effect > 5 * gap, `u: the shear mixing changed the CPU ocean by ${effect}, no more than the engines differ (${gap})`);
  }
  console.log(`21 ocean steps with shearMixing at N=8, the engines' rms relative difference against the mixing's effect over the constant drag: ${lines.join(', ')}`);
  gpuModel.destroy();
});

test('under interiorShearMixing the GPU drag between interior classes follows their Richardson number as the CPU\'s does, through twenty-one ocean steps at N=8', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const mixing = { ...OCEAN_OPTIONS, interiorShearMixing: true, closureFill: 0 };
  const { cpuModel, surfaceT0, ice, stress } = buildScenario(8);
  const mesh = cpuModel.mesh, C = mesh.nCells, dt = 1350;
  const gpuModel = await createGpuModel(new Grid(8), { topography, ocean: mixing });
  const gpuOcean = gpuModel.oceanEngine;
  const cpuOcean = (options) => { const ocean = createCpuLayeredOcean(mesh, { geography: cpuModel.geography, ...options }); ocean.initialize(surfaceT0, ice); return ocean; };
  const rms = (a, b, mask = null) => { let d = 0, r = 0; for (let x = 0; x < a.length; x++) if (!mask || mask[x]) { d += (a[x] - b[x]) ** 2; r += a[x] * a[x]; } return Math.sqrt(d / Math.max(r, 1e-300)); };
  const cpu = cpuOcean(mixing), constant = cpuOcean({ ...mixing, interiorShearMixing: false }), flux = new Float64Array(C);
  gpuOcean.initialize(surfaceT0, ice);
  for (let n = 0; n < 21; n++) {
    cpu.advance(Float64Array.from(surfaceT0), ice, flux, stress, dt);
    constant.advance(Float64Array.from(surfaceT0), ice, flux, stress, dt);
    await gpuOcean.advance(Float64Array.from(surfaceT0), ice, stress, dt);
  }
  const cpuState = cpu.serialize(), constantState = constant.serialize(), gpuState = await gpuOcean.serialize();
  const wet = Array.from(cpuState.h, (v) => v > 1);
  const lines = [];
  for (const f of ['h', 'u', 'T', 'S', 'eta']) {
    const mask = f === 'T' || f === 'S' ? wet : null, gap = rms(cpuState[f], gpuState[f], mask), effect = rms(cpuState[f], constantState[f], mask);
    lines.push(`${f} ${gap.toExponential(2)} against ${effect.toExponential(2)}`);
    assert.ok(gap < 2e-3, `${f} rms relative diff ${gap} after 21 steps`);
    if (f === 'u') assert.ok(effect > 5 * gap, `u: the interior mixing changed the CPU ocean by ${effect}, no more than the engines differ (${gap})`);
  }
  console.log(`21 ocean steps with interiorShearMixing at N=8, the engines' rms relative difference against its effect over the constant drag: ${lines.join(', ')}`);
  gpuModel.destroy();
});

test('the GPU mixed layer holds, retreats, convects, keeps to its neighbors and, under the old options, snaps back as the CPU one does, carrying the momentum of the water it moves as the CPU one does and keeping each edge\'s column momentum', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { cpuModel, surfaceT0, ice, stress } = buildScenario(8);
  const mesh = cpuModel.mesh, C = mesh.nCells, E = mesh.nEdges;
  const sea = [...Array(C).keys()].filter((i) => cpuModel.geography.land[i] === 0);
  const neutral = sea.filter((_, n) => n % 5 === 0), stable = sea.filter((_, n) => n % 5 === 1), dense = sea.filter((_, n) => n % 5 === 2);
  const oldRules = { neutralSnap: true, convectiveErosion: false, buoyancyMemory: 0, maximumMixedDepth: 200, mixedNeighbourRatio: 0, vorticityCentring: 0 };
  const unsheared = { shearEntrainment: 0 };
  for (const options of [{ ...OCEAN_OPTIONS, ...unsheared, mixedNeighbourRatio: 3, carryMomentum: true }, { ...OCEAN_OPTIONS, ...unsheared, carryMomentum: false, ...oldRules }]) {
    const cpuOcean = createCpuLayeredOcean(mesh, { geography: cpuModel.geography, ...options });
    const surfaceT = Float64Array.from(surfaceT0);
    cpuOcean.initialize(surfaceT, ice);
    const { h, Q, W, densities, layers: L } = cpuOcean;
    for (const [cells, depth, offset] of [[neutral, 300, -0.002], [stable, 600, -0.03], [dense, 50, 0.27]]) {
      for (const i of cells) {
        for (let k = 1; k < L && h[i] < depth; k++) {
          const a = k * C + i, take = Math.min(h[a] - 0.01, depth - h[i]);
          if (take <= 0) continue;
          const f = take / h[a];
          h[i] += take; Q[i] += Q[a] * f; W[i] += W[a] * f; h[a] -= take; Q[a] -= Q[a] * f; W[a] -= W[a] * f;
        }
        let below = 1;
        while (below < L - 1 && h[below * C + i] <= 5) below++;
        const t = labelTemperature(densities[below] + offset, W[i] / h[i]);
        Q[i] = h[i] * t; surfaceT[i] = t;
      }
    }
    for (let n = 0; n < cpuOcean.u.length; n++) cpuOcean.u[n] = 0.2 * Math.sin(0.7 * Math.floor(n / E) + 0.37 * (n % E));
    const saved = cpuOcean.serialize();
    cpuOcean.load(saved, surfaceT, ice);
    const still = createCpuLayeredOcean(mesh, { geography: cpuModel.geography, ...options, carryMomentum: false, shearEntrainment: 0 });
    still.load(saved, Float64Array.from(surfaceT), ice);
    const gpuModel = await createGpuModel(new Grid(8), { topography, ocean: options });
    const gpuOcean = gpuModel.oceanEngine;
    await gpuOcean.upload(saved, Float64Array.from(surfaceT), ice);
    const warmed = Float64Array.from(surfaceT, (t, i) => (neutral.includes(i) ? t + 0.05 : t));
    cpuOcean.advance(Float64Array.from(warmed), ice, new Float64Array(C), stress, 1350);
    still.advance(Float64Array.from(warmed), ice, new Float64Array(C), stress, 1350);
    gpuOcean.readSurface(Float64Array.from(warmed), ice);
    gpuOcean.setStress(stress, ice);
    gpuOcean.step(1350);
    const stepped = await gpuOcean.download();
    gpuOcean.mixedLayer(1350);
    const mixed = await gpuOcean.download();
    gpuOcean.salt(1350);
    gpuOcean.writeSurface(1350);
    await gpuModel.gpu.device.queue.onSubmittedWorkDone();
    const before = columnMomentum(mesh, gpuOcean.layers, stepped, gpuOcean.D), after = columnMomentum(mesh, gpuOcean.layers, mixed, gpuOcean.D);
    let kept = 0;
    for (let e = 0; e < E; e++) if (cpuModel.geography.edgeOcean[e]) kept = Math.max(kept, Math.abs(after.momentum[e] - before.momentum[e]) / before.scale[e]);
    const cpu = cpuOcean.serialize(), gpu = await gpuOcean.serialize(), stillState = still.serialize();
    const uGap = stats(cpu.u, gpu.u).rmsRel, uEffect = stats(cpu.u, stillState.u).rmsRel;
    console.log(`${options.neutralSnap ? 'old' : 'new'} rules${options.carryMomentum ? ', momentum carried' : ''}: the engines' u within ${uGap.toExponential(2)} rms relative, the carriage moving it by ${uEffect.toExponential(2)}; the GPU's exchanges change an edge's column momentum by ${kept.toExponential(2)} of its Σ h|u| at most`);
    assert.ok(uGap < 2e-3, `u rms relative diff ${uGap}`);
    if (options.carryMomentum) assert.ok(uEffect > 5 * uGap && kept < 1e-5, `the carriage moved u by ${uEffect} against the engines' ${uGap}; the exchanges changed a column's momentum by ${kept}`);
    else assert.ok(uEffect === 0 && kept > 1e-3, `without the carriage u moved by ${uEffect} and the columns by ${kept}`);
    let worst = 0, at = -1;
    for (let i = 0; i < C; i++) { const d = Math.abs(cpu.h[i] - gpu.h[i]); if (d > worst) { worst = d; at = i; } }
    const depths = (cells) => cells.map((i) => cpu.h[i].toFixed(0)).slice(0, 4).join(', ');
    console.log(`${options.neutralSnap ? 'old' : 'new'} rules, one step: neutral 300 m layers warmed 0.05 K → ${depths(neutral)} m, stable 600 m → ${depths(stable)} m, dense 50 m → ${depths(dense)} m; GPU mixed-layer depth within ${worst.toExponential(2)} m (cell ${at})`);
    assert.ok(worst < 0.05, `mixed-layer depth differs by ${worst} m at cell ${at}: cpu ${cpu.h[at]} gpu ${gpu.h[at]}`);
    const hStats = stats(cpu.h, gpu.h);
    assert.ok(hStats.rmsRel < 2e-3, `h rms relative diff ${hStats.rmsRel}`);
    const wet = Array.from(cpu.h, (v) => v > 1);
    for (const f of ['T', 'S']) {
      const r = stats(cpu[f].filter((_, x) => wet[x]), gpu[f].filter((_, x) => wet[x]));
      assert.ok(r.rmsRel < 2e-3, `${f} rms relative diff ${r.rmsRel}`);
    }
  }
});

test('under refillRate the GPU mixed layer refills a column in deficit by refillRate·dt and one below refillFloor to it, as the CPU one does, and reports the same thinned share', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { cpuModel, surfaceT0, ice } = buildScenario(8);
  const mesh = cpuModel.mesh, C = mesh.nCells, dt = 1350, rate = 24 / 86400;
  const options = { ...OCEAN_OPTIONS, refillRate: rate };
  const cpuOcean = createCpuLayeredOcean(mesh, { geography: cpuModel.geography, ...options });
  const surfaceT = Float64Array.from(surfaceT0);
  cpuOcean.initialize(surfaceT, ice);
  const { h, Q, W } = cpuOcean, sea = [...Array(C).keys()].filter((i) => cpuOcean.cellOcean[i]);
  const deficit = sea.filter((_, n) => n % 3 === 0), backstop = sea.filter((_, n) => n % 3 === 1);
  for (const [cells, depth] of [[deficit, 30], [backstop, 3]]) {
    for (const i of cells) {
      const amount = h[i] - depth, f = amount / h[i];
      h[C + i] += amount; Q[C + i] += Q[i] * f; W[C + i] += W[i] * f; h[i] -= amount; Q[i] -= Q[i] * f; W[i] -= W[i] * f;
    }
  }
  const saved = cpuOcean.serialize(), calm = new Float64Array(mesh.nEdges);
  const unfilled = createCpuLayeredOcean(mesh, { geography: cpuModel.geography, ...OCEAN_OPTIONS, minimumThickness: 1 });
  for (const ocean of [cpuOcean, unfilled]) {
    ocean.load(saved, Float64Array.from(surfaceT), ice);
    ocean.setStress(calm, ice);
    ocean.readSurface(Float64Array.from(surfaceT), ice);
    ocean.mixedLayer(dt);
  }
  const onGpu = async (ocean) => {
    const gpuModel = await createGpuModel(new Grid(8), { topography, ocean });
    const gpuOcean = gpuModel.oceanEngine;
    await gpuOcean.upload(saved, Float64Array.from(surfaceT), ice);
    gpuOcean.setStress(calm, ice);
    gpuOcean.readSurface(Float64Array.from(surfaceT), ice);
    gpuOcean.mixedLayer(dt);
    const out = { ...(await gpuOcean.download()), thinned: (await gpuOcean.diagnostics()).oceanThinned };
    gpuModel.destroy();
    return out;
  };
  const gpu = await onGpu(options), gpuUnfilled = await onGpu({ ...OCEAN_OPTIONS, minimumThickness: 1 }), cpuThinned = cpuOcean.diagnostics().oceanThinned;
  let cpuOff = 0, gpuOff = 0, gap = 0, capped = 0, floored = 0;
  const expected = (before) => Math.max(Math.min(50, before + rate * dt), 5);
  for (const i of [...deficit, ...backstop]) {
    const reached = unfilled.h[i] + rate * dt;
    if (reached < 5) floored++; else if (reached < 50) capped++;
    cpuOff = Math.max(cpuOff, Math.abs(cpuOcean.h[i] - expected(unfilled.h[i])));
    gpuOff = Math.max(gpuOff, Math.abs(gpu.h[i] - expected(gpuUnfilled.h[i])));
  }
  let refillGap = 0, at = -1;
  for (const i of sea) {
    refillGap = Math.max(refillGap, Math.abs((gpu.h[i] - gpuUnfilled.h[i]) - (cpuOcean.h[i] - unfilled.h[i])));
    if (Math.abs(gpu.h[i] - cpuOcean.h[i]) > gap) { gap = Math.abs(gpu.h[i] - cpuOcean.h[i]); at = i; }
  }
  console.log(`refill at ${rate * 86400} m/day: ${deficit.length} columns cut to 30 m and ${backstop.length} to 3 m; on the CPU ${capped} took refillRate·dt over what the rest of the step left and ${floored} rose to refillFloor, each to ${cpuOff.toExponential(2)} m, on the GPU to ${gpuOff.toExponential(2)} m; the engines' refills within ${refillGap.toExponential(2)} m of each other, their mixed layers within ${gap.toExponential(2)} m (cell ${at}: ${cpuOcean.h[at]} against ${gpu.h[at]} m, before the refill ${unfilled.h[at]} against ${gpuUnfilled.h[at]}); thinned ${cpuThinned.toFixed(5)} on the CPU, ${gpu.thinned.toFixed(5)} on the GPU`);
  assert.ok(capped > 100 && floored > 100, `${capped} capped and ${floored} floored columns`);
  assert.ok(cpuOff < 1e-9 && gpuOff < 1e-4, `the mixed layers are off by ${cpuOff} m on the CPU, ${gpuOff} m on the GPU`);
  assert.ok(refillGap < 1e-4 && gap < 0.2, `the engines' refills differ by ${refillGap} m, their mixed layers by ${gap} m`);
  assert.ok(cpuThinned > 0.5 && Math.abs(gpu.thinned - cpuThinned) < 1e-6, `thinned ${gpu.thinned} on the GPU against ${cpuThinned}`);
});

test('under shearEntrainment the GPU mixed layer takes in the class beneath as the CPU one does, carrying its momentum, and twenty-one steps under 10 m floors track the CPU ocean at N=8', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const ric = 0.65, sheared = { ...OCEAN_OPTIONS, mixedDepth: 10, minimumThickness: 10, shallowestMixedDepth: 10, momentumDepth: 50, carryMomentum: true, shearEntrainment: ric };
  const { cpuModel, surfaceT0, ice, stress } = buildScenario(8);
  const mesh = cpuModel.mesh, C = mesh.nCells, E = mesh.nEdges, dt = 1350, calm = new Float64Array(E);
  const cpuOcean = (options) => createCpuLayeredOcean(mesh, { geography: cpuModel.geography, ...options });
  const seeded = cpuOcean(sheared);
  seeded.initialize(Float64Array.from(surfaceT0), ice);
  for (let e = 0; e < E; e++) {
    if (!seeded.edgeOcean[e] || Math.abs(mesh.latEdge[e]) > 40 * Math.PI / 180) continue;
    const x = mesh.xEdge[3 * e], y = mesh.xEdge[3 * e + 1], r = Math.hypot(x, y) || 1e-12;
    seeded.u[e] += 0.6 * ((-y / r) * mesh.nEdge[3 * e] + (x / r) * mesh.nEdge[3 * e + 1]);
  }
  const saved = seeded.serialize();
  const L = seeded.layers, { h: h0, u: u0, T, S } = saved, rho = [1025, ...seeded.densities];
  const du = new Float64Array(E);
  for (let e = 0; e < E; e++) {
    const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1];
    if (seeded.edgeOcean[e]) for (let k = 1; k < L; k++) if (Math.min(h0[k * C + a], h0[k * C + b]) >= 5) { du[e] = u0[e] - u0[k * E + e]; break; }
  }
  const v = cellVector(mesh, du), rib = new Float64Array(C).fill(Infinity);
  for (let i = 0; i < C; i++) {
    if (!seeded.cellOcean[i]) continue;
    let below = 0;
    for (let k = 1; k < L; k++) if (h0[k * C + i] > 5) { below = k; break; }
    if (below) rib[i] = Math.max(9.81 * (rho[below] - seawaterDensity(T[i], S[i])) / 1025, 9.81 * 0.005 / 1025) * h0[i] / (v[3 * i] ** 2 + v[3 * i + 1] ** 2 + v[3 * i + 2] ** 2 + 1e-12);
  }
  const cpu = cpuOcean(sheared);
  cpu.load(saved, Float64Array.from(surfaceT0), ice);
  cpu.setStress(calm, ice);
  cpu.readSurface(Float64Array.from(surfaceT0), ice);
  cpu.mixedLayer(dt);
  const gpuModel = await createGpuModel(new Grid(8), { topography, ocean: sheared });
  const gpuOcean = gpuModel.oceanEngine;
  await gpuOcean.upload(saved, Float64Array.from(surfaceT0), ice);
  gpuOcean.setStress(calm, ice);
  gpuOcean.readSurface(Float64Array.from(surfaceT0), ice);
  gpuOcean.mixedLayer(dt);
  const gpu = await gpuOcean.download();
  let entrained = 0, near = 0, hGap = 0, hRel = 0, takenGap = 0, uGap = 0;
  for (let i = 0; i < C; i++) {
    if (!cpu.cellOcean[i]) continue;
    hGap = Math.max(hGap, Math.abs(gpu.h[i] - cpu.h[i]));
    hRel = Math.max(hRel, Math.abs(gpu.h[i] - cpu.h[i]) / cpu.h[i]);
    const taken = cpu.sheared[i];
    if (Math.abs(rib[i] / ric - 1) < 0.01) { near++; continue; }
    if (taken > 0) entrained++;
    takenGap = Math.max(takenGap, Math.abs((gpu.h[i] - h0[i]) - taken) / Math.max(taken, 1e-3));
  }
  for (let x = 0; x < L * E; x++) uGap = Math.max(uGap, Math.abs(gpu.u[x] - cpu.u[x]));
  console.log(`shear entrainment, one mixed-layer step: ${entrained} cells took in ${cpu.sheared.reduce((s, x) => s + x, 0).toFixed(1)} m on the CPU (${near} within 1 % of Ri_c set aside); the engines' mixed layers differ by ${hGap.toExponential(2)} m, ${hRel.toExponential(2)} of their depth, the entrained depth by ${takenGap.toExponential(2)} of itself, the velocities by ${uGap.toExponential(2)} m/s`);
  assert.ok(entrained > 50, `${entrained} cells entrained`);
  assert.ok(hRel < 1e-5 && takenGap < 1e-3 && uGap < 1e-5, `h ${hRel} of itself, entrained ${takenGap}, u ${uGap} m/s`);

  const heatedOptions = { ...sheared, minimumThickness: 5, shallowestMixedDepth: 5, detrainmentTime: dt, buoyancyMemory: 0 }, weak = Float64Array.from(stress, (x) => 0.02 * x);
  const heatedCpu = cpuOcean(heatedOptions);
  heatedCpu.load(saved, Float64Array.from(surfaceT0), ice);
  const warm = Float64Array.from(surfaceT0, (t, i) => t + dt * 150 / heatedCpu.capacity[i]);
  heatedCpu.setStress(weak, ice);
  heatedCpu.readSurface(Float64Array.from(warm), ice);
  heatedCpu.mixedLayer(dt);
  const heatedModel = await createGpuModel(new Grid(8), { topography, ocean: heatedOptions }), heatedGpu = heatedModel.oceanEngine;
  await heatedGpu.upload(saved, Float64Array.from(surfaceT0), ice);
  heatedGpu.setStress(weak, ice);
  heatedGpu.readSurface(Float64Array.from(warm), ice);
  heatedGpu.mixedLayer(dt);
  const heatedH = (await heatedGpu.download()).h;
  heatedModel.destroy();
  let held = 0, heldCpu = 0, heldGap = 0, retreated = 0;
  for (let i = 0; i < C; i++) {
    if (!heatedCpu.cellOcean[i] || ice[i] > 0 || !(heatedCpu.retreatDepth[i] > 0) || heatedCpu.retreatDepth[i] > 5) continue;
    if (heatedCpu.h[i] < 5 + 1e-9) retreated++;
    if (!(heatedCpu.sheared[i] > 0.1) || Math.abs(rib[i] / ric - 1) < 0.05) continue;
    held++;
    heldCpu = Math.max(heldCpu, Math.abs(heatedCpu.h[i] - (h0[i] + heatedCpu.sheared[i])));
    heldGap = Math.max(heldGap, Math.abs(heatedH[i] - heatedCpu.h[i]) / heatedCpu.h[i]);
  }
  console.log(`heated under a weak wind, ${retreated} cells retreat to the 5 m floor while ${held} that shear entrains end the step at the shear's depth on the CPU within ${heldCpu.toExponential(2)} m, the GPU's within ${heldGap.toExponential(2)} of the CPU's`);
  assert.ok(held > 20 && retreated > 20, `${held} held, ${retreated} retreated`);
  assert.ok(heldCpu < 1e-9 && heldGap < 2e-4, `the CPU ${heldCpu} m off the shear's depth, the GPU ${heldGap} of it off the CPU`);

  await gpuOcean.upload(saved, Float64Array.from(surfaceT0), ice);
  const along = cpuOcean(sheared), without = cpuOcean({ ...sheared, shearEntrainment: 0 }), flux = new Float64Array(C);
  for (const ocean of [along, without]) ocean.load(saved, Float64Array.from(surfaceT0), ice);
  for (let n = 0; n < 21; n++) {
    along.advance(Float64Array.from(surfaceT0), ice, flux, stress, dt);
    without.advance(Float64Array.from(surfaceT0), ice, flux, stress, dt);
    await gpuOcean.advance(Float64Array.from(surfaceT0), ice, stress, dt);
  }
  const rms = (a, b, mask = null) => { let d = 0, r = 0; for (let x = 0; x < a.length; x++) if (!mask || mask[x]) { d += (a[x] - b[x]) ** 2; r += a[x] * a[x]; } return Math.sqrt(d / Math.max(r, 1e-300)); };
  const cpuState = along.serialize(), withoutState = without.serialize(), gpuState = await gpuOcean.serialize(), wet = Array.from(cpuState.h, (x) => x > 1);
  const lines = [];
  for (const f of ['h', 'u', 'T', 'S', 'eta']) {
    const mask = f === 'T' || f === 'S' ? wet : null, gap = rms(cpuState[f], gpuState[f], mask), effect = rms(cpuState[f], withoutState[f], mask);
    lines.push(`${f} ${gap.toExponential(2)} against ${effect.toExponential(2)}`);
    assert.ok(gap < 2e-3, `${f} rms relative diff ${gap} after 21 steps`);
    if (f === 'h' || f === 'u') assert.ok(effect > 5 * gap, `${f}: shear entrainment changed the CPU ocean by ${effect}, no more than the engines differ (${gap})`);
  }
  console.log(`21 ocean steps under 10 m floors with shear entrainment at N=8, the engines' rms relative difference against the effect of the entrainment: ${lines.join(', ')}`);
  assert.equal((await gpuOcean.diagnostics()).oceanLimited, 0);
  gpuModel.destroy();
});

test('under the mixing package the GPU ocean takes the wind over its Ekman depth and retreats to the Langmuir-enhanced Monin–Obukhov depth as the CPU ocean does, and twenty steps under 10 m floors track it at N=8', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const packaged = { ...OCEAN_OPTIONS, mixedDepth: 10, minimumThickness: 10, shallowestMixedDepth: 10, momentumDepth: 50, carryMomentum: true, shearEntrainment: 0.65, ekmanDepth: 0.4, momentumFloor: 10, stressThickness: 5, retreatStirring: 2.5, langmuir: 2, retreatEkman: 0 };
  const { cpuModel, surfaceT0, ice, stress } = buildScenario(8);
  const mesh = cpuModel.mesh, C = mesh.nCells, E = mesh.nEdges, dt = 1350, net = 150, flux = new Float64Array(C);
  const cpuOcean = (options) => createCpuLayeredOcean(mesh, { geography: cpuModel.geography, ...options });

  const retreatOn = async (options, heat = net) => {
    const retreating = { ...options, mixedDepth: 60, detrainmentTime: dt, buoyancyMemory: 0 };
    const deep = cpuOcean(retreating);
    deep.initialize(Float64Array.from(surfaceT0), ice);
    const deepState = deep.serialize(), warm = Float64Array.from(surfaceT0, (t, i) => t + dt * heat / deep.capacity[i]);
    const cpu = cpuOcean(retreating);
    cpu.load(deepState, Float64Array.from(surfaceT0), ice);
    cpu.setStress(stress, ice);
    cpu.readSurface(Float64Array.from(warm), ice);
    const before = Float64Array.from(cpu.h.subarray(0, C));
    cpu.mixedLayer(dt);
    const model = await createGpuModel(new Grid(8), { topography, ocean: retreating }), gpu = model.oceanEngine;
    await gpu.upload(deepState, Float64Array.from(surfaceT0), ice);
    gpu.setStress(stress, ice);
    gpu.readSurface(Float64Array.from(warm), ice);
    gpu.mixedLayer(dt);
    const { h } = await gpu.download();
    model.destroy();
    const reached = Array.from(before, (b, i) => cpu.cellOcean[i] && cpu.retreatDepth[i] > 10.5 && cpu.retreatDepth[i] < b - 0.5 && h[i] < b - 0.5);
    return { L: Float64Array.from(cpu.retreatDepth), gpu: h, before, reached };
  };
  const faded = await retreatOn({ ...packaged, retreatStirring: 0, langmuir: 0 }), enhanced = await retreatOn(packaged), tauCell = cellVector(mesh, stress);
  let reached = 0, lGap = 0, offset = 0, least = Infinity, most = 0;
  for (let i = 0; i < C; i++) {
    if (!faded.reached[i] || !enhanced.reached[i]) continue;
    reached++;
    lGap = Math.max(lGap, Math.abs((enhanced.gpu[i] / faded.gpu[i]) / (enhanced.L[i] / faded.L[i]) - 1));
    offset = Math.max(offset, Math.abs(enhanced.gpu[i] / enhanced.L[i] - 1));
    const epsilon = langmuirFactor(Math.sqrt(Math.hypot(tauCell[3 * i], tauCell[3 * i + 1], tauCell[3 * i + 2]) / 1025), enhanced.before[i], 2);
    least = Math.min(least, epsilon); most = Math.max(most, epsilon);
  }
  console.log(`the Langmuir retreat: ${reached} heated cells retreat from about 60 m to their L in one step under both rules, ε from ${least.toFixed(3)} to ${most.toFixed(3)}; the GPU's ratio of the 2.5ε depth to the faded rule's within ${lGap.toExponential(2)} of the CPU's, the depths themselves within ${offset.toExponential(2)} under each engine's own buoyancy of an 8e-4 K warming, as with the options off`);
  assert.ok(reached > 20 && most > least + 0.1, `${reached} cells reached L under both rules, ε from ${least} to ${most}`);
  assert.ok(lGap < 1e-5, `the engines' retreat constants differ by ${lGap} of themselves`);
  const free = await retreatOn(packaged, 30), bounded = await retreatOn({ ...packaged, retreatEkman: 0.7 }, 30);
  let binding = 0, boundGap = 0;
  for (let i = 0; i < C; i++) {
    if (!bounded.reached[i] || !(bounded.L[i] < free.L[i] / 1.5)) continue;
    binding++;
    boundGap = Math.max(boundGap, Math.abs(bounded.gpu[i] / bounded.L[i] - 1));
  }
  console.log(`under 30 W/m² ${binding} heated cells retreat in one step to the Ekman bound 0.7u*/|f| at least a third below C_r·ε·u*³/(−B), which the bound frees of the engines' buoyancies, the GPU's depths within ${boundGap.toExponential(2)} of the CPU's`);
  assert.ok(binding > 10 && boundGap < 1e-5, `${binding} cells at the bound, the engines ${boundGap} apart`);

  const iceA = Float64Array.from({ length: C }, (_, i) => (Math.abs(mesh.latCell[i]) > 50 * Math.PI / 180 ? 1 : 0)), iceB = Float64Array.from(iceA, (x) => (x > 0 ? x - 1e-3 : 0));
  const cover = Float64Array.from({ length: C }, (_, i) => (iceA[i] > 0 ? 0.1 + 0.8 * (i % 9) / 8 : 0));
  const meltOn = async (options) => {
    const melting = { ...options, mixedDepth: 60, detrainmentTime: dt, buoyancyMemory: 0 };
    const deep = cpuOcean(melting);
    deep.initialize(Float64Array.from(surfaceT0), iceA, { concentration: cover });
    const deepState = deep.serialize(), cpu = cpuOcean(melting);
    cpu.load(deepState, Float64Array.from(surfaceT0), iceA, cover);
    cpu.advance(Float64Array.from(surfaceT0), iceB, flux, stress, dt, cover);
    const model = await createGpuModel(new Grid(8), { topography, ocean: melting }), gpu = model.oceanEngine;
    await gpu.upload(deepState, Float64Array.from(surfaceT0), iceA, cover);
    await gpu.advance(Float64Array.from(surfaceT0), iceB, stress, dt, cover);
    const { h } = await gpu.download();
    model.destroy();
    return { L: Float64Array.from(cpu.retreatDepth), cpu: Float64Array.from(cpu.h.subarray(0, C)), gpu: h };
  };
  const meltPlain = await meltOn({ ...packaged, langmuir: 0 }), meltWaved = await meltOn(packaged);
  let melting = 0, iceGap = 0, waveless = Infinity;
  for (let i = 0; i < C; i++) {
    const landed = (r) => r.L[i] > 10.5 && r.L[i] < 59.5 && Math.abs(r.cpu[i] - r.L[i]) < 1e-9 * r.L[i];
    if (!(iceA[i] > 0) || !landed(meltPlain) || !landed(meltWaved) || !(meltWaved.L[i] > 1.02 * meltPlain.L[i])) continue;
    melting++;
    iceGap = Math.max(iceGap, Math.abs((meltWaved.gpu[i] / meltPlain.gpu[i]) / (meltWaved.L[i] / meltPlain.L[i]) - 1));
    waveless = Math.min(waveless, 1 - cover[i]);
  }
  console.log(`under ice covers of 0.1 to 0.9 melting 1 mm in a step, ${melting} cells retreat to a Langmuir-enhanced depth on both engines, the GPU's enhancement within ${iceGap.toExponential(2)} of the CPU's open-water share of it`);
  assert.ok(melting > 5 && iceGap < 1e-5, `${melting} melting cells, the engines' enhancements ${iceGap} apart`);

  const seeded = cpuOcean(packaged);
  seeded.initialize(Float64Array.from(surfaceT0), ice);
  for (let e = 0; e < E; e++) {
    if (!seeded.edgeOcean[e] || Math.abs(mesh.latEdge[e]) > 40 * Math.PI / 180) continue;
    const x = mesh.xEdge[3 * e], y = mesh.xEdge[3 * e + 1], r = Math.hypot(x, y) || 1e-12;
    seeded.u[e] += 0.6 * ((-y / r) * mesh.nEdge[3 * e] + (x / r) * mesh.nEdge[3 * e + 1]);
  }
  const saved = seeded.serialize();
  const gpuModel = await createGpuModel(new Grid(8), { topography, ocean: packaged });
  const gpuOcean = gpuModel.oceanEngine;
  await gpuOcean.upload(saved, Float64Array.from(surfaceT0), ice);
  const along = cpuOcean(packaged), without = cpuOcean({ ...packaged, ekmanDepth: 0, stressThickness: 20, retreatStirring: 0, langmuir: 0 });
  for (const ocean of [along, without]) ocean.load(saved, Float64Array.from(surfaceT0), ice);
  let hGap = 0, hRel = 0, uGap = 0, depthGap = 0, between = 0;
  for (let n = 0; n < 20; n++) {
    along.advance(Float64Array.from(surfaceT0), ice, flux, stress, dt);
    without.advance(Float64Array.from(surfaceT0), ice, flux, stress, dt);
    await gpuOcean.advance(Float64Array.from(surfaceT0), ice, stress, dt);
    if (n === 0) {
      const first = await gpuOcean.download();
      for (let i = 0; i < C; i++) if (along.cellOcean[i]) { hGap = Math.max(hGap, Math.abs(first.h[i] - along.h[i])); hRel = Math.max(hRel, Math.abs(first.h[i] - along.h[i]) / along.h[i]); }
      for (let x = 0; x < first.u.length; x++) uGap = Math.max(uGap, Math.abs(first.u[x] - along.u[x]));
      const [hm] = await readRanges(gpuModel.gpu.device, gpuOcean.buffers.OD, [{ offset: gpuOcean.layout.OD.HM, length: E }]);
      for (let e = 0; e < E; e++) {
        if (!along.edgeOcean[e]) continue;
        depthGap = Math.max(depthGap, Math.abs(hm[e] - along.windDepth[e]) / along.windDepth[e]);
        if (along.windDepth[e] > 10 && along.windDepth[e] < 50) between++;
      }
    }
  }
  console.log(`one packaged step: ${between} sea edges with the Ekman depth between its bounds, the GPU's within ${depthGap.toExponential(2)} of the CPU's; the mixed layers within ${hGap.toExponential(2)} m, ${hRel.toExponential(2)} of their depth, u within ${uGap.toExponential(2)} m/s`);
  assert.ok(between > 100 && depthGap < 1e-5, `${between} edges, Ekman depths ${depthGap} apart`);
  assert.ok(hRel < 3e-5 && uGap < 1e-5, `h ${hRel} of itself, u ${uGap} m/s`);
  const rms = (a, b, mask = null) => { let d = 0, r = 0; for (let x = 0; x < a.length; x++) if (!mask || mask[x]) { d += (a[x] - b[x]) ** 2; r += a[x] * a[x]; } return Math.sqrt(d / Math.max(r, 1e-300)); };
  const cpuState = along.serialize(), withoutState = without.serialize(), gpuState = await gpuOcean.serialize(), wet = Array.from(cpuState.h, (x) => x > 1);
  const lines = [];
  for (const f of ['h', 'u', 'T', 'S', 'eta']) {
    const mask = f === 'T' || f === 'S' ? wet : null, gap = rms(cpuState[f], gpuState[f], mask), effect = rms(cpuState[f], withoutState[f], mask);
    lines.push(`${f} ${gap.toExponential(2)} against ${effect.toExponential(2)}`);
    assert.ok(gap < 2e-3, `${f} rms relative diff ${gap} after 20 steps`);
    if (f === 'h' || f === 'u') assert.ok(effect > 5 * gap, `${f}: the package's Ekman depth and retreat changed the CPU ocean by ${effect}, no more than the engines differ (${gap})`);
  }
  console.log(`20 packaged ocean steps under 10 m floors at N=8, the engines' rms relative difference against the effect of the Ekman depth, the 5 m receivers and the Langmuir retreat: ${lines.join(', ')}`);
  assert.equal((await gpuOcean.diagnostics()).oceanLimited, 0);
  gpuModel.destroy();
});

/*
 * Three mixed-layer calls at N=8 under retreatDaily and retreatLight 1,
 * from a 40 m mixed layer on a 2 m floor, each handed the same sunlight
 * below the mixed layer and the same surface, warmed by the rest of the
 * sunlight less a dark loss, the second closing a day of two calls.
 */
test('under retreatDaily and retreatLight the GPU mixed layer closes the day to the CPU\'s means, the dark loss within two f32 steps of the surface temperature, and on the same means retreats under the night gate as the CPU one does', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const dt = 1350, dark = 175, suns = [750, 125, 500];
  const night = { ...OCEAN_OPTIONS, penetration: 'IB', mixedDepth: 40, minimumThickness: 2, shallowestMixedDepth: 2, detrainmentTime: dt, retreatStirring: 2.5, retreatDaily: true, retreatLight: 1 };
  const { cpuModel, surfaceT0, ice, stress } = buildScenario(8);
  const mesh = cpuModel.mesh, C = mesh.nCells, flux = new Float64Array(C);
  const cpu = createCpuLayeredOcean(mesh, { geography: cpuModel.geography, ...night });
  cpu.initialize(Float64Array.from(surfaceT0), ice);
  const start = { ...cpu.serialize(), ...cpu.restartArrays(), Q: cpu.Q, W: cpu.W, flux };
  cpu.load(start, Float64Array.from(surfaceT0), ice);
  const model = await createGpuModel(new Grid(8), { topography, ocean: night }), gpu = model.oceanEngine, PH = model.gpu.layout.PH;
  await gpu.upload(start, Float64Array.from(surfaceT0), ice);
  const open = (i) => cpu.cellOcean[i] && !(ice[i] > 0), ulp = (x) => 2 ** (Math.floor(Math.log2(Math.abs(x))) - 23);
  const warmed = (call, T0) => Float64Array.from(T0, (t, i) => (open(i) ? t + dt * (suns[call] * (1 - cpu.transmission[i]) - dark) / cpu.capacity[i] : surfaceT0[i]));
  const call = async (n, surfaceT, gpuSurfaceT, closing) => {
    const light = Float32Array.from(cpu.transmission, (t, i) => (open(i) ? suns[n] * dt * t : 0));
    cpu.penetrated.set(light);
    model.gpu.device.queue.writeBuffer(model.gpu.buffers.PH, 4 * PH.PEN, light);
    cpu.readSurface(Float64Array.from(surfaceT), ice); cpu.setStress(stress, ice);
    gpu.readSurface(gpuSurfaceT, ice); gpu.setStress(stress, ice);
    cpu.mixedLayer(dt, closing); gpu.mixedLayer(dt, closing);
    cpu.writeSurface(surfaceT, flux, dt); gpu.writeSurface(dt);
  };
  const capacity = new Float64Array(C);
  for (let n = 0; n < 2; n++) {
    for (let i = 0; i < C; i++) capacity[i] = Math.max(capacity[i], cpu.capacity[i]);
    await call(n, warmed(n, cpu.T0), warmed(n, (await gpu.restartArrays()).previousT0), n === 1 ? 1 / (2 * dt) : 0);
  }
  const days = await gpu.restartArrays();
  let sunGap = 0, darkGap = 0, steps = 0;
  for (let i = 0; i < C; i++) {
    if (!open(i) || !(cpu.daily.sunDay[i] > 0)) continue;
    const b = 9.81 * thermalExpansion(cpu.T0[i], cpu.S0[i]) / cpu.rhoCp, step = b * capacity[i] * ulp(cpu.T0[i]) / (2 * dt);
    sunGap = Math.max(sunGap, Math.abs(days.sunDay[i] / cpu.daily.sunDay[i] - 1));
    darkGap = Math.max(darkGap, Math.abs(days.darkDay[i] - cpu.daily.darkDay[i]) / Math.abs(cpu.daily.darkDay[i]));
    steps = Math.max(steps, Math.abs(days.darkDay[i] - cpu.daily.darkDay[i]) / step);
  }
  console.log(`the day closed after two calls: the GPU's mean sunlight within ${sunGap.toExponential(2)} of the CPU's, its mean dark loss within ${darkGap.toExponential(2)}, ${steps.toFixed(2)} f32 steps of T₀ over the day`);
  assert.ok(sunGap < 1e-5 && steps < 2, `S̄ ${sunGap}, D̄ ${steps} f32 steps of T₀`);

  const means = { ...start, sunSum: new Float64Array(C), darkSum: new Float64Array(C), sunDay: Float64Array.from(cpu.daily.sunDay), darkDay: Float64Array.from(cpu.daily.darkDay), daySpan: [0] };
  cpu.load(means, Float64Array.from(surfaceT0), ice);
  await gpu.upload(means, Float64Array.from(surfaceT0), ice);
  const loaded = await gpu.restartArrays(), before = { cpu: Float64Array.from(cpu.h.subarray(0, C)), gpu: (await gpu.download()).h.slice(0, C) };
  const gate = (sunDay, darkDay, depth, i) => {
    const b = 9.81 * thermalExpansion(cpu.T0[i], cpu.S0[i]) / cpu.rhoCp, column = cpu.D[i] + cpu.eta[i];
    return darkDay[i] - b * sunDay[i] * (1 - (transmitted(depth) - transmitted(column)) - depth * deposited(depth, column));
  };
  const gates = Array.from({ length: C }, (_, i) => (open(i) ? [gate(cpu.daily.sunDay, cpu.daily.darkDay, before.cpu[i], i), gate(loaded.sunDay, loaded.darkDay, before.gpu[i], i)] : null));
  const surfaceT = warmed(2, cpu.T0);
  await call(2, surfaceT, Float64Array.from(surfaceT), 0);
  const { h } = await gpu.download();
  model.destroy();
  let signs = 0, mismatched = 0, retreated = 0, depthGap = 0, movedGap = 0;
  for (let i = 0; i < C; i++) {
    if (!gates[i]) continue;
    const [onCpu, onGpu] = gates[i];
    if (Math.abs(onCpu) > 1e-12) { signs++; if (Math.sign(onCpu) !== Math.sign(onGpu)) mismatched++; }
    depthGap = Math.max(depthGap, Math.abs(h[i] / cpu.h[i] - 1));
    const moved = before.cpu[i] - cpu.h[i];
    if (cpu.retreatDepth[i] > 0 && moved > 1) { retreated++; movedGap = Math.max(movedGap, Math.abs((before.gpu[i] - h[i]) / moved - 1)); }
  }
  console.log(`on the CPU's means, the night gate: ${signs} open cells with |B_gate| > 1e-12, ${mismatched} of a different sign on the GPU; ${retreated} cells retreated over a meter, the GPU's volumes within ${movedGap.toExponential(2)} and its depths within ${depthGap.toExponential(2)} of the CPU's`);
  assert.equal(mismatched, 0);
  assert.ok(retreated > 20 && movedGap < 1e-4 && depthGap < 1e-5, `${retreated} retreating cells, volumes ${movedGap}, depths ${depthGap} apart`);
});

test('under retreatDaily without penetration the GPU mixed layer closes its days on B alone, each engine\'s surface nudged from its own T₀, as the CPU one does', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const dt = 21600, options = { ...OCEAN_OPTIONS, penetration: null, mixedDepth: 40, minimumThickness: 2, shallowestMixedDepth: 2, retreatDaily: true };
  const { cpuModel, surfaceT0, ice, stress } = buildScenario(8);
  const C = cpuModel.mesh.nCells, cpu = createCpuLayeredOcean(cpuModel.mesh, { geography: cpuModel.geography, ...options });
  const model = await createGpuModel(new Grid(8), { topography, ocean: options }), gpu = model.oceanEngine;
  cpu.initialize(Float64Array.from(surfaceT0), ice); gpu.initialize(Float64Array.from(surfaceT0), ice);
  const flux = new Float64Array(C);
  for (let n = 0; n < 8; n++) {
    const nudged = (T0) => Float64Array.from(T0, (t, i) => (cpu.cellOcean[i] && !(ice[i] > 0) ? t + (n % 4 < 2 ? 0.2 : -0.1) : surfaceT0[i]));
    const gpuSurfaceT = nudged((await gpu.restartArrays()).previousT0);
    cpu.advance(nudged(cpu.T0), ice, flux, stress, dt, undefined, (n + 1) * dt);
    await gpu.advance(gpuSurfaceT, ice, stress, dt, undefined, (n + 1) * dt);
  }
  const [{ h }, days] = await Promise.all([gpu.download(), gpu.restartArrays()]);
  model.destroy();
  let gap = 0, depthGap = 0, closed = 0;
  for (let i = 0; i < C; i++) {
    if (!cpu.cellOcean[i] || ice[i] > 0) continue;
    if (cpu.daily.darkDay[i] !== 0) closed++;
    gap = Math.max(gap, Math.abs(days.darkDay[i] - cpu.daily.darkDay[i]) / Math.max(1e-12, Math.abs(cpu.daily.darkDay[i])));
    depthGap = Math.max(depthGap, Math.abs(h[i] / cpu.h[i] - 1));
    assert.equal(days.sunDay[i], 0);
  }
  console.log(`without penetration ${closed} open cells closed two days of B; the GPU's means within ${gap.toExponential(2)} and its depths within ${depthGap.toExponential(2)} of the CPU's`);
  assert.ok(closed > 100 && gap < 1e-3 && depthGap < 5e-3, `${closed} cells, means ${gap}, depths ${depthGap} apart`);
});

test('under retreatDaily the GPU ocean saves its daily sums, and loaded half way through a day it closes the next day as it does uninterrupted; loaded without them it starts them at zero and drops the partial day', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const dt = 10800, options = { ...OCEAN_OPTIONS, penetration: 'IB', mixedDepth: 40, minimumThickness: 2, shallowestMixedDepth: 2, retreatDaily: true, retreatLight: 1 };
  const { surfaceT0, ice, stress } = buildScenario(8);
  const C = surfaceT0.length, open = (i) => !(ice[i] > 0);
  const create = async () => { const model = await createGpuModel(new Grid(8), { topography, ocean: options }); return { model, gpu: model.oceanEngine }; };
  const call = async ({ model, gpu }, n) => {
    const hour = ((n + 0.5) * dt / 3600) % 24, sun = hour > 6 && hour < 18 ? 600 * Math.sin(Math.PI * (hour - 6) / 12) : 0;
    const { previousT0, transmission } = await gpu.restartArrays();
    model.gpu.device.queue.writeBuffer(model.gpu.buffers.PH, 4 * model.gpu.layout.PH.PEN, Float32Array.from(transmission, (t, i) => (open(i) ? sun * dt * t : 0)));
    await gpu.advance(Float64Array.from(previousT0, (t, i) => (open(i) ? t + (sun > 0 ? 0.3 : -0.2) : surfaceT0[i])), ice, stress, dt, undefined, (n + 1) * dt);
  };
  const whole = await create();
  whole.gpu.initialize(Float64Array.from(surfaceT0), ice);
  for (let n = 0; n < 12; n++) await call(whole, n);
  const saved = { ...(await whole.gpu.serialize()), ...(await whole.gpu.restartArrays()) };
  assert.equal(saved.daySpan[0], DAY_SECONDS / 2);
  const without = { ...saved };
  for (const name of [...DAILY, 'daySpan']) delete without[name];
  const resumed = await create(), forgetting = await create();
  await resumed.gpu.upload(saved, Float64Array.from(surfaceT0), ice);
  await forgetting.gpu.upload(without, Float64Array.from(surfaceT0), ice);
  const loaded = await resumed.gpu.restartArrays(), cleared = await forgetting.gpu.restartArrays();
  for (const name of DAILY) {
    assert.deepEqual(Array.from(loaded[name]), Array.from(saved[name]), `the loaded ${name}`);
    assert.ok(cleared[name].every((x) => x === 0), `${name} loads as zeros without the saved sums`);
  }
  assert.ok(saved.sunDay.some((x) => x > 0) && saved.darkDay.some((x) => x !== 0), 'the first day closed before the save');
  for (let n = 12; n < 16; n++) for (const run of [whole, resumed, forgetting]) await call(run, n);
  const [a, b, c] = await Promise.all([whole, resumed, forgetting].map(async ({ gpu }) => ({ ...(await gpu.restartArrays()), h: (await gpu.download()).h })));
  for (const run of [whole, resumed, forgetting]) run.model.destroy();
  const gap = (x, y) => { let worst = 0; for (let i = 0; i < C; i++) if (open(i)) worst = Math.max(worst, Math.abs(x[i] - y[i]) / Math.max(1e-30, Math.abs(x[i]))); return worst; };
  const gaps = Object.fromEntries(['sunDay', 'darkDay', 'h'].map((name) => [name, gap(a[name], b[name])]));
  console.log(`resumed against uninterrupted after the next close, largest relative gaps: ${Object.entries(gaps).map(([name, x]) => `${name} ${x.toExponential(1)}`).join(', ')}; without the sums the day's means stay zero and the mixed layers differ by up to ${gap(a.h, c.h).toExponential(1)}`);
  assert.ok(a.darkDay.some((x, i) => x !== saved.darkDay[i]), 'the second day closed');
  for (const [name, x] of Object.entries(gaps)) assert.ok(x < 1e-5, `${name}: ${x}`);
  assert.ok(c.sunDay.every((x) => x === 0) && c.darkDay.every((x) => x === 0), 'a load without the sums leaves the partial day unclosed');
});

test('with retreatDaily and retreatLight present but off, the GPU ocean lays out the same buffers and holds the same ones after 64 steps at N=8 as without them', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { surfaceT0, ice, stress } = buildScenario(8);
  const run = async (extra) => {
    const model = await createGpuModel(new Grid(8), { topography, ocean: { ...OCEAN_OPTIONS, penetration: 'IB', ...extra } }), gpu = model.oceanEngine;
    gpu.initialize(Float64Array.from(surfaceT0), ice);
    for (let n = 0; n < 64; n++) await gpu.advance(Float64Array.from(surfaceT0), ice, stress, 1350, undefined, (n + 1) * 1350);
    const [{ h, u, T, S, eta }, restart] = await Promise.all([gpu.download(), gpu.restartArrays()]);
    const layout = JSON.stringify(gpu.layout);
    model.destroy();
    return { layout, arrays: { h, u, T, S, eta, ...restart } };
  };
  const absent = await run({}), present = await run({ retreatDaily: false, retreatLight: 0 });
  assert.equal(present.layout, absent.layout);
  for (const [name, values] of Object.entries(absent.arrays)) assert.deepEqual(Array.from(present.arrays[name]), Array.from(values), name);
});

test('with ekmanBlock present but off, the GPU ocean lays out the same buffers, emits the same kernels and holds the same ones after 64 steps at N=8 as without it', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { surfaceT0, ice, stress } = buildScenario(8);
  const run = async (extra, steps = 64) => {
    const model = await createGpuModel(new Grid(8), { topography, ocean: { ...THIN_FLOORS, penetration: 'IB', ...extra } }), gpu = model.oceanEngine;
    gpu.initialize(Float64Array.from(surfaceT0), ice);
    for (let n = 0; n < steps; n++) await gpu.advance(Float64Array.from(surfaceT0), ice, stress, 1350, undefined, (n + 1) * 1350);
    const [{ h, u, T, S, eta }, restart] = await Promise.all([gpu.download(), gpu.restartArrays()]);
    const layout = JSON.stringify(gpu.layout);
    model.destroy();
    return { layout, arrays: { h, u, T, S, eta, ...restart } };
  };
  const absent = await run({}), present = await run({ ekmanBlock: false });
  assert.equal(present.layout, absent.layout);
  for (const [name, values] of Object.entries(absent.arrays)) assert.deepEqual(Array.from(present.arrays[name]), Array.from(values), name);
  const on = await run({ ekmanBlock: true }, 0);
  assert.notEqual(on.layout, absent.layout);
  assert.ok(JSON.parse(on.layout).OD.UEK > 0 && !('UEK' in JSON.parse(absent.layout).OD));
});

test('with upwellingSource present but zero, the GPU ocean with the Ekman block lays out the same buffers and holds the same ones after 64 steps at N=8 as without it, and with it on a coast draws its water otherwise', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { surfaceT0, ice, stress } = buildScenario(8);
  const run = async (extra, steps = 64) => {
    const model = await createGpuModel(new Grid(8), { topography, ocean: { ...THIN_FLOORS, penetration: 'IB', ekmanBlock: true, ...extra } }), gpu = model.oceanEngine;
    gpu.initialize(Float64Array.from(surfaceT0), ice);
    for (let n = 0; n < steps; n++) await gpu.advance(Float64Array.from(surfaceT0), ice, stress, 1350, undefined, (n + 1) * 1350);
    const [{ h, u, T, S, eta }, restart] = await Promise.all([gpu.download(), gpu.restartArrays()]);
    const layout = JSON.stringify(gpu.layout);
    model.destroy();
    return { layout, arrays: { h, u, T, S, eta, ...restart } };
  };
  const absent = await run({}), present = await run({ upwellingSource: 0 }), on = await run({ upwellingSource: UPWELLING });
  assert.equal(present.layout, absent.layout);
  assert.equal(on.layout, absent.layout);
  for (const [name, values] of Object.entries(absent.arrays)) assert.deepEqual(Array.from(present.arrays[name]), Array.from(values), name);
  assert.notDeepEqual(Array.from(on.arrays.h), Array.from(absent.arrays.h));
  await assert.rejects(createGpuModel(new Grid(8), { topography, ocean: { ...THIN_FLOORS, upwellingSource: 1 } }), /upwellingSource/);
});

test('the GPU ocean takes an 8-layer state that carries no class list, as the page\'s saved runs are, onto the class list as the CPU ocean does', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { cpuModel, surfaceT0, ice, stress } = buildScenario(8);
  const mesh = cpuModel.mesh, C = mesh.nCells;
  const unlisted = createCpuLayeredOcean(mesh, { geography: cpuModel.geography, ...OCEAN_OPTIONS, ...UNLISTED_OCEANS[0] });
  unlisted.initialize(Float64Array.from(surfaceT0), ice);
  for (let n = 0; n < 5; n++) unlisted.advance(Float64Array.from(surfaceT0), ice, new Float64Array(C), stress, 1350);
  const { densities, ...saved } = unlisted.serialize();
  const cpuOcean = createCpuLayeredOcean(mesh, { geography: cpuModel.geography, ...OCEAN_OPTIONS });
  cpuOcean.load(saved, Float64Array.from(surfaceT0), ice);
  const gpuModel = await createGpuModel(new Grid(8), { topography, ocean: OCEAN_OPTIONS });
  const gpuOcean = gpuModel.oceanEngine;
  await gpuOcean.upload(saved, Float64Array.from(surfaceT0), ice);
  const cpu = cpuOcean.serialize(), gpu = await gpuOcean.serialize();
  assert.equal(gpuOcean.layers, LAYER_DENSITIES.length + 1);
  assert.deepEqual(gpu.densities, LAYER_DENSITIES);
  const wet = Array.from(cpu.h, (v) => v > 1);
  for (const f of ['h', 'u', 'T', 'S', 'eta']) {
    const r = f === 'T' || f === 'S' ? stats(cpu[f].filter((_, x) => wet[x]), gpu[f].filter((_, x) => wet[x])) : stats(cpu[f], gpu[f]);
    assert.ok(r.rmsRel < 1e-6, `${f} rms relative diff ${r.rmsRel}`);
  }
  await gpuOcean.advance(Float64Array.from(surfaceT0), ice, stress, 1350);
  for (const v of (await gpuOcean.download()).h) assert.ok(Number.isFinite(v));
  gpuModel.destroy();
});

test('on the two lighter classes, a 44-class GPU state with its restart arrays loads onto both engines as a restart, and one and twenty-one GPU ocean steps from it track the CPU layered ocean at N=8', { skip: !gpuAvailable && 'webgpu not installed' }, () => lighterClasses({}, ''));

test('on the two lighter classes under the mixing package, the sunlight and a 20 m floor, a 44-class GPU state loads onto both engines and one and twenty-one GPU ocean steps from it track the CPU layered ocean at N=8', { skip: !gpuAvailable && 'webgpu not installed' }, () => lighterClasses({ minimumThickness: 20, shallowestMixedDepth: 20, penetration: 'IB', momentumDepth: 50, carryMomentum: true, refillRate: 0.00027778, refillFloor: 5, shearEntrainment: 0.65, ekmanDepth: 0.4, momentumFloor: 10, stressThickness: 20, retreatStirring: 2.5, langmuir: 2 }, ' under the package'));

async function lighterClasses(extra, name) {
  const { cpuModel, surfaceT0, ice, stress } = buildScenario(8);
  const mesh = cpuModel.mesh, C = mesh.nCells, dt = 1350, light = { ...OCEAN_OPTIONS, ...extra };
  const fortyFour = { densities: LAYER_DENSITIES.slice(2), salinities: LAYER_SALINITIES.slice(2), bottoms: LAYER_BOTTOMS.slice(2) };
  const before = await createGpuModel(new Grid(8), { topography, ocean: { ...OCEAN_OPTIONS, ...fortyFour } });
  before.oceanEngine.initialize(Float64Array.from(surfaceT0), ice);
  for (let n = 0; n < 5; n++) await before.oceanEngine.advance(Float64Array.from(surfaceT0), ice, stress, dt);
  const saved = { ...(await before.oceanEngine.serialize()), ...(await before.oceanEngine.restartArrays()) };
  before.destroy();
  assert.deepEqual(saved.densities, fortyFour.densities);
  const gpuModel = await createGpuModel(new Grid(8), { topography, ocean: light });
  const gpuOcean = gpuModel.oceanEngine;
  const cpuOcean = createCpuLayeredOcean(mesh, { geography: cpuModel.geography, ...light });
  await gpuOcean.upload(saved, Float64Array.from(surfaceT0), ice);
  cpuOcean.load(saved, Float64Array.from(surfaceT0), ice);
  assert.equal(gpuOcean.layers, LAYER_DENSITIES.length + 1);
  const restart = await gpuOcean.restartArrays();
  for (const name of ['previousT0', 'fresh', 'buoyancy', 'capacity']) assert.deepEqual(Array.from(restart[name]), Array.from(saved[name]), `the GPU's ${name}`);
  for (const name of ['fresh', 'capacity']) assert.deepEqual(Array.from(cpuOcean[name]), Array.from(saved[name]), `the CPU's ${name}`);
  const compare = async (limit, label) => {
    const cpu = cpuOcean.serialize(), gpu = await gpuOcean.serialize();
    assert.deepEqual(gpu.densities, LAYER_DENSITIES);
    const wet = Array.from(cpu.h, (v) => v > 1), lines = [];
    for (const f of ['h', 'u', 'T', 'S', 'eta']) {
      const r = f === 'T' || f === 'S' ? stats(cpu[f].filter((_, x) => wet[x]), gpu[f].filter((_, x) => wet[x])) : stats(cpu[f], gpu[f]);
      lines.push(`${f} ${r.rmsRel.toExponential(2)}`);
      assert.ok(r.rmsRel < limit, `${label}: ${f} rms relative diff ${r.rmsRel} at index ${r.at}`);
    }
    console.log(`46 classes${name}, ${label}: the engines' rms relative difference ${lines.join(', ')}`);
  };
  await compare(1e-6, 'as loaded');
  for (const steps of [1, 20]) {
    for (let n = 0; n < steps; n++) {
      cpuOcean.advance(Float64Array.from(surfaceT0), ice, new Float64Array(C), stress, dt);
      await gpuOcean.advance(Float64Array.from(surfaceT0), ice, stress, dt);
    }
    await compare(2e-3, steps === 1 ? 'after 1 step' : 'after 21 steps');
  }
  for (const v of (await gpuOcean.download()).h) assert.ok(Number.isFinite(v));
  assert.equal((await gpuOcean.diagnostics()).oceanLimited, 0);
  gpuModel.destroy();
}

test('the GPU ocean salts and freshens its mixed layer by the ice volume, cover times thickness, as the CPU one does, with the cover from the host or from the atmosphere\'s buffers, and one restarted between the two steps under the same cover takes the second as it does', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { cpuModel, cpuOcean, surfaceT0, ice, stress } = buildScenario(8);
  const mesh = cpuModel.mesh, C = mesh.nCells, land = cpuModel.geography.land, dt = 1350;
  const gpuModel = await createGpuModel(new Grid(8), { topography, ocean: OCEAN_OPTIONS });
  const gpuOcean = gpuModel.oceanEngine;
  const start = Float64Array.from(ice, (h) => (h > 0 ? 0.8 : 0));
  cpuOcean.initialize(surfaceT0, ice, { concentration: start });
  gpuOcean.initialize(surfaceT0, ice, { concentration: start });
  const grown = Float64Array.from(ice), spread = Float64Array.from(start), changed = [];
  for (let i = 0; i < C; i++) {
    if (land[i]) continue;
    if (ice[i] > 0) { grown[i] = ice[i] + 0.2; spread[i] = 0.9; changed.push(i); }
    else if (Math.abs(mesh.latCell[i]) > 50 * Math.PI / 180) { grown[i] = 0.3; spread[i] = 0.05; changed.push(i); }
  }
  const melted = Float64Array.from(grown, (h, i) => (ice[i] > 0 ? h : 0)), shrunk = Float64Array.from(spread, (a, i) => (ice[i] > 0 ? 0.5 : 0));
  const volume = (h, a) => Float64Array.from(h, (x, i) => (x > 0 ? x * a[i] : 0));
  const salinity = (state) => state.S.slice(0, C);
  async function compare(label, from, to, before) {
    const cpu = salinity(cpuOcean.serialize()), gpuS = salinity(await gpuOcean.serialize());
    const previous = (await gpuOcean.restartArrays()).previousIceVolume;
    let worst = 0, apart = Infinity;
    for (const i of changed) {
      worst = Math.max(worst, Math.abs(gpuS[i] - cpu[i]));
      const depth = cpuOcean.h[i], byThickness = (before[i] - 5) * 0.917 * (to.h[i] - from.h[i]) / depth, byVolume = (before[i] - 5) * 0.917 * (to.v[i] - from.v[i]) / depth;
      apart = Math.min(apart, Math.abs(byThickness - byVolume));
    }
    let remembered = 0;
    for (let i = 0; i < C; i++) if (!land[i]) remembered = Math.max(remembered, Math.abs(previous[i] - to.v[i]));
    console.log(`${label}: the engines' mixed-layer salinities differ by up to ${worst.toExponential(2)} psu in the ${changed.length} cells whose ice changed, where the thickness rule would differ from the volume rule by at least ${apart.toFixed(4)}; the GPU remembers the ice volume to ${remembered.toExponential(1)} m`);
    assert.ok(apart > 0.01 && worst < 5e-4, `${label}: salinity differs by ${worst} psu, the rules by ${apart}`);
    assert.ok(remembered < 1e-6, `${label}: the GPU remembers an ice volume off by ${remembered} m`);
    return cpu;
  }
  const first = { h: ice, v: volume(ice, start) }, second = { h: grown, v: volume(grown, spread) }, third = { h: melted, v: volume(melted, shrunk) };
  const before = salinity(cpuOcean.serialize());
  cpuOcean.advance(Float64Array.from(surfaceT0), grown, new Float64Array(C), stress, dt, spread);
  await gpuOcean.advance(Float64Array.from(surfaceT0), grown, stress, dt, spread);
  const middle = await compare('ice grown and spread, cover from the host', first, second, before);
  const saved = { ...(await gpuOcean.serialize()), ...(await gpuOcean.restartArrays()) };
  const restarted = await createGpuModel(new Grid(8), { topography, ocean: OCEAN_OPTIONS });
  await restarted.oceanEngine.upload(saved, Float64Array.from(surfaceT0), grown, spread);
  const loaded = await restarted.oceanEngine.restartArrays();
  assert.ok(saved.buoyancy.some((b) => b > 0), 'the first step left a buoyancy loss to carry');
  for (const name of ['previousIceVolume', 'buoyancy']) assert.ok(stats(saved[name], loaded[name]).rmsRel < 1e-6, `the restarted GPU ocean holds the ${name} it was given`);
  async function fromAtmosphere({ gpu, oceanEngine: ocean }) {
    gpu.uploadSurfaceTemperature(Float64Array.from(surfaceT0));
    gpu.device.queue.writeBuffer(gpu.buffers.S, 4 * gpu.layout.S.ICE, Float32Array.from(melted));
    gpu.uploadIce(shrunk);
    ocean.setStress(stress, melted, shrunk);
    ocean.readSurfaceFromAtmosphere();
    ocean.step(dt);
    ocean.mixedLayer(dt);
    ocean.salt(dt);
    ocean.writeSurface(dt);
    await gpu.device.queue.onSubmittedWorkDone();
  }
  cpuOcean.advance(Float64Array.from(surfaceT0), melted, new Float64Array(C), stress, dt, shrunk);
  await fromAtmosphere(gpuModel);
  await fromAtmosphere(restarted);
  await compare('ice melted back and narrowed, cover from the atmosphere\'s buffers', second, third, middle);
  const ran = { ...(await gpuOcean.serialize()), ...(await gpuOcean.restartArrays()) }, resumed = { ...(await restarted.oceanEngine.serialize()), ...(await restarted.oceanEngine.restartArrays()) };
  const wet = Array.from(ran.h, (v) => v > 1);
  const gaps = ['h', 'u', 'T', 'S', 'eta', 'previousT0', 'previousIceVolume', 'buoyancy', 'flux', 'capacity'].map((f) => [f, (f === 'T' || f === 'S' ? stats(ran[f].filter((_, x) => wet[x]), resumed[f].filter((_, x) => wet[x])) : stats(ran[f], resumed[f])).rmsRel]);
  let worst = 0;
  for (let i = 0; i < C; i++) if (!land[i]) worst = Math.max(worst, Math.abs(ran.S[i] - resumed.S[i]));
  console.log(`restarted after the first step under the same cover, the GPU ocean takes the second as the one that ran on does: mixed-layer salinity within ${worst.toExponential(1)} psu, rms relative ${gaps.map(([f, r]) => `${f} ${r.toExponential(1)}`).join(', ')}`);
  assert.ok(worst < 1e-5, `the restarted GPU ocean's mixed-layer salinity differs by up to ${worst} psu`);
  for (const [f, r] of gaps) assert.ok(r < 1e-6, `${f} rms relative diff ${r} between the GPU ocean that ran on and the restarted one`);
  restarted.destroy();
  gpuModel.destroy();
});

test('the GPU ocean starts and loads under partly covered ice without salting it, taking the cover from the sea ice it is given, as the model gives it', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const gpuModel = await createGpuModel(new Grid(4), { ocean: { everySteps: 1, thermoclineTilt: 0, salinityProfile: () => 34 } });
  const { mesh, seaIce } = gpuModel, ocean = gpuModel.oceanEngine, n = mesh.nCells, dt = 3600;
  const surfaceT = new Float64Array(n).fill(FREEZING_POINT), ice = new Float64Array(n).fill(1), calm = new Float64Array(mesh.nEdges);
  seaIce.load(ice, new Float64Array(n).fill(0.4));
  const salt = async () => {
    const [{ Wref }, { h }] = await Promise.all([ocean.restartArrays(), ocean.download()]);
    return Float64Array.from({ length: n }, (_, i) => Wref[i] + SALT_REFERENCE * h[i]);
  };
  const unchanged = async (label) => {
    const before = await salt();
    await ocean.advance(Float64Array.from(surfaceT), ice, calm, dt);
    const after = await salt();
    for (let i = 0; i < n; i++) assert.ok(Math.abs(after[i] - before[i]) < 1e-5 * before[i], `${label}: cell ${i} salt ${before[i]} -> ${after[i]}`);
  };
  ocean.initialize(Float64Array.from(surfaceT), ice);
  await unchanged('after initialize');
  await ocean.upload(await ocean.serialize(), Float64Array.from(surfaceT), ice);
  await unchanged('after upload');
  await ocean.upload(await ocean.serialize(), Float64Array.from(surfaceT), ice, new Float64Array(n).fill(1));
  const before = await salt();
  await ocean.advance(Float64Array.from(surfaceT), ice, calm, dt);
  const after = await salt();
  assert.ok(after[0] < before[0] - 10, `an ocean uploaded under full cover books the cover falling to 0.4 as melt: ${before[0]} -> ${after[0]}`);
  gpuModel.destroy();
});
