import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createOcean as createCpuLayeredOcean, LAYER_DENSITIES } from '../js/ocean/layered.module.js';
import { UNLISTED_OCEANS } from './helpers/layered.mjs';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { labelTemperature } from '../js/ocean/seawater.module.js';

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

test('under closureFill the GPU closure acts on the class flow carried onto token edges as the CPU\'s does, through twenty-one ocean steps at N=8', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const filling = { ...OCEAN_OPTIONS, closureFill: 0.5, closureHours: 1 };
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

test('the GPU mixed layer holds, retreats, convects, keeps to its neighbours and, under the old options, snaps back as the CPU one does', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { cpuModel, surfaceT0, ice, stress } = buildScenario(8);
  const mesh = cpuModel.mesh, C = mesh.nCells;
  const sea = [...Array(C).keys()].filter((i) => cpuModel.geography.land[i] === 0);
  const neutral = sea.filter((_, n) => n % 5 === 0), stable = sea.filter((_, n) => n % 5 === 1), dense = sea.filter((_, n) => n % 5 === 2);
  const oldRules = { neutralSnap: true, convectiveErosion: false, buoyancyMemory: 0, maximumMixedDepth: 200, mixedNeighbourRatio: 0, vorticityCentring: 0 };
  for (const options of [{ ...OCEAN_OPTIONS, mixedNeighbourRatio: 3 }, { ...OCEAN_OPTIONS, ...oldRules }]) {
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
    const saved = cpuOcean.serialize();
    cpuOcean.load(saved, surfaceT, ice);
    const gpuModel = await createGpuModel(new Grid(8), { topography, ocean: options });
    const gpuOcean = gpuModel.oceanEngine;
    await gpuOcean.upload(saved, Float64Array.from(surfaceT), ice);
    const warmed = Float64Array.from(surfaceT, (t, i) => (neutral.includes(i) ? t + 0.05 : t));
    cpuOcean.advance(Float64Array.from(warmed), ice, new Float64Array(C), stress, 1350);
    await gpuOcean.advance(Float64Array.from(warmed), ice, stress, 1350);
    const cpu = cpuOcean.serialize(), gpu = await gpuOcean.serialize();
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
