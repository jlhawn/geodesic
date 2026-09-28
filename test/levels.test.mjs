import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, writeFileSync, mkdtempSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { execFileSync } from 'node:child_process';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createParallelModel } from '../js/parallel.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography, topographyFromInt16 } from '../js/geography.module.js';
import { regridState, savedDeckField, DECK_FIELDS } from '../js/physics/regrid.module.js';
import { sigmaInterfaces, sigmaGridName, standardHeight, standardSigma, SIGMA_GRIDS } from '../js/dynamics/sigmaCore.module.js';
import { encodeState, decodeState, savedLevels } from '../js/stateFile.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const cam = sigmaInterfaces('cam26'), bl = sigmaInterfaces('bl34');
const topography = syntheticTopography(90, 180, (lat, lon) => ((Math.cos(lon) > 0 && Math.abs(lat) < 1.2) || lat < -1.15 ? 300 + 1500 * Math.exp(-(((lat - 0.3) / 0.3) ** 2)) : -4000));

function stats(cpu, gpu) {
  let maxDiff = 0, at = -1, sumSq = 0, sumRef = 0;
  for (let x = 0; x < cpu.length; x++) { const d = Math.abs(cpu[x] - gpu[x]); if (d > maxDiff) { maxDiff = d; at = x; } sumSq += d * d; sumRef += cpu[x] * cpu[x]; }
  return { maxDiff, at, rmsRel: Math.sqrt(sumSq / Math.max(sumRef, 1e-300)) };
}

function prepare(model) {
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i]) model.state[6][i] = 0;
  if (model.load) model.load();
  model.ocean.initialize(model.state[3], model.state[6]);
  model.land.initialize();
  return model;
}

test('the sigma grids go by name, cam26 by default, and a grid is named again from its interfaces in single precision', () => {
  assert.deepEqual(SIGMA_GRIDS, ['cam26', 'bl34']);
  assert.deepEqual(sigmaInterfaces(), cam);
  assert.equal(cam.length - 1, 27);
  assert.equal(bl.length - 1, 34);
  assert.throws(() => sigmaInterfaces('cam30'), /no sigma grid is named cam30/);
  assert.equal(sigmaGridName(cam), 'cam26');
  assert.equal(sigmaGridName(Float32Array.from(bl)), 'bl34');
  assert.equal(sigmaGridName(Float64Array.from(bl, (sigma, k) => (k === 30 ? sigma + 1e-4 : sigma))), null);
  assert.equal(sigmaGridName(bl.subarray(1)), null);
});

test('bl34 is monotone from 0 to 1 with a 40 m lowest layer, ten layers below 1.25 km, and below 2.5 km each layer at most 1.5 times as thick as the one under it, in height and in mass', () => {
  const K = bl.length - 1;
  assert.equal(bl[0], 0);
  assert.equal(bl[K], 1);
  for (let k = 1; k <= K; k++) assert.ok(bl[k] > bl[k - 1], `interface ${k}`);
  assert.ok(Math.abs(standardHeight(bl[K - 1]) - 40) < 0.5, `lowest layer ${standardHeight(bl[K - 1])} m`);
  assert.equal([...bl].filter((sigma) => sigma < 1 && standardHeight(sigma) < 1250).length, 10);
  const first = bl.findIndex((sigma) => sigma > standardSigma(2500));
  const ratios = [];
  for (let k = first - 1; k < K - 1; k++) {
    const mass = (bl[k + 1] - bl[k]) / (bl[k + 2] - bl[k + 1]);
    const height = (standardHeight(bl[k]) - standardHeight(bl[k + 1])) / (standardHeight(bl[k + 1]) - standardHeight(bl[k + 2]));
    ratios.push(height.toFixed(3));
    assert.ok(mass >= 1 && mass <= 1.5, `layer ${k} over layer ${k + 1}: ${mass} in mass`);
    assert.ok(height >= 1 && height <= 1.5, `layer ${k} over layer ${k + 1}: ${height} in height`);
  }
  console.log(`bl34 interfaces below 2.5 km: ${[...bl].slice(first).map((sigma) => `${sigma.toFixed(4)} (${standardHeight(sigma).toFixed(0)} m)`).join(', ')}; thickness ratios upward ${ratios.reverse().join(' ')}`);
});

test('above 2.5 km bl34 is cam26 exactly, interfaces and layers, down to the 2.4 km interface they share', () => {
  const shared = cam.findIndex((sigma) => standardHeight(sigma) < 2500);
  assert.equal(shared, 22);
  assert.ok(Math.abs(standardHeight(cam[shared]) - 2421) < 1);
  assert.deepEqual(bl.subarray(0, shared + 1), cam.subarray(0, shared + 1));
  const a = createModel(new Grid(2), { physics: false }), b = createModel(new Grid(2), { physics: false, levels: bl });
  for (let k = 0; k < shared; k++) {
    assert.equal(b.core.sigmaMid[k], a.core.sigmaMid[k]);
    assert.equal(b.core.diagnostics.dSigma[k], a.core.diagnostics.dSigma[k]);
  }
  assert.equal(a.core.K, 27);
  assert.equal(b.core.K, 34);
});

test('a state saved on bl34 reloads on bl34, with the deck\'s carried state, and steps as the run that saved it does; a state without levels loads on cam26', async () => {
  const model = createModel(new Grid(4), { levels: bl, ocean: false });
  const init = initializeState(model, {});
  init.forEach((values, a) => model.state[a].set(values.map(Math.fround)));
  const { mlmSubsidence, mlmHeight, mlmGate } = model.radiation;
  for (let i = 0; i < mlmHeight.length; i++) {
    mlmSubsidence[i] = Math.fround(-0.004 + 0.001 * (i % 9));
    mlmHeight[i] = Math.fround(600 + 7 * (i % 80));
    mlmGate[i] = Math.fround(0.1 + 0.1 * (i % 9));
  }
  const [pi, theta, u, surfaceT, q, qc, ice] = model.state;
  const bytes = encodeState({ N: 4, K: model.core.K, day: 0, time: 0, terrain: false, levels: model.core.levels, pi, theta, u, surfaceT, q, qc, ice, mlmSubsidence, mlmHeight, mlmGate });
  const saved = await decodeState(bytes);
  assert.ok(saved.levels instanceof Float64Array, 'the interfaces are kept in double precision');
  assert.deepEqual(savedLevels(saved), bl);
  const reloaded = createModel(new Grid(saved.N), { levels: savedLevels(saved), ocean: false });
  assert.equal(reloaded.core.K, 34);
  ['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'].forEach((name, a) => reloaded.state[a].set(saved[name]));
  for (const name of Object.keys(DECK_FIELDS)) {
    reloaded.radiation[name].set(savedDeckField(saved, name, reloaded));
    assert.deepEqual(reloaded.radiation[name], model.radiation[name], `${name} reloads as saved`);
  }
  reloaded.seaIce.load(reloaded.state[6]);
  model.seaIce.load(model.state[6]);
  for (let n = 0; n < 3; n++) { model.step(1350 * 4); reloaded.step(1350 * 4); }
  model.state.forEach((array, a) => assert.deepEqual(reloaded.state[a], array, `state array ${a} after three steps`));
  for (const name of Object.keys(DECK_FIELDS)) assert.deepEqual(reloaded.radiation[name], model.radiation[name], `${name} after three steps`);

  const legacy = await decodeState(encodeState({ N: 4, K: 27, day: 0, time: 0, pi }));
  assert.equal(legacy.levels, undefined);
  assert.deepEqual(savedLevels(legacy), cam);
  assert.deepEqual(savedLevels({ levels: Array.from(bl, Math.fround) }), bl, 'single-precision interfaces load the exact grid they name');
  const custom = Float64Array.from(cam, (sigma, k) => (k === 20 ? sigma + 0.001 : sigma));
  assert.deepEqual(savedLevels((await decodeState(encodeState({ N: 4, levels: custom })))), custom);
});

test('a bl34 state regrids across resolutions onto bl34 and not onto cam26', () => {
  const coarse = createModel(new Grid(4), { levels: bl, physics: false }), fine = createModel(new Grid(6), { levels: bl, physics: false }), other = createModel(new Grid(6), { physics: false });
  const K = coarse.core.K, C = coarse.mesh.nCells, E = coarse.mesh.nEdges;
  const pi = new Float64Array(C).fill(1e5), theta = Float64Array.from({ length: K * C }, (_, x) => 280 + 100 * (1 - coarse.core.sigmaMid[Math.floor(x / C)]));
  const [fPi, fTheta, fU] = regridState(coarse, fine, [pi, theta, new Float64Array(K * E), new Float64Array(C).fill(288)]);
  assert.equal(fTheta.length, 34 * fine.mesh.nCells);
  assert.equal(fU.length, 34 * fine.mesh.nEdges);
  for (let k = 0; k < K; k++) assert.ok(Math.abs(fTheta[k * fine.mesh.nCells] - theta[k * C]) < 1e-9);
  assert.ok(fPi.every((p) => Math.abs(p - 1e5) < 1e-6));
  assert.throws(() => regridState(coarse, other, [pi, theta, new Float64Array(K * E), new Float64Array(C)]), /layer counts differ: 34 vs 27/);
});

test('the worker-thread engine on bl34 reproduces the single-thread step bit for bit', async () => {
  const serial = createModel(new Grid(6), { levels: bl });
  const parallel = await createParallelModel(new Grid(6), { levels: bl }, 2);
  try {
    assert.equal(parallel.core.K, 34);
    const init = initializeState(serial, {});
    for (let a = 0; a < 4; a++) { serial.state[a].set(init[a]); parallel.state[a].set(init[a]); }
    for (let n = 0; n < 2; n++) { serial.step(900); parallel.step(900); }
    serial.state.forEach((array, a) => assert.deepEqual(parallel.state[a], array, `state array ${a}`));
  } finally {
    await parallel.close();
  }
});

test('the full model on bl34 steps alike on the CPU and the GPU, over a continent with its ocean', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const cpu = prepare(createModel(new Grid(6), { topography, levels: bl }));
  const gpu = prepare(await createGpuModel(new Grid(6), { topography, levels: bl }));
  assert.equal(gpu.gpu.K, 34);
  assert.equal(gpu.core.K, 34);
  for (let n = 0; n < 8; n++) { cpu.step(900); await gpu.step(900); }
  await gpu.sync();
  const d = await gpu.diagnostics(), dc = cpu.diagnostics();
  const ts = stats(cpu.state[3], gpu.state[3]), theta = stats(cpu.state[1], gpu.state[1]), q = stats(cpu.state[4], gpu.state[4]), u = stats(cpu.state[2], gpu.state[2]);
  console.log(`eight steps on bl34 at N=6 over a continent: Ts max ${ts.maxDiff.toExponential(1)} K, θ rms ${theta.rmsRel.toExponential(1)}, q rms ${q.rmsRel.toExponential(1)}, wind max ${u.maxDiff.toExponential(1)} m/s; mean Ts ${dc.meanSurfaceT.toFixed(3)} vs ${d.meanSurfaceT.toFixed(3)} K, OLR ${dc.outgoingLongwave.toFixed(2)} vs ${d.outgoingLongwave.toFixed(2)} W/m², ocean h1 ${dc.oceanUpperDepth.toFixed(2)} vs ${d.oceanUpperDepth.toFixed(2)} m`);
  for (const array of [...cpu.state, ...gpu.state]) assert.ok(array.every(Number.isFinite));
  assert.ok(ts.maxDiff < 0.02, `Ts max ${ts.maxDiff} at ${ts.at}`);
  assert.ok(theta.rmsRel < 1e-4, `θ rms ${theta.rmsRel}`);
  assert.ok(q.rmsRel < 2e-3, `q rms ${q.rmsRel}`);
  assert.ok(u.maxDiff < 0.05, `wind max ${u.maxDiff} at ${u.at}`);
  assert.ok(Math.abs(dc.meanSurfaceT - d.meanSurfaceT) < 0.01 && Math.abs(dc.outgoingLongwave - d.outgoingLongwave) < 0.5 && Math.abs(dc.oceanUpperDepth - d.oceanUpperDepth) < 1e-3);
  gpu.destroy();
});

test('spinup.mjs starts a bl34 run from a cam26 state\'s ocean and land with a fresh atmosphere at day 0, with fresh sea ice unless ICE_FROM, checkpoints it inside a day on bl34 and resumes it there alone, and refuses a state at another N', { skip: !gpuAvailable && 'webgpu not installed' }, async (t) => {
  const root = new URL('..', import.meta.url).pathname, N = 6, dir = mkdtempSync(join(tmpdir(), 'levels-'));
  t.after(() => rmSync(dir, { recursive: true, force: true }));
  const source = await createGpuModel(new Grid(N), { topography: topographyFromInt16(readFileSync(join(root, 'data/topography_0p25.bin')).buffer) });
  const { state, mesh } = source, C = mesh.nCells, land = source.geography.land;
  initializeState(source, { geostrophic: false }).forEach((values, a) => state[a].set(values));
  for (let i = 0; i < C; i++) if (land[i]) state[6][i] = 0;
  source.load();
  source.ocean.initialize(state[3], state[6]);
  source.land.initialize();
  const ocean = await source.ocean.serialize(), soil = await source.land.serialize();
  const [pi, theta, u, surfaceT, q, qc] = state, ice = Float64Array.from(state[6], (h) => (h > 0 ? 0.25 : 0)), concentration = Float64Array.from(ice, (h) => (h > 0 ? 0.6 : 0));
  for (let i = 0; i < C; i++) if (!land[i] && !(ice[i] > 0)) ocean.T[i] += 2;
  const vegetation = Float64Array.from(soil.vegetation, (v) => (v > 0 ? 0.8 : 0));
  const from = join(dir, 'src_day0100.bin');
  writeFileSync(from, encodeState({ N, K: source.core.K, day: 100, time: 100 * 86400, terrain: true, pi, theta, u, surfaceT, q, qc, ice, concentration, mlmSubsidence: source.radiation.mlmSubsidence, ocean: { h: ocean.h, u: ocean.u, T: ocean.T, S: ocean.S, eta: ocean.eta }, land: { ...soil, vegetation } }));
  source.destroy();
  const run = (env) => execFileSync(process.execPath, [join(root, 'scripts/spinup.mjs')], { cwd: root, env: { ...process.env, N: String(N), OUT: dir, FROM: from, DAYS: '1', MINUTES: '100', ...env }, encoding: 'utf8', stdio: 'pipe' });
  const read = async (name) => decodeState(new Uint8Array(readFileSync(join(dir, name))));

  const log = run({ TAG: 'fine', LEVELS: 'bl34' });
  assert.match(log, /seeded from .*src_day0100\.bin \(N=6, day 100, cam26\).*with fresh sea ice; a fresh atmosphere on bl34 \(34 layers\); the clock at day 0/);
  const fine = await read('fine_day0001.bin');
  assert.equal(fine.day, 1);
  assert.equal(fine.time, 86400);
  assert.equal(fine.K, 34);
  assert.deepEqual(savedLevels(fine), bl);
  let sea = 0, drift = 0, lands = 0, green = 0, freshIce = 0;
  for (let i = 0; i < C; i++) {
    if (!land[i] && !(fine.ice[i] > 0) && !(ice[i] > 0)) { sea++; drift += fine.ocean.T[i] - ocean.T[i]; }
    if (vegetation[i] > 0) { lands++; green += Math.abs(fine.land.vegetation[i] - 0.8); }
    freshIce = Math.max(freshIce, fine.ice[i]);
  }
  console.log(`seeded bl34 run after a day: mixed layer ${(drift / sea).toFixed(3)} K from FROM's, vegetation ${(green / lands).toFixed(4)} from FROM's on average, thickest ice ${freshIce.toFixed(2)} m (FROM had 0.25 m)`);
  assert.ok(Math.abs(drift / sea) < 0.5, `the mixed layer moved ${drift / sea} K from FROM's in a day`);
  assert.ok(green / lands < 0.02, `the vegetation is ${green / lands} from FROM's`);
  assert.ok(freshIce > 1, `the ice is fresh: at most ${freshIce} m`);

  run({ TAG: 'fine', DAYS: '2', STOP_AFTER_STEPS: '12' });
  const inside = await read('fine_day0001_step0012.bin');
  assert.equal(inside.step, 12);
  assert.deepEqual(savedLevels(inside), bl, 'the in-day checkpoint is on bl34');
  for (const field of Object.keys(DECK_FIELDS)) assert.equal(inside[field]?.length, C, `the in-day checkpoint carries ${field}`);
  assert.throws(() => run({ TAG: 'fine', DAYS: '2', LEVELS: 'cam26' }), /fine_day0001_step0012\.bin is on bl34 \(34 layers\), not cam26/);
  assert.match(run({ TAG: 'fine', DAYS: '2', LEVELS: 'bl34' }), /continuing from fine_day0001_step0012\.bin \(day 1, step 12 of 24\) on bl34 \(34 layers\)/);
  const resumed = await read('fine_day0002.bin');
  assert.equal(resumed.K, 34);
  assert.deepEqual(savedLevels(resumed), bl);

  assert.match(run({ TAG: 'iced', ICE_FROM: '1' }), /and its sea ice \(thickness, concentration, snow, skin temperature\); a fresh atmosphere on cam26 \(27 layers\)/);
  const iced = await read('iced_day0001.bin');
  assert.equal(iced.K, 27);
  let covered = 0, thickest = 0;
  for (let i = 0; i < C; i++) if (iced.ice[i] > 0) { covered++; thickest = Math.max(thickest, iced.ice[i]); assert.ok(Math.abs(iced.concentration[i] - 0.6) < 0.2, `concentration ${iced.concentration[i]} at ${i}`); }
  assert.ok(covered > 0 && thickest < 0.4, `${covered} iced cells, the thickest ${thickest} m`);

  assert.throws(() => run({ N: '8', TAG: 'coarse' }), /is N=6, not 8/);
});
