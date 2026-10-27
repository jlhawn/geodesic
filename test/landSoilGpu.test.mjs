import { test } from 'node:test';
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { createLandSurface, soilColumn } from '../js/physics/land.module.js';
import { decisionTracker, cpuDecisions, gpuDecisions, neighbourhood } from './helpers/decisions.mjs';
import { PREVIOUS_LEVELS } from './helpers/previousDefaults.mjs';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};
const { layoutFor } = gpuAvailable ? await import('../js/gpu/core.gpu.js') : {};
const { physicsSurfaceKernel } = gpuAvailable ? await import('../js/gpu/physics.gpu.js') : {};

const DT = 900, STEPS = 8;
const topography = syntheticTopography(90, 180, (lat, lon) => ((Math.cos(lon) > 0 && Math.abs(lat) < 1.2) || lat < -1.15 ? 300 : -4000));
const layers = { iceSheetBalance: true };
const control = { soil: 'slab', heatCapacity: 2.19e6 * 0.07, iceSheetBalance: true };

function prepare(model) {
  const { mesh, geography } = model, C = mesh.nCells;
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  for (let i = 0; i < C; i++) if (geography.land[i]) {
    model.state[6][i] = 0;
    if (geography.iceSheet[i] && i % 2 === 0) model.state[3][i] = 283;
    if (!geography.iceSheet[i] && mesh.latCell[i] > 1.0 && i % 2 === 0) model.state[3][i] = 276;
  }
  if (model.load) model.load();
  model.ocean.initialize(model.state[3], model.state[6]);
  model.land.initialize();
  for (let i = 0; i < C; i++) if (geography.land[i]) { model.land.soil[i] = 40; model.land.snow[i] = mesh.latCell[i] > 1.0 ? 5 : 0; }
  model.land.load({ soil: Float64Array.from(model.land.soil), snow: Float64Array.from(model.land.snow), vegetation: new Float64Array(C).fill(0.5) }, undefined, 0, model.state[3]);
  return model;
}

const stats = (a, b, skip = null) => {
  let maxDiff = 0, at = -1;
  for (let x = 0; x < a.length; x++) { if (skip && skip.has(x)) continue; const d = Math.abs(a[x] - b[x]); if (d > maxDiff) { maxDiff = d; at = x; } }
  return { maxDiff, at };
};

const pairs = new Map();
const engines = (name, land) => {
  if (!pairs.has(name)) pairs.set(name, (async () => ({
    cpu: prepare(createModel(new Grid(16), { topography, land, levels: PREVIOUS_LEVELS })),
    gpu: prepare(await createGpuModel(new Grid(16), { topography, land, levels: PREVIOUS_LEVELS })),
  }))());
  return pairs.get(name);
};

async function run({ cpu, gpu }) {
  const decisions = decisionTracker(cpu.core.K, cpu.mesh.nCells);
  for (let n = 0; n < STEPS; n++) {
    cpu.step(DT); await gpu.step(DT);
    decisions.check(cpuDecisions(cpu), gpuDecisions(await gpu.gpu.download(), await gpu.gpu.downloadPhysics()));
  }
  await gpu.sync();
  return { decisions, near: neighbourhood(cpu.mesh, decisions.parted, 1), saved: await gpu.land.serialize() };
}

async function fedStep(gpu, land) {
  const { geography } = gpu, C = gpu.mesh.nCells;
  await gpu.sync();
  const before = await gpu.land.serialize(), skinBefore = Float64Array.from(gpu.state[3]), lossBefore = Float64Array.from(gpu.land.iceLoss);
  await gpu.step(DT);
  await gpu.sync();
  const after = await gpu.land.serialize(), ph = await gpu.gpu.downloadPhysics();
  const check = createLandSurface(gpu.mesh, geography, land), T = Float64Array.from(skinBefore);
  check.load(before, null, 0);
  const fed = { skin: 0, atSkin: -1, deep: 0, atDeep: -1 };
  for (let i = 0; i < C; i++) if (geography.land[i]) check.update(i, T, ph.SFLUX, ph.EVAP[i], DT);
  for (let i = 0; i < C; i++) {
    if (!geography.land[i]) continue;
    const ds = Math.abs(T[i] - gpu.state[3][i]);
    if (ds > fed.skin) { fed.skin = ds; fed.atSkin = i; }
    for (const name of gpu.land.deepNames) { const dd = Math.abs(check[name][i] - after[name][i]); if (dd > fed.deep) { fed.deep = dd; fed.atDeep = i; } }
  }
  return { before, after, ph, skinBefore, lossBefore, fed };
}

test('eight GPU steps with the layered soil track the CPU model as closely as a slab of the skin\'s capacity does, and the three deeper layers agree to 3e-3 K at every land cell outside the decisions the engines took apart', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { cpu, gpu } = await engines('layers', layers), plain = await engines('control', control);
  const a = await run({ cpu, gpu }), b = await run(plain);
  const C = cpu.mesh.nCells, { geography } = cpu;
  assert.deepEqual([cpu.land.layered, gpu.land.layered, gpu.gpu.physics.soilColumn.depths.length, plain.gpu.gpu.physics.soilColumn], [true, true, 4, null], 'the layers are the default and the slab is asked for');
  const skin = stats(cpu.state[3], gpu.state[3], a.near), reference = stats(plain.cpu.state[3], plain.gpu.state[3], b.near);
  const deep = gpu.land.deepNames.map((name) => ({ name, ...stats(cpu.land[name], a.saved[name], a.near), all: stats(cpu.land[name], a.saved[name]) }));
  const sheets = [...geography.iceSheet].filter((v, i) => v && geography.land[i] && !a.near.has(i)).length, melted = Array.from(gpu.land.iceLoss).filter((v) => v > 0).length;
  console.log(`layered soil, ${STEPS} steps at N=16: skin max ${skin.maxDiff.toExponential(1)} K (a slab of the skin's capacity: ${reference.maxDiff.toExponential(1)} K) outside the ${a.near.size} cells within one of the ${a.decisions.parted.size} where the engines took a discrete decision apart; ${deep.map((l) => `${l.name} ${l.maxDiff.toExponential(1)} K`).join(', ')}; all cells ${deep.map((l) => l.all.maxDiff.toExponential(1)).join(', ')}; ${sheets} ice-sheet cells compared, ${melted} lost glacier ice, snow max ${stats(cpu.land.snow, a.saved.snow).maxDiff.toExponential(1)} kg/m²`);
  assert.ok(a.near.size <= 0.5 * C, `${a.near.size} cells left out`);
  assert.ok(skin.maxDiff < Math.max(0.02, 1.5 * reference.maxDiff), `skin max ${skin.maxDiff} at ${skin.at}, a slab of the same capacity ${reference.maxDiff}`);
  for (const l of deep) assert.ok(l.maxDiff < 3e-3, `${l.name} max ${l.maxDiff} at ${l.at}`);
  assert.ok(sheets > 0 && melted > 0, `${sheets} ice-sheet cells, ${melted} with glacier melt`);
  let moved = 0;
  for (let i = 0; i < C; i++) if (geography.land[i]) moved = Math.max(moved, Math.abs(cpu.land.soilT2[i] - cpu.land.soilT4[i]));
  assert.ok(moved > 0.5, 'the layers differ from one another, so the conduction is active');
  assert.ok(stats(cpu.land.snow, a.saved.snow).maxDiff < 1e-2 && stats(cpu.land.iceLoss, gpu.land.iceLoss).maxDiff < 1e-2);
  const ph = await gpu.gpu.downloadPhysics();
  assert.ok(ph.SOILT2 && ph.SOILT3 && ph.SOILT4 && !ph.SOILT5);
});

test('the WGSL for one, two and five layers gives the CPU land\'s skin and layers from the GPU\'s own fluxes, in one step and in the next', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  for (const soilDepths of [[0.2], [0.1, 0.3], [0.05, 0.1, 0.2, 0.4, 0.8]]) {
    const land = { soilDepths, iceSheetBalance: true }, gpu = prepare(await createGpuModel(new Grid(16), { topography, land, levels: PREVIOUS_LEVELS }));
    const rows = [];
    for (let step = 0; step < 2; step++) rows.push((await fedStep(gpu, land)).fed);
    console.log(`${soilDepths.length} layers: ${rows.map((fed) => `skin ${fed.skin.toExponential(1)} K, layers ${fed.deep.toExponential(1)} K`).join('; ')}`);
    assert.equal(gpu.land.deepNames.length, soilDepths.length - 1);
    for (const fed of rows) assert.ok(fed.skin < 2e-3 && fed.deep < 2e-3, `${soilDepths.length} layers: skin ${fed.skin} at ${fed.atSkin}, layers ${fed.deep} at ${fed.atDeep}`);
  }
});

test('one GPU step, fed to the CPU land with the fluxes the GPU used, gives the same skin and deep layers, and the layers\' heat content changes by the net surface flux times the step wherever the snow and the ice are unchanged', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { gpu } = await engines('layers', layers);
  const { geography } = gpu, C = gpu.mesh.nCells, column = gpu.land.column;
  const heat = (skin, arrays, i) => (geography.iceSheet[i] ? column.firn : column.soil).capacity.reduce((sum, c, k) => sum + c * (k === 0 ? skin[i] : arrays[k - 1][i]), 0);
  const { before, after, ph, skinBefore, lossBefore, fed } = await fedStep(gpu, layers);
  let compared = 0, worst = 0, worstAt = -1, sheets = 0, changed = 0;
  for (let i = 0; i < C; i++) {
    if (!geography.land[i]) continue;
    if (after.snow[i] !== before.snow[i] || gpu.land.iceLoss[i] !== lossBefore[i]) { changed++; continue; }
    const h0 = heat(skinBefore, gpu.land.deepNames.map((n) => before[n]), i), h1 = heat(gpu.state[3], gpu.land.deepNames.map((n) => after[n]), i);
    const residual = Math.abs(h1 - h0 - ph.SFLUX[i] * DT) / (Math.abs(h0) + Math.abs(h1));
    if (residual > worst) { worst = residual; worstAt = i; }
    compared++; if (geography.iceSheet[i]) sheets++;
  }
  console.log(`one step with the GPU's fluxes on the CPU land: skin ${fed.skin.toExponential(1)} K (cell ${fed.atSkin}), deep layers ${fed.deep.toExponential(1)} K (cell ${fed.atDeep}); GPU land energy row: ${compared} land cells (${sheets} on ice sheets) with the snow and ice unchanged, worst residual ${worst.toExponential(1)} of the column's heat content (cell ${worstAt}); ${changed} cells with melting or fresh snow left out`);
  assert.ok(fed.skin < 2e-3 && fed.deep < 2e-3, `skin ${fed.skin} at ${fed.atSkin}, deep ${fed.deep} at ${fed.atDeep}`);
  assert.ok(compared > 20 && sheets > 0);
  assert.ok(worst < 1e-6, `residual ${worst} at ${worstAt}`);
});

test('for soil \'slab\' the surface kernel\'s text, its constants and the PH layout are the slab\'s', { skip: !gpuAvailable && 'webgpu not installed' }, () => {
  const sha = (text) => createHash('sha256').update(text).digest('hex').slice(0, 16);
  const mesh = createModel(new Grid(16), { physics: false }).mesh;
  assert.equal(sha(physicsSurfaceKernel(false)), 'ee8de3a466192bb4');
  assert.equal(sha(physicsSurfaceKernel(true, true)), 'f3367a8e9c2b8972');
  const plain = layoutFor(mesh, 8), tall = layoutFor(mesh, 8, 1, 1, 1, null, true, 1, 1, 8, true, true, true, true);
  assert.deepEqual([plain.PH.total, sha(JSON.stringify(plain.PH))], [389407, '1f1a4f4766b1ccdd']);
  assert.deepEqual([tall.PH.total, sha(JSON.stringify(tall.PH))], [589251, 'bb66998c08f55b4e']);
  const column = soilColumn({ soil: 'layers' }), text = physicsSurfaceKernel(false, false, column);
  assert.ok(!/colC1|SOILT|SOIL_C/.test(physicsSurfaceKernel(false)));
  assert.ok(/PH_SOILT4/.test(text) && !/PH_SOILT5/.test(text) && /FIRN_G3/.test(text) && !text.includes('T += dt * net / LANDC'));
  const layered = layoutFor(mesh, 8, 0, 0, 0, null, false, 0, 0, false, false, false, false, false, 3);
  assert.equal(layered.PH.total, plain.PH.total + 3 * mesh.nCells);
  assert.deepEqual(['SOILT2', 'SOILT3', 'SOILT4'].map((name) => layered.PH[name] - plain.PH.total), [0, 1, 2].map((k) => k * mesh.nCells));
});
