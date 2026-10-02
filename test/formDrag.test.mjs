import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { topographyFromInt16, meshSubgrid } from '../js/geography.module.js';
import { formDragCoefficient } from '../js/physics/formDrag.module.js';
import { sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { cellVector } from '../js/dynamics/operators.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};
const { readRanges } = gpuAvailable ? await import('../js/gpu/device.module.js') : {};

const close = (actual, expected, tolerance, what) => assert.ok(Math.abs(actual - expected) <= tolerance * Math.abs(expected), `${what}: ${actual} against ${expected}`);
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);

function prepare(model) {
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  model.seaIce.load(model.state[6]);
  if (model.load) model.load();
  model.ocean.initialize(model.state[3], model.state[6]);
  model.land.initialize();
  return model;
}

const handCoefficient = (sigma, z) => {
  const a1 = sigma * sigma / (0.00102 * 0.00035 ** -1.9), a2 = a1 * 0.003 ** (-1.9 + 2.8);
  return 35 * 1 * 0.005 * 0.6 * 2.109 * Math.exp(-((z / 1500) ** 1.5)) * a2 * z ** -1.2;
};

test('the form drag coefficient by hand: σ_flt 100 m at 20 m, 500 m and 2 km, and the column stress at 10 m/s', () => {
  for (const [z, value] of [[20, 8.6676e-5], [500, 1.5047e-6], [2000, 7.4119e-8]]) {
    close(formDragCoefficient(100, z), handCoefficient(100, z), 1e-12, `C_tofd at ${z} m`);
    close(formDragCoefficient(100, z), value, 1e-4, `C_tofd at ${z} m against the pinned value`);
  }
  close(formDragCoefficient(200, 300), 4 * formDragCoefficient(100, 300), 1e-12, 'C_tofd ∝ σ_flt²');
  assert.equal(formDragCoefficient(0, 100), 0);
  let stress = 0;
  for (let z = 20; z < 20000; z += 1) stress += 1.2 * formDragCoefficient(100, z + 0.5) * 100;
  close(stress, 0.565, 0.01, 'ρ ∫ C_tofd U² dz from 20 m at 10 m/s, σ_flt 100 m (N/m²)');
});

function uniformWind(model, speed) {
  const { mesh, core } = model, { K, E } = core.diagnostics, u = model.state[2];
  for (let e = 0; e < E; e++) {
    const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1], lon = 0.5 * (mesh.lonCell[a] + mesh.lonCell[b]);
    const east = [-Math.sin(lon), Math.cos(lon), 0], normal = speed * (east[0] * mesh.nEdge[3 * e] + east[1] * mesh.nEdge[3 * e + 1]);
    for (let k = 0; k < K; k++) u[k * E + e] = normal;
  }
}

test('a column of the edge solve against the hand profile: σ_flt 100 m on every cell, 10 m/s, no mixing and no surface drag; none over the sea', () => {
  const N = 4, mesh = buildMesh(new Grid(N)), C = mesh.nCells;
  const sigma = 100, fields = { deviation: new Float64Array(C), anisotropy: new Float64Array(C), orientation: new Float64Array(C), slope: new Float64Array(C), filtered: new Float64Array(C).fill(sigma) };
  const model = prepare(createModel(mesh, { topography, levels: sigmaInterfaces('bl34'), subgrid: fields }));
  const { core, boundaryLayer: bl, geography } = model, { K, E, dSigma, g, geopotential } = core.diagnostics, dt = 1350 * 16 / N;
  uniformWind(model, 10);
  bl.diagnose(model.state);
  const vector = new Float64Array(3 * C);
  let checked = 0;
  for (let k = bl.kTop; k < K; k++) {
    cellVector(mesh, model.state[2].subarray(k * E, (k + 1) * E), vector);
    for (let i = 0; i < C; i++) {
      const z = (geopotential[k * C + i] - model.surfaceGeopotential[i]) / g, speed = Math.hypot(vector[3 * i], vector[3 * i + 1], vector[3 * i + 2]);
      const expected = geography.land[i] ? handCoefficient(sigma, z) * speed : 0;
      close(bl.formRate[k * C + i], expected, 1e-9, `rate at cell ${i} layer ${k} (${z.toFixed(0)} m)`);
      checked++;
    }
  }
  for (let k = 0; k < bl.kTop; k++) for (let i = 0; i < C; i++) assert.equal(bl.formRate[k * C + i], 0, 'nothing above the solve');
  bl.mixing.fill(0); bl.surfaceDrag.fill(0);
  const u = model.state[2], before = Float64Array.from(u), heat = new Float64Array(K * E);
  bl.mixEdges(model.state[0], u, 0, E, dt, heat);
  for (let e = 0; e < E; e++) {
    const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1], columnMass = 0.5 * (model.state[0][a] + model.state[0][b]);
    let stress = 0, energy = 0, warmed = 0;
    for (let k = bl.kTop; k < K; k++) {
      const rate = 0.5 * (bl.formRate[k * C + a] + bl.formRate[k * C + b]), m = columnMass * dSigma[k] / g, idx = k * E + e;
      close(u[idx], before[idx] / (1 + dt * rate), 1e-12, `edge ${e} layer ${k}: u/(1 + Δt C_tofd |U|)`);
      stress += m * rate * u[idx];
      energy += m * (before[idx] ** 2 - u[idx] ** 2); warmed += m * heat[idx];
    }
    close(bl.formStress[e], stress, 1e-12, `edge ${e} stress`);
    if (Math.abs(before[(K - 1) * E + e]) > 1) close(warmed, energy, 1e-12, `edge ${e} energy to heat`);
  }
  assert.ok(checked > 0 && geography.land.some((l) => l), 'land cells checked');
});

test('on a real state the form drag conserves each edge column’s momentum against its stresses, heats by what it removes and stays off the sea', () => {
  const model = prepare(createModel(new Grid(16), { topography, levels: sigmaInterfaces('bl34') })), dt = 1350;
  assert.ok(model.boundaryLayer.formDrag, 'the bundled N=16 fields switch the form drag on');
  for (let n = 0; n < 3; n++) model.step(dt);
  const { mesh, core, boundaryLayer: bl, state, geography } = model, { K, C, E, dSigma, g } = core.diagnostics;
  bl.diagnose(state);
  const u = Float64Array.from(state[2]), before = Float64Array.from(u), heat = new Float64Array(K * E);
  bl.mixEdges(state[0], u, 0, E, dt, heat);
  let energy = 0, warmed = 0, formed = 0;
  for (let e = 0; e < E; e++) {
    const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1], columnMass = 0.5 * (state[0][a] + state[0][b]);
    let momentum = 0;
    for (let k = 0; k < K; k++) {
      const m = columnMass * dSigma[k] / g, idx = k * E + e;
      momentum += m * (before[idx] - u[idx]);
      energy += mesh.dcEdge[e] * mesh.dvEdge[e] * m * (before[idx] ** 2 - u[idx] ** 2);
      warmed += mesh.dcEdge[e] * mesh.dvEdge[e] * m * heat[idx];
    }
    const taken = (bl.surfaceStress[e] + bl.formStress[e]) * dt;
    assert.ok(Math.abs(momentum - taken) <= 1e-11 * Math.max(1, Math.abs(momentum)), `edge ${e}: ${momentum} against ${taken}`);
    if (!geography.land[a] && !geography.land[b]) assert.equal(bl.formStress[e], 0, `sea edge ${e}`);
    formed += Math.abs(bl.formStress[e]);
  }
  for (let i = 0; i < C; i++) if (!geography.land[i]) for (let k = 0; k < K; k++) assert.equal(bl.formRate[k * C + i], 0, `sea cell ${i}`);
  assert.ok(formed > 0, 'the land takes a form stress');
  close(warmed, energy, 1e-12, 'the kinetic energy removed becomes dissipation heat');
});

test('the form drag is off without fields, without σ_flt or when asked', () => {
  assert.equal(createModel(new Grid(8), { topography }).boundaryLayer.formDrag, null, 'N=8 has no file');
  assert.equal(createModel(new Grid(16), { topography, orography: { formDrag: false } }).boundaryLayer.formDrag, null);
  assert.equal(createModel(new Grid(16), { topography, orography: false }).boundaryLayer.formDrag, null);
  const fields = meshSubgrid(buildMesh(new Grid(16)));
  assert.equal(createModel(new Grid(16), { topography, subgrid: { ...fields, filtered: undefined } }).boundaryLayer.formDrag, null);
});

test('both engines lay the same form drag and stress', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const levels = sigmaInterfaces('bl34'), dt = 1350;
  const cpu = prepare(createModel(new Grid(16), { topography, levels }));
  const gpu = prepare(await createGpuModel(new Grid(16), { topography, levels }));
  cpu.step(dt); await gpu.step(dt); await gpu.settle();
  const PH = gpu.gpu.layout.PH, S = gpu.gpu.layout.S, C = cpu.mesh.nCells, E = cpu.mesh.nEdges, K = cpu.core.K, kTop = cpu.boundaryLayer.kTop;
  const [rate, stress, surface] = await readRanges(gpu.gpu.device, gpu.gpu.buffers.PH, [['TOFD', (K - kTop) * C], ['FSTRESS', E], ['STRESS', E]].map(([name, length]) => ({ offset: PH[name], length })));
  const [u] = await readRanges(gpu.gpu.device, gpu.gpu.buffers.S, [{ offset: S.U, length: K * E }]);
  const apart = (a, b) => { let worst = 0, scale = 0; for (let n = 0; n < b.length; n++) { worst = Math.max(worst, Math.abs(a[n] - b[n])); scale = Math.max(scale, Math.abs(b[n])); } return worst / scale; };
  const report = { rate: apart(rate, cpu.boundaryLayer.formRate.subarray(kTop * C)), formStress: apart(stress, cpu.boundaryLayer.formStress), surfaceStress: apart(surface, cpu.boundaryLayer.surfaceStress), u: apart(u, cpu.state[2]) };
  console.log(`engines apart, largest difference over largest value: ${Object.entries(report).map(([k, v]) => `${k} ${v.toExponential(1)}`).join(', ')}`);
  for (const [k, v] of Object.entries(report)) assert.ok(v < 2e-3, `${k} ${v}`);
});
