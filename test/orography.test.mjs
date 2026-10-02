import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography, topographyFromInt16, subgridOrography, createGeography, surfaceGeopotential } from '../js/geography.module.js';
import { orographicColumn } from '../js/physics/orography.module.js';
import { sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};
const { readRanges } = gpuAvailable ? await import('../js/gpu/device.module.js') : {};

const close = (actual, expected, tolerance, what) => assert.ok(Math.abs(actual - expected) <= tolerance * Math.abs(expected), `${what}: ${actual} against ${expected}`);
const deg = Math.PI / 180, R = 6371220;

function ridgeFields(elevationAt) {
  const mesh = buildMesh(new Grid(8));
  const topography = syntheticTopography(720, 1440, elevationAt);
  const geography = createGeography(mesh, topography, { landBridges: {}, seaStraits: {} });
  const phis = surfaceGeopotential(mesh, geography);
  return { mesh, fields: subgridOrography(mesh, topography, Float64Array.from(phis, (p) => p / 9.80616)) };
}

test('the subgrid fields of analytic ridges: standard deviation, slope, anisotropy and orientation', () => {
  const amplitude = 400, wavelength = 4 * deg, k = 2 * Math.PI / wavelength, sampling = Math.sin(k * 0.25 * deg) / (k * 0.25 * deg);
  const cases = [
    { name: 'east–west crests', h: (lat) => 1000 + amplitude * Math.cos(k * lat), band: 50, mu: amplitude / Math.SQRT2, slope: () => amplitude * k / R * sampling / Math.SQRT2, gamma: [0, 0.1], theta: Math.PI / 2 },
    { name: 'north–south crests', h: (lat, lon) => 1000 + amplitude * Math.cos(k * lon), band: 12, mu: amplitude / Math.SQRT2, slope: (lat) => amplitude * k / (R * Math.cos(lat)) * sampling / Math.SQRT2, gamma: [0, 0.1], theta: 0 },
    { name: 'oblique crests', h: (lat, lon) => 1000 + amplitude * Math.cos(k * (lat + lon)), band: 12, mu: amplitude / Math.SQRT2, slope: () => amplitude * k * Math.SQRT2 / R * sampling / Math.SQRT2, gamma: [0, 0.15], theta: Math.PI / 4 },
    { name: 'egg crate', h: (lat, lon) => 1000 + amplitude * Math.cos(k * lat) * Math.cos(k * lon), band: 12, mu: amplitude / 2, slope: () => amplitude * k / R * sampling / 2, gamma: [0.8, 1], theta: null },
  ];
  for (const c of cases) {
    const { mesh, fields } = ridgeFields(c.h);
    let n = 0;
    for (let i = 0; i < mesh.nCells; i++) {
      const lat = mesh.latCell[i];
      if (Math.abs(lat) > c.band * deg) continue;
      n++;
      assert.ok(fields.count[i] >= 16, `${c.name}: cell ${i} holds ${fields.count[i]} points`);
      close(fields.deviation[i], c.mu, 0.1, `${c.name} μ at ${(lat / deg).toFixed(1)}°`);
      close(fields.slope[i], c.slope(lat), 0.08, `${c.name} σ at ${(lat / deg).toFixed(1)}°`);
      assert.ok(fields.anisotropy[i] >= c.gamma[0] && fields.anisotropy[i] <= c.gamma[1], `${c.name} γ ${fields.anisotropy[i]}`);
      if (c.theta !== null) assert.ok(Math.abs(Math.sin(fields.orientation[i] - c.theta)) < 0.1, `${c.name} θ ${fields.orientation[i]} against ${c.theta}`);
    }
    assert.ok(n > 20, `${c.name}: ${n} cells tested`);
  }
  const flat = ridgeFields(() => 700).fields;
  for (let i = 0; i < flat.deviation.length; i++) assert.ok(flat.deviation[i] < 1e-6 && flat.slope[i] < 1e-9, 'a plateau has no subgrid orography');
});

function uniformColumn(east = 10, north = 0) {
  const levels = sigmaInterfaces('bl34'), K = levels.length - 1, g = 9.80616, N = 0.01;
  const column = { z: new Float64Array(K), p: new Float64Array(K), rho: new Float64Array(K), theta: new Float64Array(K), east: new Float64Array(K).fill(east), north: new Float64Array(K).fill(north), pTop: new Float64Array(K), pBottom: new Float64Array(K), g };
  for (let k = 0; k < K; k++) {
    const sigma = 0.5 * (levels[k] + levels[k + 1]);
    column.pTop[k] = 1e5 * levels[k]; column.pBottom[k] = 1e5 * levels[k + 1]; column.p[k] = 1e5 * sigma;
    column.z[k] = -8000 * Math.log(sigma);
    column.rho[k] = 1.2 * Math.min(1, Math.exp(-(column.z[k] - 1000) / 8000));
    column.theta[k] = 300 * Math.exp(N * N * column.z[k] / g);
  }
  return column;
}

test('a uniform flow normal to a ridge: the blocking height, the blocked drag and the wave stress by hand', () => {
  const column = uniformColumn(), K = column.z.length, g = column.g;
  const sub = { deviation: 200, anisotropy: 0, orientation: 0, slope: 0.02 };
  const out = orographicColumn(sub, column);
  close(out.blocking, 600 - 0.5 / (0.01 / 10), 1e-6, 'Z_b = 3μ − H_n,crit U/N');
  const launch = 1.2 * (600 - out.blocking) ** 2 / 9 * (0.02 / 200) * 1 * 10 * 1 * 0.01;
  close(out.launch, launch, 1e-6, 'τ = ρ H_eff²/9 (σ/μ) G |U| N');
  close(launch, 1 / 3, 1e-6, 'τ by hand');
  for (let k = 0; k < K; k++) {
    const z = column.z[k];
    const expected = z < out.blocking ? 1 * 2 * 0.02 / (2 * 200) * Math.sqrt((out.blocking - z) / (z + 200)) * 1 * 10 / 2 : 0;
    if (expected === 0) assert.equal(out.beta[k], 0, `no blocking at ${z.toFixed(0)} m`);
    else close(out.beta[k], expected, 1e-9, `blocking rate at ${z.toFixed(0)} m`);
  }
  const critical = (Math.SQRT2 - 1) / (2 * 0.25), launchAlpha = 0.01 * (600 - out.blocking) / 10;
  let tau = 0;
  for (let k = 0; k < K; k++) {
    const expectedAbove = k === 0 ? 0 : launch * Math.min(1, 0.5 * (column.rho[k] + column.rho[k - 1]) / 1.2 * (critical / launchAlpha) ** 2);
    tau = expectedAbove - out.wave[k] * (column.pBottom[k] - column.pTop[k]) / g;
    const expectedBelow = k === K - 1 ? launch : launch * Math.min(1, 0.5 * (column.rho[k] + column.rho[k + 1]) / 1.2 * (critical / launchAlpha) ** 2);
    assert.ok(Math.abs(tau - expectedBelow) < 2e-3 * launch, `stress at the base of layer ${k} (${column.z[k].toFixed(0)} m): ${tau} against ${expectedBelow}`);
  }
  close(out.deposited, out.launch, 1e-12, 'the column takes the whole launched stress');
  assert.deepEqual(out.direction, [1, 0]);
});

test('an oblique flow over an elongated ridge: the stress turns towards the cross-ridge axis', () => {
  const angle = 30 * deg, column = uniformColumn(10 * Math.cos(angle), 10 * Math.sin(angle));
  const gamma = 0.5, out = orographicColumn({ deviation: 200, anisotropy: gamma, orientation: 0, slope: 0.02 }, column);
  const B = 1 - 0.18 * gamma - 0.04 * gamma ** 2, C = 0.48 * gamma + 0.3 * gamma ** 2, psi = -angle;
  const D1 = B * Math.cos(psi) ** 2 + C * Math.sin(psi) ** 2, D2 = (B - C) * Math.sin(psi) * Math.cos(psi), D = Math.hypot(D1, D2);
  const along = [Math.cos(angle), Math.sin(angle)], cross = [-Math.sin(angle), Math.cos(angle)];
  close(out.direction[0], (D1 * along[0] + D2 * cross[0]) / D, 1e-12, 'stress east');
  close(out.direction[1], (D1 * along[1] + D2 * cross[1]) / D, 1e-12, 'stress north');
  assert.ok(Math.atan2(out.direction[1], out.direction[0]) < angle, 'the stress lies between the wind and the cross-ridge axis');
  close(out.launch, 1.2 * (600 - out.blocking) ** 2 / 9 * (0.02 / 200) * 10 * D * 0.01, 1e-6, 'τ with √(D1² + D2²)');
  const z = column.z[column.z.length - 1], r = (Math.cos(psi) ** 2 + gamma * Math.sin(psi) ** 2) / (gamma * Math.cos(psi) ** 2 + Math.sin(psi) ** 2);
  close(out.beta[column.z.length - 1], (2 - 1 / r) * 0.02 / 400 * Math.sqrt((out.blocking - z) / (z + 200)) * (B * Math.cos(psi) ** 2 + C * Math.sin(psi) ** 2) * 10 / 2, 1e-9, 'blocking rate with r and B cos²ψ + C sin²ψ');
});

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

test('the drag conserves momentum, hands it to the ground as a stress and returns its energy as heat', () => {
  const levels = sigmaInterfaces('bl34'), model = prepare(createModel(new Grid(8), { topography, levels })), dt = 1800;
  for (let n = 0; n < 3; n++) model.step(dt);
  const { mesh, orography, state, core } = model, { K, C, E, dSigma, g } = core.diagnostics;
  let launched = 0, taken = 0;
  for (let i = 0; i < C; i++) {
    let column = 0;
    for (let k = 0; k < K; k++) column -= state[0][i] * dSigma[k] / g * orography.wave[k * C + i];
    assert.ok(column <= orography.launch[i] * (1 + 1e-9) + 1e-12, `cell ${i} takes no more than it launches`);
    launched += orography.launch[i]; taken += column;
  }
  assert.ok(launched > 0 && taken > 0.95 * launched, `the columns take ${taken} of ${launched}`);
  const u = Float64Array.from(state[2]), before = Float64Array.from(u), lost = new Float64Array(K * E);
  orography.apply(state[0], u, 0, E, dt, lost);
  let energy = 0, heat = 0, blocked = 0;
  for (let e = 0; e < E; e++) {
    const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1], columnMass = 0.5 * (state[0][a] + state[0][b]);
    let momentum = 0;
    for (let k = 0; k < K; k++) {
      const m = columnMass * dSigma[k] / g, idx = k * E + e;
      momentum += m * (before[idx] - u[idx]);
      energy += mesh.dcEdge[e] * mesh.dvEdge[e] * m * (before[idx] ** 2 - u[idx] ** 2);
      heat += mesh.dcEdge[e] * mesh.dvEdge[e] * m * lost[idx];
      if (orography.beta[k * C + a] > 0) blocked++;
    }
    assert.ok(Math.abs(momentum - orography.stress[e] * dt) <= 1e-12 * Math.max(1, Math.abs(momentum)), `edge ${e}: the column loses what the ground receives`);
  }
  assert.ok(blocked > 0, 'some flow is blocked');
  assert.ok(energy > 0, 'the drag removes kinetic energy');
  close(heat, energy, 1e-12, 'the energy it removes becomes dissipation heat');
});

test('without the scheme the model has no orographic drag', () => {
  const model = createModel(new Grid(4), { topography, orography: false });
  assert.equal(model.orography, null);
});

test('both engines lay the same blocking, wave drag and stress', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const levels = sigmaInterfaces('bl34'), dt = 1350;
  const cpu = prepare(createModel(new Grid(16), { topography, levels }));
  const gpu = prepare(await createGpuModel(new Grid(16), { topography, levels }));
  cpu.step(dt); await gpu.step(dt); await gpu.settle();
  const PH = gpu.gpu.layout.PH, C = cpu.mesh.nCells, E = cpu.mesh.nEdges, K = cpu.core.K;
  const [beta, wave, stress, launch, drag] = await readRanges(gpu.gpu.device, gpu.gpu.buffers.PH, [['OBETA', K * C], ['OWAVE', K * C], ['OSTRESS', E], ['OLAUNCH', C], ['DRAG', C]].map(([name, length]) => ({ offset: PH[name], length })));
  const apart = (a, b) => { let worst = 0, scale = 0; for (let n = 0; n < b.length; n++) { worst = Math.max(worst, Math.abs(a[n] - b[n])); scale = Math.max(scale, Math.abs(b[n])); } return worst / scale; };
  const o = cpu.orography;
  const report = { beta: apart(beta, o.beta), wave: apart(wave, o.wave), stress: apart(stress, o.stress), launch: apart(launch, o.launch), drag: apart(drag, cpu.exchange.drag) };
  console.log(`engines apart, largest difference over largest value: ${Object.entries(report).map(([k, v]) => `${k} ${v.toExponential(1)}`).join(', ')}`);
  for (const [k, v] of Object.entries(report)) assert.ok(v < 1e-3, `${k} ${v}`);
  const off = await createGpuModel(new Grid(4), { topography, orography: false });
  assert.equal(off.gpu.kernels.orography, undefined);
});
