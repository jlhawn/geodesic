import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createParallelModel } from '../js/parallel.module.js';
import { createBoundaryLayer } from '../js/physics/boundaryLayer.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { saturationHumidity } from '../js/physics/moist.module.js';

const model = createModel(new Grid(3));
const { core, mesh, boundaryLayer } = model;
const { K, C, E, dSigma, g, geopotential } = core.diagnostics;

function column(surfaceTheta, lapseAloft, windAloft, stableAbove = Infinity) {
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  const [pi, theta, u, , q] = model.state;
  core.diagnose(pi, theta, q, model.state[5]);
  for (let i = 0; i < C; i++) {
    const zb = geopotential[(K - 1) * C + i] / g;
    for (let k = 0; k < K; k++) {
      const z = geopotential[k * C + i] / g - zb;
      theta[k * C + i] = surfaceTheta + lapseAloft * Math.min(z, stableAbove) + 5e-3 * Math.max(0, z - stableAbove);
      q[k * C + i] = 1e-2 * Math.exp(-z / 2000);
    }
  }
  for (let k = 0; k < K; k++) for (let e = 0; e < E; e++) u[k * E + e] = k === K - 1 ? 0 : windAloft * mesh.nEdge[3 * e + 1];
  core.diagnose(pi, theta, q, model.state[5]);
  return model.state;
}

const columnTotal = (field, pi, i) => { let s = 0; for (let k = 0; k < K; k++) s += pi[i] * dSigma[k] / g * field[k * C + i]; return s; };

test('an unstable, sheared column gets a kilometre-deep boundary layer that mixes θ and q toward uniform and conserves both', () => {
  const state = column(300, -1e-3, 8, 1000);
  const [pi, theta, , , q, qc] = state;
  boundaryLayer.diagnose(state);
  const i = 0, zb = geopotential[(K - 1) * C + i] / g;
  const depth = boundaryLayer.depth[i] - zb;
  assert.ok(depth > 500 && depth < 4000, `boundary layer ${depth} m deep`);
  let coefficients = 0;
  for (let k = boundaryLayer.kTop; k < K - 1; k++) if (boundaryLayer.mixing[k * C + i] > 0) coefficients++;
  assert.ok(coefficients >= 2, `${coefficients} interfaces mix`);
  const thetaBefore = columnTotal(theta, pi, i), qBefore = columnTotal(q, pi, i);
  const spreadBefore = Math.abs(theta[(K - 1) * C + i] - theta[(K - 3) * C + i]);
  for (let n = 0; n < 24; n++) boundaryLayer.mixColumn(i, pi, theta, q, qc, 900);
  const spreadAfter = Math.abs(theta[(K - 1) * C + i] - theta[(K - 3) * C + i]);
  assert.ok(spreadAfter < 0.5 * spreadBefore, `θ spread ${spreadBefore} → ${spreadAfter}`);
  assert.ok(Math.abs(columnTotal(theta, pi, i) - thetaBefore) < 1e-10 * thetaBefore);
  assert.ok(Math.abs(columnTotal(q, pi, i) - qBefore) < 1e-10 * qBefore);
  console.log(`unstable column: boundary layer ${depth.toFixed(0)} m, ${coefficients} mixing interfaces, θ spread over the lowest three layers ${spreadBefore.toFixed(2)} → ${spreadAfter.toFixed(2)} K after 6 h`);
});

test('a strongly stable column stays unmixed', () => {
  const state = column(280, 2e-2, 2);
  const [pi, theta, , , q, qc] = state;
  boundaryLayer.diagnose(state);
  const i = 0, before = Float64Array.from(theta);
  let coefficients = 0;
  for (let k = boundaryLayer.kTop; k < K - 1; k++) if (boundaryLayer.mixing[k * C + i] > 0) coefficients++;
  assert.equal(coefficients, 0, 'no interface mixes');
  boundaryLayer.mixColumn(i, pi, theta, q, qc, 900);
  for (let k = 0; k < K; k++) assert.equal(theta[k * C + i], before[k * C + i]);
});

test('momentum mixing brings wind down to the surface layer and conserves each edge column\'s momentum', () => {
  const state = column(300, -1e-3, 8, 1000);
  const [pi, , u] = state;
  boundaryLayer.diagnose(state);
  const e = mesh.edgesOnCell[mesh.maxEdges * 0];
  const total = () => { let s = 0; const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1]; for (let k = 0; k < K; k++) s += 0.5 * (pi[a] + pi[b]) * dSigma[k] / g * u[k * E + e]; return s; };
  const before = total(), surfaceBefore = Math.abs(u[(K - 1) * E + e]);
  for (let n = 0; n < 24; n++) boundaryLayer.mixEdges(pi, u, e, e + 1, 900);
  assert.ok(Math.abs(u[(K - 1) * E + e]) > surfaceBefore + 0.5, `surface wind ${surfaceBefore} → ${u[(K - 1) * E + e]}`);
  assert.ok(Math.abs(total() - before) < 1e-10 * Math.abs(before) + 1e-9);
});

test('serial and parallel engines stay bit-identical with the boundary layer', async () => {
  const serial = createModel(new Grid(6));
  const parallel = await createParallelModel(new Grid(6), {}, 3);
  try {
    const init = initializeState(serial, {});
    for (let a = 0; a < init.length; a++) { serial.state[a].set(init[a]); parallel.state[a].set(init[a]); }
    serial.ocean.initialize(serial.state[3], serial.state[6]);
    parallel.ocean.initialize(parallel.state[3], parallel.state[6]);
    for (let n = 0; n < 6; n++) { serial.step(900); parallel.step(900); }
    for (let a = 0; a < serial.state.length; a++) for (let x = 0; x < serial.state[a].length; x++) assert.equal(parallel.state[a][x], serial.state[a][x], `state ${a}[${x}]`);
  } finally { await parallel.close(); }
});

test('a warm or moist sea surface deepens the momentum mixing; a cool one keeps the neutral profile', () => {
  const { exnerLayer, sigmaMid, kappa } = core.diagnostics;
  const neutralLayer = createBoundaryLayer(mesh, core, { stability: false });
  const sum = (layer) => { let s = 0; for (let k = layer.kTop; k < K - 1; k++) s += layer.mixing[k * C]; return s; };
  const compare = (offset, humidity) => {
    const state = column(300, -1e-3, 8, 1000);
    const base = (K - 1) * C;
    for (let i = 0; i < C; i++) {
      state[3][i] = state[1][base + i] * exnerLayer[base + i] * Math.pow(sigmaMid[K - 1], -kappa) + offset;
      state[4][base + i] = humidity * saturationHumidity(state[3][i], state[0][i]);
    }
    boundaryLayer.diagnose(state); neutralLayer.diagnose(state);
    return sum(boundaryLayer) / sum(neutralLayer);
  };
  const warm = compare(3, 1), cool = compare(-3, 1), moist = compare(0, 0.7);
  assert.ok(warm > 1.3 && warm < 2.4, `warm surface: ${warm} times the neutral mixing`);
  assert.ok(moist > 1.2 && moist < 2.4, `moist surface: ${moist} times the neutral mixing`);
  assert.ok(Math.abs(cool - 1) < 1e-12, `cool surface: ${cool} times the neutral mixing`);
});

test('where a mixed-layer deck runs the K-profile spans its inversion when that lies above the Richardson depth, depth itself stays the Richardson depth, and a zero deckTop changes nothing', () => {
  const state = column(290, 3e-3, 4);
  const deckTop = new Float64Array(C), plain = createBoundaryLayer(mesh, core, {}), decked = createBoundaryLayer(mesh, core, { deckTop });
  plain.diagnose(state); decked.diagnose(state);
  assert.deepEqual(decked.mixing, plain.mixing);
  assert.deepEqual(decked.depth, plain.depth);
  for (let i = 0; i < C; i++) deckTop[i] = i % 3 === 0 ? 0 : i % 3 === 1 ? 0.5 * plain.depth[i] : plain.depth[i] + 800;
  decked.diagnose(state);
  assert.deepEqual(decked.depth, plain.depth, 'the Richardson depth is the deck\'s floor and stays as found');
  let deepened = 0, added = 0;
  for (let i = 0; i < C; i++) {
    const zb = geopotential[(K - 1) * C + i] / g;
    let before = 0, after = 0;
    for (let k = decked.kTop; k < K - 1; k++) {
      const idx = k * C + i, z = 0.5 * (geopotential[idx] + geopotential[idx + C]) / g - zb;
      if (i % 3 !== 2) { assert.equal(decked.mixing[idx], plain.mixing[idx], `cell ${i} interface ${k}`); continue; }
      assert.ok(decked.mixing[idx] >= plain.mixing[idx], `cell ${i} interface ${k}: ${decked.mixing[idx]} against ${plain.mixing[idx]}`);
      assert.equal(decked.mixing[idx] > 0, z < deckTop[i] - zb, `cell ${i} interface ${k} at ${z} m under a deck at ${deckTop[i] - zb} m`);
      if (plain.mixing[idx] > 0) before++;
      if (decked.mixing[idx] > 0) after++;
    }
    if (after > before) { deepened++; added += after - before; }
  }
  console.log(`a deck 800 m above the Richardson depth adds ${added} mixing interfaces over ${deepened} of ${Math.ceil(C / 3)} decked columns`);
  assert.ok(deepened > 0.5 * Math.floor(C / 3), `${deepened} columns mix deeper`);
});

test('serial and parallel engines stay bit-identical with the mixed-layer deck carrying its height and gate', async () => {
  const radiation = { stratusSubsidence: 0, minimumInversion: 0 };
  const serial = createModel(new Grid(6), { radiation });
  const parallel = await createParallelModel(new Grid(6), { radiation }, 3);
  try {
    const init = initializeState(serial, {});
    for (let a = 0; a < init.length; a++) { serial.state[a].set(init[a]); parallel.state[a].set(init[a]); }
    serial.ocean.initialize(serial.state[3], serial.state[6]);
    parallel.ocean.initialize(parallel.state[3], parallel.state[6]);
    for (let n = 0; n < 6; n++) { serial.step(900); parallel.step(900); }
    let decked = 0;
    for (let i = 0; i < serial.mesh.nCells; i++) if (serial.radiation.mlmTop[i] > 0) decked++;
    assert.ok(decked > 0, 'the deck runs and hands its height to the boundary layer');
    for (let a = 0; a < serial.state.length; a++) for (let x = 0; x < serial.state[a].length; x++) assert.equal(parallel.state[a][x], serial.state[a][x], `state ${a}[${x}]`);
    for (const name of ['mlmHeight', 'mlmGate', 'mlmTop', 'mlmCover', 'mlmWater']) assert.deepEqual(parallel.radiation[name], serial.radiation[name], name);
  } finally { await parallel.close(); }
});
