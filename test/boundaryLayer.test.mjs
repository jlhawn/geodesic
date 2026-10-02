import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createParallelModel } from '../js/parallel.module.js';
import { createBoundaryLayer } from '../js/physics/boundaryLayer.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { saturationHumidity } from '../js/physics/moist.module.js';
import { sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { SEA_DRAG } from '../js/physics/surface.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }

const model = createModel(new Grid(3), { boundaryLayer: { turbulence: 'dry' } });
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

test('momentum mixing brings wind down to the surface layer and changes each edge column\'s momentum by the surface stress alone', () => {
  const state = column(300, -1e-3, 8, 1000);
  const [pi, , u] = state;
  boundaryLayer.diagnose(state);
  const e = mesh.edgesOnCell[mesh.maxEdges * 0];
  const total = () => { let s = 0; const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1]; for (let k = 0; k < K; k++) s += 0.5 * (pi[a] + pi[b]) * dSigma[k] / g * u[k * E + e]; return s; };
  const before = total(), surfaceBefore = Math.abs(u[(K - 1) * E + e]);
  let lost = 0;
  for (let n = 0; n < 24; n++) { boundaryLayer.mixEdges(pi, u, e, e + 1, 900); lost += 900 * boundaryLayer.surfaceStress[e]; }
  assert.ok(boundaryLayer.implicitDrag && lost !== 0);
  assert.ok(Math.abs(u[(K - 1) * E + e]) > surfaceBefore + 0.5, `surface wind ${surfaceBefore} → ${u[(K - 1) * E + e]}`);
  assert.ok(Math.abs(total() - before + lost) < 1e-10 * Math.abs(before) + 1e-9, `momentum ${before} → ${total()}, the surface took ${lost}`);
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
  const neutralLayer = createBoundaryLayer(mesh, core, { turbulence: 'dry', stability: false }), bulkLayer = createBoundaryLayer(mesh, core, { turbulence: 'dry' });
  const sum = (layer) => { let s = 0; for (let k = layer.kTop; k < K - 1; k++) s += layer.mixing[k * C]; return s; };
  const compare = (offset, humidity) => {
    const state = column(300, -1e-3, 8, 1000);
    const base = (K - 1) * C;
    for (let i = 0; i < C; i++) {
      state[3][i] = state[1][base + i] * exnerLayer[base + i] * Math.pow(sigmaMid[K - 1], -kappa) + offset;
      state[4][base + i] = humidity * saturationHumidity(state[3][i], state[0][i]);
    }
    bulkLayer.diagnose(state); neutralLayer.diagnose(state);
    return sum(bulkLayer) / sum(neutralLayer);
  };
  const warm = compare(3, 1), cool = compare(-3, 1), moist = compare(0, 0.7);
  assert.ok(warm > 1.3 && warm < 2.4, `warm surface: ${warm} times the neutral mixing`);
  assert.ok(moist > 1.2 && moist < 2.4, `moist surface: ${moist} times the neutral mixing`);
  assert.ok(Math.abs(cool - 1) < 1e-12, `cool surface: ${cool} times the neutral mixing`);
});

test('where a mixed-layer deck runs the K-profile spans its inversion when that lies above the Richardson depth, depth itself stays the Richardson depth, and a zero deckTop changes nothing', () => {
  const state = column(290, 3e-3, 4);
  const deckTop = new Float64Array(C), still = { turbulence: 'dry', entrainment: { efficiency: 0, shear: 0 } }, plain = createBoundaryLayer(mesh, core, still), decked = createBoundaryLayer(mesh, core, { turbulence: 'dry', deckTop, ...still });
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

function entrainingColumn(surfaceWarmth = 2) {
  const state = column(300, -1e-3, 8, 1000);
  const [pi, theta, , surfaceT, q, qc] = state;
  const { exnerLayer, sigmaMid, kappa } = core.diagnostics;
  const base = (K - 1) * C;
  for (let i = 0; i < C; i++) {
    for (let k = 0; k < K; k++) qc[k * C + i] = k === K - 3 ? 1e-4 : 0;
    surfaceT[i] = theta[base + i] * exnerLayer[base + i] * Math.pow(sigmaMid[K - 1], -kappa) + surfaceWarmth;
  }
  boundaryLayer.diagnose(state);
  for (let i = 0; i < C; i++) {
    const zb = geopotential[base + i] / g, h = boundaryLayer.depth[i] - zb;
    for (let k = 0; k < K - 1; k++) if (0.5 * (geopotential[k * C + i] + geopotential[(k + 1) * C + i]) / g - zb >= h) q[k * C + i] = 1e-3;
  }
  core.diagnose(pi, theta, q, qc);
  return state;
}

function expectedEntrainment(state, layer, i, { efficiency = 0.2, shear = 5, cap = 0.05, jumpFloor = 0.015, shearOnset = 5e-5 } = {}) {
  const [pi, theta, , surfaceT, q] = state;
  const { exnerLayer, sigmaMid, kappa } = core.diagnostics, { thetaV } = core.arrays;
  const base = (K - 1) * C + i, zb = geopotential[base] / g, h = layer.depth[i] - zb;
  const friction = Math.sqrt(SEA_DRAG) * 3;
  const buoyancy = g / theta[base] * SEA_DRAG * 3 * (surfaceT[i] * Math.pow(sigmaMid[K - 1], kappa) / exnerLayer[base] - theta[base] + 0.61 * theta[base] * (saturationHumidity(surfaceT[i], pi[i]) - q[base]));
  let above = -1;
  for (let k = layer.kTop; k < K - 1; k++) if (0.5 * (geopotential[k * C + i] + geopotential[(k + 1) * C + i]) / g - zb >= h) above = k;
  let weight = 0, sum = 0;
  for (let k = above + 1; k < K; k++) { weight += dSigma[k]; sum += dSigma[k] * thetaV[k * C + i]; }
  const jump = g * (thetaV[above * C + i] - sum / weight) / (sum / weight);
  const onset = shearOnset > 0 ? Math.min(1, buoyancy / shearOnset) : 1;
  return { above, buoyancy, jump, velocity: Math.min(cap, (efficiency * buoyancy + shear * onset * friction ** 3 / h) / Math.max(jump, jumpFloor)) };
}

test('a dry stable layer over a convective boundary layer is entrained at the closure\'s rate across the first interface above h, conserving the column\'s θ, water and momentum', () => {
  const state = entrainingColumn();
  const [pi, theta, u, , q, qc] = state;
  const i = 0;
  const layer = createBoundaryLayer(mesh, core, { turbulence: 'dry' }), still = createBoundaryLayer(mesh, core, { turbulence: 'dry', entrainment: { efficiency: 0, shear: 0 } });
  layer.diagnose(state); still.diagnose(state);
  const { above, buoyancy, jump, velocity } = expectedEntrainment(state, layer, i);
  assert.ok(buoyancy > 0 && jump > 0 && above >= layer.kTop, `B0 ${buoyancy}, Δb ${jump}, layer ${above}`);
  assert.ok(q[above * C + i] === 1e-3 && q[(above + 1) * C + i] > 2e-3, 'the layer above h is the dry one');
  assert.ok(Math.abs(layer.entrainment[i] - velocity) < 1e-12 * velocity, `w_e ${layer.entrainment[i]} against ${velocity}`);
  assert.equal(still.entrainment[i], 0);
  const { exnerLayer, sigmaMid, R } = core.diagnostics, idx = above * C + i;
  const rho = 0.5 * (pi[i] * sigmaMid[above] / (R * theta[idx] * exnerLayer[idx]) + pi[i] * sigmaMid[above + 1] / (R * theta[idx + C] * exnerLayer[idx + C]));
  assert.ok(Math.abs(layer.mixing[idx] - rho * velocity) < 1e-12 * rho * velocity, 'the interface above h carries ρ w_e');
  assert.equal(still.mixing[idx], 0, 'which the K-profile alone leaves shut');
  for (let k = layer.kTop; k < K - 1; k++) if (k !== above) assert.equal(layer.mixing[k * C + i], still.mixing[k * C + i], `interface ${k}`);

  const dt = 900, mass = pi[i] * dSigma[above] / g;
  const totals = () => [columnTotal(theta, pi, i), columnTotal(q, pi, i), columnTotal(qc, pi, i)];
  const before = totals(), qAbove = q[idx], thetaAbove = theta[idx];
  const plainTheta = Float64Array.from(theta), plainQ = Float64Array.from(q), plainQc = Float64Array.from(qc);
  layer.mixColumn(i, pi, theta, q, qc, dt);
  const after = totals();
  for (let n = 0; n < 3; n++) assert.ok(Math.abs(after[n] - before[n]) <= 4e-16 * Math.abs(before[n]) * K, `total ${n}: ${before[n]} → ${after[n]}`);
  assert.ok(Math.abs(mass * (q[idx] - qAbove) - dt * rho * velocity * (q[idx + C] - q[idx])) < 1e-12 * mass * Math.abs(q[idx] - qAbove), 'the layer above h gains ρ w_e dt of the boundary layer\'s water');
  assert.ok(Math.abs(mass * (theta[idx] - thetaAbove) - dt * rho * velocity * (theta[idx + C] - theta[idx])) < 1e-9 * mass * Math.abs(theta[idx] - thetaAbove), 'and loses heat at the same rate');
  still.mixColumn(i, pi, plainTheta, plainQ, plainQc, dt);
  assert.equal(plainQ[idx], qAbove, 'without entrainment the layer above h keeps its water');
  let bl = 0, blPlain = 0, blTheta = 0, blThetaPlain = 0;
  for (let k = above + 1; k < K; k++) { bl += dSigma[k] * q[k * C + i]; blPlain += dSigma[k] * plainQ[k * C + i]; blTheta += dSigma[k] * theta[k * C + i]; blThetaPlain += dSigma[k] * plainTheta[k * C + i]; }
  assert.ok(bl < blPlain && blTheta > blThetaPlain, 'the boundary layer dries and warms');

  const e = mesh.edgesOnCell[mesh.maxEdges * i], a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1];
  const momentum = () => { let s = 0; for (let k = 0; k < K; k++) s += 0.5 * (pi[a] + pi[b]) * dSigma[k] / g * u[k * E + e]; return s; };
  const uAbove = u[above * E + e], m0 = momentum(), lost = new Float64Array(K * E);
  layer.mixEdges(pi, u, e, e + 1, dt, lost);
  assert.ok(Math.abs(momentum() - m0) < 1e-13 * Math.abs(m0), `momentum ${m0} → ${momentum()}`);
  assert.ok(Math.abs(u[above * E + e]) < Math.abs(uAbove), 'the layer above h gives up momentum to the slower boundary layer');
  console.log(`entraining column: h ${(layer.depth[i] - geopotential[(K - 1) * C + i] / g).toFixed(0)} m, B0 ${buoyancy.toExponential(2)} m²/s³, Δθv ${(jump * 300 / g).toFixed(2)} K, w_e ${(1000 * velocity).toFixed(2)} mm/s; one 900 s step moves ${(1000 * (qAbove - plainQ[idx] + q[idx] - qAbove) * mass).toFixed(1)} g/m² of water above h; column θ, water and momentum change by ${((after[0] - before[0]) / before[0]).toExponential(1)}, ${((after[1] + after[2] - before[1] - before[2]) / (before[1] + before[2])).toExponential(1)}, ${((momentum() - m0) / m0).toExponential(1)}`);
});

test('the deck\'s opening and the stratiform share taper w_e: a half-open gate or a share of one half halves it, a closed gate or an EIS of 12 K stops it; a stable surface, zero coefficients, the cap and the jump floor bound it', () => {
  const deckGate = new Float64Array(C).fill(0.3), stratiform = new Float64Array(C);
  const layer = createBoundaryLayer(mesh, core, { turbulence: 'dry', deckGate, stratiform });
  const warm = entrainingColumn(2);
  layer.diagnose(warm);
  const full = expectedEntrainment(warm, layer, 0).velocity;
  assert.ok(full > 0 && Math.abs(layer.entrainment[0] - full) < 1e-12 * full);
  const at = (gate, share) => { deckGate[0] = gate; stratiform[0] = share; layer.diagnose(warm); return layer.entrainment[0]; };
  assert.equal(at(0.5, 0), full, 'a gate of one half leaves the top to the boundary layer');
  assert.ok(Math.abs(at(0.55, 0) - 0.5 * full) < 1e-12 * full, 'a half-open deck halves it');
  assert.equal(at(0.6, 0), 0, 'a closed deck owns the top');
  assert.equal(at(0.9, 0), 0);
  assert.ok(Math.abs(at(0.3, 0.5) - 0.5 * full) < 1e-12 * full, 'an EIS of 10 K halves it');
  assert.equal(at(0.3, 1), 0, 'an EIS of 12 K stops it');
  assert.ok(Math.abs(at(0.55, 0.5) - 0.25 * full) < 1e-12 * full, 'the two tapers multiply');
  const idx = expectedEntrainment(warm, layer, 0).above * C;
  const { exnerLayer, sigmaMid, R } = core.diagnostics, [pi, theta] = warm;
  const rho = 0.5 * (pi[0] * sigmaMid[idx / C] / (R * theta[idx] * exnerLayer[idx]) + pi[0] * sigmaMid[idx / C + 1] / (R * theta[idx + C] * exnerLayer[idx + C]));
  at(0.3, 0.75);
  assert.ok(Math.abs(layer.mixing[idx] - rho * 0.25 * full) < 1e-12 * rho * full, 'the interface carries ρ times the tapered w_e');
  deckGate[0] = 0.3; stratiform[0] = 0;
  const cool = entrainingColumn(-8);
  layer.diagnose(cool);
  for (let i = 0; i < C; i++) {
    assert.ok(expectedEntrainment(cool, layer, i).buoyancy < 0, `cell ${i} is stable at its surface`);
    assert.equal(layer.entrainment[i], 0, `stable surface, cell ${i}`);
  }
  const capped = createBoundaryLayer(mesh, core, { turbulence: 'dry', entrainment: { cap: 1e-4 } }), floored = createBoundaryLayer(mesh, core, { turbulence: 'dry', entrainment: { jumpFloor: 10 } });
  const again = entrainingColumn(2), off = createBoundaryLayer(mesh, core, { turbulence: 'dry', deckGate, entrainment: { efficiency: 0, shear: 0 } });
  capped.diagnose(again); floored.diagnose(again); off.diagnose(again);
  assert.ok(off.entrainment.every((x) => x === 0));
  const { buoyancy } = expectedEntrainment(again, floored, 1);
  assert.ok(buoyancy > 5e-5, 'the shear term is whole');
  assert.equal(capped.entrainment[1], 1e-4);
  assert.ok(Math.abs(floored.entrainment[1] - (0.2 * buoyancy + 5 * (Math.sqrt(SEA_DRAG) * 3) ** 3 / (floored.depth[1] - geopotential[(K - 1) * C + 1] / g)) / 10) < 1e-12);
});

test('the shear term comes in continuously with the surface buoyancy flux: w_e falls to zero as B0 falls to zero, linearly below the onset; without the onset it jumps', () => {
  const layer = createBoundaryLayer(mesh, core, { turbulence: 'dry' }), switched = createBoundaryLayer(mesh, core, { turbulence: 'dry', entrainment: { shearOnset: 0 } });
  const state = entrainingColumn(0), surfaceT = state[3], neutral = surfaceT[0];
  const flux = (warmth) => { surfaceT[0] = neutral + warmth; layer.diagnose(state); return layer.buoyancyFlux[0]; };
  let cold = -8, warm = 2;
  for (let n = 0; n < 60; n++) { const mid = 0.5 * (cold + warm); if (flux(mid) > 0) warm = mid; else cold = mid; }
  const samples = [];
  for (const target of [1e-4, 5e-5, 2.5e-5, 1e-5, 1e-6, 1e-8]) {
    let lo = cold, hi = 2;
    for (let n = 0; n < 60; n++) { const mid = 0.5 * (lo + hi); if (flux(mid) > target) hi = mid; else lo = mid; }
    flux(hi); switched.diagnose(state);
    const expected = expectedEntrainment(state, layer, 0);
    assert.ok(Math.abs(layer.entrainment[0] - expected.velocity) <= 1e-12 * expected.velocity, `B0 ${expected.buoyancy}: w_e ${layer.entrainment[0]} against ${expected.velocity}`);
    samples.push([expected.buoyancy, layer.entrainment[0], switched.entrainment[0]]);
  }
  const [, wOnset] = samples[1], [b6, w6, s6] = samples[4], [b8, w8, s8] = samples[5];
  assert.ok(w6 < 0.05 * wOnset && w8 < 1e-3 * wOnset, `w_e ${w6} at B0 ${b6} and ${w8} at ${b8} against ${wOnset} at the onset`);
  assert.ok(Math.abs(w6 / b6 - w8 / b8) < 1e-3 * (w6 / b6), 'linear in B0 below the onset');
  assert.ok(s8 > 100 * w8 && s6 > 2 * w6, `without the onset the shear term alone gives ${s8} at B0 ${b8}`);
  flux(cold);
  assert.ok(layer.buoyancyFlux[0] <= 0);
  assert.equal(layer.entrainment[0], 0);
  console.log(`shear onset: w_e ${samples.map(([b, w, x]) => `${(1000 * w).toFixed(3)} (switched ${(1000 * x).toFixed(3)}) mm/s at B0 ${b.toExponential(1)}`).join(', ')}`);
});

async function engines(entrainment) {
  const { createGpuCore } = await import('../js/gpu/core.gpu.js');
  const levels = sigmaInterfaces('bl34');
  const pair = createModel(new Grid(6), { ocean: false, levels, boundaryLayer: { entrainment, turbulence: 'dry' }, surface: { exchange: 'fixed' } });
  const { core: c, mesh: m, state, radiation, moist, boundaryLayer: layer } = pair;
  const { K: nK, C: nC, E: nE, exnerLayer, sigmaMid, kappa, geopotential: phi, g: grav } = c.diagnostics;
  const [pi, theta, u, surfaceT, q, qc] = state;
  const init = initializeState(pair, {});
  for (let a = 0; a < init.length; a++) state[a].set(init[a]);
  c.diagnose(pi, theta, q, qc);
  let seed = 2024;
  const random = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; };
  for (let i = 0; i < nC; i++) {
    const zb = phi[(nK - 1) * nC + i] / grav, surface = 286 + 14 * random(), top = 300 + 1500 * random(), below = -1e-3 + 2e-3 * random(), above = 3e-3 + 5e-3 * random(), jump = 4 * random();
    const wetBelow = 0.4 + 0.3 * random(), wetAbove = 0.1 + 0.3 * random();
    for (let k = 0; k < nK; k++) {
      const x = k * nC + i, z = phi[x] / grav - zb;
      theta[x] = z < top ? surface + below * z : surface + below * top + jump + above * (z - top);
      if (sigmaMid[k] < 0.2) theta[x] = Math.max(theta[x], init[1][x]);
      q[x] = (z < top ? wetBelow : wetAbove) * saturationHumidity(theta[x] * exnerLayer[x], pi[i] * sigmaMid[k]);
      qc[x] = z < 3000 && random() < 0.1 ? 2e-4 * random() : 0;
    }
    surfaceT[i] = theta[(nK - 1) * nC + i] * exnerLayer[(nK - 1) * nC + i] * Math.pow(sigmaMid[nK - 1], -kappa) - 3 + 6 * random();
    const gate = random();
    radiation.mlmGate[i] = gate < 0.2 ? 0.7 : gate < 0.3 ? 0.5 : gate < 0.45 ? 0.5 + 0.1 * random() : 0.3 * random();
    const share = random();
    radiation.stratiform[i] = share < 0.1 ? 1 : share < 0.35 ? random() : 0;
  }
  for (let k = 0; k < nK; k++) for (let e = 0; e < nE; e++) u[k * nE + e] = k === nK - 1 ? 2 * (random() - 0.5) : 12 * (random() - 0.5);
  for (const a of state) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  for (const a of [radiation.mlmGate, radiation.stratiform]) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  c.diagnose(pi, theta, q, qc);
  const gpu = await createGpuCore(m, { levels, physics: { entrainment, turbulence: 'dry', surfaceExchange: 'fixed' } });
  const { device, buffers, kernels } = gpu, dt = 900;
  gpu.upload(state);
  gpu.uploadPhysics({ mlmGate: radiation.mlmGate });
  device.queue.writeBuffer(buffers.PH, 4 * gpu.layout.PH.STRAT, Float32Array.from(radiation.stratiform));
  device.queue.writeBuffer(buffers.P, 0, Float32Array.from([dt, 0, 1, 0, 0, 0, 0, 0]));
  const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
  for (const [name, count] of [['pblDiagnose', nC], ['adjust', nC], ['mixMomentum', nE]]) {
    pass.setPipeline(kernels[name]);
    pass.setBindGroup(0, device.createBindGroup({ layout: kernels[name].getBindGroupLayout(0), entries: [buffers.MI, buffers.MF, buffers.LV, buffers.S, buffers.K1, buffers.D, buffers.P, buffers.PH].map((buffer, binding) => ({ binding, resource: { buffer } })) }));
    pass.dispatchWorkgroups(Math.ceil(count / 64));
  }
  pass.end();
  device.queue.submit([encoder.finish()]);
  const after = await gpu.download(), ph = await gpu.downloadPhysics();
  const before = state.map((a) => Float64Array.from(a));
  layer.diagnose(state);
  const cpu = { entrainment: Float64Array.from(layer.entrainment), mixing: Float64Array.from(layer.mixing) };
  pair.phases.adjust(0, nC, dt);
  pair.phases.mixMomentum(0, nE, dt);
  return { K: nK, C: nC, E: nE, kTop: layer.kTop, gate: radiation.mlmGate, share: radiation.stratiform, before, cpu, state, after, ph };
}

test('the dry boundary layer with entrainment matches between the engines on a random set of columns: w_e, the interface coefficients, and θ, q, qc and the wind after the step', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const on = await engines({}), off = await engines({ efficiency: 0, shear: 0 });
  const { K: nK, C: nC, E: nE, kTop } = on;
  let entraining = 0, gated = 0, stable = 0, worstW = 0, scaleW = 0, worstMix = 0, flips = 0, tapered = 0;
  for (let i = 0; i < nC; i++) {
    const cpu = on.cpu.entrainment[i], gpu = on.ph.ENTRAIN[i];
    if (cpu > 0) entraining++; else if (on.gate[i] >= 0.6 || on.share[i] >= 1) gated++; else stable++;
    if (cpu > 0 && (on.gate[i] > 0.5 || on.share[i] > 0)) tapered++;
    if ((cpu > 0) !== (gpu > 0)) { flips++; continue; }
    worstW = Math.max(worstW, Math.abs(cpu - gpu) / Math.max(cpu, 1e-4)); scaleW = Math.max(scaleW, cpu);
    let largest = 0;
    for (let k = kTop; k < nK - 1; k++) largest = Math.max(largest, on.cpu.mixing[k * nC + i]);
    for (let k = kTop; k < nK - 1; k++) { const x = k * nC + i; if (largest > 0) worstMix = Math.max(worstMix, Math.abs(on.cpu.mixing[x] - on.ph.MIX[x]) / largest); }
  }
  let theta = 0, q = 0, qc = 0, wind = 0, moved = 0, movedWind = 0;
  for (let x = 0; x < nK * nC; x++) {
    theta = Math.max(theta, Math.abs(on.state[1][x] - on.after[1][x])); q = Math.max(q, Math.abs(on.state[4][x] - on.after[4][x])); qc = Math.max(qc, Math.abs(on.state[5][x] - on.after[5][x]));
    moved = Math.max(moved, Math.abs(on.state[4][x] - off.state[4][x]));
  }
  for (let x = 0; x < nK * nE; x++) { wind = Math.max(wind, Math.abs(on.state[2][x] - on.after[2][x])); movedWind = Math.max(movedWind, Math.abs(on.state[2][x] - off.state[2][x])); }
  let offW = 0;
  for (let i = 0; i < nC; i++) offW = Math.max(offW, off.ph.ENTRAIN[i], off.cpu.entrainment[i]);
  console.log(`${nC} random columns: ${entraining} entrain (w_e up to ${(1000 * scaleW).toFixed(1)} mm/s, ${tapered} of them tapered by the deck's opening or the stratiform share), ${gated} under a closed deck or a share of one, ${stable} over a stable surface; entrainment on or off differs between the engines on ${flips}; w_e differs by ${worstW.toExponential(1)} relative, the interface coefficients by ${worstMix.toExponential(1)} of each column's largest; after the step θ by ${theta.toExponential(1)} K, q by ${q.toExponential(1)}, qc by ${qc.toExponential(1)}, the wind by ${wind.toExponential(1)} m/s, where entrainment moves q by up to ${moved.toExponential(1)} and the wind by ${movedWind.toFixed(2)} m/s`);
  assert.ok(entraining > nC / 4 && gated > nC / 10 && stable > nC / 20 && tapered > nC / 10, `${entraining} entraining, ${tapered} tapered, ${gated} gated, ${stable} stable`);
  assert.equal(offW, 0);
  assert.ok(flips <= nC / 200, `${flips} columns flip`);
  assert.ok(worstW < 1e-3 && worstMix < 1e-3, `w_e ${worstW}, coefficients ${worstMix}`);
  assert.ok(theta < 1e-3 && q < 1e-6 && qc < 1e-7 && wind < 1e-3, `θ ${theta}, q ${q}, qc ${qc}, wind ${wind}`);
  assert.ok(moved > 100 * q && movedWind > 100 * wind, `entrainment moves q by ${moved} and the wind by ${movedWind}`);
});

test('a surface parcel that tops out at the base of a stratocumulus whose descending parcel stops there is coupled to it in both engines', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { createGpuCore } = await import('../js/gpu/core.gpu.js');
  const levels = sigmaInterfaces('bl34');
  const pair = createModel(new Grid(6), { ocean: false, levels, surface: { exchange: 'fixed' } });
  const { core: c, mesh: m, state, radiation, boundaryLayer: layer } = pair;
  const { K: nK, C: nC, E: nE, exnerLayer, sigmaMid, kappa, geopotential: phi, g: grav } = c.diagnostics;
  const [pi, theta, u, surfaceT, q, qc] = state;
  const init = initializeState(pair, {});
  for (let a = 0; a < init.length; a++) state[a].set(init[a]);
  c.diagnose(pi, theta, q, qc);
  let seed = 7;
  const random = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; };
  for (let i = 0; i < nC; i++) {
    const zb = phi[(nK - 1) * nC + i] / grav, height = (k) => phi[k * nC + i] / grav - zb;
    let base = nK - 1;
    const wanted = 250 + 400 * random();
    while (base > 0 && height(base) < wanted) base--;
    const top = base - 1 - Math.floor(2 * random()), surface = 285 + 12 * random(), water = 2e-4 + 2e-4 * random();
    for (let k = 0; k < nK; k++) {
      const x = k * nC + i, p = pi[i] * sigmaMid[k];
      if (k > base) { theta[x] = surface; qc[x] = 0; }
      else if (k >= top) { theta[x] = surface + 3; qc[x] = water; }
      else { theta[x] = Math.max(surface + 10 + 4e-3 * (height(k) - height(top)), init[1][x]); qc[x] = 0; }
      q[x] = k >= top ? saturationHumidity(theta[x] * exnerLayer[x], p) : 0.3 * saturationHumidity(theta[x] * exnerLayer[x], p);
    }
    const lifted = (base + 2) * nC + i;
    for (let k = base + 1; k < nK; k++) q[k * nC + i] = 1.06 * saturationHumidity(theta[lifted] * exnerLayer[lifted], pi[i] * sigmaMid[base + 2]);
    surfaceT[i] = theta[(nK - 1) * nC + i] * exnerLayer[(nK - 1) * nC + i] * Math.pow(sigmaMid[nK - 1], -kappa) + 1 + 2 * random();
    for (let k = 0; k < nK; k++) radiation.longwave[k * nC + i] = k === top ? -60 : 0;
  }
  u.fill(0);
  for (const a of state) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  c.diagnose(pi, theta, q, qc);
  const gpu = await createGpuCore(m, { levels, physics: { surfaceExchange: 'fixed' } });
  const { device, buffers, kernels } = gpu;
  gpu.upload(state);
  gpu.uploadPhysics();
  device.queue.writeBuffer(buffers.PH, 4 * gpu.layout.PH.LWH, Float32Array.from(radiation.longwave.subarray(0, nK * nC)));
  device.queue.writeBuffer(buffers.P, 0, Float32Array.from([900, 0, 1, 0, 0, 0, 0, 0]));
  const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
  pass.setPipeline(kernels.pblDiagnose);
  pass.setBindGroup(0, device.createBindGroup({ layout: kernels.pblDiagnose.getBindGroupLayout(0), entries: [buffers.MI, buffers.MF, buffers.LV, buffers.S, buffers.K1, buffers.D, buffers.P, buffers.PH].map((buffer, binding) => ({ binding, resource: { buffer } })) }));
  pass.dispatchWorkgroups(Math.ceil(nC / 64));
  pass.end();
  device.queue.submit([encoder.finish()]);
  const ph = await gpu.downloadPhysics();
  layer.diagnose(state);
  let flips = 0, worstDepth = 0;
  for (let i = 0; i < nC; i++) {
    if (layer.regime[i] !== ph.REGIME[i]) flips++;
    worstDepth = Math.max(worstDepth, Math.abs(layer.depth[i] - ph.DEPTH[i]));
  }
  const regimes = [0, 0, 0, 0];
  for (let i = 0; i < nC; i++) regimes[layer.regime[i]]++;
  console.log(`${nC} columns of a surface layer under a stratocumulus: regimes (stable, surface, decoupled, coupled) ${regimes.join(', ')}; the engines' regimes differ in ${flips}, the depths by up to ${worstDepth.toFixed(2)} m`);
  assert.ok(regimes[3] > nC / 2, `${regimes[3]} coupled`);
  assert.equal(flips, 0);
  assert.ok(worstDepth < 0.5, `depth ${worstDepth}`);
});
