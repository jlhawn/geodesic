import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createBoundaryLayer, REGIME, ENTRAINMENT_DEFAULTS, CLOUD_TOP_DEFAULTS } from '../js/physics/boundaryLayer.module.js';
import { varianceCover } from '../js/physics/radiation.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { saturationHumidity, LATENT_HEAT, R_VAPOR, liftingCondensationLevel } from '../js/physics/moist.module.js';
import { sigmaInterfaces, VIRTUAL_FACTOR } from '../js/dynamics/sigmaCore.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }

const heights = [...Array.from({ length: 21 }, (_, n) => 100 * n), 2500, 3000, 4000, 5500, 7500, 10000, 13000, 17000, 22000, 30000];
const fine = Float64Array.from([0, ...heights.slice(1).reverse().map((z) => Math.exp(-z / 8000)), 1]);
const model = createModel(new Grid(3), { levels: fine, ocean: false });
const { core, mesh, radiation } = model;
const { K, C, E, dSigma, sigmaMid, g, cp, exnerLayer, geopotential, kappa } = core.diagnostics;
const { thetaV } = core.arrays;
const L = LATENT_HEAT;
const height = (k, i) => (geopotential[k * C + i] - geopotential[(K - 1) * C + i]) / g + (geopotential[(K - 1) * C + i] - (model.surfaceGeopotential ? model.surfaceGeopotential[i] : 0)) / g;
const interfaceZ = (k, i) => 0.5 * (geopotential[k * C + i] + geopotential[(k + 1) * C + i]) / g - geopotential[(K - 1) * C + i] / g;
const mass = (pi, k, i) => pi[i] * dSigma[k] / g;

function saturate(state) {
  const [pi, theta, , , q, qc] = state;
  for (let n = 0; n < 6; n++) {
    core.diagnose(pi, theta, q, qc);
    for (let x = 0; x < K * C; x++) {
      const k = Math.floor(x / C), i = x % C, ex = exnerLayer[x], T = theta[x] * ex, qs = saturationHumidity(T, pi[i] * sigmaMid[k]);
      let change = (q[x] - qs) / (1 + L * qs * L / (R_VAPOR * T * T) / cp);
      if (change < 0) change = Math.max(change, -qc[x]);
      q[x] -= change; qc[x] += change; theta[x] += L * change / (cp * ex);
    }
  }
  core.diagnose(pi, theta, q, qc);
}

/*
 * Every cell gets the same column: below `profile` returns θ_l and q_t by
 * height, each layer saturation-adjusted from them, the sea at `sea` K and a
 * wind rising from 1 m/s at the surface to 6 m/s at 500 m.
 */
function column(profile, sea) {
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  const state = model.state, [pi, theta, u, surfaceT, q, qc] = state;
  for (let pass = 0; pass < 3; pass++) {
    core.diagnose(pi, theta, q, qc);
    for (let i = 0; i < C; i++) for (let k = 0; k < K; k++) {
      const x = k * C + i, [level, total] = profile(height(k, i));
      theta[x] = level; q[x] = total; qc[x] = 0;
    }
  }
  saturate(state);
  for (let i = 0; i < C; i++) surfaceT[i] = sea;
  for (let k = 0; k < K; k++) for (let e = 0; e < E; e++) {
    const a = mesh.cellsOnEdge[2 * e], z = height(k, a);
    u[k * E + e] = (1 + 5 * Math.min(1, z / 500)) * mesh.nEdge[3 * e + 1];
  }
  return state;
}

const INVERSION = 1300;
const stratocumulus = (z) => (z < INVERSION ? [297, 0.0119] : [305 + 4e-3 * (z - INVERSION), 0.004 * Math.exp(-(z - INVERSION) / 2500)]);
const decoupled = (z) => (z < 500 ? [297, 0.0125] : z < INVERSION ? [297 + Math.min(1, (z - 500) / 200), 0.0125] : stratocumulus(z));

function cloudTop(qc, i) {
  for (let k = K - 1; k > 0 && height(k, i) < 3000; k--) if (qc[k * C + i] > CLOUD_TOP_DEFAULTS.threshold && !(qc[(k - 1) * C + i] > CLOUD_TOP_DEFAULTS.threshold)) return k;
  return -1;
}
function cooled(state, watts) {
  const longwave = new Float64Array(K * C);
  for (let i = 0; i < C; i++) { const top = cloudTop(state[5], i); if (top >= 0) longwave[top * C + i] = -watts; }
  return longwave;
}
const liquid = (theta, qc, x) => theta[x] - L * qc[x] / (cp * exnerLayer[x]);
const totals = (state, i) => {
  const [pi, theta, , , q, qc] = state;
  let level = 0, water = 0;
  for (let k = 0; k < K; k++) { const x = k * C + i; level += mass(pi, k, i) * liquid(theta, qc, x); water += mass(pi, k, i) * (q[x] + qc[x]); }
  return [level, water];
};

function turtonNicholls(state, i, kE, lowest, h, buoyant, sheared) {
  const [pi, theta, , , q, qc] = state, { efficiency, shear, cap, jumpFloor, evaporativeEnhancement, maximumEfficiency } = ENTRAINMENT_DEFAULTS;
  let weight = 0, sumV = 0, sumL = 0, sumQ = 0;
  for (let k = kE + 1; k <= lowest; k++) { const x = k * C + i; weight += dSigma[k]; sumV += dSigma[k] * thetaV[x]; sumL += dSigma[k] * liquid(theta, qc, x); sumQ += dSigma[k] * (q[x] + qc[x]); }
  const above = thetaV[(kE - 1) * C + i] > thetaV[kE * C + i] ? (kE - 1) * C + i : kE * C + i, mean = sumV / weight;
  const top = (kE + 1) * C + i, ex = exnerLayer[top], T = theta[top] * ex, qs = saturationHumidity(T, pi[i] * sigmaMid[kE + 1]), dqs = qs * L / (R_VAPOR * T * T);
  const gamma = L / cp * dqs, c = 1 + VIRTUAL_FACTOR * q[top] - qc[top] + (1 + VIRTUAL_FACTOR) * T * dqs;
  const jumpL = liquid(theta, qc, above) - sumL / weight, jumpQ = q[above] + qc[above] - sumQ / weight, jumpV = thetaV[above] - mean;
  const saturatedJump = c / (1 + gamma) * jumpL + (c * L / (cp * ex * (1 + gamma)) - theta[top]) * jumpQ;
  const demand = dqs * ex * jumpL - jumpQ, chi = demand > 0 ? Math.min(1, qc[top] * (1 + gamma) / demand) : 1;
  const A = Math.min(maximumEfficiency, efficiency * (1 + evaporativeEnhancement * Math.max(0, chi * (1 - saturatedJump / Math.max(jumpV, jumpFloor * mean / g)))));
  return { A, chi, jumpV, velocity: Math.min(cap, (A * buoyant + shear * sheared) / (h * Math.max(g * jumpV / mean, jumpFloor))) };
}

test('a stratocumulus-topped column over a 26 °C sea is mixed from its cloud top to the surface, entrains at the Turton–Nicholls rate, conserves θ_l, q_t and momentum, and keeps a cloud of 100–300 m at its top that the variance cover makes overcast', () => {
  const state = column(stratocumulus, 299.15), [pi, theta, u, , q, qc] = state, i = 0;
  const longwave = cooled(state, 60), top = cloudTop(qc, i);
  const layer = createBoundaryLayer(mesh, core, { longwave }), clear = createBoundaryLayer(mesh, core, {});
  layer.diagnose(state); clear.diagnose(state);
  const hc = interfaceZ(top - 1, i), zb = geopotential[(K - 1) * C + i] / g;
  assert.ok(hc > 1200 && hc < 1400, `cloud top at ${hc} m above the lowest layer`);
  assert.equal(layer.regime[i], REGIME.COUPLED);
  assert.equal(layer.cloudTopCooling[i], 60);
  assert.ok(Math.abs(layer.depth[i] - zb - hc) < 1e-9, `depth ${layer.depth[i] - zb} against the cloud top ${hc}`);
  const rho = pi[i] * sigmaMid[top] / (core.diagnostics.R * theta[top * C + i] * exnerLayer[top * C + i]);
  const velocity = Math.cbrt(g / thetaV[top * C + i] * 60 / (rho * cp) * hc);
  assert.ok(Math.abs(layer.radiativeVelocity[i] - velocity) < 1e-12 * velocity, `V ${layer.radiativeVelocity[i]} against ${velocity}`);
  let added = 0, peak = 0, peakZ = 0, inside = 0;
  for (let k = layer.kTop; k < top - 1; k++) assert.equal(layer.mixing[k * C + i], 0, `nothing mixes above the inversion's interface ${k}`);
  for (let k = top; k < K - 1; k++) {
    const extra = layer.mixing[k * C + i] - clear.mixing[k * C + i], z = interfaceZ(k, i);
    assert.ok(extra > 0, `the cloud top drives interface ${k} at ${z.toFixed(0)} m`);
    inside++; added += extra;
    if (extra > peak) { peak = extra; peakZ = z; }
  }
  assert.ok(peakZ > 0.5 * hc, `the cloud-top profile peaks at ${peakZ} m of ${hc}`);
  const { A, chi, jumpV, velocity: expected } = turtonNicholls(state, i, top - 1, K - 1, hc, velocity ** 3 + layer.buoyancyFlux[i] * hc, Math.min(1, layer.buoyancyFlux[i] / 5e-5) * layer.friction[i] ** 3);
  assert.ok(layer.buoyancyFlux[i] > 0 && A > 0.2 && chi > 0, `B0 ${layer.buoyancyFlux[i]}, A ${A}, χ* ${chi}`);
  assert.ok(Math.abs(layer.entrainment[i] - expected) < 1e-12 * expected, `w_e ${layer.entrainment[i]} against ${expected}`);
  const rhoE = 0.5 * (pi[i] * sigmaMid[top - 1] / (core.diagnostics.R * theta[(top - 1) * C + i] * exnerLayer[(top - 1) * C + i]) + rho);
  assert.ok(Math.abs(layer.mixing[(top - 1) * C + i] - rhoE * expected) < 1e-12 * rhoE * expected, 'the inversion interface carries ρ w_e');

  const before = totals(state, i), dt = 300;
  layer.mixColumn(i, pi, theta, q, qc, dt);
  const after = totals(state, i);
  assert.ok(Math.abs(after[0] - before[0]) <= 2e-16 * K * before[0] && Math.abs(after[1] - before[1]) <= 2e-16 * K * before[1], `θ_l ${(after[0] - before[0]) / before[0]}, q_t ${(after[1] - before[1]) / before[1]}`);
  const e = mesh.edgesOnCell[mesh.maxEdges * i], a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1];
  const momentum = () => { let s = 0; for (let k = 0; k < K; k++) s += 0.5 * (pi[a] + pi[b]) * dSigma[k] / g * u[k * E + e]; return s; };
  const m0 = momentum(), surface0 = u[(K - 1) * E + e];
  layer.mixEdges(pi, u, e, e + 1, dt);
  assert.ok(Math.abs(momentum() - m0) <= 1e-15 * K * Math.abs(m0), `momentum ${(momentum() - m0) / m0}`);
  assert.ok(Math.abs(u[(K - 1) * E + e]) > Math.abs(surface0), 'the surface wind gains from above');

  saturate(state);
  for (let n = 0; n < 71; n++) {
    const lw = cooled(state, 60);
    longwave.set(lw);
    for (let x = 0; x < K * C; x++) theta[x] += dt * lw[x] * g / (cp * pi[x % C] * dSigma[Math.floor(x / C)] * exnerLayer[x]);
    core.diagnose(pi, theta, q, qc);
    layer.diagnose(state);
    for (let j = 0; j < C; j++) layer.mixColumn(j, pi, theta, q, qc, dt);
    saturate(state);
  }
  layer.diagnose(state);
  const newTop = cloudTop(qc, i), levels = [], totalsQ = [];
  let base = newTop;
  while (base + 1 < K && qc[(base + 1) * C + i] > 0) base++;
  for (let k = newTop; k < K; k++) { levels.push(liquid(theta, qc, k * C + i)); totalsQ.push(q[k * C + i] + qc[k * C + i]); }
  const thickness = interfaceZ(newTop - 1, i) - interfaceZ(base, i);
  const spreadL = Math.max(...levels) - Math.min(...levels), spreadQ = Math.max(...totalsQ) - Math.min(...totalsQ);
  assert.equal(layer.regime[i], REGIME.COUPLED);
  assert.ok(spreadL < 0.4 && spreadQ < 2e-4, `after 6 h θ_l spreads by ${spreadL} K and q_t by ${spreadQ} from the surface to the cloud top`);
  assert.ok(thickness >= 100 && thickness <= 300, `cloud ${thickness} m thick`);
  const mixingDepth = layer.mixingTop[i] - geopotential[(K - 1) * C + i] / g;
  radiation.column(i, pi[i], theta, state[3][i], 5, undefined, 0, q[(K - 1) * C + i], q, qc, 0.07, 0.07, 1, 1.5e-3, 0, layer.depth[i] - geopotential[(K - 1) * C + i] / g, 0, mixingDepth);
  const covers = [];
  for (let k = newTop; k <= base; k++) covers.push(radiation.layerCover[k]);
  assert.ok(Math.max(...covers) > 0.99 && radiation.layerCover[newTop] > 0.9, `cloud layers covered ${covers.map((c) => c.toFixed(3)).join(', ')}`);
  console.log(`stratocumulus column: cloud top ${hc.toFixed(0)} m, V ${velocity.toFixed(2)} m/s, cloud-top K adds ${added.toFixed(3)} kg/m²/s over ${inside} interfaces, peaking at ${peakZ.toFixed(0)} m; A ${A.toFixed(3)} (χ* ${chi.toFixed(3)}, Δθv ${jumpV.toFixed(2)} K), w_e ${(1000 * expected).toFixed(2)} mm/s; one step keeps θ_l and q_t to ${((after[0] - before[0]) / before[0]).toExponential(1)} and ${((after[1] - before[1]) / before[1]).toExponential(1)}; after 6 h θ_l spreads by ${spreadL.toFixed(3)} K and q_t by ${(1000 * spreadQ).toFixed(3)} g/kg, the cloud is ${thickness.toFixed(0)} m thick (top ${interfaceZ(newTop - 1, i).toFixed(0)} m) and covers ${covers.map((c) => c.toFixed(3)).join(', ')}`);
});

test('a clear convective column gets the surface-driven K-profile and the M21 entrainment of the dry scheme, and the same mixing', () => {
  const state = column((z) => [300 - 1e-3 * Math.min(z, 1000) + 5e-3 * Math.max(0, z - 1000), 0.01 * Math.exp(-z / 2000)], 302);
  const [pi, theta, , , q, qc] = state;
  for (let i = 0; i < C; i++) for (let k = 0; k < K; k++) if (height(k, i) < 3000) assert.equal(qc[k * C + i], 0, 'clear');
  const moistLayer = createBoundaryLayer(mesh, core, { longwave: new Float64Array(K * C), entrainment: { jumpLayers: 1 } }), dryLayer = createBoundaryLayer(mesh, core, { turbulence: 'dry' });
  const twoLayer = createBoundaryLayer(mesh, core, { longwave: new Float64Array(K * C) });
  moistLayer.diagnose(state); dryLayer.diagnose(state); twoLayer.diagnose(state);
  let worst = 0, entraining = 0, deeper = 0;
  for (let i = 0; i < C; i++) {
    assert.equal(moistLayer.regime[i], REGIME.SURFACE);
    assert.ok(Math.abs(moistLayer.depth[i] - dryLayer.depth[i]) < 1e-9);
    if (dryLayer.entrainment[i] > 0) entraining++;
    assert.ok(Math.abs(moistLayer.entrainment[i] - dryLayer.entrainment[i]) <= 1e-12 * dryLayer.entrainment[i], `cell ${i}: w_e ${moistLayer.entrainment[i]} against ${dryLayer.entrainment[i]}`);
    if (twoLayer.entrainment[i] < dryLayer.entrainment[i]) deeper++;
    for (let k = moistLayer.kTop; k < K - 1; k++) {
      const x = k * C + i, scale = Math.max(1e-30, dryLayer.mixing[x]);
      worst = Math.max(worst, Math.abs(moistLayer.mixing[x] - dryLayer.mixing[x]) / scale);
      if (twoLayer.mixing[x] !== moistLayer.mixing[x]) assert.ok(k === [...Array(K).keys()].findLast((m) => m < K - 1 && interfaceZ(m, i) >= moistLayer.depth[i] - geopotential[(K - 1) * C + i] / g), `only the entrainment interface differs, cell ${i} interface ${k}`);
    }
  }
  assert.ok(entraining === C && worst < 1e-12, `${entraining} entraining columns, coefficients differ by ${worst}`);
  assert.equal(deeper, C, 'the jump across two layers slows entrainment under a stratified free troposphere');
  const copy = state.map((a) => Float64Array.from(a));
  for (let i = 0; i < C; i++) { moistLayer.mixColumn(i, pi, theta, q, qc, 900); dryLayer.mixColumn(i, copy[0], copy[1], copy[4], copy[5], 900); }
  let theta1 = 0, q1 = 0;
  for (let x = 0; x < K * C; x++) { theta1 = Math.max(theta1, Math.abs(theta[x] - copy[1][x])); q1 = Math.max(q1, Math.abs(q[x] - copy[4][x])); }
  assert.ok(theta1 < 1e-10 && q1 < 1e-15, `θ ${theta1}, q ${q1}`);
  console.log(`clear convective column: ${C} columns surface-driven to ${(moistLayer.depth[0] - geopotential[(K - 1) * C] / g).toFixed(0)} m, coefficients within ${worst.toExponential(1)} of the dry scheme's, w_e ${(1000 * moistLayer.entrainment[0]).toFixed(2)} mm/s (${(1000 * twoLayer.entrainment[0]).toFixed(2)} with the jump across two layers); after one step θ within ${theta1.toExponential(1)} K`);
});

test('a cloud layer over a stable subcloud layer is decoupled: the cloud-top profile stops at the decoupling height above the surface-driven layer', () => {
  const state = column(decoupled, 299.15), [, , , , , qc] = state, i = 0;
  const longwave = cooled(state, 60), layer = createBoundaryLayer(mesh, core, { longwave });
  layer.diagnose(state);
  const zb = geopotential[(K - 1) * C + i] / g, top = cloudTop(qc, i);
  const hs = layer.depth[i] - zb, base = layer.decoupling[i] - zb, hc = interfaceZ(top - 1, i);
  assert.equal(layer.regime[i], REGIME.DECOUPLED);
  assert.ok(hs < base && base < hc, `surface-driven top ${hs} m, decoupling ${base} m, cloud top ${hc} m`);
  assert.ok(base >= 500 && base <= 1000, `the parcel from the cloud top stops in the stable layer, at ${base} m`);
  assert.ok(Math.abs(layer.mixingTop[i] - zb - hc) < 1e-9);
  let between = 0;
  for (let k = layer.kTop; k < K - 1; k++) {
    const z = interfaceZ(k, i), x = k * C + i;
    if (z > hs && z <= base) { assert.ok(!(layer.mixing[x] > 0) || k === [...Array(K).keys()].findLast((m) => m < K - 1 && interfaceZ(m, i) >= hs), `interface ${k} at ${z.toFixed(0)} m between the two layers mixes only by the surface layer's entrainment`); between++; }
    if (z > base && z < hc) assert.ok(layer.mixing[x] > 0, `the cloud top mixes interface ${k} at ${z.toFixed(0)} m`);
  }
  const flat = createBoundaryLayer(mesh, core, { longwave: new Float64Array(K * C) });
  flat.diagnose(state);
  assert.notEqual(flat.regime[i], REGIME.DECOUPLED, 'without cooling there is no cloud-top layer');
  console.log(`decoupled column: surface-driven to ${hs.toFixed(0)} m, the cloud-top layer from ${base.toFixed(0)} to ${hc.toFixed(0)} m (V ${layer.radiativeVelocity[i].toFixed(2)} m/s), ${between} interfaces between them; w_e at the cloud top ${(1000 * layer.entrainment[i]).toFixed(2)} mm/s`);
});

test('the variance cover is 0 for a dry layer, 1 for one saturated beyond its spread and one half at zero saturation deficit, and the radiation uses it only inside the mixing top', () => {
  const spread = 2e-4;
  assert.ok(varianceCover(-10 * spread, spread) < 1e-20, 'a dry layer');
  assert.ok(varianceCover(-3 * spread, spread) < 2e-3);
  assert.ok(Math.abs(varianceCover(0, spread) - 0.5) < 1e-8, 'zero deficit');
  assert.ok(varianceCover(10 * spread, spread) > 1 - 1e-12, 'saturated beyond its spread');
  assert.ok(Math.abs(varianceCover(spread, spread) - 0.8413447) < 2e-7, 'one spread: the Gaussian 0.8413');
  assert.equal(varianceCover(1e-9, 0), 1); assert.equal(varianceCover(-1e-9, 0), 0); assert.equal(varianceCover(0, 0), 0.5);
  const state = column(stratocumulus, 299.15), [pi, theta, , , q, qc] = state, i = 0, top = cloudTop(qc, i);
  const x = top * C + i, ex = exnerLayer[x], total = q[x] + qc[x];
  const qsl = (level) => saturationHumidity(level * ex, pi[i] * sigmaMid[top]);
  let lo = 280, hi = 320;
  for (let n = 0; n < 80; n++) { const mid = 0.5 * (lo + hi); if (qsl(mid) < total) lo = mid; else hi = mid; }
  theta[x] = 0.5 * (lo + hi) + L * qc[x] / (cp * ex);
  core.diagnose(pi, theta, q, qc);
  const depth = interfaceZ(top - 1, i) + 50;
  const cover = (mixingDepth) => { radiation.column(i, pi[i], theta, state[3][i], 5, undefined, 0, q[(K - 1) * C + i], q, qc, 0.07, 0.07, 1, 1.5e-3, 0, depth, 0, mixingDepth); return radiation.layerCover[top]; };
  const inside = cover(depth), outside = cover(-1);
  assert.ok(Math.abs(inside - 0.5) < 1e-6, `a layer at zero saturation deficit inside the mixing top covers ${inside}`);
  assert.ok(Math.abs(outside - 0.5) > 0.01, `the humidity PDF gives it ${outside}`);
  console.log(`variance cover: ${varianceCover(-3 * spread, spread).toExponential(2)} at −3σ, ${varianceCover(0, spread)} at 0, ${varianceCover(2 * spread, spread).toFixed(4)} at 2σ; a layer at zero deficit inside the mixing top ${inside.toFixed(6)}, under the PDF cover ${outside.toFixed(3)}`);
});

test('with deckRegime \'boundaryLayer\' the deck runs where the boundary layer is coupled stratocumulus, and deckBypass leaves those columns to the resolved cloud', () => {
  const run = (radiationOptions) => {
    const m = createModel(new Grid(3), { levels: sigmaInterfaces('bl34'), ocean: false, radiation: { stratusSubsidence: -1000, ...radiationOptions } });
    const init = initializeState(m, {});
    for (let a = 0; a < init.length; a++) m.state[a].set(init[a]);
    const { K: nK, C: nC, sigmaMid: mid } = m.core.diagnostics;
    for (let k = 0; k < nK; k++) if (mid[k] < 0.86) for (let i = 0; i < nC; i++) m.state[1][k * nC + i] += 8;
    m.step(600);
    for (let i = 0; i < nC; i++) m.boundaryLayer.regime[i] = i % 2 ? REGIME.COUPLED : REGIME.SURFACE;
    m.radiation.mlmGate.fill(0.5);
    m.step(600);
    return m;
  };
  const gated = run({ deckRegime: 'boundaryLayer', gateMemory: 0 }), bypassed = run({ deckRegime: 'boundaryLayer', gateMemory: 0, deckBypass: true });
  let odd = 0, even = 0, oddBypassed = 0;
  for (let i = 0; i < gated.mesh.nCells; i++) {
    if (gated.geography && gated.geography.land[i]) continue;
    if (i % 2) { if (gated.radiation.mlmGate[i] === 1) odd++; if (bypassed.radiation.mlmGate[i] === 1 && bypassed.radiation.stratusFraction[i] === 0) oddBypassed++; }
    else if (gated.radiation.mlmGate[i] === 0 && gated.radiation.stratusFraction[i] === 0) even++;
  }
  const half = Math.floor(gated.mesh.nCells / 2);
  assert.ok(odd > 0.5 * half && even > 0.5 * half, `${odd} coupled columns open the gate, ${even} surface-driven ones keep it shut`);
  assert.equal(oddBypassed, odd, 'the bypass keeps the gate and drops the mixed-layer deck');
  for (let i = 1; i < gated.mesh.nCells; i += 2) assert.equal(bypassed.radiation.mlmCover[i], 0);
});

async function moistEngines(options = {}, moist = {}) {
  const { createGpuCore } = await import('../js/gpu/core.gpu.js');
  const levels = sigmaInterfaces('bl34');
  const pair = createModel(new Grid(6), { ocean: false, levels, boundaryLayer: options, moist });
  const { core: c, mesh: m, state, radiation: r, boundaryLayer: layer } = pair;
  const { K: nK, C: nC, E: nE, exnerLayer: ex, sigmaMid: mid, kappa: kap, geopotential: phi, g: grav } = c.diagnostics;
  const [pi, theta, u, surfaceT, q, qc] = state;
  const init = initializeState(pair, {});
  for (let a = 0; a < init.length; a++) state[a].set(init[a]);
  c.diagnose(pi, theta, q, qc);
  let seed = 4049;
  const random = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; };
  const longwave = r.longwave;
  longwave.fill(0);
  for (let i = 0; i < nC; i++) {
    const zb = phi[(nK - 1) * nC + i] / grav, surface = 286 + 14 * random(), top = 400 + 1400 * random(), split = random() < 0.35 ? top * (0.3 + 0.4 * random()) : Infinity;
    const jump = 1 + 8 * random(), wet = 0.85 + 0.2 * random(), cloudy = random() < 0.7;
    for (let k = 0; k < nK; k++) {
      const x = k * nC + i, z = phi[x] / grav - zb;
      theta[x] = z < top ? surface + (z > split ? 1.5 : 0) + 2e-4 * z : surface + 2e-4 * top + jump + 4e-3 * (z - top);
      if (mid[k] < 0.2) theta[x] = Math.max(theta[x], init[1][x]);
      const qs = saturationHumidity(theta[x] * ex[x], pi[i] * mid[k]);
      q[x] = z < top ? Math.min(wet, 0.98) * qs : 0.3 * random() * qs;
      qc[x] = cloudy && z < top && z > top - 450 ? 1e-4 + 4e-4 * random() : 0;
      if (cloudy && z < top && z > top - 450) longwave[x] = -10 - 50 * random();
    }
    surfaceT[i] = theta[(nK - 1) * nC + i] * ex[(nK - 1) * nC + i] * Math.pow(mid[nK - 1], -kap) - 2 + 5 * random();
    r.mlmGate[i] = random() < 0.2 ? 0.7 : 0.3 * random();
    r.stratiform[i] = random() < 0.2 ? random() : 0;
  }
  for (let k = 0; k < nK; k++) for (let e = 0; e < nE; e++) u[k * nE + e] = k === nK - 1 ? 2 * (random() - 0.5) : 12 * (random() - 0.5);
  for (const a of [...state, longwave, r.mlmGate, r.stratiform]) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  c.diagnose(pi, theta, q, qc);
  const gpu = await createGpuCore(m, { levels, physics: { ...options, ...moist } });
  const { device, buffers, kernels, layout } = gpu, dt = 900;
  gpu.upload(state);
  gpu.uploadPhysics({ mlmGate: r.mlmGate });
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.STRAT, Float32Array.from(r.stratiform));
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.LWH, Float32Array.from(longwave));
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
  layer.diagnose(state);
  const cpu = Object.fromEntries(['entrainment', 'mixing', 'regime', 'radiativeVelocity', 'mixingTop', 'depth', 'cloudTopCooling'].map((name) => [name, Float64Array.from(layer[name])]));
  pair.phases.adjust(0, nC, dt);
  pair.phases.mixMomentum(0, nE, dt);
  return { K: nK, C: nC, E: nE, kTop: layer.kTop, cpu, state, after, ph };
}

test('the moist boundary layer matches between the engines on a random set of columns in every regime: the regime, V, w_e, the coefficients, the depths and θ, q, qc and the wind after the step, with convection vetoed under coupled stratocumulus as well', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  for (const [label, options, moist] of [['defaults', {}], ['the jump across one layer', { entrainment: { jumpLayers: 1 } }], ['the M21 taper', { entrainment: { taper: true } }], ['no surface parcel', { cloudTop: { cumulusDepth: 0 } }], ['no convection under coupled stratocumulus', {}, { coupledVeto: true }]]) {
    const run = await moistEngines(options, moist), { K: nK, C: nC, E: nE, kTop, cpu, ph } = run;
    let coupledPlumes = 0;
    for (let i = 0; i < nC; i++) if (cpu.regime[i] === REGIME.COUPLED && ph.CUMF[i] > 0) coupledPlumes++;
    if (moist) assert.equal(coupledPlumes, 0, 'no plume rises from a coupled stratocumulus-topped layer');
    else if (label === 'defaults') assert.ok(coupledPlumes > 0, `${coupledPlumes} plumes rise from coupled layers without the veto`);
    const counts = [0, 0, 0, 0];
    let flips = 0, worstV = 0, worstW = 0, worstMix = 0, worstDepth = 0;
    for (let i = 0; i < nC; i++) {
      counts[cpu.regime[i]]++;
      if (cpu.regime[i] !== ph.REGIME[i]) { flips++; continue; }
      worstV = Math.max(worstV, Math.abs(cpu.radiativeVelocity[i] - ph.VRAD[i]) / Math.max(0.1, cpu.radiativeVelocity[i]));
      worstW = Math.max(worstW, Math.abs(cpu.entrainment[i] - ph.ENTRAIN[i]) / Math.max(1e-4, cpu.entrainment[i]));
      worstDepth = Math.max(worstDepth, Math.abs(cpu.depth[i] - ph.DEPTH[i]), Math.abs(cpu.mixingTop[i] - ph.MIXTOP[i]));
      let largest = 0;
      for (let k = kTop; k < nK - 1; k++) largest = Math.max(largest, cpu.mixing[k * nC + i]);
      for (let k = kTop; k < nK - 1; k++) if (largest > 0) worstMix = Math.max(worstMix, Math.abs(cpu.mixing[k * nC + i] - ph.MIX[k * nC + i]) / largest);
    }
    let theta = 0, q = 0, qc = 0, wind = 0, adjusted = 0;
    for (let i = 0; i < nC; i++) {
      let columnQ = 0;
      for (let k = 0; k < nK; k++) columnQ = Math.max(columnQ, Math.abs(run.state[4][k * nC + i] - run.after[4][k * nC + i]));
      if (columnQ > 1e-5) { adjusted++; continue; }
      for (let k = 0; k < nK; k++) { const x = k * nC + i; theta = Math.max(theta, Math.abs(run.state[1][x] - run.after[1][x])); q = Math.max(q, Math.abs(run.state[4][x] - run.after[4][x])); qc = Math.max(qc, Math.abs(run.state[5][x] - run.after[5][x])); }
    }
    for (let x = 0; x < nK * nE; x++) wind = Math.max(wind, Math.abs(run.state[2][x] - run.after[2][x]));
    console.log(`${label}: ${nC} columns, stable/surface/decoupled/coupled ${counts.join('/')}; the regime differs on ${flips}; V by ${worstV.toExponential(1)} relative, w_e by ${worstW.toExponential(1)}, the coefficients by ${worstMix.toExponential(1)} of each column's largest, the depths by ${worstDepth.toExponential(1)} m; after the step (but for ${adjusted} columns whose near-neutral lowest layers the dry adjustment merges in one engine only) θ by ${theta.toExponential(1)} K, q by ${q.toExponential(1)}, qc by ${qc.toExponential(1)}, the wind by ${wind.toExponential(1)} m/s`);
    assert.ok(counts[REGIME.SURFACE] > nC / 20 && counts[REGIME.DECOUPLED] > nC / 20 && counts[REGIME.COUPLED] > nC / 10 && counts[REGIME.STABLE] > nC / 20, `regimes ${counts}`);
    assert.ok(flips <= nC / 100, `${flips} columns differ in regime`);
    assert.ok(worstV < 1e-4 && worstW < 2e-3 && worstMix < 2e-3 && worstDepth < 0.05, `V ${worstV}, w_e ${worstW}, coefficients ${worstMix}, depths ${worstDepth}`);
    assert.ok(adjusted <= nC / 50, `${adjusted} columns part`);
    assert.ok(theta < 2e-3 && q < 2e-6 && qc < 2e-6 && wind < 2e-3, `θ ${theta}, q ${q}, qc ${qc}, wind ${wind}`);
  }
});
