import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { ADJUST_KERNELS } from '../js/gpu/physics.gpu.js';
import { createModel } from '../js/model.module.js';
import { createBoundaryLayer } from '../js/physics/boundaryLayer.module.js';
import { initializeState, equilibriumProfile } from '../js/physics/init.module.js';
import { saturationHumidity } from '../js/physics/moist.module.js';
import { sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { SEA_DRAG } from '../js/physics/surface.module.js';
import { withCadence } from '../js/cadence.module.js';
import { coverOf } from '../js/ocean/iceDrift.module.js';
import { PREVIOUS_LEVELS, PREVIOUS_SURFACE } from './helpers/previousDefaults.mjs';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }

const model = createModel(new Grid(16), { levels: PREVIOUS_LEVELS, boundaryLayer: { turbulence: 'dry', diagnoseFirst: false } });
const reference = equilibriumProfile(model);
const { core, mesh, boundaryLayer } = model;
const { K, C, E, dSigma, g, geopotential } = core.diagnostics;

function column(surfaceTheta, lapseAloft, windAloft, stableAbove = Infinity) {
  const init = initializeState(model, { profile: reference });
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

const nonFinite = (a, from = 0, to = a.length) => a.subarray(from, to).filter((v) => !Number.isFinite(v)).length;
const columnTotal = (field, pi, i) => { let s = 0; for (let k = 0; k < K; k++) s += pi[i] * dSigma[k] / g * field[k * C + i]; return s; };
const enthalpy = (theta, pi, i) => { const { exnerLayer, cp } = core.diagnostics; let s = 0; for (let k = 0; k < K; k++) s += pi[i] * dSigma[k] / g * cp * theta[k * C + i] * exnerLayer[k * C + i]; return s; };

test('an unstable, sheared column gets a kilometer-deep boundary layer that mixes θ and q toward uniform and conserves its enthalpy and water', () => {
  const state = column(300, -1e-3, 8, 1000);
  const [pi, theta, , , q, qc] = state;
  boundaryLayer.diagnose(state);
  const i = 0, zb = geopotential[(K - 1) * C + i] / g;
  const depth = boundaryLayer.depth[i] - zb;
  assert.ok(depth > 500 && depth < 4000, `boundary layer ${depth} m deep`);
  let coefficients = 0;
  for (let k = boundaryLayer.kTop; k < K - 1; k++) if (boundaryLayer.mixing[k * C + i] > 0) coefficients++;
  assert.ok(coefficients >= 2, `${coefficients} interfaces mix`);
  const heatBefore = enthalpy(theta, pi, i), qBefore = columnTotal(q, pi, i);
  const spreadBefore = Math.abs(theta[(K - 1) * C + i] - theta[(K - 3) * C + i]);
  for (let n = 0; n < 24; n++) boundaryLayer.mixColumn(i, pi, theta, q, qc, 900);
  const spreadAfter = Math.abs(theta[(K - 1) * C + i] - theta[(K - 3) * C + i]);
  assert.ok(spreadAfter < 0.5 * spreadBefore, `θ spread ${spreadBefore} → ${spreadAfter}`);
  assert.ok(Math.abs(enthalpy(theta, pi, i) - heatBefore) < 1e-14 * heatBefore, `enthalpy ${heatBefore} → ${enthalpy(theta, pi, i)}`);
  assert.ok(Math.abs(columnTotal(q, pi, i) - qBefore) < 1e-10 * qBefore);
  console.log(`unstable column: boundary layer ${depth.toFixed(0)} m, ${coefficients} mixing interfaces, θ spread over the lowest three layers ${spreadBefore.toFixed(2)} → ${spreadAfter.toFixed(2)} K after 6 h`);
});

function heatedColumns(boundaryLayerOptions) {
  const pair = createModel(new Grid(16), { ocean: false, boundaryLayer: { diagnoseFirst: false, ...boundaryLayerOptions } });
  const { core: c, mesh: m, radiation, state } = pair, { K: nK, C: nC, E: nE, g: grav, geopotential: phi, exnerLayer, sigmaMid, R, cp } = c.diagnostics;
  const init = initializeState(pair, { profile: equilibriumProfile(pair) });
  for (let a = 0; a < init.length; a++) state[a].set(init[a]);
  const [pi, theta, u, , q, qc] = state;
  c.diagnose(pi, theta, q, qc);
  for (let i = 0; i < nC; i++) {
    const zb = phi[(nK - 1) * nC + i] / grav;
    for (let k = 0; k < nK; k++) {
      const x = k * nC + i, z = phi[x] / grav - zb;
      theta[x] = 300 - 1e-3 * Math.min(z, 750) + 2e-3 * Math.max(0, Math.min(z, 1500) - 750) + 6e-3 * Math.max(0, z - 1500);
      q[x] = 8e-3 * Math.exp(-z / 2000);
      qc[x] = 0;
    }
  }
  for (let k = 0; k < nK; k++) for (let e = 0; e < nE; e++) u[k * nE + e] = k === nK - 1 ? 0 : 8 * m.nEdge[3 * e + 1];
  c.diagnose(pi, theta, q, qc);
  for (let i = 0; i < nC; i++) {
    const b = (nK - 1) * nC + i, rho = pi[i] * sigmaMid[nK - 1] / (R * theta[b] * exnerLayer[b]);
    radiation.sensibleHeat[i] = 300;
    radiation.evaporation[i] = 1e-4;
    radiation.surfaceBuoyancy[i] = grav / theta[b] * (300 / (rho * cp * exnerLayer[b]) + 0.61 * theta[b] * 1e-4 / rho);
  }
  return pair;
}

test('with boundaryLayer.countergradient 7.2 an unstable, heated column carries heat and water from the lower to the upper part of its mixed layer and holds its upper half more stable, as Holtslag and Boville\'s nonlocal flux does, conserving its enthalpy and water; with 0 it is bit for bit the boundary layer built without the option', () => {
  const same = (a, b) => a.length === b.length && a.every((x, n) => Object.is(x, b[n]));
  for (const turbulence of ['dry', 'moist']) {
    const runs = [7.2, 0, undefined].map((a) => {
      const pair = heatedColumns({ turbulence, ...(a === undefined ? {} : { countergradient: a }) }), { core: c, boundaryLayer: layer, state } = pair;
      const { K: nK, C: nC, dSigma: ds, g: grav, geopotential: phi, exnerLayer, cp } = c.diagnostics, [pi, theta, , , q, qc] = state;
      const heat = (i) => { let s = 0; for (let k = 0; k < nK; k++) s += pi[i] * ds[k] / grav * cp * theta[k * nC + i] * exnerLayer[k * nC + i]; return s; };
      const water = (i) => { let s = 0; for (let k = 0; k < nK; k++) s += pi[i] * ds[k] / grav * (q[k * nC + i] + qc[k * nC + i]); return s; };
      layer.diagnose(state);
      const mixing = Float64Array.from(layer.mixing), heatBefore = Float64Array.from({ length: nC }, (_, i) => heat(i)), waterBefore = Float64Array.from({ length: nC }, (_, i) => water(i));
      const zb = phi[(nK - 1) * nC] / grav, h = layer.depth[0] - zb, z = Float64Array.from({ length: nK }, (_, k) => phi[k * nC] / grav - zb);
      let carrying = 0;
      if (layer.nonlocal) for (let k = layer.kTop; k < nK - 1; k++) if (layer.nonlocal[k * nC] > 0) carrying++;
      for (let n = 0; n < 24; n++) for (let i = 0; i < nC; i++) layer.mixColumn(i, pi, theta, q, qc, 900);
      let heatMiss = 0, waterMiss = 0;
      for (let i = 0; i < nC; i++) { heatMiss = Math.max(heatMiss, Math.abs(heat(i) - heatBefore[i]) / heatBefore[i]); waterMiss = Math.max(waterMiss, Math.abs(water(i) - waterBefore[i]) / waterBefore[i]); }
      const part = (field, from, to) => { let s = 0, w = 0; for (let k = 0; k < nK; k++) if (z[k] > from && z[k] < to) { s += ds[k] * field[k * nC]; w += ds[k]; } return s / w; };
      const nearest = (height) => { let best = 0; for (let k = 0; k < nK; k++) if (Math.abs(z[k] - height) < Math.abs(z[best] - height)) best = k; return best; };
      const [middle, upper] = [nearest(0.5 * h), nearest(h)], gradient = (theta[upper * nC] - theta[middle * nC]) / (z[upper] - z[middle]);
      const result = { pair, mixing, h, carrying, heatMiss, waterMiss, gradient, lowerT: part(theta, 0.1 * h, h / 3), upperT: part(theta, h / 3, h), lowerQ: part(q, 0.1 * h, h / 3), upperQ: part(q, h / 3, h) };
      pair.step(900);
      return result;
    });
    const [on, off, none] = runs;
    assert.ok(on.h > 1000 && on.carrying >= 3, `${turbulence}: boundary layer ${on.h} m deep, ${on.carrying} interfaces carry the nonlocal flux`);
    assert.ok(on.heatMiss < 1e-13 && on.waterMiss < 1e-10, `${turbulence}: enthalpy ${on.heatMiss}, water ${on.waterMiss} relative`);
    assert.ok(on.lowerT < off.lowerT && on.upperT > off.upperT && on.lowerQ < off.lowerQ && on.upperQ > off.upperQ, `${turbulence}: θ ${off.lowerT} → ${on.lowerT} below h/3, ${off.upperT} → ${on.upperT} above; q ${off.lowerQ} → ${on.lowerQ}, ${off.upperQ} → ${on.upperQ}`);
    assert.ok(on.gradient > off.gradient, `${turbulence}: upper half ${1000 * off.gradient} → ${1000 * on.gradient} K/km`);
    assert.ok(same(off.mixing, none.mixing), `${turbulence}: coefficients`);
    for (let s = 0; s < none.pair.state.length; s++) assert.ok(same(off.pair.state[s], none.pair.state[s]), `${turbulence}: state ${s} after the columns' mixing and a step`);
    console.log(`${turbulence} heated column: h ${on.h.toFixed(0)} m, ${on.carrying} interfaces carry the nonlocal flux; after 6 h of mixing θ ${off.lowerT.toFixed(3)} → ${on.lowerT.toFixed(3)} K between 0.1 h and h/3 and ${off.upperT.toFixed(3)} → ${on.upperT.toFixed(3)} K above, q ${(1000 * off.lowerQ).toFixed(4)} → ${(1000 * on.lowerQ).toFixed(4)} and ${(1000 * off.upperQ).toFixed(4)} → ${(1000 * on.upperQ).toFixed(4)} g/kg, the upper half ${(1000 * off.gradient).toFixed(3)} → ${(1000 * on.gradient).toFixed(3)} K/km; enthalpy kept to ${on.heatMiss.toExponential(1)}, water to ${on.waterMiss.toExponential(1)}`);
  }
});

test('under 300 W/m² for 2 h velocityScale \'holtslagBoville\' spreads the heat through the mixed layer: with or without the countergradient the lowest layer warms less, and θ differs less between 0.1 h and h, than under the capped scale, heat and water mix at least as fast as momentum there, and the column gains exactly the heat and water it is given', () => {
  const settings = [{ velocityScale: 'capped' }, { velocityScale: 'capped', countergradient: 7.2 }, { velocityScale: 'holtslagBoville' }, { velocityScale: 'holtslagBoville', countergradient: 7.2 }];
  for (const turbulence of ['dry', 'moist']) {
    const [capped, cappedCg, holtslag, holtslagCg] = settings.map((options) => {
      const pair = heatedColumns({ turbulence, ...options }), { core: c, boundaryLayer: layer, radiation, state } = pair;
      const { K: nK, C: nC, dSigma: ds, g: grav, geopotential: phi, exnerLayer, sigmaMid, R, cp } = c.diagnostics, [pi, theta, , , q, qc] = state;
      const i = 0, b = (nK - 1) * nC + i, mass = (k) => pi[i] * ds[k] / grav, dt = 900;
      const heat = () => { let s = 0; for (let k = 0; k < nK; k++) s += mass(k) * cp * theta[k * nC + i] * exnerLayer[k * nC + i]; return s; };
      const water = () => { let s = 0; for (let k = 0; k < nK; k++) s += mass(k) * (q[k * nC + i] + qc[k * nC + i]); return s; };
      const start = theta[b], heatBefore = heat(), waterBefore = water();
      let asFast = true;
      for (let n = 0; n < 8; n++) {
        const rho = pi[i] * sigmaMid[nK - 1] / (R * theta[b] * exnerLayer[b]);
        radiation.surfaceBuoyancy[i] = grav / theta[b] * (300 / (rho * cp * exnerLayer[b]) + 0.61 * theta[b] * 1e-4 / rho);
        layer.diagnose(state);
        if (layer.scalarMixing) for (let k = layer.kTop; k < nK - 1; k++) if (layer.scalarMixing[k * nC + i] < layer.mixing[k * nC + i]) asFast = false;
        theta[b] += dt * 300 / (mass(nK - 1) * cp * exnerLayer[b]);
        q[b] += dt * 1e-4 / mass(nK - 1);
        layer.mixColumn(i, pi, theta, q, qc, dt);
      }
      layer.diagnose(state);
      const zb = phi[b] / grav, h = layer.depth[i] - zb;
      let low = -1, high = -1;
      for (let k = 0; k < nK; k++) { const z = phi[k * nC + i] / grav - zb; if (z > 0.1 * h && z < h) { if (high < 0) high = k; low = k; } }
      return { warming: theta[b] - start, spread: theta[high * nC + i] - theta[low * nC + i], h, asFast, heatMiss: Math.abs(heat() - heatBefore - 8 * dt * 300) / heatBefore, waterMiss: Math.abs(water() - waterBefore - 8 * dt * 1e-4) / waterBefore };
    });
    const row = (name, r) => `${name} ${r.warming.toFixed(2)} K, θ(h) − θ(0.1 h) ${r.spread.toFixed(3)} K over h ${r.h.toFixed(0)} m`;
    console.log(`${turbulence}, 2 h under 300 W/m²: the lowest layer warms ${[['capped', capped], ['capped + countergradient', cappedCg], ['holtslagBoville', holtslag], ['holtslagBoville + countergradient', holtslagCg]].map(([name, r]) => row(name, r)).join('; ')}`);
    for (const r of [capped, cappedCg, holtslag, holtslagCg]) assert.ok(r.heatMiss < 1e-12 && r.waterMiss < 1e-10 && r.asFast, `${turbulence}: enthalpy ${r.heatMiss}, water ${r.waterMiss}, heat and water at least as fast as momentum ${r.asFast}`);
    assert.ok(holtslag.warming < 0.5 * capped.warming && holtslagCg.warming < 0.5 * cappedCg.warming, `${turbulence}: ${row('capped', capped)}; ${row('holtslagBoville', holtslag)}; with the countergradient ${cappedCg.warming} and ${holtslagCg.warming} K`);
    assert.ok(Math.abs(holtslag.spread) < 0.5 * Math.abs(capped.spread) && Math.abs(holtslagCg.spread) < 0.5 * Math.abs(cappedCg.spread), `${turbulence}: θ(h) − θ(0.1 h) ${capped.spread} capped, ${holtslag.spread} holtslagBoville; with the countergradient ${cappedCg.spread} and ${holtslagCg.spread}`);
  }
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

test('with diagnoseFirst the diagnosis is the plain one of the lowest layer without the step\'s surface heat and evaporation, and the state and diagnostics are left as deposited; withhold false diagnoses the state as it stands', () => {
  const { exnerLayer, sigmaMid, kappa, cp } = core.diagnostics, { thetaV } = core.arrays;
  const heat = Float64Array.from({ length: C }, (_, i) => 150 * Math.sin(0.9 * i)), water = Float64Array.from({ length: C }, (_, i) => 6e-5 * (1 + Math.cos(1.7 * i))), dt = 900, base = (K - 1) * C;
  const first = createBoundaryLayer(mesh, core, { turbulence: 'dry', diagnoseFirst: true, depositedHeat: heat, depositedWater: water }), plain = createBoundaryLayer(mesh, core, { turbulence: 'dry', diagnoseFirst: false });
  const state = column(300, -1e-3, 8, 1000);
  for (let i = 0; i < C; i++) state[3][i] = state[1][base + i] * exnerLayer[base + i] * Math.pow(sigmaMid[K - 1], -kappa) + 2;
  const deposited = state.map((a) => Float64Array.from(a));
  plain.diagnose(state);
  const depositedDepth = Float64Array.from(plain.depth), depositedMixing = Float64Array.from(plain.mixing), depositedV = Float64Array.from(thetaV);
  for (const step of [0, dt]) {
    first.diagnose(state, 0, C, step, { withhold: false });
    assert.deepEqual(first.depth, depositedDepth);
    assert.deepEqual(first.mixing, depositedMixing);
  }
  first.diagnose(state, 0, C, dt);
  for (const bad of [0, undefined, -dt]) assert.throws(() => first.diagnose(state, 0, C, bad), /diagnoseFirst withholds/);
  for (let a = 0; a < state.length; a++) assert.deepEqual(state[a], deposited[a]);
  assert.deepEqual(Float64Array.from(thetaV), depositedV);
  const [pi, theta, , , q] = state;
  for (let i = 0; i < C; i++) {
    theta[base + i] -= dt * heat[i] / (cp * (pi[i] * dSigma[K - 1] / g)) / exnerLayer[base + i];
    q[base + i] -= dt * water[i] * g / (pi[i] * dSigma[K - 1]);
  }
  plain.diagnose(state);
  assert.deepEqual(first.depth, plain.depth);
  assert.deepEqual(first.mixing, plain.mixing);
  let deeper = 0, shallower = 0;
  for (let i = 0; i < C; i++) { if (depositedDepth[i] > first.depth[i] + 1) { assert.ok(heat[i] > 0 || water[i] > 0); deeper++; } if (depositedDepth[i] < first.depth[i] - 1) shallower++; }
  assert.ok(deeper > C / 4, `the deposit deepens ${deeper} of ${C} columns and shallows ${shallower}`);
});

test('surfaceExcess 8.5 deepens the Richardson depth of a column over a warm surface by searching from θ_v plus Troen and Mahrt\'s excess, and leaves a column over a cool one alone', () => {
  const { exnerLayer, sigmaMid, kappa } = core.diagnostics;
  const base = (K - 1) * C, plain = createBoundaryLayer(mesh, core, { diagnoseFirst: false, turbulence: 'dry' }), excessive = createBoundaryLayer(mesh, core, { diagnoseFirst: false, turbulence: 'dry', surfaceExcess: 8.5 });
  for (const [offset, deepens] of [[3, true], [-3, false]]) {
    const state = column(300, -1e-3, 8, 1000);
    for (let i = 0; i < C; i++) state[3][i] = state[1][base + i] * exnerLayer[base + i] * Math.pow(sigmaMid[K - 1], -kappa) + offset;
    plain.diagnose(state); excessive.diagnose(state);
    let deeper = 0, gain = 0;
    for (let i = 0; i < C; i++) {
      assert.ok(excessive.depth[i] >= plain.depth[i] - 1e-9, `cell ${i}: ${excessive.depth[i]} against ${plain.depth[i]}`);
      if (excessive.depth[i] > plain.depth[i] + 1) { deeper++; gain = Math.max(gain, excessive.depth[i] - plain.depth[i]); }
    }
    if (deepens) assert.ok(deeper > C / 2 && gain < 2000, `${deeper} of ${C} deeper, by up to ${gain} m`);
    else assert.deepEqual(excessive.depth, plain.depth);
  }
});

function profileColumns(extra) {
  const pair = createModel(new Grid(16), { ocean: false, levels: sigmaInterfaces('bl34'), boundaryLayer: { diagnoseFirst: false, turbulence: 'dry', ...extra }, surface: { exchange: 'fixed', ...PREVIOUS_SURFACE } });
  const { core: c, state, boundaryLayer: layer } = pair;
  const { K: nK, C: nC, E: nE, exnerLayer, sigmaMid, kappa, geopotential: phi, g: grav } = c.diagnostics;
  const [pi, theta, u, surfaceT, q, qc] = state;
  const init = initializeState(pair, {});
  for (let a = 0; a < init.length; a++) state[a].set(init[a]);
  c.diagnose(pi, theta, q, qc);
  let seed = 4099;
  const random = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; };
  for (let i = 0; i < nC; i++) {
    const zb = phi[(nK - 1) * nC + i] / grav, surface = 286 + 14 * random(), top = 300 + 1500 * random(), below = -1e-3 + 2e-3 * random(), above = 3e-3 + 5e-3 * random(), jump = 4 * random();
    for (let k = 0; k < nK; k++) {
      const x = k * nC + i, z = phi[x] / grav - zb;
      theta[x] = z < top ? surface + below * z : surface + below * top + jump + above * (z - top);
      if (sigmaMid[k] < 0.2) theta[x] = Math.max(theta[x], init[1][x]);
      q[x] = 0.5 * saturationHumidity(theta[x] * exnerLayer[x], pi[i] * sigmaMid[k]);
      qc[x] = 0;
    }
    surfaceT[i] = theta[(nK - 1) * nC + i] * exnerLayer[(nK - 1) * nC + i] * Math.pow(sigmaMid[nK - 1], -kappa) - 3 + 6 * random();
  }
  for (let k = 0; k < nK; k++) for (let e = 0; e < nE; e++) u[k * nE + e] = k === nK - 1 ? 8 * (random() - 0.5) : 12 * (random() - 0.5);
  c.diagnose(pi, theta, q, qc);
  layer.diagnose(state);
  return pair;
}

test('the K-profile options: profileHeight \'ground\' takes the profile\'s z and h from the ground, stableProfile \'holtslagBoville\' divides a stable column\'s profile by 1 + 5 z/L (5 + z/L above z/L = 1), and prandtlCoefficient is the a of Pr', () => {
  const base = profileColumns({}), ground = profileColumns({ profileHeight: 'ground' }), stable = profileColumns({ stableProfile: 'holtslagBoville' }), plain = profileColumns({ prandtlCoefficient: 0 });
  const counter = profileColumns({ countergradient: 7.2 }), counterGround = profileColumns({ countergradient: 7.2, profileHeight: 'ground' }), capped = profileColumns({ velocityScale: 'capped' }), cappedGround = profileColumns({ velocityScale: 'capped', profileHeight: 'ground' });
  const { core: c, boundaryLayer: layer } = base, { K: nK, C: nC, g: grav, cp, geopotential: phi, exnerLayer, exnerLower } = c.diagnostics, thetaV = c.arrays.thetaV;
  const shape = (z, h) => z * (1 - z / h) ** 2, mixed = (B0, ustar, h) => Math.cbrt(ustar ** 3 + 0.6 * B0 * h);
  const prandtl = (B0, ustar, h, a) => Math.pow(1 + 1.5 * 0.4 * B0 * h / ustar ** 3, -1 / 6) + a * 0.4 * 0.1 * Math.cbrt(B0 * h) / mixed(B0, ustar, h);
  const outer = (z, h, B0, ustar) => shape(z, h) * Math.cbrt(B0 * h) / (prandtl(B0, ustar, h, 7.2) * mixed(B0, ustar, h) * h);
  const cappedScale = (B0, ustar, h) => (B0 > 0 ? ustar * Math.pow(1 - 15 * Math.max(-2, -0.1 * h * 0.4 * B0 / ustar ** 3), 0.25) : ustar);
  for (const other of [ground, stable, plain]) assert.deepEqual(other.boundaryLayer.depth, layer.depth);
  assert.deepEqual(counterGround.boundaryLayer.depth, counter.boundaryLayer.depth);
  assert.deepEqual(cappedGround.boundaryLayer.depth, capped.boundaryLayer.depth);
  const columns = { stable: 0, unstable: 0 }, worst = { ground: 0, stable: 0, prandtl: 0, untouched: 0, groundScalar: 0, counter: 0, capped: 0 }, checked = { counter: 0, capped: 0 };
  let strongest = 1, lifted = 0;
  assert.equal(nonFinite(layer.depth, 0, nC) + nonFinite(layer.buoyancyFlux, 0, nC) + nonFinite(layer.friction, 0, nC), 0, 'the depths, buoyancy fluxes and friction velocities are finite');
  for (let i = 0; i < nC; i++) {
    const b = (nK - 1) * nC + i, zb = phi[b] / grav, h = layer.depth[i] - zb, B0 = layer.buoyancyFlux[i], ustar = layer.friction[i];
    if (!(h > 0)) continue;
    const lift = cp * thetaV[b] * (exnerLower[b] - exnerLayer[b]) / grav;
    lifted = Math.max(lifted, lift);
    columns[B0 < 0 ? 'stable' : 'unstable']++;
    for (let k = layer.kTop; k < nK - 1; k++) {
      const x = k * nC + i, z = 0.5 * (phi[x] / grav + phi[x + nC] / grav) - zb, m = layer.mixing[x];
      if (z >= h) continue;
      const expected = shape(z + lift, h + lift) / shape(z, h) * (B0 > 0 ? mixed(B0, ustar, h + lift) / mixed(B0, ustar, h) : 1);
      worst.ground = Math.max(worst.ground, Math.abs(ground.boundaryLayer.mixing[x] / m / expected - 1));
      worst.groundScalar = Math.max(worst.groundScalar, Math.abs(ground.boundaryLayer.scalarMixing[x] / layer.scalarMixing[x] / expected * (B0 > 0 ? prandtl(B0, ustar, h + lift, 7.2) / prandtl(B0, ustar, h, 7.2) : 1) - 1));
      if (B0 < 0) {
        const zeta = -0.4 * z * B0 / ustar ** 3, phiStable = zeta <= 1 ? 1 + 5 * zeta : 5 + zeta;
        strongest = Math.max(strongest, phiStable);
        worst.stable = Math.max(worst.stable, Math.abs(stable.boundaryLayer.mixing[x] * phiStable / m - 1));
        worst.untouched = Math.max(worst.untouched, Math.abs(plain.boundaryLayer.mixing[x] - m), Math.abs(plain.boundaryLayer.scalarMixing[x] - layer.scalarMixing[x]));
      } else {
        worst.untouched = Math.max(worst.untouched, Math.abs(stable.boundaryLayer.mixing[x] - m), Math.abs(stable.boundaryLayer.scalarMixing[x] - layer.scalarMixing[x]), Math.abs(plain.boundaryLayer.mixing[x] - m));
        if (B0 > 0) worst.prandtl = Math.max(worst.prandtl, Math.abs(plain.boundaryLayer.scalarMixing[x] / layer.scalarMixing[x] * prandtl(B0, ustar, h, 0) / prandtl(B0, ustar, h, 7.2) - 1));
      }
    }
  }
  console.log(`${columns.stable} stable and ${columns.unstable} unstable columns with h > 0; the midpoints stand up to ${lifted.toFixed(1)} m above the ground and the stable φ reaches ${strongest.toFixed(1)}; the coefficients differ from the formulas by ${worst.ground.toExponential(1)} (ground), ${worst.stable.toExponential(1)} (stable φ) and ${worst.prandtl.toExponential(1)} (Pr without a)`);
  assert.ok(columns.stable > 20 && columns.unstable > 20 && lifted > 5 && strongest > 1.5, JSON.stringify({ columns, lifted, strongest }));
  assert.equal(worst.untouched, 0, 'stableProfile leaves unstable columns alone and prandtlCoefficient the momentum and the stable columns');
  assert.ok(worst.ground < 1e-9 && worst.stable < 1e-9 && worst.prandtl < 1e-9, JSON.stringify(worst));
  for (const [run, field] of [[counter, 'nonlocal'], [capped, 'mixing']]) {
    const { depth, buoyancyFlux, friction, kTop } = run.boundaryLayer;
    assert.equal(nonFinite(depth, 0, nC) + nonFinite(buoyancyFlux, 0, nC) + nonFinite(friction, 0, nC) + nonFinite(run.boundaryLayer[field], kTop * nC, (nK - 1) * nC), 0, `the depths, buoyancy fluxes, friction velocities and ${field} are finite`);
  }
  for (let i = 0; i < nC; i++) {
    const b = (nK - 1) * nC + i, zb = phi[b] / grav, lift = cp * thetaV[b] * (exnerLower[b] - exnerLayer[b]) / grav;
    for (const [run, framed, name] of [[counter, counterGround, 'counter'], [capped, cappedGround, 'capped']]) {
      const { depth, buoyancyFlux, friction, kTop } = run.boundaryLayer, h = depth[i] - zb, B0 = buoyancyFlux[i], ustar = friction[i], field = name === 'counter' ? 'nonlocal' : 'mixing';
      if (!(h > 0) || (name === 'counter' && !(B0 > 0))) continue;
      for (let k = kTop; k < nK - 1; k++) {
        const x = k * nC + i, z = 0.5 * (phi[x] / grav + phi[x + nC] / grav) - zb;
        if (!(z < h) || !(run.boundaryLayer[field][x] > 0) || (name === 'counter' && !(z > 0.1 * h))) continue;
        const expected = name === 'counter' ? outer(z + lift, h + lift, B0, ustar) / outer(z, h, B0, ustar) : shape(z + lift, h + lift) / shape(z, h) * cappedScale(B0, ustar, h + lift) / cappedScale(B0, ustar, h);
        worst[name] = Math.max(worst[name], Math.abs(framed.boundaryLayer[field][x] / run.boundaryLayer[field][x] / expected - 1));
        checked[name]++;
      }
    }
  }
  console.log(`ground frame against the formula: the scalars ${worst.groundScalar.toExponential(1)}, the countergradient's outer layer over ${checked.counter} interfaces ${worst.counter.toExponential(1)}, the capped scale over ${checked.capped} ${worst.capped.toExponential(1)}`);
  assert.ok(checked.counter > 100 && checked.capped > 100 && worst.groundScalar < 1e-9 && worst.counter < 1e-9 && worst.capped < 1e-9, JSON.stringify({ worst, checked }));
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

test('a warm or moist sea surface deepens the momentum mixing; a cool one keeps the neutral profile', () => {
  const { exnerLayer, sigmaMid, kappa } = core.diagnostics;
  const neutralLayer = createBoundaryLayer(mesh, core, { diagnoseFirst: false, turbulence: 'dry', stability: false }), bulkLayer = createBoundaryLayer(mesh, core, { diagnoseFirst: false, turbulence: 'dry', velocityScale: 'capped' });
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
  const deckTop = new Float64Array(C), still = { turbulence: 'dry', entrainment: { efficiency: 0, shear: 0 } }, plain = createBoundaryLayer(mesh, core, { diagnoseFirst: false, ...still }), decked = createBoundaryLayer(mesh, core, { diagnoseFirst: false, turbulence: 'dry', deckTop, ...still });
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

test('a dry stable layer over a convective boundary layer is entrained at the closure\'s rate across the first interface above h, conserving the column\'s enthalpy, water and momentum', () => {
  const state = entrainingColumn();
  const [pi, theta, u, , q, qc] = state;
  const i = 0;
  const layer = createBoundaryLayer(mesh, core, { diagnoseFirst: false, turbulence: 'dry' }), still = createBoundaryLayer(mesh, core, { diagnoseFirst: false, turbulence: 'dry', entrainment: { efficiency: 0, shear: 0 } });
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

  const dt = 900, mass = pi[i] * dSigma[above] / g, { cp } = core.diagnostics;
  const totals = () => [enthalpy(theta, pi, i), columnTotal(q, pi, i), columnTotal(qc, pi, i)];
  const staticEnergy = (x) => cp * exnerLayer[x] * theta[x] + geopotential[x];
  const before = totals(), qAbove = q[idx], energyAbove = staticEnergy(idx);
  const plainTheta = Float64Array.from(theta), plainQ = Float64Array.from(q), plainQc = Float64Array.from(qc);
  layer.mixColumn(i, pi, theta, q, qc, dt);
  const after = totals();
  for (let n = 0; n < 3; n++) assert.ok(Math.abs(after[n] - before[n]) <= 4e-16 * Math.abs(before[n]) * K, `total ${n}: ${before[n]} → ${after[n]}`);
  assert.ok(Math.abs(mass * (q[idx] - qAbove) - dt * rho * velocity * (q[idx + C] - q[idx])) < 1e-12 * mass * Math.abs(q[idx] - qAbove), 'the layer above h gains ρ w_e dt of the boundary layer\'s water');
  assert.ok(Math.abs(mass * (staticEnergy(idx) - energyAbove) - dt * rho * velocity * (staticEnergy(idx + C) - staticEnergy(idx))) < 1e-9 * mass * Math.abs(staticEnergy(idx) - energyAbove), 'and loses static energy at the same rate');
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
  console.log(`entraining column: h ${(layer.depth[i] - geopotential[(K - 1) * C + i] / g).toFixed(0)} m, B0 ${buoyancy.toExponential(2)} m²/s³, Δθv ${(jump * 300 / g).toFixed(2)} K, w_e ${(1000 * velocity).toFixed(2)} mm/s; one 900 s step moves ${(1000 * (qAbove - plainQ[idx] + q[idx] - qAbove) * mass).toFixed(1)} g/m² of water above h; column enthalpy, water and momentum change by ${((after[0] - before[0]) / before[0]).toExponential(1)}, ${((after[1] + after[2] - before[1] - before[2]) / (before[1] + before[2])).toExponential(1)}, ${((momentum() - m0) / m0).toExponential(1)}`);
});

test('a mass flux sinking through the interface above h shares the entrainment: the coefficient is max(0, ρ w_e − F) + t F with F the flux and t the core\'s interface weight, ρ w_e again where the flow rises or with entrainment.subsidence off, and w_e itself stays the closure\'s', () => {
  const state = entrainingColumn();
  const [pi, theta] = state;
  const i = 0, { piSigmaDot, exnerLayer, exnerLower, sigmaMid, R } = core.diagnostics;
  const layer = createBoundaryLayer(mesh, core, { diagnoseFirst: false, turbulence: 'dry', entrainment: { subsidence: true } }), plain = createBoundaryLayer(mesh, core, { diagnoseFirst: false, turbulence: 'dry' });
  layer.diagnose(state);
  const { above, velocity } = expectedEntrainment(state, layer, i), idx = above * C + i;
  const rho = 0.5 * (pi[i] * sigmaMid[above] / (R * theta[idx] * exnerLayer[idx]) + pi[i] * sigmaMid[above + 1] / (R * theta[idx + C] * exnerLayer[idx + C]));
  const t = (exnerLower[idx] - exnerLayer[idx]) / (exnerLayer[idx + C] - exnerLayer[idx]);
  assert.ok(t > 0.3 && t < 0.7, `interface weight ${t}`);
  const saved = Float64Array.from(piSigmaDot);
  const coefficient = (flow, which = layer) => { piSigmaDot.fill(0); piSigmaDot[(above + 1) * C + i] = flow; which.diagnose(state); return which.mixing[idx]; };
  const close = (a, b, what) => assert.ok(Math.abs(a - b) < 1e-12 * Math.abs(b), `${what}: ${a} against ${b}`);
  close(coefficient(0), rho * velocity, 'no flow');
  close(coefficient(-0.5 * g * rho * velocity), rho * velocity, 'rising flow');
  close(coefficient(0.5 * g * rho * velocity), (0.5 + 0.5 * t) * rho * velocity, 'subsidence below w_e');
  assert.ok(Math.abs(layer.entrainment[i] - velocity) < 1e-12 * velocity, `w_e ${layer.entrainment[i]} against the closure's ${velocity}`);
  close(coefficient(3 * g * rho * velocity), 3 * t * rho * velocity, 'subsidence beyond w_e');
  close(coefficient(3 * g * rho * velocity, plain), rho * velocity, 'without the option');
  piSigmaDot.set(saved);
  layer.diagnose(state);
  console.log(`subsiding interface: t ${t.toFixed(3)}, ρ w_e ${(1000 * rho * velocity).toFixed(2)} g/m²/s`);
});

test('the deck\'s opening and the stratiform share taper w_e: a half-open gate or a share of one half halves it, a closed gate or an EIS of 12 K stops it; a stable surface, zero coefficients, the cap and the jump floor bound it', () => {
  const deckGate = new Float64Array(C).fill(0.3), stratiform = new Float64Array(C);
  const layer = createBoundaryLayer(mesh, core, { diagnoseFirst: false, turbulence: 'dry', deckGate, stratiform });
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
  const capped = createBoundaryLayer(mesh, core, { diagnoseFirst: false, turbulence: 'dry', entrainment: { cap: 1e-4 } }), floored = createBoundaryLayer(mesh, core, { diagnoseFirst: false, turbulence: 'dry', entrainment: { jumpFloor: 10 } });
  const again = entrainingColumn(2), off = createBoundaryLayer(mesh, core, { diagnoseFirst: false, turbulence: 'dry', deckGate, entrainment: { efficiency: 0, shear: 0 } });
  capped.diagnose(again); floored.diagnose(again); off.diagnose(again);
  assert.ok(off.entrainment.every((x) => x === 0));
  const { buoyancy } = expectedEntrainment(again, floored, 1);
  assert.ok(buoyancy > 5e-5, 'the shear term is whole');
  assert.equal(capped.entrainment[1], 1e-4);
  assert.ok(Math.abs(floored.entrainment[1] - (0.2 * buoyancy + 5 * (Math.sqrt(SEA_DRAG) * 3) ** 3 / (floored.depth[1] - geopotential[(K - 1) * C + 1] / g)) / 10) < 1e-12);
});

test('the shear term comes in continuously with the surface buoyancy flux: w_e falls to zero as B0 falls to zero, linearly below the onset; without the onset it jumps', () => {
  const layer = createBoundaryLayer(mesh, core, { diagnoseFirst: false, turbulence: 'dry' }), switched = createBoundaryLayer(mesh, core, { diagnoseFirst: false, turbulence: 'dry', entrainment: { shearOnset: 0 } });
  const state = entrainingColumn(0), surfaceT = state[3], neutral = surfaceT[0];
  const flux = (warmth) => { surfaceT[0] = neutral + warmth; layer.diagnose(state, 0, 1); return layer.buoyancyFlux[0]; };
  let cold = -8, warm = 2;
  for (let n = 0; n < 60; n++) { const mid = 0.5 * (cold + warm); if (flux(mid) > 0) warm = mid; else cold = mid; }
  const samples = [];
  for (const target of [1e-4, 5e-5, 2.5e-5, 1e-5, 1e-6, 1e-8]) {
    let lo = cold, hi = 2;
    for (let n = 0; n < 60; n++) { const mid = 0.5 * (lo + hi); if (flux(mid) > target) hi = mid; else lo = mid; }
    flux(hi); switched.diagnose(state, 0, 1);
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

async function engines(entrainment, extra = {}) {
  const { createGpuCore } = await import('../js/gpu/core.gpu.js');
  const levels = sigmaInterfaces('bl34');
  const pair = createModel(new Grid(16), { ocean: false, levels, boundaryLayer: { diagnoseFirst: false, entrainment, turbulence: 'dry', ...extra }, surface: { exchange: 'fixed', ...PREVIOUS_SURFACE } });
  const { core: c, mesh: m, state, radiation, moist, boundaryLayer: layer } = pair;
  const { K: nK, C: nC, E: nE, exnerLayer, sigmaMid, kappa, geopotential: phi, g: grav } = c.diagnostics;
  const [pi, theta, u, surfaceT, q, qc] = state;
  const init = initializeState(pair, {});
  for (let a = 0; a < init.length; a++) state[a].set(init[a]);
  c.diagnose(pi, theta, q, qc);
  let seed = 2287;
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
  const { piSigmaDot, exnerLower } = c.diagnostics;
  for (let k = 1; k < nK; k++) for (let i = 0; i < nC; i++) piSigmaDot[k * nC + i] = Math.fround(0.6 * (random() - 0.4));
  const gpu = await createGpuCore(m, { levels, physics: { diagnoseFirst: false, entrainment, turbulence: 'dry', surfaceExchange: 'fixed', exchangeOptions: PREVIOUS_SURFACE, ...extra } });
  const { device, buffers, kernels } = gpu, dt = 900;
  gpu.upload(state);
  device.queue.writeBuffer(buffers.D, 4 * gpu.layout.D.PSD, Float32Array.from(piSigmaDot));
  device.queue.writeBuffer(buffers.D, 4 * gpu.layout.D.EXL, Float32Array.from(exnerLower));
  gpu.uploadPhysics({ mlmGate: radiation.mlmGate });
  device.queue.writeBuffer(buffers.PH, 4 * gpu.layout.PH.STRAT, Float32Array.from(radiation.stratiform));
  device.queue.writeBuffer(buffers.P, 0, Float32Array.from([dt, 0, 1, 0, 0, 0, 0, 0]));
  const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
  for (const [name, count] of [['pblDiagnose', nC], ...ADJUST_KERNELS.map((name) => [name, nC]), ['mixMomentum', nE]]) {
    pass.setPipeline(kernels[name]);
    pass.setBindGroup(0, device.createBindGroup({ layout: kernels[name].getBindGroupLayout(0), entries: [buffers.MI, buffers.MF, buffers.LV, buffers.S, buffers.K1, buffers.D, buffers.P, buffers.PH].map((buffer, binding) => ({ binding, resource: { buffer } })) }));
    pass.dispatchWorkgroups(Math.ceil(count / 64));
  }
  pass.end();
  device.queue.submit([encoder.finish()]);
  const after = await gpu.download(), ph = await gpu.downloadPhysics();
  const before = state.map((a) => Float64Array.from(a));
  layer.diagnose(state);
  const cpu = { entrainment: Float64Array.from(layer.entrainment), mixing: Float64Array.from(layer.mixing), depth: Float64Array.from(layer.depth), ...(layer.nonlocal ? { nonlocal: Float64Array.from(layer.nonlocal), kinematicHeat: Float64Array.from(layer.kinematicHeat), kinematicWater: Float64Array.from(layer.kinematicWater) } : {}), ...(layer.scalarMixing ? { scalarMixing: Float64Array.from(layer.scalarMixing) } : {}) };
  pair.phases.adjust(0, nC, dt);
  pair.phases.mixMomentum(0, nE, dt);
  return { K: nK, C: nC, E: nE, kTop: layer.kTop, gate: radiation.mlmGate, share: radiation.stratiform, before, cpu, state, after, ph, mesh: m };
}

test('the dry boundary layer with entrainment matches between the engines on a random set of columns: w_e, the interface coefficients, and θ, q, qc and the wind after the step', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const on = await engines({}), off = await engines({ efficiency: 0, shear: 0 });
  const { K: nK, C: nC, E: nE, kTop, mesh: m } = on;
  let entraining = 0, gated = 0, stable = 0, worstW = 0, scaleW = 0, worstMix = 0, flips = 0, tapered = 0;
  assert.equal(nonFinite(on.cpu.entrainment, 0, nC) + nonFinite(on.ph.ENTRAIN, 0, nC) + nonFinite(on.cpu.mixing, kTop * nC, (nK - 1) * nC) + nonFinite(on.ph.MIX, kTop * nC, (nK - 1) * nC), 0, 'w_e and the interface coefficients are finite on both engines');
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
  const parted = new Uint8Array(nC);
  let theta = 0, q = 0, qc = 0, wind = 0, moved = 0, movedWind = 0, adjusted = 0, besideParted = 0;
  for (let i = 0; i < nC; i++) {
    let columnQ = 0;
    for (let k = 0; k < nK; k++) columnQ = Math.max(columnQ, Math.abs(on.state[4][k * nC + i] - on.after[4][k * nC + i]));
    if (columnQ > 1e-5) { adjusted++; parted[i] = 1; }
    for (let k = 0; k < nK; k++) {
      const x = k * nC + i;
      moved = Math.max(moved, Math.abs(on.state[4][x] - off.state[4][x]));
      if (parted[i]) continue;
      theta = Math.max(theta, Math.abs(on.state[1][x] - on.after[1][x])); q = Math.max(q, Math.abs(on.state[4][x] - on.after[4][x])); qc = Math.max(qc, Math.abs(on.state[5][x] - on.after[5][x]));
    }
  }
  for (let e = 0; e < nE; e++) {
    const beside = parted[m.cellsOnEdge[2 * e]] || parted[m.cellsOnEdge[2 * e + 1]];
    if (beside) besideParted++;
    for (let k = 0; k < nK; k++) {
      const x = k * nE + e;
      movedWind = Math.max(movedWind, Math.abs(on.state[2][x] - off.state[2][x]));
      if (!beside) wind = Math.max(wind, Math.abs(on.state[2][x] - on.after[2][x]));
    }
  }
  let offW = 0;
  for (let i = 0; i < nC; i++) offW = Math.max(offW, off.ph.ENTRAIN[i], off.cpu.entrainment[i]);
  console.log(`${nC} random columns: ${entraining} entrain (w_e up to ${(1000 * scaleW).toFixed(1)} mm/s, ${tapered} of them tapered by the deck's opening or the stratiform share), ${gated} under a closed deck or a share of one, ${stable} over a stable surface; entrainment on or off differs between the engines on ${flips}; w_e differs by ${worstW.toExponential(1)} relative, the interface coefficients by ${worstMix.toExponential(1)} of each column's largest; after the step (but for ${adjusted} columns whose near-neutral lowest layers the dry adjustment merges in one engine only) θ by ${theta.toExponential(1)} K, q by ${q.toExponential(1)}, qc by ${qc.toExponential(1)}, the wind by ${wind.toExponential(1)} m/s off the ${besideParted} edges beside such a column, where entrainment moves q by up to ${moved.toExponential(1)} and the wind by ${movedWind.toFixed(2)} m/s`);
  assert.ok(entraining > nC / 4 && gated > nC / 10 && stable > nC / 20 && tapered > nC / 10, `${entraining} entraining, ${tapered} tapered, ${gated} gated, ${stable} stable`);
  assert.equal(offW, 0);
  assert.ok(flips <= nC / 200, `${flips} columns flip`);
  assert.ok(adjusted <= nC / 50, `${adjusted} columns part`);
  assert.ok(worstW < 1e-3 && worstMix < 1e-3, `w_e ${worstW}, coefficients ${worstMix}`);
  assert.ok(theta < 1e-3 && q < 1e-6 && qc < 1e-7 && wind < 1e-3, `θ ${theta}, q ${q}, qc ${qc}, wind ${wind}`);
  assert.ok(moved > 100 * q && movedWind > 100 * wind, `entrainment moves q by ${moved} and the wind by ${movedWind}`);
});

test('the subsidence share matches between the engines on the random columns: with entrainment.subsidence the interface coefficient max(0, ρ w_e − F) + t F agrees, and it differs from the plain coefficient where the flow sinks', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const shared = await engines({ subsidence: true }), plain = await engines({});
  const { K: nK, C: nC, kTop } = shared;
  let worstMix = 0, flips = 0, differing = 0, compared = 0;
  assert.equal(nonFinite(shared.cpu.entrainment, 0, nC) + nonFinite(shared.ph.ENTRAIN, 0, nC) + nonFinite(shared.cpu.mixing, kTop * nC, (nK - 1) * nC) + nonFinite(shared.ph.MIX, kTop * nC, (nK - 1) * nC), 0, 'w_e and the interface coefficients are finite on both engines');
  for (let i = 0; i < nC; i++) {
    const cpu = shared.cpu.entrainment[i], gpu = shared.ph.ENTRAIN[i];
    if ((cpu > 0) !== (gpu > 0)) { flips++; continue; }
    if (!(cpu > 0)) continue;
    compared++;
    let largest = 0;
    for (let k = kTop; k < nK - 1; k++) largest = Math.max(largest, shared.cpu.mixing[k * nC + i]);
    for (let k = kTop; k < nK - 1; k++) { const x = k * nC + i; if (largest > 0) worstMix = Math.max(worstMix, Math.abs(shared.cpu.mixing[x] - shared.ph.MIX[x]) / largest); }
    for (let k = kTop; k < nK - 1; k++) { const x = k * nC + i; if (Math.abs(shared.cpu.mixing[x] - plain.cpu.mixing[x]) > 1e-9 * largest) { differing++; break; } }
  }
  console.log(`subsidence share: ${compared} entraining columns, the coefficient differs from the plain one on ${differing}; engines within ${worstMix.toExponential(1)} of each column's largest coefficient, ${flips} flips`);
  assert.ok(compared > nC / 4 && differing > compared / 4, `${compared} compared, ${differing} differ`);
  assert.ok(flips <= nC / 100 && worstMix < 2e-3, `${flips} flips, coefficients ${worstMix}`);
});

test('a surface parcel that tops out at the base of a stratocumulus whose descending parcel stops there is coupled to it in both engines', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { createGpuCore } = await import('../js/gpu/core.gpu.js');
  const levels = sigmaInterfaces('bl34');
  const pair = createModel(new Grid(16), { ocean: false, levels, boundaryLayer: { diagnoseFirst: false }, surface: { exchange: 'fixed' } });
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
  const gpu = await createGpuCore(m, { levels, physics: { diagnoseFirst: false, surfaceExchange: 'fixed' } });
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
  assert.equal(nonFinite(layer.depth, 0, nC) + nonFinite(ph.DEPTH, 0, nC), 0, 'the depths are finite on both engines');
  for (let i = 0; i < nC; i++) {
    if (layer.regime[i] !== ph.REGIME[i]) { flips++; continue; }
    worstDepth = Math.max(worstDepth, Math.abs(layer.depth[i] - ph.DEPTH[i]));
  }
  const regimes = [0, 0, 0, 0];
  for (let i = 0; i < nC; i++) regimes[layer.regime[i]]++;
  console.log(`${nC} columns of a surface layer under a stratocumulus: regimes (stable, surface, decoupled, coupled) ${regimes.join(', ')}; the engines' regimes differ in ${flips}, the depths of the rest by up to ${worstDepth.toFixed(2)} m`);
  assert.ok(regimes[3] > nC / 2, `${regimes[3]} coupled`);
  assert.ok(flips <= nC / 200, `${flips} columns flip`);
  assert.ok(worstDepth < 0.5, `depth ${worstDepth}`);
});

function cumulusColumns(pair) {
  const { core: c, mesh: m, state, radiation } = pair;
  const { K: nK, C: nC, exnerLayer, sigmaMid, kappa, geopotential: phi, g: grav } = c.diagnostics;
  const [pi, theta, u, surfaceT, q, qc] = state;
  const init = initializeState(pair, {});
  for (let a = 0; a < init.length; a++) state[a].set(init[a]);
  c.diagnose(pi, theta, q, qc);
  const topLayer = new Int32Array(nC).fill(-1);
  for (let i = 0; i < nC; i++) {
    const zb = phi[(nK - 1) * nC + i] / grav, height = (k) => phi[k * nC + i] / grav - zb;
    for (let k = 0; k < nK; k++) {
      const x = k * nC + i, z = height(k), moisture = z < 1800 ? 0.0153 - 3e-6 * z : 0.003;
      const virtual = z < 1800 ? 303.3 + 1e-3 * Math.max(0, z - 300) : 303.3 + 1.5 + 6 + 6e-3 * (z - 1800);
      theta[x] = Math.max(virtual / (1 + 0.61 * moisture), init[1][x]);
      q[x] = Math.min(moisture, 0.9 * saturationHumidity(theta[x] * exnerLayer[x], pi[i] * sigmaMid[k]));
      qc[x] = z > 950 && z < 1350 ? 5e-4 : 0;
      radiation.longwave[x] = 0;
      if (qc[x] > 0 && (topLayer[i] < 0 || k < topLayer[i])) topLayer[i] = k;
    }
    radiation.longwave[topLayer[i] * nC + i] = -60;
    surfaceT[i] = theta[(nK - 1) * nC + i] * exnerLayer[(nK - 1) * nC + i] * Math.pow(sigmaMid[nK - 1], -kappa) + 0.5;
  }
  u.fill(0);
  c.diagnose(pi, theta, q, qc);
  return { topLayer, nC, nK };
}

test('a cumulus-capped column keeps its cloud-top layer decoupled: with cloudTop.cumulusCapped the descending parcel stops at the cumulus base, the surface-driven depth stays the Richardson depth and no mixing crosses the layer between', () => {
  const levels = sigmaInterfaces('bl42');
  const build = (cloudTop) => createModel(new Grid(16), { ocean: false, levels, surface: { exchange: 'fixed' }, boundaryLayer: { diagnoseFirst: false, cloudTop } });
  const plain = build({ cumulusCapped: false }), capped = build({});
  const { topLayer, nC, nK } = cumulusColumns(plain);
  cumulusColumns(capped);
  plain.boundaryLayer.diagnose(plain.state);
  capped.boundaryLayer.diagnose(capped.state);
  const { geopotential: phi, g: grav } = plain.core.diagnostics;
  let coupledPlain = 0, decoupledCapped = 0, crossing = 0, inside = 0;
  for (let i = 0; i < nC; i++) {
    const zb = phi[(nK - 1) * nC + i] / grav, interfaceZ = (k) => 0.5 * (phi[k * nC + i] + phi[(k + 1) * nC + i]) / grav - zb, cloudTop = interfaceZ(topLayer[i] - 1);
    if (plain.boundaryLayer.regime[i] === 3) coupledPlain++;
    const layer = capped.boundaryLayer, depth = layer.depth[i] - zb, base = layer.decoupling[i] - zb;
    if (layer.regime[i] === 2) decoupledCapped++;
    assert.ok(depth > 50 && depth < 800, `surface-driven depth ${depth.toFixed(0)} m`);
    assert.ok(base >= depth && base < cloudTop - 100, `decoupling height ${base.toFixed(0)} m over depth ${depth.toFixed(0)} m, cloud top ${cloudTop.toFixed(0)} m`);
    assert.ok(Math.abs(layer.mixingTop[i] - zb - cloudTop) < 1, `mixing top ${(layer.mixingTop[i] - zb).toFixed(0)} m`);
    assert.ok(layer.entrainment[i] > 0, 'the cloud-top layer entrains');
    let entrained = Infinity;
    for (let k = 0; k < nK - 1; k++) if (interfaceZ(k) >= depth) entrained = interfaceZ(k);
    for (let k = 0; k < nK - 1; k++) {
      const z = interfaceZ(k), mixed = layer.mixing[k * nC + i] > 0;
      if (z > entrained + 1 && z < base - 1) { crossing++; assert.equal(mixed, false, `mixing across the cumulus layer at ${z.toFixed(0)} m`); }
      if (z > base + 1 && z < cloudTop - 1) { inside++; assert.equal(mixed, true, `no cloud-top mixing at ${z.toFixed(0)} m`); }
    }
  }
  console.log(`${nC} cumulus-capped columns under a stratus: ${coupledPlain} coupled without the cap, ${decoupledCapped} decoupled with cumulusCapped; ${crossing} interfaces between the surface-driven depth and the decoupling height unmixed, ${inside} inside the cloud-top layer mixed`);
  assert.equal(coupledPlain, nC);
  assert.equal(decoupledCapped, nC);
  assert.ok(inside > 0);
});

test('the cumulus-capped rule matches between the engines: the cumulus columns come out decoupled at the same surface-driven depth on both', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { createGpuCore } = await import('../js/gpu/core.gpu.js');
  const levels = sigmaInterfaces('bl42');
  const pair = createModel(new Grid(16), { ocean: false, levels, surface: { exchange: 'fixed' }, boundaryLayer: { diagnoseFirst: false, cloudTop: { cumulusCapped: true } } });
  const { core: c, mesh: m, state, radiation, boundaryLayer: layer } = pair;
  const { K: nK, C: nC } = c.diagnostics;
  cumulusColumns(pair);
  for (const a of state) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  c.diagnose(state[0], state[1], state[4], state[5]);
  const gpu = await createGpuCore(m, { levels, physics: { diagnoseFirst: false, surfaceExchange: 'fixed', cloudTop: { cumulusCapped: true } } });
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
  let flips = 0, worstDepth = 0, worstTop = 0;
  const regimes = [0, 0, 0, 0];
  for (let i = 0; i < nC; i++) {
    regimes[layer.regime[i]]++;
    if (layer.regime[i] !== ph.REGIME[i]) flips++;
    worstDepth = Math.max(worstDepth, Math.abs(layer.depth[i] - ph.DEPTH[i]));
    worstTop = Math.max(worstTop, Math.abs(layer.mixingTop[i] - ph.MIXTOP[i]));
  }
  console.log(`${nC} cumulus-capped columns: regimes (stable, surface, decoupled, coupled) ${regimes.join(', ')}; the engines' regimes differ in ${flips}, the depths by up to ${worstDepth.toFixed(2)} m, the mixing tops by ${worstTop.toFixed(2)} m`);
  assert.equal(regimes[2], nC);
  assert.equal(flips, 0);
  assert.ok(worstDepth < 0.5, `depth ${worstDepth}`);
  assert.ok(worstTop < 0.5, `mixing top ${worstTop}`);
});

test('the countergradient, velocityScale \'holtslagBoville\' and the K-profile options match between the engines on the random columns, alone and together: each interface\'s ρK a w* / (w_m² h), the kinematic surface fluxes, the coefficients that mix momentum and those that mix heat and water, and θ, q, qc and the wind after the step and what each option changes in them', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const holtslag = { velocityScale: 'holtslagBoville' }, built = new Map();
  const build = async (options) => { const key = JSON.stringify(options); if (!built.has(key)) built.set(key, await engines({}, options)); return built.get(key); };
  const cases = ['dry', 'moist'].flatMap((turbulence) => [['countergradient', { velocityScale: 'capped', countergradient: 7.2, turbulence }, { velocityScale: 'capped', turbulence }], ['holtslagBoville', { ...holtslag, turbulence }, { velocityScale: 'capped', turbulence }], ['holtslagBoville with the countergradient', { ...holtslag, countergradient: 7.2, turbulence }, { ...holtslag, turbulence }], ['the K-profile from the ground with the stable φ and a Prandtl a of 0', { ...holtslag, profileHeight: 'ground', stableProfile: 'holtslagBoville', prandtlCoefficient: 0, turbulence }, { ...holtslag, turbulence }], ['the K-profile from the ground with the stable φ and the countergradient', { ...holtslag, countergradient: 7.2, profileHeight: 'ground', stableProfile: 'holtslagBoville', turbulence }, { ...holtslag, countergradient: 7.2, turbulence }]].map((c) => [turbulence, ...c]));
  for (const [turbulence, name, onOptions, baseOptions] of cases) {
    const on = await build(onOptions), base = await build(baseOptions), label = `${turbulence}, ${name}`;
    const { K: nK, C: nC, E: nE, kTop, mesh: m } = on, carries = !!on.cpu.nonlocal, scalars = !!on.cpu.scalarMixing;
    const scaleOf = (a) => (a ? a.reduce((s, x) => Math.max(s, Math.abs(x)), 0) : 0);
    const heatScale = scaleOf(on.cpu.kinematicHeat), waterScale = scaleOf(on.cpu.kinematicWater), gradientScale = scaleOf(on.cpu.nonlocal);
    let carrying = 0, differing = 0, apart = 0, flips = 0, worstG = 0, worstHeat = 0, worstWater = 0, worstMix = 0, worstScalar = 0;
    const parted = new Uint8Array(nC), interfaces = (a) => nonFinite(a, kTop * nC, (nK - 1) * nC);
    const unfinished = nonFinite(on.cpu.depth, 0, nC) + nonFinite(on.ph.DEPTH, 0, nC) + interfaces(on.cpu.mixing) + interfaces(on.ph.MIX) + (carries ? interfaces(on.cpu.nonlocal) + interfaces(on.ph.CG) + nonFinite(on.cpu.kinematicHeat, 0, nC) + nonFinite(on.ph.CGH, 0, nC) + nonFinite(on.cpu.kinematicWater, 0, nC) + nonFinite(on.ph.CGW, 0, nC) : 0) + (scalars ? interfaces(on.cpu.scalarMixing) + interfaces(on.ph.MIXH) : 0);
    assert.equal(unfinished, 0, `${label}: the depths, coefficients and countergradient fluxes are finite on both engines`);
    for (let i = 0; i < nC; i++) {
      for (const run of [on, base]) for (let k = 0; k < nK; k++) if (Math.abs(run.state[4][k * nC + i] - run.after[4][k * nC + i]) > 1e-5) parted[i] = 1;
      let largest = 0, flipped = false;
      if (carries) {
        worstHeat = Math.max(worstHeat, Math.abs(on.cpu.kinematicHeat[i] - on.ph.CGH[i]) / heatScale);
        worstWater = Math.max(worstWater, Math.abs(on.cpu.kinematicWater[i] - on.ph.CGW[i]) / waterScale);
        for (let k = kTop; k < nK - 1; k++) { const x = k * nC + i; largest = Math.max(largest, on.cpu.nonlocal[x]); if ((on.cpu.nonlocal[x] > 0) !== (on.ph.CG[x] > 0)) flipped = true; }
      }
      if (largest > 0) carrying++;
      if (scalars) for (let k = kTop; k < nK - 1; k++) if (Math.abs(on.cpu.scalarMixing[k * nC + i] - on.cpu.mixing[k * nC + i]) > 1e-6 * on.cpu.mixing[k * nC + i]) { differing++; break; }
      if (!(Math.abs(on.cpu.depth[i] - on.ph.DEPTH[i]) <= 1)) { apart++; continue; }
      if (flipped) { flips++; continue; }
      let mixLargest = 0, scalarLargest = 0;
      for (let k = kTop; k < nK - 1; k++) { mixLargest = Math.max(mixLargest, on.cpu.mixing[k * nC + i]); if (scalars) scalarLargest = Math.max(scalarLargest, on.cpu.scalarMixing[k * nC + i]); }
      for (let k = kTop; k < nK - 1; k++) {
        const x = k * nC + i;
        if (largest > 0) worstG = Math.max(worstG, Math.abs(on.cpu.nonlocal[x] - on.ph.CG[x]) / Math.max(largest, 1e-2 * gradientScale));
        if (mixLargest > 0) worstMix = Math.max(worstMix, Math.abs(on.cpu.mixing[x] - on.ph.MIX[x]) / mixLargest);
        if (scalarLargest > 0) worstScalar = Math.max(worstScalar, Math.abs(on.cpu.scalarMixing[x] - on.ph.MIXH[x]) / scalarLargest);
      }
    }
    let theta = 0, q = 0, qc = 0, wind = 0, movedT = 0, movedQ = 0, termT = 0, termQ = 0;
    for (let i = 0; i < nC; i++) {
      if (parted[i]) continue;
      for (let k = 0; k < nK; k++) {
        const x = k * nC + i;
        theta = Math.max(theta, Math.abs(on.state[1][x] - on.after[1][x])); q = Math.max(q, Math.abs(on.state[4][x] - on.after[4][x])); qc = Math.max(qc, Math.abs(on.state[5][x] - on.after[5][x]));
        movedT = Math.max(movedT, Math.abs(on.state[1][x] - base.state[1][x])); movedQ = Math.max(movedQ, Math.abs(on.state[4][x] - base.state[4][x]));
        termT = Math.max(termT, Math.abs(on.state[1][x] - base.state[1][x] - (on.after[1][x] - base.after[1][x])));
        termQ = Math.max(termQ, Math.abs(on.state[4][x] - base.state[4][x] - (on.after[4][x] - base.after[4][x])));
      }
    }
    for (let e = 0; e < nE; e++) {
      if (parted[m.cellsOnEdge[2 * e]] || parted[m.cellsOnEdge[2 * e + 1]]) continue;
      for (let k = 0; k < nK; k++) wind = Math.max(wind, Math.abs(on.state[2][k * nE + e] - on.after[2][k * nE + e]));
    }
    console.log(`${label}: ${nC} random columns, ${carrying} carry the countergradient and ${differing} mix heat and water at other rates than momentum, ${apart} have depths more than 1 m apart between the engines and ${flips} more flip an interface; the countergradient's coefficients differ by ${worstG.toExponential(1)} of each column's largest or a hundredth of the largest anywhere, the kinematic heat and water fluxes by ${worstHeat.toExponential(1)} and ${worstWater.toExponential(1)} of their largest, the momentum's coefficients by ${worstMix.toExponential(1)} and the scalars' by ${worstScalar.toExponential(1)}; after the step θ ${theta.toExponential(1)} K, q ${q.toExponential(1)}, qc ${qc.toExponential(1)}, wind ${wind.toExponential(1)} m/s apart; the option moves θ by up to ${movedT.toFixed(3)} K and q by ${movedQ.toExponential(1)}, and the engines' moves differ by ${termT.toExponential(1)} K and ${termQ.toExponential(1)}`);
    if (carries) assert.ok(carrying > nC / 4, `${label}: ${carrying} columns carry the countergradient`);
    if (scalars) assert.ok(differing > nC / 4, `${label}: ${differing} columns mix heat and water at other rates than momentum`);
    assert.ok(apart + flips <= nC / 50, `${label}: ${apart} apart, ${flips} flip`);
    assert.ok(worstG < 5e-3 && worstHeat < 1e-3 && worstWater < 1e-3 && worstMix < 1e-3 && worstScalar < 1e-3, `${label}: countergradient ${worstG}, heat ${worstHeat}, water ${worstWater}, coefficients ${worstMix}, the scalars' ${worstScalar}`);
    assert.ok(termT < 1e-3 && termQ < 1e-6 && movedT > 10 * termT && movedQ > 10 * termQ, `${label}: the option moves θ ${movedT} and q ${movedQ}, the engines' moves differ by ${termT} and ${termQ}`);
    if (turbulence === 'dry') assert.ok(theta < 1e-3 && q < 1e-6 && qc < 1e-7 && wind < 1e-3, `${label}: θ ${theta}, q ${q}, qc ${qc}, wind ${wind}`);
  }
});

function relativeWindPair() {
  const pair = createModel(new Grid(16), { ocean: false, levels: sigmaInterfaces('bl34'), boundaryLayer: { turbulence: 'dry' }, surface: { relativeWind: true } });
  const { core: c, boundaryLayer: layer, state } = pair;
  const { K: nK, C: nC, E: nE } = c.diagnostics;
  const init = initializeState(pair, {});
  for (let a = 0; a < init.length; a++) state[a].set(init[a]);
  let seed = 7919;
  const random = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; };
  for (let k = 0; k < nK; k++) for (let e = 0; e < nE; e++) state[2][k * nE + e] = k === nK - 1 ? 10 * (random() - 0.5) : 24 * (random() - 0.5);
  for (const a of state) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  layer.mixing.fill(0);
  for (let k = layer.kTop; k < nK - 1; k++) for (let i = 0; i < nC; i++) layer.mixing[k * nC + i] = Math.fround(1.5 * random());
  for (let i = 0; i < nC; i++) layer.surfaceDrag[i] = Math.fround(0.004 + 0.02 * random());
  for (let e = 0; e < nE; e++) layer.seaFlow[e] = Math.fround(random() - 0.5);
  return { pair, layer, nK, nE };
}

test('with surface.relativeWind the implicit drag pulls the lowest wind toward the sea\'s: a column moving with the sea stays as it is and takes no stress, and each edge column\'s momentum changes by the stress the sea takes', () => {
  const { pair, layer, nK, nE } = relativeWindPair();
  const [pi, , u] = pair.state, { dSigma: dS, g: grav } = pair.core.diagnostics, dt = 900;
  const uniform = new Float64Array(u.length);
  for (let k = 0; k < nK; k++) for (let e = 0; e < nE; e++) uniform[k * nE + e] = layer.seaFlow[e];
  layer.mixEdges(pi, uniform, 0, nE, dt);
  let moved = 0, held = 0;
  for (let e = 0; e < nE; e++) {
    held = Math.max(held, Math.abs(layer.surfaceStress[e]));
    for (let k = 0; k < nK; k++) moved = Math.max(moved, Math.abs(uniform[k * nE + e] - layer.seaFlow[e]));
  }
  assert.ok(moved < 1e-12 && held < 1e-12, `the column moved by ${moved} m/s and took ${held} N/m²`);
  const before = Float64Array.from(u);
  layer.mixEdges(pi, u, 0, nE, dt);
  let worst = 0, scale = 0;
  for (let e = 0; e < nE; e++) {
    const a = pair.mesh.cellsOnEdge[2 * e], b = pair.mesh.cellsOnEdge[2 * e + 1];
    let change = 0;
    for (let k = 0; k < nK; k++) change += 0.5 * (pi[a] + pi[b]) * dS[k] / grav * (u[k * nE + e] - before[k * nE + e]);
    worst = Math.max(worst, Math.abs(change + dt * layer.surfaceStress[e]));
    scale = Math.max(scale, Math.abs(dt * layer.surfaceStress[e]));
  }
  assert.ok(scale > 1 && worst < 1e-9 * scale, `momentum change against the stress the sea took: ${worst} of ${scale} kg/m/s`);
});

test('with surface.relativeWind the edge solve returns as heat the kinetic energy each column loses less the work its stress does on the moving sea, and never less than none where the sea drives the air', () => {
  const { pair, layer, nK, nE } = relativeWindPair();
  const [pi, , u] = pair.state, { dSigma: dS, g: grav } = pair.core.diagnostics, dt = 900;
  for (const [label, start] of [['calm air', new Float64Array(u.length)], ['the winds', Float64Array.from(u)]]) {
    const wind = Float64Array.from(start), heat = new Float64Array(nK * nE);
    layer.mixEdges(pi, wind, 0, nE, dt, heat);
    let worst = 0, scale = 0, lowest = Infinity;
    for (let e = 0; e < nE; e++) {
      const a = pair.mesh.cellsOnEdge[2 * e], b = pair.mesh.cellsOnEdge[2 * e + 1];
      let returned = 0, lost = 0;
      for (let k = 0; k < nK; k++) {
        const m = 0.5 * (pi[a] + pi[b]) * dS[k] / grav;
        returned += m * heat[k * nE + e];
        lost += m * (start[k * nE + e] ** 2 - wind[k * nE + e] ** 2);
      }
      const work = 2 * dt * layer.surfaceStress[e] * layer.seaFlow[e];
      worst = Math.max(worst, Math.abs(returned - (lost - work)));
      scale = Math.max(scale, Math.abs(work));
      lowest = Math.min(lowest, returned);
    }
    console.log(`${label} over the moving sea: the heat returned is the loss less the sea's work to ${worst.toExponential(2)} of ${scale.toExponential(2)}; the least an edge takes is ${lowest.toExponential(2)} (2 × J/m²)`);
    assert.ok(scale > 1 && worst < 1e-9 * scale, `${label}: ${worst} of ${scale}`);
    assert.ok(lowest >= 0, `${label}: an edge takes ${lowest}`);
  }
});

test('with surface.relativeWind a coupled CPU model hands the air the sea\'s velocity from after the ice drift\'s solve, as the GPU does', () => {
  const model = createModel(new Grid(16), { levels: sigmaInterfaces('bl36'), ...withCadence({ ocean: {} }, 1350, 16), surface: { relativeWind: true } });
  const { mesh, state, seaIce, ocean, seaFlow } = model, E = mesh.nEdges, drift = ocean.iceDrift;
  initializeState(model, {}).forEach((values, s) => state[s].set(values));
  const concentration = Float64Array.from(seaIce.concentration);
  for (let i = 0; i < mesh.nCells; i++) if (Math.abs(mesh.latCell[i]) > 62 * Math.PI / 180) { state[6][i] = 1.5; concentration[i] = 0.9; state[3][i] = Math.min(state[3][i], 271.2); }
  seaIce.load(state[6], concentration);
  ocean.initialize(state[3], state[6], { concentration: seaIce.concentration });
  const formula = () => Float64Array.from({ length: E }, (_, e) => {
    const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1], area = 0.5 * (coverOf(state[6], seaIce.concentration, a) + coverOf(state[6], seaIce.concentration, b));
    return ocean.edgeOcean[e] ? (1 - area) * ocean.u[e] + (drift.valid[e] ? area * drift.velocity[e] : 0) : 0;
  });
  const solve = drift.advance, calls = [];
  drift.advance = (given) => { const before = formula(), moved = solve(given); calls.push({ before, after: formula(), handed: Float64Array.from(seaFlow) }); return moved; };
  for (let n = 0; n < 24; n++) model.step(1350);
  const apart = (x, y) => x.reduce((w, v, e) => Math.max(w, Math.abs(v - y[e])), 0);
  const lag = Math.max(...calls.map((c) => apart(c.handed, c.after))), moved = Math.max(...calls.map((c) => apart(c.before, c.after)));
  console.log(`${calls.length} drift solves move the sea's velocity under the air by up to ${moved.toExponential(2)} m/s; the air is handed it ${lag.toExponential(2)} m/s from the post-drift one`);
  assert.ok(calls.length > 1, 'the drift ran');
  assert.ok(moved > 1e-3, `the drift's solve moves the sea's velocity by ${moved} m/s`);
  assert.ok(lag < 1e-12, `the air is handed the sea's velocity ${lag} m/s from the post-drift one`);
});

test('with relativeWind both engines mix the edge winds toward the sea\'s, hand it the same stress and return the same heat', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { createGpuCore } = await import('../js/gpu/core.gpu.js');
  const { readRanges } = await import('../js/gpu/device.module.js');
  const { pair, layer, nK, nE } = relativeWindPair();
  const gpu = await createGpuCore(pair.mesh, { levels: sigmaInterfaces('bl34'), physics: { turbulence: 'dry', relativeWind: true } });
  const { device, buffers, kernels, layout } = gpu, dt = 900, [pi, , u] = pair.state, start = Float64Array.from(u), flow = Float64Array.from(layer.seaFlow);
  const run = async (seaFlow) => {
    u.set(start);
    layer.seaFlow.set(seaFlow);
    gpu.upload(pair.state);
    gpu.uploadPhysics();
    device.queue.writeBuffer(buffers.PH, 4 * layout.PH.MIX, Float32Array.from(layer.mixing));
    device.queue.writeBuffer(buffers.PH, 4 * layout.PH.SDRAG, Float32Array.from(layer.surfaceDrag));
    device.queue.writeBuffer(buffers.PH, 4 * layout.PH.UO, Float32Array.from(seaFlow));
    device.queue.writeBuffer(buffers.P, 0, Float32Array.from([dt, 0, 1, 0, 0, 0, 0, 0]));
    device.queue.writeBuffer(buffers.D, 4 * layout.D.DISS, new Float32Array(nK * nE));
    const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
    pass.setPipeline(kernels.mixMomentum);
    pass.setBindGroup(0, device.createBindGroup({ layout: kernels.mixMomentum.getBindGroupLayout(0), entries: [buffers.MI, buffers.MF, buffers.LV, buffers.S, buffers.K1, buffers.D, buffers.P, buffers.PH].map((buffer, binding) => ({ binding, resource: { buffer } })) }));
    pass.dispatchWorkgroups(Math.ceil(nE / 64));
    pass.end();
    device.queue.submit([encoder.finish()]);
    const ph = await gpu.downloadPhysics(), [gpuHeat] = await readRanges(device, buffers.D, [{ offset: layout.D.DISS, length: nK * nE }]), cpuHeat = new Float64Array(nK * nE);
    layer.mixEdges(pi, u, 0, nE, dt, cpuHeat);
    const column = (heat) => Float64Array.from({ length: nE }, (_, e) => {
      const a = pair.mesh.cellsOnEdge[2 * e], b = pair.mesh.cellsOnEdge[2 * e + 1];
      let sum = 0;
      for (let k = 0; k < nK; k++) sum += 0.5 * (pi[a] + pi[b]) * pair.core.diagnostics.dSigma[k] / pair.core.diagnostics.g * heat[k * nE + e];
      return sum;
    });
    return { cpu: Float64Array.from(layer.surfaceStress), gpu: Float64Array.from(ph.STRESS.subarray(0, nE)), cpuHeat: column(cpuHeat), gpuHeat: column(gpuHeat) };
  };
  const moving = await run(flow), still = await run(new Float64Array(nE));
  let worst = 0, scale = 0, effect = 0;
  for (let e = 0; e < nE; e++) {
    worst = Math.max(worst, Math.abs(moving.gpu[e] - moving.cpu[e]), Math.abs(still.gpu[e] - still.cpu[e]));
    scale = Math.max(scale, Math.abs(moving.cpu[e]));
    effect = Math.max(effect, Math.abs(moving.cpu[e] - still.cpu[e]));
  }
  console.log(`relative wind on ${nE} edges: the sea's motion changes the stress by up to ${effect.toExponential(2)} of ${scale.toExponential(2)} N/m²; the engines differ by up to ${worst.toExponential(2)} N/m²`);
  assert.ok(effect > 0.05 * scale, `effect ${effect} of ${scale}`);
  assert.ok(worst < 1e-4 * scale, `engines apart by ${worst} of ${scale}`);
  let heatApart = 0, heatScale = 0, heatEffect = 0;
  for (let e = 0; e < nE; e++) {
    heatApart = Math.max(heatApart, Math.abs(moving.gpuHeat[e] - moving.cpuHeat[e]), Math.abs(still.gpuHeat[e] - still.cpuHeat[e]));
    heatScale = Math.max(heatScale, Math.abs(moving.cpuHeat[e]));
    heatEffect = Math.max(heatEffect, Math.abs(moving.cpuHeat[e] - still.cpuHeat[e]));
  }
  console.log(`the heat each edge column takes: the sea's motion changes it by up to ${heatEffect.toExponential(2)} of ${heatScale.toExponential(2)}; the engines differ by up to ${heatApart.toExponential(2)} (2 × J/m²)`);
  assert.ok(heatApart < 0.01 * heatEffect, `heat apart by ${heatApart} against the sea's effect of ${heatEffect}`);
  assert.ok(heatApart < 1e-5 * heatScale, `heat apart by ${heatApart} of ${heatScale}`);
});

test('with surface.relativeWind the coupled GPU model hands its atmosphere a finite, moving sea surface velocity', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { createGpuModel } = await import('../js/gpu/model.gpu.js');
  const { readRanges } = await import('../js/gpu/device.module.js');
  const model = await createGpuModel(new Grid(16), { surface: { relativeWind: true }, ocean: { everySteps: 2 } });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  model.load();
  model.ocean.initialize(model.state[3], model.state[6]);
  for (let n = 0; n < 8; n++) await model.step(900);
  if (model.settle) await model.settle();
  const [flow] = await readRanges(model.gpu.device, model.gpu.buffers.PH, [{ offset: model.gpu.layout.PH.UO, length: model.mesh.nEdges }]);
  let largest = 0;
  for (const v of flow) { assert.ok(Number.isFinite(v)); largest = Math.max(largest, Math.abs(v)); }
  console.log(`after 8 steps the sea's surface velocity under the air reaches ${largest.toExponential(2)} m/s`);
  assert.ok(largest > 1e-4 && largest < 2, `${largest} m/s`);
});

test('both engines refuse a surface.relativeWind or a boundaryLayer.diagnoseFirst that is not true or false: the CPU model and the GPU\'s diagnosis source always, the GPU model where webgpu is installed', async () => {
  const { boundaryLayerDiagnosis } = await import('../js/gpu/physics.gpu.js');
  const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};
  for (const [group, name] of [['surface', 'relativeWind'], ['boundaryLayer', 'diagnoseFirst']]) for (const value of [1, 'true', null]) {
    const refusal = { message: `${group}.${name} is true or false, not ${JSON.stringify(value)}` }, options = { ocean: false, [group]: { [name]: value } };
    assert.throws(() => createModel(new Grid(16), options), refusal);
    if (name === 'diagnoseFirst') assert.throws(() => boundaryLayerDiagnosis('', { diagnoseFirst: value }), refusal);
    if (createGpuModel) await assert.rejects(createGpuModel(new Grid(16), options), refusal);
  }
});
