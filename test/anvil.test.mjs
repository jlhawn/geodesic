import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { MOIST_DEFAULTS, ANVIL_FIELDS, ANVIL_TRACE, LATENT_HEAT, liquidFraction, placeAnvil, saturationHumidity } from '../js/physics/moist.module.js';
import { sounding, LEVELS, SPACING, SURFACE, SURFACE_PRESSURE, STEP } from '../explain/plumeColumn.module.js';

const LIFETIME = 3 * 3600;
const near = (a, b, tolerance) => Math.abs(a - b) <= tolerance * Math.max(Math.abs(a), Math.abs(b), 1e-30);

/*
 * Cell 0 of a 42-cell mesh shrunk to the spacing of explain/plumeColumn.module.js, on the
 * bl36 layers, with Jordan's sounding and the surface of that module's plume: a deep plume
 * whose top is near 150 hPa, ready for one call of moist.plumeColumn or moist.adjust.
 */
function deepColumn(moist = {}) {
  const grid = new Grid(2), probe = buildMesh(grid), mean = Math.sqrt(probe.areaCell.reduce((a, b) => a + b, 0) / probe.nCells);
  const model = createModel(buildMesh(grid, { radius: 6.371e6 * SPACING / mean }), { ocean: false, levels: LEVELS, moist });
  const { K, C, sigmaMid, exnerLayer, geopotential, g } = model.core.diagnostics, [pi, theta, , , q, qc] = model.state, profile = sounding();
  pi[0] = SURFACE_PRESSURE;
  for (let k = 0; k < K; k++) theta[k * C] = 300;
  for (let pass = 0; pass < 4; pass++) {
    model.core.diagnoseColumn(0, pi, theta, q, qc);
    for (let k = K - 1; k >= 0; k--) { const air = profile(geopotential[k * C] / g, pi[0] * sigmaMid[k]); theta[k * C] = air.T / exnerLayer[k * C]; q[k * C] = air.q; }
  }
  model.core.diagnoseColumn(0, pi, theta, q, qc);
  model.boundaryLayer.depth[0] = geopotential[(K - 1) * C] / g + SURFACE.depthAbove;
  model.boundaryLayer.buoyancyFlux[0] = SURFACE.buoyancyFlux;
  model.boundaryLayer.friction[0] = SURFACE.friction;
  model.radiation.mlmGate[0] = SURFACE.deckGate;
  return model;
}

/*
 * A model on the default layers whose column 0 is a 6.5 K/km atmosphere from `surface` K at
 * `humidity` of saturation, with no cloud: ready for moist.autoconvertColumn.
 */
function plainColumn(surface, humidity, moist = {}) {
  const model = createModel(new Grid(16), { ocean: false, moist });
  const { K, C, sigmaMid, exnerLayer } = model.core.diagnostics, [pi, theta, , , q, qc] = model.state;
  pi.fill(101325);
  model.core.diagnose(pi, theta);
  for (let i = 0; i < C; i++) for (let k = 0; k < K; k++) {
    const idx = k * C + i, p = pi[i] * sigmaMid[k], temperature = Math.max(200, surface - 6.5e-3 * 7000 * Math.log(pi[i] / p));
    theta[idx] = temperature / exnerLayer[idx];
    q[idx] = humidity * saturationHumidity(temperature, p);
    qc[idx] = 0;
  }
  model.core.diagnoseColumn(0, pi, theta, q, qc);
  return model;
}

function totals(model, i = 0) {
  const { K, C, dSigma, cp, g, exnerLayer } = model.core.diagnostics, [pi, theta, , , q, qc] = model.state, { anvilIce } = model.moist;
  let water = 0, enthalpy = 0;
  for (let k = 0; k < K; k++) {
    const idx = k * C + i, mass = pi[i] * dSigma[k] / g;
    water += mass * (q[idx] + qc[idx] + (k < anvilIce.length / C ? anvilIce[idx] : 0));
    enthalpy += mass * (cp * theta[idx] * exnerLayer[idx] + LATENT_HEAT * q[idx]);
  }
  return { water, enthalpy };
}

test('anvilLifetime and anvilFall are off by default and anvilSpread is 1, a model without them holds no anvil and leaves the cumulus cloud above cumulusK0 empty, and bad values are refused', () => {
  assert.equal(MOIST_DEFAULTS.anvilLifetime, null);
  assert.equal(MOIST_DEFAULTS.anvilFall, null);
  assert.equal(MOIST_DEFAULTS.anvilSpread, 1);
  const model = createModel(new Grid(16), { ocean: false });
  initializeState(model).forEach((values, s) => model.state[s].set(values));
  for (const name of ANVIL_FIELDS) assert.equal(model.moist[name].length, 0, `${name} is empty`);
  const { C } = model.core.diagnostics, K0 = model.moist.cumulusK0;
  for (let n = 0; n < 6; n++) model.step(1350);
  assert.ok(model.moist.cumulusCover.subarray(0, K0 * C).every((x) => x === 0) && model.moist.cumulusWater.subarray(0, K0 * C).every((x) => x === 0), 'no cumulus cloud above cumulusK0 without the anvil');
  assert.ok(model.moist.cumulusCover.some((x) => x > 0), 'the plumes below it have cloud');
  const on = createModel(new Grid(16), { ocean: false, moist: { anvilLifetime: LIFETIME, anvilFall: 0.5, anvilSpread: 4 } });
  assert.equal(on.moist.anvilCover.length, on.moist.cumulusK0 * on.core.diagnostics.C);
  for (const bad of [0, -3600, NaN, Infinity, '3 h']) assert.throws(() => createModel(new Grid(16), { ocean: false, moist: { anvilLifetime: bad } }), /anvilLifetime/, `anvilLifetime ${bad}`);
  for (const bad of [0, -1, NaN, '1']) assert.throws(() => createModel(new Grid(16), { ocean: false, moist: { anvilLifetime: LIFETIME, anvilFall: bad } }), /anvilFall/, `anvilFall ${bad}`);
  for (const bad of [0, -1, NaN, Infinity, '2']) assert.throws(() => createModel(new Grid(16), { ocean: false, moist: { anvilLifetime: LIFETIME, anvilSpread: bad } }), /anvilSpread/, `anvilSpread ${bad}`);
});

test('the deep plume sends the ice of the air it detrains above shallowTop to the anvil: the vapor gives up exactly that ice, the layer warms by L/c_p times it, and the cover grows by g Δt D / Δp', () => {
  const dt = 4 * STEP, on = deepColumn({ anvilLifetime: LIFETIME }), off = deepColumn();
  const { K, C, cp, g, dSigma, exnerLayer } = on.core.diagnostics, K0 = on.moist.cumulusK0;
  const rain = [on, off].map((model) => { const [pi, theta, , , q, qc] = model.state; return model.moist.plumeColumn(0, pi, theta, q, qc, dt); });
  assert.ok(on.moist.deep.deep && off.moist.deep.deep, 'both columns are deep');
  assert.equal(rain[0], rain[1], 'the plume rains alike with or without the anvil');
  const { top, baseFlux } = on.moist.deep, flux = on.moist.cumulusFlux;
  assert.ok(top < K0 - 3, `the plume reaches ${K0 - top} layers above shallowTop`);
  let formed = 0, clamped = 0;
  for (let k = 0; k < K; k++) {
    const x = k * C, q = [on.state[4][x], off.state[4][x]], theta = [on.state[1][x], off.state[1][x]];
    if (k >= K0 || k < top || !(on.moist.anvilIce[x] > 0)) {
      if (k < K0) assert.ok(on.moist.anvilIce[x] === 0 && on.moist.anvilCover[x] === 0, `layer ${k} has no anvil`);
      assert.ok(q[0] === q[1] && theta[0] === theta[1], `layer ${k} is the same with and without the anvil`);
      continue;
    }
    formed++;
    const ice = on.moist.anvilIce[x], cover = on.moist.anvilCover[x], dp = on.state[0][0] * dSigma[k];
    assert.ok(near(q[1] - q[0], ice, 1e-9), `layer ${k}: the vapor gave up ${q[1] - q[0]}, the anvil holds ${ice}`);
    assert.ok(near(theta[0] - theta[1], LATENT_HEAT * ice / (cp * exnerLayer[x]), 1e-6), `layer ${k}: the layer warms by L/c_p times the ice`);
    const rate = baseFlux * (flux[k + 1] - flux[k]) * g * dt / dp;
    assert.ok(cover > 0 && cover <= 1 && near(cover, rate, 1e-12), `layer ${k}: cover ${cover} against g Δt D / Δp = ${rate}`);
    const supply = on.moist.plumeCarried[k + 1] * (1 - liquidFraction(on.moist.plumeInterfaceT[k + 1]));
    if (near(ice / cover, supply, 1e-9)) continue;
    clamped++;
    assert.ok(ice / cover < supply && q[0] === 0, `layer ${k}: the ice is limited by the vapor only when the layer is emptied`);
  }
  console.log(`the plume tops at layer ${top} of ${K0} above shallowTop: ${formed} layers of anvil, the vapor emptied in ${clamped}, cover up to ${Math.max(...on.moist.anvilCover).toFixed(3)}, in-cloud ice up to ${(1e3 * Math.max(...Array.from(on.moist.anvilIce, (x, n) => (on.moist.anvilCover[n] > 0 ? x / on.moist.anvilCover[n] : 0)))).toFixed(3)} g/kg`);
  assert.ok(formed >= 4 && clamped === 0, `${formed} anvil layers, ${clamped} limited by the vapor`);
  const a = totals(on), b = totals(off);
  assert.ok(near(a.water, b.water, 1e-13), `total water ${a.water} against ${b.water}`);
  assert.ok(near(a.enthalpy, b.enthalpy, 1e-13), `column c_p T + L q ${a.enthalpy} against ${b.enthalpy}`);
});

test('a second pass of the plume lifts the cover by (1 − a) g Δt D / Δp and adds the ice, and the cover never passes 1', () => {
  const dt = 4 * STEP, model = deepColumn({ anvilLifetime: LIFETIME });
  const { C } = model.core.diagnostics, K0 = model.moist.cumulusK0, [pi, theta, , , q, qc] = model.state;
  model.moist.plumeColumn(0, pi, theta, q, qc, dt);
  const first = Float64Array.from(model.moist.anvilCover.subarray(0, K0 * C)), firstIce = Float64Array.from(model.moist.anvilIce.subarray(0, K0 * C));
  model.core.diagnoseColumn(0, pi, theta, q, qc);
  model.moist.plumeColumn(0, pi, theta, q, qc, dt);
  let grown = 0;
  for (let k = 0; k < K0; k++) {
    const x = k * C;
    if (!(first[x] > 0)) continue;
    assert.ok(model.moist.anvilCover[x] >= first[x] && model.moist.anvilCover[x] <= 1 && model.moist.anvilIce[x] >= firstIce[x], `layer ${k} only gains`);
    if (model.moist.anvilCover[x] > first[x]) grown++;
  }
  assert.ok(grown > 0, 'a plume that fires again adds cloud');
  for (let n = 0; n < 400; n++) { model.core.diagnoseColumn(0, pi, theta, q, qc); model.moist.plumeColumn(0, pi, theta, q, qc, 8 * dt); }
  assert.ok(model.moist.anvilCover.every((x) => x >= 0 && x <= 1), 'cover stays between 0 and 1 under 400 more passes');
});

test('anvilSpread multiplies the cover the detrained air makes, to at most 1, and changes nothing else', () => {
  const dt = 4 * STEP, plain = deepColumn({ anvilLifetime: LIFETIME }), spread = deepColumn({ anvilLifetime: LIFETIME, anvilSpread: 5 }), huge = deepColumn({ anvilLifetime: LIFETIME, anvilSpread: 1e6 });
  const { C } = plain.core.diagnostics, K0 = plain.moist.cumulusK0;
  for (const model of [plain, spread, huge]) { const [pi, theta, , , q, qc] = model.state; model.moist.plumeColumn(0, pi, theta, q, qc, dt); }
  let layers = 0;
  for (let k = 0; k < K0; k++) {
    const x = k * C;
    assert.equal(spread.moist.anvilIce[x], plain.moist.anvilIce[x], `layer ${k}: the same ice`);
    assert.equal(spread.state[4][x], plain.state[4][x]);
    assert.ok(near(spread.moist.anvilCover[x], Math.min(1, 5 * plain.moist.anvilCover[x]), 1e-14), `layer ${k}: ${spread.moist.anvilCover[x]} against five times ${plain.moist.anvilCover[x]}`);
    assert.equal(huge.moist.anvilCover[x], plain.moist.anvilCover[x] > 0 ? 1 : 0, `layer ${k}: the cover is at most 1`);
    if (plain.moist.anvilCover[x] > 0) layers++;
  }
  assert.ok(layers >= 4);
});

test('the anvil decays by exp(−Δt/anvilLifetime) in cover and ice together, the ice lost sublimating into the layer\'s vapor and cooling it by L/c_p, with water and c_p T + L q exact', () => {
  const dt = 1350, model = plainColumn(300, 0.3, { anvilLifetime: 7200, iceFall: null }), { moist } = model;
  const { C, cp, exnerLayer } = model.core.diagnostics, K0 = moist.cumulusK0, [pi, theta, , , q, qc] = model.state;
  const k = 3, x = k * C;
  assert.ok(theta[x] * exnerLayer[x] < 273.15 && k < K0);
  moist.anvilCover[x] = 0.4; moist.anvilIce[x] = 2e-5;
  const before = totals(model), vapor = q[x], temperature = theta[x] * exnerLayer[x], rain = moist.autoconvertColumn(0, pi, theta, q, qc, dt);
  const decay = Math.exp(-dt / 7200), lost = 2e-5 * (1 - decay);
  assert.equal(rain, 0);
  assert.ok(near(moist.anvilCover[x], 0.4 * decay, 1e-14) && near(moist.anvilIce[x], 2e-5 * decay, 1e-14), `cover ${moist.anvilCover[x]}, ice ${moist.anvilIce[x]}`);
  assert.ok(near(q[x] - vapor, lost, 1e-9), `vapor gained ${q[x] - vapor} of ${lost}`);
  assert.ok(near(temperature - theta[x] * exnerLayer[x], LATENT_HEAT * lost / cp, 1e-9), `cooled by ${temperature - theta[x] * exnerLayer[x]}`);
  const after = totals(model);
  assert.ok(near(after.water, before.water, 1e-13) && near(after.enthalpy, before.enthalpy, 1e-13), 'water and enthalpy are exact');
  for (let m = 0; m < K0 * C; m++) if (m !== x) assert.ok(moist.anvilIce[m] === 0 && moist.anvilCover[m] === 0, 'no other layer has an anvil');
});

test('anvil ice below ANVIL_TRACE, or cover below the cumulus trace, sublimates whole and leaves no cover', () => {
  const dt = 1350, model = plainColumn(300, 0.3, { anvilLifetime: 7200, iceFall: null }), { moist } = model;
  const { C, cp, exnerLayer } = model.core.diagnostics, [pi, theta, , , q, qc] = model.state;
  const slim = 3 * C, thin = 4 * C;
  moist.anvilCover[slim] = 0.5; moist.anvilIce[slim] = 1.1 * ANVIL_TRACE;
  moist.anvilCover[thin] = 1.2e-6; moist.anvilIce[thin] = 1e-5;
  const before = totals(model), vapor = [q[slim], q[thin]], t = [theta[slim] * exnerLayer[slim], theta[thin] * exnerLayer[thin]];
  moist.autoconvertColumn(0, pi, theta, q, qc, dt);
  const decay = Math.exp(-dt / 7200);
  assert.ok(1.1 * ANVIL_TRACE * decay < ANVIL_TRACE && 1.2e-6 * decay < 1e-6, 'both fall below their traces in the step');
  for (const [n, x, ice] of [[0, slim, 1.1 * ANVIL_TRACE], [1, thin, 1e-5]]) {
    assert.ok(moist.anvilIce[x] === 0 && moist.anvilCover[x] === 0, 'the anvil is gone');
    assert.ok(near(q[x] - vapor[n], ice, 1e-9), `all ${ice} kg/kg sublimated`);
    assert.ok(near(t[n] - theta[x] * exnerLayer[x], LATENT_HEAT * ice / cp, 1e-9));
  }
  const after = totals(model);
  assert.ok(near(after.water, before.water, 1e-13) && near(after.enthalpy, before.enthalpy, 1e-13));
});

test('the anvil\'s ice falls at iceFall (ρ q_a / a)^iceFallExponent, implicitly and top down, its cover staying and reaching the layer below at least as large, and what leaves the lowest anvil layer joins the large-scale ice stream with water and energy exact', () => {
  const dt = 1350, model = plainColumn(284, 0.3, { anvilLifetime: 7200 }), { moist } = model;
  const { C, R, sigmaMid, dSigma, exnerLayer, g } = model.core.diagnostics, K0 = moist.cumulusK0, [pi, theta, , , q, qc] = model.state;
  const upper = K0 - 3, lowest = K0 - 1, cover = 0.3, ice = 4e-5;
  assert.ok(theta[lowest * C] * exnerLayer[lowest * C] < 273.15 && theta[K0 * C] * exnerLayer[K0 * C] < 273.15, 'the layers about the stream\'s entry are colder than the melting point');
  moist.anvilCover[upper * C] = cover; moist.anvilIce[upper * C] = ice;
  moist.anvilCover[lowest * C] = 0.2; moist.anvilIce[lowest * C] = 3e-5;
  const decay = Math.exp(-dt / 7200), leaving = (k, content) => {
    const T = theta[k * C] * exnerLayer[k * C], speed = 2.5 * Math.pow(pi[0] * sigmaMid[k] / (R * T) * content, 0.16), courant = speed * dt * sigmaMid[k] * g / (R * T * dSigma[k]);
    return courant / (1 + courant);
  };
  const before = totals(model), expectedUpper = leaving(upper, ice / cover), keptUpper = ice * decay * (1 - expectedUpper), fellUpper = ice * decay * expectedUpper;
  const rained = moist.autoconvertColumn(0, pi, theta, q, qc, dt);
  const middle = (upper + 1) * C;
  assert.ok(near(moist.anvilIce[upper * C], keptUpper, 1e-12) && near(moist.anvilCover[upper * C], cover * decay, 1e-14), `the upper layer keeps ${moist.anvilIce[upper * C]} of ice, ${moist.anvilCover[upper * C]} of cover`);
  assert.ok(near(moist.anvilCover[middle], cover * decay, 1e-12) && moist.anvilIce[middle] > 0, 'the next layer has the cover of the layer above');
  assert.ok(moist.anvilIce[middle] < fellUpper * 1.0000001, 'and less ice than fell, for it falls on in the same step');
  for (let k = upper + 3; k < lowest; k++) assert.ok(moist.anvilIce[k * C] === 0 || k === lowest, `layer ${k} stays empty or is fed from above`);
  assert.ok(moist.anvilIce[lowest * C] < 3e-5 * decay + 1e-30, 'the lowest anvil layer lost ice to the stream');
  const stream = qc[K0 * C] * pi[0] * dSigma[K0] / g;
  assert.ok(stream > 0, `the stream reaches layer ${K0} as cloud ice or rain: ${qc[K0 * C]} kg/kg of cloud, ${rained} kg/m² of rain at the ground`);
  const after = totals(model);
  assert.ok(near(after.water + rained, before.water, 1e-12), `water ${after.water} + rain ${rained} against ${before.water}`);
  assert.ok(near(after.enthalpy, before.enthalpy, 1e-12), 'c_p T + L q is exact');
});

test('a speed anvilFall replaces iceFall in the anvil\'s fall and in nothing else: the anvil\'s ice falls less, and the large-scale ice falls as before', () => {
  const dt = 1350, models = [plainColumn(284, 0.3, { anvilLifetime: 1e9 }), plainColumn(284, 0.3, { anvilLifetime: 1e9, anvilFall: 0.25 })];
  const { C } = models[0].core.diagnostics, K0 = models[0].moist.cumulusK0, k = K0 - 5, x = k * C, cloudy = (K0 + 6) * C;
  for (const model of models) {
    model.state[5][cloudy] = 5e-5;
    model.moist.autoconvertColumn(0, model.state[0], model.state[1], model.state[4], model.state[5], dt);
  }
  assert.ok(models[0].state[5][cloudy] < 5e-5, 'the large-scale ice falls');
  assert.deepEqual(models[0].state[5], models[1].state[5], 'the same with either speed');
  for (const model of models) {
    model.moist.anvilCover[x] = 0.3; model.moist.anvilIce[x] = 4e-5;
    model.moist.autoconvertColumn(0, model.state[0], model.state[1], model.state[4], model.state[5], dt);
  }
  assert.ok(models[1].moist.anvilIce[x] > models[0].moist.anvilIce[x], `ice kept ${models[1].moist.anvilIce[x]} against ${models[0].moist.anvilIce[x]}`);
});

test('anvil ice at or above liquidTemperature passes whole to the layer below, melting where the large-scale stream does, with water and energy exact', () => {
  const dt = 1350, model = plainColumn(312, 0.3, { anvilLifetime: 7200 }), { moist } = model;
  const { C, exnerLayer } = model.core.diagnostics, K0 = moist.cumulusK0, [pi, theta, , , q, qc] = model.state;
  const lowest = K0 - 1, cold = 2;
  assert.ok(theta[lowest * C] * exnerLayer[lowest * C] > 273.15 && theta[cold * C] * exnerLayer[cold * C] < 273.15, 'the column is warm at the bottom of the anvil layers and cold at the top');
  moist.anvilCover[lowest * C] = 0.2; moist.anvilIce[lowest * C] = 3e-5;
  moist.anvilCover[(lowest - 1) * C] = 0.2; moist.anvilIce[(lowest - 1) * C] = 3e-5;
  const before = totals(model), vapor = [q[lowest * C], q[(lowest - 1) * C]];
  const rained = moist.autoconvertColumn(0, pi, theta, q, qc, dt);
  assert.ok(moist.anvilIce[lowest * C] === 0 && moist.anvilCover[lowest * C] === 0 && moist.anvilIce[(lowest - 1) * C] === 0, 'the warm layers hold no anvil');
  assert.ok(q[lowest * C] === vapor[0] && q[(lowest - 1) * C] === vapor[1], 'and sublimated nothing');
  const after = totals(model);
  assert.ok(near(after.water + rained, before.water, 1e-12), 'water');
  assert.ok(near(after.enthalpy, before.enthalpy, 1e-12), 'energy');
});

test('a model with the anvil on, ten steps from rest, keeps cover and ice paired between 0 and 1, hands them to the radiation as the cumulus cloud above cumulusK0, and counts the anvil in the column\'s cloud water', () => {
  const model = createModel(new Grid(16), { ocean: false, moist: { anvilLifetime: LIFETIME } });
  initializeState(model).forEach((values, s) => model.state[s].set(values));
  const { C, g, dSigma } = model.core.diagnostics, K0 = model.moist.cumulusK0, { anvilCover, anvilIce, cumulusCover, cumulusWater } = model.moist;
  let formed = 0;
  for (let n = 1; n <= 10; n++) {
    model.step(1350);
    for (let x = 0; x < K0 * C; x++) {
      assert.ok(anvilCover[x] >= 0 && anvilCover[x] <= 1 && anvilIce[x] >= 0 && (anvilCover[x] > 0) === (anvilIce[x] > 0), `step ${n}, cell ${x}: cover ${anvilCover[x]}, ice ${anvilIce[x]}`);
      assert.equal(cumulusCover[x], anvilCover[x]);
      assert.ok(anvilCover[x] > 0 ? near(cumulusWater[x], anvilIce[x] / anvilCover[x], 1e-14) : cumulusWater[x] === 0);
    }
    if (n === 10) { for (let x = 0; x < K0 * C; x++) if (anvilCover[x] > 0) formed++; }
  }
  assert.ok(formed > 0, 'anvil formed');
  const diagnostics = model.diagnostics();
  let path = 0, area = 0;
  for (let i = 0; i < C; i++) { area += model.mesh.areaCell[i]; for (let k = 0; k < K0; k++) path += model.mesh.areaCell[i] * model.state[0][i] * dSigma[k] / g * anvilIce[k * C + i]; }
  assert.ok(path > 0 && diagnostics.columnCloud >= path / area, `the cloud water ${diagnostics.columnCloud} holds the anvil's ${path / area}`);
});

test('the radiation sees the anvil: with it hidden from one step\'s radiation the planet radiates more in the columns that hold it, and the cloud is back after that step\'s adjustment', () => {
  const build = () => {
    const model = createModel(new Grid(16), { ocean: false, moist: { anvilLifetime: LIFETIME } });
    initializeState(model).forEach((values, s) => model.state[s].set(values));
    for (let n = 0; n < 10; n++) model.step(1350);
    return model;
  };
  const seen = build(), hidden = build(), C = seen.mesh.nCells, K0 = seen.moist.cumulusK0;
  assert.deepEqual(Float64Array.from(seen.state[1]), Float64Array.from(hidden.state[1]), 'two builds are the same model');
  const covered = Array.from({ length: C }, (_, i) => { let cover = 0; for (let k = 0; k < K0; k++) cover = Math.max(cover, seen.moist.anvilCover[k * C + i]); return cover; });
  for (let x = 0; x < K0 * C; x++) { hidden.moist.cumulusCover[x] = 0; hidden.moist.cumulusWater[x] = 0; }
  seen.step(1350); hidden.step(1350);
  let more = 0, total = 0, count = 0, less = 0;
  for (let i = 0; i < C; i++) {
    const d = hidden.radiation.outgoing[i] - seen.radiation.outgoing[i];
    if (covered[i] > 0.01) { count++; total += d; if (d > 0.01) more++; if (d < -1e-9) less++; }
    else assert.ok(Math.abs(d) < 1e-6 || covered[i] > 0, `column ${i} without anvil radiates alike: ${d}`);
  }
  console.log(`${count} columns with more than 1 % of anvil: ${more} radiate more than 0.01 W/m² when it is hidden, by ${(total / Math.max(1, count)).toFixed(2)} W/m² on average`);
  assert.ok(count > 10 && more > count / 3 && less === 0 && total > 0, `${more} of ${count} columns, ${less} the other way, ${total} W/m²`);
  for (let x = 0; x < K0 * C; x++) assert.equal(hidden.moist.cumulusCover[x], hidden.moist.anvilCover[x], 'the cloud is back in the cumulus arrays after the step');
});

test('placeAnvil fills a model\'s anvil from a state in the K0 × C layout, leaves it zero for no state, a state without the fields or a field of another length, and gives a CPU model\'s radiation the cloud', () => {
  const model = createModel(new Grid(16), { ocean: false, moist: { anvilLifetime: LIFETIME } });
  const { C, K } = model.core.diagnostics, K0 = model.moist.cumulusK0, length = K0 * C;
  const random = (n, seed) => { let s = seed; return Float64Array.from({ length: n }, () => { s = (Math.imul(s, 1103515245) + 12345) & 0x7fffffff; return 0.05 + s / 2147483648 * 0.4; }); };
  const saved = { anvilCover: random(length, 1), anvilIce: Float32Array.from(random(length, 2), (x) => x * 1e-4) };
  model.moist.cumulusCover.fill(9); model.moist.cumulusWater.fill(9);
  placeAnvil(model, saved);
  assert.deepEqual(model.moist.anvilCover, Float64Array.from(saved.anvilCover));
  assert.deepEqual(model.moist.anvilIce, Float64Array.from(saved.anvilIce), 'a float32 field is widened');
  for (let x = 0; x < length; x++) { assert.equal(model.moist.cumulusCover[x], model.moist.anvilCover[x]); assert.equal(model.moist.cumulusWater[x], model.moist.anvilIce[x] / model.moist.anvilCover[x]); }
  assert.ok(model.moist.cumulusCover.subarray(length).every((x) => x === 9), 'the plume layers are left to placeCumulus');
  for (const [what, state] of [['no state', null], ['a state without the fields', { pi: random(C, 3) }], ['fields one layer short', { anvilCover: random(length - C, 4), anvilIce: random(length - C, 5) }], ['K × C fields', { anvilCover: random(K * C, 6), anvilIce: random(K * C, 7) }]]) {
    placeAnvil(model, state);
    assert.ok(model.moist.anvilCover.every((x) => x === 0) && model.moist.anvilIce.every((x) => x === 0), `zero for ${what}`);
    assert.ok(model.moist.cumulusCover.subarray(0, length).every((x) => x === 0), `and no cloud for the radiation for ${what}`);
  }
  const gpuLike = { mesh: { nCells: C }, core: { K }, moist: { cumulusK0: K0, cumulusCover: new Float64Array((K - K0) * C).fill(9), cumulusWater: new Float64Array((K - K0) * C).fill(9), anvilCover: new Float64Array(length), anvilIce: new Float64Array(length) } };
  placeAnvil(gpuLike, saved);
  assert.deepEqual(gpuLike.moist.anvilCover, Float64Array.from(saved.anvilCover));
  assert.ok(gpuLike.moist.cumulusCover.every((x) => x === 9), 'a GPU model\'s plume layers are left alone');
  const plain = { mesh: { nCells: C }, core: { K }, moist: { cumulusK0: K0, cumulusCover: new Float64Array(K * C), cumulusWater: new Float64Array(K * C), anvilCover: new Float64Array(0), anvilIce: new Float64Array(0) } };
  placeAnvil(plain, saved);
  assert.ok(plain.moist.cumulusCover.every((x) => x === 0), 'a model without the anvil ignores the fields');
});

test('the column\'s anvil is unchanged by the plume when it does not fire and by the moist adjustment when the layers above shallowTop are dry and still: it only decays', () => {
  const model = plainColumn(300, 0.2, { anvilLifetime: 7200, iceFall: null }), { moist } = model;
  const { C } = model.core.diagnostics, [pi, theta, u, , q, qc] = model.state;
  const x = 3 * C;
  moist.anvilCover[x] = 0.4; moist.anvilIce[x] = 2e-5;
  moist.adjust([pi, theta, u, model.state[3], q, qc, model.state[6]], 0, 1, 900);
  assert.ok(near(moist.anvilIce[x], 2e-5 * Math.exp(-900 / 7200), 1e-9), `ice ${moist.anvilIce[x]}`);
  assert.ok(near(moist.cumulusCover[x], moist.anvilCover[x], 1e-14), 'and the radiation has what is left');
});
