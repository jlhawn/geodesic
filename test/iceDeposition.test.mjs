import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, readdirSync } from 'node:fs';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { R_VAPOR, cloudSaturation, iceVaporPressure, saturationVaporPressure, criticalHumidityAt, uniformCover, liquidFraction, LATENT_HEAT, FUSION_HEAT, MOIST_DEFAULTS, CLEAR_AIR, depositionGrowth, placeCloudIce } from '../js/physics/moist.module.js';
import { ADJUST_KERNELS } from '../js/gpu/physics.gpu.js';

const root = fileURLToPath(new URL('..', import.meta.url));
const bl34 = sigmaInterfaces('bl34'), bl48 = sigmaInterfaces('bl48');
const build = (moist = {}, levels = bl34) => createModel(new Grid(16), { ocean: false, levels, moist });
const lapse = (p) => Math.max(210, 300 - 6.5e-3 * 7400 * Math.log(1e5 / p));
const saturated = (T, p) => cloudSaturation(T, p, true, 273.15, 235.15, {}).qs;

/*
 * Column 0 of the model at 1000 hPa at temperature(p), humidity(k, T, p)
 * of the layer's vapor and no cloud.
 */
function column(model, humidity, temperature = lapse) {
  const { K, C, sigmaMid, exnerLayer } = model.core.diagnostics;
  const [pi, theta, , , q, qc] = model.state;
  pi.fill(1e5);
  for (let pass = 0; pass < 3; pass++) {
    model.core.diagnose(pi, theta);
    for (let k = 0; k < K; k++) {
      const idx = k * C, p = pi[0] * sigmaMid[k], T = temperature(p);
      theta[idx] = T / exnerLayer[idx];
      q[idx] = humidity(k, T, p);
      qc[idx] = 0;
    }
  }
  model.core.diagnoseColumn(0, pi, theta, q, qc);
  return { pi, theta, q, qc };
}
function totals(model) {
  const { K, C, dSigma, g, cp, exnerLayer } = model.core.diagnostics, [pi, theta, , , q, qc] = model.state;
  let enthalpy = 0, water = 0;
  for (let k = 0; k < K; k++) { const idx = k * C, mass = pi[0] * dSigma[k] / g; enthalpy += (cp * theta[idx] * exnerLayer[idx] + LATENT_HEAT * q[idx]) * mass; water += (q[idx] + qc[idx]) * mass; }
  return { enthalpy, water };
}
const layerAt = (model, pressure) => { const { K, sigmaMid } = model.core.diagnostics; let best = 0; for (let k = 0; k < K; k++) if (Math.abs(sigmaMid[k] * 1e5 - pressure) < Math.abs(sigmaMid[best] * 1e5 - pressure)) best = k; return best; };
const halfWidth = (model, T, p) => { const out = cloudSaturation(T, p, true, 273.15, 235.15, {}), rhc = criticalHumidityAt(p, 1e5, 0.975, 0.75, 2); return (1 - rhc) * out.qs / (1 + LATENT_HEAT * out.slope / model.core.diagnostics.cp); };

/*
 * Every column of the model at 1000 hPa under the lapse rate to a 205 K
 * floor, its vapor 0.3 to 1.05 of saturation and cloud in about a third
 * of its layers, all at single precision; with convective, every fifth
 * column 5 K warmer below 800 hPa, at 0.98 of saturation there and 0.95
 * above, over a 305 K surface.
 */
function randomState(model, seed = 4242, convective = false) {
  const { core, mesh } = model, C = mesh.nCells, { K, sigmaMid, exnerLayer } = core.diagnostics, [pi, theta, , surfaceT, q, qc] = model.state;
  const random = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; };
  for (let pass = 0; pass < 3; pass++) {
    core.diagnose(pi.fill(1e5), theta, q, qc);
    for (let i = 0; i < C; i++) for (let k = 0; k < K; k++) {
      const x = k * C + i, p = 1e5 * sigmaMid[k], deep = convective && i % 5 === 0, T = Math.max(205, 300 - 6.5e-3 * 7400 * Math.log(1e5 / p)) + (deep && p > 8e4 ? 5 : 0);
      theta[x] = T / exnerLayer[x];
      if (pass === 0) { q[x] = (0.3 + 0.75 * random()) * saturated(T, p); qc[x] = random() < 0.3 ? 2e-4 * random() * (T < 250 ? 0.2 : 1) : 0; }
      if (deep) q[x] = (p > 8e4 ? 0.98 : 0.95) * saturated(T, p);
    }
  }
  surfaceT.fill(300);
  if (convective) for (let i = 0; i < C; i += 5) surfaceT[i] = 305;
  for (const a of model.state) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  core.diagnose(pi, theta, q, qc);
  return random;
}

test("iceDeposition is off by default with detrainedPhase 'ramp'; it needs iceFall, and detrainedPhase is 'ramp' or 'liquid', the latter only under it", () => {
  assert.equal(MOIST_DEFAULTS.iceDeposition, false);
  assert.equal(MOIST_DEFAULTS.detrainedPhase, 'ramp');
  assert.equal(build().moist.cloudIceShare, null);
  assert.ok(build({ iceDeposition: true }).moist.cloudIceShare instanceof Float64Array);
  assert.throws(() => build({ iceDeposition: true, iceFall: null }), /needs iceFall/);
  assert.throws(() => build({ iceDeposition: true, detrainedPhase: 'ice' }), /detrainedPhase must be/);
  assert.throws(() => build({ detrainedPhase: 'liquid' }), /needs iceDeposition/);
});

test('under iceDeposition a column whose cloud is all colder than 235.15 K or all warmer than 273.15 K steps as without it, bit for bit, and stores a glaciated fraction of 1 or 0', () => {
  for (const T of [225, 285]) {
    const run = (moist) => {
      const model = build(moist), { K, C } = model.core.diagnostics;
      const s = column(model, (j, Tj, pj) => (j % 3 === 0 ? 1.02 : 0.6) * saturated(Tj, pj), () => T);
      for (let k = 0; k < K; k += 3) s.qc[k * C] = 5e-5 * (1 + k / K);
      model.core.diagnoseColumn(0, s.pi, s.theta, s.q, s.qc);
      placeCloudIce(model, null);
      model.moist.precipitation.fill(0);
      model.moist.adjust(model.state, 0, 1, 900);
      const layers = (a) => Array.from({ length: K }, (_, k) => a[k * C]);
      return { rain: model.moist.precipitation[0], theta: layers(s.theta), q: layers(s.q), qc: layers(s.qc), share: model.moist.cloudIceShare && layers(model.moist.cloudIceShare) };
    };
    const plain = run({}), deposition = run({ iceDeposition: true });
    assert.deepEqual({ ...deposition, share: null }, { ...plain, share: null });
    assert.deepEqual(deposition.share, deposition.qc.map((c) => (c > CLEAR_AIR && T < 235.15 ? 1 : 0)));
    console.log(`${T} K: ${(1e3 * deposition.qc.reduce((sum, x) => sum + x, 0)).toFixed(4)} g/kg of cloud summed over the layers and ${(1e3 * deposition.rain).toExponential(3)} g/m² of precipitation, both as without iceDeposition`);
  }
});

test("in a mixed-phase cloud-top layer deposition grows the cloud ice from its seed by Rotstayn's form over the step, at cloud top by f_nuc + 0.1 (1 − f_nuc) of it; only that ice falls, the clear layer below keeps all it takes as ice, and the liquid converts at the plain rates whatever bergeronAutoconversion says", () => {
  const dt = 600, run = (moist) => {
    const model = build(moist), { C, sigmaMid, exnerLayer } = model.core.diagnostics, k = layerAt(model, 4e4);
    const s = column(model, (j, T, p) => (j === k ? 1 : 0.3) * saturated(T, p));
    s.qc[k * C] = 1e-4;
    model.core.diagnoseColumn(0, s.pi, s.theta, s.q, s.qc);
    model.moist.cloudIceShare.fill(0);
    const before = totals(model), T = s.theta[k * C] * exnerLayer[k * C], p = s.pi[0] * sigmaMid[k];
    const rain = model.moist.autoconvertColumn(0, s.pi, s.theta, s.q, s.qc, dt);
    return { model, k, s, before, T, p, rain };
  };
  const { model, k, s, before, T, p, rain } = run({ iceDeposition: true, bergeronAutoconversion: true }), { C, sigmaMid, dSigma, g, R } = model.core.diagnostics, ice = model.moist.cloudIceShare;
  assert.ok(T > 235.15 && T < 268.16, `the cloud layer at ${T} K is mixed and under the Bergeron factor's range`);
  const icePressure = iceVaporPressure(T), excess = (saturationVaporPressure(T) - icePressure) / icePressure, nuclei = 1000 * Math.exp(12.96 * excess - 0.639), density = p / (R * T), sublimation = LATENT_HEAT + FUSION_HEAT;
  const growth = 7.8 * Math.pow(nuclei / density, 2 / 3) * excess / (Math.cbrt(700) * (sublimation / (0.024 * T) * (sublimation / (R_VAPOR * T) - 1) + R_VAPOR * T * p / (2.21 * icePressure)));
  const cover = uniformCover(1e-4, halfWidth(model, T, p)), seed = 1e-12 * nuclei / density, active = Math.min(nuclei / 15000, 1);
  const grown = (active + (1 - active) * 0.1) * Math.min(cover * (Math.pow(2 / 3 * growth * dt + Math.pow(seed, 2 / 3), 1.5) - seed), 1e-4);
  const speed = MOIST_DEFAULTS.iceFall * Math.pow(density * grown / cover, MOIST_DEFAULTS.iceFallExponent), courant = speed * dt * sigmaMid[k] * g / (R * T * dSigma[k]), kept = grown / (1 + courant);
  const liquid = 1e-4 - grown, left = kept + liquid * Math.exp(-dt / MOIST_DEFAULTS.cloudLifetime);
  console.log(`${(p / 100).toFixed(0)} hPa, ${T.toFixed(1)} K, 0.1 g/kg of liquid over a cover of ${cover.toFixed(3)}: ${nuclei.toFixed(0)} nuclei/m³ (f_nuc ${active.toFixed(3)}) deposit ${(1e6 * grown).toFixed(3)} mg/kg of ice in ${dt} s, of which ${(1e6 * kept).toFixed(3)} stays; the layer below holds ${(1e6 * ice[(k + 1) * C]).toFixed(3)} mg/kg, all ice`);
  assert.ok(Math.abs(ice[k * C] - kept) < 1e-12 * kept, `ice ${ice[k * C]} against ${kept}`);
  assert.ok(Math.abs(s.qc[k * C] - left) < 1e-12 * left, `cloud ${s.qc[k * C]} against ${left}`);
  assert.ok(s.qc[(k + 1) * C] > 0 && ice[(k + 1) * C] === s.qc[(k + 1) * C], 'the clear layer below keeps the ice it takes as ice');
  assert.ok(Math.abs(totals(model).water + rain - before.water) < 1e-14 * before.water);
  const plain = run({ iceDeposition: true, bergeronAutoconversion: false });
  assert.deepEqual(Float64Array.from(plain.model.moist.cloudIceShare), Float64Array.from(ice));
  assert.deepEqual([...plain.s.qc, ...plain.s.q, ...plain.s.theta, plain.rain], [...s.qc, ...s.q, ...s.theta, rain]);
});

test('a layer is a cloud top when its cover is at least 0.01 under one with less: with 0.005 above it deposits at f_nuc + 0.1 (1 − f_nuc), with 0.5 above at Δz_top its own thickness; just inside 235.15 K a thin cloud glaciates whole within the step, its liquid the limit, and just inside 273.15 K hardly any', () => {
  const dt = 600, cold = 255, run = (above) => {
    const model = build({ iceDeposition: true }), { C, sigmaMid, exnerLayer } = model.core.diagnostics;
    const k = layerAt(model, 5.5e4), split = 1e5 * (sigmaMid[k - 1] + sigmaMid[k]) / 2;
    const s = column(model, (j, T, p) => 0.3 * saturated(T, p), (p) => (p < split ? 274 : cold));
    const Tabove = s.theta[(k - 1) * C] * exnerLayer[(k - 1) * C], pAbove = 1e5 * sigmaMid[k - 1], b = halfWidth(model, Tabove, pAbove);
    s.qc[(k - 1) * C] = above * above * b;
    s.qc[k * C] = 1e-4;
    model.core.diagnoseColumn(0, s.pi, s.theta, s.q, s.qc);
    model.moist.cloudIceShare.fill(0);
    const T = s.theta[k * C] * exnerLayer[k * C], p = 1e5 * sigmaMid[k];
    model.moist.autoconvertColumn(0, s.pi, s.theta, s.q, s.qc, dt);
    return { model, k, T, p, coverAbove: uniformCover(above * above * b, b), ice: model.moist.cloudIceShare[k * C] };
  };
  const thin = run(0.005), thick = run(0.5), { model, k, T, p } = thin, { C, sigmaMid, dSigma, g, R } = model.core.diagnostics;
  assert.ok(thin.coverAbove < 0.01 && thick.coverAbove >= 0.01, `covers above ${thin.coverAbove} and ${thick.coverAbove}`);
  const cover = uniformCover(1e-4, halfWidth(model, T, p)), density = p / (R * T), thickness = R * T * dSigma[k] / (g * sigmaMid[k]);
  const kept = (grown) => grown / (1 + MOIST_DEFAULTS.iceFall * Math.pow(density * grown / cover, MOIST_DEFAULTS.iceFallExponent) * dt * sigmaMid[k] * g / (R * T * dSigma[k]));
  const top = kept(depositionGrowth(T, p, 0, 1e-4, cover, dt, 0, R)), below = kept(depositionGrowth(T, p, 0, 1e-4, cover, dt, thickness, R));
  console.log(`${T.toFixed(1)} K at ${(p / 100).toFixed(0)} hPa under a warm layer of cover ${thin.coverAbove.toFixed(4)} or ${thick.coverAbove.toFixed(3)}: ${(1e6 * thin.ice).toFixed(4)} mg/kg of ice as a cloud top, ${(1e6 * thick.ice).toFixed(4)} below one ${thickness.toFixed(0)} m deep`);
  assert.ok(Math.abs(thin.ice - top) < 1e-12 * top && Math.abs(thick.ice - below) < 1e-12 * below && thick.ice > thin.ice, `${thin.ice} against ${top}; ${thick.ice} against ${below}`);
  for (const [Tlayer, water, check] of [[235.3, 2e-6, (ice, qc, made) => made === 2e-6 && ice === qc], [273.0, 1e-4, (ice, qc, made) => ice > 0 && ice <= made && made < 0.02 * 1e-4]]) {
    const model2 = build({ iceDeposition: true }), { C: C2, sigmaMid: sm, exnerLayer } = model2.core.diagnostics, k2 = layerAt(model2, 6e4);
    const s = column(model2, (j, Tj, pj) => 0.3 * saturated(Tj, pj), () => Tlayer);
    s.qc[k2 * C2] = water;
    model2.core.diagnoseColumn(0, s.pi, s.theta, s.q, s.qc);
    model2.moist.cloudIceShare.fill(0);
    const T2 = s.theta[k2 * C2] * exnerLayer[k2 * C2], p2 = 1e5 * sm[k2], a = uniformCover(water, halfWidth(model2, T2, p2));
    const made = depositionGrowth(T2, p2, 0, water, Math.max(CLEAR_AIR, a), dt, 0, model2.core.diagnostics.R);
    model2.moist.autoconvertColumn(0, s.pi, s.theta, s.q, s.qc, dt);
    const ice = model2.moist.cloudIceShare[k2 * C2], qc = s.qc[k2 * C2];
    console.log(`${Tlayer} K: ${(1e6 * water).toFixed(0)} mg/kg of liquid deposits ${(1e6 * made).toFixed(4)} mg/kg in ${dt} s; the layer ends with ${(1e6 * ice).toFixed(4)} of ice in ${(1e6 * qc).toFixed(4)} of cloud`);
    assert.ok(check(ice, qc, made), `${Tlayer} K: ice ${ice}, cloud ${qc}, deposited ${made}`);
  }
});

test('from its seed deposition grows 0.1 g/kg of in-cloud ice in 171, 68, 37, 23 and 15.6 min at 268, 263, 258, 253 and 248 K at 600 hPa and in 5.2 min at 240 K and 400 hPa; at 248 K the time grows about as p^1.5, and a cloud of 0.02 or 0.01 g/kg glaciates 2.9 or 4.6 times faster', () => {
  const { R } = build().core.diagnostics, step = 5;
  const minutes = (T, p, target) => { let ice = 0, t = 0; while (ice < target) { ice += depositionGrowth(T, p, ice, 1 - ice, 1, step, 1e4, R); t += step; } return t / 60; };
  const rows = [[268, 6e4, 171], [263, 6e4, 68], [258, 6e4, 37], [253, 6e4, 23], [248, 6e4, 15.6], [240, 4e4, 5.2], [248, 4e4, 8.7], [248, 8e4, 24.0]].map(([T, p, note]) => [T, p, note, minutes(T, p, 1e-4)]);
  console.log(rows.map(([T, p, note, m]) => `${T} K at ${p / 100} hPa: ${m.toFixed(1)} min (the note's ${note})`).join('; '));
  for (const [T, p, note, m] of rows) assert.ok(Math.abs(m - note) <= 0.03 * note + step / 60, `${T} K at ${p} Pa: ${m} min against ${note}`);
  const exponent = Math.log(rows[7][3] / rows[6][3]) / Math.log(2), thin = [2e-5, 1e-5].map((q) => rows[2][3] / minutes(258, 6e4, q));
  console.log(`p^${exponent.toFixed(2)} at 248 K; at 258 K 0.02 g/kg ${thin[0].toFixed(2)} and 0.01 g/kg ${thin[1].toFixed(2)} times faster than 0.1`);
  assert.ok(Math.abs(exponent - 1.5) < 0.1 && Math.abs(thin[0] - 2.9) < 0.15 && Math.abs(thin[1] - 4.6) < 0.2);
});

test('the memory is a fraction: with a layer\'s cloud doubled between steps, the vapor and heat that hold its equilibrium given back, the step starts from twice the ice and stores about twice the glaciated fraction', () => {
  const model = build({ iceDeposition: true }), { C, cp, exnerLayer } = model.core.diagnostics, [, theta, , , q, qc] = model.state;
  column(model, (j, T, p) => 1.05 * saturated(T, p));
  placeCloudIce(model, null);
  model.moist.adjust(model.state, 0, 1, 1350);
  const K = model.core.K;
  let k = -1;
  for (let j = 0; j < K; j++) { const T = theta[j * C] * exnerLayer[j * C]; if (k < 0 && T > 258 && T < 270 && qc[j * C] > 1e-5) k = j; }
  assert.ok(k >= 0, 'a mixed cloudy layer');
  const start = model.state.map((a) => Float64Array.from(a)), share = Float64Array.from(model.moist.cloudIceShare), x = k * C;
  share[x] = 0.3;
  const step = (doubled) => {
    const twin = build({ iceDeposition: true });
    twin.state.forEach((a, n) => a.set(start[n]));
    twin.core.diagnose(twin.state[0], twin.state[1], twin.state[4], twin.state[5]);
    twin.moist.cloudIceShare.set(share);
    if (doubled) { const added = twin.state[5][x]; twin.state[5][x] += added; twin.state[4][x] -= added; twin.state[1][x] += LATENT_HEAT * added / (cp * exnerLayer[x]); }
    twin.moist.adjust(twin.state, 0, 1, 1);
    return { qc: twin.state[5][x], share: twin.moist.cloudIceShare[x] };
  };
  const once = step(false), twice = step(true);
  console.log(`layer ${k}: ${(1e6 * once.qc).toFixed(3)} and ${(1e6 * twice.qc).toFixed(3)} mg/kg of cloud after the condensation; glaciated fraction ${once.share.toFixed(5)} and ${twice.share.toFixed(5)}`);
  assert.ok(Math.abs(twice.qc - once.qc) < 1e-3 * once.qc, 'the condensation restores the same cloud but for the extra ice that falls');
  assert.ok(Math.abs(twice.share / once.share - 2) < 0.01, `${twice.share / once.share}`);
});

test('under iceDeposition bergeronAutoconversion changes nothing: a step of every column is the same bit for bit', () => {
  const run = (bergeronAutoconversion) => {
    const model = build({ iceDeposition: true, bergeronAutoconversion, plumeEntrainmentLaw: 'gregory' });
    randomState(model);
    placeCloudIce(model, null);
    model.moist.precipitation.fill(0);
    model.phases.adjust(0, model.mesh.nCells, 1350);
    return [...model.state.flatMap((a) => [...a]), ...model.moist.cloudIceShare, ...model.moist.precipitation];
  };
  const on = run(true), off = run(false);
  assert.deepEqual(on, off);
});

test("the plumes' condensation pass caps the cloud ice the first pass left above its layer's cloud before it adds the ramp share: at 254.15 K (α = 0.5) a cloud all ice that the first pass halves and the second regrows by 0.3 of its start holds 0.8125 of its water as ice", () => {
  const model = build({ iceDeposition: true }), { moist } = model, { C, cp, sigmaMid, exnerLayer } = model.core.diagnostics;
  const { pi, theta, q, qc } = column(model, (k, T, p) => 0.3 * saturated(T, p), () => 254.15);
  const k = layerAt(model, 5e4), idx = k * C, p = pi[0] * sigmaMid[k], ice = moist.cloudIceShare;
  const vapor = (cloud, held) => {
    const liquidT = theta[idx] * exnerLayer[idx] - LATENT_HEAT * cloud / cp, out = cloudSaturation(liquidT, p, true, 273.15, 235.15, {}), b = halfWidth(model, liquidT, p);
    return (held >= b ? held : 2 * Math.sqrt(b * held) - b) * (1 + LATENT_HEAT * out.slope / cp) + out.qs - cloud;
  };
  qc[idx] = 1e-4; q[idx] = vapor(1e-4, 5e-5); ice.fill(0); ice[idx] = 1e-4;
  moist.condenseColumn(0, pi, theta, q, qc);
  assert.ok(Math.abs(qc[idx] - 5e-5) < 1e-9 && ice[idx] === 1e-4, `the first pass leaves ${qc[idx]} of cloud and ${ice[idx]} of ice`);
  theta[idx] = 254.15 / exnerLayer[idx];
  const before = qc[idx], T = theta[idx] * exnerLayer[idx];
  q[idx] = vapor(before, before + 3e-5);
  moist.condenseColumn(0, pi, theta, q, qc, true);
  const grown = qc[idx] - before, expected = Math.min(1e-4, before) + (1 - liquidFraction(T, 273.15, 235.15)) * grown;
  console.log(`cloud ${before.toExponential(4)} + ${grown.toExponential(4)} kg/kg; ice ${ice[idx].toExponential(4)}, ${(ice[idx] / qc[idx]).toFixed(6)} of the cloud`);
  assert.ok(Math.abs(grown - 3e-5) < 1e-9, `${grown}`);
  assert.ok(Math.abs(ice[idx] - expected) <= 1e-14 * expected, `${ice[idx]} against ${expected}`);
  assert.ok(Math.abs(ice[idx] / qc[idx] - 0.8125) < 1e-6, `${ice[idx] / qc[idx]}`);
});

test("placeCloudIce places a saved glaciated fraction of the model's length, and otherwise, as for a fresh start, the ramp's share 1 − α(T) where there is cloud and 0 where there is none; a model without the memory takes nothing", () => {
  const model = build({ iceDeposition: true }), { K, C } = model.core.diagnostics, [, theta, , , , qc] = model.state;
  randomState(model);
  const saved = Float64Array.from({ length: K * C }, (_, x) => (x % 7) / 7);
  placeCloudIce(model, { cloudIceShare: saved });
  assert.deepEqual(model.moist.cloudIceShare, saved);
  let ramp = null;
  for (const other of [null, {}, { cloudIceShare: saved.subarray(C) }]) {
    model.moist.cloudIceShare.fill(0.5);
    placeCloudIce(model, other);
    const exner = model.core.diagnostics.exnerLayer;
    ramp ??= Float64Array.from({ length: K * C }, (_, x) => (qc[x] > CLEAR_AIR ? 1 - liquidFraction(theta[x] * exner[x]) : 0));
    assert.deepEqual(model.moist.cloudIceShare, ramp);
  }
  assert.ok(ramp.some((v) => v > 0 && v < 1) && ramp.some((v) => v === 1) && ramp.some((v, x) => v === 0 && qc[x] > 0));
  const plain = build();
  placeCloudIce(plain, { cloudIceShare: saved });
  assert.equal(plain.moist.cloudIceShare, null);
});

test("after a plume the post-plume saturation pass caps the stale cloud ice at the cloud the first pass left, then adds the ramp's share of the regrowth under 'ramp' and nothing under 'liquid'", () => {
  const results = {};
  for (const detrainedPhase of ['ramp', 'liquid']) {
    const model = build({ iceDeposition: true, detrainedPhase }), { moist } = model, { K, C, cp, sigmaMid, exnerLayer } = model.core.diagnostics, [pi, theta, , , q, qc] = model.state;
    pi.fill(1e5);
    for (let pass = 0; pass < 3; pass++) {
      model.core.diagnose(pi, theta);
      for (let k = 0; k < K; k++) { const idx = k * C, p = pi[0] * sigmaMid[k]; theta[idx] = 254.15 / exnerLayer[idx]; q[idx] = 0.3 * saturated(254.15, p); qc[idx] = 0; }
    }
    model.core.diagnoseColumn(0, pi, theta, q, qc);
    const k = layerAt(model, 5e4), idx = k * C, p = pi[0] * sigmaMid[k], ice = moist.cloudIceShare;
    const vapor = (cloud, held) => {
      const liquidT = theta[idx] * exnerLayer[idx] - LATENT_HEAT * cloud / cp, out = cloudSaturation(liquidT, p, true, 273.15, 235.15, {}), b = halfWidth(model, liquidT, p);
      return (held >= b ? held : 2 * Math.sqrt(b * held) - b) * (1 + LATENT_HEAT * out.slope / cp) + out.qs - cloud;
    };
    qc[idx] = 1e-4; q[idx] = vapor(1e-4, 5e-5); ice.fill(0); ice[idx] = 1e-4;
    moist.condenseColumn(0, pi, theta, q, qc);
    const left = qc[idx];
    theta[idx] = 254.15 / exnerLayer[idx];
    q[idx] = vapor(left, left + 3e-5);
    moist.condenseColumn(0, pi, theta, q, qc, true);
    results[detrainedPhase] = { left, regrown: qc[idx], ice: ice[idx], alpha: liquidFraction(254.15, 273.15, 235.15) };
  }
  const { ramp, liquid } = results;
  console.log(`the first pass leaves ${ramp.left.toExponential(3)} kg/kg under 1e-4 of ice and the plume's pass regrows it to ${ramp.regrown.toExponential(3)}: ice ${ramp.ice.toExponential(4)} under 'ramp', ${liquid.ice.toExponential(4)} under 'liquid'`);
  assert.ok(ramp.left < 1e-4 && ramp.regrown > ramp.left, 'the first pass takes cloud below the stored ice and the plume regrows it');
  assert.ok(Math.abs(ramp.ice - (ramp.left + (1 - ramp.alpha) * (ramp.regrown - ramp.left))) < 1e-12 * ramp.ice, `ramp ${ramp.ice}`);
  assert.ok(Math.abs(liquid.ice - liquid.left) < 1e-12 * liquid.left, `liquid ${liquid.ice} against ${liquid.left}`);
});

test('trace.ice tallies, in kg/m² a layer, the deposition and the ice the caps remove at each visit and at the store, and changes no step', () => {
  const plain = build({ iceDeposition: true }), traced = build({ iceDeposition: true }), C = plain.mesh.nCells, K = plain.core.K;
  for (const model of [plain, traced]) { randomState(model); placeCloudIce(model, null); }
  traced.moist.trace.ice = { deposited: new Float64Array(K * C), visitCap: new Float64Array(K * C), storeCap: new Float64Array(K * C) };
  for (const model of [plain, traced]) { model.moist.precipitation.fill(0); model.phases.adjust(0, C, 1350); }
  plain.state.forEach((a, n) => assert.deepEqual(traced.state[n], a));
  assert.deepEqual(traced.moist.cloudIceShare, plain.moist.cloudIceShare);
  const t = traced.moist.trace.ice, sum = (a) => a.reduce((s, v) => s + v, 0);
  console.log(`deposited ${sum(t.deposited).toExponential(2)} kg/m², capped ${sum(t.visitCap).toExponential(2)} at the visits and ${sum(t.storeCap).toExponential(2)} at the store, summed over the layers of ${C} columns`);
  assert.ok(sum(t.deposited) > 0 && t.deposited.every((v) => v >= 0) && [t.visitCap, t.storeCap].every((a) => a.every((v) => v >= 0)));
  assert.ok(sum(t.visitCap) + sum(t.storeCap) > 0, 'the caps removed some ice');
});

test('trace.ice tallies the ice a visit caps at the cloud, in kg/m², by its hand value', () => {
  const model = build({ iceDeposition: true }), { C, dSigma, g } = model.core.diagnostics, k = layerAt(model, 4e4);
  const s = column(model, (j, T, p) => (j === k ? 1 : 0.3) * saturated(T, p));
  s.qc[k * C] = 1e-4;
  model.core.diagnoseColumn(0, s.pi, s.theta, s.q, s.qc);
  model.moist.cloudIceShare.fill(0);
  model.moist.cloudIceShare[k * C] = 1.6e-4;
  model.moist.trace.ice = { deposited: new Float64Array(model.core.K * C), visitCap: new Float64Array(model.core.K * C), storeCap: new Float64Array(model.core.K * C) };
  model.moist.autoconvertColumn(0, s.pi, s.theta, s.q, s.qc, 600);
  const mass = s.pi[0] * dSigma[k] / g, expected = (1.6e-4 - 1e-4) * mass, t = model.moist.trace.ice;
  assert.ok(Math.abs(t.visitCap[k * C] - expected) < 1e-12 * expected, `${t.visitCap[k * C]} kg/m² against ${expected}`);
  for (let x = 0; x < t.visitCap.length; x++) if (x !== k * C) assert.equal(t.visitCap[x], 0, `layer ${Math.floor(x / C)}`);
});

test('trace.ice tallies the ice the store caps at the cloud, in kg/m², by its hand value inside the mixed-phase range only, and stores the shares 1, 1 and 0 in the layers colder, inside and warmer than it', () => {
  const model = build({ iceDeposition: true }), { C, K, dSigma, g, exnerLayer } = model.core.diagnostics, k = layerAt(model, 4e4), cold = layerAt(model, 1e4), warm = layerAt(model, 9e4);
  const s = column(model, (j, T, p) => 0.3 * saturated(T, p));
  for (const at of [k, cold, warm]) s.qc[at * C] = 1e-4;
  model.core.diagnoseColumn(0, s.pi, s.theta, s.q, s.qc);
  const T = (at) => s.theta[at * C] * exnerLayer[at * C];
  assert.ok(T(k) > 235.15 && T(k) < 273.15, `the layer at ${T(k)} K is inside the mixed-phase range`);
  assert.ok(T(cold) < 235.15 && T(warm) > 273.15, `the other two layers are at ${T(cold)} and ${T(warm)} K`);
  const ice = model.moist.cloudIceShare;
  ice.fill(0);
  for (const at of [k, cold, warm]) ice[at * C] = 1.6e-4;
  model.moist.trace.ice = { deposited: new Float64Array(K * C), visitCap: new Float64Array(K * C), storeCap: new Float64Array(K * C) };
  model.moist.storeIce(0, s.pi, s.theta, s.qc);
  const expected = (1.6e-4 - 1e-4) * (s.pi[0] * dSigma[k] / g), t = model.moist.trace.ice;
  assert.ok(Math.abs(t.storeCap[k * C] - expected) < 1e-12 * expected, `${t.storeCap[k * C]} kg/m² against ${expected}`);
  for (let x = 0; x < K * C; x++) if (x !== k * C) assert.equal(t.storeCap[x], 0, `layer ${Math.floor(x / C)} of column ${x % C}`);
  assert.equal(ice[k * C], 1);
  assert.equal(ice[cold * C], 1);
  assert.equal(ice[warm * C], 0);
  assert.deepEqual(t.visitCap, new Float64Array(K * C));
});

test('every script that restores subcloudVirtual places the glaciated fraction, and every encodeState that writes subcloudVirtual writes it', () => {
  const files = [join(root, 'js/model.worker.js'), ...readdirSync(join(root, 'scripts'), { recursive: true }).filter((f) => f.endsWith('.mjs')).map((f) => join(root, 'scripts', f))];
  const missing = [];
  for (const file of files) {
    const text = readFileSync(file, 'utf8'), name = file.slice(root.length);
    if ((/subcloudVirtual\.set\(/.test(text) || text.includes("['subcloudVirtual'")) && !text.includes('placeCloudIce(')) missing.push(`${name} restores subcloudVirtual without placeCloudIce`);
    if (text.includes('placeCumulus(') && !text.includes('placeCloudIce(')) missing.push(`${name} places the cumulus without placeCloudIce`);
    for (let at = text.indexOf('encodeState('); at >= 0; at = text.indexOf('encodeState(', at + 1)) {
      let depth = 0, end = at + 'encodeState'.length;
      for (; end < text.length; end++) { if (text[end] === '(') depth++; else if (text[end] === ')' && --depth === 0) break; }
      const call = text.slice(at, end);
      if (call.includes('subcloudVirtual') && !call.includes('cloudIceShare')) missing.push(`${name} writes subcloudVirtual without cloudIceShare at ${text.slice(0, at).split('\n').length}`);
    }
  }
  assert.deepEqual(missing, []);
});

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }

test('with iceDeposition off the GPU kernels and the PH layout are those of the default; on, only the three adjust kernels change, CICE sits immediately before T2M, and every field from T2M on moves by K × C', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { createGpuCore } = await import('../js/gpu/core.gpu.js'), { mesh, core } = createModel(new Grid(16), { ocean: false, levels: bl48 });
  const plain = await createGpuCore(mesh, { levels: bl48, physics: {} }), off = await createGpuCore(mesh, { levels: bl48, physics: { iceDeposition: false } }), on = await createGpuCore(mesh, { levels: bl48, physics: { iceDeposition: true } });
  assert.deepEqual(off.kernelBodies, plain.kernelBodies);
  assert.deepEqual(off.layout.PH, plain.layout.PH);
  for (const [name, body] of Object.entries(off.kernelBodies)) for (const word of ['CICE', 'DEP_TOP_COVER', 'DEP_ACTIVE', 'dzTop', 'saturateDetrained', 'deposition(']) assert.ok(!body.includes(word), `${name} has ${word} with the option off`);
  assert.deepEqual(Object.keys(on.kernelBodies), Object.keys(off.kernelBodies));
  const changed = Object.keys(on.kernelBodies).filter((name) => on.kernelBodies[name] !== off.kernelBodies[name]);
  assert.deepEqual(changed.sort(), [...ADJUST_KERNELS].sort());
  const names = Object.keys(on.layout.PH).filter((name) => name !== 'total'), order = names.sort((a, b) => on.layout.PH[a] - on.layout.PH[b]);
  const at = order.indexOf('CICE'), length = core.diagnostics.K * mesh.nCells;
  assert.ok(at > 0);
  assert.equal(order[at + 1], 'T2M');
  assert.equal(on.layout.PH.T2M - on.layout.PH.CICE, length);
  for (const name of order.slice(0, at)) assert.equal(on.layout.PH[name], off.layout.PH[name], name);
  for (const name of [...order.slice(at + 1), 'total']) assert.equal(on.layout.PH[name], off.layout.PH[name] + length, name);
  assert.equal(off.layout.PH.CICE, undefined);
});

test("under iceDeposition the GPU's adjustPlume kernel caps the stored ice at the cloud after a plume, and adds the ramp's share of the regrowth only under detrainedPhase 'ramp'", { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { createGpuCore } = await import('../js/gpu/core.gpu.js'), { mesh } = createModel(new Grid(16), { ocean: false, levels: bl48 });
  const cap = 'PH[PH_CICE + idx] = min(max(0.0, PH[PH_CICE + idx]), IN[S_QC + idx])', share = ' + (1.0 - clamp((temperature - ICE_T) / (LIQUID_T - ICE_T), 0.0, 1.0)) * change';
  const [ramp, liquid] = [await createGpuCore(mesh, { levels: bl48, physics: { iceDeposition: true, detrainedPhase: 'ramp' } }), await createGpuCore(mesh, { levels: bl48, physics: { iceDeposition: true, detrainedPhase: 'liquid' } })];
  assert.ok(ramp.kernelBodies.adjustPlume.includes(`if (change > 0.0) { ${cap}${share}; }`), "'ramp' caps the stored ice and adds the ramp's share");
  assert.ok(liquid.kernelBodies.adjustPlume.includes(`if (change > 0.0) { ${cap}; }`), "'liquid' caps the stored ice and adds nothing");
  assert.ok(!liquid.kernelBodies.adjustPlume.includes(share), "'liquid' has no ramp share");
});

async function gpuParity(moist, { levels = bl34, dt = 1350, convective = false, reference = false, against = null, glaciated = false, exact = false } = {}) {
  const { createGpuCore } = await import('../js/gpu/core.gpu.js');
  const options = { plumeEntrainmentLaw: 'gregory', iceDeposition: true, ...moist };
  const model = createModel(new Grid(16), { ocean: false, levels, moist: options }), { core, mesh } = model, C = mesh.nCells, { K, exnerLayer } = core.diagnostics;
  const random = randomState(model, 4242, convective), [, theta, , , , qc] = model.state;
  placeCloudIce(model, null);
  const on = !!model.moist.cloudIceShare, share = model.moist.cloudIceShare ?? new Float64Array(K * C);
  if (on) for (let x = 0; x < K * C; x++) if (share[x] > 0 && share[x] < 1) share[x] = glaciated ? 1 : Math.fround(Math.min(1, Math.max(0, share[x] + 0.4 * (random() - 0.5))));
  const start = model.state.map((a) => Float64Array.from(a)), startShare = Float64Array.from(share), exner = Float64Array.from(exnerLayer);
  const gpu = await createGpuCore(mesh, { levels, physics: options });
  const { device, buffers, kernels } = gpu;
  gpu.upload(model.state);
  gpu.uploadPhysics({ mlmGate: model.radiation.mlmGate, concentration: model.seaIce.concentration, ...(on ? { cloudIceShare: share } : {}) });
  device.queue.writeBuffer(buffers.P, 0, Float32Array.from([dt, 0, 1, 0, 0, 0, 0, 0]));
  const group = device.createBindGroup({ layout: kernels.adjustMix.getBindGroupLayout(0), entries: [buffers.MI, buffers.MF, buffers.LV, buffers.S, buffers.K1, buffers.D, buffers.P, buffers.PH].map((buffer, binding) => ({ binding, resource: { buffer } })) });
  const encoder = device.createCommandEncoder(), computePass = encoder.beginComputePass();
  for (const name of ADJUST_KERNELS) { computePass.setPipeline(kernels[name]); computePass.setBindGroup(0, group); computePass.dispatchWorkgroups(Math.ceil(C / 64)); }
  computePass.end();
  device.queue.submit([encoder.finish()]);
  const after = await gpu.download(), ph = await gpu.downloadPhysics();
  model.moist.precipitation.fill(0);
  model.phases.adjust(0, C, dt);
  const gaps = { theta: new Float64Array(K * C), q: new Float64Array(K * C), qc: new Float64Array(K * C) };
  let worstTheta = 0, worstQ = 0, worstQc = 0, worstIce = 0, split = 0, mixed = 0, plumes = 0, worstQAt = 0;
  for (let i = 0; i < C; i++) if (model.moist.cumulusBaseFlux[i] > 0) plumes++;
  for (let x = 0; x < K * C; x++) {
    worstTheta = Math.max(worstTheta, Math.abs(theta[x] - after[1][x])); if (Math.abs(model.state[4][x] - after[4][x]) > worstQ) { worstQ = Math.abs(model.state[4][x] - after[4][x]); worstQAt = x; } worstQc = Math.max(worstQc, Math.abs(qc[x] - after[5][x]));
    gaps.theta[x] = theta[x] - after[1][x]; gaps.q[x] = model.state[4][x] - after[4][x]; gaps.qc[x] = qc[x] - after[5][x];
    const mine = share[x], theirs = on ? ph.CICE[x] : 0, T = theta[x] * exner[x];
    for (const [phi, cloud, who] of [[mine, qc[x], 'CPU'], [theirs, after[5][x], 'GPU']]) {
      assert.ok(phi >= 0 && phi <= 1, `${who} layer ${x}: ${phi}`);
      if (!(cloud > CLEAR_AIR)) assert.equal(phi, 0, `${who} layer ${x}: clear but ${phi}`);
    }
    if (on && qc[x] > CLEAR_AIR && T <= 235.15) assert.equal(mine, 1);
    if (T >= 273.15) assert.equal(mine, 0);
    if (mine > 0 && mine < 1) mixed++;
    const gap = Math.abs(mine * qc[x] - theirs * after[5][x]);
    if (gap > 1e-7) split++; else worstIce = Math.max(worstIce, gap);
  }
  console.log(`${JSON.stringify(moist)} at K=${K}, dt ${dt} s: ${plumes} of ${C} columns ran a plume, ${mixed} layers part ice; the engines' cloud ice differs by ${worstIce.toExponential(1)} kg/kg at most${split ? ` outside ${split} layers` : ''}; θ by ${worstTheta.toExponential(1)} K, q by ${worstQ.toExponential(1)}, q_c by ${worstQc.toExponential(1)}`);
  if (on) assert.ok(mixed > 0 && split <= (exact ? 0 : Math.max(1, 0.01 * mixed)), `${mixed} mixed layers; ${split} split`);
  if (against) {
    const moved = ['theta', 'q', 'qc'].map((name) => gaps[name].reduce((m, v, x) => Math.max(m, Math.abs(v - against.gaps[name][x])), 0));
    console.log(`beside the engines' gaps without the option, theirs with it move by θ ${moved[0].toExponential(1)} K, q ${moved[1].toExponential(1)}, q_c ${moved[2].toExponential(1)}`);
    assert.ok(moved[0] < 1e-3 && moved[1] < 1e-6 && moved[2] < 1e-7, `${moved}`);
  } else if (!reference) assert.ok(worstTheta < 1e-3 && worstQ < 1e-6 && worstQc < 1e-7, `θ ${worstTheta}, q ${worstQ}, qc ${worstQc} at layer ${Math.floor(worstQAt / C)} of column ${worstQAt % C}: ${start[4][worstQAt]} to ${model.state[4][worstQAt]} and ${after[4][worstQAt]}, T ${(start[1][worstQAt] * exner[worstQAt]).toFixed(1)} K`);
  return { start, startShare, share: Float64Array.from(share), plumes, gaps };
}

const GPU = { skip: !gpuAvailable && 'webgpu not installed' };
test('on the GPU iceDeposition matches the CPU column by column at N=16 on bl34 at 1350 s', GPU, () => gpuParity({}));
test('on the GPU iceDeposition matches the CPU column by column at N=16 on bl48, the page\'s levels, beside the engines\' gaps without it', GPU, async () => {
  const off = await gpuParity({ iceDeposition: false }, { levels: bl48, reference: true });
  await gpuParity({}, { levels: bl48, against: off });
});
test('on the GPU iceDeposition matches the CPU column by column at a 300 s step', GPU, () => gpuParity({}, { dt: 300 }));
test("on the GPU, with every mixed-phase cloud fully glaciated at the start, detrainedPhase 'liquid' matches the CPU's ice layer for layer in the plume columns", GPU, async () => {
  const run = await gpuParity({ detrainedPhase: 'liquid' }, { convective: true, glaciated: true, exact: true });
  assert.ok(run.plumes > 0, `${run.plumes} plumes`);
});
test("on the GPU iceDeposition with detrainedPhase 'liquid' matches the CPU, and the plumes' condensate is all liquid where 'ramp' gives it ice", GPU, async () => {
  const ramp = await gpuParity({}, { convective: true }), liquid = await gpuParity({ detrainedPhase: 'liquid' }, { convective: true });
  let more = 0;
  for (let x = 0; x < ramp.share.length; x++) if (ramp.share[x] > liquid.share[x] + 1e-12) more++;
  console.log(`${ramp.plumes} plumes; 'ramp' leaves more ice than 'liquid' in ${more} layers`);
  assert.ok(ramp.plumes > 0 && more > 0);
});
test('on the GPU iceDeposition with bergeronAutoconversion matches the CPU, which steps as without the factor', GPU, () => gpuParity({ bergeronAutoconversion: true }));
test("after the GPU's adjustMix and adjustPlume the cloud ice is the CPU's stage for stage where the plumes' condensation pass regrows cloud under ice the first pass left above it, and a layer set up as the hand example above a deep plume holds 0.8125 of its water as ice on both engines", GPU, async () => {
  const { createGpuCore } = await import('../js/gpu/core.gpu.js');
  const options = { plumeEntrainmentLaw: 'gregory', iceDeposition: true, cumulusMemory: 0 }, dt = 1350;
  const model = createModel(new Grid(16), { ocean: false, levels: bl48, moist: options }), { core, mesh, moist } = model, C = mesh.nCells, { K, cp, sigmaMid, exnerLayer } = core.diagnostics;
  randomState(model, 4242, true);
  placeCloudIce(model, null);
  const [pi, theta, u, , q, qc] = model.state, ice = moist.cloudIceShare;
  for (let x = 0; x < K * C; x++) if (ice[x] > 0) ice[x] = 1;
  const i0 = 0, k0 = layerAt(model, 5e3), idx0 = k0 * C + i0, p0 = pi[i0] * sigmaMid[k0], handTheta = Math.fround(254.15 / exnerLayer[idx0]), start = Math.fround(1e-4);
  const vapor = (cloud, held) => {
    const liquidT = handTheta * exnerLayer[idx0] - LATENT_HEAT * cloud / cp, out = cloudSaturation(liquidT, p0, true, 273.15, 235.15, {}), b = halfWidth(model, liquidT, p0);
    return (held >= b ? held : 2 * Math.sqrt(b * held) - b) * (1 + LATENT_HEAT * out.slope / cp) + out.qs - cloud;
  };
  theta[idx0] = handTheta; qc[idx0] = start; q[idx0] = Math.fround(vapor(start, 5e-5)); ice[idx0] = 1;

  const gpu = await createGpuCore(mesh, { levels: bl48, physics: options }), { device, buffers, kernels, layout } = gpu;
  gpu.upload(model.state);
  gpu.uploadPhysics({ mlmGate: model.radiation.mlmGate, concentration: model.seaIce.concentration, cloudIceShare: ice });
  device.queue.writeBuffer(buffers.P, 0, Float32Array.from([dt, 0, 1, 0, 0, 0, 0, 0]));
  const group = device.createBindGroup({ layout: kernels.adjustMix.getBindGroupLayout(0), entries: [buffers.MI, buffers.MF, buffers.LV, buffers.S, buffers.K1, buffers.D, buffers.P, buffers.PH].map((buffer, binding) => ({ binding, resource: { buffer } })) });
  const run = (name) => {
    const encoder = device.createCommandEncoder(), computePass = encoder.beginComputePass();
    computePass.setPipeline(kernels[name]); computePass.setBindGroup(0, group); computePass.dispatchWorkgroups(Math.ceil(C / 64)); computePass.end();
    device.queue.submit([encoder.finish()]);
  };
  run('adjustMix');
  const gpuMid = (await gpu.download())[5][idx0], gpuStale = (await gpu.downloadPhysics()).CICE[idx0];
  device.queue.writeBuffer(buffers.S, 4 * (layout.S.TH + idx0), Float32Array.of(handTheta));
  device.queue.writeBuffer(buffers.S, 4 * (layout.S.Q + idx0), Float32Array.of(vapor(gpuMid, gpuMid + 3e-5)));
  run('adjustPlume');
  const after = await gpu.download(), ph = await gpu.downloadPhysics();

  let cpuMid = 0, cpuStale = 0, exercised = 0, testable = 0, split = 0, worstOff = 0, worstIce = 0, worstTheta = 0, worstQ = 0, worstQc = 0;
  const cloudBefore = new Float64Array(K), iceBefore = new Float64Array(K), excess = new Float64Array(K * C), regrowth = new Float64Array(K * C);
  for (let i = 0; i < C; i++) {
    core.diagnoseColumn(i, pi, theta, q, qc);
    for (let k = 0; k < K; k++) ice[k * C + i] *= Math.max(0, qc[k * C + i]);
    moist.condenseColumn(i, pi, theta, q, qc);
    if (i === i0) { cpuMid = qc[idx0]; cpuStale = ice[idx0]; theta[idx0] = handTheta; q[idx0] = vapor(cpuMid, cpuMid + 3e-5); }
    const hand = i === i0 ? [theta[idx0], q[idx0], qc[idx0]] : null;
    moist.plumeColumn(i, pi, theta, q, qc, dt, u);
    if (hand) assert.deepEqual([theta[idx0], q[idx0], qc[idx0]], hand, 'the plume reaches the hand layer');
    for (let k = 0; k < K; k++) { cloudBefore[k] = qc[k * C + i]; iceBefore[k] = ice[k * C + i]; }
    if (moist.cumulusBaseFlux[i] > 0) moist.condenseColumn(i, pi, theta, q, qc, true);
    for (let k = 0; k < K; k++) {
      const x = k * C + i;
      if (!(qc[x] > cloudBefore[k] && iceBefore[k] > cloudBefore[k])) continue;
      excess[x] = iceBefore[k] - Math.max(0, cloudBefore[k]); regrowth[x] = qc[x] - cloudBefore[k]; exercised++;
    }
  }
  for (let x = 0; x < K * C; x++) {
    worstTheta = Math.max(worstTheta, Math.abs(theta[x] - after[1][x])); worstQ = Math.max(worstQ, Math.abs(q[x] - after[4][x])); worstQc = Math.max(worstQc, Math.abs(qc[x] - after[5][x]));
    if (!(excess[x] > 0)) continue;
    const gap = Math.abs(ice[x] - ph.CICE[x]);
    worstIce = Math.max(worstIce, gap);
    if (excess[x] > 1e-7 && regrowth[x] > 1e-7) { testable++; worstOff = Math.max(worstOff, gap / excess[x]); if (gap > 0.5 * excess[x]) split++; }
  }
  const shares = [ice[idx0] / qc[idx0], ph.CICE[idx0] / after[5][idx0]];
  console.log(`# ${exercised} layers regrow cloud under ice the first pass left above it, ${testable} by more than 1e-7 kg/kg under more than 1e-7 kg/kg of it, which both engines resolve; there the engines' ice differs by ${worstOff.toFixed(3)} of the excess at most, by more than half of it in ${split}; ice by ${worstIce.toExponential(1)} kg/kg at most in all ${exercised}; θ by ${worstTheta.toExponential(1)} K, q by ${worstQ.toExponential(1)}, q_c by ${worstQc.toExponential(1)}`);
  console.log(`# hand layer at ${(p0 / 100).toFixed(1)} hPa above column ${i0}'s plume: the first pass leaves ${cpuMid.toExponential(4)} (CPU) and ${gpuMid.toExponential(4)} (GPU) under ${start.toExponential(4)} of ice; ice ${shares[0].toFixed(6)} (CPU) and ${shares[1].toFixed(6)} (GPU) of the cloud`);
  assert.ok(moist.cumulusBaseFlux[i0] > 0 && ph.CUMF[i0] > 0, `plume base flux ${moist.cumulusBaseFlux[i0]} (CPU), ${ph.CUMF[i0]} (GPU)`);
  assert.ok(cpuStale === start && gpuStale === start, `stale ice ${cpuStale} (CPU), ${gpuStale} (GPU)`);
  assert.ok(Math.abs(cpuMid - 5e-5) < 1e-8 && Math.abs(gpuMid - 5e-5) < 1e-7, `${cpuMid}, ${gpuMid}`);
  assert.ok(Math.abs(shares[0] - 0.8125) < 1e-6 && Math.abs(shares[1] - 0.8125) < 1e-3, `${shares}`);
  assert.ok(testable > 100 && split <= Math.max(2, 0.01 * testable), `${exercised} regrown, ${testable} testable, ${split} split, worst ${worstOff} of the excess`);
  assert.ok(worstTheta < 1e-3 && worstQ < 1e-6 && worstQc < 1e-7, `θ ${worstTheta}, q ${worstQ}, qc ${worstQc}`);
});
const REFUSED_RAIN = [{ rainEvaporationLaw: 'kessler', rainCover: 'overlap' }, { rainEvaporationCeiling: 0.7 }];
test("iceDeposition refuses a rainEvaporationLaw other than 'deficit' and a rainEvaporationCeiling on the CPU", () => {
  for (const rain of REFUSED_RAIN) assert.throws(() => build({ iceDeposition: true, ...rain }), /rainEvaporationLaw 'deficit'/, JSON.stringify(rain));
});
test("iceDeposition refuses a rainEvaporationLaw other than 'deficit' and a rainEvaporationCeiling on the GPU", GPU, async () => {
  const { createGpuCore } = await import('../js/gpu/core.gpu.js'), { mesh } = createModel(new Grid(16), { ocean: false, levels: bl34 });
  for (const rain of REFUSED_RAIN) await assert.rejects(createGpuCore(mesh, { levels: bl34, physics: { iceDeposition: true, ...rain } }), /rainEvaporationLaw 'deficit'/, JSON.stringify(rain));
});
test('a GPU core refuses a glaciated fraction whose length is not K × C, which would run into the next field', GPU, async () => {
  const { createGpuCore } = await import('../js/gpu/core.gpu.js'), model = build({ iceDeposition: true }), { K } = model.core.diagnostics, C = model.mesh.nCells;
  const gpu = await createGpuCore(model.mesh, { levels: bl34, physics: { iceDeposition: true } }), physics = (cloudIceShare) => ({ mlmGate: model.radiation.mlmGate, concentration: model.seaIce.concentration, cloudIceShare });
  for (const length of [K * C - 1, K * C + C]) assert.throws(() => gpu.uploadPhysics(physics(new Float32Array(length))), /cloudIceShare must hold K × C/, `${length}`);
  assert.doesNotThrow(() => gpu.uploadPhysics(physics(new Float32Array(K * C))));
});
test('a GPU model under iceDeposition takes the ramp share 1 − α(T), in f32, for each cloudy layer from a fresh start, a save without the glaciated fraction and one of another length, and steps from it', GPU, async () => {
  const { createGpuModel } = await import('../js/gpu/model.gpu.js');
  const model = await createGpuModel(new Grid(16), { ocean: false, moist: { iceDeposition: true } }), { K, C, exnerLayer } = model.core.diagnostics;
  randomState(model);
  if (model.land) model.land.initialize();
  const start = model.state.map((a) => Float64Array.from(a));
  for (const [label, saved] of [['a fresh start', null], ['a save without it', { N: 16 }], ['a save of another length', { cloudIceShare: new Float32Array(K * C - C).fill(0.5) }]]) {
    start.forEach((values, a) => model.state[a].set(values));
    model.moist.cloudIceShare.fill(0.37);
    placeCloudIce(model, saved);
    const [, theta, , , , qc] = model.state;
    const ramp = Float32Array.from({ length: K * C }, (_, x) => (qc[x] > CLEAR_AIR ? 1 - liquidFraction(theta[x] * exnerLayer[x], 273.15, 235.15) : 0));
    model.load();
    const placed = (await model.gpu.downloadPhysics()).CICE;
    let differ = 0, mixed = 0;
    for (let x = 0; x < K * C; x++) { if (placed[x] !== ramp[x]) differ++; if (ramp[x] > 0 && ramp[x] < 1) mixed++; }
    await model.step(1350);
    const after = await model.gpu.download(), stepped = (await model.gpu.downloadPhysics()).CICE.subarray(0, K * C);
    const finite = [1, 4, 5].every((a) => after[a].every(Number.isFinite)), inRange = stepped.every((v) => v >= 0 && v <= 1);
    console.log(`${label}: ${mixed} layers part ice of ${K * C}; ${differ} differ from the ramp; the step leaves θ, q and q_c ${finite ? 'finite' : 'not finite'} and the fraction ${inRange ? 'in [0, 1]' : 'out of range'}`);
    assert.ok(mixed > 0 && differ === 0 && finite && inRange, label);
  }
  model.destroy();
});
