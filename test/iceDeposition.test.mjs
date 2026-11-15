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

test('every script that restores subcloudVirtual places the glaciated fraction, and every encodeState that writes subcloudVirtual writes it', () => {
  const files = [join(root, 'js/model.worker.js'), ...readdirSync(join(root, 'scripts'), { recursive: true }).filter((f) => f.endsWith('.mjs')).map((f) => join(root, 'scripts', f))];
  const missing = [];
  for (const file of files) {
    const text = readFileSync(file, 'utf8'), name = file.slice(root.length);
    if ((/subcloudVirtual\.set\(/.test(text) || text.includes("['subcloudVirtual'")) && !text.includes('placeCloudIce(')) missing.push(`${name} restores subcloudVirtual without placeCloudIce`);
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

test('with iceDeposition off the GPU kernels and the PH layout are those of the default; on, only the three adjust kernels change and CICE sits before T2M, which stays last', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { createGpuCore } = await import('../js/gpu/core.gpu.js'), { mesh } = createModel(new Grid(16), { ocean: false, levels: bl48 });
  const plain = await createGpuCore(mesh, { levels: bl48, physics: {} }), off = await createGpuCore(mesh, { levels: bl48, physics: { iceDeposition: false } }), on = await createGpuCore(mesh, { levels: bl48, physics: { iceDeposition: true } });
  assert.deepEqual(off.kernelBodies, plain.kernelBodies);
  assert.deepEqual(off.layout.PH, plain.layout.PH);
  for (const [name, body] of Object.entries(off.kernelBodies)) for (const word of ['CICE', 'DEP_TOP_COVER', 'DEP_ACTIVE', 'dzTop', 'saturateDetrained', 'deposition(']) assert.ok(!body.includes(word), `${name} has ${word} with the option off`);
  assert.deepEqual(Object.keys(on.kernelBodies), Object.keys(off.kernelBodies));
  const changed = Object.keys(on.kernelBodies).filter((name) => on.kernelBodies[name] !== off.kernelBodies[name]);
  assert.deepEqual(changed.sort(), [...ADJUST_KERNELS].sort());
  const names = Object.keys(on.layout.PH).filter((name) => name !== 'total'), order = names.sort((a, b) => on.layout.PH[a] - on.layout.PH[b]);
  assert.equal(order.at(-1), 'T2M');
  assert.equal(order.at(-2), 'CICE');
  for (const name of order.slice(0, -2)) assert.equal(on.layout.PH[name], off.layout.PH[name], name);
});

async function gpuParity(moist, { levels = bl34, dt = 1350, convective = false, reference = false, against = null } = {}) {
  const { createGpuCore } = await import('../js/gpu/core.gpu.js');
  const options = { plumeEntrainmentLaw: 'gregory', iceDeposition: true, ...moist };
  const model = createModel(new Grid(16), { ocean: false, levels, moist: options }), { core, mesh } = model, C = mesh.nCells, { K, exnerLayer } = core.diagnostics;
  const random = randomState(model, 4242, convective), [, theta, , , , qc] = model.state;
  placeCloudIce(model, null);
  const on = !!model.moist.cloudIceShare, share = model.moist.cloudIceShare ?? new Float64Array(K * C);
  if (on) for (let x = 0; x < K * C; x++) if (share[x] > 0 && share[x] < 1) share[x] = Math.fround(Math.min(1, Math.max(0, share[x] + 0.4 * (random() - 0.5))));
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
  if (on) assert.ok(mixed > 0 && split <= Math.max(1, 0.01 * mixed), `${mixed} mixed layers; ${split} split`);
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
test("on the GPU iceDeposition with detrainedPhase 'liquid' matches the CPU, and the plumes' condensate is all liquid where 'ramp' gives it ice", GPU, async () => {
  const ramp = await gpuParity({}, { convective: true }), liquid = await gpuParity({ detrainedPhase: 'liquid' }, { convective: true });
  let more = 0;
  for (let x = 0; x < ramp.share.length; x++) if (ramp.share[x] > liquid.share[x] + 1e-12) more++;
  console.log(`${ramp.plumes} plumes; 'ramp' leaves more ice than 'liquid' in ${more} layers`);
  assert.ok(ramp.plumes > 0 && more > 0);
});
test('on the GPU iceDeposition with bergeronAutoconversion matches the CPU, which steps as without the factor', GPU, () => gpuParity({ bergeronAutoconversion: true }));
