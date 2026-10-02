import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { R_VAPOR, saturationHumidity, cloudSaturation, iceVaporPressure, saturationVaporPressure, criticalHumidityAt, uniformCondensate, uniformCover, LATENT_HEAT, MOIST_DEFAULTS } from '../js/physics/moist.module.js';

const levels = sigmaInterfaces('bl34');
const build = (moist = {}) => createModel(new Grid(2), { ocean: false, levels, moist });

/*
 * Column 0 of the model at 1000 hPa under a 6.5 K/km atmosphere from
 * 300 K to a 210 K floor, humidity(k, T, p) of the layer's vapour.
 */
function column(model, humidity) {
  const { K, C, sigmaMid, exnerLayer } = model.core.diagnostics;
  const [pi, theta, , , q, qc] = model.state;
  pi.fill(1e5);
  for (let pass = 0; pass < 3; pass++) {
    model.core.diagnose(pi, theta);
    for (let k = 0; k < K; k++) {
      const idx = k * C, p = pi[0] * sigmaMid[k], T = Math.max(210, 300 - 6.5e-3 * 7400 * Math.log(1e5 / p));
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

test('saturation over ice follows the phase ramp: liquid above 273.15 K, ice below 235.15 K, the mix between, and the ice form 0.68 of water at -40 °C', () => {
  const out = { qs: 0, slope: 0, liquid: 0 };
  assert.equal(cloudSaturation(280, 8e4, true, 273.15, 235.15, out).qs, saturationHumidity(280, 8e4));
  assert.equal(out.slope, out.qs * LATENT_HEAT / (R_VAPOR * 280 * 280));
  const cold = cloudSaturation(220, 25e3, true, 273.15, 235.15, out).qs, iceE = iceVaporPressure(220);
  assert.ok(Math.abs(cold - 0.622 * iceE / (25e3 - 0.378 * iceE)) < 1e-15);
  const mid = 254.15, alpha = (mid - 235.15) / 38, e = alpha * saturationVaporPressure(mid) + (1 - alpha) * iceVaporPressure(mid);
  assert.ok(Math.abs(cloudSaturation(mid, 5e4, true, 273.15, 235.15, out).qs - 0.622 * e / (5e4 - 0.378 * e)) < 1e-15);
  const ratio = iceVaporPressure(233.15) / saturationVaporPressure(233.15);
  console.log(`e_i/e_w at -40 °C ${ratio.toFixed(3)} (0.68 by Murphy and Koop 2005), at -20 °C ${(iceVaporPressure(253.15) / saturationVaporPressure(253.15)).toFixed(3)} (0.82)`);
  assert.ok(ratio > 0.66 && ratio < 0.70);
  assert.equal(cloudSaturation(220, 25e3, false, 273.15, 235.15, out).qs, saturationHumidity(220, 25e3));
});

test('the uniform distribution condenses from the critical humidity up, its cover the square root of its condensate over the half-width, and overcast with the grid saturated', () => {
  assert.equal(criticalHumidityAt(1e5, 1e5, 0.975, 0.75, 2), 0.975);
  assert.ok(Math.abs(criticalHumidityAt(5e4, 1e5, 0.975, 0.75, 2) - (0.75 + 0.225 * Math.exp(-3))) < 1e-15);
  const b = 2e-4;
  for (const Q of [-3e-4, -2e-4, -1e-4, 0, 1e-4, 2e-4, 5e-4]) {
    const qc = uniformCondensate(Q, b), f = uniformCover(qc, b);
    assert.ok(Math.abs(f - Math.min(1, Math.max(0, (Q + b) / (2 * b)))) < 1e-12, `Q ${Q}: cover ${f}`);
  }
  assert.equal(uniformCondensate(5e-4, b), 5e-4);
  assert.equal(uniformCondensate(-3e-4, b), 0);
  const model = build(), { moist } = model, { K, C } = model.core.diagnostics;
  const k = layerAt(model, 6e4);
  const { pi, theta, q, qc } = column(model, (j, T, p) => (j === k ? 0.95 : 0.5) * saturationHumidity(T, p));
  const before = totals(model);
  moist.condenseColumn(0, pi, theta, q, qc);
  const after = totals(model);
  const T = theta[k * C] * model.core.diagnostics.exnerLayer[k * C], rhc = criticalHumidityAt(pi[0] * model.core.diagnostics.sigmaMid[k], pi[0], 0.975, 0.75, 2);
  const warmCloudy = Array.from({ length: K }, (_, j) => j).filter((j) => j !== k && model.core.diagnostics.sigmaMid[j] > 0.5 && qc[j * C] > 0).length;
  console.log(`a layer at 600 hPa (${T.toFixed(1)} K, RH_c ${rhc.toFixed(3)}) at 0.95 of water saturation holds ${(1e3 * qc[k * C]).toFixed(4)} g/kg of cloud; of the layers below 500 hPa at 0.5, ${warmCloudy} hold any`);
  assert.ok(qc[k * C] > 0 && warmCloudy === 0);
  assert.ok(Math.abs(after.water - before.water) < 1e-14 * before.water && Math.abs(after.enthalpy - before.enthalpy) < 1e-14 * before.enthalpy);
  const plain = build({ condensation: 'saturation', iceSaturation: false });
  const s = column(plain, (j, T, p) => (j === k ? 0.95 : 0.5) * saturationHumidity(T, p));
  plain.moist.condenseColumn(0, s.pi, s.theta, s.q, s.qc);
  assert.equal(s.qc[k * C], 0);
});

test('air at 220 K saturated over ice holds cloud under iceSaturation and none over water alone', () => {
  const run = (options) => {
    const model = build({ condensation: 'saturation', ...options }), k = layerAt(model, 25e3), { C } = model.core.diagnostics;
    const s = column(model, (j, T, p) => (j === k ? 1.1 : 0.3) * cloudSaturation(T, p, true, 273.15, 235.15, {}).qs);
    model.moist.condenseColumn(0, s.pi, s.theta, s.q, s.qc);
    return s.qc[k * C];
  };
  const iced = run({ iceSaturation: true }), liquid = run({ iceSaturation: false });
  console.log(`at RH_ice 1.1 near 250 hPa: ${(1e6 * iced).toFixed(2)} mg/kg of cloud over ice, ${liquid} over water`);
  assert.ok(iced > 0 && liquid === 0);
});

test('with iceNucleation clear air colder than 235 K stays supersaturated over ice below the homogeneous nucleation threshold, while the same layer with cloud in it deposits to the distribution about ice saturation', () => {
  const run = (options, seed, humidity = 1.1) => {
    const model = build(options), k = layerAt(model, 2e4), { C } = model.core.diagnostics;
    const s = column(model, (j, T, p) => (j === k ? humidity : 0.3) * cloudSaturation(T, p, true, 273.15, 235.15, {}).qs);
    s.qc[k * C] = seed;
    model.core.diagnoseColumn(0, s.pi, s.theta, s.q, s.qc);
    model.moist.condenseColumn(0, s.pi, s.theta, s.q, s.qc);
    return s.qc[k * C];
  };
  const on = { iceNucleation: true }, clear = run(on, 0), seeded = run(on, 1e-6), free = run({}, 0), freeSeeded = run({}, 1e-6), nucleated = run(on, 0, 1.6);
  console.log(`near 200 hPa at RH_ice 1.1: ${clear} kg/kg of cloud in clear air (${(1e6 * free).toFixed(2)} mg/kg without the threshold), ${(1e6 * seeded).toFixed(2)} mg/kg where 1 mg/kg was; at RH_ice 1.6 clear air nucleates ${(1e6 * nucleated).toFixed(2)} mg/kg`);
  assert.ok(clear === 0 && free > 0 && seeded === freeSeeded && nucleated > 0);
});

test('cloud ice falls in the step: the layer keeps 1/(1 + v dt/dz) of it at the speed of Heymsfield and Donner\'s form, the layers below take the rest and sublimate it into dry air, column water with the snow and moist enthalpy exact; liquid cloud and iceFall null convert as before', () => {
  const dt = 600, model = build(), { moist, core } = model, { K, C, sigmaMid, dSigma, g, R, exnerLayer } = core.diagnostics;
  const k = layerAt(model, 2e4);
  const { pi, theta, q, qc } = column(model, (j, T, p) => (j === k ? 1 : 0.3) * cloudSaturation(T, p, true, 273.15, 235.15, {}).qs);
  qc[k * C] = 3e-5;
  core.diagnoseColumn(0, pi, theta, q, qc);
  const before = totals(model), T = theta[k * C] * exnerLayer[k * C], p = pi[0] * sigmaMid[k];
  assert.ok(T < 235.15, `the cloud layer at ${T} K is all ice`);
  const out = cloudSaturation(T, p, true, 273.15, 235.15, {});
  const rhc = criticalHumidityAt(p, pi[0], 0.975, 0.75, 2), b = (1 - rhc) * out.qs / (1 + LATENT_HEAT * out.slope / core.diagnostics.cp);
  const f = uniformCover(3e-5, b), speed = MOIST_DEFAULTS.iceFall * Math.pow(p / (R * T) * 3e-5 / f, MOIST_DEFAULTS.iceFallExponent), courant = speed * dt * sigmaMid[k] * g / (R * T * dSigma[k]);
  const qBelow = q[(k + 1) * C], thetaBelow = theta[(k + 1) * C];
  const snow = moist.autoconvertColumn(0, pi, theta, q, qc, dt);
  assert.ok(Math.abs(qc[k * C] - 3e-5 / (1 + courant)) < 1e-18, `kept ${qc[k * C]} against ${3e-5 / (1 + courant)}`);
  console.log(`ice at ${(p / 100).toFixed(0)} hPa, ${T.toFixed(1)} K, 0.030 g/kg over a cover of ${f.toFixed(3)}: falls at ${speed.toFixed(3)} m/s, keeps ${(1 / (1 + courant)).toFixed(4)} of itself over ${dt} s; ${(1e3 * snow).toExponential(2)} g/m² leaves the column`);
  assert.ok(qc[(k + 1) * C] > 0, 'the layer below takes the ice');
  const mid = totals(model);
  assert.ok(Math.abs(mid.water + snow - before.water) < 1e-14 * before.water);
  assert.ok(moist.falling.moved);
  moist.condenseColumn(0, pi, theta, q, qc);
  const after = totals(model);
  assert.ok(q[(k + 1) * C] > qBelow && theta[(k + 1) * C] < thetaBelow, 'the dry layer below gains vapour and cools');
  assert.ok(Math.abs(after.water + snow - before.water) < 1e-14 * before.water && Math.abs(after.enthalpy - before.enthalpy) < 1e-14 * before.enthalpy);
  const warm = layerAt(model, 8e4);
  const runWarm = (options) => {
    const m = build(options), s = column(m, (j, Tj, pj) => (j === warm ? 1 : 0.95) * saturationHumidity(Tj, pj));
    s.qc[warm * C] = 5e-4;
    m.core.diagnoseColumn(0, s.pi, s.theta, s.q, s.qc);
    const rain = m.moist.autoconvertColumn(0, s.pi, s.theta, s.q, s.qc, dt);
    return { rain, qc: s.qc[warm * C], moved: m.moist.falling.moved };
  };
  const fall = runWarm({}), none = runWarm({ iceFall: null });
  assert.equal(fall.qc, none.qc);
  assert.equal(fall.rain, none.rain);
  assert.equal(fall.moved, false);
  const old = build({ iceFall: null }), o = column(old, (j, Tj, pj) => (j === k ? 1 : 0.3) * cloudSaturation(Tj, pj, true, 273.15, 235.15, {}).qs);
  o.qc[k * C] = 3e-5;
  old.core.diagnoseColumn(0, o.pi, o.theta, o.q, o.qc);
  old.moist.autoconvertColumn(0, o.pi, o.theta, o.q, o.qc, dt);
  assert.ok(Math.abs(o.qc[k * C] - 3e-5 * Math.exp(-dt / MOIST_DEFAULTS.cloudLifetime)) < 1e-18 && o.qc[(k + 1) * C] === 0);
});

test('a deep column of falling ice is stable at any step: the implicit fall never empties a layer it feeds nor leaves negative cloud', () => {
  const model = build(), { moist, core } = model, { K, C } = core.diagnostics;
  const { pi, theta, q, qc } = column(model, (j, T, p) => cloudSaturation(T, p, true, 273.15, 235.15, {}).qs);
  for (let k = 0; k < K; k++) if (theta[k * C] * core.diagnostics.exnerLayer[k * C] < 240) qc[k * C] = 1e-3;
  core.diagnoseColumn(0, pi, theta, q, qc);
  const before = totals(model);
  let snow = 0;
  for (const dt of [60, 600, 6000, 60000]) snow += moist.autoconvertColumn(0, pi, theta, q, qc, dt);
  const after = totals(model);
  for (let k = 0; k < K; k++) assert.ok(qc[k * C] >= 0 && Number.isFinite(qc[k * C]));
  assert.ok(Math.abs(after.water + snow - before.water) < 1e-13 * before.water);
});

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }

test('on the GPU the falling ice and the uniform condensation keep column water with the precipitation and moist enthalpy to single precision, and match the CPU column by column', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { createGpuCore } = await import('../js/gpu/core.gpu.js');
  const dt = 900, model = createModel(new Grid(6), { ocean: false, levels, moist: {} }), { core, mesh } = model, C = mesh.nCells, { K, sigmaMid, exnerLayer, dSigma, g, cp } = core.diagnostics;
  const [pi, theta, , surfaceT, q, qc] = model.state;
  let seed = 4242;
  const random = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; };
  for (let pass = 0; pass < 3; pass++) {
    core.diagnose(pi.fill(1e5), theta, q, qc);
    for (let i = 0; i < C; i++) for (let k = 0; k < K; k++) {
      const x = k * C + i, p = 1e5 * sigmaMid[k], T = Math.max(205, 300 - 6.5e-3 * 7400 * Math.log(1e5 / p));
      theta[x] = T / exnerLayer[x];
      if (pass === 0) { q[x] = (0.3 + 0.75 * random()) * cloudSaturation(T, p, true, 273.15, 235.15, {}).qs; qc[x] = random() < 0.3 ? 2e-4 * random() * (T < 250 ? 0.2 : 1) : 0; }
    }
  }
  surfaceT.fill(300);
  for (const a of model.state) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  core.diagnose(pi, theta, q, qc);
  const exner = Float64Array.from(exnerLayer);
  const budget = (th, qq, cc, i) => { let h = 0, w = 0; for (let k = 0; k < K; k++) { const x = k * C + i, m = pi[i] * dSigma[k] / g; h += (cp * th[x] * exner[x] + LATENT_HEAT * qq[x]) * m; w += (qq[x] + cc[x]) * m; } return { h, w }; };
  const before = Array.from({ length: C }, (_, i) => budget(theta, q, qc, i));
  let icy = 0;
  for (let i = 0; i < C; i++) { let cold = false; for (let k = 0; k < K; k++) if (qc[k * C + i] > 0 && theta[k * C + i] * exner[k * C + i] < 235.15) cold = true; if (cold) icy++; }
  const gpu = await createGpuCore(mesh, { levels, physics: {} });
  const { device, buffers, kernels } = gpu;
  gpu.upload(model.state);
  gpu.uploadPhysics({ mlmGate: model.radiation.mlmGate, concentration: model.seaIce.concentration });
  device.queue.writeBuffer(buffers.P, 0, Float32Array.from([dt, 0, 1, 0, 0, 0, 0, 0]));
  const group = device.createBindGroup({ layout: kernels.adjust.getBindGroupLayout(0), entries: [buffers.MI, buffers.MF, buffers.LV, buffers.S, buffers.K1, buffers.D, buffers.P, buffers.PH].map((buffer, binding) => ({ binding, resource: { buffer } })) });
  const encoder = device.createCommandEncoder(), computePass = encoder.beginComputePass();
  computePass.setPipeline(kernels.adjust); computePass.setBindGroup(0, group); computePass.dispatchWorkgroups(Math.ceil(C / 64)); computePass.end();
  device.queue.submit([encoder.finish()]);
  const after = await gpu.download(), ph = await gpu.downloadPhysics();
  model.moist.precipitation.fill(0);
  model.phases.adjust(0, C, dt);
  let worstH = 0, worstW = 0, cpuH = 0, cpuW = 0, worstTheta = 0, worstQ = 0, worstQc = 0, snowing = 0;
  for (let i = 0; i < C; i++) {
    const g1 = budget(after[1], after[4], after[5], i), c1 = budget(theta, q, qc, i);
    worstH = Math.max(worstH, Math.abs(g1.h - before[i].h) / before[i].h); worstW = Math.max(worstW, Math.abs(g1.w + ph.STEPRAIN[i] - before[i].w) / before[i].w);
    cpuH = Math.max(cpuH, Math.abs(c1.h - before[i].h) / before[i].h); cpuW = Math.max(cpuW, Math.abs(c1.w + model.moist.precipitation[i] - before[i].w) / before[i].w);
    if (model.moist.precipitation[i] > 0) snowing++;
    for (let k = 0; k < K; k++) { const x = k * C + i; worstTheta = Math.max(worstTheta, Math.abs(theta[x] - after[1][x])); worstQ = Math.max(worstQ, Math.abs(q[x] - after[4][x])); worstQc = Math.max(worstQc, Math.abs(qc[x] - after[5][x])); }
  }
  console.log(`${C} columns, ${icy} with ice cloud, ${snowing} precipitating: GPU column moist enthalpy to ${worstH.toExponential(1)} and water with the precipitation to ${worstW.toExponential(1)} relative (CPU ${cpuH.toExponential(1)} and ${cpuW.toExponential(1)}); engines differ in θ by ${worstTheta.toExponential(1)} K, q by ${worstQ.toExponential(1)}, qc by ${worstQc.toExponential(1)}`);
  assert.ok(icy > C / 2);
  assert.ok(cpuH < 1e-14 && cpuW < 1e-14, `CPU ${cpuH}, ${cpuW}`);
  assert.ok(worstH < 1e-6 && worstW < 2e-6, `GPU ${worstH}, ${worstW}`);
  assert.ok(worstTheta < 1e-3 && worstQ < 1e-6 && worstQc < 1e-7, `θ ${worstTheta}, q ${worstQ}, qc ${worstQc}`);
});
