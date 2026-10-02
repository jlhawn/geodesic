import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { saturationHumidity } from '../js/physics/moist.module.js';
import { sunDirection, SOLAR_CONSTANT, STEFAN_BOLTZMANN, GREENHOUSE_GASES } from '../js/physics/radiation.module.js';
import { LONGWAVE_TABLE, LONGWAVE_CONSTANTS, GAS_MOLAR, layerPaths, planckShare } from '../js/physics/longwave.module.js';
import { OZONE_CM_ATM } from '../js/physics/shortwaveGases.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuCore } = gpuAvailable ? await import('../js/gpu/core.gpu.js') : {};

const SPECTRAL = { longwaveScheme: 'correlated', solarGases: 'clirad' };
const GRAY_GASES = { longwaveScheme: 'gray', solarGases: 'lacisHansen' };
const TRANSPARENT = { carbonDioxide: 0, methane: 0, nitrousOxide: 0, ozoneColumn: [0, 0] };
const PLAIN = { cloudCover: 'overcast', cloudSolarAbsorption: 0, rayleighDepth: 0, landAerosol: 0, seaAerosol: 0, skylight: 0, upwardAbsorption: false, clearSkyPass: true };
const HAND = { liquid: { reflectance: 0.6430392965333067 }, thinLiquid: { emissivity: 0.7768681886822757 }, ice: { reflectance: 0.418058949528354, emissivity: 0.8998461956940361 } };
const CLOUDS = { liquid: [285, 0.1], thinLiquid: [285, 0.01], ice: [213.15, 0.02] };
const close = (a, b, tolerance) => Math.abs(a - b) <= tolerance * Math.abs(b);

function cloudColumn(radiation, { T, path, humid, sigma }) {
  const model = createModel(new Grid(4), { ocean: false, radiation: { ...PLAIN, ...radiation } });
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  const { core, radiation: rad } = model, { K, sigmaMid } = core.diagnostics, C = model.mesh.nCells, [pi, theta] = model.state;
  const i = 0, q = new Float64Array(K * C), none = new Float64Array(K * C), qc = new Float64Array(K * C);
  let k = 0;
  for (let j = 0; j < K; j++) if (Math.abs(sigmaMid[j] - sigma) < Math.abs(sigmaMid[k] - sigma)) k = j;
  core.diagnoseColumn(i, pi, theta, none, none);
  theta[k * C + i] = T / core.diagnostics.exnerLayer[k * C + i];
  core.diagnoseColumn(i, pi, theta, none, none);
  if (humid) for (let j = 0; j < K; j++) if (sigmaMid[j] > 0.3) q[j * C + i] = 0.6 * saturationHumidity(theta[j * C + i] * core.diagnostics.exnerLayer[j * C + i], pi[i] * sigmaMid[j]);
  qc[k * C + i] = path * core.diagnostics.g / (pi[i] * core.diagnostics.dSigma[k]);
  core.diagnoseColumn(i, pi, theta, q, qc);
  const beam = SOLAR_CONSTANT * 0.5, surfaceT = 290;
  const run = (water) => { rad.column(i, pi[i], theta, surfaceT, 5, undefined, beam, q[(K - 1) * C + i], q, water, 0, 0); return { ...rad.budget }; };
  const cloudy = run(qc), clear = run(none);
  return { model, i, k, q, beam, surfaceT, cloudy, clear };
}

function referenceLongwave({ model, i, k, q, surfaceT }, cloudEmissivity, ozoneProfile) {
  const { K, sigmaMid, dSigma, g, exnerLayer } = model.core.diagnostics, C = model.mesh.nCells, [pi, theta] = model.state;
  const row = new Float64Array(6), surface = STEFAN_BOLTZMANN * surfaceT ** 4;
  let outgoing = 0;
  for (const point of LONGWAVE_TABLE.points) {
    let up = planckShare(point, surfaceT) * surface;
    for (let j = K - 1; j >= 0; j--) {
      const x = j * C + i, T = theta[x] * exnerLayer[x], dry = 1 - q[x];
      layerPaths(row, pi[i] * sigmaMid[j], pi[i] * dSigma[j] / g, T, q[x], ozoneProfile[j] * OZONE_CM_ATM, GREENHOUSE_GASES.carbonDioxide * GAS_MOLAR.co2 / GAS_MOLAR.air * dry, GREENHOUSE_GASES.methane * GAS_MOLAR.ch4 / GAS_MOLAR.air * dry, GREENHOUSE_GASES.nitrousOxide * GAS_MOLAR.n2o / GAS_MOLAR.air * dry);
      let tau = 0;
      for (let n = 0; n < 6; n++) tau += point[n] * row[n];
      const eps = 1 - Math.exp(-LONGWAVE_CONSTANTS.diffusivity * tau) * (j === k ? 1 - cloudEmissivity : 1);
      up = up * (1 - eps) + eps * planckShare(point, T) * STEFAN_BOLTZMANN * T ** 4;
    }
    outgoing += up;
  }
  return outgoing;
}

test('under transparent spectral gases a liquid and an ice layer reflect and emit with the hand values of their phase optics, as under the gray gases with an open window', () => {
  const lines = [];
  for (const [name, [T, path]] of Object.entries(CLOUDS)) {
    const sigma = T > 273 ? 0.8 : 0.25;
    for (const [label, options] of [['spectral', { ...SPECTRAL, ...TRANSPARENT }], ['gray', { ...GRAY_GASES, window: 1, gasFraction: 0, ozoneAbsorption: 0, vaporAbsorption: 0 }]]) {
      const c = cloudColumn(options, { T, path, humid: false, sigma });
      const surface = STEFAN_BOLTZMANN * c.surfaceT ** 4, layer = STEFAN_BOLTZMANN * T ** 4;
      const incident = c.beam - c.cloudy.atmosphereSolar;
      const reflectance = c.cloudy.reflectedSolar / incident, emissivity = (surface - c.cloudy.outgoingLongwave) / (surface - layer);
      lines.push(`${name} ${label}: R ${reflectance.toFixed(5)} ε ${emissivity.toFixed(5)} (gases take ${(c.cloudy.atmosphereSolar / c.beam * 100).toFixed(2)} % of the beam)`);
      assert.ok(Number.isFinite(c.cloudy.outgoingLongwave) && Number.isFinite(c.cloudy.absorbedSolar), `${name} ${label}: finite fluxes`);
      if (HAND[name].reflectance) assert.ok(close(reflectance, HAND[name].reflectance, 1e-12), `${name} ${label}: reflectance ${reflectance} against ${HAND[name].reflectance}`);
      if (HAND[name].emissivity) assert.ok(close(emissivity, HAND[name].emissivity, 1e-9), `${name} ${label}: emissivity ${emissivity} against ${HAND[name].emissivity}`);
      assert.ok(Math.abs(c.cloudy.absorbedSolar + c.cloudy.reflectedSolar - c.beam) < 1e-12 * c.beam, `${name} ${label}: the shortwave closes`);
    }
  }
  console.log(lines.join('; '));
});

test('under the spectral gases in a humid column the cloud joins every g-point as 1 − (1 − ε_gas)(1 − ε_cloud) with its hand emissivity, reflects its hand share of the light the gases leave, and the gases\' overlap is all that moves its effects from the cloud alone', () => {
  const lines = [];
  for (const [name, [T, path]] of Object.entries(CLOUDS)) {
    if (!HAND[name].emissivity && !HAND[name].reflectance) continue;
    const sigma = T > 273 ? 0.8 : 0.25;
    const probe = cloudColumn({ ...SPECTRAL }, { T, path, humid: true, sigma });
    const { K } = probe.model.core.diagnostics;
    const ozoneProfile = Float64Array.from({ length: K }, (_, j) => 0.3 * (probe.model.core.diagnostics.levels[j + 1] ** 3 - probe.model.core.diagnostics.levels[j] ** 3) * (j < K / 3 ? 1 : 0.1));
    const c = cloudColumn({ ...SPECTRAL, ozoneProfile }, { T, path, humid: true, sigma }), gray = cloudColumn({ ...GRAY_GASES }, { T, path, humid: true, sigma });
    const emissivity = HAND[name].emissivity ?? (1 - Math.exp(-1.66 * 1000 * 0.090361 * path));
    const reference = referenceLongwave(c, emissivity, ozoneProfile), referenceClear = referenceLongwave(c, 0, ozoneProfile);
    assert.ok(close(c.cloudy.outgoingLongwave, reference, 1e-12) && close(c.cloudy.clearOutgoingLongwave, referenceClear, 1e-12), `${name}: OLR ${c.cloudy.outgoingLongwave} against ${reference}, clear ${c.cloudy.clearOutgoingLongwave} against ${referenceClear}`);
    const surface = STEFAN_BOLTZMANN * c.surfaceT ** 4, alone = emissivity * (surface - STEFAN_BOLTZMANN * T ** 4);
    const longwaveEffect = c.cloudy.clearOutgoingLongwave - c.cloudy.outgoingLongwave, grayEffect = gray.cloudy.clearOutgoingLongwave - gray.cloudy.outgoingLongwave;
    assert.ok(longwaveEffect > 0 && longwaveEffect < alone, `${name}: the gases mask part of the cloud's ${alone} W/m2: ${longwaveEffect}`);
    const incident = c.beam - c.cloudy.atmosphereSolar, R = c.cloudy.reflectedSolar / incident;
    if (HAND[name].reflectance) assert.ok(close(R, HAND[name].reflectance, 1e-12), `${name}: reflectance of the light the gases leave ${R}`);
    assert.ok(Math.abs(c.cloudy.atmosphereSolar - c.clear.atmosphereSolar) < 1e-9 * c.beam, `${name}: the gases take the same light with and without the cloud (the whole column's path along the beam, as the gray vapour does)`);
    const shortwaveEffect = c.cloudy.absorbedSolar - c.cloudy.clearAbsorbedSolar, grayShortwave = gray.cloudy.absorbedSolar - gray.cloudy.clearAbsorbedSolar;
    assert.ok(close(shortwaveEffect, -R * incident, 1e-12), `${name}: SWCRE ${shortwaveEffect} against −R × incident ${-R * incident}`);
    assert.ok(Math.abs(c.cloudy.absorbedSolar + c.cloudy.reflectedSolar - c.beam) < 1e-12 * c.beam, `${name}: the shortwave closes`);
    lines.push(`${name}: LWCRE ${longwaveEffect.toFixed(2)} spectral / ${grayEffect.toFixed(2)} gray / ${alone.toFixed(2)} cloud alone W/m2; SWCRE ${shortwaveEffect.toFixed(2)} spectral / ${grayShortwave.toFixed(2)} gray (incident ${incident.toFixed(1)} / ${(gray.beam - gray.cloudy.atmosphereSolar).toFixed(1)} of ${c.beam.toFixed(1)})`);
  }
  console.log(lines.join('; '));
});

function cloudyGrid() {
  const model = createModel(new Grid(6), { ocean: false, divergenceDamping: 0, radiation: { mixedLayerDeck: false } });
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  model.step(900);
  for (const a of model.state) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  const [pi, theta, , , q, qc] = model.state, { K } = model.core, C = model.mesh.nCells;
  model.core.diagnose(pi, theta, q, qc);
  const { exnerLayer, dSigma, g, sigmaMid } = model.core.diagnostics;
  for (let i = 0; i < C; i++) for (const [goal, path] of [[282, 0.1], [222, 0.02]]) {
    let best = 0;
    for (let k = 0; k < K; k++) if (Math.abs(theta[k * C + i] * exnerLayer[k * C + i] - goal) < Math.abs(theta[best * C + i] * exnerLayer[best * C + i] - goal)) best = k;
    const x = best * C + i;
    q[x] = Math.fround(0.95 * saturationHumidity(theta[x] * exnerLayer[x], pi[i] * sigmaMid[best]));
    qc[x] = Math.fround(path * g / (pi[i] * dSigma[best]));
  }
  return model;
}

async function physicsPair(base, options, dt = 864000) {
  const physics = { mixedLayerDeck: false, clearSkyPass: true, ...options };
  const model = createModel(base.mesh, { ocean: false, radiation: physics });
  model.state.forEach((a, n) => a.set(base.state[n]));
  model.seaIce.concentration.set(base.seaIce.concentration);
  model.seaIce.snow.set(base.seaIce.snow);
  model.boundaryLayer.depth.set(base.boundaryLayer.depth);
  model.time = base.time;
  const { K } = model.core, C = model.mesh.nCells;
  const reference = Float64Array.from({ length: K }, (_, k) => { let s = 0, a = 0; for (let i = 0; i < C; i++) { s += model.mesh.areaCell[i] * model.state[1][k * C + i]; a += model.mesh.areaCell[i]; } return s / a; });
  const gpu = await createGpuCore(model.mesh, { nu4: model.core.nu4, nu4Theta: model.core.nu4Theta, divergenceDamping: model.core.divergenceDamping, referenceTheta: reference, physics });
  const { device, buffers, kernels, layout } = gpu;
  gpu.upload(model.state);
  gpu.uploadPhysics({ snow: model.seaIce.snow, concentration: model.seaIce.concentration });
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.DEPTH, Float32Array.from(model.boundaryLayer.depth));
  await gpu.tendency();
  const sun = sunDirection(model.time);
  device.queue.writeBuffer(buffers.P, 0, Float32Array.from([dt, 0, sun[0], sun[1], sun[2], 0, 0, 0]));
  const group = device.createBindGroup({ layout: kernels.physics.getBindGroupLayout(0), entries: [buffers.MI, buffers.MF, buffers.LV, buffers.S, buffers.K1, buffers.D, buffers.P, buffers.PH].map((buffer, binding) => ({ binding, resource: { buffer } })) });
  const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
  pass.setPipeline(kernels.physics);
  pass.setBindGroup(0, group);
  pass.dispatchWorkgroups(Math.ceil(C / 64));
  pass.end();
  device.queue.submit([encoder.finish()]);
  const [, after] = await gpu.download(), ph = await gpu.downloadPhysics();
  const before = Float64Array.from(model.state[1]), [pi, , , , , qc] = model.state;
  model.core.diagnose(pi, model.state[1], model.state[4], qc);
  const exner = Float64Array.from(model.core.diagnostics.exnerLayer);
  model.radiation.setTime(model.time);
  model.phases.physics(0, C, dt, null);
  const heating = (theta) => Float64Array.from(theta, (t, x) => (t - before[x]) * exner[x] / dt * 86400);
  const s = model.radiation.summed, slice = (name) => Float64Array.from(ph[name].subarray(0, C));
  return {
    K, C, lit: Array.from({ length: C }, (_, i) => i).filter((i) => model.radiation.insolation(i) > 0),
    cpu: { heating: heating(model.state[1]), surface: Float64Array.from(model.radiation.surfaceFlux), down: Float64Array.from(model.radiation.surfaceShortwave), absorbed: Float64Array.from(s.absorbedSolar), outgoing: Float64Array.from(s.outgoingLongwave), shortwaveEffect: Float64Array.from(s.absorbedSolar, (x, i) => x - s.clearAbsorbedSolar[i]), longwaveEffect: Float64Array.from(s.clearOutgoingLongwave, (x, i) => x - s.outgoingLongwave[i]) },
    gpu: { heating: heating(after), surface: slice('SFLUX'), down: slice('SWDN'), absorbed: slice('ABSSUM'), outgoing: slice('OLRSUM'), shortwaveEffect: Float64Array.from(slice('ABSSUM'), (x, i) => x - ph.ABSCLRSUM[i]), longwaveEffect: Float64Array.from(slice('OLRCLRSUM'), (x, i) => x - ph.OLRSUM[i]) },
  };
}

test('with a liquid and an ice layer in every column and every default on but the mixed-layer deck (the empirical deck in its place), the engines agree on the layer heating, the surface fluxes, the top fluxes and both cloud effects under the spectral gases, and on how far the spectral gases move them from the gray', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const base = cloudyGrid();
  const spectral = await physicsPair(base, {}), gray = await physicsPair(base, GRAY_GASES);
  const { K, C, lit } = spectral, runs = [spectral, gray];
  let heating = 0, scale = 0, movedApart = 0;
  for (let i = 0; i < C; i++) for (let k = 0; k < K - 1; k++) {
    const x = k * C + i;
    heating = Math.max(heating, ...runs.map((r) => Math.abs(r.cpu.heating[x] - r.gpu.heating[x])));
    scale = Math.max(scale, Math.abs(spectral.cpu.heating[x]));
    movedApart = Math.max(movedApart, Math.abs((spectral.cpu.heating[x] - gray.cpu.heating[x]) - (spectral.gpu.heating[x] - gray.gpu.heating[x])));
  }
  const keys = ['surface', 'down', 'absorbed', 'outgoing', 'shortwaveEffect', 'longwaveEffect'];
  const worst = (key) => Math.max(...runs.map((r) => Math.max(...Array.from(r.cpu[key], (v, i) => Math.abs(v - r.gpu[key][i])))));
  const size = (key) => Math.max(...runs.map((r) => Math.max(...Array.from(r.cpu[key], Math.abs))));
  const mean = (r, side, key) => Array.from({ length: C }, (_, i) => i).reduce((s, i) => s + base.mesh.areaCell[i] * r[side][key][i], 0) / base.mesh.areaCell.reduce((a, b) => a + b, 0);
  for (const r of runs) for (const side of ['cpu', 'gpu']) for (const key of [...keys, 'heating']) assert.ok(r[side][key].every(Number.isFinite), `${side} ${key}: finite`);
  const apart = Object.fromEntries(keys.map((key) => [key, worst(key)]));
  console.log(`${C} columns at N=6 (${lit.length} lit), a 100 g/m2 layer near 282 K and a 20 g/m2 layer near 222 K in each: layer heating apart by at most ${heating.toExponential(1)} K/day of ${scale.toFixed(1)}, the spectral gases' change of it alike to ${movedApart.toExponential(1)}; ${keys.map((key) => `${key} ${apart[key].toExponential(1)} of ${size(key).toFixed(0)}`).join(', ')} W/m2; means CPU / GPU, spectral then gray: SWCRE ${runs.map((r) => `${mean(r, 'cpu', 'shortwaveEffect').toFixed(2)} / ${mean(r, 'gpu', 'shortwaveEffect').toFixed(2)}`).join(', ')}, LWCRE ${runs.map((r) => `${mean(r, 'cpu', 'longwaveEffect').toFixed(2)} / ${mean(r, 'gpu', 'longwaveEffect').toFixed(2)}`).join(', ')}`);
  assert.ok(heating < 1e-5 * scale && movedApart < 1e-5 * scale, `layer heating apart by ${heating} K/day, the gases' change of it by ${movedApart}`);
  for (const [key, d] of Object.entries(apart)) assert.ok(d < 2e-5 * size(key), `${key} apart by ${d} W/m2 of ${size(key)}`);
});
