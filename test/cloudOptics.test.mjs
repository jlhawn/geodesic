import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { saturationHumidity } from '../js/physics/moist.module.js';
import { sunDirection, cloudOptics, liquidShare, iceRadius, CLOUD_OPTICS, SOLAR_CONSTANT, STEFAN_BOLTZMANN } from '../js/physics/radiation.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuCore } = gpuAvailable ? await import('../js/gpu/core.gpu.js') : {};
const UNSCATTERED = { rayleighDepth: 0, landAerosol: 0, seaAerosol: 0, skylight: 0.15 };
const close = (a, b, tolerance = 1e-9) => Math.abs(a - b) <= tolerance * Math.abs(b);

test('the optics of a kilogram of condensate: liquid above 273.15 K with the droplet radius of sea or land, ice below 235.15 K with its radius from the temperature, an even split halfway, against hand-computed values', () => {
  assert.deepEqual([liquidShare(300), liquidShare(273.15), liquidShare(235.15), liquidShare(200)], [1, 1, 0, 0]);
  assert.ok(close(liquidShare(254.15), 0.5, 1e-12));
  assert.ok(close(iceRadius(213.15), 15.55, 1e-12) && close(iceRadius(150), 15.55, 1e-12) && close(iceRadius(253.15), 73.55, 1e-12) && close(iceRadius(290), 73.55, 1e-12));
  const cases = [
    [285, false, { liquid: 1, visible: 127.11864406779661, solar: 18.01428813559323, infrared: 149.99926 }],
    [285, true, { liquid: 1, visible: 176.47058823529412, solar: 26.453470588235295, infrared: 149.99926 }],
    [213.15, false, { liquid: 0, visible: 159.78240514469385, solar: 35.919355507703884, infrared: 115.05241157556225 }],
    [254.15, false, { liquid: 0.5, visible: 81.80949470555024, solar: 12.490479608676, infrared: 90.43447024473147 }],
  ];
  for (const [T, continental, expected] of cases) {
    const got = cloudOptics(T, continental);
    for (const key of Object.keys(expected)) assert.ok(close(got[key], expected[key]), `${T} K ${continental ? 'land' : 'sea'}: ${key} ${got[key]} against ${expected[key]}`);
  }
  assert.ok(CLOUD_OPTICS.seaDropletRadius === 11.8 && CLOUD_OPTICS.landDropletRadius === 8.5);
});

function singleLayer(radiation, T, path) {
  const model = createModel(new Grid(4), { ocean: false, radiation: { cloudCover: 'overcast', window: 1, gasFraction: 0, ozoneAbsorption: 0, cloudSolarAbsorption: 0, ...UNSCATTERED, skylight: 0, ...radiation } });
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  const { core, radiation: rad, mesh } = model, { K } = core.diagnostics, C = mesh.nCells, [pi, theta] = model.state;
  const i = 0, k = Math.floor(K / 2), dry = new Float64Array(K * C), qc = new Float64Array(K * C);
  core.diagnoseColumn(i, pi, theta, dry, dry);
  theta[k * C + i] = T / core.diagnostics.exnerLayer[k * C + i];
  qc[k * C + i] = path * core.diagnostics.g / (pi[i] * core.diagnostics.dSigma[k]);
  core.diagnoseColumn(i, pi, theta, dry, qc);
  const mu = 0.5, beam = SOLAR_CONSTANT * mu, surfaceT = 290;
  rad.column(i, pi[i], theta, surfaceT, 5, undefined, beam, null, null, qc, 0, 0);
  const b = rad.budget, surface = STEFAN_BOLTZMANN * surfaceT ** 4;
  return { reflectance: b.reflectedSolar / beam, emissivity: (surface - b.outgoingLongwave) / (surface - STEFAN_BOLTZMANN * T ** 4), closure: Math.abs(b.absorbedSolar + b.reflectedSolar - beam) / beam };
}

test('an overcast layer over a black surface reflects (1 − g)τ / ((1 − g)τ + 2μ) and emits with 1 − exp(−κ W): 100 and 10 g/m² of liquid at 285 K over sea, 20 g/m² of ice at 213.15 K, against hand-computed values; the explicit gray optics override', () => {
  const liquid = singleLayer({}, 285, 0.1), thin = singleLayer({}, 285, 0.01), ice = singleLayer({}, 213.15, 0.02), gray = singleLayer({ cloudScattering: 95, cloudAbsorption: 130 }, 213.15, 0.02);
  console.log(`liquid 100 g/m²: reflectance ${liquid.reflectance.toFixed(5)} (hand 0.64304); liquid 10 g/m²: emissivity ${thin.emissivity.toFixed(5)} (hand 0.77687); ice 20 g/m²: reflectance ${ice.reflectance.toFixed(5)} (hand 0.41806), emissivity ${ice.emissivity.toFixed(5)} (hand 0.89985); the gray optics give the ice ${gray.reflectance.toFixed(5)} and ${gray.emissivity.toFixed(5)}`);
  assert.ok(close(liquid.reflectance, 0.6430392965333067, 1e-12), `${liquid.reflectance}`);
  assert.ok(close(thin.emissivity, 0.7768681886822757, 1e-9), `${thin.emissivity}`);
  assert.ok(close(ice.reflectance, 0.418058949528354, 1e-12) && close(ice.emissivity, 0.8998461956940361, 1e-9), `${ice.reflectance} ${ice.emissivity}`);
  assert.ok(close(gray.reflectance, 1.9 / 2.9, 1e-12) && close(gray.emissivity, -Math.expm1(-2.6), 1e-9), `${gray.reflectance} ${gray.emissivity}`);
  for (const r of [liquid, thin, ice, gray]) assert.ok(r.closure < 1e-14);
});

function cloudColumns() {
  const model = createModel(new Grid(6), { ocean: false, divergenceDamping: 0, radiation: { stratus: true, mixedLayerDeck: false, exchangeCoefficient: 1.5e-3 }, boundaryLayer: { entrainment: { efficiency: 0, shear: 0 }, dragCoefficient: 1.5e-3 }, surface: { dragCoefficient: 1.5e-3 } });
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  const { K, sigmaMid } = model.core, C = model.mesh.nCells;
  for (let k = 0; k < K; k++) if (sigmaMid[k] < 0.75) for (let i = 0; i < C; i++) model.state[1][k * C + i] += 10;
  model.step(900);
  for (const a of model.state) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  const [pi, theta, , , q, qc] = model.state;
  model.core.diagnose(pi, theta, q, qc);
  const { exnerLayer, dSigma, g } = model.core.diagnostics;
  let seed = 3;
  const random = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; };
  const counts = [0, 0, 0];
  for (let i = 0; i < C; i++) {
    for (const goal of [282, 252, 222]) {
      if (random() < 0.25) continue;
      let best = -1;
      for (let k = 0; k < K; k++) { const T = theta[k * C + i] * exnerLayer[k * C + i]; if (best < 0 || Math.abs(T - goal) < Math.abs(theta[best * C + i] * exnerLayer[best * C + i] - goal)) best = k; }
      const x = best * C + i, T = theta[x] * exnerLayer[x];
      if (Math.abs(T - goal) > 12) continue;
      const qs = saturationHumidity(T, pi[i] * sigmaMid[best]);
      q[x] = Math.fround(0.95 * qs);
      qc[x] = Math.fround((0.02 + 0.28 * random()) * g / (pi[i] * dSigma[best]));
      counts[T > 273.15 ? 0 : T > 235.15 ? 1 : 2]++;
    }
  }
  return { model, counts };
}

async function physicsPair(base, options, dt = 864000) {
  const physics = { stratus: true, mixedLayerDeck: false, clearSkyPass: true, ...options };
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
    K, C, decks: model.radiation.stratusFraction.filter((f) => f > 0).length, lit: Array.from({ length: C }, (_, i) => i).filter((i) => model.radiation.insolation(i) > 0),
    cpu: { heating: heating(model.state[1]), surface: Float64Array.from(model.radiation.surfaceFlux), down: Float64Array.from(model.radiation.surfaceShortwave), shortwaveEffect: Float64Array.from(s.absorbedSolar, (x, i) => x - s.clearAbsorbedSolar[i]), longwaveEffect: Float64Array.from(s.clearOutgoingLongwave, (x, i) => x - s.outgoingLongwave[i]) },
    gpu: { heating: heating(after), surface: slice('SFLUX'), down: slice('SWDN'), shortwaveEffect: Float64Array.from(slice('ABSSUM'), (x, i) => x - ph.ABSCLRSUM[i]), longwaveEffect: Float64Array.from(slice('OLRCLRSUM'), (x, i) => x - ph.OLRSUM[i]) },
    done: () => device.destroy?.(),
  };
}

test('with warm, mixed-phase and cold cloud and the empirical deck in the columns the engines agree on the layer heating, the surface fluxes and both cloud effects under the phase optics, and on how far the phase optics move them from the gray', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { model: base, counts } = cloudColumns();
  const phase = await physicsPair(base, {}), gray = await physicsPair(base, { cloudScattering: 95, cloudAbsorption: 130 });
  const { K, C, lit } = phase;
  let heating = 0, scale = 0, moved = 0, movedApart = 0;
  for (let i = 0; i < C; i++) for (let k = 0; k < K - 1; k++) {
    const x = k * C + i;
    heating = Math.max(heating, Math.abs(phase.cpu.heating[x] - phase.gpu.heating[x]), Math.abs(gray.cpu.heating[x] - gray.gpu.heating[x]));
    scale = Math.max(scale, Math.abs(phase.cpu.heating[x]));
    moved = Math.max(moved, Math.abs(phase.cpu.heating[x] - gray.cpu.heating[x]));
    movedApart = Math.max(movedApart, Math.abs((phase.cpu.heating[x] - gray.cpu.heating[x]) - (phase.gpu.heating[x] - gray.gpu.heating[x])));
  }
  const runs = [phase, gray];
  const worst = (key) => { let d = 0; for (const r of runs) for (let i = 0; i < C; i++) d = Math.max(d, Math.abs(r.cpu[key][i] - r.gpu[key][i])); return d; };
  const size = (key) => { let d = 0; for (const r of runs) for (let i = 0; i < C; i++) d = Math.max(d, Math.abs(r.cpu[key][i])); return d; };
  const largest = (key) => { let d = 0; for (const i of lit) d = Math.max(d, Math.abs(phase.cpu[key][i] - gray.cpu[key][i])); return d; };
  const apart = { surface: worst('surface'), down: worst('down'), shortwaveEffect: worst('shortwaveEffect'), longwaveEffect: worst('longwaveEffect') };
  console.log(`${C} columns at N=6 (${lit.length} lit) with ${counts[0]} warm, ${counts[1]} mixed-phase and ${counts[2]} cold cloudy layers (${phase.decks} empirical decks): the engines' layer heating differs by at most ${heating.toExponential(1)} K/day against a largest ${scale.toFixed(1)}; the phase optics move it from the gray by up to ${moved.toFixed(2)} K/day, alike to ${movedApart.toExponential(1)}; the net surface flux differs by ${apart.surface.toExponential(1)}, the surface sunlight by ${apart.down.toExponential(1)}, the shortwave and longwave cloud effects by ${apart.shortwaveEffect.toExponential(1)} and ${apart.longwaveEffect.toExponential(1)} W/m² (the largest of each ${['surface', 'down', 'shortwaveEffect', 'longwaveEffect'].map((key) => size(key).toFixed(0)).join(', ')} W/m²), which the phase optics move by up to ${largest('shortwaveEffect').toFixed(1)} and ${largest('longwaveEffect').toFixed(1)}`);
  assert.ok(counts.every((n) => n > C / 4) && phase.decks > 10, `cloudy layers ${counts}, ${phase.decks} decks`);
  assert.ok(moved > 1 && largest('shortwaveEffect') > 10 && largest('longwaveEffect') > 5, 'the phase optics move the heating and both effects');
  assert.ok(heating < 1e-5 * scale && movedApart < 1e-5 * scale, `layer heating apart by ${heating} K/day, its change by ${movedApart}`);
  for (const [key, d] of Object.entries(apart)) assert.ok(d < 2e-5 * size(key), `${key} apart by ${d} W/m² of ${size(key)}`);
});

test('under the phase optics every sunlit column with warm, mixed-phase and cold cloud closes its shortwave: absorbed and reflected add up to the beam and the layers take what the atmosphere absorbs', () => {
  const { model } = cloudColumns(), { radiation, core } = model, { K, C, geopotential, g } = core.diagnostics, [pi, theta, , surfaceT, q, qc] = model.state;
  model.radiation.useCumulus(null, null);
  const bottom = (K - 1) * C;
  let lit = 0, closure = 0, layers = 0, cloudHeat = 0;
  for (let i = 0; i < C; i++) {
    core.diagnoseColumn(i, pi, theta, q, qc);
    const beam = radiation.insolation(i);
    radiation.column(i, pi[i], theta, surfaceT[i], 5, undefined, beam, q[bottom + i], q, qc, 0.06, 0.06, 1, undefined, 1, model.boundaryLayer.depth[i] - geopotential[bottom + i] / g, 0, 0);
    if (!(beam > 1)) continue;
    const b = radiation.budget;
    lit++;
    cloudHeat = Math.max(cloudHeat, b.cloudSolar);
    closure = Math.max(closure, Math.abs(b.absorbedSolar + b.reflectedSolar - beam) / beam);
    let sum = 0;
    for (let k = 0; k < K; k++) sum += radiation.layerFlux[k] - radiation.longwave[k * C + i] - (k === K - 1 ? b.sensibleHeat : 0);
    layers = Math.max(layers, Math.abs(sum - b.atmosphereSolar) / beam);
  }
  console.log(`${lit} columns lit by more than 1 W/m²: absorbed + reflected equals the beam to ${closure.toExponential(1)} of it, the layers' shortwave heating the atmosphere's absorption to ${layers.toExponential(1)}; the cloud absorbs up to ${cloudHeat.toFixed(1)} W/m²`);
  assert.ok(lit > 100 && cloudHeat > 1);
  assert.ok(closure < 1e-12 && layers < 1e-12, `closure ${closure}, layers ${layers}`);
});
