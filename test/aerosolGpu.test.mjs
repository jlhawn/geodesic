import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { saturationHumidity } from '../js/physics/moist.module.js';
import { DAY } from '../js/physics/radiation.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { GRAVITY } from '../js/dynamics/sigmaCore.module.js';
import { layoutFor, PHYSICS_DEFAULTS } from '../js/gpu/core.gpu.js';
import { physicsConstants, physicsFunctions, PHYSICS_FUNCTIONS, physicsKernel, radiationKernel, layeredPass, LONGWAVE_STASH_FLOATS } from '../js/gpu/physics.gpu.js';
import { syntheticAerosol } from './helpers/aerosol.mjs';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuCore } = gpuAvailable ? await import('../js/gpu/core.gpu.js') : {};
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};
const N = 16, DT = 1350, mesh = buildMesh(new Grid(N)), C = mesh.nCells;
const MAC_TEXT = /mac[A-Z]|MAC_|PH_AER|light: vec4/;

/*
 * Radiation option aerosol 'mac' on the GPU: the climatology of the model day lands in the physics buffer
 * as the CPU engine's sampler has it, the kernels take the same slabs and slices of a cell as the CPU
 * radiation, and without the option the generated WGSL has none of it.
 */

test('without \'mac\' no WGSL of the aerosol is generated and the physics buffer and the stash are the engine\'s; with it they grow by the eight fields and an array', () => {
  const base = { ...PHYSICS_DEFAULTS, R: 287 }, mac = { ...base, aerosol: 'mac' };
  assert.equal(PHYSICS_DEFAULTS.aerosol, 'flat');
  assert.equal(PHYSICS_DEFAULTS.aerosolScale, 1);
  assert.ok(!MAC_TEXT.test(physicsConstants(base)) && !MAC_TEXT.test(physicsFunctions(base)), 'the constants and functions');
  assert.equal(physicsFunctions(base), PHYSICS_FUNCTIONS);
  for (const held of [false, true]) for (const layered of [false, true]) assert.ok(!MAC_TEXT.test(radiationKernel(held, layered)) && !MAC_TEXT.test(radiationKernel(held, layered, false)), `radiation kernel held ${held} layered ${layered}`);
  for (const skin of [false, true]) assert.ok(!MAC_TEXT.test(physicsKernel(skin)) && !MAC_TEXT.test(physicsKernel(skin, false)), `physics kernel with skin ${skin}`);
  for (const held of [false, true]) assert.ok(!MAC_TEXT.test(layeredPass(base, held)), `layered pass held ${held}`);
  assert.ok(MAC_TEXT.test(physicsConstants(mac)) && /macStreams/.test(physicsConstants(mac)) && /MAC_R: array<f32, 24>/.test(physicsConstants(mac)));
  assert.notEqual(physicsFunctions(mac), PHYSICS_FUNCTIONS);
  assert.ok(!/light: vec3<f32>/.test(physicsFunctions(mac)) && /NIR_RAY \* light\.y \+ light\.w/.test(physicsFunctions(mac)));
  for (const held of [false, true]) for (const layered of [false, true]) assert.ok(/macStreams/.test(radiationKernel(held, layered, true)), `radiation kernel held ${held} layered ${layered}`);
  assert.ok(/macGasAbove|gasAbove/.test(physicsKernel(false, true)) && /macF/.test(layeredPass(mac, false)) && /mq0/.test(layeredPass(mac, true)));
  const plain = layoutFor(mesh, 8), tall = layoutFor(mesh, 8, 0, 0, 0, null, false, 0, 0, false, false, false, false, true, 0, true);
  assert.equal(plain.PH.AER, undefined);
  assert.equal(tall.PH.AER, tall.PH.total - 8 * C, 'the eight fields are the last slot of the physics buffer');
  assert.equal(tall.PH.total, layoutFor(mesh, 8, 0, 0, 0, null, false, 0, 0, false, false, false, false, true, 0, false).PH.total + 8 * C);
  for (const name of Object.keys(plain.PH)) if (name !== 'total') assert.equal(layoutFor(mesh, 8, 0, 0, 0, null, false, 0, 0, false, false, false, false, true, 0, true).PH[name], layoutFor(mesh, 8, 0, 0, 0, null, false, 0, 0, false, false, false, false, true, 0, false).PH[name], `slot ${name}`);
  assert.equal(LONGWAVE_STASH_FLOATS(8, C, false), 5 * 8 * C + 26 * C);
});

test('the GPU refuses the options the CPU radiation refuses, with the same words', () => {
  const given = (extra) => physicsConstants({ ...PHYSICS_DEFAULTS, R: 287, ...extra });
  assert.doesNotThrow(() => given({ aerosol: 'mac' }));
  assert.doesNotThrow(() => given({ aerosol: 'mac', shortwave: 'sequential', nearInfrared: 'clirad', aerosolScale: 2.5 }));
  for (const [bad, message] of [[{ aerosol: 'dust' }, /aerosol must be 'flat' or 'mac', not dust/], [{ aerosol: 'mac', solarGases: 'lacisHansen' }, /aerosol 'mac' needs solarGases 'clirad', not lacisHansen/],
    [{ aerosol: 'mac', shortwave: 'sequential', vaporShielding: 'cloudTop' }, /aerosol 'mac' does not run with vaporShielding 'cloudTop'/], [{ aerosol: 'mac', aerosolScale: -1 }, /aerosolScale must be a factor of the aerosol's depth, 0 or more, not -1/]]) assert.throws(() => given(bad), message, JSON.stringify(bad));
});

function meanTheta(model) {
  const { K } = model.core, theta = model.state[1];
  return Float64Array.from({ length: K }, (_, k) => { let s = 0, a = 0; for (let i = 0; i < C; i++) { s += model.mesh.areaCell[i] * theta[k * C + i]; a += model.mesh.areaCell[i]; } return s / a; });
}

function cloudier(model, thickness = 1) {
  const { K, sigmaMid, exnerLayer } = model.core.diagnostics, [pi, theta, , , q, qc] = model.state;
  const layers = [[0.8, 0.95, 4e-4, (i) => i % 4 === 0], [0.45, 0.6, 2e-4, (i) => i % 7 < 3], [0.2, 0.3, 6e-5, (i) => i % 5 === 1], [0.55, 0.9, 1e-4, (i) => i % 11 === 2]];
  for (const [from, to, water, test] of layers) for (let k = 0; k < K; k++) {
    if (sigmaMid[k] < from || sigmaMid[k] > to) continue;
    for (let i = 0; i < C; i++) {
      if (!test(i)) continue;
      const idx = k * C + i;
      q[idx] = saturationHumidity(theta[idx] * exnerLayer[idx], pi[i] * sigmaMid[k]);
      qc[idx] = thickness * water;
    }
  }
}

function stats(cpu, gpu) {
  let maxDiff = 0, sumSq = 0, sumRef = 0;
  const errors = new Float64Array(cpu.length);
  for (let x = 0; x < cpu.length; x++) { const d = Math.abs(cpu[x] - gpu[x]); errors[x] = d; if (d > maxDiff) maxDiff = d; sumSq += d * d; sumRef += cpu[x] * cpu[x]; }
  const scale = Math.sqrt(sumRef / cpu.length);
  return { maxDiff, rmsRel: Math.sqrt(sumSq / Math.max(sumRef, 1e-300)), p99Rel: errors.sort()[Math.floor(0.99 * (cpu.length - 1))] / Math.max(scale, 1e-300) };
}

/*
 * The CPU model and a GPU core of the same options, stepped together `steps` times from model time `time`;
 * `aerosol` false is the same without the climatology.
 */
async function pair(aerosol, radiation, { steps = 3, time = 200 * DAY - 2 * DT, thickness = 1, data = syntheticAerosol().bytes, model: modelOptions = {} } = {}) {
  const options = { stratus: true, clearSkyPass: true, ...radiation, ...(aerosol ? { aerosol: 'mac', aerosolData: data } : {}) };
  const model = createModel(new Grid(N), { ocean: false, radiation: options, surface: { exchange: 'roughness' }, ...modelOptions });
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  cloudier(model, thickness);
  model.time = time;
  const gpu = await createGpuCore(model.mesh, { nu4: model.core.nu4, nu4Theta: model.core.nu4Theta, divergenceDamping: model.core.divergenceDamping, referenceTheta: meanTheta(model), physics: options });
  gpu.upload(model.state);
  gpu.uploadPhysics({});
  for (let n = 0; n < steps; n++) { const t = model.time; model.step(DT); await gpu.stepModel(DT, t); }
  const run = { model, gpu, state: await gpu.download(), physics: await gpu.downloadPhysics() };
  for (const buffer of Object.values(gpu.buffers)) buffer.destroy();
  return run;
}

function measure(label, { model, state, physics }) {
  const r = model.radiation;
  const m = {
    absorbed: stats(r.summed.absorbedSolar, physics.ABSSUM.subarray(0, C)), atmosphere: stats(r.summed.atmosphereSolar, physics.ATMSUM.subarray(0, C)), reflected: stats(r.summed.reflectedSolar, physics.REFLSUM.subarray(0, C)),
    clear: stats(r.summed.clearAbsorbedSolar, physics.ABSCLRSUM.subarray(0, C)), olr: stats(r.outgoing, physics.OLR.subarray(0, C)), theta: stats(model.state[1], state[1]),
  };
  console.log(`${label}: against the CPU, 99th percentile of the cells (rms): ${Object.entries(m).map(([key, v]) => `${key} ${v.p99Rel.toExponential(1)} (${v.rmsRel.toExponential(1)})`).join(', ')}`);
  return m;
}

const agrees = (mac, flat) => { for (const key of ['absorbed', 'atmosphere', 'reflected', 'clear', 'olr', 'theta']) assert.ok(mac[key].p99Rel < 3 * flat[key].p99Rel + 2e-4, `${key} 99th percentile ${mac[key].p99Rel} against ${flat[key].p99Rel} without the climatology`); };
const areaMean = (values) => { let area = 0, sum = 0; for (let i = 0; i < C; i++) { area += mesh.areaCell[i]; sum += mesh.areaCell[i] * values[i]; } return sum / area; };

for (const [label, radiation, options] of [['layered', { shortwave: 'layered' }, {}], ['sequential', { shortwave: 'sequential' }, {}], ['layered with the radiation held between calls', { shortwave: 'layered', radiationEvery: 4 }, { steps: 6 }],
  ['sequential with the CLIRAD near infrared', { shortwave: 'sequential', nearInfrared: 'clirad' }, {}], ['layered with the depths doubled', { shortwave: 'layered', aerosolScale: 2 }, {}]]) {
  test(`${label}, on the GPU the climatology of the model day lands in the physics buffer as the CPU sampler has it, and the model matches the CPU model at N = 16 as it does without the aerosol`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
    const flat = await pair(false, radiation, options), mac = await pair(true, radiation, options);
    const f = measure(`${label}, without the aerosol`, flat), m = measure(`${label}, aerosol 'mac'`, mac);
    const fields = mac.model.aerosol.fields, uploaded = mac.physics.AER.subarray(0, 8 * C);
    let worst = 0;
    for (let x = 0; x < fields.length; x++) worst = Math.max(worst, Math.abs(fields[x] - uploaded[x]));
    assert.equal(worst, 0, 'the eight fields on the device are the sampler\'s');
    assert.equal(mac.model.aerosol.day, 200, 'the model crossed into day 200');
    const cpu = [areaMean(flat.model.radiation.summed.absorbedSolar), areaMean(mac.model.radiation.summed.absorbedSolar)], gpu = [areaMean(flat.physics.ABSSUM), areaMean(mac.physics.ABSSUM)];
    const atmosphere = [[flat, mac].map(({ model }) => areaMean(model.radiation.summed.atmosphereSolar)), [flat, mac].map(({ physics }) => areaMean(physics.ATMSUM))];
    console.log(`area-mean absorbed sunlight, W/m2 per step sum: CPU ${cpu[0].toFixed(2)} to ${cpu[1].toFixed(2)}, GPU ${gpu[0].toFixed(2)} to ${gpu[1].toFixed(2)}; atmosphere CPU ${atmosphere[0][0].toFixed(2)} to ${atmosphere[0][1].toFixed(2)}, GPU ${atmosphere[1][0].toFixed(2)} to ${atmosphere[1][1].toFixed(2)}`);
    assert.ok(Math.abs(atmosphere[0][1] - atmosphere[0][0]) > 1 && Math.abs(cpu[1] - cpu[0]) > 0.5, 'the climatology moves the sunlight in the CPU model');
    assert.ok(Math.abs((gpu[1] - gpu[0]) - (cpu[1] - cpu[0])) < 0.02 * Math.abs(cpu[1] - cpu[0]) + 0.02 && Math.abs((atmosphere[1][1] - atmosphere[1][0]) - (atmosphere[0][1] - atmosphere[0][0])) < 0.02 * Math.abs(atmosphere[0][1] - atmosphere[0][0]) + 0.02, 'the GPU takes the same change');
    agrees(m, f);
  });
}

test('twelve steps on the GPU across a change of day and of the month\'s weights track the CPU model, whose fields they take from the same sampler', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const data = syntheticAerosol().bytes, start = (15 - 78 + 365) * DAY - 5 * DT;
  const flat = await pair(false, { shortwave: 'layered' }, { steps: 12, time: start }), mac = await pair(true, { shortwave: 'layered' }, { steps: 12, time: start, data });
  const f = measure('twelve steps, without the aerosol', flat), m = measure('twelve steps, aerosol \'mac\'', mac);
  assert.equal(mac.model.aerosol.day, 15 - 78 + 365);
  let worst = 0;
  for (let x = 0; x < 8 * C; x++) worst = Math.max(worst, Math.abs(mac.model.aerosol.fields[x] - mac.physics.AER[x]));
  assert.equal(worst, 0);
  agrees(m, f);
});

test('the change \'mac\' makes to the shortwave heating of each layer of the sunlit columns, clear and cloudy, agrees between the engines, and so does the change in the absorbed sunlight', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const data = syntheticAerosol().bytes;
  for (const shortwave of ['sequential', 'layered']) {
    const radiation = { shortwave, radiationEvery: 4, mixedLayerDeck: false };
    const options = { steps: 1, time: 200 * DAY, data };
    const off = await pair(false, radiation, options), on = await pair(true, radiation, options), { K } = on.model.core, { cp, g, dSigma } = on.model.core.diagnostics, pi = on.model.state[0];
    const rate = (run, k, i) => run.model.radiation.held.mu[i] > 0 ? 86400 * g / (cp * pi[i] * dSigma[k]) : 0;
    let worst = 0, largest = 0, at = null, parted = 0, lit = 0, worstTop = 0, strongestTop = 0, worstBase = 0, worstLower = 0, worstFlux = 0;
    for (let i = 0; i < C; i++) {
      if (!(on.model.radiation.held.mu[i] > 0)) continue;
      lit++;
      let apart = 0;
      for (let k = 0; k < K; k++) apart = Math.max(apart, Math.abs(off.model.radiation.held.shortwave[k * C + i] - off.physics.RADSW[k * C + i]) * rate(off, k, i));
      worstBase = Math.max(worstBase, apart);
      if (apart > 1e-3) { parted++; continue; }
      for (let k = 0; k < K; k++) {
        const x = k * C + i, r = rate(on, k, i), flux = (on.model.radiation.held.shortwave[x] - off.model.radiation.held.shortwave[x]), change = flux * r, miss = Math.abs(flux - (on.physics.RADSW[x] - off.physics.RADSW[x])), d = miss * r;
        largest = Math.max(largest, Math.abs(change));
        if (d > worst) { worst = d; at = [i, k, change]; }
        worstFlux = Math.max(worstFlux, miss);
        if (on.model.core.sigmaMid[k] * pi[i] > 10000) worstLower = Math.max(worstLower, d);
      }
      const cpuTop = on.model.radiation.held.top[i] - off.model.radiation.held.top[i], gpuTop = on.physics.RADTOA[i] - off.physics.RADTOA[i];
      worstTop = Math.max(worstTop, Math.abs(cpuTop - gpuTop) / Math.max(1, on.model.radiation.held.mu[i] * 1362));
      strongestTop = Math.max(strongestTop, Math.abs(cpuTop));
    }
    console.log(`${shortwave}, one full call at N=${N} on ${lit} sunlit columns (${parted} left out, where a discrete decision parted the engines without the aerosol too; their baseline heating otherwise differs by at most ${worstBase.toExponential(1)} K/day): the aerosol changes a layer's shortwave heating by up to ${largest.toFixed(3)} K/day, and the engines' change differs by at most ${worst.toExponential(1)} K/day (cell ${at[0]} layer ${at[1]}, ${at[2].toFixed(3)} K/day), below 100 hPa by ${worstLower.toExponential(1)}, as a flux in any layer by ${worstFlux.toExponential(1)} W/m2; the change of the column's absorbed sunlight by ${worstTop.toExponential(1)} of its beam (up to ${strongestTop.toFixed(1)} W/m2)`);
    assert.ok(lit > 0.3 * C && largest > 0.02 && parted < 0.1 * lit, `${lit} sunlit columns, ${parted} parted, largest change ${largest} K/day`);
    assert.ok(worstLower < 5e-4 && worstFlux < 1e-3, `the aerosol's heating differs by ${worst} K/day at cell ${at[0]}, layer ${at[1]}, below 100 hPa by ${worstLower}`);
    assert.ok(worstTop < 1e-5, `the change in absorbed sunlight differs by ${worstTop} of the beam`);
  }
});

test('over a 4 km mountain the engines agree on what the aerosol does to the layers, those more than 5 km above the ground included, and on the sunlight each column absorbs', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const mountain = syntheticTopography(180, 360, (lat, lon) => {
    const d = Math.hypot(lat - 0.6, (lon - 0.5) * Math.cos(0.6));
    return d < 0.5 ? 1 + 4000 * Math.exp(-((d / 0.18) ** 2)) : -4000;
  }), data = syntheticAerosol().bytes, runs = [];
  for (const aerosol of [false, true]) {
    const options = { stratus: true, shortwave: 'layered', radiationEvery: 4, mixedLayerDeck: false, ...(aerosol ? { aerosol: 'mac', aerosolData: data } : {}) };
    const cpu = createModel(new Grid(N), { ocean: false, topography: mountain, radiation: options }), gpu = await createGpuModel(new Grid(N), { ocean: false, topography: mountain, radiation: options });
    const { K, R } = cpu.core.diagnostics, phis = cpu.surfaceGeopotential;
    for (let i = 0; i < C; i++) cpu.state[0][i] = 101325 * Math.exp(-phis[i] / (R * 280));
    cpu.core.diagnose(cpu.state[0], cpu.state[1], null, null);
    for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) cpu.state[1][k * C + i] = 280 / cpu.core.arrays.exnerLayer[k * C + i];
    cpu.state[3].fill(280);
    cloudier(cpu);
    cpu.state.forEach((a, n) => gpu.state[n].set(a));
    gpu.load();
    cpu.time = gpu.time = 200 * DAY;
    cpu.step(DT);
    await gpu.step(DT);
    await gpu.settle();
    const ph = await gpu.gpu.downloadPhysics(), heights = Float64Array.from({ length: K * C }, (_, x) => (cpu.core.diagnostics.geopotential[x] - phis[x % C]) / GRAVITY / 1000);
    runs.push({ held: Float64Array.from(cpu.radiation.held.shortwave), gpu: Float64Array.from(ph.RADSW.subarray(0, K * C)), mu: Float64Array.from(cpu.radiation.held.mu), heights, top: Float64Array.from(cpu.radiation.held.top), gpuTop: Float64Array.from(ph.RADTOA.subarray(0, C)), K });
    gpu.destroy();
  }
  const [off, on] = runs, { K } = on;
  let cells = 0, worst = 0, above = 0, total = 0, worstTop = 0, aboveGpu = 0, totalGpu = 0;
  const mountainCells = Array.from({ length: C }, (_, i) => i).filter((i) => on.mu[i] > 0 && mesh.latCell[i] > 0.35 && mesh.latCell[i] < 0.85 && Math.abs(mesh.lonCell[i] - 0.5) < 0.35);
  for (const i of mountainCells) {
    cells++;
    for (let k = 0; k < K; k++) {
      const x = k * C + i, cpuChange = on.held[x] - off.held[x], gpuChange = on.gpu[x] - off.gpu[x];
      worst = Math.max(worst, Math.abs(cpuChange - gpuChange));
      total += cpuChange; totalGpu += gpuChange;
      if (on.heights[x] > 5) { above += cpuChange; aboveGpu += gpuChange; }
    }
    worstTop = Math.max(worstTop, Math.abs((on.top[i] - off.top[i]) - (on.gpuTop[i] - off.gpuTop[i])));
  }
  console.log(`4 km mountain: ${cells} sunlit columns around it; the aerosol, with what it does to the light that climbs back through the gases, adds ${total.toFixed(1)} W/m2 to the layers' shortwave heating on the CPU and ${totalGpu.toFixed(1)} on the GPU, ${(100 * above / total).toFixed(1)} % and ${(100 * aboveGpu / totalGpu).toFixed(1)} % of it more than 5 km above the ground; the engines' change differs by at most ${worst.toExponential(1)} W/m2 in a layer and ${worstTop.toExponential(1)} W/m2 in the absorbed sunlight`);
  assert.ok(cells > 5 && total > 10, `${cells} columns, ${total} W/m2`);
  assert.ok(Math.abs(totalGpu - total) < 0.01 * total && Math.abs(above - aboveGpu) < 0.01 * total, 'the engines place the same heating');
  assert.ok(above > 0.02 * total && above < 0.5 * total, `${above} of ${total} W/m2 more than 5 km above the ground`);
  assert.ok(worst < 0.05 && worstTop < 0.05, `layer ${worst}, top ${worstTop} W/m2`);
});
