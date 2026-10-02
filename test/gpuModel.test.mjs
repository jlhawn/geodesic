import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { sunDirection } from '../js/physics/radiation.module.js';
import { saturationHumidity } from '../js/physics/moist.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuCore } = gpuAvailable ? await import('../js/gpu/core.gpu.js') : {};
const UNSCATTERED = { rayleighDepth: 0, landAerosol: 0, seaAerosol: 0, skylight: 0.15, upwardAbsorption: false };

function meanTheta(model) {
  const { K } = model.core, C = model.mesh.nCells, theta = model.state[1];
  return Float64Array.from({ length: K }, (_, k) => { let s = 0, a = 0; for (let i = 0; i < C; i++) { s += model.mesh.areaCell[i] * theta[k * C + i]; a += model.mesh.areaCell[i]; } return s / a; });
}
function stats(cpu, gpu) {
  let maxDiff = 0, sumSq = 0, sumRef = 0, at = -1;
  for (let x = 0; x < cpu.length; x++) { const d = Math.abs(cpu[x] - gpu[x]); if (d > maxDiff) { maxDiff = d; at = x; } sumSq += d * d; sumRef += cpu[x] * cpu[x]; }
  return { maxDiff, at, rms: Math.sqrt(sumSq / cpu.length), rmsRel: Math.sqrt(sumSq / Math.max(sumRef, 1e-300)) };
}
async function pair(N, steps, dt, inversion = 0, stratus = inversion > 0, options = {}, moist = {}, boundaryLayer = {}) {
  const model = createModel(new Grid(N), { ocean: false, radiation: { stratus, ...options }, moist, boundaryLayer });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  const { K, sigmaMid } = model.core, C = model.mesh.nCells;
  for (let k = 0; k < K; k++) if (sigmaMid[k] < 0.75) for (let i = 0; i < C; i++) model.state[1][k * C + i] += inversion;
  const gpu = await createGpuCore(model.mesh, { nu4: model.core.nu4, nu4Theta: model.core.nu4Theta, divergenceDamping: model.core.divergenceDamping, referenceTheta: meanTheta(model), physics: { stratus, ...options, ...moist, ...boundaryLayer } });
  gpu.upload(model.state);
  gpu.uploadPhysics();
  for (let n = 0; n < steps; n++) { const time = model.time; model.step(dt); await gpu.stepModel(dt, time); }
  return { model, gpu, state: await gpu.download(), physics: await gpu.downloadPhysics() };
}

test('one full GPU step with physics matches the CPU model', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { model, state, physics: after } = await pair(6, 1, 900, 10);
  const theta = stats(model.state[1], state[1]), q = stats(model.state[4], state[4]), qc = stats(model.state[5], state[5]);
  const ts = stats(model.state[3], state[3]), ice = stats(model.state[6], state[6]), u = stats(model.state[2], state[2]);
  const concentration = stats(model.seaIce.concentration, after.CONC.subarray(0, model.mesh.nCells));
  console.log(`one step at N=6: θ rms ${theta.rmsRel.toExponential(1)} max ${theta.maxDiff.toExponential(1)} K; q rms ${q.rmsRel.toExponential(1)} max ${q.maxDiff.toExponential(1)}; qc max ${qc.maxDiff.toExponential(1)}; Ts max ${ts.maxDiff.toExponential(1)} K; ice max ${ice.maxDiff.toExponential(1)} m; concentration max ${concentration.maxDiff.toExponential(1)}; wind max ${u.maxDiff.toExponential(1)} m/s`);
  assert.ok(theta.rmsRel < 1e-5, `θ rms ${theta.rmsRel}`);
  assert.ok(theta.maxDiff < 0.05, `θ max ${theta.maxDiff} K at ${theta.at}`);
  assert.ok(q.rmsRel < 5e-4, `q rms ${q.rmsRel}`);
  assert.ok(ts.maxDiff < 0.02, `Ts max ${ts.maxDiff} K at ${ts.at}`);
  assert.ok(ice.maxDiff < 1e-3, `ice max ${ice.maxDiff} m`);
  assert.ok(concentration.maxDiff < 1e-3, `concentration max ${concentration.maxDiff} at ${concentration.at}`);
  assert.ok(u.maxDiff < 1e-2, `wind max ${u.maxDiff} m/s`);
  const second = await pair(6, 2, 900, 10, true, { mixedLayerDeck: false }), C = second.model.mesh.nCells, { stratus, stratusFraction } = second.model.radiation;
  const olr = stats(second.model.radiation.outgoing, second.physics.OLR.subarray(0, C)), sw = stats(second.model.radiation.surfaceShortwave, second.physics.SWDN.subarray(0, C));
  let decked = 0, water = 0, cover = 0;
  for (let i = 0; i < C; i++) { if (stratus[i] > 0) { decked++; water += stratus[i]; } cover += stratusFraction[i]; }
  console.log(`two steps at N=6 under a 10 K inversion (the EIS deck rests on the first step's boundary layer): stratus on ${(100 * decked / C).toFixed(1)} % of sea cells, ${(1000 * water / Math.max(1, decked)).toFixed(1)} g/m² where it forms, mean cover ${(cover / C).toFixed(3)} and mean water ${(1000 * water / C).toFixed(1)} g/m² over the sea; per-cell OLR rms ${olr.rmsRel.toExponential(1)}, surface shortwave rms ${sw.rmsRel.toExponential(1)}`);
  assert.ok(decked > 0.1 * C, `stratus on ${decked} of ${C} sea cells`);
  assert.ok(olr.rmsRel < 1e-5 && sw.rmsRel < 2e-5, `per-cell OLR rms ${olr.rmsRel}, surface shortwave rms ${sw.rmsRel}; a few thin EIS decks, whose single-precision water paths differ by 4e-4, carry most of it`);
  const third = await pair(6, 2, 900, 10, true, { stratusIndex: 'ectei', mixedLayerDeck: false }), entrained = third.model.radiation;
  const olrE = stats(entrained.outgoing, third.physics.OLR.subarray(0, C)), swE = stats(entrained.surfaceShortwave, third.physics.SWDN.subarray(0, C)), coverE = stats(entrained.stratusFraction, third.physics.DECKF.subarray(0, C));
  let coverSum = 0, waterSum = 0;
  for (let i = 0; i < C; i++) { coverSum += entrained.stratusFraction[i]; waterSum += entrained.stratus[i]; }
  console.log(`the same under ECTEI: mean cover ${(coverSum / C).toFixed(3)} and mean water ${(1000 * waterSum / C).toFixed(1)} g/m²; per-cell cover max difference ${coverE.maxDiff.toExponential(1)}, OLR rms ${olrE.rmsRel.toExponential(1)}, surface shortwave rms ${swE.rmsRel.toExponential(1)}`);
  assert.ok(coverSum > 0 && coverSum < cover, `mean cover ${coverSum / C} under ECTEI against ${cover / C}`);
  assert.ok(olrE.rmsRel < 1e-5 && swE.rmsRel < 1e-5, `per-cell OLR rms ${olrE.rmsRel}, surface shortwave rms ${swE.rmsRel} under ECTEI`);
});

test('the stratiform share of the estimated inversion strength and the boundary-layer entrainment it tapers match between the engines, with the cover\'s overcast bound or without it', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  for (const [inversion, options] of [[6, {}], [9, {}], [9, { overcastWater: null }]]) {
    const { model, physics } = await pair(6, 1, 900, inversion, true, options, {}, { turbulence: 'dry' });
    const C = model.mesh.nCells, share = model.radiation.stratiform, we = model.boundaryLayer.entrainment;
    const strat = stats(share, physics.STRAT.subarray(0, C)), entrain = stats(we, physics.ENTRAIN.subarray(0, C));
    let ramp = 0, full = 0, entraining = 0, tapered = 0;
    for (let i = 0; i < C; i++) { if (share[i] > 0 && share[i] < 1) ramp++; if (share[i] >= 1) full++; if (we[i] > 0) { entraining++; if (share[i] > 0) tapered++; } }
    console.log(`one step at N=6 under a ${inversion} K inversion${options.overcastWater === null ? ' without the overcast bound' : ''}: stratiform share on the ramp in ${ramp} and whole in ${full} of ${C} columns, max engine difference ${strat.maxDiff.toExponential(1)}; ${entraining} columns entrain, ${tapered} of them tapered, w_e max difference ${(1000 * entrain.maxDiff).toExponential(1)} mm/s`);
    assert.ok(ramp > C / 20, `${ramp} columns on the ramp`);
    assert.ok(strat.maxDiff < 2e-3, `share differs by ${strat.maxDiff} at ${strat.at}`);
    assert.ok(entrain.maxDiff < 2e-5, `w_e differs by ${entrain.maxDiff} at ${entrain.at}`);
    for (let i = 0; i < C; i++) if (share[i] >= 1) assert.equal(we[i], 0, `cell ${i}: an EIS of 12 K or more entrains nothing`);
  }
});

test('with the ∇⁴ closures off, the divergence damping alone and the heat it returns match between the engines', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const run = async (divergenceDamping) => {
    const model = createModel(new Grid(6), { ocean: false, nu4Hours: Infinity, divergenceDamping });
    const init = initializeState(model, {});
    for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
    const gpu = await createGpuCore(model.mesh, { nu4: 0, nu4Theta: 0, divergenceDamping, referenceTheta: meanTheta(model) });
    gpu.upload(model.state);
    gpu.uploadPhysics();
    for (let n = 0; n < 4; n++) { const time = model.time; model.step(900); await gpu.stepModel(900, time); }
    return { model, state: await gpu.download() };
  };
  const damped = await run(0.1), free = await run(0);
  assert.equal(damped.model.core.nu4, 0);
  const theta = stats(damped.model.state[1], damped.state[1]), u = stats(damped.model.state[2], damped.state[2]);
  const moved = stats(free.model.state[2], damped.model.state[2]), heated = stats(free.model.state[1], damped.model.state[1]);
  console.log(`four steps at N=6, c = 0.1 and no ∇⁴: engines differ in θ by rms ${theta.rmsRel.toExponential(1)}, max ${theta.maxDiff.toExponential(1)} K, in wind by at most ${u.maxDiff.toExponential(1)} m/s; the damping moves the wind by up to ${moved.maxDiff.toFixed(3)} m/s and θ by up to ${heated.maxDiff.toExponential(1)} K`);
  assert.ok(theta.rmsRel < 3e-7 && theta.maxDiff < 2e-3, `θ rms ${theta.rmsRel}, max ${theta.maxDiff} K at ${theta.at}; without the damping's heat on the GPU they differ by 8e-7 and 6e-3 K`);
  assert.ok(u.maxDiff < 0.02 * moved.maxDiff, `engines differ in wind by ${u.maxDiff} m/s, the damping moved it by ${moved.maxDiff}`);
});

/*
 * One physics kernel alone against the CPU physics phase, both from the
 * same single-precision state (one step after a 10 K inversion was
 * imposed): the heating of every layer, (θ' − θ)Π/dt, over a step long
 * enough that the rounding of θ is far below the tolerance.
 */
function heatingState(radiation = {}) {
  const model = createModel(new Grid(6), { ocean: false, divergenceDamping: 0, radiation: { stratus: true, mixedLayerDeck: false, exchangeCoefficient: 1.5e-3, ...radiation }, boundaryLayer: { entrainment: { efficiency: 0, shear: 0 }, dragCoefficient: 1.5e-3 }, surface: { dragCoefficient: 1.5e-3 } });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  const { K, sigmaMid } = model.core, C = model.mesh.nCells;
  for (let k = 0; k < K; k++) if (sigmaMid[k] < 0.75) for (let i = 0; i < C; i++) model.state[1][k * C + i] += 10;
  model.step(900);
  for (const a of model.state) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  return model;
}
async function physicsHeating(base, options, dt = 864000, cumulus = null, mixingTop = null) {
  const physics = { stratus: true, mixedLayerDeck: false, ...options };
  const model = createModel(base.mesh, { ocean: false, radiation: physics, moist: physics, ice: physics });
  if (cumulus) { model.moist.cumulusCover.set(cumulus.cover); model.moist.cumulusWater.set(cumulus.water); }
  model.state.forEach((a, n) => a.set(base.state[n]));
  model.seaIce.concentration.set(base.seaIce.concentration);
  model.seaIce.snow.set(base.seaIce.snow);
  model.boundaryLayer.depth.set(base.boundaryLayer.depth);
  const buoyancy = Float64Array.from({ length: model.mesh.nCells }, (_, i) => (i % 2 ? 1e-4 : -1e-4));
  if (mixingTop) { model.boundaryLayer.mixingTop.set(mixingTop); model.boundaryLayer.buoyancyFlux.set(buoyancy); }
  model.time = base.time;
  const { K } = model.core, C = model.mesh.nCells;
  const gpu = await createGpuCore(model.mesh, { nu4: model.core.nu4, nu4Theta: model.core.nu4Theta, divergenceDamping: model.core.divergenceDamping, referenceTheta: meanTheta(model), physics });
  const { device, buffers, kernels, layout } = gpu;
  gpu.upload(model.state);
  gpu.uploadPhysics({ snow: model.seaIce.snow, concentration: model.seaIce.concentration });
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.DEPTH, Float32Array.from(model.boundaryLayer.depth));
  if (mixingTop) { device.queue.writeBuffer(buffers.PH, 4 * layout.PH.MIXTOP, Float32Array.from(mixingTop)); device.queue.writeBuffer(buffers.PH, 4 * layout.PH.BUOY, Float32Array.from(buoyancy)); }
  if (cumulus) {
    const layers = (layout.PH.CUWATER - layout.PH.CUCOVER) / C;
    device.queue.writeBuffer(buffers.PH, 4 * layout.PH.CUCOVER, Float32Array.from(cumulus.cover.subarray((K - layers) * C)));
    device.queue.writeBuffer(buffers.PH, 4 * layout.PH.CUWATER, Float32Array.from(cumulus.water.subarray((K - layers) * C)));
  }
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
  const { exnerLayer, dSigma, g, cp } = model.core.diagnostics, exner = Float64Array.from(exnerLayer);
  model.radiation.setTime(model.time);
  model.phases.physics(0, C, dt, null);
  const heating = (theta) => Float64Array.from(theta, (t, x) => (t - before[x]) * exner[x] / dt * 86400);
  const cloudy = [];
  for (let i = 0; i < C; i++) {
    let water = 0;
    for (let k = 0; k < K; k++) water += qc[k * C + i] * pi[i] * dSigma[k] / g;
    if (model.radiation.insolation(i) > 0 && (water > 1e-3 || model.radiation.stratus[i] > 0)) cloudy.push(i);
  }
  const power = (rate, i) => { let sum = 0; for (let k = 0; k < K; k++) sum += rate[k * C + i] / 86400 * cp * pi[i] * dSigma[k] / g; return sum; };
  return { K, C, cloudy, cpu: heating(model.state[1]), gpu: heating(after), power, area: model.mesh.areaCell, cpuDeck: Float64Array.from(model.radiation.stratusFraction), gpuDeck: ph.DECKF.subarray(0, C), cpuLongwave: Float64Array.from(model.radiation.longwave), gpuLongwave: ph.LWH.subarray(0, K * C) };
}

test('the heating of each layer of the sunlit cloudy columns, and the part of it the cloud water absorbs, agree between the engines', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const base = heatingState(UNSCATTERED);
  const lit = await physicsHeating(base, UNSCATTERED), scattering = await physicsHeating(base, { cloudSolarAbsorption: 0, ...UNSCATTERED });
  const { K, C, cloudy } = lit;
  let worst = 0, surface = 0, worstCloud = 0, at = null, cpuMean = 0, gpuMean = 0, area = 0, strongest = 0, decks = 0;
  for (const i of cloudy) {
    if (lit.cpuDeck[i] > 0) decks++;
    for (let k = 0; k < K; k++) {
      const x = k * C + i, d = Math.abs(lit.cpu[x] - lit.gpu[x]);
      if (k === K - 1) surface = Math.max(surface, d / Math.abs(lit.cpu[x]));
      else if (d > worst) { worst = d; at = [i, k, lit.cpu[x]]; }
      worstCloud = Math.max(worstCloud, Math.abs((lit.cpu[x] - scattering.cpu[x]) - (lit.gpu[x] - scattering.gpu[x])));
    }
  }
  for (let i = 0; i < C; i++) {
    const cpuCloud = lit.power(lit.cpu, i) - scattering.power(scattering.cpu, i), gpuCloud = lit.power(lit.gpu, i) - scattering.power(scattering.gpu, i);
    cpuMean += lit.area[i] * cpuCloud; gpuMean += lit.area[i] * gpuCloud; area += lit.area[i];
    strongest = Math.max(strongest, cpuCloud);
  }
  console.log(`physics alone at N=6 on ${cloudy.length} sunlit cloudy columns (${decks} with a deck): layer heating differs between the engines by at most ${worst.toExponential(1)} K/day (cell ${at[0]} layer ${at[1]}, ${at[2].toFixed(2)} K/day), in the lowest layer, which takes the sensible heat, by ${surface.toExponential(1)} of it; the cloud water's own heating by ${worstCloud.toExponential(1)} K/day. The cloud water absorbs a global mean ${(cpuMean / area).toFixed(3)} W/m² (GPU ${(gpuMean / area).toFixed(3)}), at most ${strongest.toFixed(1)} W/m² in a column`);
  assert.ok(cloudy.length > 0.1 * C && decks > 0, `${cloudy.length} cloudy columns, ${decks} with a deck`);
  for (let i = 0; i < C; i++) assert.ok(Math.abs(lit.cpuDeck[i] - lit.gpuDeck[i]) < 1e-5, `deck cover of cell ${i}`);
  assert.ok(worst < 1e-4, `layer heating differs by ${worst} K/day at cell ${at[0]}, layer ${at[1]}`);
  assert.ok(surface < 1e-5, `the lowest layer's heating differs by ${surface} of it`);
  assert.ok(worstCloud < 1e-4, `the cloud water's heating differs by ${worstCloud} K/day`);
  assert.ok(cpuMean / area > 0.1 && Math.abs(gpuMean - cpuMean) < 1e-3 * cpuMean, `global cloud absorption ${cpuMean / area} against ${gpuMean / area} W/m²`);
});

test('the change the Rayleigh and aerosol scattering makes to the heating of each layer of the sunlit columns agrees between the engines', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const base = heatingState(UNSCATTERED), greyIce = { iceAlbedo: 0.5, meltingIceAlbedo: 0.5, snowAgeing: false };
  const on = await physicsHeating(base, greyIce), off = await physicsHeating(base, { ...UNSCATTERED, ...greyIce });
  const { K, C } = on;
  let worst = 0, largest = 0, at = null, lit = 0;
  for (let i = 0; i < C; i++) {
    if (!(base.radiation.insolation(i) > 0)) continue;
    lit++;
    for (let k = 0; k < K - 1; k++) {
      const x = k * C + i, change = on.cpu[x] - off.cpu[x], d = Math.abs(change - (on.gpu[x] - off.gpu[x]));
      largest = Math.max(largest, Math.abs(change));
      if (d > worst) { worst = d; at = [i, k, change]; }
    }
  }
  console.log(`physics alone at N=6 on ${lit} sunlit columns: the scattering changes a layer's heating by up to ${largest.toFixed(3)} K/day, and the engines' change differs by at most ${worst.toExponential(1)} K/day (cell ${at[0]} layer ${at[1]}, ${at[2].toFixed(3)} K/day)`);
  assert.ok(lit > 0.3 * C && largest > 0.02, `${lit} sunlit columns, largest change ${largest} K/day`);
  assert.ok(worst < 2e-5, `the scattering's heating differs by ${worst} K/day at cell ${at[0]}, layer ${at[1]}`);
});

test('shallow cumulus of partial cover beside and without resolved cloud heats the layers of sunlit columns alike in both engines', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const base = heatingState(), { K, sigmaMid } = base.core, C = base.mesh.nCells;
  const cover = new Float64Array(K * C), water = new Float64Array(K * C);
  let seed = 7;
  const random = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; };
  let layers = 0;
  for (let i = 0; i < C; i++) {
    if (random() < 0.4) continue;
    for (let k = 0; k < K; k++) {
      if (!(sigmaMid[k] > 0.8 && sigmaMid[k] < 0.97)) continue;
      cover[k * C + i] = Math.fround(0.02 + 0.2 * random()); water[k * C + i] = Math.fround(2e-4 + 1e-3 * random());
      layers++;
    }
  }
  const tiny = Float64Array.from(cover, (f, x) => (f > 0 && x % 3 === 0 ? Math.fround(10 ** (-9 + 3 * random())) : f));
  const plain = await physicsHeating(base, {}), cumulus = await physicsHeating(base, {}, 864000, { cover, water }), sparse = await physicsHeating(base, {}, 864000, { cover: tiny, water });
  let worst = 0, moved = 0, at = null;
  for (let i = 0; i < C; i++) {
    if (!(base.radiation.insolation(i) > 0)) continue;
    for (let k = 0; k < K - 1; k++) {
      const x = k * C + i, d = Math.max(Math.abs(cumulus.cpu[x] - cumulus.gpu[x]), Math.abs(sparse.cpu[x] - sparse.gpu[x]));
      if (d > worst) { worst = d; at = [i, k]; }
      moved = Math.max(moved, Math.abs(cumulus.cpu[x] - plain.cpu[x]));
    }
  }
  console.log(`${layers} cumulus layers of cover 0.02–0.22, a third of them 10⁻⁹–10⁻⁶ in a second run: the engines' layer heating differs by at most ${worst.toExponential(1)} K/day (cell ${at[0]} layer ${at[1]}); the cumulus moves it by up to ${moved.toFixed(2)} K/day`);
  assert.ok(moved > 0.1, `the cumulus moves the heating by ${moved} K/day`);
  assert.ok(worst < 2e-4, `layer heating differs by ${worst} K/day`);
});

function cloudyState() {
  const base = heatingState(), { mesh, core, state } = base;
  const [pi, theta, , , q, qc] = state, C = mesh.nCells, { K } = core;
  core.diagnose(pi, theta, q, qc);
  const { exnerLayer, sigmaMid, geopotential, g } = core.diagnostics, depth = base.boundaryLayer.depth;
  let mid = 0;
  for (let k = 0; k < K; k++) if (Math.abs(sigmaMid[k] - 0.5) < Math.abs(sigmaMid[mid] - 0.5)) mid = k;
  let inside = 0, above = 0;
  for (let i = 0; i < C; i++) {
    const layers = [];
    for (let k = K - 3; k >= 0; k--) {
      const low = geopotential[k * C + i] / g < depth[i];
      if (low && layers.length === 0) layers.push(k);
      if (!low) { layers.push(k); break; }
    }
    layers.push(mid);
    layers.forEach((k, n) => {
      const x = k * C + i, qs = saturationHumidity(theta[x] * exnerLayer[x], pi[i] * sigmaMid[k]);
      q[x] = Math.fround((0.9 + 0.01 * ((i * 7 + n * 3) % 11)) * qs);
      qc[x] = Math.fround((0.002 + 0.08 * (((i * 5 + n) % 13) / 12)) * qs);
      if (geopotential[x] / g < depth[i]) inside++; else above++;
    });
  }
  let separated = 0;
  for (let i = 0; i < C; i++) {
    let runs = 0;
    for (let k = 0; k < K; k++) if (qc[k * C + i] > 0 && !(k > 0 && qc[(k - 1) * C + i] > 0)) runs++;
    if (runs > 1) separated++;
  }
  return { base, inside, above, separated };
}

test('resolved cloud of partial cover under the saturation adjustment, inside the boundary layer and above it and in separate runs of layers, heats the layers of sunlit columns alike in both engines, and the cover, its boundary-layer RHc, the overlap and the condensate bound, on its inversion ramp as well as off it, move that heating', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { base, inside, above, separated } = cloudyState(), old = { condensation: 'saturation', cloudOverlap: 'maximumRandom' };
  const pdf = await physicsHeating(base, old), moved = await physicsHeating(base, { ...old, boundaryCriticalHumidity: 0.5 }), overcast = await physicsHeating(base, { ...old, cloudCover: 'overcast' });
  const maximum = await physicsHeating(base, { ...old, cloudOverlap: 'maximum' }), unbounded = await physicsHeating(base, { ...old, overcastWater: null });
  const ramp = await physicsHeating(base, { ...old, overcastInversion: [-40, 40], overcastWater: 5e-4 }), full = await physicsHeating(base, { ...old, overcastInversion: [-1000, -999], overcastWater: 5e-4 });
  const { K, C } = pdf, lit = [];
  for (let i = 0; i < C; i++) if (base.radiation.insolation(i) > 0) lit.push(i);
  let engines = 0, scale = 0, cover = 0, boundary = 0, overlap = 0, bound = 0, belowRamp = 0, aboveRamp = 0;
  for (const i of lit) {
    for (let k = 0; k < K - 1; k++) {
      const x = k * C + i;
      engines = Math.max(engines, ...[pdf, moved, maximum, unbounded, ramp, full].map((r) => Math.abs(r.cpu[x] - r.gpu[x])));
      scale = Math.max(scale, Math.abs(pdf.cpu[x]));
      cover = Math.max(cover, Math.abs(pdf.cpu[x] - overcast.cpu[x]));
      boundary = Math.max(boundary, Math.abs(pdf.cpu[x] - moved.cpu[x]));
      overlap = Math.max(overlap, Math.abs(pdf.cpu[x] - maximum.cpu[x]));
      bound = Math.max(bound, Math.abs(pdf.cpu[x] - unbounded.cpu[x]));
      belowRamp = Math.max(belowRamp, Math.abs(ramp.cpu[x] - unbounded.cpu[x]));
      aboveRamp = Math.max(aboveRamp, Math.abs(ramp.cpu[x] - full.cpu[x]));
    }
  }
  console.log(`${lit.length} sunlit columns with ${inside} cloudy layers inside the boundary layer and ${above} above, ${separated} columns with separate runs: the engines' layer heating differs by at most ${engines.toExponential(1)} K/day against a largest ${scale.toFixed(1)}; the cover moves it by up to ${cover.toFixed(2)} K/day from overcast, the boundary layer's RHc of 0.5 by ${boundary.toFixed(2)}, maximum overlap by ${overlap.toFixed(2)}, the unbounded half-width by ${bound.toFixed(2)}, a bound of 5·10⁻⁴ kg/kg on an inversion ramp of −40 to 40 K by ${belowRamp.toFixed(2)} from the unbounded and ${aboveRamp.toFixed(2)} from the same bound everywhere`);
  assert.ok(inside > C / 4 && above > C && separated > C / 2, `${inside} cloudy layers inside, ${above} above, ${separated} columns with separate runs`);
  assert.ok(cover > 1 && boundary > 0.1 && overlap > 0.1 && bound > 0.1 && belowRamp > 0.1 && aboveRamp > 0.1, `cover ${cover}, boundary-layer RHc ${boundary}, overlap ${overlap}, bound ${bound}, ramp ${belowRamp} and ${aboveRamp} K/day`);
  assert.ok(engines < 1e-5 * scale, `layer heating differs by ${engines} K/day against ${scale}`);
});

test('the uniform condensation\'s cover of each layer\'s condensate, over ice where it is cold, and the exponential-random overlap heat the layers of sunlit columns alike in both engines and move that heating from the saturation adjustment\'s cover and maximum-random overlap', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { base } = cloudyState();
  const uniform = await physicsHeating(base, {}), liquid = await physicsHeating(base, { iceSaturation: false }), random = await physicsHeating(base, { cloudOverlap: 'maximumRandom' }), saturation = await physicsHeating(base, { condensation: 'saturation', iceSaturation: false });
  const ramp = await physicsHeating(base, { overcastInversion: [-40, 40], overcastWater: 5e-4 }), unbounded = await physicsHeating(base, { overcastWater: null });
  const { K, C } = uniform;
  let engines = 0, scale = 0, overlap = 0, cover = 0, ice = 0, blended = 0;
  for (let i = 0; i < C; i++) {
    if (!(base.radiation.insolation(i) > 0)) continue;
    for (let k = 0; k < K - 1; k++) {
      const x = k * C + i;
      engines = Math.max(engines, ...[uniform, liquid, random, ramp, unbounded].map((r) => Math.abs(r.cpu[x] - r.gpu[x])));
      scale = Math.max(scale, Math.abs(uniform.cpu[x]));
      overlap = Math.max(overlap, Math.abs(uniform.cpu[x] - random.cpu[x]));
      cover = Math.max(cover, Math.abs(random.cpu[x] - saturation.cpu[x]));
      ice = Math.max(ice, Math.abs(uniform.cpu[x] - liquid.cpu[x]));
      blended = Math.max(blended, Math.abs(ramp.cpu[x] - unbounded.cpu[x]));
    }
  }
  console.log(`the engines' layer heating differs by at most ${engines.toExponential(1)} K/day against a largest ${scale.toFixed(1)}; the exponential-random overlap moves it by up to ${overlap.toFixed(2)} K/day from maximum-random, the uniform cover by ${cover.toFixed(2)} from the saturation adjustment's, its saturation over ice by ${ice.toFixed(2)}, the stratiform blend at the cover's saturation on an inversion ramp of −40 to 40 K by ${blended.toFixed(2)}`);
  assert.ok(overlap > 1e-3 && cover > 0.1 && ice > 0.01 && blended > 0.1, `overlap ${overlap}, cover ${cover}, ice ${ice}, blend ${blended} K/day`);
  assert.ok(engines < 1e-5 * scale, `layer heating differs by ${engines} K/day against ${scale}`);
});

test('under maximum-random overlap a layer of trace cloud water joins the layers either side into one block in both engines', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { base } = cloudyState(), traced = cloudyState().base;
  const { K } = base.core, C = base.mesh.nCells, qc = traced.state[5];
  let filled = 0;
  for (let i = 0; i < C; i++) {
    let top = -1, bottom = -1;
    for (let k = 0; k < K; k++) if (qc[k * C + i] > 0) { if (top < 0) top = k; bottom = k; }
    for (let k = top + 1; k < bottom; k++) if (!(qc[k * C + i] > 0)) { qc[k * C + i] = Math.fround(1e-14); filled++; }
  }
  const options = { cloudOverlap: 'maximumRandom', condensation: 'saturation', iceSaturation: false };
  const plain = await physicsHeating(base, options), joined = await physicsHeating(traced, options);
  let engines = 0, scale = 0, moved = 0;
  for (let i = 0; i < C; i++) {
    if (!(base.radiation.insolation(i) > 0)) continue;
    for (let k = 0; k < K - 1; k++) {
      const x = k * C + i;
      engines = Math.max(engines, Math.abs(joined.cpu[x] - joined.gpu[x]));
      scale = Math.max(scale, Math.abs(joined.cpu[x]));
      moved = Math.max(moved, Math.abs(joined.cpu[x] - plain.cpu[x]));
    }
  }
  console.log(`${filled} gaps between cloud layers filled with 1e-14 of cloud water: joining the blocks moves the layer heating by up to ${moved.toFixed(2)} K/day; the engines differ by at most ${engines.toExponential(1)} K/day against a largest ${scale.toFixed(1)}`);
  assert.ok(filled > C && moved > 0.1, `${filled} gaps, heating moved ${moved} K/day`);
  assert.ok(engines < 1e-5 * scale, `layer heating differs by ${engines} K/day against ${scale}`);
});

test('the variance cover of the cloudy layers below the moist boundary layer\'s mixing top, and the longwave heating the cloud-top scheme reads, agree between the engines', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { base } = cloudyState(), { geopotential, g } = base.core.diagnostics, C = base.mesh.nCells, K = base.core.K, depth = base.boundaryLayer.depth;
  const mixingTop = Float64Array.from(depth, (d, i) => Math.fround(i % 3 === 0 ? 0 : d + (i % 3) * 400));
  const pdf = await physicsHeating(base, {}), variance = await physicsHeating(base, {}, 864000, null, mixingTop), off = await physicsHeating(base, { boundaryCover: 'pdf' }, 864000, null, mixingTop);
  const ramp = await physicsHeating(base, { overcastInversion: [-40, 40], overcastWater: 5e-4 }, 864000, null, mixingTop);
  let engines = 0, scale = 0, moved = 0, unmoved = 0, longwave = 0, longwaveScale = 0, inside = 0, blended = 0;
  for (let i = 0; i < C; i++) {
    for (let k = 0; k < K; k++) if (base.state[5][k * C + i] > 0 && geopotential[k * C + i] / g < mixingTop[i]) inside++;
    if (!(base.radiation.insolation(i) > 0)) continue;
    for (let k = 0; k < K - 1; k++) {
      const x = k * C + i;
      engines = Math.max(engines, Math.abs(variance.cpu[x] - variance.gpu[x]), Math.abs(ramp.cpu[x] - ramp.gpu[x]));
      blended = Math.max(blended, Math.abs(ramp.cpu[x] - variance.cpu[x]));
      scale = Math.max(scale, Math.abs(variance.cpu[x]));
      moved = Math.max(moved, Math.abs(variance.cpu[x] - pdf.cpu[x]));
      unmoved = Math.max(unmoved, Math.abs(off.cpu[x] - pdf.cpu[x]));
    }
  }
  for (let x = 0; x < K * C; x++) { longwave = Math.max(longwave, Math.abs(variance.cpuLongwave[x] - variance.gpuLongwave[x])); longwaveScale = Math.max(longwaveScale, Math.abs(variance.cpuLongwave[x])); }
  console.log(`${inside} cloudy layers below the mixing top: the variance cover moves the layer heating by up to ${moved.toFixed(2)} K/day against the humidity PDF, its blend into the overcast bound on an inversion ramp of −40 to 40 K by ${blended.toFixed(2)}; the engines differ by ${engines.toExponential(1)} K/day against a largest ${scale.toFixed(1)}; the longwave heating each layer keeps differs by ${longwave.toExponential(1)} W/m² against a largest ${longwaveScale.toFixed(1)}`);
  assert.ok(inside > C / 4, `${inside} cloudy layers inside`);
  assert.ok(moved > 0.1 && unmoved === 0 && blended > 0.1, `variance moves ${moved}, 'pdf' ${unmoved}, the ramp ${blended}`);
  assert.ok(engines < 5e-5 * scale, `layer heating differs by ${engines} K/day against ${scale}; the cover reads the f32 difference of q_t and q_s, as the bounded half-width does`);
  assert.ok(longwave < 1e-4 * longwaveScale, `longwave differs by ${longwave} W/m²`);
});

test('under the moist boundary layer the deck gated by its coupled stratocumulus, with the mixed-layer model or bypassed, matches between the engines', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  for (const options of [{ deckRegime: 'boundaryLayer' }, { deckRegime: 'boundaryLayer', deckBypass: true }, {}]) {
    const run = await mixedLayerPair(4, { turbulence: 'moist', ...options });
    const { K: nK } = run.model.core, keep = [];
    let coupled = 0, open = 0, parted = 0;
    for (let i = 0; i < run.C; i++) {
      if (run.model.boundaryLayer.regime[i] === 3) coupled++;
      if (run.model.radiation.mlmGate[i] > 0.5) open++;
      let worst = 0;
      for (let k = 0; k < nK; k++) worst = Math.max(worst, Math.abs(run.model.state[4][k * run.C + i] - run.state[4][k * run.C + i]));
      if (worst > 2e-5) parted++; else keep.push(i);
    }
    const pick = (a, b) => [Float64Array.from(keep.flatMap((i) => Array.from({ length: nK }, (_, k) => a[k * run.C + i]))), Float64Array.from(keep.flatMap((i) => Array.from({ length: nK }, (_, k) => b[k * run.C + i])))];
    const theta = stats(...pick(run.model.state[1], run.state[1])), q = stats(...pick(run.model.state[4], run.state[4]));
    console.log(`${JSON.stringify(options)}: four steps, ${coupled} of ${run.C} columns coupled stratocumulus, the gate open on ${open}, the deck on ${run.decked} (GPU ${run.gpuDecked}); engines differ in the gate by ${run.gate.maxDiff.toExponential(1)}, in cover by ${run.cover.maxDiff.toExponential(1)}, in deck water by rms ${run.mlmWater.rmsRel.toExponential(1)}; ${parted} columns part by more than 2·10⁻⁵ in q, where a cloud top or a parcel crosses its threshold in one engine only; elsewhere θ rms ${theta.rmsRel.toExponential(1)}, q rms ${q.rmsRel.toExponential(1)}; OLR rms ${run.olr.rmsRel.toExponential(1)}`);
    if (options.deckBypass) assert.ok(run.decked === 0 && run.gpuDecked === 0 && open > 0, `bypassed: ${run.decked} decks, gate open on ${open}`);
    else assert.ok(Math.abs(run.decked - run.gpuDecked) <= run.C / 100, `deck on ${run.decked}, GPU ${run.gpuDecked}`);
    let gateFlips = 0;
    for (let i = 0; i < run.C; i++) if (Math.abs(run.model.radiation.mlmGate[i] - run.gpuGate[i]) > 1e-3) gateFlips++;
    assert.ok(gateFlips <= run.C / 50, `${gateFlips} gates part`);
    assert.ok(parted <= run.C / 20, `${parted} columns part`);
    assert.ok(theta.rmsRel < 1e-5 && q.rmsRel < 1e-3 && run.olr.rmsRel < 1e-2, `θ ${theta.rmsRel}, q ${q.rmsRel}, OLR ${run.olr.rmsRel}`);
  }
});

test('twelve full GPU steps track the CPU model and its energy budget', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { model, state, physics } = await pair(6, 12, 900);
  const read = model.diagnostics(), d = { ...read, ...read.instantaneous };
  const C = model.mesh.nCells;
  let area = 0, ts = 0, abs = 0, atmosphere = 0, olr = 0, rain = 0;
  for (let i = 0; i < C; i++) { const a = model.mesh.areaCell[i]; area += a; ts += a * state[3][i]; abs += a * physics.ABS[i]; atmosphere += a * physics.ATMSW[i]; olr += a * physics.OLR[i]; rain += a * physics.RAIN[i]; }
  ts /= area; abs /= area; atmosphere /= area; olr /= area; rain /= area;
  const theta = stats(model.state[1], state[1]), tsStat = stats(model.state[3], state[3]);
  console.log(`twelve steps at N=6: mean Ts ${d.meanSurfaceT.toFixed(3)} vs ${ts.toFixed(3)} K; solar ${d.absorbedSolar.toFixed(2)} vs ${abs.toFixed(2)}, in the atmosphere ${d.atmosphereSolar.toFixed(2)} vs ${atmosphere.toFixed(2)}; OLR ${d.outgoingLongwave.toFixed(2)} vs ${olr.toFixed(2)} W/m²; θ rms ${theta.rmsRel.toExponential(1)}, Ts max ${tsStat.maxDiff.toExponential(1)} K`);
  assert.ok(Math.abs(d.meanSurfaceT - ts) < 0.02, `mean Ts ${d.meanSurfaceT} vs ${ts}`);
  assert.ok(Math.abs(d.absorbedSolar - abs) < 0.5, `absorbed solar ${d.absorbedSolar} vs ${abs}`);
  assert.ok(d.atmosphereSolar > 0 && Math.abs(d.atmosphereSolar - atmosphere) < 0.5, `absorbed in the atmosphere ${d.atmosphereSolar} vs ${atmosphere}`);
  assert.ok(Math.abs(d.outgoingLongwave - olr) < 0.5, `OLR ${d.outgoingLongwave} vs ${olr}`);
  assert.ok(theta.rmsRel < 1e-4, `θ rms ${theta.rmsRel}`);
  let snowOnIce = 0, worstSnow = 0;
  for (let i = 0; i < C; i++) {
    worstSnow = Math.max(worstSnow, Math.abs(model.seaIce.snow[i] - physics.SNOW[i]));
    if (physics.SNOW[i] > 0) snowOnIce++;
  }
  const ice = stats(model.state[6], state[6]), concentration = stats(model.seaIce.concentration, physics.CONC.subarray(0, C));
  console.log(`snow on ${snowOnIce} iced sea cells, engines differ by at most ${worstSnow.toExponential(1)} kg/m²; ice by ${ice.maxDiff.toExponential(1)} m, concentration by ${concentration.maxDiff.toExponential(1)}; ice fraction ${d.iceFraction.toFixed(4)}`);
  assert.ok(worstSnow < 1e-3, `snow on the surface differs between engines by ${worstSnow}`);
  assert.ok(ice.maxDiff < 1e-3 && concentration.maxDiff < 1e-3, `ice ${ice.maxDiff} m at ${ice.at}, concentration ${concentration.maxDiff} at ${concentration.at}`);
});

test('snow-ice formation matches between the engines', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const model = createModel(new Grid(6), { ocean: false });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  const C = model.mesh.nCells;
  let loaded = 0;
  for (let i = 0; i < C; i++) if (model.state[6][i] > 0) { model.seaIce.snow[i] = 100 + 200 * (i % 3); loaded++; }
  assert.ok(loaded > 0);
  const gpu = await createGpuCore(model.mesh, { nu4: model.core.nu4, nu4Theta: model.core.nu4Theta, divergenceDamping: model.core.divergenceDamping, referenceTheta: meanTheta(model) });
  gpu.upload(model.state);
  gpu.uploadPhysics();
  gpu.uploadLand({ soil: new Float64Array(C), snow: model.seaIce.snow, vegetation: new Float64Array(C) });
  const time = model.time; model.step(900); await gpu.stepModel(900, time);
  const state = await gpu.download(), physics = await gpu.downloadPhysics();
  let flooded = 0, worstSnow = 0, worstIce = 0;
  for (let i = 0; i < C; i++) if (model.state[6][i] > 0) {
    worstSnow = Math.max(worstSnow, Math.abs(model.seaIce.snow[i] - physics.SNOW[i]));
    worstIce = Math.max(worstIce, Math.abs(model.state[6][i] - state[6][i]));
    if (model.seaIce.snow[i] < 100) flooded++;
  }
  console.log(`snow-ice on ${flooded} of ${loaded} loaded cells; engines differ by ${worstSnow.toExponential(1)} kg/m² of snow and ${worstIce.toExponential(1)} m of ice`);
  assert.ok(flooded > 0, 'the heavy load floods somewhere');
  assert.ok(worstSnow < 1e-2 && worstIce < 1e-4, `snow ${worstSnow}, ice ${worstIce}`);
});

test('partly covered ice matches between the engines', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const model = createModel(new Grid(6), { ocean: false });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  const C = model.mesh.nCells;
  let seeded = 0;
  for (let i = 0; i < C; i++) if (model.state[6][i] > 0) { model.seaIce.concentration[i] = 0.5; model.seaIce.snow[i] = 10 * (i % 3); seeded++; }
  assert.ok(seeded > 0);
  const gpu = await createGpuCore(model.mesh, { nu4: model.core.nu4, nu4Theta: model.core.nu4Theta, divergenceDamping: model.core.divergenceDamping, referenceTheta: meanTheta(model) });
  gpu.upload(model.state);
  gpu.uploadPhysics();
  gpu.uploadLand({ soil: new Float64Array(C), snow: model.seaIce.snow, vegetation: new Float64Array(C) });
  gpu.uploadIce(model.seaIce.concentration);
  for (let n = 0; n < 4; n++) { const time = model.time; model.step(900); await gpu.stepModel(900, time); }
  const state = await gpu.download(), physics = await gpu.downloadPhysics();
  let worstArea = 0, worstIce = 0, worstT = 0, closed = 0, opened = 0;
  for (let i = 0; i < C; i++) if (init[6][i] > 0) {
    worstArea = Math.max(worstArea, Math.abs(model.seaIce.concentration[i] - physics.CONC[i]));
    worstIce = Math.max(worstIce, Math.abs(model.state[6][i] - state[6][i]));
    worstT = Math.max(worstT, Math.abs(model.state[3][i] - state[3][i]));
    if (model.seaIce.concentration[i] > 0.5) closed++; else if (model.seaIce.concentration[i] < 0.5) opened++;
  }
  console.log(`four steps from half cover on ${seeded} cells (${closed} closing, ${opened} opening): engines differ by ${worstArea.toExponential(1)} in concentration, ${worstIce.toExponential(1)} m of ice, ${worstT.toExponential(1)} K of skin`);
  assert.ok(closed + opened > 0, 'the concentration moved');
  assert.ok(worstArea < 1e-3 && worstIce < 1e-3 && worstT < 0.02, `concentration ${worstArea}, ice ${worstIce}, skin ${worstT}`);
});

/*
 * Every column mixed from the surface to σ 0.85 (the lowest layer's θ
 * and q throughout) under a 10 K inversion, with the running-mean
 * subsidence seeded to `seed` in both engines so the gate reads each
 * column's own inversion: the second step, the first with a diagnosed
 * boundary layer, carries the deck.
 */
async function mixedLayerPair(steps, { seed = -1e-3, height = 0, moist = { cloudLifetime: 3 * 3600, plumeCape: 70 }, step = null, turbulence = 'dry', ...options } = {}) {
  const physics = { mixedLayerDeck: true, deckRest: 'depth', minimumInversion: 2, ...options };
  const model = createModel(new Grid(6), { ocean: false, radiation: physics, moist, boundaryLayer: { turbulence } });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  const { K, sigmaMid } = model.core, C = model.mesh.nCells, theta = model.state[1], q = model.state[4];
  for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) {
    if (sigmaMid[k] < 0.85) theta[k * C + i] += 10;
    else { theta[k * C + i] = theta[(K - 1) * C + i] + (step && sigmaMid[k] < step ? 3 : 0); q[k * C + i] = q[(K - 1) * C + i]; }
  }
  model.radiation.mlmSubsidence.fill(seed);
  model.radiation.mlmHeight.fill(height);
  const gpu = await createGpuCore(model.mesh, { nu4: model.core.nu4, nu4Theta: model.core.nu4Theta, divergenceDamping: model.core.divergenceDamping, referenceTheta: meanTheta(model), physics: { ...physics, ...moist, turbulence } });
  gpu.upload(model.state);
  gpu.uploadPhysics({ mlmSubsidence: model.radiation.mlmSubsidence, mlmHeight: model.radiation.mlmHeight });
  for (let n = 0; n < steps; n++) { const time = model.time; model.step(900); await gpu.stepModel(900, time); }
  const after = await gpu.downloadPhysics(), r = model.radiation;
  const cell = (name) => after[name].subarray(0, C);
  const sunlit = [];
  for (let i = 0; i < C; i++) if (r.mlmCover[i] > 0 && r.insolation(i) > 200) sunlit.push(i);
  const pick = (values) => Float64Array.from(sunlit, (i) => values[i]);
  let decked = 0, gpuDecked = 0, partial = 0, water = 0;
  for (let i = 0; i < C; i++) {
    if (r.mlmCover[i] > 0) { decked++; water += r.mlmWater[i]; if (r.mlmCover[i] < 1) partial++; }
    if (after.MLMCOVER[i] > 0) gpuDecked++;
  }
  return {
    C, decked, gpuDecked, partial, water: water / Math.max(1, decked),
    cover: stats(r.mlmCover, cell('MLMCOVER')), mlmWater: stats(r.mlmWater, cell('MLMWATER')), entrainment: stats(r.mlmEntrainment, cell('MLMENT')),
    subsidence: stats(r.mlmSubsidence, cell('MLMSUB')), olr: stats(r.outgoing, cell('OLR')), sw: stats(r.surfaceShortwave, cell('SWDN')),
    height: stats(r.mlmHeight, cell('MLMH')), gate: stats(r.mlmGate, cell('MLMGATE')), gpuGate: Float64Array.from(cell('MLMGATE')), top: stats(r.mlmTop, cell('MLMTOP')), state: await gpu.download(), model,
    heights: Float64Array.from(r.mlmHeight), tops: Float64Array.from(r.mlmTop), depth: Float64Array.from(model.boundaryLayer.depth),
    fraction: stats(r.stratusFraction, cell('DECKF')), deck: stats(r.stratus, cell('DECK')), mean: r.mlmSubsidence,
    sunlit, sunlitWater: stats(pick(r.mlmWater), pick(cell('MLMWATER'))), sunlitDeck: stats(pick(r.stratus), pick(cell('DECK'))), waterPath: Float64Array.from(r.mlmWater),
  };
}

test('the mixed-layer deck matches between the engines: cover, water path, entrainment and the radiation they drive', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const on = await mixedLayerPair(2);
  console.log(`mixed-layer deck at N=6 under a 10 K inversion on a layer mixed to σ 0.85: on ${on.decked} of ${on.C} sea cells (GPU ${on.gpuDecked}), ${(1000 * on.water).toFixed(1)} g/m² where it forms; engines differ in cover by at most ${on.cover.maxDiff.toExponential(1)}, in water by ${on.mlmWater.maxDiff.toExponential(1)} kg/m² (rms ${on.mlmWater.rmsRel.toExponential(1)}), in entrainment by rms ${on.entrainment.rmsRel.toExponential(1)}; per-cell OLR rms ${on.olr.rmsRel.toExponential(1)}, surface shortwave rms ${on.sw.rmsRel.toExponential(1)}`);
  assert.ok(on.decked > 0.5 * on.C && on.gpuDecked === on.decked, `deck on ${on.decked} cells, ${on.gpuDecked} on the GPU`);
  assert.ok(on.cover.maxDiff < 1e-3 && on.fraction.maxDiff < 1e-3, `cover ${on.cover.maxDiff} at ${on.cover.at}`);
  assert.ok(on.mlmWater.rmsRel < 1e-4 && on.mlmWater.maxDiff < 5e-5 && on.deck.maxDiff < 5e-5, `water rms ${on.mlmWater.rmsRel}, max ${on.mlmWater.maxDiff} at ${on.mlmWater.at}`);
  assert.ok(on.entrainment.rmsRel < 1e-4, `entrainment rms ${on.entrainment.rmsRel}`);
  assert.ok(on.olr.rmsRel < 1e-5 && on.sw.rmsRel < 2e-5, `per-cell OLR rms ${on.olr.rmsRel}, surface shortwave rms ${on.sw.rmsRel}`);
  const dark = await mixedLayerPair(2, { stratusSolar: false });
  let thinned = 0;
  for (const i of on.sunlit) thinned += (dark.waterPath[i] - on.waterPath[i]) / dark.waterPath[i];
  console.log(`under the sun (${on.sunlit.length} decked cells lit by more than 200 W/m²) the cloud's absorption thins the step's water by ${(100 * thinned / on.sunlit.length).toFixed(1)} % on average; there the engines' water differs by rms ${on.sunlitWater.rmsRel.toExponential(1)}, at most ${on.sunlitWater.maxDiff.toExponential(1)} kg/m², the deck's by at most ${on.sunlitDeck.maxDiff.toExponential(1)} kg/m²`);
  assert.ok(on.sunlit.length > 50 && thinned > 0.005 * on.sunlit.length, `${on.sunlit.length} lit decks thinned by ${thinned / on.sunlit.length}`);
  assert.ok(on.sunlitWater.rmsRel < 1e-4 && on.sunlitWater.maxDiff < 5e-5 && on.sunlitDeck.maxDiff < 5e-5, `lit water rms ${on.sunlitWater.rmsRel}, max ${on.sunlitWater.maxDiff}; deck ${on.sunlitDeck.maxDiff}`);
  const split = await mixedLayerPair(2, { mixedLayer: { closure: 'buoyancy', decouplingOnset: 0, decoupledRatio: 0.02 } });
  console.log(`the buoyancy closure with decoupling from a buoyancy integral ratio of 0 to 0.02: ${split.partial} of ${split.decked} decks decoupled; cover differs by at most ${split.cover.maxDiff.toExponential(1)}, water by rms ${split.mlmWater.rmsRel.toExponential(1)}; OLR rms ${split.olr.rmsRel.toExponential(1)}, surface shortwave rms ${split.sw.rmsRel.toExponential(1)}`);
  assert.ok(split.partial > 0.2 * split.decked && split.gpuDecked === split.decked, `${split.partial} of ${split.decked} decoupled`);
  assert.ok(split.cover.maxDiff < 1e-3 && split.mlmWater.rmsRel < 1e-4, `cover ${split.cover.maxDiff}, water rms ${split.mlmWater.rmsRel}`);
  assert.ok(split.olr.rmsRel < 3e-5 && split.sw.rmsRel < 3e-5, `per-cell OLR rms ${split.olr.rmsRel}, surface shortwave rms ${split.sw.rmsRel}`);
});

test('the running mean of the subsidence at the boundary-layer top builds the same in both engines', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const run = await mixedLayerPair(4, { seed: 0 });
  let moved = 0;
  for (let i = 0; i < run.C; i++) if (run.mean[i] !== 0) moved++;
  console.log(`four steps from a zero mean: the mean moved on ${moved} of ${run.C} cells; engines differ by rms ${run.subsidence.rmsRel.toExponential(1)}, at most ${run.subsidence.maxDiff.toExponential(1)} m/s`);
  assert.ok(moved > 0.5 * run.C, `the mean moved on ${moved} cells`);
  assert.ok(run.subsidence.rmsRel < 1e-3, `mean subsidence rms ${run.subsidence.rmsRel}, max ${run.subsidence.maxDiff} at ${run.subsidence.at}`);
});

test('the deck reads the same ring-smoothed πσ̇ in both engines, and the smoothing moves it by far more than the engines differ', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const smoothed = await mixedLayerPair(2, { seed: 0, subsidenceMemory: 1e-9 }), raw = await mixedLayerPair(2, { seed: 0, subsidenceMemory: 1e-9, subsidenceSmoothing: 0 });
  const once = await mixedLayerPair(2, { seed: 0, subsidenceMemory: 1e-9, subsidenceSmoothing: 1 });
  const moved = stats(raw.mean, smoothed.mean), movedOnce = stats(raw.mean, once.mean), passes = stats(once.mean, smoothed.mean);
  console.log(`the step's subsidence at h, smoothed twice over the ring against unsmoothed: rms ${moved.rmsRel.toExponential(1)} relative, at most ${(1000 * moved.maxDiff).toFixed(3)} mm/s; once against unsmoothed rms ${movedOnce.rmsRel.toExponential(1)}, twice against once ${passes.rmsRel.toExponential(1)}; the engines differ by rms ${smoothed.subsidence.rmsRel.toExponential(1)} smoothed twice, ${once.subsidence.rmsRel.toExponential(1)} once and ${raw.subsidence.rmsRel.toExponential(1)} unsmoothed`);
  assert.ok(smoothed.subsidence.rmsRel < 1e-3 && once.subsidence.rmsRel < 1e-3 && raw.subsidence.rmsRel < 1e-3, `engines differ by ${smoothed.subsidence.rmsRel} smoothed twice, ${once.subsidence.rmsRel} once, ${raw.subsidence.rmsRel} unsmoothed`);
  assert.ok(moved.rmsRel > 30 * smoothed.subsidence.rmsRel, `smoothing moved the subsidence by ${moved.rmsRel}, the engines differ by ${smoothed.subsidence.rmsRel}`);
  assert.ok(movedOnce.rmsRel > 30 * once.subsidence.rmsRel && passes.rmsRel > 30 * Math.max(once.subsidence.rmsRel, smoothed.subsidence.rmsRel), `one pass moved the subsidence by ${movedOnce.rmsRel}, the second by ${passes.rmsRel}; the engines differ by ${once.subsidence.rmsRel}`);
});

test('over six steps the carried inversion height, the gate and the deck they give match between the engines, and the boundary layer mixes to the deck\'s height in both (under the gray gases: under the spectral ones a shallow cumulus fires on the fifth step in one column on the CPU alone)', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const run = await mixedLayerPair(6, { ...UNSCATTERED, longwaveScheme: 'gray', solarGases: 'lacisHansen' });
  const { C, model } = run, { geopotential, g } = model.core.diagnostics, K = model.core.K;
  let carried = 0, above = 0, active = 0;
  for (let i = 0; i < C; i++) {
    if (!(run.tops[i] > 0)) continue;
    active++;
    const depth = run.depth[i] - geopotential[(K - 1) * C + i] / g, h = run.heights[i];
    if (h > depth + 1) carried++;
    if (run.tops[i] > run.depth[i]) above++;
  }
  const theta = stats(model.state[1], run.state[1]), q = stats(model.state[4], run.state[4]);
  console.log(`six steps: the deck runs on ${active} of ${C} cells (GPU ${run.gpuDecked} with cover), its carried height stands above the boundary-layer top on ${carried} of them and hands the boundary layer a deeper profile on ${above}; engines differ in height by rms ${run.height.rmsRel.toExponential(1)} (max ${run.height.maxDiff.toExponential(1)} m), in mlmTop by max ${run.top.maxDiff.toExponential(1)} m, in the gate by max ${run.gate.maxDiff.toExponential(1)}, in cover by max ${run.cover.maxDiff.toExponential(1)}, in water by rms ${run.mlmWater.rmsRel.toExponential(1)} (max ${run.mlmWater.maxDiff.toExponential(1)} kg/m²); θ rms ${theta.rmsRel.toExponential(1)}, q rms ${q.rmsRel.toExponential(1)}`);
  assert.ok(active > 0.5 * C && carried > 0.5 * active && above > 0.5 * active, `active ${active}, carried ${carried}, deeper profile ${above}`);
  assert.ok(run.decked === run.gpuDecked, `deck on ${run.decked} cells, ${run.gpuDecked} on the GPU`);
  assert.ok(run.height.rmsRel < 1e-5 && run.top.rmsRel < 1e-5, `height rms ${run.height.rmsRel}, max ${run.height.maxDiff} at ${run.height.at}`);
  assert.ok(run.gate.maxDiff < 1e-6, `gate ${run.gate.maxDiff} at ${run.gate.at}`);
  assert.ok(run.cover.maxDiff < 1e-3 && run.fraction.maxDiff < 1e-3, `cover ${run.cover.maxDiff} at ${run.cover.at}`);
  assert.ok(run.mlmWater.rmsRel < 1e-4 && run.mlmWater.maxDiff < 5e-5 && run.deck.maxDiff < 5e-5, `water rms ${run.mlmWater.rmsRel}, max ${run.mlmWater.maxDiff} at ${run.mlmWater.at}`);
  assert.ok(run.entrainment.rmsRel < 1e-4, `entrainment rms ${run.entrainment.rmsRel}`);
  assert.ok(theta.rmsRel < 1e-5 && q.rmsRel < 5e-4, `θ rms ${theta.rmsRel}, q rms ${q.rmsRel}`);
  const capped = await mixedLayerPair(2, { height: 5000 });
  const ceiling = (i) => { let k = K - 1; while (model.core.sigmaMid[k] >= 0.85) k--; return geopotential[k * C + i] / g - 1; };
  let held = 0;
  for (let i = 0; i < C; i++) if (capped.tops[i] > 0) {
    assert.ok(Math.abs(capped.heights[i] - ceiling(i)) < 100, `cell ${i}: a 5000 m height held at ${capped.heights[i]} by the ceiling near ${ceiling(i)}`);
    held++;
  }
  console.log(`seeded at 5000 m, the ${held} decks start 1 m under the midpoint of the first layer above the 10 K inversion and end their step within 100 m of it; engines differ in height by rms ${capped.height.rmsRel.toExponential(1)}, in water by rms ${capped.mlmWater.rmsRel.toExponential(1)}`);
  assert.ok(held > 0.5 * C && capped.height.rmsRel < 1e-5 && capped.mlmWater.rmsRel < 1e-4 && capped.cover.maxDiff < 1e-3, `held ${held}, height rms ${capped.height.rmsRel}, water rms ${capped.mlmWater.rmsRel}, cover ${capped.cover.maxDiff}`);
});

test('with deckRest \'inversion\' the carried height starts and rests at the inversion ceiling alike in both engines', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const running = await mixedLayerPair(4, { deckRest: 'inversion' }), resting = await mixedLayerPair(4, { deckRest: 'inversion', seed: 5e-3, height: 300 });
  const { C, model } = resting, { geopotential, g } = model.core.diagnostics, K = model.core.K;
  let risen = 0;
  for (let i = 0; i < C; i++) if (!(resting.tops[i] > 0) && resting.heights[i] > 300 + 1) risen++;
  console.log(`four steps: unset heights start at the ceiling and the deck runs on ${running.decked} of ${C} cells (GPU ${running.gpuDecked}), height rms ${running.height.rmsRel.toExponential(1)}, gate max ${running.gate.maxDiff.toExponential(1)}, cover max ${running.cover.maxDiff.toExponential(1)}; under ascent ${risen} heights seeded at 300 m rise toward the ceiling while the deck rests, height rms ${resting.height.rmsRel.toExponential(1)}, max ${resting.height.maxDiff.toExponential(1)} m`);
  assert.ok(running.decked > 0.5 * C && running.decked === running.gpuDecked, `deck on ${running.decked}, GPU ${running.gpuDecked}`);
  assert.ok(running.height.rmsRel < 1e-5 && running.gate.maxDiff < 1e-6 && running.cover.maxDiff < 1e-3 && running.mlmWater.rmsRel < 1e-4, `height ${running.height.rmsRel}, gate ${running.gate.maxDiff}, cover ${running.cover.maxDiff}, water ${running.mlmWater.rmsRel}`);
  assert.ok(risen > 0.5 * C, `${risen} resting heights rose`);
  assert.ok(resting.height.rmsRel < 1e-5 && resting.gate.maxDiff < 1e-6, `height ${resting.height.rmsRel}, gate ${resting.gate.maxDiff}`);
  const ceiling = (i) => { let k = K - 1; while (model.core.sigmaMid[k] >= 0.85) k--; return geopotential[k * C + i] / g - 1; };
  for (let i = 0; i < C; i++) if (resting.tops[i] > 0 || resting.heights[i] > 300 + 1) assert.ok(resting.heights[i] <= ceiling(i) + 100, `cell ${i}: ${resting.heights[i]} against the ceiling near ${ceiling(i)}`);
});

test('with deckRest \'regime\' the deck stands down in surface-driven and decoupled columns whose inversion lies above cumulusCeiling, alike in both engines', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const shared = { deckRest: 'regime', turbulence: 'moist', longwaveScheme: 'gray', solarGases: 'lacisHansen', moist: { cloudLifetime: 3 * 3600, plumeCape: 70, condensation: 'saturation', iceFall: null } };
  const high = await mixedLayerPair(4, { ...shared, cumulusCeiling: 3000 }), low = await mixedLayerPair(4, { ...shared, cumulusCeiling: 200 });
  const regimes = (run) => { const n = [0, 0, 0, 0]; for (const r of run.model.boundaryLayer.regime) n[r]++; return n; };
  let shut = 0;
  for (let i = 0; i < low.C; i++) if (low.gpuGate[i] === 0 && low.model.radiation.mlmGate[i] === 0) shut++;
  console.log(`four moist steps under a 10 K inversion, regimes (stable, surface, decoupled, coupled) ${regimes(high).join(', ')}: with the stand-down above 3 km decks on ${high.decked} of ${high.C} cells (GPU ${high.gpuDecked}), above 200 m on ${low.decked} (GPU ${low.gpuDecked}) with ${shut} gates shut in both; height rms ${high.height.rmsRel.toExponential(1)} and ${low.height.rmsRel.toExponential(1)}, gate max ${high.gate.maxDiff.toExponential(1)} and ${low.gate.maxDiff.toExponential(1)}, water rms ${high.mlmWater.rmsRel.toExponential(1)} and ${low.mlmWater.rmsRel.toExponential(1)}`);
  assert.ok(high.decked > 0.5 * high.C && low.decked < 0.5 * high.decked && shut > 0.5 * low.C, `decks ${high.decked} and ${low.decked}, ${shut} shut`);
  for (const run of [high, low]) {
    assert.ok(run.decked === run.gpuDecked, `deck on ${run.decked}, GPU ${run.gpuDecked}`);
    assert.ok(run.height.rmsRel < 1e-5 && run.gate.maxDiff < 1e-6 && run.cover.maxDiff < 1e-3 && run.mlmWater.rmsRel < 1e-3, `height ${run.height.rmsRel}, gate ${run.gate.maxDiff}, cover ${run.cover.maxDiff}, water ${run.mlmWater.rmsRel}`);
  }
});

test('a ceilingInversion below minimumInversion holds the deck under a weaker jump, alike in both engines', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const shared = { deckRest: 'inversion', minimumInversion: 4, step: 0.92 };
  const strong = await mixedLayerPair(4, shared), weak = await mixedLayerPair(4, { ...shared, ceilingInversion: 2 });
  let lower = 0;
  for (let i = 0; i < weak.C; i++) if (weak.heights[i] < strong.heights[i] - 50) lower++;
  console.log(`four steps under a 3 K step below a 10 K inversion: with a 2 K ceiling the carried height sits lower on ${lower} of ${weak.C} cells (mean ${(weak.heights.reduce((a, b) => a + b) / weak.C).toFixed(0)} m against ${(strong.heights.reduce((a, b) => a + b) / strong.C).toFixed(0)} m), decks on ${weak.decked} and ${strong.decked}; height rms ${weak.height.rmsRel.toExponential(1)} and ${strong.height.rmsRel.toExponential(1)}, gate max ${weak.gate.maxDiff.toExponential(1)} and ${strong.gate.maxDiff.toExponential(1)}`);
  assert.ok(lower > 0.5 * weak.C, `${lower} heights lower`);
  for (const run of [weak, strong]) {
    assert.ok(run.decked === run.gpuDecked, `deck on ${run.decked}, GPU ${run.gpuDecked}`);
    assert.ok(run.height.rmsRel < 1e-5 && run.gate.maxDiff < 1e-6, `height ${run.height.rmsRel}, gate ${run.gate.maxDiff}`);
  }
});

test('the GPU model sends the deck\'s running-mean subsidence, carried height and gate and the boundary layer\'s depth, mixing top, regime and surface buoyancy flux to the device on load and reads them back on sync', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { createGpuModel } = await import('../js/gpu/model.gpu.js');
  const model = await createGpuModel(new Grid(6), { ocean: false, radiation: { mixedLayerDeck: true } });
  const C = model.mesh.nCells, init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  const loaded = Float64Array.from({ length: C }, (_, i) => -2e-3 * ((i % 7) + 1) / 7);
  const heights = Float64Array.from({ length: C }, (_, i) => (i % 3 ? 400 + 37 * (i % 11) : 0)), gates = Float64Array.from({ length: C }, (_, i) => (i % 5) / 4);
  model.radiation.mlmSubsidence.set(loaded);
  model.radiation.mlmHeight.set(heights);
  model.radiation.mlmGate.set(gates);
  const depth = Float64Array.from({ length: C }, (_, i) => 300 + 23 * (i % 13)), top = Float64Array.from(depth, (d, i) => d + 50 * (i % 4)), regimes = Float64Array.from({ length: C }, (_, i) => i % 4), buoyancy = Float64Array.from({ length: C }, (_, i) => 1e-4 * ((i % 5) - 2));
  model.boundaryLayer.depth.set(depth);
  model.boundaryLayer.mixingTop.set(top);
  model.boundaryLayer.regime.set(regimes);
  model.boundaryLayer.buoyancyFlux.set(buoyancy);
  model.load();
  const device = await model.gpu.downloadPhysics();
  for (let i = 0; i < C; i++) {
    assert.equal(device.MLMSUB[i], Math.fround(loaded[i]), `cell ${i} on the device`);
    assert.ok(device.MLMH[i] === Math.fround(heights[i]) && device.MLMGATE[i] === Math.fround(gates[i]), `cell ${i}: height and gate on the device`);
    assert.equal(device.DEPTH[i], Math.fround(depth[i]), `cell ${i}: boundary-layer depth on the device`);
    assert.ok(device.MIXTOP[i] === Math.fround(top[i]) && device.REGIME[i] === regimes[i], `cell ${i}: mixing top and regime on the device`);
    assert.equal(device.BUOY[i], Math.fround(buoyancy[i]), `cell ${i}: surface buoyancy flux on the device`);
  }
  await model.step(900); await model.step(900);
  await model.sync();
  const after = await model.gpu.downloadPhysics(), mean = model.radiation.mlmSubsidence;
  let moved = 0;
  for (let i = 0; i < C; i++) {
    assert.ok(model.radiation.mlmHeight[i] === after.MLMH[i] && model.radiation.mlmGate[i] === after.MLMGATE[i], `cell ${i}: height and gate mirrored`);
    assert.equal(model.boundaryLayer.depth[i], after.DEPTH[i], `cell ${i}: boundary-layer depth mirrored`);
    assert.ok(model.boundaryLayer.mixingTop[i] === after.MIXTOP[i] && model.boundaryLayer.regime[i] === after.REGIME[i], `cell ${i}: mixing top and regime mirrored`);
    assert.equal(model.boundaryLayer.buoyancyFlux[i], after.BUOY[i], `cell ${i}: surface buoyancy flux mirrored`);
    assert.equal(mean[i], after.MLMSUB[i], `cell ${i} mirrored`);
    assert.ok(Math.abs(mean[i] - loaded[i]) < 1e-4, `cell ${i}: ${mean[i]} from ${loaded[i]}`);
    if (mean[i] !== Math.fround(loaded[i])) moved++;
  }
  assert.ok(moved > C / 2, `the mean moved on ${moved} cells`);
});

test('the convective and large-scale rain accumulate alike in both engines, cell by cell but for the odd column whose onset falls a step apart (a plume shortens the lifetime of the cloud below its top), and add up to the precipitation', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { model, physics } = await pair(6, 24, 900, 0, false, {}, { rainEvaporation: 0 });
  const C = model.mesh.nCells, { convectivePrecipitation: convective, largeScalePrecipitation: largeScale, precipitation, rain } = model.moist;
  const largest = Math.max(...convective), largestScale = Math.max(...largeScale);
  const onset = (i) => Math.abs(convective[i] - physics.CONV[i]) > 1e-3 * largest || Math.abs(largeScale[i] - physics.COND[i]) > 1e-3 * largestScale;
  const kept = Array.from({ length: C }, (_, i) => i).filter((i) => !onset(i)), pick = (values) => Float64Array.from(kept, (i) => values[i]);
  const conv = stats(pick(convective), pick(physics.CONV)), ls = stats(pick(largeScale), pick(physics.COND)), step = stats(pick(rain), pick(physics.STEPRAIN));
  let fired = 0, rained = 0, apart = 0, area = 0, cpuMean = 0, gpuMean = 0;
  for (let i = 0; i < C; i++) {
    const a = model.mesh.areaCell[i];
    if (convective[i] > 0) fired++;
    if (largeScale[i] > 0) rained++;
    apart = Math.max(apart, Math.abs(convective[i] + largeScale[i] - precipitation[i]));
    area += a; cpuMean += a * convective[i]; gpuMean += a * physics.CONV[i];
  }
  console.log(`24 steps at N=6: convective rain on ${fired} of ${C} cells, ${(cpuMean / area).toFixed(4)} kg/m² in the mean (GPU ${(gpuMean / area).toFixed(4)}); ${C - kept.length} cells apart by more than 1e-3 of the largest cell's rain of either kind; over the rest per-cell rms ${conv.rmsRel.toExponential(1)}, large-scale on ${rained}, per-cell rms ${ls.rmsRel.toExponential(1)}; the last step's rain differs by at most ${step.maxDiff.toExponential(1)} kg/m²`);
  assert.ok(fired > C / 2 && rained > 0, `convective rain on ${fired} cells, large-scale on ${rained}`);
  assert.ok(C - kept.length <= 0.02 * C, `${C - kept.length} cells apart`);
  assert.ok(Math.abs(gpuMean - cpuMean) < 1e-3 * cpuMean, `mean convective rain ${cpuMean / area} against ${gpuMean / area}`);
  assert.ok(conv.rmsRel < 1e-3 && ls.rmsRel < 1e-3, `per-cell rms convective ${conv.rmsRel}, large-scale ${ls.rmsRel}`);
  assert.ok(step.maxDiff < 3e-4, `the last step's rain differs by ${step.maxDiff} at ${step.at}`);
  assert.ok(apart < 1e-12, `convective plus large-scale is the precipitation to ${apart}`);
});

test('both models read the rain split out at the diagnostics as means in mm/d, clearing its sums, alike but for the odd column whose onset falls a step apart, and the GPU model sends the means to the device on load and mirrors them on sync', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { createGpuModel } = await import('../js/gpu/model.gpu.js');
  const cpu = createModel(new Grid(6), { ocean: false }), gpu = await createGpuModel(new Grid(6), { ocean: false });
  const C = cpu.mesh.nCells, init = initializeState(cpu, {});
  for (let a = 0; a < init.length; a++) { cpu.state[a].set(init[a]); gpu.state[a].set(init[a]); }
  const loaded = Float64Array.from({ length: C }, (_, i) => 0.5 * (i % 9));
  gpu.moist.convectiveRain.set(loaded);
  gpu.moist.largeScaleRain.set(loaded.map((x) => x / 3));
  gpu.load();
  let device = await gpu.gpu.downloadPhysics();
  for (let i = 0; i < C; i++) assert.ok(device.CONVMEAN[i] === Math.fround(loaded[i]) && device.CONDMEAN[i] === Math.fround(loaded[i] / 3), `cell ${i}: the means on the device`);
  await gpu.diagnostics();
  device = await gpu.gpu.downloadPhysics();
  assert.equal(device.CONVMEAN[1], Math.fround(loaded[1]), 'a diagnostics frame with no time elapsed keeps the means');
  const steps = 24, seconds = steps * 900;
  for (let n = 0; n < steps; n++) { cpu.step(900); await gpu.step(900); }
  const sums = Float64Array.from(cpu.moist.convectivePrecipitation);
  device = await gpu.gpu.downloadPhysics();
  const gpuSums = Float64Array.from(device.CONV.subarray(0, C)), gpuLarge = Float64Array.from(device.COND.subarray(0, C));
  cpu.diagnostics();
  await gpu.diagnostics();
  await gpu.sync();
  device = await gpu.gpu.downloadPhysics();
  let area = 0, cpuMean = 0, gpuMean = 0;
  for (let i = 0; i < C; i++) {
    assert.ok(Math.abs(cpu.moist.convectiveRain[i] - 86400 * sums[i] / seconds) < 1e-12, `cell ${i}: the CPU mean is the sum over the interval`);
    assert.ok(Math.abs(gpu.moist.convectiveRain[i] - 86400 * gpuSums[i] / seconds) <= 1e-6 * gpu.moist.convectiveRain[i] && Math.abs(gpu.moist.largeScaleRain[i] - 86400 * gpuLarge[i] / seconds) <= 1e-6 * gpu.moist.largeScaleRain[i], `cell ${i}: the GPU mean is its sum over the interval`);
    assert.ok(gpu.moist.convectiveRain[i] === device.CONVMEAN[i] && gpu.moist.largeScaleRain[i] === device.CONDMEAN[i], `cell ${i} mirrored`);
    assert.ok(cpu.moist.convectivePrecipitation[i] === 0 && cpu.moist.largeScalePrecipitation[i] === 0 && device.CONV[i] === 0 && device.COND[i] === 0, `cell ${i}: the sums start again`);
    const a = cpu.mesh.areaCell[i];
    area += a; cpuMean += a * cpu.moist.convectiveRain[i]; gpuMean += a * gpu.moist.convectiveRain[i];
  }
  const largest = Math.max(...cpu.moist.convectiveRain), kept = Array.from({ length: C }, (_, i) => i).filter((i) => !(Math.abs(cpu.moist.convectiveRain[i] - gpu.moist.convectiveRain[i]) > 1e-3 * largest));
  const pick = (values) => Float64Array.from(kept, (i) => values[i]);
  const convective = stats(pick(cpu.moist.convectiveRain), pick(gpu.moist.convectiveRain)), largeScale = stats(pick(cpu.moist.largeScaleRain), pick(gpu.moist.largeScaleRain));
  console.log(`six hours at N=6: convective rain ${(cpuMean / area).toFixed(3)} mm/d in the mean (GPU ${(gpuMean / area).toFixed(3)}); ${C - kept.length} cells apart by more than 1e-3 of the largest cell's rain; over the rest per-cell rms ${convective.rmsRel.toExponential(1)}, large-scale per-cell rms ${largeScale.rmsRel.toExponential(1)}`);
  assert.ok(cpuMean / area > 0.1 && Math.abs(gpuMean - cpuMean) < 1e-3 * cpuMean, `mean convective rain ${cpuMean / area} against ${gpuMean / area} mm/d`);
  assert.ok(C - kept.length <= 0.01 * C, `${C - kept.length} cells apart`);
  assert.ok(convective.rmsRel < 1e-3 && largeScale.rmsRel < 1e-3, `per-cell rms convective ${convective.rmsRel}, large-scale ${largeScale.rmsRel}`);
});

test('step by step from one state with partly iced, melting polar cells, the engines agree on the ice, its snow and their albedo, the surface flux, the lowest layers and the boundary layer\'s regime', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { createGpuModel } = await import('../js/gpu/model.gpu.js');
  const { sigmaInterfaces } = await import('../js/dynamics/sigmaCore.module.js');
  const { readRanges } = await import('../js/gpu/device.module.js');
  const N = 6, dt = 1350 * 16 / N, options = { ocean: false, levels: sigmaInterfaces('bl34'), radiation: { clearSkyPass: true } };
  const cpu = createModel(new Grid(N), options), gpu = await createGpuModel(new Grid(N), options);
  const C = cpu.mesh.nCells, K = cpu.core.K, init = initializeState(cpu, {});
  const iced = [];
  for (let i = 0; i < C; i++) if (cpu.mesh.latCell[i] > 70 * Math.PI / 180) { iced.push(i); init[6][i] = 1.2; init[3][i] = 273.1; }
  for (const a of init) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  const concentration = Float64Array.from({ length: C }, (_, i) => (iced.includes(i) ? Math.fround(0.7) : 0)), snow = Float32Array.from({ length: C }, (_, i) => (iced.includes(i) ? 3 : 0));
  for (const m of [cpu, gpu]) { for (let a = 0; a < init.length; a++) m.state[a].set(init[a]); m.seaIce.load(m.state[6], concentration); m.time = 91 * 86400 + 12 * 3600; }
  cpu.seaIce.snow.set(snow);
  gpu.load();
  const { device, buffers, layout } = gpu.gpu;
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.SNOW, snow);
  const worst = { h: 0, A: 0, snow: 0, albedo: 0, skin: 0, flux: 0, iceFlux: 0, theta: 0, flips: 0 };
  let melted = 0;
  for (let n = 0; n < 8; n++) {
    cpu.step(dt);
    await gpu.step(dt);
    await gpu.sync();
    const [conc, snowNow, snowAlbedo, flux, regime] = await readRanges(device, buffers.PH, ['CONC', 'SNOW', 'SNOWALB', 'SFLUX', 'REGIME'].map((name) => ({ offset: layout.PH[name], length: C })));
    for (let i = 0; i < C; i++) {
      if (cpu.boundaryLayer.regime[i] !== regime[i]) worst.flips++;
      worst.flux = Math.max(worst.flux, Math.abs(cpu.radiation.surfaceFlux[i] - flux[i]));
      for (let k = K - 10; k < K; k++) worst.theta = Math.max(worst.theta, Math.abs(cpu.state[1][k * C + i] - gpu.state[1][k * C + i]));
    }
    for (const i of iced) {
      worst.h = Math.max(worst.h, Math.abs(cpu.state[6][i] - gpu.state[6][i])); worst.A = Math.max(worst.A, Math.abs(cpu.seaIce.concentration[i] - conc[i]));
      worst.snow = Math.max(worst.snow, Math.abs(cpu.seaIce.snow[i] - snowNow[i])); worst.albedo = Math.max(worst.albedo, Math.abs(cpu.seaIce.snowAlbedo[i] - snowAlbedo[i]));
      worst.skin = Math.max(worst.skin, Math.abs(cpu.state[3][i] - gpu.state[3][i])); worst.iceFlux = Math.max(worst.iceFlux, Math.abs(cpu.radiation.surfaceFlux[i] - flux[i]));
      if (n === 7 && cpu.state[3][i] >= 273.15 && cpu.seaIce.snow[i] < 3) melted++;
    }
  }
  console.log(`8 steps, ${iced.length} cells iced at 0.7 under 1.2 m and 3 kg/m² of snow (${melted} melting at the end): the engines differ by up to ${worst.h.toExponential(1)} m of ice, ${worst.A.toExponential(1)} of cover, ${worst.snow.toExponential(1)} kg/m² of snow, ${worst.albedo.toExponential(1)} of its albedo, ${worst.skin.toExponential(1)} K at the skin and ${worst.iceFlux.toExponential(1)} W/m² of surface flux there; everywhere ${worst.flux.toExponential(1)} W/m², ${worst.theta.toExponential(1)} K of θ in the lowest ten layers, ${worst.flips} boundary-layer regimes`);
  assert.ok(melted > 0, 'some cell melts');
  assert.ok(worst.h < 1e-4 && worst.A < 1e-4 && worst.snow < 1e-3 && worst.albedo < 1e-4 && worst.skin < 0.01 && worst.iceFlux < 0.5, JSON.stringify(worst));
  assert.ok(worst.flux < 1 && worst.theta < 0.01 && worst.flips === 0, JSON.stringify(worst));
});
