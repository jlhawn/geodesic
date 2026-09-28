import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { sunDirection } from '../js/physics/radiation.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuCore } = gpuAvailable ? await import('../js/gpu/core.gpu.js') : {};

function meanTheta(model) {
  const { K } = model.core, C = model.mesh.nCells, theta = model.state[1];
  return Float64Array.from({ length: K }, (_, k) => { let s = 0, a = 0; for (let i = 0; i < C; i++) { s += model.mesh.areaCell[i] * theta[k * C + i]; a += model.mesh.areaCell[i]; } return s / a; });
}
function stats(cpu, gpu) {
  let maxDiff = 0, sumSq = 0, sumRef = 0, at = -1;
  for (let x = 0; x < cpu.length; x++) { const d = Math.abs(cpu[x] - gpu[x]); if (d > maxDiff) { maxDiff = d; at = x; } sumSq += d * d; sumRef += cpu[x] * cpu[x]; }
  return { maxDiff, at, rms: Math.sqrt(sumSq / cpu.length), rmsRel: Math.sqrt(sumSq / Math.max(sumRef, 1e-300)) };
}
async function pair(N, steps, dt, inversion = 0, stratus = inversion > 0, options = {}) {
  const model = createModel(new Grid(N), { ocean: false, radiation: { stratus, ...options } });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  const { K, sigmaMid } = model.core, C = model.mesh.nCells;
  for (let k = 0; k < K; k++) if (sigmaMid[k] < 0.75) for (let i = 0; i < C; i++) model.state[1][k * C + i] += inversion;
  const gpu = await createGpuCore(model.mesh, { nu4: model.core.nu4, nu4Theta: model.core.nu4Theta, referenceTheta: meanTheta(model), physics: { stratus, ...options } });
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
  assert.ok(olr.rmsRel < 1e-5 && sw.rmsRel < 1e-5, `per-cell OLR rms ${olr.rmsRel}, surface shortwave rms ${sw.rmsRel}`);
  const third = await pair(6, 2, 900, 10, true, { stratusIndex: 'ectei', mixedLayerDeck: false }), entrained = third.model.radiation;
  const olrE = stats(entrained.outgoing, third.physics.OLR.subarray(0, C)), swE = stats(entrained.surfaceShortwave, third.physics.SWDN.subarray(0, C)), coverE = stats(entrained.stratusFraction, third.physics.DECKF.subarray(0, C));
  let coverSum = 0, waterSum = 0;
  for (let i = 0; i < C; i++) { coverSum += entrained.stratusFraction[i]; waterSum += entrained.stratus[i]; }
  console.log(`the same under ECTEI: mean cover ${(coverSum / C).toFixed(3)} and mean water ${(1000 * waterSum / C).toFixed(1)} g/m²; per-cell cover max difference ${coverE.maxDiff.toExponential(1)}, OLR rms ${olrE.rmsRel.toExponential(1)}, surface shortwave rms ${swE.rmsRel.toExponential(1)}`);
  assert.ok(coverSum > 0 && coverSum < cover, `mean cover ${coverSum / C} under ECTEI against ${cover / C}`);
  assert.ok(olrE.rmsRel < 1e-5 && swE.rmsRel < 1e-5, `per-cell OLR rms ${olrE.rmsRel}, surface shortwave rms ${swE.rmsRel} under ECTEI`);
});

/*
 * One physics kernel alone against the CPU physics phase, both from the
 * same single-precision state (two steps after a 10 K inversion was
 * imposed): the heating of every layer, (θ' − θ)Π/dt, over a step long
 * enough that the rounding of θ is far below the tolerance.
 */
function heatingState() {
  const model = createModel(new Grid(6), { ocean: false, radiation: { stratus: true, mixedLayerDeck: false } });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  const { K, sigmaMid } = model.core, C = model.mesh.nCells;
  for (let k = 0; k < K; k++) if (sigmaMid[k] < 0.75) for (let i = 0; i < C; i++) model.state[1][k * C + i] += 10;
  for (let n = 0; n < 2; n++) model.step(900);
  for (const a of model.state) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  return model;
}
async function physicsHeating(base, options, dt = 864000) {
  const physics = { stratus: true, mixedLayerDeck: false, ...options };
  const model = createModel(base.mesh, { ocean: false, radiation: physics });
  model.state.forEach((a, n) => a.set(base.state[n]));
  model.seaIce.concentration.set(base.seaIce.concentration);
  model.seaIce.snow.set(base.seaIce.snow);
  model.boundaryLayer.depth.set(base.boundaryLayer.depth);
  model.time = base.time;
  const { K } = model.core, C = model.mesh.nCells;
  const gpu = await createGpuCore(model.mesh, { nu4: model.core.nu4, nu4Theta: model.core.nu4Theta, referenceTheta: meanTheta(model), physics });
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
  return { K, C, cloudy, cpu: heating(model.state[1]), gpu: heating(after), power, area: model.mesh.areaCell, cpuDeck: Float64Array.from(model.radiation.stratusFraction), gpuDeck: ph.DECKF.subarray(0, C) };
}

test('the heating of each layer of the sunlit cloudy columns, and the part of it the cloud water absorbs, agree between the engines', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const base = heatingState();
  const lit = await physicsHeating(base, {}), scattering = await physicsHeating(base, { cloudSolarAbsorption: 0 });
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

test('twelve full GPU steps track the CPU model and its energy budget', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { model, state, physics } = await pair(6, 12, 900);
  const d = model.diagnostics();
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
  const gpu = await createGpuCore(model.mesh, { nu4: model.core.nu4, nu4Theta: model.core.nu4Theta, referenceTheta: meanTheta(model) });
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
  const gpu = await createGpuCore(model.mesh, { nu4: model.core.nu4, nu4Theta: model.core.nu4Theta, referenceTheta: meanTheta(model) });
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
async function mixedLayerPair(steps, { seed = -1e-3, ...options } = {}) {
  const physics = { mixedLayerDeck: true, ...options };
  const model = createModel(new Grid(6), { ocean: false, radiation: physics });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  const { K, sigmaMid } = model.core, C = model.mesh.nCells, theta = model.state[1], q = model.state[4];
  for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) {
    if (sigmaMid[k] < 0.85) theta[k * C + i] += 10;
    else { theta[k * C + i] = theta[(K - 1) * C + i]; q[k * C + i] = q[(K - 1) * C + i]; }
  }
  model.radiation.mlmSubsidence.fill(seed);
  const gpu = await createGpuCore(model.mesh, { nu4: model.core.nu4, nu4Theta: model.core.nu4Theta, referenceTheta: meanTheta(model), physics });
  gpu.upload(model.state);
  gpu.uploadPhysics({ mlmSubsidence: model.radiation.mlmSubsidence });
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
  assert.ok(on.olr.rmsRel < 1e-5 && on.sw.rmsRel < 1e-5, `per-cell OLR rms ${on.olr.rmsRel}, surface shortwave rms ${on.sw.rmsRel}`);
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

test('the GPU model sends the running-mean subsidence to the device on load and reads it back on sync', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { createGpuModel } = await import('../js/gpu/model.gpu.js');
  const model = await createGpuModel(new Grid(6), { ocean: false, radiation: { mixedLayerDeck: true } });
  const C = model.mesh.nCells, init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  const loaded = Float64Array.from({ length: C }, (_, i) => -2e-3 * ((i % 7) + 1) / 7);
  model.radiation.mlmSubsidence.set(loaded);
  model.load();
  const device = await model.gpu.downloadPhysics();
  for (let i = 0; i < C; i++) assert.equal(device.MLMSUB[i], Math.fround(loaded[i]), `cell ${i} on the device`);
  await model.step(900); await model.step(900);
  await model.sync();
  const after = await model.gpu.downloadPhysics(), mean = model.radiation.mlmSubsidence;
  let moved = 0;
  for (let i = 0; i < C; i++) {
    assert.equal(mean[i], after.MLMSUB[i], `cell ${i} mirrored`);
    assert.ok(Math.abs(mean[i] - loaded[i]) < 1e-4, `cell ${i}: ${mean[i]} from ${loaded[i]}`);
    if (mean[i] !== Math.fround(loaded[i])) moved++;
  }
  assert.ok(moved > C / 2, `the mean moved on ${moved} cells`);
});
