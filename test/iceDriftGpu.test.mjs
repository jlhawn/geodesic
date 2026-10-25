import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createOcean } from '../js/ocean/layered.module.js';
import { iceDriftMoves } from '../js/ocean/iceDrift.module.js';
import { createSeaIce, FREEZING_POINT, MINIMUM_CONCENTRATION } from '../js/physics/ice.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { DEG, streamOn, northwardOn } from './helpers/iceDrift.mjs';

/*
 * The GPU ice drift (js/gpu/iceDrift.gpu.js, merged into the hybrid ocean's
 * kernels) against the serial CPU ocean with the same drift at N = 16 with
 * everySteps 1 and DT = 1350 s, on a sea with a continent: the transport of
 * a prescribed velocity, one momentum call each of free drift and the
 * cavitating fluid, twenty coupled calls, and the single-precision
 * conservation of the volume and snow mass over a hundred calls. The wind
 * stress is handed to both engines as the same array. One GPU model at a
 * time (CLAUDE.md).
 */
let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};
const { readRanges } = gpuAvailable ? await import('../js/gpu/device.module.js') : {};

const N = 16, DT = 1350, SKIN = 2e5;
const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 300 : -4000));

function stats(reference, other) {
  let sumSq = 0, sumRef = 0, maxDiff = 0;
  for (let x = 0; x < reference.length; x++) { const d = reference[x] - other[x]; sumSq += d * d; sumRef += reference[x] * reference[x]; maxDiff = Math.max(maxDiff, Math.abs(d)); }
  return { rmsRel: Math.sqrt(sumSq / Math.max(sumRef, 1e-300)), maxDiff };
}
function windStress(mesh) {
  const stress = new Float64Array(mesh.nEdges);
  for (let e = 0; e < mesh.nEdges; e++) {
    const lon = Math.atan2(mesh.xEdge[3 * e + 1], mesh.xEdge[3 * e]);
    stress[e] = 0.08 * Math.sin(2 * mesh.latEdge[e]) * Math.cos(lon) + 0.01 * Math.sin(5 * lon);
  }
  return stress;
}

/*
 * Both engines from the same state: ice of cover 0.85 poleward of 60S with a
 * fringe of thin ice toward 55S (the whole sea at cover 0.8 with a bump when
 * `global`, with a bump of thick ice around 78S), snow on it and a cold skin,
 * the ocean initialized on the CPU and uploaded to the device with the ice in
 * the atmosphere's buffers.
 */
async function scenario(drift, { global = false, velocity = null, stress = null, keep = iceDriftMoves(drift), ocean: options = {} } = {}) {
  const gpuModel = await createGpuModel(new Grid(N), { topography, ice: { keepSmallIce: keep }, ocean: { vertical: 'hybrid', everySteps: 1, iceDrift: drift ? { ...drift, airStress: 'given' } : false, ...options } });
  const mesh = gpuModel.mesh, geography = gpuModel.geography, C = mesh.nCells, E = mesh.nEdges, gpu = gpuModel.gpu, ocean = gpuModel.oceanEngine;
  const cpu = createOcean(mesh, { geography, vertical: 'hybrid', everySteps: 1, iceDrift: drift, ...options });
  const seaIce = createSeaIce(mesh, { heatCapacity: cpu.capacity, keepSmallIce: keep });
  cpu.useSeaIce(seaIce.concentration, seaIce);
  const surfaceT = new Float64Array(C).fill(280), ice = new Float64Array(C);
  const centre = [Math.cos(-78 * DEG) * Math.cos(40 * DEG), Math.cos(-78 * DEG) * Math.sin(40 * DEG), Math.sin(-78 * DEG)];
  for (let i = 0; i < C; i++) {
    if (!cpu.cellOcean[i]) continue;
    const lat = mesh.latCell[i] / DEG;
    const bump = Math.exp(-((Math.acos(Math.min(1, mesh.xCell[3 * i] * centre[0] + mesh.xCell[3 * i + 1] * centre[1] + mesh.xCell[3 * i + 2] * centre[2])) * mesh.radius / 1.5e6) ** 2));
    if (global) { seaIce.concentration[i] = 0.4 + 0.5 * bump; ice[i] = 0.3 + 0.8 * bump; seaIce.snow[i] = 20 * bump; surfaceT[i] = FREEZING_POINT - 6 * bump; continue; }
    if (lat < -60) { seaIce.concentration[i] = 0.85; ice[i] = 0.6; seaIce.snow[i] = 12; surfaceT[i] = FREEZING_POINT - 4; } else if (lat < -55) { seaIce.concentration[i] = 0.3; ice[i] = 0.3; seaIce.snow[i] = 2; surfaceT[i] = FREEZING_POINT - 1; }
  }
  cpu.initialize(surfaceT, ice, { concentration: seaIce.concentration });
  for (let i = 0; i < C; i++) if (cpu.cellOcean[i] && !(ice[i] > 0)) surfaceT[i] = cpu.T0[i];
  seaIce.snowAlbedo.fill(0.85);
  if (velocity) cpu.iceDrift.velocity.set(velocity);
  const wind = stress ?? new Float64Array(E);
  const volume = Float64Array.from(ice, (h, i) => h * seaIce.concentration[i]);
  await ocean.upload({ ...cpu.serialize(), Q: cpu.Q, W: cpu.W, previousT0: cpu.T0, previousIceVolume: volume, fresh: new Float64Array(C), buoyancy: new Float64Array(C), flux: new Float64Array(C), capacity: cpu.capacity, ...(velocity ? { iceVelocity: velocity } : {}) }, Float64Array.from(surfaceT), ice, seaIce.concentration);
  gpu.uploadSurfaceTemperature(Float64Array.from(surfaceT));
  const write = (buffer, offset, values) => gpu.device.queue.writeBuffer(buffer, 4 * offset, Float32Array.from(values));
  write(gpu.buffers.S, gpu.layout.S.ICE, ice);
  write(gpu.buffers.PH, gpu.layout.PH.CONC, seaIce.concentration); write(gpu.buffers.PH, gpu.layout.PH.SNOW, seaIce.snow); write(gpu.buffers.PH, gpu.layout.PH.SNOWALB, seaIce.snowAlbedo);
  if (drift) ocean.setAirStress(wind);
  const s = { gpuModel, mesh, C, E, gpu, ocean, cpu, seaIce, surfaceT, ice, wind, flux: new Float64Array(C) };
  s.stepCpu = () => cpu.advance(surfaceT, ice, seaIce.oceanFlux, wind, DT, seaIce.concentration);
  s.stepGpu = async () => { await ocean.advanceCoupled(DT); await gpu.device.queue.onSubmittedWorkDone(); };
  s.readIce = async () => {
    const [h, ts] = await readRanges(gpu.device, gpu.buffers.S, [{ offset: gpu.layout.S.ICE, length: C }, { offset: gpu.layout.S.TS, length: C }]);
    const [conc, snow, albedo] = await readRanges(gpu.device, gpu.buffers.PH, ['CONC', 'SNOW', 'SNOWALB'].map((name) => ({ offset: gpu.layout.PH[name], length: C })));
    return { h, ts, conc, snow, albedo };
  };
  return s;
}
/* Area, volume, snow mass and skin enthalpy per cell area, on the CPU and from the device. */
const fieldsOf = (h, conc, snow, ts, open) => ({
  A: Float64Array.from(conc, (a, i) => (h[i] > 0 ? a : 0)), V: Float64Array.from(h, (x, i) => (x > 0 ? x * conc[i] : 0)), S: Float64Array.from(snow, (x, i) => (h[i] > 0 ? x * conc[i] : 0)),
  E: Float64Array.from(ts, (t, i) => (h[i] > 0 && open[i] ? SKIN * conc[i] * (t - FREEZING_POINT) : 0)),
});
const cpuFields = (s) => fieldsOf(s.ice, s.seaIce.concentration, s.seaIce.snow, s.surfaceT, s.cpu.cellOcean);
const gpuFields = async (s) => { const g = await s.readIce(); return fieldsOf(g.h, g.conc, g.snow, g.ts, s.cpu.cellOcean); };
const compareIce = async (label, s, limit) => {
  const a = cpuFields(s), b = await gpuFields(s), out = {};
  for (const key of ['A', 'V', 'S', 'E']) out[key] = stats(a[key], b[key]);
  console.log(`${label}: ${Object.entries(out).map(([key, r]) => `${key} rms ${r.rmsRel.toExponential(2)} max ${r.maxDiff.toExponential(2)}`).join(', ')}`);
  for (const key of ['A', 'V', 'S', 'E']) assert.ok(out[key].rmsRel < (limit[key] ?? limit.any), `${label}: ${key} rms relative ${out[key].rmsRel}`);
  return out;
};
const compareVelocity = async (label, s, limit) => {
  const { iceVelocity } = await s.ocean.restartArrays(), r = stats(s.cpu.iceDrift.velocity, iceVelocity);
  let top = 0; for (const u of s.cpu.iceDrift.velocity) top = Math.max(top, Math.abs(u));
  console.log(`${label}: ice velocity rms ${r.rmsRel.toExponential(2)} max ${r.maxDiff.toExponential(2)} of ${top.toFixed(4)} m/s`);
  assert.ok(top > 0.01 && r.rmsRel < limit, `${label}: velocity rms relative ${r.rmsRel}`);
  return r;
};

test('the GPU carries a prescribed ice velocity like the CPU: area, volume, snow and skin enthalpy after twenty calls', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const probe = await createGpuModel(new Grid(N), { topography, ocean: { vertical: 'hybrid', everySteps: 1 } }), velocity = streamOn(probe.mesh, 20 / probe.mesh.radius);
  probe.destroy();
  const s = await scenario({ rheology: 'frozen' }, { velocity });
  for (let n = 0; n < 20; n++) { s.stepCpu(); await s.stepGpu(); }
  await compareIce('transport of a prescribed velocity', s, { A: 2e-5, V: 2e-5, S: 2e-5, E: 1e-4 });
  const { iceVelocity } = await s.ocean.restartArrays();
  assert.ok(stats(velocity, iceVelocity).rmsRel < 1e-6, 'the frozen drift keeps the velocity it holds');
  s.gpuModel.destroy();
});

for (const [rheology, limit] of [['free', 2e-3], ['cavitating', 1e-2]]) {
  test(`one ${rheology} momentum call: the GPU's ice velocity follows the CPU's`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
    const probe = await createGpuModel(new Grid(N), { topography, ocean: { vertical: 'hybrid', everySteps: 1 } }), stress = windStress(probe.mesh);
    probe.destroy();
    const s = await scenario({ rheology }, { stress });
    s.stepCpu(); await s.stepGpu();
    await compareVelocity(`${rheology} drift after one call`, s, limit);
    await compareIce(`${rheology} drift after one call`, s, { any: 1e-3 });
    s.gpuModel.destroy();
  });
}

test('the hemisphere switch keeps the other hemisphere\'s ice still on the device', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const probe = await createGpuModel(new Grid(N), { topography, ocean: { vertical: 'hybrid', everySteps: 1 } }), stress = windStress(probe.mesh);
  probe.destroy();
  for (const [hemisphere, moves] of [['north', false], ['south', true]]) {
    const s = await scenario({ rheology: 'cavitating', hemisphere }, { stress });
    s.stepCpu(); await s.stepGpu();
    const { iceVelocity } = await s.ocean.restartArrays();
    let top = 0; for (const u of iceVelocity) top = Math.max(top, Math.abs(u));
    console.log(`hemisphere ${hemisphere}: the ice of the south moves at up to ${top.toFixed(4)} m/s`);
    assert.equal(top > 0.01, moves, `hemisphere ${hemisphere}`);
    if (moves) await compareVelocity(`hemisphere ${hemisphere}`, s, 1e-2);
    await compareIce(`hemisphere ${hemisphere}`, s, { any: 1e-3 });
    s.gpuModel.destroy();
  }
});

async function cavitating(options) {
  const probe = await createGpuModel(new Grid(N), { topography, ocean: { vertical: 'hybrid', everySteps: 1 } }), stress = windStress(probe.mesh);
  probe.destroy();
  const s = await scenario({ rheology: 'cavitating' }, { stress, ocean: options });
  for (let n = 0; n < 20; n++) { s.stepCpu(); await s.stepGpu(); }
  await compareVelocity('twenty coupled calls', s, 2e-2);
  await compareIce('twenty coupled calls', s, { A: 2e-3, V: 2e-3, S: 2e-3, E: 1e-2 });
  const fields = await s.ocean.download(), wet = Array.from(s.cpu.h, (v) => v > 1);
  const T = Float64Array.from(s.cpu.Q, (q, n) => q / Math.max(1e-9, s.cpu.h[n])), W = Float64Array.from(s.cpu.W, (w, n) => w / Math.max(1e-9, s.cpu.h[n]));
  const top = (a, b) => stats(a.subarray(0, s.C), b.subarray(0, s.C));
  const rT = top(T, fields.T), rS = top(W, fields.S);
  console.log(`ocean top layer after twenty coupled calls: T rms ${rT.rmsRel.toExponential(2)} max ${rT.maxDiff.toExponential(2)} K, S rms ${rS.rmsRel.toExponential(2)} max ${rS.maxDiff.toExponential(2)} psu`);
  assert.ok(rT.rmsRel < 1e-3 && rS.rmsRel < 1e-3 && wet.length > 0);
  for (const f of [fields.h, fields.u, fields.T, fields.S]) for (const v of f) assert.ok(Number.isFinite(v));
  if (options.mstarScheme === 'rh18') assert.ok(rT.maxDiff < 0.1 && rS.maxDiff < 0.02, `the top layer parts by up to ${rT.maxDiff} K and ${rS.maxDiff} psu`);
  s.gpuModel.destroy();
}
test('twenty coupled calls with the cavitating fluid: ice, ice velocity and the ocean under it follow the CPU', { skip: !gpuAvailable && 'webgpu not installed' }, () => cavitating({}));
test('twenty coupled calls with the cavitating fluid, mstarScheme rh18 and langmuirScheme rescale: ice, ice velocity and the ocean under it follow the CPU, its top layer in every cell to 0.1 K', { skip: !gpuAvailable && 'webgpu not installed' }, () => cavitating({ mstarScheme: 'rh18', langmuirScheme: 'rescale' }));

test('inside the model\'s own step every rheology takes the atmosphere\'s stress, validates on the device and keeps the ice finite and below the speed limit', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { createModel } = await import('../js/model.module.js'), { initializeState } = await import('../js/physics/init.module.js');
  const cpuModel = createModel(new Grid(N), { topography }), init = initializeState(cpuModel, {});
  for (const rheology of ['frozen', 'free', 'cavitating']) {
    const model = await createGpuModel(new Grid(N), { topography, ocean: { vertical: 'hybrid', everySteps: 1, iceDrift: { rheology } } });
    for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
    model.load();
    const dev = model.gpu.device, C = model.mesh.nCells;
    let iced = 0; for (let i = 0; i < C; i++) if (model.oceanEngine.cellOcean[i] && init[6][i] > 0) iced++;
    dev.pushErrorScope('validation');
    for (let n = 0; n < 4; n++) await model.step(DT);
    await model.settle();
    const failure = await dev.popErrorScope();
    assert.equal(failure, null, failure && failure.message);
    const { iceVelocity } = await model.oceanEngine.restartArrays();
    let top = 0; for (const u of iceVelocity) { assert.ok(Number.isFinite(u)); top = Math.max(top, Math.abs(u)); }
    const g = await readRanges(dev, model.gpu.buffers.S, [{ offset: model.gpu.layout.S.ICE, length: C }]);
    for (const h of g[0]) assert.ok(Number.isFinite(h) && h >= 0);
    console.log(`${rheology} in the model's step: ${iced} iced cells, top ice speed ${top.toFixed(4)} m/s`);
    assert.ok(top <= 2.0001 && (rheology === 'frozen' ? top === 0 : iced === 0 || top > 1e-4), `${rheology}: ${top} m/s over ${iced} iced cells`);
    model.destroy();
  }
});

test('ice that drifts into open water at a Courant number of 0.005 persists on the GPU as on the CPU, the trace cells below the minimum cover included', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const probe = await createGpuModel(new Grid(N), { topography, ocean: { vertical: 'hybrid', everySteps: 1 } }), pm = probe.mesh, toward = northwardOn(pm, (lat) => lat < -40 * DEG);
  probe.destroy();
  const velocity = Float64Array.from(toward, (d, e) => d * 0.005 * pm.dcEdge[e] / DT);
  const s = await scenario({ rheology: 'frozen' }, { velocity, keep: true });
  const before = cpuFields(s).A.filter((a) => a > 0).length;
  for (let n = 0; n < 5; n++) { s.stepCpu(); await s.stepGpu(); }
  const a = cpuFields(s), b = await gpuFields(s);
  const trace = (f) => Array.from(f.A, (v, i) => i).filter((i) => f.A[i] > 0 && f.A[i] < MINIMUM_CONCENTRATION);
  const onCpu = trace(a), onGpu = trace(b), both = onCpu.filter((i) => b.A[i] > 0);
  let worst = 0, scale = 0;
  for (const i of both) { worst = Math.max(worst, Math.abs(a.A[i] - b.A[i])); scale = Math.max(scale, a.A[i]); }
  const larger = onCpu.filter((i) => a.A[i] > 1e-5), missing = larger.filter((i) => !(b.A[i] > 0));
  console.log(`after five calls the ice covers ${a.A.filter((v) => v > 0).length} cells (${before} to start); ${onCpu.length} on the CPU and ${onGpu.length} on the device hold under 1 percent cover, ${both.length} in common to within ${worst.toExponential(2)} of ${scale.toExponential(2)}`);
  assert.ok(a.A.filter((v) => v > 0).length > before + 10 && b.A.filter((v) => v > 0).length > before + 10, 'the edge advanced on both engines');
  assert.ok(onCpu.length >= 10 && onGpu.length >= 10, `${onCpu.length} and ${onGpu.length} cells under 1 percent cover`);
  assert.equal(missing.length, 0, 'every trace cell above 1e-5 on the CPU is on the device');
  assert.ok(worst < 1e-4 * scale + 1e-8, 'their covers agree');
  await compareIce('five calls of trace arrivals', s, { A: 1e-4, V: 1e-4, S: 1e-4, E: 1e-3 });
  s.gpuModel.destroy();
});

test('over a hundred calls the GPU conserves the ice volume and snow mass to single precision, and stored ice volume follows the transport', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const probe = await createGpuModel(new Grid(N), { topography, ocean: { vertical: 'hybrid', everySteps: 1 } }), pm = probe.mesh;
  const whole = streamOn(pm, 40 / pm.radius), velocity = Float64Array.from(whole, (u, e) => (pm.latCell[pm.cellsOnEdge[2 * e]] < -71 * DEG && pm.latCell[pm.cellsOnEdge[2 * e + 1]] < -71 * DEG ? u : 0));
  probe.destroy();
  const s = await scenario({ rheology: 'frozen' }, { global: true, velocity }), mesh = s.mesh;
  const totals = async () => { const g = await s.readIce(), f = fieldsOf(g.h, g.conc, g.snow, g.ts, s.cpu.cellOcean); let V = 0, S = 0; for (let i = 0; i < s.C; i++) { V += mesh.areaCell[i] * f.V[i]; S += mesh.areaCell[i] * f.S[i]; } return { V, S }; };
  const before = await totals();
  for (let n = 0; n < 100; n++) await s.stepGpu();
  const after = await totals();
  console.log(`100 calls: volume ${(after.V / before.V - 1).toExponential(2)}, snow mass ${(after.S / before.S - 1).toExponential(2)} relative change`);
  assert.ok(Math.abs(after.V / before.V - 1) < 1e-6 && Math.abs(after.S / before.S - 1) < 1e-6);
  s.ocean.driftStep(DT); await s.gpu.device.queue.onSubmittedWorkDone();
  const { previousIceVolume } = await s.ocean.restartArrays(), g = await s.readIce();
  let worst = 0, scale = 0;
  for (let i = 0; i < s.C; i++) if (g.h[i] > 0) { worst = Math.max(worst, Math.abs(previousIceVolume[i] - g.h[i] * g.conc[i])); scale = Math.max(scale, g.h[i] * g.conc[i]); }
  console.log(`the stored ice volume differs from the measured by ${worst.toExponential(2)} m of ${scale.toFixed(3)} m after a transport alone`);
  assert.ok(worst < 1e-5 * scale, 'previousVolume follows the transport');
  s.gpuModel.destroy();
});
