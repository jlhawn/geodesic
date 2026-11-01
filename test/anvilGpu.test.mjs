import { test } from 'node:test';
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { layoutFor, createGpuCore } = gpuAvailable ? await import('../js/gpu/core.gpu.js') : {};
const { PHYSICS_KERNELS, physicsKernel, radiationKernel } = gpuAvailable ? await import('../js/gpu/physics.gpu.js') : {};
const { anvilOn, anvilConstants, anvilKernels } = gpuAvailable ? await import('../js/gpu/anvil.gpu.js') : {};
const { getDevice } = gpuAvailable ? await import('../js/gpu/device.module.js') : {};
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};
const { initializeState } = await import('../js/physics/init.module.js');

const mesh = createModel(new Grid(16), { physics: false }).mesh, C = mesh.nCells;
const sha = (text) => createHash('sha256').update(text).digest('hex').slice(0, 16);

test('the PH layout appends the anvil\'s cover and ice, K0 × C each, after every other section, and has them only with anvil layers', { skip: !gpuAvailable && 'webgpu not installed' }, () => {
  const plain = layoutFor(mesh, 8), tall = layoutFor(mesh, 8, 1, 1, 1, null, true, 1, 1, 8, true, true, true, true, 3, true), anvil = layoutFor(mesh, 8, 0, 0, 0, null, false, 0, 0, false, false, false, false, false, 0, false, 5), tallAnvil = layoutFor(mesh, 8, 1, 1, 1, null, true, 1, 1, 8, true, true, true, true, 3, true, 5);
  assert.equal(plain.PH.ANVA, undefined);
  assert.equal(plain.PH.ANVQ, undefined);
  assert.equal(anvil.PH.total, plain.PH.total + 10 * C);
  assert.equal(anvil.PH.ANVQ - anvil.PH.ANVA, 5 * C);
  assert.equal(anvil.PH.ANVA, plain.PH.total);
  for (const name of Object.keys(plain.PH)) if (name !== 'total') assert.equal(anvil.PH[name], plain.PH[name], `${name} stays where it was`);
  assert.equal(tallAnvil.PH.total, tall.PH.total + 10 * C);
  assert.equal(tallAnvil.PH.ANVA, tall.PH.total);
  assert.deepEqual(Object.keys(tallAnvil.PH).slice(0, -3), Object.keys(tall.PH).slice(0, -1), 'the sections before the anvil are the same in the same order');
});

test('the anvil\'s constants exist only with anvilLifetime, and take the fall speed from anvilFall or else iceFall and the spread from anvilSpread', { skip: !gpuAvailable && 'webgpu not installed' }, () => {
  assert.equal(anvilOn({}), false);
  assert.equal(anvilOn({ anvilLifetime: null }), false);
  assert.equal(anvilConstants({ anvilLifetime: null, iceFall: 2.5 }), '');
  assert.match(anvilConstants({ anvilLifetime: 10800, iceFall: 2.5 }), /ANVIL_LIFE: f32 = 10800;.*ANVIL_FALL_C: f32 = 2\.5;/);
  assert.match(anvilConstants({ anvilLifetime: 10800, iceFall: 2.5, anvilFall: 0.4 }), /ANVIL_FALL_C: f32 = 0\.4;/);
  assert.match(anvilConstants({ anvilLifetime: 10800, iceFall: 2.5 }), /ANVIL_SPREAD: f32 = 1;/);
  assert.match(anvilConstants({ anvilLifetime: 10800, iceFall: 2.5, anvilSpread: 4 }), /ANVIL_SPREAD: f32 = 4;/);
  assert.match(anvilConstants({ anvilLifetime: 10800, iceFall: null }), /ANVIL_FALL_C: f32 = 0;/);
  assert.throws(() => anvilConstants({ anvilLifetime: 0 }), /anvilLifetime/);
  assert.throws(() => anvilConstants({ anvilLifetime: 3600, anvilFall: -1 }), /anvilFall/);
  assert.throws(() => anvilConstants({ anvilLifetime: 3600, anvilSpread: 0 }), /anvilSpread/);
});

test('the patch puts the anvil into the plume, rain and radiation kernels at anchors that each occur once, and refuses kernels it has already patched', { skip: !gpuAvailable && 'webgpu not installed' }, () => {
  for (const layered of [false, true]) {
    const bodies = { adjustMix: PHYSICS_KERNELS.adjustMix, adjustPlume: PHYSICS_KERNELS.adjustPlume, adjustRain: PHYSICS_KERNELS.adjustRain, physics: physicsKernel(false), radiation: radiationKernel(false, layered) };
    const patched = anvilKernels(bodies);
    for (const [name, text] of Object.entries(bodies)) assert.ok(!/ANV/.test(text), `${name} has nothing of the anvil to begin with`);
    assert.ok(/iceCarried/.test(patched.adjustPlume) && /PH_ANVA/.test(patched.adjustPlume) && /PH_ANVQ/.test(patched.adjustPlume), 'the source');
    assert.ok(/anvilFall/.test(patched.adjustRain) && /ANVIL_LIFE/.test(patched.adjustRain) && /descending \+= anvilFall/.test(patched.adjustRain), 'the decay and fall');
    assert.equal(patched.adjustMix, bodies.adjustMix, 'the first adjustment kernel is not touched');
    for (const name of ['physics', 'radiation']) assert.ok(/PH_ANVA \+ k \* C \+ i/.test(patched[name]) && !/k < CU_K0\) \{ return vec2<f32>\(0\.0/.test(patched[name]), `${name} reads the anvil as cumulus cloud above cumulusK0`);
    assert.throws(() => anvilKernels(patched), /exactly once/);
  }
});

test('with anvilLifetime null the kernels the anvil patches compile to the text they always did, and with it set the model holds the anvil on the device', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { device } = await getDevice(), create = device.createShaderModule.bind(device), codes = new Map();
  device.createShaderModule = (descriptor) => { codes.set(descriptor.label, descriptor.code); return create(descriptor); };
  const build = async (physics) => {
    codes.clear();
    const gpu = await createGpuCore(mesh, { levels: sigmaInterfaces('bl48'), physics });
    const hashes = Object.fromEntries(['adjustMix', 'adjustPlume', 'adjustRain', 'radiation', 'physics', 'layered'].map((name) => [name, sha(codes.get(name))]));
    const layout = gpu.layout;
    for (const buffer of Object.values(gpu.buffers)) buffer.destroy();
    return { hashes, layout };
  };
  try {
    const off = await build({}), explicit = await build({ anvilLifetime: null }), on = await build({ anvilLifetime: 3 * 3600 });
    assert.deepEqual(off.hashes, { adjustMix: '898979d967a862c1', adjustPlume: '86c39b9f024becad', adjustRain: 'f282f9e9d102421d', radiation: 'f8a6873b16924f2e', physics: 'b05ebdef3f5d1f03', layered: 'cca38db233ea28e6' });
    assert.deepEqual(explicit.hashes, off.hashes);
    assert.equal(off.layout.PH.ANVA, undefined);
    assert.equal(off.layout.PH.total, 1811317);
    for (const name of ['adjustPlume', 'adjustRain', 'radiation', 'physics', 'layered']) assert.notEqual(on.hashes[name], off.hashes[name], `${name} changes with the anvil`);
    const K0 = on.layout.K - (on.layout.PH.CUWATER - on.layout.PH.CUCOVER) / C;
    assert.equal(on.layout.PH.ANVQ - on.layout.PH.ANVA, K0 * C);
    assert.equal(on.layout.PH.total, off.layout.PH.total + 2 * K0 * C);
    console.log(`N=16 on bl48: the adjustment, radiation and physics kernels hash as before with the anvil off; with it on they change and the PH buffer grows by 2 × ${K0} × ${C} floats`);
  } finally {
    device.createShaderModule = create;
  }
});

test('a GPU model mirrors the anvil on the host as K0 × C arrays when it has one and as empty arrays when it has none', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const off = await createGpuModel(new Grid(16), {}), on = await createGpuModel(new Grid(16), { moist: { anvilLifetime: 3 * 3600 } });
  const K0 = on.moist.cumulusK0;
  assert.equal(off.moist.anvilCover.length, 0);
  assert.equal(off.moist.anvilIce.length, 0);
  assert.equal(on.moist.anvilCover.length, K0 * C);
  assert.equal(on.moist.anvilIce.length, K0 * C);
  off.destroy(); on.destroy();
});

test('the page\'s cloud totals and cloud heights count the anvil as they count the plumes\' cumulus, in the same path as the CPU model does', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const model = await createGpuModel(new Grid(16), { ocean: false, moist: { anvilLifetime: 3 * 3600 } });
  const { state, gpu } = model, { K, sigmaMid, dSigma, g } = model.core.diagnostics, K0 = model.moist.cumulusK0, k = Math.floor(K0 / 2);
  initializeState(model, {}).forEach((values, a) => state[a].set(values));
  model.load();
  for (let n = 0; n < 4; n++) await model.step(900);
  await model.sync();
  const fields = ['cloudCumulus', 'cloudHigh', 'cloudCover', 'cloudTop', 'cloudBase'];
  const before = (await model.beginFrame({ level: 'surface', fields, diagnostics: true })).fields;
  const cover = new Float32Array(C).fill(0.3), ice = new Float32Array(C).fill(5e-5);
  gpu.device.queue.writeBuffer(gpu.buffers.PH, 4 * (gpu.layout.PH.ANVA + k * C), cover);
  gpu.device.queue.writeBuffer(gpu.buffers.PH, 4 * (gpu.layout.PH.ANVQ + k * C), ice);
  const after = (await model.beginFrame({ level: 'surface', fields, diagnostics: true })).fields;
  const [pi] = await gpu.download();
  let worst = 0;
  for (let i = 0; i < C; i++) {
    const added = pi[i] * dSigma[k] / g * 5e-5;
    worst = Math.max(worst, Math.abs(after.cloudCumulus[i] - before.cloudCumulus[i] - added) / added);
    assert.equal(after.cloudHigh[i], before.cloudHigh[i], 'the plumes\' cloud is not a class of the resolved cloud');
  }
  console.log(`one anvil layer of ${(5e-5 * 1e3).toFixed(2)} g/kg at σ ${sigmaMid[k].toFixed(3)} adds its path to cloudCumulus to ${worst.toExponential(1)} of itself; the cloud heights move in ${after.cloudTop.some((x, i) => x !== before.cloudTop[i])}`);
  assert.ok(worst < 1e-4, `the anvil's path is added to ${worst} of itself`);
  assert.ok(after.cloudTop.some((x, i) => x !== before.cloudTop[i]) || after.cloudBase.some((x, i) => x !== before.cloudBase[i]), 'the cloud top or base sees the anvil somewhere');
  assert.ok(K > K0);
  model.destroy();
});
