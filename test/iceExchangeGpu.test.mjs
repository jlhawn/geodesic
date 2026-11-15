import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { lakeArrays } from '../js/physics/lake.module.js';
import { loadSavedLand, loadSavedLake } from '../scripts/savedLand.mjs';
import { buildMesh } from '../js/mesh.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { iceExchangeOptions, FREEZING_POINT } from '../js/physics/ice.module.js';

/*
 * The GPU model with the ocean's iceExchange 'turbulent' (physics/ice.module.js): the layout of the
 * physics buffer and the kernels carry the exchange only when it is asked for, the model steps with
 * ice, the water under the ice holds heat above freezing, and the ocean's save and load carry the
 * water's heat and the conductance through the device bit for bit.
 */
let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};
const { layoutFor } = gpuAvailable ? await import('../js/gpu/core.gpu.js') : {};
const { physicsSurfaceKernel, seaSurfaceWgsl } = gpuAvailable ? await import('../js/gpu/physics.gpu.js') : {};
const { readRanges } = gpuAvailable ? await import('../js/gpu/device.module.js') : {};
const N = 16, DT = 1350;

test('the physics buffer, the surface kernel and the constants carry the turbulent exchange only when it is asked for', { skip: !gpuAvailable && 'webgpu not installed' }, () => {
  const mesh = buildMesh(new Grid(N));
  const plain = layoutFor(mesh, 8), asked = layoutFor(mesh, 8, ...new Array(22), true);
  assert.ok(plain.PH.WHEAT === undefined && plain.PH.CONDUCT === undefined && plain.PH.WET === undefined);
  assert.equal(asked.PH.WHEAT, plain.PH.total);
  assert.equal(asked.PH.CONDUCT, plain.PH.total + mesh.nCells);
  assert.equal(asked.PH.WET, plain.PH.total + 2 * mesh.nCells);
  assert.equal(asked.PH.total, plain.PH.total + 3 * mesh.nCells);
  const instant = { iceTransfer: iceExchangeOptions('instant') }, turbulent = { iceTransfer: iceExchangeOptions('turbulent') };
  assert.equal(seaSurfaceWgsl(instant), seaSurfaceWgsl({}));
  assert.equal(physicsSurfaceKernel(false), physicsSurfaceKernel(false, false, null, false, false, instant));
  assert.ok(!/PH_WHEAT|PH_CONDUCT|XFER_/.test(physicsSurfaceKernel(false)) && /PH_WHEAT/.test(physicsSurfaceKernel(false, false, null, false, false, turbulent)));
});

test('the GPU model steps with ice over a turbulent exchange: the water under the ice holds heat above freezing, nothing goes wrong, and a save and load of the ocean give the device its heat and conductance back', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const cpu = createModel(new Grid(N), { ocean: false }), init = initializeState(cpu, {});
  const model = await createGpuModel(new Grid(N), { ocean: { iceExchange: 'turbulent' } });
  const C = model.mesh.nCells, { device, buffers, layout } = model.gpu;
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  model.load();
  const read = async () => {
    const [heat, conductance, flux, cap] = await readRanges(device, buffers.PH, ['WHEAT', 'CONDUCT', 'OFLUX', 'CAP'].map((name) => ({ offset: layout.PH[name], length: C })));
    const [ice] = await readRanges(device, buffers.S, [{ offset: layout.S.ICE, length: C }]);
    return { heat: Float32Array.from(heat), conductance: Float32Array.from(conductance), flux: Float32Array.from(flux), cap: Float32Array.from(cap), ice: Float32Array.from(ice) };
  };
  const start = await read();
  assert.ok(start.ice.some((h) => h > 0), 'the start has ice');
  for (let n = 0; n < 16; n++) await model.step(DT);
  const after = await read();
  let iced = 0, held = 0, top = 0, floor = Infinity;
  for (let i = 0; i < C; i++) {
    assert.ok(Number.isFinite(after.heat[i]) && Number.isFinite(after.conductance[i]) && Number.isFinite(after.ice[i]), `cell ${i} is finite`);
    if (!(after.ice[i] > 0)) continue;
    iced++;
    assert.ok(after.heat[i] >= 0, `cell ${i}: the water under the ice holds ${after.heat[i]} J/m²`);
    if (after.heat[i] > 0) held++;
    top = Math.max(top, after.heat[i] / after.cap[i]); floor = Math.min(floor, after.conductance[i]);
  }
  console.log(`after 16 steps at N=${N}: ${iced} iced cells, ${held} with water above freezing under them (up to ${top.toFixed(3)} K), the conductance at least ${floor.toFixed(1)} W/m²/K`);
  assert.ok(iced > 0 && floor >= 0.99 * 4.08e6 * 0.006 * 5e-4, `${iced} iced cells, conductance floor ${floor}`);
  const saved = await model.ocean.serialize({ restart: true });
  assert.deepEqual(Array.from(saved.iceWater), Array.from(after.heat));
  assert.deepEqual(Array.from(saved.iceConductance), Array.from(after.conductance));
  const marked = Float64Array.from(saved.iceWater, (x, i) => (x > 0 ? 0.5 * x : x)), changed = { ...saved, iceWater: marked };
  model.ocean.load(changed, model.state[3], model.state[6]);
  const loaded = await read();
  for (let i = 0; i < C; i++) assert.equal(loaded.heat[i], Math.fround(marked[i]), `cell ${i}: the heat loaded`);
  assert.deepEqual(Array.from(loaded.conductance), Array.from(after.conductance), 'the conductance loaded');
  const instant = await createGpuModel(new Grid(N), { ocean: {} });
  assert.equal(instant.gpu.layout.PH.WHEAT, undefined, 'the default model has no water heat on the device');
  model.destroy(); instant.destroy();
});

test('winter full pack over water held at the f32 freezing point keeps a cover of exactly 1 on both engines, and over water 0.3 K above freezing it opens alike on both (iceExchange turbulent)', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const options = { ocean: false, ice: { iceExchange: 'turbulent', oceanCells: new Uint8Array(10 * N * N + 2).fill(1) } };
  const cpu = createModel(new Grid(N), options), gpu = await createGpuModel(new Grid(N), options);
  const C = cpu.mesh.nCells, init = initializeState(cpu, {}), packed = [];
  for (let i = 0; i < C; i++) if (cpu.mesh.latCell[i] > 70 * Math.PI / 180) { packed.push(i); init[6][i] = 2; init[3][i] = 250; }
  for (const a of init) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  const full = Float64Array.from({ length: C }, (_, i) => (packed.includes(i) ? 1 : 0)), warm = (i) => packed.indexOf(i) % 2 === 1;
  for (const m of [cpu, gpu]) { for (let a = 0; a < init.length; a++) m.state[a].set(init[a]); m.seaIce.load(m.state[6], full); m.time = 270 * 86400; }
  for (const i of packed) cpu.seaIce.holdWater(i, warm(i) ? FREEZING_POINT + 0.3 : Math.fround(FREEZING_POINT));
  const water = Float32Array.from(cpu.seaIce.waterHeat), conductance = Float32Array.from({ length: C }, (_, i) => (packed.includes(i) ? 500 : 0));
  cpu.seaIce.waterHeat.set(water); cpu.seaIce.conductance.set(conductance);
  gpu.load();
  const { device, buffers, layout } = gpu.gpu;
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.WHEAT, water); device.queue.writeBuffer(buffers.PH, 4 * layout.PH.CONDUCT, conductance);
  for (let n = 0; n < 4; n++) { cpu.step(DT); await gpu.step(DT); }
  await gpu.sync();
  const [conc] = await readRanges(device, buffers.PH, [{ offset: layout.PH.CONC, length: C }]);
  let held = 0, opened = 0, gap = 0, least = 1;
  for (const i of packed) {
    if (!warm(i)) { assert.equal(cpu.seaIce.concentration[i], 1, `cell ${i} on the CPU`); assert.equal(conc[i], 1, `cell ${i} on the GPU`); held++; continue; }
    assert.ok(cpu.seaIce.concentration[i] < 1 - 1e-5 && conc[i] < 1 - 1e-5, `cell ${i} opens: CPU ${cpu.seaIce.concentration[i]}, GPU ${conc[i]}`);
    gap = Math.max(gap, Math.abs(cpu.seaIce.concentration[i] - conc[i])); least = Math.min(least, cpu.seaIce.concentration[i]); opened++;
  }
  console.log(`after 4 steps at N=${N}: ${held} packed cells over water at the f32 freezing point at exactly 1 on both engines; ${opened} over water 0.3 K warmer opened to as low as ${least.toFixed(6)}, the engines within ${gap.toExponential(1)}`);
  assert.ok(held > 0 && opened > 0 && gap < 1e-6, `${held} held, ${opened} opened, gap ${gap}`);
  cpu.destroy?.(); gpu.destroy();
});

test('over cells with no ocean column under the ice the turbulent exchange holds no heat and the leads\' gain melts ice at once, alike on both engines', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const mesh = buildMesh(new Grid(N)), C = mesh.nCells, polar = [];
  for (let i = 0; i < C; i++) if (mesh.latCell[i] < -70 * Math.PI / 180) polar.push(i);
  const dry = (i) => polar.indexOf(i) % 2 === 0, oceanCells = Uint8Array.from({ length: C }, (_, i) => (polar.includes(i) && dry(i) ? 0 : 1));
  const cpu = createModel(new Grid(N), { ocean: false, ice: { iceExchange: 'turbulent', oceanCells } }), gpu = await createGpuModel(new Grid(N), { ocean: false, ice: { iceExchange: 'turbulent', oceanCells } });
  const init = initializeState(cpu, {});
  for (const i of polar) { init[6][i] = 1; init[3][i] = 268; }
  for (const a of init) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  const cover = Float64Array.from({ length: C }, (_, i) => (polar.includes(i) ? 0.8 : 0));
  for (const m of [cpu, gpu]) { for (let a = 0; a < init.length; a++) m.state[a].set(init[a]); m.seaIce.load(m.state[6], cover); m.time = 270 * 86400; }
  const conductance = Float32Array.from({ length: C }, (_, i) => (polar.includes(i) && !dry(i) ? 500 : 0));
  cpu.seaIce.conductance.set(conductance);
  gpu.load();
  const { device, buffers, layout } = gpu.gpu;
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.CONDUCT, conductance);
  const [onDevice] = await readRanges(device, buffers.PH, [{ offset: layout.PH.WET, length: C }]);
  assert.deepEqual(Array.from(onDevice), Array.from(oceanCells), 'the device takes the ice\'s oceanCells as the cells with water under the ice');
  for (let n = 0; n < 4; n++) { cpu.step(DT); await gpu.step(DT); }
  await gpu.sync();
  const [conc, wheat] = await readRanges(device, buffers.PH, [{ offset: layout.PH.CONC, length: C }, { offset: layout.PH.WHEAT, length: C }]);
  const [ice] = await readRanges(device, buffers.S, [{ offset: layout.S.ICE, length: C }]);
  let held = 0, gap = 0, dryLoss = 0, wetLoss = 0, nDry = 0;
  for (const i of polar) {
    const loss = 0.8 - cpu.seaIce.concentration[i] * cpu.state[6][i];
    if (!dry(i)) { if (cpu.seaIce.waterHeat[i] > 0 && wheat[i] > 0) held++; wetLoss += loss; continue; }
    assert.equal(cpu.seaIce.waterHeat[i], 0, `cell ${i} holds heat on the CPU`); assert.equal(wheat[i], 0, `cell ${i} holds heat on the GPU`);
    gap = Math.max(gap, Math.abs(cpu.seaIce.concentration[i] - conc[i]), Math.abs(cpu.seaIce.concentration[i] * cpu.state[6][i] - conc[i] * ice[i]));
    dryLoss += loss; nDry++;
  }
  const wet = polar.length - nDry;
  console.log(`after 4 steps at N=${N}: ${nDry} packed cells with no ocean under them hold no heat on either engine and lose ${(dryLoss / nDry).toExponential(2)} m of ice each, the engines within ${gap.toExponential(1)}; ${held} of ${wet} over the ocean hold the leads' gain and lose ${(wetLoss / wet).toExponential(2)} m`);
  assert.ok(nDry > 0 && held > 0 && dryLoss / nDry > wetLoss / wet && gap < 1e-5, `${nDry} dry, ${held} holding, losses ${dryLoss / nDry} and ${wetLoss / wet}, gap ${gap}`);
  cpu.destroy?.(); gpu.destroy();
});

const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer), TURBULENT = { mode: 'turbulent', frictionMinimum: 2e-2 };
const host = (model) => Object.entries({
  ...Object.fromEntries(model.state.map((values, a) => [`state${a}`, values])),
  concentration: model.seaIce.concentration, brine: model.seaIce.brine, windSpeed: model.surface.windSpeed, seaFlow: model.seaFlow,
  ...Object.fromEntries(['mlmSubsidence', 'mlmHeight', 'mlmGate', 'meanAbsorbedSolar', 'meanOutgoingLongwave', 'meanPlanetaryAlbedo', 'meanShortwaveCloudEffect', 'meanLongwaveCloudEffect', 'evaporation'].map((name) => [name, model.radiation[name]])),
  ...Object.fromEntries(['convectiveRain', 'largeScaleRain', 'cumulusCover', 'cumulusWater', 'subcloudVirtual'].map((name) => [name, model.moist[name]])),
  ...Object.fromEntries(['depth', 'mixingTop', 'regime', 'buoyancyFlux'].map((name) => [`boundary.${name}`, model.boundaryLayer[name]])),
  ...(model.exchange ? { exchangeWind: model.exchange.wind, ...(model.exchange.fixed ? {} : { exchangeHeat: model.exchange.heat }) } : {}),
  ...(model.skin ? { warmLayer: model.skin.warm, coolSkin: model.skin.cool } : {}),
}).filter(([, values]) => values);
const begin = async (options) => {
  const model = await createGpuModel(new Grid(N), options);
  initializeState(model, { geostrophic: !model.surfaceGeopotential }).forEach((values, a) => model.state[a].set(values));
  for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i]) model.state[6][i] = 0;
  model.load();
  model.ocean.initialize(model.state[3], model.state[6]);
  model.land.initialize();
  return model;
};
const save = async (model) => {
  await model.sync();
  const saved = { time: model.time, host: Object.fromEntries(host(model).map(([name, values]) => [name, Float64Array.from(values)])), fresh: await model.readFreshwater(), ocean: await model.ocean.serialize({ restart: true }), land: await model.land.serialize(), lake: model.lake ? lakeArrays(await model.lake.serialize()) : {} };
  model.destroy();
  return saved;
};
const restore = async (saved, options) => {
  const model = await createGpuModel(new Grid(N), options), { state } = model;
  state.forEach((values, a) => values.set(saved.host[`state${a}`]));
  model.seaIce.load(state[6], saved.host.concentration, saved.host.brine ?? null);
  for (const [name, values] of host(model)) if (!/^state|^concentration$|^brine$/.test(name)) values.set(saved.host[name]);
  model.time = saved.time;
  model.load();
  model.ocean.load(saved.ocean, state[3], state[6]);
  loadSavedLand(model, saved);
  if (model.lake) loadSavedLake(model, saved.lake, { quiet: true });
  model.writeFreshwater(saved.fresh);
  return model;
};
const fields = (value, prefix = '') => (ArrayBuffer.isView(value) || (Array.isArray(value) && typeof value[0] === 'number') ? [[prefix, value]] : value && typeof value === 'object' ? Object.entries(value).flatMap(([name, part]) => fields(part, prefix ? `${prefix}.${name}` : name)) : []);
const steps = async (model, n) => { for (let s = 0; s < n; s++) await model.step(DT); return model; };

test('a GPU run over land with iceExchange turbulent, saved by serialize({ restart: true }) on every cadence and loaded into a fresh model, steps on bit for bit with the run that did not stop: the state, the ice, the ocean with the water\'s heat and the conductance, and the land', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const options = { topography, ocean: { iceExchange: TURBULENT }, ice: { iceAlbedoPhases: true } };
  const whole = await begin(options), gcd = (a, b) => (b ? gcd(b, a % b) : a), { radiationEvery, dragEvery, oceanEvery } = whole.cadence;
  const every = [radiationEvery, dragEvery, oceanEvery].reduce((a, b) => (a * b) / gcd(a, b)), n = every * Math.ceil(16 / every);
  const ended = await save(await steps(whole, 2 * n)), saved = await save(await steps(await begin(options), n));
  const held = saved.ocean.iceWater.filter((x) => x > 0).length, stirred = saved.ocean.iceConductance.filter((x) => x > 0).length;
  assert.ok(held > 0 && stirred > 0, `the save after ${n} steps holds the water's heat in ${held} cells and a conductance in ${stirred}`);
  const after = await save(await steps(await restore(saved, options), n)), theirs = Object.fromEntries(fields(after)), differ = [];
  for (const [name, values] of fields(ended)) {
    const other = theirs[name];
    if (!other || other.length !== values.length) { differ.push(`${name} missing`); continue; }
    let count = 0, largest = 0;
    for (let i = 0; i < values.length; i++) if (!Object.is(values[i], other[i])) { count++; largest = Math.max(largest, Math.abs(values[i] - other[i])); }
    if (count) differ.push(`${name} in ${count} of ${values.length}, up to ${largest.toExponential(1)}`);
  }
  console.log(`N=${N}, ${2 * n} steps whole against ${n}, a save and a load, then ${n}: ${fields(ended).length} arrays compared; at the save ${held} cells hold the water's heat and ${stirred} a conductance`);
  assert.deepEqual(differ, [], `whole against resumed: ${differ.join('; ')}`);
});

test('an instant save with ice loaded into the turbulent model on both engines: the water under the ice starts within the f32 grain of the freezing point, as on the CPU, and over two ocean intervals the ice, the held heat and the conductance follow the CPU and no ice collapses', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const options = { topography, ocean: { iceExchange: TURBULENT } }, source = await begin({ topography }), { oceanEvery } = source.cadence;
  const saved = await save(await steps(source, 2 * oceanEvery)), C = saved.host.state6.length;
  assert.ok(!('iceWater' in saved.ocean) && saved.ocean.flux.some((x) => x > 0), 'an instant save with heat in transit to the ice');
  const place = (model) => {
    model.state.forEach((values, a) => values.set(saved.host[`state${a}`]));
    model.seaIce.load(model.state[6], saved.host.concentration, null);
    model.time = saved.time;
  };
  const gpu = await createGpuModel(new Grid(N), options), cpu = createModel(new Grid(N), options), { device, buffers, layout } = gpu.gpu;
  place(gpu); gpu.load(); gpu.ocean.load(saved.ocean, gpu.state[3], gpu.state[6]); loadSavedLand(gpu, saved);
  place(cpu); cpu.ocean.load(saved.ocean, cpu.state[3], cpu.state[6], cpu.seaIce.concentration); loadSavedLand(cpu, saved);
  const capacity = Float64Array.from(cpu.ocean.capacity), iced = [...Array(C).keys()].filter((i) => cpu.state[6][i] > 0 && cpu.ocean.cellOcean[i]);
  const onDevice = async () => { const [heat, conductance] = await readRanges(device, buffers.PH, ['WHEAT', 'CONDUCT'].map((name) => ({ offset: layout.PH[name], length: C }))); return { heat, conductance }; };
  const rise = (heat) => Math.max(...iced.map((i) => heat[i] / capacity[i]));
  const apart = (a, b, scale = () => 1) => Math.max(...iced.map((i) => Math.abs(a[i] - b[i]) / scale(i)));
  const loaded = await onDevice(), [wet] = await readRanges(device, buffers.PH, [{ offset: layout.PH.WET, length: C }]);
  assert.ok(gpu.oceanEngine.cellOcean.some((x) => !x) && wet.every((x, i) => x === gpu.oceanEngine.cellOcean[i]), 'the device\'s cells with water under the ice are the ocean\'s');
  const atLoad = { cpu: rise(cpu.seaIce.waterHeat), gpu: rise(loaded.heat), heat: apart(cpu.seaIce.waterHeat, loaded.heat, (i) => capacity[i]), conductance: apart(cpu.seaIce.conductance, loaded.conductance, (i) => cpu.seaIce.conductance[i]) };
  const volume = (model) => Float64Array.from(model.state[6], (h, i) => (h > 0 ? h * model.seaIce.concentration[i] : 0)), volume0 = volume(cpu);
  for (let n = 0; n < 2 * oceanEvery; n++) { await gpu.step(DT); cpu.step(DT); }
  await gpu.sync();
  const after = await onDevice(), [vg, vc] = [volume(gpu), volume(cpu)];
  const stepped = { heat: apart(cpu.seaIce.waterHeat, after.heat, (i) => capacity[i]), conductance: apart(cpu.seaIce.conductance, after.conductance, (i) => cpu.seaIce.conductance[i]), volume: apart(vc, vg), rise: rise(after.heat) };
  const kept = Math.min(...iced.map((i) => vg[i] / volume0[i])), total = vg.reduce((s, v, i) => s + v * gpu.mesh.areaCell[i], 0) / volume0.reduce((s, v, i) => s + v * gpu.mesh.areaCell[i], 0);
  console.log(`${iced.length} iced cells; at the load the water under the ice is up to ${atLoad.cpu.toExponential(2)} K above freezing on the CPU and ${atLoad.gpu.toExponential(2)} K on the GPU, apart by ${atLoad.heat.toExponential(2)} K, the conductance by ${atLoad.conductance.toExponential(2)}; after ${2 * oceanEvery} steps up to ${stepped.rise.toExponential(2)} K, the engines apart by ${stepped.heat.toExponential(2)} K, ${stepped.conductance.toExponential(2)} in the conductance and ${stepped.volume.toExponential(2)} m in the ice volume; the GPU ice keeps ${kept.toFixed(5)} of its volume in every cell and ${total.toFixed(6)} in all`);
  const grain = Math.fround(FREEZING_POINT) * 2 ** -23;
  assert.ok(atLoad.cpu < grain && atLoad.gpu < grain && atLoad.heat < 1e-9 && atLoad.conductance < 1e-6, `the held heat at the load beyond the f32 grain of the freezing point, ${grain} K: ${JSON.stringify(atLoad)}`);
  assert.ok(stepped.heat < 1e-3 && stepped.conductance < 1e-5 && stepped.volume < 1e-5, JSON.stringify(stepped));
  assert.ok(kept > 0.9 && total > 0.99, `kept ${kept}, ${total} in all`);
  gpu.destroy();
});
