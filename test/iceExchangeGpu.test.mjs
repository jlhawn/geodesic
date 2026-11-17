import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { Grid } from '../js/grid.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { lakeArrays } from '../js/physics/lake.module.js';
import { loadSavedLand, loadSavedLake } from '../scripts/savedLand.mjs';
import { buildMesh } from '../js/mesh.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { HYBRID_LAYERS } from '../js/ocean/column.module.js';
import { iceExchangeOptions, FREEZING_POINT, DEFAULT_ICE_EXCHANGE, DEFAULT_ICE_ALBEDO_PHASES } from '../js/physics/ice.module.js';

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
const { readRanges, getDevice } = gpuAvailable ? await import('../js/gpu/device.module.js') : {};
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

test('the GPU model at N=16 over the aquaplanet compiles its pinned WGSL, 85 modules each: under iceExchange \'instant\' of digest f99337a14076f86f and under \'turbulent\' of digest 35bb61260615ba48 without iceAlbedoPhases, and at the defaults, the turbulent exchange with frictionMinimum 2e-2 and iceAlbedoPhases, of digest b60a561d91ae451f', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const sha = (text) => createHash('sha256').update(text).digest('hex').slice(0, 16), { device } = await getDevice(), compile = device.createShaderModule;
  const compiled = async (options) => {
    const modules = {};
    let n = 0;
    device.createShaderModule = (descriptor) => { modules[`${descriptor.label ?? 'anonymous'}#${n++}`] = sha(descriptor.code); return compile.call(device, descriptor); };
    try { (await createGpuModel(new Grid(N), options)).destroy(); } finally { device.createShaderModule = compile; }
    const listed = Object.entries(modules).sort(([a], [b]) => a.localeCompare(b));
    return { modules: listed.length, digest: sha(listed.map(([name, code]) => `${name}:${code}`).join('\n')) };
  };
  const phasesOff = { iceAlbedoPhases: false };
  const pinned = { instant: await compiled({ ocean: { iceExchange: 'instant' }, ice: phasesOff }), turbulent: await compiled({ ocean: { iceExchange: 'turbulent' }, ice: phasesOff }), standard: await compiled({}) };
  assert.deepEqual(pinned, { instant: { modules: 85, digest: 'f99337a14076f86f' }, turbulent: { modules: 85, digest: '35bb61260615ba48' }, standard: { modules: 85, digest: 'b60a561d91ae451f' } }, 'a change to a kernel updates its pin');
  assert.deepEqual(pinned.standard, await compiled({ ocean: { iceExchange: DEFAULT_ICE_EXCHANGE }, ice: { iceAlbedoPhases: DEFAULT_ICE_ALBEDO_PHASES } }));
});

test('a model given no exchange takes DEFAULT_ICE_EXCHANGE on both engines, in the ocean and in the sea ice, the GPU ocean as written, and DEFAULT_ICE_ALBEDO_PHASES', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const cpu = createModel(new Grid(N), {}), gpu = await createGpuModel(new Grid(N), {}), resolved = iceExchangeOptions(DEFAULT_ICE_EXCHANGE);
  try {
    assert.deepEqual({ cpuOcean: cpu.ocean.iceExchange, cpuIce: cpu.seaIce.iceExchange, gpuOcean: gpu.oceanEngine.options.iceExchange, gpuIce: gpu.seaIce.iceExchange, cpuPhases: cpu.seaIce.iceAlbedoPhases, gpuPhases: gpu.seaIce.iceAlbedoPhases }, { cpuOcean: resolved, cpuIce: resolved, gpuOcean: DEFAULT_ICE_EXCHANGE, cpuPhases: DEFAULT_ICE_ALBEDO_PHASES, gpuPhases: DEFAULT_ICE_ALBEDO_PHASES, gpuIce: resolved });
    assert.equal(gpu.gpu.layout.PH.WHEAT === undefined, !resolved.turbulent, 'the physics buffer holds the water\'s heat only under a turbulent exchange');
  } finally { gpu.destroy(); }
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
  const instant = await createGpuModel(new Grid(N), { ocean: { iceExchange: 'instant' } });
  assert.equal(instant.gpu.layout.PH.WHEAT, undefined, 'an instant model has no water heat on the device');
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

test('thin ice over held heat melts away in a step and leaves the water no heat held, its heat in the open water, alike on both engines (iceExchange turbulent)', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const options = { ocean: false, ice: { iceExchange: 'turbulent', oceanCells: new Uint8Array(10 * N * N + 2).fill(1) } };
  const cpu = createModel(new Grid(N), options), gpu = await createGpuModel(new Grid(N), options), C = cpu.mesh.nCells, init = initializeState(cpu, {});
  const thin = [...Array(C).keys()].filter((i) => cpu.mesh.latCell[i] < -70 * Math.PI / 180);
  for (const i of thin) { init[6][i] = 0.004; init[3][i] = 268; }
  for (const a of init) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  const cover = Float64Array.from({ length: C }, (_, i) => (thin.includes(i) ? 0.5 : 0));
  for (const m of [cpu, gpu]) { for (let a = 0; a < init.length; a++) m.state[a].set(init[a]); m.seaIce.load(m.state[6], cover); m.time = 270 * 86400; }
  for (const i of thin) cpu.seaIce.holdWater(i, FREEZING_POINT + 2);
  const water = Float32Array.from(cpu.seaIce.waterHeat), conductance = Float32Array.from({ length: C }, (_, i) => (thin.includes(i) ? 2000 : 0));
  cpu.seaIce.waterHeat.set(water); cpu.seaIce.conductance.set(conductance);
  gpu.load();
  const { device, buffers, layout } = gpu.gpu;
  device.queue.writeBuffer(buffers.PH, 4 * layout.PH.WHEAT, water); device.queue.writeBuffer(buffers.PH, 4 * layout.PH.CONDUCT, conductance);
  cpu.step(DT); await gpu.step(DT); await gpu.sync();
  const [heat] = await readRanges(device, buffers.PH, [{ offset: layout.PH.WHEAT, length: C }]);
  let gap = 0;
  for (const i of thin) {
    assert.deepEqual([cpu.state[6][i], cpu.seaIce.waterHeat[i], gpu.state[6][i], heat[i]], [0, 0, 0, 0], `cell ${i}: ice and held heat on the CPU, then on the GPU`);
    assert.ok(cpu.state[3][i] > FREEZING_POINT + 1 && gpu.state[3][i] > FREEZING_POINT + 1, `cell ${i}: open water at ${cpu.state[3][i]} and ${gpu.state[3][i]} K`);
    gap = Math.max(gap, Math.abs(cpu.state[3][i] - gpu.state[3][i]));
  }
  console.log(`${thin.length} cells of 4 mm ice over water 2 K above freezing melt away in a step on both engines, the store left empty, the open water within ${gap.toExponential(2)} K`);
  assert.ok(gap < 1e-3, `the open water apart by ${gap} K`);
  gpu.destroy();
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

test('a lean turbulent save, serialize() without the restart arrays, carries the water\'s heat under the ice and the conductance, and a fresh turbulent GPU model loads them as they were saved', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const options = { topography, ocean: { iceExchange: TURBULENT }, ice: { iceAlbedoPhases: true } };
  const model = await steps(await begin(options), 16);
  await model.sync();
  const saved = { time: model.time, host: Object.fromEntries(host(model).map(([name, values]) => [name, Float64Array.from(values)])), fresh: await model.readFreshwater(), ocean: await model.ocean.serialize(), land: await model.land.serialize(), lake: model.lake ? lakeArrays(await model.lake.serialize()) : {} };
  model.destroy();
  const C = saved.host.state6.length, held = saved.ocean.iceWater?.filter((x) => x > 0).length ?? 0;
  assert.equal(saved.ocean.Qref, undefined, 'the save is lean');
  for (const name of ['iceWater', 'iceConductance']) assert.equal(saved.ocean[name]?.length, C, `the lean save holds ocean.${name}`);
  assert.ok(held > 0, `the lean save holds the water's heat in ${held} cells`);
  const back = await restore(saved, options), again = await back.ocean.serialize();
  back.destroy();
  for (const name of ['iceWater', 'iceConductance']) assert.deepEqual([...again[name].keys()].filter((i) => !Object.is(again[name][i], saved.ocean[name][i])), [], `ocean.${name} loads as saved`);
});

test('a turbulent save, restart or lean, loaded into an explicitly instant GPU model books the water\'s heat under the ice into the top layer, logged, and loads every other layer and cell as the same save without that heat and the conductance', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const options = { topography, ocean: { iceExchange: TURBULENT }, ice: { iceAlbedoPhases: true } }, instant = { ...options, ocean: { iceExchange: 'instant' } };
  const model = await steps(await begin(options), 17);
  await model.sync();
  const lean = await model.ocean.serialize(), saved = await save(model), cpu = createModel(new Grid(N), instant).ocean, C = saved.host.state6.length;
  const iced = new Set([...Array(C).keys()].filter((i) => saved.host.state6[i] > 0 && cpu.cellOcean[i]));
  for (const [kind, ocean] of [['restart', saved.ocean], ['lean', lean]]) {
    const { iceWater, iceConductance, ...without } = ocean;
    assert.ok(iceWater.some((x) => x > 0) && iceConductance.some((x) => x > 0), `the ${kind} save holds heat and a conductance under the ice`);
    const notes = [], warn = console.warn;
    console.warn = (note) => notes.push(note);
    let taken;
    try { taken = await restore({ ...saved, ocean }, instant); } finally { console.warn = warn; }
    const plain = await restore({ ...saved, ocean: without }, instant), [top, bare] = [await taken.ocean.serialize(), await plain.ocean.serialize()];
    assert.ok(!('iceWater' in top), 'the instant model saves no water heat');
    const booked = notes.map((note) => note.match(/booked the saved ocean\.iceWater, the heat the water holds above freezing under (\d+) iced cells/)).filter(Boolean);
    assert.equal(booked.length, 1, `the ${kind} load logs its booking once`);
    const changed = [...iced].filter((i) => !Object.is(top.T[i], bare.T[i])), elsewhere = [...top.T.keys()].filter((n) => !(n < C && iced.has(n)) && !Object.is(top.T[n], bare.T[n]));
    assert.deepEqual(elsewhere, [], `from the ${kind} save, only the top layer of iced cells takes the booking`);
    assert.ok(changed.length > 0 && changed.length <= Number(booked[0][1]), `${changed.length} booked cells differ from the save's layer, of ${booked[0][1]} booked`);
    const worst = Math.max(...changed.map((i) => Math.abs((top.T[i] - FREEZING_POINT) * top.h[i] * cpu.rhoCp - iceWater[i]) / Math.max(1, Math.abs(iceWater[i]))));
    assert.ok(worst < 1e-4, `the top layer holds the booked heat to ${worst} of it`);
    taken.destroy(); plain.destroy();
  }
});

const handedHeat = async (thin, bounds) => {
  const source = await begin({ topography, ocean: { iceExchange: TURBULENT } }), { oceanEvery } = source.cadence, instant = { topography, ocean: { iceExchange: 'instant' } };
  const saved = await save(await steps(source, 2 * oceanEvery)), C = saved.host.state6.length, { iceWater, iceConductance, ...without } = saved.ocean;
  assert.ok(iceWater.some((x) => x > 0), 'a turbulent save that holds heat under the ice');
  const iced = [...Array(C).keys()].filter((i) => saved.host.state6[i] > 0 && iceWater[i] !== 0), seeds = thin ? iced.filter((_, n) => n % 3 === 0) : [];
  const place = (model) => {
    model.state.forEach((values, a) => values.set(saved.host[`state${a}`]));
    model.seaIce.load(model.state[6], saved.host.concentration, null);
    model.time = saved.time;
  };
  const handed = async (ocean) => {
    const gpu = await createGpuModel(new Grid(N), instant), cpu = createModel(new Grid(N), instant), warn = console.warn, { device, buffers, layout } = gpu.gpu;
    console.warn = () => {};
    try {
      place(gpu); gpu.load(); gpu.ocean.load(ocean, gpu.state[3], gpu.state[6]); loadSavedLand(gpu, saved);
      place(cpu); cpu.ocean.load(ocean, cpu.state[3], cpu.state[6], cpu.seaIce.concentration); loadSavedLand(cpu, saved);
    } finally { console.warn = warn; }
    const holding = (values) => values.reduce((n, x) => n + (x !== 0 ? 1 : 0), 0), loaded = holding(await gpu.oceanEngine.readSection('HELD'));
    const reads = [], advance = cpu.ocean.advance;
    cpu.ocean.advance = (...given) => {
      const done = advance(...given);
      reads.push({ open: iced.filter((i) => !(cpu.state[6][i] > 0)).length, arrived: iced.filter((i) => cpu.ocean.iceDrift.arrived[i]).length });
      return done;
    };
    for (let n = 0; n < oceanEvery; n++) {
      if (n === oceanEvery - 1 && seeds.length) {
        await gpu.sync();
        for (const i of seeds) { gpu.state[6][i] = 0; cpu.state[6][i] = 0; device.queue.writeBuffer(buffers.S, 4 * (layout.S.ICE + i), Float32Array.of(0)); }
      }
      await gpu.step(DT); cpu.step(DT);
    }
    await gpu.sync();
    const [flags, arrivals, after] = await Promise.all(['ICED', 'ARRIVED', 'HELD'].map((name) => gpu.oceanEngine.readSection(name)));
    const seen = { cpu: reads, gpu: { open: iced.filter((i) => !(flags[i] > 0.5)).length, arrived: iced.filter((i) => arrivals[i] > 0.5).length } };
    const [flux] = await readRanges(device, buffers.PH, [{ offset: layout.PH.OFLUX, length: C }]), area = Float64Array.from(cpu.mesh.areaCell);
    const stored = await gpu.ocean.serialize({ restart: true }), columns = Float64Array.from(stored.Qref);
    gpu.destroy();
    return { cpu: Float64Array.from(cpu.seaIce.oceanFlux), gpu: Float64Array.from(flux), area, seen, held: { loaded, read: holding(after) }, rhoCp: cpu.ocean.rhoCp, Q: { cpu: Float64Array.from(cpu.ocean.Q), gpu: columns }, h: { cpu: Float64Array.from(cpu.ocean.h), gpu: Float64Array.from(stored.h) } };
  };
  const taken = await handed(saved.ocean), plain = await handed(without), { area } = taken, interval = oceanEvery * DT;
  const held = iced.reduce((s, i) => s + area[i] * iceWater[i], 0), beyond = (engine) => iced.reduce((s, i) => s + area[i] * (taken[engine][i] - plain[engine][i]), 0) * interval;
  assert.ok(taken.held.loaded > 0 && taken.held.read === 0 && plain.held.loaded === 0 && plain.held.read === 0, `the GPU's HELD holds heat in ${taken.held.loaded} cells after the booking load and in ${taken.held.read} after the interval, against ${plain.held.loaded} and ${plain.held.read} for a load without the booking: the first read takes the held heat once`);
  if (thin) assert.deepEqual({ gpu: taken.seen.gpu, cpu: taken.seen.cpu.at(-1) }, { gpu: { open: 4, arrived: 39 }, cpu: { open: 4, arrived: 39 } }, 'at the read 4 booked cells are open water and 39 hold ice the drift brought in, on each engine');
  console.log(`${thin ? 'bare' : 'iced'} booking, ${iced.length} booked cells, ${taken.held.loaded} of them in the GPU's HELD; open and arrived at the first read: GPU ${JSON.stringify(taken.seen.gpu)}, CPU per step ${JSON.stringify(taken.seen.cpu)}`);
  const extra = { cpu: beyond('cpu'), gpu: beyond('gpu') }, worst = Math.max(...iced.map((i) => Math.abs(taken.gpu[i] - taken.cpu[i]) / Math.max(1, Math.abs(taken.cpu[i]))));
  console.log(`${iced.length} cells held ${held.toExponential(4)} J under their ice; over the first ocean interval the first write hands the ice ${extra.cpu.toExponential(4)} J beyond a load without it on the CPU, ${extra.gpu.toExponential(4)} J on the GPU; the flux to the ice differs by up to ${worst.toExponential(2)} of it`);
  assert.ok(extra.cpu > 0 && extra.cpu <= 1.0001 * held, `the CPU hands the ice ${extra.cpu} J of the ${held} J booked`);
  assert.ok(Math.abs(extra.gpu - extra.cpu) <= 1e-3 * extra.cpu && worst < 2e-3, `the GPU hands it ${extra.gpu} J against the CPU's ${extra.cpu} J, the flux apart by up to ${worst}; bounds 1e-3 and 2e-3, against 9.6e-6 and 5.3e-4 measured with every booked cell iced and 3.1e-5 and 1.82e-3 with a third bared`);
  const heatCapacity = taken.rhoCp, layers = taken.Q.cpu.length / C;
  const moved = Object.fromEntries(['cpu', 'gpu'].map((engine) => [engine, Array.from({ length: layers }, (_, k) => { let m = 0; for (let i = 0; i < C; i++) m = Math.max(m, Math.abs(taken.h[engine][k * C + i] - plain.h[engine][k * C + i])); return m; })]));
  const firstMoved = moved.gpu.findIndex((m) => m > 0), shared = firstMoved < 0 ? layers : firstMoved, cut = { cpu: layers, gpu: shared }, deeper = moved.gpu.slice(shared), thickness = Math.max(...moved.cpu);
  console.log(`the GPU's loads hold the same thickness in layers 0 to ${shared - 1}, the layers below: ${deeper.length ? `apart by up to ${Math.max(...deeper).toExponential(2)} m` : 'none'}; the CPU's loads differ by up to ${thickness.toExponential(2)} m`);
  assert.ok(shared > HYBRID_LAYERS.fixed && thickness <= 1e-3, `the GPU's loads hold the same thickness in layers 0 to ${shared - 1}, which must take in more than the ${HYBRID_LAYERS.fixed} fixed layers; the CPU's differ by up to ${thickness} m, bound 1e-3 m, about 3 times the 3.7e-4 m measured`);
  const booked = (engine) => {
    let total = 0;
    for (let i = 0; i < C; i++) {
      let column = 0;
      for (let k = 0; k < cut[engine]; k++) column += taken.Q[engine][k * C + i] - plain.Q[engine][k * C + i];
      total += area[i] * (heatCapacity * column + (taken[engine][i] - plain[engine][i]) * interval);
    }
    return total;
  };
  const closed = { cpu: (booked('cpu') - held) / held, gpu: (booked('gpu') - held) / held };
  console.log(`over all cells, all layers on the CPU and layers 0 to ${shared - 1} on the GPU, the columns' heat change plus the heat handed to the ice, beyond a load without the booking, differs from the ${held.toExponential(4)} J held by ${closed.cpu.toExponential(4)} of it on the CPU and ${closed.gpu.toExponential(4)} on the GPU`);
  assert.ok(Math.abs(closed.cpu) < bounds.cpu && Math.abs(closed.gpu) < bounds.gpu, `the booked heat closed to ${JSON.stringify(closed)} of the held heat, the CPU over all layers and the GPU over layers 0 to ${shared - 1}; ${bounds.basis}`);
};

test('a turbulent save with ice loaded into the instant model on both engines: over the first ocean interval the first write hands the ice part of the booked heat alike on both engines, beyond what a load without that heat hands it, and the columns and the ice together close on it', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  await handedHeat(false, { cpu: 1e-5, gpu: 1e-5, basis: 'bounds 1e-5 on both: 3 times the 3.08e-6 measured on the CPU, the difference between the two loads\' hybrid rescale heat changes, 0.0781447 against 0.0781416 of the held heat; 5 times the 1.91e-6 measured on the GPU' });
});

test('a turbulent save with ice loaded into the instant model on both engines, a third of its booked cells bared in the step that reads the ocean: at the read 4 booked cells are open water and 39 hold ice the drift brought in on each engine, and over the first ocean interval the first write hands the ice part of the booked heat alike on both engines, and the columns and the ice together close on it', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  await handedHeat(true, { cpu: 1e-5, gpu: 1e-5, basis: 'bounds 1e-5 on both: 4 times the 2.29e-6 measured on the CPU, the difference between the two loads\' hybrid rescale heat changes, 0.0781406 against 0.0781383 of the held heat; 15 times the 6.5e-7 measured on the GPU' });
});

test('both engines built with no ice or ocean options step the defaults alike, the conductance at the friction floor of 490.155 W/m²/K in every iced cell so that it checks the floor only: from one state, over three ocean intervals the ice, the water\'s heat under it, the conductance and the surface temperature over the sea follow the CPU', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const source = await begin({ topography }), { oceanEvery } = source.cadence;
  const saved = await save(await steps(source, 2 * oceanEvery)), C = saved.host.state6.length;
  const gpu = await createGpuModel(new Grid(N), { topography }), cpu = createModel(new Grid(N), { topography }), { device, buffers, layout } = gpu.gpu;
  assert.ok(cpu.seaIce.turbulent && cpu.seaIce.iceAlbedoPhases === DEFAULT_ICE_ALBEDO_PHASES && layout.PH.WHEAT !== undefined, 'both engines take the default exchange and phases');
  for (const model of [gpu, cpu]) {
    model.state.forEach((values, a) => values.set(saved.host[`state${a}`]));
    model.seaIce.load(model.state[6], saved.host.concentration, null);
    model.time = saved.time;
  }
  gpu.load(); gpu.ocean.load(saved.ocean, gpu.state[3], gpu.state[6]); loadSavedLand(gpu, saved);
  cpu.ocean.load(saved.ocean, cpu.state[3], cpu.state[6], cpu.seaIce.concentration); loadSavedLand(cpu, saved);
  for (let n = 0; n < 3 * oceanEvery; n++) { await gpu.step(DT); cpu.step(DT); }
  await gpu.sync();
  const [heat, conductance] = await readRanges(device, buffers.PH, ['WHEAT', 'CONDUCT'].map((name) => ({ offset: layout.PH[name], length: C })));
  const capacity = Float64Array.from(cpu.ocean.capacity), sea = [...Array(C).keys()].filter((i) => cpu.ocean.cellOcean[i]), iced = sea.filter((i) => cpu.state[6][i] > 0);
  const volume = (model, i) => (model.state[6][i] > 0 ? model.state[6][i] * model.seaIce.concentration[i] : 0);
  const apart = {
    heat: Math.max(...iced.map((i) => Math.abs(heat[i] - cpu.seaIce.waterHeat[i]) / capacity[i])),
    conductance: Math.max(...iced.map((i) => Math.abs(conductance[i] - cpu.seaIce.conductance[i]) / Math.max(1e-9, cpu.seaIce.conductance[i]))),
    volume: Math.max(...sea.map((i) => Math.abs(volume(gpu, i) - volume(cpu, i)))),
    surface: Math.max(...sea.map((i) => Math.abs(gpu.state[3][i] - cpu.state[3][i]))),
    icedSurface: Math.max(...iced.map((i) => Math.abs(gpu.state[3][i] - cpu.state[3][i]))),
  };
  console.log(`${iced.length} iced cells of ${sea.length} sea cells after ${3 * oceanEvery} steps at the defaults; the engines differ by up to ${JSON.stringify(apart)}`);
  assert.ok(iced.length > 10 && cpu.seaIce.waterHeat.some((x) => x > 0), `${iced.length} iced cells, some holding heat`);
  const floor = cpu.ocean.rhoCp * cpu.seaIce.iceExchange.heatTransfer * cpu.seaIce.iceExchange.frictionMinimum, grain = 2 ** -23;
  const off = Math.max(...iced.map((i) => Math.abs(cpu.seaIce.conductance[i] - floor) / floor));
  assert.ok(off < 1e-14, `the CPU's conductance off the floor of ${floor} W/m²/K by ${off} of it in the worst iced cell`);
  assert.ok(apart.heat < 1e-4 && apart.conductance < grain && apart.volume < 1e-5 && apart.surface < 0.025 && apart.icedSurface < 1e-3, `${JSON.stringify(apart)}; the water-heat bound 1e-4 K is 4 times the 2.55e-5 K measured; the conductance bound is the f32 grain, 2^-23 of it, the GPU holding the floor as its f32 value, 2.49e-9 off it; volume bound 1e-5 m as in the instant-save test; surface bound 0.025 K over all sea cells, 3 times the 0.0082 K measured in one ice-free tropical cell; icedSurface bound 1e-3 K over the iced cells, 4 times the 2.5e-4 K measured`);
  gpu.destroy();
});

test('both engines built with the turbulent exchange at frictionMinimum 5e-4 and the default ice and coupled drift otherwise step alike, a real share of the iced cells above the friction floor on both engines: from one state, over three ocean intervals the ice, the water\'s heat under it, the conductance and the surface temperature over the sea follow the CPU', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const options = { topography, ocean: { iceExchange: 'turbulent' } }, source = await begin(options), { oceanEvery } = source.cadence;
  const saved = await save(await steps(source, 2 * oceanEvery)), C = saved.host.state6.length;
  const gpu = await createGpuModel(new Grid(N), options), cpu = createModel(new Grid(N), options), { device, buffers, layout } = gpu.gpu;
  assert.ok(cpu.seaIce.turbulent && layout.PH.WHEAT !== undefined && cpu.seaIce.iceExchange.frictionMinimum === 5e-4, 'both engines take the turbulent exchange at frictionMinimum 5e-4');
  for (const model of [gpu, cpu]) {
    model.state.forEach((values, a) => values.set(saved.host[`state${a}`]));
    model.seaIce.load(model.state[6], saved.host.concentration, null);
    model.time = saved.time;
  }
  gpu.load(); gpu.ocean.load(saved.ocean, gpu.state[3], gpu.state[6]); loadSavedLand(gpu, saved);
  cpu.ocean.load(saved.ocean, cpu.state[3], cpu.state[6], cpu.seaIce.concentration); loadSavedLand(cpu, saved);
  for (let n = 0; n < 3 * oceanEvery; n++) { await gpu.step(DT); cpu.step(DT); }
  await gpu.sync();
  const [heat, conductance] = await readRanges(device, buffers.PH, ['WHEAT', 'CONDUCT'].map((name) => ({ offset: layout.PH[name], length: C })));
  const capacity = Float64Array.from(cpu.ocean.capacity), sea = [...Array(C).keys()].filter((i) => cpu.ocean.cellOcean[i]), iced = sea.filter((i) => cpu.state[6][i] > 0);
  const volume = (model, i) => (model.state[6][i] > 0 ? model.state[6][i] * model.seaIce.concentration[i] : 0);
  const apart = {
    heat: Math.max(...iced.map((i) => Math.abs(heat[i] - cpu.seaIce.waterHeat[i]) / capacity[i])),
    conductance: Math.max(...iced.map((i) => Math.abs(conductance[i] - cpu.seaIce.conductance[i]) / Math.max(1e-9, cpu.seaIce.conductance[i]))),
    volume: Math.max(...sea.map((i) => Math.abs(volume(gpu, i) - volume(cpu, i)))),
    surface: Math.max(...sea.map((i) => Math.abs(gpu.state[3][i] - cpu.state[3][i]))),
    icedSurface: Math.max(...iced.map((i) => Math.abs(gpu.state[3][i] - cpu.state[3][i]))),
  };
  const floor = cpu.ocean.rhoCp * cpu.seaIce.iceExchange.heatTransfer * cpu.seaIce.iceExchange.frictionMinimum, above = (values) => iced.filter((i) => values[i] > 1.05 * floor).length / iced.length;
  const share = { cpu: above(cpu.seaIce.conductance), gpu: above(conductance) }, driven = iced.filter((i) => cpu.seaIce.conductance[i] > 1.05 * floor && conductance[i] > 1.05 * floor);
  console.log(`${iced.length} iced cells of ${sea.length} sea cells after ${3 * oceanEvery} steps at frictionMinimum 5e-4; the share of them above the friction floor ${JSON.stringify(share)}, ${driven.length} above it on both; the engines differ by up to ${JSON.stringify(apart)}`);
  assert.ok(iced.length > 10 && share.cpu >= 0.2 && share.gpu >= 0.2, `${iced.length} iced cells, the share above the floor ${JSON.stringify(share)}`);
  assert.ok(driven.length >= 0.2 * iced.length, `${driven.length} of ${iced.length} iced cells above the floor on both engines`);
  assert.ok(apart.heat < 1e-4 && apart.conductance < 2e-4 && apart.volume < 1e-5 && apart.surface < 0.025 && apart.icedSurface < 1e-3, `${JSON.stringify(apart)}; the water-heat bound 1e-4 K is 4 times the 2.51e-5 K measured; the surface bound 0.025 K over all sea cells is 3 times the 0.0082 K measured in one ice-free tropical cell; the icedSurface bound 1e-3 K is about 6 times the 1.76e-4 K measured over the iced cells; the conductance bound 2e-4 is about three times the 5.95e-5 measured over the drag and speed branch, which the f32 relative speed makes looser than the floor's f32 grain of 2^-23`);
  gpu.destroy();
});

test('both engines built with the default exchange at leadTemperature water step alike and apart from the default: from one state, over three ocean intervals the ice, the water\'s heat under it and the surface temperature follow the CPU, and the held heat parts from a CPU run at the default by far more than the engines part', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const options = { topography, ocean: { iceExchange: { ...DEFAULT_ICE_EXCHANGE, leadTemperature: 'water' } } }, source = await begin(options), { oceanEvery } = source.cadence;
  const saved = await save(await steps(source, 2 * oceanEvery)), C = saved.host.state6.length;
  const gpu = await createGpuModel(new Grid(N), options), cpu = createModel(new Grid(N), options), base = createModel(new Grid(N), { topography }), { device, buffers, layout } = gpu.gpu;
  assert.ok(cpu.seaIce.iceExchange.leadTemperature === 'water' && base.seaIce.iceExchange.leadTemperature === 'freezing' && layout.PH.WHEAT !== undefined, 'the CPU takes leadTemperature water, the base the default');
  for (const model of [gpu, cpu, base]) {
    model.state.forEach((values, a) => values.set(saved.host[`state${a}`]));
    model.seaIce.load(model.state[6], saved.host.concentration, null);
    model.time = saved.time;
  }
  gpu.load(); gpu.ocean.load(saved.ocean, gpu.state[3], gpu.state[6]); loadSavedLand(gpu, saved);
  for (const model of [cpu, base]) { model.ocean.load(saved.ocean, model.state[3], model.state[6], model.seaIce.concentration); loadSavedLand(model, saved); }
  for (let n = 0; n < 3 * oceanEvery; n++) { await gpu.step(DT); cpu.step(DT); base.step(DT); }
  await gpu.sync();
  const [heat] = await readRanges(device, buffers.PH, [{ offset: layout.PH.WHEAT, length: C }]);
  const capacity = Float64Array.from(cpu.ocean.capacity), sea = [...Array(C).keys()].filter((i) => cpu.ocean.cellOcean[i]), iced = sea.filter((i) => cpu.state[6][i] > 0 && base.state[6][i] > 0);
  const volume = (model, i) => (model.state[6][i] > 0 ? model.state[6][i] * model.seaIce.concentration[i] : 0);
  const apart = {
    heat: Math.max(...iced.map((i) => Math.abs(heat[i] - cpu.seaIce.waterHeat[i]) / capacity[i])),
    volume: Math.max(...sea.map((i) => Math.abs(volume(gpu, i) - volume(cpu, i)))),
    surface: Math.max(...sea.map((i) => Math.abs(gpu.state[3][i] - cpu.state[3][i]))),
    icedSurface: Math.max(...iced.map((i) => Math.abs(gpu.state[3][i] - cpu.state[3][i]))),
  };
  const option = Math.max(...iced.map((i) => Math.abs(cpu.seaIce.waterHeat[i] - base.seaIce.waterHeat[i]) / capacity[i]));
  console.log(`${iced.length} iced cells of ${sea.length} sea cells after ${3 * oceanEvery} steps at leadTemperature water; the engines differ by up to ${JSON.stringify(apart)}; the held heat parts from the default's by up to ${option} K`);
  assert.ok(iced.length > 10 && cpu.seaIce.waterHeat.some((x) => x > 0), `${iced.length} iced cells, some holding heat`);
  assert.ok(option > 10 * apart.heat, `the option moves the held heat by ${option} K, the engines part by ${apart.heat} K; measured 0.029 K against 2.75e-5 K`);
  assert.ok(apart.heat < 1e-4 && apart.volume < 1e-5 && apart.surface < 0.025 && apart.icedSurface < 1e-3, `${JSON.stringify(apart)}; the bounds of the default case, 4, 7, 3 and 3 times the 2.75e-5 K, 1.46e-6 m, 0.0082 K and 2.9e-4 K measured`);
  gpu.destroy();
});

const twins = async (options, prepare = () => {}) => {
  const source = await begin(options), { oceanEvery } = source.cadence;
  const saved = await save(await steps(source, 2 * oceanEvery));
  prepare(saved);
  const gpu = await createGpuModel(new Grid(N), options), cpu = createModel(new Grid(N), options), base = createModel(new Grid(N), { topography });
  for (const model of [gpu, cpu, base]) {
    model.state.forEach((values, a) => values.set(saved.host[`state${a}`]));
    model.seaIce.load(model.state[6], saved.host.concentration, null);
    model.time = saved.time;
  }
  gpu.load(); gpu.ocean.load(saved.ocean, gpu.state[3], gpu.state[6]); loadSavedLand(gpu, saved);
  for (const model of [cpu, base]) { model.ocean.load(saved.ocean, model.state[3], model.state[6], model.seaIce.concentration); loadSavedLand(model, saved); }
  for (let n = 0; n < 3 * oceanEvery; n++) { await gpu.step(DT); cpu.step(DT); base.step(DT); }
  await gpu.sync();
  const C = saved.host.state6.length, { device, buffers, layout } = gpu.gpu, [heat] = await readRanges(device, buffers.PH, [{ offset: layout.PH.WHEAT, length: C }]);
  const capacity = Float64Array.from(cpu.ocean.capacity), sea = [...Array(C).keys()].filter((i) => cpu.ocean.cellOcean[i]);
  const volume = (model, i) => (model.state[6][i] > 0 ? model.state[6][i] * model.seaIce.concentration[i] : 0);
  const apart = (cells) => ({
    heat: Math.max(0, ...cells.filter((i) => cpu.state[6][i] > 0 && gpu.state[6][i] > 0).map((i) => Math.abs(heat[i] - cpu.seaIce.waterHeat[i]) / capacity[i])),
    volume: Math.max(...sea.map((i) => Math.abs(volume(gpu, i) - volume(cpu, i)))),
    surface: Math.max(...sea.map((i) => Math.abs(gpu.state[3][i] - cpu.state[3][i]))),
  });
  return { gpu, cpu, base, saved, oceanEvery, heat, capacity, sea, volume, apart };
};

test('both engines built with the default exchange at traceMelt step alike and apart from the default: cells cut to traces of cover 0.005 over held heat melt whole on both and stay on the default, and the ice, the held heat and the surface follow the CPU', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const options = { topography, ocean: { iceExchange: { ...DEFAULT_ICE_EXCHANGE, traceMelt: true } } };
  let traces = [];
  const run = await twins(options, (saved) => {
    traces = [...saved.host.state6.keys()].filter((i) => saved.host.state6[i] > 0 && saved.ocean.iceWater[i] > 0);
    for (const i of traces) saved.host.concentration[i] = 0.005;
  });
  const { gpu, cpu, base, oceanEvery, volume } = run, gone = (model) => traces.filter((i) => !(model.state[6][i] > 0)).length, held = (model) => traces.reduce((sum, i) => sum + volume(model, i), 0);
  const apart = run.apart(run.sea), counts = { traces: traces.length, cpu: gone(cpu), gpu: gone(gpu), base: gone(base) }, melted = { option: held(base) - held(cpu), engines: Math.abs(held(gpu) - held(cpu)) };
  console.log(`traceMelt: ${JSON.stringify(counts)} traces gone after ${3 * oceanEvery} steps, their ice ${JSON.stringify(melted)} m less summed over them; the engines differ by up to ${JSON.stringify(apart)}`);
  assert.ok(cpu.seaIce.iceExchange.traceMelt && !base.seaIce.iceExchange.traceMelt, 'the CPU takes traceMelt, the base the default');
  assert.ok(counts.traces > 10 && counts.cpu > 0 && counts.base === 0 && counts.gpu === counts.cpu, `the same traces go on both engines and none on the default: ${JSON.stringify(counts)}`);
  assert.ok(melted.option > 10 * melted.engines, `the option melts ${melted.option} m of the traces, the engines part by ${melted.engines} m`);
  assert.ok(apart.heat < 1e-4 && apart.volume < 2e-4 && apart.surface < 0.08, `${JSON.stringify(apart)}; 11, 3.3 and 3.1 times the 9.4e-6 K, 6.1e-5 m and 0.0256 K measured`);
  gpu.destroy();
});

test('both engines built with leadPenetration step alike and apart from the default: the leads\' sunlight that goes below the top layer leaves the water under the ice cooler than the default\'s by far more than the engines part, and the ice and the surface follow the CPU', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const options = { topography, ice: { leadPenetration: true } };
  const run = await twins(options), { gpu, cpu, base, oceanEvery, capacity } = run;
  const iced = run.sea.filter((i) => cpu.state[6][i] > 0 && base.state[6][i] > 0 && cpu.seaIce.concentration[i] < 1), apart = run.apart(iced);
  const option = Math.max(...iced.map((i) => (base.seaIce.waterHeat[i] - cpu.seaIce.waterHeat[i]) / capacity[i]));
  console.log(`leadPenetration: ${iced.length} iced cells with leads after ${3 * oceanEvery} steps; the engines differ by up to ${JSON.stringify(apart)}; the held heat is below the default's by up to ${option} K`);
  assert.ok(iced.length > 10, `${iced.length} cells with leads`);
  const mean = iced.reduce((sum, i) => sum + (base.seaIce.waterHeat[i] - cpu.seaIce.waterHeat[i]) / capacity[i], 0) / iced.length;
  assert.ok(mean > 0, `the water under the ice is ${mean} K cooler than the default's on average`);
  assert.ok(option > 10 * apart.heat, `the option takes up to ${option} K, the engines part by ${apart.heat} K`);
  assert.ok(apart.heat < 1e-4 && apart.volume < 1e-5 && apart.surface < 0.12, `${JSON.stringify(apart)}; 3.3, 5.1 and 3.2 times the 3.0e-5 K, 1.95e-6 m and 0.0379 K measured`);
  gpu.destroy();
});

test('both engines built with traceMelt and leadPenetration together step alike and apart from the default: cells cut to traces of cover 0.005 over held heat, beside iced cells that drift feeds them from, lose far more ice than on the default and alike on both engines, the water under the other iced cells is cooler than the default\'s, and the ice, the held heat and the surface follow the CPU', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const options = { topography, ice: { leadPenetration: true }, ocean: { iceExchange: { ...DEFAULT_ICE_EXCHANGE, traceMelt: true } } };
  let traces = [];
  const run = await twins(options, (saved) => {
    traces = [...saved.host.state6.keys()].filter((i) => i % 2 === 0 && saved.host.state6[i] > 0 && saved.ocean.iceWater[i] > 0);
    for (const i of traces) saved.host.concentration[i] = 0.005;
  });
  const { gpu, cpu, base, oceanEvery, volume, capacity } = run, gone = (model) => traces.filter((i) => !(model.state[6][i] > 0)).length, held = (model) => traces.reduce((sum, i) => sum + volume(model, i), 0);
  const cut = new Set(traces), iced = run.sea.filter((i) => !cut.has(i) && cpu.state[6][i] > 0 && base.state[6][i] > 0 && cpu.seaIce.concentration[i] < 1 && Math.min(cpu.seaIce.concentration[i], base.seaIce.concentration[i]) >= 0.01);
  const apart = run.apart(run.sea), counts = { traces: traces.length, iced: iced.length, cpu: gone(cpu), gpu: gone(gpu), base: gone(base) }, melted = { option: held(base) - held(cpu), engines: Math.abs(held(gpu) - held(cpu)) };
  const cooler = iced.reduce((sum, i) => sum + (base.seaIce.waterHeat[i] - cpu.seaIce.waterHeat[i]) / capacity[i], 0) / iced.length;
  console.log(`traceMelt with leadPenetration: ${JSON.stringify(counts)} traces gone after ${3 * oceanEvery} steps, their ice ${JSON.stringify(melted)} m less summed over them; the other iced cells' water ${cooler} K cooler than the default's on average; the engines differ by up to ${JSON.stringify(apart)}`);
  assert.ok(cpu.seaIce.iceExchange.traceMelt && !base.seaIce.iceExchange.traceMelt, 'the CPU takes traceMelt, the base the default');
  assert.ok(counts.traces > 10 && counts.base === 0 && counts.gpu === counts.cpu, `the same traces go on both engines and none on the default: ${JSON.stringify(counts)}`);
  assert.ok(melted.option > 10 * melted.engines, `the options melt ${melted.option} m of the traces, the engines part by ${melted.engines} m`);
  assert.ok(iced.length > 10 && cooler > 0, `${iced.length} other iced cells, their water ${cooler} K cooler than the default's`);
  assert.ok(apart.heat < 1e-4 && apart.volume < 2e-5 && apart.surface < 0.025, `${JSON.stringify(apart)}; 12, 5.8 and 3.9 times the 8.1e-6 K, 3.5e-6 m and 0.0064 K measured`);
  gpu.destroy();
});

test('an instant save with ice loaded into the turbulent model on both engines: the water under the ice starts within the f32 grain of the freezing point, as on the CPU, and over two ocean intervals the ice, the held heat and the conductance follow the CPU and no ice collapses', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const options = { topography, ocean: { iceExchange: TURBULENT } }, source = await begin({ topography, ocean: { iceExchange: 'instant' } }), { oceanEvery } = source.cadence;
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
