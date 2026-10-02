import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createGeography, syntheticTopography } from '../js/geography.module.js';
import { createLandSurface, decompositionWarmth, decompositionMoisture, litterWarmth, dryHumusAlbedo, carbonEquilibrium, airCycle, seasonEstimate, sineSeason, SOIL_CARBON } from '../js/physics/land.module.js';
import { MELTING_POINT, SNOW_AGEING } from '../js/physics/ice.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { encodeState, decodeState } from '../js/stateFile.module.js';
import { regridLand } from '../js/physics/regrid.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel, VEGETATION_OPTIONS } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const YEAR = 365 * 86400;
const mesh = createModel(new Grid(8), { physics: false }).mesh;
const flat = () => createGeography(mesh, syntheticTopography(90, 180, () => 100), { landBridges: {}, seaStraits: {} });
const near = (a, b, tol = 1e-12) => Math.abs(a - b) <= tol * Math.max(1, Math.abs(b));
const frozen = { growthTime: 1e30, declineTime: 1e30, snowDeclineTime: 1e30, treeGrowthTime: 1e30, treeDeclineTime: 1e30 };

test('decomposition follows Lloyd and Taylor in the air and TRIFFID\'s moisture factor in the fill, litter Lieth\'s warmth in the season, and the dry soil darkens exponentially with its carbon', () => {
  assert.ok(near(decompositionWarmth(283.15), 1) && near(decompositionWarmth(273.15), 0.3021359906222061) && near(decompositionWarmth(299.15), 3.3996326829918186));
  assert.equal(decompositionWarmth(220), 0);
  assert.equal(decompositionMoisture(0.05), 0.2);
  assert.ok(near(decompositionMoisture(0.325), 0.6) && near(decompositionMoisture(0.55), 1) && near(decompositionMoisture(1), 0.64));
  assert.equal(decompositionMoisture(0.55, true), 0.2, 'frozen');
  assert.ok(near(litterWarmth(MELTING_POINT + 25), 0.8402380030563309));
  assert.equal(litterWarmth(MELTING_POINT + 0.5), 0, 'outside the growing season');
  assert.ok(near(dryHumusAlbedo(0), 0.37) && near(dryHumusAlbedo(2.6), 0.12 + 0.25 / Math.E) && near(dryHumusAlbedo(13), 0.12 + 0.25 * Math.exp(-5)));
});

function cell(options = {}) {
  const geography = flat(), land = createLandSurface(mesh, geography, { ...frozen, ...options });
  land.initialize();
  const i = [...geography.land.keys()].find((n) => geography.land[n] && !geography.iceSheet[n]);
  const surfaceT = new Float64Array(mesh.nCells).fill(290), flux = new Float64Array(mesh.nCells);
  const set = (fill, cover, trees, carbon) => { land.soil[i] = fill * land.capacity(i); land.surface[i] = 0; land.snow[i] = 0; land.vegetation[i] = cover; land.canopy[i] = trees; land.soilCarbon[i] = carbon; };
  const run = (celsius, years) => land.update(i, surfaceT, flux, 0, years * YEAR, MELTING_POINT + celsius);
  return { land, i, set, run };
}

test('a warm wet forest, a hot dry desert, a cold wet tundra and a cell losing its cover move their carbon as hand-computed', () => {
  const { land, i, set, run } = cell();
  set(0.55, 1, 1, 0);
  run(26, 0.2);
  assert.ok(near(land.soilCarbon[i], 4.014013445494061, 1e-10), `forest after 0.2 years at 26 °C: ${land.soilCarbon[i]}`);
  run(26, 10);
  assert.ok(near(land.soilCarbon[i], 6.459437778118126, 1e-10) && near(land.albedo(i), 0.13), `forest equilibrium ${land.soilCarbon[i]}`);
  assert.ok(near(land.dryAlbedo(i), 0.14084390902416827, 1e-10));
  set(0.12, 0.1, 0, 0);
  run(30, 30);
  assert.ok(near(land.soilCarbon[i], 0.5051792913500829, 1e-10) && near(land.dryAlbedo(i), 0.32585276709628985, 1e-10), `desert ${land.soilCarbon[i]}`);
  set(0.75, 0.4, 0, 0);
  run(8, 100);
  assert.ok(near(land.soilCarbon[i], 8.383854272122287, 1e-10) && near(land.dryAlbedo(i), 0.12994332610027737, 1e-10), `tundra summer ${land.soilCarbon[i]}`);
  set(0.75, 0.4, 0, 10);
  run(-5, 0.5);
  assert.ok(near(land.soilCarbon[i], 9.811185899105775, 1e-10), `tundra winter: no litter, frozen decay ${land.soilCarbon[i]}`);
  set(0.55, 0, 0, 6.459437778118126);
  run(20, 0.5);
  assert.ok(near(land.soilCarbon[i], 1.2465783320652897, 1e-10) && near(land.dryAlbedo(i), 0.2747804580534271, 1e-10), `losing its cover ${land.soilCarbon[i]}`);
  const plain = cell({ soilCarbon: false });
  plain.set(0.55, 0, 0, 3);
  plain.run(20, 1);
  assert.equal(plain.land.soilCarbon[plain.i], 3, 'without the scheme the store stands');
  assert.equal(plain.land.dryAlbedo(plain.i), 0.30);
});

test('the acceleration moves the store toward the same equilibrium, a sine year\'s stepped mean is carbonEquilibrium\'s, and the wet soil keeps half the dry albedo', () => {
  const slow = cell({ carbonAcceleration: 1 }), fast = cell();
  slow.set(0.55, 1, 1, 0); fast.set(0.55, 1, 1, 0);
  slow.run(26, 100); fast.run(26, 1);
  assert.ok(near(slow.land.soilCarbon[slow.i], fast.land.soilCarbon[fast.i], 1e-12), 'A × time is what counts');
  for (const [mean, amplitude, fill, cover] of [[-11, 16, 0.75, 0.4], [10, 11, 0.6, 0.9], [22, 8, 0.12, 0.1]]) {
    const { land, i, set } = cell({ carbonAcceleration: 1 });
    const expected = carbonEquilibrium(mean, amplitude, fill, cover);
    set(fill, cover, 0, expected);
    const surfaceT = new Float64Array(mesh.nCells).fill(290), flux = new Float64Array(mesh.nCells), steps = 365 * 4;
    let sum = 0;
    for (let n = 0; n < steps; n++) { land.update(i, surfaceT, flux, 0, YEAR / steps, MELTING_POINT + mean + amplitude * Math.cos(2 * Math.PI * (n + 0.5) / steps)); sum += land.soilCarbon[i] / steps; }
    assert.ok(Math.abs(sum / expected - 1) < 2e-3, `${mean} ± ${amplitude} °C: stepped ${sum}, equilibrium ${expected}`);
  }
  const { land, i, set } = cell();
  set(0.55, 0, 0, 2.6);
  land.surface[i] = 15;
  assert.ok(near(land.albedo(i), 0.5 * dryHumusAlbedo(2.6)), `the full surface layer halves the dry albedo: ${land.albedo(i)}`);
});

test('a fresh start and an older state start at the equilibrium of their own cover and the estimated year, a saved state keeps its carbon, the ice sheets hold none, and regridding carries it', async () => {
  const topography = syntheticTopography(180, 360, (lat, lon) => (lat < -1.2 ? 2000 : Math.cos(lon) > 0 ? 100 : -4000));
  const geography = createGeography(mesh, topography, { landBridges: {}, seaStraits: {} });
  const land = createLandSurface(mesh, geography), C = mesh.nCells;
  land.initialize();
  let checked = 0;
  for (let i = 0; i < C; i++) {
    if (!geography.land[i] || geography.iceSheet[i]) { assert.equal(land.soilCarbon[i], 0); continue; }
    const { mean, amplitude } = airCycle(mesh.latCell[i], geography.elevation[i]);
    assert.ok(near(land.soilCarbon[i], carbonEquilibrium(mean, amplitude, 0.5, land.canopy[i] + Math.max(0, 0.5 - land.canopy[i]))), `cell ${i}`);
    checked++;
  }
  assert.ok(checked > 50);
  const e = seasonEstimate(0.7, 300), c = airCycle(0.7, 300);
  assert.deepEqual(e, sineSeason(c.mean, c.amplitude, 0.9), 'the season estimate reads the same year');
  const vegetation = Float64Array.from({ length: C }, (_, i) => (i % 5) / 5), soil = Float64Array.from({ length: C }, (_, i) => 30 * (i % 10));
  const older = { soil, snow: new Float64Array(C), vegetation };
  land.load(older);
  for (let i = 0; i < C; i++) {
    if (!geography.land[i] || geography.iceSheet[i]) { assert.equal(land.soilCarbon[i], 0); continue; }
    const { mean, amplitude } = airCycle(mesh.latCell[i], geography.elevation[i]);
    const litter = Math.min(1, land.canopy[i]) + Math.max(0, land.vegetation[i] - land.canopy[i]);
    assert.ok(near(land.soilCarbon[i], carbonEquilibrium(mean, amplitude, soil[i] / 300, litter)), `older cell ${i}`);
  }
  const soilCarbon = Float64Array.from({ length: C }, (_, i) => (i % 9) * 0.7);
  land.load({ ...older, soilCarbon });
  const saved = await decodeState(encodeState({ N: 8, K: 1, day: 0, time: 0, terrain: true, land: land.serialize() }));
  const back = createLandSurface(mesh, geography);
  back.load(saved.land);
  for (let i = 0; i < C; i++) assert.ok(near(back.soilCarbon[i], geography.land[i] && !geography.iceSheet[i] ? Math.fround(soilCarbon[i]) : 0, 1e-6), `cell ${i}: kept ${back.soilCarbon[i]}`);
  const source = createModel(new Grid(8), { physics: false, topography }), target = createModel(new Grid(10), { topography });
  const carried = regridLand(source, target, { soil: new Float64Array(C).fill(150), snow: new Float64Array(C), vegetation: new Float64Array(C).fill(0.5), soilCarbon: new Float64Array(C).fill(4.2) }, null, { surfaceT: new Float64Array(C).fill(290) });
  let kept = 0, estimated = 0;
  for (let n = 0; n < target.mesh.nCells; n++) {
    if (!target.geography.land[n]) { assert.equal(carried.soilCarbon[n], 0); continue; }
    if (carried.soilCarbon[n] === 4.2) kept++; else { estimated++; assert.ok(carried.soilCarbon[n] >= 0); }
  }
  assert.ok(kept > 50, `${kept} carried, ${estimated} estimated`);
  assert.equal(regridLand(source, target, { soil: new Float64Array(C), snow: new Float64Array(C) }).soilCarbon, undefined);
});

test('every option of the land surface reaches the GPU or is refused there', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const source = createLandSurface.toString();
  const block = source.slice(source.indexOf('{', source.indexOf('geography,')) + 1, source.indexOf('} = {})'));
  const names = [...block.matchAll(/(?:^|,)\s*([A-Za-z]+)\s*(?::|=)/g)].map((m) => m[1]);
  const handled = new Set([...VEGETATION_OPTIONS, 'snowAgeing', ...Object.keys(SNOW_AGEING), 'heatCapacity', 'bucketCapacity', 'wetnessThreshold', 'albedo', 'snowAlbedo', 'fullSnow', 'latentHeatFusion', 'buffers']);
  assert.ok(names.length > 60, names.join(' '));
  assert.deepEqual(names.filter((n) => !handled.has(n)), []);
  for (const key of Object.keys(SOIL_CARBON)) assert.ok(VEGETATION_OPTIONS.includes(key), key);
  const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 ? 300 : -4000));
  await assert.rejects(createGpuModel(new Grid(4), { topography, land: { latentHeatFusion: 3e5 } }), /latentHeatFusion/);
});

const topography = syntheticTopography(90, 180, (lat, lon) => ((Math.cos(lon) > 0 && Math.abs(lat) < 1.2) || lat < -1.15 ? 300 : -4000));

test('both engines darken the bare soil alike with its carbon, wet and dry, and step the carbon alike', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const land = { carbonAcceleration: 3e5, growthTime: 4 * 3600, declineTime: 3 * 3600, snowDeclineTime: 4 * 3600 };
  const prepare = (model, s0) => {
    const C = model.mesh.nCells, init = initializeState(model, {});
    for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
    for (let i = 0; i < C; i++) if (model.geography.land[i]) { model.state[6][i] = 0; init[3][i] = MELTING_POINT - 10 + 35 * ((i * 7919) % 101) / 101; model.state[3][i] = init[3][i]; }
    model.seaIce.load(model.state[6], null);
    if (model.load) model.load();
    model.ocean.initialize(model.state[3], model.state[6]);
    let s = s0;
    const rnd = () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; };
    const snow = Float64Array.from({ length: C }, (_, i) => (Math.abs(model.mesh.latCell[i]) > 1.0 ? 30 : 0));
    const vegetation = Float64Array.from({ length: C }, () => rnd()), canopy = Float64Array.from(vegetation, (v) => v * rnd());
    const soilCarbon = Float64Array.from({ length: C }, () => 8 * rnd() * rnd());
    model.land.load({ soil: Float64Array.from({ length: C }, () => 300 * rnd()), surface: Float64Array.from({ length: C }, () => 15 * rnd() * rnd()), snow, vegetation, canopy, soilCarbon }, model.state[6]);
    return model;
  };
  const cpu = prepare(createModel(new Grid(6), { topography, land }), 11);
  const gpu = prepare(await createGpuModel(new Grid(6), { topography, land }), 11);
  const C = cpu.mesh.nCells, mask = cpu.geography.land, sheet = cpu.geography.iceSheet;
  const albedo = Float64Array.from({ length: C }, (_, i) => (mask[i] ? cpu.land.albedo(i) : NaN)), before = Float64Array.from(cpu.land.soilCarbon);
  const plain = createLandSurface(cpu.mesh, cpu.geography, { soilCarbon: false });
  plain.load({ soil: cpu.land.soil, surface: cpu.land.surface, snow: cpu.land.snow, vegetation: cpu.land.vegetation, canopy: cpu.land.canopy }, cpu.state[6]);
  await gpu.step(900);
  const ph = await gpu.gpu.downloadPhysics();
  let worstAlbedo = 0, moved = 0, cells = 0;
  for (let i = 0; i < C; i++) {
    if (!mask[i] || sheet[i]) continue;
    cells++;
    worstAlbedo = Math.max(worstAlbedo, Math.abs(ph.ADIF[i] - albedo[i]));
    moved = Math.max(moved, Math.abs(albedo[i] - plain.albedo(i)));
  }
  cpu.step(900);
  for (let n = 1; n < 24; n++) { cpu.step(900); await gpu.step(900); }
  await gpu.sync();
  const saved = await gpu.land.serialize();
  let worstCarbon = 0, scale = 0, rose = 0, fell = 0;
  for (let i = 0; i < C; i++) {
    if (!mask[i] || sheet[i]) continue;
    worstCarbon = Math.max(worstCarbon, Math.abs(cpu.land.soilCarbon[i] - saved.soilCarbon[i]));
    scale = Math.max(scale, Math.abs(cpu.land.soilCarbon[i] - before[i]));
    if (cpu.land.soilCarbon[i] > before[i] + 0.1) rose++;
    if (cpu.land.soilCarbon[i] < before[i] - 0.1) fell++;
  }
  gpu.destroy();
  console.log(`N=6, ${cells} land cells: albedo engines apart by ${worstAlbedo.toExponential(1)} (the carbon moves it by up to ${moved.toFixed(3)}); after 24 steps at A = 3e5 the carbon rose on ${rose} and fell on ${fell}, by up to ${scale.toFixed(2)} kg/m², engines apart by ${worstCarbon.toExponential(1)} kg/m²`);
  assert.ok(worstAlbedo < 1e-6 && moved > 0.1, `albedo ${worstAlbedo}, moved ${moved}`);
  assert.ok(rose > 10 && fell > 10 && worstCarbon < 1e-3 * scale, `rose ${rose}, fell ${fell}, carbon ${worstCarbon} of ${scale}`);
});
