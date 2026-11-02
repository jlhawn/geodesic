import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createModel } from '../js/model.module.js';
import { createGeography, syntheticTopography, readTopography } from '../js/geography.module.js';
import { createLandSurface } from '../js/physics/land.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { encodeState, decodeState } from '../js/stateFile.module.js';
import { options, digest, landArrays } from './helpers/savedLand.mjs';

const N = 16, mesh = buildMesh(new Grid(N)), C = mesh.nCells;
const flat = createGeography(mesh, syntheticTopography(90, 180, () => 100), { landBridges: {}, seaStraits: {} });
const geo16 = createGeography(mesh, readTopography('etopo1'));
const fresh = (geography = flat, extra = {}) => { const land = createLandSurface(mesh, geography, { start: 'neutral', deepCapacity: 300, ...extra }); land.initialize(); return land; };
const stores = (land, i) => land.soil[i] + land.surface[i] + land.snow[i] + land.deepWater[i];

const PLAIN_LAND = '8f93439d5d613b0832b50350fe8b566bd1dfedfbd061d5b83fbdf1b35285bd3a';
const PLAIN_STEPPED = '8da09236b850d009fb00ee3b7ac018decef154be168a1abfab49120f5c6dca71';

test('the option is refused below 0, without vegetation, and off by default, where the land has no deep store', () => {
  assert.throws(() => createLandSurface(mesh, flat, { deepCapacity: -1 }), /deepCapacity/);
  assert.throws(() => createLandSurface(mesh, flat, { deepCapacity: 300, vegetation: false }), /needs vegetation/);
  for (const extra of [{}, { deepCapacity: 0 }]) {
    const land = createLandSurface(mesh, flat, { start: 'neutral', ...extra });
    land.initialize();
    assert.equal(land.deepWater, null);
    assert.ok(!('deep' in land.serialize()) && !('deepWater' in land.shared));
  }
  const land = fresh();
  assert.ok(land.deepWater instanceof Float64Array && land.shared.deepWater instanceof SharedArrayBuffer);
  for (let i = 0; i < C; i++) assert.equal(land.deepWater[i], flat.land[i] && !(flat.iceSheet && flat.iceSheet[i] > 0) ? 150 : 0, `cell ${i}: a fresh start fills the store to the root zone's fill, off the ice sheets`);
});

test('what the root zone sheds fills the deep store before any of it runs off, and only the overflow of a full store is runoff', () => {
  const land = fresh(), i = 0, surfaceT = new Float64Array(C).fill(290), flux = new Float64Array(C);
  land.soil[i] = 300; land.surface[i] = 15; land.deepWater[i] = 270; land.runoff[i] = 0;
  land.deposit(i, 20, 290);
  assert.deepEqual([land.soil[i], land.surface[i], land.deepWater[i], land.runoff[i]], [300, 15, 290, 0], 'a full root zone sheds all its infiltration into the store');
  land.deposit(i, 20, 290);
  assert.deepEqual([land.soil[i], land.deepWater[i], land.runoff[i]], [300, 300, 10], 'the store fills, and the rest runs off');
  land.deepWater[i] = 280; land.soil[i] = 320; land.surface[i] = 0;
  land.update(i, surfaceT, flux, 0, 900);
  assert.deepEqual([land.soil[i], land.deepWater[i], land.runoff[i]], [300, 300, 10], 'the water above the root zone\'s capacity recharges the store');
  land.soil[i] = 330; land.surface[i] = 0;
  land.update(i, surfaceT, flux, 0, 900);
  assert.deepEqual([land.soil[i], land.deepWater[i], land.runoff[i]], [300, 300, 40], 'and runs off once it is full');
  const j = 1;
  land.soil[j] = 270; land.surface[j] = 15; land.deepWater[j] = 0; land.runoff[j] = 0;
  let before = stores(land, j), rain = 0;
  while (land.deepWater[j] < 300) {
    const shed = 30 * Math.pow(Math.min(1, land.soil[j] / 300), 4);
    land.deposit(j, 30, 290); rain += 30;
    if (land.deepWater[j] < 300) assert.equal(land.runoff[j], 0, `nothing runs off while the store has room (${land.deepWater[j]} kg/m²)`);
    assert.ok(shed > 0 && Math.abs(stores(land, j) + land.runoff[j] - before - rain) < 1e-9, 'every drop is in a store or the runoff');
  }
  assert.ok(land.runoff[j] > 0 && land.runoff[j] < 30, `the runoff ${land.runoff[j]} starts with the rain that fills the store`);
  const iced = Array.from({ length: C }, (_, n) => n).filter((n) => geo16.land[n] && geo16.iceSheet && geo16.iceSheet[n] > 0);
  assert.ok(iced.length > 0, 'ice-sheet cells at N = 16');
  const sheet = fresh(geo16), k = iced[0];
  assert.equal(sheet.deepWater[k], 0, 'an ice sheet has no deep store');
  sheet.soil[k] = 300; sheet.surface[k] = 15; sheet.runoff[k] = 0;
  sheet.deposit(k, 20, 290);
  assert.deepEqual([sheet.deepWater[k], sheet.runoff[k]], [0, 20], 'and its shed runs off');
});

test('the trees draw from the wetter store in relative terms, the root zone on a tie, and the grass and the bare ground from their own stores', () => {
  const land = fresh(), i = 0, aero = 0.01, resist = 1 + 70 * aero, dt = 3600, evaporation = 1e-4, drawn = evaporation * dt;
  const surfaceT = new Float64Array(C).fill(295), flux = new Float64Array(C);
  const set = (soil, deep, cover, trees, surface = 0) => { land.soil[i] = soil; land.deepWater[i] = deep; land.vegetation[i] = cover; land.canopy[i] = trees; land.surface[i] = surface; land.snow[i] = 0; };
  const step = () => { const before = { soil: land.soil[i], deep: land.deepWater[i], surface: land.surface[i] }; land.update(i, surfaceT, flux, evaporation, dt); return before; };
  for (const [deep, fromDeep] of [[180, true], [18, false], [45, false]]) {
    set(45, deep, 0.8, 0.4);
    const roots = 45 / 225, factor = deep / 225, drawing = 0.5 * roots + 0.5 * Math.max(roots, factor);
    assert.ok(Math.abs(land.wetness(i, aero, 295) - 0.8 * drawing / resist) < 1e-12, `deep ${deep}: the trees' factor is the larger`);
    const trees = fromDeep ? 0.5 * factor / drawing : 0;
    const before = step();
    assert.ok(Math.abs(land.deepWater[i] - (before.deep - drawn * trees)) < 1e-12, `deep ${deep}: the store gives ${before.deep - land.deepWater[i]}`);
    assert.ok(Math.abs(land.soil[i] - (before.soil - drawn * (1 - trees))) < 1e-12, `deep ${deep}: the root zone gives ${before.soil - land.soil[i]}`);
  }
  set(45, 300, 0.8, 0);
  assert.ok(Math.abs(land.wetness(i, aero, 295) - 0.8 * 0.2 / resist) < 1e-12, 'grass alone keeps the root zone\'s factor');
  let before = step();
  assert.equal(land.deepWater[i], 300, 'and leaves the store');
  assert.ok(Math.abs(land.soil[i] - (before.soil - drawn)) < 1e-12);
  set(45, 180, 0.8, 0.4, 10);
  const bare = 0.2 * 10 / 15, grass = 0.8 * 0.5 * 0.2 / resist, trees = 0.8 * 0.5 * 0.8 / resist, total = bare + grass + trees;
  assert.ok(Math.abs(land.wetness(i, aero, 295) - total) < 1e-12);
  before = step();
  const seep = (before.surface - drawn * bare / total) * (1 - Math.exp(-dt / 86400));
  assert.ok(Math.abs(land.deepWater[i] - (before.deep - drawn * trees / total)) < 1e-12, 'the trees\' share of the evaporation comes from the store');
  assert.ok(Math.abs(land.surface[i] - (before.surface - drawn * bare / total - seep)) < 1e-12, 'the bare ground\'s from the surface layer');
  assert.ok(Math.abs(land.soil[i] - (before.soil - drawn * grass / total + seep)) < 1e-12, 'the grass\'s from the root zone');
  set(0.01, 0.1, 0.8, 0.4);
  land.wetness(i, aero, 295);
  land.deepWater[i] = 0.05;
  step();
  assert.deepEqual([land.deepWater[i], land.soil[i]], [0, 0], 'an emptied store gives what it holds and the root zone the rest');
  set(45, 180, 0.8, 0.4); land.snow[i] = 0.1;
  assert.equal(land.wetness(i, aero, 295), 1);
  before = step();
  assert.equal(land.deepWater[i], 180, 'under snow the evaporation leaves the store');
  assert.ok(Math.abs(land.soil[i] - (before.soil - (drawn - 0.1))) < 1e-12);
});

test('the deep store is part of the land\'s water: over CPU model steps its change and the other stores\' close on the rain, the evaporation and the runoff', () => {
  const model = createModel(new Grid(N), { ...options, land: { soil: 'layers', deepCapacity: 300 } });
  const { geography, land } = model;
  initializeState(model, {}).forEach((a, k) => model.state[k].set(a));
  for (let i = 0; i < C; i++) if (geography.land[i]) model.state[6][i] = 0;
  if (model.load) model.load();
  land.initialize();
  const cells = Array.from({ length: C }, (_, i) => i).filter((i) => geography.land[i]), grown = cells.filter((i) => !(geography.iceSheet && geography.iceSheet[i] > 0));
  for (const i of grown) { land.vegetation[i] = 0.8; land.canopy[i] = 0.4; land.soil[i] = i % 2 ? 300 : 60; land.surface[i] = i % 2 ? 15 : 0; land.deepWater[i] = [0, 150, 280, 300][(i >> 1) % 4]; }
  const start = Float64Array.from(land.deepWater), area = (i) => model.mesh.areaCell[i];
  const held = () => cells.reduce((s, i) => s + area(i) * (land.soil[i] + land.surface[i] + land.snow[i] + land.deepWater[i] - land.iceLoss[i]), 0);
  const ran = () => cells.reduce((s, i) => s + area(i) * land.runoff[i], 0);
  const waterBefore = land.water(), heldBefore = held(), ranBefore = ran();
  let rain = 0, evaporated = 0;
  for (let n = 0; n < 4; n++) {
    model.step(1350);
    for (const i of cells) { rain += area(i) * model.moist.rain[i]; evaporated += area(i) * model.radiation.evaporation[i] * 1350; }
  }
  const change = held() - heldBefore, runoff = ran() - ranBefore;
  let drew = 0, filled = 0;
  for (const i of cells) { if (land.deepWater[i] < start[i] - 1e-9) drew++; if (land.deepWater[i] > start[i] + 1e-9) filled++; }
  console.log(`${cells.length} land cells over four steps at N=16: the store fell in ${drew} and rose in ${filled}; rain ${(rain / heldBefore).toExponential(2)}, evaporation ${(evaporated / heldBefore).toExponential(2)}, runoff ${(runoff / heldBefore).toExponential(2)} of the stores`);
  assert.ok(drew > 0.1 * grown.length && filled > 0.1 * grown.length, `the trees drew on the store in ${drew} cells, the root zone recharged it in ${filled}`);
  assert.ok(Math.abs(change - (rain - evaporated - runoff)) < 1e-9 * heldBefore, `the stores change by ${change}, the rain less the evaporation and the runoff ${rain - evaporated - runoff}`);
  assert.ok(Math.abs(land.water() - waterBefore - (rain - evaporated)) < 1e-9 * waterBefore, 'water() counts the store');
});

test('a save keeps land.deep through the land and the state file, a save without it starts the store at the root zone\'s fill, and one above the capacity is capped', async () => {
  const land = fresh(geo16, { start: 'observed' });
  const cells = Array.from({ length: C }, (_, i) => i).filter((i) => geo16.land[i] && !(geo16.iceSheet[i] > 0));
  cells.forEach((i, n) => { land.deepWater[i] = 300 * ((n * 0.37) % 1); land.soil[i] = 300 * ((n * 0.61) % 1); });
  const saved = land.serialize();
  assert.ok(saved.deep instanceof Float64Array && saved.deep.length === C);
  const back = createLandSurface(mesh, geo16, { deepCapacity: 300 });
  back.load(saved);
  assert.deepEqual(Array.from(back.deepWater), Array.from(land.deepWater));
  const decoded = await decodeState(encodeState({ N, time: 0, land: saved }));
  assert.ok(decoded.land.deep && decoded.land.deep.length === C, 'the state file carries land.deep');
  const filed = createLandSurface(mesh, geo16, { deepCapacity: 300 });
  filed.load(decoded.land);
  for (let i = 0; i < C; i++) assert.equal(filed.deepWater[i], Math.fround(land.deepWater[i]), `cell ${i} through the file`);
  const { deep, ...without } = saved, other = createLandSurface(mesh, geo16, { deepCapacity: 200 });
  other.load(without);
  for (let i = 0; i < C; i++) {
    const expected = geo16.land[i] && !(geo16.iceSheet[i] > 0) ? Math.min(1, Math.max(0, saved.soil[i] / 300)) * 200 : 0;
    assert.equal(other.deepWater[i], expected, `cell ${i}: the root zone's fill times the capacity`);
  }
  other.load(saved);
  for (const i of cells) assert.equal(other.deepWater[i], Math.min(200, saved.deep[i]), `cell ${i}: capped at the capacity`);
  const plain = createLandSurface(mesh, geo16, {});
  plain.load(saved);
  assert.ok(!('deep' in plain.serialize()), 'a land without the option drops it');
});

test('with deepCapacity 0, or without it, four steps of the CPU model give the plain land\'s digests, for the land and for the atmosphere with it', () => {
  for (const extra of [{}, { deepCapacity: 0 }]) {
    const model = createModel(new Grid(N), { ...options, land: { soil: 'layers', ...extra } });
    initializeState(model, {}).forEach((a, k) => model.state[k].set(a));
    for (let i = 0; i < C; i++) if (model.geography.land[i]) model.state[6][i] = 0;
    if (model.load) model.load();
    model.land.initialize();
    for (let n = 0; n < 4; n++) model.step(900);
    const saved = model.land.serialize();
    assert.equal(digest(landArrays(saved)), PLAIN_LAND, `${JSON.stringify(extra)}: the land`);
    assert.equal(digest([...model.state, ...landArrays(saved)]), PLAIN_STEPPED, `${JSON.stringify(extra)}: the atmosphere and the land`);
  }
});
