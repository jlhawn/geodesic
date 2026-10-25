import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { createLandSurface } from '../js/physics/land.module.js';
import { MELTING_POINT } from '../js/physics/ice.module.js';
import { saturationHumidity } from '../js/physics/moist.module.js';
import { R_DRY, CP_DRY, P0 } from '../js/dynamics/sigmaCore.module.js';
import { runoffOutlets } from '../js/ocean/layered.module.js';
import { gridMesh, stateGrid } from '../js/stateGrid.module.js';
import { DEFAULT_RELAX } from '../js/grid.module.js';
import { syntheticTopography, createGeography, encodeLandIce, decodeLandIce, landIceUrl, geographyOptions, readTopography, savedGeography, geographyHeader, pickGeography, checkGeography, stateGeography, FRESH_GEOGRAPHY, LEGACY_GEOGRAPHY, ICE_DENSITY, SEAWATER_DENSITY } from '../js/geography.module.js';
import { PREVIOUS_LEVELS } from './helpers/previousDefaults.mjs';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const topography = readTopography('etopo2022'), etopo1 = readTopography('etopo1');
const landIce = decodeLandIce(readFileSync(landIceUrl));
const observed = geographyOptions({ ...FRESH_GEOGRAPHY, iceShelves: 'land' }, landIce);
const deg = 180 / Math.PI;
const sparse = ({ rows, cols, points, count, thickness, floating = null, floatThickness = null, floatSurface = null, floatBed = null }) => {
  const index = Uint32Array.from(Array.from(count.keys()).filter((p) => count[p] > 0)), at = (field) => (p) => (field ? field[p] : 0);
  return decodeLandIce(encodeLandIce({
    rows, cols, points, index, count: Uint8Array.from(index, at(count)), thickness: Uint16Array.from(index, at(thickness)), floating: Uint8Array.from(index, at(floating)),
    floatThickness: Uint16Array.from(index, at(floatThickness)), floatSurface: Uint16Array.from(index, at(floatSurface)), floatBed: Int16Array.from(index, at(floatBed)),
  }));
};

test('the land-ice file carries the pixels holding ice: their ice points and those afloat, the thickness of each, and the surface of the floating ice and the bed under it', () => {
  const field = (f) => Array.from({ length: 7 }, (_, k) => f(k));
  const given = {
    rows: 3, cols: 5, points: 225, firnAir: 15, index: Uint32Array.from(field((k) => 2 * k + 1)), count: Uint8Array.from(field((k) => 30 * k + 15)), floating: Uint8Array.from(field((k) => 10 * k)),
    thickness: Uint16Array.from(field((k) => 4099 * k + 7)), floatThickness: Uint16Array.from(field((k) => 977 * k)), floatSurface: Uint16Array.from(field((k) => 11 * k + 3)), floatBed: Int16Array.from(field((k) => -900 * k + 50)),
  };
  const back = decodeLandIce(encodeLandIce(given));
  assert.deepEqual([back.rows, back.cols, back.points, back.firnAir], [3, 5, 225, 15]);
  for (const name of ['index', 'count', 'floating', 'thickness', 'floatThickness', 'floatSurface', 'floatBed']) assert.deepEqual(Array.from(back[name]), Array.from(given[name]), name);
  assert.throws(() => decodeLandIce(new Uint8Array(40)), /not a land-ice file/);
});

test('with land ice a land cell is ice sheet where ice covers more than half its land, with the mean thickness of that ice', () => {
  const mesh = gridMesh(new Grid(16)), rows = 90, cols = 180;
  const continent = syntheticTopography(rows, cols, (lat, lon) => (Math.cos(lon) > 0 ? 500 : -4000));
  const iced = (share, depth) => sparse({ rows, cols, points: 225, count: Uint8Array.from(continent.data, (h) => (h > 0 ? share : 0)), thickness: Uint16Array.from(continent.data, (h, n) => (h > 0 ? depth(Math.floor(n / cols)) : 0)) });
  const options = { landBridges: {}, seaStraits: {} }, rule = createGeography(mesh, continent, options);
  assert.equal(rule.iceFraction, undefined);
  const full = createGeography(mesh, continent, { ...options, landIce: iced(225, (r) => 100 + 10 * r) });
  const most = createGeography(mesh, continent, { ...options, landIce: iced(120, () => 2000) }), least = createGeography(mesh, continent, { ...options, landIce: iced(110, () => 2000) });
  let coastal = 0;
  for (let i = 0; i < mesh.nCells; i++) {
    assert.equal(full.land[i], rule.land[i]);
    assert.equal(full.shelfFraction[i], 0);
    if (!full.land[i]) { assert.equal(full.iceSheet[i], 0); assert.equal(full.iceThickness[i], 0); continue; }
    assert.equal(full.iceSheet[i], 1);
    assert.ok(Math.abs(full.iceFraction[i] - full.landFraction[i]) < 1e-12, `cell ${i}: ice ${full.iceFraction[i]}, land ${full.landFraction[i]}`);
    const row = (90 - mesh.latCell[i] * deg) / 2;
    assert.ok(Math.abs(full.iceThickness[i] - (100 + 10 * row)) < 10 * 12, `cell ${i} at row ${row.toFixed(1)}: ${full.iceThickness[i]} m`);
    assert.equal(most.iceSheet[i], 1);
    assert.equal(most.iceThickness[i], 2000);
    assert.equal(least.iceSheet[i], 0);
    assert.equal(least.iceThickness[i], 0);
    if (full.landFraction[i] < 0.9) coastal++;
  }
  assert.ok(coastal > 20, `${coastal} coastal land cells`);
  assert.throws(() => createGeography(mesh, continent, { landIce: { ...iced(225, () => 1), rows: 45 } }), /land ice is 45/);
});

test('with the ice shelves at sea the floating ice is taken away: its points count at the bed beneath it, and it is no ice sheet', () => {
  const mesh = gridMesh(new Grid(16)), rows = 90, cols = 180, freeboard = 50, H = Math.round(freeboard / (1 - ICE_DENSITY / SEAWATER_DENSITY)), bed = -800;
  const lonOf = (n) => -Math.PI + ((n % cols) + 0.5) * 2 * Math.PI / cols, grounded = (lon) => lon < 0 && lon > -Math.PI / 2, afloat = (lon) => lon >= 0 && lon < Math.PI / 2;
  const raster = syntheticTopography(rows, cols, (lat, lon) => (grounded(lon) ? 500 : afloat(lon) ? freeboard : bed));
  const where = (f, yes) => (h, n) => (f(lonOf(n)) ? yes : 0);
  const ice = sparse({
    rows, cols, points: 225, count: Uint8Array.from(raster.data, (h, n) => (grounded(lonOf(n)) || afloat(lonOf(n)) ? 225 : 0)), thickness: Uint16Array.from(raster.data, (h, n) => (grounded(lonOf(n)) ? 1500 : afloat(lonOf(n)) ? H : 0)),
    floating: Uint8Array.from(raster.data, where(afloat, 225)), floatThickness: Uint16Array.from(raster.data, where(afloat, H)), floatSurface: Uint16Array.from(raster.data, where(afloat, freeboard)),
    floatBed: Int16Array.from(raster.data, where(afloat, bed)),
  });
  const options = { landBridges: {}, seaStraits: {}, landIce: ice }, onLand = createGeography(mesh, raster, options), atSea = createGeography(mesh, raster, { ...options, iceShelves: 'sea' });
  let shelves = 0, sheets = 0;
  for (let i = 0; i < mesh.nCells; i++) {
    if (onLand.shelfFraction[i] === 1) {
      shelves++;
      assert.deepEqual([onLand.land[i], onLand.iceSheet[i], atSea.land[i], atSea.iceSheet[i]], [1, 1, 0, 0], `shelf cell ${i}`);
      assert.ok(Math.abs(onLand.elevation[i] - freeboard) < 1e-9 && Math.abs(atSea.elevation[i] - bed) < 0.1, `shelf cell ${i}: ${onLand.elevation[i]} m on land, ${atSea.elevation[i]} m at sea`);
    } else if (onLand.shelfFraction[i] === 0 && onLand.iceSheet[i]) {
      sheets++;
      assert.deepEqual([atSea.land[i], atSea.iceSheet[i], atSea.iceThickness[i], atSea.elevation[i]], [1, 1, onLand.iceThickness[i], onLand.elevation[i]], `grounded cell ${i}`);
    }
  }
  assert.ok(shelves > mesh.nCells / 32 && sheets > mesh.nCells / 32, `${shelves} shelf cells, ${sheets} grounded`);
  assert.throws(() => createGeography(mesh, raster, { iceShelves: 'sea' }), /need the land ice/);
  assert.throws(() => createGeography(mesh, raster, { ...options, iceShelves: 'lid' }), /iceShelves is land or sea or cavity, not lid/);
});

test('the repository\'s land ice puts the Antarctic and Greenland ice sheets on the ETOPO 2022 mesh at their observed area and volume', () => {
  assert.deepEqual([landIce.rows, landIce.cols, landIce.points], [topography.rows, topography.cols, 225]);
  const mesh = gridMesh(new Grid(32)), geography = createGeography(mesh, topography, observed.geography), plain = createGeography(mesh, topography);
  const sums = { south: [0, 0], north: [0, 0] };
  for (let i = 0; i < mesh.nCells; i++) {
    if (!geography.iceSheet[i]) continue;
    const lat = mesh.latCell[i] * deg, lon = mesh.lonCell[i] * deg, a = mesh.areaCell[i];
    assert.ok(lat < -50 || (lat > 58 && lon > -75 && lon < -10), `ice sheet at ${lat.toFixed(1)}, ${lon.toFixed(1)}`);
    assert.equal(geography.land[i], 1);
    assert.ok(geography.iceThickness[i] > 0 && geography.iceThickness[i] < 5000, `${geography.iceThickness[i]} m`);
    const s = sums[lat < 0 ? 'south' : 'north']; s[0] += a; s[1] += a * geography.iceThickness[i];
  }
  const [aArea, aVolume] = sums.south.map((x, k) => x / (k ? 1e15 : 1e12)), [gArea, gVolume] = sums.north.map((x, k) => x / (k ? 1e15 : 1e12));
  console.log(`N=32: Antarctic ice ${aArea.toFixed(2)} million km², ${aVolume.toFixed(1)} million km³; Greenland ${gArea.toFixed(2)} million km², ${gVolume.toFixed(2)} million km³`);
  assert.ok(aArea > 13.0 && aArea < 14.2 && aVolume > 25.5 && aVolume < 28.5, `Antarctica ${aArea} million km², ${aVolume} million km³`);
  assert.ok(gArea > 1.5 && gArea < 2.0 && gVolume > 2.6 && gVolume < 3.3, `Greenland ${gArea} million km², ${gVolume} million km³`);
  assert.deepEqual(Array.from(geography.land), Array.from(plain.land));
});

test('a state records its geography where it differs from ETOPO1 and the rule, fresh starts take ETOPO 2022 with the observed ice, and a raster other than the state\'s is refused', () => {
  assert.deepEqual(savedGeography({}), LEGACY_GEOGRAPHY);
  assert.deepEqual(savedGeography({ iceSheets: 'observed' }), { ...LEGACY_GEOGRAPHY, iceSheets: 'observed' });
  assert.deepEqual(geographyHeader(LEGACY_GEOGRAPHY), {});
  assert.deepEqual(geographyHeader(FRESH_GEOGRAPHY), { topography: 'etopo2022', iceSheets: 'observed', iceShelves: 'cavity' });
  assert.deepEqual(pickGeography({ N: 16, topography: 'etopo2022', iceSheets: 'observed', iceShelves: 'land' }), { topography: 'etopo2022', iceSheets: 'observed' });
  assert.throws(() => checkGeography({ topography: 'etopo5' }), /topography is etopo1 or etopo2022/);
  assert.throws(() => checkGeography({ iceSheets: 'bedrock' }), /iceSheets is rule or observed/);
  assert.throws(() => checkGeography({ iceSheets: 'rule', iceShelves: 'sea' }), /need the observed ice sheets/);
  assert.deepEqual(geographyOptions(LEGACY_GEOGRAPHY), { geography: {}, land: {} });
  const read = geographyOptions({ iceSheets: 'observed', iceShelves: 'sea' });
  assert.equal(read.geography.landIce.index.length, landIce.index.length);
  assert.equal(read.geography.iceShelves, 'sea');
  assert.deepEqual(read.land, { iceSheetBalance: true });
  const fresh = stateGeography(null), old = stateGeography({ N: 16 });
  assert.deepEqual([fresh.topography.name, old.topography.name], ['etopo2022', 'etopo1']);
  assert.deepEqual([old.geography, old.land], [{}, {}]);
  const recorded = { N: 16, grid: { relax: DEFAULT_RELAX }, topography: 'etopo2022', iceSheets: 'observed' };
  assert.throws(() => stateGrid(recorded, etopo1), /on the etopo2022 topography, not etopo1/);
  assert.equal(stateGrid(recorded, topography).relax, DEFAULT_RELAX);
  assert.equal(stateGrid(recorded).relax, DEFAULT_RELAX);
});

test('with the ice-sheet balance an ice sheet starts full at the cap, sheds the snow above it as discharge and melts as bare ice at the melting point, every drop counted', () => {
  const mesh = gridMesh(new Grid(16)), geography = createGeography(mesh, topography, observed.geography), C = mesh.nCells, a = mesh.areaCell;
  const land = createLandSurface(mesh, geography, observed.land), plain = createLandSurface(mesh, geography);
  land.initialize(); plain.initialize();
  let sheet = -1, ground = -1;
  for (let i = 0; i < C; i++) {
    if (!geography.land[i]) continue;
    assert.equal(land.snow[i], geography.iceSheet[i] ? 1000 : 0);
    assert.equal(plain.snow[i], 0);
    if (geography.iceSheet[i] && sheet < 0) sheet = i;
    if (!geography.iceSheet[i] && ground < 0) ground = i;
  }
  assert.ok(sheet >= 0 && ground >= 0);
  assert.equal(land.albedo(sheet), 0.8);
  const before = land.water();
  land.deposit(sheet, 30, 250, 900); land.deposit(ground, 30, 250, 900); plain.deposit(sheet, 30, 250, 900);
  assert.deepEqual([land.snow[sheet], land.runoff[sheet], land.discharge[sheet]], [1000, 30, 30]);
  assert.deepEqual([land.snow[ground], land.runoff[ground], land.discharge[ground]], [30, 0, 0]);
  assert.deepEqual([plain.snow[sheet], plain.runoff[sheet], plain.discharge[sheet]], [30, 0, 0]);
  assert.ok(Math.abs(land.water() - before - 30 * (a[sheet] + a[ground])) < 1e-12 * before, 'the discharge stays in the land\'s water until the ocean takes it');

  land.snow[sheet] = 10;
  assert.ok(Math.abs(land.albedo(sheet) - 0.65) < 1e-12, `${land.albedo(sheet)} half covered`);
  land.snow[sheet] = 0; plain.snow[sheet] = 0;
  assert.equal(land.albedo(sheet), 0.5);
  assert.equal(plain.albedo(sheet), 0.8);
  const skin = new Float64Array(C).fill(MELTING_POINT + 2), plainSkin = Float64Array.from(skin), flux = new Float64Array(C);
  const held = land.water(), ranOff = land.runoff[sheet], melt = 2 * 1e6 / 3.34e5;
  land.update(sheet, skin, flux, 0, 900); plain.update(sheet, plainSkin, flux, 0, 900);
  assert.equal(skin[sheet], MELTING_POINT);
  assert.ok(Math.abs(land.runoff[sheet] - ranOff - melt) < 1e-9 && Math.abs(land.iceLoss[sheet] - melt) < 1e-9, `melt ${land.runoff[sheet] - ranOff}, ice lost ${land.iceLoss[sheet]}, expected ${melt}`);
  assert.ok(Math.abs(land.water() - held) < 1e-12 * held, 'the melted ice comes from outside the land\'s water');
  assert.deepEqual([plainSkin[sheet], plain.iceLoss[sheet]], [MELTING_POINT + 2, 0]);

  const saved = land.serialize();
  saved.snow[sheet] = 5000; saved.snow[ground] = 5000;
  land.load(saved);
  assert.deepEqual([land.snow[sheet], land.snow[ground]], [1000, 5000]);
  assert.deepEqual([land.runoff[sheet], land.discharge[sheet], land.iceLoss[sheet]], [0, 0, 0]);
});

/*
 * Air at 1.5 times saturation in the two lowest layers over the ice
 * sheets, so that snow falls on them in the first steps.
 */
function snowOnTheSheets(model) {
  const { state, core, geography, mesh } = model, { K } = core.diagnostics, C = mesh.nCells, kappa = R_DRY / CP_DRY;
  core.diagnose(state[0], state[1]);
  const exner = core.diagnostics.exnerLayer;
  const levels = core.levels, top = PREVIOUS_LEVELS[PREVIOUS_LEVELS.length - 3];
  for (let k = 0; k < K; k++) if (levels[k] + levels[k + 1] > 2 * top) for (let i = 0; i < C; i++) if (geography.iceSheet[i]) {
    const n = k * C + i;
    state[4][n] = 1.5 * saturationHumidity(state[1][n] * exner[n], P0 * exner[n] ** (1 / kappa));
  }
}

async function prepare(model, { moisten = false, bareIce = [] } = {}) {
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i]) model.state[6][i] = 0;
  if (moisten) snowOnTheSheets(model);
  for (const i of bareIce) model.state[3][i] = MELTING_POINT + 3;
  if (model.load) model.load();
  model.ocean.initialize(model.state[3], model.state[6]);
  model.land.initialize();
  if (bareIce.length) {
    const saved = await model.land.serialize();
    for (const i of bareIce) saved.snow[i] = 0;
    model.land.load(saved);
  }
  return model;
}

test('in the CPU model the snow above a full ice sheet reaches the sea: the ocean takes every kilogram the land sheds, and the land\'s water balances its snowfall and evaporation', async () => {
  const model = await prepare(createModel(new Grid(16), { topography, ...observed }), { moisten: true });
  const { mesh, geography, land, ocean } = model, C = mesh.nCells, a = mesh.areaCell, dt = 900;
  const accumulate = ocean.accumulate;
  let delivered = 0, calved = 0;
  ocean.accumulate = (evaporation, rain, step, runoff, discharge) => {
    for (let j = 0; j < C; j++) if (geography.land[j]) { delivered += a[j] * runoff[j]; calved += a[j] * discharge[j]; }
    return accumulate(evaporation, rain, step, runoff, discharge);
  };
  const before = land.water();
  let fallen = 0, evaporated = 0;
  for (let n = 0; n < 12; n++) {
    model.step(dt);
    for (let i = 0; i < C; i++) if (geography.land[i]) { fallen += a[i] * model.moist.rain[i]; evaporated += a[i] * model.radiation.evaporation[i] * dt; }
  }
  model.phases.ocean(dt);
  let shed = 0, ran = 0;
  for (let i = 0; i < C; i++) if (geography.land[i]) { shed += a[i] * land.discharge[i]; ran += a[i] * land.runoff[i]; if (geography.iceSheet[i]) assert.ok(land.snow[i] <= 1000, `cell ${i}: ${land.snow[i]} kg/m²`); }
  console.log(`N=16, twelve steps: the sheets shed ${(shed / 1e12).toFixed(3)} Gt of the ${(fallen / 1e12).toFixed(3)} Gt that fell on land; ${(ran / 1e12).toFixed(3)} Gt ran off`);
  assert.ok(shed > 0, 'snow fell on the full sheets');
  assert.ok(Math.abs(calved - shed) <= 1e-12 * shed && Math.abs(delivered - ran) <= 1e-12 * ran, `the ocean took ${calved} of ${shed} kg discharged and ${delivered} of ${ran} kg run off`);
  const imbalance = land.water() - before - (fallen - evaporated);
  assert.ok(Math.abs(imbalance) < 1e-9 * before, `the land's water is off by ${(imbalance / before).toExponential(2)} of the store`);
});

test('both engines take the observed ice sheets and their balance alike at N=16: land codes, full sheets, the snow they shed, the bare ice they melt and the sea their discharge reaches', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const options = { topography, ...observed }, geography = createGeography(gridMesh(new Grid(16)), topography, observed.geography);
  const bareIce = [];
  for (let i = 0, n = 0; i < geography.land.length; i++) if (geography.iceSheet[i] && n++ % 3 === 0) bareIce.push(i);
  const cpu = await prepare(createModel(new Grid(16), options), { moisten: true, bareIce }), gpu = await prepare(await createGpuModel(new Grid(16), options), { moisten: true, bareIce });
  const { land, iceSheet } = cpu.geography, C = cpu.mesh.nCells, a = cpu.mesh.areaCell;
  const start = await gpu.gpu.downloadPhysics();
  let sheets = 0;
  for (let i = 0; i < C; i++) {
    assert.equal(start.LAND[i], land[i] ? (iceSheet[i] ? 2 : 1) : 0, `cell ${i}`);
    if (iceSheet[i]) { sheets++; assert.equal(start.SNOW[i], bareIce.includes(i) ? 0 : 1000, `cell ${i}`); }
  }
  assert.ok(sheets > C / 64 && bareIce.length > sheets / 4, `${sheets} ice-sheet cells, ${bareIce.length} bare`);
  const injected = Float32Array.from({ length: C }, (_, i) => (iceSheet[i] ? 20 : 0));
  for (let i = 0; i < C; i++) { cpu.land.discharge[i] += injected[i]; cpu.land.runoff[i] += injected[i]; }
  for (const name of ['DISCHARGE', 'RUNOFF']) gpu.gpu.device.queue.writeBuffer(gpu.gpu.buffers.PH, 4 * gpu.gpu.layout.PH[name], injected);
  for (let n = 0; n < 8; n++) { cpu.step(900); await gpu.step(900); }
  await gpu.sync();
  const end = await gpu.gpu.downloadPhysics();
  const worst = (pick) => { let max = 0, at = -1; for (let i = 0; i < C; i++) { const d = pick(i); if (d !== null && Math.abs(d) > max) { max = Math.abs(d); at = i; } } return [max, at]; };
  const [skin, skinAt] = worst((i) => (iceSheet[i] ? cpu.state[3][i] - gpu.state[3][i] : null));
  const [snow, snowAt] = worst((i) => (iceSheet[i] ? cpu.land.snow[i] - end.SNOW[i] : null));
  const [shed, shedAt] = worst((i) => (iceSheet[i] ? cpu.land.discharge[i] - end.DISCHARGE[i] : null));
  const [lost, lostAt] = worst((i) => (iceSheet[i] ? cpu.land.iceLoss[i] - end.ICELOSS[i] : null));
  let melted = 0, fell = 0;
  for (let i = 0; i < C; i++) if (iceSheet[i]) { melted += a[i] * cpu.land.iceLoss[i]; fell += a[i] * (cpu.land.discharge[i] - injected[i]); }
  const outlet = runoffOutlets(cpu.mesh, cpu.geography), outlets = [...new Set(Array.from(outlet).filter((o, i) => iceSheet[i] && o >= 0))];
  const [sea, seaAt] = worst((i) => (outlets.includes(i) && !(cpu.state[6][i] > 0) ? cpu.state[3][i] - gpu.state[3][i] : null));
  const [ice, iceAt] = worst((i) => (outlets.includes(i) ? cpu.state[6][i] - gpu.state[6][i] : null));
  let grown = 0;
  for (const o of outlets) grown += cpu.state[6][o] / outlets.length;
  console.log(`N=16, ${sheets} ice-sheet cells (${bareIce.length} bare), eight steps: skin ${skin.toExponential(1)} K, snow ${snow.toExponential(1)} kg/m², discharge ${shed.toExponential(1)} kg/m² (${(fell / 1e12).toFixed(3)} Gt from snowfall), ice melted ${lost.toExponential(1)} kg/m² (${(melted / 1e12).toFixed(3)} Gt); ${outlets.length} outlets: sea ${sea.toExponential(1)} K, sea ice ${ice.toExponential(1)} m (mean ${grown.toFixed(3)} m)`);
  assert.ok(skin < 0.02, `skin ${skin} K at ${skinAt}`);
  assert.ok(snow < 0.05, `snow ${snow} kg/m² at ${snowAt}`);
  assert.ok(shed < 0.05, `discharge ${shed} kg/m² at ${shedAt}`);
  assert.ok(melted > 0 && lost < 0.01, `ice melted ${lost} kg/m² apart at ${lostAt}`);
  assert.ok(fell > 0, 'snow fell on the full sheets');
  assert.ok(sea < 0.01 && ice < 0.002, `outlet sea ${sea} K at ${seaAt}, sea ice ${ice} m at ${iceAt}`);
});
