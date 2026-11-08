import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { readTopography, geographyOptions, FRESH_GEOGRAPHY } from '../js/geography.module.js';
import { sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { lakeIceAlbedo } from '../js/physics/lake.module.js';
import { cpuFrameFields } from '../js/cpuFrame.module.js';
import { encodeSeaIceStart, decodeSeaIceStart, loadSeaIceStart, SEA_ICE_START_FILE, EDGE_CONCENTRATION, NO_SEA } from '../js/physics/seaIceStart.module.js';

const d = Math.PI / 180, N = 16;
const coverGrid = { nLon: 8, nLat: 5, lon0: 0, dLon: 45, lat0: 80, dLat: -10 }, iceGrid = { nLon: 4, nLat: 3, lon0: 0, dLon: 90, lat0: 90, dLat: -45 };
const rows = [100, 50, 60, null, 0];
const synthetic = () => encodeSeaIceStart({
  coverGrid, iceGrid, source: 'synthetic',
  cover: Uint8Array.from({ length: 40 }, (_, n) => (rows[Math.floor(n / 8)] === null ? NO_SEA : rows[Math.floor(n / 8)])),
  thickness: Uint8Array.from({ length: 12 }, (_, n) => [40, 20, 0][Math.floor(n / 4)]),
  snow: Uint8Array.from({ length: 12 }, (_, n) => [50, 25, 0][Math.floor(n / 4)]),
});
const observedStart = (n) => {
  const start = decodeSeaIceStart(readFileSync(new URL(`../${SEA_ICE_START_FILE}`, import.meta.url)));
  const topography = readTopography(FRESH_GEOGRAPHY.topography), sheets = geographyOptions(FRESH_GEOGRAPHY);
  const model = createModel(new Grid(n), { topography, geography: sheets.geography, ocean: false, land: { ...sheets.land }, levels: sigmaInterfaces('bl48') });
  initializeState(model, { geostrophic: !model.surfaceGeopotential, seaIce: start }).forEach((values, a) => model.state[a].set(values));
  model.land.initialize();
  return model;
};

test('a sea-ice start decodes as written and samples bilinearly over sea, nearest within reach, nothing below the edge', () => {
  const start = decodeSeaIceStart(synthetic());
  assert.equal(start.source, 'synthetic');
  assert.deepEqual(start.coverGrid, coverGrid);
  assert.deepEqual(start.iceGrid, iceGrid);
  const between = start.at(75 * d, 0);
  assert.ok(Math.abs(between.concentration - 0.75) < 1e-12);
  assert.ok(Math.abs(between.thickness - 0.05 * (40 - 20 / 3)) < 1e-6);
  assert.ok(Math.abs(between.snow - 2 * (50 - 25 / 3)) < 1e-5);
  assert.ok(Math.abs(start.at(55 * d, 0).concentration - 0.6) < 1e-12, 'a land corner takes no weight');
  assert.equal(start.at(50 * d, 22.5 * d).concentration, 0, 'no sea within the default reach');
  assert.ok(Math.abs(decodeSeaIceStart(synthetic(), { reach: 17 }).at(50 * d, 22.5 * d).concentration - 0.6) < 1e-12, 'the nearest sea within reach');
  const thin = decodeSeaIceStart(encodeSeaIceStart({ coverGrid, iceGrid, cover: new Uint8Array(40).fill(100 * EDGE_CONCENTRATION - 1), thickness: new Uint8Array(12).fill(20), snow: new Uint8Array(12) }));
  assert.deepEqual(thin.at(70 * d, 0), { concentration: 0, thickness: 0, snow: 0 });
});

test(`${SEA_ICE_START_FILE} holds March ice where it lies at the equinox`, async () => {
  const start = await loadSeaIceStart(new URL(`../${SEA_ICE_START_FILE}`, import.meta.url).pathname);
  assert.match(start.source, /ERA5/);
  assert.match(start.source, /GIOMAS/);
  const pole = start.at(89 * d, 0), weddell = start.at(-70 * d, -50 * d);
  assert.ok(pole.concentration > 0.9 && pole.thickness > 1.5 && pole.thickness < 4 && pole.snow > 20 && pole.snow < 200, JSON.stringify(pole));
  assert.ok(weddell.concentration > 0.8 && weddell.thickness > 0.5 && weddell.thickness < 3, JSON.stringify(weddell));
  for (const [lat, lon] of [[0, 0], [40, -30], [-45, 100]]) assert.equal(start.at(lat * d, lon * d).concentration, 0, `${lat}, ${lon}`);
});

test(`a fresh start at N=${N} takes the observed March ice, its concentration and its snow, and the land leaves the sea's snow alone`, async () => {
  const start = decodeSeaIceStart(readFileSync(new URL(`../${SEA_ICE_START_FILE}`, import.meta.url)));
  const topography = readTopography(FRESH_GEOGRAPHY.topography), sheets = geographyOptions(FRESH_GEOGRAPHY);
  const model = createModel(new Grid(N), { topography, geography: sheets.geography, ocean: false, land: { ...sheets.land }, levels: sigmaInterfaces('bl48') });
  const { state, mesh } = model, C = mesh.nCells, land = model.geography.land;
  initializeState(model, { geostrophic: !model.surfaceGeopotential, seaIce: start }).forEach((values, a) => state[a].set(values));
  model.land.initialize();
  const area = [0, 0];
  let snowy = 0, iced = 0, partial = 0;
  for (let i = 0; i < C; i++) {
    const h = state[6][i], c = model.seaIce.concentration[i], observed = start.at(mesh.latCell[i], mesh.lonCell[i]);
    if (land[i]) { assert.equal(h, 0); continue; }
    assert.equal(h, observed.thickness);
    assert.equal(model.land.snow[i], observed.snow);
    if (!(h > 0)) { assert.equal(c, 0); continue; }
    assert.ok(c >= EDGE_CONCENTRATION && c <= 1 && Math.abs(c - observed.concentration) < 1e-12);
    assert.ok(Math.abs(mesh.latCell[i]) > 40 * d);
    assert.equal(state[3][i], 271.35);
    iced++;
    if (model.land.snow[i] > 0) snowy++;
    if (c < 0.95) partial++;
    area[mesh.latCell[i] > 0 ? 0 : 1] += mesh.areaCell[i] * c / 1e12;
  }
  assert.ok(area[0] > 8 && area[0] < 20, `Arctic ice area ${area[0]} Mkm²`);
  assert.ok(area[1] > 1 && area[1] < 8, `Antarctic ice area ${area[1]} Mkm²`);
  assert.ok(snowy > 0.8 * iced && partial > 0, `${snowy} of ${iced} iced cells carry snow, ${partial} partly covered`);
});

test('a fresh start at N=32 on the n1024 mask gives each lake cell the sea\'s start but its snow: the observed March ice, its cover and the sea\'s surface temperature, which FLake replaces at its first step, and no snow, which nothing would', async () => {
  const start = decodeSeaIceStart(readFileSync(new URL(`../${SEA_ICE_START_FILE}`, import.meta.url)));
  const topography = readTopography(FRESH_GEOGRAPHY.topography), sheets = geographyOptions(FRESH_GEOGRAPHY);
  const model = createModel(new Grid(32), { topography, geography: sheets.geography, ocean: false, land: { ...sheets.land }, levels: sigmaInterfaces('bl48') });
  const { state, mesh } = model, lakes = model.lake.cells, options = { geostrophic: !model.surfaceGeopotential, seaIce: start };
  const sea = initializeState({ ...model, lake: null }, options), seaCover = Float64Array.from(model.seaIce.concentration), seaSnow = Float64Array.from(model.land.snow);
  initializeState(model, options).forEach((values, a) => state[a].set(values));
  model.land.initialize();
  assert.ok(lakes.some((i) => seaSnow[i] > 0 && sea[6][i] > 0), 'the sea\'s start lays the observed ice and its snow on a lake cell');
  for (const i of lakes) assert.deepEqual([state[6][i], model.seaIce.concentration[i], state[3][i], model.land.snow[i]], [sea[6][i], seaCover[i], sea[3][i], 0], `lake cell ${i}: ice, cover, surface temperature and snow`);
});

test('at N=32 on the n1024 mask after the fresh start, the CPU\'s albedo map takes at every cell the diffuse albedo its radiation uses, on land, open sea, snowy sea ice, an open lake and a frozen one, where the ice\'s albedo is the lake ice\'s and not the sea ice\'s, and its surface-albedo diagnostic is the map\'s mean', () => {
  const model = observedStart(32), { state, mesh, seaIce } = model, land = model.geography.land, lakeCell = model.lake.isLake, h = state[6];
  const lake = (i) => lakeIceAlbedo(state[3][i]);
  const sea = (i) => seaIce.albedo(h[i], null, seaIce.snow[i], seaIce.cover(i, h[i]), state[3][i], seaIce.snowAlbedo[i]);
  const radiated = (i) => (land[i] ? model.land.albedo(i) : lakeCell[i] && h[i] > 0 ? lake(i) : sea(i));
  const cells = [...mesh.areaCell.keys()], frozen = model.lake.cells.filter((i) => h[i] > 0);
  console.log(`frozen lake cells ${frozen.join(', ')}: the lake ice's albedo ${frozen.map((i) => lake(i).toFixed(3)).join(', ')}, the sea ice's ${frozen.map((i) => sea(i).toFixed(3)).join(', ')}`);
  assert.ok(frozen.some((i) => Math.abs(lake(i) - sea(i)) > 0.05), "the start freezes a lake where the two ices' albedos differ");
  const kinds = {
    land: (i) => land[i], 'open sea': (i) => !land[i] && !lakeCell[i] && !(h[i] > 0),
    'snowy sea ice': (i) => !land[i] && !lakeCell[i] && h[i] > 0 && seaIce.snow[i] > 0, 'an open lake': (i) => lakeCell[i] && !(h[i] > 0),
  };
  for (const [kind, holds] of Object.entries(kinds)) assert.ok(cells.some(holds), `the start holds ${kind}`);
  const { albedo } = cpuFrameFields(model, { fields: ['albedo'] });
  for (const i of cells) assert.equal(albedo[i], Math.fround(radiated(i)), `cell ${i}`);
  let sum = 0, area = 0;
  for (const i of cells) { sum += mesh.areaCell[i] * albedo[i]; area += mesh.areaCell[i]; }
  const { surfaceAlbedo } = model.diagnostics();
  assert.ok(Math.abs(surfaceAlbedo - sum / area) < 1e-7, `the diagnostic ${surfaceAlbedo} against the map's mean ${sum / area}`);
});

test('at N=32 on the n1024 mask after the fresh start, the CPU\'s sea-ice diagnostics count the ice on the sea and not a frozen lake\'s', () => {
  const model = observedStart(32), { state, mesh, seaIce } = model, land = model.geography.land, lakeCell = model.lake.isLake;
  let area = 0, sea = 0, volume = 0, lakes = 0;
  for (let i = 0; i < mesh.nCells; i++) {
    const a = mesh.areaCell[i], h = state[6][i], c = seaIce.concentration[i];
    area += a;
    if (!(h > 0)) continue;
    if (lakeCell[i]) lakes += a * c;
    else if (!land[i]) { sea += a * c; volume += a * c * h; }
  }
  console.log(`sea ice ${(100 * sea / area).toFixed(3)} % of the area, ${(volume / sea).toFixed(3)} m thick; the frozen lakes' ice ${(100 * lakes / area).toFixed(4)} %`);
  assert.ok(lakes / sea > 5e-4, `the start freezes a lake: its ice ${lakes / area} of the area against the sea's ${sea / area}`);
  const d = model.diagnostics();
  assert.ok(Math.abs(d.iceFraction - sea / area) < 1e-12, `iceFraction ${d.iceFraction} against the sea's ${sea / area}`);
  assert.ok(Math.abs(d.iceThickness - volume / sea) < 1e-12, `iceThickness ${d.iceThickness} against the sea's ${volume / sea}`);
});
