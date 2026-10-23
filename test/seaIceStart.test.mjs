import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { readTopography, geographyOptions, FRESH_GEOGRAPHY } from '../js/geography.module.js';
import { sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
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
