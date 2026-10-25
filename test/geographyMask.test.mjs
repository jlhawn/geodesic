import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { gridMesh } from '../js/stateGrid.module.js';
import { checkGeography, geographyHeader, savedGeography, geographyOptions, createGeography, readTopography, subgridUrl, decodeSubgrid, landSubgrid, FRESH_GEOGRAPHY } from '../js/geography.module.js';
import { withCadence } from '../js/cadence.module.js';
import { CLIMATOLOGY_FILE, loadClimatology } from '../js/ocean/climatology.module.js';

const N = 16, DT = 1350 * 16 / N, MASKED = { ...FRESH_GEOGRAPHY, mask: 'n1024' };
const nonFinite = (fields) => Object.entries(fields).filter(([, a]) => ArrayBuffer.isView(a) && a.some((v) => !Number.isFinite(v))).map(([name]) => name);

test('the n1024 mask is a geography of ETOPO 2022 with the observed ice and its shelves as land or cavities, in the state header only when chosen, and the raster records are as they were', () => {
  assert.deepEqual(checkGeography(FRESH_GEOGRAPHY), { topography: 'etopo2022', iceSheets: 'observed', iceShelves: 'cavity', mask: 'raster' });
  const record = checkGeography(MASKED);
  assert.deepEqual(record, { ...FRESH_GEOGRAPHY, mask: 'n1024' });
  assert.deepEqual(geographyHeader(record), { topography: 'etopo2022', iceSheets: 'observed', iceShelves: 'cavity', mask: 'n1024' });
  assert.deepEqual(savedGeography(geographyHeader(record)), record);
  assert.deepEqual(geographyHeader(FRESH_GEOGRAPHY), { topography: 'etopo2022', iceSheets: 'observed', iceShelves: 'cavity' });
  assert.throws(() => checkGeography({ mask: 'n1024' }), /n1024/);
  assert.throws(() => checkGeography({ ...MASKED, iceShelves: 'sea' }), /n1024/);
  assert.throws(() => checkGeography({ ...FRESH_GEOGRAPHY, mask: 'n512' }), /mask/);
});

test(`at N=${N} the n1024 mask file gives the land and cavity geographies, its subgrid file fits it, its lakes and depressions stand at their own heights, and a fresh start on it stays finite for a day on the CPU`, async () => {
  const topography = readTopography('etopo2022'), sheets = geographyOptions(MASKED), mesh = gridMesh(new Grid(N));
  const shelves = createGeography(mesh, topography, sheets.geography), cavity = createGeography(mesh, topography, { ...sheets.geography, iceShelves: 'cavity' });
  assert.equal(shelves.mask, 'n1024');
  assert.ok(landSubgrid(decodeSubgrid(readFileSync(subgridUrl(N, 'etopo2022', 'n1024'))), shelves.land), `data/subgrid_etopo2022_n1024_N${N}.bin does not fit the mask`);
  for (let i = 0; i < mesh.nCells; i++) {
    if (cavity.cavity[i]) assert.deepEqual([shelves.land[i], cavity.land[i], cavity.ocean.land[i]], [1, 1, 0], `cavity cell ${i}`);
    else assert.equal(cavity.ocean.land[i], shelves.land[i], `cell ${i} outside the cavities`);
  }
  for (let e = 0; e < mesh.nEdges; e++) if (cavity.ocean.edgeOcean[e]) assert.ok(cavity.ocean.edgeSill[e] > 0);
  const ocean = withCadence({ ocean: {} }, DT, N).ocean, radiation = withCadence({ radiation: { clearSkyPass: true } }, DT, N).radiation;
  const climatology = await loadClimatology(new URL(`../${CLIMATOLOGY_FILE}`, import.meta.url).pathname);
  const model = createModel(new Grid(N), { topography, geography: sheets.geography, ocean: { ...ocean, climatology }, radiation, land: { ...sheets.land } });
  assert.deepEqual(model.geographyRecord, checkGeography(MASKED));
  const { state, geography, surfaceGeopotential: phis } = model, C = model.mesh.nCells, g = model.core.diagnostics.g;
  for (let i = 0; i < C; i++) {
    if (geography.lake[i]) assert.ok(Math.abs(phis[i] - g * geography.lakeLevel[i]) < 1e-9 * g * 1e4, `lake cell ${i}`);
    else if (!geography.land[i]) assert.equal(phis[i], 0);
  }
  initializeState(model, { geostrophic: !model.surfaceGeopotential }).forEach((values, a) => state[a].set(values));
  for (let i = 0; i < C; i++) if (geography.land[i]) state[6][i] = 0;
  if (model.load) model.load();
  model.ocean.initialize(state[3], state[6]);
  model.land.initialize();
  for (let n = 0; n < 86400 / DT; n++) model.step(DT);
  const atmosphere = nonFinite(Object.fromEntries(['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'].map((name, a) => [name, state[a]])));
  assert.deepEqual([...atmosphere, ...nonFinite({ h: model.ocean.h, u: model.ocean.u, eta: model.ocean.eta })], []);
});
