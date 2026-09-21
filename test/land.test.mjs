import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createGeography, syntheticTopography } from '../js/geography.module.js';
import { createLandSurface } from '../js/physics/land.module.js';
import { MELTING_POINT } from '../js/physics/ice.module.js';

const mesh = createModel(new Grid(8), { physics: false }).mesh;

test('a hemisphere continent covers half the area, its coast rings the meridian, and edges are classified consistently', () => {
  const geography = createGeography(mesh, syntheticTopography(180, 360, (lat, lon) => (Math.cos(lon) > 0 ? 500 : -4000)));
  assert.ok(Math.abs(geography.landArea - 0.5) < 0.03, `land area ${geography.landArea.toFixed(3)}`);
  for (let i = 0; i < mesh.nCells; i++) {
    const inland = Math.abs(Math.cos(mesh.lonCell[i])) > 0.2 && Math.abs(mesh.latCell[i]) < 1.4;
    if (inland) assert.equal(geography.land[i], Math.cos(mesh.lonCell[i]) > 0 ? 1 : 0, `cell ${i} at lon ${mesh.lonCell[i].toFixed(2)}`);
    if (inland && Math.abs(Math.cos(mesh.lonCell[i])) > 0.4) assert.ok(Math.abs(geography.elevation[i] - (Math.cos(mesh.lonCell[i]) > 0 ? 500 : -4000)) < 1);
  }
  for (let e = 0; e < mesh.nEdges; e++) {
    const a = geography.land[mesh.cellsOnEdge[2 * e]], b = geography.land[mesh.cellsOnEdge[2 * e + 1]];
    assert.equal(geography.edgeOcean[e], !a && !b ? 1 : 0);
  }
  for (const e of geography.coastEdges) {
    const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1];
    if (Math.abs(mesh.latCell[a]) > 1.4 || Math.abs(mesh.latCell[b]) > 1.4) continue;
    const lon = 0.5 * (mesh.lonCell[a] + mesh.lonCell[b]);
    assert.ok(Math.abs(Math.cos(lon)) < 0.25, `coast edge at lon ${lon.toFixed(2)}`);
  }
  assert.ok(geography.coastEdges.length > 0.8 * mesh.nCells / 8 && geography.coastEdges.length < 1.5 * mesh.nCells / 8);
});

test('an all-ocean raster leaves no land and no coast', () => {
  const geography = createGeography(mesh, syntheticTopography(90, 180, () => -3000));
  assert.equal(geography.landArea, 0);
  assert.equal(geography.coastEdges.length, 0);
  assert.ok(geography.edgeOcean.every((v) => v === 1));
});

test('the bucket conserves water: rain, evaporation, melt and runoff balance the stores', () => {
  const geography = createGeography(mesh, syntheticTopography(90, 180, () => 100));
  const land = createLandSurface(mesh, geography, { bucketCapacity: 150 });
  land.initialize();
  const i = 0, area = mesh.areaCell[i];
  const surfaceT = new Float64Array(mesh.nCells).fill(290), flux = new Float64Array(mesh.nCells);
  const before = land.water();
  land.deposit(i, 40, 290);
  land.update(i, surfaceT, flux, 1e-4, 3600);
  land.deposit(i, 100, 290);
  const rain = 140, evaporated = 1e-4 * 3600;
  assert.ok(Math.abs(land.water() - before - area * (rain - evaporated)) < 1e-6 * area);
  assert.equal(land.soil[i], 150);
  assert.ok(Math.abs(land.runoff[i] - (75 + rain - evaporated - 150)) < 1e-9);
  assert.ok(Math.abs(land.wetness(i) - 1) < 1e-12);
  land.soil[i] = 30;
  assert.ok(Math.abs(land.wetness(i) - 30 / 112.5) < 1e-12);
});

test('snow accumulates below freezing, raises the albedo, holds the surface at the melting point while it melts, and drains into the bucket', () => {
  const geography = createGeography(mesh, syntheticTopography(90, 180, () => 100));
  const land = createLandSurface(mesh, geography, { heatCapacity: 1e6, albedo: 0.25, snowAlbedo: 0.7, fullSnow: 20 });
  land.initialize();
  const i = 3;
  land.deposit(i, 10, MELTING_POINT - 5);
  assert.equal(land.snow[i], 10);
  assert.ok(Math.abs(land.albedo(i) - (0.25 + 0.5 * 0.45)) < 1e-12);
  const surfaceT = new Float64Array(mesh.nCells).fill(MELTING_POINT - 1), flux = new Float64Array(mesh.nCells).fill(400);
  const soilBefore = land.soil[i];
  land.update(i, surfaceT, flux, 0, 3600);
  const energy = 400 * 3600 - 1 * 1e6;
  const melted = energy / 3.34e5;
  assert.ok(melted < 10);
  assert.ok(Math.abs(land.snow[i] - (10 - melted)) < 1e-9);
  assert.ok(Math.abs(land.soil[i] - (soilBefore + melted)) < 1e-9);
  assert.ok(Math.abs(surfaceT[i] - MELTING_POINT) < 1e-9);
  flux.fill(2000);
  land.update(i, surfaceT, flux, 0, 3600);
  assert.equal(land.snow[i], 0);
  assert.ok(surfaceT[i] > MELTING_POINT);
  assert.ok(Math.abs(land.albedo(i) - 0.25) < 1e-12);
});
