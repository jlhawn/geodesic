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
  const land = createLandSurface(mesh, geography, { bucketCapacity: 150, vegetation: false });
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
  const land = createLandSurface(mesh, geography, { heatCapacity: 1e6, albedo: 0.25, snowAlbedo: 0.7, fullSnow: 20, vegetation: false });
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

test('loading a land state leaves sea cells dry and bare, except the snow on a cell that has ice to hold it', () => {
  const geography = createGeography(mesh, syntheticTopography(180, 360, (lat, lon) => (Math.cos(lon) > 0 ? 500 : -4000)));
  const land = createLandSurface(mesh, geography, { bucketCapacity: 150, vegetation: false });
  const ice = Float64Array.from({ length: mesh.nCells }, (_, i) => (i % 2 ? 1 : 0));
  land.load({ soil: new Float64Array(mesh.nCells).fill(120), snow: new Float64Array(mesh.nCells).fill(40) }, ice);
  for (let i = 0; i < mesh.nCells; i++) {
    if (geography.land[i]) { assert.equal(land.soil[i], 120); assert.equal(land.snow[i], 40); }
    else { assert.equal(land.soil[i], 0); assert.equal(land.snow[i], ice[i] > 0 ? 40 : 0); }
  }
  land.load({ soil: new Float64Array(mesh.nCells).fill(120), snow: new Float64Array(mesh.nCells).fill(40) });
  for (let i = 0; i < mesh.nCells; i++) if (!geography.land[i]) assert.equal(land.snow[i], 0);
});

const DAY = 86400;
const flat = () => createGeography(mesh, syntheticTopography(90, 180, () => 100));

test('vegetation grows over a wet bucket and dies back over a dry one on its time scales, taking the albedo and the bucket with it', () => {
  const land = createLandSurface(mesh, flat(), { growthTime: 100 * DAY, declineTime: 50 * DAY });
  land.initialize();
  const i = 0, surfaceT = new Float64Array(mesh.nCells).fill(295), flux = new Float64Array(mesh.nCells);
  assert.equal(land.vegetation[i], 1);
  assert.equal(land.capacity(i), 300);
  assert.ok(Math.abs(land.albedo(i) - 0.13) < 1e-12);
  land.vegetation[i] = 0.5; land.soil[i] = 0;
  land.update(i, surfaceT, flux, 0, 50 * DAY);
  assert.ok(Math.abs(land.vegetation[i] - 0.5 * Math.exp(-1)) < 1e-12, `dry: ${land.vegetation[i]}`);
  assert.ok(Math.abs(land.albedo(i) - (0.30 + (0.13 - 0.30) * land.vegetation[i])) < 1e-12);
  assert.equal(land.capacity(i), 300, 'the bucket does not follow the cover');
  const v0 = land.vegetation[i];
  land.soil[i] = land.capacity(i);
  land.update(i, surfaceT, flux, 0, 100 * DAY);
  assert.ok(Math.abs(land.vegetation[i] - (1 - (1 - v0) * Math.exp(-1))) < 1e-12, `wet: ${land.vegetation[i]}`);
  land.vegetation[i] = 0; land.soil[i] = 0.35 * land.capacity(i);
  land.update(i, surfaceT, flux, 0, 1e9 * DAY);
  assert.ok(Math.abs(land.vegetation[i] - 0.5) < 1e-9, `a bucket 35% full settles at half cover: ${land.vegetation[i]}`);
});

test('a browning cell keeps its bucket, and water above the root zone runs off', () => {
  const land = createLandSurface(mesh, flat(), { snowDeclineTime: 10 * DAY });
  land.initialize();
  const i = 5, surfaceT = new Float64Array(mesh.nCells).fill(MELTING_POINT - 10), flux = new Float64Array(mesh.nCells);
  land.deposit(i, 5, MELTING_POINT - 10);
  assert.equal(land.soil[i], 300);
  const before = land.water(), runoff = land.runoff[i];
  land.update(i, surfaceT, flux, 0, 50 * DAY);
  assert.ok(Math.abs(land.vegetation[i] - Math.exp(-5)) < 1e-12);
  assert.ok(Math.abs(land.soil[i] - land.capacity(i)) < 1e-9, `soil ${land.soil[i]} against capacity ${land.capacity(i)}`);
  assert.ok(Math.abs(land.runoff[i] - runoff) < 1e-9, 'browning sheds no water');
  assert.ok(Math.abs(land.water() - before) < 1e-9 * before);
  land.soil[i] = 400; land.vegetation[i] = 0.2; land.snow[i] = 0;
  land.update(i, new Float64Array(mesh.nCells).fill(MELTING_POINT + 5), flux, 0, 1);
  assert.ok(Math.abs(land.soil[i] - land.capacity(i)) < 1e-9 && Math.abs(land.runoff[i] - runoff - 100) < 1e-9, 'water above the root zone runs off');
});

test('under snow the vegetation fades over snowDeclineTime and the snow sets the albedo', () => {
  const land = createLandSurface(mesh, flat(), { snowDeclineTime: 200 * DAY, snowAlbedo: 0.55, fullSnow: 20 });
  land.initialize();
  const i = 2, surfaceT = new Float64Array(mesh.nCells).fill(MELTING_POINT - 10), flux = new Float64Array(mesh.nCells);
  land.deposit(i, 50, MELTING_POINT - 10);
  land.update(i, surfaceT, flux, 0, 200 * DAY);
  assert.ok(Math.abs(land.vegetation[i] - Math.exp(-1)) < 1e-12);
  assert.ok(Math.abs(land.albedo(i) - 0.55) < 1e-12);
});

test('a saved land state without vegetation loads green with full buckets where snow-free and bare under snow', () => {
  const geography = createGeography(mesh, syntheticTopography(180, 360, (lat, lon) => (Math.cos(lon) > 0 ? 500 : -4000)));
  const land = createLandSurface(mesh, geography);
  const snow = Float64Array.from({ length: mesh.nCells }, (_, i) => (mesh.latCell[i] > 1 ? 30 : 0));
  land.load({ soil: new Float64Array(mesh.nCells).fill(20), snow });
  for (let i = 0; i < mesh.nCells; i++) {
    if (!geography.land[i]) { assert.equal(land.vegetation[i], 0); assert.equal(land.soil[i], 0); continue; }
    if (snow[i] > 0 || geography.iceSheet[i]) { assert.equal(land.vegetation[i], 0); assert.equal(land.soil[i], 20); }
    else { assert.equal(land.vegetation[i], 1); assert.equal(land.soil[i], 300); }
  }
  land.load({ soil: new Float64Array(mesh.nCells).fill(20), snow, vegetation: new Float64Array(mesh.nCells).fill(1.5) });
  for (let i = 0; i < mesh.nCells; i++) assert.equal(land.vegetation[i], geography.land[i] && !geography.iceSheet[i] ? 1 : 0);
  assert.deepEqual(Object.keys(land.serialize()), ['soil', 'snow', 'vegetation', 'surface']);
});

test('an ice sheet keeps its albedo under anything and grows nothing', () => {
  const geography = createGeography(mesh, syntheticTopography(90, 180, (lat, lon) => (lat < -1.1 ? 2000 : Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 500 : -4000)));
  const sheet = [...geography.iceSheet].map((s, i) => (s ? i : -1)).filter((i) => i >= 0);
  assert.ok(sheet.length > 0 && sheet.every((i) => geography.land[i] && mesh.latCell[i] < -60 / 57.29578), 'the southern land is an ice sheet');
  assert.ok([...geography.land].some((l, i) => l && !geography.iceSheet[i]), 'the tropical continent is not');
  const land = createLandSurface(mesh, geography);
  land.initialize();
  const i = sheet[0], surfaceT = new Float64Array(mesh.nCells).fill(280), flux = new Float64Array(mesh.nCells);
  assert.equal(land.vegetation[i], 0);
  assert.equal(land.albedo(i), 0.8);
  land.deposit(i, 40, MELTING_POINT - 5);
  assert.equal(land.albedo(i), 0.8, 'snow does not darken it toward the snow albedo');
  land.snow[i] = 0; land.soil[i] = land.capacity(i);
  land.update(i, surfaceT, flux, 0, 400 * DAY);
  assert.equal(land.vegetation[i], 0, 'a wet bucket grows nothing on it');
  land.load({ soil: new Float64Array(mesh.nCells).fill(100), snow: new Float64Array(mesh.nCells) });
  assert.equal(land.vegetation[i], 0, 'a saved green state loads bare on it');
});

test('with vegetation the soil has a surface layer that bare ground evaporates and a root zone the cover transpires through stomata', () => {
  const land = createLandSurface(mesh, flat(), { surfaceCapacity: 15, percolationTime: 86400, stomatalResistance: 70 });
  land.initialize();
  const i = 0, surfaceT = new Float64Array(mesh.nCells).fill(295), flux = new Float64Array(mesh.nCells);
  land.vegetation[i] = 0; land.soil[i] = 20; land.surface[i] = 0;
  const before = land.water();
  land.deposit(i, 10, 290);
  assert.equal(land.surface[i], 10); assert.equal(land.soil[i], 20);
  land.deposit(i, 30, 290);
  assert.equal(land.surface[i], 15);
  const shed = 25 * Math.pow(20 / 300, 4);
  assert.ok(Math.abs(land.soil[i] - (20 + 25 - shed)) < 1e-12 && Math.abs(land.runoff[i] - shed) < 1e-12, 'the overflow infiltrates, a (soil/capacity)⁴ share running off');
  assert.ok(Math.abs(land.water() - before - mesh.areaCell[i] * 40) < 1e-9 * mesh.areaCell[i]);
  assert.ok(Math.abs(land.wetness(i, 0.01, 295) - 1) < 1e-12, 'a wet surface layer evaporates freely from bare ground');
  land.update(i, surfaceT, flux, 1e-4, 3600);
  assert.ok(Math.abs(land.surface[i] - (15 - 0.36) * Math.exp(-3600 / 86400)) < 1e-9, 'bare-ground evaporation and seepage come out of the surface layer');
  land.surface[i] = 0; land.vegetation[i] = 0;
  assert.equal(land.wetness(i, 0.01, 295), 0, 'dry bare ground evaporates nothing however wet the roots');
  land.vegetation[i] = 1; land.soil[i] = land.capacity(i);
  assert.ok(Math.abs(land.wetness(i, 0.01, 295) - 1 / (1 + 70 * 0.01)) < 1e-12, 'a full canopy transpires at the stomatal fraction');
  assert.ok(land.wetness(i, 0.01, 276) < 0.1, 'and closes its stomata in the cold');
  land.vegetation[i] = 0.3; land.soil[i] = land.capacity(i);
  const cold = new Float64Array(mesh.nCells).fill(276), warm = new Float64Array(mesh.nCells).fill(295);
  land.update(i, cold, flux, 0, 30 * DAY);
  assert.ok(Math.abs(land.vegetation[i] - 0.3) < 1e-12, 'no growth in the cold');
  land.update(i, warm, flux, 0, 30 * DAY);
  assert.ok(land.vegetation[i] > 0.35, 'growth when warm');
  const grown = land.vegetation[i];
  land.soil[i] = 0;
  land.update(i, cold, flux, 0, 30 * DAY);
  assert.ok(land.vegetation[i] < grown - 0.01, 'decline needs no warmth');
});
