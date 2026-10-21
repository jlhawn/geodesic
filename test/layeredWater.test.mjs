import { test } from 'node:test';
import assert from 'node:assert/strict';
import { createOcean, LAYER_DENSITIES, LAYER_SALINITIES, runoffOutlets } from '../js/ocean/layered.module.js';
import { seawaterDensity, labelTemperature, thermalExpansion } from '../js/ocean/seawater.module.js';
import { FREEZING_POINT } from '../js/physics/ice.module.js';
import { createGeography, syntheticTopography } from '../js/geography.module.js';
import { DEG, mesh, C, E, totalHeatSalt } from './helpers/layered.mjs';

test('the equation of state recovers each class label at its label salinity above the freezing point, and expands little near freezing', () => {
  LAYER_DENSITIES.forEach((r, k) => {
    const s = LAYER_SALINITIES[k], t = labelTemperature(r, s);
    assert.ok(Math.abs(seawaterDensity(t, s) - r) < 1e-9 && t > FREEZING_POINT, `label ${r} at ${s}: ${t} K`);
  });
  assert.ok(thermalExpansion(273.15, 35) < 7e-5 && thermalExpansion(298.15, 35) > 2.9e-4);
  assert.ok(seawaterDensity(FREEZING_POINT, 34.5) < LAYER_DENSITIES[LAYER_DENSITIES.length - 1], 'polar surface water floats on the deepest class');
});

test('runoff reaches the coastal sea beside the land it ran off, freshening it by exactly that water', () => {
  const geography = createGeography(mesh, syntheticTopography(180, 360, (lat, lon) => (Math.abs(lon) < 40 * DEG && lat > 12 * DEG && lat < 48 * DEG ? -4000 : 500)));
  const ocean = createOcean(mesh, { everySteps: 1, geography });
  const runoff = Float64Array.from(geography.land, (l) => (l ? 2 : 0));
  ocean.fresh.fill(0);
  ocean.accumulate(null, null, 0, runoff);
  let delivered = 0, ranOff = 0;
  for (let i = 0; i < C; i++) {
    ranOff += mesh.areaCell[i] * runoff[i];
    delivered -= mesh.areaCell[i] * ocean.fresh[i];
    if (ocean.fresh[i] === 0) continue;
    assert.ok(ocean.cellOcean[i], `cell ${i} receives runoff but is land`);
    let coastal = false;
    for (let m = 0; m < mesh.nEdgesOnCell[i]; m++) if (geography.land[mesh.cellsOnCell[mesh.maxEdges * i + m]]) coastal = true;
    assert.ok(coastal, `sea cell ${i} receives runoff but touches no land`);
  }
  assert.ok(ranOff > 0 && Math.abs(delivered - ranOff) < 1e-9 * ranOff, `delivered ${delivered} of ${ranOff}`);
});

test('runoff flows down the terrain to the coast the land slopes toward', () => {
  const west = -40 * DEG, east = 40 * DEG;
  const geography = createGeography(mesh, syntheticTopography(180, 360, (lat, lon) => (lon > west && lon < east && Math.abs(lat) < 50 * DEG ? 100 + 2000 * (lon - west) / (east - west) : -4000)));
  const outlet = runoffOutlets(mesh, geography);
  let land = 0, slope = 0, westward = 0;
  for (let i = 0; i < C; i++) {
    if (!geography.land[i]) { assert.equal(outlet[i], i); continue; }
    land++;
    assert.ok(outlet[i] >= 0 && !geography.land[outlet[i]], `land cell ${i} has no sea outlet`);
    if (mesh.lonCell[i] < 0 || mesh.lonCell[i] > east - 10 * DEG || Math.abs(mesh.latCell[i]) > 30 * DEG) continue;
    slope++;
    const lon = Math.atan2(mesh.xCell[3 * outlet[i] + 1], mesh.xCell[3 * outlet[i]]);
    assert.ok(lon < mesh.lonCell[i], `cell ${i} on the eastern slope drains uphill to lon ${(lon / DEG).toFixed(0)}`);
    if (lon < west + 5 * DEG) westward++;
  }
  assert.ok(land > 0.1 * C, `${land} land cells`);
  assert.ok(westward > 0.9 * slope, `${westward} of ${slope} high cells reach the western coast`);
});

test('wind stress reaches the water under sea ice, scaled by the cover and the transmission factor', () => {
  const ocean = createOcean(mesh, { iceStressTransmission: 0.8, everySteps: 1 });
  const C = mesh.nCells, E = mesh.nEdges;
  const surfaceT = new Float64Array(C).fill(FREEZING_POINT + 2), ice = new Float64Array(C), flux = new Float64Array(C);
  const total = new Float64Array(E).fill(0.1);
  ocean.advance(surfaceT, ice, flux, total, 1350);
  const open = Float64Array.from(ocean.stress);
  ice.fill(1); surfaceT.fill(FREEZING_POINT);
  ocean.advance(surfaceT, ice, flux, total, 1350);
  const full = Float64Array.from(ocean.stress);
  ocean.advance(surfaceT, ice, flux, total, 1350, new Float64Array(C).fill(0.5));
  const half = Float64Array.from(ocean.stress);
  let checked = 0;
  for (let e = 0; e < E; e++) if (open[e] !== 0) { checked++; assert.ok(Math.abs(full[e] - 0.8 * open[e]) < 1e-12 && Math.abs(half[e] - 0.9 * open[e]) < 1e-12, `edge ${e}: open ${open[e]} full ${full[e]} half ${half[e]}`); }
  assert.ok(checked > 0);
});
