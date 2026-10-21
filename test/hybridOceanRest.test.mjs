import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createOcean } from '../js/ocean/layered.module.js';
import { createGeography, topographyFromInt16 } from '../js/geography.module.js';
import { DT, DEG, mesh, C, E, hybridOcean, maxSpeed } from './helpers/layered.mjs';

test('at N=16\'s coupled ocean step of an hour and a half, an eighth of the inertial period at the poles, the sea of the real topography left alone from rest keeps its currents slow for twelve days', () => {
  const geography = createGeography(mesh, topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer));
  const ocean = createOcean(mesh, { vertical: 'hybrid', everySteps: 1, geography });
  const surfaceT = Float64Array.from({ length: C }, (_, i) => 300 - 25 * Math.sin(mesh.latCell[i]) ** 2), ice = new Float64Array(C), flux = new Float64Array(C), calm = new Float64Array(E);
  ocean.initialize(surfaceT, ice);
  let fastest = 0;
  for (let n = 0; n < 12 * 86400 / (4 * DT); n++) {
    ocean.advance(surfaceT, ice, flux, calm, 4 * DT);
    fastest = Math.max(fastest, maxSpeed(ocean));
  }
  assert.ok(fastest < 0.5, `currents reach ${fastest} m/s in twelve days`);
});

test('over a seamount the pressure force at the depths the layers share leaves a uniform stratification nearly at rest', () => {
  const { ocean, surfaceT, ice, flux, calm } = hybridOcean({}, (lat, lon) => 4000 - 3500 * Math.exp(-((lat / (8 * DEG)) ** 2 + (lon / (8 * DEG)) ** 2)));
  for (let n = 0; n < 64; n++) ocean.advance(surfaceT, ice, flux, calm, DT);
  const speed = maxSpeed(ocean);
  console.log(`seamount to 500 m under a uniform stratification: fastest current after a day ${speed.toExponential(2)} m/s`);
  assert.ok(speed < 0.01, `spurious current ${speed} m/s`);
  assert.equal(ocean.diagnostics().oceanDried, 0);
});
