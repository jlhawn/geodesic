import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createOcean } from '../js/ocean/layered.module.js';
import { createGeography, topographyFromInt16 } from '../js/geography.module.js';
import { seawaterDensity, salinityForDensity } from '../js/ocean/seawater.module.js';
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

test('with thermobaric the pressure force over a seamount still leaves a uniform stratification nearly at rest', () => {
  const { ocean, surfaceT, ice, flux, calm } = hybridOcean({ thermobaric: true }, (lat, lon) => 4000 - 3500 * Math.exp(-((lat / (8 * DEG)) ** 2 + (lon / (8 * DEG)) ** 2)));
  for (let n = 0; n < 64; n++) ocean.advance(surfaceT, ice, flux, calm, DT);
  const speed = maxSpeed(ocean);
  console.log(`seamount to 500 m under a uniform stratification, thermobaric: fastest current after a day ${speed.toExponential(2)} m/s`);
  assert.ok(speed < 0.01, `spurious current ${speed} m/s`);
});

test('below 1500 m, deep water made 2 K colder and fresher at the same surface density on one side of the globe drives a current toward the warm side only with thermobaric', () => {
  const run = (thermobaric) => {
    const { ocean, surfaceT, ice, flux, calm } = hybridOcean({ thermobaric });
    const L = ocean.layers, cold = (i) => Math.cos(mesh.lonCell[i]) < 0, { h, Q, W } = ocean;
    for (let i = 0; i < C; i++) {
      if (!cold(i)) continue;
      let top = 0;
      for (let k = 0; k < L; k++) {
        const n = k * C + i, hv = h[n];
        if (hv > 1 && top + 0.5 * hv > 1500) {
          const t = Q[n] / hv, s = W[n] / hv, colder = t - 2;
          Q[n] = hv * colder; W[n] = hv * salinityForDensity(seawaterDensity(t, s), colder);
        }
        top += hv;
      }
    }
    for (let n = 0; n < 4; n++) ocean.advance(surfaceT, ice, flux, calm, DT);
    let along = 0, agree = 0, count = 0;
    for (let e = 0; e < E; e++) {
      const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1];
      if (cold(a) === cold(b)) continue;
      let top = 0;
      for (let k = 0; k < L; k++) {
        const hv = h[k * C + a];
        if (hv > 1 && top + 0.5 * hv > 2000) { const toWarm = (cold(a) ? 1 : -1) * ocean.u[k * E + e]; along += toWarm; agree += toWarm > 0 ? 1 : 0; count++; }
        top += hv;
      }
    }
    return { mean: along / count, share: agree / count, fastest: maxSpeed(ocean) };
  };
  const off = run(false), on = run(true);
  console.log(`deep flow toward the warm side after four steps: thermobaric off ${off.mean.toExponential(2)} m/s (fastest ${off.fastest.toExponential(2)}), on ${on.mean.toExponential(2)} m/s in ${(100 * on.share).toFixed(0)} % of the front's deep edges`);
  assert.ok(off.fastest < 1e-3, `without thermobaric the compensated front moves water at ${off.fastest} m/s`);
  assert.ok(on.mean > 2e-3 && on.share > 0.9, `with thermobaric the deep flow toward the warm side averages ${on.mean} m/s in ${on.share} of the edges`);
});
