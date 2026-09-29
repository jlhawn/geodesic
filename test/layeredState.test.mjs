import { test } from 'node:test';
import assert from 'node:assert/strict';
import { createOcean } from '../js/ocean/layered.module.js';
import { createGeography, syntheticTopography } from '../js/geography.module.js';
import { RHO_AIR, DRAG, DEG, mesh, C, E, zonalWindOnEdges } from './helpers/layered.mjs';

test('load() of a saved ocean without layers starts from the climatology', () => {
  const geography = createGeography(mesh, syntheticTopography(180, 360, (lat, lon) => (Math.abs(lon) < 40 * DEG && lat > 12 * DEG && lat < 48 * DEG ? -4000 : 500)));
  const ocean = createOcean(mesh, { everySteps: 1, geography });
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => 275 + 25 * Math.cos(lat) ** 2), ice = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  const climatology = Float64Array.from(ocean.h);
  ocean.u.fill(0.3);
  ocean.load({ h1: new Float64Array(C).fill(55), u1: new Float64Array(E).fill(0.02) }, surfaceT, ice);
  for (let n = 0; n < climatology.length; n++) assert.equal(ocean.h[n], climatology[n]);
  for (let n = 0; n < ocean.u.length; n++) assert.equal(ocean.u[n], 0);
});

test('load(serialize()) reproduces h, u and eta exactly', () => {
  const ocean = createOcean(mesh, { everySteps: 1 });
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => 275 + 25 * Math.cos(lat) ** 2), ice = new Float64Array(C), flux = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  const wind = zonalWindOnEdges(mesh, (lat) => 8 * Math.cos(3 * lat));
  const stressField = Float64Array.from(wind, (w) => RHO_AIR * DRAG * 8 * w);
  for (let n = 0; n < 20; n++) ocean.advance(surfaceT, ice, flux, stressField, 1350);

  const saved = ocean.serialize();
  const reloaded = createOcean(mesh, { everySteps: 1 });
  reloaded.load(saved, Float64Array.from(surfaceT), Float64Array.from(ice));

  assert.deepEqual(Array.from(reloaded.h), Array.from(ocean.h));
  assert.deepEqual(Array.from(reloaded.u), Array.from(ocean.u));
  assert.deepEqual(Array.from(reloaded.eta), Array.from(ocean.eta));
  assert.equal(ocean.diagnostics().oceanLimited, 0);
});

test('load() fits carried-over columns to the bathymetry and fills sea cells that arrive empty', () => {
  const geography = createGeography(mesh, syntheticTopography(180, 360, (lat, lon) => (Math.abs(lon) < 40 * DEG && lat > 12 * DEG && lat < 48 * DEG ? -4000 : 500)));
  const ocean = createOcean(mesh, { geography });
  const C = mesh.nCells, surfaceT = new Float64Array(C).fill(290), ice = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  const saved = ocean.serialize();
  const L = saved.h.length / C, seaCells = [];
  for (let i = 0; i < C; i++) if (ocean.cellOcean[i]) seaCells.push(i);
  const tooDeep = seaCells[0], tooShallow = seaCells[1], empty = seaCells[2];
  for (let k = 0; k < L; k++) saved.h[k * C + tooDeep] *= 3;
  for (let k = 0; k < L; k++) saved.h[k * C + tooShallow] *= 0.2;
  for (let k = 0; k < L; k++) { saved.h[k * C + empty] = 0; saved.T[k * C + empty] = 0; saved.S[k * C + empty] = 0; }
  saved.eta[tooDeep] = 40;
  ocean.load(saved, surfaceT, ice);
  for (const i of seaCells) {
    let sum = 0;
    for (let k = 0; k < L; k++) { assert.ok(ocean.h[k * C + i] > 0, `layer ${k} at cell ${i} has water or a token`); sum += ocean.h[k * C + i]; }
    assert.ok(Math.abs(sum - ocean.D[i] - ocean.eta[i]) < 1e-6, `column ${i} sums to its depth plus sea level (${sum} vs ${ocean.D[i]} + ${ocean.eta[i]})`);
    assert.ok(Math.abs(ocean.eta[i]) <= 5, `sea level at cell ${i} is within 5 m (${ocean.eta[i]})`);
    assert.ok(ocean.h[i] >= 50 - 1e-9 || ocean.h[i] >= ocean.D[i] - 1, `mixed layer at cell ${i} keeps its floor (${ocean.h[i]})`);
    assert.ok(ocean.T0[i] > 250 && ocean.T0[i] < 320, `SST at cell ${i} is physical (${ocean.T0[i]})`);
  }
  assert.ok(Math.abs(ocean.T0[empty] - 290) < 1e-6, 'an empty sea cell takes the climatology surface temperature');
});
