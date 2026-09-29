import { test } from 'node:test';
import assert from 'node:assert/strict';
import { createOcean, rebinOcean, savedDensities, LAYER_DENSITIES, LAYER_SALINITIES, UNLISTED_LAYER_DENSITIES, EPS, THIN } from '../js/ocean/layered.module.js';
import { seawaterDensity, labelTemperature } from '../js/ocean/seawater.module.js';
import { FREEZING_POINT } from '../js/physics/ice.module.js';
import { encodeState, decodeState } from '../js/stateFile.module.js';
import { createGeography, syntheticTopography } from '../js/geography.module.js';
import { RHO_AIR, DRAG, DEG, mesh, C, E, zonalWindOnEdges, UNLISTED_OCEANS } from './helpers/layered.mjs';

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

const labels = (densities, salinities) => {
  const labelS = [35, ...salinities];
  return { labelS, labelT: [1025, ...densities].map((r, k) => Math.max(FREEZING_POINT, labelTemperature(r, labelS[k]))) };
};
const columnSums = ({ h, T, S }, L) => Array.from({ length: C }, (_, i) => {
  let water = 0, heat = 0, salt = 0;
  for (let k = 0; k < L; k++) { const n = k * C + i; water += h[n]; heat += h[n] * T[n]; salt += h[n] * S[n]; }
  return [water, heat, salt];
});
const drift = (a, b) => Math.max(...a.flatMap((column, i) => column.map((x, q) => (x === 0 ? Math.abs(b[i][q]) : Math.abs(b[i][q] - x) / Math.abs(x)))));
function unlistedState(classes) {
  const geography = createGeography(mesh, syntheticTopography(180, 360, (lat, lon) => (Math.abs(lon) < 40 * DEG && lat > 12 * DEG && lat < 48 * DEG ? 500 : -4000 + 1500 * Math.sin(3 * lon) * Math.cos(2 * lat))));
  const ocean = createOcean(mesh, { everySteps: 1, geography, ...classes });
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => Math.max(FREEZING_POINT, 303 - 35 * Math.sin(lat) ** 2)), ice = new Float64Array(C), flux = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  const stressField = Float64Array.from(zonalWindOnEdges(mesh, (lat) => 8 * Math.cos(3 * lat)), (w) => RHO_AIR * DRAG * 8 * w);
  for (let n = 0; n < 20; n++) ocean.advance(surfaceT, ice, flux, stressField, 1350);
  const { densities, ...saved } = ocean.serialize();
  return { geography, surfaceT, ice, saved };
}

test('a saved ocean carries its class list through the state file, and one on this list loads exactly as one that carries none', async () => {
  const ocean = createOcean(mesh, { everySteps: 1 });
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => 275 + 25 * Math.cos(lat) ** 2), ice = new Float64Array(C), flux = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  for (let n = 0; n < 5; n++) ocean.advance(surfaceT, ice, flux, new Float64Array(E).fill(0.05), 1350);
  const decoded = (await decodeState(encodeState({ N: 8, K: 1, day: 0, time: 0, ocean: ocean.serialize() }))).ocean;
  assert.deepEqual(Array.from(decoded.densities), LAYER_DENSITIES);
  assert.deepEqual(savedDensities(decoded, C), LAYER_DENSITIES);
  const { densities, ...unlisted } = decoded;
  const listed = createOcean(mesh, { everySteps: 1 }), plain = createOcean(mesh, { everySteps: 1 });
  listed.load(decoded, Float64Array.from(surfaceT), ice);
  plain.load(unlisted, Float64Array.from(surfaceT), ice);
  for (const name of ['h', 'u', 'Q', 'W', 'eta']) assert.deepEqual(Array.from(listed[name]), Array.from(plain[name]), name);
});

for (const classes of UNLISTED_OCEANS) {
  const count = classes.densities.length + 1;
  test(`rebinOcean carries a ${count}-layer ocean onto the class list and back, each column keeping its water, heat and salt to 1e-12, the mixed layer as it was and the classes near their labels`, () => {
    const { saved } = unlistedState(classes);
    assert.deepEqual(savedDensities(saved, C), classes.densities);
    const onto = rebinOcean(saved, classes.densities, LAYER_DENSITIES, C, labels(LAYER_DENSITIES, LAYER_SALINITIES));
    const L = LAYER_DENSITIES.length + 1, before = columnSums(saved, count);
    assert.equal(onto.h.length, L * C);
    assert.ok(drift(before, columnSums(onto, L)) < 1e-12, `column sums moved by ${drift(before, columnSums(onto, L))}`);
    for (let i = 0; i < C; i++) {
      for (const name of ['h', 'T', 'S']) assert.ok(Math.abs(onto[name][i] - saved[name][i]) <= 1e-12 * Math.abs(saved[name][i]), `the mixed layer's ${name} at ${i}`);
      if (!(saved.h[i] > 0)) continue;
      for (let k = 1; k < L; k++) assert.ok(onto.h[k * C + i] >= EPS - 1e-12, `class ${k} at ${i} holds at least its token`);
    }
    for (let e = 0; e < E; e++) assert.equal(onto.u[e], saved.u[e]);
    const off = [];
    for (let n = C; n < L * C; n++) if (onto.h[n] > THIN) off.push(Math.abs(seawaterDensity(onto.T[n], onto.S[n]) - LAYER_DENSITIES[Math.floor(n / C) - 1]));
    off.sort((a, b) => a - b);
    const median = off[off.length >> 1], ninety = off[Math.floor(0.9 * off.length)];
    const back = rebinOcean(onto, LAYER_DENSITIES, classes.densities, C, labels(classes.densities, classes.salinities));
    assert.ok(drift(before, columnSums(back, count)) < 1e-12);
    let moved = 0, total = 0;
    for (let n = C; n < saved.h.length; n++) { moved += (back.h[n] - saved.h[n]) ** 2; total += saved.h[n] ** 2; }
    console.log(`${count} layers onto the class list: classes off their labels by a median ${median.toFixed(4)}, 90% within ${ninety.toFixed(4)}, at most ${off[off.length - 1].toFixed(4)} kg/m³; there and back, the classes' thicknesses move by ${(100 * Math.sqrt(moved / total)).toFixed(1)}% rms`);
    assert.ok(median < 0.005 && ninety < 0.02, `${median}, ${ninety}`);
    assert.ok(Math.sqrt(moved / total) < 0.1);
  });

  test(`load() takes a ${count}-layer ocean that carries no class list onto the class list with the same column sums`, () => {
    const { geography, surfaceT, ice, saved } = unlistedState(classes);
    const ocean = createOcean(mesh, { everySteps: 1, geography });
    ocean.load(saved, Float64Array.from(surfaceT), ice);
    assert.equal(ocean.layers, LAYER_DENSITIES.length + 1);
    const before = columnSums(saved, count);
    for (let i = 0; i < C; i++) {
      if (!ocean.cellOcean[i]) continue;
      let water = 0, heat = 0, salt = 0;
      for (let k = 0; k < ocean.layers; k++) { const n = k * C + i; water += ocean.h[n]; heat += ocean.Q[n]; salt += ocean.W[n]; }
      assert.ok(Math.abs(water - before[i][0]) < 1e-12 * before[i][0] && Math.abs(heat - before[i][1]) < 1e-12 * before[i][1] && Math.abs(salt - before[i][2]) < 1e-12 * before[i][2], `column ${i}`);
      assert.ok(Math.abs(water - ocean.D[i] - ocean.eta[i]) < 1e-6);
      assert.equal(ocean.h[i], saved.h[i]);
    }
    assert.deepEqual(Array.from(ocean.u.subarray(0, E)), Array.from(saved.u.slice(0, E)));
  });
}
