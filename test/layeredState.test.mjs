import { test } from 'node:test';
import assert from 'node:assert/strict';
import { createOcean, rebinOcean, savedDensities, LAYER_DENSITIES, LAYER_SALINITIES, UNLISTED_LAYER_DENSITIES, THIN_LAYER_CLASSES, REBIN_CARRIED, EPS, THIN } from '../js/ocean/layered.module.js';
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
});

test('a sea cell that arrives empty takes its sea neighbors\' area-weighted column, sea level and SST, never a column at the land\'s temperature, and one ringed by empty cells takes it on the next sweep', () => {
  const geography = createGeography(mesh, syntheticTopography(180, 360, (lat, lon) => (Math.abs(lon) < 60 * DEG && Math.abs(lat) < 50 * DEG ? -4000 : 500)));
  const ocean = createOcean(mesh, { everySteps: 1, geography });
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => 275 + 25 * Math.cos(lat) ** 2), ice = new Float64Array(C), flux = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  const stressField = Float64Array.from(zonalWindOnEdges(mesh, (lat) => 8 * Math.cos(3 * lat)), (w) => RHO_AIR * DRAG * 8 * w);
  for (let n = 0; n < 20; n++) ocean.advance(surfaceT, ice, flux, stressField, 1350);
  const saved = ocean.serialize(), L = saved.h.length / C, { cellsOnCell, nEdgesOnCell, maxEdges, areaCell } = mesh;
  const neighbors = (i) => Array.from({ length: nEdgesOnCell[i] }, (_, m) => cellsOnCell[maxEdges * i + m]);
  const flat = (i) => ocean.cellOcean[i] && neighbors(i).every((j) => ocean.cellOcean[j] && ocean.D[j] === ocean.D[i]);
  const inner = Array.from({ length: C }, (_, i) => i).find((i) => flat(i) && neighbors(i).every(flat));
  assert.ok(inner !== undefined, 'the basin has a cell whose neighbors and theirs are all sea at its depth');
  const ring = neighbors(inner), outer = ring[0], emptied = [inner, ...ring];
  const empty = (i) => { for (let k = 0; k < L; k++) { saved.h[k * C + i] = 0; saved.T[k * C + i] = 0; saved.S[k * C + i] = 0; } saved.eta[i] = 0; };
  const mean = (i, wet) => {
    const from = neighbors(i).filter(wet);
    const area = from.reduce((s, j) => s + areaCell[j], 0), at = (k, j) => k * C + j;
    const h = Array.from({ length: L }, (_, k) => from.reduce((s, j) => s + areaCell[j] * saved.h[at(k, j)], 0) / area);
    const heat = from.reduce((s, j) => s + areaCell[j] * saved.h[j] * saved.T[j], 0) / area;
    return { h, eta: from.reduce((s, j) => s + areaCell[j] * saved.eta[j], 0) / area, sst: heat / h[0] };
  };
  const expected = mean(outer, (j) => !emptied.includes(j));
  const landT = 316, start = Float64Array.from(surfaceT);
  for (const i of emptied) { empty(i); start[i] = landT; }
  const reloaded = createOcean(mesh, { everySteps: 1, geography });
  reloaded.load(saved, start, ice);
  for (const i of emptied) {
    let sum = 0;
    for (let k = 0; k < L; k++) sum += reloaded.h[k * C + i];
    assert.ok(Math.abs(sum - reloaded.D[i] - reloaded.eta[i]) < 1e-6, `column ${i} sums to its depth plus sea level (${sum} vs ${reloaded.D[i]} + ${reloaded.eta[i]})`);
    assert.ok(Math.abs(reloaded.T0[i] - landT) > 5, `cell ${i} is not at the land's temperature (${reloaded.T0[i]})`);
    assert.equal(start[i], reloaded.T0[i], `the surface over cell ${i} takes its mixed layer's temperature`);
  }
  for (let k = 0; k < L; k++) assert.ok(Math.abs(reloaded.h[k * C + outer] - expected.h[k]) < 1e-9, `layer ${k} of the emptied cell is its wet neighbors' mean (${reloaded.h[k * C + outer]} vs ${expected.h[k]})`);
  assert.ok(Math.abs(reloaded.eta[outer] - expected.eta) < 0.05, `its sea level is its neighbors' (${reloaded.eta[outer]} vs ${expected.eta})`);
  assert.ok(Math.abs(reloaded.T0[outer] - expected.sst) < 1e-9, `its SST is its neighbors' (${reloaded.T0[outer]} vs ${expected.sst})`);
  const area = ring.reduce((s, j) => s + areaCell[j], 0), ringMean = (values, k) => ring.reduce((s, j) => s + areaCell[j] * values[k * C + j], 0) / area;
  for (let k = 0; k < L; k++) assert.ok(Math.abs(reloaded.h[k * C + inner] - ringMean(reloaded.h, k)) < 1e-9, `layer ${k} of the cell ringed by empty ones is the mean of the ring as the first sweep filled it`);
  assert.ok(Math.abs(reloaded.eta[inner] - ringMean(reloaded.eta, 0)) < 0.05, `the ringed cell's sea level is the ring's (${reloaded.eta[inner]})`);
});

test('sea cells that arrive empty in an enclosed basin no neighbor sweep reaches take the climatology column and sum to their depth plus sea level', () => {
  const { xCell, cellsOnCell, nEdgesOnCell, maxEdges } = mesh;
  const neighbors = (i) => Array.from({ length: nEdgesOnCell[i] }, (_, m) => cellsOnCell[maxEdges * i + m]);
  const dot = (lat, lon, i) => Math.cos(lat) * Math.cos(lon) * xCell[3 * i] + Math.cos(lat) * Math.sin(lon) * xCell[3 * i + 1] + Math.sin(lat) * xCell[3 * i + 2];
  const within = (lat, lon, i) => neighbors(i).every((j) => dot(lat, lon, i) >= dot(lat, lon, j));
  const center = Array.from({ length: C }, (_, i) => i).reduce((best, i) => (dot(0, 150 * DEG, i) > dot(0, 150 * DEG, best) ? i : best), 0);
  const lake = [center, neighbors(center)[0]];
  const geography = createGeography(mesh, syntheticTopography(180, 360, (lat, lon) => ((Math.abs(lon) < 60 * DEG && Math.abs(lat) < 50 * DEG) || lake.some((i) => within(lat, lon, i)) ? -4000 : 500)));
  const ocean = createOcean(mesh, { everySteps: 1, geography });
  assert.ok(lake.every((i) => ocean.cellOcean[i]) && lake.every((i) => neighbors(i).every((j) => lake.includes(j) || !ocean.cellOcean[j])), 'the lake is two sea cells ringed by land');
  ocean.initialize(new Float64Array(C).fill(290), new Float64Array(C));
  const saved = ocean.serialize(), L = saved.h.length / C;
  for (const i of lake) { for (let k = 0; k < L; k++) { saved.h[k * C + i] = 0; saved.T[k * C + i] = 0; saved.S[k * C + i] = 0; } saved.eta[i] = 0; }
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => 275 + 25 * Math.cos(lat) ** 2), ice = new Float64Array(C);
  const climatology = createOcean(mesh, { everySteps: 1, geography });
  climatology.initialize(Float64Array.from(surfaceT), ice);
  const reloaded = createOcean(mesh, { everySteps: 1, geography });
  reloaded.load(saved, Float64Array.from(surfaceT), ice);
  for (const i of lake) {
    let sum = 0;
    for (let k = 0; k < L; k++) {
      const n = k * C + i;
      for (const name of ['h', 'Q', 'W']) assert.equal(reloaded[name][n], climatology[name][n], `layer ${k}'s ${name} at lake cell ${i} is the climatology's`);
      sum += reloaded.h[n];
    }
    assert.ok(Math.abs(sum - reloaded.D[i] - reloaded.eta[i]) < 1e-6, `column ${i} sums to its depth plus sea level (${sum} vs ${reloaded.D[i]} + ${reloaded.eta[i]})`);
    assert.ok(Math.abs(reloaded.T0[i] - surfaceT[i]) < 1e-9, `the lake's mixed layer is at the surface temperature it was loaded with (${reloaded.T0[i]} vs ${surfaceT[i]})`);
  }
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

function withRestart(saved) {
  const arrays = {};
  REBIN_CARRIED.forEach((name, q) => { arrays[name] = Float64Array.from({ length: (name === 'remainder' ? 3 : 1) * C }, (_, i) => (name === 'capacity' ? 4e6 * (40 + (i % 11)) : name === 'previousT0' ? saved.T[i % C] : 1e-4 * (q + 1) * ((i % 7) - 3))); });
  arrays.previousIceVolume.fill(0);
  return { ...saved, ...arrays };
}

test('rebinOcean carries a stepped 44-class ocean onto the two lighter classes, the classes from 1021 to 1028 kept but for their tokens, each column\'s water, heat and salt to 1e-14, the restart arrays carried and the result restarting on the CPU loader', () => {
  const { geography, surfaceT, ice, saved: plain } = unlistedState({});
  const saved = withRestart(plain), from = savedDensities(saved, C), to = THIN_LAYER_CLASSES.densities;
  assert.deepEqual(from, LAYER_DENSITIES);
  const onto = rebinOcean(saved, from, to, C, labels(to, THIN_LAYER_CLASSES.salinities));
  const Lo = from.length + 1, Ln = to.length + 1, shift = Ln - Lo;
  assert.equal(onto.h.length, Ln * C);
  assert.ok(drift(columnSums(saved, Lo), columnSums(onto, Ln)) < 1e-14, `column sums moved by ${drift(columnSums(saved, Lo), columnSums(onto, Ln))}`);
  let tokens = 0, largest = 0;
  const short = Array.from({ length: C }, (_, i) => { let sum = 0; for (let j = 1; j < Lo; j++) sum += Math.max(0, EPS - saved.h[j * C + i]); return sum; });
  for (let j = 2; j < Lo; j++) for (let i = 0; i < C; i++) {
    const was = j * C + i, now = (j + shift) * C + i, dh = saved.h[was] - onto.h[now];
    largest = Math.max(largest, Math.abs(dh));
    assert.ok(dh >= -Math.max(0, EPS - saved.h[was]) - 1e-12 && dh <= 2 * EPS + short[i] + 1e-12, `class ${from[j - 1]} at ${i} moved by ${dh} m`);
    if (Math.abs(dh) > 1e-12) { tokens++; continue; }
    for (const name of ['T', 'S']) assert.ok(Math.abs(onto[name][now] - saved[name][was]) <= 1e-12 * Math.abs(saved[name][was]), `class ${from[j - 1]}'s ${name} at ${i}`);
  }
  for (let i = 0; i < C; i++) {
    if (!(saved.h[i] > 0)) continue;
    for (const name of ['h', 'T', 'S']) assert.ok(Math.abs(onto[name][i] - saved[name][i]) <= 1e-12 * Math.abs(saved[name][i]), `the mixed layer's ${name} at ${i}`);
    for (let k = 1; k < Ln; k++) assert.ok(onto.h[k * C + i] >= EPS - 1e-12, `class ${k} at ${i} holds at least its token`);
  }
  for (let n = 0; n < Ln * C; n++) if (onto.h[n] > 0) assert.ok(Math.abs(onto.Q[n] - onto.h[n] * onto.T[n]) <= 1e-12 * Math.abs(onto.Q[n]) && Math.abs(onto.W[n] - onto.h[n] * onto.S[n]) <= 1e-12 * Math.abs(onto.W[n]));
  for (let k = 0; k < Ln; k++) assert.deepEqual(Array.from(onto.u.subarray(k * E, (k + 1) * E)), Array.from(saved.u.slice(Math.max(k - shift, Math.min(k, 1)) * E, (Math.max(k - shift, Math.min(k, 1)) + 1) * E)), `layer ${k}'s velocity`);
  for (const name of REBIN_CARRIED) assert.deepEqual(Array.from(onto[name]), Array.from(saved[name]), name);
  console.log(`44 classes onto 46: ${tokens} class cells from 1021 to 1028 gave tokens, at most ${largest.toFixed(4)} m`);
  const ocean = createOcean(mesh, { everySteps: 1, geography, ...THIN_LAYER_CLASSES });
  ocean.load(onto, Float64Array.from(surfaceT), ice);
  assert.equal(ocean.layers, Ln);
  for (const name of ['fresh', 'capacity', 'penetrated']) assert.deepEqual(Array.from(ocean[name]), Array.from(saved[name]), `the loaded ocean's ${name}`);
});

test('a 45-layer ocean that carries no class list loads onto the two lighter classes as the 44-class list rebinned, not rebuilt', () => {
  const { geography, surfaceT, ice, saved } = unlistedState({});
  assert.deepEqual(savedDensities(saved, C), UNLISTED_LAYER_DENSITIES[2]);
  const listed = createOcean(mesh, { everySteps: 1, geography, ...THIN_LAYER_CLASSES }), unlisted = createOcean(mesh, { everySteps: 1, geography, ...THIN_LAYER_CLASSES });
  listed.load({ ...saved, densities: LAYER_DENSITIES }, Float64Array.from(surfaceT), ice);
  unlisted.load(saved, Float64Array.from(surfaceT), ice);
  assert.equal(unlisted.layers, THIN_LAYER_CLASSES.densities.length + 1);
  for (const name of ['h', 'u', 'Q', 'W', 'eta']) assert.deepEqual(Array.from(unlisted[name]), Array.from(listed[name]), name);
  const before = columnSums(saved, LAYER_DENSITIES.length + 1);
  for (let i = 0; i < C; i++) {
    if (!unlisted.cellOcean[i]) continue;
    assert.equal(unlisted.h[i], saved.h[i]);
    let water = 0;
    for (let k = 0; k < unlisted.layers; k++) water += unlisted.h[k * C + i];
    assert.ok(Math.abs(water - before[i][0]) < 1e-12 * before[i][0], `column ${i}`);
  }
});
