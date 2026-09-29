import { test } from 'node:test';
import assert from 'node:assert/strict';
import { createOcean, EPS, EDDY_SLACK, LAYER_DENSITIES } from '../js/ocean/layered.module.js';
import { seawaterDensity } from '../js/ocean/seawater.module.js';
import { mesh, C, buriedBump, interfaceRange, classContents } from './helpers/layered.mjs';

test('the eddy transport flattens a buried interface bump and keeps every class\'s volume, heat and salt; with eddyDiffusivity 0 it does nothing', () => {
  const k = LAYER_DENSITIES.indexOf(1026.86) + 1, outcropped = LAYER_DENSITIES.filter((r) => r < seawaterDensity(290, 35)).length;
  const { ocean, reshape } = buriedBump({ eddyDiffusivity: 1e6 });
  reshape(k, 200); reshape(k + 1, -200);
  const before = classContents(ocean), h0 = Float64Array.from(ocean.h), start = interfaceRange(ocean, k);
  const ranges = [start];
  for (let n = 0; n < 4; n++) { for (let s = 0; s < 25; s++) ocean.eddyTransport(2700); ranges.push(interfaceRange(ocean, k)); }
  console.log(`the ${ocean.densities[k]}/${ocean.densities[k + 1]} interface's bump over 100 steps: ${ranges.map((r) => r.toFixed(1)).join(' → ')} m`);
  assert.ok(ranges.every((r, n) => n === 0 || r < ranges[n - 1]) && ranges[4] < 0.85 * start, `range ${ranges.join(', ')}`);
  const after = classContents(ocean);
  for (let j = 0; j < ocean.layers; j++) for (const q of ['volume', 'heat', 'salt']) {
    const drift = Math.abs(after[j][q] - before[j][q]) / Math.max(Math.abs(before[j][q]), 1e-300);
    assert.ok(drift < 1e-12, `class ${j} ${q} ${before[j][q]} -> ${after[j][q]}`);
  }
  for (let n = 0; n < ocean.h.length; n++) assert.ok(ocean.h[n] >= Math.min(h0[n], EPS) - 1e-12, `layer ${Math.floor(n / C)} of cell ${n % C} thinned to ${ocean.h[n]}`);
  for (let n = 0; n < (outcropped + 1) * C; n++) assert.equal(ocean.h[n], h0[n], 'the mixed layer and the empty classes under it are untouched');

  const off = buriedBump({ eddyDiffusivity: 0 });
  off.reshape(k, 200); off.reshape(k + 1, -200);
  const untouched = [Float64Array.from(off.ocean.h), Float64Array.from(off.ocean.Q), Float64Array.from(off.ocean.W)];
  for (let s = 0; s < 25; s++) off.ocean.eddyTransport(2700);
  assert.deepEqual([off.ocean.h, off.ocean.Q, off.ocean.W].map((a) => Array.from(a)), untouched.map((a) => Array.from(a)));
});

test('a class with only its token thickness carries no eddy flux, and the interfaces around it flatten together', () => {
  const k = LAYER_DENSITIES.indexOf(1026.89) + 1;
  const { ocean, reshape } = buriedBump({ eddyDiffusivity: 1e6 });
  for (let i = 0; i < C; i++) {
    const n = k * C + i, below = (k + 1) * C + i, moved = ocean.h[n] - EPS, t = ocean.Q[below] / ocean.h[below], s = ocean.W[below] / ocean.h[below];
    ocean.h[below] += moved; ocean.Q[below] = ocean.h[below] * t; ocean.W[below] = ocean.h[below] * s;
    ocean.Q[n] *= EPS / ocean.h[n]; ocean.W[n] *= EPS / ocean.h[n]; ocean.h[n] = EPS;
  }
  reshape(k - 1, 200); reshape(k + 1, -200);
  const token = Float64Array.from(ocean.h.subarray(k * C, (k + 1) * C)), above = interfaceRange(ocean, k - 1), below = interfaceRange(ocean, k);
  for (let s = 0; s < 100; s++) ocean.eddyTransport(2700);
  assert.deepEqual(Array.from(ocean.h.subarray(k * C, (k + 1) * C)), Array.from(token));
  const aboveAfter = interfaceRange(ocean, k - 1), belowAfter = interfaceRange(ocean, k);
  console.log(`around the empty ${ocean.densities[k]} class the interfaces' bumps went ${above.toFixed(1)} → ${aboveAfter.toFixed(1)} and ${below.toFixed(1)} → ${belowAfter.toFixed(1)} m`);
  assert.ok(aboveAfter < 0.85 * above && Math.abs(aboveAfter - belowAfter) < 0.02, `${above} -> ${aboveAfter}, ${below} -> ${belowAfter}`);
});

test('the eddy transport leaves level interfaces over a bumpy bottom level, and at the stability limit stays within each class\'s water', () => {
  const bathymetry = Float64Array.from({ length: C }, (_, i) => 2500 + 2000 * Math.sin(3 * mesh.lonCell[i]) * Math.cos(2 * mesh.latCell[i]));
  const flat = createOcean(mesh, { bathymetry, thermoclineTilt: 0, salinityProfile: () => 35, eddyDiffusivity: 1e12 });
  flat.initialize(new Float64Array(C).fill(290), new Float64Array(C));
  const resting = Float64Array.from(flat.h);
  for (let s = 0; s < 50; s++) flat.eddyTransport(2700);
  let moved = 0, offset = 0;
  for (let n = 0; n < resting.length; n++) moved = Math.max(moved, Math.abs(flat.h[n] - resting[n]));
  for (let k = 1; k < flat.layers - 1; k++) {
    let low = Infinity, high = -Infinity;
    for (let i = 0; i < C; i++) {
      let z = 0, below = 0;
      for (let j = 0; j < flat.layers; j++) if (j <= k) z += resting[j * C + i]; else below += resting[j * C + i];
      if (below > 100) { low = Math.min(low, z); high = Math.max(high, z); }
    }
    if (high > low) offset = Math.max(offset, high - low);
  }
  console.log(`over the bumpy bottom the interfaces start within ${offset.toFixed(3)} m of level where they lie 100 m above the bottom, and the eddy transport at the stability limit moves them by ${moved.toFixed(3)} m in 50 steps`);
  assert.ok(moved < 0.5, `level interfaces over a bumpy bottom moved by ${moved} m`);

  const ocean = createOcean(mesh, { bathymetry, eddyDiffusivity: 1e12 });
  ocean.initialize(Float64Array.from(mesh.latCell, (lat) => 272 + 30 * Math.cos(lat) ** 2), new Float64Array(C));
  const before = classContents(ocean), h0 = Float64Array.from(ocean.h), previous = Float64Array.from(ocean.h);
  let thinnest = Infinity;
  for (let s = 0; s < 200; s++) {
    ocean.eddyTransport(2700);
    for (let n = 0; n < previous.length; n++) {
      thinnest = Math.min(thinnest, ocean.h[n] - Math.min(previous[n], EPS));
      previous[n] = ocean.h[n];
    }
  }
  assert.ok(thinnest >= -6 * EDDY_SLACK, `a class fell ${-thinnest} m below its token or its thickness before the step`);
  const after = classContents(ocean);
  let changed = 0;
  for (let n = 0; n < h0.length; n++) changed = Math.max(changed, Math.abs(ocean.h[n] - h0[n]));
  for (let j = 0; j < ocean.layers; j++) assert.ok(Math.abs(after[j].volume - before[j].volume) < 1e-12 * before[j].volume && Math.abs(after[j].heat - before[j].heat) < 1e-12 * before[j].heat, `class ${j}`);
  for (let i = 0; i < C; i++) {
    let sum = 0;
    for (let j = 0; j < ocean.layers; j++) sum += ocean.h[j * C + i];
    assert.ok(Math.abs(sum - ocean.D[i] - ocean.eta[i]) < 1e-9, `column ${i}`);
  }
  assert.ok(changed > 100, `the tilted thermocline moved by at most ${changed} m`);
});
