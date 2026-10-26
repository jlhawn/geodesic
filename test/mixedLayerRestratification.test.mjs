import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createOcean, restratificationShape, restratificationTimescale, restratificationOptions, deformationRadii, EDDY_WAVE_SPEED } from '../js/ocean/layered.module.js';
import { seawaterDensity } from '../js/ocean/seawater.module.js';

/*
 * The mixed-layer eddy restratification of js/ocean/layered.module.js
 * (mixedLayerRestratification) on the CPU at N = 16: the vertical structure
 * and the timescale, the sign of the overturning (light water over dense),
 * the fluxes of one cold column against the formula, the conservation of
 * water, heat and salt through many applications, the depth filter, the
 * restart, and the option off leaving the ocean as it was.
 */
const mesh = buildMesh(new Grid(16)), { nCells: C, nEdges: E, cellsOnEdge, areaCell, latCell, edgesOnCell, nEdgesOnCell, maxEdges, dcEdge, dvEdge, fEdge } = mesh;
const DT = 1350, GRAVITY = 9.81, RHO0 = 1025, DEPTH = 100;

function frontOcean(options = {}, temperature = (i) => 288 + 4 * Math.tanh((latCell[i] - 0.8) / 0.12)) {
  const ocean = createOcean(mesh, { everySteps: 1, eddyDiffusivity: 0, ...options }), L = ocean.layers;
  ocean.initialize(new Float64Array(C).fill(288), new Float64Array(C));
  for (let i = 0; i < C; i++) {
    let top = 0;
    for (let k = 0; k < L && top < 1.5 * DEPTH; k++) {
      const n = k * C + i;
      ocean.Q[n] = ocean.h[n] * temperature(i); ocean.W[n] = ocean.h[n] * 35;
      top += ocean.h[n];
    }
  }
  ocean.boundaryDepths.fill(DEPTH);
  return ocean;
}
const buoyancy = (ocean, n) => -GRAVITY * seawaterDensity(ocean.Q[n] / ocean.h[n], ocean.W[n] / ocean.h[n]) / RHO0;
function stratification(ocean, cells) {
  const L = ocean.layers;
  let sum = 0, count = 0;
  for (const i of cells) {
    let z = 0, bottom = -1, depthBottom = 0, depthTop = 0.5 * ocean.h[i];
    for (let k = 0; k < L; k++) { const hk = ocean.h[k * C + i]; if (z + hk <= DEPTH) { bottom = k; depthBottom = z + 0.5 * hk; } z += hk; }
    sum += (buoyancy(ocean, i) - buoyancy(ocean, bottom * C + i)) / (depthBottom - depthTop); count++;
  }
  return sum / count;
}
function potentialEnergy(ocean) {
  let total = 0;
  for (let i = 0; i < C; i++) {
    let z = 0;
    for (let k = 0; k < ocean.layers; k++) {
      const n = k * C + i, hk = ocean.h[n];
      total -= GRAVITY * areaCell[i] * seawaterDensity(ocean.Q[n] / hk, ocean.W[n] / hk) * hk * (z + 0.5 * hk);
      z += hk;
    }
  }
  return total;
}
const totals = (ocean) => {
  let volume = 0, heat = 0, salt = 0;
  for (let n = 0; n < ocean.layers * C; n++) { const a = areaCell[n % C]; volume += a * ocean.h[n]; heat += a * ocean.Q[n]; salt += a * ocean.W[n]; }
  return { volume, heat, salt };
};
const columnWater = (ocean) => Float64Array.from({ length: C }, (_, i) => { let s = 0; for (let k = 0; k < ocean.layers; k++) s += ocean.h[k * C + i]; return s; });

test('the vertical structure is zero at the surface and at the base, one at mid-depth, and averages 44/63', () => {
  assert.equal(restratificationShape(0), 0);
  assert.equal(restratificationShape(-1), 0);
  assert.ok(Math.abs(restratificationShape(-0.5) - 1) < 1e-15);
  assert.equal(restratificationShape(-1.3), 0);
  assert.equal(restratificationShape(0.2), 0);
  let sum = 0;
  const n = 20000;
  for (let m = 0; m < n; m++) sum += restratificationShape(-(m + 0.5) / n) / n;
  assert.ok(Math.abs(sum - 44 / 63) < 1e-8, `${sum}`);
  assert.ok(Math.abs(restratificationShape(-0.25) - (1 - 0.25) * (1 + 5 / 21 * 0.25)) < 1e-15);
});

test('the timescale is a sixteenth of an inertial period over 2π without wind and lengthens as the wind mixes momentum faster', () => {
  const f = 1e-4, quiet = restratificationTimescale(f, 1e-6, 100);
  assert.ok(Math.abs(quiet * 16 * f - 1) < 1e-4, `${quiet}`);
  const mixing = 0.41 * 9.8696 * 0.01 ** 2 / (f * 100 ** 2 + 4 * 100 * 0.01);
  const exact = 0.0625 * (f + 2 * mixing) / (f * f + mixing * mixing);
  assert.ok(Math.abs(restratificationTimescale(f, 0.01, 100) / exact - 1) < 1e-8);
  assert.ok(restratificationTimescale(f, 0.01, 100) > quiet && Number.isFinite(restratificationTimescale(0, 0.01, 100)));
  assert.deepEqual(restratificationOptions(true), { coefficient: 1, frontLength: 500, decayDays: 30 });
  assert.equal(restratificationOptions(false), null);
  assert.throws(() => restratificationOptions({ frontLength: 0 }), /frontLength/);
});

test('a front restratifies: light water goes over the dense, stratification grows and potential energy falls, with water, heat and salt kept', () => {
  const ocean = frontOcean({ mixedLayerRestratification: true }), cells = Array.from({ length: C }, (_, i) => i).filter((i) => Math.abs(latCell[i] - 0.8) < 0.25);
  assert.ok(cells.length > 40);
  const before = { strat: stratification(ocean, cells), energy: potentialEnergy(ocean), total: totals(ocean), column: columnWater(ocean) };
  assert.ok(Math.abs(before.strat) < 1e-9, `the mixed layer starts uniform: ${before.strat}`);
  const history = [];
  for (let n = 0; n < 240; n++) { ocean.phases.eddy(DT); if (n % 60 === 59) history.push(stratification(ocean, cells)); }
  assert.ok(history[0] > 1e-9 && history.every((v, m) => m === 0 || v > history[m - 1]), `stratification grows: ${history}`);
  assert.ok(potentialEnergy(ocean) < before.energy, 'the front\'s potential energy falls');
  const after = totals(ocean), column = columnWater(ocean);
  for (const name of ['volume', 'heat', 'salt']) assert.ok(Math.abs(after[name] / before.total[name] - 1) < 1e-13, `${name}: ${after[name] / before.total[name] - 1}`);
  let drift = 0;
  for (let i = 0; i < C; i++) drift = Math.max(drift, Math.abs(column[i] - before.column[i]));
  assert.ok(drift < 1e-8, `a column's water changes by ${drift} m`);
});

test('one cold column among warm neighbors sends its dense water out at depth and takes light water in above, the layer fluxes of the formula', () => {
  const cold = 600;
  const ocean = frontOcean({ mixedLayerRestratification: { frontLength: 500 } }, (i) => (i === cold ? 284 : 288));
  const L = ocean.layers, h0 = Float64Array.from(ocean.h), Q0 = Float64Array.from(ocean.Q);
  ocean.phases.eddy(DT);
  const H = (i) => {
    let covered = 0, mass = 0;
    for (let k = 0; k < L && covered < DEPTH; k++) { const part = Math.min(h0[k * C + i], DEPTH - covered); mass += part * seawaterDensity(Q0[k * C + i] / h0[k * C + i], 35); covered += part; }
    return { covered, b: -GRAVITY / RHO0 * mass / covered };
  };
  const mine = H(cold);
  assert.ok(Math.abs(ocean.mleThickness[cold] - mine.covered) < 1e-9 && Math.abs(ocean.mleBuoyancy[cold] / mine.b - 1) < 1e-9, 'the cell\'s mean buoyancy and thickness');
  const radii = deformationRadii(mesh, EDDY_WAVE_SPEED), expected = new Float64Array(L);
  for (let m = 0; m < nEdgesOnCell[cold]; m++) {
    const e = edgesOnCell[maxEdges * cold + m], a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1], sign = a === cold ? 1 : -1;
    const ha = H(a), hb = H(b), depth = 0.5 * (ha.covered + hb.covered);
    const u = restratificationTimescale(Math.abs(fEdge[e]), 1e-6, depth) * Math.min(dcEdge[e], radii[e]) / 500 * dvEdge[e] / dcEdge[e] * (hb.b - ha.b) * depth * depth;
    let z = 0, above = 0;
    for (let k = 0; k < L; k++) {
      z += h0[k * C + a] / 2 + h0[k * C + b] / 2;
      const below = restratificationShape(-z / depth);
      expected[k] -= sign * u * (above - below) * DT / areaCell[cold];
      above = below;
    }
  }
  let scale = 0;
  for (let k = 0; k < L; k++) scale = Math.max(scale, Math.abs(expected[k]));
  for (let k = 0; k < L; k++) assert.ok(Math.abs(ocean.h[k * C + cold] - h0[k * C + cold] - expected[k]) < 1e-9 * scale + 1e-14, `layer ${k}: ${ocean.h[k * C + cold] - h0[k * C + cold]} against ${expected[k]}`);
  assert.ok(expected[0] > 0 && expected[8] < 0 || expected.slice(0, 8).some((v) => v > 0) && expected.slice(8, 30).some((v) => v < 0), 'the surface layers gain and the deeper ones lose');
  const warmer = (ocean.Q[0 * C + cold] / ocean.h[0 * C + cold]) - (Q0[cold] / h0[cold]);
  assert.ok(warmer > 0, `the cold column's surface layer warms by ${warmer} K`);
  const k = (() => { let z = 0; for (let n = 0; n < L; n++) { z += h0[n * C + cold]; if (z > 0.6 * DEPTH) return n; } return L - 1; })();
  assert.ok(ocean.h[k * C + cold] < h0[k * C + cold], 'its deep mixed-layer water leaves');
  assert.ok(Math.abs(ocean.Q[k * C + cold] / ocean.h[k * C + cold] - Q0[k * C + cold] / h0[k * C + cold]) < 1e-12, 'carrying its own temperature');
});

test('the donors\' allowance holds the fluxes when the coefficient is far too large, with every layer kept above zero and the books balanced', () => {
  const ocean = frontOcean({ eddyDiffusivity: 1000, mixedLayerRestratification: { coefficient: 4000 } }), before = totals(ocean), column = columnWater(ocean);
  let least = Infinity;
  for (let n = 0; n < 400; n++) {
    ocean.phases.eddy(DT);
    for (let m = 0; m < ocean.layers * C; m++) least = Math.min(least, ocean.h[m]);
  }
  assert.ok(least > -1e-9, `the thinnest layer is ${least}`);
  const after = totals(ocean), now = columnWater(ocean);
  for (const name of ['volume', 'heat', 'salt']) assert.ok(Math.abs(after[name] / before[name] - 1) < 1e-12, `${name}: ${after[name] / before[name] - 1}`);
  let drift = 0;
  for (let i = 0; i < C; i++) drift = Math.max(drift, Math.abs(now[i] - column[i]));
  assert.ok(drift < 1e-7, `a column's water changes by ${drift} m`);
});

test('the depth the restratification uses jumps to a deeper boundary layer at once and relaxes to a shallower one over decayDays', () => {
  const ocean = frontOcean({ mixedLayerRestratification: { decayDays: 10 } }), i = 700, decay = 10 * 86400;
  const keep = decay / (DT + decay);
  ocean.boundaryDepths.fill(60); ocean.phases.eddy(DT);
  assert.equal(ocean.mleDepths[i], 60);
  ocean.boundaryDepths.fill(180); ocean.phases.eddy(DT);
  assert.equal(ocean.mleDepths[i], 180, 'jumps down at once');
  let expected = 180;
  ocean.boundaryDepths.fill(40);
  for (let n = 0; n < 30; n++) { ocean.phases.eddy(DT); expected = Math.max(40, (1 - keep) * 40 + keep * expected); }
  assert.ok(Math.abs(ocean.mleDepths[i] - expected) < 1e-9 && expected > 150, `${ocean.mleDepths[i]} against ${expected}`);
  const quick = frontOcean({ mixedLayerRestratification: { decayDays: 0 } });
  quick.boundaryDepths.fill(90); quick.phases.eddy(DT); quick.boundaryDepths.fill(30); quick.phases.eddy(DT);
  assert.equal(quick.mleDepths[i], 30, 'no filter without decayDays');
});

test('a state saved with the filtered depth resumes exactly, and one saved without it starts from the boundary depth', () => {
  const run = (ocean, steps) => {
    const surfaceT = new Float64Array(C).fill(288), ice = new Float64Array(C), flux = new Float64Array(C), stress = new Float64Array(E).fill(0.02);
    for (let n = 0; n < steps; n++) { surfaceT.set(ocean.T0); ocean.advance(surfaceT, ice, flux, stress, DT); }
  };
  const options = { everySteps: 1, mixedLayerRestratification: true };
  const first = createOcean(mesh, options), surfaceT = Float64Array.from(latCell, (lat) => 300 - 26 * Math.sin(lat) ** 2), ice = new Float64Array(C);
  first.initialize(surfaceT, ice);
  run(first, 12);
  const saved = { ...first.serialize(), flux: new Float64Array(C), ...first.restartArrays() };
  assert.ok(saved.mleDepths && saved.mleDepths.length === C && saved.mleDepths.some((v) => v > 0));
  const again = createOcean(mesh, options);
  again.load(saved, Float64Array.from(first.T0), ice);
  assert.deepEqual(again.mleDepths, first.mleDepths);
  run(first, 4); run(again, 4);
  for (const name of ['h', 'u', 'Q', 'W', 'mleDepths']) assert.deepEqual(again[name], first[name], `${name} resumes exactly`);
  const { mleDepths, ...older } = saved, started = createOcean(mesh, options);
  started.load(older, Float64Array.from(first.T0), ice);
  assert.deepEqual(started.mleDepths, started.boundaryDepths);
});

test('by default the restratification runs, the same bits as the option true, and false leaves it out: no restratification cells, no saved depths and other layers', () => {
  const build = (options) => {
    const ocean = createOcean(mesh, { everySteps: 1, ...options }), surfaceT = Float64Array.from(latCell, (lat) => 300 - 26 * Math.sin(lat) ** 2), ice = new Float64Array(C), flux = new Float64Array(C), stress = new Float64Array(E);
    for (let e = 0; e < E; e++) stress[e] = 0.06 * Math.sin(2 * mesh.latEdge[e]) * Math.cos(Math.atan2(mesh.xEdge[3 * e + 1], mesh.xEdge[3 * e]));
    ocean.initialize(surfaceT, ice);
    for (let n = 0; n < 16; n++) { surfaceT.set(ocean.T0); ocean.advance(surfaceT, ice, flux, stress, DT); }
    return ocean;
  };
  const plain = build({}), on = build({ mixedLayerRestratification: true });
  for (const name of ['h', 'u', 'Q', 'W', 'eta', 'boundaryDepths']) assert.deepEqual(on[name], plain[name], name);
  assert.ok(plain.restratCells && plain.restartArrays().mleDepths, 'the default saves its restratification depths');
  const off = build({ mixedLayerRestratification: false });
  assert.equal(off.restratCells, undefined);
  assert.equal(off.restartArrays().mleDepths, undefined);
  let moved = 0;
  for (let n = 0; n < off.h.length; n++) moved = Math.max(moved, Math.abs(off.h[n] - plain.h[n]));
  assert.ok(moved > 0, 'with the option off the layers differ');
});
