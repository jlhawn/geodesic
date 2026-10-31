import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createOcean, depthFields } from '../js/ocean/layered.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { upwellingTimescale, upwellingWeight, SURFACE_UPWELLING_DEPTH } from '../js/frames.module.js';
import { layerVelocities, velocityAt, worstDifference, unevenTopLayers } from './helpers/upwelling.mjs';

/*
 * The CPU ocean's upwelling memory at N=16 on a sea with a continent and
 * a shelf: after a fold its vertical velocity at any depth is the
 * divergence of the model's own layer mass fluxes, taken from the
 * ocean's tendency, interpolated within the layer that holds the depth;
 * successive folds weight each layer by the timescale at its base; and
 * the memory folds only when the ocean has stepped and restarts when time
 * goes back or a state loads.
 */
const N = 16, DT = 1350;
const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 300 + 2500 * Math.exp(-(((lat - 0.3) / 0.3) ** 2)) : -4000));

function scenario(options = {}) {
  const model = createModel(new Grid(N), { topography });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  const { mesh } = model, C = mesh.nCells, E = mesh.nEdges;
  for (let i = 0; i < C; i++) if (model.geography.land[i]) model.state[6][i] = 0;
  const surfaceT = Float64Array.from(model.state[3]), ice = Float64Array.from(model.state[6]), stress = new Float64Array(E), flux = new Float64Array(C);
  for (let e = 0; e < E; e++) { const lon = Math.atan2(mesh.xEdge[3 * e + 1], mesh.xEdge[3 * e]); stress[e] = 0.08 * Math.sin(2 * mesh.latEdge[e]) * Math.cos(lon) + 0.01 * Math.sin(5 * lon); }
  const ocean = createOcean(mesh, { geography: model.geography, vertical: 'hybrid', everySteps: 1, iceDrift: false, ...options });
  ocean.initialize(Float64Array.from(surfaceT), ice);
  const advance = (n) => { for (let m = 0; m < n; m++) ocean.advance(Float64Array.from(surfaceT), ice, flux, stress, DT); };
  const state = () => ({ h: ocean.h, u: ocean.u, eta: ocean.eta, sillEdge: ocean.sillEdge, edgeOcean: ocean.edgeOcean });
  return { model, mesh, C, E, L: ocean.layers, ocean, surfaceT, ice, advance, state };
}

// The vertical velocity at the base of each layer from the thickness tendency alone: the layers' mass fluxes are what the model moves water with, so −dh/dt summed from the top is the flow through each base.
function fromTendency(ocean, C, L) {
  const out = ocean.stages[0];
  ocean.tendency(ocean.state, out);
  const velocities = new Float64Array(L * C);
  for (let i = 0; i < C; i++) { let sum = 0; for (let k = 0; k < L; k++) { sum -= out[0][k * C + i]; velocities[k * C + i] = sum; } }
  return velocities;
}

test('after one fold the memory gives the vertical velocity of the model\'s own mass-flux divergence at any depth', () => {
  const { mesh, C, L, ocean, advance } = scenario();
  advance(8);
  const velocities = fromTendency(ocean, C, L);
  ocean.foldUpwelling(8 * DT);
  for (const depth of [20, 137, 903, 2500]) {
    const expected = velocityAt(C, L, ocean.h, ocean.cellOcean, velocities, depth), got = ocean.fields(depth).upwelling;
    const { max, scale } = worstDifference(expected, got);
    assert.ok(scale > 1e-7, `the flow through ${depth} m reaches ${scale} m/s`);
    assert.ok(max <= 1e-6 * scale, `at ${depth} m the memory parts from the model's divergence by ${max} against ${scale}`);
  }
  const temperature = () => 0;
  const mean = depthFields(mesh, L, { h: ocean.h, u: ocean.u, temperature, cellOcean: ocean.cellOcean }, 2500).upwelling;
  const { max, scale } = worstDifference(velocityAt(C, L, ocean.h, ocean.cellOcean, velocities, 2500), mean);
  assert.ok(max > 0.2 * scale, `the flow built on the mean of the two cells' thicknesses parts from the model's by only ${max} against ${scale} at 2500 m`);

  unevenTopLayers(ocean.h, C, L, ocean.cellOcean);
  const uneven = fromTendency(ocean, C, L);
  ocean.foldUpwelling(0);
  for (const depth of [5, 137]) {
    const { max: gap, scale: size } = worstDifference(velocityAt(C, L, ocean.h, ocean.cellOcean, uneven, depth), ocean.fields(depth).upwelling);
    assert.ok(size > 1e-7 && gap <= 1e-6 * size, `with uneven top layers the memory parts from the model's divergence by ${gap} against ${size} at ${depth} m`);
  }
});

test('land and cells whose floor lies above the depth are NaN, and the surface view is the flow through 50 m', () => {
  const { C, L, ocean, advance } = scenario();
  advance(4);
  const fifty = ocean.fields(SURFACE_UPWELLING_DEPTH).upwelling, surface = ocean.fields(0, SURFACE_UPWELLING_DEPTH);
  let land = 0, shown = 0;
  for (let i = 0; i < C; i++) {
    if (!ocean.cellOcean[i]) { land++; assert.ok(Number.isNaN(surface.upwelling[i]) && Number.isNaN(fifty[i]), `land at ${i}`); continue; }
    assert.ok(Number.isFinite(surface.upwelling[i]), `the surface view at sea cell ${i}`);
    if (Number.isFinite(fifty[i])) { shown++; assert.equal(surface.upwelling[i], fifty[i], `the surface view is the flow through 50 m at ${i}`); }
  }
  assert.ok(land > 0 && shown > 0);
  const deep = ocean.fields(2500);
  let dry = 0, wet = 0;
  for (let i = 0; i < C; i++) {
    assert.equal(Number.isNaN(deep.upwelling[i]), Number.isNaN(deep.temperature[i]), `the vertical velocity and the temperature are NaN in the same cells at ${i}`);
    if (!ocean.cellOcean[i]) continue;
    let floor = 0;
    for (let k = 0; k < L; k++) floor += ocean.h[k * C + i];
    if (floor <= 2500) { dry++; assert.ok(Number.isNaN(deep.upwelling[i]), `the floor at ${floor} m lies above 2500 m at ${i}`); } else { wet++; assert.ok(Number.isFinite(deep.upwelling[i])); }
  }
  assert.ok(dry > 0 && wet > 0, `${dry} sea cells with a floor above 2500 m and ${wet} below`);
});

test('the memory folds when the ocean has stepped, weights each layer by the timescale at its base, and restarts when time goes back or a state loads', () => {
  const { mesh, C, L, ocean, advance, surfaceT, ice, state } = scenario();
  const velocities = () => layerVelocities(mesh, L, state(), ocean.cellOcean);
  const check = (expected, label, depth = 500) => {
    const { max, scale } = worstDifference(velocityAt(C, L, ocean.h, ocean.cellOcean, expected, depth), ocean.fields(depth).upwelling);
    assert.ok(scale > 1e-7 && max <= 1e-6 * scale, `${label}: ${max} against ${scale} at ${depth} m`);
  };
  advance(4);
  const first = velocities();
  ocean.foldUpwelling(1e6);
  check(first, 'the first fold replaces the memory');
  ocean.foldUpwelling(2e6);
  check(first, 'a fold with no step since the last leaves the memory as it was');

  advance(1);
  const second = velocities(), mean = first.map((v, n) => 0.5 * (v + second[n]));
  ocean.foldUpwelling(1e6 + DT);
  check(mean, 'two frames a step apart average to their mean');
  advance(1);
  const third = velocities(), gap = 40 * 86400;
  const folded = Float64Array.from(mean, (m, n) => {
    const k = Math.floor(n / C), i = n % C;
    if (!ocean.cellOcean[i]) return 0;
    let base = 0;
    for (let j = 0; j <= k; j++) base += ocean.h[j * C + i];
    return m + upwellingWeight(gap, 2, upwellingTimescale(base)) * (third[n] - m);
  });
  ocean.foldUpwelling(1e6 + DT + gap);
  for (const depth of [20, 100, 500, 1500]) check(folded, `after a gap of 40 days each layer takes the weight of its own timescale`, depth);

  ocean.foldUpwelling(0);
  check(third, 'the memory restarts when time goes back');
  advance(1);
  ocean.foldUpwelling(DT);
  const loaded = ocean.serialize();
  ocean.load(loaded, Float64Array.from(surfaceT), ice);
  ocean.foldUpwelling(2 * DT);
  check(velocities(), 'the memory restarts when a state loads');
  advance(1);
  ocean.foldUpwelling(3 * DT);
  ocean.initialize(Float64Array.from(surfaceT), ice);
  for (const value of ocean.fields(500).upwelling) assert.ok(Number.isNaN(value) || value === 0, 'the memory restarts when the ocean initializes, from rest');
});

test('the model folds its ocean at the frame time, and not again until the ocean has stepped', () => {
  const { model } = scenario();
  assert.equal(typeof model.foldUpwelling, 'function');
  const { ocean } = model;
  assert.ok(ocean, 'the model carries its ocean');
  const every = ocean.everySteps;
  for (let n = 0; n < every; n++) model.step(DT);
  model.foldUpwelling(model.time);
  const folded = model.oceanFields(500).upwelling;
  assert.ok(folded.some((v) => Number.isFinite(v) && v !== 0), 'the ocean has moved water');
  for (let n = 0; n < every - 1; n++) model.step(DT);
  model.foldUpwelling(model.time);
  assert.deepEqual(model.oceanFields(500).upwelling, folded, 'no ocean step, no new fold');
  model.step(DT);
  model.foldUpwelling(model.time);
  assert.notDeepEqual(model.oceanFields(500).upwelling, folded, 'an ocean step, a new fold');
});
