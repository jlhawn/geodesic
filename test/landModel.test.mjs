import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createParallelModel } from '../js/parallel.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { initializeState } from '../js/physics/init.module.js';

const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 300 : -4000));

function prepare(model) {
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i]) model.state[6][i] = 0;
  if (model.load) model.load();
  if (model.ocean) model.ocean.initialize(model.state[3], model.state[6]);
  model.land.initialize();
  for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i]) model.land.soil[i] = 75;
  return model;
}

function landWater(model, steps = 12, dt = 1350) {
  const { mesh, geography, land } = model;
  const waterBefore = land.water();
  let rainOnLand = 0, evaporationFromLand = 0;
  for (let n = 0; n < steps; n++) {
    model.step(dt);
    for (let i = 0; i < mesh.nCells; i++) if (geography.land[i]) {
      rainOnLand += mesh.areaCell[i] * model.moist.rain[i];
      evaporationFromLand += mesh.areaCell[i] * model.radiation.evaporation[i] * dt;
    }
  }
  return { waterBefore, rainOnLand, evaporationFromLand, imbalance: land.water() - waterBefore - (rainOnLand - evaporationFromLand) };
}
const balanced = ({ waterBefore, rainOnLand, evaporationFromLand, imbalance }) => assert.ok(Math.abs(imbalance) < 1e-9 * waterBefore, `land water balance off by ${(imbalance / waterBefore).toExponential(2)} of the store (rain ${(rainOnLand / waterBefore).toExponential(2)}, evaporation ${(evaporationFromLand / waterBefore).toExponential(2)})`);

test('with a continent the ocean flows only through ocean edges, land carries no sea ice, and the bucket accounts for every drop', async () => {
  const model = prepare(createModel(new Grid(16), { topography }));
  const { mesh, geography, ocean } = model;
  assert.ok(geography.landArea > 0.35 && geography.landArea < 0.55, `land area ${geography.landArea.toFixed(2)}`);
  const water = landWater(model);
  const d = model.diagnostics();
  assert.ok(d.landFraction > 0.35 && d.landMeanT > 250 && d.landMeanT < 320, `land diagnostics ${JSON.stringify({ f: d.landFraction, T: d.landMeanT })}`);
  balanced(water);
  for (let i = 0; i < mesh.nCells; i++) if (geography.land[i]) assert.equal(model.state[6][i], 0);
  for (let k = 0; k < ocean.layers; k++) for (let e = 0; e < mesh.nEdges; e++) if (!geography.edgeOcean[e]) assert.equal(ocean.u[k * mesh.nEdges + e], 0);
  let moved = 0;
  for (let e = 0; e < mesh.nEdges; e++) if (geography.edgeOcean[e] && Math.abs(ocean.u[e]) > 0) moved++;
  assert.ok(moved > 0, 'the ocean moves somewhere');
  for (let i = 0; i < mesh.nCells; i++) if (geography.land[i]) assert.equal(ocean.h[i], 0);
});

test('with a deep store of 300 kg/m² under the trees the root zone, the surface layer, the store and the runoff account for every drop', () => {
  const model = prepare(createModel(new Grid(16), { topography, land: { deepCapacity: 300 } }));
  const { mesh, geography, land } = model, before = Float64Array.from(land.deepWater);
  const water = landWater(model);
  let moved = 0;
  for (let i = 0; i < mesh.nCells; i++) if (geography.land[i] && land.deepWater[i] !== before[i]) moved++;
  assert.ok(moved > 0, 'the trees draw on the store');
  balanced(water);
});

for (const [label, land] of [[' with the soil goal off', {}], [' with the soil goal on', { soilGoal: 'observed', soilGoalTime: 3 * 86400 }]]) test(`the change of the land's water is rain less evaporation plus the booked goal${label}`, () => {
  const model = prepare(createModel(new Grid(16), { topography, geography: { soilWaterMonthly: true }, land }));
  const { mesh, geography, land: surface } = model;
  const water = landWater(model);
  const goal = surface.budget.goal ?? 0;
  assert.ok(Math.abs(water.imbalance - goal) < 1e-9 * water.waterBefore, `change - (rain - evaporation) = ${water.imbalance}, booked ${goal}`);
  if (land.soilGoal) {
    let acted = 0;
    for (let i = 0; i < mesh.nCells; i++) if (geography.land[i] && surface.goalAdded[i] + surface.goalRemoved[i] > 0) acted++;
    assert.ok(acted > 0 && Math.abs(goal) > 1e-6 * water.waterBefore, `the goal acted on ${acted} cells and booked ${goal}`);
  } else assert.equal(goal, 0);
});

test('the parallel engine reproduces the serial model with the soil goal on, its accumulators included', async () => {
  const options = { topography, geography: { soilWaterMonthly: true }, land: { soilGoal: 'observed', coverGoal: 'observed', soilGoalTime: 3 * 86400 } };
  const serial = prepare(createModel(new Grid(16), options));
  const parallel = prepare(await createParallelModel(new Grid(16), options, 3));
  for (let n = 0; n < 4; n++) { serial.step(1350); await parallel.step(1350); }
  for (const name of ['soil', 'surface', 'goalAdded', 'goalRemoved', 'goalCover', 'goalTrees']) {
    let max = 0, size = 0;
    for (let i = 0; i < serial.mesh.nCells; i++) { max = Math.max(max, Math.abs(serial.land[name][i] - parallel.land[name][i])); size = Math.max(size, Math.abs(serial.land[name][i])); }
    assert.ok(max < 1e-9 && (name !== 'goalAdded' && name !== 'goalRemoved' || size > 0), `${name}: difference ${max}, largest ${size}`);
  }
  await parallel.close();
});

for (const [label, land] of [['', {}], [' without the deep store', { deepCapacity: 0 }]]) test(`the parallel engine reproduces the serial model over a continent${label}`, async () => {
  const serial = prepare(createModel(new Grid(16), { topography, land }));
  const parallel = prepare(await createParallelModel(new Grid(16), { topography, land }, 3));
  const dt = 1350;
  for (let n = 0; n < 4; n++) { serial.step(dt); await parallel.step(dt); }
  let maxT = 0, maxSoil = 0, maxTheta = 0, maxDeep = 0;
  for (let i = 0; i < serial.mesh.nCells; i++) {
    maxT = Math.max(maxT, Math.abs(serial.state[3][i] - parallel.state[3][i]));
    maxSoil = Math.max(maxSoil, Math.abs(serial.land.soil[i] - parallel.land.soil[i]));
    if (serial.land.deepWater) maxDeep = Math.max(maxDeep, Math.abs(serial.land.deepWater[i] - parallel.land.deepWater[i]));
  }
  for (let x = 0; x < serial.state[1].length; x++) maxTheta = Math.max(maxTheta, Math.abs(serial.state[1][x] - parallel.state[1][x]));
  assert.ok(maxT < 1e-9 && maxSoil < 1e-9 && maxTheta < 1e-9 && maxDeep < 1e-9, `serial vs parallel: Ts ${maxT} soil ${maxSoil} theta ${maxTheta} deep store ${maxDeep}`);
  await parallel.close();
});

test('under a capping inversion the stratocumulus deck forms over the warm sea and never over land', () => {
  const model = prepare(createModel(new Grid(16), { topography, radiation: { stratus: true, mixedLayerDeck: false } }));
  const { mesh, geography, core } = model, C = mesh.nCells;
  for (let k = 0; k < core.K; k++) if (core.sigmaMid[k] < 0.75) for (let i = 0; i < C; i++) model.state[1][k * C + i] += 10;
  model.step(900);
  model.step(900);
  let sea = 0, decked = 0;
  for (let i = 0; i < C; i++) {
    if (geography.land[i]) { assert.equal(model.radiation.stratus[i], 0); assert.equal(model.radiation.stratusFraction[i], 0); assert.ok(Number.isNaN(model.radiation.stabilityIndex[i])); continue; }
    sea++;
    if (!(model.radiation.stratus[i] <= 0)) {
      decked++;
      assert.ok(Number.isFinite(model.radiation.stratus[i]) && model.radiation.stratusFraction[i] <= 0.19 + 0.08 * (model.radiation.stabilityIndex[i] - 1) + 1e-12, `cover ${model.radiation.stratusFraction[i]} at EIS ${model.radiation.stabilityIndex[i]}`);
    }
  }
  assert.ok(decked > 0.1 * sea, `stratus on ${decked} of ${sea} sea cells`);
  const mixed = prepare(createModel(new Grid(16), { topography }));
  mixed.radiation.mlmSubsidence.fill(-1e-3);
  for (let k = 0; k < core.K; k++) if (core.sigmaMid[k] < 0.75) for (let i = 0; i < C; i++) mixed.state[1][k * C + i] += 10;
  mixed.step(900);
  mixed.step(900);
  let moved = 0;
  for (let i = 0; i < C; i++) {
    if (geography.land[i]) { assert.ok(mixed.radiation.stratus[i] === 0 && mixed.radiation.mlmCover[i] === 0 && mixed.radiation.mlmSubsidence[i] === -1e-3, `land cell ${i}`); continue; }
    if (mixed.radiation.mlmSubsidence[i] !== -1e-3) moved++;
  }
  assert.ok(moved > 0.5 * sea, `the mixed-layer deck's mean moved on ${moved} of ${sea} sea cells and on no land`);
});
