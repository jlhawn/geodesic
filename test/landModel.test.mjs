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

test('with a continent the ocean flows only through ocean edges, land carries no sea ice, and the bucket accounts for every drop', async () => {
  const model = prepare(createModel(new Grid(8), { topography }));
  const { mesh, geography, land, ocean } = model;
  assert.ok(geography.landArea > 0.35 && geography.landArea < 0.55, `land area ${geography.landArea.toFixed(2)}`);
  const dt = 1350 * 16 / 8;
  const waterBefore = land.water();
  let rainOnLand = 0, evaporationFromLand = 0;
  for (let n = 0; n < 12; n++) {
    model.step(dt);
    for (let i = 0; i < mesh.nCells; i++) if (geography.land[i]) {
      rainOnLand += mesh.areaCell[i] * model.moist.rain[i];
      evaporationFromLand += mesh.areaCell[i] * model.radiation.evaporation[i] * dt;
    }
  }
  const d = model.diagnostics();
  assert.ok(d.landFraction > 0.35 && d.landMeanT > 250 && d.landMeanT < 320, `land diagnostics ${JSON.stringify({ f: d.landFraction, T: d.landMeanT })}`);
  const imbalance = land.water() - waterBefore - (rainOnLand - evaporationFromLand);
  assert.ok(Math.abs(imbalance) < 1e-9 * waterBefore, `land water balance off by ${(imbalance / waterBefore).toExponential(2)} of the store (rain ${(rainOnLand / waterBefore).toExponential(2)}, evaporation ${(evaporationFromLand / waterBefore).toExponential(2)})`);
  for (let i = 0; i < mesh.nCells; i++) if (geography.land[i]) assert.equal(model.state[6][i], 0);
  for (let e = 0; e < mesh.nEdges; e++) if (!geography.edgeOcean[e]) { assert.equal(ocean.u1[e], 0); assert.equal(ocean.u2[e], 0); }
  let moved = 0;
  for (let e = 0; e < mesh.nEdges; e++) if (geography.edgeOcean[e] && Math.abs(ocean.u1[e]) > 0) moved++;
  assert.ok(moved > 0, 'the ocean moves somewhere');
  for (let i = 0; i < mesh.nCells; i++) if (geography.land[i]) assert.equal(ocean.h1[i], 50);
});

test('the parallel engine reproduces the serial model over a continent', async () => {
  const serial = prepare(createModel(new Grid(6), { topography }));
  const parallel = prepare(await createParallelModel(new Grid(6), { topography }, 3));
  const dt = 1350 * 16 / 6;
  for (let n = 0; n < 4; n++) { serial.step(dt); await parallel.step(dt); }
  let maxT = 0, maxSoil = 0, maxTheta = 0;
  for (let i = 0; i < serial.mesh.nCells; i++) {
    maxT = Math.max(maxT, Math.abs(serial.state[3][i] - parallel.state[3][i]));
    maxSoil = Math.max(maxSoil, Math.abs(serial.land.soil[i] - parallel.land.soil[i]));
  }
  for (let x = 0; x < serial.state[1].length; x++) maxTheta = Math.max(maxTheta, Math.abs(serial.state[1][x] - parallel.state[1][x]));
  assert.ok(maxT < 1e-9 && maxSoil < 1e-9 && maxTheta < 1e-9, `serial vs parallel: Ts ${maxT} soil ${maxSoil} theta ${maxTheta}`);
  await parallel.close();
});
