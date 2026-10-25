import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { gridMesh } from '../js/stateGrid.module.js';
import { createOcean } from '../js/ocean/layered.module.js';
import { createGeography, decodeLandIce, geographyOptions, landIceUrl, readTopography, FRESH_GEOGRAPHY } from '../js/geography.module.js';
import { zonalWindOnEdges } from './helpers/layered.mjs';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';

const N = 16;

test(`on ETOPO 2022 with its cavities at N=${N}, an ocean restarted from its arrays after a day steps the second day on bit for bit with one run whole`, () => {
  const mesh = gridMesh(new Grid(N)), C = mesh.nCells;
  const geography = createGeography(mesh, readTopography('etopo2022'), geographyOptions(FRESH_GEOGRAPHY, decodeLandIce(readFileSync(landIceUrl))).geography);
  const cavities = geography.cavity.reduce((n, x) => n + x, 0);
  assert.ok(cavities > 0, 'no cavity at N=16');
  const surfaceT = new Float64Array(C).fill(275), ice = new Float64Array(C);
  const wind = zonalWindOnEdges(mesh, (lat) => 0.1 * Math.cos(2 * lat)), dt = 1350, day = 64;
  const run = (ocean, steps) => { for (let n = 0; n < steps; n++) ocean.advance(ocean.surfaceT, ice, ocean.flux, wind, dt); return ocean; };
  const fresh = () => { const ocean = createOcean(mesh, { geography, everySteps: 1 }); ocean.surfaceT = Float64Array.from(surfaceT); ocean.flux = new Float64Array(C); ocean.initialize(ocean.surfaceT, ice); return ocean; };
  const whole = run(fresh(), 2 * day), first = run(fresh(), day);
  const second = createOcean(mesh, { geography, everySteps: 1 });
  second.surfaceT = Float64Array.from(first.surfaceT); second.flux = Float64Array.from(first.flux);
  second.load({ ...first.serialize(), ...first.restartArrays(), flux: Float64Array.from(first.flux) }, second.surfaceT, ice);
  run(second, day);
  const a = { ...whole.serialize(), ...whole.restartArrays() }, b = { ...second.serialize(), ...second.restartArrays() };
  for (const name of ['h', 'u', 'T', 'S', 'eta', 'shelfLoad', 'shelfMelted', 'shelfDraw', 'previousT0', 'fresh', 'buoyancy', 'boundaryDepths']) {
    assert.ok(a[name] && b[name], `${name} not saved`);
    const at = Array.from(a[name]).findIndex((x, i) => x !== b[name][i]);
    assert.equal(at, -1, `${name}[${at}]: ${a[name][at]} whole against ${b[name][at]} restarted`);
  }
  console.log(`${cavities} cavity cells; the restarted second day matches the whole run bit for bit`);
});

test(`at N=${N} a model's ocean saves the heat its top layer hands the sea ice, and a model loaded from it gives the ice the same heat before its first ocean step`, () => {
  const options = { ocean: { everySteps: 1 } };
  const prepare = (model) => { initializeState(model, {}).forEach((values, a) => model.state[a].set(values)); model.ocean.initialize(model.state[3], model.state[6]); return model; };
  const run = prepare(createModel(new Grid(N), options));
  for (let n = 0; n < 4; n++) run.step(1350);
  const saved = { ...run.ocean.serialize(), ...run.ocean.restartArrays() };
  assert.ok(saved.flux && saved.flux.some((x) => x > 0), 'the ocean hands the ice some heat');
  const resumed = createModel(new Grid(N), options);
  resumed.state.forEach((array, a) => array.set(run.state[a]));
  resumed.seaIce.load(resumed.state[6], run.seaIce.concentration);
  resumed.ocean.load(saved, resumed.state[3], resumed.state[6]);
  assert.deepEqual(Array.from(resumed.seaIce.oceanFlux), Array.from(run.seaIce.oceanFlux));
});
