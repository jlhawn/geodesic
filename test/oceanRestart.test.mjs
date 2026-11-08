import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { gridMesh } from '../js/stateGrid.module.js';
import { createOcean, HEAT_REFERENCE, SALT_REFERENCE } from '../js/ocean/layered.module.js';
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

test(`at N=${N} a model's ocean loaded from a GPU save, whose contents are stored as Qref and Wref beside float32 T and S, takes its heat and salt from Qref and Wref`, () => {
  const options = { ocean: { everySteps: 1 } };
  const prepare = (model) => { initializeState(model, {}).forEach((values, a) => model.state[a].set(values)); model.ocean.initialize(model.state[3], model.state[6]); return model; };
  const run = prepare(createModel(new Grid(N), options));
  for (let n = 0; n < 2; n++) run.step(1350);
  const { Q, W, ...rest } = run.ocean.restartArrays(), served = run.ocean.serialize(), h = served.h;
  const Qref = Float32Array.from(Q, (q, n) => q - HEAT_REFERENCE * h[n]), Wref = Float32Array.from(W, (w, n) => w - SALT_REFERENCE * h[n]);
  const saved = { ...served, T: Float32Array.from(served.T), S: Float32Array.from(served.S), ...rest, Qref, Wref };
  const resumed = createModel(new Grid(N), options);
  resumed.state.forEach((array, a) => array.set(run.state[a]));
  resumed.ocean.load(saved, resumed.state[3], resumed.state[6]);
  const at = (n) => [Qref[n] + HEAT_REFERENCE * h[n], Wref[n] + SALT_REFERENCE * h[n]];
  for (let n = 0; n < Q.length; n++) assert.deepEqual([resumed.ocean.Q[n], resumed.ocean.W[n]], at(n), `layer cell ${n}`);
});

test(`on ETOPO 2022 at N=${N}, an ocean saved between its steps as the page saves it, every array as float64, steps on in a fresh ocean bit for bit with one run whole`, () => {
  const mesh = gridMesh(new Grid(N)), C = mesh.nCells, every = 4, dt = 1350, before = 2 * every + 2, after = 2 * every - 1;
  const geography = createGeography(mesh, readTopography('etopo2022'), geographyOptions(FRESH_GEOGRAPHY, decodeLandIce(readFileSync(landIceUrl))).geography);
  const wind = zonalWindOnEdges(mesh, (lat) => 0.1 * Math.cos(2 * lat));
  const evaporation = Float64Array.from({ length: C }, (_, i) => 3e-5 * (1 + 0.5 * Math.sin(3 * i)));
  const rain = Float64Array.from({ length: C }, (_, i) => 0.04 * (1 + 0.5 * Math.cos(5 * i)));
  const discharge = Float64Array.from({ length: C }, (_, i) => (geography.land[i] ? 1e-3 : 0));
  const polar = Float64Array.from({ length: C }, (_, i) => (Math.abs(mesh.latCell[i]) > 1.1 && !geography.land[i] ? 1 : 0));
  const step = (ocean, n) => {
    for (let i = 0; i < C; i++) { ocean.ice[i] = polar[i] * (0.5 + 0.02 * n); if (ocean.transmission) ocean.penetrated[i] += polar[i] * dt * 2 * ocean.transmission[i]; }
    ocean.accumulate(evaporation, rain, dt, null, discharge);
    ocean.advance(ocean.surfaceT, ocean.ice, ocean.flux, wind, dt);
  };
  const fresh = () => {
    const ocean = createOcean(mesh, { geography, everySteps: every });
    ocean.surfaceT = new Float64Array(C).fill(275); ocean.ice = Float64Array.from(polar, (p) => 0.5 * p); ocean.flux = new Float64Array(C);
    ocean.initialize(ocean.surfaceT, ocean.ice);
    return ocean;
  };
  const whole = fresh(), first = fresh();
  for (let n = 0; n < before + after; n++) step(whole, n);
  for (let n = 0; n < before; n++) step(first, n);
  const saved = Object.fromEntries(Object.entries(first.serialize({ restart: true })).map(([name, v]) => [name, typeof v === 'object' && v.length && typeof v[0] === 'number' ? Float64Array.from(v) : v]));
  assert.equal(saved.phase, before % every, 'the save falls between the ocean\'s steps');
  for (const name of ['fresh', 'calved', 'penetrated']) assert.ok(saved[name].some((x) => x !== 0), `the save carries no ${name}`);
  const second = createOcean(mesh, { geography, everySteps: every });
  second.surfaceT = Float64Array.from(first.surfaceT); second.ice = Float64Array.from(first.ice); second.flux = Float64Array.from(saved.flux);
  second.load(saved, second.surfaceT, second.ice);
  for (let n = before; n < before + after; n++) step(second, n);
  const a = { ...whole.serialize({ restart: true }), oceanFlux: whole.flux }, b = { ...second.serialize({ restart: true }), oceanFlux: second.flux };
  assert.equal(b.phase, a.phase);
  for (const [name, values] of Object.entries(a)) if (typeof values === 'object' && values.length && typeof values[0] === 'number') {
    const at = Array.from(values).findIndex((x, i) => !Object.is(x, b[name][i]));
    assert.equal(at, -1, `${name}[${at}]: ${values[at]} whole against ${b[name][at]} restored`);
  }
});
