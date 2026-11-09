import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { gridMesh } from '../js/stateGrid.module.js';
import { createOcean, depthFields, HEAT_REFERENCE, SALT_REFERENCE } from '../js/ocean/layered.module.js';
import { cpuFrameFields } from '../js/cpuFrame.module.js';
import { createGeography, decodeLandIce, geographyOptions, landIceUrl, readTopography, FRESH_GEOGRAPHY } from '../js/geography.module.js';
import { zonalWindOnEdges } from './helpers/layered.mjs';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { thinTopLayers } from './helpers/thinLayers.mjs';
import { seawaterDensity } from '../js/ocean/seawater.module.js';

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
  second.load({ ...first.serialize(), ...first.restartArrays() }, second.surfaceT, ice);
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

test(`at N=${N} a CPU ocean loads a save with T and S and no Q or W, a top layer thinner than EPS included, with the heat and salt the save's layers hold`, () => {
  const options = { ocean: { everySteps: 1 } };
  const prepare = (model) => { initializeState(model, {}).forEach((values, a) => model.state[a].set(values)); model.ocean.initialize(model.state[3], model.state[6]); return model; };
  const run = prepare(createModel(new Grid(N), options)), C = run.mesh.nCells, column = new Float64Array(C);
  run.ocean.serialize().h.forEach((v, n) => { column[n % C] += v; });
  const thin = column.indexOf(Math.max(...column)), { Q, W, state } = thinTopLayers(run.ocean, run.state[6], [[thin, 0.005]]);
  const resumed = createModel(new Grid(N), options);
  resumed.state.forEach((array, a) => array.set(run.state[a]));
  resumed.ocean.load(state, resumed.state[3], resumed.state[6]);
  for (const [name, saved, loaded] of [['Q', Q, resumed.ocean.Q], ['W', W, resumed.ocean.W]]) {
    const at = saved.findIndex((v, n) => !(Math.abs(loaded[n] - v) <= 1e-15 * Math.abs(v)));
    assert.equal(at, -1, `${name}[${at}]${at % C === thin ? ', in the thin column,' : ''} loads as ${loaded[at]}, saved as ${saved[at]}`);
  }
});

test(`at N=${N} a CPU ocean takes a layer of 1e-6 m or thinner at its class density and a thicker one at its own water`, () => {
  const options = { ocean: { everySteps: 1 } };
  const prepare = (model) => { initializeState(model, {}).forEach((values, a) => model.state[a].set(values)); model.ocean.initialize(model.state[3], model.state[6]); return model; };
  const run = prepare(createModel(new Grid(N), options)), C = run.mesh.nCells, column = new Float64Array(C);
  run.ocean.serialize().h.forEach((v, n) => { column[n % C] += v; });
  const [above, submicron, micron] = [...column.keys()].filter((i) => column[i] > 0).sort((a, b) => column[b] - column[a]);
  const { Q, W, state } = thinTopLayers(run.ocean, run.state[6], [[above, 2e-6], [submicron, 1e-7], [micron, 1e-6]]);
  const resumed = createModel(new Grid(N), options), { ocean } = resumed;
  resumed.state.forEach((array, a) => array.set(run.state[a]));
  ocean.load(state, resumed.state[3], resumed.state[6]);
  for (const i of [submicron, micron]) assert.equal(ocean.layerRho[i], ocean.densities[0], `the ${state.h[i]} m top layer at ${i} weighs ${ocean.layerRho[i]}, not its class's ${ocean.densities[0]}`);
  const own = seawaterDensity(Q[above] / state.h[above], W[above] / state.h[above]);
  assert.ok(Math.abs(ocean.layerRho[above] - own) <= 1e-12 * own, `the 2e-6 m top layer at ${above} weighs ${ocean.layerRho[above]}, its water ${own}`);
});

test(`at N=${N} the CPU frame's SSS and the CPU ocean's mean surface salinity take each sea column's W/h at the SST's layer, the first with water, over top layers of 1e-7 m and below empty ones`, () => {
  const model = createModel(new Grid(N), { ocean: { everySteps: 1 } }), { ocean, state, mesh } = model, C = mesh.nCells, L = ocean.layers, bound = 1e-12;
  initializeState(model, {}).forEach((values, a) => state[a].set(values));
  ocean.initialize(state[3], state[6]);
  const sea = [...Array(C).keys()].filter((i) => ocean.cellOcean[i]);
  for (const cut of [1e-7, 0]) {
    ocean.initialize(state[3], state[6]);
    ocean.load(thinTopLayers(ocean, state[6], sea.map((i) => [i, cut])).state, state[3], state[6]);
    const { h, W } = ocean, { sss } = cpuFrameFields(model, { depth: 'surface', fields: ['sss'] });
    const sstLayer = depthFields(mesh, L, { h, u: new Float64Array(L * mesh.nEdges), temperature: (k) => k, cellOcean: ocean.cellOcean }, 0, { upwellingDepth: null }).temperature;
    let surface = 0, thicker = 0, area = 0;
    for (const i of sea) {
      const n = sstLayer[i] * C + i;
      let m = i;
      while (!(h[m] > 1e-6)) m += C;
      assert.ok(h[n] > 0 && Object.is(sss[i], Math.fround(W[n] / h[n])), `the CPU frame's SSS at ${i}, over a ${h[i]} m top layer: ${sss[i]}; the SST's layer, ${sstLayer[i]}, holds ${W[n] / h[n]}`);
      surface += mesh.areaCell[i] * (W[n] / h[n]); thicker += mesh.areaCell[i] * (W[m] / h[m]); area += mesh.areaCell[i];
    }
    const mean = ocean.diagnostics().oceanSalinity, expected = surface / area;
    if (cut > 0) assert.ok(Math.abs(thicker / area - expected) > 10 * bound * expected, `over top layers of ${cut} m the first layer thicker than 1e-6 m gives a mean of ${thicker / area}, within ten bounds of the SST's layer's ${expected}`);
    assert.ok(Math.abs(mean - expected) <= bound * expected, `the CPU ocean's mean surface salinity over top layers of ${cut} m is ${mean}; the SST's layer gives ${expected}`);
  }
});

test(`at N=${N} a sea column with no water, whatever heat and salt its arrays still hold, gives the CPU frame no SST or SSS and takes no weight in the CPU ocean's mean surface salinity`, () => {
  const model = createModel(new Grid(N), { ocean: { everySteps: 1 } }), { ocean, state, mesh } = model, C = mesh.nCells, L = ocean.layers;
  initializeState(model, {}).forEach((values, a) => state[a].set(values));
  ocean.initialize(state[3], state[6]);
  const sea = [...Array(C).keys()].filter((i) => ocean.cellOcean[i]), emptied = sea.filter((_, n) => n % 2 === 0);
  for (const i of emptied) for (let k = 0; k < L; k++) ocean.h[k * C + i] = 0;
  const { sst, sss } = cpuFrameFields(model, { depth: 'surface', fields: ['sst', 'sss'] });
  for (const i of emptied) assert.ok(Number.isNaN(sst[i]) && Number.isNaN(sss[i]), `the column at ${i}, with no water: SST ${sst[i]} and SSS ${sss[i]} in the CPU frame`);
  let sum = 0, area = 0;
  for (const i of sea) if (ocean.h[i] > 0) { sum += mesh.areaCell[i] * (ocean.W[i] / ocean.h[i]); area += mesh.areaCell[i]; }
  const mean = ocean.diagnostics().oceanSalinity, expected = sum / area;
  assert.ok(Math.abs(mean - expected) <= 1e-12 * expected, `with ${emptied.length} of ${sea.length} sea columns empty the CPU ocean's mean surface salinity is ${mean}; the columns with water give ${expected}`);
});

test(`at N=${N} a CPU ocean converts a save on other labels with the heat and salt its layers hold: from Q and W, or from Qref and Wref as a GPU save carries them, top layers thinner than EPS and empty included, and from h·T and h·S when it carries neither`, () => {
  const options = { ocean: { everySteps: 1 } };
  const prepare = (model) => { initializeState(model, {}).forEach((values, a) => model.state[a].set(values)); model.ocean.initialize(model.state[3], model.state[6]); return model; };
  const run = prepare(createModel(new Grid(N), options)), C = run.mesh.nCells, column = new Float64Array(C);
  run.ocean.serialize().h.forEach((v, n) => { column[n % C] += v; });
  const [thin, empty] = [...column.keys()].sort((a, b) => column[b] - column[a]), { Q, W, state } = thinTopLayers(run.ocean, run.state[6], [[thin, 0.005], [empty, 0]]);
  const other = { ...state, densities: state.densities.slice(1) }, hT = Float64Array.from(state.h, (h, n) => h * state.T[n]), hS = Float64Array.from(state.h, (h, n) => h * state.S[n]);
  const Qref = Float32Array.from(Q, (q, n) => q - HEAT_REFERENCE * state.h[n]), Wref = Float32Array.from(W, (w, n) => w - SALT_REFERENCE * state.h[n]), held = (store, reference) => Float64Array.from(store, (v, n) => v + reference * state.h[n]);
  for (const [saved, heat, salt] of [[{ ...other, Q, W }, Q, W], [{ ...other, Qref, Wref }, held(Qref, HEAT_REFERENCE), held(Wref, SALT_REFERENCE)], [other, hT, hS]]) {
    const resumed = createModel(new Grid(N), options), { ocean } = resumed, L = ocean.layers;
    resumed.state.forEach((array, a) => array.set(run.state[a]));
    ocean.load(saved, resumed.state[3], resumed.state[6]);
    for (let i = 0; i < C; i++) {
      if (!(column[i] > 0)) continue;
      let depth = 0;
      for (let k = 0; k < L; k++) depth += ocean.h[k * C + i];
      for (const [name, held, loaded] of [['heat', heat, ocean.Q], ['salt', salt, ocean.W]]) {
        let before = 0, after = 0;
        for (let k = 0; k < L; k++) { before += held[k * C + i] * depth / column[i]; after += loaded[k * C + i]; }
        assert.ok(Math.abs(after - before) <= 1e-9 * before, `column ${i}, under a ${state.h[i]} m top layer, ${saved.Q ? 'with Q and W' : saved.Qref ? 'with Qref and Wref' : 'with neither'}: ${name} ${after} loaded, ${before} saved`);
      }
    }
  }
});

