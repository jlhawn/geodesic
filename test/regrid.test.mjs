import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { interpolationWeights, regridState } from '../js/physics/regrid.module.js';
import { P0 } from '../js/dynamics/sigmaCore.module.js';

const coarse = createModel(new Grid(8));
const fine = createModel(new Grid(16));
const analytic = (x, y, z) => Math.sin(2 * Math.atan2(z, Math.hypot(x, y))) * Math.cos(Math.atan2(y, x)) + 0.5 * z;

test('interpolating a mesh onto itself is the identity', () => {
  const { cells, weights } = interpolationWeights(coarse.mesh, coarse.mesh.xCell);
  for (let i = 0; i < coarse.mesh.nCells; i++) {
    let value = 0;
    for (let m = 0; m < 3; m++) value += weights[3 * i + m] * (cells[3 * i + m] === i ? 1 : 0);
    assert.ok(Math.abs(value - 1) < 1e-9);
  }
});

test('a smooth field regridded from N=8 to N=16 matches the analytic field to second order', () => {
  const sampled = Float64Array.from({ length: coarse.mesh.nCells }, (_, i) => analytic(coarse.mesh.xCell[3 * i], coarse.mesh.xCell[3 * i + 1], coarse.mesh.xCell[3 * i + 2]));
  const { cells, weights } = interpolationWeights(coarse.mesh, fine.mesh.xCell);
  let worst = 0;
  for (let i = 0; i < fine.mesh.nCells; i++) {
    let value = 0;
    for (let m = 0; m < 3; m++) value += weights[3 * i + m] * sampled[cells[3 * i + m]];
    worst = Math.max(worst, Math.abs(value - analytic(fine.mesh.xCell[3 * i], fine.mesh.xCell[3 * i + 1], fine.mesh.xCell[3 * i + 2])));
  }
  console.log(`regrid N=8→16 max error ${worst.toExponential(2)} on a field of amplitude ~1.5`);
  assert.ok(worst < 0.06);
});

test('a state regridded to a finer mesh keeps its mass, its profile, and its solid-body rotation', () => {
  const K = coarse.core.K, C = coarse.mesh.nCells, E = coarse.mesh.nEdges;
  const omega = 1e-5;
  const pi = Float64Array.from({ length: C }, (_, i) => P0 + 1000 * coarse.mesh.xCell[3 * i + 2]);
  const theta = new Float64Array(K * C);
  for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) theta[k * C + i] = 300 + 200 * (1 - coarse.core.sigmaMid[k]) + 5 * coarse.mesh.xCell[3 * i + 2];
  const u = new Float64Array(K * E);
  for (let k = 0; k < K; k++) for (let e = 0; e < E; e++) u[k * E + e] = coarse.mesh.radius * omega * (-coarse.mesh.xEdge[3 * e + 1] * coarse.mesh.nEdge[3 * e] + coarse.mesh.xEdge[3 * e] * coarse.mesh.nEdge[3 * e + 1]);
  const surfaceT = Float64Array.from({ length: C }, (_, i) => 288 - 30 * coarse.mesh.xCell[3 * i + 2] ** 2);
  const [fPi, fTheta, fU, fSurfaceT] = regridState(coarse, fine, [pi, theta, u, surfaceT]);
  const fm = fine.mesh;
  let massCoarse = 0, massFine = 0, area = 0, worstU = 0, worstPi = 0;
  for (let i = 0; i < C; i++) massCoarse += coarse.mesh.areaCell[i] * pi[i];
  for (let i = 0; i < fm.nCells; i++) { massFine += fm.areaCell[i] * fPi[i]; area += fm.areaCell[i]; worstPi = Math.max(worstPi, Math.abs(fPi[i] - (P0 + 1000 * fm.xCell[3 * i + 2]))); }
  for (let e = 0; e < fm.nEdges; e++) {
    const exact = fm.radius * omega * (-fm.xEdge[3 * e + 1] * fm.nEdge[3 * e] + fm.xEdge[3 * e] * fm.nEdge[3 * e + 1]);
    worstU = Math.max(worstU, Math.abs(fU[(K - 1) * fm.nEdges + e] - exact));
  }
  console.log(`regrid state: mass ratio ${(massFine / massCoarse).toFixed(6)}, max ps error ${worstPi.toFixed(2)} Pa, max wind error ${worstU.toFixed(3)} m/s of ${(fm.radius * omega).toFixed(1)}`);
  assert.ok(Math.abs(massFine / massCoarse - 1) < 1e-3);
  assert.ok(worstPi < 20);
  assert.ok(worstU < 0.05 * fm.radius * omega);
  assert.ok(fTheta.every((t) => t > 250 && t < 600) && fSurfaceT.every((t) => t > 250 && t < 300));
});

test('tile sampling keeps a step field exact and, masked, takes a neighbour or the guess', async () => {
  const { sampleTiles, regridCellField, landFromSea, seaFromLand } = await import('../js/physics/regrid.module.js');
  const sm = coarse.mesh;
  const step = Float64Array.from({ length: sm.nCells }, (_, i) => (sm.latCell[i] > 0 ? 1 : 0));
  const tiles = sampleTiles(coarse, fine, step);
  assert.ok(tiles.every((v) => v === 0 || v === 1), 'a tile-sampled step field has no intermediate values');
  const south = Uint8Array.from(step, (v) => 1 - v);
  const poisoned = Float64Array.from(step, (v, i) => (south[i] ? 10 + sm.latCell[i] : -100));
  const masked = sampleTiles(coarse, fine, poisoned, south, undefined, () => -1);
  assert.ok(masked.every((v) => v > 0 || v === -1), 'masked sampling takes admitted tiles or the guess');
  assert.ok(masked.some((v) => v === -1), 'points with no admitted neighbour take the guess');
  for (let n = 0; n < masked.length; n++) if (fine.mesh.latCell[n] > 0 && fine.mesh.latCell[n] < 0.03) assert.ok(masked[n] > 0, 'points just across the boundary take the adjacent admitted tile');
  const interpolated = regridCellField(coarse, fine, poisoned, undefined, south);
  for (let n = 0; n < interpolated.length; n++) if (fine.mesh.latCell[n] < 0.15) assert.ok(interpolated[n] > 0, 'masked interpolation draws only on admitted cells within reach');
  assert.deepEqual(landFromSea(0, 290), { snow: 0, soil: 75 }, 'warm open sea implies bare land with a half-full bucket');
  assert.deepEqual(landFromSea(0.5, 260), { snow: 100, soil: 150 }, 'thick sea ice implies a deep snow cover on frozen, saturated ground');
  assert.ok(landFromSea(0.01, 265).snow >= 20, 'any sea ice implies enough snow for the full snow albedo');
  assert.equal(seaFromLand(0, 290), 0, 'bare warm land implies open water');
  assert.equal(seaFromLand(100, 265), 0.5, 'a deep snow pack implies first-year ice');
  assert.equal(seaFromLand(0, 270), 0.1, 'bare land below the seawater freezing point implies new thin ice');
});

test('land regridding carries the vegetation cover between resolutions and leaves a state without it without it', async () => {
  const { regridLand } = await import('../js/physics/regrid.module.js');
  const { syntheticTopography } = await import('../js/geography.module.js');
  const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 ? 300 : -4000));
  const source = createModel(new Grid(6), { topography }), target = createModel(new Grid(10), { topography });
  const C = source.mesh.nCells, cover = Float64Array.from({ length: C }, (_, i) => (source.mesh.latCell[i] > 0 ? 0.8 : 0.3));
  const land = { soil: new Float64Array(C).fill(100), snow: new Float64Array(C), vegetation: cover };
  const out = regridLand(source, target, land);
  let inland = 0;
  for (let n = 0; n < target.mesh.nCells; n++) {
    if (!target.geography.land[n]) { assert.equal(out.vegetation[n], 0); continue; }
    const lat = target.mesh.latCell[n];
    if (Math.abs(Math.cos(target.mesh.lonCell[n])) > 0.4 && Math.abs(lat) > 0.1) { inland++; assert.equal(out.vegetation[n], lat > 0 ? 0.8 : 0.3); }
  }
  assert.ok(inland > 50);
  assert.equal(regridLand(source, target, { soil: land.soil, snow: land.snow }).vegetation, undefined);
});

test('sea-ice concentration regrids with its ice, from the same tiles, and ice inferred from land covers its cell', async () => {
  const { regridConcentration } = await import('../js/physics/regrid.module.js');
  const { syntheticTopography } = await import('../js/geography.module.js');
  const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 300 : -4000));
  const source = createModel(new Grid(6), { topography }), target = createModel(new Grid(10), { topography });
  const { K } = source.core, C = source.mesh.nCells, E = source.mesh.nEdges;
  const sea = (i) => !source.geography.land[i];
  const ice = Float64Array.from({ length: C }, (_, i) => (sea(i) && Math.abs(source.mesh.latCell[i]) > 1.1 ? 1 + (i % 3) : 0));
  const concentration = Float64Array.from(ice, (h, i) => (h > 0 ? 0.2 + 0.2 * (i % 4) : 0));
  const surfaceT = Float64Array.from({ length: C }, (_, i) => (ice[i] > 0 ? 271.35 : 280));
  const land = { soil: new Float64Array(C), snow: Float64Array.from({ length: C }, (_, i) => (source.mesh.latCell[i] > 1.0 ? 50 : 0)) };
  const state = [new Float64Array(C).fill(P0), new Float64Array(K * C).fill(280), new Float64Array(K * E), surfaceT, new Float64Array(K * C), new Float64Array(K * C), ice];
  const outIce = regridState(source, target, state, null, { land })[6];
  const out = regridConcentration(source, target, concentration, { land, surfaceT });
  let partial = 0;
  for (let n = 0; n < target.mesh.nCells; n++) {
    if (target.geography.land[n]) continue;
    assert.equal(outIce[n] > 0, out[n] > 0, `cell ${n}: ${outIce[n]} m of ice over ${out[n]}`);
    if (out[n] > 0 && out[n] < 1) partial++;
  }
  assert.ok(partial > 0, 'the partial cover comes along');
});

test('the mixed-layer deck\'s running-mean subsidence comes back from a saved state: kept at the same resolution, interpolated over the sea at another, 0 over land and for a state saved without it', async () => {
  const { savedSubsidence } = await import('../js/physics/regrid.module.js');
  const { syntheticTopography } = await import('../js/geography.module.js');
  const { encodeState, decodeState } = await import('../js/stateFile.module.js');
  const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 300 : -4000));
  const source = createModel(new Grid(6), { topography }), target = createModel(new Grid(10), { topography });
  const C = source.mesh.nCells, sea = (model, i) => !model.geography.land[i];
  const mean = Float64Array.from({ length: C }, (_, i) => (sea(source, i) ? -1e-3 * (1 + Math.sin(source.mesh.latCell[i])) : 0));
  const saved = await decodeState(encodeState({ N: 6, K: source.core.K, day: 1, time: 86400, mlmSubsidence: mean }));
  const kept = savedSubsidence(saved, source);
  for (let i = 0; i < C; i++) assert.equal(kept[i], Math.fround(mean[i]), `cell ${i}`);
  const moved = savedSubsidence(saved, target, source);
  let seaCells = 0;
  for (let n = 0; n < target.mesh.nCells; n++) {
    if (!sea(target, n)) { assert.equal(moved[n], 0, `land cell ${n}`); continue; }
    seaCells++;
    assert.ok(moved[n] < 0 && moved[n] >= -2e-3 - 1e-9, `sea cell ${n}: ${moved[n]}`);
  }
  assert.ok(seaCells > 0);
  const legacy = savedSubsidence(await decodeState(encodeState({ N: 6, K: source.core.K, day: 1, time: 86400, pi: new Float64Array(C) })), source);
  assert.ok(legacy.length === C && legacy.every((x) => x === 0), 'a state saved without it starts from 0');
  assert.ok(savedSubsidence(null, target).every((x) => x === 0), 'a fresh start starts from 0');
});

test('the deck\'s carried height and gate survive a saved state: a model\'s fields come back through the binary file as saved, interpolated over the sea at another resolution with their starting values over land, and a state saved without them starts from an unset height and an undecided gate, never NaN', async () => {
  const { savedDeckField, DECK_FIELDS } = await import('../js/physics/regrid.module.js');
  const { syntheticTopography } = await import('../js/geography.module.js');
  const { encodeState, decodeState } = await import('../js/stateFile.module.js');
  const { initializeState } = await import('../js/physics/init.module.js');
  assert.deepEqual(DECK_FIELDS, { mlmSubsidence: 0, mlmHeight: 0, mlmGate: 0.5 });
  const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 300 : -4000));
  const source = createModel(new Grid(6), { topography, ocean: false, radiation: { stratusSubsidence: 0, minimumInversion: 0 } }), target = createModel(new Grid(10), { topography, ocean: false });
  const C = source.mesh.nCells, sea = (model, i) => !model.geography.land[i];
  assert.ok(source.radiation.mlmHeight.every((h) => h === 0) && source.radiation.mlmGate.every((x) => x === 0.5), 'a fresh model starts unset and undecided');
  initializeState(source, {}).forEach((values, a) => source.state[a].set(values));
  for (let n = 0; n < 4; n++) source.step(900);
  const { mlmSubsidence, mlmHeight, mlmGate } = source.radiation;
  let carried = 0;
  for (let i = 0; i < C; i++) if (mlmHeight[i] > 0) carried++;
  assert.ok(carried > 0.1 * C, `the deck carries a height on ${carried} of ${C} cells`);
  const [pi, theta, u, surfaceT, q, qc, ice] = source.state;
  const bytes = encodeState({ N: 6, K: source.core.K, day: 1, time: 3600, pi, theta, u, surfaceT, q, qc, ice, mlmSubsidence, mlmHeight, mlmGate });
  const saved = await decodeState(bytes);
  const restored = createModel(new Grid(6), { topography, ocean: false });
  for (const name of Object.keys(DECK_FIELDS)) restored.radiation[name].set(savedDeckField(saved, name, restored));
  for (let i = 0; i < C; i++) {
    assert.equal(restored.radiation.mlmHeight[i], Math.fround(mlmHeight[i]), `height of cell ${i}`);
    assert.equal(restored.radiation.mlmGate[i], Math.fround(mlmGate[i]), `gate of cell ${i}`);
  }
  const exact = await decodeState(encodeState({ N: 6, mlmHeight, mlmGate }, { f64: ['mlmHeight', 'mlmGate'] }));
  assert.deepEqual(savedDeckField(exact, 'mlmHeight', source), Float64Array.from(mlmHeight));
  assert.deepEqual(savedDeckField(exact, 'mlmGate', source), Float64Array.from(mlmGate));
  let highest = 0;
  for (let i = 0; i < C; i++) if (sea(source, i)) highest = Math.max(highest, Math.fround(mlmHeight[i]));
  for (const name of ['mlmHeight', 'mlmGate']) {
    const moved = savedDeckField(saved, name, target, source);
    assert.equal(moved.length, target.mesh.nCells);
    for (let n = 0; n < moved.length; n++) {
      assert.ok(Number.isFinite(moved[n]), `${name} of target cell ${n}`);
      if (!sea(target, n)) assert.equal(moved[n], DECK_FIELDS[name], `${name} of land cell ${n}`);
      else if (name === 'mlmGate') assert.ok(moved[n] >= 0 && moved[n] <= 1, `gate of sea cell ${n}: ${moved[n]}`);
      else assert.ok(moved[n] >= 0 && moved[n] <= highest * (1 + 1e-12), `height of sea cell ${n}: ${moved[n]} above the highest saved over the sea, ${highest}`);
    }
  }
  const legacy = await decodeState(encodeState({ N: 6, K: source.core.K, day: 1, time: 3600, pi, mlmSubsidence }));
  assert.ok(savedDeckField(legacy, 'mlmHeight', source).every((h) => h === 0), 'a state saved without a height starts unset');
  assert.ok(savedDeckField(legacy, 'mlmGate', source).every((x) => x === 0.5), 'a state saved without a gate starts undecided');
  assert.ok(savedDeckField(null, 'mlmGate', target).every((x) => x === 0.5) && savedDeckField(null, 'mlmHeight', target).every((h) => h === 0), 'so does a fresh start');
  const resumed = createModel(new Grid(6), { topography, ocean: false, radiation: { stratusSubsidence: 0, minimumInversion: 0 } });
  [pi, theta, u, surfaceT, q, qc, ice].forEach((values, a) => resumed.state[a].set(values));
  for (const name of Object.keys(DECK_FIELDS)) resumed.radiation[name].set(savedDeckField(legacy, name, resumed));
  resumed.step(900); resumed.step(900);
  assert.ok([...resumed.radiation.mlmHeight, ...resumed.radiation.mlmGate, ...resumed.radiation.mlmCover].every(Number.isFinite), 'a legacy state steps without NaN');
});

test('the per-cell convective and large-scale rain survive a saved state: as saved at the same resolution, interpolated and never negative at another, and zero for a state saved without them', async () => {
  const { savedRainField, RAIN_FIELDS } = await import('../js/physics/regrid.module.js');
  const { encodeState, decodeState } = await import('../js/stateFile.module.js');
  const { initializeState } = await import('../js/physics/init.module.js');
  assert.deepEqual(RAIN_FIELDS, ['convectiveRain', 'largeScaleRain']);
  const source = createModel(new Grid(6), { ocean: false }), target = createModel(new Grid(10), { ocean: false });
  initializeState(source, {}).forEach((values, a) => source.state[a].set(values));
  for (let n = 0; n < 24; n++) source.step(900);
  source.diagnostics();
  const { convectiveRain, largeScaleRain } = source.moist, C = source.mesh.nCells;
  const area = (model, values) => { let s = 0, a = 0; for (let i = 0; i < model.mesh.nCells; i++) { s += model.mesh.areaCell[i] * values[i]; a += model.mesh.areaCell[i]; } return s / a; };
  assert.ok(area(source, convectiveRain) > 0.1, `six hours rain ${area(source, convectiveRain)} mm/d convectively`);
  const saved = await decodeState(encodeState({ N: 6, K: source.core.K, day: 0, time: source.time, pi: source.state[0], convectiveRain, largeScaleRain }));
  for (const name of RAIN_FIELDS) {
    const back = savedRainField(saved, name, source), moved = savedRainField(saved, name, target, source);
    for (let i = 0; i < C; i++) assert.equal(back[i], Math.fround(source.moist[name][i]), `${name} of cell ${i}`);
    assert.equal(moved.length, target.mesh.nCells);
    assert.ok(moved.every((x) => x >= 0 && Number.isFinite(x)), `${name} at N=10 is finite and never negative`);
  }
  const before = area(source, convectiveRain), after = area(target, savedRainField(saved, 'convectiveRain', target, source));
  assert.ok(Math.abs(after - before) < 0.2 * before, `convective rain: mean ${before} mm/d at N=6, ${after} at N=10`);
  const legacy = await decodeState(encodeState({ N: 6, K: source.core.K, day: 0, time: 0, pi: source.state[0] }));
  assert.ok(RAIN_FIELDS.every((name) => savedRainField(legacy, name, source).every((x) => x === 0) && savedRainField(null, name, target).every((x) => x === 0)), 'a state saved without them, or a fresh start, starts from zero');
});

test('a state saved with the convection\'s activity loads its moist fields as saved and leaves the activity behind', async () => {
  const { savedMoistField, MOIST_FIELDS } = await import('../js/physics/regrid.module.js');
  const { encodeState, decodeState } = await import('../js/stateFile.module.js');
  const { initializeState } = await import('../js/physics/init.module.js');
  assert.deepEqual(Object.keys(MOIST_FIELDS), ['convectiveRain', 'largeScaleRain']);
  const model = createModel(new Grid(6), { ocean: false }), C = model.mesh.nCells;
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  const convectiveRain = Float64Array.from({ length: C }, (_, i) => (i % 5) / 2), convectiveActivity = Float64Array.from({ length: C }, (_, i) => (i % 7) / 6);
  const saved = await decodeState(encodeState({ N: 6, K: model.core.K, day: 0, time: 0, pi: model.state[0], convectiveRain, convectiveActivity }));
  assert.ok(saved.convectiveActivity, 'the saved state carries the activity');
  for (const name of Object.keys(MOIST_FIELDS)) model.moist[name].set(savedMoistField(saved, name, model));
  for (let i = 0; i < C; i++) assert.equal(model.moist.convectiveRain[i], Math.fround(convectiveRain[i]), `convective rain of cell ${i}`);
  assert.ok(model.moist.largeScaleRain.every((x) => x === 0), 'large-scale rain saved without it starts from zero');
  assert.equal(model.moist.convectiveActivity, undefined);
  model.step(900);
  assert.ok(model.state.every((a) => a.every(Number.isFinite)), 'the model steps');
});
