import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createOcean } from '../js/ocean/layered.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { upwellingTimescale, upwellingWeight } from '../js/frames.module.js';
import { encodeState, decodeState } from '../js/stateFile.module.js';
import { meanUpwelling } from '../js/means.module.js';
import { layerVelocities, worstDifference } from './helpers/upwelling.mjs';

/*
 * The CPU ocean's upwelling memory in a saved state at N=16 on a sea with a
 * continent and a shelf: an ocean given the time folds at each of its
 * steps, a state carries the memory with its folds and the time of the last
 * and a loaded ocean continues the mean, and the figures read a decoded
 * state at a depth as the Upwelling overlay does.
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
  const make = () => {
    const ocean = createOcean(mesh, { geography: model.geography, vertical: 'hybrid', everySteps: 1, iceDrift: false, ...options });
    ocean.initialize(Float64Array.from(surfaceT), ice);
    return ocean;
  };
  const advance = (ocean, time) => ocean.advance(Float64Array.from(surfaceT), ice, flux, stress, DT, undefined, time);
  const velocities = (ocean) => layerVelocities(mesh, ocean.layers, { h: ocean.h, u: ocean.u, eta: ocean.eta, sillEdge: ocean.sillEdge, edgeOcean: ocean.edgeOcean }, ocean.cellOcean);
  return { model, mesh, C, make, advance, velocities, surfaceT, ice };
}

// The reference mean: each layer's memory takes the velocity at its base with the weight of its own timescale.
function fold(reference, velocities, ocean, C, dt, folds) {
  const out = new Float64Array(velocities.length);
  for (let n = 0; n < out.length; n++) {
    const k = Math.floor(n / C), i = n % C;
    if (!ocean.cellOcean[i]) continue;
    let base = 0;
    for (let j = 0; j <= k; j++) base += ocean.h[j * C + i];
    const alpha = upwellingWeight(dt, folds, upwellingTimescale(base));
    out[n] = alpha === 1 ? velocities[n] : reference[n] + alpha * (velocities[n] - reference[n]);
  }
  return out;
}
const apart = (a, b) => worstDifference(a, b);

test('an ocean given the time folds at each of its steps with the weights of its layers\' timescales, and not between them', () => {
  for (const every of [1, 3]) {
    const { C, make, advance, velocities } = scenario({ everySteps: every });
    const ocean = make();
    let reference = new Float64Array(ocean.layers * C), folds = 0, last = NaN;
    for (let n = 1; n <= 7; n++) {
      const time = n * DT, stepped = advance(ocean, time);
      assert.equal(stepped, n % every === 0, `step ${n} with the ocean every ${every}`);
      if (stepped) {
        reference = fold(reference, velocities(ocean), ocean, C, folds > 0 ? time - last : 0, folds);
        folds++; last = time;
      }
      const saved = ocean.serialize();
      assert.equal(saved.upwellingFolds ?? 0, folds, `folds after step ${n} with the ocean every ${every}`);
      if (folds === 0) { assert.ok(!('upwellingMemory' in saved), 'nothing is carried before the first fold'); continue; }
      assert.equal(saved.upwellingTime, last);
      const { max, scale } = apart(reference, saved.upwellingMemory);
      assert.ok(scale > 1e-7 && max <= 1e-12, `memory after step ${n}: ${max} against ${scale}`);
    }
  }
});

test('a state carries the memory with its folds and time, a loaded ocean shows it at once and continues the mean, and one without it starts afresh', () => {
  const { C, mesh, make, advance, velocities, surfaceT, ice } = scenario();
  const source = make();
  for (let n = 1; n <= 4; n++) advance(source, n * DT);
  const saved = source.serialize();
  assert.equal(saved.upwellingFolds, 4);
  assert.equal(saved.upwellingTime, 4 * DT);

  const loaded = make();
  loaded.load(saved, Float64Array.from(surfaceT), ice);
  const again = loaded.serialize();
  assert.equal(again.upwellingFolds, 4);
  assert.equal(again.upwellingTime, 4 * DT);
  assert.equal(apart(saved.upwellingMemory, again.upwellingMemory).max, 0, 'the memory goes back as it was saved');
  for (const depth of [20, 137, 903, 2500]) {
    const { max, scale } = apart(source.fields(depth).upwelling, loaded.fields(depth).upwelling);
    assert.ok(max === 0 && scale > 1e-7, `at ${depth} m the loaded memory reads as the saved: ${max} against ${scale}`);
  }
  loaded.foldUpwelling(4 * DT + 1);
  assert.equal(loaded.serialize().upwellingFolds, 4, 'a loaded ocean folds only once it has stepped');

  advance(loaded, 5 * DT);
  const expected = fold(saved.upwellingMemory, velocities(loaded), loaded, C, DT, 4), now = loaded.serialize();
  const { max, scale } = apart(expected, now.upwellingMemory);
  assert.equal(now.upwellingFolds, 5);
  assert.ok(scale > 1e-7 && max <= 1e-12, `after a step the loaded memory continues the mean: ${max} against ${scale}`);

  const without = make(), bare = { ...saved };
  delete bare.upwellingMemory;
  without.load(bare, Float64Array.from(surfaceT), ice);
  assert.ok(!('upwellingMemory' in without.serialize()), 'a state without a memory loads with none');
  const unnamed = make();
  unnamed.load({ ...saved, upwellingTime: null }, Float64Array.from(surfaceT), ice);
  advance(unnamed, 5 * DT);
  assert.equal(unnamed.serialize().upwellingFolds, 5, 'a clock without a time still continues the mean');
  const short = make();
  short.load({ ...saved, upwellingMemory: saved.upwellingMemory.subarray(1) }, Float64Array.from(surfaceT), ice);
  assert.ok(!('upwellingMemory' in short.serialize()), 'a memory of another size is left out');
  assert.equal(mesh.nCells, C);
});

test('the figures read a decoded state\'s memory at a depth as the Upwelling overlay does', async () => {
  const { C, make, advance, surfaceT } = scenario();
  const ocean = make();
  for (let n = 1; n <= 3; n++) advance(ocean, n * DT);
  const decoded = await decodeState(encodeState({ N, surfaceT, ocean: ocean.serialize() }));
  assert.equal(decoded.ocean.upwellingFolds, 3);
  assert.equal(decoded.ocean.upwellingTime, 3 * DT);
  for (const [depth, upwellingDepth] of [[20, 20], [137, 137], [903, 903], [2500, 2500], [3900, 3900], [0, 50]]) {
    const expected = ocean.fields(depth, upwellingDepth).upwelling, got = meanUpwelling(decoded, depth, upwellingDepth);
    let seen = 0, scale = 0, max = 0;
    for (let i = 0; i < C; i++) {
      assert.equal(Number.isNaN(expected[i]), Number.isNaN(got[i]), `NaN alike at ${i}, ${depth} m`);
      if (Number.isNaN(expected[i])) continue;
      seen++; scale = Math.max(scale, Math.abs(expected[i])); max = Math.max(max, Math.abs(expected[i] - got[i]));
    }
    assert.ok(seen > 0 && max <= 1e-6 * Math.max(scale, 1e-9), `at ${depth} m: ${max} against ${scale}`);
  }
  assert.ok(meanUpwelling(decoded, 4500).every(Number.isNaN), 'no column reaches 4500 m');
  assert.throws(() => meanUpwelling({ ...decoded, ocean: { ...decoded.ocean, upwellingMemory: null } }, 100), /no upwelling memory/);
  assert.throws(() => meanUpwelling({ ...decoded, ocean: { ...decoded.ocean, upwellingMemory: decoded.ocean.upwellingMemory.subarray(1) } }, 100), /holds/);
});

test('the coupled model folds its ocean at each ocean step with the time the step ends at', () => {
  const { model } = scenario();
  const every = model.ocean.everySteps;
  assert.ok(every > 1, 'the ocean steps every few atmosphere steps');
  for (let n = 1; n <= 2 * every; n++) {
    model.step(DT);
    const saved = model.ocean.serialize(), folds = Math.floor(n / every);
    assert.equal(saved.upwellingFolds ?? 0, folds, `folds after ${n} steps`);
    if (folds) assert.equal(saved.upwellingTime, folds * every * DT);
  }
});
