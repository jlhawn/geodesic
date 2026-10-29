import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { divergence } from '../js/dynamics/operators.module.js';
import { layerVerticalVelocity, layerFieldAtLevel, smoothCells } from '../js/levels.module.js';
import { memoryWeight, VERTICAL_MEMORY } from '../js/frames.module.js';
import { encodeState, decodeState } from '../js/stateFile.module.js';
import { meanVertical } from '../js/means.module.js';

/*
 * The atmosphere's memory of its layers' vertical velocity at N=16: the
 * velocity at each layer's middle is the core's own ω there; the memory
 * folds on the model's steps with the weight of the time between folds and
 * the folds it holds; a state carries it and its clock, from which a
 * model continues the mean; the physics does not see it.
 */
const N = 16, DT = 1350;

function scenario(options = {}) {
  const model = createModel(new Grid(N), options);
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  return model;
}
const velocities = (model) => layerVerticalVelocity(model.mesh, model.core, model.state[0], model.state[2], model.state[1], new Float64Array(model.core.K * model.mesh.nCells));
const same = (a, b, label) => {
  assert.equal(a.length, b.length, `${label}: lengths`);
  for (let n = 0; n < a.length; n++) if (!Object.is(a[n], b[n])) assert.fail(`${label}: ${a[n]} against ${b[n]} at ${n}`);
};
const largest = (values) => values.reduce((m, v) => Math.max(m, Math.abs(v)), 0);
const worst = (a, b) => a.reduce((m, v, n) => Math.max(m, Math.abs(v - b[n])), 0);

test('the vertical velocity of a layer is the core\'s own ω at its middle, w = −ω R T / (p g)', () => {
  const model = scenario();
  for (let n = 0; n < 3; n++) model.step(DT);
  const { mesh, core, state } = model, { K, C, E, sigmaMid, R, g } = core.diagnostics, [pi, theta, u, , q, qc] = state;
  core.diagnose(pi, theta, q, qc);
  core.phaseFlux(state, 0, K);
  const dPi = new Float64Array(C);
  core.phaseColumn(state, [dPi], 0, C);
  const divFlux = new Float64Array(core.shared.divFlux), piSigmaDot = core.diagnostics.piSigmaDot, exner = core.diagnostics.exnerLayer, plain = new Float64Array(C);
  const got = velocities(model), expected = new Float64Array(K * C);
  for (let k = 0; k < K; k++) {
    divergence(mesh, u.subarray(k * E, (k + 1) * E), plain);
    for (let i = 0; i < C; i++) {
      const omega = sigmaMid[k] * (dPi[i] + divFlux[k * C + i] - pi[i] * plain[i]) + 0.5 * (piSigmaDot[k * C + i] + piSigmaDot[(k + 1) * C + i]);
      expected[k * C + i] = -omega * R * theta[k * C + i] * exner[k * C + i] / (pi[i] * sigmaMid[k] * g);
    }
  }
  const scale = largest(expected);
  assert.ok(scale > 1e-4, `the layers rise and sink at up to ${scale} m/s`);
  assert.ok(worst(expected, got) <= 1e-9 * scale, `the layers' velocity parts from the core's ω by ${worst(expected, got)} against ${scale}`);
});

test('a layer field at a level is interpolated in ln p between the layers\' middles and held beyond them', () => {
  const K = 5, C = 3, sigmaMid = Float64Array.from([0.1, 0.3, 0.5, 0.7, 0.95]), pi = Float64Array.from([100000, 100000, 90000]);
  const field = Float64Array.from({ length: K * C }, (_, n) => 10 * Math.floor(n / C) + (n % C));
  const at = (level) => layerFieldAtLevel(sigmaMid, K, C, pi, field, level);
  assert.deepEqual(Array.from(at(500)).slice(0, 2), [20, 21], '500 hPa is the middle of layer 2 under 1000 hPa');
  assert.ok(Math.abs(at(Math.sqrt(500 * 700))[0] - 25) < 1e-6, 'halfway in ln p between layers 2 and 3');
  assert.deepEqual(Array.from(at(1)), [0, 1, 2], 'above the top layer the top layer holds');
  assert.deepEqual(Array.from(at(2000)), [40, 41, 42], 'below the lowest layer the lowest layer holds');
  const share = Math.log(900 / 700) / Math.log(950 / 700), surface = at('surface');
  assert.ok(Math.abs(surface[0] - (30 + 10 * share)) < 1e-6, 'the surface of a 1000 hPa cell is 900 hPa');
  assert.ok(Math.abs(surface[2] - (32 + 10 * share)) < 1e-6 && Math.abs(surface[2] - at(810)[2]) < 1e-9, 'and of a 900 hPa cell 810 hPa');
});

test('the memory folds on the model\'s steps with the weight of the time since the last fold and the folds it holds', () => {
  for (const every of [1, 3]) {
    const model = scenario({ verticalEvery: every });
    assert.equal(model.verticalMemory.every, every);
    assert.equal(model.verticalMemory.folds, 0);
    const reference = new Float64Array(model.verticalMemory.values.length);
    let folds = 0, last = NaN;
    for (let n = 1; n <= 9; n++) {
      model.step(DT);
      if (folds === 0 || n % every === 0) {
        const w = velocities(model), alpha = memoryWeight(folds > 0 ? model.time - last : 0, folds, VERTICAL_MEMORY);
        for (let m = 0; m < reference.length; m++) reference[m] = alpha === 1 ? w[m] : reference[m] + alpha * (w[m] - reference[m]);
        folds++; last = model.time;
      }
      assert.equal(model.verticalMemory.folds, folds, `folds after step ${n} with a fold every ${every}`);
      assert.equal(model.verticalMemory.time, last);
      const scale = largest(reference);
      assert.ok(worst(reference, model.verticalMemory.values) <= 1e-12 * Math.max(scale, 1e-12), `the memory after step ${n} with a fold every ${every}`);
    }
  }
});

test('the first fold replaces the memory, a later fold with no time between adds the plain mean, and time going back restarts it', () => {
  const model = scenario();
  const { verticalMemory: memory } = model;
  model.step(DT);
  const first = velocities(model);
  same(memory.values, first, 'one fold holds the state\'s own velocities');
  memory.fold(model.state[0], model.state[1], model.state[2], model.time);
  assert.equal(memory.folds, 2);
  assert.ok(worst(memory.values, first) <= 1e-12 * largest(first), 'the same state again at the same time leaves the mean as it was');
  memory.fold(model.state[0], model.state[1], model.state[2], 0);
  assert.equal(memory.folds, 1, 'a fold earlier than the last restarts the count');
  assert.ok(worst(memory.values, first) <= 1e-12 * largest(first));
});

test('the memory is not seen by the physics: the prognostic fields are bit for bit the same whatever its cadence', () => {
  const a = scenario({ verticalEvery: 1 }), b = scenario({ verticalEvery: 4 });
  for (let n = 0; n < 8; n++) { a.step(DT); b.step(DT); }
  for (let s = 0; s < a.state.length; s++) same(a.state[s], b.state[s], `state array ${s}`);
  assert.notEqual(a.verticalMemory.folds, b.verticalMemory.folds);
});

test('a state carries the memory and its clock, and a model continues the mean from them', async () => {
  const model = scenario();
  for (let n = 0; n < 5; n++) model.step(DT);
  const { verticalMemory: memory } = model, carried = memory.serialize();
  assert.deepEqual(Object.keys(carried).sort(), ['verticalClock', 'verticalMemory']);
  assert.deepEqual(Array.from(carried.verticalClock), [5, 5 * DT]);

  const saved = { N, K: model.core.K, pi: model.state[0], ...carried };
  const decoded = await decodeState(encodeState(saved, { f64: ['verticalMemory'] }));
  same(decoded.verticalMemory, carried.verticalMemory, 'a float64 state holds the memory exactly');
  const other = scenario();
  assert.equal(other.verticalMemory.load(decoded), true);
  assert.equal(other.verticalMemory.folds, 5);
  assert.equal(other.verticalMemory.time, 5 * DT);
  for (let n = 0; n < 4; n++) {
    model.step(DT);
    other.verticalMemory.fold(model.state[0], model.state[1], model.state[2], model.time);
  }
  same(other.verticalMemory.values, memory.values, 'the loaded mean folds as the uninterrupted one');
  assert.equal(other.verticalMemory.folds, memory.folds);

  const rounded = await decodeState(encodeState(saved));
  assert.ok(worst(carried.verticalMemory, rounded.verticalMemory) <= 1e-6 * largest(carried.verticalMemory), 'float32 states hold the mean to single precision');
  assert.equal(rounded.verticalClock[0], 5);
  assert.equal(rounded.verticalClock[1], 5 * DT, 'the clock keeps double precision');
});

test('a state without the memory, or with one of another size, starts it afresh', () => {
  const model = scenario();
  for (let n = 0; n < 3; n++) model.step(DT);
  const carried = model.verticalMemory.serialize();
  const other = scenario();
  assert.equal(other.verticalMemory.load({}), false);
  assert.equal(other.verticalMemory.folds, 0);
  assert.equal(other.verticalMemory.load({ verticalMemory: carried.verticalMemory.subarray(1), verticalClock: carried.verticalClock }), false, 'another size');
  assert.equal(other.verticalMemory.load({ verticalMemory: carried.verticalMemory }), false, 'no clock');
  assert.equal(other.verticalMemory.load(carried), true);
  assert.equal(other.verticalMemory.load(null), false);
  assert.equal(other.verticalMemory.folds, 0);
  assert.ok(other.verticalMemory.values.every((v) => v === 0));
  assert.deepEqual(other.verticalMemory.serialize(), {}, 'a memory that holds nothing is not carried');
  other.step(DT);
  assert.equal(other.verticalMemory.folds, 1, 'the next step folds');
  assert.ok(worst(other.verticalMemory.values, velocities(other)) <= 1e-12 * largest(velocities(other)), 'and replaces the memory');
});

test('the figures read a decoded state\'s memory at a level, as the overlay does', async () => {
  const model = scenario();
  for (let n = 0; n < 4; n++) model.step(DT);
  const { mesh, core, state } = model, { K } = core.diagnostics, C = mesh.nCells;
  const saved = await decodeState(encodeState({ N, K, pi: state[0], levels: core.levels, ...model.verticalMemory.serialize() }, { f64: ['verticalMemory'] }));
  for (const level of ['surface', 850, 500, 250, 50]) {
    const expected = layerFieldAtLevel(core.diagnostics.sigmaMid, K, C, state[0], model.verticalMemory.values, level), got = meanVertical(saved, level);
    assert.ok(worst(expected, got) <= 1e-6 * Math.max(largest(expected), 1e-9), `at ${level}`);
    const smoothed = meanVertical(saved, level, { mesh });
    assert.ok(worst(smoothCells(mesh, expected), smoothed) <= 1e-6 * Math.max(largest(expected), 1e-9), `smoothed at ${level}`);
  }
  assert.throws(() => meanVertical({ ...saved, verticalMemory: null }, 500), /no vertical-velocity memory/);
  assert.throws(() => meanVertical({ ...saved, verticalMemory: saved.verticalMemory.subarray(1) }, 500), /holds/);
});
