import { test } from 'node:test';
import assert from 'node:assert/strict';
import { RADIATION_MINUTES, DRAG_MINUTES, OCEAN_MINUTES, stepsFor, withCadence } from '../js/cadence.module.js';

test('the cadences turn into steps of the resolutions\' dt, the radiation\'s and the drags\' alike, the ocean\'s step from N=64 up shortening with the cell spacing, and options given in steps win', () => {
  const dt = (N) => 1350 * 16 / N;
  assert.deepEqual([128, 64, 32, 16].map((N) => stepsFor(RADIATION_MINUTES, dt(N))), [4, 2, 1, 1]);
  assert.deepEqual([128, 64, 32, 16].map((N) => stepsFor(DRAG_MINUTES, dt(N))), [4, 2, 1, 1]);
  assert.deepEqual([128, 64, 32, 16].map((N) => withCadence({}, dt(N), N).dragEvery), [4, 2, 1, 1]);
  assert.deepEqual([128, 64].map((N) => withCadence({}, dt(N), N).ocean.everySteps), [8, 8]);
  assert.deepEqual([32, 16, 6].map((N) => withCadence({ ocean: { closureHours: 12 } }, dt(N), N).ocean), [{ closureHours: 12 }, { closureHours: 12 }, { closureHours: 12 }]);
  assert.equal(stepsFor(OCEAN_MINUTES, dt(64)), 8);
  assert.deepEqual(withCadence({}, dt(128), 128), { radiation: { radiationEvery: 4 }, ocean: { everySteps: 8 }, dragEvery: 4 });
  assert.deepEqual(withCadence({}, dt(128) / 2, 128).ocean, { everySteps: 16 });
  assert.deepEqual(withCadence({ radiation: { radiationEvery: 1, clearSkyPass: true }, ocean: { everySteps: 8 }, dragEvery: 1 }, dt(128), 128), { radiation: { radiationEvery: 1, clearSkyPass: true }, ocean: { everySteps: 8 }, dragEvery: 1 });
  assert.equal(withCadence({ dragEvery: 3 }, dt(64), 64).dragEvery, 3);
  assert.equal(withCadence({ ocean: false }, dt(64), 64).ocean, false);
});
