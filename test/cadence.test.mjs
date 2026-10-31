import { test } from 'node:test';
import assert from 'node:assert/strict';
import { RADIATION_MINUTES, DRAG_MINUTES, OCEAN_MINUTES, stepsFor, oceanMinutes, withCadence } from '../js/cadence.module.js';

test('the cadences turn into steps of the resolutions\' dt, the radiation\'s (45 minutes with the layered shortwave, 11.25 with the sequential, from N=64 up) and the drags\' alike, the ocean\'s step from N=64 up shortening with the cell spacing, and options given in steps win', () => {
  const dt = (N) => 1350 * 16 / N;
  assert.deepEqual([128, 64, 32, 16].map((N) => stepsFor(RADIATION_MINUTES.layered, dt(N))), [16, 8, 4, 2]);
  assert.deepEqual([128, 64, 32, 16].map((N) => stepsFor(RADIATION_MINUTES.sequential, dt(N))), [4, 2, 1, 1]);
  assert.deepEqual([256, 192, 128, 64, 32, 16].map((N) => withCadence({}, dt(N), N).radiation.radiationEvery), [32, 24, 16, 8, 1, 1]);
  assert.deepEqual([256, 192, 128, 64, 32, 16].map((N) => withCadence({ radiation: { shortwave: 'sequential' } }, dt(N), N).radiation.radiationEvery), [8, 6, 4, 2, 1, 1]);
  assert.deepEqual([128, 64].map((N) => withCadence({ radiation: { shortwave: 'layered' } }, dt(N), N).radiation.radiationEvery), [16, 8]);
  assert.deepEqual([128, 64, 32, 16].map((N) => stepsFor(DRAG_MINUTES, dt(N))), [4, 2, 1, 1]);
  assert.deepEqual([128, 64, 32, 16].map((N) => withCadence({}, dt(N), N).dragEvery), [4, 2, 1, 1]);
  assert.deepEqual([128, 64].map((N) => withCadence({}, dt(N), N).ocean.everySteps), [8, 8]);
  assert.deepEqual([192, 256].map((N) => withCadence({}, dt(N), N).ocean.everySteps), [4, 2]);
  assert.deepEqual([64, 128, 192, 256].map((N) => +oceanMinutes(N).toFixed(2)), [45, 22.5, 6.67, 2.81]);
  assert.deepEqual([32, 16, 6].map((N) => withCadence({ ocean: { closureHours: 12 } }, dt(N), N).ocean), [{ closureHours: 12 }, { closureHours: 12 }, { closureHours: 12 }]);
  assert.equal(stepsFor(OCEAN_MINUTES, dt(64)), 8);
  assert.deepEqual(withCadence({}, dt(128), 128), { radiation: { radiationEvery: 16 }, ocean: { everySteps: 8 }, dragEvery: 4 });
  assert.deepEqual(withCadence({ radiation: { shortwave: 'sequential' } }, dt(128), 128), { radiation: { radiationEvery: 4, shortwave: 'sequential' }, ocean: { everySteps: 8 }, dragEvery: 4 });
  assert.deepEqual(withCadence({}, dt(128) / 2, 128).ocean, { everySteps: 16 });
  assert.deepEqual(withCadence({ radiation: { radiationEvery: 1, clearSkyPass: true }, ocean: { everySteps: 8 }, dragEvery: 1 }, dt(128), 128), { radiation: { radiationEvery: 1, clearSkyPass: true }, ocean: { everySteps: 8 }, dragEvery: 1 });
  assert.equal(withCadence({ dragEvery: 3 }, dt(64), 64).dragEvery, 3);
  assert.equal(withCadence({ ocean: false }, dt(64), 64).ocean, false);
});
