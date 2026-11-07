import { test } from 'node:test';
import assert from 'node:assert/strict';
import { recordedOptions, givenOptions, headerOptions, differingGroups } from '../scripts/stateOptions.mjs';

test("a spin-up's header keeps `options`, which the page reads, to the radiation, the moist physics and the land's soil whatever other groups the run is given, and records every group in `given`", () => {
  const land = { soil: 'layers', coverGoal: 'observed', treeMoisture: false }, groups = { dynamics: { verticalAdvection: 'vanLeerWater' }, orography: false, terrain: { conserveVolume: true }, ocean: { mixing: 'kpp' } };
  const plain = headerOptions({ radiation: {}, moist: {}, land }), wide = headerOptions({ radiation: {}, moist: {}, land, ...groups });
  assert.deepEqual(wide.options, plain.options);
  assert.deepEqual(wide.options, recordedOptions({ radiation: {}, moist: {}, land: { soil: 'layers' } }));
  assert.deepEqual(Object.keys(wide.options), ['radiation', 'moist', 'land']);
  assert.deepEqual(wide.options.land, { soil: 'layers' }, 'the page sees no coverGoal');
  assert.deepEqual(Object.keys(wide.given), ['radiation', 'moist', 'land', 'dynamics', 'orography', 'terrain', 'ocean']);
  assert.deepEqual(wide.given.dynamics, { verticalAdvection: 'vanLeerWater' });
  assert.equal(wide.given.orography, false);
  assert.deepEqual(wide.given.terrain, { conserveVolume: true });
  assert.deepEqual(wide.given.land, land);
  assert.deepEqual(wide.given.radiation, wide.options.radiation);
  assert.deepEqual(Object.keys(recordedOptions(groups)), ['radiation', 'moist', 'land']);
  assert.deepEqual(givenOptions({ land: { soil: 'slab', treeMoisture: false } }).land, { soil: 'slab', treeMoisture: false });
});

test('differingGroups names the groups whose options differ, in any key order, a missing group counting', () => {
  assert.deepEqual(differingGroups({ dynamics: { a: 1, b: 2 } }, { dynamics: { b: 2, a: 1 } }), []);
  assert.deepEqual(differingGroups({ dynamics: { a: 1 }, moist: {} }, { moist: {} }), ['dynamics']);
  assert.deepEqual(differingGroups({ terrain: { conserveVolume: true } }, { terrain: { conserveVolume: false } }), ['terrain']);
});
