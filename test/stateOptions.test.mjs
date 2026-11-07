import { test } from 'node:test';
import assert from 'node:assert/strict';
import { recordedOptions, differingGroups } from '../scripts/stateOptions.mjs';

test('a header records every group a run was given over the defaults it records, and false as false', () => {
  assert.deepEqual(Object.keys(recordedOptions()), ['radiation', 'moist', 'land']);
  const recorded = recordedOptions({ radiation: {}, moist: {}, land: {}, dynamics: { verticalAdvection: 'vanLeerWater' }, orography: false, terrain: { conserveVolume: true }, ocean: { mixing: 'kpp' } });
  assert.deepEqual(Object.keys(recorded), ['radiation', 'moist', 'land', 'dynamics', 'orography', 'terrain', 'ocean']);
  assert.deepEqual(recorded.dynamics, { verticalAdvection: 'vanLeerWater' });
  assert.equal(recorded.orography, false);
  assert.deepEqual(recorded.terrain, { conserveVolume: true });
  assert.deepEqual(recordedOptions({ land: { soil: 'slab', treeMoisture: false } }).land, { soil: 'slab', treeMoisture: false });
});

test('differingGroups names the groups whose options differ, in any key order, a missing group counting', () => {
  assert.deepEqual(differingGroups({ dynamics: { a: 1, b: 2 } }, { dynamics: { b: 2, a: 1 } }), []);
  assert.deepEqual(differingGroups({ dynamics: { a: 1 }, moist: {} }, { moist: {} }), ['dynamics']);
  assert.deepEqual(differingGroups({ terrain: { conserveVolume: true } }, { terrain: { conserveVolume: false } }), ['terrain']);
});
