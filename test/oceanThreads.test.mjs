import { test } from 'node:test';
import assert from 'node:assert/strict';
import { ekmanCase, lensCase, bumpCase, basinCase } from '../explain/oceanCases.module.js';
import { threaded } from './helpers/oceanThreads.mjs';

const CASES = {
  'the walled basin': () => basinCase(),
  'the basin with no walls': () => basinCase({ walls: false }),
  'the walled basin under one f': () => basinCase({ beta: false, ocean: { closureHours: 48 } }),
  'the walled basin with mixed-layer restratification': () => basinCase({ ocean: { mixedLayerRestratification: true } }),
  'the Ekman column': () => ekmanCase({ N: 16 }),
  'the planet under the trades': () => ekmanCase({ N: 16, view: 'planet' }),
  'the lens': () => lensCase(),
  'the bump': () => bumpCase(),
};

test('every oceanCases kind steps bit for bit on three threads as on one', async () => {
  for (const [name, build] of Object.entries(CASES)) {
    const serial = build(), sim = await threaded(build(), 3);
    try {
      for (let n = 0; n < 4; n++) { serial.step(); sim.step(); }
    } finally {
      await sim.close();
    }
    for (const [key, buffer] of Object.entries(serial.ocean.shared)) assert.ok(Buffer.from(sim.ocean.shared[key]).equals(Buffer.from(buffer)), `${name}: ${key}`);
    assert.deepStrictEqual(sim.ocean.T0, serial.ocean.T0, `${name}: T0`);
    assert.deepStrictEqual(sim.ocean.S0, serial.ocean.S0, `${name}: S0`);
    assert.deepStrictEqual(sim.frame(), serial.frame(), `${name}: frame`);
  }
});
