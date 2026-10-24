import { test } from 'node:test';
import assert from 'node:assert/strict';
import { foldMemory, UPWELLING_MEMORY } from '../js/frames.module.js';

const OMEGA = 7.2921e-5, INERTIAL_30 = 2 * OMEGA * Math.sin(Math.PI / 6);

function response(frameSeconds, omega, days = 40) {
  const memory = { values: null, key: null, time: 0 };
  let low = Infinity, high = -Infinity, mean = 0, n = 0;
  for (let t = 0; t <= days * 86400; t += frameSeconds) {
    const [value] = foldMemory(memory, [1 + 50 * Math.sin(omega * t)], t, 100, UPWELLING_MEMORY);
    if (t >= (days - 3) * 86400) { low = Math.min(low, value); high = Math.max(high, value); mean += value; n++; }
  }
  return { steady: mean / n, oscillation: (high - low) / 2 / 50 };
}

test('the upwelling memory keeps a steady signal and leaves 3 percent of a near-inertial oscillation at 30 degrees', () => {
  const expected = 1 / Math.hypot(1, INERTIAL_30 * UPWELLING_MEMORY);
  assert.ok(Math.abs(expected - 0.0317) < 0.001);
  for (const frameSeconds of [506, 2025, 8100]) {
    const { steady, oscillation } = response(frameSeconds, INERTIAL_30);
    assert.ok(Math.abs(steady - 1) < 0.05, `${frameSeconds} s frames: steady ${steady}`);
    assert.ok(Math.abs(oscillation - expected) < 0.1 * expected, `${frameSeconds} s frames: oscillation ${oscillation} vs ${expected}`);
  }
});

test('a fresh upwelling memory is the running mean, within 2/(ωT) of a steady signal from day 1 to day 5', () => {
  for (const frameSeconds of [506, 2025]) {
    for (let n = 0; n < 8; n++) {
      const phase = n * Math.PI / 4, memory = { values: null, key: null, time: 0 };
      for (let t = 0; t <= 5 * 86400; t += frameSeconds) {
        const [value] = foldMemory(memory, [1 + 50 * Math.sin(INERTIAL_30 * t + phase)], t, 100, UPWELLING_MEMORY);
        if (t < 86400) continue;
        const error = Math.abs(value - 1) / 50, bound = 2 / (INERTIAL_30 * t);
        assert.ok(error <= 1.02 * bound, `${frameSeconds} s frames, phase ${phase.toFixed(2)}, day ${(t / 86400).toFixed(2)}: ${error} vs ${bound}`);
      }
    }
  }
  assert.ok(Math.abs(2 / (INERTIAL_30 * 2 * 86400) - 0.16) < 0.01 && Math.abs(2 / (INERTIAL_30 * 5 * 86400) - 0.063) < 0.005);
});

test('the upwelling memory restarts at a new depth, when time goes back, and in a cell that was not finite', () => {
  const memory = { values: null, key: null, time: 0 };
  assert.deepEqual([...foldMemory(memory, [1, 2], 0, 100, UPWELLING_MEMORY)], [1, 2]);
  const half = foldMemory(memory, [3, NaN], UPWELLING_MEMORY * Math.LN2, 100, UPWELLING_MEMORY);
  assert.ok(Math.abs(half[0] - 2) < 1e-6);
  assert.ok(Number.isNaN(half[1]));
  const wet = foldMemory(memory, [2, 5], UPWELLING_MEMORY * Math.LN2 + 600, 100, UPWELLING_MEMORY);
  assert.equal(wet[1], 5);
  assert.deepEqual([...foldMemory(memory, [7, 8], memory.time + 600, 'surface', UPWELLING_MEMORY)], [7, 8]);
  assert.deepEqual([...foldMemory(memory, [4, 4], memory.time - 600, 'surface', UPWELLING_MEMORY)], [4, 4]);
  const same = foldMemory(memory, [9, 9], memory.time, 'surface', UPWELLING_MEMORY);
  assert.deepEqual([...same], [4, 4]);
  same[0] = 0;
  assert.equal(memory.values[0], 4);
});
