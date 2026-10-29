import { test } from 'node:test';
import assert from 'node:assert/strict';
import { upwellingTimescale, upwellingWeight, UPWELLING_FAST, UPWELLING_SLOW, UPWELLING_FAST_DEPTH, UPWELLING_SLOW_DEPTH } from '../js/frames.module.js';

const OMEGA = 7.2921e-5, INERTIAL_30 = 2 * OMEGA * Math.sin(Math.PI / 6), DAY = 86400;

// One layer's memory at `depth`, folded every frameSeconds as the ocean folds it.
function memoryOf(depth, frameSeconds) {
  const timescale = upwellingTimescale(depth);
  const memory = { value: 0, folds: 0 };
  memory.fold = (sample) => {
    memory.value = memory.folds === 0 ? sample : memory.value + upwellingWeight(frameSeconds, memory.folds, timescale) * (sample - memory.value);
    memory.folds++;
    return memory.value;
  };
  return memory;
}

function response(depth, frameSeconds, omega, days) {
  const memory = memoryOf(depth, frameSeconds);
  let low = Infinity, high = -Infinity, mean = 0, n = 0;
  for (let t = 0; t <= days * DAY; t += frameSeconds) {
    const value = memory.fold(1 + 50 * Math.sin(omega * t));
    if (t >= (days - 3) * DAY) { low = Math.min(low, value); high = Math.max(high, value); mean += value; n++; }
  }
  return { steady: mean / n, oscillation: (high - low) / 2 / 50 };
}

test('the timescale is 5 days down to 50 m and 30 days from 150 m, linear in depth between', () => {
  assert.equal(UPWELLING_FAST, 5 * DAY);
  assert.equal(UPWELLING_SLOW, 30 * DAY);
  assert.equal(UPWELLING_FAST_DEPTH, 50);
  assert.equal(UPWELLING_SLOW_DEPTH, 150);
  for (const depth of [0, 10, 50]) assert.equal(upwellingTimescale(depth), 5 * DAY, `${depth} m`);
  for (const depth of [150, 300, 500, 5500]) assert.equal(upwellingTimescale(depth), 30 * DAY, `${depth} m`);
  assert.ok(Math.abs(upwellingTimescale(100) - 17.5 * DAY) < 1e-6);
  assert.ok(Math.abs(upwellingTimescale(75) - 11.25 * DAY) < 1e-6);
  let before = 0;
  for (let depth = 0; depth <= 200; depth += 5) { const t = upwellingTimescale(depth); assert.ok(t >= before, `the timescale does not fall with depth at ${depth} m`); before = t; }
});

test('the weight is the larger of the exponential\'s and the running mean\'s, and 1 for an empty memory', () => {
  const tau = 5 * DAY;
  for (const dt of [1, 506, 8100, 5 * DAY, 1e9]) assert.equal(upwellingWeight(dt, 0, tau), 1, `an empty memory takes the whole frame at dt ${dt}`);
  assert.ok(Math.abs(upwellingWeight(506, 1, tau) - 0.5) < 1e-12);
  assert.ok(Math.abs(upwellingWeight(506, 99, tau) - 0.01) < 1e-12);
  assert.ok(Math.abs(upwellingWeight(tau, 10, tau) - (1 - Math.exp(-1))) < 1e-12);
  assert.ok(Math.abs(upwellingWeight(1e9, 1000, tau) - 1) < 1e-12);
  let folds = 0;
  while (upwellingWeight(2025, folds, tau) === 1 / (folds + 1)) folds++;
  assert.ok(Math.abs(folds * 2025 - tau) < 2 * 2025, `the exponential takes over after ${folds} frames, ${(folds * 2025 / DAY).toFixed(2)} days`);
});

test('the memory keeps a steady signal and leaves a few percent of a near-inertial oscillation at 30 degrees: per 5 days at 50 m, per 30 days at 500 m', () => {
  for (const [depth, days] of [[50, 40], [500, 200]]) {
    const expected = 1 / Math.hypot(1, INERTIAL_30 * upwellingTimescale(depth));
    for (const frameSeconds of [506, 2025, 8100]) {
      const { steady, oscillation } = response(depth, frameSeconds, INERTIAL_30, days);
      assert.ok(Math.abs(steady - 1) < 0.05, `${depth} m, ${frameSeconds} s frames: steady ${steady}`);
      assert.ok(Math.abs(oscillation - expected) < 0.1 * expected, `${depth} m, ${frameSeconds} s frames: oscillation ${oscillation} vs ${expected}`);
    }
  }
  assert.ok(Math.abs(1 / Math.hypot(1, INERTIAL_30 * upwellingTimescale(50)) - 0.0317) < 0.001);
  assert.ok(Math.abs(1 / Math.hypot(1, INERTIAL_30 * upwellingTimescale(500)) - 0.00529) < 0.0002);
  const steady = memoryOf(500, 2025);
  for (let n = 0; n < 2000; n++) steady.fold(3.5);
  assert.equal(steady.value, 3.5);
});

test('a fresh memory is the running mean of its frames, within 2/(ωT) of a steady signal from day 1 to day 5, and hands over to the exponential', () => {
  const samples = Array.from({ length: 200 }, (_, n) => Math.sin(1.7 * n) + 0.01 * n);
  for (const depth of [50, 500]) {
    const memory = memoryOf(depth, 2025);
    let sum = 0;
    samples.forEach((sample, n) => {
      sum += sample;
      assert.ok(Math.abs(memory.fold(sample) - sum / (n + 1)) < 1e-12, `${depth} m: fold ${n + 1} is the mean`);
    });
  }
  for (const depth of [50, 500]) {
    for (const frameSeconds of [506, 2025]) {
      for (let n = 0; n < 8; n++) {
        const phase = n * Math.PI / 4, memory = memoryOf(depth, frameSeconds);
        for (let t = 0; t <= 5 * DAY; t += frameSeconds) {
          const value = memory.fold(1 + 50 * Math.sin(INERTIAL_30 * t + phase));
          if (t < DAY) continue;
          const error = Math.abs(value - 1) / 50, bound = 2 / (INERTIAL_30 * t);
          assert.ok(error <= 1.02 * bound, `${depth} m, ${frameSeconds} s frames, phase ${phase.toFixed(2)}, day ${(t / DAY).toFixed(2)}: ${error} vs ${bound}`);
        }
      }
    }
  }
  assert.ok(Math.abs(2 / (INERTIAL_30 * 2 * DAY) - 0.16) < 0.01 && Math.abs(2 / (INERTIAL_30 * 5 * DAY) - 0.063) < 0.005);
});
