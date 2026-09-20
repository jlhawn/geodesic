import test from 'node:test';
import assert from 'node:assert/strict';
import { createDisplayClock } from '../js/displayClock.module.js';

const STEP = 2025, INTERVAL = 450, RATE = STEP / INTERVAL, FRAME = 1000 / 120;

function run(clock, { seconds, stall = null, jitter = 50 }) {
  let target = 0, next = 0, seed = 1, last = null;
  const random = () => (seed = (seed * 48271) % 2147483647) / 2147483647;
  const speeds = [], lags = [];
  let maxAhead = -Infinity;
  for (let now = 0; now < 1000 * seconds; now += FRAME) {
    const stalled = stall && now > stall[0] && now < stall[1];
    if (!stalled && now >= next) { target += STEP; next = Math.max(next, now) + INTERVAL + (random() - 0.5) * 2 * jitter; }
    const time = clock.advance(now, target, { rate: RATE, interval: INTERVAL, running: true });
    if (last !== null && now > 20000) { speeds.push((time - last) / FRAME / RATE); lags.push((target - time) / STEP); }
    maxAhead = Math.max(maxAhead, (time - target) / STEP);
    last = time;
  }
  return { speeds, lags, maxAhead };
}

test('the display clock follows jittery model frames at a nearly constant speed, a couple of frames behind', () => {
  const { speeds, lags } = run(createDisplayClock(), { seconds: 90 });
  const mean = speeds.reduce((a, b) => a + b, 0) / speeds.length;
  assert.ok(Math.abs(mean - 1) < 0.02, `mean speed ${mean.toFixed(3)} of nominal`);
  assert.ok(Math.max(...speeds) < 1.08 && Math.min(...speeds) > 0.92, `speed ripple ${Math.min(...speeds).toFixed(3)}–${Math.max(...speeds).toFixed(3)}`);
  assert.ok(Math.min(...lags) > 0, 'never runs ahead of the model while frames keep coming');
  assert.ok(Math.max(...lags) < 4, `lag stays under four frames (max ${Math.max(...lags).toFixed(2)})`);
});

test('a stalled model leaves the display clock at most a few frames ahead, and a jump snaps it', () => {
  const { maxAhead } = run(createDisplayClock(), { seconds: 40, stall: [15000, 18000] });
  assert.ok(maxAhead > 0 && maxAhead < 3, `ran ${maxAhead.toFixed(2)} frames ahead during the stall`);
  const clock = createDisplayClock();
  clock.advance(0, 0, { rate: RATE, interval: INTERVAL });
  clock.advance(8, 100 * STEP, { rate: RATE, interval: INTERVAL });
  assert.equal(clock.time, 100 * STEP);
});

test('while paused the display clock eases forward onto the model time and never runs backwards', () => {
  const clock = createDisplayClock();
  clock.advance(0, 0, { rate: 0, interval: INTERVAL, running: false });
  clock.advance(100, STEP, { rate: 0, interval: INTERVAL, running: false });
  assert.ok(clock.time > 0 && clock.time < STEP);
  clock.advance(2000, STEP, { rate: 0, interval: INTERVAL, running: false });
  assert.equal(clock.time, STEP);
  clock.advance(2100, STEP / 2, { rate: 0, interval: INTERVAL, running: false });
  assert.equal(clock.time, STEP);
});
