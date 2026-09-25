import test from 'node:test';
import assert from 'node:assert/strict';
import { createPacer, PAUSE_LEVELS } from '../js/pace.module.js';

/*
 * Runs the pacer against a page whose late frames each second depend on
 * the pause through `lateAt(pause, second)`, and returns the pauses.
 */
function run(lateAt, seconds = 600) {
  const pacer = createPacer(), pauses = [];
  let pause = 0;
  for (let s = 0; s < seconds; s++) { pause = pacer.report(lateAt(pause, s)); pauses.push(pause); }
  return pauses;
}
const share = (pauses, predicate, from = 0) => pauses.slice(from).filter(predicate).length / (pauses.length - from);

test('no late frames, no pause', () => {
  assert.ok(run(() => 0).every((p) => p === 0));
});

test('late frames the pause does not cure leave it off nearly all the time', () => {
  const pauses = run(() => 3);
  assert.ok(share(pauses, (p) => p === 0) > 0.9, `paused ${(100 * share(pauses, (p) => p > 0)).toFixed(0)}% of the time`);
  assert.ok(Math.max(...pauses) <= PAUSE_LEVELS[1]);
});

test('late frames from GPU contention settle on the pause that cures them', () => {
  const pauses = run((pause) => (pause >= 4 ? 0 : pause >= 2 ? 3 : 6));
  assert.ok(share(pauses, (p) => p === 4, 60) > 0.8, `at 4 ms ${(100 * share(pauses, (p) => p === 4, 60)).toFixed(0)}% of the time`);
  assert.ok(share(pauses, (p) => p < 4, 60) < 0.2);
});

test('the pause lets go once the contention ends', () => {
  const pauses = run((pause, s) => (s < 200 && pause < 2 ? 5 : 0));
  assert.ok(pauses.slice(20, 200).some((p) => p >= 2));
  assert.equal(pauses.slice(400).filter((p) => p > 0).length, 0);
});

test('one late second now and then does not start a pause', () => {
  assert.ok(run((pause, s) => (s % 7 === 0 ? 2 : s % 3 === 0 ? 1 : 0)).every((p) => p === 0));
});

test('a late frame every second after the contention ends still lets the pause go', () => {
  const pauses = run((pause, s) => (s < 60 ? (pause >= 2 ? 0 : 5) : 1), 400);
  assert.ok(pauses.slice(10, 60).some((p) => p >= 2));
  assert.equal(pauses.slice(200).filter((p) => p > 0).length, 0);
});
