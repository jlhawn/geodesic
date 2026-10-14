import { test } from 'node:test';
import assert from 'node:assert/strict';
import { createNightCensus } from '../scripts/nightCensus.mjs';

const DT = 900;
// Runs one cell through the local hours from `from` to `to`, calling observe at each step's mid hour
function night(census, i, { from = 22, to = 3.75, insolation = () => 0, flux = () => -50, warm = () => 0 }) {
  for (let h = from + DT / 7200; ; h += DT / 3600) {
    const hour = ((h % 24) + 24) % 24;
    census.observe(i, hour, insolation(hour), flux(hour), DT, warm(hour));
    if (h >= to + (from > to ? 24 : 0)) break;
  }
}

test('a dark cell losing heat counts, and its warm layer is binned at 0.1 K and 1 K', () => {
  const c = createNightCensus(3);
  night(c, 0, { warm: () => 0.05 });
  night(c, 1, { warm: () => 0.4 });
  night(c, 2, { warm: () => 2.5 });
  const t = c.take();
  assert.ok(t.dark > 0 && t.darkWarm === t.dark);
  assert.equal(t.darkOver1 * 3, t.dark, 'one cell in three is over 1 K at every sample');
  assert.equal(t.darkOver01 * 3, 2 * t.dark);
  assert.equal(t.darkWarmest, 2.5);
  assert.equal(t.darkCell, 2);
});

test('polar day (the sun up at some hour of 00 to 03) and a net gain (warm air over cold water) are not counted', () => {
  const c = createNightCensus(2);
  night(c, 0, { insolation: (h) => (h < 1 ? 0 : 5), warm: () => 0.6 });
  night(c, 1, { flux: () => 30, warm: () => 0.6 });
  const t = c.take();
  assert.equal(t.dark, 0);
  assert.equal(t.darkOver01, 0);
});

test('the window is the hours since 00: a gain before it does not count and a loss after midnight does', () => {
  const c = createNightCensus(1);
  night(c, 0, { flux: (h) => (h > 12 ? 400 : -20), warm: () => 0.3 });
  const t = c.take();
  assert.ok(t.dark > 0 && t.darkOver01 === t.dark, 'the evening\'s gain is outside the window');
});

test('the counts restart after take and the accumulators restart each night', () => {
  const c = createNightCensus(1);
  night(c, 0, { warm: () => 0.3 });
  assert.ok(c.take().dark > 0);
  assert.equal(c.take().dark, 0);
  night(c, 0, { flux: () => 30, warm: () => 0.3 });
  assert.equal(c.take().dark, 0, 'a night of gain after a night of loss starts from zero');
});

test('hours outside 2.5 to 3.5 are not sampled', () => {
  const c = createNightCensus(1);
  night(c, 0, { from: 22, to: 2.2, warm: () => 3 });
  assert.equal(c.take().dark, 0);
});
