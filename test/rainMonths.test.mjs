import { test } from 'node:test';
import assert from 'node:assert/strict';
import { quarterSaves, saveDay, monthWeights, weightLabel, windowDates, areaStats } from '../scripts/figures/rainMonths.mjs';

const YEAR_ONE = [
  [30, 'Mar 20 – Apr 19', { 2: 0.383, 3: 0.617 }], [61, 'Apr 20 – May 20', { 3: 0.35, 4: 0.65 }], [91, 'May 20 – Jun 19', { 4: 0.383, 5: 0.617 }],
  [122, 'Jun 20 – Jul 20', { 5: 0.35, 6: 0.65 }], [152, 'Jul 20 – Aug 19', { 6: 0.383, 7: 0.617 }], [183, 'Aug 20 – Sep 19', { 7: 0.383, 8: 0.617 }],
  [213, 'Sep 19 – Oct 19', { 8: 0.383, 9: 0.617 }], [243, 'Oct 19 – Nov 18', { 9: 0.417, 10: 0.583 }], [274, 'Nov 19 – Dec 19', { 10: 0.383, 11: 0.617 }],
  [304, 'Dec 19 – Jan 18', { 11: 0.417, 0: 0.583 }], [335, 'Jan 19 – Feb 18', { 0: 0.417, 1: 0.583 }], [365, 'Feb 18 – Mar 20', { 1: 0.35, 2: 0.65 }],
];

test('year 1\x27s twelve save windows weigh two calendar months each, summing to 1', () => {
  YEAR_ONE.forEach(([day, dates, expected], m) => {
    assert.equal(saveDay(m + 1), day);
    const w = monthWeights(day);
    assert.ok(Math.abs(w.reduce((a, b) => a + b, 0) - 1) < 1e-12, `day ${day} sums to ${w.reduce((a, b) => a + b, 0)}`);
    w.forEach((x, month) => assert.ok(Math.abs(x - (expected[month] ?? 0)) < 1e-3, `day ${day} month ${month}: ${x}`));
    assert.equal(windowDates(day), dates);
  });
  assert.equal(weightLabel(monthWeights(30)), 'Mar 0.38, Apr 0.62');
});

test('the windows of the next year repeat the first year\x27s weights', () => {
  YEAR_ONE.forEach(([day]) => monthWeights(day + 365).forEach((x, m) => assert.ok(Math.abs(x - monthWeights(day)[m]) < 1e-12)));
});

test('a quarter report holds the three saves ending at its own day', () => {
  assert.deepEqual([1, 2, 3, 4, 5].map((k) => quarterSaves(k)), [[30, 61, 91], [122, 152, 183], [213, 243, 274], [304, 335, 365], [395, 426, 456]]);
  assert.deepEqual(quarterSaves(1, 2), [30, 61, 91, 122, 152, 183]);
  assert.deepEqual(quarterSaves(67).at(-1), 6114);
  assert.throws(() => quarterSaves(1, 5));
  assert.throws(() => quarterSaves(0));
});

test('statistics of constant fields: area-weighted means, their difference and the mean absolute difference', () => {
  const area = Float64Array.from([1, 2, 3, 4]), land = Uint8Array.from([1, 0, 0, 1]);
  const flat = areaStats(area, land, new Float64Array(4).fill(90), new Float64Array(4).fill(60));
  for (const part of ['global', 'land', 'sea']) assert.deepEqual(flat[part], { model: 90, obs: 60, diff: 30, absDiff: 30 });
  const split = areaStats(area, land, Float64Array.from([10, 20, 30, 40]), new Float64Array(4).fill(25));
  assert.equal(split.global.model, (10 + 40 + 90 + 160) / 10);
  assert.equal(split.land.model, (10 + 160) / 5);
  assert.equal(split.sea.diff, (2 * -5 + 3 * 5) / 5 );
  assert.equal(split.global.absDiff, (15 + 2 * 5 + 3 * 5 + 4 * 15) / 10);
});
