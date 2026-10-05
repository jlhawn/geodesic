import { test } from 'node:test';
import assert from 'node:assert/strict';
import { visibleCloudHeights, CLOUD_OPACITY_PATH } from '../js/frames.module.js';

// A column of K layers 1 km apart, the top one's middle at (K - 0.5) km, over sea-level ground.
const K = 10;
const z = Float64Array.from({ length: K }, (_, k) => 1000 * (K - k - 0.5));
const upper = (k) => z[k] + 500, lower = (k) => z[k] - 500;
const heights = (path, deck = 0, deckTop = 0) => visibleCloudHeights(K, z, 0, Float64Array.from(path), deck, deckTop, [0, 0]);
const column = (entries) => { const path = new Array(K).fill(0); for (const [k, p] of Object.entries(entries)) path[k] = p; return path; };

test('a single cloudy layer is seen at its own bounds', () => {
  for (const p of [0.002, 0.04, 1]) {
    const [top, base] = heights(column({ 6: p }));
    assert.ok(Math.abs(top - upper(6)) < 1e-9 && Math.abs(base - lower(6)) < 1e-9, `${p} kg/m²: ${top}, ${base}`);
  }
  const [, base] = heights(column({ 9: 0.1 }));
  assert.equal(base, 0, 'the lowest layer reaches down to the ground');
});

test('an opaque layer hides what lies beyond it', () => {
  const [top] = heights(column({ 2: 2, 7: 2 }));
  assert.ok(Math.abs(top - upper(2)) < 1e-6, `top ${top}`);
  const [, base] = heights(column({ 2: 2, 7: 2 }));
  assert.ok(Math.abs(base - lower(7)) < 1e-6, `base ${base}`);
});

test('thin high cloud over thick low cloud puts the top between them', () => {
  const [top, base] = heights(column({ 1: 0.2 * CLOUD_OPACITY_PATH, 8: 3 }));
  assert.ok(top > upper(8) + 100 && top < upper(1) - 100, `top ${top}`);
  assert.ok(Math.abs(base - lower(8)) < 1, `base ${base}`);
});

test('a clear column, or one less than 1% opaque, has no heights', () => {
  assert.deepEqual(heights(column({})), [0, 0]);
  assert.deepEqual(heights(column({ 4: 0.009 * CLOUD_OPACITY_PATH })), [0, 0]);
});

test('the deck lands in the lowest layer reaching its top', () => {
  const [top, base] = heights(column({}), 0.1, 2200);
  assert.ok(Math.abs(top - upper(7)) < 1e-9 && Math.abs(base - lower(7)) < 1e-9, `${top}, ${base}`);
  const [high] = heights(column({}), 0.1, 50000);
  assert.ok(Math.abs(high - upper(0)) < 1e-9, 'a deck above the column joins the top layer');
});

test('cloud added above never lowers the top, and the base never rises above it', () => {
  let seed = 3;
  const random = () => ((seed = (seed * 16807) % 2147483647) / 2147483647);
  for (let n = 0; n < 300; n++) {
    const path = Array.from({ length: K }, () => (random() < 0.4 ? 0.1 * random() ** 2 : 0));
    const deck = random() < 0.3 ? 0.05 * random() : 0, deckTop = 3000 * random();
    const [top, base] = heights(path, deck, deckTop);
    if (top > 0) assert.ok(base <= top + 1e-9, `base ${base} above top ${top}`);
    let highest = path.findIndex((p) => p > 0);
    if (deck > 0) { let d = K - 1; while (d > 0 && upper(d) < deckTop) d--; highest = highest < 0 ? d : Math.min(highest, d); }
    if (highest < 0) continue;
    const k = Math.floor(random() * (highest + 1));
    const more = path.slice();
    more[k] += 0.02 * random();
    const [raised] = heights(more, deck, deckTop);
    if (top > 0) assert.ok(raised >= top - 1e-6, `adding cloud to layer ${k}, at or above the highest cloud in ${highest}, lowered the top from ${top} to ${raised}`);
  }
});
