import { test } from 'node:test';
import assert from 'node:assert/strict';
import { DEFAULT_RUNS, DEFAULT_RUN, defaultRunFor } from '../js/defaultRun.module.js';

test('a resolution starts from the lowest saved run at least as fine, and above the finest from the finest', () => {
  assert.equal(DEFAULT_RUN, DEFAULT_RUNS[128]);
  for (const N of [8, 16, 32, 48, 64]) assert.equal(defaultRunFor(N), DEFAULT_RUNS[64], `N=${N}`);
  for (const N of [65, 96, 128]) assert.equal(defaultRunFor(N), DEFAULT_RUNS[128], `N=${N}`);
  for (const N of [129, 160, 192, 256]) assert.equal(defaultRunFor(N), DEFAULT_RUNS[192], `N=${N}`);
  assert.equal(defaultRunFor(48, { 32: 'a', 64: 'b', 128: 'c' }), 'b');
  assert.equal(defaultRunFor(32, { 32: 'a', 64: 'b', 128: 'c' }), 'a');
  assert.equal(defaultRunFor(300, { 32: 'a', 64: 'b', 128: 'c' }), 'c');
});
