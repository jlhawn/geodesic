import { test } from 'node:test';
import assert from 'node:assert/strict';

globalThis.IntersectionObserver ??= class { observe() {} };
const { amplification } = await import('../explain/figures/stepper.module.js');

test('forward Euler grows every oscillation and RK4 keeps it up to 2√2', () => {
  for (const x of [0.1, 0.5, 1, 2.6]) assert.ok(Math.abs(amplification('euler', x) - Math.hypot(1, x)) < 1e-12);
  for (const x of [0.5, 1, 2, 2.6, 3]) assert.ok(Math.abs(amplification('rk4', x) ** 2 - (1 - x ** 6 / 72 + x ** 8 / 576)) < 1e-12);
  assert.ok(Math.abs(amplification('rk4', 2 * Math.SQRT2) - 1) < 1e-12);
});
