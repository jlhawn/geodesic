import { test } from 'node:test';
import assert from 'node:assert/strict';
import { createSlice } from '../explain/sliceCore.module.js';

const finite = (array) => array.every(Number.isFinite);
const maxAbs = (array) => array.reduce((m, v) => Math.max(m, Math.abs(v)), 0);

test('a resting stable atmosphere stays at rest', () => {
  const slice = createSlice();
  const theta = Float64Array.from(slice.theta);
  for (let i = 0; i < 120; i++) slice.step();
  assert.equal(maxAbs(slice.u), 0);
  for (let i = 0; i < theta.length; i++) assert.ok(Math.abs(slice.theta[i] - theta[i]) < 1e-9);
  assert.ok(slice.pi.every((p) => p === 1e5));
});

test('heating a block of columns conserves mass and drives a thermally direct cell', () => {
  const slice = createSlice();
  const { M, K, heated: c, heatedHalf: h } = slice, mass = slice.mass(), right = c + h, left = c - h - 1, outside = c + h + 1;
  for (let i = 0; i < 60; i++) slice.step({ heating: 1 / 3600 });
  assert.ok(Math.abs(slice.mass() - mass) / mass < 1e-12);
  assert.ok(finite(slice.u) && finite(slice.theta) && finite(slice.pi));
  assert.ok(slice.u[1 * M + right] > 0 && slice.u[1 * M + left] < 0, 'outflow aloft at both edges of the warmed block');
  assert.ok(slice.u[(K - 1) * M + right] < 0 && slice.u[(K - 1) * M + left] > 0, 'inflow at the ground at both edges');
  assert.ok(slice.pi[c] < slice.pi[outside] && slice.pi[c] < slice.pi[c - h - 1], 'surface pressure falls under the warmed block');
  assert.ok(slice.verticalVelocity(K >> 1, c) > 0 && slice.verticalVelocity(K >> 1, outside) < 0, 'rising in the middle, sinking next door');
  assert.ok(slice.temperature(K >> 1, outside) > slice.initialTemperature(K >> 1), 'the neighbor warms aloft by sinking');
  assert.ok(maxAbs(slice.u) < 10);
});

test('ground heating with convective adjustment and friction stays bounded', () => {
  const slice = createSlice();
  const { M, K, heated: c } = slice;
  for (let i = 0; i < 360; i++) slice.step({ heating: 3 / 3600, mode: 'ground', friction: true });
  assert.ok(finite(slice.u) && finite(slice.theta));
  for (let k = 0; k < K - 1; k++) assert.ok(slice.theta[(k + 1) * M + c] <= slice.theta[k * M + c] + 1e-2, `layer ${k} is not unstable`);
  assert.ok(maxAbs(slice.u) < 30);
});
