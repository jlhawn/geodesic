import { test } from 'node:test';
import assert from 'node:assert/strict';
import { frameIsFinite, FINITE_STRIDE } from '../js/frames.module.js';

test('a frame is finite unless most of a field or the mean surface temperature is not', () => {
  const temp = new Float32Array(1000).fill(288), cloud = new Float32Array(1000).fill(0.01);
  assert.equal(frameIsFinite({ temp, cloud }, { meanSurfaceT: 288 }), true);
  assert.equal(frameIsFinite({ temp, cloud }), true);
  assert.equal(frameIsFinite({ temp, cloud }, { meanSurfaceT: NaN }), false);
  assert.equal(frameIsFinite({ temp: new Float32Array(1000).fill(NaN) }), false);
  assert.equal(frameIsFinite({ temp: new Float32Array(1000).fill(Infinity) }), false);
});

test('the ocean fields may hold NaN over land, and a third of a field being NaN is not the model gone', () => {
  const sst = new Float32Array(1000).fill(NaN), temp = new Float32Array(1000).fill(288);
  assert.equal(frameIsFinite({ sst, temp }), true);
  const third = Float32Array.from(temp);
  for (let i = 0; i < third.length; i += 3) third[i] = NaN;
  assert.equal(frameIsFinite({ temp: third }), true);
  const most = Float32Array.from(temp);
  for (let i = 0; i < most.length; i++) if (i % 4) most[i] = NaN;
  assert.equal(frameIsFinite({ temp: most }), false);
  assert.ok(FINITE_STRIDE > 1);
});
