import test from 'node:test';
import assert from 'node:assert/strict';
import { pickDevice, hoursPerMinute, TARGET_RATE } from '../js/deviceChoice.module.js';

test('an M1 Max GPU (14.5 ms a step at N=64) runs N=128', () => {
  const choice = pickDevice({ gpu: { N: 64, ms: 14.5 } }, 8);
  assert.equal(choice.engine, 'gpu');
  assert.equal(choice.N, 128);
  assert.ok(Math.abs(choice.rate - 48.5) < 0.5, `${choice.rate}`);
});

test('an iPhone GPU (165 ms a step at N=64) runs N=64', () => {
  const choice = pickDevice({ gpu: { N: 64, ms: 165 } }, 4);
  assert.deepEqual([choice.engine, choice.N], ['gpu', 64]);
  assert.ok(choice.rate >= TARGET_RATE);
});

test('a slower GPU falls to N=32', () => {
  assert.equal(pickDevice({ gpu: { N: 64, ms: 400 } }, 4).N, 32);
});

test('without the GPU, the CPU runs what its workers can carry', () => {
  const failed = { N: 64, error: 'no adapter' };
  assert.deepEqual((({ engine, N }) => [engine, N])(pickDevice({ gpu: failed, cpu: { N: 16, ms: 105 } }, 8)), ['cpu', 32]);
  assert.equal(pickDevice({ gpu: null, cpu: { N: 16, ms: 105 } }, 1).N, 16);
});

test('hours a minute from the step time', () => {
  assert.ok(Math.abs(hoursPerMinute(64, 14.5) - 388) < 1);
});

test('with neither test usable there is no choice', () => {
  assert.equal(pickDevice({ gpu: { N: 64, error: 'lost' }, cpu: null }, 4), null);
  assert.equal(pickDevice({ gpu: null, cpu: { N: 16, error: 'failed' } }, 4), null);
});
