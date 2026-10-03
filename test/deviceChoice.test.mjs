import test from 'node:test';
import assert from 'node:assert/strict';
import { pickDevice, hoursPerMinute, isMobileBrowser, probeN, TARGET_RATE, DESKTOP_MAX_N, MOBILE_MAX_N } from '../js/deviceChoice.module.js';

test('an M1 Max GPU (30 ms a step at N=64 with the whole physics) runs N=128', () => {
  const choice = pickDevice({ gpu: { N: 64, ms: 30 } }, 8);
  assert.equal(choice.engine, 'gpu');
  assert.equal(choice.N, 128);
  assert.ok(Math.abs(choice.rate - 23.4) < 0.5, `${choice.rate}`);
});

test('a GPU at 50 ms a step at N=64 runs N=64, and the test on a phone runs at N=32', () => {
  const choice = pickDevice({ gpu: { N: 64, ms: 50 } }, 4);
  assert.deepEqual([choice.engine, choice.N], ['gpu', 64]);
  assert.ok(choice.rate >= TARGET_RATE);
  assert.equal(probeN(DESKTOP_MAX_N), 64);
  assert.equal(probeN(MOBILE_MAX_N), 32);
  const phone = pickDevice({ gpu: { N: 32, ms: 20 } }, 4, { maxN: MOBILE_MAX_N });
  assert.deepEqual([phone.engine, phone.N], ['gpu', 64]);
});

test('a slower GPU falls to N=32', () => {
  assert.equal(pickDevice({ gpu: { N: 64, ms: 400 } }, 4).N, 32);
  assert.equal(pickDevice({ gpu: { N: 32, ms: 120 } }, 4, { maxN: MOBILE_MAX_N }).N, 32);
});

test('without the GPU, the CPU runs what its workers can carry', () => {
  const failed = { N: 64, error: 'no adapter' };
  assert.deepEqual((({ engine, N }) => [engine, N])(pickDevice({ gpu: failed, cpu: { N: 16, ms: 105 } }, 8)), ['cpu', 64]);
  assert.deepEqual((({ engine, N }) => [engine, N])(pickDevice({ gpu: failed, cpu: { N: 16, ms: 250 } }, 8)), ['cpu', 32]);
  assert.equal(pickDevice({ gpu: null, cpu: { N: 16, ms: 105 } }, 1).N, 32);
  assert.equal(pickDevice({ gpu: null, cpu: { N: 16, ms: 200 } }, 1).N, 16);
});

test('hours a minute from the step time', () => {
  assert.ok(Math.abs(hoursPerMinute(64, 14.5) - 388) < 1);
});

test('with neither test usable there is no choice', () => {
  assert.equal(pickDevice({ gpu: { N: 64, error: 'lost' }, cpu: null }, 4), null);
  assert.equal(pickDevice({ gpu: null, cpu: { N: 16, error: 'failed' } }, 4), null);
});

test('a phone or tablet runs N=64 at most, however fast', () => {
  assert.equal(pickDevice({ gpu: { N: 64, ms: 5 } }, 8, { maxN: MOBILE_MAX_N }).N, 64);
  assert.equal(pickDevice({ gpu: { N: 64, ms: 400 } }, 8, { maxN: MOBILE_MAX_N }).N, 32);
  assert.equal(pickDevice({ gpu: { N: 64, ms: 5 } }, 8).N, 128);
});

test('mobile browsers by their hint, their user agent or, for an iPad, their touch points', () => {
  const mac = 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/26.0 Safari/605.1.15';
  assert.equal(isMobileBrowser({ userAgent: mac, maxTouchPoints: 0 }), false);
  assert.equal(isMobileBrowser({ userAgent: mac, maxTouchPoints: 5 }), true);
  assert.equal(isMobileBrowser({ userAgent: 'Mozilla/5.0 (iPhone; CPU iPhone OS 26_0 like Mac OS X) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/26.0 Mobile/15E148 Safari/604.1' }), true);
  assert.equal(isMobileBrowser({ userAgent: 'Mozilla/5.0 (Linux; Android 15; Pixel 9) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/140.0 Mobile Safari/537.36' }), true);
  assert.equal(isMobileBrowser({ userAgent: 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/140.0 Safari/537.36', userAgentData: { mobile: false } }), false);
  assert.equal(isMobileBrowser({ userAgent: 'something', userAgentData: { mobile: true } }), true);
});
