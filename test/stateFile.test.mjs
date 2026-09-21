import { test } from 'node:test';
import assert from 'node:assert/strict';
import { gzipSync } from 'node:zlib';
import { decodeState, readState, stateName } from '../js/stateFile.module.js';

test('a saved state reads the same whether or not the file is gzipped', async () => {
  const state = { N: 4, day: 12, pi: [1e5, 99999.5], theta: [[300, 301], [302, 303]] };
  const text = JSON.stringify(state);
  assert.deepEqual(await decodeState(new TextEncoder().encode(text)), state);
  assert.deepEqual(await decodeState(new Uint8Array(gzipSync(text))), state);
  assert.deepEqual(await readState(new Response(gzipSync(text))), state);
  await assert.rejects(readState(new Response('', { status: 404 })));
});

test('state names drop the .json or .json.gz extension', () => {
  assert.equal(stateName('fresh64_state_day2190.json.gz'), 'fresh64_state_day2190');
  assert.equal(stateName('fresh64_state_day2190.json'), 'fresh64_state_day2190');
});
