import { test } from 'node:test';
import assert from 'node:assert/strict';
import { gzipSync } from 'node:zlib';
import { decodeState, fetchState, stateName } from '../js/stateFile.module.js';

const state = { N: 4, day: 12, pi: [1e5, 99999.5], theta: [[300, 301], [302, 303]] };
const text = JSON.stringify(state);

test('a saved state decodes the same whether or not its bytes are gzipped', async () => {
  assert.deepEqual(await decodeState(new TextEncoder().encode(text)), state);
  assert.deepEqual(await decodeState(new Uint8Array(gzipSync(text))), state);
});

test('a state cut into gzip parts is fetched through its manifest with a byte-accurate progress total', async () => {
  const gz = gzipSync(text);
  const cut = Math.floor(gz.length / 2);
  const files = {
    'http://h/runs/a_state_day12.json.gz.000': gz.subarray(0, cut),
    'http://h/runs/a_state_day12.json.gz.001': gz.subarray(cut),
    'http://h/runs/a_state_day12.parts.json': new TextEncoder().encode(JSON.stringify({ parts: [{ file: 'a_state_day12.json.gz.000', bytes: cut }, { file: 'a_state_day12.json.gz.001', bytes: gz.length - cut }] })),
    'http://h/runs/a_state_day12.json': new TextEncoder().encode(text),
  };
  const realFetch = globalThis.fetch;
  globalThis.fetch = async (url) => (files[url] ? new Response(files[url]) : new Response('', { status: 404 }));
  try {
    const seen = [];
    assert.deepEqual(await fetchState('http://h/runs/a_state_day12.parts.json', (received, total) => seen.push([received, total])), state);
    assert.equal(seen.at(-1)[0], gz.length);
    assert.ok(seen.every(([, total]) => total === gz.length));
    assert.deepEqual(await fetchState('http://h/runs/a_state_day12.json'), state);
    await assert.rejects(fetchState('http://h/runs/missing_state_day1.json'), /404/);
  } finally { globalThis.fetch = realFetch; }
});

test('state names drop the .json, .json.gz or .parts.json extension', () => {
  assert.equal(stateName('fresh64_state_day2190.parts.json'), 'fresh64_state_day2190');
  assert.equal(stateName('fresh64_state_day2190.json.gz'), 'fresh64_state_day2190');
  assert.equal(stateName('fresh64_state_day2190.json'), 'fresh64_state_day2190');
});
