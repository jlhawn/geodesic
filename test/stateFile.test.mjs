import { test } from 'node:test';
import assert from 'node:assert/strict';
import { gzipSync } from 'node:zlib';
import { decodeState, encodeState, fetchState, stateName, listedStates, stateDay } from '../js/stateFile.module.js';

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

test('state names drop the .json, .json.gz, .parts.json, .bin or .bin.gz extension', () => {
  assert.equal(stateName('fresh64_state_day2190.parts.json'), 'fresh64_state_day2190');
  assert.equal(stateName('fresh64_state_day2190.json.gz'), 'fresh64_state_day2190');
  assert.equal(stateName('fresh64_state_day2190.json'), 'fresh64_state_day2190');
  assert.equal(stateName('spin128c_day0678.bin'), 'spin128c_day0678');
  assert.equal(stateName('spin128c_day0678.bin.gz'), 'spin128c_day0678');
});

const binarySample = () => ({
  N: 4, K: 3, day: 12.5, time: 1080000, terrain: true,
  pi: Float64Array.from({ length: 7 }, (_, i) => 98000 + 13.25 * i),
  theta: Float64Array.from({ length: 21 }, (_, i) => 280 + Math.sin(i)),
  ice: new Float64Array(7),
  ocean: { h: Float64Array.from({ length: 14 }, (_, i) => 50 + i), eta: Float64Array.from({ length: 7 }, (_, i) => 0.01 * i - 0.03) },
  land: { soil: Float64Array.from({ length: 7 }, (_, i) => 10 * i), snow: new Float64Array(7) },
});

function assertClose(actual, expected, label) {
  assert.equal(actual.length, expected.length, `${label} length`);
  for (let i = 0; i < expected.length; i++) assert.ok(Math.abs(actual[i] - expected[i]) <= 1e-7 * Math.max(1, Math.abs(expected[i])), `${label}[${i}]: ${actual[i]} against ${expected[i]}`);
}

test('a binary state decodes to the JSON state shape, float32 by default and float64 where asked, gzipped or not', async () => {
  const saved = binarySample();
  const encoded = encodeState(saved, { f64: ['pi', 'ocean.eta'] });
  for (const bytes of [encoded, new Uint8Array(gzipSync(encoded))]) {
    const decoded = await decodeState(bytes);
    assert.deepEqual({ N: decoded.N, K: decoded.K, day: decoded.day, time: decoded.time, terrain: decoded.terrain }, { N: 4, K: 3, day: 12.5, time: 1080000, terrain: true });
    assert.ok(decoded.pi instanceof Float64Array && decoded.theta instanceof Float32Array && decoded.ocean.eta instanceof Float64Array && decoded.land.soil instanceof Float32Array);
    assert.deepEqual(Array.from(decoded.pi), Array.from(saved.pi));
    assert.deepEqual(Array.from(decoded.ocean.eta), Array.from(saved.ocean.eta));
    assertClose(decoded.theta, saved.theta, 'theta');
    assertClose(decoded.ocean.h, saved.ocean.h, 'ocean.h');
    assertClose(decoded.land.soil, saved.land.soil, 'land.soil');
    assert.equal(decoded.arrays, undefined);
  }
});

test("a binary state in the spin-up driver's layout, at a misaligned offset, or cut short", async () => {
  const header = new TextEncoder().encode(JSON.stringify({ N: 2, day: 3, arrays: [{ name: 'pi', length: 3, offset: 0 }, { name: 'ocean.h', length: 2, offset: 12 }] }) + ' ');
  assert.notEqual(header.length % 8, 0);
  const start = 8 + Math.ceil(header.length / 4) * 4;
  const bytes = new Uint8Array(start + 20);
  bytes.set([0x47, 0x43, 0x4d, 0x53]);
  new DataView(bytes.buffer).setUint32(4, header.length, true);
  bytes.set(header, 8);
  new Float32Array(bytes.buffer, start, 5).set([1.5, 2.5, 3.5, 60, 70]);
  const decoded = await decodeState(bytes);
  assert.deepEqual(Array.from(decoded.pi), [1.5, 2.5, 3.5]);
  assert.deepEqual(Array.from(decoded.ocean.h), [60, 70]);

  const encoded = encodeState(binarySample(), { f64: ['pi'] });
  const shifted = new Uint8Array(encoded.length + 4);
  shifted.set(encoded, 4);
  const fromShifted = await decodeState(shifted.subarray(4));
  assert.deepEqual(Array.from(fromShifted.pi), Array.from(binarySample().pi));
  assertClose(fromShifted.theta, binarySample().theta, 'shifted theta');

  await assert.rejects(decodeState(encoded.subarray(0, encoded.length - 4)), /runs past the end/);
});

test('a directory index lists JSON, parts and binary states, and their days', () => {
  const html = ['layered64_state_day910.json', 'ocean64_state_day930.parts.json', 'ocean64_state_day930.json.gz.000', 'spin128c_day0678.bin', 'spin128c_day0650.bin.gz', 'layered64_day030.json', 'spin128c.log']
    .map((file) => `<li><a href="${encodeURIComponent(file)}">${file}</a></li>`).join('\n');
  assert.deepEqual(listedStates(html), ['layered64_state_day910.json', 'ocean64_state_day930.parts.json', 'spin128c_day0650.bin.gz', 'spin128c_day0678.bin']);
  assert.equal(stateDay('runs/spin128c_day0678.bin'), 678);
  assert.equal(stateDay('runs/layered64_state_day910.json'), 910);
  assert.equal(stateDay('runs/notes.txt'), null);
});
