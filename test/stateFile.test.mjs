import { test } from 'node:test';
import assert from 'node:assert/strict';
import { gzipSync, gunzipSync } from 'node:zlib';
import { readFileSync, readdirSync } from 'node:fs';
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

/*
 * fetchState against decodeState of the whole file, through a fetch whose
 * bodies arrive in chunks of `size` bytes (or with no body at all).
 */
function sameBits(actual, expected, path = 'state') {
  if (ArrayBuffer.isView(expected)) {
    assert.ok(ArrayBuffer.isView(actual) && actual.constructor === expected.constructor && actual.length === expected.length, `${path}: ${actual?.constructor?.name} of ${actual?.length} against ${expected.constructor.name} of ${expected.length}`);
    assert.ok(Buffer.from(actual.buffer, actual.byteOffset, actual.byteLength).equals(Buffer.from(expected.buffer, expected.byteOffset, expected.byteLength)), `${path} differs`);
  } else if (expected && typeof expected === 'object') {
    assert.deepEqual(Object.keys(actual).sort(), Object.keys(expected).sort(), `${path} keys`);
    for (const key of Object.keys(expected)) sameBits(actual[key], expected[key], `${path}.${key}`);
  } else assert.equal(actual, expected, path);
}

function chunkedFetch(files, size, { bodyless = false } = {}) {
  return async (url) => {
    const bytes = files[url];
    if (!bytes) return new Response('', { status: 404 });
    if (bodyless) return { ok: true, status: 200, body: null, headers: new Headers({ 'content-length': String(bytes.length) }), arrayBuffer: async () => bytes.slice().buffer, json: async () => JSON.parse(new TextDecoder().decode(bytes)) };
    let at = 0;
    return new Response(new ReadableStream({
      pull(controller) {
        if (at >= bytes.length) { controller.close(); return; }
        controller.enqueue(bytes.slice(at, at + size));
        at += size;
      },
    }), { headers: { 'content-length': String(bytes.length) } });
  };
}

async function withFetch(fetcher, work) {
  const realFetch = globalThis.fetch;
  globalThis.fetch = fetcher;
  try { return await work(); } finally { globalThis.fetch = realFetch; }
}

function partsOf(name, gz, cuts) {
  const files = {}, parts = [];
  let from = 0;
  for (const [n, to] of [...cuts, gz.length].entries()) {
    const file = `${name}.${String(n).padStart(3, '0')}`;
    files[`http://h/runs/${file}`] = gz.subarray(from, to);
    parts.push({ file, bytes: to - from });
    from = to;
  }
  files[`http://h/runs/${name.replace(/\.(?:bin|json)\.gz$/, '')}.parts.json`] = new TextEncoder().encode(JSON.stringify({ parts }));
  return files;
}

test('fetchState decodes small binary and JSON states, gzipped or not and in parts, as decodeState does, through chunks of any size', async () => {
  const encoded = encodeState(binarySample(), { f64: ['pi', 'ocean.eta'] });
  const gz = new Uint8Array(gzipSync(encoded)), jsonGz = new Uint8Array(gzipSync(text));
  const files = {
    'http://h/runs/s_day12.bin': encoded,
    'http://h/runs/s_day12.bin.gz': gz,
    'http://h/runs/a_state_day12.json': new TextEncoder().encode(text),
    'http://h/runs/a_state_day12.json.gz': jsonGz,
    ...partsOf('s_day12.bin.gz', gz, [5, 6, 40]),
    ...partsOf('a_state_day12.json.gz', jsonGz, [1, 30]),
  };
  const expected = { bin: await decodeState(encoded), json: await decodeState(new TextEncoder().encode(text)) };
  for (const size of [1, 2, 3, 7, 8, 9, 64, 4096]) {
    for (const bodyless of size === 4096 ? [false, true] : [false]) {
      await withFetch(chunkedFetch(files, size, { bodyless }), async () => {
        for (const name of ['s_day12.bin', 's_day12.bin.gz', 's_day12.parts.json']) {
          const seen = [];
          sameBits(await fetchState(`http://h/runs/${name}`, (received, total) => seen.push([received, total])), expected.bin, `${name} in chunks of ${size}${bodyless ? ', no body' : ''}`);
          assert.equal(seen.at(-1)[0], name === 's_day12.bin' ? encoded.length : gz.length);
          assert.ok(seen.every(([, total]) => total === seen.at(-1)[0]));
        }
        for (const name of ['a_state_day12.json', 'a_state_day12.json.gz', 'a_state_day12.parts.json']) sameBits(await fetchState(`http://h/runs/${name}`), expected.json, `${name} in chunks of ${size}`);
      });
    }
  }
});

test('fetchState rejects a missing part, a state cut short and one longer than its arrays', async () => {
  const encoded = encodeState(binarySample(), { f64: ['pi'] });
  const gz = new Uint8Array(gzipSync(encoded));
  const files = partsOf('s_day12.bin.gz', gz, [10, 20]);
  const missing = { ...files };
  delete missing['http://h/runs/s_day12.bin.gz.001'];
  const longer = new Uint8Array(encoded.length + 4);
  longer.set(encoded);
  const shortGz = new Uint8Array(gzipSync(encoded.subarray(0, encoded.length - 4)));
  await withFetch(chunkedFetch({ ...missing, 'http://h/runs/short.bin': encoded.subarray(0, encoded.length - 4), 'http://h/runs/short.bin.gz': shortGz, 'http://h/runs/long.bin': longer, 'http://h/runs/long.bin.gz': new Uint8Array(gzipSync(longer)) }, 7), async () => {
    await assert.rejects(fetchState('http://h/runs/s_day12.parts.json'), /s_day12\.bin\.gz\.001: 404/);
    await assert.rejects(fetchState('http://h/runs/short.bin'), /runs past the end/);
    await assert.rejects(fetchState('http://h/runs/short.bin.gz'), /runs past the end/);
    await assert.rejects(decodeState(encoded.subarray(0, encoded.length - 4)), /runs past the end/);
    await assert.rejects(fetchState('http://h/runs/long.bin'), /runs past the end/);
    await assert.rejects(fetchState('http://h/runs/long.bin.gz'), /runs past the end/);
  });
});

test("fetchState reads the page's default states in parts as decodeState reads the whole file, and eleven64's as one .bin.gz and as a .bin", async () => {
  const runs = new URL('../runs/', import.meta.url);
  const manifests = readdirSync(runs).filter((file) => file.endsWith('.parts.json')).sort();
  assert.ok(manifests.length >= 2, `parts in runs/: ${manifests}`);
  for (const manifest of manifests) {
    const { parts } = JSON.parse(readFileSync(new URL(manifest, runs)));
    const files = { [`http://h/runs/${manifest}`]: readFileSync(new URL(manifest, runs)) };
    for (const { file } of parts) files[`http://h/runs/${file}`] = new Uint8Array(readFileSync(new URL(file, runs)));
    const whole = () => Buffer.concat(parts.map(({ file }) => files[`http://h/runs/${file}`]));
    const expected = await decodeState(whole());
    assert.ok(expected.N >= 64 && expected.theta instanceof Float32Array);
    await withFetch(chunkedFetch(files, 1 << 20), async () => sameBits(await fetchState(`http://h/runs/${manifest}`), expected, manifest));
    if (!manifest.startsWith('eleven64')) continue;
    const name = manifest.replace('.parts.json', ''), gz = whole();
    const single = { [`http://h/runs/${name}.bin.gz`]: gz, [`http://h/runs/${name}.bin`]: new Uint8Array(gunzipSync(gz)) };
    await withFetch(chunkedFetch(single, 65536 + 3), async () => {
      sameBits(await fetchState(`http://h/runs/${name}.bin.gz`), expected, `${name}.bin.gz`);
      sameBits(await fetchState(`http://h/runs/${name}.bin`), expected, `${name}.bin`);
    });
  }
});

test('a directory index lists JSON, parts and binary states, and their days', () => {
  const html = ['layered64_state_day910.json', 'ocean64_state_day930.parts.json', 'ocean64_state_day930.json.gz.000', 'spin128c_day0678.bin', 'spin128c_day0650.bin.gz', 'layered64_day030.json', 'spin128c.log']
    .map((file) => `<li><a href="${encodeURIComponent(file)}">${file}</a></li>`).join('\n');
  assert.deepEqual(listedStates(html), ['layered64_state_day910.json', 'ocean64_state_day930.parts.json', 'spin128c_day0650.bin.gz', 'spin128c_day0678.bin']);
  assert.equal(stateDay('runs/spin128c_day0678.bin'), 678);
  assert.equal(stateDay('runs/layered64_state_day910.json'), 910);
  assert.equal(stateDay('runs/notes.txt'), null);
});
