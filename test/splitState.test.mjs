import { test } from 'node:test';
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { mkdtempSync, writeFileSync, readFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { createHash } from 'node:crypto';
import { gunzipSync } from 'node:zlib';
import { runId, encodeState, decodeState } from '../js/stateFile.module.js';
import { SST_FIELDS } from '../js/sstMeans.module.js';

const SCRIPT = new URL('../scripts/splitState.mjs', import.meta.url).pathname;
const hex = (bytes) => createHash('sha256').update(bytes).digest('hex');

test('splitState writes hashed parts with a manifest naming the run by its content, and rewrites a manifest in place', async () => {
  const dir = mkdtempSync(join(tmpdir(), 'split-'));
  const state = { N: 4, day: 7, pi: Array.from({ length: 4000 }, (_, i) => Math.sin(i)) };
  writeFileSync(join(dir, 'c_state_day7.json'), JSON.stringify(state));
  execFileSync('node', [SCRIPT, join(dir, 'c_state_day7.json'), String(1 / 1024)]);
  const manifest = JSON.parse(readFileSync(join(dir, 'c_state_day7.parts.json'), 'utf8'));
  assert.equal(manifest.name, 'c_state_day7');
  assert.ok(manifest.parts.length > 1);
  for (const part of manifest.parts) {
    const bytes = readFileSync(join(dir, part.file));
    assert.equal(part.bytes, bytes.length);
    assert.equal(part.sha256, hex(bytes));
  }
  assert.equal(manifest.id, await runId(manifest.parts));
  assert.deepEqual(JSON.parse(gunzipSync(Buffer.concat(manifest.parts.map((part) => readFileSync(join(dir, part.file))))).toString()), state);
  const written = readFileSync(join(dir, 'c_state_day7.parts.json'), 'utf8');
  execFileSync('node', [SCRIPT, join(dir, 'c_state_day7.parts.json')]);
  assert.equal(readFileSync(join(dir, 'c_state_day7.parts.json'), 'utf8'), written);
});

test('splitState leaves a binary state\'s SST ring out of its parts and the rest as it was saved, and splits a state without one byte for byte', async () => {
  const dir = mkdtempSync(join(tmpdir(), 'split-'));
  const C = 42, values = (length, from) => Float64Array.from({ length }, (_, i) => from + i / 3);
  const state = {
    N: 2, K: 3, day: 12, time: 12 * 86400, terrain: true, grid: { relax: 30 },
    pi: values(C, 1e5), theta: Float32Array.from(values(3 * C, 290)), levels: Float64Array.of(0, 0.3, 0.7, 1), energyRecord: values(6, 1),
    ocean: { vertical: 'hybrid', generator: 'hybrid', h: values(2 * C, 10), T: Float32Array.from(values(2 * C, 280)), densities: Float64Array.of(1024.5, 1027.25) },
    land: { snow: Float32Array.from(values(C, 0)), record: values(4, 2) },
  };
  const f64 = ['pi', 'ocean.h'], ring = { sstMeans: Float32Array.from(values(2 * C, 280)), sstMeanDays: Float64Array.of(11, 12), sstSum: values(C, 0), sstSamples: 5 };
  const split = (name, bytes) => {
    writeFileSync(join(dir, `${name}.bin`), bytes);
    execFileSync(process.execPath, [SCRIPT, join(dir, `${name}.bin`), String(1 / 1024)]);
    const { parts } = JSON.parse(readFileSync(join(dir, `${name}.parts.json`), 'utf8'));
    return gunzipSync(Buffer.concat(parts.map((part) => readFileSync(join(dir, part.file)))));
  };
  const plain = encodeState(state, { f64 });
  assert.ok(split('plain_day0012', plain).equals(plain));
  const parted = await decodeState(new Uint8Array(split('ringed_day0012', encodeState({ ...state, ...ring }, { f64 }))));
  for (const name of SST_FIELDS) assert.ok(!(name in parted), name);
  assert.deepEqual(parted, await decodeState(plain));
});

test('splitState records RADIATION and MOIST in the header of a binary state\'s parts as its options, and keeps the arrays as they were saved', async () => {
  const dir = mkdtempSync(join(tmpdir(), 'split-'));
  const C = 42, values = (length, from) => Float64Array.from({ length }, (_, i) => from + i / 3);
  const state = { N: 2, K: 3, day: 12, time: 12 * 86400, grid: { relax: 30 }, pi: values(C, 1e5), theta: Float32Array.from(values(3 * C, 290)), ocean: { vertical: 'hybrid', h: values(2 * C, 10) } };
  const f64 = ['pi', 'ocean.h'], radiation = { cloudColumns: true, cloudSpread: 0.75 }, moist = { stratocumulusInversion: 7, cumulusAscent: 'kinetic' };
  const plain = encodeState(state, { f64 });
  writeFileSync(join(dir, 'f_day0012.bin'), plain);
  execFileSync(process.execPath, [SCRIPT, join(dir, 'f_day0012.bin'), String(1 / 1024)], { env: { ...process.env, RADIATION: JSON.stringify(radiation), MOIST: JSON.stringify(moist) } });
  const { parts } = JSON.parse(readFileSync(join(dir, 'f_day0012.parts.json'), 'utf8'));
  const parted = await decodeState(new Uint8Array(gunzipSync(Buffer.concat(parts.map((part) => readFileSync(join(dir, part.file)))))));
  assert.deepEqual(parted.options, { radiation, moist });
  delete parted.options;
  assert.deepEqual(parted, await decodeState(plain));
});
