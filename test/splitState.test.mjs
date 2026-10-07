import { test } from 'node:test';
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { mkdtempSync, writeFileSync, readFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { createHash } from 'node:crypto';
import { gunzipSync } from 'node:zlib';
import { runId } from '../js/stateFile.module.js';

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
