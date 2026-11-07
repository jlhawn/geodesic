import { test } from 'node:test';
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';

test('the device test answers with a measured step on a browser-like worker', () => {
  const output = execFileSync(process.execPath, [fileURLToPath(new URL('./helpers/workerProbe.mjs', import.meta.url))], { encoding: 'utf8', stdio: ['ignore', 'pipe', 'ignore'] });
  const reply = JSON.parse(output.trim().split('\n').pop());
  assert.equal(reply.error, undefined);
  assert.ok(reply.result.cpu.ms > 0, JSON.stringify(reply.result));
});
