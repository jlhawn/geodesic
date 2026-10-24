import { test } from 'node:test';
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { mkdtempSync, readFileSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { decodeState } from '../js/stateFile.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }

const root = fileURLToPath(new URL('..', import.meta.url)), N = 16;
const arms = {
  cavity: { ICE_SHELVES: 'cavity', OCEAN: JSON.stringify({ iceDrift: false }) },
  drift: { ICE_SHELVES: 'land', OCEAN: JSON.stringify({ iceDrift: { rheology: 'cavitating', hemisphere: 'both' } }) },
};

const fields = (state, prefix = '') => Object.entries(state).flatMap(([name, value]) => (ArrayBuffer.isView(value) ? [[prefix + name, value]] : value && typeof value === 'object' && !Array.isArray(value) ? fields(value, `${prefix}${name}.`) : []));

for (const [arm, env] of Object.entries(arms)) {
  test(`at N=${N} spinup.mjs with ${arm === 'cavity' ? 'the ice-shelf cavities' : 'the cavitating ice drift'} run 2 days whole and run 1 day, saved, then resumed for the second, ends on the same state bit for bit`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
    const out = mkdtempSync(join(tmpdir(), `restart-${arm}-`));
    try {
      const run = (extra) => execFileSync(process.execPath, [join(root, 'scripts/spinup.mjs')], { cwd: root, env: { ...process.env, N: String(N), MINUTES: '100', OUT: out, ...env, ...extra }, encoding: 'utf8', stdio: ['ignore', 'pipe', 'pipe'] });
      run({ TAG: 'whole', DAYS: '2' });
      run({ TAG: 'split', DAYS: '1' });
      run({ TAG: 'split', DAYS: '2' });
      const [whole, split] = await Promise.all(['whole', 'split'].map((tag) => decodeState(new Uint8Array(readFileSync(join(out, `${tag}_day0002.bin`))))));
      const other = Object.fromEntries(fields(split)), apart = [];
      for (const [name, values] of fields(whole)) {
        const theirs = other[name];
        if (!theirs || theirs.length !== values.length) { apart.push(`${name} missing`); continue; }
        let count = 0, largest = 0;
        for (let i = 0; i < values.length; i++) if (values[i] !== theirs[i] && !(Number.isNaN(values[i]) && Number.isNaN(theirs[i]))) { count++; largest = Math.max(largest, Math.abs(values[i] - theirs[i])); }
        if (count) apart.push(`${name} ${largest.toExponential(1)} (${count})`);
      }
      assert.deepEqual(apart, [], `day 2, whole against split: ${apart.join(', ')}`);
    } finally {
      rmSync(out, { recursive: true, force: true });
    }
  });
}
