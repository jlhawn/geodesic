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

test(`at N=${N} spinup.mjs resumed from a state saved under OCEAN_CD keeps the C-D Coriolis scheme without the switch and its relaxation time, and v stays nearer R(u) than u is large`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const out = mkdtempSync(join(tmpdir(), 'restart-cd-'));
  try {
    const run = (extra) => execFileSync(process.execPath, [join(root, 'scripts/spinup.mjs')], { cwd: root, env: { ...process.env, N: String(N), MINUTES: '100', OUT: out, TAG: 'cd', ICE_SHELVES: 'land', OCEAN: JSON.stringify({ iceDrift: false }), ...extra }, encoding: 'utf8', stdio: ['ignore', 'pipe', 'pipe'] });
    run({ DAYS: '1', OCEAN_CD: '1', OCEAN_CD_TAU: '200000' });
    const first = decodeState(new Uint8Array(readFileSync(join(out, 'cd_day0001.bin'))));
    const resumed = run({ DAYS: '2' });
    const [day1, day2] = await Promise.all([first, decodeState(new Uint8Array(readFileSync(join(out, 'cd_day0002.bin'))))]);
    for (const state of [day1, day2]) {
      assert.equal(state.ocean.coriolis, 'cd');
      assert.equal(state.ocean.cdTau, 200000);
    }
    const lines = resumed.split('\n').filter((line) => line.startsWith('C-D Coriolis'));
    assert.ok(lines.length >= 1, 'the resumed segment reports the scheme');
    for (const line of lines) {
      const [misfit, speed] = [...line.matchAll(/(\d\.\d+e[-+]\d+) m\/s/g)].map((m) => Number(m[1]));
      assert.ok(misfit < speed, `${line}`);
    }
  } finally {
    rmSync(out, { recursive: true, force: true });
  }
});
