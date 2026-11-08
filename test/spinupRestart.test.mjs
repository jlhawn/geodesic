import { test } from 'node:test';
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { copyFileSync, mkdtempSync, readFileSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { decodeState } from '../js/stateFile.module.js';
import { withoutMeans } from '../js/intervalMeans.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }

const root = fileURLToPath(new URL('..', import.meta.url)), N = 16;
const arms = {
  cavity: { ICE_SHELVES: 'cavity', OCEAN: JSON.stringify({ iceDrift: false }) },
  drift: { ICE_SHELVES: 'land', OCEAN: JSON.stringify({ iceDrift: { rheology: 'cavitating', hemisphere: 'both' } }) },
};

const decoded = (file) => decodeState(new Uint8Array(readFileSync(file)));
const fields = (state, prefix = '') => Object.entries(state).flatMap(([name, value]) => (ArrayBuffer.isView(value) ? [[prefix + name, value]] : value && typeof value === 'object' && !Array.isArray(value) ? fields(value, `${prefix}${name}.`) : []));
function apart(whole, split, except = []) {
  const other = Object.fromEntries(fields(withoutMeans(split))), out = [];
  for (const [name, values] of fields(withoutMeans(whole))) {
    if (except.includes(name)) continue;
    const theirs = other[name];
    if (!theirs || theirs.length !== values.length) { out.push(`${name} missing`); continue; }
    let count = 0, largest = 0;
    for (let i = 0; i < values.length; i++) if (values[i] !== theirs[i] && !(Number.isNaN(values[i]) && Number.isNaN(theirs[i]))) { count++; largest = Math.max(largest, Math.abs(values[i] - theirs[i])); }
    if (count) out.push(`${name} ${largest.toExponential(1)} (${count})`);
  }
  return out;
}

for (const [arm, env] of Object.entries(arms)) {
  test(`at N=${N} spinup.mjs with ${arm === 'cavity' ? 'the ice-shelf cavities' : 'the cavitating ice drift'} run 2 days whole and run 1 day, saved, then resumed for the second, ends on the same state bit for bit`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
    const out = mkdtempSync(join(tmpdir(), `restart-${arm}-`));
    try {
      const run = (extra) => execFileSync(process.execPath, [join(root, 'scripts/spinup.mjs')], { cwd: root, env: { ...process.env, N: String(N), MINUTES: '100', OUT: out, ...env, ...extra }, encoding: 'utf8', stdio: ['ignore', 'pipe', 'pipe'] });
      run({ TAG: 'whole', DAYS: '2' });
      run({ TAG: 'split', DAYS: '1' });
      run({ TAG: 'split', DAYS: '2' });
      const [whole, split] = await Promise.all(['whole', 'split'].map((tag) => decoded(join(out, `${tag}_day0002.bin`))));
      assert.deepEqual([whole.means.days, whole.means.samples, split.means.days, split.means.samples], [2, 128, 1, 64], 'each save\'s means cover its own segment');
      const parted = apart(whole, split);
      assert.deepEqual(parted, [], `day 2, whole against split: ${parted.join(', ')}`);
    } finally {
      rmSync(out, { recursive: true, force: true });
    }
  });
}

test('at N=32 spinup.mjs on the n1024 mask with FLake in its lakes, resumed inside a day and run to its end whole, and stopped again inside it with lake outflow not yet drained to the sea, saved, then resumed, ends on the same state bit for bit but for the day\'s rain and radiation means and its energy record, which cover the steps each run took', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const out = mkdtempSync(join(tmpdir(), 'restart-lakes-'));
  try {
    const run = (extra) => execFileSync(process.execPath, [join(root, 'scripts/spinup.mjs')], { cwd: root, env: { ...process.env, N: '32', MINUTES: '100', OUT: out, DAYS: '1', ...extra }, encoding: 'utf8', stdio: ['ignore', 'pipe', 'pipe'] });
    run({ TAG: 'whole', STOP_AFTER_STEPS: '32' });
    copyFileSync(join(out, 'whole_day0000_step0032.bin'), join(out, 'split_day0000_step0032.bin'));
    run({ TAG: 'whole' });
    run({ TAG: 'split', STOP_AFTER_STEPS: '32' });
    const inside = await decoded(join(out, 'split_day0000_step0064.bin'));
    const lakes = inside.lakeMean.reduce((n, v) => n + (v !== 0), 0), pending = inside.lakePending.reduce((n, v) => n + (v !== 0), 0);
    assert.ok(lakes > 0 && pending > 0, `${pending} of ${lakes} lake cells hold outflow not yet drained at the stop`);
    run({ TAG: 'split' });
    const [whole, split] = await Promise.all(['whole', 'split'].map((tag) => decoded(join(out, `${tag}_day0001.bin`))));
    const parted = apart(whole, split, ['energyRecord', 'convectiveRain', 'largeScaleRain', 'meanAbsorbedSolar', 'meanOutgoingLongwave', 'meanPlanetaryAlbedo', 'meanShortwaveCloudEffect', 'meanLongwaveCloudEffect']);
    assert.deepEqual(parted, [], `day 1, whole against split: ${parted.join(', ')}`);
  } finally {
    rmSync(out, { recursive: true, force: true });
  }
});
