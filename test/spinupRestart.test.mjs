import { test } from 'node:test';
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { copyFileSync, mkdtempSync, readFileSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { decodeState } from '../js/stateFile.module.js';
import { withoutMeans } from '../js/intervalMeans.module.js';
import { withoutOptionGroups } from '../scripts/optionEnvironment.mjs';
import { loadSeaIceStart, SEA_ICE_START_FILE } from '../js/physics/seaIceStart.module.js';
import { stateGrid } from '../js/stateGrid.module.js';
import { stateGeography, createGeography } from '../js/geography.module.js';
import { seaCells, seaIceExtent } from '../js/seaIceExtent.module.js';

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
      const run = (extra) => execFileSync(process.execPath, [join(root, 'scripts/spinup.mjs')], { cwd: root, env: { ...withoutOptionGroups(), N: String(N), MINUTES: '100', OUT: out, ...env, ...extra }, encoding: 'utf8', stdio: ['ignore', 'pipe', 'pipe'] });
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

test('at N=32 spinup.mjs on the n1024 mask with FLake in its lakes, resumed inside a day and run to its end whole, and stopped again inside it with lake outflow not yet drained to the sea, saved, then resumed, ends on the same state bit for bit but for the day\'s rain and radiation means, which cover the steps each run took; and a run from a fresh start, which lays no observed sea-ice snow on a lake, ends the day on the state of its own restart from step 32 but for those means and the energy record, which holds day 0 only for the run that took it whole', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const out = mkdtempSync(join(tmpdir(), 'restart-lakes-'));
  try {
    const run = (extra) => execFileSync(process.execPath, [join(root, 'scripts/spinup.mjs')], { cwd: root, env: { ...withoutOptionGroups(), N: '32', MINUTES: '100', OUT: out, DAYS: '1', ...extra }, encoding: 'utf8', stdio: ['ignore', 'pipe', 'pipe'] });
    run({ TAG: 'whole', STOP_AFTER_STEPS: '32' });
    copyFileSync(join(out, 'whole_day0000_step0032.bin'), join(out, 'split_day0000_step0032.bin'));
    run({ TAG: 'whole' });
    run({ TAG: 'split', STOP_AFTER_STEPS: '32' });
    const inside = await decoded(join(out, 'split_day0000_step0064.bin'));
    const lakes = inside.lakeMean.reduce((n, v) => n + (v !== 0), 0), pending = inside.lakePending.reduce((n, v) => n + (v !== 0), 0);
    assert.ok(lakes > 0 && pending > 0, `${pending} of ${lakes} lake cells hold outflow not yet drained at the stop`);
    run({ TAG: 'split' });
    run({ TAG: 'fresh' });
    const [fresh, whole, split] = await Promise.all(['fresh', 'whole', 'split'].map((tag) => decoded(join(out, `${tag}_day0001.bin`))));
    const means = ['convectiveRain', 'largeScaleRain', 'meanAbsorbedSolar', 'meanOutgoingLongwave', 'meanPlanetaryAlbedo', 'meanShortwaveCloudEffect', 'meanLongwaveCloudEffect'];
    const parted = apart(whole, split, means);
    assert.deepEqual(parted, [], `day 1, whole against split: ${parted.join(', ')}`);
    const start = await loadSeaIceStart(fileURLToPath(new URL(`../${SEA_ICE_START_FILE}`, import.meta.url))), { mesh } = stateGrid(whole);
    const snowed = [...whole.lakeMean.keys()].filter((i) => whole.lakeMean[i] !== 0 && start.at(mesh.latCell[i], mesh.lonCell[i]).snow > 0);
    assert.ok(snowed.length > 0, 'the observed sea ice holds snow at a lake cell');
    assert.deepEqual(snowed.map((i) => fresh.land.snow[i]), snowed.map((i) => whole.land.snow[i]), `day 1, land.snow at lake cells ${snowed.join(', ')}: the fresh run against its restart`);
    const restarted = apart(fresh, whole, [...means, 'energyRecord']);
    assert.deepEqual(restarted, [], `day 1, fresh against whole: ${restarted.join(', ')}`);
  } finally {
    rmSync(out, { recursive: true, force: true });
  }
});

test("at N=64 on the n1024 mask spinup.mjs's daily line gives the sea-ice extents of the day's save, with FLake over the sea cells and saying its frozen lakes are left out, and with the lake model off over every cell off land and saying nothing of lakes", { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const out = mkdtempSync(join(tmpdir(), 'spinup-extents-'));
  try {
    for (const [tag, lake] of [['flake', {}], ['off', { LAKE: '{"model":"off"}' }]]) {
      execFileSync(process.execPath, [join(root, 'scripts/spinup.mjs')], { cwd: root, env: { ...withoutOptionGroups(), N: '64', MINUTES: '100', OUT: out, DAYS: '1', TAG: tag, ...lake }, encoding: 'utf8', stdio: ['ignore', 'pipe', 'pipe'] });
      const line = readFileSync(join(out, `${tag}.log`), 'utf8').split('\n').find((l) => l.startsWith('day 1 ('));
      const printed = line?.match(/, ice [\d.]+% \(N ([\d.]+) S ([\d.]+) Mkm²(, lakes left out)?\), albedo /);
      assert.ok(printed, `${tag}: the day-1 line gives the extents: ${line}`);
      const saved = await decoded(join(out, `${tag}_day0001.bin`)), { mesh } = stateGrid(saved), g = stateGeography(saved), geography = createGeography(mesh, g.topography, g.geography);
      const left = seaIceExtent(mesh, seaCells({ geography, lake: { isLake: geography.lake } }), saved.concentration);
      const counted = seaIceExtent(mesh, seaCells({ geography, lake: null }), saved.concentration);
      console.log(`${tag}: day 1 line N ${printed[1]} S ${printed[2]} Mkm²${printed[3] ?? ''}; the save, lakes left out N ${left[0].toFixed(3)} S ${left[1].toFixed(3)}, counted N ${counted[0].toFixed(3)} S ${counted[1].toFixed(3)}`);
      assert.ok(counted[0] - left[0] >= 0.1, `${tag}: the frozen lakes hold ${(counted[0] - left[0]).toFixed(3)} Mkm² of the northern extent, enough to show in its tenths`);
      const want = tag === 'flake' ? [...left.map((x) => x.toFixed(1)), ', lakes left out'] : [...counted.map((x) => x.toFixed(1)), ''];
      assert.deepEqual([printed[1], printed[2], printed[3] ?? ''], want, `${tag}: the extents and the label`);
    }
  } finally {
    rmSync(out, { recursive: true, force: true });
  }
});
