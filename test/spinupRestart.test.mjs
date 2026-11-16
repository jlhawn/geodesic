import { test } from 'node:test';
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { copyFileSync, mkdtempSync, readFileSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { decodeState } from '../js/stateFile.module.js';
import { createEnergyRecord, FIELDS } from '../js/physics/energyRecord.module.js';
import { withoutMeans } from '../js/intervalMeans.module.js';
import { withoutOptionGroups } from '../scripts/optionEnvironment.mjs';
import { loadSeaIceStart, SEA_ICE_START_FILE } from '../js/physics/seaIceStart.module.js';
import { stateGrid } from '../js/stateGrid.module.js';
import { stateGeography, createGeography } from '../js/geography.module.js';
import { seaCells, seaIceExtent, polarIcedCells } from '../js/seaIceExtent.module.js';
import { DEFAULT_ICE_EXCHANGE, DEFAULT_ICE_ALBEDO_PHASES } from '../js/physics/ice.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }

const root = fileURLToPath(new URL('..', import.meta.url)), N = 16;
const arms = {
  cavity: { ICE_SHELVES: 'cavity', OCEAN: JSON.stringify({ iceDrift: false }) },
  drift: { ICE_SHELVES: 'land', OCEAN: JSON.stringify({ iceDrift: { rheology: 'cavitating', hemisphere: 'both' } }) },
  deposition: { MOIST: JSON.stringify({ iceDeposition: true }) },
  memory: { MOIST: JSON.stringify({ iceDeposition: true }), RADIATION: JSON.stringify({ cloudPhase: 'memory' }) },
};

const decoded = (file) => decodeState(new Uint8Array(readFileSync(file)));
const fields = (state, prefix = '') => Object.entries(state).flatMap(([name, value]) => (ArrayBuffer.isView(value) ? [[prefix + name, value]] : value && typeof value === 'object' && !Array.isArray(value) ? fields(value, `${prefix}${name}.`) : []));
const gapRows = (length, days) => { const record = createEnergyRecord((length - 2) / FIELDS.length); for (const day of days) record.add(day, null); return record.values.map(Number.isNaN); };
function apart(whole, split, except = [], gapDays = []) {
  const other = Object.fromEntries(fields(withoutMeans(split))), out = [];
  for (const [name, values] of fields(withoutMeans(whole))) {
    if (except.includes(name)) continue;
    const theirs = other[name];
    if (!theirs || theirs.length !== values.length) { out.push(`${name} missing`); continue; }
    const gaps = name === 'energyRecord' ? gapRows(values.length, gapDays) : null;
    let count = 0, largest = 0;
    for (let i = 0; i < values.length; i++) if ((values[i] !== theirs[i] || !Number.isFinite(values[i])) && !(gaps && gaps[i] && Number.isNaN(values[i]) && Number.isNaN(theirs[i]))) { count++; largest = Math.max(largest, Math.abs(values[i] - theirs[i])); }
    if (count) out.push(`${name} ${largest.toExponential(1)} (${count})`);
  }
  return out;
}

for (const [arm, env] of Object.entries(arms)) {
  test(`at N=${N} spinup.mjs with ${{ cavity: 'the ice-shelf cavities', drift: 'the cavitating ice drift', deposition: 'the deposition-limited ice fall', memory: "the deposition-limited ice fall under cloudPhase 'memory'" }[arm]} run 2 days whole and run 1 day, saved, then resumed for the second, ends on the same state bit for bit`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
    const out = mkdtempSync(join(tmpdir(), `restart-${arm}-`));
    try {
      const run = (extra) => execFileSync(process.execPath, [join(root, 'scripts/spinup.mjs')], { cwd: root, env: { ...withoutOptionGroups(), N: String(N), MINUTES: '100', OUT: out, ...env, ...extra }, encoding: 'utf8', stdio: ['ignore', 'pipe', 'pipe'] });
      run({ TAG: 'whole', DAYS: '2' });
      run({ TAG: 'split', DAYS: '1' });
      run({ TAG: 'split', DAYS: '2' });
      const [whole, split] = await Promise.all(['whole', 'split'].map((tag) => decoded(join(out, `${tag}_day0002.bin`))));
      assert.deepEqual([whole.means.days, whole.means.samples, split.means.days, split.means.samples], [2, 128, 1, 64], 'each save\'s means cover its own segment');
      if (arm === 'deposition' || arm === 'memory') assert.ok(whole.cloudIceShare && split.cloudIceShare && whole.cloudIceShare.some((v) => v > 0 && v < 1), 'both saves hold a glaciated fraction with mixed layers');
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
    const parted = apart(whole, split, means, [1]);
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

test("at N=16 spinup.mjs records the default ice exchange in its header, and a run that resumes a save under another exchange says so and takes its own: the default save resumed under OCEAN iceExchange instant holds no water heat and keeps the save's iceAlbedoPhases, and that save resumed under OCEAN iceExchange turbulent, as confirm128 started, holds the water's heat by the end of day 2", { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const out = mkdtempSync(join(tmpdir(), 'restart-switch-')), turbulent = { mode: 'turbulent', frictionMinimum: 2e-2 };
  try {
    const run = (extra) => execFileSync(process.execPath, [join(root, 'scripts/spinup.mjs')], { cwd: root, env: Object.fromEntries(Object.entries({ ...withoutOptionGroups(), N: '16', MINUTES: '100', OUT: out, TAG: 'switch', ...extra }).filter(([, value]) => value !== undefined)), encoding: 'utf8', stdio: ['ignore', 'pipe', 'pipe'] });
    run({ DAYS: '1', STOP_AFTER_STEPS: '32' });
    const made = await decoded(join(out, 'switch_day0000_step0032.bin')), C = made.surfaceT.length;
    assert.deepEqual([made.given.ocean.iceExchange, made.given.ice.iceAlbedoPhases], [DEFAULT_ICE_EXCHANGE, DEFAULT_ICE_ALBEDO_PHASES], 'the header records the default exchange and iceAlbedoPhases');
    assert.equal(made.ocean.iceWater?.length, C);
    assert.match(run({ DAYS: '1', OCEAN: JSON.stringify({ iceExchange: 'instant' }) }), /iceExchange "instant" where switch_day0000_step0032\.bin was made under \{"mode":"turbulent","frictionMinimum":0\.02\}/);
    const day1 = await decoded(join(out, 'switch_day0001.bin'));
    assert.deepEqual([day1.given.ocean.iceExchange, day1.given.ice.iceAlbedoPhases], ['instant', DEFAULT_ICE_ALBEDO_PHASES], "the resumed run takes its own exchange and its save's iceAlbedoPhases");
    assert.equal(day1.ocean.iceWater, undefined, 'the instant run holds no water heat');
    assert.match(run({ DAYS: '2', OCEAN: JSON.stringify({ iceExchange: turbulent }) }), /iceExchange \{"mode":"turbulent","frictionMinimum":0\.02\} where switch_day0001\.bin was made under "instant"/);
    const day2 = await decoded(join(out, 'switch_day0002.bin'));
    assert.deepEqual(day2.given.ocean.iceExchange, turbulent);
    const held = day2.ocean.iceWater?.length === C ? day2.ocean.iceWater.filter((x) => x > 0).length : 0;
    assert.ok(held > 0, `the turbulent run holds the water's heat under its ice in ${held} cells`);
    assert.ok(day2.surfaceT.every(Number.isFinite) && day2.ocean.T.every(Number.isFinite));
  } finally {
    rmSync(out, { recursive: true, force: true });
  }
});

test("at N=64 on the n1024 mask spinup.mjs's daily line gives the sea-ice extents of the day's save, with FLake over the sea cells and saying its frozen lakes are left out, and with the lake model off over every cell off land and saying nothing of lakes, and a surface shortwave over the save's iced sea cells poleward of 60° within a factor of 2 of N 17.3 and S 41.9 W/m²",{ skip: !gpuAvailable && 'webgpu not installed' }, async () => {
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
      const iced = polarIcedCells(mesh, seaCells({ geography, lake: tag === 'flake' ? { isLake: geography.lake } : null }), saved.ice), north = iced.filter((i) => mesh.latCell[i] > 0).length;
      const shortwave = line.match(/surface shortwave ([-\d.]+) net longwave [-\d.]+ W\/m² \(iced cells poleward of 60°: N ([-\d.]+) S ([-\d.]+)\)/);
      console.log(`${tag}: iced sea cells poleward of 60° in the save N ${north} S ${iced.length - north}; the line gives ${shortwave?.[0]}`);
      assert.ok(north > 0 && iced.length > north, `${tag}: the save holds sea ice poleward of 60° in each hemisphere`);
      assert.ok(shortwave && [[shortwave[2], 17.3], [shortwave[3], 41.9]].every(([x, near]) => +x > near / 2 && +x < near * 2), `${tag}: the line gives the iced cells poleward of 60° a surface shortwave within a factor of 2 of N 17.3 and S 41.9 W/m²: ${shortwave?.[0]}`);
    }
  } finally {
    rmSync(out, { recursive: true, force: true });
  }
});

test("at N=16 spinup.mjs with LAND soilGoal 'observed', stopped inside a day and resumed, carries the goal's window and ends day 2 on the land of the run that did not stop, bit for bit, whose window opened on day 0; a segment that ends at day 1 closes its window, and the next one's holds the second day alone and says so in its header", { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const out = mkdtempSync(join(tmpdir(), 'restart-goal-'));
  try {
    const run = (extra) => execFileSync(process.execPath, [join(root, 'scripts/spinup.mjs')], { cwd: root, env: { ...withoutOptionGroups(), N: '16', MINUTES: '100', OUT: out, LAND: JSON.stringify({ soilGoal: 'observed' }), ...extra }, encoding: 'utf8', stdio: ['ignore', 'pipe', 'pipe'] });
    run({ TAG: 'whole', DAYS: '2' });
    run({ TAG: 'halt', DAYS: '2', STOP_AFTER_STEPS: '96' });
    const inside = await decoded(join(out, 'halt_day0001_step0032.bin'));
    assert.deepEqual([inside.goalSince, inside.goalOpen], [0, true], 'a save inside a day leaves the window open');
    assert.ok(inside.land.goalAdded.some((x) => x > 0), 'and carries its sums');
    run({ TAG: 'halt', DAYS: '2' });
    run({ TAG: 'close', DAYS: '1' });
    run({ TAG: 'close', DAYS: '2' });
    const [whole, halt, close] = await Promise.all(['whole', 'halt', 'close'].map((tag) => decoded(join(out, `${tag}_day0002.bin`))));
    assert.deepEqual([whole.goalSince, whole.goalOpen, halt.goalSince, halt.goalOpen], [0, false, 0, false]);
    const land = (state) => Object.fromEntries(fields(state.land));
    const differing = Object.entries(land(whole)).filter(([name, values]) => !land(halt)[name] || values.some((x, i) => x !== land(halt)[name][i])).map(([name]) => name);
    assert.deepEqual(differing, [], `day 2 land, whole against halted: ${differing.join(', ')}`);
    assert.deepEqual([close.goalSince, close.goalOpen], [1, false], 'the closed window of day 1 is not carried');
    const sum = (state, name) => state.land[name].reduce((s, x) => s + x, 0);
    assert.ok(sum(close, 'goalAdded') > 0 && sum(close, 'goalAdded') < 0.8 * sum(whole, 'goalAdded'), `one day of ${sum(close, 'goalAdded').toFixed(0)} against two of ${sum(whole, 'goalAdded').toFixed(0)}`);
    assert.deepEqual(whole.given.land.soilGoal, 'observed');
    assert.ok(!('soilGoal' in whole.options.land), 'the page never sees the goal');
  } finally {
    rmSync(out, { recursive: true, force: true });
  }
});

test("at N=16 spinup.mjs with OCEAN iceExchange turbulent and ICE iceAlbedoPhases, stopped inside a day, saves the water's heat under the ice and the conductance, and resumed with no OCEAN takes the exchange its save records and ends day 2 on the state of the run that did not stop, bit for bit but for the day's means", { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const out = mkdtempSync(join(tmpdir(), 'restart-exchange-')), iceExchange = { mode: 'turbulent', frictionMinimum: 2e-2 };
  try {
    const run = (extra) => execFileSync(process.execPath, [join(root, 'scripts/spinup.mjs')], { cwd: root, env: Object.fromEntries(Object.entries({ ...withoutOptionGroups(), N: '16', MINUTES: '100', OUT: out, DAYS: '2', OCEAN: JSON.stringify({ iceExchange }), ICE: JSON.stringify({ iceAlbedoPhases: true }), ...extra }).filter(([, value]) => value !== undefined)), encoding: 'utf8', stdio: ['ignore', 'pipe', 'pipe'] });
    run({ TAG: 'whole' });
    run({ TAG: 'halt', STOP_AFTER_STEPS: '96' });
    const inside = await decoded(join(out, 'halt_day0001_step0032.bin')), C = inside.surfaceT.length;
    assert.deepEqual(inside.given.ocean.iceExchange, iceExchange);
    for (const name of ['iceWater', 'iceConductance']) assert.equal(inside.ocean[name]?.length, C, `the save inside the day holds ocean.${name}`);
    const held = inside.ocean.iceWater.filter((x) => x > 0).length, stirred = inside.ocean.iceConductance.filter((x) => x > 0).length;
    assert.ok(held > 0 && stirred > 0, `the save inside the day holds the water's heat in ${held} cells and a conductance in ${stirred}`);
    assert.match(run({ TAG: 'halt', OCEAN: undefined }), /iceExchange \{"mode":"turbulent","frictionMinimum":0\.02\}, as halt_day0001_step0032\.bin records/);
    const [whole, halt] = await Promise.all(['whole', 'halt'].map((tag) => decoded(join(out, `${tag}_day0002.bin`))));
    assert.deepEqual(halt.given.ocean, { iceExchange }, 'the resumed run records the exchange it took');
    const means = ['convectiveRain', 'largeScaleRain', 'meanAbsorbedSolar', 'meanOutgoingLongwave', 'meanPlanetaryAlbedo', 'meanShortwaveCloudEffect', 'meanLongwaveCloudEffect'];
    const parted = apart(whole, halt, [...means, 'energyRecord']), gaps = gapRows(whole.energyRecord.length, [2]);
    assert.deepEqual(parted, [], `day 2, whole against halted: ${parted.join(', ')}`);
    const recorded = [...whole.energyRecord.keys()].filter((i) => (gaps[i] ? !Number.isNaN(halt.energyRecord[i]) : !Object.is(whole.energyRecord[i], halt.energyRecord[i])));
    assert.deepEqual(recorded, [], 'the energy record is the whole run\'s but for day 2, which the stop leaves unrecorded');
    assert.ok(whole.ocean.iceWater.some((x) => x > 0), 'the day-2 save holds the water\'s heat');
  } finally {
    rmSync(out, { recursive: true, force: true });
  }
});
