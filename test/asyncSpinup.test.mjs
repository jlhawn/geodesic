import { test, before, after } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, writeFileSync, mkdtempSync, mkdirSync, rmSync, readdirSync, copyFileSync, existsSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { spawn, spawnSync } from 'node:child_process';
import { Grid } from '../js/grid.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { decodeState, encodeState } from '../js/stateFile.module.js';
import { decodeForcing, forcingName } from '../js/forcing.module.js';
import { withOceanOf } from '../js/oceanHandOff.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};
const skip = !gpuAvailable && 'webgpu not installed';

const root = new URL('..', import.meta.url).pathname;
const N = 6, PER_DAY = 24;
let dir, fresh, forcing;

function run(script, env) {
  const result = spawnSync(process.execPath, [join(root, 'scripts', script)], { cwd: root, env: { ...process.env, N: String(N), ...env }, encoding: 'utf8' });
  if (result.status !== 0) console.log(`${script} exited ${result.status}: ${result.stderr.slice(-2000)}`);
  return result.status;
}
function driver(env) {
  const result = spawnSync('bash', [join(root, 'scripts/asyncSpinup.sh')], { cwd: root, env: { ...process.env, N: String(N), YEAR_DAYS: '2', COUPLED_YEARS: '1', PER_YEAR: '1', SNAPSHOT_DAYS: '0', RETRY_WAIT: '1', ...env }, encoding: 'utf8' });
  if (result.status !== 0) console.log(`asyncSpinup.sh exited ${result.status}: ${result.stderr.slice(-2000)}`);
  return result.status;
}
const load = async (file) => decodeState(new Uint8Array(readFileSync(file)));
const text = (file) => readFileSync(file, 'utf8');
async function waitFor(check, seconds = 60) {
  const until = Date.now() + seconds * 1000;
  while (!check()) {
    if (Date.now() > until) throw new Error('timed out');
    await new Promise((resolve) => setTimeout(resolve, 5));
  }
}
const exited = (child) => new Promise((resolve) => child.on('exit', (code, signal) => resolve(code ?? signal)));
function outDir(name, withFresh = null) {
  const out = join(dir, name);
  mkdirSync(out);
  if (withFresh) copyFileSync(fresh, join(out, withFresh));
  return out;
}

before(async () => {
  dir = mkdtempSync(join(tmpdir(), 'asyncSpinup-'));
  if (!gpuAvailable) return;
  const topography = topographyFromInt16(readFileSync(join(root, 'data/topography_0p25.bin')).buffer);
  const model = await createGpuModel(new Grid(N), { topography });
  const { state, mesh } = model;
  initializeState(model, { geostrophic: !model.surfaceGeopotential }).forEach((values, a) => state[a].set(values));
  for (let i = 0; i < mesh.nCells; i++) if (model.geography.land[i]) state[6][i] = 0;
  model.load();
  model.ocean.initialize(state[3], state[6]);
  model.land.initialize();
  const ocean = await model.ocean.serialize(), land = await model.land.serialize();
  const [pi, theta, u, surfaceT, q, qc, ice] = state;
  fresh = join(dir, 'fresh_day0000.bin');
  writeFileSync(fresh, encodeState({ N, K: model.core.K, day: 0, time: 0, terrain: !!model.surfaceGeopotential, pi, theta, u, surfaceT, q, qc, ice, concentration: model.seaIce.concentration, mlmSubsidence: model.radiation.mlmSubsidence, ocean: { h: ocean.h, u: ocean.u, T: ocean.T, S: ocean.S, eta: ocean.eta }, land }));
  model.destroy();
  const recorded = outDir('recorded', 'rec_day0000.bin');
  forcing = join(recorded, 'forcing');
  assert.equal(run('spinup.mjs', { OUT: recorded, TAG: 'rec', RECORD: forcing, DAYS: '2', MINUTES: '100' }), 0);
});
after(() => rmSync(dir, { recursive: true, force: true }));

test('the hand-off takes the sea surface, ice and ocean and keeps the atmosphere, the land and the clock', () => {
  const C = 5, land = Uint8Array.from([1, 0, 0, 1, 0]);
  const cells = (base) => Float32Array.from({ length: C }, (_, i) => base + i);
  const coupled = { N: 1, day: 730, time: 730 * 86400, oceanYears: 100, levels: Float64Array.from([0, 0.5, 1]), pi: cells(1e5), surfaceT: cells(280), ice: cells(0), concentration: cells(0.1), mlmSubsidence: cells(0.01), mlmHeight: cells(800), mlmGate: cells(0.2), ocean: { h: cells(1), T: cells(275) }, land: { soil: cells(10), snow: cells(1) } };
  const alone = { N: 1, day: 37230, time: 37230 * 86400, oceanYears: 100, levels: Float64Array.from([0, 0.4, 1]), pi: cells(9e4), surfaceT: cells(270), ice: cells(2), concentration: cells(0.5), mlmSubsidence: cells(0.02), mlmHeight: cells(900), mlmGate: cells(0.7), ocean: { h: cells(5), T: cells(271), Q: cells(7) }, land: { soil: cells(20), snow: cells(3) } };
  const merged = withOceanOf(coupled, alone, land);
  for (let i = 0; i < C; i++) {
    const from = land[i] ? coupled : alone;
    assert.equal(merged.surfaceT[i], from.surfaceT[i]);
    assert.equal(merged.ice[i], from.ice[i]);
    assert.equal(merged.concentration[i], from.concentration[i]);
    assert.equal(merged.land.snow[i], from.land.snow[i]);
  }
  assert.equal(merged.ocean.h, alone.ocean.h);
  assert.equal(merged.ocean.Q, alone.ocean.Q);
  assert.equal(merged.pi, coupled.pi);
  for (const field of ['levels', 'mlmSubsidence', 'mlmHeight', 'mlmGate']) assert.equal(merged[field], coupled[field], field);
  assert.equal(merged.land.soil, coupled.land.soil);
  assert.deepEqual([merged.day, merged.time, merged.oceanYears], [730, 730 * 86400, 200]);
  assert.throws(() => withOceanOf(coupled, { ...alone, N: 2 }, land), /N=2/);
});

test('an ocean-only run stopped inside a day and continued ends byte for byte where the uninterrupted run does', { skip }, async () => {
  const env = { TAG: 'alone', STATE: join(dir, 'recorded/rec_day0002.bin'), FORCING: forcing, YEARS: '3', DAYS_PER_YEAR: '2', SNAPSHOT_DAYS: '1' };
  const whole = outDir('whole'), parts = outDir('parts');
  assert.equal(run('oceanSpinup.mjs', { ...env, OUT: whole }), 0);
  assert.equal(run('oceanSpinup.mjs', { ...env, OUT: parts, STOP_AFTER_STEPS: String(PER_DAY + 4) }), 0);
  assert.deepEqual(readdirSync(parts).filter((f) => f.endsWith('.bin')), ['alone_year0000_day001_step0004.bin']);
  const inside = await load(join(parts, 'alone_year0000_day001_step0004.bin'));
  assert.deepEqual([inside.oceanYears, inside.oceanDays, inside.oceanStep, inside.day], [0, 1, 4, 3]);
  assert.ok(inside.ocean.Q && inside.ocean.flux, 'the in-day file carries the restart arrays');
  for (const field of ['mlmSubsidence', 'mlmHeight', 'mlmGate']) assert.ok(inside[field], `the in-day file carries ${field}`);
  assert.equal(run('oceanSpinup.mjs', { ...env, OUT: parts, STOP_AFTER_STEPS: String(PER_DAY + 16) }), 0);
  assert.deepEqual(readdirSync(parts).filter((f) => f.endsWith('.bin')).sort(), ['alone_year0001.bin', 'alone_year0001_day000_step0020.bin']);
  assert.equal(run('oceanSpinup.mjs', { ...env, OUT: parts }), 0);
  assert.deepEqual(readdirSync(parts).filter((f) => f.endsWith('.bin')).sort(), ['alone_year0002.bin', 'alone_year0003.bin']);
  assert.ok(readFileSync(join(whole, 'alone_year0003.bin')).equals(readFileSync(join(parts, 'alone_year0003.bin'))), 'the continued run ends where the uninterrupted one does');
  const log = text(join(parts, 'alone.log'));
  assert.match(log, /stopped by STOP_AFTER_STEPS=28 after 0 years 1 days 4 of 24 steps alone/);
  assert.match(log, /year 3: .*global SST -?[0-9.]+ °C, drift -?[0-9.]+ K since the start, -?[0-9.]+ K over the last 3 years;/);
});

test('SIGTERM makes an ocean-only run save and exit 0, and the run continues to the same end', { skip }, async () => {
  const env = { TAG: 'alone', STATE: join(dir, 'recorded/rec_day0002.bin'), FORCING: forcing, YEARS: '20', DAYS_PER_YEAR: '2', SNAPSHOT_DAYS: '0' };
  const whole = outDir('signalWhole'), parts = outDir('signalParts');
  assert.equal(run('oceanSpinup.mjs', { ...env, OUT: whole }), 0);
  const child = spawn(process.execPath, [join(root, 'scripts/oceanSpinup.mjs')], { cwd: root, env: { ...process.env, N: String(N), ...env, OUT: parts }, stdio: 'ignore' });
  await waitFor(() => existsSync(join(parts, 'alone.log')) && text(join(parts, 'alone.log')).includes('saved alone_year0002.bin'));
  const sent = Date.now();
  child.kill('SIGTERM');
  assert.equal(await exited(child), 0);
  assert.ok(Date.now() - sent < 20000, `stopped ${Date.now() - sent} ms after SIGTERM`);
  const log = text(join(parts, 'alone.log'));
  assert.match(log, /SIGTERM at .*: saving after the ocean step in progress/);
  assert.match(log, /stopped by SIGTERM after \d+ years/);
  assert.ok(!existsSync(join(parts, 'alone_year0020.bin')), 'the run stopped before its end');
  assert.equal(run('oceanSpinup.mjs', { ...env, OUT: parts }), 0);
  assert.ok(readFileSync(join(whole, 'alone_year0020.bin')).equals(readFileSync(join(parts, 'alone_year0020.bin'))));
});

test('a coupled segment stopped inside a day carries the recorded day across the interruption', { skip }, async () => {
  const whole = outDir('coupledWhole', 'c_day0000.bin'), parts = outDir('coupledParts', 'c_day0000.bin');
  const env = { TAG: 'c', DAYS: '3', MINUTES: '100' };
  assert.equal(run('spinup.mjs', { ...env, OUT: whole, RECORD: join(whole, 'forcing') }), 0);
  assert.equal(run('spinup.mjs', { ...env, OUT: parts, RECORD: join(parts, 'forcing'), STOP_AFTER_STEPS: String(PER_DAY + 12) }), 0);
  assert.deepEqual(readdirSync(parts).filter((f) => f.endsWith('.bin')).sort(), ['c_day0000.bin', 'c_day0001_step0012.bin']);
  const inside = await load(join(parts, 'c_day0001_step0012.bin'));
  assert.deepEqual([inside.day, inside.step, inside.forcingSteps, inside.forcingOceanSteps, inside.forcingSeconds], [1, 12, 12, 3, 12 * 3600]);
  for (const field of ['mlmSubsidence', 'mlmHeight', 'mlmGate', 'boundaryDepth', 'mixingTop', 'boundaryRegime', 'boundaryBuoyancy']) assert.ok(inside[field], `the in-day checkpoint carries ${field}`);
  assert.match(text(join(parts, 'c.log')), /stopped by STOP_AFTER_STEPS=36 at day 1 and 12 of 24 steps/);
  assert.equal(run('spinup.mjs', { ...env, OUT: parts, RECORD: join(parts, 'forcing') }), 0);
  assert.deepEqual(readdirSync(parts).filter((f) => f.endsWith('.bin')).sort(), ['c_day0000.bin', 'c_day0003.bin']);
  assert.ok(readFileSync(join(whole, 'forcing', forcingName(1))).equals(readFileSync(join(parts, 'forcing', forcingName(1)))));
  for (const day of [2, 3]) {
    const [a, b] = await Promise.all([whole, parts].map(async (out) => decodeForcing(new Uint8Array(readFileSync(join(out, 'forcing', forcingName(day)))))));
    assert.deepEqual([b.seconds, b.steps, b.oceanSteps], [a.seconds, a.steps, a.oceanSteps]);
    for (const name of ['stress', 'netFlux', 'shortwave', 'shortwaveDown', 'evaporation', 'rain']) {
      let off = 0, size = 0;
      for (let i = 0; i < a.fields[name].length; i++) { off += Math.abs(a.fields[name][i] - b.fields[name][i]); size += Math.abs(a.fields[name][i]); }
      assert.ok(off < 0.05 * size, `day ${day} ${name}: the continued run's mean is ${(100 * off / size).toFixed(1)}% off the uninterrupted one's`);
    }
  }
});

test('a coupled run split at each day\'s end ends byte for byte where the uninterrupted run does', { skip }, async () => {
  const whole = outDir('splitWhole', 's_day0000.bin'), parts = outDir('splitParts', 's_day0000.bin');
  const env = { TAG: 's', MINUTES: '100' };
  assert.equal(run('spinup.mjs', { ...env, OUT: whole, DAYS: '3' }), 0);
  assert.equal(run('spinup.mjs', { ...env, OUT: parts, DAYS: '1' }), 0);
  const day = await load(join(parts, 's_day0001.bin'));
  for (const field of ['windSpeed', 'evaporation', 'cumulusCover', 'cumulusWater']) assert.ok(day[field], `the day's file carries ${field}`);
  assert.ok(day.ocean.Q && day.ocean.flux && day.ocean.capacity, 'the day\'s file carries the ocean\'s restart arrays');
  assert.equal(run('spinup.mjs', { ...env, OUT: parts, DAYS: '2' }), 0);
  assert.equal(run('spinup.mjs', { ...env, OUT: parts, DAYS: '3' }), 0);
  assert.ok(readFileSync(join(whole, 's_day0003.bin')).equals(readFileSync(join(parts, 's_day0003.bin'))), 'the split run ends where the uninterrupted one does');
});

test('a coupled run split inside a day ends the day with the uninterrupted run\'s state', { skip }, async () => {
  const whole = outDir('halfWhole', 'h_day0000.bin'), parts = outDir('halfParts', 'h_day0000.bin');
  const env = { TAG: 'h', MINUTES: '100', DAYS: '1' };
  assert.equal(run('spinup.mjs', { ...env, OUT: whole }), 0);
  assert.equal(run('spinup.mjs', { ...env, OUT: parts, STOP_AFTER_STEPS: '12' }), 0);
  const inside = await load(join(parts, 'h_day0000_step0012.bin'));
  for (const field of ['rainTotal', 'runoffTotal', 'rainSeen', 'runoffSeen']) assert.ok(inside[field], `the in-day file carries ${field}`);
  assert.equal(run('spinup.mjs', { ...env, OUT: parts }), 0);
  const [a, b] = await Promise.all([whole, parts].map((out) => load(join(out, 'h_day0001.bin'))));
  const means = ['convectiveRain', 'largeScaleRain', 'meanAbsorbedSolar', 'meanOutgoingLongwave', 'meanPlanetaryAlbedo', 'meanShortwaveCloudEffect', 'meanLongwaveCloudEffect'];
  const apart = [];
  for (const [name, x] of [...Object.entries(a), ...Object.entries(a.ocean).map(([k, v]) => [`ocean.${k}`, v]), ...Object.entries(a.land).map(([k, v]) => [`land.${k}`, v])]) {
    if (!ArrayBuffer.isView(x) || means.includes(name)) continue;
    const y = name.startsWith('ocean.') ? b.ocean[name.slice(6)] : name.startsWith('land.') ? b.land[name.slice(5)] : b[name];
    if (!x.every((v, n) => Object.is(v, y[n]))) apart.push(name);
  }
  assert.deepEqual(apart, []);
});

test('the asynchronous driver alternates coupled and ocean-only phases, hands the ocean over and picks up where it stands', { skip }, async () => {
  const out = outDir('driver', 'async6_day0000.bin');
  assert.equal(driver({ OUT: out, CYCLES: '2', OCEAN_YEARS: '2' }), 0);
  assert.deepEqual(readdirSync(out).filter((f) => f.endsWith('.bin')).sort(), ['async6_c01_year0002.bin', 'async6_c02_year0002.bin', 'async6_day0000.bin', 'async6_day0002.bin', 'async6_day0004.bin']);
  assert.ok(!existsSync(join(out, 'async6_c01_forcing')) && !existsSync(join(out, 'async6_c02_forcing')), 'the records are deleted once used');
  assert.equal(text(join(out, 'async6_cycles.txt')), 'start=0\ndone=1\ndone=2\n');
  const summaries = text(join(out, 'async6_cycles.log')).split('\n').filter((line) => / cycle \d+: /.test(line));
  assert.equal(summaries.length, 2);
  assert.match(summaries[1], /cycle 2: 2 coupled years to day 4, 4 ocean-only years; warm pool .* °C, cold tongue .*; Southern Ocean 60–70S 0–60 m .*; global SST .* °C, drift .* K since the start, .* K over the last 2 years; ice extent N .* Mkm²/);
  const [handed, alone, second] = await Promise.all(['async6_day0004.bin', 'async6_c01_year0002.bin', 'async6_c02_year0002.bin'].map((f) => load(join(out, f))));
  assert.deepEqual([handed.day, handed.oceanFrom, handed.oceanYears], [4, 'async6_c01_year0002.bin', 2]);
  assert.deepEqual([alone.day, alone.oceanYears, second.day], [6, 2, 8]);
  assert.match(text(join(out, 'async6.log')), /ocean, sea ice and sea surface of async6_day0002\.bin replaced by those of .*async6_c01_year0002\.bin/);

  assert.equal(driver({ OUT: out, CYCLES: '2', OCEAN_YEARS: '2' }), 0);
  assert.equal(readdirSync(out).filter((f) => f.endsWith('.bin')).length, 5, 'a finished schedule runs nothing more');
  assert.equal(run('spinup.mjs', { OUT: out, TAG: 'async6', DAYS: '5', MINUTES: '100', OCEAN_FROM: join(out, 'async6_c01_year0002.bin') }), 0);
  assert.match(text(join(out, 'async6.log')), /async6_day0004\.bin already carries the ocean of .*async6_c01_year0002\.bin/);

  writeFileSync(join(out, 'STOP_async'), '');
  assert.equal(driver({ OUT: out, CYCLES: '3', OCEAN_YEARS: '2' }), 0);
  assert.match(text(join(out, 'async6_cycles.log')), /stopped at STOP_async or STOP_async6 at coupled day 5/);
  assert.ok(!existsSync(join(out, 'async6_c03_forcing')));
});

test('SIGTERM to the driver reaches its running phase, and a second start finishes as an uninterrupted one does', { skip }, async () => {
  const env = { N: String(N), YEAR_DAYS: '2', COUPLED_YEARS: '1', PER_YEAR: '1', SNAPSHOT_DAYS: '0', RETRY_WAIT: '1', CYCLES: '1', OCEAN_YEARS: '20' };
  const whole = outDir('driverWhole', 'async6_day0000.bin'), parts = outDir('driverParts', 'async6_day0000.bin');
  assert.equal(driver({ ...env, OUT: whole }), 0);
  const child = spawn('bash', [join(root, 'scripts/asyncSpinup.sh')], { cwd: root, env: { ...process.env, ...env, OUT: parts }, stdio: 'ignore' });
  await waitFor(() => existsSync(join(parts, 'async6_c01.log')) && text(join(parts, 'async6_c01.log')).includes('saved async6_c01_year0002.bin'));
  assert.equal(driver({ ...env, OUT: parts }), 0);
  assert.match(text(join(parts, 'async6_cycles.log')), new RegExp(`another driver \\(pid ${child.pid}\\) holds .*async6_cycles\\.lock; not starting`));
  child.kill('SIGTERM');
  assert.equal(await exited(child), 0);
  assert.match(text(join(parts, 'async6_cycles.log')), /stopped by a signal at coupled day 2/);
  assert.match(text(join(parts, 'async6_c01.log')), /stopped by SIGTERM after \d+ years/);
  assert.ok(!existsSync(join(parts, 'async6_c01_year0020.bin')));
  assert.equal(driver({ ...env, OUT: parts }), 0);
  assert.ok(readFileSync(join(whole, 'async6_c01_year0020.bin')).equals(readFileSync(join(parts, 'async6_c01_year0020.bin'))));
  assert.equal(text(join(parts, 'async6_cycles.txt')), 'start=0\ndone=1\n');
  assert.ok(!existsSync(join(parts, 'async6_cycles.lock')), 'the lock goes with the driver');
});

test('SYNC_CMD receives every file the phases save, and RESTORE_CMD lets another machine carry the schedule on', { skip }, async () => {
  const first = outDir('syncFirst', 'async6_day0000.bin'), second = outDir('syncSecond'), bucket = outDir('bucket');
  const env = { OCEAN_YEARS: '2', SNAPSHOT_DAYS: '1', BUCKET: bucket, SYNC_CMD: 'mkdir -p "$BUCKET/$(dirname "${1#$OUT/}")" && cp "$1" "$BUCKET/${1#$OUT/}"' };
  assert.equal(driver({ ...env, OUT: first, CYCLES: '1' }), 0);
  const held = spawnSync('find', [bucket, '-type', 'f'], { encoding: 'utf8' }).stdout.split('\n').filter(Boolean).map((f) => f.slice(bucket.length + 1)).sort();
  for (const name of ['async6_day0002.bin', 'async6_c01_forcing/forcing-0001.bin', 'async6_c01_forcing/forcing-0002.bin', 'async6_c01_year0000_day001.bin', 'async6_c01_year0002.bin', 'async6.log', 'async6_c01.log', 'async6_cycles.log', 'async6_cycles.txt']) assert.ok(held.includes(name), `${name} reached the bucket: ${held.join(', ')}`);
  assert.equal(driver({ ...env, OUT: second, CYCLES: '2', RESTORE_CMD: 'cp -R "$BUCKET/." "$OUT/"' }), 0);
  for (const log of ['async6.log', 'async6_c01.log', 'async6_c02.log']) assert.doesNotMatch(text(join(second, log)), /SYNC_CMD failed/, log);
  assert.equal(text(join(second, 'async6_cycles.txt')), 'start=0\ndone=1\ndone=2\n');
  const handed = await load(join(second, 'async6_day0004.bin'));
  assert.deepEqual([handed.oceanFrom, handed.oceanYears], ['async6_c01_year0002.bin', 2]);
  assert.ok(existsSync(join(second, 'async6_c02_year0002.bin')));
});
