import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, writeFileSync, mkdtempSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { spawnSync } from 'node:child_process';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { inLongitudes, convectionLine, equatorialOcean, equatorLine, heatingProfile, tropicalBoxOf, TROPICAL_BOXES } from '../js/audit.module.js';

const mesh = buildMesh(new Grid(16));
const C = mesh.nCells, E = mesh.nEdges;

test('a box runs east from its western longitude, across the date line when it has to', () => {
  for (const lon of [160, 180, -180, -150, -100]) assert.ok(inLongitudes(lon, 160, -100), `${lon} lies in 160E-100W`);
  for (const lon of [0, -90, 150]) assert.ok(!inLongitudes(lon, 160, -100), `${lon} lies outside 160E-100W`);
  for (const lon of [-180, -1, 0, 90, 179.9]) assert.ok(inLongitudes(lon, -180, 180), `${lon} lies in the whole circle`);
});

test('a heating profile peaks in the 50 hPa bin of largest mass-weighted mean, not at the largest thin layer, and its centroid weighs only the heating', () => {
  const p = [150e2, 420e2, 440e2, 610e2, 960e2, 970e2], dp = [50e2, 20e2, 20e2, 100e2, 10e2, 10e2], Q = [5, 3, 2, 1, 6, -4];
  const peak = heatingProfile(Q, p, dp, { top: 200e2 });
  assert.equal(peak.layer, 960e2);
  assert.equal(peak.value, 6);
  assert.deepEqual(peak.bin, [400e2, 450e2]);
  assert.equal(peak.binValue, (3 * 20e2 + 2 * 20e2) / 40e2);
  assert.equal(peak.centroid, (420e2 * 3 * 20e2 + 440e2 * 2 * 20e2 + 610e2 * 1 * 100e2 + 960e2 * 6 * 10e2) / (3 * 20e2 + 2 * 20e2 + 100e2 + 6 * 10e2));
  assert.ok(Number.isNaN(heatingProfile([-1, -2], [500e2, 600e2], [1, 1]).centroid), 'no heating, no centroid');
});

test('each tropical cell counts in the first box that holds it, the sea boxes keep sea and the Amazon land', () => {
  const land = Uint8Array.from({ length: C }, (_, i) => (mesh.lonCell[i] * 180 / Math.PI < -55 ? 1 : 0)), boxOf = tropicalBoxOf(mesh, land);
  const deg = 180 / Math.PI;
  for (let i = 0; i < C; i++) {
    const lat = mesh.latCell[i] * deg, lon = mesh.lonCell[i] * deg;
    const first = TROPICAL_BOXES.findIndex(([, [s, n, w, e], surface]) => lat >= s && lat <= n && inLongitudes(lon, w, e) && (surface === 'all' || (surface === 'land') === !!land[i]));
    assert.equal(boxOf[i], first);
    if (boxOf[i] >= 0) assert.ok(TROPICAL_BOXES[boxOf[i]][2] === 'all' || (TROPICAL_BOXES[boxOf[i]][2] === 'land') === !!land[i]);
  }
  assert.ok([0, 1, 4].every((b) => boxOf.includes(b)), 'the ITCZ, the warm pool and the Amazon hold cells');
});

test('the convection line gives the convective share, the SE Pacific rain and how often it rains convectively there, and the Pacific ITCZ rain', () => {
  const split = { convective: new Float64Array(C).fill(3), largeScale: new Float64Array(C).fill(1), wet: new Float64Array(C).fill(0.5) };
  const line = convectionLine(mesh, { land: new Uint8Array(C), ice: new Float64Array(C) }, split, 30);
  assert.equal(line, 'convection after 30 days: convective share global 0.75, 15S-15N 0.75; SE Pacific 10-30S 110-80W sea 4.00 mm/d (convective 3.00), convective rain on 0.50 of its column-days; Pacific ITCZ 5-12N 160E-100W 4.00 mm/d');
});

test('the equator line reads the mixed layer\'s current, the fastest eastward class, the stress, the mixed-layer depth and the thermocline from a layered ocean', () => {
  const eastward = (speed) => Float64Array.from({ length: E }, (_, e) => {
    const x = mesh.xEdge[3 * e], y = mesh.xEdge[3 * e + 1], r = Math.hypot(x, y);
    return speed * (-y * mesh.nEdge[3 * e] + x * mesh.nEdge[3 * e + 1]) / r;
  });
  const speeds = [0.3, 0.8, -0.1], thickness = [50, 100, 200];
  const ocean = {
    h: Float64Array.from({ length: 3 * C }, (_, n) => thickness[Math.floor(n / C)]),
    u: Float64Array.from(speeds.flatMap((speed) => Array.from(eastward(speed)))),
    densities: [1023, 1025],
  };
  const o = equatorialOcean(mesh, new Uint8Array(C), ocean, eastward(-0.05));
  const near = (value, expected, what) => assert.ok(Math.abs(value - expected) < 0.03 * Math.abs(expected), `${what} ${value}, not ${expected}`);
  near(o.surface, 0.3, 'surface current');
  near(o.surfaceEast, 0.3, 'eastern surface current');
  near(o.undercurrent, 0.8, 'undercurrent');
  near(o.stress, -0.05, 'stress');
  assert.equal(o.undercurrentClass, 1023);
  assert.ok(Math.abs(o.undercurrentDepth - 100) < 1e-9 && Math.abs(o.mixedEast - 50) < 1e-9, `undercurrent at ${o.undercurrentDepth} m, mixed layer ${o.mixedEast} m`);
  assert.ok(Math.abs(o.thermoclineWest - 150) < 1e-9 && Math.abs(o.thermoclineEast - 150) < 1e-9, `the 1024 class top at ${o.thermoclineWest} and ${o.thermoclineEast} m`);
  assert.match(equatorLine(mesh, new Uint8Array(C), ocean, eastward(-0.05), 30), /^equator after 30 days \(2S-2N, eastward \+\): surface current 160E-100W \+0\.3\d m\/s, 140W-100W \+0\.3\d m\/s; undercurrent \+0\.8\d m\/s at 100 m \(class 1023, 180-100W\); stress 160E-100W -0\.05\d N\/m²; mixed layer 140W-100W 50 m; 1024 class top 150E-180 150 m, 120W-90W 150 m$/);
});

test('scripts/verticalAudit.mjs prints every headline number of a saved state with its Earth reference and a verdict', async () => {
  const { createModel } = await import('../js/model.module.js');
  const { topographyFromInt16 } = await import('../js/geography.module.js');
  const { initializeState } = await import('../js/physics/init.module.js');
  const { encodeState } = await import('../js/stateFile.module.js');
  const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
  const model = createModel(new Grid(12), { topography, ocean: false });
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  model.land.initialize();
  const [pi, theta, u, surfaceT, q, qc, ice] = model.state, { mlmSubsidence, mlmHeight, mlmGate } = model.radiation;
  const dir = mkdtempSync(join(tmpdir(), 'verticalAudit-')), file = join(dir, 'audit12_day0000.bin');
  writeFileSync(file, encodeState({ N: 12, K: model.core.K, day: 0, time: 0, terrain: true, levels: model.core.levels, pi, theta, u, surfaceT, q, qc, ice, concentration: model.seaIce.concentration, mlmSubsidence, mlmHeight, mlmGate, land: model.land.serialize() }));
  const withEffects = join(dir, 'audit12_day0001.bin'), C = model.mesh.nCells;
  writeFileSync(withEffects, encodeState({ N: 12, K: model.core.K, day: 1, time: 0, terrain: true, levels: model.core.levels, pi, theta, u, surfaceT, q, qc, ice, concentration: model.seaIce.concentration, mlmSubsidence, mlmHeight, mlmGate, land: model.land.serialize(), meanShortwaveCloudEffect: new Float64Array(C).fill(-40), meanLongwaveCloudEffect: new Float64Array(C).fill(25) }));
  const audit = (state) => spawnSync(process.execPath, [new URL('../scripts/verticalAudit.mjs', import.meta.url).pathname, state], { env: { ...process.env, STEPS: '2' }, encoding: 'utf8' });
  const run = audit(file), saved = audit(withEffects);
  rmSync(dir, { recursive: true, force: true });
  assert.equal(run.status, 0, run.stderr);
  assert.equal(saved.status, 0, saved.stderr);
  const lines = run.stdout.split('\n').filter((line) => line.includes(' Earth '));
  const classLines = lines.filter((line) => line.includes('clear sky, '));
  assert.equal(lines.length - classLines.length, 86, `${lines.length - classLines.length} rows besides the surface classes`);
  assert.ok(classLines.some((line) => /clear sky, open sea 0-30: albedo at the top +0\.[0-9]{3}  Earth 0\.080\.\.0\.100  -> /.test(line)), 'a row for the tropical open sea at the top');
  assert.ok(classLines.some((line) => /clear sky, partly vegetated \(0\.2-0\.7\): surface albedo +0\.[0-9]{3}  Earth 0\.180\.\.0\.250  -> /.test(line)), 'a row for partly vegetated land');
  assert.equal(classLines.filter((line) => line.includes('direct-beam surface albedo')).length, 5, 'five rows of the open sea by the sun\'s cosine');
  assert.match(run.stdout, /global clear-sky albedo, mean over the window's 2 steps +0\.[0-9]{3}  Earth n\/a\.\.n\/a  -> n\/a  \[the window's 2 steps 0\.[0-9]{3}; an outcome of the surface classes below, Earth about 0\.15 globally\]/);
  assert.match(run.stdout, /30S-30N clear-sky albedo, mean over the window's 2 steps +0\.[0-9]{3}  Earth n\/a\.\.n\/a  -> n\/a  \[the window's 2 steps 0\.[0-9]{3}; an outcome/);
  for (const name of ['global shortwave cloud effect', 'global longwave cloud effect', '30S-30N shortwave cloud effect', '30S-30N longwave cloud effect']) assert.ok(lines.some((line) => line.includes(`${name}, mean over the window's 2 steps (W/m2)`)), `a row for the window's ${name}`);
  assert.match(saved.stdout, /global shortwave cloud effect, day mean of the state's last day \(W\/m2\) +-40\.0  Earth -43\.0\.\.-51\.0  -> too weak by x1\.08  \[the window's 2 steps -?[0-9.]+\]/);
  assert.match(saved.stdout, /30S-30N longwave cloud effect, day mean of the state's last day \(W\/m2\) +25\.0  Earth n\/a\.\.n\/a  -> n\/a  \[the window's 2 steps -?[0-9.]+; Earth \+26 \+- 3 globally\]/);
  for (const name of ['SE Pacific 10-30S 110-80W: rain', 'Peru 5-20S 90-75W: columns firing a step', "deck's virtual jump above h", 'estimated inversion strength', 'saved running-mean deck sink', 'deck-height sink now, as the dynamics leaves it', 'deck-height sink now, as the deck reads it, 2 ring passes', 'omega700 (Pa/s)', 'low cloud', "deck's start height h", 'deck runs, share of column-steps', 'deck height where it runs', "deck's cloud-layer thickness where it runs", "deck's liquid water path where it runs", 'global evaporation (mm/d)', 'resolved inversion (m)', "resolved inversion's thetaV jump", 'Pacific ITCZ 5-12N 160E-100W: omega500', 'zonal-mean rain peak', 'omega700 grid-scale share', 'SH Hadley peak', 'NH Hadley peak', 'Pacific ITCZ 5-12N 160E-100W: Q1-QR peak over 50 hPa bins (hPa)', 'Pacific ITCZ 5-12N 160E-100W: Q1-QR centroid (hPa)', 'warm pool 10S-10N 120-170E sea: Q1-QR peak over 50 hPa bins (hPa)', 'warm pool 10S-10N 120-170E sea: Q1-QR centroid (hPa)', 'large-scale heating below 1 km', 'global rain (mm/d)', 'convective share of the rain, 15S-15N', 'low-cloud cover, radiative', 'low-cloud water path in cloud', 'coupled stratocumulus share', 'cloud-top cooling', 'Namibia 10-20S 0-10E: rain', 'California 20-30N 130-120W: low-cloud cover']) {
    assert.ok(lines.some((line) => line.includes(name)), `a row for ${name}`);
  }
  for (const line of lines) assert.match(line, /Earth (-?[0-9.]+\.\.-?[0-9.]+|n\/a\.\.n\/a)  -> (matches|too weak by x[0-9.]+|too strong by x[0-9.]+|wrong sign|too noisy to tell|none|n\/a)/);
  assert.match(run.stdout, /SH Hadley peak \(1e9 kg\/s\) +-?[0-9]+\.[0-9]/);
  assert.match(run.stdout, /replica over [1-9][0-9]* column-steps: running mean \|diff\| 0\.0e\+0 m\/s, gate \|diff\| 0\.0e\+0, run decisions differing 0/);
  assert.match(run.stdout, /off: subsidence [0-9.]+, off: jump [0-9.]+, off: gate memory [0-9.]+, stood down [0-9.]+; failing the subsidence test [0-9.]+, the jump test [0-9.]+/);
});

test('scripts/tropicalHeating.mjs replays the moist step exactly, closes the heat and total-water budgets, matches the model\'s sensible heat and prints a summary of every tropical box', async () => {
  const { createModel } = await import('../js/model.module.js');
  const { topographyFromInt16 } = await import('../js/geography.module.js');
  const { initializeState } = await import('../js/physics/init.module.js');
  const { encodeState } = await import('../js/stateFile.module.js');
  const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
  const model = createModel(new Grid(12), { topography, ocean: false });
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  model.land.initialize();
  for (let n = 0; n < 6; n++) model.step(1800);
  const [pi, theta, u, surfaceT, q, qc, ice] = model.state, { mlmSubsidence, mlmHeight, mlmGate } = model.radiation;
  const dir = mkdtempSync(join(tmpdir(), 'tropicalHeating-')), file = join(dir, 'heat12_day0000.bin');
  writeFileSync(file, encodeState({ N: 12, K: model.core.K, day: 0, time: model.time, terrain: true, levels: model.core.levels, pi, theta, u, surfaceT, q, qc, ice, concentration: model.seaIce.concentration, mlmSubsidence, mlmHeight, mlmGate, land: model.land.serialize() }));
  const run = spawnSync(process.execPath, [new URL('../scripts/tropicalHeating.mjs', import.meta.url).pathname, file], { env: { ...process.env, STEPS: '3' }, encoding: 'utf8' });
  rmSync(dir, { recursive: true, force: true });
  assert.equal(run.status, 0, run.stderr);
  const [, mismatch, columns, heat, water, sensible, own] = run.stdout.match(/differs from the model's in (\d+) of (\d+) column-steps.*close on the temperature change to ([0-9.e+-]+) K a step and on the change of q_t to ([0-9.e+-]+) kg\/kg a step.*global sensible heat (-?[0-9.]+) against the model's (-?[0-9.]+) W\/m2/);
  assert.equal(Number(mismatch), 0, `${mismatch} of ${columns} column-steps`);
  assert.ok(Number(columns) > 0 && Number(heat) < 1e-12 && Number(water) < 1e-16, `heat ${heat} K, water ${water} kg/kg a step`);
  assert.equal(sensible, own, 'the global sensible heat');
  const summaries = run.stdout.split('\n').filter((line) => line.startsWith('summary '));
  assert.deepEqual(summaries.map((line) => line.slice(8, line.indexOf(': {'))), TROPICAL_BOXES.map(([name]) => name));
  for (const line of summaries) assert.ok('q1rCentroid' in JSON.parse(line.slice(line.indexOf(': {') + 2)), line);
});
