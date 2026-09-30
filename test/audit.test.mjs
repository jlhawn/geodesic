import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, writeFileSync, mkdtempSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { spawnSync } from 'node:child_process';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { inLongitudes, convectionLine, equatorialOcean, equatorLine } from '../js/audit.module.js';

const mesh = buildMesh(new Grid(16));
const C = mesh.nCells, E = mesh.nEdges;

test('a box runs east from its western longitude, across the date line when it has to', () => {
  for (const lon of [160, 180, -180, -150, -100]) assert.ok(inLongitudes(lon, 160, -100), `${lon} lies in 160E-100W`);
  for (const lon of [0, -90, 150]) assert.ok(!inLongitudes(lon, 160, -100), `${lon} lies outside 160E-100W`);
  for (const lon of [-180, -1, 0, 90, 179.9]) assert.ok(inLongitudes(lon, -180, 180), `${lon} lies in the whole circle`);
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
  const run = spawnSync(process.execPath, [new URL('../scripts/verticalAudit.mjs', import.meta.url).pathname, file], { env: { ...process.env, STEPS: '2' }, encoding: 'utf8' });
  rmSync(dir, { recursive: true, force: true });
  assert.equal(run.status, 0, run.stderr);
  const lines = run.stdout.split('\n').filter((line) => line.includes(' Earth '));
  assert.equal(lines.length, 26, `${lines.length} rows`);
  for (const name of ['SE Pacific 10-30S 110-80W: rain', 'Peru 5-20S 90-75W: columns firing a step', "deck's virtual jump above h", 'estimated inversion strength', 'saved 10-day deck sink', 'omega700 (Pa/s)', 'low cloud', 'deck height h', 'resolved inversion', 'Pacific ITCZ 5-12N 160E-100W: omega500', 'zonal-mean rain peak', 'omega700 grid-scale share', 'SH Hadley peak', 'NH Hadley peak']) {
    assert.ok(lines.some((line) => line.includes(name)), `a row for ${name}`);
  }
  for (const line of lines) assert.match(line, /Earth (-?[0-9.]+\.\.-?[0-9.]+|n\/a\.\.n\/a)  -> (matches|too weak by x[0-9.]+|too strong by x[0-9.]+|wrong sign|too noisy to tell|none|n\/a)/);
  assert.match(run.stdout, /SH Hadley peak \(1e9 kg\/s\) +-?[0-9]+\.[0-9]/);
});
