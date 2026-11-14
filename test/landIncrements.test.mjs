import { test } from 'node:test';
import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { mkdtempSync, writeFileSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { encodeState } from '../js/stateFile.module.js';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createGeography, readTopography } from '../js/geography.module.js';
import { parseArguments } from '../scripts/landIncrements.mjs';

const SCRIPT = fileURLToPath(new URL('../scripts/landIncrements.mjs', import.meta.url));
const N = 16, mesh = buildMesh(new Grid(N)), C = mesh.nCells, geo = createGeography(mesh, readTopography('etopo1'));
const open = (i) => geo.land[i] && !geo.iceSheet[i];

function save(dir, day, since, open_, scale) {
  const field = (value) => Float32Array.from({ length: C }, (_, i) => (open(i) ? value * scale : 0));
  const land = { soil: new Float32Array(C), goalAdded: field(3), goalRemoved: field(1), goalCover: field(0.01), goalTrees: field(-0.002) };
  writeFileSync(join(dir, `syn_day${String(day).padStart(4, '0')}.bin`), encodeState({ N, K: 2, day, time: day * 86400, goalSince: since, goalOpen: open_, land }));
}
const run = (dir, day, ...extra) => spawnSync(process.execPath, [SCRIPT, join(dir, `syn_day${String(day).padStart(4, '0')}.bin`), ...extra], { encoding: 'utf8' });

test('landIncrements sums the three monthly windows of a quarter and sets the same quarter of earlier years beside it', () => {
  const dir = mkdtempSync(join(tmpdir(), 'increments-'));
  try {
    for (const [day, since] of [[30, 0], [61, 30], [91, 61]]) save(dir, day, since, false, 1);
    for (const [day, since] of [[395, 365], [426, 395], [456, 426]]) save(dir, day, since, false, 0.5);
    const out = run(dir, 456);
    assert.equal(out.status, 0, out.stderr);
    assert.match(out.stdout, /days 365–456 \(3 saves\)/);
    const row = (text, title) => text.split(title)[1].split('\n').find((l) => l.startsWith('| Open land')).split('|').map((x) => x.trim().replace('−', '-'));
    const water = row(out.stdout, '### Water'), cover = row(out.stdout, '### Vegetation'), trees = row(out.stdout, '### Tree');
    const cells = [...Array(C).keys()].filter(open).length;
    assert.equal(Number(water[2]), cells);
    assert.ok(Math.abs(Number(water[3]) - 3) < 0.06 && Math.abs(Number(water[4]) - 3) < 0.06, `3 windows of (3 - 1) x 0.5 kg/m2: ${water}`);
    const area = [...Array(C).keys()].filter(open).reduce((sum, i) => sum + mesh.areaCell[i], 0);
    assert.ok(Math.abs(Number(water[5]) - 3 * area / 1e12) < 0.02 * 3 * area / 1e12, `Gt ${water[5]}`);
    assert.ok(water[6].startsWith('+6.0 / 6.0'), `the year before sums twice as much: ${water[6]}`);
    assert.ok(Math.abs(Number(cover[3]) - 0.015) < 1e-4 && Math.abs(Number(trees[3]) + 0.003) < 1e-4, `${cover[3]} ${trees[3]}`);
    const first = run(dir, 91);
    assert.equal(first.status, 0);
    assert.match(first.stdout, /days 0–91/);
    assert.ok(!/Net \/ abs/.test(first.stdout), 'the first year has no earlier year');
    const one = run(dir, 456, '--years', '0');
    assert.ok(!/Net \/ abs/.test(one.stdout));
  } finally { rmSync(dir, { recursive: true, force: true }); }
});

test('landIncrements exits 3 and prints nothing for a state without increments or with an open window', () => {
  const dir = mkdtempSync(join(tmpdir(), 'increments-'));
  try {
    save(dir, 30, 0, true, 1);
    const halted = run(dir, 30);
    assert.deepEqual([halted.status, halted.stdout], [3, '']);
    writeFileSync(join(dir, 'syn_day0060.bin'), encodeState({ N, K: 2, day: 60, time: 60 * 86400, land: { soil: new Float32Array(C) } }));
    const plain = run(dir, 60);
    assert.deepEqual([plain.status, plain.stdout], [3, '']);
  } finally { rmSync(dir, { recursive: true, force: true }); }
});

test('landIncrements takes --days and --years and refuses other arguments', () => {
  assert.deepEqual(parseArguments(['a_day0001.bin', '--days', '30', '--years', '2']), { file: 'a_day0001.bin', days: 30, years: 2 });
  assert.throws(() => parseArguments([]), /usage/);
  assert.throws(() => parseArguments(['a', 'b']), /usage/);
});
