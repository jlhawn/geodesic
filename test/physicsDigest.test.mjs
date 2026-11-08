import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { gridMesh } from '../js/stateGrid.module.js';
import { areaMean, bandMask, hemisphereArea, compareDigests, runDigest } from '../scripts/physicsDigest.mjs';

const mesh = gridMesh(new Grid(16)), C = mesh.nCells;
const sphere = () => { let a = 0; for (let i = 0; i < C; i++) a += mesh.areaCell[i]; return a; };

test('the digest\x27s means weigh cells by area, over a mask or a latitude band', () => {
  assert.ok(Math.abs(areaMean(mesh, new Float64Array(C).fill(7.5)) - 7.5) < 1e-12);
  const sinLat = Float64Array.from({ length: C }, (_, i) => Math.sin(mesh.latCell[i]));
  assert.ok(Math.abs(areaMean(mesh, sinLat)) < 1e-6, `global mean of sin(lat) ${areaMean(mesh, sinLat)}`);
  const north = Uint8Array.from({ length: C }, (_, i) => (mesh.latCell[i] > 0 ? 1 : 0));
  assert.ok(areaMean(mesh, sinLat, north) > 0.45 && areaMean(mesh, sinLat, north) < 0.55, `northern mean of sin(lat) ${areaMean(mesh, sinLat, north)}`);
  const tropics = bandMask(mesh, -30, 30), share = areaMean(mesh, Float64Array.from(tropics));
  assert.ok(Math.abs(share - 0.5) < 0.02, `30S-30N holds ${share} of the sphere`);
  const inner = bandMask(mesh, -30, 30, north);
  for (let i = 0; i < C; i++) assert.equal(inner[i], tropics[i] && north[i] ? 1 : 0);
  assert.ok(Number.isNaN(areaMean(mesh, sinLat, new Uint8Array(C))));
});

test('sea-ice area sums cover times cell area over the sea cells of one hemisphere, in Mkm²', () => {
  const full = new Float64Array(C).fill(1), all = new Uint8Array(C).fill(1), half = sphere() / 2 / 1e12;
  const n = hemisphereArea(mesh, full, all, true), s = hemisphereArea(mesh, full, all, false);
  assert.ok(Math.abs(n + s - 2 * half) < 1e-9 * half && Math.abs(n - half) < 0.03 * half, `${n} and ${s} Mkm² against ${half}`);
  assert.ok(Math.abs(half - 255.0) < 1, `a hemisphere is ${half} Mkm²`);
  const sea = Uint8Array.from({ length: C }, (_, i) => i % 3), cover = Float64Array.from({ length: C }, (_, i) => (i % 7) / 7);
  let expected = 0;
  for (let i = 0; i < C; i++) if (sea[i] && mesh.latCell[i] > 0) expected += mesh.areaCell[i] * cover[i] / 1e12;
  assert.ok(Math.abs(hemisphereArea(mesh, cover, sea, true) - expected) < 1e-12 * half);
});

test('compare lists the rows that moved, largest relative change first, and the rows added or dropped', () => {
  const earlier = { stateHash: 'a', rows: [{ name: 'x', units: 'K', value: 1 }, { name: 'y', units: 'W/m²', value: 200 }, { name: 'z', units: '', value: 3 }, { name: 'gone', units: '', value: 0 }] };
  assert.deepEqual(compareDigests(earlier, earlier), { moved: [], added: [], removed: [], sameState: true });
  const later = { stateHash: 'b', rows: [{ name: 'x', units: 'K', value: 1.5 }, { name: 'y', units: 'W/m²', value: 202 }, { name: 'z', units: '', value: 3 + 1e-12 }, { name: 'new', units: '', value: 1 }] };
  const { moved, added, removed, sameState } = compareDigests(earlier, later);
  assert.deepEqual(moved.map((m) => m.name), ['x', 'y']);
  assert.ok(Math.abs(moved[0].delta - 0.5) < 1e-15 && Math.abs(moved[1].relative - 2 / 202) < 1e-15);
  assert.deepEqual(added, ['new']);
  assert.deepEqual(removed, ['gone']);
  assert.equal(sameState, false);
});

test('one step of the default-physics digest gives every row a finite value, the same twice', async () => {
  const [a, b] = [await runDigest({ steps: 1 }), await runDigest({ steps: 1 })];
  console.log(`${a.rows.length} rows, state ${a.stateHash}`);
  assert.ok(a.rows.length >= 30 && a.rows.every((r) => Number.isFinite(r.value)), a.rows.filter((r) => !Number.isFinite(r.value)).map((r) => r.name).join(', '));
  assert.equal(a.stateHash, b.stateHash);
  assert.equal(compareDigests(a, b).moved.length, 0);
});

test('the 16-step digest reproduces the digest of record: the same state, no row moved, added or dropped', async () => {
  const record = JSON.parse(readFileSync(new URL('./helpers/physicsDigest.json', import.meta.url), 'utf8'));
  const now = await runDigest({ steps: record.steps }), { moved, added, removed, sameState } = compareDigests(record, now);
  const made = `the record was made on ${record.commit} (data/ ${record.data}) under Node ${record.node}; this run is ${now.commit} (data/ ${now.data}) under Node ${now.node}`;
  const groups = Object.keys(now.groupHashes).filter((name) => now.groupHashes[name] !== record.groupHashes?.[name]);
  assert.deepEqual([moved.map((m) => `${m.name} ${m.before} → ${m.after}`), added, removed], [[], [], []], made);
  assert.ok(sameState, `state ${now.stateHash} against the record's ${record.stateHash}, differing in ${groups.join(', ') || 'no group'}; ${made}`);
  assert.deepEqual(now.hashed, record.hashed);
});
