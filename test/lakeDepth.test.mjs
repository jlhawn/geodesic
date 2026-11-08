import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createGeography, readTopography, geographyOptions, FRESH_GEOGRAPHY } from '../js/geography.module.js';
import { createLake, decodeLakeDepths, encodeLakeDepths, readLakeDepths, lakeDepthsUrl, sourceLakeDepths, FLAKE } from '../js/physics/lake.module.js';

const DEG = 180 / Math.PI;
const mesh = buildMesh(new Grid(64)), geography = createGeography(mesh, readTopography('etopo2022'), geographyOptions(FRESH_GEOGRAPHY).geography);
const near = (lat, lon) => { let best = -1, d = Infinity; for (let i = 0; i < mesh.nCells; i++) if (geography.lake[i]) { const e = (mesh.latCell[i] * DEG - lat) ** 2 + ((((mesh.lonCell[i] * DEG - lon) % 360) + 540) % 360 - 180) ** 2; if (e < d) { d = e; best = i; } } return best; };

test('data/lakedepth.bin holds the HydroLAKES lakes of at least 100 km², largest first, and round-trips', () => {
  const lakes = decodeLakeDepths(readFileSync(lakeDepthsUrl));
  assert.ok(lakes.length > 1500 && lakes.length < 2000, `${lakes.length} lakes`);
  assert.ok(lakes.every((l, n) => l.area >= 100 && l.depth > 0 && [1, 2, 3].includes(l.source) && (n === 0 || l.area <= lakes[n - 1].area)));
  const victoria = lakes.find((l) => Math.abs(l.lat - 0.4) < 0.5 && Math.abs(l.lon - 33.2) < 0.5 && l.area > 6e4);
  assert.ok(victoria && Math.abs(victoria.depth - 41) < 1, 'Victoria, 41 m');
  assert.deepEqual(decodeLakeDepths(encodeLakeDepths(lakes.slice(0, 5))), lakes.slice(0, 5));
});

test("missingBed 'hydrolakes' gives the lakes without an ETOPO bed HydroLAKES' depth and leaves every other lake cell as it was", () => {
  const plain = createLake(mesh, geography), sourced = createLake(mesh, geography, { missingBed: 'hydrolakes' });
  const bedless = plain.cells.filter((i) => !(-geography.elevation[i] >= 1));
  assert.ok(bedless.length > 5, `${bedless.length} lake cells without a bed`);
  for (const i of plain.cells) if (!bedless.includes(i)) assert.equal(sourced.depth[i], plain.depth[i]);
  const changed = bedless.filter((i) => sourced.depth[i] !== plain.depth[i]);
  assert.ok(changed.length >= 0.8 * bedless.length, `${changed.length} of ${bedless.length} changed`);
  for (const i of plain.cells) assert.ok(sourced.depth[i] >= FLAKE.minimumDepth && sourced.depth[i] <= FLAKE.maximumDepth);
  assert.ok(Math.abs(sourced.depth[near(-1, 33)] - 41.1) < 0.1, 'Victoria');
  assert.equal(sourced.depth[near(-6, 29.5)], FLAKE.maximumDepth, 'Tanganyika at the IFS bound');
  assert.equal(plain.depth[near(-1, 33)], 10);
  assert.throws(() => createLake(mesh, geography, { missingBed: 'gldb' }), /missingBed/);
});

test('sourceLakeDepths reads only lake cells and finds nothing without lakes near them', () => {
  const depth = sourceLakeDepths(mesh, geography.lake, readLakeDepths());
  for (let i = 0; i < mesh.nCells; i++) if (!geography.lake[i]) assert.ok(Number.isNaN(depth[i]));
  assert.ok(sourceLakeDepths(mesh, geography.lake, [{ lon: -150, lat: 0, area: 1e5, depth: 30, source: 1 }]).every((d) => Number.isNaN(d)));
});

test('a lake group that reaches no pour point takes a big lake whose extent covers it while the groups at its pour point leave enough of it unclaimed', () => {
  const mesh = buildMesh(new Grid(16)), C = mesh.nCells, D = 180 / Math.PI;
  const unit = (i) => [Math.cos(mesh.latCell[i]) * Math.cos(mesh.lonCell[i]), Math.cos(mesh.latCell[i]) * Math.sin(mesh.lonCell[i]), Math.sin(mesh.latCell[i])];
  const km = (i, j) => { const p = unit(i), q = unit(j); return 6371 * Math.acos(Math.min(1, p[0] * q[0] + p[1] * q[1] + p[2] * q[2])); };
  const far = 0, near = Array.from({ length: C }, (_, j) => j).find((j) => km(far, j) > 1250 && km(far, j) < 1450);
  const big = { lon: mesh.lonCell[near] * D, lat: mesh.latCell[near] * D, area: 2e5, depth: 33 };
  const lone = new Uint8Array(C);
  lone[far] = 1;
  assert.ok(Number.isNaN(sourceLakeDepths(mesh, lone, [big], { extent: 0 })[far]), 'beyond its reach the lake is no candidate without the extent');
  assert.equal(sourceLakeDepths(mesh, lone, [big])[far], 33, 'its extent covers the far group');
  const both = Uint8Array.from(lone);
  both[near] = 1;
  const claimed = sourceLakeDepths(mesh, both, [big]);
  assert.equal(claimed[near], 33);
  assert.ok(Number.isNaN(claimed[far]), 'the group at its pour point claims all of it');
  assert.equal(sourceLakeDepths(mesh, both, [{ ...big, area: 5e5 }])[far], 33, 'a lake larger than the group at its pour point is left over');
});
