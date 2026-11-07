import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { createModel, dragConstants } from '../js/model.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { topographyFromInt16, surfaceGeopotential, terrainOptions, terrainKey, createGeography, readTopography, geographyOptions, FRESH_GEOGRAPHY } from '../js/geography.module.js';
import { orographicColumn, namedConstants, OROGRAPHY_DEFAULTS, OROGRAPHY_CONSTANTS, LOTT_MILLER } from '../js/physics/orography.module.js';
import { orographyConstants } from '../js/gpu/orography.gpu.js';

const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);

function ridgeColumn() {
  const K = 30, z = Float64Array.from({ length: K }, (_, k) => 150 * (K - k) - 75), theta = Float64Array.from(z, (h) => 290 * Math.exp(1e-4 * h / 9.80616));
  const p = Float64Array.from(z, (h) => 1e5 * Math.exp(-h / 8000)), edge = (h) => 1e5 * Math.exp(-h / 8000);
  return { z, p, rho: Float64Array.from(p, () => 1.2), theta, east: new Float64Array(K).fill(10), north: new Float64Array(K), pTop: Float64Array.from(z, (h) => edge(h + 75)), pBottom: Float64Array.from(z, (h) => edge(Math.max(0, h - 75))) };
}

test('the orography constants \'cy41r2\' take the IFS\'s operational G, 0.3, over the defaults; constants given by name win, and the raster\'s fields take the set over LOTT_MILLER', () => {
  assert.deepEqual(namedConstants({}), {});
  assert.deepEqual(namedConstants({ constants: 'cy41r2' }), { ...OROGRAPHY_DEFAULTS, waveDrag: 0.3 });
  assert.equal(namedConstants({ constants: 'cy41r2', waveDrag: 0.5 }).waveDrag, 0.5);
  assert.throws(() => namedConstants({ constants: 'cy47r3' }), /cy41r2/);
  assert.deepEqual(dragConstants({ raster: true }, { constants: 'cy41r2' }), OROGRAPHY_CONSTANTS.cy41r2);
  assert.deepEqual(dragConstants({ raster: true }, {}), LOTT_MILLER);
  assert.deepEqual(dragConstants({}, { waveDrag: 2 }), { waveDrag: 2 });
  assert.equal(createModel(new Grid(16), { topography, orography: { constants: 'cy41r2' } }).orography.options.waveDrag, 0.3);
  assert.match(orographyConstants({ orography: dragConstants({}, { constants: 'cy41r2' }) }), /const ORO_G: f32 = 0\.3;/);
});

test('with \'cy41r2\' an unblocked ridge launches 0.3/1.23 of the default\'s wave stress and the same blocking', () => {
  const sub = { deviation: 200, anisotropy: 1, orientation: Math.PI / 2, slope: 0.02 };
  const base = orographicColumn(sub, ridgeColumn(), OROGRAPHY_DEFAULTS), ifs = orographicColumn(sub, ridgeColumn(), OROGRAPHY_CONSTANTS.cy41r2);
  assert.ok(base.launch > 0, `the default launches ${base.launch} N/m²`);
  assert.ok(Math.abs(ifs.launch / base.launch - 0.3 / 1.23) < 1e-12, `launch ratio ${ifs.launch / base.launch}`);
  assert.equal(ifs.blocking, base.blocking);
  assert.deepEqual(Array.from(ifs.beta), Array.from(base.beta));
});

test('the terrain options take conserveVolume and passes; true keeps the defaults', () => {
  assert.deepEqual(terrainOptions(true), {});
  assert.deepEqual(terrainOptions({ conserveVolume: true }), { conserveVolume: true });
  assert.deepEqual(terrainOptions({ passes: 0 }), { passes: 0 });
  assert.throws(() => terrainOptions({ smoothing: 'volume' }), /smoothing/);
  assert.throws(() => terrainOptions({ conserveVolume: 1 }), /conserveVolume/);
  assert.throws(() => terrainOptions({ passes: 1.5 }), /passes/);
});

const landMasses = ({ nCells: C, cellsOnCell, nEdgesOnCell, maxEdges }, land) => {
  const mass = new Int32Array(C).fill(-1), masses = [];
  for (let seed = 0; seed < C; seed++) {
    if (!land[seed] || mass[seed] >= 0) continue;
    const queue = [seed]; mass[seed] = masses.length;
    for (let n = 0; n < queue.length; n++) for (let m = 0; m < nEdgesOnCell[queue[n]]; m++) { const j = cellsOnCell[maxEdges * queue[n] + m]; if (land[j] && mass[j] < 0) { mass[j] = masses.length; queue.push(j); } }
    masses.push(queue);
  }
  return masses;
};

function keptVolume(mesh, geography, label) {
  const { land, lake } = geography, { areaCell, nCells: C } = mesh;
  const raw = surfaceGeopotential(mesh, geography, { passes: 0 }), kept = surfaceGeopotential(mesh, geography, { conserveVolume: true }), masses = landMasses(mesh, land);
  let worst = 0, islands = 0;
  for (const cells of masses) {
    const volume = (phi) => cells.reduce((s, i) => s + areaCell[i] * phi[i], 0), scale = cells.reduce((s, i) => s + areaCell[i] * Math.abs(raw[i]), 0);
    const error = Math.abs(volume(kept) - volume(raw));
    assert.ok(error <= 1e-12 * scale, `${label}: a land mass of ${cells.length} cells holds ${volume(kept)} against ${volume(raw)}`);
    if (scale > 0) worst = Math.max(worst, error / scale);
    if (cells.length === 1) { islands += raw[cells[0]] > 0; assert.equal(kept[cells[0]], raw[cells[0]], `${label}: the one-cell island ${cells[0]} keeps its height`); }
  }
  for (let i = 0; i < C; i++) if (!land[i]) assert.equal(kept[i], lake?.[i] ? raw[i] : 0, `${label}: water cell ${i} stays as it was`);
  return { kept, masses: masses.length, islands, worst };
}

for (const N of [16, 24]) test(`with conserveVolume each land mass at N=${N} keeps its volume to 1e-12 through passes that trade it between neighbors, the sea stays at 0 and the models take it`, () => {
  const model = createModel(new Grid(N), { topography, physics: false, terrain: { conserveVolume: true } });
  const { mesh, geography } = model, { kept, masses, worst } = keptVolume(mesh, geography, `N=${N}`);
  assert.deepEqual(Array.from(createModel(new Grid(N), { topography, physics: false }).surfaceGeopotential), Array.from(surfaceGeopotential(mesh, geography)), 'terrain true is the default smoothing');
  assert.deepEqual(Array.from(model.surfaceGeopotential), Array.from(kept), 'the model takes the option');
  console.log(`N=${N}: ${masses} land masses, each volume kept to ${worst.toExponential(1)}`);
});

for (const N of [64, 128]) test(`with conserveVolume each land mass of the page grid's geography at N=${N} keeps its volume to 1e-12, the sea stays at 0, lakes at their levels, and one-cell islands keep their heights`, () => {
  const mesh = buildMesh(new Grid(N), { omega: 7.2921e-5 });
  const geography = createGeography(mesh, readTopography(FRESH_GEOGRAPHY.topography), geographyOptions(FRESH_GEOGRAPHY).geography);
  assert.equal(geography.mask, 'n1024');
  const { masses, islands, worst } = keptVolume(mesh, geography, `N=${N}`);
  assert.ok(islands > 0, `N=${N} has a raised one-cell island`);
  console.log(`N=${N}: ${masses} land masses, ${islands} raised one-cell islands, each volume kept to ${worst.toExponential(1)}`);
});

test('with conserveVolume no land cell depends on a water cell: a one-cell island keeps its height, a flat plateau stays flat to its coast, and raising the sea beside real coasts moves no land cell', () => {
  const mesh = buildMesh(new Grid(16), { omega: 7.2921e-5 }), { nCells: C, latCell, lonCell } = mesh;
  const distance = (i, lat, lon) => Math.acos(Math.min(1, Math.sin(latCell[i]) * Math.sin(lat) + Math.cos(latCell[i]) * Math.cos(lat) * Math.cos(lonCell[i] - lon)));
  const land = Uint8Array.from({ length: C }, (_, i) => (distance(i, 0.5, 0) < 0.5 ? 1 : 0)), elevation = Float64Array.from(land, (l) => 2000 * l);
  let island = 0;
  for (let i = 0; i < C; i++) if (distance(i, -0.5, Math.PI) < distance(island, -0.5, Math.PI)) island = i;
  land[island] = 1; elevation[island] = 1000;
  const flat = { land, elevation }, raw = surfaceGeopotential(mesh, flat, { passes: 0 }), kept = surfaceGeopotential(mesh, flat, { conserveVolume: true });
  assert.ok(surfaceGeopotential(mesh, flat)[island] < raw[island], 'the default smoothing lowers the island');
  assert.deepEqual(Array.from(kept), Array.from(raw), 'the island and the plateau keep their heights');
  const { geography } = createModel(new Grid(16), { topography, physics: false }), g = 9.80616;
  const shore = Float64Array.from(surfaceGeopotential(mesh, geography, { passes: 0 }), (phi) => phi / g);
  let seed = 7;
  const random = () => (seed = (seed * 16807) % 2147483647) / 2147483647;
  const raised = Float64Array.from(shore, (h, i) => (geography.land[i] ? h : 500 + 1000 * random()));
  const before = surfaceGeopotential(mesh, { ...geography, surface: shore }, { conserveVolume: true }), after = surfaceGeopotential(mesh, { ...geography, surface: raised }, { conserveVolume: true });
  const plain = surfaceGeopotential(mesh, { ...geography, surface: raised });
  let coasts = 0;
  for (let i = 0; i < C; i++) {
    if (!geography.land[i]) { assert.equal(after[i], g * raised[i], `water cell ${i} stays as it was`); continue; }
    assert.equal(after[i], before[i], `land cell ${i} does not see the water`);
    coasts += plain[i] !== surfaceGeopotential(mesh, { ...geography, surface: shore })[i];
  }
  assert.ok(coasts > 0, 'the default smoothing sees the raised water');
});

test('with conserveVolume the passes take at least as much grid-scale roughness off a random field over land without coasts as the default smoothing', () => {
  const mesh = buildMesh(new Grid(16), { omega: 7.2921e-5 }), { nCells: C, areaCell, cellsOnCell, nEdgesOnCell, maxEdges } = mesh;
  let seed = 20261107;
  const random = () => (seed = (seed * 16807) % 2147483647) / 2147483647;
  const geography = { land: new Uint8Array(C).fill(1), elevation: Float64Array.from({ length: C }, () => 1000 * random()) };
  const roughness = (phi) => {
    let sum = 0, area = 0;
    for (let i = 0; i < C; i++) {
      let around = 0;
      for (let m = 0; m < nEdgesOnCell[i]; m++) around += phi[cellsOnCell[maxEdges * i + m]];
      const d = phi[i] - around / nEdgesOnCell[i];
      sum += areaCell[i] * d * d; area += areaCell[i];
    }
    return Math.sqrt(sum / area);
  };
  const before = roughness(surfaceGeopotential(mesh, geography, { passes: 0 })), smoothed = roughness(surfaceGeopotential(mesh, geography)), kept = roughness(surfaceGeopotential(mesh, geography, { conserveVolume: true }));
  console.log(`rms of h less its neighbors' mean: ${(smoothed / before).toFixed(5)} of the field's after the default passes, ${(kept / before).toFixed(5)} with conserveVolume`);
  assert.ok(kept <= smoothed, `${kept} against ${smoothed}`);
});

test('terrain options that make the same surface share a key', () => {
  assert.equal(terrainKey(true), terrainKey({}));
  assert.equal(terrainKey(true), terrainKey({ conserveVolume: false, passes: 2 }));
  assert.equal(terrainKey({ passes: 3, conserveVolume: true }), terrainKey({ conserveVolume: true, passes: 3 }));
  assert.notEqual(terrainKey({ conserveVolume: true }), terrainKey(true));
  assert.notEqual(terrainKey({ passes: 1 }), terrainKey(true));
});

test('a run over terrain of other options keeps the bundled subgrid file, its σ_flt and the named drag constants', () => {
  const plain = createModel(new Grid(16), { topography }), kept = createModel(new Grid(16), { topography, terrain: { conserveVolume: true } });
  assert.notDeepEqual(Array.from(kept.surfaceGeopotential), Array.from(plain.surfaceGeopotential));
  for (const model of [plain, kept]) {
    assert.equal(model.orography.fields.raster, undefined);
    assert.ok(model.orography.fields.filtered.some((value) => value > 0));
    assert.deepEqual(model.orography.options, OROGRAPHY_DEFAULTS);
  }
  assert.deepEqual(Array.from(kept.orography.fields.deviation), Array.from(plain.orography.fields.deviation));
});
