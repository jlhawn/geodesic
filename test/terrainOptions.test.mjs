import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { createModel, dragConstants, orographyFields } from '../js/model.module.js';
import { topographyFromInt16, surfaceGeopotential, terrainOptions, terrainKey, SUBGRID_FIELDS } from '../js/geography.module.js';
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

for (const N of [16, 24]) test(`with conserveVolume each land mass at N=${N} keeps the volume the smoothing takes from it, the sea stays at 0 and the models take it`, () => {
  const model = createModel(new Grid(N), { topography, physics: false, terrain: { conserveVolume: true } });
  const { mesh, geography } = model, { land } = geography, { areaCell, cellsOnCell, nEdgesOnCell, maxEdges, nCells: C } = mesh;
  const raw = surfaceGeopotential(mesh, geography, { passes: 0 }), smoothed = surfaceGeopotential(mesh, geography), kept = surfaceGeopotential(mesh, geography, { conserveVolume: true });
  assert.deepEqual(Array.from(createModel(new Grid(N), { topography, physics: false }).surfaceGeopotential), Array.from(smoothed), 'terrain true is the default smoothing');
  assert.deepEqual(Array.from(model.surfaceGeopotential), Array.from(kept), 'the model takes the option');
  const mass = new Int32Array(C).fill(-1);
  let masses = 0, lost = 0, total = 0;
  for (let seed = 0; seed < C; seed++) {
    if (!land[seed] || mass[seed] >= 0) continue;
    const queue = [seed]; mass[seed] = masses;
    for (let n = 0; n < queue.length; n++) for (let m = 0; m < nEdgesOnCell[queue[n]]; m++) { const j = cellsOnCell[maxEdges * queue[n] + m]; if (land[j] && mass[j] < 0) { mass[j] = masses; queue.push(j); } }
    const volume = (phi) => queue.reduce((s, i) => s + areaCell[i] * phi[i], 0);
    if (queue.some((i) => smoothed[i] > 0)) assert.ok(Math.abs(volume(kept) - volume(raw)) <= 1e-9 * Math.max(1, Math.abs(volume(raw))), `land mass ${masses} of ${queue.length} cells: ${volume(kept)} against ${volume(raw)}`);
    lost += volume(raw) - volume(smoothed); total += volume(raw);
    masses++;
  }
  for (let i = 0; i < C; i++) if (!land[i]) assert.equal(kept[i], smoothed[i], `cell ${i} off the land`);
  console.log(`N=${N}: ${masses} land masses; the default smoothing loses ${(100 * lost / total).toFixed(2)} % of the land's volume`);
  assert.ok(lost > 0, 'the default smoothing loses volume');
});

test('terrain options that make the same surface share a key', () => {
  assert.equal(terrainKey(true), terrainKey({}));
  assert.equal(terrainKey(true), terrainKey({ conserveVolume: false, passes: 2 }));
  assert.equal(terrainKey({ passes: 3, conserveVolume: true }), terrainKey({ conserveVolume: true, passes: 3 }));
  assert.notEqual(terrainKey({ conserveVolume: true }), terrainKey(true));
  assert.notEqual(terrainKey({ passes: 1 }), terrainKey(true));
});

test('a subgrid file is refused over terrain of other options than its own, and the raster\'s fields are computed against the run\'s surface', () => {
  const model = createModel(new Grid(16), { topography, physics: false }), { mesh, geography } = model, g = 9.80616;
  const file = Object.fromEntries(SUBGRID_FIELDS.map(([name]) => [name, new Float64Array(mesh.nCells).fill(0.5)]));
  const kept = surfaceGeopotential(mesh, geography, { conserveVolume: true });
  assert.equal(orographyFields(mesh, topography, geography, model.surfaceGeopotential, file, g).raster, undefined);
  assert.equal(orographyFields(mesh, topography, geography, model.surfaceGeopotential, file, g, { conserveVolume: false, passes: 2 }).raster, undefined);
  const refused = orographyFields(mesh, topography, geography, kept, file, g, { conserveVolume: true });
  assert.equal(refused.raster, true);
  const computed = orographyFields(mesh, topography, geography, kept, false, g, { conserveVolume: true });
  assert.deepEqual(Array.from(refused.deviation), Array.from(computed.deviation));
  assert.equal(orographyFields(mesh, topography, geography, kept, { ...file, terrain: { conserveVolume: true } }, g, { conserveVolume: true }).raster, undefined);
});
