import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid, DEFAULT_RELAX } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createOcean } from '../js/ocean/layered.module.js';
import { createGeography, syntheticTopography } from '../js/geography.module.js';
import { encodeState, decodeState } from '../js/stateFile.module.js';
import { stateGrid, stateRelax, gridHeader, gridMesh, savedDryColumns, maskMismatches, STATE_RELAX } from '../js/stateGrid.module.js';

const N = 16, DEG = Math.PI / 180;
const coast = syntheticTopography(360, 720, (lat, lon) => (lat > 5 * DEG * Math.sin(7 * lon) + 3 * DEG * Math.cos(11 * lon) && Math.abs(lon) < 100 * DEG ? 800 : -4000));

function savedOn(relax) {
  const mesh = gridMesh(new Grid(N, { relax }));
  const geography = createGeography(mesh, coast);
  const ocean = createOcean(mesh, { geography });
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => 275 + 25 * Math.cos(lat) ** 2);
  ocean.initialize(surfaceT, new Float64Array(mesh.nCells));
  return { N, K: 1, day: 0, time: 0, surfaceT, ocean: ocean.serialize(), land: geography.land };
}

test('a grid keeps its Lloyd passes, its mesh carries them, and a state written with them reads them back from its header', async () => {
  assert.equal(new Grid(N).relaxPasses, DEFAULT_RELAX);
  const grid = new Grid(N, { relax: 8 }), mesh = gridMesh(grid);
  assert.equal(grid.relaxPasses, 8);
  assert.equal(mesh.relax, 8);
  assert.deepEqual(gridHeader(mesh), { grid: { relax: 8 } });
  assert.deepEqual(gridHeader(grid), { grid: { relax: 8 } });
  assert.deepEqual(gridHeader({ nCells: 1 }), {});
  const model = createModel(mesh, { physics: false });
  const bytes = encodeState({ N, ...gridHeader(model.mesh), K: model.core.K, day: 3, time: 3 * 86400, pi: model.state[0] });
  const saved = await decodeState(bytes);
  assert.deepEqual(saved.grid, { relax: 8 });
  assert.equal(saved.pi.length, model.state[0].length);
  const found = stateGrid(saved, null);
  assert.equal(found.relax, 8);
  assert.equal(found.source, 'header');
  assert.equal(found.grid.relaxPasses, 8);
  assert.deepEqual(Array.from(found.mesh.xCell), Array.from(mesh.xCell));
  const json = await decodeState(new TextEncoder().encode(JSON.stringify({ N, grid: { relax: 0 } })));
  assert.equal(stateGrid(json, null).relax, 0);
});

test('a state saved without the header is matched to the grid whose land mask is its ocean\'s dry columns, on a synthetic N = 16 coast', () => {
  const [first, second] = STATE_RELAX.map((relax) => savedOn(relax));
  const disagree = maskMismatches(first.land, second.land);
  assert.ok(disagree > 0, `the two grids put the coast in different cells (${disagree})`);
  for (const [saved, relax] of [[first, STATE_RELAX[0]], [second, STATE_RELAX[1]]]) {
    assert.deepEqual(Array.from(savedDryColumns(saved)), Array.from(saved.land, (l) => (l ? 1 : 0)));
    const found = stateGrid(saved, coast);
    assert.equal(found.relax, relax);
    assert.equal(found.source, 'mask');
    assert.equal(found.mismatches, 0);
    assert.equal(found.grid.relaxPasses, relax);
    assert.equal(found.mesh.relax, relax);
    assert.deepEqual(stateRelax(saved, coast), { relax, source: 'mask', mismatches: 0 });
    const loaded = createGeography(found.mesh, coast).land;
    assert.equal(maskMismatches(savedDryColumns(saved), loaded), 0);
  }
  const water = { ...second, ocean: { ...second.ocean, h: Float64Array.from(second.ocean.h, (h, n) => (h > 0 ? h : 1)) } };
  assert.deepEqual(stateRelax(water, coast), { relax: DEFAULT_RELAX, source: 'default', mismatches: null }, 'an ocean with no dry column is on the default grid');
  assert.deepEqual(stateRelax(second, null), { relax: DEFAULT_RELAX, source: 'default', mismatches: null }, 'nothing to match against');
  assert.deepEqual(stateRelax({ N, K: 1 }, coast), { relax: DEFAULT_RELAX, source: 'default', mismatches: null }, 'no ocean');
  assert.equal(stateRelax({ ...second, grid: { relax: STATE_RELAX[0] } }, coast).relax, STATE_RELAX[0], 'the header wins over the coast');
});
