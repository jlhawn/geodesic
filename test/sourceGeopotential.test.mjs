import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { topographyFromInt16, surfaceGeopotential, decodeSubgrid } from '../js/geography.module.js';
import { sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { withCadence } from '../js/cadence.module.js';
import { savedTerrain } from '../js/stateFile.module.js';

/*
 * The page's worker rebalances a saved run's surface pressure from the
 * terrain a physics-free model on the saved mesh with the current
 * topography has (js/model.worker.js, sourceFor) to the model's own. At
 * the same resolution it takes that terrain from the model it built, so
 * the two must be the same geopotential, cell for cell, whether the model
 * runs with terrain or without.
 */
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const levels = sigmaInterfaces('bl36');

for (const N of [16, 32]) {
  test(`at N=${N} the page's model has the terrain of the saved mesh's physics-free model, with terrain on and off`, () => {
    const source = createModel(new Grid(N), { physics: false, levels, topography });
    assert.ok(source.surfaceGeopotential && source.surfaceGeopotential.some((phi) => phi > 0));
    const subgrid = decodeSubgrid(readFileSync(new URL(`../data/subgrid_N${N}.bin`, import.meta.url)));
    for (const terrain of [true, false]) {
      const options = { topography, subgrid, terrain, levels };
      const model = createModel(new Grid(N), { ...options, ...withCadence(options, 1350 * 16 / N, N) });
      const own = model.surfaceGeopotential ?? surfaceGeopotential(model.mesh, model.geography);
      assert.equal(model.surfaceGeopotential === null, !terrain);
      assert.equal(own.length, source.surfaceGeopotential.length);
      for (let i = 0; i < own.length; i++) assert.equal(own[i], source.surfaceGeopotential[i], `terrain ${terrain}, cell ${i}`);
    }
  });
}

test('a state saved over terrain of other options starts on it: the page\'s model and the saved mesh\'s physics-free model take the terrain its header records, so at the same resolution the rebalance has nothing to move', () => {
  const terrain = { conserveVolume: true }, header = { N: 16, terrain: true, terrainOptions: terrain };
  assert.deepEqual(savedTerrain(header), terrain);
  assert.deepEqual(savedTerrain({ terrain: true, given: { terrain } }), terrain);
  assert.equal(savedTerrain({ terrain: true }), undefined);
  const recorded = savedTerrain(header) ?? true, subgrid = decodeSubgrid(readFileSync(new URL('../data/subgrid_N16.bin', import.meta.url)));
  const source = createModel(new Grid(16), { physics: false, levels, topography, terrain: recorded }), plain = createModel(new Grid(16), { physics: false, levels, topography });
  const options = { topography, subgrid, terrain: recorded, levels }, model = createModel(new Grid(16), { ...options, ...withCadence(options, 1350, 16) });
  assert.deepEqual(Array.from(model.surfaceGeopotential), Array.from(source.surfaceGeopotential));
  assert.ok(model.surfaceGeopotential.some((phi, i) => phi !== plain.surfaceGeopotential[i]), 'the recorded terrain is not the default one');
});
