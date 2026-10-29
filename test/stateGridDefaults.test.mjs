import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createOcean } from '../js/ocean/layered.module.js';
import { createGeography, stateGeography } from '../js/geography.module.js';
import { lakeWater, lakeModelOf } from '../js/physics/lake.module.js';
import { decodeState } from '../js/stateFile.module.js';
import { stateGrid, savedDryColumns, maskMismatches } from '../js/stateGrid.module.js';
import { DEFAULT_RUNS } from '../js/defaultRun.module.js';

const runs = new URL('../runs/', import.meta.url);
const readParts = async (manifest) => decodeState(Buffer.concat(JSON.parse(readFileSync(new URL(manifest, runs))).parts.map(({ file }) => readFileSync(new URL(file, runs)))));

for (const N of [64, 128]) {
  const manifest = DEFAULT_RUNS[N].replace(/^runs\//, '');
  test(`the N=${N} default, ${manifest}, loads on the grid it was made on with no cell off its coast, and on the raster mask without a grid header its coast finds that grid`, async () => {
    const saved = await readParts(manifest);
    const { topography, geography } = stateGeography(saved);
    const found = stateGrid(saved, topography, { geography });
    const made = createGeography(found.mesh, topography, geography), water = made.ocean ?? made;
    assert.equal(maskMismatches(savedDryColumns(saved), Uint8Array.from(water.land, (land, i) => land || made.lake?.[i] || 0)), 0);
    if (!geography.mask) {
      const bare = stateGrid({ ...saved, grid: undefined }, topography, { geography });
      assert.deepEqual({ relax: bare.relax, source: bare.source, mismatches: bare.mismatches }, { relax: found.relax, source: 'mask', mismatches: 0 });
    }
    if (saved.N > 64) return;
    const ocean = createOcean(found.mesh, { geography: lakeModelOf(made) !== 'off' ? lakeWater(found.mesh, made) : made });
    const surfaceT = Float64Array.from(saved.surfaceT), ice = Float64Array.from(saved.ice);
    ocean.load(saved.ocean, surfaceT, ice, Float64Array.from(saved.concentration));
    const C = found.mesh.nCells, L = ocean.layers, Ls = saved.ocean.h.length / C;
    for (let i = 0; i < C; i++) {
      if (!ocean.cellOcean[i]) continue;
      let now = 0, was = 0;
      for (let k = 0; k < L; k++) now += ocean.h[k * C + i];
      for (let k = 0; k < Ls; k++) was += saved.ocean.h[k * C + i];
      assert.ok(Math.abs(now - was) < 1, `cell ${i} keeps its water through the conversion (${now} against ${was})`);
      assert.ok(Math.abs(ocean.T0[i] - saved.ocean.T[i]) < 0.02, `cell ${i} keeps its saved SST (${ocean.T0[i]} vs ${saved.ocean.T[i]})`);
      assert.equal(surfaceT[i], saved.surfaceT[i]);
    }
  });
}
