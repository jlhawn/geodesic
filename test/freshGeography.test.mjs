import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid, DEFAULT_RELAX } from '../js/grid.module.js';
import { gridMesh } from '../js/stateGrid.module.js';
import { createGeography, readTopography, geographyOptions, subgridUrl, decodeSubgrid, landSubgrid, FRESH_GEOGRAPHY } from '../js/geography.module.js';

const etopo1 = readTopography('etopo1'), etopo2022 = readTopography('etopo2022'), observed = geographyOptions({ ...FRESH_GEOGRAPHY, mask: 'raster' });

test('ETOPO 2022 moves the coast of 632 cells at N=128, its sub-grid files fit its land mask on every mesh, and ETOPO1\'s still fit ETOPO1\'s', () => {
  for (const N of [16, 32, 64, 128, 192]) {
    const mesh = gridMesh(new Grid(N)), old = createGeography(mesh, etopo1), fresh = createGeography(mesh, etopo2022, observed.geography);
    const fitted = decodeSubgrid(readFileSync(subgridUrl(N, 'etopo2022'))), kept = decodeSubgrid(readFileSync(subgridUrl(N)));
    assert.equal(fitted.relax, DEFAULT_RELAX);
    assert.ok(landSubgrid(fitted, fresh.land), `N=${N}: data/subgrid_etopo2022_N${N}.bin does not fit the ETOPO 2022 land`);
    assert.ok(landSubgrid(kept, old.land), `N=${N}: data/subgrid_N${N}.bin no longer fits the ETOPO1 land`);
    let toSea = 0, toLand = 0;
    for (let i = 0; i < mesh.nCells; i++) {
      if (old.land[i] && !fresh.land[i]) toSea++;
      if (!old.land[i] && fresh.land[i]) toLand++;
    }
    console.log(`N=${N}: ETOPO 2022 turns ${toSea} land cells to sea and ${toLand} sea cells to land`);
    if (N === 128) assert.deepEqual([toSea, toLand], [370, 262]);
  }
});
