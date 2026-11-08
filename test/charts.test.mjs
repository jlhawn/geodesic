import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { stateGeography, FRESH_GEOGRAPHY } from '../js/geography.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { gridHeader } from '../js/stateGrid.module.js';
import { levelFields } from '../js/levels.module.js';
import { chartModel, chartFor } from '../js/charts.module.js';

test('the charts diagnose a state on the terrain it records: over Tibet the 500 hPa surface stands above the ground, where the run itself puts it', () => {
  const mesh = new Grid(16), header = { N: 16, ...gridHeader(mesh), ...FRESH_GEOGRAPHY, terrain: true, day: 0 };
  const { topography, geography } = stateGeography(header), DEG = 180 / Math.PI;
  for (const terrain of [true, { conserveVolume: true }]) {
    const run = createModel(mesh, { topography, geography, physics: false, terrain }), phis = run.surfaceGeopotential;
    const [pi, theta, u] = initializeState(run, {});
    const state = { ...header, ...(terrain === true ? {} : { terrainOptions: terrain }), levels: run.core.levels, pi, theta, u };
    const { model } = chartModel(state, topography, geography);
    assert.deepEqual(Array.from(model.surfaceGeopotential), Array.from(phis));
    const heights = (core) => {
      const p = Float64Array.from(pi), t = Float64Array.from(theta);
      core.diagnose(p, t);
      return Array.from(levelFields(core, p, t, () => new Float64Array(3 * run.mesh.nCells), 500).height);
    };
    const height = heights(model.core);
    assert.deepEqual(height, heights(run.core));
    let tibet = -1;
    for (let i = 0; i < run.mesh.nCells; i++) {
      const lat = run.mesh.latCell[i] * DEG, lon = run.mesh.lonCell[i] * DEG;
      if (lat > 27 && lat < 40 && lon > 75 && lon < 105 && (tibet < 0 || phis[i] > phis[tibet])) tibet = i;
    }
    const ground = phis[tibet] / run.core.diagnostics.g;
    assert.ok(ground > 2000, `Tibet's highest cell stands ${ground.toFixed(0)} m up`);
    assert.ok(height[tibet] > ground + 1000, `500 hPa at ${height[tibet].toFixed(0)} m over ground at ${ground.toFixed(0)} m`);
    assert.match(chartFor(state, '500', topography, geography), /500 hPa geopotential height/);
  }
});
