import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { geographyOptions, readTopography, FRESH_GEOGRAPHY } from '../js/geography.module.js';

const N = 32, DT = 1350 * 16 / N, topography = readTopography('etopo2022'), sheets = geographyOptions({ ...FRESH_GEOGRAPHY });

function prepare(passed) {
  const model = createModel(new Grid(N), { topography, geography: sheets.geography, land: sheets.land, lake: { model: 'flake' } });
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i] || model.geography.lake?.[i]) model.state[6][i] = 0;
  model.ocean.initialize(model.state[3], model.state[6]);
  model.land.initialize();
  const saved = model.lake.serialize();
  for (const i of model.lake.cells) saved.passed[i] = passed;
  model.lake.load(saved);
  return model;
}

test(`at N=${N} a CPU model whose lakes load a running total of passed water hands the sea only what they pass after the load`, () => {
  const fresh = (m) => { let t = 0; for (let i = 0; i < m.mesh.nCells; i++) t += m.ocean.fresh[i] * m.mesh.areaCell[i]; return t; };
  const base = prepare(0), loaded = prepare(1000);
  let lakeArea = 0;
  for (const i of base.lake.cells) lakeArea += base.mesh.areaCell[i];
  for (let n = 0; n < base.ocean.everySteps - 1; n++) { base.step(DT); loaded.step(DT); }
  const gap = Math.abs(fresh(loaded) - fresh(base)), load = 1000 * lakeArea;
  console.log(`${base.lake.cells.length} lake cells: the sea's fresh water differs by ${(gap / load).toExponential(1)} of the loaded total`);
  assert.ok(fresh(base) !== 0, 'the lakes pass water');
  assert.ok(gap < 1e-9 * load, `the sea's fresh water differs by ${gap} kg against a load of ${load} kg`);
});

test(`at N=${N} a saved lake carries its column depths, and a run whose lakes have other depths refuses it`, () => {
  const model = prepare(0), saved = model.lake.serialize(), i = model.lake.cells[0];
  assert.equal(saved.depth[i], model.lake.depth[i]);
  model.lake.load(saved);
  saved.depth[i] += 5;
  assert.throws(() => model.lake.load(saved), /other depths/);
});
