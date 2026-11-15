import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createModel } from '../js/model.module.js';
import { createOcean } from '../js/ocean/layered.module.js';
import { iceExchangeOptions, DEFAULT_ICE_EXCHANGE, DEFAULT_ICE_ALBEDO_PHASES } from '../js/physics/ice.module.js';
import { PUBLISHED_OCEAN, PUBLISHED_ICE } from '../explain/published.module.js';

test('the explainer builds its oceans and sea ice as its shipped data and figures were made, the exchange instant and no iceAlbedoPhases where the defaults differ, in every builder that can form sea ice', () => {
  const instant = iceExchangeOptions('instant');
  assert.notDeepEqual(iceExchangeOptions(DEFAULT_ICE_EXCHANGE), instant);
  assert.notEqual(DEFAULT_ICE_ALBEDO_PHASES, false);
  assert.deepEqual(createOcean(buildMesh(new Grid(16)), PUBLISHED_OCEAN).iceExchange, instant);
  const model = createModel(new Grid(16), { ocean: false, ice: PUBLISHED_ICE });
  assert.deepEqual([model.seaIce.iceExchange, model.seaIce.iceAlbedoPhases], [instant, false]);
  const source = (path) => readFileSync(new URL(`../${path}`, import.meta.url), 'utf8');
  for (const [path, pattern] of [
    ['explain/oceanColumn.module.js', /createOcean\(m, \{ \.\.\.PUBLISHED_OCEAN, /],
    ['explain/oceanCases.module.js', /const START = \{ \.\.\.PUBLISHED_OCEAN, /],
    ['explain/seaClasses.module.js', /createOcean\(mesh, \{ iceExchange: PUBLISHED_OCEAN\.iceExchange \}\)/],
    ['explain/plumeColumn.module.js', /createModel\([^;]*\{ ocean: false, ice: PUBLISHED_ICE, /],
    ['explain/seaIceColumn.module.js', /iceExchange: 'instant', iceAlbedoPhases: false/],
    ['scripts/explainLayers.mjs', /ice: \{ \.\.\.PUBLISHED_ICE, \.\.\.ICE \}/],
    ['scripts/explainOcean.mjs', /const madeUnder = resumed \? savedIceExchange\(\{ options: base \}\) : PUBLISHED_ICE\.iceExchange;/],
  ]) assert.match(source(path), pattern, path);
});
