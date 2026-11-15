import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, readdirSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createModel } from '../js/model.module.js';
import { createOcean } from '../js/ocean/layered.module.js';
import { iceExchangeOptions, DEFAULT_ICE_EXCHANGE, DEFAULT_ICE_ALBEDO_PHASES } from '../js/physics/ice.module.js';
import { PUBLISHED_OCEAN, PUBLISHED_ICE } from '../explain/published.module.js';
import { seaIceColumn, iceGrowth } from '../explain/seaIceColumn.module.js';
import { ekmanCase, lensCase, bumpCase, basinCase } from '../explain/oceanCases.module.js';

const instant = iceExchangeOptions('instant');
const source = (path) => readFileSync(new URL(`../${path}`, import.meta.url), 'utf8');

test('the published pins differ from the model defaults', () => {
  assert.notDeepEqual(iceExchangeOptions(DEFAULT_ICE_EXCHANGE), instant);
  assert.notEqual(DEFAULT_ICE_ALBEDO_PHASES, false);
});

test('the explainer builds the objects it exposes as its shipped data and figures were made: the exchange instant, the albedo phases off', () => {
  assert.deepEqual(createOcean(buildMesh(new Grid(16)), PUBLISHED_OCEAN).iceExchange, instant);
  const model = createModel(new Grid(16), { ocean: false, ice: PUBLISHED_ICE });
  assert.deepEqual([model.seaIce.iceExchange, model.seaIce.iceAlbedoPhases], [instant, false]);
  const seas = [
    ['iceGrowth', iceGrowth({ days: 1 }).sea],
    ['seaIceColumn', seaIceColumn({ years: 1 }).sea],
    ['seaIceColumn dark', seaIceColumn({ years: 1, ice: 'dark' }).sea],
  ];
  for (const [name, sea] of seas) assert.deepEqual([name, sea.iceExchange, sea.iceAlbedoPhases], [name, instant, false]);
  const sims = [
    ['ekmanCase', ekmanCase({ N: 16 })], ['lensCase', lensCase()], ['bumpCase', bumpCase()], ['basinCase', basinCase()],
  ];
  for (const [name, sim] of sims) assert.deepEqual([name, sim.ocean.iceExchange], [name, instant]);
});

test('the builders that run below N = 16, or that cannot be imported without running, take the pins in their source', () => {
  for (const [path, pattern, count] of [
    ['explain/oceanColumn.module.js', /createOcean\(m, \{ \.\.\.PUBLISHED_OCEAN, /, 1],
    ['explain/oceanCases.module.js', /const START = \{ \.\.\.PUBLISHED_OCEAN, /, 1],
    ['explain/seaClasses.module.js', /createOcean\(mesh, \{ iceExchange: PUBLISHED_OCEAN\.iceExchange \}\)/, 1],
    ['explain/plumeColumn.module.js', /createModel\([^;]*\{ ocean: false, ice: PUBLISHED_ICE, /, 1],
    ['explain/seaIceColumn.module.js', /iceExchange: 'instant', iceAlbedoPhases: false/g, 2],
    ['scripts/explainLayers.mjs', /ice: \{ \.\.\.PUBLISHED_ICE, \.\.\.ICE \}/, 1],
    ['scripts/explainOcean.mjs', /const madeUnder = resumed \? savedIceExchange\(\{ options: base \}\) : PUBLISHED_ICE\.iceExchange;/, 1],
    ['scripts/explainOcean.mjs', /\{ ocean: \{ everySteps: 4 \}, ice: DARK_ICE \}/, 1],
    ['scripts/explainOcean.mjs', /const DARK_ICE = \{[^}]*iceAlbedoPhases: false \}/, 1],
    ['scripts/explainOcean.mjs', /\{ \.\.\.GIVEN_ICE, iceAlbedoPhases: resumed \? savedIceAlbedoPhases\(\{ options: base \}\) : PUBLISHED_ICE\.iceAlbedoPhases \}/, 1],
    ['scripts/explainOcean.mjs', /const OCEAN = \{ iceExchange: madeUnder, \.\.\.base\.ocean, /, 1],
  ]) assert.equal(source(path).match(new RegExp(pattern.source, pattern.flags.replace('g', '') + 'g'))?.length ?? 0, count, `${path} ${pattern}`);
});

test('every call that builds a sea ice, ocean or model in the explainer is one the pins above cover', () => {
  const listed = {
    'explain/oceanCases.module.js': 3,
    'explain/oceanColumn.module.js': 1,
    'explain/plumeColumn.module.js': 1,
    'explain/seaClasses.module.js': 1,
    'explain/seaIceColumn.module.js': 2,
    'scripts/explainLayers.mjs': 1,
    'scripts/explainOcean.mjs': 1,
    'scripts/explainSurfaceBudget.mjs': 1,
  };
  const files = [
    ...readdirSync(new URL('../explain/', import.meta.url)).filter((f) => f.endsWith('.js')).map((f) => `explain/${f}`),
    ...readdirSync(new URL('../explain/figures/', import.meta.url)).filter((f) => f.endsWith('.js')).map((f) => `explain/figures/${f}`),
    ...readdirSync(new URL('../scripts/', import.meta.url)).filter((f) => /^explain.*\.mjs$/.test(f)).map((f) => `scripts/${f}`),
  ];
  const found = {};
  for (const path of files) {
    const n = source(path).match(/\b(?:createSeaIce|createOcean|createModel|createParallelModel|createGpuModel)\(/g)?.length ?? 0;
    if (n) found[path] = n;
  }
  assert.deepEqual(found, listed);
});
