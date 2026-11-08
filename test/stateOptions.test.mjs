import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { topographyFromInt16, surfaceGeopotential, terrainOptions, terrainKey, terrainCache, rebalanceSurfacePressure } from '../js/geography.module.js';
import { encodeState, decodeState } from '../js/stateFile.module.js';
import { recordedOptions, givenOptions, headerOptions, differingGroups, rebuiltOptions, runTerrain, moveOntoTerrain, requireOwnSubgrid, stateTerrain } from '../scripts/stateOptions.mjs';

const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);

test("a spin-up's header keeps `options`, which the page reads, to the radiation, the moist physics and the land's soil whatever other groups the run is given, and records every group in `given`", () => {
  const land = { soil: 'layers', coverGoal: 'observed', treeMoisture: false }, groups = { dynamics: { verticalAdvection: 'vanLeerWater' }, orography: false, terrain: { conserveVolume: true }, ocean: { mixing: 'kpp' } };
  const plain = headerOptions({ radiation: {}, moist: {}, land }), wide = headerOptions({ radiation: {}, moist: {}, land, ...groups });
  assert.deepEqual(wide.options, plain.options);
  assert.deepEqual(wide.options, recordedOptions({ radiation: {}, moist: {}, land: { soil: 'layers' } }));
  assert.deepEqual(Object.keys(wide.options), ['radiation', 'moist', 'land']);
  assert.deepEqual(wide.options.land, { soil: 'layers' }, 'the page sees no coverGoal');
  assert.deepEqual(Object.keys(wide.given), ['radiation', 'moist', 'land', 'dynamics', 'orography', 'terrain', 'ocean']);
  assert.deepEqual(wide.given.dynamics, { verticalAdvection: 'vanLeerWater' });
  assert.equal(wide.given.orography, false);
  assert.deepEqual(wide.given.terrain, { conserveVolume: true });
  assert.deepEqual(wide.given.land, land);
  assert.deepEqual(wide.given.radiation, wide.options.radiation);
  assert.deepEqual(Object.keys(recordedOptions(groups)), ['radiation', 'moist', 'land']);
  assert.deepEqual(givenOptions({ land: { soil: 'slab', treeMoisture: false } }).land, { soil: 'slab', treeMoisture: false });
});

test('differingGroups names the groups whose options differ, in any key order, a missing group counting', () => {
  assert.deepEqual(differingGroups({ dynamics: { a: 1, b: 2 } }, { dynamics: { b: 2, a: 1 } }), []);
  assert.deepEqual(differingGroups({ dynamics: { a: 1 }, moist: {} }, { moist: {} }), ['dynamics']);
  assert.deepEqual(differingGroups({ terrain: { conserveVolume: true } }, { terrain: { conserveVolume: false } }), ['terrain']);
});

test('a run that rebuilds a saved state takes the groups its header records under given, each group it is given over the recorded one key by key and false in its place, and refuses a header without them unless a stand-in is given, and a stand-in beside a header that records them unless the two agree', () => {
  const { options, given } = headerOptions({ radiation: { aerosol: 'flat' }, moist: { pcapeMinimum: 0 }, land: { soil: 'slab' }, dynamics: { verticalAdvection: 'vanLeerWater' }, orography: false, gravityWaves: {} });
  const header = { options, given };
  assert.deepEqual(rebuiltOptions(header), given);
  const rebuilt = rebuiltOptions(header, { moist: { bergeronAutoconversion: false }, orography: { constants: 'cy41r2' }, gravityWaves: false, ocean: { mixing: 'kpp' } });
  assert.deepEqual(rebuilt.moist, { ...given.moist, bergeronAutoconversion: false });
  assert.equal(rebuilt.moist.pcapeMinimum, 0);
  assert.deepEqual(rebuilt.orography, { constants: 'cy41r2' });
  assert.equal(rebuilt.gravityWaves, false);
  assert.deepEqual(rebuilt.ocean, { mixing: 'kpp' });
  assert.deepEqual([rebuilt.radiation, rebuilt.land, rebuilt.dynamics], [given.radiation, given.land, given.dynamics]);
  assert.throws(() => rebuiltOptions({ options }, { moist: {} }), /records no given options/);
  assert.deepEqual(rebuiltOptions({ options }, { moist: { pcapeMinimum: 0.05 } }, {}), { moist: { pcapeMinimum: 0.05 } });
  assert.throws(() => rebuiltOptions(header, {}, { moist: {} }), /this header records them, and GIVEN differs from them in radiation, moist, land, dynamics, orography, gravityWaves:/);
  assert.throws(() => rebuiltOptions(header, {}, {}), /GIVEN differs from them in/);
  const agreeing = Object.fromEntries(Object.entries(given).reverse().map(([group, values]) => [group, values && Object.fromEntries(Object.entries(values).reverse())]));
  assert.deepEqual(rebuiltOptions(header, { moist: { pcapeMinimum: 0.05 } }, agreeing), { ...given, moist: { ...given.moist, pcapeMinimum: 0.05 } });
});

test("a state saved over terrain of other options keeps its options, given and terrain options through the file and resumes on that terrain through spinup.mjs's runTerrain and moveOntoTerrain, with TERRAIN or without, its surface pressure unmoved; a state over the default terrain moves onto it once", async () => {
  const terrain = { conserveVolume: true }, model = createModel(new Grid(16), { topography, terrain, physics: false }), { core, mesh, geography } = model;
  const [pi, theta] = initializeState(model, {});
  const save = (pressure, records) => decodeState(encodeState({ N: 16, K: core.K, day: 0, time: 0, terrain: true, ...records, pi: pressure, theta }));
  const resume = (saved, given) => {
    const run = runTerrain(given, saved) ?? true, state = [Float64Array.from(saved.pi), Float64Array.from(saved.theta)];
    return { run, from: moveOntoTerrain(saved, run, { core, mesh, geography }, state, model.surfaceGeopotential), moved: state[0] };
  };
  const records = { terrainOptions: terrain, ...headerOptions({ radiation: {}, moist: {}, land: { soil: 'layers', coverGoal: 'observed' }, dynamics: { verticalAdvection: 'vanLeerWater' }, terrain }) };
  const kept = await save(pi, records);
  assert.deepEqual(kept.options, records.options);
  assert.deepEqual(kept.given, records.given);
  assert.deepEqual(kept.terrainOptions, terrain);
  for (const given of [undefined, terrain, { conserveVolume: true, passes: 2 }]) {
    const { run, from, moved } = resume(kept, given);
    assert.equal(terrainKey(run), terrainKey(terrain), `TERRAIN ${JSON.stringify(given)}`);
    assert.equal(from, null, `TERRAIN ${JSON.stringify(given)}`);
    assert.deepEqual(Array.from(moved), Array.from(kept.pi));
  }
  const once = resume(await save(pi, headerOptions({ radiation: {}, moist: {}, land: {} })), terrain);
  assert.equal(once.from, true);
  assert.ok(once.moved.some((p, i) => p !== Math.fround(pi[i])), 'the default terrain\'s state moves onto conserveVolume\'s surface');
  const again = resume(await save(once.moved, records), undefined);
  assert.equal(again.from, null);
  assert.deepEqual(Array.from(again.moved), Array.from(once.moved, Math.fround));
});

test('a state saved without terrain resumes onto terrain with its surface pressure moved up from sea level, one saved over terrain resumes without terrain with it moved down to sea level, and a state stays where it stands when the run matches it', async () => {
  const model = createModel(new Grid(16), { topography, physics: false }), { core, mesh, geography } = model;
  const [pi, theta] = initializeState(model, {});
  const [flat, raised] = await Promise.all([false, true].map((terrain) => decodeState(encodeState({ N: 16, K: core.K, day: 0, time: 0, terrain, pi, theta }))));
  assert.deepEqual([stateTerrain(flat), stateTerrain(raised)], [false, true]);
  const resume = (saved, terrain, toPhi) => {
    const state = [Float64Array.from(saved.pi), Float64Array.from(saved.theta)];
    return { from: moveOntoTerrain(saved, terrain, { core, mesh, geography }, state, toPhi), moved: state[0] };
  };
  const rebalanced = (saved, fromPhi, toPhi) => {
    const p = Float64Array.from(saved.pi), t = Float64Array.from(saved.theta);
    core.diagnose(p, t, null, null);
    rebalanceSurfacePressure(core, p, t, fromPhi, toPhi);
    return p;
  };
  const up = resume(flat, true, model.surfaceGeopotential), down = resume(raised, false, null);
  assert.equal(up.from, false);
  assert.deepEqual(up.moved, rebalanced(flat, null, model.surfaceGeopotential));
  assert.ok(up.moved.some((p, i) => p !== flat.pi[i]), 'the pressure moves up onto the terrain');
  assert.equal(down.from, true);
  assert.deepEqual(down.moved, rebalanced(raised, model.surfaceGeopotential, null));
  for (const [saved, terrain, toPhi] of [[flat, false, null], [raised, true, model.surfaceGeopotential]]) {
    const kept = resume(saved, terrain, toPhi);
    assert.equal(kept.from, null);
    assert.deepEqual(kept.moved, Float64Array.from(saved.pi));
  }
});

test('a run over terrain of other options than the default refuses to start with the fine grids out of reach, and starts with them in reach, on the default terrain or without orography', () => {
  const before = process.env.TERRAIN_CACHE, other = [{ conserveVolume: true }, { passes: 1 }];
  try {
    if (terrainCache()) for (const terrain of other) assert.doesNotThrow(() => requireOwnSubgrid(terrain, {}));
    process.env.TERRAIN_CACHE = 'none';
    for (const terrain of other) assert.throws(() => requireOwnSubgrid(terrain, {}), /rebuilt below its own surface/);
    for (const [terrain, orography] of [[undefined, {}], [true, {}], [{ passes: 2 }, {}], [{ conserveVolume: false }, {}], [false, {}], [{ conserveVolume: true }, false]]) assert.doesNotThrow(() => requireOwnSubgrid(terrain, orography));
  } finally {
    if (before === undefined) delete process.env.TERRAIN_CACHE; else process.env.TERRAIN_CACHE = before;
  }
});
