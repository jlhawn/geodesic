import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { topographyFromInt16, surfaceGeopotential, terrainOptions, terrainKey, terrainCache, rebalanceSurfacePressure } from '../js/geography.module.js';
import { encodeState, decodeState } from '../js/stateFile.module.js';
import { startOptions } from '../js/snapshots.module.js';
import { DIAGNOSE_FIRST } from '../js/physics/boundaryLayer.module.js';
import { recordedOptions, givenOptions, headerOptions, differingGroups, rebuiltOptions, resumeNotes, runTerrain, moveOntoTerrain, requireOwnSubgrid, stateTerrain } from '../scripts/stateOptions.mjs';

const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);

test("a spin-up's header keeps `options`, which the page reads, to the radiation, the moist physics, the land's soil and the boundary layer whatever other groups the run is given, and records every group in `given`", () => {
  const land = { soil: 'layers', coverGoal: 'observed', treeMoisture: false }, groups = { dynamics: { verticalAdvection: 'vanLeerWater' }, orography: false, terrain: { conserveVolume: true }, ocean: { mixing: 'kpp' } };
  const plain = headerOptions({ radiation: {}, moist: {}, land }), wide = headerOptions({ radiation: {}, moist: {}, land, ...groups });
  assert.deepEqual(wide.options, plain.options);
  assert.deepEqual(wide.options, recordedOptions({ radiation: {}, moist: {}, land: { soil: 'layers' } }));
  assert.deepEqual(Object.keys(wide.options), ['radiation', 'moist', 'land', 'boundaryLayer']);
  assert.deepEqual(wide.options.land, { soil: 'layers' }, 'the page sees no coverGoal');
  assert.deepEqual(Object.keys(wide.given), ['radiation', 'moist', 'land', 'boundaryLayer', 'dynamics', 'orography', 'terrain', 'ocean']);
  assert.deepEqual(wide.given.dynamics, { verticalAdvection: 'vanLeerWater' });
  assert.equal(wide.given.orography, false);
  assert.deepEqual(wide.given.terrain, { conserveVolume: true });
  assert.deepEqual(wide.given.land, land);
  assert.deepEqual(wide.given.radiation, wide.options.radiation);
  assert.deepEqual(Object.keys(recordedOptions(groups)), ['radiation', 'moist', 'land', 'boundaryLayer']);
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
  assert.throws(() => rebuiltOptions(header, {}, { moist: {} }), /this header records them, and GIVEN differs from them in radiation, moist, land, boundaryLayer, dynamics, orography, gravityWaves:/);
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

const blank = () => {
  const model = createModel(new Grid(16), { topography, physics: false }), [pi, theta] = initializeState(model, {});
  return (records) => decodeState(encodeState({ N: 16, K: model.core.K, day: 0, time: 0, terrain: true, ...records, pi, theta }));
};

test('a header records the boundary layer\'s order its run was given in `options` and `given`, for both values, through the file, the default where the run is given none, and false in both for a run given no boundary layer', async () => {
  const save = blank();
  for (const diagnoseFirst of [true, false]) {
    const kept = await save(headerOptions({ boundaryLayer: { diagnoseFirst } }));
    assert.deepEqual([kept.options.boundaryLayer, kept.given.boundaryLayer], [{ diagnoseFirst }, { diagnoseFirst }]);
  }
  const unset = await save(headerOptions({}));
  assert.deepEqual([unset.options.boundaryLayer, unset.given.boundaryLayer], [{ diagnoseFirst: DIAGNOSE_FIRST }, { diagnoseFirst: DIAGNOSE_FIRST }]);
  const off = await save(headerOptions({ boundaryLayer: false }));
  assert.deepEqual([off.options.boundaryLayer, off.given.boundaryLayer], [false, false]);
});

test('a `given` without boundaryLayer.diagnoseFirst reads as the plain order: a year-2 control\'s header rebuilds and records so, as in fixedSst, a resume at the default names boundaryLayer, a GIVEN equal to the record as made stands in, a resume of a state that records the order in its options alone, as a page export does, names boundaryLayer where the run differs from it, and one that records it nowhere says the order it was made in cannot be told, with the order the run\'s boundary layer runs, and nothing for a run without one', async () => {
  const save = blank(), strip = ({ boundaryLayer, ...rest }) => rest;
  const fresh = headerOptions({ land: { coverGoal: 'observed' } }), plainGiven = headerOptions({ land: { coverGoal: 'observed' }, boundaryLayer: { diagnoseFirst: false } }).given;
  const control = await save({ options: strip(fresh.options), given: strip(fresh.given) });
  assert.deepEqual(Object.keys(control.given), ['radiation', 'moist', 'land']);
  const run = rebuiltOptions(control), rerecorded = headerOptions(run);
  assert.deepEqual([run.boundaryLayer, rerecorded.options.boundaryLayer, rerecorded.given.boundaryLayer], [{ diagnoseFirst: false }, { diagnoseFirst: false }, { diagnoseFirst: false }]);
  assert.deepEqual(resumeNotes(fresh.given, control, 'ctl'), ["options other than ctl's header records: boundaryLayer"]);
  assert.deepEqual(resumeNotes(plainGiven, control, 'ctl'), []);
  assert.deepEqual(resumeNotes(fresh.given, await save(fresh), 'c'), []);
  const early = await save({ options: strip(fresh.options) });
  assert.deepEqual([resumeNotes(fresh.given, early, 'seventeen64'), resumeNotes(plainGiven, early, 'seventeen64'), resumeNotes({ ...fresh.given, boundaryLayer: false }, early, 'seventeen64')], [["seventeen64 records no boundary-layer order, so the order it was made in cannot be told: this run's boundary layer runs diagnoseFirst"], ["seventeen64 records no boundary-layer order, so the order it was made in cannot be told: this run's boundary layer runs in the plain order"], []]);
  const exported = await save({ options: fresh.options }), exportedPlain = await save({ options: headerOptions({ boundaryLayer: { diagnoseFirst: false } }).options });
  assert.deepEqual([resumeNotes(fresh.given, exported, 'page'), resumeNotes(plainGiven, exported, 'page')], [[], ["options other than page's header records in its options: boundaryLayer"]]);
  assert.deepEqual([resumeNotes(fresh.given, exportedPlain, 'page'), resumeNotes(plainGiven, exportedPlain, 'page')], [["options other than page's header records in its options: boundaryLayer"], []]);
  assert.deepEqual(rebuiltOptions(control, {}, run), run);
  assert.deepEqual(resumeNotes(fresh.given, null, 'none'), []);
});

test('the page starts a state whose options record diagnoseFirst false in the plain order, one whose options record no boundaryLayer group at today\'s default, records the order it starts, and takes a start message\'s boundaryLayer false over the record', async () => {
  const save = blank();
  const plain = startOptions((await save(headerOptions({ boundaryLayer: { diagnoseFirst: false } }))).options, {}), bare = startOptions((await save({})).options, {});
  assert.deepEqual([plain.boundaryLayer, bare.boundaryLayer], [{ diagnoseFirst: false }, { diagnoseFirst: DIAGNOSE_FIRST }]);
  assert.equal(startOptions((await save(headerOptions({ boundaryLayer: { diagnoseFirst: false } }))).options, { boundaryLayer: false }).boundaryLayer, false);
  assert.equal(createModel(new Grid(16), { topography, ...plain }).boundaryLayer.diagnoseFirst, false);
  assert.equal(createModel(new Grid(16), { topography, ...bare }).boundaryLayer.diagnoseFirst, DIAGNOSE_FIRST);
});
