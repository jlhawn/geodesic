import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel, STATE_NAMES } from '../js/model.module.js';
import { encodeState, decodeState } from '../js/stateFile.module.js';
import { geographyOptions, readTopography, FRESH_GEOGRAPHY } from '../js/geography.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { restoreSurface, savedIceOptions, resumedIceExchange, resumedOcean, resumedIce } from '../scripts/savedLand.mjs';

const N = 16, DT = 1350, ICE = { meltingIceAlbedo: 0.4, iceMeltingRange: 3 };
const topography = readTopography('etopo2022'), sheets = geographyOptions({ ...FRESH_GEOGRAPHY });
const OPTIONS = { topography, geography: sheets.geography, land: { ...sheets.land }, ocean: false };

const restored = (saved, ice) => {
  const model = createModel(new Grid(N), { ...OPTIONS, ...(ice ? { ice } : {}) });
  STATE_NAMES.forEach((name, a) => model.state[a].set(saved[name]));
  return model;
};

test('a model restored from a save by restoreSurface has the running model\'s surface albedo', async () => {
  const run = createModel(new Grid(N), { ...OPTIONS, ice: ICE });
  const { mesh, state, seaIce, land, geography } = run, C = mesh.nCells;
  initializeState(run, { geostrophic: !run.surfaceGeopotential }).forEach((values, a) => state[a].set(values));
  for (let i = 0; i < C; i++) if (geography.land[i]) state[6][i] = 0;
  land.initialize();
  for (let s = 0; s < 6; s++) run.step(DT);
  assert.ok(state.every((values) => values.every(Number.isFinite)), 'the stepped model stays finite');
  const polar = [];
  for (let i = 0; i < C; i++) if (!geography.land[i] && Math.abs(mesh.latCell[i]) > 1.1) polar.push(i);
  for (const i of polar) {
    state[6][i] = 0.3 + 0.4 * (i % 4); state[3][i] = 268 + (i % 7);
    seaIce.snow[i] = (i % 3) * 0.05; seaIce.snowAlbedo[i] = 0.66 + 0.04 * (i % 4);
  }
  seaIce.load(state[6], Float64Array.from({ length: C }, (_, i) => (state[6][i] > 0 ? 0.7 + 0.04 * (i % 5) : 0)));
  const saved = await decodeState(encodeState({
    N, K: run.core.K, day: 0, time: run.time, given: { ice: ICE },
    ...Object.fromEntries(STATE_NAMES.map((name, a) => [name, state[a]])), concentration: seaIce.concentration, land: land.serialize(),
  }));
  assert.deepEqual(savedIceOptions(saved), { iceExchange: 'instant', ...ICE });

  const back = restored(saved, savedIceOptions(saved));
  assert.deepEqual(restoreSurface(back, saved, () => {}), []);
  const iced = polar.filter((i) => state[6][i] > 0);
  assert.ok(iced.length > 20, `${iced.length} cells keep their ice`);
  let worst = 0;
  for (let i = 0; i < C; i++) for (const mu of [null, 0.2, 0.7]) worst = Math.max(worst, Math.abs(back.surfaceAlbedoAt(i, mu) - run.surfaceAlbedoAt(i, mu)));
  assert.ok(worst < 1e-5, `largest albedo difference ${worst}`);

  const plain = restored(saved, null);
  restoreSurface(plain, saved, () => {});
  const apart = Math.max(...iced.map((i) => Math.abs(plain.diffuseAlbedoAt(i) - run.diffuseAlbedoAt(i))));
  assert.ok(apart > 1e-3, `without the save's ice options the albedo differs by up to ${apart}`);

  const { snowAlbedo, ...withoutSnowAlbedo } = saved.land;
  assert.ok(snowAlbedo);
  assert.deepEqual(restoreSurface(restored(saved, ICE), { ...saved, land: withoutSnowAlbedo }, () => {}), ['land.snowAlbedo']);
});

test('restoreSurface hands a turbulent save\'s held heat and conductance to the sea ice, and savedIceOptions carries the ocean\'s iceExchange', async () => {
  const run = createModel(new Grid(N), { ...OPTIONS, ice: { iceExchange: 'turbulent' } });
  const { mesh, state, seaIce, land } = run, C = mesh.nCells;
  initializeState(run, { geostrophic: !run.surfaceGeopotential }).forEach((values, a) => state[a].set(values));
  land.initialize();
  const iceWater = Float64Array.from({ length: C }, (_, i) => (i % 3) * 1e6), iceConductance = Float64Array.from({ length: C }, (_, i) => 20 + (i % 5));
  const saved = await decodeState(encodeState({
    N, K: run.core.K, day: 0, time: run.time, given: { ocean: { iceExchange: 'turbulent' } },
    ...Object.fromEntries(STATE_NAMES.map((name, a) => [name, state[a]])), concentration: seaIce.concentration, land: land.serialize(), ocean: { iceWater, iceConductance },
  }));
  assert.deepEqual(savedIceOptions(saved), { iceExchange: 'turbulent' });
  const back = restored(saved, savedIceOptions(saved));
  assert.ok(back.seaIce.turbulent);
  assert.deepEqual(restoreSurface(back, saved, () => {}), []);
  assert.deepEqual(Array.from(back.seaIce.waterHeat), Array.from(iceWater));
  assert.deepEqual(Array.from(back.seaIce.conductance), Array.from(iceConductance));
  const { ocean, ...withoutExchange } = saved;
  assert.ok(ocean.iceWater);
  assert.deepEqual(restoreSurface(restored(saved, savedIceOptions(saved)), withoutExchange, () => {}), ['ocean.iceWater', 'ocean.iceConductance']);
});

test('a resumed run whose ocean options name no iceExchange takes the one its save was made under: the one it records, in given or else in options, or the instant exchange where it records none; a save that holds the held heat but records none is refused', () => {
  const turbulent = { mode: 'turbulent', frictionMinimum: 2e-2 }, held = { ocean: { iceWater: new Float64Array(4) } };
  assert.equal(resumedIceExchange(null, {}), undefined);
  assert.equal(resumedIceExchange({ given: { ocean: { iceExchange: turbulent } } }, { iceExchange: 'instant' }), undefined, 'the run\'s own option stands');
  assert.deepEqual(resumedIceExchange({ given: { ocean: { iceExchange: turbulent } }, ...held }, {}), turbulent);
  assert.deepEqual(resumedIceExchange({ options: { ocean: { iceExchange: turbulent } }, ...held }, { closureHours: 3 }), turbulent);
  assert.deepEqual(savedIceOptions({ options: { ocean: { iceExchange: turbulent } } }), { iceExchange: turbulent });
  assert.equal(resumedIceExchange({ given: {} }, {}), 'instant', 'a save that records no exchange was made under the instant one');
  assert.equal(resumedIceExchange({}, {}), 'instant');
  assert.deepEqual(savedIceOptions({ given: { ice: { iceAlbedoPhases: true } } }), { iceExchange: 'instant', iceAlbedoPhases: true });
  assert.throws(() => resumedIceExchange({ given: {}, ...held }, {}), /iceWater/);
  assert.equal(resumedIceExchange({ given: { ocean: { iceExchange: turbulent } } }, false), undefined, 'an ocean given as false takes no exchange');
  assert.equal(resumedOcean({ given: { ocean: { iceExchange: turbulent } } }, false), false);
  assert.deepEqual(resumedIce({ given: { ocean: { iceExchange: turbulent }, ice: { iceAlbedoPhases: true, meltOpening: 0.5 } } }, { meltOpening: 0.25 }), { iceAlbedoPhases: true, meltOpening: 0.25 }, 'the save\'s ice options under the run\'s, the exchange left to the ocean');
});
