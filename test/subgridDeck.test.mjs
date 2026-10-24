import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { sigmaInterfaces, SIGMA_GRIDS } from '../js/dynamics/sigmaCore.module.js';
import { physicsConstants } from '../js/gpu/physics.gpu.js';
import { LATENT_HEAT } from '../js/physics/moist.module.js';
import { deckProfile, SHARE } from './helpers/deckProfile.mjs';

const N = 16;
const MOIST = { cloudLifetime: 3 * 3600, plumeCape: 70 };
const DECK = { mixedLayerDeck: true, deckRest: 'depth', minimumInversion: 2 };

function deckRun(levels, radiation = {}, turbulence = 'dry') {
  const model = createModel(new Grid(N), { ocean: false, levels: sigmaInterfaces(levels), radiation: { ...DECK, ...radiation }, moist: MOIST, boundaryLayer: { turbulence }, surface: { exchange: 'roughness' } });
  const profile = deckProfile(model), saved = model.state.map((a) => Float64Array.from(a));
  model.step(900);
  model.state.forEach((a, n) => a.set(saved[n]));
  model.step(900);
  const C = model.mesh.nCells;
  return { model, ...profile, saved, C, decked: Array.from({ length: C }, (_, i) => i).filter((i) => model.radiation.mlmTop[i] > 0) };
}

test('the sub-grid inversion finds the share of free troposphere mixed into the layer above the mixed layer, on bl60, cam26 and bl48 where it is built, and puts the deck\'s top inside that layer', () => {
  for (const [levels, mean] of [['bl60', 0.005], ['cam26', 0.01], ...(SIGMA_GRIDS.includes('bl48') ? [['bl48', 0.01]] : [])]) {
    const { model, ambiguous: a, C, decked } = deckRun(levels, { deckInversion: 'subgrid' });
    const r = model.radiation, { geopotential, g } = model.core.diagnostics;
    let worst = 0, sum = 0, outside = 0;
    for (const i of decked) {
      assert.equal(r.mlmLayer[i], a, `${levels} cell ${i}: inversion in layer ${r.mlmLayer[i]}, the mix is in ${a}`);
      const miss = Math.abs(r.mlmShare[i] - SHARE(i));
      worst = Math.max(worst, miss); sum += miss;
      const upper = 0.5 * (geopotential[(a - 1) * C + i] + geopotential[a * C + i]) / g, lower = 0.5 * (geopotential[a * C + i] + geopotential[(a + 1) * C + i]) / g;
      if (!(r.mlmTop[i] > lower - 5 && r.mlmTop[i] < upper + 5)) outside++;
    }
    console.log(`${levels}: the deck runs on ${decked.length} of ${C} cells, each on an inversion in layer ${a}; the share above it misses the share mixed in (0.1 to 0.9) by ${(sum / decked.length).toFixed(4)} on average, at most ${worst.toFixed(4)}; ${outside} tops outside the layer`);
    assert.ok(decked.length > 0.5 * C, `${levels}: deck on ${decked.length} of ${C} cells`);
    assert.ok(sum / decked.length < mean && worst < 0.03, `${levels}: share missed by ${sum / decked.length} on average, at most ${worst}`);
    assert.equal(outside, 0, `${levels}: ${outside} tops outside layer ${a}`);
  }
});

test('where the deck runs on its sub-grid inversion the boundary layer mixes only below the inversion\'s layer, and the deck\'s entrainment, ρ w_e / m at that layer\'s lower interface, is the only exchange across it, in the dry and the moist scheme', () => {
  for (const turbulence of ['dry', 'moist']) {
    const { model, ambiguous: a, C, decked } = deckRun('bl60', { deckInversion: 'subgrid' }, turbulence);
    const r = model.radiation, bl = model.boundaryLayer;
    let sealed = 0, lightest = Infinity, heaviest = 0;
    for (const i of decked) {
      if (bl.depth[i] > r.mlmTop[i] || !(r.mlmEntrainment[i] > 0)) continue;
      sealed++;
      for (let k = bl.kTop; k < a; k++) assert.equal(bl.mixing[k * C + i], 0, `${turbulence} cell ${i}: interface ${k} above the inversion's layer mixes`);
      const density = bl.mixing[a * C + i] * Math.max(r.mlmShare[i], 1e-3) / r.mlmEntrainment[i];
      lightest = Math.min(lightest, density); heaviest = Math.max(heaviest, density);
      assert.equal(bl.entrainment[i], r.mlmEntrainment[i], `${turbulence} cell ${i}: entrainment ${bl.entrainment[i]} against the deck's ${r.mlmEntrainment[i]}`);
      assert.ok(bl.mixing[(a + 1) * C + i] > 0, `${turbulence} cell ${i}: the profile does not reach the inversion's layer`);
    }
    console.log(`${turbulence}: ${sealed} of ${decked.length} decks sealed (the rest under a deeper boundary layer); the coefficient at the inversion's layer times its share over the deck's w_e, the air's density, from ${lightest.toFixed(3)} to ${heaviest.toFixed(3)} kg/m³`);
    assert.ok(sealed > 0.5 * decked.length, `${turbulence}: ${sealed} of ${decked.length} sealed`);
    assert.ok(lightest > 0.9 && heaviest < 1.3, `${turbulence}: density ${lightest} to ${heaviest}`);
  }
});

test('with deckCloud \'top\' the deck\'s water sits in the layer that holds its top, the sub-grid inversion or the carried height; the inversion\'s layer hands its radiative heating under a deck to the mixed layer below it, so the deck\'s longwave cooling moves there from σ 0.92', () => {
  const sigma = deckRun('bl60', { deckInversion: 'subgrid' }), top = deckRun('bl60', { deckInversion: 'subgrid', deckCloud: 'top' }), grid = deckRun('bl60', { deckCloud: 'top' });
  const { ambiguous: a, C } = top, stratus = top.model.radiation.stratusLayer;
  for (const i of sigma.decked) assert.equal(sigma.model.radiation.deckLayer[i], stratus, `cell ${i}`);
  for (const i of top.decked) assert.equal(top.model.radiation.deckLayer[i], top.model.radiation.mlmShare[i] < 1 ? a : a + 1, `cell ${i}`);
  const both = top.decked.filter((i) => top.model.radiation.stratus[i] > 0 && sigma.model.radiation.stratus[i] > 0);
  for (const run of [sigma, top]) for (const i of both) assert.equal(run.model.radiation.longwave[a * C + i], 0, `cell ${i}: the inversion's layer kept longwave heating`);
  const mean = (run, k) => both.reduce((s, i) => s + run.model.radiation.longwave[k * C + i], 0) / both.length;
  console.log(`over ${both.length} cells clouded both ways, the longwave heating of the mixed layer's top layer ${a + 1} goes from ${mean(sigma, a + 1).toFixed(2)} to ${mean(top, a + 1).toFixed(2)}, of the σ 0.92 layer ${stratus} from ${mean(sigma, stratus).toFixed(2)} to ${mean(top, stratus).toFixed(2)}`);
  assert.ok(both.length > 0.5 * C, `${both.length} cells clouded both ways`);
  assert.ok(mean(top, a + 1) < mean(sigma, a + 1) - 20 && mean(top, stratus) > mean(sigma, stratus), 'the cooling stayed at σ 0.92');
  const { geopotential, g } = grid.model.core.diagnostics;
  let away = 0;
  for (const i of grid.decked) {
    const k = grid.model.radiation.deckLayer[i], h = grid.model.radiation.mlmTop[i];
    if (!(h > geopotential[(k + 1) * C + i] / g && h < geopotential[(k - 1) * C + i] / g)) away++;
  }
  console.log(`on the grid ceiling the deck's layer brackets its height on all but ${away} of ${grid.decked.length} cells`);
  assert.ok(grid.decked.length > 0.5 * C && away === 0, `${away} of ${grid.decked.length} away from their layer`);
});

test('the inversion\'s layer hands its heating down inside its column: each column\'s longwave heating sums to its outgoing longwave less the surface\'s net, and with the radiation held its shortwave heating to the sunlight its air absorbs', () => {
  for (const levels of ['bl60', ...(SIGMA_GRIDS.includes('bl48') ? ['bl48'] : [])]) for (const radiationEvery of [1, 2]) {
    const { model, C } = deckRun(levels, { deckInversion: 'subgrid', deckCloud: 'top', radiationEvery }, 'moist'), r = model.radiation, K = model.core.K;
    let worst = 0, sealed = 0;
    for (let i = 0; i < C; i++) {
      if (r.mlmLayer[i] > 0) sealed++;
      let sum = 0;
      if (radiationEvery === 1) {
        for (let k = 0; k < K; k++) sum += r.longwave[k * C + i];
        worst = Math.max(worst, Math.abs(sum + r.outgoing[i] + r.surfaceFlux[i] - r.surfaceSolar[i] + r.sensibleHeat[i] + LATENT_HEAT * r.evaporation[i]));
      } else {
        for (let k = 0; k < K; k++) sum += r.held.shortwave[k * C + i];
        worst = Math.max(worst, Math.abs(sum - r.held.atmosphere[i]));
      }
    }
    const label = `${levels}, ${radiationEvery === 1 ? 'longwave' : 'held shortwave'}`;
    console.log(`${label}: ${sealed} of ${C} columns sealed, worst column off its boundary fluxes by ${worst.toExponential(2)} W/m²`);
    assert.ok(sealed > 0.5 * C && worst < 1e-9, `${label}: ${sealed} sealed, ${worst} W/m² off`);
  }
});

test('with deckResolved \'above\' the deck replaces the resolved cloud of the layer holding its top and every layer below it inside its fraction: columns without a deck are untouched, each column keeps its energy, and the decked columns carry less water and pass more sunlight', () => {
  const options = { deckInversion: 'subgrid', deckCloud: 'top', cloudColumns: true, cloudSpread: 0.75 };
  for (const radiationEvery of [1, 2]) {
    const added = deckRun('bl48', { ...options, radiationEvery, deckResolved: 'added' }, 'moist'), above = deckRun('bl48', { ...options, radiationEvery, deckResolved: 'above' }, 'moist');
    if (radiationEvery > 1) for (const { model, saved } of [added, above]) { model.state.forEach((a, n) => a.set(saved[n])); model.step(900); }
    const a = added.model.radiation, r = above.model.radiation, C = above.C, K = above.model.core.K;
    let bare = 0, decked = 0, affected = 0, wetter = 0, darker = 0, lighter = 0, gain = 0, worst = 0, fromMiss = 0;
    for (let i = 0; i < C; i++) {
      let sum = 0;
      if (radiationEvery === 1) {
        for (let k = 0; k < K; k++) sum += r.longwave[k * C + i];
        worst = Math.max(worst, Math.abs(sum + r.outgoing[i] + r.surfaceFlux[i] - r.surfaceSolar[i] + r.sensibleHeat[i] + LATENT_HEAT * r.evaporation[i]));
      } else {
        for (let k = 0; k < K; k++) sum += r.held.shortwave[k * C + i];
        worst = Math.max(worst, Math.abs(sum - r.held.atmosphere[i]));
      }
      if (!(r.stratus[i] > 0 && r.stratusFraction[i] > 0)) {
        bare++;
        assert.ok(r.outgoing[i] === a.outgoing[i] && r.surfaceShortwave[i] === a.surfaceShortwave[i] && r.lowWater[i] === a.lowWater[i], `cell ${i}: changed without a deck`);
        for (let k = 0; k < K; k++) assert.equal(r.longwave[k * C + i], a.longwave[k * C + i], `cell ${i} layer ${k}`);
        continue;
      }
      decked++;
      if (r.deckFrom[i] !== r.deckLayer[i]) fromMiss++;
      if (r.lowWater[i] > a.lowWater[i]) wetter++;
      if (r.surfaceShortwave[i] < a.surfaceShortwave[i] - 1e-9) darker++;
      if (!(a.lowWater[i] - r.lowWater[i] > 1e-6 && a.surfaceShortwave[i] > 1)) continue;
      affected++;
      if (r.surfaceShortwave[i] > a.surfaceShortwave[i]) lighter++;
      gain += r.surfaceShortwave[i] - a.surfaceShortwave[i];
    }
    const label = radiationEvery === 1 ? 'radiation every step' : 'radiation held';
    console.log(`${label}: ${decked} decked and ${bare} bare columns; the replaced layers start at the deck's own layer on all but ${fromMiss}; of ${affected} sunlit columns losing water, ${lighter} lighter at the surface, by ${(gain / affected).toFixed(2)} W/m² on average; ${darker} decked columns darker; worst column off its boundary fluxes by ${worst.toExponential(2)} W/m²`);
    assert.ok(decked > 0.5 * C && fromMiss === 0 && wetter === 0, `${label}: ${decked} decked, ${fromMiss} replaced from another layer, ${wetter} wetter`);
    assert.ok(affected > 0.1 * decked && lighter === affected && darker === 0 && worst < 1e-9, `${label}: ${lighter} of ${affected} lighter, ${darker} darker, ${worst} W/m² off`);
  }
});

test('the deck options refuse other values and need the mixed-layer deck, in both engines', async () => {
  const { PHYSICS_DEFAULTS } = await import('../js/gpu/core.gpu.js');
  for (const [options, message] of [[{ deckInversion: 'cell' }, /deckInversion must be/], [{ deckCloud: 'layer' }, /deckCloud must be/], [{ deckResolved: 'twice' }, /deckResolved must be/], [{ deckInversion: 'subgrid', prognosticHeight: false }, /needs the mixed-layer deck/], [{ deckCloud: 'top', mixedLayerDeck: false }, /needs the mixed-layer deck/], [{ deckResolved: 'above', prognosticHeight: false }, /needs the mixed-layer deck/]]) {
    assert.throws(() => createModel(new Grid(N), { ocean: false, radiation: options }), message, `CPU: ${JSON.stringify(options)}`);
    assert.throws(() => physicsConstants({ ...PHYSICS_DEFAULTS, R: 287, ...options }), message, `GPU: ${JSON.stringify(options)}`);
  }
});
