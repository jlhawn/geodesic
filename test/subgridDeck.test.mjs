import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { sigmaInterfaces, SIGMA_GRIDS } from '../js/dynamics/sigmaCore.module.js';
import { physicsConstants } from '../js/gpu/physics.gpu.js';
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
  return { model, ...profile, C, decked: Array.from({ length: C }, (_, i) => i).filter((i) => model.radiation.mlmTop[i] > 0) };
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

test('with deckCloud \'top\' the deck\'s water sits in the layer that holds its top, the sub-grid inversion or the carried height, and its longwave cooling moves there from σ 0.92', () => {
  const sigma = deckRun('bl60', { deckInversion: 'subgrid' }), top = deckRun('bl60', { deckInversion: 'subgrid', deckCloud: 'top' }), grid = deckRun('bl60', { deckCloud: 'top' });
  const { ambiguous: a, C } = top, stratus = top.model.radiation.stratusLayer;
  for (const i of sigma.decked) assert.equal(sigma.model.radiation.deckLayer[i], stratus, `cell ${i}`);
  for (const i of top.decked) assert.equal(top.model.radiation.deckLayer[i], top.model.radiation.mlmShare[i] < 1 ? a : a + 1, `cell ${i}`);
  const both = top.decked.filter((i) => top.model.radiation.stratus[i] > 0 && sigma.model.radiation.stratus[i] > 0);
  const mean = (run, k) => both.reduce((s, i) => s + run.model.radiation.longwave[k * C + i], 0) / both.length;
  console.log(`over ${both.length} cells clouded both ways, the longwave heating of the inversion's layer ${a} goes from ${mean(sigma, a).toFixed(2)} to ${mean(top, a).toFixed(2)}, of the σ 0.92 layer ${stratus} from ${mean(sigma, stratus).toFixed(2)} to ${mean(top, stratus).toFixed(2)}`);
  assert.ok(both.length > 0.5 * C, `${both.length} cells clouded both ways`);
  assert.ok(mean(top, a) < mean(sigma, a) && mean(top, stratus) > mean(sigma, stratus), 'the cooling stayed at σ 0.92');
  const { geopotential, g } = grid.model.core.diagnostics;
  let away = 0;
  for (const i of grid.decked) {
    const k = grid.model.radiation.deckLayer[i], h = grid.model.radiation.mlmTop[i];
    if (!(h > geopotential[(k + 1) * C + i] / g && h < geopotential[(k - 1) * C + i] / g)) away++;
  }
  console.log(`on the grid ceiling the deck's layer brackets its height on all but ${away} of ${grid.decked.length} cells`);
  assert.ok(grid.decked.length > 0.5 * C && away === 0, `${away} of ${grid.decked.length} away from their layer`);
});

test('the deck options refuse other values and need the mixed-layer deck, in both engines', async () => {
  const { PHYSICS_DEFAULTS } = await import('../js/gpu/core.gpu.js');
  for (const [options, message] of [[{ deckInversion: 'cell' }, /deckInversion must be/], [{ deckCloud: 'layer' }, /deckCloud must be/], [{ deckInversion: 'subgrid', prognosticHeight: false }, /needs the mixed-layer deck/], [{ deckCloud: 'top', mixedLayerDeck: false }, /needs the mixed-layer deck/]]) {
    assert.throws(() => createModel(new Grid(N), { ocean: false, radiation: options }), message, `CPU: ${JSON.stringify(options)}`);
    assert.throws(() => physicsConstants({ ...PHYSICS_DEFAULTS, R: 287, ...options }), message, `GPU: ${JSON.stringify(options)}`);
  }
});
