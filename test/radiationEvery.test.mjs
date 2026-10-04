import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { STEFAN_BOLTZMANN } from '../js/physics/radiation.module.js';

function heldModel(radiationEvery, options = {}) {
  const model = createModel(new Grid(6), { ocean: false, radiation: { stratus: true, clearSkyPass: true, radiationEvery, ...options }, surface: { exchange: 'roughness' } });
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  return model;
}

test('radiationEvery must be a whole number of steps, 1 or more', () => {
  for (const radiationEvery of [0, 2.5, -1, NaN]) assert.throws(() => heldModel(radiationEvery), /radiationEvery/);
});

test('a full call keeps fluxes that close: the layers\' shortwave and the surface\'s make the absorbed sunlight, the layers\' longwave and the surface\'s make the OLR, and a change in the surface\'s emission is absorbed or escapes whole', () => {
  for (const options of [{}, { longwaveOverlap: 'random' }, { longwaveScheme: 'gray' }]) {
    const model = heldModel(4, options), { held, longwave } = model.radiation, { K } = model.core, C = model.mesh.nCells;
    model.step(900);
    let worst = { shortwave: 0, longwave: 0, share: 0 }, lit = 0;
    for (let i = 0; i < C; i++) {
      let shortwave = held.absorbed[i], emitted = held.back[i] - held.emission[i], share = held.escape[i];
      for (let k = 0; k < K; k++) { shortwave += held.shortwave[k * C + i]; emitted += longwave[k * C + i]; share += held.share[k * C + i]; }
      if (held.mu[i] > 0) lit++;
      worst = { shortwave: Math.max(worst.shortwave, Math.abs(shortwave - held.top[i])), longwave: Math.max(worst.longwave, Math.abs(emitted + held.outgoing[i])), share: Math.max(worst.share, Math.abs(share - 1)) };
      assert.ok(held.escape[i] > 0 && held.escape[i] < 1 && held.clearEscape[i] >= held.escape[i] - 1e-12, `cell ${i}: escape ${held.escape[i]}, clear ${held.clearEscape[i]}`);
    }
    console.log(`${JSON.stringify(options)}: ${lit} of ${C} columns lit over the call; worst closure: shortwave ${worst.shortwave.toExponential(1)}, longwave ${worst.longwave.toExponential(1)} W/m², shares ${worst.share.toExponential(1)}`);
    assert.ok(lit > C / 3, `${lit} lit columns`);
    assert.ok(worst.shortwave < 1e-9 && worst.longwave < 1e-9 && worst.share < 1e-12, JSON.stringify(worst));
  }
});

test('between calls the sunlight scales with the cosine of the zenith angle, a cell whose sun rises between calls takes sunlight, and the OLR follows the surface\'s emission', () => {
  const model = heldModel(4), radiation = model.radiation, { held } = radiation, C = model.mesh.nCells, dt = 900;
  model.step(dt);
  const rising = [];
  for (let i = 0; i < C; i++) if (radiation.cosZenith(i) === 0 && held.mu[i] > 0) rising.push(i);
  const before = Float64Array.from(radiation.summed.absorbedSolar), insolationBefore = Float64Array.from(radiation.summed.insolation);
  const top = Float64Array.from(held.top), outgoing = Float64Array.from(held.outgoing), escape = Float64Array.from(held.escape), emission = Float64Array.from(held.emission);
  for (let n = 1; n < 4; n++) {
    const time = model.time, skin = Float64Array.from(model.state[3]), ice = Float64Array.from(model.state[6]);
    model.step(dt);
    radiation.setTime(time);
    for (let i = 0; i < C; i++) {
      const scale = held.mu[i] > 0 ? radiation.cosZenith(i) / held.mu[i] : 0, gained = radiation.summed.absorbedSolar[i] - before[i];
      assert.ok(Math.abs(gained - scale * top[i]) < 1e-9 * Math.max(1, top[i]), `step ${n + 1} cell ${i}: absorbed ${gained} against ${scale} × ${top[i]}`);
      assert.ok(Math.abs(radiation.summed.insolation[i] - insolationBefore[i] - radiation.insolation(i)) < 1e-9, `step ${n + 1} cell ${i}: insolation`);
      if (ice[i] > 0) continue;
      const lift = STEFAN_BOLTZMANN * skin[i] ** 4 - emission[i];
      assert.ok(Math.abs(radiation.outgoing[i] - (outgoing[i] + escape[i] * lift)) < 1e-9, `step ${n + 1} cell ${i}: OLR ${radiation.outgoing[i]} against ${outgoing[i]} + ${escape[i]} × ${lift}`);
    }
    before.set(radiation.summed.absorbedSolar); insolationBefore.set(radiation.summed.insolation);
    radiation.setTime(model.time);
  }
  let sunrise = 0;
  for (const i of rising) if (radiation.cosZenith(i) > 0 || radiation.summed.absorbedSolar[i] > 0) sunrise++;
  console.log(`${rising.length} cells dark at the call whose sun rises before the next take the call's sunlight; ${sunrise} of them lit by the fourth step`);
  assert.ok(rising.length > 0 && sunrise > 0, `${rising.length} rising cells`);
});

test('with radiationEvery 1 the model steps as without the option', () => {
  const plain = createModel(new Grid(6), { ocean: false, radiation: { stratus: true }, surface: { exchange: 'roughness' } });
  initializeState(plain, {}).forEach((values, a) => plain.state[a].set(values));
  const once = heldModel(1, { clearSkyPass: false });
  for (let n = 0; n < 3; n++) { plain.step(900); once.step(900); }
  for (let a = 0; a < plain.state.length; a++) assert.deepEqual(once.state[a], plain.state[a], `state ${a}`);
});
