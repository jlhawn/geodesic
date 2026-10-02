import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, writeFileSync, mkdtempSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { spawnSync } from 'node:child_process';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { saturationHumidity } from '../js/physics/moist.module.js';
import { STEFAN_BOLTZMANN } from '../js/physics/radiation.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { encodeState } from '../js/stateFile.module.js';

function partlyCloudy() {
  const model = createModel(new Grid(6), { ocean: false, divergenceDamping: 0, radiation: { mixedLayerDeck: false } });
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  model.step(900);
  const [pi, theta, , , q, qc] = model.state, { K } = model.core, C = model.mesh.nCells;
  model.core.diagnose(pi, theta, q, qc);
  const { exnerLayer, sigmaMid } = model.core.diagnostics;
  for (let i = 0; i < C; i++) {
    for (let k = 0; k < K; k++) {
      if (!(sigmaMid[k] > 0.3 && sigmaMid[k] < 0.9) || (k + i) % 5 === 0) continue;
      const x = k * C + i, qs = saturationHumidity(theta[x] * exnerLayer[x], pi[i] * sigmaMid[k]);
      q[x] = (0.85 + 0.01 * ((3 * i + k) % 11)) * qs;
      qc[x] = (0.002 + 0.05 * (((i + 2 * k) % 9) / 8)) * qs;
    }
  }
  model.core.diagnose(pi, theta, q, qc);
  return model;
}

function columns(model, radiation) {
  const run = createModel(model.mesh, { ocean: false, radiation: { mixedLayerDeck: false, clearSkyPass: true, ...radiation } });
  run.state.forEach((a, n) => a.set(model.state[n]));
  const [pi, theta, , surfaceT, q, qc] = run.state, { K } = run.core, C = run.mesh.nCells;
  run.core.diagnose(pi, theta, q, qc);
  run.radiation.setTime(model.time);
  const rad = run.radiation, out = { outgoing: new Float64Array(C), back: new Float64Array(C), clear: new Float64Array(C), longwave: new Float64Array(K * C), closure: 0, scale: 0 };
  for (let i = 0; i < C; i++) {
    rad.column(i, pi[i], theta, surfaceT[i], 5, undefined, rad.insolation(i), q[(K - 1) * C + i], q, qc, 0.07, 0.07);
    out.outgoing[i] = rad.budget.outgoingLongwave; out.back[i] = rad.budget.downwardLongwave; out.clear[i] = rad.budget.clearOutgoingLongwave;
    let sum = 0;
    for (let k = 0; k < K; k++) { out.longwave[k * C + i] = rad.longwave[k * C + i]; sum += rad.longwave[k * C + i]; }
    const difference = STEFAN_BOLTZMANN * surfaceT[i] ** 4 - out.back[i] - out.outgoing[i];
    out.closure = Math.max(out.closure, Math.abs(sum - difference));
    out.scale = Math.max(out.scale, out.outgoing[i]);
  }
  return out;
}

test('under exponential-random overlap in the longwave the layers\' longwave heating of every column sums to the surface\'s emission less the back radiation and the OLR, as under random overlap; the overlap raises the OLR and lowers the back radiation of partly cloudy columns, and leaves the clear-sky OLR alone', () => {
  const model = partlyCloudy();
  const exponential = columns(model, {}), random = columns(model, { longwaveOverlap: 'random' });
  const C = model.mesh.nCells;
  let raised = 0, lowered = 0, moved = 0;
  for (let i = 0; i < C; i++) {
    raised += exponential.outgoing[i] - random.outgoing[i]; lowered += random.back[i] - exponential.back[i];
    moved = Math.max(moved, Math.abs(exponential.outgoing[i] - random.outgoing[i]));
    assert.equal(exponential.clear[i], random.clear[i], `the clear-sky OLR of cell ${i}`);
  }
  console.log(`${C} partly cloudy columns at N=6: the layers' longwave closes on the flux difference to ${exponential.closure.toExponential(1)} W/m2 (random ${random.closure.toExponential(1)}); exponential-random less random OLR ${(raised / C).toFixed(3)}, back radiation ${(-lowered / C).toFixed(3)} W/m2 in the mean, at most ${moved.toFixed(2)} in a column`);
  assert.ok(exponential.closure < 1e-11 * exponential.scale && random.closure < 1e-11 * random.scale, `closure ${exponential.closure} and ${random.closure} W/m2`);
  assert.ok(raised > 0 && lowered > 0 && moved > 0.5, `the overlap raises the OLR by ${raised / C} and lowers the back radiation by ${lowered / C} W/m2, at most ${moved}`);
});

test('with a decorrelation length so short that adjacent cloudy layers overlap at random (α = 0), the two-region longwave gives the random overlap\'s OLR, back radiation and layer heating to rounding', () => {
  const model = partlyCloudy(), short = { decorrelationLength: 1e-6, decorrelationSlope: 0 };
  const chain = columns(model, short), random = columns(model, { ...short, longwaveOverlap: 'random' });
  const C = model.mesh.nCells, { K } = model.core;
  let flux = 0, heating = 0, scale = 0;
  for (let i = 0; i < C; i++) {
    flux = Math.max(flux, Math.abs(chain.outgoing[i] - random.outgoing[i]) / random.outgoing[i], Math.abs(chain.back[i] - random.back[i]) / random.back[i]);
    for (let k = 0; k < K; k++) { heating = Math.max(heating, Math.abs(chain.longwave[k * C + i] - random.longwave[k * C + i])); scale = Math.max(scale, Math.abs(random.longwave[k * C + i])); }
  }
  assert.ok(flux < 1e-13 && heating < 1e-11 * scale, `OLR and back radiation apart by ${flux} relative, layer longwave by ${heating} W/m2 of ${scale}`);
});

test('longwave overlap rejects an unknown overlap', () => {
  const model = createModel(new Grid(4), { ocean: false });
  assert.throws(() => createModel(model.mesh, { ocean: false, radiation: { longwaveOverlap: 'maximum' } }), /longwaveOverlap must be 'exponentialRandom' or 'random'/);
});

test('scripts/longwaveOverlap.mjs repeats the radiation\'s OLR and surface downward longwave under its own overlap, exponential-random by default and random on request, and finds the exponential-random OLR above the random', () => {
  const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
  const model = createModel(new Grid(12), { topography, ocean: false });
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  model.land.initialize();
  for (let n = 0; n < 6; n++) model.step(1800);
  const [pi, theta, u, surfaceT, q, qc, ice] = model.state, { mlmSubsidence, mlmHeight, mlmGate } = model.radiation;
  const dir = mkdtempSync(join(tmpdir(), 'longwaveOverlap-')), file = join(dir, 'overlap12_day0000.bin');
  writeFileSync(file, encodeState({ N: 12, K: model.core.K, day: 0, time: model.time, terrain: true, levels: model.core.levels, pi, theta, u, surfaceT, q, qc, ice, concentration: model.seaIce.concentration, mlmSubsidence, mlmHeight, mlmGate, land: model.land.serialize() }));
  const script = (radiation) => spawnSync(process.execPath, [new URL('../scripts/longwaveOverlap.mjs', import.meta.url).pathname, file], { env: { ...process.env, RADIATION: JSON.stringify(radiation) }, encoding: 'utf8' });
  const runs = [script({}), script({ longwaveOverlap: 'random' })];
  rmSync(dir, { recursive: true, force: true });
  for (const [n, run] of runs.entries()) {
    assert.equal(run.status, 0, run.stderr);
    const [, olr, down] = run.stdout.match(/repeats the radiation's OLR to ([0-9.e+-]+) and its surface downward longwave to ([0-9.e+-]+) relative/);
    assert.ok(Number(olr) < 1e-6 && Number(down) < 1e-6, `${n ? 'random' : 'exponential-random'}: OLR to ${olr}, downward to ${down}`);
    assert.match(run.stdout, n ? /the random recomputation \(the radiation's longwaveOverlap\)/ : /the exponential-random recomputation \(the radiation's longwaveOverlap\)/);
    const [, raised] = run.stdout.match(/^global: .*exponential-random less random: OLR (-?[0-9.]+),/m);
    assert.ok(Number(raised) > 0, `exponential-random less random OLR ${raised} W/m2`);
  }
});
