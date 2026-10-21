import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createModel } from '../js/model.module.js';
import { createSigmaCore, sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { createRadiation } from '../js/physics/radiation.module.js';
import { cloudSaturation } from '../js/physics/moist.module.js';
import { makeColumn, readColumn, K, LEVELS } from '../explain/radiativeColumn.module.js';
import { createSweep, zeroCrossings, uniformCloud, cloudFromHumidity, LAYERS, PATHS, SUNS, HIGHEST } from '../explain/cloudEffects.module.js';

const MLS = readColumn(JSON.parse(readFileSync(new URL('../explain/data/mls_bl36.json', import.meta.url))));
const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const cell = (sweep, layer, grams) => LAYERS.indexOf(layer) * sweep.cols + PATHS.findIndex((w) => Math.abs(w * 1000 - grams) < 1e-9);

const sweep = createSweep(MLS);
sweep.advance();

test('the clouds a reader can choose: the 23 layers from the ground to 141 hPa and 30 water paths from 5 to 400 g/m²', () => {
  assert.equal(LAYERS.length, K - HIGHEST);
  assert.deepEqual([LAYERS[0], LAYERS.at(-1)], [K - 1, HIGHEST]);
  close(MLS.pressure[HIGHEST] / 100, 140.9, 0.05, 'the highest cloud\'s layer, hPa');
  assert.equal(PATHS.length, 30);
  for (const grams of [5, 20, 100, 400]) assert.ok(PATHS.some((w) => Math.abs(w * 1000 - grams) < 1e-9), `${grams} g/m² is a stop`);
  for (let c = 1; c < PATHS.length; c++) assert.ok(PATHS[c] > PATHS[c - 1]);
});

test('the sweep holds the pinned clouds of the midlatitude summer under the average sun, and every cell is the column\'s own', () => {
  assert.ok(sweep.done);
  const pinned = [[26, 100, -116.09, 6.0], [22, 100, -114.75, 36.58], [17, 20, -38.15, 79.37], [17, 5, -11.01, 28.88]];
  for (const [layer, grams, sw, lw] of pinned) {
    const n = cell(sweep, layer, grams);
    close(sweep.shortwave[n], sw, 0.005, `SWCRE of ${grams} g/m² at layer ${layer}`);
    close(sweep.longwave[n], lw, 0.005, `LWCRE of ${grams} g/m² at layer ${layer}`);
  }
  const column = makeColumn({ profile: MLS, humidity: MLS.humidity, cloud: { pressure: MLS.pressure[19], path: 0.14, cover: 0.35 } });
  const n = cell(sweep, 19, 140);
  close(column.budget.shortwaveCloudEffect, 0.35 * sweep.shortwave[n], 1e-9, 'a cover of 0.35 scales the sweep\'s SWCRE');
  close(column.budget.longwaveCloudEffect, 0.35 * sweep.longwave[n], 1e-9, 'and its LWCRE');
  assert.ok(sweep.largest() > 150 && sweep.largest() < 200, `the largest corner effect ${sweep.largest()}`);
});

test('low clouds cool whatever their water, thin high clouds warm, and the line of no net effect climbs toward thinner cloud as the cloud sinks', () => {
  const net = Float64Array.from(sweep.shortwave, (sw, n) => sw + sweep.longwave[n]), lines = zeroCrossings(net, sweep.rows, sweep.cols);
  const row = (layer) => LAYERS.indexOf(layer);
  for (const layer of [K - 1, 30, 26, 24]) {
    assert.equal(lines[row(layer)].length, 0, `layer ${layer} has no crossing`);
    for (let c = 0; c < sweep.cols; c++) assert.ok(net[row(layer) * sweep.cols + c] < 0, `layer ${layer} cools at ${PATHS[c] * 1000} g/m²`);
  }
  const at = (layer) => lines[row(layer)][0] * 1000;
  for (const [layer, lo, hi] of [[22, 12, 20], [20, 50, 60], [17, 100, 120], [HIGHEST, 120, 200]]) {
    assert.equal(lines[row(layer)].length, 1, `one crossing at layer ${layer}`);
    assert.ok(at(layer) > lo && at(layer) < hi, `layer ${layer} crosses at ${at(layer).toFixed(1)} g/m²`);
  }
  assert.ok(at(22) < at(20) && at(20) < at(17) && at(17) < at(HIGHEST));
  console.log(`no net effect at ${[22, 20, 17, HIGHEST].map((layer) => `${(MLS.pressure[layer] / 100).toFixed(0)} hPa ${at(layer).toFixed(1)} g/m²`).join(', ')}`);
});

test('at night every cloud warms, and in dry air the low cloud warms three times as much', () => {
  const night = createSweep(MLS);
  night.start({ sun: SUNS.night });
  night.advance();
  for (let n = 0; n < night.rows * night.cols; n++) { assert.equal(night.shortwave[n], 0); assert.ok(night.longwave[n] >= 0); }
  close(night.longwave[cell(night, 26, 100)], sweep.longwave[cell(sweep, 26, 100)], 1e-9, 'the low cloud\'s LWCRE does not depend on the sun');
  const dry = createSweep(MLS);
  dry.start({ humidity: null });
  dry.advance();
  close(dry.longwave[cell(dry, 26, 100)], 18.13, 0.005, 'the low cloud\'s LWCRE in dry air');
  close(dry.longwave[cell(dry, 17, 20)], 129.25, 0.01, 'the high cloud\'s LWCRE in dry air');
});

test('zeroCrossings interpolates in the logarithm of the water path and skips rows without a change of sign', () => {
  const paths = [0.01, 0.1, 1], values = Float64Array.from([3, 1, -1, -2, -1, -3, 2, NaN, -2]);
  const lines = zeroCrossings(values, 3, 3, paths);
  assert.equal(lines[0].length, 1);
  close(Math.log10(lines[0][0]), -0.5, 1e-12, 'midway between 0.1 and 1 in the logarithm');
  assert.equal(lines[1].length, 0);
  assert.equal(lines[2].length, 0);
});

test('uniformCloud is the model\'s own condensation and the radiation\'s own cover, layer by layer', () => {
  const model = createModel(new Grid(16), { ocean: false, levels: sigmaInterfaces('bl36') });
  const { C, sigmaMid, exnerLayer, dSigma, g } = model.core.diagnostics, [pi, theta, , , q, qc] = model.state;
  for (const k of [26, 22, 17, HIGHEST]) for (const humidity of [0.8, 0.95, 1, 1.1, 1.3]) {
    pi.fill(MLS.surfacePressure);
    for (let pass = 0; pass < 3; pass++) {
      model.core.diagnose(pi, theta);
      for (let j = 0; j < K; j++) { theta[j * C] = MLS.temperature[j] / exnerLayer[j * C]; q[j * C] = 0; qc[j * C] = 0; }
    }
    const temperature = theta[k * C] * exnerLayer[k * C], pressure = pi[0] * sigmaMid[k], before = theta[k * C];
    q[k * C] = humidity * cloudSaturation(temperature, pressure, true, 273.15, 235.15, {}).qs;
    model.core.diagnoseColumn(0, pi, theta, q, qc);
    model.moist.condenseColumn(0, pi, theta, q, qc);
    model.radiation.column(0, pi[0], theta, MLS.groundTemperature, 5, undefined, 340.5, null, q, qc, 0.07, 0.06, 1, 0, 0, 0, 0, 0, 'radiation');
    const cloud = uniformCloud({ temperature, pressure, surfacePressure: pi[0], humidity, mass: pi[0] * dSigma[k] / g });
    close(cloud.condensate, qc[k * C], 1e-18, `condensate at layer ${k}, humidity ${humidity}`);
    close(cloud.vapor, q[k * C], 1e-18, `vapor left at layer ${k}, humidity ${humidity}`);
    close(cloud.warming, (theta[k * C] - before) * exnerLayer[k * C], 1e-12, `latent heat at layer ${k}, humidity ${humidity}`);
    if (cloud.condensate > 0) close(cloud.cover, model.radiation.layerCover[k], 1e-12, `cover at layer ${k}, humidity ${humidity}`);
  }
});

test('the uniform distribution: cloud begins at the critical humidity, covers half the cell at saturation and all of it at 2 − RH_c', () => {
  const at = (humidity, temperature = 285, pressure = 850e2, surfacePressure = 1000e2) => uniformCloud({ temperature, pressure, surfacePressure, humidity });
  const { critical } = at(1);
  close(critical, 0.903, 0.0005, 'RH_c at 850 hPa');
  assert.equal(at(critical - 1e-6).condensate, 0);
  assert.ok(at(critical + 1e-3).condensate > 0);
  close(at(1).cover, 0.5, 1e-12, 'saturated on average, half the cell is cloud');
  close(at(0.95).cover, (0.95 - critical) / (2 * (1 - critical)), 1e-12, 'the cover is the share of the cell above saturation');
  close(at(2 - critical).cover, 1, 1e-9, 'overcast from 2 − RH_c');
  close(at(0.95).condensate * 1000, 0.021, 0.0005, 'condensate at 0.95, g/kg');
  close(at(0.95).cover, 0.24, 0.005, 'cover at 0.95');
  close(at(1).condensate * 1000, 0.092, 0.0005, 'condensate at saturation, g/kg');
  const low = MLS.surfacePressure * 0.999;
  close(uniformCloud({ temperature: 294, pressure: low, surfacePressure: MLS.surfacePressure, humidity: 1 }).critical, 0.975, 0.001, 'RH_c at the ground');
  close(uniformCloud({ temperature: 230, pressure: 250e2, surfacePressure: MLS.surfacePressure, humidity: 1 }).critical, 0.75, 1e-6, 'RH_c aloft');
});

test('the cloud a layer\'s humidity makes, blended by its cover, is the model\'s own partial cover: the same infrared, and the same sunlight once the cloud holds more than a few g/m²', () => {
  const mesh = buildMesh(new Grid(16)), C = mesh.nCells, core = createSigmaCore(mesh, { levels: LEVELS }), { exnerLayer } = core.diagnostics, ps = MLS.surfacePressure;
  const pi = new Float64Array(C).fill(ps), theta = new Float64Array(K * C), q = new Float64Array(K * C), qc = new Float64Array(K * C);
  const pdf = createRadiation(mesh, core, { longwaveLines: 'scaled', solarConstant: 340.5 * 1.5, ozoneProfile: MLS.ozone, nearInfrared: 'clirad' });
  pdf.useCondensation(createModel(new Grid(16), { ocean: false, levels: LEVELS }).moist.condensation);
  const column = makeColumn({ profile: MLS, humidity: MLS.humidity });
  const call = (water) => {
    qc.fill(0);
    if (water) qc[water.k * C] = water.condensate;
    pdf.column(0, ps, theta, MLS.groundTemperature, 5, undefined, 340.5, null, q, qc, 0.07, 0.07, 1, 0, 0, 0, 0, 0, 'radiation');
    return { ...pdf.budget, cover: water ? pdf.layerCover[water.k] : 0 };
  };
  for (const [k, humidity] of [[26, 0.95], [26, 1], [26, 1.1], [22, 1], [17, 1], [17, 1.12], [HIGHEST, 1]]) {
    const made = cloudFromHumidity(MLS, k, humidity);
    core.diagnoseColumn(0, pi, theta);
    for (let j = 0; j < K; j++) { theta[j * C] = made.temperature[j] / exnerLayer[j * C]; q[j * C] = made.humidity[j]; }
    core.diagnoseColumn(0, pi, theta);
    const cloudy = call({ k, condensate: made.made.condensate }), clear = call(null);
    column.set({ humidity: made.humidity, temperature: made.temperature, cloud: made.cloud });
    const b = column.budget, what = `layer ${k} at ${humidity}`;
    close(made.cloud.cover, cloudy.cover, 1e-12, `${what}: cover`);
    close(b.olr, cloudy.outgoingLongwave, 1e-9, `${what}: OLR`);
    close(b.dlr, cloudy.downwardLongwave, 1e-9, `${what}: DLR`);
    close(b.clearOlr, clear.outgoingLongwave, 1e-9, `${what}: clear OLR`);
    close(b.clearAbsorbed, clear.absorbedSolar, 1e-9, `${what}: clear sunlight`);
    close(b.absorbed, cloudy.absorbedSolar, made.made.gridPath > 0.01 ? 1e-6 : 0.1, `${what}: sunlight absorbed`);
  }
});
