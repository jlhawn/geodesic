import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { DEFAULT_SHORTWAVE, SOLAR_CONSTANT } from '../js/physics/radiation.module.js';
import { MOIST_DEFAULTS } from '../js/physics/moist.module.js';
import { withCadence } from '../js/cadence.module.js';
import { PREVIOUS_SHORTWAVE, PREVIOUS_BERGERON, PREVIOUS_RADIATION, PREVIOUS_MOIST } from './helpers/previousDefaults.mjs';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }

function cloudyColumn(radiation) {
  const model = createModel(new Grid(16), { ocean: false, radiation: { cloudCover: 'overcast', cloudColumns: false, clearSkyPass: true, ...radiation } });
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  const { core, radiation: rad, mesh } = model, { K, g, dSigma } = core.diagnostics, C = mesh.nCells, [pi, theta] = model.state;
  const i = 0, k = Math.floor(K / 2), dry = new Float64Array(K * C), qc = new Float64Array(K * C);
  core.diagnoseColumn(i, pi, theta, dry, dry);
  qc[k * C + i] = 0.1 * g / (pi[i] * dSigma[k]);
  core.diagnoseColumn(i, pi, theta, dry, qc);
  const mu = 0.5;
  rad.column(i, pi[i], theta, 290, 5, undefined, { beam: SOLAR_CONSTANT * mu, mu }, 0, dry, qc, 0.1, 0.1);
  return { scheme: rad.shortwave, absorbed: rad.budget.absorbedSolar, reflected: rad.budget.reflectedSolar, atmosphere: rad.budget.atmosphereSolar };
}

function coldLiquid(moist) {
  const model = createModel(new Grid(16), { ocean: false, moist });
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  const { core, moist: physics } = model, { K, C, sigmaMid, exnerLayer } = core.diagnostics, [pi, theta, , , q, qc] = model.state;
  let k = 0;
  for (let j = 0; j < K; j++) if (Math.abs(sigmaMid[j] - 0.45) < Math.abs(sigmaMid[k] - 0.45)) k = j;
  theta[k * C] = 255 / exnerLayer[k * C];
  qc[k * C] = 3e-4;
  core.diagnoseColumn(0, pi, theta, q, qc);
  physics.autoconvertColumn(0, pi, theta, q, qc, 600);
  return qc[k * C];
}

test('the shortwave defaults to the layered scheme and moist.bergeronAutoconversion to true, and the previous defaults are the sequential scheme and false', () => {
  assert.equal(DEFAULT_SHORTWAVE, 'layered');
  assert.equal(MOIST_DEFAULTS.bergeronAutoconversion, true);
  assert.deepEqual([PREVIOUS_SHORTWAVE, PREVIOUS_BERGERON], [{ shortwave: 'sequential' }, { bergeronAutoconversion: false }]);
  assert.equal(PREVIOUS_RADIATION.shortwave, 'sequential');
  assert.equal(PREVIOUS_MOIST.bergeronAutoconversion, false);
  const plain = cloudyColumn({}), layered = cloudyColumn({ shortwave: 'layered' }), sequential = cloudyColumn(PREVIOUS_SHORTWAVE);
  assert.deepEqual([plain.scheme, layered.scheme, sequential.scheme], ['layered', 'layered', 'sequential']);
  assert.deepEqual(plain, layered, 'the default is the layered scheme');
  assert.ok(Math.abs(plain.reflected - sequential.reflected) > 1 && Math.abs(plain.atmosphere - sequential.atmosphere) > 0.1, `the schemes differ over a cloud: reflected ${plain.reflected} against ${sequential.reflected}, atmosphere ${plain.atmosphere} against ${sequential.atmosphere}`);
});

test('the drivers hold the default radiation 45 minutes and the sequential one 11.25 from N = 64 up, and a model built with the previous options is the sequential, Bergeron-free one', () => {
  const dt = (N) => 1350 * 16 / N;
  for (const [N, layered, sequential] of [[128, 16, 4], [64, 8, 2]]) {
    assert.equal(withCadence({}, dt(N), N).radiation.radiationEvery, layered);
    assert.equal(withCadence({ radiation: PREVIOUS_SHORTWAVE }, dt(N), N).radiation.radiationEvery, sequential);
    assert.equal(withCadence({ radiation: PREVIOUS_RADIATION }, dt(N), N).radiation.radiationEvery, sequential);
  }
  assert.deepEqual([32, 16].map((N) => withCadence({}, dt(N), N).radiation.radiationEvery), [1, 1]);
  const previous = createModel(new Grid(16), { ocean: false, radiation: PREVIOUS_SHORTWAVE, moist: PREVIOUS_BERGERON });
  assert.equal(previous.radiation.shortwave, 'sequential');
  const [kept, previousKept, forced] = [coldLiquid({}), coldLiquid(PREVIOUS_BERGERON), coldLiquid({ bergeronAutoconversion: true })];
  assert.equal(kept, forced, 'the default converts the cold liquid with the factor');
  assert.ok(kept < previousKept, `the cold liquid left in the layer, ${kept} against ${previousKept} with the factor off`);
});

test('the GPU engine takes the same defaults: PHYSICS_DEFAULTS names the layered shortwave and the Bergeron factor, and a model asked for the previous options keeps them', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { PHYSICS_DEFAULTS } = await import('../js/gpu/core.gpu.js');
  assert.equal(PHYSICS_DEFAULTS.shortwave, DEFAULT_SHORTWAVE);
  assert.equal(PHYSICS_DEFAULTS.bergeronAutoconversion, true);
  const { createGpuModel } = await import('../js/gpu/model.gpu.js');
  const plain = await createGpuModel(new Grid(16), { ocean: false }), previous = await createGpuModel(new Grid(16), { ocean: false, radiation: PREVIOUS_SHORTWAVE, moist: PREVIOUS_BERGERON });
  assert.deepEqual([plain.gpu.physics.shortwave, plain.gpu.physics.bergeronAutoconversion], ['layered', true]);
  assert.deepEqual([previous.gpu.physics.shortwave, previous.gpu.physics.bergeronAutoconversion], ['sequential', false]);
});
