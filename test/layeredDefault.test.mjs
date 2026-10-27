import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { gunzipSync, constants } from 'node:zlib';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { DEFAULT_SHORTWAVE, SOLAR_CONSTANT } from '../js/physics/radiation.module.js';
import { MOIST_DEFAULTS } from '../js/physics/moist.module.js';
import { withCadence } from '../js/cadence.module.js';
import { DEFAULT_RUNS } from '../js/defaultRun.module.js';
import { DEFAULT_SOIL } from '../js/physics/land.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { PREVIOUS_SHORTWAVE, PREVIOUS_BERGERON, PREVIOUS_RADIATION, PREVIOUS_MOIST, PREVIOUS_LAND } from './helpers/previousDefaults.mjs';

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

test('the page\'s default states record no shortwave scheme and no Bergeron option, so the defaults apply on load and set the drivers\' cadence', () => {
  const root = new URL('..', import.meta.url).pathname;
  for (const [N, manifest] of Object.entries(DEFAULT_RUNS)) {
    const { parts } = JSON.parse(readFileSync(root + manifest, 'utf8'));
    const head = gunzipSync(readFileSync(root + manifest.replace(/[^/]*$/, '') + parts[0].file), { finishFlush: constants.Z_SYNC_FLUSH });
    assert.equal(String.fromCharCode(...head.subarray(0, 4)), 'GCMS');
    const { options = {} } = JSON.parse(new TextDecoder().decode(head.subarray(8, 8 + head.readUInt32LE(4))));
    assert.equal(options.radiation?.shortwave, undefined, `${manifest} records no shortwave scheme`);
    assert.equal(options.moist?.bergeronAutoconversion, undefined, `${manifest} records no Bergeron option`);
    assert.equal(options.land?.soil, undefined, `${manifest} records no soil option`);
    const every = withCadence({ radiation: options.radiation ?? {} }, 1350 * 16 / Number(N), Number(N)).radiation.radiationEvery;
    assert.equal(every, { 64: 8, 128: 16, 192: 24 }[N], `the page's radiation cadence at N = ${N}`);
  }
});

test('land.soil defaults to the layers and the previous default is the slab: a model built without the option has the layers and one asked for the previous option has the slab', () => {
  assert.equal(DEFAULT_SOIL, 'layers');
  assert.deepEqual(PREVIOUS_LAND, { soil: 'slab' });
  const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 ? 300 : -4000));
  const plain = createModel(new Grid(16), { ocean: false, topography }), previous = createModel(new Grid(16), { ocean: false, topography, land: PREVIOUS_LAND });
  assert.deepEqual([plain.land.layered, plain.land.deepNames], [true, ['soilT2', 'soilT3', 'soilT4']]);
  assert.deepEqual([previous.land.layered, previous.land.deepNames], [false, []]);
});

test('the GPU engine takes the same defaults: PHYSICS_DEFAULTS names the layered shortwave and the Bergeron factor, and a model asked for the previous options keeps them', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { PHYSICS_DEFAULTS } = await import('../js/gpu/core.gpu.js');
  assert.equal(PHYSICS_DEFAULTS.shortwave, DEFAULT_SHORTWAVE);
  assert.equal(PHYSICS_DEFAULTS.bergeronAutoconversion, true);
  const { createGpuModel } = await import('../js/gpu/model.gpu.js');
  const plain = await createGpuModel(new Grid(16), { ocean: false }), previous = await createGpuModel(new Grid(16), { ocean: false, radiation: PREVIOUS_SHORTWAVE, moist: PREVIOUS_BERGERON });
  assert.deepEqual([plain.gpu.physics.shortwave, plain.gpu.physics.bergeronAutoconversion], ['layered', true]);
  assert.deepEqual([previous.gpu.physics.shortwave, previous.gpu.physics.bergeronAutoconversion], ['sequential', false]);
  const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 ? 300 : -4000));
  const land = await createGpuModel(new Grid(16), { ocean: false, topography }), slab = await createGpuModel(new Grid(16), { ocean: false, topography, land: PREVIOUS_LAND });
  assert.deepEqual([land.gpu.physics.soilColumn.depths.length, land.land.deepNames], [4, ['soilT2', 'soilT3', 'soilT4']]);
  assert.deepEqual([slab.gpu.physics.soilColumn, slab.land.deepNames, slab.gpu.layout.PH.SOILT2], [null, [], undefined]);
});
