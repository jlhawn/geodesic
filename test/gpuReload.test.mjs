import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { initializeState } from '../js/physics/init.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};
const { readBuffer } = gpuAvailable ? await import('../js/gpu/device.module.js') : {};

const N = 16, DT = 1350 * 16 / N;
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const build = () => createGpuModel(new Grid(N), { topography, radiation: { radiationEvery: 2 }, dragEvery: 2 });

const mirrors = (model) => Object.fromEntries(Object.entries({
  state: model.state, concentration: model.seaIce.concentration, windSpeed: model.surface.windSpeed,
  ...Object.fromEntries(['mlmSubsidence', 'mlmHeight', 'mlmGate', 'meanAbsorbedSolar', 'meanOutgoingLongwave', 'meanPlanetaryAlbedo', 'meanShortwaveCloudEffect', 'meanLongwaveCloudEffect', 'evaporation'].map((name) => [name, model.radiation[name]])),
  ...Object.fromEntries(['convectiveRain', 'largeScaleRain', 'cumulusCover', 'cumulusWater', 'subcloudVirtual'].map((name) => [name, model.moist[name]])),
  ...Object.fromEntries(['depth', 'mixingTop', 'regime', 'buoyancyFlux'].map((name) => [`boundary.${name}`, model.boundaryLayer[name]])),
  ...(model.exchange ? { exchangeWind: model.exchange.wind, ...(model.exchange.fixed ? {} : { exchangeHeat: model.exchange.heat }) } : {}),
}).filter(([, value]) => value));
const copy = (arrays) => Object.fromEntries(Object.entries(arrays).map(([name, value]) => [name, Array.isArray(value) ? value.map((a) => Float64Array.from(a)) : Float64Array.from(value)]));

/*
 * A saved state with every field a load sends to the device, taken from a
 * model that has stepped, so that none of them is a start's value.
 */
async function savedState() {
  const model = await build();
  initializeState(model, { geostrophic: !model.surfaceGeopotential }).forEach((values, a) => model.state[a].set(values));
  for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i]) model.state[6][i] = 0;
  model.load();
  model.ocean.initialize(model.state[3], model.state[6]);
  model.land.initialize();
  for (let n = 0; n < 9; n++) await model.step(DT);
  await model.sync();
  const saved = { mirrors: copy(mirrors(model)), ocean: await model.ocean.serialize({ restart: true }), land: await model.land.serialize() };
  model.destroy();
  return saved;
}

function place(model, saved, loadOptions) {
  for (const [name, target] of Object.entries(mirrors(model))) {
    if (name === 'state') target.forEach((array, a) => array.set(saved.mirrors.state[a]));
    else target.set(saved.mirrors[name]);
  }
  model.load(loadOptions);
  model.ocean.load(saved.ocean, model.state[3], model.state[6]);
  model.land.load(saved.land);
}

async function dump(model) {
  const { device, buffers } = model.gpu, ocean = model.oceanEngine.buffers, out = {};
  for (const [name, buffer] of [...Object.entries(buffers), ...Object.entries(ocean).map(([name, buffer]) => [`ocean ${name}`, buffer])]) {
    if (['P', 'FP', 'PR', 'ocean OF'].includes(name)) continue;
    out[name] = await readBuffer(device, buffer, buffer.size, Uint32Array);
  }
  return out;
}

function assertSame(label, reference, got) {
  assert.deepEqual(Object.keys(got), Object.keys(reference));
  for (const [name, words] of Object.entries(reference)) {
    const other = got[name];
    assert.equal(other.length, words.length, `${label}: ${name} length`);
    let differ = 0, first = -1;
    for (let x = 0; x < words.length; x++) if (words[x] !== other[x]) { differ++; if (first < 0) first = x; }
    assert.equal(differ, 0, `${label}: ${name} differs in ${differ} of ${words.length} words, first at ${first}`);
  }
}

test('a load leaves every device buffer as a fresh model\'s load does, after steps have dirtied them', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const saved = await savedState();
  const model = await build();
  place(model, saved, { ocean: false });
  const fresh = await dump(model);
  for (let n = 0; n < 6; n++) await model.step(DT);
  await model.diagnostics();
  const stepped = await dump(model);
  const dirtied = Object.keys(fresh).filter((name) => stepped[name].some((word, x) => word !== fresh[name][x]));
  for (const name of ['S', 'T', 'K1', 'K4', 'D', 'PH', 'FR', 'ocean S', 'ocean T', 'ocean K1', 'ocean K4', 'ocean OD']) assert.ok(dirtied.includes(name), `six steps and a frame change ${name}`);
  place(model, saved, { ocean: false });
  assertSame('loaded again', fresh, await dump(model));
  for (let n = 0; n < 6; n++) await model.step(DT);
  await model.diagnostics();
  place(model, saved);
  assertSame('loaded again with the ocean initialized first', fresh, await dump(model));
  console.log(`N=${N}: six steps and a frame change ${dirtied.length} of ${Object.keys(fresh).length} buffers (${dirtied.join(', ')}); loading the same state again restores all ${Object.keys(fresh).length} word for word, with or without the ocean's initialize`);
  model.destroy();
});
