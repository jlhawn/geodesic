import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { shelfBands } from './helpers/iceShelf.mjs';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};
const { readRanges } = gpuAvailable ? await import('../js/gpu/device.module.js') : {};

/*
 * The GPU ocean packed onto the sea's cells, edges and vertices
 * (js/gpu/oceanPacking.gpu.js) against the same ocean with every row kept
 * (ocean option packed: false), from one state and in lockstep. The packed
 * ocean runs the same arithmetic on the same values, so the ocean's stores,
 * everything it saves and restores, the page's frames and the atmosphere's S
 * and PH agree to the bit after two simulated days with a continent, the ice
 * drift on, and with ice shelves melting from beneath and calving.
 * One GPU model at a time (CLAUDE.md).
 */
const N = 16, DT = 900, STEPS = 192;
const continent = syntheticTopography(90, 180, (lat, lon) => ((Math.cos(lon) > 0 && Math.abs(lat) < 1.2) || lat < -1.15 ? 300 + 2500 * Math.exp(-(((lat - 0.3) / 0.3) ** 2)) : -4000));
const shelves = shelfBands();

async function prepared(topography, geography, ocean) {
  const model = await createGpuModel(new Grid(N), { topography, geography, ocean: { everySteps: 2, ...ocean }, land: { growthTime: 3 * 3600, declineTime: 2 * 3600, snowDeclineTime: 4 * 3600 } });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i]) model.state[6][i] = 0;
  model.load();
  model.ocean.initialize(model.state[3], model.state[6]);
  model.land.initialize();
  const C = model.mesh.nCells, soil = new Float64Array(C), snow = new Float64Array(C);
  for (let i = 0; i < C; i++) if (model.geography.land[i]) { soil[i] = 40; snow[i] = model.mesh.latCell[i] > 1.0 ? 5 : 0; }
  model.land.load({ soil, snow, vegetation: new Float64Array(C).fill(0.5) });
  return model;
}

function differing(label, a, b) {
  assert.equal(a.length, b.length, `${label}: length`);
  const x = new Uint32Array(Float32Array.from(a).buffer), y = new Uint32Array(Float32Array.from(b).buffer);
  let count = 0, first = -1;
  for (let n = 0; n < x.length; n++) if (x[n] !== y[n] && !(Number.isNaN(a[n]) && Number.isNaN(b[n]))) { count++; if (first < 0) first = n; }
  assert.equal(count, 0, `${label}: ${count} of ${x.length} differ, the first at ${first} (${a[first]} against ${b[first]})`);
}

async function compare(whole, packed) {
  await Promise.all([whole.settle(), packed.settle()]);
  const [w, p] = [whole.oceanEngine, packed.oceanEngine];
  assert.ok(p.packing.packed && !w.packing.packed);
  const stores = await Promise.all([w.readStores(), p.readStores()]);
  for (const name of ['h', 'u', 'q', 'w']) differing(`store ${name}`, stores[0][name], stores[1][name]);
  const downloads = await Promise.all([w.download(), p.download()]);
  for (const name of ['h', 'u', 'T', 'S', 'eta', 'h1', 'thermoclineDepth']) differing(`download ${name}`, downloads[0][name], downloads[1][name]);
  const restarts = await Promise.all([w.restartArrays(), p.restartArrays()]);
  for (const name of Object.keys(restarts[0])) differing(`restart ${name}`, restarts[0][name], restarts[1][name]);
  const fields = ['sst', 'sss', 'layerDepth', 'thermocline', 'ssh', 'currents', 'current', 'upwelling'];
  const frames = await Promise.all([w, p].map((o, n) => o.frame({ fields, diagnostics: true, depth: 0, upwellingDepth: 100, time: [whole, packed][n].time })));
  for (const name of fields) differing(`frame ${name}`, frames[0].fields[name], frames[1].fields[name]);
  for (const name of ['S', 'PH']) {
    const views = await Promise.all([whole, packed].map((m) => readRanges(m.gpu.device, m.gpu.buffers[name], [{ offset: 0, length: m.gpu.layout[name].total }])));
    differing(`atmosphere ${name}`, views[0][0], views[1][0]);
  }
  for (const name of ['UPWM', 'STRESS', 'T0', 'ETA', 'MELT', 'DRAW']) {
    const sections = await Promise.all([w.readSection(name), p.readSection(name)]);
    differing(`section ${name}`, sections[0], sections[1]);
  }
}

for (const [label, topography, geography, ocean] of [
  ['a continent with the ice drift on', continent, {}, {}],
  ['ice shelves with a cavity beneath', shelves.raster, { ...shelves.options, iceShelves: 'cavity' }, { iceDrift: false }],
]) {
  test(`the packed GPU ocean keeps every field of the unpacked one to the bit, ${label}`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
    const whole = await prepared(topography, geography, { ...ocean, packed: false }), packed = await prepared(topography, geography, ocean);
    assert.ok(packed.oceanEngine.packing.C < whole.oceanEngine.packing.C, 'the packed ocean keeps fewer cells');
    await compare(whole, packed);
    for (let n = 0; n < STEPS; n++) {
      await whole.step(DT); await packed.step(DT);
      if (n % 48 === 47) await compare(whole, packed);
    }
    whole.destroy(); packed.destroy();
  });
}
