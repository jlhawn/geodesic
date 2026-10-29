import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createOcean } from '../js/ocean/layered.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { upwellingTimescale, upwellingWeight } from '../js/frames.module.js';
import { layerVelocities, worstDifference } from './helpers/upwelling.mjs';

/*
 * The GPU ocean's upwelling memory in a saved state at N=16 on a sea with a
 * continent and a shelf: a coupled model folds at each ocean step with the
 * weights of its layers' timescales, serialize carries the memory with its
 * folds and the time of the last, an upload restores it and the next fold
 * continues the mean, a CPU ocean loaded from the device's state reads the
 * same field, and the prognostic fields do not see the memory. One GPU model
 * at a time (CLAUDE.md).
 */
let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const N = 16, DT = 900, HYBRID = { vertical: 'hybrid', everySteps: 1, iceDrift: false };
const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 300 + 2500 * Math.exp(-(((lat - 0.3) / 0.3) ** 2)) : -4000));
const skip = !gpuAvailable && 'webgpu not installed';

async function build() {
  const model = await createGpuModel(new Grid(N), { topography, ocean: HYBRID });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i]) model.state[6][i] = 0;
  model.load();
  model.ocean.initialize(model.state[3], model.state[6]);
  model.land.initialize();
  return model;
}
const velocities = (model, saved) => layerVelocities(model.mesh, model.oceanEngine.layers, { h: saved.h, u: saved.u, eta: saved.eta, sillEdge: model.oceanEngine.sillEdge, edgeOcean: model.oceanEngine.edgeOcean }, model.oceanEngine.cellOcean);
function fold(previous, w, h, cellOcean, C, dt, folds) {
  return Float64Array.from(w, (x, n) => {
    const k = Math.floor(n / C), i = n % C;
    if (!cellOcean[i]) return 0;
    let base = 0;
    for (let j = 0; j <= k; j++) base += h[j * C + i];
    const alpha = upwellingWeight(dt, folds, upwellingTimescale(base));
    return alpha === 1 ? x : previous[n] + alpha * (x - previous[n]);
  });
}
const same = (a, b, label) => {
  assert.equal(a.length, b.length, `${label}: lengths`);
  for (let n = 0; n < a.length; n++) if (!Object.is(a[n], b[n])) assert.fail(`${label}: ${a[n]} against ${b[n]} at ${n}`);
};

test('a coupled model folds its ocean at each ocean step with the time the step ends at and the weights of its layers\' timescales', { skip }, async () => {
  const model = await build();
  try {
    const { mesh } = model, C = mesh.nCells, { cellOcean } = model.oceanEngine;
    let reference = new Float64Array(model.oceanEngine.layers * C);
    for (let n = 1; n <= 5; n++) {
      await model.step(DT);
      const saved = await model.oceanEngine.serialize();
      reference = fold(reference, velocities(model, saved), saved.h, cellOcean, C, n > 1 ? DT : 0, n - 1);
      assert.equal(saved.upwellingFolds, n, `folds after step ${n}`);
      assert.equal(saved.upwellingTime, n * DT, `time of the last fold after step ${n}`);
      const { max, scale } = worstDifference(reference, saved.upwellingMemory);
      if (n === 1 || n === 5) console.log(`ocean step ${n}: largest |w| ${scale.toExponential(2)} m/s, GPU parts from the CPU fold by ${max.toExponential(2)}`);
      assert.ok(scale > 1e-7 && max <= 1e-4 * scale, `memory after step ${n}: ${max} against ${scale}`);
    }
  } finally { model.destroy(); }
});

test('an uploaded state restores the memory with its folds and time and the next fold continues the mean; one without it starts afresh', { skip }, async () => {
  const model = await build();
  try {
    const C = model.mesh.nCells, gpu = model.oceanEngine, surfaceT = Float64Array.from(model.state[3]), ice = Float64Array.from(model.state[6]);
    for (let n = 0; n < 3; n++) await model.step(DT);
    const saved = await gpu.serialize();
    assert.equal(saved.upwellingFolds, 3);

    const bare = { ...saved };
    delete bare.upwellingMemory;
    gpu.upload(bare, Float64Array.from(surfaceT), ice);
    assert.ok(!('upwellingMemory' in await gpu.serialize()), 'a state without a memory uploads with none');
    assert.ok((await gpu.readSection('UPWM')).every((v) => v === 0));

    gpu.upload(saved, Float64Array.from(surfaceT), ice);
    const again = await gpu.serialize();
    assert.equal(again.upwellingFolds, 3);
    assert.equal(again.upwellingTime, 3 * DT);
    same(again.upwellingMemory, saved.upwellingMemory, 'the memory goes back as it was saved');

    await model.step(DT);
    const next = await gpu.serialize(), expected = fold(saved.upwellingMemory, velocities(model, next), next.h, gpu.cellOcean, C, DT, 3);
    const { max, scale } = worstDifference(expected, next.upwellingMemory);
    assert.equal(next.upwellingFolds, 4);
    assert.equal(next.upwellingTime, 4 * DT);
    assert.ok(scale > 1e-7 && max <= 1e-4 * scale, `after a step the uploaded memory continues the mean: ${max} against ${scale}`);
  } finally { model.destroy(); }
});

test('a CPU ocean loaded from the device\'s state reads the same field at any depth, and the prognostic fields are the same with the memory and without it', { skip }, async () => {
  const runs = [];
  for (const withMemory of [true, false]) {
    const model = await build();
    try {
      const gpu = model.oceanEngine, surfaceT = Float64Array.from(model.state[3]), ice = Float64Array.from(model.state[6]);
      for (let n = 0; n < 4; n++) await model.step(DT);
      const saved = await gpu.serialize();
      if (withMemory) {
        const cpu = createOcean(model.mesh, { geography: model.geography, ...HYBRID });
        cpu.initialize(Float64Array.from(surfaceT), ice);
        cpu.load(saved, Float64Array.from(surfaceT), ice);
        for (const depth of [20, 137, 903, 2500]) {
          const expected = cpu.fields(depth).upwelling, got = (await gpu.frame({ fields: ['upwelling'], depth, upwellingDepth: depth })).fields.upwelling;
          const { max, scale } = worstDifference(expected, got);
          assert.ok(scale > 1e-7 && max <= 1e-5 * scale, `at ${depth} m the CPU reads the device's memory as the device does: ${max} against ${scale}`);
        }
      } else delete saved.upwellingMemory;
      gpu.upload(saved, Float64Array.from(surfaceT), ice);
      for (let n = 0; n < 3; n++) await model.step(DT);
      const after = await gpu.serialize();
      runs.push(['h', 'u', 'T', 'S', 'eta'].map((name) => [name, Float64Array.from(after[name])]));
    } finally { model.destroy(); }
  }
  runs[0].forEach(([name, values], n) => same(values, runs[1][n][1], `${name} after a load with the memory and without it`));
});
