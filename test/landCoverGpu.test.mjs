import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { readTopography } from '../js/geography.module.js';
import { startPlaceholders } from '../js/physics/land.module.js';
import { PREVIOUS_LEVELS, PREVIOUS_PCAPE_MINIMUM, PREVIOUS_BOUNDARY_LAYER } from './helpers/previousDefaults.mjs';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const topography = readTopography('etopo1');
const FAST = { growthTime: 3 * 3600, declineTime: 2 * 3600, snowDeclineTime: 4 * 3600, treeGrowthTime: 2 * 3600, treeDeclineTime: 3600, treeMoisture: false };
const neutral = startPlaceholders().neutral.carbon;

function prepare(model) {
  const C = model.mesh.nCells, land = model.geography.land, init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  for (let i = 0; i < C; i++) if (land[i]) model.state[6][i] = 0;
  if (model.load) model.load();
  model.ocean.initialize(model.state[3], model.state[6]);
  model.land.initialize();
  const each = (value) => Float64Array.from({ length: C }, (_, i) => (land[i] ? value(i) : 0));
  model.land.load({
    soil: each(() => 40), snow: each((i) => (model.mesh.latCell[i] > 1.0 ? 5 : 0)), vegetation: each(() => 0.5), canopy: each(() => 0.2), soilCarbon: each(() => neutral),
    seasonLength: each(() => 0.8), seasonWarmth: each(() => 4), rainMean: each(() => 3), demandMean: each(() => 4), litterMean: each(() => 0.3), decayMean: each(() => 0.2),
    record: Float64Array.of(0, 4, 4),
  });
  return model;
}

const apart = (a, b, keep) => { let worst = 0; for (let i = 0; i < a.length; i++) if (keep(i)) worst = Math.max(worst, Math.abs(a[i] - b[i])); return worst; };

for (const [coverGoal, snowHold] of [['water', false], ['observed', false], ['observed', true]]) {
  test(`from an observed start both engines relax the cover and the trees alike under coverGoal '${coverGoal}'${snowHold ? ' with snowHold' : ''}, hold the carbon at the placeholder, and a jump leaves the trees${coverGoal === 'observed' ? ', with the boundary layer diagnosed after the deposit' : ''}`, { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
    const land = { ...FAST, coverGoal, snowHold };
    const cpu = prepare(createModel(new Grid(16), { topography, land, levels: PREVIOUS_LEVELS, moist: PREVIOUS_PCAPE_MINIMUM, boundaryLayer: { ...(coverGoal === 'observed' ? PREVIOUS_BOUNDARY_LAYER : {}), velocityScale: 'capped' } })), gpu = prepare(await createGpuModel(new Grid(16), { topography, land, levels: PREVIOUS_LEVELS, moist: PREVIOUS_PCAPE_MINIMUM, boundaryLayer: { ...(coverGoal === 'observed' ? PREVIOUS_BOUNDARY_LAYER : {}), velocityScale: 'capped' } }));
    const C = cpu.mesh.nCells, grown = (i) => cpu.geography.land[i] && !cpu.geography.iceSheet[i], count = Array.from({ length: C }, (_, i) => i).filter(grown).length;
    assert.deepEqual([...cpu.land.record], [0, 4, 4]);
    for (let n = 0; n < 16; n++) { cpu.step(900); await gpu.step(900); }
    await gpu.sync();
    const saved = await gpu.land.serialize(), vegetation = apart(cpu.land.vegetation, saved.vegetation, grown), trees = apart(cpu.land.canopy, saved.canopy, grown), free = apart(cpu.land.snowFreeCover, saved.snowFreeCover, grown);
    let moved = 0, grew = 0;
    for (let i = 0; i < C; i++) if (grown(i)) { if (Math.abs(cpu.land.vegetation[i] - 0.5) > 0.02) moved++; if (Math.abs(cpu.land.canopy[i] - 0.2) > 0.02) grew++; }
    const carbon = apart(cpu.land.soilCarbon, saved.soilCarbon, grown), held = apart(saved.soilCarbon, new Float64Array(C).fill(neutral), grown);
    console.log(`coverGoal '${coverGoal}', ${count} vegetated cells at N=16 after 16 steps: cover differs by ${vegetation.toExponential(1)}, trees ${trees.toExponential(1)}, snow-free cover ${free.toExponential(1)}, carbon ${carbon.toExponential(1)} (held to ${held.toExponential(1)} of the placeholder); the cover moved in ${moved} cells and the trees in ${grew}`);
    assert.ok(moved > 0.2 * count && grew > 0.2 * count, `${moved} cells with moving cover and ${grew} with moving trees of ${count}`);
    assert.ok(vegetation < 1e-4 && trees < 1e-4 && free < 1e-4, `cover ${vegetation}, trees ${trees}, snow-free cover ${free}`);
    assert.ok(carbon < 1e-6 && held < 1e-6, `carbon ${carbon}, held to ${held}`);
    const snowy = Array.from({ length: C }, (_, i) => i).filter((i) => grown(i) && cpu.land.snow[i] > 0 && saved.snow[i] > 0);
    const kept = snowy.filter((i) => Math.abs(cpu.land.vegetation[i] - 0.5) < 1e-6 && Math.abs(saved.vegetation[i] - 0.5) < 1e-6).length;
    if (snowHold) assert.ok(snowy.length > 0.05 * count && kept > 0.5 * snowy.length, `${kept} of ${snowy.length} snowy cells keep their cover`);
    else assert.ok(kept < 0.1 * snowy.length, `under snow the cover decays: ${kept} of ${snowy.length} snowy cells keep it`);
    const before = Float64Array.from(saved.canopy);
    await gpu.land.jump();
    cpu.land.jump(Math.fround);
    const after = await gpu.land.serialize();
    let changed = 0;
    for (let i = 0; i < C; i++) { assert.equal(after.canopy[i], before[i], `GPU trees of cell ${i}`); if (grown(i) && Math.abs(after.soilCarbon[i] - neutral) > 1e-3) changed++; }
    const target = apart(cpu.land.soilCarbon, after.soilCarbon, grown), scale = Math.max(...after.soilCarbon);
    console.log(`jump: the carbon of ${changed} of ${count} cells leaves the placeholder, engines ${target.toExponential(1)} apart (of up to ${scale.toFixed(1)} kg/m²)`);
    assert.ok(changed > 0.9 * count && target < 1e-3 * scale, `${changed} cells changed, engines ${target} apart`);
    gpu.destroy();
  });
}
