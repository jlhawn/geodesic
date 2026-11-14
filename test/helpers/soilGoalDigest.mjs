import { createHash } from 'node:crypto';
import { Grid } from '../../js/grid.module.js';
import { buildMesh } from '../../js/mesh.module.js';
import { createModel } from '../../js/model.module.js';
import { initializeState } from '../../js/physics/init.module.js';
import { readTopography } from '../../js/geography.module.js';
import { startPlaceholders } from '../../js/physics/land.module.js';
import { physicsSurfaceKernel } from '../../js/gpu/physics.gpu.js';
import { layoutFor } from '../../js/gpu/core.gpu.js';
import { PREVIOUS_LEVELS, PREVIOUS_PCAPE_MINIMUM, PREVIOUS_BOUNDARY_LAYER, PREVIOUS_AIR_SEA } from './previousDefaults.mjs';

const hash = (parts) => { const h = createHash('sha256'); for (const part of parts) h.update(part); return h.digest('hex').slice(0, 32); };

export function shaderDigest() {
  const texts = [];
  for (const held of [false, true]) for (const cover of [false, true]) for (const deep of [false, true]) texts.push(physicsSurfaceKernel(held, false, null, cover, deep, { treeline: true, vegetation: true }), physicsSurfaceKernel(held, true, null, cover, deep, {}));
  return hash(texts);
}

export function layoutDigest() {
  const mesh = buildMesh(new Grid(16)), args = [8, 0, 0, 0, null, false, 0, 0, false, false, false, false, false, 3, false, true];
  return hash([JSON.stringify(layoutFor(mesh, ...args)), JSON.stringify(layoutFor(mesh, ...args, true, true))]);
}

const FAST = { growthTime: 3 * 3600, declineTime: 2 * 3600, snowDeclineTime: 4 * 3600, treeGrowthTime: 2 * 3600, treeDeclineTime: 3600, treeMoisture: false };

export function prepareLand(model, soil = 40) {
  const neutral = startPlaceholders().neutral.carbon;
  const C = model.mesh.nCells, land = model.geography.land, init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  for (let i = 0; i < C; i++) if (land[i]) model.state[6][i] = 0;
  if (model.load) model.load();
  model.ocean.initialize(model.state[3], model.state[6]);
  model.land.initialize();
  const each = (value) => Float64Array.from({ length: C }, (_, i) => (land[i] ? value(i) : 0));
  model.land.load({
    soil: each((i) => (i % 3 ? soil : 280)), surface: each((i) => (i % 2 ? 15 : 0)), snow: each((i) => (model.mesh.latCell[i] > 1.0 ? 5 : 0)), vegetation: each(() => 0.5), canopy: each(() => 0.2), soilCarbon: each(() => neutral),
    seasonLength: each(() => 0.8), seasonWarmth: each(() => 4), rainMean: each(() => 3), demandMean: each(() => 4), litterMean: each(() => 0.3), decayMean: each(() => 0.2),
    record: Float64Array.of(0, 4, 4),
  });
  return model;
}

export async function landDigest(createGpuModel) {
  const topography = readTopography('etopo1'), land = { ...FAST, coverGoal: 'observed' };
  const gpu = prepareLand(await createGpuModel(new Grid(16), { topography, land, levels: PREVIOUS_LEVELS, surface: PREVIOUS_AIR_SEA, moist: PREVIOUS_PCAPE_MINIMUM, boundaryLayer: { ...PREVIOUS_BOUNDARY_LAYER, velocityScale: 'capped' } }));
  for (let n = 0; n < 4; n++) await gpu.step(900);
  await gpu.sync();
  const saved = await gpu.land.serialize();
  gpu.destroy();
  return hash(Object.keys(saved).sort().map((name) => new Uint8Array(saved[name].buffer, saved[name].byteOffset, saved[name].byteLength)));
}
