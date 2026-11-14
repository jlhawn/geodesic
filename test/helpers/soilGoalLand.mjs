import { initializeState } from '../../js/physics/init.module.js';
import { startPlaceholders } from '../../js/physics/land.module.js';

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
