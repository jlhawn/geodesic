import { test } from 'node:test';
import assert from 'node:assert/strict';
import { stateGeography, surfaceGeopotential, FRESH_GEOGRAPHY } from '../js/geography.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { decodeState, encodeState } from '../js/stateFile.module.js';
import { gridHeader } from '../js/stateGrid.module.js';
import { Grid } from '../js/grid.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const skip = !gpuAvailable && 'webgpu not installed';

test('a figure tool\'s GPU model stands on the terrain its state records, and on none for a state saved without terrain', { skip }, async () => {
  const { createGpuModel } = await import('../js/gpu/model.gpu.js');
  const { gpuModelFrom } = await import('../scripts/figures/figureState.mjs');
  const mesh = new Grid(16);
  for (const terrain of [{ conserveVolume: true }, false]) {
    const header = { N: 16, ...gridHeader(mesh), ...FRESH_GEOGRAPHY, terrain: terrain !== false, ...(terrain ? { terrainOptions: terrain } : {}) };
    const { topography, geography } = stateGeography(header);
    const source = await createGpuModel(mesh, { topography, geography, terrain });
    const { state } = source;
    initializeState(source, { geostrophic: !source.surfaceGeopotential }).forEach((values, a) => state[a].set(values));
    for (let i = 0; i < source.mesh.nCells; i++) if (source.geography.land[i]) state[6][i] = 0;
    source.load();
    source.ocean.initialize(state[3], state[6]);
    source.land.initialize();
    const ocean = await source.ocean.serialize(), land = await source.land.serialize();
    const [pi, theta, u, surfaceT, q, qc, ice] = state;
    const saved = await decodeState(encodeState({
      ...header, K: source.core.K, levels: source.core.levels, day: 0, time: 0, pi, theta, u, surfaceT, q, qc, ice,
      concentration: source.seaIce.concentration, mlmSubsidence: source.radiation.mlmSubsidence, ocean: { h: ocean.h, u: ocean.u, T: ocean.T, S: ocean.S, eta: ocean.eta }, land,
    }));
    const recorded = source.surfaceGeopotential && Array.from(source.surfaceGeopotential);
    if (terrain) assert.notDeepEqual(recorded, Array.from(surfaceGeopotential(mesh, source.geography, {})), 'conserveVolume moves the surface off the default');
    else assert.equal(recorded, null);
    source.destroy();
    const model = await gpuModelFrom(saved);
    try {
      assert.deepEqual(model.surfaceGeopotential && Array.from(model.surfaceGeopotential), recorded, JSON.stringify(terrain));
    } finally { model.destroy(); }
  }
});
