import { createHash } from 'node:crypto';
import { Grid } from '../../js/grid.module.js';
import { createModel } from '../../js/model.module.js';
import { initializeState } from '../../js/physics/init.module.js';
import { syntheticTopography } from '../../js/geography.module.js';
import { createLandSurface } from '../../js/physics/land.module.js';
import { PREVIOUS_LEVELS } from './previousDefaults.mjs';

export const topography = syntheticTopography(90, 180, (lat, lon) => ((Math.cos(lon) > 0 && Math.abs(lat) < 1.2) || lat < -1.15 ? 300 : -4000));
export const options = { topography, ocean: false, levels: PREVIOUS_LEVELS, land: { soil: 'slab' } };

export const digest = (arrays) => {
  const hash = createHash('sha256');
  for (const a of arrays) hash.update(Buffer.from(Float64Array.from(a).buffer));
  return hash.digest('hex');
};
export const landArrays = (land) => Object.keys(land).sort().map((key) => land[key]);

/*
 * A saved state at N = 16: land with snow on the northern part, sea ice
 * 1.5 m thick poleward of 1.2 rad carrying 12 kg/m² of snow in the land's
 * snow array, and the model's initial atmosphere. `iced` lists the sea-ice
 * cells.
 */
export function savedFixture() {
  const model = createModel(new Grid(16), options);
  const { mesh, geography } = model, C = mesh.nCells;
  const state = initializeState(model, {}).map((a) => Float64Array.from(a));
  const land = createLandSurface(mesh, geography, { soil: 'slab' });
  land.initialize();
  const iced = [];
  for (let i = 0; i < C; i++) {
    if (geography.land[i]) { state[6][i] = 0; land.soil[i] = 40; land.snow[i] = mesh.latCell[i] > 1.0 ? 5 : 0; }
    else if (mesh.latCell[i] > 1.2) { state[6][i] = 1.5; state[3][i] = 271.35; land.snow[i] = 12; iced.push(i); }
  }
  return { mesh, geography, state, iced, saved: { N: 16, time: 40 * 86400, land: land.serialize() } };
}

/*
 * The digest of a CPU model that has loaded the fixture's state, then its
 * land by `load(model, saved)`, and stepped once.
 */
export function cpuStepped(fixture, load) {
  const model = createModel(new Grid(16), options);
  fixture.state.forEach((a, k) => model.state[k].set(a));
  if (model.load) model.load();
  load(model, fixture.saved);
  model.step(900);
  return digest([...model.state, ...landArrays(model.land.serialize())]);
}
