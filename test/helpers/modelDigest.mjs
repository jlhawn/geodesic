import { createHash } from 'node:crypto';
import { Grid } from '../../js/grid.module.js';
import { createModel } from '../../js/model.module.js';
import { initializeState } from '../../js/physics/init.module.js';
import { PREVIOUS_LEVELS, PREVIOUS_RADIATION, PREVIOUS_MOIST, PREVIOUS_OCEAN } from './previousDefaults.mjs';

export const OVERCAST = { cloudCover: 'overcast' };
export const UNSCATTERED = { rayleighDepth: 0, nearInfraredRayleigh: 0, landAerosol: 0, seaAerosol: 0, skylight: 0.15, upwardAbsorption: false };
export const GRAY_GASES = { longwaveScheme: 'gray', solarGases: 'lacisHansen' };
export const REDIAGNOSED = { prognosticHeight: false, gateMemory: 0 };
export const GREY_ICE = { iceAlbedo: 0.5, meltingIceAlbedo: 0.5, snowAgeing: false };
export const RAYLEIGH_TOP = { surface: { topDragDays: 5, spongeDays: 0 }, gravityWaves: false };
export function modelDigest(radiation, moist = {}, ice = {}, { surface = {}, gravityWaves = {} } = {}) {
  const model = createModel(new Grid(16), { levels: PREVIOUS_LEVELS, ocean: { ...PREVIOUS_OCEAN, eddyDiffusivity: 0, closureFill: 0 }, divergenceDamping: 0, moist: { ...PREVIOUS_MOIST, cloudLifetime: 3 * 3600, capeClosure: 'threshold', plumeSourceDepth: 'boundaryLayer', cumulusClosure: 0.06, convectionType: 'top', plumeEntrainmentLaw: 'gregory', plumePhase: 'liquid', plumeConversion: 'zhangMcFarlane', plumeCape: 70, stratiformLifetime: null, condensation: 'saturation', iceSaturation: false, iceFall: null, fogDroplets: null, fogDeposition: 0, ...moist }, boundaryLayer: { turbulence: 'dry', entrainment: { efficiency: 0, shear: 0 }, dragCoefficient: 1.5e-3 }, surface: { dragCoefficient: 1.5e-3, ...surface }, radiation: { ...PREVIOUS_RADIATION, exchangeCoefficient: 1.5e-3, longwaveLines: 'scaled', ...radiation }, ice, gravityWaves });
  initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
  for (let n = 0; n < 12; n++) model.step(900);
  const hash = createHash('sha256');
  for (const a of [...model.state, model.radiation.stratus, model.radiation.stratusFraction, model.radiation.surfaceFlux, model.radiation.outgoing]) hash.update(new Uint8Array(a.buffer, a.byteOffset, a.byteLength));
  return { digest: hash.digest('hex').slice(0, 32), radiation: model.radiation };
}
