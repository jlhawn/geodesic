import { createSeaIce, FREEZING_POINT } from '../../js/physics/ice.module.js';
import { hybridOcean, DT } from './layered.mjs';

export const POLAR = 62 * Math.PI / 180;

export function rig(iceExchange, { thickness = 1.2, cover = 0.8, flux = 0, global = false, geography = null, iceDrift = null, water = 290, keepSmallIce = false } = {}) {
  const r = hybridOcean({ iceExchange, everySteps: 1, ...(geography ? { geography } : {}), ...(iceDrift ? { iceDrift } : {}) }, () => 4000, water), { mesh, C, ocean, surfaceT, ice } = r;
  r.seaIce = createSeaIce(mesh, { heatCapacity: ocean.capacity, oceanCells: ocean.cellOcean, iceExchange: ocean.iceExchange, ...(keepSmallIce ? { keepSmallIce } : {}) });
  ocean.useSeaIce(r.seaIce.concentration, r.seaIce);
  for (let i = 0; i < C; i++) if (global || Math.abs(mesh.latCell[i]) > POLAR) { ice[i] = thickness; surfaceT[i] = FREEZING_POINT - 3; r.seaIce.concentration[i] = cover; }
  ocean.initialize(surfaceT, ice, { concentration: r.seaIce.concentration });
  for (let i = 0; i < C; i++) if (!(ice[i] > 0)) surfaceT[i] = ocean.T0[i];
  r.heat = new Float64Array(C).fill(flux);
  r.volume = () => { let v = 0; for (let i = 0; i < C; i++) v += mesh.areaCell[i] * ice[i] * (ice[i] > 0 ? r.seaIce.concentration[i] : 0); return v; };
  r.cellEnergy = () => Float64Array.from({ length: C }, (_, i) => {
    let q = 0;
    for (let k = 1; k < ocean.layers; k++) q += ocean.Q[k * C + i] - 272 * ocean.h[k * C + i];
    const iced = ice[i] > 0, top = (iced ? ocean.h[i] * FREEZING_POINT + (r.seaIce.turbulent ? r.seaIce.waterHeat[i] : 0) / ocean.rhoCp : ocean.h[i] * surfaceT[i]) - 272 * ocean.h[i];
    return ocean.rhoCp * (q + top) + (iced ? r.seaIce.energy(surfaceT[i], ice[i], r.seaIce.snow[i], r.seaIce.cover(i, ice[i])) : 0);
  });
  r.cycle = () => {
    ocean.readSurface(surfaceT, ice, r.seaIce.concentration);
    ocean.salt();
    ocean.writeSurface(surfaceT, r.seaIce.oceanFlux, DT);
    for (let i = 0; i < C; i++) r.seaIce.update(surfaceT, ice, r.heat, i, DT);
  };
  r.seed = (rise = 1, cold = 0.2) => {
    let n = 0;
    for (let i = 0; i < C; i++) if (ice[i] > 0) r.seaIce.waterHeat[i] = ocean.capacity[i] * (n++ % 7 === 3 ? -cold : rise);
  };
  r.step = () => {
    ocean.advance(surfaceT, ice, r.seaIce.oceanFlux, r.calm, DT, r.seaIce.concentration);
    for (let i = 0; i < C; i++) r.seaIce.update(surfaceT, ice, r.heat, i, DT);
  };
  return r;
}
