/*
 * The seastore figure's year: oceanColumn's column of the model's own ocean,
 * started as the model's analytic start leaves the sea at a latitude no
 * further from the equator than 45°, where the sea carries no ice: the
 * start's surface temperature from initializeState, its salinity from
 * createOcean's profile and its e-folding depth stretched by createOcean's
 * thermocline tilt, all at that latitude. The heating swings with the
 * annual harmonic of the sunlight a sea absorbs there under a sky that
 * passes TRANSMISSION of it, peaking at the June solstice, and the sea
 * loses DAMPING W/m² for each kelvin it is warmer than the start.
 */
import { createSeasonalColumn, RHO_CP, AIR_COLUMN, SLAB, STORE, DAY, YEAR_DAYS } from './oceanColumn.module.js';
import { dailySunlight, TOY } from './seaIceColumn.module.js';
import { FREEZING_POINT } from '../js/physics/ice.module.js';

export { RHO_CP, AIR_COLUMN, SLAB, STORE, YEAR_DAYS };
export const DAMPING = 10, STEP = 10800, TILT = 0.3, TRANSMISSION = TOY.transmission, POLEWARD = 45;
const RAD = Math.PI / 180;

export function startingSea(latitude) {
  const phi = Math.abs(latitude), polar = Math.max(0, Math.min(1, (phi - 50) / 20));
  return {
    temperature: Math.max(FREEZING_POINT + 0.3, 273.15 + 30 * Math.cos(latitude * RAD) ** 2 - 2 - 4 * polar),
    salinity: 34 + 2 * Math.exp(-(((phi - 25) / 20) ** 2)),
    stretch: 1 - TILT + 2 * TILT * Math.cos(latitude * RAD) ** 2,
  };
}

export function sunSwing(latitude, transmission = TRANSMISSION) {
  const { mean, water } = dailySunlight(latitude);
  let sum = 0;
  for (let d = 0; d < YEAR_DAYS; d++) sum += transmission * mean[d] * (1 - (Number.isNaN(water[d]) ? 0 : water[d])) * Math.sin(2 * Math.PI * (d + 0.5) / YEAR_DAYS);
  return 2 * sum / YEAR_DAYS;
}

export function columnOptions({ latitude = POLEWARD, heating = sunSwing(latitude), wind = 0.1, days = YEAR_DAYS } = {}) {
  if (!(Math.abs(latitude) <= POLEWARD)) throw new Error(`the column carries no ice, so it holds for latitudes within ${POLEWARD}° of the equator, not ${latitude}`);
  return { latitude, heating, wind, damping: DAMPING, dt: STEP, days, mean: startingSea(latitude).temperature, every: DAY / STEP, storeDepth: STORE, light: true };
}

export function createSeaStore(settings = {}) {
  const options = columnOptions(settings), column = createSeasonalColumn(options), L = column.layers, r = column.record;
  let sent = 0;
  return {
    column,
    layout: { heating: options.heating, wind: options.wind, mean: options.mean, layers: L, samples: Math.floor(column.steps / column.every) + 1, minimum: column.minimum },
    run(steps) { column.run(steps); },
    get done() { return column.done; },
    batch() {
      const from = sent, to = column.samples, slice = (a) => a.slice(from, to), wide = (a) => a.slice(from * L, to * L);
      sent = to;
      return {
        from, time: slice(r.time), forcing: slice(r.forcing), sst: slice(r.seas.model.sst), slab: slice(r.seas.slab.sst), depth: slice(r.mixedDepth), boundary: slice(r.boundaryDepth),
        salinity: slice(r.salinity), stored: slice(r.stored), heat: slice(r.seas.model.heat), slabHeat: slice(r.seas.slab.heat), bases: wide(r.bases), temperatures: wide(r.temperatures),
      };
    },
  };
}

export const airWarming = (joules) => joules / AIR_COLUMN;
