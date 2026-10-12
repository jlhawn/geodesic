/*
 * CONDUCTIVITY, SNOW_CONDUCTIVITY, SNOW_DENSITY and MINIMUM_THICKNESS copy
 * createSeaIce's defaults; test/seaIceYears.test.mjs checks them against its update.
 */
import { FREEZING_POINT, MELTING_POINT } from '../js/physics/ice.module.js';
import { salinityForDensity } from '../js/ocean/seawater.module.js';
import { seaIceColumn, RHO_CP, DEEP_CLASS, DEPTHS } from './seaIceColumn.module.js';

export const START_DAY = 183, START_T = MELTING_POINT, YEARS = 20, SNOWFALL = 100;
export const CONDUCTIVITY = 2.0, SNOW_CONDUCTIVITY = 0.31, SNOW_DENSITY = 300, MINIMUM_THICKNESS = 0.1;
export const SINKING_SALINITY = salinityForDensity(DEEP_CLASS, FREEZING_POINT);
export { DEPTHS, DEEP_CLASS, FREEZING_POINT };

export function conduction(thickness, snow, skin) {
  if (!(thickness > 0)) return 0;
  return (FREEZING_POINT - skin) / (Math.max(thickness, MINIMUM_THICKNESS) / CONDUCTIVITY + snow / (SNOW_DENSITY * SNOW_CONDUCTIVITY));
}

export const snowDepth = (snow) => snow / SNOW_DENSITY;
export const iceEquivalent = (snow) => snowDepth(snow) * CONDUCTIVITY / SNOW_CONDUCTIVITY;
export const heatAboveFreezing = (temperature, depth) => RHO_CP * depth * Math.max(0, temperature - FREEZING_POINT);

export function iceSeason(present, from, to) {
  let best = null, start = -1;
  for (let d = from; d <= to; d++) {
    const iced = d < to && present[d] > 0;
    if (iced && start < 0) start = d;
    if (!iced && start >= 0) { if (!best || d - start > best[1] - best[0]) best = [start, d]; start = -1; }
  }
  if (!best) return { forms: null, melts: null, days: 0, allYear: false, none: true };
  const [first, end] = best;
  const carried = first === from && from > 0 && present[from - 1] > 0, lasting = end === to && (to >= present.length || present[to] > 0);
  return { forms: carried ? null : first, melts: lasting ? null : end, days: end - first, allYear: end - first === to - from, none: false };
}

export function seaIceYears({ latitude = 75, depth = DEPTHS.floor, ice = 'bright', convergence = 85, snow = false, years = YEARS } = {}) {
  const run = seaIceColumn({ latitude, depth, ice, convergence, snowfall: snow ? SNOWFALL : 0, years, startDay: START_DAY, start: START_T });
  const r = run.record, skin = new Float32Array(run.days), conducted = new Float32Array(run.days);
  for (let d = 0; d < run.days; d++) {
    const A = r.cover[d];
    skin[d] = r.thickness[d] > 0 && A > 0 && A < 1 ? (r.surface[d] - (1 - A) * FREEZING_POINT) / A : r.surface[d];
    conducted[d] = conduction(r.thickness[d], r.snow[d], skin[d]);
  }
  function summary(y = years - 1) {
    const from = y * 365, to = from + 365, season = iceSeason(r.coverMax, from, to), base = run.year(y);
    const waterHigh = y === 0 ? Math.max(base.waterHigh, START_T) : base.waterHigh;
    return { ...base, ...season, waterHigh, brine: base.salinityHigh - base.salinityLow, summerHeat: heatAboveFreezing(waterHigh, depth) };
  }
  return { ...run, convergence, ice, snow, skin, conducted, summary };
}
