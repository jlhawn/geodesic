import { OZONE_GRID, OZONE_PROFILES } from './ozoneTable.module.js';

/*
 * The ozone climatology: the AFGL atmospheres' ozone (Anderson et al. 1986;
 * js/physics/ozoneTable.module.js, scripts/ozoneTable.mjs) as the column
 * above each pressure, blended by latitude and season. The tropical profile
 * holds equatorward of 15°, the midlatitude ones at 45° and the subarctic
 * ones poleward of 60°, linear in latitude between; each hemisphere's summer
 * profile holds on 15 July and its winter one on 15 January (the AFGL
 * months), blended by the cosine of the time of year between. A layer's
 * ozone is the difference of the blended column above its two interfaces,
 * linear in pressure between the grid's pressures (log-spaced, 0.1 Pa to
 * 1100 hPa) and to zero above the first.
 */
export const OZONE_ROWS = ['tropical', 'midlatitudeSummer', 'midlatitudeWinter', 'subarcticSummer', 'subarcticWinter'];
export const SUMMER_DAY = 117;
const STEP = Math.log(OZONE_GRID.highest / OZONE_GRID.lowest) / (OZONE_GRID.points - 1);
const PRESSURES = Float64Array.from({ length: OZONE_GRID.points }, (_, j) => OZONE_GRID.lowest * Math.exp(STEP * j));

// The five profiles' weights at latitude lat (radians) and year fraction f (0 at the March equinox).
export function ozoneWeights(lat, f, out = new Float64Array(5)) {
  const a = Math.abs(lat) * 180 / Math.PI, season = Math.cos(2 * Math.PI * (f - SUMMER_DAY / 365)) * (lat < 0 ? -1 : 1), summer = 0.5 * (1 + season);
  const middle = Math.min(1, Math.max(0, (a - 15) / 30)) * (1 - Math.min(1, Math.max(0, (a - 45) / 15))), polar = Math.min(1, Math.max(0, (a - 45) / 15));
  out[0] = 1 - middle - polar;
  out[1] = middle * summer; out[2] = middle * (1 - summer);
  out[3] = polar * summer; out[4] = polar * (1 - summer);
  return out;
}

// The blended column (cm-atm) above pressure p (Pa).
export function ozoneAbove(p, weights) {
  if (!(p > 0)) return 0;
  const x = Math.log(p / OZONE_GRID.lowest) / STEP;
  const j = Math.min(OZONE_GRID.points - 2, Math.floor(x));
  let column = 0;
  if (j < 0) { for (let r = 0; r < 5; r++) column += weights[r] * OZONE_PROFILES[r][0]; return column * p / PRESSURES[0]; }
  const t = (p - PRESSURES[j]) / (PRESSURES[j + 1] - PRESSURES[j]);
  for (let r = 0; r < 5; r++) { const row = OZONE_PROFILES[r]; column += weights[r] * (row[j] + (row[j + 1] - row[j]) * t); }
  return column;
}
