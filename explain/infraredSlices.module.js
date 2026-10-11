/*
 * The model's longwave g-points on a dry column (q = 0): the gases' paths
 * exactly as radiation.column builds them, each slice's share of the ground's
 * σTs⁴, its optical depth from the top, the pressure where that depth
 * reaches 1 (null where the whole column stays below 1), and what it sends
 * to space. The slices' fluxes sum to clearLongwave's.
 */
import { gasPaths, opticalDepth, planckShare, clearLongwave, longwaveTableFor, GAS_MOLAR, LONGWAVE_CONSTANTS } from '../js/physics/longwave.module.js';
import { GREENHOUSE_GASES, STEFAN_BOLTZMANN } from '../js/physics/radiation.module.js';
import { OZONE_CM_ATM } from '../js/physics/shortwaveGases.module.js';

export const GASES = ['carbonDioxide', 'ozone', 'methane', 'nitrousOxide'];
export const CLEAR = 0.05;
const PATHS = { carbonDioxide: 'co2', ozone: 'o3', methane: 'ch4', nitrousOxide: 'n2o' };
const ROW = { carbonDioxide: 2, ozone: 3, methane: 4, nitrousOxide: 5 };

export function dryColumn(profile, { carbonDioxide = 1, gases = {} } = {}) {
  const on = { ...Object.fromEntries(GASES.map((g) => [g, true])), ...gases };
  const { levels, surfacePressure: ps, temperature: T, groundTemperature: Ts } = profile, K = T.length;
  const mass = Float64Array.from({ length: K }, (_, k) => (levels[k + 1] - levels[k]) * ps / LONGWAVE_CONSTANTS.gravity);
  const mixed = (amount, molar) => new Float64Array(K).fill(amount * molar / GAS_MOLAR.air);
  return {
    levels, ps, T, Ts, q: new Float64Array(K),
    o3: Float64Array.from(profile.ozone, (x, k) => (on.ozone ? x * OZONE_CM_ATM / mass[k] : 0)),
    co2: mixed(on.carbonDioxide ? carbonDioxide * GREENHOUSE_GASES.carbonDioxide : 0, GAS_MOLAR.co2),
    ch4: mixed(on.methane ? GREENHOUSE_GASES.methane : 0, GAS_MOLAR.ch4),
    n2o: mixed(on.nitrousOxide ? GREENHOUSE_GASES.nitrousOxide : 0, GAS_MOLAR.n2o),
  };
}

function atLogPressure(p, pressures, values, bottom = null) {
  const K = pressures.length;
  if (p <= pressures[0]) return values[0];
  if (p >= pressures[K - 1]) return bottom ? values[K - 1] + (bottom[1] - values[K - 1]) * Math.log(p / pressures[K - 1]) / Math.log(bottom[0] / pressures[K - 1]) : values[K - 1];
  let k = 1;
  while (pressures[k] < p) k++;
  const f = Math.log(p / pressures[k - 1]) / Math.log(pressures[k] / pressures[k - 1]);
  return values[k - 1] + f * (values[k] - values[k - 1]);
}

export function slices(profile, options = {}) {
  const column = dryColumn(profile, options), table = longwaveTableFor(profile.levels), paths = gasPaths(column, table);
  const { levels, ps, T, Ts } = column, K = T.length, emission = STEFAN_BOLTZMANN * Ts ** 4;
  const tau = new Float64Array(K), planck = new Float64Array(K);
  const points = table.points.map((row, j) => {
    const gas = Object.fromEntries(GASES.map((g) => [g, LONGWAVE_CONSTANTS.diffusivity * row[ROW[g]] * paths[PATHS[g]].reduce((a, b) => a + b, 0)]));
    let depth = 0, level = null, layer = -1;
    for (let k = 0; k < K; k++) {
      tau[k] = opticalDepth(row, paths, k);
      planck[k] = planckShare(row, T[k]) * STEFAN_BOLTZMANN * T[k] ** 4;
      if (level === null && depth + tau[k] >= 1) {
        layer = k;
        level = ps * (levels[k] + (1 - depth) / tau[k] * (levels[k + 1] - levels[k]));
      }
      depth += tau[k];
    }
    let down = 0;
    for (let k = 0; k < K; k++) down = down * Math.exp(-tau[k]) - planck[k] * Math.expm1(-tau[k]);
    const share = planckShare(row, Ts), ground = share * emission;
    let up = ground;
    for (let k = K - 1; k >= 0; k--) up = up * Math.exp(-tau[k]) - planck[k] * Math.expm1(-tau[k]);
    const strongest = GASES.reduce((best, g) => (gas[g] > gas[best] ? g : best), GASES[0]);
    return {
      index: j, share, ground, up, down, depth, gas, layer, level,
      dominant: depth >= CLEAR ? strongest : null,
      escaping: share * Math.exp(-depth),
      temperature: level === null ? Ts : atLogPressure(level, profile.pressure, T),
      height: level === null ? 0 : atLogPressure(level, profile.pressure, profile.height, [ps, 0]),
    };
  });
  const totals = clearLongwave(column, { table });
  const olr = totals.up[0], dlr = totals.down[K];
  const order = points.map((p) => p.index).sort((a, b) => {
    const p = points[a], q = points[b];
    if ((p.level === null) !== (q.level === null)) return p.level === null ? -1 : 1;
    return p.level === null ? p.depth - q.depth : q.level - p.level;
  });
  return { points, order, groundEmission: emission, olr, dlr, greenhouse: emission - olr, escaping: points.reduce((a, p) => a + p.escaping, 0), groundTemperature: Ts, surfacePressure: ps };
}
