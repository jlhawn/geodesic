/*
 * The model's longwave g-points on a column of the given air: the gases'
 * paths exactly as radiation.column builds them, each slice's share of the
 * ground's σTs⁴, its optical depth from the top split by absorber, the
 * pressure where that depth reaches 1 (null where the whole column stays
 * below 1), and what it sends to space. The air is dry unless `humidity`
 * scales the profile's vapor (or gives q per layer); vapor absorbs through
 * its lines and its self-continuum, each of which can be switched off by
 * zeroing its coefficients. The slices' fluxes sum to clearLongwave's on
 * the same table.
 */
import { gasPaths, opticalDepth, planckShare, clearLongwave, longwaveTableFor, GAS_MOLAR, LONGWAVE_CONSTANTS } from '../js/physics/longwave.module.js';
import { GREENHOUSE_GASES, STEFAN_BOLTZMANN } from '../js/physics/radiation.module.js';
import { OZONE_CM_ATM } from '../js/physics/shortwaveGases.module.js';

export const GASES = ['carbonDioxide', 'ozone', 'methane', 'nitrousOxide'];
export const VAPOR = ['vaporLines', 'continuum'];
export const ABSORBERS = [...GASES, ...VAPOR];
export const CLEAR = 0.05;
const PATHS = { carbonDioxide: 'co2', ozone: 'o3', methane: 'ch4', nitrousOxide: 'n2o', vaporLines: 'line', continuum: 'continuum' };
const ROW = { carbonDioxide: 2, ozone: 3, methane: 4, nitrousOxide: 5, vaporLines: 0, continuum: 1 };
const variants = new WeakMap();

export function airColumn(profile, { carbonDioxide = 1, gases = {}, humidity = 0 } = {}) {
  const on = { ...Object.fromEntries(GASES.map((g) => [g, true])), ...gases };
  const { levels, surfacePressure: ps, temperature: T, groundTemperature: Ts } = profile, K = T.length;
  const mass = Float64Array.from({ length: K }, (_, k) => (levels[k + 1] - levels[k]) * ps / LONGWAVE_CONSTANTS.gravity);
  const q = typeof humidity === 'number' ? Float64Array.from({ length: K }, (_, k) => (humidity > 0 ? humidity * profile.humidity[k] : 0)) : Float64Array.from(humidity);
  const mixed = (amount, molar) => Float64Array.from(q, (x) => amount * molar / GAS_MOLAR.air * Math.max(0, 1 - Math.max(0, x)));
  return {
    levels, ps, T, Ts, q,
    o3: Float64Array.from(profile.ozone, (x, k) => (on.ozone ? x * OZONE_CM_ATM / mass[k] : 0)),
    co2: mixed(on.carbonDioxide ? carbonDioxide * GREENHOUSE_GASES.carbonDioxide : 0, GAS_MOLAR.co2),
    ch4: mixed(on.methane ? GREENHOUSE_GASES.methane : 0, GAS_MOLAR.ch4),
    n2o: mixed(on.nitrousOxide ? GREENHOUSE_GASES.nitrousOxide : 0, GAS_MOLAR.n2o),
  };
}

export function dryColumn(profile, options = {}) {
  return airColumn(profile, { ...options, humidity: 0 });
}

export function vaporTable(levels, { vaporLines = true, continuum = true } = {}) {
  const table = longwaveTableFor(levels), key = `${vaporLines} ${continuum}`;
  if (vaporLines && continuum) return table;
  if (!variants.has(table)) variants.set(table, new Map());
  const made = variants.get(table);
  if (!made.has(key)) made.set(key, { ...table, points: table.points.map((row) => Object.assign(row.slice(), vaporLines ? {} : { [ROW.vaporLines]: 0 }, continuum ? {} : { [ROW.continuum]: 0 })) });
  return made.get(key);
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
  const column = airColumn(profile, options), table = vaporTable(profile.levels, options.gases), paths = gasPaths(column, table);
  const { levels, ps, T, Ts } = column, K = T.length, emission = STEFAN_BOLTZMANN * Ts ** 4;
  const tau = new Float64Array(K), planck = new Float64Array(K);
  const points = table.points.map((row, j) => {
    const gas = Object.fromEntries(ABSORBERS.map((g) => [g, LONGWAVE_CONSTANTS.diffusivity * row[ROW[g]] * paths[PATHS[g]].reduce((a, b) => a + b, 0)]));
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
    const strongest = ABSORBERS.reduce((best, g) => (gas[g] > gas[best] ? g : best), ABSORBERS[0]);
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
  const precipitableWater = column.q.reduce((sum, q, k) => sum + Math.max(0, q) * (levels[k + 1] - levels[k]) * ps / LONGWAVE_CONSTANTS.gravity, 0);
  return { points, order, groundEmission: emission, olr, dlr, greenhouse: emission - olr, escaping: points.reduce((a, p) => a + p.escaping, 0), precipitableWater, groundTemperature: Ts, surfacePressure: ps };
}
