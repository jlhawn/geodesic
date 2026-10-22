import { LONGWAVE_SPECTRAL_MODEL, LONGWAVE_POINTS } from './longwaveTable.module.js';
import { LONGWAVE_SPECTRAL_MODEL as BL36_SPECTRAL_MODEL, LONGWAVE_POINTS as BL36_POINTS } from './longwaveTableBl36.module.js';
import { LONGWAVE_SPECTRAL_MODEL as BL48_SPECTRAL_MODEL, LONGWAVE_POINTS as BL48_POINTS } from './longwaveTableBl48.module.js';
import { LONGWAVE_SPECTRAL_MODEL as BL60_SPECTRAL_MODEL, LONGWAVE_POINTS as BL60_POINTS } from './longwaveTableBl60.module.js';
import { LONGWAVE_SPECTRAL_MODEL as BL36_VOIGT_SPECTRAL_MODEL, LONGWAVE_POINTS as BL36_VOIGT_POINTS } from './longwaveTableBl36Voigt.module.js';
import { LONGWAVE_SPECTRAL_MODEL as BL36_BY_GAS_SPECTRAL_MODEL, LONGWAVE_POINTS as BL36_BY_GAS_POINTS } from './longwaveTableBl36ByGas.module.js';
import { sigmaGridName } from '../dynamics/sigmaCore.module.js';

/*
 * The clear-sky longwave gas optics: g-points of a simple spectral model
 * (scripts/longwaveFit.mjs, which fits it to RRTMG over the standard
 * atmospheres and writes the table). A layer's optical depth in g-point g is
 * 1.66 (k_line,g a + k_self,g b + k_CO2,g c + k_O3,g d + k_CH4,g m + k_N2O,g n)
 * with the layer's paths: a the vapour mass times p_V/p_ref (p_ref 500 hPa),
 * b the vapour mass times its vapour pressure (Pa) times
 * exp(tSelf (1/T - 1/296)), c the CO2 mass times p_V/p_ref exp(tCo2 (T - 250)),
 * d the ozone mass times (p_V/p_ref)^nO3, m and n the methane and nitrous
 * oxide masses times p/p_ref; p is the layer's mid-pressure and p_V, for the
 * vapour lines, CO2 and ozone, sqrt(p^2 + p_D^2): the Voigt half-width in its
 * root-sum-square form as a pressure, p_D (the table's dopplerH2o, dopplerCo2
 * and dopplerO3) the pressure at which the gas's Lorentz half-width equals
 * its Doppler half-width at 250 K. The g-point emits the share of sigma T^4
 * its row's quartic in (T - 250)/100 gives.
 *
 * A table with Doppler cores (dopplerCoreH2o, dopplerCoreCo2, dopplerCoreO3)
 * has twelve more coefficients at the end of each row (LONGWAVE_CORES): for
 * the vapor, CO2 and ozone, one per node of a line's core, x = |nu -
 * nu_line|/alpha_D between consecutive CORE_EDGES, alpha_D the Doppler 1/e
 * half-width. A node's coefficient is its mean absorption in the Doppler
 * limit, and its path the gas's own (the vapor mass, the CO2 mass times
 * exp(tCo2 (T - 250)), the ozone mass) times coreFactor: the node's mean of
 * the Voigt function H(a, x) over its mean of exp(-x^2), a = sqrt(ln 2) p/p_D
 * the Lorentz half-width over alpha_D. Such a table's dopplerH2o, dopplerCo2
 * and dopplerO3 are 0: the far wings scale with the Lorentz width alone, as a
 * Voigt line's do.
 */
export const LONGWAVE_CONSTANTS = { pRef: 50000, diffusivity: 1.66, gravity: 9.806 };

/*
 * The rows with their quartics' coefficients shifted, in proportion to each
 * row's share at 250 K, so that the shares sum to 1 at every temperature
 * and the surface emits σT⁴ in all: the table spans 10-3250 cm-1 only.
 */
export function normalizedPoints(points) {
  const rows = points.map((row) => row.slice());
  for (let n = 0; n < 5; n++) {
    const excess = rows.reduce((s, row) => s + row[6 + n], 0) - (n === 0 ? 1 : 0), weight = rows.reduce((s, row) => s + row[6], 0);
    for (const row of rows) row[6 + n] -= excess * row[6] / weight;
  }
  return rows;
}
export const LONGWAVE_TABLE = { ...LONGWAVE_SPECTRAL_MODEL, points: normalizedPoints(LONGWAVE_POINTS) };
const BL36_TABLE = { ...BL36_SPECTRAL_MODEL, points: normalizedPoints(BL36_POINTS) };
const BL48_TABLE = { ...BL48_SPECTRAL_MODEL, points: normalizedPoints(BL48_POINTS) };
const BL60_TABLE = { ...BL60_SPECTRAL_MODEL, points: normalizedPoints(BL60_POINTS) };
const BL36_VOIGT_TABLE = { ...BL36_VOIGT_SPECTRAL_MODEL, points: normalizedPoints(BL36_VOIGT_POINTS) };
const BL36_BY_GAS_TABLE = { ...BL36_BY_GAS_SPECTRAL_MODEL, points: normalizedPoints(BL36_BY_GAS_POINTS) };
export const LONGWAVE_TABLE_FILES = { bl34: 'longwaveTable.module.js', bl36: 'longwaveTableBl36.module.js', bl36Voigt: 'longwaveTableBl36Voigt.module.js', bl36ByGas: 'longwaveTableBl36ByGas.module.js', bl48: 'longwaveTableBl48.module.js', bl60: 'longwaveTableBl60.module.js' };

/*
 * The table fitted on a grid's own columns: bl36's for bl36 and for bl42,
 * which shares its top layers, bl48's for bl48, bl60's for bl60, the one
 * fitted on bl34 for every other grid. With lines 'voigt', the table with
 * Doppler cores fitted on bl36's columns; with 'scaledByGas', bl36's scaled
 * lines refitted against each gas's RRTMG bands as well. bl48 and bl60 keep
 * their own tables whatever the lines: on their layers the bl36 'voigt'
 * fit's heating errors at 3-30 hPa are 1.4 to four times larger.
 */
export function longwaveTableFor(levels, lines = 'scaled') {
  if (lines !== 'scaled' && lines !== 'scaledByGas' && lines !== 'voigt') throw new Error(`longwaveLines must be 'scaled', 'scaledByGas' or 'voigt', not ${lines}`);
  const name = levels ? sigmaGridName(levels) : null;
  if (name === 'bl48') return BL48_TABLE;
  if (name === 'bl60') return BL60_TABLE;
  if (lines === 'voigt') return BL36_VOIGT_TABLE;
  if (lines === 'scaledByGas') return BL36_BY_GAS_TABLE;
  return name === 'bl36' || name === 'bl42' ? BL36_TABLE : LONGWAVE_TABLE;
}
export const GAS_MOLAR = { air: 28.964, h2o: 18.015, co2: 44.01, o3: 47.997, ch4: 16.04, n2o: 44.013 };
const STEFAN = 5.670374419e-8;

export function vaporPressure(q, p) {
  return q * p / (0.622 + 0.378 * q);
}

export function planckShare(row, T) {
  const t = (T - 250) / 100;
  return row[6] + t * (row[7] + t * (row[8] + t * (row[9] + t * row[10])));
}

/*
 * The paths of one layer into out[0..5] (line, continuum, CO2, ozone,
 * methane, nitrous oxide): p its mid-pressure (Pa), mass its air mass
 * (kg/m2), q its specific humidity, ozone its ozone mass (kg/m2), and the
 * well-mixed gases' mass mixing ratios.
 */
export function voigtScale(p, doppler) {
  return (doppler > 0 ? Math.sqrt(p * p + doppler * doppler) : p) / LONGWAVE_CONSTANTS.pRef;
}

export const longwaveCores = (model) => model.dopplerCoreCo2 > 0;
export const LONGWAVE_ABSORBERS = ['line', 'continuum', 'co2', 'o3', 'ch4', 'n2o'];
export const CORE_EDGES = [0, 0.6, 1.2, 1.8, 2.6], CORE_GASES = ['lineCore', 'co2Core', 'o3Core'];
export const LONGWAVE_CORES = CORE_GASES.flatMap((n) => CORE_EDGES.slice(1).map((_, j) => `${n}${j + 1}`));

/*
 * The core nodes' factors at log2 a = -10, -10 + 1/3, ... 10 (CORE_FACTOR_POINTS
 * each, node by node), with H(a, x) = (1/pi) ∫ exp(-(x + a tan t)^2) dt over
 * |t| < pi/2 and the node means taken at the midpoints of 32 equal parts.
 */
export const CORE_FACTOR_POINTS = 61;
let coreFactorTable = null;
export function coreFactors() {
  if (coreFactorTable) return coreFactorTable;
  const nodes = CORE_EDGES.length - 1, parts = 32, steps = 400, out = new Float64Array(nodes * CORE_FACTOR_POINTS);
  const voigt = (a, x) => {
    const lo = Math.atan((-x - 7) / a), h = (Math.atan((7 - x) / a) - lo) / steps;
    let s = 0;
    for (let i = 0; i < steps; i++) { const u = x + a * Math.tan(lo + (i + 0.5) * h); s += Math.exp(-u * u); }
    return s * h / Math.PI;
  };
  for (let j = 0; j < nodes; j++) {
    const xs = Array.from({ length: parts }, (_, i) => CORE_EDGES[j] + (i + 0.5) * (CORE_EDGES[j + 1] - CORE_EDGES[j]) / parts);
    const gauss = xs.reduce((s, x) => s + Math.exp(-x * x), 0);
    for (let t = 0; t < CORE_FACTOR_POINTS; t++) {
      const a = 2 ** (-10 + t / 3);
      out[j * CORE_FACTOR_POINTS + t] = xs.reduce((s, x) => s + voigt(a, x), 0) / gauss;
    }
  }
  return (coreFactorTable = out);
}

export function coreFactor(table, node, logA) {
  const u = Math.min(Math.max((logA + 10) * 3, 0), CORE_FACTOR_POINTS - 1), i = Math.min(Math.floor(u), CORE_FACTOR_POINTS - 2), w = u - i, r = node * CORE_FACTOR_POINTS + i;
  return table[r] * (1 - w) + table[r + 1] * w;
}

export function layerPaths(out, p, mass, T, q, ozone, co2, ch4, n2o, model = LONGWAVE_TABLE) {
  const scale = p / LONGWAVE_CONSTANTS.pRef, vapour = Math.max(0, q) * mass;
  out[0] = vapour * voigtScale(p, model.dopplerH2o ?? 0);
  out[1] = vapour * vaporPressure(Math.max(0, q), p) * Math.exp(model.tSelf * (1 / T - 1 / 296));
  out[2] = co2 * mass * voigtScale(p, model.dopplerCo2 ?? 0) * Math.exp(model.tCo2 * (T - 250));
  out[3] = ozone * Math.pow(voigtScale(p, model.dopplerO3 ?? 0), model.nO3);
  out[4] = ch4 * mass * scale;
  out[5] = n2o * mass * scale;
  if (model.dopplerCoreCo2 > 0) {
    const table = coreFactors(), nodes = CORE_EDGES.length - 1, carbon = co2 * mass * Math.exp(model.tCo2 * (T - 250));
    const logH2o = Math.log2(SQRT_LN2 * p / model.dopplerCoreH2o), logCo2 = Math.log2(SQRT_LN2 * p / model.dopplerCoreCo2), logO3 = Math.log2(SQRT_LN2 * p / model.dopplerCoreO3);
    for (let j = 0; j < nodes; j++) {
      out[6 + j] = vapour * coreFactor(table, j, logH2o);
      out[6 + nodes + j] = carbon * coreFactor(table, j, logCo2);
      out[6 + 2 * nodes + j] = ozone * coreFactor(table, j, logO3);
    }
  }
  return out;
}
const SQRT_LN2 = Math.sqrt(Math.LN2);

/*
 * The paths of a whole column given as layer means (column.levels sigma
 * interfaces, column.ps, and per-layer T, q and mass mixing ratios o3, co2,
 * ch4, n2o), for scripts.
 */
export function gasPaths(column, model = LONGWAVE_TABLE) {
  const K = column.T.length, names = longwaveCores(model) ? [...LONGWAVE_ABSORBERS, ...LONGWAVE_CORES] : LONGWAVE_ABSORBERS;
  const out = Object.fromEntries(names.map((n) => [n, new Float64Array(K)])), row = new Float64Array(names.length);
  for (let k = 0; k < K; k++) {
    const p = 0.5 * (column.levels[k] + column.levels[k + 1]) * column.ps, mass = (column.levels[k + 1] - column.levels[k]) * column.ps / LONGWAVE_CONSTANTS.gravity;
    layerPaths(row, p, mass, column.T[k], column.q[k], column.o3[k] * mass, column.co2[k], column.ch4[k], column.n2o[k], model);
    names.forEach((n, j) => { out[n][k] = row[j]; });
  }
  return out;
}

export function opticalDepth(row, paths, k) {
  let cores = 0;
  if (row.length > 11) for (let m = 0; m < LONGWAVE_CORES.length; m++) if (row[11 + m]) cores += row[11 + m] * paths[LONGWAVE_CORES[m]][k];
  return LONGWAVE_CONSTANTS.diffusivity * (row[0] * paths.line[k] + row[1] * paths.continuum[k] + row[2] * paths.co2[k] + row[3] * paths.o3[k] + row[4] * paths.ch4[k] + row[5] * paths.n2o[k] + cores);
}

/*
 * Clear-sky fluxes of a column over a black surface at column.Ts, at the
 * interfaces top first: up, down and net (up less down), W/m2.
 */
export function clearLongwave(column, { table = LONGWAVE_TABLE } = {}) {
  const K = column.T.length, paths = gasPaths(column, table);
  const up = new Float64Array(K + 1), down = new Float64Array(K + 1), tr = new Float64Array(K), B = new Float64Array(K);
  for (const row of table.points) {
    for (let k = 0; k < K; k++) {
      tr[k] = Math.exp(-opticalDepth(row, paths, k));
      B[k] = planckShare(row, column.T[k]) * STEFAN * column.T[k] ** 4;
    }
    let d = 0;
    for (let k = 0; k < K; k++) { d = d * tr[k] + B[k] * (1 - tr[k]); down[k + 1] += d; }
    let u = planckShare(row, column.Ts) * STEFAN * column.Ts ** 4;
    up[K] += u;
    for (let k = K - 1; k >= 0; k--) { u = u * tr[k] + B[k] * (1 - tr[k]); up[k] += u; }
  }
  return { up, down, net: Array.from(up, (x, k) => x - down[k]) };
}
