/*
 * Prescribed seasons over one column of the model's own ocean: createOcean
 * on the 42 cells of Grid(2), with all its classes, a flat 4000 m floor and
 * the same start in every column, a mixed layer FLOOR meters deep at `mean`
 * and 35 psu over the classes at their subtropical depths. Before every
 * step each cell's surface temperature takes F = F₁ sin(2πt/year) −
 * λ(T − mean), t at the middle of the step from the March equinox, at the
 * capacity the ocean publishes, as the coupled model's surface update
 * does, and every edge carries the same eastward stress. With `currents`
 * false the edge velocities are set back to rest after every step, so
 * each column keeps to itself and only the mixed layer's exchanges and
 * the restoring of the classes change it; the column read is the cell
 * nearest `latitude`. The slabs take the same F at a fixed capacity from
 * the same start. Temperatures are kelvin, depths meters below the
 * surface, heat J/m² relative to the start, times model days.
 */
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { cellVector } from '../js/dynamics/operators.module.js';
import { CP_DRY, GRAVITY } from '../js/dynamics/sigmaCore.module.js';
import { createOcean, EPS } from '../js/ocean/layered.module.js';

export const DAY = 86400, YEAR_DAYS = 365, STEP = 5400, START_T = 288.15, FLOOR = 50, CAP = 600;
export const RHO_CP = 1025 * 3985;
export const AIR_COLUMN = CP_DRY * 1e5 / GRAVITY;
export const SLAB = 4.2e7;
export const SEAS = { slab: SLAB, fixed: FLOOR * RHO_CP };
export const HEATING_PEAK = YEAR_DAYS / 4;
const OMEGA = 2 * Math.PI / (YEAR_DAYS * DAY);

let mesh = null;
export const columnMesh = () => (mesh ??= buildMesh(new Grid(2)));

export function eastwardStress(m, tau) {
  const stress = new Float64Array(m.nEdges);
  for (let e = 0; e < m.nEdges; e++) {
    const x = m.xEdge[3 * e], y = m.xEdge[3 * e + 1], r = Math.hypot(x, y);
    stress[e] = r > 0 ? tau * (-y * m.nEdge[3 * e] + x * m.nEdge[3 * e + 1]) / r : 0;
  }
  return stress;
}

export function nearestCell(m, latitude) {
  let best = 0;
  for (let i = 0; i < m.nCells; i++) if (Math.abs(m.latCell[i] * 180 / Math.PI - latitude) < Math.abs(m.latCell[best] * 180 / Math.PI - latitude) - 1e-9) best = i;
  return best;
}

export function slabCycle({ capacity, heating = 150, damping = 20 }) {
  const amplitude = heating / Math.hypot(damping, capacity * OMEGA);
  return { amplitude, range: 2 * amplitude, lag: Math.atan2(capacity * OMEGA, damping) / OMEGA / DAY };
}

export function createSeasonalColumn({
  latitude = 32, heating = 150, wind = 0.05, damping = 20, sea = 'model', dt = STEP, days = YEAR_DAYS, start = 0, mean = START_T,
  currents = false, every = 1, storeDepth = CAP, ocean: oceanOptions = {},
} = {}) {
  if (sea !== 'model' && !(sea in SEAS)) throw new Error(`sea is 'model', 'slab' or 'fixed', not ${sea}`);
  const steps = Math.round(days * DAY / dt), samples = Math.floor(steps / every) + 1;
  const m = columnMesh(), C = m.nCells, cell = nearestCell(m, latitude), live = sea === 'model';
  const ocean = live ? createOcean(m, { everySteps: 1, thermoclineTilt: 0, salinityProfile: () => 35, mixedDepth: FLOOR, ...oceanOptions }) : null;
  const L = live ? ocean.layers : 0;
  const surfaceT = new Float64Array(C).fill(mean), ice = new Float64Array(C), flux = new Float64Array(C);
  const stress = eastwardStress(m, wind);
  if (live) ocean.initialize(surfaceT, ice);
  const names = [...Object.keys(SEAS), ...(live ? ['model'] : [])];
  const slabT = Object.fromEntries(Object.keys(SEAS).map((name) => [name, mean]));
  const heatIn = Object.fromEntries(names.map((name) => [name, 0])), applied = Object.fromEntries(names.map((name) => [name, 0]));

  const record = {
    time: new Float64Array(samples), forcing: new Float32Array(samples),
    seas: Object.fromEntries(names.map((name) => [name, { sst: new Float32Array(samples), heat: new Float64Array(samples), flux: new Float32Array(samples) }])),
    ...(live ? {
      mixedDepth: new Float32Array(samples), salinity: new Float32Array(samples), stored: new Float64Array(samples),
      exchange: new Float32Array(samples), exchangeClass: new Int8Array(samples).fill(-1),
      bases: new Float32Array(samples * L), temperatures: new Float32Array(samples * L),
    } : {}),
  };
  const columnHeat = () => { let q = 0; for (let k = 0; k < L; k++) q += ocean.Q[k * C + cell]; return RHO_CP * q; };
  const startHeat = live ? columnHeat() : 0;
  let step = 0, sample = 0, forcing = 0;

  function content(n, depth) {
    let q = 0, top = 0;
    for (let k = 0; k < L && top < depth; k++) {
      const bottom = record.bases[n * L + k];
      q += (Math.min(bottom, depth) - top) * record.temperatures[n * L + k];
      top = bottom;
    }
    return RHO_CP * q;
  }
  const heatAbove = (n, depth) => content(n, depth) - content(0, depth);
  const thickness = (n, k) => record.bases[n * L + k] - (k > 0 ? record.bases[n * L + k - 1] : 0);

  function write() {
    const n = sample++;
    record.time[n] = start + step * dt / DAY;
    record.forcing[n] = forcing;
    for (const name of Object.keys(SEAS)) {
      const s = record.seas[name];
      s.sst[n] = slabT[name]; s.heat[n] = SEAS[name] * (slabT[name] - mean); s.flux[n] = applied[name];
    }
    if (!live) return;
    const s = record.seas.model;
    s.sst[n] = ocean.T0[cell]; s.heat[n] = columnHeat() - startHeat; s.flux[n] = applied.model;
    record.mixedDepth[n] = ocean.h[cell];
    record.salinity[n] = ocean.S0[cell];
    let z = 0;
    for (let k = 0; k < L; k++) {
      const h = ocean.h[k * C + cell];
      z += h;
      record.bases[n * L + k] = z;
      record.temperatures[n * L + k] = k === 0 ? ocean.T0[cell] : ocean.Q[k * C + cell] / h;
    }
    record.stored[n] = heatAbove(n, storeDepth);
    if (n === 0) return;
    const change = record.mixedDepth[n] - record.mixedDepth[n - 1];
    let most = -1, size = 0;
    for (let k = 1; k < L; k++) {
      const d = (thickness(n - 1, k) - thickness(n, k)) * Math.sign(change);
      if (d > size) { size = d; most = k; }
    }
    record.exchange[n] = change;
    record.exchangeClass[n] = most;
  }

  function advance() {
    if (step >= steps) return false;
    forcing = heating * Math.sin(OMEGA * (start * DAY + (step + 0.5) * dt));
    for (const name of Object.keys(SEAS)) {
      applied[name] = forcing - damping * (slabT[name] - mean);
      slabT[name] += applied[name] * dt / SEAS[name];
      heatIn[name] += applied[name] * dt;
    }
    if (live) {
      applied.model = forcing - damping * (surfaceT[cell] - mean);
      for (let i = 0; i < C; i++) surfaceT[i] += (forcing - damping * (surfaceT[i] - mean)) * dt / ocean.capacity[i];
      heatIn.model += applied.model * dt;
      ocean.advance(surfaceT, ice, flux, stress, dt);
      if (!currents) ocean.u.fill(0);
    }
    step++;
    if (step % every === 0) write();
    return step < steps;
  }

  function layersAt(n) {
    const out = [];
    let top = 0;
    for (let k = 0; k < L; k++) {
      const bottom = record.bases[n * L + k];
      if (bottom - top > 2 * EPS) out.push({ layer: k, density: k === 0 ? null : ocean.densities[k], top, bottom, temperature: record.temperatures[n * L + k] });
      top = bottom;
    }
    return out;
  }

  function temperatureAt(n, depth) {
    for (const layer of layersAt(n)) if (depth < layer.bottom) return layer.temperature;
    return NaN;
  }

  function events({ threshold = 0.05 * dt * every / 3600, upTo = sample - 1 } = {}) {
    const out = [];
    for (let n = 1; n <= upTo; n++) {
      const change = record.exchange[n];
      if (Math.abs(change) < threshold) continue;
      const kind = change > 0 ? 'entrains' : 'detrains', last = out[out.length - 1];
      if (last && last.kind === kind && last.to === record.time[n - 1]) {
        last.to = record.time[n]; last.meters += Math.abs(change);
        last.classes.add(record.exchangeClass[n]);
      } else out.push({ kind, from: record.time[n - 1], to: record.time[n], meters: Math.abs(change), classes: new Set([record.exchangeClass[n]]) });
    }
    return out.map((e) => ({ ...e, classes: [...e.classes].filter((k) => k > 0).map((k) => ocean.densities[k]) }));
  }

  function readouts(n = sample - 1, name = sea) {
    const s = record.seas[name], perYear = Math.round(YEAR_DAYS * DAY / (dt * every));
    const from = Math.max(0, n - perYear + 1);
    let low = Infinity, high = -Infinity, warmest = from;
    for (let j = from; j <= n; j++) {
      if (s.sst[j] < low) low = s.sst[j];
      if (s.sst[j] > high) { high = s.sst[j]; warmest = j; }
    }
    const peak = Math.floor((record.time[warmest] - HEATING_PEAK) / YEAR_DAYS) * YEAR_DAYS + HEATING_PEAK;
    const lag = warmest < n && peak >= record.time[from] ? record.time[warmest] - peak : null;
    const spring = Math.floor(record.time[n] / YEAR_DAYS) * YEAR_DAYS;
    let since = 0;
    while (since < n && record.time[since] < spring) since++;
    return {
      day: record.time[n], sst: s.sst[n], low, high, range: high - low, warmestDay: record.time[warmest], lag,
      storedSinceSpring: s.heat[n] - s.heat[since],
      mixedDepth: name === 'model' ? record.mixedDepth[n] : SEAS[name] / RHO_CP,
    };
  }

  write();
  return {
    cell, latitude: m.latCell[cell] * 180 / Math.PI, steps, dt, every, ocean, layers: L, mean, record, heatIn,
    densities: live ? Array.from(ocean.densities) : [],
    stirring: Math.hypot(...cellVector(m, stress).subarray(3 * cell, 3 * cell + 3)),
    step: advance,
    run(n = steps) { for (let j = 0; j < n && advance(); j++); return this; },
    get done() { return step >= steps; },
    get samples() { return sample; },
    layersAt, temperatureAt, heatAbove, events, readouts,
  };
}

export function seasonalColumn(options) {
  return createSeasonalColumn(options).run();
}
