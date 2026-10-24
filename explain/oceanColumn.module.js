/*
 * Prescribed seasons over one column of the model's own ocean: createOcean
 * on the 42 cells of Grid(2), with all its layers, a flat 4000 m floor and
 * the same start in every column, the model's analytic start at `mean`,
 * every cell placed at `latitude` (its latitude and Coriolis parameter, the
 * only things the start and the column physics read of a cell's place) and
 * the column read the cell of Grid(2) nearest CELL_LATITUDE, away from the
 * poles, where the eastward stress is zero.
 * Before every step each cell's surface temperature takes F =
 * F₁ sin(2πt/year) − λ(T − T_start), t at the middle of the step from the
 * March equinox, at the capacity the ocean publishes, as the coupled
 * model's surface update does, and every edge carries the same eastward
 * stress. With `light` the heating is sunlight: the share the ocean
 * publishes as passing its top layer is handed to it as `penetrated`, to
 * be deposited over the layers beneath, and the rest warms the surface.
 * With `currents` false the edge velocities are set back to rest after
 * every step, so each column keeps to itself and only the column physics
 * change it. The slab takes the same F at a fixed capacity from the same
 * start. Temperatures are kelvin, depths meters below the surface, heat
 * J/m² relative to the start, times model days.
 */
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { cellVector } from '../js/dynamics/operators.module.js';
import { CP_DRY, GRAVITY } from '../js/dynamics/sigmaCore.module.js';
import { createOcean } from '../js/ocean/layered.module.js';
import { PUBLISHED_OCEAN } from './published.module.js';
import { MODE_DEPTH } from '../js/ocean/column.module.js';

export const DAY = 86400, YEAR_DAYS = 365, STEP = 5400, START_T = 288.15, STORE = 100, CELL_LATITUDE = 32;
export const RHO_CP = 1025 * 3985;
export const AIR_COLUMN = CP_DRY * 1e5 / GRAVITY;
export const SLAB = 4.2e7;
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

export function atLatitude(m, latitude) {
  const phi = latitude * Math.PI / 180;
  return { ...m, latCell: new Float64Array(m.nCells).fill(phi), fCell: new Float64Array(m.nCells).fill(2 * m.omega * Math.sin(phi)) };
}

export function slabCycle({ capacity, heating = 150, damping = 20 }) {
  const amplitude = heating / Math.hypot(damping, capacity * OMEGA);
  return { amplitude, range: 2 * amplitude, lag: Math.atan2(capacity * OMEGA, damping) / OMEGA / DAY };
}

export function createSeasonalColumn({
  latitude = 32, heating = 150, wind = 0.05, damping = 20, sea = 'model', dt = STEP, days = YEAR_DAYS, start = 0, mean = START_T,
  currents = false, every = 1, storeDepth = STORE, ocean: oceanOptions = {}, light = false,
} = {}) {
  if (sea !== 'model' && sea !== 'slab') throw new Error(`sea is 'model' or 'slab', not ${sea}`);
  const steps = Math.round(days * DAY / dt), samples = Math.floor(steps / every) + 1;
  const grid = columnMesh(), cell = nearestCell(grid, CELL_LATITUDE), m = atLatitude(grid, latitude), C = m.nCells, live = sea === 'model';
  const ocean = live ? createOcean(m, { ...PUBLISHED_OCEAN, everySteps: 1, ...oceanOptions }) : null;
  const L = live ? ocean.layers : 0, minimum = live ? Float32Array.from(ocean.table.minimum) : new Float32Array(0);
  const surfaceT = new Float64Array(C).fill(mean), ice = new Float64Array(C), flux = new Float64Array(C);
  const stress = eastwardStress(m, wind);
  if (live) { ocean.initialize(surfaceT, ice); surfaceT.set(ocean.T0); }
  const reference = Float64Array.from(surfaceT), names = ['slab', ...(live ? ['model'] : [])], slabStart = reference[cell];
  let slab = slabStart;
  const heatIn = Object.fromEntries(names.map((name) => [name, 0])), applied = Object.fromEntries(names.map((name) => [name, 0]));

  const record = {
    time: new Float64Array(samples), forcing: new Float32Array(samples),
    seas: Object.fromEntries(names.map((name) => [name, { sst: new Float32Array(samples), heat: new Float64Array(samples), flux: new Float32Array(samples) }])),
    ...(live ? {
      mixedDepth: new Float32Array(samples), boundaryDepth: new Float32Array(samples), salinity: new Float32Array(samples), stored: new Float64Array(samples),
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

  function write() {
    const n = sample++;
    record.time[n] = start + step * dt / DAY;
    record.forcing[n] = forcing;
    const s = record.seas.slab;
    s.sst[n] = slab; s.heat[n] = SLAB * (slab - slabStart); s.flux[n] = applied.slab;
    if (!live) return;
    const o = record.seas.model;
    o.sst[n] = ocean.T0[cell]; o.heat[n] = columnHeat() - startHeat; o.flux[n] = applied.model;
    record.mixedDepth[n] = ocean.mixedDepths[cell];
    record.boundaryDepth[n] = ocean.boundaryDepths[cell];
    record.salinity[n] = ocean.S0[cell];
    let z = 0;
    for (let k = 0; k < L; k++) {
      const h = ocean.h[k * C + cell];
      z += h;
      record.bases[n * L + k] = z;
      record.temperatures[n * L + k] = ocean.Q[k * C + cell] / h;
    }
    record.stored[n] = heatAbove(n, storeDepth);
  }

  function advance() {
    if (step >= steps) return false;
    forcing = heating * Math.sin(OMEGA * (start * DAY + (step + 0.5) * dt));
    applied.slab = forcing - damping * (slab - slabStart);
    slab += applied.slab * dt / SLAB;
    heatIn.slab += applied.slab * dt;
    if (live) {
      applied.model = forcing - damping * (surfaceT[cell] - reference[cell]);
      for (let i = 0; i < C; i++) {
        const passing = light ? forcing * ocean.transmission[i] : 0;
        surfaceT[i] += (forcing - damping * (surfaceT[i] - reference[i]) - passing) * dt / ocean.capacity[i];
        ocean.penetrated[i] += passing * dt;
      }
      heatIn.model += applied.model * dt;
      ocean.advance(surfaceT, ice, flux, stress, dt, undefined, start * DAY + (step + 1) * dt);
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
      const bottom = record.bases[n * L + k], thickness = bottom - top;
      out.push({ layer: k, top, bottom, thickness, minimum: minimum[k], mode: Math.min(1, Math.max(0, (thickness - minimum[k]) / MODE_DEPTH)), temperature: record.temperatures[n * L + k] });
      top = bottom;
    }
    return out;
  }

  function temperatureAt(n, depth) {
    for (const layer of layersAt(n)) if (depth < layer.bottom) return layer.temperature;
    return NaN;
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
      mixedDepth: name === 'model' ? record.mixedDepth[n] : SLAB / RHO_CP,
    };
  }

  write();
  return {
    cell, latitude, steps, dt, every, ocean, layers: L, minimum, mean, record, heatIn,
    stirring: Math.hypot(...cellVector(m, stress).subarray(3 * cell, 3 * cell + 3)),
    step: advance,
    run(n = steps) { for (let j = 0; j < n && advance(); j++); return this; },
    get done() { return step >= steps; },
    get samples() { return sample; },
    layersAt, temperatureAt, heatAbove, readouts,
  };
}

export function seasonalColumn(options) {
  return createSeasonalColumn(options).run();
}
