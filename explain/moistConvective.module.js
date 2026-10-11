/*
 * One column of radiativeColumn.module.js, stepped every 3 hours with its
 * vapor none, frozen or at a fixed relative humidity, and its convection
 * off, dry or moist. Radiation alone takes three 1-hour steps per 3 hours:
 * with tens of kg/m² of vapor held in the air, the 40 m lowest layer relaxes
 * by its own infrared in about an hour, and a 3-hour step oscillates without
 * bound. "Frozen" copies the vapor the radiation sees; the convecting air
 * keeps the fixed relative humidity of the lowest layer. Layers run top
 * first; moistAdiabat takes pressures from the ground up.
 */
import { makeColumn, K } from './radiativeColumn.module.js';
import { STEP, STEPS_PER_DAY, BALANCED, CALM_DAYS, lapseRate } from './radiativeConvective.module.js';
import { adiabatThrough } from './moistParcel.module.js';
import { saturationHumidity } from '../js/physics/moist.module.js';

export const SUBSTEPS = 3;

export function createWetColumn({ vapor = 'relative', convection = 'moist', ...options } = {}) {
  const column = makeColumn({ ...options, humidity: vapor === 'none' ? null : 'relative' });
  const run = { column, vapor, convection, since: 0, calm: 0, last: null, before: null, set, restart, step, settled };

  function changed() {
    if (run.last) run.before = run.last;
    run.last = null;
    run.since = 0;
    run.calm = 0;
  }

  function set({ vapor: next, convection: mode, ...changes } = {}) {
    if (mode !== undefined) run.convection = mode;
    if (next !== undefined) {
      run.vapor = next;
      changes.humidity = next === 'relative' ? 'relative' : next === 'frozen' ? Float64Array.from(column.q) : null;
    }
    if (Object.keys(changes).length) column.set(changes);
    changed();
    return run;
  }

  function restart() {
    column.restart();
    changed();
    return run;
  }

  function step() {
    if (run.convection === 'off') for (let s = 0; s < SUBSTEPS; s++) column.step(STEP / SUBSTEPS, { convection: false });
    else column.step(STEP, { convection: run.convection });
    run.since += 1 / STEPS_PER_DAY;
    const { imbalance, groundNet } = column.budget;
    run.calm = Math.abs(imbalance) < BALANCED && Math.abs(groundNet) < BALANCED ? run.calm + 1 : 0;
    if (!run.last && settled()) run.last = snapshot(column);
    return run;
  }

  function settled() {
    return run.calm >= CALM_DAYS * STEPS_PER_DAY;
  }

  if (vapor === 'frozen') column.set({ humidity: Float64Array.from(column.q) });
  return run;
}

function snapshot(column) {
  const copy = (a) => Float64Array.from(a);
  return { T: copy(column.T), height: copy(column.height), pressure: copy(column.pressure), Ts: column.Ts };
}

export function nearestLayer(column, pressure) {
  let best = 0;
  for (let k = 1; k < K; k++) if (Math.abs(column.pressure[k] - pressure) < Math.abs(column.pressure[best] - pressure)) best = k;
  return best;
}

export const troposphericLapse = (column, lower = 850e2, upper = 300e2) => lapseRate(column, nearestLayer(column, upper), nearestLayer(column, lower));

export function convectionTop(column) {
  const k = column.convectionTop;
  return k < 0 ? null : { k, pressure: column.pressure[k], height: column.height[k], temperature: column.T[k] };
}

export function parcelHumidity(column) {
  const { humidity, relativeHumidity, surfacePressure } = column.settings, p = column.pressure[K - 1];
  return humidity ? relativeHumidity(p / surfacePressure) * saturationHumidity(column.T[K - 1], p) : 0;
}

export function moistAdiabat(column, pressures, out = new Float64Array(pressures.length)) {
  return adiabatThrough(column.T[K - 1], parcelHumidity(column), column.pressure[K - 1], pressures, out);
}
