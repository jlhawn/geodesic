/*
 * Two columns of radiativeColumn.module.js under the same sunlight, gases
 * and ground, stepped together every 3 hours from the same start: one with
 * the ground's sensible heat and the model's dry convective adjustment, the
 * other with radiation alone. Layers run top first, as in the column.
 */
import { makeColumn, K } from './radiativeColumn.module.js';
import { GREENHOUSE_GASES } from '../js/physics/radiation.module.js';

export const STEP = 3 * 3600;
export const STEPS_PER_DAY = 86400 / STEP;
export const PPM = GREENHOUSE_GASES.carbonDioxide * 1e6;
export const BALANCED = 0.005;
export const CALM_DAYS = 30;

export function createPair({ convection = true, ...options } = {}) {
  const column = makeColumn(options), ghost = makeColumn(options);
  const pair = { column, ghost, convection, since: 0, set, restart, step, settled };

  function set({ convection, ...changes } = {}) {
    if (convection !== undefined) pair.convection = convection;
    if (Object.keys(changes).length) { column.set(changes); ghost.set(changes); }
    pair.since = 0;
    return pair;
  }

  function restart() {
    column.restart();
    ghost.restart();
    pair.since = 0;
    return pair;
  }

  function step() {
    column.step(STEP, { convection: pair.convection });
    ghost.step(STEP, { convection: false });
    pair.since += 1 / STEPS_PER_DAY;
    return pair;
  }

  function settled() {
    return pair.since >= CALM_DAYS && [column, ghost].every(({ budget }) => Math.abs(budget.imbalance) < BALANCED && Math.abs(budget.groundNet) < BALANCED);
  }

  return pair;
}

export function tropopause(column, margin = 1) {
  const { k, temperature } = column.coldPoint();
  if (k < 0) return null;
  let above = -Infinity, below = -Infinity;
  for (let j = 0; j < k; j++) above = Math.max(above, column.T[j]);
  for (let j = k + 1; j < K; j++) below = Math.max(below, column.T[j]);
  return above - temperature >= margin && below - temperature >= margin ? column.coldPoint() : null;
}

export const lapseRate = (column, upper, lower = K - 1) => (column.T[lower] - column.T[upper]) / (column.height[upper] - column.height[lower]) * 1000;

export const logScale = (v, small, large) => Math.sign(v) * Math.min(1, Math.log1p(Math.abs(v) / small) / Math.log1p(large / small));

export function pickLayers(centers, gap, first = []) {
  const chosen = [];
  const fits = (k) => chosen.every((j) => Math.abs(centers[j] - centers[k]) >= gap);
  for (const k of first) if (k >= 0 && k < centers.length && Number.isFinite(centers[k]) && !chosen.includes(k) && fits(k)) chosen.push(k);
  for (let k = centers.length - 1; k >= 0; k--) if (Number.isFinite(centers[k]) && !chosen.includes(k) && fits(k)) chosen.push(k);
  return chosen.sort((a, b) => a - b);
}

export const daysPerSecond = (since) => Math.max(25, since / 1.5);
