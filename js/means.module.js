import { memoryWeight, upwellingAtDepth, VERTICAL_MEMORY } from './frames.module.js';
import { verticalWork, verticalColumn, layerFieldAtLevel, smoothCells } from './levels.module.js';
import { savedLevels } from './stateFile.module.js';

/*
 * The running means a saved state carries, which the Vertical motion and
 * Upwelling overlays read and the figures can use at once:
 *
 *   verticalMemory (K × C) the vertical velocity at the middle of each
 *     atmosphere layer, m/s, upward positive (layerVerticalVelocity),
 *     layer k of cell i at k * C + i, a running mean over VERTICAL_MEMORY;
 *   verticalClock [folds, time] how many states the mean holds and the
 *     model time in seconds of the last, from which the mean continues;
 *   ocean.upwellingMemory (L × C) the vertical velocity at the base of each
 *     ocean layer (cellUpwelling), a running mean over upwellingTimescale
 *     of the base's depth, with ocean.upwellingFolds and ocean.upwellingTime.
 *
 * A state without them starts the means afresh.
 */

/*
 * The CPU model's memory of its layers' vertical velocity: fold(pi, theta,
 * u, time) takes the state at `time` into the mean with the weight
 * memoryWeight(Δt, folds, VERTICAL_MEMORY), Δt being the time since the last
 * fold; step(…) folds on every `every`-th call and on the first after a
 * restart. The mean restarts when time goes back and when load(saved) finds
 * no mean of this model's size in the state.
 */
export function createVerticalMemory(mesh, core, { every = 1 } = {}) {
  const { K, C } = core.diagnostics, values = new Float64Array(K * C), work = verticalWork(core), column = new Float64Array(K);
  let calls = 0;
  const memory = {
    values, every, folds: 0, time: NaN,
    restart() { memory.folds = 0; memory.time = NaN; values.fill(0); },
    fold(pi, theta, u, time) {
      if (time < memory.time) memory.folds = 0;
      const dt = memory.folds > 0 && Number.isFinite(memory.time) ? time - memory.time : 0;
      const alpha = memoryWeight(dt, memory.folds, VERTICAL_MEMORY);
      for (let i = 0; i < C; i++) {
        verticalColumn(mesh, core, pi, u, theta, i, work, column);
        for (let k = 0; k < K; k++) { const n = k * C + i; values[n] = alpha === 1 ? column[k] : values[n] + alpha * (column[k] - values[n]); }
      }
      memory.folds++; memory.time = time;
    },
    step(pi, theta, u, time) {
      calls++;
      if (memory.folds === 0 || calls % every === 0) memory.fold(pi, theta, u, time);
    },
    load(saved) {
      memory.restart();
      const kept = saved && saved.verticalMemory, clock = saved && saved.verticalClock;
      if (!kept || kept.length !== K * C || !clock || !(clock[0] > 0)) return false;
      values.set(kept);
      memory.folds = clock[0]; memory.time = Number.isFinite(clock[1]) ? clock[1] : NaN;
      return true;
    },
    serialize() { return memory.folds > 0 ? { verticalMemory: Float64Array.from(values), verticalClock: Float64Array.of(memory.folds, memory.time) } : {}; },
  };
  return memory;
}

/*
 * The ocean's memory of a decoded state (decodeState of
 * stateFile.module.js) at `depth` metres, upward positive in m/s, per cell
 * (NaN on land and where the sea floor lies above the depth): linear in
 * depth between the velocities at the base of the layers, as the Upwelling
 * overlay reads it. With `upwellingDepth` the velocity is taken there in the
 * cells whose water column reaches `depth`.
 */
export function meanUpwelling(saved, depth, upwellingDepth = depth) {
  const ocean = saved && saved.ocean;
  if (!ocean || !ocean.upwellingMemory || !(ocean.upwellingFolds > 0)) throw new Error('the state carries no upwelling memory');
  const C = saved.surfaceT.length, L = ocean.h.length / C;
  if (ocean.upwellingMemory.length !== L * C) throw new Error(`the upwelling memory holds ${ocean.upwellingMemory.length} values, not ${L} layers of ${C} cells`);
  return upwellingAtDepth(L, C, ocean.h, ocean.upwellingMemory, depth, upwellingDepth);
}

/*
 * The atmosphere's memory of a decoded state at a level, 'surface' or a
 * pressure in hPa, upward positive in m/s, per cell: interpolated in ln p
 * between the layers' middles as the Vertical motion overlay reads it. The
 * overlay also averages each cell with its neighbours; given the state's
 * `mesh`, so does this.
 */
export function meanVertical(saved, level, { mesh = null } = {}) {
  if (!saved || !saved.verticalMemory || !saved.verticalClock || !(saved.verticalClock[0] > 0)) throw new Error('the state carries no vertical-velocity memory');
  const sigma = savedLevels(saved), K = sigma.length - 1, C = saved.pi.length;
  if (saved.verticalMemory.length !== K * C) throw new Error(`the vertical-velocity memory holds ${saved.verticalMemory.length} values, not ${K} layers of ${C} cells`);
  const sigmaMid = Float64Array.from({ length: K }, (_, k) => 0.5 * (sigma[k] + sigma[k + 1]));
  const field = layerFieldAtLevel(sigmaMid, K, C, saved.pi, saved.verticalMemory, level);
  return mesh ? smoothCells(mesh, field) : field;
}
