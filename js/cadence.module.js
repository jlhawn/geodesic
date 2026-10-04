/*
 * How often the model's slower parts run, as intervals of model time that
 * the drivers (the page's worker and the scripts) turn into steps of their
 * dt: the full radiation every RADIATION_MINUTES (the radiation's
 * radiationEvery; the radiation held between calls, see applyHeld in
 * js/physics/radiation.module.js), the gravity-wave and orographic drags'
 * profiles every DRAG_MINUTES (the models' dragEvery; the profiles held
 * and applied every step between, see js/model.module.js), and from N=64
 * up the ocean's step (its everySteps) OCEAN_MINUTES at N=64 and shorter
 * in proportion to the cell spacing at finer N (22.5 minutes at N=128),
 * the cadence the spin-ups ran: 45 minutes at N=128 drives its currents
 * to the 5 m/s cap within a day. Options given in steps win; coarser grids
 * and the tests keep the engines' own defaults, the radiation and the
 * drags every step and the ocean every 4 steps. See "M24 —
 * Performance" in docs/c-grid-dynamical-core.md for what each interval
 * costs in realism and saves in time.
 */
export const RADIATION_MINUTES = 11.25;
export const DRAG_MINUTES = 11.25;
export const OCEAN_MINUTES = 45;

export const stepsFor = (minutes, dt) => Math.max(1, Math.round(minutes * 60 / dt));

export function withCadence({ radiation = {}, ocean = {}, dragEvery } = {}, dt, N) {
  return {
    radiation: { radiationEvery: stepsFor(RADIATION_MINUTES, dt), ...radiation },
    ocean: ocean === false || N < 64 ? ocean : { everySteps: stepsFor(OCEAN_MINUTES * 64 / N, dt), ...ocean },
    dragEvery: dragEvery ?? stepsFor(DRAG_MINUTES, dt),
  };
}
