import { DEFAULT_SHORTWAVE } from './physics/radiation.module.js';

/*
 * How often the model's slower parts run, as intervals of model time that
 * the drivers (the page's worker and the scripts) turn into steps of their
 * dt: from N=64 up the full radiation every RADIATION_MINUTES of its
 * shortwave scheme (the radiation's radiationEvery; the radiation held
 * between calls, see applyHeld in js/physics/radiation.module.js: 45
 * minutes with the layered shortwave, every 16th step at N=128 and every
 * 8th at N=64, and 11.25 minutes with the sequential one), the
 * gravity-wave and orographic drags' profiles every DRAG_MINUTES (the
 * models' dragEvery; the profiles held and applied every step between,
 * see js/model.module.js), and from N=64 up the ocean's step (its
 * everySteps, oceanMinutes()) OCEAN_MINUTES at
 * N=64 and shorter in proportion to the cell spacing up to N=128 (22.5
 * minutes there, the cadence the spin-ups ran: 45 minutes at N=128 drives
 * its currents to the 5 m/s cap within a day), then with the cube of the
 * spacing: at N=192 a 15-minute step pins 700,000 edges at the cap on its
 * first day and is NaN on its fifth where 7.5 minutes (every 4 steps)
 * holds five days with currents under 2.6 m/s, and at N=256 5.6 minutes
 * clamps 1.3 million edges on the first day where 2.8 (every 2) holds.
 * Options given in steps win; coarser grids
 * and the tests keep the engines' own defaults, the radiation and the
 * drags every step and the ocean every 4 steps. See "M24 —
 * Performance" in docs/c-grid-dynamical-core.md for what each interval
 * costs in realism and saves in time, and docs/cloud-gas-overlap.md for
 * the layered shortwave's cost and the error of holding it.
 */
export const RADIATION_MINUTES = { layered: 45, sequential: 11.25 };
export const DRAG_MINUTES = 11.25;
export const OCEAN_MINUTES = 45;

export const stepsFor = (minutes, dt) => Math.max(1, Math.round(minutes * 60 / dt));
export const oceanMinutes = (N) => (N <= 128 ? OCEAN_MINUTES * 64 / N : OCEAN_MINUTES * 0.5 * (128 / N) ** 3);

export function withCadence({ radiation = {}, ocean = {}, dragEvery } = {}, dt, N) {
  return {
    radiation: { radiationEvery: N < 64 ? 1 : stepsFor(RADIATION_MINUTES[radiation.shortwave ?? DEFAULT_SHORTWAVE] ?? RADIATION_MINUTES[DEFAULT_SHORTWAVE], dt), ...radiation },
    ocean: ocean === false || N < 64 ? ocean : { everySteps: stepsFor(oceanMinutes(N), dt), ...ocean },
    dragEvery: dragEvery ?? stepsFor(DRAG_MINUTES, dt),
  };
}
