/*
 * How often the model's slower parts run, as intervals of model time that
 * the drivers (the page's worker and the scripts) turn into steps of their
 * dt: the full radiation every RADIATION_MINUTES (the radiation's
 * radiationEvery; the radiation held between calls, see applyHeld in
 * js/physics/radiation.module.js). Options given in steps win; the
 * engines' own default, the radiation every step, is for the tests. See
 * "Performance" in docs/c-grid-dynamical-core.md for what the interval
 * costs in realism and saves in time.
 */
export const RADIATION_MINUTES = 11.25;

export const stepsFor = (minutes, dt) => Math.max(1, Math.round(minutes * 60 / dt));

export function withCadence({ radiation = {} } = {}, dt) {
  return { radiation: { radiationEvery: stepsFor(RADIATION_MINUTES, dt), ...radiation } };
}
