import { DEFAULT_SHORTWAVE } from '../js/physics/radiation.module.js';
import { MOIST_DEFAULTS } from '../js/physics/moist.module.js';

/*
 * The options a state's header records (js/stateFile.module.js), by group:
 * those the run was given, over the radiation's shortwave scheme and the
 * moist physics' Bergeron-Findeisen conversion at their defaults, so that
 * the page starts the state as it was made when those defaults change.
 * Later layers win.
 */
export function recordedOptions(...layers) {
  const defaults = { radiation: { shortwave: DEFAULT_SHORTWAVE }, moist: { bergeronAutoconversion: MOIST_DEFAULTS.bergeronAutoconversion } };
  return Object.fromEntries(Object.entries(defaults).map(([group, values]) => [group, Object.assign({}, values, ...layers.map((layer) => layer?.[group]))]));
}
