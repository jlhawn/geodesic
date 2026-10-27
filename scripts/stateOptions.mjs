import { DEFAULT_SHORTWAVE } from '../js/physics/radiation.module.js';
import { MOIST_DEFAULTS } from '../js/physics/moist.module.js';

const SOIL_OPTIONS = ['soil', 'soilDepths', 'soilHeatCapacity', 'soilConductivity', 'firnDensity'];
export const soilOptions = (land = {}) => Object.fromEntries(SOIL_OPTIONS.filter((name) => land[name] !== undefined).map((name) => [name, land[name]]));

/*
 * The options a state's header records (js/stateFile.module.js), by group:
 * those the run was given, over the radiation's shortwave scheme, the
 * moist physics' Bergeron-Findeisen conversion and the land's soil at
 * their defaults, so that the page starts the state as it was made when
 * those defaults change. Later layers win.
 */
export function recordedOptions(...layers) {
  const defaults = { radiation: { shortwave: DEFAULT_SHORTWAVE }, moist: { bergeronAutoconversion: MOIST_DEFAULTS.bergeronAutoconversion }, land: { soil: 'slab' } };
  return Object.fromEntries(Object.entries(defaults).map(([group, values]) => [group, Object.assign({}, values, ...layers.map((layer) => layer?.[group]))]));
}
