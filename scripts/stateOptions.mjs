import { DEFAULT_SHORTWAVE } from '../js/physics/radiation.module.js';
import { MOIST_DEFAULTS } from '../js/physics/moist.module.js';
import { DEFAULT_SOIL } from '../js/physics/land.module.js';
import { DEFAULT_AEROSOL } from '../js/aerosol.module.js';

const SOIL_OPTIONS = ['soil', 'soilDepths', 'soilHeatCapacity', 'soilConductivity', 'firnDensity'];
export const soilOptions = (land = {}) => Object.fromEntries(SOIL_OPTIONS.filter((name) => land[name] !== undefined).map((name) => [name, land[name]]));

/*
 * The options a state's header records (js/stateFile.module.js), by group:
 * those the run was given, over the radiation's shortwave scheme and
 * aerosol, the moist physics' Bergeron-Findeisen conversion and PCAPE
 * minimum and the land's soil at their defaults, so that the page starts
 * the state as it was made when those defaults change. Later layers win,
 * and a group given as false is recorded as false.
 */
export function recordedOptions(...layers) {
  const defaults = { radiation: { shortwave: DEFAULT_SHORTWAVE, aerosol: DEFAULT_AEROSOL }, moist: { bergeronAutoconversion: MOIST_DEFAULTS.bergeronAutoconversion, pcapeMinimum: MOIST_DEFAULTS.pcapeMinimum }, land: { soil: DEFAULT_SOIL } };
  const groups = new Set([...Object.keys(defaults), ...layers.flatMap((layer) => Object.keys(layer ?? {}))]);
  return Object.fromEntries([...groups].map((group) => {
    const given = layers.map((layer) => layer?.[group]).filter((values) => values !== undefined);
    return [group, given.at(-1) === false ? false : Object.assign({}, defaults[group], ...given)];
  }));
}

const canonical = (value) => (value && typeof value === 'object' && !Array.isArray(value) ? `{${Object.keys(value).sort().map((key) => `${JSON.stringify(key)}:${canonical(value[key])}`).join(',')}}` : JSON.stringify(value));

/*
 * The groups whose options differ between two records, a group either
 * lacks counting as different.
 */
export const differingGroups = (a = {}, b = {}) => [...new Set([...Object.keys(a), ...Object.keys(b)])].filter((group) => canonical(a[group]) !== canonical(b[group]));
