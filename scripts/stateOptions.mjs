import { DEFAULT_SHORTWAVE } from '../js/physics/radiation.module.js';
import { MOIST_DEFAULTS } from '../js/physics/moist.module.js';
import { DEFAULT_SOIL } from '../js/physics/land.module.js';
import { DEFAULT_AEROSOL } from '../js/aerosol.module.js';

const SOIL_OPTIONS = ['soil', 'soilDepths', 'soilHeatCapacity', 'soilConductivity', 'firnDensity'];
export const soilOptions = (land = {}) => Object.fromEntries(SOIL_OPTIONS.filter((name) => land[name] !== undefined).map((name) => [name, land[name]]));

const DEFAULTS = { radiation: { shortwave: DEFAULT_SHORTWAVE, aerosol: DEFAULT_AEROSOL }, moist: { bergeronAutoconversion: MOIST_DEFAULTS.bergeronAutoconversion, pcapeMinimum: MOIST_DEFAULTS.pcapeMinimum }, land: { soil: DEFAULT_SOIL } };

/*
 * The options a state's header records (js/stateFile.module.js), by group:
 * those the run was given, over the radiation's shortwave scheme and
 * aerosol, the moist physics' Bergeron-Findeisen conversion and PCAPE
 * minimum and the land's soil at their defaults, so that the page starts
 * the state as it was made when those defaults change. Later layers win.
 */
export function recordedOptions(...layers) {
  return Object.fromEntries(Object.entries(DEFAULTS).map(([group, values]) => [group, Object.assign({}, values, ...layers.map((layer) => layer?.[group]))]));
}

/*
 * Every option group a run was given, over the defaults recordedOptions
 * records, a group given as false recorded as false: a header's `given`,
 * which the scripts read and the page does not. Later layers win.
 */
export function givenOptions(...layers) {
  const groups = new Set([...Object.keys(DEFAULTS), ...layers.flatMap((layer) => Object.keys(layer ?? {}))]);
  return Object.fromEntries([...groups].map((group) => {
    const given = layers.map((layer) => layer?.[group]).filter((values) => values !== undefined);
    return [group, given.at(-1) === false ? false : Object.assign({}, DEFAULTS[group], ...given)];
  }));
}

/*
 * A spin-up header's two records of the groups it was given: `options`,
 * the radiation, the moist physics and the land's soil, which the page
 * starts the state with, and `given`, every group.
 */
export const headerOptions = ({ radiation, moist, land, ...groups } = {}) => ({ options: recordedOptions({ radiation, moist, land: soilOptions(land) }), given: givenOptions({ radiation, moist, land, ...groups }) });

const canonical = (value) => (value && typeof value === 'object' && !Array.isArray(value) ? `{${Object.keys(value).sort().map((key) => `${JSON.stringify(key)}:${canonical(value[key])}`).join(',')}}` : JSON.stringify(value));

/*
 * The groups whose options differ between two records, a group either
 * lacks counting as different.
 */
export const differingGroups = (a = {}, b = {}) => [...new Set([...Object.keys(a), ...Object.keys(b)])].filter((group) => canonical(a[group]) !== canonical(b[group]));
