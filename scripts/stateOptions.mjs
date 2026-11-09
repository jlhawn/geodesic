import { DEFAULT_SHORTWAVE } from '../js/physics/radiation.module.js';
import { MOIST_DEFAULTS } from '../js/physics/moist.module.js';
import { DEFAULT_SOIL } from '../js/physics/land.module.js';
import { DEFAULT_AEROSOL } from '../js/aerosol.module.js';
import { DIAGNOSE_FIRST } from '../js/physics/boundaryLayer.module.js';
import { terrainKey, terrainOptions, terrainCache, surfaceGeopotential, rebalanceSurfacePressure } from '../js/geography.module.js';
import { savedTerrain, stateTerrain } from '../js/stateFile.module.js';

const SOIL_OPTIONS = ['soil', 'soilDepths', 'soilHeatCapacity', 'soilConductivity', 'firnDensity'];
export const soilOptions = (land = {}) => Object.fromEntries(SOIL_OPTIONS.filter((name) => land[name] !== undefined).map((name) => [name, land[name]]));

const DEFAULTS = { radiation: { shortwave: DEFAULT_SHORTWAVE, aerosol: DEFAULT_AEROSOL }, moist: { bergeronAutoconversion: MOIST_DEFAULTS.bergeronAutoconversion, pcapeMinimum: MOIST_DEFAULTS.pcapeMinimum }, land: { soil: DEFAULT_SOIL }, boundaryLayer: { diagnoseFirst: DIAGNOSE_FIRST } };

/*
 * The options a state's header records (js/stateFile.module.js), by group:
 * those the run was given, over the radiation's shortwave scheme and
 * aerosol, the moist physics' Bergeron-Findeisen conversion and PCAPE
 * minimum, the land's soil and the boundary layer's diagnoseFirst at their
 * defaults, so that the page starts the state as it was made when those
 * defaults change. Later layers win.
 */
export function recordedOptions(...layers) {
  return Object.fromEntries(Object.entries(DEFAULTS).map(([group, values]) => [group, Object.assign({}, values, ...layers.map((layer) => layer?.[group]))]));
}

/*
 * Every option group a run was given, over the defaults recordedOptions
 * records, a group given as false recorded as false: a header's `given`.
 * Later layers win.
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
 * the radiation, the moist physics, the land's soil and the boundary
 * layer, which the page starts the state with, and `given`, every group.
 */
export const headerOptions = ({ radiation, moist, land, boundaryLayer, ...groups } = {}) => {
  const options = recordedOptions({ radiation, moist, land: soilOptions(land), boundaryLayer });
  return { options: boundaryLayer === false ? { ...options, boundaryLayer: false } : options, given: givenOptions({ radiation, moist, land, boundaryLayer, ...groups }) };
};

const canonical = (value) => (value && typeof value === 'object' && !Array.isArray(value) ? `{${Object.keys(value).sort().map((key) => `${JSON.stringify(key)}:${canonical(value[key])}`).join(',')}}` : JSON.stringify(value));

/*
 * The groups whose options differ between two records, a group either
 * lacks counting as different.
 */
export const differingGroups = (a = {}, b = {}) => [...new Set([...Object.keys(a), ...Object.keys(b)])].filter((group) => canonical(a[group]) !== canonical(b[group]));

const isGroup = (value) => !!value && typeof value === 'object' && !Array.isArray(value);

/*
 * A header's `given` as its run was made: a given written while the plain
 * order was the default holds no boundaryLayer.diagnoseFirst, and reads as
 * false.
 */
export function asMade(given) {
  const group = given?.boundaryLayer;
  if (!isGroup(given) || group === false || (isGroup(group) && group.diagnoseFirst !== undefined)) return given;
  return { ...given, boundaryLayer: { ...(isGroup(group) ? group : {}), diagnoseFirst: false } };
}

/*
 * The lines a run that resumes `base` under the record `given` logs: the
 * groups whose options differ from those base's `given` records as made,
 * or, for a base that records no `given` (saved before Nov 7, 2026), that
 * the boundary layer runs diagnoseFirst over a state made in the plain
 * order.
 */
export function resumeNotes(given, base, name) {
  if (!base) return [];
  if (base.given === undefined) return given.boundaryLayer?.diagnoseFirst === true ? [`${name} records no given options (saved before Nov 7, 2026): the boundary layer now runs diagnoseFirst over a state made in the plain order`] : [];
  const other = differingGroups(given, asMade(base.given));
  return other.length ? [`options other than ${name}'s header records: ${other.join(', ')}`] : [];
}

/*
 * The option groups of a run that rebuilds a saved state: the groups its
 * header records under `given`, or `standIn`'s for a header without them,
 * with each group of `groups` over the recorded one, key by key, and one
 * given as false or as no object in its place. Throws when neither holds
 * a record, and when both do and `standIn` is not the header's record
 * group for group: the record holds the radiation, moist and land groups
 * with their defaults filled in, and the ocean's when the run had one.
 */
export function rebuiltOptions(header, groups = {}, standIn) {
  if (header.given !== undefined && standIn !== undefined) {
    const differing = differingGroups(header.given, standIn);
    if (differing.length) throw new Error(`GIVEN stands in for the given options only where the state's header records none; this header records them, and GIVEN differs from them in ${differing.join(', ')}: unset GIVEN, and give the groups to change in their own variables`);
  }
  const recorded = header.given === undefined ? standIn : asMade(header.given);
  if (!isGroup(recorded)) throw new Error('the state\'s header records no given options (no state saved before Nov 7, 2026 does), so its run cannot be rebuilt as it was made: GIVEN (JSON, {} for today\'s defaults) stands in for them');
  return Object.fromEntries([...new Set([...Object.keys(recorded), ...Object.keys(groups)])].map((group) => {
    const own = recorded[group], value = groups[group];
    return [group, value === undefined ? own : isGroup(value) && isGroup(own) ? { ...own, ...value } : value];
  }));
}

export { savedTerrain, stateTerrain };

/*
 * The terrain whose surface a state's pressure stands on (stateTerrain,
 * false for none), when a run over `terrain` resumes it and that is
 * another surface; null when the state stands on `terrain`'s own.
 */
export function terrainMove(header, terrain = true) {
  const from = stateTerrain(header);
  return terrainKey(from) !== terrainKey(terrain) ? from : null;
}

/*
 * Throws for a run over terrain of other options than the default whose
 * orography the bundled sub-grid files would serve: their fields lie below
 * the default terrain's surface, and orographyFields rebuilds them below
 * the run's own only from the fine grids of terrainCache().
 */
export function requireOwnSubgrid(terrain, orography) {
  if (terrain === false || orography === false || terrainKey(terrain ?? true) === terrainKey(true) || terrainCache()) return;
  throw new Error(`terrain ${JSON.stringify(terrain)} needs its sub-grid orography rebuilt below its own surface from fine_h5.f32 and fine_flt2.f32, which neither TERRAIN_CACHE nor a geodesic-terrain-cache beside or above the repository or the one this worktree belongs to holds; run on the default terrain, or give the grids`);
}

/*
 * A run's terrain: TERRAIN when given, else the terrain the state it
 * resumes records, undefined for the default.
 */
export const runTerrain = (given, saved) => given ?? (saved ? savedTerrain(saved) : undefined);

/*
 * Moves a resumed state's surface pressure, state[0] over state[1]'s θ,
 * from the surface it stands on onto `toPhi`, the surface of a run over
 * `terrain` (false for none), when that is another surface; returns the
 * terrain it moved from (false for none), or null when nothing moved.
 */
export function moveOntoTerrain(header, terrain, { core, mesh, geography }, state, toPhi) {
  const from = terrainMove(header, terrain);
  if (from === null) return null;
  core.diagnose(state[0], state[1], null, null);
  rebalanceSurfacePressure(core, state[0], state[1], from === false ? null : surfaceGeopotential(mesh, geography, terrainOptions(from)), toPhi);
  return from;
}
