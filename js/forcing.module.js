import { encodeState, decodeState } from './stateFile.module.js';

/*
 * One model day of what the ocean and its sea ice receive from the
 * atmosphere in a coupled run: recorded by scripts/spinup.mjs (RECORD)
 * and replayed over the ocean alone by scripts/oceanSpinup.mjs. A day is
 * one file, forcing-DDDD.bin for the model day DDDD it ends, in the
 * binary state format of stateFile.module.js: 'GCMS', the JSON header
 * { kind: 'forcing', version, N, day, time, seconds, steps, oceanSteps,
 * arrays }, then float32 arrays, each a mean over the day's `seconds`:
 *
 *   stress         edges  N/m²     normal wind stress on the ocean's top
 *                                  layer after the ice transmission, as
 *                                  each ocean step received it (zero off
 *                                  the ocean)
 *   netFlux        cells  W/m²     net heat flux into the surface (PH SFLUX)
 *   shortwave      cells  W/m²     shortwave absorbed at the surface; the
 *                                  longwave and turbulent part of netFlux
 *                                  is netFlux − shortwave
 *   shortwaveDown  cells  W/m²     shortwave reaching the surface
 *   sensible       cells  W/m²     sensible heat flux, upward
 *   evaporation    cells  kg/m²/s
 *   rain           cells  kg/m²/s  precipitation, liquid or frozen
 *   snowfall       cells  kg/m²/s  the part of it that fell as snow on
 *                                  sea cells, judged by the lowest air
 *                                  temperature after each step (zero on
 *                                  land)
 *   runoff         cells  kg/m²/s  runoff leaving each land cell, before
 *                                  the ocean routes it to its outlet
 *   surfaceT       cells  K        skin temperature
 *   sst            cells  K        the ocean's surface temperature:
 *                                  surfaceT over open water, freezing
 *                                  under ice
 *   ice            cells  m        sea-ice thickness over the iced part
 *   concentration  cells  1        iced fraction of the cell
 *
 * `steps` atmosphere steps went into the cell means and `oceanSteps`
 * ocean steps into the stress.
 */
export const FORCING_VERSION = 1;
export const FORCING_FIELDS = [
  ['stress', 'edges'], ['netFlux', 'cells'], ['shortwave', 'cells'], ['shortwaveDown', 'cells'], ['sensible', 'cells'],
  ['evaporation', 'cells'], ['rain', 'cells'], ['snowfall', 'cells'], ['runoff', 'cells'],
  ['surfaceT', 'cells'], ['sst', 'cells'], ['ice', 'cells'], ['concentration', 'cells'],
];

export const forcingName = (day) => `forcing-${String(day).padStart(4, '0')}.bin`;
export function forcingDay(name) {
  const match = name.match(/^forcing-(\d+)\.bin$/);
  return match ? Number(match[1]) : null;
}

const sizes = (N) => ({ cells: 10 * N * N + 2, edges: 30 * N * N });

function check(N, fields) {
  const expected = sizes(N);
  for (const [name, where] of FORCING_FIELDS) {
    if (!fields[name]) throw new Error(`forcing lacks ${name}`);
    if (fields[name].length !== expected[where]) throw new Error(`forcing ${name} has ${fields[name].length} values, not the ${expected[where]} ${where} of N=${N}`);
  }
}

export function encodeForcing({ N, day, time, seconds, steps, oceanSteps }, fields) {
  check(N, fields);
  const arrays = Object.fromEntries(FORCING_FIELDS.map(([name]) => [name, fields[name]]));
  return encodeState({ kind: 'forcing', version: FORCING_VERSION, N, day, time, seconds, steps, oceanSteps, ...arrays });
}

export async function decodeForcing(bytes) {
  const { kind, version, N, day, time, seconds, steps, oceanSteps, ...rest } = await decodeState(bytes);
  if (kind !== 'forcing') throw new Error('not a forcing file');
  if (version !== FORCING_VERSION) throw new Error(`forcing version ${version}, not ${FORCING_VERSION}`);
  const fields = Object.fromEntries(FORCING_FIELDS.map(([name]) => [name, rest[name]]));
  check(N, fields);
  return { N, day, time, seconds, steps, oceanSteps, fields };
}
