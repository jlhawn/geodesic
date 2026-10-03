// What the figure dumps share: reading a saved state, its tag and the
// season of its day, and a GPU model holding the state as a continuing
// spin-up segment would (scripts/spinup.mjs).
import { readFileSync } from 'node:fs';
import { basename } from 'node:path';
import { Grid } from '../../js/grid.module.js';
import { topographyFromInt16 } from '../../js/geography.module.js';
import { createGpuModel } from '../../js/gpu/model.gpu.js';
import { decodeState, savedLevels, stateName } from '../../js/stateFile.module.js';
import { savedDeckField, DECK_FIELDS, savedMoistField, MOIST_FIELDS, savedRadiationField, RADIATION_FIELDS } from '../../js/physics/regrid.module.js';

export const readTopography = () => topographyFromInt16(readFileSync(new URL('../../data/topography_0p25.bin', import.meta.url)).buffer);
export const readState = async (file) => decodeState(new Uint8Array(readFileSync(file)));
export const tagOf = (file) => basename(stateName(file)).replace(/_day\d+(?:_step\d+)?$/, '');

/*
 * The model calendar: day 0 is the March equinox, 91 the June solstice,
 * 183 the September equinox, 274 the December solstice, a year 365 days.
 */
const MARKERS = [[0, 'March equinox'], [91, 'June solstice'], [183, 'September equinox'], [274, 'December solstice'], [365, 'March equinox']];
const SEASONS = ['northern spring', 'northern summer', 'northern autumn', 'northern winter'];
export function season(day) {
  const d = ((Math.round(day) % 365) + 365) % 365, year = Math.floor(Math.round(day) / 365) + 1;
  const quarter = MARKERS.findIndex(([start], n) => d >= start && d < MARKERS[n + 1][0]);
  const [at, name] = MARKERS.reduce((best, m) => (Math.abs(d - m[0]) < Math.abs(d - best[0]) ? m : best));
  const offset = d - at;
  return `${SEASONS[quarter]}, ${name}${offset ? ` ${offset > 0 ? '+' : '−'} ${Math.abs(offset)} d` : ''}${year > 1 ? `, year ${year}` : ''}`;
}

export const figureHeader = (file, saved) => ({ tag: tagOf(file), day: saved.day, N: saved.N, season: season(saved.day) });

/*
 * A GPU model on the state's grid and level set with the run's options
 * (OCEAN, RADIATION, MOIST, BOUNDARY_LAYER, SURFACE, LAND as JSON in the
 * environment; the radiation's clear-sky pass on unless RADIATION turns
 * it off), loaded with everything a whole-day snapshot carries.
 */
export async function gpuModelFrom(saved) {
  const env = (name) => JSON.parse(process.env[name] ?? '{}');
  const model = await createGpuModel(new Grid(saved.N), {
    topography: readTopography(), levels: savedLevels(saved), ocean: env('OCEAN'), radiation: { clearSkyPass: true, ...env('RADIATION') },
    moist: env('MOIST'), boundaryLayer: env('BOUNDARY_LAYER'), surface: env('SURFACE'), land: env('LAND'),
  });
  const { state } = model;
  ['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'].forEach((name, a) => state[a].set(saved[name]));
  model.seaIce.load(state[6], saved.concentration ?? null);
  for (const field of Object.keys(DECK_FIELDS)) model.radiation[field].set(savedDeckField(saved, field, model));
  for (const field of Object.keys(MOIST_FIELDS)) model.moist[field].set(savedMoistField(saved, field, model));
  for (const field of Object.keys(RADIATION_FIELDS)) model.radiation[field].set(savedRadiationField(saved, field, model));
  if (saved.boundaryDepth) model.boundaryLayer.depth.set(saved.boundaryDepth);
  if (saved.mixingTop) model.boundaryLayer.mixingTop.set(saved.mixingTop);
  if (saved.boundaryRegime) model.boundaryLayer.regime.set(saved.boundaryRegime);
  if (saved.boundaryBuoyancy) model.boundaryLayer.buoyancyFlux.set(saved.boundaryBuoyancy);
  if (saved.windSpeed) model.surface.windSpeed.set(saved.windSpeed);
  if (saved.evaporation) model.radiation.evaporation.set(saved.evaporation);
  if (saved.exchangeHeat && model.exchange && !model.exchange.fixed) model.exchange.heat.set(saved.exchangeHeat);
  if (saved.exchangeWind && model.exchange) model.exchange.wind.set(saved.exchangeWind);
  for (const field of ['cumulusCover', 'cumulusWater']) if (saved[field] && saved[field].length === model.moist[field].length) model.moist[field].set(saved[field]);
  if (saved.subcloudVirtual && saved.subcloudVirtual.length === model.moist.subcloudVirtual.length) model.moist.subcloudVirtual.set(saved.subcloudVirtual); else model.moist.subcloudVirtual.fill(0);
  model.time = saved.time;
  model.load();
  model.ocean.load(saved.ocean, state[3], state[6]);
  model.land.load(saved.land);
  return model;
}

export const dt = (N) => 1350 * 16 / N;
export const DEG = 180 / Math.PI;

/* East and north components of the cell vectors v (three per cell) at cell i. */
export function eastNorth(mesh, v, i) {
  const lon = mesh.lonCell[i], lat = mesh.latCell[i];
  const e = [-Math.sin(lon), Math.cos(lon), 0], n = [-Math.sin(lat) * Math.cos(lon), -Math.sin(lat) * Math.sin(lon), Math.cos(lat)];
  return [v[3 * i] * e[0] + v[3 * i + 1] * e[1] + v[3 * i + 2] * e[2], v[3 * i] * n[0] + v[3 * i + 1] * n[1] + v[3 * i + 2] * n[2]];
}

export const rounded = (v, n) => (v === null || v === undefined || !Number.isFinite(v) ? null : +v.toFixed(n));
