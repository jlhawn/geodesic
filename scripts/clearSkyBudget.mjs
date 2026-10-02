// The clear-sky shortwave budget of a saved state on the CPU engine, term by
// term, by 10-degree band, by surface type and by surface class:
//   node scripts/clearSkyBudget.mjs <state.bin>
// Each column with its cloud, deck and cumulus taken away is lit at TIMES
// (48) instants spread over day DAY (the day ending DAY days after the
// equinox, as the spin-up's day lines count; the state's day by default),
// with the state's vapour, sea ice and snow held fixed, and its sunlight
// split into what leaves the top (atmos, the same column's reflection over
// a black surface, and surf, the rest: the surface's reflection as seen
// from the top), what ozone, vapour and aerosol absorb, and what the
// surface absorbs (sfcabs). Fluxes are day means over the row's area
// (W/m2); albedo is the clear-sky reflected over the incoming, surfAlb the
// sunlight the surface reflects over what reaches it, and the last column
// the surface's reflection seen from the top over surfAlb times the
// incoming. The surface types and classes are columns over that surface
// alone (scripts/clearSkyClasses.mjs: a sea cell's open water and ice are
// lit separately), with each class's reference range and verdict.
// RADIATION (JSON) passes options to the radiation, LAND (JSON) to the land.
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { createModel, STATE_NAMES } from '../js/model.module.js';
import { decodeState, savedLevels } from '../js/stateFile.module.js';
import { clearSkyClasses, classTable, TERMS, REFERENCE_SOURCES } from './clearSkyClasses.mjs';

const FILE = process.argv[2];
if (!FILE) { console.error('usage: node scripts/clearSkyBudget.mjs <state.bin>'); process.exit(1); }
const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
const TIMES = Number(process.env.TIMES ?? 48), RADIATION = JSON.parse(process.env.RADIATION ?? '{}'), LAND = JSON.parse(process.env.LAND ?? '{}'), AT = Number(process.env.DAY ?? saved.day);
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const model = createModel(new Grid(saved.N), { topography, levels: savedLevels(saved), ocean: false, radiation: { clearSkyPass: true, ...RADIATION }, land: LAND });
const { mesh, core, state, seaIce, land } = model;
const { K, C } = core.diagnostics;
STATE_NAMES.forEach((name, a) => state[a].set(saved[name]));
seaIce.load(state[6], saved.concentration ?? null);
land.load(saved.land, state[6]);
const result = clearSkyClasses(model, { day: AT, times: TIMES, ozone: RADIATION.ozoneAbsorption ?? 0.03 });

const deg = 180 / Math.PI;
const rows = new Map();
const add = (key, area, values) => {
  if (!rows.has(key)) rows.set(key, { area: 0, ...Object.fromEntries(TERMS.map((t) => [t, 0])) });
  const r = rows.get(key);
  r.area += area;
  for (const t of TERMS) r[t] += area * values[t];
};
for (let i = 0; i < C; i++) {
  const a = mesh.areaCell[i], s = Object.fromEntries(TERMS.map((t) => [t, result.cells[i][t] / TIMES]));
  add('global', a, s);
  const lat = mesh.latCell[i] * deg, b = Math.min(80, Math.floor(lat / 10) * 10);
  add(`${b >= 0 ? `${b}N` : `${-b}S`}..${b + 10 >= 0 ? `${b + 10}N` : `${-(b + 10)}S`}`.replace('0N..', '0..').replace(/^0\.\./, 'EQ..'), a, s);
  if (lat >= -30 && lat <= 30) add('30S-30N', a, s);
}
const f = (x, d = 1) => (Number.isFinite(x) ? x.toFixed(d) : 'n/a');
const line = (key, m, surfAlb, extra = '') => console.log(`${key.padEnd(11)} ${f(m('insolation')).padStart(6)} ${f(m('reflected')).padStart(5)} ${f(m('atmosphere')).padStart(6)} ${f(m('reflected') - m('atmosphere')).padStart(5)} ${f(m('ozone')).padStart(6)} ${f(m('vapour')).padStart(6)} ${f(m('aerosol')).padStart(5)} ${f(m('absorbedSurface')).padStart(6)}  ${f(m('reflected') / m('insolation'), 3)}   ${f(m('atmosphere') / m('insolation'), 3)}       ${f(surfAlb, 3)}    ${f((m('reflected') - m('atmosphere')) / (surfAlb * m('insolation')), 3)}${extra}`);
console.log(`clear-sky shortwave budget of ${FILE.split('/').pop()} (N=${saved.N}, K=${K}) lit over day ${AT} at ${TIMES} instants; RADIATION ${JSON.stringify(RADIATION)}; LAND ${JSON.stringify(LAND)}`);
console.log('row          insol  refl  atmos  surf  ozone vapour aeros sfcabs  albedo  atmos/insol  surfAlb  surf/(surfAlb*insol)');
const order = ['global', '30S-30N', ...[...rows.keys()].filter((k) => k.includes('..')).sort((x, y) => parseFloat(x) * (x.includes('S') ? -1 : 1) - parseFloat(y) * (y.includes('S') ? -1 : 1))];
for (const key of order) {
  const r = rows.get(key);
  const m = (t) => r[t] / r.area;
  line(key, m, 1 - m('absorbedSurface') / m('down'));
}
for (const r of result.typeRows) line(r.name, (t) => r.raw[t], r.surfaceAlbedo, `   (area share ${f(r.areaShare, 3)})`);
console.log('');
for (const l of classTable(result)) console.log(l);
console.log(`references: ${REFERENCE_SOURCES}; the open sea's direct-beam albedo against the range spanned by Taylor et al. (1996) and Fresnel reflection from flat water at the bin's mean cosine, widened by 0.01 on each side`);
