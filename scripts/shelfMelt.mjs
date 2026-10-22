// The ice shelves' basal melt in states saved with iceShelves 'cavity',
// against Rignot et al. (2013, Science 341, 266, Table 1), over each
// interval between the states given (the first from day 0):
// - the Antarctic total and cavity area, south of 60°S, against the
//   1,325 ± 235 Gt/yr and 1.554 million km² surveyed;
// - each shelf of SHELVES (the cavity cells whose centers lie in its
//   box): cells, area, Gt/yr and meters of water a year against Table 1
//   (errors of joined shelves in quadrature);
// - the melt of each shelf's deepest third of cells by draft (by its
//   grounding lines) against its shallowest (by its front), the
//   Antarctic melt by draft, and the cells where water froze onto the base;
// - the area mean of the water 20 m below the base at the interval's end,
//   its temperature above the interface's freezing point (the thermal
//   driving);
// - the sheets' store (shelfDraw) at the interval's end, its growth and
//   its share of the melt, which is melt beyond the discharge.
//   node scripts/shelfMelt.mjs runs/shelves/shelf64_day0365.bin runs/shelves/shelf64_day0730.bin
import { readFileSync } from 'node:fs';
import { decodeState } from '../js/stateFile.module.js';
import { stateGrid } from '../js/stateGrid.module.js';
import { stateGeography, createGeography } from '../js/geography.module.js';
import { baseFreezing, SHELF_MELT } from '../js/ocean/iceShelfMelt.module.js';

export const SHELVES = [
  { name: 'Ross', lat: [-86, -77.3], lon: [[155, 180], [-180, -140]], area: 500809, melt: 47.7, error: 24.4 },
  { name: 'Filchner–Ronne', lat: [-84, -74.3], lon: [[-85, -29]], area: 443140, melt: 155.4, error: 36.4 },
  { name: 'Amery', lat: [-73.6, -68.5], lon: [[66, 75]], area: 60654, melt: 35.5, error: 23 },
  { name: 'Larsen C', lat: [-69.6, -66], lon: [[-66, -59]], area: 46465, melt: 20.7, error: 67 },
  { name: 'Getz', lat: [-75.5, -73.4], lon: [[-136, -114]], area: 34018, melt: 144.9, error: 14 },
  { name: 'Amundsen (Pine Island, Thwaites, Crosson, Dotson)', lat: [-76, -74], lon: [[-114, -98]], area: 20780, melt: 282.4, error: 12 },
];
const SURVEYED = { area: 1553978, melt: 1325, error: 235 }, YEAR = 365, LAYER = 20, GRAVITY = 9.81, DRAFTS = [0, 200, 400, 600, 800, Infinity];

const files = process.argv.slice(2);
if (!files.length) throw new Error('usage: node scripts/shelfMelt.mjs STATE...');
let mesh = null, geography = null, C = 0, cells = null, south = null, regions = null;
const states = [];
for (const path of files) {
  const st = await decodeState(new Uint8Array(readFileSync(path)));
  if (!st.ocean?.shelfMelted) throw new Error(`${path} holds no shelf melt: it was not run with iceShelves 'cavity'`);
  if (!mesh) {
    const { topography, geography: options } = stateGeography(st);
    mesh = stateGrid(st, topography).mesh; C = mesh.nCells;
    geography = createGeography(mesh, topography, options);
    cells = [...geography.cavity.keys()].filter((i) => geography.cavity[i]);
    south = cells.filter((i) => mesh.latCell[i] < -Math.PI / 3);
    const where = (i) => { const lat = mesh.latCell[i] * 180 / Math.PI, lon = ((mesh.lonCell[i] * 180 / Math.PI + 540) % 360) - 180; return { lat, lon }; };
    regions = SHELVES.map((shelf) => ({ ...shelf, cells: south.filter((i) => { const { lat, lon } = where(i); return lat >= shelf.lat[0] && lat <= shelf.lat[1] && shelf.lon.some(([a, b]) => lon >= a && lon <= b); }) }));
  } else if (st.ocean.shelfMelted.length !== C) throw new Error(`${path} is on another grid`);
  const { h, T, S, shelfLoad } = st.ocean, L = h.length / C, driving = new Float64Array(C);
  for (const i of cells) {
    let slab = 0, t = 0, s = 0;
    for (let k = 0; k < L && slab < LAYER; k++) { const n = k * C + i, part = Math.min(h[n], LAYER - slab); if (part > 0) { t += part * T[n]; s += part * S[n]; slab += part; } }
    driving[i] = t / slab - (SHELF_MELT.liquidusSlope * s / slab + baseFreezing(GRAVITY * shelfLoad[i]));
  }
  states.push({ path, day: st.day, melted: Float64Array.from(st.ocean.shelfMelted), draw: Float64Array.from(st.ocean.shelfDraw), driving });
}

const area = (list) => list.reduce((sum, i) => sum + mesh.areaCell[i], 0);
const sum = (list, field) => list.reduce((total, i) => total + mesh.areaCell[i] * field[i], 0);
const fixed = (x, digits = 1) => (Number.isFinite(x) ? x.toFixed(digits) : '–');
console.log(`${files.length} state(s), N=${Math.round(Math.sqrt((C - 2) / 10))}: ${cells.length} cavity cells, ${south.length} south of 60°S over ${fixed(area(south) / 1e12, 3)} million km² (${fixed(SURVEYED.area / 1e6, 3)} surveyed)`);
for (const { name, cells: list, area: observed } of regions) console.log(`  ${name}: ${list.length} cells, ${fixed(area(list) / 1e9, 0)} thousand km² (${fixed(observed / 1e3, 0)} observed)`);

let before = { day: 0, melted: new Float64Array(C), draw: new Float64Array(C) };
for (const after of states) {
  const years = (after.day - before.day) / YEAR, delta = new Float64Array(C);
  for (const i of cells) delta[i] = after.melted[i] - before.melted[i];
  const gt = (list) => sum(list, delta) / 1e12 / years, water = (list) => sum(list, delta) / area(list) / 1000 / years;
  const total = gt(south), drawn = (sum(south, after.draw) - sum(south, before.draw)) / 1e12;
  console.log(`\ndays ${before.day}–${after.day} (${after.path}):`);
  console.log(`  Antarctic melt ${fixed(total, 0)} Gt/yr against ${SURVEYED.melt} ± ${SURVEYED.error}, ${fixed(water(south), 2)} m/yr of water against 0.85; north of 60°S ${fixed(gt(cells.filter((i) => !south.includes(i))), 1)} Gt/yr`);
  console.log(`  store drawn ${fixed(sum(south, after.draw) / 1e12, 0)} Gt in all, ${fixed(drawn / years, 0)} Gt/yr over the interval, ${fixed(100 * drawn / (total * years), 0)} % of the melt`);
  console.log('  shelf | Gt/yr (Rignot) | m/yr water (Rignot) | deepest third, shallowest third m/yr | freezing cells | thermal driving K');
  for (const { name, cells: list, area: observed, melt, error } of regions) {
    if (!list.length) { console.log(`  ${name} | no cavity cells`); continue; }
    const byDraft = [...list].sort((a, b) => geography.iceDraft[b] - geography.iceDraft[a]), third = Math.max(1, Math.floor(list.length / 3));
    const deep = byDraft.slice(0, third), shallow = byDraft.slice(-third), freezing = list.filter((i) => delta[i] < 0);
    console.log(`  ${name} | ${fixed(gt(list))} (${melt} ± ${error}) | ${fixed(water(list), 2)} (${fixed(melt * 1e3 / observed, 2)}) | ${fixed(water(deep), 2)} at ${fixed(geography.iceDraft[deep[0]], 0)}–${fixed(geography.iceDraft[deep[deep.length - 1]], 0)} m, ${fixed(water(shallow), 2)} at ${fixed(geography.iceDraft[shallow[0]], 0)}–${fixed(geography.iceDraft[shallow[shallow.length - 1]], 0)} m | ${freezing.length} of ${list.length} | ${fixed(sum(list, after.driving) / area(list), 3)}`);
  }
  const bins = DRAFTS.slice(0, -1).map((low, b) => south.filter((i) => geography.iceDraft[i] >= low && geography.iceDraft[i] < DRAFTS[b + 1]));
  console.log(`  by draft: ${bins.map((list, b) => `${DRAFTS[b]}${Number.isFinite(DRAFTS[b + 1]) ? `–${DRAFTS[b + 1]}` : '+'} m ${list.length ? `${fixed(water(list), 2)} m/yr over ${fixed(area(list) / 1e9, 0)}k km²` : 'none'}`).join('; ')}`);
  const freezing = south.filter((i) => delta[i] < 0);
  console.log(`  freezing onto the base: ${freezing.length} cells, ${fixed(area(freezing) / 1e9, 0)}k km², ${fixed(gt(freezing))} Gt/yr, mean draft ${fixed(freezing.reduce((s, i) => s + geography.iceDraft[i], 0) / Math.max(1, freezing.length), 0)} m; in ${regions.filter((r) => r.cells.some((i) => delta[i] < 0)).map((r) => r.name.split(' (')[0]).join(', ') || 'no named shelf'}`);
  before = after;
}
