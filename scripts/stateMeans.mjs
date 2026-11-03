// The 30-day means ending at its save that a state of scripts/spinup.mjs
// carries (js/intervalMeans.module.js): the header's summary, then each field's
// area-weighted mean globally, over land, over sea and over 30S–30N, over
// the cells where it is a number:
//   node scripts/stateMeans.mjs STATE.bin [field…]
import { readFileSync } from 'node:fs';
import { decodeState } from '../js/stateFile.module.js';
import { stateGrid } from '../js/stateGrid.module.js';
import { stateGeography, createGeography } from '../js/geography.module.js';

const [file, ...asked] = process.argv.slice(2);
if (!file) { console.error('usage: node scripts/stateMeans.mjs STATE.bin [field…]'); process.exit(1); }
const saved = await decodeState(new Uint8Array(readFileSync(file)));
const group = saved.means;
if (!group || !Array.isArray(group.fields)) { console.log(`${file}: day ${saved.day}, no 30-day means`); process.exit(0); }
const { topography, geography } = stateGeography(saved), mesh = stateGrid(saved, topography).mesh, land = createGeography(mesh, topography, geography).land;
const C = mesh.nCells, deg = 180 / Math.PI, end = saved.step ? Math.floor(saved.time / 86400) : saved.time / 86400;
console.log(`${file}: day ${saved.day}${saved.step ? `, step ${saved.step}` : ''}; ${group.days ? `${group.days}-day means ending at the save (model days ${(end - group.days).toFixed(3)} to ${end.toFixed(3)}), ${group.samples} samples every ${group.every} steps` : 'no whole day of means'}${Number.isInteger(group.daySamples) ? `, with the ring and ${group.daySamples} samples of the day under way to carry on` : ''}`);
const regions = { global: () => true, land: (i) => land[i], sea: (i) => !land[i], '30S-30N': (i) => Math.abs(mesh.latCell[i] * deg) <= 30 };
console.log(`| field | ${Object.keys(regions).join(' | ')} |\n|---|${Object.keys(regions).map(() => '---|').join('')}`);
for (const name of asked.length ? asked : group.fields) {
  const values = group[name];
  if (!values) { console.log(`| ${name} | not in the state |`); continue; }
  if (values.length !== C) { console.log(`| ${name} | ${values.length / C} values a cell |`); continue; }
  const cells = Object.values(regions).map((within) => {
    let sum = 0, area = 0;
    for (let i = 0; i < C; i++) if (within(i) && Number.isFinite(values[i])) { sum += mesh.areaCell[i] * values[i]; area += mesh.areaCell[i]; }
    return area > 0 ? (sum / area).toPrecision(6) : '—';
  });
  console.log(`| ${name} | ${cells.join(' | ')} |`);
}
