#!/usr/bin/env node
/*
 * Writes data/references/gpcp-precip-1991-2020.bin and its .json: GPCP v2.3's monthly precipitation climatology over
 * 1991–2020 on its 2.5° grid, read from NOAA PSL's long-term-mean file by OPeNDAP. The .bin is float32 little-endian
 * [12][72][144], mm/day, rows from 88.75°S northward and columns from 1.25°E eastward (cell centers), no missing values.
 *
 *   node scripts/referenceGpcp.mjs [out.bin]
 */
import { writeFileSync } from 'node:fs';
import { createHash } from 'node:crypto';

const OUT = process.argv[2] ?? new URL('../data/references/gpcp-precip-1991-2020.bin', import.meta.url).pathname;
const URL_BASE = 'https://psl.noaa.gov/thredds/dodsC/Datasets/gpcp/precip.mon.ltm.1991-2020.nc';
const NY = 72, NX = 144;
const das = await (await fetch(`${URL_BASE}.das`)).text();
const range = das.match(/ltm_range ([\d.]+), ([\d.]+)/).slice(1).map((d) => new Date(Date.UTC(1800, 0, 1) + Number(d) * 864e5).getUTCFullYear());
if (range[0] !== 1991 || range[1] !== 2020) throw new Error(`ltm_range ${range}`);
const dods = async (query) => { const b = Buffer.from(await (await fetch(`${URL_BASE}.dods?${query}`)).arrayBuffer()), at = b.indexOf('\nData:\n') + 7; return Array.from({ length: b.readUInt32BE(at) }, (_, k) => b.readFloatBE(at + 8 + 4 * k)); };
const values = await dods(`precip.precip%5B0:1:11%5D%5B0:1:${NY - 1}%5D%5B0:1:${NX - 1}%5D`), lat = await dods('lat'), n = values.length;
if (n !== 12 * NY * NX || lat.length !== NY) throw new Error(`${n} values, ${lat.length} latitudes`);
if (lat[0] !== -88.75 || lat[NY - 1] !== 88.75) throw new Error(`lat ${lat[0]}..${lat[NY - 1]}`);
const grid = new Float32Array(n);
for (let m = 0; m < 12; m++) for (let j = 0; j < NY; j++) for (let i = 0; i < NX; i++) {
  const x = values[(m * NY + j) * NX + i];
  if (!(x >= 0 && x < 1e3)) throw new Error(`month ${m + 1} row ${j} column ${i}: ${x}`);
  grid[(m * NY + j) * NX + i] = x;
}
const bytes = Buffer.from(grid.buffer);
writeFileSync(OUT, bytes);
let s = 0, w = 0;
for (let m = 0; m < 12; m++) for (let j = 0; j < NY; j++) { const c = Math.cos((-88.75 + 2.5 * j) * Math.PI / 180); for (let i = 0; i < NX; i++) { s += c * grid[(m * NY + j) * NX + i]; w += c; } }
writeFileSync(OUT.replace(/\.bin$/, '.json'), JSON.stringify({
  name: 'Precipitation, GPCP v2.3 1991–2020 monthly climatology on 2.5°',
  source: { dataset: 'GPCP Version 2.3 Combined Precipitation Data Set, monthly long-term mean 1991–2020 (NOAA PSL precip.mon.ltm.1991-2020.nc)', doi: '10.7289/V56971M6', url: 'https://psl.noaa.gov/data/gridded/data.gpcp.html', references: ['Adler et al. 2018, Atmosphere 9, 138'], license: 'public (NOAA PSL); cite Adler et al. and acknowledge NOAA PSL' },
  period: '1991-01 to 2020-12', box: 'global 2.5° grid', mask: 'none: a box score samples the grid on the model\'s own land or sea cells of the box (boxes of the scorecard)',
  definition: 'mm/day, float32 little-endian [month 12][row 72][column 144], rows from 88.75°S northward, columns from 1.25°E eastward. The model\'s matching quantity is the day-mean total precipitation over the box\'s model cells, the reference interpolated bilinearly to each cell center and area-weighted the same way.',
  script: 'scripts/referenceGpcp.mjs', sha256: createHash('sha256').update(bytes).digest('hex'), globalMean: Number((s / w).toFixed(3)),
}, null, 1) + '\n');
console.log(`wrote ${OUT}: ${bytes.length} bytes, global mean ${(s / w).toFixed(3)} mm/day`);
