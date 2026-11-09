#!/usr/bin/env node
/*
 * Writes data/references/era5-fields-1991-2020.bin and its .json: ERA5's 1991–2020 monthly climatology of the fields the
 * scorecard's box rows compare with, from ERA5's monthly means (NCAR RDA ds633.1) read by OPeNDAP at every fourth 0.25°
 * point, a 1° grid. Monthly means of the reanalysis are means over every hour, so they match the run's day means.
 *
 *   t2m    2 m air temperature (K), 30°N–70°N: the 40–60N land row against the model's 2 m temperature
 *   msl    mean sea-level pressure (Pa), 10°S–10°N: the equatorial SLP contrast
 *   metss  mean eastward turbulent surface stress (N/m²), 10°S–10°N: the equatorial stress row
 *   lsm    ERA5's land-sea mask over both bands: a box row on the model's land cells weighs the four points around a
 *          cell by lsm, and on its sea cells by 1 − lsm, so a coastal cell takes the points of its own kind
 *
 * The .bin is float32 little-endian, the fields in the order of the .json's `fields`, each [month 12][row][column 360],
 * rows from the band's southern edge northward at 1° and columns from 0°E eastward.
 *
 *   node scripts/referenceEra5Fields.mjs [out.bin]
 */
import { writeFileSync } from 'node:fs';
import { createHash } from 'node:crypto';

const OUT = process.argv[2] ?? new URL('../data/references/era5-fields-1991-2020.bin', import.meta.url).pathname;
const BASE = 'https://thredds.rda.ucar.edu/thredds/dodsC/files/g/d633001';
const FIELDS = [
  { name: 't2m', path: 'e5.moda.an.sfc', code: '128_167_2t', v: 'VAR_2T', south: 30, north: 70, units: 'K' },
  { name: 'msl', path: 'e5.moda.an.sfc', code: '128_151_msl', v: 'MSL', south: -10, north: 10, units: 'Pa' },
  { name: 'metss', path: 'e5.moda.fc.sfc.meanflux', code: '235_041_metss', v: 'METSS', south: -10, north: 10, units: 'N/m²' },
];
const YEARS = Array.from({ length: 30 }, (_, k) => 1991 + k), NX = 360;
let bytes = 0;

async function read(path, code, v, year, south, north, months = 12) {
  const rowFrom = 4 * (90 - north), rowTo = 4 * (90 - south);
  const url = `${BASE}/${path}/${year}/${path}.${code}.ll025sc.${year}010100_${year}120100.nc.dods?${v}.${v}%5B0:1:${months - 1}%5D%5B${rowFrom}:4:${rowTo}%5D%5B0:4:1439%5D`;
  for (let attempt = 1; ; attempt++) {
    try {
      const res = await fetch(url);
      if (!res.ok) throw new Error(`${res.status}`);
      const buf = Buffer.from(await res.arrayBuffer()), at = buf.indexOf('\nData:\n') + 7, n = buf.readUInt32BE(at), rows = north - south + 1;
      if (at < 7 || n !== months * rows * NX) throw new Error(`${n} values`);
      bytes += buf.length;
      const out = new Float64Array(n);
      for (let m = 0; m < months; m++) for (let r = 0; r < rows; r++) for (let i = 0; i < NX; i++) {
        const x = buf.readFloatBE(at + 8 + 4 * ((m * rows + r) * NX + i));
        if (!(Number.isFinite(x) && Math.abs(x) < 1e19)) throw new Error(`${v} ${year} month ${m + 1}: missing value`);
        out[(m * rows + (rows - 1 - r)) * NX + i] = x;
      }
      return out;
    } catch (e) { if (attempt === 4) throw new Error(`${url}: ${e.message}`); await new Promise((r) => setTimeout(r, 5000 * attempt)); }
  }
}

const parts = [], fields = [];
for (const f of FIELDS) {
  const rows = f.north - f.south + 1, sum = new Float64Array(12 * rows * NX);
  for (const year of YEARS) { const x = await read(f.path, f.code, f.v, year, f.south, f.north); for (let k = 0; k < sum.length; k++) sum[k] += x[k]; }
  parts.push(Float32Array.from(sum, (s) => s / YEARS.length));
  fields.push({ name: f.name, units: f.units, south: f.south, north: f.north, rows, months: 12, columns: NX, source: `${f.path}.${f.code}` });
  console.log(`${f.name}: 30 years, ${(bytes / 1e6).toFixed(0)} MB so far`);
}
const lsm = await read('e5.moda.an.sfc', '128_172_lsm', 'LSM', 2000, -10, 70, 1);
parts.push(Float32Array.from(lsm));
fields.push({ name: 'lsm', units: 'fraction', south: -10, north: 70, rows: 81, months: 1, columns: NX, source: 'e5.moda.an.sfc.128_172_lsm, January 2000' });
const bytesOut = Buffer.concat(parts.map((p) => Buffer.from(p.buffer)));
writeFileSync(OUT, bytesOut);
writeFileSync(OUT.replace(/\.bin$/, '.json'), JSON.stringify({
  name: 'ERA5 1991–2020 monthly climatology of 2 m temperature (30–70N), sea-level pressure and eastward surface stress (10S–10N) on 1°',
  source: { dataset: 'ERA5, C3S/ECMWF, via NCAR RDA ds633.1: monthly means of daily means (e5.moda.an.sfc, e5.moda.fc.sfc.meanflux)', citation: 'ERA5 Reanalysis (Monthly Mean 0.25 Degree Latitude-Longitude Grid), ECMWF, NSF NCAR 2019', doi: '10.5065/P8GT-0R61', url: 'https://thredds.rda.ucar.edu/thredds/catalog/files/g/d633001/catalog.html', references: ['Hersbach et al. 2020, Q. J. R. Meteorol. Soc. 146, 1999–2049'], license: 'CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/), ECMWF\'s license for Copernicus data since 2 July 2025, which replaced the Copernicus license', changes: 'a 1991-2020 monthly climatology from a 1° subsample' },
  period: '1991-01 to 2020-12, each calendar month the mean of 30 years', box: 'per field: south to north, all longitudes', mask: 'lsm (ERA5\'s land-sea mask, January 2000) is the last field: a box row takes the model\'s own land or sea cells and samples each at its center bilinearly, the four points\' weights multiplied by lsm on a land cell and by 1 − lsm on a sea cell and renormalized (plain bilinear where they sum to zero), then averages the cells by area (scripts/references.mjs, sampleEra5)',
  definition: 'float32 little-endian, fields in this order, each [month][row][column], rows from `south` northward at 1°, columns from 0°E at 1°. Matching model quantities, all from the run\'s day means: t2m the model\'s air temperature at 2 m above its surface, not the skin and not the lowest layer (about 20 m up): the model keeps no 2 m diagnostic yet (one is being added with the surface layer), so until it lands a 40–60N land row compares its lowest layer\'s air at that layer\'s own height and is not definition-matched; msl its mean sea-level pressure; metss its eastward turbulent stress on the surface (the run\'s exchange, not a rebuilt one).',
  script: 'scripts/referenceEra5Fields.mjs', sha256: createHash('sha256').update(bytesOut).digest('hex'), fields,
}, null, 1) + '\n');
console.log(`wrote ${OUT}: ${bytesOut.length} bytes, ${(bytes / 1e6).toFixed(0)} MB read`);
