// ERA5's zonal-mean zonal wind over 5S-5N (41 rows, cos-weighted) on its
// pressure levels at 00 UTC on DAYS of each month of YEARS, as JSON
// {levels, series: [{date, doy, u: [per level]}]} for qboCompare.mjs.
//   YEARS=2007-2019 DAYS=1,8,15,22 LEVELS=1,2,3,5,7,10,20,30,50 OUT=data/era5_equatorial_wind.json node scripts/era5EquatorialWind.mjs
import { writeFileSync } from 'node:fs';
import { fieldRows, pool, ERA5 } from './era5Raw.mjs';
const [y0, y1] = (process.env.YEARS ?? '2007-2019').split('-').map(Number);
const DAYS = (process.env.DAYS ?? '1,8,15,22').split(',').map(Number);
const LEVELS = (process.env.LEVELS ?? '1,2,3,5,7,10,20,30,50').split(',').map(Number);
const r0 = 340, r1 = 381, NX = ERA5.NX, out = [], t0 = Date.now();
const tasks = [];
for (let y = y0; y <= y1; y++) for (let m = 1; m <= 12; m++) for (const d of DAYS) {
  const date = `${y}/${String(m).padStart(2, '0')}/${String(d).padStart(2, '0')}`;
  tasks.push(async () => {
    const u = [];
    for (const level of LEVELS) {
      const x = await fieldRows(date, 'u_component_of_wind', level, 0, r0, r1);
      let s = 0, w = 0;
      for (let r = 0; r < r1 - r0; r++) { const c = Math.cos(ERA5.latOfRow(r0 + r) * Math.PI / 180); for (let i = 0; i < NX; i++) s += c * x[r * NX + i]; w += c * NX; }
      u.push(s / w);
    }
    const doy = Math.round((Date.UTC(y, m - 1, d) - Date.UTC(y, 0, 1)) / 86400e3);
    return { date, doy, u };
  });
}
const series = await pool(tasks, 12);
process.stderr.write(`${series.length} dates in ${((Date.now() - t0) / 1000).toFixed(0)} s\n`);
writeFileSync(process.env.OUT ?? 'era5_equatorial_wind.json', JSON.stringify({ levels: LEVELS, series }));
