// ERA5's eastward turbulent surface stress and mean sea-level pressure over
// its sea points within 3 deg of the equator (25 rows; land-sea mask below
// 0.5; cos-weighted) in the 10-degree sections of equatorialSections.mjs,
// 120E to 80W, on DAYS of each month of YEARS, all 24 hours, as JSON
// {years, days, latitude, sections, seaPoints, tauX, slp, tauXByHour,
// slpByHour} for armMetrics.mjs. tauX is in mN/m2, the raw files' hourly
// accumulation (N m-2 s, the hour ending at the time stamp) over 3600 s,
// slp in hPa; tauX and slp are [month][section] means over the days and
// hours, the ByHour arrays [month][hour UTC][section] means over the days.
//   YEARS=2007-2019 DAYS=1,8,15,22 OUT=data/era5_equatorial_stress.json node scripts/era5EquatorialStress.mjs
import { writeFileSync } from 'node:fs';
import { surfaceRows, pool, ERA5 } from './era5Raw.mjs';
import { SECTIONS, sectionOf } from './equatorialSections.mjs';
const [y0, y1] = (process.env.YEARS ?? '2007-2019').split('-').map(Number);
const DAYS = (process.env.DAYS ?? '1,8,15,22').split(',').map(Number);
const LAT = 3, r0 = (90 - LAT) / 0.25, r1 = (90 + LAT) / 0.25 + 1, NX = ERA5.NX, NS = SECTIONS.length, t0 = Date.now();
const lsm = await surfaceRows(`${y0}/01/01`, 'land_sea_mask', 0, r0, r1);
const section = new Int8Array(lsm.length).fill(-1), weight = new Float64Array(lsm.length), total = new Float64Array(NS), seaPoints = new Array(NS).fill(0);
for (let r = 0; r < r1 - r0; r++) for (let c = 0; c < NX; c++) {
  const n = r * NX + c, s = sectionOf(0.25 * c);
  if (s < 0 || !(lsm[n] < 0.5)) continue;
  section[n] = s; weight[n] = Math.cos(ERA5.latOfRow(r0 + r) * Math.PI / 180); total[s] += weight[n]; seaPoints[s]++;
}
const sectionMeans = (x, scale) => { const sum = new Float64Array(NS); for (let n = 0; n < x.length; n++) if (section[n] >= 0) sum[section[n]] += weight[n] * x[n]; return Array.from(sum, (v, s) => scale * v / total[s]); };
const VARS = { tauX: ['eastward_turbulent_surface_stress', 1000 / 3600], slp: ['mean_sea_level_pressure', 0.01] };
const tasks = [];
for (let y = y0; y <= y1; y++) for (let m = 1; m <= 12; m++) for (const d of DAYS) for (const [key, [variable, scale]] of Object.entries(VARS)) {
  const date = `${y}/${String(m).padStart(2, '0')}/${String(d).padStart(2, '0')}`;
  tasks.push(async () => { const hours = []; for (let h = 0; h < 24; h++) hours.push(sectionMeans(await surfaceRows(date, variable, h, r0, r1), scale)); return { key, month: m - 1, hours }; });
}
const results = await pool(tasks, 32);
const round = (x) => Math.round(100 * x) / 100, out = { years: [y0, y1], days: DAYS, latitude: LAT, sections: SECTIONS, seaPoints };
for (const key of Object.keys(VARS)) {
  const byHour = Array.from({ length: 12 }, (_, m) => {
    const samples = results.filter((r) => r.key === key && r.month === m);
    return Array.from({ length: 24 }, (_, h) => Array.from({ length: NS }, (_, s) => samples.reduce((a, r) => a + r.hours[h][s], 0) / samples.length));
  });
  out[key] = byHour.map((hours) => Array.from({ length: NS }, (_, s) => round(hours.reduce((a, v) => a + v[s], 0) / 24)));
  out[`${key}ByHour`] = byHour.map((hours) => hours.map((v) => v.map(round)));
}
process.stderr.write(`${tasks.length / 2} dates in ${((Date.now() - t0) / 1000).toFixed(0)} s\n`);
writeFileSync(process.env.OUT ?? 'era5_equatorial_stress.json', JSON.stringify(out));
