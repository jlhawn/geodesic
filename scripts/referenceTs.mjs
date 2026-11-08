#!/usr/bin/env node
/*
 * Writes data/references/global-ts-era5.json: the global mean surface temperature by calendar month over 1991–2020 in
 * the model's own definition, from ERA5's monthly means (NCAR RDA ds633.1, read by OPeNDAP at every fourth 0.25° point,
 * a 1° grid of 181 × 360). The model's Ts (surfaceT) is the land's surface temperature on land, the ice's surface
 * temperature on sea ice and the top ocean layer's temperature on open sea, so a point takes ERA5's skin temperature
 * (skt) where its land-sea mask is at least 0.5 or it holds any sea ice (ci > 0), and its bulk SST (sstk) elsewhere, or skt at a sea point without one.
 * ERA5's skt is the grid box's skin, its tiles weighted by cover, so at a point with ice it is the ice's skin and the
 * open water of the leads together: the model's matching quantity there is the same cover-weighted mean, gridBoxSkin of
 * js/intervalMeans.module.js, not the ice's own surfaceT. The annual means weigh each month by its mean length over
 * 1991–2020 (February 28 + 8/30 days).
 * For comparison it also writes the same mean with 2 m air temperature (2t) in place of skt, the definition of the
 * station and Berkeley Earth records, and the gap between the two by month.
 *
 *   node scripts/referenceTs.mjs [out.json]
 */
import { writeFileSync } from 'node:fs';

const OUT = process.argv[2] ?? new URL('../data/references/global-ts-era5.json', import.meta.url).pathname;
const BASE = 'https://thredds.rda.ucar.edu/thredds/dodsC/files/g/d633001/e5.moda.an.sfc';
const CODES = { skt: ['128_235_skt', 'SKT'], t2m: ['128_167_2t', 'VAR_2T'], ci: ['128_031_ci', 'CI'], sst: ['128_034_sstk', 'SSTK'], lsm: ['128_172_lsm', 'LSM'] };
const NY = 181, NX = 360, N = NY * NX, YEARS = Array.from({ length: 30 }, (_, k) => 1991 + k);
let bytes = 0;

async function field(name, year) {
  const [code, v] = CODES[name], url = `${BASE}/${year}/e5.moda.an.sfc.${code}.ll025sc.${year}010100_${year}120100.nc.dods?${v}.${v}%5B0:1:11%5D%5B0:4:720%5D%5B0:4:1439%5D`;
  for (let attempt = 1; ; attempt++) {
    try {
      const res = await fetch(url);
      if (!res.ok) throw new Error(`${res.status}`);
      const buf = Buffer.from(await res.arrayBuffer());
      let at = buf.indexOf('\nData:\n'); if (at < 0) throw new Error('no data marker'); at += 7;
      const n = buf.readUInt32BE(at); if (n !== 12 * N || buf.length !== at + 8 + 4 * n) throw new Error(`${n} values, ${buf.length - at} bytes`);
      bytes += buf.length;
      const out = new Float32Array(n);
      for (let k = 0; k < n; k++) { const x = buf.readFloatBE(at + 8 + 4 * k); out[k] = Number.isFinite(x) && Math.abs(x) < 1e19 ? x : NaN; }
      return out;
    } catch (e) { if (attempt === 4) throw new Error(`${url}: ${e.message}`); await new Promise((r) => setTimeout(r, 5000 * attempt)); }
  }
}

const weight = Float64Array.from({ length: NY }, (_, j) => Math.cos((90 - j) * Math.PI / 180));
const lsm = (await field('lsm', 2000)).subarray(0, N);
const sums = { skin: new Float64Array(12), air: new Float64Array(12), landSkin: new Float64Array(12), landAir: new Float64Array(12), iceSkin: new Float64Array(12), iceAir: new Float64Array(12) };
const areas = { all: new Float64Array(12), land: new Float64Array(12), ice: new Float64Array(12) }, unsampled = new Float64Array(12);
for (const year of YEARS) {
  const [skt, t2m, ci, sst] = await Promise.all(['skt', 't2m', 'ci', 'sst'].map((name) => field(name, year)));
  for (let m = 0; m < 12; m++) for (let p = 0; p < N; p++) {
    const k = m * N + p, w = weight[Math.floor(p / NX)], land = lsm[p] >= 0.5, iced = !land && ci[k] > 0;
    if (!Number.isFinite(skt[k]) || !Number.isFinite(t2m[k])) throw new Error(`${year} month ${m + 1} point ${p}: no skt or 2t`);
    const open = !land && !iced && Number.isFinite(sst[k]);
    if (!land && !iced && !open) unsampled[m] += w;
    sums.skin[m] += w * (open ? sst[k] : skt[k]); sums.air[m] += w * (open ? sst[k] : t2m[k]); areas.all[m] += w;
    if (land) { sums.landSkin[m] += w * skt[k]; sums.landAir[m] += w * t2m[k]; areas.land[m] += w; }
    if (iced) { sums.iceSkin[m] += w * skt[k]; sums.iceAir[m] += w * t2m[k]; areas.ice[m] += w; }
  }
  console.log(`${year}: ${(bytes / 1e6).toFixed(0)} MB so far`);
}
const C = (x) => Number((x - 273.15).toFixed(3)), K = (x) => Number(x.toFixed(3)), months = (f) => Array.from({ length: 12 }, (_, m) => f(m));
const monthly = months((m) => C(sums.skin[m] / areas.all[m])), air = months((m) => C(sums.air[m] / areas.all[m]));
const DAYS = [31, 28 + 8 / 30, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], YEAR = DAYS.reduce((a, b) => a + b, 0);
const annual = (xs) => K(xs.reduce((a, x, m) => a + x * DAYS[m], 0) / YEAR);
writeFileSync(OUT, JSON.stringify({
  name: 'Global mean surface temperature in the model\'s definition, 1991–2020 monthly climatology',
  source: { dataset: 'ERA5, C3S/ECMWF, via NCAR RDA ds633.1: monthly means of daily means on single levels (e5.moda.an.sfc: skt, 2t, ci, sstk, lsm)', citation: 'ERA5 Reanalysis (Monthly Mean 0.25 Degree Latitude-Longitude Grid), ECMWF, NSF NCAR 2019', doi: '10.5065/P8GT-0R61', url: 'https://thredds.rda.ucar.edu/thredds/catalog/files/g/d633001/e5.moda.an.sfc/catalog.html', references: ['Hersbach et al. 2020, Q. J. R. Meteorol. Soc. 146, 1999–2049'], license: 'CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/), ECMWF\'s license for Copernicus data since 2 July 2025, which replaced the Copernicus license', changes: 'a 1991-2020 monthly climatology from a 1° subsample' },
  period: '1991-01 to 2020-12', box: 'global', mask: 'ERA5\'s land-sea mask (lsm ≥ 0.5 is land) and sea-ice cover (ci > 0 is ice) at each point',
  definition: '°C, month m at index m − 1, cos-latitude weighted over every fourth 0.25° point. monthly: skt on land and at sea points with ice, sstk on open sea, area-weighted over the globe; skt is the grid box\'s skin (at an ice point the ice and its leads by cover), so the model\'s matching quantity is means.gridBoxSkin of the run\'s day means (c·surfaceT + (1 − c)·sst on sea cells with ice, surfaceT elsewhere), not surfaceT. annual and airAnnual weigh the months by their mean lengths over 1991–2020. air: the same with 2t in place of skt (the station records\' definition), for comparison only. gap: monthly − air, the skin less the air over land and ice in the global mean. seaIce: skt and 2t at sea points with ice, both grid-box values (ice and leads together).',
  script: 'scripts/referenceTs.mjs',
  monthly, annual: annual(monthly), air, airAnnual: annual(air), gap: months((m) => K(monthly[m] - air[m])),
  land: { skin: months((m) => C(sums.landSkin[m] / areas.land[m])), air: months((m) => C(sums.landAir[m] / areas.land[m])), share: months((m) => K(areas.land[m] / areas.all[m])) },
  seaWithoutSst: months((m) => Number((unsampled[m] / areas.all[m]).toExponential(2))),
  seaIce: { skin: months((m) => C(sums.iceSkin[m] / areas.ice[m])), air: months((m) => C(sums.iceAir[m] / areas.ice[m])), share: months((m) => K(areas.ice[m] / areas.all[m])) },
}, null, 1) + '\n');
console.log(`skin ${monthly.join(' ')}; annual ${annual(monthly)} °C, air ${annual(air)}, gap ${(annual(monthly) - annual(air)).toFixed(3)} K; ${(bytes / 1e6).toFixed(0)} MB`);
