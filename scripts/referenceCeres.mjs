#!/usr/bin/env node
/*
 * Writes data/references/ceres-ebaf-global-monthly.json: CERES-EBAF Ed4.2.1's global monthly climatology at the top of
 * the atmosphere, 2001–2024, read from the obs4MIPs files on ESGF ORNL by OPeNDAP. For each field and calendar month the
 * 24-year mean of each 1° cell is area-weighted (cos latitude) to a global mean. Fields already in the output are kept.
 *
 *   [VARS=rsdt,rsut,...] node scripts/referenceCeres.mjs [out.json]
 */
import { writeFileSync, existsSync, readFileSync } from 'node:fs';

const OUT = process.argv[2] ?? new URL('../data/references/ceres-ebaf-global-monthly.json', import.meta.url).pathname;
const U = 'https://esgf-node.ornl.gov/thredds/dodsC/user_pub_work/obs4MIPs/NASA-LaRC/CERES-EBAF-4-2-1/mon';
const VARS = (process.env.VARS ?? 'rsdt,rsut,rsutcs,rlut,rlutcs,rstcre,rltcre').split(',');
const file = (v) => `${U}/${v}/gn/v20260220/${v}_mon_CERES-EBAF-4-2-1_NASA-GISS_gn_200003-202509.nc`;
const NY = 24, NR = 180, NC = 360, NP = NR * NC;
const axis = (text, name) => { const lines = text.split('\n'), at = lines.findIndex((l) => l.startsWith(`${name}[`)); return lines[at + 1].split(',').map(Number); };
const dateOf = (days) => new Date(Date.UTC(2000, 2, 1) + days * 86400e3);
const get = async (url) => {
  for (let attempt = 1; ; attempt++) {
    try { const res = await fetch(url); if (!res.ok) throw new Error(`${res.status}`); return res; }
    catch (e) { if (attempt === 3) throw new Error(`${url}: ${e.message}`); await new Promise((r) => setTimeout(r, 5000 * attempt)); }
  }
};
const out = existsSync(OUT) ? JSON.parse(readFileSync(OUT, 'utf8')) : {
  name: 'CERES-EBAF Ed4.2.1 TOA, global monthly climatology',
  source: { dataset: 'CERES-EBAF-TOA Edition 4.2.1 (obs4MIPs CERES-EBAF-4-2-1, NASA-GISS gn, 200003-202509)', doi: '10.5067/TERRA-AQUA-NOAA20/CERES/EBAF_L3B004.2.1', url: 'https://esgf-node.ornl.gov/thredds/catalog/user_pub_work/obs4MIPs/NASA-LaRC/CERES-EBAF-4-2-1/mon/catalog.html', references: ['Loeb et al. 2018, J. Climate 31, 895–918', 'Loeb et al. 2020, J. Climate 33, 61–75, doi:10.1175/JCLI-D-19-0381.1'], license: 'CC BY 4.0; acknowledge per https://ceres.larc.nasa.gov/' },
  period: '2001-01 to 2024-12, each calendar month the mean of 24 years', box: 'global', mask: 'none (every 1° cell)',
  definition: 'W/m², month m at index m − 1. rsdt incoming solar, rsut and rlut all-sky outgoing, rsutcs and rlutcs EBAF\'s clear-sky fluxes of cloud-free areas (clr_c), rstcre and rltcre EBAF\'s cloud radiative effects (see clearSky)',
  clearSky: 'rstcre and rltcre are EBAF\'s cloud radiative effects against its total-region clear sky (clr_t), the definition of the model\'s clear-sky pass: their annual means, -45.13 and 25.70 W/m², sit +0.48 and -2.10 from the cloud-free-area rsutcs - rsut and rlutcs - rlut (-45.61, 27.80), the adjustments Loeb et al. (2020) give (+0.5, -2.2). Score the model\'s SWCRE and LWCRE against rstcre and rltcre.',
  script: 'scripts/referenceCeres.mjs', fields: {},
};
let total = 0;
for (const v of VARS) {
  if (out.fields[v]?.length === 12 && out.fields[v].every(Number.isFinite)) continue;
  const ax = await (await get(`${file(v)}.ascii?time,lat,lon`)).text(), time = axis(ax, 'time'), lat = axis(ax, 'lat'), lon = axis(ax, 'lon');
  if (lat.length !== NR || lon.length !== NC || lat[0] !== -89.5 || lon[0] !== 0.5) throw new Error(`${v}: grid`);
  const weight = lat.map((la) => Math.cos(la * Math.PI / 180)), means = [];
  for (let m = 1; m <= 12; m++) {
    const A = 9 + m, Z = A + 12 * (NY - 1);
    for (let t = A; t <= Z; t += 12) { const d = dateOf(time[t]); if (d.getUTCMonth() !== m - 1 || d.getUTCFullYear() !== 2001 + (t - A) / 12) throw new Error(`${v}: time ${t} is ${d.toISOString()}`); }
    const buf = Buffer.from(await (await get(`${file(v)}.dods?${v}.${v}%5B${A}:12:${Z}%5D%5B0:1:${NR - 1}%5D%5B0:1:${NC - 1}%5D`)).arrayBuffer());
    let at = buf.indexOf('\nData:\n'); if (at < 0) throw new Error(`${v} ${m}: no data marker`); at += 7;
    const n = buf.readUInt32BE(at); if (n !== NY * NP || buf.length !== at + 8 + 4 * n) throw new Error(`${v} ${m}: ${n} values, ${buf.length - at} bytes`);
    let g = 0, gw = 0, minN = NY;
    for (let p = 0; p < NP; p++) {
      let s = 0, c = 0;
      for (let y = 0; y < NY; y++) { const x = buf.readFloatBE(at + 8 + 4 * (y * NP + p)); if (Number.isFinite(x) && Math.abs(x) < 1e19) { s += x; c++; } }
      minN = Math.min(minN, c);
      if (c) { const w = weight[Math.floor(p / NC)]; g += w * s / c; gw += w; }
    }
    if (minN < NY) throw new Error(`${v} ${m}: a cell has only ${minN} years`);
    means.push(g / gw); total += buf.length;
    console.log(v, m, (g / gw).toFixed(3), `${(total / 1e6).toFixed(0)} MB so far`);
  }
  out.fields[v] = means;
  writeFileSync(OUT, JSON.stringify(out, null, 1) + '\n');
}
const annual = (v) => out.fields[v] ? out.fields[v].reduce((a, b) => a + b, 0) / 12 : NaN;
console.log(`annual: rsdt ${annual('rsdt').toFixed(2)}, rsut ${annual('rsut').toFixed(2)}, rlut ${annual('rlut').toFixed(2)}, cloud-free LWCRE ${(annual('rlutcs') - annual('rlut')).toFixed(2)}, rltcre ${annual('rltcre').toFixed(2)}, cloud-free SWCRE ${(annual('rsutcs') - annual('rsut')).toFixed(2)}, rstcre ${annual('rstcre').toFixed(2)}`);
