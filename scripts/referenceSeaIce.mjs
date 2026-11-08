#!/usr/bin/env node
/*
 * Writes data/references/nsidc-sea-ice-monthly.json: the Arctic and Antarctic sea-ice extent and area by calendar month
 * over 1991–2020 from NSIDC's Sea Ice Index, version 4 (G02135), fetched from its monthly CSV files (or read from DIR,
 * where they are named as on the server). A month missing in a year (December 1987 to January 1988) is left out of its
 * mean, so each month averages 29 or 30 years; the counts are recorded.
 *
 *   node scripts/referenceSeaIce.mjs [DIR] [out.json]
 */
import { readFileSync, writeFileSync, existsSync } from 'node:fs';
import { join } from 'node:path';

const [DIR = null, OUT = new URL('../data/references/nsidc-sea-ice-monthly.json', import.meta.url).pathname] = process.argv.slice(2);
const URL_OF = (h, m) => `https://noaadata.apps.nsidc.org/NOAA/G02135/${h === 'N' ? 'north' : 'south'}/monthly/data/${h}_${String(m).padStart(2, '0')}_extent_v4.0.csv`;
const text = async (h, m) => {
  const local = DIR && join(DIR, `${h}_${String(m).padStart(2, '0')}_extent_v4.0.csv`);
  if (local && existsSync(local)) return readFileSync(local, 'utf8');
  const res = await fetch(URL_OF(h, m));
  if (!res.ok) throw new Error(`${URL_OF(h, m)}: ${res.status}`);
  return res.text();
};
const out = {
  name: 'Sea-ice extent and area by hemisphere, 1991–2020 monthly climatology',
  source: { dataset: 'NSIDC Sea Ice Index, Version 4 (G02135), monthly extent and area CSV files', doi: '10.7265/a98x-0f50', url: 'https://noaadata.apps.nsidc.org/NOAA/G02135/', references: ['Fetterer, F., K. Knowles, W. N. Meier, M. Savoie, A. K. Windnagel 2025, Sea Ice Index, Version 4, NSIDC'], license: 'no restrictions; cite the data set (https://nsidc.org/data/g02135)' },
  period: '1991-01 to 2020-12', box: 'northern and southern hemispheres', mask: 'sea ice of the ocean only: the passive-microwave product masks lakes and land, so the model\'s lake ice must be left out of its count',
  definition: 'Mkm², month m at index m − 1. Extent: the area of cells whose concentration is at least 15 %. Area: the sum of concentration times cell area (the pole hole counted as ice in the Arctic). The model\'s matching quantity is the area of sea (not lake) cells whose 30-day or monthly mean concentration is at least 15 %, from the run\'s day means.',
  script: 'scripts/referenceSeaIce.mjs', years: {}, extent: {}, area: {},
};
for (const h of ['N', 'S']) {
  out.extent[h] = []; out.area[h] = []; out.years[h] = [];
  for (let m = 1; m <= 12; m++) {
    const rows = (await text(h, m)).trim().split('\n').slice(1).map((l) => l.split(',').map((x) => x.trim())).filter((r) => +r[0] >= 1991 && +r[0] <= 2020 && +r[4] > 0 && +r[5] > 0);
    if (rows.length < 29) throw new Error(`${h} month ${m}: ${rows.length} years`);
    const mean = (k) => Number((rows.reduce((s, r) => s + +r[k], 0) / rows.length).toFixed(3));
    out.extent[h].push(mean(4)); out.area[h].push(mean(5)); out.years[h].push(rows.length);
  }
}
writeFileSync(OUT, JSON.stringify(out, null, 1) + '\n');
console.log(`Arctic March ${out.extent.N[2]}, September ${out.extent.N[8]}; Antarctic February ${out.extent.S[1]}, September ${out.extent.S[8]} Mkm² extent`);
