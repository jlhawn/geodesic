#!/usr/bin/env node
/*
 * Writes data/references/global-ts-berkeley.json: the global mean surface temperature by calendar month over 1991–2020,
 * from Berkeley Earth's Land/Ocean record (Rohde and Hausfather 2020), the variant with air temperature over sea ice.
 * Each month is the record's absolute Jan 1951–Dec 1980 climatology for it plus the mean of its 30 monthly anomalies.
 *
 *   curl -o Land_and_Ocean_complete.txt https://berkeley-earth-temperature.s3.us-west-1.amazonaws.com/Global/Land_and_Ocean_complete.txt
 *   node scripts/referenceTs.mjs Land_and_Ocean_complete.txt [out.json]
 */
import { readFileSync, writeFileSync } from 'node:fs';

const [input, OUT = new URL('../data/references/global-ts-berkeley.json', import.meta.url).pathname] = process.argv.slice(2);
if (!input) throw new Error('usage: node scripts/referenceTs.mjs Land_and_Ocean_complete.txt [out.json]');
const lines = readFileSync(input, 'latin1').split('\n');
const head = lines.findIndex((l) => l.startsWith('% Estimated Jan 1951-Dec 1980 monthly absolute temperature:'));
const climatology = lines[head + 2].replace('%', '').trim().split(/\s+/).map(Number);
const global = Number(lines.find((l) => l.includes('Using air temperature above sea ice:') && /\d/.test(l)).match(/([\d.]+) \+\/-/)[1]);
const end = lines.findIndex((l) => l.includes('Sea Ice Temperature Inferred from Water Temperatures'));
if (climatology.length !== 12 || !climatology.every(Number.isFinite) || head < 0 || end < head) throw new Error('unexpected layout');
const sums = new Float64Array(12), counts = new Float64Array(12);
for (const l of lines.slice(head, end)) {
  const [year, month, anomaly] = l.trim().split(/\s+/).map(Number);
  if (year >= 1991 && year <= 2020 && month >= 1 && month <= 12 && Number.isFinite(anomaly)) { sums[month - 1] += anomaly; counts[month - 1]++; }
}
if (!counts.every((c) => c === 30)) throw new Error(`months counted ${counts.join(' ')}`);
const monthly = climatology.map((c, m) => Number((c + sums[m] / counts[m]).toFixed(3)));
const annual = monthly.reduce((a, b) => a + b, 0) / 12;
writeFileSync(OUT, JSON.stringify({
  name: 'Global mean surface temperature, 1991–2020 monthly climatology',
  source: { dataset: 'Berkeley Earth Land/Ocean Temperature Record (Land_and_Ocean_complete.txt; land analysis 04-Jan-2025, HadSST4 ocean)', doi: '10.5194/essd-12-3469-2020', url: 'https://berkeley-earth-temperature.s3.us-west-1.amazonaws.com/Global/Land_and_Ocean_complete.txt', references: ['Rohde and Hausfather 2020, Earth Syst. Sci. Data 12, 3469–3479'], license: 'CC BY-NC 4.0 (https://berkeleyearth.org/data/); attribute to Berkeley Earth' },
  period: '1991-01 to 2020-12', box: 'global', mask: 'none: land air temperature, sea-surface temperature over open ocean, air temperature over sea ice',
  definition: `°C, month m at index m − 1: the record's 1951–1980 absolute climatology (global ${global} ± 0.019 °C) plus the month's mean anomaly over 1991–2020. The model's matching quantity is the surface skin temperature (surfaceT) from the run's day means, area-weighted over the globe.`,
  script: 'scripts/referenceTs.mjs', monthly, annual: Number(annual.toFixed(3)),
}, null, 1) + '\n');
console.log(`1991–2020: ${monthly.join(' ')}; annual ${annual.toFixed(3)} °C (1951–1980 ${global})`);
