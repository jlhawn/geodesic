#!/usr/bin/env node
/*
 * Writes data/references/tao-adcp-equatorial-undercurrent.json: the mean zonal current profile and the Equatorial
 * Undercurrent's core at the four equatorial TAO/TRITON ADCP moorings (0°, 165°E, 170°W, 140°W, 110°W), from PMEL's
 * daily ADCP data (ERDDAP pmelTaoDyAdcp, u_1205 in cm/s), quality 1–3 only. At each depth a calendar month's mean is the mean of its
 * days over the whole record, and the profile is the mean of the twelve months, so the record's uneven seasons do not
 * bias it; a depth enters only where every month has at least MIN_DAYS days from at least MIN_YEARS distinct years, so
 * every depth kept is a climatology of a decade or more and none rests on a few years' anomalies. The core is the profile's eastward
 * maximum, its depth and speed refined by a parabola through the maximum and its two neighbors.
 *
 *   node scripts/referenceTaoAdcp.mjs [out.json]
 */
import { writeFileSync } from 'node:fs';

const OUT = process.argv[2] ?? new URL('../data/references/tao-adcp-equatorial-undercurrent.json', import.meta.url).pathname;
const ERDDAP = 'https://data.pmel.noaa.gov/pmel/erddap/tabledap/pmelTaoDyAdcp.csv';
const STATIONS = [['0n165e', 165], ['0n170w', -170], ['0n140w', -140], ['0n110w', -110]], MIN_DAYS = 30, MIN_YEARS = 10;
let bytes = 0;

async function rows(station) {
  const url = `${ERDDAP}?time,depth,u_1205&station=%22${station}%22&QU_5205%3E=1&QU_5205%3C=3`;
  for (let attempt = 1; ; attempt++) {
    try {
      const res = await fetch(url);
      if (!res.ok) throw new Error(`${res.status}`);
      const text = await res.text();
      bytes += text.length;
      return text.trim().split('\n').slice(2).map((l) => l.split(',')).map(([t, z, u]) => ({ month: Number(t.slice(5, 7)) - 1, year: Number(t.slice(0, 4)), depth: Number(z), u: Number(u) / 100 })).filter((r) => Number.isFinite(r.u) && Math.abs(r.u) < 5);
    } catch (e) { if (attempt === 4) throw new Error(`${url}: ${e.message}`); await new Promise((r) => setTimeout(r, 10000 * attempt)); }
  }
}

const out = {
  name: 'Equatorial Undercurrent core and mean zonal current at the TAO/TRITON ADCP moorings',
  source: { dataset: 'TAO/TRITON daily ADCP currents (NOAA PMEL GTMBA), ERDDAP dataset pmelTaoDyAdcp', url: 'https://data.pmel.noaa.gov/pmel/erddap/tabledap/pmelTaoDyAdcp.html', references: ['McPhaden et al. 1998, J. Geophys. Res. 103, 14169–14240', 'Johnson et al. 2002, Prog. Oceanogr. 52, 31–61'], license: 'public (NOAA); acknowledge the GTMBA Project Office of NOAA/PMEL' },
  period: 'each mooring\'s whole ADCP record (years listed per station)', box: 'the four moorings on the equator', mask: 'quality flags 1–3 (highest, standard, lower)',
  definition: `m/s eastward and m below the surface. profile: at each depth, the mean of the twelve calendar-month means of the daily values, a depth kept only where every calendar month has at least ${MIN_DAYS} days from at least ${MIN_YEARS} distinct years, so every depth listed is climatological (each depth's \`days\` and \`years\` are its sparsest month's). core: the profile's eastward maximum, depth and speed from a parabola through it and its neighbors. The model's matching quantity is the day-mean zonal velocity of its ocean at the equator and the mooring's longitude, averaged over a whole number of years (or the same calendar months), interpolated to these depths.`,
  script: 'scripts/referenceTaoAdcp.mjs', stations: {},
};
for (const [station, longitude] of STATIONS) {
  const data = await rows(station), bins = new Map();
  for (const r of data) {
    if (!bins.has(r.depth)) bins.set(r.depth, Array.from({ length: 12 }, () => ({ s: 0, n: 0, years: new Set() })));
    const b = bins.get(r.depth)[r.month]; b.s += r.u; b.n++; b.years.add(r.year);
  }
  const profile = [...bins.entries()].filter(([, ms]) => ms.every((b) => b.n >= MIN_DAYS && b.years.size >= MIN_YEARS)).map(([depth, ms]) => ({ depth, u: Number((ms.reduce((a, b) => a + b.s / b.n, 0) / 12).toFixed(4)), days: Math.min(...ms.map((b) => b.n)), years: Math.min(...ms.map((b) => b.years.size)) })).sort((a, b) => a.depth - b.depth);
  if (profile.length < 5) throw new Error(`${station}: ${profile.length} depths`);
  let k = 0;
  for (let i = 1; i < profile.length; i++) if (profile[i].u > profile[k].u) k = i;
  let depth = profile[k].depth, speed = profile[k].u;
  if (k > 0 && k < profile.length - 1) {
    const [a, b, c] = [profile[k - 1], profile[k], profile[k + 1]], h = (c.depth - a.depth) / 2;
    if (Math.abs(b.depth - a.depth - h) < 1e-6 && Math.abs(c.depth - b.depth - h) < 1e-6) {
      const curvature = a.u - 2 * b.u + c.u, shift = curvature < 0 ? 0.5 * (a.u - c.u) / curvature : 0;
      depth = b.depth + shift * h; speed = b.u - 0.25 * (a.u - c.u) * shift;
    }
  }
  const years = [...new Set(data.map((r) => r.year))].sort((a, b) => a - b);
  out.stations[station] = { longitude, latitude: 0, years: [years[0], years.at(-1)], yearsWithData: years.length, core: { speed: Number(speed.toFixed(3)), depth: Number(depth.toFixed(1)) }, profile };
  console.log(`${station}: ${data.length} values ${years[0]}–${years.at(-1)}, ${profile.length} depths ${profile[0].depth}–${profile.at(-1).depth} m; core ${speed.toFixed(3)} m/s at ${depth.toFixed(1)} m; ${(bytes / 1e6).toFixed(0)} MB so far`);
}
writeFileSync(OUT, JSON.stringify(out, null, 1) + '\n');
