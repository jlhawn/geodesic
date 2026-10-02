// The tuning score of a run and the readers it takes its numbers from:
// a spin-up log (scripts/spinup.mjs, whose day line is the day's mean), a vertical audit
// (scripts/verticalAudit.mjs) and a state's 60-90N sea-ice volume.
// Score = Σ w_i e_i², e_i = (x_i − target_i)/tolerance_i, and 10 for a
// term the run did not give (a NaN, a missing log line, state or audit row).
// TERMS are the first sweep's, TERMS2 the second's (sweep2.mjs).
import { readFileSync } from 'node:fs';
import { Grid } from '../../js/grid.module.js';
import { buildMesh } from '../../js/mesh.module.js';
import { decodeState } from '../../js/stateFile.module.js';

export const TERMS = [
  { key: 'balance', label: 'ASR - OLR, day mean (W/m2)', target: 0, tolerance: 3, weight: 4 },
  { key: 'albedo', label: 'albedo, day mean', target: 0.30, tolerance: 0.015, weight: 2 },
  { key: 'rain', label: 'global rain, audit (mm/d)', target: 2.7, tolerance: 0.2, weight: 1 },
  { key: 'sepLow', label: 'SE Pacific low cloud, radiative', target: 0.6, tolerance: 0.1, weight: 1 },
  { key: 'peruLow', label: 'Peru low cloud, radiative', target: 0.6, tolerance: 0.1, weight: 1 },
  { key: 'sepLwp', label: 'SE Pacific deck LWP as the radiation takes it (g/m2)', target: 100, tolerance: 50, weight: 1 },
  { key: 'peruLwp', label: 'Peru deck LWP as the radiation takes it (g/m2)', target: 100, tolerance: 50, weight: 1 },
  { key: 'sepRain', label: 'SE Pacific rain (mm/d)', target: 0.2, tolerance: 0.2, weight: 1 },
  { key: 'peruRain', label: 'Peru rain (mm/d)', target: 0.2, tolerance: 0.2, weight: 1 },
  { key: 'itczRain', label: 'Pacific ITCZ rain (mm/d)', target: 7.5, tolerance: 1.5, weight: 1 },
  { key: 'itczPeak', label: 'Pacific ITCZ heating peak (hPa)', target: 450, tolerance: 50, weight: 0.5 },
  { key: 'zonalPeak', label: 'zonal-mean rain peak latitude (deg N)', target: 7.5, tolerance: 2.5, weight: 0.5 },
  { key: 'stress', label: 'equatorial stress 2S-2N 160E-100W (N/m2)', target: -0.05, tolerance: 0.015, weight: 1 },
  { key: 'arctic', label: '60-90N ice loss (1e3 km3/day)', target: 0.15, tolerance: 0.03, weight: 2 },
];

export const TERMS2 = [
  { key: 'fBalance', label: 'fresh start: ASR - OLR, days 6-10 (W/m2)', target: 0, tolerance: 3, weight: 4 },
  { key: 'fAlbedo', label: 'fresh start: albedo, days 6-10', target: 0.30, tolerance: 0.015, weight: 3 },
  { key: 'fOlr', label: 'fresh start: OLR, days 6-10 (W/m2)', target: 240, tolerance: 4, weight: 1 },
  { key: 'fSwcre', label: 'fresh start: shortwave cloud effect, days 6-10 (W/m2)', target: -47, tolerance: 4, weight: 2 },
  { key: 'fLwcre', label: 'fresh start: longwave cloud effect, days 6-10 (W/m2)', target: 27, tolerance: 3, weight: 2 },
  { key: 'fRain', label: 'fresh start: global rain, days 6-10 (mm/d)', target: 2.7, tolerance: 0.2, weight: 1 },
  { key: 'balance', label: 'eight64: ASR - OLR, day 186 (W/m2)', target: 0, tolerance: 3, weight: 2 },
  { key: 'clearAlbedo', label: 'eight64: clear-sky albedo, day 186', target: 0.15, tolerance: 0.01, weight: 3 },
  { key: 'sepLow', label: 'SE Pacific low cloud, radiative', target: 0.6, tolerance: 0.1, weight: 1 },
  { key: 'peruLow', label: 'Peru low cloud, radiative', target: 0.6, tolerance: 0.1, weight: 1 },
  { key: 'sepDeckWater', label: "SE Pacific deck's own water path where it runs (g/m2)", target: 100, tolerance: 50, weight: 1 },
  { key: 'peruDeckWater', label: "Peru deck's own water path where it runs (g/m2)", target: 100, tolerance: 50, weight: 1 },
  { key: 'sepThickness', label: "SE Pacific deck's cloud-layer thickness (m)", target: 300, tolerance: 150, weight: 0.5 },
  { key: 'peruThickness', label: "Peru deck's cloud-layer thickness (m)", target: 300, tolerance: 150, weight: 0.5 },
  { key: 'sepRain', label: 'SE Pacific rain (mm/d)', target: 0.2, tolerance: 0.2, weight: 1 },
  { key: 'itczRain', label: 'Pacific ITCZ rain (mm/d)', target: 7.5, tolerance: 1.5, weight: 1 },
  { key: 'itczPeak', label: 'Pacific ITCZ heating peak (hPa)', target: 450, tolerance: 50, weight: 0.5 },
  { key: 'stress', label: 'equatorial stress 2S-2N 160E-100W (N/m2)', target: -0.05, tolerance: 0.015, weight: 0.5 },
  { key: 'arctic', label: '60-90N ice loss (1e3 km3/day)', target: 0.15, tolerance: 0.03, weight: 2 },
];

export function score(values, terms = TERMS) {
  const errors = {}, parts = {};
  let total = 0;
  for (const { key, target, tolerance, weight } of terms) {
    const x = values[key];
    const e = Number.isFinite(x) ? (x - target) / tolerance : 10;
    errors[key] = e;
    parts[key] = weight * e * e;
    total += parts[key];
  }
  return { errors, parts, total };
}

// The clear-sky albedo of a day from its ASR, albedo and shortwave cloud
// effect: the day's insolation is ASR / (1 - albedo), its clear-sky ASR
// ASR - SWCRE.
export function clearAlbedo(asr, albedo, swcre) {
  return 1 - (asr - swcre) * (1 - albedo) / asr;
}

export function readLog(file) {
  const days = [];
  let stress = null, nan = false;
  for (const line of readFileSync(file, 'utf8').split('\n')) {
    const d = line.match(/^day (\d+) .*ASR ([-\d.]+) \(atmosphere [-\d.]+\) OLR ([-\d.]+) W.*precip ([-\d.]+) mm\/d.*albedo ([-\d.]+),.*clamped (\d+)/);
    if (d) {
      const row = { day: +d[1], asr: +d[2], olr: +d[3], precip: +d[4], albedo: +d[5], clamped: +d[6] };
      const c = line.match(/SWCRE ([-\d.]+) LWCRE ([-\d.]+)/);
      if (c) Object.assign(row, { swcre: +c[1], lwcre: +c[2], clearAlbedo: clearAlbedo(row.asr, row.albedo, +c[1]) });
      const r = line.match(/clear-sky reflectance ([-\d.]+)/);
      if (r) row.clearAlbedo = +r[1];
      Object.assign(row, { meanAlbedo: row.albedo, meanAsr: row.asr, meanOlr: row.olr, meanPrecip: row.precip });
      days.push(row);
    }
    const s = line.match(/stress 160E-100W ([-\d.]+) N\/m²/);
    if (s) stress = +s[1];
    if (/NaN on day/.test(line)) nan = true;
  }
  return { days, stress, nan };
}

const AUDIT = {
  sepLow: ['SE Pacific', 'low-cloud cover, radiative'],
  peruLow: ['Peru', 'low-cloud cover, radiative'],
  sepRain: ['SE Pacific', 'rain (mm/d)'],
  peruRain: ['Peru', 'rain (mm/d)'],
  itczRain: ['Pacific ITCZ', 'rain (mm/d)'],
  itczPeak: ['Pacific ITCZ', "firing columns' convective heating peak (hPa)"],
  rain: ['', 'global rain (mm/d)'],
  evaporation: ['', 'global evaporation (mm/d)'],
  sepRuns: ['SE Pacific', 'deck runs, share of column-steps'],
  peruRuns: ['Peru', 'deck runs, share of column-steps'],
  sepInversion: ['SE Pacific', 'resolved inversion (m)'],
  peruInversion: ['Peru', 'resolved inversion (m)'],
  namibiaLow: ['Namibia', 'low-cloud cover, radiative'],
};

export function readAudit(file) {
  const lines = readFileSync(file, 'utf8').split('\n'), out = {};
  const find = (box, metric) => lines.find((l) => l.startsWith(`  ${box}`) && l.includes(`: ${metric} `) || (box === '' && l.startsWith(`  ${metric} `)));
  const first = (line, metric) => {
    if (!line) return NaN;
    const token = line.slice(line.indexOf(metric) + metric.length).trim().split(/\s+/)[0];
    return token === 'n/a' ? NaN : Number(token);
  };
  for (const [key, [box, metric]] of Object.entries(AUDIT)) out[key] = first(find(box, metric), metric);
  for (const [key, box] of [['sepLwp', 'SE Pacific'], ['peruLwp', 'Peru']]) {
    const line = find(box, "deck's liquid water path where it runs (g/m2)");
    const m = line && line.match(/capped at \d+: ([-\d.]+|n\/a)/);
    out[key] = !m || m[1] === 'n/a' ? 0 : Number(m[1]);
  }
  for (const [key, box, metric] of [['sepDeckWater', 'SE Pacific', "deck's liquid water path where it runs (g/m2)"], ['peruDeckWater', 'Peru', "deck's liquid water path where it runs (g/m2)"], ['sepThickness', 'SE Pacific', "deck's cloud-layer thickness where it runs (m)"], ['peruThickness', 'Peru', "deck's cloud-layer thickness where it runs (m)"]]) out[key] = first(find(box, metric), metric);
  const zonal = lines.find((l) => l.startsWith('  zonal-mean rain peak (mm/d)'));
  const z = zonal && zonal.match(/\[at ([-\d.]+) deg/);
  out.zonalPeak = z ? Number(z[1]) : NaN;
  out.zonalPeakRain = first(zonal, 'zonal-mean rain peak (mm/d)');
  if (!Number.isFinite(out.itczPeak)) out.itczPeak = 1000;
  for (const key of ['sepRain', 'peruRain']) if (!Number.isFinite(out[key])) out[key] = 0;
  return out;
}

// The means over days first..last of a log's day lines (balance ASR - OLR),
// NaN when a day is missing.
export function dayMeans(log, first, last) {
  const days = log.days.filter((d) => d.day >= first && d.day <= last), whole = days.length === last - first + 1 && !log.nan;
  const mean = (f) => (whole ? days.reduce((s, d) => s + f(d), 0) / days.length : NaN);
  return { balance: mean((d) => d.asr - d.olr), asr: mean((d) => d.asr), olr: mean((d) => d.olr), albedo: mean((d) => d.albedo), clearAlbedo: mean((d) => d.clearAlbedo), swcre: mean((d) => d.swcre), lwcre: mean((d) => d.lwcre), rain: mean((d) => d.precip) };
}

// The mean of several audits' readings, each key over the windows that gave
// it; the deck's own water path and thickness are where it runs, 0 when it
// ran in none of the windows.
export function meanAudits(audits) {
  const out = {};
  for (const key of new Set(audits.flatMap(Object.keys))) {
    const xs = audits.map((a) => a[key]).filter(Number.isFinite);
    out[key] = xs.length ? xs.reduce((s, x) => s + x, 0) / xs.length : /DeckWater|Thickness/.test(key) ? 0 : NaN;
  }
  return out;
}

const meshes = new Map();
export async function arcticVolume(file) {
  const s = await decodeState(new Uint8Array(readFileSync(file)));
  if (!meshes.has(s.N)) meshes.set(s.N, buildMesh(new Grid(s.N)));
  const mesh = meshes.get(s.N), deg = 180 / Math.PI;
  let volume = 0;
  for (let i = 0; i < mesh.nCells; i++) if (mesh.latCell[i] * deg >= 60 && s.ice[i] > 0) volume += mesh.areaCell[i] * s.concentration[i] * s.ice[i];
  return volume / 1e12;
}
