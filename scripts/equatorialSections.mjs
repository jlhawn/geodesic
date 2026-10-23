// The equatorial Pacific in 10-degree sections from 120E to 80W, where
// era5EquatorialStress.mjs averages ERA5 (data/era5_equatorial_stress.json)
// and armMetrics.mjs the model, ERA5's values in a model day's month, as
// day means and at the day's UTC hour, and a pair's table from armMetrics
// outputs: each section's stress and sea-level pressure as the mean over a
// side's runs ± its standard error, and the difference ± its own.
//   node scripts/equatorialSections.mjs base-s0.txt base-s1.txt base-s2.txt -- fix-s0.txt fix-s1.txt fix-s2.txt
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dayOfYear } from './woaSstBias.mjs';

export const SECTIONS = Array.from({ length: 16 }, (_, s) => [120 + 10 * s, 130 + 10 * s]);
export const sectionName = ([a, b]) => `${a <= 180 ? a + 'E' : 360 - a + 'W'}-${b <= 180 ? b + 'E' : 360 - b + 'W'}`.replace('180E', '180').replace('180W', '180');
export const sectionOf = (lonDeg) => { const lon = ((lonDeg % 360) + 360) % 360; return lon >= 120 && lon < 280 ? Math.floor((lon - 120) / 10) : -1; };

export const MONTHS = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'];
const MONTH_START = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
export function monthOfDay(modelDay) {
  if (!Number.isFinite(modelDay)) throw new Error(`model day ${modelDay}`);
  const d = dayOfYear(Math.floor(modelDay));
  let m = 11;
  while (MONTH_START[m] > d) m--;
  return m;
}

// The model's day begins with the sun over Greenwich, at 12 UTC.
export function utcHourOfDay(modelDay) {
  if (!Number.isFinite(modelDay)) throw new Error(`model day ${modelDay}`);
  return Math.round(12 + 24 * (modelDay - Math.floor(modelDay))) % 24;
}

export function era5Sections(era5, modelDay) {
  const m = monthOfDay(modelDay), h = utcHourOfDay(modelDay), next = (h + 1) % 24;
  return { era5Tau: era5.tauX[m], era5Slp: era5.slp[m], era5TauAt: era5.tauXByHour[m][h].map((v, s) => (v + era5.tauXByHour[m][next][s]) / 2), era5SlpAt: era5.slpByHour[m][h] };
}

const f = (x, d) => (Number.isFinite(x) ? x.toFixed(d) : '—');
export function sectionTable({ tau, slp, era5Tau, era5Slp, era5TauAt, era5SlpAt }, hours) {
  const rows = SECTIONS.map((s, k) => `| ${sectionName(s)} | ${f(tau[k], 2)} | ${f(era5Tau[k], 2)} | ${f(tau[k] / era5Tau[k], 2)} | ${f(era5TauAt[k], 2)} | ${f(slp[k], 2)} | ${f(era5Slp[k], 2)} | ${f(era5SlpAt[k], 2)} |`);
  return [`| section | tau_x (mN/m2) | ERA5 tau_x | ratio | ERA5 tau_x at ${hours} | SLP (hPa) | ERA5 SLP | ERA5 SLP at ${hours} |`, '|---|---|---|---|---|---|---|---|', ...rows].join('\n');
}

export function readSections(text) {
  return SECTIONS.map((s) => {
    const line = text.split('\n').find((l) => l.startsWith(`| ${sectionName(s)} |`));
    if (!line) throw new Error(`no ${sectionName(s)} row`);
    const [tau, era5Tau, , era5TauAt, slp, era5Slp, era5SlpAt] = line.split('|').slice(2, -1).map(Number);
    return { tau, era5Tau, era5TauAt, slp, era5Slp, era5SlpAt };
  });
}

const stats = (xs) => { const n = xs.length, mean = xs.reduce((a, x) => a + x, 0) / n; return { mean, se: n > 1 ? Math.sqrt(xs.reduce((a, x) => a + (x - mean) ** 2, 0) / (n - 1) / n) : NaN }; };
export function pairTable(baseTexts, fixTexts) {
  const [base, fix] = [baseTexts, fixTexts].map((texts) => texts.map(readSections));
  const pm = ({ mean, se }, d) => `${f(mean, d)} ± ${f(se, d)}`, at = baseTexts[0].match(/\| ERA5 tau_x at ([^|]*?) \|/)[1];
  const table = (key, obs, d, ratio) => SECTIONS.map((s, k) => {
    const b = stats(base.map((run) => run[k][key])), x = stats(fix.map((run) => run[k][key])), diff = { mean: x.mean - b.mean, se: Math.hypot(b.se, x.se) };
    const of = (v) => (ratio ? ` | ${f(v.mean / base[0][k][obs], 2)}` : '');
    return `| ${sectionName(s)} | ${f(base[0][k][obs], d)} | ${f(base[0][k][`${obs}At`], d)} | ${pm(b, d)}${of(b)} | ${pm(x, d)}${of(x)} | ${pm(diff, d)} |`;
  });
  return [
    `tau_x (mN/m2), mean over ${base.length} base and ${fix.length} feature runs ± its standard error; ERA5's day mean, then at ${at}`,
    `| section | ERA5 | ERA5 at ${at} | base | base/ERA5 | feature | feature/ERA5 | feature - base |`, '|---|---|---|---|---|---|---|---|', ...table('tau', 'era5Tau', 1, true),
    '', 'SLP (hPa)',
    `| section | ERA5 | ERA5 at ${at} | base | feature | feature - base |`, '|---|---|---|---|---|---|', ...table('slp', 'era5Slp', 2, false),
  ].join('\n');
}

if (process.argv[1] && fileURLToPath(import.meta.url) === process.argv[1]) {
  const args = process.argv.slice(2), cut = args.indexOf('--');
  if (cut < 1 || cut === args.length - 1) throw new Error('usage: equatorialSections.mjs base.txt... -- feature.txt...');
  const read = (files) => files.map((file) => readFileSync(file, 'utf8'));
  console.log(pairTable(read(args.slice(0, cut)), read(args.slice(cut + 1))));
}
