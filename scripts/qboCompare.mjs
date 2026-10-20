// ERA5 against the runs on ERA5's levels: the 5S-5N wind's mean, the sd of its monthly-mean
// anomalies from the monthly composite, the periodogram's strongest period (12-60 months), the
// semiannual harmonic (amplitude, day of the westerly maximum from the March equinox), and the
// annual-mean w* over 15S-15N from the TEM means; each run over days FIRST-LAST of its fixedSst.mjs
// WINDS=1 log and from its TEM_HOURS means.
//   node scripts/qboCompare.mjs data/era5_equatorial_wind.json data/era5_tem_annual.json FIRST LAST NAME=RUN.log[,RUN_tem.json] ...
import { readFileSync } from 'node:fs';
const [saoFile, temFile, first, last, ...runs] = process.argv.slice(2);
const LEVELS = [1, 2, 3, 5, 7, 10, 20, 30, 50], W_LEVELS = [1, 1.5, 2, 3, 5, 7, 10, 20, 30, 50, 70];
const YEAR = 365.2425, monthOf = (t) => Math.floor((((t % YEAR) + YEAR) % YEAR) / YEAR * 12), stamp = (t) => Math.floor(t / YEAR) * 12 + monthOf(t);
function era5Samples(file) {
  const { levels, series } = JSON.parse(readFileSync(file, 'utf8'));
  return series.map((s) => { const [y, m, d] = s.date.split('/').map(Number); return { t: (Date.UTC(y, m - 1, d) - Date.UTC(2007, 2, 20)) / 86400e3, u: LEVELS.map((p) => s.u[levels.indexOf(p)]) }; });
}
function logSamples(file, from, to) {
  const out = [];
  for (const line of readFileSync(file, 'utf8').split('\n')) {
    const m = line.match(/^upper winds day (\d+) \([^)]*\): (.*)$/);
    if (!m || Number(m[1]) < from || Number(m[1]) > to) continue;
    const layers = m[2].split('; ').map((s) => { const [p, , w] = s.split(' '); return [Number(p), Number(w.split('/')[2])]; }).sort((a, b) => a[0] - b[0]);
    out.push({ t: Number(m[1]), u: LEVELS.map((p) => {
      let j = 0;
      while (j < layers.length - 2 && layers[j + 1][0] < p) j++;
      const w = Math.log(p / layers[j][0]) / Math.log(layers[j + 1][0] / layers[j][0]);
      return (1 - w) * layers[j][1] + w * layers[j + 1][1];
    }) });
  }
  return out;
}
function stats(samples) {
  return LEVELS.map((p, l) => {
    const sums = new Float64Array(12), counts = new Float64Array(12);
    for (const s of samples) { sums[monthOf(s.t)] += s.u[l]; counts[monthOf(s.t)]++; }
    const comp = Array.from(sums, (x, m) => x / counts[m]), mean = comp.reduce((a, b) => a + b) / 12;
    let c = 0, sn = 0;
    comp.forEach((x, m) => { const a = 2 * Math.PI * 2 * (m + 0.5) / 12; c += x * Math.cos(a); sn += x * Math.sin(a); });
    const sao = 2 * Math.hypot(c, sn) / 12, phase = ((Math.atan2(sn, c) / (2 * Math.PI) * 6 + 6) % 6) / 12 * YEAR;
    const monthly = new Map();
    for (const s of samples) { const k = stamp(s.t), e = monthly.get(k) ?? [0, 0]; e[0] += s.u[l] - comp[monthOf(s.t)]; e[1]++; monthly.set(k, e); }
    const keys = [...monthly.keys()].sort((a, b) => a - b), anomalies = keys.map((k) => monthly.get(k)[0] / monthly.get(k)[1]);
    const sd = Math.sqrt(anomalies.reduce((a, x) => a + x * x, 0) / anomalies.length);
    let best = [0, 0];
    for (let period = 12; period <= 60; period += 0.25) {
      let cc = 0, ss = 0;
      anomalies.forEach((x, n) => { cc += x * Math.cos(2 * Math.PI * keys[n] / period); ss += x * Math.sin(2 * Math.PI * keys[n] / period); });
      if (cc * cc + ss * ss > best[1]) best = [period, cc * cc + ss * ss];
    }
    return { mean, sd, period: best[0], sao, phase };
  });
}
function era5W(file) {
  const t = JSON.parse(readFileSync(file, 'utf8')), B = t.lat.length;
  return W_LEVELS.map((p) => {
    const l = t.levels.indexOf(p);
    if (l < 0) return NaN;
    let s = 0, w = 0;
    for (let b = 0; b < B; b++) if (Math.abs(t.lat[b]) < 15) { const c = Math.cos(t.lat[b] * Math.PI / 180); s += c * t.wStar[l * B + b]; w += c; }
    return 1000 * s / w;
  });
}
function modelW(file) {
  const t = JSON.parse(readFileSync(file, 'utf8')), B = t.lat.length;
  let belt = 0;
  for (let b = 0; b < B; b++) if (Math.abs(t.lat[b]) < 15) belt += t.area[b];
  return W_LEVELS.map((p) => {
    const j = t.pressures.findIndex((x) => Math.abs(x - 100 * p) < 1e-6 * x);
    if (j < 0) return NaN;
    let m = 0;
    for (let b = 0; b < B; b++) if (Math.abs(t.lat[b]) < 15) m += -t.omegaStar[j * B + b] * t.area[b];
    return 1000 * m * 7000 / (100 * p * belt);
  });
}
const columns = [{ name: 'ERA5 2007-19', s: stats(era5Samples(saoFile)), w: era5W(temFile) }];
for (const run of runs) {
  const [name, files] = run.split('='), [log, tem] = files.split(',');
  const samples = logSamples(log, Number(first), Number(last));
  columns.push({ name: `${name} (${samples.length} d)`, s: stats(samples), w: tem ? modelW(tem) : W_LEVELS.map(() => NaN) });
}
const f = (x, d = 1) => (Number.isFinite(x) ? x.toFixed(d) : '-').padStart(7);
const table = (title, pick, d = 1) => { console.log(`${title}\n   hPa ${columns.map((c) => c.name.padStart(20)).join('')}`); LEVELS.forEach((p, l) => console.log(`${String(p).padStart(6)} ${columns.map((c) => f(pick(c.s[l]), d).padStart(20)).join('')}`)); };
table('5S-5N zonal wind, mean (m/s)', (s) => s.mean);
table('sd of monthly-mean anomalies (m/s)', (s) => s.sd);
table('strongest period of the anomalies, 12-60 months', (s) => s.period);
console.log(`semiannual harmonic, amplitude m/s / day of westerly maximum\n   hPa ${columns.map((c) => c.name.padStart(20)).join('')}`);
LEVELS.slice(0, 5).forEach((p, l) => console.log(`${String(p).padStart(6)} ${columns.map((c) => `${c.s[l].sao.toFixed(1)} / ${c.s[l].phase.toFixed(0)}`.padStart(20)).join('')}`));
console.log(`w* over 15S-15N, annual mean (mm/s)\n   hPa ${columns.map((c) => c.name.padStart(20)).join('')}`);
W_LEVELS.forEach((p, l) => console.log(`${String(p).padStart(6)} ${columns.map((c) => f(c.w[l], 2).padStart(20)).join('')}`));
