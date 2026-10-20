// ERA5's transformed Eulerian mean on its pressure levels from range reads
// of the raw files: per sample (date, hour) zonal means and eddy
// covariances over 1440 longitudes on each 0.25-degree row within LAT of
// the equator, then psi = [v'θ']/θ_p, ω* = ω + (psi cos φ)_φ/(a cos φ)
// (φ-derivative over ±DROWS rows), w* = -H ω*/p (H = 7 km) and the
// adiabatic cooling ω* θ_p T/θ, averaged over the samples.
//   YEARS=2010-2019 DATES=03-23,03-29 HOURS=0,12 LEVELS=1,2,3 OUT=x.json node scripts/era5Tem.mjs
import { writeFileSync } from 'node:fs';
import { fieldRows, pool, ERA5 } from './era5Raw.mjs';

const [y0, y1] = (process.env.YEARS ?? '2010-2019').split('-').map(Number);
const DATES = (process.env.DATES ?? '03-23,03-29,04-04,04-10,04-16').split(',');
const HOURS = (process.env.HOURS ?? '0,12').split(',').map(Number);
const LEVELS = (process.env.LEVELS ?? '1,2,3,5,7,10,20,30,50,70').split(',').map(Number);
const LAT = Number(process.env.LAT ?? 20), DROWS = Number(process.env.DROWS ?? 4), TROPICS = Number(process.env.TROPICS ?? 15);
const VARS = (process.env.VARS ?? 'temperature,u_component_of_wind,v_component_of_wind,vertical_velocity').split(',');
const KAPPA = 287.06 / 1004.7, H = 7000, a = 6371220, NX = ERA5.NX;
const r0 = Math.round((90 - LAT) / 0.25), r1 = Math.round((90 + LAT) / 0.25) + 1, R = r1 - r0, L = LEVELS.length;
const lat = Float64Array.from({ length: R }, (_, r) => ERA5.latOfRow(r0 + r) * Math.PI / 180);
const names = ['u', 'T', 'omega', 'wBar', 'wStar', 'cooling', 'vTheta', 'uv', 'uOmega'];
const sums = Object.fromEntries(names.map((n) => [n, new Float64Array(L * R)]));
let samples = 0;
const t0 = Date.now();
for (let y = y0; y <= y1; y++) for (const md of DATES) {
  const date = `${y}/${md.replace('-', '/')}`;
  for (const hour of HOURS) {
    const tasks = [];
    for (const level of LEVELS) for (const v of VARS) tasks.push(() => fieldRows(date, v, level, hour, r0, r1));
    const got = await pool(tasks, 16), field = (l, v) => got[l * VARS.length + VARS.indexOf(v)];
    const z = Object.fromEntries(['u', 'v', 'omega', 'T', 'vT', 'uv', 'uOmega'].map((n) => [n, new Float64Array(L * R)]));
    for (let l = 0; l < L; l++) {
      const T = field(l, 'temperature'), u = field(l, 'u_component_of_wind'), v = field(l, 'v_component_of_wind'), w = field(l, 'vertical_velocity');
      for (let r = 0; r < R; r++) {
        let su = 0, sv = 0, sw = 0, sT = 0, svT = 0, suv = 0, suw = 0;
        for (let x = r * NX; x < (r + 1) * NX; x++) {
          su += u ? u[x] : 0; sv += v ? v[x] : 0; sw += w ? w[x] : 0; sT += T[x];
          if (v) { svT += v[x] * T[x]; if (u) suv += u[x] * v[x]; }
          if (u && w) suw += u[x] * w[x];
        }
        const n = l * R + r;
        z.u[n] = su / NX; z.v[n] = sv / NX; z.omega[n] = sw / NX; z.T[n] = sT / NX;
        z.vT[n] = svT / NX - z.v[n] * z.T[n]; z.uv[n] = suv / NX - z.u[n] * z.v[n]; z.uOmega[n] = suw / NX - z.u[n] * z.omega[n];
      }
    }
    const theta = (l, r) => z.T[l * R + r] * (1000 / LEVELS[l]) ** KAPPA;
    const thetaP = (l, r) => { const lo = Math.max(0, l - 1), hi = Math.min(L - 1, l + 1); return (theta(hi, r) - theta(lo, r)) / Math.log(LEVELS[hi] / LEVELS[lo]) / (100 * LEVELS[l]); };
    const psi = (l, r) => z.vT[l * R + r] * (1000 / LEVELS[l]) ** KAPPA / thetaP(l, r);
    for (let l = 0; l < L; l++) for (let r = 0; r < R; r++) {
      const n = l * R + r, lo = Math.max(0, r - DROWS), hi = Math.min(R - 1, r + DROWS), p = 100 * LEVELS[l];
      const omegaStar = z.omega[n] + (psi(l, hi) * Math.cos(lat[hi]) - psi(l, lo) * Math.cos(lat[lo])) / (lat[hi] - lat[lo]) / (a * Math.cos(lat[r]));
      sums.u[n] += z.u[n]; sums.T[n] += z.T[n]; sums.omega[n] += z.omega[n];
      sums.wBar[n] += -H * z.omega[n] / p; sums.wStar[n] += -H * omegaStar / p;
      sums.cooling[n] += 86400 * omegaStar * thetaP(l, r) * z.T[n] / theta(l, r);
      sums.vTheta[n] += z.vT[n] * (1000 / LEVELS[l]) ** KAPPA; sums.uv[n] += z.uv[n]; sums.uOmega[n] += z.uOmega[n];
    }
    samples++;
    process.stderr.write(`${date} ${hour}Z done (${samples}, ${((Date.now() - t0) / 1000).toFixed(0)} s)\n`);
  }
}
for (const n of names) sums[n] = sums[n].map((x) => x / samples);
const mean = (x, l, from, to) => { let s = 0, w = 0; for (let r = 0; r < R; r++) { const d = lat[r] * 180 / Math.PI; if (d < from || d > to) continue; s += Math.cos(lat[r]) * x[l * R + r]; w += Math.cos(lat[r]); } return s / w; };
console.log(`ERA5 ${y0}-${y1}, ${DATES.join(' ')}, ${HOURS.join('/')} UTC: ${samples} samples; ${TROPICS}S-${TROPICS}N means: T K; u at 5S-5N m/s; w-bar, w*, eddy part mm/s; cooling K/day`);
console.log('    hPa      T   u_eq    w-bar     w*   eddy  cooling');
LEVELS.forEach((p, l) => {
  const w = mean(sums.wStar, l, -TROPICS, TROPICS), wb = mean(sums.wBar, l, -TROPICS, TROPICS);
  console.log(`  ${String(p).padStart(5)} ${mean(sums.T, l, -TROPICS, TROPICS).toFixed(1).padStart(6)} ${mean(sums.u, l, -5, 5).toFixed(1).padStart(6)} ${(1e3 * wb).toFixed(3).padStart(8)} ${(1e3 * w).toFixed(3).padStart(6)} ${(1e3 * (w - wb)).toFixed(3).padStart(6)} ${mean(sums.cooling, l, -TROPICS, TROPICS).toFixed(2).padStart(8)}`);
});
if (process.env.OUT) writeFileSync(process.env.OUT, JSON.stringify({ years: [y0, y1], dates: DATES, hours: HOURS, levels: LEVELS, samples, lat: Array.from(lat, (x) => x * 180 / Math.PI), ...Object.fromEntries(names.map((n) => [n, Array.from(sums[n])])) }));
