// A spectral reference for the molecular atmosphere over a black surface in
// the radiation's visible band, and the grey or few-band depths that follow it:
//   node scripts/rayleighReference.mjs
// The band: of a 5778 K Planck spectrum, the share VISIBLE (the radiation's
// VISIBLE_FRACTION) at the short end, less the share OZONE (0.03) that ozone
// takes, removed from the shortest wavelengths (the Hartley-Huggins bands);
// WAVELENGTHS (40, or 400 with NEAR_INFRARED, where 40 put the reference
// 2.5 % low) equal intervals in wavelength, each weighted by its Planck
// energy, with the Rayleigh depth of Hansen & Travis (1974) at 1013.25 hPa
// at its midpoint,
// tau(l) = 0.008569 l^-4 (1 + 0.0113 l^-2 + 0.00013 l^-4), l in um.
// Each wavelength is reflected by the model's two-stream tau/(tau + 2 mu); the
// band-mean reflectance against mu is the reference the model's depths are
// fitted to. With NEAR_INFRARED=1 the band is instead the rest of the
// spectrum, from VISIBLE's edge to the 0.995 share, carrying 1 - VISIBLE of
// the beam. Beside it, the same spectrum by doubling-adding with the
// azimuth-averaged Rayleigh phase function (scalar, unpolarised; GAUSS nodes),
// which shows the two-stream's own error. The global-mean reflected flux of a
// full atmosphere is (S0/2) x band share x the integral of mu R(mu) over mu,
// since the sunlit hemisphere's area is uniform in mu at every instant.
import { VISIBLE_FRACTION } from '../js/physics/radiation.module.js';

const VISIBLE = Number(process.env.VISIBLE ?? VISIBLE_FRACTION), OZONE = Number(process.env.OZONE ?? 0.03), N = Number(process.env.WAVELENGTHS ?? (process.env.NEAR_INFRARED === '1' ? 400 : 40));
const S0 = 1362, T_SUN = 5778, GAUSS = Number(process.env.GAUSS ?? 24), GREY = Number(process.env.GREY ?? 0.18);
const HC_K = 14387.77;

const planck = (l) => 1 / (l ** 5 * Math.expm1(HC_K / (l * T_SUN)));
const grid = (a, b, n) => Array.from({ length: n + 1 }, (_, k) => a * (b / a) ** (k / n));
function cumulative() {
  const edges = grid(0.05, 200, 400000), cum = [0];
  for (let k = 1; k < edges.length; k++) {
    const m = 0.5 * (edges[k] + edges[k - 1]);
    cum.push(cum[k - 1] + planck(m) * (edges[k] - edges[k - 1]));
  }
  const total = cum[cum.length - 1];
  return { edges, share: cum.map((c) => c / total), total };
}
const spectrum = cumulative();
function wavelengthAt(share) {
  const { edges, share: s } = spectrum;
  let lo = 0, hi = s.length - 1;
  while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (s[mid] < share) lo = mid; else hi = mid; }
  return edges[lo] + (edges[hi] - edges[lo]) * (share - s[lo]) / (s[hi] - s[lo]);
}
const shareBelow = (l) => {
  const { edges, share: s } = spectrum;
  let lo = 0, hi = edges.length - 1;
  while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (edges[mid] < l) lo = mid; else hi = mid; }
  return s[lo] + (s[hi] - s[lo]) * (l - edges[lo]) / (edges[hi] - edges[lo]);
};
const rayleigh = (l) => 0.008569 * l ** -4 * (1 + 0.0113 * l ** -2 + 0.00013 * l ** -4);

const NEAR_INFRARED = process.env.NEAR_INFRARED === '1';
const low = NEAR_INFRARED ? wavelengthAt(VISIBLE) : wavelengthAt(OZONE), high = NEAR_INFRARED ? wavelengthAt(0.995) : wavelengthAt(VISIBLE), SHARE = NEAR_INFRARED ? 1 - VISIBLE : VISIBLE - OZONE;
const bands = [];
for (let k = 0; k < N; k++) {
  const a = low + (high - low) * k / N, b = low + (high - low) * (k + 1) / N, l = 0.5 * (a + b);
  bands.push({ l, weight: shareBelow(b) - shareBelow(a), tau: rayleigh(l) });
}
const bandShare = bands.reduce((s, x) => s + x.weight, 0);
for (const x of bands) x.weight /= bandShare;

const twoStream = (tau, mu) => tau / (tau + 2 * mu);
const mixture = (set, mu) => set.reduce((s, { weight, tau }) => s + weight * twoStream(tau, mu), 0);

function gaussLegendre(n) {
  const x = [], w = [];
  for (let i = 1; i <= n; i++) {
    let z = Math.cos(Math.PI * (i - 0.25) / (n + 0.5)), dp = 0;
    for (let it = 0; it < 100; it++) {
      let p1 = 1, p2 = 0;
      for (let j = 1; j <= n; j++) { const p3 = p2; p2 = p1; p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j; }
      dp = n * (z * p1 - p2) / (z * z - 1);
      const dz = p1 / dp; z -= dz;
      if (Math.abs(dz) < 1e-15) break;
    }
    x.push(0.5 * (1 - z)); w.push(1 / ((1 - z * z) * dp * dp));
  }
  return { x, w };
}
const { x: nodes, w: weights } = gaussLegendre(GAUSS);
const phase = (a, b) => 0.375 * (3 - a * a - b * b + 3 * a * a * b * b);
const matmul = (A, B) => A.map((row) => B[0].map((_, j) => row.reduce((s, v, k) => s + v * B[k][j], 0)));
const matvec = (A, v) => A.map((row) => row.reduce((s, x, k) => s + x * v[k], 0));
function inverse(A) {
  const n = A.length, M = A.map((row, i) => [...row, ...Array.from({ length: n }, (_, j) => (i === j ? 1 : 0))]);
  for (let c = 0; c < n; c++) {
    let p = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
    [M[c], M[p]] = [M[p], M[c]];
    const d = M[c][c];
    for (let j = 0; j < 2 * n; j++) M[c][j] /= d;
    for (let r = 0; r < n; r++) if (r !== c) { const f = M[r][c]; if (f) for (let j = 0; j < 2 * n; j++) M[r][j] -= f * M[c][j]; }
  }
  return M.map((row) => row.slice(n));
}
function doublingAdding(tau, mus) {
  const doublings = 30, d = tau / 2 ** doublings;
  let R = nodes.map((mi) => nodes.map((mj, j) => d / mi * 0.5 * phase(mi, mj) * weights[j]));
  let T = nodes.map((mi, i) => nodes.map((mj, j) => (i === j ? 1 - d / mi : 0) + d / mi * 0.5 * phase(mi, mj) * weights[j]));
  let rho = mus.map((m0) => nodes.map((mi) => d / mi * 0.5 * phase(mi, m0)));
  let sigma = rho.map((v) => [...v]);
  let e = mus.map((m0) => Math.exp(-d / m0));
  for (let s = 0; s < doublings; s++) {
    const M = inverse(matmul(R, R).map((row, i) => row.map((v, j) => (i === j ? 1 : 0) - v)));
    const TM = matmul(T, M);
    for (let b = 0; b < mus.length; b++) {
      const Rr = matvec(R, rho[b]), D = matvec(M, sigma[b].map((v, i) => v + e[b] * Rr[i]));
      const RD = matvec(R, D), U = rho[b].map((v, i) => e[b] * v + RD[i]);
      const TU = matvec(T, U), TD = matvec(T, D);
      rho[b] = rho[b].map((v, i) => v + TU[i]); sigma[b] = sigma[b].map((v, i) => e[b] * v + TD[i]); e[b] *= e[b];
    }
    const RT = matmul(matmul(TM, R), T);
    R = R.map((row, i) => row.map((v, j) => v + RT[i][j]));
    T = matmul(TM, T);
  }
  return mus.map((m0, b) => rho[b].reduce((s, v, i) => s + weights[i] * nodes[i] * v, 0) / m0);
}

const MUS = [0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];
const FLUX_MUS = Array.from({ length: 200 }, (_, k) => (k + 0.5) / 200);
const reference = (mu) => mixture(bands, mu);
const exactByBand = bands.map((x) => doublingAdding(x.tau, [...MUS, ...FLUX_MUS]));
const exact = (index) => bands.reduce((s, x, b) => s + x.weight * exactByBand[b][index], 0);
const flux = (R) => S0 / 2 * SHARE * FLUX_MUS.reduce((s, mu) => s + mu * R(mu), 0) / FLUX_MUS.length;
const exactFlux = S0 / 2 * SHARE * FLUX_MUS.reduce((s, mu, k) => s + mu * exact(MUS.length + k), 0) / FLUX_MUS.length;

const FIT_MUS = Array.from({ length: 91 }, (_, k) => 0.1 + 0.01 * k);
const worstRelative = (R) => Math.max(...FIT_MUS.map((mu) => Math.abs(R(mu) / reference(mu) - 1)));
function nelderMead(f, x0, step, iterations = 4000) {
  const n = x0.length;
  let simplex = [x0, ...x0.map((_, i) => x0.map((v, j) => v + (i === j ? step[j] : 0)))].map((x) => ({ x, f: f(x) }));
  for (let it = 0; it < iterations; it++) {
    simplex.sort((a, b) => a.f - b.f);
    const centroid = x0.map((_, j) => simplex.slice(0, n).reduce((s, p) => s + p.x[j], 0) / n);
    const worst = simplex[n], along = (t) => centroid.map((c, j) => c + t * (worst.x[j] - c));
    const r = along(-1), fr = f(r);
    if (fr < simplex[0].f) { const e = along(-2), fe = f(e); simplex[n] = fe < fr ? { x: e, f: fe } : { x: r, f: fr }; }
    else if (fr < simplex[n - 1].f) simplex[n] = { x: r, f: fr };
    else {
      const c = along(0.5), fc = f(c);
      if (fc < worst.f) simplex[n] = { x: c, f: fc };
      else simplex = simplex.map((p, i) => (i === 0 ? p : { x: p.x.map((v, j) => simplex[0].x[j] + 0.5 * (v - simplex[0].x[j])), f: 0 })).map((p, i) => (i === 0 ? p : { x: p.x, f: f(p.x) }));
    }
  }
  simplex.sort((a, b) => a.f - b.f);
  return simplex[0];
}
const softmax = (z) => { const e = z.map(Math.exp), s = e.reduce((a, b) => a + b, 0); return e.map((v) => v / s); };
function subBands(n) {
  const decode = (x) => { const w = softmax([0, ...x.slice(0, n - 1)]); return w.map((weight, k) => ({ weight, tau: Math.exp(x[n - 1 + k]) })); };
  const objective = (x) => worstRelative((mu) => mixture(decode(x), mu));
  const sorted = [...bands].sort((a, b) => a.tau - b.tau);
  const start = [];
  let acc = 0, k = 0, groups = Array.from({ length: n }, () => ({ w: 0, t: 0 }));
  for (const x of sorted) { const g = Math.min(n - 1, Math.floor(acc * n)); groups[g].w += x.weight; groups[g].t += x.weight * x.tau; acc += x.weight; k++; }
  groups = groups.map((g) => ({ w: g.w, t: g.t / g.w }));
  for (let j = 1; j < n; j++) start.push(Math.log(groups[j].w / groups[0].w));
  for (let j = 0; j < n; j++) start.push(Math.log(groups[j].t));
  let best = nelderMead(objective, start, start.map(() => 0.3));
  for (let r = 0; r < 6; r++) best = nelderMead(objective, best.x, best.x.map(() => 0.05));
  return decode(best.x).sort((a, b) => a.tau - b.tau);
}

const meanTau = bands.reduce((s, x) => s + x.weight * x.tau, 0);
const greyFit = nelderMead((x) => worstRelative((mu) => twoStream(Math.exp(x[0]), mu)), [Math.log(0.15)], [0.2]);
const greyMinimax = Math.exp(greyFit.x[0]);
const referenceFlux = flux(reference);
let lo = 0.01, hi = 1;
for (let it = 0; it < 100; it++) { const m = 0.5 * (lo + hi); if (flux((mu) => twoStream(m, mu)) < referenceFlux) lo = m; else hi = m; }
const greyFlux = 0.5 * (lo + hi);
const fits = [2, 3].map((n) => ({ n, set: subBands(n) }));

const f3 = (x) => x.toFixed(3), f4 = (x) => x.toFixed(4), pct = (x) => `${(100 * x >= 0 ? '+' : '')}${(100 * x).toFixed(1)}%`;
console.log(NEAR_INFRARED ? `near-infrared band ${f3(low)}-${f3(high)} um of a ${T_SUN} K Planck spectrum (share ${f3(bandShare)} of the beam, above the ${VISIBLE} below ${f3(low)} um), ${N} wavelengths` : `visible band ${f3(low)}-${f3(high)} um of a ${T_SUN} K Planck spectrum (share ${f3(bandShare)} of the beam: ${VISIBLE} below ${f3(high)} um less ozone's ${OZONE} below ${f3(low)} um), ${N} wavelengths`);
console.log(`Rayleigh depth at 1013.25 hPa: ${f3(rayleigh(high))} at ${f3(high)} um, ${f4(rayleigh(0.55))} at 0.55 um, ${f3(rayleigh(low))} at ${f3(low)} um; band mean (energy-weighted) ${f4(meanTau)}`);
console.log(`grey depth minimising the largest relative error over mu 0.1-1: ${f4(greyMinimax)} (largest error ${pct(greyFit.f)}); grey depth giving the reference's global-mean flux: ${f4(greyFlux)}`);
for (const { n, set } of fits) console.log(`${n} sub-bands (weight, depth): ${set.map((x) => `(${f4(x.weight)}, ${f4(x.tau)})`).join(' ')}; largest relative error over mu 0.1-1 ${pct(worstRelative((mu) => mixture(set, mu)))}`);
console.log('');
console.log(`   mu  reference  exact(DA)  2s/exact  grey ${GREY}  vs ref   grey ${f3(greyFlux)}  vs ref  ${fits.map(({ n }) => `${n}-band  vs ref`).join('  ')}`);
MUS.forEach((mu, m) => {
  const r = reference(mu), ex = exact(m), g = twoStream(GREY, mu), gf = twoStream(greyFlux, mu);
  console.log(`${mu.toFixed(2).padStart(5)}  ${f4(r).padStart(9)}  ${f4(ex).padStart(9)}  ${pct(r / ex - 1).padStart(8)}  ${f4(g).padStart(9)}  ${pct(g / r - 1).padStart(7)}  ${f4(gf).padStart(10)}  ${pct(gf / r - 1).padStart(6)}  ${fits.map(({ set }) => { const v = mixture(set, mu); return `${f4(v)}  ${pct(v / r - 1).padStart(6)}`; }).join('  ')}`);
});
console.log('');
console.log(`global-mean reflected flux of a full atmosphere over a black surface, W/m2 (S0 ${S0}, band share ${SHARE.toFixed(4)}):`);
console.log(`  reference (two-stream per wavelength) ${referenceFlux.toFixed(2)}; exact (doubling-adding) ${exactFlux.toFixed(2)}`);
console.log(`  grey ${GREY}: ${flux((mu) => twoStream(GREY, mu)).toFixed(2)}; grey ${f4(greyMinimax)} (minimax): ${flux((mu) => twoStream(greyMinimax, mu)).toFixed(2)}; ${fits.map(({ n, set }) => `${n} sub-bands: ${flux((mu) => mixture(set, mu)).toFixed(2)}`).join('; ')}`);
