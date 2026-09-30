// The equatorial ocean's response to the wind in one saved state, on the CPU:
//   [LAT=2] [OCEAN='{"closureHours":48}'] node scripts/equatorialOcean.mjs <state.bin>
// OCEAN adds to everySteps 8 and must match the run's ocean options.
// By 20-degree longitude bin along the equator: the lowest-layer wind and the
// stress it makes, the mixed-layer depth, the depths of the thermocline classes
// and of the 20 C isotherm, the zonal current by class and by depth, the
// momentum budget of the mixed layer and of the water above the 1024 class,
// term by term as the ocean's tendency computes it, the occupancy of the
// thermocline classes with the closure's pull near their token edges, and the
// change of the mixed layer's current over one full ocean step.
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { createModel, STATE_NAMES } from '../js/model.module.js';
import { decodeState, savedLevels } from '../js/stateFile.module.js';
import { cellVector, gradient, curl, kineticEnergy, laplacianVelocity } from '../js/dynamics/operators.module.js';
import { seawaterDensity } from '../js/ocean/seawater.module.js';
import { THIN, PV_FLOOR, THERMOCLINE_DENSITY, closureCoefficient, closureVelocity } from '../js/ocean/layered.module.js';

const FILE = process.argv[2];
if (!FILE) throw new Error('usage: node scripts/equatorialOcean.mjs <state.bin>');
const LAT = Number(process.env.LAT ?? 2);
const OCEAN = { everySteps: 8, ...JSON.parse(process.env.OCEAN ?? '{}') };
const DEFAULTS = { interfacialDrag: 2e-4, bottomDrag: 3e-3, minimumThickness: 50, vorticityCentring: 0.5, closureHours: 12, closureFill: false, density: 1025, gravity: 9.81 };
const opt = { ...DEFAULTS, ...OCEAN };
const say = (s = '') => console.log(s);
const f = (x, d = 1) => (Number.isFinite(x) ? x.toFixed(d) : '—');
const col = (s, n = 8) => String(s).padStart(n);

const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
const N = saved.N, levels = savedLevels(saved);
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const model = createModel(new Grid(N), { topography, levels, ocean: OCEAN });
const { mesh, core, state, surface, seaIce, ocean } = model;
STATE_NAMES.forEach((name, a) => state[a].set(saved[name]));
seaIce.load(state[6], saved.concentration ?? null);
model.time = saved.time;
ocean.load(saved.ocean, state[3], state[6]);

const { K, C, E, sigmaMid, R, g: gAir, exnerLayer } = core.diagnostics;
const { nEdgesOnEdge, edgesOnEdge, weightsOnEdge, maxEdgesOnEdge, dcEdge, dvEdge, cellsOnEdge, verticesOnEdge, fVertex, latCell, lonCell } = mesh;
const L = ocean.layers, rho = ocean.densities, rho0 = opt.density, g = opt.gravity;
const { h, u, Q, eta, cellOcean, edgeOcean } = ocean;
const hEdge = new Float64Array(ocean.shared.hEdge), pressure = new Float64Array(ocean.shared.pressure);
const gradEta = new Float64Array(ocean.shared.gradEta), gradRho = new Float64Array(ocean.shared.gradRho);
const at = (k, i) => k * C + i, ae = (k, e) => k * E + e;
const deg = 180 / Math.PI, dt = 1350 * 16 / N, dtOcean = OCEAN.everySteps * dt;
const upperLayers = rho.filter((r, k) => k > 0 && r < THERMOCLINE_DENSITY).length;
let spacing = 0;
for (let e = 0; e < E; e++) spacing += dcEdge[e];
spacing /= E;
const nu4 = closureCoefficient(spacing, opt.closureHours);

const BINS = [];
for (let west = 120; west < 280; west += 20) BINS.push([west, west + 20]);
const binName = ([a, b]) => `${a <= 180 ? `${a}E` : `${360 - a}W`}-${b <= 180 ? `${b}E` : `${360 - b}W`}`;
const binOf = new Int8Array(C).fill(-1);
for (let i = 0; i < C; i++) {
  if (!cellOcean[i] || Math.abs(latCell[i] * deg) > LAT) continue;
  const lon = ((lonCell[i] * deg) % 360 + 360) % 360;
  const b = BINS.findIndex(([a, c]) => lon >= a && lon < c);
  binOf[i] = b;
}
const east = (i, v) => -Math.sin(lonCell[i]) * v[3 * i] + Math.cos(lonCell[i]) * v[3 * i + 1];
const zonal = (edgeField) => { const v = new Float64Array(3 * C); cellVector(mesh, edgeField, v); return Float64Array.from({ length: C }, (_, i) => east(i, v)); };
function binMean(value, weight = null) {
  const sum = new Float64Array(BINS.length), w = new Float64Array(BINS.length);
  for (let i = 0; i < C; i++) {
    const b = binOf[i];
    if (b < 0) continue;
    const x = typeof value === 'function' ? value(i) : value[i];
    const wi = mesh.areaCell[i] * (weight ? (typeof weight === 'function' ? weight(i) : weight[i]) : 1);
    if (!Number.isFinite(x) || !(wi > 0)) continue;
    sum[b] += wi * x; w[b] += wi;
  }
  return Array.from(sum, (s, b) => (w[b] > 0 ? s / w[b] : NaN));
}
const row = (label, values, d = 1, scale = 1) => say(`${label.padEnd(34)}${values.map((x) => col(f(x * scale, d))).join('')}`);
const header = () => say(`${''.padEnd(34)}${BINS.map((b) => col(binName(b), 8)).join('')}`);
const count = binMean(() => 1).map((_, b) => binOf.reduce((n, x) => n + (x === b ? 1 : 0), 0));

say(`== ${FILE.split('/').pop()}: day ${saved.day}, N=${N}, ${L} ocean layers, ${LAT}S-${LAT}N, ocean step ${dtOcean} s (every ${OCEAN.everySteps} of ${dt} s)`);
say(`   ocean options ${JSON.stringify(OCEAN)}; nu4 ${nu4.toExponential(2)} m^4/s at mean spacing ${(spacing / 1e3).toFixed(1)} km`);
header();
row('sea cells', count, 0);

// ---------------- 1. wind and stress ----------------
core.phaseFlux(state, 0, K);
core.phaseColumn(state, [new Float64Array(C)], 0, C);
surface.lowestWindSpeed(state[2]);
const atmosphereStress = surface.stress(state, new Float64Array(E));
ocean.setStress(atmosphereStress, state[6], seaIce.concentration);
const bottom = K - 1;
const wind = zonal(state[2].subarray(bottom * E, K * E));
const windHeight = Float64Array.from({ length: C }, (_, i) => {
  const T = state[1][bottom * C + i] * exnerLayer[bottom * C + i];
  return R * T / gAir * Math.log(1 / sigmaMid[bottom]);
});
const tauAtm = zonal(atmosphereStress), tauOcean = zonal(ocean.stress);
say('\n-- 1. wind and stress');
row('lowest-layer height (m)', binMean(windHeight), 0);
row('lowest-layer zonal wind (m/s)', binMean(wind), 2);
row('10 m wind, log law z0=1e-4 (m/s)', binMean((i) => wind[i] * Math.log(10 / 1e-4) / Math.log(windHeight[i] / 1e-4)), 2);
row('wind speed |V| (m/s)', binMean(surface.windSpeed), 2);
row('zonal stress, atmosphere (N/m2)', binMean(tauAtm), 4);
row('zonal stress, ocean (N/m2)', binMean(tauOcean), 4);
const wind10 = Float64Array.from({ length: C }, (_, i) => wind[i] * Math.log(10 / 1e-4) / Math.log(windHeight[i] / 1e-4));
const airDensity = Float64Array.from({ length: C }, (_, i) => state[0][i] * sigmaMid[bottom] / (R * state[1][bottom * C + i] * exnerLayer[bottom * C + i]));
const bulk = binMean((i) => airDensity[i] * Math.abs(wind10[i]) * wind10[i]);
row('Cd at 10 m, bin stress / bin rho|u|u (1e-3)', binMean(tauAtm).map((t, b) => t / bulk[b]), 2, 1e3);
row('Large-Pond Cd at that 10 m wind (1e-3)', binMean((i) => (Math.abs(wind10[i]) < 11 ? 1.2 : 0.49 + 0.065 * Math.abs(wind10[i]))), 2);

// ---------------- 2. mixed layer and thermocline ----------------
const layerT = (k, i) => Q[at(k, i)] / Math.max(h[at(k, i)], 1e-9) - 273.15;
function isothermDepth(i, target) {
  let top = 0, prevMid = null, prevT = null;
  for (let k = 0; k < L; k++) {
    const hk = h[at(k, i)];
    if (k > 0 && hk <= THIN) { top += hk; continue; }
    const t = layerT(k, i), mid = k === 0 ? 0 : top + 0.5 * hk;
    if (k === 0 && t < target) return 0;
    if (prevT !== null && t < target) return prevMid + (prevT - target) / (prevT - t) * (mid - prevMid);
    prevMid = k === 0 ? hk : mid; prevT = t; top += hk;
  }
  return NaN;
}
const classTop = (i, density) => { let z = 0; for (let k = 0; k < L; k++) { if (k > 0 && rho[k] >= density) return z; z += h[at(k, i)]; } return z; };
say('\n-- 2. mixed layer and thermocline');
row('SST (C)', binMean((i) => layerT(0, i)), 2);
row('mixed-layer depth h0 (m)', binMean((i) => h[i]), 1);
row(`share of cells with h0 <= ${opt.minimumThickness + 0.5} m (%)`, binMean((i) => (h[i] <= opt.minimumThickness + 0.5 ? 100 : 0)), 0);
row('min h0 in bin (m)', BINS.map((_, b) => { let m = Infinity; for (let i = 0; i < C; i++) if (binOf[i] === b) m = Math.min(m, h[i]); return m; }), 1);
row('max h0 in bin (m)', BINS.map((_, b) => { let m = -Infinity; for (let i = 0; i < C; i++) if (binOf[i] === b) m = Math.max(m, h[i]); return m; }), 1);
row('20 C isotherm depth (m)', binMean((i) => isothermDepth(i, 20)), 1);
row('24 C isotherm depth (m)', binMean((i) => isothermDepth(i, 24)), 1);
for (const d of [1023, 1024, 1025, 1026]) row(`top of the ${d} class (m)`, binMean((i) => classTop(i, d)), 1);
row('sea level eta (cm)', binMean(eta), 1, 100);
row('mixed-layer density - 1000', binMean((i) => seawaterDensity(Q[i] / h[i], ocean.W[i] / h[i]) - 1000), 2);
row('first class under ML (density-1000)', binMean((i) => { for (let k = 1; k < L; k++) if (h[at(k, i)] > THIN) return rho[k] - 1000; return NaN; }), 2);
row('its thickness (m)', binMean((i) => { for (let k = 1; k < L; k++) if (h[at(k, i)] > THIN) return h[at(k, i)]; return NaN; }), 1);

// ---------------- 3. zonal current by class and by depth ----------------
const uZonal = [];
for (let k = 0; k < L; k++) uZonal.push(zonal(u.subarray(k * E, (k + 1) * E)));
say('\n-- 3. zonal current (cm/s), thickness-weighted over the cells where the class is thicker than THIN; [mean depth of its middle, m]');
header();
for (let k = 0; k < L; k++) {
  const w = (i) => (k === 0 || h[at(k, i)] > THIN ? h[at(k, i)] : 0);
  const cover = binMean((i) => (k === 0 || h[at(k, i)] > THIN ? 1 : 0));
  if (Math.max(...cover) < 0.25) continue;
  const mid = binMean((i) => { let z = 0; for (let j = 0; j < k; j++) z += h[at(j, i)]; return z + 0.5 * h[at(k, i)]; }, w);
  const um = binMean(uZonal[k], w);
  say(`${(k === 0 ? 'mixed layer' : `class ${rho[k].toFixed(2)}`).padEnd(18)}${''.padEnd(16)}${BINS.map((_, b) => col(cover[b] >= 0.25 ? `${f(um[b] * 100, 0)}[${f(mid[b], 0)}]` : '', 8)).join('')}`);
}
function currentAtDepth(i, z) {
  let top = 0;
  for (let k = 0; k < L; k++) { const hk = h[at(k, i)]; if (k > 0 && hk <= THIN) { top += hk; continue; } if (z < top + hk) return uZonal[k][i]; top += hk; }
  return NaN;
}
say('\nzonal current by depth (cm/s)');
header();
for (const z of [5, 25, 50, 75, 100, 125, 150, 175, 200, 250, 300]) row(`  ${z} m`, binMean((i) => currentAtDepth(i, z)), 0, 100);
row('  max eastward in 0-300 m', binMean((i) => { let m = -Infinity; for (let z = 5; z <= 300; z += 5) m = Math.max(m, currentAtDepth(i, z)); return m; }), 0, 100);
row('  its depth (m)', binMean((i) => { let m = -Infinity, zm = NaN; for (let z = 5; z <= 300; z += 5) { const c = currentAtDepth(i, z); if (c > m) { m = c; zm = z; } } return zm; }), 0);

// ---------------- 4. momentum budget by term ----------------
const stage = ocean.stages[0];
ocean.tendency(ocean.state, stage);
const du = stage[1];
const terms = ['coriolis f', 'rel. vorticity', '-grad K', '-g grad eta', 'baroclinic PGF', 'stress', 'drag above', 'drag below', 'bottom drag', 'nu4 closure', 'thin-layer relax'];
const budget = terms.map(() => new Float64Array(L * E));
const zeta = new Float64Array(C > 0 ? mesh.nVertices : 0), qf = new Float64Array(E), qz = new Float64Array(E), fluxPV = new Float64Array(E);
const Kc = new Float64Array(C), grad = new Float64Array(E), lap = new Float64Array(E), lap2 = new Float64Array(E), divS = new Float64Array(C), curlS = new Float64Array(mesh.nVertices), phi = new Float64Array(C), closureU = new Float64Array(E);
const relax = Math.min(1 / 3600, 1 / dtOcean);
for (let k = 0; k < L; k++) {
  const oe = k * E, uk = u.subarray(oe, oe + E);
  for (let e = 0; e < E; e++) fluxPV[e] = edgeOcean[e] ? hEdge[oe + e] * uk[e] : 0;
  curl(mesh, uk, zeta);
  for (let e = 0; e < E; e++) {
    const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1], va = verticesOnEdge[2 * e], vb = verticesOnEdge[2 * e + 1];
    const denom = Math.max(hEdge[oe + e], k > 0 ? opt.vorticityCentring * 0.5 * (h[at(k, a)] + h[at(k, b)]) : 0, PV_FLOOR);
    qf[e] = 0.5 * (fVertex[va] + fVertex[vb]) / denom;
    qz[e] = 0.5 * (zeta[va] + zeta[vb]) / denom;
  }
  kineticEnergy(mesh, uk, Kc);
  gradient(mesh, Kc, grad);
  for (let i = 0; i < C; i++) phi[i] = k === 0 ? 0 : g * pressure[at(k, i)] / rho0;
  const gradP = new Float64Array(E);
  gradient(mesh, phi, gradP);
  if (nu4 > 0) { laplacianVelocity(mesh, k === 0 || !opt.closureFill ? uk : closureVelocity(mesh, uk, hEdge.subarray(oe, oe + E), edgeOcean, closureU), lap, divS, curlS); laplacianVelocity(mesh, lap, lap2, divS, curlS); }
  for (let e = 0; e < E; e++) {
    if (!edgeOcean[e]) continue;
    let sf = 0, sz = 0;
    for (let s = 0; s < nEdgesOnEdge[e]; s++) {
      const o = edgesOnEdge[maxEdgesOnEdge * e + s], w = weightsOnEdge[maxEdgesOnEdge * e + s] * dvEdge[o] * fluxPV[o];
      sf += w * 0.5 * (qf[e] + qf[o]); sz += w * 0.5 * (qz[e] + qz[o]);
    }
    const n = oe + e, he = Math.max(hEdge[n], opt.minimumThickness), hd = k === 0 ? he : Math.max(hEdge[n], THIN);
    budget[0][n] = sf / dcEdge[e];
    budget[1][n] = sz / dcEdge[e];
    budget[2][n] = -grad[e];
    budget[3][n] = -g * gradEta[e];
    budget[4][n] = k === 0 ? -g / rho0 * 0.5 * hEdge[e] * gradRho[e] : -gradP[e];
    if (k === 0) budget[5][n] = ocean.stress[e] / rho0 / he;
    if (k > 0) { let j = k - 1; while (j > 0 && hEdge[ae(j, e)] < THIN) j--; budget[6][n] = opt.interfacialDrag * (u[ae(j, e)] - u[n]) / hd; }
    if (k < L - 1) { let j = k + 1; while (j < L - 1 && hEdge[ae(j, e)] < THIN) j++; if (hEdge[ae(j, e)] >= THIN) budget[7][n] = -opt.interfacialDrag * (u[n] - u[ae(j, e)]) / hd; }
    let isBottom = k === L - 1;
    if (!isBottom) { isBottom = true; for (let j = k + 1; j < L; j++) if (hEdge[ae(j, e)] >= THIN) { isBottom = false; break; } }
    if (isBottom) budget[8][n] = -opt.bottomDrag * Math.abs(u[n]) * u[n] / he;
    budget[9][n] = nu4 > 0 ? -nu4 * lap2[e] : 0;
    if (k > 0 && hEdge[n] < THIN) {
      for (let t = 0; t < 10; t++) budget[t][n] = 0;
      budget[10][n] = (u[ae(k - 1, e)] - u[n]) * relax;
    }
  }
}
let worst = 0, scale = 0;
for (let n = 0; n < L * E; n++) { let s = 0; for (const t of budget) s += t[n]; worst = Math.max(worst, Math.abs(s - du[n])); scale = Math.max(scale, Math.abs(du[n])); }
say(`\n-- 4. momentum budget (1e-7 m/s2, zonal, cell-reconstructed); terms sum to the ocean's own tendency within ${worst.toExponential(2)} of max ${scale.toExponential(2)} m/s2`);
header();
say('mixed layer (k = 0), per unit mass:');
for (let t = 0; t < terms.length; t++) if (t !== 6 && t !== 8 && t !== 10) row(`  ${terms[t]}`, binMean(zonal(budget[t].subarray(0, E))), 2, 1e7);
row('  sum = du0/dt', binMean(zonal(du.subarray(0, E))), 2, 1e7);
row('  -g d(eta)/dx alone as eta slope (cm/1000km)', binMean(zonal(budget[3].subarray(0, E))), 1, -1e8 / g);

const integrated = terms.map(() => new Float64Array(E)), integratedDu = new Float64Array(E), weightE = new Float64Array(E);
for (let e = 0; e < E; e++) {
  for (let k = 0; k <= upperLayers; k++) {
    const n = ae(k, e);
    for (let t = 0; t < terms.length; t++) integrated[t][e] += hEdge[n] * budget[t][n];
    integratedDu[e] += hEdge[n] * du[n];
    weightE[e] += hEdge[n];
  }
}
say(`\nlayers above the ${THERMOCLINE_DENSITY} class (k = 0..${upperLayers}), depth-integrated h*term (1e-5 m2/s2, h-weighted; the stress row is tau/rho0 where h0 >= minimumThickness):`);
for (let t = 0; t < terms.length; t++) if (t !== 10 && t !== 8) row(`  ${terms[t]}`, binMean(zonal(integrated[t])), 2, 1e5);
row('  sum', binMean(zonal(Float64Array.from({ length: E }, (_, e) => integrated.reduce((s, x) => s + x[e], 0)))), 2, 1e5);
const pgf = binMean(zonal(Float64Array.from({ length: E }, (_, e) => integrated[3][e] + integrated[4][e])));
const tx = binMean(zonal(integrated[5]));
row('  pressure force / stress (-1 = balance)', pgf.map((p, b) => p / tx[b]), 2);

say('\nper class (1e-7 m/s2 per unit mass, thickness-weighted where thicker than THIN): u cm/s | PGF = -g grad eta + baroclinic | drag above + below | nu4 | Coriolis + rel. vorticity | sum');
const termZonal = budget.map((t) => { const out = []; for (let k = 0; k < L; k++) out.push(zonal(t.subarray(k * E, (k + 1) * E))); return out; });
const duZonal = []; for (let k = 0; k < L; k++) duZonal.push(zonal(du.subarray(k * E, (k + 1) * E)));
for (const b of [2, 3, 4, 5, 6]) {
  say(`  ${binName(BINS[b])}`);
  for (let k = 0; k < L; k++) {
    const w = (i) => (binOf[i] === b && (k === 0 || h[at(k, i)] > THIN) ? h[at(k, i)] : 0);
    let cover = 0, n = 0; for (let i = 0; i < C; i++) if (binOf[i] === b) { n++; if (k === 0 || h[at(k, i)] > THIN) cover++; }
    if (cover < 0.25 * n) continue;
    const m = (v) => binMean(v, w)[b];
    const mid = m((i) => { let z = 0; for (let j = 0; j < k; j++) z += h[at(j, i)]; return z + 0.5 * h[at(k, i)]; });
    if (mid > 400) break;
    const pg = m((i) => termZonal[3][k][i] + termZonal[4][k][i]), dr = m((i) => termZonal[6][k][i] + termZonal[7][k][i]), st = m(termZonal[5][k]);
    say(`    ${(k === 0 ? 'ML' : rho[k].toFixed(2)).padEnd(8)} z ${f(mid, 0).padStart(4)} h ${f(m((i) => h[at(k, i)]), 0).padStart(3)}  u ${f(m(uZonal[k]) * 100, 1).padStart(6)} | stress ${f(st * 1e7, 2).padStart(6)} PGF ${f(pg * 1e7, 2).padStart(6)} (eta ${f(m(termZonal[3][k]) * 1e7, 2)}) drag ${f(dr * 1e7, 2).padStart(6)} nu4 ${f(m(termZonal[9][k]) * 1e7, 2).padStart(6)} cor ${f(m((i) => termZonal[0][k][i] + termZonal[1][k][i] + termZonal[2][k][i]) * 1e7, 2).padStart(6)} relax ${f(m(termZonal[10][k]) * 1e7, 2).padStart(5)} | sum ${f(m(duZonal[k]) * 1e7, 2).padStart(6)}`);
  }
}

const columns = [];
const column = (t) => { if (columns[t]) return columns[t]; const out = new Float64Array(E); for (let e = 0; e < E; e++) for (let k = 0; k < L; k++) out[e] += hEdge[ae(k, e)] * budget[t][ae(k, e)]; return (columns[t] = out); };
say('\nwhole column, depth-integrated h*term (1e-5 m2/s2): interfacial drag sums to zero where the mixed layer is at least minimumThickness thick at the edge');
header();
row('  stress', binMean(zonal(column(5))), 2, 1e5);
row('  -g grad eta + baroclinic PGF', binMean(zonal(Float64Array.from(column(3), (x, e) => x + column(4)[e]))), 2, 1e5);
row('  interfacial drag (net)', binMean(zonal(Float64Array.from(column(6), (x, e) => x + column(7)[e]))), 2, 1e5);
row('  bottom drag', binMean(zonal(column(8))), 2, 1e5);
row('  nu4 closure', binMean(zonal(column(9))), 2, 1e5);
row('  Coriolis + rel. vorticity + grad K', binMean(zonal(Float64Array.from(column(0), (x, e) => x + column(1)[e] + column(2)[e]))), 2, 1e5);

say('\nmeridional section of the zonal current (cm/s) over 140W-110W, by latitude');
const LATS = []; for (let a = -8; a < 8; a += 1) LATS.push(a);
say(`${'depth'.padEnd(10)}${LATS.map((a) => col(`${a}..${a + 1}`, 7)).join('')}`);
const latBin = Int8Array.from({ length: C }, (_, i) => { const lon = ((lonCell[i] * deg) % 360 + 360) % 360, la = latCell[i] * deg; return cellOcean[i] && lon >= 220 && lon < 250 && la >= -8 && la < 8 ? Math.floor(la + 8) : -1; });
const latMean = (fn) => { const s = new Float64Array(LATS.length), w = new Float64Array(LATS.length); for (let i = 0; i < C; i++) { const b = latBin[i]; if (b < 0) continue; const x = fn(i); if (!Number.isFinite(x)) continue; s[b] += mesh.areaCell[i] * x; w[b] += mesh.areaCell[i]; } return Array.from(s, (x, b) => x / w[b]); };
for (const z of [5, 50, 75, 100, 125, 150, 200]) say(`${`${z} m`.padEnd(10)}${latMean((i) => currentAtDepth(i, z)).map((x) => col(f(x * 100, 0), 7)).join('')}`);
say(`${'h0 (m)'.padEnd(10)}${latMean((i) => h[i]).map((x) => col(f(x, 0), 7)).join('')}`);
say(`${'20C (m)'.padEnd(10)}${latMean((i) => isothermDepth(i, 20)).map((x) => col(f(x, 0), 7)).join('')}`);
say(`${'SST (C)'.padEnd(10)}${latMean((i) => layerT(0, i)).map((x) => col(f(x, 1), 7)).join('')}`);
say(`${'taux'.padEnd(10)}${latMean((i) => tauOcean[i]).map((x) => col(f(x * 1000, 0), 7)).join('')}  (1e-3 N/m2)`);
const tauNorth = (() => { const v = new Float64Array(3 * C); cellVector(mesh, ocean.stress, v); return Float64Array.from({ length: C }, (_, i) => -Math.sin(latCell[i]) * (Math.cos(lonCell[i]) * v[3 * i] + Math.sin(lonCell[i]) * v[3 * i + 1]) + Math.cos(latCell[i]) * v[3 * i + 2]); })();
say(`${'tauy'.padEnd(10)}${latMean((i) => tauNorth[i]).map((x) => col(f(x * 1000, 0), 7)).join('')}  (1e-3 N/m2)`);

// ---------------- 5. class occupancy and the closure near token edges ----------------
const { edgesOnCell, nEdgesOnCell, maxEdges, edgesOnVertex, nEdge, latEdge } = mesh;
const eastOfEdge = new Float64Array(E), bandEdge = new Uint8Array(E), bandCell = (i) => cellOcean[i] && Math.abs(latCell[i] * deg) <= LAT && ((lonCell[i] * deg) % 360 + 360) % 360 >= 180 && ((lonCell[i] * deg) % 360 + 360) % 360 < 250;
for (let e = 0; e < E; e++) {
  const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
  const lon = Math.atan2(Math.sin(lonCell[a]) + Math.sin(lonCell[b]), Math.cos(lonCell[a]) + Math.cos(lonCell[b]));
  eastOfEdge[e] = -Math.sin(lon) * nEdge[3 * e] + Math.cos(lon) * nEdge[3 * e + 1];
  bandEdge[e] = edgeOcean[e] && bandCell(a) && bandCell(b) && Math.abs(latEdge[e] * deg) <= LAT ? 1 : 0;
}
const ring = (mark) => {
  const out = new Uint8Array(E);
  for (let e = 0; e < E; e++) {
    if (!mark[e]) continue;
    for (const c of [cellsOnEdge[2 * e], cellsOnEdge[2 * e + 1]]) for (let m = 0; m < nEdgesOnCell[c]; m++) out[edgesOnCell[maxEdges * c + m]] = 1;
    for (const v of [verticesOnEdge[2 * e], verticesOnEdge[2 * e + 1]]) for (let m = 0; m < 3; m++) { const o = edgesOnVertex[3 * v + m]; if (o >= 0) out[o] = 1; }
  }
  return out;
};
say(`\n-- 5. thermocline classes on the equator, ${LAT}S-${LAT}N 180-110W: cells holding the class (> THIN), edges with it on both sides, token edges (hEdge < THIN);`);
say('   the closure (1e-7 m/s2, zonal, least squares over edge normals) on edges within two rings of a token edge and on the others, and the zonal current there (cm/s)');
say('class     cells  both  token   nu4 near  nu4 clean   u near  u clean');
for (let k = 1; k < L && rho[k] < 1025.6; k++) {
  let cells = 0, present = 0, edges = 0, both = 0, tokens = 0;
  for (let i = 0; i < C; i++) if (bandCell(i)) { cells++; if (h[at(k, i)] > THIN) present++; }
  const token = Uint8Array.from({ length: E }, (_, e) => (hEdge[ae(k, e)] < THIN || !edgeOcean[e] ? 1 : 0));
  const near = ring(ring(token));
  const acc = [[0, 0, 0], [0, 0, 0]];
  for (let e = 0; e < E; e++) {
    if (!bandEdge[e]) continue;
    edges++;
    if (h[at(k, cellsOnEdge[2 * e])] > THIN && h[at(k, cellsOnEdge[2 * e + 1])] > THIN) both++;
    if (token[e]) { tokens++; continue; }
    const a = acc[near[e] ? 0 : 1], w = eastOfEdge[e];
    a[0] += budget[9][ae(k, e)] * w; a[1] += w * w; a[2] += u[ae(k, e)] * w;
  }
  if (present < 0.02 * cells) continue;
  const z = (a, j) => (a[1] > 0 ? a[j] / a[1] : NaN);
  say(`${rho[k].toFixed(2).padEnd(9)}${col(f(100 * present / cells, 0), 5)}%${col(f(100 * both / edges, 0), 5)}%${col(f(100 * tokens / edges, 0), 5)}%${col(f(z(acc[0], 0) * 1e7, 2), 11)}${col(f(z(acc[1], 0) * 1e7, 2), 11)}${col(f(z(acc[0], 2) * 100, 1), 9)}${col(f(z(acc[1], 2) * 100, 1), 9)}`);
}

// ---------------- 6. one full ocean step ----------------
const u0 = Float64Array.from(u.subarray(0, E)), uAll = Float64Array.from(u);
const surfaceT = Float64Array.from(state[3]), iceCopy = Float64Array.from(state[6]), oceanFlux = new Float64Array(C);
const h0Before = Float64Array.from(h.subarray(0, C));
let stepped = false;
for (let s = 0; s < OCEAN.everySteps && !stepped; s++) stepped = ocean.advance(surfaceT, iceCopy, oceanFlux, ocean.stress.slice(), dt, seaIce.concentration);
const tauUsed = zonal(ocean.stress);
say(`\n-- 6. one full ocean step of ${dtOcean} s from this state (stress as sampled at this atmosphere step; the ice factor re-applied: zonal stress ${f(binMean(tauUsed)[6], 4)} N/m2 in 240-260E)`);
header();
row('du0/dt over the step (1e-7 m/s2)', binMean(zonal(Float64Array.from({ length: E }, (_, e) => (u[e] - u0[e]) / dtOcean))), 2, 1e7);
row('tendency at the start (1e-7 m/s2)', binMean(zonal(du.subarray(0, E))), 2, 1e7);
row('depth-mean shift and rules (1e-7)', binMean(zonal(Float64Array.from({ length: E }, (_, e) => (u[e] - u0[e]) / dtOcean - du[e]))), 2, 1e7);
row('h0 change over the step (m)', binMean((i) => h[i] - h0Before[i]), 2);
row('stress alone on u0 over 1 day (m/s)', binMean(zonal(budget[5].subarray(0, E))), 2, 86400);
