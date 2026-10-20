// The transformed Eulerian mean of saved states on pressure surfaces
// (Andrews, Holton and Leovy 1987), in pressure coordinates: each
// column interpolated in ln p to PRESSURES (hPa), zonal means and eddy
// covariances over BAND-degree latitude bands, and per band and pressure
// psi = [v'θ']/θ_p, v* = v - psi_p, ω* = ω + (psi cos φ)_φ / (a cos φ),
// w* = -H ω*/p with H = 7 km, and the residual ascent's adiabatic cooling
// ω* θ_p T/θ (K/day), each the mean over the states. Printed as means over
// the bands within TROPICS degrees of the equator; OUT names a JSON file
// for the band tables.
//   PRESSURES=1,2,3,5,7,10 node scripts/residualCirculation.mjs STATE.bin [STATE.bin ...]
import { readFileSync, writeFileSync } from 'node:fs';
import { cellVector } from '../js/dynamics/operators.module.js';
import { R_DRY, CP_DRY, P0, GRAVITY, createSigmaCore } from '../js/dynamics/sigmaCore.module.js';
import { decodeState, savedLevels } from '../js/stateFile.module.js';
import { stateGrid, stateRelax } from '../js/stateGrid.module.js';

const KAPPA = R_DRY / CP_DRY, DEG = 180 / Math.PI, H = 7000;
export const TEM_PRESSURES = [0.2, 0.3, 0.4, 0.5, 0.7, 1, 1.5, 2, 3, 4, 5, 7, 10, 15, 20, 30, 50, 70, 100];

// Per cell and pressure (Pa): θ, T, the east and north wind and ω (Pa/s),
// interpolated in ln p between the layers' midpoints σ π, ω from the
// interfaces' π σ-dot and σ (π_t + v·∇π) as js/levels.module.js's
// verticalVelocity forms them.
export function pressureColumns(mesh, core, state, pressures) {
  const { K, C, E, sigmaMid, sigmaUpper, sigmaLower, dSigma } = core.diagnostics;
  const { maxEdges, nEdgesOnCell, edgesOnCell, edgeSignOnCell, cellsOnEdge, dvEdge, areaCell, latCell, lonCell } = mesh;
  const [pi, theta, u] = state, J = pressures.length;
  const shape = Float64Array.from({ length: K }, (_, k) => (sigmaLower[k] ** (1 + KAPPA) - sigmaUpper[k] ** (1 + KAPPA)) / ((1 + KAPPA) * dSigma[k]));
  const east = new Float64Array(K * C), north = new Float64Array(K * C), vector = new Float64Array(3 * C);
  for (let k = 0; k < K; k++) {
    cellVector(mesh, u.subarray(k * E, (k + 1) * E), vector);
    for (let i = 0; i < C; i++) {
      const lat = latCell[i], lon = lonCell[i], x = vector[3 * i], y = vector[3 * i + 1], z = vector[3 * i + 2];
      east[k * C + i] = -Math.sin(lon) * x + Math.cos(lon) * y;
      north[k * C + i] = -Math.sin(lat) * Math.cos(lon) * x - Math.sin(lat) * Math.sin(lon) * y + Math.cos(lat) * z;
    }
  }
  const out = Object.fromEntries(['theta', 'T', 'east', 'north', 'omega'].map((n) => [n, new Float64Array(J * C)]));
  const massDivergence = new Float64Array(K), divergence = new Float64Array(K), flows = new Float64Array(K + 1);
  for (let i = 0; i < C; i++) {
    massDivergence.fill(0); divergence.fill(0);
    for (let m = 0; m < nEdgesOnCell[i]; m++) {
      const e = edgesOnCell[maxEdges * i + m], f = edgeSignOnCell[maxEdges * i + m] * dvEdge[e];
      const piEdge = 0.5 * (pi[cellsOnEdge[2 * e]] + pi[cellsOnEdge[2 * e + 1]]);
      for (let k = 0; k < K; k++) { const flow = f * u[k * E + e]; massDivergence[k] += piEdge * flow; divergence[k] += flow; }
    }
    let dPi = 0;
    for (let k = 0; k < K; k++) { massDivergence[k] /= areaCell[i]; divergence[k] /= areaCell[i]; dPi -= massDivergence[k] * dSigma[k]; }
    for (let k = 0, cumulative = 0; k < K; k++) { cumulative += massDivergence[k] * dSigma[k]; flows[k + 1] = k === K - 1 ? 0 : -cumulative - sigmaLower[k] * dPi; }
    for (let j = 0; j < J; j++) {
      const p = pressures[j], sigma = p / pi[i];
      let k = 0;
      while (k < K - 2 && pi[i] * sigmaMid[k + 1] < p) k++;
      const t = Math.min(1, Math.max(0, Math.log(p / (pi[i] * sigmaMid[k])) / Math.log(sigmaMid[k + 1] / sigmaMid[k])));
      const at = (x) => x[k * C + i] + t * (x[(k + 1) * C + i] - x[k * C + i]);
      let layer = 0;
      while (layer < K - 1 && sigmaLower[layer] < sigma) layer++;
      const s = Math.min(1, Math.max(0, (sigma - sigmaUpper[layer]) / dSigma[layer]));
      const advection = (n) => massDivergence[n] - pi[i] * divergence[n];
      const th = at(theta), temperature = (1 - t) * theta[k * C + i] * (pi[i] / P0) ** KAPPA * shape[k] + t * theta[(k + 1) * C + i] * (pi[i] / P0) ** KAPPA * shape[k + 1];
      out.theta[j * C + i] = th;
      out.T[j * C + i] = temperature;
      out.east[j * C + i] = at(east);
      out.north[j * C + i] = at(north);
      out.omega[j * C + i] = flows[layer] + s * (flows[layer + 1] - flows[layer]) + sigma * (dPi + advection(k) + t * (advection(k + 1) - advection(k)));
    }
  }
  return out;
}

// Per band b and pressure j (index j * bands + b): the zonal means, eddy
// covariances and the residual circulation.
export function residualCirculation(mesh, columns, pressures, bandWidth = 2.5) {
  const C = mesh.nCells, J = pressures.length, bands = Math.round(180 / bandWidth), a = mesh.radius ?? 6371220;
  const band = Int32Array.from(mesh.latCell, (lat) => Math.min(bands - 1, Math.floor((lat * DEG + 90) / bandWidth)));
  const lat = Float64Array.from({ length: bands }, (_, b) => (-90 + (b + 0.5) * bandWidth) / DEG);
  const area = new Float64Array(bands);
  for (let i = 0; i < C; i++) area[band[i]] += mesh.areaCell[i];
  const mean = (x, j) => { const m = new Float64Array(bands); for (let i = 0; i < C; i++) m[band[i]] += mesh.areaCell[i] * x[j * C + i]; return m.map((v, b) => v / area[b]); };
  const covariance = (x, y, j, mx, my) => { const m = new Float64Array(bands); for (let i = 0; i < C; i++) m[band[i]] += mesh.areaCell[i] * x[j * C + i] * y[j * C + i]; return m.map((v, b) => v / area[b] - mx[b] * my[b]); };
  const field = () => new Float64Array(J * bands);
  const r = Object.fromEntries(['u', 'v', 'omega', 'theta', 'T', 'vTheta', 'uv', 'uOmega', 'thetaP', 'psi', 'vStar', 'omegaStar', 'wBar', 'wStar', 'cooling'].map((n) => [n, field()]));
  for (let j = 0; j < J; j++) {
    const u = mean(columns.east, j), v = mean(columns.north, j), w = mean(columns.omega, j), th = mean(columns.theta, j), T = mean(columns.T, j);
    const vTheta = covariance(columns.north, columns.theta, j, v, th), uv = covariance(columns.east, columns.north, j, u, v), uOmega = covariance(columns.east, columns.omega, j, u, w);
    for (let b = 0; b < bands; b++) {
      const n = j * bands + b;
      r.u[n] = u[b]; r.v[n] = v[b]; r.omega[n] = w[b]; r.theta[n] = th[b]; r.T[n] = T[b]; r.vTheta[n] = vTheta[b]; r.uv[n] = uv[b]; r.uOmega[n] = uOmega[b];
    }
  }
  const dp = (x, j, b) => {
    const lo = Math.max(0, j - 1), hi = Math.min(J - 1, j + 1);
    return (x[hi * bands + b] - x[lo * bands + b]) / Math.log(pressures[hi] / pressures[lo]) / pressures[j];
  };
  for (let j = 0; j < J; j++) for (let b = 0; b < bands; b++) {
    const n = j * bands + b;
    r.thetaP[n] = dp(r.theta, j, b);
    r.psi[n] = r.vTheta[n] / r.thetaP[n];
  }
  for (let j = 0; j < J; j++) for (let b = 0; b < bands; b++) {
    const n = j * bands + b, lo = Math.max(0, b - 1), hi = Math.min(bands - 1, b + 1);
    const flux = (r.psi[j * bands + hi] * Math.cos(lat[hi]) - r.psi[j * bands + lo] * Math.cos(lat[lo])) / (lat[hi] - lat[lo]);
    r.omegaStar[n] = r.omega[n] + flux / (a * Math.cos(lat[b]));
    r.vStar[n] = r.v[n] - dp(r.psi, j, b);
    r.wBar[n] = -H * r.omega[n] / pressures[j];
    r.wStar[n] = -H * r.omegaStar[n] / pressures[j];
    r.cooling[n] = 86400 * r.omegaStar[n] * r.thetaP[n] * r.T[n] / r.theta[n];
  }
  return { bands, lat, area, ...r };
}

export function bandMean(tem, x, from, to, j) {
  let sum = 0, weight = 0;
  for (let b = 0; b < tem.bands; b++) {
    const l = tem.lat[b] * DEG;
    if (l < from || l > to) continue;
    sum += tem.area[b] * x[j * tem.bands + b]; weight += tem.area[b];
  }
  return sum / weight;
}

export function temSums() {
  let sums = null, meta = null, count = 0;
  return {
    add(tem) {
      if (!sums) sums = Object.fromEntries(Object.entries(tem).filter(([n, x]) => x instanceof Float64Array && !['lat', 'area'].includes(n)).map(([n, x]) => [n, new Float64Array(x.length)]));
      for (const n of Object.keys(sums)) for (let m = 0; m < sums[n].length; m++) sums[n][m] += tem[n][m];
      meta ??= { bands: tem.bands, lat: tem.lat, area: tem.area };
      count++;
    },
    mean() { return { ...meta, count, ...Object.fromEntries(Object.entries(sums).map(([n, x]) => [n, x.map((v) => v / count)])) }; },
  };
}

export function tropicalTable(tem, pressures, tropics = 15, bandWidth = 2.5) {
  const lines = [`means of ${tem.count} samples over ${tropics}°S-${tropics}°N (${bandWidth}° bands): T, K; u at 5°S-5°N, m/s; w-bar, w* and its eddy part, mm/s (H = 7 km); v* at ${tropics}°N and ${tropics}°S, mm/s; the ascent's adiabatic cooling, K/day`, '    hPa      T   u_eq    w-bar     w*  eddy   v*N    v*S  cooling'];
  pressures.forEach((p, j) => {
    const m = (x, from = -tropics, to = tropics) => bandMean(tem, x, from, to, j);
    const w = m(tem.wStar), wb = m(tem.wBar);
    lines.push(`  ${(p / 100).toPrecision(3).padStart(5)} ${m(tem.T).toFixed(1).padStart(6)} ${m(tem.u, -5, 5).toFixed(1).padStart(6)} ${(1e3 * wb).toFixed(3).padStart(8)} ${(1e3 * w).toFixed(3).padStart(6)} ${(1e3 * (w - wb)).toFixed(3).padStart(6)} ${(1e3 * m(tem.vStar, tropics - bandWidth, tropics + bandWidth)).toFixed(0).padStart(6)} ${(1e3 * m(tem.vStar, -tropics - bandWidth, -tropics + bandWidth)).toFixed(0).padStart(6)} ${m(tem.cooling).toFixed(2).padStart(8)}`);
  });
  return lines.join('\n');
}

export function temJson(tem, extra = {}) {
  return JSON.stringify({ ...extra, count: tem.count, lat: Array.from(tem.lat, (x) => x * DEG), ...Object.fromEntries(Object.entries(tem).filter(([n, x]) => x instanceof Float64Array && n !== 'lat').map(([n, x]) => [n, Array.from(x)])) });
}

async function main(files) {
  const pressures = (process.env.PRESSURES ? process.env.PRESSURES.split(',').map(Number) : TEM_PRESSURES).map((p) => 100 * p);
  const bandWidth = Number(process.env.BAND ?? 2.5), tropics = Number(process.env.TROPICS ?? 15), sums = temSums();
  let mesh = null, core = null, N = 0, relax = null;
  for (const file of files) {
    const saved = await decodeState(new Uint8Array(readFileSync(file)));
    const found = stateRelax(saved).relax;
    if (!mesh) { N = saved.N; relax = found; mesh = stateGrid(saved).mesh; core = createSigmaCore(mesh, { levels: savedLevels(saved) }); }
    else if (saved.N !== N || found !== relax) throw new Error(`${file} is not on ${files[0]}'s grid`);
    const state = [Float64Array.from(saved.pi), Float64Array.from(saved.theta), Float64Array.from(saved.u)];
    sums.add(residualCirculation(mesh, pressureColumns(mesh, core, state, pressures), pressures, bandWidth));
    console.log(`${file.replace(/.*\//, '')}: day ${saved.day}`);
  }
  const tem = sums.mean();
  console.log(tropicalTable(tem, pressures, tropics, bandWidth));
  if (process.env.OUT) writeFileSync(process.env.OUT, temJson(tem, { files, pressures, bandWidth }));
}

if (import.meta.url === `file://${process.argv[1]}`) await main(process.argv.slice(2));
