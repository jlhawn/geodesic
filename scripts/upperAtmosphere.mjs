// The layers above 200 hPa of a state: zonal means by latitude band, the
// zonally asymmetric part, the largest winds and the step's Courant numbers.
//   node scripts/upperAtmosphere.mjs <state.bin> [<state.bin> ...]
// prints, for each state, every such layer's zonal-mean zonal wind and
// temperature in 10-degree bands, its eddy kinetic energy and eddy
// temperature (the departures from the band means), the largest edge wind,
// the rms divergence, and the horizontal and vertical Courant numbers at the
// N's time step (dt = 1350 * 16 / N); with several states, the means over
// them. SURFACE and GRAVITY_WAVES (JSON, as for scripts/spinup.mjs; the
// model's defaults where unset) name the top treatment, whose zonal force
// the tables give, by layer and band, from the CPU engine's own operators
// on each state, with the gravity waves' absolute momentum flux by layer
// and band (unless GRAVITY_WAVES is false) and the mass flux the force drives
// poleward of 50 degrees (f v* = -F, the force's downward control,
// Haynes et al. 1991) and the adiabatic warming its descent gives each
// layer of the 70-90 degree caps.
// spinup.mjs's STRATOSPHERE line takes upperWindLine from here.
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { cellVector, divergence } from '../js/dynamics/operators.module.js';
import { R_DRY, CP_DRY, P0, GRAVITY, createSigmaCore, sigmaGridName } from '../js/dynamics/sigmaCore.module.js';
import { decodeState, savedLevels } from '../js/stateFile.module.js';
import { SIDEREAL_DAY } from '../js/model.module.js';
import { createSurface, TOP_DRAG } from '../js/physics/surface.module.js';
import { spongeGeometry, spongeRates, dampEddies, spongeSigmaFor, SPONGE } from '../js/dynamics/sponge.module.js';
import { createGravityWaveDrag } from '../js/physics/gravityWaves.module.js';

const KAPPA = R_DRY / CP_DRY, DEG = 180 / Math.PI;

export function upperLayers(levels, below = 0.2) {
  const out = [];
  for (let k = 0; k < levels.length - 1 && 0.5 * (levels[k] + levels[k + 1]) < below; k++) out.push(k);
  return out;
}

// Per layer: the cells' east and north wind and temperature, their band
// means, the eddy kinetic energy and temperature variance about them.
export function upperAtmosphere(mesh, levels, state, { bandWidth = 10, dt = null, layers = upperLayers(levels) } = {}) {
  const C = mesh.nCells, E = mesh.nEdges, K = levels.length - 1;
  const [pi, theta, u] = state;
  const bands = Math.round(180 / bandWidth), band = Int32Array.from(mesh.latCell, (lat) => Math.min(bands - 1, Math.floor((lat * DEG + 90) / bandWidth)));
  const east = new Float64Array(C), north = new Float64Array(C), vector = new Float64Array(3 * C);
  const flux = new Float64Array(E), div = new Float64Array(C), columnDiv = new Float64Array(C), cumulative = new Float64Array(C);
  let minDc = Infinity;
  for (let e = 0; e < E; e++) minDc = Math.min(minDc, mesh.dcEdge[e]);
  const allDiv = [];
  for (let k = 0; k < K; k++) {
    for (let e = 0; e < E; e++) flux[e] = 0.5 * (pi[mesh.cellsOnEdge[2 * e]] + pi[mesh.cellsOnEdge[2 * e + 1]]) * u[k * E + e];
    divergence(mesh, flux, div);
    allDiv.push(Float64Array.from(div));
    for (let i = 0; i < C; i++) columnDiv[i] += div[i] * (levels[k + 1] - levels[k]);
  }
  const rows = [];
  for (let k = 0, kk = 0; k < K && kk < layers.length; k++) {
    const dSigma = levels[k + 1] - levels[k];
    let worstW = 0;
    for (let i = 0; i < C; i++) {
      cumulative[i] += allDiv[k][i] * dSigma;
      const above = -(cumulative[i] - allDiv[k][i] * dSigma) + levels[k] * columnDiv[i], belowFlow = -cumulative[i] + levels[k + 1] * columnDiv[i];
      worstW = Math.max(worstW, Math.max(Math.abs(above), Math.abs(belowFlow)) / (pi[i] * dSigma));
    }
    if (k !== layers[kk]) continue;
    kk++;
    const shape = (levels[k + 1] ** (1 + KAPPA) - levels[k] ** (1 + KAPPA)) / ((1 + KAPPA) * dSigma);
    cellVector(mesh, u.subarray(k * E, (k + 1) * E), vector);
    const area = new Float64Array(bands), uMean = new Float64Array(bands), vMean = new Float64Array(bands), tMean = new Float64Array(bands);
    const T = new Float64Array(C);
    for (let i = 0; i < C; i++) {
      const lat = mesh.latCell[i], lon = mesh.lonCell[i], x = vector[3 * i], y = vector[3 * i + 1], z = vector[3 * i + 2];
      east[i] = -Math.sin(lon) * x + Math.cos(lon) * y;
      north[i] = -Math.sin(lat) * Math.cos(lon) * x - Math.sin(lat) * Math.sin(lon) * y + Math.cos(lat) * z;
      T[i] = theta[k * C + i] * (pi[i] / P0) ** KAPPA * shape;
      const a = mesh.areaCell[i], b = band[i];
      area[b] += a; uMean[b] += a * east[i]; vMean[b] += a * north[i]; tMean[b] += a * T[i];
    }
    for (let b = 0; b < bands; b++) { uMean[b] /= area[b]; vMean[b] /= area[b]; tMean[b] /= area[b]; }
    let total = 0, eke = 0, tVar = 0, divSq = 0;
    for (let i = 0; i < C; i++) {
      const a = mesh.areaCell[i], b = band[i];
      total += a;
      eke += a * 0.5 * ((east[i] - uMean[b]) ** 2 + (north[i] - vMean[b]) ** 2);
      tVar += a * (T[i] - tMean[b]) ** 2;
      divSq += a * allDiv[k][i] ** 2 / pi[i] ** 2;
    }
    let maxWind = 0;
    for (let e = 0; e < E; e++) maxWind = Math.max(maxWind, Math.abs(u[k * E + e]));
    rows.push({
      k, pressure: 500 * (levels[k] + levels[k + 1]), bandWidth, uMean, vMean, tMean, east: Float64Array.from(east), north: Float64Array.from(north), T,
      eke: eke / total, tRms: Math.sqrt(tVar / total), divRms: Math.sqrt(divSq / total), maxWind,
      ...(dt ? { courant: maxWind * dt / minDc, verticalCourant: worstW * dt } : {}),
    });
  }
  return { rows, bands, band };
}

const jet = (row, south) => {
  const n = row.uMean.length;
  let best = -Infinity, at = 0;
  for (let b = 0; b < n; b++) {
    const lat = -90 + (b + 0.5) * row.bandWidth;
    if ((south ? lat > -20 : lat < 20) || row.uMean[b] <= best) continue;
    best = row.uMean[b]; at = lat;
  }
  return `${best.toFixed(0)}@${at.toFixed(0)}`;
};

const bandMean = (row, from, to) => {
  let sum = 0, n = 0;
  row.uMean.forEach((u, b) => { const lat = -90 + (b + 0.5) * row.bandWidth; if (lat > from && lat < to) { sum += u; n++; } });
  return (sum / n).toFixed(0);
};

// One line: per layer above 200 hPa, the strongest zonal-mean westerly
// poleward of 20 degrees in each hemisphere (m/s at the band's latitude,
// 2.5-degree bands), the zonal-mean wind within 5 degrees of the
// equator and over 57.5-62.5 S and N, the largest edge wind, the eddy
// kinetic energy, the eddy temperature's rms, the rms divergence and the
// Courant numbers.
export function upperWindLine(mesh, levels, state, dt, day) {
  const { rows } = upperAtmosphere(mesh, levels, state, { bandWidth: 2.5, dt });
  return `upper winds day ${day} (jet N / jet S / 5S-5N / 60S / 60N, m/s@lat; max wind m/s; eddy KE m²/s²; eddy T rms K; divergence rms 1e-6/s; Courant horizontal/vertical): ` + rows.map((r) => `${r.pressure.toPrecision(3)} hPa ${jet(r, false)}/${jet(r, true)}/${bandMean(r, -5, 5)}/${bandMean(r, -62.5, -57.5)}/${bandMean(r, 57.5, 62.5)} ${r.maxWind.toFixed(0)} ${r.eke.toFixed(0)} ${r.tRms.toFixed(1)} ${(1e6 * r.divRms).toFixed(1)} ${r.courant.toFixed(2)}/${r.verticalCourant.toFixed(2)}`).join('; ');
}

// The zonal mean of the east component of a treatment's wind tendency.
function zonalForce(mesh, tendency, layers, band, bands) {
  const C = mesh.nCells, E = mesh.nEdges, vector = new Float64Array(3 * C);
  return layers.map((k) => {
    cellVector(mesh, tendency.subarray(k * E, (k + 1) * E), vector);
    const sum = new Float64Array(bands), area = new Float64Array(bands);
    for (let i = 0; i < C; i++) {
      const lon = mesh.lonCell[i], a = mesh.areaCell[i];
      sum[band[i]] += a * (-Math.sin(lon) * vector[3 * i] + Math.cos(lon) * vector[3 * i + 1]);
      area[band[i]] += a;
    }
    return sum.map((s, b) => s / area[b]);
  });
}

// The wind tendency of the top treatment on the CPU engine: the Rayleigh drag
// and the sponge that SURFACE names and the gravity-wave drag of
// GRAVITY_WAVES, each with the model's defaults where they are not named.
function topTendency(mesh, levels, state, options, waves) {
  const core = createSigmaCore(mesh, { levels });
  core.diagnose(state[0], state[1]);
  const surface = createSurface(mesh, core, { topSigma: TOP_DRAG.sigma, topDragDays: TOP_DRAG.days, ...options });
  const K = levels.length - 1, E = mesh.nEdges, out = [new Float64Array(mesh.nCells), new Float64Array(K * mesh.nCells), new Float64Array(K * E)];
  surface.applyTop(state, out);
  const rates = spongeRates(core.sigmaMid, options.spongeSigma ?? spongeSigmaFor(sigmaGridName(levels)), options.spongeDays ?? SPONGE.days);
  if (rates.some((r) => r > 0)) {
    const geometry = spongeGeometry(mesh), means = new Float64Array(2 * geometry.bands), layer = new Float64Array(E);
    rates.forEach((rate, k) => {
      if (!(rate > 0)) return;
      layer.set(state[2].subarray(k * E, (k + 1) * E));
      dampEddies(mesh, geometry, layer, rate, 1, means);
      for (let e = 0; e < E; e++) out[2][k * E + e] += layer[e] - state[2][k * E + e];
    });
  }
  if (waves === false) return { tendency: out[2], absoluteFlux: null };
  const drag = createGravityWaveDrag(mesh, core, { ...waves, diagnose: true });
  drag.compute(state);
  drag.applyEdges(out[2], 0, E, 1);
  return { tendency: out[2], absoluteFlux: drag.absoluteFlux };
}

async function main(files) {
  const SURFACE = JSON.parse(process.env.SURFACE ?? '{}'), WAVES = JSON.parse(process.env.GRAVITY_WAVES ?? '{}');
  let mesh = null, N = 0, sums = null, count = 0, levels = null, fluxSums = null;
  for (const file of files) {
    const saved = await decodeState(new Uint8Array(readFileSync(file)));
    if (!mesh) { N = saved.N; mesh = buildMesh(new Grid(N), { omega: 2 * Math.PI / SIDEREAL_DAY }); levels = savedLevels(saved); }
    const dt = 1350 * 16 / N, state = [Float64Array.from(saved.pi), Float64Array.from(saved.theta), Float64Array.from(saved.u)];
    const { rows, bands, band } = upperAtmosphere(mesh, levels, state, { dt });
    const top = topTendency(mesh, levels, state, SURFACE, WAVES), force = zonalForce(mesh, top.tendency, rows.map((r) => r.k), band, bands);
    if (top.absoluteFlux) {
      fluxSums ??= rows.map(() => new Float64Array(bands));
      rows.forEach((r, j) => {
        const sum = new Float64Array(bands), area = new Float64Array(bands);
        for (let i = 0; i < mesh.nCells; i++) { sum[band[i]] += mesh.areaCell[i] * top.absoluteFlux[r.k * mesh.nCells + i]; area[band[i]] += mesh.areaCell[i]; }
        for (let b = 0; b < bands; b++) fluxSums[j][b] += sum[b] / area[b];
      });
    }
    console.log(`${file}: day ${saved.day}, N=${N}, ${levels.length - 1} layers; max wind / eddy KE / eddy T rms / divergence rms / Courant`);
    for (const r of rows) console.log(`  ${r.pressure.toPrecision(3)} hPa: ${r.maxWind.toFixed(0)} m/s / ${r.eke.toFixed(0)} m²/s² / ${r.tRms.toFixed(2)} K / ${(1e6 * r.divRms).toFixed(2)}e-6 /s / ${r.courant.toFixed(2)}, ${r.verticalCourant.toFixed(3)}`);
    if (!sums) sums = rows.map(() => ({ u: new Float64Array(bands), T: new Float64Array(bands), F: new Float64Array(bands), eke: 0, tRms: 0 }));
    rows.forEach((r, j) => { for (let b = 0; b < bands; b++) { sums[j].u[b] += r.uMean[b]; sums[j].T[b] += r.tMean[b]; sums[j].F[b] += force[j][b]; } sums[j].eke += r.eke; sums[j].tRms += r.tRms; sums[j].p = r.pressure; sums[j].k = r.k; });
    count++;
  }
  const bands = sums[0].u.length, lats = Array.from({ length: bands }, (_, b) => -90 + (b + 0.5) * 180 / bands);
  const head = lats.map((l) => l.toFixed(0).padStart(6)).join('');
  const table = (title, pick, digits) => {
    console.log(`${title}${count > 1 ? ` (mean of ${count} states)` : ''}\n  hPa     ${head}`);
    for (const s of sums) console.log(`  ${s.p.toPrecision(3).padStart(6)}  ${Array.from(pick(s), (x) => (x / count).toFixed(digits).padStart(6)).join('')}`);
  };
  table('zonal-mean zonal wind, m/s, by band centre', (s) => s.u, 1);
  table('zonal-mean temperature, K', (s) => s.T, 1);
  if (fluxSums) table('the gravity waves\' absolute momentum flux rising through each layer, both directions, mPa', (s) => fluxSums[sums.indexOf(s)].map((x) => 1e3 * x), 2);
  table(`top treatment's zonal force, m/s/day (SURFACE ${JSON.stringify(SURFACE)}, GRAVITY_WAVES ${JSON.stringify(WAVES)})`, (s) => s.F.map((x) => 86400 * x), 2);
  const omega = 2 * Math.PI / SIDEREAL_DAY, a = mesh.radius ?? 6371220;
  for (const [name, edge, sign] of [['north', 50, 1], ['south', -50, -1]]) {
    const capT = (x) => { let t = 0, w = 0; lats.forEach((lat, b) => { if (sign * lat > 70) { t += x.T[b] / count; w++; } }); return t / w; };
    let above = 0;
    const cells = sums.map((s, j) => {
      const k = s.k, dp = (levels[k + 1] - levels[k]) * P0, p = 0.5 * (levels[k] + levels[k + 1]) * P0;
      let force = 0, weight = 0;
      lats.forEach((lat, b) => { if (Math.abs(lat - edge) < 5.01) { force += s.F[b] / count; weight++; } });
      const layerFlux = -sign * (force / weight) / (2 * omega * Math.sin(edge / DEG)) * dp / GRAVITY;
      const poleward = above + 0.5 * layerFlux;
      above += layerFlux;
      const descent = poleward * 2 * Math.PI * a * Math.cos(edge / DEG) / (2 * Math.PI * a * a * (1 - Math.sin(Math.abs(edge) / DEG)));
      const up = sums[Math.max(0, j - 1)], down = sums[Math.min(sums.length - 1, j + 1)];
      const pUp = 0.5 * (levels[up.k] + levels[up.k + 1]) * P0, pDown = 0.5 * (levels[down.k] + levels[down.k + 1]) * P0;
      const stability = KAPPA * capT(s) / p - (capT(down) - capT(up)) / (pDown - pUp);
      return `${s.p.toPrecision(3)} hPa ${poleward.toFixed(0)} kg/m/s, ${(86400 * GRAVITY * descent * stability).toFixed(2)} K/day`;
    });
    console.log(`the force's poleward mass flux across ${Math.abs(edge)}°${name === 'north' ? 'N' : 'S'} (f v* = -F over the layers above each layer's middle) and the adiabatic warming of its descent over the cap at the 70-90° cap's static stability: ${cells.join('; ')}`);
  }
}

if (import.meta.url === `file://${process.argv[1]}`) await main(process.argv.slice(2));
