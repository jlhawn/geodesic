// The atmosphere's boundary-layer zonal momentum budget on the equator over
// the Pacific by 20 degrees of longitude, from one saved state, on the CPU
// with the ocean off:
//   node scripts/equatorialMomentum.mjs <state.bin>
// SPIN (1) full steps are taken first and the budget is the next step.
// Terms are the east component of each edge tendency reconstructed at the
// cells (1e-5 m/s2), area means over the sea cells within LAT (3) degrees:
// the dynamics' terms from the step's starting state (pgf the sigma
// pressure-gradient force, corF and corZeta the f and relative-vorticity
// parts of the PV flux, gradK, vertical advection, drag the surface drag on
// the lowest layer; sumK1 their sum, beside the core's own tendency), then
// the step's actual increments over dt: rk4, the closure split into nu4 and
// divDamp, blMix the boundary layer's implicit mixing, plume the plume's
// momentum transport, total. BOUNDARY_LAYER and SURFACE (JSON) pass options.
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { createModel, STATE_NAMES } from '../js/model.module.js';
import { SEA_DRAG } from '../js/physics/surface.module.js';
import { decodeState, savedLevels } from '../js/stateFile.module.js';
import { savedDeckField, DECK_FIELDS } from '../js/physics/regrid.module.js';
import { createRK4Arrays } from '../js/dynamics/integrators.module.js';
import { divergence, gradient, curl, kineticEnergy, laplacianVelocity, cellVector } from '../js/dynamics/operators.module.js';

const FILE = process.argv[2], LAT = Number(process.env.LAT ?? 3);
const BOUNDARY_LAYER = JSON.parse(process.env.BOUNDARY_LAYER ?? '{}'), SURFACE = JSON.parse(process.env.SURFACE ?? '{}');
const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const model = createModel(new Grid(saved.N), { topography, levels: savedLevels(saved), ocean: false, boundaryLayer: BOUNDARY_LAYER, surface: SURFACE });
const { mesh, core, state, phases, radiation, boundaryLayer: bl, moist, seaIce, land, surface } = model;
const { K, C, E, V, sigmaMid, sigmaLower, dSigma, R, g, exnerLayer, geopotential, piSigmaDot } = core.diagnostics;
const { thetaV } = core.arrays;
STATE_NAMES.forEach((name, a) => state[a].set(saved[name]));
seaIce.load(state[6], saved.concentration ?? null);
for (const field of Object.keys(DECK_FIELDS)) radiation[field].set(savedDeckField(saved, field, model));
land.load({ soil: Float64Array.from(saved.land.soil), snow: Float64Array.from(saved.land.snow), ...(saved.land.vegetation ? { vegetation: Float64Array.from(saved.land.vegetation) } : {}), ...(saved.land.surface ? { surface: Float64Array.from(saved.land.surface) } : {}) }, state[6]);
model.time = saved.time;
const dt = 1350 * 16 / saved.N, deg = 180 / Math.PI;
const { latCell, lonCell, cellsOnEdge, verticesOnEdge, nEdgesOnEdge, edgesOnEdge, weightsOnEdge, maxEdgesOnEdge, dcEdge, dvEdge, fVertex, areaCell } = mesh;
const landMask = model.geography.land;
const rk4 = createRK4Arrays(STATE_NAMES.map((_, a) => state[a].length));

function step() {
  radiation.setTime(model.time);
  rk4(model.tendency, state, dt);
  phases.physics(0, C, dt, model.totals);
  phases.closure(0, K, dt);
  phases.adjust(0, C, dt);
  phases.mixMomentum(0, E, dt);
  phases.dissipate(0, C);
  model.time += dt;
}
for (let n = 0; n < Number(process.env.SPIN ?? 1); n++) step();

const BINS = [[140, 160], [160, 180], [180, 200], [200, 220], [220, 240], [240, 260], [260, 280]];
const binName = ([a, b]) => `${a <= 180 ? a + 'E' : 360 - a + 'W'}-${b <= 180 ? b + 'E' : 360 - b + 'W'}`.replace('180E', '180').replace('180W', '180');
const binOf = new Int8Array(C).fill(-1);
for (let i = 0; i < C; i++) {
  if (landMask[i] || state[6][i] > 0 || Math.abs(latCell[i] * deg) > LAT) continue;
  const lon = ((lonCell[i] * deg) % 360 + 360) % 360;
  binOf[i] = BINS.findIndex(([a, b]) => lon >= a && lon < b);
}
const cells = []; for (let i = 0; i < C; i++) if (binOf[i] >= 0) cells.push(i);
const NB = BINS.length;
const east = (vec, i) => -Math.sin(lonCell[i]) * vec[3 * i] + Math.cos(lonCell[i]) * vec[3 * i + 1];
const north = (vec, i) => -Math.sin(latCell[i]) * (Math.cos(lonCell[i]) * vec[3 * i] + Math.sin(lonCell[i]) * vec[3 * i + 1]) + Math.cos(latCell[i]) * vec[3 * i + 2];
const vec = new Float64Array(3 * C);
function eastAt(edgeField, out) {
  for (const i of cells) {
    let x = 0, y = 0, z = 0;
    for (let m = 0; m < mesh.nEdgesOnCell[i]; m++) {
      const e = mesh.edgesOnCell[mesh.maxEdges * i + m], w = 0.5 * dcEdge[e] * dvEdge[e] * edgeField[e];
      x += w * mesh.nEdge[3 * e]; y += w * mesh.nEdge[3 * e + 1]; z += w * mesh.nEdge[3 * e + 2];
    }
    vec[3 * i] = x / areaCell[i]; vec[3 * i + 1] = y / areaCell[i]; vec[3 * i + 2] = z / areaCell[i];
    out[i] = east(vec, i);
  }
  return out;
}
function northAt(edgeField, out) { eastAt(edgeField, out); for (const i of cells) out[i] = north(vec, i); return out; }

const [pi, theta, u] = state;
const surfaceT = Float64Array.from(state[3]);
const u0 = Float64Array.from(u);
core.phaseFlux(state, 0, K);
core.phaseColumn(state, [new Float64Array(C)], 0, C);
core.phaseVertex(state, 0, V);
const flux = new Float64Array(core.shared.flux), piVertex = new Float64Array(core.shared.piVertex);
const psd = Float64Array.from(piSigmaDot), geo = Float64Array.from(geopotential), thv = Float64Array.from(thetaV), exL = Float64Array.from(exnerLayer);
const k1 = STATE_NAMES.map((_, a) => new Float64Array(state[a].length));
core.phaseLayer(state, k1, 0, K, 'momentum');
const dragOut = [null, null, new Float64Array(K * E)];
surface.lowestWindSpeed(u);
surface.applyLayers(state, dragOut, 0, K);
const stressEdge = surface.stress(state, new Float64Array(E));

const TERMS = ['pgf', 'corF', 'corZeta', 'gradK', 'vertical', 'drag', 'sumK1', 'phaseLayer+drag', 'rk4', 'nu4', 'divDamp', 'closure', 'blMix', 'plume', 'total'];
const T = Object.fromEntries(TERMS.map((t) => [t, new Float64Array(K * C)]));
const uEast = new Float64Array(K * C), vNorth = new Float64Array(K * C), zMid = new Float64Array(K * C);
const piEdge = Float64Array.from({ length: E }, (_, e) => 0.5 * (pi[cellsOnEdge[2 * e]] + pi[cellsOnEdge[2 * e + 1]]));
const lnPi = Float64Array.from(pi, Math.log), gradLnPi = gradient(mesh, lnPi);
const zeta = new Float64Array(V), qF = new Float64Array(V), qZ = new Float64Array(V), qEF = new Float64Array(E), qEZ = new Float64Array(E), kin = new Float64Array(C), gK = new Float64Array(E), gPhi = new Float64Array(E);
const fields = Object.fromEntries(['pgf', 'corF', 'corZeta', 'gradK', 'vertical'].map((t) => [t, new Float64Array(E)]));
const scratch = new Float64Array(C), tmpE = new Float64Array(E), sumE = new Float64Array(E);
for (let k = 0; k < K; k++) {
  const off = k * C, uk = u0.subarray(k * E, (k + 1) * E), fk = flux.subarray(k * E, (k + 1) * E);
  curl(mesh, uk, zeta);
  for (let v = 0; v < V; v++) { qF[v] = fVertex[v] / piVertex[v]; qZ[v] = zeta[v] / piVertex[v]; }
  for (let e = 0; e < E; e++) { qEF[e] = 0.5 * (qF[verticesOnEdge[2 * e]] + qF[verticesOnEdge[2 * e + 1]]); qEZ[e] = 0.5 * (qZ[verticesOnEdge[2 * e]] + qZ[verticesOnEdge[2 * e + 1]]); }
  kineticEnergy(mesh, uk, kin); gradient(mesh, kin, gK);
  gradient(mesh, geo.subarray(off, off + C), gPhi);
  for (let e = 0; e < E; e++) {
    const i = cellsOnEdge[2 * e], j = cellsOnEdge[2 * e + 1];
    let pf = 0, pz = 0;
    for (let s = 0; s < nEdgesOnEdge[e]; s++) {
      const slot = maxEdgesOnEdge * e + s, o = edgesOnEdge[slot], w = weightsOnEdge[slot] * dvEdge[o] * fk[o];
      pf += w * 0.5 * (qEF[e] + qEF[o]); pz += w * 0.5 * (qEZ[e] + qEZ[o]);
    }
    fields.corF[e] = pf / dcEdge[e]; fields.corZeta[e] = pz / dcEdge[e];
    fields.pgf[e] = -gPhi[e] - R * 0.5 * (thv[off + i] * exL[off + i] + thv[off + j] * exL[off + j]) * gradLnPi[e];
    fields.gradK[e] = -gK[e];
    const lowerFlow = 0.5 * (psd[(k + 1) * C + i] + psd[(k + 1) * C + j]), upperFlow = 0.5 * (psd[k * C + i] + psd[k * C + j]);
    const lowerU = k === K - 1 ? 0 : 0.5 * (uk[e] + u0[(k + 1) * E + e]), upperU = k === 0 ? 0 : 0.5 * (uk[e] + u0[(k - 1) * E + e]);
    fields.vertical[e] = -(lowerFlow * lowerU - upperFlow * upperU - uk[e] * (lowerFlow - upperFlow)) / (piEdge[e] * dSigma[k]);
  }
  const target = new Float64Array(C);
  for (const t of Object.keys(fields)) { eastAt(fields[t], target); for (const i of cells) T[t][off + i] = target[i]; }
  for (let e = 0; e < E; e++) sumE[e] = fields.pgf[e] + fields.corF[e] + fields.corZeta[e] + fields.gradK[e] + fields.vertical[e] + dragOut[2][k * E + e];
  eastAt(sumE, target); for (const i of cells) T.sumK1[off + i] = target[i];
  for (let e = 0; e < E; e++) tmpE[e] = k1[2][k * E + e] + dragOut[2][k * E + e];
  eastAt(tmpE, target); for (const i of cells) T['phaseLayer+drag'][off + i] = target[i];
  eastAt(dragOut[2].subarray(k * E, (k + 1) * E), target); for (const i of cells) T.drag[off + i] = target[i];
  eastAt(uk, target); for (const i of cells) uEast[off + i] = target[i];
  northAt(uk, target); for (const i of cells) vNorth[off + i] = target[i];
  for (const i of cells) zMid[off + i] = geo[off + i] / g;
}
const zs = Float64Array.from({ length: C }, (_, i) => (model.surfaceGeopotential ? model.surfaceGeopotential[i] / g : 0));
const lowRho = new Float64Array(C);
for (const i of cells) { const idx = (K - 1) * C + i; lowRho[i] = pi[i] * sigmaMid[K - 1] / (R * theta[idx] * exL[idx]); }
const tauEast = eastAt(stressEdge, new Float64Array(C));
const speed = Float64Array.from(surface.windSpeed);

radiation.setTime(model.time);
rk4(model.tendency, state, dt);
const u1 = Float64Array.from(u);
phases.physics(0, C, dt, model.totals);
const depth = Float64Array.from(bl.depth), mixing = Float64Array.from(bl.mixing), went = Float64Array.from(bl.entrainment), deckTop = Float64Array.from(radiation.mlmTop);
let spacing = 0; for (let e = 0; e < E; e++) spacing += dcEdge[e]; spacing /= E;
const nu4 = core.nu4, divStep = core.divergenceDamping * spacing * spacing;
const lap = new Float64Array(E), lap2 = new Float64Array(E), dS = new Float64Array(C), cS = new Float64Array(V), nuE = new Float64Array(E), ddE = new Float64Array(E);
const target = new Float64Array(C);
for (let k = 0; k < K; k++) {
  const uk = u1.subarray(k * E, (k + 1) * E);
  laplacianVelocity(mesh, uk, lap, dS, cS); laplacianVelocity(mesh, lap, lap2, dS, cS);
  for (let e = 0; e < E; e++) { nuE[e] = -nu4 * lap2[e]; tmpE[e] = uk[e] + dt * nuE[e]; }
  divergence(mesh, tmpE, dS); gradient(mesh, dS, ddE);
  for (let e = 0; e < E; e++) ddE[e] *= divStep / dt;
  eastAt(nuE, target); for (const i of cells) T.nu4[k * C + i] = target[i];
  eastAt(ddE, target); for (const i of cells) T.divDamp[k * C + i] = target[i];
}
phases.closure(0, K, dt);
const u2 = Float64Array.from(u);
phases.adjust(0, C, dt);
const u2b = Float64Array.from(u);
bl.mixEdges(state[0], u, 0, E, dt, core.arrays.dissipation);
const u3 = Float64Array.from(u);
moist.transportMomentum(state[0], u, 0, E, dt, core.arrays.dissipation);
const u4 = Float64Array.from(u);
for (let k = 0; k < K; k++) {
  const s = k * E;
  const put = (name, a, b) => { for (let e = 0; e < E; e++) tmpE[e] = (b[s + e] - a[s + e]) / dt; eastAt(tmpE, target); for (const i of cells) T[name][k * C + i] = target[i]; };
  put('rk4', u0, u1); put('closure', u1, u2); put('blMix', u2b, u3); put('plume', u3, u4); put('total', u0, u4);
}
let adjustChange = 0; for (let x = 0; x < u2.length; x++) adjustChange = Math.max(adjustChange, Math.abs(u2b[x] - u2[x]));

const f = (x, d = 2) => (Number.isFinite(x) ? x.toFixed(d) : '—');
const s5 = (x) => f(x * 1e5);
const binMean = (fn) => { const s = new Float64Array(NB), a = new Float64Array(NB); for (const i of cells) { const v = fn(i); if (!Number.isFinite(v)) continue; s[binOf[i]] += areaCell[i] * v; a[binOf[i]] += areaCell[i]; } return Array.from(s, (x, b) => x / a[b]); };
const hAbove = (i) => (deckTop[i] > 0 ? Math.max(depth[i], deckTop[i]) : depth[i]) - zs[i];
const blLayers = (i) => { const h = hAbove(i), ks = []; for (let k = K - 1; k >= 0; k--) if (zMid[k * C + i] - zs[i] < h) ks.push(k); return ks.length ? ks : [K - 1]; };
const blMass = (i) => blLayers(i).reduce((m, k) => m + pi[i] * dSigma[k] / g, 0);
const blMean = (arr, i) => { let s = 0, w = 0; for (const k of blLayers(i)) { s += dSigma[k] * arr[k * C + i]; w += dSigma[k]; } return s / w; };
const mixTop = (i) => { let top = 0; for (let k = 0; k < K - 1; k++) if (mixing[k * C + i] > 0) { top = 0.5 * (zMid[k * C + i] + zMid[(k + 1) * C + i]) - zs[i]; break; } return top; };
const atHeight = (arr, i, z) => { for (let k = K - 1; k > 0; k--) { const za = zMid[(k - 1) * C + i] - zs[i], zb = zMid[k * C + i] - zs[i]; if (z <= za) { const t = Math.max(0, Math.min(1, (z - zb) / (za - zb))); return arr[k * C + i] + t * (arr[(k - 1) * C + i] - arr[k * C + i]); } } return arr[i]; };
const slp = (i) => pi[i] * Math.exp(g * zs[i] / (R * theta[(K - 1) * C + i] * exL[(K - 1) * C + i]));
const count = binMean(() => 1).map((_, b) => cells.filter((i) => binOf[i] === b).length);
const head = '| | ' + BINS.map(binName).join(' | ') + ' |\n|---|' + BINS.map(() => '---|').join('');
const row = (label, vals, fmt = s5) => `| ${label} | ${vals.map(fmt).join(' | ')} |`;
console.log(`== ${FILE.split('/').pop()} day ${saved.day}, budget on step ${Number(process.env.SPIN ?? 1) + 1}, N=${saved.N}, K=${K}, dt ${dt} s, |lat| <= ${LAT}, sea cells ${count.join(', ')}; adjust phase max |du| ${adjustChange.toExponential(1)}`);
console.log(`nu4 ${nu4.toExponential(3)} m4/s, divergence damping nu_d ${(divStep / dt).toExponential(3)} m2/s, mean spacing ${(spacing / 1e3).toFixed(1)} km`);

console.log('\n-- state');
console.log(head);
const slpBins = binMean((i) => slp(i) / 100);
console.log(row('sea-level pressure (hPa)', slpBins, (x) => f(x, 2)));
console.log(row('surface temperature (C)', binMean((i) => surfaceT[i] - 273.15), (x) => f(x, 2)));
console.log(row('surface pressure (hPa)', binMean((i) => pi[i] / 100), (x) => f(x, 2)));
console.log(row('u lowest layer (m/s)', binMean((i) => uEast[(K - 1) * C + i]), (x) => f(x)));
console.log(row('v lowest layer (m/s)', binMean((i) => vNorth[(K - 1) * C + i]), (x) => f(x)));
console.log(row('|V| lowest layer (m/s)', binMean((i) => speed[i]), (x) => f(x)));
console.log(row('share of cells at the gust floor', binMean((i) => (speed[i] < 3 ? 1 : 0)), (x) => f(x)));
console.log(row('lowest-layer mid height (m)', binMean((i) => zMid[(K - 1) * C + i] - zs[i]), (x) => f(x, 0)));
console.log(row('u BL mean (m/s)', binMean((i) => blMean(uEast, i)), (x) => f(x)));
console.log(row('BL depth h, Richardson (m)', binMean((i) => depth[i] - zs[i]), (x) => f(x, 0)));
console.log(row('BL depth with the deck (m)', binMean(hAbove), (x) => f(x, 0)));
console.log(row('top of momentum mixing incl. entrainment (m)', binMean(mixTop), (x) => f(x, 0)));
console.log(row('entrainment w_e (mm/s)', binMean((i) => went[i] * 1e3), (x) => f(x)));
console.log(row('u at BL top (m/s)', binMean((i) => atHeight(uEast, i, hAbove(i))), (x) => f(x)));
console.log(row('u at 1.5 km (m/s)', binMean((i) => atHeight(uEast, i, 1500)), (x) => f(x)));
console.log(row('u at 3 km (m/s)', binMean((i) => atHeight(uEast, i, 3000)), (x) => f(x)));
console.log(row('stress tau_x (N/m2)', binMean((i) => tauEast[i]), (x) => f(x, 4)));
console.log(row('rho lowest layer (kg/m3)', binMean((i) => lowRho[i]), (x) => f(x, 3)));
console.log(row('BL mass (kg/m2)', binMean(blMass), (x) => f(x, 0)));

console.log('\n-- zonal acceleration, 1e-5 m/s2: lowest layer');
console.log(head);
for (const t of TERMS) console.log(row(t, binMean((i) => T[t][(K - 1) * C + i])));
console.log('\n-- zonal acceleration, 1e-5 m/s2: boundary-layer mass mean (layers whose mid lies below h)');
console.log(head);
for (const t of TERMS) console.log(row(t, binMean((i) => blMean(T[t], i))));
console.log(row('tau_x / BL mass', binMean((i) => tauEast[i] / blMass(i))));
console.log(row('Cd max(|V|,3) |u| u / h (lowest rho)', binMean((i) => (SURFACE.dragCoefficient ?? SEA_DRAG) * Math.max(speed[i], 3) * uEast[(K - 1) * C + i] / hAbove(i))));

console.log('\n-- zonal SLP gradient between bin centres and the acceleration -(1/rho) dp/dx it implies (1e-5 m/s2)');
const R_E = mesh.radius ?? 6.371e6;
const dx = 20 * Math.PI / 180 * R_E;
console.log('| ' + BINS.slice(0, -1).map((b, j) => `${binName(b)} to ${binName(BINS[j + 1])}`).join(' | ') + ' |\n|' + BINS.slice(1).map(() => '---|').join(''));
console.log('| ' + BINS.slice(0, -1).map((_, j) => s5(-(slpBins[j + 1] - slpBins[j]) * 100 / (1.18 * dx))).join(' | ') + ' |');
console.log(`basin 140E-80W: SLP ${f(slpBins[0])} -> ${f(slpBins[NB - 1])} hPa, max ${f(Math.max(...slpBins))}, min ${f(Math.min(...slpBins))}`);

console.log('\n-- profiles by layer (bins pooled 180-100W and 140E-180)');
for (const [label, set] of [['180-100W', [2, 3, 4, 5]], ['140E-180', [0, 1]]]) {
  const sel = cells.filter((i) => set.includes(binOf[i]));
  const m = (fn) => { let s = 0, a = 0; for (const i of sel) { s += areaCell[i] * fn(i); a += areaCell[i]; } return s / a; };
  console.log(`\n${label}\n| k | z (m) | u | v | pgf | corF | corZeta | gradK | vertical | drag | nu4 | divDamp | blMix | plume | rk4 | total | rhoK/dz at the layer's top interface (kg/m2/s) | mixing's flux at the layer's base (N/m2, + carries easterly momentum down) |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|`);
  let fluxAbove = 0;
  for (let k = 14; k < K; k++) {
    const massK = m((i) => pi[i] * dSigma[k] / g);
    fluxAbove += massK * m((i) => T.blMix[k * C + i]);
    console.log(`| ${k} | ${f(m((i) => zMid[k * C + i] - zs[i]), 0)} | ${f(m((i) => uEast[k * C + i]))} | ${f(m((i) => vNorth[k * C + i]))} | ${['pgf', 'corF', 'corZeta', 'gradK', 'vertical', 'drag', 'nu4', 'divDamp', 'blMix', 'plume', 'rk4', 'total'].map((t) => s5(m((i) => T[t][k * C + i]))).join(' | ')} | ${k > 0 ? f(m((i) => mixing[(k - 1) * C + i]), 3) : '—'} | ${f(fluxAbove, 4)} |`);
  }
}
