// The thermocline classes' ∇⁴ closure in one saved state, on the CPU, with the
// token edges' input as the ocean takes it and as candidate fills would:
//   [OCEAN='{"closureHours":12}'] [LAT=2] node scripts/closureCoupling.mjs <state.bin>
// For each treatment of the token edges: the depth integral of the interior
// classes' closure by 20-degree bin along the equator beside the stress
// (1e-5 m2/s2); per class over 180-110W, the closure on the edges within two
// rings of a token edge and on the others (1e-7 m/s2, zonal, least squares over
// the edge normals); and the rate at which the closure on a class's thick edges
// follows the mixed layer, d(closure)/d(u_ML) for 0.1 m/s eastward added to the
// mixed layer within 10 degrees of the equator with the tokens following it
// (1e-6 1/s), beside the interfacial drag's r/h.
// The treatments: 'tokens', the layer above's velocity the token edges carry;
// 'fill w', closureVelocity's least-squares fit with weight w; 'interior', the
// fit only on token edges with water denser than the class beneath both of
// their cells, so a class under the sea floor keeps its tokens; 'interior2',
// the same with the token edges next to the filled ones fitted from them;
// 'extended', the ocean's closure: interior2 with closureAdjoint.
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { createModel, STATE_NAMES } from '../js/model.module.js';
import { decodeState, savedLevels } from '../js/stateFile.module.js';
import { cellVector, laplacianVelocity } from '../js/dynamics/operators.module.js';
import { THIN, CLOSURE_RIDGE, closureCoefficient, closureVelocity, closureAdjoint } from '../js/ocean/layered.module.js';

const FILE = process.argv[2];
if (!FILE) throw new Error('usage: node scripts/closureCoupling.mjs <state.bin>');
const OCEAN = { everySteps: 8, ...JSON.parse(process.env.OCEAN ?? '{}') };
const LAT = Number(process.env.LAT ?? 2);
const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
const N = saved.N;
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const model = createModel(new Grid(N), { topography, levels: savedLevels(saved), ocean: OCEAN });
const { mesh, state, surface, seaIce, ocean, core } = model;
STATE_NAMES.forEach((name, a) => state[a].set(saved[name]));
seaIce.load(state[6], saved.concentration ?? null);
model.time = saved.time;
ocean.load(saved.ocean, state[3], state[6]);
const C = mesh.nCells, E = mesh.nEdges, L = ocean.layers, rho = ocean.densities, rho0 = 1025;
const { h, u, edgeOcean, cellOcean } = ocean;
const { dcEdge, cellsOnEdge, latCell, lonCell, latEdge, nEdgesOnCell, edgesOnCell, maxEdges, verticesOnEdge, edgesOnVertex, nEdge, tEdge, nEdgesOnEdge, edgesOnEdge, maxEdgesOnEdge } = mesh;
const deg = 180 / Math.PI;
new Float64Array(ocean.shared.params)[1] = OCEAN.everySteps * 1350 * 16 / N;
ocean.tendency(ocean.state, ocean.stages[0]);
const hEdge = new Float64Array(ocean.shared.hEdge);
let spacing = 0;
for (let e = 0; e < E; e++) spacing += dcEdge[e];
spacing /= E;
const nu4 = closureCoefficient(spacing, OCEAN.closureHours ?? 12);
core.phaseFlux(state, 0, core.K);
core.phaseColumn(state, [new Float64Array(C)], 0, C);
surface.lowestWindSpeed(state[2]);
ocean.setStress(surface.stress(state, new Float64Array(E)), state[6], seaIce.concentration);

const BINS = [];
for (let w = 160; w < 280; w += 20) BINS.push([w, w + 20]);
const binName = ([a, b]) => `${a <= 180 ? `${a}E` : `${360 - a}W`}-${b <= 180 ? `${b}E` : `${360 - b}W`}`;
const lonOf = (i) => ((lonCell[i] * deg) % 360 + 360) % 360;
const binOf = Int8Array.from({ length: C }, (_, i) => (cellOcean[i] && Math.abs(latCell[i] * deg) <= LAT ? BINS.findIndex(([a, b]) => lonOf(i) >= a && lonOf(i) < b) : -1));
const zonal = (f) => {
  const v = new Float64Array(3 * C);
  cellVector(mesh, f, v);
  return Float64Array.from({ length: C }, (_, i) => -Math.sin(lonCell[i]) * v[3 * i] + Math.cos(lonCell[i]) * v[3 * i + 1]);
};
const binMean = (x) => {
  const s = new Float64Array(BINS.length), w = new Float64Array(BINS.length);
  for (let i = 0; i < C; i++) { const b = binOf[i]; if (b < 0 || !Number.isFinite(x[i])) continue; s[b] += mesh.areaCell[i] * x[i]; w[b] += mesh.areaCell[i]; }
  return Array.from(s, (a, b) => a / w[b]);
};
const row = (label, xs, scale) => console.log(`${label.padEnd(34)}${xs.map((x) => (x * scale).toFixed(2).padStart(10)).join('')}`);
const ae = (k, e) => k * E + e;

const interior = (k, e) => {
  let a = false, b = false;
  for (let j = k + 1; j < L; j++) { a ||= h[j * C + cellsOnEdge[2 * e]] > THIN; b ||= h[j * C + cellsOnEdge[2 * e + 1]] > THIN; }
  return a && b;
};
function fit(e, valid, values) {
  const nx = nEdge[3 * e], ny = nEdge[3 * e + 1], nz = nEdge[3 * e + 2], tx = tEdge[3 * e], ty = tEdge[3 * e + 1], tz = tEdge[3 * e + 2];
  let saa = 0, sab = 0, sbb = 0, sau = 0, sbu = 0, found = false;
  for (let m = 0; m < nEdgesOnEdge[e]; m++) {
    const o = edgesOnEdge[maxEdgesOnEdge * e + m];
    if (!valid[o]) continue;
    const a = nEdge[3 * o] * nx + nEdge[3 * o + 1] * ny + nEdge[3 * o + 2] * nz, b = nEdge[3 * o] * tx + nEdge[3 * o + 1] * ty + nEdge[3 * o + 2] * tz;
    saa += a * a; sab += a * b; sbb += b * b; sau += a * values[o]; sbu += b * values[o]; found = true;
  }
  if (!found) return null;
  const p = saa + CLOSURE_RIDGE, q = sbb + CLOSURE_RIDGE;
  return (q * sau - sab * sbu) / (p * q - sab * sab);
}
const filled = new Float64Array(E), second = new Float64Array(E), valid = new Uint8Array(E);
function closureInput(k, uAll, mode) {
  const uk = uAll.subarray(k * E, (k + 1) * E);
  if (!mode.fill) return uk;
  closureVelocity(mesh, uk, hEdge.subarray(k * E, (k + 1) * E), edgeOcean, mode.fill, filled);
  if (!mode.interior) return filled;
  for (let e = 0; e < E; e++) if (filled[e] !== uk[e] && !interior(k, e)) filled[e] = uk[e];
  if (!mode.rings) return filled;
  for (let e = 0; e < E; e++) valid[e] = edgeOcean[e] && (hEdge[ae(k, e)] >= THIN || filled[e] !== uk[e]) ? 1 : 0;
  second.set(filled);
  for (let e = 0; e < E; e++) {
    if (valid[e] || !edgeOcean[e] || !interior(k, e)) continue;
    const value = fit(e, valid, filled);
    if (value !== null) second[e] = mode.fill * value + (1 - mode.fill) * uk[e];
  }
  return second;
}
const lap = new Float64Array(E), lap2 = new Float64Array(E), divS = new Float64Array(C), curlS = new Float64Array(mesh.nVertices);
const rings = { deepest: new Float64Array(ocean.shared.deepestEdge), k: 0, valid: new Uint8Array(E), second: new Float64Array(E) };
function closure(k, uAll, mode, out) {
  if (mode.extended) {
    rings.k = k;
    laplacianVelocity(mesh, closureVelocity(mesh, uAll.subarray(k * E, (k + 1) * E), hEdge.subarray(k * E, (k + 1) * E), edgeOcean, 1, filled, rings), lap, divS, curlS);
    laplacianVelocity(mesh, lap, lap2, divS, curlS);
    closureAdjoint(mesh, lap2, 1, rings);
  } else {
    laplacianVelocity(mesh, closureInput(k, uAll, mode), lap, divS, curlS);
    laplacianVelocity(mesh, lap, lap2, divS, curlS);
  }
  for (let e = 0; e < E; e++) out[e] = !edgeOcean[e] || hEdge[ae(k, e)] < THIN ? 0 : -nu4 * lap2[e];
  return out;
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
const eastOfEdge = new Float64Array(E), band = new Uint8Array(E);
for (let e = 0; e < E; e++) {
  const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
  const lon = Math.atan2(Math.sin(lonCell[a]) + Math.sin(lonCell[b]), Math.cos(lonCell[a]) + Math.cos(lonCell[b]));
  eastOfEdge[e] = -Math.sin(lon) * nEdge[3 * e] + Math.cos(lon) * nEdge[3 * e + 1];
  const lo = ((lon * deg) % 360 + 360) % 360;
  band[e] = edgeOcean[e] && Math.abs(latEdge[e] * deg) <= LAT && lo >= 180 && lo < 250 ? 1 : 0;
}
const MODES = [{ name: 'tokens' }, { name: 'fill 0.25', fill: 0.25 }, { name: 'fill 0.5', fill: 0.5 }, { name: 'fill 1', fill: 1 }, { name: 'interior', fill: 1, interior: true }, { name: 'interior2', fill: 1, interior: true, rings: 2 }, { name: 'extended', extended: true }];
const classes = [];
for (let k = 1; k < L && rho[k] < 1025.6; k++) if (rho[k] >= 1021.4) classes.push(k);

console.log(`== ${FILE.split('/').pop()}: day ${saved.day}, N=${N}, nu4 ${nu4.toExponential(2)} m^4/s, ${LAT}S-${LAT}N`);
console.log(`${''.padEnd(34)}${BINS.map((b) => binName(b).padStart(10)).join('')}`);
row('stress / rho0 (1e-5 m2/s2)', binMean(zonal(Float64Array.from(ocean.stress, (t) => t / rho0))), 1e5);
const term = new Float64Array(E);
for (const mode of MODES) {
  const column = new Float64Array(E), rows = [];
  for (let k = 1; k < L; k++) {
    closure(k, u, mode, term);
    for (let e = 0; e < E; e++) column[e] += hEdge[ae(k, e)] * term[e];
    if (!classes.includes(k)) continue;
    const token = Uint8Array.from({ length: E }, (_, e) => (hEdge[ae(k, e)] < THIN || !edgeOcean[e] ? 1 : 0)), near = ring(ring(token));
    let a = 0, w = 0, ac = 0, wc = 0;
    for (let e = 0; e < E; e++) {
      if (!band[e] || token[e]) continue;
      const x = eastOfEdge[e];
      if (near[e]) { a += term[e] * x; w += x * x; } else { ac += term[e] * x; wc += x * x; }
    }
    rows.push(`${rho[k].toFixed(2)} ${w > 0 ? (a / w * 1e7).toFixed(1) : '—'}/${wc > 0 ? (ac / wc * 1e7).toFixed(1) : '—'}`);
  }
  row(`${mode.name}: column closure (1e-5)`, binMean(zonal(column)), 1e5);
  console.log(`   near/clean 180-110W (1e-7 m/s2): ${rows.join(', ')}`);
}

const delta = 0.1, uP = Float64Array.from(u);
for (let e = 0; e < E; e++) if (edgeOcean[e] && Math.abs(latEdge[e] * deg) <= 10) uP[e] += delta * eastOfEdge[e];
for (let e = 0; e < E; e++) for (let k = 1; k < L; k++) if (hEdge[ae(k, e)] < THIN) uP[ae(k, e)] = uP[ae(k - 1, e)];
const r = OCEAN.interfacialDrag ?? 2e-4;
console.log(`\nd(closure)/d(u_ML) on the thick edges of each class, 2S-2N 180-110W (1e-6 1/s); the interfacial drag's r/h is ${(r / 50 * 1e6).toFixed(1)} at 50 m, ${(r / 20 * 1e6).toFixed(1)} at 20 m`);
const base = new Float64Array(E), pert = new Float64Array(E);
for (const mode of MODES) {
  const rows = [];
  for (const k of classes) {
    closure(k, u, mode, base);
    closure(k, uP, mode, pert);
    let a = 0, w = 0;
    for (let e = 0; e < E; e++) { if (!band[e] || hEdge[ae(k, e)] < THIN) continue; const x = eastOfEdge[e]; a += hEdge[ae(k, e)] * (pert[e] - base[e]) * x; w += hEdge[ae(k, e)] * delta * x * x; }
    if (w > 0) rows.push(`${rho[k].toFixed(2)} ${(a / w * 1e6).toFixed(1)}`);
  }
  console.log(`${mode.name.padEnd(10)} ${rows.join(', ')}`);
}
