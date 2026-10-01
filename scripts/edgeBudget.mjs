// The ocean's momentum budget on chosen edges of one saved state, layer by
// layer, term by term as the tendency computes it, on the CPU (1e-5 m/s2;
// 'drag+rest' is the tendency less the other terms):
//   EDGES=20091,20281 node scripts/edgeBudget.mjs <state.bin>
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { createModel, STATE_NAMES } from '../js/model.module.js';
import { decodeState, savedLevels } from '../js/stateFile.module.js';
import { gradient, curl, kineticEnergy, laplacianVelocity } from '../js/dynamics/operators.module.js';
import { THIN, PV_FLOOR, closureCoefficient } from '../js/ocean/layered.module.js';
const FILE = process.argv[2], EDGES = (process.env.EDGES ?? '').split(',').map(Number);
const OCEAN = { everySteps: 8, ...JSON.parse(process.env.OCEAN ?? '{}') };
const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
const N = saved.N, topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const model = createModel(new Grid(N), { topography, levels: savedLevels(saved), ocean: OCEAN });
const { mesh, core, state, surface, seaIce, ocean } = model;
STATE_NAMES.forEach((name, a) => state[a].set(saved[name]));
seaIce.load(state[6], saved.concentration ?? null);
model.time = saved.time;
ocean.load(saved.ocean, state[3], state[6]);
core.phaseFlux(state, 0, core.K); core.phaseColumn(state, [new Float64Array(mesh.nCells)], 0, mesh.nCells); surface.lowestWindSpeed(state[2]);
ocean.setStress(surface.stress(state, new Float64Array(mesh.nEdges)), state[6], seaIce.concentration);
const C = mesh.nCells, E = mesh.nEdges, L = ocean.layers, rho = ocean.densities, g = 9.81, rho0 = 1025, dtOcean = OCEAN.everySteps * 1350 * 16 / N;
new Float64Array(ocean.shared.params)[1] = dtOcean;
const stage = ocean.stages[0];
ocean.tendency(ocean.state, stage);
const du = stage[1], { h, u, edgeOcean } = ocean;
const hEdge = new Float64Array(ocean.shared.hEdge), pressure = new Float64Array(ocean.shared.pressure), gradEta = new Float64Array(ocean.shared.gradEta), gradRho = new Float64Array(ocean.shared.gradRho);
const { nEdgesOnEdge, edgesOnEdge, weightsOnEdge, maxEdgesOnEdge, dcEdge, dvEdge, cellsOnEdge, verticesOnEdge, fVertex } = mesh;
let spacing = 0; for (let e = 0; e < E; e++) spacing += dcEdge[e]; spacing /= E;
const nu4 = closureCoefficient(spacing, OCEAN.closureHours ?? 12), centring = OCEAN.vorticityCentring ?? 0.5, floor = OCEAN.minimumThickness ?? 50;
const zeta = new Float64Array(mesh.nVertices), Kc = new Float64Array(C), gK = new Float64Array(E), phi = new Float64Array(C), gP = new Float64Array(E), lap = new Float64Array(E), lap2 = new Float64Array(E), dS = new Float64Array(C), cS = new Float64Array(mesh.nVertices), fluxPV = new Float64Array(E), q = new Float64Array(E);
const ae = (k, e) => k * E + e, at = (k, i) => k * C + i;
console.log(`== ${FILE.split('/').pop()} day ${saved.day} step ${saved.step ?? 0}; terms in 1e-5 m/s2`);
for (let k = 0; k < L; k++) {
  const show = EDGES.filter((e) => k === 0 || hEdge[ae(k, e)] >= THIN);
  if (!show.length) continue;
  const uk = u.subarray(k * E, (k + 1) * E);
  for (let e = 0; e < E; e++) fluxPV[e] = edgeOcean[e] ? hEdge[ae(k, e)] * uk[e] : 0;
  curl(mesh, uk, zeta);
  for (let e = 0; e < E; e++) { const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1]; q[e] = 0.5 * (zeta[verticesOnEdge[2 * e]] + fVertex[verticesOnEdge[2 * e]] + zeta[verticesOnEdge[2 * e + 1]] + fVertex[verticesOnEdge[2 * e + 1]]) / Math.max(hEdge[ae(k, e)], k > 0 ? centring * 0.5 * (h[at(k, a)] + h[at(k, b)]) : 0, PV_FLOOR); }
  kineticEnergy(mesh, uk, Kc); gradient(mesh, Kc, gK);
  for (let i = 0; i < C; i++) phi[i] = k === 0 ? 0 : g * pressure[at(k, i)] / rho0;
  gradient(mesh, phi, gP);
  laplacianVelocity(mesh, uk, lap, dS, cS); laplacianVelocity(mesh, lap, lap2, dS, cS);
  for (const e of show) {
    let cor = 0; for (let s = 0; s < nEdgesOnEdge[e]; s++) { const o = edgesOnEdge[maxEdgesOnEdge * e + s]; cor += weightsOnEdge[maxEdgesOnEdge * e + s] * dvEdge[o] * fluxPV[o] * 0.5 * (q[e] + q[o]); }
    cor /= dcEdge[e];
    const eta = -g * gradEta[e], bc = k === 0 ? -g / rho0 * 0.5 * hEdge[e] * gradRho[e] : -gP[e];
    const he = Math.max(hEdge[ae(k, e)], floor), st = k === 0 ? ocean.stress[e] / rho0 / he : 0;
    const closure = -nu4 * lap2[e];
    const rest = du[ae(k, e)] - (cor - gK[e] + eta + bc + st + closure);
    console.log(`edge ${e} ${k === 0 ? 'ML' : rho[k].toFixed(2)} hEdge ${hEdge[ae(k, e)].toFixed(0)} u ${uk[e].toFixed(2)} | -g grad eta ${(eta * 1e5).toFixed(2)} baroclinic ${(bc * 1e5).toFixed(2)} Coriolis ${(cor * 1e5).toFixed(2)} -grad K ${(-gK[e] * 1e5).toFixed(2)} stress ${(st * 1e5).toFixed(2)} nu4 ${(closure * 1e5).toFixed(2)} drag+rest ${(rest * 1e5).toFixed(2)} | du/dt ${(du[ae(k, e)] * 1e5).toFixed(2)}`);
  }
}
