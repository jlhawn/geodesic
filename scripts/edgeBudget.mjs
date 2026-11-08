// The ocean's momentum budget on chosen edges of one saved state, layer by
// layer, term by term as the tendency computes it, on the CPU (1e-5 m/s2;
// 'rest' is the tendency less the printed terms). A row the tendency pins
// shows only its relaxation toward the layer above, which replaces its
// forces; below the top layer, rows thinner than twice minimumLayer are
// pinned and left out:
//   EDGES=20091,20281 node scripts/edgeBudget.mjs <state.bin>
import { readFileSync } from 'node:fs';
import { stateGeography } from '../js/geography.module.js';
import { createModel, STATE_NAMES } from '../js/model.module.js';
import { decodeState, savedLevels, stateTerrain } from '../js/stateFile.module.js';
import { gradient, curl, kineticEnergy, laplacianVelocity } from '../js/dynamics/operators.module.js';
import { PV_FLOOR, closureCoefficient, oceanCells, inertialStepLimit } from '../js/ocean/layered.module.js';
import { H_MIN } from '../js/ocean/column.module.js';
import { withCadence } from '../js/cadence.module.js';
import { stateGrid } from '../js/stateGrid.module.js';
import { diagnoseExchange } from './surfaceExchange.mjs';
const FILE = process.argv[2], EDGES = (process.env.EDGES ?? '').split(',').map(Number), OCEAN = JSON.parse(process.env.OCEAN ?? '{}');
const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
const N = saved.N, dt = 1350 * 16 / N, ground = stateGeography(saved), { topography } = ground;
const model = createModel(stateGrid(saved, topography).mesh, { terrain: stateTerrain(saved), topography, geography: ground.geography, land: ground.land, levels: savedLevels(saved), ocean: withCadence({ ocean: OCEAN }, dt, N).ocean });
const { mesh, core, state, surface, seaIce, ocean } = model;
STATE_NAMES.forEach((name, a) => state[a].set(saved[name]));
seaIce.load(state[6], saved.concentration ?? null);
model.time = saved.time;
ocean.load(saved.ocean, state[3], state[6]);
core.phaseFlux(state, 0, core.K); core.phaseColumn(state, [new Float64Array(mesh.nCells)], 0, mesh.nCells); diagnoseExchange(model);
ocean.setStress(surface.stress(state, new Float64Array(mesh.nEdges)), state[6], seaIce.concentration);
const C = mesh.nCells, E = mesh.nEdges, L = ocean.layers, rho = ocean.densities, g = OCEAN.gravity ?? 9.81, rho0 = OCEAN.density ?? 1025;
const dtOcean = ocean.everySteps * dt, part = dtOcean / Math.max(1, Math.ceil(dtOcean / inertialStepLimit(mesh))), rate = Math.min(1 / 3600, 1 / part);
const params = new Float64Array(ocean.shared.params);
params[0] = rate;
params[1] = part;
const stage = ocean.stages[0];
ocean.tendency(ocean.state, stage);
const du = stage[1], { h, u, stress, edgeOcean } = ocean, shared = (name) => new Float64Array(ocean.shared[name]);
const hEdge = shared('hEdge'), pressure = shared('pressure'), gradEta = shared('gradEta'), layerRho = shared('layerRho'), columnAbove = shared('columnAbove'), edgeAbove = shared('edgeAbove'), edgeBottom = shared('edgeBottom');
const { nEdgesOnEdge, edgesOnEdge, weightsOnEdge, maxEdgesOnEdge, dcEdge, dvEdge, cellsOnEdge, verticesOnEdge, fVertex } = mesh;
let spacing = 0; for (let e = 0; e < E; e++) spacing += dcEdge[e]; spacing /= E;
const nu4 = closureCoefficient(spacing, OCEAN.closureHours ?? 12, OCEAN.closureSpacing), centring = OCEAN.vorticityCentring ?? 0.5, minimumLayer = OCEAN.minimumLayer ?? H_MIN;
const pvFloor = Float64Array.from(ocean.table.minimum, (m) => Math.min(PV_FLOOR, m)), { shelf, draft } = oceanCells(C, model.geography);
const topDrag = OCEAN.topDrag ?? 2.5e-3, bottomDrag = OCEAN.bottomDrag ?? 3e-3, momentumDepth = OCEAN.momentumDepth ?? 50;
const thinEdge = 2 * minimumLayer, tokenEdge = Float64Array.from(ocean.table.minimum, (m) => 2 * Math.min(m, ocean.table.spacing ?? m));
const zeta = new Float64Array(mesh.nVertices), Kc = new Float64Array(C), gK = new Float64Array(E), lap = new Float64Array(E), lap2 = new Float64Array(E), dS = new Float64Array(C), cS = new Float64Array(mesh.nVertices), fluxPV = new Float64Array(E), q = new Float64Array(E);
const ae = (k, e) => k * E + e, at = (k, i) => k * C + i, f = (x) => (x * 1e5).toFixed(2), fine = (x) => (x * 1e5).toFixed(4);
console.log(`== ${FILE.split('/').pop()} day ${saved.day} step ${saved.step ?? 0}; terms in 1e-5 m/s2; pinned rows relax at 1/${(1 / rate).toFixed(0)} s`);
for (const e of EDGES) if (!edgeOcean[e]) console.log(`edge ${e} is not an ocean edge: its tendency is zero`);
for (let k = 0; k < L; k++) {
  const rows = EDGES.filter((e) => edgeOcean[e] && (k === 0 || hEdge[ae(k, e)] > thinEdge));
  if (!rows.length) continue;
  const uk = u.subarray(k * E, (k + 1) * E);
  for (let e = 0; e < E; e++) fluxPV[e] = edgeOcean[e] ? hEdge[ae(k, e)] * uk[e] : 0;
  curl(mesh, uk, zeta);
  for (let e = 0; e < E; e++) { const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1]; q[e] = 0.5 * (zeta[verticesOnEdge[2 * e]] + fVertex[verticesOnEdge[2 * e]] + zeta[verticesOnEdge[2 * e + 1]] + fVertex[verticesOnEdge[2 * e + 1]]) / Math.max(hEdge[ae(k, e)], centring * 0.5 * (h[at(k, a)] + h[at(k, b)]), pvFloor[k]); }
  kineticEnergy(mesh, uk, Kc); gradient(mesh, Kc, gK);
  laplacianVelocity(mesh, uk, lap, dS, cS); laplacianVelocity(mesh, lap, lap2, dS, cS);
  for (const e of rows) {
    const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1], na = at(k, a), nb = at(k, b), Ha = columnAbove[na], Hb = columnAbove[nb], hk = hEdge[ae(k, e)], tendency = du[ae(k, e)];
    const head = `edge ${e} layer ${k} ${rho[k].toFixed(2)} hEdge ${hk.toFixed(1)} u ${uk[e].toFixed(2)}`;
    const ra = Ha - draft[a], rb = Hb - draft[b];
    if (k > 0 && hk <= tokenEdge[k] && !(Math.min(ra + h[na], rb + h[nb]) > Math.max(ra, rb))) {
      const relaxation = (u[ae(k - 1, e)] - uk[e]) * rate;
      console.log(`${head} | pinned, relaxing toward layer ${k - 1} ${f(relaxation)} rest ${fine(tendency - relaxation)} | du/dt ${f(tendency)}`);
      continue;
    }
    let cor = 0; for (let s = 0; s < nEdgesOnEdge[e]; s++) { const o = edgesOnEdge[maxEdgesOnEdge * e + s]; cor += weightsOnEdge[maxEdgesOnEdge * e + s] * dvEdge[o] * fluxPV[o] * 0.5 * (q[e] + q[o]); }
    cor /= dcEdge[e];
    const top = Math.max(Ha, Hb), bottom = Math.min(Ha + h[na], Hb + h[nb]), zR = bottom > top ? 0.5 * (top + bottom) : 0.5 * (Ha + 0.5 * h[na] + Hb + 0.5 * h[nb]);
    const eta = -g * gradEta[e], bc = -g / rho0 * (pressure[nb] + layerRho[nb] * (zR - Hb) - pressure[na] - layerRho[na] * (zR - Ha)) / dcEdge[e];
    const floored = Math.max(hk, minimumLayer), st = k === 0 ? stress[e] / rho0 / floored : 0;
    const shelfDrag = k === 0 ? -topDrag * 0.5 * (shelf[a] + shelf[b]) * Math.abs(uk[e]) * uk[e] / floored : 0;
    const layerTop = edgeAbove[ae(k, e)], floor = edgeBottom[e], overlap = Math.min(layerTop + hk, floor) - Math.max(layerTop, floor - momentumDepth);
    const drag = overlap > 0 ? -bottomDrag * Math.abs(uk[e]) * uk[e] * (overlap / momentumDepth) / floored : 0;
    const closure = -nu4 * lap2[e];
    const rest = tendency - (cor - gK[e] + eta + bc + st + shelfDrag + drag + closure);
    console.log(`${head} | -g grad eta ${f(eta)} baroclinic ${f(bc)} Coriolis ${f(cor)} -grad K ${f(-gK[e])} stress ${f(st)} shelf drag ${f(shelfDrag)} bottom drag ${f(drag)} nu4 ${f(closure)} rest ${fine(rest)} | du/dt ${f(tendency)}`);
  }
}
