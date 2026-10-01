// One ocean cell and its neighbours in saved states: depth, sea level, the
// layers thicker than THIN with their temperatures, and each edge's length,
// whether it is sea, and the velocity of every layer it carries.
//   CELL=6777 node scripts/cellColumn.mjs <state.bin>...
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { createModel, STATE_NAMES } from '../js/model.module.js';
import { decodeState, savedLevels } from '../js/stateFile.module.js';
import { THIN } from '../js/ocean/layered.module.js';
const CELL = Number(process.env.CELL ?? 6777);
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
let model = null, modelKey = null;
for (const file of process.argv.slice(2)) {
  const saved = await decodeState(new Uint8Array(readFileSync(file)));
  const key = `${saved.N} ${savedLevels(saved).join(',')}`;
  if (key !== modelKey) { model = createModel(new Grid(saved.N), { topography, levels: savedLevels(saved), ocean: { everySteps: 8, ...JSON.parse(process.env.OCEAN ?? '{}') } }); modelKey = key; }
  const { mesh, ocean, state } = model;
  STATE_NAMES.forEach((name, a) => state[a].set(saved[name]));
  ocean.load(saved.ocean, state[3], state[6]);
  const C = mesh.nCells, E = mesh.nEdges, L = ocean.layers, rho = ocean.densities, { h, u, Q, D, cellOcean, edgeOcean, eta } = ocean, deg = 180 / Math.PI;
  const at = (k, i) => k * C + i;
  const where = (i) => `${(mesh.latCell[i] * deg).toFixed(2)},${(mesh.lonCell[i] * deg).toFixed(2)}`;
  const layers = (i) => { const out = []; for (let k = 0; k < L; k++) if (k === 0 || h[at(k, i)] > THIN) out.push(`${k === 0 ? 'ML' : rho[k].toFixed(2)}:${h[at(k, i)].toFixed(0)}m/${(Q[at(k, i)] / h[at(k, i)] - 273.15).toFixed(1)}C`); return out.join(' '); };
  console.log(`== ${file.split('/').pop()} day ${saved.day}: cell ${CELL} at ${where(CELL)} D ${D[CELL].toFixed(0)} m, eta ${(eta[CELL] * 100).toFixed(1)} cm, area ${(mesh.areaCell[CELL] / 1e6).toFixed(0)} km2`);
  console.log(`   layers ${layers(CELL)}`);
  for (let m = 0; m < mesh.nEdgesOnCell[CELL]; m++) {
    const e = mesh.edgesOnCell[mesh.maxEdges * CELL + m], o = mesh.cellsOnCell[mesh.maxEdges * CELL + m];
    const speeds = []; for (let k = 0; k < L; k++) { const he = Math.min(k === 0 ? Infinity : h[at(k, CELL)], k === 0 ? Infinity : h[at(k, o)]); if (k === 0 || he > THIN) speeds.push(`${k === 0 ? 'ML' : rho[k].toFixed(2)}:${u[k * E + e].toFixed(2)}`); }
    console.log(`   edge ${e} (sea ${edgeOcean[e]}, dc ${(mesh.dcEdge[e] / 1e3).toFixed(0)} km, dv ${(mesh.dvEdge[e] / 1e3).toFixed(0)} km) to cell ${o} at ${where(o)} ${cellOcean[o] ? `sea D ${D[o].toFixed(0)} m eta ${(eta[o] * 100).toFixed(1)} cm` : 'LAND'}; u ${speeds.join(' ')}`);
    if (cellOcean[o]) console.log(`      ${o} layers ${layers(o)}`);
  }
}
