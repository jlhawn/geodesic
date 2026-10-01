// The grid-scale share of each class's thickness-flux divergence in saved
// states, the measure of a checkerboard: over cells within LAT degrees of the
// equator where the class is thicker than THIN, the area-weighted mean square
// of the cell's div(h u) minus the mean over its neighbours holding the class,
// over its variance (white noise about 1.17, a smooth field near 0), with the
// rms divergence, for the mixed layer and three groups of classes, split by
// the cell's edges: no token edge ('clean'), token edges with a denser class
// beneath both cells, which closureVelocity fits, or any token edge over the
// sea floor.
//   [LAT=10] node scripts/divergenceNoise.mjs <state.bin>...
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { topographyFromInt16, createGeography } from '../js/geography.module.js';
import { decodeState } from '../js/stateFile.module.js';
import { divergence } from '../js/dynamics/operators.module.js';
import { createOcean, THIN } from '../js/ocean/layered.module.js';

const files = process.argv.slice(2);
if (!files.length) throw new Error('usage: node scripts/divergenceNoise.mjs <state.bin>...');
const LAT = Number(process.env.LAT ?? 10);
let cache = null;
for (const FILE of files) {
  const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
  const N = saved.N;
  if (!cache || cache.N !== N) {
    const mesh = buildMesh(new Grid(N));
    const geography = createGeography(mesh, topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer));
    cache = { N, mesh, geography };
  }
  const { mesh, geography } = cache, C = mesh.nCells, E = mesh.nEdges;
  const ocean = createOcean(mesh, { geography, everySteps: 8 });
  ocean.load(saved.ocean, Float64Array.from(saved.surfaceT), Float64Array.from(saved.ice));
  new Float64Array(ocean.shared.params)[1] = 8 * 1350 * 16 / N;
  ocean.tendency(ocean.state, ocean.stages[0]);
  const L = ocean.layers, rho = ocean.densities, { h, u, edgeOcean, cellOcean } = ocean;
  const hEdge = new Float64Array(ocean.shared.hEdge), deepest = new Float64Array(ocean.shared.deepestEdge);
  const { nEdgesOnCell, edgesOnCell, cellsOnCell, maxEdges, areaCell, latCell } = mesh;
  const deg = 180 / Math.PI;
  const groups = { ML: [0], 'thermocline 1022-1024.75': [], '1025-1026.5': [], 'deep >1026.5': [] };
  for (let k = 1; k < L; k++) {
    if (rho[k] >= 1022 && rho[k] <= 1024.75 + 1e-9) groups['thermocline 1022-1024.75'].push(k);
    else if (rho[k] > 1024.75 + 1e-9 && rho[k] <= 1026.5 + 1e-9) groups['1025-1026.5'].push(k);
    else if (rho[k] > 1026.5 + 1e-9) groups['deep >1026.5'].push(k);
  }
  const F = new Float64Array(E), div = new Float64Array(C);
  console.log(`== ${FILE.split('/').pop()} N=${N} day ${saved.day}; |lat| <= ${LAT}; grid-scale share of div(h u) (cell minus mean of its thick neighbours, over the variance; white noise ~1.17), rms div in 1e-6 m/s, cells`);
  for (const [name, ks] of Object.entries(groups)) {
    const acc = ['clean', 'fitted', 'excluded'].map(() => ({ r: 0, v: 0, s: 0, a: 0, n: 0, m: 0 }));
    const pass = [];
    for (const k of ks) {
      for (let e = 0; e < E; e++) F[e] = edgeOcean[e] ? hEdge[k * E + e] * u[k * E + e] : 0;
      divergence(mesh, F, div);
      const thick = (i) => cellOcean[i] && h[k * C + i] > THIN;
      for (let i = 0; i < C; i++) {
        if (!thick(i) || Math.abs(latCell[i] * deg) > LAT) continue;
        let kind = 0, nb = 0, sum = 0;
        for (let m = 0; m < nEdgesOnCell[i]; m++) {
          const e = edgesOnCell[maxEdges * i + m], j = cellsOnCell[maxEdges * i + m];
          if (k > 0 && edgeOcean[e] && hEdge[k * E + e] < THIN) kind = Math.max(kind, deepest[e] > k ? 1 : 2);
          if (thick(j)) { nb++; sum += div[j]; }
        }
        if (!nb) continue;
        pass.push([kind, i, div[i], sum / nb]);
        acc[kind].m += areaCell[i] * div[i]; acc[kind].a += areaCell[i];
      }
    }
    for (const a of acc) a.m /= a.a || 1;
    for (const [kind, i, d, nm] of pass) { const a = acc[kind]; a.v += areaCell[i] * (d - a.m) ** 2; a.r += areaCell[i] * (d - nm) ** 2; a.s += areaCell[i] * d * d; a.n++; }
    const show = (a) => (a.n ? `${(a.r / a.v).toFixed(3)} rms ${(1e6 * Math.sqrt(a.s / a.a)).toFixed(2)} n ${a.n}` : '—');
    console.log(`  ${name.padEnd(26)} clean ${show(acc[0])} | beside fitted tokens ${show(acc[1])} | beside floor tokens ${show(acc[2])}`);
  }
}
