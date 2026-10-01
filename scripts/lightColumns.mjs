// Ocean columns filled with light water to depth in saved states: cells within
// LAT degrees of the equator, deeper than DEPTH metres, whose mixed layer and
// classes lighter than DENSITY hold more than SHARE of the column.
//   [LAT=25] [DEPTH=300] [DENSITY=1025] [SHARE=0.7] node scripts/lightColumns.mjs <state.bin>...
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { decodeState } from '../js/stateFile.module.js';

const LAT = Number(process.env.LAT ?? 25), DEPTH = Number(process.env.DEPTH ?? 300), DENSITY = Number(process.env.DENSITY ?? 1025), SHARE = Number(process.env.SHARE ?? 0.7);
const meshes = {}, deg = 180 / Math.PI;
for (const file of process.argv.slice(2)) {
  const s = await decodeState(new Uint8Array(readFileSync(file)));
  const mesh = meshes[s.N] ?? (meshes[s.N] = buildMesh(new Grid(s.N)));
  const C = mesh.nCells, L = s.ocean.h.length / C, rho = s.ocean.densities;
  const found = [];
  for (let i = 0; i < C; i++) {
    if (Math.abs(mesh.latCell[i] * deg) >= LAT) continue;
    let total = 0, light = 0;
    for (let k = 0; k < L; k++) { const hk = s.ocean.h[k * C + i]; total += hk; if (k === 0 || rho[k - 1] < DENSITY) light += hk; }
    if (total > DEPTH && light > SHARE * total) found.push(`${i}@${(mesh.latCell[i] * deg).toFixed(1)},${(mesh.lonCell[i] * deg).toFixed(1)} ${light.toFixed(0)}/${total.toFixed(0)} m`);
  }
  console.log(`${file.split('/').pop()}, day ${s.day}: ${found.length} columns; ${found.join('; ')}`);
}
