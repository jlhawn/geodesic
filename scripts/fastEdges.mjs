// The edges where a layer runs faster than a threshold in saved states, counted
// by layer and by 5-degree box, the TOP fastest listed with their cells and the
// layer's thickness in each; a token layer follows the one above it, so one
// fast edge counts once for each token beneath it.
//   [TOP=6] node scripts/fastEdges.mjs <threshold m/s> <state.bin>...
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { decodeState } from '../js/stateFile.module.js';
const limit = Number(process.argv[2]);
const meshes = {};
for (const file of process.argv.slice(3)) {
  const s = await decodeState(new Uint8Array(readFileSync(file)));
  const mesh = meshes[s.N] ?? (meshes[s.N] = buildMesh(new Grid(s.N)));
  const C = mesh.nCells, E = mesh.nEdges, L = s.ocean.h.length / C, deg = 180 / Math.PI;
  const rows = [], perLayer = new Map();
  for (let k = 0; k < L; k++) for (let e = 0; e < E; e++) {
    const v = s.ocean.u[k * E + e];
    if (!(Math.abs(v) > limit)) continue;
    const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1];
    perLayer.set(k, (perLayer.get(k) ?? 0) + 1);
    rows.push({ k, e, v, lat: mesh.latEdge[e] * deg, lon: Math.atan2(Math.sin(mesh.lonCell[a]) + Math.sin(mesh.lonCell[b]), Math.cos(mesh.lonCell[a]) + Math.cos(mesh.lonCell[b])) * deg, ha: s.ocean.h[k * C + a], hb: s.ocean.h[k * C + b], a, b });
  }
  rows.sort((x, y) => Math.abs(y.v) - Math.abs(x.v));
  const region = new Map();
  for (const r of rows) { const key = `${Math.round(r.lat / 5) * 5},${Math.round(r.lon / 5) * 5}`; region.set(key, (region.get(key) ?? 0) + 1); }
  console.log(`== ${file.split('/').pop()} day ${s.day}: ${rows.length} layer edges above ${limit} m/s; by layer ${[...perLayer].map(([k, n]) => `${k}:${n}`).join(' ')}; by 5-degree box (lat,lon) ${[...region].sort((x, y) => y[1] - x[1]).slice(0, 8).map(([k, n]) => `${k}:${n}`).join(' ')}`);
  for (const r of rows.slice(0, Number(process.env.TOP ?? 6))) console.log(`   k ${r.k} edge ${r.e} u ${r.v.toFixed(2)} at ${r.lat.toFixed(1)},${r.lon.toFixed(1)} cells ${r.a}/${r.b} h ${r.ha.toFixed(1)}/${r.hb.toFixed(1)}`);
}
