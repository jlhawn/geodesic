// The equatorial Pacific's temperature section from a saved state: the
// area-weighted mean over 2S–2N in 2° longitude bins from 120E to 80W at
// every 5 m from 0 to 300 m, each column's temperature interpolated
// between its layers' mid-depths (the mixed layer and every interior
// class thicker than a token), on the CPU.
//   node scripts/figures/eqsection.mjs <state.bin> <out.json>
import { writeFileSync } from 'node:fs';
import { Grid } from '../../js/grid.module.js';
import { buildMesh } from '../../js/mesh.module.js';
import { createGeography } from '../../js/geography.module.js';
import { EPS } from '../../js/ocean/layered.module.js';
import { readState, readTopography, figureHeader, DEG, rounded } from './figureState.mjs';

const [file, out] = process.argv.slice(2);
if (!file || !out) { console.error('usage: node scripts/figures/eqsection.mjs <state.bin> <out.json>'); process.exit(2); }
const s = await readState(file);
const mesh = buildMesh(new Grid(s.N)), C = mesh.nCells, L = s.ocean.h.length / C;
const geo = createGeography(mesh, readTopography());
const DEPTHS = Array.from({ length: 61 }, (_, i) => 5 * i);
const LONS = Array.from({ length: 81 }, (_, i) => 120 + 2 * i);
const sums = LONS.map(() => DEPTHS.map(() => [0, 0]));
let columns = 0;
for (let i = 0; i < C; i++) {
  if (geo.land[i] || !(s.ocean.h[i] > 0)) continue;
  if (Math.abs(mesh.latCell[i] * DEG) > 2) continue;
  let lon = mesh.lonCell[i] * DEG; if (lon < 0) lon += 360;
  const b = Math.round((lon - 120) / 2); if (b < 0 || b >= LONS.length) continue;
  const mids = [], temps = []; let z = 0;
  for (let k = 0; k < L; k++) {
    const hk = s.ocean.h[k * C + i];
    if (k > 0 && hk <= 1.1 * EPS) { z += Math.max(0, hk); continue; }
    mids.push(z + hk / 2); temps.push(s.ocean.T[k * C + i] - 273.15); z += hk;
  }
  columns++;
  for (let d = 0; d < DEPTHS.length && DEPTHS[d] <= z; d++) {
    const depth = DEPTHS[d];
    let t;
    if (depth <= mids[0]) t = temps[0];
    else if (depth >= mids[mids.length - 1]) t = temps[temps.length - 1];
    else { let j = 0; while (mids[j + 1] < depth) j++; const w = (depth - mids[j]) / (mids[j + 1] - mids[j]); t = temps[j] + w * (temps[j + 1] - temps[j]); }
    sums[b][d][0] += mesh.areaCell[i] * t; sums[b][d][1] += mesh.areaCell[i];
  }
}
const T = sums.map((col) => col.map(([q, w]) => (w > 0 ? rounded(q / w, 3) : null)));
writeFileSync(out, JSON.stringify({ ...figureHeader(file, s), layers: L, lons: LONS, depths: DEPTHS, T }));
console.log(`wrote ${out} (${columns} columns in ${LONS.length} bins, ${L} layers)`);
