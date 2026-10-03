// The per-cell fields of the four-panel state map (stateMaps.py): the
// mixed layer's temperature (sea), the sea ice's thickness, concentration
// and snow, the land's vegetation cover and trees, and the surface
// temperature, from a saved state on the CPU.
//   node scripts/figures/stateMaps.mjs <state.bin> <out.json>
import { writeFileSync } from 'node:fs';
import { Grid } from '../../js/grid.module.js';
import { buildMesh } from '../../js/mesh.module.js';
import { createGeography } from '../../js/geography.module.js';
import { readState, readTopography, figureHeader, DEG, rounded } from './figureState.mjs';

const [file, out] = process.argv.slice(2);
if (!file || !out) { console.error('usage: node scripts/figures/stateMaps.mjs <state.bin> <out.json>'); process.exit(2); }
const saved = await readState(file);
const mesh = buildMesh(new Grid(saved.N)), C = mesh.nCells;
const geo = createGeography(mesh, readTopography());
const { ocean, land } = saved, zero = new Float32Array(C);
const vegetation = land.vegetation ?? zero, trees = land.canopy ?? zero, snow = land.snow ?? zero, concentration = saved.concentration ?? zero;
const columns = { lat: [], lon: [], land: [], iceSheet: [], ice: [], concentration: [], snow: [], sst: [], ts: [], vegetation: [], trees: [] };
for (let i = 0; i < C; i++) {
  const onLand = !!geo.land[i], sea = !onLand && ocean.h[i] > 0;
  columns.lat.push(rounded(mesh.latCell[i] * DEG, 2));
  columns.lon.push(rounded(mesh.lonCell[i] * DEG, 2));
  columns.land.push(onLand ? 1 : 0);
  columns.iceSheet.push(geo.iceSheet[i] ? 1 : 0);
  columns.ice.push(onLand ? 0 : rounded(saved.ice[i], 2));
  columns.concentration.push(onLand ? 0 : rounded(saved.ice[i] > 0 ? (concentration[i] > 0 ? concentration[i] : 1) : 0, 3));
  columns.snow.push(rounded(snow[i], 1));
  columns.sst.push(sea ? rounded(ocean.T[i] - 273.15, 2) : null);
  columns.ts.push(rounded(saved.surfaceT[i] - 273.15, 1));
  columns.vegetation.push(onLand ? rounded(vegetation[i], 3) : 0);
  columns.trees.push(onLand ? rounded(trees[i], 3) : 0);
}
writeFileSync(out, JSON.stringify({ ...figureHeader(file, saved), ...columns }));
console.log(`wrote ${out} (${C} cells)`);
