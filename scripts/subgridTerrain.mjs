#!/usr/bin/env node
/*
 * Writes data/subgrid_N<N>.bin, each land cell's subgrid orography from
 * GMTED2010's 30″ mean elevation as the IFS defines its fields (Cy47r3
 * Part IV §11.3.3–11.3.4), from the 2′30″ grids of
 * scripts/subgridTerrain.py's `filter` in the download cache:
 *  - μ, γ, θ, σ (geography.module.js's subgridOrography) of h5, the 30″
 *    orography smoothed at 5 km, less the orography the model resolves
 *    (the surface geopotential of the 0.25° raster on this mesh over g),
 *    so the fields hold the scales between 5 km and what the mesh carries;
 *  - σ_flt, the square root of the cell's area-weighted mean of
 *    (h₂ − h₂₀)², the 3–22 km band the form drag takes.
 * Every 2′30″ point counts for the cell nearest it; sea cells hold zeros.
 *
 *   node scripts/subgridTerrain.mjs [CACHE] [N ...]
 */
import { readFileSync, writeFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { topographyFromInt16, subgridOrography, encodeSubgrid, decodeSubgrid, subgridUrl } from '../js/geography.module.js';
import { createModel } from '../js/model.module.js';

export const FINE_ROWS = 4320, FINE_COLS = 8640;

export function readFine(cache, name, rows = FINE_ROWS, cols = FINE_COLS) {
  const bytes = readFileSync(`${cache}/fine_${name}.f32`);
  if (bytes.byteLength !== 4 * rows * cols) throw new Error(`fine_${name}.f32: ${bytes.byteLength} bytes`);
  return new Float32Array(bytes.buffer, bytes.byteOffset, rows * cols);
}

/*
 * The fields for a mesh from the fine grids `h5` and `flt2` (rows from
 * the north, columns from −180°) and the model of that mesh (its mesh,
 * land mask and surface geopotential).
 */
export function meshFields(model, h5, flt2, rows = FINE_ROWS, cols = FINE_COLS) {
  const g = model.core.diagnostics.g, phis = model.surfaceGeopotential;
  return subgridOrography(model.mesh, { rows, cols, data: h5 }, phis ? Float64Array.from(phis, (p) => p / g) : null, model.geography.land, { filtered: { data: flt2 } });
}

function summary(N, model, fields) {
  const { land, iceSheet } = model.geography, { areaCell } = model.mesh;
  const mean = (name, keep = () => true) => {
    let w = 0, s = 0;
    for (let i = 0; i < land.length; i++) if (land[i] && keep(i)) { w += areaCell[i]; s += areaCell[i] * fields[name][i]; }
    return s / w;
  };
  const bins = [50, 100, 200, 400, Infinity];
  const share = bins.map(() => 0);
  let total = 0, points = 0, cells = 0;
  for (let i = 0; i < land.length; i++) {
    if (!land[i]) continue;
    total += areaCell[i]; points += fields.count[i]; cells++;
    share[bins.findIndex((b) => fields.deviation[i] < b)] += areaCell[i];
  }
  const fmt = (x, d = 3) => x.toFixed(d);
  console.log(`N=${N}: ${cells} land cells, ${fmt(points / cells, 1)} fine points a cell; land-mean μ ${fmt(mean('deviation'), 1)} m, σ ${fmt(mean('slope'), 4)}, γ ${fmt(mean('anisotropy'))}, σ_flt ${fmt(mean('filtered'), 1)} m (off the ice sheets ${fmt(mean('filtered', (i) => !iceSheet[i]), 1)} m)`);
  console.log(`  μ < 50 / 50–100 / 100–200 / 200–400 / ≥ 400 m: ${share.map((s) => fmt(s / total, 2)).join(' / ')} of the land`);
}

if (import.meta.url === `file://${process.argv[1]}`) {
  const cache = process.argv[2] ?? new URL('../../geodesic-terrain-cache', import.meta.url).pathname;
  const Ns = process.argv.length > 3 ? process.argv.slice(3).map(Number) : [16, 32, 64, 128];
  const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
  const h5 = readFine(cache, 'h5'), flt2 = readFine(cache, 'flt2');
  for (const N of Ns) {
    const start = performance.now();
    const model = createModel(new Grid(N), { physics: false, topography });
    const fields = meshFields(model, h5, flt2);
    const bytes = encodeSubgrid(fields);
    writeFileSync(subgridUrl(N), new Uint8Array(bytes));
    const back = decodeSubgrid(bytes);
    const worst = Object.fromEntries([['deviation', 1], ['anisotropy', 0.01], ['orientation', 0.01], ['slope', 1e-3], ['filtered', 1]].map(([name, floor]) => {
      let w = 0;
      for (let i = 0; i < back[name].length; i++) if (model.geography.land[i]) w = Math.max(w, Math.abs(back[name][i] - fields[name][i]) / Math.max(floor, Math.abs(fields[name][i])));
      return [name, w.toExponential(1)];
    }));
    summary(N, model, fields);
    console.log(`  ${subgridUrl(N).pathname.split('/').slice(-2).join('/')}: ${bytes.byteLength} bytes, ${((performance.now() - start) / 1000).toFixed(1)} s, quantization's worst relative error ${JSON.stringify(worst)} (floors 1 m, 0.01, 0.01 rad, 10⁻³, 1 m)`);
  }
}
