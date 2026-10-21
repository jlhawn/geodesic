// The per-cell fields of the four-panel ocean figure (oceanPanels.py)
// from a saved state on the CPU: the sea surface temperature at the
// state's instant, the SST and WOA23's monthly climatology that
// scripts/woaSstBias.mjs compares (the state's window mean over the
// window's days when it carries one, as `sstMean` and `sstWindow`, else
// the instant on the state's date), the current averaged over the top
// 15 m (east and north components) and the mixed-layer depth (where the
// density first exceeds the surface layer's by 0.03 kg/m³), for a slab
// or a hybrid ocean.
//   node scripts/figures/oceanPanels.mjs <state.bin> <out.json>
import { readFileSync, writeFileSync } from 'node:fs';
import { createGeography } from '../../js/geography.module.js';
import { cellVector } from '../../js/dynamics/operators.module.js';
import { THIN } from '../../js/ocean/layered.module.js';
import { layerDensities, mixedLayerDepth } from '../../js/ocean/column.module.js';
import { decodeMonthlySst, climatologyAt, MONTHLY_SST_FILE } from '../../js/ocean/climatology.module.js';
import { stateSst } from '../woaSstBias.mjs';
import { readState, readGeography, figureHeader, DEG, eastNorth, rounded } from './figureState.mjs';
import { stateGrid } from '../../js/stateGrid.module.js';

const [file, out] = process.argv.slice(2);
if (!file || !out) { console.error('usage: node scripts/figures/oceanPanels.mjs <state.bin> <out.json>'); process.exit(2); }
const KELVIN = 273.15, SURFACE_DEPTH = 15;
const saved = await readState(file);
const { topography, geography } = readGeography(saved);
const mesh = stateGrid(saved, topography).mesh, C = mesh.nCells, E = mesh.nEdges;
const geo = createGeography(mesh, topography, geography);
const woa = decodeMonthlySst(readFileSync(new URL(`../../${MONTHLY_SST_FILE}`, import.meta.url)));
const judged = stateSst(saved, woa), { months } = judged, mean = judged.kind === 'mean';
const { ocean } = saved, hybrid = ocean.vertical === 'hybrid', L = ocean.densities.length, { h, u, T, S } = ocean;
const at = (k, i) => k * C + i;
const visible = (k, i) => hybrid || k === 0 || h[at(k, i)] > THIN;
const velocity = [];
for (let k = 0; k < L; k++) velocity.push(cellVector(mesh, u.subarray(k * E, (k + 1) * E)));
const surfaceCurrent = (i) => {
  let top = 0, ue = 0, vn = 0, w = 0;
  for (let k = 0; k < L && top < SURFACE_DEPTH; k++) {
    const hk = h[at(k, i)];
    if (visible(k, i)) {
      const part = Math.min(hk, SURFACE_DEPTH - top), [e, n] = eastNorth(mesh, velocity[k], i);
      ue += part * e; vn += part * n; w += part;
    }
    top += hk;
  }
  return w > 0 ? [ue / w, vn / w] : [null, null];
};
const hc = new Float64Array(L), Tc = new Float64Array(L), Sc = new Float64Array(L), rhoc = new Float64Array(L);
const mixedDepth = (i) => {
  let K = 0;
  for (let k = 0; k < L; k++) {
    if (!visible(k, i) || h[at(k, i)] <= 0.01) continue;
    hc[K] = h[at(k, i)]; Tc[K] = T[at(k, i)]; Sc[K] = S[at(k, i)]; K++;
  }
  return K ? mixedLayerDepth(hc, layerDensities(Tc, Sc, K, rhoc), K) : null;
};
const columns = { lat: [], lon: [], land: [], ice: [], sst: [], ...(mean ? { sstMean: [] } : {}), woa: [], ue: [], vn: [], mld: [] };
for (let i = 0; i < C; i++) {
  let column = 0; for (let k = 0; k < L; k++) column += h[at(k, i)];
  const sea = !geo.land[i] && column > 1;
  columns.lat.push(rounded(mesh.latCell[i] * DEG, 2));
  columns.lon.push(rounded(mesh.lonCell[i] * DEG, 2));
  columns.land.push(sea ? 0 : 1);
  columns.ice.push(sea ? rounded(saved.ice[i], 2) : 0);
  const clim = sea ? climatologyAt(woa, months, mesh.latCell[i], mesh.lonCell[i]) : NaN;
  columns.sst.push(sea ? rounded(T[i] - KELVIN, 2) : null);
  if (mean) columns.sstMean.push(sea ? rounded(judged.model[i] - KELVIN, 2) : null);
  columns.woa.push(rounded(clim - KELVIN, 2));
  const [ue, vn] = sea ? surfaceCurrent(i) : [null, null];
  columns.ue.push(rounded(ue, 3)); columns.vn.push(rounded(vn, 3));
  columns.mld.push(sea ? rounded(mixedDepth(i), 1) : null);
}
const sstWindow = mean ? { kind: 'mean', days: [judged.days[0], judged.days[judged.days.length - 1]], length: judged.days.length } : { kind: 'instant', day: judged.day };
writeFileSync(out, JSON.stringify({ ...figureHeader(file, saved), vertical: hybrid ? 'hybrid' : 'slab', surfaceDepth: SURFACE_DEPTH, sstWindow, ...columns }));
console.log(`${out}: ${C} cells, ${columns.sst.filter((v) => v !== null).length} at sea; ${mean ? `the ${sstWindow.length}-day mean SST of days ${sstWindow.days.join('–')}` : 'the instant SST'} against WOA23`);
