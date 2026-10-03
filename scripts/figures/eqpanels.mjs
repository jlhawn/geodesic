// The equatorial Pacific's atmosphere and ocean from a saved state, for
// eqpanels.py: the state loaded on the GPU as a continuing segment, four
// steps, then a surface frame of the lowest layer's wind, temperature and
// the sea-level pressure and a frame of the mixed layer's current and
// temperature, with the top of the THERMOCLINE_DENSITY class from the
// saved layers (as spinup.mjs's "ocean after N days" line takes it), per
// cell within 16° of the equator from 89E to 69W, with the cells'
// polygons. OCEAN, RADIATION and the rest as for spinup.mjs; the run uses
// OCEAN='{"everySteps":8}'.
//   node scripts/figures/eqpanels.mjs <state.bin> <out.json>
import { writeFileSync } from 'node:fs';
import { THERMOCLINE_DENSITY, LAYER_DENSITIES, savedDensities } from '../../js/ocean/layered.module.js';
import { readState, figureHeader, gpuModelFrom, dt as stepOf, DEG, eastNorth, rounded } from './figureState.mjs';

const [file, out] = process.argv.slice(2);
if (!file || !out) { console.error('usage: node scripts/figures/eqpanels.mjs <state.bin> <out.json>'); process.exit(2); }
const saved = await readState(file);
const model = await gpuModelFrom(saved);
const { mesh } = model, C = mesh.nCells, dt = stepOf(saved.N);
for (let n = 0; n < 4; n++) await model.step(dt);
await model.settle();
const atm = (await model.beginFrame({ level: 'surface', fields: ['wind', 'temperature', 'mslp'] })).fields;
const sea = (await model.beginFrame({ level: 'surface', depth: 0, fields: ['currents', 'sst'] })).fields;
const classes = savedDensities(saved.ocean, C) ?? LAYER_DENSITIES;
const kT = classes.filter((r) => r < THERMOCLINE_DENSITY).length;
const rows = [], polys = [];
const { verticesOnCell, nEdgesOnCell, latVertex, lonVertex, maxEdges } = mesh;
for (let i = 0; i < C; i++) {
  const lat = mesh.latCell[i] * DEG; let lon = mesh.lonCell[i] * DEG; if (lon < 0) lon += 360;
  if (Math.abs(lat) > 16 || lon < 89 || lon > 291) continue;
  const [wu, wv] = eastNorth(mesh, atm.wind, i), land = !!model.geography.land[i];
  let top = null, cu = null, cv = null, sst = null;
  if (!land && saved.ocean.h[i] > 0) {
    top = 0; for (let k = 0; k <= kT; k++) top += saved.ocean.h[k * C + i];
    [cu, cv] = eastNorth(mesh, sea.currents, i); sst = sea.sst[i] - 273.15;
  }
  const poly = [];
  for (let k = 0; k < nEdgesOnCell[i]; k++) {
    const v = verticesOnCell[maxEdges * i + k]; let vlon = lonVertex[v] * DEG;
    if (vlon < 0) vlon += 360; if (vlon - lon > 180) vlon -= 360; if (lon - vlon > 180) vlon += 360;
    poly.push([rounded(vlon, 3), rounded(latVertex[v] * DEG, 3)]);
  }
  polys.push(poly);
  rows.push([rounded(lon, 2), rounded(lat, 2), land ? 1 : 0, rounded(wu, 2), rounded(wv, 2), rounded(atm.temperature[i] - 273.15, 2), rounded(atm.mslp[i] / 100, 2), rounded(cu, 3), rounded(cv, 3), rounded(sst, 2), rounded(top, 0)]);
}
writeFileSync(out, JSON.stringify({ ...figureHeader(file, saved), steps: 4, thermoclineDensity: THERMOCLINE_DENSITY, columns: ['lon', 'lat', 'land', 'u10', 'v10', 'tair', 'mslp', 'cu', 'cv', 'sst', 'thermocline'], rows, polys }));
console.log(`wrote ${out} (${rows.length} cells)`);
process.exit(0);
