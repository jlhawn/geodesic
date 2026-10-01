// The equatorial undercurrent's latitude band in saved states: the class with
// the strongest thickness-weighted eastward flow over 180-100W within 2S-2N
// (js/audit.module.js equatorialOcean), its eastward flow and mean depth by
// degree of latitude from 8S to 8N over the same longitudes, and the band
// around the equator where that flow exceeds half its peak.
//   node scripts/undercurrentBand.mjs <state.bin>...
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { decodeState } from '../js/stateFile.module.js';
import { cellVector } from '../js/dynamics/operators.module.js';
import { equatorialOcean, inLongitudes } from '../js/audit.module.js';

const meshes = {}, deg = 180 / Math.PI;
for (const file of process.argv.slice(2)) {
  const s = await decodeState(new Uint8Array(readFileSync(file)));
  const mesh = meshes[s.N] ?? (meshes[s.N] = buildMesh(new Grid(s.N)));
  const C = mesh.nCells, E = mesh.nEdges, ocean = s.ocean, densities = Array.from(ocean.densities);
  const land = Uint8Array.from({ length: C }, (_, i) => (ocean.h[i] > 0 ? 0 : 1));
  const o = equatorialOcean(mesh, land, { h: ocean.h, u: ocean.u, densities }, new Float64Array(E));
  const k = densities.indexOf(o.undercurrentClass) + 1;
  const vector = cellVector(mesh, Float64Array.from(ocean.u.slice(k * E, (k + 1) * E)));
  const rows = new Map();
  for (let i = 0; i < C; i++) {
    const lat = mesh.latCell[i] * deg, lon = mesh.lonCell[i] * deg;
    if (land[i] || Math.abs(lat) > 8.5 || !inLongitudes(lon, -180, -100)) continue;
    const row = Math.round(lat), a = mesh.areaCell[i], thick = ocean.h[k * C + i];
    let above = 0;
    for (let j = 0; j < k; j++) above += ocean.h[j * C + i];
    const east = -Math.sin(mesh.lonCell[i]) * vector[3 * i] + Math.cos(mesh.lonCell[i]) * vector[3 * i + 1];
    const r = rows.get(row) ?? { mass: 0, flow: 0, depth: 0 };
    r.mass += a * thick; r.flow += a * thick * east; r.depth += a * thick * (above + thick / 2);
    rows.set(row, r);
  }
  const profile = [...rows].sort((x, y) => x[0] - y[0]).map(([lat, r]) => ({ lat, u: r.mass > 0 ? r.flow / r.mass : NaN, depth: r.mass > 0 ? r.depth / r.mass : NaN }));
  const peak = profile.filter((p) => Math.abs(p.lat) <= 2).reduce((best, p) => (p.u > best.u ? p : best), { u: -Infinity });
  let south = peak.lat, north = peak.lat;
  const at = (lat) => profile.find((p) => p.lat === lat);
  while (at(south - 1) && at(south - 1).u > 0.5 * peak.u) south--;
  while (at(north + 1) && at(north + 1).u > 0.5 * peak.u) north++;
  console.log(`${file.split('/').pop()} day ${s.day}: class ${o.undercurrentClass} ${o.undercurrent.toFixed(2)} m/s at ${o.undercurrentDepth.toFixed(0)} m (2S-2N, 180-100W); peak ${peak.u.toFixed(2)} m/s at ${peak.lat}°, above half of it ${south}° to ${north}°; by latitude ${profile.map((p) => `${p.lat}:${p.u.toFixed(2)}@${p.depth.toFixed(0)}`).join(' ')}`);
}
