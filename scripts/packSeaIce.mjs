// Writes data/seaice_march.bin, the sea ice a fresh start takes at the
// March equinox (js/physics/seaIceStart.module.js reads it).
//
// Concentration: ERA5's sea_ice_cover (Hersbach et al. 2020, Q. J. R.
// Meteorol. Soc. 146, 1999-2049; read from ARCO-ERA5, public), the mean of
// 00 UTC on March 17, 20 and 23 of 1991-2020 on ERA5's 0.25-degree grid
// poleward of 35 degrees.
//
// Thickness of the covered part and the snow on it, in water equivalent:
// GIOMAS (Zhang and Rothrock 2003, Mon. Wea. Rev. 131, 845-861; Polar
// Science Center, APL/UW, public), sum(heff)/sum(area) and
// sum(snow)/sum(area) over the Marches of 1991-2020 at each sea point whose
// March concentration averages EDGE_CONCENTRATION or more; each point of a
// 1-degree grid poleward of 30 degrees takes the nearest of them in its
// hemisphere.
//
//   node scripts/packSeaIce.mjs [out.bin]
import { writeFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { surfaceRows, pool, ERA5, parseHeader } from './era5Raw.mjs';
import { encodeSeaIceStart, SEA_ICE_START_FILE, EDGE_CONCENTRATION, NO_SEA } from '../js/physics/seaIceStart.module.js';

const YEARS = Array.from({ length: 30 }, (_, n) => 1991 + n), DAYS = [17, 20, 23], POLEWARD = 35, FILLED = 30, MARCH = 2;
const GIOMAS = 'https://pscfiles.apl.washington.edu/zhang/Global_seaice';
const THICKNESS_STEP = 0.05, SNOW_STEP = 2, R = 6371;
const out = process.argv[2] ?? new URL(`../${SEA_ICE_START_FILE}`, import.meta.url).pathname;
const { NX, NY, latOfRow } = ERA5, deg = Math.PI / 180;
const unit = (lat, lon) => [Math.cos(lat * deg) * Math.cos(lon * deg), Math.cos(lat * deg) * Math.sin(lon * deg), Math.sin(lat * deg)];

const edgeRows = Math.round((90 - POLEWARD) / 0.25) + 1, bands = [[0, edgeRows], [NY - edgeRows, NY]];
const sum = new Float64Array(NX * NY), count = new Uint16Array(NX * NY);
await pool(YEARS.flatMap((year) => DAYS.flatMap((day) => bands.map(([from, to]) => async () => {
  const rows = await surfaceRows(`${year}/03/${day}`, 'sea_ice_cover', 0, from, to);
  for (let m = 0; m < rows.length; m++) if (!Number.isNaN(rows[m])) { sum[from * NX + m] += rows[m]; count[from * NX + m]++; }
}))), 8);
const cover = new Uint8Array(NX * NY);
for (const [from, to] of bands) for (let n = from * NX; n < to * NX; n++) cover[n] = count[n] ? Math.round(100 * Math.min(1, Math.max(0, sum[n] / count[n]))) : NO_SEA;

async function giomas(name, year) {
  for (let t = 0; ; t++) {
    try {
      const response = await fetch(`${GIOMAS}/${name}.H${year}.nc.gz`);
      if (!response.ok) throw new Error(`${response.status} ${name}.H${year}.nc.gz`);
      let bytes = new Uint8Array(await response.arrayBuffer());
      if (bytes[0] === 0x1f && bytes[1] === 0x8b) bytes = gunzipSync(bytes);
      return new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
    } catch (e) { if (t + 1 >= 5) throw e; await new Promise((s) => setTimeout(s, 3000 * (t + 1))); }
  }
}
function variable(view, name, record = 0) {
  const header = parseHeader(view), meta = header.vars[name];
  if (!meta) throw new Error(`no ${name} among ${Object.keys(header.vars)}`);
  const length = (dim) => header.dims.find(([d]) => d === dim)[1], recordDim = header.dims.find(([, n]) => n === 0)?.[0];
  const size = meta.dims.filter((d) => d !== recordDim).reduce((n, d) => n * length(d), 1), width = { 3: 2, 4: 4, 5: 4, 6: 8 }[meta.type];
  const records = Object.values(header.vars).filter((v) => v.dims[0] === recordDim);
  const stride = records.length === 1 ? size * width : records.reduce((n, v) => n + v.vsize, 0), at = meta.begin + (meta.dims[0] === recordDim ? record * stride : 0);
  const get = { 3: (o) => view.getInt16(o), 4: (o) => view.getInt32(o), 5: (o) => view.getFloat32(o), 6: (o) => view.getFloat64(o) }[meta.type];
  return { values: Float64Array.from({ length: size }, (_, n) => get(at + width * n)), units: meta.at.units };
}
let grid = null;
const totals = {}, units = {};
await pool(YEARS.flatMap((year) => ['heff', 'area', 'snow'].map((name) => async () => {
  const view = await giomas(name, year);
  if (variable(view, 'month', MARCH).values[0] !== 3) throw new Error(`${name}.H${year}: record ${MARCH} is not March`);
  grid ??= { lat: variable(view, 'lat_scaler').values, lon: variable(view, 'lon_scaler').values, kmt: variable(view, 'kmt').values };
  const { values, units: given } = variable(view, name, MARCH), total = (totals[name] ??= new Float64Array(values.length));
  units[name] = given;
  for (let n = 0; n < values.length; n++) if (values[n] < 1e3) total[n] += values[n];
})), 6);
const points = [];
for (let n = 0; n < grid.lat.length; n++) {
  if (!(grid.kmt[n] > 0) || totals.area[n] < EDGE_CONCENTRATION * YEARS.length) continue;
  points.push({ x: unit(grid.lat[n], grid.lon[n]), north: grid.lat[n] > 0, thickness: totals.heff[n] / totals.area[n], snow: 1000 * totals.snow[n] / totals.area[n] });
}
const iceGrid = { nLon: 360, nLat: 181, lon0: 0, dLon: 1, lat0: 90, dLat: -1 };
const thickness = new Uint8Array(iceGrid.nLon * iceGrid.nLat), snow = new Uint8Array(iceGrid.nLon * iceGrid.nLat);
for (let r = 0; r < iceGrid.nLat; r++) {
  const lat = iceGrid.lat0 + r * iceGrid.dLat;
  if (Math.abs(lat) < FILLED) continue;
  for (let c = 0; c < iceGrid.nLon; c++) {
    const x = unit(lat, iceGrid.lon0 + c * iceGrid.dLon);
    let best = null, bestDot = -2;
    for (const p of points) {
      if (p.north !== lat > 0) continue;
      const dot = x[0] * p.x[0] + x[1] * p.x[1] + x[2] * p.x[2];
      if (dot > bestDot) { bestDot = dot; best = p; }
    }
    thickness[r * iceGrid.nLon + c] = Math.max(1, Math.min(255, Math.round(best.thickness / THICKNESS_STEP)));
    snow[r * iceGrid.nLon + c] = Math.min(255, Math.round(best.snow / SNOW_STEP));
  }
}
const span = `${YEARS[0]}-${YEARS.at(-1)}`;
const source = `ERA5 sea_ice_cover (Copernicus C3S, ARCO-ERA5), mean of 00 UTC March ${DAYS.join(', ')} of ${span}, percent, poleward of ${POLEWARD} degrees; `
  + `thickness of the covered part and snow water equivalent from GIOMAS (Zhang and Rothrock 2003; APL/UW PSC) March ${span}, heff/area and snow/area at the nearest point of mean concentration ${EDGE_CONCENTRATION} or more`;
writeFileSync(out, encodeSeaIceStart({ coverGrid: { nLon: NX, nLat: NY, lon0: 0, dLon: 0.25, lat0: 90, dLat: -0.25 }, iceGrid, cover, thickness, snow, thicknessStep: THICKNESS_STEP, snowStep: SNOW_STEP, source }));

const extent = [0, 0];
for (let r = 0; r < NY; r++) {
  const cell = R * R * (0.25 * deg) ** 2 * Math.cos(latOfRow(r) * deg);
  for (let c = 0; c < NX; c++) if (cover[r * NX + c] !== NO_SEA && cover[r * NX + c] >= 100 * EDGE_CONCENTRATION) extent[latOfRow(r) > 0 ? 0 : 1] += cell;
}
const mean = (north, key) => { const chosen = points.filter((p) => p.north === north); return chosen.reduce((s, p) => s + p[key], 0) / chosen.length; };
console.log(`${out}: ERA5 extent of the mean concentration ${(extent[0] / 1e6).toFixed(2)} (north), ${(extent[1] / 1e6).toFixed(2)} Mkm² (south); GIOMAS units heff ${units.heff}, area ${units.area}, snow ${units.snow}; `
  + `${points.length} GIOMAS ice points, mean thickness ${mean(true, 'thickness').toFixed(2)} / ${mean(false, 'thickness').toFixed(2)} m, snow ${mean(true, 'snow').toFixed(0)} / ${mean(false, 'snow').toFixed(0)} kg/m² (north / south)`);
