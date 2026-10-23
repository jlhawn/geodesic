/*
 * The sea ice a fresh start takes at the March equinox, day 0
 * (data/seaice_march.bin, written by scripts/packSeaIce.mjs, which names
 * the sources): the observed concentration on one regular
 * latitude–longitude grid, and on another the thickness of the covered
 * part and the snow on it in water equivalent. `at(lat, lon)`, in
 * radians, gives { concentration, thickness, snow } at a point, thickness
 * in m and snow in kg/m² per unit ice area: the concentration bilinear
 * between the four surrounding grid points that hold sea, or that of the
 * nearest within `reach` degrees when none does, none below
 * EDGE_CONCENTRATION, the ice edge of the passive-microwave records; the
 * thickness and snow bilinear on their grid.
 *
 * Layout, little-endian:
 *   0   char[4]  'SIC1'
 *   4   uint32   byte offset of the data
 *   8   uint16   nLon, nLat of the concentration grid, then of the thickness grid
 *   16  float32  lon0, dLon, lat0, dLat of each (degrees; centre of column 0 / row 0),
 *                then the thickness step (m) and the snow step (kg/m²)
 *   56  uint16   n, then n bytes of ASCII provenance; zero padding to a multiple of 8
 *   data: uint8 concentration[nLat][nLon] in percent, NO_SEA where the source has
 *         no sea; uint8 thickness[nLat][nLon], then snow[nLat][nLon], in their steps.
 */
export const SEA_ICE_START_FILE = 'data/seaice_march.bin';
export const EDGE_CONCENTRATION = 0.15;
export const NO_SEA = 255;
const MAGIC = 'SIC1', deg = 180 / Math.PI;

export function encodeSeaIceStart({ coverGrid, iceGrid, cover, thickness, snow, thicknessStep = 0.05, snowStep = 2, source = '' }) {
  const text = Uint8Array.from(source, (ch) => ch.charCodeAt(0)), offset = Math.ceil((58 + text.length) / 8) * 8;
  const bytes = new Uint8Array(offset + cover.length + thickness.length + snow.length), view = new DataView(bytes.buffer);
  [...MAGIC].forEach((ch, n) => view.setUint8(n, ch.charCodeAt(0)));
  view.setUint32(4, offset, true);
  [coverGrid.nLon, coverGrid.nLat, iceGrid.nLon, iceGrid.nLat].forEach((v, n) => view.setUint16(8 + 2 * n, v, true));
  [coverGrid.lon0, coverGrid.dLon, coverGrid.lat0, coverGrid.dLat, iceGrid.lon0, iceGrid.dLon, iceGrid.lat0, iceGrid.dLat, thicknessStep, snowStep].forEach((v, n) => view.setFloat32(16 + 4 * n, v, true));
  view.setUint16(56, text.length, true);
  bytes.set(text, 58);
  bytes.set(cover, offset);
  bytes.set(thickness, offset + cover.length);
  bytes.set(snow, offset + cover.length + thickness.length);
  return bytes;
}

const unit = (latDeg, lonDeg) => [Math.cos(latDeg / deg) * Math.cos(lonDeg / deg), Math.cos(latDeg / deg) * Math.sin(lonDeg / deg), Math.sin(latDeg / deg)];

function bilinear({ nLon, nLat, lon0, dLon, lat0, dLat }, values, held, latDeg, lonDeg) {
  const y = (latDeg - lat0) / dLat, x = (lonDeg - lon0) / dLon, r0 = Math.floor(y), c0 = Math.floor(x), fy = y - r0, fx = x - c0;
  let sum = 0, weight = 0;
  for (const [r, c, w] of [[r0, c0, (1 - fy) * (1 - fx)], [r0, c0 + 1, (1 - fy) * fx], [r0 + 1, c0, fy * (1 - fx)], [r0 + 1, c0 + 1, fy * fx]]) {
    if (r < 0 || r >= nLat || !(w > 0)) continue;
    const n = r * nLon + (((c % nLon) + nLon) % nLon);
    if (held(n)) { sum += w * values[n]; weight += w; }
  }
  return weight > 0 ? sum / weight : null;
}

function nearest({ nLon, nLat, lon0, dLon, lat0, dLat }, held, latDeg, lonDeg, reach) {
  const rows = Math.ceil(reach / Math.abs(dLat)), r0 = Math.round((latDeg - lat0) / dLat), c0 = Math.round((lonDeg - lon0) / dLon);
  const cols = Math.min(Math.ceil(nLon / 2), Math.ceil(reach / Math.abs(dLon) / Math.max(0.05, Math.cos(latDeg / deg))));
  const x = unit(latDeg, lonDeg);
  let best = -1, bestDot = Math.cos(reach / deg);
  for (let r = Math.max(0, r0 - rows); r <= Math.min(nLat - 1, r0 + rows); r++) {
    for (let dc = -cols; dc <= cols; dc++) {
      const c = (((c0 + dc) % nLon) + nLon) % nLon, n = r * nLon + c;
      if (!held(n)) continue;
      const p = unit(lat0 + r * dLat, lon0 + c * dLon), dot = x[0] * p[0] + x[1] * p[1] + x[2] * p[2];
      if (dot >= bestDot) { bestDot = dot; best = n; }
    }
  }
  return best;
}

export function decodeSeaIceStart(bytes, { reach = 1 } = {}) {
  const data = ArrayBuffer.isView(bytes) ? new Uint8Array(bytes.buffer, bytes.byteOffset, bytes.byteLength) : new Uint8Array(bytes);
  const view = new DataView(data.buffer, data.byteOffset, data.byteLength);
  const magic = String.fromCharCode(...data.subarray(0, 4));
  if (magic !== MAGIC) throw new Error(`not a sea-ice start (${magic})`);
  const offset = view.getUint32(4, true), [nLon, nLat, iceLon, iceLat] = [0, 1, 2, 3].map((n) => view.getUint16(8 + 2 * n, true));
  const f = Array.from({ length: 10 }, (_, n) => view.getFloat32(16 + 4 * n, true));
  const coverGrid = { nLon, nLat, lon0: f[0], dLon: f[1], lat0: f[2], dLat: f[3] }, iceGrid = { nLon: iceLon, nLat: iceLat, lon0: f[4], dLon: f[5], lat0: f[6], dLat: f[7] };
  const [thicknessStep, snowStep] = [f[8], f[9]], plane = iceLon * iceLat;
  const source = String.fromCharCode(...data.subarray(58, 58 + view.getUint16(56, true)));
  const cover = data.slice(offset, offset + nLon * nLat);
  const thickness = data.slice(offset + cover.length, offset + cover.length + plane), snow = data.slice(offset + cover.length + plane, offset + cover.length + 2 * plane);
  const sea = (n) => cover[n] !== NO_SEA, everywhere = () => true;
  function at(lat, lon) {
    const latDeg = lat * deg, lonDeg = lon * deg;
    let percent = bilinear(coverGrid, cover, sea, latDeg, lonDeg);
    if (percent === null) { const n = nearest(coverGrid, sea, latDeg, lonDeg, reach); percent = n < 0 ? 0 : cover[n]; }
    const concentration = Math.min(1, percent / 100), h = thicknessStep * (bilinear(iceGrid, thickness, everywhere, latDeg, lonDeg) ?? 0);
    if (concentration < EDGE_CONCENTRATION || !(h > 0)) return { concentration: 0, thickness: 0, snow: 0 };
    return { concentration, thickness: h, snow: snowStep * (bilinear(iceGrid, snow, everywhere, latDeg, lonDeg) ?? 0) };
  }
  return { source, coverGrid, iceGrid, thicknessStep, snowStep, cover, thickness, snow, at };
}

export async function loadSeaIceStart(source, options = {}) {
  if (!source) return null;
  if (typeof source === 'object' && typeof source.at === 'function') return source;
  if (typeof source === 'string' && !/^[a-z]+:\/\//i.test(source) && typeof process !== 'undefined' && process.versions && process.versions.node) {
    const { readFileSync } = await import('node:fs');
    return decodeSeaIceStart(readFileSync(source), options);
  }
  const response = await fetch(source);
  if (!response.ok) throw new Error(`sea-ice start ${source}: ${response.status}`);
  return decodeSeaIceStart(await response.arrayBuffer(), options);
}
