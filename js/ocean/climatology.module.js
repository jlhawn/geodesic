/*
 * An ocean climatology on a regular latitude–longitude grid at a few
 * standard depths: potential temperature and practical salinity, stored
 * as int16 with a scale and offset (data/woa_annual_1deg.bin, the World
 * Ocean Atlas 2023 annual mean written by scripts/packWoa.py, which
 * documents the layout). `columnAt(lat, lon)`, in radians, gives the
 * column at a point as { depths, T, S }, T in kelvin, over the levels
 * that hold water there: each level interpolated bilinearly between the
 * four surrounding grid points when all four hold water at that level,
 * otherwise taken from the nearest of them that does, the column ending
 * at the first level none of them reaches. A point whose four
 * surrounding grid points are all land takes the nearest column within
 * `reach` degrees whole; null when there is none.
 */
export const CLIMATOLOGY_FILE = 'data/woa_annual_1deg.bin';
const MAGIC = 'WOA1', KELVIN = 273.15;

export function decodeClimatology(bytes, { reach = 5 } = {}) {
  const view = ArrayBuffer.isView(bytes) ? new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength) : new DataView(bytes);
  const magic = String.fromCharCode(...[0, 1, 2, 3].map((n) => view.getUint8(n)));
  if (magic !== MAGIC) throw new Error(`not an ocean climatology (${magic})`);
  const offset = view.getUint32(4, true), nLon = view.getUint16(8, true), nLat = view.getUint16(10, true), nDepth = view.getUint16(12, true), missing = view.getInt16(14, true);
  const [lon0, dLon, lat0, dLat, tOffset, tScale, sOffset, sScale] = [0, 1, 2, 3, 4, 5, 6, 7].map((n) => view.getFloat32(16 + 4 * n, true));
  const depths = Float64Array.from({ length: nDepth }, (_, j) => view.getFloat32(48 + 4 * j, true));
  const textAt = 48 + 4 * nDepth, textLength = view.getUint16(textAt, true);
  const source = String.fromCharCode(...Array.from({ length: textLength }, (_, n) => view.getUint8(textAt + 2 + n)));
  const count = nDepth * nLat * nLon;
  let T, S;
  if (new Uint8Array(Uint16Array.of(1).buffer)[0] === 1) {
    const data = new Uint8Array(view.buffer, view.byteOffset + offset, 4 * count).slice().buffer;
    T = new Int16Array(data, 0, count); S = new Int16Array(data, 2 * count, count);
  } else {
    T = new Int16Array(count); S = new Int16Array(count);
    for (let n = 0; n < count; n++) { T[n] = view.getInt16(offset + 2 * n, true); S[n] = view.getInt16(offset + 2 * (count + n), true); }
  }
  return climatology({ nLon, nLat, nDepth, lon0, dLon, lat0, dLat, tOffset, tScale, sOffset, sScale, missing, depths, T, S, source }, reach);
}

/*
 * The file for a grid of values: T in °C and S in psu, each
 * [depth][lat][lon] flattened, NaN where missing.
 */
export function encodeClimatology({ nLon, nLat, lon0, dLon, lat0, dLat, depths, T, S, source = '', tOffset = 15, tScale = 0.001, sOffset = 20, sScale = 0.001, missing = -32768 }) {
  const nDepth = depths.length, count = nDepth * nLat * nLon;
  let head = 48 + 4 * nDepth + 2 + source.length;
  head += (8 - (head % 8)) % 8;
  const buffer = new ArrayBuffer(head + 4 * count), view = new DataView(buffer);
  for (let n = 0; n < 4; n++) view.setUint8(n, MAGIC.charCodeAt(n));
  view.setUint32(4, head, true); view.setUint16(8, nLon, true); view.setUint16(10, nLat, true); view.setUint16(12, nDepth, true); view.setInt16(14, missing, true);
  [lon0, dLon, lat0, dLat, tOffset, tScale, sOffset, sScale].forEach((value, n) => view.setFloat32(16 + 4 * n, value, true));
  depths.forEach((d, j) => view.setFloat32(48 + 4 * j, d, true));
  view.setUint16(48 + 4 * nDepth, source.length, true);
  for (let n = 0; n < source.length; n++) view.setUint8(50 + 4 * nDepth + n, source.charCodeAt(n));
  const pack = (value, off, scale) => (Number.isFinite(value) ? Math.round((value - off) / scale) : missing);
  for (let n = 0; n < count; n++) { view.setInt16(head + 2 * n, pack(T[n], tOffset, tScale), true); view.setInt16(head + 2 * (count + n), pack(S[n], sOffset, sScale), true); }
  return new Uint8Array(buffer);
}

/*
 * The climatology from a file path (Node) or URL (fetch), an already
 * decoded one as it is, or null for none.
 */
export async function loadClimatology(source, options = {}) {
  if (!source) return null;
  if (typeof source.columnAt === 'function') return source;
  if (typeof source === 'string' && !/^[a-z]+:\/\//i.test(source) && typeof process !== 'undefined' && process.versions && process.versions.node) {
    const { readFileSync } = await import('node:fs');
    return decodeClimatology(readFileSync(source), options);
  }
  const response = await fetch(source);
  if (!response.ok) throw new Error(`ocean climatology ${source}: ${response.status}`);
  return decodeClimatology(await response.arrayBuffer(), options);
}

/*
 * Temperature and salinity at depth z in a column, linear between its
 * levels and constant beyond its first and last.
 */
export function profileAt({ depths, T, S }, z) {
  const n = depths.length;
  if (z <= depths[0]) return [T[0], S[0]];
  if (z >= depths[n - 1]) return [T[n - 1], S[n - 1]];
  let j = 1;
  while (depths[j] < z) j++;
  const f = (z - depths[j - 1]) / (depths[j] - depths[j - 1]);
  return [T[j - 1] + f * (T[j] - T[j - 1]), S[j - 1] + f * (S[j] - S[j - 1])];
}

function climatology(grid, reach) {
  const { nLon, nLat, nDepth, lon0, dLon, lat0, dLat, tOffset, tScale, sOffset, sScale, missing, depths, T, S } = grid;
  const plane = nLat * nLon, deg = 180 / Math.PI;
  const at = (j, r, c) => j * plane + r * nLon + c;
  const wrap = (c) => ((c % nLon) + nLon) % nLon;
  const valid = (j, r, c) => T[at(j, r, c)] !== missing && S[at(j, r, c)] !== missing;
  const levels = (r, c) => { let j = 0; while (j < nDepth && valid(j, r, c)) j++; return j; };
  const column = (n) => ({ depths: depths.slice(0, n), T: new Float64Array(n), S: new Float64Array(n) });
  function nearestColumn(latDeg, lonDeg) {
    const rows = Math.ceil(reach / Math.abs(dLat)), r0 = Math.round((latDeg - lat0) / dLat), c0 = Math.round((lonDeg - lon0) / dLon);
    const cosLat = Math.max(0.05, Math.cos(latDeg / deg)), cols = Math.min(Math.ceil(nLon / 2), Math.ceil(reach / Math.abs(dLon) / cosLat));
    const x = [Math.cos(latDeg / deg) * Math.cos(lonDeg / deg), Math.cos(latDeg / deg) * Math.sin(lonDeg / deg), Math.sin(latDeg / deg)];
    let best = null, bestDistance = Infinity;
    for (let r = Math.max(0, r0 - rows); r <= Math.min(nLat - 1, r0 + rows); r++) {
      for (let dc = -cols; dc <= cols; dc++) {
        const c = wrap(c0 + dc);
        if (!valid(0, r, c)) continue;
        const la = (lat0 + r * dLat) / deg, lo = (lon0 + c * dLon) / deg;
        const distance = Math.acos(Math.max(-1, Math.min(1, x[0] * Math.cos(la) * Math.cos(lo) + x[1] * Math.cos(la) * Math.sin(lo) + x[2] * Math.sin(la)))) * deg;
        if (distance <= reach && distance < bestDistance) { bestDistance = distance; best = [r, c]; }
      }
    }
    if (!best) return null;
    const [r, c] = best, out = column(levels(r, c));
    for (let j = 0; j < out.depths.length; j++) { const n = at(j, r, c); out.T[j] = tOffset + tScale * T[n] + KELVIN; out.S[j] = sOffset + sScale * S[n]; }
    return out;
  }
  function columnAt(lat, lon) {
    const latDeg = lat * deg, lonDeg = lon * deg;
    const y = Math.max(0, Math.min(nLat - 1, (latDeg - lat0) / dLat)), r0 = Math.min(nLat - 2, Math.floor(y)), fy = y - r0;
    const x = (lonDeg - lon0) / dLon, cx = Math.floor(x), fx = x - cx, c0 = wrap(cx), c1 = wrap(cx + 1);
    const corners = [[r0, c0, (1 - fy) * (1 - fx)], [r0, c1, (1 - fy) * fx], [r0 + 1, c0, fy * (1 - fx)], [r0 + 1, c1, fy * fx]];
    const depthOf = corners.map(([r, c]) => levels(r, c));
    const n = Math.max(...depthOf);
    if (n === 0) return nearestColumn(latDeg, lonDeg);
    const nearest = corners.map((corner, m) => [corner, depthOf[m]]).sort((a, b) => b[0][2] - a[0][2]);
    const out = column(n);
    for (let j = 0; j < n; j++) {
      let t = 0, s = 0;
      if (depthOf.every((d) => d > j)) {
        for (const [r, c, w] of corners) { const m = at(j, r, c); t += w * T[m]; s += w * S[m]; }
      } else {
        const [[r, c]] = nearest.find(([, d]) => d > j);
        const m = at(j, r, c); t = T[m]; s = S[m];
      }
      out.T[j] = tOffset + tScale * t + KELVIN;
      out.S[j] = sOffset + sScale * s;
    }
    return out;
  }
  return { ...grid, reach, columnAt };
}
