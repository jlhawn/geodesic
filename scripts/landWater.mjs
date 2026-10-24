/*
 * The land of GSHHG (Wessel and Smith 1996, J. Geophys. Res. 101, 8741;
 * version 2.3.7, LGPL) on a latitude–longitude raster (rows from the
 * north, columns from −180°). readGshhg reads a gshhs_?.b file: per
 * polygon a header of eleven big-endian int32 (id, points, flag, west,
 * east, south, north, area, full area, container, ancestor; the flag's
 * low byte is the level) and its points as int32 pairs of microdegrees.
 * Levels 1 to 4 nest (land, lake, island in a lake, pond on that island),
 * so a point is land where it lies inside an odd number of them; levels 5
 * and 6, the Antarctic ice front and grounding line, are left out.
 */
import { readFileSync } from 'node:fs';

export function readGshhg(path) {
  const bytes = readFileSync(path), view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength), polygons = [];
  for (let at = 0; at + 44 <= bytes.byteLength;) {
    const n = view.getInt32(at + 4), flag = view.getInt32(at + 8);
    const xy = new Float64Array(2 * n);
    for (let p = 0; p < n; p++) { xy[2 * p] = view.getInt32(at + 44 + 8 * p) * 1e-6; xy[2 * p + 1] = view.getInt32(at + 48 + 8 * p) * 1e-6; }
    polygons.push({ id: view.getInt32(at), level: flag & 255, xy });
    at += 44 + 8 * n;
  }
  return polygons;
}

/*
 * 1 for each raster point whose whole cell lies on land: the sub × sub
 * sample points spread across it are all inside an odd number of the
 * polygons (each {level, xy}, xy alternating longitude and latitude in
 * degrees, closed or not), so a channel through any part of the cell keeps
 * it from counting.
 */
export function landRaster(polygons, rows, cols, { levels = [1, 2, 3, 4], sub = 3, band = 120 } = {}) {
  const R = rows * sub, W = cols * sub, dlat = 180 / R, dlon = 360 / W, bands = Math.ceil(R / band);
  const counts = new Int32Array(bands + 1);
  const unwrap = (polygon) => {
    const xy = polygon.xy, n = xy.length / 2, x = new Float64Array(n + 1), y = new Float64Array(n + 1);
    for (let p = 0; p < n; p++) {
      let v = xy[2 * p];
      if (p > 0) { while (v - x[p - 1] > 180) v -= 360; while (v - x[p - 1] < -180) v += 360; }
      x[p] = v; y[p] = xy[2 * p + 1];
    }
    x[n] = x[0]; y[n] = y[0];
    if (Math.abs(x[n] - x[n - 1]) > 180) x[n] = x[n - 1] + ((x[0] - x[n - 1] + 540) % 360) - 180;
    return { x, y };
  };
  const rowOf = (lat) => (90 - lat) / dlat - 0.5;
  const kept = polygons.filter((p) => levels.includes(p.level)).map(unwrap);
  const edgesIn = Array.from({ length: bands }, () => []);
  kept.forEach(({ x, y }, k) => {
    for (let p = 0; p + 1 < x.length; p++) {
      const a = rowOf(y[p]), b = rowOf(y[p + 1]), r0 = Math.ceil(Math.min(a, b)), r1 = Math.floor(Math.max(a, b) - 1e-12);
      if (r1 < r0 || r1 < 0 || r0 >= R) continue;
      for (let s = Math.max(0, Math.floor(r0 / band)); s <= Math.min(bands - 1, Math.floor(r1 / band)); s++) edgesIn[s].push(k, p);
    }
  });
  const land = new Uint8Array(rows * cols);
  const marks = new Uint8Array(band * (W + 1)), hits = [];
  for (let s = 0; s < bands; s++) {
    marks.fill(0); hits.length = 0;
    const top = s * band, bottom = Math.min(R, top + band), list = edgesIn[s];
    for (let q = 0; q < list.length; q += 2) {
      const { x, y } = kept[list[q]], p = list[q + 1];
      const a = rowOf(y[p]), b = rowOf(y[p + 1]);
      const r0 = Math.max(top, Math.ceil(Math.min(a, b))), r1 = Math.min(bottom - 1, Math.floor(Math.max(a, b) - 1e-12));
      for (let r = r0; r <= r1; r++) {
        const t = (r - a) / (b - a);
        hits.push(r - top, list[q], x[p] + t * (x[p + 1] - x[p]));
      }
    }
    const order = Array.from({ length: hits.length / 3 }, (_, h) => h);
    order.sort((u, v) => hits[3 * u] - hits[3 * v] || hits[3 * u + 1] - hits[3 * v + 1] || hits[3 * u + 2] - hits[3 * v + 2]);
    const toggle = (row, from, to) => {
      const c0 = Math.ceil((from + 180) / dlon - 0.5), c1 = Math.ceil((to + 180) / dlon - 0.5);
      if (c1 <= c0) return;
      const span = (u, v) => { marks[row * (W + 1) + u] ^= 1; marks[row * (W + 1) + v] ^= 1; };
      let u = ((c0 % W) + W) % W, length = Math.min(W, c1 - c0);
      while (length > 0) { const v = Math.min(W, u + length); span(u, v); length -= v - u; u = 0; }
    };
    for (let h = 0; h + 1 < order.length;) {
      const first = order[h], next = order[h + 1];
      if (hits[3 * first] === hits[3 * next] && hits[3 * first + 1] === hits[3 * next + 1]) { toggle(hits[3 * first], hits[3 * first + 2], hits[3 * next + 2]); h += 2; } else h += 1;
    }
    const inside = new Uint8Array(W);
    for (let r = top; r < bottom; r++) {
      let on = 0;
      const base = (r - top) * (W + 1);
      for (let c = 0; c < W; c++) { on ^= marks[base + c]; inside[c] = on; }
      const row = Math.floor(r / sub), first = r % sub === 0;
      for (let c = 0; c < cols; c++) {
        let all = 1;
        for (let k = 0; k < sub; k++) all &= inside[c * sub + k];
        const n = row * cols + c;
        land[n] = first ? all : land[n] & all;
      }
    }
  }
  return land;
}
