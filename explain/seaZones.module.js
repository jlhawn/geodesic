/* Grids here are indexed i * nj + j, with i running north along the circles of latitude and j down through the
   boundaries between classes. contours() orients every line so that larger values lie on its right in a frame
   drawn with i to the right and j down the page; for the overturning streamfunction, the northward flow above
   each boundary, that orientation is the direction the water moves. */
import { classOverturning, classSection } from './oceanRun.module.js';

export const SHALLOW = 500, FLOOR = 4000, THICK = 1;

export function overturningGrid(header, month = null) {
  const flow = classOverturning(header, month), base = classSection(header, month).base, L = header.layers, M = L + 1, caps = flow.caps;
  const lats = [-90, ...caps, 90], ni = lats.length, psi = new Float64Array(ni * M), depth = new Float64Array(ni * M);
  const pole = (i, b) => { for (let k = 1; k <= L; k++) depth[i * M + k] = base[b * L + k - 1]; };
  pole(0, 0);
  pole(ni - 1, header.section.bins - 1);
  caps.forEach((_, c) => { psi.set(flow.psi.subarray(c * M, (c + 1) * M), (c + 1) * M); depth.set(flow.depth.subarray(c * M, (c + 1) * M), (c + 1) * M); });
  return { lats, ni, nj: M, psi, depth, flow };
}

export function northwardByClass(header, month = null) {
  const flow = month ? header.overturning.months.find((m) => m.name === month) : header.overturning, base = classSection(header, month).base;
  const L = header.layers, caps = header.overturning.caps, R = header.constants.radius, values = [];
  caps.forEach((lat, c) => {
    const span = 2 * Math.PI * R * Math.cos(lat * Math.PI / 180);
    for (let k = 0; k < L; k++) {
      const thick = (b) => base[b * L + k] - (k ? base[b * L + k - 1] : 0), h = 0.5 * (thick(c) + thick(c + 1));
      values.push(h < THICK ? null : 100 * flow.mass[c * L + k] * 1e6 / (span * h));
    }
  });
  return { lats: caps, values };
}

export function overturningCells(header, month = null) {
  const { flow } = overturningGrid(header, month), M = header.layers + 1, caps = flow.caps;
  const mixed = (lat) => flow.psi[caps.indexOf(lat) * M + 1];
  const deep = (sign, from, to) => {
    let best = { psi: 0 };
    caps.forEach((lat, c) => { if (lat >= from && lat <= to) for (let k = 0; k <= header.layers; k++) if (flow.depth[c * M + k] >= 1000 && sign * flow.psi[c * M + k] > sign * best.psi) best = { psi: flow.psi[c * M + k], latitude: lat, depth: flow.depth[c * M + k] }; });
    return best;
  };
  return { mixed, south: mixed(-10), north: mixed(10), deepNorth: deep(1, 40, 90), deepSouth: deep(-1, -90, -40) };
}

export function contours(values, ni, nj, level) {
  const at = (i, j) => values[i * nj + j], above = (i, j) => at(i, j) > level;
  const cross = (i0, j0, i1, j1) => { const t = (level - at(i0, j0)) / (at(i1, j1) - at(i0, j0)); return [i0 + t * (i1 - i0), j0 + t * (j1 - j0)]; };
  const segments = [];
  for (let i = 0; i + 1 < ni; i++) for (let j = 0; j + 1 < nj; j++) {
    const corners = [[i, j], [i + 1, j], [i + 1, j + 1], [i, j + 1]], edges = [];
    for (let e = 0; e < 4; e++) {
      const [a, b] = [corners[e], corners[(e + 1) % 4]];
      if (above(...a) !== above(...b)) edges.push({ key: a[1] === b[1] ? `h${Math.min(a[0], b[0])},${a[1]}` : `v${a[0]},${Math.min(a[1], b[1])}`, point: cross(...a, ...b), e });
    }
    if (edges.length === 2) segments.push(orient(edges[0], edges[1], corners, above));
    else if (edges.length === 4) {
      const center = 0.25 * (at(i, j) + at(i + 1, j) + at(i + 1, j + 1) + at(i, j + 1)) > level, high = above(i, j);
      const pairs = center === high ? [[0, 1], [2, 3]] : [[0, 3], [1, 2]];
      for (const [p, q] of pairs) segments.push(orient(edges[p], edges[q], corners, above));
    }
  }
  const starts = new Map(segments.map((s) => [s.from, s])), ends = new Map(segments.map((s) => [s.to, s])), used = new Set(), lines = [];
  for (const segment of segments) {
    if (used.has(segment)) continue;
    let first = segment;
    while (ends.has(first.from) && ends.get(first.from) !== segment && !used.has(ends.get(first.from))) first = ends.get(first.from);
    const points = [first.a];
    let s = first, closed = false;
    while (s && !used.has(s)) {
      used.add(s);
      points.push(s.b);
      s = starts.get(s.to);
      if (s === first) closed = true;
    }
    lines.push({ points, closed });
  }
  return lines;
}

function orient(p, q, corners, above) {
  const mid = [(p.point[0] + q.point[0]) / 2, (p.point[1] + q.point[1]) / 2];
  const nearest = corners.reduce((best, c) => (Math.hypot(c[0] - mid[0], c[1] - mid[1]) < Math.hypot(best[0] - mid[0], best[1] - mid[1]) ? c : best));
  const d = [q.point[0] - p.point[0], q.point[1] - p.point[1]], right = (nearest[0] - p.point[0]) * -d[1] + (nearest[1] - p.point[1]) * d[0] > 0;
  return right === above(...nearest) ? { from: p.key, to: q.key, a: p.point, b: q.point } : { from: q.key, to: p.key, a: q.point, b: p.point };
}
