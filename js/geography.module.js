/*
 * Land and elevation on the mesh from a regular latitude–longitude
 * elevation raster: rows from north to south, columns from west to
 * east starting at −180°, each value the elevation in metres, negative
 * under the sea. Every raster point is assigned to the nearest cell by
 * walking the cell neighbour graph from the previous point's cell; a
 * cell's elevation is the mean of its points and its land fraction the
 * share of them above sea level. A cell is land when that fraction
 * exceeds landThreshold, and along the polylines of landBridges (the
 * Panama isthmus by default) when the mask would otherwise let a sea
 * path cross them; and the cells along the polylines of seaStraits
 * (Hormuz, Bab-el-Mandeb and Gibraltar by default) are sea at least
 * their sill deep when the mask would otherwise close them. A mesh
 * that resolves a bridge or strait on its own is left untouched. The
 * NARROW_STRAITS (Bosporus, the Danish straits) are defined but not
 * applied: a cell of 75–150 km opening a strait a few km wide would
 * exchange water far too freely. edgeOcean marks the edges between two ocean
 * cells, the only ones the ocean flows through; coastEdges lists the
 * edges between land and ocean.
 */
export function topographyFromInt16(buffer, rows = 720, cols = 1440) {
  return { rows, cols, data: new Int16Array(buffer) };
}

export function syntheticTopography(rows, cols, elevationAt) {
  const data = new Int16Array(rows * cols);
  for (let r = 0; r < rows; r++) {
    const lat = Math.PI / 2 - (r + 0.5) * Math.PI / rows;
    for (let c = 0; c < cols; c++) data[r * cols + c] = elevationAt(lat, -Math.PI + (c + 0.5) * 2 * Math.PI / cols);
  }
  return { rows, cols, data };
}

export const LAND_BRIDGES = {
  panama: [[9.8, -84.5], [9.0, -82.5], [8.8, -80.5], [9.2, -79.2], [8.5, -77.8], [7.5, -77.0]],
};

export const SEA_STRAITS = {
  hormuz: { sill: 90, points: [[26.6, 53.5], [26.4, 55.5], [26.2, 56.6], [25.5, 57.2]] },
  babElMandeb: { sill: 150, points: [[14.0, 42.0], [13.0, 43.0], [12.6, 43.6], [12.3, 44.5]] },
  gibraltar: { sill: 300, points: [[36.1, -6.5], [35.95, -5.6], [36.0, -4.8]] },
};

export const NARROW_STRAITS = {
  bosporus: { sill: 50, points: [[41.6, 29.0], [41.1, 29.1], [40.6, 28.0], [40.2, 26.5], [39.9, 25.8]] },
  danish: { sill: 50, points: [[55.5, 12.7], [56.1, 12.3], [56.8, 11.5], [57.6, 10.5]] },
};

export function createGeography(mesh, topography, { landThreshold = 0.5, landBridges = LAND_BRIDGES, seaStraits = SEA_STRAITS } = {}) {
  const { nCells: C, nEdges: E, xCell, cellsOnCell, nEdgesOnCell, cellsOnEdge, areaCell, latCell, lonCell } = mesh;
  const { rows, cols, data } = topography;
  const count = new Int32Array(C), landCount = new Int32Array(C), sum = new Float64Array(C);
  const nearest = (x, y, z, start) => {
    let best = start, bestDot = x * xCell[3 * start] + y * xCell[3 * start + 1] + z * xCell[3 * start + 2];
    for (;;) {
      let next = best;
      for (let m = 0; m < nEdgesOnCell[best]; m++) {
        const j = cellsOnCell[mesh.maxEdges * best + m];
        const dot = x * xCell[3 * j] + y * xCell[3 * j + 1] + z * xCell[3 * j + 2];
        if (dot > bestDot) { bestDot = dot; next = j; }
      }
      if (next === best) return best;
      best = next;
    }
  };
  let cell = 0;
  for (let r = 0; r < rows; r++) {
    const lat = Math.PI / 2 - (r + 0.5) * Math.PI / rows, cosLat = Math.cos(lat), sinLat = Math.sin(lat);
    for (let c = 0; c < cols; c++) {
      const lon = -Math.PI + (c + 0.5) * 2 * Math.PI / cols;
      cell = nearest(cosLat * Math.cos(lon), cosLat * Math.sin(lon), sinLat, cell);
      const elevation = data[r * cols + c];
      count[cell]++;
      sum[cell] += elevation;
      if (elevation > 0) landCount[cell]++;
    }
  }
  const land = new Uint8Array(C), landFraction = new Float64Array(C), elevation = new Float64Array(C);
  for (let i = 0; i < C; i++) {
    if (count[i] === 0) {
      const r = Math.min(rows - 1, Math.max(0, Math.floor((Math.PI / 2 - latCell[i]) / (Math.PI / rows))));
      const c = ((Math.floor((lonCell[i] + Math.PI) / (2 * Math.PI / cols)) % cols) + cols) % cols;
      const value = data[r * cols + c];
      count[i] = 1; sum[i] = value; landCount[i] = value > 0 ? 1 : 0;
    }
    landFraction[i] = landCount[i] / count[i];
    elevation[i] = sum[i] / count[i];
    land[i] = landFraction[i] > landThreshold ? 1 : 0;
  }
  const trace = (points, visit) => {
    let start = 0;
    for (let p = 0; p + 1 < points.length; p++) {
      const [la0, lo0] = points[p], [la1, lo1] = points[p + 1];
      const steps = Math.max(1, Math.ceil(Math.hypot(la1 - la0, (lo1 - lo0) * Math.cos(la0 * Math.PI / 180)) / 0.2));
      for (let n = 0; n <= steps; n++) {
        const la = (la0 + (la1 - la0) * n / steps) * Math.PI / 180, lo = (lo0 + (lo1 - lo0) * n / steps) * Math.PI / 180;
        start = nearest(Math.cos(la) * Math.cos(lo), Math.cos(la) * Math.sin(lo), Math.sin(la), start);
        visit(start);
      }
    }
  };
  const deg = 180 / Math.PI;
  const seaPath = (points) => {
    const lats = points.map((q) => q[0]), lons = points.map((q) => q[1]);
    const box = [Math.min(...lats) - 4, Math.max(...lats) + 4, Math.min(...lons) - 4, Math.max(...lons) + 4];
    const inBox = (i) => latCell[i] * deg >= box[0] && latCell[i] * deg <= box[1] && lonCell[i] * deg >= box[2] && lonCell[i] * deg <= box[3];
    const cells = []; trace(points, (i) => cells.push(i));
    const ends = [cells[0], cells[cells.length - 1]].map((c) => {
      if (!land[c]) return c;
      for (let m = 0; m < nEdgesOnCell[c]; m++) { const j = cellsOnCell[mesh.maxEdges * c + m]; if (!land[j]) return j; }
      return -1;
    });
    if (ends[0] < 0 || ends[1] < 0) return false;
    const seen = new Uint8Array(C), queue = [ends[0]]; seen[ends[0]] = 1;
    while (queue.length) {
      const i = queue.pop(); if (i === ends[1]) return true;
      for (let m = 0; m < nEdgesOnCell[i]; m++) { const j = cellsOnCell[mesh.maxEdges * i + m]; if (j < 0 || seen[j] || land[j] || !inBox(j)) continue; seen[j] = 1; queue.push(j); }
    }
    return false;
  };
  for (const points of Object.values(landBridges)) {
    if (!seaPath(points)) continue;
    trace(points, (i) => { land[i] = 1; elevation[i] = Math.max(elevation[i], 1); });
  }
  for (const { sill, points } of Object.values(seaStraits)) {
    if (seaPath(points)) continue;
    trace(points, (i) => { land[i] = 0; elevation[i] = Math.min(elevation[i], -sill); });
  }
  const edgeOcean = new Uint8Array(E);
  const coast = [];
  for (let e = 0; e < E; e++) {
    const a = land[cellsOnEdge[2 * e]], b = land[cellsOnEdge[2 * e + 1]];
    edgeOcean[e] = !a && !b ? 1 : 0;
    if (a !== b) coast.push(e);
  }
  // The topography carries no ice mask: the ice sheets are Antarctica's land and Greenland's interior above 800 m.
  const iceSheet = new Uint8Array(C);
  for (let i = 0; i < C; i++) {
    if (!land[i]) continue;
    const lat = latCell[i] * 180 / Math.PI, lon = lonCell[i] * 180 / Math.PI;
    if (lat < -60 || (lat > 60 && lat < 84 && lon > -73 && lon < -12 && elevation[i] > 800)) iceSheet[i] = 1;
  }
  let landArea = 0, total = 0;
  for (let i = 0; i < C; i++) { total += areaCell[i]; if (land[i]) landArea += areaCell[i]; }
  return { land, landFraction, elevation, iceSheet, edgeOcean, coastEdges: Int32Array.from(coast), landArea: landArea / total };
}

/*
 * The surface geopotential for the dynamical core: the land elevation,
 * clamped at sea level and smoothed on the mesh so that nothing sits at
 * the grid scale, times g. Each pass relaxes every land cell toward the
 * mean of its neighbours; ocean cells stay at sea level and act as the
 * boundary, so coasts slope down to the sea over `passes` cells.
 */
export function surfaceGeopotential(mesh, geography, { g = 9.80616, passes = 2, weight = 0.5 } = {}) {
  const { nCells: C, maxEdges, nEdgesOnCell, cellsOnCell } = mesh;
  const { land, elevation } = geography;
  let height = Float64Array.from(elevation, (h, i) => (land[i] ? Math.max(0, h) : 0));
  for (let pass = 0; pass < passes; pass++) {
    const next = new Float64Array(C);
    for (let i = 0; i < C; i++) {
      if (!land[i]) continue;
      let sum = 0;
      for (let m = 0; m < nEdgesOnCell[i]; m++) sum += height[cellsOnCell[maxEdges * i + m]];
      next[i] = (1 - weight) * height[i] + weight * sum / nEdgesOnCell[i];
    }
    height = next;
  }
  return Float64Array.from(height, (h) => g * h);
}

/*
 * The subgrid orography of each cell after Baines and Palmer (1990), as
 * the IFS computes its fields for Lott and Miller's (1997) scheme (IFS
 * Cy47r3 Part IV §11.3.4): the raster, clamped at sea level, less the
 * orography the dynamics resolves (`resolved`, per cell in metres,
 * interpolated linearly to each raster point on the triangle of cell
 * centres that holds it); every raster point's gradient by central
 * differences, K = ½(h_x² + h_y²), L = ½(h_x² − h_y²), M = h_x h_y, h
 * and h² averaged over the points nearest each cell with area weights;
 * then the standard deviation μ² = ⟨h²⟩ − ⟨h⟩², the anisotropy
 * γ² = (K − √(L² + M²))/(K + √(L² + M²)), the orientation θ = ½ atan2(M, L)
 * of the principal axis (radians from east towards north, the direction
 * of the steepest mean-square slope) and the slope σ² = K + √(L² + M²).
 * x is east and y north; `count` is the number of raster points a cell
 * holds. With `land`, a cell the mask makes sea has no subgrid orography.
 * With `filtered`, a raster on the same grid of the band-passed variance
 * (σ_flt² per point, scripts/subgridTerrain.py), `filtered` per cell is
 * the square root of its area-weighted mean. With `spacing` (metres) the
 * central differences reach that far east, west, north and south
 * (bilinear on the raster, at least one point), as on the IFS's 5 km grid.
 */
export function subgridOrography(mesh, topography, resolved = null, land = null, { filtered: band = null, spacing = 0 } = {}) {
  const { nCells: C, xCell, cellsOnCell, nEdgesOnCell, maxEdges, verticesOnCell, cellsOnVertex, radius } = mesh;
  const { rows, cols, data } = topography;
  const R = radius ?? 6371220;
  const residual = new Float64Array(rows * cols), owner = new Int32Array(rows * cols);
  const nearest = (x, y, z, start) => {
    let best = start, bestDot = x * xCell[3 * start] + y * xCell[3 * start + 1] + z * xCell[3 * start + 2];
    for (;;) {
      let next = best;
      for (let m = 0; m < nEdgesOnCell[best]; m++) {
        const j = cellsOnCell[maxEdges * best + m];
        const dot = x * xCell[3 * j] + y * xCell[3 * j + 1] + z * xCell[3 * j + 2];
        if (dot > bestDot) { bestDot = dot; next = j; }
      }
      if (next === best) return best;
      best = next;
    }
  };
  const height = (i) => (resolved ? resolved[i] : 0);
  const at = (x, y, z, i) => {
    if (!resolved) return 0;
    let bestWeight = -Infinity, value = height(i);
    for (let m = 0; m < nEdgesOnCell[i]; m++) {
      const v = verticesOnCell[maxEdges * i + m];
      if (v < 0) continue;
      const a = cellsOnVertex[3 * v], b = cellsOnVertex[3 * v + 1], c = cellsOnVertex[3 * v + 2];
      const ax = xCell[3 * a], ay = xCell[3 * a + 1], az = xCell[3 * a + 2], bx = xCell[3 * b], by = xCell[3 * b + 1], bz = xCell[3 * b + 2], cx = xCell[3 * c], cy = xCell[3 * c + 1], cz = xCell[3 * c + 2];
      const det = ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx);
      const wa = (x * (by * cz - bz * cy) - y * (bx * cz - bz * cx) + z * (bx * cy - by * cx)) / det;
      const wb = (ax * (y * cz - z * cy) - ay * (x * cz - z * cx) + az * (x * cy - y * cx)) / det;
      const wc = (ax * (by * z - bz * y) - ay * (bx * z - bz * x) + az * (bx * y - by * x)) / det;
      const sum = wa + wb + wc, low = Math.min(wa, wb, wc) / sum;
      if (low > bestWeight) { bestWeight = low; value = (wa * height(a) + wb * height(b) + wc * height(c)) / sum; }
    }
    return value;
  };
  let cell = 0;
  for (let r = 0; r < rows; r++) {
    const lat = Math.PI / 2 - (r + 0.5) * Math.PI / rows, cosLat = Math.cos(lat), sinLat = Math.sin(lat);
    for (let c = 0; c < cols; c++) {
      const lon = -Math.PI + (c + 0.5) * 2 * Math.PI / cols, x = cosLat * Math.cos(lon), y = cosLat * Math.sin(lon);
      cell = nearest(x, y, sinLat, cell);
      owner[r * cols + c] = cell;
      residual[r * cols + c] = Math.max(0, data[r * cols + c]) - at(x, y, sinLat, cell);
    }
  }
  const weight = new Float64Array(C), sum = new Float64Array(C), square = new Float64Array(C), kSum = new Float64Array(C), lSum = new Float64Array(C), mSum = new Float64Array(C), bandSum = new Float64Array(C), count = new Int32Array(C);
  const dLat = Math.PI / rows, dLon = 2 * Math.PI / cols;
  const wrap = (c) => ((c % cols) + cols) % cols;
  const sample = (row, column) => {
    const r0 = Math.floor(row), c0 = Math.floor(column), tr = row - r0, tc = column - c0, r1 = Math.min(rows - 1, r0 + 1), a = wrap(c0), b = wrap(c0 + 1);
    return (1 - tr) * ((1 - tc) * residual[r0 * cols + a] + tc * residual[r0 * cols + b]) + tr * ((1 - tc) * residual[r1 * cols + a] + tc * residual[r1 * cols + b]);
  };
  for (let r = 0; r < rows; r++) {
    const lat = Math.PI / 2 - (r + 0.5) * dLat, w = Math.cos(lat);
    const reachX = spacing > 0 ? Math.min(cols / 4, Math.max(1, spacing / (R * w * dLon))) : 1, dx = 2 * reachX * R * w * dLon;
    const reachY = spacing > 0 ? Math.max(1, spacing / (R * dLat)) : 1;
    const north = Math.max(0, r - reachY), south = Math.min(rows - 1, r + reachY), dy = (south - north) * R * dLat;
    for (let c = 0; c < cols; c++) {
      const n = r * cols + c, i = owner[n], h = residual[n];
      const hx = spacing > 0 ? (sample(r, c + reachX) - sample(r, c - reachX)) / dx : (residual[r * cols + (c + 1) % cols] - residual[r * cols + (c + cols - 1) % cols]) / dx;
      const hy = spacing > 0 ? (sample(north, c) - sample(south, c)) / dy : (residual[north * cols + c] - residual[south * cols + c]) / dy;
      weight[i] += w; sum[i] += w * h; square[i] += w * h * h; count[i]++;
      if (band) bandSum[i] += w * band.data[n];
      kSum[i] += w * 0.5 * (hx * hx + hy * hy); lSum[i] += w * 0.5 * (hx * hx - hy * hy); mSum[i] += w * hx * hy;
    }
  }
  const deviation = new Float64Array(C), anisotropy = new Float64Array(C), orientation = new Float64Array(C), slope = new Float64Array(C), filtered = band ? new Float64Array(C) : null;
  for (let i = 0; i < C; i++) {
    if (!(weight[i] > 0) || (land && !land[i])) continue;
    if (band) filtered[i] = Math.sqrt(Math.max(0, bandSum[i] / weight[i]));
    const mean = sum[i] / weight[i], kk = kSum[i] / weight[i], ll = lSum[i] / weight[i], mm = mSum[i] / weight[i], spread = Math.hypot(ll, mm);
    deviation[i] = Math.sqrt(Math.max(0, square[i] / weight[i] - mean * mean));
    anisotropy[i] = kk + spread > 0 ? Math.sqrt(Math.max(0, (kk - spread) / (kk + spread))) : 1;
    orientation[i] = 0.5 * Math.atan2(mm, ll);
    slope[i] = Math.sqrt(kk + spread);
  }
  return { deviation, anisotropy, orientation, slope, count, ...(band ? { filtered } : {}) };
}

/*
 * The per-mesh subgrid orography files data/subgrid_N<N>.bin that
 * scripts/subgridTerrain.mjs writes from GMTED2010: a 16-byte header
 * ('SGO1', the cell count, the field count, 0) and per field one uint16
 * per cell, value = q × scale + offset, in SUBGRID_FIELDS' order.
 */
export const SUBGRID_FIELDS = [
  ['deviation', 0.1, 0],
  ['anisotropy', 1 / 65535, 0],
  ['orientation', Math.PI / 65535, -Math.PI / 2],
  ['slope', 1e-5, 0],
  ['filtered', 0.05, 0],
];

export function encodeSubgrid(fields) {
  const C = fields.deviation.length, out = new ArrayBuffer(16 + 2 * C * SUBGRID_FIELDS.length);
  const head = new DataView(out);
  [83, 71, 79, 49].forEach((b, j) => head.setUint8(j, b));
  head.setUint32(4, C, true); head.setUint32(8, SUBGRID_FIELDS.length, true);
  SUBGRID_FIELDS.forEach(([name, scale, offset], f) => {
    const q = new Uint16Array(out, 16 + 2 * C * f, C), values = fields[name];
    for (let i = 0; i < C; i++) q[i] = Math.min(65535, Math.max(0, Math.round((values[i] - offset) / scale)));
  });
  return out;
}

export function decodeSubgrid(buffer) {
  const bytes = buffer instanceof ArrayBuffer ? buffer : buffer.buffer.slice(buffer.byteOffset, buffer.byteOffset + buffer.byteLength);
  const head = new DataView(bytes);
  if (String.fromCharCode(...new Uint8Array(bytes, 0, 4)) !== 'SGO1') throw new Error('not a subgrid orography file');
  const C = head.getUint32(4, true), count = head.getUint32(8, true);
  if (count !== SUBGRID_FIELDS.length || bytes.byteLength !== 16 + 2 * C * count) throw new Error('subgrid orography file of the wrong size');
  const fields = {};
  SUBGRID_FIELDS.forEach(([name, scale, offset], f) => {
    const q = new Uint16Array(bytes, 16 + 2 * C * f, C);
    fields[name] = Float64Array.from(q, (v) => v * scale + offset);
  });
  for (let i = 0; i < C; i++) if (!(fields.deviation[i] > 0)) fields.orientation[i] = 0;
  return fields;
}

export const subgridUrl = (N) => new URL(`../data/subgrid_N${N}.bin`, import.meta.url);

/*
 * The fields for a mesh: `given` when it has the mesh's cell count, else
 * under node the file for N if there is one, else null (the caller
 * computes them from the 0.25° raster and says so).
 */
export function meshSubgrid(mesh, given = undefined) {
  if (given === false) return null;
  if (given) return given.deviation.length === mesh.nCells ? given : null;
  const N = Math.round(Math.sqrt((mesh.nCells - 2) / 10));
  const fs = globalThis.process?.getBuiltinModule?.('node:fs');
  if (!fs || 10 * N * N + 2 !== mesh.nCells) return null;
  const url = subgridUrl(N);
  if (!fs.existsSync(url)) return null;
  return decodeSubgrid(fs.readFileSync(url));
}

/*
 * `fields` with the sea cells of `land` zeroed, or null when more than
 * 0.5 % of the land cells hold nothing: fields made for another land mask.
 */
export function landSubgrid(fields, land) {
  const masked = {};
  for (const [name] of SUBGRID_FIELDS) if (fields[name]) masked[name] = Float64Array.from(fields[name], (v, i) => (land[i] ? v : 0));
  let cells = 0, bare = 0;
  for (let i = 0; i < land.length; i++) {
    if (!land[i]) continue;
    cells++;
    if (!(masked.deviation[i] > 0 || masked.slope[i] > 0 || (masked.filtered && masked.filtered[i] > 0))) bare++;
  }
  return bare <= 0.005 * cells ? masked : null;
}

const SUBGRID_FALLBACK = {
  missing: (N) => `no data/subgrid_N${N}.bin for this mesh`,
  terrain: () => 'the terrain is off, and the files hold only the scales below the resolved terrain',
  land: (N) => `data/subgrid_N${N}.bin was made for another land mask`,
};
const fallbackSaid = new Set();
export function subgridFallbackNote(mesh, reason = 'missing') {
  const N = Math.round(Math.sqrt((mesh.nCells - 2) / 10));
  if (fallbackSaid.has(`${N} ${reason}`)) return;
  fallbackSaid.add(`${N} ${reason}`);
  console.log(`subgrid orography: ${SUBGRID_FALLBACK[reason](N)}; computed from the topography's raster (at 0.25°, σ about a third of the 5 km fields') with Lott and Miller's constants, and no form drag`);
}

/*
 * Moves a column's surface pressure between surface geopotentials
 * hydrostatically, so a state saved over one terrain starts balanced
 * over another: π scales by exp(−Δφ / (R T)) with T the lowest layer's
 * temperature.
 */
export function rebalanceSurfacePressure(core, pi, theta, fromPhi, toPhi, R = 287.04) {
  const { K, C, exnerLayer } = core.diagnostics;
  const bottom = (K - 1) * C;
  for (let i = 0; i < C; i++) {
    const delta = (toPhi ? toPhi[i] : 0) - (fromPhi ? fromPhi[i] : 0);
    if (delta !== 0) pi[i] *= Math.exp(-delta / (R * theta[bottom + i] * exnerLayer[bottom + i]));
  }
  return pi;
}
