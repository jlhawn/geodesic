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
 * edges between land and ocean. edgeSill is each ocean edge's sill: the
 * greatest depth to which the two cells' deepest raster points stay
 * joined through the raster points of those two cells (eight-connected),
 * so the deepest channel the two cells hold between them, 0 where raster
 * land parts them; where the two cells' points never touch (an edge
 * shorter than the raster spacing) the shallowest point along the great
 * circle between the centers, sampled every 0.1°. Along a seaStraits
 * polyline the edges between its cells are at least the strait's sill
 * deep, whether or not the mesh resolves it. With `landIce`
 * (decodeLandIce, on the topography's grid), a land cell is ice sheet
 * where ice covers more than half of its land; iceFraction is the share
 * of its raster points under ice, iceThickness the mean thickness of that
 * ice (m), shelfFraction the share afloat, shelfThickness the floating
 * ice's thickness spread over the whole cell (m of ice) and wetDepth the
 * depth of the sea beneath its points, land counting as none, averaged
 * over the cell: the water a cell holds per unit area. Without it the
 * ice sheets are Antarctica's land and Greenland's interior above 800 m.
 * With `iceShelves` 'sea' the floating ice is taken away: its points
 * count at the depth of the bed beneath them and hold no ice. With
 * 'cavity' the atmosphere's geography is that of 'land', and `ocean` is
 * the geography of the water: that of 'sea', except that a cell 'sea'
 * turns to sea stays land where less than cavityDepth of water lies under
 * its shelf. The others are cavities, land above and sea below, under a
 * shelf of draft iceDraft = (ρi/ρw)·shelfThickness, their bed at wetDepth
 * so that the column holds the water of the cell's points.
 */
export function topographyFromInt16(buffer, rows = 720, cols = 1440, name = 'etopo1') {
  return { rows, cols, data: new Int16Array(buffer), name };
}

/*
 * The rasters under data/ a geography is built on, 0.25° block means of
 * the elevation in meters: ETOPO1's ice surface (scripts/topography.mjs)
 * and ETOPO 2022's (scripts/landIce.mjs), the bed plus the ice above it.
 * A state records its geography where it differs from LEGACY_GEOGRAPHY,
 * so one without the record is on ETOPO1 with the latitude and elevation
 * rule; fresh starts take FRESH_GEOGRAPHY. `iceShelves` is the floating
 * ice: 'land', an ice surface at its freeboard with the coastline at the
 * shelf front; 'sea', open water down to the bed; or 'cavity', the ice
 * surface above and the sea beneath it.
 */
export const TOPOGRAPHIES = { etopo1: 'topography_0p25.bin', etopo2022: 'topography2022_0p25.bin' };
export const ICE_SHEETS = ['rule', 'observed'], ICE_SHELVES = ['land', 'sea', 'cavity'];
export const LEGACY_GEOGRAPHY = { topography: 'etopo1', iceSheets: 'rule', iceShelves: 'land' };
export const FRESH_GEOGRAPHY = { topography: 'etopo2022', iceSheets: 'observed', iceShelves: 'land' };
export const topographyUrl = (name = 'etopo1') => new URL(`../data/${TOPOGRAPHIES[name]}`, import.meta.url);

export function checkGeography(record) {
  const { topography, iceSheets, iceShelves } = { ...LEGACY_GEOGRAPHY, ...record };
  if (!TOPOGRAPHIES[topography]) throw new Error(`topography is ${Object.keys(TOPOGRAPHIES).join(' or ')}, not ${topography}`);
  if (!ICE_SHEETS.includes(iceSheets)) throw new Error(`iceSheets is ${ICE_SHEETS.join(' or ')}, not ${iceSheets}`);
  if (!ICE_SHELVES.includes(iceShelves)) throw new Error(`iceShelves is ${ICE_SHELVES.join(' or ')}, not ${iceShelves}`);
  if (iceShelves !== 'land' && iceSheets !== 'observed') throw new Error(`ice shelves at ${iceShelves} need the observed ice sheets`);
  return { topography, iceSheets, iceShelves };
}

export const savedGeography = (saved) => checkGeography(Object.fromEntries(Object.keys(LEGACY_GEOGRAPHY).map((key) => [key, saved?.[key] ?? LEGACY_GEOGRAPHY[key]])));
export const geographyHeader = (record) => Object.fromEntries(Object.entries(checkGeography(record)).filter(([key, value]) => value !== LEGACY_GEOGRAPHY[key]));
export const pickGeography = (source) => geographyHeader(savedGeography(source));
export const geographyRecordOf = (topography, options = {}) => (topography?.name ? checkGeography({ topography: topography.name, iceSheets: options.landIce ? 'observed' : 'rule', iceShelves: options.iceShelves ?? 'land' }) : null);

export function readTopography(name = 'etopo1') {
  const fs = globalThis.process?.getBuiltinModule?.('node:fs');
  if (!fs) throw new Error(`data/${TOPOGRAPHIES[name]} is read under node only`);
  const bytes = fs.readFileSync(topographyUrl(name));
  return topographyFromInt16(bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength), 720, 1440, name);
}

/*
 * Under node, what a model of a saved state (or of a fresh start, with
 * null) needs of its geography: the record, the raster, and the
 * geography and land options.
 */
export function stateGeography(saved) {
  const record = saved ? savedGeography(saved) : { ...FRESH_GEOGRAPHY };
  return { record, topography: readTopography(record.topography), ...geographyOptions(record) };
}

/*
 * data/landice_0p25.bin, which scripts/landIce.mjs writes from ETOPO
 * 2022 on the topography's grid: a 24-byte header ('LIC3', rows,
 * columns, raster points per pixel, the n pixels holding ice, and the
 * firn air the flotation took off the floating ice's surface, float32
 * meters), then for those pixels their indices (uint32), the mean
 * thickness of their ice points and of those afloat (ice equivalent),
 * the mean surface of those afloat and the mean bed under them, in
 * meters (uint16, uint16, uint16, int16), and the counts of their ice
 * points and of those afloat (uint8, uint8).
 */
export const ICE_DENSITY = 917, SEAWATER_DENSITY = 1028;
export const landIceUrl = new URL('../data/landice_0p25.bin', import.meta.url);

export function encodeLandIce({ rows, cols, points, firnAir = 0, index, count, floating, thickness, floatThickness, floatSurface, floatBed }) {
  const n = index.length, out = new ArrayBuffer(24 + 14 * n), head = new DataView(out);
  [76, 73, 67, 51].forEach((b, j) => head.setUint8(j, b));
  head.setUint32(4, rows, true); head.setUint32(8, cols, true); head.setUint32(12, points, true); head.setUint32(16, n, true); head.setFloat32(20, firnAir, true);
  new Uint32Array(out, 24, n).set(index);
  new Uint16Array(out, 24 + 4 * n, n).set(thickness);
  new Uint16Array(out, 24 + 6 * n, n).set(floatThickness);
  new Uint16Array(out, 24 + 8 * n, n).set(floatSurface);
  new Int16Array(out, 24 + 10 * n, n).set(floatBed);
  new Uint8Array(out, 24 + 12 * n, n).set(count);
  new Uint8Array(out, 24 + 13 * n, n).set(floating);
  return out;
}

export function decodeLandIce(buffer) {
  const bytes = buffer instanceof ArrayBuffer ? buffer : buffer.buffer.slice(buffer.byteOffset, buffer.byteOffset + buffer.byteLength);
  const head = new DataView(bytes);
  if (bytes.byteLength < 24 || String.fromCharCode(...new Uint8Array(bytes, 0, 4)) !== 'LIC3') throw new Error('not a land-ice file');
  const rows = head.getUint32(4, true), cols = head.getUint32(8, true), points = head.getUint32(12, true), n = head.getUint32(16, true), firnAir = head.getFloat32(20, true);
  if (bytes.byteLength !== 24 + 14 * n) throw new Error('land-ice file of the wrong size');
  return {
    rows, cols, points, firnAir, index: new Uint32Array(bytes, 24, n), thickness: new Uint16Array(bytes, 24 + 4 * n, n), floatThickness: new Uint16Array(bytes, 24 + 6 * n, n),
    floatSurface: new Uint16Array(bytes, 24 + 8 * n, n), floatBed: new Int16Array(bytes, 24 + 10 * n, n), count: new Uint8Array(bytes, 24 + 12 * n, n), floating: new Uint8Array(bytes, 24 + 13 * n, n),
  };
}

/*
 * The geography and land options of a geography record: the observed
 * ice sheets take the ice of data/landice_0p25.bin (`landIce` when
 * given, else read from the file under node) and the land's
 * iceSheetBalance.
 */
export function geographyOptions(record = {}, landIce = null) {
  const { iceSheets, iceShelves } = checkGeography(record);
  if (iceSheets === 'rule') return { geography: {}, land: {} };
  if (!landIce) {
    const fs = globalThis.process?.getBuiltinModule?.('node:fs');
    if (!fs) throw new Error('the observed ice sheets need the land ice of data/landice_0p25.bin');
    landIce = decodeLandIce(fs.readFileSync(landIceUrl));
  }
  return { geography: { landIce, ...(iceShelves !== 'land' ? { iceShelves } : {}) }, land: { iceSheetBalance: true } };
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

export function createGeography(mesh, topography, { landThreshold = 0.5, landBridges = LAND_BRIDGES, seaStraits = SEA_STRAITS, landIce = null, iceShelves = 'land', cavityDepth = 50 } = {}) {
  const { nCells: C, nEdges: E, xCell, cellsOnCell, nEdgesOnCell, edgesOnCell, cellsOnEdge, areaCell, latCell, lonCell } = mesh;
  const { rows, cols } = topography;
  if (landIce && (landIce.rows !== rows || landIce.cols !== cols)) throw new Error(`the land ice is ${landIce.rows} × ${landIce.cols}, the topography ${rows} × ${cols}`);
  if (!ICE_SHELVES.includes(iceShelves)) throw new Error(`iceShelves is ${ICE_SHELVES.join(' or ')}, not ${iceShelves}`);
  if (iceShelves !== 'land' && !landIce) throw new Error(`ice shelves at ${iceShelves} need the land ice`);
  if (iceShelves === 'cavity') {
    const options = { landThreshold, landBridges, seaStraits, landIce };
    return cavityGeography(mesh, createGeography(mesh, topography, options), createGeography(mesh, topography, { ...options, iceShelves: 'sea' }), cavityDepth);
  }
  const slot = landIce ? new Int32Array(rows * cols).fill(-1) : null, afloat = iceShelves === 'sea';
  if (landIce) landIce.index.forEach((n, k) => { slot[n] = k; });
  const data = !afloat ? topography.data : Float64Array.from(topography.data, (value, n) => {
    const k = slot[n];
    if (k < 0 || !landIce.floating[k]) return value;
    return value + landIce.floating[k] / landIce.points * (landIce.floatBed[k] - landIce.floatSurface[k]);
  });
  const count = new Int32Array(C), landCount = new Int32Array(C), sum = new Float64Array(C), owner = new Int32Array(rows * cols);
  const iceShare = landIce ? new Float64Array(C) : null, icePoints = landIce ? new Float64Array(C) : null, iceVolume = landIce ? new Float64Array(C) : null, shelfShare = landIce ? new Float64Array(C) : null, shelfVolume = landIce ? new Float64Array(C) : null, wetSum = landIce ? new Float64Array(C) : null;
  const addIce = (i, n) => {
    const k = slot[n];
    if (k < 0) return;
    const { count: ice, floating, thickness, floatThickness } = landIce, held = afloat ? ice[k] - floating[k] : ice[k];
    iceShare[i] += held / landIce.points; icePoints[i] += held; shelfShare[i] += floating[k] / landIce.points;
    shelfVolume[i] += floating[k] * floatThickness[k] / landIce.points;
    iceVolume[i] += Math.max(0, ice[k] * thickness[k] - (afloat ? floating[k] * floatThickness[k] : 0));
  };
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
      owner[r * cols + c] = cell;
      count[cell]++;
      sum[cell] += elevation;
      if (elevation > 0) landCount[cell]++;
      if (landIce) { addIce(cell, r * cols + c); wetSum[cell] += Math.max(0, -elevation); }
    }
  }
  const land = new Uint8Array(C), landFraction = new Float64Array(C), elevation = new Float64Array(C);
  for (let i = 0; i < C; i++) {
    if (count[i] === 0) {
      const r = Math.min(rows - 1, Math.max(0, Math.floor((Math.PI / 2 - latCell[i]) / (Math.PI / rows))));
      const c = ((Math.floor((lonCell[i] + Math.PI) / (2 * Math.PI / cols)) % cols) + cols) % cols;
      const value = data[r * cols + c];
      count[i] = 1; sum[i] = value; landCount[i] = value > 0 ? 1 : 0;
      if (landIce) { addIce(i, r * cols + c); wetSum[i] = Math.max(0, -value); }
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
  const edgeSill = new Float64Array(E);
  const pointStart = new Int32Array(C + 1), pointList = new Int32Array(rows * cols);
  for (let n = 0; n < rows * cols; n++) pointStart[owner[n] + 1]++;
  for (let i = 0; i < C; i++) pointStart[i + 1] += pointStart[i];
  {
    const fill = Int32Array.from(pointStart.subarray(0, C));
    for (let n = 0; n < rows * cols; n++) pointList[fill[owner[n]]++] = n;
  }
  const local = new Int32Array(rows * cols).fill(-1), parent = new Int32Array(rows * cols);
  const find = (p) => { while (parent[p] !== p) { parent[p] = parent[parent[p]]; p = parent[p]; } return p; };
  const depthOf = (n) => Math.max(0, -data[n]);
  const depthAt = (x, y, z) => {
    const r = Math.min(rows - 1, Math.max(0, Math.floor((Math.PI / 2 - Math.asin(Math.max(-1, Math.min(1, z)))) / (Math.PI / rows))));
    const c = ((Math.floor((Math.atan2(y, x) + Math.PI) / (2 * Math.PI / cols)) % cols) + cols) % cols;
    return depthOf(r * cols + c);
  };
  const sampleStep = 0.1 * Math.PI / 180;
  const alongCenters = (a, b) => {
    const ax = xCell[3 * a], ay = xCell[3 * a + 1], az = xCell[3 * a + 2], bx = xCell[3 * b], by = xCell[3 * b + 1], bz = xCell[3 * b + 2];
    const theta = Math.acos(Math.max(-1, Math.min(1, ax * bx + ay * by + az * bz))), sinTheta = Math.sin(theta), n = Math.max(1, Math.ceil(theta / sampleStep));
    let sill = Infinity;
    for (let s = 0; s <= n && sill > 0; s++) {
      const t = s / n, wa = sinTheta > 0 ? Math.sin((1 - t) * theta) / sinTheta : 1 - t, wb = sinTheta > 0 ? Math.sin(t * theta) / sinTheta : t;
      sill = Math.min(sill, depthAt(wa * ax + wb * bx, wa * ay + wb * by, wa * az + wb * bz));
    }
    return sill;
  };
  const deepestOf = (i) => { let best = -1, depth = -1; for (let p = pointStart[i]; p < pointStart[i + 1]; p++) { const d = depthOf(pointList[p]); if (d > depth) { depth = d; best = pointList[p]; } } return best; };
  for (let e = 0; e < E; e++) {
    if (!edgeOcean[e]) continue;
    const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1], da = deepestOf(a), db = deepestOf(b);
    if (da < 0 || db < 0) { edgeSill[e] = alongCenters(a, b); continue; }
    const points = [];
    for (const i of [a, b]) for (let p = pointStart[i]; p < pointStart[i + 1]; p++) points.push(pointList[p]);
    points.sort((m, n) => depthOf(n) - depthOf(m));
    let sill = -1;
    for (let q = 0; q < points.length; q++) {
      const n = points[q], r = Math.floor(n / cols), c = n - r * cols;
      local[n] = q; parent[n] = n;
      for (let dr = -1; dr <= 1; dr++) {
        const rr = r + dr;
        if (rr < 0 || rr >= rows) continue;
        for (let dc = -1; dc <= 1; dc++) {
          const m = rr * cols + ((c + dc + cols) % cols);
          if (m === n || local[m] < 0) continue;
          const x = find(n), y = find(m);
          if (x !== y) parent[x] = y;
        }
      }
      if (local[da] >= 0 && local[db] >= 0 && find(da) === find(db)) { sill = depthOf(n); break; }
    }
    for (const n of points) local[n] = -1;
    edgeSill[e] = sill >= 0 ? sill : alongCenters(a, b);
  }
  for (const { sill, points } of Object.values(seaStraits)) {
    let previous = -1;
    trace(points, (i) => {
      if (previous >= 0 && previous !== i) {
        for (let m = 0; m < nEdgesOnCell[previous]; m++) {
          if (cellsOnCell[mesh.maxEdges * previous + m] !== i) continue;
          const e = edgesOnCell[mesh.maxEdges * previous + m];
          if (edgeOcean[e]) edgeSill[e] = Math.max(edgeSill[e], sill);
        }
      }
      previous = i;
    });
  }
  const iceSheet = new Uint8Array(C), iceFraction = landIce ? new Float64Array(C) : null, iceThickness = landIce ? new Float64Array(C) : null, shelfFraction = landIce ? new Float64Array(C) : null, shelfThickness = landIce ? new Float64Array(C) : null, wetDepth = landIce ? new Float64Array(C) : null;
  for (let i = 0; i < C; i++) {
    if (landIce) { shelfFraction[i] = shelfShare[i] / count[i]; shelfThickness[i] = shelfVolume[i] / count[i]; wetDepth[i] = wetSum[i] / count[i]; }
    if (!land[i]) continue;
    if (landIce) {
      iceFraction[i] = iceShare[i] / count[i];
      if (iceFraction[i] > 0.5 * landFraction[i]) { iceSheet[i] = 1; iceThickness[i] = iceVolume[i] / icePoints[i]; }
      continue;
    }
    const lat = latCell[i] * 180 / Math.PI, lon = lonCell[i] * 180 / Math.PI;
    if (lat < -60 || (lat > 60 && lat < 84 && lon > -73 && lon < -12 && elevation[i] > 800)) iceSheet[i] = 1;
  }
  let landArea = 0, total = 0;
  for (let i = 0; i < C; i++) { total += areaCell[i]; if (land[i]) landArea += areaCell[i]; }
  return { land, landFraction, elevation, iceSheet, ...(landIce ? { iceFraction, iceThickness, shelfFraction, shelfThickness, wetDepth } : {}), edgeOcean, edgeSill, coastEdges: Int32Array.from(coast), landArea: landArea / total };
}

function cavityGeography(mesh, shelf, sea, cavityDepth) {
  const { nCells: C, nEdges: E, cellsOnEdge } = mesh;
  const cavity = new Uint8Array(C), iceDraft = new Float64Array(C), land = Uint8Array.from(shelf.land), elevation = Float64Array.from(sea.elevation);
  for (let i = 0; i < C; i++) {
    if (!shelf.land[i] || sea.land[i]) continue;
    const draft = shelf.shelfThickness[i] * ICE_DENSITY / SEAWATER_DENSITY;
    if (sea.wetDepth[i] - draft < cavityDepth) continue;
    cavity[i] = 1; iceDraft[i] = draft; land[i] = 0; elevation[i] = -sea.wetDepth[i];
  }
  for (let i = 0; i < C; i++) if (land[i]) elevation[i] = shelf.elevation[i];
  const edgeOcean = new Uint8Array(E), edgeSill = new Float64Array(E), coast = [];
  for (let e = 0; e < E; e++) {
    const a = land[cellsOnEdge[2 * e]], b = land[cellsOnEdge[2 * e + 1]];
    edgeOcean[e] = !a && !b ? 1 : 0;
    if (edgeOcean[e]) edgeSill[e] = sea.edgeSill[e];
    if (a !== b) coast.push(e);
  }
  return { ...shelf, cavity, iceDraft, ocean: { land, elevation, edgeOcean, edgeSill, coastEdges: Int32Array.from(coast) } };
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
 * ('SGO2', the cell count, the field count, the Lloyd passes of the mesh
 * the fields were computed on) and per field one uint16 per cell,
 * value = q × scale + offset, in SUBGRID_FIELDS' order. Fields without an
 * integer `relax` are written as 'SGO1' with 0 in its place, and 'SGO1'
 * files decode with relax null.
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
  const relax = Number.isInteger(fields.relax) && fields.relax >= 0 ? fields.relax : null;
  [83, 71, 79, relax === null ? 49 : 50].forEach((b, j) => head.setUint8(j, b));
  head.setUint32(4, C, true); head.setUint32(8, SUBGRID_FIELDS.length, true); head.setUint32(12, relax ?? 0, true);
  SUBGRID_FIELDS.forEach(([name, scale, offset], f) => {
    const q = new Uint16Array(out, 16 + 2 * C * f, C), values = fields[name];
    for (let i = 0; i < C; i++) q[i] = Math.min(65535, Math.max(0, Math.round((values[i] - offset) / scale)));
  });
  return out;
}

export function decodeSubgrid(buffer) {
  const bytes = buffer instanceof ArrayBuffer ? buffer : buffer.buffer.slice(buffer.byteOffset, buffer.byteOffset + buffer.byteLength);
  const head = new DataView(bytes);
  const magic = String.fromCharCode(...new Uint8Array(bytes, 0, 4));
  if (magic !== 'SGO1' && magic !== 'SGO2') throw new Error('not a subgrid orography file');
  const C = head.getUint32(4, true), count = head.getUint32(8, true);
  if (count !== SUBGRID_FIELDS.length || bytes.byteLength !== 16 + 2 * C * count) throw new Error('subgrid orography file of the wrong size');
  const fields = { relax: magic === 'SGO2' ? head.getUint32(12, true) : null };
  SUBGRID_FIELDS.forEach(([name, scale, offset], f) => {
    const q = new Uint16Array(bytes, 16 + 2 * C * f, C);
    fields[name] = Float64Array.from(q, (v) => v * scale + offset);
  });
  for (let i = 0; i < C; i++) if (!(fields.deviation[i] > 0)) fields.orientation[i] = 0;
  return fields;
}

export const subgridUrl = (N, topography = 'etopo1') => new URL(`../data/subgrid_${topography === 'etopo1' ? '' : `${topography}_`}N${N}.bin`, import.meta.url);

/*
 * The fields for a mesh: `given` when it has the mesh's cell count, else
 * under node the file for N and the topography if there is one, else null
 * (the caller computes them from the 0.25° raster and says so).
 */
export function meshSubgrid(mesh, given = undefined, topography = 'etopo1') {
  if (given === false) return null;
  if (given) return given.deviation.length === mesh.nCells ? given : null;
  const N = Math.round(Math.sqrt((mesh.nCells - 2) / 10));
  const fs = globalThis.process?.getBuiltinModule?.('node:fs');
  if (!fs || 10 * N * N + 2 !== mesh.nCells) return null;
  const url = subgridUrl(N, topography);
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
  missing: (file) => `no ${file} for this mesh`,
  off: () => 'the subgrid files are turned off (subgrid false)',
  terrain: () => 'the terrain is off, and the files hold only the scales below the resolved terrain',
  land: (file) => `${file} was made for another land mask`,
  relax: (file, mesh, fields) => `${file} was made for the mesh of ${fields.relax} Lloyd passes, not this one of ${mesh.relax}`,
};
const fallbackSaid = new Set();
export function subgridFallbackNote(mesh, reason = 'missing', fields = null, topography = 'etopo1') {
  const N = Math.round(Math.sqrt((mesh.nCells - 2) / 10)), file = `data/${subgridUrl(N, topography).pathname.split('/').pop()}`;
  if (fallbackSaid.has(`${file} ${reason}`)) return;
  fallbackSaid.add(`${file} ${reason}`);
  console.log(`subgrid orography: ${SUBGRID_FALLBACK[reason](file, mesh, fields)}; computed from the topography's raster (at 0.25°, σ about a third of the 5 km fields') with Lott and Miller's constants, and no form drag`);
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
