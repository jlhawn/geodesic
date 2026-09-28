/*
 * Land and elevation on the mesh from a regular latitude–longitude
 * elevation raster: rows from north to south, columns from west to
 * east starting at −180°, each value the elevation in metres, negative
 * under the sea. Every raster point is assigned to the nearest cell by
 * walking the cell neighbour graph from the previous point's cell; a
 * cell's elevation is the mean of its points and its land fraction the
 * share of them above sea level. A cell is land when that fraction
 * exceeds landThreshold, and along the polylines of landBridges (the
 * Panama isthmus by default) regardless, since a coarse mesh would
 * otherwise open a seaway where a narrow isthmus falls below the
 * threshold. edgeOcean marks the edges between two ocean
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

export function createGeography(mesh, topography, { landThreshold = 0.5, landBridges = LAND_BRIDGES } = {}) {
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
  for (const points of Object.values(landBridges)) {
    let start = 0;
    for (let p = 0; p + 1 < points.length; p++) {
      const [la0, lo0] = points[p], [la1, lo1] = points[p + 1];
      const steps = Math.max(1, Math.ceil(Math.hypot(la1 - la0, (lo1 - lo0) * Math.cos(la0 * Math.PI / 180)) / 0.2));
      for (let n = 0; n <= steps; n++) {
        const la = (la0 + (la1 - la0) * n / steps) * Math.PI / 180, lo = (lo0 + (lo1 - lo0) * n / steps) * Math.PI / 180;
        const i = nearest(Math.cos(la) * Math.cos(lo), Math.cos(la) * Math.sin(lo), Math.sin(la), start);
        start = i; land[i] = 1; elevation[i] = Math.max(elevation[i], 1);
      }
    }
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
