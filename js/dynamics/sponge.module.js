/*
 * The top sponge on the zonally asymmetric wind (Shepherd, Semeniuk &
 * Koshyk 1996): in each sponge layer the departure of the normal velocity
 * from its zonal mean relaxes at the layer's rate, and the zonal mean is
 * left alone. The zonal mean of the cells' east and north wind (Section
 * 3.7's reconstruction) in latitude bands as wide as the mesh spacing,
 * the area-weighted sum over the band's area-weighted cosine of
 * latitude, is interpolated linearly in latitude between the bands' mean
 * latitudes to each edge, times the cosine of the edge's latitude, so
 * that solid-body rotation is kept to the poles, and projected onto the
 * edge's normal. The cell on a pole, whose east is undefined, is left
 * out of its band.
 */
export function spongeGeometry(mesh) {
  const { nCells: C, nEdges: E, maxEdges, nEdgesOnCell, edgesOnCell, dcEdge, dvEdge, nEdge, areaCell, latCell, lonCell, xEdge } = mesh;
  let spacing = 0;
  for (let e = 0; e < E; e++) spacing += dcEdge[e];
  spacing /= E * mesh.radius;
  const bands = Math.max(2, Math.round(Math.PI / spacing)), width = Math.PI / bands;
  const cellBand = Int32Array.from(latCell, (lat) => Math.min(bands - 1, Math.floor((lat + Math.PI / 2) / width)));
  const bandArea = new Float64Array(bands), bandCosine = new Float64Array(bands), bandLat = new Float64Array(bands);
  const onPole = (i) => Math.abs(latCell[i]) > Math.PI / 2 - 1e-9;
  for (let i = 0; i < C; i++) {
    if (onPole(i)) continue;
    const b = cellBand[i];
    bandArea[b] += areaCell[i]; bandCosine[b] += areaCell[i] * Math.cos(latCell[i]); bandLat[b] += areaCell[i] * latCell[i];
  }
  for (let b = 0; b < bands; b++) bandLat[b] = bandArea[b] > 0 ? bandLat[b] / bandArea[b] : -Math.PI / 2 + (b + 0.5) * width;
  const eastWeight = new Float64Array(maxEdges * C), northWeight = new Float64Array(maxEdges * C);
  for (let i = 0; i < C; i++) {
    if (onPole(i)) continue;
    const lat = latCell[i], lon = lonCell[i];
    const east = [-Math.sin(lon), Math.cos(lon), 0], north = [-Math.sin(lat) * Math.cos(lon), -Math.sin(lat) * Math.sin(lon), Math.cos(lat)];
    for (let m = 0; m < nEdgesOnCell[i]; m++) {
      const e = edgesOnCell[maxEdges * i + m], b = cellBand[i], w = bandCosine[b] > 0 ? 0.5 * dcEdge[e] * dvEdge[e] / bandCosine[b] : 0;
      eastWeight[maxEdges * i + m] = w * (nEdge[3 * e] * east[0] + nEdge[3 * e + 1] * east[1]);
      northWeight[maxEdges * i + m] = w * (nEdge[3 * e] * north[0] + nEdge[3 * e + 1] * north[1] + nEdge[3 * e + 2] * north[2]);
    }
  }
  const edgeBand = new Int32Array(E), edgeShare = new Float64Array(E), edgeEast = new Float64Array(E), edgeNorth = new Float64Array(E);
  for (let e = 0; e < E; e++) {
    const x = xEdge[3 * e], y = xEdge[3 * e + 1], z = xEdge[3 * e + 2], r = Math.hypot(x, y, z), lat = Math.asin(z / r), lon = Math.atan2(y, x);
    let lower = Math.max(0, Math.min(bands - 2, Math.floor((lat + Math.PI / 2) / width - 0.5)));
    while (lower > 0 && bandLat[lower] > lat) lower--;
    while (lower < bands - 2 && bandLat[lower + 1] <= lat) lower++;
    edgeBand[e] = lower;
    edgeShare[e] = Math.max(0, Math.min(1, (lat - bandLat[lower]) / (bandLat[lower + 1] - bandLat[lower])));
    edgeEast[e] = Math.cos(lat) * (-Math.sin(lon) * nEdge[3 * e] + Math.cos(lon) * nEdge[3 * e + 1]);
    edgeNorth[e] = Math.cos(lat) * (-Math.sin(lat) * Math.cos(lon) * nEdge[3 * e] - Math.sin(lat) * Math.sin(lon) * nEdge[3 * e + 1] + Math.cos(lat) * nEdge[3 * e + 2]);
  }
  const order = Int32Array.from({ length: C }, (_, i) => i).sort((a, b) => cellBand[a] - cellBand[b] || a - b);
  const bandStart = new Int32Array(bands + 1);
  for (let i = 0; i < C; i++) bandStart[cellBand[i] + 1]++;
  for (let b = 0; b < bands; b++) bandStart[b + 1] += bandStart[b];
  return { bands, cellBand, eastWeight, northWeight, edgeBand, edgeShare, edgeEast, edgeNorth, order, bandStart };
}

export const SPONGE = { sigma: 0.005, days: 1 };

/*
 * The sponge's onset on a grid by name: on bl36 the 78 Pa (over p0) at
 * which the IFS's begins (IFS CY48r1 Part III, 2.2.11), so that the
 * 0.3-1 hPa layer keeps its equatorial waves; SPONGE.sigma elsewhere.
 */
export function spongeSigmaFor(gridName) {
  return gridName === 'bl36' ? 78 / 101325 : SPONGE.sigma;
}

/*
 * The sponge's rate (1/s) in each layer: 1/days at the top, falling
 * linearly in σ to zero at `sigma`.
 */
export function spongeRates(sigmaMid, sigma, days) {
  return Float64Array.from(sigmaMid, (s) => (days > 0 && s < sigma ? (sigma - s) / sigma / (days * 86400) : 0));
}

/*
 * One step of the sponge on one layer's normal velocity, in place:
 * u ← ū + (u − ū)/(1 + r dt), the zonal mean ū from `means` (2 per band,
 * filled here).
 */
export function dampEddies(mesh, geometry, u, rate, dt, means) {
  const { nCells: C, nEdges: E, maxEdges, nEdgesOnCell, edgesOnCell } = mesh;
  const { bands, cellBand, eastWeight, northWeight, edgeBand, edgeShare, edgeEast, edgeNorth } = geometry;
  means.fill(0, 0, 2 * bands);
  for (let i = 0; i < C; i++) {
    let east = 0, north = 0;
    for (let m = 0; m < nEdgesOnCell[i]; m++) {
      const slot = maxEdges * i + m, ue = u[edgesOnCell[slot]];
      east += eastWeight[slot] * ue;
      north += northWeight[slot] * ue;
    }
    means[2 * cellBand[i]] += east;
    means[2 * cellBand[i] + 1] += north;
  }
  const keep = 1 / (1 + rate * dt);
  for (let e = 0; e < E; e++) {
    const b = edgeBand[e], s = edgeShare[e];
    const east = (1 - s) * means[2 * b] + s * means[2 * b + 2], north = (1 - s) * means[2 * b + 1] + s * means[2 * b + 3];
    const mean = east * edgeEast[e] + north * edgeNorth[e];
    u[e] = mean + (u[e] - mean) * keep;
  }
}
