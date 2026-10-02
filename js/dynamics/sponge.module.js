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
 * Rayleigh friction on the zonal-mean wind of the top layers, standing in
 * for the mesospheric gravity-wave drag above a low lid: profiles of the decay time (days) against
 * log-pressure height z = H ln(p0/p), H = 7 km, as pairs [z (m), days],
 * the rate's logarithm linear in z between them, the last held above
 * and none below the first. `holtonWehrbein`: Holton & Wehrbein (1980,
 * PAGEOPH 118, 284), 5 days at 65 km to 2 days at 75 km, as Rind et al.
 * (1984) quote them; `rind`: the GISS 21-layer model's drag in its layers
 * of 65-75 km (Rind, Suozzo, Lacis, Russell & Hansen 1984, NASA
 * TM-86183), its winter decay times 2, 1 and 0.5 days at 65, 70 and 75 km.
 */
export const LID_FRICTION = {
  holtonWehrbein: [[65e3, 5], [75e3, 2]],
  rind: [[65e3, 2], [70e3, 1], [75e3, 0.5]],
};
const FRICTION_SCALE_HEIGHT = 7e3;

/*
 * Each layer's friction rate (1/s): the profile's rate averaged over the
 * layer's mass at the reference surface pressure p0. `profile` is a name
 * in LID_FRICTION, an array of [z, days] pairs, or null for none.
 */
export function lidFrictionRates(levels, profile, p0 = 101325) {
  const K = levels.length - 1, rates = new Float64Array(K), points = typeof profile === 'string' ? LID_FRICTION[profile] : profile;
  if (typeof profile === 'string' && !points) throw new Error(`no lid friction is named ${profile}; the profiles are ${Object.keys(LID_FRICTION).join(', ')}`);
  if (!points || !points.length) return rates;
  const rateAt = (p) => {
    const z = FRICTION_SCALE_HEIGHT * Math.log(p0 / p);
    if (z < points[0][0]) return 0;
    let n = 0;
    while (n < points.length - 1 && points[n + 1][0] <= z) n++;
    if (n === points.length - 1) return 1 / (points[n][1] * 86400);
    const [z0, d0] = points[n], [z1, d1] = points[n + 1], t = (z - z0) / (z1 - z0);
    return Math.exp((1 - t) * Math.log(1 / d0) + t * Math.log(1 / d1)) / 86400;
  };
  const steps = 2000;
  for (let k = 0; k < K; k++) {
    const top = levels[k] * p0, bottom = levels[k + 1] * p0;
    let sum = 0;
    for (let n = 0; n < steps; n++) sum += rateAt(top + (n + 0.5) * (bottom - top) / steps);
    rates[k] = sum / steps;
  }
  return rates;
}

/*
 * The lid friction of a grid by name: none.
 */
export function lidFrictionFor() {
  return null;
}

/*
 * One step of the sponge on one layer's normal velocity, in place:
 * u ← ū/(1 + r̄ dt) + (u − ū)/(1 + r dt), the zonal mean ū from `means`
 * (2 per band, filled here), r the eddies' rate and r̄ the zonal mean's
 * (the lid friction's).
 */
export function dampEddies(mesh, geometry, u, rate, dt, means, meanRate = 0) {
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
  const keep = 1 / (1 + rate * dt), meanKeep = 1 / (1 + meanRate * dt);
  for (let e = 0; e < E; e++) {
    const b = edgeBand[e], s = edgeShare[e];
    const east = (1 - s) * means[2 * b] + s * means[2 * b + 2], north = (1 - s) * means[2 * b + 1] + s * means[2 * b + 3];
    const mean = east * edgeEast[e] + north * edgeNorth[e];
    u[e] = mean * meanKeep + (u[e] - mean) * keep;
  }
}
