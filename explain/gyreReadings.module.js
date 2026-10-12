/*
 * Readings of oceanCases' basin frames for the gyre figure. Speeds are means
 * from the surface to the case's flat 4000 m floor, in m/s, positive north or
 * east; widths are meters along 30°N; transports are Sv (10⁶ m³/s). Each cell
 * of the frame's 30°N row is as wide as half the gaps to its neighbors, the
 * rule the frame's own Sv totals use.
 */
import { EARTH_RADIUS } from '../js/mesh.module.js';

export const DEPTH = 4000, ROW = 30;
const DEG = Math.PI / 180, PER_DEGREE = EARTH_RADIUS * DEG * Math.cos(ROW * DEG);

export function rowWidths(longitude) {
  const n = longitude.length;
  return Float64Array.from(longitude, (lon, s) => {
    const left = s > 0 ? 0.5 * (lon - longitude[s - 1]) : 0.5 * (longitude[s + 1] - lon);
    const right = s < n - 1 ? 0.5 * (longitude[s + 1] - lon) : left;
    return left + right;
  });
}

export function westernCurrent({ longitude, north }) {
  const n = longitude.length;
  if (n < 2 || !(north[0] > 0)) return null;
  const widths = rowWidths(longitude), coast = longitude[0] - 0.5 * (longitude[1] - longitude[0]);
  let k = 0, peak = north[0], transport = north[0] * widths[0];
  while (k + 1 < n && north[k + 1] > 0) { k++; peak = Math.max(peak, north[k]); transport += north[k] * widths[k]; }
  const end = k + 1 < n ? longitude[k] + (longitude[k + 1] - longitude[k]) * north[k] / (north[k] - north[k + 1]) : longitude[n - 1] + 0.5 * widths[n - 1];
  return { coast, end, width: (end - coast) * PER_DEGREE, peak: peak / DEPTH, transport: transport * PER_DEGREE / 1e6 };
}

export function zonalExtremes({ latitude, east }, south = 10, north = 50) {
  let westward = null, eastward = null;
  for (let b = 0; b < latitude.length; b++) {
    if (latitude[b] < south || latitude[b] > north || !Number.isFinite(east[b])) continue;
    const speed = east[b] / DEPTH;
    if (!westward || speed < westward.speed) westward = { speed, latitude: latitude[b] };
    if (!eastward || speed > eastward.speed) eastward = { speed, latitude: latitude[b] };
  }
  return { westward, eastward };
}

export function gyreReadings(frame) {
  return {
    day: frame.day,
    west: frame.west / DEPTH,
    interior: frame.interior / DEPTH,
    east: frame.east / DEPTH,
    ratio: frame.interior < 0 && frame.west > 0 ? frame.west / -frame.interior : null,
    current: frame.walls ? westernCurrent(frame.parallel) : null,
    sverdrup: frame.sverdrup === null ? null : { speed: frame.sverdrup / DEPTH, transport: frame.sverdrupTotal },
    across: { north: frame.northward, south: frame.southward },
    zonal: zonalExtremes(frame.zonal),
  };
}
