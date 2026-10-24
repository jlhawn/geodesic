/*
 * Lakes on a latitude–longitude raster (rows from the north, columns from
 * −180°) from GeoJSON in the form of Natural Earth's ne_10m_lakes:
 * features with properties featurecla, name and scalerank, and Polygon or
 * MultiPolygon geometry. Features whose featurecla names an intermittent
 * water body are left out, as are any whose `intermittent` property is
 * set, and with `pans` those it holds, whose featurecla and the center of
 * whose bounds it matches. A point belongs to a lake where its center lies inside an odd
 * number of the lake's rings.
 */
export const INTERMITTENT = /playa|intermittent|ephemeral/i;

export const AUSTRALIAN_PANS = { featurecla: 'Alkaline Lake', box: [110, 155, -45, -10] };

const center = (polygons) => {
  let west = 180, east = -180, south = 90, north = -90;
  for (const rings of polygons) for (const ring of rings) for (const [lon, lat] of ring) { west = Math.min(west, lon); east = Math.max(east, lon); south = Math.min(south, lat); north = Math.max(north, lat); }
  return [(west + east) / 2, (south + north) / 2];
};

export function lakeFeatures(geojson, { exclude = INTERMITTENT, pans = null } = {}) {
  const lakes = [], left = [];
  for (const feature of geojson.features ?? []) {
    const { properties = {}, geometry } = feature;
    if (!geometry || !['Polygon', 'MultiPolygon'].includes(geometry.type)) continue;
    const record = { name: properties.name ?? null, featurecla: properties.featurecla ?? null, scalerank: properties.scalerank ?? null };
    if ((record.featurecla && exclude.test(record.featurecla)) || properties.intermittent) { left.push(record); continue; }
    const polygons = geometry.type === 'Polygon' ? [geometry.coordinates] : geometry.coordinates;
    if (pans && record.featurecla === pans.featurecla) {
      const [lon, lat] = center(polygons), [w, e, s, n] = pans.box;
      if (lon >= w && lon <= e && lat >= s && lat <= n) { left.push({ ...record, pan: true }); continue; }
    }
    lakes.push({ ...record, polygons });
  }
  return { lakes, left };
}

export function rasterizeLakes(lakes, rows, cols) {
  if (lakes.length >= 65535) throw new Error(`${lakes.length} lakes, more than a uint16 raster holds`);
  const id = new Uint16Array(rows * cols), dlat = 180 / rows, dlon = 360 / cols, crossings = [];
  lakes.forEach((lake, k) => {
    for (const rings of lake.polygons) {
      let south = 90, north = -90;
      for (const ring of rings) for (const [, lat] of ring) { south = Math.min(south, lat); north = Math.max(north, lat); }
      const r0 = Math.max(0, Math.floor((90 - north) / dlat - 0.5)), r1 = Math.min(rows - 1, Math.ceil((90 - south) / dlat - 0.5));
      for (let r = r0; r <= r1; r++) {
        const lat = 90 - (r + 0.5) * dlat;
        crossings.length = 0;
        for (const ring of rings) {
          for (let p = 0; p + 1 < ring.length; p++) {
            const [x0, y0] = ring[p], [x1, y1] = ring[p + 1];
            if ((y0 > lat) !== (y1 > lat)) crossings.push(x0 + (lat - y0) * (x1 - x0) / (y1 - y0));
          }
        }
        crossings.sort((a, b) => a - b);
        for (let q = 0; q + 1 < crossings.length; q += 2) {
          const c0 = Math.max(0, Math.ceil((crossings[q] + 180) / dlon - 0.5)), c1 = Math.min(cols - 1, Math.ceil((crossings[q + 1] + 180) / dlon - 0.5) - 1);
          for (let c = c0; c <= c1; c++) id[r * cols + c] = k + 1;
        }
      }
    }
  });
  return id;
}

/*
 * Each lake's level from the raster's surface (`surface`, meters, on the
 * same grid): where at least `flat` of the lake's points share one value
 * to the meter, the raster holds the water's surface and that is the
 * level; otherwise, the raster holding the lake's bed, the level is the
 * `low` quantile of the surface on the points just outside the lake (its
 * four-neighbors in no lake of its own), the shore at the water's edge.
 */
export function lakeLevels(id, surface, rows, cols, count, { flat = 0.3, low = 0.1 } = {}) {
  const inside = Array.from({ length: count }, () => new Map()), shore = Array.from({ length: count }, () => []);
  for (let r = 0; r < rows; r++) {
    for (let c = 0; c < cols; c++) {
      const n = r * cols + c, k = id[n];
      if (!k) continue;
      const v = Math.round(surface[n]), m = inside[k - 1];
      m.set(v, (m.get(v) ?? 0) + 1);
      for (const [rr, cc] of [[r - 1, c], [r + 1, c], [r, c - 1], [r, c + 1]]) {
        if (rr < 0 || rr >= rows) continue;
        const q = rr * cols + ((cc + cols) % cols);
        if (id[q] !== k) shore[k - 1].push(surface[q]);
      }
    }
  }
  return inside.map((m, k) => {
    let points = 0, mode = 0, most = 0;
    for (const [v, n] of m) { points += n; if (n > most) { most = n; mode = v; } }
    if (!points) return { level: NaN, points: 0, flat: false };
    if (most >= flat * points) return { level: mode, points, flat: true };
    const s = Float64Array.from(shore[k]).sort();
    return { level: s.length ? s[Math.min(s.length - 1, Math.floor(low * s.length))] : mode, points, flat: false };
  });
}
