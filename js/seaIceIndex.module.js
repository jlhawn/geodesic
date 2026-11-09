/*
 * The model's sea-ice extent and area as NSIDC's Sea Ice Index counts them: a day's concentration cut at 15 % over sea
 * cells (neither land nor lake) inside the index's polar stereographic grids, so closed basins off those grids (the
 * Caspian, the Black Sea and Azov) stay out as they do from the index. The grids (NSIDC-0051 and G02135 user guides):
 * north, true at 70°N about 45°W, x from −3850 to 3750 km and y from −5350 to 5850 km; south, true at 70°S about 0°,
 * x from −3950 to 3950 km and y from −3950 to 4350 km; on a sphere of the Hughes 1980 radius, a few km from the
 * ellipsoid's positions, which matters only for cells on a grid's edge, where no ice forms.
 */
const RADIUS = 6378.273, TRUE_LATITUDE = 70 * Math.PI / 180, SCALE = (1 + Math.sin(TRUE_LATITUDE)) / 2;
export const INDEX_GRIDS = {
  N: { x: [-3850, 3750], y: [-5350, 5850], central: -45 * Math.PI / 180 },
  S: { x: [-3950, 3950], y: [-3950, 4350], central: 0 },
};

export function inIndexGrid(lat, lon) {
  const north = lat > 0, g = INDEX_GRIDS[north ? 'N' : 'S'];
  const r = 2 * RADIUS * SCALE * Math.tan(Math.PI / 4 - Math.abs(lat) / 2), d = lon - g.central;
  const x = r * Math.sin(d), y = north ? -r * Math.cos(d) : r * Math.cos(d);
  return x >= g.x[0] && x <= g.x[1] && y >= g.y[0] && y <= g.y[1];
}

export function indexCells(mesh, geography) {
  const { land, lake } = geography;
  return Uint8Array.from({ length: mesh.nCells }, (_, i) => (!land[i] && !(lake && lake[i]) && inIndexGrid(mesh.latCell[i], mesh.lonCell[i]) ? 1 : 0));
}

export function seaIceIndexCut(mesh, cells, concentration) {
  const out = { N: { extent: 0, area: 0 }, S: { extent: 0, area: 0 } };
  for (let i = 0; i < mesh.nCells; i++) {
    const c = concentration[i];
    if (!cells[i] || !(c >= 0.15)) continue;
    const h = out[mesh.latCell[i] > 0 ? 'N' : 'S'];
    h.extent += mesh.areaCell[i] / 1e12; h.area += c * mesh.areaCell[i] / 1e12;
  }
  return out;
}
