/*
 * The sea-ice extent and the iced polar cells that spinup.mjs logs each day, over the cells the model runs as sea:
 * neither land nor a lake the lake model holds. With the lake model off, a lake cell is sea and its ice is the sea's.
 */
const DEG = 180 / Math.PI;

export function seaCells(model) {
  const land = model.geography.land, lake = model.lake ? model.lake.isLake : null;
  return Uint8Array.from(land, (l, i) => (l || (lake && lake[i]) ? 0 : 1));
}

export function seaIceExtent(mesh, sea, concentration) {
  let north = 0, south = 0;
  for (let i = 0; i < mesh.nCells; i++) {
    if (!sea[i] || !(concentration[i] >= 0.15)) continue;
    if (mesh.latCell[i] > 0) north += mesh.areaCell[i]; else south += mesh.areaCell[i];
  }
  return [north / 1e12, south / 1e12];
}

export function polarIcedCells(mesh, sea, ice) {
  const cells = [];
  for (let i = 0; i < mesh.nCells; i++) if (sea[i] && ice[i] > 0 && Math.abs(mesh.latCell[i] * DEG) >= 60) cells.push(i);
  return cells;
}
