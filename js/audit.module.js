import { cellVector } from './dynamics/operators.module.js';

/*
 * The boxes of the vertical-motion and convection audit (M21 in
 * docs/c-grid-dynamical-core.md) and the two lines the spin-up log
 * prints from them. A box is [south, north, west, east] in degrees,
 * east of west across the date line when west > east (160, -100 is
 * 160E-100W).
 */
export const BOXES = {
  sePacific: [-30, -10, -110, -80], peru: [-20, -5, -90, -75], itcz: [5, 12, 160, -100],
};
const DEG = 180 / Math.PI;

export function inLongitudes(lon, west, east) {
  if (east - west >= 360) return true;
  const span = (((east - west) % 360) + 360) % 360, offset = (((lon - west) % 360) + 360) % 360;
  return offset <= span;
}

export function boxCells(mesh, [south, north, west, east], keep = () => true) {
  const cells = [];
  for (let i = 0; i < mesh.nCells; i++) {
    const lat = mesh.latCell[i] * DEG, lon = mesh.lonCell[i] * DEG;
    if (lat >= south && lat <= north && inLongitudes(lon, west, east) && keep(i)) cells.push(i);
  }
  return cells;
}

export function areaMean(mesh, cells, value) {
  let sum = 0, area = 0;
  for (const i of cells) { const a = mesh.areaCell[i]; sum += a * value(i); area += a; }
  return area > 0 ? sum / area : NaN;
}

/*
 * One line on the rain split over `days`: `convective` and `largeScale`
 * are each cell's mean rain in mm/d, `wet` the fraction of its days with
 * any convective rain; the SE Pacific box keeps its ice-free sea cells.
 */
export function convectionLine(mesh, { land, ice }, { convective, largeScale, wet }, days) {
  const all = boxCells(mesh, [-90, 90, -180, 180]), tropics = boxCells(mesh, [-15, 15, -180, 180]);
  const share = (cells) => areaMean(mesh, cells, (i) => convective[i]) / areaMean(mesh, cells, (i) => convective[i] + largeScale[i]);
  const deck = boxCells(mesh, BOXES.sePacific, (i) => !land[i] && !(ice[i] > 0)), itcz = boxCells(mesh, BOXES.itcz);
  const total = (cells) => areaMean(mesh, cells, (i) => convective[i] + largeScale[i]);
  const f = (x, d = 2) => x.toFixed(d);
  return `convection after ${days} days: convective share global ${f(share(all))}, 15S-15N ${f(share(tropics))}; SE Pacific 10-30S 110-80W sea ${f(total(deck))} mm/d `
    + `(convective ${f(areaMean(mesh, deck, (i) => convective[i]))}), convective rain on ${f(areaMean(mesh, deck, (i) => wet[i]))} of its column-days; `
    + `Pacific ITCZ 5-12N 160E-100W ${f(total(itcz))} mm/d`;
}

/*
 * The equatorial Pacific ocean, 2S-2N, from a serialized ocean (h and u
 * layer-major with the mixed layer first, and the interior class
 * densities) and the wind stress on its edges (N/m²): the mixed layer's
 * eastward current, the strongest thickness-weighted mean eastward
 * current of the classes to 1026.0 over 180-100W with the class's mean
 * depth, the eastward stress, the mixed-layer depth in the east and the
 * depth of the top of the class `thermocline` west and east.
 */
export function equatorialOcean(mesh, land, ocean, stress, thermocline = 1024) {
  const C = mesh.nCells, E = mesh.nEdges, L = ocean.h.length / C, densities = ocean.densities;
  const band = (west, east) => boxCells(mesh, [-2, 2, west, east], (i) => !land[i]);
  const eastward = (i, vector) => -Math.sin(mesh.lonCell[i]) * vector[3 * i] + Math.cos(mesh.lonCell[i]) * vector[3 * i + 1];
  const layerVector = (k) => cellVector(mesh, Float64Array.from(ocean.u.slice(k * E, (k + 1) * E)));
  const h = (k, i) => ocean.h[k * C + i];
  const surface = layerVector(0), wind = cellVector(mesh, Float64Array.from(stress));
  const pacific = band(160, -100), eastern = band(-140, -100), central = band(-180, -100);
  let under = { u: -Infinity, depth: NaN, density: NaN };
  for (let k = 1; k < L && densities[k - 1] <= 1026.0; k++) {
    const vector = layerVector(k);
    let mass = 0, flow = 0, depth = 0, area = 0;
    for (const i of central) {
      let above = 0;
      for (let j = 0; j < k; j++) above += h(j, i);
      const a = mesh.areaCell[i], thick = h(k, i);
      mass += a * thick; flow += a * thick * eastward(i, vector); depth += a * thick * (above + thick / 2); area += a;
    }
    if (mass / area < 5) continue;
    if (flow / mass > under.u) under = { u: flow / mass, depth: depth / mass, density: densities[k - 1] };
  }
  const lighter = densities.filter((r) => r < thermocline).length;
  const classTop = (cells) => areaMean(mesh, cells, (i) => { let d = 0; for (let k = 0; k <= lighter; k++) d += h(k, i); return d; });
  return {
    surface: areaMean(mesh, pacific, (i) => eastward(i, surface)), surfaceEast: areaMean(mesh, eastern, (i) => eastward(i, surface)),
    undercurrent: under.u, undercurrentDepth: under.depth, undercurrentClass: under.density,
    stress: areaMean(mesh, pacific, (i) => eastward(i, wind)), mixedEast: areaMean(mesh, eastern, (i) => h(0, i)),
    thermoclineWest: classTop(band(150, 180)), thermoclineEast: classTop(band(-120, -90)), thermocline,
  };
}

export function equatorLine(mesh, land, ocean, stress, days) {
  const o = equatorialOcean(mesh, land, ocean, stress);
  const v = (x) => `${x >= 0 ? '+' : ''}${x.toFixed(2)}`;
  return `equator after ${days} days (2S-2N, eastward +): surface current 160E-100W ${v(o.surface)} m/s, 140W-100W ${v(o.surfaceEast)} m/s; `
    + `undercurrent ${v(o.undercurrent)} m/s at ${o.undercurrentDepth.toFixed(0)} m (class ${o.undercurrentClass}, 180-100W); stress 160E-100W ${o.stress.toFixed(3)} N/m²; `
    + `mixed layer 140W-100W ${o.mixedEast.toFixed(0)} m; ${o.thermocline} class top 150E-180 ${o.thermoclineWest.toFixed(0)} m, 120W-90W ${o.thermoclineEast.toFixed(0)} m`;
}
