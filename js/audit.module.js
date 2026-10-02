import { cellVector } from './dynamics/operators.module.js';

/*
 * The boxes of the vertical-motion and convection audit (M21 in
 * docs/c-grid-dynamical-core.md) and the two lines the spin-up log
 * prints from them. A box is [south, north, west, east] in degrees,
 * east of west across the date line when west > east (160, -100 is
 * 160E-100W).
 */
export const BOXES = {
  sePacific: [-30, -10, -110, -80], peru: [-20, -5, -90, -75], itcz: [5, 12, 160, -100], namibia: [-20, -10, 0, 10], california: [20, 30, -130, -120],
};
const DEG = 180 / Math.PI;

export function inLongitudes(lon, west, east) {
  if (east - west >= 360) return true;
  const span = (((east - west) % 360) + 360) % 360, offset = (((lon - west) % 360) + 360) % 360;
  return offset <= span;
}

/*
 * The tropical boxes of scripts/tropicalHeating.mjs and of the audit's
 * heating rows, [name, box, surface] with surface 'all', 'sea' or 'land';
 * tropicalBoxOf gives each cell the index of the first box that holds it
 * (-1 for none), so that a cell in two boxes counts in the first.
 */
export const TROPICAL_BOXES = [
  ['Pacific ITCZ 5-12N 160E-100W', BOXES.itcz, 'all'],
  ['warm pool 10S-10N 120-170E sea', [-10, 10, 120, 170], 'sea'],
  ['SPCZ 20-5S 160E-150W sea', [-20, -5, 160, -150], 'sea'],
  ['N Pacific trades 15-25N 170-130W sea', [15, 25, -170, -130], 'sea'],
  ['Amazon 10S-2N 70-50W land', [-10, 2, -70, -50], 'land'],
];

export function tropicalBoxOf(mesh, land, boxes = TROPICAL_BOXES) {
  const boxOf = new Int8Array(mesh.nCells).fill(-1);
  for (let i = 0; i < mesh.nCells; i++) {
    const lat = mesh.latCell[i] * DEG, lon = mesh.lonCell[i] * DEG;
    boxes.forEach(([, [south, north, west, east], surface], b) => {
      const kept = surface === 'all' || (surface === 'land') === !!land[i];
      if (boxOf[i] < 0 && lat >= south && lat <= north && inLongitudes(lon, west, east) && kept) boxOf[i] = b;
    });
  }
  return boxOf;
}

/*
 * A box's apparent heat source or moisture sink Q (K/day per layer, top
 * first) at the layers' mean pressures p and thicknesses dp (Pa), over the
 * layers below `top`: the layer of largest Q, the `width` (50 hPa) bin of
 * largest mass-weighted mean Σ Q dp / Σ dp over the layers whose
 * midpoints fall in it, and the centroid Σ p Q dp / Σ Q dp over the
 * layers where Q > 0.
 */
export function heatingProfile(Q, p, dp, { top = 100e2, width = 50e2 } = {}) {
  let best = -1, weighted = 0, weight = 0;
  const bins = new Map();
  for (let k = 0; k < Q.length; k++) {
    if (!(p[k] >= top)) continue;
    if (best < 0 || Q[k] > Q[best]) best = k;
    const key = Math.floor(p[k] / width), bin = bins.get(key) ?? [0, 0];
    bin[0] += Q[k] * dp[k]; bin[1] += dp[k];
    bins.set(key, bin);
    if (Q[k] > 0) { weighted += p[k] * Q[k] * dp[k]; weight += Q[k] * dp[k]; }
  }
  let bin = null, binValue = -Infinity;
  for (const [key, [sum, mass]] of bins) if (sum / mass > binValue) { binValue = sum / mass; bin = key; }
  return { layer: best >= 0 ? p[best] : NaN, value: best >= 0 ? Q[best] : NaN, bin: bin === null ? [NaN, NaN] : [bin * width, (bin + 1) * width], binValue, centroid: weight > 0 ? weighted / weight : NaN };
}

/*
 * The bulk sensible heat (W/m²) the CPU model's radiation step gives cell
 * i from the lowest air temperature airT and the surface temperature
 * surfaceT it read, and the ice thickness before the step: over partly
 * iced sea the skin blends in open water at the freezing point.
 */
export function bulkSensible(model, i, airT, surfaceT, ice, { seaDrag, landDrag, freezing, gustiness = 3 }) {
  const pi = model.state[0], { sigmaMid, K, R, cp } = model.core.diagnostics;
  let skin = surfaceT;
  const land = model.geography.land[i];
  if (!land) { const cover = model.seaIce.cover(i, ice); if (ice > 0 && cover < 1) skin = cover * surfaceT + (1 - cover) * freezing; }
  const density = pi[i] * sigmaMid[K - 1] / (R * airT);
  return density * (land ? landDrag : seaDrag) * Math.max(model.surface.windSpeed[i], gustiness) * cp * (skin - airT);
}

/*
 * The layer Exner functions of a column of surface pressure ps, as the
 * core forms them: the mass-weighted mean of (p/p0)^κ over each layer.
 */
export function layerExner(ps, { K, sigmaLower, sigmaUpper, dSigma, kappa, p0 }, out = new Float64Array(K)) {
  let upper = 0;
  for (let k = 0; k < K; k++) {
    const lower = Math.pow(ps * sigmaLower[k] / p0, kappa);
    out[k] = (lower * sigmaLower[k] - upper * sigmaUpper[k]) / ((1 + kappa) * dSigma[k]);
    upper = lower;
  }
  return out;
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
