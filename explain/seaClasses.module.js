/*
 * The model's seawater and its interior classes, for Part 8's seawater
 * figure. Temperatures are in kelvin, salinities in psu and densities in
 * kg/m³ at the surface. startSurface is a copy of initializeState's sea
 * (js/physics/init.module.js), which cannot be called without building a
 * whole model; test/seaClasses.test.mjs holds the two equal.
 */
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createOcean, LAYER_DENSITIES, LAYER_SALINITIES, THIN } from '../js/ocean/layered.module.js';
import { SEAWATER, seawaterDensity, thermalExpansion, halineContraction, labelTemperature } from '../js/ocean/seawater.module.js';
import { FREEZING_POINT } from '../js/physics/ice.module.js';
import { FOUR_CLASSES } from './oceanCases.module.js';

export const ALL_CLASSES = { densities: LAYER_DENSITIES, salinities: LAYER_SALINITIES };
export { FOUR_CLASSES, FREEZING_POINT, THIN };
export const CELSIUS = 273.15;

export function classPoints({ densities, salinities } = ALL_CLASSES) {
  return densities.map((rho, k) => ({ rho, s: salinities[k], t: Math.max(FREEZING_POINT, labelTemperature(rho, salinities[k])) }));
}

export function water(t, s) {
  const alpha = thermalExpansion(t, s), beta = halineContraction(t, s);
  return { t, s, rho: seawaterDensity(t, s), alpha, beta, perKelvin: -alpha * SEAWATER.rho0, perPsu: beta * SEAWATER.rho0 };
}

export function nearestClass(rho, points) {
  let best = 0;
  for (let k = 1; k < points.length; k++) if (Math.abs(points[k].rho - rho) < Math.abs(points[best].rho - rho)) best = k;
  return best;
}

export function mixture(a, b, share = 0.5) {
  const t = a.t + (b.t - a.t) * share, s = a.s + (b.s - a.s) * share;
  const rho = seawaterDensity(t, s), ra = seawaterDensity(a.t, a.s), rb = seawaterDensity(b.t, b.s);
  return { t, s, rho, average: ra + (rb - ra) * share, excess: rho - (ra + (rb - ra) * share), overDenser: rho - Math.max(ra, rb) };
}

export function densityMaximum(s) {
  const { a0, lambda1, nu, t0, s0 } = SEAWATER;
  return t0 - (1 + nu * (s - s0) / a0) / lambda1;
}

export function startSurface(mesh) {
  const C = mesh.nCells, surfaceT = new Float64Array(C), ice = new Float64Array(C), deg = 180 / Math.PI;
  for (let i = 0; i < C; i++) {
    const lat = mesh.latCell[i] * deg, polar = Math.max(0, Math.min(1, (Math.abs(lat) - 50) / 20));
    ice[i] = lat > 72 ? 1.5 : lat < -68 ? 0.7 : 0;
    surfaceT[i] = ice[i] > 0 ? FREEZING_POINT : Math.max(FREEZING_POINT + 0.3, CELSIUS + 30 * Math.cos(mesh.latCell[i]) ** 2 - 2 - 4 * polar);
  }
  return { surfaceT, ice };
}

function columns(mesh, ocean, picks) {
  const C = mesh.nCells, L = ocean.layers, n = picks.length;
  const out = {
    layers: L, densities: Float64Array.from(ocean.densities), latitude: new Float64Array(n), eta: new Float64Array(n), iced: new Uint8Array(n),
    thickness: new Float64Array(n * L), temperature: new Float64Array(n * L), salinity: new Float64Array(n * L),
  };
  picks.forEach((i, j) => {
    out.latitude[j] = mesh.latCell[i] * 180 / Math.PI;
    out.eta[j] = ocean.eta[i];
    for (let k = 0; k < L; k++) {
      const h = ocean.h[k * C + i];
      out.thickness[j * L + k] = h;
      out.temperature[j * L + k] = ocean.Q[k * C + i] / h;
      out.salinity[j * L + k] = ocean.W[k * C + i] / h;
    }
  });
  return out;
}

export function classesByLatitude(N = 16, lists = [ALL_CLASSES, FOUR_CLASSES]) {
  const mesh = buildMesh(new Grid(N)), C = mesh.nCells;
  const order = Array.from({ length: C }, (_, i) => i).sort((a, b) => mesh.latCell[a] - mesh.latCell[b]);
  const picks = [];
  for (const i of order) if (!picks.length || mesh.latCell[i] - mesh.latCell[picks.at(-1)] > 1e-9) picks.push(i);
  return lists.map((classes) => {
    const { surfaceT, ice } = startSurface(mesh);
    const ocean = createOcean(mesh, classes);
    ocean.initialize(surfaceT, ice);
    const out = columns(mesh, ocean, picks);
    picks.forEach((i, j) => { out.iced[j] = ice[i] > 0 ? 1 : 0; });
    return out;
  });
}

export function columnAt(section, j) {
  const L = section.layers, base = new Float64Array(L);
  let z = 0;
  for (let k = 0; k < L; k++) { z += section.thickness[j * L + k]; base[k] = z; }
  return base;
}

export function layerAt(section, j, depth) {
  const L = section.layers;
  let top = 0;
  for (let k = 0; k < L; k++) {
    const h = section.thickness[j * L + k];
    if (k > 0 && h <= THIN) continue;
    if (depth < top + h) return { k, top, base: top + h };
    top += h;
  }
  return null;
}

export function thickClasses(section, j) {
  let count = 0;
  for (let k = 1; k < section.layers; k++) if (section.thickness[j * section.layers + k] > THIN) count++;
  return count;
}

export function nearestLatitude(section, latitude) {
  let best = 0;
  for (let j = 1; j < section.latitude.length; j++) if (Math.abs(section.latitude[j] - latitude) < Math.abs(section.latitude[best] - latitude)) best = j;
  return best;
}
