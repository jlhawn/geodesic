/*
 * The model's seawater and its layers, for Part 8's seawater figure.
 * Temperatures are in kelvin, salinities in psu and densities in kg/m³
 * at the surface. The ocean's vertical coordinate is the hybrid:
 * HYBRID_LAYERS gives each of its 58 layers a label density and a
 * minimum thickness: a layer is fixed-depth where its isopycnals lie
 * closer than its minimum and isopycnal where they lie further apart.
 * startSurface is a copy of initializeState's sea (js/physics/init.module.js),
 * which cannot be called without building a whole model;
 * test/seaClasses.test.mjs holds the two equal.
 */
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createOcean } from '../js/ocean/layered.module.js';
import { HYBRID_LAYERS, labelBoundaries, H_MIN, MODE_DEPTH } from '../js/ocean/column.module.js';
import { SEAWATER, seawaterDensity, thermalExpansion, halineContraction, salinityForDensity } from '../js/ocean/seawater.module.js';
import { FREEZING_POINT } from '../js/physics/ice.module.js';

export { HYBRID_LAYERS, FREEZING_POINT, H_MIN, MODE_DEPTH };
export const CELSIUS = 273.15;
export const SURFACE_LAYERS = HYBRID_LAYERS.fixed;
export const SURFACE_DEPTH = HYBRID_LAYERS.depths[SURFACE_LAYERS];
export const DEPTH_SCALE = 20, SECTION_FLOOR = 4000;
const BOUNDARIES = labelBoundaries(HYBRID_LAYERS.densities), STRETCH = Math.log1p(SECTION_FLOOR / DEPTH_SCALE);

export function depthFraction(z) {
  return Math.log1p(Math.min(Math.max(z, 0), SECTION_FLOOR) / DEPTH_SCALE) / STRETCH;
}

export function depthAt(fraction) {
  return DEPTH_SCALE * Math.expm1(Math.min(Math.max(fraction, 0), 1) * STRETCH);
}

export function water(t, s) {
  const alpha = thermalExpansion(t, s), beta = halineContraction(t, s);
  return { t, s, rho: seawaterDensity(t, s), alpha, beta, perKelvin: -alpha * SEAWATER.rho0, perPsu: beta * SEAWATER.rho0 };
}

export function labelOf(rho) {
  let k = 0;
  while (k < HYBRID_LAYERS.K - 1 && rho >= BOUNDARIES[k + 1]) k++;
  return k;
}

export function labelRange(k) {
  return { below: BOUNDARIES[k], above: BOUNDARIES[k + 1] };
}

export function isopycnalLine(rho, tLow, tHigh, step = 0.25) {
  const out = [];
  for (let t = tLow; t <= tHigh + 1e-9; t += step) out.push([salinityForDensity(rho, t + CELSIUS), t]);
  return out;
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

export function layerWeight(h, k) {
  return Math.min(1, Math.max(0, (h - HYBRID_LAYERS.minimum[k]) / MODE_DEPTH));
}

export function layerMode(h, k) {
  if (h <= H_MIN + 1e-6) return 'empty';
  const w = layerWeight(h, k);
  return w < 1e-6 ? 'fixed' : w > 1 - 1e-6 ? 'isopycnal' : 'between';
}

export function layersByLatitude(N = 16) {
  const mesh = buildMesh(new Grid(N)), C = mesh.nCells, { surfaceT, ice } = startSurface(mesh);
  const order = Array.from({ length: C }, (_, i) => i).sort((a, b) => mesh.latCell[a] - mesh.latCell[b]);
  const picks = [];
  for (const i of order) if (!picks.length || mesh.latCell[i] - mesh.latCell[picks.at(-1)] > 1e-9) picks.push(i);
  const ocean = createOcean(mesh);
  ocean.initialize(surfaceT, ice);
  const K = ocean.layers, n = picks.length;
  const out = {
    layers: K, latitude: new Float64Array(n), iced: new Uint8Array(n), floor: new Float64Array(n),
    interface: new Float64Array(n * (K + 1)), temperature: new Float64Array(n * K), salinity: new Float64Array(n * K), density: new Float64Array(n * K), weight: new Float64Array(n * K),
  };
  picks.forEach((i, j) => {
    out.latitude[j] = mesh.latCell[i] * 180 / Math.PI;
    out.iced[j] = ice[i] > 0 ? 1 : 0;
    let z = 0;
    for (let k = 0; k < K; k++) {
      const h = ocean.h[k * C + i], m = j * K + k;
      out.interface[j * (K + 1) + k] = z;
      z += h;
      out.temperature[m] = ocean.Q[k * C + i] / h;
      out.salinity[m] = ocean.W[k * C + i] / h;
      out.density[m] = seawaterDensity(out.temperature[m], out.salinity[m]);
      out.weight[m] = layerWeight(h, k);
    }
    out.interface[j * (K + 1) + K] = z;
    out.floor[j] = z;
  });
  return out;
}

export function thickness(section, j, k) {
  const base = j * (section.layers + 1);
  return section.interface[base + k + 1] - section.interface[base + k];
}

export function modeAt(section, j, k) {
  return layerMode(thickness(section, j, k), k);
}

export function layerAt(section, j, depth) {
  const K = section.layers, base = j * (K + 1);
  if (depth >= section.floor[j]) return null;
  for (let k = 0; k < K; k++) {
    const top = section.interface[base + k], bottom = section.interface[base + k + 1];
    if (depth < bottom) return { k, top, base: bottom, mode: layerMode(bottom - top, k), weight: section.weight[j * K + k] };
  }
  return null;
}

export function modeCounts(section, j) {
  const counts = { fixed: 0, between: 0, isopycnal: 0, empty: 0 };
  for (let k = 0; k < section.layers; k++) counts[modeAt(section, j, k)]++;
  return counts;
}

export function emptyEverywhere(section) {
  let count = 0;
  for (let k = section.layers - 1; k >= 0; k--) {
    let empty = true;
    for (let j = 0; j < section.latitude.length && empty; j++) empty = modeAt(section, j, k) === 'empty';
    if (!empty) break;
    count++;
  }
  return count;
}

export function nearestLatitude(section, latitude) {
  let best = 0;
  for (let j = 1; j < section.latitude.length; j++) if (Math.abs(section.latitude[j] - latitude) < Math.abs(section.latitude[best] - latitude)) best = j;
  return best;
}
