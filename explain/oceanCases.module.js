/*
 * Small cases of the model's own ocean for Part 8's figures: createOcean,
 * unchanged, on buildMesh(new Grid(N)), or on a copy of that mesh whose
 * fEdge and fVertex hold one value ("constant f"), started from the
 * analytic profile with a salinity of 35 and no tilt of the thermocline, so
 * that every column is the same, and with no heat, fresh water or ice. The
 * layers are the model's 58 (ocean.table). Cell i of every array is cell i
 * of new Grid(N), the order the explainer's globes draw. Units are SI (m, s,
 * m/s, m²/s); vectors come as east and north components; fields are
 * Float32Array and every array a frame returns is a fresh copy. A case
 * factory returns a simulation that steps the ocean and reports what a
 * figure draws; ekman, lens, bump and basin run one to its end and return
 * its last frame.
 *
 * The Ekman layer of a column is read from its layers' ring-mean currents:
 * the reference is the thickness-weighted current of the layers between
 * EKMAN_REACH and 1.5 EKMAN_REACH boundary-layer depths down (at least
 * EKMAN_FLOOR meters), and the Ekman transport is the water above
 * EKMAN_REACH boundary-layer depths times its current relative to that.
 */
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { cellVector, gradient, tangential } from '../js/dynamics/operators.module.js';
import { createOcean, depthFields } from '../js/ocean/layered.module.js';
import { seawaterDensity } from '../js/ocean/seawater.module.js';
import { createGeography, syntheticTopography } from '../js/geography.module.js';

export const RHO = 1025, G = 9.81, DEG = Math.PI / 180, DAY = 86400;
export const EKMAN_REACH = 2, EKMAN_FLOOR = 40, PUMP_DEPTH = 100;

const meshes = new Map();
export function caseMesh(N) {
  if (!meshes.has(N)) meshes.set(N, buildMesh(new Grid(N)));
  return meshes.get(N);
}

export function zonalWindOnEdges(m, speedAt) {
  const u = new Float64Array(m.nEdges);
  for (let e = 0; e < m.nEdges; e++) {
    const x = m.xEdge[3 * e], y = m.xEdge[3 * e + 1], r = Math.hypot(x, y);
    const east = r > 0 ? [-y / r, x / r, 0] : [0, 0, 0];
    u[e] = speedAt(m.latEdge[e]) * (east[0] * m.nEdge[3 * e] + east[1] * m.nEdge[3 * e + 1] + east[2] * m.nEdge[3 * e + 2]);
  }
  return u;
}

const bases = new WeakMap();
export function cellBasis(mesh) {
  if (bases.has(mesh)) return bases.get(mesh);
  const C = mesh.nCells, east = new Float64Array(3 * C), north = new Float64Array(3 * C);
  for (let i = 0; i < C; i++) {
    const x = mesh.xCell[3 * i], y = mesh.xCell[3 * i + 1], z = mesh.xCell[3 * i + 2], r = Math.hypot(x, y);
    if (r > 1e-12) { east.set([-y / r, x / r, 0], 3 * i); north.set([-z * x / r, -z * y / r, r], 3 * i); }
    else { east.set([0, 1, 0], 3 * i); north.set([-Math.sign(z), 0, 0], 3 * i); }
  }
  const basis = { east, north };
  bases.set(mesh, basis);
  return basis;
}
const dot = (v, basis, i) => v[3 * i] * basis[3 * i] + v[3 * i + 1] * basis[3 * i + 1] + v[3 * i + 2] * basis[3 * i + 2];

export function distancesFrom(mesh, latitude, longitude, points = 'cell') {
  const p = [Math.cos(latitude * DEG) * Math.cos(longitude * DEG), Math.cos(latitude * DEG) * Math.sin(longitude * DEG), Math.sin(latitude * DEG)];
  const x = points === 'edge' ? mesh.xEdge : mesh.xCell, n = x.length / 3, out = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const c = (x[3 * i] * p[0] + x[3 * i + 1] * p[1] + x[3 * i + 2] * p[2]) / Math.hypot(x[3 * i], x[3 * i + 1], x[3 * i + 2]);
    out[i] = Math.acos(Math.max(-1, Math.min(1, c))) * mesh.radius;
  }
  return out;
}

export function meanSpacing(mesh) {
  let sum = 0;
  for (let e = 0; e < mesh.nEdges; e++) sum += mesh.dcEdge[e];
  return sum / mesh.nEdges;
}

const START = { everySteps: 1, thermoclineTilt: 0, salinityProfile: () => 35 };
const recipe = (extra = {}) => ({ ...START, ...extra });

function layoutOf(kind, N, mesh, ocean, extra = {}) {
  return {
    kind, N, nCells: mesh.nCells, layers: ocean.layers, fixed: ocean.table.fixed, densities: Float64Array.from(ocean.table.densities), minimum: Float64Array.from(ocean.table.minimum),
    latitude: Float32Array.from(mesh.latCell, (v) => v / DEG), longitude: Float32Array.from(mesh.lonCell, (v) => v / DEG),
    sea: Uint8Array.from(ocean.cellOcean), ...extra,
  };
}

function components(mesh, edgeValues, east, north, vector = new Float64Array(3 * mesh.nCells)) {
  const basis = cellBasis(mesh);
  cellVector(mesh, edgeValues, vector);
  for (let i = 0; i < mesh.nCells; i++) { east[i] = dot(vector, basis.east, i); north[i] = dot(vector, basis.north, i); }
}

function layerSpeeds(mesh, ocean, cells) {
  const { nCells: C, nEdges: E } = mesh, L = ocean.layers, vector = new Float64Array(3 * C), out = new Float32Array(L);
  for (let k = 0; k < L; k++) {
    cellVector(mesh, ocean.u.subarray(k * E, (k + 1) * E), vector);
    let sum = 0, area = 0;
    for (const i of cells) {
      if (!(ocean.h[k * C + i] > 0)) continue;
      sum += mesh.areaCell[i] * Math.hypot(vector[3 * i], vector[3 * i + 1], vector[3 * i + 2]); area += mesh.areaCell[i];
    }
    out[k] = area > 0 ? sum / area : NaN;
  }
  return out;
}

function fastestCurrent(mesh, ocean) {
  const { nCells: C, nEdges: E } = mesh, vector = new Float64Array(3 * C);
  let fastest = 0;
  for (let k = 0; k < ocean.layers; k++) {
    cellVector(mesh, ocean.u.subarray(k * E, (k + 1) * E), vector);
    for (let i = 0; i < C; i++) if (ocean.cellOcean[i] && ocean.h[k * C + i] > 0) fastest = Math.max(fastest, Math.hypot(vector[3 * i], vector[3 * i + 1], vector[3 * i + 2]));
  }
  return fastest;
}

export function ekmanLayer(thick, east, north, base, reach = EKMAN_REACH) {
  const depth = Math.max(reach * base, EKMAN_FLOOR), L = thick.length;
  let z = 0, weight = 0, refEast = 0, refNorth = 0;
  for (let k = 0; k < L; k++) {
    const middle = z + 0.5 * thick[k];
    z += thick[k];
    if (middle > depth && middle < 1.5 * depth) { weight += thick[k]; refEast += thick[k] * east[k]; refNorth += thick[k] * north[k]; }
  }
  if (weight > 0) { refEast /= weight; refNorth /= weight; }
  let transportEast = 0, transportNorth = 0;
  z = 0;
  for (let k = 0; k < L; k++) {
    const middle = z + 0.5 * thick[k];
    z += thick[k];
    if (middle < depth) { transportEast += thick[k] * (east[k] - refEast); transportNorth += thick[k] * (north[k] - refNorth); }
  }
  return { depth, east: refEast, north: refNorth, transportEast, transportNorth };
}

export function meanWithin(thick, values, depth) {
  let z = 0, sum = 0, within = 0;
  for (let k = 0; k < thick.length && z < depth; k++) {
    const part = Math.min(thick[k], depth - z);
    sum += part * values[k]; within += part; z += thick[k];
  }
  return within > 0 ? sum / within : 0;
}

function drift(east, north, theory = 1) {
  const magnitude = Math.hypot(east, north);
  return { east, north, magnitude, angle: Math.atan2(-north, east) / DEG, ratio: magnitude / theory };
}

function fit(points) {
  const n = points.length;
  if (n < 2) return null;
  let mt = 0, mr = 0;
  for (const [t, r] of points) { mt += t / n; mr += r / n; }
  let num = 0, den = 0;
  for (const [t, r] of points) { num += (t - mt) * (r - mr); den += (t - mt) ** 2; }
  return den > 0 ? num / den : null;
}

export const PLANET_WINDS = {
  westerlies: (lat) => Math.exp(-(((Math.abs(lat) / DEG - 45) / 10) ** 2)),
  trades: (lat) => {
    const l = Math.abs(lat) / DEG;
    if (l < 30) return -0.5 * Math.cos(Math.PI * l / 60);
    if (l < 65) return Math.sin(Math.PI * (l - 30) / 35);
    return -0.2 * Math.sin(Math.PI * (l - 65) / 25);
  },
};
export const columnWind = (latitude) => (lat) => Math.exp(-(((lat / DEG - latitude) / 10) ** 2));

export function ekmanCase({ N = 8, tau = 0.1, latitude = 45, days = 4, dt = 1350, view = 'column', wind = 'trades', averaging = 2 * DAY, sst = 290 } = {}) {
  const mesh = caseMesh(N), { nCells: C, nEdges: E, areaCell } = mesh;
  const ocean = createOcean(mesh, recipe()), L = ocean.layers;
  const surfaceT = new Float64Array(C).fill(sst), ice = new Float64Array(C), flux = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  if (view !== 'column' && view !== 'planet') throw new Error(`the Ekman view is 'column' or 'planet', not ${view}`);
  if (view === 'planet' && !PLANET_WINDS[wind]) throw new Error(`the planet's wind is ${Object.keys(PLANET_WINDS).join(' or ')}, not ${wind}`);
  const shape = view === 'planet' ? PLANET_WINDS[wind] : columnWind(latitude);
  const stress = zonalWindOnEdges(mesh, (lat) => tau * shape(lat));
  const basis = cellBasis(mesh), half = 0.5 * meanSpacing(mesh) / mesh.radius / DEG;
  const coriolis = (lat) => 2 * mesh.omega * Math.sin(lat);
  const layerVectors = Array.from({ length: L }, () => new Float64Array(3 * C));
  const allLayers = () => { for (let k = 0; k < L; k++) cellVector(mesh, ocean.u.subarray(k * E, (k + 1) * E), layerVectors[k]); };
  const sim = { kind: 'ekman', view, mesh, ocean, time: 0, end: days * DAY, nextDt: () => dt };
  let advanceExtra = () => {}, frameExtra = () => ({}), layoutExtra = () => ({});

  if (view === 'column') {
    const ring = [];
    let ringArea = 0, probe = 0, best = Infinity, theory = 0, meanStress = 0, meanF = 0;
    for (let i = 0; i < C; i++) {
      const lat = mesh.latCell[i] / DEG, lon = mesh.lonCell[i] / DEG;
      if (Math.abs(lat - latitude) < half) { ring.push(i); ringArea += areaCell[i]; }
      const miss = Math.abs(lat - latitude) + 0.01 * Math.abs(lon);
      if (miss < best) { best = miss; probe = i; }
    }
    for (const i of ring) {
      const w = areaCell[i] / ringArea;
      theory += w * tau * shape(mesh.latCell[i]) / (RHO * Math.abs(coriolis(mesh.latCell[i])));
      meanStress += w * tau * shape(mesh.latCell[i]); meanF += w * coriolis(mesh.latCell[i]);
    }
    const f = coriolis(latitude * DEG), period = 2 * Math.PI / Math.abs(f);
    const snaps = [], series = {};
    for (const name of ['time', 'topEast', 'topNorth', 'mixedEast', 'mixedNorth', 'deepEast', 'deepNorth', 'mixedDepth', 'boundaryDepth']) series[name] = [];
    const record = () => {
      allLayers();
      const snap = { time: sim.time, east: new Float64Array(L), north: new Float64Array(L), thick: new Float64Array(L), nu: new Float64Array(L), mixed: 0, boundary: 0 };
      for (const i of ring) {
        const w = areaCell[i] / ringArea;
        for (let k = 0; k < L; k++) {
          snap.east[k] += w * dot(layerVectors[k], basis.east, i); snap.north[k] += w * dot(layerVectors[k], basis.north, i);
          snap.thick[k] += w * ocean.h[k * C + i]; snap.nu[k] += w * ocean.viscosity[k * C + i];
        }
        snap.mixed += w * ocean.mixedDepths[i]; snap.boundary += w * ocean.boundaryDepths[i];
      }
      const layer = ekmanLayer(snap.thick, snap.east, snap.north, snap.boundary);
      snaps.push(snap);
      const row = [sim.time, snap.east[0], snap.north[0], meanWithin(snap.thick, snap.east, snap.mixed), meanWithin(snap.thick, snap.north, snap.mixed), layer.east, layer.north, snap.mixed, snap.boundary];
      Object.keys(series).forEach((name, j) => series[name].push(row[j]));
    };
    record();
    advanceExtra = record;
    frameExtra = () => {
      const settled = sim.time >= 2 * period, from = settled ? sim.time - 2 * period + 0.5 * dt : sim.time;
      const used = snaps.filter((s) => s.time >= from), n = used.length;
      const mean = { east: new Float64Array(L), north: new Float64Array(L), thick: new Float64Array(L), nu: new Float64Array(L), mixed: 0, boundary: 0 };
      for (const s of used) {
        for (let k = 0; k < L; k++) { mean.east[k] += s.east[k] / n; mean.north[k] += s.north[k] / n; mean.thick[k] += s.thick[k] / n; mean.nu[k] += s.nu[k] / n; }
        mean.mixed += s.mixed / n; mean.boundary += s.boundary / n;
      }
      const layer = ekmanLayer(mean.thick, mean.east, mean.north, mean.boundary);
      const center = new Float32Array(L), thickness = new Float32Array(L), east = new Float32Array(L), north = new Float32Array(L), viscosity = new Float32Array(L);
      let z = 0, weight = 0, nuSum = 0;
      for (let k = 0; k < L; k++) {
        center[k] = z + 0.5 * mean.thick[k]; thickness[k] = mean.thick[k];
        east[k] = mean.east[k] - layer.east; north[k] = mean.north[k] - layer.north; viscosity[k] = mean.nu[k];
        if (center[k] < mean.boundary) { weight += mean.thick[k]; nuSum += mean.thick[k] * mean.nu[k]; }
        z += mean.thick[k];
      }
      const mixedEast = meanWithin(mean.thick, east, mean.mixed), mixedNorth = meanWithin(mean.thick, north, mean.mixed);
      const copy = (name, type = Float32Array) => type.from(series[name]);
      return {
        latitude, cellLatitude: mesh.latCell[probe] / DEG, f, inertialPeriod: period, theory, ringCells: ring.length, stressMean: meanStress, fMean: meanF,
        mixedDepth: mean.mixed, boundaryDepth: mean.boundary, averagedOver: used.length > 1 ? sim.time - used[0].time + dt : 0, settled,
        reach: layer.depth, reference: { east: layer.east, north: layer.north },
        transport: drift(layer.transportEast, layer.transportNorth, theory), top: drift(east[0], north[0]), mixed: drift(mixedEast, mixedNorth),
        meanViscosity: weight > 0 ? nuSum / weight : 0,
        profile: { center, thickness, east, north, viscosity },
        series: { time: copy('time', Float64Array), topEast: copy('topEast'), topNorth: copy('topNorth'), mixedEast: copy('mixedEast'), mixedNorth: copy('mixedNorth'), deepEast: copy('deepEast'), deepNorth: copy('deepNorth'), mixedDepth: copy('mixedDepth'), boundaryDepth: copy('boundaryDepth') },
      };
    };
    layoutExtra = () => ({ ring: Int32Array.from(ring), probe });
  } else {
    const capacity = Math.max(1, Math.ceil(averaging / dt)), width = 5 * C, snapshots = new Float32Array(capacity * width), sums = new Float64Array(width);
    const thick = new Float64Array(L), east = new Float64Array(L), north = new Float64Array(L);
    let count = 0, next = 0;
    advanceExtra = () => {
      allLayers();
      const w = depthFields(mesh, L, { h: ocean.h, u: ocean.u, temperature: () => 0, cellOcean: ocean.cellOcean }, PUMP_DEPTH, { thin: 0 }).upwelling;
      const slot = snapshots.subarray(next * width, (next + 1) * width);
      if (count === capacity) for (let j = 0; j < width; j++) sums[j] -= slot[j];
      for (let i = 0; i < C; i++) {
        for (let k = 0; k < L; k++) { thick[k] = ocean.h[k * C + i]; east[k] = dot(layerVectors[k], basis.east, i); north[k] = dot(layerVectors[k], basis.north, i); }
        const layer = ekmanLayer(thick, east, north, ocean.boundaryDepths[i]);
        slot[i] = layer.transportEast; slot[C + i] = layer.transportNorth; slot[2 * C + i] = w[i]; slot[3 * C + i] = ocean.boundaryDepths[i]; slot[4 * C + i] = ocean.mixedDepths[i];
      }
      for (let j = 0; j < width; j++) sums[j] += slot[j];
      next = (next + 1) % capacity; count = Math.min(capacity, count + 1);
    };
    const bins = 36;
    frameExtra = () => {
      const mean = (part) => Float32Array.from({ length: C }, (_, i) => (count ? sums[part * C + i] / count : 0));
      const ekmanEast = mean(0), ekmanNorth = mean(1), upwelling = mean(2), boundaryDepth = mean(3), mixedDepth = mean(4);
      const zonal = { latitude: new Float32Array(bins), stress: new Float32Array(bins), theory: new Float32Array(bins), ekmanNorth: new Float32Array(bins), upwelling: new Float32Array(bins) };
      const area = new Float64Array(bins);
      let equator = 0, equatorArea = 0;
      for (let i = 0; i < C; i++) {
        const lat = mesh.latCell[i] / DEG, b = Math.min(bins - 1, Math.floor((lat + 90) / 5)), a = areaCell[i];
        zonal.ekmanNorth[b] += a * ekmanNorth[i]; zonal.upwelling[b] += a * upwelling[i]; area[b] += a;
        if (Math.abs(lat) < half) { equator += a * upwelling[i]; equatorArea += a; }
      }
      for (let b = 0; b < bins; b++) {
        const lat = (-87.5 + 5 * b) * DEG;
        zonal.latitude[b] = lat / DEG; zonal.stress[b] = tau * shape(lat); zonal.theory[b] = -tau * shape(lat) / (RHO * coriolis(lat));
        zonal.ekmanNorth[b] /= area[b]; zonal.upwelling[b] /= area[b];
      }
      return {
        wind, ekmanEast, ekmanNorth, upwelling, boundaryDepth, mixedDepth, zonal, equatorUpwelling: equator / equatorArea, pumpDepth: PUMP_DEPTH,
        eta: Float32Array.from(ocean.eta), stress: Float32Array.from(mesh.latCell, (lat) => tau * shape(lat)), averagedOver: count * dt,
      };
    };
  }

  sim.step = () => {
    ocean.advance(surfaceT, ice, flux, stress, dt);
    sim.time += dt;
    advanceExtra();
  };
  sim.layout = () => layoutOf('ekman', N, mesh, ocean, { view, ...layoutExtra() });
  sim.frame = () => ({ kind: 'ekman', view, tau, time: sim.time, day: sim.time / DAY, done: sim.time >= sim.end, ...frameExtra() });
  return sim;
}

export function internalWaveSpeed(thickness, density, { gravity = G, reference = RHO, merge = 1e-4 } = {}) {
  const h = [], rho = [];
  for (let k = 0; k < thickness.length; k++) {
    if (!(thickness[k] > 0)) continue;
    if (h.length && Math.abs(density[k] - rho[rho.length - 1]) < merge) { const n = h.length - 1; rho[n] = (rho[n] * h[n] + density[k] * thickness[k]) / (h[n] + thickness[k]); h[n] += thickness[k]; }
    else { h.push(thickness[k]); rho.push(density[k]); }
  }
  const n = h.length - 1;
  if (n < 1) return NaN;
  const g = Array.from({ length: n }, (_, m) => Math.max(1e-9, gravity * (rho[m + 1] - rho[m]) / reference));
  const diag = Array.from({ length: n }, (_, m) => 1 / h[m] + 1 / h[m + 1]), off = Array.from({ length: n - 1 }, (_, m) => -1 / h[m + 1]);
  let w = Float64Array.from({ length: n }, () => 1), lambda = 0;
  const solve = (b) => {
    const c = new Float64Array(n), d = new Float64Array(n);
    c[0] = off[0] / diag[0]; d[0] = b[0] / diag[0];
    for (let m = 1; m < n; m++) { const den = diag[m] - off[m - 1] * c[m - 1]; c[m] = m < n - 1 ? off[m] / den : 0; d[m] = (b[m] - off[m - 1] * d[m - 1]) / den; }
    const x = new Float64Array(n);
    x[n - 1] = d[n - 1];
    for (let m = n - 2; m >= 0; m--) x[m] = d[m] - c[m] * x[m + 1];
    return x;
  };
  for (let iteration = 0; iteration < 400; iteration++) {
    const x = solve(Float64Array.from(w, (v, m) => g[m] * v));
    let top = 0, bottom = 0;
    for (let m = 0; m < n; m++) { top += w[m] * g[m] * w[m]; bottom += w[m] * g[m] * x[m]; }
    lambda = top / bottom;
    const norm = Math.hypot(...x);
    w = x.map((v) => v / norm);
  }
  return 1 / Math.sqrt(lambda);
}

function pileCase(kind, { N, latitude, width, days, dt, sst, disturb }) {
  const mesh = caseMesh(N), { nCells: C, nEdges: E, areaCell } = mesh;
  const ocean = createOcean(mesh, recipe());
  const surfaceT = new Float64Array(C).fill(sst), ice = new Float64Array(C), flux = new Float64Array(C), calm = new Float64Array(E);
  ocean.initialize(surfaceT, ice);
  const distance = distancesFrom(mesh, latitude, 0), edgeDistance = distancesFrom(mesh, latitude, 0, 'edge');
  const shape = Float64Array.from(distance, (d) => Math.exp(-((d / width) ** 2)));
  const reshape = (k, amountAt) => {
    for (let i = 0; i < C; i++) {
      const n = k * C + i, t = ocean.Q[n] / ocean.h[n], s = ocean.W[n] / ocean.h[n];
      ocean.h[n] += amountAt(i); ocean.Q[n] = ocean.h[n] * t; ocean.W[n] = ocean.h[n] * s;
    }
  };
  const marker = disturb({ ocean, C, shape, reshape });
  for (let i = 0; i < C; i++) { let sum = 0; for (let k = 0; k < ocean.layers; k++) sum += ocean.h[k * C + i]; ocean.eta[i] = sum - ocean.D[i]; }

  let centre = 0;
  for (let i = 0; i < C; i++) if (distance[i] < distance[centre]) centre = i;
  const far = [], ring = [];
  for (let i = 0; i < C; i++) {
    if (distance[i] > Math.max(4e6, 2.5 * width)) far.push(i);
    if (distance[i] > 0.6 * width && distance[i] < 1.4 * width) ring.push(i);
  }
  const ringEdges = [];
  for (let e = 0; e < E; e++) if (edgeDistance[e] > 0.6 * width && edgeDistance[e] < 1.4 * width) ringEdges.push(e);
  const spacing = meanSpacing(mesh), basis = cellBasis(mesh);
  const axis = [Math.cos(latitude * DEG), 0, Math.sin(latitude * DEG)], normal = [-Math.sin(latitude * DEG), 0, Math.cos(latitude * DEG)];
  const slice = [];
  for (let i = 0; i < C; i++) {
    const x = mesh.xCell[3 * i], y = mesh.xCell[3 * i + 1], z = mesh.xCell[3 * i + 2];
    const along = Math.atan2(y, x * axis[0] + z * axis[2]) * mesh.radius;
    const across = Math.asin(Math.max(-1, Math.min(1, x * normal[0] + z * normal[2]))) * mesh.radius;
    if (Math.abs(across) < 0.5 * spacing && Math.abs(along) < Math.max(4e6, 3 * width)) slice.push([along, i]);
  }
  slice.sort((a, b) => a[0] - b[0]);
  const sliceCells = Int32Array.from(slice, ([, i]) => i), sliceDistance = Float32Array.from(slice, ([d]) => d);

  const meanOf = (field, cells) => { let s = 0, a = 0; for (const i of cells) { s += areaCell[i] * field(i); a += areaCell[i]; } return s / a; };
  const thermocline = (i) => { let z = 0; for (let k = 0; k <= marker; k++) z += ocean.h[k * C + i]; return z; };
  const columnMass = (i) => { let mass = 0; for (let k = 0; k < ocean.layers; k++) { const n = k * C + i, hk = ocean.h[n]; if (hk > 0) mass += seawaterDensity(ocean.Q[n] / hk, ocean.W[n] / hk) * hk; } return mass; };
  const gradEta = new Float64Array(E), tangent = new Float64Array(E), vector = new Float64Array(3 * C);
  const farColumn = far[0], farThickness = Float64Array.from({ length: ocean.layers }, (_, k) => ocean.h[k * C + farColumn]);
  const farDensity = Float64Array.from(farThickness, (hk, k) => seawaterDensity(ocean.Q[k * C + farColumn] / hk, ocean.W[k * C + farColumn] / hk));
  const waveSpeed = internalWaveSpeed(farThickness, farDensity), f = 2 * mesh.omega * Math.sin(latitude * DEG);
  let fixedDepth = 0;
  while (fixedDepth < ocean.layers && farThickness[fixedDepth] - ocean.table.minimum[fixedDepth] < 0.5) fixedDepth++;

  const sim = { kind, mesh, ocean, time: 0, end: days * DAY, nextDt: typeof dt === 'function' ? () => dt(sim.time) : () => dt, centre, distance, marker, waveSpeed };
  sim.step = () => {
    const step = sim.nextDt();
    ocean.advance(surfaceT, ice, flux, calm, step);
    sim.time += step;
    if (sim.afterStep) sim.afterStep();
  };
  sim.geostrophic = () => {
    gradient(mesh, ocean.eta, gradEta);
    tangential(mesh, ocean.u.subarray(0, E), tangent);
    let ab = 0, bb = 0, aa = 0;
    for (const e of ringEdges) { const a = mesh.fEdge[e] * tangent[e], b = G * gradEta[e]; ab += a * b; bb += b * b; aa += a * a; }
    return { ratio: bb > 0 ? ab / bb : NaN, correlation: aa > 0 && bb > 0 ? ab / Math.sqrt(aa * bb) : NaN };
  };
  sim.common = () => {
    const L = ocean.layers, n = sliceCells.length;
    const currentEast = new Float32Array(C), currentNorth = new Float32Array(C);
    components(mesh, ocean.u.subarray(0, E), currentEast, currentNorth, vector);
    const interfaces = new Float32Array(n * L), across = new Float32Array(n * L), eta = new Float32Array(n), thickness = new Float32Array(n * L);
    for (let k = 0; k < L; k++) {
      cellVector(mesh, ocean.u.subarray(k * E, (k + 1) * E), vector);
      sliceCells.forEach((i, s) => { across[s * L + k] = dot(vector, basis.north, i); });
    }
    sliceCells.forEach((i, s) => {
      eta[s] = ocean.eta[i];
      let z = ocean.eta[i];
      for (let k = 0; k < L; k++) { z -= ocean.h[k * C + i]; interfaces[s * L + k] = z; thickness[s * L + k] = ocean.h[k * C + i]; }
    });
    const etaFar = meanOf((i) => ocean.eta[i], far);
    return {
      eta: Float32Array.from(ocean.eta), thermocline: Float32Array.from({ length: C }, (_, i) => thermocline(i)), currentEast, currentNorth,
      rise: ocean.eta[centre] - etaFar, depression: thermocline(centre) - meanOf(thermocline, far), bottomPressure: (columnMass(centre) - meanOf(columnMass, far)) / RHO,
      geostrophic: sim.geostrophic(), layerSpeeds: layerSpeeds(mesh, ocean, ring),
      slice: { distance: Float32Array.from(sliceDistance), eta, interfaces, across, thickness },
    };
  };
  sim.layoutExtra = () => ({ centre, latitudeOfCentre: mesh.latCell[centre] / DEG, width, marker, fixedDepth, waveSpeed, deformationRadius: waveSpeed / Math.abs(f), slice: { cells: Int32Array.from(sliceCells), distance: Float32Array.from(sliceDistance) } });
  return sim;
}

export function lensCase({ N = 16, latitude = 45, depth = 60, width = 1.5e6, days = 30, dt = 3600, sst = 300, light = 1025.25, dense = 1026.77 } = {}) {
  let up = -1, down = -1, pair = null;
  const sim = pileCase('lens', {
    N, latitude, width, days, dt, sst,
    disturb: ({ ocean, C, shape, reshape }) => {
      const labels = Array.from(ocean.table.densities);
      up = labels.indexOf(light); down = labels.indexOf(dense);
      if (up < 0 || down <= up) throw new Error(`the lens moves water from a layer labeled ${dense} to a shallower one labeled ${light}, both among ${labels.join(', ')}`);
      const amount = Float64Array.from(shape, (s, i) => Math.min(depth * s, 0.5 * ocean.h[down * C + i]));
      const own = (k) => seawaterDensity(ocean.Q[k * C] / ocean.h[k * C], ocean.W[k * C] / ocean.h[k * C]);
      pair = [own(up), own(down)];
      reshape(up, (i) => amount[i]);
      reshape(down, (i) => -amount[i]);
      return up;
    },
  });
  const { mesh, ocean } = sim;
  sim.layout = () => layoutOf('lens', N, mesh, ocean, { light: ocean.table.densities[up], dense: ocean.table.densities[down], lens: up, donor: down, ...sim.layoutExtra() });
  sim.frame = () => {
    const common = sim.common();
    return { kind: 'lens', time: sim.time, day: sim.time / DAY, done: sim.time >= sim.end, depth, width, latitude, ...common, lensDensity: pair[0], donorDensity: pair[1], expected: (pair[1] - pair[0]) / RHO * common.depression };
  };
  return sim;
}

export function bumpCase({ N = 16, latitude = 45, height = 1, width = 1.5e6, days = 5, dt = (t) => (t < DAY ? 600 : 3600), sst = 300, threshold = 0.02 } = {}) {
  let deepest = -1;
  const sim = pileCase('bump', {
    N, latitude, width, days, dt, sst,
    disturb: ({ ocean, C, shape, reshape }) => {
      deepest = ocean.layers - 1;
      while (deepest > 0 && ocean.h[deepest * C] < 10) deepest--;
      reshape(deepest, (i) => height * shape[i]);
      return ocean.table.densities.indexOf(1025.25);
    },
  });
  const { mesh, ocean, distance, centre } = sim, { nCells: C, areaCell } = mesh;
  const bin = 0.8 * meanSpacing(mesh), bins = Math.ceil(Math.PI * mesh.radius / bin), profile = new Float64Array(bins), area = new Float64Array(bins);
  let total = 0, sea = 0, depth = 0;
  for (let i = 0; i < C; i++) { area[Math.min(bins - 1, Math.floor(distance[i] / bin))] += areaCell[i]; sea += areaCell[i]; depth += areaCell[i] * ocean.D[i]; }
  for (let i = 0; i < C; i++) total += areaCell[i] * ocean.eta[i];
  const level = total / sea, core = [];
  for (let i = 0; i < C; i++) if (distance[i] < 2 * width) core.push(i);
  const heldIn = () => { let s = 0; for (const i of core) s += areaCell[i] * (ocean.eta[i] - level); return s; };
  const initial = heldIn(), fronts = [];
  const radial = () => {
    profile.fill(0);
    for (let i = 0; i < C; i++) profile[Math.min(bins - 1, Math.floor(distance[i] / bin))] += areaCell[i] * ocean.eta[i];
    for (let b = 0; b < bins; b++) profile[b] = area[b] > 0 ? profile[b] / area[b] : NaN;
  };
  const front = () => {
    radial();
    const edge = threshold * height;
    for (let b = bins - 2; b >= 0; b--) {
      if (!(profile[b] > edge)) continue;
      const beyond = Number.isFinite(profile[b + 1]) ? profile[b + 1] : 0;
      return (b + 0.5 + (profile[b] - edge) / Math.max(1e-12, profile[b] - beyond)) * bin;
    }
    return 0;
  };
  fronts.push([0, front()]);
  sim.afterStep = () => { if (sim.time <= 2 * DAY) fronts.push([sim.time, front()]); };
  const gravitySpeed = Math.sqrt(G * depth / sea);
  sim.layout = () => layoutOf('bump', N, mesh, ocean, { height, bumped: deepest, ...sim.layoutExtra(), bin });
  sim.frame = () => {
    radial();
    const latest = fronts[fronts.length - 1], tracking = latest[1] < 0.8 * Math.PI * mesh.radius && sim.time <= 2 * DAY;
    const frontSpeed = tracking && sim.time >= 3 * 3600 ? fit(fronts.filter(([t]) => t >= 2 * 3600)) : null;
    return {
      kind: 'bump', time: sim.time, day: sim.time / DAY, done: sim.time >= sim.end, height, width, latitude, ...sim.common(),
      centreHeight: ocean.eta[centre] - level, kept: heldIn() / initial, level,
      profile: { distance: Float32Array.from({ length: bins }, (_, b) => (b + 0.5) * bin), eta: Float32Array.from(profile) },
      front: tracking ? latest[1] : null, frontSpeed, gravitySpeed,
    };
  };
  return sim;
}

const BASIN = { west: -40, east: 40, south: 12, north: 48, depth: 4000 };
export const basinElevation = (lat, lon) => (Math.abs(lon) < BASIN.east * DEG && lat > BASIN.south * DEG && lat < BASIN.north * DEG ? -BASIN.depth : 500);
export const basinWind = (lat) => { const l = lat / DEG; return l > 15 && l < 45 ? -Math.cos(Math.PI * (l - 15) / 30) : 0; };

export function basinMesh(N, beta = true, latitude = 30) {
  const mesh = caseMesh(N);
  if (beta) return mesh;
  const f = 2 * mesh.omega * Math.sin(latitude * DEG);
  return { ...mesh, fEdge: new Float64Array(mesh.nEdges).fill(f), fVertex: new Float64Array(mesh.nVertices).fill(f) };
}

export function basinCase({ N = 16, dt = 5400, walls = true, beta = true, tau0 = 0.1, days = 100, sst = 300, ocean: options = {} } = {}) {
  const mesh = basinMesh(N, beta), { nCells: C, nEdges: E, areaCell } = mesh;
  const geography = walls ? createGeography(mesh, syntheticTopography(180, 360, basinElevation), { landBridges: {}, seaStraits: {} }) : null;
  const ocean = createOcean(mesh, recipe({ geography, ...options }));
  const surfaceT = Float64Array.from(mesh.latCell, typeof sst === 'function' ? sst : () => sst), ice = new Float64Array(C), flux = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  const stress = zonalWindOnEdges(mesh, (lat) => tau0 * basinWind(lat));
  const half = 0.5 * meanSpacing(mesh) / mesh.radius / DEG, L = ocean.layers;
  const lat = (i) => mesh.latCell[i] / DEG, lon = (i) => mesh.lonCell[i] / DEG;
  const parallel = [];
  for (let i = 0; i < C; i++) if (ocean.cellOcean[i] && Math.abs(lat(i) - 30) < half && lon(i) >= BASIN.west && lon(i) <= BASIN.east) parallel.push(i);
  parallel.sort((a, b) => mesh.lonCell[a] - mesh.lonCell[b]);
  const beta30 = beta ? 2 * mesh.omega * Math.cos(30 * DEG) / mesh.radius : 0;
  const sverdrup = beta ? -tau0 * Math.PI / (30 * DEG * mesh.radius) / (RHO * beta30) : null;
  const basinWidth = (BASIN.east - BASIN.west) * DEG * mesh.radius * Math.cos(30 * DEG);
  const total = new Float64Array(E), vector = new Float64Array(3 * C);
  const sim = { kind: 'basin', mesh, ocean, time: 0, end: days * DAY, nextDt: () => dt };
  sim.step = () => { ocean.advance(surfaceT, ice, flux, stress, dt); sim.time += dt; };
  const transports = () => {
    for (let e = 0; e < E; e++) {
      let t = 0;
      const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1];
      for (let k = 0; k < L; k++) t += 0.5 * (ocean.h[k * C + a] + ocean.h[k * C + b]) * ocean.u[k * E + e];
      total[e] = ocean.edgeOcean[e] ? t : 0;
    }
    const east = new Float32Array(C), north = new Float32Array(C);
    components(mesh, total, east, north, vector);
    for (let i = 0; i < C; i++) if (!ocean.cellOcean[i]) { east[i] = 0; north[i] = 0; }
    return { east, north };
  };
  const band = (north, lo, hi) => {
    let v = 0, area = 0;
    for (let i = 0; i < C; i++) {
      if (!ocean.cellOcean[i] || lat(i) < 25 || lat(i) > 35 || lon(i) < lo || lon(i) >= hi) continue;
      v += areaCell[i] * north[i]; area += areaCell[i];
    }
    return area > 0 ? v / area : NaN;
  };
  sim.layout = () => layoutOf('basin', N, mesh, ocean, { walls, beta, tau0, parallel: Int32Array.from(parallel), stress: Float32Array.from(mesh.latCell, (l) => tau0 * basinWind(l)) });
  sim.frame = () => {
    const { east: transportEast, north: transportNorth } = transports();
    const currentEast = new Float32Array(C), currentNorth = new Float32Array(C);
    components(mesh, ocean.u.subarray(0, E), currentEast, currentNorth, vector);
    let strongest = null, strongestValue = -Infinity;
    for (let centre = -38; centre <= 38; centre += 4) { const v = band(transportNorth, centre - 2, centre + 2); if (Number.isFinite(v) && v > strongestValue) { strongestValue = v; strongest = centre; } }
    const across = Float32Array.from(parallel, (i) => transportNorth[i]), longitude = Float32Array.from(parallel, (i) => lon(i));
    let northward = 0, southward = 0;
    for (let s = 0; s < parallel.length; s++) {
      const left = s > 0 ? 0.5 * (longitude[s] - longitude[s - 1]) : 0.5 * (longitude[s + 1] - longitude[s]);
      const right = s < parallel.length - 1 ? 0.5 * (longitude[s + 1] - longitude[s]) : left;
      const dx = (left + right) * DEG * mesh.radius * Math.cos(30 * DEG), sv = across[s] * dx / 1e6;
      if (sv > 0) northward += sv; else southward += sv;
    }
    const eta = Float32Array.from(ocean.eta, (v, i) => (ocean.cellOcean[i] ? v : NaN));
    const zonal = { latitude: Float32Array.from({ length: 36 }, (_, b) => -87.5 + 5 * b), east: new Float32Array(36) }, zonalArea = new Float64Array(36);
    for (let i = 0; i < C; i++) {
      if (!ocean.cellOcean[i]) continue;
      const b = Math.min(35, Math.floor((lat(i) + 90) / 5));
      zonal.east[b] += areaCell[i] * transportEast[i]; zonalArea[b] += areaCell[i];
    }
    for (let b = 0; b < 36; b++) zonal.east[b] = zonalArea[b] > 0 ? zonal.east[b] / zonalArea[b] : NaN;
    const diagnostics = ocean.diagnostics();
    return {
      kind: 'basin', time: sim.time, day: sim.time / DAY, done: sim.time >= sim.end, walls, beta, tau0,
      eta, transportEast, transportNorth, currentEast, currentNorth,
      west: band(transportNorth, -40, -32), interior: band(transportNorth, -20, 20), east: band(transportNorth, 32, 40), strongest,
      parallel: { longitude, north: across }, northward, southward, zonal,
      sverdrup, sverdrupTotal: sverdrup === null ? null : sverdrup * basinWidth / 1e6,
      fastest: fastestCurrent(mesh, ocean), limited: diagnostics.oceanLimited,
    };
  };
  return sim;
}

export const CASES = { ekman: ekmanCase, lens: lensCase, bump: bumpCase, basin: basinCase };

export function runCase(sim, days = sim.end / DAY) {
  while (sim.time < days * DAY - 1e-6) sim.step();
  return sim.frame();
}

export const ekman = (options = {}) => runCase(ekmanCase(options));
export const lens = (options = {}) => runCase(lensCase(options));
export const bump = (options = {}) => runCase(bumpCase(options));
export const basin = (options = {}) => runCase(basinCase(options));
