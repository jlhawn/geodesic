/*
 * Small cases of the model's own ocean for Part 8's figures: createOcean,
 * unchanged, on buildMesh(new Grid(N)), or on a copy of that mesh whose
 * fEdge and fVertex hold one value ("constant f"), with the four interior
 * classes of FOUR_CLASSES unless `classes` says otherwise, and with no heat,
 * fresh water or ice. Cell i of every array is cell i of new Grid(N), the
 * order the explainer's globes draw. Units are SI (m, s, m/s, m²/s);
 * vectors come as east and north components; fields are Float32Array and
 * every array a frame returns is a fresh copy. A case factory returns a
 * simulation that steps the ocean and reports what a figure draws; ekman,
 * lens, bump and basin run one to its end and return its last frame.
 */
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { cellVector, gradient, tangential } from '../js/dynamics/operators.module.js';
import { createOcean, depthFields } from '../js/ocean/layered.module.js';
import { labelTemperature, seawaterDensity } from '../js/ocean/seawater.module.js';
import { createGeography, syntheticTopography } from '../js/geography.module.js';

export const RHO = 1025, G = 9.81, DEG = Math.PI / 180, DAY = 86400;
export const FOUR_CLASSES = { densities: [1023, 1025, 1026.5, 1026.95], salinities: [35, 35, 34.86, 34.8], bottoms: [170, 500, 1004, 12000] };

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

export function totalHeatSalt(ocean, m) {
  let heat = 0, salt = 0;
  for (let i = 0; i < m.nCells; i++) {
    if (!ocean.cellOcean[i]) continue;
    for (let k = 0; k < ocean.layers; k++) { heat += m.areaCell[i] * ocean.Q[k * m.nCells + i]; salt += m.areaCell[i] * ocean.W[k * m.nCells + i]; }
  }
  return { heat, salt };
}

export function northwardMixedTransport(ocean, m) {
  const vector = cellVector(m, ocean.u.subarray(0, m.nEdges)), { north } = cellBasis(m);
  const out = new Float64Array(m.nCells);
  for (let i = 0; i < m.nCells; i++) out[i] = ocean.h[i] * dot(vector, north, i);
  return out;
}

export const UNIFORM = { everySteps: 1, thermoclineTilt: 0, salinityProfile: () => 34.5, mixedDepth: 50 };
export function uniformOcean(m, options = {}) {
  const ocean = createOcean(m, { ...UNIFORM, ...options });
  const surfaceT = new Float64Array(m.nCells).fill(278), ice = new Float64Array(m.nCells), flux = new Float64Array(m.nCells);
  ocean.initialize(surfaceT, ice);
  return { ocean, surfaceT, ice, flux, calm: new Float64Array(m.nEdges) };
}

export function deepen({ ocean, surfaceT, ice }, cells, depth, offset) {
  const { h, Q, W, densities, layers: L } = ocean, n = h.length / L;
  for (const i of cells) {
    for (let k = 1; k < L && h[i] < depth; k++) {
      const a = k * n + i, take = Math.min(h[a] - 0.01, depth - h[i]);
      if (take <= 0) continue;
      const f = take / h[a];
      h[i] += take; Q[i] += Q[a] * f; W[i] += W[a] * f;
      h[a] -= take; Q[a] -= Q[a] * f; W[a] -= W[a] * f;
    }
    let below = 1;
    while (below < L - 1 && h[below * n + i] <= 5) below++;
    const t = labelTemperature(densities[below] + offset, W[i] / h[i]);
    Q[i] = h[i] * t; surfaceT[i] = t;
  }
  ocean.load(ocean.serialize(), surfaceT, ice);
}

export function buriedBump(m, options) {
  const C = m.nCells;
  const ocean = createOcean(m, { bathymetry: new Float64Array(C).fill(4000), thermoclineTilt: 0, salinityProfile: () => 35, ...options });
  ocean.initialize(new Float64Array(C).fill(290), new Float64Array(C));
  const distance = distancesFrom(m, -45, 0);
  const reshape = (k, amount) => {
    for (let i = 0; i < C; i++) {
      const n = k * C + i, t = ocean.Q[n] / ocean.h[n], s = ocean.W[n] / ocean.h[n];
      ocean.h[n] += amount * Math.exp(-((distance[i] / 2e6) ** 2));
      ocean.Q[n] = ocean.h[n] * t; ocean.W[n] = ocean.h[n] * s;
    }
  };
  return { ocean, reshape };
}

export function interfaceRange(ocean, k) {
  const C = ocean.h.length / ocean.layers;
  let low = Infinity, high = -Infinity;
  for (let i = 0; i < C; i++) {
    let z = 0;
    for (let j = 0; j <= k; j++) z += ocean.h[j * C + i];
    low = Math.min(low, z); high = Math.max(high, z);
  }
  return high - low;
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

function recipe(classes, extra = {}) {
  return { everySteps: 1, thermoclineTilt: 0, salinityProfile: () => 35, ...(classes ?? {}), ...extra };
}

function layoutOf(kind, N, mesh, ocean, extra = {}) {
  return {
    kind, N, nCells: mesh.nCells, layers: ocean.layers, densities: Float64Array.from(ocean.densities),
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

// The first layer below the mixed layer: the one holding the mixed-layer depth, or the next if that one ends there.
function firstBelow(ocean, C, i) {
  const depth = ocean.mixedDepths[i];
  let k = 0, top = 0;
  while (k < ocean.layers - 1 && top + ocean.h[k * C + i] <= depth + 1e-9) { top += ocean.h[k * C + i]; k++; }
  return Math.max(1, k);
}
// The mixed layer's mean velocity: the thickness-weighted mean of the layers within its depth.
function mixedVectors(mesh, ocean, layerVectors, out) {
  const C = mesh.nCells;
  out.fill(0);
  for (let i = 0; i < C; i++) {
    const depth = ocean.mixedDepths[i];
    let top = 0, within = 0;
    for (let k = 0; k < ocean.layers && top < depth; k++) {
      const part = Math.min(ocean.h[k * C + i], depth - top);
      for (let c = 0; c < 3; c++) out[3 * i + c] += part * layerVectors[k][3 * i + c];
      within += part; top += ocean.h[k * C + i];
    }
    if (within > 0) for (let c = 0; c < 3; c++) out[3 * i + c] /= within;
  }
}

function drift(east, north, theory) {
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

export function ekmanCase({ N = 8, tau = 0.1, latitude = 45, days = 4, dt = 1350, view = 'column', wind = 'trades', averaging = 2 * DAY, sst = 290, classes = FOUR_CLASSES } = {}) {
  const mesh = caseMesh(N), { nCells: C, nEdges: E, areaCell } = mesh;
  const ocean = createOcean(mesh, recipe(classes));
  const surfaceT = new Float64Array(C).fill(sst), ice = new Float64Array(C), flux = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  if (view !== 'column' && view !== 'planet') throw new Error(`the Ekman view is 'column' or 'planet', not ${view}`);
  if (view === 'planet' && !PLANET_WINDS[wind]) throw new Error(`the planet's wind is ${Object.keys(PLANET_WINDS).join(' or ')}, not ${wind}`);
  const shape = view === 'planet' ? PLANET_WINDS[wind] : columnWind(latitude);
  const stress = zonalWindOnEdges(mesh, (lat) => tau * shape(lat));
  const basis = cellBasis(mesh), half = 0.5 * meanSpacing(mesh) / mesh.radius / DEG;
  const coriolis = (lat) => 2 * mesh.omega * Math.sin(lat);
  const ml = new Float64Array(3 * C), beneath = new Float64Array(3 * C), layerVectors = Array.from({ length: ocean.layers }, () => new Float64Array(3 * C));
  const allLayers = () => { for (let k = 0; k < ocean.layers; k++) cellVector(mesh, ocean.u.subarray(k * E, (k + 1) * E), layerVectors[k]); };
  const sim = { kind: 'ekman', view, mesh, ocean, time: 0, end: days * DAY, nextDt: () => dt };
  let advanceExtra = () => {}, frameExtra = () => ({}), layoutExtra = () => ({});

  if (view === 'column') {
    const ring = [];
    let ringArea = 0, probe = 0, best = Infinity, theory = 0;
    for (let i = 0; i < C; i++) {
      const lat = mesh.latCell[i] / DEG, lon = mesh.lonCell[i] / DEG;
      if (Math.abs(lat - latitude) < half) { ring.push(i); ringArea += areaCell[i]; }
      const miss = Math.abs(lat - latitude) + 0.01 * Math.abs(lon);
      if (miss < best) { best = miss; probe = i; }
    }
    for (const i of ring) theory += areaCell[i] / ringArea * tau * shape(mesh.latCell[i]) / (RHO * Math.abs(coriolis(mesh.latCell[i])));
    const f = coriolis(latitude * DEG), period = 2 * Math.PI / Math.abs(f);
    const series = { time: [0], east: [0], north: [0], relativeEast: [0], relativeNorth: [0], cellEast: [0], cellNorth: [0], mixedDepth: [], transportEast: [0], transportNorth: [0], ekmanEast: [0], ekmanNorth: [0] };
    const record = () => {
      allLayers();
      mixedVectors(mesh, ocean, layerVectors, ml);
      let e0 = 0, n0 = 0, re = 0, rn = 0, depth = 0, te = 0, tn = 0, ke = 0, kn = 0;
      for (const i of ring) {
        const w = areaCell[i] / ringArea, h = ocean.mixedDepths[i], k = firstBelow(ocean, C, i);
        for (let c = 0; c < 3; c++) beneath[3 * i + c] = layerVectors[k][3 * i + c];
        const ue = dot(ml, basis.east, i), un = dot(ml, basis.north, i), be = dot(beneath, basis.east, i), bn = dot(beneath, basis.north, i);
        e0 += w * ue; n0 += w * un; re += w * (ue - be); rn += w * (un - bn); depth += w * h;
        te += w * h * ue; tn += w * h * un; ke += w * h * (ue - be); kn += w * h * (un - bn);
      }
      return { e0, n0, re, rn, depth, te, tn, ke, kn, ce: dot(ml, basis.east, probe), cn: dot(ml, basis.north, probe) };
    };
    const first = record();
    series.mixedDepth.push(first.depth);
    advanceExtra = () => {
      const r = record();
      series.time.push(sim.time); series.east.push(r.e0); series.north.push(r.n0); series.relativeEast.push(r.re); series.relativeNorth.push(r.rn);
      series.cellEast.push(r.ce); series.cellNorth.push(r.cn); series.mixedDepth.push(r.depth);
      series.transportEast.push(r.te); series.transportNorth.push(r.tn); series.ekmanEast.push(r.ke); series.ekmanNorth.push(r.kn);
    };
    const below = () => {
      const east = new Float32Array(ocean.layers), north = new Float32Array(ocean.layers), vector = new Float64Array(3 * C);
      for (let k = 0; k < ocean.layers; k++) {
        cellVector(mesh, ocean.u.subarray(k * E, (k + 1) * E), vector);
        let se = 0, sn = 0, area = 0;
        for (const i of ring) if (ocean.h[k * C + i] > 0) { se += areaCell[i] * dot(vector, basis.east, i); sn += areaCell[i] * dot(vector, basis.north, i); area += areaCell[i]; }
        east[k] = area > 0 ? se / area : NaN; north[k] = area > 0 ? sn / area : NaN;
      }
      return { east, north };
    };
    frameExtra = () => {
      const n = series.time.length, from = sim.time - 2 * period;
      let te = 0, tn = 0, ke = 0, kn = 0, count = 0, earliest = sim.time;
      for (let s = n - 1; s >= 1 && series.time[s] > from; s--) { te += series.transportEast[s]; tn += series.transportNorth[s]; ke += series.ekmanEast[s]; kn += series.ekmanNorth[s]; earliest = series.time[s - 1]; count++; }
      const mean = (v) => (count ? v / count : 0);
      const copy = (name, type = Float32Array) => type.from(series[name]);
      return {
        latitude, cellLatitude: mesh.latCell[probe] / DEG, f, inertialPeriod: period, theory, ringCells: ring.length,
        mixedDepth: series.mixedDepth[n - 1], averagedOver: sim.time - earliest, settled: sim.time >= 2 * period,
        transport: drift(mean(te), mean(tn), theory), ekman: drift(mean(ke), mean(kn), theory), below: below(),
        series: { time: copy('time', Float64Array), east: copy('east'), north: copy('north'), relativeEast: copy('relativeEast'), relativeNorth: copy('relativeNorth'), cellEast: copy('cellEast'), cellNorth: copy('cellNorth'), mixedDepth: copy('mixedDepth') },
      };
    };
    layoutExtra = () => ({ ring: Int32Array.from(ring), probe });
  } else {
    const capacity = Math.max(1, Math.ceil(averaging / dt)), width = 5 * C, snapshots = new Float32Array(capacity * width), sums = new Float64Array(width);
    let count = 0, next = 0;
    const meanDepth = () => { let s = 0, a = 0; for (let i = 0; i < C; i++) { s += areaCell[i] * ocean.mixedDepths[i]; a += areaCell[i]; } return s / a; };
    advanceExtra = () => {
      allLayers();
      mixedVectors(mesh, ocean, layerVectors, ml);
      const w = depthFields(mesh, ocean.layers, { h: ocean.h, u: ocean.u, temperature: () => 0, cellOcean: ocean.cellOcean }, meanDepth(), { thin: 0 }).upwelling;
      const slot = snapshots.subarray(next * width, (next + 1) * width);
      if (count === capacity) for (let j = 0; j < width; j++) sums[j] -= slot[j];
      for (let i = 0; i < C; i++) {
        const h = ocean.mixedDepths[i], k = firstBelow(ocean, C, i);
        const ue = dot(ml, basis.east, i), un = dot(ml, basis.north, i);
        slot[i] = h * ue; slot[C + i] = h * un;
        slot[2 * C + i] = h * (ue - dot(layerVectors[k], basis.east, i)); slot[3 * C + i] = h * (un - dot(layerVectors[k], basis.north, i));
        slot[4 * C + i] = w[i];
      }
      for (let j = 0; j < width; j++) sums[j] += slot[j];
      next = (next + 1) % capacity; count = Math.min(capacity, count + 1);
    };
    const bins = 36;
    frameExtra = () => {
      const mean = (part) => Float32Array.from({ length: C }, (_, i) => (count ? sums[part * C + i] / count : 0));
      const transportEast = mean(0), transportNorth = mean(1), ekmanEast = mean(2), ekmanNorth = mean(3), upwelling = mean(4);
      const zonal = { latitude: new Float32Array(bins), stress: new Float32Array(bins), theory: new Float32Array(bins), transportNorth: new Float32Array(bins), ekmanNorth: new Float32Array(bins), upwelling: new Float32Array(bins) };
      const area = new Float64Array(bins);
      let equator = 0, equatorArea = 0;
      for (let i = 0; i < C; i++) {
        const lat = mesh.latCell[i] / DEG, b = Math.min(bins - 1, Math.floor((lat + 90) / 5)), a = areaCell[i];
        zonal.transportNorth[b] += a * transportNorth[i]; zonal.ekmanNorth[b] += a * ekmanNorth[i]; zonal.upwelling[b] += a * upwelling[i]; area[b] += a;
        if (Math.abs(lat) < half) { equator += a * upwelling[i]; equatorArea += a; }
      }
      for (let b = 0; b < bins; b++) {
        const lat = (-87.5 + 5 * b) * DEG;
        zonal.latitude[b] = lat / DEG; zonal.stress[b] = tau * shape(lat); zonal.theory[b] = -tau * shape(lat) / (RHO * coriolis(lat));
        zonal.transportNorth[b] /= area[b]; zonal.ekmanNorth[b] /= area[b]; zonal.upwelling[b] /= area[b];
      }
      return {
        wind, transportEast, transportNorth, ekmanEast, ekmanNorth, upwelling, zonal, equatorUpwelling: equator / equatorArea,
        eta: Float32Array.from(ocean.eta), mixedDepth: Float32Array.from(ocean.mixedDepths),
        stress: Float32Array.from(mesh.latCell, (lat) => tau * shape(lat)), averagedOver: count * dt,
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

function pileCase(kind, { N, latitude, width, days, dt, sst, classes, disturb }) {
  const mesh = caseMesh(N), { nCells: C, nEdges: E, areaCell } = mesh;
  const ocean = createOcean(mesh, recipe(classes));
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
  const thermoclineClass = disturb({ ocean, C, shape, reshape }) ?? 1;
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
  const boundary = 0.5 * (ocean.densities[thermoclineClass] + ocean.densities[Math.min(ocean.layers - 1, thermoclineClass + 1)]);
  const thermocline = (i) => { let z = 0; for (let k = 0; k < ocean.layers; k++) { const n = k * C + i, hk = ocean.h[n]; if (hk > 1e-6 && seawaterDensity(ocean.Q[n] / hk, ocean.W[n] / hk) >= boundary) return z; z += hk; } return z; };
  const gradEta = new Float64Array(E), tangent = new Float64Array(E), vector = new Float64Array(3 * C);

  const sim = { kind, mesh, ocean, time: 0, end: days * DAY, nextDt: typeof dt === 'function' ? () => dt(sim.time) : () => dt, centre, distance };
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
    const interfaces = new Float32Array(n * L), across = new Float32Array(n * L), eta = new Float32Array(n);
    for (let k = 0; k < L; k++) {
      cellVector(mesh, ocean.u.subarray(k * E, (k + 1) * E), vector);
      sliceCells.forEach((i, s) => { across[s * L + k] = dot(vector, basis.north, i); });
    }
    sliceCells.forEach((i, s) => {
      eta[s] = ocean.eta[i];
      let z = ocean.eta[i];
      for (let k = 0; k < L; k++) { z -= ocean.h[k * C + i]; interfaces[s * L + k] = z; }
    });
    const etaFar = meanOf((i) => ocean.eta[i], far);
    return {
      eta: Float32Array.from(ocean.eta), thermocline: Float32Array.from({ length: C }, (_, i) => thermocline(i)), currentEast, currentNorth,
      rise: ocean.eta[centre] - etaFar, depression: thermocline(centre) - meanOf(thermocline, far),
      geostrophic: sim.geostrophic(), layerSpeeds: layerSpeeds(mesh, ocean, ring),
      slice: { distance: Float32Array.from(sliceDistance), eta, interfaces, across },
    };
  };
  sim.layoutExtra = () => ({ centre, latitudeOfCentre: mesh.latCell[centre] / DEG, width, slice: { cells: Int32Array.from(sliceCells), distance: Float32Array.from(sliceDistance) } });
  return sim;
}

export function lensCase({ N = 16, latitude = 45, depth = 60, width = 1.5e6, days = 30, dt = 3600, sst = 300, classes = FOUR_CLASSES, light = 1023, dense = 1026.77 } = {}) {
  let up = -1, down = -1;
  const sim = pileCase('lens', {
    N, latitude, width, days, dt, sst, classes,
    disturb: ({ ocean, C, shape, reshape }) => {
      up = ocean.densities.indexOf(light, 1); down = ocean.densities.indexOf(dense, 1);
      if (up < 1 || down < 1) throw new Error(`the lens moves water between two of the classes ${ocean.densities.slice(1).join(', ')}`);
      const amount = Float64Array.from(shape, (s, i) => Math.min(depth * s, ocean.h[down * C + i] - 1));
      reshape(up, (i) => amount[i]);
      reshape(down, (i) => -amount[i]);
      for (let i = 0; i < C; i++) { const n = up * C + i; if (amount[i] > 0) ocean.Q[n] = ocean.h[n] * labelTemperature(light, ocean.W[n] / ocean.h[n]); }
      return up;
    },
  });
  const { mesh, ocean } = sim, step = ocean.densities[down] - ocean.densities[up];
  sim.layout = () => layoutOf('lens', N, mesh, ocean, { light, dense, ...sim.layoutExtra() });
  sim.frame = () => {
    const common = sim.common();
    return { kind: 'lens', time: sim.time, day: sim.time / DAY, done: sim.time >= sim.end, depth, width, latitude, ...common, expected: step / RHO * common.depression };
  };
  return sim;
}

export function bumpCase({ N = 16, latitude = 45, height = 1, width = 1.5e6, days = 5, dt = (t) => (t < DAY ? 600 : 3600), sst = 300, classes = FOUR_CLASSES, threshold = 0.02 } = {}) {
  const sim = pileCase('bump', {
    N, latitude, width, days, dt, sst, classes,
    disturb: ({ ocean, shape, reshape }) => { reshape(ocean.layers - 1, (i) => height * shape[i]); },
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
  sim.layout = () => layoutOf('bump', N, mesh, ocean, { height, ...sim.layoutExtra(), bin });
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

export function basinCase({ N = 16, dt = 5400, walls = true, beta = true, tau0 = 0.1, days = 100, classes = FOUR_CLASSES, sst = 300 } = {}) {
  const mesh = basinMesh(N, beta), { nCells: C, nEdges: E, areaCell } = mesh;
  const geography = walls ? createGeography(mesh, syntheticTopography(180, 360, basinElevation), { landBridges: {}, seaStraits: {} }) : null;
  const ocean = createOcean(mesh, recipe(classes, { geography }));
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
