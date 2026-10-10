// Precomputes the layered-model frames that explain/sphere.html and
// explain/primitive.html play back:
//   node scripts/explainLayers.mjs jw06|heldsuarez|tracks|state [OUT]
// writing explain/data/<case>_N16.bin, jw06state_N16.bin for state.
import { writeFileSync, mkdirSync } from 'node:fs';
import { dirname } from 'node:path';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createSigmaCore } from '../js/dynamics/sigmaCore.module.js';
import { createRK4Arrays } from '../js/dynamics/integrators.module.js';
import { EARTH, edgeNormalVelocity, hyperdiffusion, cellVelocity } from '../explain/swCases.module.js';

const which = process.argv[2], N = 16, DAY = 86400, dt = 450 * 16 / N;
const out = process.argv[3] ?? new URL(`../explain/data/${which === 'state' ? 'jw06state' : which}_N${N}.bin`, import.meta.url).pathname;
const { a, omega, g } = EARTH, R = 287, cp = 1004.5, p0 = 1e5;
const mesh = buildMesh(new Grid(N), { radius: a, omega }), C = mesh.nCells, E = mesh.nEdges;

function jw06() {
  const eta0 = 0.252, u0 = 35, T0 = 288, gamma = 0.005, deltaT = 4.8e5, etaT = 0.2;
  const etaV = (eta) => (eta - eta0) * Math.PI / 2;
  const shape = (lat) => { const s = Math.sin(lat), c = Math.cos(lat); return { A: -2 * s ** 6 * (c * c + 1 / 3) + 10 / 63, B: 8 / 5 * c ** 3 * (s * s + 2 / 3) - Math.PI / 4 }; };
  const zonalWind = (eta, lat) => u0 * Math.cos(etaV(eta)) ** 1.5 * Math.sin(2 * lat) ** 2;
  const geopotential = (eta, lat) => {
    const exp = R * gamma / g;
    let mean = T0 * g / gamma * (1 - eta ** exp);
    if (eta < etaT) mean -= R * deltaT * ((Math.log(eta / etaT) + 137 / 60) * etaT ** 5 - 5 * etaT ** 4 * eta + 5 * etaT ** 3 * eta ** 2 - (10 / 3) * etaT ** 2 * eta ** 3 + (5 / 4) * etaT * eta ** 4 - eta ** 5 / 5);
    const cv = Math.cos(etaV(eta)) ** 1.5, { A, B } = shape(lat);
    return mean + u0 * cv * (A * u0 * cv + B * a * omega);
  };
  const perturbation = (lon, lat) => {
    const lonC = Math.PI / 9, latC = 2 * Math.PI / 9, radius = a / 10;
    const r = a * Math.acos(Math.min(1, Math.sin(latC) * Math.sin(lat) + Math.cos(latC) * Math.cos(lat) * Math.cos(lon - lonC)));
    return Math.exp(-((r / radius) ** 2));
  };
  const surface = Float64Array.from(mesh.latCell, (lat) => geopotential(1, lat));
  const core = createSigmaCore(mesh, { g, cp, R, p0, surfaceGeopotential: surface, nu4: hyperdiffusion(mesh, 3), nu4Theta: hyperdiffusion(mesh, 3) });
  const K = core.K, pi = new Float64Array(C).fill(p0), theta = new Float64Array(K * C);
  core.diagnose(pi, theta);
  const { exnerLayer, exnerLower } = core.arrays;
  for (let i = 0; i < C; i++) {
    let below = surface[i], thetaBelow = 0;
    for (let k = K - 1; k >= 0; k--) {
      const idx = k * C + i, target = geopotential(core.sigmaMid[k], mesh.latCell[i]);
      const fromBelow = k === K - 1 ? 0 : cp * thetaBelow * (exnerLayer[idx + C] - exnerLower[idx]);
      theta[idx] = (target - below - fromBelow) / (cp * (exnerLower[idx] - exnerLayer[idx]));
      below = target; thetaBelow = theta[idx];
    }
  }
  const u = new Float64Array(K * E);
  for (let k = 0; k < K; k++) u.set(edgeNormalVelocity(mesh, (lon, lat) => ({ zonal: zonalWind(core.sigmaMid[k], lat) + perturbation(lon, lat), meridional: 0 })), k * E);
  return { core, state: [pi, theta, u], tendency: core.tendency, surface };
}

function heldSuarez() {
  const core = createSigmaCore(mesh, { g, cp, R, p0, nu4: hyperdiffusion(mesh, 3), nu4Theta: hyperdiffusion(mesh, 3) });
  const { K, sigmaMid, kappa, exnerLayer } = core.diagnostics;
  const kf = 1 / DAY, ka = 1 / (40 * DAY), ks = 1 / (4 * DAY), sigmaB = 0.7, deltaTy = 60, deltaThetaZ = 10, tMin = 200, tMax = 315;
  const cosLat = Float64Array.from(mesh.latCell, Math.cos), sinLat = Float64Array.from(mesh.latCell, Math.sin);
  const forcing = (state, outs) => {
    const [pi, theta, u] = state, [, dTheta, dU] = outs;
    for (let k = 0; k < K; k++) {
      const s = sigmaMid[k], weight = Math.max(0, (s - sigmaB) / (1 - sigmaB));
      for (let i = 0; i < C; i++) {
        const idx = k * C + i, p = pi[i] * s, c2 = cosLat[i] * cosLat[i];
        const tEq = Math.max(tMin, (tMax - deltaTy * sinLat[i] * sinLat[i] - deltaThetaZ * Math.log(p / p0) * c2) * Math.pow(p / p0, kappa));
        dTheta[idx] -= (ka + (ks - ka) * weight * c2 * c2) * (theta[idx] - tEq / exnerLayer[idx]);
      }
      if (weight > 0) for (let e = 0; e < E; e++) dU[k * E + e] -= kf * weight * u[k * E + e];
    }
  };
  const tendency = (state, outs) => { core.tendency(state, outs); forcing(state, outs); };
  const pi = new Float64Array(C).fill(p0), theta = new Float64Array(K * C);
  core.diagnose(pi, theta);
  let seed = 12345;
  const noise = () => { seed = (seed * 1664525 + 1013904223) >>> 0; return seed / 2 ** 32 - 0.5; };
  for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) {
    const p = pi[i] * sigmaMid[k], c2 = cosLat[i] ** 2;
    const tEq = Math.max(tMin, (tMax - deltaTy * sinLat[i] ** 2 - deltaThetaZ * Math.log(p / p0) * c2) * (p / p0) ** kappa);
    theta[k * C + i] = tEq / exnerLayer[k * C + i] + 0.1 * noise();
  }
  return { core, state: [pi, theta, new Float64Array(K * E)], tendency, surface: null };
}

const run = which === 'jw06' || which === 'tracks' || which === 'state' ? jw06() : which === 'heldsuarez' ? heldSuarez() : null;
if (!run) { console.error('usage: node scripts/explainLayers.mjs jw06|heldsuarez|tracks|state [OUT]'); process.exit(1); }
const { core, state, tendency } = run, K = core.K, step = createRK4Arrays([C, K * C, K * E]);
const [pi, theta, u] = state, { exnerLayer } = core.arrays, sigmaMid = core.sigmaMid;

function columnAt(i, pressure, values) {
  for (let k = K - 1; k > 0; k--) {
    const pk = sigmaMid[k] * pi[i], pu = sigmaMid[k - 1] * pi[i];
    if (pressure <= pk && pressure >= pu) { const w = Math.log(pk / pressure) / Math.log(pk / pu); return values(k) * (1 - w) + values(k - 1) * w; }
  }
  return values(pressure > sigmaMid[K - 1] * pi[i] ? K - 1 : 0);
}
const velocity = [0, 0];
function frame() {
  core.diagnose(pi, theta);
  const ps = new Int16Array(C), t850 = new Int16Array(C), wind250 = new Int16Array(C);
  for (let i = 0; i < C; i++) {
    ps[i] = Math.round((pi[i] / 100 - 1000) * 20);
    t850[i] = Math.round((columnAt(i, 85000, (k) => theta[k * C + i] * exnerLayer[k * C + i]) - 260) * 50);
    wind250[i] = Math.round(columnAt(i, 25000, (k) => { cellVelocity(mesh, u, i, velocity, 0, k * E); return Math.hypot(velocity[0], velocity[1]); }) * 100);
  }
  return [ps, t850, wind250];
}

const BINS = 36;
const means = { u: new Float64Array(BINS * K), T: new Float64Array(BINS * K), count: new Float64Array(BINS * K), samples: 0 };
function accumulate() {
  core.diagnose(pi, theta);
  for (let i = 0; i < C; i++) {
    const bin = Math.min(BINS - 1, Math.floor((mesh.latCell[i] + Math.PI / 2) / Math.PI * BINS));
    for (let k = 0; k < K; k++) {
      cellVelocity(mesh, u, i, velocity, 0, k * E);
      means.u[bin * K + k] += velocity[0] * mesh.areaCell[i]; means.T[bin * K + k] += theta[k * C + i] * exnerLayer[k * C + i] * mesh.areaCell[i]; means.count[bin * K + k] += mesh.areaCell[i];
    }
  }
  means.samples++;
}

if (which === 'tracks') {
  const { maxEdges, nEdgesOnCell, cellsOnCell, verticesOnCell, cellsOnVertex, cellsOnEdge, dcEdge, dvEdge, areaCell, xCell } = mesh, dSigma = core.diagnostics.dSigma, levels = core.levels;
  const surface = run.surface;
  const velocities = new Float64Array(2 * K * C), sigmaDot = new Float64Array((K + 1) * C), v2 = [0, 0];
  const scratch = [new Float64Array(C), new Float64Array(K * C), new Float64Array(K * E)];
  const parcels = [];
  for (const sigma of [0.9, 0.7]) for (let latDeg = 32; latDeg <= 62; latDeg += 6) for (let lonDeg = 80; lonDeg <= 230; lonDeg += 10) parcels.push({ lat: latDeg * Math.PI / 180, lon: lonDeg * Math.PI / 180, sigma, cell: 0 });
  const energy = [], times = [], START = 6 * DAY, END = 11 * DAY;
  const nearest = (x, y, z, start) => { let best = start, bestDot = x * xCell[3 * start] + y * xCell[3 * start + 1] + z * xCell[3 * start + 2]; for (;;) { let next = best; for (let m = 0; m < nEdgesOnCell[best]; m++) { const j = cellsOnCell[maxEdges * best + m], d = x * xCell[3 * j] + y * xCell[3 * j + 1] + z * xCell[3 * j + 2]; if (d > bestDot) { bestDot = d; next = j; } } if (next === best) return best; best = next; } };
  const triple = (x, y, z, i, j) => { const ax = xCell[3 * i], ay = xCell[3 * i + 1], az = xCell[3 * i + 2], bx = xCell[3 * j], by = xCell[3 * j + 1], bz = xCell[3 * j + 2]; return x * (ay * bz - az * by) + y * (az * bx - ax * bz) + z * (ax * by - ay * bx); };
  function weights(p) {
    const x = Math.cos(p.lat) * Math.cos(p.lon), y = Math.cos(p.lat) * Math.sin(p.lon), z = Math.sin(p.lat);
    p.cell = nearest(x, y, z, p.cell);
    let best = null, bestMin = -Infinity;
    for (let m = 0; m < nEdgesOnCell[p.cell]; m++) {
      const v = verticesOnCell[maxEdges * p.cell + m], [i, j, k] = [cellsOnVertex[3 * v], cellsOnVertex[3 * v + 1], cellsOnVertex[3 * v + 2]];
      const w = [triple(x, y, z, j, k), triple(x, y, z, k, i), triple(x, y, z, i, j)], total = w[0] + w[1] + w[2];
      const smallest = Math.min(...w) / total;
      if (smallest > bestMin) { bestMin = smallest; best = [[i, w[0] / total], [j, w[1] / total], [k, w[2] / total]]; }
    }
    return best;
  }
  function sample(p) {
    const ws = weights(p), mids = core.sigmaMid;
    let k = 0; while (k < K - 2 && p.sigma > mids[k + 1]) k++;
    const f = Math.min(1, Math.max(0, (p.sigma - mids[k]) / (mids[k + 1] - mids[k])));
    let u = 0, v = 0, th = 0, ps = 0, sd = 0;
    let kk = 0; while (kk < K - 1 && p.sigma > levels[kk + 1]) kk++;
    const g2 = Math.min(1, Math.max(0, (p.sigma - levels[kk]) / (levels[kk + 1] - levels[kk])));
    for (const [i, w] of ws) {
      u += w * ((1 - f) * velocities[2 * (k * C + i)] + f * velocities[2 * ((k + 1) * C + i)]);
      v += w * ((1 - f) * velocities[2 * (k * C + i) + 1] + f * velocities[2 * ((k + 1) * C + i) + 1]);
      th += w * ((1 - f) * theta[k * C + i] + f * theta[(k + 1) * C + i]);
      ps += w * pi[i];
      sd += w * ((1 - g2) * sigmaDot[kk * C + i] + g2 * sigmaDot[(kk + 1) * C + i]);
    }
    return { u, v, theta: th, ps, sigmaDot: sd };
  }
  function measure(time) {
    let kinetic = 0, internal = 0, mass = 0, thetaMass = 0;
    core.diagnose(pi, theta);
    for (let k = 0; k < K; k++) for (let e = 0; e < E; e++) { const i = cellsOnEdge[2 * e], j = cellsOnEdge[2 * e + 1], uu = u[k * E + e] ** 2; kinetic += 0.25 * dcEdge[e] * dvEdge[e] * uu * dSigma[k] * (pi[i] + pi[j]) / g; }
    for (let i = 0; i < C; i++) { mass += areaCell[i] * pi[i] / g; internal += areaCell[i] * surface[i] * pi[i] / g; for (let k = 0; k < K; k++) { const layer = areaCell[i] * pi[i] * dSigma[k] / g; internal += layer * cp * theta[k * C + i] * exnerLayer[k * C + i]; thetaMass += layer * theta[k * C + i]; } }
    const area = 4 * Math.PI * a * a;
    energy.push({ day: +(time / DAY).toFixed(3), kinetic: kinetic / area, internal: internal / area, mass: mass / area, thetaMass: thetaMass / area });
  }
  const steps = Math.round(12 * DAY / dt), records = [], t0 = performance.now();
  for (let n = 0; n <= steps; n++) {
    const time = n * dt;
    if (n % Math.round(3 * 3600 / dt) === 0) measure(time);
    if (time >= START - 1 && time <= END + 1) {
      core.diagnose(pi, theta);
      for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) { cellVelocity(mesh, u, i, v2, 0, k * E); velocities[2 * (k * C + i)] = v2[0]; velocities[2 * (k * C + i) + 1] = v2[1]; }
      core.tendency(state, scratch);
      for (let k = 0; k <= K; k++) for (let i = 0; i < C; i++) sigmaDot[k * C + i] = core.arrays.piSigmaDot[k * C + i] / pi[i];
      if (Math.round((time - START) / dt) % Math.round(3600 / dt) === 0) {
        times.push(time);
        records.push(parcels.map((p) => { const s = sample(p); return [p.lat * 180 / Math.PI, p.lon * 180 / Math.PI, p.sigma * s.ps / 100, s.theta]; }));
      }
      if (time < END) for (const p of parcels) {
        const s = sample(p);
        p.lat += s.v / a * dt; p.lon += s.u / (a * Math.cos(p.lat)) * dt; p.sigma = Math.min(0.995, Math.max(0.05, p.sigma + s.sigmaDot * dt));
      }
    }
    if (n % Math.round(DAY / dt) === 0) console.error(`tracks day ${(time / DAY).toFixed(0)}: ${((performance.now() - t0) / 1000).toFixed(0)} s`);
    if (n < steps) step(tendency, state, dt);
  }
  const header = { version: 1, case: 'tracks', N, K, parcels: parcels.length, times, fields: ['lat', 'lon', 'pressure', 'theta'], energy };
  const json = Buffer.from(JSON.stringify(header)), pad = (4 - ((8 + json.length) % 4)) % 4, data = new Float32Array(records.length * parcels.length * 4);
  records.forEach((r, h) => r.forEach((v, q) => data.set(v, (h * parcels.length + q) * 4)));
  const head = Buffer.alloc(8); head.write('EXT1', 0); head.writeUInt32LE(json.length + pad, 4);
  mkdirSync(dirname(out), { recursive: true });
  writeFileSync(out, Buffer.concat([head, json, Buffer.alloc(pad, 32), Buffer.from(data.buffer)]));
  console.error(`wrote ${out}: ${parcels.length} parcels, ${records.length} hours, ${energy.length} energy samples`);
  process.exit(0);
}

if (which === 'state') {
  const DAYS = 9, steps = Math.round(DAYS * DAY / dt), t0 = performance.now();
  for (let n = 1; n <= steps; n++) {
    step(tendency, state, dt);
    if (n % Math.round(DAY / dt) === 0) console.error(`state day ${(n * dt / DAY).toFixed(0)}: ${((performance.now() - t0) / 1000).toFixed(0)} s`);
  }
  const arrays = [['pi', 'Pa', ['C'], 'i', pi], ['theta', 'K', ['K', 'C'], 'k * C + i', theta], ['u', 'm/s', ['K', 'E'], 'k * E + e', u], ['surfaceGeopotential', 'm²/s²', ['C'], 'i', run.surface]];
  let offset = 0;
  const layout = arrays.map(([name, unit, shape, index, values]) => { const entry = { name, unit, shape, index, offset, length: values.length }; offset += values.length; return entry; });
  const header = {
    version: 1, case: 'jw06', N, K, C, E, day: DAYS, time: steps * dt, dt, levels: Array.from(core.levels), R, cp, p0, g, a, omega, nu4: core.nu4, nu4Theta: core.nu4Theta,
    data: 'float32, little-endian, the arrays one after another; offset and length count floats from the start of the data, after the header; layer k = 0 is the top; cells and edges in buildMesh(new Grid(N)) order; u is the wind across each edge, along its normal nEdge', arrays: layout,
  };
  const json = Buffer.from(JSON.stringify(header)), pad = (4 - ((8 + json.length) % 4)) % 4, data = new Float32Array(offset);
  for (const [n, entry] of layout.entries()) data.set(arrays[n][4], entry.offset);
  const head = Buffer.alloc(8); head.write('EXS1', 0); head.writeUInt32LE(json.length + pad, 4);
  mkdirSync(dirname(out), { recursive: true });
  writeFileSync(out, Buffer.concat([head, json, Buffer.alloc(pad, 32), Buffer.from(data.buffer)]));
  let low = 0; for (let i = 0; i < C; i++) if (pi[i] < pi[low]) low = i;
  console.error(`wrote ${out}: day ${DAYS}, ${offset} floats, deepest low ${(pi[low] / 100).toFixed(1)} hPa at ${(mesh.latCell[low] * 180 / Math.PI).toFixed(1)}° ${(mesh.lonCell[low] * 180 / Math.PI).toFixed(1)}°, ${((performance.now() - t0) / 1000).toFixed(0)} s`);
  process.exit(0);
}

const plan = which === 'jw06'
  ? { days: 12, sequences: [{ name: 'wave', every: 3 * 3600, from: 0, to: 12 * DAY }] }
  : { days: 200, sequences: [{ name: 'spinup', every: DAY, from: 0, to: 120 * DAY }, { name: 'settled', every: 6 * 3600, from: 170 * DAY, to: 200 * DAY }], meanFrom: 100 * DAY };
const sequences = plan.sequences.map((s) => ({ ...s, times: [], frames: [] }));
const t0 = performance.now(), steps = Math.round(plan.days * DAY / dt), mass0 = core.mass(pi);
for (let n = 0; n <= steps; n++) {
  const time = n * dt;
  for (const s of sequences) if (time >= s.from && time <= s.to + 1 && Math.abs(time / s.every - Math.round(time / s.every)) < 1e-6) { s.times.push(time); s.frames.push(frame()); }
  if (plan.meanFrom !== undefined && time >= plan.meanFrom && Math.round(time / dt) % Math.round(6 * 3600 / dt) === 0) accumulate();
  if (n % Math.round(DAY / dt) === 0) {
    let min = Infinity, max = -Infinity; for (const p of pi) { min = Math.min(min, p); max = Math.max(max, p); }
    console.error(`${which} day ${(time / DAY).toFixed(0)}: surface pressure ${(min / 100).toFixed(1)}–${(max / 100).toFixed(1)} hPa, mass drift ${((core.mass(pi) - mass0) / mass0).toExponential(1)}, ${((performance.now() - t0) / 1000).toFixed(0)} s`);
  }
  if (n < steps) step(tendency, state, dt);
}

const header = {
  version: 1, case: which, N, K, dt, levels: Array.from(core.levels), fields: [{ name: 'surfacePressure', unit: 'hPa', scale: 1 / 20, offset: 1000 }, { name: 'temperature850', unit: 'K', scale: 1 / 50, offset: 260 }, { name: 'wind250', unit: 'm/s', scale: 1 / 100, offset: 0 }],
  sequences: sequences.map((s) => ({ name: s.name, times: s.times })),
};
if (plan.meanFrom !== undefined) {
  header.zonalMean = { bins: BINS, from: plan.meanFrom, to: plan.days * DAY, samples: means.samples, sigma: Array.from(sigmaMid), u: Array.from(means.u, (v, n) => +(v / means.count[n]).toFixed(3)), T: Array.from(means.T, (v, n) => +(v / means.count[n]).toFixed(3)) };
}
const json = Buffer.from(JSON.stringify(header)), pad = (4 - ((8 + json.length) % 4)) % 4;
const blocks = [];
for (const s of sequences) for (const f of s.frames) for (const field of f) blocks.push(Buffer.from(field.buffer));
const head = Buffer.alloc(8); head.write('EXF1', 0); head.writeUInt32LE(json.length + pad, 4);
mkdirSync(new URL('../explain/data/', import.meta.url), { recursive: true });
writeFileSync(out, Buffer.concat([head, json, Buffer.alloc(pad, 32), ...blocks]));
console.error(`wrote ${out}: ${sequences.map((s) => `${s.name} ${s.frames.length} frames`).join(', ')}, ${((performance.now() - t0) / 60000).toFixed(1)} min`);
