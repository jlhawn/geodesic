// Precomputes the layered-model frames that explain/sphere.html plays back:
//   node scripts/explainLayers.mjs jw06|heldsuarez [OUT]
// writing explain/data/<case>_N16.bin.
import { writeFileSync, mkdirSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createSigmaCore } from '../js/dynamics/sigmaCore.module.js';
import { createRK4Arrays } from '../js/dynamics/integrators.module.js';
import { EARTH, edgeNormalVelocity, hyperdiffusion, cellVelocity } from '../explain/swCases.module.js';

const which = process.argv[2], N = 16, DAY = 86400, dt = 450 * 16 / N;
const out = process.argv[3] ?? new URL(`../explain/data/${which}_N${N}.bin`, import.meta.url).pathname;
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

const run = which === 'jw06' ? jw06() : which === 'heldsuarez' ? heldSuarez() : null;
if (!run) { console.error('usage: node scripts/explainLayers.mjs jw06|heldsuarez [OUT]'); process.exit(1); }
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
