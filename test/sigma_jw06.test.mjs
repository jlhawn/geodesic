import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createSigmaCore, sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { createRK4Arrays } from '../js/dynamics/integrators.module.js';
import { EARTH, DAY, edgeNormalVelocity, hyperdiffusion } from './helpers/sphere.mjs';

const N = +(process.env.SIGMA_TEST_N ?? 16);
const STEADY_DAYS = +(process.env.SIGMA_STEADY_DAYS ?? 5);
const WAVE_DAYS = +(process.env.SIGMA_WAVE_DAYS ?? 10);
const { a, omega, g } = EARTH;
const R = 287, cp = 1004.5, p0 = 1e5;
const eta0 = 0.252, u0 = 35, T0 = 288, gamma = 0.005, deltaT = 4.8e5, etaT = 0.2;

function etaV(eta) { return (eta - eta0) * Math.PI / 2; }
function shape(lat) {
  const s = Math.sin(lat), c = Math.cos(lat);
  return { A: -2 * s ** 6 * (c * c + 1 / 3) + 10 / 63, B: 8 / 5 * c ** 3 * (s * s + 2 / 3) - Math.PI / 4 };
}
function zonalWind(eta, lat) {
  return u0 * Math.cos(etaV(eta)) ** 1.5 * Math.sin(2 * lat) ** 2;
}
function geopotentialJW(eta, lat) {
  const exp = R * gamma / g;
  let mean = T0 * g / gamma * (1 - eta ** exp);
  if (eta < etaT) {
    mean -= R * deltaT * ((Math.log(eta / etaT) + 137 / 60) * etaT ** 5 - 5 * etaT ** 4 * eta + 5 * etaT ** 3 * eta ** 2 - (10 / 3) * etaT ** 2 * eta ** 3 + (5 / 4) * etaT * eta ** 4 - eta ** 5 / 5);
  }
  const cv = Math.cos(etaV(eta)) ** 1.5;
  const { A, B } = shape(lat);
  return mean + u0 * cv * (A * u0 * cv + B * a * omega);
}
function temperatureJW(eta, lat) {
  const exp = R * gamma / g;
  let mean = T0 * eta ** exp;
  if (eta < etaT) mean += deltaT * (etaT - eta) ** 5;
  const ev = etaV(eta);
  const { A, B } = shape(lat);
  return mean + 0.75 * (eta * Math.PI * u0 / R) * Math.sin(ev) * Math.sqrt(Math.cos(ev)) * (A * 2 * u0 * Math.cos(ev) ** 1.5 + B * a * omega);
}
function perturbation(lon, lat) {
  const lonC = Math.PI / 9, latC = 2 * Math.PI / 9, radius = a / 10;
  const r = a * Math.acos(Math.min(1, Math.sin(latC) * Math.sin(lat) + Math.cos(latC) * Math.cos(lat) * Math.cos(lon - lonC)));
  return Math.exp(-((r / radius) ** 2));
}

function initialize(mesh, core, perturbed) {
  const C = mesh.nCells, E = mesh.nEdges, K = core.K;
  const pi = new Float64Array(C).fill(p0);
  const surface = core.surface;
  const theta = new Float64Array(K * C);
  core.diagnose(pi, theta);
  const { exnerLayer, exnerLower } = core.arrays;
  for (let i = 0; i < C; i++) {
    const lat = mesh.latCell[i];
    let below = surface[i], thetaBelow = 0;
    for (let k = K - 1; k >= 0; k--) {
      const idx = k * C + i;
      const target = geopotentialJW(core.sigmaMid[k], lat);
      const fromBelow = k === K - 1 ? 0 : cp * thetaBelow * (exnerLayer[idx + C] - exnerLower[idx]);
      theta[idx] = (target - below - fromBelow) / (cp * (exnerLower[idx] - exnerLayer[idx]));
      below = target;
      thetaBelow = theta[idx];
    }
  }
  const u = new Float64Array(K * E);
  for (let k = 0; k < K; k++) {
    const eta = core.sigmaMid[k];
    const uk = edgeNormalVelocity(mesh, (lon, lat) => ({ zonal: zonalWind(eta, lat) + (perturbed ? perturbation(lon, lat) : 0), meridional: 0 }));
    u.set(uk, k * E);
  }
  return [pi, theta, u];
}

function build(N) {
  const grid = new Grid(N);
  const mesh = buildMesh(grid, { radius: a, omega });
  const surface = Float64Array.from(mesh.latCell, (lat) => geopotentialJW(1, lat));
  const core = createSigmaCore(mesh, { g, cp, R, p0, surfaceGeopotential: surface, nu4: hyperdiffusion(mesh, 3), nu4Theta: hyperdiffusion(mesh, 3) });
  core.surface = surface;
  const step = createRK4Arrays([mesh.nCells, core.K * mesh.nCells, core.K * mesh.nEdges]);
  return { mesh, core, step };
}

function surfacePressureStats(pi) {
  let min = Infinity, max = -Infinity;
  for (const p of pi) { min = Math.min(min, p); max = Math.max(max, p); }
  return { min: min / 100, max: max / 100 };
}

function integrate(core, step, state, dt, days, onDay) {
  const perDay = Math.round(DAY / dt);
  for (let d = 1; d <= days; d++) {
    for (let n = 0; n < perDay; n++) step(core.tendency, state, dt);
    onDay(d);
  }
}

test(`JW06 at N=${N}: the balanced baroclinic base state stays steady for ${STEADY_DAYS} days`, () => {
  const { mesh, core, step } = build(N);
  const state = initialize(mesh, core, false);
  const u0 = Float64Array.from(state[2]);
  const dt = 450 * 16 / N;
  const mass0 = core.mass(state[0]);
  const rows = [];
  const t0 = performance.now();
  integrate(core, step, state, dt, STEADY_DAYS, (d) => {
    const ps = surfacePressureStats(state[0]);
    let du = 0, scale = 0;
    for (let x = 0; x < u0.length; x++) { du = Math.max(du, Math.abs(state[2][x] - u0[x])); scale = Math.max(scale, Math.abs(u0[x])); }
    rows.push(`d${d} ps ${ps.min.toFixed(2)}–${ps.max.toFixed(2)} hPa, max|Δu| ${du.toFixed(2)} m/s`);
  });
  const seconds = (performance.now() - t0) / 1000;
  const ps = surfacePressureStats(state[0]);
  console.log(`JW06 steady N=${N} K=${core.K} dt=${dt}s: ${rows.join(' | ')}; mass drift ${((core.mass(state[0]) - mass0) / mass0).toExponential(1)}; ${seconds.toFixed(0)} s`);
  assert.ok(Number.isFinite(ps.min));
  assert.ok(Math.abs(core.mass(state[0]) - mass0) / mass0 < 1e-12);
  assert.ok(1000 - ps.min < 2 && ps.max - 1000 < 2);
  let du = 0;
  for (let x = 0; x < u0.length; x++) du = Math.max(du, Math.abs(state[2][x] - u0[x]));
  assert.ok(du < 2);
});

test(`JW06 at N=${N}: the perturbed base state grows a baroclinic wave within ${WAVE_DAYS} days`, () => {
  const { mesh, core, step } = build(N);
  const state = initialize(mesh, core, true);
  const dt = 450 * 16 / N;
  const rows = [];
  const t0 = performance.now();
  integrate(core, step, state, dt, WAVE_DAYS, (d) => {
    const ps = surfacePressureStats(state[0]);
    rows.push(`d${d} ps ${ps.min.toFixed(1)}–${ps.max.toFixed(1)}`);
  });
  const seconds = (performance.now() - t0) / 1000;
  const ps = surfacePressureStats(state[0]);
  console.log(`JW06 wave N=${N}: ${rows.join(' | ')}; ${seconds.toFixed(0)} s`);
  assert.ok(Number.isFinite(ps.min));
  assert.ok(ps.min < 990);
});
