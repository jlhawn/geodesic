import { Grid } from '../../js/grid.module.js';
import { buildMesh } from '../../js/mesh.module.js';
import { createSigmaCore } from '../../js/dynamics/sigmaCore.module.js';
import { createRK4Arrays } from '../../js/dynamics/integrators.module.js';
import { EARTH, DAY, edgeNormalVelocity, hyperdiffusion } from './sphere.mjs';

export const N = +(process.env.SIGMA_TEST_N ?? 16);
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
function perturbation(lon, lat) {
  const lonC = Math.PI / 9, latC = 2 * Math.PI / 9, radius = a / 10;
  const r = a * Math.acos(Math.min(1, Math.sin(latC) * Math.sin(lat) + Math.cos(latC) * Math.cos(lat) * Math.cos(lon - lonC)));
  return Math.exp(-((r / radius) ** 2));
}

export function initialize(mesh, core, perturbed) {
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

export function build(N) {
  const grid = new Grid(N);
  const mesh = buildMesh(grid, { radius: a, omega });
  const surface = Float64Array.from(mesh.latCell, (lat) => geopotentialJW(1, lat));
  const core = createSigmaCore(mesh, { g, cp, R, p0, surfaceGeopotential: surface, nu4: hyperdiffusion(mesh, 3), nu4Theta: hyperdiffusion(mesh, 3) });
  core.surface = surface;
  const step = createRK4Arrays([mesh.nCells, core.K * mesh.nCells, core.K * mesh.nEdges]);
  return { mesh, core, step };
}

export function surfacePressureStats(pi) {
  let min = Infinity, max = -Infinity;
  for (const p of pi) { min = Math.min(min, p); max = Math.max(max, p); }
  return { min: min / 100, max: max / 100 };
}

export function integrate(core, step, state, dt, days, onDay) {
  const perDay = Math.round(DAY / dt);
  for (let d = 1; d <= days; d++) {
    for (let n = 0; n < perDay; n++) step(core.tendency, state, dt);
    onDay(d);
  }
}
