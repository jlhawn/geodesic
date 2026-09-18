import { Grid } from '../../js/grid.module.js';
import { buildMesh } from '../../js/mesh.module.js';
import { createShallowWater } from '../../js/dynamics/shallowWater.module.js';
import { createRK4 } from '../../js/dynamics/integrators.module.js';

export const EARTH = { a: 6.37122e6, omega: 7.292e-5, g: 9.80616 };
export const DAY = 86400;

export function setup(N, options = {}) {
  const grid = new Grid(N, { relax: options.relax ?? 10 });
  const mesh = buildMesh(grid, { radius: EARTH.a, omega: EARTH.omega });
  const modelOptions = typeof options.model === 'function' ? options.model(mesh) : options.model;
  const model = createShallowWater(mesh, { g: EARTH.g, ...modelOptions });
  const step = createRK4(mesh.nCells, mesh.nEdges);
  return { grid, mesh, model, step };
}

export function meanSpacing(mesh) {
  let sum = 0;
  for (let e = 0; e < mesh.nEdges; e++) sum += mesh.dcEdge[e];
  return sum / mesh.nEdges;
}

export function hyperdiffusion(mesh, hours) {
  return (meanSpacing(mesh) / Math.PI) ** 4 / (hours * 3600);
}

export function lonLat(x, y, z) {
  return { lon: Math.atan2(y, x), lat: Math.atan2(z, Math.hypot(x, y)) };
}

export function edgeNormalVelocity(mesh, windAt) {
  const { nEdges, xEdge, nEdge } = mesh;
  const u = new Float64Array(nEdges);
  for (let e = 0; e < nEdges; e++) {
    const x = xEdge[3 * e], y = xEdge[3 * e + 1], z = xEdge[3 * e + 2];
    const { lon, lat } = lonLat(x, y, z);
    const { zonal, meridional } = windAt(lon, lat);
    const cosLon = Math.cos(lon), sinLon = Math.sin(lon), cosLat = Math.cos(lat), sinLat = Math.sin(lat);
    const vx = -sinLon * zonal - sinLat * cosLon * meridional;
    const vy = cosLon * zonal - sinLat * sinLon * meridional;
    const vz = cosLat * meridional;
    u[e] = vx * nEdge[3 * e] + vy * nEdge[3 * e + 1] + vz * nEdge[3 * e + 2];
  }
  return u;
}

export function cellField(mesh, valueAt) {
  const { nCells, xCell } = mesh;
  const h = new Float64Array(nCells);
  for (let i = 0; i < nCells; i++) {
    const { lon, lat } = lonLat(xCell[3 * i], xCell[3 * i + 1], xCell[3 * i + 2]);
    h[i] = valueAt(lon, lat);
  }
  return h;
}

export function cellNorms(mesh, field, reference) {
  let num2 = 0, den2 = 0, numInf = 0, denInf = 0;
  for (let i = 0; i < mesh.nCells; i++) {
    const d = field[i] - reference[i];
    num2 += mesh.areaCell[i] * d * d;
    den2 += mesh.areaCell[i] * reference[i] * reference[i];
    numInf = Math.max(numInf, Math.abs(d));
    denInf = Math.max(denInf, Math.abs(reference[i]));
  }
  return { l2: Math.sqrt(num2 / den2), linf: numInf / denInf };
}

export function edgeNorms(mesh, field, reference) {
  let num2 = 0, den2 = 0, numInf = 0, denInf = 0;
  for (let e = 0; e < mesh.nEdges; e++) {
    const w = mesh.dcEdge[e] * mesh.dvEdge[e];
    const d = field[e] - reference[e];
    num2 += w * d * d;
    den2 += w * reference[e] * reference[e];
    numInf = Math.max(numInf, Math.abs(d));
    denInf = Math.max(denInf, Math.abs(reference[e]));
  }
  return { l2: Math.sqrt(num2 / den2), linf: numInf / denInf };
}

export function run(model, step, h, u, dt, seconds, onDay) {
  const steps = Math.round(seconds / dt);
  const perDay = Math.round(DAY / dt);
  for (let n = 1; n <= steps; n++) {
    step(model.tendency, h, u, dt);
    if (onDay && n % perDay === 0) onDay(n / perDay);
  }
}
