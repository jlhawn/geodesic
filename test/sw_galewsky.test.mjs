import { test } from 'node:test';
import assert from 'node:assert/strict';
import { curl } from '../js/dynamics/operators.module.js';
import { EARTH, DAY, setup, edgeNormalVelocity, cellField, cellNorms, run, hyperdiffusion } from './helpers/sphere.mjs';

const N = +(process.env.SW_TEST_N ?? 32);
const { a, omega, g } = EARTH;
const phi0 = Math.PI / 7, phi1 = Math.PI / 2 - phi0, uMax = 80, en = Math.exp(-4 / (phi1 - phi0) ** 2);
const hMean = 1e4;

function jet(lat) {
  if (lat <= phi0 || lat >= phi1) return 0;
  return (uMax / en) * Math.exp(1 / ((lat - phi0) * (lat - phi1)));
}

function balancedHeight() {
  const steps = 20000;
  const dphi = Math.PI / steps;
  const table = new Float64Array(steps + 1);
  let integral = 0;
  for (let k = 1; k <= steps; k++) {
    const mid = -Math.PI / 2 + (k - 0.5) * dphi;
    const u = jet(mid);
    integral += a * u * (2 * omega * Math.sin(mid) + Math.tan(mid) * u / a) * dphi;
    table[k] = -integral / g;
  }
  let weighted = 0, area = 0;
  for (let k = 0; k <= steps; k++) {
    const w = Math.cos(-Math.PI / 2 + k * dphi);
    weighted += w * table[k];
    area += w;
  }
  const offset = hMean - weighted / area;
  return (lat) => {
    const x = (lat + Math.PI / 2) / dphi;
    const k = Math.min(steps - 1, Math.max(0, Math.floor(x)));
    const f = x - k;
    return offset + table[k] * (1 - f) + table[k + 1] * f;
  };
}

function perturbation(lon, lat) {
  const alpha = 1 / 3, beta = 1 / 15, phi2 = Math.PI / 4, amplitude = 120;
  return amplitude * Math.cos(lat) * Math.exp(-((lon / alpha) ** 2)) * Math.exp(-(((phi2 - lat) / beta) ** 2));
}

function maxVorticity(mesh, u) {
  const zeta = curl(mesh, u);
  let m = 0;
  for (let v = 0; v < mesh.nVertices; v++) m = Math.max(m, Math.abs(zeta[v]));
  return m;
}

function eddyEnergy(mesh, model, h, u, uZonal) {
  let e = 0;
  for (let k = 0; k < mesh.nEdges; k++) {
    const he = 0.5 * (h[mesh.cellsOnEdge[2 * k]] + h[mesh.cellsOnEdge[2 * k + 1]]);
    e += 0.5 * mesh.dcEdge[k] * mesh.dvEdge[k] * he * (u[k] - uZonal[k]) ** 2;
  }
  return e;
}

const closure = (mesh) => ({ nu4: hyperdiffusion(mesh, 3) });

test(`Galewsky jet at N=${N}: the balanced zonal jet holds for 3 days, before grid noise seeds the instability`, () => {
  const { mesh, model, step } = setup(N, { model: closure });
  const height = balancedHeight();
  const uExact = edgeNormalVelocity(mesh, (lon, lat) => ({ zonal: jet(lat), meridional: 0 }));
  const hExact = cellField(mesh, (lon, lat) => height(lat));
  const h = Float64Array.from(hExact), u = Float64Array.from(uExact);
  const dt = 240 * 32 / N;
  const trace = [];
  run(model, step, h, u, dt, 3 * DAY, (day) => trace.push(cellNorms(mesh, h, hExact).l2));
  const eh = cellNorms(mesh, h, hExact);
  const zeta0 = maxVorticity(mesh, uExact), zeta = maxVorticity(mesh, u);
  console.log(`Galewsky unperturbed N=${N}: h l2 by day ${trace.map((x) => x.toExponential(2)).join(', ')}; linf ${eh.linf.toExponential(2)}; max |ζ| ×${(zeta / zeta0).toFixed(2)}`);
  if (N >= 32) assert.ok(eh.l2 < 1e-3);
  assert.ok(zeta > 0.8 * zeta0 && zeta < 1.2 * zeta0);
});

test(`Galewsky jet at N=${N}: the perturbed jet goes barotropically unstable within 6 days`, () => {
  const { mesh, model, step } = setup(N, { model: closure });
  const height = balancedHeight();
  const uZonal = edgeNormalVelocity(mesh, (lon, lat) => ({ zonal: jet(lat), meridional: 0 }));
  const h = cellField(mesh, (lon, lat) => height(lat) + perturbation(lon, lat));
  const u = Float64Array.from(uZonal);
  const dt = 240 * 32 / N;
  const zeta0 = maxVorticity(mesh, u);
  const start = model.diagnostics(h, u);
  const growth = [];
  const t0 = performance.now();
  run(model, step, h, u, dt, 6 * DAY, (day) => {
    growth.push({ day, eke: eddyEnergy(mesh, model, h, u, uZonal) / start.kinetic, zeta: maxVorticity(mesh, u) / zeta0 });
  });
  const seconds = (performance.now() - t0) / 1000;
  const end = model.diagnostics(h, u);
  console.log(`Galewsky perturbed N=${N}: ` + growth.map((r) => `day ${r.day}: EKE/KE ${r.eke.toExponential(2)}, max|ζ| ×${r.zeta.toFixed(2)}`).join('; ') + `; mass drift ${((end.mass - start.mass) / start.mass).toExponential(1)}; energy drift ${((end.energy - start.energy) / start.energy).toExponential(2)}; ${seconds.toFixed(0)} s`);
  assert.ok(Number.isFinite(end.energy));
  assert.ok(growth[5].eke > 10 * growth[0].eke);
  assert.ok(growth[5].eke > 1e-3);
  assert.ok(Math.abs(end.mass - start.mass) / start.mass < 1e-12);
});
