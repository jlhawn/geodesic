import { test } from 'node:test';
import assert from 'node:assert/strict';
import { EARTH, DAY, setup, edgeNormalVelocity, cellField, lonLat, run } from './helpers/sphere.mjs';

const N = +(process.env.SW_TEST_N ?? 32);
const DAYS = +(process.env.SW_TEST_DAYS ?? 14);
const { a, omega, g } = EARTH;
const w = 7.848e-6, K = 7.848e-6, R = 4, h0 = 8e3;

function wind(lon, lat) {
  const c = Math.cos(lat), s = Math.sin(lat);
  return {
    zonal: a * w * c + a * K * c ** (R - 1) * (R * s * s - c * c) * Math.cos(R * lon),
    meridional: -a * K * R * c ** (R - 1) * s * Math.sin(R * lon),
  };
}

function height(lon, lat) {
  const c = Math.cos(lat);
  const A = (w / 2) * (2 * omega + w) * c * c + (K * K / 4) * c ** (2 * R) * ((R + 1) * c * c + (2 * R * R - R - 2) - 2 * R * R / (c * c));
  const B = (2 * (omega + w) * K / ((R + 1) * (R + 2))) * c ** R * ((R * R + 2 * R + 2) - (R + 1) ** 2 * c * c);
  const C = (K * K / 4) * c ** (2 * R) * ((R + 1) * c * c - (R + 2));
  return h0 + (a * a / g) * (A + B * Math.cos(R * lon) + C * Math.cos(2 * R * lon));
}

function wavePhase(mesh, h) {
  let cs = 0, sn = 0;
  for (let i = 0; i < mesh.nCells; i++) {
    const { lon, lat } = lonLat(mesh.xCell[3 * i], mesh.xCell[3 * i + 1], mesh.xCell[3 * i + 2]);
    if (Math.abs(lat) > 0.6) continue;
    const weight = mesh.areaCell[i] * Math.cos(lat) ** R;
    cs += weight * h[i] * Math.cos(R * lon);
    sn += weight * h[i] * Math.sin(R * lon);
  }
  return Math.atan2(sn, cs) / R;
}

test(`Williamson TC6 at N=${N}: Rossby–Haurwitz wave 4 propagates for ${DAYS} days conserving mass and energy`, () => {
  const { mesh, model, step } = setup(N);
  const u = edgeNormalVelocity(mesh, wind);
  const h = cellField(mesh, height);
  const dt = 300 * 32 / N;
  const start = model.diagnostics(h, u);
  const phase0 = wavePhase(mesh, h);
  const analyticRate = (R * (R + 3) * w - 2 * omega) / ((R + 1) * (R + 2));
  let unwrapped = 0, previous = phase0;
  const t0 = performance.now();
  run(model, step, h, u, dt, DAYS * DAY, (day) => {
    const phase = wavePhase(mesh, h);
    let delta = phase - previous;
    while (delta > Math.PI / R) delta -= 2 * Math.PI / R;
    while (delta < -Math.PI / R) delta += 2 * Math.PI / R;
    unwrapped += delta;
    previous = phase;
    if (day % 7 === 0) {
      const d = model.diagnostics(h, u);
      console.log(`  day ${day}: phase ${(unwrapped * 180 / Math.PI).toFixed(1)}° (analytic ${(analyticRate * day * DAY * 180 / Math.PI).toFixed(1)}°), energy drift ${((d.energy - start.energy) / start.energy).toExponential(2)}, enstrophy drift ${((d.enstrophy - start.enstrophy) / start.enstrophy).toExponential(2)}`);
    }
  });
  const seconds = (performance.now() - t0) / 1000;
  const end = model.diagnostics(h, u);
  const measuredRate = unwrapped / (DAYS * DAY);
  let hMin = Infinity, hMax = -Infinity;
  for (let i = 0; i < mesh.nCells; i++) { hMin = Math.min(hMin, h[i]); hMax = Math.max(hMax, h[i]); }
  console.log(`TC6 N=${N}: phase speed ${(measuredRate * DAY * 180 / Math.PI).toFixed(2)}°/day vs analytic ${(analyticRate * DAY * 180 / Math.PI).toFixed(2)}°/day; h range ${hMin.toFixed(0)}–${hMax.toFixed(0)} m; mass drift ${((end.mass - start.mass) / start.mass).toExponential(1)}; energy drift ${((end.energy - start.energy) / start.energy).toExponential(2)}; enstrophy drift ${((end.enstrophy - start.enstrophy) / start.enstrophy).toExponential(2)}; ${seconds.toFixed(0)} s`);
  assert.ok(Number.isFinite(end.energy));
  assert.ok(Math.abs(end.mass - start.mass) / start.mass < 1e-12);
  assert.ok(Math.abs(end.energy - start.energy) / start.energy < 1e-3);
  assert.ok(Math.abs(end.enstrophy - start.enstrophy) / start.enstrophy < 1e-2);
  assert.ok(Math.abs(measuredRate - analyticRate) / analyticRate < 0.15);
  assert.ok(hMin > 7000 && hMax < 11000);
});
