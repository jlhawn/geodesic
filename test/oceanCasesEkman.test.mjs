import { test } from 'node:test';
import assert from 'node:assert/strict';
import { caseMesh, ekman, EKMAN_REACH, RHO, DEG } from '../explain/oceanCases.module.js';
import { ekmanSpiral, fitViscosity, columnDays } from '../explain/ekmanTheory.module.js';
import { HYBRID_LAYERS } from '../js/ocean/column.module.js';

const N = 16;
const mesh = caseMesh(N);
const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const between = (v, lo, hi, what) => assert.ok(v >= lo && v <= hi, `${what}: ${v} not within ${lo} to ${hi}`);
const angles = (profile, upTo) => Array.from({ length: upTo }, (_, k) => Math.atan2(-profile.north[k], profile.east[k]) / DEG);

test('a band of westerlies at 45°N drives a spiral of currents whose water carries τ/(ρ₀f) at right angles to the wind, to the right, and to the left at 45°S', () => {
  const north = ekman({ N, latitude: 45, days: columnDays(45) }), south = ekman({ N, latitude: -45, days: columnDays(-45) });
  const show = (r) => `${r.transport.magnitude.toFixed(4)} m²/s at ${r.transport.angle.toFixed(1)}° (${r.transport.ratio.toFixed(3)} of τ/(ρ₀f) = ${r.theory.toFixed(4)}); top layer ${(100 * r.top.magnitude).toFixed(2)} cm/s at ${r.top.angle.toFixed(1)}°; mixed layer ${(100 * r.mixed.magnitude).toFixed(2)} cm/s at ${r.mixed.angle.toFixed(1)}°, ${r.mixedDepth.toFixed(1)} m deep, boundary layer ${r.boundaryDepth.toFixed(1)} m, summed over ${r.reach.toFixed(1)} m; mean viscosity ${r.meanViscosity.toFixed(4)} m²/s`;
  console.log(`45°N: ${show(north)}; inertial period ${(north.inertialPeriod / 3600).toFixed(2)} h, ${north.ringCells} cells`);
  console.log(`45°S: ${show(south)}`);
  close(north.inertialPeriod / 3600, 16.92, 0.01, 'inertial period at 45° in hours');
  close(north.theory, 0.1 / (RHO * 2 * mesh.omega * Math.sin(45 * DEG)), 0.025, 'the ring’s τ/(ρ₀f) is near the one at 45°');
  assert.ok(north.settled && north.averagedOver >= 2 * north.inertialPeriod - 1350 && north.averagedOver <= 2 * north.inertialPeriod + 1350);
  close(north.transport.ratio, 1, 0.05, 'Ekman transport against τ/(ρ₀f) at 45°N');
  between(north.transport.angle, 86, 94, 'the transport’s angle to the right of the wind');
  close(south.transport.ratio, north.transport.ratio, 0.01, 'the same transport at 45°S');
  close(south.transport.angle, -north.transport.angle, 0.5, 'to the left of the wind at 45°S');
  between(north.top.angle, 45, 60, 'the top layer drifts about 45° to the right'); between(100 * north.top.magnitude, 5, 6.2, 'top layer cm/s');
  between(north.mixed.angle, 80, 92, 'the mixed layer’s drift is at right angles'); between(100 * north.mixed.magnitude, 2.3, 2.9, 'mixed layer cm/s');
  close(north.mixed.magnitude * north.mixedDepth, north.transport.magnitude, 0.03 * north.transport.magnitude, 'the mixed layer carries the transport');
  between(north.mixedDepth, 32, 38, 'mixed layer depth'); between(north.boundaryDepth, 35, 41, 'boundary layer depth');
  assert.ok(north.mixedDepth < north.boundaryDepth && north.boundaryDepth * EKMAN_REACH === north.reach);
  between(north.meanViscosity, 0.011, 0.017, 'viscosity averaged over the boundary layer');

  const p = north.profile, a = angles(p, 12);
  assert.equal(p.center.length, HYBRID_LAYERS.K);
  between(a[0], 45, 60, 'the top layer’s angle'); assert.ok(a[8] > a[0] + 90, `the current turns more than 90° to the right between the top layer and ${p.center[8].toFixed(0)} m: ${a[0]} to ${a[8]}`);
  for (let k = 1; k <= 7; k++) assert.ok(a[k] > a[k - 1], `the current turns further right at layer ${k}`);
  assert.ok(Math.hypot(p.east[9], p.north[9]) < 0.05 * Math.hypot(p.east[0], p.north[0]), 'the current has died away by 40 m');
  const nu = Array.from(p.viscosity);
  between(Math.max(...nu.slice(0, 9)), 0.02, 0.03, 'the largest viscosity within the boundary layer'); assert.ok(nu[12] < 5e-4, `viscosity below the boundary layer: ${nu[12]}`);

  const n = north.series.time.length;
  assert.equal(n, columnDays(45) * 64 + 1);
  for (const name of Object.keys(north.series)) assert.equal(north.series[name].length, n, name);
  let fastest = 0;
  for (let s = 0; s < n; s++) fastest = Math.max(fastest, Math.hypot(north.series.mixedEast[s] - north.series.deepEast[s], north.series.mixedNorth[s] - north.series.deepNorth[s]));
  const deepFit = (r) => { const q = fitViscosity({ tau: r.stressMean, f: r.fMean, center: r.profile.center, thickness: r.profile.thickness, east: r.profile.east, north: r.profile.north, reach: r.reach }); return { ...q, surface: Math.abs(ekmanSpiral({ tau: r.stressMean, f: r.fMean, viscosity: q.viscosity }).surface) }; };
  const fit45 = deepFit(north);
  console.log(`transport ${north.transport.magnitude.toFixed(5)} m²/s against ${north.theory.toFixed(5)}; fitted viscosity ${fit45.viscosity.toFixed(4)} m²/s, misfit ${(100 * fit45.rms).toFixed(2)} cm/s of ${(100 * fit45.surface).toFixed(2)}, ${(100 * fit45.rms / fit45.surface).toFixed(1)} percent`);
  between(fit45.rms / fit45.surface, 0.04, 0.09, 'the spiral fits the 45°N layers to within a few percent of its surface current');
  let loops = 0;
  for (let s = 0; s < north.series.time.length; s++) if (north.series.time[s] > north.series.time[north.series.time.length - 1] - 17 * 3600) loops = Math.max(loops, Math.hypot(north.series.mixedEast[s] - north.series.deepEast[s] - north.mixed.east, north.series.mixedNorth[s] - north.series.deepNorth[s] - north.mixed.north));
  console.log(`on day 4 the mixed layer’s loops reach ${(100 * loops).toFixed(2)} cm/s around its drift of ${(100 * north.mixed.magnitude).toFixed(2)} cm/s`);
  between(100 * loops, 0.33, 0.57, 'the inertial loops on day 4, cm/s');
  console.log(`the mixed layer’s speed against the water below peaks at ${(100 * fastest).toFixed(2)} cm/s while it circles and drifts at ${(100 * north.mixed.magnitude).toFixed(2)} cm/s`);
  assert.ok(fastest > 2 * north.mixed.magnitude, 'the mixed layer overshoots in an inertial circle before it settles');
  console.log(`the mixed layer deepens from ${north.series.mixedDepth[0].toFixed(1)} m to ${north.series.mixedDepth[n - 1].toFixed(1)} m and the boundary layer from ${north.series.boundaryDepth[0].toFixed(1)} to ${north.series.boundaryDepth[n - 1].toFixed(1)} m; layer angles ${a.slice(0, 9).map((v) => v.toFixed(0)).join(', ')}° at ${Array.from(p.center.slice(0, 9), (z) => z.toFixed(0)).join(', ')} m; speeds ${Array.from({ length: 9 }, (_, k) => (100 * Math.hypot(p.east[k], p.north[k])).toFixed(2)).join(', ')} cm/s; viscosity ${nu.slice(0, 10).map((v) => v.toFixed(4)).join(', ')} m²/s`);
  between(north.series.mixedDepth[0], 4.5, 6, 'the analytic start’s mixed layer'); between(north.series.mixedDepth[n - 1], 34, 38, 'the mixed layer after four days');
  between(north.series.boundaryDepth[0], 0, 0.01, 'no boundary layer before the first step'); between(north.series.boundaryDepth[n - 1], 37, 41, 'the boundary layer after four days');
  between(a[4], 88, 100, 'the current at 16 m points at right angles to the wind'); between(Math.hypot(p.east[7], p.north[7]) / Math.hypot(p.east[0], p.north[0]), 0.25, 0.42, 'a third of the top layer’s speed remains at 29 m');
});
