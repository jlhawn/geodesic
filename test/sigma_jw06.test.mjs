import { test } from 'node:test';
import assert from 'node:assert/strict';
import { N, build, initialize, integrate, surfacePressureStats } from './helpers/jw06.mjs';

const STEADY_DAYS = +(process.env.SIGMA_STEADY_DAYS ?? 5);

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
