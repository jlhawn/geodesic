import { test } from 'node:test';
import assert from 'node:assert/strict';
import { N, build, initialize, integrate, surfacePressureStats } from './helpers/jw06.mjs';

const WAVE_DAYS = +(process.env.SIGMA_WAVE_DAYS ?? 10);

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
