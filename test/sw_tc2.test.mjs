import { test } from 'node:test';
import assert from 'node:assert/strict';
import { EARTH, DAY, setup, edgeNormalVelocity, cellField, cellNorms, edgeNorms, run } from './helpers/sphere.mjs';

const SIZES = (process.env.SW_TEST_N ?? '16,32').split(',').map(Number);
const DAYS = +(process.env.SW_TEST_DAYS ?? 5);
const { a, omega, g } = EARTH;
const u0 = 2 * Math.PI * a / (12 * DAY);
const h0 = 2.94e4 / g;

const errors = {};
for (const N of SIZES) {
  test(`Williamson TC2 at N=${N}: steady geostrophic zonal flow stays steady for ${DAYS} days`, () => {
    const { mesh, model, step } = setup(N);
    const wind = () => ({ zonal: 0, meridional: 0 });
    const uExact = edgeNormalVelocity(mesh, (lon, lat) => ({ zonal: u0 * Math.cos(lat), meridional: 0 }));
    const hExact = cellField(mesh, (lon, lat) => h0 - (a * omega * u0 + u0 * u0 / 2) * Math.sin(lat) ** 2 / g);
    const h = Float64Array.from(hExact);
    const u = Float64Array.from(uExact);
    const dt = 300 * 32 / N;
    const before = model.diagnostics(h, u);
    const t0 = performance.now();
    run(model, step, h, u, dt, DAYS * DAY);
    const seconds = (performance.now() - t0) / 1000;
    const after = model.diagnostics(h, u);
    const eh = cellNorms(mesh, h, hExact);
    const eu = edgeNorms(mesh, u, uExact);
    errors[N] = eh.l2;
    console.log(`TC2 N=${N} dt=${dt}s: h l2 ${eh.l2.toExponential(2)} linf ${eh.linf.toExponential(2)}; u l2 ${eu.l2.toExponential(2)} linf ${eu.linf.toExponential(2)}; mass drift ${((after.mass - before.mass) / before.mass).toExponential(1)}; energy drift ${((after.energy - before.energy) / before.energy).toExponential(1)}; ${seconds.toFixed(1)} s`);
    assert.ok(Number.isFinite(eh.l2) && Number.isFinite(eu.l2));
    assert.ok(Math.abs(after.mass - before.mass) / before.mass < 1e-12);
    if (N >= 32) assert.ok(eh.l2 < 1e-3);
  });
}

test('TC2 error decreases with resolution', () => {
  const Ns = SIZES.slice().sort((p, q) => p - q);
  for (let k = 1; k < Ns.length; k++) {
    assert.ok(errors[Ns[k]] < errors[Ns[k - 1]], `N=${Ns[k - 1]}: ${errors[Ns[k - 1]]} vs N=${Ns[k]}: ${errors[Ns[k]]}`);
  }
});
