import { test } from 'node:test';
import assert from 'node:assert/strict';
import { ekmanSpiral, stressOfWind, windForStress, tracks, coriolis, inertialPeriod, WIND_HEIGHT, RHO } from '../explain/ekmanTheory.module.js';
import { ekmanCase, runCase } from '../explain/oceanCases.module.js';

const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const degrees = ([east, north]) => Math.atan2(north, east) * 180 / Math.PI;

test('the spiral carries τ/(ρ₀f) at right angles to the wind, and turns 45° at the surface', () => {
  for (const [latitude, sign] of [[45, -1], [-30, 1]]) {
    const tau = 0.1, f = coriolis(latitude), spiral = ekmanSpiral({ tau, f, viscosity: 0.02 });
    const steps = 80000, dz = 40 * spiral.depth / steps;
    let east = 0, north = 0;
    for (let s = 0; s < steps; s++) { const [u, v] = spiral.at(-(s + 0.5) * dz); east += u * dz; north += v * dz; }
    close(east, 0, 1e-6, 'eastward transport');
    close(north, sign * tau / (RHO * Math.abs(f)), 1e-6, 'northward transport');
    assert.deepEqual(spiral.transport, [0, -tau / (RHO * f)]);
    close(degrees(spiral.at(0)), sign * 45, 1e-9, 'surface angle');
    close(Math.hypot(...spiral.at(0)), tau / (RHO * Math.sqrt(0.02 * Math.abs(f))), 1e-12, 'surface speed');
    const deep = spiral.at(-spiral.reversal), top = spiral.at(0);
    close((deep[0] * top[0] + deep[1] * top[1]) / (Math.hypot(...deep) * Math.hypot(...top)), -1, 1e-12, 'the current opposes the surface at πd');
    close(Math.hypot(...deep) / Math.hypot(...top), Math.exp(-Math.PI), 1e-12, 'and is e^−π as fast');
  }
});

test('the spiral solves A u″ = −f ẑ × u with the wind stress at the top', () => {
  const tau = 0.15, f = coriolis(50), A = 0.04, spiral = ekmanSpiral({ tau, f, viscosity: A }), h = 1e-3;
  for (const z of [-3, -20, -55]) {
    const [u0, v0] = spiral.at(z), [up, vp] = spiral.at(z + h), [um, vm] = spiral.at(z - h);
    close(A * (up - 2 * u0 + um) / (h * h), -f * v0, 1e-9, `east at ${z} m`);
    close(A * (vp - 2 * v0 + vm) / (h * h), f * u0, 1e-9, `north at ${z} m`);
  }
  const [u1, v1] = spiral.at(-h), [u2, v2] = spiral.at(-2 * h), [u0, v0] = spiral.at(0);
  close(RHO * A * (3 * u0 - 4 * u1 + u2) / (2 * h), tau, 1e-6, 'stress east at the surface');
  close(RHO * A * (3 * v0 - 4 * v1 + v2) / (2 * h), 0, 1e-6, 'no stress north at the surface');
});

test('the wind stress follows the model’s own drag over the sea', () => {
  close(WIND_HEIGHT, 20, 0.5, 'the lowest bl36 level');
  for (const [speed, tau] of [[5, 0.0247], [10, 0.1348], [15, 0.4043]]) close(stressOfWind(speed), tau, 5e-4, `${speed} m/s`);
  for (const tau of [0.02, 0.1, 0.3]) close(stressOfWind(windForStress(tau)), tau, 1e-9, `the wind for ${tau} N/m²`);
});

test('the tracks integrate the series, and the period is the one the case uses', () => {
  const n = 101, time = Float64Array.from({ length: n }, (_, s) => 100 * s), w = 2 * Math.PI / 10000;
  const circling = { time, east: time.map((t) => Math.cos(w * t)), north: time.map((t) => Math.sin(w * t)), relativeEast: time.map(() => -0.5), relativeNorth: time.map(() => 0) };
  const { slab, beneath } = tracks(circling);
  close(Math.hypot(slab[2 * n - 2], slab[2 * n - 1]), 0, 2e-2, 'a whole inertial circle closes');
  close(beneath[2 * n - 2], slab[2 * n - 2] + 0.5 * 10000, 1e-9, 'the water beneath runs half a meter a second faster east');

  const sim = ekmanCase({ latitude: 45, days: 1 }), frame = runCase(sim);
  close(frame.inertialPeriod, inertialPeriod(45), 1e-9, 'inertial period');
  close(inertialPeriod(45) / 3600, 16.92, 0.005, 'inertial period in hours');
  const path = tracks(frame.series), last = frame.series.time.length - 1;
  let east = 0;
  for (let s = 1; s <= last; s++) east += 0.5 * (frame.series.time[s] - frame.series.time[s - 1]) * (frame.series.east[s] + frame.series.east[s - 1]);
  close(path.slab[2 * last], east, 1e-6, 'the slab’s eastward run');
  assert.ok(path.slab[2 * last + 1] < -300, `the slab drifts south in the north: ${path.slab[2 * last + 1]} m`);
});
