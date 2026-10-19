import { test } from 'node:test';
import assert from 'node:assert/strict';
import { ekmanSpiral, fitViscosity, stressOfWind, windForStress, tracks, coriolis, inertialPeriod, WIND_HEIGHT, RHO } from '../explain/ekmanTheory.module.js';
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

test('the fit recovers the viscosity of a spiral sampled on the model’s layers, and of the spiral with a current added to it far below', () => {
  const thickness = Float64Array.from({ length: 30 }, (_, k) => (k < 12 ? 3 * 1.08 ** k : 10)), center = new Float64Array(30);
  for (let k = 0, z = 0; k < 30; k++) { center[k] = z + 0.5 * thickness[k]; z += thickness[k]; }
  for (const [latitude, viscosity] of [[45, 0.018], [-30, 0.03], [15, 0.01]]) {
    const tau = 0.1, f = coriolis(latitude), spiral = ekmanSpiral({ tau, f, viscosity });
    const east = Float64Array.from(center, (z) => spiral.at(-z)[0]), north = Float64Array.from(center, (z) => spiral.at(-z)[1]);
    const fit = fitViscosity({ tau, f, center, thickness, east, north, reach: 200 });
    close(fit.viscosity, viscosity, 1e-3 * viscosity, `viscosity at ${latitude}°`);
    close(fit.rms, 0, 1e-6, 'misfit');
  }
  const f = coriolis(45), spiral = ekmanSpiral({ tau: 0.1, f, viscosity: 0.02 });
  const east = Float64Array.from(center, (z, k) => spiral.at(-z)[0] * (1 + 0.1 * Math.sin(k))), north = Float64Array.from(center, (z, k) => spiral.at(-z)[1] * (1 - 0.1 * Math.sin(k)));
  const noisy = fitViscosity({ tau: 0.1, f, center, thickness, east, north, reach: 200 });
  close(noisy.viscosity, 0.02, 0.1 * 0.02, 'a 10 percent ripple moves the fitted viscosity by less than 10 percent');
  assert.ok(noisy.rms > 0 && noisy.rms < 0.1 * Math.abs(spiral.surface));
});

test('the tracks integrate the series, the period is the one the case uses, and the model’s own Ekman layer is a spiral that a constant viscosity fits', () => {
  const n = 101, time = Float64Array.from({ length: n }, (_, s) => 100 * s), w = 2 * Math.PI / 10000;
  const circling = { time, topEast: time.map((t) => Math.cos(w * t)), topNorth: time.map((t) => Math.sin(w * t)), mixedEast: time.map(() => 0.25), mixedNorth: time.map(() => 0), deepEast: time.map(() => -0.5), deepNorth: time.map(() => 0) };
  const { top, mixed, deep } = tracks(circling);
  close(Math.hypot(top[2 * n - 2], top[2 * n - 1]), 0, 2e-2, 'a whole inertial circle closes');
  close(mixed[2 * n - 2], 0.25 * 10000, 1e-9, 'the mixed layer runs a quarter meter a second east');
  close(deep[2 * n - 2], -0.5 * 10000, 1e-9, 'and the water below half a meter a second west');

  const sim = ekmanCase({ latitude: 45, days: 4 }), frame = runCase(sim);
  close(frame.inertialPeriod, inertialPeriod(45), 1e-9, 'inertial period');
  close(inertialPeriod(45) / 3600, 16.92, 0.005, 'inertial period in hours');
  const path = tracks(frame.series), last = frame.series.time.length - 1;
  let east = 0;
  for (let s = 1; s <= last; s++) east += 0.5 * (frame.series.time[s] - frame.series.time[s - 1]) * (frame.series.mixedEast[s] + frame.series.mixedEast[s - 1]);
  close(path.mixed[2 * last], east, 1e-6, 'the mixed layer’s eastward run');
  assert.ok(path.mixed[2 * last + 1] < -300, `the mixed layer drifts south in the north: ${path.mixed[2 * last + 1]} m`);
  assert.ok(Math.hypot(path.top[2 * last], path.top[2 * last + 1]) > Math.hypot(path.mixed[2 * last], path.mixed[2 * last + 1]), 'the top layer goes farther than the mixed layer');

  const p = frame.profile, fit = fitViscosity({ tau: frame.stressMean, f: frame.fMean, center: p.center, thickness: p.thickness, east: p.east, north: p.north, reach: frame.reach });
  const spiral = ekmanSpiral({ tau: frame.stressMean, f: frame.fMean, viscosity: fit.viscosity });
  console.log(`fitted viscosity ${fit.viscosity.toFixed(4)} m²/s, rms ${(100 * fit.rms).toFixed(2)} cm/s, surface current ${(100 * spiral.surface).toFixed(2)} cm/s in theory against ${(100 * Math.hypot(p.east[0], p.north[0])).toFixed(2)} cm/s in the top layer; model's mean viscosity ${frame.meanViscosity.toFixed(4)}; reversal ${spiral.reversal.toFixed(1)} m, e-folding ${spiral.depth.toFixed(1)} m`);
  assert.ok(fit.viscosity > 0.015 && fit.viscosity < 0.022, `fitted viscosity ${fit.viscosity}`);
  assert.ok(fit.rms < 0.1 * Math.abs(spiral.surface), `the spiral fits to ${fit.rms} m/s against its ${spiral.surface} m/s surface current`);
  close(Math.hypot(p.east[0], p.north[0]), Math.abs(spiral.surface), 0.25 * Math.abs(spiral.surface), 'the top layer’s current is near the spiral’s surface current');
  assert.ok(fit.viscosity > frame.meanViscosity, 'the fit exceeds the boundary-layer mean, which includes the layer’s weakly mixed base');
});
