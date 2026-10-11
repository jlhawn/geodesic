import { test } from 'node:test';
import assert from 'node:assert/strict';
import { liftParcel } from '../explain/moistParcel.module.js';

globalThis.IntersectionObserver ??= class { observe() {} };
const { reading } = await import('../explain/figures/moistlift.module.js');

const C0 = 273.15;
const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const path = (ground, rh) => liftParcel({ T0: C0 + ground, rh: rh / 100 }, 190e2);

test('the figure opens on a parcel at 15 °C and 70 percent holding 7.47 of the 10.67 g/kg it can hold, with its cloud base at 678 m', () => {
  const P = path(15, 70), s = reading(P, 0);
  close(1000 * s.q, 7.47, 0.005, 'vapor');
  close(1000 * s.limit, 10.67, 0.005, 'the most it can hold');
  close(s.rh, 0.7, 1e-12, 'relative humidity');
  close(s.dew - C0, 9.6, 0.05, 'dew point');
  close(1000 * s.lapse, 9.77, 0.005, 'dry below cloud base');
  close(P.lcl.height, 678, 0.5, 'cloud base');
  assert.equal(reading(P, 650).cloudy, false);
  close(reading(P, 650).rh, 0.985, 0.0005, 'nearly saturated just below cloud base');
  assert.equal(reading(P, 700).cloudy, true);
});

test('lifted to 5.5 km the cloudy parcel is 14.0 K warmer than a dry one and its θ has risen to 305.2 K; at 12 km it cools at nearly the dry rate', () => {
  const P = path(15, 70);
  const mid = reading(P, 5500), top = reading(P, 12000);
  close(mid.p / 100, 498.5, 0.05, 'pressure at 5.5 km');
  close(mid.T - C0, -23.04, 0.005, 'parcel at 5.5 km');
  close(mid.dry - C0, -37.03, 0.005, 'dry parcel at 5.5 km');
  close(1000 * mid.condensate, 6.26, 0.005, 'cloud water made by 5.5 km');
  close(1000 * mid.q, 1.20, 0.005, 'vapor left at 5.5 km');
  close(mid.theta, 305.2, 0.05, 'θ at 5.5 km');
  close(1000 * mid.lapse, 8.08, 0.005, 'cooling rate at 5.5 km');
  close((P.lcl.temperature - mid.T) / (5.5 - P.lcl.height / 1000), 6.5, 0.05, 'mean cooling from cloud base to 5.5 km, K per km');
  for (const [z, rate] of [[1000, 5.42], [3000, 6.41]]) close(1000 * reading(P, z).lapse, rate, 0.005, `cooling rate at ${z} m`);
  close(top.p / 100, 190.8, 0.05, 'pressure at 12 km');
  close(top.T - C0, -80.81, 0.005, 'parcel at 12 km');
  close(top.dry - C0, -93.75, 0.005, 'dry parcel at 12 km');
  close(1000 * top.condensate, 7.46, 0.005, 'cloud water made by 12 km');
  close(top.theta, 308.9, 0.05, 'θ at 12 km');
  close(1000 * top.lapse, 9.76, 0.005, 'cooling rate at 12 km');
});

test('warm damp air from a tropical sea and cold air at freezing', () => {
  const warm = path(27, 80), cold = path(0, 50);
  close(1000 * reading(warm, 0).q, 17.99, 0.005, 'tropical vapor');
  close(1000 * reading(warm, 0).limit, 22.48, 0.005, 'tropical limit');
  close(warm.lcl.height, 452, 0.5, 'tropical cloud base');
  const mid = reading(warm, 5500), top = reading(warm, 12000);
  close(mid.T - C0, -1.63, 0.005, 'tropical parcel at 5.5 km');
  close(mid.T - mid.dry, 25.57, 0.005, 'tropical warmth over a dry parcel at 5.5 km');
  close(1000 * mid.condensate, 11.19, 0.005, 'tropical cloud water by 5.5 km');
  close(top.T - top.dry, 32.92, 0.005, 'tropical warmth over a dry parcel at 12 km');
  close(top.theta, 353.0, 0.05, 'tropical θ at 12 km');
  close(1000 * reading(cold, 0).q, 1.91, 0.005, 'vapor at 0 °C and 50 percent');
  close(cold.lcl.height, 1197, 0.5, 'cold cloud base');
  close(reading(cold, 12000).T - reading(cold, 12000).dry, 3.18, 0.005, 'cold warmth over a dry parcel at 12 km');
  for (const [ground, limit] of [[0, 3.81], [15, 10.67], [30, 26.84], [35, 35.79]]) close(1000 * reading(path(ground, 50), 0).limit, limit, 0.005, `most vapor at ${ground} °C`);
});

test('at 100 percent the cloud base is the ground, and perfectly dry air never makes cloud', () => {
  const fog = path(15, 100), dry = path(15, 0);
  assert.equal(fog.lcl.height, 0);
  assert.equal(reading(fog, 0).cloudy, true);
  assert.equal(dry.lcl, null);
  const s = reading(dry, 5500);
  assert.equal(s.cloudy, false);
  assert.equal(s.condensate, 0);
  assert.ok(Number.isNaN(s.dew));
  assert.equal(s.rh, 0);
});
