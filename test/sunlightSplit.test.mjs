import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { readColumn, K } from '../explain/radiativeColumn.module.js';
import { createSplit, lapseTropopause, GROUNDS } from '../explain/sunlightSplit.module.js';

const MLS = readColumn(JSON.parse(readFileSync(new URL('../explain/data/mls_bl36.json', import.meta.url))));
const split = createSplit(MLS);
const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const MU_QUARTER = Math.asin(0.25) * 180 / Math.PI;

test('the six parts add up to the sunlight arriving, for every ground and sun', () => {
  for (const ground of [...Object.keys(GROUNDS), 0, 1]) for (const elevation of [90, 47, 30, 12, 5]) for (const ozone of [true, false]) {
    const s = split({ elevation, ground, ozone }), what = JSON.stringify({ elevation, ground, ozone });
    close(s.arriving, 1362 * Math.sin(elevation * Math.PI / 180), 1e-9, `arriving ${what}`);
    close(s.byAir + s.byGround + s.ozone + s.gases + s.aerosol + s.ground, s.arriving, 1e-9 * s.arriving, `sum ${what}`);
    close(s.ozone + s.gases + s.aerosol, s.air, 1e-9 * s.arriving, `the air's absorbers ${what}`);
    assert.ok(s.byGround >= -1e-9, `the ground adds reflection ${what}`);
    if (!ozone) assert.equal(s.ozone, 0, `no ozone ${what}`);
    if (ground === 0) close(s.byGround, 0, 1e-9, `a black ground reflects nothing ${what}`);
  }
});

test('reflected shares and the air’s absorption over sea with the sun overhead, at 30° and at μ = 0.25', () => {
  const cases = [[90, 0.0707, 0.0490, 49.05, 35.75], [30, 0.1383, 0.0831, 35.29, 25.29], [MU_QUARTER, 0.2438, 0.1294, 26.39, 18.84]];
  for (const [elevation, reflected, byAir, air, ozone] of cases) {
    const s = split({ elevation });
    close(s.reflected / s.arriving, reflected, 5e-4, `reflected share at ${elevation.toFixed(1)}°`);
    close(s.byAir / s.arriving, byAir, 5e-4, `reflected by the air at ${elevation.toFixed(1)}°`);
    close(s.air, air, 0.05, `absorbed in the air at ${elevation.toFixed(1)}°`);
    close(s.ozone, ozone, 0.05, `by ozone at ${elevation.toFixed(1)}°`);
  }
  const overhead = split();
  close(overhead.ground, 1216.6, 0.1, 'absorbed by the ground, sun overhead');
  close(overhead.albedo[0], 0.024, 0.0005, 'the sea\'s direct albedo, sun overhead');
  close(split({ elevation: 5 }).albedo[0], 0.321, 0.001, 'the sea\'s direct albedo, sun 5° up');
});

test('ozone heats the air near 1.6 hPa at 26.5 K/day and the lowest 10 km at 0.12 to 0.16', () => {
  const overhead = split();
  close(overhead.peak.rate, 26.49, 0.01, 'peak heating');
  close(overhead.peak.pressure, 161, 1, 'peak pressure, Pa');
  close(overhead.lower.low, 0.122, 0.001, 'least heating below 10 km');
  close(overhead.lower.high, 0.158, 0.001, 'most heating below 10 km');
  close(split({ elevation: 30 }).peak.pressure, 65, 1, 'a lower sun heats higher up');
  assert.ok(split({ ozone: false }).peak.rate < 0.3, 'without ozone the stratosphere barely warms');
  for (let k = 0; k < K; k++) assert.ok(overhead.heating[k] > 0, `layer ${k} absorbs some sunlight`);
});

test('the tropopause is where the air stops cooling by 2 K per km', () => {
  const height = MLS.height, pressure = MLS.pressure;
  for (const top of [9000, 11000, 16000]) {
    const T = Float64Array.from(height, (z) => 288 - 6.5e-3 * Math.min(z, top));
    close(lapseTropopause(T, height, pressure).height, top, 600, `a 6.5 K/km troposphere ${top / 1000} km deep`);
  }
  const mls = lapseTropopause(MLS.temperature, height, pressure);
  assert.ok(mls.height > 12000 && mls.height < 15000, `midlatitude summer at ${mls.height} m`);
  assert.ok(mls.pressure > 12000 && mls.pressure < 20000, `midlatitude summer at ${mls.pressure} Pa`);
  assert.equal(lapseTropopause(Float64Array.from(height, (z) => 288 - 9.8e-3 * z), height, pressure), null, 'a column cooling all the way up has none');
});
