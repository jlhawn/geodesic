import { test } from 'node:test';
import assert from 'node:assert/strict';
import { caseMesh, ekman, PUMP_DEPTH, DEG, DAY } from '../explain/oceanCases.module.js';

const N = 16;
const mesh = caseMesh(N), C = mesh.nCells;
const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const between = (v, lo, hi, what) => assert.ok(v >= lo && v <= hi, `${what}: ${v} not within ${lo} to ${hi}`);

test('over the planet the trades lift water at the equator and sink it in the belts near 15° and 30°, and the westerlies’ Ekman transport is τ/(ρ₀f)', () => {
  const trades = ekman({ N, view: 'planet', wind: 'trades' }), westerlies = ekman({ N, view: 'planet', wind: 'westerlies' });
  const bin = (r, lat) => r.zonal.latitude.findIndex((l) => l === lat);
  const w = (r, lat) => r.zonal.upwelling[bin(r, lat)] * DAY;
  const belt = (r, from, to) => { let sum = 0, area = 0; for (let i = 0; i < C; i++) { const lat = Math.abs(mesh.latCell[i]) / DEG; if (lat >= from && lat < to) { sum += mesh.areaCell[i] * r.upwelling[i]; area += mesh.areaCell[i]; } } return sum / area * DAY; };
  const ratio = (r, lat) => r.zonal.ekmanNorth[bin(r, lat)] / r.zonal.theory[bin(r, lat)];
  console.log(`trades, as the readout averages them over the belts of both hemispheres: ${belt(trades, 10, 20).toFixed(3)} m/day at 10–20°, ${belt(trades, 20, 30).toFixed(3)} at 20–30°, ${belt(trades, 30, 40).toFixed(3)} at 30–40°; Ekman transport against τ/(ρ₀f) at 42.5° and 47.5°: ${[42.5, 47.5].map((lat) => ratio(trades, lat).toFixed(3)).join(', ')}`);
  between(belt(trades, 10, 20), -0.12, -0.03, 'sinking at 10–20°, m/day'); between(belt(trades, 20, 30), -0.07, -0.03, 'sinking at 20–30°, m/day'); between(belt(trades, 30, 40), -0.08, -0.04, 'sinking at 30–40°, m/day');
  for (const lat of [42.5, 47.5]) close(ratio(trades, lat), 1, 0.05, `Ekman transport under the trades and westerlies at ${lat}°`);
  console.log(`trades: ${PUMP_DEPTH} m down the equator rises ${(trades.equatorUpwelling * DAY).toFixed(3)} m/day; 12.5° ${w(trades, 12.5).toFixed(3)}, 17.5° ${w(trades, 17.5).toFixed(3)}, 22.5° ${w(trades, 22.5).toFixed(3)}, 27.5° ${w(trades, 27.5).toFixed(3)}, 32.5° ${w(trades, 32.5).toFixed(3)}; westerlies only: ${(westerlies.equatorUpwelling * DAY).toFixed(3)} m/day`);
  for (const lat of [42.5, 47.5, -42.5, -47.5]) {
    const b = bin(westerlies, lat), got = westerlies.zonal.ekmanNorth[b], want = westerlies.zonal.theory[b];
    console.log(`westerlies at ${lat}°: ${got.toFixed(3)} m²/s north against ${want.toFixed(3)}`);
    close(got / want, 1, 0.06, `zonal-mean Ekman transport at ${lat}°`);
  }
  between(trades.equatorUpwelling * DAY, 1.4, 2.1, 'equatorial upwelling under the trades, m/day');
  assert.ok(Math.abs(westerlies.equatorUpwelling * DAY) < 0.02, 'no wind at the equator, no upwelling');
  close(trades.averagedOver, 2 * DAY, 1350, 'the window');
  assert.ok(Math.max(...trades.boundaryDepth) < PUMP_DEPTH / 2, 'the boundary layer is far above the depth the pumping is read at');
  for (const name of ['ekmanEast', 'ekmanNorth', 'upwelling', 'eta', 'mixedDepth', 'boundaryDepth', 'stress']) assert.equal(trades[name].length, C, name);
});
