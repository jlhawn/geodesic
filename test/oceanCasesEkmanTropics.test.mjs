import { test } from 'node:test';
import assert from 'node:assert/strict';
import { ekman } from '../explain/oceanCases.module.js';
import { ekmanSpiral, fitViscosity, columnDays } from '../explain/ekmanTheory.module.js';

const N = 16;
const between = (v, lo, hi, what) => assert.ok(v >= lo && v <= hi, `${what}: ${v} not within ${lo} to ${hi}`);

test('at 15°N the Ekman transport is about τ/(ρ₀f) and turns about a right angle, the boundary layer is 40 m deep and a constant viscosity fits the layers to within a sixth of the surface current', () => {
  const mid = ekman({ N, latitude: 15, days: columnDays(15) });
  console.log(`15°N: ${mid.transport.magnitude.toFixed(3)} m²/s at ${mid.transport.angle.toFixed(1)}° (${mid.transport.ratio.toFixed(3)} of τ/(ρ₀f) = ${mid.theory.toFixed(3)}); top layer ${mid.top.angle.toFixed(1)}°; mixed ${mid.mixedDepth.toFixed(1)} m, boundary ${mid.boundaryDepth.toFixed(1)} m; cells ${mid.ringCells}`);
  const share = (r) => { const q = fitViscosity({ tau: r.stressMean, f: r.fMean, center: r.profile.center, thickness: r.profile.thickness, east: r.profile.east, north: r.profile.north, reach: r.reach }); return q.rms / Math.abs(ekmanSpiral({ tau: r.stressMean, f: r.fMean, viscosity: q.viscosity }).surface); };
  console.log(`the best constant viscosity misses the layers by ${(100 * share(mid)).toFixed(1)} percent of its surface current at 15°N; boundary layer ${mid.boundaryDepth.toFixed(2)} m`);
  between(share(mid), 0.1, 0.18, 'misfit at 15°N');
  between(mid.boundaryDepth, 38, 43, 'boundary layer at 15°N');
  between(mid.transport.ratio, 0.9, 1.1, 'Ekman transport against τ/(ρ₀f) at 15°N'); between(mid.transport.angle, 85, 105, 'angle at 15°N');
});
