import { test } from 'node:test';
import assert from 'node:assert/strict';
import { caseMesh, ekman, DEG } from '../explain/oceanCases.module.js';
import { ekmanSpiral, fitViscosity, columnDays } from '../explain/ekmanTheory.module.js';

const N = 16;
const mesh = caseMesh(N), C = mesh.nCells;
const between = (v, lo, hi, what) => assert.ok(v >= lo && v <= hi, `${what}: ${v} not within ${lo} to ${hi}`);

test('at 5°N the Ekman balance has weakened: the transport turns well short of a right angle and runs a third above the ring’s τ/(ρ₀f), and the boundary layer is deeper', () => {
  const low = ekman({ N, latitude: 5, days: columnDays(5) });
  console.log(`5°N: ${low.transport.magnitude.toFixed(3)} m²/s at ${low.transport.angle.toFixed(1)}° (${low.transport.ratio.toFixed(3)} of τ/(ρ₀f) = ${low.theory.toFixed(3)}); top layer ${low.top.angle.toFixed(1)}°; mixed ${low.mixedDepth.toFixed(1)} m, boundary ${low.boundaryDepth.toFixed(1)} m; cells ${low.ringCells}`);
  const share = (r) => { const q = fitViscosity({ tau: r.stressMean, f: r.fMean, center: r.profile.center, thickness: r.profile.thickness, east: r.profile.east, north: r.profile.north, reach: r.reach }); return q.rms / Math.abs(ekmanSpiral({ tau: r.stressMean, f: r.fMean, viscosity: q.viscosity }).surface); };
  console.log(`the best constant viscosity misses the layers by ${(100 * share(low)).toFixed(1)} percent of its surface current at 5°N; boundary layer ${low.boundaryDepth.toFixed(2)} m`);
  between(share(low), 0.15, 0.25, 'misfit at 5°N');
  between(low.boundaryDepth, 49, 55, 'boundary layer at 5°N');
  between(low.transport.ratio, 1.25, 1.5, 'Ekman transport against τ/(ρ₀f) at 5°N'); between(low.transport.angle, 62, 78, 'angle at 5°N');
  const lats = Array.from(mesh.latCell, (v) => v / DEG), ring = [];
  for (let i = 0; i < C; i++) if (Math.abs(lats[i] - 5) < 0.5 * (mesh.dcEdge.reduce((s, d) => s + d, 0) / mesh.nEdges) / mesh.radius / DEG) ring.push(lats[i]);
  between(Math.min(...ring), 2.5, 3.5, 'the 5°N ring reaches down to about 3°N'); between(Math.max(...ring), 6.5, 7.5, 'and up to about 7°N');
});
