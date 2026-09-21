import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { interpolationWeights, regridState } from '../js/physics/regrid.module.js';
import { P0 } from '../js/dynamics/sigmaCore.module.js';

const coarse = createModel(new Grid(8));
const fine = createModel(new Grid(16));
const analytic = (x, y, z) => Math.sin(2 * Math.atan2(z, Math.hypot(x, y))) * Math.cos(Math.atan2(y, x)) + 0.5 * z;

test('interpolating a mesh onto itself is the identity', () => {
  const { cells, weights } = interpolationWeights(coarse.mesh, coarse.mesh.xCell);
  for (let i = 0; i < coarse.mesh.nCells; i++) {
    let value = 0;
    for (let m = 0; m < 3; m++) value += weights[3 * i + m] * (cells[3 * i + m] === i ? 1 : 0);
    assert.ok(Math.abs(value - 1) < 1e-9);
  }
});

test('a smooth field regridded from N=8 to N=16 matches the analytic field to second order', () => {
  const sampled = Float64Array.from({ length: coarse.mesh.nCells }, (_, i) => analytic(coarse.mesh.xCell[3 * i], coarse.mesh.xCell[3 * i + 1], coarse.mesh.xCell[3 * i + 2]));
  const { cells, weights } = interpolationWeights(coarse.mesh, fine.mesh.xCell);
  let worst = 0;
  for (let i = 0; i < fine.mesh.nCells; i++) {
    let value = 0;
    for (let m = 0; m < 3; m++) value += weights[3 * i + m] * sampled[cells[3 * i + m]];
    worst = Math.max(worst, Math.abs(value - analytic(fine.mesh.xCell[3 * i], fine.mesh.xCell[3 * i + 1], fine.mesh.xCell[3 * i + 2])));
  }
  console.log(`regrid N=8→16 max error ${worst.toExponential(2)} on a field of amplitude ~1.5`);
  assert.ok(worst < 0.06);
});

test('a state regridded to a finer mesh keeps its mass, its profile, and its solid-body rotation', () => {
  const K = coarse.core.K, C = coarse.mesh.nCells, E = coarse.mesh.nEdges;
  const omega = 1e-5;
  const pi = Float64Array.from({ length: C }, (_, i) => P0 + 1000 * coarse.mesh.xCell[3 * i + 2]);
  const theta = new Float64Array(K * C);
  for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) theta[k * C + i] = 300 + 200 * (1 - coarse.core.sigmaMid[k]) + 5 * coarse.mesh.xCell[3 * i + 2];
  const u = new Float64Array(K * E);
  for (let k = 0; k < K; k++) for (let e = 0; e < E; e++) u[k * E + e] = coarse.mesh.radius * omega * (-coarse.mesh.xEdge[3 * e + 1] * coarse.mesh.nEdge[3 * e] + coarse.mesh.xEdge[3 * e] * coarse.mesh.nEdge[3 * e + 1]);
  const surfaceT = Float64Array.from({ length: C }, (_, i) => 288 - 30 * coarse.mesh.xCell[3 * i + 2] ** 2);
  const [fPi, fTheta, fU, fSurfaceT] = regridState(coarse, fine, [pi, theta, u, surfaceT]);
  const fm = fine.mesh;
  let massCoarse = 0, massFine = 0, area = 0, worstU = 0, worstPi = 0;
  for (let i = 0; i < C; i++) massCoarse += coarse.mesh.areaCell[i] * pi[i];
  for (let i = 0; i < fm.nCells; i++) { massFine += fm.areaCell[i] * fPi[i]; area += fm.areaCell[i]; worstPi = Math.max(worstPi, Math.abs(fPi[i] - (P0 + 1000 * fm.xCell[3 * i + 2]))); }
  for (let e = 0; e < fm.nEdges; e++) {
    const exact = fm.radius * omega * (-fm.xEdge[3 * e + 1] * fm.nEdge[3 * e] + fm.xEdge[3 * e] * fm.nEdge[3 * e + 1]);
    worstU = Math.max(worstU, Math.abs(fU[(K - 1) * fm.nEdges + e] - exact));
  }
  console.log(`regrid state: mass ratio ${(massFine / massCoarse).toFixed(6)}, max ps error ${worstPi.toFixed(2)} Pa, max wind error ${worstU.toFixed(3)} m/s of ${(fm.radius * omega).toFixed(1)}`);
  assert.ok(Math.abs(massFine / massCoarse - 1) < 1e-3);
  assert.ok(worstPi < 20);
  assert.ok(worstU < 0.05 * fm.radius * omega);
  assert.ok(fTheta.every((t) => t > 250 && t < 600) && fSurfaceT.every((t) => t > 250 && t < 300));
});

test('tile sampling keeps a step field exact and, masked, never draws on excluded cells', async () => {
  const { sampleTiles, regridCellField } = await import('../js/physics/regrid.module.js');
  const sm = coarse.mesh;
  const step = Float64Array.from({ length: sm.nCells }, (_, i) => (sm.latCell[i] > 0 ? 1 : 0));
  const tiles = sampleTiles(coarse, fine, step);
  assert.ok(tiles.every((v) => v === 0 || v === 1), 'a tile-sampled step field has no intermediate values');
  const south = Uint8Array.from(step, (v) => 1 - v);
  const poisoned = Float64Array.from(step, (v, i) => (south[i] ? 10 + sm.latCell[i] : -100));
  const masked = sampleTiles(coarse, fine, poisoned, south, undefined, { reach: 3e6, fill: -1 });
  assert.ok(masked.every((v) => v > 0 || v === -1), 'masked sampling takes admitted tiles or the fill');
  assert.ok(masked.some((v) => v === -1), 'points far from every admitted tile take the fill');
  const near = fine.mesh.latCell.map((lat, n) => (lat > 0 && lat < 0.2 ? masked[n] : 5));
  assert.ok(near.every((v) => v > 0), 'excluded points near the boundary take the nearest admitted tile');
  const interpolated = regridCellField(coarse, fine, poisoned, undefined, south);
  for (let n = 0; n < interpolated.length; n++) if (fine.mesh.latCell[n] < 0.15) assert.ok(interpolated[n] > 0, "masked interpolation draws only on admitted cells within reach");
});
