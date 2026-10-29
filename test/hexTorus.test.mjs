import { test } from 'node:test';
import assert from 'node:assert/strict';
import { gradient, divergence, tangential } from '../js/dynamics/operators.module.js';
import { hexTorus, threeColors, broadHigh, linearWater } from '../explain/hexTorus.module.js';

const SPACING = 60e3, F = 2 * 7.2921e-5 * Math.SQRT1_2;

test('on the wrapped patch the model’s sum rebuilds a uniform wind along every edge exactly, and its weights hand the wind no energy', () => {
  const { mesh } = hexTorus(9, 6, SPACING), u = new Float64Array(mesh.nEdges);
  for (let e = 0; e < mesh.nEdges; e++) u[e] = 3 * mesh.nEdge[3 * e] + 4 * mesh.nEdge[3 * e + 1];
  const along = tangential(mesh, u);
  for (let e = 0; e < mesh.nEdges; e++) assert.ok(Math.abs(along[e] - (4 * mesh.nEdge[3 * e] - 3 * mesh.nEdge[3 * e + 1])) < 1e-12, `edge ${e}: ${along[e]}`);
  const r = Float64Array.from({ length: mesh.nEdges }, (_, e) => Math.sin(1.7 * e) + Math.cos(0.3 * e * e)), turned = tangential(mesh, r);
  assert.ok(Math.abs(r.reduce((a, v, e) => a + v * turned[e], 0)) < 1e-12 * r.reduce((a, v) => a + v * v, 0));
});

test('the three-color pattern has one cell of each color at every corner, drives a wind across every edge and none along any, and spreads at −6/d²', () => {
  for (const [cols, rows] of [[9, 6], [6, 8]]) {
    const torus = hexTorus(cols, rows, SPACING), { mesh, color, neighbor } = torus;
    for (let i = 0; i < mesh.nCells; i++) for (let k = 0; k < 6; k++) {
      const a = neighbor(i, k), b = neighbor(i, (k + 1) % 6);
      assert.equal(new Set([color[i], color[a], color[b]]).size, 3, `corner ${k} of cell ${i}`);
    }
    const p = threeColors(torus), push = gradient(mesh, p), along = tangential(mesh, push), spread = divergence(mesh, push);
    assert.ok(Math.min(...Array.from(push, Math.abs)) > 0.9 / SPACING);
    assert.ok(Math.max(...Array.from(along, Math.abs)) < 1e-12 / SPACING);
    for (let i = 0; i < mesh.nCells; i++) assert.ok(Math.abs(spread[i] + 6 / SPACING ** 2 * p[i]) < 1e-12 / SPACING ** 2);
  }
});

test('let go in the sea, the three-color pattern sloshes every 2πd/(√6 c) with the planet turning or stopped, and a broad high holds its pressure only while the planet turns', () => {
  const torus = hexTorus(9, 6, SPACING), { mesh } = torus, dt = 600, c = 2;
  const run = (start, f, hours) => {
    const step = linearWater(mesh, { c, f, rho: 1025 }), state = { p: Float64Array.from(start), u: new Float64Array(mesh.nEdges) }, norm = start.reduce((a, v) => a + v * v, 0);
    return Array.from({ length: hours * 3600 / dt }, () => { step(state, dt); return state.p.reduce((a, v, i) => a + v * start[i], 0) / norm; });
  };
  const three = threeColors(torus), turning = run(three, F, 48), stopped = run(three, 0, 48), period = 2 * Math.PI * SPACING / (Math.sqrt(6) * c);
  turning.forEach((a, k) => {
    assert.ok(Math.abs(a - stopped[k]) < 1e-12);
    assert.ok(Math.abs(a - Math.cos(2 * Math.PI * (k + 1) * dt / period)) < 1e-4);
  });
  const high = broadHigh(torus, 1.2 * SPACING), held = run(high, F, 96).slice(-288), loose = run(high, 0, 96).slice(-288);
  assert.ok(Math.min(...held) > 0.6, `turning: ${Math.min(...held)}`);
  assert.ok(Math.min(...loose) < 0, `stopped: ${Math.min(...loose)}`);
});
