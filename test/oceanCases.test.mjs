import { test } from 'node:test';
import assert from 'node:assert/strict';
import * as original from './helpers/layered.mjs';
import { caseMesh, zonalWindOnEdges, internalWaveSpeed, ekmanLayer, meanWithin, EKMAN_REACH, RHO, G, DAY } from '../explain/oceanCases.module.js';
import { ekmanSpiral, coriolis } from '../explain/ekmanTheory.module.js';
import { HYBRID_LAYERS } from '../js/ocean/column.module.js';

const N = 16;
const mesh = caseMesh(N), C = mesh.nCells;
const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const sameArray = (a, b, what) => { assert.equal(a.length, b.length, `${what} length`); for (let n = 0; n < a.length; n++) assert.equal(a[n], b[n], `${what}[${n}]`); };

test('the case mesh and wind helper give what test/helpers/layered.mjs gives', () => {
  const speed = (lat) => 8 * Math.cos(3 * lat);
  sameArray(zonalWindOnEdges(mesh, speed), original.zonalWindOnEdges(original.mesh, speed), 'zonal wind');
  for (const name of ['xCell', 'xEdge', 'areaCell', 'dcEdge', 'fEdge', 'latCell']) sameArray(mesh[name], original.mesh[name], name);
});

test('the Ekman layer is read against the water below it, and the stack’s gravity-wave speed follows its stratification', () => {
  const tau = 0.1, f = coriolis(45), spiral = ekmanSpiral({ tau, f, viscosity: 0.02 });
  const thick = Float64Array.from({ length: 80 }, () => 5), east = new Float64Array(80), north = new Float64Array(80);
  for (let k = 0; k < 80; k++) { const [u, v] = spiral.at(-(k + 0.5) * 5); east[k] = 0.004 + u; north[k] = -0.002 + v; }
  const layer = ekmanLayer(thick, east, north, 40);
  close(layer.depth, EKMAN_REACH * 40, 1e-12, 'depth the transport is summed over');
  close(layer.east, 0.004, 5e-4, 'the reference current east'); close(layer.north, -0.002, 5e-4, 'the reference current north');
  close(layer.transportNorth / spiral.transport[1], 1, 0.06, 'the truncated spiral carries nearly all of τ/(ρ₀f)');
  close(layer.transportEast, 0, 0.05 * Math.abs(spiral.transport[1]), 'and no transport along the wind');
  close(meanWithin([10, 10, 10], [1, 2, 3], 15), 4 / 3, 1e-12, 'mean within 15 m of three 10 m layers');

  const c2 = internalWaveSpeed([100, 300], [1024, 1026]);
  close(c2, Math.sqrt(G * (2 / RHO) * 100 * 300 / 400), 1e-9, 'two layers: √(g′ h₁h₂/(h₁+h₂))');
  const n = 400, uniform = internalWaveSpeed(new Array(n).fill(10), Array.from({ length: n }, (_, k) => 1024 + 1e-3 * (k + 0.5) * 10));
  close(uniform, Math.sqrt(G / RHO * 1e-3) * 4000 / Math.PI, 0.01 * uniform, 'uniform stratification: N D/π');
  close(internalWaveSpeed([100, 300, 0.5, 0.5], [1024, 1026, 1026.00001, 1026.00002]), c2, 2e-3 * c2, 'a meter of vanished layers of the same water joins the layer above it');
});

test('the worker streams a layout and frames whose transferred buffers are fresh, and reports an unknown case', async () => {
  const messages = [];
  globalThis.self = { postMessage: (message, transfer = []) => messages.push({ message, transfer }) };
  await import('../explain/ocean.worker.js');
  const send = (data) => self.onmessage({ data });
  send({ type: 'start', run: 7, kind: 'ekman', options: { days: 0.25, N } });
  assert.deepEqual(messages.map((m) => m.message.type), ['layout', 'frame']);
  const [layout, first] = messages.map((m) => m.message);
  assert.equal(layout.run, 7); assert.equal(layout.nCells, C); assert.equal(layout.latitude.length, C); assert.equal(layout.layers, HYBRID_LAYERS.K);
  assert.equal(first.time, 0);
  for (let n = 0; n < 60 && !messages[messages.length - 1].message.done; n++) send({ type: 'advance', seconds: 3 * 3600 });
  const last = messages[messages.length - 1].message;
  assert.equal(last.type, 'frame'); assert.ok(last.done); close(last.time, 0.25 * DAY, 1350, 'stops at the end');
  for (const { transfer } of messages) {
    assert.equal(new Set(transfer).size, transfer.length, 'each buffer once');
    for (const buffer of transfer) assert.ok(buffer instanceof ArrayBuffer && !(buffer instanceof SharedArrayBuffer), 'only plain buffers move');
  }
  send({ type: 'start', run: 8, kind: 'gyre' });
  const error = messages[messages.length - 1].message;
  assert.equal(error.type, 'error'); assert.equal(error.run, 8); assert.match(error.message, /no ocean case called gyre/);
  delete globalThis.self;
});
