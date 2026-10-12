import { test } from 'node:test';
import assert from 'node:assert/strict';
import * as original from './helpers/layered.mjs';
import {
  caseMesh, zonalWindOnEdges, totalHeatSalt, northwardMixedTransport, uniformOcean, deepen, buriedBump, interfaceRange,
  ekman, lensCase, bumpCase, runCase, RHO, G, DEG, DAY,
} from '../explain/oceanCases.module.js';

const mesh = caseMesh(8), C = mesh.nCells;
const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const sameArray = (a, b, what) => { assert.equal(a.length, b.length, `${what} length`); for (let n = 0; n < a.length; n++) assert.equal(a[n], b[n], `${what}[${n}]`); };

test('the copied helpers give what test/helpers/layered.mjs gives on the same mesh', () => {
  const speed = (lat) => 8 * Math.cos(3 * lat);
  sameArray(zonalWindOnEdges(mesh, speed), original.zonalWindOnEdges(original.mesh, speed), 'zonal wind');

  const ours = uniformOcean(mesh, { densities: [1026, 1027], bottoms: [400], salinities: [34.6, 34.7] });
  const theirs = original.uniformOcean({ densities: [1026, 1027], bottoms: [400], salinities: [34.6, 34.7] });
  sameArray(ours.ocean.h, theirs.ocean.h, 'uniform ocean h');
  sameArray(ours.surfaceT, theirs.surfaceT, 'uniform surface temperature');
  assert.deepEqual(totalHeatSalt(ours.ocean, mesh), totalHeatSalt(theirs.ocean, original.mesh));
  assert.deepEqual(totalHeatSalt(ours.ocean, mesh), original.totalHeatSalt(theirs.ocean, original.mesh));

  const cells = [5, 77, 300];
  deepen(ours, cells, 180, 0.05);
  original.deepen(theirs, cells, 180, 0.05);
  sameArray(ours.ocean.h, theirs.ocean.h, 'deepened h');
  sameArray(ours.ocean.Q, theirs.ocean.Q, 'deepened heat');
  for (const k of [0, 1, 2]) assert.equal(interfaceRange(ours.ocean, k), original.interfaceRange(theirs.ocean, k));

  const wind = zonalWindOnEdges(mesh, (lat) => 0.1 * Math.cos(lat));
  for (let n = 0; n < 6; n++) { ours.ocean.advance(ours.surfaceT, ours.ice, ours.flux, wind, 1350); theirs.ocean.advance(theirs.surfaceT, theirs.ice, theirs.flux, wind, 1350); }
  const a = northwardMixedTransport(ours.ocean, mesh), b = original.northwardMixedTransport(theirs.ocean, original.mesh);
  let polar = 0;
  for (let i = 0; i < C; i++) {
    assert.ok(Number.isFinite(a[i]), `cell ${i}`);
    if (Math.hypot(mesh.xCell[3 * i], mesh.xCell[3 * i + 1]) < 1e-12) { polar++; continue; }
    close(a[i], b[i], 1e-12 * (1 + Math.abs(b[i])), `northward transport in cell ${i}`);
  }
  assert.equal(polar, 2, 'the two cells centered on the poles, where north has no direction');

  const mine = buriedBump(mesh, { everySteps: 1 }), yours = original.buriedBump({ everySteps: 1 });
  mine.reshape(2, 40); yours.reshape(2, 40);
  sameArray(mine.ocean.h, yours.ocean.h, 'buried bump h');
});

test('a band of westerlies at 45°N drives a mixed-layer drift of τ/(ρ₀f) at right angles to the wind, to the right, and to the left at 45°S', () => {
  const north = ekman({ latitude: 45 }), south = ekman({ latitude: -45 });
  const describe = (r) => `${r.ekman.magnitude.toFixed(3)} m²/s at ${r.ekman.angle.toFixed(1)}° (${r.ekman.ratio.toFixed(3)} of τ/(ρ₀f) = ${r.theory.toFixed(3)}); with the water beneath ${r.transport.magnitude.toFixed(3)} m²/s at ${r.transport.angle.toFixed(1)}°`;
  console.log(`45°N: ${describe(north)}; inertial period ${(north.inertialPeriod / 3600).toFixed(2)} h, ${north.ringCells} cells, mixed layer ${north.mixedDepth.toFixed(1)} m`);
  console.log(`45°S: ${describe(south)}`);
  close(north.inertialPeriod / 3600, 16.92, 0.01, 'inertial period at 45° in hours');
  close(north.theory, 0.1 / (RHO * 2 * mesh.omega * Math.sin(45 * DEG)), 0.05, 'the ring\'s τ/(ρ₀f) is near the one at 45°');
  close(north.ekman.ratio, 1, 0.05, 'Ekman transport against τ/(ρ₀f) at 45°N');
  assert.ok(north.ekman.north < 0 && north.ekman.angle > 75 && north.ekman.angle < 95, `the drift is to the right of the wind: ${north.ekman.angle}°`);
  close(south.ekman.ratio, 1, 0.05, 'Ekman transport against τ/(ρ₀f) at 45°S');
  assert.ok(south.ekman.north > 0 && south.ekman.angle < -75 && south.ekman.angle > -95, `the drift is to the left of the wind: ${south.ekman.angle}°`);
  close(north.transport.ratio, 1, 0.1, 'the whole mixed layer\'s transport');
  assert.ok(north.settled && north.averagedOver >= 2 * north.inertialPeriod - 1350 && north.averagedOver <= 2 * north.inertialPeriod + 1350);
  const n = north.series.time.length;
  assert.equal(n, 4 * 64 + 1);
  for (const name of ['east', 'north', 'relativeEast', 'relativeNorth', 'cellEast', 'cellNorth', 'mixedDepth']) assert.equal(north.series[name].length, n, name);
  let fastest = 0;
  for (let s = 0; s < n; s++) fastest = Math.max(fastest, Math.hypot(north.series.relativeEast[s], north.series.relativeNorth[s]));
  const settled = Math.hypot(north.series.relativeEast[n - 1], north.series.relativeNorth[n - 1]);
  console.log(`the slab's speed peaks at ${(100 * fastest).toFixed(2)} cm/s while it circles and drifts at ${(100 * settled).toFixed(2)} cm/s after 4 days`);
  assert.ok(fastest > 1.5 * settled, 'the layer overshoots in an inertial circle before it settles');
  assert.ok(Math.abs(north.below.east[2]) < 0.2 * Math.abs(north.series.east[n - 1]) + 0.01, 'the water beneath moves far more slowly');
});

test('over the planet the trades lift water at the equator and the westerlies push it equatorward at τ/(ρ₀f)', () => {
  const trades = ekman({ view: 'planet', wind: 'trades' }), westerlies = ekman({ view: 'planet', wind: 'westerlies' });
  const bin = (r, lat) => r.zonal.latitude.findIndex((l) => l === lat);
  console.log(`trades: upwelling at the equator ${(trades.equatorUpwelling * DAY).toFixed(2)} m/day; westerlies: ${(westerlies.equatorUpwelling * DAY).toFixed(3)} m/day`);
  for (const lat of [42.5, 47.5, -42.5, -47.5]) {
    const b = bin(westerlies, lat), got = westerlies.zonal.transportNorth[b], want = westerlies.zonal.theory[b];
    console.log(`westerlies at ${lat}°: ${got.toFixed(3)} m²/s north against ${want.toFixed(3)}`);
    close(got / want, 1, 0.1, `zonal-mean Ekman transport at ${lat}°`);
  }
  assert.ok(trades.equatorUpwelling * DAY > 0.3, `equatorial upwelling ${trades.equatorUpwelling * DAY} m/day`);
  assert.ok(Math.abs(westerlies.equatorUpwelling * DAY) < 0.02, 'no wind at the equator, no upwelling');
  close(trades.averagedOver, 2 * DAY, 1350, 'the window');
  for (const name of ['transportEast', 'transportNorth', 'ekmanEast', 'ekmanNorth', 'upwelling', 'eta', 'mixedDepth', 'stress']) assert.equal(trades[name].length, C, name);
});

test('a lens of light water 60 m thick raises the sea surface by Δρ/ρ₀ of its thickness, and the surface current around it is geostrophic while the deep water barely moves', () => {
  const sim = lensCase({ days: 10 });
  const layout = sim.layout(), frames = [];
  for (const day of [1, 5, 10]) { runCase(sim, day); frames.push(sim.frame()); }
  for (const f of frames) console.log(`day ${f.day}: surface ${(100 * f.rise).toFixed(1)} cm up over a thermocline ${f.depression.toFixed(1)} m down (expected ${(100 * f.expected).toFixed(1)} cm); f·u against g∂η/∂n ${f.geostrophic.ratio.toFixed(3)}, correlation ${f.geostrophic.correlation.toFixed(3)}; speeds ${Array.from(f.layerSpeeds, (v) => (100 * v).toFixed(2)).join(', ')} cm/s`);
  const [, day5, day10] = frames;
  close(day10.depression, 59.4, 0.5, 'the thermocline at the centre cell');
  for (const f of [day5, day10]) {
    assert.ok(f.rise > 0.85 * f.expected && f.rise < 1.02 * f.expected, `rise ${f.rise} against ${f.expected}`);
    assert.ok(f.rise > 0.16 && f.rise < 0.21, `rise ${f.rise} m`);
    close(f.geostrophic.ratio, 1, 0.1, 'geostrophic ratio');
    assert.ok(f.geostrophic.correlation > 0.99, `correlation ${f.geostrophic.correlation}`);
    assert.ok(f.layerSpeeds[4] < 0.15 * f.layerSpeeds[0], 'the deepest class barely moves');
  }
  assert.equal(day10.slice.distance.length, layout.slice.cells.length);
  assert.equal(day10.slice.interfaces.length, layout.slice.cells.length * layout.layers);
  assert.equal(sim.ocean.diagnostics().oceanLimited, 0);
});

test('a 1 m bump of sea surface races away at about √(gD) and leaves a small geostrophic eddy in every layer', () => {
  const sim = bumpCase({ days: 1 }), depth = 4000;
  runCase(sim, 8 / 24);
  const wave = sim.frame();
  console.log(`after 8 h the edge of the ring is ${(wave.front / 1e3).toFixed(0)} km out, moving at ${wave.frontSpeed.toFixed(1)} m/s against √(gD) ${wave.gravitySpeed.toFixed(1)} m/s`);
  close(wave.gravitySpeed, Math.sqrt(G * depth), 1e-9, '√(gD)');
  close(wave.frontSpeed / wave.gravitySpeed, 1, 0.1, 'front speed');
  runCase(sim);
  const eddy = sim.frame(), speeds = Array.from(eddy.layerSpeeds);
  console.log(`after a day the centre stands ${(100 * eddy.centreHeight).toFixed(1)} cm up, ${(100 * eddy.kept).toFixed(0)} % of the bump's volume within twice its width; f·u against g∂η/∂n ${eddy.geostrophic.ratio.toFixed(3)}; speeds ${speeds.map((v) => (100 * v).toFixed(2)).join(', ')} cm/s`);
  assert.ok(eddy.centreHeight > 0.1 && eddy.centreHeight < 0.25, `centre ${eddy.centreHeight} m`);
  close(eddy.geostrophic.ratio, 1, 0.1, 'geostrophic eddy');
  assert.ok(Math.max(...speeds) < 1.02 * Math.min(...speeds), 'the same speed in every layer');
  close(eddy.level * 4 * Math.PI * mesh.radius ** 2, Math.PI * 1.5e6 ** 2, 0.05 * Math.PI * 1.5e6 ** 2, 'the volume of the bump spread over the sea');
});

test('the worker streams a layout and frames whose transferred buffers are fresh, and reports an unknown case', async () => {
  const messages = [];
  globalThis.self = { postMessage: (message, transfer = []) => messages.push({ message, transfer }) };
  await import('../explain/ocean.worker.js');
  const send = (data) => self.onmessage({ data });
  send({ type: 'start', run: 7, kind: 'ekman', options: { days: 0.25 } });
  assert.deepEqual(messages.map((m) => m.message.type), ['layout', 'frame']);
  const [layout, first] = messages.map((m) => m.message);
  assert.equal(layout.run, 7); assert.equal(layout.nCells, C); assert.equal(layout.latitude.length, C);
  assert.equal(first.time, 0);
  for (let n = 0; n < 10 && !messages[messages.length - 1].message.done; n++) send({ type: 'advance', seconds: 3 * 3600 });
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
