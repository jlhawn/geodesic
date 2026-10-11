import { test } from 'node:test';
import assert from 'node:assert/strict';
import { OMEGA, RADIUS, cross, coriolis, frame, launch, parts, step, latitude, longitude } from '../explain/coriolisSphere.module.js';

const D = Math.PI / 180;
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const norm = (a) => Math.hypot(...a);
const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const heading = (s) => { const { east, north } = frame(latitude(s.x), longitude(s.x)); return Math.atan2(dot(s.v, east), dot(s.v, north)); };
const loops = (lat, omega = OMEGA) => {
  const s = launch(lat * D, 0, 0, 1, omega), expected = 2 * Math.PI / Math.abs(2 * omega * Math.sin(lat * D)), crossings = [];
  let previous = dot(s.v, frame(lat * D, 0).north);
  for (let n = 0; n < Math.ceil(4.5 * expected / 60); n++) {
    const t = s.t;
    step(s, 60);
    const now = dot(s.v, frame(latitude(s.x), longitude(s.x)).north);
    if (previous > 0 && now <= 0) crossings.push(t + 60 * previous / (previous - now));
    previous = now;
  }
  return { count: crossings.length, period: (crossings[crossings.length - 1] - crossings[0]) / (crossings.length - 1), expected };
};

test('the Coriolis acceleration is −2Ω × v, at right angles to the axis and to the velocity', () => {
  for (const v of [[10, 0, 0], [0, 10, 0], [3, -4, 12], [-7, 2, -5], [0, 0, 9]]) {
    const a = coriolis(v), expected = cross([0, 0, OMEGA], v).map((c) => -2 * c);
    for (let i = 0; i < 3; i++) close(a[i], expected[i], 1e-15, `component ${i}`);
    close(dot(a, [0, 0, 1]), 0, 1e-18, 'along the axis');
    close(dot(a, v), 0, 1e-15, 'along the velocity');
  }
});

test('at the equator the force has no part along the surface, and its vertical part is 2Ωu', () => {
  for (let h = 0; h < 360; h += 15) {
    const s = launch(0, 0.7, h * D, 25), { surface, vertical } = parts(s.x, s.v), u = 25 * Math.sin(h * D);
    close(norm(surface), 0, 1e-15, `surface part heading ${h}°`);
    close(vertical, 2 * OMEGA * u, 1e-15, `vertical part heading ${h}°`);
  }
});

test('at the pole the force has no vertical part and lies wholly along the surface', () => {
  for (const pole of [90, -90]) for (let h = 0; h < 360; h += 30) {
    const s = launch(pole * D, 1.2, h * D, 25), { surface, vertical } = parts(s.x, s.v);
    close(vertical, 0, 1e-15, `vertical part at ${pole}° heading ${h}°`);
    close(norm(surface), 2 * OMEGA * 25, 1e-15, `surface part at ${pole}° heading ${h}°`);
  }
});

test('a 10 m/s east wind at 45° north is pushed 1.031e-3 m/s² south and 1.031e-3 m/s² up', () => {
  const s = launch(45 * D, -0.4, 90 * D, 10), { up, surface, vertical } = parts(s.x, s.v), { north } = frame(45 * D, -0.4);
  close(norm(surface), 1.031e-3, 5e-7, 'surface part');
  close(dot(surface, north), -norm(surface), 1e-15, 'surface part points south');
  close(vertical, 1.031e-3, 5e-7, 'vertical part, upward');
  close(dot(up, s.x) / RADIUS, 1, 1e-15, 'up is the local vertical');
});

test('the surface part turns the parcel to the right in the north and to the left in the south', () => {
  for (const lat of [10, 30, 45, 60, 80, -10, -30, -45, -60, -80]) for (let h = 0; h < 360; h += 45) {
    const s = launch(lat * D, 2, h * D, 20), { up, surface } = parts(s.x, s.v), f = 2 * OMEGA * Math.sin(lat * D);
    close(norm(surface), Math.abs(f) * 20, 1e-15, `size at ${lat}° heading ${h}°`);
    close(dot(surface, s.v), 0, 1e-15, `at right angles to the velocity at ${lat}° heading ${h}°`);
    assert.equal(Math.sign(dot(cross(s.v, surface), up)), lat > 0 ? -1 : 1, `side at ${lat}° heading ${h}°`);
    const before = heading(s);
    for (let n = 0; n < 10; n++) step(s, 60);
    let turned = heading(s) - before;
    turned -= 2 * Math.PI * Math.round(turned / (2 * Math.PI));
    const geodesic = 20 * Math.sin(h * D) * Math.tan(lat * D) / RADIUS;
    close(turned, (f + geodesic) * 600, Math.abs(f) * 600 * 0.01, `heading change at ${lat}° heading ${h}°`);
  }
});

test('the speed is conserved and the parcel stays on the sphere for three days', () => {
  for (const [lat, h, speed] of [[45, 90, 10], [45, 90, 100], [20, 0, 300], [0, 90, 150], [-60, 200, 300], [89, 45, 50]]) {
    const s = launch(lat * D, 0.3, h * D, speed);
    let worst = 0;
    for (let t = 0; t < 3 * 86400; t += 60) { step(s, 60); worst = Math.max(worst, Math.abs(norm(s.v) / speed - 1)); }
    assert.ok(worst < 1e-6, `speed drift ${worst} at ${lat}°, ${speed} m/s`);
    close(norm(s.x) / RADIUS, 1, 1e-12, 'radius');
    close(dot(s.x, s.v) / (RADIUS * speed), 0, 1e-12, 'velocity along the surface');
  }
});

test('a slow parcel goes round once every 2π/f', () => {
  for (const lat of [45, 30, -60]) {
    const { count, period, expected } = loops(lat);
    assert.ok(count >= 3, `${count} loops at ${lat}°`);
    assert.ok(Math.abs(period / expected - 1) < 1e-3, `period ${period} s against ${expected} s at ${lat}°`);
  }
});

test('the period scales as 1/|Ω|: at twice Earth\'s spin a slow parcel at 45° goes round in 8.46 hours', () => {
  const earth = loops(45).period;
  for (const k of [2, 0.5, -1, -2, 0.1]) for (const lat of [45, -60]) {
    const { count, period, expected } = loops(lat, k * OMEGA);
    assert.ok(count >= 3, `${count} loops at ${lat}°, ${k} × Earth's spin`);
    assert.ok(Math.abs(period / expected - 1) < 1e-3, `period ${period} s against ${expected} s at ${lat}°, ${k} × Earth's spin`);
    if (lat === 45) close(period * Math.abs(k), earth, earth * 1e-3, `period times |k| at ${k} × Earth's spin`);
  }
  close(loops(45, 2 * OMEGA).period / 3600, 8.46, 0.005, 'hours per loop at 45°, twice Earth\'s spin');
});

test('a parcel launched on a planet spinning at another rate keeps that rate when stepped', () => {
  const a = launch(45 * D, 0, 90 * D, 100, 2 * OMEGA), b = launch(45 * D, 0, 90 * D, 100);
  for (let n = 0; n < 600; n++) { step(a, 60); step(b, 60, 2 * OMEGA); }
  assert.deepEqual(a.x, b.x);
  assert.deepEqual(a.v, b.v);
});

test('with no spin nothing turns the parcel: it runs along a great circle at a steady speed', () => {
  for (const [lat, h, speed] of [[45, 90, 100], [20, 0, 300], [-60, 200, 50], [0, 45, 150], [89, 270, 80]]) {
    const s = launch(lat * D, 0.3, h * D, speed, 0), x0 = s.x.map((c) => c / RADIUS), v0 = s.v.map((c) => c / speed), normal = cross(x0, v0);
    assert.equal(norm(parts(s.x, s.v, 0).full), 0, `force at ${lat}° heading ${h}°`);
    let drift = 0, sideways = 0, off = 0;
    for (let n = 1; n <= 3 * 1440; n++) {
      step(s, 60);
      const theta = speed * s.t / RADIUS, along = x0.map((c, i) => RADIUS * (c * Math.cos(theta) + v0[i] * Math.sin(theta)));
      drift = Math.max(drift, Math.abs(norm(s.v) / speed - 1));
      sideways = Math.max(sideways, Math.abs(dot(s.x, normal)) / RADIUS, Math.abs(dot(s.v, normal)) / speed);
      off = Math.max(off, norm(s.x.map((c, i) => c - along[i])) / RADIUS);
    }
    assert.ok(drift < 1e-9, `speed drift ${drift} at ${lat}° heading ${h}°`);
    assert.ok(sideways < 1e-9, `sideways ${sideways} at ${lat}° heading ${h}°`);
    assert.ok(off < 1e-8, `distance from the great circle position ${off} at ${lat}° heading ${h}°`);
  }
});

test('spun backward, the planet turns the parcel to the left in the north and to the right in the south', () => {
  const omega = -OMEGA;
  for (const lat of [10, 30, 45, 60, 80, -10, -30, -45, -60, -80]) for (let h = 0; h < 360; h += 45) {
    const s = launch(lat * D, 2, h * D, 20, omega), { up, surface } = parts(s.x, s.v, omega), f = 2 * omega * Math.sin(lat * D);
    close(norm(surface), Math.abs(f) * 20, 1e-15, `size at ${lat}° heading ${h}°`);
    assert.equal(Math.sign(dot(cross(s.v, surface), up)), lat > 0 ? 1 : -1, `side at ${lat}° heading ${h}°`);
    const before = heading(s);
    for (let n = 0; n < 10; n++) step(s, 60);
    let turned = heading(s) - before;
    turned -= 2 * Math.PI * Math.round(turned / (2 * Math.PI));
    const geodesic = 20 * Math.sin(h * D) * Math.tan(lat * D) / RADIUS;
    close(turned, (f + geodesic) * 600, Math.abs(f) * 600 * 0.01, `heading change at ${lat}° heading ${h}°`);
    assert.equal(Math.sign(turned - geodesic * 600), lat > 0 ? -1 : 1, `turned left in the north at ${lat}° heading ${h}°`);
  }
});

test('the force and both its parts change sign with Ω and grow in proportion to it', () => {
  for (const lat of [0, 20, 45, 70, 90, -35, -90]) for (let h = 0; h < 360; h += 30) {
    const s = launch(lat * D, 0.9, h * D, 30), earth = parts(s.x, s.v);
    for (const k of [2, 0.5, 0, -0.5, -1, -2]) {
      const p = parts(s.x, s.v, k * OMEGA);
      close(p.vertical, k * earth.vertical, 1e-15, `vertical part at ${lat}° heading ${h}°, ${k} × Earth's spin`);
      for (let i = 0; i < 3; i++) {
        close(p.full[i], k * earth.full[i], 1e-15, `force component ${i} at ${lat}° heading ${h}°, ${k} × Earth's spin`);
        close(p.surface[i], k * earth.surface[i], 1e-15, `surface component ${i} at ${lat}° heading ${h}°, ${k} × Earth's spin`);
      }
    }
  }
  const s = launch(45 * D, -0.4, 90 * D, 10);
  close(parts(s.x, s.v, -OMEGA).vertical, -1.031e-3, 5e-7, 'a 10 m/s east wind at 45° north is pushed down when the planet spins backward');
});

test('launched east at 45° north, the loops reach as far south, last as long and drift as far west as the prose says', () => {
  const loops = (speed) => {
    const s = launch(45 * D, 0, 90 * D, speed), tops = [];
    let lowest = 90, lon = 0, last = 0, previous = dot(s.v, frame(45 * D, 0).north);
    for (let n = 0; n < 6 * 1440; n++) {
      const t = s.t;
      step(s, 60);
      const phi = latitude(s.x), here = longitude(s.x), now = dot(s.v, frame(phi, here).north);
      lon += here - last - 2 * Math.PI * Math.round((here - last) / (2 * Math.PI));
      last = here;
      lowest = Math.min(lowest, phi / D);
      if (previous > 0 && now <= 0) tops.push([t + 60 * previous / (previous - now), lon]);
      previous = now;
    }
    const [t0, l0] = tops[0], [t1, l1] = tops[tops.length - 1];
    return { lowest, hours: (t1 - t0) / (tops.length - 1) / 3600, west: -(l1 - l0) / D / ((t1 - t0) / 86400) };
  };
  const slow = loops(10), fast = loops(100), km = slow.west * D * RADIUS * Math.cos(45 * D) / 1000;
  assert.ok(km >= 5 && km <= 8, `10 m/s drifts west ${km} km a day`);
  close(fast.lowest, 22.7, 0.2, 'southernmost latitude at 100 m/s');
  close(fast.hours, 22.4, 0.3, 'hours per loop at 100 m/s');
  close(fast.west, 13.75, 0.3, 'degrees of longitude a day westward at 100 m/s');
});
