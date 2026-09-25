import test from 'node:test';
import assert from 'node:assert/strict';
import { twoFingerMotion } from '../js/gestures.module.js';

const close = (actual, expected, tolerance = 1e-12) => assert.ok(Math.abs(actual - expected) <= tolerance, `${actual} vs ${expected}`);

test('a parallel slide pans by the midpoint shift alone', () => {
  const m = twoFingerMotion({ x: 100, y: 100 }, { x: 200, y: 100 }, { x: 130, y: 80 }, { x: 230, y: 80 });
  close(m.dx, 30); close(m.dy, -20); close(m.scale, 1); close(m.turn, 0);
  close(m.x, 180); close(m.y, 80);
});

test('spreading about a fixed midpoint zooms without panning', () => {
  const m = twoFingerMotion({ x: 150, y: 200 }, { x: 250, y: 200 }, { x: 100, y: 200 }, { x: 300, y: 200 });
  close(m.dx, 0); close(m.dy, 0); close(m.scale, 2); close(m.turn, 0);
});

test('a clockwise turn on screen is positive and wraps across ±π', () => {
  const quarter = twoFingerMotion({ x: 0, y: 0 }, { x: 100, y: 0 }, { x: 0, y: 0 }, { x: 0, y: 100 });
  close(quarter.turn, Math.PI / 2);
  close(quarter.scale, 1);
  const across = twoFingerMotion({ x: 0, y: 0 }, { x: -100, y: 1 }, { x: 0, y: 0 }, { x: -100, y: -1 });
  close(across.turn, 2 * Math.atan2(1, 100));
});

test('coincident fingers leave the zoom alone', () => {
  const m = twoFingerMotion({ x: 5, y: 5 }, { x: 5, y: 5 }, { x: 6, y: 5 }, { x: 8, y: 5 });
  close(m.scale, 1);
});
