// A maximin Latin hypercube of `points` points over `parameters` (500
// seeded trials, the one whose least distance in the unit cube is largest),
// with the defaults as point 0, written to `file` once and read back after.
import { existsSync, readFileSync, writeFileSync } from 'node:fs';

function mulberry32(a) {
  return () => { a |= 0; a = (a + 0x6d2b79f5) | 0; let t = Math.imul(a ^ (a >>> 15), 1 | a); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
}

function latinHypercube(n, d, random) {
  const columns = Array.from({ length: d }, () => {
    const order = [...Array(n).keys()];
    for (let i = n - 1; i > 0; i--) { const j = Math.floor(random() * (i + 1)); [order[i], order[j]] = [order[j], order[i]]; }
    return order.map((k) => (k + random()) / n);
  });
  return Array.from({ length: n }, (_, i) => columns.map((c) => c[i]));
}

const minimumDistance = (u) => {
  let least = Infinity;
  for (let i = 0; i < u.length; i++) for (let j = i + 1; j < u.length; j++) least = Math.min(least, u[i].reduce((s, x, k) => s + (x - u[j][k]) ** 2, 0));
  return Math.sqrt(least);
};

export function readDesign(file, parameters, points, seed) {
  if (!existsSync(file)) {
    const random = mulberry32(seed);
    let best = null, bestDistance = -1;
    for (let trial = 0; trial < 500; trial++) {
      const u = latinHypercube(points, parameters.length, random), distance = minimumDistance(u);
      if (distance > bestDistance) { best = u; bestDistance = distance; }
    }
    const rows = [parameters.map((p) => p.base), ...best.map((u) => parameters.map((p, k) => Number((p.low + u[k] * (p.high - p.low)).toPrecision(4))))];
    writeFileSync(file, ['point,' + parameters.map((p) => p.key).join(','), ...rows.map((r, i) => `${i},${r.join(',')}`)].join('\n') + '\n');
    console.log(`design: ${points} points, seed ${seed}, least distance in the unit cube ${bestDistance.toFixed(3)}`);
  }
  return readFileSync(file, 'utf8').trim().split('\n').slice(1).map((line) => {
    const [point, ...values] = line.split(',').map(Number);
    return { point, ...Object.fromEntries(parameters.map((p, k) => [p.key, values[k]])) };
  });
}
