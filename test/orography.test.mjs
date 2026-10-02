import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { syntheticTopography, subgridOrography, createGeography, surfaceGeopotential } from '../js/geography.module.js';

const close = (actual, expected, tolerance, what) => assert.ok(Math.abs(actual - expected) <= tolerance * Math.abs(expected), `${what}: ${actual} against ${expected}`);
const deg = Math.PI / 180, R = 6371220;

function ridgeFields(elevationAt) {
  const mesh = buildMesh(new Grid(8));
  const topography = syntheticTopography(720, 1440, elevationAt);
  const geography = createGeography(mesh, topography, { landBridges: {}, seaStraits: {} });
  const phis = surfaceGeopotential(mesh, geography);
  return { mesh, fields: subgridOrography(mesh, topography, Float64Array.from(phis, (p) => p / 9.80616)) };
}

test('the subgrid fields of analytic ridges: standard deviation, slope, anisotropy and orientation', () => {
  const amplitude = 400, wavelength = 4 * deg, k = 2 * Math.PI / wavelength, sampling = Math.sin(k * 0.25 * deg) / (k * 0.25 * deg);
  const cases = [
    { name: 'east–west crests', h: (lat) => 1000 + amplitude * Math.cos(k * lat), band: 50, mu: amplitude / Math.SQRT2, slope: () => amplitude * k / R * sampling / Math.SQRT2, gamma: [0, 0.1], theta: Math.PI / 2 },
    { name: 'north–south crests', h: (lat, lon) => 1000 + amplitude * Math.cos(k * lon), band: 12, mu: amplitude / Math.SQRT2, slope: (lat) => amplitude * k / (R * Math.cos(lat)) * sampling / Math.SQRT2, gamma: [0, 0.1], theta: 0 },
    { name: 'oblique crests', h: (lat, lon) => 1000 + amplitude * Math.cos(k * (lat + lon)), band: 12, mu: amplitude / Math.SQRT2, slope: () => amplitude * k * Math.SQRT2 / R * sampling / Math.SQRT2, gamma: [0, 0.15], theta: Math.PI / 4 },
    { name: 'egg crate', h: (lat, lon) => 1000 + amplitude * Math.cos(k * lat) * Math.cos(k * lon), band: 12, mu: amplitude / 2, slope: () => amplitude * k / R * sampling / 2, gamma: [0.8, 1], theta: null },
  ];
  for (const c of cases) {
    const { mesh, fields } = ridgeFields(c.h);
    let n = 0;
    for (let i = 0; i < mesh.nCells; i++) {
      const lat = mesh.latCell[i];
      if (Math.abs(lat) > c.band * deg) continue;
      n++;
      assert.ok(fields.count[i] >= 16, `${c.name}: cell ${i} holds ${fields.count[i]} points`);
      close(fields.deviation[i], c.mu, 0.1, `${c.name} μ at ${(lat / deg).toFixed(1)}°`);
      close(fields.slope[i], c.slope(lat), 0.08, `${c.name} σ at ${(lat / deg).toFixed(1)}°`);
      assert.ok(fields.anisotropy[i] >= c.gamma[0] && fields.anisotropy[i] <= c.gamma[1], `${c.name} γ ${fields.anisotropy[i]}`);
      if (c.theta !== null) assert.ok(Math.abs(Math.sin(fields.orientation[i] - c.theta)) < 0.1, `${c.name} θ ${fields.orientation[i]} against ${c.theta}`);
    }
    assert.ok(n > 20, `${c.name}: ${n} cells tested`);
  }
  const flat = ridgeFields(() => 700).fields;
  for (let i = 0; i < flat.deviation.length; i++) assert.ok(flat.deviation[i] < 1e-6 && flat.slope[i] < 1e-9, 'a plateau has no subgrid orography');
});
