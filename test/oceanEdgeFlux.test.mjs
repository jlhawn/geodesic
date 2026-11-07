import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createOcean, SHEAR_DIFFUSIVITY, HYCOM_DIFFUSION_SPEED, shearResolution } from '../js/ocean/layered.module.js';
import { edgeLayerFluxes } from '../scripts/oceanEdgeFlux.mjs';

test('the heat the edge fluxes carry across a latitude is the heat that leaves the cells south of it in a step', () => {
  const CELSIUS = 273.15, mesh = buildMesh(new Grid(16)), { nCells: C, nEdges: E, cellsOnEdge, latCell, areaCell } = mesh, dt = 1350;
  const ocean = createOcean(mesh, { everySteps: 1, eddyDiffusivity: 0 }), L = ocean.layers;
  const surfaceT = Float64Array.from(latCell, (lat) => 300 - 26 * Math.sin(lat) ** 2), ice = new Float64Array(C), oceanFlux = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  const stress = new Float64Array(E);
  for (let e = 0; e < E; e++) {
    const x = mesh.xEdge[3 * e], y = mesh.xEdge[3 * e + 1], r = Math.hypot(x, y);
    stress[e] = r > 0 ? -0.1 * Math.cos(3 * mesh.latEdge[e]) * (-y * mesh.nEdge[3 * e] + x * mesh.nEdge[3 * e + 1]) / r : 0;
  }
  const step = () => { surfaceT.set(ocean.T0); ocean.advance(surfaceT, ice, oceanFlux, stress, dt); };
  for (let n = 0; n < 2 * 86400 / dt; n++) step();
  const fluxes = edgeLayerFluxes(mesh, { layers: L, fixed: ocean.table.fixed, options: { tracerAdvection: 'third', shearDiffusivity: SHEAR_DIFFUSIVITY, shearDiffusivityCap: 4000, lateralDiffusionSpeed: HYCOM_DIFFUSION_SPEED }, cellOcean: ocean.cellOcean, depth: ocean.D, sill: ocean.sillEdge });
  const out = { volume: new Float64Array(L), net: new Float64Array(L), heat: new Float64Array(L), diffusion: new Float64Array(L), donor: new Int32Array(L) };
  const caps = [-30, -10, -5, 5, 10, 30].map((degrees) => {
    const phi = degrees * Math.PI / 180, edges = [];
    for (let e = 0; e < E; e++) { const south = latCell[cellsOnEdge[2 * e]] < phi; if (south !== latCell[cellsOnEdge[2 * e + 1]] < phi) edges.push([e, south ? 1 : -1]); }
    return { degrees, phi, edges };
  });
  const snapshot = () => {
    const h = Float64Array.from(ocean.h), T = Float64Array.from(ocean.Q, (q, n) => q / h[n]);
    fluxes.load({ h, u: Float64Array.from(ocean.u), T, eta: Float64Array.from(ocean.eta) });
    return caps.map(({ phi, edges }) => {
      let held = 0, water = 0, model = 0, donor = 0, volume = 0, column = 0, columnHeat = 0;
      for (let i = 0; i < C; i++) if (latCell[i] < phi) for (let k = 0; k < L; k++) { const n = k * C + i; held += areaCell[i] * (ocean.Q[n] - CELSIUS * ocean.h[n]); water += areaCell[i] * ocean.h[n]; }
      for (const [e, sign] of edges) {
        fluxes.edge(e, out);
        for (let k = 0; k < L; k++) {
          model += sign * (out.heat[k] - CELSIUS * out.volume[k] + out.diffusion[k]); donor += sign * out.volume[k] * (T[k * C + out.donor[k]] - CELSIUS); volume += sign * out.volume[k];
          const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1], thick = 0.5 * (h[k * C + a] + h[k * C + b]);
          column += thick; columnHeat += thick * (0.5 * (T[k * C + a] + T[k * C + b]) - CELSIUS);
        }
      }
      return { held, water, model, donor, volume, mean: columnHeat / column };
    });
  };
  const total = () => { let q = 0; for (let n = 0; n < L * C; n++) q += areaCell[n % C] * (ocean.Q[n] - CELSIUS * ocean.h[n]); return q; };
  const areaSouth = caps.map(({ phi }) => { let a = 0, all = 0; for (let i = 0; i < C; i++) { all += areaCell[i]; if (latCell[i] < phi) a += areaCell[i]; } return a / all; });
  const steps = 12 * 3600 / dt, start = total(), first = snapshot(), sums = caps.map(() => ({ model: 0, donor: 0, volume: 0, mean: 0 }));
  let previous = first;
  for (let n = 0; n < steps; n++) {
    step();
    const now = snapshot();
    caps.forEach((_, c) => { for (const name of ['model', 'donor', 'volume', 'mean']) sums[c][name] += 0.5 * (previous[c][name] + now[c][name]) / steps; });
    previous = now;
  }
  const global = (total() - start) / (steps * dt);
  console.log(`the whole sea's heat changes by ${(global / 1e9).toFixed(4)}e9 K·m³/s through no edge`);
  caps.forEach(({ degrees }, c) => {
    const span = steps * dt, actual = (previous[c].held - first[c].held) / span, filled = (previous[c].water - first[c].water) / span;
    const missing = filled + sums[c].volume, model = -sums[c].model + missing * sums[c].mean + areaSouth[c] * global, donor = -sums[c].donor + missing * sums[c].mean + areaSouth[c] * global;
    console.log(`${degrees}°: water south rises ${(filled / 1e6).toFixed(3)} Sv against ${(-sums[c].volume / 1e6).toFixed(3)} Sv through the edges; heat south rises ${(actual / 1e9).toFixed(4)}e9 K·m³/s, the model's edge fluxes ${(model / 1e9).toFixed(4)}e9, the donor rule ${(donor / 1e9).toFixed(4)}e9`);
    assert.ok(Math.abs(model - actual) < 0.01 * Math.abs(actual) + 1e7, `${degrees}°: ${model} against ${actual}`);
  });
});

test('under shearDiffusivityResolution the edge fluxes scale the shear-scaled diffusion of the fixed-depth layers by the resolution function, and nothing else', () => {
  const mesh = buildMesh(new Grid(16)), { nCells: C, nEdges: E, latCell, latEdge } = mesh, L = 3, fixed = 2, resolve = { scale: 4, power: 2 };
  const h = new Float64Array(L * C), u = new Float64Array(L * E), T = new Float64Array(L * C);
  for (let i = 0; i < C; i++) for (let k = 0; k < L; k++) { h[k * C + i] = [10, 20, 3970][k]; T[k * C + i] = 300 - 20 * Math.sin(latCell[i]) ** 2 - 2 * k; }
  for (let k = 0; k < L; k++) for (let e = 0; e < E; e++) u[k * E + e] = 0.4 * Math.exp(-k) * Math.cos(latEdge[e]) * Math.sin(3 * e);
  const options = { tracerAdvection: 'upwind', shearDiffusivity: SHEAR_DIFFUSIVITY, shearDiffusivityCap: 4000, lateralDiffusionSpeed: 0 };
  const shape = { layers: L, fixed, cellOcean: new Uint8Array(C).fill(1), depth: new Float64Array(C).fill(4000) };
  const whole = edgeLayerFluxes(mesh, { ...shape, options }), scaled = edgeLayerFluxes(mesh, { ...shape, options: { ...options, shearDiffusivityResolution: resolve } });
  const state = { h, u, T, eta: new Float64Array(C) }, F = shearResolution(mesh, resolve);
  whole.load(state); scaled.load(state);
  const out = () => ({ volume: new Float64Array(L), net: new Float64Array(L), heat: new Float64Array(L), diffusion: new Float64Array(L), donor: new Int32Array(L) });
  const a = out(), b = out();
  let moved = 0;
  for (let e = 0; e < E; e++) {
    whole.edge(e, a); scaled.edge(e, b);
    for (let k = 0; k < L; k++) {
      assert.equal(b.volume[k], a.volume[k]); assert.equal(b.heat[k], a.heat[k]);
      const expected = k < fixed ? F[e] * a.diffusion[k] : 0;
      assert.ok(Math.abs(b.diffusion[k] - expected) <= 1e-12 * Math.abs(a.diffusion[k]), `edge ${e} layer ${k}: ${b.diffusion[k]} against ${expected}`);
      if (k < fixed && F[e] < 0.9 && a.diffusion[k] !== 0) moved++;
    }
  }
  assert.ok(moved > 0, 'no edge carries tapered diffusion');
});
