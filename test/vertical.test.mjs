import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { cellVector } from '../js/dynamics/operators.module.js';
import { levelFields, verticalVelocity } from '../js/levels.module.js';
import { createOcean as createLayeredOcean, depthFields } from '../js/ocean/layered.module.js';

const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 300 + 2500 * Math.exp(-(((lat - 0.3) / 0.3) ** 2)) : -4000));

test('the vertical velocity at a level rests on the same πσ̇ the core diagnoses, and vanishes at the top and the ground', () => {
  const model = createModel(new Grid(6), { topography });
  const { mesh, core, state } = model;
  const C = mesh.nCells, E = mesh.nEdges, K = core.K;
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) state[a].set(init[a]);
  for (let n = 0; n < 6; n++) model.step(900);
  const [pi, theta, u, , q] = state;
  core.diagnose(pi, theta, q, state[5]);
  core.tendency(state, state.map((x) => new Float64Array(x.length)));
  const own = core.arrays.piSigmaDot;
  const layerWind = (k) => cellVector(mesh, u.subarray(k * E, (k + 1) * E));
  const piSigmaDot = new Float64Array((K + 1) * C);
  const { temperature } = levelFields(core, pi, theta, layerWind, 500, q);
  const w = verticalVelocity(mesh, core, pi, u, 500, temperature, new Float32Array(C), piSigmaDot);
  let worst = 0, scale = 0;
  for (let n = 0; n < piSigmaDot.length; n++) { worst = Math.max(worst, Math.abs(piSigmaDot[n] - own[n])); scale = Math.max(scale, Math.abs(own[n])); }
  assert.ok(scale > 1e-4, `the flow has vertical motion: ${scale}`);
  assert.ok(worst <= 1e-9 * scale, `πσ̇ differs from the core's by ${worst} against ${scale}`);
  for (let i = 0; i < C; i++) { assert.equal(piSigmaDot[i], 0); assert.equal(piSigmaDot[K * C + i], 0); }
  let rms = 0, top = 0;
  for (let i = 0; i < C; i++) { rms += w[i] * w[i]; top = Math.max(top, Math.abs(w[i])); }
  rms = Math.sqrt(rms / C);
  assert.ok(rms > 1e-5 && top < 5, `500 hPa vertical velocity rms ${rms} m/s, largest ${top}`);
  const lowest = verticalVelocity(mesh, core, pi, u, 'surface', levelFields(core, pi, theta, layerWind, 'surface', q).temperature);
  const underground = verticalVelocity(mesh, core, pi, u, 1000, levelFields(core, pi, theta, layerWind, 1000, q).temperature);
  for (let i = 0; i < C; i++) if (pi[i] < 1000e2) assert.ok(Math.abs(underground[i] - lowest[i]) < 1e-12, `a level under the ground shows the lowest layer at ${i}`);
});

test('the fields at a depth pick the layer holding it, mask the sea floor, and integrate the transport above it to a divergence-free whole', () => {
  const model = createModel(new Grid(5), { topography });
  const { mesh } = model, C = mesh.nCells, E = mesh.nEdges;
  const ocean = createLayeredOcean(mesh, { geography: model.geography });
  const init = initializeState(model, {});
  ocean.initialize(init[3], init[6]);
  const L = ocean.layers, { h, u, cellOcean } = ocean;
  for (let e = 0; e < E; e++) u[e] = 0.3 * Math.sin(3 * mesh.latEdge[e]) * Math.cos(2 * Math.atan2(mesh.xEdge[3 * e + 1], mesh.xEdge[3 * e]));
  for (let k = 1; k < L; k++) for (let e = 0; e < E; e++) u[k * E + e] = u[e] / (k + 1);
  for (let e = 0; e < E; e++) if (!ocean.edgeOcean[e]) for (let k = 0; k < L; k++) u[k * E + e] = 0;
  const temperature = (k, i) => 300 - 3 * k - 1e-3 * i;
  const surface = depthFields(mesh, L, { h, u, temperature, cellOcean }, 0);
  const top = cellVector(mesh, u.subarray(0, E));
  for (let i = 0; i < C; i++) {
    if (!cellOcean[i]) { assert.ok(Number.isNaN(surface.temperature[i]) && Number.isNaN(surface.upwelling[i])); continue; }
    assert.equal(surface.temperature[i], Math.fround(temperature(0, i)));
    assert.equal(surface.upwelling[i], 0, 'nothing lies above the surface');
    for (let c = 0; c < 3; c++) assert.ok(Math.abs(surface.current[3 * i + c] - top[3 * i + c]) < 1e-6);
  }
  const deep = depthFields(mesh, L, { h, u, temperature, cellOcean }, 250);
  let wet = 0, dry = 0, sum = 0, magnitude = 0;
  for (let i = 0; i < C; i++) {
    if (!cellOcean[i]) continue;
    let above = 0, layer = -1;
    for (let k = 0; k < L; k++) { if (250 < above + h[k * C + i]) { layer = k; break; } above += h[k * C + i]; }
    if (layer < 0) { dry++; assert.ok(Number.isNaN(deep.temperature[i]) && Number.isNaN(deep.upwelling[i]) && deep.current[3 * i] === 0); continue; }
    wet++;
    assert.ok(layer >= 1, `250 m lies below the mixed layer at ${i}`);
    assert.equal(deep.temperature[i], Math.fround(temperature(layer, i)));
    sum += mesh.areaCell[i] * deep.upwelling[i]; magnitude += mesh.areaCell[i] * Math.abs(deep.upwelling[i]);
  }
  assert.ok(wet > 0.9 * (wet + dry), `most ocean columns reach 250 m: ${wet} of ${wet + dry}`);
  assert.ok(magnitude > 0);
  assert.ok(Math.abs(sum) < 1e-6 * magnitude, `the upwelling through 250 m sums to ${sum} against ${magnitude} in magnitude`);
  const abyss = depthFields(mesh, L, { h, u, temperature, cellOcean }, 20000);
  for (let i = 0; i < C; i++) assert.ok(Number.isNaN(abyss.temperature[i]));
});
