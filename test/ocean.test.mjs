import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { cellVector } from '../js/dynamics/operators.module.js';
import { createOcean } from '../js/ocean/reducedGravity.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';

const mesh = buildMesh(new Grid(6));
const { nCells: C, nEdges: E } = mesh;
const RHO_AIR = 1.2, DRAG = 1.5e-3, RHO = 1025;

function zonalWindOnEdges(speedAt) {
  const u = new Float64Array(E);
  for (let e = 0; e < E; e++) {
    const x = mesh.xEdge[3 * e], y = mesh.xEdge[3 * e + 1], r = Math.hypot(x, y);
    const east = r > 0 ? [-y / r, x / r, 0] : [0, 0, 0];
    u[e] = speedAt(mesh.latEdge[e]) * (east[0] * mesh.nEdge[3 * e] + east[1] * mesh.nEdge[3 * e + 1] + east[2] * mesh.nEdge[3 * e + 2]);
  }
  return u;
}

function northwardTransport(ocean) {
  const vector = cellVector(mesh, ocean.u1);
  const out = new Float64Array(C);
  for (let i = 0; i < C; i++) {
    const x = mesh.xCell[3 * i], y = mesh.xCell[3 * i + 1], z = mesh.xCell[3 * i + 2], r = Math.hypot(x, y);
    const north = [-z * x / r, -z * y / r, r];
    out[i] = ocean.h1[i] * (vector[3 * i] * north[0] + vector[3 * i + 1] * north[1] + vector[3 * i + 2] * north[2]);
  }
  return out;
}

test('the ocean at rest under no wind stays at rest, bit for bit', () => {
  const ocean = createOcean(mesh, { everySteps: 1 });
  const surfaceT = new Float64Array(C).fill(290), ice = new Float64Array(C), flux = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  const calm = new Float64Array(E);
  for (let n = 0; n < 20; n++) ocean.advance(surfaceT, ice, flux, calm, 3600);
  assert.ok(ocean.u1.every((x) => x === 0) && ocean.u2.every((x) => x === 0));
  assert.ok(ocean.h1.every((x) => x === 50) && ocean.h2.every((x) => x === 350));
  assert.ok(surfaceT.every((x) => x === 290) && flux.every((x) => x === 0));
});

test('westerlies drive an equatorward Ekman transport of τ/(ρf) in both hemispheres', () => {
  const ocean = createOcean(mesh, { everySteps: 1 });
  const surfaceT = new Float64Array(C).fill(290), ice = new Float64Array(C), flux = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  const speed = 10;
  const wind = zonalWindOnEdges((lat) => speed * Math.exp(-(((Math.abs(lat) * 180 / Math.PI - 45) / 10) ** 2)));
  const stressField = Float64Array.from(wind, (w) => RHO_AIR * DRAG * speed * w);
  for (let n = 0; n < 24 * 12; n++) ocean.advance(surfaceT, ice, flux, stressField, 3600);
  const transport = northwardTransport(ocean);
  let north = 0, south = 0, an = 0, as = 0;
  for (let i = 0; i < C; i++) {
    const lat = mesh.latCell[i] * 180 / Math.PI;
    if (lat > 42 && lat < 48) { north += mesh.areaCell[i] * transport[i]; an += mesh.areaCell[i]; }
    if (lat < -42 && lat > -48) { south += mesh.areaCell[i] * transport[i]; as += mesh.areaCell[i]; }
  }
  north /= an; south /= as;
  const stress = RHO_AIR * DRAG * speed * speed;
  const f = 2 * mesh.omega * Math.sin(45 * Math.PI / 180);
  const ekman = stress / (RHO * f);
  assert.ok(north < 0 && Math.abs(-north - ekman) < 0.25 * ekman, `NH transport ${north} vs Ekman ${-ekman} m²/s`);
  assert.ok(south > 0 && Math.abs(south - ekman) < 0.25 * ekman, `SH transport ${south} vs Ekman ${ekman} m²/s`);
  console.log(`Ekman transport at 45°: ${(-north).toFixed(3)} (N) ${south.toFixed(3)} (S) vs τ/ρf ${ekman.toFixed(3)} m²/s; max |u1| ${Math.max(...ocean.u1.map(Math.abs)).toFixed(3)} m/s`);
});

test('wind-driven flow moves heat around but conserves it', () => {
  const ocean = createOcean(mesh, { everySteps: 1 });
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => 300 - 25 * Math.sin(lat) ** 2), ice = new Float64Array(C), flux = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  const wind = zonalWindOnEdges((lat) => 8 * Math.cos(3 * lat));
  const stressField = Float64Array.from(wind, (w) => RHO_AIR * DRAG * 8 * w);
  const heat = () => { let h = 0; for (let i = 0; i < C; i++) h += mesh.areaCell[i] * (ocean.H1[i] + ocean.H2[i]); return h; };
  const before = heat(), sstBefore = Float64Array.from(surfaceT);
  for (let n = 0; n < 24 * 5; n++) ocean.advance(surfaceT, ice, flux, stressField, 3600);
  assert.ok(Math.abs(heat() - before) < 1e-11 * before, `heat ${before} → ${heat()}`);
  let moved = 0;
  for (let i = 0; i < C; i++) moved = Math.max(moved, Math.abs(surfaceT[i] - sstBefore[i]));
  assert.ok(moved > 1e-3, `the surface temperature changed by only ${moved} K`);
});

test('heat converged under ice goes to the ice base as a flux and the water stays at the freezing point', () => {
  const ocean = createOcean(mesh, { everySteps: 1, diffusivity: 2 });
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => Math.abs(lat) > 1.2 ? 260 : 295);
  const ice = Float64Array.from(mesh.latCell, (lat) => (Math.abs(lat) > 1.2 ? 1 : 0)), flux = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  for (let n = 0; n < 6; n++) ocean.advance(surfaceT, ice, flux, new Float64Array(E), 3600);
  let edgeFlux = 0, interior = 0;
  for (let i = 0; i < C; i++) {
    if (ice[i] > 0) {
      assert.equal(ocean.T1[i], 271.35);
      let openNeighbour = false;
      for (let k = 0; k < mesh.nEdgesOnCell[i]; k++) if (ice[mesh.cellsOnCell[mesh.maxEdges * i + k]] === 0) openNeighbour = true;
      if (openNeighbour) edgeFlux = Math.max(edgeFlux, flux[i]); else interior = Math.max(interior, Math.abs(flux[i]));
    } else assert.equal(flux[i], 0);
  }
  assert.ok(edgeFlux > 1, `ice at the edge receives ${edgeFlux} W/m²`);
  assert.ok(interior < 0.1 * edgeFlux, `ice interior receives ${interior} W/m² against ${edgeFlux} at the edge`);
});

test('the coupled model steps with the ocean and stays bounded', () => {
  const m = createModel(new Grid(4));
  const init = initializeState(m, {});
  for (let a = 0; a < init.length; a++) m.state[a].set(init[a]);
  m.ocean.initialize(m.state[3], m.state[6]);
  for (let n = 0; n < 48 * 3; n++) m.step(1800);
  const d = m.diagnostics();
  assert.ok(Number.isFinite(d.maxWind) && d.maxWind < 80);
  assert.ok(d.iceFraction >= 0 && d.iceFraction <= 1);
  assert.ok(d.oceanUpperDepth > 20 && d.oceanUpperDepth < 100, `upper layer ${d.oceanUpperDepth} m`);
  assert.ok(d.oceanSpeed > 0 && d.oceanSpeed < 2, `max current ${d.oceanSpeed} m/s`);
  assert.ok(Number.isFinite(d.oceanHeat) && d.oceanThermoclineT > 270 && d.oceanThermoclineT < 300, `thermocline ${d.oceanThermoclineT} K`);
  for (let i = 0; i < m.mesh.nCells; i++) assert.ok(m.state[3][i] > 200 && m.state[3][i] < 320, `surface temperature ${m.state[3][i]} at cell ${i} (ice ${m.state[6][i]} m)`);
  console.log(`N=4 three days with the ocean: h1 ${d.oceanUpperDepth.toFixed(2)} m, max current ${d.oceanSpeed.toFixed(3)} m/s, thermocline ${d.oceanThermoclineT.toFixed(1)} K, ocean heat ${(d.oceanHeat / 1e9).toFixed(2)} GJ/m²`);
});
