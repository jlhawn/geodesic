import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createOcean } from '../js/ocean/layered.module.js';
import { HYBRID_LAYERS, H_MIN } from '../js/ocean/column.module.js';
import { seawaterDensity } from '../js/ocean/seawater.module.js';
import { zonalWindOnEdges, totalHeatSalt, uniformOcean } from './helpers/layered.mjs';

/*
 * The hybrid ocean (vertical: 'hybrid') on Grid(8) and Grid(16): flat
 * floors, an analytic start with no tilt, 34.5 psu everywhere, driven
 * through advance with the ocean's own surface temperature handed back.
 */
const DT = 1350, DEG = Math.PI / 180;
const meshes = new Map();
const meshOf = (N) => { if (!meshes.has(N)) meshes.set(N, buildMesh(new Grid(N))); return meshes.get(N); };
const FLAT = { vertical: 'hybrid', everySteps: 1, thermoclineTilt: 0, salinityProfile: () => 34.5 };
function hybridOcean(N, options = {}, depth = () => 4000, surface = 290) {
  const mesh = meshOf(N), C = mesh.nCells, E = mesh.nEdges;
  const bathymetry = Float64Array.from({ length: C }, (_, i) => depth(mesh.latCell[i], mesh.lonCell[i], i));
  const ocean = createOcean(mesh, { ...FLAT, bathymetry, ...options });
  const surfaceT = new Float64Array(C).fill(surface), ice = new Float64Array(C), flux = new Float64Array(C), calm = new Float64Array(E);
  ocean.initialize(surfaceT, ice);
  return { mesh, C, E, ocean, surfaceT, ice, flux, calm };
}
const columnSums = (ocean, C) => Float64Array.from({ length: C }, (_, i) => { let s = 0; for (let k = 0; k < ocean.layers; k++) s += ocean.h[k * C + i]; return s; });
const maxSpeed = (ocean) => { let m = 0; for (let n = 0; n < ocean.u.length; n++) m = Math.max(m, Math.abs(ocean.u[n])); return m; };
const volumeOf = (ocean, mesh) => { let v = 0; const C = mesh.nCells; for (let i = 0; i < C; i++) if (ocean.cellOcean[i]) for (let k = 0; k < ocean.layers; k++) v += mesh.areaCell[i] * ocean.h[k * C + i]; return v; };

test('the hybrid start puts every column on the depth table clipped to its floor, sums to the depth, keeps the minimum thickness, and leaves a uniform stratification exactly at rest', () => {
  const { mesh, C, ocean, surfaceT, ice, flux, calm } = hybridOcean(8);
  const sums = columnSums(ocean, C);
  for (let i = 0; i < C; i++) assert.ok(Math.abs(sums[i] - ocean.D[i] - ocean.eta[i]) < 1e-9, `column ${i} sums to ${sums[i]}`);
  let thinnest = Infinity;
  for (let n = 0; n < ocean.h.length; n++) thinnest = Math.min(thinnest, ocean.h[n]);
  assert.ok(thinnest >= H_MIN - 1e-12, `thinnest layer ${thinnest} m`);
  assert.equal(ocean.layers, HYBRID_LAYERS.K);
  for (let k = 0; k < 12; k++) assert.ok(Math.abs(ocean.h[k * C] - HYBRID_LAYERS.minimum[k]) < 1e-9, `surface layer ${k} at its minimum`);
  for (let n = 0; n < 10; n++) ocean.advance(surfaceT, ice, flux, calm, DT);
  assert.ok(maxSpeed(ocean) < 1e-12, `at rest: ${maxSpeed(ocean)} m/s`);
  assert.ok(Math.max(...ocean.eta.map(Math.abs)) < 1e-9, 'the free surface stays level');
});

test('through dynamics, regrid and column physics the hybrid ocean conserves water, heat and salt to rounding, with no dried layer', () => {
  const { mesh, C, ocean, surfaceT, ice, flux, calm } = hybridOcean(8);
  const wind = zonalWindOnEdges(mesh, (lat) => 0.1 * Math.cos(2 * lat));
  ocean.advance(surfaceT, ice, flux, wind, DT);
  const before = totalHeatSalt(ocean, mesh), volume = volumeOf(ocean, mesh);
  for (let n = 0; n < 24; n++) ocean.advance(surfaceT, ice, flux, wind, DT);
  const after = totalHeatSalt(ocean, mesh);
  assert.ok(Math.abs(after.heat - before.heat) < 1e-10 * before.heat, `heat ${before.heat} → ${after.heat}`);
  assert.ok(Math.abs(after.salt - before.salt) < 1e-10 * before.salt, `salt ${before.salt} → ${after.salt}`);
  assert.ok(Math.abs(volumeOf(ocean, mesh) - volume) < 1e-10 * volume, 'volume');
  const d = ocean.diagnostics();
  assert.equal(d.oceanDried, 0);
  assert.ok(d.oceanSpeed > 0.01 && d.oceanSpeed < 1, `speed ${d.oceanSpeed}`);
  assert.ok(d.oceanUpperDepth > 5 && d.oceanUpperDepth < 100, `mixed layer ${d.oceanUpperDepth} m`);
});

test('tropical easterlies drive equatorial upwelling and a poleward surface drift, and the mixed layer deepens under the stronger wind', () => {
  const { mesh, C, ocean, surfaceT, ice, flux } = hybridOcean(16);
  const wind = zonalWindOnEdges(mesh, (lat) => -0.08 * Math.exp(-((lat / (15 * DEG)) ** 2)));
  for (let n = 0; n < 5 * 64; n++) ocean.advance(surfaceT, ice, flux, wind, DT);
  const f = ocean.fields(50);
  let equator = 0, nEq = 0, side = 0, nSide = 0, mixedEq = 0, mixedPole = 0, nPole = 0;
  for (let i = 0; i < C; i++) {
    const lat = Math.abs(mesh.latCell[i]);
    if (lat < 3 * DEG) { equator += f.upwelling[i]; nEq++; mixedEq += f.h1[i]; }
    else if (lat > 6 * DEG && lat < 12 * DEG) { side += f.upwelling[i]; nSide++; }
    else if (lat > 60 * DEG) { mixedPole += f.h1[i]; nPole++; }
  }
  console.log(`after five days of easterlies: equatorial upwelling ${(equator / nEq * 86400).toFixed(2)} m/day, 6–12° ${(side / nSide * 86400).toFixed(2)} m/day; mixed layer at the equator ${(mixedEq / nEq).toFixed(1)} m, poleward of 60° ${(mixedPole / nPole).toFixed(1)} m`);
  assert.ok(equator / nEq > 0, 'upwelling at the equator');
  assert.ok(side / nSide < equator / nEq, 'less off the equator');
  assert.ok(mixedEq / nEq > mixedPole / nPole, 'the wind mixes deeper');
});

test('a slab ocean converts into the hybrid layers with every column\'s heat and salt kept and the barotropic transport within rounding', () => {
  const mesh = meshOf(8), C = mesh.nCells, E = mesh.nEdges;
  const slab = uniformOcean({ bathymetry: new Float64Array(C).fill(3000) }, mesh), wind = zonalWindOnEdges(mesh, (lat) => 0.1 * Math.cos(2 * lat));
  for (let n = 0; n < 12; n++) slab.ocean.advance(slab.surfaceT, slab.ice, slab.flux, wind, DT);
  const saved = slab.ocean.serialize(), Ls = slab.ocean.layers;
  const hybrid = createOcean(mesh, { ...FLAT, bathymetry: new Float64Array(C).fill(3000) });
  hybrid.load(saved, Float64Array.from(slab.surfaceT), slab.ice);
  const L = hybrid.layers;
  for (let i = 0; i < C; i++) {
    let hs = 0, qs = 0, ws = 0, hh = 0, qh = 0, wh = 0;
    for (let k = 0; k < Ls; k++) { const n = k * C + i; hs += saved.h[n]; qs += saved.h[n] * saved.T[n]; ws += saved.h[n] * saved.S[n]; }
    for (let k = 0; k < L; k++) { const n = k * C + i; hh += hybrid.h[n]; qh += hybrid.Q[n]; wh += hybrid.W[n]; }
    assert.ok(Math.abs(hh - hs) < 1e-9 * hs, `column ${i} water ${hs} → ${hh}`);
    assert.ok(Math.abs(qh - qs) < 1e-9 * qs, `column ${i} heat ${qs} → ${qh}`);
    assert.ok(Math.abs(wh - ws) < 1e-9 * ws, `column ${i} salt ${ws} → ${wh}`);
  }
  let worst = 0;
  for (let e = 0; e < E; e++) {
    const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1];
    let was = 0, is = 0, scale = 0;
    for (let k = 0; k < Ls; k++) { const he = 0.5 * (saved.h[k * C + a] + saved.h[k * C + b]); was += he * saved.u[k * E + e]; scale += he * Math.abs(saved.u[k * E + e]); }
    for (let k = 0; k < L; k++) { const he = 0.5 * (hybrid.h[k * C + a] + hybrid.h[k * C + b]); is += he * hybrid.u[k * E + e]; }
    if (scale > 0) worst = Math.max(worst, Math.abs(is - was) / scale);
  }
  assert.ok(worst < 1e-9, `transport changed by ${worst} of its scale`);
  assert.ok(hybrid.u.every(Number.isFinite));
  const serialized = hybrid.serialize();
  assert.equal(serialized.vertical, 'hybrid');
  assert.equal(serialized.densities.length, L);
  const again = createOcean(mesh, { ...FLAT, bathymetry: new Float64Array(C).fill(3000) });
  again.load(serialized, Float64Array.from(slab.surfaceT), slab.ice);
  for (let n = 0; n < L * C; n++) assert.ok(Math.abs(again.h[n] - hybrid.h[n]) < 1e-9, `thickness ${n} reloads`);
  for (let n = 0; n < 10; n++) { hybrid.advance(slab.surfaceT, slab.ice, slab.flux, wind, DT); again.advance(Float64Array.from(slab.surfaceT), slab.ice, slab.flux, wind, DT); }
  assert.ok(hybrid.u.every(Number.isFinite) && hybrid.Q.every(Number.isFinite));
});

test('the fields at a depth take the hybrid layer that holds it, however thin, with the mixed-layer depth as h1', () => {
  const { mesh, C, ocean } = hybridOcean(8);
  for (let i = 0; i < C; i++) ocean.Q[C + i] = ocean.h[C + i] * 280;
  const top = ocean.h[C], surface = ocean.fields(0), within = ocean.fields(0.5 * (ocean.h[0] + ocean.h[C]));
  assert.ok(top < 5, `the second layer is ${top} m thick`);
  assert.ok(Math.abs(surface.temperature[0] - ocean.T0[0]) < 1e-3 && Math.abs(within.temperature[0] - 280) < 1e-9, `${surface.temperature[0]} at the surface against T0 ${ocean.T0[0]}, ${within.temperature[0]} within the second layer`);
  assert.equal(surface.h1, ocean.mixedDepths);
  assert.ok(surface.h1.every((d, i) => d > 0 && d <= 4000.5) && Number.isFinite(surface.thermoclineDepth[0]), `mixed layer ${surface.h1[0]} m, thermocline ${surface.thermoclineDepth[0]} m`);
  assert.equal(mesh.nCells, C);
});

test('over a seamount the pressure force at the depths the layers share leaves a uniform stratification nearly at rest', () => {
  const { mesh, C, ocean, surfaceT, ice, flux, calm } = hybridOcean(16, {}, (lat, lon) => 4000 - 3500 * Math.exp(-((lat / (8 * DEG)) ** 2 + (lon / (8 * DEG)) ** 2)));
  for (let n = 0; n < 64; n++) ocean.advance(surfaceT, ice, flux, calm, DT);
  const speed = maxSpeed(ocean);
  console.log(`seamount to 500 m under a uniform stratification: fastest current after a day ${speed.toExponential(2)} m/s`);
  assert.ok(speed < 0.01, `spurious current ${speed} m/s`);
  assert.equal(ocean.diagnostics().oceanDried, 0);
});

test('the z* and isopycnal generators run the same ocean and conserve it, and a 50 m shelf keeps its layers', () => {
  for (const generator of ['zstar', 'isopycnal']) {
    const { mesh, C, ocean, surfaceT, ice, flux } = hybridOcean(8, { generator });
    const wind = zonalWindOnEdges(mesh, (lat) => 0.1 * Math.cos(2 * lat));
    ocean.advance(surfaceT, ice, flux, wind, DT);
    const before = totalHeatSalt(ocean, mesh);
    for (let n = 0; n < 12; n++) ocean.advance(surfaceT, ice, flux, wind, DT);
    const after = totalHeatSalt(ocean, mesh);
    assert.ok(Math.abs(after.heat - before.heat) < 1e-10 * before.heat, `${generator}: heat`);
    assert.ok(Math.abs(after.salt - before.salt) < 1e-10 * before.salt, `${generator}: salt`);
    assert.ok(ocean.u.every(Number.isFinite), `${generator}: finite`);
    const sums = columnSums(ocean, C);
    for (let i = 0; i < C; i++) assert.ok(Math.abs(sums[i] - ocean.D[i] - ocean.eta[i]) < 1e-9, `${generator}: column ${i}`);
  }
  const shelf = hybridOcean(8, {}, (lat) => (Math.abs(lat) < 20 * DEG ? 50 : 4000));
  const wind = zonalWindOnEdges(shelf.mesh, (lat) => 0.1 * Math.cos(2 * lat));
  for (let n = 0; n < 12; n++) shelf.ocean.advance(shelf.surfaceT, shelf.ice, shelf.flux, wind, DT);
  assert.ok(shelf.ocean.u.every(Number.isFinite) && shelf.ocean.Q.every(Number.isFinite), 'the shelf runs');
  const sums = columnSums(shelf.ocean, shelf.C);
  for (let i = 0; i < shelf.C; i++) assert.ok(Math.abs(sums[i] - shelf.ocean.D[i] - shelf.ocean.eta[i]) < 1e-9, `shelf column ${i}`);
  assert.equal(shelf.ocean.diagnostics().oceanDried, 0);
});

test('the slab options that have no hybrid meaning are refused', () => {
  const mesh = meshOf(8), C = mesh.nCells;
  assert.throws(() => createOcean(mesh, { ...FLAT, bathymetry: new Float64Array(C).fill(4000), ekmanBlock: true }), /slab/);
  assert.throws(() => createOcean(mesh, { ...FLAT, bathymetry: new Float64Array(C).fill(4000), vertical: 'sigma' }), /vertical/);
});

test('the worker engine computes the hybrid ocean bit for bit as the serial one, tendencies, regrid and column physics alike', async () => {
  const { createParallelModel } = await import('../js/parallel.module.js');
  const options = { vertical: 'hybrid', everySteps: 1, thermoclineTilt: 0, salinityProfile: () => 34.5 };
  const parallel = await createParallelModel(new Grid(8), { ocean: options }, 3), mesh = parallel.mesh, C = mesh.nCells, E = mesh.nEdges;
  const serial = createOcean(mesh, options);
  const surfaceT = new Float64Array(C).fill(290), ice = new Float64Array(C), flux = new Float64Array(C), wind = zonalWindOnEdges(mesh, (lat) => 0.1 * Math.cos(2 * lat));
  const surfaceP = Float64Array.from(surfaceT);
  serial.initialize(surfaceT, ice);
  parallel.ocean.initialize(surfaceP, ice);
  try {
    for (let n = 0; n < 6; n++) { serial.advance(surfaceT, ice, flux, wind, DT); parallel.ocean.advance(surfaceP, ice, flux, wind, DT); }
    for (const name of ['h', 'u', 'Q', 'W', 'eta']) {
      const a = serial[name], b = parallel.ocean[name];
      let differ = 0;
      for (let n = 0; n < a.length; n++) if (a[n] !== b[n]) differ++;
      assert.equal(differ, 0, `${name}: ${differ} of ${a.length} values differ`);
    }
    assert.ok(parallel.ocean.diagnostics().oceanSpeed > 0.01, 'the ocean moved');
    assert.ok(parallel.ocean.mixedDepths.some((d) => d > 5), 'the column physics ran');
  } finally {
    await parallel.close();
  }
});
