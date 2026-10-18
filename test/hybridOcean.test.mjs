import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createOcean, beyondCells, limitedEdgeValue } from '../js/ocean/layered.module.js';
import { createGeography, syntheticTopography } from '../js/geography.module.js';
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

test('a brackish sea keeps its fixed-depth surface layers, hands back its first layer\'s temperature, and keeps its heat once stepping', () => {
  const { mesh, C, ocean, surfaceT, ice, flux, calm } = hybridOcean(8, { salinityProfile: () => 12 }, () => 200, 283);
  const atMinimum = (tolerance, when) => { for (let i = 0; i < C; i++) for (let k = 0; k < HYBRID_LAYERS.fixed; k++) assert.ok(Math.abs(ocean.h[k * C + i] - HYBRID_LAYERS.minimum[k]) < tolerance, `surface layer ${k} of column ${i} ${when}: ${ocean.h[k * C + i]} m`); };
  atMinimum(1e-9, 'at the start');
  for (let i = 0; i < C; i++) assert.ok(Math.abs(surfaceT[i] - ocean.Q[i] / ocean.h[i]) < 0.1, `the surface temperature handed back is within the analytic profile's first-layer mean at ${i}`);
  ocean.advance(surfaceT, ice, flux, calm, DT);
  const before = totalHeatSalt(ocean, mesh);
  for (let n = 0; n < 4; n++) ocean.advance(surfaceT, ice, flux, calm, DT);
  const after = totalHeatSalt(ocean, mesh);
  assert.ok(Math.abs(after.heat - before.heat) < 1e-9 * Math.abs(before.heat), `heat kept: ${before.heat} → ${after.heat}`);
  atMinimum(1e-6, 'after five steps');
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

test('the beyond cell of each edge lies on the far side of its cell along the edge, joined to it by the recorded edge', () => {
  const mesh = meshOf(8), { nEdges: E, cellsOnEdge, xCell } = mesh, { cell, edge } = beyondCells(mesh);
  let worst = 1;
  for (let e = 0; e < E; e++) {
    const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
    const d = [0, 1, 2].map((c) => xCell[3 * b + c] - xCell[3 * a + c]);
    for (const side of [0, 1]) {
      const from = side ? b : a, j = cell[2 * e + side], ee = edge[2 * e + side];
      assert.ok(j >= 0 && j !== a && j !== b, `edge ${e} side ${side} has a beyond cell`);
      assert.ok((cellsOnEdge[2 * ee] === from && cellsOnEdge[2 * ee + 1] === j) || (cellsOnEdge[2 * ee + 1] === from && cellsOnEdge[2 * ee] === j), 'the beyond edge joins the cell to its beyond cell');
      const v = [0, 1, 2].map((c) => xCell[3 * j + c] - xCell[3 * from + c]);
      const cos = (side ? 1 : -1) * (v[0] * d[0] + v[1] * d[1] + v[2] * d[2]) / Math.hypot(...v) / Math.hypot(...d);
      worst = Math.min(worst, cos);
    }
  }
  console.log(`beyond cells at N=8: the worst alignment with the edge's line is ${worst.toFixed(3)}`);
  assert.ok(worst > 0.5, `every beyond cell lies within 60° of the edge's line, worst cosine ${worst}`);
});

test('the limited edge value is the midpoint on a linear field, the third-order upwind value on a smooth monotone one, bounded by its cells, and the upwind value at an extremum', () => {
  const at = (values) => limitedEdgeValue(Float64Array.from(values), 1, 2, 0, 1, 1);
  assert.ok(Math.abs(at([-1, 0, 1]) - 0.5) < 1e-15, 'linear');
  const smooth = [Math.exp(-0.1), 1, Math.exp(0.1)], quick = (2 * smooth[2] + 5 * smooth[1] - smooth[0]) / 6;
  assert.ok(Math.abs(at(smooth) - quick) < 1e-15, `third-order upwind ${quick} against ${at(smooth)}`);
  assert.equal(at([1, 0, 1]), 0, 'a minimum takes the upwind value');
  assert.equal(at([0, 1, 0]), 1, 'a maximum takes the upwind value');
  assert.equal(at([5, 5, 5]), 5, 'uniform');
  for (const [uu, u, d] of [[0, 1, 10], [0, 9, 10], [-50, 1, 2], [3, 1, 2], [1, 1, 2]]) {
    const v = at([uu, u, d]);
    assert.ok(v >= Math.min(u, d) - 1e-15 && v <= Math.max(u, d) + 1e-15, `${[uu, u, d]} gives ${v} within its cells`);
  }
  assert.ok(Math.abs(limitedEdgeValue(Float64Array.from([-2, 0, 1]), 1, 2, 0, 2, 1) - 0.5) < 1e-15, 'a linear field over a doubled upwind spacing is still the midpoint');
});

test('on a smooth field under solid-body rotation the third-order tracer tendency is far closer to the exact one than upwind, and conserves heat', () => {
  const run = (tracerAdvection) => {
    const { mesh, C, E, ocean } = hybridOcean(8, { tracerAdvection });
    const { xCell, xEdge, nEdge, areaCell, radius } = mesh, omega = 1e-7, h = ocean.h, L = ocean.layers;
    for (let e = 0; e < E; e++) {
      const r = Math.hypot(xEdge[3 * e], xEdge[3 * e + 1], xEdge[3 * e + 2]), y = xEdge[3 * e + 1] / r, z = xEdge[3 * e + 2] / r;
      const ue = omega * radius * (-z * nEdge[3 * e + 1] + y * nEdge[3 * e + 2]);
      for (let k = 0; k < L; k++) ocean.u[k * E + e] = ue;
    }
    const field = (i) => { const r = Math.hypot(xCell[3 * i], xCell[3 * i + 1], xCell[3 * i + 2]); return 290 + 5 * (xCell[3 * i + 2] / r) ** 2; };
    for (let i = 0; i < C; i++) for (let k = 0; k < L; k++) ocean.Q[k * C + i] = h[k * C + i] * field(i);
    ocean.prepare(ocean.state);
    const [dh, , dQ] = ocean.stages[0];
    ocean.tendencyLayers(ocean.state, ocean.stages[0], 0, 1, 'tracers');
    let err = 0, ref = 0, net = 0, gross = 0;
    for (let i = 0; i < C; i++) {
      const r = Math.hypot(xCell[3 * i], xCell[3 * i + 1], xCell[3 * i + 2]), y = xCell[3 * i + 1] / r, z = xCell[3 * i + 2] / r;
      const exact = -h[i] * 10 * omega * y * z, got = dQ[i] - field(i) * dh[i];
      err += areaCell[i] * (got - exact) ** 2; ref += areaCell[i] * exact ** 2; net += areaCell[i] * dQ[i]; gross += areaCell[i] * Math.abs(dQ[i]);
    }
    return { error: Math.sqrt(err / ref), net: Math.abs(net) / gross };
  };
  const upwind = run('upwind'), third = run('third');
  console.log(`tracer tendency on a smooth field at N=8: rms error ${(100 * third.error).toFixed(1)} % third-order against ${(100 * upwind.error).toFixed(1)} % upwind; net heat tendency ${third.net.toExponential(1)} of the gross`);
  assert.ok(third.net < 1e-12 && upwind.net < 1e-12, 'the tendency moves heat without making it');
  assert.ok(third.error < 0.5 * upwind.error, `the third-order error ${third.error} is well under the upwind error ${upwind.error}`);
  assert.ok(third.error < 0.25, 'and small in itself');
});

test('the shear-scaled diffusivity moves heat down the gradient in the fixed-depth layers only, under the cap, and makes none, and the background diffusion reaches every layer', () => {
  const run = (shearDiffusivity, lateralDiffusionSpeed = 0) => {
    const { mesh, C, E, ocean } = hybridOcean(8, { shearDiffusivity, lateralDiffusionSpeed, tracerAdvection: 'upwind' });
    const { xCell, xEdge, nEdge, areaCell, radius, cellsOnCell, nEdgesOnCell, maxEdges } = mesh, omega = 1e-7, h = ocean.h, L = ocean.layers;
    for (let e = 0; e < E; e++) {
      const r = Math.hypot(xEdge[3 * e], xEdge[3 * e + 1], xEdge[3 * e + 2]), y = xEdge[3 * e + 1] / r, z = xEdge[3 * e + 2] / r;
      for (let k = 0; k < L; k++) ocean.u[k * E + e] = omega * radius * (-z * nEdge[3 * e + 1] + y * nEdge[3 * e + 2]);
    }
    const field = (i) => { const r = Math.hypot(xCell[3 * i], xCell[3 * i + 1], xCell[3 * i + 2]); return 290 + 5 * (xCell[3 * i + 2] / r) ** 2; };
    for (let i = 0; i < C; i++) for (let k = 0; k < L; k++) ocean.Q[k * C + i] = h[k * C + i] * field(i);
    ocean.prepare(ocean.state);
    ocean.tendencyLayers(ocean.state, ocean.stages[0], 0, L, 'tracers');
    const dQ = Float64Array.from(ocean.stages[0][2]);
    return { dQ, C, L, areaCell, field, cellsOnCell, nEdgesOnCell, maxEdges, fixed: ocean.table.fixed, h };
  };
  const off = run(0), on = run(0.2), { C, L, areaCell, field, cellsOnCell, nEdgesOnCell, maxEdges, fixed } = on;
  let net = 0, gross = 0, agree = 0, counted = 0, deepChanged = 0;
  for (let i = 0; i < C; i++) {
    const d = on.dQ[i] - off.dQ[i];
    net += areaCell[i] * d; gross += areaCell[i] * Math.abs(d);
    let nbr = 0; for (let m = 0; m < nEdgesOnCell[i]; m++) nbr += field(cellsOnCell[maxEdges * i + m]); nbr = nbr / nEdgesOnCell[i] - field(i);
    if (Math.abs(nbr) > 1e-3) { counted++; if (Math.sign(d) === Math.sign(nbr)) agree++; }
    for (let k = fixed; k < L; k++) if (on.dQ[k * C + i] !== off.dQ[k * C + i]) deepChanged++;
  }
  const maxRate = Math.max(...Array.from({ length: C }, (_, i) => Math.abs(on.dQ[i] - off.dQ[i]) / on.h[i]));
  console.log(`shear diffusivity at N=8: down-gradient on ${agree} of ${counted} cells, net ${(net / gross).toExponential(1)} of the gross, the largest surface rate ${(maxRate * 86400).toExponential(2)} K/day, no change in ${L - fixed} deep layers`);
  assert.ok(Math.abs(net) < 1e-12 * gross, 'conserves heat');
  assert.ok(agree > 0.95 * counted, 'down the gradient');
  assert.equal(deepChanged, 0, 'the layers below the fixed-depth ones are untouched');
  assert.ok(gross > 0, 'it acts');
  assert.ok(maxRate * 86400 < 0.2, `under the cap the surface warms or cools by less than 0.2 K/day on a 5 K field at N=8, not ${maxRate * 86400}`);
  const back = run(0, 0.005);
  let backNet = 0, backGross = 0, deepAgree = 0, deepCounted = 0;
  for (let i = 0; i < C; i++) {
    let nbr = 0; for (let m = 0; m < nEdgesOnCell[i]; m++) nbr += field(cellsOnCell[maxEdges * i + m]); nbr = nbr / nEdgesOnCell[i] - field(i);
    for (let k = 0; k < L; k++) { const d = back.dQ[k * C + i] - off.dQ[k * C + i]; backNet += areaCell[i] * d; backGross += areaCell[i] * Math.abs(d); if (k >= fixed && back.h[k * C + i] > 1 && Math.abs(nbr) > 1e-3) { deepCounted++; if (Math.sign(d) === Math.sign(nbr)) deepAgree++; } }
  }
  console.log(`background diffusion at 5 mm/s: down-gradient on ${deepAgree} of ${deepCounted} deep cell-layers, net ${(backNet / backGross).toExponential(1)} of the gross`);
  assert.ok(Math.abs(backNet) < 1e-12 * backGross && deepAgree > 0.95 * deepCounted && deepCounted > 0, 'the background diffusion is conservative and down-gradient in the deep layers');
});

test('a class held at its floor on one side of an edge, with no depth in common with the layer across it, follows the layer above while the front between two stratifications adjusts', () => {
  const mesh = meshOf(8), C = mesh.nCells, E = mesh.nEdges, deep = (i) => Math.cos(mesh.lonCell[i]) > 0;
  const bathymetry = Float64Array.from({ length: C }, (_, i) => (deep(i) ? 4000 : 1500));
  const ocean = createOcean(mesh, { ...FLAT, bathymetry });
  const surfaceT = Float64Array.from({ length: C }, (_, i) => (deep(i) ? 300 : 280)), ice = new Float64Array(C), flux = new Float64Array(C), calm = new Float64Array(E);
  ocean.initialize(surfaceT, ice);
  const L = ocean.layers, { cellsOnEdge } = mesh, tokenEdge = Array.from(HYBRID_LAYERS.minimum, (m) => 2 * Math.min(m, HYBRID_LAYERS.spacing));
  const tokens = () => {
    const Ha = new Float64Array(E), Hb = new Float64Array(E), list = [];
    for (let k = 0; k < L; k++) for (let e = 0; e < E; e++) {
      const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1], ha = ocean.h[k * C + a], hb = ocean.h[k * C + b], he = Math.min(ha, hb);
      if (k >= HYBRID_LAYERS.fixed && ocean.edgeOcean[e] && he > 2 * H_MIN && he <= tokenEdge[k] && Math.min(Ha[e] + ha, Hb[e] + hb) <= Math.max(Ha[e], Hb[e])) list.push([k, e]);
      Ha[e] += ha; Hb[e] += hb;
    }
    return list;
  };
  let shear = 0, count = Infinity;
  for (let n = 0; n < 64; n++) {
    ocean.advance(surfaceT, ice, flux, calm, DT);
    const list = tokens();
    count = Math.min(count, list.length);
    for (const [k, e] of list) shear = Math.max(shear, Math.abs(ocean.u[k * E + e] - ocean.u[(k - 1) * E + e]));
  }
  assert.ok(count > 500, `${count} token edge-layers`);
  assert.ok(maxSpeed(ocean) > 0.5, `the front adjusts: ${maxSpeed(ocean)} m/s`);
  assert.ok(shear < 1e-2, `a token's shear to the layer above reached ${shear} m/s`);
  assert.equal(ocean.diagnostics().oceanLimited, 0);
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

test('an edge takes the deepest channel its two cells hold between them as its sill, the hybrid ocean cuts the edge column there so the classes below carry nothing across, and the slab ocean is untouched', () => {
  const mesh = meshOf(8), C = mesh.nCells, E = mesh.nEdges, { cellsOnEdge, lonCell, latCell } = mesh;
  const wrap = (x) => Math.atan2(Math.sin(x), Math.cos(x));
  const ridge = (lon) => Math.abs(lon) < 1.5 * DEG, strip = (lon) => Math.abs(wrap(lon - Math.PI)) < 1.5 * DEG, within = (lat) => Math.abs(lat) < 40 * DEG;
  const topography = syntheticTopography(180, 360, (lat, lon) => (within(lat) && ridge(lon) ? -300 : within(lat) && strip(lon) ? 100 : -4000));
  const geography = createGeography(mesh, topography);
  const across = (e, lon0) => { const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1], la = wrap(lonCell[a] - lon0), lb = wrap(lonCell[b] - lon0); return geography.edgeOcean[e] && la * lb < 0 && Math.abs(la) > 3.3 * DEG && Math.abs(lb) > 3.3 * DEG && Math.abs(la) < 15 * DEG && Math.abs(lb) < 15 * DEG && Math.abs(latCell[a]) < 30 * DEG && Math.abs(latCell[b]) < 30 * DEG; };
  const ridgeEdges = [], stripEdges = [], openEdges = [];
  for (let e = 0; e < E; e++) {
    if (across(e, 0)) ridgeEdges.push(e);
    else if (across(e, Math.PI)) stripEdges.push(e);
    else if (geography.edgeOcean[e] && Math.abs(lonCell[cellsOnEdge[2 * e]]) > 20 * DEG && Math.abs(lonCell[cellsOnEdge[2 * e]]) < 70 * DEG) openEdges.push(e);
  }
  assert.ok(ridgeEdges.length >= 4 && stripEdges.length >= 4 && openEdges.length > 20, `${ridgeEdges.length} ridge, ${stripEdges.length} strip, ${openEdges.length} open edges`);
  for (const e of ridgeEdges) assert.equal(geography.edgeSill[e], 300, `ridge edge ${e}`);
  for (const e of stripEdges) assert.equal(geography.edgeSill[e], 0, `land-strip edge ${e}`);
  for (const e of openEdges) assert.equal(geography.edgeSill[e], 4000, `open edge ${e}`);
  const slab = createOcean(mesh, { everySteps: 1, geography }), hybrid = createOcean(mesh, { ...FLAT, geography });
  const floor = (e) => Math.min(hybrid.D[cellsOnEdge[2 * e]], hybrid.D[cellsOnEdge[2 * e + 1]]);
  for (const e of ridgeEdges) { assert.equal(hybrid.sillEdge[e], 300); assert.ok(slab.sillEdge[e] > 1000 && slab.sillEdge[e] === floor(e), `the slab keeps min(D) at ${e}: ${slab.sillEdge[e]}`); }
  for (const e of stripEdges) { assert.equal(hybrid.sillEdge[e], 50); assert.equal(slab.sillEdge[e], floor(e)); }
  for (const e of openEdges) { assert.equal(hybrid.sillEdge[e], floor(e)); assert.equal(slab.sillEdge[e], floor(e)); }
  const front = Float64Array.from({ length: C }, (_, i) => (lonCell[i] < 0 ? 280 : 300)), ice = new Float64Array(C), flux = new Float64Array(C), calm = new Float64Array(E);
  const L = hybrid.layers;
  const transport = (ocean) => { let through = 0, below = 0; for (const e of ridgeEdges) { let z = 0; for (let k = 0; k < L; k++) { const t = Math.abs(ocean.hEdge[k * E + e] * ocean.u[k * E + e]) * mesh.dvEdge[e]; through += t; if (z >= 300) below += t; z += ocean.hEdge[k * E + e]; } } return { through, below }; };
  const open = createOcean(mesh, { ...FLAT, geography, edgeSill: null });
  hybrid.initialize(front, ice); open.initialize(front, ice);
  let sillThrough = 0, sillBelow = 0, openThrough = 0;
  for (let n = 0; n < 48; n++) {
    hybrid.advance(front, ice, flux, calm, DT); open.advance(front, ice, flux, calm, DT);
    const a = transport(hybrid), b = transport(open);
    sillThrough += a.through; sillBelow += a.below; openThrough += b.through;
  }
  for (const e of ridgeEdges) {
    let z = 0, sum = 0;
    for (let k = 0; k < L; k++) { const he = hybrid.hEdge[k * E + e]; if (z >= 300) assert.equal(he, 0, `layer ${k} of ridge edge ${e} at ${z} m`); z += Math.min(hybrid.h[k * C + cellsOnEdge[2 * e]], hybrid.h[k * C + cellsOnEdge[2 * e + 1]]); sum += he; }
    assert.ok(Math.abs(sum - 300) < 1, `ridge edge ${e} column ${sum} m`);
  }
  assert.ok(maxSpeed(hybrid) > 0.1 && maxSpeed(open) > 0.1, 'the front adjusts');
  assert.equal(sillBelow, 0);
  assert.ok(openThrough > 0 && sillThrough < 0.5 * openThrough, `transport through the ridge edges over 48 steps: ${(sillThrough / 48 / 1e6).toFixed(1)} Sv with the sill, ${(openThrough / 48 / 1e6).toFixed(1)} Sv without`);
  assert.equal(hybrid.diagnostics().oceanLimited, 0);
});

test('the slab options that have no hybrid meaning are refused', () => {
  const mesh = meshOf(8), C = mesh.nCells;
  assert.throws(() => createOcean(mesh, { ...FLAT, bathymetry: new Float64Array(C).fill(4000), ekmanBlock: true }), /slab/);
  assert.throws(() => createOcean(mesh, { ...FLAT, bathymetry: new Float64Array(C).fill(4000), vertical: 'sigma' }), /vertical/);
  assert.throws(() => createOcean(mesh, { ...FLAT, bathymetry: new Float64Array(C).fill(4000), shearDiffusivity: 0 }), /shearDiffusivity/);
  assert.throws(() => createOcean(mesh, { bathymetry: new Float64Array(C).fill(4000), tracerAdvection: 'third' }), /hybrid/);
  assert.throws(() => createOcean(mesh, { bathymetry: new Float64Array(C).fill(4000), shearDiffusivity: 0.2 }), /hybrid/);
  assert.throws(() => createOcean(mesh, { bathymetry: new Float64Array(C).fill(4000), lateralDiffusionSpeed: 0.005 }), /hybrid/);
  createOcean(mesh, { ...FLAT, bathymetry: new Float64Array(C).fill(4000), tracerAdvection: 'upwind', shearDiffusivity: 0 });
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
