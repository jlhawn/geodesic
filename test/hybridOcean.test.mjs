import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createOcean, beyondCells, limitedEdgeValue, inertialStepLimit, EDGE_TRANSPORTS } from '../js/ocean/layered.module.js';
import { HYBRID_LAYERS, H_MIN } from '../js/ocean/column.module.js';
import { DT, DEG, FLAT, mesh, C, E, hybridOcean, columnSums, maxSpeed, volumeOf, zonalWindOnEdges, totalHeatSalt, uniformOcean } from './helpers/layered.mjs';

test('the hybrid start puts every column on the depth table clipped to its floor, sums to the depth, keeps the minimum thickness, and leaves a uniform stratification exactly at rest', () => {
  const { mesh, C, ocean, surfaceT, ice, flux, calm } = hybridOcean();
  const sums = columnSums(ocean);
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
  const { mesh, C, ocean, surfaceT, ice, flux, calm } = hybridOcean({ salinityProfile: () => 12 }, () => 200, 283);
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

test('through dynamics, regrid and column physics the hybrid ocean conserves water, heat and salt to rounding, with no dried layer, whichever edge transport it keeps', () => {
  for (const edgeTransport of EDGE_TRANSPORTS) {
    const { mesh, C, ocean, surfaceT, ice, flux, calm } = hybridOcean({ edgeTransport });
    const wind = zonalWindOnEdges(mesh, (lat) => 0.1 * Math.cos(2 * lat));
    ocean.advance(surfaceT, ice, flux, wind, DT);
    const before = totalHeatSalt(ocean, mesh), volume = volumeOf(ocean);
    for (let n = 0; n < 24; n++) ocean.advance(surfaceT, ice, flux, wind, DT);
    const after = totalHeatSalt(ocean, mesh);
    assert.ok(Math.abs(after.heat - before.heat) < 1e-10 * before.heat, `${edgeTransport}: heat ${before.heat} → ${after.heat}`);
    assert.ok(Math.abs(after.salt - before.salt) < 1e-10 * before.salt, `${edgeTransport}: salt ${before.salt} → ${after.salt}`);
    assert.ok(Math.abs(volumeOf(ocean) - volume) < 1e-10 * volume, `${edgeTransport}: volume`);
    const d = ocean.diagnostics();
    assert.equal(d.oceanDried, 0);
    assert.ok(d.oceanSpeed > 0.01 && d.oceanSpeed < 1, `${edgeTransport}: speed ${d.oceanSpeed}`);
    assert.ok(d.oceanUpperDepth > 5 && d.oceanUpperDepth < 100, `${edgeTransport}: mixed layer ${d.oceanUpperDepth} m`);
  }
});

const tangentialMisfit = (ocean) => {
  const reconstructed = new Float64Array(E);
  let misfit = 0, speed = 0;
  for (let k = 0; k < ocean.layers; k++) {
    ocean.reconstruct(ocean.u, k * E, reconstructed);
    for (let e = 0; e < E; e++) { if (!ocean.edgeOcean[e]) continue; misfit += (ocean.v[k * E + e] - reconstructed[e]) ** 2; speed += ocean.u[k * E + e] ** 2; }
  }
  return Math.sqrt(misfit / speed);
};

test('with coriolis cd the hybrid ocean conserves water, heat and salt to rounding as with TRiSK, v carries none of them and stays near the reconstruction of u', () => {
  const { mesh, ocean, surfaceT, ice, flux } = hybridOcean({ coriolis: 'cd' });
  assert.equal(ocean.state.length, 5);
  assert.equal(hybridOcean().ocean.state.length, 4);
  const wind = zonalWindOnEdges(mesh, (lat) => 0.1 * Math.cos(2 * lat));
  ocean.advance(surfaceT, ice, flux, wind, DT);
  const before = totalHeatSalt(ocean, mesh), volume = volumeOf(ocean);
  for (let n = 0; n < 24; n++) ocean.advance(surfaceT, ice, flux, wind, DT);
  const after = totalHeatSalt(ocean, mesh);
  assert.ok(Math.abs(after.heat - before.heat) < 1e-10 * before.heat, `heat ${before.heat} → ${after.heat}`);
  assert.ok(Math.abs(after.salt - before.salt) < 1e-10 * before.salt, `salt ${before.salt} → ${after.salt}`);
  assert.ok(Math.abs(volumeOf(ocean) - volume) < 1e-10 * volume, 'volume');
  const d = ocean.diagnostics();
  assert.equal(d.oceanDried, 0);
  assert.ok(d.oceanSpeed > 0.01 && d.oceanSpeed < 1, `speed ${d.oceanSpeed}`);
  for (const x of ocean.v) assert.ok(Number.isFinite(x) && Math.abs(x) <= 5, 'v is finite and within the speed limit');
  const misfit = tangentialMisfit(ocean);
  console.log(`coriolis cd after 25 steps: rms(v − R(u)) / rms(u) = ${misfit.toExponential(2)}`);
  assert.ok(misfit < 0.1, `v is ${misfit} of u away from R(u)`);
});

test('coriolis cd is recorded in the saved ocean, an unknown scheme or a relaxation time that is not positive is refused, and a load starts v as the reconstruction of u', () => {
  const { mesh, ocean, surfaceT, ice, flux } = hybridOcean({ coriolis: 'cd', cdTau: 1e5 });
  const wind = zonalWindOnEdges(mesh, (lat) => 0.1 * Math.cos(2 * lat));
  for (let n = 0; n < 6; n++) ocean.advance(surfaceT, ice, flux, wind, DT);
  const saved = ocean.serialize();
  assert.equal(saved.coriolis, 'cd');
  assert.equal(saved.cdTau, 1e5);
  assert.equal(hybridOcean().ocean.serialize().coriolis, undefined);
  assert.equal(hybridOcean().ocean.v, null);
  const twin = hybridOcean({ coriolis: 'cd', cdTau: 1e5 });
  twin.ocean.load(saved, twin.surfaceT, twin.ice);
  assert.equal(tangentialMisfit(twin.ocean), 0);
  assert.ok(Math.abs(tangentialMisfit(ocean) - tangentialMisfit(twin.ocean)) > 0, 'a running v is not exactly R(u)');
  assert.throws(() => hybridOcean({ coriolis: 'd' }), /coriolis is one of/);
  assert.throws(() => hybridOcean({ coriolis: 'cd', cdTau: 0 }), /cdTau/);
});

test('a step longer than a radian of the fastest inertial rotation takes its dynamics in equal parts that keep water, heat and salt to rounding', () => {
  const { mesh, ocean, surfaceT, ice, flux } = hybridOcean();
  const wind = zonalWindOnEdges(mesh, (lat) => 0.1 * Math.cos(2 * lat)), step = 4 * 3600;
  assert.equal(Math.ceil(step / inertialStepLimit(mesh)), 3);
  ocean.advance(surfaceT, ice, flux, wind, step);
  const before = totalHeatSalt(ocean, mesh), volume = volumeOf(ocean);
  for (let n = 0; n < 6; n++) ocean.advance(surfaceT, ice, flux, wind, step);
  const after = totalHeatSalt(ocean, mesh);
  assert.ok(Math.abs(after.heat - before.heat) < 1e-10 * before.heat, `heat ${before.heat} → ${after.heat}`);
  assert.ok(Math.abs(after.salt - before.salt) < 1e-10 * before.salt, `salt ${before.salt} → ${after.salt}`);
  assert.ok(Math.abs(volumeOf(ocean) - volume) < 1e-10 * volume, 'volume');
  assert.equal(ocean.diagnostics().oceanDried, 0);
});

test('the beyond cell of each edge lies on the far side of its cell along the edge, joined to it by the recorded edge', () => {
  const { cellsOnEdge, xCell } = mesh, { cell, edge } = beyondCells(mesh);
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
  console.log(`beyond cells at N=16: the worst alignment with the edge's line is ${worst.toFixed(3)}`);
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
    const { mesh, C, E, ocean } = hybridOcean({ tracerAdvection });
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
  console.log(`tracer tendency on a smooth field at N=16: rms error ${(100 * third.error).toFixed(1)} % third-order against ${(100 * upwind.error).toFixed(1)} % upwind; net heat tendency ${third.net.toExponential(1)} of the gross`);
  assert.ok(third.net < 1e-12 && upwind.net < 1e-12, 'the tendency moves heat without making it');
  assert.ok(third.error < 0.5 * upwind.error, `the third-order error ${third.error} is well under the upwind error ${upwind.error}`);
  assert.ok(third.error < 0.05, `and small in itself: ${third.error}`);
});

test('the shear-scaled diffusivity moves heat down the gradient in the fixed-depth layers only, under the cap, and makes none, and the background diffusion reaches every layer', () => {
  const run = (shearDiffusivity, lateralDiffusionSpeed = 0) => {
    const { mesh, C, E, ocean } = hybridOcean({ shearDiffusivity, lateralDiffusionSpeed, tracerAdvection: 'upwind' });
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
  console.log(`shear diffusivity at N=16: down-gradient on ${agree} of ${counted} cells, net ${(net / gross).toExponential(1)} of the gross, the largest surface rate ${(maxRate * 86400).toExponential(2)} K/day, no change in ${L - fixed} deep layers`);
  assert.ok(Math.abs(net) < 1e-12 * gross, 'conserves heat');
  assert.ok(agree > 0.95 * counted, 'down the gradient');
  assert.equal(deepChanged, 0, 'the layers below the fixed-depth ones are untouched');
  assert.ok(gross > 0, 'it acts');
  assert.ok(maxRate * 86400 < 1e-3, `under the cap the surface warms or cools by less than 0.001 K/day on a 5 K field at N=16, not ${maxRate * 86400}`);
  const back = run(0, 0.005);
  let backNet = 0, backGross = 0, deepAgree = 0, deepCounted = 0;
  for (let i = 0; i < C; i++) {
    let nbr = 0; for (let m = 0; m < nEdgesOnCell[i]; m++) nbr += field(cellsOnCell[maxEdges * i + m]); nbr = nbr / nEdgesOnCell[i] - field(i);
    for (let k = 0; k < L; k++) { const d = back.dQ[k * C + i] - off.dQ[k * C + i]; backNet += areaCell[i] * d; backGross += areaCell[i] * Math.abs(d); if (k >= fixed && back.h[k * C + i] > 1 && Math.abs(nbr) > 1e-3) { deepCounted++; if (Math.sign(d) === Math.sign(nbr)) deepAgree++; } }
  }
  console.log(`background diffusion at 5 mm/s: down-gradient on ${deepAgree} of ${deepCounted} deep cell-layers, net ${(backNet / backGross).toExponential(1)} of the gross`);
  assert.ok(Math.abs(backNet) < 1e-12 * backGross && deepAgree > 0.95 * deepCounted && deepCounted > 0, 'the background diffusion is conservative and down-gradient in the deep layers');
});

test('a slab ocean converts into the hybrid layers with every column\'s heat and salt kept and the barotropic transport within rounding', () => {
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
  const { mesh, C, ocean } = hybridOcean();
  for (let i = 0; i < C; i++) ocean.Q[C + i] = ocean.h[C + i] * 280;
  const top = ocean.h[C], surface = ocean.fields(0), within = ocean.fields(0.5 * (ocean.h[0] + ocean.h[C]));
  assert.ok(top < 5, `the second layer is ${top} m thick`);
  assert.ok(Math.abs(surface.temperature[0] - ocean.T0[0]) < 1e-3 && Math.abs(within.temperature[0] - 280) < 1e-9, `${surface.temperature[0]} at the surface against T0 ${ocean.T0[0]}, ${within.temperature[0]} within the second layer`);
  assert.equal(surface.h1, ocean.mixedDepths);
  assert.ok(surface.h1.every((d, i) => d > 0 && d <= 4000.5) && Number.isFinite(surface.thermoclineDepth[0]), `mixed layer ${surface.h1[0]} m, thermocline ${surface.thermoclineDepth[0]} m`);
  assert.equal(mesh.nCells, C);
});

test('the z* and isopycnal generators run the same ocean and conserve it, and a 50 m shelf keeps its layers', () => {
  for (const generator of ['zstar', 'isopycnal']) {
    const { mesh, C, ocean, surfaceT, ice, flux } = hybridOcean({ generator });
    const wind = zonalWindOnEdges(mesh, (lat) => 0.1 * Math.cos(2 * lat));
    ocean.advance(surfaceT, ice, flux, wind, DT);
    const before = totalHeatSalt(ocean, mesh);
    for (let n = 0; n < 12; n++) ocean.advance(surfaceT, ice, flux, wind, DT);
    const after = totalHeatSalt(ocean, mesh);
    assert.ok(Math.abs(after.heat - before.heat) < 1e-10 * before.heat, `${generator}: heat`);
    assert.ok(Math.abs(after.salt - before.salt) < 1e-10 * before.salt, `${generator}: salt`);
    assert.ok(ocean.u.every(Number.isFinite), `${generator}: finite`);
    const sums = columnSums(ocean);
    for (let i = 0; i < C; i++) assert.ok(Math.abs(sums[i] - ocean.D[i] - ocean.eta[i]) < 1e-9, `${generator}: column ${i}`);
  }
  const shelf = hybridOcean({}, (lat) => (Math.abs(lat) < 20 * DEG ? 50 : 4000));
  const wind = zonalWindOnEdges(shelf.mesh, (lat) => 0.1 * Math.cos(2 * lat));
  for (let n = 0; n < 12; n++) shelf.ocean.advance(shelf.surfaceT, shelf.ice, shelf.flux, wind, DT);
  assert.ok(shelf.ocean.u.every(Number.isFinite) && shelf.ocean.Q.every(Number.isFinite), 'the shelf runs');
  const sums = columnSums(shelf.ocean);
  for (let i = 0; i < shelf.C; i++) assert.ok(Math.abs(sums[i] - shelf.ocean.D[i] - shelf.ocean.eta[i]) < 1e-9, `shelf column ${i}`);
  assert.equal(shelf.ocean.diagnostics().oceanDried, 0);
});

test('the slab\'s coordinate and the options without a hybrid meaning are refused', () => {
  assert.throws(() => createOcean(mesh, { ...FLAT, bathymetry: new Float64Array(C).fill(4000), vertical: 'slab' }), /vertical/);
  assert.throws(() => createOcean(mesh, { ...FLAT, bathymetry: new Float64Array(C).fill(4000), vertical: 'sigma' }), /vertical/);
  assert.throws(() => createOcean(mesh, { ...FLAT, bathymetry: new Float64Array(C).fill(4000), shearDiffusivity: 0 }), /shearDiffusivity/);
  assert.throws(() => createOcean(mesh, { ...FLAT, bathymetry: new Float64Array(C).fill(4000), tracerAdvection: 'central' }), /tracerAdvection/);
  createOcean(mesh, { ...FLAT, bathymetry: new Float64Array(C).fill(4000), tracerAdvection: 'upwind', shearDiffusivity: 0 });
});

test('the worker engine computes the hybrid ocean bit for bit as the serial one, tendencies, regrid and column physics alike', async () => {
  const { createParallelModel } = await import('../js/parallel.module.js');
  const options = { vertical: 'hybrid', everySteps: 1, thermoclineTilt: 0, salinityProfile: () => 34.5 };
  const parallel = await createParallelModel(new Grid(16), { ocean: options }, 3), mesh = parallel.mesh, C = mesh.nCells, E = mesh.nEdges;
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
