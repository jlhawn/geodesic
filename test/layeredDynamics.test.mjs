import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { cellVector, laplacianVelocity } from '../js/dynamics/operators.module.js';
import { createOcean, closureVelocity, THIN, EPS } from '../js/ocean/layered.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { seawaterDensity } from '../js/ocean/seawater.module.js';
import { RHO_AIR, DRAG, RHO, DEG, mesh, C, E, zonalWindOnEdges, totalHeatSalt, northwardMixedTransport, slowOcean } from './helpers/layered.mjs';

test('the ocean at rest under zero stress on a bumpy bottom stays at rest and conserves heat, salt and surface temperature', () => {
  const bathymetry = new Float64Array(C);
  for (let i = 0; i < C; i++) bathymetry[i] = 4000 + 1500 * Math.sin(3 * mesh.lonCell[i]) * Math.cos(2 * mesh.latCell[i]);
  const ocean = createOcean(mesh, { everySteps: 1, bathymetry, closureHours: 12, thermoclineTilt: 0, salinityProfile: () => 35 });
  const surfaceT = new Float64Array(C).fill(290), ice = new Float64Array(C), flux = new Float64Array(C), calm = new Float64Array(E);
  ocean.initialize(surfaceT, ice);
  const surfaceBefore = Float64Array.from(surfaceT);
  const before = totalHeatSalt(ocean, mesh);
  for (let n = 0; n < 24; n++) ocean.advance(surfaceT, ice, flux, calm, 1350);
  let maxU = 0, maxEta = 0;
  for (const x of ocean.u) maxU = Math.max(maxU, Math.abs(x));
  for (const x of ocean.eta) maxEta = Math.max(maxEta, Math.abs(x));
  assert.ok(maxU < 1e-5, `|u| max ${maxU} m/s`);
  assert.ok(maxEta < 1e-4, `|eta| max ${maxEta} m`);
  const after = totalHeatSalt(ocean, mesh);
  assert.ok(Math.abs(after.heat - before.heat) < 1e-10 * Math.abs(before.heat), `heat ${before.heat} -> ${after.heat}`);
  assert.ok(Math.abs(after.salt - before.salt) < 1e-10 * Math.abs(before.salt), `salt ${before.salt} -> ${after.salt}`);
  let maxDrift = 0;
  for (let i = 0; i < C; i++) maxDrift = Math.max(maxDrift, Math.abs(surfaceT[i] - surfaceBefore[i]));
  assert.ok(maxDrift < 1e-9, `surface temperature drifted by ${maxDrift} K at rest`);
  assert.equal(ocean.diagnostics().oceanLimited, 0);
});

test('the sum of layer thicknesses minus the bathymetry equals the free surface in every ocean cell', () => {
  const ocean = createOcean(mesh, { everySteps: 1 });
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => 275 + 25 * Math.cos(lat) ** 2), ice = new Float64Array(C), flux = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  const wind = zonalWindOnEdges(mesh, (lat) => 8 * Math.cos(3 * lat));
  const stressField = Float64Array.from(wind, (w) => RHO_AIR * DRAG * 8 * w);
  for (let n = 0; n < 40; n++) ocean.advance(surfaceT, ice, flux, stressField, 1350);
  let maxMismatch = 0;
  for (let i = 0; i < C; i++) {
    if (!ocean.cellOcean[i]) continue;
    let sum = 0;
    for (let k = 0; k < ocean.layers; k++) sum += ocean.h[k * C + i];
    maxMismatch = Math.max(maxMismatch, Math.abs(sum - ocean.D[i] - ocean.eta[i]));
  }
  assert.ok(maxMismatch < 1e-9, `max |sum h - D - eta| = ${maxMismatch} m`);
  assert.equal(ocean.diagnostics().oceanLimited, 0);
});

test('westerlies at 45N and 45S drive an equatorward mixed-layer Ekman transport near tau/(rho f)', async () => {
  const ocean = await slowOcean(mesh, { ocean: { everySteps: 1 }, surfaceT: new Float64Array(C).fill(290) });
  console.log(`Ekman transport on the ${ocean.engine} ocean`);
  const speed = 10;
  const wind = zonalWindOnEdges(mesh, (lat) => speed * Math.exp(-(((Math.abs(lat) * 180 / Math.PI - 45) / 10) ** 2)));
  const stressField = Float64Array.from(wind, (w) => RHO_AIR * DRAG * speed * w);
  const dt = 1350;
  const f = 2 * mesh.omega * Math.sin(45 * DEG);
  const inertialPeriod = 2 * Math.PI / f;
  await ocean.advance(stressField, dt, 12 * 64);
  // A suddenly-imposed wind excites a ~17h inertial oscillation in the mixed
  // layer that an instantaneous snapshot aliases against daily sampling;
  // average over two inertial periods to recover the mean Ekman balance.
  const windowSteps = Math.round(2 * inertialPeriod / dt);
  let north = 0, south = 0, samples = 0;
  for (let n = 0; n < windowSteps; n++) {
    await ocean.advance(stressField, dt);
    const transport = northwardMixedTransport(await ocean.download(), mesh);
    let nSum = 0, sSum = 0, nArea = 0, sArea = 0;
    for (let i = 0; i < C; i++) {
      const lat = mesh.latCell[i] * 180 / Math.PI;
      if (lat > 42 && lat < 48) { nSum += mesh.areaCell[i] * transport[i]; nArea += mesh.areaCell[i]; }
      if (lat < -42 && lat > -48) { sSum += mesh.areaCell[i] * transport[i]; sArea += mesh.areaCell[i]; }
    }
    north += nSum / nArea; south += sSum / sArea; samples++;
  }
  north /= samples; south /= samples;
  const stress = RHO_AIR * DRAG * speed * speed;
  const ekman = stress / (RHO * f);
  console.log(`mixed-layer transport averaged over two inertial periods: ${north.toPrecision(4)} m²/s at 45°N, ${south.toPrecision(4)} m²/s at 45°S, against an Ekman ${ekman.toPrecision(4)} m²/s`);
  assert.ok(north < 0 && Math.abs(-north - ekman) < 0.25 * ekman, `NH transport ${north} vs Ekman ${-ekman} m^2/s`);
  assert.ok(south > 0 && Math.abs(south - ekman) < 0.25 * ekman, `SH transport ${south} vs Ekman ${ekman} m^2/s`);
  assert.equal((await ocean.diagnostics()).oceanLimited, 0);
  await ocean.close();
});

test('a subtropical wind over a closed 80-degree basin drives a Stommel western boundary current', async () => {
  const basin = buildMesh(new Grid(16));
  const { nCells: bC, nEdges: bE } = basin;
  const topography = syntheticTopography(180, 360, (lat, lon) => (Math.abs(lon) < 40 * DEG && lat > 12 * DEG && lat < 48 * DEG ? -4000 : 500));
  const ocean = await slowOcean(basin, { ocean: { everySteps: 1 }, topography, surfaceT: Float64Array.from(basin.latCell, (lat) => 275 + 25 * Math.cos(lat) ** 2) });
  console.log(`Stommel basin on the ${ocean.engine} ocean`);
  const wind = zonalWindOnEdges(basin, (lat) => {
    const l = lat / DEG;
    return l > 15 && l < 45 ? -Math.cos(Math.PI * (l - 15) / 30) : 0;
  });
  const stressField = Float64Array.from(wind, (w) => 0.1 * w);
  // The impulsively-started gyre sheds a near-symmetric pair of boundary
  // currents that only separates into the western-intensified Stommel
  // balance slowly: a day-by-day run shows the west/east pooled transport
  // ratio crossing 1 only after day ~90, and still a coin-flip at day 80;
  // 100 days gives a comfortable, non-marginal margin on every check below
  // while keeping the whole file's runtime well under the budget.
  const days = 100;
  await ocean.advance(stressField, 1350, days * 64);

  const { h, u } = await ocean.download();
  const Ubt = new Float64Array(bE);
  for (let e = 0; e < bE; e++) {
    let t = 0;
    const a = basin.cellsOnEdge[2 * e], b = basin.cellsOnEdge[2 * e + 1];
    for (let k = 0; k < ocean.layers; k++) t += 0.5 * (h[k * bC + a] + h[k * bC + b]) * u[k * bE + e];
    Ubt[e] = t;
  }
  const vector = cellVector(basin, Ubt);
  function band(lonLo, lonHi) {
    let v = 0, area = 0;
    for (let i = 0; i < bC; i++) {
      const lat = basin.latCell[i] / DEG, lon = basin.lonCell[i] / DEG;
      if (!ocean.cellOcean[i] || lat < 25 || lat > 35 || lon < lonLo || lon >= lonHi) continue;
      const x = basin.xCell[3 * i], y = basin.xCell[3 * i + 1], z = basin.xCell[3 * i + 2], r = Math.hypot(x, y);
      v += basin.areaCell[i] * (vector[3 * i] * (-z * x / r) + vector[3 * i + 1] * (-z * y / r) + vector[3 * i + 2] * r);
      area += basin.areaCell[i];
    }
    return area > 0 ? v / area : NaN;
  }
  const west = band(-40, -32);
  const interior = band(-20, 20);
  let maxBin = -Infinity, maxCenter = null;
  for (let center = -38; center <= 38; center += 4) {
    const v = band(center - 2, center + 2);
    if (Number.isFinite(v) && v > maxBin) { maxBin = v; maxCenter = center; }
  }
  console.log(`after ${days} days the western band carries ${west.toPrecision(4)} m²/s north, the interior ${interior.toPrecision(4)} m²/s, the strongest 4-degree band centred at ${maxCenter}°`);
  assert.ok(west > 0, `western boundary transport ${west} m^2/s should be northward`);
  assert.ok(interior < 0, `interior transport ${interior} m^2/s should be southward`);
  assert.ok(west >= 5 * Math.abs(interior), `western transport ${west} should be at least 5x the interior's ${interior}`);
  assert.ok(maxCenter === -38 || maxCenter === -34, `the largest 4-degree band is centred at ${maxCenter}, not one of the westernmost two`);
  assert.equal((await ocean.diagnostics()).oceanLimited, 0);
  await ocean.close();
});

test('wind-driven advection over 100 steps conserves total heat and salt', () => {
  const ocean = createOcean(mesh, { everySteps: 1 });
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => 300 - 25 * Math.sin(lat) ** 2), ice = new Float64Array(C), flux = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  const wind = zonalWindOnEdges(mesh, (lat) => 8 * Math.cos(3 * lat));
  const stressField = Float64Array.from(wind, (w) => RHO_AIR * DRAG * 8 * w);
  const before = totalHeatSalt(ocean, mesh);
  for (let n = 0; n < 100; n++) ocean.advance(surfaceT, ice, flux, stressField, 1350);
  const after = totalHeatSalt(ocean, mesh);
  // Under this wind the free surface swings by tens of centimetres, and the
  // per-cell rescale to the sub-stepped barotropic eta at the end of step()
  // is not exactly heat/salt-conservative (about 1e-10 relative over these
  // 100 steps).
  assert.ok(Math.abs(after.heat - before.heat) < 1e-9 * Math.abs(before.heat), `heat ${before.heat} -> ${after.heat}`);
  assert.ok(Math.abs(after.salt - before.salt) < 1e-9 * Math.abs(before.salt), `salt ${before.salt} -> ${after.salt}`);
  assert.equal(ocean.diagnostics().oceanLimited, 0);
});

test('interfacial drag, constant or from the shear, moves momentum between the layers at an edge without changing the column\'s, the classes a few metres thick included', () => {
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => 300 - 25 * Math.sin(lat) ** 2), ice = new Float64Array(C);
  const build = (options) => {
    const ocean = createOcean(mesh, { everySteps: 1, ...options });
    ocean.initialize(Float64Array.from(surfaceT), ice);
    for (let n = 0; n < ocean.u.length; n++) ocean.u[n] = 0.2 * Math.sin(0.7 * Math.floor(n / E) + 0.37 * (n % E));
    ocean.tendency(ocean.state, ocean.stages[0]);
    return ocean;
  };
  const still = build({ shearMixing: false, interfacialDrag: 0 }), du0 = still.stages[0][1];
  for (const options of [{ shearMixing: false, interfacialDrag: 5e-5 }, { shearMixing: true, interfacialDrag: 5e-5 }]) {
    const dragging = build(options), L = dragging.layers;
    const hEdge = new Float64Array(dragging.shared.hEdge), du = dragging.stages[0][1];
    let worst = 0, largest = 0, thin = 0;
    for (let e = 0; e < E; e++) {
      if (!dragging.edgeOcean[e]) continue;
      let net = 0;
      for (let k = 0; k < L; k++) {
        const n = k * E + e;
        if (k > 0 && hEdge[n] < THIN) continue;
        if (k > 0 && hEdge[n] < 50) thin++;
        const moved = hEdge[n] * (du[n] - du0[n]);
        net += moved; largest = Math.max(largest, Math.abs(moved));
      }
      worst = Math.max(worst, Math.abs(net));
    }
    console.log(`${JSON.stringify(options)}: ${thin} layer edges between ${THIN} and 50 m thick; the largest momentum moved by the drag ${largest.toExponential(2)} m²/s², the largest column imbalance ${worst.toExponential(2)}`);
    assert.ok(thin > 100, `only ${thin} thin layer edges`);
    assert.ok(largest > 1e-7 && worst < 1e-12 * largest, `the drag changes a column's momentum by ${worst} against ${largest} moved`);
  }
});

test('under shearMixing the drag between the mixed layer and the class beneath follows the Pacanowski–Philander viscosity of their Richardson number from the full velocity difference, no less than interfacialDrag', () => {
  const ocean = createOcean(mesh, { everySteps: 1, shearMixing: true, interfacialDrag: 0 }), floored = createOcean(mesh, { everySteps: 1, shearMixing: true, interfacialDrag: 5e-5 });
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => 300 - 25 * Math.sin(lat) ** 2);
  ocean.initialize(Float64Array.from(surfaceT), new Float64Array(C));
  floored.initialize(Float64Array.from(surfaceT), new Float64Array(C));
  const L = ocean.layers, rho = ocean.densities, hEdge = new Float64Array(ocean.shared.hEdge), g = 9.81, rho0 = 1025;
  let lowest = Infinity;
  const mixed = (i) => seawaterDensity(ocean.Q[i] / ocean.h[i], ocean.W[i] / ocean.h[i]);
  const errors = [], richardsons = [];
  for (const speed of [0.05, 0.1, 0.2, 0.4, 0.8]) {
    ocean.u.fill(0);
    for (let e = 0; e < E; e++) ocean.u[e] = speed * (-mesh.xEdge[3 * e + 1] * mesh.nEdge[3 * e] + mesh.xEdge[3 * e] * mesh.nEdge[3 * e + 1]);
    ocean.tendency(ocean.state, ocean.stages[0]);
    floored.u.set(ocean.u);
    floored.tendency(floored.state, floored.stages[0]);
    for (let e = 0; e < E; e++) {
      if (!ocean.edgeOcean[e] || Math.abs(mesh.latEdge[e]) > 40 * DEG) continue;
      let k = 1;
      while (k < L && hEdge[k * E + e] < THIN) k++;
      if (k === L) continue;
      const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1];
      const dz = Math.max(THIN, 0.5 * (hEdge[e] + hEdge[k * E + e])), buoyancy = Math.max(0, g * (rho[k] - 0.5 * (mixed(a) + mixed(b))) / rho0);
      const shear = speed ** 2 * (mesh.xEdge[3 * e] ** 2 + mesh.xEdge[3 * e + 1] ** 2), richardson = buoyancy * dz / shear;
      const expected = (1e-2 / (1 + 5 * richardson) ** 2 + 1e-4) / dz, rate = ocean.interfaceRate(e, 0, k);
      lowest = Math.min(lowest, floored.interfaceRate(e, 0, k) / Math.max(5e-5, rate));
      if (richardson < 2 && rate < 0.99 * 0.5 * Math.min(Math.max(hEdge[e], 20), Math.max(hEdge[k * E + e], THIN)) / 3600) { errors.push(Math.abs(rate / expected - 1)); richardsons.push(richardson); }
    }
  }
  errors.sort((p, q) => p - q); richardsons.sort((p, q) => p - q);
  const at = (list, p) => list[Math.floor(p * (list.length - 1))];
  console.log(`${errors.length} mixed-layer bases at Ri below 2 (${at(richardsons, 0.1).toFixed(2)}–${at(richardsons, 0.9).toFixed(2)}, 10th–90th percentile): the drag coefficient within ${(100 * at(errors, 0.5)).toFixed(1)}% of the Pacanowski–Philander value at the median, ${(100 * at(errors, 0.9)).toFixed(1)}% at the 90th percentile`);
  assert.ok(errors.length > 100 && at(richardsons, 0.1) < 0.3, `${errors.length} interfaces over Ri ${at(richardsons, 0.1)}–${at(richardsons, 0.9)}`);
  assert.ok(at(errors, 0.5) < 0.05 && at(errors, 0.9) < 0.2, `median ${at(errors, 0.5)}, 90th percentile ${at(errors, 0.9)}`);
  assert.ok(Math.abs(lowest - 1) < 1e-12, `with interfacialDrag 5e-5 m/s the coefficient is the larger of that and the viscosity's (${lowest})`);
});

test('closureVelocity gives the token edges beside a class a weighted share of the class\'s own flow, and the closure then pulls the class less toward the flow of the layer above them', () => {
  const m = buildMesh(new Grid(16)), mE = m.nEdges;
  const rotation = (e, sign) => sign * 0.3 * (-m.xEdge[3 * e + 1] * m.nEdge[3 * e] + m.xEdge[3 * e] * m.nEdge[3 * e + 1]);
  const inside = (i) => Math.sin(3 * m.lonCell[i]) + 0.5 * Math.cos(5 * m.latCell[i]) > 0.2;
  const hEdge = new Float64Array(mE), sea = new Uint8Array(mE).fill(1), exact = new Float64Array(mE), slaved = new Float64Array(mE);
  for (let e = 0; e < mE; e++) {
    const present = inside(m.cellsOnEdge[2 * e]) && inside(m.cellsOnEdge[2 * e + 1]);
    hEdge[e] = present ? 12 : EPS;
    exact[e] = rotation(e, 1);
    slaved[e] = present ? exact[e] : rotation(e, -1);
  }
  const beside = new Uint8Array(mE), near = new Uint8Array(mE);
  for (let e = 0; e < mE; e++) {
    if (hEdge[e] >= THIN) continue;
    for (let s = 0; s < m.nEdgesOnEdge[e]; s++) { const o = m.edgesOnEdge[m.maxEdgesOnEdge * e + s]; if (hEdge[o] >= THIN) { beside[e] = 1; near[o] = 1; } }
  }
  const pull = (field) => {
    const lap2 = laplacianVelocity(m, laplacianVelocity(m, field));
    let sum = 0, n = 0;
    for (let e = 0; e < mE; e++) if (near[e]) { sum += lap2[e] ** 2; n++; }
    return Math.sqrt(sum / n);
  };
  const whole = closureVelocity(m, slaved, hEdge, sea, 1, new Float64Array(mE)), half = closureVelocity(m, slaved, hEdge, sea, 0.5, new Float64Array(mE));
  let error = 0, size = 0, count = 0;
  for (let e = 0; e < mE; e++) {
    if (!beside[e]) { assert.equal(whole[e], slaved[e]); assert.equal(half[e], slaved[e]); continue; }
    assert.ok(Math.abs(half[e] - 0.5 * (whole[e] + slaved[e])) < 1e-12);
    error += (whole[e] - exact[e]) ** 2; size += exact[e] ** 2; count++;
  }
  const before = pull(slaved), halfPull = pull(half), wholePull = pull(whole);
  console.log(`${count} token edges beside the class: the fit ${(100 * Math.sqrt(error / size)).toFixed(1)}% rms from its flow; the closure on the class edges beside them ${before.toExponential(2)} on the velocity above, ${halfPull.toExponential(2)} with half the fit, ${wholePull.toExponential(2)} with all of it`);
  assert.ok(count > 100 && Math.sqrt(error / size) < 0.1, `the fit is ${Math.sqrt(error / size)} rms from the class's flow`);
  assert.ok(wholePull < halfPull && halfPull < 0.75 * before, `the closure's pull ${before}, ${halfPull}, ${wholePull}`);
});
