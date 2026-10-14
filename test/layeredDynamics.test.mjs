import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { cellVector, laplacianVelocity } from '../js/dynamics/operators.module.js';
import { createOcean, closureVelocity, closureAdjoint, THIN, EPS, PV_FLOOR } from '../js/ocean/layered.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { seawaterDensity, labelTemperature } from '../js/ocean/seawater.module.js';
import { RHO_AIR, DRAG, RHO, DEG, mesh, C, E, zonalWindOnEdges, totalHeatSalt, northwardMixedTransport, slowOcean, deepen } from './helpers/layered.mjs';

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
  const still = build({ shearMixing: false, interiorShearMixing: false, interfacialDrag: 0 }), du0 = still.stages[0][1];
  for (const options of [{ shearMixing: false, interiorShearMixing: false, interfacialDrag: 5e-5 }, { shearMixing: false, interiorShearMixing: true, interfacialDrag: 5e-5 }, { shearMixing: true, interfacialDrag: 5e-5 }]) {
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

test('the wind gives every sea edge\'s column τ/ρ₀, over the top momentumDepth on the mixed layer and the classes at least PV_FLOOR thick, and over the mixed layer alone where it is deeper than that', () => {
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => 300 - 25 * Math.sin(lat) ** 2), ice = new Float64Array(C);
  const column = { ocean: createOcean(mesh, { everySteps: 1, mixedDepth: 10, minimumThickness: 10, shallowestMixedDepth: 10 }), surfaceT, ice };
  const { ocean } = column, L = ocean.layers;
  ocean.initialize(surfaceT, ice);
  deepen(column, [...Array(C).keys()].filter((i) => i % 3 !== 0), 80, -0.02);
  for (let n = 0; n < ocean.u.length; n++) ocean.u[n] = 0.2 * Math.sin(0.7 * Math.floor(n / E) + 0.37 * (n % E));
  const tau = Float64Array.from(zonalWindOnEdges(mesh, (lat) => 8 * Math.cos(3 * lat)), (w) => RHO_AIR * DRAG * 8 * w);
  ocean.setStress(new Float64Array(E), ice);
  ocean.tendency(ocean.state, ocean.stages[0]);
  const still = Float64Array.from(ocean.stages[0][1]);
  ocean.setStress(tau, ice);
  ocean.tendency(ocean.state, ocean.stages[0]);
  const du = ocean.stages[0][1], hEdge = new Float64Array(ocean.shared.hEdge);
  let worst = 0, depthOff = 0, alone = 0, spread = 0, deep = 0;
  for (let e = 0; e < E; e++) {
    if (!ocean.edgeOcean[e] || Math.abs(tau[e]) < 1e-6) continue;
    const push = tau[e] / RHO, h0 = hEdge[e];
    let total = 0, receivers = 0;
    for (let k = 0; k < L; k++) {
      const n = k * E + e, gain = du[n] - still[n];
      total += hEdge[n] * gain;
      if (k > 0 && gain !== 0) { receivers++; assert.ok(hEdge[n] >= PV_FLOOR, `edge ${e}: a class ${hEdge[n]} m thick took the stress`); }
    }
    worst = Math.max(worst, Math.abs(total - push) / Math.abs(push));
    const depth = push / (du[e] - still[e]);
    if (h0 >= 50) {
      deep++;
      assert.equal(receivers, 0, `edge ${e}: a ${h0} m mixed layer shares the stress`);
      depthOff = Math.max(depthOff, Math.abs(depth - h0) / h0);
    } else {
      if (receivers) spread++; else alone++;
      assert.ok(depth > h0 * (1 - 1e-12) && depth < 50 * (1 + 1e-12), `edge ${e}: the mixed layer, ${h0} m at the edge, takes the stress over ${depth} m`);
    }
  }
  console.log(`${deep} sea edges with the mixed layer at least 50 m thick, ${spread} thinner ones sharing the stress with the classes beneath, ${alone} taking it alone: the column's gain within ${worst.toExponential(2)} of τ/ρ₀, the deep mixed layer's divisor within ${depthOff.toExponential(2)} of its thickness`);
  assert.ok(deep > 50 && spread > 50, `${deep} deep and ${spread} shared edges`);
  assert.ok(worst < 1e-12, `a column gains ${worst} off τ/ρ₀`);
  assert.ok(depthOff < 1e-12);
});

test('under ekmanDepth the wind gives every sea edge\'s column τ/ρ₀ over the larger of the mixed layer and the wind\'s Ekman depth min(momentumDepth, max(momentumFloor, c_E u*/|f|)), on the mixed layer and the classes at least stressThickness thick', () => {
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => 300 - 25 * Math.sin(lat) ** 2), ice = new Float64Array(C);
  const tau = Float64Array.from(zonalWindOnEdges(mesh, (lat) => 8 * Math.cos(3 * lat)), (w) => RHO_AIR * DRAG * 8 * w);
  for (const stressThickness of [PV_FLOOR, THIN]) {
    const column = { ocean: createOcean(mesh, { everySteps: 1, mixedDepth: 10, minimumThickness: 10, shallowestMixedDepth: 10, ekmanDepth: 0.4, momentumFloor: 10, stressThickness }), surfaceT, ice };
    const { ocean } = column, L = ocean.layers;
    ocean.initialize(Float64Array.from(surfaceT), ice);
    deepen(column, [...Array(C).keys()].filter((i) => i % 3 !== 0), 80, -0.02);
    for (let n = 0; n < ocean.u.length; n++) ocean.u[n] = 0.2 * Math.sin(0.7 * Math.floor(n / E) + 0.37 * (n % E));
    ocean.setStress(new Float64Array(E), ice);
    ocean.tendency(ocean.state, ocean.stages[0]);
    const still = Float64Array.from(ocean.stages[0][1]);
    ocean.setStress(tau, ice);
    ocean.tendency(ocean.state, ocean.stages[0]);
    const du = ocean.stages[0][1], hEdge = new Float64Array(ocean.shared.hEdge);
    let worst = 0, depthWorst = 0, alone = 0, spread = 0, between = 0, thinReceivers = 0;
    for (let e = 0; e < E; e++) {
      if (!ocean.edgeOcean[e] || Math.abs(tau[e]) < 1e-6) continue;
      let along = 0;
      for (let s = 0; s < mesh.nEdgesOnEdge[e]; s++) { const o = mesh.edgesOnEdge[mesh.maxEdgesOnEdge * e + s]; along += mesh.weightsOnEdge[mesh.maxEdgesOnEdge * e + s] * mesh.dvEdge[o] * ocean.stress[o]; }
      const ustar = Math.sqrt(Math.hypot(ocean.stress[e], along / mesh.dcEdge[e]) / RHO), hm = Math.min(50, Math.max(10, 0.4 * ustar / Math.abs(mesh.fEdge[e])));
      depthWorst = Math.max(depthWorst, Math.abs(ocean.windDepth[e] - hm) / hm);
      if (hm > 10 && hm < 50) between++;
      const push = tau[e] / RHO, h0 = hEdge[e];
      let total = 0, receivers = 0;
      for (let k = 0; k < L; k++) {
        const n = k * E + e, gain = du[n] - still[n];
        total += hEdge[n] * gain;
        if (k > 0 && gain !== 0) { receivers++; assert.ok(hEdge[n] >= stressThickness, `edge ${e}: a class ${hEdge[n]} m thick took the stress`); if (hEdge[n] < PV_FLOOR) thinReceivers++; }
      }
      worst = Math.max(worst, Math.abs(total - push) / Math.abs(push));
      const depth = push / (du[e] - still[e]);
      if (h0 >= hm) {
        alone++;
        assert.equal(receivers, 0, `edge ${e}: a ${h0} m mixed layer over its ${hm} m Ekman depth shares the stress`);
        assert.ok(Math.abs(depth - h0) <= 1e-12 * h0, `edge ${e}: the ${h0} m mixed layer takes the stress over ${depth} m`);
      } else {
        if (receivers) spread++;
        assert.ok(depth > h0 * (1 - 1e-12) && depth < hm * (1 + 1e-12), `edge ${e}: the mixed layer, ${h0} m at the edge, takes the stress over ${depth} m within an Ekman depth of ${hm} m`);
      }
    }
    console.log(`stressThickness ${stressThickness}: ${between} sea edges with the Ekman depth between its bounds, ${alone} with the mixed layer as deep, ${spread} sharing the stress with the classes beneath (${thinReceivers} classes under PV_FLOOR among the receivers); the column's gain within ${worst.toExponential(2)} of τ/ρ₀, the Ekman depth within ${depthWorst.toExponential(2)} of min(50, max(10, 0.4u*/|f|))`);
    assert.ok(between > 100 && alone > 50 && spread > 50, `${between} edges between, ${alone} alone, ${spread} shared`);
    assert.ok(worst < 1e-12, `a column gains ${worst} off τ/ρ₀`);
    assert.ok(depthWorst < 1e-12, `the Ekman depth is ${depthWorst} off`);
    if (stressThickness === THIN) assert.ok(thinReceivers > 0);
  }
});

test('the interfacial drag keeps each column\'s momentum under a 10 m mixed layer, whatever the refill floor', () => {
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => 300 - 25 * Math.sin(lat) ** 2), ice = new Float64Array(C);
  for (const floors of [{}, { minimumThickness: 10, shallowestMixedDepth: 10 }]) {
    const build = (options) => {
      const ocean = createOcean(mesh, { everySteps: 1, mixedDepth: 10, ...floors, ...options });
      ocean.initialize(Float64Array.from(surfaceT), ice);
      for (let n = 0; n < ocean.u.length; n++) ocean.u[n] = 0.2 * Math.sin(0.7 * Math.floor(n / E) + 0.37 * (n % E));
      ocean.tendency(ocean.state, ocean.stages[0]);
      return ocean;
    };
    const du0 = build({ interfacialDrag: 0 }).stages[0][1];
    for (const options of [{ interfacialDrag: 5e-5 }, { shearMixing: true, interfacialDrag: 5e-5 }, { blockDrag: true, interfacialDrag: 5e-5 }]) {
      const dragging = build(options), L = dragging.layers;
      const hEdge = new Float64Array(dragging.shared.hEdge), du = dragging.stages[0][1];
      let worst = 0, largest = 0, thin = 0;
      for (let e = 0; e < E; e++) {
        if (!dragging.edgeOcean[e]) continue;
        if (hEdge[e] < 50) thin++;
        let net = 0;
        for (let k = 0; k < L; k++) {
          const n = k * E + e;
          if (k > 0 && hEdge[n] < THIN) continue;
          const moved = hEdge[n] * (du[n] - du0[n]);
          net += moved; largest = Math.max(largest, Math.abs(moved));
        }
        worst = Math.max(worst, Math.abs(net));
      }
      console.log(`floors ${JSON.stringify(floors)}, ${JSON.stringify(options)}: ${thin} sea edges with the mixed layer under 50 m; the largest momentum moved by the drag ${largest.toExponential(2)} m²/s², the largest column imbalance ${worst.toExponential(2)}`);
      assert.ok(thin > 1000, `only ${thin} thin mixed-layer edges`);
      assert.ok(largest > 1e-7 && worst < 1e-12 * largest, `the drag changes a column's momentum by ${worst} against ${largest} moved`);
    }
  }
});

/*
 * A 10 m mixed layer over 25, 25, 30 and 410 m classes and a deep one, at
 * 45°N under a zonal stress of 0.1 N/m² peaking there, with no stirring
 * and no surface flux: the transport of everything above the deepest
 * class relative to it, and the mixed layer's own velocity relative to
 * it, averaged over two inertial periods after four days.
 */
async function ekmanStack(stacked, latitude = 45) {
  const classes = { densities: [1023, 1024, 1025, 1026, 1026.95], salinities: [35, 35, 35, 34.9, 34.8], bottoms: [35, 60, 90, 500, 12000] };
  const options = { everySteps: 1, thermoclineTilt: 0, salinityProfile: () => 34, mixedDepth: 10, minimumThickness: 10, shallowestMixedDepth: 10, stirring: 0, ...stacked, ...classes };
  const ocean = await slowOcean(mesh, { ocean: options, surfaceT: new Float64Array(C).fill(303) });
  const tau = 0.1, dt = 1350, shape = (lat) => Math.exp(-(((lat / DEG - latitude) / 10) ** 2));
  const stress = zonalWindOnEdges(mesh, (lat) => tau * shape(lat));
  let spacing = 0;
  for (let e = 0; e < E; e++) spacing += mesh.dcEdge[e];
  const half = 0.5 * spacing / E / mesh.radius / DEG, ring = [];
  for (let i = 0; i < C; i++) if (Math.abs(mesh.latCell[i] / DEG - latitude) < half) ring.push(i);
  const area = ring.reduce((sum, i) => sum + mesh.areaCell[i], 0), f = (lat) => 2 * mesh.omega * Math.sin(lat);
  const theory = ring.reduce((sum, i) => sum + mesh.areaCell[i] / area * tau * shape(mesh.latCell[i]) / (RHO * f(mesh.latCell[i])), 0);
  const L = ocean.layers, deepest = L - 1, period = 2 * Math.PI / f(latitude * DEG), steps = Math.round(2 * period / dt);
  await ocean.advance(stress, dt, Math.round(4 * 86400 / dt));
  const vectors = Array.from({ length: L }, () => new Float64Array(3 * C));
  let east = 0, north = 0, mlEast = 0, mlNorth = 0, depth = 0;
  for (let n = 0; n < steps; n++) {
    await ocean.advance(stress, dt);
    const { h, u } = await ocean.download();
    for (let k = 0; k < L; k++) cellVector(mesh, u.subarray(k * E, (k + 1) * E), vectors[k]);
    for (const i of ring) {
      const x = mesh.xCell[3 * i], y = mesh.xCell[3 * i + 1], z = mesh.xCell[3 * i + 2], r = Math.hypot(x, y), w = mesh.areaCell[i] / area / steps;
      const toEast = [-y / r, x / r, 0], toNorth = [-z * x / r, -z * y / r, r];
      const along = (k, axis) => (vectors[k][3 * i] - vectors[deepest][3 * i]) * axis[0] + (vectors[k][3 * i + 1] - vectors[deepest][3 * i + 1]) * axis[1] + (vectors[k][3 * i + 2] - vectors[deepest][3 * i + 2]) * axis[2];
      for (let k = 0; k < deepest; k++) { east += w * h[k * C + i] * along(k, toEast); north += w * h[k * C + i] * along(k, toNorth); }
      mlEast += w * along(0, toEast); mlNorth += w * along(0, toNorth); depth += w * h[i];
    }
  }
  const limited = (await ocean.diagnostics()).oceanLimited;
  await ocean.close();
  return { engine: ocean.engine, ratio: Math.hypot(east, north) / theory, angle: Math.atan2(-north, east) / DEG, speed: Math.hypot(mlEast, mlNorth), theory, depth, limited, f: f(latitude * DEG) };
}

test('under a 10 m mixed layer the wind spread over momentumDepth drives the Ekman transport at 45°N within 7 % of τ/(ρ₀f), the mixed layer moving as a 50 m one would, and at momentumDepth 10 the mixed layer feels its full Coriolis force', async () => {
  const spread = await ekmanStack({ momentumDepth: 50, ekmanDepth: 0 }), alone = await ekmanStack({ momentumDepth: 10, ekmanDepth: 0 });
  const bound = 1.3 * 0.1 / (RHO * spread.f * 50);
  for (const [label, run] of [['momentumDepth 50', spread], ['momentumDepth 10', alone]]) console.log(`${label} on the ${run.engine} ocean, a ${run.depth.toFixed(1)} m mixed layer: the transport above the deepest class ${run.ratio.toFixed(3)} of τ/(ρ₀f) (${run.theory.toFixed(3)} m²/s) at ${run.angle.toFixed(1)}° to the right of the wind; the mixed layer at ${run.speed.toFixed(4)} m/s relative to the deepest class`);
  for (const run of [spread, alone]) {
    assert.ok(Math.abs(run.ratio - 1) < 0.07, `the stack carries ${run.ratio} of τ/(ρ₀f)`);
    assert.ok(Math.abs(run.depth - 10) < 0.5, `the mixed layer is ${run.depth} m`);
    assert.equal(run.limited, 0);
  }
  assert.ok(Math.abs(spread.angle - 90) < 5, `the transport turns ${spread.angle}° from the wind`);
  assert.ok(spread.speed < bound, `the mixed layer runs at ${spread.speed} m/s, more than ${bound}`);
  assert.ok(alone.speed > 3 * spread.speed, `with momentumDepth 10 the mixed layer runs at ${alone.speed} m/s against ${spread.speed}`);
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

test('under interiorShearMixing the drag between two interior classes follows the Pacanowski–Philander viscosity of their Richardson number down to backgroundViscosity over their spacing, while the mixed layer\'s base keeps interfacialDrag', () => {
  const ocean = createOcean(mesh, { everySteps: 1, interiorShearMixing: true }), constant = createOcean(mesh, { everySteps: 1 });
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => 300 - 25 * Math.sin(lat) ** 2);
  for (const o of [ocean, constant]) o.initialize(Float64Array.from(surfaceT), new Float64Array(C));
  const L = ocean.layers, rho = ocean.densities, hEdge = new Float64Array(ocean.shared.hEdge), g = 9.81, rho0 = 1025;
  const errors = [], quiet = [];
  let base = 0, bases = 0, unchanged = 0, pairs = 0;
  for (const speed of [0, 0.05, 0.2, 0.8]) {
    for (let n = 0; n < ocean.u.length; n++) { const k = Math.floor(n / E), e = n % E; ocean.u[n] = (k % 2 ? speed : 0) * (-mesh.xEdge[3 * e + 1] * mesh.nEdge[3 * e] + mesh.xEdge[3 * e] * mesh.nEdge[3 * e + 1]); }
    constant.u.set(ocean.u);
    ocean.tendency(ocean.state, ocean.stages[0]);
    constant.tendency(constant.state, constant.stages[0]);
    for (let e = 0; e < E; e++) {
      if (!ocean.edgeOcean[e] || Math.abs(mesh.latEdge[e]) > 40 * DEG) continue;
      const present = [];
      for (let k = 1; k < L; k++) if (hEdge[k * E + e] >= THIN) present.push(k);
      if (present.length < 2) continue;
      base += Math.abs(ocean.interfaceRate(e, 0, present[0]) - 2e-4); bases++;
      for (let m = 0; m + 1 < present.length; m++) {
        const up = present[m], down = present[m + 1];
        pairs++;
        if (constant.interfaceRate(e, up, down) === 2e-4) unchanged++;
        const dz = Math.max(THIN, 0.5 * (hEdge[up * E + e] + hEdge[down * E + e])), buoyancy = g * (rho[down] - rho[up]) / rho0;
        const shear = (up % 2 === down % 2 ? 0 : speed ** 2 * (mesh.xEdge[3 * e] ** 2 + mesh.xEdge[3 * e + 1] ** 2)) + 1e-12;
        const expected = Math.min((1e-2 / (1 + 5 * buoyancy * dz / shear) ** 2 + 1e-4) / dz, 0.5 * Math.min(hEdge[up * E + e], hEdge[down * E + e]) / 3600), rate = ocean.interfaceRate(e, up, down);
        if (speed === 0) quiet.push(rate * dz / 1e-4);
        else if (Math.abs(mesh.latEdge[e]) < 5 * DEG) errors.push(Math.abs(rate / expected - 1));
      }
    }
  }
  errors.sort((p, q) => p - q); quiet.sort((p, q) => p - q);
  const at = (list, p) => list[Math.floor(p * (list.length - 1))];
  console.log(`${quiet.length} interior interfaces at rest: r Δz / backgroundViscosity ${at(quiet, 0).toFixed(3)}–${at(quiet, 1).toFixed(3)}; ${errors.length} sheared ones within 5° of the equator within ${(100 * at(errors, 0.5)).toFixed(1)}% of the Pacanowski–Philander value at the median; ${bases} mixed-layer bases at interfacialDrag`);
  assert.ok(quiet.length > 1000 && Math.abs(at(quiet, 0) - 1) < 1e-9 && Math.abs(at(quiet, 1) - 1) < 1e-9);
  assert.ok(errors.length > 100 && at(errors, 0.5) < 0.05, `median ${at(errors, 0.5)}`);
  assert.ok(bases > 100 && base === 0 && unchanged === pairs, `mixed-layer bases off interfacialDrag by ${base}; ${unchanged} of ${pairs} interior interfaces at interfacialDrag without the option`);
});

const solidBody = (e) => -mesh.xEdge[3 * e + 1] * mesh.nEdge[3 * e] + mesh.xEdge[3 * e] * mesh.nEdge[3 * e + 1];
const eastward = zonalWindOnEdges(mesh, () => 1);

function dragColumn(ocean, { h0, up, hUp, down, lighter = 0.33 }) {
  const L = ocean.layers, rho = ocean.densities;
  for (let i = 0; i < C; i++) {
    const was = [];
    for (let k = 0; k < L; k++) was.push([ocean.Q[k * C + i] / ocean.h[k * C + i], ocean.W[k * C + i] / ocean.h[k * C + i]]);
    for (let k = 0; k < L; k++) {
      const a = k * C + i, thickness = k === 0 ? h0 : k === up ? hUp : k === down ? ocean.D[i] - h0 - hUp - (L - 3) * EPS : EPS;
      const [T, S] = k === 0 ? [labelTemperature(rho[up] - lighter, was[0][1]), was[0][1]] : was[k];
      ocean.h[a] = thickness; ocean.Q[a] = thickness * T; ocean.W[a] = thickness * S;
    }
  }
}

test('under blockDrag the drag is interfacialDrag where the interface lies a layer spacing or more above the wind depth, the unfloored Pacanowski–Philander viscosity from the wind depth down, the mixed layer\'s base included, and the two blended linearly between', () => {
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => 300 - 25 * Math.sin(lat) ** 2), ice = new Float64Array(C);
  const stress = zonalWindOnEdges(mesh, (lat) => 0.1 * (0.05 + 0.95 * Math.cos(lat) ** 2));
  const build = (options) => {
    const ocean = createOcean(mesh, { everySteps: 1, mixedDepth: 20, ...options });
    ocean.initialize(Float64Array.from(surfaceT), ice);
    ocean.setStress(stress, ice);
    for (let n = 0; n < ocean.u.length; n++) ocean.u[n] = (Math.floor(n / E) % 2 ? 0.2 : 0) * solidBody(n % E);
    ocean.tendency(ocean.state, ocean.stages[0]);
    return ocean;
  };
  const block = build({ blockDrag: true }), twin = build({ shearMixing: true, interfacialDrag: 0 });
  const L = block.layers, hEdge = new Float64Array(block.shared.hEdge);
  const count = { baseHeld: 0, baseBlend: 0, baseFree: 0, interiorHeld: 0, interiorBlend: 0, interiorFree: 0 };
  let differs = 0, wrong = 0, shallow = Infinity, deep = 0;
  for (let e = 0; e < E; e++) {
    if (!block.edgeOcean[e]) continue;
    const present = [0];
    for (let k = 1; k < L; k++) if (hEdge[k * E + e] >= THIN) present.push(k);
    shallow = Math.min(shallow, block.windDepth[e]); deep = Math.max(deep, block.windDepth[e]);
    for (let m = 0; m + 1 < present.length; m++) {
      const up = present[m], down = present[m + 1];
      let depth = 0;
      for (let k = 0; k < down; k++) depth += hEdge[k * E + e];
      const held = Math.min(1, (block.windDepth[e] - depth) / Math.max(THIN, 0.5 * (hEdge[up * E + e] + hEdge[down * E + e])));
      const rate = block.interfaceRate(e, up, down), viscous = twin.interfaceRate(e, up, down);
      count[`${up === 0 ? 'base' : 'interior'}${held >= 1 ? 'Held' : held > 0 ? 'Blend' : 'Free'}`]++;
      if (held >= 1) { if (rate !== 2e-4) wrong++; if (viscous !== 2e-4) differs++; } else if (held <= 0) { if (rate !== viscous || !(rate < 2e-4)) wrong++; } else if (Math.abs(rate - (viscous + (2e-4 - viscous) * held)) > 1e-15) wrong++;
    }
  }
  console.log(`wind depths ${shallow.toFixed(1)}–${deep.toFixed(1)} m over a 20 m mixed layer: ${JSON.stringify(count)} interfaces; ${differs} of those held would differ under the viscosity alone`);
  assert.ok(shallow < 20 && deep > 49, 'the wind depth straddles the mixed layer');
  for (const [name, n] of Object.entries(count)) assert.ok(n > 20, `${name}: ${n} interfaces`);
  assert.ok(differs > 20, 'the held drag is not the viscosity\'s');
  assert.equal(wrong, 0, `${wrong} interfaces off interfacialDrag where held, off the unfloored viscosity from the wind depth down, or off the blend between`);
});

test('under blockDrag the viscosity beneath the Ekman block takes the design\'s reference values over 50 m, from the Richardson number of the full velocity difference, where interfacialDrag over the same distance is the 10⁻² m²/s of a wind-mixed layer', () => {
  const surfaceT = new Float64Array(C).fill(300), ice = new Float64Array(C);
  const lighter = 0.33, g = 9.81, rho0 = 1025, dz = 50, buoyancy = g * lighter / rho0;
  const ocean = createOcean(mesh, { everySteps: 1, blockDrag: true, bathymetry: new Float64Array(C).fill(180) });
  ocean.initialize(Float64Array.from(surfaceT), ice);
  dragColumn(ocean, { h0: 60, up: 8, hUp: 40, down: 12, lighter });
  ocean.setStress(zonalWindOnEdges(mesh, () => 0.05), ice);
  ocean.tendency(ocean.state, ocean.stages[0]);
  const edges = [];
  for (let e = 0; e < E; e++) if (ocean.edgeOcean[e] && Math.abs(mesh.latEdge[e]) < 20 * DEG) edges.push(e);
  const tangent = (field, e) => {
    let sum = 0;
    for (let s = 0; s < mesh.nEdgesOnEdge[e]; s++) { const o = mesh.edgesOnEdge[mesh.maxEdgesOnEdge * e + s]; sum += mesh.weightsOnEdge[mesh.maxEdgesOnEdge * e + s] * mesh.dvEdge[o] * field[o]; }
    return sum / mesh.dcEdge[e];
  };
  const shear = edges.map((e) => eastward[e] ** 2 + tangent(eastward, e) ** 2).sort((p, q) => p - q), typical = Math.sqrt(shear[Math.floor(shear.length / 2)]);
  const reference = [[0.1, 9e-5], [0.25, 4.2e-5], [0.5, 1.8e-5], [1, 7.6e-6], [5, 2.3e-6]], lines = [];
  for (const [target, expected] of reference) {
    const speed = Math.sqrt(buoyancy * dz / target) / typical;
    ocean.u.fill(0);
    for (let e = 0; e < E; e++) ocean.u[8 * E + e] = speed * eastward[e];
    const near = [];
    for (const e of edges) {
      const richardson = buoyancy * dz / (speed * speed * (eastward[e] ** 2 + tangent(eastward, e) ** 2) + 1e-12);
      const rate = ocean.interfaceRate(e, 0, 8), viscosity = (1e-2 / (1 + 5 * richardson) ** 2 + 1e-4) / dz;
      assert.ok(Math.abs(rate / viscosity - 1) < 1e-9, `Ri ${richardson}: r ${rate} against ${viscosity}`);
      if (Math.abs(richardson / target - 1) < 0.02) near.push(rate);
    }
    near.sort((p, q) => p - q);
    const median = near[Math.floor(near.length / 2)];
    lines.push(`Ri ${target}: ${median.toExponential(2)} m/s over ${near.length} edges`);
    assert.ok(near.length >= 5 && Math.abs(median / expected - 1) < 0.05, `Ri ${target}: r ${median} against the design's ${expected} (${near.length} edges)`);
  }
  console.log(`r over 50 m at the design's Richardson numbers: ${lines.join('; ')}`);
  assert.ok(Math.abs(2e-4 * dz / (1e-2 + 1e-4) - 1) < 0.01, 'interfacialDrag over 50 m is the viscosity at Ri 0');
});

test('blockDrag excludes shearMixing and interiorShearMixing, and with an Ekman block deeper than the ocean it leaves the constant drag\'s run unchanged to the bit', () => {
  assert.throws(() => createOcean(mesh, { blockDrag: true, shearMixing: true }), /blockDrag/);
  assert.throws(() => createOcean(mesh, { blockDrag: true, interiorShearMixing: true }), /blockDrag/);
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => 300 - 25 * Math.sin(lat) ** 2), ice = new Float64Array(C), flux = new Float64Array(C);
  const stress = zonalWindOnEdges(mesh, (lat) => 0.1 * Math.cos(3 * lat));
  const run = (options) => {
    const ocean = createOcean(mesh, { everySteps: 1, ekmanDepth: 0, momentumDepth: 6000, ...options });
    ocean.initialize(Float64Array.from(surfaceT), ice);
    for (let n = 0; n < 6; n++) ocean.advance(Float64Array.from(surfaceT), ice, flux, stress, 1350);
    return ocean;
  };
  const constant = run({}), block = run({ blockDrag: true });
  for (const name of ['u', 'h', 'Q', 'W', 'eta']) assert.deepEqual(Array.from(block[name]), Array.from(constant[name]), name);
  let top = 0;
  for (const x of block.u) top = Math.max(top, Math.abs(x));
  assert.ok(top > 0.01, `the run moved (${top} m/s)`);
});

test('under blockDrag the drag across the Ekman block\'s base falls continuously from interfacialDrag to the Richardson-number viscosity as the mixed layer deepens through the wind depth', () => {
  const surfaceT = new Float64Array(C).fill(300), ice = new Float64Array(C);
  const stress = zonalWindOnEdges(mesh, () => 0.08);
  const hUp = 10, depths = [10, 25, 35, 49.9, 50, 50.1, 70, 80];
  const rates = (h0, options) => {
    const ocean = createOcean(mesh, { everySteps: 1, bathymetry: new Float64Array(C).fill(400), ...options });
    ocean.initialize(Float64Array.from(surfaceT), ice);
    dragColumn(ocean, { h0, up: 8, hUp, down: 12 });
    ocean.setStress(stress, ice);
    for (let e = 0; e < E; e++) { ocean.u[e] = -0.3 * eastward[e]; ocean.u[8 * E + e] = 0.4 * eastward[e]; }
    ocean.tendency(ocean.state, ocean.stages[0]);
    const hEdge = new Float64Array(ocean.shared.hEdge), seen = new Map();
    for (let e = 0; e < E; e++) {
      if (!ocean.edgeOcean[e] || Math.abs(mesh.latEdge[e]) > 20 * DEG || ocean.windDepth[e] !== 50) continue;
      assert.ok(Math.abs(hEdge[e] - h0) < 1e-9 && hEdge[8 * E + e] === hUp && hEdge[9 * E + e] === EPS, `edge column ${hEdge[e]}, ${hEdge[8 * E + e]}, ${hEdge[9 * E + e]}`);
      seen.set(e, ocean.interfaceRate(e, 0, 8));
    }
    return seen;
  };
  const rows = depths.map((h0) => ({ h0, block: rates(h0, { blockDrag: true }), viscous: rates(h0, { shearMixing: true, interfacialDrag: 0 }) }));
  for (const { h0, block, viscous } of rows) {
    const held = Math.min(1, (50 - h0 - 7 * EPS) / (0.5 * (h0 + hUp)));
    let worst = 0, low = Infinity, high = 0;
    for (const [e, rate] of block) {
      const expected = held >= 1 ? 2e-4 : viscous.get(e) + (2e-4 - viscous.get(e)) * Math.max(0, held);
      worst = Math.max(worst, Math.abs(rate - expected));
      low = Math.min(low, rate); high = Math.max(high, rate);
    }
    console.log(`mixed layer ${h0} m over the 50 m wind depth: ${block.size} equatorial bases, held ${Math.max(0, held).toFixed(3)}, r from ${low.toExponential(2)} to ${high.toExponential(2)} m/s`);
    assert.ok(block.size > 50 && worst <= 1e-15, `h0 ${h0}: ${block.size} edges, ${worst} off the blend`);
  }
  const first = rows[0], before = rows.find((r) => r.h0 === 49.9), after = rows.find((r) => r.h0 === 50.1), last = rows[rows.length - 1];
  assert.ok([...first.block.values()].every((r) => r === 2e-4), 'a 10 m mixed layer keeps interfacialDrag');
  let jump = 0;
  for (const [e, rate] of before.block) jump = Math.max(jump, Math.abs(after.block.get(e) - rate));
  assert.ok(jump < 0.01 * 2e-4, `across 49.9 to 50.1 m the drag moves by ${jump}`);
  for (const [e, rate] of last.block) assert.equal(rate, last.viscous.get(e), 'a mixed layer deeper than the wind depth and a layer spacing is the viscosity alone');
});

test('under blockDrag the undercurrent slab sets free beneath the Ekman block: a 100 m slab between a 0.3 m/s westward current above and still water below, driven by 0.05 N/m² over 175 m, balances eastward at 0.7 m/s with the mixed layer\'s base below the wind depth, and westward under a 10 m mixed layer as under the constant drag', () => {
  const surfaceT = new Float64Array(C).fill(300), ice = new Float64Array(C);
  const stress = zonalWindOnEdges(mesh, () => 0.05);
  const hSlab = 100, push = 0.05 / (1025 * 175), upper = -0.3;
  const speeds = Array.from({ length: 111 }, (_, n) => -0.5 + 0.02 * n);
  const balance = (options, h0) => {
    const ocean = createOcean(mesh, { everySteps: 1, bathymetry: new Float64Array(C).fill(h0 + hSlab + 80), ...options });
    ocean.initialize(Float64Array.from(surfaceT), ice);
    dragColumn(ocean, { h0, up: 8, hUp: hSlab, down: 12, lighter: 0.33 });
    ocean.setStress(stress, ice);
    const edges = [];
    for (let e = 0; e < E; e++) if (ocean.edgeOcean[e] && Math.abs(mesh.latEdge[e]) < 8 * DEG && Math.abs(eastward[e]) > 0.7 && ocean.windDepth[e] === 50) edges.push(e);
    ocean.tendency(ocean.state, ocean.stages[0]);
    const at = (speed) => {
      for (let e = 0; e < E; e++) { ocean.u[e] = upper * eastward[e]; ocean.u[8 * E + e] = speed * eastward[e]; }
      return edges.map((e) => {
        const r1 = ocean.interfaceRate(e, 0, 8), r2 = ocean.interfaceRate(e, 8, 12);
        return { r1, r2, residual: push + (r1 * (upper - speed) - r2 * speed) / hSlab };
      });
    };
    const roots = [];
    let last = at(speeds[0]);
    for (let n = 1; n < speeds.length; n++) {
      const now = at(speeds[n]);
      for (let m = 0; m < edges.length; m++) if (last[m].residual > 0 && now[m].residual <= 0) roots.push({ speed: speeds[n - 1] + 0.02 * last[m].residual / (last[m].residual - now[m].residual), r1: now[m].r1, r2: now[m].r2 });
      last = now;
    }
    roots.sort((p, q) => p.speed - q.speed);
    return { roots, edges: edges.length, median: roots[Math.floor(roots.length / 2)] };
  };
  const control = balance({}, 60), under = balance({ blockDrag: true }, 60), within = balance({ blockDrag: true }, 10);
  const show = (name, b) => console.log(`${name}: ${b.roots.length} of ${b.edges} equatorial edges balance, median ${b.median.speed.toFixed(3)} m/s at r ${b.median.r1.toExponential(2)} above and ${b.median.r2.toExponential(2)} below`);
  show('constant drag', control); show('blockDrag, mixed layer 60 m', under); show('blockDrag, mixed layer 10 m', within);
  const balanced = (hP, r1, r2) => (hP + r1 * upper) / (r1 + r2);
  assert.ok(control.roots.length > 20 && control.median.speed > -0.095 && control.median.speed < -0.06, `the slab's balance at the constant 2e-4 m/s is ${control.median.speed} m/s`);
  assert.ok(under.roots.length > 20 && under.median.speed > 0.3 && under.median.speed < 1, `the slab's balance beneath the block is ${under.median.speed} m/s`);
  assert.ok(within.roots.length > 20 && within.median.speed < 0 && within.median.r1 > 1.2e-4, `above the wind depth the slab stays coupled to the current above: ${within.median.speed} m/s at r ${within.median.r1}`);
  for (const b of [control, under, within]) assert.ok(Math.abs(b.median.speed / balanced(push * hSlab, b.median.r1, b.median.r2) - 1) < 0.1, 'the balance is the slab formula');
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

test('with rings, closureVelocity fits a class\'s token edges where a thicker class lies beneath both cells, in two rings, and leaves the tokens over the sea floor on the velocity above', () => {
  const m = buildMesh(new Grid(16)), mE = m.nEdges, k = 5;
  const rotation = (e, sign) => sign * 0.3 * (-m.xEdge[3 * e + 1] * m.nEdge[3 * e] + m.xEdge[3 * e] * m.nEdge[3 * e + 1]);
  const inside = (i) => Math.sin(3 * m.lonCell[i]) + 0.5 * Math.cos(5 * m.latCell[i]) > 0.2;
  const east = (e) => Math.cos(m.lonCell[m.cellsOnEdge[2 * e]]) > 0 && Math.cos(m.lonCell[m.cellsOnEdge[2 * e + 1]]) > 0;
  const hEdge = new Float64Array(mE), sea = new Uint8Array(mE).fill(1), exact = new Float64Array(mE), slaved = new Float64Array(mE), deepest = new Float64Array(mE);
  for (let e = 0; e < mE; e++) {
    const present = inside(m.cellsOnEdge[2 * e]) && inside(m.cellsOnEdge[2 * e + 1]);
    hEdge[e] = present ? 12 : EPS;
    exact[e] = rotation(e, 1);
    slaved[e] = present ? exact[e] : rotation(e, -1);
    deepest[e] = east(e) ? k + 3 : k - 1;
  }
  const thick = (e) => hEdge[e] >= THIN, neighbours = (e) => Array.from({ length: m.nEdgesOnEdge[e] }, (_, s) => m.edgesOnEdge[m.maxEdgesOnEdge * e + s]);
  const first = Uint8Array.from({ length: mE }, (_, e) => (!thick(e) && neighbours(e).some(thick) ? 1 : 0));
  const second = Uint8Array.from({ length: mE }, (_, e) => (!thick(e) && !first[e] && neighbours(e).some((o) => first[o] && deepest[o] > k) ? 1 : 0));
  const rings = { deepest, k, valid: new Uint8Array(mE), second: new Float64Array(mE) };
  const out = closureVelocity(m, slaved, hEdge, sea, 1, new Float64Array(mE), rings);
  const error = [0, 0], size = [0, 0], count = [0, 0];
  let floor = 0;
  for (let e = 0; e < mE; e++) {
    if (thick(e) || deepest[e] <= k || !(first[e] || second[e])) { assert.equal(out[e], slaved[e]); if (!thick(e) && (first[e] || second[e])) floor++; continue; }
    const r = first[e] ? 0 : 1;
    error[r] += (out[e] - exact[e]) ** 2; size[r] += exact[e] ** 2; count[r]++;
  }
  const [one, two] = [0, 1].map((r) => Math.sqrt(error[r] / size[r]));
  console.log(`${count[0]} token edges beside the class fitted ${(100 * one).toFixed(1)}% rms from its flow, ${count[1]} next to them ${(100 * two).toFixed(1)}%; ${floor} over the sea floor left on the velocity above`);
  assert.ok(count[0] > 30 && count[1] > 30 && floor > 30, `${count}, ${floor}`);
  assert.ok(one < 0.1 && two < 0.3, `the fits are ${one} and ${two} rms from the class's flow`);
});

test('closureAdjoint returns the closure on the fitted token edges to the edges they were fitted from, so the closure of the filled flow does work on the class\'s own edges alone and never adds energy', () => {
  const m = buildMesh(new Grid(16)), mE = m.nEdges, k = 5;
  const inside = (i) => Math.sin(3 * m.lonCell[i]) + 0.5 * Math.cos(5 * m.latCell[i]) > 0.2;
  const hEdge = new Float64Array(mE), sea = new Uint8Array(mE).fill(1), deepest = new Float64Array(mE);
  for (let e = 0; e < mE; e++) { hEdge[e] = inside(m.cellsOnEdge[2 * e]) && inside(m.cellsOnEdge[2 * e + 1]) ? 12 : EPS; deepest[e] = Math.cos(m.lonCell[m.cellsOnEdge[2 * e]]) > -0.5 ? k + 3 : k - 1; }
  const rings = { deepest, k, valid: new Uint8Array(mE), second: new Float64Array(mE) }, weight = (e) => m.dcEdge[e] * m.dvEdge[e];
  let worst = 0, least = Infinity, unextended = Infinity;
  for (let trial = 0; trial < 5; trial++) {
    const u = Float64Array.from({ length: mE }, (_, e) => (hEdge[e] >= THIN ? Math.sin(12.9898 * e * (trial + 1)) : 0));
    const filled = closureVelocity(m, u, hEdge, sea, 1, new Float64Array(mE), rings);
    const closure = laplacianVelocity(m, laplacianVelocity(m, filled));
    let whole = 0;
    for (let e = 0; e < mE; e++) whole += weight(e) * filled[e] * closure[e];
    let own = 0;
    for (let e = 0; e < mE; e++) if (rings.valid[e] === 1) own += weight(e) * u[e] * closure[e];
    closureAdjoint(m, closure, 1, rings);
    let back = 0;
    for (let e = 0; e < mE; e++) if (rings.valid[e] === 1) back += weight(e) * u[e] * closure[e];
    worst = Math.max(worst, Math.abs(back - whole) / whole); least = Math.min(least, whole); unextended = Math.min(unextended, own / whole);
  }
  console.log(`the work of the closure on the thick edges with the adjoint matches the filled flow's ∇⁴ energy to ${worst.toExponential(2)}; without it as little as ${unextended.toFixed(2)} of it`);
  assert.ok(worst < 1e-9 && least > 0);
});

test('under ekmanDepth 0.4 the wind spread over its own Ekman depth drives the transport of a stack under a 10 m mixed layer within 5 % of τ/(ρ₀f) at about 90° at 45°N, where that depth is 38 m, and within 5 % at 15°N, where it is momentumDepth', async () => {
  const north = await ekmanStack({ momentumDepth: 50, ekmanDepth: 0.4 }), south = await ekmanStack({ momentumDepth: 50, ekmanDepth: 0.4 }, 15);
  const hm = 0.4 * Math.sqrt(0.1 / RHO) / north.f, spread = 0.1 / (RHO * north.f * hm);
  for (const [label, run] of [['45°N', north], ['15°N', south]]) console.log(`ekmanDepth 0.4 at ${label} on the ${run.engine} ocean, a ${run.depth.toFixed(1)} m mixed layer: the transport above the deepest class ${run.ratio.toFixed(3)} of τ/(ρ₀f) (${run.theory.toFixed(3)} m²/s) at ${run.angle.toFixed(1)}° to the right of the wind; the mixed layer at ${run.speed.toFixed(4)} m/s relative to the deepest class`);
  console.log(`at 45°N the peak stress's Ekman depth is ${hm.toFixed(1)} m, over which the stress alone would move the mixed layer at ${spread.toFixed(4)} m/s`);
  assert.ok(Math.abs(hm - 38) < 1);
  for (const run of [north, south]) {
    assert.ok(Math.abs(run.ratio - 1) < 0.05, `the stack carries ${run.ratio} of τ/(ρ₀f)`);
    assert.ok(Math.abs(run.depth - 10) < 0.5, `the mixed layer is ${run.depth} m`);
    assert.equal(run.limited, 0);
  }
  assert.ok(Math.abs(north.angle - 90) < 5, `the transport turns ${north.angle}° from the wind`);
  assert.ok(north.speed > 0.8 * spread && north.speed < 1.3 * spread, `the mixed layer runs at ${north.speed} m/s against ${spread} over the Ekman depth`);
});
