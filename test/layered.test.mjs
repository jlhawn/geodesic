import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { cellVector } from '../js/dynamics/operators.module.js';
import { createOcean } from '../js/ocean/layered.module.js';
import { createGeography, syntheticTopography } from '../js/geography.module.js';

const RHO_AIR = 1.2, DRAG = 1.5e-3, RHO = 1025, DEG = Math.PI / 180;
const mesh = buildMesh(new Grid(8));
const { nCells: C, nEdges: E } = mesh;

function zonalWindOnEdges(m, speedAt) {
  const u = new Float64Array(m.nEdges);
  for (let e = 0; e < m.nEdges; e++) {
    const x = m.xEdge[3 * e], y = m.xEdge[3 * e + 1], r = Math.hypot(x, y);
    const east = r > 0 ? [-y / r, x / r, 0] : [0, 0, 0];
    u[e] = speedAt(m.latEdge[e]) * (east[0] * m.nEdge[3 * e] + east[1] * m.nEdge[3 * e + 1] + east[2] * m.nEdge[3 * e + 2]);
  }
  return u;
}

// Area-weighted sums of h·T and h·S over every layer: what the module's own
// flux-divergence advection and mixing conserve exactly, independent of the
// rhoCp-scaled diagnostics().oceanHeat proxy.
function totalHeatSalt(ocean, m) {
  let heat = 0, salt = 0;
  for (let i = 0; i < m.nCells; i++) {
    if (!ocean.cellOcean[i]) continue;
    for (let k = 0; k < ocean.layers; k++) { heat += m.areaCell[i] * ocean.Q[k * m.nCells + i]; salt += m.areaCell[i] * ocean.W[k * m.nCells + i]; }
  }
  return { heat, salt };
}

function northwardMixedTransport(ocean, m) {
  const vector = cellVector(m, ocean.u.subarray(0, m.nEdges));
  const out = new Float64Array(m.nCells);
  for (let i = 0; i < m.nCells; i++) {
    const x = m.xCell[3 * i], y = m.xCell[3 * i + 1], z = m.xCell[3 * i + 2], r = Math.hypot(x, y);
    const north = [-z * x / r, -z * y / r, r];
    out[i] = ocean.h[i] * (vector[3 * i] * north[0] + vector[3 * i + 1] * north[1] + vector[3 * i + 2] * north[2]);
  }
  return out;
}

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

test('westerlies at 45N and 45S drive an equatorward mixed-layer Ekman transport near tau/(rho f)', () => {
  const ocean = createOcean(mesh, { everySteps: 1 });
  const surfaceT = new Float64Array(C).fill(290), ice = new Float64Array(C), flux = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  const speed = 10;
  const wind = zonalWindOnEdges(mesh, (lat) => speed * Math.exp(-(((Math.abs(lat) * 180 / Math.PI - 45) / 10) ** 2)));
  const stressField = Float64Array.from(wind, (w) => RHO_AIR * DRAG * speed * w);
  const dt = 1350;
  const f = 2 * mesh.omega * Math.sin(45 * DEG);
  const inertialPeriod = 2 * Math.PI / f;
  for (let n = 0; n < 12 * 64; n++) ocean.advance(surfaceT, ice, flux, stressField, dt);
  // A suddenly-imposed wind excites a ~17h inertial oscillation in the mixed
  // layer that an instantaneous snapshot aliases against daily sampling;
  // average over two inertial periods to recover the mean Ekman balance.
  const windowSteps = Math.round(2 * inertialPeriod / dt);
  let north = 0, south = 0, samples = 0;
  for (let n = 0; n < windowSteps; n++) {
    ocean.advance(surfaceT, ice, flux, stressField, dt);
    const transport = northwardMixedTransport(ocean, mesh);
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
  assert.ok(north < 0 && Math.abs(-north - ekman) < 0.25 * ekman, `NH transport ${north} vs Ekman ${-ekman} m^2/s`);
  assert.ok(south > 0 && Math.abs(south - ekman) < 0.25 * ekman, `SH transport ${south} vs Ekman ${ekman} m^2/s`);
  assert.equal(ocean.diagnostics().oceanLimited, 0);
});

test('a subtropical wind over a closed 80-degree basin drives a Stommel western boundary current', () => {
  const basin = buildMesh(new Grid(16));
  const { nCells: bC, nEdges: bE } = basin;
  const geography = createGeography(basin, syntheticTopography(180, 360, (lat, lon) => (Math.abs(lon) < 40 * DEG && lat > 12 * DEG && lat < 48 * DEG ? -4000 : 500)));
  const ocean = createOcean(basin, { everySteps: 1, geography });
  const surfaceT = Float64Array.from(basin.latCell, (lat) => 275 + 25 * Math.cos(lat) ** 2), ice = new Float64Array(bC), flux = new Float64Array(bC);
  ocean.initialize(surfaceT, ice);
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
  for (let n = 0; n < days * 64; n++) ocean.advance(surfaceT, ice, flux, stressField, 1350);

  const Ubt = new Float64Array(bE);
  for (let e = 0; e < bE; e++) {
    let t = 0;
    const a = basin.cellsOnEdge[2 * e], b = basin.cellsOnEdge[2 * e + 1];
    for (let k = 0; k < ocean.layers; k++) t += 0.5 * (ocean.h[k * bC + a] + ocean.h[k * bC + b]) * ocean.u[k * bE + e];
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
  assert.ok(west > 0, `western boundary transport ${west} m^2/s should be northward`);
  assert.ok(interior < 0, `interior transport ${interior} m^2/s should be southward`);
  assert.ok(west >= 5 * Math.abs(interior), `western transport ${west} should be at least 5x the interior's ${interior}`);

  let maxBin = -Infinity, maxCenter = null;
  for (let center = -38; center <= 38; center += 4) {
    const v = band(center - 2, center + 2);
    if (Number.isFinite(v) && v > maxBin) { maxBin = v; maxCenter = center; }
  }
  assert.ok(maxCenter === -38 || maxCenter === -34, `the largest 4-degree band is centred at ${maxCenter}, not one of the westernmost two`);
  assert.equal(ocean.diagnostics().oceanLimited, 0);
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
  // At rest (eta ~ microns) the module conserves heat to ~1e-14 relative (see
  // the first test above); under this wind the free surface swings by tens
  // of centimetres, and the per-cell rescale to the sub-stepped barotropic
  // eta at the end of step() is not perfectly heat/salt-conservative, giving
  // a measured drift around 2-9e-10 relative after 100 steps. 1e-9 keeps a
  // comfortable margin above that measured drift while still catching a
  // real conservation break (which is orders of magnitude larger).
  assert.ok(Math.abs(after.heat - before.heat) < 1e-9 * Math.abs(before.heat), `heat ${before.heat} -> ${after.heat}`);
  assert.ok(Math.abs(after.salt - before.salt) < 1e-9 * Math.abs(before.salt), `salt ${before.salt} -> ${after.salt}`);
  assert.equal(ocean.diagnostics().oceanLimited, 0);
});

test('surface warming shallows the mixed layer by detrainment; surface cooling deepens it by entrainment', () => {
  const ocean = createOcean(mesh, { everySteps: 1 });
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => 275 + 25 * Math.cos(lat) ** 2), ice = new Float64Array(C), flux = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  const wind = zonalWindOnEdges(mesh, () => 6);
  const stressField = Float64Array.from(wind, (w) => RHO_AIR * DRAG * 6 * w);
  const dt = 1350;
  const d0 = ocean.diagnostics().oceanUpperDepth;
  for (let n = 0; n < 10; n++) {
    for (let i = 0; i < C; i++) surfaceT[i] += 1;
    ocean.advance(surfaceT, ice, flux, stressField, dt);
  }
  const d1 = ocean.diagnostics().oceanUpperDepth;
  assert.ok(d1 < d0, `mixed layer depth ${d0} -> ${d1} m under warming should decrease`);
  for (let n = 0; n < 10; n++) {
    for (let i = 0; i < C; i++) surfaceT[i] -= 1;
    ocean.advance(surfaceT, ice, flux, stressField, dt);
  }
  const d2 = ocean.diagnostics().oceanUpperDepth;
  assert.ok(d2 > d1, `mixed layer depth ${d1} -> ${d2} m under cooling should increase`);
  assert.equal(ocean.diagnostics().oceanLimited, 0);
});

test('load() from the old two-layer format seeds the climatology and honours h1 within [minimumThickness, D] and u1 on ocean edges', () => {
  const geography = createGeography(mesh, syntheticTopography(180, 360, (lat, lon) => (Math.abs(lon) < 40 * DEG && lat > 12 * DEG && lat < 48 * DEG ? -4000 : 500)));
  const ocean = createOcean(mesh, { everySteps: 1, geography });
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => 275 + 25 * Math.cos(lat) ** 2), ice = new Float64Array(C);

  // The climatology's own mixed layer is ~60 m deep everywhere (mixedDepth),
  // so a shallowing request to 55 m only ever gives mass to the interior
  // (unbounded above) and is honoured exactly, unlike a deepening request
  // (see below).
  const u1 = Float64Array.from({ length: E }, (_, e) => 0.02 * Math.sin(e));
  const saved = { h1: new Float64Array(C).fill(55), h2: new Float64Array(C).fill(900), u1, u2: new Float64Array(E), T2: new Float64Array(C).fill(280) };
  ocean.load(saved, surfaceT, ice);
  for (let i = 0; i < C; i++) {
    if (!ocean.cellOcean[i]) continue;
    assert.ok(Math.abs(ocean.h[i] - 55) < 1e-6, `mixed layer depth ${ocean.h[i]} at cell ${i} should honour the requested 55 m`);
  }
  for (let e = 0; e < E; e++) {
    if (ocean.edgeOcean[e]) assert.equal(ocean.u[e], u1[e]);
    else assert.equal(ocean.u[e], 0);
  }

  // Extreme requests still land within [minimumThickness, D]: the loader
  // only exchanges mass with the single interior layer directly below the
  // mixed layer, so a huge h1 request is capped by that layer's own
  // thickness long before it reaches the seafloor.
  const extreme = { h1: new Float64Array(C).fill(1e6), h2: new Float64Array(C).fill(1), u1: new Float64Array(E).fill(0.3), u2: new Float64Array(E), T2: new Float64Array(C) };
  ocean.load(extreme, surfaceT, ice);
  for (let i = 0; i < C; i++) {
    if (!ocean.cellOcean[i]) continue;
    assert.ok(ocean.h[i] >= 10 - 1e-9 && ocean.h[i] <= ocean.D[i] + 1e-9, `mixed layer depth ${ocean.h[i]} at cell ${i} should stay within [minimumThickness, D=${ocean.D[i]}]`);
  }

  const shallow = { h1: new Float64Array(C).fill(0.001), h2: new Float64Array(C).fill(900), u1: new Float64Array(E).fill(-0.1), u2: new Float64Array(E), T2: new Float64Array(C) };
  ocean.load(shallow, surfaceT, ice);
  for (let i = 0; i < C; i++) {
    if (!ocean.cellOcean[i]) continue;
    assert.ok(ocean.h[i] >= 10 - 1e-9, `mixed layer depth ${ocean.h[i]} at cell ${i} should not go below minimumThickness`);
  }
});

test('load(serialize()) reproduces h, u and eta exactly', () => {
  const ocean = createOcean(mesh, { everySteps: 1 });
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => 275 + 25 * Math.cos(lat) ** 2), ice = new Float64Array(C), flux = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  const wind = zonalWindOnEdges(mesh, (lat) => 8 * Math.cos(3 * lat));
  const stressField = Float64Array.from(wind, (w) => RHO_AIR * DRAG * 8 * w);
  for (let n = 0; n < 20; n++) ocean.advance(surfaceT, ice, flux, stressField, 1350);

  const saved = ocean.serialize();
  const reloaded = createOcean(mesh, { everySteps: 1 });
  reloaded.load(saved, Float64Array.from(surfaceT), Float64Array.from(ice));

  assert.deepEqual(Array.from(reloaded.h), Array.from(ocean.h));
  assert.deepEqual(Array.from(reloaded.u), Array.from(ocean.u));
  assert.deepEqual(Array.from(reloaded.eta), Array.from(ocean.eta));
  assert.equal(ocean.diagnostics().oceanLimited, 0);
});

test('load() fits carried-over columns to the bathymetry and fills sea cells that arrive empty', () => {
  const geography = createGeography(mesh, syntheticTopography(180, 360, (lat, lon) => (Math.abs(lon) < 40 * DEG && lat > 12 * DEG && lat < 48 * DEG ? -4000 : 500)));
  const ocean = createOcean(mesh, { geography });
  const C = mesh.nCells, surfaceT = new Float64Array(C).fill(290), ice = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  const saved = ocean.serialize();
  const L = saved.h.length / C, seaCells = [];
  for (let i = 0; i < C; i++) if (ocean.cellOcean[i]) seaCells.push(i);
  const tooDeep = seaCells[0], tooShallow = seaCells[1], empty = seaCells[2];
  for (let k = 0; k < L; k++) saved.h[k * C + tooDeep] *= 3;
  for (let k = 0; k < L; k++) saved.h[k * C + tooShallow] *= 0.2;
  for (let k = 0; k < L; k++) { saved.h[k * C + empty] = 0; saved.T[k * C + empty] = 0; saved.S[k * C + empty] = 0; }
  saved.eta[tooDeep] = 40;
  ocean.load(saved, surfaceT, ice);
  for (const i of seaCells) {
    let sum = 0;
    for (let k = 0; k < L; k++) { assert.ok(ocean.h[k * C + i] > 0, `layer ${k} at cell ${i} has water or a token`); sum += ocean.h[k * C + i]; }
    assert.ok(Math.abs(sum - ocean.D[i] - ocean.eta[i]) < 1e-6, `column ${i} sums to its depth plus sea level (${sum} vs ${ocean.D[i]} + ${ocean.eta[i]})`);
    assert.ok(Math.abs(ocean.eta[i]) <= 5, `sea level at cell ${i} is within 5 m (${ocean.eta[i]})`);
    assert.ok(ocean.h[i] >= 50 - 1e-9 || ocean.h[i] >= ocean.D[i] - 1, `mixed layer at cell ${i} keeps its floor (${ocean.h[i]})`);
    assert.ok(ocean.T0[i] > 250 && ocean.T0[i] < 320, `SST at cell ${i} is physical (${ocean.T0[i]})`);
  }
  assert.ok(Math.abs(ocean.T0[empty] - 290) < 1e-6, 'an empty sea cell takes the climatology surface temperature');
});
