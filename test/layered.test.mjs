import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { cellVector } from '../js/dynamics/operators.module.js';
import { createOcean, LAYER_DENSITIES, runoffOutlets, EPS, EDDY_SLACK } from '../js/ocean/layered.module.js';
import { seawaterDensity, labelTemperature, thermalExpansion } from '../js/ocean/seawater.module.js';
import { FREEZING_POINT } from '../js/physics/ice.module.js';
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
  // a measured drift around 2.5e-9 relative after 100 steps with the
  // 24 classes and the −1 °C polar interior. 5e-9 keeps a margin above
  // that measured drift while still catching a real conservation break
  // (orders of magnitude larger).
  assert.ok(Math.abs(after.heat - before.heat) < 5e-9 * Math.abs(before.heat), `heat ${before.heat} -> ${after.heat}`);
  assert.ok(Math.abs(after.salt - before.salt) < 5e-9 * Math.abs(before.salt), `salt ${before.salt} -> ${after.salt}`);
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

test('load() of a saved ocean without layers starts from the climatology', () => {
  const geography = createGeography(mesh, syntheticTopography(180, 360, (lat, lon) => (Math.abs(lon) < 40 * DEG && lat > 12 * DEG && lat < 48 * DEG ? -4000 : 500)));
  const ocean = createOcean(mesh, { everySteps: 1, geography });
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => 275 + 25 * Math.cos(lat) ** 2), ice = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  const climatology = Float64Array.from(ocean.h);
  ocean.u.fill(0.3);
  ocean.load({ h1: new Float64Array(C).fill(55), u1: new Float64Array(E).fill(0.02) }, surfaceT, ice);
  for (let n = 0; n < climatology.length; n++) assert.equal(ocean.h[n], climatology[n]);
  for (let n = 0; n < ocean.u.length; n++) assert.equal(ocean.u[n], 0);
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

test('the equation of state recovers each class label and expands little near freezing', () => {
  for (const r of LAYER_DENSITIES) assert.ok(Math.abs(seawaterDensity(labelTemperature(r), 35) - r) < 1e-9, `label ${r}`);
  assert.ok(thermalExpansion(273.15, 35) < 7e-5 && thermalExpansion(298.15, 35) > 2.9e-4);
  assert.ok(seawaterDensity(FREEZING_POINT, 34.5) < LAYER_DENSITIES[LAYER_DENSITIES.length - 1], 'polar surface water floats on the deepest class');
});

test('every initial column is statically stable, polar columns under ice included', () => {
  const ocean = createOcean(mesh, { everySteps: 1 });
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => Math.max(FREEZING_POINT, 302 - 35 * Math.sin(lat) ** 2));
  const ice = Float64Array.from(surfaceT, (t) => (t <= FREEZING_POINT ? 1 : 0));
  ocean.initialize(surfaceT, ice);
  const L = ocean.layers;
  for (let i = 0; i < C; i++) {
    let above = seawaterDensity(ocean.Q[i] / ocean.h[i], ocean.W[i] / ocean.h[i]);
    for (let k = 1; k < L; k++) {
      const n = k * C + i;
      if (ocean.h[n] <= 5) continue;
      const r = seawaterDensity(ocean.Q[n] / ocean.h[n], ocean.W[n] / ocean.h[n]);
      assert.ok(Math.abs(r - ocean.densities[k]) < 1e-9, `layer ${k} of cell ${i} starts at ${r} against its label ${ocean.densities[k]}`);
      assert.ok(r >= above - 1e-9, `cell ${i} at ${(mesh.latCell[i] / DEG).toFixed(0)}°: layer ${k} (${r.toFixed(3)}) lies under denser water (${above.toFixed(3)})`);
      above = r;
    }
  }
});

test('interior layers relax back to their label densities without losing heat or salt', () => {
  const ocean = createOcean(mesh, { everySteps: 1 });
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => 300 - 25 * Math.sin(lat) ** 2), ice = new Float64Array(C), flux = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  const L = ocean.layers;
  const error = () => {
    let sum = 0, count = 0;
    for (let i = 0; i < C; i++) for (let k = 2; k < L - 1; k++) {
      const n = k * C + i;
      if (ocean.h[n] <= 20) continue;
      sum += Math.abs(seawaterDensity(ocean.Q[n] / ocean.h[n], ocean.W[n] / ocean.h[n]) - ocean.densities[k]); count++;
    }
    return sum / count;
  };
  for (let i = 0; i < C; i++) for (let k = 2; k < L - 1; k++) {
    const n = k * C + i;
    if (ocean.h[n] > 20) ocean.Q[n] += ocean.h[n] * ((i + k) % 2 ? 0.3 : -0.3);
  }
  const before = totalHeatSalt(ocean, mesh), start = error();
  for (let n = 0; n < 400; n++) ocean.advance(surfaceT, ice, flux, new Float64Array(E), 1350);
  const after = totalHeatSalt(ocean, mesh), end = error();
  assert.ok(Math.abs(after.heat - before.heat) < 1e-8 * before.heat && Math.abs(after.salt - before.salt) < 1e-8 * before.salt, `heat ${before.heat} -> ${after.heat}, salt ${before.salt} -> ${after.salt}`);
  assert.ok(start > 0.04 && end < 0.5 * start, `mean distance from the labels ${start} -> ${end} kg/m³`);
});

test('runoff reaches the coastal sea beside the land it ran off, freshening it by exactly that water', () => {
  const geography = createGeography(mesh, syntheticTopography(180, 360, (lat, lon) => (Math.abs(lon) < 40 * DEG && lat > 12 * DEG && lat < 48 * DEG ? -4000 : 500)));
  const ocean = createOcean(mesh, { everySteps: 1, geography });
  const runoff = Float64Array.from(geography.land, (l) => (l ? 2 : 0));
  ocean.fresh.fill(0);
  ocean.accumulate(null, null, 0, runoff);
  let delivered = 0, ranOff = 0;
  for (let i = 0; i < C; i++) {
    ranOff += mesh.areaCell[i] * runoff[i];
    delivered -= mesh.areaCell[i] * ocean.fresh[i];
    if (ocean.fresh[i] === 0) continue;
    assert.ok(ocean.cellOcean[i], `cell ${i} receives runoff but is land`);
    let coastal = false;
    for (let m = 0; m < mesh.nEdgesOnCell[i]; m++) if (geography.land[mesh.cellsOnCell[mesh.maxEdges * i + m]]) coastal = true;
    assert.ok(coastal, `sea cell ${i} receives runoff but touches no land`);
  }
  assert.ok(ranOff > 0 && Math.abs(delivered - ranOff) < 1e-9 * ranOff, `delivered ${delivered} of ${ranOff}`);
});

test('runoff flows down the terrain to the coast the land slopes toward', () => {
  const west = -40 * DEG, east = 40 * DEG;
  const geography = createGeography(mesh, syntheticTopography(180, 360, (lat, lon) => (lon > west && lon < east && Math.abs(lat) < 50 * DEG ? 100 + 2000 * (lon - west) / (east - west) : -4000)));
  const outlet = runoffOutlets(mesh, geography);
  let land = 0, westward = 0;
  for (let i = 0; i < C; i++) {
    if (!geography.land[i]) { assert.equal(outlet[i], i); continue; }
    land++;
    assert.ok(outlet[i] >= 0 && !geography.land[outlet[i]], `land cell ${i} has no sea outlet`);
    if (mesh.lonCell[i] < 0 || mesh.lonCell[i] > east - 10 * DEG || Math.abs(mesh.latCell[i]) > 30 * DEG) continue;
    const lon = Math.atan2(mesh.xCell[3 * outlet[i] + 1], mesh.xCell[3 * outlet[i]]);
    assert.ok(lon < mesh.lonCell[i], `cell ${i} on the eastern slope drains uphill to lon ${(lon / DEG).toFixed(0)}`);
    if (lon < west + 5 * DEG) westward++;
  }
  assert.ok(land > 50 && westward > 5, `${westward} high cells reach the western coast`);
});

test('wind stress reaches the water under sea ice, scaled by the cover and the transmission factor', () => {
  const ocean = createOcean(mesh, { iceStressTransmission: 0.8, everySteps: 1 });
  const C = mesh.nCells, E = mesh.nEdges;
  const surfaceT = new Float64Array(C).fill(FREEZING_POINT + 2), ice = new Float64Array(C), flux = new Float64Array(C);
  const total = new Float64Array(E).fill(0.1);
  ocean.advance(surfaceT, ice, flux, total, 1350);
  const open = Float64Array.from(ocean.stress);
  ice.fill(1); surfaceT.fill(FREEZING_POINT);
  ocean.advance(surfaceT, ice, flux, total, 1350);
  const full = Float64Array.from(ocean.stress);
  ocean.advance(surfaceT, ice, flux, total, 1350, new Float64Array(C).fill(0.5));
  const half = Float64Array.from(ocean.stress);
  let checked = 0;
  for (let e = 0; e < E; e++) if (open[e] !== 0) { checked++; assert.ok(Math.abs(full[e] - 0.8 * open[e]) < 1e-12 && Math.abs(half[e] - 0.9 * open[e]) < 1e-12, `edge ${e}: open ${open[e]} full ${full[e]} half ${half[e]}`); }
  assert.ok(checked > 0);
});

function buriedBump(options) {
  const ocean = createOcean(mesh, { bathymetry: new Float64Array(C).fill(4000), thermoclineTilt: 0, salinityProfile: () => 35, ...options });
  ocean.initialize(new Float64Array(C).fill(290), new Float64Array(C));
  const centre = [Math.cos(-45 * DEG), 0, Math.sin(-45 * DEG)];
  const reshape = (k, amount) => {
    for (let i = 0; i < C; i++) {
      const x = mesh.xCell[3 * i], y = mesh.xCell[3 * i + 1], z = mesh.xCell[3 * i + 2], r = Math.hypot(x, y, z);
      const distance = Math.acos(Math.min(1, (x * centre[0] + y * centre[1] + z * centre[2]) / r)) * mesh.radius;
      const n = k * C + i, t = ocean.Q[n] / ocean.h[n], s = ocean.W[n] / ocean.h[n];
      ocean.h[n] += amount * Math.exp(-((distance / 2e6) ** 2));
      ocean.Q[n] = ocean.h[n] * t; ocean.W[n] = ocean.h[n] * s;
    }
  };
  return { ocean, reshape };
}

function interfaceRange(ocean, k) {
  let low = Infinity, high = -Infinity;
  for (let i = 0; i < C; i++) {
    let z = 0;
    for (let j = 0; j <= k; j++) z += ocean.h[j * C + i];
    low = Math.min(low, z); high = Math.max(high, z);
  }
  return high - low;
}

function classContents(ocean) {
  return Array.from({ length: ocean.layers }, (_, k) => {
    let volume = 0, heat = 0, salt = 0;
    for (let i = 0; i < C; i++) { const n = k * C + i; volume += mesh.areaCell[i] * ocean.h[n]; heat += mesh.areaCell[i] * ocean.Q[n]; salt += mesh.areaCell[i] * ocean.W[n]; }
    return { volume, heat, salt };
  });
}

test('the eddy transport flattens a buried interface bump and keeps every class\'s volume, heat and salt; with eddyDiffusivity 0 it does nothing', () => {
  const k = 20;
  const { ocean, reshape } = buriedBump({ eddyDiffusivity: 1e6 });
  reshape(k, 200); reshape(k + 1, -200);
  const before = classContents(ocean), h0 = Float64Array.from(ocean.h), start = interfaceRange(ocean, k);
  const ranges = [start];
  for (let n = 0; n < 4; n++) { for (let s = 0; s < 25; s++) ocean.eddyTransport(2700); ranges.push(interfaceRange(ocean, k)); }
  console.log(`the ${ocean.densities[k]}/${ocean.densities[k + 1]} interface's bump over 100 steps: ${ranges.map((r) => r.toFixed(1)).join(' → ')} m`);
  assert.ok(ranges.every((r, n) => n === 0 || r < ranges[n - 1]) && ranges[4] < 0.85 * start, `range ${ranges.join(', ')}`);
  const after = classContents(ocean);
  for (let j = 0; j < ocean.layers; j++) for (const q of ['volume', 'heat', 'salt']) {
    const drift = Math.abs(after[j][q] - before[j][q]) / Math.max(Math.abs(before[j][q]), 1e-300);
    assert.ok(drift < 1e-12, `class ${j} ${q} ${before[j][q]} -> ${after[j][q]}`);
  }
  for (let n = 0; n < ocean.h.length; n++) assert.ok(ocean.h[n] >= Math.min(h0[n], EPS) - 1e-12, `layer ${Math.floor(n / C)} of cell ${n % C} thinned to ${ocean.h[n]}`);
  for (let n = 0; n < 12 * C; n++) assert.equal(ocean.h[n], h0[n], 'the mixed layer and the empty classes under it are untouched');

  const off = buriedBump({ eddyDiffusivity: 0 });
  off.reshape(k, 200); off.reshape(k + 1, -200);
  const untouched = [Float64Array.from(off.ocean.h), Float64Array.from(off.ocean.Q), Float64Array.from(off.ocean.W)];
  for (let s = 0; s < 25; s++) off.ocean.eddyTransport(2700);
  assert.deepEqual([off.ocean.h, off.ocean.Q, off.ocean.W].map((a) => Array.from(a)), untouched.map((a) => Array.from(a)));
});

test('a class with only its token thickness carries no eddy flux, and the interfaces around it flatten together', () => {
  const k = 19;
  const { ocean, reshape } = buriedBump({ eddyDiffusivity: 1e6 });
  for (let i = 0; i < C; i++) {
    const n = k * C + i, below = (k + 1) * C + i, moved = ocean.h[n] - EPS, t = ocean.Q[below] / ocean.h[below], s = ocean.W[below] / ocean.h[below];
    ocean.h[below] += moved; ocean.Q[below] = ocean.h[below] * t; ocean.W[below] = ocean.h[below] * s;
    ocean.Q[n] *= EPS / ocean.h[n]; ocean.W[n] *= EPS / ocean.h[n]; ocean.h[n] = EPS;
  }
  reshape(k - 1, 200); reshape(k + 1, -200);
  const token = Float64Array.from(ocean.h.subarray(k * C, (k + 1) * C)), above = interfaceRange(ocean, k - 1), below = interfaceRange(ocean, k);
  for (let s = 0; s < 100; s++) ocean.eddyTransport(2700);
  assert.deepEqual(Array.from(ocean.h.subarray(k * C, (k + 1) * C)), Array.from(token));
  const aboveAfter = interfaceRange(ocean, k - 1), belowAfter = interfaceRange(ocean, k);
  console.log(`around the empty ${ocean.densities[k]} class the interfaces' bumps went ${above.toFixed(1)} → ${aboveAfter.toFixed(1)} and ${below.toFixed(1)} → ${belowAfter.toFixed(1)} m`);
  assert.ok(aboveAfter < 0.85 * above && Math.abs(aboveAfter - belowAfter) < 0.02, `${above} -> ${aboveAfter}, ${below} -> ${belowAfter}`);
});

test('the eddy transport leaves level interfaces over a bumpy bottom level, and at the stability limit stays within each class\'s water', () => {
  const bathymetry = Float64Array.from({ length: C }, (_, i) => 2500 + 2000 * Math.sin(3 * mesh.lonCell[i]) * Math.cos(2 * mesh.latCell[i]));
  const flat = createOcean(mesh, { bathymetry, thermoclineTilt: 0, salinityProfile: () => 35, eddyDiffusivity: 1e12 });
  flat.initialize(new Float64Array(C).fill(290), new Float64Array(C));
  const resting = Float64Array.from(flat.h);
  for (let s = 0; s < 50; s++) flat.eddyTransport(2700);
  let moved = 0, offset = 0;
  for (let n = 0; n < resting.length; n++) moved = Math.max(moved, Math.abs(flat.h[n] - resting[n]));
  for (let k = 1; k < flat.layers - 1; k++) {
    let low = Infinity, high = -Infinity;
    for (let i = 0; i < C; i++) {
      let z = 0, below = 0;
      for (let j = 0; j < flat.layers; j++) if (j <= k) z += resting[j * C + i]; else below += resting[j * C + i];
      if (below > 100) { low = Math.min(low, z); high = Math.max(high, z); }
    }
    if (high > low) offset = Math.max(offset, high - low);
  }
  console.log(`over the bumpy bottom the interfaces start within ${offset.toFixed(3)} m of level where they lie 100 m above the bottom, and the eddy transport at the stability limit moves them by ${moved.toFixed(3)} m in 50 steps`);
  assert.ok(moved < 0.5, `level interfaces over a bumpy bottom moved by ${moved} m`);

  const ocean = createOcean(mesh, { bathymetry, eddyDiffusivity: 1e12 });
  ocean.initialize(Float64Array.from(mesh.latCell, (lat) => 272 + 30 * Math.cos(lat) ** 2), new Float64Array(C));
  const before = classContents(ocean), h0 = Float64Array.from(ocean.h), previous = Float64Array.from(ocean.h);
  let thinnest = Infinity;
  for (let s = 0; s < 200; s++) {
    ocean.eddyTransport(2700);
    for (let n = 0; n < previous.length; n++) {
      thinnest = Math.min(thinnest, ocean.h[n] - Math.min(previous[n], EPS));
      previous[n] = ocean.h[n];
    }
  }
  assert.ok(thinnest >= -6 * EDDY_SLACK, `a class fell ${-thinnest} m below its token or its thickness before the step`);
  const after = classContents(ocean);
  let changed = 0;
  for (let n = 0; n < h0.length; n++) changed = Math.max(changed, Math.abs(ocean.h[n] - h0[n]));
  for (let j = 0; j < ocean.layers; j++) assert.ok(Math.abs(after[j].volume - before[j].volume) < 1e-12 * before[j].volume && Math.abs(after[j].heat - before[j].heat) < 1e-12 * before[j].heat, `class ${j}`);
  for (let i = 0; i < C; i++) {
    let sum = 0;
    for (let j = 0; j < ocean.layers; j++) sum += ocean.h[j * C + i];
    assert.ok(Math.abs(sum - ocean.D[i] - ocean.eta[i]) < 1e-9, `column ${i}`);
  }
  assert.ok(changed > 100, `the tilted thermocline moved by at most ${changed} m`);
});
