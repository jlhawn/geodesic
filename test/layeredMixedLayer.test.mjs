import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createOcean, LAYER_BOTTOMS, transmitted, THIN, langmuirFactor, stokesSurfaceLayer, FRICTION_PER_WIND } from '../js/ocean/layered.module.js';
import { SEAWATER, seawaterDensity, thermalExpansion, labelTemperature } from '../js/ocean/seawater.module.js';
import { cellVector } from '../js/dynamics/operators.module.js';
import { RHO_AIR, DRAG, RHO, DEG, mesh, C, E, zonalWindOnEdges, UNIFORM, uniformOcean, deepen, slowOcean, columnMomentum } from './helpers/layered.mjs';

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

test('a convectively neutral mixed layer holds its depth; under neutralSnap it returns to shallowestMixedDepth', () => {
  for (const neutralSnap of [false, true]) {
    const column = uniformOcean({ neutralSnap });
    const { ocean, surfaceT, ice, flux, calm } = column;
    deepen(column, [...Array(C).keys()], 300, -0.002);
    for (let n = 0; n < 8; n++) ocean.advance(surfaceT, ice, flux, calm, 1350);
    const expected = neutralSnap ? 50 : 300;
    for (let i = 0; i < C; i++) assert.ok(Math.abs(ocean.h[i] - expected) < 1e-6, `neutralSnap ${neutralSnap}: cell ${i} is ${ocean.h[i]} m deep, not ${expected}`);
  }
});

test('a deep mixed layer retreats toward shallowestMixedDepth over detrainmentTime while the surface gains buoyancy, and holds without a flux', () => {
  const dt = 1350, steps = 4;
  for (const warming of [0.05, 0]) {
    const column = uniformOcean();
    const { ocean, surfaceT, ice, flux, calm } = column;
    deepen(column, [...Array(C).keys()], 300, -0.002);
    for (let n = 0; n < steps; n++) {
      for (let i = 0; i < C; i++) surfaceT[i] += warming;
      ocean.advance(surfaceT, ice, flux, calm, dt);
    }
    const expected = warming > 0 ? 50 + 250 * (1 - dt / 86400) ** steps : 300;
    for (let i = 0; i < C; i++) assert.ok(Math.abs(ocean.h[i] - expected) < 1e-4, `warming ${warming} K a step: cell ${i} is ${ocean.h[i]} m deep, not ${expected}`);
  }
});

/*
 * Open water under a steady 0.05 N/m² wind gaining 50 W/m² net, 200 W/m²
 * of it sunlight, its surface updated as the atmosphere's step does it:
 * with penetration the share the ocean publishes goes to the waiting light
 * rather than the SST. The mixed layer retreats toward the Monin–Obukhov
 * depth L = 2 m u*³/B of the heating it keeps, so with penetration toward
 * the depth h that solves h = L(h), S·(transmitted(h) − transmitted(D))
 * of the sunlight passing below it.
 */
test('with the same net flux, part of it sunlight, the mixed layer retreats to a deeper Monin–Obukhov depth with penetration than without, at the depth where L(h) = h', () => {
  const dt = 2700, g = 9.81, tau = 0.05, sunlight = 200, net = 50, stirring = 0.8, stirringDepth = 100, days = 4;
  const band = [...Array(C).keys()].filter((i) => Math.abs(mesh.latCell[i]) < 50 * DEG), runs = {};
  for (const penetration of [null, 'IB']) {
    const { ocean, surfaceT, ice, flux } = uniformOcean({ penetration, minimumThickness: 5, shallowestMixedDepth: 5, buoyancyMemory: 0, stirring, stirringDepth });
    const stress = zonalWindOnEdges(mesh, () => tau);
    for (let n = 0; n < days * 86400 / dt; n++) {
      for (let i = 0; i < C; i++) {
        const passed = penetration ? sunlight * ocean.transmission[i] : 0;
        surfaceT[i] += dt * (net - passed) / ocean.capacity[i];
        if (penetration) ocean.penetrated[i] += dt * passed;
      }
      ocean.advance(surfaceT, ice, flux, stress, dt);
    }
    let depth = 0, target = 0, worst = 0;
    for (const i of band) {
      const alpha = thermalExpansion(surfaceT[i], ocean.W[i] / ocean.h[i]), bottom = ocean.D[i] + ocean.eta[i];
      const L = (h) => 2 * stirring * Math.exp(-h / stirringDepth) * (tau / RHO) ** 1.5 * ocean.rhoCp / (g * alpha * (net - (penetration ? sunlight * (transmitted(h) - transmitted(bottom)) : 0)));
      let h = 10;
      for (let it = 0; it < 100; it++) h = L(h);
      depth += ocean.h[i] / band.length; target += h / band.length; worst = Math.max(worst, Math.abs(ocean.h[i] / h - 1));
    }
    runs[penetration ?? 'none'] = { depth, target, worst };
    console.log(`penetration ${penetration}: after ${days} days the mixed layer within 50° of the equator is ${depth.toFixed(2)} m deep against a mean L(h) = h of ${target.toFixed(2)} m, each column within ${(100 * worst).toFixed(1)} % of its own`);
  }
  const { none, IB } = runs;
  assert.ok(Math.abs(IB.depth / IB.target - 1) < 0.03 && IB.worst < 0.05, `with penetration the layer is ${IB.depth} m against ${IB.target}, a column ${IB.worst} off`);
  assert.ok(IB.target > none.target + 3 && IB.depth - none.depth > 0.75 * (IB.target - none.target), `penetration deepens the layer by ${IB.depth - none.depth} m against ${IB.target - none.target}`);
});

test('a mixed layer within the density range of the class beneath erodes it at B/(h N²) while it loses buoyancy, N² the class\'s remaining density span over its thickness', () => {
  const dt = 1350, g = 9.81, rho0 = 1025;
  const column = uniformOcean({ buoyancyMemory: 0, convectiveRate: 1e5 / 86400 });
  const { ocean, surfaceT, ice, flux, calm } = column;
  deepen(column, [...Array(C).keys()], 300, -0.002);
  const i = 0, L = ocean.layers;
  let k = 1;
  while (ocean.h[k * C + i] <= 5) k++;
  const H = ocean.h[k * C + i], s = ocean.W[i] / ocean.h[i], rho = ocean.densities;
  for (let j = 0; j < C; j++) surfaceT[j] -= 0.05;
  const t = surfaceT[i], rm = seawaterDensity(t, s), densest = rho[k] + 0.5 * (k < L - 1 ? rho[k + 1] - rho[k] : rho[k] - rho[k - 1]);
  const loss = g * thermalExpansion(t, s) * 0.05 * 300 / dt;
  const expected = 300 + loss * dt * rho0 * H / (g * 300 * (densest - rm));
  ocean.advance(surfaceT, ice, flux, calm, dt);
  console.log(`a 300 m layer ${(rm - rho[k]).toFixed(4)} kg/m³ from the label of the ${H.toFixed(0)} m class beneath, cooled 0.05 K in ${dt} s (B = ${loss.toExponential(2)} m²/s³): ${ocean.h[i].toFixed(4)} m against ${expected.toFixed(4)}`);
  assert.ok(expected > 300.05 && Math.abs(ocean.h[i] - expected) < 1e-6 * expected, `${ocean.h[i]} m against ${expected}`);
});

test('a mixed layer deeper than mixedNeighbourRatio times its neighbours\' mean detrains the excess over detrainmentTime, and convection stops at that depth', () => {
  const dt = 1350, i = 0, around = [];
  for (let m = 0; m < mesh.nEdgesOnCell[i]; m++) around.push(mesh.cellsOnCell[mesh.maxEdges * i + m]);
  const meanAround = (ocean) => around.reduce((sum, j) => sum + ocean.h[j], 0) / around.length;
  for (const ratio of [3, 0]) {
    const column = uniformOcean({ mixedNeighbourRatio: ratio });
    deepen(column, [i], 600, -0.05);
    column.ocean.advance(column.surfaceT, column.ice, column.flux, column.calm, dt);
    const reach = 3 * meanAround(column.ocean), expected = ratio ? 600 - (600 - reach) * dt / 86400 : 600;
    assert.ok(Math.abs(column.ocean.h[i] - expected) < 0.5, `ratio ${ratio}: a stable 600 m layer among ${reach / 3} m ones is ${column.ocean.h[i]} m after a step, not ${expected}`);

    const dense = uniformOcean({ mixedNeighbourRatio: ratio, convectiveRate: 1e5 / 86400 });
    deepen(dense, [i], 50, 0.5);
    dense.ocean.advance(dense.surfaceT, dense.ice, dense.flux, dense.calm, dt);
    const cap = 3 * meanAround(dense.ocean);
    if (ratio) assert.ok(Math.abs(dense.ocean.h[i] - cap) < 0.1, `convection took the dense layer to ${dense.ocean.h[i]} m, not ${cap}`);
    else assert.ok(dense.ocean.h[i] > cap + 50, `without the guard convection takes the dense layer past ${cap} m (${dense.ocean.h[i]})`);
  }
});

test('a neutral 600 m mixed layer beside 50 m ones at 60°S moves the free surface by about its steric deficit and no water faster than 1 m/s over five days at N=16', async () => {
  const coarse = buildMesh(new Grid(16));
  const offset = (i) => Math.abs(coarse.latCell[i] + 60 * DEG) + Math.abs(coarse.lonCell[i]);
  let centre = 0;
  for (let i = 1; i < coarse.nCells; i++) if (offset(i) < offset(centre)) centre = i;
  const beside = coarse.cellsOnCell[coarse.maxEdges * centre];
  let eta0, contrast, steric;
  const ocean = await slowOcean(coarse, {
    ocean: { ...UNIFORM, mixedNeighbourRatio: 0, bottoms: LAYER_BOTTOMS.map((z) => 0.6 * z) },
    surfaceT: new Float64Array(coarse.nCells).fill(278),
    start(column) {
      deepen(column, [centre], 600, -0.002);
      const { ocean } = column, rhoMl = (i) => seawaterDensity(ocean.Q[i] / ocean.h[i], ocean.W[i] / ocean.h[i]);
      const upperMass = (i) => {
        let z = 0, mass = 0;
        for (let k = 0; k < ocean.layers && z < 1000; k++) { const part = Math.min(ocean.h[k * coarse.nCells + i], 1000 - z); mass += part * (k ? ocean.densities[k] : rhoMl(i)); z += part; }
        return mass;
      };
      eta0 = Float64Array.from(ocean.eta);
      contrast = rhoMl(centre) - rhoMl(beside);
      steric = (upperMass(centre) - upperMass(beside)) / RHO;
    },
  });
  console.log(`the neutral 600 m mixed layer at N=16 on the ${ocean.engine} ocean`);
  const calm = new Float64Array(coarse.nEdges);
  let surface = 0, fastest = 0, state;
  for (let n = 0; n < 5 * 16; n++) {
    await ocean.advance(calm, 5400);
    state = await ocean.download();
    for (let i = 0; i < coarse.nCells; i++) surface = Math.max(surface, Math.abs(state.eta[i] - eta0[i]));
    for (const v of state.u) fastest = Math.max(fastest, Math.abs(v));
  }
  const { eta, h } = state;
  console.log(`600 m mixed layer ${contrast.toFixed(2)} kg/m³ denser than its 50 m neighbours at 60°S, N=16: over five days |Δη| ≤ ${(100 * surface).toFixed(1)} cm (at the column ${(100 * (eta[centre] - eta0[centre])).toFixed(1)} cm, its steric deficit ${(100 * steric).toFixed(1)} cm), |u| ≤ ${fastest.toFixed(3)} m/s; ${h[centre].toFixed(0)} m deep at the end`);
  assert.ok(contrast > 0.2, `the column is only ${contrast} kg/m³ denser at the surface`);
  assert.ok(surface < 0.15 && surface < 1.5 * steric, `free surface moved ${surface} m against a steric deficit of ${steric} m`);
  assert.ok(fastest < 1, `fastest water ${fastest} m/s`);
  assert.ok(h[centre] > 500, `the neutral layer holds (${h[centre]} m)`);
  assert.equal((await ocean.diagnostics()).oceanLimited, 0);
  await ocean.close();
});

test('the layer left a few metres thick under a 600 m mixed layer at 45°S stays in balance at N=32, where potential vorticity on its smaller edge thickness alone runs it half as fast again', async () => {
  const fine = buildMesh(new Grid(32)), n = fine.nCells, calm = new Float64Array(fine.nEdges);
  const offset = (i) => Math.abs(fine.latCell[i] + 45 * DEG) + Math.abs(fine.lonCell[i]);
  let centre = 0;
  for (let i = 1; i < n; i++) if (offset(i) < offset(centre)) centre = i;
  const fastest = {};
  for (const vorticityCentring of [0.5, 0]) {
    const ocean = await slowOcean(fine, {
      ocean: { everySteps: 1, vorticityCentring },
      surfaceT: Float64Array.from(fine.latCell, (lat) => Math.max(272, 302 - 32 * Math.sin(lat) ** 2)),
      start: (column) => deepen(column, [centre], 600, -0.05),
    });
    console.log(`the remnant at N=32, vorticityCentring ${vorticityCentring}, on the ${ocean.engine} ocean`);
    const { h } = await ocean.download();
    let remnant = 1;
    while (h[remnant * n + centre] <= 5) remnant++;
    assert.ok(h[remnant * n + centre] < 20, `the layer under the mixed layer is ${h[remnant * n + centre]} m thick`);
    fastest[vorticityCentring] = 0;
    for (let step = 0; step < 2 * 32; step++) {
      await ocean.advance(calm, 2700);
      for (const v of (await ocean.download()).u) fastest[vorticityCentring] = Math.max(fastest[vorticityCentring], Math.abs(v));
    }
    await ocean.close();
  }
  console.log(`a 600 m mixed layer at 45°S over a thin remnant, N=32, two days: fastest water ${fastest[0.5].toFixed(2)} m/s with the potential vorticity's thickness at least half the centred one, ${fastest[0].toFixed(2)} m/s on the edge thickness alone`);
  assert.ok(fastest[0.5] < 0.3, `${fastest[0.5]} m/s`);
  assert.ok(fastest[0] > 1.5 * fastest[0.5], `the edge thickness no longer drives the remnant (${fastest[0]} against ${fastest[0.5]} m/s)`);
});

test('the thinned share counts the sea whose mixed layer is below its floor while the water beneath could fill it, and the refill tally takes what the floor drew', () => {
  const ocean = createOcean(mesh, { everySteps: 1 });
  const surfaceT = Float64Array.from(mesh.latCell, (lat) => 275 + 25 * Math.cos(lat) ** 2), ice = new Float64Array(C), flux = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  assert.equal(ocean.diagnostics().oceanThinned, 0);
  const { h, Q, W } = ocean, L = ocean.layers, sea = [];
  for (let i = 0; i < C; i++) if (ocean.cellOcean[i]) sea.push(i);
  const [short, emptied] = sea;
  const move = (i, amount) => { const f = amount / h[i]; h[C + i] += amount; Q[C + i] += Q[i] * f; W[C + i] += W[i] * f; h[i] -= amount; Q[i] -= Q[i] * f; W[i] -= W[i] * f; };
  move(short, 20);
  for (let k = 1; k < L; k++) { const n = k * C + emptied; h[n] = 0.01; }
  h[emptied] = 30;
  const seaArea = sea.reduce((s, i) => s + mesh.areaCell[i], 0);
  assert.ok(Math.abs(ocean.diagnostics().oceanThinned - mesh.areaCell[short] / seaArea) < 1e-12, 'only the column that could refill counts');
  ocean.refilled.water.fill(0); ocean.refilled.heat.fill(0);
  ocean.advance(surfaceT, ice, flux, new Float64Array(mesh.nEdges), 1350);
  assert.ok(ocean.refilled.water[short] > 0, `the floor drew water into the shortened column (${ocean.refilled.water[short]} m)`);
  assert.ok(ocean.refilled.water.every((x) => x >= 0));
  assert.ok(h[short] >= 50 - 1e-9);
  assert.ok(Math.abs(ocean.diagnostics().oceanThinned) < 1e-12, 'the instant refill leaves nothing thinned that could be filled');
});

test('at 10 m floors a 600 m mixed layer beside 10 m ones at 60°S, and the remnant under one at 45°S, run no faster than 1 m/s and reach no clamp', async () => {
  const floors = { mixedDepth: 10, minimumThickness: 10, shallowestMixedDepth: 10 };
  const runs = [];
  const coarse = buildMesh(new Grid(16)), fine = buildMesh(new Grid(32));
  for (const [m, latitude, options, depth, offset, steps, dt] of [
    [coarse, -60, { ...UNIFORM, mixedNeighbourRatio: 0, bottoms: LAYER_BOTTOMS.map((z) => 0.6 * z) }, 600, -0.002, 5 * 16, 5400],
    [fine, -45, { everySteps: 1, vorticityCentring: 0.5 }, 600, -0.05, 2 * 32, 2700],
  ]) {
    const offsetFrom = (i) => Math.abs(m.latCell[i] - latitude * DEG) + Math.abs(m.lonCell[i]);
    let center = 0;
    for (let i = 1; i < m.nCells; i++) if (offsetFrom(i) < offsetFrom(center)) center = i;
    const ocean = await slowOcean(m, {
      ocean: { ...options, ...floors },
      surfaceT: latitude === -60 ? new Float64Array(m.nCells).fill(278) : Float64Array.from(m.latCell, (lat) => Math.max(272, 302 - 32 * Math.sin(lat) ** 2)),
      start: (column) => deepen(column, [center], depth, offset),
    });
    const calm = new Float64Array(m.nEdges);
    let fastest = 0, thinnest = Infinity;
    for (let n = 0; n < steps; n++) {
      await ocean.advance(calm, dt);
      const { h, u } = await ocean.download();
      for (const v of u) fastest = Math.max(fastest, Math.abs(v));
      for (let i = 0; i < m.nCells; i++) if (ocean.cellOcean[i]) thinnest = Math.min(thinnest, h[i]);
    }
    const limited = (await ocean.diagnostics()).oceanLimited, { h } = await ocean.download();
    runs.push({ latitude, engine: ocean.engine, fastest, thinnest, limited, deep: h[center] });
    await ocean.close();
  }
  for (const run of runs) console.log(`a 600 m mixed layer at ${-run.latitude}°S among 10 m ones on the ${run.engine} ocean: fastest water ${run.fastest.toFixed(3)} m/s, the thinnest mixed layer ${run.thinnest.toFixed(1)} m, the deep one ${run.deep.toFixed(0)} m, ${run.limited} clamped edges`);
  for (const run of runs) {
    assert.ok(run.fastest < 1, `${run.latitude}: ${run.fastest} m/s`);
    assert.equal(run.limited, 0);
    assert.ok(run.thinnest > 9.9 && run.thinnest < 11, `${run.latitude}: the thin mixed layers are ${run.thinnest} m`);
  }
});

/*
 * A 50 m mixed layer at 303 K over classes 30, 70, 150, 300 and 3400 m
 * thick on a flat 4000 m sea, each layer given its own velocities, the
 * classes kept on their labels.
 */
const STACK = { everySteps: 1, thermoclineTilt: 0, salinityProfile: () => 34, mixedDepth: 50, restoreTime: 1e30, densities: [1026.0, 1026.1, 1026.2, 1026.3, 1026.95], salinities: [34.9, 34.9, 34.9, 34.9, 34.8], bottoms: [80, 150, 300, 600, 12000] };
function stacked(options = {}, mixedDensity = null) {
  const ocean = createOcean(mesh, { ...STACK, carryMomentum: true, ...options });
  const surfaceT = new Float64Array(C).fill(303), ice = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  if (mixedDensity) {
    const t = labelTemperature(mixedDensity, 34);
    for (let i = 0; i < C; i++) { ocean.Q[i] = ocean.h[i] * t; surfaceT[i] = t; }
    ocean.load(ocean.serialize(), surfaceT, ice);
  }
  for (let k = 0; k < ocean.layers; k++) for (let e = 0; e < E; e++) ocean.u[k * E + e] = (0.3 - 0.04 * k) * Math.sin(0.37 * e + 1.3 * k) + 0.02 * k;
  const start = { h: Float64Array.from(ocean.h), u: Float64Array.from(ocean.u), eta: Float64Array.from(ocean.eta) };
  return { ocean, surfaceT, ice, start };
}
function columnChange({ ocean, start }) {
  const before = columnMomentum(mesh, ocean.layers, start, ocean.D), after = columnMomentum(mesh, ocean.layers, ocean, ocean.D);
  let worst = 0;
  for (let e = 0; e < E; e++) worst = Math.max(worst, Math.abs(after.momentum[e] - before.momentum[e]) / before.scale[e]);
  return worst;
}

test('with carryMomentum every edge\'s column keeps its momentum through the mixed layer\'s retreat, refill and convection: the water a retreat leaves keeps the mixed layer\'s velocity, a refill takes the mean of what the layer now holds and convection the mean of the water it swallows; without it the velocities are left alone', () => {
  const dt = 1350, cells = [...Array(C).keys()], { cellsOnEdge } = mesh, a = (e) => cellsOnEdge[2 * e], b = (e) => cellsOnEdge[2 * e + 1];
  const warmed = (column) => { column.ocean.readSurface(Float64Array.from(column.surfaceT, (t) => t + 0.05), column.ice); column.ocean.mixedLayer(dt); };

  const retreat = stacked({ shallowestMixedDepth: 10, minimumThickness: 10, detrainmentTime: dt, buoyancyMemory: 0 });
  warmed(retreat);
  const { u: u0 } = retreat.start, { h, u } = retreat.ocean;
  assert.ok(cells.every((i) => Math.abs(h[i] - 10) < 1e-9 && Math.abs(h[C + i] - 70) < 1e-9), 'every cell detrains 40 m into the first class');
  let off = 0;
  for (let e = 0; e < E; e++) off = Math.max(off, Math.abs(u[e] - u0[e]), Math.abs(u[E + e] - (40 * u0[e] + 30 * u0[E + e]) / 70));
  const retreatKept = columnChange(retreat);
  console.log(`retreat 50 → 10 m: the class takes (40 u₀ + 30 u₁)/70 and the mixed layer keeps u₀ to ${off.toExponential(2)} m/s; the columns keep their momentum to ${retreatKept.toExponential(2)}`);
  assert.ok(off < 1e-12 && retreatKept < 1e-12);

  const refill = stacked({ minimumThickness: 10 });
  const patch = new Set(cells.filter((i) => Math.abs(mesh.latCell[i] - 30 * DEG) < 15 * DEG && Math.abs(mesh.lonCell[i]) < 20 * DEG));
  const { ocean: r } = refill;
  for (const i of patch) { const f = 47 / r.h[i]; r.h[C + i] += 47; r.Q[C + i] += r.Q[i] * f; r.W[C + i] += r.W[i] * f; r.h[i] -= 47; r.Q[i] -= r.Q[i] * f; r.W[i] -= r.W[i] * f; }
  refill.start.h.set(r.h);
  r.readSurface(refill.surfaceT, refill.ice);
  r.mixedLayer(dt);
  let inside = 0, border = 0, mixedOff = 0;
  for (let e = 0; e < E; e++) {
    const ina = patch.has(a(e)), inb = patch.has(b(e));
    if (ina && inb) { inside++; mixedOff = Math.max(mixedOff, Math.abs(r.u[e] - (3 * refill.start.u[e] + 7 * refill.start.u[E + e]) / 10)); } else if (ina || inb) border++;
  }
  const refillKept = columnChange(refill);
  console.log(`refill 3 → 10 m over ${patch.size} cells: inside the patch the mixed layer takes (3 u₀ + 7 u₁)/10 to ${mixedOff.toExponential(2)} m/s on ${inside} edges; the columns keep their momentum to ${refillKept.toExponential(2)}, the ${border} edges on the patch's border included`);
  assert.ok([...patch].every((i) => Math.abs(r.h[i] - 10) < 1e-9));
  assert.ok(inside > 5 && border > 5 && mixedOff < 1e-12 && refillKept < 1e-12);

  const convection = stacked({ maximumMixedDepth: 300, convectiveRate: 1e5 / 86400 }, 1026.8);
  const { ocean: c, start } = convection;
  const before = Float64Array.from(c.h);
  c.readSurface(convection.surfaceT, convection.ice);
  c.mixedLayer(dt);
  assert.ok(cells.every((i) => c.h[i] > 299.9), `convection deepens every mixed layer to 300 m (${Math.min(...cells.map((i) => c.h[i]))})`);
  let meanOff = 0, moved = 0;
  for (let e = 0; e < E; e++) {
    const reach = Math.max(0.5 * (before[a(e)] + before[b(e)]), 0.5 * (c.h[a(e)] + c.h[b(e)]));
    let top = 0, sum = 0;
    for (let k = 0; k < c.layers && top < reach; k++) { const hk = Math.min(before[k * C + a(e)], before[k * C + b(e)]), part = Math.min(k ? hk : 0.5 * (before[a(e)] + before[b(e)]), reach - top); sum += part * start.u[k * E + e]; top += k ? hk : 0.5 * (before[a(e)] + before[b(e)]); }
    meanOff = Math.max(meanOff, Math.abs(c.u[e] - sum / reach)); moved = Math.max(moved, Math.abs(c.u[e] - start.u[e]));
  }
  const convectionKept = columnChange(convection);
  console.log(`convection 50 → 300 m: the mixed layer takes the mean of the old top 300 m to ${meanOff.toExponential(2)} m/s, moving by up to ${moved.toFixed(3)} m/s; the columns keep their momentum to ${convectionKept.toExponential(2)}`);
  assert.ok(meanOff < 1e-5 && moved > 0.1 && convectionKept < 1e-12);

  const still = stacked({ shallowestMixedDepth: 10, minimumThickness: 10, detrainmentTime: dt, buoyancyMemory: 0, carryMomentum: false });
  warmed(still);
  assert.ok(Math.abs(still.ocean.h[0] - 10) < 1e-9);
  assert.deepEqual(still.ocean.u, still.start.u);
});

/*
 * stacked()'s sea with every third mixed layer cut to 30, 3 and 49.9 m,
 * the water put into the first class, and one mixed-layer step taken.
 */
const CUTS = [30, 3, 49.9];
function cut(options, dt) {
  const column = stacked(options), { ocean } = column, { h, Q, W } = ocean;
  for (let i = 0; i < C; i++) {
    const amount = h[i] - CUTS[i % 3], f = amount / h[i];
    h[C + i] += amount; Q[C + i] += Q[i] * f; W[C + i] += W[i] * f; h[i] -= amount; Q[i] -= Q[i] * f; W[i] -= W[i] * f;
  }
  column.start.h.set(h);
  const sums = (store, i) => { let s = 0, size = 0; for (let k = 0; k < ocean.layers; k++) { s += store[k * C + i]; size += Math.abs(store[k * C + i]); } return [s, size]; };
  const before = [h, Q, W].map((store) => Array.from({ length: C }, (_, i) => sums(store, i)));
  ocean.readSurface(column.surfaceT, column.ice);
  ocean.refilled.water.fill(0);
  ocean.mixedLayer(dt);
  column.kept = Math.max(...[h, Q, W].flatMap((store, n) => Array.from({ length: C }, (_, i) => Math.abs(sums(store, i)[0] - before[n][i][0]) / before[n][i][1])));
  return column;
}

test('under refillRate a mixed layer below its floor deepens by refillRate·dt and never stays thinner than refillFloor, keeping each column\'s water, heat, salt and each edge\'s momentum, and counts as thinned; at refillRate 0 it refills at once, bit for bit as by default', () => {
  const dt = 1350, rate = 24 / 86400, expected = [30 + rate * dt, 5, 50];
  const capped = cut({ refillRate: rate }, dt), { ocean } = capped;
  let off = 0, tally = 0;
  for (let i = 0; i < C; i++) {
    off = Math.max(off, Math.abs(ocean.h[i] - expected[i % 3]));
    tally = Math.max(tally, Math.abs(ocean.refilled.water[i] - (expected[i % 3] - CUTS[i % 3])));
  }
  let area = 0, short = 0;
  for (let i = 0; i < C; i++) { area += mesh.areaCell[i]; if (i % 3 < 2) short += mesh.areaCell[i]; }
  const momentum = columnChange(capped), thinned = ocean.diagnostics().oceanThinned;
  console.log(`refill at ${rate * 86400} m/day over ${dt} s: 30 m → ${ocean.h[0].toFixed(6)}, 3 m → ${ocean.h[1].toFixed(6)}, 49.9 m → ${ocean.h[2].toFixed(6)} m, each to ${off.toExponential(2)} m; tallied to ${tally.toExponential(2)} m; columns keep water, heat and salt to ${capped.kept.toExponential(2)} and each edge its momentum to ${momentum.toExponential(2)}; thinned ${thinned.toFixed(4)} of the sea against ${(short / area).toFixed(4)}`);
  assert.ok(off < 1e-12 && tally < 1e-12, `the mixed layers are off by ${off} m, the tally by ${tally} m`);
  assert.ok(capped.kept < 1e-12 && momentum < 1e-12, `the columns changed by ${capped.kept}, an edge's momentum by ${momentum}`);
  assert.ok(Math.abs(thinned - short / area) < 1e-12, `thinned ${thinned} against ${short / area}`);

  const instant = cut({ refillRate: 0, refillFloor: 7 }, dt), plain = cut({}, dt);
  assert.ok(Array.from({ length: C }, (_, i) => instant.ocean.h[i]).every((x) => Math.abs(x - 50) < 1e-9), 'at refillRate 0 every mixed layer is back on its floor');
  for (const name of ['h', 'u', 'Q', 'W', 'eta']) assert.deepEqual(instant.ocean[name], plain.ocean[name], `${name} at refillRate 0 is the default's`);
  assert.equal(instant.ocean.diagnostics().oceanThinned, 0);
});

/*
 * The shear across the mixed layer's base as the CPU ocean takes it:
 * on each sea edge the mixed layer's velocity less that of the first
 * class at least THIN thick there, reconstructed to the cells.
 */
function baseShear(ocean) {
  const { h, u, layers: L } = ocean, { cellsOnEdge } = mesh, du = new Float64Array(E);
  for (let e = 0; e < E; e++) {
    const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
    for (let k = 1; k < L; k++) if (Math.min(h[k * C + a], h[k * C + b]) >= THIN) { du[e] = u[e] - u[k * E + e]; break; }
  }
  const v = cellVector(mesh, du);
  return Float64Array.from({ length: C }, (_, i) => v[3 * i] ** 2 + v[3 * i + 1] ** 2 + v[3 * i + 2] ** 2);
}

/*
 * A 20 m mixed layer moving east at 0.5 m/s over a still 60 m class,
 * its buoyancy step set cell by cell so that the bulk Richardson number
 * Δb·h/|Δu|² is 0.4 wherever the shear is at least 0.2 m²/s², and far
 * above Ri_c elsewhere. The two waters differ along the one direction
 * in which the quadratic equation of state is linear, so mixing them
 * keeps the buoyancy step in proportion to the mixed layer's share.
 */
const RIC = 0.65, PWP = { everySteps: 1, thermoclineTilt: 0, salinityProfile: () => 35, densities: [1026.5, 1027.5], salinities: [35, 35], bottoms: [80, 12000], restoreTime: 1e30, carryMomentum: true, shearEntrainment: RIC, minimumThickness: 5, shallowestMixedDepth: 5, mixedDepth: 20 };
function pwpColumns(options = {}) {
  const { a0, b0, lambda1, lambda2, nu, t0, s0 } = SEAWATER, g = 9.81, rho0 = 1025;
  const ocean = createOcean(mesh, { ...PWP, ...options }), { h, Q, W, u, eta, D } = ocean;
  const t1 = labelTemperature(1026.5, 35), ta = t1 - t0, sa = 35 - s0;
  const rhoT = -a0 * (1 + lambda1 * ta) - nu * sa, rhoS = b0 * (1 - lambda2 * sa) - nu * ta;
  const r = (-nu + Math.sqrt(nu * nu - a0 * lambda1 * b0 * lambda2)) / (a0 * lambda1);
  const east = zonalWindOnEdges(mesh, () => 0.5), v = cellVector(mesh, east), shear = Float64Array.from({ length: C }, (_, i) => v[3 * i] ** 2 + v[3 * i + 1] ** 2 + v[3 * i + 2] ** 2);
  const surfaceT = new Float64Array(C), ice = new Float64Array(C);
  ocean.initialize(new Float64Array(C).fill(t1), ice);
  for (let i = 0; i < C; i++) {
    const db = shear[i] >= 0.2 ? 0.4 * (shear[i] + 1e-12) / 20 : 40 * (shear[i] + 1e-12) / 20 + 0.01, x = db * rho0 / g / (r * rhoT + rhoS);
    const tm = t1 - r * x, sm = 35 - x;
    h[i] = 20; h[C + i] = 60; h[2 * C + i] = D[i] - 80;
    Q[i] = 20 * tm; W[i] = 20 * sm; Q[C + i] = 60 * t1; W[C + i] = 60 * 35;
    Q[2 * C + i] = h[2 * C + i] * labelTemperature(1027.5, 35); W[2 * C + i] = h[2 * C + i] * 35;
    eta[i] = 0; surfaceT[i] = tm;
  }
  ocean.load(ocean.serialize(), surfaceT, ice);
  u.fill(0);
  u.set(east);
  ocean.setStress(new Float64Array(E), ice);
  ocean.readSurface(surfaceT, ice);
  const entraining = new Set([...Array(C).keys()].filter((i) => shear[i] >= 0.2));
  const inner = [...entraining].filter((i) => [...Array(mesh.nEdgesOnCell[i]).keys()].every((m) => entraining.has(mesh.cellsOnCell[mesh.maxEdges * i + m])));
  return { ocean, surfaceT, ice, shear, inner, rib: (i) => g * (1026.5 - seawaterDensity(ocean.Q[i] / ocean.h[i], ocean.W[i] / ocean.h[i])) / rho0 * ocean.h[i] / (baseShear(ocean)[i] + 1e-12) };
}
const columnSums = (ocean) => [ocean.h, ocean.Q, ocean.W].map((store) => Array.from({ length: C }, (_, i) => { let s = 0, size = 0; for (let k = 0; k < ocean.layers; k++) { s += store[k * C + i]; size += Math.abs(store[k * C + i]); } return [s, size]; }));

test('under shearEntrainment a mixed layer whose bulk Richardson number is below Ri_c takes in the class beneath to h·√(Ri_c/Ri_b) in one step, which the exchange of momentum brings to Ri_c, so that a second step takes nothing; columns keep their water, heat and salt and edges their momentum', () => {
  const column = pwpColumns(), { ocean, inner } = column, dt = 1350;
  const before = inner.map((i) => column.rib(i)), sums = columnSums(ocean), start = { h: Float64Array.from(ocean.h), u: Float64Array.from(ocean.u), eta: Float64Array.from(ocean.eta) };
  ocean.mixedLayer(dt);
  const expected = 20 * Math.sqrt(RIC / 0.4);
  let closed = 0, ideal = 0, ri = 0, taken = 0;
  inner.forEach((i, n) => {
    closed = Math.max(closed, Math.abs(ocean.h[i] - 20 * Math.sqrt(RIC / before[n])));
    ideal = Math.max(ideal, Math.abs(ocean.h[i] - expected));
    ri = Math.max(ri, Math.abs(column.rib(i) - RIC));
    taken = Math.max(taken, Math.abs(ocean.sheared[i] - (ocean.h[i] - 20)));
  });
  const after = columnSums(ocean);
  const kept = Math.max(...after.flatMap((cells, n) => cells.map(([s], i) => Math.abs(s - sums[n][i][0]) / sums[n][i][1])));
  const before2 = Float64Array.from(ocean.h);
  ocean.mixedLayer(dt);
  const second = Math.max(...inner.map((i) => ocean.h[i] - before2[i]));
  const momentum = (() => { const a = columnMomentum(mesh, ocean.layers, start, ocean.D), b = columnMomentum(mesh, ocean.layers, ocean, ocean.D); let worst = 0; for (let e = 0; e < E; e++) worst = Math.max(worst, Math.abs(b.momentum[e] - a.momentum[e]) / a.scale[e]); return worst; })();
  console.log(`PWP: ${inner.length} cells at Ri_b 0.4 deepen 20 → ${ocean.h[inner[0]].toFixed(6)} m, to the closed form h₀√(Ri_c/Ri_b) within ${closed.toExponential(2)} m (${expected.toFixed(6)} m within ${ideal.toExponential(2)}); Ri_b after the exchange of momentum ${RIC} within ${ri.toExponential(2)}; a second step takes ${second.toExponential(2)} m; columns keep water, heat and salt to ${kept.toExponential(2)}, edges their momentum to ${momentum.toExponential(2)} over both steps; sheared tallies the entrainment to ${taken.toExponential(2)} m`);
  assert.ok(inner.length > 300, `${inner.length} cells`);
  assert.ok(closed < 1e-12 && ideal < 1e-9, `h off the closed form by ${closed} m, off ${expected} by ${ideal} m`);
  assert.ok(ri < 1e-9, `Ri_b after the exchange off Ri_c by ${ri}`);
  assert.ok(second < 1e-9, `a second step took ${second} m`);
  assert.ok(kept < 1e-12 && momentum < 1e-12 && taken < 1e-12, `columns ${kept}, momentum ${momentum}, tally ${taken}`);
});

test('under shearEntrainment a heated mixed layer whose bulk Richardson number is below Ri_c ends the step at the shear\'s depth h·√(Ri_c/Ri_b), its retreat in the same step held there however shallow its Monin–Obukhov depth', () => {
  const dt = 1350, net = 150;
  const column = pwpColumns({ buoyancyMemory: 0, detrainmentTime: dt }), { ocean, surfaceT, ice, inner } = column;
  for (let i = 0; i < C; i++) surfaceT[i] += dt * net / ocean.capacity[i];
  ocean.setStress(zonalWindOnEdges(mesh, () => 0.002), ice);
  ocean.readSurface(surfaceT, ice);
  const before = inner.map((i) => column.rib(i));
  ocean.mixedLayer(dt);
  let held = 0, worst = 0, deepest = 0;
  inner.forEach((i, n) => {
    if (!(ocean.retreatDepth[i] > 0)) return;
    held++;
    worst = Math.max(worst, Math.abs(ocean.h[i] - 20 * Math.sqrt(RIC / before[n])));
    deepest = Math.max(deepest, ocean.retreatDepth[i]);
  });
  console.log(`${held} heated cells at Ri_b about 0.4 with Monin–Obukhov depths of at most ${deepest.toFixed(3)} m end the step at the shear's depth within ${worst.toExponential(2)} m`);
  assert.ok(held > 300 && deepest < 5, `${held} heated cells, L up to ${deepest} m`);
  assert.ok(worst < 1e-9, `the heated layer ended ${worst} m off h·√(Ri_c/Ri_b)`);
});

test('shearEntrainment without carryMomentum throws: the base\'s Richardson number would never rise', () => {
  assert.throws(() => createOcean(mesh, { shearEntrainment: RIC }), /carryMomentum/);
  assert.throws(() => pwpColumns({ carryMomentum: false }), /carryMomentum/);
});

test('under shearEntrainment a 5 m mixed layer at 45°N over a stratification of N = 0.01 /s, the wind on it alone switched on and no interfacial drag, as the theory has none, deepens in half an inertial period to Pollard, Rhines and Thompson\'s 1.7u*/√(Nf), within a quarter, and without it stays 5 m deep', () => {
  const step = 0.1, N2 = 1e-4, thick = step * 9.81 / (1025 * N2), n = 10, tau = 0.1, f45 = 2 * 7.2921e-5 * Math.SQRT1_2, steps = 24, dt = Math.PI / f45 / steps;
  const densities = [...Array(n).keys()].map((k) => 1025.5 + step * k).concat([1027.5]);
  const run = (shearEntrainment) => {
    const ocean = createOcean(mesh, { everySteps: 1, thermoclineTilt: 0, salinityProfile: () => 35, densities, salinities: densities.map(() => 35), bottoms: densities.map((_, k) => 100 + 100 * k), restoreTime: 1e30, carryMomentum: true, shearEntrainment, interfacialDrag: 0, momentumDepth: 5, minimumThickness: 5, shallowestMixedDepth: 5, stirring: 0, mixedDepth: 5 });
    const { h, Q, W, eta, D } = ocean, L = ocean.layers, tm = labelTemperature(1025.5 - step / 2, 35);
    const surfaceT = new Float64Array(C).fill(tm), ice = new Float64Array(C), flux = new Float64Array(C);
    ocean.initialize(surfaceT, ice);
    for (let i = 0; i < C; i++) {
      h[i] = 5; Q[i] = 5 * tm; W[i] = 5 * 35; eta[i] = 0;
      for (let k = 1; k < L; k++) { h[k * C + i] = k < L - 1 ? thick : D[i] - 5 - n * thick; Q[k * C + i] = h[k * C + i] * labelTemperature(densities[k - 1], 35); W[k * C + i] = h[k * C + i] * 35; }
    }
    ocean.load(ocean.serialize(), surfaceT, ice);
    const stress = zonalWindOnEdges(mesh, () => tau);
    for (let s = 0; s < steps; s++) ocean.advance(surfaceT, ice, flux, stress, dt);
    const band = [...Array(C).keys()].filter((i) => Math.abs(mesh.latCell[i] / DEG - 45) < 5).map((i) => ocean.h[i]).sort((a, b) => a - b);
    return band[band.length >> 1];
  };
  const prt = 1.7 * Math.sqrt(tau / 1025) / Math.sqrt(Math.sqrt(N2) * f45), deepened = run(RIC), still = run(0);
  console.log(`PRT: under a ${tau} N/m² wind over N = ${Math.sqrt(N2)} /s, half an inertial period at 45°N takes the 5 m mixed layer to ${deepened.toFixed(2)} m against 1.7u*/√(Nf) = ${prt.toFixed(2)} m; without shearEntrainment it is ${still.toFixed(2)} m`);
  assert.ok(Math.abs(deepened / prt - 1) < 0.25, `${deepened} m against ${prt} m`);
  assert.ok(Math.abs(still - 5) < 0.5, `without shear entrainment the layer is ${still} m`);
});

test('with shearEntrainment, ekmanDepth, retreatStirring and langmuir 0 and stressThickness PV_FLOOR the ocean steps as without the options, bit for bit, over 48 steps from the analytic start', () => {
  const run = (options) => {
    const { ocean, surfaceT, ice, flux } = uniformOcean({ penetration: 'IB', carryMomentum: true, minimumThickness: 10, shallowestMixedDepth: 10, refillRate: 24 / 86400, ...options });
    const stress = zonalWindOnEdges(mesh, (lat) => -0.1 * Math.cos(3 * lat));
    for (let n = 0; n < 48; n++) { for (let i = 0; i < C; i++) surfaceT[i] = 278 + 2 * Math.cos(mesh.latCell[i]) - 0.01 * n; ocean.advance(surfaceT, ice, flux, stress, 1350); }
    return ocean;
  };
  const off = run({ shearEntrainment: 0, ekmanDepth: 0, momentumFloor: 10, stressThickness: 20, retreatStirring: 0, langmuir: 0 }), plain = run({});
  for (const name of ['h', 'u', 'Q', 'W', 'eta']) assert.deepEqual(off[name], plain[name], `${name} with the options off is the default's`);
});

test('under retreatStirring and langmuir a heated mixed layer retreats toward C_r·ε·u*³/(−B), unfaded, with ε the Langmuir factor at the depth the step began with, between 1 and its cap', () => {
  const dt = 1350, net = 150, cap = 2;
  const stress = zonalWindOnEdges(mesh, (lat) => 0.01 + 0.2 * Math.cos(lat) ** 4);
  const tauCell = cellVector(mesh, stress);
  const run = (mixedDepth, options) => {
    const { ocean, surfaceT, ice } = uniformOcean({ mixedDepth, minimumThickness: 5, shallowestMixedDepth: 5, buoyancyMemory: 0, ...options });
    for (let i = 0; i < C; i++) surfaceT[i] += dt * net / ocean.capacity[i];
    ocean.setStress(stress, ice);
    ocean.readSurface(surfaceT, ice);
    const began = Float64Array.from(ocean.h.subarray(0, C)), alpha = Float64Array.from(began, (h, i) => thermalExpansion(surfaceT[i], ocean.W[i] / h));
    ocean.mixedLayer(dt);
    return { L: Float64Array.from(ocean.retreatDepth), began, alpha, rhoCp: ocean.rhoCp };
  };
  let cells = 0, worst = 0, closed = 0, least = Infinity, most = 0, capped = 0;
  for (const depth of [6, 80]) for (let faded = run(depth, {}), constant = run(depth, { retreatStirring: 2.5 }), langmuir = run(depth, { retreatStirring: 2.5, langmuir: cap }), i = 0; i < C; i++) {
    if (!(faded.L[i] > 0)) continue;
    cells++;
    const h = faded.began[i], tau = Math.hypot(tauCell[3 * i], tauCell[3 * i + 1], tauCell[3 * i + 2]), ustar = Math.sqrt(tau / RHO);
    const scale = faded.L[i] / (2 * 0.8 * Math.exp(-h / 100)), epsilon = langmuirFactor(ustar, h, cap);
    assert.equal(constant.began[i], h);
    worst = Math.max(worst, Math.abs(constant.L[i] / (2.5 * scale) - 1), Math.abs(langmuir.L[i] / (2.5 * epsilon * scale) - 1));
    closed = Math.max(closed, Math.abs(constant.L[i] / (2.5 * ustar ** 3 * faded.rhoCp / (9.81 * faded.alpha[i] * net)) - 1));
    const ratio = langmuir.L[i] / constant.L[i];
    least = Math.min(least, ratio); most = Math.max(most, ratio);
    if (ratio > cap - 1e-12) capped++;
  }
  console.log(`${cells} heated cells, 6 and 80 m deep: L at 2.5 and at 2.5ε within ${worst.toExponential(2)} of the faded rule's u*³/(−B) times their constants, and the 2.5 rule within ${closed.toExponential(2)} of 2.5u*³ρc_p/(gα·${net} W/m²); ε from ${least.toFixed(3)} to ${most.toFixed(3)}, ${capped} cells at the cap ${cap}`);
  assert.ok(cells > 600, `${cells} cells retreated`);
  assert.ok(worst < 1e-12, `the retreat depth is ${worst} off its closed form`);
  assert.ok(closed < 1e-9, `the 2.5 rule is ${closed} off the flux's Monin–Obukhov depth`);
  assert.ok(least >= 1 && most <= cap * (1 + 1e-15) && least < 1.5 && capped > 0, `ε from ${least} to ${most}`);

  const { ocean, surfaceT, ice } = uniformOcean({ mixedDepth: 30, minimumThickness: 5, shallowestMixedDepth: 5, buoyancyMemory: 0, carryMomentum: true, shearEntrainment: 0.65, ekmanDepth: 0.4, stressThickness: THIN, retreatStirring: 2.5, langmuir: cap });
  for (let e = 0; e < E; e++) {
    if (Math.abs(mesh.latEdge[e]) > 30 * DEG) continue;
    const x = mesh.xEdge[3 * e], y = mesh.xEdge[3 * e + 1], r = Math.hypot(x, y) || 1e-12;
    ocean.u[e] += 0.4 * ((-y / r) * mesh.nEdge[3 * e] + (x / r) * mesh.nEdge[3 * e + 1]);
  }
  for (let i = 0; i < C; i++) surfaceT[i] += dt * net / ocean.capacity[i];
  ocean.setStress(stress, ice);
  ocean.readSurface(surfaceT, ice);
  const sums = columnSums(ocean), start = { h: Float64Array.from(ocean.h), u: Float64Array.from(ocean.u), eta: Float64Array.from(ocean.eta) };
  ocean.mixedLayer(dt);
  const after = columnSums(ocean), books = Math.max(...after.flatMap((cells, n) => cells.map(([s], i) => Math.abs(s - sums[n][i][0]) / sums[n][i][1])));
  const a = columnMomentum(mesh, ocean.layers, start, ocean.D), b = columnMomentum(mesh, ocean.layers, ocean, ocean.D);
  let momentum = 0, retreated = 0, entrained = 0;
  for (let e = 0; e < E; e++) momentum = Math.max(momentum, Math.abs(b.momentum[e] - a.momentum[e]) / Math.max(a.scale[e], 1e-30));
  for (let i = 0; i < C; i++) { if (ocean.h[i] < start.h[i] - 1e-9) retreated++; if (ocean.sheared[i] > 0) entrained++; }
  console.log(`with the whole package, ${retreated} cells retreat and ${entrained} entrain by shear in one step; columns keep water, heat and salt within ${books.toExponential(2)} and edges their momentum within ${momentum.toExponential(2)}`);
  assert.ok(retreated > 100 && entrained > 100, `${retreated} retreated, ${entrained} entrained`);
  assert.ok(books < 1e-12 && momentum < 1e-12, `columns ${books}, momentum ${momentum}`);
});

test('under retreatEkman the heated mixed layer retreats toward no more than retreatEkman·u*/|f|: a 300 m layer under a strong wind and a weak gain lands on that Ekman bound, while the unbounded rule holds it deeper, and wherever the bound lies deeper than C_r·ε·u*³/(−B) the depth keeps that closed form', () => {
  const dt = 1350, net = 50, c = 0.7;
  const stress = zonalWindOnEdges(mesh, (lat) => 0.01 + 0.2 * Math.cos(lat) ** 4);
  const tauCell = cellVector(mesh, stress);
  const run = (mixedDepth, options) => {
    const { ocean, surfaceT, ice } = uniformOcean({ mixedDepth, minimumThickness: 5, shallowestMixedDepth: 5, buoyancyMemory: 0, detrainmentTime: dt, retreatStirring: 2.5, langmuir: 2, ...options });
    for (let i = 0; i < C; i++) surfaceT[i] += dt * net / ocean.capacity[i];
    ocean.setStress(stress, ice);
    ocean.readSurface(surfaceT, ice);
    ocean.mixedLayer(dt);
    return { L: Float64Array.from(ocean.retreatDepth), h: Float64Array.from(ocean.h.subarray(0, C)) };
  };
  const lines = [];
  for (const depth of [6, 300]) {
    const free = run(depth, {}), bounded = run(depth, { retreatEkman: c });
    let heated = 0, binding = 0, worst = 0, landed = 0, deeper = 0;
    for (let i = 0; i < C; i++) {
      if (!(free.L[i] > 0)) continue;
      heated++;
      const ustar = Math.sqrt(Math.hypot(tauCell[3 * i], tauCell[3 * i + 1], tauCell[3 * i + 2]) / RHO), limit = mesh.fCell[i] !== 0 ? c * ustar / Math.abs(mesh.fCell[i]) : Infinity;
      worst = Math.max(worst, Math.abs(bounded.L[i] / Math.min(free.L[i], limit) - 1));
      if (limit >= free.L[i]) continue;
      binding++;
      if (limit < depth && Math.abs(bounded.h[i] - Math.max(5, limit)) < 1e-9) landed++;
      if (free.h[i] > bounded.h[i] + 1) deeper++;
    }
    lines.push(`${depth} m: ${heated} heated cells, ${binding} with the bound below C_r·ε·u*³/(−B), ${landed} landing on it, ${deeper} held over a meter deeper without it`);
    assert.ok(worst < 1e-12, `the retreat depth is ${worst} off min(C_r·ε·u*³/(−B), ${c}u*/|f|) at ${depth} m`);
    assert.ok(binding > 300 && heated - binding > 100, `${binding} of ${heated} heated cells bound at ${depth} m`);
    if (depth === 300) assert.ok(landed > 300 && deeper > 300, `${landed} landed on the bound, ${deeper} held deeper without it`);
  }
  console.log(`the Ekman bound ${c}u*/|f| under ${net} W/m²: ${lines.join('; ')}`);
});

test('under langmuir the Langmuir factor\'s excess over 1 scales with the open water\'s share of the cell, so a melting layer under full pack ice retreats toward C_r·u*³/(−B) and one under partial cover toward C_r·(1 + (ε − 1)(1 − a))·u*³/(−B)', () => {
  const dt = 1350, cap = 2;
  const stress = zonalWindOnEdges(mesh, (lat) => 0.05 + 0.2 * Math.cos(lat) ** 4);
  const ice = Float64Array.from({ length: C }, (_, i) => (Math.abs(mesh.latCell[i]) > 40 * DEG ? 1 : 0));
  const concentration = Float64Array.from({ length: C }, (_, i) => (i % 3 === 0 ? 0 : 0.1 + 0.8 * (i % 7) / 6));
  const run = (options) => {
    const ocean = createOcean(mesh, { ...UNIFORM, mixedDepth: 30, minimumThickness: 5, shallowestMixedDepth: 5, buoyancyMemory: 0, retreatStirring: 2.5, ...options });
    const surfaceT = new Float64Array(C).fill(278);
    ocean.initialize(surfaceT, ice, { concentration });
    for (let i = 0; i < C; i++) if (!ice[i]) surfaceT[i] += dt * 150 / ocean.capacity[i];
    ocean.setStress(stress, ice, concentration);
    ocean.readSurface(surfaceT, ice, concentration);
    const began = Float64Array.from(ocean.h.subarray(0, C));
    ocean.mixedLayer(dt);
    return { L: Float64Array.from(ocean.retreatDepth), began, tauCell: cellVector(mesh, ocean.stress) };
  };
  const plain = run({}), waved = run({ langmuir: cap });
  let full = 0, partial = 0, open = 0, damped = 0, worst = 0;
  for (let i = 0; i < C; i++) {
    if (!(plain.L[i] > 0)) continue;
    const cover = ice[i] > 0 ? (concentration[i] > 0 ? concentration[i] : 1) : 0, t = plain.tauCell;
    const epsilon = langmuirFactor(Math.sqrt(Math.hypot(t[3 * i], t[3 * i + 1], t[3 * i + 2]) / RHO), plain.began[i], cap);
    worst = Math.max(worst, Math.abs(waved.L[i] / plain.L[i] - (1 + (epsilon - 1) * (1 - cover))));
    if (cover === 1) full++; else if (cover > 0) partial++; else open++;
    if (cover > 0 && epsilon > 1.05) damped++;
  }
  console.log(`${full} heated cells under full ice, ${partial} under partial cover and ${open} in open water: the Langmuir retreat over the plain one within ${worst.toExponential(2)} of 1 + (ε − 1)(1 − a), ${damped} iced cells with ε above 1.05`);
  assert.ok(full > 30 && partial > 50 && open > 100 && damped > 50, `${full} full, ${partial} partial, ${open} open, ${damped} damped`);
  assert.ok(worst < 1e-12, `the Langmuir factor is ${worst} off its open-water share`);
});

/*
 * CVMix's cvmix_kpp_ustokes_SL_model, the surface-layer Stokes drift of
 * the theory wave, written out with erfc from its Maclaurin series.
 */
function publishedStokes(u10, h) {
  const erfc = (x) => { let term = x, sum = x; for (let n = 1; n < 80; n++) { term *= -x * x / n; sum += term / (2 * n + 1); } return 1 - 2 / Math.sqrt(Math.PI) * sum; };
  const us = 0.0162 * u10, hm0 = 0.0246 * u10 ** 2, fp = 0.877 * 9.81 / (2 * Math.PI * 1.075 * u10), fm = 1.296 * fp;
  const vstokes = 0.125 * Math.PI * 0.667 * fm * hm0 ** 2, kphil = 0.176 * us / vstokes, kstar = 2.56 * kphil, z0 = 0.2 * h;
  const r1 = (0.151 / kphil / z0 - 0.84) * (1 - Math.exp(-2 * kphil * z0));
  const r2 = -(0.84 + 0.0591 / kphil / z0) * Math.sqrt(2 * Math.PI * kphil * z0) * erfc(Math.sqrt(2 * kphil * z0));
  const r3 = (0.0632 / kstar / z0 + 0.125) * (1 - Math.exp(-2 * kstar * z0));
  const r4 = (0.125 + 0.0946 / kstar / z0) * Math.sqrt(2 * Math.PI * kstar * z0) * erfc(Math.sqrt(2 * kstar * z0));
  return us * (0.715 + r1 + r2 + r3 + r4);
}

test('the surface-layer Langmuir number of a 10 m/s theory wave over a 10 m layer is the published formula\'s, and the factor is 1 in a calm and never below 1 nor above its cap', () => {
  const ustar = 10 * FRICTION_PER_WIND, published = publishedStokes(10, 10), model = stokesSurfaceLayer(10, 10);
  const la = Math.sqrt(ustar / model), laPublished = Math.sqrt(ustar / published);
  console.log(`U10 10 m/s, h 10 m: u_SL ${model.toFixed(6)} m/s against the formula's ${published.toFixed(6)}, La_SL ${la.toFixed(4)} against ${laPublished.toFixed(4)}, ε ${langmuirFactor(ustar, 10, 5).toFixed(4)}`);
  assert.ok(Math.abs(la / laPublished - 1) < 0.01 && Math.abs(model / 0.07161749 - 1) < 0.01, `La_SL ${la} against ${laPublished}`);
  assert.ok(Math.abs(langmuirFactor(ustar, 10, 5) - Math.sqrt(1 + la ** -2 / 1.5 ** 2 + la ** -4 / 5.4 ** 4)) < 1e-12);
  assert.equal(langmuirFactor(0, 10, 2), 1);
  let previous = Infinity;
  for (const h of [1, 3, 10, 30, 100, 300, 600]) {
    const at = langmuirFactor(ustar, h, 5);
    assert.ok(at <= previous + 1e-12, `the factor grows from ${previous} to ${at} at ${h} m`);
    previous = at;
    for (const u10 of [0.002, 0.1, 1, 3, 10, 20, 40]) { const e = langmuirFactor(u10 * FRICTION_PER_WIND, h, 2); assert.ok(e >= 1 && e <= 2 && Number.isFinite(e), `ε ${e} at ${u10} m/s over ${h} m`); }
  }
});
