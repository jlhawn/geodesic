import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createOcean, LAYER_BOTTOMS } from '../js/ocean/layered.module.js';
import { seawaterDensity, thermalExpansion, labelTemperature } from '../js/ocean/seawater.module.js';
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
    let centre = 0;
    for (let i = 1; i < m.nCells; i++) if (offsetFrom(i) < offsetFrom(centre)) centre = i;
    const ocean = await slowOcean(m, {
      ocean: { ...options, ...floors },
      surfaceT: latitude === -60 ? new Float64Array(m.nCells).fill(278) : Float64Array.from(m.latCell, (lat) => Math.max(272, 302 - 32 * Math.sin(lat) ** 2)),
      start: (column) => deepen(column, [centre], depth, offset),
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
    runs.push({ latitude, engine: ocean.engine, fastest, thinnest, limited, deep: h[centre] });
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
