import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createOcean, LAYER_BOTTOMS } from '../js/ocean/layered.module.js';
import { seawaterDensity, thermalExpansion } from '../js/ocean/seawater.module.js';
import { RHO_AIR, DRAG, RHO, DEG, mesh, C, zonalWindOnEdges, UNIFORM, uniformOcean, deepen, slowOcean } from './helpers/layered.mjs';

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
