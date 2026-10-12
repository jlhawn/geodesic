import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createModel } from '../js/model.module.js';
import { createOcean, EPS } from '../js/ocean/layered.module.js';
import { seawaterDensity, thermalExpansion, halineContraction } from '../js/ocean/seawater.module.js';
import { createSeaIce, FREEZING_POINT } from '../js/physics/ice.module.js';

const mesh = buildMesh(new Grid(2)), C = mesh.nCells, E = mesh.nEdges;
const ICE_DENSITY = 917, ICE_SALINITY = 5, HOUR = 3600, GRAVITY = 9.81;

function frozenColumn(options = {}) {
  const ocean = createOcean(mesh, { everySteps: 1, thermoclineTilt: 0, salinityProfile: () => 34, mixedDepth: 50, convectiveRate: 0, stirring: 0, ...options });
  const surfaceT = new Float64Array(C).fill(FREEZING_POINT + 0.001), ice = new Float64Array(C);
  const sea = createSeaIce(mesh, { heatCapacity: ocean.capacity });
  ocean.useSeaIce(sea.concentration);
  ocean.initialize(surfaceT, ice);
  return { ocean, sea, surfaceT, ice, flux: new Float64Array(C), calm: new Float64Array(E) };
}

test('new ice over a small part of a cell salts the mixed layer by its volume, cover times thickness, not by its thickness', () => {
  const { ocean, sea, surfaceT, ice, flux, calm } = frozenColumn();
  const i = 0, S0 = ocean.W[i] / ocean.h[i], depth = ocean.h[i], salt = Float64Array.from(ocean.W.subarray(0, C));
  flux[i] = -400;
  for (let n = 0; n < 2; n++) sea.update(surfaceT, ice, flux, i, HOUR);
  const cover = sea.concentration[i], volume = cover * ice[i];
  ocean.advance(surfaceT, ice, sea.oceanFlux, calm, HOUR, sea.concentration);
  const gained = ocean.W[i] / ocean.h[i] - S0;
  const byVolume = (S0 - ICE_SALINITY) * ICE_DENSITY * volume / 1000 / depth, byThickness = (S0 - ICE_SALINITY) * ICE_DENSITY * ice[i] / 1000 / depth;
  console.log(`ice ${ice[i].toFixed(3)} m thick over ${(100 * cover).toFixed(1)}% of the cell salts the ${depth} m mixed layer by ${gained.toFixed(5)} psu; by its thickness it would be ${byThickness.toFixed(4)}`);
  assert.ok(cover > 0.02 && cover < 0.04 && ice[i] > 0.29, `the ice covers ${cover} of the cell, ${ice[i]} m thick`);
  assert.ok(Math.abs(gained - byVolume) < 1e-6 * byVolume, `salinity gained ${gained} against ${byVolume} by volume`);
  assert.ok(gained > 0.0044 && gained < 0.0048 && byThickness > 30 * gained, `${gained} psu by volume, ${byThickness} by thickness`);
  for (let j = 1; j < C; j++) assert.ok(Math.abs(ocean.W[j] - salt[j]) < 1e-6, `cell ${j} without ice keeps its salt: ${salt[j]} -> ${ocean.W[j]}`);
});

test('the brine of new ice over a small part of a cell erodes the denser class under the mixed layer at the buoyancy of its volume, cover times thickness, not of its thickness', () => {
  const { ocean, sea, surfaceT, ice, flux, calm } = frozenColumn({ buoyancyMemory: 0, convectiveRate: 0.1 });
  const i = 0, rho = ocean.densities, depth = ocean.h[i], s0 = ocean.W[i] / depth, previousT0 = ocean.T0[i], rm = seawaterDensity(FREEZING_POINT, s0);
  let k = 1;
  while (ocean.h[k * C + i] <= EPS) k++;
  const below = ocean.h[k * C + i], others = Float64Array.from(ocean.h.subarray(0, C));
  flux[i] = -400;
  for (let n = 0; n < 2; n++) sea.update(surfaceT, ice, flux, i, HOUR);
  const cover = sea.concentration[i], volume = cover * ice[i];
  ocean.advance(surfaceT, ice, sea.oceanFlux, calm, HOUR, sea.concentration);
  const deepened = ocean.h[i] - depth, eroded = below - ocean.h[k * C + i];
  const buoyancy = (grown) => GRAVITY * (thermalExpansion(FREEZING_POINT, s0) * (previousT0 - FREEZING_POINT) * depth + halineContraction(FREEZING_POINT, s0) * (s0 - ICE_SALINITY) * ICE_DENSITY * grown / 1000) / HOUR;
  const N2 = GRAVITY * (rho[k] + 0.5 * (rho[k + 1] - rho[k]) - rm) / (rho[0] * below);
  const byVolume = buoyancy(volume) * HOUR / (depth * N2), byThickness = buoyancy(ice[i]) * HOUR / (depth * N2);
  console.log(`ice ${ice[i].toFixed(3)} m thick over ${(100 * cover).toFixed(1)}% of the cell deepens the ${depth} m mixed layer into the ${rho[k]} class beneath by ${deepened.toFixed(4)} m, where the buoyancy of its volume gives ${byVolume.toFixed(4)} m and of its thickness ${byThickness.toFixed(1)} m`);
  assert.ok(cover > 0.02 && cover < 0.04 && ice[i] > 0.29, `the ice covers ${cover} of the cell, ${ice[i]} m thick`);
  assert.ok(rho[k] > rm && rm > rho[k] - 0.5 * (rho[k] - rho[k - 1]) && below > 100, `the mixed layer at ${rm} lies in the lighter half of the ${below} m class labelled ${rho[k]}`);
  assert.ok(Math.abs(deepened - byVolume) < 1e-6 * byVolume, `the mixed layer deepened by ${deepened} m against ${byVolume} by volume`);
  assert.ok(Math.abs(eroded - deepened) < 1e-6 * deepened, `the class beneath gave ${eroded} m for the mixed layer's ${deepened}`);
  assert.ok(deepened > 10 && byThickness > 30 * deepened, `${deepened} m by volume, ${byThickness} by thickness`);
  for (let j = 1; j < C; j++) assert.ok(Math.abs(ocean.h[j] - others[j]) < 1e-6, `cell ${j} without ice keeps its mixed layer: ${others[j]} -> ${ocean.h[j]}`);
});

test('the ocean salts or freshens by exactly the ice volume the ice step froze or melted, and a freeze–melt cycle gives the salt back', () => {
  const { ocean, sea, surfaceT, ice, flux, calm } = frozenColumn();
  let area = 0;
  for (let i = 0; i < C; i++) area += mesh.areaCell[i];
  const frozen = () => (sea.budget.frozen - sea.budget.melted + sea.budget.snowIce / ICE_DENSITY) / area;
  const S0 = ocean.W[0] / ocean.h[0];
  let worst = 0, throughput = 0, saltiest = S0, hours = 0;
  function hour(heating, snowfall = 0) {
    flux.fill(heating);
    const before = frozen(), w = ocean.W[0], s = w / ocean.h[0];
    for (let i = 0; i < C; i++) {
      if (snowfall > 0) sea.deposit(i, snowfall, 260, ice, surfaceT);
      sea.update(surfaceT, ice, flux, i, HOUR);
    }
    ocean.advance(surfaceT, ice, sea.oceanFlux, calm, HOUR, sea.concentration);
    const booked = (s - ICE_SALINITY) * ICE_DENSITY * (frozen() - before) / 1000;
    worst = Math.max(worst, Math.abs(ocean.W[0] - w - booked));
    throughput += Math.abs(booked);
    saltiest = Math.max(saltiest, ocean.W[0] / ocean.h[0]);
    hours++;
  }
  for (let n = 0; n < 120; n++) hour(-250, n >= 48 && n < 60 ? 6 : 0);
  const grown = { thickness: ice[0], cover: sea.concentration[0] };
  while (ice[0] > 0 && hours < 1000) hour(300);
  const S1 = ocean.W[0] / ocean.h[0];
  console.log(`ice ${grown.thickness.toFixed(3)} m over ${(100 * grown.cover).toFixed(0)}% of the cells, then melted, over ${hours} hourly steps: the salt booked matches the ice's frozen volume to ${worst.toExponential(1)} psu·m; the mixed layer went from ${S0} to ${saltiest.toFixed(4)} psu and back to ${S1.toFixed(6)}`);
  assert.ok(sea.budget.leadFrozen > 0 && sea.budget.snowIce > 0 && sea.budget.lateralMelted > 0 && sea.budget.melted > 0, 'the cycle freezes in the leads, floods snow into ice, melts at the edges and melts away');
  assert.ok(grown.cover < 1 && ice[0] === 0, `the ice covered ${grown.cover} of each cell and is gone after ${hours} hours`);
  for (let i = 1; i < C; i++) assert.equal(ocean.W[i], ocean.W[0], `cell ${i} stays as cell 0 does`);
  assert.ok(ocean.u.every((v) => v === 0), 'level columns do not flow');
  assert.ok(worst < 1e-9 * ocean.W[0], `the salt each ocean step books differs from the ice's frozen volume by up to ${worst} psu·m`);
  assert.ok(saltiest - S0 > 0.1 && Math.abs(S1 - S0) < 1e-3 * (saltiest - S0), `salinity ${S0} -> ${saltiest} -> ${S1}`);
  assert.ok(throughput > 10, `${throughput} psu·m of brine booked and returned`);
});

test('the ocean starts and loads under partly covered ice without salting it, taking the cover from the sea ice it is given, as the model gives it', () => {
  const model = createModel(new Grid(2), { ocean: { everySteps: 1, thermoclineTilt: 0, salinityProfile: () => 34 } });
  const { ocean, seaIce } = model, n = model.mesh.nCells;
  const surfaceT = new Float64Array(n).fill(FREEZING_POINT), ice = new Float64Array(n).fill(1), calm = new Float64Array(model.mesh.nEdges);
  seaIce.load(ice, new Float64Array(n).fill(0.4));
  const unchanged = (label) => {
    const before = Float64Array.from(ocean.W.subarray(0, n));
    ocean.advance(surfaceT, ice, seaIce.oceanFlux, calm, HOUR, seaIce.concentration);
    for (let i = 0; i < n; i++) assert.ok(Math.abs(ocean.W[i] - before[i]) < 1e-12 * before[i], `${label}: cell ${i} salt ${before[i]} -> ${ocean.W[i]}`);
  };
  ocean.initialize(surfaceT, ice);
  unchanged('after initialize');
  ocean.load(ocean.serialize(), surfaceT, ice);
  unchanged('after load');
  const alone = createOcean(model.mesh, { everySteps: 1, thermoclineTilt: 0, salinityProfile: () => 34 });
  alone.initialize(surfaceT, ice);
  const before = alone.W[0];
  alone.advance(surfaceT, ice, new Float64Array(n), calm, HOUR, seaIce.concentration);
  assert.ok(alone.W[0] < before - 10, `an ocean that started under full cover books the cover falling to 0.4 as melt: ${before} -> ${alone.W[0]}`);
});
