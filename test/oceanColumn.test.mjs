import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createOcean } from '../js/ocean/layered.module.js';
import { HYBRID_LAYERS, JERLOV, transmitted, MODE_DEPTH, boundaryLayer, interiorMixing, mixedLayerDepth, columnScratch } from '../js/ocean/column.module.js';
import { seawaterDensity } from '../js/ocean/seawater.module.js';
import { seasonalColumn, createSeasonalColumn, slabCycle, columnMesh, nearestCell, atLatitude, RHO_CP, AIR_COLUMN, SLAB, STORE, YEAR_DAYS, DAY } from '../explain/oceanColumn.module.js';
import { columnOptions } from '../explain/seaStore.module.js';

const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const C = (k) => k - 273.15;

test('a cubic meter of seawater holds 4.08 MJ per kelvin, and the whole column of air above a square meter holds as much as the top 2.5 m of the sea', () => {
  const ocean = createOcean(buildMesh(new Grid(16)), { everySteps: 1 });
  assert.equal(RHO_CP, ocean.rhoCp, 'the ocean\'s own ρc_p');
  close(RHO_CP, 4.08e6, 0.005e6, 'ρc_p');
  close(AIR_COLUMN, 1.02e7, 0.005e7, 'c_p·1000 hPa/g');
  close(AIR_COLUMN / RHO_CP, 2.5, 0.01, 'meters of sea per air column');
  close(SLAB / RHO_CP, 10.28, 0.005, 'Part 7\'s slab in meters');
  close(SLAB / AIR_COLUMN, 4.1, 0.05, 'Part 7\'s slab in air columns');
  close(HYBRID_LAYERS.minimum[0] * RHO_CP / AIR_COLUMN, 1.2, 0.01, 'the top layer in air columns');
  close(HYBRID_LAYERS.depths[12], 57, 0.5, 'the twelve surface layers, meters');
  close(HYBRID_LAYERS.minimum[11], 7, 0.05, 'the twelfth layer, meters');
  close(HYBRID_LAYERS.depths[21], 144, 0.5, 'the next nine layers end at'); close(HYBRID_LAYERS.depths[44], 374, 0.5, 'the 23 layers of 10 m end at'); close(HYBRID_LAYERS.depths[58], 3019, 1, 'the minimum thicknesses add to');
  close(250 * DAY / (3 * RHO_CP), 1.76, 0.01, 'kelvins a day of 250 W/m² would warm 3 m of water');
  close(4000 * RHO_CP, 1.63e10, 0.005e10, 'the whole 4000 m column');
  console.log(`ρc_p ${(RHO_CP / 1e6).toFixed(4)} MJ/m³/K; the air column ${(AIR_COLUMN / 1e7).toFixed(4)}e7 J/m²/K = ${(AIR_COLUMN / RHO_CP).toFixed(3)} m of sea; Part 7's slab ${(SLAB / AIR_COLUMN).toFixed(2)} air columns, the 3 m top layer ${(3 * RHO_CP / AIR_COLUMN).toFixed(2)}`);
});

test('every cell of the column is at the latitude asked for, the column read is a cell away from the poles, and the stress it stirs with is the wind\'s', () => {
  const m = columnMesh();
  assert.equal(m.nCells, 42);
  const column = createSeasonalColumn({ days: 0, latitude: 40 });
  close(m.latCell[column.cell] * 180 / Math.PI, 31.72, 0.01, 'the cell read');
  assert.equal(column.cell, nearestCell(m, 32));
  close(column.stirring, 0.05, 0.0015, 'the stress at the column');
  assert.equal(column.layers, HYBRID_LAYERS.K, 'all the hybrid layers');
  assert.equal(column.layers, 58);
  const placed = atLatitude(m, 40);
  assert.ok(placed.latCell.every((v) => v === 40 * Math.PI / 180) && placed.fCell.every((v) => v === 2 * m.omega * Math.sin(40 * Math.PI / 180)));
  const o = column.ocean, n = m.nCells;
  for (let k = 0; k < 12; k++) for (let i = 1; i < n; i++) { assert.equal(o.h[k * n + i], o.h[k * n]); close(o.Q[k * n + i], o.Q[k * n], 1e-9, `cell ${i} layer ${k}`); }
  assert.ok(column.stirring > 0.04);
});

test('the slabs answer a year of heating as the linear response says: 4.6 K and 23 days behind for 10 m at λ = 20, 7.7 K and 40 days at λ = 10, 2.2 K and 65 days for 50 m', () => {
  const per100 = (capacity, damping) => slabCycle({ capacity, heating: 100, damping });
  const cases = [[SLAB, 20, 4.6, 23], [SLAB, 10, 7.7, 40], [50 * RHO_CP, 20, 2.2, 65], [50 * RHO_CP, 10, 2.4, 77], [100 * RHO_CP, 20, 1.2, 77], [100 * RHO_CP, 10, 1.2, 84]];
  for (const [capacity, damping, amplitude, lag] of cases) {
    const cycle = per100(capacity, damping);
    close(cycle.amplitude, amplitude, 0.05, `amplitude for ${(capacity / RHO_CP).toFixed(0)} m at λ ${damping}`);
    close(cycle.lag, lag, 0.6, `lag for ${(capacity / RHO_CP).toFixed(0)} m at λ ${damping}`);
  }
  const slabs = seasonalColumn({ sea: 'slab', days: 2 * YEAR_DAYS });
  const second = slabs.readouts(slabs.samples - 1, 'slab'), cycle = slabCycle({ capacity: SLAB, heating: 150, damping: 20 });
  close(second.range, cycle.range, 0.03 * cycle.range, 'the second year\'s range');
  close(second.lag, cycle.lag, 1, 'the second year\'s lag');
});

test('the sunlight the heating carries passes the 3 m top layer by Jerlov type IB, 31 percent of it, and is deposited over the layers beneath', () => {
  const heating = 100, dt = 10800;
  const run = (light) => createSeasonalColumn({ latitude: 45, heating, damping: 0, wind: 0, mean: 286.15, start: YEAR_DAYS / 4, days: dt / DAY, dt, light }).run();
  const lit = run(true), dark = run(false), o = lit.ocean, n = columnMesh().nCells, cell = lit.cell;
  const share = transmitted(3, 'IB');
  close(share, 0.31, 0.005, 'the share that passes 3 m');
  assert.equal(o.transmission[cell], share);
  close(transmitted(50, 'IB'), 0.0173, 0.0005, 'the share that passes 50 m');
  const start = createSeasonalColumn({ latitude: 45, mean: 286.15, days: 0 }).ocean;
  const top = (r) => RHO_CP * (r.ocean.Q[cell] - start.Q[cell]);
  let below = 0;
  for (let k = 1; k < lit.layers; k++) below += RHO_CP * (o.Q[k * n + cell] - start.Q[k * n + cell]);
  const sun = (heating * Math.sin(2 * Math.PI * (YEAR_DAYS / 4 * DAY + 0.5 * dt) / (YEAR_DAYS * DAY))) * dt;
  console.log(`one step of ${heating} W/m² of sunlight: the top layer gains ${(top(lit) / sun).toFixed(3)} of it, the layers beneath ${(below / sun).toFixed(3)}; without the light option the top layer gains ${(top(dark) / sun).toFixed(3)}`);
  close((top(dark) - top(lit)) / sun, share, 0.005, 'the share the light option moves below the top layer');
  close(top(dark) / sun, 1, 0.03, 'without the option the top layer takes nearly all of it');
  close(top(lit) + below, sun, 1e-6 * sun, 'all of it is in the column');
});

test('with no wind the column holds exactly the heat prescribed, and with wind to within 6 percent of the year\'s net', () => {
  const calm = seasonalColumn({ latitude: 45, heating: 133, damping: 10, wind: 0, mean: 286.15, dt: 10800, every: 8, days: 120, light: true });
  const heat = calm.record.seas.model.heat[calm.samples - 1];
  console.log(`120 days without wind: prescribed ${(calm.heatIn.model / 1e9).toFixed(4)} GJ/m², the column holds ${(heat / 1e9).toFixed(4)}`);
  close(heat, calm.heatIn.model, 1e-6 * Math.abs(calm.heatIn.model), 'heat in the column, without wind');
  const windy = seasonalColumn({ latitude: 45, heating: 133, damping: 10, wind: 0.1, mean: 286.15, dt: 10800, every: 8, days: 120, light: true });
  const held = windy.record.seas.model.heat[windy.samples - 1];
  console.log(`120 days at 0.1 N/m²: prescribed ${(windy.heatIn.model / 1e9).toFixed(4)} GJ/m², the column holds ${(held / 1e9).toFixed(4)}`);
  close(held, windy.heatIn.model, 0.06 * Math.abs(windy.heatIn.model), 'heat in the column, with wind');
});

test('at 15°N a year under 0.1 N/m² swings the surface by 2.2 K where Part 7\'s slab swings 7.7 K, and with no wind 4.6 K; the layers of the thermocline are isopycnal from 264 m', () => {
  const windy = seasonalColumn(columnOptions({ latitude: 15 })), calm = seasonalColumn(columnOptions({ latitude: 15, wind: 0 }));
  const swing = (run, name) => { const s = run.record.seas[name].sst; return Math.max(...s) - Math.min(...s); };
  console.log(`15°N: surface swing ${swing(windy, 'model').toFixed(2)} K with wind, ${swing(calm, 'model').toFixed(2)} K without; slab ${swing(windy, 'slab').toFixed(2)} K; mixed layer ${Math.min(...windy.record.mixedDepth.slice(20)).toFixed(1)}–${Math.max(...windy.record.mixedDepth).toFixed(1)} m with wind`);
  close(swing(windy, 'model'), 2.25, 0.05, 'surface swing with wind');
  close(swing(calm, 'model'), 4.64, 0.05, 'surface swing without wind');
  close(swing(windy, 'slab'), 7.66, 0.05, 'the slab\'s swing');
  const layers = windy.layersAt(120), iso = layers.filter((l) => l.mode > 0.01 && l.top > 100);
  close(Math.min(...windy.record.mixedDepth.slice(20)), 53.9, 1, 'the shallowest mixed layer'); close(Math.max(...windy.record.mixedDepth), 119.1, 2, 'the deepest');
  assert.equal(iso[0].layer, 33);
  close(iso[0].top, 264, 2, 'top of the first isopycnal layer');
  assert.ok(iso.every((l) => l.mode > 0.99), 'isopycnal layers sit well clear of their minimum');
  assert.ok(layers.slice(0, 33).every((l) => l.mode < 0.01), 'the 33 layers above are held at their minimum thicknesses');
  assert.ok(MODE_DEPTH === 5);
});

test('the column physics the prose gives: the boundary layer is where the work of mixing equals ρ₀ m* u*³ Δt + n* E_c, the interior mixing is Pacanowski and Philander\'s, and the mixed layer is where the density first passes the top layer\'s by 0.03 kg/m³', () => {
  const K = 30, h = new Float64Array(K).fill(5), S = new Float64Array(K).fill(35), T = Float64Array.from({ length: K }, (_, k) => 290 - 0.01 * 5 * k), kappa = new Float64Array(K), scratch = columnScratch(K);
  const ustar = 0.0098, f = 1e-4, dt = 10800, rho0 = 1026;
  const calm = boundaryLayer({ h, T, S, K, ustar, f, dt, langmuir: 0 }, kappa, scratch);
  const rotation = 1 + (calm.depth * f / (0.5 * ustar)) ** 2;
  close(0.5 * Math.sqrt(0.1 / 1025) / (2 * 7.2921159e-5 * Math.sin(Math.PI / 4)), 47.9, 0.1, 'the rotation limit at 45°N under 0.1 N/m²');
  console.log(`u* ${ustar} m/s, f ${f}: boundary layer ${calm.depth.toFixed(2)} m, energy ${calm.energy.toFixed(3)} J/m², work ${calm.work.toFixed(3)} J/m²`);
  assert.ok(calm.depth > 5 && calm.depth < 150, `depth ${calm.depth}`);
  close(calm.energy, rho0 * dt * 1.25 * ustar ** 3 / rotation, 1e-9 * calm.energy, 'energy from the wind');
  close(calm.work, calm.energy, 1e-4 * calm.energy, 'the work of mixing equals the energy');
  const released = 3, convecting = boundaryLayer({ h, T, S, K, ustar, f, dt, langmuir: 0, released }, kappa, scratch);
  close(convecting.energy, rho0 * dt * 1.25 * ustar ** 3 / (1 + (convecting.depth * f / (0.5 * ustar)) ** 2) + 0.2 * released, 1e-9 * convecting.energy, 'energy with convection');
  close(convecting.velocity, Math.cbrt(ustar ** 3 + released / (rho0 * dt)), 1e-12, 'v*');
  assert.ok(convecting.depth > calm.depth, 'convection deepens the layer');
  let z = 0;
  boundaryLayer({ h, T, S, K, ustar, f, dt, langmuir: 0 }, kappa, scratch);
  for (let k = 0; k < 3; k++) { z += h[k]; close(kappa[k], 0.4 * ustar * z * (1 - z / calm.depth) ** 2, 1e-9, `diffusivity at ${z} m`); }
  const warm = Float64Array.from({ length: K }, (_, k) => (k < 2 ? 291 : 290)), du = [0.3, 0.05], n = new Float64Array(K), kp = new Float64Array(K), shear = new Float64Array(K);
  shear[0] = du[0] ** 2; shear[1] = du[1] ** 2;
  interiorMixing({ h, T: warm, S, shear, K }, kp, n, scratch);
  const g = 9.81, b = (k) => g * (seawaterDensity(warm[k + 1], 35) - seawaterDensity(warm[k], 35)) / rho0;
  for (const k of [0, 1]) {
    const Ri = b(k) * 5 / shear[k], nu = 1e-2 / (1 + 5 * Ri) ** 2 + 1e-4;
    console.log(`Ri ${Ri.toFixed(3)} at interface ${k}`);
    close(n[k], nu, 1e-9, `viscosity ${k}`); close(kp[k], nu / (1 + 5 * Ri) + 1e-5, 1e-9, `diffusivity ${k}`);
  }
  const unstable = Float64Array.from({ length: K }, (_, k) => 290 + (k === 3 ? 1 : 0)), nu2 = new Float64Array(K), ka2 = new Float64Array(K);
  interiorMixing({ h, T: unstable, S, K }, ka2, nu2, scratch);
  assert.equal(ka2[2], 0.1); assert.equal(ka2[3] <= 0.1, true);
  const rho = Float64Array.from({ length: 6 }, (_, k) => 1025 + (k < 3 ? 0 : 0.06 + 0.01 * (k - 3))), hh = new Float64Array(6).fill(10);
  close(mixedLayerDepth(hh, rho, 6), 30, 1e-9, 'the mixed layer, midway between the centers of the layers that straddle the excess');
  assert.equal(JERLOV.IB.join(' '), '0.67 1 17');
});
