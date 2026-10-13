import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createOcean } from '../js/ocean/layered.module.js';
import { seasonalColumn, createSeasonalColumn, slabCycle, columnMesh, nearestCell, RHO_CP, AIR_COLUMN, SLAB, SEAS, FLOOR, CAP, YEAR_DAYS, DAY } from '../explain/oceanColumn.module.js';

const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const C = (k) => k - 273.15;
const PER_DAY = 16;

test('a cubic meter of seawater holds 4.08 MJ per kelvin, and the whole column of air above a square meter holds as much as the top 2.5 m of the sea', () => {
  const ocean = createOcean(buildMesh(new Grid(2)), { everySteps: 1 });
  assert.equal(RHO_CP, ocean.rhoCp, 'the ocean\'s own ρc_p');
  close(RHO_CP, 4.08e6, 0.005e6, 'ρc_p');
  close(AIR_COLUMN, 1.02e7, 0.005e7, 'c_p·1000 hPa/g');
  close(AIR_COLUMN / RHO_CP, 2.5, 0.01, 'meters of sea per air column');
  close(SLAB / RHO_CP, 10.28, 0.005, 'Part 7\'s slab in meters');
  close(SLAB / AIR_COLUMN, 4.1, 0.05, 'Part 7\'s slab in air columns');
  close(SEAS.fixed / AIR_COLUMN, 20, 0.05, 'the 50 m floor in air columns');
  close(CAP * RHO_CP, 2.45e9, 0.005e9, 'the 600 m cap');
  close(4000 * RHO_CP, 1.63e10, 0.005e10, 'the whole 4000 m column');
  console.log(`ρc_p ${(RHO_CP / 1e6).toFixed(4)} MJ/m³/K; the air column ${(AIR_COLUMN / 1e7).toFixed(4)}e7 J/m²/K = ${(AIR_COLUMN / RHO_CP).toFixed(3)} m of sea; Part 7's slab ${(SLAB / AIR_COLUMN).toFixed(2)} air columns, the 50 m floor ${(SEAS.fixed / AIR_COLUMN).toFixed(2)}`);
});

test('the column is the cell of Grid(2) nearest the latitude asked for, and the stress it stirs with is the wind\'s', () => {
  const m = columnMesh();
  assert.equal(m.nCells, 42);
  const column = createSeasonalColumn({ days: 0 });
  close(column.latitude, 31.72, 0.01, 'the default column');
  close(column.stirring, 0.05, 0.0015, 'the stress at the column');
  assert.equal(column.layers, 47, 'the mixed layer and all 46 classes');
  close(m.latCell[nearestCell(m, 60)] * 180 / Math.PI, 58.28, 0.01, 'the cell nearest 60°');
});

test('the slabs answer a year of heating as the linear response says: 4.6 K and 23 days behind for 10 m at λ = 20, 7.7 K and 40 days at λ = 10, 2.2 K and 65 days for 50 m', () => {
  const per100 = (capacity, damping) => slabCycle({ capacity, heating: 100, damping });
  const cases = [[SLAB, 20, 4.6, 23], [SLAB, 10, 7.7, 40], [SEAS.fixed, 20, 2.2, 65], [SEAS.fixed, 10, 2.4, 77], [100 * RHO_CP, 20, 1.2, 77], [100 * RHO_CP, 10, 1.2, 84]];
  for (const [capacity, damping, amplitude, lag] of cases) {
    const cycle = per100(capacity, damping);
    close(cycle.amplitude, amplitude, 0.05, `amplitude for ${(capacity / RHO_CP).toFixed(0)} m at λ ${damping}`);
    close(cycle.lag, lag, 0.6, `lag for ${(capacity / RHO_CP).toFixed(0)} m at λ ${damping}`);
  }
  const slabs = seasonalColumn({ sea: 'slab', days: 2 * YEAR_DAYS });
  for (const name of ['slab', 'fixed']) {
    const second = slabs.readouts(slabs.samples - 1, name), cycle = slabCycle({ capacity: SEAS[name], heating: 150, damping: 20 });
    close(second.range, cycle.range, 0.03 * cycle.range, `${name}: the second year's range`);
    close(second.lag, cycle.lag, 1, `${name}: the second year's lag`);
  }
});

test('heated without damping, the model\'s mixed layer sits on its 50 m floor from day 30 of the half-year of heating while its surface warms from 15.0 to 21.4 °C, and the column stores the 1.51 GJ/m² it is given', () => {
  const t0 = performance.now();
  const column = seasonalColumn({ damping: 0, days: YEAR_DAYS / 2 });
  const seconds = (performance.now() - t0) / 1000, { record: r } = column, N = column.samples - 1;
  let low = Infinity, high = 0;
  for (let n = 30 * PER_DAY; n <= N; n++) { low = Math.min(low, r.mixedDepth[n]); high = Math.max(high, r.mixedDepth[n]); }
  const given = 150 * YEAR_DAYS * DAY / Math.PI;
  console.log(`half a year in ${seconds.toFixed(1)} s (${(1000 * seconds / column.steps).toFixed(2)} ms a step): mixed layer ${low.toFixed(2)}–${high.toFixed(2)} m from day 10; SST ${C(r.seas.model.sst[0]).toFixed(2)} → ${C(r.seas.model.sst[N]).toFixed(2)} °C, a fixed 50 m layer's ${C(r.seas.fixed.sst[N]).toFixed(2)}; the column gained ${(r.seas.model.heat[N] / 1e9).toFixed(4)} GJ/m² of the ${(given / 1e9).toFixed(4)} given, ${(r.stored[N] / 1e9).toFixed(4)} of it above 600 m`);
  close(r.time[N], 182.5, 1e-9, 'half a year');
  assert.ok(low >= FLOOR && high < FLOOR + 0.6, `the mixed layer stays within 0.6 m of its floor: ${low}–${high} m`);
  close(C(r.seas.model.sst[0]), 15, 1e-4, 'the start');
  close(C(r.seas.model.sst[N]), 21.38, 0.01, 'the surface after the heating');
  close(C(r.seas.fixed.sst[N]), 15 + given / SEAS.fixed, 0.01, 'a fixed 50 m layer warms by 7.4 K');
  close(given / SEAS.fixed, 7.37, 0.01, 'the warming of the fixed layer');
  close(column.heatIn.model, given, 1e-4 * given, 'the heat put in');
  close(r.seas.model.heat[N], given, 0.005 * given, 'the column keeps the heat put into it, within what its neighbors exchange with it');
  const seasonal = column.layersAt(N).filter((layer) => layer.layer > 0 && layer.bottom < 70);
  assert.ok(seasonal.length >= 5 && seasonal.every((layer) => layer.density < 1025.25 && layer.temperature < r.seas.model.sst[N]), `the summer leaves thin, light classes under the mixed layer: ${seasonal.map((layer) => layer.density).join(', ')}`);
});

test('a year with damping: the floor from day 30 while the column gains heat, convection to 356 m by the end of winter, and a spring retreat that halves in 0.75 days, the model\'s one-day detrainment, and stands 15 m above the floor at the March equinox', () => {
  const t0 = performance.now();
  const column = seasonalColumn();
  const seconds = (performance.now() - t0) / 1000, { record: r } = column, N = column.samples - 1;
  let warming = 0, high = 0;
  for (let n = 30 * PER_DAY; n <= N && r.time[n] < 200; n++) if (r.seas.model.flux[n] > 0) { warming = r.time[n]; high = Math.max(high, r.mixedDepth[n]); }
  let deepest = 0;
  for (let n = 0; n <= N; n++) if (r.mixedDepth[n] > r.mixedDepth[deepest]) deepest = n;
  let half = deepest;
  while (r.mixedDepth[half] - FLOOR > (r.mixedDepth[deepest] - FLOOR) / 2) half++;
  const at = (day) => r.mixedDepth[day * PER_DAY];
  const events = column.events().filter((e) => e.meters > 10);
  const model = column.readouts(N, 'model'), slab = column.readouts(N, 'slab'), fixed = column.readouts(N, 'fixed');
  console.log(`a year in ${seconds.toFixed(1)} s (${(1000 * seconds / column.steps).toFixed(2)} ms a step); the column gains heat until day ${warming.toFixed(1)}, the mixed layer within ${(high - FLOOR).toFixed(2)} m of its floor from day 30; ${at(240).toFixed(1)}, ${at(270).toFixed(1)}, ${at(300).toFixed(1)} and ${at(330).toFixed(1)} m on days 240, 270, 300 and 330; deepest ${r.mixedDepth[deepest].toFixed(1)} m on day ${r.time[deepest].toFixed(2)}; half-way back in ${(r.time[half] - r.time[deepest]).toFixed(2)} days, ${at(YEAR_DAYS).toFixed(1)} m at the March equinox; ${events.map((e) => `${e.kind} ${e.meters.toFixed(0)} m, days ${e.from.toFixed(1)}–${e.to.toFixed(1)}, class ${e.classes.join('/')}`).join('; ')}`);
  console.log(`SST range and lag behind the strongest heating: 10 m slab ${slab.range.toFixed(2)} K, ${slab.lag.toFixed(1)} d; 50 m ${fixed.range.toFixed(2)} K, ${fixed.lag.toFixed(1)} d; the model's mixed layer ${model.range.toFixed(2)} K, ${model.lag.toFixed(1)} d; the column ends ${(r.seas.model.heat[N] / 1e9).toFixed(3)} GJ/m² against ${(column.heatIn.model / 1e9).toFixed(3)} put in`);
  close(warming, 153.0, 0.1, 'the last day the column gains heat');
  assert.ok(high < FLOOR + 1, `the floor while the column gains heat: ${high} m`);
  close(at(240), 66.9, 0.1, 'day 240'); close(at(270), 128.7, 0.1, 'day 270'); close(at(300), 279.6, 0.1, 'day 300'); close(at(330), 341.7, 0.1, 'day 330');
  close(r.mixedDepth[deepest], 356.5, 0.1, 'the deepest mixed layer');
  close(r.time[deepest], 356.19, 0.01, 'on day');
  close(r.time[half] - r.time[deepest], 0.75, 0.01, 'half the excess detrained');
  close(at(YEAR_DAYS), 64.8, 0.05, 'above the floor at the March equinox');
  assert.deepEqual(events.map((e) => [e.kind, e.classes]), [['entrains', [1025.25]], ['detrains', [1025.25]]], 'winter swallows the class beneath and spring gives its water back');
  close(events[0].meters, 265.9, 0.1, 'entrained from day 265 to 330'); close(events[1].meters, 290.6, 0.1, 'detrained in spring');
  const winter = column.layersAt(deepest), spring = column.layersAt(N);
  close(winter[0].bottom, r.mixedDepth[deepest], 1e-3, 'the profile\'s mixed layer');
  close(C(spring[0].temperature), 13.69, 0.01, 'the spring surface'); close(C(spring[1].temperature), 13.71, 0.01, 'the winter water left in the class beneath');
  close(C(column.temperatureAt(0, 100)), 13.97, 0.01, 'the start at 100 m'); close(C(column.temperatureAt(deepest, 100)), 13.64, 0.01, 'the winter at 100 m');
  close(slab.range, 13.86, 0.01, 'the slab\'s range'); close(slab.lag, 22.4, 0.1, 'the slab\'s lag');
  close(fixed.range, 7.27, 0.01, 'the fixed layer\'s range'); close(fixed.lag, 57.3, 0.1, 'the fixed layer\'s lag');
  close(model.range, 5.03, 0.01, 'the mixed layer\'s range'); close(model.lag, 55.5, 0.1, 'the mixed layer\'s lag');
  close(r.seas.model.heat[N], column.heatIn.model, 0.02 * 1.5e9, 'the column keeps what it is given, within 2 percent of a season\'s heat');
  close(r.stored[N], r.seas.model.heat[N], 0.01e9, 'the year\'s changes lie above 600 m');
});
