import { test } from 'node:test';
import { HYBRID_LAYERS } from '../js/ocean/column.module.js';
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
  assert.equal(column.layers, HYBRID_LAYERS.K, 'all the hybrid layers');
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


