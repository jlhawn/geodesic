import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createGeography, readTopography } from '../js/geography.module.js';
import { createLandSurface, TRACE_SNOW } from '../js/physics/land.module.js';
import { monthWeights } from '../js/monthCalendar.module.js';

const DAY = 86400, DT = 900;
const mesh = buildMesh(new Grid(16)), C = mesh.nCells;
const geo = createGeography(mesh, readTopography('etopo1'), { soilWaterMonthly: true });
const open = (i) => geo.land[i] && !geo.iceSheet[i];
const constant = (value) => Array.from({ length: 12 }, () => new Float64Array(C).fill(value));
const withMonths = (months) => ({ ...geo, wetnessMonthly: months });

function build(geography, options = {}) {
  const land = createLandSurface(mesh, geography, options);
  land.initialize(0);
  return land;
}
function step(land, dt = DT, count = 1, time = 0, kelvin = 285) {
  const T = new Float64Array(C).fill(kelvin), flux = new Float64Array(C);
  for (let n = 0; n < count; n++) { land.setTime(time + n * dt); for (let i = 0; i < C; i++) if (geo.land[i]) land.update(i, T, flux, 0, dt); }
}
const pick = (test) => { for (let i = 0; i < C; i++) if (test(i)) return i; throw new Error('no such cell'); };

test('the soil goal is refused for other values, a bad time, a land without vegetation and a geography without the monthly field', () => {
  assert.throws(() => createLandSurface(mesh, geo, { soilGoal: 'wet' }), /soilGoal is 'free' or 'observed'/);
  assert.throws(() => createLandSurface(mesh, geo, { soilGoal: 'observed', soilGoalTime: 0 }), /soilGoalTime/);
  assert.throws(() => createLandSurface(mesh, geo, { soilGoal: 'observed', soilGoalTime: NaN }), /soilGoalTime/);
  assert.throws(() => createLandSurface(mesh, geo, { soilGoal: 'observed', vegetation: false }), /needs vegetation/);
  const { wetnessMonthly, ...without } = geo;
  assert.throws(() => createLandSurface(mesh, without, { soilGoal: 'observed' }), /monthly soil water/);
  assert.equal(createLandSurface(mesh, geo, { soilGoal: 'free' }).goaled, false);
});

test("with soilGoal 'free' the land has no goal arrays and no goal in its budget or save, and steps as without the option", () => {
  const a = build(geo), b = build(geo, { soilGoal: 'free' });
  for (const land of [a, b]) { for (let i = 0; i < C; i++) if (open(i)) land.soil[i] = 90; step(land, DT, 6); }
  assert.ok(!('goal' in b.budget) && b.goalAdded === null);
  assert.deepEqual(Object.keys(b.serialize()).filter((k) => k.startsWith('goal')), []);
  for (const name of ['soil', 'surface', 'vegetation', 'canopy']) assert.deepEqual(Array.from(a[name]), Array.from(b[name]), name);
  assert.doesNotThrow(() => b.setTime(1e7));
  assert.doesNotThrow(() => b.resetIncrements());
});

test('the root zone and the surface layer relax to fill x capacity with the e-folding time, and the water is booked as the goal', () => {
  const tau = 10 * DAY, land = build(withMonths(constant(0.3)), { soilGoal: 'observed', soilGoalTime: tau });
  const i = pick(open), k = pick((n) => open(n) && n !== i);
  land.soil[i] = 240; land.surface[i] = 0; land.soil[k] = 10; land.surface[k] = 15;
  const cap = land.capacity(i), steps = Math.round(tau / DT);
  step(land, DT, steps);
  const left = Math.exp(-1);
  assert.ok(Math.abs(land.soil[i] - (0.3 * cap + (240 - 0.3 * cap) * left)) < 5, `soil ${land.soil[i]}`);
  assert.ok(land.goalRemoved[i] > 50 && land.goalAdded[i] > 0, 'water out of the root zone, into the surface layer');
  assert.ok(land.soil[k] > 10 + 0.3 * cap * 0.5 && land.surface[k] < 15, 'the dry cell fills and the full surface layer drains');
  let booked = 0;
  for (let n = 0; n < C; n++) if (geo.land[n]) booked += mesh.areaCell[n] * (land.goalAdded[n] - land.goalRemoved[n]);
  assert.ok(Math.abs(land.budget.goal - booked) < 1e-9 * Math.abs(booked) + 1, `budget ${land.budget.goal} vs ${booked}`);
});

test('goalAdded less goalRemoved is the change of the root zone and the surface layer when nothing else moves water', () => {
  const land = build(withMonths(constant(0.45)), { soilGoal: 'observed' });
  for (let n = 0; n < C; n++) if (geo.land[n]) { land.soil[n] = (n % 7) * 40; land.surface[n] = (n % 3) * 7; }
  const before = Float64Array.from({ length: C }, (_, n) => land.soil[n] + land.surface[n]);
  step(land, DT, 8);
  let moved = 0;
  for (let n = 0; n < C; n++) {
    if (!open(n)) continue;
    const change = land.soil[n] + land.surface[n] - before[n];
    assert.ok(Math.abs(change - (land.goalAdded[n] - land.goalRemoved[n])) < 1e-9, `cell ${n}`);
    assert.ok(land.goalAdded[n] >= 0 && land.goalRemoved[n] >= 0);
    if (land.goalAdded[n] + land.goalRemoved[n] > 0) moved++;
  }
  assert.ok(moved > 0.5 * geo.land.reduce((s, v, n) => s + (open(n) ? 1 : 0), 0));
});

test('the goal leaves ice sheets, snow-covered cells, the sea and cells without a goal alone', () => {
  const months = constant(0.3), hole = pick((n) => open(n) && n > 5), snowy = pick((n) => open(n) && n > hole);
  for (const m of months) m[hole] = NaN;
  const land = build(withMonths(months), { soilGoal: 'observed' });
  for (let n = 0; n < C; n++) { land.soil[n] = 120; land.surface[n] = 7; }
  land.snow[snowy] = 10 * TRACE_SNOW * 1e6;
  const sheet = geo.land.findIndex((v, n) => v && geo.iceSheet[n]);
  assert.ok(sheet >= 0, 'N = 16 has an ice-sheet cell');
  const before = { soil: Float64Array.from(land.soil), surface: Float64Array.from(land.surface) };
  step(land, DT, 4, 0, 260);
  for (const n of [hole, snowy, sheet]) assert.equal(land.goalAdded[n] + land.goalRemoved[n], 0, `cell ${n}`);
  for (const n of [hole, snowy]) assert.equal(land.soil[n] + land.surface[n], before.soil[n] + before.surface[n], `cell ${n}`);
  for (let n = 0; n < C; n++) if (!geo.land[n]) { assert.equal(land.goalAdded[n] + land.goalRemoved[n], 0); assert.equal(land.soil[n], before.soil[n]); }
  assert.ok(land.goalRemoved[pick((n) => open(n) && n !== hole && n !== snowy)] > 0);
});

test('the goal fill is the mid-month values interpolated on the model calendar, cyclic across December and January', () => {
  const months = Array.from({ length: 12 }, (_, m) => new Float64Array(C).fill(0.05 * (m + 1)));
  const i = pick(open);
  for (const day of [0, 15, 45, 79, 300, 330, 400, -20]) {
    const land = build(withMonths(months), { soilGoal: 'observed', soilGoalTime: DT });
    land.soil[i] = 0; land.surface[i] = 0;
    land.setTime(day * DAY);
    const T = new Float64Array(C).fill(285), flux = new Float64Array(C);
    land.update(i, T, flux, 0, DT);
    const [m0, m1, f] = monthWeights(day), fill = months[m0][i] * (1 - f) + months[m1][i] * f, a = 1 - Math.exp(-1);
    assert.ok(Math.abs(land.goalAdded[i] - a * fill * (land.capacity(i) + 15)) < 1e-9 * 300, `day ${day}: ${land.goalAdded[i]}`);
  }
});

test('the cover and tree increments are the observed goal\'s change less the model\'s own, and zero under snow or with the cover goal at water', () => {
  const months = constant(0.35), geography = withMonths(months);
  const A = build(geography, { soilGoal: 'observed', coverGoal: 'observed' }), B = build(geography, { soilGoal: 'observed', coverGoal: 'water' });
  const cells = [];
  for (let n = 0; n < C; n++) if (open(n) && geo.greenMax[n] > 0.2) cells.push(n);
  assert.ok(cells.length > 20);
  for (const land of [A, B]) for (const n of cells) { land.vegetation[n] = 0.2 + 0.5 * ((n % 5) / 5); land.canopy[n] = 0.1; land.soil[n] = 100 + (n % 4) * 40; land.snow[n] = 0; }
  const snowy = cells[3];
  for (const land of [A, B]) land.snow[snowy] = 20;
  const before = { vegetation: Float64Array.from(A.vegetation), canopy: Float64Array.from(A.canopy) };
  const tau = 3600;
  const T = new Float64Array(C).fill(288), flux = new Float64Array(C);
  for (const land of [A, B]) { land.setTime(0); for (const n of cells) land.update(n, T.slice(), flux, 0, tau); }
  let some = 0;
  for (const n of cells) {
    const g = (land, name) => land[name][n];
    if (n === snowy) { assert.equal(A.goalCover[n], 0); assert.equal(A.goalTrees[n], 0); continue; }
    const own = g(B, 'vegetation') - before.vegetation[n];
    assert.ok(Math.abs(A.goalCover[n] - (g(A, 'vegetation') - before.vegetation[n] - own)) < 1e-12, `cover ${n}`);
    const f = A.treeFactor(n) * g(A, 'vegetation'), was = before.canopy[n], rate = 1 - Math.exp(-tau / (f > was ? 10 * 365 * DAY : 3 * 365 * DAY));
    assert.ok(Math.abs(A.goalTrees[n] - (g(A, 'canopy') - was - (f - was) * rate)) < 1e-12, `trees ${n}`);
    assert.equal(B.goalCover[n], 0); assert.equal(B.goalTrees[n], 0);
    if (Math.abs(A.goalCover[n]) > 1e-9) some++;
  }
  assert.ok(some > 0, 'the two goals differ somewhere');
});

test('the four arrays reach a save and load back, a save without them loads zeros, and resetIncrements clears them', () => {
  const land = build(withMonths(constant(0.2)), { soilGoal: 'observed', coverGoal: 'observed' });
  for (let n = 0; n < C; n++) if (open(n)) land.soil[n] = 200;
  step(land, DT, 3);
  const saved = land.serialize();
  for (const name of ['goalAdded', 'goalRemoved', 'goalCover', 'goalTrees']) assert.ok(saved[name] instanceof Float64Array && saved[name].length === C, name);
  assert.ok(saved.goalRemoved.some((v) => v > 0));
  const again = build(withMonths(constant(0.2)), { soilGoal: 'observed', coverGoal: 'observed' });
  again.load(saved);
  for (const name of ['goalAdded', 'goalRemoved', 'goalCover', 'goalTrees']) assert.deepEqual(Array.from(again[name]), Array.from(land[name]), name);
  const { goalAdded, goalRemoved, goalCover, goalTrees, ...plain } = saved;
  again.load(plain);
  for (const name of ['goalAdded', 'goalRemoved', 'goalCover', 'goalTrees']) assert.ok(again[name].every((v) => v === 0), name);
  land.resetIncrements();
  for (const name of ['goalAdded', 'goalRemoved', 'goalCover', 'goalTrees']) assert.ok(land[name].every((v) => v === 0), name);
});

test('regridLand keeps no goal increments, at the same cells or on another mesh', async () => {
  const { regridLand } = await import('../js/physics/regrid.module.js');
  const land = build(geo, { soilGoal: 'observed', soilGoalTime: 3600 });
  for (let i = 0; i < C; i++) if (open(i)) land.soil[i] = 250;
  step(land, DT, 4);
  const saved = land.serialize();
  assert.ok(saved.goalRemoved.some((x) => x > 0));
  const target = { mesh: buildMesh(new Grid(32)) };
  target.geography = createGeography(target.mesh, readTopography('etopo1'));
  const same = regridLand({ mesh, geography: geo }, { mesh, geography: geo, land }, saved);
  assert.deepEqual(Object.keys(same).filter((name) => name.startsWith('goal')), []);
  const finer = regridLand({ mesh, geography: geo }, { mesh: target.mesh, geography: target.geography, land: createLandSurface(target.mesh, target.geography, {}) }, saved);
  assert.deepEqual(Object.keys(finer).filter((name) => name.startsWith('goal')), []);
});
