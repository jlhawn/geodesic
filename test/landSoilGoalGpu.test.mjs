import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createModel } from '../js/model.module.js';
import { readTopography } from '../js/geography.module.js';
import { MID_MONTH_DAYS } from '../js/monthCalendar.module.js';
import { PREVIOUS_LEVELS, PREVIOUS_PCAPE_MINIMUM, PREVIOUS_BOUNDARY_LAYER, PREVIOUS_AIR_SEA } from './helpers/previousDefaults.mjs';
import { physicsSurfaceKernel } from '../js/gpu/physics.gpu.js';
import { prepareLand } from './helpers/soilGoalLand.mjs';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};
const { layoutFor } = gpuAvailable ? await import('../js/gpu/core.gpu.js') : {};

const topography = readTopography('etopo1');
const FAST = { growthTime: 3 * 3600, declineTime: 2 * 3600, snowDeclineTime: 4 * 3600, treeGrowthTime: 2 * 3600, treeDeclineTime: 3600, treeMoisture: false };
const GOAL = { ...FAST, coverGoal: 'observed', soilGoal: 'observed', soilGoalTime: 12 * 3600 };
const options = (land) => ({ topography, geography: { soilWaterMonthly: true }, land, levels: PREVIOUS_LEVELS, surface: PREVIOUS_AIR_SEA, moist: PREVIOUS_PCAPE_MINIMUM, boundaryLayer: { ...PREVIOUS_BOUNDARY_LAYER, velocityScale: 'capped' } });
const START = MID_MONTH_DAYS[3] * 86400 - 8 * 900;
const FIELDS = ['soil', 'surface', 'deepWater', 'goalAdded', 'goalRemoved', 'goalCover', 'goalTrees'];
const SAVED = { deepWater: 'deep' };

test('the soil goal adds twelve monthly fills and four accumulators of C floats before the 2 m temperature and moves no other slot', { skip: !gpuAvailable && 'webgpu not installed' }, () => {
  const mesh = buildMesh(new Grid(16)), C = mesh.nCells, args = [8, 0, 0, 0, null, false, 0, 0, false, false, false, false, false, 3, false, true, true, false, false, false, false, false];
  const without = layoutFor(mesh, ...args), plain = layoutFor(mesh, ...args, false), goal = layoutFor(mesh, ...args, true);
  assert.equal(without.PH.GOALFILL, undefined);
  assert.deepEqual(plain.PH, without.PH);
  assert.equal(goal.PH.GOALFILL, without.PH.T2M);
  assert.equal(goal.PH.GOALADD, goal.PH.GOALFILL + 12 * C);
  assert.deepEqual([goal.PH.GOALREM, goal.PH.GOALVEG, goal.PH.GOALTREE], [1, 2, 3].map((k) => goal.PH.GOALADD + k * C));
  assert.equal(goal.PH.T2M, goal.PH.total - C);
  assert.equal(goal.PH.total, without.PH.total + 16 * C);
  for (const name of Object.keys(without.PH)) if (name !== 'total') assert.equal(goal.PH[name], without.PH[name] + (name === 'T2M' ? 16 * C : 0), name);
});

test("with soilGoal 'free' every surface kernel is the kernel of a model without the option and carries no goal code", () => {
  for (const held of [false, true]) for (const skin of [false, true]) for (const cover of [false, true]) for (const deep of [false, true]) {
    const phys = skin ? {} : { treeline: true, vegetation: true };
    const plain = physicsSurfaceKernel(held, skin, null, cover, deep, phys);
    const free = physicsSurfaceKernel(held, skin, null, cover, deep, { ...phys, soilGoal: 'free', soilGoalTime: 30 * 86400 });
    const goal = physicsSurfaceKernel(held, skin, null, cover, deep, { ...phys, soilGoal: 'observed', soilGoalTime: 30 * 86400 });
    const which = JSON.stringify({ held, skin, cover, deep });
    assert.equal(free, plain, which);
    assert.doesNotMatch(free, /PH_GOAL/, which);
    assert.match(goal, /PH_GOALFILL/, which);
  }
});

test("the GPU refuses soilGoal 'observed' without vegetation and values other than 'free' and 'observed'", { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  await assert.rejects(createGpuModel(new Grid(16), options({ ...GOAL, vegetation: false })), /vegetation/);
  await assert.rejects(createGpuModel(new Grid(16), options({ ...GOAL, soilGoal: 'wet' })), /soilGoal/);
});

test('from a dry start across a month boundary both engines relax the soil and the surface layer to the observed fill and book the increments alike, and resetIncrements zeroes both', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const cpu = prepareLand(createModel(new Grid(16), options(GOAL))), gpu = prepareLand(await createGpuModel(new Grid(16), options(GOAL)));
  cpu.time = START; gpu.time = START;
  const C = cpu.mesh.nCells, land = cpu.geography.land, open = (i) => land[i] && !cpu.geography.iceSheet[i], cells = Array.from({ length: C }, (_, i) => i).filter(open);
  assert.ok(gpu.gpu.layout.PH.GOALFILL !== undefined, 'the GPU allocates the fills');
  for (let n = 0; n < 16; n++) { cpu.step(900); await gpu.step(900); }
  await gpu.sync();
  const saved = await gpu.land.serialize(), cpuSaved = cpu.land.serialize();
  const report = [];
  for (const name of FIELDS) {
    const a = cpuSaved[SAVED[name] ?? name], b = saved[SAVED[name] ?? name];
    let scale = 0, worst = 0, sum = 0, at = -1;
    for (const i of cells) { scale = Math.max(scale, Math.abs(a[i])); const gap = Math.abs(a[i] - b[i]); if (gap > worst) { worst = gap; at = i; } sum += b[i] - a[i]; }
    report.push({ name, scale, worst, mean: Math.abs(sum / cells.length), at });
  }
  console.log(report.map((r) => `${r.name} scale ${r.scale.toPrecision(3)} worst ${r.worst.toExponential(1)} at ${r.at}, mean gap ${r.mean.toExponential(1)}`).join('\n'));
  for (const r of report) {
    assert.ok(r.mean < 1e-5 * Math.max(r.scale, 1), `${r.name} means apart by ${r.mean} of ${r.scale}`);
    assert.ok(r.worst < 2e-4 * Math.max(r.scale, 1), `${r.name} cell ${r.at} apart by ${r.worst} of ${r.scale}`);
  }
  const moved = (name) => cells.filter((i) => Math.abs(saved[name][i]) > 1e-3).length;
  assert.ok(moved('goalAdded') > 0.5 * cells.length && moved('goalRemoved') > 0.1 * cells.length, `${moved('goalAdded')} cells gained and ${moved('goalRemoved')} lost water of ${cells.length}`);
  assert.ok(moved('goalCover') > 0.1 * cells.length && moved('goalTrees') > 0.05 * cells.length, `${moved('goalCover')} cells with a cover increment and ${moved('goalTrees')} with a tree increment`);
  const net = (budget) => budget.goal, relative = Math.abs(net(gpu.land.budget) - net(cpu.land.budget)) / Math.max(1, Math.abs(net(cpu.land.budget)));
  assert.ok(relative < 1e-4, `the booked goal water differs by ${relative}: GPU ${net(gpu.land.budget)}, CPU ${net(cpu.land.budget)}`);
  const booked = net(gpu.land.budget);

  cpu.land.resetIncrements(); await gpu.land.resetIncrements();
  const reset = await gpu.land.serialize(), cpuReset = cpu.land.serialize();
  for (const name of FIELDS.slice(3)) { assert.ok(Array.from(reset[name]).every((v) => v === 0), `GPU ${name} zero`); assert.ok(Array.from(cpuReset[name]).every((v) => v === 0), `CPU ${name} zero`); }
  assert.equal(reset.soil[cells[0]], saved.soil[cells[0]], 'the stores are kept');
  assert.equal(gpu.land.budget.goal, booked, 'the booked water is kept');
  await gpu.step(900); await gpu.sync();
  const after = await gpu.land.serialize();
  assert.ok(cells.some((i) => after.goalAdded[i] > 0), 'the window starts again');
  gpu.destroy();
});

test('across the year end with the moisture gate on the trees, both engines agree on the stores and the four accumulators', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const land = { ...GOAL, treeMoisture: true }, start = (365 - 79) * 86400 - 8 * 900;
  const cpu = prepareLand(createModel(new Grid(16), options(land))), gpu = prepareLand(await createGpuModel(new Grid(16), options(land)));
  cpu.time = start; gpu.time = start;
  const C = cpu.mesh.nCells, cells = Array.from({ length: C }, (_, i) => i).filter((i) => cpu.geography.land[i] && !cpu.geography.iceSheet[i]);
  for (let n = 0; n < 16; n++) { cpu.step(900); await gpu.step(900); }
  await gpu.sync();
  const saved = await gpu.land.serialize(), cpuSaved = cpu.land.serialize();
  for (const name of FIELDS) {
    const a = cpuSaved[SAVED[name] ?? name], b = saved[SAVED[name] ?? name];
    let scale = 0, worst = 0;
    for (const i of cells) { scale = Math.max(scale, Math.abs(a[i])); worst = Math.max(worst, Math.abs(a[i] - b[i])); }
    assert.ok(worst < (name === 'goalTrees' ? 3e-3 : 2e-4) * Math.max(scale, 1), `${name} apart by ${worst} of ${scale}`);
  }
  assert.ok(cells.some((i) => saved.goalAdded[i] > 1e-3), 'the goal acted');
  gpu.destroy();
});

test('a save reloads the accumulators into the GPU and the window goes on from them', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const gpu = prepareLand(await createGpuModel(new Grid(16), options(GOAL)));
  gpu.time = START;
  for (let n = 0; n < 4; n++) await gpu.step(900);
  await gpu.sync();
  const first = await gpu.land.serialize();
  gpu.land.load(first);
  const second = await gpu.land.serialize();
  assert.equal(gpu.land.budget.goal, 0, 'the budget counts from the load');
  for (const name of FIELDS.slice(3)) assert.deepEqual(Array.from(second[name]), Array.from(first[name]), name);
  gpu.destroy();
});
