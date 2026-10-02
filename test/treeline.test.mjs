import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createGeography, syntheticTopography } from '../js/geography.module.js';
import { createLandSurface, sineSeason, insolationCycle, seasonEstimate, treelineFactor, SEASON_ESTIMATE } from '../js/physics/land.module.js';
import { MELTING_POINT } from '../js/physics/ice.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { encodeState, decodeState } from '../js/stateFile.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const DAY = 86400, YEAR = 365 * DAY;
const mesh = createModel(new Grid(8), { physics: false }).mesh;
const flat = (height = 100) => createGeography(mesh, syntheticTopography(90, 180, () => height), { landBridges: {}, seaStraits: {} });
const near = (a, b, tol = 1e-12) => Math.abs(a - b) < tol;
const deg = Math.PI / 180;

test('a sine year gives the share above the threshold and the mean excess over it in closed form, and the insolation its annual mean and harmonic', () => {
  const s = sineSeason(-4, 15, 0.9);
  assert.ok(near(s.length, 0.39407455499635624) && near(s.warmth, 2.58174392223435), JSON.stringify(s));
  assert.deepEqual(sineSeason(20, 5, 0.9), { length: 1, warmth: 19.1 });
  assert.deepEqual(sineSeason(-20, 5, 0.9), { length: 0, warmth: 0 });
  const equator = insolationCycle(0, 365);
  assert.ok(Math.abs(equator.mean - 415.58) < 0.02 && equator.amplitude < 1e-9, `equator ${JSON.stringify(equator)}: S₀/π J₀(ε) and no annual harmonic`);
  const north = insolationCycle(65 * deg), south = insolationCycle(-65 * deg);
  assert.ok(near(north.mean, south.mean, 1e-9) && near(north.amplitude, south.amplitude, 1e-9));
  const [a, b] = SEASON_ESTIMATE.mean, [k, c0, c1] = SEASON_ESTIMATE.amplitude;
  for (const [lat, z] of [[30, 0], [65, 0], [65, 1000], [10, 4000]]) {
    const q = insolationCycle(lat * deg), expected = sineSeason(a + b * q.mean - 0.0065 * z, Math.min(k * q.amplitude, c0 + c1 * q.amplitude), 0.9);
    assert.deepEqual(seasonEstimate(lat * deg, z), expected);
  }
  console.log(`estimated season at 0 m: ${[40, 50, 60, 65, 70, 75, 80].map((lat) => { const e = seasonEstimate(lat * deg); return `${lat}N ${Math.round(365 * e.length)} d at ${(0.9 + e.warmth / Math.max(e.length, 94 / 365)).toFixed(1)} °C`; }).join(', ')}`);
});

test('the treeline factor reads the season\'s mean over at least 94 days against 6.4–8.0 °C, and the season means follow the lowest air', () => {
  const land = createLandSurface(mesh, flat());
  land.initialize();
  const i = 5;
  land.seasonLength[i] = 0.4; land.seasonWarmth[i] = 2.8;
  assert.ok(near(land.treeFactor(i), 0.9374999999999994), `a 146-day season at 7.9 °C: ${land.treeFactor(i)}`);
  land.seasonLength[i] = 0.2; land.seasonWarmth[i] = 1.5;
  assert.ok(near(land.treeFactor(i), 0.20279255319148926), `a 73-day season counts as 94 days: ${land.treeFactor(i)}`);
  land.seasonLength[i] = 0; land.seasonWarmth[i] = 0;
  assert.equal(land.treeFactor(i), 0, 'no season, no trees');
  land.seasonLength[i] = 0.4; land.seasonWarmth[i] = 2.8;
  const surfaceT = new Float64Array(mesh.nCells).fill(285), flux = new Float64Array(mesh.nCells), half = 3 * YEAR * Math.LN2;
  land.update(i, surfaceT, flux, 0, half, MELTING_POINT + 10);
  assert.ok(near(land.seasonLength[i], 0.7) && near(land.seasonWarmth[i], 5.95), `a warm half-life: ${land.seasonLength[i]}, ${land.seasonWarmth[i]}`);
  land.update(i, surfaceT, flux, 0, half, MELTING_POINT - 3);
  assert.ok(near(land.seasonLength[i], 0.35) && near(land.seasonWarmth[i], 2.975), `a cold one: ${land.seasonLength[i]}, ${land.seasonWarmth[i]}`);
  land.update(i, surfaceT, flux, 0, half);
  assert.ok(near(land.seasonLength[i], 0.35), 'an update without the air leaves them');
});

test('trees grow toward the cover times the treeline factor over ten years of snow-free time, die back over three, and under snow hold unless the warmth fails', () => {
  const land = createLandSurface(mesh, flat());
  land.initialize();
  const i = 6, cap = land.capacity(i), flux = new Float64Array(mesh.nCells);
  const warm = new Float64Array(mesh.nCells).fill(295), cold = new Float64Array(mesh.nCells).fill(MELTING_POINT - 10);
  const set = (length, warmth, cover, trees, snow = 0) => { land.seasonLength[i] = length; land.seasonWarmth[i] = warmth; land.vegetation[i] = cover; land.canopy[i] = trees; land.snow[i] = snow; land.soil[i] = cap; land.surface[i] = 0; };
  set(1, 20, 1, 0.2);
  land.update(i, warm, flux, 0, 10 * YEAR);
  assert.ok(near(land.canopy[i], 0.7056964470628462), `growth: ${land.canopy[i]}`);
  set(1, 0, 1, 0.6);
  land.update(i, warm, flux, 0, 3 * YEAR);
  assert.ok(near(land.canopy[i], 0.2207276647028654), `die-back where the season is 0.9 °C: ${land.canopy[i]}`);
  set(0.5, 0.5 * (7.2 - 0.9), 1, 0.8, 100);
  assert.ok(near(land.treeFactor(i), 0.5));
  land.update(i, cold, flux, 0, 3 * YEAR);
  assert.ok(near(land.canopy[i], 0.6103638323514327), `under snow, down toward the factor: ${land.canopy[i]}`);
  set(1, 20, 0.8, 0.8, 100);
  land.update(i, cold, flux, 0, 200 * DAY);
  assert.ok(land.vegetation[i] < 0.7 && land.canopy[i] === 0.8, `the trees stand while the cover decays under snow: cover ${land.vegetation[i]}, trees ${land.canopy[i]}`);
  set(1, 20, 0.3, 0.8);
  land.soil[i] = 0.2 * cap;
  land.update(i, warm, flux, 0, 3 * YEAR);
  assert.ok(land.canopy[i] < 0.8 * Math.exp(-1) + 0.3, `a drought thins them toward the cover: ${land.canopy[i]}`);
  const plain = createLandSurface(mesh, flat(), { treeline: false });
  plain.initialize();
  plain.vegetation[i] = 0.6; plain.canopy[i] = 0.2; plain.soil[i] = cap; plain.seasonLength[i] = 0; plain.seasonWarmth[i] = 0;
  plain.update(i, warm, flux, 0, DAY);
  assert.ok(plain.canopy[i] >= plain.vegetation[i], 'without the treeline the standing cover follows the cover up at once');
});

test('a fresh start and an older state start the trees at the cover times the estimated season\'s factor, a saved state keeps its season means and trees, and the ice sheets grow none', async () => {
  const topography = syntheticTopography(180, 360, (lat, lon) => (lat < -1.2 ? 2000 : Math.cos(lon) > 0 ? 100 : -4000));
  const geography = createGeography(mesh, topography, { landBridges: {}, seaStraits: {} });
  const land = createLandSurface(mesh, geography), C = mesh.nCells;
  land.initialize();
  let polar = 0, temperate = 0, sheet = 0;
  for (let i = 0; i < C; i++) {
    if (!geography.land[i]) { assert.equal(land.seasonLength[i], 0); assert.equal(land.canopy[i], 0); continue; }
    const e = seasonEstimate(mesh.latCell[i], geography.elevation[i]);
    assert.ok(near(land.seasonLength[i], e.length) && near(land.seasonWarmth[i], e.warmth), `cell ${i}`);
    if (geography.iceSheet[i]) { sheet++; assert.equal(land.canopy[i], 0); continue; }
    assert.ok(near(land.canopy[i], 0.5 * land.treeFactor(i)), `cell ${i}: trees ${land.canopy[i]}`);
    const lat = Math.abs(mesh.latCell[i]) / deg;
    if (lat > 78) { polar++; assert.equal(land.canopy[i], 0, `no trees at ${lat.toFixed(1)}°`); }
    if (lat < 45) { temperate++; assert.equal(land.canopy[i], 0.5, `full potential at ${lat.toFixed(1)}°`); }
  }
  assert.ok(polar >= 2 && temperate > 50 && sheet > 3, `${polar} polar, ${temperate} temperate, ${sheet} ice-sheet cells`);
  const vegetation = Float64Array.from({ length: C }, (_, i) => (i % 5) / 5), snow = Float64Array.from({ length: C }, (_, i) => (i % 3 === 0 ? 30 : 0));
  const older = { soil: new Float64Array(C).fill(50), snow, vegetation, canopy: Float64Array.from(vegetation, (v) => Math.min(1, v + 0.25)) };
  land.load(older);
  for (let i = 0; i < C; i++) if (geography.land[i] && !geography.iceSheet[i]) assert.ok(near(land.canopy[i], land.treeFactor(i) * vegetation[i]), `an older state's standing cover is not a tree cover: cell ${i}`);
  const seasonLength = Float64Array.from({ length: C }, (_, i) => (i % 7) / 7), seasonWarmth = Float64Array.from({ length: C }, (_, i) => (i % 11) / 2);
  const canopy = Float64Array.from({ length: C }, (_, i) => (i % 13) / 13);
  land.load({ ...older, seasonLength, seasonWarmth, canopy });
  const saved = await decodeState(encodeState({ N: 8, K: 1, day: 0, time: 0, terrain: true, land: land.serialize() }));
  const back = createLandSurface(mesh, geography);
  back.load(saved.land);
  for (let i = 0; i < C; i++) {
    const kept = geography.land[i] ? 1 : 0;
    assert.ok(near(back.seasonLength[i], kept * Math.fround(seasonLength[i]), 1e-7) && near(back.seasonWarmth[i], kept * Math.fround(seasonWarmth[i]), 1e-6), `cell ${i}: season`);
    assert.ok(near(back.canopy[i], kept && !geography.iceSheet[i] ? Math.fround(canopy[i]) : 0, 1e-7), `cell ${i}: trees ${back.canopy[i]}`);
  }
});

test('the model hands the land its lowest air temperature each step', () => {
  const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 300 : -4000));
  const model = createModel(new Grid(6), { topography, land: { seasonMemory: 2 * 900 } });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i]) model.state[6][i] = 0;
  model.land.initialize();
  const { K, C, exnerLayer } = model.core.diagnostics, start = Float64Array.from(model.land.seasonWarmth);
  model.step(900);
  let toward = 0, away = 0;
  for (let i = 0; i < C; i++) {
    if (!model.geography.land[i]) continue;
    const excess = Math.max(0, model.state[1][(K - 1) * C + i] * exnerLayer[(K - 1) * C + i] - MELTING_POINT - 0.9), moved = model.land.seasonWarmth[i] - start[i];
    if (Math.abs(excess - start[i]) < 2) continue;
    if (Math.sign(moved) === Math.sign(excess - start[i]) && Math.abs(moved) > 0.2 * Math.abs(excess - start[i])) toward++; else away++;
  }
  assert.ok(toward > 50 && away === 0, `${toward} land cells moved toward their air, ${away} did not`);
});

const topography = syntheticTopography(90, 180, (lat, lon) => ((Math.cos(lon) > 0 && Math.abs(lat) < 1.2) || lat < -1.15 ? 300 : -4000));

test('over 48 GPU steps the season means and the tree cover evolve as on the CPU', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const land = { seasonMemory: 6 * 3600, treeGrowthTime: 3 * 3600, treeDeclineTime: 2 * 3600, treelineWarmth: [6, 22], growthTime: 4 * 3600, declineTime: 3 * 3600, snowDeclineTime: 4 * 3600 };
  const prepare = (model) => {
    const C = model.mesh.nCells, init = initializeState(model, {});
    for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
    for (let i = 0; i < C; i++) if (model.geography.land[i]) model.state[6][i] = 0;
    model.seaIce.load(model.state[6], null);
    if (model.load) model.load();
    model.ocean.initialize(model.state[3], model.state[6]);
    let s = 3;
    const rnd = () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; };
    const snow = Float64Array.from({ length: C }, (_, i) => (Math.abs(model.mesh.latCell[i]) > 0.9 ? 30 : 0));
    const vegetation = Float64Array.from({ length: C }, () => rnd()), canopy = Float64Array.from({ length: C }, () => rnd());
    const seasonLength = Float64Array.from({ length: C }, () => rnd()), seasonWarmth = Float64Array.from({ length: C }, () => 10 * rnd());
    model.land.load({ soil: Float64Array.from({ length: C }, () => 300 * rnd()), snow, vegetation, canopy, seasonLength, seasonWarmth }, model.state[6]);
    return model;
  };
  const cpu = prepare(createModel(new Grid(6), { topography, land }));
  const gpu = prepare(await createGpuModel(new Grid(6), { topography, land }));
  const C = cpu.mesh.nCells, before = Float64Array.from(cpu.land.canopy), lengthBefore = Float64Array.from(cpu.land.seasonLength);
  for (let n = 0; n < 48; n++) { cpu.step(900); await gpu.step(900); }
  await gpu.sync();
  const saved = await gpu.land.serialize();
  const { K, exnerLayer } = cpu.core.diagnostics;
  let worstLength = 0, worstWarmth = 0, worstTrees = 0, worstAir = 0, grew = 0, died = 0, snowy = 0, inSeason = 0, cells = 0;
  for (let i = 0; i < C; i++) {
    if (!cpu.geography.land[i] || cpu.geography.iceSheet[i]) continue;
    cells++;
    worstAir = Math.max(worstAir, Math.abs(cpu.state[1][(K - 1) * C + i] - gpu.state[1][(K - 1) * C + i]) * exnerLayer[(K - 1) * C + i]);
    worstLength = Math.max(worstLength, Math.abs(cpu.land.seasonLength[i] - saved.seasonLength[i]));
    worstWarmth = Math.max(worstWarmth, Math.abs(cpu.land.seasonWarmth[i] - saved.seasonWarmth[i]));
    worstTrees = Math.max(worstTrees, Math.abs(cpu.land.canopy[i] - saved.canopy[i]));
    if (cpu.land.canopy[i] > before[i] + 0.05) grew++;
    if (cpu.land.canopy[i] < before[i] - 0.05) died++;
    if (cpu.land.snow[i] > 0) snowy++;
    if (cpu.land.seasonLength[i] > lengthBefore[i]) inSeason++;
  }
  gpu.destroy();
  console.log(`N=6, 48 steps over ${cells} land cells (${snowy} under snow, ${inSeason} in season): trees grew on ${grew} and died back on ${died}; engines apart by ${worstLength.toExponential(1)} in season length, ${worstWarmth.toExponential(1)} K in season warmth (the lowest air by ${worstAir.toExponential(1)} K), ${worstTrees.toExponential(1)} in tree cover`);
  assert.ok(grew > 10 && died > 10 && snowy > 5 && inSeason > 10 && cells - inSeason > 5, `${grew} grew, ${died} died, ${snowy} snowy, ${inSeason} in season`);
  assert.ok(worstLength < 1e-4 && worstWarmth <= worstAir && worstTrees < 1e-4, `season length ${worstLength}, warmth ${worstWarmth} under air ${worstAir}, trees ${worstTrees}`);
});

test('land regridding carries the season means and fills land the source lacks from the estimate, its trees at the guessed cover times the estimate\'s factor', async () => {
  const { regridLand } = await import('../js/physics/regrid.module.js');
  const make = (N, relief) => { const m = createModel(new Grid(N), { physics: false, topography: relief }); return { mesh: m.mesh, geography: createGeography(m.mesh, relief, { landBridges: {}, seaStraits: {} }) }; };
  const islands = syntheticTopography(90, 180, (lat, lon) => (lat > 1.3 && Math.cos(lon) < 0 ? 50 : 0) || (((Math.cos(lon) > 0 && Math.abs(lat) < 1.2) || lat < -1.15) ? 300 : -4000));
  const source = make(6, topography), target = make(8, islands), C = source.mesh.nCells;
  const land = { soil: new Float64Array(C).fill(100), snow: new Float64Array(C), vegetation: new Float64Array(C).fill(0.5), canopy: new Float64Array(C).fill(0.3), seasonLength: new Float64Array(C).fill(0.45), seasonWarmth: new Float64Array(C).fill(4) };
  const out = regridLand(source, target, land);
  let carried = 0, filled = 0, partial = 0;
  for (let n = 0; n < target.mesh.nCells; n++) {
    if (!target.geography.land[n]) { assert.equal(out.seasonLength[n], 0); assert.equal(out.seasonWarmth[n], 0); continue; }
    if (out.seasonLength[n] === 0.45 && out.seasonWarmth[n] === 4) { carried++; continue; }
    assert.ok(out.seasonLength[n] > 0 && out.seasonLength[n] <= 1, `cell ${n}: estimated ${out.seasonLength[n]}`);
    assert.ok(near(out.canopy[n], 0.5 * treelineFactor(out.seasonLength[n], out.seasonWarmth[n])), `cell ${n}: trees ${out.canopy[n]} on new land`);
    filled++; if (out.canopy[n] < 0.5) partial++;
  }
  assert.ok(carried > 0.8 * [...target.geography.land].filter(Boolean).length && filled > 5 && partial > 5, `${carried} carried, ${filled} filled, ${partial} of them under the full cover`);
  const same = regridLand(source, source, land);
  assert.deepEqual([...same.seasonWarmth], [...land.seasonWarmth]);
  assert.equal(regridLand(source, target, { soil: land.soil, snow: land.snow, vegetation: land.vegetation }).seasonLength, undefined);
});
