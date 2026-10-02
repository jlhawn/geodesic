import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createGeography, syntheticTopography } from '../js/geography.module.js';
import { createLandSurface, aridityFactor, moistureEstimate, insolationCycle, treelineFactor, FOREST_ARIDITY, MOISTURE_ESTIMATE } from '../js/physics/land.module.js';
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

test('the moisture factor ramps with the aridity index P/PET between the forest limits, and the moisture means follow the rain and the potential evaporation', () => {
  const [dry, wet] = FOREST_ARIDITY;
  assert.equal(aridityFactor(0.5, 5), Math.max(0, Math.min(1, (0.1 - dry) / (wet - dry))));
  assert.ok(near(aridityFactor(1.2, 3, [0.2, 1.0]), 0.25), 'P/PET 0.4 a quarter of the way from 0.2 to 1.0');
  assert.equal(aridityFactor(3, 3, [0.2, 1.0]), 1);
  assert.equal(aridityFactor(0.1, 3, [0.2, 1.0]), 0);
  assert.equal(aridityFactor(0.5, 0, [0.2, 1.0]), 1, 'no evaporative demand, no moisture limit');
  const equator = moistureEstimate(0, 0.5, 0.5), q = insolationCycle(0).mean;
  assert.ok(near(equator.demand, -2.14 + 0.0134 * q) && near(equator.rain, equator.demand * (0.01 + 0.79 * 0.5 + 0.63 * 0.5)), `the estimate at the equator: ${JSON.stringify(equator)}`);
  assert.ok(Math.abs(equator.demand - 3.43) < 0.01 && Math.abs(equator.rain - 2.47) < 0.01, 'about 1250 mm/yr of demand and P/PET 0.72 half full and half covered');
  assert.deepEqual(MOISTURE_ESTIMATE, { demand: [-2.14, 0.0134], aridity: [0.01, 0.79, 0.63] });
  const land = createLandSurface(mesh, flat());
  land.initialize();
  const i = 7, surfaceT = new Float64Array(mesh.nCells).fill(290), flux = new Float64Array(mesh.nCells), half = 3 * YEAR * Math.LN2;
  land.rainMean[i] = 1; land.demandMean[i] = 4;
  land.deposit(i, 2 * half / DAY, 290, half);
  assert.ok(near(land.rainMean[i], 1.5), `rain halfway to 2 mm/d over a half-life: ${land.rainMean[i]}`);
  land.deposit(i, 5, 290);
  assert.ok(near(land.rainMean[i], 1.5), 'a deposit without its step leaves the mean');
  land.deposit(i, 6 * half / DAY, MELTING_POINT - 5, half);
  assert.ok(near(land.rainMean[i], 3.75), `snowfall counts: ${land.rainMean[i]}`);
  land.update(i, surfaceT, flux, 0, half, null, 6 / DAY);
  assert.ok(near(land.demandMean[i], 5), `the potential evaporation (6 mm/d) halfway: ${land.demandMean[i]}`);
  land.update(i, surfaceT, flux, 0, half);
  assert.ok(near(land.demandMean[i], 5), 'an update without it leaves the mean');
});

test('trees grow toward the treeline factor times the moisture factor times the cover, die back toward it where the climate dries, and the rest of the cover is grass', () => {
  const land = createLandSurface(mesh, flat(), { forestAridity: [0.2, 1.0] });
  land.initialize();
  const i = 6, cap = land.capacity(i), flux = new Float64Array(mesh.nCells), warm = new Float64Array(mesh.nCells).fill(295), cold = new Float64Array(mesh.nCells).fill(MELTING_POINT - 10);
  const set = (rain, demand, cover, trees, snow = 0) => { land.seasonLength[i] = 1; land.seasonWarmth[i] = 20; land.rainMean[i] = rain; land.demandMean[i] = demand; land.vegetation[i] = cover; land.canopy[i] = trees; land.snow[i] = snow; land.soil[i] = cap; land.surface[i] = 0; };
  set(1.8, 3, 1, 0.1);
  assert.ok(near(land.moistureFactor(i), 0.5) && near(land.treeFactor(i), 0.5), `P/PET 0.6: ${land.moistureFactor(i)}`);
  land.update(i, warm, flux, 0, 10 * YEAR);
  assert.ok(near(land.canopy[i], 0.5 - 0.4 * Math.exp(-1)), `growth toward half: ${land.canopy[i]}`);
  set(0.3, 3, 1, 0.8);
  land.update(i, warm, flux, 0, 3 * YEAR);
  assert.ok(near(land.canopy[i], 0.8 * Math.exp(-1)), `die-back in the arid class: ${land.canopy[i]}`);
  set(0.9, 3, 1, 0.8, 100);
  land.update(i, cold, flux, 0, 3 * YEAR);
  assert.ok(near(land.canopy[i], 0.125 + 0.675 * Math.exp(-1)), `under snow down toward the potential 0.125 at P/PET 0.3: ${land.canopy[i]}`);
  const wet = createLandSurface(mesh, flat(), { treeMoisture: false });
  wet.initialize();
  wet.seasonLength[i] = 1; wet.seasonWarmth[i] = 20; wet.rainMean[i] = 0; wet.demandMean[i] = 5;
  assert.equal(wet.treeFactor(i), 1, 'without the moisture gate the season alone sets the potential');
});

test('the vegetated albedo blends forest and grass by their shares of the cover, snow buries grass less 0.06 times its share and trees mask it as before', () => {
  const land = createLandSurface(mesh, flat(), { soilCarbon: false });
  land.initialize();
  const i = 3;
  land.vegetation[i] = 0.6; land.canopy[i] = 0.2; land.surface[i] = 0; land.snow[i] = 0;
  const cover = 0.20 + (0.13 - 0.20) / 3;
  assert.ok(near(land.albedo(i), 0.30 + (cover - 0.30) * 0.6), `a third forest: ${land.albedo(i)}`);
  land.canopy[i] = 0.6;
  assert.ok(near(land.albedo(i), 0.30 + (0.13 - 0.30) * 0.6), 'all forest');
  land.canopy[i] = 0.9;
  assert.ok(near(land.albedo(i), 0.30 + (0.13 - 0.30) * 0.6), 'trees standing above the cover read as forest');
  land.canopy[i] = 0;
  assert.ok(near(land.albedo(i), 0.30 + (0.20 - 0.30) * 0.6), 'all grass');
  land.vegetation[i] = 0;
  assert.ok(near(land.albedo(i), 0.30), 'bare');
  land.vegetation[i] = 0.6; land.canopy[i] = 0.2; land.snow[i] = 40; land.snowAlbedo[i] = 0.85;
  const own = 0.85 - 0.06 * 0.4;
  assert.ok(near(land.albedo(i), own + (0.27 - own) * 0.2 / 0.7), `snow on grass under sparse trees: ${land.albedo(i)}`);
  land.canopy[i] = 0;
  assert.ok(near(land.albedo(i), 0.85 - 0.06 * 0.6), `snow on grassland: ${land.albedo(i)}`);
  land.canopy[i] = 0.6;
  assert.ok(near(land.albedo(i), 0.85 + (0.27 - 0.85) * 0.6 / 0.7), 'snow under trees alone, as before');
  const plain = createLandSurface(mesh, flat(), { grassland: false, soilCarbon: false });
  plain.initialize();
  plain.vegetation[i] = 0.6; plain.canopy[i] = 0.2; plain.surface[i] = 0; plain.snow[i] = 0;
  assert.ok(near(plain.albedo(i), 0.30 + (0.13 - 0.30) * 0.6), 'without grassland one vegetated albedo');
});

test('a fresh start and an older state take the moisture means from the estimate and their trees at most at the cover times the potential; a saved state keeps them through a state file', async () => {
  const topography = syntheticTopography(180, 360, (lat, lon) => (lat < -1.2 ? 2000 : Math.cos(lon) > 0 ? 100 : -4000));
  const geography = createGeography(mesh, topography, { landBridges: {}, seaStraits: {} });
  const land = createLandSurface(mesh, geography), C = mesh.nCells;
  land.initialize();
  for (let i = 0; i < C; i++) {
    if (!geography.land[i]) { assert.equal(land.rainMean[i], 0); assert.equal(land.demandMean[i], 0); continue; }
    const e = moistureEstimate(mesh.latCell[i], 0.5, geography.iceSheet[i] ? 0 : 0.5);
    assert.ok(near(land.rainMean[i], e.rain) && near(land.demandMean[i], e.demand), `cell ${i}`);
    if (!geography.iceSheet[i]) assert.ok(near(land.canopy[i], 0.5 * land.treeFactor(i)), `cell ${i}: trees ${land.canopy[i]}`);
  }
  const soil = Float64Array.from({ length: C }, (_, i) => 300 * ((i % 9) / 8)), vegetation = Float64Array.from({ length: C }, (_, i) => (i % 5) / 5), snow = new Float64Array(C);
  const seasonLength = new Float64Array(C).fill(0.6), seasonWarmth = new Float64Array(C).fill(8);
  land.load({ soil, snow, vegetation, canopy: new Float64Array(C).fill(0.9), seasonLength, seasonWarmth });
  let gated = 0;
  for (let i = 0; i < C; i++) {
    if (!geography.land[i] || geography.iceSheet[i]) continue;
    const e = moistureEstimate(mesh.latCell[i], soil[i] / 300, vegetation[i]);
    assert.ok(near(land.rainMean[i], e.rain) && near(land.demandMean[i], e.demand), `cell ${i}: estimated from its own fill`);
    assert.ok(near(land.canopy[i], vegetation[i] * treelineFactor(0.6, 8) * aridityFactor(e.rain, e.demand)), `cell ${i}: an older state's trees restart at the potential`);
    if (aridityFactor(e.rain, e.demand) < 1 && vegetation[i] > 0) gated++;
  }
  assert.ok(gated > 20, `${gated} cells held below the season's potential by their moisture`);
  const young = Float64Array.from({ length: C }, (_, i) => (i % 3) / 20);
  land.load({ soil, snow, vegetation, canopy: young, seasonLength, seasonWarmth });
  let kept = 0;
  for (let i = 0; i < C; i++) {
    if (!geography.land[i] || geography.iceSheet[i]) continue;
    const e = moistureEstimate(mesh.latCell[i], soil[i] / 300, vegetation[i]), potential = vegetation[i] * treelineFactor(0.6, 8) * aridityFactor(e.rain, e.demand);
    assert.ok(near(land.canopy[i], Math.min(young[i], potential)), `cell ${i}: an older state's trees below the potential stand`);
    if (young[i] < potential) kept++;
  }
  assert.ok(kept > 20, `${kept} cells keep their younger trees`);
  const rainMean = Float64Array.from({ length: C }, (_, i) => (i % 7) / 2), demandMean = Float64Array.from({ length: C }, (_, i) => 1 + (i % 11) / 3), canopy = Float64Array.from({ length: C }, (_, i) => (i % 13) / 13);
  land.load({ soil, snow, vegetation, canopy, seasonLength, seasonWarmth, rainMean, demandMean });
  const saved = await decodeState(encodeState({ N: 8, K: 1, day: 0, time: 0, terrain: true, land: land.serialize() }));
  const back = createLandSurface(mesh, geography);
  back.load(saved.land);
  for (let i = 0; i < C; i++) {
    const kept = geography.land[i] ? 1 : 0;
    assert.ok(near(back.rainMean[i], kept * Math.fround(rainMean[i]), 1e-6) && near(back.demandMean[i], kept * Math.fround(demandMean[i]), 1e-6), `cell ${i}: moisture means`);
    assert.ok(near(back.canopy[i], kept && !geography.iceSheet[i] ? Math.fround(canopy[i]) : 0, 1e-7), `cell ${i}: trees kept`);
  }
});

test('the model hands the land each step\'s rain and potential evaporation', () => {
  const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 300 : -4000));
  const model = createModel(new Grid(6), { topography, land: { moistureMemory: 900 } });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i]) model.state[6][i] = 0;
  model.land.initialize();
  for (let n = 0; n < 4; n++) model.step(900);
  const keep = 1 - Math.exp(-1);
  let checked = 0, rained = 0;
  for (let i = 0; i < model.mesh.nCells; i++) {
    if (!model.geography.land[i]) continue;
    const before = model.land.rainMean[i];
    model.land.deposit(i, 0, 290, 900);
    assert.ok(near(model.land.rainMean[i], before * (1 - keep), 1e-9));
    if (before > 0.1) rained++;
    if (model.land.demandMean[i] > 0.1) checked++;
  }
  assert.ok(rained > 10 && checked > 50, `${rained} land cells with rain, ${checked} with a potential evaporation`);
});

test('land regridding carries the moisture means and fills new land from the estimate', async () => {
  const { regridLand } = await import('../js/physics/regrid.module.js');
  const make = (N, relief) => { const m = createModel(new Grid(N), { physics: false, topography: relief }); return { mesh: m.mesh, geography: createGeography(m.mesh, relief, { landBridges: {}, seaStraits: {} }) }; };
  const continent = syntheticTopography(90, 180, (lat, lon) => ((Math.cos(lon) > 0 && Math.abs(lat) < 1.2) || lat < -1.15 ? 300 : -4000));
  const islands = syntheticTopography(90, 180, (lat, lon) => (lat > 1.3 && Math.cos(lon) < 0 ? 50 : 0) || (((Math.cos(lon) > 0 && Math.abs(lat) < 1.2) || lat < -1.15) ? 300 : -4000));
  const source = make(6, continent), target = make(8, islands), C = source.mesh.nCells;
  const land = { soil: new Float64Array(C).fill(150), snow: new Float64Array(C), vegetation: new Float64Array(C).fill(0.5), canopy: new Float64Array(C).fill(0.3), seasonLength: new Float64Array(C).fill(0.45), seasonWarmth: new Float64Array(C).fill(4), rainMean: new Float64Array(C).fill(2.5), demandMean: new Float64Array(C).fill(3.5) };
  const out = regridLand(source, target, land);
  let carried = 0, filled = 0;
  for (let n = 0; n < target.mesh.nCells; n++) {
    if (!target.geography.land[n]) { assert.equal(out.rainMean[n], 0); assert.equal(out.demandMean[n], 0); continue; }
    if (out.rainMean[n] === 2.5 && out.demandMean[n] === 3.5) { carried++; continue; }
    assert.ok(out.demandMean[n] > 0, `cell ${n}: estimated demand ${out.demandMean[n]}`);
    filled++;
  }
  assert.ok(carried > 0.8 * [...target.geography.land].filter(Boolean).length && filled > 5, `${carried} carried, ${filled} filled`);
  assert.equal(regridLand(source, target, { soil: land.soil, snow: land.snow, vegetation: land.vegetation }).rainMean, undefined);
});

const topography = syntheticTopography(90, 180, (lat, lon) => ((Math.cos(lon) > 0 && Math.abs(lat) < 1.2) || lat < -1.15 ? 300 : -4000));

test('both engines give land the same albedo from random cover, trees, snow and snow albedo, forest and grass blended', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const cpu = createModel(new Grid(6), { topography }), gpu = await createGpuModel(new Grid(6), { topography });
  const C = cpu.mesh.nCells, landMask = cpu.geography.land;
  let s = 5;
  const rnd = () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; };
  const init = initializeState(cpu, {});
  for (let i = 0; i < C; i++) if (landMask[i]) { init[6][i] = 0; init[3][i] = MELTING_POINT - 15 + 20 * rnd(); }
  const soil = Float64Array.from({ length: C }, (_, i) => (landMask[i] ? 300 * rnd() : 0)), surface = Float64Array.from({ length: C }, () => 15 * rnd());
  const snow = Float64Array.from({ length: C }, () => (rnd() < 0.4 ? 40 * rnd() : 0));
  const vegetation = Float64Array.from({ length: C }, () => rnd()), canopy = Float64Array.from(vegetation, (v) => Math.min(1, v * 1.3 * rnd()));
  const snowAlbedo = Float64Array.from({ length: C }, () => 0.5 + 0.35 * rnd());
  const seasonLength = new Float64Array(C).fill(0.6), seasonWarmth = new Float64Array(C).fill(8), rainMean = Float64Array.from({ length: C }, () => 4 * rnd()), demandMean = Float64Array.from({ length: C }, () => 1 + 4 * rnd());
  for (const m of [cpu, gpu]) {
    for (let a = 0; a < init.length; a++) m.state[a].set(init[a]);
    m.seaIce.load(m.state[6], null);
    m.land.load({ soil, surface, snow, vegetation, canopy, snowAlbedo, seasonLength, seasonWarmth, rainMean, demandMean }, m.state[6]);
  }
  gpu.load();
  const expected = Float64Array.from({ length: C }, (_, i) => (landMask[i] ? cpu.land.albedo(i) : NaN));
  const plain = createModel(new Grid(6), { topography, land: { grassland: false } });
  plain.land.load({ soil, surface, snow, vegetation, canopy, snowAlbedo, seasonLength, seasonWarmth, rainMean, demandMean }, init[6]);
  await gpu.step(900);
  const ph = await gpu.gpu.downloadPhysics();
  let worst = 0, cells = 0, grassy = 0, buried = 0, most = 0;
  for (let i = 0; i < C; i++) {
    if (!landMask[i] || cpu.geography.iceSheet[i]) continue;
    cells++;
    worst = Math.max(worst, Math.abs(ph.ADIF[i] - expected[i]));
    if (vegetation[i] - cpu.land.canopy[i] > 0.2) { grassy++; if (snow[i] > 0) buried++; }
    most = Math.max(most, Math.abs(expected[i] - plain.land.albedo(i)));
  }
  gpu.destroy();
  console.log(`N=6: ${cells} land cells (${grassy} with grass above 0.2, ${buried} of them under snow): albedo engines apart by ${worst.toExponential(1)}; grass moves it by up to ${most.toFixed(3)}`);
  assert.ok(grassy > 20 && buried > 5 && most > 0.03, `${grassy} grassy, ${buried} buried, most ${most}`);
  assert.ok(worst < 1e-6, `engines apart by ${worst}`);
});

test('over 48 GPU steps the moisture means and the gated trees evolve as on the CPU', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const land = { seasonMemory: 6 * 3600, moistureMemory: 6 * 3600, treeGrowthTime: 3 * 3600, treeDeclineTime: 2 * 3600, treelineWarmth: [6, 22], forestAridity: [0.2, 1.0], growthTime: 4 * 3600, declineTime: 3 * 3600, snowDeclineTime: 4 * 3600 };
  const prepare = (model) => {
    const C = model.mesh.nCells, init = initializeState(model, {});
    for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
    for (let i = 0; i < C; i++) if (model.geography.land[i]) model.state[6][i] = 0;
    model.seaIce.load(model.state[6], null);
    if (model.load) model.load();
    model.ocean.initialize(model.state[3], model.state[6]);
    let s = 9;
    const rnd = () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; };
    const snow = Float64Array.from({ length: C }, (_, i) => (Math.abs(model.mesh.latCell[i]) > 0.9 ? 30 : 0));
    const vegetation = Float64Array.from({ length: C }, () => rnd()), canopy = Float64Array.from({ length: C }, () => rnd());
    const seasonLength = Float64Array.from({ length: C }, () => 0.5 + 0.5 * rnd()), seasonWarmth = Float64Array.from({ length: C }, () => 10 * rnd());
    const rainMean = Float64Array.from({ length: C }, () => 3 * rnd()), demandMean = Float64Array.from({ length: C }, () => 0.5 + 4 * rnd());
    model.land.load({ soil: Float64Array.from({ length: C }, () => 300 * rnd()), snow, vegetation, canopy, seasonLength, seasonWarmth, rainMean, demandMean }, model.state[6]);
    return model;
  };
  const cpu = prepare(createModel(new Grid(6), { topography, land }));
  const gpu = prepare(await createGpuModel(new Grid(6), { topography, land }));
  const C = cpu.mesh.nCells, rainBefore = Float64Array.from(cpu.land.rainMean), demandBefore = Float64Array.from(cpu.land.demandMean), treesBefore = Float64Array.from(cpu.land.canopy);
  for (let n = 0; n < 48; n++) { cpu.step(900); await gpu.step(900); }
  await gpu.sync();
  const saved = await gpu.land.serialize();
  let worstRain = 0, worstDemand = 0, worstTrees = 0, rainMoved = 0, demandMoved = 0, gated = 0, scaleRain = 0, scaleDemand = 0, cells = 0;
  for (let i = 0; i < C; i++) {
    if (!cpu.geography.land[i] || cpu.geography.iceSheet[i]) continue;
    cells++;
    worstRain = Math.max(worstRain, Math.abs(cpu.land.rainMean[i] - saved.rainMean[i]));
    worstDemand = Math.max(worstDemand, Math.abs(cpu.land.demandMean[i] - saved.demandMean[i]));
    worstTrees = Math.max(worstTrees, Math.abs(cpu.land.canopy[i] - saved.canopy[i]));
    scaleRain = Math.max(scaleRain, cpu.land.rainMean[i]); scaleDemand = Math.max(scaleDemand, cpu.land.demandMean[i]);
    if (Math.abs(cpu.land.rainMean[i] - rainBefore[i]) > 0.3) rainMoved++;
    if (Math.abs(cpu.land.demandMean[i] - demandBefore[i]) > 0.3) demandMoved++;
    if (cpu.land.moistureFactor(i) < 0.9 && cpu.land.canopy[i] < treesBefore[i] - 0.05) gated++;
  }
  gpu.destroy();
  console.log(`N=6, 48 steps over ${cells} land cells: rain means moved on ${rainMoved} (to ${scaleRain.toFixed(1)} mm/d), demand on ${demandMoved} (to ${scaleDemand.toFixed(1)} mm/d), trees thinned by moisture on ${gated}; engines apart by ${worstRain.toExponential(1)} and ${worstDemand.toExponential(1)} mm/d and ${worstTrees.toExponential(1)} in tree cover`);
  assert.ok(rainMoved > 10 && demandMoved > 10 && gated > 10, `${rainMoved}, ${demandMoved}, ${gated}`);
  assert.ok(worstRain < 0.1 * scaleRain && worstDemand < 0.1 * scaleDemand && worstTrees < 0.02, `rain ${worstRain}, demand ${worstDemand}, trees ${worstTrees}`);
});
