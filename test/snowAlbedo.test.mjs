import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createGeography, syntheticTopography } from '../js/geography.module.js';
import { createLandSurface } from '../js/physics/land.module.js';
import { createSeaIce, agedSnowAlbedo, refreshedSnowAlbedo, SNOW_AGEING, MELTING_POINT, FREEZING_POINT } from '../js/physics/ice.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { encodeState, decodeState } from '../js/stateFile.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const DAY = 86400;
const mesh = createModel(new Grid(8), { physics: false }).mesh;
const flat = () => createGeography(mesh, syntheticTopography(90, 180, () => 100), { landBridges: {}, seaStraits: {} });
const near = (a, b, tol = 1e-12) => Math.abs(a - b) < tol;

test('snow ages as Douville et al. (1995): cold snow loses 0.008 a day slowed by the cold, wet snow relaxes to its floor at 0.24 a day, and snowfall refreshes it', () => {
  const plain = { ...SNOW_AGEING, ageingActivation: 0 };
  assert.ok(near(agedSnowAlbedo(0.85, MELTING_POINT - 20, DAY, 0.5, plain), 0.842));
  assert.ok(near(agedSnowAlbedo(0.85, MELTING_POINT - 20, DAY, 0.5), 0.8481162461858663), 'at -20 C grain growth runs at exp(-1.446) = 0.235 of its pace at the melting point');
  assert.ok(near(agedSnowAlbedo(0.80, MELTING_POINT - 10, 2 * DAY, 0.5), 0.792019673307976), 'at -10 C at 0.499 of it');
  assert.ok(near(agedSnowAlbedo(0.505, MELTING_POINT - 30, 30 * DAY, 0.5, plain), 0.5), 'not below the floor');
  assert.ok(near(agedSnowAlbedo(0.85, MELTING_POINT, DAY, 0.5), 0.7753197513732937));
  assert.ok(near(agedSnowAlbedo(0.85, MELTING_POINT - 1.5, DAY, 0.5), 0.7753197513732937), 'within 2 K of melting the snow is wet');
  assert.ok(near(agedSnowAlbedo(0.85, MELTING_POINT, DAY, 0.7), 0.817994179159983), 'toward the floor it is given');
  assert.ok(near(refreshedSnowAlbedo(0.6, 5), 0.725) && near(refreshedSnowAlbedo(0.6, 20), 0.85) && refreshedSnowAlbedo(0.6, 0) === 0.6);
});

test('the land ages its snow after each update, refreshes it with snowfall and starts the next snow fresh once the ground is bare', () => {
  const land = createLandSurface(mesh, flat());
  land.initialize();
  const i = 3, surfaceT = new Float64Array(mesh.nCells).fill(MELTING_POINT - 20), flux = new Float64Array(mesh.nCells);
  assert.equal(land.snowAlbedo[i], 0.85);
  land.deposit(i, 50, MELTING_POINT - 20);
  land.update(i, surfaceT, flux, 0, DAY);
  assert.ok(near(land.snowAlbedo[i], 0.8481162461858663), `cold day ${land.snowAlbedo[i]}`);
  surfaceT.fill(MELTING_POINT - 1);
  land.snowAlbedo[i] = 0.85;
  land.update(i, surfaceT, flux, 0, DAY);
  assert.ok(near(land.snowAlbedo[i], 0.7753197513732937), `wet day ${land.snowAlbedo[i]}`);
  land.deposit(i, 5, MELTING_POINT - 5);
  assert.ok(near(land.snowAlbedo[i], 0.7753197513732937 + 0.5 * (0.85 - 0.7753197513732937)));
  land.deposit(i, 5, MELTING_POINT + 5);
  assert.ok(near(land.snowAlbedo[i], 0.7753197513732937 + 0.5 * (0.85 - 0.7753197513732937)), 'rain does not refresh it');
  land.snowAlbedo[i] = 0.6;
  flux.fill(5000);
  land.update(i, surfaceT, flux, 0, DAY);
  assert.equal(land.snow[i], 0);
  assert.equal(land.snowAlbedo[i], 0.85, 'bare ground holds the fresh albedo for the next snow');
});

test('trees standing above the snow darken it linearly to forestSnowAlbedo at closedCanopy, and without the treeline the standing cover outlasts the cover\'s decay under snow', () => {
  const land = createLandSurface(mesh, flat(), { treeline: false, grassland: false }), open = createLandSurface(mesh, flat(), { snowMasking: false, grassland: false });
  land.initialize(); open.initialize();
  const i = 4;
  for (const m of [land, open]) { m.snow[i] = 40; m.snowAlbedo[i] = 0.8; }
  for (const [standing, expected] of [[0, 0.8], [0.35, 0.8 + (0.27 - 0.8) * 0.5], [0.7, 0.27], [1, 0.27]]) {
    land.canopy[i] = standing; open.canopy[i] = standing;
    assert.ok(near(land.albedo(i), expected), `standing ${standing}: ${land.albedo(i)}`);
    assert.ok(near(open.albedo(i), 0.8), 'without masking the snow is its own');
  }
  land.snow[i] = 0; land.canopy[i] = 0.35;
  const ground = land.albedo(i);
  land.snow[i] = 10;
  assert.ok(near(land.albedo(i), ground + 0.5 * ((0.8 + (0.27 - 0.8) * 0.5) - ground)), 'thin snow blends with the ground as before');
  const surfaceT = new Float64Array(mesh.nCells).fill(MELTING_POINT - 10), flux = new Float64Array(mesh.nCells);
  land.snow[i] = 100; land.vegetation[i] = 0.8; land.canopy[i] = 0.8;
  land.update(i, surfaceT, flux, 0, 100 * DAY);
  assert.ok(near(land.vegetation[i], 0.6962597806667126), `cover ${land.vegetation[i]}`);
  assert.ok(near(land.canopy[i], 0.7751389579641572), `standing ${land.canopy[i]}`);
  land.snow[i] = 0; land.soil[i] = land.capacity(i);
  const warm = new Float64Array(mesh.nCells).fill(295);
  land.update(i, warm, flux, 0, 300 * DAY);
  assert.equal(land.canopy[i], land.vegetation[i], 'regrowth past the standing cover carries it');
});

test('bare sea ice is 0.62 cold and darkens to 0.48 over the last kelvin below melting; the snow on it ages toward 0.70 and is refreshed by snowfall', () => {
  const sea = createSeaIce(mesh);
  for (const [T, bare] of [[MELTING_POINT - 10, 0.62], [MELTING_POINT - 1, 0.62], [MELTING_POINT - 0.5, 0.55], [MELTING_POINT, 0.48]]) {
    assert.ok(near(sea.albedo(1, null, 0, 1, T), bare), `T ${T}: ${sea.albedo(1, null, 0, 1, T)}`);
    assert.ok(near(sea.albedo(0.25, null, 0, 1, T), 0.06 + (bare - 0.06) * 0.5), 'thin ice keeps its ramp from the water');
    assert.ok(near(sea.albedo(1, null, 10, 1, T, 0.8), bare + 0.5 * (0.8 - bare)), 'the snow at its own albedo over 20 kg/m2');
  }
  const fixed = createSeaIce(mesh, { snowAgeing: false });
  assert.ok(near(fixed.albedo(1, null, 20, 1, MELTING_POINT, 0.6), 0.75), 'without ageing the snow is iceSnowAlbedo');
  const surfaceT = new Float64Array([MELTING_POINT - 20]), ice = new Float64Array([1]), flux = new Float64Array([0]);
  sea.deposit(0, 30, MELTING_POINT - 20, ice, surfaceT);
  assert.equal(sea.snowAlbedo[0], 0.85);
  sea.snowAlbedo[0] = 0.75;
  sea.deposit(0, 2, MELTING_POINT - 20, ice, surfaceT);
  assert.ok(near(sea.snowAlbedo[0], 0.77));
  surfaceT[0] = MELTING_POINT - 20;
  sea.update(surfaceT, ice, flux, 0, 3600);
  assert.ok(near(sea.snowAlbedo[0], agedSnowAlbedo(0.77, surfaceT[0], 3600, 0.7)), 'aged at the skin it ends the step with');
  surfaceT[0] = MELTING_POINT; flux[0] = 50;
  for (let n = 0; n < 24; n++) sea.update(surfaceT, ice, flux, 0, 3600);
  assert.ok(sea.snowAlbedo[0] < 0.76 && sea.snowAlbedo[0] > 0.7, `melting snow ${sea.snowAlbedo[0]}`);
  flux[0] = 2e4;
  for (let n = 0; n < 24 && sea.snow[0] > 0; n++) sea.update(surfaceT, ice, flux, 0, 3600);
  assert.equal(sea.snow[0], 0);
  assert.equal(sea.snowAlbedo[0], 0.85, 'bare ice holds the fresh albedo for the next snow');
  const water = new Float64Array([FREEZING_POINT + 3]), none = new Float64Array([0]);
  sea.snowAlbedo[0] = 0.6;
  sea.update(water, none, new Float64Array([0]), 0, 900);
  assert.equal(sea.snowAlbedo[0], 0.85, 'open water too');
});

test('the land shares its snow albedo with the sea ice, saves it and the standing cover, and an older state loads them fresh and at the cover (without the treeline)', async () => {
  const geography = createGeography(mesh, syntheticTopography(180, 360, (lat, lon) => (Math.cos(lon) > 0 ? 500 : -4000)), { landBridges: {}, seaStraits: {} });
  const model = createModel(new Grid(8), { physics: true, topography: syntheticTopography(180, 360, (lat, lon) => (Math.cos(lon) > 0 ? 500 : -4000)), geography: { landBridges: {}, seaStraits: {} } });
  assert.equal(model.seaIce.snowAlbedo.buffer, model.land.snowAlbedo.buffer);
  const land = createLandSurface(mesh, geography, { treeline: false }), C = mesh.nCells;
  land.initialize();
  const snow = Float64Array.from({ length: C }, (_, i) => (i % 3 === 0 ? 30 : 0)), vegetation = Float64Array.from({ length: C }, (_, i) => (i % 5) / 5);
  const ice = Float64Array.from({ length: C }, (_, i) => (geography.land[i] ? 0 : 1));
  land.load({ soil: new Float64Array(C).fill(50), snow, vegetation }, ice);
  for (let i = 0; i < C; i++) {
    assert.equal(land.snowAlbedo[i], 0.85);
    assert.equal(land.canopy[i], land.vegetation[i]);
  }
  const snowAlbedo = Float64Array.from({ length: C }, (_, i) => 0.5 + 0.3 * ((i * 7) % 11) / 10), canopy = Float64Array.from({ length: C }, (_, i) => Math.min(1, vegetation[i] + 0.25));
  land.load({ soil: new Float64Array(C).fill(50), snow, vegetation, snowAlbedo, canopy }, ice);
  const saved = await decodeState(encodeState({ N: 8, K: 1, day: 0, time: 0, terrain: true, land: land.serialize() }));
  const back = createLandSurface(mesh, geography, { treeline: false });
  back.load(saved.land, ice);
  for (let i = 0; i < C; i++) {
    const kept = snow[i] > 0 && (geography.land[i] || ice[i] > 0);
    assert.ok(near(back.snowAlbedo[i], kept ? Math.fround(snowAlbedo[i]) : 0.85, 1e-7), `cell ${i}: snow albedo ${back.snowAlbedo[i]}`);
    const standing = geography.land[i] && !geography.iceSheet[i] && vegetation[i] > 0 ? Math.max(Math.fround(vegetation[i]), Math.fround(canopy[i])) : back.vegetation[i];
    assert.ok(near(back.canopy[i], standing, 1e-7), `cell ${i}: standing ${back.canopy[i]} vs ${standing}`);
  }
});

const topography = syntheticTopography(90, 180, (lat, lon) => ((Math.cos(lon) > 0 && Math.abs(lat) < 1.2) || lat < -1.15 ? 300 : -4000));

test('both engines give land, ice and the snow on them the same albedo from random snow, snow albedo, standing cover and skin temperature', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const cpu = createModel(new Grid(6), { topography }), gpu = await createGpuModel(new Grid(6), { topography });
  const C = cpu.mesh.nCells, landMask = cpu.geography.land;
  let s = 11;
  const rnd = () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; };
  const init = initializeState(cpu, {});
  for (let i = 0; i < C; i++) {
    if (landMask[i]) { init[6][i] = 0; init[3][i] = MELTING_POINT - 15 + 16 * rnd(); continue; }
    if (rnd() < 0.5) { init[6][i] = 0.1 + 2 * rnd(); init[3][i] = MELTING_POINT - 4 + 4 * rnd(); }
  }
  const soil = Float64Array.from({ length: C }, (_, i) => (landMask[i] ? 300 * rnd() : 0));
  const snow = Float64Array.from({ length: C }, () => (rnd() < 0.6 ? 40 * rnd() : 0));
  const vegetation = Float64Array.from({ length: C }, () => rnd()), canopy = Float64Array.from(vegetation, (v) => Math.min(1, v + 0.3 * rnd()));
  const snowAlbedo = Float64Array.from({ length: C }, () => 0.5 + 0.35 * rnd());
  for (const m of [cpu, gpu]) {
    for (let a = 0; a < init.length; a++) m.state[a].set(init[a]);
    m.seaIce.load(m.state[6], null);
    m.land.load({ soil, snow, vegetation, canopy, snowAlbedo }, m.state[6]);
  }
  gpu.load();
  const expected = Float64Array.from({ length: C }, (_, i) => (landMask[i] ? cpu.land.albedo(i) : cpu.seaIce.albedo(cpu.state[6][i], null, cpu.seaIce.snow[i], cpu.seaIce.cover(i, cpu.state[6][i]), cpu.state[3][i], cpu.seaIce.snowAlbedo[i])));
  await gpu.step(900);
  const ph = await gpu.gpu.downloadPhysics();
  let landWorst = 0, iceWorst = 0, iced = 0, masked = 0, melting = 0;
  for (let i = 0; i < C; i++) {
    const d = Math.abs(ph.ADIF[i] - expected[i]);
    if (landMask[i]) { landWorst = Math.max(landWorst, d); if (snow[i] > 0 && canopy[i] > 0.1 && !(cpu.geography.iceSheet && cpu.geography.iceSheet[i])) masked++; }
    else if (init[6][i] > 0) { iceWorst = Math.max(iceWorst, d); iced++; if (init[3][i] > MELTING_POINT - 1) melting++; }
  }
  gpu.destroy();
  console.log(`N=6: land albedo engines apart by ${landWorst.toExponential(1)} (${masked} snow cells under a standing cover), ice by ${iceWorst.toExponential(1)} over ${iced} iced cells (${melting} within 1 K of melting)`);
  assert.ok(masked > 20 && iced > 20 && melting > 5, `${masked} masked, ${iced} iced, ${melting} melting`);
  assert.ok(landWorst < 1e-6 && iceWorst < 1e-5, `engines apart by ${landWorst} over land and ${iceWorst} over ice`);
});

test('over 48 GPU steps the snow albedo and the standing cover evolve as on the CPU, on land and on the ice', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const options = { land: { treeline: false, coldSnowAgeing: 2, meltingSnowAgeing: 12, canopyMemory: 6 * 3600, snowDeclineTime: 4 * 3600, refreshSnowfall: 0.5, wetSnowRange: 10 }, ice: { coldSnowAgeing: 2, meltingSnowAgeing: 12, refreshSnowfall: 0.5, wetSnowRange: 10 } };
  const prepare = (model) => {
    const C = model.mesh.nCells, init = initializeState(model, {});
    for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
    for (let i = 0; i < C; i++) if (model.geography.land[i]) model.state[6][i] = 0;
    for (let x = 0; x < model.state[4].length; x++) if (Math.abs(model.mesh.latCell[x % C]) > 0.9) model.state[4][x] *= 3;
    model.seaIce.load(model.state[6], null);
    if (model.load) model.load();
    model.ocean.initialize(model.state[3], model.state[6]);
    const snow = Float64Array.from({ length: C }, (_, i) => (Math.abs(model.mesh.latCell[i]) > 0.9 || Math.abs(model.state[3][i] - MELTING_POINT) < 4 ? 30 : 0));
    model.land.load({ soil: new Float64Array(C).fill(60), snow, vegetation: new Float64Array(C).fill(0.6), snowAlbedo: new Float64Array(C).fill(0.8) }, model.state[6]);
    return model;
  };
  const cpu = prepare(createModel(new Grid(6), { topography, ...options }));
  const gpu = prepare(await createGpuModel(new Grid(6), { topography, ...options }));
  for (let n = 0; n < 48; n++) { cpu.step(900); await gpu.step(900); }
  await gpu.sync();
  const saved = await gpu.land.serialize(), C = cpu.mesh.nCells;
  let worst = 0, sum = 0, n = 0, moved = 0, refreshed = 0, wet = 0, canopyWorst = 0, standing = 0, at = -1;
  for (let i = 0; i < C; i++) {
    if (!(cpu.land.snow[i] > 0)) continue;
    const d = Math.abs(cpu.land.snowAlbedo[i] - saved.snowAlbedo[i]);
    if (d > worst) { worst = d; at = i; }
    sum += d * d; n++;
    if (Math.abs(cpu.land.snowAlbedo[i] - 0.8) > 0.01) moved++;
    if (cpu.land.snow[i] > 30.01) refreshed++;
    if (cpu.state[3][i] >= MELTING_POINT - 10) wet++;
    if (cpu.geography.land[i]) { canopyWorst = Math.max(canopyWorst, Math.abs(cpu.land.canopy[i] - saved.canopy[i])); if (cpu.land.canopy[i] > cpu.land.vegetation[i] + 0.05) standing++; }
  }
  const rms = Math.sqrt(sum / n);
  gpu.destroy();
  console.log(`N=6, 48 steps: ${n} snow cells, ${moved} moved from 0.8 (${refreshed} snowed on, ${wet} wet), snow albedo engines rms ${rms.toExponential(1)}, max ${worst.toExponential(1)} at ${at}; standing cover above the cover on ${standing} cells, engines max ${canopyWorst.toExponential(1)}`);
  assert.ok(moved > n / 2 && refreshed > 3 && standing > 5 && wet > 5 && n - wet > 5, `${moved} moved, ${refreshed} snowed on, ${wet} of ${n} wet, ${standing} standing`);
  assert.ok(rms < 1e-4 && worst < 1e-3, `snow albedo rms ${rms}, max ${worst}`);
  assert.ok(canopyWorst < 1e-4, `standing cover max ${canopyWorst}`);
});

test('the GPU, whose land and sea ice share one snow ageing, refuses an ageing option given to one of them alone or differently to each', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  for (const options of [{ land: { snowAgeing: false } }, { ice: { snowAgeing: false } }, { ice: { coldSnowAgeing: 0.01 } }, { land: { refreshSnowfall: 5 }, ice: { refreshSnowfall: 8 } }]) {
    await assert.rejects(createGpuModel(new Grid(4), { topography, ...options }), /share/, JSON.stringify(options));
  }
  const both = await createGpuModel(new Grid(4), { topography, land: { coldSnowAgeing: 0.01 }, ice: { coldSnowAgeing: 0.01 } });
  assert.equal(both.gpu.physics.coldSnowAgeing, 0.01);
  both.destroy();
  const seaOnly = await createGpuModel(new Grid(4), { ice: { snowAgeing: false } });
  assert.equal(seaOnly.gpu.physics.snowAgeing, false, 'without land the sea ice sets it');
  seaOnly.destroy();
});
