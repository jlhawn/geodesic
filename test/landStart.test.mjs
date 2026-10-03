import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createGeography, syntheticTopography } from '../js/geography.module.js';
import { createLandSurface, treelineFactor, aridityFactor, litterWarmth, decompositionWarmth, decompositionMoisture, carbonEquilibrium, carbonRecord, recordWeight, freshJumpDue, startPlaceholders, dryHumusAlbedo, START_CODES, FRESH_JUMPS } from '../js/physics/land.module.js';
import { MELTING_POINT } from '../js/physics/ice.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { encodeState, decodeState } from '../js/stateFile.module.js';
import { regridLand } from '../js/physics/regrid.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const DAY = 86400, YEAR = 365 * DAY;
const mesh = createModel(new Grid(8), { physics: false }).mesh;
const flat = () => createGeography(mesh, syntheticTopography(90, 180, () => 100), { landBridges: {}, seaStraits: {} });
const near = (a, b, tol = 1e-12) => Math.abs(a - b) <= tol * Math.max(1, Math.abs(b));
const topography = syntheticTopography(90, 180, (lat, lon) => ((Math.cos(lon) > 0 && Math.abs(lat) < 1.2) || lat < -1.15 ? 300 : -4000));

test('the start options set the cover, the held trees and the held topsoil carbon, and empty the record at age 0', () => {
  const p = startPlaceholders();
  assert.deepEqual(p, { neutral: { cover: 0.5, share: 0.5, carbon: 2.6 * Math.LN2 }, bare: { cover: 0, share: 0, carbon: 0 }, green: { cover: 1, share: null, carbon: 13 } });
  assert.ok(near(dryHumusAlbedo(p.neutral.carbon), (0.37 + 0.12) / 2), 'the neutral soil is midway between mineral and humus');
  assert.ok(Math.abs(dryHumusAlbedo(p.green.carbon) - 0.12) < 0.002, 'the green soil is within e⁻⁵ of humus');
  const geography = createGeography(mesh, syntheticTopography(180, 360, (lat, lon) => (lat < -1.2 ? 2000 : Math.cos(lon) > 0 ? 100 : -4000)), { landBridges: {}, seaStraits: {} });
  for (const start of ['neutral', 'bare', 'green']) {
    const land = createLandSurface(mesh, geography, { start });
    land.initialize();
    assert.deepEqual([...land.record], [0, START_CODES[start], START_CODES[start]]);
    let cells = 0;
    for (let i = 0; i < mesh.nCells; i++) {
      for (const field of ['seasonLength', 'seasonWarmth', 'rainMean', 'demandMean', 'litterMean', 'decayMean']) assert.equal(land[field][i], 0, `${start} ${field} ${i}`);
      if (!geography.land[i] || geography.iceSheet[i]) { assert.equal(land.vegetation[i], 0); assert.equal(land.canopy[i], 0); assert.equal(land.soilCarbon[i], 0); continue; }
      cells++;
      assert.equal(land.vegetation[i], p[start].cover);
      assert.equal(land.snowFreeCover[i], p[start].cover);
      assert.equal(land.canopy[i], start === 'neutral' ? 0.25 : 0, `${start}: the green trees wait for a record`);
      assert.equal(land.soilCarbon[i], p[start].carbon);
      assert.equal(land.soil[i], 150);
    }
    assert.ok(cells > 50);
  }
  assert.throws(() => createLandSurface(mesh, geography, { start: 'atlas' }), /start/);
  const old = createLandSurface(mesh, geography);
  old.load({ soil: new Float64Array(mesh.nCells).fill(100), snow: new Float64Array(mesh.nCells) });
  assert.deepEqual([...old.record], [-1, 0, 0], 'a state without a record keeps exponential means and holds nothing');
});

test('a fresh record is the plain average of every step while it is within its memory and the exponential mean after, for the season, moisture and carbon means alike', () => {
  const dt = 900, land = createLandSurface(mesh, flat(), { seasonMemory: 6 * dt, moistureMemory: 4 * dt });
  land.initialize();
  const i = 5, j = 9, cap = land.capacity(i), flux = new Float64Array(mesh.nCells), surfaceT = new Float64Array(mesh.nCells).fill(280);
  land.soil[i] = 0.4 * cap;
  const air = (k) => MELTING_POINT - 3 + 1.7 * k, potential = (k) => (1 + 0.3 * k) / DAY, rain = (k) => (k % 3 === 0 ? 0 : 0.5 + 0.1 * k);
  const inputs = { seasonLength: [], seasonWarmth: [], demandMean: [], litterMean: [], decayMean: [], rainMean: [] };
  const hand = Object.fromEntries(Object.keys(inputs).map((k) => [k, 0]));
  for (let k = 0; k < 10; k++) {
    const T = air(k), fill = land.soil[i] / cap;
    inputs.seasonLength.push(T >= MELTING_POINT + 0.9 ? 1 : 0);
    inputs.seasonWarmth.push(Math.max(0, T - MELTING_POINT - 0.9));
    inputs.demandMean.push(86400 * potential(k));
    inputs.litterMean.push(Math.min(1, fill / 0.75) * litterWarmth(T));
    inputs.decayMean.push(decompositionWarmth(T) * decompositionMoisture(fill, T < MELTING_POINT));
    inputs.rainMean.push(86400 * rain(k) / dt);
    land.update(i, surfaceT, flux, 0, dt, T, potential(k));
    land.deposit(j, rain(k), 290, dt);
    land.advance(dt);
    for (const [name, values] of Object.entries(inputs)) {
      const memory = name === 'demandMean' || name === 'rainMean' ? 4 * dt : 6 * dt, n = values.length;
      if (n * dt <= memory) hand[name] = values.reduce((s, x) => s + x, 0) / n;
      else hand[name] += (values[n - 1] - hand[name]) * (1 - Math.exp(-dt / memory));
      const cell = name === 'rainMean' ? j : i;
      assert.ok(near(land[name][cell], hand[name], 1e-14), `${name} after ${n} steps: ${land[name][cell]} against ${hand[name]}`);
    }
    assert.equal(land.record[0], (k + 1) * dt);
    assert.equal(land.canopy[i], 0.5 * land.vegetation[i], 'the trees are held at half the cover');
    assert.equal(land.soilCarbon[i], 2.6 * Math.LN2, 'the carbon is held');
  }
  assert.ok(inputs.seasonLength.includes(0) && inputs.seasonLength.includes(1) && inputs.decayMean.some((x, k) => x > 0 && air(k) < MELTING_POINT), 'the steps cross the season threshold and freezing');
  assert.equal(recordWeight(0, dt, 6 * dt), 1);
  assert.equal(recordWeight(2 * dt, dt, 6 * dt), 1 / 3);
  assert.equal(recordWeight(5 * dt, dt, 6 * dt), 1 / 6);
  assert.equal(recordWeight(6 * dt, dt, 6 * dt), 1 - Math.exp(-1 / 6));
  assert.equal(recordWeight(-1, dt, 6 * dt), 1 - Math.exp(-1 / 6));
});

test('in float32 a plain average over three years of steps stays within 3·10⁻⁴ of the exact mean', () => {
  const f = Math.fround;
  let s = 7;
  const rnd = () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; };
  const air = (t, mean, amplitude) => MELTING_POINT + mean + amplitude * Math.cos(2 * Math.PI * t / YEAR) + 6 * Math.cos(2 * Math.PI * t / DAY) + 2 * (rnd() - 0.5);
  const signals = {
    'season length, boreal': (t) => (air(t, -3, 18) >= MELTING_POINT + 0.9 ? 1 : 0),
    'season warmth, boreal': (t) => Math.max(0, air(t, -3, 18) - MELTING_POINT - 0.9),
    'rain in showers': () => (rnd() < 0.03 ? 100 * rnd() : 0),
    'decomposition, tundra': (t) => { const T = air(t, -11, 16); return decompositionWarmth(T) * (T < MELTING_POINT ? 0.2 : 0.64); },
  };
  const lines = [];
  for (const dt of [337.5, 168.75]) {
    for (const [name, signal] of Object.entries(signals)) {
      s = 7;
      let m = 0, sum = 0, carry = 0, age = 0, worst = 0, scale = 0;
      const steps = Math.round(3 * YEAR / dt);
      for (let k = 0; k < steps; k++) {
        const x = f(signal(k * dt)), w = f(recordWeight(age, dt, 3 * YEAR));
        m = f(m + f(f(x - m) * w));
        const y = x - carry, t = sum + y;
        carry = (t - sum) - y; sum = t; age += dt;
        const exact = sum / (k + 1);
        scale = Math.max(scale, Math.abs(exact));
        worst = Math.max(worst, Math.abs(m - exact));
      }
      lines.push(`dt ${dt} s ${name}: ${steps} steps, worst ${(worst / scale).toExponential(1)} of the mean's scale`);
      assert.ok(worst < 3e-4 * scale, lines[lines.length - 1]);
    }
  }
  console.log(lines.join('\n'));
});

test('the jump sets the trees to the cover times the record\'s factors and the carbon to its equilibrium, releases the hold, and changes nothing when repeated', () => {
  const geography = createGeography(mesh, syntheticTopography(180, 360, (lat, lon) => (lat < -1.2 ? 2000 : Math.cos(lon) > 0 ? 100 : -4000)), { landBridges: {}, seaStraits: {} });
  const land = createLandSurface(mesh, geography), C = mesh.nCells;
  land.initialize();
  const cells = [...Array(C).keys()].filter((i) => geography.land[i] && !geography.iceSheet[i]), sheet = [...Array(C).keys()].find((i) => geography.iceSheet[i]);
  const [a, b, c] = cells;
  land.seasonLength[a] = 0.5; land.seasonWarmth[a] = 3.25; land.rainMean[a] = 2; land.demandMean[a] = 4; land.vegetation[a] = 0.8; land.litterMean[a] = 0.3; land.decayMean[a] = 0.6;
  land.seasonLength[b] = 0.2; land.seasonWarmth[b] = 0.5; land.rainMean[b] = 5; land.demandMean[b] = 2; land.vegetation[b] = 0.4; land.litterMean[b] = 0.1; land.decayMean[b] = 0.05;
  for (const field of ['seasonLength', 'seasonWarmth', 'rainMean', 'demandMean', 'litterMean', 'decayMean']) land[field][c] = land[field][a];
  land.snow[c] = 30; land.vegetation[c] = 0.3; land.snowFreeCover[c] = 0.6;
  const report = land.jump();
  assert.ok(near(land.canopy[a], 0.625 * 0.375 * 0.8), `a 183-day season at 7.4 °C and P/PET 0.5: trees ${land.canopy[a]}`);
  assert.ok(near(land.soilCarbon[a], 0.5 * 70 * 0.8 * 0.3 / 0.6), `S* = 35 kg/m² × cover × 0.3/0.6: ${land.soilCarbon[a]}`);
  assert.equal(land.canopy[b], 0, 'a season too short for trees');
  assert.ok(near(land.canopy[c], 0.625 * 0.375 * 0.6), `under snow the trees take the cover of the last snow-free step: ${land.canopy[c]}`);
  assert.ok(near(land.soilCarbon[c], 35 * 0.6 * 0.3 / 0.6), `and so does the litter: ${land.soilCarbon[c]}`);
  assert.ok(near(land.soilCarbon[b], 35 * 0.4 * 2), `all grass: ${land.soilCarbon[b]}`);
  assert.equal(land.soilCarbon[sheet], 0); assert.equal(land.canopy[sheet], 0);
  assert.equal(land.record[1], 0, 'the hold is released');
  for (const i of cells.slice(3, 40)) {
    assert.ok(near(land.canopy[i], treelineFactor(land.seasonLength[i], land.seasonWarmth[i]) * aridityFactor(land.rainMean[i], land.demandMean[i]) * land.vegetation[i]));
    assert.equal(land.soilCarbon[i], 0, 'an empty record holds no carbon');
  }
  assert.ok(report.global.trees.length === 2 && report.bands.length > 5 && near(report.bands.reduce((s, x) => s + x.share, 0), 1), JSON.stringify(report.global));
  const again = { canopy: Float64Array.from(land.canopy), soilCarbon: Float64Array.from(land.soilCarbon) };
  const second = land.jump();
  assert.deepEqual([...land.canopy], [...again.canopy]);
  assert.deepEqual([...land.soilCarbon], [...again.soilCarbon]);
  assert.deepEqual(second.global.carbon, [second.global.carbon[1], second.global.carbon[1]], 'the repeat reports no change');
  const fill = 0.55, sine = carbonRecord(10, 11, fill);
  for (const i of cells) { land.soil[i] = fill * 300; land.vegetation[i] = 0.9; land.snow[i] = 0; land.litterMean[i] = sine.litter; land.decayMean[i] = sine.decay; }
  land.jump();
  for (const i of cells) assert.ok(near(land.soilCarbon[i], carbonEquilibrium(10, 11, fill, Math.min(1, land.canopy[i]) + Math.max(0, 0.9 - land.canopy[i])), 1e-13), `cell ${i}: the jump from a sine year's record lands on carbonEquilibrium`);
});

test('the snow-free cover follows the cover while the cell is free of snow and keeps it under snow', () => {
  const land = createLandSurface(mesh, flat());
  land.initialize();
  const i = [...Array(mesh.nCells).keys()].find((n) => land.land[n]), C = mesh.nCells, flux = new Float64Array(C), surfaceT = new Float64Array(C).fill(290);
  land.soil[i] = 0.9 * land.capacity(i);
  for (let k = 0; k < 4; k++) { land.update(i, surfaceT, flux, 0, 86400, MELTING_POINT + 15, 2 / DAY); land.advance(86400); assert.equal(land.snowFreeCover[i], land.vegetation[i]); }
  const before = land.vegetation[i];
  assert.ok(before > 0.5, `the cover grew: ${before}`);
  land.deposit(i, 40, MELTING_POINT - 5, 86400);
  surfaceT[i] = MELTING_POINT - 5;
  for (let k = 0; k < 20; k++) { land.update(i, surfaceT, flux, 0, 86400, MELTING_POINT - 5, 0); land.advance(86400); surfaceT[i] = MELTING_POINT - 5; }
  assert.ok(land.snow[i] > 0 && land.vegetation[i] < before * 0.98, `the cover decays under snow: ${land.vegetation[i]} from ${before}`);
  assert.equal(land.snowFreeCover[i], before, 'the snow-free cover keeps the cover the snow came on');
  const saved = land.serialize(), back = createLandSurface(mesh, flat());
  back.load(saved);
  assert.equal(back.snowFreeCover[i], before);
  const old = createLandSurface(mesh, flat()), { snowFreeCover, ...older } = saved;
  old.load(older);
  assert.equal(old.snowFreeCover[i], old.vegetation[i], 'a state without one takes its cover');
});

test('a state saved after the jump reloads its record exactly and its fields as float32', async () => {
  const land = createLandSurface(mesh, flat(), { start: 'green' });
  land.initialize();
  const C = mesh.nCells, flux = new Float64Array(C), surfaceT = new Float64Array(C).fill(285);
  for (let k = 0; k < 5; k++) { for (let i = 0; i < C; i++) if (land.land[i]) { land.update(i, surfaceT, flux, 0, 3600, MELTING_POINT + 4 + (i % 9), (i % 5) / DAY); land.deposit(i, (i % 4) * 0.3, 285, 3600); } land.advance(3600); }
  land.jump();
  const direct = createLandSurface(mesh, flat());
  direct.load(land.serialize());
  for (const name of ['canopy', 'soilCarbon', 'seasonLength', 'seasonWarmth', 'rainMean', 'demandMean', 'litterMean', 'decayMean', 'vegetation', 'snowFreeCover', 'soil']) assert.deepEqual([...direct[name]], [...land[name]], name);
  assert.deepEqual([...direct.record], [5 * 3600, 0, 3]);
  const file = await decodeState(encodeState({ N: 8, K: 1, day: 0, time: 0, terrain: true, land: land.serialize() }));
  assert.ok(file.land.record instanceof Float64Array);
  const back = createLandSurface(mesh, flat());
  back.load(file.land);
  assert.deepEqual([...back.record], [...land.record]);
  for (const name of ['canopy', 'soilCarbon', 'litterMean', 'decayMean']) for (let i = 0; i < C; i++) assert.equal(back[name][i], Math.fround(land[name][i]), `${name} ${i}`);
});

test('a fresh record jumps at 365 and 730 days of its age, and a state without one never', () => {
  assert.deepEqual(FRESH_JUMPS, [365, 730]);
  const fresh = [364 * DAY, 1, 1];
  assert.equal(freshJumpDue(fresh, 364 * DAY, 365 * DAY), true);
  assert.equal(freshJumpDue(fresh, 365 * DAY, 366 * DAY), false, 'a segment that starts at the jump does not repeat it');
  assert.equal(freshJumpDue(fresh, 729.5 * DAY, 730 * DAY), true);
  assert.equal(freshJumpDue(fresh, 731 * DAY, 800 * DAY), false);
  assert.equal(freshJumpDue([-1, 0, 0], -1, -1), false);
  assert.equal(freshJumpDue([400 * DAY, 0, 0], 0, 800 * DAY), false, 'a land that did not start fresh');
});

test('land regridding carries the record and the carbon\'s record and fills new land from the sine year', () => {
  const make = (N, relief) => { const m = createModel(new Grid(N), { physics: false, topography: relief }); return { mesh: m.mesh, geography: createGeography(m.mesh, relief, { landBridges: {}, seaStraits: {} }) }; };
  const islands = syntheticTopography(90, 180, (lat, lon) => (lat > 1.3 && Math.cos(lon) < 0 ? 50 : 0) || (((Math.cos(lon) > 0 && Math.abs(lat) < 1.2) || lat < -1.15) ? 300 : -4000));
  const source = make(6, topography), target = make(8, islands), C = source.mesh.nCells;
  const land = { soil: new Float64Array(C).fill(150), snow: new Float64Array(C), vegetation: new Float64Array(C).fill(0.5), snowFreeCover: new Float64Array(C).fill(0.4), soilCarbon: new Float64Array(C).fill(3), litterMean: new Float64Array(C).fill(0.25), decayMean: new Float64Array(C).fill(0.75), record: Float64Array.of(1e6, 1, 1) };
  const out = regridLand(source, target, land);
  assert.deepEqual([...out.record], [1e6, 1, 1]);
  let carried = 0, filled = 0;
  for (let n = 0; n < target.mesh.nCells; n++) {
    if (!target.geography.land[n]) { assert.equal(out.litterMean[n], 0); assert.equal(out.decayMean[n], 0); assert.equal(out.snowFreeCover[n], 0); continue; }
    assert.ok(out.snowFreeCover[n] === 0.4 || out.snowFreeCover[n] === 0.5 || out.snowFreeCover[n] === 0, `cell ${n}: snow-free cover ${out.snowFreeCover[n]}`);
    if (out.litterMean[n] === 0.25 && out.decayMean[n] === 0.75) { carried++; continue; }
    filled++;
    const e = carbonRecord(0, 0, 0.5);
    assert.ok(out.decayMean[n] > 0 && out.litterMean[n] >= 0 && e.decay > 0, `cell ${n}: ${out.litterMean[n]}, ${out.decayMean[n]}`);
  }
  assert.ok(carried > 50 && filled > 3, `${carried} carried, ${filled} filled`);
  assert.deepEqual([...regridLand(source, source, land).record], [1e6, 1, 1]);
});

const RAINING_START = { capeClosure: 'threshold', plumeSourceDepth: 'boundaryLayer', cumulusClosure: 0.06 };

test('a CPU model\'s fresh land holds its trees and carbon through its steps, and after the jump they run free', () => {
  const model = createModel(new Grid(6), { topography, land: { start: 'neutral', treeGrowthTime: 3600, treeDeclineTime: 3600 }, moist: RAINING_START });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i]) model.state[6][i] = 0;
  model.land.initialize();
  for (let n = 0; n < 6; n++) model.step(900);
  const { land } = model, cells = [...Array(model.mesh.nCells).keys()].filter((i) => model.geography.land[i] && !model.geography.iceSheet[i]);
  assert.equal(land.record[0], 6 * 900);
  for (const i of cells) { assert.equal(land.canopy[i], 0.5 * land.vegetation[i]); assert.equal(land.soilCarbon[i], 2.6 * Math.LN2); }
  land.jump();
  const jumped = Float64Array.from(land.canopy);
  model.step(900);
  const moved = cells.filter((i) => Math.abs(land.canopy[i] - jumped[i]) > 1e-6).length, off = cells.filter((i) => Math.abs(land.canopy[i] - 0.5 * land.vegetation[i]) > 1e-3).length;
  assert.ok(moved > 10 && off > 10, `${moved} cells' trees moved after the jump, ${off} away from half the cover`);
});

const prepareGpu = async (land, start = 'neutral') => {
  const model = await createGpuModel(new Grid(6), { topography, land: { start, ...land }, moist: RAINING_START });
  const C = model.mesh.nCells, init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  for (let i = 0; i < C; i++) if (model.geography.land[i]) model.state[6][i] = 0;
  model.seaIce.load(model.state[6], null);
  model.load();
  model.ocean.initialize(model.state[3], model.state[6]);
  model.land.initialize();
  return model;
};
const RECORDS = ['seasonLength', 'seasonWarmth', 'rainMean', 'demandMean', 'litterMean', 'decayMean'];

test('on the GPU the fresh record is the plain average of each step\'s own values within its memory and the exponential mean after', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const dt = 900, memory = 6 * dt, options = { seasonMemory: memory, moistureMemory: memory };
  const averaged = await prepareGpu(options), single = await prepareGpu(options);
  const C = averaged.mesh.nCells, cells = [...Array(C).keys()].filter((i) => averaged.geography.land[i] && !averaged.geography.iceSheet[i]);
  const hand = Object.fromEntries(RECORDS.map((name) => [name, new Float64Array(C)])), sums = Object.fromEntries(RECORDS.map((name) => [name, new Float64Array(C)]));
  const worst = Object.fromEntries(RECORDS.map((name) => [name, 0])), scale = Object.fromEntries(RECORDS.map((name) => [name, 0]));
  for (let n = 1; n <= 10; n++) {
    single.land.record[0] = 0;
    await averaged.step(dt); await single.step(dt);
    const [mean, own] = [await averaged.land.serialize(), await single.land.serialize()];
    for (const i of cells) assert.equal(mean.vegetation[i], own.vegetation[i], 'the record feeds nothing back while the land is held');
    for (const name of RECORDS) {
      for (const i of cells) {
        sums[name][i] += own[name][i];
        hand[name][i] = n * dt <= memory ? sums[name][i] / n : hand[name][i] + (own[name][i] - hand[name][i]) * Math.fround(1 - Math.exp(-dt / memory));
        worst[name] = Math.max(worst[name], Math.abs(mean[name][i] - hand[name][i]));
        scale[name] = Math.max(scale[name], Math.abs(hand[name][i]));
      }
    }
    if (n === 10) for (const i of cells) { assert.equal(mean.canopy[i], 0.5 * mean.vegetation[i]); assert.equal(mean.soilCarbon[i], Math.fround(2.6 * Math.LN2)); }
  }
  averaged.destroy(); single.destroy();
  console.log(`N=6, 10 steps (6 within the memory) over ${cells.length} land cells, worst against the mean of each step's own values: ${RECORDS.map((name) => `${name} ${worst[name].toExponential(1)} of ${scale[name].toFixed(2)}`).join(', ')}`);
  for (const name of RECORDS) assert.ok(scale[name] > 0 && worst[name] <= 2e-6 * scale[name], `${name}: ${worst[name]} of ${scale[name]}`);
});

test('the GPU holds the start\'s placeholders, jumps to the hand-computed equilibria of its own record, repeats the jump to the bit and reloads a state saved after it to the bit', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const options = { seasonMemory: 6 * 3600, moistureMemory: 6 * 3600, treeGrowthTime: 3600, treeDeclineTime: 3600 };
  const bare = await prepareGpu(options, 'bare');
  for (let n = 0; n < 4; n++) await bare.step(900);
  const bareState = await bare.land.serialize(), C = bare.mesh.nCells;
  const cells = [...Array(C).keys()].filter((i) => bare.geography.land[i] && !bare.geography.iceSheet[i]);
  for (const i of cells) { assert.equal(bareState.canopy[i], 0); assert.equal(bareState.soilCarbon[i], 0); }
  assert.ok(cells.every((i) => bareState.vegetation[i] >= 0) && bare.land.record[1] === 2);
  bare.destroy();
  const green = await prepareGpu(options, 'green');
  let before = await green.land.serialize();
  for (let n = 0; n < 4; n++) {
    await green.step(900);
    const now = await green.land.serialize();
    for (const i of cells) {
      if (now.snowFreeCover[i] !== before.snowFreeCover[i]) assert.equal(now.snowFreeCover[i], now.vegetation[i], `step ${n + 1} cell ${i}: the snow-free cover is the cover when it moves`);
      if (!(before.snow[i] > 0) && !(now.snow[i] > 0)) assert.equal(now.snowFreeCover[i], now.vegetation[i], `step ${n + 1} cell ${i}: free of snow`);
      const expected = treelineFactor(now.seasonLength[i], now.seasonWarmth[i]) * aridityFactor(before.rainMean[i], now.demandMean[i]) * now.vegetation[i];
      assert.ok(Math.abs(now.canopy[i] - expected) < 1e-5, `step ${n + 1} cell ${i}: green trees ${now.canopy[i]} against ${expected}`);
      assert.equal(now.soilCarbon[i], 13);
    }
    before = now;
  }
  assert.ok(cells.some((i) => before.canopy[i] > 0.1) && cells.some((i) => before.vegetation[i] < 1), 'some green trees stand and the cover runs free');
  const buried = cells.filter((i) => before.canopy[i] > 0.05).slice(0, 12);
  for (const i of buried) { before.snow[i] = 200; before.vegetation[i] = 0.7; before.snowFreeCover[i] = 0.95; }
  green.land.load(before, green.state[6]);
  await green.step(900);
  before = await green.land.serialize();
  for (const i of buried) assert.ok(before.snow[i] > 0 && before.snowFreeCover[i] === Math.fround(0.95) && before.vegetation[i] < 0.7, `cell ${i} under snow: cover ${before.vegetation[i]}, snow-free cover ${before.snowFreeCover[i]}`);
  const report = await green.land.jump();
  const jumped = await green.land.serialize();
  let treed = 0, snowy = 0;
  for (const i of cells) {
    const cover = before.snow[i] > 0 ? before.snowFreeCover[i] : before.vegetation[i];
    if (before.snow[i] > 0 && before.snowFreeCover[i] !== before.vegetation[i]) snowy++;
    const trees = Math.fround(Math.min(1, treelineFactor(before.seasonLength[i], before.seasonWarmth[i]) * aridityFactor(before.rainMean[i], before.demandMean[i]) * cover));
    const litter = Math.min(1, trees) + Math.max(0, cover - trees);
    const carbon = before.decayMean[i] > 0 ? Math.fround(0.5 / YEAR * 70 * YEAR * litter * before.litterMean[i] / before.decayMean[i]) : 0;
    assert.equal(jumped.canopy[i], trees, `cell ${i}: trees`);
    assert.ok(near(jumped.soilCarbon[i], carbon, 1e-7), `cell ${i}: carbon ${jumped.soilCarbon[i]} against ${carbon}`);
    if (trees > 0.05) treed++;
  }
  assert.ok(treed > 10 && green.land.record[1] === 0 && report.global.carbon[1] < 13, `${treed} cells with trees; carbon ${report.global.carbon}`);
  console.log(`GPU jump: ${treed} cells with trees, ${snowy} under snow whose snow-free cover differs from their cover`);
  await green.land.jump();
  const repeated = await green.land.serialize();
  for (const name of ['canopy', 'soilCarbon']) assert.deepEqual([...repeated[name]], [...jumped[name]], `${name} repeated`);
  const file = await decodeState(encodeState({ N: 6, K: 1, day: 0, time: 0, terrain: true, land: repeated }));
  const reloaded = await prepareGpu(options, 'neutral');
  reloaded.land.load(file.land);
  const back = await reloaded.land.serialize();
  for (const name of ['canopy', 'soilCarbon', 'vegetation', 'snowFreeCover', 'soil', 'surface', 'snow', ...RECORDS]) assert.deepEqual(cells.map((i) => back[name][i]), cells.map((i) => repeated[name][i]), `${name} reloaded`);
  assert.deepEqual([...back.record], [...repeated.record]);
  await green.step(900); await reloaded.step(900);
  const freed = await green.land.serialize();
  assert.ok(cells.filter((i) => Math.abs(freed.canopy[i] - jumped.canopy[i]) > 1e-6).length > 10, 'after the jump the trees run free');
  green.destroy(); reloaded.destroy();
});
