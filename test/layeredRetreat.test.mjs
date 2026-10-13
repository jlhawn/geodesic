import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { withCadence } from '../js/cadence.module.js';
import { createSeasonalColumn, columnMesh, DAY } from '../explain/oceanColumn.module.js';
import { createOcean, THIN_LAYER_CLASSES, DAILY, transmitted, deposited } from '../js/ocean/layered.module.js';
import { thermalExpansion } from '../js/ocean/seawater.module.js';
import { pwpColumn } from './helpers/pwp.mjs';

/*
 * The model's own ocean on Grid(2) (explain/oceanColumn.module.js) with
 * the currents held at rest, the mixing package over the
 * two lighter classes, a 1 m floor and IB light, from a 50 m mixed layer:
 * a 12 h half-sine of absorbed sunlight with a daily mean of 250 W/m²
 * (local noon at 12 UTC) and a constant dark loss of N/S times that, a
 * different N/S on three cells of the 58.3°N ring, under a constant wind
 * U, for 120 days at 45 min steps. The dawn depth is the mixed layer's
 * at 06 UTC.
 */
const MEAN_SUN = 250, STEP = 2700, DAYS = 120, GRAVITY = 9.81, RHO_CP = 1025 * 3985;
const PACKAGE = {
  ...THIN_LAYER_CLASSES, penetration: 'IB', carryMomentum: true, momentumDepth: 50, refillRate: 1 / 3600, refillFloor: 5, minimumThickness: 1, shallowestMixedDepth: 1,
  shearEntrainment: 0.65, ekmanDepth: 0.4, momentumFloor: 10, stressThickness: 20, retreatStirring: 2.5, langmuir: 2, retreatEkman: 0.7,
};
const MODES = { off: {}, daily: { retreatDaily: true }, light: { retreatDaily: true, retreatLight: 1 } };
const CELLS = [1, 5, 9], RATIOS = [0.6, 0.7, 0.8];
const sunlight = (t) => {
  const hour = (t / 3600) % 24;
  return hour > 6 && hour < 18 ? MEAN_SUN * Math.PI * Math.sin(Math.PI * (hour - 6) / 12) : 0;
};
const loss = (i) => MEAN_SUN * (RATIOS[CELLS.indexOf(i)] ?? 0.7);
const stressOf = (U) => 1.2 * 1.2e-3 * U * U;
const mean = (values) => values.reduce((s, x) => s + x, 0) / values.length;

function column(U, mode, { days = DAYS, start = 0, restart = null } = {}) {
  return createSeasonalColumn({ latitude: 45, heating: 0, damping: 0, wind: stressOf(U), dt: STEP, days, start, sunlight, loss, restart, ocean: { ...PACKAGE, ...MODES[mode] } });
}

const runs = new Map();
function dawns(U, mode) {
  const key = `${U} ${mode}`;
  if (runs.has(key)) return runs.get(key);
  const col = column(U, mode), depths = CELLS.map(() => []);
  while (col.step()) {
    const seconds = Math.round(col.record.time[col.samples - 1] * DAY);
    if (seconds % DAY === DAY / 4) CELLS.forEach((cell, n) => depths[n].push(col.ocean.h[cell]));
  }
  const result = { col, depths, dawn: depths.map((d) => mean(d.slice(-10))), drift: depths.map((d) => Math.abs(mean(d.slice(-5)) - mean(d.slice(-20, -15))) / mean(d.slice(-10))) };
  runs.set(key, result);
  return result;
}

const passing = (h) => transmitted(h) - transmitted(4000);
const absorbedAbove = (h) => 1 - passing(h) - h * deposited(h, 4000);
function trappingDepth(ratio) {
  let low = 0.01, high = 300;
  for (let n = 0; n < 100; n++) { const mid = Math.sqrt(low * high); if (absorbedAbove(mid) > ratio) high = mid; else low = mid; }
  return Math.sqrt(low * high);
}

const show = (values) => values.map((x) => x.toFixed(2)).join(' / ');

test('at 2 m/s the night gate holds the dawn depth within a quarter of the zero-wind trapping depth and within a factor 1.5 of a resolved PWP column, and no shallower than without it', () => {
  const latitude = columnMesh().latCell[CELLS[0]] * 180 / Math.PI;
  const off = dawns(2, 'off'), daily = dawns(2, 'daily'), light = dawns(2, 'light');
  const trap = RATIOS.map(trappingDepth);
  const pwp = RATIOS.map((ratio) => pwpColumn({ wind: 2, sun: sunlight, loss: ratio * MEAN_SUN, latitude }).dawn);
  console.log(`U 2 m/s at ${latitude.toFixed(1)}°N, N/S ${RATIOS.join(' / ')}: dawn depth off ${show(off.dawn)}, retreatDaily ${show(daily.dawn)}, with retreatLight 1 ${show(light.dawn)} m; trapping depth ${show(trap)}, PWP ${show(pwp)} m`);
  RATIOS.forEach((ratio, n) => {
    assert.ok(Math.abs(light.dawn[n] / trap[n] - 1) < 0.25, `N/S ${ratio}: ${light.dawn[n]} m against a trapping depth of ${trap[n]} m`);
    assert.ok(light.dawn[n] < 1.5 * pwp[n] && light.dawn[n] > pwp[n] / 1.5, `N/S ${ratio}: ${light.dawn[n]} m against PWP's ${pwp[n]} m`);
    assert.ok(light.dawn[n] >= off.dawn[n], `N/S ${ratio}: ${light.dawn[n]} m with the gate, ${off.dawn[n]} m without retreatDaily`);
  });
});

test('at 8 m/s the wind holds the dawn depth and the night gate changes it by under a tenth', () => {
  const daily = dawns(8, 'daily'), light = dawns(8, 'light');
  console.log(`U 8 m/s, N/S ${RATIOS.join(' / ')}: dawn depth with retreatDaily ${show(daily.dawn)}, with retreatLight 1 ${show(light.dawn)} m`);
  RATIOS.forEach((ratio, n) => assert.ok(Math.abs(light.dawn[n] / daily.dawn[n] - 1) < 0.1, `N/S ${ratio}: ${light.dawn[n]} against ${daily.dawn[n]} m`));
});

test('under retreatDaily the dawn depth drifts by under 5 % over days 100 to 120', () => {
  for (const [U, mode] of [[2, 'daily'], [2, 'light'], [8, 'daily'], [8, 'light']]) {
    const { drift } = dawns(U, mode);
    console.log(`U ${U} ${mode}: drift ${drift.map((x) => (100 * x).toFixed(1)).join(' / ')} %`);
    drift.forEach((x, n) => assert.ok(x < 0.05, `U ${U} ${mode}, N/S ${RATIOS[n]}: ${x}`));
  }
});

test('the day\'s means are the previous UTC day\'s absorbed sunlight and its dark loss in buoyancy units, b N', () => {
  const { col } = dawns(2, 'light'), { ocean } = col;
  let sum = 0;
  for (let s = 0; s < DAY / STEP; s++) sum += sunlight((s + 0.5) * STEP);
  const expected = sum / (DAY / STEP);
  CELLS.forEach((cell, n) => {
    const b = GRAVITY * thermalExpansion(ocean.T0[cell], ocean.S0[cell]) / RHO_CP;
    const sun = ocean.daily.sunDay[cell], dark = ocean.daily.darkDay[cell] / (b * RATIOS[n] * MEAN_SUN);
    console.log(`N/S ${RATIOS[n]}: S̄ ${sun.toFixed(6)} W/m² against the forcing's ${expected.toFixed(6)}; D̄/(b N) ${dark.toFixed(4)}`);
    assert.ok(Math.abs(sun / expected - 1) < 1e-12, `${sun} against ${expected}`);
    assert.ok(Math.abs(dark - 1) < 0.02, `D̄/(b N) ${dark}`);
  });
});

test('retreatLight needs retreatDaily and penetration and lies in [0, 1], and retreatDaily needs the model time', () => {
  const mesh = columnMesh();
  assert.throws(() => createOcean(mesh, { penetration: 'IB', retreatLight: 1 }), /needs retreatDaily/);
  assert.throws(() => createOcean(mesh, { retreatDaily: true, retreatLight: 1 }), /needs penetration/);
  assert.throws(() => createOcean(mesh, { penetration: 'IB', retreatDaily: true, retreatLight: 2 }), /from 0 \(no gate\) to 1/);
  const ocean = createOcean(mesh, { everySteps: 1, retreatDaily: true }), C = mesh.nCells;
  const surfaceT = new Float64Array(C).fill(288), ice = new Float64Array(C);
  ocean.initialize(surfaceT, ice);
  assert.throws(() => ocean.advance(surfaceT, ice, new Float64Array(C), new Float64Array(mesh.nEdges), 3600), /needs the model time/);
});

test('saved with its daily sums and loaded half way through a day, the column continues across the next two UTC days as it does uninterrupted, and without them it drops the partial day', () => {
  const C = columnMesh().nCells, whole = column(2, 'light', { days: 3 }).run(), first = column(2, 'light', { days: 1.5 }).run();
  const saved = { ...first.ocean.serialize(), ...first.ocean.restartArrays(), Q: first.ocean.Q, W: first.ocean.W, flux: new Float64Array(C) };
  const resumed = column(2, 'light', { start: 1.5, days: 1.5, restart: saved }).run();
  const without = { ...saved };
  for (const name of [...DAILY, 'daySpan']) delete without[name];
  const forgetting = column(2, 'light', { start: 1.5, days: 1.5, restart: without }).run();
  const gap = (a, b) => { let worst = 0; for (let n = 0; n < a.length; n++) worst = Math.max(worst, Math.abs(a[n] - b[n]) / Math.max(1e-30, Math.abs(a[n]))); return worst; };
  const gaps = Object.fromEntries([['h', (o) => o.h], ['Q', (o) => o.Q], ['W', (o) => o.W], ...DAILY.map((name) => [name, (o) => o.daily[name]])].map(([name, pick]) => [name, gap(pick(whole.ocean), pick(resumed.ocean))]));
  console.log(`resumed against uninterrupted, largest relative gaps: ${Object.entries(gaps).map(([name, x]) => `${name} ${x.toExponential(1)}`).join(', ')}; without the sums the mixed layers differ by up to ${gap(whole.ocean.h.subarray(0, C), forgetting.ocean.h.subarray(0, C)).toExponential(1)}`);
  assert.equal(saved.daySpan[0], 0.5 * DAY);
  assert.ok(CELLS.every((cell) => whole.ocean.daily.sunDay[cell] > 0 && resumed.ocean.daily.sunDay[cell] > 0));
  for (const [name, x] of Object.entries(gaps)) assert.ok(x < 1e-9, `${name}: ${x}`);
  assert.ok(gap(whole.ocean.h.subarray(0, C), forgetting.ocean.h.subarray(0, C)) > 1e-3, 'a load without the sums drops the partial day and waits a whole day for its first means');
});

test('with retreatDaily and retreatLight present but off, two days of the CPU model at N = 8 match the model without them bit for bit', () => {
  const N = 8, DT = 900, ocean = { ...PACKAGE, minimumThickness: 10, shallowestMixedDepth: 10 };
  const build = (extra) => {
    const model = createModel(new Grid(N), { levels: sigmaInterfaces('bl36'), ...withCadence({ ocean: { ...ocean, ...extra } }, DT, N) });
    initializeState(model).forEach((values, s) => model.state[s].set(values));
    model.ocean.initialize(model.state[3], model.state[6]);
    for (let n = 0; n < 2 * DAY / DT; n++) model.step(DT);
    return model;
  };
  const absent = build({}), present = build({ retreatDaily: false, retreatLight: 0 });
  const same = (a, b) => Buffer.compare(Buffer.from(a.buffer, a.byteOffset, a.byteLength), Buffer.from(b.buffer, b.byteOffset, b.byteLength)) === 0;
  absent.state.forEach((values, s) => assert.ok(same(values, present.state[s]), `state ${s}`));
  for (const name of ['h', 'u', 'Q', 'W', 'eta', 'capacity', 'transmission']) assert.ok(same(absent.ocean[name], present.ocean[name]), name);
});
