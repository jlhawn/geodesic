import { test } from 'node:test';
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createSeaIce, FREEZING_POINT, MELTING_POINT } from '../js/physics/ice.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { hybridOcean, DT } from './helpers/layered.mjs';

/*
 * The ice and the hybrid ocean under it with iceExchange 'turbulent' (physics/ice.module.js): the water under
 * the ice keeps the heat it holds above freezing, the ice draws on it at McPhee's bottom flux and Maykut and
 * Perovich's lateral rate, and the leads' gain warms it. The ocean and the ice stepped alone, with no air
 * above, are a closed system whose heat, ice volume and salt are checked step by step; a restart carries
 * what the exchange adds; and the default exchange is the model it was.
 */
const POLAR = 62 * Math.PI / 180;

function rig(iceExchange, { thickness = 1.2, cover = 0.8, flux = 0, global = false } = {}) {
  const r = hybridOcean({ iceExchange, everySteps: 1 }, () => 4000, 290), { mesh, C, ocean, surfaceT, ice } = r;
  r.seaIce = createSeaIce(mesh, { heatCapacity: ocean.capacity, iceExchange: ocean.iceExchange });
  ocean.useSeaIce(r.seaIce.concentration, r.seaIce);
  for (let i = 0; i < C; i++) if (global || Math.abs(mesh.latCell[i]) > POLAR) { ice[i] = thickness; surfaceT[i] = FREEZING_POINT - 3; r.seaIce.concentration[i] = cover; }
  ocean.initialize(surfaceT, ice, { concentration: r.seaIce.concentration });
  for (let i = 0; i < C; i++) if (!(ice[i] > 0)) surfaceT[i] = ocean.T0[i];
  r.heat = new Float64Array(C).fill(flux);
  r.volume = () => { let v = 0; for (let i = 0; i < C; i++) v += mesh.areaCell[i] * ice[i] * (ice[i] > 0 ? r.seaIce.concentration[i] : 0); return v; };
  r.cellEnergy = () => Float64Array.from({ length: C }, (_, i) => {
    let q = 0;
    for (let k = 1; k < ocean.layers; k++) q += ocean.Q[k * C + i] - 272 * ocean.h[k * C + i];
    const iced = ice[i] > 0, top = (iced ? ocean.h[i] * FREEZING_POINT + (r.seaIce.turbulent ? r.seaIce.waterHeat[i] : 0) / ocean.rhoCp : ocean.h[i] * surfaceT[i]) - 272 * ocean.h[i];
    return ocean.rhoCp * (q + top) + (iced ? r.seaIce.energy(surfaceT[i], ice[i], r.seaIce.snow[i], r.seaIce.cover(i, ice[i])) : 0);
  });
  r.cycle = () => {
    ocean.readSurface(surfaceT, ice, r.seaIce.concentration);
    ocean.salt();
    ocean.writeSurface(surfaceT, r.seaIce.oceanFlux, DT);
    for (let i = 0; i < C; i++) r.seaIce.update(surfaceT, ice, r.heat, i, DT);
  };
  r.seed = (rise = 1, cold = 0.2) => {
    let n = 0;
    for (let i = 0; i < C; i++) if (ice[i] > 0) r.seaIce.waterHeat[i] = ocean.capacity[i] * (n++ % 7 === 3 ? -cold : rise);
  };
  r.step = () => {
    ocean.advance(surfaceT, ice, r.seaIce.oceanFlux, r.calm, DT, r.seaIce.concentration);
    for (let i = 0; i < C; i++) r.seaIce.update(surfaceT, ice, r.heat, i, DT);
  };
  return r;
}

test('the ocean\'s surface and the ice conserve heat exactly through bottom melt, lateral melt, the leads\' gain and the freezing of supercooled water, and the ice volume closes with its budget', () => {
  const r = rig('turbulent', { flux: 40, global: true }), { mesh, C, ocean, seaIce, ice } = r;
  const area = mesh.areaCell.reduce((s, a) => s + a, 0), input = DT * 40 * area;
  r.seed(1.5);
  let before = r.cellEnergy(), moved = 0, warm = 0, worst = 0, supercooled = 0;
  const volume0 = r.volume();
  for (let n = 0; n < 48; n++) {
    r.cycle();
    const after = r.cellEnergy();
    let change = 0;
    for (let i = 0; i < C; i++) change += mesh.areaCell[i] * (after[i] - before[i]);
    before = after;
    worst = Math.max(worst, Math.abs(change - input) / input);
    assert.ok(Math.abs(change - input) < 1e-9 * input, `cycle ${n}: the heat changed by ${change / input} of the air's input`);
    for (let i = 0; i < C; i++) if (ice[i] > 0) { moved += mesh.areaCell[i] * seaIce.waterHeat[i]; warm = Math.max(warm, ocean.T0[i] - FREEZING_POINT); }
    supercooled += seaIce.oceanFlux.reduce((s, x) => s + (x < 0 ? 1 : 0), 0);
  }
  const grown = seaIce.budget.frozen - seaIce.budget.melted;
  assert.ok(Math.abs(r.volume() - volume0 - grown) < 1e-9 * volume0, `the ice volume went from ${volume0} to ${r.volume()} and the budget says ${grown}`);
  console.log(`48 cycles of the ocean's surface and its ice under 40 W/m²: heat closes to ${worst.toExponential(1)} of a cycle's input; the water under the ice reached ${warm.toFixed(2)} K above freezing; ice volume ${(volume0 / 1e12).toFixed(2)} → ${(r.volume() / 1e12).toFixed(2)} thousand km³ (${(grown / 1e12).toFixed(3)} by the budget, ${(seaIce.budget.lateralMelted / 1e12).toFixed(3)} melted at the sides or by leads); ${supercooled} cell-cycles freezing supercooled water`);
  assert.ok(moved > 0 && warm > 0 && supercooled > 0, 'the water under the ice held heat above freezing and some was supercooled');
  assert.ok(seaIce.budget.lateralMelted > 0 && seaIce.budget.melted + seaIce.budget.frozen !== 0);
});

test('the salt the ocean takes from the ice follows the ice volume that the exchange changed, cell by cell', () => {
  const r = rig('turbulent', { flux: 30, global: true }), { mesh, C, ocean, seaIce, ice } = r;
  const columnSalt = (i) => { let w = 0; for (let k = 0; k < ocean.layers; k++) w += ocean.W[k * C + i]; return w; };
  const volumeAt = (i) => (ice[i] > 0 ? seaIce.concentration[i] : 0) * ice[i];
  r.seed(1.5);
  let previous = Float64Array.from({ length: C }, (_, i) => volumeAt(i)), worst = 0, moved = 0;
  for (let n = 0; n < 24; n++) {
    const volume = Float64Array.from({ length: C }, (_, i) => volumeAt(i)), before = Float64Array.from({ length: C }, (_, i) => columnSalt(i));
    ocean.readSurface(r.surfaceT, ice, seaIce.concentration);
    const salinity = Float64Array.from({ length: C }, (_, i) => ocean.W[i] / ocean.h[i]);
    ocean.salt();
    for (let i = 0; i < C; i++) {
      const expected = (salinity[i] - 5) * 917 * (volume[i] - previous[i]) / 1000, got = columnSalt(i) - before[i];
      if (n > 0 && i === 0) assert.ok(Math.abs(volume[i] - previous[i]) > 0, 'the polar cell\'s ice changed');
      worst = Math.max(worst, Math.abs(got - expected)); moved += Math.abs(expected);
    }
    previous = volume;
    ocean.writeSurface(r.surfaceT, seaIce.oceanFlux, DT);
    for (let i = 0; i < C; i++) seaIce.update(r.surfaceT, ice, r.heat, i, DT);
  }
  console.log(`24 cycles: the column salt followed (S − 5) ρ_i ΔV/1000 cell by cell to ${worst.toExponential(1)} psu m, against ${moved.toExponential(2)} psu m moved in all and the columns' 1.4e5`);
  assert.ok(moved > 0 && worst < 1e-9, `salt differs by ${worst} psu m`);
});

test('the turbulent exchange keeps the water warm under the ice and the ice thicker than the instant exchange, which melts what the layer holds at once', () => {
  const instant = rig('instant', { flux: 40 }), turbulent = rig('turbulent', { flux: 40 });
  for (let n = 0; n < 24; n++) { instant.step(); turbulent.step(); }
  let warm = 0, flat = 0;
  for (let i = 0; i < turbulent.C; i++) {
    if (turbulent.ice[i] > 0) warm = Math.max(warm, turbulent.ocean.T0[i] - FREEZING_POINT);
    if (instant.ice[i] > 0) flat = Math.max(flat, instant.ocean.T0[i] - FREEZING_POINT);
  }
  console.log(`after 24 steps the ice volume is ${(instant.volume() / 1e12).toFixed(3)} (instant) and ${(turbulent.volume() / 1e12).toFixed(3)} (turbulent) thousand km³; the water under the ice is ${flat.toExponential(1)} and ${warm.toFixed(2)} K above freezing`);
  assert.equal(flat, 0, 'the instant exchange holds the water at freezing');
  assert.ok(warm > 0.01, `the turbulent exchange leaves it ${warm} K above freezing`);
  assert.ok(turbulent.volume() > instant.volume(), `ice volume ${turbulent.volume()} against ${instant.volume()}`);
});

test('a restart carries what the exchange adds: the water\'s heat and the conductance resume the run bit for bit, and without them the ice resumes differently', () => {
  const run = (r, steps) => { for (let n = 0; n < steps; n++) r.step(); return r; };
  const whole = run(rig('turbulent', { flux: 40 }), 12), first = run(rig('turbulent', { flux: 40 }), 6);
  assert.ok(first.seaIce.waterHeat.some((x) => x > 0), 'the water holds heat at the restart');
  const resume = (omit) => {
    const second = rig('turbulent', { flux: 40 }), saved = { ...first.ocean.serialize(), ...first.ocean.restartArrays() };
    for (const name of omit) delete saved[name];
    second.surfaceT.set(first.surfaceT); second.ice.set(first.ice);
    for (const name of ['concentration', 'snow', 'snowAlbedo', 'oceanFlux']) second.seaIce[name].set(first.seaIce[name]);
    second.ocean.load(saved, second.surfaceT, second.ice, second.seaIce.concentration);
    return second;
  };
  const same = (a, b, what) => {
    const at = Array.from(a).findIndex((x, i) => x !== b[i]);
    assert.equal(at, -1, `${what}[${at}]: ${a[at]} whole against ${b[at]} restarted`);
  };
  const compare = (b) => {
    for (const [name, pick] of [['surfaceT', (r) => r.surfaceT], ['ice', (r) => r.ice], ['concentration', (r) => r.seaIce.concentration], ['snow', (r) => r.seaIce.snow], ['waterHeat', (r) => r.seaIce.waterHeat], ['oceanFlux', (r) => r.seaIce.oceanFlux], ['Q', (r) => r.ocean.Q], ['W', (r) => r.ocean.W], ['h', (r) => r.ocean.h]]) same(pick(whole), pick(b), name);
  };
  compare(run(resume([]), 6));
  console.log('the restarted second half matches the whole run bit for bit');
  const lost = run(resume(['iceWater', 'iceConductance']), 6);
  assert.throws(() => compare(lost), /whole against/, 'dropping the water\'s heat changes the run');
});

test('the ocean and the sea ice agree on the exchange, and the model builds the ice from the ocean\'s option', () => {
  const grid = new Grid(16), model = createModel(grid, { ocean: { iceExchange: 'turbulent' } });
  assert.equal(model.seaIce.turbulent, true);
  assert.equal(model.ocean.iceExchange.mode, 'turbulent');
  assert.equal(createModel(grid, { ocean: {} }).seaIce.turbulent, false);
  assert.throws(() => createModel(grid, { ocean: { iceExchange: 'turbulent' }, ice: { iceExchange: 'instant' } }), /iceExchange/);
  assert.throws(() => createModel(grid, { ocean: { iceExchange: 'gradual' } }), /iceExchange/);
});

test('with the default exchange the coupled model and the ice are the model before iceExchange, bit for bit', () => {
  const sha = (arrays) => { const h = createHash('sha256'); for (const a of arrays) h.update(Buffer.from(Float64Array.from(a).buffer)); return h.digest('hex').slice(0, 16); };
  const m = createModel(new Grid(16), { ocean: { everySteps: 4 } });
  initializeState(m, {}).forEach((v, a) => m.state[a].set(v));
  m.ocean.initialize(m.state[3], m.state[6]);
  for (let n = 0; n < 48; n++) m.step(1350);
  const coupled = sha([m.state[3], m.state[6], m.seaIce.concentration, m.seaIce.snow, m.seaIce.snowAlbedo, m.seaIce.oceanFlux, m.ocean.Q, m.ocean.W, m.ocean.h, m.ocean.u]);
  const sea = createSeaIce(m.mesh), surfaceT = new Float64Array([FREEZING_POINT - 5]), ice = new Float64Array([0.6]), flux = new Float64Array(1);
  sea.concentration[0] = 0.5; sea.snow[0] = 5; sea.oceanFlux[0] = 4;
  const seen = [];
  for (let n = 0; n < 200 * 96; n++) {
    const phase = n < 40 * 96 ? 'winter' : n < 120 * 96 ? 'summer' : 'autumn';
    flux[0] = phase === 'summer' ? 160 : phase === 'winter' ? -120 : -200;
    if (phase === 'winter' && n % 96 === 0) sea.deposit(0, 1.5, MELTING_POINT - 10, ice, surfaceT);
    sea.update(surfaceT, ice, flux, 0, 900, phase === 'summer' ? 120 : 0);
    if (n % 500 === 0) seen.push(surfaceT[0], ice[0], sea.concentration[0], sea.snow[0]);
  }
  const column = sha([seen, [surfaceT[0], ice[0], sea.concentration[0], sea.snow[0], sea.budget.leadFrozen, sea.budget.lateralMelted, sea.budget.melted, sea.budget.frozen]]);
  assert.equal(coupled, '9e6191dca0b1752c');
  assert.equal(column, '3819c4d3c156057d');
});
