import { test } from 'node:test';
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { readFileSync } from 'node:fs';
import { createSeaIce, FREEZING_POINT, MELTING_POINT, COVER_GRAIN, DEFAULT_ICE_EXCHANGE, frictionMinimumNote } from '../js/physics/ice.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { DT, mesh as oceanMesh } from './helpers/layered.mjs';
import { rig, POLAR } from './helpers/iceExchange.mjs';

/*
 * The ice and the hybrid ocean under it with iceExchange 'turbulent' (physics/ice.module.js): the water under
 * the ice keeps the heat it holds above freezing, the ice draws on it at McPhee's bottom flux and Maykut and
 * Perovich's lateral rate, and the leads' gain warms it. The ocean and the ice stepped alone, with no air
 * above, are a closed system whose heat, ice volume and salt are checked step by step; a restart carries
 * what the exchange adds; and the default exchange keeps the default model's pinned hashes.
 */

test('the ocean\'s surface and the ice conserve heat exactly through bottom melt, lateral melt, the lateral melt dropped below the f32 grain over water barely above freezing, the leads\' gain and the freezing of supercooled water, and the ice volume closes with its budget', () => {
  const r = rig('turbulent', { flux: 40, global: true }), { mesh, C, ocean, seaIce, ice } = r;
  const area = mesh.areaCell.reduce((s, a) => s + a, 0), input = DT * 40 * area;
  r.seed(1.5);
  const { lateralRate, lateralExponent } = seaIce.iceExchange;
  let fine = 0;
  for (let i = 0; i < C; i++) {
    if (!(ice[i] > 0) || i % 3 !== 1) continue;
    seaIce.waterHeat[i] = ocean.capacity[i] * 0.005;
    if (lateralRate * Math.pow(0.005, lateralExponent) * DT < COVER_GRAIN) fine++;
  }
  assert.ok(fine > 0, 'some cells hold water whose lateral melt falls below the grain');
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
  console.log(`48 cycles of the ocean's surface and its ice under 40 W/m², ${fine} iced cells starting 0.005 K above freezing where the lateral melt drops out: heat closes to ${worst.toExponential(1)} of a cycle's input; the water under the ice reached ${warm.toFixed(2)} K above freezing; ice volume ${(volume0 / 1e12).toFixed(2)} → ${(r.volume() / 1e12).toFixed(2)} thousand km³ (${(grown / 1e12).toFixed(3)} by the budget, ${(seaIce.budget.lateralMelted / 1e12).toFixed(3)} melted at the sides or by leads); ${supercooled} cell-cycles freezing supercooled water`);
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

test('a lean turbulent save, serialize() without the restart arrays, carries the water\'s heat under the ice and the conductance, and a turbulent ocean loads them as they were saved; an instant save carries neither', () => {
  const run = rig('turbulent', { flux: 40, global: true }), { C } = run;
  for (let n = 0; n < 4; n++) run.cycle();
  const lean = run.ocean.serialize();
  assert.equal(lean.Q, undefined, 'the save is lean');
  for (const name of ['iceWater', 'iceConductance']) assert.equal(lean[name]?.length, C, `the lean save holds ${name}`);
  assert.ok(lean.iceWater.some((x) => x > 0) && lean.iceConductance.some((x) => x > 0), 'the lean save holds heat and a conductance under the ice');
  const back = rig('turbulent', { flux: 40, global: true });
  back.ice.set(run.ice); back.surfaceT.set(run.surfaceT); back.seaIce.concentration.set(run.seaIce.concentration);
  back.ocean.load(lean, back.surfaceT, back.ice, back.seaIce.concentration);
  assert.deepEqual(Array.from(back.seaIce.waterHeat), Array.from(lean.iceWater));
  assert.deepEqual(Array.from(back.seaIce.conductance), Array.from(lean.iceConductance));
  const plain = rig('instant', { flux: 40, global: true });
  plain.cycle();
  assert.ok(!('iceWater' in plain.ocean.serialize()) && !('iceConductance' in plain.ocean.serialize()));
});

test('a turbulent save, restart or lean, loaded into an explicitly instant ocean books the water\'s heat under the ice into the top layer, logged, and its first read keeps that heat where an instant read sets the freezing point under ice; held heat over a cell without ice is dropped, logged', () => {
  const run = rig('turbulent', { flux: 40, global: true });
  for (let n = 0; n < 4; n++) run.cycle();
  for (const restart of [true, false]) {
    const saved = run.ocean.serialize({ restart }), { iceWater, iceConductance, ...without } = saved;
    assert.ok(iceWater.some((x) => x > 0) && iceConductance.some((x) => x > 0), `the ${restart ? 'restart' : 'lean'} save holds heat and a conductance under the ice`);
    const [taken, plain] = [saved, without].map((given) => {
      const r = rig('instant', { flux: 40, global: true });
      r.ice.set(run.ice); r.surfaceT.set(run.surfaceT); r.seaIce.concentration.set(run.seaIce.concentration);
      r.ice[[...Array(r.C).keys()].find((i) => r.ocean.cellOcean[i] && iceWater[i] > 0)] = 0;
      const notes = [], warn = console.warn, booked = [...Array(r.C).keys()].filter((i) => r.ocean.cellOcean[i] && r.ice[i] > 0);
      console.warn = (note) => notes.push(note);
      try { r.ocean.load(given, r.surfaceT, r.ice, r.seaIce.concentration); } finally { console.warn = warn; }
      if (given === saved) {
        assert.equal(notes.length, 2, 'the booking and the drop are logged once each');
        assert.match(notes[0], new RegExp(`under ${booked.length} iced cells`));
        assert.match(notes[1], /dropped \d+\.\d+ EJ of the saved ocean\.iceWater, held over 1 cells without ice over an ocean column/);
        for (const i of booked) {
          assert.equal(r.ocean.Q[i], r.ocean.h[i] * FREEZING_POINT + iceWater[i] / r.ocean.rhoCp);
          assert.equal(r.ocean.T0[i], r.ocean.Q[i] / r.ocean.h[i]);
        }
      } else assert.deepEqual(notes, [], 'a save without iceWater books nothing');
      r.ocean.readSurface(r.surfaceT, r.ice, r.seaIce.concentration);
      return { r, booked };
    });
    const top = (r, i) => r.ocean.Q[i] - r.ocean.h[i] * FREEZING_POINT;
    const kept = Math.max(...taken.booked.map((i) => Math.abs(top(taken.r, i) - iceWater[i] / taken.r.ocean.rhoCp) / Math.max(1, Math.abs(taken.r.ocean.Q[i])))), plainTop = Math.max(...plain.booked.map((i) => Math.abs(top(plain.r, i))));
    assert.ok(kept < 1e-14 && plainTop === 0, `the first read keeps the booked heat to ${kept} of the layer's heat, and a load without it reads the freezing point to ${plainTop}`);
    taken.r.ocean.readSurface(taken.r.surfaceT, taken.r.ice, taken.r.seaIce.concentration);
    assert.equal(Math.max(...taken.booked.map((i) => Math.abs(top(taken.r, i)))), 0, 'a second read with no write between reads the freezing point: the booked heat is read once');
    assert.ok(!('iceWater' in taken.r.ocean.serialize({ restart })), 'the instant ocean saves no water heat');
  }
});

test('ice on a cell with no ocean under it, as on a lake, holds no heat under the turbulent exchange and takes the leads\' gain as the instant exchange does, and the ocean stays finite there', () => {
  const { nCells, nEdges, latCell, cellsOnEdge } = oceanMesh, land = Uint8Array.from({ length: nCells }, (_, i) => (latCell[i] > POLAR && i % 4 === 0 ? 1 : 0));
  const geography = { land, elevation: Float64Array.from(land, (l) => (l ? 100 : -4000)), edgeOcean: Uint8Array.from({ length: nEdges }, (_, e) => (land[cellsOnEdge[2 * e]] || land[cellsOnEdge[2 * e + 1]] ? 0 : 1)) };
  const [instant, turbulent] = ['instant', 'turbulent'].map((mode) => {
    const r = rig(mode, { geography, flux: 40 }), { C, ocean, ice, surfaceT, seaIce } = r;
    for (let n = 0; n < 8; n++) r.step();
    const dry = [];
    let held = 0;
    for (let i = 0; i < C; i++) if (ice[i] > 0) { if (!ocean.cellOcean[i]) dry.push(i); else if (seaIce.turbulent && seaIce.waterHeat[i] > 0) held++; }
    assert.ok(ocean.serialize().T.every(Number.isFinite), `${mode}: every layer's temperature is finite`);
    const at = (a) => dry.map((i) => a[i]);
    return { held, dry: { cells: dry, ice: at(ice), cover: at(seaIce.concentration), surfaceT: at(surfaceT), snow: at(seaIce.snow), water: seaIce.turbulent ? at(seaIce.waterHeat) : dry.map(() => 0), T0: at(ocean.T0), Q: at(ocean.Q) } };
  });
  assert.ok(turbulent.dry.cells.length > 0 && turbulent.held > 0, `${turbulent.dry.cells.length} iced cells with no ocean under them; ${turbulent.held} over the ocean hold heat`);
  assert.deepEqual(turbulent.dry, instant.dry);
  for (const t of turbulent.dry.T0) assert.equal(t, FREEZING_POINT);
});

test('the ocean and the sea ice agree on the exchange: the model builds the ice from the ocean\'s option, refuses an ice option of the other mode, and gives the ice the ocean\'s constants', () => {
  const grid = new Grid(16), model = createModel(grid, { ocean: { iceExchange: 'turbulent' } });
  assert.equal(model.seaIce.turbulent, true);
  assert.equal(model.ocean.iceExchange.mode, 'turbulent');
  assert.equal(createModel(grid, { ocean: { iceExchange: 'instant' } }).seaIce.turbulent, false);
  assert.throws(() => createModel(grid, { ocean: { iceExchange: 'turbulent' }, ice: { iceExchange: 'instant' } }), /the ocean's iceExchange is turbulent and the sea ice's is instant/);
  assert.throws(() => createModel(grid, { ocean: { iceExchange: 'gradual' } }), /iceExchange/);
  const constants = createModel(grid, { ocean: { iceExchange: { mode: 'turbulent', frictionMinimum: 2e-2 } }, ice: { iceExchange: 'turbulent' } });
  assert.deepEqual([constants.seaIce.iceExchange.frictionMinimum, constants.ocean.iceExchange.frictionMinimum], [2e-2, 2e-2]);
});

test('with the instant exchange the coupled model over 48 steps and a 200-day ice column keep their pinned hashes', () => {
  const sha = (arrays) => { const h = createHash('sha256'); for (const a of arrays) h.update(Buffer.from(Float64Array.from(a).buffer)); return h.digest('hex').slice(0, 16); };
  const m = createModel(new Grid(16), { ocean: { everySteps: 4, iceExchange: 'instant' } });
  initializeState(m, {}).forEach((v, a) => m.state[a].set(v));
  m.ocean.initialize(m.state[3], m.state[6]);
  for (let n = 0; n < 48; n++) m.step(1350);
  const coupled = sha([m.state[3], m.state[6], m.seaIce.concentration, m.seaIce.snow, m.seaIce.snowAlbedo, m.seaIce.oceanFlux, m.ocean.Q, m.ocean.W, m.ocean.h, m.ocean.u]);
  const sea = createSeaIce(m.mesh, { iceExchange: 'instant' }), surfaceT = new Float64Array([FREEZING_POINT - 5]), ice = new Float64Array([0.6]), flux = new Float64Array(1);
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
  assert.equal(coupled, 'd25e8222cf1e1ff8');
  assert.equal(column, '3819c4d3c156057d');
});

test('an instant save taken as its ocean hands the ice heat loads into the turbulent exchange holding no heat under the ice and keeping the heat in transit, which the ice takes as the instant ice does; the heat then closes step by step and no ice collapses', () => {
  const instant = rig('instant', { flux: 40, global: true }), turbulent = rig('turbulent', { flux: 40, global: true }), { mesh, C } = instant;
  const area = mesh.areaCell.reduce((s, a) => s + a, 0), input = DT * 40 * area;
  for (let n = 0; n < 4; n++) instant.cycle();
  instant.ocean.readSurface(instant.surfaceT, instant.ice, instant.seaIce.concentration);
  instant.ocean.salt();
  for (let i = 0; i < C; i++) if (instant.ice[i] > 0) instant.ocean.Q[i] += instant.ocean.h[i] * (0.2 + 0.1 * (i % 3));
  instant.ocean.writeSurface(instant.surfaceT, instant.seaIce.oceanFlux, DT);
  const saved = instant.ocean.serialize({ restart: true });
  assert.ok(!('iceWater' in saved) && saved.flux.every((x, i) => !(instant.ice[i] > 0) || x > 0), 'an instant save whose every iced cell has heat in transit to the ice');
  turbulent.ice.set(instant.ice); turbulent.surfaceT.set(instant.surfaceT); turbulent.seaIce.concentration.set(instant.seaIce.concentration);
  turbulent.ocean.load(saved, turbulent.surfaceT, turbulent.ice, turbulent.seaIce.concentration);
  const held = [...turbulent.seaIce.waterHeat.keys()].filter((i) => turbulent.ice[i] > 0 && turbulent.seaIce.waterHeat[i] !== 0);
  assert.deepEqual(held, [], 'iced cells holding heat at the load');
  assert.deepEqual(Array.from(turbulent.seaIce.oceanFlux), Array.from(saved.flux), 'the load keeps the heat in transit');
  const stored = (r) => r.cellEnergy().reduce((s, e, i) => s + mesh.areaCell[i] * e, 0), transit = saved.flux.reduce((s, x, i) => s + mesh.areaCell[i] * x * DT, 0);
  const atSave = stored(instant) + transit, volume0 = Float64Array.from(instant.ice, (h, i) => h * instant.seaIce.concentration[i]);
  assert.ok(Math.abs(stored(turbulent) + transit - atSave) < 1e-12 * atSave, 'the load holds the save\'s heat');
  let before = atSave, worst = 0;
  for (let n = 0; n < 9; n++) {
    for (const r of [instant, turbulent]) if (n === 0) { for (let i = 0; i < C; i++) r.seaIce.update(r.surfaceT, r.ice, r.heat, i, DT); } else r.cycle();
    const after = stored(turbulent);
    worst = Math.max(worst, Math.abs(after - before - input) / transit);
    before = after;
  }
  const volume = (r) => Float64Array.from(r.ice, (h, i) => h * r.seaIce.concentration[i]), [vi, vt] = [volume(instant), volume(turbulent)];
  const kept = Math.min(...[...vt.keys()].filter((i) => volume0[i] > 0).map((i) => vt[i] / volume0[i])), ahead = Math.min(...[...vt.keys()].map((i) => vt[i] - vi[i]));
  console.log(`transit ${(transit / area).toExponential(3)} J/m², the heat closed to ${worst.toExponential(2)} of it over 9 steps; the turbulent ice keeps at least ${kept.toFixed(6)} of its volume at the save, at least ${ahead.toExponential(2)} m above the instant ice's in every cell`);
  assert.ok(worst < 1e-9, `the heat closes to ${worst} of the transit`);
  assert.ok(kept > 0.95 && ahead >= 0, `kept ${kept}, above the instant ice by ${ahead} m`);
});

test('a turbulent exchange given without its own frictionMinimum, as the bare string or an object, is noted while its 5e-4 is not the default\'s, and spinup, oceanSpinup, explainOcean and the page worker log that note', () => {
  assert.equal(frictionMinimumNote('turbulent'), 'iceExchange turbulent takes frictionMinimum 5e-4; the default is 2e-2');
  assert.equal(frictionMinimumNote({ mode: 'turbulent', heatTransfer: 0.006 }), frictionMinimumNote('turbulent'));
  for (const exchange of [undefined, 'instant', { mode: 'turbulent', frictionMinimum: 5e-4 }, DEFAULT_ICE_EXCHANGE]) assert.equal(frictionMinimumNote(exchange), null, JSON.stringify(exchange));
  for (const path of ['scripts/spinup.mjs', 'scripts/oceanSpinup.mjs', 'scripts/explainOcean.mjs', 'js/model.worker.js']) assert.match(readFileSync(new URL(`../${path}`, import.meta.url), 'utf8'), /frictionMinimumNote\(/, path);
});
