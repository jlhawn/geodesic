import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createOcean } from '../js/ocean/layered.module.js';
import { iceDriftMoves, ICE_REMNANT } from '../js/ocean/iceDrift.module.js';
import { createSeaIce, seaIceRemoval, FREEZING_POINT, MINIMUM_CONCENTRATION, MINIMUM_VOLUME, NEGLIGIBLE_CONCENTRATION, NEGLIGIBLE_VOLUME } from '../js/physics/ice.module.js';
import { DEG, meshOf, rig, northwardOn, nearestCell, windOn, unit, dot } from './helpers/iceDrift.mjs';

/*
 * Ice that the drift carries into open cells arrives a few thousandths of a
 * cell at a time (a Courant number of 0.005 per ocean call at N = 128).
 * With the sea ice's `keepSmallIce` it persists, freezes the leads around it
 * and melts only in the water's heat; without it the removal floors of the
 * thermodynamics (MINIMUM_CONCENTRATION, MINIMUM_VOLUME) melt it on arrival.
 */
const near = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const DT = 3600;

function pack(keepSmallIce, { courant = 0.005, water = 0, edge = -65, dt = DT } = {}) {
  const m = meshOf(16), r = rig(m, { rheology: 'frozen' }, null, {}, { keepSmallIce }), toward = northwardOn(m, (lat) => lat < -50 * DEG);
  r.dt = dt;
  for (let i = 0; i < r.C; i++) {
    r.surfaceT[i] = FREEZING_POINT + water;
    if (m.latCell[i] < edge * DEG) { r.thickness[i] = 1; r.concentration[i] = 0.9; r.snow[i] = 10; r.surfaceT[i] = FREEZING_POINT - 4; }
  }
  for (let e = 0; e < r.E; e++) r.drift.velocity[e] = toward[e] * courant * m.dcEdge[e] / dt;
  r.initial = Uint8Array.from(r.thickness, (h) => (h > 0 ? 1 : 0));
  r.volume = () => { let v = 0; for (let i = 0; i < r.C; i++) v += m.areaCell[i] * r.thickness[i] * r.concentration[i]; return v; };
  r.cells = (test) => { let n = 0; for (let i = 0; i < r.C; i++) if (r.thickness[i] > 0 && test(i)) n++; return n; };
  r.carry = () => { r.drift.arrived.fill(0); r.drift.transport(r.env()); };
  return r;
}

test('the removal floors are the minimum ones unless the drift moves the ice, in the sea ice and in the model that builds it', () => {
  assert.deepEqual(seaIceRemoval(false), { concentration: MINIMUM_CONCENTRATION, volume: MINIMUM_VOLUME });
  assert.deepEqual(seaIceRemoval(true), { concentration: NEGLIGIBLE_CONCENTRATION, volume: NEGLIGIBLE_VOLUME });
  assert.ok(NEGLIGIBLE_CONCENTRATION < 1e-3 * MINIMUM_CONCENTRATION && NEGLIGIBLE_VOLUME < MINIMUM_VOLUME);
  assert.deepEqual(createSeaIce(meshOf(16)).removal, seaIceRemoval(false));
  assert.deepEqual(createSeaIce(meshOf(16), { keepSmallIce: true }).removal, seaIceRemoval(true));
  for (const option of [true, {}, { rheology: 'free' }, { rheology: 'cavitating', airStress: 'given' }]) assert.ok(iceDriftMoves(option), JSON.stringify(option));
  for (const option of [false, undefined, null, { rheology: 'frozen' }]) assert.ok(!iceDriftMoves(option), JSON.stringify(option));
  const grid = new Grid(16), removal = (ocean, ice = {}) => createModel(grid, { ocean, ice }).seaIce.removal;
  assert.deepEqual(removal({ iceDrift: { rheology: 'cavitating' } }), seaIceRemoval(true));
  assert.deepEqual(removal({ iceDrift: true }), seaIceRemoval(true));
  assert.deepEqual(removal({ iceDrift: { rheology: 'frozen' } }), seaIceRemoval(false));
  assert.deepEqual(removal({ iceDrift: false }), seaIceRemoval(false));
  assert.deepEqual(removal({}), seaIceRemoval(true));
  assert.deepEqual(removal(false), seaIceRemoval(false));
  assert.deepEqual(removal({ iceDrift: true }, { keepSmallIce: false }), seaIceRemoval(false));
  assert.throws(() => { const r = rig(meshOf(16), { rheology: 'frozen' }); r.drift.attach({ snow: r.seaIce.snow, meltAway: r.seaIce.meltAway }); }, /needs the sea ice/);
});

test('ice drifting into open water at a Courant number of 0.005 advances the edge and the dynamic term sums to zero, where the minimum floors melt it on arrival', () => {
  const kept = pack(true), melted = pack(false), V0 = kept.volume(), n0 = kept.cells(() => true);
  assert.equal(melted.cells(() => true), n0);
  const calls = 40;
  let balance = 0, removedKept = 0, removedMelted = 0;
  for (const [r, name] of [[kept, 'kept'], [melted, 'melted']]) {
    let removed = 0;
    for (let k = 0; k < calls; k++) {
      const before = r.volume();
      r.carry();
      removed += r.drift.count.removed;
      balance = Math.max(balance, Math.abs(r.volume() + r.drift.count.removed - before) / V0);
    }
    if (name === 'kept') removedKept = removed; else removedMelted = removed;
  }
  const fringe = (r) => r.cells((i) => !r.initial[i]), trace = (r) => r.cells((i) => r.concentration[i] < MINIMUM_CONCENTRATION);
  console.log(`${calls} calls at a Courant number of 0.005: ${fringe(kept)} cells beyond the pack of ${n0} with ${trace(kept)} under 1 percent cover and ${(removedKept / V0).toExponential(2)} of the volume removed; the minimum floors leave ${fringe(melted)} and melt ${(removedMelted / V0).toExponential(2)} on arrival`);
  assert.ok(balance < 1e-13, `volume balance of a call ${balance}`);
  assert.ok(fringe(kept) >= 20 && trace(kept) >= 10, `${fringe(kept)} cells beyond the pack, ${trace(kept)} under 1 percent`);
  assert.equal(fringe(melted), 0, 'the minimum floors keep the edge where it was');
  assert.ok(removedKept < 1e-6 * V0, `removed ${removedKept / V0}`);
  assert.ok(removedMelted > 1e-2 * V0, `melted on arrival ${removedMelted / V0}`);
  assert.ok(kept.concentration.every((a) => a === 0 || (a >= NEGLIGIBLE_CONCENTRATION && a <= 1)));
  near(kept.volume() / V0, 1, 1e-6, 'the ice volume of the globe');
  near(kept.seaIce.budget.melted, removedKept, 1e-9 * removedKept + 1e-6, 'what the budget counts as melted is what the drift removed');
});

test('the thermodynamics keep the arrived ice and freeze the leads around it, where the minimum floors remove ice below 1 percent cover in a step', () => {
  const m = meshOf(16), seed = nearestCell(m, -62 * DEG, 30 * DEG);
  const step = (keep, flux) => {
    const s = createSeaIce(m, { keepSmallIce: keep }), surfaceT = new Float64Array(m.nCells).fill(FREEZING_POINT), h = new Float64Array(m.nCells), fluxes = new Float64Array(m.nCells).fill(flux);
    h[seed] = 1; s.concentration[seed] = 0.005; s.snow[seed] = 8; surfaceT[seed] = FREEZING_POINT - 4;
    s.update(surfaceT, h, fluxes, seed, DT);
    return { h: h[seed], A: s.concentration[seed], snow: s.snow[seed] };
  };
  const cold = step(true, -40), warm = step(true, 100), old = step(false, -40);
  console.log(`a cell of 0.005 cover and 1 m: after a cold step ${cold.A.toFixed(5)} (${cold.h.toFixed(3)} m), after a warm step ${warm.A.toFixed(5)}; with the minimum floors ${old.A}`);
  assert.ok(cold.h > 0 && cold.A > 0.005 && cold.snow > 0, 'the leads freeze onto it');
  assert.ok(warm.h > 0 && warm.A < 0.005 && warm.A > NEGLIGIBLE_CONCENTRATION, 'it melts from the sides by the lead heat and stays');
  assert.equal(old.h, 0); assert.equal(old.A, 0);

  const kept = pack(true, { water: 0.6 }), melted = pack(false, { water: 0.6 }), flux = new Float64Array(kept.C).fill(-30), growth = [];
  const beyond = (r) => Array.from({ length: r.C }, (_, i) => i).filter((i) => !r.initial[i] && r.thickness[i] > 0);
  for (let n = 0; n < 6; n++) {
    for (const r of [kept, melted]) {
      for (let i = 0; i < r.C; i++) r.seaIce.update(r.surfaceT, r.thickness, flux, i, DT);
      r.carry();
    }
    const fringe = beyond(kept), area = fringe.reduce((s, i) => s + m.areaCell[i] * kept.concentration[i], 0);
    growth.push(area);
    if (n === 0) assert.ok(fringe.length >= 20, `${fringe.length} cells beyond the pack after the first call`);
    assert.equal(beyond(melted).length, 0, `the minimum floors at step ${n}`);
  }
  const fringe = beyond(kept), top = Math.max(...fringe.map((i) => kept.concentration[i]));
  console.log(`ice area beyond the pack over six steps of cooling: ${growth.map((a) => (a / m.radius ** 2).toExponential(2)).join(', ')} (unit sphere); the cover reaches ${top.toFixed(4)} over ${fringe.length} cells`);
  assert.ok(growth.every((a, n) => n === 0 || a > growth[n - 1]), 'the area beyond the pack grows at every step');
  assert.ok(growth[5] > 4 * growth[0] && top > 2 * MINIMUM_CONCENTRATION, `${growth}`);
});

function coupled(iceDrift, setup, iceOptions = {}, water = 2) {
  const m = meshOf(16), C = m.nCells, E = m.nEdges;
  const ocean = createOcean(m, { everySteps: 1, thermoclineTilt: 0, salinityProfile: () => 34.5, iceDrift });
  const seaIce = createSeaIce(m, { heatCapacity: ocean.capacity, ...iceOptions });
  ocean.useSeaIce(seaIce.concentration, seaIce);
  const surfaceT = new Float64Array(C).fill(FREEZING_POINT + water), ice = new Float64Array(C), o = { m, C, E, ocean, seaIce, surfaceT, ice, concentration: seaIce.concentration, calm: new Float64Array(E) };
  setup(o);
  ocean.initialize(surfaceT, ice, { concentration: seaIce.concentration });
  for (let i = 0; i < C; i++) if (!(ice[i] > 0)) surfaceT[i] = FREEZING_POINT + water;
  o.env = (dt = 1350) => ({ surfaceT, ice, concentration: seaIce.concentration, stress: o.calm, waterU: ocean.u, eta: ocean.eta, dt });
  o.step = (stress = o.calm, dt = 1350) => ocean.advance(surfaceT, ice, seaIce.oceanFlux, stress, dt, seaIce.concentration);
  return o;
}

test('ice that drifts into water above freezing at a Courant number of 0.005 persists, then melts in exactly the heat the top layer hands it', () => {
  const m = meshOf(16), A = nearestCell(m, -60 * DEG, 0), B = m.cellsOnCell[m.maxEdges * A], edge = m.edgesOnCell[m.maxEdges * A], sign = m.cellsOnEdge[2 * edge] === A ? 1 : -1;
  const run = (keep) => {
    const o = coupled({ rheology: 'frozen' }, (r) => { r.ice[A] = 1; r.concentration[A] = 0.9; r.seaIce.snow[A] = 10; r.surfaceT[A] = FREEZING_POINT - 5; }, { keepSmallIce: keep }), { ocean, seaIce } = o, drift = ocean.iceDrift;
    o.step();
    assert.equal(o.ice[B], 0);
    const water = o.surfaceT[B];
    drift.velocity[edge] = sign * 0.005 * m.dcEdge[edge] / DT;
    drift.transport(o.env(DT));
    return { o, ocean, seaIce, drift, water };
  };
  const control = run(false);
  assert.equal(control.o.ice[B], 0, 'the minimum floors melt the arriving ice at once'); assert.equal(control.drift.arrived[B], 0);

  const { o, ocean, seaIce, drift, water } = run(true);
  assert.ok(o.ice[B] > 0 && o.concentration[B] >= NEGLIGIBLE_CONCENTRATION && o.concentration[B] < MINIMUM_CONCENTRATION, `B holds ${o.concentration[B]}`);
  assert.equal(drift.arrived[B], 1); assert.equal(drift.waterT[B], water);
  const inflow = seaIce.energy(o.surfaceT[B], o.ice[B], seaIce.snow[B], o.concentration[B]), volume = m.areaCell[B] * o.ice[B] * o.concentration[B];
  ocean.readSurface(o.surfaceT, o.ice, o.concentration);
  assert.equal(ocean.T0[B], water); assert.equal(ocean.Q[B], ocean.h[B] * water);
  ocean.writeSurface(o.surfaceT, seaIce.oceanFlux, DT);
  const handed = seaIce.oceanFlux[B] * DT / (ocean.rhoCp * ocean.h[B]), capacity = ocean.capacity[B];
  near(handed, water - FREEZING_POINT, 1e-9, 'the heat the ice base receives');
  const melted = seaIce.budget.melted;
  seaIce.update(o.surfaceT, o.ice, new Float64Array(o.C), B, DT);
  assert.equal(o.ice[B], 0, 'the water\'s heat melts the arrived ice');
  near(seaIce.budget.melted - melted, volume, 1e-12 * volume, 'the ice that melted is the ice that arrived');
  const expected = water + inflow / capacity;
  console.log(`B arrives at ${(water - FREEZING_POINT).toFixed(1)} K above freezing with ${(volume / m.areaCell[B] * 1e3).toFixed(2)} mm of ice and ends ${(o.surfaceT[B] - expected).toExponential(2)} K from the water's heat less the ice's latent heat`);
  near(o.surfaceT[B], expected, 1e-8, 'the water cools by exactly the ice\'s energy');
});

test('a cell below the minimum cover whose thickness, snow or skin temperature has left the bounds melts away with its energy and mass, a consistent one stays', () => {
  const m = meshOf(16), r = rig(m, { rheology: 'frozen' }, null, {}, { keepSmallIce: true }), dt = DT;
  r.dt = dt; r.surfaceT.fill(FREEZING_POINT + 0.5);
  const cases = [['consistent', 1, 10, FREEZING_POINT - 3, true], ['thick', 2 * ICE_REMNANT.thickness, 10, FREEZING_POINT - 3, false], ['snowy', 1, 2 * ICE_REMNANT.snow, FREEZING_POINT - 3, false], ['cold', 1, 10, ICE_REMNANT.skinLow - 50, false], ['warm', 1, 10, ICE_REMNANT.skinHigh + 20, false]];
  const sites = [[-60, 0], [-40, 72], [-20, 144], [0, 216], [20, 288]].map(([lat, lon]) => nearestCell(m, lat * DEG, lon * DEG));
  const pairs = cases.map(([name, h, s, T, ok], n) => {
    const c = sites[n], d = m.cellsOnCell[m.maxEdges * c], e = m.edgesOnCell[m.maxEdges * c], sign = m.cellsOnEdge[2 * e] === c ? 1 : -1;
    r.thickness[c] = h; r.concentration[c] = 0.004; r.snow[c] = s; r.surfaceT[c] = T;
    r.drift.velocity[e] = sign * 0.001 * m.dcEdge[e] / dt;
    return { name, c, d, ok };
  });
  const energy = () => { let s = 0; for (let i = 0; i < r.C; i++) s += m.areaCell[i] * r.seaIce.energy(r.surfaceT[i], r.thickness[i], r.snow[i], r.concentration[i]); return s; };
  const volume = () => { let v = 0; for (let i = 0; i < r.C; i++) v += m.areaCell[i] * r.thickness[i] * r.concentration[i]; return v; };
  const before = energy(), V0 = volume(), melted = r.seaIce.budget.melted;
  r.drift.arrived.fill(0);
  r.drift.transport(r.env());
  let held = 0;
  for (let i = 0; i < r.C; i++) if (r.drift.arrived[i]) held += m.areaCell[i] * r.seaIce.slabHeatCapacity * (r.drift.waterT[i] - FREEZING_POINT);
  for (const { name, c, d, ok } of pairs) {
    assert.equal(r.thickness[c] > 0, ok, `${name}: the source cell`); assert.equal(r.thickness[d] > 0, ok, `${name}: the cell the ice reached`);
    if (!ok) assert.ok(r.concentration[c] === 0 && r.concentration[d] === 0 && r.snow[c] === 0 && r.snow[d] === 0);
    else assert.ok(r.concentration[c] > 0.0039 && r.concentration[d] >= NEGLIGIBLE_CONCENTRATION && r.concentration[d] < 1e-5);
  }
  const gone = r.seaIce.budget.melted - melted;
  near((volume() + gone) / V0, 1, 1e-13, 'volume plus what melted away');
  near(gone, r.drift.count.removed, 1e-12 * gone, 'what the budget counts is what the drift reports');
  near((energy() + held) / before, 1, 1e-12, 'energy, with the water the arrived ice holds off the surface');
  for (let i = 0; i < r.C; i++) if (r.thickness[i] > 0) {
    assert.ok(r.thickness[i] <= ICE_REMNANT.thickness && r.snow[i] <= ICE_REMNANT.snow && r.surfaceT[i] >= ICE_REMNANT.skinLow && r.surfaceT[i] <= ICE_REMNANT.skinHigh, `cell ${i}`);
  }
});

function drain(seed, bounded, { courant = 0.19, calls = 150, centers = 6 } = {}) {
  const saved = { ...ICE_REMNANT };
  if (!bounded) Object.assign(ICE_REMNANT, { thickness: Infinity, snow: Infinity, skinLow: -Infinity, skinHigh: Infinity });
  try {
    let state = seed;
    const rand = () => { state = (state + 0x6D2B79F5) | 0; let t = Math.imul(state ^ (state >>> 15), 1 | state); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
    const m = meshOf(16), r = rig(m, { rheology: 'frozen' }, null, {}, { keepSmallIce: true });
    r.dt = DT;
    for (let i = 0; i < r.C; i++) {
      r.thickness[i] = rand() < 0.5 ? 0.1 + 0.3 * rand() : 2 + 3 * rand(); r.concentration[i] = 0.3 + 0.7 * rand();
      r.snow[i] = rand() < 0.5 ? 2 * rand() : 40 + 40 * rand(); r.surfaceT[i] = rand() < 0.5 ? FREEZING_POINT - rand() : FREEZING_POINT - 25 - 5 * rand();
      r.seaIce.snowAlbedo[i] = rand() < 0.5 ? 0.6 : 0.85;
    }
    const sites = Array.from({ length: centers }, () => Math.floor(rand() * r.C)), closeness = Float64Array.from({ length: r.C }, (_, i) => Math.max(...sites.map((c) => dot(unit(m.xCell, i), unit(m.xCell, c)))));
    for (let e = 0; e < r.E; e++) r.drift.velocity[e] = (closeness[m.cellsOnEdge[2 * e]] > closeness[m.cellsOnEdge[2 * e + 1]] ? 1 : -1) * courant * m.dcEdge[e] / DT;
    const out = { cells: 0, trace: 0, worst: [], nonfinite: 0 };
    for (let k = 0; k < calls; k++) {
      r.drift.arrived.fill(0);
      r.drift.transport(r.env());
      for (let i = 0; i < r.C; i++) {
        if (!(r.thickness[i] > 0)) continue;
        out.cells++;
        if (![r.thickness[i], r.concentration[i], r.snow[i], r.surfaceT[i], r.seaIce.snowAlbedo[i]].every(Number.isFinite)) out.nonfinite++;
        if (r.concentration[i] < MINIMUM_CONCENTRATION) {
          out.trace++;
          if (r.thickness[i] > saved.thickness || r.snow[i] > saved.snow || r.surfaceT[i] < saved.skinLow || r.surfaceT[i] > saved.skinHigh) out.worst.push(i);
        }
        assert.ok(r.seaIce.snowAlbedo[i] >= 0 && r.seaIce.snowAlbedo[i] <= 1, 'snow albedo');
      }
    }
    return out;
  } finally { Object.assign(ICE_REMNANT, saved); }
}

test('under strongly divergent flows no cell below the minimum cover keeps a thickness, snow or skin temperature outside the bounds, which the same flows give it without them', () => {
  let guarded = 0, free = 0, trace = 0, nonfinite = 0;
  for (let seed = 1; seed <= 8; seed++) {
    const on = drain(seed, true), off = drain(seed, false);
    guarded += on.worst.length; free += off.worst.length; trace += on.trace; nonfinite += on.nonfinite;
  }
  console.log(`${trace} cell-calls below 1 percent cover over eight draining runs; outside the bounds: ${guarded} with them, ${free} without`);
  assert.ok(trace > 1000);
  assert.equal(nonfinite, 0);
  assert.equal(guarded, 0);
  assert.ok(free > 0, 'the bounds are doing work');
});

test('over eighty coupled calls with the cavitating fluid and the thermodynamics the ice volume closes with its budget, the drift removes next to none of it, and the stored volume follows the transport', () => {
  const m = meshOf(16), C = m.nCells;
  const patch = (o) => {
    const c = [Math.cos(-70 * DEG), 0, Math.sin(-70 * DEG)];
    for (let i = 0; i < C; i++) if (Math.acos(Math.min(1, dot(unit(m.xCell, i), c))) < 25 * DEG) { o.ice[i] = 0.9; o.concentration[i] = 0.8; o.seaIce.snow[i] = 15; o.surfaceT[i] = FREEZING_POINT - 3; }
  };
  const wind = windOn(m, [0.08, 0, 0]);
  const history = (keep) => {
    const o = coupled({ rheology: 'cavitating' }, patch, { keepSmallIce: keep }, 0.2), { ocean, seaIce } = o, flux = new Float64Array(C), start = { ...seaIce.budget };
    const volume = () => { let v = 0; for (let i = 0; i < C; i++) v += m.areaCell[i] * o.ice[i] * o.concentration[i]; return v; };
    const V0 = volume();
    let removed = 0, dynamic = 0, trace = 0, follows = 0;
    for (let n = 0; n < 80; n++) {
      flux.fill(n < 40 ? -30 : 60);
      for (let i = 0; i < C; i++) seaIce.update(o.surfaceT, o.ice, flux, i, 1350);
      const before = volume();
      o.step(wind);
      dynamic += Math.abs(volume() - before); removed += ocean.iceDrift.count.removed;
      for (let i = 0; i < C; i++) if (o.ice[i] > 0 && o.concentration[i] < MINIMUM_CONCENTRATION) trace++;
    }
    const kept = ocean.restartArrays().previousIceVolume;
    for (let i = 0; i < C; i++) follows = Math.max(follows, Math.abs(kept[i] - o.ice[i] * o.concentration[i]));
    const d = (key) => seaIce.budget[key] - start[key], ledger = d('frozen') + d('snowIce') / 917 - d('melted');
    return { V0, V1: volume(), ledger, removed, dynamic, trace, follows };
  };
  const on = history(true), off = history(false);
  console.log(`80 calls: closure ${((on.V1 - on.V0 - on.ledger) / on.V0).toExponential(2)} of the volume; the drift removes ${(on.removed / on.V0).toExponential(2)} (${(off.removed / off.V0).toExponential(2)} with the minimum floors), moves ${(on.dynamic / on.V0).toExponential(2)} net per call summed; ${on.trace} trace cell-calls`);
  for (const run of [on, off]) near(run.V1 - run.V0, run.ledger, 1e-12 * run.V0, 'the volume against the budget');
  assert.ok(on.removed < 1e-6 * on.V0, `removed ${on.removed / on.V0}`);
  assert.ok(off.removed > 1e-3 * off.V0, `melted on arrival with the minimum floors ${off.removed / off.V0}`);
  assert.ok(on.dynamic < 2 * on.removed + 1e-12 * on.V0, 'the dynamic term of every call is the volume the drift removes and nothing else');
  assert.ok(on.trace > 100 && off.trace === 0);
  assert.ok(on.follows < 1e-12 && off.follows < 1e-12, `previousVolume differs from the measured volume by ${on.follows}`);
});

test('the surface energy changes by exactly the surface flux while ice drifts into slightly warm water, freezes on and melts away', () => {
  const m = meshOf(16), r = pack(true, { water: 0.3, edge: -60 }), flux = new Float64Array(r.C);
  const energy = () => { let s = 0; for (let i = 0; i < r.C; i++) s += m.areaCell[i] * r.seaIce.energy(r.surfaceT[i], r.thickness[i], r.snow[i], r.concentration[i]); return s; };
  const scale = m.areaCell.reduce((s, a, i) => s + (r.initial[i] ? a * r.seaIce.latent : 0), 0);
  let expected = energy(), held = 0, worst = 0, arrivals = 0, trace = 0;
  for (let n = 0; n < 60; n++) {
    flux.fill(n < 20 ? -30 : 80);
    for (let i = 0; i < r.C; i++) r.seaIce.update(r.surfaceT, r.thickness, flux, i, DT);
    expected += DT * m.areaCell.reduce((s, a, i) => s + a * flux[i], 0);
    r.carry();
    for (let i = 0; i < r.C; i++) {
      if (r.drift.arrived[i]) { held += m.areaCell[i] * r.seaIce.slabHeatCapacity * (r.drift.waterT[i] - FREEZING_POINT); arrivals++; }
      if (r.thickness[i] > 0 && r.concentration[i] < MINIMUM_CONCENTRATION) trace++;
    }
    worst = Math.max(worst, Math.abs(energy() + held - expected) / scale);
  }
  console.log(`60 steps with ${arrivals} arrivals in water 0.3 K above freezing and ${trace} trace cell-steps: the energy follows the flux to ${worst.toExponential(2)} of the pack's latent heat`);
  assert.ok(arrivals >= 20 && trace >= 20 && worst < 1e-11, `${arrivals} ${trace} ${worst}`);
});
