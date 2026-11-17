import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createSeaIce as defaultSeaIce, openWaterAlbedo, iceExchangeOptions, ICE_EXCHANGE, DEFAULT_ICE_EXCHANGE, FREEZING_POINT, MELTING_POINT } from '../js/physics/ice.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { PREVIOUS_SEA_ICE } from './helpers/previousDefaults.mjs';

const createSeaIce = (mesh, options = {}) => defaultSeaIce(mesh, { ...PREVIOUS_SEA_ICE, ...options });
const model = createModel(new Grid(16));
const seaIce = createSeaIce(model.mesh);
const energyOf = (sea, surfaceT, ice, i = 0) => sea.energy(surfaceT[i], ice[i], sea.snow[i], sea.cover(i, ice[i]));

test('the surface energy changes by exactly the surface flux through freezing, growth, melting and thaw', () => {
  const surfaceT = new Float64Array([FREEZING_POINT + 2]), ice = new Float64Array([0]);
  const flux = new Float64Array(1);
  const dt = 900;
  let energy = energyOf(seaIce, surfaceT, ice);
  const scale = seaIce.slabHeatCapacity * 2;
  let froze = false, melted = false;
  for (let n = 0; n < 4 * 96 * 60; n++) {
    flux[0] = n < 96 * 90 ? -150 : 250;
    seaIce.update(surfaceT, ice, flux, 0, dt);
    energy += dt * flux[0];
    assert.ok(Math.abs(energyOf(seaIce, surfaceT, ice) - energy) < 1e-9 * scale, `step ${n}: energy ${energyOf(seaIce, surfaceT, ice)} vs ${energy}`);
    if (ice[0] > 0) { froze = true; assert.ok(surfaceT[0] <= MELTING_POINT + 1e-9 && surfaceT[0] < FREEZING_POINT + 20); }
    if (froze && ice[0] === 0) melted = true;
  }
  assert.ok(froze && melted, 'expected the cell to freeze over and later thaw');
  assert.ok(surfaceT[0] > FREEZING_POINT, 'open water ends above freezing');
});

test('a cold skin over thick ice conducts heat up from the base and grows the ice', () => {
  const surfaceT = new Float64Array([FREEZING_POINT - 20]), ice = new Float64Array([1]);
  const flux = new Float64Array([0]);
  seaIce.update(surfaceT, ice, flux, 0, 900);
  assert.ok(ice[0] > 1, 'ice grew');
  assert.ok(surfaceT[0] > FREEZING_POINT - 20, 'the skin warmed');
});

test('open water is dark under a high sun and bright near the horizon', () => {
  assert.ok(openWaterAlbedo(1) < 0.03);
  assert.ok(Math.abs(openWaterAlbedo(0.5) - 0.07) < 0.01);
  assert.ok(openWaterAlbedo(0.1) > 0.25 && openWaterAlbedo(0.1) < 0.4);
  for (let mu = 0.05; mu < 0.8; mu += 0.05) assert.ok(openWaterAlbedo(mu) > openWaterAlbedo(mu + 0.05));
  assert.ok(openWaterAlbedo(1) < openWaterAlbedo(0.5));
  assert.ok(seaIce.albedo(0, 0.1) > seaIce.albedo(0, 0.9));
  assert.ok(Math.abs(seaIce.albedo(0) - 0.06) < 1e-12, 'diffuse light sees the diffuse albedo');
  assert.equal(seaIce.albedo(1, 0.1), seaIce.albedo(1, 0.9));
});

test('albedo rises from open water to thick ice', () => {
  assert.ok(seaIce.albedo(0) < 0.1);
  assert.ok(seaIce.albedo(0.25) > seaIce.albedo(0) && seaIce.albedo(0.25) < seaIce.albedo(1));
  assert.equal(seaIce.albedo(1), seaIce.albedo(3));
});

test('the initial state has ice at the poles and the model steps with it', () => {
  const m = createModel(new Grid(16));
  const init = initializeState(m, {});
  for (let a = 0; a < init.length; a++) m.state[a].set(init[a]);
  const d0 = m.diagnostics();
  assert.ok(d0.iceFraction > 0.02 && d0.iceFraction < 0.3, `initial ice fraction ${d0.iceFraction}`);
  for (let n = 0; n < 32; n++) m.step(1350);
  const d = m.diagnostics();
  assert.ok(Number.isFinite(d.maxWind) && d.maxWind < 80);
  assert.ok(d.planetaryAlbedo > 0.05 && d.planetaryAlbedo < 0.6, `planetary albedo ${d.planetaryAlbedo}`);
  assert.ok(d.iceFraction > 0);
  console.log(`N=16 after 12 h: ice fraction ${d.iceFraction.toFixed(3)}, mean thickness ${d.iceThickness.toFixed(2)} m, surface albedo ${d.surfaceAlbedo.toFixed(3)}, planetary albedo ${d.planetaryAlbedo.toFixed(3)}, TCW ${(1000 * d.columnCloud).toFixed(0)} g/m²`);
});

test('snow on the ice brightens it, insulates it, melts before it, and goes into the water when the ice is gone', () => {
  const snowy = createSeaIce(model.mesh), bare = createSeaIce(model.mesh);
  const flux = (w) => new Float64Array([w]);
  assert.ok(Math.abs(snowy.albedo(1, null, 0) - 0.62) < 1e-12 && Math.abs(snowy.albedo(1, null, 20) - 0.85) < 1e-12 && Math.abs(snowy.albedo(1, null, 10) - 0.735) < 1e-12);
  assert.equal(snowy.albedo(0, null, 20), snowy.albedo(0), 'open water is not brightened');
  const ice = new Float64Array([1, 0]), water = new Float64Array([FREEZING_POINT - 3, FREEZING_POINT + 1]);
  assert.equal(snowy.deposit(0, 5, MELTING_POINT - 5, ice, water), true);
  assert.equal(snowy.snow[0], 5);
  assert.equal(snowy.deposit(0, 5, MELTING_POINT + 5, ice, water), false, 'rain on ice is not snow');
  assert.equal(snowy.deposit(1, 5, MELTING_POINT - 5, ice, water), true, 'snow on open water melts at once');
  assert.equal(snowy.snow[1], 0);
  assert.ok(Math.abs((FREEZING_POINT + 1 - water[1]) - 5 * snowy.latentHeatFusion / snowy.slabHeatCapacity) < 1e-12, 'and cools the water by its latent heat');
  assert.equal(snowy.snow[0], 5);
  snowy.snow[1] = 4; water[1] = FREEZING_POINT + 1;
  snowy.update(water, ice, new Float64Array([0, 0]), 1, 900);
  assert.equal(snowy.snow[1], 0, 'stray snow on open water melts in the update');
  assert.ok(Math.abs((FREEZING_POINT + 1 - water[1]) - 4 * snowy.latentHeatFusion / snowy.slabHeatCapacity) < 1e-12);

  const iceB = new Float64Array([1]), tB = new Float64Array([FREEZING_POINT - 20]);
  bare.update(tB, iceB, flux(0), 0, 900);
  snowy.snow[0] = 30;
  const iceS = new Float64Array([1]), tS = new Float64Array([FREEZING_POINT - 20]);
  snowy.update(tS, iceS, flux(0), 0, 900);
  assert.ok(iceS[0] > 1 && iceS[0] - 1 < 0.7 * (iceB[0] - 1) && iceS[0] - 1 > 0.5 * (iceB[0] - 1), `snow slows the growth: ${iceS[0] - 1} against ${iceB[0] - 1}`);

  snowy.snow[0] = 10;
  const iceM = new Float64Array([1]), tM = new Float64Array([MELTING_POINT]);
  let steps = 0;
  while (snowy.snow[0] > 0 && steps++ < 1000) snowy.update(tM, iceM, flux(300), 0, 900);
  assert.ok(steps < 1000 && iceM[0] > 0.999, `the ice waits for the snow: ${iceM[0]} m left after ${steps} steps`);

  const surfaceT = new Float64Array([FREEZING_POINT - 5]), thin = new Float64Array([0.3]);
  snowy.snow[0] = 0; snowy.concentration[0] = 1;
  let energy = energyOf(snowy, surfaceT, thin);
  const scale = snowy.slabHeatCapacity * 2, dt = 900;
  for (let n = 0; n < 96 * 40; n++) {
    if (n % 96 === 0 && n < 96 * 10) { snowy.deposit(0, 2, MELTING_POINT - 5, thin); energy -= snowy.latentHeatFusion * 2; }
    const w = n < 96 * 10 ? -50 : 300;
    snowy.update(surfaceT, thin, flux(w), 0, dt);
    energy += dt * w;
    assert.ok(Math.abs(energyOf(snowy, surfaceT, thin) - energy) < 1e-9 * scale, `step ${n}: energy ${energyOf(snowy, surfaceT, thin)} vs ${energy}`);
  }
  assert.equal(thin[0], 0, 'the ice melted away');
  assert.equal(snowy.snow[0], 0, 'and took its snow into the water');

  const tA = new Float64Array([MELTING_POINT]), iceA = new Float64Array([0.05]), tC = new Float64Array([MELTING_POINT]), iceC = new Float64Array([0.05]);
  snowy.snow[0] = 10;
  for (let n = 0; n < 96; n++) { snowy.update(tA, iceA, flux(400), 0, dt); bare.update(tC, iceC, flux(400), 0, dt); }
  assert.ok(iceA[0] === 0 && iceC[0] === 0);
  assert.ok(Math.abs((tC[0] - tA[0]) - 10 * snowy.latentHeatFusion / snowy.slabHeatCapacity) < 1e-9, `the snow's latent heat came out of the water: ${tC[0] - tA[0]} K`);
});

test('snow heavier than the freeboard floods into snow-ice, restoring the freeboard and conserving energy', () => {
  const sea = createSeaIce(model.mesh), flux = new Float64Array([0]), freeboard = 1026 - 917;
  const ice = new Float64Array([0.5]), t = new Float64Array([FREEZING_POINT - 10]);
  sea.snow[0] = 20;
  sea.update(t, ice, flux, 0, 900);
  assert.equal(sea.snow[0], 20, 'light snow rides above the water line');
  assert.equal(sea.budget.snowIce, 0);
  sea.snow[0] = 200;
  const before = sea.energy(t[0], ice[0], sea.snow[0]), mass = 917 * ice[0] + sea.snow[0], grown = ice[0];
  sea.update(t, ice, flux, 0, 900);
  assert.ok(sea.snow[0] < 200 && sea.snow[0] > 0 && ice[0] > grown + 0.05, `snow ${sea.snow[0]} kg/m² on ${ice[0]} m of ice`);
  assert.ok(Math.abs(freeboard * ice[0] - sea.snow[0]) < 1e-9, `the flooded snow brings the ice back to the water line: ${freeboard * ice[0] - sea.snow[0]}`);
  assert.ok(Math.abs(sea.energy(t[0], ice[0], sea.snow[0]) - before) < 1e-6, 'no energy is made or lost in the conversion');
  assert.ok(sea.budget.snowIce > 0 && Math.abs(917 * ice[0] + sea.snow[0] - mass - 917 * (ice[0] - grown - sea.budget.snowIce / model.mesh.areaCell[0] / 917)) < 1e-9, 'mass moves from the snow to the ice');
});

test('partly covered ice conserves energy exactly through partial melt, lead freezing, thaw and refreezing', () => {
  const sea = createSeaIce(model.mesh);
  const surfaceT = new Float64Array([FREEZING_POINT - 5]), ice = new Float64Array([0.6]), flux = new Float64Array(1);
  sea.concentration[0] = 0.5; sea.snow[0] = 5; sea.oceanFlux[0] = 4;
  const dt = 900, day = 96, scale = sea.slabHeatCapacity * 2;
  let energy = energyOf(sea, surfaceT, ice), previous = sea.concentration[0];
  const seen = { closing: false, opening: false, thawed: false, refrozen: false, partialRefreeze: false };
  for (let n = 0; n < 200 * day; n++) {
    const phase = n < 40 * day ? 'winter' : n < 120 * day ? 'summer' : 'autumn';
    flux[0] = phase === 'summer' ? 160 : phase === 'winter' ? -120 : -200;
    const contrast = phase === 'summer' ? 120 : 0;
    if (phase === 'winter' && n % day === 0 && sea.deposit(0, 1.5, MELTING_POINT - 10, ice, surfaceT)) energy -= sea.latentHeatFusion * 1.5;
    sea.update(surfaceT, ice, flux, 0, dt, contrast);
    energy += dt * (flux[0] + sea.oceanFlux[0]);
    assert.ok(Math.abs(energyOf(sea, surfaceT, ice) - energy) < 1e-9 * scale, `step ${n} (${phase}): energy ${energyOf(sea, surfaceT, ice)} vs ${energy}`);
    const A = sea.concentration[0];
    assert.ok((ice[0] > 0) === (A > 0) && A <= 1, `step ${n}: ${ice[0]} m of ice over ${A} of the cell`);
    if (phase === 'winter' && A > previous) seen.closing = true;
    if (phase === 'summer' && A > 0 && A < previous) seen.opening = true;
    if (phase === 'summer' && previous > 0 && A === 0) seen.thawed = true;
    if (phase === 'autumn' && ice[0] > 0) { seen.refrozen = true; if (A < 1) seen.partialRefreeze = true; }
    previous = A;
  }
  assert.deepEqual(seen, { closing: true, opening: true, thawed: true, refrozen: true, partialRefreeze: true });
  assert.ok(sea.budget.leadFrozen > 0 && sea.budget.lateralMelted > 0, `leads froze ${sea.budget.leadFrozen} and melted ${sea.budget.lateralMelted}`);
});

test('a half-covered cell in the sun loses its area faster than a full one loses thickness, and is gone sooner', () => {
  const half = createSeaIce(model.mesh), full = createSeaIce(model.mesh), dt = 900;
  half.concentration[0] = 0.5; full.concentration[0] = 1;
  const tHalf = new Float64Array([MELTING_POINT]), hHalf = new Float64Array([1]), tFull = new Float64Array([MELTING_POINT]), hFull = new Float64Array([0.5]);
  const onIce = 200, contrast = 340 * (0.5 - 0.06);
  const fluxHalf = new Float64Array([onIce + 0.5 * contrast]), fluxFull = new Float64Array([onIce]);
  for (let n = 0; n < 96; n++) { half.update(tHalf, hHalf, fluxHalf, 0, dt, contrast); full.update(tFull, hFull, fluxFull, 0, dt, contrast); }
  const areaLost = 1 - half.concentration[0] / 0.5, thicknessLost = 1 - hFull[0] / 0.5;
  console.log(`after a day in the sun: the half-covered cell lost ${(100 * areaLost).toFixed(1)}% of its area, the full cell ${(100 * thicknessLost).toFixed(1)}% of its thickness and ${(100 * (1 - full.concentration[0])).toFixed(1)}% of its area`);
  assert.ok(areaLost > thicknessLost && areaLost > 1 - full.concentration[0], `area ${areaLost} against thickness ${thicknessLost}`);
  let stepsHalf = 96, stepsFull = 96;
  while (hHalf[0] > 0 && stepsHalf < 96 * 100) { half.update(tHalf, hHalf, fluxHalf, 0, dt, contrast); stepsHalf++; }
  while (hFull[0] > 0 && stepsFull < 96 * 100) { full.update(tFull, hFull, fluxFull, 0, dt, contrast); stepsFull++; }
  console.log(`the half-covered cell melted away after ${(stepsHalf / 96).toFixed(1)} days, the full cell of the same volume after ${(stepsFull / 96).toFixed(1)}`);
  assert.ok(half.concentration[0] === 0 && full.concentration[0] === 0 && stepsHalf < stepsFull);
});

test('areaLossFactor scales the area melting takes, conserves energy, and 1 is Hibler\'s rule', () => {
  const run = (options) => {
    const sea = createSeaIce(model.mesh, options), surfaceT = new Float64Array([MELTING_POINT]), ice = new Float64Array([1]), flux = new Float64Array([200]);
    sea.concentration[0] = 0.8;
    let energy = energyOf(sea, surfaceT, ice);
    for (let n = 0; n < 96; n++) {
      sea.update(surfaceT, ice, flux, 0, 900, 150);
      energy += 900 * flux[0];
      assert.ok(Math.abs(energyOf(sea, surfaceT, ice) - energy) < 1e-9 * sea.slabHeatCapacity, `step ${n}: energy ${energyOf(sea, surfaceT, ice)} vs ${energy}`);
    }
    return { area: sea.concentration[0], volume: sea.concentration[0] * ice[0] };
  };
  const hibler = run({}), one = run({ areaLossFactor: 1 }), half = run({ areaLossFactor: 0.5 }), none = run({ areaLossFactor: 0 });
  assert.deepEqual(one, hibler);
  console.log(`a day melting from 0.8 cover of 1 m: area ${hibler.area.toFixed(4)} (Hibler), ${half.area.toFixed(4)} (0.5), ${none.area.toFixed(4)} (0); volume ${hibler.volume.toFixed(4)}, ${half.volume.toFixed(4)}, ${none.volume.toFixed(4)}`);
  assert.ok(hibler.area < half.area && half.area < none.area && none.area === 0.8, `areas ${hibler.area}, ${half.area}, ${none.area}`);
  const lost = (r) => 0.8 - r.area;
  assert.ok(Math.abs(lost(half) / lost(hibler) - 0.5) < 0.1, `half the factor loses ${lost(half) / lost(hibler)} of Hibler's area`);
});

test("with areaLossFactor a melting cell gives up that share of Hibler's opening in one step and thins by the rest, and its energy closes over a day at 0.5, 0.25 and 0", () => {
  const dt = 900, contrast = 340 * (0.5 - 0.06), scale = seaIce.slabHeatCapacity * 2;
  const run = (options, steps) => {
    const sea = createSeaIce(model.mesh, options);
    sea.concentration[0] = 0.5;
    const surfaceT = new Float64Array([MELTING_POINT]), ice = new Float64Array([1]), flux = new Float64Array([200 + 0.5 * contrast]);
    let energy = energyOf(sea, surfaceT, ice), worst = 0;
    for (let n = 0; n < steps; n++) {
      sea.update(surfaceT, ice, flux, 0, dt, contrast);
      energy += dt * (flux[0] + sea.oceanFlux[0]);
      worst = Math.max(worst, Math.abs(energyOf(sea, surfaceT, ice) - energy));
    }
    return { A: sea.concentration[0], h: ice[0], worst };
  };
  const hibler = run({}, 1), half = run({ areaLossFactor: 0.5 }, 1), none = run({ areaLossFactor: 0 }, 1);
  console.log(`one step in the sun from cover 0.5: Hibler's opening ${(0.5 - hibler.A).toExponential(4)}, areaLossFactor 0.5 ${(0.5 - half.A).toExponential(4)}, 0 ${(0.5 - none.A).toExponential(4)}; thickness ${hibler.h.toFixed(6)}, ${half.h.toFixed(6)}, ${none.h.toFixed(6)} m`);
  assert.ok(hibler.A < 0.5 && Math.abs(0.5 - half.A - 0.5 * (0.5 - hibler.A)) < 1e-15, `opening ${0.5 - half.A} against half of ${0.5 - hibler.A}`);
  assert.equal(none.A, 0.5);
  assert.ok(none.h < half.h && half.h < hibler.h, 'the melt that opens no water thins the ice');
  for (const areaLossFactor of [0.5, 0.25, 0]) {
    const day = run({ areaLossFactor }, 96);
    assert.ok(day.worst < 1e-9 * scale, `areaLossFactor ${areaLossFactor}: energy off by ${day.worst} J/m²`);
  }
});

test('a lead under a cold sky closes: the area rises toward full cover as the volume grows by what freezes', () => {
  const sea = createSeaIce(model.mesh, { leadExchange: 0 }), dt = 900;
  const surfaceT = new Float64Array([FREEZING_POINT - 10]), ice = new Float64Array([1]), flux = new Float64Array([-100]);
  sea.concentration[0] = 0.4;
  const volume = () => sea.concentration[0] * ice[0], start = volume(), area = model.mesh.areaCell[0];
  let previous = sea.concentration[0];
  for (let n = 0; n < 96 * 60; n++) {
    sea.update(surfaceT, ice, flux, 0, dt);
    assert.ok(sea.concentration[0] >= previous, `step ${n}: the lead reopened, ${sea.concentration[0]} after ${previous}`);
    previous = sea.concentration[0];
  }
  const closing = 100 * dt / (sea.latent * sea.leadClosing), expected = 1 - 1 / (1 / 0.6 + closing * 96 * 60);
  console.log(`sixty days at −100 W/m²: concentration 0.4 → ${sea.concentration[0].toFixed(3)} (${expected.toFixed(3)} for leads that close as their open area squared), volume ${start.toFixed(3)} → ${volume().toFixed(3)} m, ${(sea.budget.leadFrozen / area).toFixed(3)} m of it frozen in the leads`);
  assert.ok(sea.concentration[0] > 0.85 && Math.abs(sea.concentration[0] - expected) < 0.02, `concentration ${sea.concentration[0]} against ${expected}`);
  assert.ok(Math.abs(volume() - start - sea.budget.frozen / area) < 1e-9, 'the volume grows by exactly what froze');
  assert.ok(sea.budget.leadFrozen > 0 && sea.budget.leadFrozen < sea.budget.frozen, 'part of it in the leads, the rest under the ice');
});

test("with leadClosure 'hibler' a lead under a cold sky closes as its open area, Hibler's h0 law, and the volume still grows by what freezes", () => {
  const sea = createSeaIce(model.mesh, { leadExchange: 0, leadClosure: 'hibler' }), dt = 900, steps = 96 * 20;
  const surfaceT = new Float64Array([FREEZING_POINT - 10]), ice = new Float64Array([1]), flux = new Float64Array([-100]);
  sea.concentration[0] = 0.4;
  const volume = () => sea.concentration[0] * ice[0], start = volume(), area = model.mesh.areaCell[0];
  for (let n = 0; n < steps; n++) sea.update(surfaceT, ice, flux, 0, dt);
  const closing = 100 * dt / (sea.latent * sea.leadClosing), expected = 1 - 0.6 * (1 - closing) ** steps;
  console.log(`twenty days at −100 W/m²: concentration 0.4 → ${sea.concentration[0].toFixed(4)} (${expected.toFixed(4)} for leads that close as their open area)`);
  assert.ok(Math.abs(sea.concentration[0] - expected) < 1e-6, `concentration ${sea.concentration[0]} against ${expected}`);
  assert.ok(Math.abs(volume() - start - sea.budget.frozen / area) < 1e-9, 'the volume grows by exactly what froze');
  assert.throws(() => createSeaIce(model.mesh, { leadClosure: 'linear' }), /leadClosure/);
});

test('a state saved without concentration loads as full cover wherever it has ice', () => {
  const sea = createSeaIce(model.mesh), C = model.mesh.nCells;
  const ice = Float64Array.from({ length: C }, (_, i) => (i % 3 === 0 ? 0 : 0.5 + i % 2));
  sea.load(ice);
  for (let i = 0; i < C; i++) assert.equal(sea.concentration[i], ice[i] > 0 ? 1 : 0);
  const saved = Float64Array.from({ length: C }, (_, i) => (i % 4) / 3);
  sea.load(ice, saved);
  for (let i = 0; i < C; i++) assert.equal(sea.concentration[i], ice[i] > 0 ? (saved[i] > 0 ? Math.min(1, saved[i]) : 1) : 0);
  assert.equal(sea.cover(0, 0), 0);
  sea.concentration[1] = 0;
  assert.equal(sea.cover(1, 0.5), 1, 'ice without a concentration covers its cell');
});

test('the albedo blends linearly in the concentration between open water and full cover', () => {
  const sea = createSeaIce(model.mesh);
  const before = (h, mu, snow) => {
    const water = mu === null ? 0.06 : openWaterAlbedo(mu);
    if (h <= 0) return water;
    const bare = water + (0.62 - water) * Math.min(1, h / 0.5);
    return bare + (0.85 - bare) * Math.min(1, snow / 20);
  };
  for (const mu of [null, 0.2, 0.9]) for (const [h, snow] of [[0.1, 0], [0.3, 8], [2, 30]]) {
    assert.equal(sea.albedo(h, mu, snow, 1), before(h, mu, snow), 'full cover is the ice');
    assert.equal(sea.albedo(h, mu, snow), before(h, mu, snow), 'ice covers its cell unless told otherwise');
    assert.equal(sea.albedo(h, mu, snow, 0), before(0, mu, snow), 'no cover is open water');
    assert.equal(sea.albedoContrast(h, mu, snow), before(h, mu, snow) - before(0, mu, snow));
    for (const A of [0.1, 0.5, 0.85]) assert.ok(Math.abs(sea.albedo(h, mu, snow, A) - (A * before(h, mu, snow) + (1 - A) * before(0, mu, snow))) < 1e-15);
  }
});

test('the leads at the freezing point lose the heat the colder ice would otherwise lose, and freeze', () => {
  const run = (leadExchange) => {
    const sea = createSeaIce(model.mesh, { leadExchange });
    sea.concentration[0] = 0.5;
    const surfaceT = new Float64Array([FREEZING_POINT - 20]), ice = new Float64Array([1]), flux = new Float64Array([0]);
    const before = energyOf(sea, surfaceT, ice), volume = 0.5;
    for (let n = 0; n < 96; n++) sea.update(surfaceT, ice, flux, 0, 900);
    const area = model.mesh.areaCell[0];
    return { sea, skin: surfaceT[0], grown: sea.concentration[0] * ice[0] - volume, leads: sea.budget.leadFrozen / area, energy: energyOf(sea, surfaceT, ice) - before };
  };
  const split = run(10), shared = run(0);
  console.log(`a day at zero net flux, half covered, skin 20 K below freezing: with the split the concentration rose to ${split.sea.concentration[0].toFixed(4)} with ${split.leads.toExponential(2)} m frozen in the leads and ${(split.grown - split.leads).toExponential(2)} m under the ice (skin ${(split.skin - FREEZING_POINT).toFixed(2)} K); without it ${shared.sea.concentration[0]} and ${shared.grown.toExponential(2)} m under the ice (skin ${(shared.skin - FREEZING_POINT).toFixed(2)} K)`);
  assert.ok(split.sea.concentration[0] > 0.5 && split.leads > 0, 'the leads freeze and close');
  assert.equal(shared.sea.concentration[0], 0.5);
  assert.equal(shared.leads, 0, 'without the split the leads share the cell\'s zero flux');
  assert.ok(Math.abs(split.grown - split.leads - (split.sea.budget.frozen / model.mesh.areaCell[0] - split.leads)) < 1e-12, 'the volume grows by what froze');
  assert.ok(split.skin > shared.skin && split.grown - split.leads < shared.grown, 'the ice, relieved of the leads\' loss, warms and grows less at its base');
  for (const { energy } of [split, shared]) assert.ok(Math.abs(energy) < 1e-9 * 2 * split.sea.slabHeatCapacity, `energy moved by ${energy} J/m² at zero flux`);
});

test('with a transmission, open water keeps the passing share of its absorbed sunlight out of its temperature and adds it to the pending store, the two changing by exactly the surface and ocean fluxes; ice and its leads pass none', () => {
  const capacity = new Float64Array([4.1e6 * 20, 4.1e6 * 50, 4.1e6 * 50]), transmission = new Float64Array([0.102, 0.0174, 0.0174]), penetrated = new Float64Array([5, 7, 11]);
  const sea = createSeaIce(model.mesh, { heatCapacity: capacity, transmission, penetrated });
  sea.load(new Float64Array([0, 1, 0.5]), new Float64Array([0, 1, 0.6]));
  sea.oceanFlux.set([0, 3, 2]);
  const surfaceT = new Float64Array([FREEZING_POINT + 12, FREEZING_POINT - 5, FREEZING_POINT - 1]), ice = new Float64Array([0, 1, 0.5]);
  const flux = new Float64Array([180, 120, 140]), solar = 640, dt = 900;
  const before = Float64Array.from(penetrated), water = surfaceT[0];
  sea.update(surfaceT, ice, flux, 0, dt, 0, solar);
  const pen = solar * transmission[0], stored = capacity[0] * (surfaceT[0] - water) + (penetrated[0] - before[0]);
  assert.equal(penetrated[0] - before[0], dt * pen);
  assert.ok(Math.abs(stored - dt * flux[0]) < 4 * Number.EPSILON * capacity[0] * surfaceT[0], `open water: capacity·ΔT + Δpenetrated = ${stored} J/m², the fluxes ${dt * flux[0]}, to the rounding of the temperature`);
  assert.ok(capacity[0] * (surfaceT[0] - water) < dt * flux[0] - 0.9 * dt * pen, 'the passing light does not warm the mixed layer');
  for (const i of [1, 2]) {
    sea.update(surfaceT, ice, flux, i, dt, i === 2 ? 30 : 0, solar);
    assert.equal(penetrated[i], before[i], `cell ${i}, ${i === 1 ? 'full ice' : 'ice with leads'}, passes no light`);
  }
  const plain = createSeaIce(model.mesh, { heatCapacity: capacity }), alone = new Float64Array([water]);
  plain.update(alone, new Float64Array(1), flux, 0, dt, 0, solar);
  assert.equal(alone[0], water + dt * flux[0] / capacity[0], 'without a transmission all of the sunlight stays in the mixed layer');
});

for (const brineReservoir of [false, true]) test(`with iceLight${brineReservoir ? ' and brineReservoir' : ''} the surface energy${brineReservoir ? ', the brine counted,' : ''} changes by exactly the surface flux through a sunlit melt and a cold autumn${brineReservoir ? ', the brine ending each step at most half the latent heat of the ice, reaching it, and drained in the cold' : ''}`, () => {
  const sea = createSeaIce(model.mesh, { iceLight: true, brineReservoir });
  const surfaceT = new Float64Array([FREEZING_POINT - 1]), ice = new Float64Array([1]), flux = new Float64Array(1), dt = 900;
  sea.load(ice);
  const total = () => sea.energy(surfaceT[0], ice[0], sea.snow[0], sea.cover(0, ice[0]), brineReservoir ? sea.brine[0] : 0);
  let energy = total(), capped = false, cold = 0, spent = false;
  for (let n = 0; n < 96 * 90; n++) {
    const sunny = n < 96 * 23, solar = sunny ? 400 : 0;
    flux[0] = sunny ? 120 : -120;
    const Q0 = brineReservoir ? sea.brine[0] : 0, T0 = surfaceT[0];
    sea.update(surfaceT, ice, flux, 0, dt, 0, solar, 0.5);
    energy += dt * flux[0];
    assert.ok(Math.abs(total() - energy) < 1e-9 * sea.latent * 2, `step ${n}: energy ${total()} vs ${energy}`);
    assert.ok(ice[0] > 0, `step ${n}: the ice is still there`);
    if (!brineReservoir) continue;
    const Q = sea.brine[0], cap = 0.5 * sea.latent * ice[0];
    assert.ok(Q >= 0 && Q <= cap * (1 + 1e-12), `step ${n}: brine ${Q} over ${cap} for ${ice[0]} m`);
    if (Q >= cap * (1 - 1e-9)) capped = true;
    if (T0 < Math.fround(MELTING_POINT) && Q0 === 0) assert.equal(Q, 0, `step ${n}: an empty brine stays empty under a skin below the melting point`);
    if (!sunny) {
      assert.ok(Q <= Q0, `step ${n}: the brine only drains in the cold`);
      if (Q > 0) cold++;
      if (Q0 > 0 && Q === 0) spent = true;
    }
  }
  if (brineReservoir) assert.ok(capped && cold > 0 && spent, `capped ${capped}, ${cold} cold steps with brine, spent ${spent}`);
});

test('in autumn the brine holds an upper half-layer of the ice at the melting point: the skin takes its conduction from there, the base the downward conduction below it, and the brine pays the difference', () => {
  const k = 2, rs = 300 * 0.31, dt = 1350;
  const build = (brineReservoir) => {
    const sea = createSeaIce(model.mesh, { iceLight: true, brineReservoir, conductivity: k, minimumThickness: 0.1, snowDensity: 300, snowConductivity: 0.31 });
    const surfaceT = new Float64Array([MELTING_POINT]), ice = new Float64Array([1.5]), flux = new Float64Array([-90]);
    sea.load(ice, new Float64Array([1]));
    sea.snow[0] = 10;
    if (brineReservoir) sea.brine[0] = 3e7;
    return { sea, surfaceT, ice, flux };
  };
  const held = build(true), free = build(false);
  let steps = 0, drained = false;
  for (let n = 0; n < 64 * 40 && !drained; n++) {
    const { sea, surfaceT, ice, flux } = held, T = surfaceT[0], h = ice[0], s = sea.snow[0], Q = sea.brine[0], A = sea.cover(0, h), volume = sea.concentration[0] * h;
    sea.update(surfaceT, ice, flux, 0, dt);
    free.sea.update(free.surfaceT, free.ice, free.flux, 0, dt);
    if (!(Q > 0) || T >= MELTING_POINT) continue;
    const upper = (MELTING_POINT - T) / (h / (4 * k) + s / rs), lower = (FREEZING_POINT - MELTING_POINT) / (3 * h / (4 * k)), cost = dt * (upper - lower);
    if (cost < Q) {
      assert.ok(Math.abs(Q - sea.brine[0] - cost) < 1e-9 * Q, `step ${n}: the brine pays ${Q - sea.brine[0]} for ${cost}`);
      const iceFlux = flux[0] + (1 - A) * sea.leadExchange * (FREEZING_POINT - T);
      assert.ok(Math.abs(surfaceT[0] - T - dt * (iceFlux + upper) / sea.skinHeatCapacity) < 1e-9, `step ${n}: the skin takes ${upper} W/m²`);
      const leadIce = Math.max(0, -(1 - A) * (flux[0] - A * sea.leadExchange * (FREEZING_POINT - T)) * dt / sea.latent);
      assert.ok(sea.concentration[0] * ice[0] <= volume + leadIce * (1 + 1e-9), `step ${n}: the ice grows only in its leads while the brine lasts`);
      steps++;
    } else { assert.equal(sea.brine[0], 0, `step ${n}: the last of the brine`); drained = true; }
  }
  console.log(`the brine held the upper half-layer for ${(steps / 64).toFixed(1)} days, the skin at ${held.surfaceT[0].toFixed(2)} K against ${free.surfaceT[0].toFixed(2)} K without it, the ice ${held.ice[0].toFixed(3)} m against ${free.ice[0].toFixed(3)} m`);
  assert.ok(steps > 64 && drained, `${steps} steps`);
  assert.ok(held.surfaceT[0] > free.surfaceT[0] && held.ice[0] < free.ice[0]);
});

test('over ten years a perennial column with the brine reservoir ends thicker than with iceLight alone, as FMM97 find for latent storage', () => {
  const sigma = 5.67e-8, dt = 1350, perDay = 64;
  const column = (brineReservoir) => {
    const sea = createSeaIce(model.mesh, { iceLight: true, brineReservoir }), surfaceT = new Float64Array([250]), ice = new Float64Array([2.5]), flux = new Float64Array(1);
    sea.load(ice);
    sea.oceanFlux[0] = 2;
    let low = Infinity, mean = 0;
    for (let y = 0; y < 10; y++) for (let n = 0; n < 365 * perDay; n++) {
      const d = n / perDay, air = Math.min(273.2, 255 + 18 * Math.cos(2 * Math.PI * (d - 200) / 365)), sw = 0.6 * Math.max(0, 340 * Math.cos(2 * Math.PI * (d - 172) / 365) - 40);
      const solar = sw * (1 - sea.albedo(ice[0], null, sea.snow[0], sea.cover(0, ice[0]), surfaceT[0], sea.snowAlbedo[0])), T = surfaceT[0];
      flux[0] = solar + 180 + 120 * (air - 237) / 36.2 - 0.98 * sigma * T ** 4 + 12 * (air - T);
      sea.update(surfaceT, ice, flux, 0, dt, 0, solar, 0.5);
      const volume = sea.concentration[0] * ice[0];
      if (y === 9) { low = Math.min(low, volume); mean += volume / (365 * perDay); }
    }
    return { low, mean };
  };
  const alone = column(false), stored = column(true);
  console.log(`year 10 mean ice volume: ${alone.mean.toFixed(3)} m with iceLight alone, ${stored.mean.toFixed(3)} m with the brine reservoir (least ${alone.low.toFixed(2)} and ${stored.low.toFixed(2)} m)`);
  assert.ok(alone.low > 0.5 && stored.low > 0.5, 'both columns are perennial');
  assert.ok(stored.mean > alone.mean + 0.5, 'the brine thickens the perennial ice');
});

for (const melt of [20, 60]) test(`with the ice melting at ${melt} W/m² under a steady sun the brine stores ${melt < 40 ? 'most of' : 'all'} the light the ice holds, a skin at the melting point releasing none of it, as eq. 9 stores its source term`, () => {
  const dt = 1350, sun = 150, i0 = 0.18 + (0.35 - 0.18) * 0.5, h0 = 2;
  const sea = createSeaIce(model.mesh, { iceLight: true, brineReservoir: true });
  const surfaceT = new Float64Array([MELTING_POINT]), ice = new Float64Array([h0]), flux = new Float64Array(1);
  sea.load(ice, new Float64Array([1]));
  const base = melt + sun * i0 * Math.exp(-1.5 * (h0 - 0.1)) + 1.8 / (h0 / 2);
  const total = () => sea.energy(surfaceT[0], ice[0], sea.snow[0], sea.cover(0, ice[0]), sea.brine[0]);
  let inside = 0, atMelting = 0, both = 0, energy = total();
  for (let n = 0; n < 64 * 30; n++) {
    const T = surfaceT[0], h = ice[0], Q = sea.brine[0], passed = sun * i0, held = passed - passed * Math.exp(-1.5 * (h - 0.1));
    flux[0] = base - 15 * (T - MELTING_POINT);
    sea.update(surfaceT, ice, flux, 0, dt, 0, sun, 0.5);
    energy += dt * flux[0];
    assert.ok(Math.abs(total() - energy) < 1e-9 * sea.latent * 2, `step ${n}: energy ${total()} vs ${energy}`);
    if (T < Math.fround(MELTING_POINT) && Q > 0) both++;
    inside += dt * held;
    if (T >= MELTING_POINT) { assert.ok(Math.abs(sea.brine[0] - Q - dt * held) < 1e-9 * (Q + dt * held), `step ${n}: from a skin at the melting point the brine gains ${sea.brine[0] - Q} for ${dt * held}`); atMelting++; }
  }
  console.log(`melting at ${melt} W/m² for 30 days: the brine holds ${(sea.brine[0] / 1e6).toFixed(1)} MJ/m² of the ${(inside / 1e6).toFixed(1)} the ice held, the skin at the melting point in ${atMelting} of ${64 * 30} steps, filling and releasing in ${both}`);
  if (melt < 40) assert.ok(sea.brine[0] > 0.75 * inside && both > 64 * 20, `the brine keeps most of the light the ice holds, filling and releasing in ${both} steps`);
  else assert.ok(atMelting === 64 * 30 && sea.brine[0] > (1 - 1e-9) * inside, 'the brine keeps all the light the ice holds');
});

test('a melting skin saved at single precision reloads into the same step: the brine fills from the light whether it held any or not', async () => {
  const { encodeState, decodeState } = await import('../js/stateFile.module.js');
  const C = model.mesh.nCells, dt = 1350, sun = 300;
  const fresh = () => {
    const sea = createSeaIce(model.mesh, { iceLight: true, brineReservoir: true }), surfaceT = new Float64Array(C).fill(FREEZING_POINT), ice = new Float64Array(C);
    ice[0] = 2; ice[1] = 2; surfaceT[0] = MELTING_POINT; surfaceT[1] = MELTING_POINT;
    return { sea, surfaceT, ice };
  };
  const twin = fresh();
  twin.sea.load(twin.ice, Float64Array.from(twin.ice, (h) => (h > 0 ? 1 : 0)), Float64Array.from({ length: C }, (_, i) => (i === 1 ? 1e7 : 0)));
  const saved = await decodeState(encodeState({ N: 16, K: model.core.K, day: 1, time: 86400, surfaceT: twin.surfaceT, ice: twin.ice, concentration: twin.sea.concentration, brine: twin.sea.brine }));
  assert.ok(saved.surfaceT[0] < MELTING_POINT, 'the saved skin rounds below the melting point');
  const loaded = fresh();
  loaded.surfaceT.set(saved.surfaceT); loaded.ice.set(saved.ice);
  loaded.sea.load(loaded.ice, saved.concentration, saved.brine);
  for (const { sea, surfaceT, ice } of [twin, loaded]) for (const i of [0, 1]) sea.update(surfaceT, ice, new Float64Array(C).fill(50), i, dt, 0, sun, 0.5);
  for (const i of [0, 1]) {
    const start = i === 1 ? 1e7 : 0;
    assert.ok(twin.sea.brine[i] > start, `cell ${i}: the brine fills`);
    assert.ok(Math.abs(loaded.sea.brine[i] - twin.sea.brine[i]) < 1e-6 * twin.sea.brine[i], `cell ${i}: reloaded ${loaded.sea.brine[i]} against ${twin.sea.brine[i]}`);
  }
});

const TURBULENT = (extra = {}) => ({ iceExchange: { mode: 'turbulent', ...extra } });
const volumeOf = (sea, ice) => sea.concentration[0] * ice[0];

test('iceExchange is the mode name or an object of the exchange constants that names its mode, CICE\'s where it names none; the default is turbulent with a friction-velocity floor of 2 cm/s', () => {
  assert.equal(iceExchangeOptions('instant').mode, 'instant');
  assert.deepEqual(iceExchangeOptions(), iceExchangeOptions(DEFAULT_ICE_EXCHANGE));
  assert.equal(iceExchangeOptions('turbulent').turbulent, true);
  assert.deepEqual(Object.fromEntries(Object.keys(ICE_EXCHANGE).map((key) => [key, iceExchangeOptions('turbulent')[key]])), { heatTransfer: 0.006, frictionMinimum: 5e-4, lateralMelt: 1.6e-6, lateralExponent: 1.36, floeDiameter: 300, floeShape: 0.66, leadTemperature: 'freezing' });
  assert.ok(Math.abs(iceExchangeOptions('turbulent').lateralRate - 1.6e-6 * Math.PI / (0.66 * 300)) < 1e-20);
  assert.equal(iceExchangeOptions({ mode: 'turbulent', frictionMinimum: 5e-3 }).frictionMinimum, 5e-3);
  const resolved = iceExchangeOptions({ mode: 'turbulent', floeDiameter: 100 });
  assert.deepEqual(iceExchangeOptions(resolved), resolved, 'a resolved set resolves to itself');
  for (const bad of [{ frictionMinimum: 5e-3 }, {}, { ...resolved, mode: 'instant' }, { ...resolved, lateralRate: 2 * resolved.lateralRate }, { mode: 'instant', turbulent: true }]) assert.throws(() => iceExchangeOptions(bad), /iceExchange/, JSON.stringify(bad));
  for (const bad of ['gradual', { mode: 'turbulent', speed: 1 }, { heatTransfer: 0 }, { frictionMinimum: -1 }, null, 3]) assert.throws(() => iceExchangeOptions(bad), /iceExchange/, JSON.stringify(bad));
  assert.deepEqual(defaultSeaIce(model.mesh).iceExchange, iceExchangeOptions({ mode: 'turbulent', frictionMinimum: 2e-2 }));
  assert.equal(defaultSeaIce(model.mesh, { iceExchange: 'instant' }).turbulent, false);
  assert.equal(createSeaIce(model.mesh, TURBULENT()).turbulent, true);
});

test('under full ice the water gives the base ρ c_p c_h u* (T_w − T_f) as turbulence carries it, no more than the layer holds, and the ice and the water trade exactly that heat', () => {
  const rise = 0.4, conductance = 80, dt = 900;
  const run = (c, heat) => {
    const sea = createSeaIce(model.mesh, TURBULENT({ lateralMelt: 0 })), surfaceT = new Float64Array([FREEZING_POINT]), ice = new Float64Array([1]), flux = new Float64Array(1);
    sea.concentration[0] = 1; sea.conductance[0] = c; sea.waterHeat[0] = heat;
    const before = energyOf(sea, surfaceT, ice) + sea.waterHeat[0];
    sea.update(surfaceT, ice, flux, 0, dt);
    return { sea, volume: volumeOf(sea, ice), water: sea.waterHeat[0], energy: energyOf(sea, surfaceT, ice) + sea.waterHeat[0] - before };
  };
  const stored = sea0().slabHeatCapacity * rise;
  function sea0() { return createSeaIce(model.mesh); }
  const free = run(conductance, stored);
  const taken = conductance * rise * dt;
  assert.ok(Math.abs(stored - free.water - taken) < 1e-6, `the water gave ${stored - free.water} J/m², McPhee's flux times the step is ${taken}`);
  assert.ok(Math.abs((1 - free.volume) * free.sea.latent - taken) < 1e-3, `the ice melted ${(1 - free.volume) * free.sea.latent} J/m² of its base`);
  assert.ok(Math.abs(free.energy) < 1e-9 * free.sea.slabHeatCapacity, `energy moved by ${free.energy} J/m²`);
  const limited = run(1e5, stored);
  assert.ok(limited.water < 1e-6 && Math.abs((1 - limited.volume) * limited.sea.latent - stored) < 1e-3, `a conductance of 1e5 empties the layer into ${(1 - limited.volume) * limited.sea.latent} J/m² of ice, of ${stored} held`);
  assert.ok(Math.abs(limited.energy) < 1e-9 * limited.sea.slabHeatCapacity);
  const cold = run(conductance, 0);
  assert.equal(cold.volume, 1, 'water at the freezing point melts nothing');
  const weak = run(conductance / 4, stored);
  assert.ok(Math.abs((1 - weak.volume) / (1 - free.volume) - 0.25) < 1e-9, 'the flux is proportional to the conductance');
});

test('the leads\' net gain, sunlight included, warms the water under the ice and melts none of it; a net loss takes that heat first and freezes new ice only for the rest', () => {
  const dt = 900, make = (mode) => {
    const sea = createSeaIce(model.mesh, mode === 'turbulent' ? TURBULENT({ lateralMelt: 0 }) : {}), surfaceT = new Float64Array([FREEZING_POINT]), ice = new Float64Array([1]);
    sea.concentration[0] = 0.5; sea.conductance[0] = 0;
    return { sea, surfaceT, ice };
  };
  const step = (r, w, contrast = 0, heat = 0) => { r.sea.waterHeat[0] = heat; r.sea.update(r.surfaceT, r.ice, new Float64Array([w]), 0, dt, contrast); return r; };
  const gain = 0.5 * 100 * dt;
  const warm = step(make('turbulent'), 100), instant = step(make('instant'), 100);
  assert.ok(Math.abs(warm.sea.waterHeat[0] - gain) < 1e-9 && Math.abs(volumeOf(warm.sea, warm.ice) - 0.5) < 1e-15, `the water gained ${warm.sea.waterHeat[0]} of ${gain} J/m² and the ice kept ${volumeOf(warm.sea, warm.ice)} m`);
  assert.ok(volumeOf(instant.sea, instant.ice) < 0.5 - 0.99 * gain / instant.sea.latent, 'the instant exchange melts ice with the same heat');
  const sun = step(make('turbulent'), 100, 150);
  assert.ok(Math.abs(sun.sea.waterHeat[0] - 0.5 * (100 + 0.5 * 150) * dt) < 1e-9, 'sunlight absorbed in the leads is part of the gain');
  const loss = 0.5 * 100 * dt;
  const covered = step(make('turbulent'), -100, 0, 1.6 * loss);
  assert.ok(Math.abs(covered.sea.waterHeat[0] - 0.6 * loss) < 1e-9 && covered.sea.budget.leadFrozen === 0, 'a loss the water\'s heat covers freezes nothing');
  const short = step(make('turbulent'), -100, 0, 0.4 * loss);
  assert.ok(short.sea.waterHeat[0] < 1e-9, 'a loss beyond it leaves the water at freezing');
  assert.ok(Math.abs(short.sea.budget.leadFrozen / model.mesh.areaCell[0] * short.sea.latent - 0.6 * loss) < 1e-6, `new ice in the leads took ${short.sea.budget.leadFrozen / model.mesh.areaCell[0] * short.sea.latent} J/m² of the ${0.6 * loss} the water could not give`);
});

test('iceExchange leadTemperature water, turbulent only, has the air see the leads at the temperature of the water under the ice: for the grid-box flux those leads give, the ice takes what it takes with the leads at freezing and the water loses leadExchange more for each kelvin it holds above freezing', () => {
  assert.equal(iceExchangeOptions().leadTemperature, 'freezing');
  assert.equal(iceExchangeOptions({ mode: 'turbulent', leadTemperature: 'water' }).leadTemperature, 'water');
  for (const bad of [{ mode: 'instant', leadTemperature: 'water' }, { mode: 'turbulent', leadTemperature: 'sea' }]) assert.throws(() => iceExchangeOptions(bad), /leadTemperature/, JSON.stringify(bad));
  const dt = 900, A = 0.2, rise = 4, F = -50;
  const run = (leadTemperature) => {
    const sea = createSeaIce(model.mesh, TURBULENT({ lateralMelt: 0, leadTemperature })), surfaceT = new Float64Array([MELTING_POINT - 3]), ice = new Float64Array([1]);
    sea.concentration[0] = A; sea.conductance[0] = 0; sea.waterHeat[0] = sea.slabHeatCapacity * rise;
    const lead = sea.leadTemperature(0);
    sea.update(surfaceT, ice, new Float64Array([F - (1 - A) * sea.leadExchange * (lead - FREEZING_POINT)]), 0, dt);
    return { sea, lead, skin: surfaceT[0], volume: volumeOf(sea, ice) };
  };
  const cold = run('freezing'), warm = run('water');
  assert.equal(cold.lead, FREEZING_POINT);
  assert.ok(Math.abs(warm.lead - FREEZING_POINT - rise) < 1e-12, `the leads at ${warm.lead} K`);
  assert.ok(Math.abs(warm.skin - cold.skin) < 1e-12 && Math.abs(warm.volume - cold.volume) < 1e-15, `the ice skin ${warm.skin} against ${cold.skin} K, its volume ${warm.volume} against ${cold.volume} m`);
  const lost = cold.sea.waterHeat[0] - warm.sea.waterHeat[0], expected = (1 - A) * cold.sea.leadExchange * rise * dt;
  assert.ok(Math.abs(lost - expected) < 1e-9 * expected, `the water lost ${lost} J/m² more, against ${expected}`);
});

test('ice that melts away over held heat under the turbulent exchange leaves the water no heat held, its heat in the open water; over a cell with no ocean column the water under the ice is at the freezing point whatever heat was loaded there', () => {
  const sea = createSeaIce(model.mesh, TURBULENT()), surfaceT = new Float64Array([FREEZING_POINT]), ice = new Float64Array([0.004]);
  sea.concentration[0] = 0.5; sea.conductance[0] = 2000; sea.waterHeat[0] = sea.slabHeatCapacity * 2;
  sea.update(surfaceT, ice, new Float64Array(1), 0, 1350);
  assert.deepEqual([ice[0], sea.concentration[0], sea.waterHeat[0], sea.waterTemperature(0)], [0, 0, 0, FREEZING_POINT]);
  assert.ok(surfaceT[0] > FREEZING_POINT + 1.9, `the open water at ${surfaceT[0] - FREEZING_POINT} K above freezing`);
  const dry = createSeaIce(model.mesh, { ...TURBULENT(), oceanCells: new Uint8Array(model.mesh.nCells) });
  dry.waterHeat[0] = dry.slabHeatCapacity * 2;
  assert.equal(dry.waterTemperature(0), FREEZING_POINT);
});

test('floes melt at their sides as Maykut and Perovich (1987) measured: w_lat = m1 (T_w − T_f)^m2 over the floe\'s perimeter, area and volume in proportion, the thickness kept, the heat from the water', () => {
  const dt = 900, ice0 = 1, A0 = 0.8, run = (rise, options = {}) => {
    const sea = createSeaIce(model.mesh, TURBULENT(options)), surfaceT = new Float64Array([FREEZING_POINT]), ice = new Float64Array([ice0]);
    sea.concentration[0] = A0; sea.conductance[0] = 0; sea.waterHeat[0] = sea.slabHeatCapacity * rise;
    const stored = sea.waterHeat[0];
    sea.update(surfaceT, ice, new Float64Array(1), 0, dt);
    return { sea, area: sea.concentration[0], thickness: ice[0], taken: stored - sea.waterHeat[0], r: 1 - sea.concentration[0] / A0 };
  };
  const two = run(2), one = run(1);
  const expected = 1.6e-6 * 2 ** 1.36 * Math.PI / (0.66 * 300) * dt;
  assert.ok(Math.abs(two.r - expected) < 1e-12, `the cover lost ${two.r} of itself, w_lat π/(α D) dt is ${expected}`);
  assert.ok(Math.abs(two.thickness - ice0) < 1e-12, `the floes kept their thickness: ${two.thickness}`);
  assert.ok(Math.abs(two.taken - two.sea.latent * expected * A0 * ice0) < 1e-3, `the water gave ${two.taken} J/m² for ${two.sea.latent * expected * A0 * ice0}`);
  assert.ok(Math.abs(two.r / one.r - 2 ** 1.36) < 1e-9, `doubling T_w − T_f multiplies the loss by ${two.r / one.r}, 2^1.36 is ${2 ** 1.36}`);
  assert.equal(run(2, { lateralMelt: 0 }).r, 0, 'a zero coefficient melts nothing at the sides');
  assert.ok(Math.abs(run(2, { floeDiameter: 150 }).r / two.r - 2) < 1e-9, 'smaller floes melt faster by their perimeter');
  const big = run(1e4);
  assert.ok(big.taken <= big.sea.slabHeatCapacity * 1e4 + 1e-6 && big.r <= 1, 'a layer\'s lateral melt cannot exceed the floe or the heat held');
});

test('with the turbulent exchange the ice and the water under it conserve energy exactly through partial melt, lead freezing, thaw and refreezing, and the ice outlasts the instant exchange\'s', () => {
  const run = (options, conductance) => {
    const sea = createSeaIce(model.mesh, options), surfaceT = new Float64Array([FREEZING_POINT - 5]), ice = new Float64Array([0.6]), flux = new Float64Array(1);
    sea.concentration[0] = 0.5; sea.snow[0] = 5; sea.conductance[0] = conductance;
    const dt = 900, day = 96, scale = sea.slabHeatCapacity * 2, total = () => energyOf(sea, surfaceT, ice) + (ice[0] > 0 ? sea.waterHeat[0] : 0);
    let energy = total(), thawed = null, warm = 0, lateral = 0;
    for (let n = 0; n < 200 * day; n++) {
      const phase = n < 40 * day ? 'winter' : n < 120 * day ? 'summer' : 'autumn';
      flux[0] = phase === 'summer' ? 160 : phase === 'winter' ? -120 : -200;
      if (phase === 'winter' && n % day === 0 && sea.deposit(0, 1.5, MELTING_POINT - 10, ice, surfaceT)) energy -= sea.latentHeatFusion * 1.5;
      sea.update(surfaceT, ice, flux, 0, dt, phase === 'summer' ? 120 : 0);
      energy += dt * flux[0];
      assert.ok(Math.abs(total() - energy) < 1e-9 * scale, `step ${n} (${phase}): energy ${total()} vs ${energy}`);
      if (ice[0] > 0) warm = Math.max(warm, sea.waterHeat[0] / sea.slabHeatCapacity);
      else if (thawed === null) thawed = n / day;
      lateral = sea.budget.lateralMelted;
    }
    return { thawed, warm, lateral, sea };
  };
  const instant = run({}, 0), turbulent = run(TURBULENT(), 60);
  console.log(`a half-covered cell through a summer: the instant exchange thaws it on day ${instant.thawed?.toFixed(1)}, the turbulent exchange on day ${turbulent.thawed?.toFixed(1)}, the water under the ice warming to ${turbulent.warm.toFixed(2)} K above freezing`);
  assert.ok(turbulent.warm > 0.05, `the water under the ice reached ${turbulent.warm} K above freezing`);
  assert.ok(turbulent.thawed === null || instant.thawed === null || turbulent.thawed > instant.thawed, `the turbulent exchange thaws on day ${turbulent.thawed} against ${instant.thawed}`);
  assert.ok(turbulent.sea.budget.leadFrozen > 0 && turbulent.lateral > 0);
});

test('ice that melts away leaves the water with the heat it held above freezing, and ice that forms in open water leaves it at freezing', () => {
  const sea = createSeaIce(model.mesh, TURBULENT({ lateralMelt: 0 })), dt = 900, flux = new Float64Array([300]);
  const surfaceT = new Float64Array([MELTING_POINT]), ice = new Float64Array([0.0003]);
  sea.concentration[0] = 1; sea.waterHeat[0] = 3e6; sea.conductance[0] = 20;
  const before = energyOf(sea, surfaceT, ice) + sea.waterHeat[0];
  sea.update(surfaceT, ice, flux, 0, dt);
  assert.equal(ice[0], 0, 'the thin ice melted away');
  assert.ok(Math.abs(sea.slabHeatCapacity * (surfaceT[0] - FREEZING_POINT) - (before + dt * flux[0])) < 1e-9 * sea.slabHeatCapacity, `the open water holds ${sea.slabHeatCapacity * (surfaceT[0] - FREEZING_POINT)} J/m² above freezing, ${before + dt * flux[0]} came in`);
  const open = new Float64Array([FREEZING_POINT + 0.01]), none = new Float64Array([0]);
  sea.waterHeat[0] = 7e5;
  for (let n = 0; n < 400 && none[0] === 0; n++) sea.update(open, none, new Float64Array([-250]), 0, dt);
  assert.ok(none[0] > 0 && sea.waterHeat[0] === 0 && open[0] === FREEZING_POINT, `ice formed at ${none[0]} m with the water at ${sea.waterHeat[0]} J/m²`);
});
