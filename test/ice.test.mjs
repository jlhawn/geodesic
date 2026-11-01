import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createSeaIce, openWaterAlbedo, iceExchangeOptions, ICE_EXCHANGE, FREEZING_POINT, MELTING_POINT } from '../js/physics/ice.module.js';
import { initializeState } from '../js/physics/init.module.js';

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

const TURBULENT = (extra = {}) => ({ iceExchange: { mode: 'turbulent', ...extra } });
const volumeOf = (sea, ice) => sea.concentration[0] * ice[0];

test('iceExchange is the mode name or an object of the exchange constants, CICE\'s by default', () => {
  assert.equal(iceExchangeOptions().mode, 'instant');
  assert.equal(iceExchangeOptions('turbulent').turbulent, true);
  assert.deepEqual(Object.fromEntries(Object.keys(ICE_EXCHANGE).map((key) => [key, iceExchangeOptions('turbulent')[key]])), { heatTransfer: 0.006, frictionMinimum: 5e-4, lateralMelt: 1.6e-6, lateralExponent: 1.36, floeDiameter: 300, floeShape: 0.66 });
  assert.ok(Math.abs(iceExchangeOptions('turbulent').lateralRate - 1.6e-6 * Math.PI / (0.66 * 300)) < 1e-20);
  assert.equal(iceExchangeOptions({ frictionMinimum: 5e-3 }).mode, 'turbulent', 'an object of constants without a mode is the turbulent exchange');
  assert.deepEqual(iceExchangeOptions(iceExchangeOptions({ floeDiameter: 100 })).floeDiameter, 100, 'a resolved set resolves to itself');
  for (const bad of ['gradual', { mode: 'turbulent', speed: 1 }, { heatTransfer: 0 }, { frictionMinimum: -1 }, null, 3]) assert.throws(() => iceExchangeOptions(bad), /iceExchange/, JSON.stringify(bad));
  assert.equal(createSeaIce(model.mesh).turbulent, false);
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
