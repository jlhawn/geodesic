import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createSeaIce, openWaterAlbedo, FREEZING_POINT, MELTING_POINT } from '../js/physics/ice.module.js';
import { initializeState } from '../js/physics/init.module.js';

const model = createModel(new Grid(2));
const seaIce = createSeaIce(model.mesh);

test('the surface energy changes by exactly the surface flux through freezing, growth, melting and thaw', () => {
  const surfaceT = new Float64Array([FREEZING_POINT + 2]), ice = new Float64Array([0]);
  const flux = new Float64Array(1);
  const dt = 900;
  let energy = seaIce.energy(surfaceT[0], ice[0]);
  const scale = seaIce.slabHeatCapacity * 2;
  let froze = false, melted = false;
  for (let n = 0; n < 4 * 96 * 60; n++) {
    flux[0] = n < 96 * 90 ? -150 : 250;
    seaIce.update(surfaceT, ice, flux, 0, dt);
    energy += dt * flux[0];
    assert.ok(Math.abs(seaIce.energy(surfaceT[0], ice[0]) - energy) < 1e-9 * scale, `step ${n}: energy ${seaIce.energy(surfaceT[0], ice[0])} vs ${energy}`);
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
  const m = createModel(new Grid(4));
  const init = initializeState(m, {});
  for (let a = 0; a < init.length; a++) m.state[a].set(init[a]);
  const d0 = m.diagnostics();
  assert.ok(d0.iceFraction > 0.02 && d0.iceFraction < 0.3, `initial ice fraction ${d0.iceFraction}`);
  for (let n = 0; n < 24; n++) m.step(1800);
  const d = m.diagnostics();
  assert.ok(Number.isFinite(d.maxWind) && d.maxWind < 80);
  assert.ok(d.planetaryAlbedo > 0.05 && d.planetaryAlbedo < 0.6, `planetary albedo ${d.planetaryAlbedo}`);
  assert.ok(d.iceFraction > 0);
  console.log(`N=4 after 12 h: ice fraction ${d.iceFraction.toFixed(3)}, mean thickness ${d.iceThickness.toFixed(2)} m, surface albedo ${d.surfaceAlbedo.toFixed(3)}, planetary albedo ${d.planetaryAlbedo.toFixed(3)}, TCW ${(1000 * d.columnCloud).toFixed(0)} g/m²`);
});

test('snow on the ice brightens it, insulates it, melts before it, and goes into the water when the ice is gone', () => {
  const snowy = createSeaIce(model.mesh), bare = createSeaIce(model.mesh);
  const flux = (w) => new Float64Array([w]);
  assert.ok(Math.abs(snowy.albedo(1, null, 0) - 0.5) < 1e-12 && Math.abs(snowy.albedo(1, null, 20) - 0.75) < 1e-12 && Math.abs(snowy.albedo(1, null, 10) - 0.625) < 1e-12);
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
  snowy.snow[0] = 0;
  let energy = snowy.energy(surfaceT[0], thin[0], snowy.snow[0]);
  const scale = snowy.slabHeatCapacity * 2, dt = 900;
  for (let n = 0; n < 96 * 40; n++) {
    if (n % 96 === 0 && n < 96 * 10) { snowy.deposit(0, 2, MELTING_POINT - 5, thin); energy -= snowy.latentHeatFusion * 2; }
    const w = n < 96 * 10 ? -50 : 300;
    snowy.update(surfaceT, thin, flux(w), 0, dt);
    energy += dt * w;
    assert.ok(Math.abs(snowy.energy(surfaceT[0], thin[0], snowy.snow[0]) - energy) < 1e-9 * scale, `step ${n}: energy ${snowy.energy(surfaceT[0], thin[0], snowy.snow[0])} vs ${energy}`);
  }
  assert.equal(thin[0], 0, 'the ice melted away');
  assert.equal(snowy.snow[0], 0, 'and took its snow into the water');

  const tA = new Float64Array([MELTING_POINT]), iceA = new Float64Array([0.05]), tC = new Float64Array([MELTING_POINT]), iceC = new Float64Array([0.05]);
  snowy.snow[0] = 10;
  for (let n = 0; n < 96; n++) { snowy.update(tA, iceA, flux(400), 0, dt); bare.update(tC, iceC, flux(400), 0, dt); }
  assert.ok(iceA[0] === 0 && iceC[0] === 0);
  assert.ok(Math.abs((tC[0] - tA[0]) - 10 * snowy.latentHeatFusion / snowy.slabHeatCapacity) < 1e-9, `the snow's latent heat came out of the water: ${tC[0] - tA[0]} K`);
});
