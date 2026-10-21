import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createOcean } from '../js/ocean/layered.module.js';
import { seawaterDensity } from '../js/ocean/seawater.module.js';
import { createSeaIce, openWaterAlbedo, FREEZING_POINT } from '../js/physics/ice.module.js';
import { SOLAR_CONSTANT } from '../js/physics/radiation.module.js';
import { dailyMean, sunTrack } from '../explain/insolation.module.js';
import { seaIceColumn, iceGrowth, RHO_CP, brine, iceForDensity, dailySunlight, DEEP_DENSITY, DEEP_WATER, DEPTHS, ICE_SALINITY, ICE_DENSITY } from '../explain/seaIceColumn.module.js';

const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const C = (k) => k - 273.15;

test('a meter of ice whose top is 20 K below freezing conducts 40 W/m² and grows 1.1 cm a day, and the Stefan law gives 1.4, 2.0 and 2.5 m in 180 days at 10, 20 and 30 K', () => {
  const day = iceGrowth();
  console.log(`1 m of ice across 20 K conducts ${day.conduction.toFixed(2)} W/m² and grows ${(100 * day.perDay).toFixed(3)} cm in a day`);
  close(day.conduction, 40, 1e-9, 'conduction');
  close(100 * day.perDay, 1.12, 0.005, 'growth in cm a day');
  const stefan = [10, 20, 30].map((below) => iceGrowth({ thickness: 0.1, below, days: 180 }).thickness);
  console.log(`from 0.1 m, 180 days at 10, 20 and 30 K below freezing grow ${stefan.map((h) => h.toFixed(2)).join(', ')} m`);
  close(stefan[0], 1.43, 0.005, '10 K'); close(stefan[1], 2.02, 0.005, '20 K'); close(stefan[2], 2.47, 0.005, '30 K');
});

test('seawater freezes at −1.8 °C, cools only a few thousandths of a kelvin below it before ice forms, and sits at it under the ice', () => {
  close(C(FREEZING_POINT), -1.8, 1e-9, 'the freezing point');
  const sea = createSeaIce({ nCells: 1, areaCell: [1] }, { heatCapacity: Float64Array.of(50 * 1025 * 3985) });
  const T = Float64Array.of(FREEZING_POINT + 0.1), ice = new Float64Array(1), cooling = Float64Array.of(-100);
  let coldest = T[0], hours = 0;
  while (ice[0] === 0 && hours < 1000) { sea.update(T, ice, cooling, 0, 3600); hours++; if (ice[0] === 0) coldest = Math.min(coldest, T[0]); }
  console.log(`50 m of water losing 100 W/m² reaches ${(FREEZING_POINT - coldest).toFixed(4)} K below freezing before its first ice, after ${hours} hours: ${ice[0].toFixed(3)} m over ${(100 * sea.concentration[0]).toFixed(1)} percent of the cell`);
  assert.ok(coldest < FREEZING_POINT && coldest > FREEZING_POINT - 0.005 && ice[0] > 0, `coldest open water ${coldest}`);
  close(T[0], FREEZING_POINT, 1e-9, 'the new ice\'s water');
  const run = seaIceColumn({ ice: 'dark', years: 4 }), r = run.record;
  let iced = 0;
  for (let d = 0; d < run.days; d++) if (r.thickness[d] > 0) { iced++; assert.equal(r.water[d], Math.fround(FREEZING_POINT), `day ${d}: water under ice`); }
  assert.ok(iced > 100, `${iced} days with ice`);
  close(seawaterDensity(FREEZING_POINT, 34), 1026.47, 0.005, '34 psu at freezing');
  close(DEEP_DENSITY, 1026.77, 0.005, 'the deep water');
});

test('the deep water the gauge compares with is the bottom of the model\'s own start: 2 °C at 34.7 psu, 1026.77 kg/m³, at every latitude', () => {
  const m = buildMesh(new Grid(16)), ocean = createOcean(m, { everySteps: 1 }), surfaceT = new Float64Array(m.nCells).fill(288), n = m.nCells, k = ocean.layers - 1;
  ocean.initialize(surfaceT, new Float64Array(n));
  let worst = 0;
  for (let i = 0; i < n; i++) {
    const h = ocean.h[k * n + i], T = ocean.Q[k * n + i] / h, S = ocean.W[k * n + i] / h;
    close(T, DEEP_WATER.T, 0.05, `cell ${i}: bottom layer temperature`); close(S, DEEP_WATER.S, 0.01, `cell ${i}: bottom layer salinity`);
    worst = Math.max(worst, Math.abs(seawaterDensity(T, S) - DEEP_DENSITY));
  }
  console.log(`the bottom layer of the start is within ${worst.toFixed(4)} kg/m³ of ${DEEP_DENSITY.toFixed(3)}`);
  assert.ok(worst < 0.02, `${worst}`);
  close(DEEP_WATER.T, 273.15 + 2, 1e-9, '2 °C'); close(seawaterDensity(FREEZING_POINT, 34.38), DEEP_DENSITY, 0.005, 'water at freezing reaches it near 34.38 psu');
  close(iceForDensity({ depth: 50 }), 0.711, 0.002, 'meters of ice over 50 m that salt 34 psu water at freezing to it');
  close(iceForDensity({ depth: DEPTHS.slab }), 0.146, 0.002, '... over the slab\'s 10.28 m'); close(iceForDensity({ depth: DEPTHS.deep }), 2.846, 0.005, '... over 200 m');
});


test('the feedback: over 50 m at 75°N with 85 W/m² brought in, the model\'s bright ice lasts the year 3.6 to 5.3 m thick, where Part 7\'s dark ice is gone for 274 days; deeper water stays open longer', () => {
  const bright = seaIceColumn().year(), dark = seaIceColumn({ ice: 'dark' }).year();
  const describe = (y) => `${y.thinnest.toFixed(2)}–${y.thickest.toFixed(2)} m of ice, ${y.iceFreeDays} days without any, water ${C(y.waterLow).toFixed(2)} to ${C(y.waterHigh).toFixed(2)} °C, surface down to ${C(y.surfaceLow).toFixed(1)} °C, salinity ${y.salinityLow.toFixed(3)}–${y.salinityHigh.toFixed(3)}, densest ${(y.againstDeep >= 0 ? '+' : '') + y.againstDeep.toFixed(3)} kg/m³ against the deep water`;
  console.log(`bright ice: ${describe(bright)}`);
  console.log(`dark ice:   ${describe(dark)}`);
  close(bright.thinnest, 3.55, 0.01, 'bright: thinnest'); close(bright.thickest, 5.29, 0.01, 'bright: thickest');
  assert.equal(bright.iceFreeDays, 0, 'the bright ice never goes');
  close(dark.thickest, 0.80, 0.01, 'dark: thickest');
  assert.equal(dark.iceFreeDays, 274, 'days without dark ice');
  close(C(dark.waterHigh), 5.35, 0.01, 'the open water under dark ice\'s summer');
  close(dark.salinityGained, 0.428, 0.001, 'the brine of a winter of dark ice');
  close(dark.againstDeep, 0.039, 0.002, 'dark ice: the densest water is just denser than the deep water');
  close(bright.againstDeep, 2.05, 0.02, 'the brine of perennial ice over a layer that never mixes down');
  const open = [DEPTHS.slab, 50, 200].map((depth) => seaIceColumn({ ice: 'dark', depth }).year().iceFreeDays);
  console.log(`dark ice is gone for ${open.join(', ')} days over 10.28, 50 and 200 m of water`);
  assert.deepEqual(open, [194, 274, 365], 'the deeper the water, the longer it stays open');
});

test('the sunlight is the model\'s daily mean, and the open water\'s albedo is its own, weighted by the day\'s light', () => {
  const light = dailySunlight(75);
  for (const day of [0, 45, 91, 183, 274]) {
    close(light.mean[day], dailyMean(75, day), 1e-9, `day ${day}: daily mean`);
    if (light.mean[day] === 0) continue;
    const { x, z } = sunTrack(day, 288), c = Math.cos(75 * Math.PI / 180), s = Math.sin(75 * Math.PI / 180);
    let absorbed = 0;
    for (let k = 0; k < 288; k++) { const mu = c * x[k] + s * z[k]; if (mu > 0) absorbed += SOLAR_CONSTANT * mu * (1 - openWaterAlbedo(mu)) / 288; }
    close(light.mean[day] * (1 - openWaterAlbedo(light.sun[day])), absorbed, 1e-9 * absorbed, `day ${day}: open water absorbs`);
  }
  console.log(`at 75°N the June solstice brings ${light.mean[91].toFixed(1)} W/m² to the top of the air, and open water reflects ${(100 * light.water[91]).toFixed(1)} percent of it, ${(100 * light.water[0]).toFixed(1)} at the March equinox`);
  close(light.water[91], 0.095, 0.0005, 'the solstice\'s open-water albedo');
});

test('twenty-four model years of hourly steps take well under a second', () => {
  seaIceColumn({ years: 1 });
  const t0 = performance.now();
  const run = seaIceColumn({ years: 24 });
  const ms = performance.now() - t0;
  console.log(`24 years, ${24 * 365 * 24} hourly steps, in ${ms.toFixed(0)} ms`);
  assert.equal(run.days, 24 * 365);
  assert.ok(ms < 500, `${ms} ms`);
});

test('every meter of ice frozen out of 50 m at 34 psu adds 0.53 psu and 0.42 kg/m³; water at freezing reaches the deep water\'s density after 0.71 m of ice, at 34.38 psu; the first winter ends 0.77 kg/m³ denser than the deep water under bright ice and 0.67 under dark', () => {
  close(brine(34, 1, 50), 0.532, 0.001, 'psu per meter of ice');
  const per = seawaterDensity(FREEZING_POINT, 34 + brine(34, 1, 50)) - seawaterDensity(FREEZING_POINT, 34);
  close(per, 0.423, 0.001, 'kg/m³ per meter of ice');
  const ice = iceForDensity({ depth: 50 });
  close(ice, 0.711, 0.001, 'ice for the deep water\'s density'); close(34 + brine(34, ice, 50), 34.378, 0.001, 'its salinity');
  const first = (ice) => seaIceColumn({ ice, startDay: 183, start: 273.15 }).year(0).againstDeep;
  close(first('bright'), 0.772, 0.005, 'bright, first winter'); close(first('dark'), 0.673, 0.005, 'dark, first winter');
  console.log(`${per.toFixed(3)} kg/m³ per meter of ice over 50 m; the deep water's density after ${ice.toFixed(3)} m of ice`);
});

test('the ocean under ice: the top layer is 3 m, a centimeter of ice salts it by 0.09 psu and makes it denser than the layer beneath, and the wind\'s stress reaches the water at 0.8 under full cover', () => {
  const m = buildMesh(new Grid(16)), ocean = createOcean(m, { everySteps: 1 }), n = m.nCells, surfaceT = new Float64Array(n).fill(274), ice = new Float64Array(n);
  ocean.initialize(surfaceT, ice);
  close(ocean.h[0], 3, 1e-9, 'the top layer, meters'); close(ocean.capacity[0], 3 * RHO_CP, 1e-3 * RHO_CP, 'its capacity');
  const S0 = ocean.W[0] / ocean.h[0], perCentimeter = (S0 - ICE_SALINITY) * ICE_DENSITY / 1000 * 0.01 / 3;
  close(perCentimeter, 0.09, 0.01, 'psu a centimeter of ice adds');
  assert.ok(seawaterDensity(ocean.Q[0] / 3, S0 + perCentimeter) > seawaterDensity(ocean.Q[n] / ocean.h[n], ocean.W[n] / ocean.h[n]), 'salted, the top layer is denser than the layer beneath');
  console.log(`a centimeter of ice salts the 3 m top layer (${S0.toFixed(2)} psu) by ${perCentimeter.toFixed(3)} psu`);
  const stress = new Float64Array(m.nEdges).fill(1);
  ocean.setStress(stress, new Float64Array(n).fill(1));
  assert.ok(Array.from(ocean.stress).every((v) => Math.abs(v - 0.8) < 1e-12), 'under full ice cover 0.8 of the stress');
  ocean.setStress(stress, new Float64Array(n));
  assert.ok(Array.from(ocean.stress).every((v) => v === 1), 'open water 1');
});
