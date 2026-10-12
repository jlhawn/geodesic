import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createOcean } from '../js/ocean/layered.module.js';
import { seawaterDensity } from '../js/ocean/seawater.module.js';
import { createSeaIce, openWaterAlbedo, FREEZING_POINT } from '../js/physics/ice.module.js';
import { SOLAR_CONSTANT } from '../js/physics/radiation.module.js';
import { dailyMean, sunTrack } from '../explain/insolation.module.js';
import { seaIceColumn, iceGrowth, brine, iceForDensity, dailySunlight, DEEP_CLASS, DEPTHS, ICE_SALINITY, ICE_DENSITY } from '../explain/seaIceColumn.module.js';

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
  assert.equal(DEEP_CLASS, 1026.95, 'the class beneath');
});

test('brine: a meter of ice frozen out of a 50 m layer adds 0.532 psu and 0.42 kg/m³, 1.13 m makes it as dense as the 1026.95 class at 34.6 psu, and the ocean books brine the same way, by volume', () => {
  close(brine(34, 1, 50), 0.532, 0.0005, 'psu from a meter of ice');
  close(seawaterDensity(FREEZING_POINT, 34 + brine(34, 1, 50)) - seawaterDensity(FREEZING_POINT, 34), 0.42, 0.005, 'density from a meter of ice');
  const needed = iceForDensity();
  close(needed, 1.13, 0.005, 'ice to reach the class'); close(34 + brine(34, needed, 50), 34.6, 0.005, 'salinity at the class');
  close(iceForDensity({ depth: DEPTHS.slab }), 0.232, 0.001, 'from Part 7\'s slab'); close(iceForDensity({ depth: 200 }), 4.52, 0.005, 'from 200 m');
  console.log(`a meter of ice salts 50 m of 34 psu water by ${brine(34, 1, 50).toFixed(4)} psu; ${needed.toFixed(3)} m brings it to the ${DEEP_CLASS} class (${iceForDensity({ depth: DEPTHS.slab }).toFixed(3)} m from a 10.28 m slab, ${iceForDensity({ depth: 200 }).toFixed(2)} m from 200 m)`);

  const mesh = buildMesh(new Grid(2)), n = mesh.nCells;
  const ocean = createOcean(mesh, { everySteps: 1, thermoclineTilt: 0, salinityProfile: () => 34, mixedDepth: 50, convectiveRate: 0, stirring: 0 });
  const surfaceT = new Float64Array(n).fill(FREEZING_POINT + 0.001), ice = new Float64Array(n), flux = new Float64Array(n);
  const sea = createSeaIce(mesh, { heatCapacity: ocean.capacity });
  ocean.useSeaIce(sea.concentration);
  ocean.initialize(surfaceT, ice);
  const S0 = ocean.W[0] / ocean.h[0], depth = ocean.h[0];
  flux[0] = -400;
  for (let k = 0; k < 2; k++) sea.update(surfaceT, ice, flux, 0, 3600);
  ocean.advance(surfaceT, ice, sea.oceanFlux, new Float64Array(mesh.nEdges), 3600, sea.concentration);
  const volume = sea.concentration[0] * ice[0];
  close(ocean.W[0] / ocean.h[0] - S0, brine(S0, volume, depth), 1e-9, `the ocean's brine for ${volume.toFixed(4)} m of ice over ${(100 * sea.concentration[0]).toFixed(1)} percent of the cell`);
  assert.equal(ICE_SALINITY, 5); assert.equal(ICE_DENSITY, 917);

  const run = seaIceColumn({ years: 6 }), r = run.record;
  let worst = 0;
  for (let d = 0; d < run.days; d++) worst = Math.max(worst, Math.abs(r.salinity[d] - (ICE_SALINITY + (34 - ICE_SALINITY) * Math.exp(ICE_DENSITY * r.volume[d] / 1000 / 50))));
  assert.ok(worst < 1e-3, `the column's salinity follows its ice volume to ${worst} psu`);
  const seasonal = seaIceColumn({ ice: 'dark' }).year();
  close(seasonal.salinityDrift, 0, 1e-4, 'a year of seasonal ice gives its salt back');
});

test('the feedback: over 50 m at 75°N with 85 W/m² brought in, the model\'s bright ice lasts the year 3.6 to 5.3 m thick, where Part 7\'s dark ice is gone for 274 days; deeper water stays open longer', () => {
  const bright = seaIceColumn().year(), dark = seaIceColumn({ ice: 'dark' }).year();
  const describe = (y) => `${y.thinnest.toFixed(2)}–${y.thickest.toFixed(2)} m of ice, ${y.iceFreeDays} days without any, water ${C(y.waterLow).toFixed(2)} to ${C(y.waterHigh).toFixed(2)} °C, surface down to ${C(y.surfaceLow).toFixed(1)} °C, salinity ${y.salinityLow.toFixed(3)}–${y.salinityHigh.toFixed(3)}, densest ${(y.againstClass >= 0 ? '+' : '') + y.againstClass.toFixed(3)} kg/m³ against the class`;
  console.log(`bright ice: ${describe(bright)}`);
  console.log(`dark ice:   ${describe(dark)}`);
  close(bright.thinnest, 3.55, 0.01, 'bright: thinnest'); close(bright.thickest, 5.29, 0.01, 'bright: thickest');
  assert.equal(bright.iceFreeDays, 0, 'the bright ice never goes');
  close(dark.thickest, 0.80, 0.01, 'dark: thickest');
  assert.equal(dark.iceFreeDays, 274, 'days without dark ice');
  close(C(dark.waterHigh), 5.35, 0.01, 'the open water under dark ice\'s summer');
  close(dark.salinityGained, 0.428, 0.001, 'the brine of a winter of dark ice');
  close(dark.againstClass, -0.137, 0.001, 'dark ice: the densest water is still lighter than the class');
  assert.ok(bright.againstClass > 1.8, `the brine of perennial ice over a layer that never mixes down: ${bright.againstClass}`);
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
