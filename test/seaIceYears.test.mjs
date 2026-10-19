import { test } from 'node:test';
import assert from 'node:assert/strict';
import { createSeaIce, FREEZING_POINT, MELTING_POINT } from '../js/physics/ice.module.js';
import { seawaterDensity } from '../js/ocean/seawater.module.js';
import { seasonPhrase } from '../explain/insolation.module.js';
import { seaIceYears, conduction, iceSeason, heatAboveFreezing, iceEquivalent, snowDepth, SINKING_SALINITY, DEEP_DENSITY, DEPTHS, START_DAY, YEARS } from '../explain/seaIceYears.module.js';

const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const C = (k) => k - 273.15;
const phrase = (run, d) => seasonPhrase(run.record.day[d]);

test('the conduction through ice and snow is the model\'s: the ice its update freezes on or melts off in one second, times the latent heat', () => {
  for (const [thickness, snow, below] of [[1, 0, 20], [2.5, 0, 35], [0.05, 0, 10], [1.2, 100, 30], [3, 300, 25], [2, 0, -1.8]]) {
    const sea = createSeaIce({ nCells: 1, areaCell: [1] }, { heatCapacity: Float64Array.of(50 * 1025 * 3985) });
    const ice = Float64Array.of(thickness), skin = Float64Array.of(FREEZING_POINT - below);
    sea.load(ice);
    sea.snow[0] = snow;
    const expected = conduction(thickness, snow, skin[0]);
    sea.update(skin, ice, new Float64Array(1), 0, 1);
    close((sea.cover(0, ice[0]) * ice[0] - thickness) * sea.latent, expected, 1e-3, `${thickness} m of ice under ${snow} kg/m² of snow, top ${below} K below freezing`);
  }
  close(conduction(1, 0, FREEZING_POINT - 20), 40, 1e-12, 'a meter across 20 K');
  close(snowDepth(100), 1 / 3, 1e-12, 'a year of snow, in meters');
  close(iceEquivalent(100), 2.15, 0.005, 'a year of snow insulates like this much ice');
  assert.equal(conduction(0, 0, FREEZING_POINT - 20), 0);
});

test('a year\'s ice season is its longest run of days with ice, and says when it carries over', () => {
  const days = (pattern) => Float32Array.from(pattern, (c) => (c === '#' ? 1 : 0));
  assert.deepEqual(iceSeason(days('..####...'), 0, 9), { forms: 2, melts: 6, days: 4, allYear: false, none: false });
  assert.deepEqual(iceSeason(days('.#..####.'), 0, 9), { forms: 4, melts: 8, days: 4, allYear: false, none: false });
  assert.deepEqual(iceSeason(days('###..###...'), 3, 11), { forms: 5, melts: 8, days: 3, allYear: false, none: false });
  assert.deepEqual(iceSeason(days('######...'), 3, 9), { forms: null, melts: 6, days: 3, allYear: false, none: false });
  assert.deepEqual(iceSeason(days('...######'), 0, 6), { forms: 3, melts: null, days: 3, allYear: false, none: false });
  assert.deepEqual(iceSeason(days('#########'), 3, 6), { forms: null, melts: null, days: 3, allYear: true, none: false });
  assert.deepEqual(iceSeason(days('.........'), 0, 9), { forms: null, melts: null, days: 0, allYear: false, none: true });
});

test('the figure\'s column: open water at 0 °C on the September equinox freezes 45 days later, and the brine of 0.72 m of ice makes it as dense as the deep water', () => {
  const run = seaIceYears(), r = run.record;
  assert.equal(r.day[0], START_DAY);
  assert.equal(run.days, YEARS * 365);
  close(heatAboveFreezing(MELTING_POINT, DEPTHS.fifty) / 1e6, 367.6, 0.05, 'MJ/m² to lose before freezing');
  close(run.summary(0).summerHeat / 1e6, 367.6, 0.05, 'the first year\'s most heat above freezing is the start\'s');
  let first = 0;
  while (!(r.coverMax[first] > 0)) first++;
  assert.equal(phrase(run, first), '45 days after the September equinox');
  for (let d = first; d < run.days; d++) if (r.thickness[d] > 0) assert.equal(r.water[d], Math.fround(FREEZING_POINT), `day ${d}: the water under the ice`);
  const solstice = 91;
  assert.equal(phrase(run, solstice), 'the December solstice');
  const flow = run.conducted[solstice], growth = 100 * 86400 * flow / run.sea.latent;
  console.log(`on the first December solstice: ${r.thickness[solstice].toFixed(2)} m of ice over ${(100 * r.cover[solstice]).toFixed(0)}% of the sea, top at ${C(run.skin[solstice]).toFixed(1)} °C, conducting ${flow.toFixed(1)} W/m² and freezing ${growth.toFixed(2)} cm a day onto its base`);
  close(flow, 40.3, 0.1, 'the solstice\'s conduction'); close(growth, 1.14, 0.01, 'cm a day');
  let crossed = 0;
  while (!(r.density[crossed] > DEEP_DENSITY)) crossed++;
  console.log(`the water first passes the deep water's density, ${DEEP_DENSITY.toFixed(3)} kg/m³, ${phrase(run, crossed)}, under ${r.volume[crossed].toFixed(3)} m of ice, at ${r.salinity[crossed].toFixed(3)} psu (water at freezing is as dense as the deep water at ${SINKING_SALINITY.toFixed(3)} psu)`);
  assert.equal(phrase(run, crossed), '19 days before the December solstice');
  close(r.volume[crossed], 0.72, 0.01, 'ice at the crossing');
  close(SINKING_SALINITY, 34.38, 0.005, 'the salinity of the deep water\'s density at freezing');
  close(seawaterDensity(FREEZING_POINT, SINKING_SALINITY), DEEP_DENSITY, 1e-9, 'its density');
});

test('the feedback: the model\'s bright ice survives every summer after the first and thickens for years, while Part 7\'s dark ice settles into a season of three months', () => {
  const bright = seaIceYears(), dark = seaIceYears({ ice: 'dark' });
  const b1 = bright.summary(0), d1 = dark.summary(0), b = bright.summary(), d = dark.summary();
  console.log(`first summer: the bright ice is gone ${phrase(bright, b1.melts)}, the dark ice ${phrase(dark, d1.melts)}`);
  console.log(`year ${YEARS}: bright ${b.thinnest.toFixed(2)}–${b.thickest.toFixed(2)} m all year, salinity ${b.salinityLow.toFixed(2)}–${b.salinityHigh.toFixed(2)}, ${b.againstDeep.toFixed(2)} kg/m³ past the deep water`);
  console.log(`year ${YEARS}: dark from ${phrase(dark, d.forms)} to ${phrase(dark, d.melts)}, up to ${d.thickest.toFixed(2)} m, water up to ${C(d.waterHigh).toFixed(2)} °C holding ${(d.summerHeat / 1e9).toFixed(2)} GJ/m² (${(d.summerHeat / dark.sea.latent).toFixed(1)} m of ice), brine ${d.brine.toFixed(2)} psu, densest ${d.againstDeep.toFixed(3)} kg/m³ against the deep water`);
  assert.equal(phrase(bright, b1.melts), '26 days before the September equinox');
  assert.equal(phrase(dark, d1.melts), '11 days before the June solstice');
  for (let y = 1; y < YEARS; y++) assert.ok(bright.summary(y).allYear, `bright ice all year in year ${y + 1}`);
  close(b.thinnest, 3.58, 0.01, 'bright: thinnest'); close(b.thickest, 5.32, 0.01, 'bright: thickest');
  close(b.againstDeep, 2.06, 0.02, 'perennial ice keeps its brine in the layer');
  assert.equal(phrase(dark, d.forms), '38 days before the March equinox');
  assert.equal(phrase(dark, d.melts), '38 days before the June solstice');
  close(d.thickest, 0.80, 0.005, 'dark: thickest'); close(C(d.waterHigh), 5.35, 0.01, 'dark: warmest water');
  close(d.summerHeat / 1e9, 1.46, 0.005, 'GJ/m² above freezing'); close(d.brine, 0.43, 0.005, 'psu of brine'); close(d.againstDeep, 0.04, 0.005, 'kg/m³ against the deep water');
});

test('the deeper the mixed layer, the more heat it holds above freezing and the later it freezes; snow keeps either ice from melting', () => {
  const slab = seaIceYears({ ice: 'dark', depth: DEPTHS.slab }), deep = seaIceYears({ ice: 'dark', depth: DEPTHS.deep });
  const s = slab.summary(), last = deep.summary();
  console.log(`dark ice over ${DEPTHS.slab.toFixed(2)} m: from ${phrase(slab, s.forms)} to ${phrase(slab, s.melts)}, water up to ${C(s.waterHigh).toFixed(1)} °C holding ${(s.summerHeat / 1e9).toFixed(2)} GJ/m²; over 200 m: none from year ${[0, 1, 2, 3, 4].findIndex((y) => deep.summary(y).none) + 1} on, water ${C(last.waterLow).toFixed(2)} to ${C(last.waterHigh).toFixed(2)} °C`);
  assert.equal(phrase(slab, s.forms), '8 days before the December solstice');
  close(s.summerHeat / 1e9, 1.00, 0.005, 'the slab\'s heat above freezing');
  for (let y = 3; y < YEARS; y++) assert.ok(deep.summary(y).none, `no ice over 200 m in year ${y + 1}`);
  assert.ok(last.summerHeat > 2e9 && last.waterLow > FREEZING_POINT, '200 m keeps its heat through the winter');
  const bright = seaIceYears({ snow: true }), dark = seaIceYears({ ice: 'dark', snow: true });
  assert.deepEqual(Array.from(dark.record.volume), Array.from(bright.record.volume), 'under snow, dark and bright ice are alike');
  let previous = 0;
  for (let y = 1; y < YEARS; y++) {
    const year = bright.summary(y);
    assert.ok(year.allYear && year.thinnest > previous, `year ${y + 1}: snowy ice never melts and thickens`);
    previous = year.thinnest;
  }
  console.log(`with snow, year ${YEARS}: ${bright.summary().thinnest.toFixed(1)}–${bright.summary().thickest.toFixed(1)} m of ice under up to ${snowDepth(bright.summary().deepestSnow).toFixed(1)} m of snow`);
});

test('the controls: the dark ice of year 10 is 0.94 m thick over 85 percent of the sea for three months, 60 W/m² builds the bright ice to 15.7 m, 92 W/m² melts it each summer and ends the dark ice after four winters, and snow piles the ice to 12 m', () => {
  const dark = seaIceYears({ ice: 'dark' }), year = dark.year(9), r = dark.record;
  let cover = 0;
  for (let d = 9 * 365; d < 10 * 365; d++) cover = Math.max(cover, r.cover[d]);
  console.log(`dark ice, year 10: ${year.thickestIce.toFixed(3)} m thick over ${(100 * cover).toFixed(1)} percent, ${dark.summary(9).days} days, water up to ${C(year.waterHigh).toFixed(2)} °C`);
  close(year.thickestIce, 0.94, 0.005, 'thickest ice'); close(cover, 0.85, 0.005, 'its cover'); assert.equal(dark.summary(9).days, 91);
  const heavy = seaIceYears({ convergence: 60 }).summary(), light = seaIceYears({ convergence: 92 }), lightDark = seaIceYears({ ice: 'dark', convergence: 92 });
  console.log(`60 W/m²: bright ice ${heavy.thinnest.toFixed(2)}–${heavy.thickest.toFixed(2)} m in year 20; 92 W/m²: bright ice-free ${light.summary().iceFreeDays} days, dark ice gone from year ${[...Array(YEARS).keys()].find((y) => lightDark.summary(y).none) + 1}`);
  close(heavy.thickest, 15.65, 0.05, 'bright ice at 60 W/m²');
  assert.ok(heavy.allYear && !light.summary().allYear && light.summary().iceFreeDays > 0, 'at 92 W/m² the bright ice melts each summer');
  assert.equal([...Array(YEARS).keys()].find((y) => lightDark.summary(y).none), 4, 'no dark ice from the fifth year');
  const snowy = seaIceYears({ snow: true }).summary();
  close(snowy.thinnest, 11.8, 0.1, 'snowy ice, thinnest'); close(snowy.thickest, 12.1, 0.1, 'snowy ice, thickest');
});
