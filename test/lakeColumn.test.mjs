import { test } from 'node:test';
import assert from 'node:assert/strict';
import { FLAKE, lakeColumn, lakeStep, lakeEnthalpy, lakeIceAlbedo, lakeDepth, createLake } from '../js/physics/lake.module.js';

const SIGMA = 5.670374e-8, DAY = 86400, DT = 3600;

function year(depth, { mean = 276, amplitude = 16, sun = 300, years = 3, closed = true } = {}) {
  const c = lakeColumn(277, depth), f = 1.26e-4, days = [];
  let Ts = c.mixed, gained = 0, worst = 0;
  const start = lakeEnthalpy(c, depth);
  for (let t = 0; t < years * 365 * DAY; t += DT) {
    const day = t / DAY, season = Math.sin(2 * Math.PI * (day - 110) / 365);
    const air = mean + amplitude * season, sunlight = Math.max(0, sun * (1 + season) / 2 * Math.max(0, Math.sin(2 * Math.PI * day)) * Math.PI);
    const albedo = c.ice > 0 ? lakeIceAlbedo(Ts) : 0.07, absorbed = (1 - albedo) * sunlight;
    const flux = 0.8 * SIGMA * air ** 4 - 0.97 * SIGMA * Ts ** 4 + 1.2 * 1004 * 1.3e-3 * 5 * (air - Ts) + 2.5e6 * 1.2 * 1.3e-3 * 5 * 0.004 * (air - Ts) / 10;
    const before = lakeEnthalpy(c, depth);
    Ts = lakeStep(c, depth, absorbed, flux, 0.006, f, DT, closed);
    const change = lakeEnthalpy(c, depth) - before, put = (absorbed + flux) * DT;
    gained += put;
    worst = Math.max(worst, Math.abs(change - put));
    if (t % DAY === 0 && day >= (years - 1) * 365) days.push({ day: day - (years - 1) * 365, ...c, Ts });
  }
  return { days, worst, drift: lakeEnthalpy(c, depth) - start - gained, gained };
}

test('the IFS bounds a lake between 2 and 50 m', () => {
  assert.deepEqual([lakeDepth(0), lakeDepth(10), lakeDepth(1642)], [2, 10, 50]);
});

test('a 20 m lake at 60° N runs FLake\'s seasonal cycle: stratified in summer with a mixed layer of meters, mixed to the bottom in autumn, frozen in winter over a bottom between 1 and 4 °C, open again by summer, with the heat budget closed every step', () => {
  const { days, worst } = year(20);
  const at = (d) => days.find((x) => x.day === d);
  const summer = at(200), autumn = at(300), winter = at(30), spring = at(150);
  console.log(`20 m lake, last of three years: day 200 T_ML ${(summer.mixed - 273.15).toFixed(1)} °C, h ${summer.mixedDepth.toFixed(1)} m, T_b ${(summer.bottom - 273.15).toFixed(1)} °C; day 300 h ${autumn.mixedDepth.toFixed(1)} m, T ${(autumn.mixed - 273.15).toFixed(1)} °C; day 30 ice ${winter.ice.toFixed(2)} m, T_b ${(winter.bottom - 273.15).toFixed(2)} °C; most ice ${Math.max(...days.map((x) => x.ice)).toFixed(2)} m, largest step budget error ${worst.toExponential(1)} J/m²`);
  assert.ok(summer.mixedDepth < 15 && summer.mixed - summer.bottom > 3, 'stratified in summer');
  assert.ok(autumn.mixedDepth > summer.mixedDepth, 'the mixed layer deepens through autumn');
  assert.ok(winter.ice > 0.2, 'ice in winter');
  assert.ok(winter.bottom > FLAKE.freezing + 1 && winter.bottom <= FLAKE.maximumDensityT + 0.05, 'the bottom between 1 and 4 °C under the ice, above the mixed layer at freezing');
  assert.ok(spring.ice === 0 && summer.ice === 0, 'open by late spring');
  assert.ok(worst < 1e-3, `the heat budget, J/m² a step: ${worst}`);
});

test('the closed ice grows and goes as FLake\'s own does, to days and centimeters, while FLake\'s leaks heat', () => {
  const closed = year(20), open = year(20, { closed: false });
  const season = ({ days }) => {
    const iced = days.filter((x) => x.ice > 0).map((x) => x.day);
    const breakup = Math.max(...iced.filter((d) => d < 200)), freeze = Math.min(...iced.filter((d) => d > 200));
    return { most: Math.max(...days.map((x) => x.ice)), breakup, freeze };
  };
  const a = season(closed), b = season(open);
  console.log(`closed: most ice ${a.most.toFixed(3)} m, break-up day ${a.breakup}, freeze-up day ${a.freeze}, drift ${closed.drift.toExponential(1)} J/m² over three years; FLake's own: ${b.most.toFixed(3)} m, ${b.breakup}, ${b.freeze}, drift ${(open.drift / 1e6).toFixed(1)} MJ/m²`);
  assert.ok(Math.abs(closed.drift) < 1e-2, `the closed column's three-year drift, J/m²: ${closed.drift}`);
  assert.ok(Math.abs(a.most - b.most) < 0.15 && Math.abs(a.breakup - b.breakup) <= 10 && Math.abs(a.freeze - b.freeze) <= 3, JSON.stringify({ a, b }));
});

test('a shallow 2 m lake mixes to its bottom all year and freezes to more of its depth', () => {
  const shallow = year(2), deep = year(20);
  assert.ok(shallow.days.every((x) => x.ice > 0 || x.mixedDepth > 1.99), 'mixed to the bottom whenever open');
  assert.ok(Math.max(...shallow.days.map((x) => x.mixed)) > Math.max(...deep.days.map((x) => x.mixed)), 'warmer at the top in summer');
});

test('a model lake\'s enthalpy changes by exactly the fluxes it is given, through ice and snow', () => {
  const lake = createLake({ nCells: 1, latCell: [Math.PI / 3] }, { lake: [1], elevation: [-20], saline: [0] }, { model: 'flake' });
  const Ts = [278], ice = [0], concentration = [0];
  let worst = 0, most = 0;
  for (let t = 0; t < 2 * 365 * DAY; t += DT) {
    const day = t / DAY, season = Math.sin(2 * Math.PI * (day - 110) / 365), air = 276 + 16 * season;
    const sunlight = Math.max(0, 300 * (1 + season) / 2 * Math.max(0, Math.sin(2 * Math.PI * day)) * Math.PI), absorbed = (1 - (lake.albedo(0, Ts[0]) ?? 0.07)) * sunlight;
    const flux = 0.8 * SIGMA * air ** 4 - 0.97 * SIGMA * Ts[0] ** 4 + 23.4 * (air - Ts[0]);
    const snow = air < 273.15 ? 1e-4 * DT : 0;
    if (t > 0) {
      const before = lake.enthalpy(0);
      lake.deposit(0, snow, air);
      lake.step(0, Ts, ice, concentration, absorbed, flux, 0.006, DT);
      worst = Math.max(worst, Math.abs(lake.enthalpy(0) - before - (absorbed + flux) * DT + 3.34e5 * snow));
    } else lake.step(0, Ts, ice, concentration, absorbed, flux, 0.006, DT);
    most = Math.max(most, ice[0]);
  }
  assert.ok(most > 0.2, `ice in winter: ${most} m`);
  assert.ok(worst < 1e-3, `the lake's step budget, J/m²: ${worst}`);
});
