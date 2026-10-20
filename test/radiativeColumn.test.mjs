import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createSigmaCore } from '../js/dynamics/sigmaCore.module.js';
import { createRadiation, DAY } from '../js/physics/radiation.module.js';
import { makeColumn, readColumn, ozoneAt, meanOzone, meanSunlight, LEVELS, K } from '../explain/radiativeColumn.module.js';
import { insolationMap } from '../explain/insolation.module.js';

const MLS = readColumn(JSON.parse(readFileSync(new URL('../explain/data/mls_bl36.json', import.meta.url))));
const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const sum = (values) => values.reduce((a, b) => a + b, 0);

test('the standard column is the midlatitude summer on the 36 bl36 layers', () => {
  assert.equal(MLS.temperature.length, K);
  assert.equal(K, 36);
  close(MLS.groundTemperature, 294.2, 0, 'ground temperature');
  close(MLS.surfacePressure, 101300, 0, 'surface pressure');
  close(MLS.ozoneColumn, 331.6, 0.05, 'ozone column, DU');
  close(sum(MLS.ozone) * 1000, MLS.ozoneColumn, 0.05, 'ozone layers sum to the column');
  for (let k = 1; k < K; k++) assert.ok(MLS.pressure[k] > MLS.pressure[k - 1] && MLS.height[k] < MLS.height[k - 1], `layer ${k} lies below layer ${k - 1}`);
});

test('every column call closes: absorbed sunlight − OLR = the layers\' terms + the ground\'s net flux', () => {
  const cases = [
    { profile: MLS, beam: 1362, mu: 1, albedo: 'sea' },
    { profile: MLS, beam: 340.5, mu: 0.25, albedo: 0.3 },
    { profile: MLS, beam: 0 },
    { profile: MLS, beam: 681, mu: 0.5, albedo: 0.85, carbonDioxide: 4, gases: { ozone: false } },
    {},
  ];
  for (const options of cases) for (const sensible of [true, false]) {
    const { budget, terms } = makeColumn(options).evaluate({ sensible });
    const layers = sum(terms.sunlight) + sum(terms.infrared) + sum(terms.ground);
    const scale = Math.max(budget.absorbed, budget.olr);
    close(budget.absorbed - budget.olr, layers + budget.groundNet, 1e-9 * scale, `closure for ${JSON.stringify({ ...options, profile: undefined, sensible })}`);
    close(sum(terms.sunlight), budget.air, 1e-9 * scale, 'the layers take the air\'s sunlight');
    close(budget.ozone + budget.oxygen + budget.carbonDioxide + budget.vapor + budget.aerosol + budget.upwardNearInfrared, budget.air, 1e-9 * scale, 'the air\'s sunlight by absorber');
    close(budget.insolation, budget.reflected + budget.absorbed, 1e-9 * Math.max(1, budget.insolation), 'sunlight is reflected or absorbed');
    close(terms.ground[K - 1], sensible ? budget.sensible : 0, 0, 'heat from the ground enters the lowest layer');
  }
});

test('the dry midlatitude summer: OLR 356.2, back radiation 127.6, greenhouse 68.6 W/m² at Ts 294.2 K', () => {
  const { budget } = makeColumn({ profile: MLS, beam: 0 });
  close(budget.groundEmission, 424.8, 0.05, 'σTs⁴');
  close(budget.olr, 356.2, 1, 'OLR');
  close(budget.dlr, 127.6, 1, 'DLR');
  close(budget.greenhouse, 68.6, 1, 'G = σTs⁴ − OLR');
  close(budget.escaping, 0.653, 0.005, 'share of the ground\'s emission reaching space');
});

test('the sun overhead on the dry midlatitude summer over sea: 7.1 % reflected, 49.1 W/m² taken by the air', () => {
  const { budget, heating, pressure } = makeColumn({ profile: MLS, beam: 1362, mu: 1, albedo: 'sea' });
  close(budget.reflected, 96.3, 1, 'reflected');
  close(budget.albedo, 0.071, 0.001, 'reflected share');
  close(budget.air, 49.1, 0.05, 'absorbed in the air');
  close(budget.ozone, 35.8, 0.5, 'by ozone');
  close(budget.ground, 1216.6, 0.05, 'by the ground');
  let peak = 0;
  for (let k = 1; k < K; k++) if (heating.sunlight[k] > heating.sunlight[peak]) peak = k;
  close(heating.sunlight[peak], 26.5, 0.5, 'peak heating, K/day');
  close(pressure[peak], 160, 5, 'peak heating pressure, Pa');
});

test('the global-mean column settles at Ts 283.3 K with 298.0 W/m² in and out', () => {
  const column = makeColumn();
  for (let n = 0; n < 1000 * 8; n++) column.step(3 * 3600);
  const { budget } = column;
  close(column.Ts, 283.3, 0.3, 'Ts');
  close(budget.absorbed, 298.0, 0.05, 'ASR');
  close(budget.olr, 298.0, 0.05, 'OLR');
  close(budget.imbalance, 0, 0.01, 'ASR − OLR');
  close(budget.albedo, 0.125, 0.002, 'planetary albedo');
  close(column.T[K - 1] - column.T[K - 8], 9.77e-3 * (column.height[K - 8] - column.height[K - 1]), 0.05, 'the lowest kilometer is dry adiabatic');
  close(sum(column.terms.convection), 0, 1e-9, 'convection moves heat without making it');
});

test('without convection or heat from the ground the column settles at Ts 284.9 K, 33.7 K/km in its lowest layers', () => {
  const column = makeColumn();
  for (let n = 0; n < 1000 * 8; n++) column.step(3 * 3600, { convection: false });
  close(column.Ts, 284.9, 0.1, 'Ts');
  close(column.budget.sensible, 0, 0, 'no heat from the ground');
  close(column.budget.imbalance, 0, 0.01, 'ASR − OLR');
  close((column.T[K - 1] - column.T[K - 8]) / (column.height[K - 8] - column.height[K - 1]) * 1000, 33.7, 0.2, 'lapse rate, K/km');
});

test('restart returns the column to its start', () => {
  const column = makeColumn({ profile: MLS });
  const T = Float64Array.from(column.T), Ts = column.Ts;
  for (let n = 0; n < 50; n++) column.step(3600, { convection: false });
  assert.ok(Math.abs(column.Ts - Ts) > 0.1);
  column.restart();
  for (let k = 0; k < K; k++) close(column.T[k], T[k], 1e-9, `T at layer ${k}`);
  assert.equal(column.Ts, Ts);
  assert.equal(column.time, 0);
});

test('ozoneAt is the model\'s own climatology at that latitude and day', () => {
  const mesh = buildMesh(new Grid(2)), core = createSigmaCore(mesh, { levels: LEVELS }), { C, exnerLayer } = core.diagnostics;
  const radiation = createRadiation(mesh, core, { longwaveLines: 'scaled' });
  const pi = new Float64Array(C).fill(101325), theta = new Float64Array(K * C), q = new Float64Array(K * C);
  for (const [i, day] of [[2, 0], [9, 120], [24, 120], [41, 250]]) {
    radiation.setTime(day * DAY);
    core.diagnoseColumn(i, pi, theta);
    for (let k = 0; k < K; k++) theta[k * C + i] = MLS.temperature[k] / exnerLayer[k * C + i];
    core.diagnoseColumn(i, pi, theta);
    radiation.column(i, 101325, theta, 290, 5, undefined, 800, null, q, null, 0.1, 0.1);
    const own = { ...radiation.budget };
    const ours = makeColumn({ temperature: MLS.temperature, groundTemperature: 290, ozoneProfile: ozoneAt(mesh.latCell[i] * 180 / Math.PI, day), beam: 800, mu: 800 / 1362, albedo: 0.1 }).budget;
    close(ours.ozone, own.ozoneSolar, 1e-9, `ozone sunlight at cell ${i}, day ${day}`);
    close(ours.olr, own.outgoingLongwave, 1e-9, `OLR at cell ${i}, day ${day}`);
  }
  close(sum(meanOzone()) * 1000, 320.6, 0.1, 'annual, area mean ozone column, DU');
  const days = new Float64Array(K);
  for (let d = 0; d < 365; d++) ozoneAt(45, d + 0.5).forEach((x, k) => { days[k] += x / 365; });
  meanOzone(45).forEach((x, k) => close(x, days[k], 1e-6 * sum(days), `annual mean ozone at 45°, layer ${k}`));
});

test('mean sunlight: 415.6 W/m² a year at the equator, 174.1 at the pole, 541.8 on the pole\'s June solstice', () => {
  close(meanSunlight(0).insolation, 415.6, 0.01, 'equator, annual');
  close(meanSunlight(45).insolation, 307.1, 0.1, '45°, annual');
  close(meanSunlight(90).insolation, 174.1, 0.1, 'pole, annual');
  close(meanSunlight(90, { from: 91, days: 1 }).insolation, 541.8, 0.1, 'pole, June solstice');
  close(meanSunlight(45, { from: 274, days: 1 }).insolation, 116.6, 0.1, '45°, December solstice');
  close(meanSunlight(0).mu, 0.753, 0.001, 'equator, insolation-weighted μ');
  const { annual } = insolationMap();
  for (const latitude of [0, 20, 45, 70, 90]) close(meanSunlight(latitude).insolation, annual[90 + latitude], 0.01, `${latitude}°, annual, against the seasons map`);
});
