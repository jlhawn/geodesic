import { test } from 'node:test';
import assert from 'node:assert/strict';
import { K, makeColumn, manabeWetherald } from '../explain/radiativeColumn.module.js';
import { tropopause, STEP, STEPS_PER_DAY } from '../explain/radiativeConvective.module.js';
import { createWetColumn, convectionTop, troposphericLapse, nearestLayer, parcelHumidity, moistAdiabat, SUBSTEPS } from '../explain/moistConvective.module.js';
import { saturationHumidity } from '../js/physics/moist.module.js';
import { KAPPA } from '../explain/physics.module.js';

const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const until = (run, days = 4000) => { for (let n = 0; n < days * STEPS_PER_DAY && !run.settled(); n++) run.step(); assert.ok(run.settled(), 'the run settles'); return run; };

const wet = until(createWetColumn({ albedo: 0.3 }));

test('the wet column settles on day 1,507 at 286.6 K, with convection up to 318 hPa and the cold point at 102 hPa', () => {
  const { column } = wet, { budget } = column;
  close(column.time / 86400, 1507.1, 0.2, 'settling day');
  close(column.Ts, 286.64, 0.01, 'ground');
  close(budget.absorbed, 249.63, 0.05, 'absorbed sunlight');
  close(budget.olr, 249.63, 0.05, 'OLR');
  close(budget.albedo, 0.267, 0.001, 'planetary albedo');
  close(budget.greenhouse, 133.16, 0.05, 'G');
  close(budget.precipitableWater, 12.50, 0.01, 'precipitable water, kg/m²');
  close(budget.evaporation * 86400, 2.07, 0.005, 'rain = evaporation, mm/day');
  close(budget.latent, 59.89, 0.05, 'L·E');
  close(budget.sensible, 17.03, 0.05, 'sensible heat');
  const top = convectionTop(column), cold = tropopause(column);
  close(top.pressure / 100, 317.7, 0.1, 'top of convection, hPa');
  close(top.height, 8541, 5, 'top of convection, m');
  close(cold.pressure / 100, 101.8, 0.1, 'cold point, hPa');
  close(cold.temperature, 202.66, 0.02, 'cold point, K');
  close(column.pressure[nearestLayer(column, 850e2)] / 100, 850.5, 0.1, 'the layer nearest 850 hPa');
  close(column.pressure[nearestLayer(column, 300e2)] / 100, 317.7, 0.1, 'the layer nearest 300 hPa');
  close(troposphericLapse(column), 8.12, 0.005, 'lapse rate from 850 to 300 hPa, K/km');
});

test('the moist adiabat guide is the curve the convecting layers lie on, and in dry air it is the dry adiabat', () => {
  const { column } = wet, top = column.convectionTop, up = Array.from({ length: K }, (_, n) => column.pressure[K - 1 - n]);
  close(parcelHumidity(column), manabeWetherald(column.pressure[K - 1] / 101325) * saturationHumidity(column.T[K - 1], column.pressure[K - 1]), 1e-15, 'the parcel carries the lowest layer\'s relative humidity');
  const guide = moistAdiabat(column, up);
  for (let k = top; k < K; k++) close(guide[K - 1 - k], column.T[k], 1e-9, `layer ${k}`);
  assert.ok(column.T[top - 1] - guide[K - top] > 1, 'the layer above the top of convection is warmer than the adiabat');
  const dry = makeColumn({ albedo: 0.3, temperature: column.T, groundTemperature: column.Ts });
  const plain = moistAdiabat(dry, up);
  for (let n = 0; n < K; n++) if (up[n] >= 100e2) close(plain[n], dry.T[K - 1] * (up[n] / up[0]) ** KAPPA, 0.01, `dry air follows the dry adiabat at ${up[n].toFixed(0)} Pa`);
});

test('a change keeps the last settled state; frozen vapor is a copy, and frozen with radiation alone settles at 308.5 K', () => {
  const run = createWetColumn({ albedo: 0.3 });
  run.column.set({ temperature: Float64Array.from(wet.column.T), groundTemperature: wet.column.Ts });
  run.set({ vapor: 'frozen', convection: 'off' });
  assert.equal(run.before, null, 'no settled state yet');
  const held = Float64Array.from(run.column.q);
  for (let k = 0; k < K; k++) close(held[k], wet.column.q[k], 1e-15, `layer ${k} is frozen at its amount`);
  until(run);
  for (let k = 0; k < K; k++) assert.equal(run.column.q[k], held[k], `layer ${k} stays frozen`);
  close(run.column.Ts, 308.52, 0.01, 'radiation alone with the vapor frozen');
  close(run.column.Ts - wet.column.Ts, 21.88, 0.01, 'warmer than with moist convection');
  assert.equal(run.column.convectionTop, -1, 'no convection');
  close(run.last.Ts, run.column.Ts, 0, 'the settled state is kept');
  run.set({ convection: 'dry' });
  close(run.before.Ts, 308.52, 0.01, 'the change keeps it as the state before');
  assert.equal(run.last, null, 'and waits for the next');
  assert.ok(!run.settled(), 'a change unsettles the run');
  const once = makeColumn({ albedo: 0.3, humidity: held, temperature: wet.column.T, groundTemperature: wet.column.Ts });
  for (let n = 0; n < 20000; n++) once.step(STEP, { convection: false });
  close(once.Ts, 308.52, 0.01, 'the 1-hour substeps settle where 3-hour steps do');
});

test('radiation alone takes 1-hour substeps: frozen with a hot column\'s 66 kg/m² of vapor, 3-hour steps run away within days', () => {
  const warm = Float64Array.from(wet.column.T, (T) => T + 25), q = Float64Array.from(warm, (T, k) => (wet.column.pressure[k] < 100e2 ? 3e-6 : manabeWetherald(wet.column.pressure[k] / 101325) * saturationHumidity(T, wet.column.pressure[k])));
  const start = { albedo: 0.07, temperature: warm, groundTemperature: wet.column.Ts + 25 };
  const direct = makeColumn({ ...start, humidity: q });
  assert.ok(direct.budget.precipitableWater > 60, `${direct.budget.precipitableWater} kg/m²`);
  let ran = false;
  for (let n = 0; n < 10 * STEPS_PER_DAY && !ran; n++) { direct.step(STEP, { convection: false }); ran = !(Math.abs(direct.T[K - 1] - 300) < 100); }
  assert.ok(ran, 'a 3-hour step runs away');
  const run = createWetColumn({ ...start, vapor: 'none', convection: 'off' });
  run.column.set({ humidity: q });
  let swing = 0, previous = run.column.T[K - 1], change = 0;
  for (let n = 0; n < 60 * STEPS_PER_DAY; n++) {
    run.step();
    const next = run.column.T[K - 1] - previous;
    if (n > 8 && Math.sign(next) * Math.sign(change) < 0 && Math.abs(next) > 0.5) swing++;
    change = next; previous = run.column.T[K - 1];
  }
  assert.equal(SUBSTEPS, 3);
  assert.ok(Number.isFinite(run.column.Ts) && run.column.T.every(Number.isFinite), 'the substeps stay finite');
  assert.equal(swing, 0, 'and the lowest layer does not oscillate');
});
