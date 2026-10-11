import { test } from 'node:test';
import assert from 'node:assert/strict';
import { K, makeColumn } from '../explain/radiativeColumn.module.js';
import { createPair, tropopause, lapseRate, logScale, pickLayers, daysPerSecond, STEPS_PER_DAY, BALANCED, CALM_DAYS } from '../explain/radiativeConvective.module.js';

const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);
const sum = (values) => values.reduce((a, b) => a + b, 0);
const until = (pair, days) => { for (let n = 0; n < days * STEPS_PER_DAY && !pair.settled(); n++) pair.step(); return pair; };

const settled = until(createPair(), 1500);

test('from 288 K everywhere the pair settles within 1,000 days: 283.3 K with convection, 284.9 K with radiation alone', () => {
  const { column, ghost } = settled;
  assert.ok(settled.settled(), 'settled');
  assert.ok(settled.since < 1000, `settled on day ${settled.since}`);
  for (const c of [column, ghost]) {
    assert.ok(Math.abs(c.budget.imbalance) < BALANCED && Math.abs(c.budget.groundNet) < BALANCED, 'both columns balance');
    close(c.time, column.time, 0, 'the two columns keep the same clock');
  }
  close(column.Ts, 283.31, 0.02, 'ground with convection');
  close(ghost.Ts, 284.85, 0.02, 'ground with radiation alone');
  close(ghost.Ts - column.Ts, 1.54, 0.03, 'convection cools the ground by only 1.5 K in dry air');
  close(column.budget.absorbed, 298.0, 0.1, 'absorbed sunlight');
  close(column.budget.albedo, 0.125, 0.001, 'planetary albedo');
  close(column.budget.greenhouse, 67.3, 0.1, 'G');
});

test('convection holds the lowest 2 km on the dry adiabat; radiation alone leaves its lowest half kilometer at 33.7 K/km', () => {
  const { column, ghost } = settled;
  let top = K - 1;
  while (Math.abs(column.theta[top - 1] - column.theta[K - 1]) < 1e-6) top--;
  close(column.pressure[top] / 100, 790, 1, 'the mixed layer reaches 790 hPa');
  close(column.height[top], 1980, 30, 'about 2 km up');
  close(lapseRate(column, top), 9.77, 0.01, 'dry adiabatic, g/c_p');
  close(lapseRate(ghost, K - 8), 33.7, 0.2, 'radiation alone, lowest 8 layers');
  close(ghost.height[K - 8], 530, 10, 'which reach about 530 m');
  close(sum(column.terms.convection), 0, 1e-9, 'convection only moves heat');
  close(column.heating.ground[K - 1] + column.heating.convection[K - 1] + column.heating.infrared[K - 1] + column.heating.sunlight[K - 1], 0, 0.01, 'the lowest layer balances');
  assert.equal(sum(ghost.terms.ground), 0, 'no heat from the ground without convection');
  const tp = tropopause(column);
  close(tp.pressure / 100, 53.8, 0.1, 'tropopause pressure');
  close(tp.temperature, 208.1, 0.1, 'tropopause temperature');
});

test('doubling carbon dioxide keeps the state, cools the upper stratosphere within days and warms the ground by 1.3 K over months', () => {
  const pair = until(createPair(), 1500), before = { Ts: pair.column.Ts, top: pair.column.T[1] };
  pair.set({ carbonDioxide: 2 });
  close(pair.column.Ts, before.Ts, 0, 'set keeps the ground');
  close(pair.column.T[1], before.top, 0, 'set keeps the air');
  assert.ok(!pair.settled() && pair.since === 0, 'a change starts the clock again');
  assert.ok(pair.column.budget.imbalance > 3, 'more carbon dioxide traps infrared at once');
  let stratosphere = null, ground = null;
  for (let n = 1; n <= 2000 * STEPS_PER_DAY && !pair.settled(); n++) {
    pair.step();
    if (stratosphere === null && before.top - pair.column.T[1] > 17.1 / 2) stratosphere = n / STEPS_PER_DAY;
    if (ground === null && pair.column.Ts - before.Ts > 1.28 / 2) ground = n / STEPS_PER_DAY;
  }
  assert.ok(pair.settled(), 'settles again');
  close(pair.column.Ts - before.Ts, 1.28, 0.03, 'ground warming');
  close(pair.column.T[1] - before.top, -17.1, 0.2, 'cooling at 0.65 hPa');
  assert.ok(stratosphere < 5, `the upper stratosphere goes halfway in ${stratosphere} days`);
  assert.ok(ground > 40 && ground < 70, `the ground goes halfway in ${ground} days`);
});

test('a change of convection alone restarts the clock without touching the radiation', () => {
  const pair = createPair({ convection: false });
  assert.equal(pair.convection, false);
  assert.equal(pair.column.settings.convection, undefined, 'convection is the pair\'s, not the column\'s');
  for (let n = 0; n < 10; n++) pair.step();
  const T = Float64Array.from(pair.column.T);
  pair.set({ convection: true });
  assert.equal(pair.since, 0);
  assert.deepEqual(Float64Array.from(pair.column.T), T, 'the state is kept');
  pair.restart();
  assert.equal(pair.column.time, 0);
  close(pair.column.Ts, 288, 0, 'restart returns to the start');
  assert.ok(!pair.settled() && CALM_DAYS > 0, 'no column is settled at the start');
});

test('the tropopause is the coldest layer with warmer air both above and below it', () => {
  assert.equal(tropopause(makeColumn()), null, 'an isothermal column has none');
  const falling = Array.from({ length: K }, (_, k) => 200 + 2.5 * k);
  assert.equal(tropopause(makeColumn({ temperature: falling })), null, 'cooling all the way up has none');
  const bent = falling.map((t, k) => (k < 10 ? 225 - 2 * k : t));
  const tp = tropopause(makeColumn({ temperature: bent }));
  assert.equal(tp.k, 9);
  close(tp.temperature, 207, 1e-9, 'its temperature');
});

test('arrow lengths: zero at zero, odd, growing by equal steps per factor of ten, clipped at the largest', () => {
  assert.equal(logScale(0, 0.03, 40), 0);
  for (const v of [0.01, 0.1, 1, 10]) close(logScale(-v, 0.03, 40), -logScale(v, 0.03, 40), 1e-15, `odd at ${v}`);
  close(logScale(40, 0.03, 40), 1, 1e-15, 'the largest fills the half');
  assert.equal(logScale(400, 0.03, 40), 1);
  const steps = [0.1, 1, 10, 100].map((v) => logScale(v * 1e-3, 1e-9, 1));
  close(steps[1] - steps[0], steps[2] - steps[1], 1e-3, 'equal steps far above the small scale');
  for (let v = 0.001; v < 40; v *= 1.5) assert.ok(logScale(v * 1.5, 0.03, 40) > logScale(v, 0.03, 40), `grows at ${v}`);
});

test('arrow rows: the listed layers first, then from the ground up, none closer than the gap', () => {
  const centers = [10, 30, 40, 52, NaN, 90, 95, 100];
  const chosen = pickLayers(centers, 15, [6, 7]);
  assert.deepEqual(chosen, [0, 1, 3, 6]);
  for (const a of chosen) for (const b of chosen) if (a !== b) assert.ok(Math.abs(centers[a] - centers[b]) >= 15);
  assert.deepEqual(pickLayers([0, 5, 10], 6), [0, 2]);
});

test('the run starts at 25 days a second and speeds up as the column settles', () => {
  assert.equal(daysPerSecond(0), 25);
  assert.equal(daysPerSecond(30), 25);
  assert.equal(daysPerSecond(300), 200);
});
