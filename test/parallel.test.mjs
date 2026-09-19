import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createParallelModel } from '../js/parallel.module.js';
import { initializeState } from '../js/physics/init.module.js';

test('the worker-thread step reproduces the single-thread step bit for bit', async () => {
  const serial = createModel(new Grid(8));
  const parallel = await createParallelModel(new Grid(8), {}, 3);
  const init = initializeState(serial, {});
  for (let a = 0; a < 4; a++) { serial.state[a].set(init[a]); parallel.state[a].set(init[a]); }
  try {
    for (let n = 0; n < 4; n++) { serial.step(600); parallel.step(600); }
    for (let a = 0; a < serial.state.length; a++) {
      const s = serial.state[a], p = parallel.state[a];
      for (let i = 0; i < s.length; i++) if (s[i] !== p[i]) assert.fail(`array ${a} differs at ${i}: ${s[i]} vs ${p[i]}`);
    }
    const ds = serial.diagnostics(), dp = parallel.diagnostics();
    assert.equal(ds.mass, dp.mass);
    assert.ok(Math.abs(ds.absorbedSolar - dp.absorbedSolar) < 1e-9 * ds.absorbedSolar);
    assert.ok(Math.abs(ds.outgoingLongwave - dp.outgoingLongwave) < 1e-9 * ds.outgoingLongwave);
    assert.equal(serial.time, parallel.time);
  } finally {
    await parallel.close();
  }
});

test('worker threads speed up an N=16 step', async () => {
  const N = +(process.env.PARALLEL_TEST_N ?? 16);
  const serial = createModel(new Grid(N));
  const parallel = await createParallelModel(new Grid(N));
  const init = initializeState(serial, {});
  for (let a = 0; a < 4; a++) { serial.state[a].set(init[a]); parallel.state[a].set(init[a]); }
  serial.step(450); parallel.step(450);
  const steps = 6;
  let t0 = performance.now();
  for (let n = 0; n < steps; n++) serial.step(450);
  const serialMs = (performance.now() - t0) / steps;
  t0 = performance.now();
  for (let n = 0; n < steps; n++) parallel.step(450);
  const parallelMs = (performance.now() - t0) / steps;
  console.log(`N=${N}: serial ${serialMs.toFixed(0)} ms/step, ${parallel.workers} workers ${parallelMs.toFixed(0)} ms/step, speedup ${(serialMs / parallelMs).toFixed(1)}×`);
  try {
    for (let a = 0; a < serial.state.length; a++) assert.deepEqual(parallel.state[a], serial.state[a]);
    assert.ok(parallelMs < serialMs);
  } finally {
    await parallel.close();
  }
});
