import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Worker } from 'node:worker_threads';
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

/*
 * Wall time (ms) of `threads` threads each doing the same fixed burn at
 * once, so that threads × burn(1) / burn(threads) is the speedup the
 * machine can give `threads` workers at this moment.
 */
function burn(threads, iterations = 3e7) {
  const source = "const { parentPort, workerData } = require('node:worker_threads'); let x = 0; for (let i = 0; i < workerData; i++) x += Math.sqrt(i); parentPort.postMessage(x);";
  const t0 = performance.now();
  return Promise.all(Array.from({ length: threads }, () => new Promise((resolve, reject) => {
    const worker = new Worker(source, { eval: true, workerData: iterations });
    worker.once('message', resolve);
    worker.once('error', reject);
  }))).then(() => performance.now() - t0);
}

test('worker threads speed up an N=16 atmosphere step by a fair share of the cores free right now', async (t) => {
  const N = +(process.env.PARALLEL_TEST_N ?? 16);
  const serial = createModel(new Grid(N), { ocean: false });
  const parallel = await createParallelModel(new Grid(N), { ocean: false });
  try {
    const headroom = parallel.workers * (await burn(1)) / (await burn(parallel.workers));
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
    const speedup = serialMs / parallelMs;
    console.log(`N=${N}: serial ${serialMs.toFixed(0)} ms/step, ${parallel.workers} workers ${parallelMs.toFixed(0)} ms/step, speedup ${speedup.toFixed(1)}× with ${headroom.toFixed(1)}× of headroom`);
    for (let a = 0; a < serial.state.length; a++) assert.deepEqual(parallel.state[a], serial.state[a]);
    if (headroom < 1.5) { t.skip(`only ${headroom.toFixed(1)}× of headroom for ${parallel.workers} workers: the cores are busy`); return; }
    assert.ok(speedup >= 0.3 * headroom, `speedup ${speedup.toFixed(2)}× against ${headroom.toFixed(1)}× of headroom`);
  } finally {
    await parallel.close();
  }
});
