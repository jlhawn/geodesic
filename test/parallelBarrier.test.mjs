import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createParallelModel, PHASE } from '../js/parallel.module.js';
import { initializeState } from '../js/physics/init.module.js';

const WORKERS = 8, ROUNDS = 2000, LATE_MS = 0.05, SETTLE_MS = 20;
const GENERATION = 0, PHASE_SLOT = 1, DONE = 2, CLAIMED = 6;

function spin(ms) { for (const end = performance.now() + ms; performance.now() < end;); }

function controlOf(model) {
  let ctrl = null;
  const { notify } = Atomics;
  Atomics.notify = (array, index, count) => { ctrl ??= array; return notify(array, index, count); };
  try {
    model.step(600);
  } finally {
    Atomics.notify = notify;
  }
  return ctrl;
}

test('each worker passes each phase barrier once when the main thread announces the next phase at once and wakes the workers late', async () => {
  const parallel = await createParallelModel(new Grid(16), {}, WORKERS);
  const init = initializeState(parallel, {});
  for (let a = 0; a < 4; a++) parallel.state[a].set(init[a]);
  let passes = 0;
  try {
    const ctrl = controlOf(parallel);
    spin(SETTLE_MS);
    Atomics.store(ctrl, DONE, 0);
    for (let round = 0; round < ROUNDS; round++) {
      Atomics.store(ctrl, CLAIMED, 0);
      passes += Atomics.exchange(ctrl, DONE, 0);
      Atomics.store(ctrl, PHASE_SLOT, PHASE.IDLE);
      Atomics.add(ctrl, GENERATION, 1);
      spin(LATE_MS);
      Atomics.notify(ctrl, GENERATION);
      const deadline = performance.now() + 10000;
      while (Atomics.load(ctrl, DONE) < WORKERS) if (performance.now() > deadline) assert.fail(`round ${round}: ${Atomics.load(ctrl, DONE)} of ${WORKERS} workers passed`);
    }
    spin(SETTLE_MS);
    passes += Atomics.exchange(ctrl, DONE, 0);
  } finally {
    await parallel.close();
  }
  assert.equal(passes, WORKERS * ROUNDS, `${passes - WORKERS * ROUNDS} passes beyond one per worker in ${ROUNDS} rounds`);
});
