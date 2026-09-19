import { createModel, STATE_NAMES, stateLengths } from './model.module.js';
import { shareMesh } from './mesh.module.js';
import { parallelism, spawn } from './threads.module.js';

export const PHASE = { IDLE: 0, FLUX: 1, COLUMN: 2, LAYER: 3, CELL: 4, ADVANCE: 5, COMBINE: 6, ADJUST: 7, EXIT: 8 };

function split(n, index, workers) {
  return [Math.floor(n * index / workers), Math.floor(n * (index + 1) / workers)];
}

export function workerRanges(index, workers, { K, C, E, V }) {
  const lengths = stateLengths({ K, C, E });
  return {
    layers: split(K, index, workers),
    cells: split(C, index, workers),
    vertices: split(V, index, workers),
    arrays: STATE_NAMES.map((name) => split(lengths[name], index, workers)),
  };
}

export const TOTALS = ['absorbedSolar', 'outgoingLongwave', 'sensibleHeat', 'evaporation'];

/*
 * The model with its time step computed by worker threads. The mesh and
 * every array that crosses a phase boundary live once in shared memory;
 * each worker owns a block of layers, a block of cells, a block of
 * vertices and a block of each state array, and the main thread drives
 * the RK4 phases through a shared counter. Sums are done in the same
 * order as on one thread, so the result is bit-identical to createModel.
 */
export async function createParallelModel(grid, options = {}, workers = null) {
  workers ??= Math.max(1, await parallelism() - 2);
  const model = createModel(grid, options);
  const { K, C, E } = model.core.diagnostics;
  const lengths = stateLengths({ K, C, E });
  const allocate = () => Object.fromEntries(Object.entries(lengths).map(([name, n]) => [name, new SharedArrayBuffer(8 * n)]));
  const trial = allocate();
  const stages = [allocate(), allocate(), allocate(), allocate()];
  const control = { ints: new SharedArrayBuffer(4 * 8), floats: new SharedArrayBuffer(8 * 4), totals: new SharedArrayBuffer(8 * TOTALS.length * workers) };
  const ctrl = new Int32Array(control.ints);
  const params = new Float64Array(control.floats);
  const totals = new Float64Array(control.totals);
  const meshShared = shareMesh(model.mesh);
  const buffers = { ...model.shared, trial, stages };
  const workerOptions = { ...options, nu4Hours: options.nu4Hours, physics: options.physics ?? true };
  delete workerOptions.buffers;

  const threads = [];
  const failures = [];
  const ready = [];
  for (let index = 0; index < workers; index++) {
    let resolveReady, rejectReady;
    ready.push(new Promise((resolve, reject) => { resolveReady = resolve; rejectReady = reject; }));
    threads.push(await spawn(new URL('./parallel.worker.js', import.meta.url), { index, workers, meshShared, buffers, options: workerOptions, control }, {
      onMessage(message) {
        if (message.type === 'ready') resolveReady();
        if (message.type === 'error') failures.push(`worker ${message.index} phase ${message.phase}: ${message.message}`);
      },
      onError(error) { failures.push(`worker ${index}: ${error && error.stack ? error.stack : error}`); rejectReady(new Error(failures.join('\n'))); },
    }));
  }
  await Promise.all(ready);
  for (const thread of threads) thread.unref();

  function run(phase, { useTrial = 0, stage = 0, dt = 0, factor = 0 } = {}) {
    Atomics.store(ctrl, 3, useTrial);
    Atomics.store(ctrl, 4, stage);
    params[0] = dt;
    params[1] = factor;
    params[2] = model.time;
    Atomics.store(ctrl, 2, 0);
    Atomics.store(ctrl, 1, phase);
    Atomics.add(ctrl, 0, 1);
    Atomics.notify(ctrl, 0);
    let done;
    while ((done = Atomics.load(ctrl, 2)) < workers) {
      if (Atomics.wait(ctrl, 2, done, 30000) === 'timed-out') throw new Error(`phase ${phase}: ${done} of ${workers} workers finished after 30 s${failures.length ? '\n' + failures.join('\n') : ''}`);
    }
    if (Atomics.load(ctrl, 5)) {
      const detail = failures.join('\n') || 'a worker reported an error';
      throw new Error(`phase ${phase} failed:\n${detail}`);
    }
  }

  function tendencyPhases(useTrial, stage) {
    run(PHASE.FLUX, { useTrial, stage });
    run(PHASE.COLUMN, { useTrial, stage });
    run(PHASE.LAYER, { useTrial, stage });
    run(PHASE.CELL, { useTrial, stage });
  }

  model.step = function step(dt) {
    tendencyPhases(0, 0);
    run(PHASE.ADVANCE, { stage: 0, factor: dt / 2 });
    tendencyPhases(1, 1);
    run(PHASE.ADVANCE, { stage: 1, factor: dt / 2 });
    tendencyPhases(1, 2);
    run(PHASE.ADVANCE, { stage: 2, factor: dt });
    tendencyPhases(1, 3);
    run(PHASE.COMBINE, { dt });
    run(PHASE.ADJUST, { dt });
    model.time += dt;
  };

  const serialDiagnostics = model.diagnostics;
  model.diagnostics = function diagnostics() {
    const sums = Object.fromEntries(TOTALS.map((name) => [name, 0]));
    for (let w = 0; w < workers; w++) TOTALS.forEach((name, t) => { sums[name] += totals[TOTALS.length * w + t]; });
    return serialDiagnostics(sums);
  };

  model.workers = workers;
  model.close = async function close() {
    Atomics.store(ctrl, 1, PHASE.EXIT);
    Atomics.add(ctrl, 0, 1);
    Atomics.notify(ctrl, 0);
    await Promise.all(threads.map((thread) => thread.terminate()));
  };
  return model;
}
