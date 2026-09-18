import { availableParallelism } from 'node:os';
import { Worker } from 'node:worker_threads';
import { createModel } from './model.module.js';
import { shareMesh } from './mesh.module.js';

export const PHASE = { IDLE: 0, FLUX: 1, COLUMN: 2, LAYER: 3, CELL: 4, ADVANCE: 5, COMBINE: 6, ADJUST: 7, EXIT: 8 };

function split(n, index, workers) {
  return [Math.floor(n * index / workers), Math.floor(n * (index + 1) / workers)];
}

export function workerRanges(index, workers, { K, C, E, V }) {
  return {
    layers: split(K, index, workers),
    cells: split(C, index, workers),
    vertices: split(V, index, workers),
    arrays: [split(C, index, workers), split(K * C, index, workers), split(K * E, index, workers), split(C, index, workers)],
  };
}

/*
 * The model with its time step computed by worker threads. The mesh and
 * every array that crosses a phase boundary live once in shared memory;
 * each worker owns a block of layers, a block of cells, a block of
 * vertices and a block of each state array, and the main thread drives
 * the RK4 phases through a shared counter. Sums are done in the same
 * order as on one thread, so the result is bit-identical to createModel.
 */
export async function createParallelModel(grid, options = {}, workers = Math.max(1, availableParallelism() - 2)) {
  const model = createModel(grid, options);
  const { K, C, E, V } = model.core.diagnostics;
  const lengths = { pi: C, theta: K * C, u: K * E, surfaceT: C };
  const allocate = () => Object.fromEntries(Object.entries(lengths).map(([name, n]) => [name, new SharedArrayBuffer(8 * n)]));
  const trial = allocate();
  const stages = [allocate(), allocate(), allocate(), allocate()];
  const control = { ints: new SharedArrayBuffer(4 * 8), floats: new SharedArrayBuffer(8 * 4), totals: new SharedArrayBuffer(8 * 3 * workers) };
  const ctrl = new Int32Array(control.ints);
  const params = new Float64Array(control.floats);
  const totals = new Float64Array(control.totals);
  const meshShared = shareMesh(model.mesh);
  const buffers = { ...model.shared, trial, stages };
  const workerOptions = { ...options, nu4Hours: options.nu4Hours, physics: options.physics ?? true };
  delete workerOptions.buffers;

  const threads = [];
  const failures = [];
  await Promise.all(Array.from({ length: workers }, (_, index) => new Promise((resolve, reject) => {
    const thread = new Worker(new URL('./parallel.worker.js', import.meta.url), { workerData: { index, workers, meshShared, buffers, options: workerOptions, control } });
    thread.on('message', (message) => {
      if (message.type === 'ready') resolve();
      if (message.type === 'error') failures.push(`worker ${message.index} phase ${message.phase}: ${message.message}`);
    });
    thread.on('error', (error) => { failures.push(`worker ${index}: ${error && error.stack ? error.stack : error}`); reject(error); });
    thread.unref();
    threads.push(thread);
  })));

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
    run(PHASE.ADJUST);
    model.time += dt;
  };

  const serialDiagnostics = model.diagnostics;
  model.diagnostics = function diagnostics() {
    const sums = { absorbedSolar: 0, outgoingLongwave: 0, sensibleHeat: 0 };
    for (let w = 0; w < workers; w++) { sums.absorbedSolar += totals[3 * w]; sums.outgoingLongwave += totals[3 * w + 1]; sums.sensibleHeat += totals[3 * w + 2]; }
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
