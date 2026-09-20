import { createModel, STATE_NAMES, stateLengths } from './model.module.js';
import { shareMesh } from './mesh.module.js';
import { parallelism, spawn } from './threads.module.js';

export const PHASE = { IDLE: 0, FLUX: 1, COLUMN: 2, LAYER: 3, PHYSICS: 4, ADVANCE: 5, COMBINE: 6, CLOSURE: 7, ADJUST: 8, EXIT: 9 };

function blocks(kind, n, size, extra = {}) {
  const chunks = [];
  for (let from = 0; from < n; from += size) chunks.push({ kind, from, to: Math.min(n, from + size), ...extra });
  return chunks;
}

/*
 * The work units of each phase. Workers claim units through a shared
 * counter, so the split adapts to however fast each core is: a layer at
 * a time for the layer-partitioned phases, and for the rest blocks of
 * cells, vertices or array elements sized to give every worker several
 * units whatever the resolution.
 */
export function phaseChunks({ K, C, E, V }, workers = 8) {
  const lengths = stateLengths({ K, C, E });
  const size = (n, floor) => Math.max(floor, Math.ceil(n / (4 * workers)));
  const arrays = STATE_NAMES.flatMap((name, a) => blocks('array', lengths[name], size(lengths[name], 4096), { a }));
  return {
    [PHASE.FLUX]: blocks('layers', K, 1),
    [PHASE.COLUMN]: [...blocks('cells', C, size(C, 32)), ...blocks('vertices', V, size(V, 64))],
    [PHASE.LAYER]: [...blocks('momentum', K, 1), ...blocks('tracers', K, 1)],
    [PHASE.PHYSICS]: blocks('cells', C, size(C, 32)),
    [PHASE.ADVANCE]: arrays,
    [PHASE.COMBINE]: arrays,
    [PHASE.CLOSURE]: [...blocks('momentum', K, 1), ...blocks('tracers', K, 1)],
    [PHASE.ADJUST]: [...blocks('cells', C, size(C, 16)), ...blocks('edges', E, size(E, 64))],
  };
}

export const TOTALS = ['absorbedSolar', 'outgoingLongwave', 'sensibleHeat', 'evaporation', 'insolation', 'reflectedSolar'];

/*
 * The model with its time step computed by worker threads. The mesh and
 * every array that crosses a phase boundary live once in shared memory;
 * the main thread drives the RK4 phases through a shared generation
 * counter and the workers claim each phase's work units from a shared
 * chunk counter. Every array element is computed by exactly one worker
 * with the single-thread arithmetic, so the state is bit-identical to
 * createModel's; only the radiation totals are summed in a different
 * order.
 */
export async function createParallelModel(grid, options = {}, workers = null) {
  workers ??= Math.max(1, await parallelism());
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
  const workerOptions = { ...options, nu4Hours: options.nu4Hours, physics: options.physics ?? true, ocean: false };
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

  const phaseTime = new Float64Array(PHASE.EXIT + 1);
  function run(phase, { useTrial = 0, stage = 0, dt = 0, factor = 0 } = {}) {
    const started = performance.now();
    Atomics.store(ctrl, 3, useTrial);
    Atomics.store(ctrl, 4, stage);
    params[0] = dt;
    params[1] = factor;
    params[2] = model.time;
    Atomics.store(ctrl, 6, 0);
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
    phaseTime[phase] += performance.now() - started;
  }

  function tendencyPhases(useTrial, stage) {
    run(PHASE.FLUX, { useTrial, stage });
    run(PHASE.COLUMN, { useTrial, stage });
    run(PHASE.LAYER, { useTrial, stage });
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
    model.phases.ocean(dt);
    run(PHASE.PHYSICS, { dt });
    run(PHASE.CLOSURE, { dt });
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
  model.phaseTime = phaseTime;
  model.close = async function close() {
    Atomics.store(ctrl, 1, PHASE.EXIT);
    Atomics.add(ctrl, 0, 1);
    Atomics.notify(ctrl, 0);
    await Promise.all(threads.map((thread) => thread.terminate()));
  };
  return model;
}
