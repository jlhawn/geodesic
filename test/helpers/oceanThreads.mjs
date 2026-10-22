/*
 * An oceanCases simulation with its ocean spread over worker threads. Each
 * worker builds a twin of the ocean from sim.mesh and sim.oceanOptions on
 * the ocean's shared arrays, and the ocean hands every pass (the layer
 * tendencies and the regrid, column, eddy and barotropic passes) to the
 * workers and the calling thread in blocks of layers, cells or edges, as
 * createParallelModel's workers take the coupled model's ocean. Every array
 * element is computed by one thread with the single-thread arithmetic, so a
 * threaded run is bit-identical to runCase on one thread. OCEAN_THREADS sets
 * the number of threads, the calling one included (default the core count,
 * at most 8); 1 leaves the sim as it is.
 */
import { Worker } from 'node:worker_threads';
import { availableParallelism } from 'node:os';
import { shareMesh } from '../../js/mesh.module.js';
import { PHASE, OCEAN_PASSES, OCEAN_PASS_FUNCTIONS, phaseChunks } from '../../js/parallel.module.js';
import { runCase } from '../../explain/oceanCases.module.js';

export const SLOT = { generation: 0, phase: 1, done: 2, trial: 3, stage: 4, claimed: 5 };
const PATIENCE_MS = 60000, FAILURE_BYTES = 4096;

export function oceanThreads() {
  const asked = Math.floor(Number(process.env.OCEAN_THREADS));
  return asked >= 1 ? asked : Math.min(8, availableParallelism());
}

export const passChunks = (mesh, ocean, threads) => phaseChunks({ K: 0, C: mesh.nCells, E: mesh.nEdges, V: mesh.nVertices, L: ocean.layers }, threads);

export function claim(ctrl, ocean, chunks) {
  const phase = Atomics.load(ctrl, SLOT.phase), stage = Atomics.load(ctrl, SLOT.stage);
  const units = chunks[phase], input = Atomics.load(ctrl, SLOT.trial) ? ocean.trial : ocean.state;
  for (let c; (c = Atomics.add(ctrl, SLOT.claimed, 1)) < units.length;) {
    const { from, to, kind } = units[c];
    if (phase === PHASE.OCEAN) ocean.tendencyLayers(input, ocean.stages[stage], from, to, kind);
    else if (phase === PHASE.OCEAN_PREPARE_CELLS) ocean.prepareCells(from, to, input);
    else if (phase === PHASE.OCEAN_PREPARE_EDGES) ocean.prepareEdges(from, to, input);
    else ocean[OCEAN_PASS_FUNCTIONS[phase]](from, to);
  }
}

export function report(failure, error) {
  const length = new Int32Array(failure, 0, 1);
  if (Atomics.compareExchange(length, 0, 0, -1) !== 0) return;
  const text = new TextEncoder().encode(String(error && error.stack ? error.stack : error)).subarray(0, FAILURE_BYTES);
  new Uint8Array(failure, 4).set(text);
  Atomics.store(length, 0, text.length);
}

export async function threaded(sim, threads = oceanThreads()) {
  const { mesh, ocean } = sim;
  sim.close = async () => {};
  if (threads < 2) return sim;
  const workers = threads - 1, control = new SharedArrayBuffer(4 * 6), failure = new SharedArrayBuffer(4 + FAILURE_BYTES);
  const ctrl = new Int32Array(control), failed = new Int32Array(failure, 0, 1);
  const oceanOptions = Object.fromEntries(Object.entries(sim.oceanOptions).filter(([, value]) => typeof value !== 'function'));
  const workerData = { mesh: shareMesh(mesh), oceanOptions, buffers: ocean.shared, control, failure, threads };
  const pool = [];
  try {
    await Promise.all(Array.from({ length: workers }, () => new Promise((resolve, reject) => {
      const worker = new Worker(new URL('./oceanThreads.worker.mjs', import.meta.url), { workerData });
      pool.push(worker);
      worker.once('message', resolve).once('error', reject).once('exit', (code) => reject(new Error(`an ocean thread exited with code ${code} before it was ready`)));
    })));
  } catch (error) {
    await Promise.all(pool.map((worker) => worker.terminate()));
    throw error;
  }
  for (const worker of pool) worker.unref();
  const chunks = passChunks(mesh, ocean, threads);

  function run(phase, useTrial = 0, stage = 0) {
    Atomics.store(ctrl, SLOT.phase, phase);
    Atomics.store(ctrl, SLOT.trial, useTrial);
    Atomics.store(ctrl, SLOT.stage, stage);
    Atomics.store(ctrl, SLOT.claimed, 0);
    Atomics.store(ctrl, SLOT.done, 0);
    Atomics.add(ctrl, SLOT.generation, 1);
    Atomics.notify(ctrl, SLOT.generation);
    let own = null;
    try { claim(ctrl, ocean, chunks); } catch (error) { own = error; }
    for (let done; (done = Atomics.load(ctrl, SLOT.done)) < workers;) {
      if (Atomics.wait(ctrl, SLOT.done, done, PATIENCE_MS) === 'timed-out') throw new Error(`${done} of ${workers} ocean threads finished pass ${phase} within ${PATIENCE_MS / 1000} s`);
    }
    if (own) throw own;
    const length = Atomics.load(failed, 0);
    if (length !== 0) throw new Error(`an ocean thread failed in pass ${phase}:\n${new TextDecoder().decode(new Uint8Array(failure, 4, Math.max(0, length)).slice())}`);
  }
  ocean.setLayerRunner((useTrial, stage) => run(PHASE.OCEAN, useTrial, stage));
  ocean.setColumnRunner((pass, useTrial = 0) => run(OCEAN_PASSES[pass], useTrial));
  sim.close = () => {
    ocean.setLayerRunner(null);
    ocean.setColumnRunner(null);
    Atomics.store(ctrl, SLOT.phase, PHASE.EXIT);
    Atomics.add(ctrl, SLOT.generation, 1);
    Atomics.notify(ctrl, SLOT.generation);
    return Promise.all(pool.map((worker) => worker.terminate()));
  };
  return sim;
}

export async function runThreaded(sim, days) {
  await threaded(sim);
  try {
    return runCase(sim, days);
  } finally {
    await sim.close();
  }
}
