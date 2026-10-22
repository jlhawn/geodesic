import { workerData, parentPort } from 'node:worker_threads';
import { meshFromShared } from '../../js/mesh.module.js';
import { createOcean } from '../../js/ocean/layered.module.js';
import { PHASE } from '../../js/parallel.module.js';
import { SLOT, passChunks, claim, report } from './oceanThreads.mjs';

const { mesh: shared, oceanOptions, buffers, control, failure, threads } = workerData;
const mesh = meshFromShared(shared);
const ocean = createOcean(mesh, { ...oceanOptions, buffers });
const ctrl = new Int32Array(control), chunks = passChunks(mesh, ocean, threads);
parentPort.postMessage('ready');
for (let generation = 0; ;) {
  // the notify that announces a generation can wake a thread that has already done it
  while (Atomics.load(ctrl, SLOT.generation) === generation) Atomics.wait(ctrl, SLOT.generation, generation);
  generation = Atomics.load(ctrl, SLOT.generation);
  if (Atomics.load(ctrl, SLOT.phase) === PHASE.EXIT) break;
  try { claim(ctrl, ocean, chunks); } catch (error) { report(failure, error); }
  Atomics.add(ctrl, SLOT.done, 1);
  Atomics.notify(ctrl, SLOT.done);
}
