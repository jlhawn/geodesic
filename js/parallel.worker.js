import { parentPort, workerData } from 'node:worker_threads';
import { meshFromShared } from './mesh.module.js';
import { createModel } from './model.module.js';
import { PHASE, workerRanges } from './parallel.module.js';

const { index, workers, meshShared, buffers, options, control } = workerData;
const mesh = meshFromShared(meshShared);
const model = createModel(mesh, { ...options, buffers });
const { K, C, E, V } = model.core.diagnostics;
const ctrl = new Int32Array(control.ints);
const params = new Float64Array(control.floats);
const totals = new Float64Array(control.totals);
const state = model.state;
const trial = ['pi', 'theta', 'u', 'surfaceT'].map((name) => new Float64Array(buffers.trial[name]));
const stages = buffers.stages.map((stage) => ['pi', 'theta', 'u', 'surfaceT'].map((name) => new Float64Array(stage[name])));
const ranges = workerRanges(index, workers, { K, C, E, V });
const sums = { absorbedSolar: 0, outgoingLongwave: 0, sensibleHeat: 0 };

function run(phase) {
  const input = Atomics.load(ctrl, 3) ? trial : state;
  const out = stages[Atomics.load(ctrl, 4)];
  switch (phase) {
    case PHASE.FLUX: model.phases.flux(input, ranges.layers[0], ranges.layers[1]); break;
    case PHASE.COLUMN:
      model.phases.column(input, out, ranges.cells[0], ranges.cells[1]);
      model.phases.vertex(input, ranges.vertices[0], ranges.vertices[1]);
      break;
    case PHASE.LAYER: model.phases.layer(input, out, ranges.layers[0], ranges.layers[1]); break;
    case PHASE.CELL:
      model.radiation.setTime(params[2]);
      model.phases.cell(input, out, ranges.cells[0], ranges.cells[1], sums);
      totals[3 * index] = sums.absorbedSolar;
      totals[3 * index + 1] = sums.outgoingLongwave;
      totals[3 * index + 2] = sums.sensibleHeat;
      break;
    case PHASE.ADVANCE: {
      const factor = params[1];
      for (let a = 0; a < 4; a++) {
        const [from, to] = ranges.arrays[a];
        const s = state[a], k = out[a], t = trial[a];
        for (let i = from; i < to; i++) t[i] = s[i] + factor * k[i];
      }
      break;
    }
    case PHASE.COMBINE: {
      const w = params[0] / 6;
      for (let a = 0; a < 4; a++) {
        const [from, to] = ranges.arrays[a];
        const s = state[a], k1 = stages[0][a], k2 = stages[1][a], k3 = stages[2][a], k4 = stages[3][a];
        for (let i = from; i < to; i++) s[i] += w * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
      }
      break;
    }
    case PHASE.ADJUST: model.phases.adjust(ranges.cells[0], ranges.cells[1]); break;
    default: break;
  }
}

parentPort.postMessage({ type: 'ready', index });
let generation = 0;
for (;;) {
  Atomics.wait(ctrl, 0, generation);
  generation = Atomics.load(ctrl, 0);
  const phase = Atomics.load(ctrl, 1);
  if (phase === PHASE.EXIT) break;
  try {
    run(phase);
  } catch (error) {
    Atomics.store(ctrl, 5, 1);
    parentPort.postMessage({ type: 'error', index, phase, message: error && error.stack ? error.stack : String(error) });
  }
  Atomics.add(ctrl, 2, 1);
  Atomics.notify(ctrl, 2);
}
