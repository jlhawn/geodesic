import { emptyBuffer, readRanges, createKernel } from './device.module.js';
import { SST_REFERENCE, sstMeanOf } from '../sstMeans.module.js';

/*
 * The GPU side of ../sstMeans.module.js for a model of model.gpu.js:
 * after each model step that ended an ocean step, one dispatch adds each
 * ocean cell's OD.T0, the top-layer temperature the ocean's writeSurface
 * left, less SST_REFERENCE into a sum of its own. `step` records through
 * the model's core, so inside model.stepBatch it belongs in the per-step
 * callback; `day` and `checkpoint` belong outside a batch, after the
 * day's steps are queued, and `day` queues its read and the clearing of
 * the sum before it returns, so the next day's steps may follow at once.
 */
export function createGpuSstRecorder(model, from = null) {
  const { gpu, oceanEngine: ocean } = model, { device } = gpu, C = model.mesh.nCells, { OD } = ocean.layout;
  const sum = emptyBuffer(device, 4 * C);
  if (from) device.queue.writeBuffer(sum, 0, Float32Array.from(from.sum));
  const kernel = createKernel(device, `
@group(0) @binding(0) var<storage, read_write> OD: array<f32>;
@group(0) @binding(1) var<storage, read_write> SUM: array<f32>;
@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x);
  if (i >= ${C} || OD[${OD.CMASK} + i] < 0.5) { return; }
  SUM[i] += OD[${OD.T0} + i] - ${SST_REFERENCE};
}`, 'sstSum');
  let counted = 0, samples = from ? from.samples : 0;
  return {
    step() {
      if (++counted % ocean.everySteps !== 0) return;
      gpu.compute((pass) => kernel.run(pass, [ocean.buffers.OD, sum], C));
      samples++;
    },
    day() {
      const read = readRanges(device, sum, [{ offset: 0, length: C }]), taken = samples;
      gpu.clearBuffer(sum, 0, 4 * C);
      samples = 0;
      return read.then(([sums]) => ({ mean: sstMeanOf(sums, taken, ocean.cellOcean), samples: taken }));
    },
    async checkpoint() {
      const [sums] = await readRanges(device, sum, [{ offset: 0, length: C }]);
      return { sum: Float64Array.from(sums), samples };
    },
    destroy() { sum.destroy(); },
  };
}
