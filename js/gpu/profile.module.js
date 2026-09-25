import { getDevice } from './device.module.js';

/*
 * A profile of the GPU model, for finding what makes a device slow. Runs
 * `steps` steps one at a time, each timed from its submission until the
 * device has finished it. Where the device offers 'timestamp-query', GPU
 * timestamps around every compute pass are charged to the kernels the
 * pass dispatched; a pass whose two timestamps are equal (browsers may
 * round them, Chrome to 0.1 ms) counts as unresolved rather than as
 * time. Also measures an empty round trip to the device. The
 * model advances by those steps, so the caller stops its own loop first.
 */
const MAX_PASSES = 2048;

export async function profileGpu(model, { steps = 16, dt }) {
  const { device } = model.gpu, { adapter } = await getDevice();
  const settle = () => device.queue.onSubmittedWorkDone();
  await settle();
  const trips = [];
  for (let n = 0; n < 20; n++) { const start = performance.now(); await settle(); trips.push(performance.now() - start); }

  const timed = device.features.has('timestamp-query');
  const passes = [], profiled = new WeakMap(), undo = [];
  device.pushErrorScope('validation');
  const querySet = timed ? device.createQuerySet({ type: 'timestamp', count: 2 * MAX_PASSES }) : null;
  if (timed) {
    const begin = GPUCommandEncoder.prototype.beginComputePass, setPipeline = GPUComputePassEncoder.prototype.setPipeline;
    GPUCommandEncoder.prototype.beginComputePass = function (descriptor = {}) {
      if (passes.length >= MAX_PASSES || descriptor.timestampWrites) return begin.call(this, descriptor);
      const record = { kernels: [] }, index = passes.push(record) - 1;
      const pass = begin.call(this, { ...descriptor, timestampWrites: { querySet, beginningOfPassWriteIndex: 2 * index, endOfPassWriteIndex: 2 * index + 1 } });
      profiled.set(pass, record);
      return pass;
    };
    GPUComputePassEncoder.prototype.setPipeline = function (pipeline) {
      profiled.get(this)?.kernels.push(pipeline.label || 'unnamed');
      return setPipeline.call(this, pipeline);
    };
    undo.push(() => { GPUCommandEncoder.prototype.beginComputePass = begin; GPUComputePassEncoder.prototype.setPipeline = setPipeline; });
  }
  const stepTimes = [];
  let failure = null;
  try {
    for (let n = 0; n < steps; n++) {
      const start = performance.now();
      await model.step(dt);
      await settle();
      stepTimes.push(performance.now() - start);
    }
  } catch (error) {
    failure = error;
  } finally {
    for (const restore of undo) restore();
  }
  const invalid = await device.popErrorScope();
  if (failure) throw failure;
  if (invalid) throw new Error(`the profiled steps failed validation: ${invalid.message}`);

  let kernels = null, gpuMs = null;
  if (timed && passes.length) {
    const bytes = 16 * passes.length;
    const resolved = device.createBuffer({ size: bytes, usage: GPUBufferUsage.QUERY_RESOLVE | GPUBufferUsage.COPY_SRC });
    const read = device.createBuffer({ size: bytes, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST });
    const encoder = device.createCommandEncoder();
    encoder.resolveQuerySet(querySet, 0, 2 * passes.length, resolved, 0);
    encoder.copyBufferToBuffer(resolved, 0, read, 0, bytes);
    device.queue.submit([encoder.finish()]);
    await read.mapAsync(GPUMapMode.READ);
    const stamps = new BigUint64Array(read.getMappedRange().slice(0));
    read.unmap();
    read.destroy(); resolved.destroy(); querySet.destroy();
    const totals = new Map();
    let total = 0;
    passes.forEach(({ kernels: names }, i) => {
      const resolved = stamps[2 * i + 1] > stamps[2 * i], ms = resolved ? Number(stamps[2 * i + 1] - stamps[2 * i]) / 1e6 : 0;
      const label = [...new Set(names)].join(' + ') || 'no kernel';
      const entry = totals.get(label) ?? totals.set(label, { name: label, ms: 0, passes: 0, unresolved: 0 }).get(label);
      entry.ms += ms / steps; entry.passes += 1 / steps; total += ms;
      if (!resolved) entry.unresolved += 1 / steps;
    });
    kernels = [...totals.values()].sort((a, b) => b.ms - a.ms);
    gpuMs = total / steps;
  } else if (querySet) querySet.destroy();

  const sorted = (values) => [...values].sort((a, b) => a - b), median = (values) => sorted(values)[values.length >> 1];
  const info = adapter.info ?? {};
  return {
    device: [info.vendor, info.architecture, info.device, info.description].filter(Boolean).join(' / ') || 'unknown',
    timestamps: timed, steps, stepMedian: median(stepTimes), stepMin: sorted(stepTimes)[0], stepMax: sorted(stepTimes)[steps - 1],
    roundTrip: median(trips), gpuMs, kernels,
  };
}
