import { getDevice } from './device.module.js';

/*
 * A profile of the GPU model, for finding what makes a device slow. Runs
 * `steps` steps one at a time, each timed from its submission until the
 * device has finished it. Where the device offers 'timestamp-query', GPU
 * timestamps around every compute pass are charged to the kernels the
 * pass dispatched; a pass whose two timestamps are equal (browsers may
 * round them, Chrome to 0.1 ms) counts as unresolved rather than as
 * time. With `split` every dispatch is recorded into a timed pass of its
 * own, the pipeline and bind groups set before it carried over, so that
 * each kernel's time is known alone: the dispatches run in the same order
 * on the same buffers, and only the passes' boundaries differ. Each row
 * carries its dispatches' resolved times, sorted, with their median and
 * least, which tell a kernel's own time from time another process's work
 * took on a shared device. Also measures an empty round trip to the
 * device. The model advances by those steps, so the caller stops its own
 * loop first.
 */
const QUERIES = 4096, MAX_PASSES = 32 * QUERIES / 2;

export async function profileGpu(model, { steps = 16, dt, split = false }) {
  const { device } = model.gpu, { adapter } = await getDevice();
  const settle = () => device.queue.onSubmittedWorkDone();
  await settle();
  const trips = [];
  for (let n = 0; n < 20; n++) { const start = performance.now(); await settle(); trips.push(performance.now() - start); }

  const timed = device.features.has('timestamp-query');
  const passes = [], profiled = new WeakMap(), undo = [], querySets = [];
  device.pushErrorScope('validation');
  const timestampWrites = (index) => {
    const set = Math.floor(2 * index / QUERIES), slot = 2 * index % QUERIES;
    if (set === querySets.length) querySets.push(device.createQuerySet({ type: 'timestamp', count: QUERIES }));
    return { querySet: querySets[set], beginningOfPassWriteIndex: slot, endOfPassWriteIndex: slot + 1 };
  };
  if (timed) {
    const begin = GPUCommandEncoder.prototype.beginComputePass, setPipeline = GPUComputePassEncoder.prototype.setPipeline;
    const timedPass = (encoder, descriptor, record) => {
      const index = passes.push(record) - 1;
      return begin.call(encoder, { ...descriptor, timestampWrites: timestampWrites(index) });
    };
    GPUCommandEncoder.prototype.beginComputePass = function (descriptor = {}) {
      if (descriptor.timestampWrites) return begin.call(this, descriptor);
      if (split) return splitPass(this, descriptor, timedPass, () => passes.length < MAX_PASSES);
      if (passes.length >= MAX_PASSES) return begin.call(this, descriptor);
      const record = { kernels: [] };
      const pass = timedPass(this, descriptor, record);
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
    const bytes = 8 * QUERIES * querySets.length;
    const resolved = device.createBuffer({ size: bytes, usage: GPUBufferUsage.QUERY_RESOLVE | GPUBufferUsage.COPY_SRC });
    const read = device.createBuffer({ size: bytes, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST });
    const encoder = device.createCommandEncoder();
    querySets.forEach((set, n) => encoder.resolveQuerySet(set, 0, Math.min(QUERIES, 2 * passes.length - n * QUERIES), resolved, 8 * QUERIES * n));
    encoder.copyBufferToBuffer(resolved, 0, read, 0, bytes);
    device.queue.submit([encoder.finish()]);
    await read.mapAsync(GPUMapMode.READ);
    const stamps = new BigUint64Array(read.getMappedRange().slice(0));
    read.unmap();
    read.destroy(); resolved.destroy();
    const totals = new Map();
    let total = 0;
    passes.forEach(({ kernels: names, groups, within }, i) => {
      const resolved = stamps[2 * i + 1] > stamps[2 * i], ms = resolved ? Number(stamps[2 * i + 1] - stamps[2 * i]) / 1e6 : 0;
      const label = [...new Set(names)].join(' + ') || 'no kernel';
      const key = within ? `${label}|${groups}|${within.kernels.join(' + ')}` : label;
      const entry = totals.get(key) ?? totals.set(key, { name: label, ms: 0, passes: 0, unresolved: 0, times: [], ...(within ? { groups, within: [...new Set(within.kernels)].join(' + ') } : {}) }).get(key);
      entry.ms += ms / steps; entry.passes += 1 / steps; total += ms;
      if (resolved) entry.times.push(ms); else entry.unresolved += 1 / steps;
    });
    kernels = [...totals.values()].sort((a, b) => b.ms - a.ms);
    for (const entry of kernels) {
      const times = entry.times.sort((a, b) => a - b);
      entry.median = times.length ? times[times.length >> 1] : 0; entry.least = times.length ? times[0] : 0;
    }
    gpuMs = total / steps;
  }
  for (const set of querySets) set.destroy();

  const sorted = (values) => [...values].sort((a, b) => a - b), median = (values) => sorted(values)[values.length >> 1];
  const info = adapter.info ?? {};
  return {
    device: [info.vendor, info.architecture, info.device, info.description].filter(Boolean).join(' / ') || 'unknown',
    timestamps: timed, split, steps, stepMedian: median(stepTimes), stepMin: sorted(stepTimes)[0], stepMax: sorted(stepTimes)[steps - 1],
    roundTrip: median(trips), gpuMs, kernels,
  };
}

/*
 * What the engine records into one compute pass, recorded instead as one
 * timed pass per dispatch: the pipeline and the bind groups last set are
 * set again on each new pass. `within` names the kernels the engine's own
 * pass held, for telling apart one kernel's dispatches from different
 * passes.
 */
function splitPass(encoder, descriptor, timedPass, room) {
  const within = { kernels: [] }, bound = [];
  let pipeline = null;
  const run = (groups, dispatch) => {
    const label = pipeline?.label || 'unnamed';
    within.kernels.push(label);
    if (!room()) throw new Error('the split profile ran out of timestamp queries; profile fewer steps');
    const pass = timedPass(encoder, descriptor, { kernels: [label], groups, within });
    pass.setPipeline(pipeline);
    bound.forEach((args, index) => args && pass.setBindGroup(index, ...args));
    dispatch(pass);
    pass.end();
  };
  return {
    label: descriptor.label ?? '',
    setPipeline(next) { pipeline = next; },
    setBindGroup(index, ...args) { bound[index] = args; },
    dispatchWorkgroups(x, y = 1, z = 1) { run(x * y * z, (pass) => pass.dispatchWorkgroups(x, y, z)); },
    dispatchWorkgroupsIndirect(buffer, offset) { run(0, (pass) => pass.dispatchWorkgroupsIndirect(buffer, offset)); },
    pushDebugGroup() {}, popDebugGroup() {}, insertDebugMarker() {},
    end() {},
  };
}
