/*
 * WebGPU device access shared by Node (Google's Dawn through the
 * `webgpu` package) and the browser (navigator.gpu), with the small set
 * of buffer and pipeline helpers the engine needs. Everything on the
 * GPU is single precision: the engine uploads Float64Array state as
 * f32 and reads it back the same way.
 */
let cached = null;

export async function getDevice() {
  if (cached) return cached;
  let gpu = globalThis.navigator && globalThis.navigator.gpu;
  if (!gpu) {
    const { create, globals } = await import('webgpu');
    Object.assign(globalThis, globals);
    gpu = create([]);
  }
  const adapter = await gpu.requestAdapter();
  if (!adapter) throw new Error('no WebGPU adapter');
  const wanted = ['maxStorageBuffersPerShaderStage', 'maxStorageBufferBindingSize', 'maxBufferSize', 'maxComputeWorkgroupsPerDimension'];
  const requiredLimits = {};
  for (const name of wanted) if (adapter.limits[name] !== undefined) requiredLimits[name] = adapter.limits[name];
  const requiredFeatures = adapter.features.has('timestamp-query') ? ['timestamp-query'] : [];
  const device = await adapter.requestDevice({ requiredLimits, requiredFeatures });
  cached = { gpu, adapter, device };
  return cached;
}

export function storageBuffer(device, data, extraUsage = 0) {
  const buffer = device.createBuffer({ size: Math.max(4, Math.ceil(data.byteLength / 4) * 4), usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST | GPUBufferUsage.COPY_SRC | extraUsage, mappedAtCreation: true });
  new (data.constructor)(buffer.getMappedRange()).set(data);
  buffer.unmap();
  return buffer;
}

export function emptyBuffer(device, byteLength) {
  return device.createBuffer({ size: Math.max(4, Math.ceil(byteLength / 4) * 4), usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST | GPUBufferUsage.COPY_SRC });
}

export async function readBuffer(device, buffer, byteLength, Type = Float32Array) {
  const staging = device.createBuffer({ size: byteLength, usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ });
  const encoder = device.createCommandEncoder();
  encoder.copyBufferToBuffer(buffer, 0, staging, 0, byteLength);
  device.queue.submit([encoder.finish()]);
  await staging.mapAsync(GPUMapMode.READ);
  const out = new Type(staging.getMappedRange().slice(0));
  staging.unmap();
  staging.destroy();
  return out;
}

/*
 * Several ranges of one buffer, given in floats, read back through one
 * staging buffer and one map.
 */
export function readRanges(device, buffer, ranges) {
  const total = ranges.reduce((sum, r) => sum + r.length, 0);
  const staging = device.createBuffer({ size: 4 * total, usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ });
  const encoder = device.createCommandEncoder();
  let at = 0;
  for (const r of ranges) { encoder.copyBufferToBuffer(buffer, 4 * r.offset, staging, 4 * at, 4 * r.length); at += r.length; }
  device.queue.submit([encoder.finish()]);
  return staging.mapAsync(GPUMapMode.READ).then(() => {
    const all = new Float32Array(staging.getMappedRange().slice(0));
    staging.unmap();
    staging.destroy();
    let from = 0;
    return ranges.map((r) => { const view = all.subarray(from, from + r.length); from += r.length; return view; });
  });
}

/*
 * A compute pipeline with bind groups cached per buffer set: `run`
 * records one dispatch of `count` invocations at workgroup size 64.
 */
export function createKernel(device, code, label = 'kernel') {
  const module = device.createShaderModule({ code, label });
  const pipeline = device.createComputePipeline({ label, layout: 'auto', compute: { module, entryPoint: 'main' } });
  const groups = new Map();
  function bindGroup(buffers) {
    const key = buffers.map((b) => b.label ?? b.__id ?? (b.__id = Math.random())).join('|');
    let group = groups.get(key);
    if (!group) {
      group = device.createBindGroup({ layout: pipeline.getBindGroupLayout(0), entries: buffers.map((buffer, binding) => ({ binding, resource: { buffer } })) });
      groups.set(key, group);
    }
    return group;
  }
  function run(pass, buffers, count) {
    pass.setPipeline(pipeline);
    pass.setBindGroup(0, bindGroup(buffers));
    pass.dispatchWorkgroups(Math.ceil(count / 64));
  }
  return { run, pipeline };
}

export const REDUCTION_WORKGROUP = 64;
export const reductionGroups = (count) => Math.ceil(count / REDUCTION_WORKGROUP);

/*
 * A reduction over `count` items in one dispatch: each invocation
 * evaluates one expression per quantity for its item after `setup`,
 * each workgroup combines them by sum, min or max in shared memory, and
 * invocation 0 writes the partials to OUT[base + q·groups + workgroup],
 * with workgroups past 65535 in the dispatch's second dimension.
 * finishReduction combines the partials on the host in double precision.
 */
export function reductionKernel(quantities, { count, base, setup = '' }) {
  const W = REDUCTION_WORKGROUP, groups = reductionGroups(count);
  const identity = (kind) => (kind === 'min' ? '3.0e38' : kind === 'max' ? '-3.0e38' : '0.0');
  const combine = (kind, a, b) => (kind === 'sum' ? `${a} + ${b}` : `${kind}(${a}, ${b})`);
  const cell = (q) => `acc[${q * W} + slot]`;
  return `var<workgroup> acc: array<f32, ${quantities.length * W}>;
@compute @workgroup_size(${W}) fn main(@builtin(workgroup_id) wg: vec3<u32>, @builtin(local_invocation_index) li: u32) {
  let slot = i32(li);
  let wgIndex = i32(wg.x) + i32(wg.y) * 65535;
  let i = wgIndex * ${W} + slot;
${quantities.map(([, kind], q) => `  ${cell(q)} = ${identity(kind)};`).join('\n')}
  if (i < ${count}) {
${setup}
${quantities.map(([, , expr], q) => `    ${cell(q)} = ${expr};`).join('\n')}
  }
  workgroupBarrier();
  for (var s = ${W / 2}; s > 0; s = s / 2) {
    if (slot < s) {
${quantities.map(([, kind], q) => `      ${cell(q)} = ${combine(kind, cell(q), `acc[${q * W} + slot + s]`)};`).join('\n')}
    }
    workgroupBarrier();
  }
  if (slot == 0 && wgIndex < ${groups}) {
${quantities.map((_, q) => `    OUT[${base} + ${q * groups} + wgIndex] = acc[${q * W}];`).join('\n')}
  }
}`;
}

export function finishReduction(quantities, partials, count) {
  const groups = reductionGroups(count), out = {};
  quantities.forEach(([name, kind], q) => {
    let value = kind === 'min' ? Infinity : kind === 'max' ? -Infinity : 0;
    for (let w = 0; w < groups; w++) {
      const x = partials[q * groups + w];
      value = kind === 'sum' ? value + x : kind === 'min' ? Math.min(value, x) : Math.max(value, x);
    }
    out[name] = value;
  });
  return out;
}
