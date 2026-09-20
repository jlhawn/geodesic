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
  const device = await adapter.requestDevice({ requiredLimits });
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
