// Probe: with eight storage bindings of varied sizes, does arrayLength() report the true length,
// and do clamped reads of each buffer's last element (direct, and inside a helper function) return it?
import { createRequire } from 'node:module';
const require = createRequire('/Users/jlhawn/git_repos/jlhawn/geodesic/package.json');
const { create, globals } = await import(require.resolve('webgpu'));
Object.assign(globalThis, globals);
const flags = (process.env.FLAGS ?? '').split(' ').filter(Boolean);
const gpu = create(flags); const adapter = await gpu.requestAdapter();
const device = await adapter.requestDevice({ requiredLimits: { maxStorageBuffersPerShaderStage: adapter.limits.maxStorageBuffersPerShaderStage, maxStorageBufferBindingSize: adapter.limits.maxStorageBufferBindingSize, maxBufferSize: adapter.limits.maxBufferSize } });
const sizes = [1000, 3000, 500, 20000, 20000, 200000, 8, 1500000];
const decl = sizes.map((_, b) => `@group(0) @binding(${b}) var<storage, read> B${b}: array<f32>;`).join('\n');
const code = `${decl}
@group(0) @binding(8) var<storage, read_write> OUT: array<f32>;
${sizes.map((_, b) => `fn last${b}() -> f32 { return B${b}[i32(arrayLength(&B${b})) - 1]; }`).join('\n')}
@compute @workgroup_size(1) fn main() {
${sizes.map((_, b) => `  OUT[${b}] = f32(arrayLength(&B${b})); OUT[${8 + b}] = B${b}[i32(arrayLength(&B${b})) - 1]; OUT[${16 + b}] = last${b}();`).join('\n')}
}`;
const bufs = sizes.map((n, b) => { const a = new Float32Array(n); a[n - 1] = 1000 + b; const buf = device.createBuffer({ size: a.byteLength, usage: 0x88 }); device.queue.writeBuffer(buf, 0, a); return buf; });
const out = device.createBuffer({ size: 4 * 24, usage: 0x84 }), rd = device.createBuffer({ size: 4 * 24, usage: 0x9 });
const module = device.createShaderModule({ code }); for (const m of (await module.getCompilationInfo()).messages) console.log(m.type, m.message);
const pipeline = device.createComputePipeline({ layout: 'auto', compute: { module, entryPoint: 'main' } });
const bind = device.createBindGroup({ layout: pipeline.getBindGroupLayout(0), entries: [...bufs, out].map((buffer, binding) => ({ binding, resource: { buffer } })) });
const enc = device.createCommandEncoder(); const pass = enc.beginComputePass(); pass.setPipeline(pipeline); pass.setBindGroup(0, bind); pass.dispatchWorkgroups(1); pass.end();
enc.copyBufferToBuffer(out, 0, rd, 0, 4 * 24); device.queue.submit([enc.finish()]);
await rd.mapAsync(1); const o = new Float32Array(rd.getMappedRange().slice(0));
console.log('binding: true length / arrayLength / last element direct / last element via helper');
sizes.forEach((n, b) => console.log(`  B${b}: ${n} / ${o[b]} / ${o[8 + b]} / ${o[16 + b]}  ${o[b] === n && o[8 + b] === 1000 + b && o[16 + b] === 1000 + b ? 'ok' : 'WRONG'}`));
process.exit(0);
