// Probe: a storage buffer indexed inside helper functions (which receive Tint's array-length struct by value)
// against direct indexing in main. DECL=runtime (array<f32>) | fixed (array<f32, N>). FLAGS passes Dawn toggles.
import { createRequire } from 'node:module';
const require = createRequire('/Users/jlhawn/git_repos/jlhawn/geodesic/package.json');
const { create, globals } = await import(require.resolve('webgpu'));
Object.assign(globalThis, globals);
const flags = (process.env.FLAGS ?? '').split(' ').filter(Boolean);
const gpu = create(flags); const adapter = await gpu.requestAdapter(); const device = await adapter.requestDevice();
const N = 300000, DECL = process.env.DECL ?? 'runtime';
const arr = DECL === 'fixed' ? `array<f32, ${N}>` : 'array<f32>';
const code = `
@group(0) @binding(0) var<storage, read> IN: ${arr};
@group(0) @binding(1) var<storage, read_write> OUT: array<f32>;
fn deep(idx: i32, k: i32) -> f32 { return IN[idx + k * 1000] * 0.5 + IN[idx]; }
fn mid(idx: i32, k: i32) -> f32 { var acc = 0.0; for (var j = 0; j < 3; j++) { acc += deep(idx, k + j); } return acc; }
@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= ${N} - 4000) { return; }
  var direct = 0.0;
  for (var j = 0; j < 3; j++) { direct += IN[i + (1 + j) * 1000] * 0.5 + IN[i]; }
  OUT[i] = mid(i, 1) - direct;
}`;
const input = new Float32Array(N); for (let i = 0; i < N; i++) input[i] = Math.fround(Math.sin(i * 0.001) * 100 + i * 1e-3);
const bI = device.createBuffer({ size: input.byteLength, usage: 0x88 }); device.queue.writeBuffer(bI, 0, input);
const bO = device.createBuffer({ size: input.byteLength, usage: 0x84 }); const rd = device.createBuffer({ size: input.byteLength, usage: 0x9 });
const module = device.createShaderModule({ code }); for (const m of (await module.getCompilationInfo()).messages) console.log(m.type, m.message);
const pipeline = device.createComputePipeline({ layout: 'auto', compute: { module, entryPoint: 'main' } });
const bind = device.createBindGroup({ layout: pipeline.getBindGroupLayout(0), entries: [{ binding: 0, resource: { buffer: bI } }, { binding: 1, resource: { buffer: bO } }] });
const enc = device.createCommandEncoder(); const pass = enc.beginComputePass(); pass.setPipeline(pipeline); pass.setBindGroup(0, bind); pass.dispatchWorkgroups(Math.ceil(N / 64)); pass.end();
enc.copyBufferToBuffer(bO, 0, rd, 0, input.byteLength); device.queue.submit([enc.finish()]);
await rd.mapAsync(1); const out = new Float32Array(rd.getMappedRange().slice(0));
let bad = 0, first = -1, worst = 0; for (let i = 0; i < N - 4000; i++) { const d = Math.abs(out[i]); if (d > 1e-3) { bad++; if (first < 0) first = i; } worst = Math.max(worst, d); }
console.log(`DECL=${DECL} FLAGS=${flags.join(' ') || '(none)'}: helper-function reads differ from direct reads in ${bad} of ${N - 4000} elements (first at ${first}, worst ${worst.toExponential(2)})`);
process.exit(0);
