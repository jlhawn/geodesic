// Minimal reproduction: Apple M5 (Metal) miscompiles a down-counting loop over a private array of more than 32 floats
// when Tint's robustness clamps its indices and the loop bounds become constants only after inlining.
// node repro_min.mjs            -> robustness on (the default, as in Chrome) is wrong; disable_robustness is right
// env K (array length, default 34), M (rows used, default 15), FORM=down|up, DUMP=1 prints the generated MSL
import { createRequire } from 'node:module';
const require = createRequire('/Users/jlhawn/git_repos/jlhawn/geodesic/package.json');
const { create, globals } = await import(require.resolve('webgpu'));
Object.assign(globalThis, globals);
const K = +(process.env.K ?? 34), M = +(process.env.M ?? 15), N = 4096, FORM = process.env.FORM ?? 'down', f = Math.fround;
const readBack = FORM === 'up'
  ? 'for (var q = 0; q < m; q++) { let j = m - 1 - q; OUT[j * N + i] = gain[j]; }'
  : 'for (var j = m - 1; j >= 0; j--) { OUT[j * N + i] = gain[j]; }';
const code = `const K: i32 = ${K}; const M: i32 = ${M}; const N: i32 = ${N};
@group(0) @binding(0) var<storage, read> L: array<f32>;
@group(0) @binding(1) var<storage, read_write> OUT: array<f32>;
fn solve(i: i32, m: i32) {
  var gain: array<f32, K>;
  gain[0] = L[i];
  for (var j = 1; j < m; j++) { gain[j] = L[j * N + i] * gain[j - 1]; }
  ${readBack}
}
@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= N) { return; }
  solve(i, M);
}`;
let seed = 12345; const rand = () => { seed = (seed * 1664525 + 1013904223) >>> 0; return seed / 4294967296; };
const L = Float32Array.from({ length: M * N }, () => f(0.5 + rand()));
const expected = new Float32Array(M * N);
for (let i = 0; i < N; i++) { let g = L[i]; expected[i] = g; for (let j = 1; j < M; j++) { g = f(L[j * N + i] * g); expected[j * N + i] = g; } }
const keep = [];
for (const flags of [[], ['enable-dawn-features=disable_robustness']]) {
  const features = [...flags.map((x) => x.split('=')[1]), ...(process.env.DUMP ? ['dump_shaders', 'disable_symbol_renaming'] : [])];
  const gpu = create(features.length ? [`enable-dawn-features=${features.join(',')}`] : []); keep.push(gpu);
  const adapter = await gpu.requestAdapter(), device = await adapter.requestDevice(); keep.push(device);
  const bL = device.createBuffer({ size: L.byteLength, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST }); device.queue.writeBuffer(bL, 0, L);
  const bO = device.createBuffer({ size: L.byteLength, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC });
  const read = device.createBuffer({ size: L.byteLength, usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ });
  const pipeline = device.createComputePipeline({ layout: 'auto', compute: { module: device.createShaderModule({ code }), entryPoint: 'main' } });
  const bind = device.createBindGroup({ layout: pipeline.getBindGroupLayout(0), entries: [bL, bO].map((buffer, binding) => ({ binding, resource: { buffer } })) });
  const enc = device.createCommandEncoder(), pass = enc.beginComputePass();
  pass.setPipeline(pipeline); pass.setBindGroup(0, bind); pass.dispatchWorkgroups(N / 64); pass.end();
  enc.copyBufferToBuffer(bO, 0, read, 0, L.byteLength); device.queue.submit([enc.finish()]);
  await read.mapAsync(GPUMapMode.READ); const out = new Float32Array(read.getMappedRange().slice(0));
  let wrong = 0, first = -1;
  for (let x = 0; x < out.length; x++) if (out[x] !== expected[x]) { wrong++; if (first < 0) first = x; }
  const where = first < 0 ? '' : `; first at thread ${first % N}, j ${Math.floor(first / N)}: ${out[first]} for ${expected[first]}`;
  console.log(`${adapter.info.vendor} ${adapter.info.architecture} K=${K} M=${M} ${FORM} robustness ${flags.length ? 'off' : 'on '}: ${wrong} of ${out.length} wrong${where}`);
}
process.exit(0);
