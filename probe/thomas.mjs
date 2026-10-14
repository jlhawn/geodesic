// Probe: the physics kernel's thomas() solve on this GPU against a float32 JS solve.
// FORM=ptr (verbatim pointer-to-array form) | inline (solve written inline on local arrays) | fwd (backward sweep as a forward loop)
import { getDevice } from '/Users/jlhawn/git_repos/jlhawn/geodesic/js/gpu/device.module.js';
const K = 36, KTOP = 2, n = K - KTOP, N = 8192, FORM = process.env.FORM ?? 'ptr';
const f = Math.fround;
const { device } = await getDevice();
const body = {
  ptr: `
fn thomas(n: i32, upper: ptr<function, array<f32, K>>, lower: ptr<function, array<f32, K>>, rhs: ptr<function, array<f32, K>>) {
  var gain: array<f32, K>;
  var denominator = 1.0 + (*upper)[0] + (*lower)[0];
  gain[0] = -(*lower)[0] / denominator;
  (*rhs)[0] = (*rhs)[0] / denominator;
  for (var j = 1; j < n; j++) {
    denominator = 1.0 + (*upper)[j] + (*lower)[j] + (*upper)[j] * gain[j - 1];
    gain[j] = -(*lower)[j] / denominator;
    (*rhs)[j] = ((*rhs)[j] + (*upper)[j] * (*rhs)[j - 1]) / denominator;
  }
  for (var j = n - 2; j >= 0; j--) { (*rhs)[j] -= gain[j] * (*rhs)[j + 1]; }
}
@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= N) { return; }
  var upper: array<f32, K>; var lower: array<f32, K>; var rhs: array<f32, K>;
  for (var j = 0; j < n; j++) { upper[j] = U[j * N + i]; lower[j] = L[j * N + i]; rhs[j] = R[j * N + i]; }
  thomas(n, &upper, &lower, &rhs);
  for (var j = 0; j < n; j++) { OUT[j * N + i] = rhs[j]; }
}`,
  inline: `
@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= N) { return; }
  var gain: array<f32, K>; var rhs: array<f32, K>;
  var denominator = 1.0 + U[i] + L[i];
  gain[0] = -L[i] / denominator;
  rhs[0] = R[i] / denominator;
  for (var j = 1; j < n; j++) {
    let u = U[j * N + i]; let l = L[j * N + i];
    denominator = 1.0 + u + l + u * gain[j - 1];
    gain[j] = -l / denominator;
    rhs[j] = (R[j * N + i] + u * rhs[j - 1]) / denominator;
  }
  for (var j = n - 2; j >= 0; j--) { rhs[j] -= gain[j] * rhs[j + 1]; }
  for (var j = 0; j < n; j++) { OUT[j * N + i] = rhs[j]; }
}`,
  fwd: `
fn thomas(n: i32, upper: ptr<function, array<f32, K>>, lower: ptr<function, array<f32, K>>, rhs: ptr<function, array<f32, K>>) {
  var gain: array<f32, K>;
  var denominator = 1.0 + (*upper)[0] + (*lower)[0];
  gain[0] = -(*lower)[0] / denominator;
  (*rhs)[0] = (*rhs)[0] / denominator;
  for (var j = 1; j < n; j++) {
    denominator = 1.0 + (*upper)[j] + (*lower)[j] + (*upper)[j] * gain[j - 1];
    gain[j] = -(*lower)[j] / denominator;
    (*rhs)[j] = ((*rhs)[j] + (*upper)[j] * (*rhs)[j - 1]) / denominator;
  }
  for (var m = 0; m < n - 1; m++) { let j = n - 2 - m; (*rhs)[j] -= gain[j] * (*rhs)[j + 1]; }
}
@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= N) { return; }
  var upper: array<f32, K>; var lower: array<f32, K>; var rhs: array<f32, K>;
  for (var j = 0; j < n; j++) { upper[j] = U[j * N + i]; lower[j] = L[j * N + i]; rhs[j] = R[j * N + i]; }
  thomas(n, &upper, &lower, &rhs);
  for (var j = 0; j < n; j++) { OUT[j * N + i] = rhs[j]; }
}`,
}[FORM];
const code = `const K: i32 = ${K}; const n: i32 = ${n}; const N: i32 = ${N};
@group(0) @binding(0) var<storage, read> U: array<f32>; @group(0) @binding(1) var<storage, read> L: array<f32>;
@group(0) @binding(2) var<storage, read> R: array<f32>; @group(0) @binding(3) var<storage, read_write> OUT: array<f32>;
${body}`;
let seed = 12345; const rand = () => { seed = (seed * 1664525 + 1013904223) >>> 0; return seed / 4294967296; };
const U = new Float32Array(n * N), L = new Float32Array(n * N), R = new Float32Array(n * N);
for (let i = 0; i < N; i++) for (let j = 0; j < n; j++) {
  const scale = Math.exp(6 * rand() - 3);
  U[j * N + i] = j > 0 ? f(scale * rand()) : 0; L[j * N + i] = j < n - 1 ? f(scale * rand()) : 0; R[j * N + i] = f(280 + 40 * rand());
}
const buf = (arr, usage) => { const b = device.createBuffer({ size: arr.byteLength, usage }); device.queue.writeBuffer(b, 0, arr); return b; };
const bU = buf(U, 0x88), bL = buf(L, 0x88), bR = buf(R, 0x88);
const bO = device.createBuffer({ size: R.byteLength, usage: 0x84 });
const read = device.createBuffer({ size: R.byteLength, usage: 0x9 });
device.pushErrorScope('validation');
const module = device.createShaderModule({ code });
const info = await module.getCompilationInfo(); for (const m of info.messages) console.log(m.type, m.lineNum, m.message);
const pipeline = device.createComputePipeline({ layout: 'auto', compute: { module, entryPoint: 'main' } });
const bind = device.createBindGroup({ layout: pipeline.getBindGroupLayout(0), entries: [bU, bL, bR, bO].map((b, k) => ({ binding: k, resource: { buffer: b } })) });
const enc = device.createCommandEncoder(); const pass = enc.beginComputePass(); pass.setPipeline(pipeline); pass.setBindGroup(0, bind); pass.dispatchWorkgroups(Math.ceil(N / 64)); pass.end();
enc.copyBufferToBuffer(bO, 0, read, 0, R.byteLength); device.queue.submit([enc.finish()]);
const err = await device.popErrorScope(); if (err) console.log('validation:', err.message);
await read.mapAsync(1); const out = new Float32Array(read.getMappedRange().slice(0));
// float32 JS twin
const ref = new Float32Array(n * N); let worst = 0, worstI = -1, worstJ = -1, bad = 0;
for (let i = 0; i < N; i++) {
  const gain = new Float32Array(n), rhs = new Float32Array(n);
  let den = f(1 + U[i] + L[i]); gain[0] = f(-L[i] / den); rhs[0] = f(R[i] / den);
  for (let j = 1; j < n; j++) { const u = U[j * N + i], l = L[j * N + i]; den = f(f(f(1 + u) + l) + f(u * gain[j - 1])); gain[j] = f(-l / den); rhs[j] = f(f(R[j * N + i] + f(u * rhs[j - 1])) / den); }
  for (let j = n - 2; j >= 0; j--) rhs[j] = f(rhs[j] - f(gain[j] * rhs[j + 1]));
  for (let j = 0; j < n; j++) { ref[j * N + i] = rhs[j]; const d = Math.abs(out[j * N + i] - rhs[j]); if (d > worst) { worst = d; worstI = i; worstJ = j; } if (d > 1e-3) bad++; }
}
console.log(`FORM=${FORM}: ${N} systems of ${n}: largest |gpu - js| ${worst.toExponential(2)} K at system ${worstI} layer ${worstJ}; values over 1e-3 apart: ${bad} of ${n * N}`);
if (worstI >= 0) console.log('  gpu', Array.from({ length: n }, (_, j) => out[j * N + worstI].toFixed(3)).join(' '), '\n  js ', Array.from({ length: n }, (_, j) => ref[j * N + worstI].toFixed(3)).join(' '));
process.exit(0);
