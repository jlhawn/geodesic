// node pagecmp.mjs <repo> <state.bin> [steps]: steps the page's GPU model from a saved state and prints hashes of its buffers.
// FLAGS sets Dawn toggles (e.g. enable-dawn-features=disable_robustness); OUT writes the S buffer as raw f32.
import { createHash } from 'node:crypto';
import { writeFileSync } from 'node:fs';
import { createRequire } from 'node:module';
const REPO = process.argv[2].replace(/\/$/, ''), FILE = process.argv[3], STEPS = Number(process.argv[4] ?? 32);
const require = createRequire('/Users/jlhawn/git_repos/jlhawn/geodesic/package.json');
const { create, globals } = await import(require.resolve('webgpu'));
Object.assign(globalThis, globals);
const flags = (process.env.FLAGS ?? '').split(' ').filter(Boolean);
Object.defineProperty(globalThis, 'navigator', { value: { gpu: create(flags) }, configurable: true });
const { gpuModelFrom, readState } = await import(`${REPO}/scripts/figures/figureState.mjs`);
const { readBuffer } = await import(`${REPO}/js/gpu/device.module.js`);
const saved = await readState(FILE);
const model = await gpuModelFrom(saved);
const dt = 1350 * 16 / saved.N;
for (let n = 0; n < STEPS; n++) { await model.step(dt); await model.settle(); }
const { device, buffers } = model.gpu;
const out = [];
for (const name of ['S', 'PH', 'D']) {
  const data = await readBuffer(device, buffers[name], buffers[name].size, Float32Array);
  out.push(`${name} ${createHash('sha256').update(Buffer.from(data.buffer)).digest('hex').slice(0, 12)}`);
  if (process.env.OUT) writeFileSync(`${process.env.OUT}.${name}`, Buffer.from(data.buffer));
}
if (process.env.OUT) writeFileSync(`${process.env.OUT}.layout.json`, JSON.stringify(model.gpu.layout));
const d = await model.diagnostics();
console.log(`${REPO.replace(/.*\//, '')} flags=[${flags.join(' ')}] N=${saved.N} after ${STEPS} steps: ${out.join(', ')}`);
process.exit(0);
