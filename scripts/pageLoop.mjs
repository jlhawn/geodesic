// The page's model rate in Node: the GPU loop of js/model.worker.js
// (loop() and yieldToPage(), at most QUEUE_DEPTH 2 submissions in flight,
// no pace pauses as no page reports late frames) from a saved state with
// the drivers' defaults (scripts/figures/figureState.mjs), each frame's
// fields (FIELDS, wind and speed at the surface) queued ahead of its
// stepsPerFrame steps (the worker's default for the resolution) and read
// back after them, for SECONDS (120) of wall time after a warm-up frame.
// MODE 'batch' queues each frame's steps as one batch (model.stepBatch),
// as the worker does; 'steps' queues and yields after every step, as it
// did before. Prints simulated hours per wall minute and seconds per
// model day.
//   MODE=steps node scripts/pageLoop.mjs runs/eleven128_day1825.bin
import { gpuModelFrom, readState } from './figures/figureState.mjs';

const file = process.argv[2] ?? new URL('../runs/eleven128_day1825.bin', import.meta.url).pathname;
const MODE = process.env.MODE ?? 'batch', SECONDS = Number(process.env.SECONDS ?? 120), FIELDS = (process.env.FIELDS ?? 'wind,speed').split(',');
const saved = await readState(file);
const model = await gpuModelFrom(saved);
const N = saved.N, dt = 1350 * 16 / N, stepsPerFrame = Math.max(2, Math.round(24 * 16 / N)), QUEUE_DEPTH = 2;
const subscription = { level: 'surface', depth: 'surface', fields: FIELDS, diagnostics: false }, inFlight = [];
async function yieldToPage() {
  inFlight.push(model.settle());
  if (inFlight.length >= QUEUE_DEPTH) await inFlight.shift();
}
async function frame() {
  const capturing = model.beginFrame(subscription);
  if (MODE === 'batch') { await model.stepBatch(stepsPerFrame, dt, null, false); await yieldToPage(); }
  else for (let n = 0; n < stepsPerFrame; n++) { await model.step(dt); await yieldToPage(); }
  await capturing;
  await new Promise((resolve) => setTimeout(resolve, 0));
}
await frame();
await Promise.all(inFlight.splice(0));
const start = performance.now(), time0 = model.time;
let frames = 0;
while (performance.now() - start < 1000 * SECONDS) { await frame(); frames++; }
await Promise.all(inFlight.splice(0));
const seconds = (performance.now() - start) / 1000, simulated = model.time - time0;
console.log(`${file.replace(/.*\//, '')} N=${N} ${MODE}: ${frames} frames of ${stepsPerFrame} steps in ${seconds.toFixed(1)} s, ${(simulated / 3600 / (seconds / 60)).toFixed(1)} simulated h/min, ${(seconds / (simulated / 86400)).toFixed(1)} s per model day`);
process.exit(0);
