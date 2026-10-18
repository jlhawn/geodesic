// The GPU ocean's time per ocean step, alone and by phase, from a saved
// state: loads the state as scripts/profileGpu.mjs does (OCEAN as JSON in
// the environment picks the ocean, {"vertical":"hybrid"} for the hybrid on
// a slab state), steps the coupled model for two ocean intervals so the
// stress and surface are the atmosphere's, then times STEPS (32) ocean
// steps (advanceCoupled with the ocean's own step) each awaited alone, and
// each phase of the hybrid's step the same way, REPS (3) times over. SPLIT=1
// adds js/gpu/profile.module.js's per-kernel split of the ocean step, ROWS
// (12) rows of it. One GPU model at a time (CLAUDE.md); timings under
// contention are worthless, so alternate base and candidate runs.
//   OCEAN='{"vertical":"hybrid"}' node scripts/paceOcean.mjs runs/eleven128h_day1825.bin
import { gpuModelFrom, readState, dt as dtOf } from './figures/figureState.mjs';
import { profileGpu } from '../js/gpu/profile.module.js';

const file = process.argv[2], STEPS = Number(process.argv[3] ?? process.env.STEPS ?? 32), REPS = Number(process.env.REPS ?? 3), ROWS = Number(process.env.ROWS ?? 12);
if (!file) { console.error('usage: [OCEAN=json] [REPS=3] [SPLIT=1] node scripts/paceOcean.mjs <state.bin> [STEPS=32]'); process.exit(1); }
const saved = await readState(file);
const model = await gpuModelFrom(saved);
const ocean = model.oceanEngine, { device } = model.gpu, dt = dtOf(saved.N), every = ocean.everySteps ?? 1, dtO = dt * every;
const settle = () => device.queue.onSubmittedWorkDone();
for (let n = 0; n < 2 * every; n++) await model.step(dt);
await settle();
const median = (a) => [...a].sort((x, y) => x - y)[a.length >> 1];
const time = async (label, run) => {
  const times = [];
  for (let n = 0; n < STEPS; n++) { const start = performance.now(); run(); await settle(); times.push(performance.now() - start); }
  console.log(`  ${label.padEnd(10)} ${median(times).toFixed(2).padStart(8)} ms (min ${Math.min(...times).toFixed(2)})`);
  return median(times);
};
console.log(`${file.replace(/.*\//, '')} N=${saved.N}, ${ocean.vertical} ocean, step ${dtO} s, ${STEPS} steps awaited one at a time`);
for (let r = 0; r < REPS; r++) {
  const full = await time('advance', () => ocean.advanceCoupled(dtO));
  const names = ['dynamics', 'rescale', 'eddy', 'shift', 'regrid', 'column'].filter((p) => ocean.phases[p]);
  let sum = 0;
  for (const p of names) sum += await time(p, () => (p === 'rescale' || p === 'shift' ? ocean.phases[p]() : ocean.phases[p](dtO)));
  console.log(`  phases ${sum.toFixed(1)} ms of ${full.toFixed(1)} per ocean step; ${(full * 86400 / dtO / 1000).toFixed(2)} s a simulated day for the ocean alone`);
}
if (process.env.SPLIT === '1') {
  const p = await profileGpu({ gpu: model.gpu, step: (d) => ocean.advanceCoupled(d) }, { steps: STEPS, dt: dtO, split: true });
  console.log(`split: step median ${p.stepMedian.toFixed(1)} ms, GPU ${p.gpuMs?.toFixed(1)} ms`);
  for (const k of (p.kernels ?? []).slice(0, ROWS)) console.log(`  ${k.ms.toFixed(3).padStart(8)} ms ${k.passes.toFixed(1).padStart(5)} dispatches ${(1000 * k.median).toFixed(0).padStart(7)} µs median  ${k.name}`);
}
process.exit(0);
