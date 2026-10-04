// The GPU model's time per step and where it goes, kernel by kernel, from a
// saved state: loads the state as a continuing spin-up would (the options
// as JSON in OCEAN, RADIATION, MOIST, BOUNDARY_LAYER, SURFACE and LAND, see
// scripts/figures/figureState.mjs), steps WARM steps (8) to settle, then
// profiles STEPS steps (32: two ocean calls and eight full radiation
// calls at N=128 under js/cadence.module.js) with
// js/gpu/profile.module.js: each step submitted alone and awaited, its
// wall time from submission to completion, and with the device's
// timestamp queries the GPU time of every compute pass, charged to the
// kernels it dispatched and averaged per step (so the ocean's passes are
// amortised over its call interval, the radiation over its). Prints the
// step median, the empty round trip, the summed GPU time and the ROWS (28)
// costliest passes. SPLIT=1 times every dispatch in a pass of its own and
// prints each kernel's time with its dispatches a step, the mean, median
// and least of its dispatches' times, its workgroups a dispatch and the
// kernels of the engine's pass it belongs to.
//   node scripts/profileGpu.mjs [state.bin]
import { gpuModelFrom, readState } from './figures/figureState.mjs';
import { profileGpu } from '../js/gpu/profile.module.js';

const file = process.argv[2] ?? new URL('../runs/eleven128_day1825.bin', import.meta.url).pathname;
const WARM = Number(process.env.WARM ?? 8), STEPS = Number(process.env.STEPS ?? 32), ROWS = Number(process.env.ROWS ?? 28), SPLIT = process.env.SPLIT === '1';
const saved = await readState(file);
const model = await gpuModelFrom(saved);
const dt = 1350 * 16 / saved.N;
for (let n = 0; n < WARM; n++) await model.step(dt);
await model.settle();
const p = await profileGpu(model, { steps: STEPS, dt, split: SPLIT });
console.log(`${file.replace(/.*\//, '')} N=${saved.N}, ${STEPS} steps${SPLIT ? ', each dispatch its own pass' : ''}; ${p.device}: step median ${p.stepMedian.toFixed(1)} ms (min ${p.stepMin.toFixed(1)}, max ${p.stepMax.toFixed(1)}), round trip ${p.roundTrip.toFixed(2)} ms, GPU time ${p.gpuMs?.toFixed(1)} ms/step over ${p.timestamps ? 'timestamps' : 'no timestamps'}`);
for (const k of (p.kernels ?? []).slice(0, ROWS)) {
  const unresolved = k.unresolved ? ` (unresolved ${k.unresolved.toFixed(1)})` : '';
  if (SPLIT) console.log(`  ${k.ms.toFixed(3).padStart(8)} ms  ${k.passes.toFixed(3).padStart(7)} dispatches  ${(1000 * k.ms / k.passes).toFixed(1).padStart(8)} µs each (median ${(1000 * k.median).toFixed(1)}, least ${(1000 * k.least).toFixed(1)})  ${String(k.groups).padStart(7)} workgroups  ${k.name}  in [${k.within}]${unresolved}`);
  else console.log(`  ${k.ms.toFixed(2).padStart(7)} ms  ${k.passes.toFixed(1).padStart(5)} passes  ${k.name}${unresolved}`);
}
process.exit(0);
