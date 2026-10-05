// The GPU model's pace in wall seconds per model day from a saved state,
// stepped as the page and the spin-up step it: loads the state as
// scripts/profileGpu.mjs does, steps WARM steps (8), then DAYS (1) model
// days, either each step queued and awaited before the next (BATCH 1, the
// page's worker) or BATCH steps recorded into one submission and awaited
// together (model.stepBatch, the spin-up's BATCH), with one diagnostics
// frame a day. Prints the seconds per model day and the milliseconds per
// step.
//   BATCH=8 node scripts/paceGpu.mjs runs/eleven128_day1825.bin
import { gpuModelFrom, readState } from './figures/figureState.mjs';

const file = process.argv[2] ?? new URL('../runs/eleven128_day1825.bin', import.meta.url).pathname;
const WARM = Number(process.env.WARM ?? 8), DAYS = Number(process.env.DAYS ?? 1), BATCH = Math.max(1, Math.round(Number(process.env.BATCH ?? 1)));
const saved = await readState(file);
const model = await gpuModelFrom(saved);
const dt = 1350 * 16 / saved.N, perDay = Math.round(86400 / dt), total = Math.round(DAYS * perDay);
for (let n = 0; n < WARM; n++) await model.step(dt);
await model.settle();
const start = performance.now();
let d = null;
for (let done = 0; done < total;) {
  const count = Math.min(BATCH, total - done, perDay - (done % perDay));
  if (BATCH === 1) { await model.step(dt); await model.settle(); } else await model.stepBatch(count, dt);
  done += count;
  if (done % perDay === 0) d = await model.diagnostics();
}
const seconds = (performance.now() - start) / 1000;
d ??= await model.diagnostics();
console.log(`${file.replace(/.*\//, '')} N=${saved.N}, ${total} steps (${DAYS} d) ${BATCH === 1 ? 'awaited one at a time' : `queued ${BATCH} to a submission`}: ${(seconds / DAYS).toFixed(1)} s per model day, ${(1000 * seconds / total).toFixed(1)} ms per step; ${total % perDay ? 'mean' : 'day mean'} ASR ${d.absorbedSolar.toFixed(2)} OLR ${d.outgoingLongwave.toFixed(2)} W/m², Ts ${(d.meanSurfaceT - 273.15).toFixed(3)} °C`);
process.exit(0);
