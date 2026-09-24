// One spin-up segment on the GPU engine with the page's configuration:
// continue from the newest runs/<TAG>_dayNNNN.bin (or start from the fresh
// initial state when there is none), step until MINUTES of wall time have
// passed or DAYS is reached, finishing the simulated day, save a binary
// snapshot and keep the two newest. Logs one line a day to runs/<TAG>.log
// and exits with 2 on NaN.
//
// Environment: N (128), TAG (spin<N>), MINUTES (15), DAYS (none), OUT
// (runs/), OCEAN (JSON options for the ocean, e.g. '{"closureHours":3}').
// scripts/spinup.sh runs segments back to back.
import { readFileSync, writeFileSync, readdirSync, renameSync, unlinkSync, appendFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { createGpuModel } from '../js/gpu/model.gpu.js';
import { decodeState, encodeState } from '../js/stateFile.module.js';

const N = Number(process.env.N ?? 128), TAG = process.env.TAG ?? `spin${N}`, MINUTES = Number(process.env.MINUTES ?? 15), DAYS = Number(process.env.DAYS ?? Infinity), KEEP = 2;
const OUT = process.env.OUT ?? new URL('../runs/', import.meta.url).pathname;
const OCEAN = JSON.parse(process.env.OCEAN ?? '{}');
const log = (line) => { console.log(line); appendFileSync(`${OUT}/${TAG}.log`, line + '\n'); };
const dayOf = (file) => Number(file.match(/_day(\d+)\.bin$/)[1]);
const snapshots = () => readdirSync(OUT).filter((f) => f.startsWith(`${TAG}_day`) && /_day\d+\.bin$/.test(f)).sort((a, b) => dayOf(a) - dayOf(b));

const t0 = performance.now();
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const model = await createGpuModel(new Grid(N), { topography, ocean: OCEAN });
const { mesh, core, state } = model;
const C = mesh.nCells, dt = 1350 * 16 / N, perDay = Math.round(86400 / dt);
const existing = snapshots();
if (existing.length) {
  const file = existing[existing.length - 1];
  const saved = await decodeState(new Uint8Array(readFileSync(`${OUT}/${file}`)));
  if (saved.N !== N) throw new Error(`${file} is N=${saved.N}`);
  ['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'].forEach((name, a) => state[a].set(saved[name]));
  model.time = saved.time;
  model.load();
  model.ocean.load(saved.ocean, state[3], state[6]);
  model.land.load({ soil: Float64Array.from(saved.land.soil), snow: Float64Array.from(saved.land.snow) });
  log(`--- ${new Date().toISOString()} continuing from ${file} (day ${saved.day}) after ${((performance.now() - t0) / 1000).toFixed(0)} s of setup`);
} else {
  initializeState(model, { geostrophic: !model.surfaceGeopotential }).forEach((values, a) => state[a].set(values));
  for (let i = 0; i < C; i++) if (model.geography.land[i]) state[6][i] = 0;
  model.load();
  model.ocean.initialize(state[3], state[6]);
  model.land.initialize();
  log(`--- ${new Date().toISOString()} fresh start at N=${N} (${C} cells, dt ${dt} s, ${perDay} steps a day) after ${((performance.now() - t0) / 1000).toFixed(0)} s of setup`);
}

await model.diagnostics();
const start = performance.now();
let day = Math.round(model.time / 86400);
for (;;) {
  for (let n = 0; n < perDay; n++) { await model.step(dt); if (n % 8 === 7) await model.settle(); }
  day++;
  const d = await model.diagnostics();
  const minutes = (performance.now() - start) / 60000;
  log(`day ${day} (${minutes.toFixed(1)} min): Ts ${(d.meanSurfaceT - 273.15).toFixed(2)} °C, ASR ${d.absorbedSolar.toFixed(1)} OLR ${d.outgoingLongwave.toFixed(1)} W/m², ps ${(d.piMin / 100).toFixed(0)}–${(d.piMax / 100).toFixed(0)} hPa, max wind ${d.maxWind.toFixed(1)} m/s, precip ${(86400 * d.precipitation).toFixed(2)} mm/d, ice ${(100 * d.iceFraction).toFixed(1)}%, ocean h1 ${d.oceanUpperDepth.toFixed(0)} m, interior ${(d.oceanInteriorT - 273.15).toFixed(2)} °C, currents ≤ ${d.oceanSpeed.toFixed(2)} m/s, transport ${d.oceanTransport.toFixed(0)} Sv, clamped ${d.oceanLimited}`);
  if (!Number.isFinite(d.meanSurfaceT) || !Number.isFinite(d.maxWind) || !Number.isFinite(d.oceanSpeed)) { log(`NaN on day ${day}; stopping`); process.exit(2); }
  if (minutes >= MINUTES || day >= DAYS) break;
}

await model.sync();
const ocean = await model.ocean.serialize(), land = await model.land.serialize();
const name = `${TAG}_day${String(day).padStart(4, '0')}.bin`;
const [pi, theta, u, surfaceT, q, qc, ice] = state;
writeFileSync(`${OUT}/${name}.partial`, encodeState({ N, K: core.K, day, time: model.time, terrain: !!model.surfaceGeopotential, pi, theta, u, surfaceT, q, qc, ice, ocean: { h: ocean.h, u: ocean.u, T: ocean.T, S: ocean.S, eta: ocean.eta }, land: { soil: land.soil, snow: land.snow } }));
renameSync(`${OUT}/${name}.partial`, `${OUT}/${name}`);
const kept = snapshots();
for (const old of kept.slice(0, Math.max(0, kept.length - KEEP))) unlinkSync(`${OUT}/${old}`);
log(`saved ${name} after ${((performance.now() - t0) / 60000).toFixed(1)} min; keeping ${snapshots().join(', ')}`);
process.exit(0);
