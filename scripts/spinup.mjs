// One spin-up segment on the GPU engine with the page's configuration:
// continue from the newest <TAG>_dayDDDD.bin in OUT (or start from the
// fresh initial state when there is none), step until MINUTES of wall
// time have passed or DAYS is reached, finishing the simulated day, save a
// binary snapshot and keep the KEEP newest. Logs one line a day to
// <TAG>.log, with the sea-ice extent of each hemisphere (the area of the
// cells at least 15% covered), and at the end of the segment the rain,
// vegetation and surface temperature of the regions in BOXES, and exits
// with 2 on NaN.
//
// SIGTERM or SIGINT stops the segment after the ocean step in progress
// and exits 0: at a day's end it saves <TAG>_dayDDDD.bin as usual, inside
// a day <TAG>_dayDDDD_stepSSSS.bin (DDDD days and SSSS steps done, with
// the ocean's restart arrays and, when RECORD is set, the forcing
// recorder's part of the day), which the next segment continues from and
// deletes once it has saved a whole day.
// Whole-day names alone count towards KEEP, and only whole-day names match
// the day patterns of the shell drivers.
//
// Environment: N (128), TAG (spin<N>), MINUTES (15), DAYS (none), KEEP (2), OUT
// (runs/), OCEAN (JSON options for the ocean, e.g. '{"closureHours":3}'),
// RADIATION (JSON options for the radiation, e.g. '{"cloudSolarAbsorption":0}'),
// RECORD (a directory to write each day's ocean and sea-ice forcing into as
// forcing-DDDD.bin, see js/forcing.module.js), BATCH (1: the steps are queued
// one at a time, waiting every eighth; more: that many steps go to the
// GPU in one submission, byte-identical, for drivers where submitting
// costs more than it does on Metal), LAND_FROM (a state whose land seeds a
// fresh start), OCEAN_FROM (a state of the same N whose ocean, sea ice and
// sea surface replace the snapshot's, see js/oceanHandOff.module.js; the
// snapshots name it as oceanFrom, and a snapshot that already does is
// continued without replacing again), SYNC_CMD (a shell command run after
// every snapshot, forcing file and log update with the file's path as $1,
// see scripts/runControl.mjs), STOP_AFTER_STEPS (for tests: stop as on
// SIGTERM once this many steps have run).
// scripts/spinup.sh runs segments back to back, scripts/asyncSpinup.sh
// alternates them with ocean-only spin-ups.
import { readFileSync, writeFileSync, readdirSync, renameSync, unlinkSync, appendFileSync, mkdirSync } from 'node:fs';
import { basename } from 'node:path';
import { Grid } from '../js/grid.module.js';
import { topographyFromInt16, createGeography } from '../js/geography.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { createGpuModel } from '../js/gpu/model.gpu.js';
import { decodeState, encodeState } from '../js/stateFile.module.js';
import { savedSubsidence, regridLand } from '../js/physics/regrid.module.js';
import { readRanges } from '../js/gpu/device.module.js';
import { LAYER_DENSITIES, THERMOCLINE_DENSITY } from '../js/ocean/layered.module.js';
import { createForcingRecorder } from '../js/gpu/forcing.gpu.js';
import { forcingName } from '../js/forcing.module.js';
import { withOceanOf } from '../js/oceanHandOff.module.js';
import { stopOnSignal, syncAfterSave } from './runControl.mjs';

const BOXES = {
  sahara: [16, 30, -10, 32], arabia: [16, 30, 38, 55], sahel: [8, 16, -15, 35], india: [15, 28, 72, 88], congo: [-5, 5, 12, 30], amazon: [-10, 3, -70, -50],
  seAsia: [10, 25, 95, 110], borneo: [-4, 7, 108, 119], europe: [45, 55, 0, 30], eastUS: [32, 45, -95, -75], siberia: [55, 65, 60, 120],
  ausInterior: [-30, -20, 120, 145], kalahari: [-27, -20, 17, 25], gobi: [38, 46, 90, 110], usSouthwest: [30, 37, -117, -106], cerrado: [-20, -10, -55, -42],
  ausNorth: [-18, -11, 125, 145], ausEast: [-37, -25, 148, 154], ausSoutheast: [-43, -34, 140, 150], ausWest: [-30, -20, 114, 120], newGuinea: [-11, -1, 130, 151],
};

const N = Number(process.env.N ?? 128), TAG = process.env.TAG ?? `spin${N}`, MINUTES = Number(process.env.MINUTES ?? 15), DAYS = Number(process.env.DAYS ?? Infinity), KEEP = Number(process.env.KEEP ?? 2);
const OUT = process.env.OUT ?? new URL('../runs/', import.meta.url).pathname;
const OCEAN = JSON.parse(process.env.OCEAN ?? '{}');
const RADIATION = JSON.parse(process.env.RADIATION ?? '{}');
const OCEAN_FROM = process.env.OCEAN_FROM, STOP_AFTER_STEPS = Number(process.env.STOP_AFTER_STEPS ?? Infinity);
const log = (line) => { console.log(line); appendFileSync(`${OUT}/${TAG}.log`, line + '\n'); };
const stop = stopOnSignal(log), hook = syncAfterSave(process.env.SYNC_CMD, log);
const positionOf = (file) => { const [, day, step] = file.match(/_day(\d+)(?:_step(\d+))?\.bin$/); return Number(day) * 1e6 + Number(step ?? 0); };
const inDay = (file) => /_step\d+\.bin$/.test(file);
const snapshots = () => readdirSync(OUT).filter((f) => f.startsWith(`${TAG}_day`) && /_day\d+(?:_step\d+)?\.bin$/.test(f)).sort((a, b) => positionOf(a) - positionOf(b));

const t0 = performance.now();
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const model = await createGpuModel(new Grid(N), { topography, ocean: OCEAN, radiation: RADIATION });
const { mesh, core, state } = model;
const C = mesh.nCells, dt = 1350 * 16 / N, perDay = Math.round(86400 / dt), BATCH = Math.max(1, Math.round(Number(process.env.BATCH ?? 1)));
function loadSaved(saved) {
  ['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'].forEach((name, a) => state[a].set(saved[name]));
  model.seaIce.load(state[6], saved.concentration ?? null);
  model.radiation.mlmSubsidence.set(savedSubsidence(saved, model));
  model.time = saved.time;
  model.load();
  model.ocean.load(saved.ocean, state[3], state[6]);
  model.land.load({ soil: Float64Array.from(saved.land.soil), snow: Float64Array.from(saved.land.snow), ...(saved.land.vegetation ? { vegetation: Float64Array.from(saved.land.vegetation) } : {}), ...(saved.land.surface ? { surface: Float64Array.from(saved.land.surface) } : {}) });
}
const existing = snapshots();
let saved = null, startStep = 0, oceanYears = 0, oceanFrom = null;
if (existing.length) {
  const file = existing[existing.length - 1];
  saved = await decodeState(new Uint8Array(readFileSync(`${OUT}/${file}`)));
  if (saved.N !== N) throw new Error(`${file} is N=${saved.N}`);
  startStep = saved.step ?? 0;
  oceanYears = saved.oceanYears ?? 0;
  oceanFrom = saved.oceanFrom ?? null;
  if (OCEAN_FROM && oceanFrom === basename(OCEAN_FROM)) log(`${file} already carries the ocean of ${OCEAN_FROM}`);
  else if (OCEAN_FROM) {
    const alone = await decodeState(new Uint8Array(readFileSync(OCEAN_FROM)));
    saved = withOceanOf(saved, alone, model.geography.land);
    oceanYears = saved.oceanYears;
    oceanFrom = basename(OCEAN_FROM);
    log(`ocean, sea ice and sea surface of ${file} replaced by those of ${OCEAN_FROM} (day ${alone.day}, ${alone.oceanYears ?? 0} years alone; ${oceanYears} in all)`);
  }
  loadSaved(saved);
  log(`--- ${new Date().toISOString()} continuing from ${file} (day ${saved.day}${startStep ? `, step ${startStep} of ${perDay}` : ''}) after ${((performance.now() - t0) / 1000).toFixed(0)} s of setup`);
} else {
  if (OCEAN_FROM) throw new Error(`OCEAN_FROM needs a ${TAG}_dayDDDD.bin in ${OUT} to hand the ocean to`);
  initializeState(model, { geostrophic: !model.surfaceGeopotential }).forEach((values, a) => state[a].set(values));
  for (let i = 0; i < C; i++) if (model.geography.land[i]) state[6][i] = 0;
  model.load();
  model.ocean.initialize(state[3], state[6]);
  model.land.initialize();
  if (process.env.LAND_FROM) {
    const seed = await decodeState(new Uint8Array(readFileSync(process.env.LAND_FROM)));
    const sourceMesh = seed.N === N ? mesh : buildMesh(new Grid(seed.N));
    const source = { mesh: sourceMesh, geography: seed.N === N ? model.geography : createGeography(sourceMesh, topography) };
    model.land.load(regridLand(source, { mesh, geography: model.geography, land: model.land }, seed.land, null, { ice: seed.ice, surfaceT: seed.surfaceT }), state[6]);
    log(`land seeded from ${process.env.LAND_FROM} (N=${seed.N}, day ${seed.day})`);
  }
  log(`--- ${new Date().toISOString()} fresh start at N=${N} (${C} cells, dt ${dt} s, ${perDay} steps a day) after ${((performance.now() - t0) / 1000).toFixed(0)} s of setup`);
}
const everySteps = model.oceanEngine ? model.oceanEngine.everySteps : 1;
if (startStep % everySteps) throw new Error(`the snapshot stopped at step ${startStep}, off the ocean's ${everySteps}-step cadence`);

const deg = 180 / Math.PI, land = model.geography.land;
const inBox = Object.fromEntries(Object.entries(BOXES).map(([k, [a, b, c, e]]) => [k, [...Array(C).keys()].filter((i) => land[i] && mesh.latCell[i] * deg >= a && mesh.latCell[i] * deg <= b && mesh.lonCell[i] * deg >= c && mesh.lonCell[i] * deg <= e)]));
const PH = model.gpu.layout.PH;
const readRain = async () => { const [a, b] = await readRanges(model.gpu.device, model.gpu.buffers.PH, [{ offset: PH.CONV, length: C }, { offset: PH.COND, length: C }]); return Float64Array.from(a, (x, i) => x + b[i]); };
const iceArea = async () => { const { fields } = await model.beginFrame({ fields: ['concentration'] }); let north = 0, south = 0; for (let i = 0; i < C; i++) if (fields.concentration[i] >= 0.15) { if (mesh.latCell[i] > 0) north += mesh.areaCell[i]; else south += mesh.areaCell[i]; } return [north / 1e12, south / 1e12]; };
const rain0 = await readRain();
await model.diagnostics();
const RECORD = process.env.RECORD;
if (RECORD) mkdirSync(RECORD, { recursive: true });
const partOfDay = RECORD && startStep && saved.forcingSums ? { sums: saved.forcingSums, steps: saved.forcingSteps, oceanSteps: saved.forcingOceanSteps, seconds: saved.forcingSeconds, rain: saved.forcingRain, runoff: saved.forcingRunoff } : null;
if (RECORD && startStep && !partOfDay) log(`the snapshot holds no recorded part of its day, so forcing day ${saved.day + 1} covers its last ${perDay - startStep} steps alone`);
const recorder = RECORD ? await createForcingRecorder(model, partOfDay) : null;
if (stop.requested) { log(`stopped by ${stop.requested} before the first step; nothing to save`); await hook.drain(); process.exit(0); }
const start = performance.now();
const day0 = Math.round((model.time - startStep * dt) / 86400);
let day = day0, step = startStep, taken = 0, iceNorth = 0, iceSouth = 0;
const halted = () => stop.requested || taken >= STOP_AFTER_STEPS;
for (;;) {
  if (BATCH === 1) {
    while (step < perDay) {
      await model.step(dt); recorder?.step(); step++; taken++;
      if (step % 8 === 0) await model.settle();
      if (halted() && step % everySteps === 0) break;
    }
  } else {
    while (step < perDay) {
      const count = Math.min(BATCH, perDay - step);
      await model.stepBatch(count, dt, recorder ? () => recorder.step() : null);
      step += count; taken += count;
      if (halted() && step % everySteps === 0) break;
    }
  }
  if (step < perDay) break;
  step = 0;
  day++;
  const d = await model.diagnostics();
  if (recorder) {
    const file = `${RECORD}/${forcingName(day)}`;
    writeFileSync(`${file}.partial`, await recorder.day(day));
    renameSync(`${file}.partial`, file);
    hook.after(file);
  }
  const [north, south] = await iceArea();
  iceNorth += north; iceSouth += south;
  const minutes = (performance.now() - start) / 60000;
  log(`day ${day} (${minutes.toFixed(1)} min): Ts ${(d.meanSurfaceT - 273.15).toFixed(2)} °C, ASR ${d.absorbedSolar.toFixed(1)} (atmosphere ${d.atmosphereSolar.toFixed(1)}) OLR ${d.outgoingLongwave.toFixed(1)} W/m², ps ${(d.piMin / 100).toFixed(0)}–${(d.piMax / 100).toFixed(0)} hPa, max wind ${d.maxWind.toFixed(1)} m/s, precip ${(86400 * d.precipitation).toFixed(2)} mm/d, ice ${(100 * d.iceFraction).toFixed(1)}% (N ${north.toFixed(1)} S ${south.toFixed(1)} Mkm²), albedo ${d.planetaryAlbedo.toFixed(3)}, ocean h1 ${d.oceanUpperDepth.toFixed(0)} m, interior ${(d.oceanInteriorT - 273.15).toFixed(2)} °C, currents ≤ ${d.oceanSpeed.toFixed(2)} m/s, transport ${d.oceanTransport.toFixed(0)} Sv, clamped ${d.oceanLimited}`);
  if (!Number.isFinite(d.meanSurfaceT) || !Number.isFinite(d.maxWind) || !Number.isFinite(d.oceanSpeed)) { log(`NaN on day ${day}; stopping`); await hook.drain(); process.exit(2); }
  if (minutes >= MINUTES || day >= DAYS || halted()) break;
}

const partial = {};
if (step && recorder) {
  await model.diagnostics();
  const part = await recorder.checkpoint();
  Object.assign(partial, { forcingSteps: part.steps, forcingOceanSteps: part.oceanSteps, forcingSeconds: part.seconds, forcingSums: part.sums, forcingRain: part.rain, forcingRunoff: part.runoff });
}
await model.sync();
const rain1 = await readRain(), days = day - day0;
const ocean = await model.ocean.serialize({ restart: step > 0 }), landState = await model.land.serialize();
const kT = LAYER_DENSITIES.findIndex((r) => r >= THERMOCLINE_DENSITY);
const sea = (a, b, c, e) => [...Array(C).keys()].filter((i) => !land[i] && mesh.latCell[i] * deg >= a && mesh.latCell[i] * deg <= b && (c <= e ? mesh.lonCell[i] * deg >= c && mesh.lonCell[i] * deg <= e : mesh.lonCell[i] * deg >= c || mesh.lonCell[i] * deg <= e));
const meanOf = (cells, f) => cells.reduce((s, i) => s + f(i), 0) / Math.max(1, cells.length);
const sst = (cells) => meanOf(cells, (i) => ocean.T[i] - 273.15), classTop = (cells) => meanOf(cells, (i) => { let d = 0; for (let k = 0; k < kT; k++) d += ocean.h[k * C + i]; return d; });
const warmPool = sea(-10, 10, 120, 160), coldTongue = sea(-2, 2, -110, -90), westPacific = sea(-5, 5, 140, 170), eastPacific = sea(-5, 5, -120, -90);
if (days > 0 && !step) {
  log(`regions after ${days} days (rain mm/d / vegetation / surface °C): ` + Object.entries(inBox).map(([name, cells]) => {
    let r = 0, v = 0, t = 0;
    for (const i of cells) { r += rain1[i] - rain0[i]; v += landState.vegetation ? landState.vegetation[i] : 0; t += state[3][i]; }
    const n = Math.max(1, cells.length);
    return `${name} ${(r / n / days).toFixed(1)}/${(v / n).toFixed(2)}/${(t / n - 273.15).toFixed(0)}`;
  }).join(', ') + `; sea ice mean N ${(iceNorth / days).toFixed(1)} S ${(iceSouth / days).toFixed(1)} Mkm²`);
  log(`ocean after ${days} days: warm pool ${sst(warmPool).toFixed(1)} °C, cold tongue ${sst(coldTongue).toFixed(1)} °C (W−E ${(sst(westPacific) - sst(eastPacific)).toFixed(1)} K), ${THERMOCLINE_DENSITY} class top W Pac ${classTop(westPacific).toFixed(0)} m, E Pac ${classTop(eastPacific).toFixed(0)} m`);
}
const name = `${TAG}_day${String(day).padStart(4, '0')}${step ? `_step${String(step).padStart(4, '0')}` : ''}.bin`;
const [pi, theta, u, surfaceT, q, qc, ice] = state, { concentration } = model.seaIce, { mlmSubsidence } = model.radiation;
const header = { N, K: core.K, day, ...(step ? { step } : {}), time: model.time, terrain: !!model.surfaceGeopotential, ...(oceanYears ? { oceanYears } : {}), ...(oceanFrom ? { oceanFrom } : {}) };
writeFileSync(`${OUT}/${name}.partial`, encodeState({ ...header, pi, theta, u, surfaceT, q, qc, ice, concentration, mlmSubsidence, ocean: step ? ocean : { h: ocean.h, u: ocean.u, T: ocean.T, S: ocean.S, eta: ocean.eta }, land: landState, ...partial }, { f64: ['forcingRain', 'forcingRunoff'] }));
renameSync(`${OUT}/${name}.partial`, `${OUT}/${name}`);
const kept = snapshots(), whole = kept.filter((f) => !inDay(f));
for (const old of whole.slice(0, Math.max(0, whole.length - KEEP))) unlinkSync(`${OUT}/${old}`);
for (const old of kept) if (inDay(old) && old !== name) unlinkSync(`${OUT}/${old}`);
log(`saved ${name} after ${((performance.now() - t0) / 60000).toFixed(1)} min; keeping ${snapshots().join(', ')}`);
if (halted()) log(`stopped by ${stop.requested ?? `STOP_AFTER_STEPS=${STOP_AFTER_STEPS}`} at day ${day}${step ? ` and ${step} of ${perDay} steps` : ''}`);
hook.after(`${OUT}/${name}`);
hook.after(`${OUT}/${TAG}.log`);
await hook.drain();
process.exit(0);
