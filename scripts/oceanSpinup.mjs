// The ocean and its sea ice spun up alone under a recorded year of surface
// forcing (spinup.mjs RECORD, js/forcing.module.js), looped for YEARS
// years, so that a coupled run can start from a spun-up ocean. The GPU
// model is built as spinup.mjs builds it and the state loaded into it, but
// its atmosphere never steps: each day of FORCING in turn, in day order,
// drives the sea surface, the ice and the ocean at the coupled model's dt
// and ocean cadence (createForcedOcean in js/gpu/forcing.gpu.js), with the
// open water restored to the recorded SST at RESTORE W/m²/K. At each
// year's end it saves OUT/<TAG>_yearYYYY.bin, a full state whose
// atmosphere and land are those loaded (the land keeps the snow the ice
// carries), its time advanced by the days run, keeps the KEEP newest,
// and logs to OUT/<TAG>.log the ocean line of spinup.mjs, the Southern
// Ocean column at 60–70S, the SST and sea-ice extent against the
// recorded last day, the global mean SST and its drift since the start
// and over the last ten years (or all of them when fewer), the ocean's
// interior temperature, fastest current, largest transport and clamped
// count, and the wall time; NaN in the SST or the ocean exits with 2.
// Every SNAPSHOT_DAYS days of the year between those it saves
// <TAG>_yearYYYY_dayDDD.bin (YYYY years and DDD days done), and SIGTERM
// or SIGINT stops it after the ocean step in progress, saving
// <TAG>_yearYYYY_dayDDD_stepSSSS.bin when that falls inside a day, and
// exits 0. Only the newest of those in-year files is kept. It continues
// from the newest of all these files in OUT when there is one, and
// otherwise starts from STATE, which should be a snapshot at the start or
// the end of the recorded cycle so that its season matches. Each file
// carries the global mean SST at STATE and at each year's end since
// (sstByYear), its position (oceanYears whole years, oceanDays days in
// all and oceanStep steps into the next) and the ocean's restart arrays
// (restartArrays in js/gpu/layeredOcean.gpu.js), so that a run continued
// from any of them follows the uninterrupted one bit for bit.
//
// Environment: N (64), STATE, FORCING, YEARS (1, the year to stop after),
// DAYS_PER_YEAR (365), LOOP_DAYS (DAYS_PER_YEAR: the first that many days
// of FORCING are looped, a whole number of years), RESTORE (30), TAG
// (ocean<N>), OUT (runs/), KEEP (2), SNAPSHOT_DAYS (30; 0 saves at the
// years' ends alone), OCEAN and RADIATION (JSON options as in spinup.mjs;
// OCEAN's everySteps sets the ocean step), SYNC_CMD and STOP_AFTER_STEPS
// (as in spinup.mjs).
import { readFileSync, writeFileSync, readdirSync, renameSync, unlinkSync, appendFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { createGpuModel } from '../js/gpu/model.gpu.js';
import { decodeState, encodeState } from '../js/stateFile.module.js';
import { savedSubsidence } from '../js/physics/regrid.module.js';
import { readRanges } from '../js/gpu/device.module.js';
import { LAYER_DENSITIES, THERMOCLINE_DENSITY } from '../js/ocean/layered.module.js';
import { FREEZING_POINT } from '../js/physics/ice.module.js';
import { createForcedOcean } from '../js/gpu/forcing.gpu.js';
import { decodeForcing, forcingDay } from '../js/forcing.module.js';
import { stopOnSignal, syncAfterSave } from './runControl.mjs';

const N = Number(process.env.N ?? 64), TAG = process.env.TAG ?? `ocean${N}`, YEARS = Number(process.env.YEARS ?? 1), KEEP = Number(process.env.KEEP ?? 2);
const DAYS_PER_YEAR = Number(process.env.DAYS_PER_YEAR ?? 365), LOOP_DAYS = Number(process.env.LOOP_DAYS ?? DAYS_PER_YEAR), RESTORE = Number(process.env.RESTORE ?? 30);
const SNAPSHOT_DAYS = Number(process.env.SNAPSHOT_DAYS ?? 30), STOP_AFTER_STEPS = Number(process.env.STOP_AFTER_STEPS ?? Infinity);
const OUT = process.env.OUT ?? new URL('../runs/', import.meta.url).pathname, FORCING = process.env.FORCING;
const OCEAN = JSON.parse(process.env.OCEAN ?? '{}'), RADIATION = JSON.parse(process.env.RADIATION ?? '{}');
if (!FORCING) throw new Error('FORCING must name a directory of recorded forcing');
if (LOOP_DAYS % DAYS_PER_YEAR) throw new Error(`LOOP_DAYS ${LOOP_DAYS} is not a whole number of ${DAYS_PER_YEAR}-day years`);
const log = (line) => { console.log(line); appendFileSync(`${OUT}/${TAG}.log`, line + '\n'); };
const stop = stopOnSignal(log), hook = syncAfterSave(process.env.SYNC_CMD, log);
const pad = (x, width) => String(x).padStart(width, '0');
const FILE = /_year(\d+)(?:_day(\d+)(?:_step(\d+))?)?\.bin$/;
const positionOf = (file) => { const [, year, day, step] = file.match(FILE); return (Number(year) * 1e4 + Number(day ?? 0)) * 1e5 + Number(step ?? 0); };
const inYear = (file) => /_day\d+(?:_step\d+)?\.bin$/.test(file);
const snapshots = () => readdirSync(OUT).filter((f) => f.startsWith(`${TAG}_year`) && FILE.test(f)).sort((a, b) => positionOf(a) - positionOf(b));
const nameAt = (done, step) => {
  const year = Math.floor(done / DAYS_PER_YEAR), day = done % DAYS_PER_YEAR;
  return `${TAG}_year${pad(year, 4)}${day || step ? `_day${pad(day, 3)}` : ''}${step ? `_step${pad(step, 4)}` : ''}.bin`;
};

const t0 = performance.now();
const days = readdirSync(FORCING).filter((f) => forcingDay(f) !== null).sort((a, b) => forcingDay(a) - forcingDay(b)).slice(0, LOOP_DAYS);
if (days.length < LOOP_DAYS) throw new Error(`${FORCING} holds ${days.length} days of forcing, not ${LOOP_DAYS}`);
const existing = snapshots();
const from = existing.length ? `${OUT}/${existing[existing.length - 1]}` : process.env.STATE;
if (!from) throw new Error(`no ${TAG}_yearYYYY.bin in ${OUT} and no STATE to start from`);
const saved = await decodeState(new Uint8Array(readFileSync(from)));
if (saved.N !== N) throw new Error(`${from} is N=${saved.N}`);

const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const model = await createGpuModel(new Grid(N), { topography, ocean: OCEAN, radiation: RADIATION });
const { mesh, core, state, gpu } = model;
const C = mesh.nCells, dt = 1350 * 16 / N, perDay = Math.round(86400 / dt), everySteps = model.oceanEngine.everySteps, oceanDt = everySteps * dt;
const doneBefore = existing.length ? saved.oceanDays ?? (saved.oceanYears ?? 0) * DAYS_PER_YEAR : 0;
const stepBefore = existing.length ? saved.oceanStep ?? 0 : 0;
if (stepBefore % everySteps) throw new Error(`${from} stopped at step ${stepBefore}, off the ocean's ${everySteps}-step cadence`);
const baseDay = saved.day - doneBefore, baseTime = saved.time - doneBefore * 86400 - stepBefore * dt;
['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'].forEach((name, a) => state[a].set(saved[name]));
model.seaIce.load(state[6], saved.concentration ?? null);
model.radiation.mlmSubsidence.set(savedSubsidence(saved, model));
model.time = saved.time;
model.load();
model.ocean.load(saved.ocean, state[3], state[6]);
model.land.load({ soil: Float64Array.from(saved.land.soil), snow: Float64Array.from(saved.land.snow), ...(saved.land.vegetation ? { vegetation: Float64Array.from(saved.land.vegetation) } : {}), ...(saved.land.surface ? { surface: Float64Array.from(saved.land.surface) } : {}) });
const forced = createForcedOcean(model);

const deg = 180 / Math.PI, land = model.geography.land;
const kT = LAYER_DENSITIES.findIndex((r) => r >= THERMOCLINE_DENSITY);
const sea = (a, b, c, e) => [...Array(C).keys()].filter((i) => !land[i] && mesh.latCell[i] * deg >= a && mesh.latCell[i] * deg <= b && (c <= e ? mesh.lonCell[i] * deg >= c && mesh.lonCell[i] * deg <= e : mesh.lonCell[i] * deg >= c || mesh.lonCell[i] * deg <= e));
const meanOf = (cells, f) => cells.reduce((s, i) => s + f(i), 0) / Math.max(1, cells.length);
const warmPool = sea(-10, 10, 120, 160), coldTongue = sea(-2, 2, -110, -90), westPacific = sea(-5, 5, 140, 170), eastPacific = sea(-5, 5, -120, -90);
const southern = sea(-70, -60, -180, 180), ocean = sea(-90, 90, -180, 180);
const BANDS = [[0, 60], [60, 200], [200, 500], [500, 1000]];
const S = gpu.layout.S, PH = gpu.layout.PH;
const globalSST = (surfaceT, ice) => {
  let area = 0, sum = 0;
  for (const i of ocean) { area += mesh.areaCell[i]; sum += mesh.areaCell[i] * (ice[i] > 0 ? FREEZING_POINT : surfaceT[i]); }
  return sum / area;
};
const history = existing.length && saved.sstByYear ? Array.from(saved.sstByYear) : [];
if (!history.length) history[Math.floor(doneBefore / DAYS_PER_YEAR)] = globalSST(saved.surfaceT, saved.ice);
const fixed = (x, digits) => (Number.isFinite(x) ? x.toFixed(digits) : '—');

const where = (done, step) => `${Math.floor(done / DAYS_PER_YEAR)} years${done % DAYS_PER_YEAR || step ? ` ${done % DAYS_PER_YEAR} days` : ''}${step ? ` ${step} of ${perDay} steps` : ''} alone`;
const phase = ((forcingDay(days[doneBefore % LOOP_DAYS]) - 1 - (baseDay + doneBefore)) % DAYS_PER_YEAR + DAYS_PER_YEAR) % DAYS_PER_YEAR;
log(`--- ${new Date().toISOString()} ocean spin-up at N=${N} from ${from} (day ${saved.day}, ${where(doneBefore, stepBefore)}) under ${FORCING} days ${forcingDay(days[0])}–${forcingDay(days[days.length - 1])}, restoring ${RESTORE} W/m²/K, dt ${dt} s, ocean step ${oceanDt} s, after ${((performance.now() - t0) / 1000).toFixed(0)} s of setup${phase ? `; the forcing is ${phase} days out of season with the state` : ''}`);

async function readSurface() {
  const [surfaceT, ice] = await readRanges(gpu.device, gpu.buffers.S, [{ offset: S.TS, length: C }, { offset: S.ICE, length: C }]);
  const [concentration] = await readRanges(gpu.device, gpu.buffers.PH, [{ offset: PH.CONC, length: C }]);
  return { o: await model.ocean.serialize({ restart: true }), surfaceT, ice, concentration };
}

async function save(done, step, read = null) {
  await model.settle();
  const { o, surfaceT, ice, concentration } = read ?? await readSurface();
  const name = nameAt(done, step), landState = await model.land.serialize();
  writeFileSync(`${OUT}/${name}.partial`, encodeState({
    N, K: core.K, day: baseDay + done, time: baseTime + done * 86400 + step * dt, terrain: !!model.surfaceGeopotential,
    oceanYears: Math.floor(done / DAYS_PER_YEAR), oceanDays: done, ...(step ? { oceanStep: step } : {}),
    pi: saved.pi, theta: saved.theta, u: saved.u, surfaceT, q: saved.q, qc: saved.qc, ice, concentration, mlmSubsidence: model.radiation.mlmSubsidence,
    ocean: o, land: landState, sstByYear: Float64Array.from(history, (x) => x ?? NaN),
  }, { f64: ['sstByYear'] }));
  renameSync(`${OUT}/${name}.partial`, `${OUT}/${name}`);
  const kept = snapshots(), years = kept.filter((f) => !inYear(f));
  for (const old of years.slice(0, Math.max(0, years.length - KEEP))) unlinkSync(`${OUT}/${old}`);
  for (const old of kept) if (inYear(old) && old !== name) unlinkSync(`${OUT}/${old}`);
  log(`saved ${name} after ${((performance.now() - t0) / 60000).toFixed(1)} min; keeping ${snapshots().join(', ')}`);
  hook.after(`${OUT}/${name}`);
  hook.after(`${OUT}/${TAG}.log`);
}

async function yearEnd(done, last, seconds, daysRun) {
  await model.settle();
  const year = done / DAYS_PER_YEAR, read = await readSurface(), { o, surfaceT, ice, concentration } = read;
  const d = await model.oceanEngine.diagnostics();
  const sst = (cells) => meanOf(cells, (i) => o.T[i] - 273.15), classTop = (cells) => meanOf(cells, (i) => { let depth = 0; for (let k = 0; k < kT; k++) depth += o.h[k * C + i]; return depth; });
  log(`ocean after year ${year}: warm pool ${sst(warmPool).toFixed(1)} °C, cold tongue ${sst(coldTongue).toFixed(1)} °C (W−E ${(sst(westPacific) - sst(eastPacific)).toFixed(1)} K), ${THERMOCLINE_DENSITY} class top W Pac ${classTop(westPacific).toFixed(0)} m, E Pac ${classTop(eastPacific).toFixed(0)} m`);
  const column = BANDS.map(([top, bottom]) => {
    let sum = 0, weight = 0;
    for (const i of southern) {
      let above = 0;
      for (let k = 0; k < model.oceanEngine.layers; k++) {
        const h = o.h[k * C + i], overlap = Math.max(0, Math.min(bottom, above + h) - Math.max(top, above));
        sum += mesh.areaCell[i] * overlap * o.T[k * C + i]; weight += mesh.areaCell[i] * overlap; above += h;
      }
    }
    return `${top}–${bottom} m ${weight > 0 ? (sum / weight - 273.15).toFixed(2) : '—'}`;
  });
  let area = 0, drift = 0, square = 0;
  const extent = { north: 0, south: 0, recordedNorth: 0, recordedSouth: 0 };
  for (const i of ocean) {
    const a = mesh.areaCell[i], difference = (ice[i] > 0 ? FREEZING_POINT : surfaceT[i]) - last.fields.sst[i];
    area += a; drift += a * difference; square += a * difference * difference;
    const cover = ice[i] > 0 ? (concentration[i] > 0 ? concentration[i] : 1) : 0, north = mesh.latCell[i] > 0;
    if (cover >= 0.15) extent[north ? 'north' : 'south'] += a;
    if (last.fields.concentration[i] >= 0.15) extent[north ? 'recordedNorth' : 'recordedSouth'] += a;
  }
  const mean = globalSST(surfaceT, ice), back = Math.min(10, year);
  history[year] = mean;
  log(`year ${year}: Southern Ocean 60–70S ${column.join(', ')} °C; SST − recorded day ${last.day} mean ${(drift / area).toFixed(3)} K, rms ${Math.sqrt(square / area).toFixed(3)} K; global SST ${fixed(mean - 273.15, 2)} °C, drift ${fixed(mean - history[0], 3)} K since the start, ${fixed(mean - history[year - back], 3)} K over the last ${back} years; ice extent N ${(extent.north / 1e12).toFixed(1)} (recorded ${(extent.recordedNorth / 1e12).toFixed(1)}) S ${(extent.south / 1e12).toFixed(1)} (recorded ${(extent.recordedSouth / 1e12).toFixed(1)}) Mkm²; interior ${fixed(d.oceanInteriorT - 273.15, 2)} °C, currents ≤ ${fixed(d.oceanSpeed, 2)} m/s, transport ${fixed(d.oceanTransport, 0)} Sv, clamped ${d.oceanLimited}; ${seconds.toFixed(0)} s (${(seconds / daysRun).toFixed(2)} s a day)`);
  if (!Number.isFinite(mean) || !Number.isFinite(d.oceanSpeed) || !Number.isFinite(d.oceanInteriorT)) { log(`NaN in year ${year}; stopping`); await hook.drain(); process.exit(2); }
  await save(done, 0, read);
}

if (stop.requested) { log(`stopped by ${stop.requested} before the first step; nothing to save`); await hook.drain(); process.exit(0); }
let done = doneBefore, step = stepBefore, taken = 0, yearStart = performance.now(), daysRun = 0;
const halted = () => stop.requested || taken >= STOP_AFTER_STEPS;
while (done < YEARS * DAYS_PER_YEAR) {
  const last = await decodeForcing(new Uint8Array(readFileSync(`${FORCING}/${days[done % LOOP_DAYS]}`)));
  if (last.N !== N) throw new Error(`${days[done % LOOP_DAYS]} is N=${last.N}`);
  forced.setDay(last.fields, oceanDt);
  let oceanSteps = 0;
  while (step < perDay) {
    const stepped = forced.step(dt, RESTORE);
    step++; taken++;
    if (stepped && ++oceanSteps % 8 === 0) await model.settle();
    if (stepped && halted()) break;
  }
  if (step < perDay) break;
  step = 0;
  done++;
  daysRun++;
  if (done % DAYS_PER_YEAR === 0) {
    await yearEnd(done, last, (performance.now() - yearStart) / 1000, daysRun);
    yearStart = performance.now(); daysRun = 0;
  } else if ((SNAPSHOT_DAYS > 0 && (done % DAYS_PER_YEAR) % SNAPSHOT_DAYS === 0) || halted()) await save(done, 0);
  if (halted()) break;
}
if (step) await save(done, step);
if (halted()) log(`stopped by ${stop.requested ?? `STOP_AFTER_STEPS=${STOP_AFTER_STEPS}`} after ${where(done, step)}`);
await hook.drain();
process.exit(0);
