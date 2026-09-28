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
// carries), its time advanced by the years run, keeps the KEEP newest,
// and logs to OUT/<TAG>.log the ocean line of spinup.mjs, the Southern
// Ocean column at 60–70S, the SST and sea-ice extent against the
// recorded last day, and the wall time. It continues from the newest
// OUT/<TAG>_yearYYYY.bin when there is one, and otherwise starts from
// STATE, which should be a snapshot at the start or the end of the
// recorded cycle so that its season matches.
//
// Environment: N (64), STATE, FORCING, YEARS (1, the year to stop after),
// DAYS_PER_YEAR (365, the first that many days of FORCING make the year),
// RESTORE (30), TAG (ocean<N>), OUT (runs/), KEEP (2), OCEAN and RADIATION
// (JSON options as in spinup.mjs).
import { readFileSync, writeFileSync, readdirSync, renameSync, unlinkSync, appendFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { createGpuModel } from '../js/gpu/model.gpu.js';
import { decodeState, encodeState, savedLevels } from '../js/stateFile.module.js';
import { savedSubsidence } from '../js/physics/regrid.module.js';
import { readRanges } from '../js/gpu/device.module.js';
import { LAYER_DENSITIES, THERMOCLINE_DENSITY } from '../js/ocean/layered.module.js';
import { FREEZING_POINT } from '../js/physics/ice.module.js';
import { createForcedOcean } from '../js/gpu/forcing.gpu.js';
import { decodeForcing, forcingDay } from '../js/forcing.module.js';

const N = Number(process.env.N ?? 64), TAG = process.env.TAG ?? `ocean${N}`, YEARS = Number(process.env.YEARS ?? 1), KEEP = Number(process.env.KEEP ?? 2);
const DAYS_PER_YEAR = Number(process.env.DAYS_PER_YEAR ?? 365), RESTORE = Number(process.env.RESTORE ?? 30);
const OUT = process.env.OUT ?? new URL('../runs/', import.meta.url).pathname, FORCING = process.env.FORCING;
const OCEAN = JSON.parse(process.env.OCEAN ?? '{}'), RADIATION = JSON.parse(process.env.RADIATION ?? '{}');
if (!FORCING) throw new Error('FORCING must name a directory of recorded forcing');
const log = (line) => { console.log(line); appendFileSync(`${OUT}/${TAG}.log`, line + '\n'); };
const yearOf = (file) => Number(file.match(/_year(\d+)\.bin$/)[1]);
const snapshots = () => readdirSync(OUT).filter((f) => f.startsWith(`${TAG}_year`) && /_year\d+\.bin$/.test(f)).sort((a, b) => yearOf(a) - yearOf(b));

const t0 = performance.now();
const days = readdirSync(FORCING).filter((f) => forcingDay(f) !== null).sort((a, b) => forcingDay(a) - forcingDay(b)).slice(0, DAYS_PER_YEAR);
if (days.length < DAYS_PER_YEAR) throw new Error(`${FORCING} holds ${days.length} days of forcing, not ${DAYS_PER_YEAR}`);
const existing = snapshots();
const from = existing.length ? `${OUT}/${existing[existing.length - 1]}` : process.env.STATE;
if (!from) throw new Error(`no ${TAG}_yearYYYY.bin in ${OUT} and no STATE to start from`);
const saved = await decodeState(new Uint8Array(readFileSync(from)));
if (saved.N !== N) throw new Error(`${from} is N=${saved.N}`);
const firstYear = existing.length ? yearOf(existing[existing.length - 1]) : 0;

const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const model = await createGpuModel(new Grid(N), { topography, ocean: OCEAN, radiation: RADIATION, levels: savedLevels(saved) });
const { mesh, core, state, gpu } = model;
const C = mesh.nCells, dt = 1350 * 16 / N, perDay = Math.round(86400 / dt), oceanDt = model.oceanEngine.everySteps * dt;
['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'].forEach((name, a) => state[a].set(saved[name]));
model.seaIce.load(state[6], saved.concentration ?? null);
model.radiation.mlmSubsidence.set(savedSubsidence(saved, model));
model.time = saved.time;
model.load();
model.ocean.load(saved.ocean, state[3], state[6]);
model.land.load({ soil: Float64Array.from(saved.land.soil), snow: Float64Array.from(saved.land.snow), ...(saved.land.vegetation ? { vegetation: Float64Array.from(saved.land.vegetation) } : {}), ...(saved.land.surface ? { surface: Float64Array.from(saved.land.surface) } : {}) });
const forced = createForcedOcean(model);
const phase = ((forcingDay(days[0]) - 1 - saved.day) % DAYS_PER_YEAR + DAYS_PER_YEAR) % DAYS_PER_YEAR;
log(`--- ${new Date().toISOString()} ocean spin-up at N=${N} from ${from} (day ${saved.day}, year ${firstYear}) under ${FORCING} days ${forcingDay(days[0])}–${forcingDay(days[days.length - 1])}, restoring ${RESTORE} W/m²/K, dt ${dt} s, ocean step ${oceanDt} s, after ${((performance.now() - t0) / 1000).toFixed(0)} s of setup${phase ? `; the forcing starts ${phase} days out of season with the state` : ''}`);

const deg = 180 / Math.PI, land = model.geography.land;
const kT = LAYER_DENSITIES.findIndex((r) => r >= THERMOCLINE_DENSITY);
const sea = (a, b, c, e) => [...Array(C).keys()].filter((i) => !land[i] && mesh.latCell[i] * deg >= a && mesh.latCell[i] * deg <= b && (c <= e ? mesh.lonCell[i] * deg >= c && mesh.lonCell[i] * deg <= e : mesh.lonCell[i] * deg >= c || mesh.lonCell[i] * deg <= e));
const meanOf = (cells, f) => cells.reduce((s, i) => s + f(i), 0) / Math.max(1, cells.length);
const warmPool = sea(-10, 10, 120, 160), coldTongue = sea(-2, 2, -110, -90), westPacific = sea(-5, 5, 140, 170), eastPacific = sea(-5, 5, -120, -90);
const southern = sea(-70, -60, -180, 180), ocean = sea(-90, 90, -180, 180);
const BANDS = [[0, 60], [60, 200], [200, 500], [500, 1000]];
const S = gpu.layout.S, PH = gpu.layout.PH;

let year = firstYear, last = null;
while (year < YEARS) {
  const start = performance.now();
  for (const file of days) {
    last = await decodeForcing(new Uint8Array(readFileSync(`${FORCING}/${file}`)));
    if (last.N !== N) throw new Error(`${file} is N=${last.N}`);
    forced.setDay(last.fields, oceanDt);
    let oceanSteps = 0;
    for (let n = 0; n < perDay; n++) if (forced.step(dt, RESTORE) && ++oceanSteps % 8 === 0) await model.settle();
  }
  await model.settle();
  year++;
  const seconds = (performance.now() - start) / 1000;
  model.time = saved.time + (year - firstYear) * DAYS_PER_YEAR * 86400;
  const o = await model.ocean.serialize();
  const [surfaceT, ice] = await readRanges(gpu.device, gpu.buffers.S, [{ offset: S.TS, length: C }, { offset: S.ICE, length: C }]);
  const [concentration] = await readRanges(gpu.device, gpu.buffers.PH, [{ offset: PH.CONC, length: C }]);
  const sst = (cells) => meanOf(cells, (i) => o.T[i] - 273.15), classTop = (cells) => meanOf(cells, (i) => { let d = 0; for (let k = 0; k < kT; k++) d += o.h[k * C + i]; return d; });
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
  log(`year ${year}: Southern Ocean 60–70S ${column.join(', ')} °C; SST − recorded day ${last.day} mean ${(drift / area).toFixed(3)} K, rms ${Math.sqrt(square / area).toFixed(3)} K; ice extent N ${(extent.north / 1e12).toFixed(1)} (recorded ${(extent.recordedNorth / 1e12).toFixed(1)}) S ${(extent.south / 1e12).toFixed(1)} (recorded ${(extent.recordedSouth / 1e12).toFixed(1)}) Mkm²; ${seconds.toFixed(0)} s (${(seconds / DAYS_PER_YEAR).toFixed(2)} s a day)`);
  const name = `${TAG}_year${String(year).padStart(4, '0')}.bin`;
  const landState = await model.land.serialize();
  writeFileSync(`${OUT}/${name}.partial`, encodeState({
    N, K: core.K, day: Math.round(model.time / 86400), time: model.time, terrain: !!model.surfaceGeopotential, levels: core.levels, oceanYears: year,
    pi: saved.pi, theta: saved.theta, u: saved.u, surfaceT, q: saved.q, qc: saved.qc, ice, concentration, mlmSubsidence: model.radiation.mlmSubsidence,
    ocean: { h: o.h, u: o.u, T: o.T, S: o.S, eta: o.eta }, land: landState,
  }));
  renameSync(`${OUT}/${name}.partial`, `${OUT}/${name}`);
  const kept = snapshots();
  for (const old of kept.slice(0, Math.max(0, kept.length - KEEP))) unlinkSync(`${OUT}/${old}`);
  log(`saved ${name} after ${((performance.now() - t0) / 60000).toFixed(1)} min; keeping ${snapshots().join(', ')}`);
}
process.exit(0);
