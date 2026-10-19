// The hybrid ocean (vertical: 'hybrid' in js/ocean/layered.module.js)
// alone on the CPU under a recorded year of surface forcing
// (js/forcing.module.js, recorded by scripts/spinup.mjs RECORD at its own
// N), the fields interpolated onto this N's mesh (js/physics/regrid.module.js),
// from the WOA23 climatology (data/woa_annual_1deg.bin). The sea ice and
// its concentration are prescribed from the recording and kept from the
// ocean, which would hold its surface at freezing against an unlimited
// sink: under ice the surface is instead relaxed toward freezing at
// RESTORE_ICE W/m²/K, the ice-ocean exchange of a mixed layer under pack,
// and the recorded ice's growth and melt salt and freshen the surface as
// freshwater taken or given at the ice's density less its salinity. The
// recorded stress already carries the ice transmission. Over open
// water each ocean step applies the recorded net flux less the sunlight
// that passes the top layer (handed below as the coupled model hands it),
// the restoring −RESTORE·(SST − target) toward WOA23's monthly SST on the
// forcing day's date (RESTORE_TO woa, the default; the recorded SST where
// WOA23 has none) or toward the recorded SST (RESTORE_TO recorded), and
// the recorded evaporation, rain and runoff, and holds the surface no
// colder than freezing. The recording holds daily means; with DIURNAL (on by
// default) each cell's shortwave is replayed as a half-sine over its
// local day, 12 hours long and of the recorded mean, and the rest of the
// net flux as the mean, so the wind has the day's stratification to
// erode as it does under the coupled model's sun. The recorded stress, which carries the ice transmission
// already, is undone for it so the ocean applies its own.
//
// Logs OUT/<TAG>.log every LOG_DAYS days: the global and tropical SST
// against the recording, the mixed-layer depths, the 140W equatorial
// column, the Pacific undercurrent's core (the fastest eastward current
// between 30 and 400 m within 2° of the equator), the Southern Ocean column, the ocean
// diagnostics and the pace. Saves OUT/<TAG>_dayDDD.bin (stateFile.module.js
// format, kind 'hybridOcean': the ocean's serialize() arrays, the surface
// temperature and the sunlight not yet deposited) every SAVE_DAYS days
// and at the end, continuing from the newest such file in OUT when there
// is one.
//
// Environment: N (32), FORCING (runs/twelve-local/forcing), OUT
// (runs/hybrid), TAG (hybrid<N>), DAYS (365), STEP (the ocean step in
// seconds; the coupled cadence of js/cadence.module.js by default),
// RESTORE (30 W/m²/K), RESTORE_ICE (100 W/m²/K), RESTORE_TO (woa), DIURNAL (1), LOG_DAYS (10),
// SAVE_DAYS (30), OCEAN (JSON options for createOcean), WORKERS (the
// machine's parallelism less one; 0 or 1 runs the ocean on this thread),
// ENGINE (cpu; gpu runs the same replay on the GPU hybrid ocean of
// js/gpu/layeredOcean.gpu.js, the surface budget still on this thread,
// with the surface temperature and the sunlight not yet deposited
// handed to the device and read back each ocean step; one GPU model at
// a time on this Mac), STRESS_SCALE (1; multiplies the recorded stress,
// an instrument for the recording's wind bias).
import { readFileSync, writeFileSync, readdirSync, appendFileSync, mkdirSync, existsSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { gridMesh, stateGrid } from '../js/stateGrid.module.js';
import { createGeography, topographyFromInt16 } from '../js/geography.module.js';
import { createOcean, THERMOCLINE_DENSITY } from '../js/ocean/layered.module.js';
import { layerDensities, mixedLayerDepth } from '../js/ocean/column.module.js';
import { loadClimatology, CLIMATOLOGY_FILE, decodeMonthlySst, monthWeights, climatologyAt, MONTHLY_SST_FILE } from '../js/ocean/climatology.module.js';
import { decodeForcing, forcingDay } from '../js/forcing.module.js';
import { decodeState, encodeState } from '../js/stateFile.module.js';
import { interpolationWeights, regridCellField, regridEdgeField } from '../js/physics/regrid.module.js';
import { FREEZING_POINT } from '../js/physics/ice.module.js';
import { oceanMinutes } from '../js/cadence.module.js';
import { stopOnSignal } from './runControl.mjs';
import { parallelism } from '../js/threads.module.js';
import { createParallelModel } from '../js/parallel.module.js';

const N = Number(process.env.N ?? 32), FORCING = process.env.FORCING ?? 'runs/twelve-local/forcing', OUT = process.env.OUT ?? 'runs/hybrid', TAG = process.env.TAG ?? `hybrid${N}`;
const DAYS = Number(process.env.DAYS ?? 365), RESTORE = Number(process.env.RESTORE ?? 30), RESTORE_ICE = Number(process.env.RESTORE_ICE ?? 100), LOG_DAYS = Number(process.env.LOG_DAYS ?? 10), SAVE_DAYS = Number(process.env.SAVE_DAYS ?? 30);
const STEP = Number(process.env.STEP ?? oceanMinutes(N) * 60), OCEAN = JSON.parse(process.env.OCEAN ?? '{}'), DIURNAL = (process.env.DIURNAL ?? '1') !== '0';
const STRESS_SCALE = Number(process.env.STRESS_SCALE ?? 1);
const RESTORE_TO = process.env.RESTORE_TO ?? 'woa', EQUINOX_DAY_OF_YEAR = 78, WORKERS = Number(process.env.WORKERS ?? Math.max(1, await parallelism() - 1));
if (!['woa', 'recorded'].includes(RESTORE_TO)) throw new Error(`RESTORE_TO is woa or recorded, not ${RESTORE_TO}`);
if (86400 % STEP) throw new Error(`STEP ${STEP} s does not divide a day`);
if (!existsSync(OUT)) mkdirSync(OUT, { recursive: true });
const log = (line) => { console.log(line); appendFileSync(`${OUT}/${TAG}.log`, line + '\n'); };
const stop = stopOnSignal(log);
const fixed = (x, n) => (Number.isFinite(x) ? x.toFixed(n) : String(x));

const t0 = performance.now();
const files = readdirSync(FORCING).filter((f) => forcingDay(f) !== null).sort((a, b) => forcingDay(a) - forcingDay(b));
if (!files.length) throw new Error(`${FORCING} holds no forcing days`);
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const first = await decodeForcing(new Uint8Array(readFileSync(`${FORCING}/${files[0]}`)));
const source = { mesh: stateGrid(first, topography).mesh };
source.geography = createGeography(source.mesh, topography);
const mesh = source.mesh.nCells === 10 * N * N + 2 ? source.mesh : gridMesh(new Grid(N)), C = mesh.nCells, E = mesh.nEdges;
const geography = mesh === source.mesh ? source.geography : createGeography(mesh, topography), target = { mesh, geography };
const same = mesh === source.mesh;
const atCells = same ? null : interpolationWeights(source.mesh, mesh.xCell), atEdges = same ? null : interpolationWeights(source.mesh, mesh.xEdge);
const seaMask = Uint8Array.from(source.geography.land, (l) => (l ? 0 : 1)), landMask = Uint8Array.from(source.geography.land, (l) => (l ? 1 : 0));
const SEA_FIELDS = ['netFlux', 'shortwave', 'sst', 'ice', 'concentration'], ANY_FIELDS = ['evaporation', 'rain'];
function regrid(fields) {
  if (same) return Object.fromEntries(Object.entries(fields).map(([k, v]) => [k, Float64Array.from(v)]));
  const out = {};
  for (const name of SEA_FIELDS) out[name] = regridCellField(source, target, Float64Array.from(fields[name]), atCells, seaMask);
  for (const name of ANY_FIELDS) out[name] = regridCellField(source, target, Float64Array.from(fields[name]), atCells);
  out.runoff = regridCellField(source, target, Float64Array.from(fields.runoff), atCells, landMask);
  out.stress = regridEdgeField(source, target, Float64Array.from(fields.stress), atEdges);
  for (let i = 0; i < C; i++) { if (geography.land[i]) { out.ice[i] = 0; out.concentration[i] = 0; out.runoff[i] = Math.max(0, out.runoff[i]); } else { out.runoff[i] = 0; out.ice[i] = Math.max(0, out.ice[i]); out.concentration[i] = Math.min(1, Math.max(0, out.concentration[i])); if (out.ice[i] > 0 && out.concentration[i] <= 0) out.concentration[i] = 1; } }
  for (let e = 0; e < E; e++) if (!geography.edgeOcean[e]) out.stress[e] = 0;
  return out;
}

const climatologyPath = new URL(`../${CLIMATOLOGY_FILE}`, import.meta.url).pathname, oceanOptions = { vertical: 'hybrid', everySteps: 1, penetration: 'IB', ...OCEAN };
const ENGINE = process.env.ENGINE ?? 'cpu';
if (ENGINE !== 'cpu' && ENGINE !== 'gpu') throw new Error(`ENGINE is cpu or gpu, not ${ENGINE}`);
const parallel = ENGINE === 'cpu' && WORKERS > 1 ? await createParallelModel(mesh, { topography, ocean: { climatology: climatologyPath, ...oceanOptions } }, WORKERS) : null;
if (parallel && parallel.mesh !== mesh) throw new Error('the parallel model built its own mesh');
const gpuModel = ENGINE === 'gpu' ? await (await import('../js/gpu/model.gpu.js')).createGpuModel(mesh, { topography, ocean: { climatology: climatologyPath, ...oceanOptions } }) : null;
if (gpuModel && gpuModel.mesh !== mesh) throw new Error('the GPU model built its own mesh');
const readRanges = gpuModel ? (await import('../js/gpu/device.module.js')).readRanges : null;
const ocean = gpuModel ? gpuModel.oceanEngine : parallel ? parallel.ocean : createOcean(mesh, { geography, climatology: await loadClimatology(climatologyPath), ...oceanOptions });
const woa = RESTORE_TO === 'woa' ? decodeMonthlySst(readFileSync(new URL(`../${MONTHLY_SST_FILE}`, import.meta.url))) : null;
function targetSst(modelDay, recorded) {
  const out = Float64Array.from(recorded);
  if (!woa) return out;
  const months = monthWeights(woa, (((modelDay + EQUINOX_DAY_OF_YEAR) % 365) + 365) % 365);
  for (let i = 0; i < C; i++) { if (geography.land[i]) continue; const t = climatologyAt(woa, months, mesh.latCell[i], mesh.lonCell[i]); if (Number.isFinite(t)) out[i] = t; }
  return out;
}
const L = ocean.layers, at = (k, i) => k * C + i;
const surfaceT = new Float64Array(C), flux = new Float64Array(C), stress = new Float64Array(E), rainAmount = new Float64Array(C), runoffAmount = new Float64Array(C), noIce = new Float64Array(C), evaporation = new Float64Array(C), previousVolume = new Float64Array(C);
const ICE_DENSITY = 917, ICE_SALINITY = 5;
const deg = 180 / Math.PI, land = geography.land;

/*
 * The ocean behind one interface for both engines: the CPU ocean's arrays
 * are read in place; the GPU's are mirrored at each day's end, and each
 * ocean step hands the device the surface temperature and the sunlight
 * not yet deposited and reads them back with the capacity and the
 * transmission of its top layer.
 */
const mirror = gpuModel ? { h: new Float64Array(L * C), u: new Float64Array(L * E), T: new Float64Array(L * C), S: new Float64Array(L * C), capacity: new Float64Array(C), transmission: new Float64Array(C), penetrated: new Float64Array(C), none: new Float64Array(C), hc: new Float64Array(L), Tc: new Float64Array(L), Sc: new Float64Array(L), rho: new Float64Array(L) } : null;
const PHL = gpuModel ? gpuModel.gpu.layout.PH : null, SL = gpuModel ? gpuModel.gpu.layout.S : null, device = gpuModel ? gpuModel.gpu.device : null, buffers = gpuModel ? gpuModel.gpu.buffers : null;
async function pullSurface() {
  const [capacity, transmission, penetrated] = await readRanges(device, buffers.PH, [{ offset: PHL.CAP, length: C }, { offset: PHL.TRANS, length: C }, { offset: PHL.PEN, length: C }]);
  mirror.capacity.set(capacity); mirror.transmission.set(transmission); mirror.penetrated.set(penetrated);
}
async function pullState() {
  const s = await ocean.serialize();
  mirror.h.set(s.h); mirror.u.set(s.u); mirror.T.set(s.T); mirror.S.set(s.S);
}
const io = gpuModel ? {
  h: mirror.h, u: mirror.u, capacity: mirror.capacity, transmission: mirror.transmission, penetrated: mirror.penetrated,
  layerT: (k, i) => mirror.T[at(k, i)], surfaceS: (i) => mirror.S[i],
  mixedDepths() {
    const out = new Float64Array(C);
    for (let i = 0; i < C; i++) {
      if (land[i]) continue;
      let K = 0;
      for (let k = 0; k < L; k++) { const hk = mirror.h[at(k, i)]; if (!(hk > 0.01)) continue; mirror.hc[K] = hk; mirror.Tc[K] = mirror.T[at(k, i)]; mirror.Sc[K] = mirror.S[at(k, i)]; K++; }
      layerDensities(mirror.Tc, mirror.Sc, K, mirror.rho);
      out[i] = mixedLayerDepth(mirror.hc, mirror.rho, K);
    }
    return out;
  },
  async initialize() { gpuModel.state[3].set(surfaceT); gpuModel.load({ ocean: false }); const start = ocean.initialize(surfaceT, noIce); await pullSurface(); await pullState(); return start; },
  async load(saved) {
    gpuModel.state[3].set(surfaceT); gpuModel.load({ ocean: false });
    ocean.upload({ h: saved.h, u: saved.u, T: saved.T, S: saved.S, eta: saved.eta, densities: Array.from(saved.densities), vertical: 'hybrid' }, surfaceT, noIce);
    if (saved.penetrated) device.queue.writeBuffer(buffers.PH, 4 * PHL.PEN, Float32Array.from(saved.penetrated));
    await pullSurface(); await pullState();
  },
  setDay() {
    device.queue.writeBuffer(buffers.PH, 4 * PHL.EVAP, Float32Array.from(evaporation));
    device.queue.writeBuffer(buffers.PH, 4 * PHL.RAIN, Float32Array.from(rainAmount));
    device.queue.writeBuffer(buffers.PH, 4 * PHL.RUNOFF, Float32Array.from(runoffAmount));
  },
  async step(time) {
    ocean.forgetAccumulated('rain'); ocean.forgetAccumulated('runoff');
    ocean.accumulateFreshwater(STEP);
    device.queue.writeBuffer(buffers.PH, 4 * PHL.PEN, Float32Array.from(mirror.penetrated));
    await ocean.advance(surfaceT, noIce, stress, STEP, mirror.none, time);
    const [[ts], [pen, capacity, transmission]] = await Promise.all([readRanges(device, buffers.S, [{ offset: SL.TS, length: C }]), readRanges(device, buffers.PH, [{ offset: PHL.PEN, length: C }, { offset: PHL.CAP, length: C }, { offset: PHL.TRANS, length: C }])]);
    surfaceT.set(ts); mirror.penetrated.set(pen); mirror.capacity.set(capacity); mirror.transmission.set(transmission);
  },
  async endDay() { await pullState(); },
  diagnostics: () => ocean.diagnostics(),
  async serialize() { return { ...(await ocean.serialize()), penetrated: mirror.penetrated }; },
  close: () => gpuModel.settle(),
  engine: 'the GPU',
} : {
  h: ocean.h, u: ocean.u, capacity: ocean.capacity, transmission: ocean.transmission, penetrated: ocean.penetrated,
  layerT: (k, i) => ocean.Q[at(k, i)] / ocean.h[at(k, i)], surfaceS: (i) => ocean.W[i] / ocean.h[i],
  mixedDepths: () => ocean.fields().h1,
  async initialize() { return ocean.initialize(surfaceT, noIce); },
  async load(saved) { ocean.load({ h: saved.h, u: saved.u, T: saved.T, S: saved.S, eta: saved.eta, densities: Array.from(saved.densities) }, surfaceT, noIce); if (saved.penetrated) ocean.penetrated.set(saved.penetrated); },
  setDay() {},
  async step(time) { ocean.accumulate(evaporation, rainAmount, STEP, runoffAmount); ocean.advance(surfaceT, noIce, flux, stress, STEP, undefined, time); },
  async endDay() {},
  diagnostics: async () => ocean.diagnostics(),
  async serialize() { return { ...ocean.serialize(), penetrated: ocean.penetrated }; },
  async close() { if (parallel) await parallel.close(); },
  engine: parallel ? `${WORKERS} worker threads` : 'one thread',
};

const sea = (a, b, c, e) => [...Array(C).keys()].filter((i) => !land[i] && mesh.latCell[i] * deg >= a && mesh.latCell[i] * deg <= b && mesh.lonCell[i] * deg >= c && mesh.lonCell[i] * deg <= e);
const all = sea(-90, 90, -180, 180), tropics = sea(-20, 20, -180, 180), southern = sea(-70, -60, -180, 180), band = [...sea(-2, 2, 140, 180), ...sea(-2, 2, -180, -85)];
const areaOf = (cells) => cells.reduce((s, i) => s + mesh.areaCell[i], 0);
const meanOf = (cells, f) => cells.reduce((s, i) => s + mesh.areaCell[i] * f(i), 0) / areaOf(cells);
const nearest = (lat, lon) => all.reduce((best, i) => { const d = Math.hypot(mesh.latCell[i] * deg - lat, mesh.lonCell[i] * deg - lon); return d < best.d ? { i, d } : best; }, { i: -1, d: Infinity }).i;
const station = nearest(0, -140);
const eastAt = (i) => { const x = mesh.xCell[3 * i], y = mesh.xCell[3 * i + 1], r = Math.hypot(x, y); return [-y / r, x / r, 0]; };
function eastward(k, i) {
  const east = eastAt(i);
  let x = 0, y = 0, z = 0;
  for (let m = 0; m < mesh.nEdgesOnCell[i]; m++) {
    const e = mesh.edgesOnCell[mesh.maxEdges * i + m], w = 0.5 * mesh.dcEdge[e] * mesh.dvEdge[e] * io.u[k * E + e];
    x += w * mesh.nEdge[3 * e]; y += w * mesh.nEdge[3 * e + 1]; z += w * mesh.nEdge[3 * e + 2];
  }
  return (x * east[0] + y * east[1] + z * east[2]) / mesh.areaCell[i];
}
const column = (i, every) => { const parts = []; let z = 0; for (let k = 0; k < L; k++) { const hk = io.h[at(k, i)]; if (k % every === 0 || k < 4) parts.push(`${(z + 0.5 * hk).toFixed(0)}m ${(io.layerT(k, i) - 273.15).toFixed(2)}`); z += hk; if (z > 1500) break; } return parts.join(', '); };

const saves = readdirSync(OUT).filter((f) => f.startsWith(`${TAG}_day`) && f.endsWith('.bin')).sort();
let dayDone = 0, time = 0;
const ice = new Float64Array(C), concentration = new Float64Array(C);
let day0 = regrid(first.fields);
day0.target = targetSst(first.day, day0.sst);
ice.set(day0.ice); concentration.set(day0.concentration);
if (saves.length) {
  const saved = await decodeState(new Uint8Array(readFileSync(`${OUT}/${saves[saves.length - 1]}`)));
  if (saved.N !== N) throw new Error(`${saves[saves.length - 1]} is N=${saved.N}`);
  dayDone = saved.day; time = saved.time;
  const forcingToday = await decodeForcing(new Uint8Array(readFileSync(`${FORCING}/${files[dayDone % files.length]}`)));
  day0 = regrid(forcingToday.fields);
  day0.target = targetSst(forcingToday.day, day0.sst);
  ice.set(day0.ice); concentration.set(day0.concentration);
  surfaceT.set(saved.surfaceT);
  await io.load(saved);
  log(`continuing from ${saves[saves.length - 1]} at day ${dayDone} on ${io.engine}`);
} else {
  surfaceT.set(day0.sst);
  const start = await io.initialize();
  log(`N=${N}: ${C} cells, ${E} edges, ${L} layers, ocean step ${STEP} s, forcing ${files.length} days from ${FORCING}${same ? '' : ` (N=${first.N}, interpolated)`}, started from the atlas (${start ? `${start.atlas} cells, ${start.analytic} from their neighbors` : 'no atlas'}), restoring ${RESTORE} W/m²/K to ${RESTORE_TO === 'woa' ? "WOA23's monthly SST on the forcing day's date" : 'the recorded SST'}, ice prescribed with ${RESTORE_ICE} W/m²/K toward freezing beneath it, shortwave ${DIURNAL ? 'as a half-sine over each local day' : 'as the daily mean'}, ${io.engine}`);
}

async function save(day) {
  const s = await io.serialize();
  const bytes = encodeState({ kind: 'hybridOcean', N, day, time, h: Float32Array.from(s.h), u: Float32Array.from(s.u), T: Float32Array.from(s.T), S: Float32Array.from(s.S), eta: Float32Array.from(s.eta), densities: Float64Array.from(s.densities), surfaceT: Float32Array.from(surfaceT), penetrated: Float32Array.from(s.penetrated) });
  const name = `${TAG}_day${String(day).padStart(3, '0')}.bin`;
  writeFileSync(`${OUT}/${name}`, bytes);
  return name;
}

async function report(day, fields, seconds) {
  const d = await io.diagnostics(), h1 = io.mixedDepths();
  const sst = meanOf(all, (i) => surfaceT[i]), drift = meanOf(all, (i) => surfaceT[i] - fields.target[i]), tropicDrift = meanOf(tropics, (i) => surfaceT[i] - fields.target[i]);
  const rms = Math.sqrt(meanOf(all, (i) => (surfaceT[i] - fields.target[i]) ** 2)), recorded = meanOf(all, (i) => surfaceT[i] - fields.sst[i]), tropicRecorded = meanOf(tropics, (i) => surfaceT[i] - fields.sst[i]);
  let core = 0, coreDepth = 0, coreCell = -1;
  for (const i of band) { let z = 0; for (let k = 0; k < L; k++) { const hk = io.h[at(k, i)], u = eastward(k, i); if (z + 0.5 * hk > 30 && u > core) { core = u; coreDepth = z + 0.5 * hk; coreCell = i; } z += hk; if (z > 400) break; } }
  const mixed = (cells) => meanOf(cells, (i) => h1[i]);
  log(`day ${day}: SST ${fixed(sst - 273.15, 2)} °C, − ${RESTORE_TO === 'woa' ? 'WOA' : 'recorded'} ${fixed(drift, 3)} K (20S–20N ${fixed(tropicDrift, 3)}), rms ${fixed(rms, 3)}${RESTORE_TO === 'woa' ? `, − recorded ${fixed(recorded, 3)} K (20S–20N ${fixed(tropicRecorded, 3)})` : ''}; mixed layer ${fixed(mixed(all), 1)} m (tropics ${fixed(mixed(tropics), 1)}, 60–70S ${fixed(mixed(southern), 1)}); undercurrent core ${fixed(core, 3)} m/s at ${coreDepth.toFixed(0)} m${coreCell >= 0 ? ` (${(mesh.lonCell[coreCell] * deg).toFixed(0)}E)` : ''}; thermocline ${fixed(d.oceanThermoclineDepth, 0)} m; interior ${fixed(d.oceanInteriorT - 273.15, 2)} °C; currents ≤ ${fixed(d.oceanSpeed, 2)} m/s, transport ${fixed(d.oceanTransport, 0)} Sv, clamped ${d.oceanLimited}, dried ${d.oceanDried}; ${(seconds / Math.max(1, LOG_DAYS)).toFixed(0)} s a day`);
  log(`  140W 0N: ${column(station, 4)}`);
  log(`  60–70S column: ${column(southern[Math.floor(southern.length / 2)] ?? station, 6)}`);
}

const perDay = Math.round(86400 / STEP);
let lastLog = performance.now();
for (let day = dayDone; day < DAYS; day++) {
  let forcing = day0;
  if (day !== dayDone) {
    const decoded = await decodeForcing(new Uint8Array(readFileSync(`${FORCING}/${files[day % files.length]}`)));
    forcing = regrid(decoded.fields);
    forcing.target = targetSst(decoded.day, forcing.sst);
  }
  ice.set(forcing.ice); concentration.set(forcing.concentration);
  for (let i = 0; i < C; i++) {
    rainAmount[i] = forcing.rain[i] * STEP; runoffAmount[i] = forcing.runoff[i] * STEP;
    const volume = land[i] ? 0 : forcing.ice[i] * forcing.concentration[i], growth = day === dayDone ? 0 : (volume - previousVolume[i]) / 86400;
    previousVolume[i] = volume;
    const salinity = land[i] ? 35 : Math.max(1, io.surfaceS(i));
    evaporation[i] = forcing.evaporation[i] + growth * ICE_DENSITY * (1 - ICE_SALINITY / salinity);
  }
  for (let e = 0; e < E; e++) stress[e] = STRESS_SCALE * forcing.stress[e];
  io.setDay();
  for (let s = 0; s < perDay; s++) {
    const middle = ((time + 0.5 * STEP) % 86400) / 3600;
    for (let i = 0; i < C; i++) {
      if (land[i]) continue;
      if (ice[i] > 0) { surfaceT[i] = Math.max(FREEZING_POINT, surfaceT[i] - STEP * RESTORE_ICE * (surfaceT[i] - FREEZING_POINT) / io.capacity[i]); continue; }
      let sun = forcing.shortwave[i];
      if (DIURNAL) { const local = (middle + mesh.lonCell[i] * deg / 15 + 48) % 24; sun *= local > 6 && local < 18 ? Math.PI * Math.sin(Math.PI * (local - 6) / 12) : 0; }
      const pen = sun * io.transmission[i];
      io.penetrated[i] += STEP * pen;
      let T = surfaceT[i];
      T += STEP * (forcing.netFlux[i] - forcing.shortwave[i] + sun - RESTORE * (T - forcing.target[i]) - pen) / io.capacity[i];
      surfaceT[i] = Math.max(FREEZING_POINT, T);
    }
    await io.step(time + STEP);
    time += STEP;
  }
  await io.endDay();
  if (!all.every((i) => Number.isFinite(surfaceT[i]))) { log(`NaN in the SST on day ${day + 1}`); await save(day + 1); process.exit(2); }
  const done = day + 1;
  if (done % LOG_DAYS === 0 || done === DAYS) { await report(done, forcing, (performance.now() - lastLog) / 1000); lastLog = performance.now(); }
  if (done % SAVE_DAYS === 0 || done === DAYS || stop.requested) { log(`saved ${await save(done)} after ${((performance.now() - t0) / 60000).toFixed(1)} min`); }
  if (stop.requested) break;
}
await io.close();
