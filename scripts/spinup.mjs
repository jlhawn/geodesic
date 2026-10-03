// One spin-up segment on the GPU engine with the page's configuration:
// continue from the newest <TAG>_dayDDDD.bin in OUT (or start from the
// fresh initial state when there is none), step until MINUTES of wall
// time have passed or DAYS is reached, finishing the simulated day, save a
// binary snapshot and keep the KEEP newest. Logs one line a day to
// <TAG>.log, with the sea-ice extent of each hemisphere (the area of the
// cells at least 15% covered), and at the end of the segment the rain,
// vegetation and surface temperature of the regions in BOXES, the
// equatorial Pacific's surface temperatures and thermocline, and the
// convection and equator lines of js/audit.module.js, and exits with 2 on
// NaN. The daily line's ASR (with its part absorbed in the atmosphere),
// OLR and albedo are means over the day's steps, the albedo the day's
// reflected over its incoming sunlight: the last step alone sees the sun
// at the same UTC hour every day, over whatever cloud lies under it then,
// which between runs with different cloud geography moves the albedo by
// up to 0.04 and the fluxes by up to 15 W/m². The radiation's clear-sky
// pass is on unless RADIATION sets clearSkyPass false, and the line gives
// the day-mean shortwave and longwave cloud effects after the albedo
// (SWCRE, ASR less clear-sky ASR; LWCRE, clear-sky OLR less OLR) and the
// clear-sky reflectance (the day's clear-sky reflected over its incoming
// sunlight), then the day-mean sunlight absorbed and net longwave
// (downward less emitted) at the surface of the sea cells (sea ice
// included), and the sunlight absorbed at the surface of the sea cells
// poleward of 60° that hold ice at the day's end, leads included, in each
// hemisphere. The
// state saved carries the last day's per-cell convective and large-scale
// rain, absorbed sunlight, outgoing longwave, albedo and cloud effects,
// the boundary layer's depth, mixing top, regime and surface buoyancy
// flux, what the next step reads from the last before it computes its
// own (the lowest layer's wind speed, the evaporation and the shallow
// cumulus) and the ocean's restart arrays, so that a run continued from
// a day's snapshot steps on bit for bit as the run that did not stop.
//
// SIGTERM or SIGINT stops the segment after the ocean step in progress
// and exits 0: at a day's end it saves <TAG>_dayDDDD.bin as usual, inside
// a day <TAG>_dayDDDD_stepSSSS.bin (DDDD days and SSSS steps done, with
// the running rain and runoff totals the ocean takes its freshwater from
// as differences and, when RECORD is set, the forcing recorder's part of
// the day), which the next segment continues from and deletes once it has
// saved a whole day.
// Whole-day names alone count towards KEEP, and only whole-day names match
// the day patterns of the shell drivers.
//
// Environment: N (128), TAG (spin<N>), MINUTES (15), DAYS (none), KEEP (2), OUT
// (runs/), OCEAN (JSON options for the ocean, e.g. '{"closureHours":3}'),
// RADIATION (JSON options for the radiation, e.g. '{"cloudSolarAbsorption":0}'),
// MOIST (JSON options for the moist physics, e.g. '{"plumeEntrainment":0.15}'),
// BOUNDARY_LAYER (JSON options for the boundary layer, e.g.
// '{"entrainment":{"efficiency":0.3}}'), SURFACE (JSON options for the
// surface layer, e.g. '{"exchange":"fixed"}' for constant coefficients or
// '{"dragCoefficient":1.3e-3}', a fixed sea drag coefficient, which its
// heat and vapour exchange share, and the model top's sponge and Rayleigh
// drag), GRAVITY_WAVES (JSON options for the gravity-wave drag of
// js/physics/gravityWaves.module.js, e.g. '{"flux":3e-3}', or false),
// LAND (JSON options for the land surface, e.g. '{"treeMoisture":false}';
// '{"start":"bare"}' or '{"start":"green"}' make a fresh start's land
// bare or green in place of 'neutral', see START_CODES in
// js/physics/land.module.js), LAND_JUMPS (the
// model days at whose end land.jump() sets every land cell's trees and
// topsoil carbon to the equilibrium of its own record and ends a fresh
// start's hold, logging the change by land area globally and by 10° band,
// e.g. '365,730'; 'fresh', the default, jumps as the land's record passes
// 365 and 730 days when the land started fresh, which a state saved
// without a record never does; 'none' never),
// DIVERGENCE_DAMPING (the model's DIVERGENCE_DAMPING: the coefficient c of the
// core's divergence damping, the tendency c d²/dt ∇δ with d the mean
// distance between cell centres),
// RECORD (a directory to write each day's ocean and sea-ice forcing into as
// forcing-DDDD.bin, see js/forcing.module.js), BATCH (1: the steps are queued
// one at a time, waiting every eighth; more: that many steps go to the
// GPU in one submission, byte-identical, for drivers where submitting
// costs more than it does on Metal), LEVELS (bl36: the sigma grid of a
// fresh start, one of SIGMA_GRIDS in js/dynamics/sigmaCore.module.js; a
// run continues on its snapshot's grid), OCEAN_FROM (a state of the same
// N whose ocean, sea ice and sea surface replace the snapshot's, see
// js/oceanHandOff.module.js; the snapshots name it as oceanFrom, and a
// snapshot that already does is continued without replacing again),
// SYNC_CMD (a shell command run after every snapshot, forcing file and log
// update with the file's path as $1, see scripts/runControl.mjs),
// STOP_AFTER_STEPS (for tests: stop as on SIGTERM once this many steps
// have run), STRATOSPHERE (1: two more daily lines, the layer-mean
// temperature of every layer above 200 hPa, global and by zone, and those
// layers' winds, see scripts/upperAtmosphere.mjs), TOP_BUDGET (a number of
// steps a day, dividing the day's, at which the zonal-mean budgets of the top
// six layers' wind and temperature are sampled, read out daily as the 'top
// budget' line, see createTopBudget in scripts/upperAtmosphere.mjs; with
// BATCH 1 only).
// A fresh start can take from saved states: FROM, a state at
// the same N, gives the ocean, the land, the sea-surface temperature of
// its mixed layer and the land-surface temperature and, with
// ATMOSPHERE=carry (the default), its atmosphere (pi, theta, u, q, qc and
// the surface temperature, which fresh sea ice takes as its skin) and the
// mixed-layer deck's carried state, remapped onto the run's sigma grid
// when FROM is on another (remapLevels in js/physics/regrid.module.js);
// ATMOSPHERE=fresh starts the atmosphere and the deck from the initial
// state on the run's grid. The sea ice starts fresh and the clock at day
// 0; ICE_FROM=1 takes FROM's sea ice (thickness, concentration, snow, skin
// temperature) too; LAND_FROM, a state at any N, then gives the land.
// Without FROM, the ocean starts from the climatology CLIMATOLOGY (a file
// in js/ocean/climatology.module.js's format; data/woa_annual_1deg.bin
// when it exists, the World Ocean Atlas) and the sea surface from its
// mixed layer; CLIMATOLOGY=none starts it from the analytic climatology
// under the atmosphere's initial sea surface.
// scripts/spinup.sh runs segments back to back, scripts/asyncSpinup.sh
// alternates them with ocean-only spin-ups.
import { readFileSync, writeFileSync, readdirSync, renameSync, unlinkSync, appendFileSync, mkdirSync, existsSync } from 'node:fs';
import { basename } from 'node:path';
import { Grid } from '../js/grid.module.js';
import { topographyFromInt16, createGeography } from '../js/geography.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { createGpuModel } from '../js/gpu/model.gpu.js';
import { decodeState, encodeState, savedLevels } from '../js/stateFile.module.js';
import { sigmaInterfaces, sigmaGridName, R_DRY, CP_DRY, P0 } from '../js/dynamics/sigmaCore.module.js';
import { savedDeckField, DECK_FIELDS, savedMoistField, MOIST_FIELDS, savedRadiationField, RADIATION_FIELDS, regridLand, remapLevels } from '../js/physics/regrid.module.js';
import { readRanges } from '../js/gpu/device.module.js';
import { LAYER_DENSITIES, THERMOCLINE_DENSITY } from '../js/ocean/layered.module.js';
import { createForcingRecorder } from '../js/gpu/forcing.gpu.js';
import { forcingName } from '../js/forcing.module.js';
import { withOceanOf } from '../js/oceanHandOff.module.js';
import { CLIMATOLOGY_FILE } from '../js/ocean/climatology.module.js';
import { stopOnSignal, syncAfterSave } from './runControl.mjs';
import { freshJumpDue } from '../js/physics/land.module.js';
import { convectionLine, equatorLine } from '../js/audit.module.js';
import { upperWindLine, createTopBudget } from './upperAtmosphere.mjs';

const BOXES = {
  sahara: [16, 30, -10, 32], arabia: [16, 30, 38, 55], sahel: [8, 16, -15, 35], india: [15, 28, 72, 88], congo: [-5, 5, 12, 30], amazon: [-10, 3, -70, -50],
  seAsia: [10, 25, 95, 110], borneo: [-4, 7, 108, 119], europe: [45, 55, 0, 30], eastUS: [32, 45, -95, -75], siberia: [55, 65, 60, 120],
  ausInterior: [-30, -20, 120, 145], kalahari: [-27, -20, 17, 25], gobi: [38, 46, 90, 110], usSouthwest: [30, 37, -117, -106], cerrado: [-20, -10, -55, -42],
  ausNorth: [-18, -11, 125, 145], ausEast: [-37, -25, 148, 154], ausSoutheast: [-43, -34, 140, 150], ausWest: [-30, -20, 114, 120], newGuinea: [-11, -1, 130, 151],
};

const N = Number(process.env.N ?? 128), TAG = process.env.TAG ?? `spin${N}`, MINUTES = Number(process.env.MINUTES ?? 15), DAYS = Number(process.env.DAYS ?? Infinity), KEEP = Number(process.env.KEEP ?? 2);
const OUT = process.env.OUT ?? new URL('../runs/', import.meta.url).pathname;
const OCEAN = JSON.parse(process.env.OCEAN ?? '{}');
const RADIATION = { clearSkyPass: true, ...JSON.parse(process.env.RADIATION ?? '{}') }, MOIST = JSON.parse(process.env.MOIST ?? '{}'), BOUNDARY_LAYER = JSON.parse(process.env.BOUNDARY_LAYER ?? '{}'), SURFACE = JSON.parse(process.env.SURFACE ?? '{}'), LAND = JSON.parse(process.env.LAND ?? '{}'), DAMPING = process.env.DIVERGENCE_DAMPING === undefined ? {} : { divergenceDamping: Number(process.env.DIVERGENCE_DAMPING) };
const GRAVITY_WAVES = process.env.GRAVITY_WAVES === undefined ? {} : { gravityWaves: JSON.parse(process.env.GRAVITY_WAVES) };
const OCEAN_FROM = process.env.OCEAN_FROM, STOP_AFTER_STEPS = Number(process.env.STOP_AFTER_STEPS ?? Infinity);
const ATMOSPHERE = process.env.ATMOSPHERE ?? 'carry', STRATOSPHERE = process.env.STRATOSPHERE === '1';
const LAND_JUMPS = process.env.LAND_JUMPS ?? 'fresh', jumpDays = LAND_JUMPS === 'fresh' || LAND_JUMPS === 'none' ? [] : LAND_JUMPS.split(',').map(Number);
if (jumpDays.some((d) => !(Number.isInteger(d) && d > 0))) throw new Error(`LAND_JUMPS is 'fresh', 'none' or model days such as 365,730, not ${LAND_JUMPS}`);
if (ATMOSPHERE !== 'carry' && ATMOSPHERE !== 'fresh') throw new Error(`ATMOSPHERE is carry or fresh, not ${ATMOSPHERE}`);
const log = (line) => { console.log(line); appendFileSync(`${OUT}/${TAG}.log`, line + '\n'); };
const stop = stopOnSignal(log), hook = syncAfterSave(process.env.SYNC_CMD, log);
const positionOf = (file) => { const [, day, step] = file.match(/_day(\d+)(?:_step(\d+))?\.bin$/); return Number(day) * 1e6 + Number(step ?? 0); };
const inDay = (file) => /_step\d+\.bin$/.test(file);
const snapshots = () => readdirSync(OUT).filter((f) => f.startsWith(`${TAG}_day`) && /_day\d+(?:_step\d+)?\.bin$/.test(f)).sort((a, b) => positionOf(a) - positionOf(b));

const t0 = performance.now();
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const existing = snapshots(), file = existing[existing.length - 1];
let saved = file ? await decodeState(new Uint8Array(readFileSync(`${OUT}/${file}`))) : null;
if (saved && saved.N !== N) throw new Error(`${file} is N=${saved.N}`);
const levels = saved ? savedLevels(saved) : sigmaInterfaces(process.env.LEVELS ?? 'bl36');
const grid = `${sigmaGridName(levels) ?? 'a saved grid'} (${levels.length - 1} layers)`;
if (saved && process.env.LEVELS && sigmaGridName(levels) !== process.env.LEVELS) throw new Error(`${file} is on ${grid}, not ${process.env.LEVELS}`);
const fallback = new URL(`../${CLIMATOLOGY_FILE}`, import.meta.url).pathname, chosen = process.env.CLIMATOLOGY ?? (existsSync(fallback) ? fallback : 'none');
const climatology = !saved && !process.env.FROM && chosen !== 'none' ? chosen : null;
const model = await createGpuModel(new Grid(N), { topography, ocean: climatology ? { ...OCEAN, climatology } : OCEAN, radiation: RADIATION, moist: MOIST, boundaryLayer: BOUNDARY_LAYER, surface: SURFACE, land: LAND, ...DAMPING, ...GRAVITY_WAVES, levels });
const { mesh, core, state } = model;
const C = mesh.nCells, dt = 1350 * 16 / N, perDay = Math.round(86400 / dt), BATCH = Math.max(1, Math.round(Number(process.env.BATCH ?? 1)));
function loadSaved(saved) {
  ['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'].forEach((name, a) => state[a].set(saved[name]));
  model.seaIce.load(state[6], saved.concentration ?? null);
  for (const field of Object.keys(DECK_FIELDS)) model.radiation[field].set(savedDeckField(saved, field, model));
  for (const field of Object.keys(MOIST_FIELDS)) model.moist[field].set(savedMoistField(saved, field, model));
  for (const field of Object.keys(RADIATION_FIELDS)) model.radiation[field].set(savedRadiationField(saved, field, model));
  if (saved.boundaryDepth) model.boundaryLayer.depth.set(saved.boundaryDepth);
  if (saved.mixingTop) model.boundaryLayer.mixingTop.set(saved.mixingTop);
  if (saved.boundaryRegime) model.boundaryLayer.regime.set(saved.boundaryRegime);
  if (saved.boundaryBuoyancy) model.boundaryLayer.buoyancyFlux.set(saved.boundaryBuoyancy);
  if (saved.windSpeed) model.surface.windSpeed.set(saved.windSpeed);
  if (saved.evaporation) model.radiation.evaporation.set(saved.evaporation);
  if (saved.exchangeHeat && model.exchange && !model.exchange.fixed) model.exchange.heat.set(saved.exchangeHeat);
  if (saved.exchangeWind && model.exchange) model.exchange.wind.set(saved.exchangeWind);
  for (const field of ['cumulusCover', 'cumulusWater']) if (saved[field] && saved[field].length === model.moist[field].length) model.moist[field].set(saved[field]);
  if (saved.subcloudVirtual && saved.subcloudVirtual.length === model.moist.subcloudVirtual.length) model.moist.subcloudVirtual.set(saved.subcloudVirtual); else model.moist.subcloudVirtual.fill(0);
  model.time = saved.time;
  model.load();
  model.ocean.load(saved.ocean, state[3], state[6]);
  model.land.load(saved.land);
}
let startStep = 0, oceanYears = 0, oceanFrom = null;
if (saved) {
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
  log(`--- ${new Date().toISOString()} continuing from ${file} (day ${saved.day}${startStep ? `, step ${startStep} of ${perDay}` : ''}) on ${grid} after ${((performance.now() - t0) / 1000).toFixed(0)} s of setup`);
} else {
  if (OCEAN_FROM) throw new Error(`OCEAN_FROM needs a ${TAG}_dayDDDD.bin in ${OUT} to hand the ocean to`);
  initializeState(model, { geostrophic: !model.surfaceGeopotential }).forEach((values, a) => state[a].set(values));
  for (let i = 0; i < C; i++) if (model.geography.land[i]) state[6][i] = 0;
  const from = process.env.FROM ? await decodeState(new Uint8Array(readFileSync(process.env.FROM))) : null, iceFrom = !!from && process.env.ICE_FROM === '1';
  const fromLevels = from ? savedLevels(from) : null, fromGrid = from ? sigmaGridName(fromLevels) ?? 'a saved grid' : null;
  let atmosphere = `a fresh atmosphere on ${grid}`;
  if (from) {
    if (from.N !== N) throw new Error(`FROM ${process.env.FROM} is N=${from.N}, not ${N}`);
    if (ATMOSPHERE === 'carry') {
      const { pi, theta, u, q, qc } = remapLevels(fromLevels, levels, from, mesh);
      [pi, theta, u, from.surfaceT, q, qc].forEach((values, a) => { if (values) state[a].set(values); });
      for (const field of Object.keys(DECK_FIELDS)) model.radiation[field].set(savedDeckField(from, field, model));
      for (const field of Object.keys(MOIST_FIELDS)) model.moist[field].set(savedMoistField(from, field, model));
      for (const field of Object.keys(RADIATION_FIELDS)) model.radiation[field].set(savedRadiationField(from, field, model));
      const same = fromLevels.length === levels.length && fromLevels.every((sigma, k) => sigma === levels[k]);
      atmosphere = `its atmosphere${same ? '' : ` remapped from ${fromGrid}`} on ${grid} and its deck`;
    }
    const land = model.geography.land;
    for (let i = 0; i < C; i++) {
      if (land[i]) { state[3][i] = from.surfaceT[i]; continue; }
      if (iceFrom) state[6][i] = from.ice[i];
      if (iceFrom && state[6][i] > 0) state[3][i] = from.surfaceT[i];
      else if (!(state[6][i] > 0) && from.ocean.h[i] > 1) state[3][i] = from.ocean.T[i];
    }
    if (iceFrom) model.seaIce.load(state[6], from.concentration ?? null);
  }
  model.load();
  if (from) {
    model.ocean.load(from.ocean, state[3], state[6]);
    model.land.load(from.land, iceFrom ? state[6] : null);
    log(`seeded from ${process.env.FROM} (N=${from.N}, day ${from.day}, ${fromGrid}): the ocean, the land (soil, snow, ${from.land.snowAlbedo ? 'snow albedo, ' : ''}${from.land.vegetation ? 'vegetation, ' : ''}${from.land.canopy ? 'standing cover, ' : ''}${from.land.seasonLength ? 'season means, ' : ''}${from.land.surface ? 'surface water, ' : ''}surface temperature) and the sea-surface temperature of its mixed layer, ${iceFrom ? 'and its sea ice (thickness, concentration, snow, skin temperature)' : 'with fresh sea ice'}; ${atmosphere}; the clock at day 0`);
  } else {
    const started = model.ocean.initialize(state[3], state[6]);
    model.land.initialize();
    log(started ? `ocean and sea surface from ${climatology}: ${started.atlas} sea cells, ${started.analytic} from the analytic climatology where it has no water` : 'ocean from the analytic climatology');
  }
  if (process.env.LAND_FROM) {
    const seed = await decodeState(new Uint8Array(readFileSync(process.env.LAND_FROM)));
    const sourceMesh = seed.N === N ? mesh : buildMesh(new Grid(seed.N));
    const source = { mesh: sourceMesh, geography: seed.N === N ? model.geography : createGeography(sourceMesh, topography) };
    model.land.load(regridLand(source, { mesh, geography: model.geography, land: model.land }, seed.land, null, { ice: seed.ice, surfaceT: seed.surfaceT }), state[6]);
    log(`land seeded from ${process.env.LAND_FROM} (N=${seed.N}, day ${seed.day})`);
  }
  const { record } = model.land, held = ['', 'the trees at half the cover and the topsoil carbon at', 'no trees and the topsoil carbon at', 'the trees at the cover times their record\'s factor and the topsoil carbon at'][record[1]];
  log(`land: record ${record[0] < 0 ? 'taken over (exponential means)' : `${(record[0] / 86400).toFixed(0)} days old`}${held ? `, ${['', 'neutral', 'bare', 'green'][record[2]]} start holding ${held} ${Object.values(model.land.placeholders)[record[1] - 1].carbon.toFixed(3)} kg/m²` : ''}; jumps ${LAND_JUMPS === 'fresh' ? (record[2] >= 1 ? 'as the record passes 365 and 730 days' : 'none (the land did not start fresh)') : LAND_JUMPS}`);
  log(`--- ${new Date().toISOString()} fresh start at N=${N} on ${grid} (${C} cells, dt ${dt} s, ${perDay} steps a day) after ${((performance.now() - t0) / 1000).toFixed(0)} s of setup`);
}
const everySteps = model.oceanEngine ? model.oceanEngine.everySteps : 1;
if (startStep % everySteps) throw new Error(`the snapshot stopped at step ${startStep}, off the ocean's ${everySteps}-step cadence`);

const deg = 180 / Math.PI, land = model.geography.land;
const inBox = Object.fromEntries(Object.entries(BOXES).map(([k, [a, b, c, e]]) => [k, [...Array(C).keys()].filter((i) => land[i] && mesh.latCell[i] * deg >= a && mesh.latCell[i] * deg <= b && mesh.lonCell[i] * deg >= c && mesh.lonCell[i] * deg <= e)]));
const PH = model.gpu.layout.PH;
const split = { convective: new Float64Array(C), largeScale: new Float64Array(C), wet: new Float64Array(C), seconds: 0 };
let splitTime = model.time;
const readSplit = async () => {
  const [convective, largeScale] = await readRanges(model.gpu.device, model.gpu.buffers.PH, [{ offset: PH.CONVMEAN, length: C }, { offset: PH.CONDMEAN, length: C }]);
  const seconds = model.time - splitTime, days = seconds / 86400;
  splitTime = model.time;
  if (!(seconds > 0)) return;
  for (let i = 0; i < C; i++) { split.convective[i] += days * convective[i]; split.largeScale[i] += days * largeScale[i]; if (convective[i] > 0) split.wet[i] += days; }
  split.seconds += seconds;
};
const iceArea = async () => { const { fields } = await model.beginFrame({ fields: ['concentration', 'ice'] }); let north = 0, south = 0; for (let i = 0; i < C; i++) if (fields.concentration[i] >= 0.15) { if (mesh.latCell[i] > 0) north += mesh.areaCell[i]; else south += mesh.areaCell[i]; } return [north / 1e12, south / 1e12, fields.ice]; };
await model.diagnostics();
if (startStep && saved.rainTotal) {
  const queue = model.gpu.device.queue;
  queue.writeBuffer(model.gpu.buffers.PH, 4 * PH.RAIN, Float32Array.from(saved.rainTotal));
  queue.writeBuffer(model.gpu.buffers.PH, 4 * PH.RUNOFF, Float32Array.from(saved.runoffTotal));
  if (model.oceanEngine && saved.rainSeen) {
    const { OD } = model.oceanEngine.layout;
    queue.writeBuffer(model.oceanEngine.buffers.OD, 4 * OD.RAINSEEN, Float32Array.from(saved.rainSeen));
    queue.writeBuffer(model.oceanEngine.buffers.OD, 4 * OD.RUNOFFSEEN, Float32Array.from(saved.runoffSeen));
  }
}
const RECORD = process.env.RECORD;
if (RECORD) mkdirSync(RECORD, { recursive: true });
const partOfDay = RECORD && startStep && saved.forcingSums ? { sums: saved.forcingSums, steps: saved.forcingSteps, oceanSteps: saved.forcingOceanSteps, seconds: saved.forcingSeconds, rain: saved.forcingRain, runoff: saved.forcingRunoff } : null;
if (RECORD && startStep && !partOfDay) log(`the snapshot holds no recorded part of its day, so forcing day ${saved.day + 1} covers its last ${perDay - startStep} steps alone`);
const recorder = RECORD ? await createForcingRecorder(model, partOfDay) : null;
if (stop.requested) { log(`stopped by ${stop.requested} before the first step; nothing to save`); await hook.drain(); process.exit(0); }
const start = performance.now();
const day0 = Math.round((model.time - startStep * dt) / 86400);
let day = day0, step = startStep, dayStart = startStep, taken = 0, iceNorth = 0, iceSouth = 0;
log(`ASR, atmosphere, OLR and albedo below are day means over the day's ${perDay} steps${startStep ? ` (day ${day0 + 1}'s over its last ${perDay - startStep})` : ''}, the albedo the day's reflected over its incoming sunlight${RADIATION.clearSkyPass ? ', and so are SWCRE and LWCRE, the shortwave (ASR less clear-sky ASR) and longwave (clear-sky OLR less OLR) cloud effects' : ''}`);
const halted = () => stop.requested || taken >= STOP_AFTER_STEPS;
const TOP_BUDGET = Number(process.env.TOP_BUDGET ?? 0), budget = TOP_BUDGET ? createTopBudget(mesh, levels, { surface: SURFACE }) : null;
let budgetArmed = false;
if (budget) {
  if (BATCH !== 1 || perDay % TOP_BUDGET) throw new Error(`TOP_BUDGET needs BATCH 1 and a number of samples dividing the day's ${perDay} steps`);
  const S = model.gpu.layout.S, top = 6 * C, ocean = model.gpu.hooks.beforePhysics;
  let before = null;
  model.gpu.hooks.beforePhysics = async (dt, n) => {
    if (ocean) await ocean(dt, n);
    if (budgetArmed) [before] = await readRanges(model.gpu.device, model.gpu.buffers.S, [{ offset: S.TH, length: top }]);
  };
  model.gpu.hooks.afterPhysics = async (dt) => {
    if (!budgetArmed) return;
    const [after, pi] = await readRanges(model.gpu.device, model.gpu.buffers.S, [{ offset: S.TH, length: top }, { offset: S.PI, length: C }]);
    const [east] = await readRanges(model.gpu.device, model.gpu.buffers.PH, [{ offset: PH.GWE, length: top }]);
    budget.radiation(before, after, pi, dt);
    budget.waves(east);
  };
  await model.sync();
  budget.start(state, model.time);
}
const ZONES = [['global', -90, 90], ['20S-20N', -20, 20], ['35-55N', 35, 55], ['35-55S', -55, -35], ['70-90N', 70, 90], ['70-90S', -90, -70]];
function stratosphereLine(day) {
  const kappa = R_DRY / CP_DRY, K = levels.length - 1, rows = [];
  for (let k = 0; k < K && 0.5 * (levels[k] + levels[k + 1]) < 0.2; k++) {
    const shape = (levels[k + 1] ** (1 + kappa) - levels[k] ** (1 + kappa)) / ((1 + kappa) * (levels[k + 1] - levels[k]));
    const means = ZONES.map(([, a, b]) => {
      let sum = 0, area = 0;
      for (let i = 0; i < C; i++) { const lat = mesh.latCell[i] * deg; if (lat < a || lat > b) continue; sum += mesh.areaCell[i] * state[1][k * C + i] * (state[0][i] / P0) ** kappa * shape; area += mesh.areaCell[i]; }
      return (sum / area).toFixed(1);
    });
    rows.push(`${(500 * (levels[k] + levels[k + 1])).toPrecision(3)} hPa ${means.join('/')}`);
  }
  return `stratosphere day ${day}, layer-mean temperature (K; ${ZONES.map(([name]) => name).join('/')}), layers at their pressure under 1000 hPa: ${rows.join('; ')}`;
}
const band = (b) => `${Math.abs(b.from)}${b.from > 0 ? 'N' : b.from < 0 ? 'S' : ''}–${Math.abs(b.to)}${b.to > 0 ? 'N' : b.to < 0 ? 'S' : ''}`;
const change = ([before, after], digits) => `${before.toFixed(digits)} → ${after.toFixed(digits)}`;
const jumpLine = (day, r) => `land jump at the end of day ${day} (record ${r.age < 0 ? 'taken over (exponential means)' : `${(r.age / 86400).toFixed(2)} days`}): land means trees ${change(r.global.trees, 3)}, topsoil carbon ${change(r.global.carbon, 2)} kg/m², dry soil albedo ${change(r.global.dry, 3)}, land albedo ${change(r.global.albedo, 3)}; by band (share of the land: trees, carbon, land albedo) ${r.bands.map((b) => `${band(b)} (${b.share.toFixed(3)}: ${change(b.trees, 2)}, ${change(b.carbon, 1)}, ${change(b.albedo, 3)})`).join(', ')}`;
let landAge = model.land.record[0];
for (;;) {
  if (BATCH === 1) {
    while (step < perDay) {
      budgetArmed = !!budget && (step + 1) % (perDay / TOP_BUDGET) === 0;
      await model.step(dt); recorder?.step(); step++; taken++;
      if (budgetArmed) { await model.sync(); budget.sample(state); budgetArmed = false; }
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
  const [absorbedSum, atmosphereSum, longwaveSum] = await readRanges(model.gpu.device, model.gpu.buffers.PH, [{ offset: PH.ABSSUM, length: C }, { offset: PH.ATMSUM, length: C }, { offset: PH.LWSFCSUM, length: C }]);
  let seaSolar = 0, seaLongwave = 0, seaArea = 0;
  for (let i = 0; i < C; i++) if (!land[i]) { seaSolar += mesh.areaCell[i] * (absorbedSum[i] - atmosphereSum[i]); seaLongwave += mesh.areaCell[i] * longwaveSum[i]; seaArea += mesh.areaCell[i]; }
  const stepsToday = perDay - dayStart;
  seaSolar /= seaArea * stepsToday;
  seaLongwave /= seaArea * stepsToday;
  step = 0; dayStart = 0;
  day++;
  const d = await model.diagnostics();
  await readSplit();
  if (recorder) {
    const file = `${RECORD}/${forcingName(day)}`;
    writeFileSync(`${file}.partial`, await recorder.day(day));
    renameSync(`${file}.partial`, file);
    hook.after(file);
  }
  const [north, south, iceNow] = await iceArea();
  const iceSolar = [0, 0], icedArea = [0, 0];
  for (let i = 0; i < C; i++) if (!land[i] && iceNow[i] > 0 && Math.abs(mesh.latCell[i] * deg) >= 60) { const s = mesh.latCell[i] > 0 ? 0 : 1; iceSolar[s] += mesh.areaCell[i] * (absorbedSum[i] - atmosphereSum[i]); icedArea[s] += mesh.areaCell[i]; }
  const iceSurface = iceSolar.map((x, s) => (icedArea[s] > 0 ? x / (icedArea[s] * stepsToday) : 0));
  iceNorth += north; iceSouth += south;
  const minutes = (performance.now() - start) / 60000;
  log(`day ${day} (${minutes.toFixed(1)} min): Ts ${(d.meanSurfaceT - 273.15).toFixed(2)} °C, ASR ${d.absorbedSolar.toFixed(1)} (atmosphere ${d.atmosphereSolar.toFixed(1)}) OLR ${d.outgoingLongwave.toFixed(1)} W/m², ps ${(d.piMin / 100).toFixed(0)}–${(d.piMax / 100).toFixed(0)} hPa, max wind ${d.maxWind.toFixed(1)} m/s, precip ${(86400 * d.precipitation).toFixed(2)} mm/d, ice ${(100 * d.iceFraction).toFixed(1)}% (N ${north.toFixed(1)} S ${south.toFixed(1)} Mkm²), albedo ${d.planetaryAlbedo.toFixed(3)}, ${d.shortwaveCloudEffect === undefined ? '' : `SWCRE ${d.shortwaveCloudEffect.toFixed(1)} LWCRE ${d.longwaveCloudEffect.toFixed(1)}, clear-sky reflectance ${(1 - d.clearAbsorbedSolar * (1 - d.planetaryAlbedo) / d.absorbedSolar).toFixed(4)}, `}sea surface shortwave ${seaSolar.toFixed(1)} net longwave ${seaLongwave.toFixed(1)} W/m² (iced cells poleward of 60°: N ${iceSurface[0].toFixed(1)} S ${iceSurface[1].toFixed(1)}), ocean h1 ${d.oceanUpperDepth.toFixed(0)} m, interior ${(d.oceanInteriorT - 273.15).toFixed(2)} °C, currents ≤ ${d.oceanSpeed.toFixed(2)} m/s, transport ${d.oceanTransport.toFixed(0)} Sv, clamped ${d.oceanLimited}`);
  if (STRATOSPHERE) { await model.sync(); log(stratosphereLine(day)); log(upperWindLine(mesh, levels, state, dt, day)); }
  if (budget) { await model.sync(); log(budget.line(day, state, model.time)); }
  if (!Number.isFinite(d.meanSurfaceT) || !Number.isFinite(d.maxWind) || !Number.isFinite(d.oceanSpeed)) { log(`NaN on day ${day}; stopping`); await hook.drain(); process.exit(2); }
  if (jumpDays.includes(day) || (LAND_JUMPS === 'fresh' && freshJumpDue(model.land.record, landAge, model.land.record[0]))) log(jumpLine(day, await model.land.jump()));
  landAge = model.land.record[0];
  if (minutes >= MINUTES || day >= DAYS || halted()) break;
}

const partial = {};
if (step) {
  const [rainTotal, runoffTotal] = await readRanges(model.gpu.device, model.gpu.buffers.PH, [{ offset: PH.RAIN, length: C }, { offset: PH.RUNOFF, length: C }]);
  Object.assign(partial, { rainTotal, runoffTotal });
  if (model.oceanEngine) {
    const { OD } = model.oceanEngine.layout, [rainSeen, runoffSeen] = await readRanges(model.gpu.device, model.oceanEngine.buffers.OD, [{ offset: OD.RAINSEEN, length: C }, { offset: OD.RUNOFFSEEN, length: C }]);
    Object.assign(partial, { rainSeen, runoffSeen });
  }
}
if (step && recorder) {
  await model.diagnostics();
  const part = await recorder.checkpoint();
  Object.assign(partial, { forcingSteps: part.steps, forcingOceanSteps: part.oceanSteps, forcingSeconds: part.seconds, forcingSums: part.sums, forcingRain: part.rain, forcingRunoff: part.runoff });
}
await model.sync();
const days = day - day0;
const ocean = await model.ocean.serialize({ restart: true }), landState = await model.land.serialize();
const kT = LAYER_DENSITIES.filter((r) => r < THERMOCLINE_DENSITY).length;
const sea = (a, b, c, e) => [...Array(C).keys()].filter((i) => !land[i] && mesh.latCell[i] * deg >= a && mesh.latCell[i] * deg <= b && (c <= e ? mesh.lonCell[i] * deg >= c && mesh.lonCell[i] * deg <= e : mesh.lonCell[i] * deg >= c || mesh.lonCell[i] * deg <= e));
const meanOf = (cells, f) => cells.reduce((s, i) => s + f(i), 0) / Math.max(1, cells.length);
const sst = (cells) => meanOf(cells, (i) => ocean.T[i] - 273.15), classTop = (cells) => meanOf(cells, (i) => { let d = 0; for (let k = 0; k <= kT; k++) d += ocean.h[k * C + i]; return d; });
const warmPool = sea(-10, 10, 120, 160), coldTongue = sea(-2, 2, -110, -90), westPacific = sea(-5, 5, 140, 170), eastPacific = sea(-5, 5, -120, -90);
if (days > 0 && !step) {
  const splitDays = split.seconds / 86400, rainOf = (i) => (split.convective[i] + split.largeScale[i]) / splitDays;
  log(`regions after ${days} days (rain mm/d / vegetation / surface °C): ` + Object.entries(inBox).map(([name, cells]) => {
    let r = 0, v = 0, t = 0;
    for (const i of cells) { r += rainOf(i); v += landState.vegetation ? landState.vegetation[i] : 0; t += state[3][i]; }
    const n = Math.max(1, cells.length);
    return `${name} ${(r / n).toFixed(1)}/${(v / n).toFixed(2)}/${(t / n - 273.15).toFixed(0)}`;
  }).join(', ') + `; sea ice mean N ${(iceNorth / days).toFixed(1)} S ${(iceSouth / days).toFixed(1)} Mkm²`);
  log(`ocean after ${days} days: warm pool ${sst(warmPool).toFixed(1)} °C, cold tongue ${sst(coldTongue).toFixed(1)} °C (W−E ${(sst(westPacific) - sst(eastPacific)).toFixed(1)} K), ${THERMOCLINE_DENSITY} class top W Pac ${classTop(westPacific).toFixed(0)} m, E Pac ${classTop(eastPacific).toFixed(0)} m`);
  const daily = (sums) => Float64Array.from(sums, (x) => x / splitDays);
  log(convectionLine(mesh, { land, ice: state[6] }, { convective: daily(split.convective), largeScale: daily(split.largeScale), wet: daily(split.wet) }, days));
  const [stress] = await readRanges(model.gpu.device, model.oceanEngine.buffers.OD, [{ offset: model.oceanEngine.layout.OD.STRESS, length: mesh.nEdges }]);
  log(equatorLine(mesh, land, ocean, stress, days));
}
const name = `${TAG}_day${String(day).padStart(4, '0')}${step ? `_step${String(step).padStart(4, '0')}` : ''}.bin`;
const [pi, theta, u, surfaceT, q, qc, ice] = state, { concentration } = model.seaIce, { mlmSubsidence, mlmHeight, mlmGate, meanAbsorbedSolar, meanOutgoingLongwave, meanPlanetaryAlbedo, meanShortwaveCloudEffect, meanLongwaveCloudEffect, evaporation } = model.radiation, { convectiveRain, largeScaleRain, cumulusCover, cumulusWater, subcloudVirtual } = model.moist, { windSpeed } = model.surface, boundaryDepth = model.boundaryLayer.depth, mixingTop = model.boundaryLayer.mixingTop, boundaryRegime = model.boundaryLayer.regime, boundaryBuoyancy = model.boundaryLayer.buoyancyFlux;
const exchangeFields = model.exchange ? { ...(model.exchange.fixed ? {} : { exchangeHeat: model.exchange.heat }), exchangeWind: model.exchange.wind } : {};
const header = { N, K: core.K, day, ...(step ? { step } : {}), time: model.time, terrain: !!model.surfaceGeopotential, levels: core.levels, ...(oceanYears ? { oceanYears } : {}), ...(oceanFrom ? { oceanFrom } : {}) };
writeFileSync(`${OUT}/${name}.partial`, encodeState({ ...header, pi, theta, u, surfaceT, q, qc, ice, concentration, mlmSubsidence, mlmHeight, mlmGate, convectiveRain, largeScaleRain, meanAbsorbedSolar, meanOutgoingLongwave, meanPlanetaryAlbedo, ...(RADIATION.clearSkyPass ? { meanShortwaveCloudEffect, meanLongwaveCloudEffect } : {}), boundaryDepth, mixingTop, boundaryRegime, boundaryBuoyancy, windSpeed, evaporation, ...exchangeFields, cumulusCover, cumulusWater, subcloudVirtual, ocean, land: landState, ...partial }, { f64: ['forcingRain', 'forcingRunoff'] }));
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
