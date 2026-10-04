import { Grid } from './grid.module.js';
import { createModel } from './model.module.js';
import { createParallelModel } from './parallel.module.js';
import { createGpuModel } from './gpu/model.gpu.js';
import { withCadence } from './cadence.module.js';
import { initializeState } from './physics/init.module.js';
import { cellVector } from './dynamics/operators.module.js';
import { regridState, regridOcean, regridLand, regridConcentration, savedDeckField, DECK_FIELDS, savedMoistField, MOIST_FIELDS, savedRadiationField, RADIATION_FIELDS } from './physics/regrid.module.js';
import { topographyFromInt16, rebalanceSurfacePressure, surfaceGeopotential, decodeSubgrid, subgridUrl } from './geography.module.js';
import { decodeClimatology, CLIMATOLOGY_FILE } from './ocean/climatology.module.js';
import { regridCellField } from './physics/regrid.module.js';
import { levelFields, dewPoint, wetBulb, miseryIndex, verticalVelocity, smoothCells } from './levels.module.js';
import { initialHumidity } from './physics/init.module.js';
import { fetchState, stateName, savedLevels } from './stateFile.module.js';
import { sigmaInterfaces } from './dynamics/sigmaCore.module.js';
import { LEVEL_FIELDS, OCEAN_FIELDS, CLOUD_TYPES, RAIN_MEMORY, VERTICAL_MEMORY } from './frames.module.js';
import { createPacer } from './pace.module.js';
import { profileGpu } from './gpu/profile.module.js';
import { freshJumpDue } from './physics/land.module.js';

const FREEZING = 273.15;
let model = null, serving = false, running = false, dt = 450, stepsPerFrame = 24, frame = 0;
let subscription = { level: 'surface', depth: 'surface', fields: [], diagnostics: false }, layerWinds = [];
const rain = { total: null, time: 0 };
const vertical = { values: null, level: null, time: 0 };

function restartRain() { rain.total = null; rain.time = model.time; if (model.restartPrecipitation) model.restartPrecipitation(); }

/*
 * A frame for the CPU engines, from the model's own arrays: the rain
 * folds into its three-hour memory every frame, and only the subscribed
 * fields and diagnostics are built. The cell-center winds of a layer are
 * reconstructed only when that layer bounds the level somewhere.
 */
async function cpuFrame({ level, depth, fields, diagnostics: summarize }) {
  const { mesh, core, state } = model;
  const C = mesh.nCells, E = mesh.nEdges, [pi, theta, u, , q, , iceField] = state;
  const time = model.time, interval = time - rain.time;
  if (!rain.total || interval < 0) rain.total = new Float32Array(C);
  const keep = Math.exp(-Math.max(0, interval) / RAIN_MEMORY), fallen = model.moist.precipitation;
  for (let i = 0; i < C; i++) rain.total[i] = rain.total[i] * keep + fallen[i];
  rain.time = time;
  let diagnostics = null;
  if (summarize) {
    diagnostics = await model.diagnostics();
    if (interval <= 0) {
      let area = 0, recent = 0;
      for (let i = 0; i < C; i++) { area += mesh.areaCell[i]; recent += mesh.areaCell[i] * rain.total[i]; }
      diagnostics.precipitation = recent / area / RAIN_MEMORY;
    }
  } else {
    model.restartPrecipitation();
  }
  const want = new Set(fields), out = {};
  const onLand = model.geography ? model.geography.land : null;
  const sea = (source) => Float32Array.from({ length: C }, (_, i) => (onLand && onLand[i] ? NaN : source[i]));
  if (fields.some((name) => LEVEL_FIELDS.has(name))) {
    layerWinds.fill(null);
    const layerWind = (k) => layerWinds[k] ??= cellVector(mesh, u.subarray(k * E, (k + 1) * E), new Float64Array(3 * C));
    const f = levelFields(core, pi, theta, layerWind, level, q);
    for (const [name, values] of Object.entries({ temperature: f.temperature, height: f.height, humidity: f.humidity, speed: f.speed, wind: f.vector })) if (want.has(name)) out[name] = values;
    const comfort = (fn) => Float32Array.from(f.temperature, (t, i) => fn(t - FREEZING, f.humidity[i], f.speed[i]) + FREEZING);
    if (want.has('dewPoint')) out.dewPoint = comfort(dewPoint);
    if (want.has('wetBulb')) out.wetBulb = comfort(wetBulb);
    if (want.has('misery')) out.misery = comfort(miseryIndex);
    if (want.has('vertical')) {
      const smoothed = smoothCells(mesh, verticalVelocity(mesh, core, pi, u, level, f.temperature));
      const keep = vertical.values && vertical.level === level ? Math.exp(-Math.max(0, time - vertical.time) / VERTICAL_MEMORY) : 0;
      if (!vertical.values || vertical.values.length !== C) vertical.values = new Float32Array(C);
      for (let i = 0; i < C; i++) vertical.values[i] = keep * vertical.values[i] + (1 - keep) * smoothed[i];
      vertical.level = level; vertical.time = time;
      out.vertical = Float32Array.from(vertical.values);
    } else vertical.level = null;
  }
  if (want.has('ps')) out.ps = Float32Array.from(pi);
  if (want.has('mslp')) {
    out.mslp = Float32Array.from(pi);
    if (model.surfaceGeopotential) {
      const phis = model.surfaceGeopotential, { R, g, exnerLayer } = core.diagnostics, bottom = (core.K - 1) * C;
      for (let i = 0; i < C; i++) out.mslp[i] = pi[i] * Math.exp(phis[i] / (R * (theta[bottom + i] * exnerLayer[bottom + i] + 0.00325 * phis[i] / g)));
    }
  }
  if (want.has('water')) out.water = Float32Array.from({ length: C }, (_, i) => model.moist.columnWater(pi, q, i));
  if (want.has('cloud')) out.cloud = Float32Array.from({ length: C }, (_, i) => model.cloudWater(i));
  if (CLOUD_TYPES.some((name) => want.has(name))) {
    const parts = Array.from({ length: C }, (_, i) => model.cloudParts(i));
    for (const name of CLOUD_TYPES) if (want.has(name)) out[name] = Float32Array.from(parts, (part) => part[name]);
  }
  if (want.has('rain')) out.rain = Float32Array.from(rain.total);
  if (want.has('ice')) out.ice = Float32Array.from(iceField);
  if (want.has('concentration')) out.concentration = Float32Array.from(model.seaIce.concentration);
  if (want.has('albedo')) out.albedo = Float32Array.from(iceField, (h, i) => (onLand && onLand[i] ? model.land.albedo(i) : model.seaIce.albedo(h, null, model.seaIce.snow[i], model.seaIce.cover(i, h), model.state[3][i], model.seaIce.snowAlbedo[i])));
  if (want.has('shortwave')) out.shortwave = Float32Array.from(model.radiation.surfaceShortwave);
  if (want.has('longwave')) out.longwave = Float32Array.from(model.radiation.outgoing);
  if (model.land && want.has('soil')) out.soil = Float32Array.from(model.land.soil);
  if (model.land && want.has('snow')) out.snow = Float32Array.from(model.land.snow);
  if (model.land && want.has('vegetation')) out.vegetation = Float32Array.from(model.land.vegetation);
  const ocean = model.oceanFields && fields.some((name) => OCEAN_FIELDS.has(name)) ? model.oceanFields(depth === 'surface' ? 0 : Number(depth)) : null;
  if (ocean) {
    for (const [name, values] of Object.entries({ sst: ocean.temperature, sss: ocean.S1, layerDepth: ocean.h1, thermocline: ocean.thermoclineDepth, ssh: ocean.eta, current: ocean.speed, upwelling: ocean.upwelling })) if (want.has(name)) out[name] = sea(values);
    if (want.has('currents')) {
      const vector = ocean.current;
      if (onLand) for (let i = 0; i < C; i++) if (onLand[i]) vector.fill(0, 3 * i, 3 * i + 3);
      out.currents = Float32Array.from(vector);
    }
  }
  return { time, level, depth, fields: out, diagnostics };
}

const captureFrame = () => (model.beginFrame ? model.beginFrame(subscription) : cpuFrame(subscription));

/*
 * The energy record takes a day from the frames' diagnostics: each frame's
 * means, over the steps since the diagnostics before it, are summed into
 * the day its last step falls in, and the day goes in once its sum holds a
 * whole day's steps to within a frame (the frames do not end on the day).
 * A frame whose means span more than a frame (the panel was closed, so no
 * diagnostics were taken) restarts the sum, so the record gains days only
 * while the panel stays open; the state's own record stands otherwise.
 */
let energyDay = { day: 0, steps: 0, asr: 0, olr: 0, ts: 0 };
function feedEnergy(d, time) {
  if (!(d.meanSteps > 0)) return;
  const perDay = Math.round(86400 / dt), day = Math.ceil(time / 86400 - 1e-9);
  if (d.meanSteps > stepsPerFrame || day !== energyDay.day) {
    if (day === energyDay.day + 1 && Math.abs(energyDay.steps - perDay) <= stepsPerFrame && d.meanSteps <= stepsPerFrame) model.energyRecord.add(energyDay.day, { asr: energyDay.asr / energyDay.steps, olr: energyDay.olr / energyDay.steps, ts: energyDay.ts });
    energyDay = { day, steps: 0, asr: 0, olr: 0, ts: 0 };
    if (d.meanSteps > stepsPerFrame) return;
  }
  energyDay.steps += d.meanSteps; energyDay.asr += d.meanSteps * d.absorbedSolar; energyDay.olr += d.meanSteps * d.outgoingLongwave; energyDay.ts = d.meanSurfaceT;
}
function placeEnergy(model, saved) {
  if (!model.energyRecord) return;
  model.energyRecord.clear();
  if (saved && saved.energyRecord && saved.energyRecord.length === model.energyRecord.values.length) model.energyRecord.load(saved.energyRecord);
  energyDay = { day: 0, steps: 0, asr: 0, olr: 0, ts: 0 };
}

function postFrame({ time, level, depth, fields, diagnostics }) {
  if (diagnostics && model.energyRecord) { feedEnergy(diagnostics, time); diagnostics.energy = model.energyRecord.windows(); diagnostics.energyDay = model.energyRecord.newest; }
  const transfer = [...new Set(Object.values(fields).map((values) => values.buffer))];
  self.postMessage({ type: 'frame', frame: frame++, time, day: time / 86400, level, depth, engine: model.engine ?? 'cpu', pause: model.beginFrame ? pace.pause : 0, fields, diagnostics }, transfer);
}

async function sendFrame() { postFrame(await captureFrame()); }

/*
 * A frame outside the loop, reporting any failure on the status line.
 * While a batch of steps is still finishing after a pause, it waits for
 * that batch's own frame, so that the page's last frame carries the
 * latest subscription. halt() pauses and resolves once no batch is
 * stepping.
 */
let stepping = false, resend = false, batchDone = Promise.resolve();
async function halt() { running = false; resend = false; await batchDone; }
const report = (error) => status(`error: ${error && error.stack ? error.stack : error}`);
function refresh() {
  if (stepping) { resend = true; return; }
  sendFrame().catch(report);
}

/*
 * Work on the model outside the loop (the start, a snapshot, a restore, a
 * profile, the device test) holds it, one holder at a time: the loop
 * stops first, pause and resume requests that arrive meanwhile are kept
 * in `held.resume` and take effect when the holder is done, and a frame
 * asked for meanwhile follows then. `resume` is whether the loop runs
 * afterwards unless told otherwise. A restore that comes before the first
 * start waits in `pendingRestore` and follows it.
 */
let held = null, pendingRestore = null;
async function hold(work, { resume = false } = {}) {
  while (held) await held.done;
  let release;
  held = { resume, done: new Promise((resolve) => { release = resolve; }) };
  let finish;
  try {
    await halt();
    batchDone = new Promise((resolve) => { finish = resolve; });
    stepping = true;
    await work();
  } finally {
    stepping = false;
    if (finish) finish();
    const again = held.resume;
    held = null;
    release();
    if (again && model) { running = true; loop(); } else if (resend) { resend = false; refresh(); }
  }
}

/*
 * The GPU draws the page's globe too, and it takes queued work in order,
 * so the worker keeps at most QUEUE_DEPTH frames' steps in flight: after
 * queuing a frame's batch it waits for the one before to finish, which keeps the device
 * busy without letting the queue run ahead of the page's frames. The page
 * reports once a second how many of its frames came late for want of the
 * GPU, and the pacer (js/pace.module.js) turns that into an idle pause:
 * the worker lets the device drain and then waits that long. Without
 * reports (the page hidden) there is no pause.
 */
const QUEUE_DEPTH = 2;
const pacer = createPacer(), inFlight = [];
const pace = { pause: 0, heard: -Infinity };
function adjustPace({ late }) {
  pace.heard = performance.now();
  pace.pause = pacer.report(late);
}
const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));
async function yieldToPage() {
  inFlight.push(model.settle());
  if (pace.pause > 0 && performance.now() - pace.heard < 3000) { await Promise.all(inFlight.splice(0)); await sleep(pace.pause); }
  else if (inFlight.length >= QUEUE_DEPTH) await inFlight.shift();
}

/*
 * On the GPU a step only queues work: the frame's kernels and read-backs
 * are queued ahead of the frame's steps, which go to the device together
 * as one batch (model.stepBatch) yielding to the page once, and the frame
 * is posted after them. The CPU engines step their arrays in place, so they build the
 * frame after the steps. A land that started fresh jumps once its
 * record passes each of FRESH_JUMPS, as the spin-up's default does.
 */
async function loop() {
  if (!running) return;
  stepping = true;
  let finish;
  batchDone = new Promise((resolve) => { finish = resolve; });
  try {
    const landAge = model.land ? model.land.record[0] : -1;
    if (model.beginFrame) {
      const capturing = model.beginFrame(subscription);
      capturing.catch(() => {}); // when a step fails first, its error is the one reported
      await model.stepBatch(stepsPerFrame, dt, null, false);
      await yieldToPage();
      postFrame(await capturing);
    } else {
      for (let n = 0; n < stepsPerFrame; n++) await model.step(dt);
      await sendFrame();
    }
    if (model.land && freshJumpDue(model.land.record, landAge, model.land.record[0])) await model.land.jump();
  } catch (error) {
    running = false;
    report(error);
  } finally {
    stepping = false;
    finish();
  }
  if (running) setTimeout(loop, 0);
  else if (resend) refresh();
  resend = false;
}

const status = (text, fraction = null) => self.postMessage({ type: 'status', text, fraction });

/*
 * Fetches a saved state while reporting the bytes received against the
 * total, which is most of the wait for a large snapshot.
 */
async function fetchWithProgress(url, from, to) {
  const name = stateName(url.replace(/.*\//, '').replace(/[?#].*/, ''));
  let reported = -1;
  status(`loading ${name}…`, from);
  return fetchState(url, (received, total) => {
    if (!total) return;
    if (received >= total) { status(`parsing ${name}…`, to); return; }
    const percent = Math.floor(100 * received / total);
    if (percent !== reported) { reported = percent; status(`loading ${name}: ${(received / 1048576).toFixed(0)} of ${(total / 1048576).toFixed(0)} MB`, from + (to - from) * received / total); }
  });
}

/*
 * The initial state comes from a saved run (a *_state_*.json written by
 * the emergence driver) when the start message names one, regridded if
 * it was saved at another resolution and otherwise written straight into
 * the model's own state arrays; without one from initializeState. A saved
 * run's surface pressure is rebalanced from its terrain, taken to be the
 * current topography's, to the model's; the CPU engines' first frame
 * reads the core's diagnosis, so they diagnose even when the two match.
 */
function initialState(model, saved, N) {
  if (!saved) {
    const fresh = initializeState(model, { geostrophic: !model.surfaceGeopotential });
    if (model.geography) for (let i = 0; i < fresh[6].length; i++) if (model.geography.land[i]) fresh[6][i] = 0;
    return fresh;
  }
  const kept = [saved.pi, saved.theta, saved.u, saved.surfaceT];
  if (saved.q) kept.push(saved.q);
  if (saved.q && saved.qc) kept.push(saved.qc);
  if (saved.q && saved.qc && saved.ice) kept.push(saved.ice);
  model.time = saved.time;
  let carried, fromPhi = null;
  if (saved.N === N) {
    carried = kept.map((values, a) => { model.state[a].set(values); return model.state[a]; });
    if (saved.terrain) fromPhi = model.surfaceGeopotential ?? (model.geography ? surfaceGeopotential(model.mesh, model.geography) : null);
  } else {
    const source = sourceFor(saved);
    carried = regridState(source, model, kept.map((a) => Float64Array.from(a)), (fraction, text) => status(`regridding day ${saved.day} from N=${saved.N} to N=${N}: ${text}…`, 0.8 + 0.12 * fraction), { land: saved.land ?? null });
    if (saved.terrain) fromPhi = regridCellField(source, model, source.surfaceGeopotential);
  }
  if (fromPhi !== model.surfaceGeopotential || (fromPhi && !model.beginFrame)) {
    model.core.diagnose(carried[0], carried[1], null, null);
    rebalanceSurfacePressure(model.core, carried[0], carried[1], fromPhi, model.surfaceGeopotential);
  }
  if (carried.length < 5) carried.push(initialHumidity(model, carried[0], carried[1]));
  if (carried.length < 6) carried.push(new Float64Array(carried[1].length));
  if (carried.length < 7) carried.push(Float64Array.from(carried[3], (t) => (t < 271.35 ? 0.5 : 0)));
  if (model.geography) for (let i = 0; i < carried[6].length; i++) if (model.geography.land[i]) carried[6][i] = 0;
  return carried;
}

/*
 * The sea-ice concentration that goes with the placed ice: the saved
 * run's, regridded if it was saved at another resolution, or full cover
 * wherever a run saved without one has ice.
 */
function placeIce(model, saved, N) {
  const kept = !saved ? model.seaIce.concentration : !saved.concentration ? null : saved.N === N ? saved.concentration : regridConcentration(sourceFor(saved), model, saved.concentration, { land: saved.land ?? null, surfaceT: saved.surfaceT });
  model.seaIce.load(model.state[6], kept);
}

/*
 * The mixed-layer deck's carried state (running-mean subsidence,
 * inversion height and gate), the last means of the convective and
 * large-scale rain and those of the absorbed sunlight, outgoing longwave
 * and planetary albedo: the saved run's,
 * regridded if it was saved at another resolution, or each field's
 * starting value for a fresh start or a run saved without it.
 */
function placeDeck(model, saved, N) {
  for (const name of Object.keys(DECK_FIELDS)) model.radiation[name].set(savedDeckField(saved, name, model, saved && saved[name] && saved.N !== N ? sourceFor(saved) : null));
  for (const name of Object.keys(MOIST_FIELDS)) model.moist[name].set(savedMoistField(saved, name, model, saved && saved[name] && saved.N !== N ? sourceFor(saved) : null));
  for (const name of Object.keys(RADIATION_FIELDS)) model.radiation[name].set(savedRadiationField(saved, name, model, saved && saved[name] && saved.N !== N ? sourceFor(saved) : null));
  placeCumulus(model, saved, N);
}

/*
 * The cumulus cloud's cover and water are saved in the GPU's layout, the
 * plume layers from cumulusK0 down, (K − K0) × C; the CPU keeps K × C.
 * A run at another resolution starts without them, as a fresh start does.
 */
const CUMULUS_FIELDS = ['cumulusCover', 'cumulusWater'];
function placeCumulus(model, saved, N) {
  const C = model.mesh.nCells, K = model.core.K, K0 = model.moist.cumulusK0 ?? 0;
  for (const name of CUMULUS_FIELDS) {
    const target = model.moist[name];
    if (!target) continue;
    target.fill(0);
    const values = saved && saved.N === N ? saved[name] : null;
    if (!values) continue;
    if (values.length === target.length) target.set(values);
    else if (values.length === (K - K0) * C && target.length === K * C) target.set(values, K0 * C);
  }
}
function savedCumulus(model) {
  const C = model.mesh.nCells, K = model.core.K, K0 = model.moist.cumulusK0 ?? 0, out = {};
  for (const name of CUMULUS_FIELDS) {
    const values = model.moist[name];
    if (values) out[name] = Float64Array.from(values.length === K * C ? values.subarray(K0 * C) : values).buffer;
  }
  return out;
}

/*
 * A physics-free model on the saved run's mesh and sigma grid, with the
 * current topography so its land mask can steer the regrid; kept for
 * the ocean and land that follow the state, until the start is done.
 */
let sourceModel = null;
const sameLevels = (a, b) => a.length === b.length && a.every((sigma, k) => sigma === b[k]);
function sourceFor(saved) {
  const levels = savedLevels(saved);
  if (!sourceModel || sourceModel.N !== saved.N || sourceModel.topography !== currentTopography || !sameLevels(sourceModel.model.core.levels, levels)) {
    sourceModel = { N: saved.N, topography: currentTopography, model: createModel(new Grid(saved.N), { physics: false, levels, ...(currentTopography ? { topography: currentTopography } : {}) }) };
  }
  return sourceModel.model;
}

/*
 * The land state comes with a saved run when it has one, regridded if
 * needed; otherwise it starts as `initialize` makes it.
 */
function placeLand(model, saved, N) {
  if (!model.land) return;
  if (saved && saved.land) model.land.load(saved.N === N ? saved.land : regridLand(sourceFor(saved), model, saved.land, (fraction, text) => status(`regridding ${text}…`, 0.94), { ice: saved.ice ?? null, surfaceT: saved.surfaceT ?? null }), model.state[6]);
  else model.land.initialize();
}

let currentTopography = null;
const topographies = new Map();
async function loadTopography(url) {
  if (topographies.has(url)) return topographies.get(url);
  const response = await fetch(url);
  if (!response.ok) throw new Error(`topography ${url}: ${response.status}`);
  const topography = topographyFromInt16(await response.arrayBuffer());
  topographies.set(url, topography);
  return topography;
}

const subgrids = new Map();
async function loadSubgrid(N) {
  if (subgrids.has(N)) return subgrids.get(N);
  const response = await fetch(subgridUrl(N)).catch(() => null);
  const fields = response && response.ok ? await response.arrayBuffer().then(decodeSubgrid).catch(() => null) : null;
  subgrids.set(N, fields);
  return fields;
}

/*
 * The ocean climatology a fresh start takes its ocean and sea surface
 * from, fetched with the bytes received on the status line. Without a
 * `url` it is the repository's World Ocean Atlas file, and a start goes on
 * from the analytic climatology when that cannot be had.
 */
const climatologies = new Map();
async function loadOceanClimatology(url, from, to) {
  const source = url ?? new URL(`../${CLIMATOLOGY_FILE}`, import.meta.url).href;
  if (climatologies.has(source)) return climatologies.get(source);
  try {
    status('loading the ocean climatology…', from);
    const response = await fetch(source);
    if (!response.ok) throw new Error(`ocean climatology ${source}: ${response.status}`);
    const total = Number(response.headers.get('content-length')) || 0, parts = [];
    let received = 0;
    for (const reader = response.body.getReader(); ;) {
      const { done, value } = await reader.read();
      if (done) break;
      parts.push(value); received += value.length;
      if (total) status(`loading the ocean climatology: ${(received / 1048576).toFixed(1)} of ${(total / 1048576).toFixed(1)} MB`, from + (to - from) * Math.min(1, received / total));
    }
    const bytes = new Uint8Array(received);
    let at = 0;
    for (const part of parts) { bytes.set(part, at); at += part.length; }
    const climatology = decodeClimatology(bytes);
    climatologies.set(source, climatology);
    return climatology;
  } catch (error) {
    if (url) throw error;
    console.warn(`starting from the analytic ocean: ${error}`);
    return null;
  }
}

/*
 * The geography the page draws once: the land mask, the land fraction
 * and elevation of every cell, and the coast as segments between the
 * vertices of each edge that separates land from ocean, each with its
 * land cell for the projection.
 */
function geographyMessage(model) {
  const { geography, mesh } = model;
  if (!geography) return { land: null };
  const { coastEdges } = geography;
  const coast = new Float32Array(6 * coastEdges.length), coastCells = new Int32Array(coastEdges.length);
  for (let n = 0; n < coastEdges.length; n++) {
    const e = coastEdges[n];
    for (let side = 0; side < 2; side++) {
      const v = mesh.verticesOnEdge[2 * e + side];
      const r = Math.hypot(mesh.xVertex[3 * v], mesh.xVertex[3 * v + 1], mesh.xVertex[3 * v + 2]);
      for (let axis = 0; axis < 3; axis++) coast[6 * n + 3 * side + axis] = mesh.xVertex[3 * v + axis] / r;
    }
    const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1];
    coastCells[n] = geography.land[a] ? a : b;
  }
  return { land: Uint8Array.from(geography.land), landFraction: Float32Array.from(geography.landFraction), elevation: Float32Array.from(geography.elevation), coast, coastCells, terrain: !!model.surfaceGeopotential };
}

self.onmessage = async (event) => {
  const message = event.data;
  if (message.type === 'start') {
    await hold(() => start(message), { resume: !message.paused }).catch(report);
    if (pendingRestore) {
      const saved = pendingRestore;
      pendingRestore = null;
      if (model) await hold(() => restore(saved), { resume: running }).catch(report);
      else report('the snapshot was not restored: the model failed to start');
    }
  } else if (message.type === 'pause') {
    if (held) held.resume = false;
    else running = false;
  } else if (message.type === 'resume') {
    if (held) held.resume = true;
    else if (model && !running) { running = true; if (!stepping) loop(); }
  } else if (message.type === 'pace') {
    adjustPace(message);
  } else if (message.type === 'subscribe') {
    subscription = { level: 'surface', depth: 'surface', fields: [], diagnostics: false, ...message.subscription };
    if (serving && !running) refresh();
  } else if (message.type === 'snapshot') {
    if (serving) await hold(snapshot).catch(report);
    else report('the model is not ready yet');
  } else if (message.type === 'restore') {
    if (lastStart) await hold(() => restore(message.snapshot)).catch(report);
    else pendingRestore = message.snapshot;
  } else if (message.type === 'profile') {
    if (serving) await hold(profile, { resume: running }).catch((error) => self.postMessage({ type: 'profile', error: String(error && error.stack ? error.stack : error) }));
    else self.postMessage({ type: 'profile', error: 'the model is not ready yet' });
  } else if (message.type === 'probe') {
    await hold(() => probe(message)).catch((error) => self.postMessage({ type: 'probe', error: String(error && error.stack ? error.stack : error) }));
  }
};

/*
 * The device test the page runs before its first start on a device: a
 * model on the GPU at the message's N (64 on desktops, 32 on phones) from
 * the fresh initial state, stepped as the loop steps it after a short
 * warm-up, reporting milliseconds a step. Without a GPU, or when that
 * fails, one thread of the CPU engine at N=16 instead. It runs alone:
 * anything else on the worker's thread slows the steps.
 */
const PROBE = { gpuN: 64, cpuN: 16, warmup: 4, steps: 12, cpuSteps: 3 };
async function probe(message) {
  const result = { gpu: null, cpu: null }, gpuN = message.N ?? PROBE.gpuN;
  const topography = message.land === false ? null : await loadTopography(message.topography ?? new URL('../data/topography_0p25.bin', import.meta.url).href);
  const options = { ...(topography ? { topography } : {}), terrain: message.terrain !== false };
  const prepare = (test, N) => {
    for (const [a, values] of initialState(test, null, N).entries()) test.state[a].set(values);
    placeIce(test, null, N);
    placeDeck(test, null, N);
    if (test.load) test.load();
    if (test.ocean) test.ocean.initialize(test.state[3], test.state[6]);
    if (test.land) test.land.initialize();
  };
  if (message.engine === 'gpu' && typeof navigator !== 'undefined' && navigator.gpu) {
    let test = null;
    try {
      status('testing the GPU…', 0.1);
      test = await createGpuModel(new Grid(gpuN), { ...options, ...withCadence({}, 1350 * 16 / gpuN, gpuN), ...(topography ? { subgrid: await loadSubgrid(gpuN) } : {}) });
      prepare(test, gpuN);
      const step = 1350 * 16 / gpuN, queued = [];
      for (let n = 0; n < PROBE.warmup; n++) await test.step(step);
      await test.settle();
      const begin = performance.now();
      for (let n = 0; n < PROBE.steps; n++) {
        await test.step(step);
        queued.push(test.settle());
        if (queued.length >= QUEUE_DEPTH) await queued.shift();
      }
      await Promise.all(queued);
      result.gpu = { N: gpuN, ms: (performance.now() - begin) / PROBE.steps };
    } catch (error) {
      result.gpu = { N: gpuN, error: String(error) };
    } finally {
      if (test && test.destroy) test.destroy();
    }
  }
  if (!result.gpu || result.gpu.error) {
    status('testing the CPU…', 0.3);
    const test = createModel(new Grid(PROBE.cpuN), { ...options, ...withCadence({}, 1350 * 16 / PROBE.cpuN, PROBE.cpuN), ...(topography ? { subgrid: await loadSubgrid(PROBE.cpuN) } : {}) });
    prepare(test, PROBE.cpuN);
    const step = 1350 * 16 / PROBE.cpuN;
    await test.step(step);
    const begin = performance.now();
    for (let n = 0; n < PROBE.cpuSteps; n++) await test.step(step);
    result.cpu = { N: PROBE.cpuN, ms: (performance.now() - begin) / PROBE.cpuSteps };
  }
  self.postMessage({ type: 'probe', result });
}

/*
 * A GPU profile for the page (js/gpu/profile.module.js), with the time
 * one frame of the current subscription takes to compute and read back.
 * It runs while holding the model.
 */
async function profile() {
  if (!model.beginFrame) { self.postMessage({ type: 'profile', error: 'the profile needs the GPU engine' }); return; }
  const result = await profileGpu(model, { steps: 16, dt });
  const start = performance.now();
  const captured = await model.beginFrame(subscription);
  result.frame = performance.now() - start;
  postFrame(captured);
  self.postMessage({ type: 'profile', result: { ...result, N: currentN, dt, queueDepth: QUEUE_DEPTH, pause: pace.pause } });
}

let lastStart = null, currentN = null;

/*
 * While holding the model: refreshes the mirrors from the device if the
 * engine keeps them there, and hands the page a copy of the state and the
 * ocean as transferable buffers.
 */
async function snapshot() {
  if (model.sync) await model.sync();
  const names = ['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'];
  const arrays = Object.fromEntries(names.map((name, a) => [name, Float64Array.from(model.state[a]).buffer]));
  arrays.concentration = Float64Array.from(model.seaIce.concentration).buffer;
  for (const name of Object.keys(DECK_FIELDS)) arrays[name] = Float64Array.from(model.radiation[name]).buffer;
  for (const name of Object.keys(MOIST_FIELDS)) arrays[name] = Float64Array.from(model.moist[name]).buffer;
  Object.assign(arrays, savedCumulus(model));
  for (const name of Object.keys(RADIATION_FIELDS)) arrays[name] = Float64Array.from(model.radiation[name]).buffer;
  arrays.levels = Float64Array.from(model.core.levels).buffer;
  let ocean = null, land = null;
  if (model.ocean) {
    const o = await model.ocean.serialize();
    ocean = Object.fromEntries(Object.entries(o).map(([k, v]) => [k, Float64Array.from(v).buffer]));
  }
  if (model.land) {
    const l = await model.land.serialize();
    land = Object.fromEntries(Object.entries(l).map(([k, v]) => [k, Float64Array.from(v).buffer]));
  }
  const transfer = [...Object.values(arrays), ...(ocean ? Object.values(ocean) : []), ...(land ? Object.values(land) : [])];
  self.postMessage({ type: 'snapshotData', N: currentN, K: model.core.K, day: model.time / 86400, time: model.time, terrain: !!model.surfaceGeopotential, arrays, ocean, land }, transfer);
  refresh();
}

/*
 * Restores a snapshot while holding the model: into the running model
 * when the resolution and the sigma grid match, otherwise by starting
 * over with the snapshot as the saved state. The model stays paused
 * afterwards unless the page asked meanwhile for it to run.
 */
async function restore(snapshot) {
  const saved = { N: snapshot.N, K: snapshot.K, day: snapshot.day, time: snapshot.time, terrain: !!snapshot.terrain };
  for (const [name, buffer] of Object.entries(snapshot.arrays)) saved[name] = new Float64Array(buffer);
  if (snapshot.ocean) saved.ocean = Object.fromEntries(Object.entries(snapshot.ocean).map(([k, buffer]) => [k, new Float64Array(buffer)]));
  if (snapshot.land) saved.land = Object.fromEntries(Object.entries(snapshot.land).map(([k, buffer]) => [k, new Float64Array(buffer)]));
  if (!model || saved.N !== currentN || !sameLevels(savedLevels(saved), model.core.levels)) { await start({ ...lastStart, saved, paused: true }); return; }
  serving = false;
  status('restoring the snapshot…', 0.8);
  const init = initialState(model, saved, currentN);
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  placeIce(model, saved, currentN);
  placeDeck(model, saved, currentN);
  if (model.load) model.load();
  if (model.dragsDue) model.dragsDue();
  if (model.ocean) { if (saved.ocean) model.ocean.load(saved.ocean, model.state[3], model.state[6]); else model.ocean.initialize(model.state[3], model.state[6]); }
  placeLand(model, saved, currentN);
  placeEnergy(model, saved);
  model.time = saved.time;
  restartRain();
  serving = true;
  self.postMessage({ type: 'ready', N: currentN, cells: model.mesh.nCells, layers: model.core.K, levels: Array.from(model.core.levels), dt, day: model.time / 86400, workers: lastStart?.workers ?? 1, ocean: !!model.ocean, ...geographyMessage(model) });
  await sendFrame();
  sourceModel = null;
}

async function start(message) {
  lastStart = { ...message, saved: null };
  let saved = message.saved ?? null;
  if (!saved && message.from) {
    saved = await fetchWithProgress(message.from, 0, 0.5);
  }
  const N = message.N ?? saved?.N ?? 16;
  const gpuWanted = message.engine === 'gpu' && typeof navigator !== 'undefined' && navigator.gpu;
  const options = { ...(message.options ?? {}) };
  if (message.land !== false) { status('loading the topography…', 0.52); options.topography = await loadTopography(message.topography ?? new URL('../data/topography_0p25.bin', import.meta.url).href); options.subgrid = await loadSubgrid(N); }
  if (!saved && message.land !== false && message.climatology !== false) {
    const climatology = await loadOceanClimatology(message.climatology ?? null, 0.53, 0.55);
    if (climatology) options.ocean = { ...(options.ocean ?? {}), climatology };
  }
  options.terrain = message.terrain !== false;
  if (saved) options.levels = savedLevels(saved);
  else if (message.levels) options.levels = sigmaInterfaces(message.levels);
  const workers = message.workers ?? 1;
  const step = message.dt ?? 1350 * 16 / N;
  Object.assign(options, withCadence(options, step, N));
  status(`building the N=${N} grid…`, 0.55);
  const grid = new Grid(N);
  status(gpuWanted ? 'compiling the GPU model…' : workers > 1 ? `starting ${workers} workers…` : 'building the model…', 0.65);
  const built = gpuWanted ? await createGpuModel(grid, options) : workers > 1 ? await createParallelModel(grid, options, workers) : createModel(grid, options);
  serving = false;
  model = built;
  currentN = N;
  currentTopography = options.topography ?? null;
  dt = step;
  stepsPerFrame = message.stepsPerFrame ?? Math.max(2, Math.round((gpuWanted ? 24 : 8) * 16 / N));
  status(saved ? 'placing the saved state…' : 'building the initial state…', 0.8);
  const init = initialState(model, saved, N);
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  placeIce(model, saved, N);
  placeDeck(model, saved, N);
  status('uploading the state…', 0.92);
  if (model.load) model.load();
  if (model.ocean) {
    if (saved && saved.ocean) model.ocean.load(saved.N === N ? saved.ocean : regridOcean(sourceFor(saved), model, saved.ocean, (fraction, text) => status(`regridding ${text}…`, 0.93)), model.state[3], model.state[6]);
    else model.ocean.initialize(model.state[3], model.state[6]);
  }
  placeLand(model, saved, N);
  placeEnergy(model, saved);
  layerWinds = new Array(model.core.K).fill(null);
  if (message.subscription) subscription = { ...subscription, ...message.subscription };
  restartRain();
  serving = true;
  self.postMessage({ type: 'ready', N, cells: model.mesh.nCells, layers: model.core.K, levels: Array.from(model.core.levels), dt, day: model.time / 86400, workers, ocean: !!model.ocean, ...geographyMessage(model) });
  refresh();
  sourceModel = null;
}
