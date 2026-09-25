import { Grid } from './grid.module.js';
import { createModel } from './model.module.js';
import { createParallelModel } from './parallel.module.js';
import { createGpuModel } from './gpu/model.gpu.js';
import { initializeState } from './physics/init.module.js';
import { cellVector } from './dynamics/operators.module.js';
import { regridState, regridOcean, regridLand } from './physics/regrid.module.js';
import { topographyFromInt16, rebalanceSurfacePressure } from './geography.module.js';
import { regridCellField } from './physics/regrid.module.js';
import { levelFields, dewPoint, wetBulb, miseryIndex } from './levels.module.js';
import { initialHumidity } from './physics/init.module.js';
import { fetchState, stateName } from './stateFile.module.js';
import { LEVEL_FIELDS, OCEAN_FIELDS, RAIN_MEMORY } from './frames.module.js';
import { createPacer } from './pace.module.js';

const FREEZING = 273.15;
let model = null, serving = false, running = false, dt = 450, stepsPerFrame = 24, frame = 0;
let subscription = { level: 'surface', fields: [], diagnostics: false }, layerWinds = [];
const rain = { total: null, time: 0 };

function restartRain() { rain.total = null; rain.time = model.time; if (model.restartPrecipitation) model.restartPrecipitation(); }

/*
 * A frame for the CPU engines, from the model's own arrays: the rain
 * folds into its three-hour memory every frame, and only the subscribed
 * fields and diagnostics are built. The cell-center winds of a layer are
 * reconstructed only when that layer bounds the level somewhere.
 */
async function cpuFrame({ level, fields, diagnostics: summarize }) {
  const { mesh, core, state } = model;
  const C = mesh.nCells, E = mesh.nEdges, [pi, theta, u, , q, qc, iceField] = state;
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
  if (want.has('cloud')) out.cloud = Float32Array.from({ length: C }, (_, i) => model.moist.columnWater(pi, qc, i));
  if (want.has('rain')) out.rain = Float32Array.from(rain.total);
  if (want.has('ice')) out.ice = Float32Array.from(iceField);
  if (want.has('albedo')) out.albedo = Float32Array.from(iceField, (h, i) => (onLand && onLand[i] ? model.land.albedo(i) : model.seaIce.albedo(h)));
  if (want.has('shortwave')) out.shortwave = Float32Array.from(model.radiation.surfaceShortwave);
  if (want.has('longwave')) out.longwave = Float32Array.from(model.radiation.outgoing);
  if (model.land && want.has('soil')) out.soil = Float32Array.from(model.land.soil);
  if (model.land && want.has('snow')) out.snow = Float32Array.from(model.land.snow);
  const ocean = model.oceanFields && fields.some((name) => OCEAN_FIELDS.has(name)) ? model.oceanFields() : null;
  if (ocean) {
    for (const [name, values] of Object.entries({ sst: ocean.T1, sss: ocean.S1, layerDepth: ocean.h1, thermocline: ocean.thermoclineDepth, ssh: ocean.eta })) if (want.has(name)) out[name] = sea(values);
    if (want.has('current') || want.has('currents')) {
      const vector = cellVector(mesh, ocean.u1, new Float64Array(3 * C));
      if (onLand) for (let i = 0; i < C; i++) if (onLand[i]) vector.fill(0, 3 * i, 3 * i + 3);
      if (want.has('currents')) out.currents = Float32Array.from(vector);
      if (want.has('current')) out.current = sea(Float32Array.from({ length: C }, (_, i) => Math.hypot(vector[3 * i], vector[3 * i + 1], vector[3 * i + 2])));
    }
  }
  return { time, level, fields: out, diagnostics };
}

const captureFrame = () => (model.beginFrame ? model.beginFrame(subscription) : cpuFrame(subscription));

function postFrame({ time, level, fields, diagnostics }) {
  const transfer = [...new Set(Object.values(fields).map((values) => values.buffer))];
  self.postMessage({ type: 'frame', frame: frame++, time, day: time / 86400, level, engine: model.engine ?? 'cpu', pause: model.beginFrame ? pace.pause : 0, fields, diagnostics }, transfer);
}

async function sendFrame() { postFrame(await captureFrame()); }

/*
 * A frame outside the loop, reporting any failure on the status line.
 * While a batch of steps is still finishing after a pause, it waits for
 * that batch's own frame, so that the page's last frame carries the
 * latest subscription. halt() pauses and resolves once no batch is
 * stepping, before anything reads or replaces the model's state.
 */
let stepping = false, resend = false, batchDone = Promise.resolve();
async function halt() { running = false; resend = false; await batchDone; }
const report = (error) => status(`error: ${error && error.stack ? error.stack : error}`);
function refresh() {
  if (stepping) { resend = true; return; }
  sendFrame().catch(report);
}

/*
 * The GPU draws the page's globe too, and it takes queued work in order,
 * so the worker keeps at most QUEUE_DEPTH steps in flight: after queuing
 * a step it waits for the one before to finish, which keeps the device
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
 * are queued ahead of the batch of steps and the frame is posted after
 * it. The CPU engines step their arrays in place, so they build the
 * frame after the steps.
 */
async function loop() {
  if (!running) return;
  stepping = true;
  let finish;
  batchDone = new Promise((resolve) => { finish = resolve; });
  try {
    if (model.beginFrame) {
      const capturing = model.beginFrame(subscription);
      capturing.catch(() => {}); // when a step fails first, its error is the one reported
      for (let n = 0; n < stepsPerFrame; n++) { await model.step(dt); await yieldToPage(); }
      postFrame(await capturing);
    } else {
      for (let n = 0; n < stepsPerFrame; n++) await model.step(dt);
      await sendFrame();
    }
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
 * it was saved at another resolution; otherwise from initializeState.
 */
function initialState(model, saved, N) {
  if (!saved) {
    const fresh = initializeState(model, { geostrophic: !model.surfaceGeopotential });
    if (model.geography) for (let i = 0; i < fresh[6].length; i++) if (model.geography.land[i]) fresh[6][i] = 0;
    return fresh;
  }
  const arrays = [saved.pi, saved.theta, saved.u, saved.surfaceT].map((a) => Float64Array.from(a));
  if (saved.q) arrays.push(Float64Array.from(saved.q));
  if (saved.q && saved.qc) arrays.push(Float64Array.from(saved.qc));
  if (saved.q && saved.qc && saved.ice) arrays.push(Float64Array.from(saved.ice));
  model.time = saved.time;
  const source = saved.N !== N || saved.terrain ? sourceFor(saved) : null;
  const carried = saved.N === N ? arrays : regridState(source, model, arrays, (fraction, text) => status(`regridding day ${saved.day} from N=${saved.N} to N=${N}: ${text}…`, 0.8 + 0.12 * fraction), { land: saved.land ?? null });
  const fromPhi = saved.terrain ? (saved.N === N ? source.surfaceGeopotential : regridCellField(source, model, source.surfaceGeopotential)) : null;
  if (fromPhi || model.surfaceGeopotential) {
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
 * A physics-free model on the saved run's mesh, with the current
 * topography so its land mask can steer the regrid; kept for the ocean
 * and land that follow the state.
 */
let sourceModel = null;
function sourceFor(saved) {
  if (!sourceModel || sourceModel.N !== saved.N || sourceModel.topography !== currentTopography) {
    sourceModel = { N: saved.N, topography: currentTopography, model: createModel(new Grid(saved.N), { physics: false, ...(currentTopography ? { topography: currentTopography } : {}) }) };
  }
  return sourceModel.model;
}

/*
 * The land state comes with a saved run when it has one, regridded if
 * needed; otherwise the buckets start half full and bare.
 */
function placeLand(model, saved, N) {
  if (!model.land) return;
  if (saved && saved.land) model.land.load(saved.N === N ? { soil: Float64Array.from(saved.land.soil), snow: Float64Array.from(saved.land.snow) } : regridLand(sourceFor(saved), model, saved.land, (fraction, text) => status(`regridding ${text}…`, 0.94), { ice: saved.ice ?? null, surfaceT: saved.surfaceT ?? null }));
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
    try { await start(message); } catch (error) { status(`error: ${error && error.stack ? error.stack : error}`); }
  } else if (message.type === 'pause') {
    running = false;
  } else if (message.type === 'resume') {
    if (model && !running) { running = true; if (!stepping) loop(); }
  } else if (message.type === 'pace') {
    adjustPace(message);
  } else if (message.type === 'subscribe') {
    subscription = { level: 'surface', fields: [], diagnostics: false, ...message.subscription };
    if (serving && !running) refresh();
  } else if (message.type === 'snapshot') {
    if (serving) await snapshot();
  } else if (message.type === 'restore') {
    try { await restore(message.snapshot); } catch (error) { status(`error: ${error && error.stack ? error.stack : error}`); }
  }
};

let lastStart = null, currentN = null;

/*
 * Pauses, refreshes the mirrors from the device if the engine keeps
 * them there, and hands the page a copy of the state and the ocean as
 * transferable buffers.
 */
async function snapshot() {
  await halt();
  if (model.sync) await model.sync();
  const names = ['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'];
  const arrays = Object.fromEntries(names.map((name, a) => [name, Float64Array.from(model.state[a]).buffer]));
  let ocean = null, land = null;
  if (model.ocean) {
    const o = await model.ocean.serialize();
    ocean = Object.fromEntries(Object.entries(o).map(([k, v]) => [k, Float64Array.from(v).buffer]));
  }
  if (model.land) {
    const l = await model.land.serialize();
    land = { soil: Float64Array.from(l.soil).buffer, snow: Float64Array.from(l.snow).buffer };
  }
  const transfer = [...Object.values(arrays), ...(ocean ? Object.values(ocean) : []), ...(land ? Object.values(land) : [])];
  self.postMessage({ type: 'snapshotData', N: currentN, K: model.core.K, day: model.time / 86400, time: model.time, terrain: !!model.surfaceGeopotential, arrays, ocean, land }, transfer);
  refresh();
}

/*
 * Restores a snapshot: into the running model when the resolution
 * matches, otherwise by starting over with the snapshot as the saved
 * state. The model stays paused afterwards.
 */
async function restore(snapshot) {
  const saved = { N: snapshot.N, K: snapshot.K, day: snapshot.day, time: snapshot.time, terrain: !!snapshot.terrain };
  for (const [name, buffer] of Object.entries(snapshot.arrays)) saved[name] = new Float64Array(buffer);
  if (snapshot.ocean) saved.ocean = Object.fromEntries(Object.entries(snapshot.ocean).map(([k, buffer]) => [k, new Float64Array(buffer)]));
  if (snapshot.land) saved.land = Object.fromEntries(Object.entries(snapshot.land).map(([k, buffer]) => [k, new Float64Array(buffer)]));
  await halt();
  if (!model || saved.N !== currentN) { await start({ ...lastStart, saved, paused: true }); return; }
  serving = false;
  status('restoring the snapshot…', 0.8);
  const init = initialState(model, saved, currentN);
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  if (model.load) model.load();
  if (model.ocean) { if (saved.ocean) model.ocean.load(saved.ocean, model.state[3], model.state[6]); else model.ocean.initialize(model.state[3], model.state[6]); }
  placeLand(model, saved, currentN);
  model.time = saved.time;
  restartRain();
  serving = true;
  self.postMessage({ type: 'ready', N: currentN, cells: model.mesh.nCells, layers: model.core.K, dt, day: model.time / 86400, workers: lastStart?.workers ?? 1, ocean: !!model.ocean, ...geographyMessage(model) });
  await sendFrame();
}

async function start(message) {
  lastStart = message;
  serving = false;
  let saved = message.saved ?? null;
  if (!saved && message.from) {
    saved = await fetchWithProgress(message.from, 0, 0.5);
  }
  const N = message.N ?? saved?.N ?? 16;
  currentN = N;
  const gpuWanted = message.engine === 'gpu' && typeof navigator !== 'undefined' && navigator.gpu;
  const options = { ...(message.options ?? {}) };
  if (message.land !== false) { status('loading the topography…', 0.52); options.topography = await loadTopography(message.topography ?? new URL('../data/topography_0p25.bin', import.meta.url).href); }
  currentTopography = options.topography ?? null;
  options.terrain = message.terrain !== false;
  dt = message.dt ?? 1350 * 16 / N;
  stepsPerFrame = message.stepsPerFrame ?? Math.max(2, Math.round((gpuWanted ? 24 : 8) * 16 / N));
  const workers = message.workers ?? 1;
  status(`building the N=${N} grid…`, 0.55);
  const grid = new Grid(N);
  status(gpuWanted ? 'compiling the GPU model…' : workers > 1 ? `starting ${workers} workers…` : 'building the model…', 0.65);
  model = gpuWanted ? await createGpuModel(grid, options) : workers > 1 ? await createParallelModel(grid, options, workers) : createModel(grid, options);
  status(saved ? 'placing the saved state…' : 'building the initial state…', 0.8);
  const init = initialState(model, saved, N);
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  status('uploading the state…', 0.92);
  if (model.load) model.load();
  if (model.ocean) {
    if (saved && saved.ocean) model.ocean.load(saved.N === N ? saved.ocean : regridOcean(sourceFor(saved), model, saved.ocean, (fraction, text) => status(`regridding ${text}…`, 0.93)), model.state[3], model.state[6]);
    else model.ocean.initialize(model.state[3], model.state[6]);
  }
  placeLand(model, saved, N);
  layerWinds = new Array(model.core.K).fill(null);
  if (message.subscription) subscription = { ...subscription, ...message.subscription };
  restartRain();
  serving = true;
  self.postMessage({ type: 'ready', N, cells: model.mesh.nCells, layers: model.core.K, dt, day: model.time / 86400, workers, ocean: !!model.ocean, ...geographyMessage(model) });
  refresh();
  running = !message.paused;
  if (running) loop();
}
