import { Grid } from './grid.module.js';
import { createModel } from './model.module.js';
import { createParallelModel } from './parallel.module.js';
import { createGpuModel } from './gpu/model.gpu.js';
import { initializeState } from './physics/init.module.js';
import { cellVector } from './dynamics/operators.module.js';
import { regridState, regridOcean } from './physics/regrid.module.js';
import { levelFields } from './levels.module.js';
import { initialHumidity } from './physics/init.module.js';

let model = null, running = false, dt = 450, stepsPerFrame = 24, frame = 0;
let level = 'surface', layerWinds = [], lastFrameTime = 0;

/*
 * A frame carries the fields of the selected pressure level (or the
 * lowest layer) plus surface pressure; the cell-center winds of a layer
 * are reconstructed only when that layer bounds the level somewhere.
 */
async function postFrame() {
  if (model.sync) await model.sync();
  const diagnostics = await model.diagnostics();
  const { mesh, core, state, time } = model;
  const [pi, theta, u, surfaceT] = state;
  const E = mesh.nEdges;
  layerWinds.fill(null);
  const layerWind = (k) => layerWinds[k] ??= cellVector(mesh, u.subarray(k * E, (k + 1) * E), new Float64Array(3 * mesh.nCells));
  const q = state[4], qc = state[5];
  const ice = Float32Array.from(state[6]), albedo = Float32Array.from(state[6], (h) => model.seaIce.albedo(h));
  const shortwave = Float32Array.from(model.radiation.surfaceShortwave), longwave = Float32Array.from(model.radiation.outgoing);
  const precipitation = Float32Array.from(model.moist.precipitation);
  const fields = levelFields(core, pi, theta, layerWind, level, q);
  const ps = Float32Array.from(pi), ts = Float32Array.from(surfaceT);
  const water = new Float32Array(mesh.nCells), cloud = new Float32Array(mesh.nCells);
  for (let i = 0; i < mesh.nCells; i++) { water[i] = model.moist.columnWater(pi, q, i); cloud[i] = model.moist.columnWater(pi, qc, i); }
  const interval = time - lastFrameTime;
  lastFrameTime = time;
  for (let i = 0; i < mesh.nCells; i++) precipitation[i] = interval > 0 ? precipitation[i] / interval * 86400 : 0;
  const message = { type: 'frame', frame: frame++, time, day: time / 86400, level, ps, ts, ...fields, precipitation, water, cloud, ice, albedo, shortwave, longwave, diagnostics, engine: model.engine ?? 'cpu' };
  self.postMessage(message, [ps.buffer, ts.buffer, fields.speed.buffer, fields.vector.buffer, fields.temperature.buffer, fields.height.buffer, fields.humidity.buffer, precipitation.buffer, water.buffer, cloud.buffer, ice.buffer, albedo.buffer, shortwave.buffer, longwave.buffer]);
}

async function loop() {
  if (!running) return;
  for (let n = 0; n < stepsPerFrame; n++) await model.step(dt);
  await postFrame();
  setTimeout(loop, 0);
}

const status = (text, fraction = null) => self.postMessage({ type: 'status', text, fraction });

/*
 * Fetches a JSON file while reporting the bytes received against the
 * Content-Length, which is most of the wait for a large snapshot.
 */
async function fetchWithProgress(url, from, to) {
  const response = await fetch(url);
  const total = Number(response.headers.get('content-length')) || 0;
  const name = url.replace(/.*\//, '');
  if (!response.body || !total) { status(`loading ${name}…`, from); return response.json(); }
  const reader = response.body.getReader();
  const chunks = [];
  let received = 0, reported = -1;
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    chunks.push(value);
    received += value.length;
    const percent = Math.floor(100 * received / total);
    if (percent !== reported) { reported = percent; status(`loading ${name}: ${(received / 1048576).toFixed(0)} of ${(total / 1048576).toFixed(0)} MB`, from + (to - from) * received / total); }
  }
  const bytes = new Uint8Array(received);
  let at = 0;
  for (const chunk of chunks) { bytes.set(chunk, at); at += chunk.length; }
  status(`parsing ${name}…`, to);
  return JSON.parse(new TextDecoder().decode(bytes));
}

/*
 * The initial state comes from a saved run (a *_state_*.json written by
 * the emergence driver) when the start message names one, regridded if
 * it was saved at another resolution; otherwise from initializeState.
 */
function initialState(model, saved, N) {
  if (!saved) return initializeState(model, {});
  const arrays = [saved.pi, saved.theta, saved.u, saved.surfaceT].map((a) => Float64Array.from(a));
  if (saved.q) arrays.push(Float64Array.from(saved.q));
  if (saved.q && saved.qc) arrays.push(Float64Array.from(saved.qc));
  if (saved.q && saved.qc && saved.ice) arrays.push(Float64Array.from(saved.ice));
  model.time = saved.time;
  const carried = saved.N === N ? arrays : regridState(createModel(new Grid(saved.N)), model, arrays, (fraction, text) => status(`regridding day ${saved.day} from N=${saved.N} to N=${N}: ${text}…`, 0.8 + 0.12 * fraction));
  if (carried.length < 5) carried.push(initialHumidity(model, carried[0], carried[1]));
  if (carried.length < 6) carried.push(new Float64Array(carried[1].length));
  if (carried.length < 7) carried.push(Float64Array.from(carried[3], (t) => (t < 271.35 ? 0.5 : 0)));
  return carried;
}

self.onmessage = async (event) => {
  const message = event.data;
  if (message.type === 'start') {
    try { await start(message); } catch (error) { status(`error: ${error && error.stack ? error.stack : error}`); }
  } else if (message.type === 'pause') {
    running = false;
  } else if (message.type === 'resume') {
    if (model && !running) { running = true; loop(); }
  } else if (message.type === 'level') {
    level = message.level;
    if (model) postFrame();
  } else if (message.type === 'snapshot') {
    if (model) await snapshot();
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
  running = false;
  if (model.sync) await model.sync();
  const names = ['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'];
  const arrays = Object.fromEntries(names.map((name, a) => [name, Float64Array.from(model.state[a]).buffer]));
  let ocean = null;
  if (model.ocean) {
    const o = await model.ocean.serialize();
    ocean = Object.fromEntries(Object.entries(o).map(([k, v]) => [k, Float64Array.from(v).buffer]));
  }
  const transfer = [...Object.values(arrays), ...(ocean ? Object.values(ocean) : [])];
  self.postMessage({ type: 'snapshotData', N: currentN, K: model.core.K, day: model.time / 86400, time: model.time, arrays, ocean }, transfer);
  postFrame();
}

/*
 * Restores a snapshot: into the running model when the resolution
 * matches, otherwise by starting over with the snapshot as the saved
 * state. The model stays paused afterwards.
 */
async function restore(snapshot) {
  const saved = { N: snapshot.N, K: snapshot.K, day: snapshot.day, time: snapshot.time };
  for (const [name, buffer] of Object.entries(snapshot.arrays)) saved[name] = new Float64Array(buffer);
  if (snapshot.ocean) saved.ocean = Object.fromEntries(Object.entries(snapshot.ocean).map(([k, buffer]) => [k, new Float64Array(buffer)]));
  running = false;
  if (!model || saved.N !== currentN) { await start({ ...lastStart, saved, paused: true }); return; }
  status('restoring the snapshot…', 0.8);
  const init = initialState(model, saved, currentN);
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  if (model.load) model.load();
  if (model.ocean) { if (saved.ocean) model.ocean.load(saved.ocean, model.state[3], model.state[6]); else model.ocean.initialize(model.state[3], model.state[6]); }
  model.time = saved.time;
  lastFrameTime = model.time;
  layerWinds.fill(null);
  self.postMessage({ type: 'ready', N: currentN, cells: model.mesh.nCells, layers: model.core.K, dt, day: model.time / 86400, workers: lastStart?.workers ?? 1 });
  await postFrame();
}

async function start(message) {
  lastStart = message;
  let saved = message.saved ?? null;
  if (!saved && message.from) {
    saved = await fetchWithProgress(message.from, 0, 0.5);
  }
  const N = message.N ?? saved?.N ?? 16;
  currentN = N;
  const gpuWanted = message.engine === 'gpu' && typeof navigator !== 'undefined' && navigator.gpu;
  dt = message.dt ?? 1350 * 16 / N;
  stepsPerFrame = message.stepsPerFrame ?? Math.max(2, Math.round((gpuWanted ? 24 : 8) * 16 / N));
  const workers = message.workers ?? 1;
  status(`building the N=${N} grid…`, 0.55);
  const grid = new Grid(N);
  status(gpuWanted ? 'compiling the GPU model…' : workers > 1 ? `starting ${workers} workers…` : 'building the model…', 0.65);
  model = gpuWanted ? await createGpuModel(grid, message.options ?? {}) : workers > 1 ? await createParallelModel(grid, message.options ?? {}, workers) : createModel(grid, message.options ?? {});
  status(saved ? 'placing the saved state…' : 'building the initial state…', 0.8);
  const init = initialState(model, saved, N);
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  status('uploading the state…', 0.92);
  if (model.load) model.load();
  if (model.ocean) {
    if (saved && saved.ocean) model.ocean.load(saved.N === N ? saved.ocean : regridOcean(createModel(new Grid(saved.N)), model, saved.ocean, (fraction, text) => status(`regridding ${text}…`, 0.93)), model.state[3], model.state[6]);
    else model.ocean.initialize(model.state[3], model.state[6]);
  }
  layerWinds = new Array(model.core.K).fill(null);
  level = message.level ?? 'surface';
  lastFrameTime = model.time;
  self.postMessage({ type: 'ready', N, cells: model.mesh.nCells, layers: model.core.K, dt, day: model.time / 86400, workers });
  postFrame();
  running = !message.paused;
  if (running) loop();
}
