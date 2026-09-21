import { Grid } from './grid.module.js';
import { createModel } from './model.module.js';
import { createParallelModel } from './parallel.module.js';
import { createGpuModel } from './gpu/model.gpu.js';
import { initializeState } from './physics/init.module.js';
import { cellVector } from './dynamics/operators.module.js';
import { regridState, regridOcean, regridLand } from './physics/regrid.module.js';
import { topographyFromInt16, rebalanceSurfacePressure } from './geography.module.js';
import { regridCellField } from './physics/regrid.module.js';
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
  const onLand = model.geography ? model.geography.land : null;
  const ice = Float32Array.from(state[6]), albedo = Float32Array.from(state[6], (h, i) => (onLand && onLand[i] ? model.land.albedo(i) : model.seaIce.albedo(h)));
  const soil = model.land ? Float32Array.from(model.land.soil) : new Float32Array(0), snow = model.land ? Float32Array.from(model.land.snow) : new Float32Array(0);
  const oceanFields = model.oceanFields ? model.oceanFields() : null;
  const sst = new Float32Array(oceanFields ? mesh.nCells : 0), current = new Float32Array(sst.length), currentVector = new Float32Array(3 * sst.length), layerDepth = new Float32Array(sst.length), thermocline = new Float32Array(sst.length);
  if (oceanFields) {
    const vector = cellVector(mesh, oceanFields.u1, new Float64Array(3 * mesh.nCells));
    for (let i = 0; i < mesh.nCells; i++) {
      const sea = !onLand || !onLand[i];
      sst[i] = sea ? oceanFields.T1[i] : NaN; layerDepth[i] = sea ? oceanFields.h1[i] : NaN; thermocline[i] = sea ? oceanFields.T2[i] : NaN;
      currentVector[3 * i] = vector[3 * i]; currentVector[3 * i + 1] = vector[3 * i + 1]; currentVector[3 * i + 2] = vector[3 * i + 2];
      current[i] = sea ? Math.hypot(vector[3 * i], vector[3 * i + 1], vector[3 * i + 2]) : NaN;
    }
  }
  const shortwave = Float32Array.from(model.radiation.surfaceShortwave), longwave = Float32Array.from(model.radiation.outgoing);
  const precipitation = Float32Array.from(model.moist.precipitation);
  const fields = levelFields(core, pi, theta, layerWind, level, q);
  const ps = Float32Array.from(pi), ts = Float32Array.from(surfaceT);
  const mslp = Float32Array.from(pi);
  if (model.surfaceGeopotential) {
    const phis = model.surfaceGeopotential, { R, g, exnerLayer } = core.diagnostics, bottom = (core.K - 1) * mesh.nCells;
    for (let i = 0; i < mesh.nCells; i++) mslp[i] = pi[i] * Math.exp(phis[i] / (R * (theta[bottom + i] * exnerLayer[bottom + i] + 0.00325 * phis[i] / g)));
  }
  const water = new Float32Array(mesh.nCells), cloud = new Float32Array(mesh.nCells);
  for (let i = 0; i < mesh.nCells; i++) { water[i] = model.moist.columnWater(pi, q, i); cloud[i] = model.moist.columnWater(pi, qc, i); }
  const interval = time - lastFrameTime;
  lastFrameTime = time;
  for (let i = 0; i < mesh.nCells; i++) precipitation[i] = interval > 0 ? precipitation[i] / interval * 86400 : 0;
  const message = { type: 'frame', frame: frame++, time, day: time / 86400, level, ps, mslp, ts, ...fields, precipitation, water, cloud, ice, albedo, shortwave, longwave, soil, snow, sst, current, currentVector, layerDepth, thermocline, diagnostics, engine: model.engine ?? 'cpu' };
  self.postMessage(message, [ps.buffer, mslp.buffer, ts.buffer, fields.speed.buffer, fields.vector.buffer, fields.temperature.buffer, fields.height.buffer, fields.humidity.buffer, precipitation.buffer, water.buffer, cloud.buffer, ice.buffer, albedo.buffer, shortwave.buffer, longwave.buffer, soil.buffer, snow.buffer, sst.buffer, current.buffer, currentVector.buffer, layerDepth.buffer, thermocline.buffer]);
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
  const source = saved.N === N && !saved.terrain ? null : createModel(new Grid(saved.N), { physics: false, ...(saved.terrain ? { topography: currentTopography } : {}) });
  const carried = saved.N === N ? arrays : regridState(source, model, arrays, (fraction, text) => status(`regridding day ${saved.day} from N=${saved.N} to N=${N}: ${text}…`, 0.8 + 0.12 * fraction));
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
 * The land state comes with a saved run when it has one, regridded if
 * needed; otherwise the buckets start half full and bare.
 */
function placeLand(model, saved, N) {
  if (!model.land) return;
  if (saved && saved.land) model.land.load(saved.N === N ? { soil: Float64Array.from(saved.land.soil), snow: Float64Array.from(saved.land.snow) } : regridLand(createModel(new Grid(saved.N)), model, saved.land, (fraction, text) => status(`regridding ${text}…`, 0.94)));
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
  postFrame();
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
  running = false;
  if (!model || saved.N !== currentN) { await start({ ...lastStart, saved, paused: true }); return; }
  status('restoring the snapshot…', 0.8);
  const init = initialState(model, saved, currentN);
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  if (model.load) model.load();
  if (model.ocean) { if (saved.ocean) model.ocean.load(saved.ocean, model.state[3], model.state[6]); else model.ocean.initialize(model.state[3], model.state[6]); }
  placeLand(model, saved, currentN);
  model.time = saved.time;
  lastFrameTime = model.time;
  layerWinds.fill(null);
  self.postMessage({ type: 'ready', N: currentN, cells: model.mesh.nCells, layers: model.core.K, dt, day: model.time / 86400, workers: lastStart?.workers ?? 1, ...geographyMessage(model) });
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
    if (saved && saved.ocean) model.ocean.load(saved.N === N ? saved.ocean : regridOcean(createModel(new Grid(saved.N)), model, saved.ocean, (fraction, text) => status(`regridding ${text}…`, 0.93)), model.state[3], model.state[6]);
    else model.ocean.initialize(model.state[3], model.state[6]);
  }
  placeLand(model, saved, N);
  layerWinds = new Array(model.core.K).fill(null);
  level = message.level ?? 'surface';
  lastFrameTime = model.time;
  self.postMessage({ type: 'ready', N, cells: model.mesh.nCells, layers: model.core.K, dt, day: model.time / 86400, workers, ...geographyMessage(model) });
  postFrame();
  running = !message.paused;
  if (running) loop();
}
