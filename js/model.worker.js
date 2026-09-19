import { Grid } from './grid.module.js';
import { createModel } from './model.module.js';
import { createParallelModel } from './parallel.module.js';
import { initializeState } from './physics/init.module.js';
import { cellVector } from './dynamics/operators.module.js';
import { regridState } from './physics/regrid.module.js';
import { levelFields } from './levels.module.js';

let model = null, running = false, dt = 450, stepsPerFrame = 24, frame = 0;
let level = 'surface', layerWinds = [];

/*
 * A frame carries the fields of the selected pressure level (or the
 * lowest layer) plus surface pressure; the cell-center winds of a layer
 * are reconstructed only when that layer bounds the level somewhere.
 */
function postFrame() {
  const { mesh, core, state, time } = model;
  const [pi, theta, u, surfaceT] = state;
  const E = mesh.nEdges;
  layerWinds.fill(null);
  const layerWind = (k) => layerWinds[k] ??= cellVector(mesh, u.subarray(k * E, (k + 1) * E), new Float64Array(3 * mesh.nCells));
  const fields = levelFields(core, pi, theta, layerWind, level);
  const ps = Float32Array.from(pi), ts = Float32Array.from(surfaceT);
  const diagnostics = model.diagnostics();
  const message = { type: 'frame', frame: frame++, time, day: time / 86400, level, ps, ts, ...fields, diagnostics };
  self.postMessage(message, [ps.buffer, ts.buffer, fields.speed.buffer, fields.vector.buffer, fields.temperature.buffer, fields.height.buffer]);
}

function loop() {
  if (!running) return;
  for (let n = 0; n < stepsPerFrame; n++) model.step(dt);
  postFrame();
  setTimeout(loop, 0);
}

const status = (text) => self.postMessage({ type: 'status', text });

/*
 * The initial state comes from a saved run (a *_state_*.json written by
 * the emergence driver) when the start message names one, regridded if
 * it was saved at another resolution; otherwise from initializeState.
 */
function initialState(model, saved, N) {
  if (!saved) return initializeState(model, {});
  const arrays = [saved.pi, saved.theta, saved.u, saved.surfaceT].map((a) => Float64Array.from(a));
  model.time = saved.time;
  if (saved.N === N) return arrays;
  status(`regridding day ${saved.day} from N=${saved.N} to N=${N}…`);
  return regridState(createModel(new Grid(saved.N)), model, arrays);
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
  }
};

async function start(message) {
  let saved = null;
  if (message.from) {
    status(`loading ${message.from.replace(/.*\//, '')}…`);
    saved = await (await fetch(message.from)).json();
  }
  const N = message.N ?? saved?.N ?? 16;
  dt = message.dt ?? 450 * 16 / N;
  stepsPerFrame = message.stepsPerFrame ?? Math.max(2, Math.round(24 * 16 / N));
  const workers = message.workers ?? 1;
  status(`building the N=${N} grid${workers > 1 ? ` and ${workers} workers` : ''}…`);
  const grid = new Grid(N);
  model = workers > 1 ? await createParallelModel(grid, message.options ?? {}, workers) : createModel(grid, message.options ?? {});
  const init = initialState(model, saved, N);
  for (let a = 0; a < 4; a++) model.state[a].set(init[a]);
  layerWinds = new Array(model.core.K).fill(null);
  level = message.level ?? 'surface';
  self.postMessage({ type: 'ready', N, cells: model.mesh.nCells, layers: model.core.K, dt, day: model.time / 86400, workers });
  postFrame();
  running = !message.paused;
  if (running) loop();
}
