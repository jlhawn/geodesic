import { Grid } from './grid.module.js';
import { createModel } from './model.module.js';
import { createParallelModel } from './parallel.module.js';
import { initializeState } from './physics/init.module.js';
import { cellVector } from './dynamics/operators.module.js';
import { regridState } from './physics/regrid.module.js';

let model = null, running = false, dt = 450, stepsPerFrame = 24, frame = 0;
let vector = null, jetLayer = 0;

function layerWind(k) {
  const { mesh, core, state } = model;
  const C = mesh.nCells, E = mesh.nEdges;
  cellVector(mesh, state[2].subarray(k * E, (k + 1) * E), vector);
  const speed = new Float32Array(C);
  for (let i = 0; i < C; i++) speed[i] = Math.hypot(vector[3 * i], vector[3 * i + 1], vector[3 * i + 2]);
  return { speed, vector: Float32Array.from(vector) };
}

function postFrame() {
  const { core, state, time } = model;
  const [pi, , , surfaceT] = state;
  const surface = layerWind(core.K - 1), jet = layerWind(jetLayer);
  const ps = Float32Array.from(pi), ts = Float32Array.from(surfaceT);
  const diagnostics = model.diagnostics();
  const message = { type: 'frame', frame: frame++, time, day: time / 86400, ps, ts, wind: surface.speed, windVector: surface.vector, jet: jet.speed, jetVector: jet.vector, diagnostics };
  self.postMessage(message, [ps.buffer, ts.buffer, surface.speed.buffer, surface.vector.buffer, jet.speed.buffer, jet.vector.buffer]);
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
async function initialState(model, message) {
  if (!message.from) return initializeState(model, message.init ?? {});
  status(`loading ${message.from.replace(/.*\//, '')}…`);
  const saved = await (await fetch(message.from)).json();
  const arrays = [saved.pi, saved.theta, saved.u, saved.surfaceT].map((a) => Float64Array.from(a));
  model.time = saved.time;
  if (saved.N === message.N) return arrays;
  status(`regridding day ${saved.day} from N=${saved.N} to N=${message.N}…`);
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
  }
};

async function start(message) {
  const N = message.N ?? 16;
  dt = message.dt ?? 450 * 16 / N;
  stepsPerFrame = message.stepsPerFrame ?? Math.max(2, Math.round(24 * 16 / N));
  const workers = message.workers ?? 1;
  status(`building the N=${N} grid${workers > 1 ? ` and ${workers} workers` : ''}…`);
  const grid = new Grid(N);
  model = workers > 1 ? await createParallelModel(grid, message.options ?? {}, workers) : createModel(grid, message.options ?? {});
  const init = await initialState(model, message);
  for (let a = 0; a < 4; a++) model.state[a].set(init[a]);
  vector = new Float64Array(3 * model.mesh.nCells);
  jetLayer = model.core.sigmaMid.findIndex((s) => s > 0.25);
  self.postMessage({ type: 'ready', N, cells: model.mesh.nCells, layers: model.core.K, dt, day: model.time / 86400, workers });
  postFrame();
  running = true;
  loop();
}
