import { Grid } from './grid.module.js';
import { createModel } from './model.module.js';
import { initializeState } from './physics/init.module.js';
import { cellVector } from './dynamics/operators.module.js';

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

self.onmessage = (event) => {
  const message = event.data;
  if (message.type === 'start') {
    const N = message.N ?? 16;
    dt = message.dt ?? 450 * 16 / N;
    stepsPerFrame = message.stepsPerFrame ?? 24;
    const grid = new Grid(N);
    model = createModel(grid, message.options ?? {});
    const init = initializeState(model, message.init ?? {});
    for (let a = 0; a < 4; a++) model.state[a].set(init[a]);
    vector = new Float64Array(3 * model.mesh.nCells);
    jetLayer = model.core.sigmaMid.findIndex((s) => s > 0.25);
    self.postMessage({ type: 'ready', N, cells: model.mesh.nCells, layers: model.core.K, dt });
    postFrame();
    running = true;
    loop();
  } else if (message.type === 'pause') {
    running = false;
  } else if (message.type === 'resume') {
    if (model && !running) { running = true; loop(); }
  }
};
