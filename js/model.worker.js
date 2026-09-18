import { Grid } from './grid.module.js';
import { createModel } from './model.module.js';
import { initializeState } from './physics/init.module.js';
import { cellVector } from './dynamics/operators.module.js';

let model = null, running = false, dt = 450, stepsPerFrame = 24, frame = 0;
let surfaceWind = null, vector = null;

function postFrame() {
  const { mesh, core, state, time } = model;
  const [pi, , u, surfaceT] = state;
  const C = mesh.nCells, E = mesh.nEdges, K = core.K;
  cellVector(mesh, u.subarray((K - 1) * E, K * E), vector);
  for (let i = 0; i < C; i++) surfaceWind[i] = Math.hypot(vector[3 * i], vector[3 * i + 1], vector[3 * i + 2]);
  const ps = Float32Array.from(pi), ts = Float32Array.from(surfaceT), wind = Float32Array.from(surfaceWind);
  const diagnostics = model.diagnostics();
  self.postMessage({ type: 'frame', frame: frame++, time, day: time / 86400, ps, ts, wind, diagnostics }, [ps.buffer, ts.buffer, wind.buffer]);
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
    surfaceWind = new Float64Array(model.mesh.nCells);
    vector = new Float64Array(3 * model.mesh.nCells);
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
