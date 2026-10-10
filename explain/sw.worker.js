import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createShallowWater } from '../js/dynamics/shallowWater.module.js';
import { createRK4 } from '../js/dynamics/integrators.module.js';
import { curl } from '../js/dynamics/operators.module.js';
import { EARTH, edgeNormalVelocity, cellField, hyperdiffusion, bump, galewsky, haurwitz, cellVelocity } from './swCases.module.js';

const CASES = { bump, galewsky, haurwitz };
const BUDGET_MS = 40;
let sim = null;

function arrowCells(mesh, coarse) {
  const picks = [];
  for (const cell of new Grid(coarse, { relax: 0 })) {
    const c = cell.centerVertex;
    let best = -1, bestDot = -2;
    for (let i = 0; i < mesh.nCells; i++) { const d = c.x * mesh.xCell[3 * i] + c.y * mesh.xCell[3 * i + 1] + c.z * mesh.xCell[3 * i + 2]; if (d > bestDot) { bestDot = d; best = i; } }
    picks.push(best);
  }
  return Int32Array.from(picks);
}

function start({ N, kind, options = {}, spin = 1, closureHours = 0, coarse = 8, run = 0 }) {
  const mesh = buildMesh(new Grid(N), { radius: EARTH.a, omega: EARTH.omega * spin });
  const model = createShallowWater(mesh, { g: EARTH.g, nu4: closureHours ? hyperdiffusion(mesh, closureHours) : 0 });
  const step = createRK4(mesh.nCells, mesh.nEdges);
  const setup = CASES[kind](options);
  const h = cellField(mesh, setup.height), u = edgeNormalVelocity(mesh, setup.wind);
  const picks = arrowCells(mesh, coarse);
  sim = { mesh, model, step, h, u, mean: setup.mean, dt: 240 * 32 / N, time: 0, owed: 0, picks, kind, run, initialMass: model.diagnostics(h, u).mass };
  send();
}

function send() {
  const { mesh, h, u, mean, picks, kind, time } = sim;
  const field = new Float32Array(mesh.nCells);
  if (kind === 'galewsky') {
    const zeta = curl(mesh, u);
    for (let i = 0; i < mesh.nCells; i++) { let s = 0; const n = mesh.nEdgesOnCell[i]; for (let m = 0; m < n; m++) s += zeta[mesh.verticesOnCell[mesh.maxEdges * i + m]]; field[i] = s / n; }
  } else for (let i = 0; i < mesh.nCells; i++) field[i] = h[i] - mean;
  const arrows = new Float32Array(3 * picks.length);
  for (let k = 0; k < picks.length; k++) cellVelocity(mesh, u, picks[k], arrows, k);
  let phase = null;
  if (kind === 'haurwitz') {
    let cs = 0, sn = 0;
    for (let i = 0; i < mesh.nCells; i++) {
      const z = mesh.xCell[3 * i + 2], lat = Math.asin(Math.max(-1, Math.min(1, z)));
      if (Math.abs(lat) > 0.6) continue;
      const lon = Math.atan2(mesh.xCell[3 * i + 1], mesh.xCell[3 * i]), w = mesh.areaCell[i] * Math.cos(lat) ** 4;
      cs += w * h[i] * Math.cos(4 * lon); sn += w * h[i] * Math.sin(4 * lon);
    }
    phase = Math.atan2(sn, cs) / 4;
  }
  const mass = sim.model.diagnostics(h, u).mass;
  self.postMessage({ type: 'frame', run: sim.run, time, field, arrows, picks, phase, massDrift: (mass - sim.initialMass) / sim.initialMass }, [field.buffer, arrows.buffer]);
}

function advance(seconds) {
  sim.owed += seconds;
  const t0 = performance.now();
  while (sim.owed >= sim.dt && performance.now() - t0 < BUDGET_MS) { sim.step(sim.model.tendency, sim.h, sim.u, sim.dt); sim.owed -= sim.dt; sim.time += sim.dt; }
  if (sim.owed > sim.dt) sim.owed = sim.dt;
  send();
}

self.onmessage = ({ data }) => {
  if (data.type === 'start') start(data);
  else if (data.type === 'advance' && sim) advance(data.seconds);
};
