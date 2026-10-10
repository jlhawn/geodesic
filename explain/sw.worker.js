import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createShallowWater } from '../js/dynamics/shallowWater.module.js';
import { createRK4 } from '../js/dynamics/integrators.module.js';
import { curl } from '../js/dynamics/operators.module.js';
import { EARTH, edgeNormalVelocity, cellField, hyperdiffusion, bump, galewsky } from './swCases.module.js';

const CASES = { bump, galewsky };
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

function cellVelocity(mesh, u, i, out, k) {
  const { maxEdges, nEdgesOnCell, edgesOnCell, nEdge, xCell } = mesh;
  const x = xCell[3 * i], y = xCell[3 * i + 1], z = xCell[3 * i + 2];
  const rho = Math.hypot(x, y) || 1e-12, east = [-y / rho, x / rho, 0], north = [-z * x / rho, -z * y / rho, rho];
  let saa = 0, sab = 0, sbb = 0, sau = 0, sbu = 0;
  for (let m = 0; m < nEdgesOnCell[i]; m++) {
    const e = edgesOnCell[maxEdges * i + m], n = [nEdge[3 * e], nEdge[3 * e + 1], nEdge[3 * e + 2]];
    const a = n[0] * east[0] + n[1] * east[1] + n[2] * east[2], b = n[0] * north[0] + n[1] * north[1] + n[2] * north[2];
    saa += a * a; sab += a * b; sbb += b * b; sau += a * u[e]; sbu += b * u[e];
  }
  const det = saa * sbb - sab * sab, ue = (sau * sbb - sbu * sab) / det, un = (saa * sbu - sab * sau) / det;
  for (let c = 0; c < 3; c++) out[3 * k + c] = ue * east[c] + un * north[c];
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
  const mass = sim.model.diagnostics(h, u).mass;
  self.postMessage({ type: 'frame', run: sim.run, time, field, arrows, picks, massDrift: (mass - sim.initialMass) / sim.initialMass }, [field.buffer, arrows.buffer]);
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
