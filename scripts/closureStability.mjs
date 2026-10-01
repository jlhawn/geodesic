// The growth of the ∇⁴ closure's grid-scale modes over one ocean step in a
// saved state, by class, for each treatment of the token edges: power
// iteration on one RK4 step of du/dt = −ν₄ ∇⁴ F(u) over a class's thick edges,
// F filling its token edges as the ocean's closureVelocity does and the token
// edges' own velocity held at zero, as a velocity carried from the layer above
// is held through the step. A growth above 1 is a mode the step amplifies.
// Treatments: 'tokens' as they are, 'beside' and 'interior' the fills of
// closureVelocity, 'extended' the interior fill with closureAdjoint.
//   [OCEAN='{"closureHours":12}'] [ITERATIONS=60] node scripts/closureStability.mjs <state.bin>
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { topographyFromInt16, createGeography } from '../js/geography.module.js';
import { decodeState } from '../js/stateFile.module.js';
import { laplacianVelocity } from '../js/dynamics/operators.module.js';
import { createOcean, closureCoefficient, closureVelocity, closureAdjoint, THIN } from '../js/ocean/layered.module.js';

const FILE = process.argv[2];
if (!FILE) throw new Error('usage: node scripts/closureStability.mjs <state.bin>');
const OCEAN = { everySteps: 8, ...JSON.parse(process.env.OCEAN ?? '{}') }, ITERATIONS = Number(process.env.ITERATIONS ?? 60);
const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
const N = saved.N, mesh = buildMesh(new Grid(N)), C = mesh.nCells, E = mesh.nEdges;
const geography = createGeography(mesh, topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer));
const ocean = createOcean(mesh, { geography, ...OCEAN });
ocean.load(saved.ocean, Float64Array.from(saved.surfaceT), Float64Array.from(saved.ice));
const dt = OCEAN.everySteps * 1350 * 16 / N;
new Float64Array(ocean.shared.params)[1] = dt;
ocean.tendency(ocean.state, ocean.stages[0]);
const hEdge = new Float64Array(ocean.shared.hEdge), deepest = new Float64Array(ocean.shared.deepestEdge), { edgeOcean } = ocean, L = ocean.layers, rho = ocean.densities;
let spacing = 0;
for (let e = 0; e < E; e++) spacing += mesh.dcEdge[e];
const nu4 = closureCoefficient(spacing / E, OCEAN.closureHours ?? 12);
const MODES = (process.env.MODES ?? 'tokens,beside,interior,extended').split(',');
const filled = new Float64Array(E), lap = new Float64Array(E), lap2 = new Float64Array(E), divS = new Float64Array(C), curlS = new Float64Array(mesh.nVertices);
const rings = { deepest, k: 0, valid: new Uint8Array(E), second: new Float64Array(E) };
function rate(k, mode, v, out) {
  const he = hEdge.subarray(k * E, (k + 1) * E);
  const input = mode === 'tokens' ? v : closureVelocity(mesh, v, he, edgeOcean, OCEAN.closureFill ?? 1, filled, mode === 'interior' || mode === 'extended' ? Object.assign(rings, { k }) : null);
  laplacianVelocity(mesh, input, lap, divS, curlS);
  laplacianVelocity(mesh, lap, lap2, divS, curlS);
  if (mode === 'extended') closureAdjoint(mesh, lap2, OCEAN.closureFill ?? 1, rings);
  for (let e = 0; e < E; e++) out[e] = edgeOcean[e] && he[e] >= THIN ? -nu4 * lap2[e] : 0;
}
const s = [0, 1, 2, 3].map(() => new Float64Array(E)), trial = new Float64Array(E), next = new Float64Array(E);
function step(k, mode, v) {
  rate(k, mode, v, s[0]);
  for (let e = 0; e < E; e++) trial[e] = v[e] + 0.5 * dt * s[0][e];
  rate(k, mode, trial, s[1]);
  for (let e = 0; e < E; e++) trial[e] = v[e] + 0.5 * dt * s[1][e];
  rate(k, mode, trial, s[2]);
  for (let e = 0; e < E; e++) trial[e] = v[e] + dt * s[2][e];
  rate(k, mode, trial, s[3]);
  for (let e = 0; e < E; e++) next[e] = v[e] + dt / 6 * (s[0][e] + 2 * s[1][e] + 2 * s[2][e] + s[3][e]);
  return next;
}
console.log(`== ${FILE.split('/').pop()}: day ${saved.day}, N=${N}, ν₄ ${nu4.toExponential(2)} m⁴/s, step ${dt} s; growth of the fastest mode over one step after ${ITERATIONS} iterations, and where it lies`);
const classes = (process.env.CLASSES ?? '1022,1022.5,1023,1023.5,1024,1024.5,1025,1026,1026.5,1026.8,1026.9').split(',').map(Number).map((r) => rho.findIndex((x) => Math.abs(x - r) < 1e-6)).filter((k) => k > 0);
for (const k of classes) {
  const parts = [];
  for (const mode of MODES) {
    let v = Float64Array.from({ length: E }, (_, e) => (edgeOcean[e] && hEdge[k * E + e] >= THIN ? Math.sin(12.9898 * e) : 0));
    let growth = 0;
    for (let n = 0; n < ITERATIONS; n++) {
      const norm0 = Math.sqrt(v.reduce((a, x) => a + x * x, 0));
      const w = step(k, mode, v);
      const norm1 = Math.sqrt(w.reduce((a, x) => a + x * x, 0));
      growth = norm1 / norm0;
      v = Float64Array.from(w, (x) => x / norm1);
    }
    let at = 0;
    for (let e = 0; e < E; e++) if (Math.abs(v[e]) > Math.abs(v[at])) at = e;
    parts.push(`${mode} ${growth.toFixed(3)} at ${(mesh.latEdge[at] * 180 / Math.PI).toFixed(1)},${(Math.atan2(mesh.xEdge[3 * at + 1], mesh.xEdge[3 * at]) * 180 / Math.PI).toFixed(1)}`);
  }
  console.log(`  ${rho[k].toFixed(2)}: ${parts.join('; ')}`);
}
