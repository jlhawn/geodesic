// The ∇⁴ closure's work on each class's own flow in a saved state, as a
// damping rate −Σ w u·(−ν₄ ∇⁴ F(u)) / Σ w u² over the class's thick edges
// in a box (w = dc·dv; positive damps), for the token edges as they are
// ('tokens', carrying the layer above), the ocean's closure ('default': the
// interior fill in two rings with closureAdjoint) and the same fill on the
// first ring alone with its transpose ('one ring').
//   [BOX=-2,2,180,260] [CLASSES=1022.5,1023,...] node scripts/closureDrag.mjs <state.bin>...
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { topographyFromInt16, createGeography } from '../js/geography.module.js';
import { decodeState } from '../js/stateFile.module.js';
import { laplacianVelocity } from '../js/dynamics/operators.module.js';
import { createOcean, closureCoefficient, closureVelocity, closureAdjoint, THIN } from '../js/ocean/layered.module.js';

const FILES = process.argv.slice(2);
if (!FILES.length) throw new Error('usage: node scripts/closureDrag.mjs <state.bin>...');
const [south, north, west, east] = (process.env.BOX ?? '-2,2,180,260').split(',').map(Number);
const CLASSES = (process.env.CLASSES ?? '1022.25,1022.5,1022.75,1023,1023.25,1023.5,1023.75,1024,1024.5').split(',').map(Number);
const deg = 180 / Math.PI;
let grid = null;
for (const FILE of FILES) {
  const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
  const N = saved.N;
  if (!grid || grid.N !== N) {
    const mesh = buildMesh(new Grid(N));
    grid = { N, mesh, geography: createGeography(mesh, topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer)) };
  }
  const { mesh, geography } = grid, E = mesh.nEdges, { dcEdge, dvEdge, latEdge, xEdge } = mesh;
  const ocean = createOcean(mesh, { geography, everySteps: 8 });
  ocean.load(saved.ocean, Float64Array.from(saved.surfaceT), Float64Array.from(saved.ice));
  new Float64Array(ocean.shared.params)[1] = 8 * 1350 * 16 / N;
  ocean.tendency(ocean.state, ocean.stages[0]);
  const hEdge = new Float64Array(ocean.shared.hEdge), deepest = new Float64Array(ocean.shared.deepestEdge), { edgeOcean, u } = ocean, rho = ocean.densities;
  let spacing = 0;
  for (let e = 0; e < E; e++) spacing += dcEdge[e];
  const nu4 = closureCoefficient(spacing / E, 12);
  const inBox = Uint8Array.from({ length: E }, (_, e) => {
    const lat = latEdge[e] * deg, lon = ((Math.atan2(xEdge[3 * e + 1], xEdge[3 * e]) * deg) % 360 + 360) % 360;
    return lat >= south && lat <= north && lon >= west && lon < east ? 1 : 0;
  });
  const rings = { deepest, k: 0, valid: new Uint8Array(E), second: new Float64Array(E) }, filled = new Float64Array(E);
  function closure(k, v, mode) {
    let input = v;
    if (mode !== 'tokens') {
      rings.k = k;
      input = closureVelocity(mesh, v, hEdge.subarray(k * E, (k + 1) * E), edgeOcean, 1, filled, rings);
      if (mode === 'one ring') for (let e = 0; e < E; e++) if (rings.valid[e] === 3) { rings.valid[e] = 0; input[e] = v[e]; }
    }
    const lap2 = laplacianVelocity(mesh, laplacianVelocity(mesh, input));
    if (mode !== 'tokens') closureAdjoint(mesh, lap2, 1, rings);
    return lap2;
  }
  console.log(`== ${FILE.split('/').pop()}: day ${saved.day}, N=${N}; the closure's damping rate of each class's own flow over its thick edges in ${south}..${north} lat, ${west}..${east} E (1e-6 1/s, positive damps)`);
  for (const r of CLASSES) {
    const k = rho.findIndex((x) => Math.abs(x - r) < 1e-6);
    if (k <= 0) continue;
    const v = u.slice(k * E, (k + 1) * E), thick = (e) => inBox[e] && edgeOcean[e] && hEdge[k * E + e] >= THIN;
    let count = 0;
    for (let e = 0; e < E; e++) if (thick(e)) count++;
    const parts = ['tokens', 'default', 'one ring'].map((mode) => {
      const lap2 = closure(k, v, mode);
      let work = 0, energy = 0;
      for (let e = 0; e < E; e++) if (thick(e)) { const w = dcEdge[e] * dvEdge[e]; work += w * v[e] * nu4 * lap2[e]; energy += w * v[e] * v[e]; }
      return `${mode} ${(1e6 * work / energy).toFixed(2)}`;
    });
    console.log(`  ${r.toFixed(2)} (${count} edges): ${parts.join(', ')}`);
  }
}
