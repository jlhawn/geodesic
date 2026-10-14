// node cpu20.mjs <repo> [real]  -> hash of the CPU ocean state after 20 steps at N=8
import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
const REPO = process.argv[2], REAL = process.argv[3] === 'real';
const js = (p) => import(`${REPO}/js/${p}`);
const { Grid } = await js('grid.module.js');
const { createModel } = await js('model.module.js');
const { createOcean } = await js('ocean/layered.module.js');
const { initializeState } = await js('physics/init.module.js');
const { syntheticTopography, topographyFromInt16 } = await js('geography.module.js');
const topography = REAL ? topographyFromInt16(readFileSync(`${REPO}/data/topography_0p25.bin`).buffer) : syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 300 : -4000));
const options = JSON.parse(process.env.OCEAN ?? '{}');
const cpuModel = createModel(new Grid(8), { topography });
const init = initializeState(cpuModel, {});
for (let a = 0; a < init.length; a++) cpuModel.state[a].set(init[a]);
for (let i = 0; i < cpuModel.mesh.nCells; i++) if (cpuModel.geography.land[i]) cpuModel.state[6][i] = 0;
const mesh = cpuModel.mesh, surfaceT0 = Float64Array.from(cpuModel.state[3]), ice = Float64Array.from(cpuModel.state[6]);
const stress = new Float64Array(mesh.nEdges);
for (let e = 0; e < mesh.nEdges; e++) { const lon = Math.atan2(mesh.xEdge[3 * e + 1], mesh.xEdge[3 * e]); stress[e] = 0.08 * Math.sin(2 * mesh.latEdge[e]) * Math.cos(lon) + 0.01 * Math.sin(5 * lon); }
const ocean = createOcean(mesh, { geography: cpuModel.geography, everySteps: 1, ...options });
ocean.initialize(Float64Array.from(surfaceT0), ice);
const flux = new Float64Array(mesh.nCells);
let thin = 0;
for (let n = 0; n < 20; n++) {
  ocean.advance(Float64Array.from(surfaceT0), ice, flux, stress, 1350);
  const hEdge = new Float64Array(ocean.shared.hEdge);
  for (let e = 0; e < mesh.nEdges; e++) if (ocean.edgeOcean[e] && hEdge[e] < 50) thin++;
}
const hash = createHash('sha256');
for (const a of [ocean.h, ocean.u, ocean.Q, ocean.W, ocean.eta]) hash.update(Buffer.from(a.buffer, a.byteOffset, a.byteLength));
console.log(`${REAL ? 'real' : 'synthetic'} ${hash.digest('hex').slice(0, 16)} edge-steps with h0e<50 after a step: ${thin}`);
