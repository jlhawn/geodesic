// node gpuIdentity.mjs <repo>: hashes of (1) the coupled N=6 model after stepBatch(33) as test/stepBatch.test.mjs
// prepares it, and (2) the GPU ocean of test/layeredGpu.test.mjs's N=8 scenario after 1 + 20 advances.
import { createHash } from 'node:crypto';
const REPO = process.argv[2].replace(/\/$/, '');
const { create, globals } = await import('/Users/jlhawn/git_repos/jlhawn/geodesic/node_modules/webgpu/index.js');
Object.assign(globalThis, globals);
Object.defineProperty(globalThis, 'navigator', { value: { gpu: create([]) }, configurable: true });
const js = (p) => import(`${REPO}/js/${p}`);
const { Grid } = await js('grid.module.js');
const { initializeState } = await js('physics/init.module.js');
const { syntheticTopography } = await js('geography.module.js');
const { createGpuModel } = await js('gpu/model.gpu.js');
const { readBuffer } = await js('gpu/device.module.js');
const { createModel } = await js('model.module.js');
const OCEAN = JSON.parse(process.env.OCEAN ?? '{}');
const hashOf = async (device, buffer) => createHash('sha256').update(Buffer.from((await readBuffer(device, buffer, buffer.size, Uint32Array)).buffer)).digest('hex').slice(0, 12);

{
  const topography = syntheticTopography(90, 180, (lat, lon) => ((Math.cos(lon) > 0 && Math.abs(lat) < 1.2) || lat < -1.15 ? 300 + 2500 * Math.exp(-(((lat - 0.3) / 0.3) ** 2)) : -4000));
  const model = await createGpuModel(new Grid(6), { topography, ocean: OCEAN, land: { growthTime: 3 * 3600, declineTime: 2 * 3600, snowDeclineTime: 4 * 3600 } });
  const init = initializeState(model, {});
  for (let a = 0; a < init.length; a++) model.state[a].set(init[a]);
  for (let i = 0; i < model.mesh.nCells; i++) if (model.geography.land[i]) model.state[6][i] = 0;
  model.load();
  model.ocean.initialize(model.state[3], model.state[6]);
  model.land.initialize();
  const C = model.mesh.nCells, soil = new Float64Array(C), snow = new Float64Array(C);
  for (let i = 0; i < C; i++) if (model.geography.land[i]) { soil[i] = 40; snow[i] = model.mesh.latCell[i] > 1.0 ? 5 : 0; }
  model.land.load({ soil, snow, vegetation: new Float64Array(C).fill(0.5) });
  await model.stepBatch(33, 900);
  await model.settle();
  const { device, buffers } = model.gpu, ocean = model.oceanEngine.buffers;
  const named = { S: buffers.S, PH: buffers.PH, D: buffers.D, OS: ocean.S, OD: ocean.OD };
  const out = [];
  for (const [name, buffer] of Object.entries(named)) out.push(`${name} ${await hashOf(device, buffer)}`);
  console.log(`coupled N=6, stepBatch(33): ${out.join(', ')}`);
  model.destroy();
}
{
  const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 300 : -4000));
  const cpuModel = createModel(new Grid(8), { topography });
  const init = initializeState(cpuModel, {});
  for (let a = 0; a < init.length; a++) cpuModel.state[a].set(init[a]);
  for (let i = 0; i < cpuModel.mesh.nCells; i++) if (cpuModel.geography.land[i]) cpuModel.state[6][i] = 0;
  const mesh = cpuModel.mesh, surfaceT0 = Float64Array.from(cpuModel.state[3]), ice = Float64Array.from(cpuModel.state[6]);
  const stress = new Float64Array(mesh.nEdges);
  for (let e = 0; e < mesh.nEdges; e++) { const lon = Math.atan2(mesh.xEdge[3 * e + 1], mesh.xEdge[3 * e]); stress[e] = 0.08 * Math.sin(2 * mesh.latEdge[e]) * Math.cos(lon) + 0.01 * Math.sin(5 * lon); }
  const gpuModel = await createGpuModel(new Grid(8), { topography, ocean: { everySteps: 1, ...OCEAN } });
  const ocean = gpuModel.oceanEngine;
  ocean.initialize(surfaceT0, ice);
  for (let n = 0; n < 21; n++) await ocean.advance(Float64Array.from(surfaceT0), ice, stress, 1350);
  const { device } = gpuModel.gpu;
  console.log(`ocean N=8, 21 advances: OS ${await hashOf(device, ocean.buffers.S)}, OD ${await hashOf(device, ocean.buffers.OD)}, PH ${await hashOf(device, gpuModel.gpu.buffers.PH)}`);
  gpuModel.destroy();
}
process.exit(0);
