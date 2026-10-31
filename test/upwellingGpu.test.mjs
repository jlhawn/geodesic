import { test, after } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createOcean } from '../js/ocean/layered.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { upwellingTimescale, upwellingWeight, SURFACE_UPWELLING_DEPTH } from '../js/frames.module.js';
import { layerVelocities, velocityAt, worstDifference, unevenTopLayers } from './helpers/upwelling.mjs';

/*
 * The GPU ocean's upwelling memory against the CPU's at N=16 on a sea
 * with a continent and a shelf: the fold of one state, the same state on
 * both engines, then both stepped through a few ocean steps with a frame
 * after each, at several depths and in the surface view; and the GPU's
 * memory against the layer weights computed from its own state, folded
 * only when the ocean has stepped and restarted when time goes back or a
 * state uploads. One GPU model at a time (CLAUDE.md).
 */
let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const N = 16, DT = 1350, HYBRID = { vertical: 'hybrid', everySteps: 1, iceDrift: false };
const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 300 + 2500 * Math.exp(-(((lat - 0.3) / 0.3) ** 2)) : -4000));
const DEPTHS = [10, 50, 137, 450, 903, 2500, 4500];

let built = null;
function fixture() {
  built ??= (async () => {
    const cpuModel = createModel(new Grid(N), { topography });
    const init = initializeState(cpuModel, {});
    for (let a = 0; a < init.length; a++) cpuModel.state[a].set(init[a]);
    const { mesh, geography } = cpuModel, C = mesh.nCells, E = mesh.nEdges;
    for (let i = 0; i < C; i++) if (geography.land[i]) cpuModel.state[6][i] = 0;
    const surfaceT = Float64Array.from(cpuModel.state[3]), ice = Float64Array.from(cpuModel.state[6]), stress = new Float64Array(E);
    for (let e = 0; e < E; e++) { const lon = Math.atan2(mesh.xEdge[3 * e + 1], mesh.xEdge[3 * e]); stress[e] = 0.08 * Math.sin(2 * mesh.latEdge[e]) * Math.cos(lon) + 0.01 * Math.sin(5 * lon); }
    const model = await createGpuModel(new Grid(N), { topography, ocean: HYBRID });
    model.oceanEngine.initialize(Float64Array.from(surfaceT), ice);
    return { mesh, geography, C, surfaceT, ice, stress, model, gpu: model.oceanEngine, device: model.gpu.device };
  })();
  return built;
}
after(async () => { if (built) (await built).model.destroy(); });

// A CPU ocean stepped `steps` times, its state sent to the device, and a step that advances both.
async function paired(steps, uneven = false) {
  const { mesh, geography, C, surfaceT, ice, stress, gpu, device } = await fixture();
  const cpu = createOcean(mesh, { geography, ...HYBRID });
  cpu.initialize(Float64Array.from(surfaceT), ice);
  const flux = new Float64Array(C), zero = new Float64Array(C);
  for (let n = 0; n < steps; n++) cpu.advance(Float64Array.from(surfaceT), ice, flux, stress, DT);
  if (uneven) unevenTopLayers(cpu.h, C, cpu.layers, cpu.cellOcean);
  const sync = () => gpu.upload({ ...cpu.serialize(), Q: cpu.Q, W: cpu.W, previousT0: cpu.T0, previousIceVolume: ice, fresh: zero, buoyancy: zero, flux: zero, capacity: cpu.capacity }, Float64Array.from(surfaceT), ice);
  await sync();
  const step = async () => {
    cpu.advance(Float64Array.from(surfaceT), ice, flux, stress, DT);
    await gpu.advance(Float64Array.from(surfaceT), ice, stress, DT);
  };
  return { cpu, gpu, device, sync, step, C, mesh };
}

const frameOf = async (gpu, depth, time, upwellingDepth = depth) => (await gpu.frame({ fields: ['upwelling'], depth, upwellingDepth, time })).fields.upwelling;

test('the GPU fold of one state matches the CPU ocean\'s at every depth, in the surface view and where the floor lies above the depth', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { cpu, gpu } = await paired(3, true);
  const time = 3 * DT;
  cpu.foldUpwelling(time);
  for (const depth of DEPTHS) {
    const expected = cpu.fields(depth).upwelling, got = await frameOf(gpu, depth, time);
    const { max, scale, rms, rmsScale } = worstDifference(expected, got);
    console.log(`state fold at ${depth} m: largest |w| ${scale.toExponential(2)} m/s, largest difference ${max.toExponential(2)}, rms ${rms.toExponential(2)} of ${rmsScale.toExponential(2)}`);
    if (depth < 4000) assert.ok(scale > 1e-7, `the flow through ${depth} m reaches ${scale} m/s`);
    assert.ok(max <= 1e-4 * Math.max(scale, 1e-6), `at ${depth} m the GPU parts from the CPU by ${max} against ${scale}`);
  }
  const surface = worstDifference(cpu.fields(0, SURFACE_UPWELLING_DEPTH).upwelling, await frameOf(gpu, 0, time, SURFACE_UPWELLING_DEPTH));
  assert.ok(surface.scale > 1e-7 && surface.max <= 1e-4 * surface.scale, `the surface view parts by ${surface.max} against ${surface.scale}`);
});

test('stepped through ocean steps with a frame after each, the GPU memory tracks the CPU\'s at several depths', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { cpu, gpu, step } = await paired(2);
  let time = 2 * DT;
  for (let n = 0; n < 5; n++) {
    if (n > 0) { await step(); time += DT; }
    cpu.foldUpwelling(time);
    for (const depth of [50, 450, 903, 2500]) {
      const expected = cpu.fields(depth).upwelling, got = await frameOf(gpu, depth, time);
      const { max, scale, rms, rmsScale } = worstDifference(expected, got);
      console.log(`frame ${n + 1} at ${depth} m: largest |w| ${scale.toExponential(2)} m/s, largest difference ${max.toExponential(2)}, rms ${rms.toExponential(2)} of ${rmsScale.toExponential(2)}`);
      assert.ok(max <= 2e-3 * scale && rms <= 1e-3 * rmsScale, `frame ${n + 1} at ${depth} m: the GPU parts from the CPU by ${max} (rms ${rms}) against ${scale} (rms ${rmsScale})`);
    }
  }
});

test('the GPU memory is the layer-weighted fold of the velocities computed from its own state, folded only when the ocean has stepped and restarted by a backward time or an upload', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { cpu, gpu, device, sync, step, C, mesh } = await paired(3);
  const { cellOcean, sillEdge, edgeOcean } = cpu, L = gpu.layers;
  const sea = (values) => values.filter((_, n) => cellOcean[n % C]);
  const own = async () => {
    const state = await gpu.download(), memory = await gpu.readSection('UPWM');
    return { h: state.h, W: layerVelocities(mesh, L, { h: state.h, u: state.u, eta: state.eta, sillEdge, edgeOcean }, cellOcean), memory: Float64Array.from(memory) };
  };
  const compare = (expected, got, label) => {
    const { max, scale } = worstDifference(sea(expected), sea(got));
    console.log(`${label}: largest difference ${max.toExponential(2)} against ${scale.toExponential(2)}`);
    assert.ok(scale > 1e-7 && max <= 1e-4 * scale, `${label}: ${max} against ${scale}`);
  };
  const start = 1e6, gap = 40 * 86400;
  const first = await frameOf(gpu, 500, start), again = await own();
  compare(again.W, again.memory, 'the first fold replaces the memory');
  assert.deepEqual(await frameOf(gpu, 500, start + 5e5), first, 'a frame with no ocean step since the last fold leaves the memory as it was');

  await step();
  await frameOf(gpu, 500, start + gap);
  const second = await own();
  const folded = Float64Array.from(again.W, (previous, n) => {
    const k = Math.floor(n / C), i = n % C;
    if (!cellOcean[i]) return 0;
    let base = 0;
    for (let j = 0; j <= k; j++) base += second.h[j * C + i];
    return previous + upwellingWeight(gap, 1, upwellingTimescale(base)) * (second.W[n] - previous);
  });
  compare(folded, second.memory, 'after a gap of 40 days each layer takes the weight of its own timescale');
  for (const depth of [20, 450, 2500]) {
    const { max, scale } = worstDifference(velocityAt(C, L, second.h, cellOcean, second.memory, depth), await frameOf(gpu, depth, start + gap));
    assert.ok(max <= 1e-5 * Math.max(scale, 1e-6), `the frame at ${depth} m reads the memory by interpolation: ${max} against ${scale}`);
  }

  await frameOf(gpu, 500, 0);
  const restarted = await own();
  compare(restarted.W, restarted.memory, 'the memory restarts when time goes back');
  await step();
  await frameOf(gpu, 500, DT);
  await sync();
  assert.ok((await own()).memory.every((v) => v === 0), 'an upload clears the memory');
  await frameOf(gpu, 500, 2 * DT);
  const uploaded = await own();
  compare(uploaded.W, uploaded.memory, 'the first fold after an upload replaces the memory');
});
