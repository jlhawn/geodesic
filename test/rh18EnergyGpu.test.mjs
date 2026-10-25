import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createOcean as createCpuOcean } from '../js/ocean/layered.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';

/*
 * The GPU hybrid ocean's boundary layer with mstarScheme 'rh18'
 * (js/gpu/layeredOcean.gpu.js: oReadSurface keeps the heat and salt the
 * surface gave the step, oColumnCells forms the surface buoyancy flux from
 * them, blEnergy takes Reichl and Hallberg's m*N C + m*S) against the CPU's
 * at N=16, under a surface that warms the north, cools the south and
 * grows ice at the south pole. One GPU model at a time (CLAUDE.md).
 */
let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const topography = syntheticTopography(90, 180, (lat, lon) => (Math.cos(lon) > 0 && Math.abs(lat) < 1.2 ? 300 : -4000));
const HYBRID = { vertical: 'hybrid', everySteps: 1, iceDrift: false };
const N = 16, DT = 1350;

function stats(cpu, gpu, mask = null) {
  let maxDiff = 0, at = -1, sumSq = 0, sumRef = 0;
  for (let x = 0; x < cpu.length; x++) {
    if (mask && !mask[x]) continue;
    const d = Math.abs(cpu[x] - gpu[x]);
    if (d > maxDiff) { maxDiff = d; at = x; }
    sumSq += d * d; sumRef += cpu[x] * cpu[x];
  }
  return { maxDiff, at, rmsRel: Math.sqrt(sumSq / Math.max(sumRef, 1e-300)) };
}
function scenario(n) {
  const cpuModel = createModel(new Grid(n), { topography });
  const init = initializeState(cpuModel, {});
  for (let a = 0; a < init.length; a++) cpuModel.state[a].set(init[a]);
  for (let i = 0; i < cpuModel.mesh.nCells; i++) if (cpuModel.geography.land[i]) cpuModel.state[6][i] = 0;
  const mesh = cpuModel.mesh, stress = new Float64Array(mesh.nEdges);
  for (let e = 0; e < mesh.nEdges; e++) {
    const lon = Math.atan2(mesh.xEdge[3 * e + 1], mesh.xEdge[3 * e]);
    stress[e] = 0.08 * Math.sin(2 * mesh.latEdge[e]) * Math.cos(lon) + 0.01 * Math.sin(5 * lon);
  }
  return { cpuModel, mesh, surfaceT: Float64Array.from(cpuModel.state[3]), ice: Float64Array.from(cpuModel.state[6]), stress };
}
const cpuFields = (ocean) => ({
  h: ocean.h, u: ocean.u, eta: ocean.eta,
  T: Float64Array.from(ocean.Q, (q, n) => q / Math.max(1e-9, ocean.h[n])), S: Float64Array.from(ocean.W, (w, n) => w / Math.max(1e-9, ocean.h[n])),
});

async function parity(options, reference, moved) {
  const { cpuModel, mesh, surfaceT, ice, stress } = scenario(N), C = mesh.nCells, chosen = { ...HYBRID, ...options };
  const cpu = createCpuOcean(mesh, { geography: cpuModel.geography, ...chosen }), other = createCpuOcean(mesh, { geography: cpuModel.geography, ...HYBRID, ...reference });
  const gpuModel = await createGpuModel(new Grid(N), { topography, ocean: chosen });
  const gpu = gpuModel.oceanEngine, flux = new Float64Array(C);
  for (const ocean of [cpu, other, gpu]) ocean.initialize(Float64Array.from(surfaceT), ice);
  const forced = (n) => Float64Array.from(surfaceT, (t, i) => t + (mesh.latCell[i] > 0 ? 0.01 : -0.01) * (n + 1));
  const grown = (n) => Float64Array.from(ice, (v, i) => (mesh.latCell[i] < -1.2 ? 0.02 * (n + 1) : v));
  for (let n = 0; n < 12; n++) {
    cpu.advance(forced(n), grown(n), flux, stress, DT);
    other.advance(forced(n), grown(n), flux, stress, DT);
    await gpu.advance(forced(n), grown(n), stress, DT);
  }
  const final = await gpu.download(), { boundaryDepths } = await gpu.restartArrays(), wet = Array.from(cpu.h, (v) => v > 1), a = cpuFields(cpu), limit = { h: 1e-2, u: 2e-2, T: 1e-3, S: 1e-4, eta: 1e-2 };
  for (const f of ['h', 'u', 'T', 'S', 'eta']) {
    const r = stats(a[f], final[f], f === 'T' || f === 'S' ? wet : null);
    console.log(`${f} rms ${r.rmsRel.toExponential(2)} max ${r.maxDiff.toExponential(2)}`);
    assert.ok(r.rmsRel < limit[f], `${f} rms relative diff ${r.rmsRel} (largest ${r.maxDiff} at ${r.at})`);
  }
  let engines = 0, schemes = 0;
  for (let i = 0; i < C; i++) if (wet[i]) { engines += Math.abs(boundaryDepths[i] - cpu.boundaryDepths[i]); schemes += Math.abs(other.boundaryDepths[i] - cpu.boundaryDepths[i]); }
  console.log(`boundary depths: GPU against CPU ${engines.toFixed(2)} m summed, ${JSON.stringify(options)} against ${JSON.stringify(reference)} on the CPU ${schemes.toFixed(2)} m`);
  assert.ok(schemes > moved, `the scheme moves the boundary layer: ${schemes} m`);
  assert.ok(engines < 0.1 * schemes, `the engines part by ${engines} m, the two schemes by ${schemes} m`);
  gpuModel.destroy();
}

test('with mstarScheme rh18 the GPU hybrid ocean forms the surface buoyancy flux from the heat the surface gave it and the brine of the ice grown, and its boundary layer tracks the CPU\'s over twelve steps at N=16', { skip: !gpuAvailable && 'webgpu not installed' }, () => parity({ mstarScheme: 'rh18' }, {}, 10));

test('with langmuirScheme rescale the GPU hybrid ocean multiplies the m* of rh18 by the capped Langmuir enhancement of the wind\'s surface-layer Langmuir number, tapered by the ice, and its boundary layer tracks the CPU\'s over twelve steps at N=16', { skip: !gpuAvailable && 'webgpu not installed' }, () => parity({ mstarScheme: 'rh18', langmuirScheme: 'rescale' }, { mstarScheme: 'rh18' }, 5));
