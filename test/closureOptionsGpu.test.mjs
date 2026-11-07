import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { createSigmaCore, sigmaInterfaces, DIVERGENCE_DAMPING, DIVERGENCE_SURFACES, TEMPERATURE_CLOSURES } from '../js/dynamics/sigmaCore.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuCore } = gpuAvailable ? await import('../js/gpu/core.gpu.js') : {};
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};
const { readRanges } = gpuAvailable ? await import('../js/gpu/device.module.js') : {};
const skip = !gpuAvailable && 'webgpu not installed';

const N = 16, dt = 1350, levels = sigmaInterfaces('bl48'), K = levels.length - 1;
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const model = createModel(new Grid(N), { topography, levels, physics: false });
const { mesh } = model, { spacing } = model.core, C = mesh.nCells, E = mesh.nEdges;
const start = initializeState(model, {});
const nu4 = (spacing / Math.PI) ** 4 / (3 * 3600);
const closureOnly = { levels, nu4, nu4Theta: nu4, spongeRates: null, spongeMeanRates: null, surfaceGeopotential: model.surfaceGeopotential };
const apart = (a, b) => { let diff = 0, scale = 0; for (let n = 0; n < b.length; n++) { diff = Math.max(diff, Math.abs(a[n] - b[n])); scale = Math.max(scale, Math.abs(b[n])); } return { diff, scale, relative: diff / scale }; };
const minus = (a, b) => Float64Array.from(a, (x, n) => x - b[n]);
const noise = (n) => { const x = Math.sin(12.9898 * n + 78.233) * 43758.5453; return x - Math.floor(x) - 0.5; };
const said = ({ diff, scale, relative }) => `${diff.toExponential(2)} of ${scale.toExponential(2)} (${relative.toExponential(2)})`;

async function closureStep(options, initial) {
  const cpu = createSigmaCore(mesh, { ...closureOnly, ...options, splitClosure: true });
  const state = initial.map((a) => Float64Array.from(a));
  cpu.diagnose(state[0], state[1], state[4], state[5]);
  const exner = Float64Array.from(cpu.diagnostics.exnerLayer);
  cpu.phaseClosure(state, 0, K, dt);
  const gpu = await createGpuCore(mesh, { ...closureOnly, ...options, dragCoefficient: 0, topDragDays: 0, gravityWaves: false });
  gpu.upload(initial);
  await gpu.tendency();
  gpu.closure(dt);
  return { cpu: state, gpu: await gpu.download(), exner };
}

function prepare(m) {
  const init = initializeState(m, {});
  for (let a = 0; a < init.length; a++) m.state[a].set(init[a]);
  m.seaIce.load(m.state[6]);
  if (m.load) m.load();
  m.ocean.initialize(m.state[3], m.state[6]);
  m.land.initialize();
  return m;
}

test(`one closure step with divergenceViscosity on the GPU matches the CPU at N=${N} on bl48 along σ and the hybrid surfaces, and twice c d²/dt damps twice the coefficient's step`, { skip }, async () => {
  const initial = start.map((a) => Float64Array.from(a)), equivalent = DIVERGENCE_DAMPING * spacing * spacing / dt;
  for (let k = 0; k < K; k++) for (let e = 0; e < E; e++) initial[2][k * E + e] += 2 * (noise(k * C + mesh.cellsOnEdge[2 * e + 1]) - noise(k * C + mesh.cellsOnEdge[2 * e]));
  const wind32 = Float64Array.from(initial[2], Math.fround);
  for (const divergenceSurfaces of DIVERGENCE_SURFACES) {
    const none = await closureStep({ divergenceSurfaces }, initial);
    const coefficient = await closureStep({ divergenceSurfaces, divergenceDamping: DIVERGENCE_DAMPING }, initial);
    const viscosity = await closureStep({ divergenceSurfaces, divergenceDamping: DIVERGENCE_DAMPING, divergenceViscosity: 2 * equivalent }, initial);
    const damping = (run, engine) => minus(run[engine][2], none[engine][2]);
    const twice = apart(damping(viscosity, 'gpu'), Float64Array.from(damping(coefficient, 'gpu'), (x) => 2 * x));
    const engines = apart(damping(viscosity, 'gpu'), damping(viscosity, 'cpu'));
    const step = apart(minus(viscosity.gpu[2], wind32), minus(viscosity.cpu[2], initial[2]));
    console.log(`${divergenceSurfaces}: ν_d ${(2 * equivalent).toExponential(3)} m²/s; on the GPU its damping against twice the coefficient's ${said(twice)} m/s; GPU against CPU, the damping ${said(engines)}, the step's wind change ${said(step)}`);
    assert.ok(engines.scale > 0.05, `the damping moves the wind ${engines.scale} m/s at most`);
    assert.ok(twice.relative < 1e-3, `${divergenceSurfaces}: twice the step ${twice.diff} m/s against ${twice.scale}`);
    assert.ok(engines.relative < 1e-2, `${divergenceSurfaces}: damping ${engines.diff} m/s against ${engines.scale}`);
    assert.ok(step.relative < 1e-3, `${divergenceSurfaces}: wind change ${step.diff} m/s against ${step.scale}`);
  }
});

test(`one closure step with temperatureClosure 'standard' on the GPU matches the CPU at N=${N} on bl48 over the real terrain`, { skip }, async () => {
  const runs = {};
  for (const temperatureClosure of TEMPERATURE_CLOSURES) runs[temperatureClosure] = await closureStep({ divergenceDamping: DIVERGENCE_DAMPING, temperatureClosure }, start);
  const { exner } = runs.sigma, theta32 = Float64Array.from(start[1], Math.fround);
  const kelvin = (a) => Float64Array.from(a, (x, n) => x * exner[n]);
  const shift = apart(kelvin(minus(runs.standard.gpu[1], runs.sigma.gpu[1])), kelvin(minus(runs.standard.cpu[1], runs.sigma.cpu[1])));
  const change = Object.fromEntries(TEMPERATURE_CLOSURES.map((t) => [t, apart(kelvin(minus(runs[t].gpu[1], theta32)), kelvin(minus(runs[t].cpu[1], start[1])))]));
  console.log(`GPU against CPU in K: the standard reference's shift ${said(shift)}; the step's change along σ ${said(change.sigma)}, less the standard atmosphere ${said(change.standard)}`);
  assert.ok(shift.scale > 1e-3, `the reference moves T ${shift.scale} K at most`);
  assert.ok(shift.relative < 1e-2, `shift ${shift.diff} K against ${shift.scale}`);
  for (const t of TEMPERATURE_CLOSURES) assert.ok(change[t].relative < 1e-2, `${t}: ${change[t].diff} K against ${change[t].scale}`);
  assert.equal(apart(runs.standard.gpu[2], runs.sigma.gpu[2]).diff, 0, 'the wind does not see the temperature closure');
});

test(`the GPU model takes terrain { conserveVolume } and orography { constants: 'cy41r2' }: the CPU model's surface geopotential on the device, and the same blocking, wave drag and stress after a step`, { skip }, async () => {
  const options = { topography, levels: sigmaInterfaces('bl34'), terrain: { conserveVolume: true }, orography: { constants: 'cy41r2' } };
  const cpu = prepare(createModel(new Grid(N), options));
  const gpu = prepare(await createGpuModel(new Grid(N), options));
  assert.equal(cpu.orography.options.waveDrag, 0.3);
  assert.notDeepEqual(Array.from(cpu.surfaceGeopotential), Array.from(model.surfaceGeopotential));
  assert.deepEqual(Array.from(gpu.surfaceGeopotential), Array.from(cpu.surfaceGeopotential));
  const [phis] = await readRanges(gpu.gpu.device, gpu.gpu.buffers.MF, [{ offset: gpu.gpu.layout.MF.PHIS, length: C }]);
  assert.deepEqual(Array.from(phis), Array.from(cpu.surfaceGeopotential, Math.fround));
  cpu.step(dt); await gpu.step(dt); await gpu.settle();
  const PH = gpu.gpu.layout.PH, Kc = cpu.core.K, o = cpu.orography;
  const [beta, wave, stress, launch] = await readRanges(gpu.gpu.device, gpu.gpu.buffers.PH, [['OBETA', Kc * C], ['OWAVE', Kc * C], ['OSTRESS', E], ['OLAUNCH', C]].map(([name, length]) => ({ offset: PH[name], length })));
  const report = { beta: apart(beta, o.beta), wave: apart(wave, o.wave), stress: apart(stress, o.stress), launch: apart(launch, o.launch) };
  console.log(`engines apart with cy41r2 over the volume-kept terrain: ${Object.entries(report).map(([k, v]) => `${k} ${v.relative.toExponential(1)}`).join(', ')}`);
  for (const [k, v] of Object.entries(report)) assert.ok(v.scale > 0 && v.relative < 1e-3, `${k} ${v.relative} of ${v.scale}`);
});
