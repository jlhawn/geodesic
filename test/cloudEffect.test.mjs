import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { CLEAR_SUMMED } from '../js/physics/radiation.module.js';
import { encodeState, decodeState } from '../js/stateFile.module.js';
import { RADIATION_FIELDS, savedRadiationField } from '../js/physics/regrid.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const CLEAR_SLOTS = { clearAbsorbedSolar: 'ABSCLRSUM', clearOutgoingLongwave: 'OLRCLRSUM' };
const EFFECT_SLOTS = { meanShortwaveCloudEffect: 'SWCREMEAN', meanLongwaveCloudEffect: 'LWCREMEAN' };
const DT = 900;
const UNSCATTERED = { rayleighDepth: 0, landAerosol: 0, seaAerosol: 0, skylight: 0.15, upwardAbsorption: false };
const GRAY = { cloudScattering: 95, cloudAbsorption: 130 };
const STEPS = 16;

function stats(cpu, gpu) {
  let maxDiff = 0, sumSq = 0, sumRef = 0;
  for (let x = 0; x < cpu.length; x++) { const d = Math.abs(cpu[x] - gpu[x]); maxDiff = Math.max(maxDiff, d); sumSq += d * d; sumRef += cpu[x] * cpu[x]; }
  return { maxDiff, rmsRel: Math.sqrt(sumSq / Math.max(sumRef, 1e-300)) };
}

test('the clear-sky pass gives each column the top-of-atmosphere fluxes of the same column with its resolved cloud taken away', () => {
  const model = createModel(new Grid(6), { ocean: false, radiation: { clearSkyPass: true } });
  const init = initializeState(model, {});
  init.forEach((values, a) => model.state[a].set(values));
  const { radiation, core } = model, { K, C } = core.diagnostics;
  const [pi, theta, , surfaceT, q] = model.state;
  const qc = new Float64Array(K * C), none = new Float64Array(K * C);
  for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) qc[k * C + i] = k % 3 === 0 ? 2e-4 : 0;
  for (let i = 0; i < C; i++) core.diagnoseColumn(i, pi, theta, q, qc);
  radiation.setTime(0);
  const bottom = (K - 1) * C;
  let checked = 0, lit = 0, worst = 0;
  for (let i = 0; i < C; i++) {
    const args = (water) => [i, pi[i], theta, surfaceT[i], 5, undefined, radiation.insolation(i), q[bottom + i], q, water, 0.07, 0.06, 1, undefined, 0, 0, 0, 0];
    radiation.column(...args(none));
    const clear = { absorbedSolar: radiation.budget.absorbedSolar, outgoingLongwave: radiation.budget.outgoingLongwave };
    assert.equal(radiation.budget.clearAbsorbedSolar, clear.absorbedSolar, `cell ${i}: a cloudless column's clear-sky ASR is its ASR`);
    assert.ok(Math.abs(radiation.budget.clearOutgoingLongwave - clear.outgoingLongwave) <= 1e-12 * clear.outgoingLongwave, `cell ${i}: a cloudless column's clear-sky OLR is its OLR`);
    radiation.column(...args(qc));
    const b = radiation.budget;
    worst = Math.max(worst, Math.abs(b.clearAbsorbedSolar - clear.absorbedSolar), Math.abs(b.clearOutgoingLongwave - clear.outgoingLongwave) / clear.outgoingLongwave);
    assert.ok(Math.abs(b.clearAbsorbedSolar - clear.absorbedSolar) <= 1e-12 * Math.max(1, clear.absorbedSolar), `cell ${i}: the cloudy column's clear-sky ASR ${b.clearAbsorbedSolar} against ${clear.absorbedSolar}`);
    assert.ok(Math.abs(b.clearOutgoingLongwave - clear.outgoingLongwave) <= 1e-12 * clear.outgoingLongwave, `cell ${i}: the cloudy column's clear-sky OLR ${b.clearOutgoingLongwave} against ${clear.outgoingLongwave}`);
    assert.ok(b.outgoingLongwave < b.clearOutgoingLongwave, `cell ${i}: cloud lowers the OLR`);
    if (b.insolation > 0) { lit++; assert.ok(b.absorbedSolar < b.clearAbsorbedSolar, `cell ${i}: cloud lowers the ASR`); }
    checked++;
  }
  console.log(`${checked} columns (${lit} lit) with cloud in every third layer: clear-sky ASR and OLR equal the cloudless column's to ${worst.toExponential(1)}`);
});

async function engines(radiation, moist = {}, cloud = 0) {
  const cpu = createModel(new Grid(6), { ocean: false, radiation, moist }), gpu = await createGpuModel(new Grid(6), { ocean: false, radiation, moist });
  const init = initializeState(cpu, {}), { K, C, sigmaMid } = cpu.core.diagnostics;
  for (let k = 0; k < K; k++) if (sigmaMid[k] > 0.5 && sigmaMid[k] < 0.9) for (let i = 0; i < C; i++) init[5][k * C + i] = cloud;
  for (let a = 0; a < init.length; a++) { cpu.state[a].set(init[a]); gpu.state[a].set(init[a]); }
  gpu.load();
  return { cpu, gpu, C: cpu.mesh.nCells };
}

test('both engines sum the clear-sky fluxes per cell alike and read out the day-mean cloud effects, mirrored per cell and carried in a saved state (under the gray optics and the saturation adjustment, whose trajectories stay together cell by cell)', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { cpu, gpu, C } = await engines({ clearSkyPass: true, cloudOverlap: 'maximumRandom', ...UNSCATTERED, ...GRAY }, { condensation: 'saturation', iceSaturation: false, iceFall: null }, 3e-4);
  for (let s = 0; s < STEPS; s++) { cpu.step(DT); await gpu.step(DT); }
  const device = await gpu.gpu.downloadPhysics();
  const cpuSums = Object.fromEntries(CLEAR_SUMMED.map((name) => [name, Float64Array.from(cpu.radiation.summed[name])]));
  const gpuSums = Object.fromEntries(CLEAR_SUMMED.map((name) => [name, Float64Array.from(device[CLEAR_SLOTS[name]].subarray(0, C))]));
  const parity = Object.fromEntries(CLEAR_SUMMED.map((name) => [name, stats(cpuSums[name], gpuSums[name])]));
  console.log(`${STEPS} steps at N=6, per-cell clear-sky sums CPU against GPU: ${CLEAR_SUMMED.map((name) => `${name} rms ${parity[name].rmsRel.toExponential(1)} (max ${parity[name].maxDiff.toExponential(1)} W/m² summed)`).join(', ')}`);
  for (const name of CLEAR_SUMMED) assert.ok(parity[name].rmsRel < 1e-4, `${name}: per-cell rms ${parity[name].rmsRel}`);

  const dc = cpu.diagnostics(), dg = await gpu.diagnostics();
  await gpu.sync();
  console.log(`day means over the ${STEPS} steps (CPU / GPU): ASR ${dc.absorbedSolar.toFixed(3)} / ${dg.absorbedSolar.toFixed(3)}, clear ${dc.clearAbsorbedSolar.toFixed(3)} / ${dg.clearAbsorbedSolar.toFixed(3)}; OLR ${dc.outgoingLongwave.toFixed(3)} / ${dg.outgoingLongwave.toFixed(3)}, clear ${dc.clearOutgoingLongwave.toFixed(3)} / ${dg.clearOutgoingLongwave.toFixed(3)}; SWCRE ${dc.shortwaveCloudEffect.toFixed(3)} / ${dg.shortwaveCloudEffect.toFixed(3)}, LWCRE ${dc.longwaveCloudEffect.toFixed(3)} / ${dg.longwaveCloudEffect.toFixed(3)} W/m²`);
  for (const d of [dc, dg]) {
    assert.ok(Math.abs(d.shortwaveCloudEffect - (d.absorbedSolar - d.clearAbsorbedSolar)) < 1e-3 && Math.abs(d.longwaveCloudEffect - (d.clearOutgoingLongwave - d.outgoingLongwave)) < 1e-3, 'each effect is the difference of its means');
    assert.ok(d.shortwaveCloudEffect < -1 && d.longwaveCloudEffect > 1, `the run has cloud to see: SWCRE ${d.shortwaveCloudEffect}, LWCRE ${d.longwaveCloudEffect}`);
  }
  for (const name of ['clearAbsorbedSolar', 'clearOutgoingLongwave']) assert.ok(Math.abs(dg[name] - dc[name]) < 1e-4 * dc[name], `${name}: the engines agree`);
  for (const name of ['shortwaveCloudEffect', 'longwaveCloudEffect']) assert.ok(Math.abs(dg[name] - dc[name]) < 0.05, `${name}: the engines agree to 0.05 W/m², ${dc[name]} against ${dg[name]}`);

  const after = await gpu.gpu.downloadPhysics();
  for (let i = 0; i < C; i++) {
    for (const name of CLEAR_SUMMED) assert.ok(cpu.radiation.summed[name][i] === 0 && after[CLEAR_SLOTS[name]][i] === 0, `cell ${i}: ${name} starts again`);
    assert.ok(Math.abs(cpu.radiation.meanShortwaveCloudEffect[i] - (cpu.radiation.meanAbsorbedSolar[i] - cpuSums.clearAbsorbedSolar[i] / STEPS)) <= 1e-9, `cell ${i}: the CPU's per-cell shortwave effect`);
    assert.ok(Math.abs(cpu.radiation.meanLongwaveCloudEffect[i] - (cpuSums.clearOutgoingLongwave[i] / STEPS - cpu.radiation.meanOutgoingLongwave[i])) <= 1e-9, `cell ${i}: the CPU's per-cell longwave effect`);
    for (const [name, slot] of Object.entries(EFFECT_SLOTS)) assert.equal(gpu.radiation[name][i], after[slot][i], `cell ${i}: ${name} mirrored`);
  }
  const effects = Object.fromEntries(Object.keys(EFFECT_SLOTS).map((name) => [name, stats(cpu.radiation[name], gpu.radiation[name])]));
  console.log(`per-cell cloud effects CPU against GPU: ${Object.entries(effects).map(([name, s]) => `${name} rms ${s.rmsRel.toExponential(1)} (max ${s.maxDiff.toExponential(1)} W/m²)`).join(', ')}`);
  for (const s of Object.values(effects)) assert.ok(s.rmsRel < 2e-3 && s.maxDiff < 0.5, `per-cell effects apart by ${s.maxDiff} W/m² (rms ${s.rmsRel})`);

  const saved = await decodeState(encodeState({ N: 6, K: gpu.core.K, day: 0, time: gpu.time, pi: gpu.state[0], ...Object.fromEntries(Object.keys(RADIATION_FIELDS).map((name) => [name, gpu.radiation[name]])) }));
  const fresh = await createGpuModel(new Grid(6), { ocean: false, radiation: { clearSkyPass: true } });
  fresh.state.forEach((values, a) => values.set(gpu.state[a]));
  for (const name of Object.keys(RADIATION_FIELDS)) fresh.radiation[name].set(savedRadiationField(saved, name, fresh));
  fresh.load();
  const placed = await fresh.gpu.downloadPhysics();
  for (let i = 0; i < C; i++) for (const [name, slot] of Object.entries(EFFECT_SLOTS)) assert.equal(placed[slot][i], Math.fround(gpu.radiation[name][i]), `cell ${i}: the saved ${name} on the device`);
  fresh.destroy(); gpu.destroy();
});

test('without the clear-sky pass the read-out has no cloud effects and the clear-sky sums stay empty on both engines', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { cpu, gpu, C } = await engines({});
  for (let s = 0; s < 4; s++) { cpu.step(DT); await gpu.step(DT); }
  const device = await gpu.gpu.downloadPhysics();
  for (let i = 0; i < C; i++) for (const name of CLEAR_SUMMED) assert.ok(cpu.radiation.summed[name][i] === 0 && device[CLEAR_SLOTS[name]][i] === 0, `cell ${i}: ${name}`);
  const dc = cpu.diagnostics(), dg = await gpu.diagnostics();
  for (const d of [dc, dg]) for (const name of ['shortwaveCloudEffect', 'longwaveCloudEffect', 'clearAbsorbedSolar', 'clearOutgoingLongwave']) assert.ok(!(name in d), `no ${name}`);
  gpu.destroy();
});
