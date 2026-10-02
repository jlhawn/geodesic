import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { clearLongwave } from '../js/physics/longwave.module.js';
import { BENCHMARK, modelColumn, referenceAt, interfaceAt, MOLAR } from '../scripts/standardAtmospheres.mjs';
import { runColumn } from '../scripts/radiationBenchmark.mjs';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};

const SPECTRAL = { longwaveScheme: 'correlated', solarGases: 'clirad' };
const DT = 900;

test('the longwave g-points follow RRTMG over the standard atmospheres and LBLRTM for doubled CO2', () => {
  const lines = [];
  for (const a of ['TROP', 'MLS', 'MLW', 'SAW']) {
    const column = modelColumn(BENCHMARK.atmospheres[a]), ref = BENCHMARK.rrtmgLongwave[a].levels, r = clearLongwave(column), K = column.T.length;
    const olr = r.up[0] - ref[ref.length - 1].up, dlr = r.down[K] - ref[0].down, net = interfaceAt(column, r.net, 20000) - referenceAt(ref, 20000);
    lines.push(`${a} ${olr.toFixed(2)}/${dlr.toFixed(2)}/${net.toFixed(2)}`);
    for (const x of [olr, dlr, net]) assert.ok(Math.abs(x) < 3, `${a}: OLR, DLR and 200 hPa misses ${olr}, ${dlr}, ${net} W/m2`);
  }
  const mls = modelColumn(BENCHMARK.atmospheres.MLS), K = mls.T.length;
  const withCo2 = (vmr) => ({ ...mls, co2: mls.co2.map((_, k) => vmr * MOLAR.co2 / MOLAR.air * (1 - mls.q[k])), ch4: mls.ch4.map((_, k) => 806e-9 * MOLAR.ch4 / MOLAR.air * (1 - mls.q[k])), n2o: mls.n2o.map((_, k) => 275e-9 * MOLAR.n2o / MOLAR.air * (1 - mls.q[k])) });
  const one = clearLongwave(withCo2(287e-6)), two = clearLongwave(withCo2(574e-6)), ref = BENCHMARK.iacono.co2Doubling.longwave;
  const forcing = { toa: one.net[0] - two.net[0], p20000: interfaceAt(mls, one.net, 20000) - interfaceAt(mls, two.net, 20000), surface: two.down[K] - one.down[K] };
  console.log(`OLR/DLR/net-200-hPa misses against RRTMG: ${lines.join(', ')} W/m2; doubled CO2 TOA ${forcing.toa.toFixed(2)} 200 hPa ${forcing.p20000.toFixed(2)} surface ${forcing.surface.toFixed(2)} W/m2`);
  for (const level of ['toa', 'p20000', 'surface']) assert.ok(Math.abs(forcing[level] / ref[level] - 1) < 0.1, `doubled CO2 at ${level}: ${forcing[level]} against ${ref[level]}`);
});

test('the radiation column carries the g-points: a cloudless column at night has the clear-sky longwave of the same layers', () => {
  const column = modelColumn(BENCHMARK.atmospheres.MLW);
  const r = runColumn(SPECTRAL, column), clear = clearLongwave(column);
  assert.ok(Math.abs(r.olr - clear.up[0]) < 0.2 && Math.abs(r.dlr - clear.down[column.T.length]) < 0.2, `OLR ${r.olr} against ${clear.up[0]}, DLR ${r.dlr} against ${clear.down[column.T.length]}`);
  assert.ok(Math.abs(r.budget.clearOutgoingLongwave - r.olr) <= 1e-12 * r.olr, 'no cloud: the clear-sky OLR is the OLR');
});

test('the solar gases absorb what RRTMG gives within 2 per cent at an overhead sun, and close the column budget', () => {
  for (const a of ['TROP', 'MLS', 'MLW', 'SAW']) {
    const column = modelColumn(BENCHMARK.atmospheres[a]), ref = BENCHMARK.rrtmgShortwave[a].levels, top = ref[ref.length - 1], sfc = ref[0];
    const r = runColumn(SPECTRAL, column, { beam: top.down, albedo: 0.2, solarConstant: top.down });
    const atmosphere = top.net - sfc.net, miss = r.budget.atmosphereSolar / atmosphere - 1;
    assert.ok(Math.abs(miss) < 0.02, `${a}: atmosphere absorbs ${r.budget.atmosphereSolar} against ${atmosphere}`);
    const layers = r.sw.reduce((s, x) => s + x, 0);
    assert.ok(Math.abs(layers - r.budget.atmosphereSolar) < 1e-9 * top.down, `${a}: the layers take the atmosphere's absorption, ${layers} against ${r.budget.atmosphereSolar}`);
    assert.ok(r.budget.oxygenSolar > 0 && r.budget.carbonDioxideSolar > 0, `${a}: O2 and CO2 absorb`);
  }
});

async function engines(radiation) {
  const cpu = createModel(new Grid(6), { ocean: false, radiation }), gpu = await createGpuModel(new Grid(6), { ocean: false, radiation });
  const init = initializeState(cpu, {});
  for (let a = 0; a < init.length; a++) { cpu.state[a].set(init[a]); gpu.state[a].set(init[a]); }
  gpu.load();
  return { cpu, gpu, C: cpu.mesh.nCells };
}

function stats(cpu, gpu) {
  let maxDiff = 0, sumSq = 0, sumRef = 0;
  for (let x = 0; x < cpu.length; x++) { const d = Math.abs(cpu[x] - gpu[x]); maxDiff = Math.max(maxDiff, d); sumSq += d * d; sumRef += cpu[x] * cpu[x]; }
  return { maxDiff, rmsRel: Math.sqrt(sumSq / Math.max(sumRef, 1e-300)) };
}

test('both engines carry the spectral gases alike: per-cell sums of absorbed and outgoing radiation, clear and all-sky, at N=6', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const { cpu, gpu, C } = await engines({ ...SPECTRAL, clearSkyPass: true });
  const slots = { absorbedSolar: 'ABSSUM', atmosphereSolar: 'ATMSUM', outgoingLongwave: 'OLRSUM', clearAbsorbedSolar: 'ABSCLRSUM', clearOutgoingLongwave: 'OLRCLRSUM' };
  const compare = async (names, limit) => {
    const device = await gpu.gpu.downloadPhysics();
    return names.map((name) => {
      const s = stats(cpu.radiation.summed[name], device[slots[name]].subarray(0, C));
      assert.ok(s.rmsRel < limit, `${name}: per-cell rms ${s.rmsRel}`);
      return `${name} rms ${s.rmsRel.toExponential(1)} (max ${s.maxDiff.toExponential(1)})`;
    });
  };
  for (let s = 0; s < 4; s++) { cpu.step(DT); await gpu.step(DT); }
  const early = await compare(Object.keys(slots), 1e-5);
  for (let s = 4; s < 24; s++) { cpu.step(DT); await gpu.step(DT); }
  const later = await compare(['clearAbsorbedSolar', 'clearOutgoingLongwave'], 1e-4);
  console.log(`spectral gases, CPU against GPU per-cell sums at N=6, 4 steps: ${early.join(', ')}; 24 steps: ${later.join(', ')}`);
  gpu.destroy();
});
