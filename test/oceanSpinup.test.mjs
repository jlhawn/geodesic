import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, writeFileSync, mkdtempSync, rmSync, readdirSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { execFileSync } from 'node:child_process';
import { Grid } from '../js/grid.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { decodeState, encodeState } from '../js/stateFile.module.js';
import { decodeForcing, forcingName } from '../js/forcing.module.js';
import { savedSubsidence } from '../js/physics/regrid.module.js';
import { FREEZING_POINT } from '../js/physics/ice.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};
const { createForcingRecorder } = gpuAvailable ? await import('../js/gpu/forcing.gpu.js') : {};
const { readRanges } = gpuAvailable ? await import('../js/gpu/device.module.js') : {};

const root = new URL('..', import.meta.url).pathname;
const topography = topographyFromInt16(readFileSync(join(root, 'data/topography_0p25.bin')).buffer);
const N = 6;

async function freshModel() {
  const model = await createGpuModel(new Grid(N), { topography });
  const { state, mesh } = model;
  initializeState(model, { geostrophic: !model.surfaceGeopotential }).forEach((values, a) => state[a].set(values));
  for (let i = 0; i < mesh.nCells; i++) if (model.geography.land[i]) state[6][i] = 0;
  model.load();
  model.ocean.initialize(state[3], state[6]);
  model.land.initialize();
  return model;
}

async function loadInto(model, saved) {
  const { state } = model;
  ['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'].forEach((name, a) => state[a].set(saved[name]));
  model.seaIce.load(state[6], saved.concentration ?? null);
  model.radiation.mlmSubsidence.set(savedSubsidence(saved, model));
  model.time = saved.time;
  model.load();
  model.ocean.load(saved.ocean, state[3], state[6]);
  model.land.load({ soil: Float64Array.from(saved.land.soil), snow: Float64Array.from(saved.land.snow), vegetation: Float64Array.from(saved.land.vegetation) });
}

test('the recorded stress and fluxes are the means of what each step handed the ocean', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const model = await freshModel();
  const { mesh, gpu, oceanEngine: ocean } = model, C = mesh.nCells, E = mesh.nEdges;
  await model.diagnostics();
  const recorder = await createForcingRecorder(model);
  const dt = 1350 * 16 / N, perDay = Math.round(86400 / dt);
  const stress = new Float64Array(E), net = new Float64Array(C), evaporation = new Float64Array(C);
  let oceanSteps = 0;
  for (let n = 1; n <= perDay; n++) {
    await model.step(dt);
    recorder.step();
    const [[s], [f, e]] = await Promise.all([readRanges(gpu.device, ocean.buffers.OD, [{ offset: ocean.layout.OD.STRESS, length: E }]), readRanges(gpu.device, gpu.buffers.PH, [{ offset: gpu.layout.PH.SFLUX, length: C }, { offset: gpu.layout.PH.EVAP, length: C }])]);
    if (n % ocean.everySteps === 0) { oceanSteps++; for (let x = 0; x < E; x++) stress[x] += s[x]; }
    for (let i = 0; i < C; i++) { net[i] += f[i]; evaporation[i] += e[i]; }
  }
  await model.diagnostics();
  const day = await decodeForcing(await recorder.day(1));
  assert.equal(day.steps, perDay);
  assert.equal(day.oceanSteps, oceanSteps);
  let largest = 0;
  for (let x = 0; x < E; x++) largest = Math.max(largest, Math.abs(stress[x] / oceanSteps));
  assert.ok(largest > 1e-3, `the stress is ${largest} N/m² at most`);
  for (let x = 0; x < E; x++) assert.ok(Math.abs(day.fields.stress[x] - stress[x] / oceanSteps) <= 1e-6 * largest, `stress on edge ${x}: ${day.fields.stress[x]} recorded, ${stress[x] / oceanSteps} used`);
  for (let i = 0; i < C; i++) {
    assert.ok(Math.abs(day.fields.netFlux[i] - net[i] / perDay) <= 1e-3, `net flux at ${i}`);
    assert.ok(Math.abs(day.fields.evaporation[i] - evaporation[i] / perDay) <= 1e-9, `evaporation at ${i}`);
  }
  model.destroy();
});

test('two recorded days replayed over the ocean alone keep its SST with the coupled run and load back into the coupled model', { skip: !gpuAvailable && 'webgpu not installed' }, async (t) => {
  const dir = mkdtempSync(join(tmpdir(), 'oceanSpinup-'));
  t.after(() => rmSync(dir, { recursive: true, force: true }));
  const fresh = await freshModel();
  const ocean = await fresh.ocean.serialize(), land = await fresh.land.serialize();
  const [pi, theta, u, surfaceT, q, qc, ice] = fresh.state;
  writeFileSync(join(dir, 'rec_day0000.bin'), encodeState({ N, K: fresh.core.K, day: 0, time: 0, terrain: !!fresh.surfaceGeopotential, pi, theta, u, surfaceT, q, qc, ice, concentration: fresh.seaIce.concentration, mlmSubsidence: fresh.radiation.mlmSubsidence, ocean: { h: ocean.h, u: ocean.u, T: ocean.T, S: ocean.S, eta: ocean.eta }, land }));
  fresh.destroy();
  const forcing = join(dir, 'forcing');
  const run = (script, env) => execFileSync(process.execPath, [join(root, 'scripts', script)], { cwd: root, env: { ...process.env, N: String(N), OUT: dir, ...env }, encoding: 'utf8' });
  run('spinup.mjs', { TAG: 'rec', RECORD: forcing, DAYS: '2', MINUTES: '100' });
  assert.deepEqual(readdirSync(forcing).sort(), [forcingName(1), forcingName(2)]);
  run('oceanSpinup.mjs', { TAG: 'alone', STATE: join(dir, 'rec_day0000.bin'), FORCING: forcing, YEARS: '1', DAYS_PER_YEAR: '2' });

  const [alone, coupled] = await Promise.all(['alone_year0001.bin', 'rec_day0002.bin'].map(async (name) => decodeState(new Uint8Array(readFileSync(join(dir, name))))));
  const recorded = await decodeForcing(new Uint8Array(readFileSync(join(forcing, forcingName(2)))));
  const model = await createGpuModel(new Grid(N), { topography });
  const { mesh } = model, sea = model.geography.land.map((l) => !l), C = mesh.nCells;
  const sst = (s, i) => (s.ice[i] > 0 ? FREEZING_POINT : s.surfaceT[i]);
  let area = 0, square = 0, worst = 0;
  for (let i = 0; i < C; i++) {
    if (!sea[i]) continue;
    const off = sst(alone, i) - recorded.fields.sst[i];
    area += mesh.areaCell[i]; square += mesh.areaCell[i] * off * off;
    worst = Math.max(worst, Math.abs(sst(alone, i) - sst(coupled, i)));
  }
  console.log(`after two days alone: SST rms ${Math.sqrt(square / area).toFixed(3)} K from the recorded day, at most ${worst.toFixed(3)} K from the coupled run`);
  assert.ok(Math.sqrt(square / area) < 0.5, `SST rms ${Math.sqrt(square / area)} K from the recorded day`);
  assert.ok(worst < 0.5, `SST up to ${worst} K from the coupled run`);

  await loadInto(model, alone);
  await model.step(1350 * 16 / N);
  await model.sync();
  const d = await model.diagnostics();
  for (const array of model.state) for (const x of array) assert.ok(Number.isFinite(x));
  for (const key of ['meanSurfaceT', 'maxWind', 'oceanSpeed', 'oceanInteriorT']) assert.ok(Number.isFinite(d[key]), key);
  model.destroy();
});
