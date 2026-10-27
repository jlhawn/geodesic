import { test } from 'node:test';
import assert from 'node:assert/strict';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { createSigmaCore, sigmaInterfaces, VIRTUAL_FACTOR } from '../js/dynamics/sigmaCore.module.js';
import { createRK4Arrays } from '../js/dynamics/integrators.module.js';
import { divergence } from '../js/dynamics/operators.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { withCadence } from '../js/cadence.module.js';
import { FREEZING_POINT } from '../js/physics/ice.module.js';
import { geographyOptions, readTopography, FRESH_GEOGRAPHY } from '../js/geography.module.js';
import { CLIMATOLOGY_FILE, loadClimatology } from '../js/ocean/climatology.module.js';
import { FLAKE } from '../js/physics/lake.module.js';

let gpuAvailable = true;
try { await import('webgpu'); } catch { gpuAvailable = false; }
const { createGpuModel } = gpuAvailable ? await import('../js/gpu/model.gpu.js') : {};
const { createGpuCore } = gpuAvailable ? await import('../js/gpu/core.gpu.js') : {};

const N = 16, DT = 1350 * 16 / N, SPINUP_STEPS = 16, MEASURED_STEPS = 6;
const cadence = withCadence({ ocean: false }, DT, N);
const OPTIONS = { ocean: false, levels: sigmaInterfaces('bl36'), ...cadence };

/*
 * An aquaplanet over a slab six hours from initializeState, so that the
 * boundary layer, the clouds and the convection are all under way.
 */
let spun = null;
function spunUp() {
  if (spun) return spun;
  const model = createModel(new Grid(N), OPTIONS);
  initializeState(model).forEach((values, s) => model.state[s].set(values));
  for (let n = 0; n < SPINUP_STEPS; n++) model.step(DT);
  spun = model;
  return model;
}

let spunBergeron = null;
function spunUpBergeron() {
  if (spunBergeron) return spunBergeron;
  const model = createModel(new Grid(N), { ...OPTIONS, moist: { bergeronAutoconversion: true } });
  initializeState(model).forEach((values, s) => model.state[s].set(values));
  for (let n = 0; n < SPINUP_STEPS; n++) model.step(DT);
  spunBergeron = model;
  return model;
}

/*
 * The spun-up air over a layered ocean started under it.
 */
function overOcean(surface = {}) {
  const slab = spunUp(), model = createModel(new Grid(N), { ...OPTIONS, ocean: {}, surface });
  slab.state.forEach((values, s) => model.state[s].set(values));
  model.seaIce.load(model.state[6], slab.seaIce.concentration);
  model.ocean.initialize(model.state[3], model.state[6]);
  model.time = slab.time;
  return model;
}

/*
 * The same with the sea's skin, its warm layer started at 1 K over the open
 * water so that the air sees T_s ≠ T_slab from the first step.
 */
function overOceanSkin() {
  const model = overOcean({ warmLayer: true, coolSkin: true });
  for (let i = 0; i < model.mesh.nCells; i++) if (!(model.state[6][i] > 0)) model.skin.warm[i] = 1;
  return model;
}

/*
 * The air's enthalpy c_p T, latent heat L q and kinetic energy, and the
 * surface's store (seaIce.energy, with an ocean's heat capacity over open
 * water), as global means in J/m², summed in double precision with the
 * Exner function of each state's own π.
 */
function bookkeeper(model) {
  const { mesh, core, seaIce, ocean } = model, { K, C, E, dSigma, cp, g } = core.diagnostics, L = model.moist.latentHeat;
  const { areaCell, cellsOnEdge, dcEdge, dvEdge } = mesh, area = areaCell.reduce((s, a) => s + a, 0);
  const look = createSigmaCore(mesh, { levels: core.levels }), exner = look.diagnostics.exnerLayer;
  return (state) => {
    const [pi, theta, u, surfaceT, q, , ice] = state;
    look.diagnose(pi, theta);
    let H = 0, Lq = 0, KE = 0, S = 0;
    for (let i = 0; i < C; i++) {
      for (let k = 0; k < K; k++) { const n = k * C + i, m = areaCell[i] * pi[i] * dSigma[k] / g; H += m * cp * theta[n] * exner[n]; Lq += m * L * q[n]; }
      S += areaCell[i] * (ocean && !(ice[i] > 0) ? ocean.capacity[i] * (surfaceT[i] - FREEZING_POINT) : seaIce.energy(surfaceT[i], ice[i], seaIce.snow[i], seaIce.cover(i, ice[i])));
    }
    for (let k = 0; k < K; k++) for (let e = 0; e < E; e++) {
      const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
      KE += 0.5 * dcEdge[e] * dvEdge[e] * u[k * E + e] ** 2 * 0.5 * (pi[a] + pi[b]) * dSigma[k] / g;
    }
    return { H: H / area, Lq: Lq / area, KE: KE / area, S: S / area };
  };
}

const largest = (rows) => Math.max(...rows.map(Math.abs));
const mean = (rows) => rows.reduce((s, x) => s + x, 0) / rows.length;
const show = (x) => x.toExponential(1);

for (const [label, build] of [['over a slab', spunUp], ['over a slab with bergeronAutoconversion', spunUpBergeron], ['over a layered ocean', () => overOcean()], ['over a layered ocean with the warm layer and the cool skin', overOceanSkin]]) test(`in double precision the CPU model's step ${label} keeps the energy of the air and its surface through the ∇⁴ closure, the boundary layer's mixing, the moist physics and the momentum mixing with the heat it returns, and the RK4 dynamics within its truncation`, () => {
  const model = build(), { core, boundaryLayer, phases, state, radiation } = model, { K, C, E } = core.diagnostics;
  const books = bookkeeper(model), rk4 = createRK4Arrays(state.map((a) => a.length));
  const shadow = state.map((a) => Float64Array.from(a));
  const record = { dynamics: [], closure: [], boundaryLayer: [], adjust: [], momentum: [] };
  for (let n = 0; n < MEASURED_STEPS; n++) {
    radiation.setTime(model.time);
    let a = books(state);
    rk4(model.tendency, state, DT);
    let b = books(state);
    record.dynamics.push((b.H + b.Lq + b.KE - a.H - a.Lq - a.KE) / DT);
    const slabBefore = Float64Array.from(state[3]), iceBefore = Float64Array.from(state[6]), penBefore = model.ocean && model.ocean.penetrated ? Float64Array.from(model.ocean.penetrated) : null;
    phases.physics(0, C, DT, model.totals);
    if (model.skin) {
      let seen = 0;
      for (let i = 0; i < C; i++) if (!(iceBefore[i] > 0) && !(state[6][i] > 0)) {
        const paid = (DT * (radiation.surfaceFlux[i] + model.seaIce.oceanFlux[i]) - (penBefore ? model.ocean.penetrated[i] - penBefore[i] : 0)) / model.ocean.capacity[i];
        assert.ok(Math.abs(state[3][i] - slabBefore[i] - paid) <= 1e-9, `cell ${i}: the slab gains ${state[3][i] - slabBefore[i]} K for ${paid} K of flux at the skin`);
        if (model.skin.offset(i, model.ocean.capacity[i]) !== 0) seen++;
      }
      assert.ok(seen > 0, 'the air sees the skin');
    }
    a = books(state);
    phases.closure(0, K, DT);
    b = books(state);
    record.closure.push((b.H + b.Lq - a.H - a.Lq) / DT);
    let momentum = b.KE - a.KE;
    shadow.forEach((copy, s) => copy.set(state[s]));
    const unmixed = books(shadow);
    for (let i = 0; i < C; i++) boundaryLayer.mixColumn(i, shadow[0], shadow[1], shadow[4], shadow[5], DT);
    const mixed = books(shadow);
    record.boundaryLayer.push((mixed.H + mixed.Lq - unmixed.H - unmixed.Lq) / DT);
    a = b;
    phases.adjust(0, C, DT);
    b = books(state);
    record.adjust.push((b.H + b.Lq + b.S - a.H - a.Lq - a.S) / DT);
    a = b;
    phases.mixMomentum(0, E, DT);
    phases.dissipate(0, C);
    b = books(state);
    momentum += b.H + b.Lq + b.KE - a.H - a.Lq - a.KE;
    record.momentum.push(momentum / DT);
    model.time += DT;
    model.radiationSteps++;
  }
  for (const [name, rows] of Object.entries(record)) console.log(`CPU ${label} ${name}: ${show(mean(rows))} W/m² on average, ${show(largest(rows))} at most`);
  for (const name of ['closure', 'boundaryLayer', 'adjust', 'momentum']) assert.ok(largest(record[name]) < 1e-6, `${name}: ${largest(record[name])} W/m²`);
  assert.ok(Math.abs(mean(record.dynamics)) < 0.02, `dynamics: ${mean(record.dynamics)} W/m²`);
});

test("in double precision the CPU's adjustment keeps the energy of the air and its water while plumes root above the 50 hPa stack under plumeDeparture 'elevated', over a spun-up slab with its lowest 60 hPa cooled by 6 K in two columns of five", () => {
  const model = createModel(new Grid(N), { ...OPTIONS, moist: { plumeDeparture: 'elevated', stratocumulusInversion: null } });
  initializeState(model).forEach((values, s) => model.state[s].set(values));
  for (let n = 0; n < SPINUP_STEPS; n++) model.step(DT);
  const { K, C, sigmaMid, exnerLayer } = model.core.diagnostics, [pi, theta] = model.state, { phases, moist } = model;
  for (let i = 0; i < C; i++) {
    if (i % 5 > 1) continue;
    model.core.diagnoseColumn(i, pi, theta, model.state[4], model.state[5]);
    for (let k = 0; k < K; k++) if (sigmaMid[k] > 1 - 60e2 / pi[i]) theta[k * C + i] -= 6 / exnerLayer[k * C + i];
  }
  const books = bookkeeper(model), rows = [];
  let elevated = 0;
  for (let n = 0; n < MEASURED_STEPS; n++) {
    model.step(DT);
    const a = books(model.state);
    phases.adjust(0, C, DT);
    const b = books(model.state);
    rows.push((b.H + b.Lq + b.S - a.H - a.Lq - a.S) / DT);
    for (let i = 0; i < C; i++) {
      let stack = K - 1;
      while (stack > 0 && pi[i] * sigmaMid[stack - 1] >= pi[i] - 50e2) stack--;
      if (moist.momentum.source[i] < stack) elevated++;
    }
  }
  console.log(`CPU over a slab with plumeDeparture elevated, ${elevated} plume-steps rooted above the stack: adjust ${show(mean(rows))} W/m² on average, ${show(largest(rows))} at most`);
  assert.ok(elevated > 0, 'plumes rooted above the stack');
  assert.ok(largest(rows) < 1e-6, `adjust: ${largest(rows)} W/m²`);
});

/*
 * The energy tendency (W/m²) of one evaluation of the dynamics, from the
 * tendencies of π, θ and u in `out`: Σ c_p Π (π dθ/dt + (1 + κ) θ dπ/dt)
 * Δσ/g over the cells and layers (Ππ goes as π^(1+κ)) and the kinetic
 * energy's change on the edges. Water's L q is carried in flux form and
 * keeps its total, so it is left out.
 */
function energyTendency(model, state, out) {
  const { mesh, core } = model, { K, C, E, dSigma, cp, g, kappa } = core.diagnostics;
  const { areaCell, cellsOnEdge, dcEdge, dvEdge } = mesh;
  const look = createSigmaCore(mesh, { levels: core.levels }), exner = look.diagnostics.exnerLayer;
  look.diagnose(state[0], state[1]);
  let rate = 0, area = 0;
  for (let i = 0; i < C; i++) {
    area += areaCell[i];
    for (let k = 0; k < K; k++) {
      const n = k * C + i;
      rate += areaCell[i] * dSigma[k] / g * cp * exner[n] * (out[1][n] * state[0][i] + (1 + kappa) * state[1][n] * out[0][i]);
    }
  }
  for (let k = 0; k < K; k++) for (let e = 0; e < E; e++) {
    const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1], x = k * E + e, u = state[2][x];
    rate += 0.5 * dcEdge[e] * dvEdge[e] * dSigma[k] / g * (u * out[2][x] * 0.5 * (state[0][a] + state[0][b]) + 0.25 * u * u * (out[0][a] + out[0][b]));
  }
  return rate / area;
}

/*
 * −∫ ω (α_v − α) dm (W/m²), the work moist buoyancy does, with
 * ω = σ(∂π/∂t + v·∇π) + πσ̇ from the CPU core's own fluxes.
 */
function virtualWork(model, state) {
  const { mesh, core } = model, { K, C, E, dSigma, sigmaMid, g, R } = core.diagnostics;
  const out = state.map((a) => new Float64Array(a.length)), spread = new Float64Array(C);
  core.tendency(state, out);
  const { piSigmaDot, exnerLayer } = core.diagnostics, divFlux = new Float64Array(core.shared.divFlux);
  let work = 0, area = 0;
  for (let k = 0; k < K; k++) {
    divergence(mesh, state[2].subarray(k * E, (k + 1) * E), spread);
    for (let i = 0; i < C; i++) {
      const n = k * C + i, pi = state[0][i];
      const omega = sigmaMid[k] * (out[0][i] + divFlux[n] - pi * spread[i]) + 0.5 * (piSigmaDot[n] + piSigmaDot[n + C]);
      work -= mesh.areaCell[i] * pi * dSigma[k] / g * omega * R * state[1][n] * exnerLayer[n] * (VIRTUAL_FACTOR * state[4][n] - state[5][n]) / (pi * sigmaMid[k]);
    }
  }
  for (let i = 0; i < C; i++) area += mesh.areaCell[i];
  return work / area;
}

/*
 * The spun-up state rounded to single precision with its water ×WETTER,
 * which raises the virtual work tenfold and keeps it well above single
 * precision's floor, and its dry twin: θ_v as θ and no water, which
 * the dynamics moves as it moves the moist state, to rounding.
 */
const WETTER = 10;
function twins(model) {
  const { K, C } = model.core.diagnostics;
  const moist = model.state.map((a) => Float64Array.from(a, Math.fround));
  for (let n = 0; n < K * C; n++) { moist[4][n] = Math.fround(WETTER * moist[4][n]); moist[5][n] = Math.fround(WETTER * moist[5][n]); }
  const dry = moist.map((a) => Float64Array.from(a));
  for (let n = 0; n < K * C; n++) { dry[1][n] = Math.fround(moist[1][n] * (1 + VIRTUAL_FACTOR * moist[4][n] - moist[5][n])); dry[4][n] = 0; dry[5][n] = 0; }
  return { moist, dry };
}

test('the CPU dynamics changes the energy of moist air as it changes its dry twin\'s, so the work its θ_v buoyancy does is paid from c_p T', () => {
  const model = spunUp(), { moist, dry } = twins(model);
  const out = moist.map((a) => new Float64Array(a.length));
  model.core.tendency(moist, out);
  const wet = energyTendency(model, moist, out);
  model.core.tendency(dry, out);
  const plain = energyTendency(model, dry, out), work = virtualWork(model, moist);
  console.log(`CPU: energy tendency ${wet.toExponential(3)} W/m² moist, ${plain.toExponential(3)} its dry twin; the virtual work ${work.toExponential(3)}`);
  assert.ok(work > 0.02, `virtual work ${work} W/m²`);
  assert.ok(Math.abs(wet - plain) < 0.2 * work, `moist ${wet} against dry ${plain} W/m², virtual work ${work}`);
});

test('the GPU dynamics changes the energy of moist air as it changes its dry twin\'s, to single precision', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const model = spunUp(), { moist, dry } = twins(model), { K, C } = model.core.diagnostics;
  const referenceTheta = Float64Array.from({ length: K }, (_, k) => { let s = 0, a = 0; for (let i = 0; i < C; i++) { s += model.mesh.areaCell[i] * moist[1][k * C + i]; a += model.mesh.areaCell[i]; } return s / a; });
  const gpu = await createGpuCore(model.mesh, { levels: model.core.levels, dragCoefficient: 0, topDragDays: 0, spongeRates: null, spongeMeanRates: null, gravityWaves: false, referenceTheta });
  const rates = [];
  for (const state of [moist, dry]) {
    gpu.upload(state);
    await gpu.tendency();
    rates.push(energyTendency(model, state, await gpu.download(gpu.buffers.K1)));
  }
  for (const buffer of Object.values(gpu.buffers)) buffer.destroy();
  const [wet, plain] = rates, work = virtualWork(model, moist);
  console.log(`GPU: energy tendency ${wet.toExponential(3)} W/m² moist, ${plain.toExponential(3)} its dry twin; the virtual work ${work.toExponential(3)}`);
  assert.ok(Math.abs(wet - plain) < 0.2 * work, `moist ${wet} against dry ${plain} W/m², virtual work ${work}`);
});

/*
 * The GPU model holding the spun-up state after one step of its own, so
 * that its physics has laid down the boundary layer's inputs; then the
 * ∇⁴ closure's θ, q and qc passes, and the boundary layer's diagnosis
 * and mixing with the saturation adjustment that follows it, each
 * dispatched alone on the state as it stands. The bounds are two to
 * four times what an M1 Max leaves here.
 */
test('in single precision the GPU\'s ∇⁴ closure and boundary-layer mixing keep the air\'s c_p T + L q to their rounding', { skip: !gpuAvailable && 'webgpu not installed' }, async () => {
  const model = spunUp(), books = bookkeeper(model);
  const gpuModel = await createGpuModel(model.mesh, OPTIONS);
  model.state.forEach((a, s) => gpuModel.state[s].set(a));
  gpuModel.seaIce.concentration.set(model.seaIce.concentration);
  gpuModel.time = model.time;
  gpuModel.load();
  const { gpu } = gpuModel, { device, buffers, kernels } = gpu, { C } = model.core.diagnostics;
  const run = async (names, parameters) => {
    device.queue.writeBuffer(buffers.P, 0, Float32Array.from(parameters));
    const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
    for (const name of names) {
      pass.setPipeline(kernels[name]);
      pass.setBindGroup(0, device.createBindGroup({ layout: kernels[name].getBindGroupLayout(0), entries: [buffers.MI, buffers.MF, buffers.LV, buffers.S, buffers.K1, buffers.D, buffers.P, buffers.PH].map((buffer, binding) => ({ binding, resource: { buffer } })) }));
      pass.dispatchWorkgroups(Math.ceil(C / 64));
    }
    pass.end();
    device.queue.submit([encoder.finish()]);
  };
  const record = { closure: [], boundaryLayer: [] };
  for (let n = 0; n < 8; n++) {
    await gpuModel.step(DT);
    let before = books(await gpu.download());
    await run(['lapScalar1', 'lapScalar2'], [DT * model.core.nu4Theta, 0, 0, 0, 0, 0, 0, 0]);
    let after = books(await gpu.download());
    record.closure.push((after.H + after.Lq - before.H - before.Lq) / DT);
    before = after;
    await run(['pblDiagnose', 'adjustMix'], [DT, 0, 1, 0, 0, 0, 0, 0]);
    after = books(await gpu.download());
    record.boundaryLayer.push((after.H + after.Lq - before.H - before.Lq) / DT);
  }
  gpuModel.destroy();
  for (const [name, rows] of Object.entries(record)) {
    console.log(`GPU ${name}: ${show(mean(rows))} W/m² on average, ${show(largest(rows))} at most`);
    assert.ok(Math.abs(mean(rows)) < 2e-3 && largest(rows) < 5e-3, `${name}: ${mean(rows)} W/m² on average, ${largest(rows)} at most`);
  }
});

test('at N=32 FLake holds the n1024 lakes: dry to the ocean, finite for two days from a fresh start, each lake\'s heat changing by the surface flux it is given where nothing falls on it', async () => {
  const n = 32, dt = 1350 * 16 / n, topography = readTopography('etopo2022'), sheets = geographyOptions({ ...FRESH_GEOGRAPHY, mask: 'n1024' });
  const ocean = withCadence({ ocean: {} }, dt, n).ocean, radiation = withCadence({ radiation: { clearSkyPass: true } }, dt, n).radiation;
  const climatology = await loadClimatology(new URL(`../${CLIMATOLOGY_FILE}`, import.meta.url).pathname);
  const model = createModel(new Grid(n), { topography, geography: sheets.geography, ocean: { ...ocean, climatology }, radiation, land: { ...sheets.land }, lake: { model: 'flake' } });
  const { state, geography, lake } = model, C = model.mesh.nCells;
  assert.ok(lake.cells.length >= 2, `${lake.cells.length} lake cells`);
  for (const i of lake.cells) assert.deepEqual([geography.land[i], model.ocean.cellOcean[i]], [0, 0], `lake cell ${i}`);
  initializeState(model, { geostrophic: !model.surfaceGeopotential }).forEach((values, a) => state[a].set(values));
  for (let i = 0; i < C; i++) if (geography.land[i]) state[6][i] = 0;
  if (model.load) model.load();
  model.ocean.initialize(state[3], state[6]);
  model.land.initialize();
  let checked = 0, worst = 0;
  for (let s = 0; s < 2 * 86400 / dt; s++) {
    const heat = lake.cells.map((i) => (lake.mean[i] === 0 ? null : lake.enthalpy(i))), water = lake.cells.map((i) => lake.passed[i]);
    model.step(dt);
    lake.cells.forEach((i, k) => {
      if (heat[k] === null || lake.passed[i] - water[k] + model.radiation.evaporation[i] * dt > 1e-12) return;
      checked++;
      const given = model.radiation.surfaceFlux[i] * dt, grain = FLAKE.rhoWater * FLAKE.cWater * lake.depth[i] * lake.mean[i] * Number.EPSILON;
      worst = Math.max(worst, Math.abs(lake.enthalpy(i) - heat[k] - given) / (1e-9 * Math.max(1, Math.abs(given)) + 4 * grain));
    });
  }
  const nonFinite = (fields) => Object.entries(fields).filter(([, v]) => ArrayBuffer.isView(v) && v.some((x) => !Number.isFinite(x))).map(([name]) => name);
  assert.deepEqual([...nonFinite(Object.fromEntries(['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'].map((name, k) => [name, state[k]]))), ...nonFinite(lake.serialize())], []);
  for (const i of lake.cells) assert.ok(state[3][i] > 240 && state[3][i] < 315 && lake.mixedDepth[i] > 0 && lake.mixedDepth[i] <= lake.depth[i], `lake cell ${i}: ${state[3][i]} K, mixed ${lake.mixedDepth[i]} of ${lake.depth[i]} m`);
  assert.ok(checked > 0 && worst < 1, `${checked} lake steps, worst heat-budget miss ${worst} of its bound`);
});
