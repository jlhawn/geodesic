import { emptyBuffer, storageBuffer, readRanges } from './device.module.js';
import { physicsConstants, PHYSICS_FUNCTIONS, SEA_SURFACE_WGSL, snowOnSea } from './physics.gpu.js';
import { nearestLayer, STABILITY_SIGMA } from '../physics/radiation.module.js';
import { encodeForcing } from '../forcing.module.js';

/*
 * The GPU side of forcing.module.js for a model of model.gpu.js.
 *
 * createForcingRecorder sums, after each model step, what the step
 * handed the ocean and sea ice, and `day` turns the sums into one day's
 * forcing file. It counts the model's steps from its own creation to
 * know which were ocean steps, so it must be created before the model's
 * first step and see every step after; `step` records through the
 * model's core, so inside model.stepBatch it belongs in the per-step
 * callback, where it lands after its step. Rain and runoff come from the
 * running totals rather than the sums: rain from PH CONV + COND, runoff
 * from model.land.runoff, which the model's diagnostics frame advances,
 * so `day` belongs after the day's model.diagnostics().
 *
 * createForcedOcean steps the ocean and sea ice alone under a recorded
 * day. Each `step` is one atmosphere step of the physics kernel's
 * sea-cell surface (SEA_SURFACE_WGSL) under the recorded net flux plus
 * the restoring −restore·(SST − recorded SST) on the open water, and of
 * the recorded snowfall; every everySteps-th, counted from creation as
 * the coupled model counts, also takes the recorded freshwater and steps
 * the ocean as advanceCoupled does, with the recorded stress in place of
 * the lowest wind's, and returns true.
 * The lead/ice split of the flux takes the recorded shortwave reaching
 * the surface times the diffuse albedo contrast of the ice as it stands.
 * Surface temperature, ice and concentration stay on the device in the
 * atmosphere's buffers, where the coupled model keeps them.
 */
const WORKGROUP = 64;

function seq(names) { const out = {}; let off = 0; for (const [name, n] of names) { out[name] = off; off += n; } out.total = off; return out; }

function compile(model, name, body, prefix, offsets, out) {
  const { gpu, oceanEngine: ocean } = model;
  const { device, buffers } = gpu;
  const phys = gpu.physics;
  const constants = physicsConstants({ ...phys, kTop: gpu.kTop, stratusLayer: nearestLayer(gpu.sigmaMid, phys.stratusSigma), stabilityLayer: nearestLayer(gpu.sigmaMid, STABILITY_SIGMA) });
  const names = ['MI', 'MF', 'LV', 'IN', 'OUT', 'D', 'P', 'PH', 'OD'];
  const code = `${gpu.preludeConstants}
${names.map((n, b) => `@group(0) @binding(${b}) var<storage, read_write> ${n}: array<${n === 'MI' ? 'i32' : 'f32'}>;`).join('\n')}
${constants}
${PHYSICS_FUNCTIONS}
const O_STRESS: i32 = ${ocean.layout.OD.STRESS};
${Object.entries(offsets).filter(([k]) => k !== 'total').map(([k, v]) => `const ${prefix}${k}: i32 = ${v};`).join('\n')}
@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let n = i32(id.x) + i32(id.y) * ${65535 * WORKGROUP};
${body}
}`;
  const layout = device.createBindGroupLayout({ entries: names.map((_, binding) => ({ binding, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } })) });
  const pipeline = device.createComputePipeline({ label: name, layout: device.createPipelineLayout({ bindGroupLayouts: [layout] }), compute: { module: device.createShaderModule({ code, label: name }), entryPoint: 'main' } });
  const params = storageBuffer(device, new Float32Array(8));
  const group = device.createBindGroup({ layout, entries: [buffers.MI, buffers.MF, buffers.LV, buffers.S, out, buffers.D, params, buffers.PH, ocean.buffers.OD].map((buffer, binding) => ({ binding, resource: { buffer } })) });
  return function run(values, count) {
    gpu.writeParams(Float32Array.from({ length: 8 }, (_, k) => values[k] ?? 0), params);
    gpu.compute((pass) => {
      pass.setPipeline(pipeline);
      pass.setBindGroup(0, group);
      const groups = Math.ceil(count / WORKGROUP);
      pass.dispatchWorkgroups(Math.min(groups, 65535), Math.ceil(groups / 65535));
    });
  };
}

export async function createForcingRecorder(model) {
  const { gpu, mesh, oceanEngine: ocean } = model;
  const { device } = gpu;
  const C = mesh.nCells, E = mesh.nEdges, PH = gpu.layout.PH;
  const SUMMED = ['netFlux', 'shortwave', 'shortwaveDown', 'sensible', 'evaporation', 'surfaceT', 'sst', 'ice', 'concentration', 'snowfall'];
  const A = seq([...SUMMED.map((name) => [name, C]), ['seen', C], ['stress', E]]);
  const acc = emptyBuffer(device, 4 * A.total);
  const run = compile(model, 'recordForcing', `
  if (n < E && P[0] > 0.5) { OUT[A_stress + n] += OD[O_STRESS + n]; }
  if (n >= C) { return; }
  let i = n; let ice = IN[S_ICE + i]; let conc = PH[PH_CONC + i];
  OUT[A_netFlux + i] += PH[PH_SFLUX + i];
  OUT[A_shortwave + i] += PH[PH_ABS + i] - PH[PH_ATMSW + i];
  OUT[A_shortwaveDown + i] += PH[PH_SWDN + i];
  OUT[A_sensible + i] += PH[PH_SH + i];
  OUT[A_evaporation + i] += PH[PH_EVAP + i];
  OUT[A_surfaceT + i] += IN[S_TS + i];
  OUT[A_sst + i] += select(IN[S_TS + i], FREEZING, ice > 0.0);
  OUT[A_ice + i] += ice;
  OUT[A_concentration + i] += select(0.0, select(conc, 1.0, conc <= 0.0), ice > 0.0);
  let fallen = PH[PH_CONV + i] + PH[PH_COND + i];
  let amount = select(fallen, fallen - OUT[A_seen + i], fallen >= OUT[A_seen + i]);
  OUT[A_seen + i] = fallen;
  let bottom = (K - 1) * C + i;
  if (PH[PH_LAND + i] < 0.5 && amount > 0.0 && IN[S_TH + bottom] * D[D_EXM + bottom] < MELTING) { OUT[A_snowfall + i] += amount; }`, 'A_', A, acc);

  const readRain = async () => { const [conv, cond] = await readRanges(device, gpu.buffers.PH, [{ offset: PH.CONV, length: C }, { offset: PH.COND, length: C }]); return Float64Array.from(conv, (x, i) => x + cond[i]); };
  let rainBefore = await readRain(), runoffBefore = model.land ? Float64Array.from(model.land.runoff) : null, timeBefore = model.time;
  device.queue.writeBuffer(acc, 4 * A.seen, Float32Array.from(rainBefore));
  let counted = 0, steps = 0, oceanSteps = 0;

  return {
    step() {
      const oceanStep = ++counted % ocean.everySteps === 0;
      steps++;
      if (oceanStep) oceanSteps++;
      run([oceanStep ? 1 : 0], Math.max(C, E));
    },
    async day(day) {
      const [views, rain] = await Promise.all([readRanges(device, acc, [...SUMMED.map((name) => ({ offset: A[name], length: C })), { offset: A.stress, length: E }]), readRain()]);
      const seconds = model.time - timeBefore;
      const fields = Object.fromEntries(SUMMED.map((name, n) => [name, Float32Array.from(views[n], (x) => x / (name === 'snowfall' ? seconds : steps))]));
      fields.stress = Float32Array.from(views[SUMMED.length], (x) => x / Math.max(1, oceanSteps));
      fields.rain = Float32Array.from(rain, (x, i) => (x - rainBefore[i]) / seconds);
      fields.runoff = model.land ? Float32Array.from(model.land.runoff, (x, i) => (x - runoffBefore[i]) / seconds) : new Float32Array(C);
      const bytes = encodeForcing({ N: Math.round(Math.sqrt((C - 2) / 10)), day, time: model.time, seconds, steps, oceanSteps }, fields);
      device.queue.writeBuffer(acc, 0, new Float32Array(A.seen));
      device.queue.writeBuffer(acc, 4 * A.stress, new Float32Array(E));
      rainBefore = rain; runoffBefore = model.land ? Float64Array.from(model.land.runoff) : null; timeBefore = model.time;
      steps = 0; oceanSteps = 0;
      return bytes;
    },
  };
}

export function createForcedOcean(model) {
  const { gpu, mesh, oceanEngine: ocean } = model;
  const { device } = gpu;
  const C = mesh.nCells, E = mesh.nEdges, PH = gpu.layout.PH;
  const F = seq([['netFlux', C], ['shortwaveDown', C], ['sst', C], ['snowfall', C], ['stress', E]]);
  const forcing = emptyBuffer(device, 4 * F.total);
  const run = compile(model, 'forcedSurface', `
  let i = n; if (i >= C) { return; }
  if (PH[PH_LAND + i] > 0.5) { return; }
  let dt = P[0];
  let skin = IN[S_TS + i]; let ice = IN[S_ICE + i];
  let conc0 = PH[PH_CONC + i]; let snow0 = PH[PH_SNOW + i];
  let cover = select(0.0, select(conc0, 1.0, conc0 <= 0.0), ice > 0.0);
  let pull = -P[1] * (select(skin, FREEZING, ice > 0.0) - OUT[FC_sst + i]);
  let net = OUT[FC_netFlux + i] + (1.0 - cover) * pull;
  let contrast = OUT[FC_shortwaveDown + i] * (surfaceAlbedo(ice, ALB_DIF_WATER, snow0) - ALB_DIF_WATER) + pull;
  let ocean = PH[PH_OFLUX + i]; let capacity = PH[PH_CAP + i];
  var T = skin; var h = ice;
  ${SEA_SURFACE_WGSL}
  IN[S_TS + i] = T; IN[S_ICE + i] = h;
  let snowfall = OUT[FC_snowfall + i] * dt;
  if (snowfall > 0.0) {
    ${snowOnSea('snowfall')}
  }`, 'FC_', F, forcing);
  let counted = 0, perInterval = 0;

  return {
    /*
     * The day's forcing (forcing.module.js fields), with the freshwater
     * written as the amounts of one ocean step of `oceanDt` seconds.
     */
    setDay(fields, oceanDt) {
      for (const name of ['netFlux', 'shortwaveDown', 'sst', 'snowfall', 'stress']) device.queue.writeBuffer(forcing, 4 * F[name], Float32Array.from(fields[name]));
      device.queue.writeBuffer(gpu.buffers.PH, 4 * PH.EVAP, Float32Array.from(fields.evaporation));
      device.queue.writeBuffer(gpu.buffers.PH, 4 * PH.RAIN, Float32Array.from(fields.rain, (x) => x * oceanDt));
      device.queue.writeBuffer(gpu.buffers.PH, 4 * PH.RUNOFF, Float32Array.from(fields.runoff, (x) => x * oceanDt));
      perInterval = oceanDt;
    },
    step(dt, restore) {
      run([dt, restore], C);
      if (++counted % ocean.everySteps !== 0) return false;
      const oceanDt = ocean.everySteps * dt;
      if (Math.abs(oceanDt - perInterval) > 1e-6 * oceanDt) throw new Error(`setDay was given an ocean step of ${perInterval} s, not ${oceanDt} s`);
      ocean.forgetAccumulated('rain');
      ocean.forgetAccumulated('runoff');
      ocean.accumulateFreshwater(oceanDt);
      ocean.readSurfaceFromAtmosphere();
      const encoder = device.createCommandEncoder();
      encoder.copyBufferToBuffer(forcing, 4 * F.stress, ocean.buffers.OD, 4 * ocean.layout.OD.STRESS, 4 * E);
      device.queue.submit([encoder.finish()]);
      ocean.step(oceanDt);
      ocean.mixedLayer(oceanDt);
      ocean.salt(oceanDt);
      ocean.writeSurface(oceanDt);
      return true;
    },
  };
}
