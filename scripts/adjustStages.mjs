// The engines' adjust step and boundary-layer diagnosis, stage by stage,
// from one saved state. The CPU model takes the state (with the land and
// terrain) through its physics phase; that state, rounded to single
// precision, and the phase's boundary-layer fields go to a GPU core, which
// runs the adjust kernel alone, cut short after the mixing, after the
// first condensation, after the plumes, or whole; the CPU runs the same
// stages column by column. Then both diagnose the boundary layer from the
// CPU's adjusted state, the GPU's friction velocity and surface buoyancy
// flux taken from the CPU. Per stage it prints the layers below and above
// the mixing top whose θ (1e-3 K), q or qc (1e-7) part, the plumes whose
// base flux parts by 1 %, and the regimes and mixing tops that part.
// The GPU's heights are above the ground, the CPU's above sea level: the
// mixing top and depth are moved by the terrain on the way.
//   node scripts/adjustStages.mjs runs/eleven64_day1825.bin
//   BC=cloudLayer DUMP=1234 node scripts/adjustStages.mjs runs/eleven64_day1825.bin
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { decodeState, savedLevels } from '../js/stateFile.module.js';
import { getDevice } from '../js/gpu/device.module.js';
import { createGpuCore } from '../js/gpu/core.gpu.js';
import { readTopography } from './figures/figureState.mjs';

const BC = process.env.BC ?? 'uniform', DUMP = process.env.DUMP ? Number(process.env.DUMP) : -1;
const saved = await decodeState(new Uint8Array(readFileSync(process.argv[2])));
const levels = savedLevels(saved);
const { device } = await getDevice();
const plain = device.createShaderModule.bind(device);
let stage = 3;
device.createShaderModule = (descriptor) => {
  if (descriptor.label !== 'adjust' && descriptor.label !== 'pblDiagnose') return plain(descriptor);
  let code = descriptor.code;
  const swap = (from, to) => { if (code.split(from).length !== 2) throw new Error(`${descriptor.label} kernel text moved: ${from}`); code = code.replace(from, to); };
  if (descriptor.label === 'adjust') {
    swap('  let pi = IN[S_PI + i]; let dt = P[0]; let bottom = K - 1;\n  var mixes = false;\n', '  let pi = IN[S_PI + i]; let dt = P[0]; let bottom = K - 1;\n  diagnoseColumn(i);\n  var mixes = false;\n');
    swap('  diagnoseColumn(i);\n  saturateColumn(i, pi);\n  let produced = plumeColumn(i, pi, dt);\n  if (PH[PH_CUMF + i] > 0.0) { saturateColumn(i, pi); }\n',
      `  if (${stage} == 0) { return; }\n  diagnoseColumn(i);\n  saturateColumn(i, pi);\n  if (${stage} == 1) { return; }\n  let produced = plumeColumn(i, pi, dt);\n  if (PH[PH_CUMF + i] > 0.0) { saturateColumn(i, pi); }\n  if (${stage} == 2) { return; }\n`);
  } else {
    swap('  let friction = sqrt(PH[PH_DRAG + i]) * xWind(i, speed);', '  let friction = PH[PH_USTAR + i];');
    const buoyancy = code.match(/ {2}let buoyancy = select\(GRAV \/ IN\[S_TH \+ base\][^\n]*\n/);
    if (!buoyancy) throw new Error('pblDiagnose kernel text moved: buoyancy');
    swap(buoyancy[0], '  let buoyancy = PH[PH_BUOY + i];\n');
  }
  return plain({ ...descriptor, code });
};

const moistOptions = { boundaryCondensation: BC };
const model = createModel(new Grid(saved.N), { ocean: false, levels, moist: moistOptions, topography: readTopography() });
const { core, moist, boundaryLayer: bl, radiation, mesh } = model, C = mesh.nCells, { K, geopotential, g } = core.diagnostics, dt = 1350 * 16 / saved.N;
['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'].forEach((name, a) => model.state[a].set(saved[name]));
if (saved.boundaryDepth) bl.depth.set(saved.boundaryDepth);
model.time = saved.time;
model.seaIce.load(model.state[6], saved.concentration ?? null);
if (model.land && saved.land) model.land.load(saved.land);
if (saved.windSpeed) model.surface.windSpeed.set(saved.windSpeed);
if (saved.exchangeWind && model.exchange) model.exchange.wind.set(saved.exchangeWind);
if (saved.mlmGate) radiation.mlmGate.set(saved.mlmGate);
if (saved.mlmHeight) radiation.mlmHeight.set(saved.mlmHeight);
if (saved.subcloudVirtual && saved.subcloudVirtual.length === moist.subcloudVirtual.length) moist.subcloudVirtual.set(saved.subcloudVirtual);
for (let i = 0; i < C; i++) core.diagnoseColumn(i, model.state[0], model.state[1], model.state[4], model.state[5]);
model.phases.physics(0, C, dt, model.totals);
const single = (a) => { for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]); };
model.state.forEach(single);
[bl.cloudLayer, bl.mixing, bl.mixingTop, bl.depth, bl.regime, bl.buoyancyFlux, bl.friction, radiation.stratiform, radiation.mlmGate, radiation.sensibleHeat, radiation.evaporation, moist.subcloudVirtual, model.seaIce.concentration].forEach(single);
const start = model.state.map((a) => Float64Array.from(a)), startSub = Float64Array.from(moist.subcloudVirtual), startGeo = Float64Array.from(geopotential);
const lift = (i) => (model.surfaceGeopotential ? model.surfaceGeopotential[i] / g : 0);
const ground = (a) => Float64Array.from(a, (v, i) => v - lift(i));
const mixed = (k, i) => startGeo[k * C + i] / g < bl.mixingTop[i];
const landCode = Float32Array.from(model.geography.land, (l, i) => (l ? (model.geography.iceSheet && model.geography.iceSheet[i] ? 2 : 1) : 0));

async function gpuCore() {
  const gpu = await createGpuCore(mesh, { levels, surfaceGeopotential: model.surfaceGeopotential, physics: { ...moistOptions, landed: true } });
  const put = (name, values) => device.queue.writeBuffer(gpu.buffers.PH, 4 * gpu.layout.PH[name], Float32Array.from(values));
  const run = (kernel) => {
    device.queue.writeBuffer(gpu.buffers.P, 0, Float32Array.from([dt, 0, 1, 0, 0, 0, 0, 0]));
    const { buffers } = gpu, group = device.createBindGroup({ layout: gpu.kernels[kernel].getBindGroupLayout(0), entries: [buffers.MI, buffers.MF, buffers.LV, buffers.S, buffers.K1, buffers.D, buffers.P, buffers.PH].map((buffer, binding) => ({ binding, resource: { buffer } })) });
    const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
    pass.setPipeline(gpu.kernels[kernel]); pass.setBindGroup(0, group); pass.dispatchWorkgroups(Math.ceil(C / 64)); pass.end();
    device.queue.submit([encoder.finish()]);
  };
  return { gpu, put, run };
}

const names = ['after the mixing', 'after the condensation', 'after the plumes', 'after the whole adjust step'];
for (stage = 0; stage < 4; stage++) {
  model.state.forEach((a, n) => a.set(start[n])); moist.subcloudVirtual.set(startSub);
  const [pi, theta, u, , q, qc] = model.state;
  if (stage === 3) model.phases.adjust(0, C, dt);
  else {
    for (let i = 0; i < C; i++) {
      bl.mixColumn(i, pi, theta, q, qc, dt);
      if (stage === 0) continue;
      core.diagnoseColumn(i, pi, theta, q, qc);
      moist.condenseColumn(i, pi, theta, q, qc);
      if (stage === 1) continue;
      moist.plumeColumn(i, pi, theta, q, qc, dt, u);
      if (moist.cumulusBaseFlux[i] > 0) moist.condenseColumn(i, pi, theta, q, qc);
    }
  }
  const { gpu, put, run } = await gpuCore();
  gpu.upload(start);
  gpu.uploadPhysics({ mlmGate: radiation.mlmGate, concentration: model.seaIce.concentration, land: landCode });
  put('MIX', bl.mixing); put('REGIME', bl.regime); put('MIXTOP', ground(bl.mixingTop)); put('STRAT', radiation.stratiform); put('DEPTH', ground(bl.depth));
  put('BUOY', bl.buoyancyFlux); put('USTAR', bl.friction); put('SUBTV', startSub); put('SH', radiation.sensibleHeat); put('EVAP', radiation.evaporation); put('CLOUDK', bl.cloudLayer);
  run('adjust');
  const after = await gpu.download(), ph = await gpu.downloadPhysics();
  const tally = () => ({ n: 0, theta: 0, q: 0, qc: 0, worstT: 0, worstQ: 0, worstQc: 0 });
  const below = tally(), above = tally(), thetaApart = new Set(), qcApart = new Set();
  let plumesApart = 0;
  for (let i = 0; i < C; i++) {
    if (stage >= 2 && Math.abs(moist.cumulusBaseFlux[i] - ph.CUMF[i]) > 1e-2 * Math.max(moist.cumulusBaseFlux[i], ph.CUMF[i])) plumesApart++;
    for (let k = 0; k < K; k++) {
      const x = k * C + i, t = mixed(k, i) ? below : above;
      const dT = Math.abs(theta[x] - after[1][x]), dQ = Math.abs(q[x] - after[4][x]), dQc = Math.abs(qc[x] - after[5][x]);
      t.n++;
      if (dT > 1e-3) { t.theta++; thetaApart.add(i); }
      if (dQ > 1e-7) t.q++;
      if (dQc > 1e-7) { t.qc++; qcApart.add(i); }
      t.worstT = Math.max(t.worstT, dT); t.worstQ = Math.max(t.worstQ, dQ); t.worstQc = Math.max(t.worstQc, dQc);
    }
  }
  const line = (name, t) => `${name} the mixing top, ${t.n} layers: θ apart in ${t.theta} (max ${t.worstT.toExponential(1)} K), q in ${t.q} (max ${t.worstQ.toExponential(1)}), qc in ${t.qc} (max ${t.worstQc.toExponential(1)})`;
  console.log(`${BC}, ${names[stage]}: ${thetaApart.size} of ${C} columns with θ apart, ${qcApart.size} with qc apart${stage >= 2 ? `, ${plumesApart} plumes' base flux apart by 1 %` : ''}\n  ${line('below', below)}\n  ${line('above', above)}`);
  if (DUMP >= 0) {
    const i = DUMP;
    console.log(`  cell ${i}: mixing top ${bl.mixingTop[i].toFixed(1)} m, cloud run bottom ${bl.cloudLayer[i]}, regime ${bl.regime[i]}, base flux ${moist.cumulusBaseFlux[i].toExponential(4)}/${ph.CUMF[i].toExponential(4)}`);
    for (let k = K - 1; k >= 0; k--) {
      const x = k * C + i;
      if (start[5][x] === 0 && qc[x] === 0 && after[5][x] === 0 && Math.abs(theta[x] - after[1][x]) < 1e-4) continue;
      console.log(`    k${k} ${(startGeo[x] / g).toFixed(0)} m ${mixed(k, i) ? 'mixed' : 'above'}: θ ${start[1][x].toFixed(4)} → ${theta[x].toFixed(5)}/${after[1][x].toFixed(5)}, q ${start[4][x].toExponential(4)} → ${q[x].toExponential(5)}/${after[4][x].toExponential(5)}, qc ${start[5][x].toExponential(3)} → ${qc[x].toExponential(4)}/${after[5][x].toExponential(4)}`);
    }
  }
}

const adjusted = model.state.map((a) => Float64Array.from(a, Math.fround));
model.state.forEach((a, n) => a.set(adjusted[n]));
bl.diagnose(model.state);
[bl.friction, bl.buoyancyFlux, radiation.longwave, radiation.stratiform, radiation.mlmTop].forEach(single);
const { gpu, put, run } = await gpuCore();
gpu.upload(adjusted);
gpu.uploadPhysics({ mlmGate: radiation.mlmGate, concentration: model.seaIce.concentration, land: landCode });
put('LWH', radiation.longwave); put('BUOY', bl.buoyancyFlux); put('USTAR', bl.friction); put('STRAT', radiation.stratiform);
put('MLMTOP', Float64Array.from(radiation.mlmTop, (v, i) => (v > 0 ? v - lift(i) : 0)));
run('pblDiagnose');
const ph = await gpu.downloadPhysics();
let cloudTopped = 0, regimes = 0, tops = 0, worstTop = 0, depths = 0, cooling = 0, entrainment = 0;
for (let i = 0; i < C; i++) {
  if (bl.regime[i] >= 2) cloudTopped++;
  if (bl.regime[i] !== ph.REGIME[i]) regimes++;
  const top = Math.abs(bl.mixingTop[i] - (ph.MIXTOP[i] + lift(i)));
  if (top > 1) tops++;
  worstTop = Math.max(worstTop, top);
  if (Math.abs(bl.depth[i] - (ph.DEPTH[i] + lift(i))) > 1) depths++;
  if (Math.abs(bl.cloudTopCooling[i] - ph.CTCOOL[i]) > 1e-3 * Math.max(1, Math.abs(bl.cloudTopCooling[i]))) cooling++;
  if (Math.abs(bl.entrainment[i] - ph.ENTRAIN[i]) > 1e-3 * Math.max(1e-3, bl.entrainment[i])) entrainment++;
}
console.log(`${BC}, the boundary layer's diagnosis of the adjusted state: ${cloudTopped} of ${C} columns cloud-topped; regime apart in ${regimes}, mixing top by over 1 m in ${tops} (largest ${worstTop.toFixed(1)} m), depth in ${depths}, cloud-top cooling in ${cooling}, entrainment by 1e-3 in ${entrainment}`);
process.exit(0);
