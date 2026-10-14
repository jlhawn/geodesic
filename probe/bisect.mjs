import { Grid } from '/Users/jlhawn/git_repos/jlhawn/geodesic/js/grid.module.js';
import { ADJUST_KERNELS } from '/Users/jlhawn/git_repos/jlhawn/geodesic/js/gpu/physics.gpu.js';
import { createModel } from '/Users/jlhawn/git_repos/jlhawn/geodesic/js/model.module.js';
import { createParallelModel } from '/Users/jlhawn/git_repos/jlhawn/geodesic/js/parallel.module.js';
import { createBoundaryLayer } from '/Users/jlhawn/git_repos/jlhawn/geodesic/js/physics/boundaryLayer.module.js';
import { initializeState } from '/Users/jlhawn/git_repos/jlhawn/geodesic/js/physics/init.module.js';
import { saturationHumidity } from '/Users/jlhawn/git_repos/jlhawn/geodesic/js/physics/moist.module.js';
import { sigmaInterfaces } from '/Users/jlhawn/git_repos/jlhawn/geodesic/js/dynamics/sigmaCore.module.js';
import { SEA_DRAG } from '/Users/jlhawn/git_repos/jlhawn/geodesic/js/physics/surface.module.js';
async function engines(entrainment) {
  const { createGpuCore } = await import('/Users/jlhawn/git_repos/jlhawn/geodesic/js/gpu/core.gpu.js');
  const levels = sigmaInterfaces('bl34');
  const pair = createModel(new Grid(6), { ocean: false, levels, boundaryLayer: { entrainment, turbulence: 'dry' }, surface: { exchange: 'fixed' } });
  const { core: c, mesh: m, state, radiation, moist, boundaryLayer: layer } = pair;
  const { K: nK, C: nC, E: nE, exnerLayer, sigmaMid, kappa, geopotential: phi, g: grav } = c.diagnostics;
  const [pi, theta, u, surfaceT, q, qc] = state;
  const init = initializeState(pair, {});
  for (let a = 0; a < init.length; a++) state[a].set(init[a]);
  c.diagnose(pi, theta, q, qc);
  let seed = 2287;
  const random = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; };
  for (let i = 0; i < nC; i++) {
    const zb = phi[(nK - 1) * nC + i] / grav, surface = 286 + 14 * random(), top = 300 + 1500 * random(), below = -1e-3 + 2e-3 * random(), above = 3e-3 + 5e-3 * random(), jump = 4 * random();
    const wetBelow = 0.4 + 0.3 * random(), wetAbove = 0.1 + 0.3 * random();
    for (let k = 0; k < nK; k++) {
      const x = k * nC + i, z = phi[x] / grav - zb;
      theta[x] = z < top ? surface + below * z : surface + below * top + jump + above * (z - top);
      if (sigmaMid[k] < 0.2) theta[x] = Math.max(theta[x], init[1][x]);
      q[x] = (z < top ? wetBelow : wetAbove) * saturationHumidity(theta[x] * exnerLayer[x], pi[i] * sigmaMid[k]);
      qc[x] = z < 3000 && random() < 0.1 ? 2e-4 * random() : 0;
    }
    surfaceT[i] = theta[(nK - 1) * nC + i] * exnerLayer[(nK - 1) * nC + i] * Math.pow(sigmaMid[nK - 1], -kappa) - 3 + 6 * random();
    const gate = random();
    radiation.mlmGate[i] = gate < 0.2 ? 0.7 : gate < 0.3 ? 0.5 : gate < 0.45 ? 0.5 + 0.1 * random() : 0.3 * random();
    const share = random();
    radiation.stratiform[i] = share < 0.1 ? 1 : share < 0.35 ? random() : 0;
  }
  for (let k = 0; k < nK; k++) for (let e = 0; e < nE; e++) u[k * nE + e] = k === nK - 1 ? 2 * (random() - 0.5) : 12 * (random() - 0.5);
  for (const a of state) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  for (const a of [radiation.mlmGate, radiation.stratiform]) for (let x = 0; x < a.length; x++) a[x] = Math.fround(a[x]);
  c.diagnose(pi, theta, q, qc);
  const gpu = await createGpuCore(m, { levels, physics: { entrainment, turbulence: 'dry', surfaceExchange: 'fixed' } });
  const { device, buffers, kernels } = gpu, dt = 900;
  gpu.upload(state);
  gpu.uploadPhysics({ mlmGate: radiation.mlmGate });
  device.queue.writeBuffer(buffers.PH, 4 * gpu.layout.PH.STRAT, Float32Array.from(radiation.stratiform));
  device.queue.writeBuffer(buffers.P, 0, Float32Array.from([dt, 0, 1, 0, 0, 0, 0, 0]));
  const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
  for (const [name, count] of (process.env.COEF ? [['pblDiagnose', nC], ['adjustMix', nC]] : [['pblDiagnose', nC], ...ADJUST_KERNELS.map((name) => [name, nC]), ['mixMomentum', nE]])) {
    pass.setPipeline(kernels[name]);
    pass.setBindGroup(0, device.createBindGroup({ layout: kernels[name].getBindGroupLayout(0), entries: [buffers.MI, buffers.MF, buffers.LV, buffers.S, buffers.K1, buffers.D, buffers.P, buffers.PH].map((buffer, binding) => ({ binding, resource: { buffer } })) }));
    pass.dispatchWorkgroups(Math.ceil(count / 64));
  }
  pass.end();
  device.queue.submit([encoder.finish()]);
  const after = await gpu.download(), ph = await gpu.downloadPhysics();
  const before = state.map((a) => Float64Array.from(a));
  layer.diagnose(state);
  const cpu = { entrainment: Float64Array.from(layer.entrainment), mixing: Float64Array.from(layer.mixing) };
  pair.phases.adjust(0, nC, dt);
  pair.phases.mixMomentum(0, nE, dt);
  return { K: nK, C: nC, E: nE, kTop: layer.kTop, gate: radiation.mlmGate, share: radiation.stratiform, before, cpu, state, after, ph };
}

const on = await engines({});
if (process.env.COEF) { await import('./coef.mjs'); globalThis.__coefHook(on); process.exit(0); }
const { K: nK, C: nC, E: nE, kTop } = on;
const names = ['pi', 'theta', 'u', 'Ts', 'q', 'qc'];
for (const a of [1, 4, 5]) {
  const byLayer = Array.from({ length: nK }, () => 0); let cols = 0, worst = 0, wi = -1;
  for (let i = 0; i < nC; i++) { let m = 0; for (let k = 0; k < nK; k++) { const x = k * nC + i, d = Math.abs(on.state[a][x] - on.after[a][x]); byLayer[k] = Math.max(byLayer[k], d); m = Math.max(m, d); } if (m > (a === 1 ? 1e-3 : 1e-6)) cols++; if (m > worst) { worst = m; wi = i; } }
  console.log(`${names[a]}: ${cols} of ${nC} columns differ; worst ${worst.toExponential(2)} at column ${wi}; by layer (kTop=${kTop}): ${byLayer.map((v, k) => k >= kTop - 1 ? v.toExponential(0) : '').filter(Boolean).join(' ')}`);
  if (a === 1 || a === 4) {
    const i = wi; const wet = Array.from({ length: nK }, (_, k) => on.before[5][k * nC + i] > 0).some(Boolean);
    console.log(`  column ${i}: wet ${wet}, gate ${on.gate[i].toFixed(2)}, share ${on.share[i].toFixed(2)}, w_e cpu ${on.cpu.entrainment[i].toExponential(2)} gpu ${on.ph.ENTRAIN[i].toExponential(2)}`);
    const rows = [];
    for (let k = kTop; k < nK; k++) { const x = k * nC + i; rows.push(`k${k} before ${on.before[a][x].toFixed(a === 1 ? 2 : 5)} cpu ${on.state[a][x].toFixed(a === 1 ? 2 : 5)} gpu ${on.after[a][x].toFixed(a === 1 ? 2 : 5)} mix ${on.cpu.mixing[k * nC + i].toExponential(1)}/${on.ph.MIX[k * nC + i].toExponential(1)}`); }
    console.log('  ' + rows.join('\n  '));
  }
}
let wind = 0, we = -1; for (let x = 0; x < nK * nE; x++) { const d = Math.abs(on.state[2][x] - on.after[2][x]); if (d > wind) { wind = d; we = x; } }
console.log(`wind: worst ${wind.toExponential(2)} at k ${Math.floor(we / nE)} edge ${we % nE}; edges differing > 1e-3: ${Array.from({ length: nK * nE }, (_, x) => Math.abs(on.state[2][x] - on.after[2][x]) > 1e-3).filter(Boolean).length} of ${nK * nE}`);
process.exit(0);
