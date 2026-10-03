// How far the CPU model's rain moves under θ perturbations near the size
// of the engines' single-precision differences, on the setup of the rain
// accumulation parity test in test/gpuModel.test.mjs (24 steps of 900 s at
// N=6). One run is plain, the other has ±AMP K of uniform noise added to
// every layer's θ (from the generator seeded with SEED) before the first
// step and, unless EACH=0, after every step; the count of parted cells
// depends on the seed, so compare settings over several. It prints the
// cells whose convective or large-scale rain parts by more than 1e-3 of
// the largest cell's (the test's outlier rule) and, for
// each, the first discrete decision that parted between the runs: the
// regime, the cloud-top run's lowest layer, the count of mixed interfaces,
// whether the plume fires, its base flux by 1 % or its top, a layer's
// cloud against the cloud-top threshold, or a merge of the dry adjustment.
//   node scripts/perturbedRain.mjs
//   SEED=2024 node scripts/perturbedRain.mjs
//   BC=cloudLayer AMP=1e-3 EACH=0 node scripts/perturbedRain.mjs
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { CLOUD_TOP_DEFAULTS } from '../js/physics/boundaryLayer.module.js';

const BC = process.env.BC ?? 'uniform', STEPS = Number(process.env.STEPS ?? 24), AMP = Number(process.env.AMP ?? 1e-4), EACH = process.env.EACH !== '0';
const moist = { rainEvaporation: 0, excessVelocity: 'convective', plumePhase: 'liquid', plumeConversion: 'zhangMcFarlane', boundaryCondensation: BC };
const make = () => createModel(new Grid(6), { ocean: false, radiation: { stratus: false, surfaceExchange: 'fixed' }, moist, surface: { exchange: 'fixed' } });
const a = make(), b = make();
const init = initializeState(a, {});
for (let n = 0; n < init.length; n++) { a.state[n].set(init[n]); b.state[n].set(init[n]); }
let seed = Number(process.env.SEED ?? 99);
const random = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648 - 0.5; };
const perturb = () => { for (let x = 0; x < b.state[1].length; x++) b.state[1][x] += 2 * AMP * random(); };
perturb();
const { K } = a.core, C = a.mesh.nCells, threshold = CLOUD_TOP_DEFAULTS.threshold;
const decisions = new Map();
const note = (i, step, what) => { if (!decisions.has(i)) decisions.set(i, `step ${step}: ${what}`); };
const mixedInterfaces = (m, i) => { let n = 0; for (let k = 0; k < K; k++) if (m.boundaryLayer.mixing[k * C + i] > 0) n++; return n; };
for (let n = 1; n <= STEPS; n++) {
  a.step(900); b.step(900);
  if (EACH) perturb();
  for (let i = 0; i < C; i++) {
    const A = a.boundaryLayer, B = b.boundaryLayer, fa = a.moist.cumulusBaseFlux[i], fb = b.moist.cumulusBaseFlux[i];
    if (A.regime[i] !== B.regime[i]) note(i, n, `regime ${A.regime[i]}/${B.regime[i]}`);
    if (A.cloudLayer[i] !== B.cloudLayer[i]) note(i, n, `cloud-top run to layer ${A.cloudLayer[i]}/${B.cloudLayer[i]}`);
    if (mixedInterfaces(a, i) !== mixedInterfaces(b, i)) note(i, n, `mixed interfaces ${mixedInterfaces(a, i)}/${mixedInterfaces(b, i)}`);
    if ((fa > 0) !== (fb > 0)) note(i, n, `the plume fires ${fa > 0}/${fb > 0}`);
    else if (fa > 0 && Math.abs(fa - fb) > 1e-2 * Math.max(fa, fb)) note(i, n, `plume base flux ${fa.toExponential(3)}/${fb.toExponential(3)} kg/m²/s`);
    else if (fa > 0 && Math.abs(a.moist.cumulusTop[i] - b.moist.cumulusTop[i]) > 1) note(i, n, `plume top ${(a.moist.cumulusTop[i] / 100).toFixed(0)}/${(b.moist.cumulusTop[i] / 100).toFixed(0)} hPa`);
    for (let k = 0; k < K; k++) { const x = k * C + i; if ((a.state[5][x] > threshold) !== (b.state[5][x] > threshold)) { note(i, n, `cloud in layer ${k} ${a.state[5][x].toExponential(2)}/${b.state[5][x].toExponential(2)}`); break; } }
    for (let k = 0; k < K - 1; k++) {
      const x = k * C + i, y = x + C;
      if ((Math.abs(a.state[1][x] - a.state[1][y]) <= 1e-6 * a.state[1][y]) !== (Math.abs(b.state[1][x] - b.state[1][y]) <= 1e-6 * b.state[1][y])) { note(i, n, `dry adjustment merges layer ${k}`); break; }
    }
  }
}
const convective = [a.moist.convectivePrecipitation, b.moist.convectivePrecipitation], largeScale = [a.moist.largeScalePrecipitation, b.moist.largeScalePrecipitation];
const largest = Math.max(...convective[0]), largestScale = Math.max(...largeScale[0]);
const apart = [];
for (let i = 0; i < C; i++) if (Math.abs(convective[0][i] - convective[1][i]) > 1e-3 * largest || Math.abs(largeScale[0][i] - largeScale[1][i]) > 1e-3 * largestScale) apart.push(i);
console.log(`${BC}, ±${AMP} K of θ ${EACH ? 'before every step' : 'once'}, ${STEPS} steps at N=6: ${apart.length} of ${C} cells' rain apart by 1e-3 of the largest (convective ${largest.toFixed(3)}, large-scale ${largestScale.toFixed(3)} kg/m²); a discrete decision parted in ${decisions.size} cells`);
for (const i of apart) console.log(`  cell ${i}: convective ${(Math.abs(convective[0][i] - convective[1][i]) / largest).toExponential(1)} and large-scale ${(Math.abs(largeScale[0][i] - largeScale[1][i]) / largestScale).toExponential(1)} of the largest; ${decisions.get(i) ?? 'no discrete decision parted'}`);
