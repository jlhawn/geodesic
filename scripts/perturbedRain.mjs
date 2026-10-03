// How far the CPU model's rain moves under θ perturbations near the size
// of the engines' single-precision differences, on the setup of the rain
// accumulation parity test in test/gpuModel.test.mjs (24 steps of 900 s at
// N=6). One run is plain, the other has ±AMP K of uniform noise added to
// every layer's θ (from the generator seeded with SEED) before the first
// step and, unless EACH=0, after every step; the count of parted cells
// depends on the seed, so compare settings over several. It prints the
// cells whose convective or large-scale rain parts by more than 1e-3 of
// the largest cell's (the test's outlier rule) and, for
// each, the first discrete decision that parted between the runs (the
// regime, whether the plume fires, its base flux by 1 % or its top, a
// layer's cloud against the cloud-top threshold, or a merge of the dry
// adjustment, by the rule of test/helpers/decisions.mjs that the parity
// tests use), and how many cells those decisions cover with their
// neighbours within one and two cells.
//   node scripts/perturbedRain.mjs
//   SEED=2024 node scripts/perturbedRain.mjs
//   BC=cloudLayer AMP=1e-3 EACH=0 node scripts/perturbedRain.mjs
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { decisionTracker, cpuDecisions, neighbourhood, DECISION_KINDS } from '../test/helpers/decisions.mjs';

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
const { K } = a.core, C = a.mesh.nCells;
const decisions = decisionTracker(K, C);
for (let n = 1; n <= STEPS; n++) {
  a.step(900); b.step(900);
  if (EACH) perturb();
  decisions.check(cpuDecisions(a), cpuDecisions(b));
}
const convective = [a.moist.convectivePrecipitation, b.moist.convectivePrecipitation], largeScale = [a.moist.largeScalePrecipitation, b.moist.largeScalePrecipitation];
const largest = Math.max(...convective[0]), largestScale = Math.max(...largeScale[0]);
const apart = [];
for (let i = 0; i < C; i++) if (Math.abs(convective[0][i] - convective[1][i]) > 1e-3 * largest || Math.abs(largeScale[0][i] - largeScale[1][i]) > 1e-3 * largestScale) apart.push(i);
console.log(`${BC}, ±${AMP} K of θ ${EACH ? 'before every step' : 'once'}, ${STEPS} steps at N=6: ${apart.length} of ${C} cells' rain apart by 1e-3 of the largest (convective ${largest.toFixed(3)}, large-scale ${largestScale.toFixed(3)} kg/m²); a discrete decision parted in ${decisions.parted.size} cells (${DECISION_KINDS.map((kind) => `${kind} ${decisions.kinds[kind].size}`).join(', ')}), ${neighbourhood(a.mesh, decisions.parted, 1).size} with their neighbours, ${neighbourhood(a.mesh, decisions.parted, 2).size} with those within two cells; ${apart.filter((i) => !neighbourhood(a.mesh, decisions.parted, 2).has(i)).length} of the cells apart lie outside those`);
for (const i of apart) console.log(`  cell ${i}: convective ${(Math.abs(convective[0][i] - convective[1][i]) / largest).toExponential(1)} and large-scale ${(Math.abs(largeScale[0][i] - largeScale[1][i]) / largestScale).toExponential(1)} of the largest; ${decisions.first.get(i) ?? 'no discrete decision parted'}`);
