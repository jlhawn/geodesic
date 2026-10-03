// How many columns' discrete decisions the CPU model takes apart from
// itself under θ perturbations near the size of the engines'
// single-precision differences, on the setup of a parity test: SETUP
// 'cloudEffect' (test/cloudEffect.test.mjs, 16 steps), 'dayMeans'
// (test/dayMeans.test.mjs, 24 steps) or 'treeline'
// (test/treeline.test.mjs, 48 steps). One run is plain, the other has ±AMP
// K of uniform noise on every layer's θ (seeded with SEED) before the first
// step and after every step. It counts the columns whose decisions parted
// (test/helpers/decisions.mjs and the tests' own: OLR or absorbed sunlight
// by 1 W/m² at a step; for the cloud effects the bottom block of uniform q
// and the layers holding cloud water; on land whether the step is in
// season and whether there is snow) and their neighbours within one and
// two cells, and prints the test's measures between the runs over every
// column and over the columns kept: for the treeline those outside two
// cells of a land decision that parted, for the others outside two cells
// of any.
//   SETUP=cloudEffect node scripts/perturbedDecisions.mjs
//   SETUP=treeline SEED=7 node scripts/perturbedDecisions.mjs
import { Grid } from '../js/grid.module.js';
import { createModel } from '../js/model.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { syntheticTopography } from '../js/geography.module.js';
import { decisionTracker, cpuDecisions, neighbourhood } from '../test/helpers/decisions.mjs';

const SETUP = process.env.SETUP ?? 'cloudEffect', AMP = Number(process.env.AMP ?? 1e-4);
let seed = Number(process.env.SEED ?? 99);
const random = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648 - 0.5; };

const UNSCATTERED = { rayleighDepth: 0, nearInfraredRayleigh: 0, landAerosol: 0, seaAerosol: 0, skylight: 0.15, upwardAbsorption: false, longwaveScheme: 'gray', solarGases: 'lacisHansen' };
const SETUPS = {
  cloudEffect: {
    steps: 16,
    make() {
      const model = createModel(new Grid(6), { ocean: false, radiation: { clearSkyPass: true }, surface: { exchange: 'fixed' }, moist: { excessVelocity: 'convective', plumeEntrainmentLaw: 'gregory' } });
      const init = initializeState(model, {}), { K, C, sigmaMid } = model.core.diagnostics;
      for (let k = 0; k < K; k++) if (sigmaMid[k] > 0.5 && sigmaMid[k] < 0.9) for (let i = 0; i < C; i++) init[5][k * C + i] = 3e-4;
      init.forEach((values, a) => model.state[a].set(values));
      return model;
    },
  },
  dayMeans: {
    steps: 24,
    make() {
      const model = createModel(new Grid(6), { ocean: false, radiation: UNSCATTERED, surface: { exchange: 'fixed' }, moist: { excessVelocity: 'convective', convectionType: 'cloudDepth', plumePhase: 'liquid', plumeConversion: 'zhangMcFarlane' } });
      initializeState(model, {}).forEach((values, a) => model.state[a].set(values));
      return model;
    },
  },
  treeline: {
    steps: 48,
    make() {
      const topography = syntheticTopography(90, 180, (lat, lon) => ((Math.cos(lon) > 0 && Math.abs(lat) < 1.2) || lat < -1.15 ? 300 : -4000));
      const land = { soilCarbon: false, seasonMemory: 6 * 3600, treeGrowthTime: 3 * 3600, treeDeclineTime: 2 * 3600, treelineWarmth: [6, 22], growthTime: 4 * 3600, declineTime: 3 * 3600, snowDeclineTime: 4 * 3600 };
      const model = createModel(new Grid(6), { topography, land, surface: { exchange: 'fixed' } });
      const C = model.mesh.nCells, init = initializeState(model, {});
      init.forEach((values, a) => model.state[a].set(values));
      for (let i = 0; i < C; i++) if (model.geography.land[i]) model.state[6][i] = 0;
      model.seaIce.load(model.state[6], null);
      model.ocean.initialize(model.state[3], model.state[6]);
      let s = 3;
      const rnd = () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; };
      const snow = Float64Array.from({ length: C }, (_, i) => (Math.abs(model.mesh.latCell[i]) > 0.9 ? 30 : 0));
      const vegetation = Float64Array.from({ length: C }, () => rnd()), canopy = Float64Array.from({ length: C }, () => rnd());
      const seasonLength = Float64Array.from({ length: C }, () => rnd()), seasonWarmth = Float64Array.from({ length: C }, () => 10 * rnd());
      model.land.load({ soil: Float64Array.from({ length: C }, () => 300 * rnd()), snow, vegetation, canopy, seasonLength, seasonWarmth }, model.state[6]);
      return model;
    },
  },
};

const setup = SETUPS[SETUP];
const a = setup.make(), b = setup.make(), { K } = a.core, C = a.mesh.nCells, mesh = a.mesh;
const perturb = () => { for (let x = 0; x < b.state[1].length; x++) b.state[1][x] += 2 * AMP * random(); };
perturb();
const decisions = decisionTracker(K, C), own = { merged: new Set(), condensed: new Set(), radiation: new Set(), season: new Set(), snow: new Set() };
const mixedTop = (q, i) => { let k = K - 1; while (k > 0 && Math.abs(q[(k - 1) * C + i] - q[(K - 1) * C + i]) <= 1e-6 * q[(K - 1) * C + i]) k--; return k; };
const cloudyLayers = (qc, i) => { let bits = ''; for (let k = 0; k < K; k++) bits += qc[k * C + i] > 0 ? '1' : '0'; return bits; };
const keep = 900 / (6 * 3600), inSeason = (after, before, i) => Math.round((after[i] - (1 - keep) * before[i]) / keep);
const { exnerLayer } = a.core.diagnostics, airGap = new Float64Array(C);
for (let n = 0; n < setup.steps; n++) {
  const absorbed = [Float64Array.from(a.radiation.summed.absorbedSolar), Float64Array.from(b.radiation.summed.absorbedSolar)];
  const seasons = SETUP === 'treeline' ? [Float64Array.from(a.land.seasonLength), Float64Array.from(b.land.seasonLength)] : null;
  a.step(900); b.step(900);
  perturb();
  decisions.check(cpuDecisions(a), cpuDecisions(b));
  for (let i = 0; i < C; i++) {
    if (SETUP === 'cloudEffect' && mixedTop(a.state[4], i) !== mixedTop(b.state[4], i)) own.merged.add(i);
    if (SETUP === 'cloudEffect' && cloudyLayers(a.state[5], i) !== cloudyLayers(b.state[5], i)) own.condensed.add(i);
    if (Math.abs(a.radiation.outgoing[i] - b.radiation.outgoing[i]) > 1 || Math.abs(a.radiation.summed.absorbedSolar[i] - absorbed[0][i] - (b.radiation.summed.absorbedSolar[i] - absorbed[1][i])) > 1) own.radiation.add(i);
    if (seasons && a.geography.land[i]) {
      const x = (K - 1) * C + i;
      airGap[i] = Math.max(airGap[i], Math.abs(a.state[1][x] - b.state[1][x]) * exnerLayer[x]);
      if (inSeason(a.land.seasonLength, seasons[0], i) !== inSeason(b.land.seasonLength, seasons[1], i)) own.season.add(i);
      if ((a.land.snow[i] > 0) !== (b.land.snow[i] > 0)) own.snow.add(i);
    }
  }
}
const parted = new Set([...decisions.parted, ...Object.values(own).flatMap((cells) => [...cells])]);
const near1 = neighbourhood(mesh, parted, 1), near2 = neighbourhood(mesh, parted, 2);
console.log(`${SETUP}, ±${AMP} K of θ before every step, ${setup.steps} steps at N=6: a decision parted in ${parted.size} of ${C} columns (${Object.entries(decisions.kinds).map(([kind, cells]) => `${kind} ${cells.size}`).join(', ')}; ${Object.entries(own).map(([kind, cells]) => `${kind} ${cells.size}`).join(', ')}), ${near1.size} with their neighbours, ${near2.size} within two cells`);

const rmsRel = (x, y, cells) => { let d = 0, r = 0; for (const i of cells) { d += (x[i] - y[i]) ** 2; r += x[i] ** 2; } return Math.sqrt(d / Math.max(r, 1e-300)); };
const maxDiff = (x, y, cells) => cells.reduce((m, i) => Math.max(m, Math.abs(x[i] - y[i])), 0);
const all = Array.from({ length: C }, (_, i) => i), kept = all.filter((i) => !near2.has(i));
const meanOver = (values, cells) => { let s = 0, w = 0; for (const i of cells) { s += mesh.areaCell[i] * values[i]; w += mesh.areaCell[i]; } return s / w; };
if (SETUP === 'dayMeans') {
  const sums = ['absorbedSolar', 'reflectedSolar', 'outgoingLongwave'].map((name) => [name, Float64Array.from(a.radiation.summed[name]), Float64Array.from(b.radiation.summed[name])]);
  const da = a.diagnostics(), db = b.diagnostics();
  console.log(`  day means over every column: ASR ${(Math.abs(db.absorbedSolar - da.absorbedSolar) / da.absorbedSolar).toExponential(1)} of itself, OLR ${(Math.abs(db.outgoingLongwave - da.outgoingLongwave) / da.outgoingLongwave).toExponential(1)}, albedo ${Math.abs(db.planetaryAlbedo - da.planetaryAlbedo).toExponential(1)}`);
  for (const [label, cells] of [['every column', all], [`the ${kept.length} kept`, kept]]) {
    const asr = [meanOver(a.radiation.meanAbsorbedSolar, cells), meanOver(b.radiation.meanAbsorbedSolar, cells)], olr = [meanOver(a.radiation.meanOutgoingLongwave, cells), meanOver(b.radiation.meanOutgoingLongwave, cells)];
    console.log(`  over ${label}: per-cell sums ${sums.map(([name, x, y]) => `${name} rms ${rmsRel(x, y, cells).toExponential(1)}`).join(', ')}; area mean of the per-cell ASR ${(Math.abs(asr[1] - asr[0]) / asr[0]).toExponential(1)} of itself, OLR ${(Math.abs(olr[1] - olr[0]) / olr[0]).toExponential(1)}`);
  }
}
if (SETUP === 'cloudEffect') {
  a.diagnostics(); b.diagnostics();
  for (const [label, cells] of [['every column', all], [`the ${kept.length} kept`, kept]]) console.log(`  over ${label}: ${['meanShortwaveCloudEffect', 'meanLongwaveCloudEffect'].map((name) => `${name} rms ${rmsRel(a.radiation[name], b.radiation[name], cells).toExponential(1)} (max ${maxDiff(a.radiation[name], b.radiation[name], cells).toExponential(1)} W/m²)`).join(', ')}`);
}
if (SETUP === 'treeline') {
  const landCells = all.filter((i) => a.geography.land[i] && !a.geography.iceSheet[i]), landParted = new Set([...own.season, ...own.snow]);
  const landNear = neighbourhood(mesh, landParted, 2), landKept = landCells.filter((i) => !landNear.has(i));
  console.log(`  the land's own decisions (in season, under snow) parted in ${landParted.size} land cells, ${landCells.length - landKept.length} of ${landCells.length} with the land cells within two cells; the lowest air apart by up to ${Math.max(...landCells.map((i) => airGap[i])).toExponential(1)} K, by more than 1e-3 K on ${landCells.filter((i) => airGap[i] > 1e-3).length} land cells`);
  for (const [label, cells] of [[`all ${landCells.length} land cells`, landCells], [`the ${landKept.length} land cells kept`, landKept]]) console.log(`  over ${label}: season length ${maxDiff(a.land.seasonLength, b.land.seasonLength, cells).toExponential(1)}, warmth ${maxDiff(a.land.seasonWarmth, b.land.seasonWarmth, cells).toExponential(1)} K, trees ${maxDiff(a.land.canopy, b.land.canopy, cells).toExponential(1)}`);
}
