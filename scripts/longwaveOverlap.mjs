// The longwave of a saved state under three cloud overlaps, on the
// single-thread CPU engine with the ocean off:
//   node scripts/longwaveOverlap.mjs <state.bin>
// One model step is taken first (the plumes' cumulus is not saved); then
// every column's radiation runs once at the state's time (dt 0, so the deck's
// carried state does not move) and its longwave is recomputed from the same
// gases (the correlated g-points), temperatures, cloud water, covers and
// infrared optics in three ways: as the radiation does it, each layer's cover
// entering through its emissivity f (1 - exp(-kappa W/f)), which is the
// expectation under random overlap; and as the expectation over the two-state
// Markov chain of exponential-random overlap (adjacent layers overlapping with
// alpha = exp(-dz/z0), z0 the radiation's decorrelation length, the chain of
// the radiation's own shortwave cover) and of maximum-random overlap
// (alpha 1), each layer cloudy over f with its in-cloud emissivity
// 1 - exp(-kappa W/f) or clear. The deck's layer keeps the radiation's own
// emissivity, outside the chain. Printed: how closely the random recomputation
// repeats the radiation's OLR, and the OLR, the longwave cloud effect and the
// downward longwave at the surface under each overlap, global, 30S-30N and in
// the tropical boxes of scripts/tropicalHeating.mjs.
// RADIATION, MOIST, BOUNDARY_LAYER and SURFACE (JSON) pass options as to
// scripts/verticalAudit.mjs; the radiation must use longwaveScheme
// 'correlated'.
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { createModel, STATE_NAMES } from '../js/model.module.js';
import { decodeState, savedLevels } from '../js/stateFile.module.js';
import { savedDeckField, DECK_FIELDS } from '../js/physics/regrid.module.js';
import { DECORRELATION_LENGTH, DECORRELATION_SLOPE, GREENHOUSE_GASES, OZONE_COLUMN, STEFAN_BOLTZMANN, YEAR } from '../js/physics/radiation.module.js';
import { ozoneWeights, ozoneAbove as climatologyAbove } from '../js/physics/ozone.module.js';
import { LONGWAVE_TABLE, LONGWAVE_CONSTANTS, GAS_MOLAR, layerPaths, planckShare } from '../js/physics/longwave.module.js';
import { OZONE_CM_ATM } from '../js/physics/shortwaveGases.module.js';
import { SEA_DRAG, LAND_DRAG } from '../js/physics/surface.module.js';
import { FREEZING_POINT } from '../js/physics/ice.module.js';
import { BOXES, inLongitudes } from '../js/audit.module.js';

const FILE = process.argv[2];
if (!FILE) { console.error('usage: node scripts/longwaveOverlap.mjs <state.bin>'); process.exit(1); }
const RADIATION = JSON.parse(process.env.RADIATION ?? '{}'), MOIST = JSON.parse(process.env.MOIST ?? '{}'), BOUNDARY_LAYER = JSON.parse(process.env.BOUNDARY_LAYER ?? '{}'), SURFACE = JSON.parse(process.env.SURFACE ?? '{}');
if ((RADIATION.longwaveScheme ?? 'correlated') !== 'correlated') throw new Error('the recomputation follows the correlated longwave');
const t0 = performance.now();
const say = (s = '') => console.log(s);

const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const model = createModel(new Grid(saved.N), { topography, levels: savedLevels(saved), ocean: false, radiation: { clearSkyPass: true, ...RADIATION }, moist: MOIST, boundaryLayer: BOUNDARY_LAYER, surface: SURFACE });
const { mesh, core, state, radiation, boundaryLayer: bl, moist, seaIce, land, surface } = model;
const { K, C, levels, sigmaMid, dSigma, g, exnerLayer, geopotential } = core.diagnostics;
STATE_NAMES.forEach((name, a) => state[a].set(saved[name]));
seaIce.load(state[6], saved.concentration ?? null);
for (const field of Object.keys(DECK_FIELDS)) radiation[field].set(savedDeckField(saved, field, model));
land.load(saved.land, state[6]);
model.time = saved.time;
const zs = Float64Array.from({ length: C }, (_, i) => (model.surfaceGeopotential ? model.surfaceGeopotential[i] / g : 0));
if (saved.boundaryDepth) bl.depth.set(Float64Array.from(saved.boundaryDepth, (z, i) => z + zs[i]));
if (saved.mixingTop) bl.mixingTop.set(Float64Array.from(saved.mixingTop, (z, i) => z + zs[i]));
if (saved.boundaryRegime) bl.regime.set(saved.boundaryRegime);
if (saved.boundaryBuoyancy) bl.buoyancyFlux.set(saved.boundaryBuoyancy);
const dt = 1350 * 16 / saved.N, bottom = K - 1;
model.step(dt);
for (let i = 0; i < C; i++) core.diagnoseColumn(i, state[0], state[1], state[4], state[5]);
radiation.setTime(model.time);
surface.lowestWindSpeed(state[2]);

const [pi, theta, , surfaceT, q, qc, ice] = state;
const deg = 180 / Math.PI, area = mesh.areaCell, landMask = model.geography.land;
const lat = Float64Array.from(mesh.latCell, (x) => x * deg), lon = Float64Array.from(mesh.lonCell, (x) => x * deg);
const R = { ...{ carbonDioxide: GREENHOUSE_GASES.carbonDioxide, methane: GREENHOUSE_GASES.methane, nitrousOxide: GREENHOUSE_GASES.nitrousOxide, ozone: 'afgl', ozoneColumn: OZONE_COLUMN, ozoneHeight: 25e3, ozoneWidth: 5e3, scaleHeight: 7e3, decorrelationLength: DECORRELATION_LENGTH, decorrelationSlope: DECORRELATION_SLOPE, gustiness: 3, cumulusCloud: true }, ...RADIATION };
const wellMixed = [R.carbonDioxide * GAS_MOLAR.co2 / GAS_MOLAR.air, R.methane * GAS_MOLAR.ch4 / GAS_MOLAR.air, R.nitrousOxide * GAS_MOLAR.n2o / GAS_MOLAR.air];
const ozoneAbove = (sigma) => (sigma <= 0 ? 0 : (1 + Math.exp(-R.ozoneHeight / R.ozoneWidth)) / (1 + Math.exp((-R.scaleHeight * Math.log(sigma) - R.ozoneHeight) / R.ozoneWidth)));
const ozoneShare = Float64Array.from({ length: K }, (_, k) => ozoneAbove(levels[k + 1]) - ozoneAbove(levels[k]));
const layerOzone = new Float64Array(K), weights = new Float64Array(5);
const points = LONGWAVE_TABLE.points, NG = points.length;
const paths = Array.from({ length: 6 }, () => new Float64Array(K)), pathRow = new Float64Array(6);
const T = new Float64Array(K), water = new Float64Array(K), cover = new Float64Array(K), inCloud = new Float64Array(K), effective = new Float64Array(K), outside = new Float64Array(K), gas = new Float64Array(NG * K), planck = new Float64Array(NG * K), surfaceUp = new Float64Array(NG);
const alpha = new Float64Array(K);

// Expected upward flux at the top and downward flux at the surface of one
// g-point over the Markov chain whose adjacent layers overlap with alpha
// (alpha[k] between layers k - 1 and k); `random` takes the effective
// emissivities instead.
function expected(gp, overlap) {
  const E = (k, s) => { const eg = gas[gp * K + k]; return s ? 1 - (1 - eg) * (1 - inCloud[k]) : 1 - (1 - eg) * (1 - outside[k]); };
  const joint = (k) => {
    const a = cover[k - 1], b = cover[k], al = overlap === 'maximum' ? 1 : overlap === 'random' ? 0 : alpha[k];
    const pair = al * Math.max(a, b) + (1 - al) * (a + b - a * b), both = a + b - pair;
    return [1 - a - b + both, b - both, a - both, both];
  };
  let up0 = 0, up1 = 0;
  {
    const k = bottom, f = cover[k];
    up0 = (1 - f) * (surfaceUp[gp] * (1 - E(k, 0)) + E(k, 0) * planck[gp * K + k]);
    up1 = f * (surfaceUp[gp] * (1 - E(k, 1)) + E(k, 1) * planck[gp * K + k]);
  }
  for (let k = bottom - 1; k >= 0; k--) {
    const [p00, p01, p10, p11] = joint(k + 1), f = cover[k], fb = cover[k + 1];
    const from0 = fb < 1 ? up0 / (1 - fb) : 0, from1 = fb > 0 ? up1 / fb : 0;
    const in0 = from0 * p00 + from1 * p01, in1 = from0 * p10 + from1 * p11;
    up0 = in0 * (1 - E(k, 0)) + E(k, 0) * planck[gp * K + k] * (1 - f);
    up1 = in1 * (1 - E(k, 1)) + E(k, 1) * planck[gp * K + k] * f;
  }
  let down0 = 0, down1 = 0;
  {
    const f = cover[0];
    down0 = (1 - f) * E(0, 0) * planck[gp * K];
    down1 = f * E(0, 1) * planck[gp * K];
  }
  for (let k = 1; k < K; k++) {
    const [p00, p01, p10, p11] = joint(k), f = cover[k], fa = cover[k - 1];
    const from0 = fa < 1 ? down0 / (1 - fa) : 0, from1 = fa > 0 ? down1 / fa : 0;
    const in0 = from0 * p00 + from1 * p10, in1 = from0 * p01 + from1 * p11;
    down0 = in0 * (1 - E(k, 0)) + E(k, 0) * planck[gp * K + k] * (1 - f);
    down1 = in1 * (1 - E(k, 1)) + E(k, 1) * planck[gp * K + k] * f;
  }
  return [up0 + up1, down0 + down1];
}
function randomFlux(gp) {
  let up = surfaceUp[gp], down = 0;
  for (let k = bottom; k >= 0; k--) { const e = 1 - (1 - gas[gp * K + k]) * (1 - effective[k]); up = up * (1 - e) + e * planck[gp * K + k]; }
  for (let k = 0; k < K; k++) { const e = 1 - (1 - gas[gp * K + k]) * (1 - effective[k]); down = down * (1 - e) + e * planck[gp * K + k]; }
  return [up, down];
}

const NAMES = ['random (the radiation)', 'exponential-random', 'maximum-random'];
const out = { olr: NAMES.map(() => new Float64Array(C)), down: NAMES.map(() => new Float64Array(C)), model: new Float64Array(C), modelDown: new Float64Array(C), clear: new Float64Array(C) };
let worst = 0, worstDown = 0, worstChain = 0, decks = 0;
for (let i = 0; i < C; i++) {
  let skin = surfaceT[i], openSea = 0;
  const h = ice[i];
  if (!landMask[i]) { const cov = seaIce.cover(i, h); openSea = 1 - cov; if (h > 0 && cov < 1) skin = cov * surfaceT[i] + (1 - cov) * FREEZING_POINT; }
  const drag = landMask[i] ? LAND_DRAG : SURFACE.dragCoefficient ?? SEA_DRAG;
  const albedo = model.surfaceAlbedo[i];
  radiation.column(i, pi[i], theta, skin, surface.windSpeed[i], undefined, radiation.insolation(i), q[bottom * C + i], q, qc, albedo, albedo, 1, drag, openSea, bl.depth[i] - geopotential[bottom * C + i] / g, 0, bl.mixingTop[i] - geopotential[bottom * C + i] / g);
  const budget = radiation.budget;
  out.model[i] = budget.outgoingLongwave; out.modelDown[i] = budget.downwardLongwave; out.clear[i] = budget.clearOutgoingLongwave;
  const ozoneCell = R.ozoneColumn[0] + (R.ozoneColumn[1] - R.ozoneColumn[0]) * Math.sin(mesh.latCell[i]) ** 2;
  if (R.ozone === 'afgl') { ozoneWeights(mesh.latCell[i], (model.time % YEAR) / YEAR, weights); let above = 0; for (let k = 0; k < K; k++) { const below = climatologyAbove(pi[i] * levels[k + 1], weights); layerOzone[k] = below - above; above = below; } }
  else for (let k = 0; k < K; k++) layerOzone[k] = ozoneCell * ozoneShare[k];
  const z0 = R.decorrelationLength - R.decorrelationSlope * Math.abs(lat[i]);
  const deck = budget.stratus, fraction = budget.stratusFraction;
  if (deck > 0) decks++;
  for (let k = 0; k < K; k++) {
    const idx = k * C + i, mass = pi[i] * dSigma[k] / g, kappa = radiation.infrared[k];
    T[k] = theta[idx] * exnerLayer[idx];
    water[k] = Math.max(0, qc[idx]) * mass;
    const cumulus = R.cumulusCloud ? moist.cumulusCover[idx] * moist.cumulusWater[idx] * mass : 0;
    if (cumulus > 0) water[k] += cumulus;
    cover[k] = water[k] > 0 ? radiation.layerCover[k] : 0;
    inCloud[k] = water[k] > 0 ? 1 - Math.exp(-kappa * water[k] / cover[k]) : 0;
    effective[k] = cover[k] * inCloud[k];
    outside[k] = 0;
    if (deck > 0 && k === radiation.stratusLayer) {
      effective[k] = fraction * (1 - Math.exp(-kappa * (water[k] + deck))) + (1 - fraction) * effective[k];
      outside[k] = effective[k]; cover[k] = 0; inCloud[k] = 0;
    }
    alpha[k] = k > 0 ? Math.exp(-(geopotential[(k - 1) * C + i] - geopotential[idx]) / (g * z0)) : 0;
    const dry = Math.max(0, 1 - Math.max(0, q[idx]));
    layerPaths(pathRow, pi[i] * sigmaMid[k], mass, T[k], q[idx], layerOzone[k] * OZONE_CM_ATM, wellMixed[0] * dry, wellMixed[1] * dry, wellMixed[2] * dry);
    for (let j = 0; j < 6; j++) paths[j][k] = pathRow[j];
  }
  const surfaceEmission = STEFAN_BOLTZMANN * skin ** 4;
  points.forEach((row, gp) => {
    surfaceUp[gp] = planckShare(row, skin) * surfaceEmission;
    for (let k = 0; k < K; k++) {
      const tau = LONGWAVE_CONSTANTS.diffusivity * (row[0] * paths[0][k] + row[1] * paths[1][k] + row[2] * paths[2][k] + row[3] * paths[3][k] + row[4] * paths[4][k] + row[5] * paths[5][k]);
      gas[gp * K + k] = -Math.expm1(-tau);
      planck[gp * K + k] = planckShare(row, T[k]) * STEFAN_BOLTZMANN * T[k] ** 4;
    }
  });
  const sums = [[0, 0], [0, 0], [0, 0]];
  let chainRandom = 0, chainRandomDown = 0;
  for (let gp = 0; gp < NG; gp++) {
    const r = randomFlux(gp), e = expected(gp, 'exponential'), m = expected(gp, 'maximum'), c = expected(gp, 'random');
    sums[0][0] += r[0]; sums[0][1] += r[1]; sums[1][0] += e[0]; sums[1][1] += e[1]; sums[2][0] += m[0]; sums[2][1] += m[1];
    chainRandom += c[0]; chainRandomDown += c[1];
  }
  sums.forEach(([up, down], n) => { out.olr[n][i] = up; out.down[n][i] = down; });
  worst = Math.max(worst, Math.abs(sums[0][0] - budget.outgoingLongwave) / budget.outgoingLongwave);
  worstDown = Math.max(worstDown, Math.abs(sums[0][1] - budget.downwardLongwave) / budget.downwardLongwave);
  worstChain = Math.max(worstChain, Math.abs(chainRandom - sums[0][0]) / sums[0][0], Math.abs(chainRandomDown - sums[0][1]) / sums[0][1]);
}
const sea = (i) => !landMask[i];
const regions = [
  ['global', () => true], ['30S-30N', (i) => Math.abs(lat[i]) <= 30],
  ['Pacific ITCZ 5-12N 160E-100W', (i) => lat[i] >= BOXES.itcz[0] && lat[i] <= BOXES.itcz[1] && inLongitudes(lon[i], BOXES.itcz[2], BOXES.itcz[3])],
  ['warm pool 10S-10N 120-170E sea', (i) => sea(i) && Math.abs(lat[i]) <= 10 && inLongitudes(lon[i], 120, 170)],
  ['SPCZ 20-5S 160E-150W sea', (i) => sea(i) && lat[i] >= -20 && lat[i] <= -5 && inLongitudes(lon[i], 160, -150)],
  ['N Pacific trades 15-25N 170-130W sea', (i) => sea(i) && lat[i] >= 15 && lat[i] <= 25 && inLongitudes(lon[i], -170, -130)],
  ['60-90N', (i) => lat[i] >= 60], ['60-90S', (i) => lat[i] <= -60],
];
const f = (x, d) => (Number.isFinite(x) ? x.toFixed(d) : 'n/a');
say(`longwave overlap of ${FILE.split('/').pop()}: day ${saved.day}, N=${saved.N}, one CPU step, then every column at t = ${f(model.time / 86400, 3)} d; ${C} columns, ${decks} with a deck`);
say(`the random recomputation repeats the radiation's OLR to ${worst.toExponential(1)} and its surface downward longwave to ${worstDown.toExponential(1)} relative in every column; the chain with alpha 0 repeats the random recomputation to ${worstChain.toExponential(1)}`);
for (const [name, keep] of regions) {
  let A = 0, model0 = 0, clear = 0;
  const olr = [0, 0, 0], down = [0, 0, 0];
  for (let i = 0; i < C; i++) {
    if (!keep(i)) continue;
    const a = area[i];
    A += a; model0 += a * out.model[i]; clear += a * out.clear[i];
    for (let n = 0; n < 3; n++) { olr[n] += a * out.olr[n][i]; down[n] += a * out.down[n][i]; }
  }
  say(`${name}: clear-sky OLR ${f(clear / A, 2)}; ${NAMES.map((label, n) => `${label} OLR ${f(olr[n] / A, 2)}, LWCRE ${f((clear - olr[n]) / A, 2)}, surface downward ${f(down[n] / A, 2)}`).join('; ')} W/m2; exponential-random less random: OLR ${f((olr[1] - olr[0]) / A, 2)}, LWCRE ${f((olr[0] - olr[1]) / A, 2)}, surface downward ${f((down[1] - down[0]) / A, 2)} W/m2`);
}
say(`(${f((performance.now() - t0) / 1000, 0)} s)`);
