// The headline numbers of the vertical-motion and convection audit (M21 in
// docs/c-grid-dynamical-core.md) from one saved state, on the single-thread
// CPU engine with the ocean off, each beside its Earth reference and a
// verdict:
//   node scripts/verticalAudit.mjs <state.bin>
// RADIATION (JSON) passes options to the model's radiation, e.g.
// '{"subsidenceSmoothing":0}' for a deck that reads the flux unsmoothed.
// Signs: omega (Pa/s) > 0 and sink (mm/s) > 0 are descent.
//
// From the state's own winds (stage 0 of the next step): omega at the
// model level nearest 700 and 500 hPa, the saved running-mean deck subsidence,
// the share of omega700's variance at the neighbouring-cell scale (the
// mean square of cell minus neighbour mean over the variance: white noise
// 1.167, a wave six cells long 0.065), the Hadley cells' peaks of the
// meridional mass streamfunction in its exact discrete form, and the sink
// at the deck's height h as the deck interpolates it, from pi*sigma-dot as
// the dynamics leaves it and as the deck reads it after its ring means,
// with its per-cell spread (the standard deviation over the box) and
// grid-scale share. Over a
// window of STEPS (8) steps: the rain and its convective share, the
// fraction of columns firing (convective rain above 1 mm/d in a step),
// the low cloud (a layer below 680 hPa with more than 1e-5 kg/kg of cloud
// water, or the deck's cover), the estimated inversion strength, and the
// mixed-layer deck's height and virtual potential temperature jump above
// it as the deck sees them at the start of each step's physics, and its
// gate replicated from the same inputs: each column-step runs the deck,
// or is off because the running-mean subsidence fails its test ("off:
// subsidence"), because the subsidence passes and the jump is below the
// minimum ("off: jump"), or because both pass but the gate's running mean
// has not yet turned ("off: gate memory"); beside the attribution, the
// share failing each test on its own, and the deck's height after the
// steps on which it ran. After the window, the resolved inversion (the
// interface of largest dθv/dz between 100 m and 3 km) and its θv jump.
// "u" is twice the grid-noise standard error of a
// snapshot's box mean; a verdict is "too noisy to tell" when u exceeds
// half the value.
//
// Earth references (September): deck boxes 0.1-0.3 mm/d with little
// convective rain and low cloud 0.6-0.7 (Wood 2012), omega700 0.03-0.05
// Pa/s and a boundary-layer-top divergence of 3-5e-6 /s (ERA-Interim),
// an inversion of 6-12 K at 1.0-1.5 km along 20S (VOCALS-REx) and an
// estimated inversion strength of 5-8 K (Wood and Bretherton 2006); the
// central and east Pacific ITCZ at 5-12N rising at 0.05-0.10 Pa/s at 500
// hPa with 6-9 mm/d, the zonal-mean rain peaking at 6-7 mm/d near 8N
// (GPCP); the winter Hadley cell 100-200e9 kg/s and the summer one
// 10-50e9 kg/s.
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { createModel, STATE_NAMES } from '../js/model.module.js';
import { decodeState, savedLevels } from '../js/stateFile.module.js';
import { savedDeckField, DECK_FIELDS } from '../js/physics/regrid.module.js';
import { divergence } from '../js/dynamics/operators.module.js';
import { createMixedLayer, dycomsLongwave } from '../js/physics/mixedLayer.module.js';
import { DECK_CLOUD_LEVELS, ringMean } from '../js/physics/radiation.module.js';
import { LATENT_HEAT } from '../js/physics/moist.module.js';
import { BOXES, inLongitudes } from '../js/audit.module.js';

const FILE = process.argv[2];
if (!FILE) { console.error('usage: node scripts/verticalAudit.mjs <state.bin>'); process.exit(1); }
const STEPS = Number(process.env.STEPS ?? 8), RADIATION = JSON.parse(process.env.RADIATION ?? '{}');
const t0 = performance.now();
const say = (s = '') => console.log(s);

const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const model = createModel(new Grid(saved.N), { topography, levels: savedLevels(saved), ocean: false, radiation: RADIATION });
const { mesh, core, state, phases, radiation, boundaryLayer: bl, moist, seaIce, land } = model;
const { K, C, E, levels, sigmaMid, sigmaLower, dSigma, R, g, cp, p0, exnerLayer, exnerLower, geopotential, piSigmaDot } = core.diagnostics;
const { thetaV } = core.arrays;
STATE_NAMES.forEach((name, a) => state[a].set(saved[name]));
seaIce.load(state[6], saved.concentration ?? null);
for (const field of Object.keys(DECK_FIELDS)) radiation[field].set(savedDeckField(saved, field, model));
land.load({ soil: Float64Array.from(saved.land.soil), snow: Float64Array.from(saved.land.snow), ...(saved.land.vegetation ? { vegetation: Float64Array.from(saved.land.vegetation) } : {}), ...(saved.land.surface ? { surface: Float64Array.from(saved.land.surface) } : {}) }, state[6]);
model.time = saved.time;
const [pi, theta, u, , q, qc, ice] = state;
const deg = 180 / Math.PI, area = mesh.areaCell, landMask = model.geography.land, dt = 1350 * 16 / saved.N, bottom = (K - 1) * C;
const lat = Float64Array.from(mesh.latCell, (x) => x * deg), lon = Float64Array.from(mesh.lonCell, (x) => x * deg);
const { maxEdges, nEdgesOnCell, cellsOnCell } = mesh;
const zs = Float64Array.from({ length: C }, (_, i) => (model.surfaceGeopotential ? model.surfaceGeopotential[i] / g : 0));
const gates = radiation.deckGates;

const sea = Uint8Array.from({ length: C }, (_, i) => (!landMask[i] && !(ice[i] > 0) ? 1 : 0));
const everywhere = new Uint8Array(C).fill(1);
const boxMask = ([south, north, west, east], base) => Uint8Array.from({ length: C }, (_, i) => (base[i] && lat[i] >= south && lat[i] <= north && inLongitudes(lon[i], west, east) ? 1 : 0));
const sePacific = boxMask(BOXES.sePacific, sea), peru = boxMask(BOXES.peru, sea), itcz = boxMask(BOXES.itcz, everywhere);

function mean(x, mask) {
  let s = 0, a = 0;
  for (let i = 0; i < C; i++) if (mask[i] && Number.isFinite(x[i])) { s += area[i] * x[i]; a += area[i]; }
  return a > 0 ? s / a : NaN;
}
function spread(x, mask) {
  const m = mean(x, mask);
  let s = 0, a = 0;
  for (let i = 0; i < C; i++) if (mask[i] && Number.isFinite(x[i])) { s += area[i] * (x[i] - m) ** 2; a += area[i]; }
  return Math.sqrt(s / a);
}
const neighbourMean = (x, i) => { let s = 0; for (let m = 0; m < nEdgesOnCell[i]; m++) s += x[cellsOnCell[maxEdges * i + m]]; return s / nEdgesOnCell[i]; };
function noise(x, mask) {
  let A = 0, m = 0; const idx = [];
  for (let i = 0; i < C; i++) {
    if (!mask[i] || !Number.isFinite(x[i])) continue;
    let ok = true; for (let n = 0; n < nEdgesOnCell[i]; n++) if (!Number.isFinite(x[cellsOnCell[maxEdges * i + n]])) ok = false;
    if (!ok) continue;
    idx.push(i); A += area[i]; m += area[i] * x[i];
  }
  m /= A;
  let v = 0, r = 0;
  for (const i of idx) { v += area[i] * (x[i] - m) ** 2; r += area[i] * (x[i] - neighbourMean(x, i)) ** 2; }
  return { ratio: r / v, u: 2 * Math.sqrt(r / A) / Math.sqrt(idx.length) };
}

// ---------------- stage 0: omega, the saved deck sink, the streamfunction ----------------
const dPi = new Float64Array(C);
core.phaseFlux(state, 0, K);
core.phaseColumn(state, [dPi], 0, C);
const D = Float64Array.from(new Float64Array(core.shared.divFlux)), psd = Float64Array.from(piSigmaDot), thv = Float64Array.from(thetaV), exM = Float64Array.from(exnerLayer);
const delta = new Float64Array(K * C), omega = new Float64Array(K * C);
for (let k = 0; k < K; k++) {
  divergence(mesh, u.subarray(k * E, (k + 1) * E), delta.subarray(k * C, (k + 1) * C));
  for (let i = 0; i < C; i++) { const idx = k * C + i; omega[idx] = 0.5 * (psd[idx] + psd[idx + C]) + sigmaMid[k] * (dPi[i] + D[idx] - pi[i] * delta[idx]); }
}
const shadow = createMixedLayer({ cp, R, g, latentHeat: LATENT_HEAT, referencePressure: p0, cloudLevels: DECK_CLOUD_LEVELS, ...(RADIATION.mixedLayer ?? {}) });
function deckGeometry(i, mixedDepth) {
  const b = bottom + i, surface = geopotential[b] - cp * thetaV[b] * (exnerLower[b] - exnerLayer[b]);
  const depth = mixedDepth + (geopotential[b] - surface) / g;
  let ceiling = Infinity;
  for (let k = K - 2; k >= 1; k--) {
    const upper = (geopotential[k * C + i] - surface) / g;
    if ((geopotential[(k + 1) * C + i] - surface) / g >= shadow.maximumHeight) break;
    if (upper > depth && thetaV[k * C + i] - thetaV[(k + 1) * C + i] >= gates.minimumInversion) { ceiling = upper - 1; break; }
  }
  const h = radiation.mlmHeight[i] > 0 ? shadow.bound(radiation.mlmHeight[i], depth, ceiling) : depth;
  let k = K - 1;
  for (; k >= 0 && geopotential[k * C + i] - surface < g * h; k--);
  if (k < 1) return null;
  const interfaceHeight = (m) => (geopotential[m * C + i] + cp * thetaV[m * C + i] * (exnerLayer[m * C + i] - exnerLower[(m - 1) * C + i]) - surface) / g;
  let lowerHeight = 0, lower = K, m = K - 1;
  for (; m > k && interfaceHeight(m) < h; m--) { lowerHeight = interfaceHeight(m); lower = m; }
  const upperHeight = interfaceHeight(m), density = pi[i] * sigmaMid[m] / (R * thetaV[m * C + i] * exnerLayer[m * C + i]);
  const sink = (passes) => {
    const lowerFlow = lower < K ? ringMean(mesh, piSigmaDot, lower * C, i, passes) : 0;
    return (lowerFlow + (ringMean(mesh, piSigmaDot, m * C, i, passes) - lowerFlow) * (h - lowerHeight) / (upperHeight - lowerHeight)) / (density * g);
  };
  return { surface, depth, ceiling, h, k, sink };
}
const deckSink = [0, gates.subsidenceSmoothing].map(() => new Float64Array(C).fill(NaN));
{
  const depthBefore = Float64Array.from(bl.depth);
  bl.diagnose(state, 0, C);
  for (let i = 0; i < C; i++) {
    const mixedDepth = bl.depth[i] - geopotential[bottom + i] / g;
    if (landMask[i] || !(1 - seaIce.cover(i, ice[i]) > 0) || !(mixedDepth > 0)) continue;
    const column = deckGeometry(i, mixedDepth);
    if (column) [0, gates.subsidenceSmoothing].forEach((passes, n) => { deckSink[n][i] = 1000 * column.sink(passes); });
  }
  bl.depth.set(depthBefore);
}
function omegaAt(P) {
  const om = new Float64Array(C).fill(NaN);
  for (let i = 0; i < C; i++) {
    if (100 * P > pi[i]) continue;
    let best = 0; for (let k = 1; k < K; k++) if (Math.abs(pi[i] * sigmaMid[k] - 100 * P) < Math.abs(pi[i] * sigmaMid[best] - 100 * P)) best = k;
    om[i] = omega[best * C + i];
  }
  return om;
}
const omega700 = omegaAt(700), omega500 = omegaAt(500);
const sinkSaved = Float64Array.from(radiation.mlmSubsidence, (v, i) => (landMask[i] ? NaN : -1000 * v));
const savedHeight = (mask) => { let a = 0, h = 0; for (let i = 0; i < C; i++) if (mask[i] && radiation.mlmHeight[i] > 0) { a += area[i]; h += area[i] * radiation.mlmHeight[i]; } return h / a; };

const NB = 180, binOf = (i) => Math.min(NB - 1, Math.max(0, Math.floor(lat[i] + 90)));
const capSum = new Float64Array(NB * K);
for (let i = 0; i < C; i++) { let cum = 0; const b = binOf(i); for (let k = 0; k < K; k++) { cum += D[k * C + i] * dSigma[k]; capSum[b * K + k] += area[i] * cum; } }
const psi = [];
for (let b = 0; b < NB; b++) { const row = new Float64Array(K); for (let bb = b; bb < NB; bb++) for (let k = 0; k < K; k++) row[k] -= capSum[bb * K + k] / g; psi.push(row); }
function hadley(south, north, sign) {
  let best = 0, at = NaN;
  for (let b = 0; b < NB; b++) {
    const L = b - 90; if (L < south || L > north) continue;
    for (let k = 0; k < K; k++) if (sigmaLower[k] >= 0.1 && sigmaLower[k] <= 0.95 && sign * psi[b][k] > sign * best) { best = psi[b][k]; at = L; }
  }
  return { value: best, lat: at };
}
const stage0 = (performance.now() - t0) / 1000;

// ---------------- the window: rain, firing, low cloud, the deck's start ----------------
const longwave = dycomsLongwave();
const acc = Object.fromEntries(['fire', 'low', 'cloudy', 'eis', 'eisN', 'attempt', 'h', 'jump', 'runs', 'offSubsidence', 'offJump', 'offMemory', 'failSubsidence', 'failJump', 'ran', 'ranH'].map((name) => [name, new Float64Array(C)]));
const predicted = { mean: new Float64Array(C).fill(NaN), gate: new Float64Array(C).fill(NaN), runs: new Uint8Array(C) };
const replica = { mean: 0, gate: 0, decisions: 0, checked: 0 };
function deckStart() {
  predicted.mean.fill(NaN);
  for (let i = 0; i < C; i++) {
    if (landMask[i] || !(sePacific[i] || peru[i])) continue;
    const b = bottom + i, mixedDepth = bl.depth[i] - geopotential[b] / g;
    if (!(1 - seaIce.cover(i, ice[i]) > 0) || !(mixedDepth > 0)) continue;
    const column = deckGeometry(i, mixedDepth);
    if (!column) continue;
    const { surface, h } = column;
    let weight = 0, heat = 0, water = 0, k = K - 1;
    for (; k >= 0 && geopotential[k * C + i] - surface < g * h; k--) {
      const idx = k * C + i, cloud = Math.max(0, qc[idx]);
      heat += dSigma[k] * (theta[idx] - LATENT_HEAT * cloud / (cp * exnerLayer[idx]));
      water += dSigma[k] * (Math.max(0, q[idx]) + cloud);
      weight += dSigma[k];
    }
    if (k < 1) continue;
    const above = k * C + i, aboveCloud = Math.max(0, qc[above]);
    const forcing = { surfacePressure: pi[i], sensibleHeat: 0, evaporation: 0, radiation: longwave, subsidence: () => 0,
      thetaLAbove: theta[above] - LATENT_HEAT * aboveCloud / (cp * exnerLayer[above]), qtAbove: Math.max(0, q[above]) + aboveCloud };
    const d = shadow.diagnose({ h, thetaL: heat / weight, qt: water / weight }, forcing);
    acc.attempt[i] += 1; acc.h[i] += h; acc.jump[i] += d.virtualJump;
    const keep = Math.exp(-dt / gates.subsidenceMemory), mean = radiation.mlmSubsidence[i] * keep - column.sink(gates.subsidenceSmoothing) * (1 - keep);
    const sinking = !(mean > -gates.stratusSubsidence), pass = sinking && d.virtualJump >= gates.minimumInversion ? 1 : 0;
    const gate = gates.gateMemory > 0 ? radiation.mlmGate[i] - (pass - radiation.mlmGate[i]) * Math.expm1(-dt / gates.gateMemory) : pass;
    const runs = gate > 0.5 || (gate === 0.5 && pass === 1);
    predicted.mean[i] = mean; predicted.gate[i] = gate; predicted.runs[i] = runs ? 1 : 0;
    if (runs) acc.runs[i] += 1;
    else if (pass) acc.offMemory[i] += 1;
    else if (!sinking) acc.offSubsidence[i] += 1;
    else acc.offJump[i] += 1;
    if (!sinking) acc.failSubsidence[i] += 1;
    if (!(d.virtualJump >= gates.minimumInversion)) acc.failJump[i] += 1;
  }
}
function deckAfter() {
  for (let i = 0; i < C; i++) {
    if (!Number.isFinite(predicted.mean[i])) continue;
    replica.checked++;
    replica.mean = Math.max(replica.mean, Math.abs(predicted.mean[i] - radiation.mlmSubsidence[i]));
    replica.gate = Math.max(replica.gate, Math.abs(predicted.gate[i] - radiation.mlmGate[i]));
    const ran = radiation.mlmTop[i] > 0;
    if (ran !== !!predicted.runs[i]) replica.decisions++;
    if (ran) { acc.ran[i] += 1; acc.ranH[i] += radiation.mlmHeight[i]; }
  }
}
const physicsPhase = phases.physics;
phases.physics = (...args) => {
  deckStart();
  physicsPhase(...args);
  deckAfter();
  for (let i = 0; i < C; i++) if (Number.isFinite(radiation.stabilityIndex[i])) { acc.eis[i] += radiation.stabilityIndex[i]; acc.eisN[i] += 1; }
};
model.restartPrecipitation();
const before = new Float64Array(C);
for (let n = 0; n < STEPS; n++) {
  before.set(moist.convectivePrecipitation);
  model.step(dt);
  for (let i = 0; i < C; i++) {
    if ((moist.convectivePrecipitation[i] - before[i]) * 86400 / dt > 1) acc.fire[i] += 1;
    let low = false, any = false;
    for (let k = 0; k < K; k++) if (qc[k * C + i] > 1e-5) { any = true; if (pi[i] * sigmaMid[k] > 680e2) low = true; }
    if (any) acc.cloudy[i] += 1;
    acc.low[i] += low ? 1 : Math.min(1, radiation.stratusFraction[i]);
  }
}
const perDay = 86400 / (STEPS * dt);
const convective = Float64Array.from(moist.convectivePrecipitation, (x) => x * perDay), rain = Float64Array.from(convective, (x, i) => x + moist.largeScalePrecipitation[i] * perDay);
const ratio = (sum, count, mask) => { let s = 0, w = 0; for (let i = 0; i < C; i++) if (mask[i] && count[i] > 0) { s += area[i] * sum[i]; w += area[i] * count[i]; } return s / w; };
for (let i = 0; i < C; i++) core.diagnoseColumn(i, pi, theta, q, qc);
const inversionZ = new Float64Array(C).fill(NaN), inversionJump = new Float64Array(C).fill(NaN);
for (let i = 0; i < C; i++) {
  if (landMask[i]) continue;
  let best = -Infinity;
  for (let k = K - 2; k >= 1; k--) {
    const zu = geopotential[k * C + i] / g - zs[i], zl = geopotential[(k + 1) * C + i] / g - zs[i];
    if (zu > 3000) break;
    if (zl < 100) continue;
    const jump = thetaV[k * C + i] - thetaV[(k + 1) * C + i], gradient = jump / (zu - zl);
    if (gradient > best) { best = gradient; inversionZ[i] = 0.5 * (zu + zl); inversionJump[i] = jump; }
  }
}
const zonalPeak = (() => {
  const s = new Float64Array(90), a = new Float64Array(90);
  for (let i = 0; i < C; i++) if (Math.abs(lat[i]) < 45) { const b = Math.min(89, Math.floor(lat[i] + 45)); s[b] += area[i] * rain[i]; a[b] += area[i]; }
  let best = 0; for (let b = 0; b < 90; b++) if (s[b] / a[b] > s[best] / a[best]) best = b;
  return { value: s[best] / a[best], lat: best - 45 + 0.5 };
})();

// ---------------- the table ----------------
const f = (x, d) => (Number.isFinite(x) ? x.toFixed(d) : 'n/a');
function verdict(value, lo, hi, u = 0) {
  if (![value, lo, hi].every(Number.isFinite)) return 'n/a';
  if (u > 0 && u >= Math.abs(value) / 2) return 'too noisy to tell';
  if (lo !== 0 && value !== 0 && Math.sign(value) !== Math.sign(lo)) return 'wrong sign';
  const m = Math.abs(value), low = Math.min(Math.abs(lo), Math.abs(hi)), high = Math.max(Math.abs(lo), Math.abs(hi));
  if (m >= low && m <= high) return 'matches';
  if (m === 0) return 'none';
  const factor = m < low ? low / m : m / high;
  return `${m < low ? 'too weak' : 'too strong'} by x${f(factor, factor < 2 ? 2 : 1)}`;
}
const rows = [];
const row = (name, value, digits, lo, hi, u = 0, note = '') => rows.push([name, value, digits, lo, hi, u, note]);
for (const [name, mask] of [['SE Pacific 10-30S 110-80W', sePacific], ['Peru 5-20S 90-75W', peru]]) {
  const total = mean(rain, mask), conv = mean(convective, mask), h = savedHeight(mask), n700 = noise(omega700, mask), nSink = noise(sinkSaved, mask);
  row(`${name}: rain (mm/d)`, total, 2, 0.1, 0.3);
  row(`${name}: convective share of the rain`, conv / total, 2, 0, 0.1);
  row(`${name}: columns firing a step`, mean(Float64Array.from(acc.fire, (x) => x / STEPS), mask), 3, 0, 0.01);
  row(`${name}: deck's virtual jump above h (K)`, ratio(acc.jump, acc.attempt, mask), 2, 6, 12);
  row(`${name}: estimated inversion strength (K)`, ratio(acc.eis, acc.eisN, mask), 2, 5, 8);
  row(`${name}: saved running-mean deck sink at h=${f(h, 0)} m (mm/s)`, mean(sinkSaved, mask), 2, 3e-3 * h, 5e-3 * h, nSink.u, `spread ${f(spread(sinkSaved, mask), 2)}, grid-scale share ${f(nSink.ratio, 3)}`);
  for (const [n, label] of [[0, 'as the dynamics leaves it'], [1, `as the deck reads it, ${gates.subsidenceSmoothing} ring passes`]]) {
    const x = deckSink[n], noiseX = noise(x, mask);
    row(`${name}: deck-height sink now, ${label} (mm/s)`, mean(x, mask), 2, 3e-3 * h, 5e-3 * h, noiseX.u, `spread ${f(spread(x, mask), 2)}, grid-scale share ${f(noiseX.ratio, 3)}`);
  }
  row(`${name}: omega700 (Pa/s)`, mean(omega700, mask), 4, 0.03, 0.05, n700.u);
  row(`${name}: low cloud`, mean(Float64Array.from(acc.low, (x) => x / STEPS), mask), 3, 0.6, 0.7, 0, `resolved cloud in any layer ${f(mean(Float64Array.from(acc.cloudy, (x) => x / STEPS), mask), 3)}`);
  row(`${name}: deck's start height h (m)`, ratio(acc.h, acc.attempt, mask), 0, 1000, 1500);
  const share = (x) => ratio(x, acc.attempt, mask);
  row(`${name}: deck runs, share of column-steps`, share(acc.runs), 3, 0.6, 1, 0, `off: subsidence ${f(share(acc.offSubsidence), 3)}, off: jump ${f(share(acc.offJump), 3)}, off: gate memory ${f(share(acc.offMemory), 3)}; failing the subsidence test ${f(share(acc.failSubsidence), 3)}, the jump test ${f(share(acc.failJump), 3)}`);
  row(`${name}: deck height where it runs (m)`, ratio(acc.ranH, acc.ran, mask), 0, 1000, 1500);
  row(`${name}: resolved inversion (m)`, mean(inversionZ, mask), 0, 1000, 1500);
  row(`${name}: resolved inversion's thetaV jump (K)`, mean(inversionJump, mask), 2, 6, 12);
}
row('Pacific ITCZ 5-12N 160E-100W: rain (mm/d)', mean(rain, itcz), 2, 6, 9, 0, `convective share ${f(mean(convective, itcz) / mean(rain, itcz), 2)}`);
row('Pacific ITCZ 5-12N 160E-100W: omega500 (Pa/s)', mean(omega500, itcz), 4, -0.05, -0.10, noise(omega500, itcz).u);
row('zonal-mean rain peak (mm/d)', zonalPeak.value, 2, 6, 7, 0, `at ${f(zonalPeak.lat, 1)} deg; Earth near 8N`);
const grid = noise(omega700, everywhere);
row('omega700 grid-scale share, global (white noise 1.167)', grid.ratio, 3, 0, 0.065, 0, `SE Pacific ${f(noise(omega700, sePacific).ratio, 3)}`);
const sh = hadley(-35, 15, -1), nh = hadley(0, 35, 1);
row('SH Hadley peak (1e9 kg/s)', sh.value / 1e9, 1, -100, -200, 0, `at ${sh.lat} deg`);
row('NH Hadley peak (1e9 kg/s)', nh.value / 1e9, 1, 10, 50, 0, `at ${nh.lat} deg`);

say(`vertical audit of ${FILE.split('/').pop()}: day ${saved.day}, N=${saved.N}, K=${K}; CPU, ocean off; ${STEPS} steps of ${dt} s after the stage-0 snapshot`);
say(`deck gates: ${gates.subsidenceSmoothing} ring passes, subsidence memory ${f(gates.subsidenceMemory / 86400, 2)} d, sink at least ${f(1000 * gates.stratusSubsidence, 2)} mm/s, jump at least ${f(gates.minimumInversion, 1)} K, gate memory ${f(gates.gateMemory / 86400, 2)} d; replica over ${replica.checked} column-steps: running mean |diff| ${replica.mean.toExponential(1)} m/s, gate |diff| ${replica.gate.toExponential(1)}, run decisions differing ${replica.decisions}`);
const width = Math.max(...rows.map(([name]) => name.length));
for (const [name, value, digits, lo, hi, u, note] of rows) {
  const range = `${f(lo, digits)}..${f(hi, digits)}`;
  say(`  ${name.padEnd(width)}  ${f(value, digits).padStart(9)}${u ? ` u ${f(u, digits)}` : ''}  Earth ${range}  -> ${verdict(value, lo, hi, u)}${note ? `  [${note}]` : ''}`);
}
say(`(${f(stage0, 0)} s to the snapshot, ${f((performance.now() - t0) / 1000, 0)} s in all)`);
