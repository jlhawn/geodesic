// The headline numbers of the vertical-motion and convection audit (M21 in
// docs/c-grid-dynamical-core.md) from one saved state, on the single-thread
// CPU engine with the ocean off, each beside its Earth reference and a
// verdict:
//   node scripts/verticalAudit.mjs <state.bin>
// RADIATION (JSON) passes options to the model's radiation, e.g.
// '{"subsidenceSmoothing":0}' for a deck that reads the flux unsmoothed,
// MOIST (JSON) to its moist physics, BOUNDARY_LAYER (JSON) to its
// boundary layer and SURFACE (JSON) to its surface.
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
// grid-scale share. SKIP (0) steps are run first and the snapshot taken
// after them, so windows can start at different hours. Over a
// window of STEPS (8) steps: the rain and its convective share, the
// fraction of columns firing (convective rain above 1 mm/d in a step)
// and of those whose temperature convection changes at all, the global
// rain and the convective share of 15S-15N, the convective heating
// profile of the Pacific ITCZ's firing columns (the pressure of its
// maximum, and its mean over the lowest 100 m) and the box's large-scale
// heating (condensation, autoconversion and rain evaporation) below 1 km,
// the low cloud (a layer below 680 hPa with more than 1e-5 kg/kg of cloud
// water, or the deck's cover), the estimated inversion strength, and the
// mixed-layer deck's height and virtual potential temperature jump above
// it as the deck sees them at the start of each step's physics, and its
// gate replicated from the same inputs: each column-step runs the deck,
// or is off because the running-mean subsidence fails its test ("off:
// subsidence"), because the subsidence passes and the jump is below the
// minimum ("off: jump"), because both pass but the gate's running mean
// has not yet turned ("off: gate memory"), or because its regime stands it
// down ("stood down"); beside the attribution, the
// share failing each test on its own, and the deck's height after the
// steps on which it ran, with its cloud layer's thickness at the step's
// start and its liquid water path, uncapped and as the radiation caps it,
// with the share of its running steps on which the cap binds.
// After the window, the resolved inversion (the
// interface of largest dθv/dz between 100 m and 3 km) and its θv jump.
// The cloud-radiative effects, global and over 30S-30N: shortwave (ASR
// less clear-sky ASR) and longwave (clear-sky OLR less OLR), the day means
// of the state's last day where the state carries them, else the
// window's means (the radiation runs its clear-sky pass here unless
// RADIATION sets clearSkyPass false), with the window's in the note.
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
// 10-50e9 kg/s; the global cloud-radiative effects -47 +- 4 W/m2
// shortwave and +26 +- 3 W/m2 longwave (CERES EBAF).
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
import { REGIME } from '../js/physics/boundaryLayer.module.js';

const FILE = process.argv[2];
if (!FILE) { console.error('usage: node scripts/verticalAudit.mjs <state.bin>'); process.exit(1); }
const STEPS = Number(process.env.STEPS ?? 8), SKIP = Number(process.env.SKIP ?? 0), RADIATION = JSON.parse(process.env.RADIATION ?? '{}'), MOIST = JSON.parse(process.env.MOIST ?? '{}'), BOUNDARY_LAYER = JSON.parse(process.env.BOUNDARY_LAYER ?? '{}'), SURFACE = JSON.parse(process.env.SURFACE ?? '{}');
const t0 = performance.now();
const say = (s = '') => console.log(s);

const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const model = createModel(new Grid(saved.N), { topography, levels: savedLevels(saved), ocean: false, radiation: { clearSkyPass: true, ...RADIATION }, moist: MOIST, boundaryLayer: BOUNDARY_LAYER, surface: SURFACE });
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
// the GPU engine that saves states measures heights from the surface, this engine from sea level
if (saved.boundaryDepth) bl.depth.set(Float64Array.from(saved.boundaryDepth, (z, i) => z + zs[i]));
if (saved.mixingTop) bl.mixingTop.set(Float64Array.from(saved.mixingTop, (z, i) => z + zs[i]));
if (saved.boundaryRegime) bl.regime.set(saved.boundaryRegime);
if (saved.boundaryBuoyancy) bl.buoyancyFlux.set(saved.boundaryBuoyancy);
for (let n = 0; n < SKIP; n++) model.step(dt);
if (SKIP > 0) radiation.restartSums();
const gates = radiation.deckGates;

const sea = Uint8Array.from({ length: C }, (_, i) => (!landMask[i] && !(ice[i] > 0) ? 1 : 0));
const everywhere = new Uint8Array(C).fill(1);
const boxMask = ([south, north, west, east], base) => Uint8Array.from({ length: C }, (_, i) => (base[i] && lat[i] >= south && lat[i] <= north && inLongitudes(lon[i], west, east) ? 1 : 0));
const sePacific = boxMask(BOXES.sePacific, sea), peru = boxMask(BOXES.peru, sea), itcz = boxMask(BOXES.itcz, everywhere), namibia = boxMask(BOXES.namibia, sea), california = boxMask(BOXES.california, sea);

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
  const interfaceHeight = (m) => (geopotential[m * C + i] + cp * thetaV[m * C + i] * (exnerLayer[m * C + i] - exnerLower[(m - 1) * C + i]) - surface) / g;
  let ceiling = Infinity, capping = -1;
  for (let k = K - 2; k >= 1; k--) {
    if ((geopotential[(k + 1) * C + i] - surface) / g >= shadow.maximumHeight) break;
    if ((geopotential[k * C + i] - surface) / g > depth && thetaV[k * C + i] - thetaV[(k + 1) * C + i] >= gates.ceilingInversion) { capping = k; ceiling = (geopotential[k * C + i] - surface) / g - 1; break; }
  }
  const regime = gates.deckRest === 'regime' ? bl.regime[i] : -1;
  const standDown = (regime === REGIME.SURFACE || regime === REGIME.DECOUPLED) && !(capping >= 0 && interfaceHeight(capping + 1) <= gates.cumulusCeiling);
  const resting = gates.deckRest !== 'depth' && regime !== REGIME.COUPLED && !standDown && ceiling < shadow.maximumHeight ? ceiling : depth;
  const h = radiation.mlmHeight[i] > 0 ? shadow.bound(radiation.mlmHeight[i], depth, ceiling) : resting;
  let k = K - 1;
  for (; k >= 0 && geopotential[k * C + i] - surface < g * h; k--);
  if (k < 1) return null;
  let lowerHeight = 0, lower = K, m = K - 1;
  for (; m > k && interfaceHeight(m) < h; m--) { lowerHeight = interfaceHeight(m); lower = m; }
  const upperHeight = interfaceHeight(m), density = pi[i] * sigmaMid[m] / (R * thetaV[m * C + i] * exnerLayer[m * C + i]);
  const sink = (passes) => {
    const lowerFlow = lower < K ? ringMean(mesh, piSigmaDot, lower * C, i, passes) : 0;
    return (lowerFlow + (ringMean(mesh, piSigmaDot, m * C, i, passes) - lowerFlow) * (h - lowerHeight) / (upperHeight - lowerHeight)) / (density * g);
  };
  return { surface, depth, ceiling, h, k, sink, standDown };
}
const deckSink = [0, gates.subsidenceSmoothing].map(() => new Float64Array(C).fill(NaN));
{
  const before = ['depth', 'regime', 'mixingTop', 'buoyancyFlux'].map((name) => Float64Array.from(bl[name]));
  bl.diagnose(state, 0, C);
  for (let i = 0; i < C; i++) {
    const mixedDepth = bl.depth[i] - geopotential[bottom + i] / g;
    if (landMask[i] || !(1 - seaIce.cover(i, ice[i]) > 0) || !(mixedDepth > 0)) continue;
    const column = deckGeometry(i, mixedDepth);
    if (column) [0, gates.subsidenceSmoothing].forEach((passes, n) => { deckSink[n][i] = 1000 * column.sink(passes); });
  }
  ['depth', 'regime', 'mixingTop', 'buoyancyFlux'].forEach((name, n) => bl[name].set(before[n]));
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
const acc = Object.fromEntries(['fire', 'convect', 'low', 'cloudy', 'eis', 'eisN', 'attempt', 'h', 'jump', 'runs', 'offSubsidence', 'offJump', 'offMemory', 'offStood', 'failSubsidence', 'failJump', 'ran', 'ranH', 'lowCover', 'lowWater', 'ranThick', 'ranCloudy', 'ranPath', 'ranPathSeen', 'ranCapped', 'evaporation', 'stable', 'surface', 'decoupled', 'coupled', 'cooling', 'velocity', 'entrainment', 'mixingTop'].map((name) => [name, new Float64Array(C)]));
const predicted = { mean: new Float64Array(C).fill(NaN), gate: new Float64Array(C).fill(NaN), runs: new Uint8Array(C) };
const replica = { mean: 0, gate: 0, decisions: 0, checked: 0 };
function deckStart() {
  predicted.mean.fill(NaN);
  for (let i = 0; i < C; i++) {
    if (landMask[i] || !(sePacific[i] || peru[i] || namibia[i] || california[i])) continue;
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
    const regimeTest = gates.deckRegime === 'boundaryLayer' ? bl.regime[i] === REGIME.COUPLED : d.virtualJump >= gates.minimumInversion;
    const sinking = !(mean > -gates.stratusSubsidence), pass = sinking && regimeTest && !column.standDown ? 1 : 0;
    const gate = column.standDown ? 0 : gates.gateMemory > 0 ? radiation.mlmGate[i] - (pass - radiation.mlmGate[i]) * Math.expm1(-dt / gates.gateMemory) : pass;
    const runs = (gate > 0.5 || (gate === 0.5 && pass === 1)) && !gates.deckBypass;
    predicted.mean[i] = mean; predicted.gate[i] = gate; predicted.runs[i] = runs ? 1 : 0;
    if (runs) {
      acc.runs[i] += 1;
      if (d.cloudy) { acc.ranCloudy[i] += 1; acc.ranThick[i] += h - d.cloudBase; }
    } else if (column.standDown) acc.offStood[i] += 1;
    else if (pass) acc.offMemory[i] += 1;
    else if (!sinking) acc.offSubsidence[i] += 1;
    else acc.offJump[i] += 1;
    if (!sinking) acc.failSubsidence[i] += 1;
    if (!regimeTest) acc.failJump[i] += 1;
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
    if (ran) {
      acc.ran[i] += 1; acc.ranH[i] += radiation.mlmHeight[i];
      acc.ranPath[i] += radiation.mlmWater[i]; acc.ranPathSeen[i] += Math.min(gates.stratusWaterMax, radiation.mlmWater[i]);
      if (radiation.mlmWater[i] > gates.stratusWaterMax) acc.ranCapped[i] += 1;
    }
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
moist.trace.convection = new Float64Array(K * C);
moist.trace.largeScale = new Float64Array(K * C);
const heating = { convection: new Float64Array(K), largeScale: new Float64Array(K), pressure: new Float64Array(K), height: new Float64Array(K), fired: 0, area: 0 };
const tropicsMask = boxMask([-15, 15, -180, 180], everywhere), TOP_BINS = 10;
const plumes = { tops: new Float64Array(TOP_BINS), deep: 0, shallow: 0, area: 0, deepFlux: 0, shallowFlux: 0 };
for (let n = 0; n < STEPS; n++) {
  before.set(moist.convectivePrecipitation);
  moist.trace.convection.fill(0);
  moist.trace.largeScale.fill(0);
  model.step(dt);
  for (let i = 0; i < C; i++) {
    const fired = (moist.convectivePrecipitation[i] - before[i]) * 86400 / dt > 1;
    if (fired) acc.fire[i] += 1;
    let convected = false;
    for (let k = 0; k < K; k++) if (moist.trace.convection[k * C + i] !== 0) convected = true;
    if (convected) acc.convect[i] += 1;
    if (itcz[i]) {
      const a = area[i];
      heating.area += a;
      if (fired) heating.fired += a;
      for (let k = 0; k < K; k++) {
        const idx = k * C + i;
        heating.largeScale[k] += a * moist.trace.largeScale[idx];
        heating.pressure[k] += a * pi[i] * sigmaMid[k];
        heating.height[k] += a * (geopotential[idx] / g - zs[i]);
        if (fired) heating.convection[k] += a * moist.trace.convection[idx];
      }
    }
    if (tropicsMask[i]) {
      plumes.area += area[i];
      if (moist.cumulusBaseFlux[i] > 0) {
        const top = moist.cumulusTop[i];
        plumes.tops[Math.min(TOP_BINS - 1, Math.floor(top / 100e2))] += area[i];
        if (top < pi[i] * moist.deepSigma) { plumes.deep += area[i]; plumes.deepFlux += area[i] * moist.cumulusBaseFlux[i]; } else { plumes.shallow += area[i]; plumes.shallowFlux += area[i] * moist.cumulusBaseFlux[i]; }
      }
    }
    let low = false, any = false;
    for (let k = 0; k < K; k++) if (qc[k * C + i] > 1e-5) { any = true; if (pi[i] * sigmaMid[k] > 680e2) low = true; }
    if (any) acc.cloudy[i] += 1;
    acc.low[i] += low ? 1 : Math.min(1, radiation.stratusFraction[i]);
    acc.lowCover[i] += radiation.lowCover[i]; acc.lowWater[i] += radiation.lowWater[i];
    acc.evaporation[i] += 86400 * radiation.evaporation[i];
    acc[['stable', 'surface', 'decoupled', 'coupled'][bl.regime[i]]][i] += 1;
    acc.cooling[i] += bl.cloudTopCooling[i]; acc.velocity[i] += bl.radiativeVelocity[i]; acc.entrainment[i] += bl.entrainment[i]; acc.mixingTop[i] += bl.mixingTop[i] - geopotential[bottom + i] / g;
  }
}
const perDay = 86400 / (STEPS * dt);
if (radiation.clearSkyPass) radiation.readMeans(STEPS);
const itczProfile = (() => {
  const perStep = 86400 / dt, conv = Float64Array.from(heating.convection, (x) => x / heating.fired * perStep), large = Float64Array.from(heating.largeScale, (x) => x / heating.area * perStep);
  const p = Float64Array.from(heating.pressure, (x) => x / heating.area / 100), z = Float64Array.from(heating.height, (x) => x / heating.area);
  let peak = 0, low = 0, lowMass = 0, most = 0, least = 0;
  for (let k = 0; k < K; k++) {
    if (conv[k] > conv[peak]) peak = k;
    if (z[k] < 100) { low += conv[k] * dSigma[k]; lowMass += dSigma[k]; }
    if (z[k] < 1000) { if (large[k] > large[most]) most = k; if (large[k] < large[least]) least = k; }
  }
  return { peakP: heating.fired > 0 ? p[peak] : NaN, peak: conv[peak], low: low / lowMass, fired: heating.fired / heating.area, most: large[most], mostZ: z[most], least: large[least], leastZ: z[least] };
})();
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
function deckRows(name, mask) {
  const share = (x) => ratio(x, acc.attempt, mask);
  row(`${name}: deck's start height h (m)`, share(acc.h), 0, 1000, 1500);
  row(`${name}: deck runs, share of column-steps`, share(acc.runs), 3, 0.6, 1, 0, `off: subsidence ${f(share(acc.offSubsidence), 3)}, off: jump ${f(share(acc.offJump), 3)}, off: gate memory ${f(share(acc.offMemory), 3)}, stood down ${f(share(acc.offStood), 3)}; failing the subsidence test ${f(share(acc.failSubsidence), 3)}, the jump test ${f(share(acc.failJump), 3)}`);
  row(`${name}: deck height where it runs (m)`, ratio(acc.ranH, acc.ran, mask), 0, 1000, 1500);
  row(`${name}: deck's cloud-layer thickness where it runs (m)`, ratio(acc.ranThick, acc.ranCloudy, mask), 0, 200, 400, 0, `cloudy on ${f(ratio(acc.ranCloudy, acc.runs, mask), 3)} of its running steps, at the step's start`);
  row(`${name}: deck's liquid water path where it runs (g/m2)`, 1000 * ratio(acc.ranPath, acc.ran, mask), 1, 50, 150, 0, `as the radiation takes it, capped at ${f(1000 * gates.stratusWaterMax, 0)}: ${f(1000 * ratio(acc.ranPathSeen, acc.ran, mask), 1)}, the cap binding on ${f(ratio(acc.ranCapped, acc.ran, mask), 3)} of its running steps`);
}
const windowMean = (x, mask) => mean(Float64Array.from(x, (v) => v / STEPS), mask);
function boundaryRows(name, mask) {
  const cover = windowMean(acc.lowCover, mask);
  row(`${name}: low-cloud cover, radiative`, cover, 3, 0.6, 0.7, 0, `presence ${f(windowMean(acc.low, mask), 3)}`);
  row(`${name}: low-cloud water path in cloud (g/m2)`, 1000 * windowMean(acc.lowWater, mask) / cover, 1, 50, 150, 0, `grid mean ${f(1000 * windowMean(acc.lowWater, mask), 1)}`);
  row(`${name}: boundary-layer regime, coupled stratocumulus share`, windowMean(acc.coupled, mask), 3, NaN, NaN, 0, `decoupled ${f(windowMean(acc.decoupled, mask), 3)}, surface-driven ${f(windowMean(acc.surface, mask), 3)}, stable ${f(windowMean(acc.stable, mask), 3)}`);
  row(`${name}: cloud-top cooling (W/m2)`, windowMean(acc.cooling, mask), 1, NaN, NaN, 0, `V ${f(windowMean(acc.velocity, mask), 2)} m/s, w_e ${f(1000 * windowMean(acc.entrainment, mask), 2)} mm/s, mixing top ${f(windowMean(acc.mixingTop, mask), 0)} m`);
}
for (const [name, mask] of [['SE Pacific 10-30S 110-80W', sePacific], ['Peru 5-20S 90-75W', peru]]) {
  const total = mean(rain, mask), conv = mean(convective, mask), h = savedHeight(mask), n700 = noise(omega700, mask), nSink = noise(sinkSaved, mask);
  row(`${name}: rain (mm/d)`, total, 2, 0.1, 0.3);
  row(`${name}: convective share of the rain`, conv / total, 2, 0, 0.1);
  row(`${name}: columns firing a step`, mean(Float64Array.from(acc.fire, (x) => x / STEPS), mask), 3, 0, 0.01, 0, `convecting ${f(mean(Float64Array.from(acc.convect, (x) => x / STEPS), mask), 3)}`);
  row(`${name}: deck's virtual jump above h (K)`, ratio(acc.jump, acc.attempt, mask), 2, 6, 12);
  row(`${name}: estimated inversion strength (K)`, ratio(acc.eis, acc.eisN, mask), 2, 5, 8);
  row(`${name}: saved running-mean deck sink at h=${f(h, 0)} m (mm/s)`, mean(sinkSaved, mask), 2, 3e-3 * h, 5e-3 * h, nSink.u, `spread ${f(spread(sinkSaved, mask), 2)}, grid-scale share ${f(nSink.ratio, 3)}`);
  for (const [n, label] of [[0, 'as the dynamics leaves it'], [1, `as the deck reads it, ${gates.subsidenceSmoothing} ring passes`]]) {
    const x = deckSink[n], noiseX = noise(x, mask);
    row(`${name}: deck-height sink now, ${label} (mm/s)`, mean(x, mask), 2, 3e-3 * h, 5e-3 * h, noiseX.u, `spread ${f(spread(x, mask), 2)}, grid-scale share ${f(noiseX.ratio, 3)}`);
  }
  row(`${name}: omega700 (Pa/s)`, mean(omega700, mask), 4, 0.03, 0.05, n700.u);
  row(`${name}: low cloud`, mean(Float64Array.from(acc.low, (x) => x / STEPS), mask), 3, 0.6, 0.7, 0, `resolved cloud in any layer ${f(mean(Float64Array.from(acc.cloudy, (x) => x / STEPS), mask), 3)}`);
  deckRows(name, mask);
  row(`${name}: resolved inversion (m)`, mean(inversionZ, mask), 0, 1000, 1500);
  row(`${name}: resolved inversion's thetaV jump (K)`, mean(inversionJump, mask), 2, 6, 12);
  boundaryRows(name, mask);
}
for (const [name, mask] of [['Namibia 10-20S 0-10E', namibia], ['California 20-30N 130-120W', california]]) {
  row(`${name}: rain (mm/d)`, mean(rain, mask), 2, 0.1, 0.3);
  deckRows(name, mask);
  row(`${name}: resolved inversion (m)`, mean(inversionZ, mask), 0, 1000, 1500);
  row(`${name}: resolved inversion's thetaV jump (K)`, mean(inversionJump, mask), 2, 6, 12);
  boundaryRows(name, mask);
}
row('Pacific ITCZ 5-12N 160E-100W: rain (mm/d)', mean(rain, itcz), 2, 6, 9, 0, `convective share ${f(mean(convective, itcz) / mean(rain, itcz), 2)}`);
row('Pacific ITCZ 5-12N 160E-100W: omega500 (Pa/s)', mean(omega500, itcz), 4, -0.05, -0.10, noise(omega500, itcz).u);
row('Pacific ITCZ 5-12N 160E-100W: firing columns\' convective heating peak (hPa)', itczProfile.peakP, 0, 400, 500, 0, `${f(itczProfile.peak, 2)} K/d; ${f(itczProfile.low, 2)} K/d over the lowest 100 m; firing ${f(itczProfile.fired, 3)} of the column-steps`);
row('Pacific ITCZ 5-12N 160E-100W: large-scale heating below 1 km, largest |K/d|', Math.max(itczProfile.most, -itczProfile.least), 2, NaN, NaN, 0, `${f(itczProfile.most, 2)} at ${f(itczProfile.mostZ, 0)} m, ${f(itczProfile.least, 2)} at ${f(itczProfile.leastZ, 0)} m`);
const tropics = boxMask([-15, 15, -180, 180], everywhere);
row('global rain (mm/d)', mean(rain, everywhere), 2, 2.6, 2.8, 0, `convective share ${f(mean(convective, everywhere) / mean(rain, everywhere), 2)}`);
row('global evaporation (mm/d)', windowMean(acc.evaporation, everywhere), 2, 2.6, 2.8);
row('convective share of the rain, 15S-15N', mean(convective, tropics) / mean(rain, tropics), 2, NaN, NaN, 0, `15S-15N rain ${f(mean(rain, tropics), 2)} mm/d`);
row('zonal-mean rain peak (mm/d)', zonalPeak.value, 2, 6, 7, 0, `at ${f(zonalPeak.lat, 1)} deg; Earth near 8N`);
{
  const cloudBand = boxMask([-30, 30, -180, 180], everywhere), dayMeans = !!saved.meanShortwaveCloudEffect;
  const field = (name) => (dayMeans ? Float64Array.from(saved[name]) : radiation[name]);
  const windowNote = (name, mask) => (radiation.clearSkyPass ? `the window's ${STEPS} steps ${f(mean(radiation[name], mask), 1)}` : 'no clear-sky pass in the window');
  const source = dayMeans ? 'day mean of the state\'s last day' : `mean over the window's ${STEPS} steps`;
  if (dayMeans || radiation.clearSkyPass) for (const [name, mask, global] of [['global', everywhere, true], ['30S-30N', cloudBand, false]]) {
    row(`${name} shortwave cloud effect, ${source} (W/m2)`, mean(field('meanShortwaveCloudEffect'), mask), 1, global ? -43 : NaN, global ? -51 : NaN, 0, `${windowNote('meanShortwaveCloudEffect', mask)}${global ? '' : '; Earth -47 +- 4 globally'}`);
    row(`${name} longwave cloud effect, ${source} (W/m2)`, mean(field('meanLongwaveCloudEffect'), mask), 1, global ? 23 : NaN, global ? 29 : NaN, 0, `${windowNote('meanLongwaveCloudEffect', mask)}${global ? '' : '; Earth +26 +- 3 globally'}`);
  }
}
const grid = noise(omega700, everywhere);
row('omega700 grid-scale share, global (white noise 1.167)', grid.ratio, 3, 0, 0.065, 0, `SE Pacific ${f(noise(omega700, sePacific).ratio, 3)}`);
const sh = hadley(-35, 15, -1), nh = hadley(0, 35, 1);
row('SH Hadley peak (1e9 kg/s)', sh.value / 1e9, 1, -100, -200, 0, `at ${sh.lat} deg`);
row('NH Hadley peak (1e9 kg/s)', nh.value / 1e9, 1, 10, 50, 0, `at ${nh.lat} deg`);

say(`vertical audit of ${FILE.split('/').pop()}: day ${saved.day}, N=${saved.N}, K=${K}; CPU, ocean off; ${SKIP ? `${SKIP} steps, then ` : ''}${STEPS} steps of ${dt} s after the stage-0 snapshot`);
say(`deck gates: ${gates.subsidenceSmoothing} ring passes, subsidence memory ${f(gates.subsidenceMemory / 86400, 2)} d, sink at least ${f(1000 * gates.stratusSubsidence, 2)} mm/s, jump at least ${f(gates.minimumInversion, 1)} K, gate memory ${f(gates.gateMemory / 86400, 2)} d; replica over ${replica.checked} column-steps: running mean |diff| ${replica.mean.toExponential(1)} m/s, gate |diff| ${replica.gate.toExponential(1)}, run decisions differing ${replica.decisions}`);
const width = Math.max(...rows.map(([name]) => name.length));
for (const [name, value, digits, lo, hi, u, note] of rows) {
  const range = `${f(lo, digits)}..${f(hi, digits)}`;
  say(`  ${name.padEnd(width)}  ${f(value, digits).padStart(9)}${u ? ` u ${f(u, digits)}` : ''}  Earth ${range}  -> ${verdict(value, lo, hi, u)}${note ? `  [${note}]` : ''}`);
}
say(`plumes over 15S-15N (share of the column-steps): deep ${f(plumes.deep / plumes.area, 3)} (mean base flux ${f(plumes.deepFlux / plumes.deep, 4)} kg/m2/s), shallow ${f(plumes.shallow / plumes.area, 3)} (${f(plumes.shallowFlux / plumes.shallow, 4)}); tops by 100 hPa of their top interface: ${Array.from(plumes.tops, (x, b) => `${b * 100}-${b * 100 + 100} ${f(x / plumes.area, 3)}`).join(', ')}`);
{
  const bins = [];
  for (let west = -180; west < 180; west += 20) {
    const mask = boxMask([-5, 5, west, west + 20], everywhere);
    bins.push(`${west < 0 ? `${-west}W` : `${west}E`} ${f(mean(omega500, mask), 3)}`);
  }
  say(`equatorial (5S-5N) omega500 by 20 degrees of longitude from the west edge (Pa/s, + descent): ${bins.join(', ')}`);
}
say(`(${f(stage0, 0)} s to the snapshot, ${f((performance.now() - t0) / 1000, 0)} s in all)`);
