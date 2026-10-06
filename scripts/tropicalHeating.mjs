// The heating and the rain of the tropical rain bands by process and layer,
// on the single-thread CPU engine with the ocean off:
//   node scripts/tropicalHeating.mjs <state.bin>
// STEPS (one day) steps from the state; RADIATION, MOIST, BOUNDARY_LAYER and
// SURFACE (JSON) pass options as to scripts/verticalAudit.mjs. Boxes
// (TROPICAL_BOXES of js/audit.module.js): the Pacific ITCZ (5-12N 160E-100W,
// all cells, the audit's), the warm pool (10S-10N 120-170E sea), the SPCZ
// (20-5S 160E-150W sea), the north Pacific trades (15-25N 170-130W sea) and
// the Amazon (10S-2N 70-50W land).
//
// Every phase of every step is measured in each box column as it runs: the
// dynamics (the RK4 step, the closure, the kinetic energy returned as heat),
// the radiation (its longwave from the radiation's own per-layer record, the
// lowest layer's sensible heat from the bulk formula with the radiation's
// inputs, the shortwave the rest), the boundary layer's mixing (its change of
// T_l = T - L q_c/c_p; the cloud it evaporates is counted with the
// condensation), the moist step, the fusion heat of snow and the dry
// adjustment. The moist step is replayed column by column on copies with
// instances of the moist physics built as the model builds its own, and the
// replay must reproduce the model's columns bit for bit: the uniform
// condensation, the deep plume (with plumeClosure 'cape', which runs the same
// deep plume without the shallow one: the latent heat of the rain it makes,
// the evaporation into its downdraft, and the rest, the two drafts' transport
// of s_l), the shallow plume (the difference), the condensation after the
// plumes, the rain's fall and evaporation (a copy of autoconvertColumn that
// also tags the rain by origin: liquid converted below or above 700 hPa, and
// falling ice), and the adjustment after the ice fall. The sum of the terms
// closes on each column's temperature change. The total water q_t = q + q_c
// is charged by process at the same points (the surface evaporation, the
// boundary layer's mixing, the condensation, the deep and shallow plumes,
// the evaporation of large-scale and deep-plume rain, the conversion to rain,
// the ice fall with its melting, the filler and the dry adjustment), and its
// terms close on each column's q_t change.
//
// Printed: the checks; per box the surface rain by source with the
// stratiform share (melted falling ice and large-scale conversion above 700
// hPa over all the rain), the profiles in K/day and W/m2 per layer, the
// moistening by process in K/day (L/cp dq_t/dt), the apparent heat source Q1
// (all physics), Q1 without radiation (Q1R) and the moisture sink
// Q2 = -(L/cp) dq_t/dt of the physics, each with its peak by layer and over
// 50 hPa bins and its centroid (heatingProfile of js/audit.module.js), and
// the column integrals of Q1R and Q2 by process in W/m2; over land boxes, the
// share of the convective rain by local solar hour; the firing columns'
// convective heating peak over the day (the layer of largest mean heating in
// the moist physics' convective trace, without mass weights, over the ITCZ
// column-steps raining more than 1 mm/d of convective rain)
// with its terms; the deep plume's fate in each box (deck veto, no condensation
// level, plume topping below 700 hPa, CAPE at most plumeCape, closure or
// inhibition, fired), its CAPE, inhibition, base flux and top, and the level
// at which it stops being buoyant against the same plume without entrainment;
// the plume of the same column under saturation adjustment (condensation
// 'saturation', from the same total water and liquid-water temperature) and
// the condensate the uniform condensation forms and converts below 700 hPa in
// the columns where only that plume fires; and the box's temperature and
// relative humidity against Jordan's (1958) mean West Indies sounding.
import { readFileSync } from 'node:fs';
import { Grid } from '../js/grid.module.js';
import { topographyFromInt16 } from '../js/geography.module.js';
import { createModel, STATE_NAMES } from '../js/model.module.js';
import { decodeState, savedLevels } from '../js/stateFile.module.js';
import { savedDeckField, DECK_FIELDS } from '../js/physics/regrid.module.js';
import { createMoistPhysics, MOIST_DEFAULTS, SOURCE_EXCESS, FUSION_HEAT, IFS_PRECIPITATION, surfaceLayerVelocity, LATENT_HEAT, R_VAPOR, CLEAR_AIR, DECK_OPEN, DECK_CLOSED, COUPLED_REGIME, DEEP_REFERENCE, DEEP_CLOUD_DEPTH, IFS_ENTRAINMENT, saturationHumidity, saturationVaporPressure, cloudSaturation, criticalHumidityAt, uniformCover, fogConstants, liquidFraction, liftingCondensationLevel } from '../js/physics/moist.module.js';
import { VIRTUAL_FACTOR } from '../js/dynamics/sigmaCore.module.js';
import { SEA_DRAG, LAND_DRAG } from '../js/physics/surface.module.js';
import { FREEZING_POINT } from '../js/physics/ice.module.js';
import { TROPICAL_BOXES, tropicalBoxOf, heatingProfile, bulkSensible, lowestHeight, layerExner } from '../js/audit.module.js';

const FILE = process.argv[2];
if (!FILE) { console.error('usage: node scripts/tropicalHeating.mjs <state.bin>'); process.exit(1); }
const RADIATION = JSON.parse(process.env.RADIATION ?? '{}'), MOIST = JSON.parse(process.env.MOIST ?? '{}'), BOUNDARY_LAYER = JSON.parse(process.env.BOUNDARY_LAYER ?? '{}'), SURFACE = JSON.parse(process.env.SURFACE ?? '{}');
const t0 = performance.now();
const say = (s = '') => console.log(s);

const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const model = createModel(new Grid(saved.N), { topography, levels: savedLevels(saved), ocean: false, radiation: RADIATION, moist: MOIST, boundaryLayer: BOUNDARY_LAYER, surface: SURFACE });
const { mesh, core, state, phases, radiation, boundaryLayer: bl, moist, seaIce, land, surface } = model;
const { K, C, levels, sigmaMid, sigmaLower, sigmaUpper, dSigma, R, g, cp, p0, kappa, exnerLayer, exnerLower, geopotential } = core.diagnostics;
const { thetaV } = core.arrays;
STATE_NAMES.forEach((name, a) => state[a].set(saved[name]));
seaIce.load(state[6], saved.concentration ?? null);
for (const field of Object.keys(DECK_FIELDS)) radiation[field].set(savedDeckField(saved, field, model));
land.load(saved.land, state[6]);
model.time = saved.time;
const [pi, theta, wind, surfaceT, q, qc, ice] = state;
const deg = 180 / Math.PI, area = mesh.areaCell, landMask = model.geography.land, dt = 1350 * 16 / saved.N, bottom = K - 1;
const STEPS = Number(process.env.STEPS ?? Math.round(86400 / dt));
const zs = Float64Array.from({ length: C }, (_, i) => (model.surfaceGeopotential ? model.surfaceGeopotential[i] / g : 0));
if (saved.boundaryDepth) bl.depth.set(Float64Array.from(saved.boundaryDepth, (z, i) => z + zs[i]));
if (saved.mixingTop) bl.mixingTop.set(Float64Array.from(saved.mixingTop, (z, i) => z + zs[i]));
if (saved.boundaryRegime) bl.regime.set(saved.boundaryRegime);
if (saved.boundaryBuoyancy) bl.buoyancyFlux.set(saved.boundaryBuoyancy);
if (saved.subcloudVirtual && saved.subcloudVirtual.length === moist.subcloudVirtual.length) moist.subcloudVirtual.set(saved.subcloudVirtual);

const O = { ...MOIST_DEFAULTS, ...MOIST }, bechtold = O.capeClosure === 'bechtold', mixedPlume = O.plumePhase === 'mixed';
if (O.plumeClosure === 'maximum') throw new Error('the deep and shallow split needs plumeClosure separate or cape');
const L = O.latentHeat;
const lon = Float64Array.from(mesh.lonCell, (x) => x * deg);
const BOX_LIST = TROPICAL_BOXES;
const boxOf = tropicalBoxOf(mesh, landMask);
const traced = [];
for (let i = 0; i < C; i++) if (boxOf[i] >= 0) traced.push(i);
const NB = BOX_LIST.length;

const twinOptions = (extra) => ({
  boundaryDepth: bl.depth, boundaryRegime: bl.regime, deckGate: radiation.mlmGate,
  boundaryTop: bl.turbulence === 'moist' ? bl.mixingTop : null, boundaryCloudLayer: bl.turbulence === 'moist' ? bl.cloudLayer : null, stratiform: radiation.stratiform,
  surfaceBuoyancy: bl.buoyancyFlux, frictionVelocity: bl.friction, land: model.geography.land, iceSheet: model.geography.iceSheet ?? null, surfaceDrag: bl.implicitDrag ? bl.surfaceDrag : null, surfaceSensible: radiation.sensibleHeat, surfaceEvaporation: radiation.evaporation, buffers: { subcloudVirtual: moist.shared.subcloudVirtual },
  ...Object.fromEntries(['liquidTemperature', 'iceTemperature'].filter((key) => key in RADIATION).map((key) => [key, RADIATION[key]])),
  ...MOIST, ...extra,
});
const replica = createMoistPhysics(mesh, core, twinOptions({}));
const deepOnly = createMoistPhysics(mesh, core, twinOptions({ plumeClosure: 'cape' }));
const adjusted = createMoistPhysics(mesh, core, twinOptions({ condensation: 'saturation', plumeClosure: 'cape' }));
for (const twin of [replica, deepOnly, adjusted]) twin.useSeaIce(seaIce.concentration);
const boundaryTop = bl.turbulence === 'moist' ? bl.mixingTop : null;
const gustiness = RADIATION.gustiness ?? 3;
const FOG = fogConstants(O.fogDroplets), foggy = O.fogDroplets !== null, depositing = O.fogDeposition > 0 && bl.implicitDrag;
const continental = Uint8Array.from({ length: C }, (_, i) => (landMask[i] && !(model.geography.iceSheet && model.geography.iceSheet[i]) ? 1 : 0));

const TERMS = ['dynamics', 'closure', 'dissipation', 'shortwave', 'longwave', 'sensible', 'boundaryLayer', 'condensationMixed', 'condensationUniform', 'deepRain', 'deepFreezing', 'deepMelting', 'deepDowndraftEvaporation', 'deepTransport', 'shallow', 'recondensation', 'largeScaleEvaporation', 'convectiveEvaporation', 'afterFallPositive', 'afterFallNegative', 'snow', 'dryAdjustment'];
const SHORT = ['dyn', 'clos', 'diss', 'SW', 'LW', 'sens', 'BL', 'condBL', 'condU', 'dRain', 'dFrz', 'dMelt', 'dDDev', 'dTrans', 'shal', 'recond', 'LSev', 'CVev', 'fall+', 'fall-', 'snow', 'dry'];
const DYNAMICS = new Set(['dynamics', 'closure']);
const T_ = Object.fromEntries(TERMS.map((t, n) => [t, n]));
const NT = TERMS.length;
const WATER = ['dynamics', 'closure', 'evaporation', 'boundaryLayer', 'condensation', 'deep', 'shallow', 'largeScaleEvaporation', 'convectiveEvaporation', 'conversion', 'iceFall', 'filler', 'dryAdjustment'];
const WATER_SHORT = ['dyn', 'clos', 'evap', 'BL', 'cond', 'deep', 'shal', 'LSev', 'CVev', 'conv', 'ice', 'fill', 'dry'];
const W_ = Object.fromEntries(WATER.map((t, n) => [t, n]));
const NW = WATER.length, HOURS = 24;
const acc = Array.from({ length: NB }, () => ({
  heat: new Float64Array(NT * K), energy: new Float64Array(NT * K), total: new Float64Array(K), water: new Float64Array(NW * K), waterMass: new Float64Array(NW * K), totalWater: new Float64Array(K), pressure: new Float64Array(K), thickness: new Float64Array(K), height: new Float64Array(K), areaSteps: 0,
  temperature: new Float64Array(K), humidity: new Float64Array(K), localHour: new Float64Array(HOURS),
  rain: { deep: 0, shallow: 0, liquidLow: 0, liquidHigh: 0, iceMelted: 0, iceGround: 0, model: 0, modelConvective: 0 },
  detrained: new Float64Array(K), detrainedPressure: new Float64Array(K), deepMade: 0,
  plume: { typed: 0, deck: 0, noLcl: 0, notCloudy: 0, shallowTop: new Float64Array(10), weakCape: 0, weakCapeSum: 0, closed: 0, closedCape: 0, closedInhibition: 0, fired: 0, firedCape: 0, firedInhibition: 0, firedFlux: 0, firedTop: new Float64Array(10), limited: 0, candidateCape: new Float64Array(8), consumption: 0, cloudBase: 0, cloudBaseN: 0, pcape: 0, pcapeBoundary: 0, tau: 0, speed: 0, taus: [], closedByForcing: 0, closedByForcingSum: 0,
    dilute: new Float64Array(10), undilute: new Float64Array(10), diluteN: 0, neverBuoyant: 0, pureNeverBuoyant: 0, diluteCape: 0, pureCape: 0 },
  counter: { subcloud: 0, cloudLayer: 0, fires: 0, firesOnly: 0, both: 0, actualOnly: 0, flux: 0, fluxActual: 0, low: new Float64Array(4), lowOnly: new Float64Array(4), temperature: new Float64Array(K), humidity: new Float64Array(K), pressure: new Float64Array(K), weight: 0 },
}));
const audit = { conv: new Float64Array(K), pressure: new Float64Array(K), fired: 0, area: 0, parts: new Float64Array(5 * K) };
const checks = { replicaColumns: 0, replicaMismatch: 0, shadowMismatch: 0, probeChecked: 0, probeMismatch: 0, deepFluxMismatch: 0, traceOff: 0, closure: 0, waterClosure: 0, waterUnassigned: 0, sensibleGlobal: 0, sensibleModel: 0, precipOff: 0 };

const ex = new Float64Array(K);
const exnerOf = (i, out) => layerExner(pi[i], core.diagnostics, out);
const nT = traced.length;
const stepStart = new Float64Array(nT * K), stepStartQt = new Float64Array(nT * K), mark = new Float64Array(nT * K), markQ = new Float64Array(nT * K), markQc = new Float64Array(nT * K), exTrue = new Float64Array(nT * K);
const stepHeat = new Float64Array(nT * NT * K), stepWater = new Float64Array(nT * NW * K);
const slot = new Int32Array(C).fill(-1);
traced.forEach((i, n) => { slot[i] = n; });
const water = (n, term, k, dq) => { stepWater[(n * NW + W_[term]) * K + k] += dq; };
function charge(n, i, term, waterTerm) {
  for (let k = 0; k < K; k++) {
    const idx = k * C + i, s = n * K + k, dq = q[idx] + qc[idx] - (markQ[s] + markQc[s]);
    stepHeat[(n * NT + T_[term]) * K + k] += (theta[idx] - mark[s]) * exTrue[s];
    if (waterTerm) water(n, waterTerm, k, dq); else checks.waterUnassigned = Math.max(checks.waterUnassigned, Math.abs(dq));
    mark[s] = theta[idx]; markQ[s] = q[idx]; markQc[s] = qc[idx];
  }
}
const add = (n, term, k, dT) => { stepHeat[(n * NT + T_[term]) * K + k] += dT; };

const physicsPhase = phases.physics, closurePhase = phases.closure, dissipatePhase = phases.dissipate, oceanPhase = phases.ocean;
const mixColumn = bl.mixColumn, moistAdjust = moist.adjust, dryAdjust = surface.convectiveAdjustment;
phases.ocean = (...args) => {
  traced.forEach((i, n) => {
    exnerOf(i, ex);
    for (let k = 0; k < K; k++) {
      const s = n * K + k, idx = k * C + i;
      exTrue[s] = ex[k]; stepHeat[(n * NT + T_.dynamics) * K + k] += theta[idx] * ex[k] - stepStart[s];
      water(n, 'dynamics', k, q[idx] + qc[idx] - stepStartQt[s]);
      mark[s] = theta[idx]; markQ[s] = q[idx]; markQc[s] = qc[idx];
    }
  });
  oceanPhase(...args);
};
const staleExner = new Float64Array(nT * K), physicsAir = new Float64Array(C), physicsSurface = new Float64Array(C), physicsExner = new Float64Array(C);
const sensibleOptions = { seaDrag: SURFACE.dragCoefficient ?? SEA_DRAG, landDrag: LAND_DRAG, freezing: FREEZING_POINT, gustiness };
phases.physics = (iFrom, iTo, step, sums) => {
  for (let i = 0; i < C; i++) { physicsExner[i] = exnerLayer[bottom * C + i]; physicsAir[i] = theta[bottom * C + i] * physicsExner[i]; physicsSurface[i] = surfaceT[i]; }
  traced.forEach((i, n) => { for (let k = 0; k < K; k++) staleExner[n * K + k] = exnerLayer[k * C + i]; });
  const coverBefore = Float64Array.from({ length: C }, (_, i) => seaIce.cover(i, ice[i])), heightBefore = Float64Array.from({ length: C }, (_, i) => lowestHeight(model, i));
  physicsPhase(iFrom, iTo, step, sums);
  let sensibleSum = 0;
  const sensibleOf = (i) => bulkSensible(model, i, physicsAir[i], physicsSurface[i], coverBefore[i], heightBefore[i], sensibleOptions);
  for (let i = 0; i < C; i++) sensibleSum += area[i] * sensibleOf(i);
  checks.sensibleGlobal = sensibleSum; checks.sensibleModel = sums.sensibleHeat;
  traced.forEach((i, n) => {
    const ratio = (k) => exTrue[n * K + k] / staleExner[n * K + k];
    for (let k = 0; k < K; k++) {
      const s = n * K + k, mass = pi[i] * dSigma[k] / g, total = (theta[k * C + i] - mark[s]) * exTrue[s];
      const longwave = radiation.longwave[k * C + i] * step / (cp * mass) * ratio(k);
      const sensible = k === bottom ? sensibleOf(i) * step / (cp * mass) * ratio(k) : 0;
      add(n, 'longwave', k, longwave); add(n, 'sensible', k, sensible); add(n, 'shortwave', k, total - longwave - sensible);
      water(n, 'evaporation', k, q[k * C + i] + qc[k * C + i] - (markQ[s] + markQc[s]));
      mark[s] = theta[k * C + i]; markQ[s] = q[k * C + i]; markQc[s] = qc[k * C + i];
    }
  });
};
phases.closure = (...args) => { closurePhase(...args); traced.forEach((i, n) => charge(n, i, 'closure', 'closure')); };
phases.dissipate = (...args) => { dissipatePhase(...args); traced.forEach((i, n) => charge(n, i, 'dissipation', null)); };
const preMixQc = new Float64Array(nT * K), mixedPhase = new Float64Array(nT * K);
bl.mixColumn = (i, piArray, thetaArray, qArray, qcArray, step) => {
  const n = slot[i];
  if (n >= 0) for (let k = 0; k < K; k++) preMixQc[n * K + k] = qc[k * C + i];
  mixColumn(i, piArray, thetaArray, qArray, qcArray, step);
  if (n < 0) return;
  for (let k = 0; k < K; k++) {
    const s = n * K + k, idx = k * C + i, dT = (theta[idx] - mark[s]) * exTrue[s], phase = L / cp * (qc[idx] - markQc[s]);
    add(n, 'boundaryLayer', k, dT - phase); mixedPhase[n * K + k] += phase;
    water(n, 'boundaryLayer', k, q[idx] + qc[idx] - (markQ[s] + markQc[s]));
    mark[s] = theta[idx]; markQ[s] = q[idx]; markQc[s] = qc[idx];
  }
};
surface.convectiveAdjustment = (piArray, thetaArray, iFrom, iTo, qArray, qcArray) => {
  traced.forEach((i, n) => charge(n, i, 'snow', 'dryAdjustment'));
  dryAdjust(piArray, thetaArray, iFrom, iTo, qArray, qcArray);
  traced.forEach((i, n) => charge(n, i, 'dryAdjustment', 'dryAdjustment'));
};

const scratch = Array.from({ length: 4 }, () => ({ theta: new Float64Array(K * C), q: new Float64Array(K * C), qc: new Float64Array(K * C) }));
function copyColumn(i, from, to) { for (let k = 0; k < K; k++) { const idx = k * C + i; to.theta[idx] = from.theta[idx]; to.q[idx] = from.q[idx]; to.qc[idx] = from.qc[idx]; } }
const live = { theta, q, qc };
const columnT = new Float64Array(K), columnQ = new Float64Array(K), columnQc = new Float64Array(K);
function hold(i, s) { for (let k = 0; k < K; k++) { const idx = k * C + i; columnT[k] = s.theta[idx]; columnQ[k] = s.q[idx]; columnQc[k] = s.qc[idx]; } }
function heated(i, s, n, term) {
  for (let k = 0; k < K; k++) {
    const idx = k * C + i, dT = (s.theta[idx] - columnT[k]) * exnerLayer[idx];
    add(n, term, k, dT);
    columnT[k] = s.theta[idx]; columnQ[k] = s.q[idx]; columnQc[k] = s.qc[idx];
  }
}
const saturated = { qs: 0, slope: 0, liquid: 1 };
const saturation = (T, p) => cloudSaturation(T, p, O.iceSaturation, O.liquidTemperature, O.iceTemperature, saturated);
const halfWidth = (p, ps) => (1 - criticalHumidityAt(p, ps, O.surfaceCriticalHumidity, O.topCriticalHumidity, O.criticalExponent)) * saturated.qs / (1 + L * saturated.slope / cp);
const upperInterface = (i, k) => (geopotential[k * C + i] + cp * thetaV[k * C + i] * (exnerLayer[k * C + i] - exnerLower[(k - 1) * C + i])) / g;
function longCloudShare(twin, i, k, iced) {
  if (twin.cumulusBaseFlux[i] > 0 && pi[i] * sigmaMid[k] >= twin.cumulusTop[i]) return 0;
  let share = 0;
  if (boundaryTop !== null) {
    if (!(geopotential[k * C + i] / g < boundaryTop[i])) share = radiation.stratiform[i];
    else if (bl.regime[i] === COUPLED_REGIME) share = 1;
  }
  return Math.max(share, iced);
}
const fall = { evaporated: new Float64Array(K), convective: new Float64Array(K), sublimated: new Float64Array(K), refunded: new Float64Array(K), converted: new Float64Array(K), melted: new Float64Array(K), tags: { low: 0, high: 0, ice: 0 }, ground: 0, rain: 0, streamed: 0, convectiveLeft: 0 };
const reserve = new Float64Array(K);
function shadowFall(twin, i, s, stream, iced, frozenStream = null) {
  const th = s.theta, qq = s.q, cc = s.qc;
  fall.evaporated.fill(0); fall.convective.fill(0); fall.converted.fill(0); fall.melted.fill(0); fall.sublimated.fill(0); fall.refunded.fill(0);
  let frozen = 0;
  const tags = fall.tags; tags.low = 0; tags.high = 0; tags.ice = 0;
  let rain = 0, convective = 0, streamed = 0, descending = 0, settling = 0;
  const floor = O.autoconversionFloor === 'boundaryLayer' && bl.depth ? bl.depth[i] : null;
  for (let k = 0; k < K; k++) {
    const idx = k * C + i;
    if (rain > 0 && (O.evaporationInCloud || !(cc[idx] > CLEAR_AIR)) && O.rainEvaporation > 0) {
      const exl = exnerLayer[idx], mass = pi[i] * dSigma[k] / g;
      const temperature = th[idx] * exl;
      const qs = O.iceSaturation ? saturation(temperature, pi[i] * sigmaMid[k]).qs : saturationHumidity(temperature, pi[i] * sigmaMid[k]);
      const slope = O.iceSaturation ? saturated.slope : qs * L / (R_VAPOR * temperature * temperature);
      const deficit = Math.max(0, (qs - qq[idx]) / (1 + L * slope / cp)) * mass;
      const evaporated = Math.min(rain, O.rainEvaporation * deficit);
      if (evaporated > 0) {
        const before = rain;
        rain = Math.max(0, rain - evaporated);
        const keep = before > 0 ? rain / before : 0;
        tags.low *= keep; tags.high *= keep; tags.ice *= keep;
        fall.evaporated[k] += evaporated;
        qq[idx] += evaporated / mass;
        th[idx] -= L * evaporated / (mass * cp * exl);
      }
    }
    if (stream) {
      convective = Math.max(0, convective + stream[k]);
      if (frozenStream) {
        const arriving = frozen + frozenStream[k];
        frozen = Math.min(convective, Math.max(0, arriving));
        const fusion = FUSION_HEAT * (Math.max(0, arriving) - frozen + Math.min(0, arriving)) / (pi[i] * dSigma[k] / g * cp);
        if (fusion !== 0) { th[idx] -= fusion / exnerLayer[idx]; fall.refunded[k] -= fusion; }
      }
      const spare = convective - reserve[k];
      if (spare > 0 && k > twin.deep.base && O.plumeRainEvaporation > 0 && O.rainEvaporation > 0 && (O.evaporationInCloud || !(cc[idx] > CLEAR_AIR))) {
        const exl = exnerLayer[idx], mass = pi[i] * dSigma[k] / g;
        const temperature = th[idx] * exl;
        const qs = saturationHumidity(temperature, pi[i] * sigmaMid[k]);
        const slope = qs * L / (R_VAPOR * temperature * temperature);
        const airborne = -convective * Math.expm1(-O.plumeRainEvaporation * Math.max(0, 1 - qq[idx] / qs) * R * temperature * dSigma[k] / (sigmaMid[k] * g));
        const evaporated = Math.min(spare, airborne, O.rainEvaporation * Math.max(0, (qs - qq[idx]) / (1 + L * slope / cp)) * mass);
        if (evaporated > 0) {
          const sublimated = frozen > 0 ? evaporated * frozen / convective : 0;
          convective -= evaporated;
          streamed += evaporated;
          fall.convective[k] += evaporated;
          qq[idx] += evaporated / mass;
          if (frozenStream) {
            frozen = Math.max(0, frozen - sublimated);
            fall.sublimated[k] += sublimated;
            th[idx] -= (L * evaporated + FUSION_HEAT * sublimated) / (mass * cp * exl);
          } else th[idx] -= L * evaporated / (mass * cp * exl);
        }
      }
    }
    let liquid = cc[idx];
    if (O.iceFall !== null) {
      const mass = pi[i] * dSigma[k] / g;
      if (descending > 0) {
        const melted = descending * liquidFraction(th[idx] * exnerLayer[idx], O.liquidTemperature, O.iceTemperature);
        cc[idx] += (descending - melted) / mass;
        rain += melted;
        tags.ice += melted;
        fall.melted[k] += melted;
        descending = 0;
      }
      liquid = cc[idx];
      if (cc[idx] > 0) {
        const temperature = th[idx] * exnerLayer[idx], pressure = pi[i] * sigmaMid[k];
        const share = liquidFraction(temperature, O.liquidTemperature, O.iceTemperature), iceWater = (1 - share) * cc[idx];
        liquid = share * cc[idx];
        if (iceWater > 0) {
          saturation(temperature, pressure);
          const cover = Math.max(CLEAR_AIR, uniformCover(cc[idx], halfWidth(pressure, pi[i])));
          const speed = O.iceFall * Math.pow(pressure / (R * temperature) * iceWater / cover, O.iceFallExponent);
          const courant = speed * dt * sigmaMid[k] * g / (R * temperature * dSigma[k]), leaving = iceWater * courant / (1 + courant);
          cc[idx] -= leaving;
          descending = leaving * mass;
        }
      }
    }
    if (settling > 0) {
      const added = settling * g / (pi[i] * dSigma[k]);
      cc[idx] += added;
      liquid += added;
      settling = 0;
    }
    if (!(cc[idx] > 0)) continue;
    const floored = O.autoconversionFloor !== 'none' && (floor === null ? k >= K - 2 : k > 0 && upperInterface(i, k) < floor);
    if (floored) {
      if (liquid > 0 && (foggy || (depositing && k === K - 1))) {
        const temperature = th[idx] * exnerLayer[idx], pressure = pi[i] * sigmaMid[k], mass = pi[i] * dSigma[k] / g;
        let courant = 0;
        if (foggy) {
          const speed = (continental[i] ? FOG.settleLand : FOG.settleSea) * Math.pow(pressure / (R * temperature) * liquid, 2 / 3);
          courant = speed * dt * sigmaMid[k] * g / (R * temperature * dSigma[k]);
        }
        if (depositing && k === K - 1) courant += Math.min(O.fogDeposition * bl.surfaceDrag[i], O.fogDepositionLimit * pressure / (R * temperature)) * dt * g / (pi[i] * dSigma[k]);
        const leaving = liquid * courant / (1 + courant);
        cc[idx] -= leaving;
        if (k === K - 1) {
          rain += leaving * mass;
          if (pressure > 700e2) tags.low += leaving * mass; else tags.high += leaving * mass;
        } else settling = leaving * mass;
        if (foggy) {
          const kept = liquid - leaving;
          if (kept > 0) {
            const converted = -kept * Math.expm1(-(continental[i] ? FOG.drizzleLand : FOG.drizzleSea) * Math.pow(kept, 1.47) * dt);
            cc[idx] -= converted;
            const made = mass * converted;
            rain += made;
            fall.converted[k] += made;
            if (pressure > 700e2) tags.low += made; else tags.high += made;
          }
        }
      }
      continue;
    }
    const excess = Math.max(0, liquid - O.autoconversionThreshold);
    let lifetime = O.upperCloudLifetime !== null && pi[i] * sigmaMid[k] < O.shallowTop ? O.upperCloudLifetime : O.cloudLifetime;
    if (O.stratiformLifetime !== null) lifetime += longCloudShare(twin, i, k, iced) * (O.stratiformLifetime - lifetime);
    const converted = Math.min(liquid, excess * (1 - Math.exp(-O.autoconversionRate * dt)) + liquid * (1 - Math.exp(-dt / lifetime)));
    cc[idx] -= converted;
    const made = pi[i] * dSigma[k] / g * converted;
    rain += made;
    fall.converted[k] += made;
    if (pi[i] * sigmaMid[k] > 700e2) tags.low += made; else tags.high += made;
  }
  fall.ground = descending; fall.rain = rain; fall.streamed = streamed; fall.convectiveLeft = convective;
  return rain + descending;
}

// The deep plume's ascent as plumeColumn makes it, with where it stops being
// buoyant, and the same plume without entrainment.
const env = { T: new Float64Array(K), p: new Float64Array(K), dp: new Float64Array(K), z: new Float64Array(K), s: new Float64Array(K), q: new Float64Array(K) };
const plumeAir = { T: 0, liquid: 0, ice: 0 };
function saturatedTemperature(energy, pressure, guess) {
  let t = guess;
  for (let n = 0; n < 4; n++) {
    const qs = saturationHumidity(t, pressure);
    t -= (cp * t + L * qs - energy) / (cp + L * L * qs / (R_VAPOR * t * t));
  }
  return t;
}
const plumeSaturated = { qs: 0, slope: 0, liquid: 1 };
function plumeState(energy, water, height, pressure, guess) {
  const dry = (energy - g * height) / cp;
  if (mixedPlume) {
    const liquidT = RADIATION.liquidTemperature ?? O.liquidTemperature, iceT = RADIATION.iceTemperature ?? O.iceTemperature;
    plumeAir.ice = 0;
    if (!(water > cloudSaturation(dry, pressure, true, liquidT, iceT, plumeSaturated).qs)) { plumeAir.T = dry; plumeAir.liquid = 0; return; }
    const target = energy - g * height + L * water, span = 1 / (liquidT - iceT);
    let t = Math.max(dry, guess);
    for (let n = 0; n < 4; n++) {
      const { qs, slope, liquid } = cloudSaturation(t, pressure, true, liquidT, iceT, plumeSaturated), held = water - qs, frozen = 1 - liquid;
      t -= (cp * t + L * qs - FUSION_HEAT * frozen * held - target) / (cp + (L + FUSION_HEAT * frozen) * slope + (liquid > 0 && liquid < 1 ? FUSION_HEAT * span * held : 0));
    }
    cloudSaturation(t, pressure, true, liquidT, iceT, plumeSaturated);
    plumeAir.T = t; plumeAir.liquid = Math.max(0, water - plumeSaturated.qs); plumeAir.ice = (1 - plumeSaturated.liquid) * plumeAir.liquid;
    return;
  }
  if (!(water > saturationHumidity(dry, pressure))) { plumeAir.T = dry; plumeAir.liquid = 0; return; }
  const t = saturatedTemperature(energy - g * height + L * water, pressure, Math.max(dry, guess));
  plumeAir.T = t;
  plumeAir.liquid = Math.max(0, water - saturationHumidity(t, pressure));
}
const ramp = (x) => Math.min(1, Math.max(0, x));
function ascend(i, s, entraining) {
  const { T, p, dp, z } = env, envS = env.s, envQ = env.q;
  const bottomK = K - 1, depthBL = bl.depth[i];
  let mass = 0, energy = 0, water = 0, source = bottomK;
  const surface = O.plumeSourceDepth === 'surface50';
  for (let k = bottomK; k >= 0; k--) {
    if (k < bottomK && (((surface || !(upperInterface(i, k + 1) < depthBL)) && !(p[k] >= pi[i] - O.cumulusSourceDepth)) || !(p[k] > O.shallowTop))) break;
    mass += dp[k]; energy += dp[k] * envS[k]; water += dp[k] * envQ[k]; source = k;
  }
  let sourceS = O.plumeSource === 'lowest' ? envS[bottomK] : energy / mass, sourceQ = O.plumeSource === 'lowest' ? envQ[bottomK] : water / mass;
  if (surface) {
    const buoyancy = bl.buoyancyFlux[i], density = p[bottomK] / (R * T[bottomK]), layer = O.excessVelocity === 'surfaceLayer', b = bottomK * C + i;
    const velocity = layer ? surfaceLayerVelocity(radiation.sensibleHeat[i], radiation.evaporation[i], density, T[bottomK], cp * thetaV[b] * (exnerLower[b] - exnerLayer[b]) / g, cp, g)
      : Math.max(buoyancy > 0 ? Math.cbrt(buoyancy * Math.max(0, depthBL - z[bottomK])) : 0, bl.friction[i]);
    if (velocity > 0) {
      const dT = Math.min(SOURCE_EXCESS.temperature, SOURCE_EXCESS.coefficient * radiation.sensibleHeat[i] / (density * cp * velocity));
      const dq = Math.min(SOURCE_EXCESS.humidity, SOURCE_EXCESS.coefficient * radiation.evaporation[i] / (density * velocity));
      sourceS += cp * (layer ? Math.max(0, dT) : dT); sourceQ += layer ? Math.max(0, dq) : dq;
    }
  }
  const lcl = liftingCondensationLevel((sourceS - g * z[bottomK]) / cp, sourceQ, p[bottomK], kappa);
  if (!lcl || !(lcl.pressure > pi[i] * levels[1])) return { status: 'noLcl' };
  const virtual = O.virtualBuoyancy ? VIRTUAL_FACTOR : 0, loading = O.virtualBuoyancy ? 1 : 0;
  let sp = sourceS, wq = sourceQ, w2 = 0, below = 0, inhibition = 0, cloudy = false, started = false, top = -1, guess = 0, cape = 0, base = -1, unbuoyant = NaN, buoyant = false, baseSaturation = 0;
  const ifs = O.plumeEntrainmentLaw === 'ifs', liquidT = RADIATION.liquidTemperature ?? O.liquidTemperature, iceT = RADIATION.iceTemperature ?? O.iceTemperature, sat = { qs: 0, slope: 0, liquid: 1 };
  for (let k = source - 1; k >= 0; k--) {
    const lower = upperInterface(i, k + 1), upper = k > 0 ? upperInterface(i, k) : Infinity, depth = upper - lower;
    const mixes = pi[i] * levels[k + 1] <= lcl.pressure;
    if (mixes && !started) { started = true; base = k + 1; w2 = O.plumeVelocity * O.plumeVelocity; if (ifs) baseSaturation = cloudSaturation(T[k], p[k], O.iceSaturation, liquidT, iceT, sat).qs; }
    let epsilon = 0, mixing = 0;
    if (ifs) {
      const humidity = Math.min(1, Math.max(0, s.q[k * C + i]) / cloudSaturation(T[k], p[k], O.iceSaturation, liquidT, iceT, sat).qs);
      const detrained = mixes ? IFS_ENTRAINMENT.detrainment * (IFS_ENTRAINMENT.detrainmentHumidity - humidity) : 0;
      epsilon = entraining && mixes && below > 0 ? IFS_ENTRAINMENT.entrainment * (IFS_ENTRAINMENT.humidity - humidity) * (sat.qs / baseSaturation) ** IFS_ENTRAINMENT.scaleExponent : 0;
      mixing = entraining ? IFS_ENTRAINMENT.drag * (epsilon > 0 ? epsilon : detrained) : 0;
    } else {
      epsilon = !entraining ? 0 : mixes ? Math.max(O.plumeEntrainmentFloor, O.plumeEntrainment * Math.max(0, below) / w2) : 0;
      mixing = O.plumeDrag * epsilon;
    }
    const half = Math.exp(-epsilon * (z[k] - lower));
    const midS = envS[k] + (sp - envS[k]) * half, midQ = envQ[k] + (wq - envQ[k]) * half;
    plumeState(midS, midQ, z[k], p[k], guess);
    guess = plumeAir.T;
    const idx = k * C + i, air = Math.max(0, s.q[idx]), cloud = Math.max(0, s.qc[idx]);
    const environment = T[k] * (1 + virtual * air - loading * cloud), rising = plumeAir.T * (1 + virtual * (midQ - plumeAir.liquid) - loading * plumeAir.liquid);
    let work = R * (rising - environment) * dp[k] / p[k];
    const buoyancy = g * (rising - environment) / environment;
    if (plumeAir.liquid > 0) cloudy = true;
    if (!cloudy && work < 0) inhibition -= work;
    if (cloudy && buoyancy > 0) buoyant = true;
    else if (buoyant && !Number.isFinite(unbuoyant)) unbuoyant = p[k];
    if (O.plumeCapeParcel === 'undilute' && entraining) {
      plumeState(sourceS, sourceQ, z[k], p[k], plumeAir.T);
      work = R * (plumeAir.T * (1 + virtual * (sourceQ - plumeAir.liquid)) - environment) * dp[k] / p[k];
    }
    if (mixes) {
      if (!(depth < Infinity)) { top = k; break; }
      const x = 2 * mixing * depth, decay = Math.exp(-x);
      w2 = w2 * decay + 2 * O.plumeAcceleration * buoyancy * depth * (x > 0 ? -Math.expm1(-x) / x : 1);
      if (!(w2 > 0)) { top = k; break; }
    }
    if (cloudy && work > 0) cape += work;
    below = buoyancy;
    const full = Math.exp(-epsilon * depth);
    sp = envS[k] + (sp - envS[k]) * full; wq = envQ[k] + (wq - envQ[k]) * full;
    if (mixes && O.plumeConversion === 'sundqvist') {
      plumeState(sp, wq, upper, pi[i] * levels[k], guess);
      if (plumeAir.liquid > (landMask[i] ? IFS_PRECIPITATION.landThreshold : IFS_PRECIPITATION.seaThreshold)) {
        const alpha = liquidFraction(plumeAir.T, liquidT, iceT), P = IFS_PRECIPITATION;
        const bergeron = plumeAir.T < P.bergeron ? 1 + 0.5 * Math.sqrt(Math.min(P.bergeron - plumeAir.T, P.bergeron - P.ice)) : 1;
        const speed = Math.min(P.speed, Math.max(O.plumeVelocity, Math.sqrt(w2))), critical = P.critical / bergeron;
        const rate = P.conversion * (P.liquidFactor * alpha + 1 - alpha) * bergeron / (P.velocityScale * speed) * -Math.expm1(-((plumeAir.liquid / critical) ** 2));
        const fallen = -plumeAir.liquid * Math.expm1(-rate * depth); wq -= fallen; sp += L * fallen;
        if (mixedPlume && plumeAir.ice > 0) { const frozen = fallen * plumeAir.ice / plumeAir.liquid; sp += FUSION_HEAT * frozen; }
      }
    } else if (mixes) {
      plumeState(sp, wq, upper, pi[i] * levels[k], guess);
      const excess = plumeAir.liquid - O.plumeRainThreshold;
      if (excess > 0) { const fallen = -excess * Math.expm1(-O.plumeRainRate * depth); wq -= fallen; sp += L * fallen; if (mixedPlume && plumeAir.ice > 0) { const frozen = fallen * plumeAir.ice / plumeAir.liquid; sp += FUSION_HEAT * frozen; } }
    }
  }
  if (top === 0) top = 1;
  const deepTop = cloudy && top >= 1 && (O.convectionType === 'testParcel' || (O.convectionType === 'cloudDepth' ? pi[i] * (levels[base] - levels[top]) > DEEP_CLOUD_DEPTH : levels[top] * DEEP_REFERENCE < O.shallowTop));
  return { status: !cloudy ? 'notCloudy' : deepTop ? 'candidate' : 'shallow', source, top, cape, inhibition, lclPressure: lcl.pressure, base, buoyant, unbuoyant, topPressure: top >= 0 ? pi[i] * levels[top] : NaN };
}
function fillEnvironment(i, s) {
  for (let k = 0; k < K; k++) {
    const idx = k * C + i, cloud = Math.max(0, s.qc[idx]);
    env.T[k] = s.theta[idx] * exnerLayer[idx]; env.p[k] = pi[i] * sigmaMid[k]; env.dp[k] = pi[i] * dSigma[k]; env.z[k] = geopotential[idx] / g;
    env.s[k] = cp * env.T[k] + g * env.z[k] - L * cloud; env.q[k] = Math.max(0, s.q[idx]) + cloud;
  }
}
const openOf = (i) => (O.deckVeto ? ramp((DECK_CLOSED - radiation.mlmGate[i]) / (DECK_CLOSED - DECK_OPEN)) : 1) * (O.coupledVeto && bl.regime[i] === COUPLED_REGIME ? 0 : 1);
const bin100 = (p) => Math.min(9, Math.max(0, Math.floor(p / 100e2)));

const deepHeat = new Float64Array(K), deepWater = new Float64Array(K), lsCondensed = new Float64Array(K), uniformLayer = new Uint8Array(K), replayQt = new Float64Array(K);
function replayWater(n, i, s, term) {
  for (let k = 0; k < K; k++) { const idx = k * C + i, now = s.q[idx] + s.qc[idx]; water(n, term, k, now - replayQt[k]); replayQt[k] = now; }
}
moist.trace.convection = new Float64Array(K * C);
moist.trace.largeScale = new Float64Array(K * C);
const before = { precipitation: new Float64Array(C), convective: new Float64Array(C) };
const cellRain = new Float64Array(nT);
const replayed = Array.from({ length: nT }, () => ({ theta: new Float64Array(K), q: new Float64Array(K), qc: new Float64Array(K), rained: 0, convected: 0 }));
moist.adjust = (st, iFrom, iTo, step) => {
  const [sa, sb, sc, sd] = scratch;
  traced.forEach((i, n) => {
    const b = boxOf[i], A = acc[b], a = area[i];
    copyColumn(i, live, sa);
    copyColumn(i, live, sd);
    core.diagnoseColumn(i, pi, sa.theta, sa.q, sa.qc);
    hold(i, sa);
    for (let k = 0; k < K; k++) replayQt[k] = sa.q[k * C + i] + sa.qc[k * C + i];
    replica.condenseColumn(i, pi, sa.theta, sa.q, sa.qc);
    replayWater(n, i, sa, 'condensation');
    for (let k = 0; k < K; k++) {
      const idx = k * C + i, below = boundaryTop !== null && geopotential[idx] / g < boundaryTop[i];
      lsCondensed[k] = sa.qc[idx] - columnQc[k];
      const held = O.boundaryCondensation === 'uniform' || (O.boundaryCondensation === 'cloudLayer' && boundaryTop !== null && k <= bl.cloudLayer[i]);
      uniformLayer[k] = O.condensation === 'uniform' && (held || !below) ? 1 : 0;
      add(n, uniformLayer[k] ? 'condensationUniform' : 'condensationMixed', k, (sa.theta[idx] - columnT[k]) * exnerLayer[idx] + mixedPhase[n * K + k]);
      mixedPhase[n * K + k] = 0;
      columnT[k] = sa.theta[idx]; columnQ[k] = sa.q[idx]; columnQc[k] = sa.qc[idx];
    }
    copyColumn(i, sa, sb);
    fillEnvironment(i, sa);
    const open = openOf(i);
    const dilute = open > 0 ? ascend(i, sa, true) : { status: 'deck' };
    const pure = dilute.status !== 'deck' && dilute.status !== 'noLcl' ? ascend(i, sa, false) : null;
    deepOnly.plumeColumn(i, pi, sb.theta, sb.q, sb.qc, step, wind);
    if (O.convectionType === 'testParcel' && dilute.status === 'candidate' && !deepOnly.deep.parcelDeep) dilute.status = 'shallow';
    if (deepOnly.deep.parcelDeep) A.plume.typed += a;
    for (let k = 0; k < K; k++) {
      const idx = k * C + i;
      deepHeat[k] = deepOnly.deep.deep ? (sb.theta[idx] - sa.theta[idx]) * exnerLayer[idx] : 0;
      deepWater[k] = deepOnly.deep.deep ? sb.q[idx] + sb.qc[idx] - (sa.q[idx] + sa.qc[idx]) : 0;
    }
    const produced = replica.plumeColumn(i, pi, sa.theta, sa.q, sa.qc, step, wind);
    const deep = replica.deep.deep;
    if (deep !== deepOnly.deep.deep || (deep && replica.deep.baseFlux !== deepOnly.deep.baseFlux)) checks.deepFluxMismatch++;
    if (deep) {
      const baseFlux = deepOnly.deep.baseFlux, dd = deepOnly.deep;
      for (let k = dd.top; k < dd.base; k++) {
        const leaving = Math.max(0, deepOnly.cumulusFlux[k + 1] - deepOnly.cumulusFlux[k]) * deepOnly.plumeCarried[k + 1] * baseFlux * step;
        A.detrained[k] += a * leaving; A.detrainedPressure[k] += a * leaving * pi[i] * sigmaMid[k];
        if (k > dd.top && k < dilute.source) A.deepMade += a * deepOnly.cumulusFlux[k] * deepOnly.plumeRain[k] * baseFlux * step;
      }
      for (let k = 0; k < K; k++) {
        const mass = pi[i] * dSigma[k] / g, made = k > deepOnly.deep.top && k < dilute.source ? deepOnly.cumulusFlux[k] * deepOnly.plumeRain[k] * baseFlux * step : 0;
        const evaporated = made - deepOnly.convectiveFall[k];
        const freezing = mixedPlume && k > deepOnly.deep.top && k < dilute.source ? FUSION_HEAT * deepOnly.cumulusFlux[k] * deepOnly.plumeFrozen[k] * baseFlux * step : 0;
        const melting = mixedPlume ? -FUSION_HEAT * deepOnly.convectiveMelted[k] * baseFlux * step : 0;
        add(n, 'deepRain', k, L * made / (cp * mass));
        add(n, 'deepFreezing', k, freezing / (cp * mass));
        add(n, 'deepMelting', k, melting / (cp * mass));
        add(n, 'deepDowndraftEvaporation', k, -L * evaporated / (cp * mass));
        add(n, 'deepTransport', k, deepHeat[k] - (L * (made - evaporated) + freezing + melting) / (cp * mass));
      }
    }
    for (let k = 0; k < K; k++) {
      const idx = k * C + i, dT = (sa.theta[idx] - columnT[k]) * exnerLayer[idx];
      add(n, 'shallow', k, dT - deepHeat[k]);
      columnT[k] = sa.theta[idx]; columnQ[k] = sa.q[idx]; columnQc[k] = sa.qc[idx];
      const now = sa.q[idx] + sa.qc[idx];
      water(n, 'deep', k, deepWater[k]); water(n, 'shallow', k, now - replayQt[k] - deepWater[k]); replayQt[k] = now;
    }
    if (replica.cumulusBaseFlux[i] > 0) { replica.condenseColumn(i, pi, sa.theta, sa.q, sa.qc); heated(i, sa, n, 'recondensation'); replayWater(n, i, sa, 'condensation'); }
    const iced = ice && ice[i] > 0 ? (seaIce.concentration[i] > 0 ? seaIce.concentration[i] : 1) : 0;
    const stream = deep ? replica.convectiveFall : null;
    reserve.fill(0);
    if (deep) for (let k = K - 2; k >= 0; k--) reserve[k] = Math.max(0, reserve[k + 1] - replica.convectiveFall[k + 1]);
    copyColumn(i, sa, sc);
    const frozenStream = deep && mixedPlume ? replica.convectiveFrozen : null;
    const shadowRained = shadowFall(replica, i, sc, stream ? Float64Array.from(stream) : null, iced, frozenStream ? Float64Array.from(frozenStream) : null);
    const rained = replica.autoconvertColumn(i, pi, sa.theta, sa.q, sa.qc, step, stream, iced, frozenStream);
    let differs = shadowRained !== rained;
    for (let k = 0; k < K && !differs; k++) { const idx = k * C + i; if (sc.theta[idx] !== sa.theta[idx] || sc.q[idx] !== sa.q[idx] || sc.qc[idx] !== sa.qc[idx]) differs = true; }
    if (differs) checks.shadowMismatch++;
    for (let k = 0; k < K; k++) {
      const idx = k * C + i, mass = pi[i] * dSigma[k] / g, now = sa.q[idx] + sa.qc[idx];
      add(n, 'largeScaleEvaporation', k, -L * fall.evaporated[k] / (cp * mass));
      add(n, 'convectiveEvaporation', k, -(L * fall.convective[k] + FUSION_HEAT * fall.sublimated[k]) / (cp * mass));
      add(n, 'deepMelting', k, fall.refunded[k]);
      const evaporated = fall.evaporated[k] / mass, streamed = fall.convective[k] / mass, converted = -fall.converted[k] / mass;
      water(n, 'largeScaleEvaporation', k, evaporated); water(n, 'convectiveEvaporation', k, streamed); water(n, 'conversion', k, converted);
      water(n, 'iceFall', k, now - replayQt[k] - evaporated - streamed - converted); replayQt[k] = now;
    }
    hold(i, sa);
    const convected = deep ? replica.falling.convective + replica.deep.shallowRain : produced;
    if (replica.falling.moved) {
      replica.condenseColumn(i, pi, sa.theta, sa.q, sa.qc);
      for (let k = 0; k < K; k++) { const idx = k * C + i, dT = (sa.theta[idx] - columnT[k]) * exnerLayer[idx]; add(n, dT > 0 ? 'afterFallPositive' : 'afterFallNegative', k, dT); }
      replayWater(n, i, sa, 'condensation');
    }
    replica.fillColumn(i, pi, sa.q);
    replica.fillColumn(i, pi, sa.qc);
    replayWater(n, i, sa, 'filler');
    const out = replayed[n];
    for (let k = 0; k < K; k++) { const idx = k * C + i; out.theta[k] = sa.theta[idx]; out.q[k] = sa.q[idx]; out.qc[k] = sa.qc[idx]; }
    out.rained = rained; out.convected = convected;
    A.rain.deep += a * (deep ? replica.falling.convective : 0);
    A.rain.shallow += a * (deep ? replica.deep.shallowRain : produced);
    A.rain.liquidLow += a * fall.tags.low; A.rain.liquidHigh += a * fall.tags.high; A.rain.iceMelted += a * fall.tags.ice; A.rain.iceGround += a * fall.ground;

    const P = A.plume;
    if (dilute.status === 'deck') P.deck += a;
    else if (dilute.status === 'noLcl') P.noLcl += a;
    else if (dilute.status === 'notCloudy') P.notCloudy += a;
    else if (dilute.status === 'shallow') P.shallowTop[bin100(dilute.topPressure)] += a;
    if (dilute.lclPressure) { P.cloudBase += a * dilute.lclPressure; P.cloudBaseN += a; }
    if (dilute.status === 'candidate' || dilute.status === 'shallow') {
      P.diluteN += a; P.diluteCape += a * dilute.cape; P.pureCape += a * (pure ? pure.cape : 0);
      if (dilute.buoyant) P.dilute[bin100(Number.isFinite(dilute.unbuoyant) ? dilute.unbuoyant : dilute.topPressure)] += a; else P.neverBuoyant += a;
      if (pure && pure.buoyant) P.undilute[bin100(Number.isFinite(pure.unbuoyant) ? pure.unbuoyant : pure.topPressure)] += a; else P.pureNeverBuoyant += a;
    }
    if (dilute.status === 'candidate') {
      checks.probeChecked++;
      const d = deepOnly.deep;
      if (d.top !== dilute.top || d.cape !== dilute.cape || d.inhibition !== dilute.inhibition) checks.probeMismatch++;
      P.candidateCape[Math.min(7, Math.floor(d.cape / 60))] += a;
      if (d.deep) {
        P.fired += a; P.firedCape += a * d.cape; P.firedInhibition += a * d.inhibition; P.firedFlux += a * d.baseFlux; P.firedTop[bin100(dilute.topPressure)] += a; P.consumption += a * d.consumption;
        const relaxed = bechtold ? (d.consumptionP > 0 ? Math.max(0, d.pcape - d.pcapeBoundary) / (d.tau * d.consumptionP) : 0) : d.consumption > 0 && d.cape > O.plumeCape ? (d.cape - O.plumeCape) / (O.plumeRelaxation * d.consumption) : 0;
        const gate = ramp(0.5 + (O.inhibitionThreshold - d.inhibition) / Math.max(1, O.inhibitionThreshold));
        if (d.baseFlux < open * gate * relaxed * (1 - 1e-12)) P.limited += a;
        if (bechtold) { P.pcape += a * d.pcape; P.pcapeBoundary += a * d.pcapeBoundary; P.tau += a * d.tau; P.speed += a * d.speed; P.taus.push([d.tau, a]); }
      } else if (bechtold && d.consumptionP > 0 && !(d.pcape > d.pcapeBoundary)) { P.closedByForcing += a; P.closedByForcingSum += a * d.cape; }
      else if (!bechtold && d.cape <= O.plumeCape) { P.weakCape += a; P.weakCapeSum += a * d.cape; }
      else { P.closed += a; P.closedCape += a * d.cape; P.closedInhibition += a * d.inhibition; }
    }

    const Q = A.counter;
    adjusted.condenseColumn(i, pi, sd.theta, sd.q, sd.qc);
    adjusted.plumeColumn(i, pi, sd.theta, sd.q, sd.qc, step, wind);
    const counterFires = adjusted.deep.deep;
    const low = [0, 0, 0, 0];
    for (let k = 0; k < K; k++) {
      if (!(pi[i] * sigmaMid[k] > 700e2)) continue;
      const mass = pi[i] * dSigma[k] / g, u = uniformLayer[k] ? 0 : 2;
      low[u] += (uniformLayer[k] ? lsCondensed[k] : lsCondensed[k] + markQc[n * K + k] - preMixQc[n * K + k]) * mass;
      low[u + 1] += fall.converted[k];
    }
    for (let m = 0; m < 4; m++) Q.low[m] += a * low[m];
    if (dilute.lclPressure) for (let k = 0; k < K; k++) {
      if (!uniformLayer[k] || !(lsCondensed[k] > 0)) continue;
      const pk = pi[i] * sigmaMid[k], water = lsCondensed[k] * pi[i] * dSigma[k] / g;
      if (pk > dilute.lclPressure) Q.subcloud += a * water; else if (pk > 700e2) Q.cloudLayer += a * water;
    }
    if (counterFires) { Q.fires += a; Q.flux += a * adjusted.deep.baseFlux; }
    if (deep) Q.fluxActual += a * replica.deep.baseFlux;
    if (counterFires && deep) Q.both += a;
    if (deep && !counterFires) Q.actualOnly += a;
    if (counterFires && !deep) {
      Q.firesOnly += a; Q.weight += a;
      for (let m = 0; m < 4; m++) Q.lowOnly[m] += a * low[m];
      for (let k = 0; k < K; k++) {
        const idx = k * C + i, T = live.theta[idx] * exnerLayer[idx], p = pi[i] * sigmaMid[k];
        Q.temperature[k] += a * T; Q.humidity[k] += a * Math.max(0, live.q[idx]) / saturationHumidity(T, p); Q.pressure[k] += a * p;
      }
    }
  });
  traced.forEach((i) => { before.precipitation[i] = moist.precipitation[i]; before.convective[i] = moist.convectivePrecipitation[i]; });
  moist.trace.convection.fill(0); moist.trace.largeScale.fill(0);
  moistAdjust(st, iFrom, iTo, step);
  traced.forEach((i, n) => {
    const out = replayed[n];
    checks.replicaColumns++;
    let differs = false;
    for (let k = 0; k < K; k++) { const idx = k * C + i; if (out.theta[k] !== theta[idx] || out.q[k] !== q[idx] || out.qc[k] !== qc[idx]) differs = true; }
    if (differs) checks.replicaMismatch++;
    const fell = moist.precipitation[i] - before.precipitation[i];
    if (Math.abs(fell - (out.rained + out.convected)) > 1e-12 * Math.max(1, fell)) checks.precipOff++;
    const b = boxOf[i], a = area[i];
    const convectiveFell = moist.convectivePrecipitation[i] - before.convective[i];
    cellRain[n] += fell;
    acc[b].rain.model += a * fell; acc[b].rain.modelConvective += a * convectiveFell;
    acc[b].localHour[Math.floor((((12 + ((model.time + step / 2) % 86400) / 3600 + lon[i] / 15) % HOURS) + HOURS) % HOURS)] += a * convectiveFell;
    let traceSum = 0, partsSum = 0;
    for (let k = 0; k < K; k++) {
      const s = (n * NT) * K;
      const parts = ['deepRain', 'deepFreezing', 'deepMelting', 'deepDowndraftEvaporation', 'deepTransport', 'shallow', 'convectiveEvaporation'].reduce((x, t) => x + stepHeat[s + T_[t] * K + k], 0);
      traceSum += Math.abs(moist.trace.convection[k * C + i]); partsSum += Math.abs(moist.trace.convection[k * C + i] - parts);
    }
    if (partsSum > 1e-9 * Math.max(1, traceSum)) checks.traceOff++;
    if (b === 0) {
      const fired = (moist.convectivePrecipitation[i] - before.convective[i]) * 86400 / step > 1;
      audit.area += a;
      if (fired) audit.fired += a;
      for (let k = 0; k < K; k++) {
        audit.pressure[k] += a * pi[i] * sigmaMid[k];
        if (!fired) continue;
        audit.conv[k] += a * moist.trace.convection[k * C + i];
        const s = (n * NT) * K;
        ['deepRain', 'deepDowndraftEvaporation', 'deepTransport', 'shallow', 'convectiveEvaporation'].forEach((t, m) => { audit.parts[m * K + k] += a * stepHeat[s + T_[t] * K + k]; });
      }
    }
    for (let k = 0; k < K; k++) { const s = n * K + k, idx = k * C + i; mark[s] = theta[idx]; markQ[s] = q[idx]; markQc[s] = qc[idx]; }
  });
};

model.restartPrecipitation();
for (let step = 0; step < STEPS; step++) {
  stepHeat.fill(0); stepWater.fill(0);
  traced.forEach((i, n) => {
    exnerOf(i, ex);
    const b = boxOf[i], A = acc[b], a = area[i];
    for (let k = 0; k < K; k++) {
      const idx = k * C + i, T = theta[idx] * ex[k], p = pi[i] * sigmaMid[k];
      stepStart[n * K + k] = T; stepStartQt[n * K + k] = q[idx] + qc[idx];
      A.temperature[k] += a * T; A.humidity[k] += a * Math.max(0, q[idx]) / saturationHumidity(T, p);
    }
  });
  model.step(dt);
  traced.forEach((i, n) => {
    exnerOf(i, ex);
    const b = boxOf[i], A = acc[b], a = area[i];
    A.areaSteps += a;
    for (let k = 0; k < K; k++) {
      const idx = k * C + i, s = n * K + k, mass = pi[i] * dSigma[k] / g, total = theta[idx] * ex[k] - stepStart[s];
      let sum = 0;
      for (let t = 0; t < NT; t++) {
        const dT = stepHeat[(n * NT + t) * K + k];
        sum += dT;
        A.heat[t * K + k] += a * dT;
        A.energy[t * K + k] += a * cp * mass * dT;
      }
      checks.closure = Math.max(checks.closure, Math.abs(sum - total));
      A.total[k] += a * total;
      let waterSum = 0;
      for (let w = 0; w < NW; w++) {
        const dq = stepWater[(n * NW + w) * K + k];
        waterSum += dq;
        A.water[w * K + k] += a * dq;
        A.waterMass[w * K + k] += a * mass * dq;
      }
      const waterTotal = q[idx] + qc[idx] - stepStartQt[s];
      checks.waterClosure = Math.max(checks.waterClosure, Math.abs(waterSum - waterTotal));
      A.totalWater[k] += a * waterTotal;
      A.pressure[k] += a * pi[i] * sigmaMid[k];
      A.thickness[k] += a * pi[i] * dSigma[k];
      A.height[k] += a * (geopotential[idx] / g - zs[i]);
    }
  });
  if ((step + 1) % 32 === 0) console.error(`step ${step + 1} of ${STEPS} (${((performance.now() - t0) / 1000).toFixed(0)} s)`);
}

const f = (x, d) => (Number.isFinite(x) ? x.toFixed(d) : 'n/a');
const globalArea = area.reduce((x, y) => x + y, 0);
const perDay = 86400 / dt;
say(`tropical heating of ${FILE.split('/').pop()}: day ${saved.day}, N=${saved.N}, K=${K}; CPU, ocean off; ${STEPS} steps of ${dt} s${Object.keys(MOIST).length ? `; MOIST ${JSON.stringify(MOIST)}` : ''}${Object.keys(RADIATION).length ? `; RADIATION ${JSON.stringify(RADIATION)}` : ''}${Object.keys(BOUNDARY_LAYER).length ? `; BOUNDARY_LAYER ${JSON.stringify(BOUNDARY_LAYER)}` : ''}${Object.keys(SURFACE).length ? `; SURFACE ${JSON.stringify(SURFACE)}` : ''}`);
say(`checks: the replayed moist step differs from the model's in ${checks.replicaMismatch} of ${checks.replicaColumns} column-steps (precipitation in ${checks.precipOff}); the copy of the rain's fall in ${checks.shadowMismatch}; the deep-only plume's base flux in ${checks.deepFluxMismatch}; the probe's top, CAPE or inhibition in ${checks.probeMismatch} of ${checks.probeChecked} candidates; the convective trace against its terms in ${checks.traceOff}; the terms close on the temperature change to ${checks.closure.toExponential(1)} K a step and on the change of q_t to ${checks.waterClosure.toExponential(1)} kg/kg a step (q_t changed outside a water term by ${checks.waterUnassigned.toExponential(1)}); last step's global sensible heat ${f(checks.sensibleGlobal / globalArea, 4)} against the model's ${f(checks.sensibleModel / globalArea, 4)} W/m2`);
const jordan = [[1015, 26.3, 18.9], [1000, 25.4, 18.4], [950, 22.2, 16.6], [900, 19.6, 14.4], [850, 17.3, 12.4], [800, 14.8, 10.5], [750, 12.0, 8.8], [700, 8.9, 7.4], [650, 5.6, 6.0], [600, 1.8, 4.7], [550, -2.3, 3.6], [500, -6.9, 2.6], [450, -11.9, 1.8], [400, -17.3, 1.2], [350, -24.4, 0.66], [300, -32.4, 0.35], [250, -42.1, 0.14], [200, -54.0, 0.04], [175, -60.2, 0.02], [150, -67.2, 0.01], [125, -73.0, 0.005], [100, -73.3, 0.003]];
function jordanAt(p) {
  const hPa = p / 100;
  if (hPa >= jordan[0][0] || hPa < jordan[jordan.length - 1][0]) return null;
  for (let n = 0; n < jordan.length - 1; n++) {
    const [p0j, t0j, r0j] = jordan[n], [p1j, t1j, r1j] = jordan[n + 1];
    if (hPa <= p0j && hPa >= p1j) {
      const w = Math.log(hPa / p0j) / Math.log(p1j / p0j), T = 273.15 + t0j + w * (t1j - t0j), r = 1e-3 * (r0j + w * (r1j - r0j));
      return { T, rh: r * p / (0.622 + r) / saturationVaporPressure(T) };
    }
  }
  return null;
}
function weightedMedian(pairs) {
  if (!pairs.length) return NaN;
  const sorted = [...pairs].sort((x, y) => x[0] - y[0]), half = sorted.reduce((s, [, w]) => s + w, 0) / 2;
  let run = 0;
  for (const [v, w] of sorted) { run += w; if (run >= half) return v; }
  return sorted[sorted.length - 1][0];
}
const PHYSICS_WATER = WATER.filter((t) => t !== 'dynamics' && t !== 'closure');
acc.forEach((A, b) => {
  const [name, , surfaceKind] = BOX_LIST[b], S = A.areaSteps;
  const p = Float64Array.from(A.pressure, (x) => x / S), z = Float64Array.from(A.height, (x) => x / S), dp = Float64Array.from(A.thickness, (x) => x / S);
  const kday = (t, k) => A.heat[t * K + k] / S * perDay, wm2 = (t, k) => A.energy[t * K + k] / S / dt;
  const moistening = (w, k) => L / cp * A.water[w * K + k] / S * perDay, waterWm2 = (w) => Array.from({ length: K }, (_, k) => L * A.waterMass[w * K + k] / S / dt).reduce((x, y) => x + y, 0);
  say(`\n=== ${name} (${f(S / STEPS / 1e12, 2)} Mkm2) ===`);
  const r = A.rain, mm = (x) => x / S * perDay;
  const ls = mm(r.liquidLow + r.liquidHigh + r.iceMelted + r.iceGround);
  say(`surface rain, mm/d: ${f(mm(r.model), 2)} (the model's; convective ${f(mm(r.modelConvective), 2)}, share ${f(r.modelConvective / r.model, 2)}): deep plume ${f(mm(r.deep), 2)}, shallow plume ${f(mm(r.shallow), 2)}, large-scale ${f(ls, 2)} = liquid converted below 700 hPa ${f(mm(r.liquidLow), 2)} + above 700 hPa ${f(mm(r.liquidHigh), 2)} + falling ice melted ${f(mm(r.iceMelted), 2)} + reaching the ground as ice ${f(mm(r.iceGround), 2)}; stratiform share (melted falling ice and conversion above 700 hPa over all the rain) ${f((r.iceMelted + r.liquidHigh) / r.model, 3)}`);
  const physicsTerms = TERMS.filter((t) => !DYNAMICS.has(t));
  const Q1 = Float64Array.from({ length: K }, (_, k) => physicsTerms.reduce((s, t) => s + kday(T_[t], k), 0));
  const QR = Float64Array.from({ length: K }, (_, k) => kday(T_.shortwave, k) + kday(T_.longwave, k));
  const Q1R = Float64Array.from(Q1, (x, k) => x - QR[k]);
  const deepAll = Float64Array.from({ length: K }, (_, k) => ['deepRain', 'deepFreezing', 'deepMelting', 'deepDowndraftEvaporation', 'deepTransport', 'convectiveEvaporation'].reduce((s, t) => s + kday(T_[t], k), 0));
  const Q2 = Float64Array.from({ length: K }, (_, k) => -PHYSICS_WATER.reduce((s, t) => s + moistening(W_[t], k), 0));
  say('profile, K/day (box means; Q1 all physics, Q1-QR without radiation, Q2 = -(L/cp) the physics\' change of q_t; dyn the RK4 step, clos the closure, diss the kinetic energy returned as heat; condBL the condensation in the layers below the mixing top (to saturation), condU the uniform condensation in those above it, each with the cloud the boundary layer\'s mixing evaporated there, dRain the latent heat of the deep plume\'s rain, dDDev its downdraft\'s evaporation, dTrans its two drafts\' transport of s_l, shal the shallow plume, recond the condensation after the plumes, LSev and CVev the evaporation of large-scale and deep-plume rain, fall+ and fall- the adjustment after the ice fall, snow its fusion heat, dry the dry adjustment):');
  say(`  ${'hPa'.padStart(5)} ${'z m'.padStart(6)} ${'Q1'.padStart(6)} ${'Q1-QR'.padStart(6)} ${'Q2'.padStart(6)} ${SHORT.map((s) => s.padStart(6)).join(' ')} ${'total'.padStart(6)}`);
  for (let k = 0; k < K; k++) {
    if (p[k] < 70e2) continue;
    say(`  ${f(p[k] / 100, 0).padStart(5)} ${f(z[k], 0).padStart(6)} ${f(Q1[k], 2).padStart(6)} ${f(Q1R[k], 2).padStart(6)} ${f(Q2[k], 2).padStart(6)} ${TERMS.map((t) => f(kday(T_[t], k), 2).padStart(6)).join(' ')} ${f(A.total[k] / S * perDay, 2).padStart(6)}`);
  }
  say('profile, W/m2 per layer:');
  say(`  ${'hPa'.padStart(5)} ${SHORT.map((s) => s.padStart(6)).join(' ')}`);
  for (let k = 0; k < K; k++) {
    if (p[k] < 70e2) continue;
    say(`  ${f(p[k] / 100, 0).padStart(5)} ${TERMS.map((t) => f(wm2(T_[t], k), 1).padStart(6)).join(' ')}`);
  }
  say(`  ${'sum'.padStart(5)} ${TERMS.map((t) => f(Array.from({ length: K }, (_, k) => wm2(T_[t], k)).reduce((x, y) => x + y, 0), 1).padStart(6)).join(' ')}`);
  say('moistening by process, K/day (L/cp dq_t/dt, positive moistens; dyn the RK4 step, clos the closure, evap the surface evaporation, BL the boundary layer\'s mixing, cond the condensation (no change of q_t), deep and shal the deep and shallow plumes, LSev and CVev the evaporation of large-scale and deep-plume rain, conv the condensate converted to rain, ice the ice fall with its melting, fill the filler, dry the dry adjustment):');
  say(`  ${'hPa'.padStart(5)} ${WATER_SHORT.map((s) => s.padStart(6)).join(' ')} ${'total'.padStart(6)}`);
  for (let k = 0; k < K; k++) {
    if (p[k] < 70e2) continue;
    say(`  ${f(p[k] / 100, 0).padStart(5)} ${WATER.map((t) => f(moistening(W_[t], k), 2).padStart(6)).join(' ')} ${f(L / cp * A.totalWater[k] / S * perDay, 2).padStart(6)}`);
  }
  say(`  column, mm/d: ${WATER.map((t, w) => `${WATER_SHORT[w]} ${f(waterWm2(w) / L * 86400, 2)}`).join(', ')}`);
  const column = (terms) => terms.map((t) => [t, Array.from({ length: K }, (_, k) => wm2(T_[t], k)).reduce((x, y) => x + y, 0)]);
  const q1rParts = column(physicsTerms.filter((t) => t !== 'shortwave' && t !== 'longwave')), q2Parts = PHYSICS_WATER.map((t) => [t, -waterWm2(W_[t])]);
  const listed = (parts) => parts.filter(([, v]) => Math.abs(v) >= 0.05).map(([t, v]) => `${SHORT[T_[t]] ?? WATER_SHORT[W_[t]]} ${f(v, 1)}`).join(', ');
  say(`column integrals, W/m2: Q1 ${f(column(physicsTerms).reduce((x, [, v]) => x + v, 0), 1)} (radiation ${f(column(['shortwave', 'longwave']).reduce((x, [, v]) => x + v, 0), 1)}); Q1-QR ${f(q1rParts.reduce((x, [, v]) => x + v, 0), 1)} = ${listed(q1rParts)}; Q2 ${f(q2Parts.reduce((x, [, v]) => x + v, 0), 1)} = ${listed(q2Parts.map(([t, v]) => [t, v]))}`);
  for (const [label, profile] of [['Q1', Q1], ['Q1-QR', Q1R], ['deep plume with its rain\'s evaporation', deepAll], ['Q2', Q2]]) {
    const pk = heatingProfile(profile, p, dp);
    say(`${label} peak below 100 hPa: layer at ${f(pk.layer / 100, 0)} hPa (${f(pk.value, 2)} K/d); over 50 hPa bins ${f(pk.bin[0] / 100, 0)}-${f(pk.bin[1] / 100, 0)} hPa (${f(pk.binValue, 2)} K/d); centroid of its positive part ${f(pk.centroid / 100, 0)} hPa`);
  }
  if (surfaceKind === 'land') {
    const total = A.localHour.reduce((x, y) => x + y, 0);
    let peak = 0;
    for (let h = 1; h < HOURS; h++) if (A.localHour[h] > A.localHour[peak]) peak = h;
    say(`convective rain by local solar hour (share of the box's ${f(mm(r.modelConvective), 2)} mm/d; the step's middle at each cell's longitude): ${Array.from(A.localHour, (x, h) => `${h} ${f(x / total, 3)}`).join(', ')}; largest at ${peak}-${peak + 1} h`);
  }
  {
    const P = A.plume, total = S;
    const tops = Array.from(P.shallowTop, (x, n) => `${n * 100}-${n * 100 + 100} ${f(x / total, 3)}`).filter((x) => !x.endsWith(' 0.000')).join(', ');
    const firedTops = Array.from(P.firedTop, (x, n) => `${n * 100}-${n * 100 + 100} ${f(x / P.fired, 3)}`).filter((x) => !x.endsWith(' 0.000')).join(', ');
    const capes = Array.from(P.candidateCape, (x, n) => `${n * 60}-${n === 7 ? 'up' : n * 60 + 60} ${f(x / total, 3)}`).join(', ');
    say(`deep plume, share of the column-steps: deck veto ${f(P.deck / total, 3)}, no condensation level ${f(P.noLcl / total, 3)}, never cloudy ${f(P.notCloudy / total, 3)}, ${O.convectionType === 'testParcel' ? 'cloudy but its test parcel\'s cloud no deeper than 200 hPa' : O.convectionType === 'cloudDepth' ? 'cloudy but its cloud no deeper than 200 hPa' : 'cloudy but topping below 700 hPa'} ${f(P.shallowTop.reduce((x, y) => x + y, 0) / total, 3)} (tops by hPa: ${tops}); ${bechtold ? `deep candidates whose PCAPE the boundary-layer part PCAPE_bl takes up ${f(P.closedByForcing / total, 3)} (mean CAPE ${f(P.closedByForcingSum / P.closedByForcing, 0)} J/kg)` : `deep candidates with CAPE at most plumeCape ${O.plumeCape} ${f(P.weakCape / total, 3)} (mean CAPE ${f(P.weakCapeSum / P.weakCape, 0)} J/kg)`}, otherwise closed by inhibition or consumption ${f(P.closed / total, 3)} (CAPE ${f(P.closedCape / P.closed, 0)}, inhibition ${f(P.closedInhibition / P.closed, 1)} J/kg), fired ${f(P.fired / total, 3)}`);
    say(`  candidates' CAPE (J/kg, share of all column-steps): ${capes}`);
    say(`  fired: mean CAPE ${f(P.firedCape / P.fired, 0)} J/kg, inhibition ${f(P.firedInhibition / P.fired, 1)} J/kg, base flux ${f(P.firedFlux / P.fired, 4)} kg/m2/s (box mean ${f(P.firedFlux / total, 5)}), consumption F ${f(P.consumption / P.fired, 4)} J/kg per s per kg/m2/s, held below its closure by the boundary-loss or Courant limit on ${f(P.limited / P.fired, 3)}; tops by hPa: ${firedTops}`);
    if (bechtold) say(`  fired, the closure of Bechtold et al. (2014): mean PCAPE ${f(P.pcape / P.fired, 1)} Pa, PCAPE_bl ${f(P.pcapeBoundary / P.fired, 1)} Pa, w-bar ${f(P.speed / P.fired, 2)} m/s, tau ${f(P.tau / P.fired / 60, 1)} min (median ${f(weightedMedian(P.taus) / 60, 1)} min)`);
    say(`  mean condensation level ${f(P.cloudBase / P.cloudBaseN / 100, 0)} hPa; where the cloudy plume, once buoyant, first stops being buoyant (or its top), hPa, share of the cloudy plumes: entraining ${Array.from(P.dilute, (x, n) => `${n * 100}-${n * 100 + 100} ${f(x / P.diluteN, 3)}`).filter((x) => !x.endsWith(' 0.000')).join(', ')}, never buoyant ${f(P.neverBuoyant / P.diluteN, 3)}; the same plume without entrainment ${Array.from(P.undilute, (x, n) => `${n * 100}-${n * 100 + 100} ${f(x / P.diluteN, 3)}`).filter((x) => !x.endsWith(' 0.000')).join(', ')}, never buoyant ${f(P.pureNeverBuoyant / P.diluteN, 3)}; mean CAPE of the cloudy plumes ${f(P.diluteCape / P.diluteN, 0)} J/kg, without entrainment ${f(P.pureCape / P.diluteN, 0)} J/kg (rain-out at plumeRainRate, condensate loading)`);
  }
  {
    const Q = A.counter, total = S;
    say(`under saturation adjustment of the same column the deep plume fires on ${f(Q.fires / total, 3)} of the column-steps (mean base flux ${f(Q.flux / total, 5)} kg/m2/s over the box, against ${f(Q.fluxActual / total, 5)}): on ${f(Q.both / total, 3)} both fire, on ${f(Q.firesOnly / total, 3)} only it, on ${f(Q.actualOnly / total, 3)} only the model's`);
    const day = (x, w) => f(x / w * perDay, 2);
    say(`  below 700 hPa, kg/m2/day over the box: the uniform condensation (layers above the mixing top) condenses net ${day(Q.low[0], total)} and its cloud converts ${day(Q.low[1], total)} to rain; the layers below the mixing top (saturation adjustment, net of the boundary layer's mixing) ${day(Q.low[2], total)} and ${day(Q.low[3], total)}; of the uniform condensation's gross condensate ${f(Q.subcloud / total * perDay, 3)} forms below the plume's condensation level and ${day(Q.cloudLayer, total)} between it and 700 hPa`);
    say(`  in the columns where only the adjusted plume fires, per unit of their area: uniform ${day(Q.lowOnly[0], Q.firesOnly)} condensed, ${day(Q.lowOnly[1], Q.firesOnly)} converted; below the mixing top ${day(Q.lowOnly[2], Q.firesOnly)} and ${day(Q.lowOnly[3], Q.firesOnly)}`);
  }
  say('sounding against Jordan (1958): hPa, box T - Jordan (K), box RH (over water) / Jordan\'s; the columns where only the adjusted plume fires in brackets');
  const rows = [];
  for (let k = 0; k < K; k++) {
    const J = jordanAt(p[k]);
    if (!J) continue;
    const T = A.temperature[k] / S, rh = A.humidity[k] / S, Q = A.counter;
    const sub = Q.weight > 0 ? ` (${f(Q.temperature[k] / Q.weight - J.T, 1)}, ${f(Q.humidity[k] / Q.weight, 2)})` : '';
    rows.push(`${f(p[k] / 100, 0)} ${f(T - J.T, 1)} ${f(rh, 2)}/${f(J.rh, 2)}${sub}`);
  }
  say(`  ${rows.join('; ')}`);
  const nearest = (hPa) => { let best = 0; for (let k = 1; k < K; k++) if (Math.abs(p[k] - 100 * hPa) < Math.abs(p[best] - 100 * hPa)) best = k; return best; };
  const sounding = Object.fromEntries([848, 704, 516, 439].flatMap((hPa) => { const k = nearest(hPa), J = jordanAt(p[k]); return [[`dT${hPa}`, A.temperature[k] / S - J.T], [`rh${hPa}`, A.humidity[k] / S]]; }));
  for (const hPa of [946, 963]) sounding[`rh${hPa}`] = A.humidity[nearest(hPa)] / S;
  sounding.rhLowest = A.humidity[K - 1] / S;
  const shallowExport = -Array.from({ length: K }, (_, k) => (p[k] > 950e2 ? A.waterMass[W_.shallow * K + k] : 0)).reduce((x, y) => x + y, 0) / S / dt * 86400;
  let wettest = 0;
  traced.forEach((i, n) => { if (boxOf[i] === b) wettest = Math.max(wettest, cellRain[n] / (STEPS * dt) * 86400); });
  const hourly = A.localHour.reduce((x, y) => x + y, 0) > 0 ? A.localHour.indexOf(Math.max(...A.localHour)) : NaN;
  const layerMean = (profile, top, bottomP) => { let sum = 0, mass = 0; for (let k = 0; k < K; k++) if (p[k] >= top && p[k] <= bottomP) { sum += profile[k] * dp[k]; mass += dp[k]; } return sum / mass; };
  const longwave = Float64Array.from({ length: K }, (_, k) => kday(T_.longwave, k)), q1r = heatingProfile(Q1R, p, dp), P = A.plume;
  const summary = { rain: mm(r.model), convectiveShare: r.modelConvective / r.model, firing: P.fired / S, q1rLayer: q1r.layer / 100, q1rBin: (q1r.bin[0] + q1r.bin[1]) / 200, q1rBinValue: q1r.binValue, q1rCentroid: q1r.centroid / 100, largeScaleBelow700: mm(r.liquidLow), stratiformShare: (r.iceMelted + r.liquidHigh) / r.model, longwave300to500: layerMean(longwave, 300e2, 500e2), dilutedCape: P.diluteCape / P.diluteN, undilutedCape: P.pureCape / P.diluteN, ...sounding, undiluteStop700: P.undilute[7] / P.diluteN, shallowExport950: shallowExport, wettestCell: wettest, convectivePeakHour: hourly, firedAbove300: P.firedTop.slice(0, 3).reduce((x, y) => x + y, 0) / P.fired, typedDeep: P.typed / S, ...(bechtold ? { tauMedian: weightedMedian(P.taus) / 60, pcape: P.pcape / P.fired, pcapeBoundary: P.pcapeBoundary / P.fired } : {}) };
  {
    const out = A.detrained.reduce((x, y) => x + y, 0), at = A.detrainedPressure.reduce((x, y) => x + y, 0), above = Array.from(A.detrained).reduce((x, y, k) => x + (p[k] < 400e2 ? y : 0), 0);
    summary.detrainedPerRain = out / r.modelConvective; summary.detrainedPerMade = out / A.deepMade; summary.detrainedCentroid = out > 0 ? at / out / 100 : NaN; summary.detrainedAbove400 = out > 0 ? above / out : NaN;
    say(`deep plume's detrained condensate (sum over layers of the mass flux's fall times the condensate it carried in, over the day): ${f(mm(out), 3)} mm/d, ${f(out / r.modelConvective, 3)} of the convective rain, ${f(out / A.deepMade, 3)} of the rain the plume made, its mass-weighted pressure ${f(summary.detrainedCentroid, 0)} hPa, ${f(summary.detrainedAbove400, 3)} of it above 400 hPa; by layer (hPa mm/d): ${Array.from(A.detrained, (x, k) => [p[k], mm(x)]).filter(([, v]) => v > 5e-4).map(([pk, v]) => `${f(pk / 100, 0)} ${f(v, 3)}`).join(', ')}`);
  }
  say(`summary ${name}: ${JSON.stringify(Object.fromEntries(Object.entries(summary).map(([key, v]) => [key, Number(v.toPrecision(6))])))}`);
});
{
  const p = Float64Array.from(audit.pressure, (x) => x / audit.area / 100);
  const conv = Float64Array.from(audit.conv, (x) => x / audit.fired * perDay);
  let peak = 0;
  for (let k = 0; k < K; k++) if (conv[k] > conv[peak]) peak = k;
  const part = (m, k) => audit.parts[m * K + k] / audit.fired * perDay;
  let upper = 0;
  for (let k = 0; k < K; k++) if (p[k] >= 400 && p[k] <= 520 && conv[k] > conv[upper]) upper = k;
  say(`\nthe firing columns' convective heating over the day (the convective trace, no mass weights), Pacific ITCZ: firing on ${f(audit.fired / audit.area, 3)} of the column-steps; peak at ${f(p[peak], 0)} hPa, ${f(conv[peak], 2)} K/d = deep rain ${f(part(0, peak), 2)} + downdraft evaporation ${f(part(1, peak), 2)} + deep transport ${f(part(2, peak), 2)} + shallow ${f(part(3, peak), 2)} + rain evaporation below cloud base ${f(part(4, peak), 2)}; largest at 400-520 hPa ${f(p[upper], 0)} hPa, ${f(conv[upper], 2)} K/d = ${f(part(0, upper), 2)} + ${f(part(1, upper), 2)} + ${f(part(2, upper), 2)} + ${f(part(3, upper), 2)} + ${f(part(4, upper), 2)}`);
  say(`  firing-column profile (hPa K/d): ${Array.from(conv, (x, k) => `${f(p[k], 0)} ${f(x, 1)}`).filter((_, k) => p[k] > 150).join(', ')}`);
  const deepOnlyProfile = Float64Array.from({ length: K }, (_, k) => part(0, k) + part(1, k) + part(2, k) + part(4, k));
  let deepPeak = -1;
  for (let k = 0; k < K; k++) if (p[k] >= 100 && (deepPeak < 0 || deepOnlyProfile[k] > deepOnlyProfile[deepPeak])) deepPeak = k;
  say(`  the same columns without the shallow plume (the deep plume and its rain's evaporation): peak at ${f(p[deepPeak], 0)} hPa, ${f(deepOnlyProfile[deepPeak], 2)} K/d; profile ${Array.from(deepOnlyProfile, (x, k) => `${f(p[k], 0)} ${f(x, 1)}`).filter((_, k) => p[k] > 150).join(', ')}`);
}
{
  let sum = 0, convective = 0, wettest = 0, at = -1;
  for (let i = 0; i < C; i++) {
    sum += area[i] * moist.precipitation[i]; convective += area[i] * moist.convectivePrecipitation[i];
    if (moist.precipitation[i] > wettest) { wettest = moist.precipitation[i]; at = i; }
  }
  const scale = 86400 / (STEPS * dt);
  say(`global rain over the ${STEPS} steps ${f(sum / globalArea * scale, 3)} mm/d (convective ${f(convective / globalArea * scale, 3)}); wettest cell ${f(wettest * scale, 1)} mm/d at ${f(mesh.latCell[at] * deg, 1)}, ${f(lon[at], 1)}`);
  say(`global summary: ${JSON.stringify({ rain: Number((sum / globalArea * scale).toPrecision(6)), convective: Number((convective / globalArea * scale).toPrecision(6)), wettestCell: Number((wettest * scale).toPrecision(6)) })}`);
}
say(`(${f((performance.now() - t0) / 1000, 0)} s)`);
