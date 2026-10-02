import { getDevice, storageBuffer, emptyBuffer, readBuffer, readRanges, reductionKernel, finishReduction, reductionGroups as groupsOf } from './device.module.js';
import { sigmaInterfaces, R_DRY, CP_DRY, P0, GRAVITY, VIRTUAL_FACTOR } from '../dynamics/sigmaCore.module.js';
import { sunDirection, nearestLayer, STABILITY_SIGMA, UNDECIDED, RAYLEIGH_BANDS, LAND_AEROSOL, SEA_AEROSOL, GREENHOUSE_GASES, OZONE_COLUMN, YEAR } from '../physics/radiation.module.js';
import { VAPOR_STRENGTH } from '../physics/shortwaveGases.module.js';
import { physicsConstants, PHYSICS_FUNCTIONS, PHYSICS_KERNELS } from './physics.gpu.js';
import { MOIST_DEFAULTS } from '../physics/moist.module.js';
import { SEA_DRAG } from '../physics/surface.module.js';

const MAX_EDGES = 6, MAX_EDGES_ON_EDGE = 10, WORKGROUP = 64, RING_SLOTS = 16384, MAXIMUM_SURFACE_PRESSURE = 110000;

/*
 * The hydrostatic sigma-coordinate core of sigmaCore.module.js on the
 * GPU, in single precision. The mesh lives in one integer and one float
 * buffer, the level constants in a third; the state (pi, theta, u,
 * surfaceT, q, qc, ice) is one f32 buffer with the CPU layout inside it,
 * and the RK4 stages, the trial state and the tendency use the same
 * layout. One tendency evaluation is eight dispatches — edge mass
 * fluxes, layer divergences, the column diagnosis with dπ/dt and πσ̇,
 * kite-weighted π on vertices, potential vorticity on vertices and on
 * edges, the cell tendencies (kinetic energy, geopotential plus kinetic
 * energy, flux-form transport of θ, q and qc), and the momentum
 * tendency with the RTSK PV flux, the two pressure-gradient terms, the
 * vertical advection, the aerodynamic drag on the lowest layer and the
 * top sponge — then a fused advance, and after the four stages a fused
 * combine. The ∇⁴ closures follow once per step as in the split CPU
 * core. Every kernel binds the same seven buffers in the same order.
 */
export function layoutFor(mesh, K, cumulusLayers = 0, momentumLayers = 0) {
  const C = mesh.nCells, E = mesh.nEdges, V = mesh.nVertices;
  const KC = K * C, KE = K * E, KV = K * V;
  const seq = (names) => { const out = {}; let off = 0; for (const [name, n] of names) { out[name] = off; off += n; } out.total = off; return out; };
  const MI = seq([['COE', 2 * E], ['VOE', 2 * E], ['EOC', MAX_EDGES * C], ['ESC', MAX_EDGES * C], ['COC', MAX_EDGES * C], ['NEC', C], ['COV', 3 * V], ['EOV', 3 * V], ['ESV', 3 * V], ['EOE', MAX_EDGES_ON_EDGE * E], ['NEE', E]]);
  const MF = seq([['AREA', C], ['ATRI', V], ['DC', E], ['DV', E], ['FV', V], ['KAV', 3 * V], ['PVW', MAX_EDGES_ON_EDGE * E], ['NEDGE', 3 * E], ['LAT', C], ['XC', 3 * C], ['GPHIS', E], ['PHIS', C]]);
  const LV = seq([['SL', K], ['SU', K], ['DS', K], ['SM', K], ['TOP', K], ['CL', K], ['CM', K], ['CD', K], ['CA', K], ['CB', K], ['CT', K], ['GR', K], ['GABS', K], ['SHAPE', K], ['OZ', K], ['GASE', K], ['AER', K], ['OZS', K]]);
  const S = seq([['PI', C], ['TH', KC], ['U', KE], ['TS', C], ['Q', KC], ['QC', KC], ['ICE', C]]);
  const D = seq([['FLUX', KE], ['DIV', KC], ['PSD', (K + 1) * C], ['EXL', KC], ['EXM', KC], ['DEX', KC], ['THL', KC], ['QL', KC], ['QCL', KC], ['THV', KC], ['GEO', KC], ['PIV', V], ['QV', KV], ['QE', KE], ['PHI', KC], ['DRAG', C], ['WIND', C], ['LAPA', KE], ['LAPB', KE], ['DIVS', KC], ['CURLS', KV], ['LAP1', 3 * KC], ['LNPI', C], ['DISS', KE]]);
  const PH = seq([['SFLUX', C], ['OFLUX', C], ['CAP', C], ['ADIF', C], ['MIX', KC], ['DEPTH', C], ['RAIN', C], ['ABS', C], ['OLR', C], ['SH', C], ['EVAP', C], ['INS', C], ['REFL', C], ['TAU', C], ['CONV', C], ['COND', C], ['SWDN', C], ['LAND', C], ['DRAG', C], ['SOIL', C], ['SNOW', C], ['CONC', C], ['RUNOFF', C], ['VEG', C], ['SURF', C], ['DECK', C], ['DECKF', C], ['MLMSUB', C], ['MLMCOVER', C], ['MLMWATER', C], ['MLMENT', C], ['MLMH', C], ['MLMGATE', C], ['MLMTOP', C], ['ATMSW', C], ['CONVMEAN', C], ['CONDMEAN', C], ['STEPRAIN', C], ['ENTRAIN', C], ['BUOY', C], ['USTAR', C], ['STRAT', C], ['REGIME', C], ['MIXTOP', C], ['VRAD', C], ['CTCOOL', C], ['LWH', KC], ['CUMF', C], ['CUTOP', C], ['CUCOVER', cumulusLayers * C], ['CUWATER', cumulusLayers * C], ['MOMU', (momentumLayers + 1) * C], ['MOMK', momentumLayers * C], ['MOMD', (momentumLayers + 1) * C], ['MOMKD', momentumLayers * C], ['MOMS', C], ['ABSSUM', C], ['ATMSUM', C], ['OLRSUM', C], ['INSSUM', C], ['REFLSUM', C], ['ASRMEAN', C], ['OLRMEAN', C], ['ALBMEAN', C], ['ABSCLRSUM', C], ['OLRCLRSUM', C], ['SWCREMEAN', C], ['LWCREMEAN', C], ['LWSFCSUM', C]]);
  const FR = seq([['T', C], ['Z', C], ['RH', C], ['SPD', C], ['WIND', 3 * C], ['DP', C], ['WB', C], ['MI', C], ['W', C], ['WM', C], ['TPW', C], ['TCW', C], ['MSLP', C], ['RAIN', C], ['RUNOFF', C], ['RDONE', C], ['PART', REDUCED.length * groupsOf(C)]]);
  return { C, E, V, K, KC, KE, KV, MI, MF, LV, S, D, PH, FR };
}

function prelude(L, constants) {
  const { C, E, V, K } = L;
  const consts = Object.entries({ ...L.MI, ...Object.fromEntries(Object.entries(L.MF).map(([k, v]) => ['F_' + k, v])), ...Object.fromEntries(Object.entries(L.LV).map(([k, v]) => ['L_' + k, v])), ...Object.fromEntries(Object.entries(L.S).map(([k, v]) => ['S_' + k, v])), ...Object.fromEntries(Object.entries(L.D).map(([k, v]) => ['D_' + k, v])), ...Object.fromEntries(Object.entries(L.PH).map(([k, v]) => ['PH_' + k, v])), ...Object.fromEntries(Object.entries(L.FR).map(([k, v]) => ['FR_' + k, v])), GROUPS: groupsOf(C) })
    .filter(([k]) => k !== 'total').map(([k, v]) => `const ${k}: i32 = ${v};`).join('\n');
  return `
const C: i32 = ${C}; const E: i32 = ${E}; const V: i32 = ${V}; const K: i32 = ${K};
const MAXE: i32 = ${MAX_EDGES}; const MAXEE: i32 = ${MAX_EDGES_ON_EDGE};
const KAPPA: f32 = ${constants.kappa}; const CP: f32 = ${constants.cp}; const P0: f32 = ${constants.p0}; const GRAV: f32 = ${constants.g}; const RGAS: f32 = ${constants.R};
const VIRT: f32 = ${VIRTUAL_FACTOR}; const CDRAG: f32 = ${constants.dragCoefficient}; const GUST: f32 = ${constants.gustiness};
${consts}
@group(0) @binding(0) var<storage, read_write> MI: array<i32>;
@group(0) @binding(1) var<storage, read_write> MF: array<f32>;
@group(0) @binding(2) var<storage, read_write> LV: array<f32>;
@group(0) @binding(3) var<storage, read_write> IN: array<f32>;
@group(0) @binding(4) var<storage, read_write> OUT: array<f32>;
@group(0) @binding(5) var<storage, read_write> D: array<f32>;
@group(0) @binding(6) var<storage, read_write> P: array<f32>;
@group(0) @binding(7) var<storage, read_write> PH: array<f32>;
${constants.physics}
fn diagnoseColumn(i: i32) {
  let pi = IN[S_PI + i];
  let exner0 = pow(pi / P0, KAPPA);
  var thBelow = 0.0; var qBelow = 0.0; var qcBelow = 0.0; var thvBelow = 0.0; var geoBelow = 0.0;
  for (var k = K - 1; k >= 0; k--) {
    let idx = k * C + i;
    let th = IN[S_TH + idx]; let q = IN[S_Q + idx]; let qc = IN[S_QC + idx];
    D[D_EXL + idx] = exner0 * LV[L_CL + k];
    D[D_EXM + idx] = exner0 * LV[L_CM + k];
    D[D_DEX + idx] = exner0 / pi * LV[L_CD + k];
    let thv = th * (1.0 + VIRT * q - qc);
    D[D_THV + idx] = thv;
    var geo = 0.0;
    if (k == K - 1) { geo = CP * exner0 * thv * LV[L_CB + K - 1] - LV[L_GR + K - 1]; }
    else {
      let t = LV[L_CT + k];
      D[D_THL + idx] = th + t * (thBelow - th);
      D[D_QL + idx] = q + t * (qBelow - q);
      D[D_QCL + idx] = qc + t * (qcBelow - qc);
      geo = geoBelow + CP * exner0 * (thvBelow * LV[L_CA + k] + thv * LV[L_CB + k]) - LV[L_GR + k];
    }
    D[D_GEO + idx] = geo;
    thBelow = th; qBelow = q; qcBelow = qc; thvBelow = thv; geoBelow = geo;
  }
}
${PHYSICS_FUNCTIONS}
`;
}

/*
 * The global sums behind the model's diagnostics, one partial per
 * workgroup: [name, how the partials combine, the cell's term].
 */
const REDUCED = [
  ['area', 'sum', 'a'], ['mass', 'sum', 'a * pi'], ['surfaceT', 'sum', 'a * IN[S_TS + i]'],
  ['piMin', 'min', 'pi'], ['piMax', 'max', 'pi'], ['maxWind', 'max', 'wind'],
  ['water', 'sum', 'a * water'], ['cloud', 'sum', 'a * cloud'],
  ['rain', 'sum', 'a * PH[PH_RAIN + i]'], ['recentRain', 'sum', 'a * OUT[FR_RAIN + i]'],
  ['iceArea', 'sum', 'a * cover'], ['iceVolume', 'sum', 'a * cover * max(0.0, IN[S_ICE + i])'], ['albedo', 'sum', 'a * PH[PH_ADIF + i]'],
  ['absorbedSolar', 'sum', 'a * PH[PH_ABS + i]'], ['atmosphereSolar', 'sum', 'a * PH[PH_ATMSW + i]'], ['outgoingLongwave', 'sum', 'a * PH[PH_OLR + i]'], ['sensibleHeat', 'sum', 'a * PH[PH_SH + i]'],
  ['evaporation', 'sum', 'a * PH[PH_EVAP + i]'], ['insolation', 'sum', 'a * PH[PH_INS + i]'], ['reflectedSolar', 'sum', 'a * PH[PH_REFL + i]'],
  ['absorbedSum', 'sum', 'a * PH[PH_ABSSUM + i]'], ['atmosphereSum', 'sum', 'a * PH[PH_ATMSUM + i]'], ['outgoingSum', 'sum', 'a * PH[PH_OLRSUM + i]'],
  ['insolationSum', 'sum', 'a * PH[PH_INSSUM + i]'], ['reflectedSum', 'sum', 'a * PH[PH_REFLSUM + i]'],
  ['clearAbsorbedSum', 'sum', 'a * PH[PH_ABSCLRSUM + i]'], ['clearOutgoingSum', 'sum', 'a * PH[PH_OLRCLRSUM + i]'],
  ['landArea', 'sum', 'land * a'], ['landT', 'sum', 'land * a * IN[S_TS + i]'], ['snowArea', 'sum', 'select(0.0, land * a, PH[PH_SNOW + i] > 1.0)'],
  ['soil', 'sum', 'land * a * PH[PH_SOIL + i]'],
];
const REDUCED_SETUP = `    let a = MF[F_AREA + i]; let pi = IN[S_PI + i];
    let conc = PH[PH_CONC + i]; let cover = select(0.0, select(conc, 1.0, conc <= 0.0), IN[S_ICE + i] > 0.0);
    let land = select(0.0, 1.0, PH[PH_LAND + i] > 0.5);
    var water = 0.0; var cloud = 0.0; var wind = 0.0;
    for (var k = 0; k < K; k++) {
      let d = pi * LV[L_DS + k] / GRAV;
      water += d * IN[S_Q + k * C + i]; cloud += d * IN[S_QC + k * C + i];
      for (var m = 0; m < MAXE; m++) { wind = max(wind, abs(IN[S_U + k * E + MI[EOC + MAXE * i + m]])); }
    }`;

/*
 * What the page draws, computed where the state lives. frameFields
 * interpolates in ln p to the pressure P[0] (0 for the lowest layer) as
 * levels.module.js does, with the column's geopotential integrated in
 * registers, derives the comfort measures there as levels.module.js
 * does, adds the column water, cloud and sea-level pressure, and when
 * P[3] asks for it the vertical velocity at the level as
 * levels.module.js derives it from the layers' mass-flux divergences;
 * frameVertical then averages it with the neighbours' and folds it into
 * its memory, WM ← WM·P[4] + smoothed·(1 − P[4]).
 * frameRain folds the step accumulators into the three-hour rain,
 * S ← S·P[1] + rain, and the running runoff, which moves to RDONE for
 * the host to count when the diagnostics are taken (P[2] = 1); when
 * P[5] is 86400 over the seconds since the last frame, it turns the
 * convective and large-scale sums into their means in mm/d, CONVMEAN and
 * CONDMEAN, and every frame clears the sums with the rain. Alike, when
 * P[6] is one over the steps since the last frame, the sums of the
 * radiation over those steps become the means ASRMEAN and OLRMEAN and
 * the albedo ALBMEAN (reflected over incoming summed, 0 where no sun
 * rose), with clearSkyPass the cloud effects SWCREMEAN (absorbed less
 * clear-sky absorbed) and LWCREMEAN (clear-sky less all-sky outgoing),
 * and every frame clears those sums too.
 */
const COMFORT_WGSL = `
fn dewPointC(t: f32, rh: f32) -> f32 {
  let g = log(clamp(rh, 1e-3, 1.0)) + 17.625 * t / (243.04 + t);
  return 243.04 * g / (17.625 - g);
}
fn wetBulbC(t: f32, rh: f32) -> f32 {
  let p = 100.0 * clamp(rh, 0.05, 0.99);
  return t * atan(0.151977 * sqrt(p + 8.313659)) + atan(t + p) - atan(p - 1.676331) + 0.00391838 * pow(p, 1.5) * atan(0.023101 * p) - 4.686035;
}
fn heatIndexC(t: f32, rh: f32) -> f32 {
  let f = t * 9.0 / 5.0 + 32.0; let r = 100.0 * clamp(rh, 0.0, 1.0);
  var hi = -42.379 + 2.04901523 * f + 10.14333127 * r - 0.22475541 * f * r - 6.83783e-3 * f * f - 5.481717e-2 * r * r + 1.22874e-3 * f * f * r + 8.5282e-4 * f * r * r - 1.99e-6 * f * f * r * r;
  if (r < 13.0 && f <= 112.0) { hi -= ((13.0 - r) / 4.0) * sqrt(max(0.0, (17.0 - abs(f - 95.0)) / 17.0)); }
  else if (r > 85.0 && f <= 87.0) { hi += ((r - 85.0) / 10.0) * ((87.0 - f) / 5.0); }
  return (hi - 32.0) * 5.0 / 9.0;
}
fn windChillC(t: f32, v: f32) -> f32 {
  let k = 3.6 * v;
  if (k < 4.8) { return t; }
  let p = pow(k, 0.16);
  return 13.12 + 0.6215 * t - 11.37 * p + 0.3965 * t * p;
}
fn miseryC(t: f32, rh: f32, v: f32) -> f32 {
  if (t >= 26.7) { return max(t, heatIndexC(t, rh)); }
  if (t <= 10.0) { return min(t, windChillC(t, v)); }
  return t;
}
`;

const FRAME_KERNELS = {
  frameFields: `${COMFORT_WGSL}
@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= C) { return; }
  let pi = IN[S_PI + i];
  let exner0 = pow(pi / P0, KAPPA);
  let pressure = P[0];
  var k = K - 2;
  if (pressure > 0.0) { k = 0; while (k < K - 2 && pi * LV[L_SM + k + 1] < pressure) { k++; } }
  var geo = 0.0; var thvBelow = 0.0; var geoK = 0.0; var geoK1 = 0.0;
  for (var j = K - 1; j >= k; j--) {
    let n = j * C + i;
    let thv = IN[S_TH + n] * (1.0 + VIRT * IN[S_Q + n] - IN[S_QC + n]);
    if (j == K - 1) { geo = CP * exner0 * thv * LV[L_CB + j] - LV[L_GR + j]; }
    else { geo = geo + CP * exner0 * (thvBelow * LV[L_CA + j] + thv * LV[L_CB + j]) - LV[L_GR + j]; }
    if (j == k + 1) { geoK1 = geo + LV[L_GABS + j]; }
    if (j == k) { geoK = geo + LV[L_GABS + j]; }
    thvBelow = thv;
  }
  let pk = pi * LV[L_SM + k]; let pk1 = pi * LV[L_SM + k + 1];
  var t = 1.0;
  if (pressure > 0.0) { t = (log(pressure) - log(pk)) / (log(pk1) - log(pk)); }
  let tw = clamp(t, 0.0, 1.0);
  let wind = cellWind(i, k) + tw * (cellWind(i, k + 1) - cellWind(i, k));
  let n = k * C + i; let n1 = n + C;
  let tk = IN[S_TH + n] * exner0 * LV[L_CM + k]; let tk1 = IN[S_TH + n1] * exner0 * LV[L_CM + k + 1];
  let temperature = tk + tw * (tk1 - tk);
  let lowest = pi * LV[L_SM + K - 1];
  var here = lowest;
  if (pressure > 0.0) { here = min(pressure, lowest); }
  let q = IN[S_Q + n] + tw * (IN[S_Q + n1] - IN[S_Q + n]);
  let phis = MF[F_PHIS + i];
  var height = (geoK + t * (geoK1 - geoK) + phis) / GRAV;
  if (t > 1.0) { height = (geoK1 + phis - RGAS * tk1 * log(pressure / lowest)) / GRAV; }
  let rh = min(1.5, q / qsat(temperature, here)); let speed = length(wind); let celsius = temperature - 273.15;
  OUT[FR_T + i] = temperature; OUT[FR_Z + i] = height; OUT[FR_RH + i] = rh;
  OUT[FR_SPD + i] = speed;
  OUT[FR_DP + i] = dewPointC(celsius, rh) + 273.15; OUT[FR_WB + i] = wetBulbC(celsius, rh) + 273.15; OUT[FR_MI + i] = miseryC(celsius, rh, speed) + 273.15;
  OUT[FR_WIND + 3 * i] = wind.x; OUT[FR_WIND + 3 * i + 1] = wind.y; OUT[FR_WIND + 3 * i + 2] = wind.z;
  var water = 0.0; var cloud = 0.0;
  for (var j = 0; j < K; j++) {
    let d = pi * LV[L_DS + j] / GRAV;
    water += d * IN[S_Q + j * C + i]; cloud += d * IN[S_QC + j * C + i];
  }
  OUT[FR_TPW + i] = water; OUT[FR_TCW + i] = cloud + PH[PH_DECKF + i] * PH[PH_DECK + i];
  let tb = IN[S_TH + (K - 1) * C + i] * exner0 * LV[L_CM + K - 1];
  OUT[FR_MSLP + i] = pi * exp(phis / (RGAS * (tb + 0.00325 * phis / GRAV)));
  if (P[3] > 0.5) {
    var massDivergence: array<f32, K>; var divergence: array<f32, K>;
    for (var j = 0; j < K; j++) { massDivergence[j] = 0.0; divergence[j] = 0.0; }
    for (var m = 0; m < MAXE; m++) {
      let e = MI[EOC + MAXE * i + m]; let f = f32(MI[ESC + MAXE * i + m]) * MF[F_DV + e];
      let piEdge = 0.5 * (IN[S_PI + MI[COE + 2 * e]] + IN[S_PI + MI[COE + 2 * e + 1]]);
      for (var j = 0; j < K; j++) { let flow = f * IN[S_U + j * E + e]; massDivergence[j] += piEdge * flow; divergence[j] += flow; }
    }
    var dPi = 0.0;
    for (var j = 0; j < K; j++) { massDivergence[j] /= MF[F_AREA + i]; divergence[j] /= MF[F_AREA + i]; dPi -= massDivergence[j] * LV[L_DS + j]; }
    let sigma = here / pi;
    let advectionK = massDivergence[k] - pi * divergence[k]; let advectionK1 = massDivergence[k + 1] - pi * divergence[k + 1];
    let piAdvection = advectionK + tw * (advectionK1 - advectionK);
    var layer = 0;
    while (layer < K - 1 && LV[L_SL + layer] < sigma) { layer++; }
    let s = clamp((sigma - LV[L_SU + layer]) / LV[L_DS + layer], 0.0, 1.0);
    var cumulative = 0.0; var upper = 0.0; var flow = 0.0;
    for (var j = 0; j < K; j++) {
      cumulative += massDivergence[j] * LV[L_DS + j];
      var lower = -cumulative - LV[L_SL + j] * dPi;
      if (j == K - 1) { lower = 0.0; }
      if (j == layer) { flow = upper + s * (lower - upper); }
      upper = lower;
    }
    let omega = flow + sigma * (dPi + piAdvection);
    OUT[FR_W + i] = -omega * RGAS * temperature / (here * GRAV);
  }
}`,
  frameVertical: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= C) { return; }
  var sum = OUT[FR_W + i]; let n = MI[NEC + i];
  for (var m = 0; m < MAXE; m++) { if (m < n) { sum += OUT[FR_W + MI[COC + MAXE * i + m]]; } }
  OUT[FR_WM + i] = OUT[FR_WM + i] * P[4] + (sum / f32(n + 1)) * (1.0 - P[4]);
}`,
  frameRain: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= C) { return; }
  OUT[FR_RAIN + i] = OUT[FR_RAIN + i] * P[1] + PH[PH_RAIN + i];
  let tally = OUT[FR_RUNOFF + i] + PH[PH_RUNOFF + i];
  if (P[2] > 0.5) { OUT[FR_RDONE + i] = tally; OUT[FR_RUNOFF + i] = 0.0; } else { OUT[FR_RUNOFF + i] = tally; }
  if (P[5] > 0.0) { PH[PH_CONVMEAN + i] = PH[PH_CONV + i] * P[5]; PH[PH_CONDMEAN + i] = PH[PH_COND + i] * P[5]; }
  PH[PH_RAIN + i] = 0.0; PH[PH_RUNOFF + i] = 0.0; PH[PH_CONV + i] = 0.0; PH[PH_COND + i] = 0.0;
  if (P[6] > 0.0) {
    PH[PH_ASRMEAN + i] = PH[PH_ABSSUM + i] * P[6]; PH[PH_OLRMEAN + i] = PH[PH_OLRSUM + i] * P[6];
    PH[PH_ALBMEAN + i] = select(0.0, PH[PH_REFLSUM + i] / PH[PH_INSSUM + i], PH[PH_INSSUM + i] > 0.0);
    if (CLEAR_SKY) { PH[PH_SWCREMEAN + i] = (PH[PH_ABSSUM + i] - PH[PH_ABSCLRSUM + i]) * P[6]; PH[PH_LWCREMEAN + i] = (PH[PH_OLRCLRSUM + i] - PH[PH_OLRSUM + i]) * P[6]; }
  }
  PH[PH_ABSSUM + i] = 0.0; PH[PH_ATMSUM + i] = 0.0; PH[PH_OLRSUM + i] = 0.0; PH[PH_INSSUM + i] = 0.0; PH[PH_REFLSUM + i] = 0.0; PH[PH_ABSCLRSUM + i] = 0.0; PH[PH_OLRCLRSUM + i] = 0.0; PH[PH_LWSFCSUM + i] = 0.0;
}`,
};

const KERNELS = {
  flux: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let n = i32(id.x); if (n >= K * E) { return; }
  let e = n % E;
  let piEdge = 0.5 * (IN[S_PI + MI[COE + 2 * e]] + IN[S_PI + MI[COE + 2 * e + 1]]);
  D[D_FLUX + n] = piEdge * IN[S_U + n];
}`,
  divergence: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let n = i32(id.x); if (n >= K * C) { return; }
  let k = n / C; let i = n % C;
  var sum = 0.0;
  for (var m = 0; m < MAXE; m++) {
    let e = MI[EOC + MAXE * i + m];
    sum += f32(MI[ESC + MAXE * i + m]) * D[D_FLUX + k * E + e] * MF[F_DV + e];
  }
  D[D_DIV + n] = sum / MF[F_AREA + i];
}`,
  column: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= C) { return; }
  let pi = IN[S_PI + i];
  var sum = 0.0;
  for (var k = 0; k < K; k++) { sum += D[D_DIV + k * C + i] * LV[L_DS + k]; }
  let dPi = -sum;
  OUT[S_PI + i] = dPi;
  var cumulative = 0.0;
  D[D_PSD + i] = 0.0;
  for (var k = 0; k < K; k++) {
    cumulative += D[D_DIV + k * C + i] * LV[L_DS + k];
    D[D_PSD + (k + 1) * C + i] = -cumulative - LV[L_SL + k] * dPi;
  }
  D[D_PSD + K * C + i] = 0.0;
  D[D_LNPI + i] = log(pi);
  diagnoseColumn(i);
  let bottom = (K - 1) * C + i;
  let airT = IN[S_TH + bottom] * D[D_EXM + bottom];
  let rho = pi * LV[L_SM + K - 1] / (RGAS * airT);
  let mass = pi * LV[L_DS + K - 1] / GRAV;
  D[D_DRAG + i] = PH[PH_DRAG + i] * rho * max(D[D_WIND + i], GUST) / mass;
}`,
  vertexPi: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let v = i32(id.x); if (v >= V) { return; }
  var sum = 0.0;
  for (var m = 0; m < 3; m++) { sum += MF[F_KAV + 3 * v + m] * IN[S_PI + MI[COV + 3 * v + m]]; }
  D[D_PIV + v] = sum / MF[F_ATRI + v];
}`,
  pvVertex: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let n = i32(id.x); if (n >= K * V) { return; }
  let k = n / V; let v = n % V;
  var zeta = 0.0;
  for (var m = 0; m < 3; m++) {
    let e = MI[EOV + 3 * v + m];
    zeta += f32(MI[ESV + 3 * v + m]) * IN[S_U + k * E + e] * MF[F_DC + e];
  }
  zeta = zeta / MF[F_ATRI + v];
  D[D_QV + n] = (zeta + MF[F_FV + v]) / D[D_PIV + v];
}`,
  pvEdge: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let n = i32(id.x); if (n >= K * E) { return; }
  let k = n / E; let e = n % E;
  D[D_QE + n] = 0.5 * (D[D_QV + k * V + MI[VOE + 2 * e]] + D[D_QV + k * V + MI[VOE + 2 * e + 1]]);
}`,
  cellTendency: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= C) { return; }
  let pi = IN[S_PI + i]; let area = MF[F_AREA + i]; let dragHere = D[D_DRAG + i];
  var edges: array<i32, MAXE>; var neighbours: array<i32, MAXE>; var signs: array<f32, MAXE>; var dcs: array<f32, MAXE>; var dvs: array<f32, MAXE>; var dragThere: array<f32, MAXE>;
  for (var m = 0; m < MAXE; m++) {
    let slot = MAXE * i + m;
    let e = MI[EOC + slot]; let j = MI[COC + slot]; let sign = f32(MI[ESC + slot]);
    edges[m] = e; neighbours[m] = j; signs[m] = sign;
    dcs[m] = MF[F_DC + e]; dvs[m] = MF[F_DV + e] * abs(sign); dragThere[m] = D[D_DRAG + j];
  }
  var upperFlow = D[D_PSD + i];
  for (var k = 0; k < K; k++) {
    let idx = k * C + i; let row = k * C; let flux0 = k * E;
    let fT = IN[S_TH + idx]; let fQ = IN[S_Q + idx]; let fC = IN[S_QC + idx];
    let top = LV[L_TOP + k]; let bottomLayer = k == K - 1;
    var kinetic = 0.0; var dragPower = 0.0; var divT = 0.0; var divQ = 0.0; var divC = 0.0;
    for (var m = 0; m < MAXE; m++) {
      let e = edges[m]; let j = neighbours[m];
      let dc = dcs[m]; let dv = dvs[m];
      let u = IN[S_U + flux0 + e];
      kinetic += 0.25 * dc * dv * u * u;
      let rate = top + select(0.0, 0.5 * (dragHere + dragThere[m]), bottomLayer);
      dragPower += 0.5 * dc * dv * rate * u * u;
      let carried = signs[m] * D[D_FLUX + flux0 + e] * 0.5;
      divT += carried * (fT + IN[S_TH + row + j]) * dv;
      divQ += carried * (fQ + IN[S_Q + row + j]) * dv;
      divC += carried * (fC + IN[S_QC + row + j]) * dv;
    }
    D[D_PHI + idx] = D[D_GEO + idx] + kinetic / area;
    let divergence = D[D_DIV + idx];
    let lowerFlow = D[D_PSD + (k + 1) * C + i];
    let layerMass = pi * LV[L_DS + k];
    var lowerT = 0.0; var lowerQ = 0.0; var lowerC = 0.0; var upperT = 0.0; var upperQ = 0.0; var upperC = 0.0;
    if (k < K - 1) { lowerT = D[D_THL + idx]; lowerQ = D[D_QL + idx]; lowerC = D[D_QCL + idx]; }
    if (k > 0) { upperT = D[D_THL + idx - C]; upperQ = D[D_QL + idx - C]; upperC = D[D_QCL + idx - C]; }
    let verticalT = (lowerFlow * lowerT - upperFlow * upperT - fT * (lowerFlow - upperFlow)) / layerMass;
    let verticalQ = (lowerFlow * lowerQ - upperFlow * upperQ - fQ * (lowerFlow - upperFlow)) / layerMass;
    let verticalC = (lowerFlow * lowerC - upperFlow * upperC - fC * (lowerFlow - upperFlow)) / layerMass;
    OUT[S_TH + idx] = -(divT / area - fT * divergence) / pi - verticalT + dragPower / area / (CP * D[D_EXM + idx]);
    OUT[S_Q + idx] = -(divQ / area - fQ * divergence) / pi - verticalQ;
    OUT[S_QC + idx] = -(divC / area - fC * divergence) / pi - verticalC;
    upperFlow = lowerFlow;
  }
}`,
  momentum: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let e = i32(id.x); if (e >= E) { return; }
  let i = MI[COE + 2 * e]; let j = MI[COE + 2 * e + 1];
  let dc = MF[F_DC + e];
  var others: array<i32, MAXEE>; var weights: array<f32, MAXEE>;
  for (var s = 0; s < MAXEE; s++) { others[s] = MI[EOE + MAXEE * e + s]; weights[s] = MF[F_PVW + MAXEE * e + s]; }
  let piEdge = 0.5 * (IN[S_PI + i] + IN[S_PI + j]);
  let lnPiStep = D[D_LNPI + j] - D[D_LNPI + i];
  let gphis = MF[F_GPHIS + e];
  var upperFlow = 0.5 * (D[D_PSD + i] + D[D_PSD + j]);
  var uAbove = 0.0;
  var u = IN[S_U + e];
  for (var k = 0; k < K; k++) {
    let n = k * E + e; let flux0 = k * E; let off = k * C;
    let qHere = 0.5 * D[D_QE + n];
    var pv = 0.0;
    for (var s = 0; s < MAXEE; s++) {
      let other = others[s];
      pv += weights[s] * D[D_FLUX + flux0 + other] * (qHere + 0.5 * D[D_QE + flux0 + other]);
    }
    let gradPhi = (D[D_PHI + off + j] - D[D_PHI + off + i]) / dc;
    let pgfPi = RGAS * 0.5 * (D[D_THV + off + i] * D[D_EXM + off + i] + D[D_THV + off + j] * D[D_EXM + off + j]) * lnPiStep / dc;
    let lowerFlow = 0.5 * (D[D_PSD + (k + 1) * C + i] + D[D_PSD + (k + 1) * C + j]);
    var uBelow = 0.0; var lowerU = 0.0; var upperU = 0.0;
    if (k < K - 1) { uBelow = IN[S_U + n + E]; lowerU = 0.5 * (u + uBelow); }
    if (k > 0) { upperU = 0.5 * (u + uAbove); }
    let vertical = (lowerFlow * lowerU - upperFlow * upperU - u * (lowerFlow - upperFlow)) / (piEdge * LV[L_DS + k]);
    var du = pv / dc - gradPhi - gphis - pgfPi - vertical;
    if (k == K - 1) { du -= 0.5 * (D[D_DRAG + i] + D[D_DRAG + j]) * u; }
    du -= LV[L_TOP + k] * u;
    OUT[S_U + n] = du;
    uAbove = u; u = uBelow; upperFlow = lowerFlow;
  }
}`,
  advance: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let n = i32(id.x); if (n >= ${'S_TOTAL'}) { return; }
  OUT[n] = IN[n] + P[0] * D[n];
}`,
  combine: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let n = i32(id.x); if (n >= ${'S_TOTAL'}) { return; }
  IN[n] += P[0] * (OUT[n] + 2.0 * D[n] + 2.0 * MF[n] + LV[n]);
}`,
  lapScalar1: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= C) { return; }
  let area = MF[F_AREA + i]; let piHere = IN[S_PI + i];
  var neighbours: array<i32, MAXE>; var dvs: array<f32, MAXE>; var dcs: array<f32, MAXE>; var piThere: array<f32, MAXE>;
  for (var m = 0; m < MAXE; m++) {
    let e = MI[EOC + MAXE * i + m]; let j = MI[COC + MAXE * i + m];
    neighbours[m] = j; dvs[m] = MF[F_DV + e]; dcs[m] = MF[F_DC + e]; piThere[m] = IN[S_PI + j];
  }
  for (var k = 0; k < K; k++) {
    let idx = k * C + i; let row = k * C;
    let hereT = IN[S_TH + idx] * D[D_EXM + idx]; let hereQ = IN[S_Q + idx] * piHere; let hereC = IN[S_QC + idx] * piHere;
    var sumT = 0.0; var sumQ = 0.0; var sumC = 0.0;
    for (var m = 0; m < MAXE; m++) {
      let j = neighbours[m]; let dv = dvs[m]; let dc = dcs[m];
      sumT += dv * (IN[S_TH + row + j] * D[D_EXM + row + j] - hereT) / dc;
      sumQ += dv * (IN[S_Q + row + j] * piThere[m] - hereQ) / dc;
      sumC += dv * (IN[S_QC + row + j] * piThere[m] - hereC) / dc;
    }
    D[D_LAP1 + idx] = sumT / area; D[D_LAP1 + K * C + idx] = sumQ / area; D[D_LAP1 + 2 * K * C + idx] = sumC / area;
  }
}`,
  lapScalar2: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= C) { return; }
  let area = MF[F_AREA + i]; let piHere = IN[S_PI + i];
  var neighbours: array<i32, MAXE>; var dvs: array<f32, MAXE>; var dcs: array<f32, MAXE>;
  for (var m = 0; m < MAXE; m++) {
    let e = MI[EOC + MAXE * i + m];
    neighbours[m] = MI[COC + MAXE * i + m]; dvs[m] = MF[F_DV + e]; dcs[m] = MF[F_DC + e];
  }
  for (var k = 0; k < K; k++) {
    let idx = k * C + i; let row = k * C;
    for (var field = 0; field < 3; field++) {
      let base = field * K * C;
      let here = D[D_LAP1 + base + idx];
      var sum = 0.0;
      for (var m = 0; m < MAXE; m++) { sum += dvs[m] * (D[D_LAP1 + base + row + neighbours[m]] - here) / dcs[m]; }
      let lap2 = sum / area;
      if (field == 0) { IN[S_TH + idx] -= P[0] * lap2 / D[D_EXM + idx]; }
      else if (field == 1) { IN[S_Q + idx] -= P[0] * lap2 / piHere; }
      else { IN[S_QC + idx] -= P[0] * lap2 / piHere; }
    }
  }
}`,
  divCurl: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let n = i32(id.x);
  let fromLap = P[1] > 0.5;
  if (n < C) {
    let i = n;
    let area = MF[F_AREA + i];
    var edges: array<i32, MAXE>; var weights: array<f32, MAXE>;
    for (var m = 0; m < MAXE; m++) { let e = MI[EOC + MAXE * i + m]; edges[m] = e; weights[m] = f32(MI[ESC + MAXE * i + m]); }
    for (var k = 0; k < K; k++) {
      var sum = 0.0;
      for (var m = 0; m < MAXE; m++) {
        let e = edges[m];
        let u = select(IN[S_U + k * E + e], D[D_LAPA + k * E + e], fromLap);
        sum += weights[m] * u * MF[F_DV + e];
      }
      D[D_DIVS + k * C + i] = sum / area;
    }
  } else if (n < C + V) {
    let v = n - C;
    let area = MF[F_ATRI + v];
    var edges: array<i32, 3>; var weights: array<f32, 3>;
    for (var m = 0; m < 3; m++) { let e = MI[EOV + 3 * v + m]; edges[m] = e; weights[m] = f32(MI[ESV + 3 * v + m]); }
    for (var k = 0; k < K; k++) {
      var sum = 0.0;
      for (var m = 0; m < 3; m++) {
        let e = edges[m];
        let u = select(IN[S_U + k * E + e], D[D_LAPA + k * E + e], fromLap);
        sum += weights[m] * u * MF[F_DC + e];
      }
      D[D_CURLS + k * V + v] = sum / area;
    }
  }
}`,
  lapVelocity: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let e = i32(id.x); if (e >= E) { return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1]; let va = MI[VOE + 2 * e]; let vb = MI[VOE + 2 * e + 1];
  let dc = MF[F_DC + e]; let dv = MF[F_DV + e];
  let second = P[1] > 0.5;
  for (var k = 0; k < K; k++) {
    let n = k * E + e;
    let lap = (D[D_DIVS + k * C + b] - D[D_DIVS + k * C + a]) / dc
      - (D[D_CURLS + k * V + vb] - D[D_CURLS + k * V + va]) / dv;
    if (second) {
      let before = IN[S_U + n]; let after = before - P[0] * lap;
      IN[S_U + n] = after;
      D[D_DISS + n] = before * before - after * after;
    } else { D[D_LAPA + n] = lap; }
  }
}`,
  divergenceDamp: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let e = i32(id.x); if (e >= E) { return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  let dc = MF[F_DC + e];
  let accumulate = P[2] > 0.5;
  for (var k = 0; k < K; k++) {
    let n = k * E + e;
    let before = IN[S_U + n]; let after = before + P[0] * (D[D_DIVS + k * C + b] - D[D_DIVS + k * C + a]) / dc;
    IN[S_U + n] = after;
    let lost = before * before - after * after;
    if (accumulate) { D[D_DISS + n] += lost; } else { D[D_DISS + n] = lost; }
  }
}`,
  dissipationHeat: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let n = i32(id.x); if (n >= K * C) { return; }
  let k = n / C; let i = n % C;
  var sum = 0.0;
  for (var m = 0; m < MAXE; m++) { let e = MI[EOC + MAXE * i + m]; sum += abs(f32(MI[ESC + MAXE * i + m])) * MF[F_DC + e] * MF[F_DV + e] * D[D_DISS + k * E + e]; }
  IN[S_TH + n] += 0.25 * sum / MF[F_AREA + i] / (CP * D[D_EXM + n]);
}`,
  dissipationClear: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let n = i32(id.x); if (n >= K * E) { return; }
  D[D_DISS + n] = 0.0;
}`,
};

export const PHYSICS_DEFAULTS = {
  solarConstant: 1362, cloudAbsorption: 130, cloudScattering: 95, cloudSolarAbsorption: 0.4, stratus: true, stratusIndex: 'eis', stratusScale: 0.15, stratusWaterMax: 0.15, stratusSigma: 0.92,
  mixedLayerDeck: true, mixedLayer: {}, stratusSubsidence: -1e-3, minimumInversion: 4, ceilingInversion: null, subsidenceMemory: 2 * 86400, subsidenceSmoothing: 2, cloudCover: 'pdf', criticalHumidity: 0.8, boundaryCriticalHumidity: 0.85, coverFloor: 0.01, overcastWater: 5e-5, overcastInversion: [8, 12], cloudOverlap: 'maximumRandom', prognosticHeight: true, deckRest: 'regime', cumulusCeiling: 2000, gateMemory: 86400, stratusSolar: true, window: 0.25, tauEquator: 5.3, tauPole: 1.325, linearFraction: 0.1,
  gasFraction: 0.2, gasOpticalDepth: 7, ozoneAbsorption: 0.03, vaporAbsorption: 1, ozoneHeight: 25e3, ozoneWidth: 5e3, ozoneOpacity: 4, scaleHeight: 7e3,
  exchangeCoefficient: SEA_DRAG, latentHeat: 2.5e6, vaporCoupling: 0.55, skylight: 0, clearSkyPass: false,
  longwaveScheme: 'correlated', solarGases: 'clirad', ...GREENHOUSE_GASES, ozone: 'afgl', ozoneColumn: OZONE_COLUMN, ozoneProfile: null, vaporStrength: VAPOR_STRENGTH,
  rayleighBands: RAYLEIGH_BANDS, rayleighDepth: null, upwardAbsorption: true, visibleFraction: 0.5, landAerosol: LAND_AEROSOL, seaAerosol: SEA_AEROSOL, aerosolAlbedo: 0.95, aerosolAsymmetry: 0.7, aerosolHeight: 2000,
  slabHeatCapacity: 2.1e7, skinHeatCapacity: 2e5, conductivity: 2, minimumThickness: 0.1, iceDensity: 917, latentHeatFusion: 3.34e5, leadClosing: 0.3, leadExchange: 10,
  diffuseWaterAlbedo: 0.06, iceAlbedo: 0.5, fullAlbedoThickness: 0.5, iceSnowAlbedo: 0.75, iceFullSnow: 20, snowConductivity: 0.31, snowDensity: 300, waterDensity: 1026,
  ...MOIST_DEFAULTS,
  richardsonCritical: 0.5, vonKarman: 0.4, searchTop: 0.5, stability: true, turbulence: 'moist', cloudTop: {},
  boundaryCover: 'variance', varianceFloor: 0.002, varianceScale: 5, mixingLength: 300, stableMixingLength: 30, deckRegime: 'inversion', deckBypass: false,
  landed: false, landHeatCapacity: 1e6, bucketCapacity: 150, wetnessThreshold: 0.75, landAlbedo: 0.2, snowAlbedo: 0.55, fullSnow: 20,
  vegetation: true, bareAlbedo: 0.30, vegetatedAlbedo: 0.13, soilDarkening: true, wetSoilAlbedo: 0.15, darkeningWetness: [0.2, 0.5], rootZoneCapacity: 300, dryWetness: 0.1, wetWetness: 0.6, iceSheetAlbedo: 0.8,
  growthTime: 180 * 86400, declineTime: 365 * 86400, snowDeclineTime: 720 * 86400, surfaceCapacity: 15, percolationTime: 86400, stomatalResistance: 70, growthColdest: 278.15, growthWarmest: 288.15,
};

export async function createGpuCore(mesh, {
  levels = sigmaInterfaces(), g = GRAVITY, cp = CP_DRY, R = R_DRY, p0 = P0, nu4 = 0, nu4Theta = 0, divergenceDamping = 0,
  dragCoefficient = SEA_DRAG, gustiness = 3, topSigma = 0.02, topDragDays = 5, referenceTheta = null, surfaceGeopotential = null, physics: physicsOptions = {},
} = {}) {
  const phys = { ...PHYSICS_DEFAULTS, ...physicsOptions, R };
  const { device } = await getDevice();
  const K = levels.length - 1;
  let cumulusK0 = K;
  while (cumulusK0 > 0 && 0.5 * (levels[cumulusK0 - 1] + levels[cumulusK0]) * MAXIMUM_SURFACE_PRESSURE > phys.shallowTop) cumulusK0--;
  phys.cumulusK0 = cumulusK0;
  const L = layoutFor(mesh, K, K - cumulusK0, phys.plumeMomentum ? K : 0);
  const { C, E, V } = L;
  const kappa = R / cp;
  const sigmaUpper = levels.subarray(0, K), sigmaLower = levels.subarray(1, K + 1);
  const dSigma = Float64Array.from(sigmaLower, (s, k) => s - sigmaUpper[k]);
  const sigmaMid = Float64Array.from(sigmaLower, (s, k) => 0.5 * (s + sigmaUpper[k]));
  const topRate = Float64Array.from(sigmaMid, (s) => (topDragDays > 0 && s < topSigma ? (topSigma - s) / topSigma / (topDragDays * 86400) : 0));
  const shape = Float64Array.from({ length: K }, (_, k) => phys.linearFraction * (levels[k + 1] - levels[k]) + (1 - phys.linearFraction) * (levels[k + 1] ** 4 - levels[k] ** 4));
  const ozoneAbove = (sigma) => (sigma <= 0 ? 0 : (1 + Math.exp(-phys.ozoneHeight / phys.ozoneWidth)) / (1 + Math.exp((-phys.scaleHeight * Math.log(sigma) - phys.ozoneHeight) / phys.ozoneWidth)));
  const beamLeft = (sigma) => Math.exp(-phys.ozoneOpacity * ozoneAbove(sigma));
  const ozoneFraction = Float64Array.from({ length: K }, (_, k) => (beamLeft(levels[k]) - beamLeft(levels[k + 1])) / (1 - Math.exp(-phys.ozoneOpacity)));
  const gasEmissivity = Float64Array.from({ length: K }, (_, k) => 1 - Math.exp(-phys.gasOpticalDepth * (levels[k + 1] - levels[k])));
  const aerosolFraction = Float64Array.from({ length: K }, (_, k) => levels[k + 1] ** (phys.scaleHeight / phys.aerosolHeight) - levels[k] ** (phys.scaleHeight / phys.aerosolHeight));
  const ozoneShare = Float64Array.from({ length: K }, (_, k) => ozoneAbove(levels[k + 1]) - ozoneAbove(levels[k]));
  let kTop = 0;
  while (kTop < K - 1 && sigmaMid[kTop] <= phys.searchTop) kTop++;

  const mi = new Int32Array(L.MI.total);
  const put = (arr, off, src) => arr.set(src, off);
  const padded = (source, fill) => Int32Array.from({ length: MAX_EDGES * C }, (_, slot) => {
    const i = Math.floor(slot / MAX_EDGES), m = slot % MAX_EDGES;
    return m < mesh.nEdgesOnCell[i] ? source[mesh.maxEdges * i + m] : fill(i);
  });
  put(mi, L.MI.COE, mesh.cellsOnEdge); put(mi, L.MI.VOE, mesh.verticesOnEdge);
  put(mi, L.MI.EOC, padded(mesh.edgesOnCell, (i) => mesh.edgesOnCell[mesh.maxEdges * i])); put(mi, L.MI.ESC, padded(mesh.edgeSignOnCell, () => 0));
  put(mi, L.MI.COC, padded(mesh.cellsOnCell, (i) => i)); put(mi, L.MI.NEC, mesh.nEdgesOnCell); put(mi, L.MI.COV, mesh.cellsOnVertex); put(mi, L.MI.EOV, mesh.edgesOnVertex);
  put(mi, L.MI.ESV, mesh.edgeSignOnVertex); put(mi, L.MI.NEE, mesh.nEdgesOnEdge);
  const edgesOnEdge = Int32Array.from({ length: MAX_EDGES_ON_EDGE * E }, (_, slot) => (slot % MAX_EDGES_ON_EDGE < mesh.nEdgesOnEdge[Math.floor(slot / MAX_EDGES_ON_EDGE)] ? mesh.edgesOnEdge[mesh.maxEdgesOnEdge * Math.floor(slot / MAX_EDGES_ON_EDGE) + slot % MAX_EDGES_ON_EDGE] : Math.floor(slot / MAX_EDGES_ON_EDGE)));
  put(mi, L.MI.EOE, edgesOnEdge);
  const mf = new Float32Array(L.MF.total);
  put(mf, L.MF.AREA, mesh.areaCell); put(mf, L.MF.ATRI, mesh.areaTriangle); put(mf, L.MF.DC, mesh.dcEdge); put(mf, L.MF.DV, mesh.dvEdge); put(mf, L.MF.FV, mesh.fVertex);
  put(mf, L.MF.KAV, mesh.kiteAreasOnVertex);
  put(mf, L.MF.PVW, Float64Array.from({ length: MAX_EDGES_ON_EDGE * E }, (_, slot) => {
    const e = Math.floor(slot / MAX_EDGES_ON_EDGE), s = slot % MAX_EDGES_ON_EDGE;
    if (s >= mesh.nEdgesOnEdge[e]) return 0;
    const source = mesh.maxEdgesOnEdge * e + s;
    return mesh.weightsOnEdge[source] * mesh.dvEdge[mesh.edgesOnEdge[source]];
  }));
  put(mf, L.MF.NEDGE, mesh.nEdge); put(mf, L.MF.LAT, mesh.latCell); put(mf, L.MF.XC, mesh.xCell);
  if (surfaceGeopotential) put(mf, L.MF.PHIS, surfaceGeopotential);
  if (surfaceGeopotential) put(mf, L.MF.GPHIS, Float64Array.from({ length: E }, (_, e) => (surfaceGeopotential[mesh.cellsOnEdge[2 * e + 1]] - surfaceGeopotential[mesh.cellsOnEdge[2 * e]]) / mesh.dcEdge[e]));
  const lv = new Float32Array(L.LV.total);
  put(lv, L.LV.SL, sigmaLower); put(lv, L.LV.SU, sigmaUpper); put(lv, L.LV.DS, dSigma); put(lv, L.LV.SM, sigmaMid); put(lv, L.LV.TOP, topRate);
  const cl = Float64Array.from(sigmaLower, (s) => Math.pow(s, kappa));
  const cm = Float64Array.from(sigmaLower, (s, k) => (Math.pow(s, 1 + kappa) - Math.pow(sigmaUpper[k], 1 + kappa)) / ((1 + kappa) * dSigma[k]));
  const cd = Float64Array.from(sigmaLower, (s, k) => (kappa / (1 + kappa)) * (Math.pow(s, 1 + kappa) - Math.pow(sigmaUpper[k], 1 + kappa)) / dSigma[k]);
  const ca = Float64Array.from(sigmaLower, (s, k) => (k < K - 1 ? cm[k + 1] - cl[k] : 0));
  const cb = Float64Array.from(sigmaLower, (s, k) => cl[k] - cm[k]);
  const ct = Float64Array.from(sigmaLower, (s, k) => (k < K - 1 ? cb[k] / (cm[k + 1] - cm[k]) : 0));
  put(lv, L.LV.CL, cl); put(lv, L.LV.CM, cm); put(lv, L.LV.CD, cd); put(lv, L.LV.CA, ca); put(lv, L.LV.CB, cb); put(lv, L.LV.CT, ct);
  const thetaRef = referenceTheta ? Float64Array.from(referenceTheta) : new Float64Array(K);
  const gr = new Float64Array(K), gabs = new Float64Array(K);
  gr[K - 1] = cp * thetaRef[K - 1] * cb[K - 1];
  for (let k = K - 2; k >= 0; k--) gr[k] = cp * (thetaRef[k + 1] * ca[k] + thetaRef[k] * cb[k]);
  gabs[K - 1] = gr[K - 1];
  for (let k = K - 2; k >= 0; k--) gabs[k] = gabs[k + 1] + gr[k];
  put(lv, L.LV.GR, gr); put(lv, L.LV.GABS, gabs); put(lv, L.LV.SHAPE, shape); put(lv, L.LV.OZ, ozoneFraction); put(lv, L.LV.GASE, gasEmissivity); put(lv, L.LV.AER, aerosolFraction); put(lv, L.LV.OZS, ozoneShare);

  const buffers = {
    MI: storageBuffer(device, mi), MF: storageBuffer(device, mf), LV: storageBuffer(device, lv),
    S: emptyBuffer(device, 4 * L.S.total), T: emptyBuffer(device, 4 * L.S.total),
    K1: emptyBuffer(device, 4 * L.S.total), K2: emptyBuffer(device, 4 * L.S.total), K3: emptyBuffer(device, 4 * L.S.total), K4: emptyBuffer(device, 4 * L.S.total),
    D: emptyBuffer(device, 4 * L.D.total), P: storageBuffer(device, new Float32Array(8)), PH: emptyBuffer(device, 4 * L.PH.total),
    FR: emptyBuffer(device, 4 * L.FR.total), FP: storageBuffer(device, new Float32Array(8)), PR: emptyBuffer(device, 32 * RING_SLOTS),
  };
  for (const [name, b] of Object.entries(buffers)) b.label = name;

  const layout = device.createBindGroupLayout({ entries: Array.from({ length: 8 }, (_, binding) => ({ binding, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } })) });
  const pipelineLayout = device.createPipelineLayout({ bindGroupLayouts: [layout] });
  const head = prelude(L, { kappa, cp, p0, g, R, dragCoefficient, gustiness, physics: physicsConstants({ ...phys, kTop, stratusLayer: nearestLayer(sigmaMid, phys.stratusSigma), stabilityLayer: nearestLayer(sigmaMid, STABILITY_SIGMA) }) });
  const preludeConstants = head.slice(0, head.indexOf('@group(0) @binding(0)'));
  let meshSpacing = 0;
  for (let e = 0; e < E; e++) meshSpacing += mesh.dcEdge[e];
  meshSpacing /= E;
  const divergenceStep = divergenceDamping * meshSpacing * meshSpacing;
  const kernels = {};
  for (const [name, body] of Object.entries({ ...KERNELS, ...PHYSICS_KERNELS, ...FRAME_KERNELS, frameReduce: reductionKernel(REDUCED, { count: C, base: 'FR_PART', setup: REDUCED_SETUP }) })) {
    const code = head + body.replaceAll('S_TOTAL', String(L.S.total)).replaceAll('i32(id.x)', '(i32(id.x) + i32(id.y) * 4194240)');
    const module = device.createShaderModule({ code, label: name });
    kernels[name] = device.createComputePipeline({ label: name, layout: pipelineLayout, compute: { module, entryPoint: 'main' } });
  }
  const groups = new Map();
  function group(IN, OUT, D = buffers.D, P = buffers.P, MF = buffers.MF, LV = buffers.LV) {
    const key = [IN.label, OUT.label, D.label, P.label, MF.label, LV.label].join('|');
    let g = groups.get(key);
    if (!g) {
      g = device.createBindGroup({ layout, entries: [buffers.MI, MF, LV, IN, OUT, D, P, buffers.PH].map((buffer, binding) => ({ binding, resource: { buffer } })) });
      groups.set(key, g);
    }
    return g;
  }
  function dispatch(pass, name, bindGroup, count) {
    pass.setPipeline(kernels[name]);
    pass.setBindGroup(0, bindGroup);
    const groups = Math.ceil(count / WORKGROUP);
    pass.dispatchWorkgroups(Math.min(groups, 65535), Math.ceil(groups / 65535));
  }

  function tendencyPasses(pass, IN, OUT) {
    const g = group(IN, OUT);
    dispatch(pass, 'flux', g, L.KE);
    dispatch(pass, 'divergence', g, L.KC);
    dispatch(pass, 'column', g, C);
    dispatch(pass, 'vertexPi', g, V);
    dispatch(pass, 'pvVertex', g, L.KV);
    dispatch(pass, 'pvEdge', g, L.KE);
    dispatch(pass, 'cellTendency', g, C);
    dispatch(pass, 'momentum', g, E);
  }

  /*
   * Step work is recorded through `encode`/`compute`, and the eight
   * parameters of P (or of another kernel's parameter buffer) through
   * `writeParams`: the core's, the ocean's and the forcing recorder's.
   * Outside a batch each record is its own submission after a queue
   * write of the parameters. Inside `batched` every record goes into one
   * encoder, and each parameter set is staged in a host array and copied
   * into its buffer by a copy command at its place in that encoder. The
   * staged sets reach the ring buffer PR in one queue write just before
   * the encoder is submitted, so the dispatches see exactly the
   * parameters they would have seen stepping one at a time. A full ring
   * submits what it holds and starts again at its first slot, which the
   * queue orders after the work already submitted.
   */
  const staged = new Float32Array(8 * RING_SLOTS);
  let batch = null;
  function encode(record) {
    const encoder = batch ? batch.encoder : device.createCommandEncoder();
    record(encoder);
    if (!batch) device.queue.submit([encoder.finish()]);
  }
  function compute(record) {
    encode((encoder) => { const pass = encoder.beginComputePass(); record(pass); pass.end(); });
  }
  function flush(more) {
    if (batch.slots) device.queue.writeBuffer(buffers.PR, 0, staged.subarray(0, 8 * batch.slots));
    device.queue.submit([batch.encoder.finish()]);
    batch.submissions++;
    batch.slots = 0;
    batch.encoder = more ? device.createCommandEncoder() : null;
  }
  function writeParams(values, target = buffers.P) {
    if (!batch) { device.queue.writeBuffer(target, 0, values); return; }
    if (batch.slots === RING_SLOTS) flush(true);
    staged.set(values, 8 * batch.slots);
    batch.encoder.copyBufferToBuffer(buffers.PR, 32 * batch.slots, target, 0, 32);
    batch.slots++;
  }
  function clearBuffer(buffer, offset, size) {
    if (batch) batch.encoder.clearBuffer(buffer, offset, size);
    else device.queue.writeBuffer(buffer, offset, new Uint8Array(size));
  }
  async function batched(run) {
    if (batch) throw new Error('a batch is already being recorded');
    batch = { encoder: device.createCommandEncoder(), slots: 0, submissions: 0 };
    try {
      await run();
      flush(false);
      return batch.submissions;
    } finally {
      batch = null;
    }
  }

  const params = new Float32Array(8);
  function setParams(values) { params.set(values); writeParams(params); }

  function rungeKutta(dt) {
    const stage = (IN, OUT, factor, next) => {
      setParams([factor]);
      compute((pass) => {
        tendencyPasses(pass, IN, OUT);
        if (next) dispatch(pass, 'advance', group(buffers.S, next, OUT), L.S.total);
      });
    };
    stage(buffers.S, buffers.K1, dt / 2, buffers.T);
    stage(buffers.T, buffers.K2, dt / 2, buffers.T);
    stage(buffers.T, buffers.K3, dt, buffers.T);
    stage(buffers.T, buffers.K4, 0, null);
    setParams([dt / 6]);
    compute((pass) => dispatch(pass, 'combine', group(buffers.S, buffers.K1, buffers.K2, buffers.P, buffers.K3, buffers.K4), L.S.total));
  }

  async function step(dt) {
    rungeKutta(dt);
    closurePasses(dt);
    await device.queue.onSubmittedWorkDone();
  }

  function closurePasses(dt) {
    const g = group(buffers.S, buffers.K1);
    if (nu4Theta > 0) {
      setParams([dt * nu4Theta, 0]);
      compute((pass) => {
        dispatch(pass, 'lapScalar1', g, C);
        dispatch(pass, 'lapScalar2', g, C);
      });
    }
    if (nu4 > 0) {
      for (const second of [0, 1]) {
        setParams([dt * nu4, second]);
        compute((pass) => {
          dispatch(pass, 'divCurl', g, C + V);
          dispatch(pass, 'lapVelocity', g, E);
        });
      }
    }
    if (divergenceStep > 0) {
      setParams([divergenceStep, 0, nu4 > 0 ? 1 : 0]);
      compute((pass) => {
        dispatch(pass, 'divCurl', g, C);
        dispatch(pass, 'divergenceDamp', g, E);
      });
    }
  }

  let stepCount = 0;
  const hooks = { beforePhysics: null };
  async function stepModel(dt, time) {
    const sun = sunDirection(time);
    rungeKutta(dt);
    stepCount++;
    if (hooks.beforePhysics) await hooks.beforePhysics(dt, stepCount);
    const g = group(buffers.S, buffers.K1);
    setParams([dt, 0, sun[0], sun[1], sun[2], (time % YEAR) / YEAR]);
    compute((pass) => {
      dispatch(pass, 'physics', g, C);
      dispatch(pass, 'pblDiagnose', g, C);
    });
    closurePasses(dt);
    setParams([dt, 0, sun[0], sun[1], sun[2], (time % YEAR) / YEAR]);
    compute((pass) => {
      dispatch(pass, 'adjust', g, C);
      dispatch(pass, 'mixMomentum', g, E);
      dispatch(pass, 'dissipationHeat', g, L.KC);
      dispatch(pass, 'dissipationClear', g, L.KE);
    });
  }

  const retained = { land: null, drag: null, soil: null, snow: null, vegetation: null, concentration: null, mlmSubsidence: null, mlmHeight: null, mlmGate: null, convectiveRain: null, largeScaleRain: null, meanAbsorbedSolar: null, meanOutgoingLongwave: null, meanPlanetaryAlbedo: null, meanShortwaveCloudEffect: null, meanLongwaveCloudEffect: null, boundaryDepth: null, mixingTop: null, regime: null, buoyancyFlux: null };
  function uploadPhysics({ capacity = null, oceanFlux = null, land, drag, soil, snow, vegetation, concentration, mlmSubsidence, mlmHeight, mlmGate, convectiveRain, largeScaleRain, meanAbsorbedSolar, meanOutgoingLongwave, meanPlanetaryAlbedo, meanShortwaveCloudEffect, meanLongwaveCloudEffect, boundaryDepth, mixingTop, regime, buoyancyFlux } = {}) {
    for (const [name, value] of Object.entries({ land, drag, soil, snow, vegetation, concentration, mlmSubsidence, mlmHeight, mlmGate, convectiveRain, largeScaleRain, meanAbsorbedSolar, meanOutgoingLongwave, meanPlanetaryAlbedo, meanShortwaveCloudEffect, meanLongwaveCloudEffect, boundaryDepth, mixingTop, regime, buoyancyFlux })) if (value !== undefined) retained[name] = value;
    const ph = new Float32Array(L.PH.total);
    for (let i = 0; i < C; i++) {
      const lat = mesh.latCell[i];
      ph[L.PH.TAU + i] = phys.tauEquator + (phys.tauPole - phys.tauEquator) * Math.sin(lat) ** 2;
      ph[L.PH.CAP + i] = capacity ? capacity[i] : phys.slabHeatCapacity;
      ph[L.PH.OFLUX + i] = oceanFlux ? oceanFlux[i] : 0;
      ph[L.PH.LAND + i] = retained.land ? retained.land[i] : 0;
      ph[L.PH.DRAG + i] = retained.drag ? retained.drag[i] : dragCoefficient;
      ph[L.PH.SOIL + i] = retained.soil ? retained.soil[i] : 0;
      ph[L.PH.SNOW + i] = retained.snow ? retained.snow[i] : 0;
      ph[L.PH.CONC + i] = retained.concentration ? retained.concentration[i] : 0;
      ph[L.PH.VEG + i] = retained.vegetation ? retained.vegetation[i] : 0;
      ph[L.PH.SURF + i] = retained.surface ? retained.surface[i] : 0;
      ph[L.PH.MLMSUB + i] = retained.mlmSubsidence ? retained.mlmSubsidence[i] : 0;
      ph[L.PH.MLMH + i] = retained.mlmHeight ? retained.mlmHeight[i] : 0;
      ph[L.PH.MLMGATE + i] = retained.mlmGate ? retained.mlmGate[i] : UNDECIDED;
      ph[L.PH.CONVMEAN + i] = retained.convectiveRain ? retained.convectiveRain[i] : 0;
      ph[L.PH.CONDMEAN + i] = retained.largeScaleRain ? retained.largeScaleRain[i] : 0;
      ph[L.PH.ASRMEAN + i] = retained.meanAbsorbedSolar ? retained.meanAbsorbedSolar[i] : 0;
      ph[L.PH.OLRMEAN + i] = retained.meanOutgoingLongwave ? retained.meanOutgoingLongwave[i] : 0;
      ph[L.PH.ALBMEAN + i] = retained.meanPlanetaryAlbedo ? retained.meanPlanetaryAlbedo[i] : 0;
      ph[L.PH.SWCREMEAN + i] = retained.meanShortwaveCloudEffect ? retained.meanShortwaveCloudEffect[i] : 0;
      ph[L.PH.LWCREMEAN + i] = retained.meanLongwaveCloudEffect ? retained.meanLongwaveCloudEffect[i] : 0;
      ph[L.PH.DEPTH + i] = retained.boundaryDepth ? retained.boundaryDepth[i] : 0;
      ph[L.PH.MIXTOP + i] = retained.mixingTop ? retained.mixingTop[i] : 0;
      ph[L.PH.REGIME + i] = retained.regime ? retained.regime[i] : 0;
      ph[L.PH.BUOY + i] = retained.buoyancyFlux ? retained.buoyancyFlux[i] : 0;
    }
    device.queue.writeBuffer(buffers.PH, 0, ph);
  }
  function uploadLand({ soil, snow, vegetation, surface = null }) {
    retained.soil = soil; retained.snow = snow; retained.vegetation = vegetation; retained.surface = surface;
    device.queue.writeBuffer(buffers.PH, 4 * L.PH.SOIL, Float32Array.from(soil));
    device.queue.writeBuffer(buffers.PH, 4 * L.PH.SNOW, Float32Array.from(snow));
    device.queue.writeBuffer(buffers.PH, 4 * L.PH.VEG, Float32Array.from(vegetation));
    device.queue.writeBuffer(buffers.PH, 4 * L.PH.SURF, surface ? Float32Array.from(surface) : new Float32Array(C));
    device.queue.writeBuffer(buffers.PH, 4 * L.PH.RUNOFF, new Float32Array(C));
  }
  function uploadIce(concentration) {
    retained.concentration = concentration;
    device.queue.writeBuffer(buffers.PH, 4 * L.PH.CONC, Float32Array.from(concentration));
  }
  function downloadPhysics() {
    return readRanges(device, buffers.PH, [{ offset: 0, length: L.PH.total }]).then(([ph]) => Object.fromEntries(Object.entries(L.PH).filter(([k]) => k !== 'total').map(([k, off]) => [k, ph.subarray(off)])));
  }

  async function tendency(IN = buffers.S, OUT = buffers.K1) {
    const encoder = device.createCommandEncoder();
    const pass = encoder.beginComputePass();
    tendencyPasses(pass, IN, OUT);
    pass.end();
    device.queue.submit([encoder.finish()]);
    await device.queue.onSubmittedWorkDone();
  }

  const names = ['PI', 'TH', 'U', 'TS', 'Q', 'QC', 'ICE'];
  function upload(state) {
    const packed = new Float32Array(L.S.total);
    state.forEach((array, a) => packed.set(array, L.S[names[a]]));
    device.queue.writeBuffer(buffers.S, 0, packed);
    device.queue.writeBuffer(buffers.T, 0, packed);
    const zero = new Float32Array(L.S.total);
    for (const b of [buffers.K1, buffers.K2, buffers.K3, buffers.K4]) device.queue.writeBuffer(b, 0, zero);
    device.queue.writeBuffer(buffers.D, 0, new Float32Array(L.D.total));
  }
  function uploadSurfaceTemperature(surfaceT) {
    const values = Float32Array.from(surfaceT);
    device.queue.writeBuffer(buffers.S, 4 * L.S.TS, values);
    device.queue.writeBuffer(buffers.T, 4 * L.S.TS, values);
  }
  function download(buffer = buffers.S) {
    const lengths = [C, L.KC, L.KE, C, L.KC, L.KC, C];
    return readRanges(device, buffer, [{ offset: 0, length: L.S.total }]).then(([packed]) => names.map((name, a) => Float64Array.from(packed.subarray(L.S[name], L.S[name] + lengths[a]))));
  }
  async function downloadDiagnostics() {
    const d = await readBuffer(device, buffers.D, 4 * L.D.total);
    const out = Object.fromEntries(Object.entries(L.D).filter(([k]) => k !== 'total').map(([k, off]) => [k, d.subarray(off)]));
    for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) out.GEO[k * C + i] += gabs[k];
    return out;
  }

  /*
   * Queues the frame kernels and the read-backs of exactly the named
   * fields, plus the diagnostics' partial sums and each cell's runoff
   * since the last diagnostics when asked; the promise resolves to
   * { fields, sums, runoff } once the copies land.
   */
  const FIELDS = {
    temperature: ['FR', 'T', 1], height: ['FR', 'Z', 1], humidity: ['FR', 'RH', 1], speed: ['FR', 'SPD', 1], wind: ['FR', 'WIND', 3],
    dewPoint: ['FR', 'DP', 1], wetBulb: ['FR', 'WB', 1], misery: ['FR', 'MI', 1], vertical: ['FR', 'WM', 1],
    water: ['FR', 'TPW', 1], cloud: ['FR', 'TCW', 1], mslp: ['FR', 'MSLP', 1], rain: ['FR', 'RAIN', 1], ps: ['S', 'PI', 1], ice: ['S', 'ICE', 1],
    albedo: ['PH', 'ADIF', 1], shortwave: ['PH', 'SWDN', 1], longwave: ['PH', 'OLR', 1], soil: ['PH', 'SOIL', 1], snow: ['PH', 'SNOW', 1], vegetation: ['PH', 'VEG', 1],
    concentration: ['PH', 'CONC', 1],
  };
  const frameParams = new Float32Array(8);
  function frame({ pressure = 0, keep = 1, keepVertical = 0, fields = [], diagnostics = false, rainScale = 0, meanScale = 0 } = {}) {
    const wanted = fields.filter((name) => name in FIELDS), vertical = wanted.includes('vertical');
    frameParams.set([pressure, keep, diagnostics ? 1 : 0, vertical ? 1 : 0, keepVertical, rainScale, meanScale]);
    device.queue.writeBuffer(buffers.FP, 0, frameParams);
    const g = group(buffers.S, buffers.FR, buffers.D, buffers.FP);
    const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
    if (wanted.some((name) => FIELDS[name][0] === 'FR' && name !== 'rain')) dispatch(pass, 'frameFields', g, C);
    if (vertical) dispatch(pass, 'frameVertical', g, C);
    if (diagnostics) dispatch(pass, 'frameReduce', g, C);
    dispatch(pass, 'frameRain', g, C);
    pass.end();
    device.queue.submit([encoder.finish()]);
    const reads = {};
    for (const name of wanted) { const [buffer, part, n] = FIELDS[name]; (reads[buffer] ??= []).push({ name, offset: L[buffer][part], length: n * C }); }
    if (diagnostics) (reads.FR ??= []).push({ name: 'sums', offset: L.FR.PART, length: REDUCED.length * groupsOf(C) }, { name: 'runoff', offset: L.FR.RDONE, length: C });
    return Promise.all(Object.entries(reads).map(([buffer, ranges]) => readRanges(device, buffers[buffer], ranges).then((views) => views.map((view, n) => [ranges[n].name, view])))).then((lists) => {
      const out = { fields: {}, sums: null, runoff: null };
      for (const [name, view] of lists.flat()) {
        if (name === 'sums') out.sums = finishReduction(REDUCED, view, C);
        else if (name === 'runoff') out.runoff = view;
        else out.fields[name] = view;
      }
      return out;
    });
  }
  function clearFrame() { device.queue.writeBuffer(buffers.FR, 0, new Float32Array(L.FR.total)); }
  function setWindSpeed(windSpeed) { device.queue.writeBuffer(buffers.D, 4 * L.D.WIND, Float32Array.from(windSpeed)); }

  return { device, mesh, meshSpacing, preludeConstants, layout: L, buffers, kernels, step, stepModel, hooks, batched, encode, compute, writeParams, clearBuffer, get stepCount() { return stepCount; }, tendency, upload, uploadSurfaceTemperature, download, downloadDiagnostics, frame, clearFrame, uploadPhysics, uploadLand, uploadIce, downloadPhysics, setWindSpeed, K, C, E, V, kTop, dSigma, sigmaMid, physics: phys };
}
