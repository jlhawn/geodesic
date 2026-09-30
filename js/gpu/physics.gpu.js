import { MINIMUM_CONCENTRATION, MINIMUM_VOLUME } from '../physics/ice.module.js';
import { MIXED_LAYER_DEFAULTS, DYCOMS_LONGWAVE } from '../physics/mixedLayer.module.js';
import { DECK_CLOUD_LEVELS, UNDECIDED, VISIBLE_PATH } from '../physics/radiation.module.js';
import { CLEAR_AIR } from '../physics/moist.module.js';

/*
 * The column physics of the model as WGSL, one thread per column (or per
 * edge for momentum mixing), sharing the core's bindings and layouts:
 * the three-band gray radiation with clouds, the stratocumulus deck
 * (the EIS fit or the mixed-layer model) and the zenith/diffuse
 * surface reflection, bulk surface fluxes, the zero-layer sea ice and
 * its concentration, the boundary-layer diagnosis, and the adjustment
 * phase — boundary-layer mixing by the implicit tridiagonal solve,
 * saturation adjustment, the triggered, entraining Betts–Miller
 * convection with its anvil and downdraft, autoconversion and the rain's
 * fall, the filler and the dry convective adjustment. Each is a line-by-line port of the JavaScript module it
 * names; the physics reads the Exner ratios the last RK4 stage left in
 * the diagnostic buffer, as the CPU does, and the adjustment
 * re-diagnoses the column first. The cloud water's shortwave absorption
 * sits behind the constant CLOUD_SW, in the two-stream and in a pass of
 * its own for the mixed layer's sunlight, so that with
 * cloudSolarAbsorption 0 the kernel compiles to the purely scattering
 * one bit for bit: the shader compiler reassociates, and sharing terms
 * between the two changes the rounding.
 */
export function physicsConstants(o) {
  if (o.stratusIndex !== 'eis' && o.stratusIndex !== 'ectei') throw new Error(`stratusIndex must be 'eis' or 'ectei', not ${o.stratusIndex}`);
  const m = { ...MIXED_LAYER_DEFAULTS, cloudLevels: DECK_CLOUD_LEVELS, ...o.mixedLayer };
  if (m.closure !== 'radiative' && m.closure !== 'buoyancy') throw new Error(`closure must be 'radiative' or 'buoyancy', not ${m.closure}`);
  if (m.drizzle) throw new Error('the GPU mixed-layer deck runs without drizzle');
  if (![0, 1, 2].includes(o.subsidenceSmoothing)) throw new Error(`subsidenceSmoothing must be 0, 1 or 2, not ${o.subsidenceSmoothing}`);
  return `
const S0: f32 = ${o.solarConstant}; const STEFAN: f32 = 5.670374419e-8; const LHEAT: f32 = ${o.latentHeat}; const EPSILON: f32 = 0.622; const RVAP: f32 = ${o.R / 0.622};
const PDF_COVER: bool = ${o.cloudCover === 'pdf'}; const VISIBLE_PATH: f32 = ${VISIBLE_PATH}; const RHC: f32 = ${o.criticalHumidity}; const RHC_BL: f32 = ${o.boundaryCriticalHumidity}; const COVER_FLOOR: f32 = ${o.coverFloor ?? 0.01};
const CLOUD_ABS: f32 = ${o.cloudAbsorption}; const CLOUD_SCAT: f32 = ${o.cloudScattering}; const CLOUD_SW: f32 = ${o.cloudSolarAbsorption}; const WINDOW: f32 = ${o.window}; const GAS_FRAC: f32 = ${o.gasFraction};
const STRATUS: bool = ${!!o.stratus}; const ECTEI: bool = ${o.stratusIndex === 'ectei'}; const STRATUS_SCALE: f32 = ${o.stratusScale}; const STRATUS_MAX: f32 = ${o.stratusWaterMax}; const STRATUS_K: i32 = ${o.stratusLayer}; const STABILITY_K: i32 = ${o.stabilityLayer};
const VAPOR_FRAC: f32 = ${1 - o.window - o.gasFraction}; const OZONE_ABS: f32 = ${o.ozoneAbsorption}; const VAPOR_ABS: f32 = ${o.vaporAbsorption}; const CEX: f32 = ${o.exchangeCoefficient};
const VCOUP: f32 = ${o.vaporCoupling}; const COUPLED: bool = ${o.vaporCoupling > 0}; const SKYLIGHT: f32 = ${o.skylight}; const DIFFUSE_MU: f32 = 0.6;
const ALB_ICE: f32 = ${o.iceAlbedo}; const FULLALB: f32 = ${o.fullAlbedoThickness}; const ALB_DIF_WATER: f32 = ${o.diffuseWaterAlbedo};
const ALB_ICESNOW: f32 = ${o.iceSnowAlbedo}; const FULLSNOW_ICE: f32 = ${o.iceFullSnow}; const KSNOW: f32 = ${o.snowConductivity}; const RHOSNOW: f32 = ${o.snowDensity}; const RHOICE: f32 = ${o.iceDensity}; const RHOWATER: f32 = ${o.waterDensity};
const FREEZING: f32 = 271.35; const MELTING: f32 = 273.15; const SKINC: f32 = ${o.skinHeatCapacity}; const COND: f32 = ${o.conductivity}; const HMIN: f32 = ${o.minimumThickness}; const LATENT_ICE: f32 = ${o.iceDensity * o.latentHeatFusion};
const LEADC: f32 = ${o.leadClosing}; const LEADX: f32 = ${o.leadExchange}; const MIN_CONC: f32 = ${MINIMUM_CONCENTRATION}; const MIN_VOLUME: f32 = ${MINIMUM_VOLUME};
const RELAX: f32 = ${o.relaxationTime}; const RH_REF: f32 = ${o.referenceHumidity}; const AUTO_T: f32 = ${o.autoconversionThreshold}; const AUTO_R: f32 = ${o.autoconversionRate}; const CLOUD_LIFE: f32 = ${o.cloudLifetime};
const DETRAIN: f32 = ${o.detrainment}; const ANVIL: f32 = ${o.anvilDepth}; const RAIN_EVAP: f32 = ${o.rainEvaporation};
const AUTO_BL: bool = ${o.autoconversionFloor === 'boundaryLayer'}; const CLEAR_AIR: f32 = ${CLEAR_AIR}; const PARCEL_DEPTH: f32 = ${o.parcelDepth}; const ENTRAIN: f32 = ${o.entrainmentRate}; const CAPE_MIN: f32 = ${o.capeThreshold}; const CIN_MAX: f32 = ${o.inhibitionThreshold}; const ACT_MEM: f32 = ${o.activityMemory}; const SHALLOW_TOP: f32 = ${o.shallowTop}; const DOWNDRAFT: f32 = ${o.downdraftEvaporation}; const SHALLOW_RH: f32 = ${o.shallowHumidity};
const BL_PARCEL: bool = ${o.boundaryParcel !== false}; const FROM_SURFACE: bool = ${o.adjustFrom === 'surface'}; const DECK_VETO: bool = ${o.deckVeto !== false}; const EVAP_IN_CLOUD: bool = ${!!o.evaporationInCloud}; const AUTO_NONE: bool = ${o.autoconversionFloor === 'none'};
const EAGER: bool = ${o.shallowTrigger != null}; const EAGER_CAPE: f32 = ${o.shallowTrigger ? o.shallowTrigger.cape : 0}; const EAGER_CIN: f32 = ${o.shallowTrigger ? o.shallowTrigger.inhibition : 0}; const EAGER_TOP: f32 = ${o.shallowTrigger ? o.shallowTrigger.top ?? o.shallowTop : 0}; const DRAFT_MASS: bool = ${o.downdraftSpread === 'mass'}; const PARCEL_VIRT: f32 = ${o.virtualBuoyancy === false ? 0 : 'VIRT'};
const BL_WE: f32 = ${o.entrainmentVelocity ?? 0}; const RIC: f32 = ${o.richardsonCritical}; const KARMAN: f32 = ${o.vonKarman}; const STABILITY: bool = ${o.stability ? 'true' : 'false'}; const KTOP: i32 = ${o.kTop};
const LANDED: bool = ${!!o.landed}; const LANDC: f32 = ${o.landHeatCapacity}; const BUCKET: f32 = ${o.bucketCapacity}; const WETT: f32 = ${o.wetnessThreshold}; const ALB_LAND: f32 = ${o.landAlbedo}; const VEGETATED: bool = ${!!o.vegetation}; const ALB_BARE: f32 = ${o.bareAlbedo}; const ALB_VEG: f32 = ${o.vegetatedAlbedo}; const ROOTCAP: f32 = ${o.rootZoneCapacity};
const MLM_DECK: bool = ${!!o.mixedLayerDeck}; const STRATUS_SOLAR: bool = ${!!o.stratusSolar}; const MLM_SUBSIDENCE: f32 = ${o.stratusSubsidence}; const MLM_MININV: f32 = ${o.minimumInversion}; const MLM_MEMORY: f32 = ${o.subsidenceMemory};
const MLM_LEVELS: i32 = ${m.cloudLevels}; const MLM_NODES: i32 = ${m.cloudLevels + 1}; const MLM_BUOYANCY: bool = ${m.closure === 'buoyancy'}; const MLM_DELTA: f32 = 1.0 / EPSILON - 1.0; const MLM_LC: f32 = LHEAT / CP;
const MLM_A1: f32 = ${m.entrainmentEfficiency}; const MLM_A2: f32 = ${m.evaporativeEnhancement}; const MLM_AMAX: f32 = ${m.maximumEfficiency}; const MLM_WEMAX: f32 = ${m.maximumEntrainment}; const MLM_MINJUMP: f32 = ${m.minimumJump};
const MLM_ONSET: f32 = ${m.decouplingOnset}; const MLM_DRATIO: f32 = ${m.decoupledRatio}; const MLM_DCOVER: f32 = ${m.decoupledCover}; const DYC_F0: f32 = ${DYCOMS_LONGWAVE.F0}; const DYC_F1: f32 = ${DYCOMS_LONGWAVE.F1}; const DYC_K: f32 = ${DYCOMS_LONGWAVE.kappa};
const MLM_PASSES: i32 = ${o.subsidenceSmoothing}; const MLM_PROGNOSTIC: bool = ${o.prognosticHeight ? 'true' : 'false'}; const MLM_GATEMEM: f32 = ${o.gateMemory}; const MLM_UNDECIDED: f32 = ${UNDECIDED}; const MLM_HMEM: f32 = ${m.heightMemory}; const MLM_HMAX: f32 = ${m.maximumHeight};
const ALB_ICESHEET: f32 = ${o.iceSheetAlbedo}; const SURFCAP: f32 = ${o.surfaceCapacity}; const PERCT: f32 = ${o.percolationTime}; const RSTOM: f32 = ${o.stomatalResistance}; const GROWCOLD: f32 = ${o.growthColdest}; const GROWWARM: f32 = ${o.growthWarmest}; const VEG_DRY: f32 = ${o.dryWetness}; const VEG_WET: f32 = ${o.wetWetness}; const VEG_GROW: f32 = ${o.growthTime}; const VEG_DECLINE: f32 = ${o.declineTime}; const VEG_SNOW: f32 = ${o.snowDeclineTime}; const ALB_SNOW: f32 = ${o.snowAlbedo}; const FULLSNOW: f32 = ${o.fullSnow}; const LFUS: f32 = ${o.latentHeatFusion};
`;
}

export const PHYSICS_FUNCTIONS = `
fn esat(T: f32) -> f32 { return 611.2 * exp(17.67 * (T - 273.15) / (T - 29.65)); }
fn qsat(T: f32, p: f32) -> f32 { let es = esat(T); let dry = p - (1.0 - EPSILON) * es; return select(1.0, EPSILON * es / dry, dry > 0.0); }
fn openWaterAlbedo(mu: f32) -> f32 { return 0.026 / (pow(mu, 1.7) + 0.065) + 0.15 * (mu - 0.1) * (mu - 0.5) * (mu - 1.0); }
fn surfaceAlbedo(h: f32, water: f32, snow: f32) -> f32 {
  if (h <= 0.0) { return water; }
  let bare = water + (ALB_ICE - water) * min(1.0, h / FULLALB);
  return bare + (ALB_ICESNOW - bare) * min(1.0, snow / FULLSNOW_ICE);
}
fn cellWind(i: i32, k: i32) -> vec3<f32> {
  var w = vec3<f32>(0.0, 0.0, 0.0);
  for (var m = 0; m < MAXE; m++) {
    let e = MI[EOC + MAXE * i + m];
    let s = abs(f32(MI[ESC + MAXE * i + m])) * 0.5 * MF[F_DC + e] * MF[F_DV + e] * IN[S_U + k * E + e];
    w += s * vec3<f32>(MF[F_NEDGE + 3 * e], MF[F_NEDGE + 3 * e + 1], MF[F_NEDGE + 3 * e + 2]);
  }
  return w / MF[F_AREA + i];
}
fn band(fraction: f32, eps: ptr<function, array<f32, K>>, temperature: ptr<function, array<f32, K>>, netFlux: ptr<function, array<f32, K>>, surfaceEmission: f32) -> vec2<f32> {
  var down = 0.0;
  for (var k = 0; k < K; k++) {
    let t = (*temperature)[k]; let e = (*eps)[k];
    let emitted = fraction * e * STEFAN * t * t * t * t;
    (*netFlux)[k] += e * down - 2.0 * emitted;
    down = down * (1.0 - e) + emitted;
  }
  var up = fraction * surfaceEmission;
  for (var k = K - 1; k >= 0; k--) {
    let t = (*temperature)[k]; let e = (*eps)[k];
    let emitted = fraction * e * STEFAN * t * t * t * t;
    (*netFlux)[k] += e * up;
    up = up * (1.0 - e) + emitted;
  }
  return vec2<f32>(up, down);
}
fn condensationLevel(T: f32, q: f32, p: f32) -> vec2<f32> {
  let e = q * p / (EPSILON + (1.0 - EPSILON) * q);
  let y = log(e / 611.2);
  let dewPoint = (273.15 * 17.67 - 29.65 * y) / (17.67 - y);
  if (dewPoint >= T) { return vec2<f32>(T, p); }
  let lclT = 1.0 / (1.0 / (dewPoint - 56.0) + log(T / dewPoint) / 800.0) + 56.0;
  return vec2<f32>(lclT, p * pow(lclT / T, 1.0 / KAPPA));
}
fn moistAdiabat(T: f32, qs: f32) -> f32 { return GRAV / CP * (1.0 + LHEAT * qs / (RGAS * T)) / (1.0 + LHEAT * LHEAT * qs / (CP * RVAP * T * T)); }
fn inversionStrength(stability: f32, lowerT: f32, upperT: f32, depth: f32) -> f32 {
  let T = 0.5 * (lowerT + upperT);
  return stability - (GRAV / CP - moistAdiabat(T, qsat(T, 85000.0))) * depth;
}
fn entrainmentIndex(inversion: f32, lowerQ: f32, upperQ: f32) -> f32 { return inversion - 0.23 * LHEAT / CP * (lowerQ - upperQ); }
fn deckWater(lclT: f32, lclP: f32, thickness: f32) -> f32 {
  if (thickness <= 0.0) { return 0.0; }
  let qs = qsat(lclT, lclP);
  let lapse = lclP / (RGAS * lclT) * qs * (LHEAT * moistAdiabat(lclT, qs) / (RVAP * lclT * lclT) - GRAV / (RGAS * lclT));
  return min(STRATUS_MAX, STRATUS_SCALE * 0.5 * lapse * thickness * thickness);
}
fn cloudKeep(path: f32) -> f32 {
  if (CLOUD_SW > 0.0) { return exp(-CLOUD_SW * path); }
  return 1.0;
}
fn shortwave(cloudDepth: f32, keep: f32, mu: f32, adir: f32, adif: f32) -> vec4<f32> {
  let reflectance = select(0.0, cloudDepth / (cloudDepth + 2.0 * mu), mu > 0.0 && cloudDepth > 0.0);
  var direct = (1.0 - SKYLIGHT) * select(1.0, exp(-cloudDepth / mu), cloudDepth > 0.0 && mu > 0.0);
  var diffuse = 1.0 - reflectance - direct;
  var returned = select(0.0, cloudDepth / (cloudDepth + 2.0 * DIFFUSE_MU), cloudDepth > 0.0);
  if (CLOUD_SW > 0.0) { direct *= keep; diffuse *= keep; returned *= keep; }
  let upward = adir * direct + adif * diffuse;
  let reflections = returned * upward / (1.0 - adif * returned);
  return vec4<f32>((1.0 - adir) * direct + (1.0 - adif) * (diffuse + reflections), direct + diffuse + reflections, direct, (1.0 - keep) * (1.0 + upward / (1.0 - adif * returned)));
}
fn moistLapse(temperature: f32, pressure: f32) -> f32 {
  let qs = qsat(temperature, pressure);
  return (RGAS * temperature + LHEAT * qs) / (CP + LHEAT * LHEAT * EPSILON * qs / (RGAS * temperature * temperature));
}
fn saturatedTemperature(energy: f32, pressure: f32, guess: f32, humidity: f32) -> f32 {
  var t = guess;
  for (var n = 0; n < 4; n++) {
    let qs = humidity * qsat(t, pressure);
    t -= (CP * t + LHEAT * qs - energy) / (CP + LHEAT * LHEAT * qs / (RVAP * t * t));
  }
  return t;
}
fn relaxedFraction(x: f32) -> f32 {
  if (x < 1e-2) { return x * (1.0 - 0.5 * x * (1.0 - x / 3.0 * (1.0 - 0.25 * x))); }
  return 1.0 - exp(-x);
}
fn thomas(n: i32, upper: ptr<function, array<f32, K>>, lower: ptr<function, array<f32, K>>, rhs: ptr<function, array<f32, K>>) {
  var gain: array<f32, K>;
  var denominator = 1.0 + (*upper)[0] + (*lower)[0];
  gain[0] = -(*lower)[0] / denominator;
  (*rhs)[0] = (*rhs)[0] / denominator;
  for (var j = 1; j < n; j++) {
    denominator = 1.0 + (*upper)[j] + (*lower)[j] + (*upper)[j] * gain[j - 1];
    gain[j] = -(*lower)[j] / denominator;
    (*rhs)[j] = ((*rhs)[j] + (*upper)[j] * (*rhs)[j - 1]) / denominator;
  }
  for (var j = n - 2; j >= 0; j--) { (*rhs)[j] -= gain[j] * (*rhs)[j + 1]; }
}
`;

/*
 * The mixed-layer deck of radiation.module.js and mixedLayer.module.js
 * (drizzle off), line by line except where single precision needs
 * another form: the Newton iterations stop at a single-precision
 * tolerance, cloud base solves ln(q_s/q_t) = 0, the step advances θ_l and
 * q_t by the jump form θ_l' = θ_l + dt (w_e Δθ_l + H)/h' of the same flux
 * update, and the running means of the subsidence and the gate and the
 * carried height's relaxation move by 1 − e^(−dt/memory) from its series
 * when that is small. The heights are above the surface, as the core's
 * geopotential here excludes the terrain's, so the deck's mlmTop is its
 * h itself.
 */
const MIXED_LAYER_WGSL = `
struct MlmAir { T: f32, ql: f32, qv: f32, qs: f32, dqs: f32 }
struct MlmState { h: f32, thetaL: f32, qt: f32 }
struct MlmOut { lwp: f32, cover: f32, entrainment: f32, jump: f32, heat: f32, water: f32 }
struct MlmDeck { ok: bool, cover: f32, water: f32, entrainment: f32, top: f32 }
struct MlmSun { incident: f32, mu: f32, adir: f32, adif: f32, path: f32, layer: f32, clear: f32 }
fn mlmFinite(x: f32) -> bool { return (bitcast<u32>(x) & 0x7f800000u) != 0x7f800000u; }
fn mlmFresh(x: f32) -> f32 {
  if (x < 1e-2) { return x * (1.0 - 0.5 * x * (1.0 - x / 3.0 * (1.0 - 0.25 * x))); }
  return 1.0 - exp(-x);
}
fn mlmRest(i: i32, depth: f32, dt: f32) {
  let h = PH[PH_MLMH + i];
  if (!(h > 0.0)) { return; }
  var settled = depth;
  if (MLM_HMEM > 0.0) { settled = h + (depth - h) * mlmFresh(dt / MLM_HMEM); }
  PH[PH_MLMH + i] = settled;
}
fn mlmPressure(x: f32) -> f32 { return P0 * pow(x, 1.0 / KAPPA); }
fn mlmSlope(T: f32, p: f32) -> f32 {
  let e = esat(T); let dry = p - (1.0 - EPSILON) * e;
  return EPSILON * p * e * 17.67 * 243.5 / ((T - 29.65) * (T - 29.65) * dry * dry);
}
fn mlmSaturate(thetaL: f32, qt: f32, x: f32) -> MlmAir {
  let p = mlmPressure(x); let Tl = thetaL * x;
  var qs = qsat(Tl, p);
  if (qt <= qs) { return MlmAir(Tl, 0.0, qt, qs, mlmSlope(Tl, p)); }
  var T = Tl + MLM_LC * (qt - qs) / (1.0 + MLM_LC * mlmSlope(Tl, p));
  for (var it = 0; it < 30; it++) {
    qs = qsat(T, p);
    let step = (T - Tl - MLM_LC * (qt - qs)) / (1.0 + MLM_LC * mlmSlope(T, p));
    T -= step;
    if (abs(step) < 1e-6 * T) { break; }
  }
  qs = qsat(T, p);
  let ql = max(0.0, qt - qs);
  return MlmAir(T, ql, qt - ql, qs, mlmSlope(T, p));
}
fn mlmVirtual(a: MlmAir, x: f32) -> f32 { return a.T / x * (1.0 + MLM_DELTA * a.qv - a.ql); }
fn mlmCloudBase(thetaL: f32, qt: f32, piS: f32) -> f32 {
  var x = piS;
  for (var it = 0; it < 50; it++) {
    let T = thetaL * x; let p = mlmPressure(x); let e = esat(T); let dry = p - (1.0 - EPSILON) * e;
    let qs = EPSILON * e / dry;
    let dlnT = 17.67 * 243.5 / ((T - 29.65) * (T - 29.65)) * p / dry;
    let derivative = dlnT * thetaL - p / (KAPPA * x) / dry;
    let step = log(qs / qt) / derivative;
    x -= step;
    if (abs(step) < 1e-6) { break; }
  }
  return x;
}
fn mlmRadiation(below: f32, above: f32) -> f32 { return DYC_F0 * exp(-DYC_K * above) + DYC_F1 * exp(-DYC_K * below); }
fn mlmAbsorbed(sun: MlmSun, lwp: f32) -> f32 {
  let water = min(STRATUS_MAX, lwp);
  if (!(water > 0.0) || !(sun.incident > 0.0)) { return 0.0; }
  let total = sun.path + water;
  let sw = shortwave(CLOUD_SCAT * total, cloudKeep(total), sun.mu, sun.adir, sun.adif);
  return sun.incident * sw.w / total * (sun.layer + water) - sun.clear;
}
fn mlmDiagnose(s: MlmState, ps: f32, sensible: f32, evaporation: f32, thetaAbove: f32, qtAbove: f32, sun: MlmSun) -> MlmOut {
  var zNode: array<f32, MLM_NODES>; var qlNode: array<f32, MLM_NODES>; var rhoNode: array<f32, MLM_NODES>;
  var aNode: array<f32, MLM_NODES>; var bNode: array<f32, MLM_NODES>; var pathNode: array<f32, MLM_NODES>;
  let h = s.h; let thetaL = s.thetaL; let qt = s.qt;
  let piS = pow(ps / P0, KAPPA);
  let thetaVDry = thetaL * (1.0 + MLM_DELTA * qt);
  var zb = 0.0; var piB = piS;
  if (!(qt >= qsat(thetaL * piS, ps))) { piB = mlmCloudBase(thetaL, qt, piS); zb = CP * thetaVDry * (piS - piB) / GRAV; }
  let cloudy = zb < h;
  let dz = (h - zb) / f32(MLM_LEVELS);
  var lwp = 0.0; var piH = 0.0; var thetaVTop = 0.0; var qlTop = 0.0; var topA = 0.0; var topB = 0.0; var topGamma = 0.0; var topSlope = 0.0;
  if (cloudy) {
    var x = piB;
    var air = mlmSaturate(thetaL, qt, x);
    var tv = mlmVirtual(air, x);
    for (var j = 0; j <= MLM_LEVELS; j++) {
      if (j > 0) {
        let predicted = x - GRAV * dz / (CP * tv);
        let tvPredicted = mlmVirtual(mlmSaturate(thetaL, qt, predicted), predicted);
        x -= 0.5 * GRAV * dz / CP * (1.0 / tv + 1.0 / tvPredicted);
        air = mlmSaturate(thetaL, qt, x);
        tv = mlmVirtual(air, x);
      }
      let p = mlmPressure(x); let gamma = MLM_LC * air.dqs; let theta = air.T / x;
      let c = 1.0 + MLM_DELTA * air.qv - air.ql + (1.0 + MLM_DELTA) * air.T * air.dqs;
      zNode[j] = zb + f32(j) * dz; qlNode[j] = air.ql;
      rhoNode[j] = p / (RGAS * air.T * (1.0 + MLM_DELTA * air.qv - air.ql));
      aNode[j] = c / (1.0 + gamma);
      bNode[j] = c * MLM_LC / (x * (1.0 + gamma)) - theta;
      if (j == 0) { pathNode[j] = 0.0; } else { pathNode[j] = pathNode[j - 1] + 0.5 * dz * (rhoNode[j - 1] * qlNode[j - 1] + rhoNode[j] * qlNode[j]); }
      if (j == MLM_LEVELS) { topGamma = gamma; topSlope = air.dqs; }
    }
    lwp = pathNode[MLM_LEVELS]; piH = x; thetaVTop = tv; qlTop = qlNode[MLM_LEVELS]; topA = aNode[MLM_LEVELS]; topB = bNode[MLM_LEVELS];
  } else {
    piH = piS - GRAV * h / (CP * thetaVDry);
    thetaVTop = thetaVDry;
  }
  let pH = mlmPressure(piH);
  let density = (ps - pH) / (GRAV * h);
  let jumpTheta = thetaAbove - thetaL; let jumpQ = qtAbove - qt;
  let jumpVirtual = mlmVirtual(mlmSaturate(thetaAbove, qtAbove, piH), piH) - thetaVTop;
  let capped = jumpVirtual > 0.0; let jump = max(MLM_MINJUMP, jumpVirtual);
  let fluxSurface = mlmRadiation(0.0, lwp); let fluxTop = mlmRadiation(lwp, 0.0);
  let divergence = fluxTop - fluxSurface;
  let absorbed = mlmAbsorbed(sun, lwp);
  let netDivergence = divergence - absorbed;
  var efficiency = MLM_A1;
  if (cloudy && qlTop > 0.0 && capped) {
    let saturatedJump = topA * jumpTheta + topB * jumpQ;
    let demand = topSlope * piH * jumpTheta - jumpQ;
    var chi = 1.0;
    if (demand > 0.0) { chi = min(1.0, qlTop * (1.0 + topGamma) / demand); }
    efficiency = min(MLM_AMAX, MLM_A1 * (1.0 + MLM_A2 * max(0.0, chi * (1.0 - saturatedJump / jump))));
  }
  let heat0 = sensible / (density * CP); let water0 = evaporation / density;
  let heatRate = (heat0 - netDivergence / (density * CP)) / h;
  let waterRate = water0 / h;
  let aDry = 1.0 + MLM_DELTA * qt; let bDry = MLM_DELTA * thetaL; let top = min(zb, h);
  let dryJump = aDry * jumpTheta + bDry * jumpQ;
  let dry0Low = aDry * heat0 + bDry * water0; let dry0High = aDry * (heat0 - top * heatRate) + bDry * (water0 - top * waterRate);
  let dry1High = -(top / h) * dryJump;
  var I0 = 0.5 * top * (dry0Low + dry0High); var I1 = 0.5 * top * dry1High;
  if (cloudy) {
    var before0 = 0.0; var before1 = 0.0;
    for (var j = 0; j <= MLM_LEVELS; j++) {
      let z = zNode[j];
      var sunBelow = 0.0;
      if (absorbed > 0.0) { sunBelow = absorbed * pathNode[j] / lwp; }
      let flux0 = aNode[j] * (heat0 - z * heatRate - (mlmRadiation(pathNode[j], lwp - pathNode[j]) - fluxSurface - sunBelow) / (density * CP)) + bNode[j] * (water0 - z * waterRate);
      let flux1 = -(z / h) * (aNode[j] * jumpTheta + bNode[j] * jumpQ);
      if (j > 0) { I0 += 0.5 * dz * (before0 + flux0); I1 += 0.5 * dz * (before1 + flux1); }
      before0 = flux0; before1 = flux1;
    }
  }
  var entrainment = 0.0;
  if (capped) {
    if (MLM_BUOYANCY) {
      let denominator = h * jump - 2.5 * efficiency * I1;
      entrainment = MLM_WEMAX;
      if (denominator > 0.0) { entrainment = min(MLM_WEMAX, max(0.0, 2.5 * efficiency * I0 / denominator)); }
    } else {
      entrainment = min(MLM_WEMAX, max(0.0, efficiency * divergence / (density * CP * jump)));
    }
  }
  let integral = I0 + entrainment * I1;
  let dryLow = dry0Low; let dryHigh = dry0High + entrainment * dry1High;
  var negative = 0.0;
  if (dryLow < 0.0 && dryHigh < 0.0) { negative = 0.5 * top * (dryLow + dryHigh); }
  else if (dryLow < 0.0) { negative = 0.5 * dryLow * top * dryLow / (dryLow - dryHigh); }
  else if (dryHigh < 0.0) { negative = 0.5 * dryHigh * top * dryHigh / (dryHigh - dryLow); }
  let rest = integral - negative;
  var decoupled = 1.0;
  if (negative < 0.0) {
    decoupled = MLM_DCOVER;
    if (rest > 0.0) {
      let ratio = -negative / rest;
      if (ratio <= MLM_ONSET) { decoupled = 1.0; }
      else if (ratio < MLM_DRATIO) { decoupled = 1.0 - (1.0 - MLM_DCOVER) * (ratio - MLM_ONSET) / (MLM_DRATIO - MLM_ONSET); }
    }
  }
  let cover = select(0.0, decoupled, cloudy && lwp > 0.0 && capped);
  return MlmOut(lwp, cover, entrainment, jumpVirtual, entrainment * jumpTheta + heat0 - netDivergence / (density * CP), entrainment * jumpQ + water0);
}
fn mlmInterface(i: i32, m: i32) -> f32 {
  let idx = m * C + i;
  return (D[D_GEO + idx] + LV[L_GABS + m] + CP * D[D_THV + idx] * (D[D_EXM + idx] - D[D_EXL + idx - C])) / GRAV;
}
fn mlmRing(i: i32, m: i32) -> f32 {
  let n = MI[NEC + i];
  var sum = D[D_PSD + m * C + i];
  for (var s = 0; s < MAXE; s++) { if (s < n) { sum += D[D_PSD + m * C + MI[COC + MAXE * i + s]]; } }
  return sum / f32(n + 1);
}
fn mlmFlow(i: i32, m: i32) -> f32 {
  if (MLM_PASSES == 0) { return D[D_PSD + m * C + i]; }
  if (MLM_PASSES == 1) { return mlmRing(i, m); }
  let n = MI[NEC + i];
  var sum = mlmRing(i, m);
  for (var s = 0; s < MAXE; s++) { if (s < n) { sum += mlmRing(MI[COC + MAXE * i + s], m); } }
  return sum / f32(n + 1);
}
fn mlmCeiling(i: i32, floor: f32) -> f32 {
  for (var k = K - 2; k >= 1; k--) {
    let upper = (D[D_GEO + k * C + i] + LV[L_GABS + k]) / GRAV;
    if ((D[D_GEO + (k + 1) * C + i] + LV[L_GABS + k + 1]) / GRAV >= MLM_HMAX) { break; }
    if (upper > floor && D[D_THV + k * C + i] - D[D_THV + (k + 1) * C + i] >= MLM_MININV) { return upper - 1.0; }
  }
  return MLM_HMAX;
}
fn mlmColumn(i: i32, pi: f32, mixedDepth: f32, sensible: f32, evaporation: f32, dt: f32, sun: MlmSun) -> MlmDeck {
  let none = MlmDeck(false, 0.0, 0.0, 0.0, 0.0);
  let bottom = (K - 1) * C + i;
  let depth = mixedDepth + (D[D_GEO + bottom] + LV[L_GABS + K - 1]) / GRAV;
  var ceiling = MLM_HMAX;
  if (MLM_PROGNOSTIC) { ceiling = min(MLM_HMAX, mlmCeiling(i, depth)); }
  var h = depth;
  if (MLM_PROGNOSTIC && PH[PH_MLMH + i] > 0.0) { h = max(depth, min(ceiling, PH[PH_MLMH + i])); }
  var weight = 0.0; var heat = 0.0; var water = 0.0; var k = K - 1;
  for (; k >= 0; k--) {
    let idx = k * C + i;
    if (!(D[D_GEO + idx] + LV[L_GABS + k] < GRAV * h)) { break; }
    let cloud = max(0.0, IN[S_QC + idx]);
    heat += LV[L_DS + k] * (IN[S_TH + idx] - LHEAT * cloud / (CP * D[D_EXM + idx]));
    water += LV[L_DS + k] * (max(0.0, IN[S_Q + idx]) + cloud);
    weight += LV[L_DS + k];
  }
  if (k < 1) { mlmRest(i, depth, dt); return none; }
  let above = k * C + i; let aboveCloud = max(0.0, IN[S_QC + above]);
  var lowerHeight = 0.0; var lower = K; var m = K - 1;
  for (; m > k; m--) {
    let z = mlmInterface(i, m);
    if (!(z < h)) { break; }
    lowerHeight = z; lower = m;
  }
  let upperHeight = mlmInterface(i, m);
  var lowerFlow = 0.0;
  if (lower < K) { lowerFlow = mlmFlow(i, lower); }
  let flow = lowerFlow + (mlmFlow(i, m) - lowerFlow) * (h - lowerHeight) / (upperHeight - lowerHeight);
  let density = pi * LV[L_SM + m] / (RGAS * D[D_THV + m * C + i] * D[D_EXM + m * C + i]);
  let subsidence = -flow / (density * GRAV);
  let x = dt / MLM_MEMORY;
  var fresh = 1.0 - exp(-x);
  if (x < 1e-2) { fresh = x * (1.0 - 0.5 * x * (1.0 - x / 3.0 * (1.0 - 0.25 * x))); }
  let mean = PH[PH_MLMSUB + i] + (subsidence - PH[PH_MLMSUB + i]) * fresh;
  PH[PH_MLMSUB + i] = mean;
  let sinking = !(mean > -MLM_SUBSIDENCE);
  let thetaAbove = IN[S_TH + above] - LHEAT * aboveCloud / (CP * D[D_EXM + above]);
  let qtAbove = max(0.0, IN[S_Q + above]) + aboveCloud;
  let start = MlmState(h, heat / weight, water / weight);
  var now = MlmOut(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  var passed = 0.0;
  if (sinking) {
    now = mlmDiagnose(start, pi, sensible, evaporation, thetaAbove, qtAbove, sun);
    if (now.jump >= MLM_MININV) { passed = 1.0; }
  }
  var gate = passed;
  if (MLM_GATEMEM > 0.0) { gate = PH[PH_MLMGATE + i] + (passed - PH[PH_MLMGATE + i]) * mlmFresh(dt / MLM_GATEMEM); }
  PH[PH_MLMGATE + i] = gate;
  if (!(gate > MLM_UNDECIDED || (gate == MLM_UNDECIDED && passed > 0.0))) { mlmRest(i, depth, dt); return none; }
  if (!sinking) { now = mlmDiagnose(start, pi, sensible, evaporation, thetaAbove, qtAbove, sun); }
  var next = now; var top = h;
  if (dt > 0.0) {
    let deepened = h + dt * (now.entrainment + subsidence);
    top = deepened;
    if (MLM_PROGNOSTIC) { top = max(depth, min(ceiling, deepened)); }
    next = mlmDiagnose(MlmState(top, start.thetaL + dt * now.heat / deepened, start.qt + dt * now.water / deepened), pi, sensible, evaporation, thetaAbove, qtAbove, sun);
  }
  if (!(mlmFinite(next.lwp) && mlmFinite(next.cover) && mlmFinite(now.entrainment))) { mlmRest(i, depth, dt); return none; }
  PH[PH_MLMH + i] = top;
  return MlmDeck(true, next.cover, next.lwp, now.entrainment, select(0.0, top, MLM_PROGNOSTIC));
}
`;

/*
 * The sea-cell branch of the physics kernel's surface update, shared as
 * text so that a kernel with a prescribed net flux steps the same ice;
 * it expects h, T, net, contrast, ocean, capacity, cover, snow0, dt and i.
 */
export const SEA_SURFACE_WGSL = `if (h <= 0.0) {
    if (snow0 > 0.0) { T -= LFUS * snow0 / capacity; PH[PH_SNOW + i] = 0.0; }
    T += dt * (net + ocean) / capacity;
    var fresh = 0.0;
    if (T < FREEZING) {
      let volume = (FREEZING - T) * capacity / LATENT_ICE; let area = min(1.0, volume / LEADC);
      if (area >= MIN_CONC && volume >= MIN_VOLUME) { h = volume / area; T = FREEZING; fresh = area; }
    }
    PH[PH_CONC + i] = fresh;
  } else {
    var snow = snow0;
    let split = contrast - LEADX * (FREEZING - T);
    let iceFlux = net - (1.0 - cover) * split; let waterFlux = net + cover * split;
    let conduction = (FREEZING - T) / (max(h, HMIN) / COND + snow / (RHOSNOW * KSNOW));
    T += dt * (iceFlux + conduction) / SKINC;
    var thickness = h + dt * (conduction - ocean) / LATENT_ICE;
    if (T > MELTING) {
      var excess = (T - MELTING) * SKINC; T = MELTING;
      let fromSnow = min(snow, excess / LFUS);
      snow -= fromSnow; excess -= fromSnow * LFUS;
      thickness -= excess / LATENT_ICE;
    }
    let leadHeat = (1.0 - cover) * (waterFlux + ocean) * dt;
    let melted = max(0.0, -(cover * (thickness - h))) + max(0.0, leadHeat / LATENT_ICE);
    let leadIce = max(0.0, -leadHeat / LATENT_ICE);
    let remaining = cover * SKINC * (T - FREEZING) - LATENT_ICE * (cover * thickness - leadHeat / LATENT_ICE) - LFUS * cover * snow;
    var volume = cover * thickness - leadHeat / LATENT_ICE;
    let area = min(1.0, cover - melted / (2.0 * h) + (1.0 - cover) * leadIce / LEADC);
    if (area < cover) { volume += (cover - area) * (snow / RHOICE - SKINC * (T - FREEZING) / LATENT_ICE); }
    else if (area > cover) { T = FREEZING + cover * (T - FREEZING) / area; snow = cover * snow / area; }
    if (area < MIN_CONC || volume < MIN_VOLUME) { T = FREEZING + remaining / capacity; h = 0.0; snow = 0.0; PH[PH_CONC + i] = 0.0; } else {
      let spread = volume / area;
      let flooded = max(0.0, snow - (RHOWATER - RHOICE) * spread) * RHOICE / RHOWATER;
      snow -= flooded; h = spread + flooded / RHOICE; PH[PH_CONC + i] = area;
    }
    PH[PH_SNOW + i] = snow;
  }`;

/*
 * Snow of `amount` kg/m² falling on sea cell i, on its ice or into its water.
 */
export const snowOnSea = (amount) => `if (IN[S_ICE + i] > 0.0) {
      let conc = PH[PH_CONC + i]; let cover = select(conc, 1.0, conc <= 0.0);
      PH[PH_SNOW + i] += ${amount};
      IN[S_ICE + i] += (1.0 - cover) * (${amount}) / (RHOICE * cover);
    } else { IN[S_TS + i] -= LFUS * (${amount}) / PH[PH_CAP + i]; }`;

export const PHYSICS_KERNELS = {
  physics: `${MIXED_LAYER_WGSL}@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= C) { return; }
  let pi = IN[S_PI + i]; let skin = IN[S_TS + i]; let ice = IN[S_ICE + i]; let bottom = (K - 1) * C + i;
  let wind = cellWind(i, K - 1);
  let ws = length(wind);
  D[D_WIND + i] = ws;
  let sun = vec3<f32>(P[2], P[3], P[4]);
  let mu = max(0.0, MF[F_XC + 3 * i] * sun.x + MF[F_XC + 3 * i + 1] * sun.y + MF[F_XC + 3 * i + 2] * sun.z);
  let beam = S0 * mu;
  let onLand = PH[PH_LAND + i] > 0.5; let onIceSheet = PH[PH_LAND + i] > 1.5;
  let soil0 = PH[PH_SOIL + i]; let snow0 = PH[PH_SNOW + i]; let veg0 = PH[PH_VEG + i]; let surf0 = PH[PH_SURF + i];
  let bucket = select(BUCKET, ROOTCAP, VEGETATED);
  let bareAlbedo = select(ALB_LAND, ALB_BARE + (ALB_VEG - ALB_BARE) * veg0, VEGETATED);
  let landAlbedo = select(bareAlbedo + min(1.0, snow0 / FULLSNOW) * (ALB_SNOW - bareAlbedo), ALB_ICESHEET, onIceSheet);
  let conc0 = PH[PH_CONC + i];
  let cover = select(0.0, select(conc0, 1.0, conc0 <= 0.0), ice > 0.0);
  let waterDir = openWaterAlbedo(mu);
  let iceDif = surfaceAlbedo(ice, ALB_DIF_WATER, snow0); let iceDir = surfaceAlbedo(ice, waterDir, snow0);
  let adif = select(cover * iceDif + (1.0 - cover) * ALB_DIF_WATER, landAlbedo, onLand);
  let adir = select(cover * iceDir + (1.0 - cover) * waterDir, landAlbedo, onLand);
  let ts = select(skin, cover * skin + (1.0 - cover) * FREEZING, !onLand && ice > 0.0 && cover < 1.0);
  let warmth = clamp((ts - GROWCOLD) / (GROWWARM - GROWCOLD), 0.0, 1.0);
  let aero = select(CEX, PH[PH_DRAG + i], LANDED) * max(ws, GUST);
  let roots = min(1.0, soil0 / (WETT * bucket));
  let bareWet = (1.0 - veg0) * min(1.0, surf0 / SURFCAP);
  let canopyWet = veg0 * roots / (1.0 + RSTOM * aero / max(0.05, warmth));
  let landWet = select(roots, bareWet + canopyWet, VEGETATED);
  let wetness = select(1.0, select(landWet, 1.0, snow0 > 0.0), onLand);
  let bareShare = select(0.0, bareWet / max(1e-12, bareWet + canopyWet), VEGETATED && snow0 <= 0.0);
  let ozoneHeating = beam * OZONE_ABS;
  let surfaceEmission = STEFAN * ts * ts * ts * ts;
  var vaporE: array<f32, K>; var mixedE: array<f32, K>; var cloudE: array<f32, K>; var temperature: array<f32, K>; var netFlux: array<f32, K>;
  var cloudPath = 0.0; var columnCover = 0.0;
  let tau0 = PH[PH_TAU + i];
  var deck = 0.0; var fraction = 0.0; var mlmCover = 0.0; var mlmWater = 0.0; var mlmEntrainment = 0.0; var mlmTop = 0.0;
  let mixedDepth = PH[PH_DEPTH + i] - (D[D_GEO + bottom] + LV[L_GABS + K - 1]) / GRAV;
  if (STRATUS && MLM_DECK && !onLand && 1.0 - cover > 0.0 && mixedDepth > 0.0) {
    let airT = IN[S_TH + bottom] * D[D_EXM + bottom];
    let rho = pi * LV[L_SM + K - 1] / (RGAS * airT);
    let exchange = rho * select(CEX, PH[PH_DRAG + i], LANDED) * max(ws, GUST);
    let sensible = exchange * CP * (ts - airT);
    let evap = wetness * max(0.0, exchange * (qsat(ts, pi) - IN[S_Q + bottom]));
    var deckSun = MlmSun(0.0, mu, adir, adif, 0.0, 0.0, 0.0);
    if (STRATUS_SOLAR && CLOUD_SW > 0.0) {
      var path = 0.0; var layer = 0.0;
      for (var k = 0; k < K; k++) {
        let water = max(0.0, IN[S_QC + k * C + i]) * (pi * LV[L_DS + k] / GRAV);
        path += water;
        if (k == STRATUS_K) { layer = water; }
      }
      var lit = beam - ozoneHeating;
      if (VAPOR_ABS > 0.0 && mu > 0.0) {
        let magnification = 35.0 / sqrt(1224.0 * mu * mu + 1.0);
        var vapor = 0.0;
        for (var k = 0; k < K; k++) { vapor += max(0.0, IN[S_Q + k * C + i]) * pi * LV[L_DS + k] / GRAV * sqrt(LV[L_SM + k]) * 0.1 * magnification; }
        lit -= lit * (VAPOR_ABS * 2.9 * vapor / (pow(1.0 + 141.5 * vapor, 0.635) + 5.925 * vapor));
      }
      let sky = shortwave(CLOUD_SCAT * path, cloudKeep(path), mu, adir, adif);
      deckSun = MlmSun(lit, mu, adir, adif, path, layer, select(0.0, lit * sky.w / path, path > 0.0) * layer);
    }
    let mixed = mlmColumn(i, pi, mixedDepth, sensible, evap, P[0], deckSun);
    if (mixed.ok) {
      mlmCover = mixed.cover; mlmWater = mixed.water; mlmEntrainment = mixed.entrainment; mlmTop = mixed.top;
      fraction = mixed.cover * (1.0 - cover);
      if (fraction > 0.0) { deck = min(STRATUS_MAX, mixed.water); }
    }
    if (deck <= 0.0) { fraction = 0.0; }
  }
  if (STRATUS && !MLM_DECK && !onLand && mixedDepth > 0.0 && IN[S_Q + bottom] > 0.0) {
    let upper = STABILITY_K * C + i; let lowerT = IN[S_TH + bottom] * D[D_EXM + bottom];
    let lcl = condensationLevel(lowerT, IN[S_Q + bottom], pi * LV[L_SM + K - 1]);
    let base = max(0.0, CP * (lowerT - lcl.x) / GRAV);
    let height = (D[D_GEO + upper] - D[D_GEO + bottom] + LV[L_GABS + STABILITY_K] - LV[L_GABS + K - 1]) / GRAV;
    let inversion = inversionStrength(IN[S_TH + upper] - IN[S_TH + bottom], lowerT, IN[S_TH + upper] * D[D_EXM + upper], height - base);
    let index = select(inversion, entrainmentIndex(inversion, IN[S_Q + bottom], IN[S_Q + upper]), ECTEI);
    fraction = clamp(0.19 + 0.08 * (index - 1.0), 0.0, 1.0) * clamp((ts - 278.15) / 5.0, 0.0, 1.0) * (1.0 - cover);
    if (fraction > 0.0) { deck = deckWater(lcl.x, lcl.y, mixedDepth - base); }
    if (deck <= 0.0) { fraction = 0.0; }
  }
  PH[PH_DECK + i] = deck; PH[PH_DECKF + i] = fraction;
  PH[PH_MLMCOVER + i] = mlmCover; PH[PH_MLMWATER + i] = mlmWater; PH[PH_MLMENT + i] = mlmEntrainment; PH[PH_MLMTOP + i] = mlmTop;
  for (var k = 0; k < K; k++) {
    let idx = k * C + i;
    let mass = pi * LV[L_DS + k] / GRAV;
    var eps = 1.0 - exp(-tau0 * LV[L_SHAPE + k]);
    if (COUPLED) { eps = 1.0 - exp(-VCOUP * max(0.0, IN[S_Q + idx]) * mass); }
    let water = max(0.0, IN[S_QC + idx]) * mass;
    cloudPath += water;
    var f = 1.0;
    if (PDF_COVER && water > 0.0) {
      let inside = (D[D_GEO + idx] + LV[L_GABS + k] - D[D_GEO + bottom] - LV[L_GABS + K - 1]) / GRAV < mixedDepth;
      let qsl = qsat(IN[S_TH + idx] * D[D_EXM + idx], pi * LV[L_SM + k]); let width = (1.0 - select(RHC, RHC_BL, inside)) * qsl;
      f = clamp((max(0.0, IN[S_Q + idx]) + max(0.0, IN[S_QC + idx]) - qsl + width) / (2.0 * width), COVER_FLOOR, 1.0);
    }
    if (PDF_COVER) { columnCover = max(columnCover, f * (1.0 - exp(-water / VISIBLE_PATH))); }
    cloudE[k] = select(0.0, f * (1.0 - exp(-CLOUD_ABS * water / f)), water > 0.0);
    if (STRATUS && k == STRATUS_K && deck > 0.0) { cloudE[k] = fraction * (1.0 - exp(-CLOUD_ABS * (water + deck))) + (1.0 - fraction) * cloudE[k]; }
    let clear = 1.0 - cloudE[k];
    vaporE[k] = 1.0 - (1.0 - eps) * clear;
    mixedE[k] = 1.0 - (1.0 - LV[L_GASE + k]) * clear;
    temperature[k] = IN[S_TH + idx] * D[D_EXM + idx];
    netFlux[k] = ozoneHeating * LV[L_OZ + k];
  }
  var incident = beam - ozoneHeating;
  var vaporHeating = 0.0;
  if (VAPOR_ABS > 0.0 && mu > 0.0) {
    let magnification = 35.0 / sqrt(1224.0 * mu * mu + 1.0);
    var path = 0.0; var taken = 0.0;
    for (var k = 0; k < K; k++) {
      path += max(0.0, IN[S_Q + k * C + i]) * pi * LV[L_DS + k] / GRAV * sqrt(LV[L_SM + k]) * 0.1 * magnification;
      let through = VAPOR_ABS * 2.9 * path / (pow(1.0 + 141.5 * path, 0.635) + 5.925 * path);
      netFlux[k] += incident * (through - taken);
      taken = through;
    }
    vaporHeating = incident * taken;
    incident -= vaporHeating;
  }
  if (!(columnCover > 0.0)) { columnCover = 1.0; }
  var sw = shortwave(CLOUD_SCAT * cloudPath / columnCover, cloudKeep(cloudPath / columnCover), mu, adir, adif);
  if (PDF_COVER && columnCover < 1.0) { sw = columnCover * sw + (1.0 - columnCover) * shortwave(0.0, 1.0, mu, adir, adif); }
  var clearShare = select(0.0, incident * sw.w / cloudPath, cloudPath > 0.0);
  var deckShare = 0.0;
  if (deck > 0.0) {
    let decked = shortwave(CLOUD_SCAT * (cloudPath + deck), cloudKeep(cloudPath + deck), mu, adir, adif);
    deckShare = fraction * incident * decked.w / (cloudPath + deck);
    clearShare *= 1.0 - fraction;
    sw = fraction * decked + (1.0 - fraction) * sw;
  }
  var cloudHeating = 0.0;
  if (CLOUD_SW > 0.0 && incident > 0.0 && sw.w > 0.0) {
    for (var k = 0; k < K; k++) {
      let water = max(0.0, IN[S_QC + k * C + i]) * (pi * LV[L_DS + k] / GRAV);
      let share = clearShare * water + deckShare * (water + select(0.0, deck, k == STRATUS_K));
      netFlux[k] += share;
      cloudHeating += share;
    }
  }
  let absorbed = incident * sw.x;
  let v = band(VAPOR_FRAC, &vaporE, &temperature, &netFlux, surfaceEmission);
  let g = band(GAS_FRAC, &mixedE, &temperature, &netFlux, surfaceEmission);
  let w = band(WINDOW, &cloudE, &temperature, &netFlux, surfaceEmission);
  let outgoing = v.x + g.x + w.x; let back = v.y + g.y + w.y;
  let airT = temperature[K - 1];
  let rho = pi * LV[L_SM + K - 1] / (RGAS * airT);
  let exchange = rho * select(CEX, PH[PH_DRAG + i], LANDED) * max(ws, GUST);
  let sensible = exchange * CP * (ts - airT);
  let evap = wetness * max(0.0, exchange * (qsat(ts, pi) - IN[S_Q + bottom]));
  netFlux[K - 1] += sensible;
  let net = absorbed - surfaceEmission + back - sensible - LHEAT * evap;
  let dt = P[0];
  for (var k = 0; k < K; k++) {
    let idx = k * C + i;
    let mass = pi * LV[L_DS + k] / GRAV;
    IN[S_TH + idx] += dt * netFlux[k] / (CP * mass) / D[D_EXM + idx];
  }
  IN[S_Q + bottom] += dt * evap * GRAV / (pi * LV[L_DS + K - 1]);
  let swdn = incident * sw.y;
  PH[PH_SWDN + i] = swdn;
  let directDown = incident * sw.z;
  let contrast = directDown * (iceDir - waterDir) + (swdn - directDown) * (iceDif - ALB_DIF_WATER);
  PH[PH_SFLUX + i] = net; PH[PH_ABS + i] = absorbed + ozoneHeating + vaporHeating + cloudHeating; PH[PH_ATMSW + i] = ozoneHeating + vaporHeating + cloudHeating; PH[PH_OLR + i] = outgoing; PH[PH_SH + i] = sensible; PH[PH_EVAP + i] = evap; PH[PH_INS + i] = beam; PH[PH_REFL + i] = incident - absorbed - cloudHeating; PH[PH_ADIF + i] = adif;
  let ocean = PH[PH_OFLUX + i]; let capacity = PH[PH_CAP + i];
  var T = skin; var h = ice;
  if (onLand) {
    var soil = soil0; var snow = snow0; var surf = surf0;
    T += dt * net / LANDC;
    var left = evap * dt;
    let fromSnow = min(snow, left);
    snow -= fromSnow; left -= fromSnow;
    T -= LFUS * fromSnow / LANDC;
    if (VEGETATED) {
      let fromSurface = min(surf, left * bareShare);
      surf -= fromSurface; left -= fromSurface;
      let seep = surf * (1.0 - exp(-dt / PERCT));
      surf -= seep; soil += seep;
    }
    soil = max(0.0, soil - left);
    if (snow > 0.0 && T > MELTING) {
      let energy = (T - MELTING) * LANDC;
      let melt = min(snow, energy / LFUS);
      snow -= melt; soil += melt;
      T = MELTING + (energy - melt * LFUS) / LANDC;
    }
    var cap = bucket;
    if (VEGETATED) {
      var veg = veg0 * exp(-dt / VEG_SNOW);
      if (snow <= 0.0) {
        let goal = clamp((min(soil, cap) / cap - VEG_DRY) / (VEG_WET - VEG_DRY), 0.0, 1.0);
        veg = veg0 + (goal - veg0) * (1.0 - exp(-dt * select(1.0, warmth, goal > veg0) / select(VEG_DECLINE, VEG_GROW, goal > veg0)));
      }
      if (onIceSheet) { veg = 0.0; }
      PH[PH_VEG + i] = veg;
      cap = ROOTCAP;
    }
    if (soil > cap) { PH[PH_RUNOFF + i] += soil - cap; soil = cap; }
    PH[PH_SOIL + i] = soil; PH[PH_SNOW + i] = snow; PH[PH_SURF + i] = surf;
  } else ${SEA_SURFACE_WGSL}
  IN[S_TS + i] = T; IN[S_ICE + i] = h;
}`,
  pblDiagnose: `@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= C) { return; }
  diagnoseColumn(i);
  let pi = IN[S_PI + i]; let base = (K - 1) * C + i;
  let bottomWind = cellWind(i, K - 1);
  let speed = length(bottomWind);
  let friction = sqrt(PH[PH_DRAG + i]) * max(speed, GUST);
  let zb = (D[D_GEO + base] + LV[L_GABS + K - 1]) / GRAV;
  var found = false; var riPrev = 0.0; var zPrev = zb; var depth = zb;
  for (var k = K - 2; k >= KTOP; k--) {
    if (found) { continue; }
    let idx = k * C + i;
    let z = (D[D_GEO + idx] + LV[L_GABS + k]) / GRAV;
    let dw = cellWind(i, k) - bottomWind;
    let shear = dot(dw, dw) + 100.0 * friction * friction;
    let ri = GRAV * (D[D_THV + idx] - D[D_THV + base]) * (z - zb) / (D[D_THV + base] * shear);
    if (ri > RIC) { depth = zPrev + (z - zPrev) * (RIC - riPrev) / (ri - riPrev); found = true; }
    else { riPrev = ri; zPrev = z; if (k == KTOP) { depth = z; } }
  }
  PH[PH_DEPTH + i] = depth;
  var top = depth;
  if (MLM_PROGNOSTIC && PH[PH_MLMTOP + i] > 0.0) { top = max(depth, PH[PH_MLMTOP + i]); }
  let h = top - zb;
  for (var k = KTOP; k < K; k++) { PH[PH_MIX + k * C + i] = 0.0; }
  if (h <= 0.0) { return; }
  var scale = friction;
  if (STABILITY) {
    let moisture = select(0.61 * IN[S_TH + base] * (qsat(IN[S_TS + i], pi) - IN[S_Q + base]), 0.0, PH[PH_LAND + i] > 0.5);
    let buoyancy = GRAV / IN[S_TH + base] * PH[PH_DRAG + i] * max(speed, GUST) * (IN[S_TS + i] * pow(LV[L_SM + K - 1], KAPPA) / D[D_EXM + base] - IN[S_TH + base] + moisture);
    if (buoyancy > 0.0) { scale = friction * pow(1.0 - 15.0 * max(-2.0, -0.1 * h * KARMAN * buoyancy / (friction * friction * friction)), 0.25); }
  }
  var entrainK = -1;
  for (var k = KTOP; k < K - 1; k++) {
    let idx = k * C + i; let below = idx + C;
    let zAbove = (D[D_GEO + idx] + LV[L_GABS + k]) / GRAV; let zBelow = (D[D_GEO + below] + LV[L_GABS + k + 1]) / GRAV;
    let z = 0.5 * (zAbove + zBelow) - zb;
    if (z >= h) { entrainK = k; continue; }
    let diffusivity = KARMAN * scale * z * (1.0 - z / h) * (1.0 - z / h);
    let rhoAbove = pi * LV[L_SM + k] / (RGAS * IN[S_TH + idx] * D[D_EXM + idx]);
    let rhoBelow = pi * LV[L_SM + k + 1] / (RGAS * IN[S_TH + below] * D[D_EXM + below]);
    PH[PH_MIX + idx] = 0.5 * (rhoAbove + rhoBelow) * diffusivity / (zAbove - zBelow);
  }
  if (BL_WE > 0.0 && entrainK >= KTOP) {
    let idx = entrainK * C + i; let below = idx + C;
    PH[PH_MIX + idx] = 0.5 * (pi * LV[L_SM + entrainK] / (RGAS * IN[S_TH + idx] * D[D_EXM + idx]) + pi * LV[L_SM + entrainK + 1] / (RGAS * IN[S_TH + below] * D[D_EXM + below])) * BL_WE;
  }
}`,
  adjust: `fn upperInterface(i: i32, k: i32) -> f32 {
  let idx = k * C + i;
  return (D[D_GEO + idx] + LV[L_GABS + k] + CP * D[D_THV + idx] * (D[D_EXM + idx] - D[D_EXL + idx - C])) / GRAV;
}
fn mixField(fieldOff: i32, i: i32, pi: f32, dt: f32) {
  var upper: array<f32, K>; var lower: array<f32, K>; var rhs: array<f32, K>;
  let n = K - KTOP;
  for (var j = 0; j < n; j++) {
    let k = KTOP + j;
    let mass = pi * LV[L_DS + k] / GRAV;
    upper[j] = select(0.0, dt * PH[PH_MIX + (k - 1) * C + i] / mass, j > 0);
    lower[j] = select(0.0, dt * PH[PH_MIX + k * C + i] / mass, k < K - 1);
    rhs[j] = IN[fieldOff + k * C + i];
  }
  thomas(n, &upper, &lower, &rhs);
  for (var j = 0; j < n; j++) { IN[fieldOff + (KTOP + j) * C + i] = rhs[j]; }
}
@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= C) { return; }
  let pi = IN[S_PI + i]; let dt = P[0]; let bottom = K - 1;
  var mixes = false;
  for (var k = KTOP; k < K - 1; k++) { if (PH[PH_MIX + k * C + i] > 0.0) { mixes = true; } }
  if (mixes) { mixField(S_TH, i, pi, dt); mixField(S_Q, i, pi, dt); mixField(S_QC, i, pi, dt); }
  diagnoseColumn(i);
  // saturation adjustment
  for (var k = 0; k < K; k++) {
    let idx = k * C + i;
    let ex = D[D_EXM + idx];
    let temperature = IN[S_TH + idx] * ex;
    let pressure = pi * LV[L_SM + k];
    let qs = qsat(temperature, pressure);
    let slope = qs * LHEAT / (RVAP * temperature * temperature);
    var change = (IN[S_Q + idx] - qs) / (1.0 + LHEAT * slope / CP);
    if (change < 0.0) { change = max(change, -IN[S_QC + idx]); }
    if (change == 0.0) { continue; }
    IN[S_Q + idx] -= change; IN[S_QC + idx] += change; IN[S_TH + idx] += LHEAT * change / (CP * ex);
  }
  // convection: the triggered, entraining Betts–Miller of moist.module.js
  var T: array<f32, K>; var p: array<f32, K>; var dp: array<f32, K>; var z: array<f32, K>; var Tref: array<f32, K>; var qref: array<f32, K>;
  for (var k = 0; k < K; k++) { let idx = k * C + i; T[k] = IN[S_TH + idx] * D[D_EXM + idx]; p[k] = pi * LV[L_SM + k]; dp[k] = pi * LV[L_DS + k]; z[k] = (D[D_GEO + idx] + LV[L_GABS + k]) / GRAV; }
  var produced = 0.0; var cloudBase = K; var downdraft = 0.0;
  let decked = DECK_VETO && PH[PH_MLMGATE + i] > 0.5;
  var top = -1; var base = -1; var cape = 0.0; var inhibition = 0.0; var lclP = 0.0; var thetaP = 0.0; var qP = 0.0; var energy0 = 0.0;
  if (!decked) {
    let depthTop = select(-1e30, PH[PH_DEPTH + i], BL_PARCEL);
    var weight = 0.0; var heat = 0.0; var water = 0.0; var source = bottom;
    for (var k = bottom; k >= 0; k--) {
      if (k < bottom && !(upperInterface(i, k + 1) < depthTop) && !(p[k] >= pi - PARCEL_DEPTH)) { break; }
      let idx = k * C + i;
      weight += LV[L_DS + k]; heat += LV[L_DS + k] * IN[S_TH + idx]; water += LV[L_DS + k] * max(0.0, IN[S_Q + idx]);
      source = k;
    }
    thetaP = heat / weight; qP = water / weight;
    let T0 = thetaP * D[D_EXM + bottom * C + i]; let pb = p[bottom];
    if (qP > 0.0) {
      let e = qP * pb / (EPSILON + (1.0 - EPSILON) * qP);
      let y = log(e / 611.2);
      let dewPoint = (273.15 * 17.67 - 29.65 * y) / (17.67 - y);
      var lclT = T0; lclP = pb;
      if (dewPoint < T0) { lclT = 1.0 / (1.0 / (dewPoint - 56.0) + log(T0 / dewPoint) / 800.0) + 56.0; lclP = pb * pow(lclT / T0, 1.0 / KAPPA); }
      if (lclP >= p[0]) {
        base = 0;
        for (var k = bottom; k >= 0; k--) { if (pi * LV[L_SU + k] < lclP) { base = k; break; } }
        var height = z[bottom] + CP * (T0 - lclT) / GRAV;
        var energy = CP * lclT + GRAV * height + LHEAT * qP;
        energy0 = energy;
        var temperature = lclT; var free = false;
        for (var k = bottom; k >= 0; k--) {
          let idx = k * C + i; let air = max(0.0, IN[S_Q + idx]);
          let saturated = p[k] < lclP;
          var vapour = qP;
          if (!saturated) { Tref[k] = thetaP * D[D_EXM + idx]; }
          else {
            let environment = CP * T[k] + GRAV * z[k] + LHEAT * air;
            energy = environment + (energy - environment) * exp(-ENTRAIN * (z[k] - height));
            height = z[k];
            temperature = saturatedTemperature(energy - GRAV * z[k], p[k], temperature, 1.0);
            Tref[k] = temperature;
            vapour = qsat(temperature, p[k]);
          }
          if (k >= source) { continue; }
          let work = RGAS * (Tref[k] * (1.0 + PARCEL_VIRT * vapour) - T[k] * (1.0 + PARCEL_VIRT * air)) * dp[k] / p[k];
          if (!free && saturated && work > 0.0) { free = true; }
          if (!free) { if (work < 0.0) { inhibition -= work; } }
          else if (work > 0.0) { cape += work; top = k; }
          if (saturated && T[k] - Tref[k] > 10.0) { break; }
        }
      }
    }
  }
  var passed = 0.0;
  if (top >= 0) { passed = clamp(0.5 + (cape - CAPE_MIN) / max(1.0, CAPE_MIN), 0.0, 1.0) * clamp(0.5 + (CIN_MAX - inhibition) / max(1.0, CIN_MAX), 0.0, 1.0); }
  var activity = passed;
  if (ACT_MEM > 0.0) { activity = PH[PH_CONVACT + i] + (passed - PH[PH_CONVACT + i]) * relaxedFraction(dt / ACT_MEM); }
  PH[PH_CONVACT + i] = activity;
  let eager = EAGER && top > 0 && p[top] > EAGER_TOP && cape > EAGER_CAPE && inhibition < EAGER_CIN;
  if (top >= 0 && (eager || activity > 0.5 || (activity == 0.5 && passed > 0.5))) {
    let parcelBase = base;
    if (FROM_SURFACE) { base = bottom; }
    let shallow = top > 0 && p[top] > SHALLOW_TOP;
    if (shallow) {
      let above = top - 1; let aboveQ = max(0.0, IN[S_Q + above * C + i]);
      let aboveEnergy = CP * T[above] + GRAV * z[above] + LHEAT * aboveQ;
      let span = lclP - p[above];
      for (var k = top; k <= base; k++) {
        let chi = clamp((lclP - p[k]) / span, 0.0, 1.0);
        let energy = energy0 + chi * (aboveEnergy - energy0) - GRAV * z[k];
        let water = qP + chi * (aboveQ - qP);
        var t = (energy - LHEAT * water) / CP;
        if (water > SHALLOW_RH * qsat(t, p[k])) { t = saturatedTemperature(energy, p[k], t, SHALLOW_RH); qref[k] = SHALLOW_RH * qsat(t, p[k]); }
        else { qref[k] = water; }
        Tref[k] = t;
      }
    } else {
      for (var k = top; k <= base; k++) { qref[k] = RH_REF * qsat(Tref[k], p[k]); }
    }
    var heating = 0.0; var drying = 0.0; var depth = 0.0;
    for (var k = top; k <= base; k++) { heating += CP * (Tref[k] - T[k]) * dp[k]; drying -= (qref[k] - IN[S_Q + k * C + i]) * dp[k]; depth += dp[k]; }
    if (shallow || heating > 0.0) {
      let rate = dt / RELAX;
      var rain = 0.0;
      if (!shallow && drying > 0.0) {
        let shift = (LHEAT * drying - heating) / (CP * depth);
        for (var k = top; k <= base; k++) { Tref[k] += shift; }
        rain = drying / GRAV * rate;
      } else {
        let shiftQ = drying / depth; let shiftT = -heating / (CP * depth);
        for (var k = top; k <= base; k++) { qref[k] += shiftQ; Tref[k] += shiftT; }
      }
      for (var k = top; k <= base; k++) {
        let idx = k * C + i;
        IN[S_TH + idx] += (Tref[k] - T[k]) * rate / D[D_EXM + idx];
        IN[S_Q + idx] += (qref[k] - IN[S_Q + idx]) * rate;
      }
      if (rain > 0.0) {
        if (DETRAIN > 0.0) {
          var anvilMass = 0.0; var anvilBottom = top;
          for (var k = top; k <= base; k++) { if (k == top || anvilMass < ANVIL) { anvilMass += dp[k]; anvilBottom = k; } else { break; } }
          let detrained = DETRAIN * rain;
          for (var k = top; k <= anvilBottom; k++) { IN[S_QC + k * C + i] += detrained * GRAV / anvilMass; }
          rain -= detrained;
        }
        produced = rain; cloudBase = parcelBase; downdraft = DOWNDRAFT * rain;
      }
    }
  }
  // autoconversion, and the rain and the downdraft's share evaporating as they fall
  var rained = 0.0; var left = downdraft;
  let floor = PH[PH_DEPTH + i];
  var subcloud = 0.0;
  if (DRAFT_MASS) { for (var k = cloudBase + 1; k < K; k++) { subcloud += LV[L_DS + k]; } }
  for (var k = 0; k < K; k++) {
    let idx = k * C + i;
    let mass = pi * LV[L_DS + k] / GRAV;
    var offer = left;
    if (DRAFT_MASS && k > cloudBase) { offer = min(left, downdraft * LV[L_DS + k] / subcloud); }
    let open = rained > 0.0 && (EVAP_IN_CLOUD || !(IN[S_QC + idx] > CLEAR_AIR)); let draft = offer > 0.0 && k > cloudBase;
    if ((open || draft) && RAIN_EVAP > 0.0) {
      let ex = D[D_EXM + idx];
      let temperature = IN[S_TH + idx] * ex;
      let qs = qsat(temperature, pi * LV[L_SM + k]);
      let slope = qs * LHEAT / (RVAP * temperature * temperature);
      let deficit = max(0.0, (qs - IN[S_Q + idx]) / (1.0 + LHEAT * slope / CP)) * mass;
      let available = select(0.0, rained, open) + select(0.0, offer, draft);
      let evaporated = min(available, RAIN_EVAP * deficit);
      if (evaporated > 0.0) {
        var fromDraft = 0.0; var fromRain = 0.0;
        if (evaporated >= available) { fromDraft = select(0.0, offer, draft); fromRain = select(0.0, rained, open); }
        else { fromDraft = select(0.0, evaporated * offer / available, draft); fromRain = evaporated - fromDraft; }
        left -= fromDraft;
        rained = max(0.0, rained - fromRain);
        IN[S_Q + idx] += (fromRain + fromDraft) / mass;
        IN[S_TH + idx] -= LHEAT * (fromRain + fromDraft) / (mass * CP * ex);
      }
    }
    let qc = IN[S_QC + idx];
    if (!(qc > 0.0)) { continue; }
    if (AUTO_NONE) { } else if (AUTO_BL) { if (k > 0 && upperInterface(i, k) < floor) { continue; } } else if (k >= K - 2) { continue; }
    let excess = max(0.0, qc - AUTO_T);
    let converted = min(qc, excess * (1.0 - exp(-AUTO_R * dt)) + qc * (1.0 - exp(-dt / CLOUD_LIFE)));
    IN[S_QC + idx] = qc - converted;
    rained += mass * converted;
  }
  let convected = produced - (downdraft - left);
  // filler
  for (var f = 0; f < 2; f++) {
    let off = select(S_Q, S_QC, f == 1);
    for (var k = 0; k < K - 1; k++) {
      let idx = k * C + i;
      if (IN[off + idx] < 0.0) { IN[off + idx + C] += IN[off + idx] * LV[L_DS + k] / LV[L_DS + k + 1]; IN[off + idx] = 0.0; }
    }
    if (IN[off + bottom * C + i] < 0.0) { IN[off + bottom * C + i] = 0.0; }
  }
  PH[PH_RAIN + i] += rained + convected; PH[PH_COND + i] += rained; PH[PH_CONV + i] += convected; PH[PH_STEPRAIN + i] = rained + convected;
  let airT = IN[S_TH + bottom * C + i] * D[D_EXM + bottom * C + i];
  if (PH[PH_LAND + i] < 0.5 && airT < MELTING && rained + convected > 0.0) {
    ${snowOnSea('rained + convected')}
    IN[S_TH + bottom * C + i] += LFUS * (rained + convected) * GRAV / (CP * pi * LV[L_DS + K - 1] * D[D_EXM + bottom * C + i]);
  }
  if (PH[PH_LAND + i] > 0.5) {
    if (airT < MELTING) {
      PH[PH_SNOW + i] += rained + convected;
      IN[S_TH + bottom * C + i] += LFUS * (rained + convected) * GRAV / (CP * pi * LV[L_DS + K - 1] * D[D_EXM + bottom * C + i]);
    }
    else {
      let cap = select(BUCKET, ROOTCAP, VEGETATED);
      var soil = PH[PH_SOIL + i];
      if (VEGETATED) {
        var surf = PH[PH_SURF + i] + rained + convected;
        if (surf > SURFCAP) {
          let infiltration = surf - SURFCAP; surf = SURFCAP;
          let shed = infiltration * pow(min(1.0, soil / cap), 4.0);
          PH[PH_RUNOFF + i] += shed; soil += infiltration - shed;
        }
        PH[PH_SURF + i] = surf;
      } else { soil += rained + convected; }
      if (soil > cap) { PH[PH_RUNOFF + i] += soil - cap; soil = cap; }
      PH[PH_SOIL + i] = soil;
    }
  }
  // dry convective adjustment
  var blockTop: array<i32, K>; var blockHeat: array<f32, K>; var blockWeight: array<f32, K>; var blockMass: array<f32, K>; var blockQ: array<f32, K>; var blockQc: array<f32, K>;
  var nb = 0; var merges = 0;
  for (var k = K - 1; k >= 0; k--) {
    let idx = k * C + i; let w = D[D_EXM + idx] * LV[L_DS + k];
    blockTop[nb] = k; blockHeat[nb] = IN[S_TH + idx] * w; blockWeight[nb] = w; blockMass[nb] = LV[L_DS + k];
    blockQ[nb] = IN[S_Q + idx] * LV[L_DS + k]; blockQc[nb] = IN[S_QC + idx] * LV[L_DS + k];
    nb++;
    loop {
      if (nb < 2) { break; }
      if (!(blockHeat[nb - 2] / blockWeight[nb - 2] > (blockHeat[nb - 1] / blockWeight[nb - 1]) * (1.0 + 1e-6))) { break; }
      blockHeat[nb - 2] += blockHeat[nb - 1]; blockWeight[nb - 2] += blockWeight[nb - 1]; blockMass[nb - 2] += blockMass[nb - 1];
      blockQ[nb - 2] += blockQ[nb - 1]; blockQc[nb - 2] += blockQc[nb - 1]; blockTop[nb - 2] = blockTop[nb - 1];
      nb--; merges++;
    }
  }
  if (merges > 0) {
    var lowest = K - 1;
    for (var b = 0; b < nb; b++) {
      let top = blockTop[b];
      if (top < lowest) {
        let mixed = blockHeat[b] / blockWeight[b]; let mq = blockQ[b] / blockMass[b]; let mc = blockQc[b] / blockMass[b];
        for (var k = top; k <= lowest; k++) { let idx = k * C + i; IN[S_TH + idx] = mixed; IN[S_Q + idx] = mq; IN[S_QC + idx] = mc; }
      }
      lowest = top - 1;
    }
  }
}`,
  mixMomentum: `@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let e = i32(id.x); if (e >= E) { return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  var mixes = false;
  for (var k = KTOP; k < K - 1; k++) { if (PH[PH_MIX + k * C + a] + PH[PH_MIX + k * C + b] > 0.0) { mixes = true; } }
  if (!mixes) { return; }
  let columnMass = 0.5 * (IN[S_PI + a] + IN[S_PI + b]); let dt = P[0];
  var upper: array<f32, K>; var lower: array<f32, K>; var rhs: array<f32, K>;
  let n = K - KTOP;
  for (var j = 0; j < n; j++) {
    let k = KTOP + j;
    let mass = columnMass * LV[L_DS + k] / GRAV;
    upper[j] = select(0.0, dt * 0.5 * (PH[PH_MIX + (k - 1) * C + a] + PH[PH_MIX + (k - 1) * C + b]) / mass, j > 0);
    lower[j] = select(0.0, dt * 0.5 * (PH[PH_MIX + k * C + a] + PH[PH_MIX + k * C + b]) / mass, k < K - 1);
    rhs[j] = IN[S_U + k * E + e];
  }
  var before: array<f32, K>;
  for (var j = 0; j < n; j++) { before[j] = rhs[j]; }
  thomas(n, &upper, &lower, &rhs);
  for (var j = 0; j < n; j++) { IN[S_U + (KTOP + j) * E + e] = rhs[j]; }
  var share: array<f32, K>;
  var loss = 0.0; var total = 0.0;
  for (var j = 0; j < n; j++) {
    let mass = columnMass * LV[L_DS + KTOP + j] / GRAV; let change = rhs[j] - before[j];
    loss += mass * (before[j] * before[j] - rhs[j] * rhs[j]);
    share[j] = mass * change * change;
  }
  for (var j = 0; j < n - 1; j++) {
    let k = KTOP + j; let shear = rhs[j] - rhs[j + 1];
    let part = dt * 0.5 * (PH[PH_MIX + k * C + a] + PH[PH_MIX + k * C + b]) * shear * shear;
    share[j] += part; share[j + 1] += part;
  }
  for (var j = 0; j < n; j++) { total += share[j]; }
  if (total <= 0.0) { return; }
  for (var j = 0; j < n; j++) {
    let k = KTOP + j;
    D[D_DISS + k * E + e] += loss * share[j] / (total * columnMass * LV[L_DS + k] / GRAV);
  }
}`,
};
