import { MINIMUM_CONCENTRATION, MINIMUM_VOLUME, MELTING_POINT } from '../physics/ice.module.js';
import { DARKENING_WETNESS, TRACE_SNOW, LLOYD_TAYLOR, MIAMI, startPlaceholders } from '../physics/land.module.js';
import { MIXED_LAYER_DEFAULTS, DYCOMS_LONGWAVE } from '../physics/mixedLayer.module.js';
import { DECK_CLOUD_LEVELS, UNDECIDED, VISIBLE_PATH, REFERENCE_PRESSURE, REFERENCE_RESISTANCE } from '../physics/radiation.module.js';
import { CLEAR_AIR, DECK_OPEN, DECK_CLOSED, CUMULUS_FLOOR, CUMULUS_TRACE, DEEP_REFERENCE, DEEP_CLOUD_DEPTH, IFS_ENTRAINMENT, TEST_PARCEL, IFS_PRECIPITATION, RETIRED_OPTIONS, FUSION_HEAT, BECHTOLD, SOURCE_EXCESS, fogConstants } from '../physics/moist.module.js';
import { ENTRAINMENT_DEFAULTS, CLOUD_TOP_DEFAULTS } from '../physics/boundaryLayer.module.js';
import { LONGWAVE_TABLE as DEFAULT_LONGWAVE_TABLE, LONGWAVE_CONSTANTS, GAS_MOLAR } from '../physics/longwave.module.js';
import { OZONE_GRID, OZONE_PROFILES } from '../physics/ozoneTable.module.js';
import { SUMMER_DAY } from '../physics/ozone.module.js';
import { OZONE_SHARES, OZONE_COEFFICIENTS, VISIBLE_VAPOR, VAPOR_COEFFICIENTS, VAPOR_WEIGHTS, OXYGEN, CO2_COEFFICIENT, STP_DEPTH, OZONE_CM_ATM, SCALING_PRESSURE, SCALING_EXPONENT } from '../physics/shortwaveGases.module.js';
import { exchangeConstants, EXCHANGE_WGSL } from './exchange.gpu.js';

/*
 * The column physics of the model as WGSL, one thread per column (or per
 * edge for momentum mixing), sharing the core's bindings and layouts:
 * the three-band gray radiation with clouds, the stratocumulus deck
 * (the EIS fit or the mixed-layer model) and the zenith/diffuse
 * surface reflection, bulk surface fluxes, the zero-layer sea ice and
 * its concentration, the boundary-layer diagnosis, and the adjustment
 * phase — boundary-layer mixing by the implicit tridiagonal solve,
 * saturation adjustment, the shallow and deep convective plume,
 * autoconversion and the rain's fall, the filler and the dry convective adjustment. Each is a line-by-line port of the JavaScript module it
 * names; the physics reads the Exner ratios the last RK4 stage left in
 * the diagnostic buffer, as the CPU does, and the adjustment
 * re-diagnoses the column first. The adjustment is three kernels
 * (ADJUST_KERNELS): the mixing and first saturation, the plumes, and the
 * rest, which takes the plume's rain from PH's CU fields. The cloud water's shortwave absorption
 * sits behind the constant CLOUD_SW, in the two-stream and in a pass of
 * its own for the mixed layer's sunlight, so that with
 * cloudSolarAbsorption 0 the kernel compiles to the purely scattering
 * one bit for bit: the shader compiler reassociates, and sharing terms
 * between the two changes the rounding.
 */
export function physicsConstants(o) {
  if (o.plumeEntrainmentLaw !== 'ifs' && o.plumeEntrainmentLaw !== 'gregory') throw new Error(`plumeEntrainmentLaw must be 'ifs' or 'gregory', not ${o.plumeEntrainmentLaw}`);
  if (o.plumeConversion !== 'sundqvist' && o.plumeConversion !== 'zhangMcFarlane') throw new Error(`plumeConversion must be 'sundqvist' or 'zhangMcFarlane', not ${o.plumeConversion}`);
  if (o.plumePhase !== 'mixed' && o.plumePhase !== 'liquid') throw new Error(`plumePhase must be 'mixed' or 'liquid', not ${o.plumePhase}`);
  if (o.convectionType !== 'testParcel' && o.convectionType !== 'cloudDepth' && o.convectionType !== 'top') throw new Error(`convectionType must be 'testParcel', 'cloudDepth' or 'top', not ${o.convectionType}`);
  if (o.stratusIndex !== 'eis' && o.stratusIndex !== 'ectei') throw new Error(`stratusIndex must be 'eis' or 'ectei', not ${o.stratusIndex}`);
  const m = { ...MIXED_LAYER_DEFAULTS, cloudLevels: DECK_CLOUD_LEVELS, ...o.mixedLayer };
  if (m.closure !== 'radiative' && m.closure !== 'buoyancy') throw new Error(`closure must be 'radiative' or 'buoyancy', not ${m.closure}`);
  if (m.drizzle) throw new Error('the GPU mixed-layer deck runs without drizzle');
  if (o.cloudOverlap !== 'maximum' && o.cloudOverlap !== 'maximumRandom' && o.cloudOverlap !== 'exponentialRandom') throw new Error(`cloudOverlap must be 'maximum', 'maximumRandom' or 'exponentialRandom', not ${o.cloudOverlap}`);
  const darkening = o.soilDarkening === true ? 'surface' : o.soilDarkening;
  if (darkening !== false && !DARKENING_WETNESS[darkening]) throw new Error(`soilDarkening is 'surface', 'rootZone' or false, not ${o.soilDarkening}`);
  const wetting = o.darkeningWetness ?? DARKENING_WETNESS[darkening || 'surface'];
  if (!(o.treelineWarmth[1] > o.treelineWarmth[0])) throw new Error(`treelineWarmth must rise from its first to its second temperature, not ${o.treelineWarmth}`);
  if (!(o.forestAridity[1] > o.forestAridity[0])) throw new Error(`forestAridity must rise from its first to its second index, not ${o.forestAridity}`);
  if (!(o.overcastInversion?.[1] > o.overcastInversion?.[0])) throw new Error(`overcastInversion must rise from its first to its second EIS, not ${o.overcastInversion}`);
  if (o.deckSlab !== 'fraction' && o.deckSlab !== 'midpoint') throw new Error(`deckSlab must be 'fraction' or 'midpoint', not ${o.deckSlab}`);
  if (typeof o.cumulusMemory !== 'number' || !(o.cumulusMemory >= 0 && o.cumulusMemory < Infinity)) throw new Error(`cumulusMemory must be a time in seconds, 0 or more, not ${o.cumulusMemory}`);
  if (o.deckReference !== 'interpolate' && o.deckReference !== 'layer') throw new Error(`deckReference must be 'interpolate' or 'layer', not ${o.deckReference}`);
  if (o.deckRest !== 'depth' && o.deckRest !== 'inversion' && o.deckRest !== 'regime') throw new Error(`deckRest must be 'depth', 'inversion' or 'regime', not ${o.deckRest}`);
  if (![0, 1, 2].includes(o.subsidenceSmoothing)) throw new Error(`subsidenceSmoothing must be 0, 1 or 2, not ${o.subsidenceSmoothing}`);
  for (const retired of RETIRED_OPTIONS) if (retired in o) throw new Error(`${retired} belongs to the retired Betts–Miller convection; the plume is the only scheme`);
  if (o.cumulusSource !== 'mean' && o.cumulusSource !== 'lowest') throw new Error(`cumulusSource must be 'mean' or 'lowest', not ${o.cumulusSource}`);
  if (o.plumeClosure !== 'maximum' && o.plumeClosure !== 'separate' && o.plumeClosure !== 'cape') throw new Error(`plumeClosure must be 'maximum', 'separate' or 'cape', not ${o.plumeClosure}`);
  if (o.plumeCapeParcel !== 'plume' && o.plumeCapeParcel !== 'undilute') throw new Error(`plumeCapeParcel must be 'plume' or 'undilute', not ${o.plumeCapeParcel}`);
  if (o.plumeSource !== 'mean' && o.plumeSource !== 'lowest') throw new Error(`plumeSource must be 'mean' or 'lowest', not ${o.plumeSource}`);
  if (o.plumeConsumption !== 'all' && o.plumeConsumption !== 'buoyant') throw new Error(`plumeConsumption must be 'all' or 'buoyant', not ${o.plumeConsumption}`);
  if (o.plumeSourceDepth !== 'surface50' && o.plumeSourceDepth !== 'boundaryLayer') throw new Error(`plumeSourceDepth must be 'surface50' or 'boundaryLayer', not ${o.plumeSourceDepth}`);
  if (o.excessVelocity !== 'surfaceLayer' && o.excessVelocity !== 'convective') throw new Error(`excessVelocity must be 'surfaceLayer' or 'convective', not ${o.excessVelocity}`);
  if (o.capeClosure !== 'bechtold' && o.capeClosure !== 'threshold') throw new Error(`capeClosure must be 'bechtold' or 'threshold', not ${o.capeClosure}`);
  if (o.pcapeBoundary !== 'positive' && o.pcapeBoundary !== 'signed') throw new Error(`pcapeBoundary must be 'positive' or 'signed', not ${o.pcapeBoundary}`);
  if (o.condensation !== 'uniform' && o.condensation !== 'saturation') throw new Error(`condensation must be 'uniform' or 'saturation', not ${o.condensation}`);
  if (o.autoconversionFloor !== 'lowest' && o.autoconversionFloor !== 'boundaryLayer' && o.autoconversionFloor !== 'none') throw new Error(`autoconversionFloor must be 'lowest', 'boundaryLayer' or 'none', not ${o.autoconversionFloor}`);
  if (o.fogDroplets !== null && !(Array.isArray(o.fogDroplets) && o.fogDroplets.length === 2 && o.fogDroplets.every((n) => Number.isFinite(n) && n > 0))) throw new Error(`fogDroplets must be null or [continental, sea] droplet numbers above 0 per cm³, not ${JSON.stringify(o.fogDroplets)}`);
  if (!(Number.isFinite(o.fogDeposition) && o.fogDeposition >= 0)) throw new Error(`fogDeposition must be an efficiency of 0 or more, not ${o.fogDeposition}`);
  if (!(Number.isFinite(o.fogDepositionLimit) && o.fogDepositionLimit > 0)) throw new Error(`fogDepositionLimit must be a speed above 0 in m/s, not ${o.fogDepositionLimit}`);
  const fog = fogConstants(o.fogDroplets);
  if (o.boundaryCondensation !== 'cloudLayer' && o.boundaryCondensation !== 'uniform' && o.boundaryCondensation !== 'saturation') throw new Error(`boundaryCondensation must be 'cloudLayer', 'uniform' or 'saturation', not ${o.boundaryCondensation}`);
  const entrainment = { ...ENTRAINMENT_DEFAULTS, ...o.entrainment };
  const cloudTop = { ...CLOUD_TOP_DEFAULTS, ...o.cloudTop };
  if (o.turbulence !== 'moist' && o.turbulence !== 'dry') throw new Error(`turbulence must be 'moist' or 'dry', not ${o.turbulence}`);
  if (o.boundaryCover !== 'variance' && o.boundaryCover !== 'pdf') throw new Error(`boundaryCover must be 'variance' or 'pdf', not ${o.boundaryCover}`);
  if (o.deckRegime !== 'inversion' && o.deckRegime !== 'boundaryLayer') throw new Error(`deckRegime must be 'inversion' or 'boundaryLayer', not ${o.deckRegime}`);
  const moistTurbulence = o.turbulence === 'moist';
  if (o.longwaveScheme !== 'correlated' && o.longwaveScheme !== 'gray') throw new Error(`longwaveScheme must be 'correlated' or 'gray', not ${o.longwaveScheme}`);
  if (o.longwaveOverlap !== 'exponentialRandom' && o.longwaveOverlap !== 'random') throw new Error(`longwaveOverlap must be 'exponentialRandom' or 'random', not ${o.longwaveOverlap}`);
  if (o.solarGases !== 'clirad' && o.solarGases !== 'lacisHansen') throw new Error(`solarGases must be 'clirad' or 'lacisHansen', not ${o.solarGases}`);
  if (o.ozoneProfile) throw new Error('the GPU radiation takes its ozone from its climatology, not an ozoneProfile');
  if (o.ozone !== 'afgl' && o.ozone !== 'idealized') throw new Error(`ozone must be 'afgl' or 'idealized', not ${o.ozone}`);
  const LONGWAVE_TABLE = o.longwaveTable ?? DEFAULT_LONGWAVE_TABLE, points = LONGWAVE_TABLE.points, f32 = (values) => `array<f32, ${values.length}>(${values.map((v) => v.toPrecision(9)).join(', ')})`;
  const rayleigh = o.rayleighDepth != null ? [[1, o.rayleighDepth]] : o.rayleighBands;
  if (!(rayleigh.length >= 1 && rayleigh.length <= 3 && Math.abs(rayleigh.reduce((sum, [w]) => sum + w, 0) - 1) < 1e-9 && rayleigh.every(([w, tau]) => w > 0 && tau >= 0))) throw new Error(`rayleighBands must be one to three [weight, depth] pairs whose weights sum to 1, not ${JSON.stringify(rayleigh)}`);
  const visibleSum = (term) => rayleigh.map(([w, tau]) => `${w} * ${term(`${tau / REFERENCE_PRESSURE} * light.y + light.z`)}`).join(' + ');
  return `
fn visibleStreams(cloudDepth: f32, keep: f32, mu: f32, adir: f32, adif: f32, light: vec3<f32>) -> vec4<f32> { return ${visibleSum((depth) => `stream(cloudDepth + ${depth}, keep, mu, adir, adif)`)}; }
fn visibleEscape(cloudDepth: f32, keep: f32, mu: f32, adir: f32, adif: f32, light: vec3<f32>) -> f32 { return ${visibleSum((depth) => `streamEscape(cloudDepth + ${depth}, keep, mu, adir, adif)`)}; }
fn visibleReflection(cloudDepth: f32, keep: f32, mu: f32, light: vec3<f32>) -> f32 { return ${visibleSum((depth) => `reflection(cloudDepth + ${depth}, keep, mu)`)}; }
const MOIST_BL: bool = ${moistTurbulence}; const CT_THRESH: f32 = ${cloudTop.threshold}; const CT_HMAX: f32 = ${cloudTop.maximumHeight}; const CT_PERT: f32 = ${cloudTop.perturbation}; const CT_PROFILE: f32 = ${cloudTop.profile}; const CT_EXCESS: f32 = ${cloudTop.excess}; const CT_TOLERANCE: f32 = ${cloudTop.tolerance}; const CT_CUMULUS: f32 = ${cloudTop.cumulusDepth};
const BL_A2: f32 = ${entrainment.evaporativeEnhancement}; const BL_AMAX: f32 = ${entrainment.maximumEfficiency}; const BL_TAPER: bool = ${!!entrainment.taper}; const BL_JUMP2: bool = ${entrainment.jumpLayers > 1};
const VARIANCE_COVER: bool = ${moistTurbulence && o.boundaryCover === 'variance'}; const VAR_FLOOR: f32 = ${o.varianceFloor}; const VAR_SCALE: f32 = ${o.varianceScale}; const MIX_LENGTH: f32 = ${o.mixingLength}; const STABLE_LENGTH: f32 = ${o.stableMixingLength};
const MLM_BL_GATE: bool = ${o.deckRegime === 'boundaryLayer'}; const MLM_BYPASS: bool = ${!!o.deckBypass};
const S0: f32 = ${o.solarConstant}; const STEFAN: f32 = 5.670374419e-8; const LHEAT: f32 = ${o.latentHeat}; const EPSILON: f32 = 0.622; const RVAP: f32 = ${o.R / 0.622};
const PDF_COVER: bool = ${o.cloudCover === 'pdf'}; const VISIBLE_PATH: f32 = ${VISIBLE_PATH}; const RHC: f32 = ${o.criticalHumidity}; const RHC_BL: f32 = ${o.boundaryCriticalHumidity}; const COVER_FLOOR: f32 = ${o.coverFloor ?? 0.01}; const BOUND_WIDTH: bool = ${o.overcastWater != null}; const OVERCAST_WATER: f32 = ${o.overcastWater ?? 0}; const OVERCAST_EIS: f32 = ${o.overcastInversion[0]}; const OVERCAST_RAMP: f32 = ${o.overcastInversion[1] - o.overcastInversion[0]}; const RANDOM_OVERLAP: bool = ${o.cloudOverlap === 'maximumRandom'}; const EXP_OVERLAP: bool = ${o.cloudOverlap === 'exponentialRandom'}; const DECOR_LEN: f32 = ${o.decorrelationLength}; const DECOR_SLOPE: f32 = ${o.decorrelationSlope};
const CLOUD_ABS: f32 = ${o.cloudAbsorption ?? 0}; const CLOUD_SCAT: f32 = ${o.cloudScattering ?? 0}; const GRAY_LW: bool = ${o.cloudAbsorption != null}; const GRAY_SW: bool = ${o.cloudScattering != null};
const LIQUID_T: f32 = ${o.liquidTemperature}; const ICE_T: f32 = ${o.iceTemperature}; const DROP_SEA: f32 = ${o.seaDropletRadius}; const DROP_LAND: f32 = ${o.landDropletRadius}; const LIQUID_IR: f32 = ${o.liquidInfrared}; const IR_DIFFUSIVITY: f32 = ${o.diffusivity}; const ICE_WARMEST: f32 = ${o.iceFitWarmest}; const ICE_COLDEST: f32 = ${o.iceFitColdest};
const CLOUD_SW: f32 = ${o.cloudSolarAbsorption}; const WINDOW: f32 = ${o.window}; const GAS_FRAC: f32 = ${o.gasFraction};
const STRATUS: bool = ${!!o.stratus}; const ECTEI: bool = ${o.stratusIndex === 'ectei'}; const STRATUS_SCALE: f32 = ${o.stratusScale}; const STRATUS_MAX: f32 = ${o.stratusWaterMax}; const STRATUS_K: i32 = ${o.stratusLayer}; const STABILITY_K: i32 = ${o.stabilityLayer};
const VAPOR_FRAC: f32 = ${1 - o.window - o.gasFraction}; const OZONE_ABS: f32 = ${o.ozoneAbsorption}; const VAPOR_ABS: f32 = ${o.vaporAbsorption}; const CEX: f32 = ${o.exchangeCoefficient};
const VCOUP: f32 = ${o.vaporCoupling}; const COUPLED: bool = ${o.vaporCoupling > 0}; const SKYLIGHT: f32 = ${o.skylight}; const DIFFUSE_MU: f32 = 0.6; const CLEAR_SKY: bool = ${!!o.clearSkyPass};
const LW_CORRELATED: bool = ${o.longwaveScheme === 'correlated'}; const LW_CHAIN: bool = ${o.longwaveOverlap === 'exponentialRandom'}; const SOLAR_CLIRAD: bool = ${o.solarGases === 'clirad'}; const NG: i32 = ${points.length}; const LW_D: f32 = ${LONGWAVE_CONSTANTS.diffusivity}; const LW_PREF: f32 = ${LONGWAVE_CONSTANTS.pRef};
const LW_TSELF: f32 = ${LONGWAVE_TABLE.tSelf}; const LW_TCO2: f32 = ${LONGWAVE_TABLE.tCo2}; const LW_NO3: f32 = ${LONGWAVE_TABLE.nO3};
const LW_DH2O: f32 = ${LONGWAVE_TABLE.dopplerH2o ?? 0}; const LW_DCO2: f32 = ${LONGWAVE_TABLE.dopplerCo2 ?? 0}; const LW_DO3: f32 = ${LONGWAVE_TABLE.dopplerO3 ?? 0};
const CO2_MASS: f32 = ${o.carbonDioxide * GAS_MOLAR.co2 / GAS_MOLAR.air}; const CH4_MASS: f32 = ${o.methane * GAS_MOLAR.ch4 / GAS_MOLAR.air}; const N2O_MASS: f32 = ${o.nitrousOxide * GAS_MOLAR.n2o / GAS_MOLAR.air};
const LW_K: array<f32, ${6 * points.length}> = ${f32(points.flatMap((row) => row.slice(0, 6)))};
const LW_PLANCK: array<f32, ${5 * points.length}> = ${f32(points.flatMap((row) => row.slice(6, 11)))};
const OZ_AFGL: bool = ${o.ozone === 'afgl'}; const OZ_N: i32 = ${OZONE_GRID.points}; const OZ_LOW: f32 = ${OZONE_GRID.lowest}; const OZ_STEP: f32 = ${Math.log(OZONE_GRID.highest / OZONE_GRID.lowest) / (OZONE_GRID.points - 1)}; const OZ_SUMMER: f32 = ${SUMMER_DAY / 365};
const OZ_TABLE: array<f32, ${5 * OZONE_GRID.points}> = ${f32(OZONE_PROFILES.flat())};
const OZ_P: array<f32, ${OZONE_GRID.points}> = ${f32(Array.from({ length: OZONE_GRID.points }, (_, j) => OZONE_GRID.lowest * Math.exp(Math.log(OZONE_GRID.highest / OZONE_GRID.lowest) / (OZONE_GRID.points - 1) * j)))};
const OZ_EQ: f32 = ${o.ozoneColumn[0]}; const OZ_POLE: f32 = ${o.ozoneColumn[1]}; const OZ_KG: f32 = ${OZONE_CM_ATM}; const VIS_O3: f32 = ${(OZONE_SHARES[6] * OZONE_COEFFICIENTS[6] + OZONE_SHARES[7] * OZONE_COEFFICIENTS[7]) / (OZONE_SHARES[6] + OZONE_SHARES[7])};
const O3_SHARE: array<f32, 8> = ${f32(OZONE_SHARES)}; const O3_COEF: array<f32, 8> = ${f32(OZONE_COEFFICIENTS)};
const H2O_K: array<f32, 10> = ${f32(VAPOR_COEFFICIENTS)}; const H2O_W: array<f32, 10> = ${f32(VAPOR_WEIGHTS)};
const VIS_H2O_S: f32 = ${VISIBLE_VAPOR.share}; const VIS_H2O_K: f32 = ${VISIBLE_VAPOR.coefficient}; const VAP_STRENGTH: f32 = ${o.vaporStrength};
const O2_SHARE: f32 = ${OXYGEN.share}; const O2_K: f32 = ${OXYGEN.coefficient}; const O2_PATH: f32 = ${OXYGEN.mixingRatio * STP_DEPTH}; const CO2_PATH: f32 = ${o.carbonDioxide * STP_DEPTH}; const CO2_SW_K: f32 = ${CO2_COEFFICIENT};
const SCALE_P: f32 = ${SCALING_PRESSURE}; const SCALE_N: f32 = ${SCALING_EXPONENT};
const SCATTER: bool = ${rayleigh.some(([, tau]) => tau > 0) || o.nearInfraredRayleigh > 0 || o.landAerosol > 0 || o.seaAerosol > 0}; const NIR_RAY: f32 = ${o.nearInfraredRayleigh / REFERENCE_PRESSURE}; const UPWARD: bool = ${!!o.upwardAbsorption}; const DIFFUSE_PATH: f32 = ${5 / 3}; const VIS_FRAC: f32 = ${o.visibleFraction}; const LAND_AER: f32 = ${o.landAerosol}; const SEA_AER: f32 = ${o.seaAerosol}; const AER_ABS: f32 = ${1 - o.aerosolAlbedo}; const AER_SCAT: f32 = ${(1 - o.aerosolAsymmetry) * o.aerosolAlbedo};
const ALB_ICE: f32 = ${o.iceAlbedo}; const FULLALB: f32 = ${o.fullAlbedoThickness}; const ALB_DIF_WATER: f32 = ${o.diffuseWaterAlbedo}; const ALB_ICEMELT: f32 = ${o.meltingIceAlbedo}; const ICE_MELTRANGE: f32 = ${o.iceMeltingRange};
const AGEING: bool = ${!!o.snowAgeing}; const ALB_FRESH: f32 = ${o.freshSnowAlbedo}; const SNOW_COLDAGE: f32 = ${o.coldSnowAgeing}; const SNOW_MELTAGE: f32 = ${o.meltingSnowAgeing}; const SNOW_REFRESH: f32 = ${o.refreshSnowfall}; const WETSNOW: f32 = ${o.wetSnowRange}; const AGE_ACT: f32 = ${o.ageingActivation}; const ICE_SNOWFLOOR: f32 = ${o.iceSnowFloor};
const ALB_ICESNOW: f32 = ${o.iceSnowAlbedo}; const FULLSNOW_ICE: f32 = ${o.iceFullSnow}; const KSNOW: f32 = ${o.snowConductivity}; const RHOSNOW: f32 = ${o.snowDensity}; const RHOICE: f32 = ${o.iceDensity}; const RHOWATER: f32 = ${o.waterDensity};
const FREEZING: f32 = 271.35; const MELTING: f32 = 273.15; const SKINC: f32 = ${o.skinHeatCapacity}; const COND: f32 = ${o.conductivity}; const HMIN: f32 = ${o.minimumThickness}; const LATENT_ICE: f32 = ${o.iceDensity * o.latentHeatFusion};
const LEADC: f32 = ${o.leadClosing}; const LEADX: f32 = ${o.leadExchange}; const MIN_CONC: f32 = ${MINIMUM_CONCENTRATION}; const MIN_VOLUME: f32 = ${MINIMUM_VOLUME};
const AUTO_T: f32 = ${o.autoconversionThreshold}; const AUTO_R: f32 = ${o.autoconversionRate}; const CLOUD_LIFE: f32 = ${o.cloudLifetime}; const UPPER_LIFE: f32 = ${o.upperCloudLifetime ?? o.cloudLifetime}; const UPPER_SPLIT: bool = ${o.upperCloudLifetime != null}; const STRAT_LIFE: f32 = ${o.stratiformLifetime ?? o.cloudLifetime}; const STRAT_SPLIT: bool = ${o.stratiformLifetime != null};
const RAIN_EVAP: f32 = ${o.rainEvaporation};
const UNIFORM: bool = ${o.condensation === 'uniform'}; const BL_UNIFORM: bool = ${o.boundaryCondensation === 'uniform'}; const BL_CLOUDLAYER: bool = ${o.boundaryCondensation === 'cloudLayer' && moistTurbulence}; const ICE_SAT: bool = ${!!o.iceSaturation}; const NUCLEATION: bool = ${!!o.iceNucleation && !!o.iceSaturation}; const LFUSION: f32 = ${FUSION_HEAT}; const RHC_SURF: f32 = ${o.surfaceCriticalHumidity}; const RHC_TOP: f32 = ${o.topCriticalHumidity}; const RHC_EXP: f32 = ${o.criticalExponent};
const ICE_FALL: bool = ${o.iceFall != null}; const FALL_C: f32 = ${o.iceFall ?? 0}; const FALL_EXP: f32 = ${o.iceFallExponent};
const FOG: bool = ${o.fogDroplets !== null}; const FOG_SETTLE_LAND: f32 = ${fog.settleLand}; const FOG_SETTLE_SEA: f32 = ${fog.settleSea}; const FOG_KK_LAND: f32 = ${fog.drizzleLand}; const FOG_KK_SEA: f32 = ${fog.drizzleSea}; const FOG_THIRDS: f32 = 0.6666666666666666;
const FOG_DEP: bool = ${o.fogDeposition > 0 && o.surfaceExchange === 'roughness' && !!o.implicitDrag}; const FOG_E: f32 = ${o.fogDeposition}; const FOG_VMAX: f32 = ${o.fogDepositionLimit};
const AUTO_BL: bool = ${o.autoconversionFloor === 'boundaryLayer'}; const CLEAR_AIR: f32 = ${CLEAR_AIR}; const CIN_MAX: f32 = ${o.inhibitionThreshold}; const SHALLOW_TOP: f32 = ${o.shallowTop};
const DECK_VETO: bool = ${o.deckVeto !== false}; const COUPLED_VETO: bool = ${!!o.coupledVeto && o.turbulence !== 'dry'}; const EVAP_IN_CLOUD: bool = ${!!o.evaporationInCloud}; const AUTO_NONE: bool = ${o.autoconversionFloor === 'none'};
const DECK_OPEN: f32 = ${DECK_OPEN}; const DECK_CLOSED: f32 = ${DECK_CLOSED}; const PARCEL_VIRT: f32 = ${o.virtualBuoyancy === false ? 0 : 'VIRT'};
const CU_FLOOR: f32 = ${CUMULUS_FLOOR}; const CU_K0: i32 = ${o.cumulusK0 ?? 0}; const CU_C: f32 = ${o.cumulusClosure}; const CU_EPS: f32 = ${o.cumulusEntrainment}; const CU_DEL: f32 = ${o.cumulusDetrainment}; const CU_SOURCE: f32 = ${o.cumulusSourceDepth}; const CU_LOSS: f32 = ${o.cumulusBoundaryLoss};
const CU_FRIC: f32 = ${o.cumulusFriction}; const CU_OVER: f32 = ${o.cumulusOvershoot}; const CU_WU: f32 = ${o.cumulusUpdraft}; const CU_RAIN: bool = ${o.cumulusRain != null}; const CU_RAIN_Q: f32 = ${o.cumulusRain ?? 0}; const CU_LOWEST: bool = ${o.cumulusSource === 'lowest'};
const CU_LOADING: f32 = ${o.virtualBuoyancy === false ? 0 : 1}; const CU_CLOUD: bool = ${o.cumulusCloud !== false && o.cloudCover === 'pdf'}; const CU_MEMORY: f32 = ${o.cumulusMemory}; const CU_TRACE: f32 = ${CUMULUS_TRACE};
const PL_SEPARATE: bool = ${o.plumeClosure === 'separate' && o.convectionType === 'top'}; const PL_BY_DEPTH: bool = ${o.convectionType === 'cloudDepth'}; const PL_PARCEL: bool = ${o.convectionType === 'testParcel'}; const PL_MIXED: bool = ${o.plumePhase === 'mixed'}; const PL_SUNDQVIST: bool = ${o.plumeConversion === 'sundqvist'}; const SQ_C00: f32 = ${IFS_PRECIPITATION.conversion}; const SQ_LIQ: f32 = ${IFS_PRECIPITATION.liquidFactor}; const SQ_CRIT: f32 = ${IFS_PRECIPITATION.critical}; const SQ_SEA: f32 = ${IFS_PRECIPITATION.seaThreshold}; const SQ_LAND: f32 = ${IFS_PRECIPITATION.landThreshold}; const SQ_BF: f32 = ${IFS_PRECIPITATION.bergeron}; const SQ_ICE: f32 = ${IFS_PRECIPITATION.ice}; const SQ_SPEED: f32 = ${IFS_PRECIPITATION.speed}; const SQ_VSCALE: f32 = ${IFS_PRECIPITATION.velocityScale}; const TP_EPS: f32 = ${TEST_PARCEL.entrainment}; const TP_REMOVE: f32 = ${TEST_PARCEL.removal}; const PL_IFS: bool = ${o.plumeEntrainmentLaw === 'ifs'}; const IFS_EPS: f32 = ${IFS_ENTRAINMENT.entrainment}; const IFS_RH: f32 = ${IFS_ENTRAINMENT.humidity}; const IFS_DEL: f32 = ${IFS_ENTRAINMENT.detrainment}; const IFS_DRH: f32 = ${IFS_ENTRAINMENT.detrainmentHumidity}; const IFS_DRAG: f32 = ${IFS_ENTRAINMENT.drag}; const DEEP_DEPTH: f32 = ${DEEP_CLOUD_DEPTH}; const PL_RELAXED: bool = ${o.plumeClosure !== 'maximum'}; const PL_LOWEST: bool = ${o.plumeSource === 'lowest'}; const PL_UNDILUTE: bool = ${o.plumeCapeParcel === 'undilute'}; const PL_W0: f32 = ${o.plumeVelocity}; const PL_ACC: f32 = ${o.plumeAcceleration}; const PL_DRAG: f32 = ${o.plumeDrag}; const PL_EPS: f32 = ${o.plumeEntrainment}; const PL_FLOOR: f32 = ${o.plumeEntrainmentFloor}; const PL_GROWTH: f32 = ${o.plumeMassGrowth};
const PL_MOMENTUM: bool = ${!!o.plumeMomentum}; const PL_BUOYANT_F: bool = ${o.plumeConsumption === 'buoyant'}; const PL_RAIN_RATE: f32 = ${o.plumeRainRate}; const PL_RAIN_Q: f32 = ${o.plumeRainThreshold}; const PL_EVAP: f32 = ${o.plumeRainEvaporation}; const DD_SHARE: f32 = ${o.downdraftShare}; const DD_EPS: f32 = ${o.downdraftEntrainment}; const PL_CAPE: f32 = ${o.plumeCape}; const PL_TAU: f32 = ${o.plumeRelaxation}; const DEEP_REFERENCE: f32 = ${DEEP_REFERENCE};
const PL_SURFACE50: bool = ${o.plumeSourceDepth === 'surface50'}; const EX_COEF: f32 = ${SOURCE_EXCESS.coefficient}; const EX_T: f32 = ${SOURCE_EXCESS.temperature}; const EX_Q: f32 = ${SOURCE_EXCESS.humidity};
const EX_LAYER: bool = ${o.excessVelocity === 'surfaceLayer'}; const EX_SCALE: f32 = ${SOURCE_EXCESS.scale}; const EX_STAB: f32 = ${SOURCE_EXCESS.layer}; const EX_KARMAN: f32 = ${SOURCE_EXCESS.karman}; const EX_USTAR: f32 = ${SOURCE_EXCESS.friction}; const EX_VIRT: f32 = ${SOURCE_EXCESS.virtual};
const PL_BECHTOLD: bool = ${o.capeClosure === 'bechtold'}; const BT_POSITIVE: bool = ${o.pcapeBoundary === 'positive'}; const BT_SCALE: f32 = ${BECHTOLD.resolution / BECHTOLD.reference}; const BT_SHORT: f32 = ${BECHTOLD.shortest}; const BT_LONG: f32 = ${BECHTOLD.longest}; const BT_WIND: f32 = ${BECHTOLD.boundaryWind}; const BT_TSTAR: f32 = ${BECHTOLD.temperatureScale}; const SUB_K: i32 = ${o.subcloudLayers ?? 0};
const BL_ENTRAIN: bool = ${entrainment.efficiency > 0 || entrainment.shear > 0}; const BL_A: f32 = ${entrainment.efficiency}; const BL_AS: f32 = ${entrainment.shear}; const BL_WEMAX: f32 = ${entrainment.cap}; const BL_BMIN: f32 = ${entrainment.jumpFloor}; const BL_ONSET: f32 = ${entrainment.shearOnset}; const RIC: f32 = ${o.richardsonCritical}; const KARMAN: f32 = ${o.vonKarman}; const STABILITY: bool = ${o.stability ? 'true' : 'false'}; const KTOP: i32 = ${o.kTop};
const LANDED: bool = ${!!o.landed}; const LANDC: f32 = ${o.landHeatCapacity}; const BUCKET: f32 = ${o.bucketCapacity}; const WETT: f32 = ${o.wetnessThreshold}; const ALB_LAND: f32 = ${o.landAlbedo}; const VEGETATED: bool = ${!!o.vegetation}; const ALB_BARE: f32 = ${o.bareAlbedo}; const ALB_VEG: f32 = ${o.vegetatedAlbedo}; const DARKENING: bool = ${darkening !== false}; const DARK_SURFACE: bool = ${darkening === 'surface'}; const ALB_WETSOIL: f32 = ${o.wetSoilAlbedo}; const DARK_FROM: f32 = ${wetting[0]}; const DARK_SPAN: f32 = ${wetting[1] - wetting[0]}; const ROOTCAP: f32 = ${o.rootZoneCapacity};
const MLM_DECK: bool = ${!!o.mixedLayerDeck}; const STRATUS_SOLAR: bool = ${!!o.stratusSolar}; const MLM_SUBSIDENCE: f32 = ${o.stratusSubsidence}; const MLM_MININV: f32 = ${o.minimumInversion}; const MLM_CEILINV: f32 = ${o.ceilingInversion ?? o.minimumInversion}; const MLM_MEMORY: f32 = ${o.subsidenceMemory};
const MLM_LEVELS: i32 = ${m.cloudLevels}; const MLM_NODES: i32 = ${m.cloudLevels + 1}; const MLM_BUOYANCY: bool = ${m.closure === 'buoyancy'}; const MLM_DELTA: f32 = 1.0 / EPSILON - 1.0; const MLM_LC: f32 = LHEAT / CP;
const MLM_A1: f32 = ${m.entrainmentEfficiency}; const MLM_A2: f32 = ${m.evaporativeEnhancement}; const MLM_AMAX: f32 = ${m.maximumEfficiency}; const MLM_WEMAX: f32 = ${m.maximumEntrainment}; const MLM_MINJUMP: f32 = ${m.minimumJump};
const MLM_ONSET: f32 = ${m.decouplingOnset}; const MLM_DRATIO: f32 = ${m.decoupledRatio}; const MLM_DCOVER: f32 = ${m.decoupledCover}; const DYC_F0: f32 = ${DYCOMS_LONGWAVE.F0}; const DYC_F1: f32 = ${DYCOMS_LONGWAVE.F1}; const DYC_K: f32 = ${DYCOMS_LONGWAVE.kappa};
const MLM_PASSES: i32 = ${o.subsidenceSmoothing}; const MLM_FRACTION: bool = ${o.deckSlab === 'fraction'}; const MLM_INTERPOLATE: bool = ${o.deckReference === 'interpolate'}; const MLM_PROGNOSTIC: bool = ${o.prognosticHeight ? 'true' : 'false'}; const MLM_GATEMEM: f32 = ${o.gateMemory}; const MLM_UNDECIDED: f32 = ${UNDECIDED}; const MLM_HMEM: f32 = ${m.heightMemory}; const MLM_HMAX: f32 = ${m.maximumHeight}; const MLM_REST_INVERSION: bool = ${o.deckRest !== 'depth'}; const MLM_REST_REGIME: bool = ${o.deckRest === 'regime' && moistTurbulence}; const MLM_CUCEIL: f32 = ${o.cumulusCeiling};
const ALB_ICESHEET: f32 = ${o.iceSheetAlbedo}; const SURFCAP: f32 = ${o.surfaceCapacity}; const PERCT: f32 = ${o.percolationTime}; const RSTOM: f32 = ${o.stomatalResistance}; const GROWCOLD: f32 = ${o.growthColdest}; const GROWWARM: f32 = ${o.growthWarmest}; const VEG_DRY: f32 = ${o.dryWetness}; const VEG_WET: f32 = ${o.wetWetness}; const VEG_GROW: f32 = ${o.growthTime}; const VEG_DECLINE: f32 = ${o.declineTime}; const VEG_SNOW: f32 = ${o.snowDeclineTime}; const ALB_SNOW: f32 = ${o.snowAlbedo}; const FULLSNOW: f32 = ${o.fullSnow}; const LFUS: f32 = ${o.latentHeatFusion};
const ALB_OLDSNOW: f32 = ${o.oldSnowAlbedo}; const MASKED: bool = ${!!o.snowMasking && !!o.vegetation}; const ALB_FOREST: f32 = ${o.forestSnowAlbedo}; const CLOSED_CANOPY: f32 = ${o.closedCanopy}; const CANOPY_MEM: f32 = ${o.canopyMemory};
const TREELINE: bool = ${!!o.treeline}; const SEASON_C: f32 = ${o.seasonThreshold}; const SEASON_K: f32 = ${MELTING_POINT + o.seasonThreshold}; const SEASON_SHORTEST: f32 = ${o.minimumSeason / 365}; const SEASON_MEM: f32 = ${o.seasonMemory}; const TREE_LO: f32 = ${o.treelineWarmth[0]}; const TREE_SPAN: f32 = ${o.treelineWarmth[1] - o.treelineWarmth[0]}; const TREE_GROW: f32 = ${o.treeGrowthTime}; const TREE_DECLINE: f32 = ${o.treeDeclineTime};
const REF_RESIST: f32 = ${REFERENCE_RESISTANCE}; const GATED: bool = ${!!o.treeline && !!o.treeMoisture}; const MOIST_MEM: f32 = ${o.moistureMemory}; const ARID_LO: f32 = ${o.forestAridity[0]}; const ARID_SPAN: f32 = ${o.forestAridity[1] - o.forestAridity[0]};
const GRASSY: bool = ${!!o.grassland && !!o.vegetation}; const ALB_FORESTV: f32 = ${o.forestAlbedo}; const ALB_GRASS: f32 = ${o.grassAlbedo}; const GRASS_SNOW: f32 = ${o.grassSnowDarkening};
const HUMIC: bool = ${!!o.soilCarbon && !!o.vegetation}; const ALB_MINERAL: f32 = ${o.mineralAlbedo}; const ALB_HUMUS: f32 = ${o.humusAlbedo}; const HUMUS_SCALE: f32 = ${100 / (o.topsoilMass * o.organicScale)};
const LITTER_IN: f32 = ${o.litterInput}; const LITTER_TREE: f32 = ${o.treeLitter}; const LITTER_GRASS: f32 = ${o.grassLitter}; const DECAY_RATE: f32 = ${1 / o.soilTurnover}; const DECAY_WILT: f32 = ${o.decompositionWilting}; const DECAY_OPT: f32 = ${0.5 * (1 + o.decompositionWilting)}; const CARBON_ACC: f32 = ${o.carbonAcceleration};
const LT_E: f32 = ${LLOYD_TAYLOR.activation}; const LT_REF: f32 = ${1 / LLOYD_TAYLOR.reference}; const LT_T0: f32 = ${LLOYD_TAYLOR.offset}; const MIAMI_A: f32 = ${MIAMI[0]}; const MIAMI_B: f32 = ${MIAMI[1]};
const HELD_SHARE: f32 = ${startPlaceholders(o).neutral.share}; const HELD_NEUTRAL_C: f32 = ${startPlaceholders(o).neutral.carbon}; const HELD_GREEN_C: f32 = ${startPlaceholders(o).green.carbon};
${exchangeConstants(o)}`;
}

export const PHYSICS_FUNCTIONS = `
fn esat(T: f32) -> f32 { return 611.2 * exp(17.67 * (T - 273.15) / (T - 29.65)); }
fn qsat(T: f32, p: f32) -> f32 { let es = esat(T); let dry = p - (1.0 - EPSILON) * es; return select(1.0, EPSILON * es / dry, dry > 0.0); }
fn smallRate(x: f32) -> f32 { return select(1.0 - exp(-x), x * (1.0 - 0.5 * x), x < 1e-3); }
fn cloudSat(T: f32, p: f32) -> vec2<f32> {
  let alpha = select(1.0, clamp((T - ICE_T) / (LIQUID_T - ICE_T), 0.0, 1.0), ICE_SAT);
  var es = esat(T);
  if (alpha < 1.0) { es = alpha * es + (1.0 - alpha) * 611.21 * exp(22.587 * (T - 273.16) / (T + 0.7)); }
  let dry = p - (1.0 - EPSILON) * es;
  let qs = select(1.0, EPSILON * es / dry, dry > 0.0);
  return vec2<f32>(qs, qs * (LHEAT + (1.0 - alpha) * LFUSION) / (RVAP * T * T));
}
fn criticalProfile(p: f32, ps: f32) -> f32 { return RHC_TOP + (RHC_SURF - RHC_TOP) * exp(1.0 - pow(ps / p, RHC_EXP)); }
fn uniformWidth(s: vec2<f32>, p: f32, ps: f32) -> f32 {
  let rhc = criticalProfile(p, ps);
  return (1.0 - rhc) * s.x / (1.0 + LHEAT * s.y / CP);
}
fn plumeSeen(seen: f32, cu: vec2<f32>, mass: f32) -> f32 {
  if (!PDF_COVER || !(cu.y * mass > 0.0)) { return seen; }
  return max(seen, cu.x * smallRate(cu.y * mass / (cu.x * VISIBLE_PATH)));
}
fn uniformCover(qc: f32, b: f32) -> f32 {
  if (qc >= b) { return 1.0; }
  if (qc > 0.0) { return sqrt(qc / b); }
  return 0.0;
}
fn openWaterAlbedo(mu: f32) -> f32 { return 0.026 / (pow(mu, 1.7) + 0.065) + 0.15 * (mu - 0.1) * (mu - 0.5) * (mu - 1.0); }
fn surfaceAlbedo(h: f32, water: f32, snow: f32, T: f32, snowAlbedo: f32) -> f32 {
  if (h <= 0.0) { return water; }
  let bareIce = ALB_ICE + (ALB_ICEMELT - ALB_ICE) * clamp((T - MELTING + ICE_MELTRANGE) / ICE_MELTRANGE, 0.0, 1.0);
  let bare = water + (bareIce - water) * min(1.0, h / FULLALB);
  return bare + (select(ALB_ICESNOW, snowAlbedo, AGEING) - bare) * min(1.0, snow / FULLSNOW_ICE);
}
fn agedSnow(albedo: f32, T: f32, dt: f32, floor: f32) -> f32 {
  let days = dt / 86400.0;
  if (T >= MELTING - WETSNOW) { return floor + (albedo - floor) * exp(-SNOW_MELTAGE * days); }
  var pace = 1.0;
  if (AGE_ACT > 0.0) { pace = min(1.0, exp(AGE_ACT * (T - MELTING) / (MELTING * T))); }
  return max(floor, albedo - SNOW_COLDAGE * pace * days);
}
fn refreshedSnow(albedo: f32, fall: f32) -> f32 { return albedo + min(1.0, fall / SNOW_REFRESH) * (ALB_FRESH - albedo); }
fn cellWind(i: i32, k: i32) -> vec3<f32> {
  var w = vec3<f32>(0.0, 0.0, 0.0);
  for (var m = 0; m < MAXE; m++) {
    let e = MI[EOC + MAXE * i + m];
    let s = abs(f32(MI[ESC + MAXE * i + m])) * 0.5 * MF[F_DC + e] * MF[F_DV + e] * IN[S_U + k * E + e];
    w += s * vec3<f32>(MF[F_NEDGE + 3 * e], MF[F_NEDGE + 3 * e + 1], MF[F_NEDGE + 3 * e + 2]);
  }
  return w / MF[F_AREA + i];
}
fn edgesOf(i: i32) -> array<i32, MAXE> {
  var e: array<i32, MAXE>;
  for (var m = 0; m < MAXE; m++) { e[m] = MI[EOC + MAXE * i + m]; }
  return e;
}
fn edgeWind(edges: ptr<function, array<i32, MAXE>>, i: i32, k: i32) -> vec3<f32> {
  var w = vec3<f32>(0.0, 0.0, 0.0);
  for (var m = 0; m < MAXE; m++) {
    let e = (*edges)[m];
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
fn ozoneAbsorb(path: f32) -> f32 {
  var a = 0.0;
  for (var b = 0; b < 8; b++) { a += O3_SHARE[b] * relaxedFraction(O3_COEF[b] * path); }
  return a;
}
fn visibleVapour(path: f32) -> f32 { return VIS_H2O_S * relaxedFraction(VIS_H2O_K * VAP_STRENGTH * path); }
fn nearInfraredGases(vapour: f32, oxygen: f32, co2: f32) -> f32 {
  var a = 0.0;
  for (var j = 0; j < 10; j++) { a += H2O_W[j] * relaxedFraction(H2O_K[j] * VAP_STRENGTH * vapour); }
  return VAPOR_ABS * a + O2_SHARE * relaxedFraction(O2_K * sqrt(oxygen)) + CO2_SW_K * sqrt(co2);
}
fn ozoneColumnAt(i: i32) -> f32 { let s = sin(MF[F_LAT + i]); return OZ_EQ + (OZ_POLE - OZ_EQ) * s * s; }
fn ozoneWeightsAt(i: i32) -> array<f32, 5> {
  let lat = MF[F_LAT + i]; let a = abs(lat) * 57.29577951308232;
  let summer = 0.5 * (1.0 + cos(6.283185307179586 * (P[1] - OZ_SUMMER)) * select(1.0, -1.0, lat < 0.0));
  let polar = clamp((a - 45.0) / 15.0, 0.0, 1.0); let middle = clamp((a - 15.0) / 30.0, 0.0, 1.0) * (1.0 - polar);
  return array<f32, 5>(1.0 - middle - polar, middle * summer, middle * (1.0 - summer), polar * summer, polar * (1.0 - summer));
}
fn ozoneAboveAt(p: f32, w: array<f32, 5>) -> f32 {
  if (!(p > 0.0)) { return 0.0; }
  let j = min(OZ_N - 2, i32(floor(log(p / OZ_LOW) / OZ_STEP)));
  var column = 0.0;
  if (j < 0) { for (var r = 0; r < 5; r++) { column += w[r] * OZ_TABLE[r * OZ_N]; } return column * p / OZ_LOW; }
  let t = (p - OZ_P[j]) / (OZ_P[j + 1] - OZ_P[j]);
  for (var r = 0; r < 5; r++) { let b = r * OZ_N + j; column += w[r] * (OZ_TABLE[b] + (OZ_TABLE[b + 1] - OZ_TABLE[b]) * t); }
  return column;
}
fn columnOzone(i: i32, pi: f32, layers: ptr<function, array<f32, K>>) -> f32 {
  if (!OZ_AFGL) { let column = ozoneColumnAt(i); for (var k = 0; k < K; k++) { (*layers)[k] = column * LV[L_OZS + k]; } return column; }
  let w = ozoneWeightsAt(i); var above = 0.0;
  for (var k = 0; k < K; k++) { let below = ozoneAboveAt(pi * LV[L_SL + k], w); (*layers)[k] = below - above; above = below; }
  return above;
}
fn solarGases(i: i32, pi: f32, mu: f32, layerOzone: ptr<function, array<f32, K>>, ozoneTaken: ptr<function, array<f32, K>>, gasTaken: ptr<function, array<f32, K>>) -> vec4<f32> {
  let magnification = 35.0 / sqrt(1224.0 * mu * mu + 1.0);
  var ozone = 0.0; var vapour = 0.0; var oxygen = 0.0; var co2 = 0.0;
  for (var k = 0; k < K; k++) {
    let idx = k * C + i; let p = pi * LV[L_SM + k]; let mass = pi * LV[L_DS + k] / GRAV; let q = max(0.0, IN[S_Q + idx]);
    let scaling = pow(p / SCALE_P, SCALE_N) * magnification; let dry = mass * max(0.0, 1.0 - q);
    ozone += (*layerOzone)[k] * magnification;
    vapour += q * mass * 0.1 * scaling * (1.0 + 0.00135 * (IN[S_TH + idx] * D[D_EXM + idx] - 240.0));
    oxygen += O2_PATH * dry * scaling; co2 += CO2_PATH * dry * scaling;
    (*ozoneTaken)[k] = ozoneAbsorb(ozone);
    (*gasTaken)[k] = VAPOR_ABS * visibleVapour(vapour) + nearInfraredGases(vapour, oxygen, co2);
  }
  return vec4<f32>(VAPOR_ABS * visibleVapour(vapour), vapour, oxygen, co2);
}
fn solarGasTotals(i: i32, pi: f32, mu: f32, layerOzone: ptr<function, array<f32, K>>) -> vec3<f32> {
  let magnification = 35.0 / sqrt(1224.0 * mu * mu + 1.0);
  var ozone = 0.0; var vapour = 0.0; var oxygen = 0.0; var co2 = 0.0;
  for (var k = 0; k < K; k++) {
    let idx = k * C + i; let p = pi * LV[L_SM + k]; let mass = pi * LV[L_DS + k] / GRAV; let q = max(0.0, IN[S_Q + idx]);
    let scaling = pow(p / SCALE_P, SCALE_N) * magnification; let dry = mass * max(0.0, 1.0 - q);
    ozone += (*layerOzone)[k] * magnification;
    vapour += q * mass * 0.1 * scaling * (1.0 + 0.00135 * (IN[S_TH + idx] * D[D_EXM + idx] - 240.0));
    oxygen += O2_PATH * dry * scaling; co2 += CO2_PATH * dry * scaling;
  }
  return vec3<f32>(VAPOR_ABS * visibleVapour(vapour), ozoneAbsorb(ozone), VAPOR_ABS * visibleVapour(vapour) + nearInfraredGases(vapour, oxygen, co2));
}
fn nearInfraredUpward(i: i32, pi: f32, down: vec3<f32>, upward: ptr<function, array<f32, K>>) -> f32 {
  var vapour = down.x; var oxygen = down.y; var co2 = down.z;
  var before = nearInfraredGases(vapour, oxygen, co2); var loss = 0.0;
  for (var k = K - 1; k >= 0; k--) {
    let idx = k * C + i; let p = pi * LV[L_SM + k]; let mass = pi * LV[L_DS + k] / GRAV; let q = max(0.0, IN[S_Q + idx]);
    let scaling = pow(p / SCALE_P, SCALE_N) * DIFFUSE_PATH; let dry = mass * max(0.0, 1.0 - q);
    vapour += q * mass * 0.1 * scaling * (1.0 + 0.00135 * (IN[S_TH + idx] * D[D_EXM + idx] - 240.0));
    oxygen += O2_PATH * dry * scaling; co2 += CO2_PATH * dry * scaling;
    let through = nearInfraredGases(vapour, oxygen, co2);
    (*upward)[k] = through - before;
    loss += through - before;
    before = through;
  }
  return loss;
}
fn longwavePaths(i: i32, k: i32, pi: f32, T: f32, ozone: f32) -> array<f32, 6> {
  let idx = k * C + i; let p = pi * LV[L_SM + k]; let mass = pi * LV[L_DS + k] / GRAV; let q = max(0.0, IN[S_Q + idx]);
  let scale = p / LW_PREF; let vapour = q * mass; let dry = max(0.0, 1.0 - q) * mass * scale;
  let co2 = max(0.0, 1.0 - q) * mass * sqrt(p * p + LW_DCO2 * LW_DCO2) / LW_PREF;
  return array<f32, 6>(vapour * sqrt(p * p + LW_DH2O * LW_DH2O) / LW_PREF, vapour * (q * p / (0.622 + 0.378 * q)) * exp(LW_TSELF * (1.0 / T - 1.0 / 296.0)), CO2_MASS * co2 * exp(LW_TCO2 * (T - 250.0)),
    ozone * OZ_KG * pow(sqrt(p * p + LW_DO3 * LW_DO3) / LW_PREF, LW_NO3), CH4_MASS * dry, N2O_MASS * dry);
}
fn longwaveDepth(g: i32, row: array<f32, 6>) -> f32 {
  let r = 6 * g;
  return LW_D * (LW_K[r] * row[0] + LW_K[r + 1] * row[1] + LW_K[r + 2] * row[2] + LW_K[r + 3] * row[3] + LW_K[r + 4] * row[4] + LW_K[r + 5] * row[5]);
}
fn planckShare(g: i32, t: f32) -> f32 {
  let b = 5 * g;
  return LW_PLANCK[b] + t * (LW_PLANCK[b + 1] + t * (LW_PLANCK[b + 2] + t * (LW_PLANCK[b + 3] + t * LW_PLANCK[b + 4])));
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
fn cloudOptics(T: f32, continental: bool) -> vec2<f32> {
  let liquid = clamp((T - ICE_T) / (LIQUID_T - ICE_T), 0.0, 1.0);
  let droplet = select(DROP_SEA, DROP_LAND, continental);
  let c = clamp(T - 273.15, ICE_COLDEST, ICE_WARMEST);
  let crystal = 0.5 * (326.3 + c * (12.42 + c * (0.197 + c * 0.0012)));
  let liquidVisible = 1500.0 / droplet; let iceVisible = 1000.0 * (3.448e-3 + 2.431 / crystal);
  let solar = liquid * (1.0 - (0.829 + 2.482e-3 * droplet)) * liquidVisible + (1.0 - liquid) * (1.0 - (0.7661 + 5.851e-4 * crystal)) * iceVisible;
  let infrared = IR_DIFFUSIVITY * 1000.0 * (liquid * LIQUID_IR + (1.0 - liquid) * (0.005 + 1.0 / crystal));
  return vec2<f32>(select(solar, CLOUD_SCAT, GRAY_SW), select(infrared, CLOUD_ABS, GRAY_LW));
}
fn cloudKeep(path: f32) -> f32 {
  if (CLOUD_SW > 0.0) { return exp(-CLOUD_SW * path); }
  return 1.0;
}
fn erfApprox(x: f32) -> f32 {
  let t = 1.0 / (1.0 + 0.3275911 * abs(x));
  let y = 1.0 - t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429)))) * exp(-x * x);
  return select(y, -y, x < 0.0);
}
fn varianceCover(deficit: f32, spread: f32) -> f32 {
  if (!(spread > 0.0)) { return select(select(0.5, 0.0, deficit < 0.0), 1.0, deficit > 0.0); }
  return 0.5 * (1.0 + erfApprox(deficit / (1.4142135623730951 * spread)));
}
fn turbulentCover(i: i32, k: i32, pi: f32, mixTop: f32) -> f32 {
  let bottom = (K - 1) * C + i; let idx = k * C + i; let ex = D[D_EXM + idx];
  let water = max(0.0, IN[S_QC + idx]); let total = max(0.0, IN[S_Q + idx]) + water;
  let level = IN[S_TH + idx] - LHEAT * water / (CP * ex); let liquidT = level * ex;
  let saturated = cloudSat(liquidT, pi * LV[L_SM + k]); let qs = saturated.x; let slope = saturated.y;
  let a = 1.0 / (1.0 + LHEAT / CP * slope);
  let zk = D[D_GEO + idx] + LV[L_GABS + k]; let zs = D[D_GEO + bottom] + LV[L_GABS + K - 1];
  var gradientQ = 0.0; var gradientL = 0.0; var n = 0.0;
  if (k > 0) {
    let jdx = idx - C; let zj = D[D_GEO + jdx] + LV[L_GABS + k - 1];
    if (!((zj - zs) / GRAV >= mixTop)) {
      let jWater = max(0.0, IN[S_QC + jdx]); let dz = (zj - zk) / GRAV;
      gradientQ += (max(0.0, IN[S_Q + jdx]) + jWater - total) / dz;
      gradientL += (IN[S_TH + jdx] - LHEAT * jWater / (CP * D[D_EXM + jdx]) - level) / dz;
      n += 1.0;
    }
  }
  if (k < K - 1) {
    let jdx = idx + C; let zj = D[D_GEO + jdx] + LV[L_GABS + k + 1];
    let jWater = max(0.0, IN[S_QC + jdx]); let dz = (zj - zk) / GRAV;
    gradientQ += (max(0.0, IN[S_Q + jdx]) + jWater - total) / dz;
    gradientL += (IN[S_TH + jdx] - LHEAT * jWater / (CP * D[D_EXM + jdx]) - level) / dz;
    n += 1.0;
  }
  if (n > 1.0) { gradientQ /= n; gradientL /= n; }
  let z = zk / GRAV;
  let asymptote = select(STABLE_LENGTH, MIX_LENGTH, PH[PH_BUOY + i] > 0.0);
  let length = 0.4 * z / (1.0 + 0.4 * z / asymptote);
  let spread = max(VAR_FLOOR * qs, VAR_SCALE * length * a * abs(gradientQ - ex * slope * gradientL));
  return varianceCover(a * (total - qs), spread);
}
fn layerCover(idx: i32, k: i32, bottom: i32, pi: f32, water: f32, mixedDepth: f32, stratiform: f32, mixTop: f32) -> f32 {
  if (!PDF_COVER || !(water > 0.0)) { return 1.0; }
  let height = (D[D_GEO + idx] + LV[L_GABS + k] - D[D_GEO + bottom] - LV[L_GABS + K - 1]) / GRAV;
  let inside = height < mixedDepth;
  if (VARIANCE_COVER && height < mixTop) {
    let fv = clamp(turbulentCover(idx - k * C, k, pi, mixTop), COVER_FLOOR, 1.0);
    if (!(stratiform > 0.0)) { return fv; }
    return (1.0 - stratiform) * fv + stratiform * overcastCover(idx, k, pi, RHC_BL);
  }
  let qsl = qsat(IN[S_TH + idx] * D[D_EXM + idx], pi * LV[L_SM + k]); let condensate = max(0.0, IN[S_QC + idx]);
  let excess = max(0.0, IN[S_Q + idx]) + condensate - qsl;
  var rhc = select(RHC, RHC_BL, inside);
  let width = (1.0 - rhc) * qsl;
  var f = clamp((excess + width) / (2.0 * width), COVER_FLOOR, 1.0);
  if (UNIFORM) {
    let pressure = pi * LV[L_SM + k];
    f = clamp(uniformCover(condensate, uniformWidth(cloudSat(IN[S_TH + idx] * D[D_EXM + idx] - LHEAT * condensate / CP, pressure), pressure, pi)), COVER_FLOOR, 1.0);
    rhc = criticalProfile(pressure, pi);
  }
  if (!(stratiform > 0.0)) { return f; }
  return (1.0 - stratiform) * f + stratiform * overcastCover(idx, k, pi, rhc);
}
fn overcastCover(idx: i32, k: i32, pi: f32, rhc: f32) -> f32 {
  let pressure = pi * LV[L_SM + k]; let condensate = max(0.0, IN[S_QC + idx]); let total = max(0.0, IN[S_Q + idx]) + condensate;
  if (UNIFORM) {
    let s = cloudSat(IN[S_TH + idx] * D[D_EXM + idx] - LHEAT * condensate / CP, pressure);
    let a = 1.0 / (1.0 + LHEAT * s.y / CP);
    let bound = min(a * (1.0 - rhc) * s.x, max(condensate, OVERCAST_WATER));
    return clamp((a * (total - s.x) + bound) / (2.0 * bound), COVER_FLOOR, 1.0);
  }
  let qs = qsat(IN[S_TH + idx] * D[D_EXM + idx], pressure);
  let bound = min((1.0 - rhc) * qs, max(condensate, OVERCAST_WATER));
  return clamp((total - qs + bound) / (2.0 * bound), COVER_FLOOR, 1.0);
}
fn stratiformShare(i: i32, bottom: i32, pi: f32) -> f32 {
  if (!(IN[S_Q + bottom] > 0.0)) { return 0.0; }
  let upper = STABILITY_K * C + i; let lowerT = IN[S_TH + bottom] * D[D_EXM + bottom];
  let lcl = condensationLevel(lowerT, IN[S_Q + bottom], pi * LV[L_SM + K - 1]);
  let base = max(0.0, CP * (lowerT - lcl.x) / GRAV);
  let height = (D[D_GEO + upper] - D[D_GEO + bottom] + LV[L_GABS + STABILITY_K] - LV[L_GABS + K - 1]) / GRAV;
  let inversion = inversionStrength(IN[S_TH + upper] - IN[S_TH + bottom], lowerT, IN[S_TH + upper] * D[D_EXM + upper], height - base);
  return clamp((inversion - OVERCAST_EIS) / OVERCAST_RAMP, 0.0, 1.0);
}
fn overlap(seen: f32, last: bool, blocks: ptr<function, vec3<f32>>) {
  if (seen > 0.0) { (*blocks).x = max((*blocks).x, seen); }
  if ((*blocks).x > 0.0 && (!(seen > 0.0) || last)) { (*blocks).y *= 1.0 - (*blocks).x; (*blocks).z = max((*blocks).z, (*blocks).x); (*blocks).x = 0.0; }
}
fn overlapCover(blocks: vec3<f32>, layered: vec2<f32>) -> f32 {
  if (!PDF_COVER) { return 0.0; }
  if (EXP_OVERLAP) { return layered.x; }
  return select(blocks.z, 1.0 - blocks.y, RANDOM_OVERLAP);
}
fn overlapAlpha(i: i32, k: i32) -> f32 {
  if (k == 0) { return 0.0; }
  let dz = (D[D_GEO + (k - 1) * C + i] + LV[L_GABS + k - 1] - D[D_GEO + k * C + i] - LV[L_GABS + k]) / GRAV;
  return exp(-dz / (DECOR_LEN - DECOR_SLOPE * abs(MF[F_LAT + i]) * 57.29577951308232));
}
fn overlapLayer(seen: f32, alpha: f32, layered: ptr<function, vec2<f32>>) {
  let previous = (*layered).y;
  if (previous < 1.0) {
    let pair = alpha * max(previous, seen) + (1.0 - alpha) * (previous + seen - previous * seen);
    (*layered).x = 1.0 - (1.0 - (*layered).x) * (1.0 - pair) / (1.0 - previous);
  } else { (*layered).x = 1.0; }
  (*layered).y = seen;
}
fn shortwave(cloudDepth: f32, keep: f32, mu: f32, adir: f32, adif: f32, light: vec3<f32>) -> vec4<f32> {
  if (!SCATTER) { return stream(cloudDepth, keep, mu, adir, adif); }
  return light.x * visibleStreams(cloudDepth, keep, mu, adir, adif, light) + (1.0 - light.x) * stream(cloudDepth + NIR_RAY * light.y, keep, mu, adir, adif);
}
fn streamEscape(cloudDepth: f32, keep: f32, mu: f32, adir: f32, adif: f32) -> f32 {
  let s = stream(cloudDepth, keep, mu, adir, adif);
  let reflectance = select(0.0, cloudDepth / (cloudDepth + 2.0 * mu), mu > 0.0 && cloudDepth > 0.0);
  return 1.0 - s.x - keep * reflectance - s.w;
}
fn reflection(cloudDepth: f32, keep: f32, mu: f32) -> f32 { return select(0.0, keep * cloudDepth / (cloudDepth + 2.0 * mu), mu > 0.0 && cloudDepth > 0.0); }
fn escapes(cloudDepth: f32, keep: f32, mu: f32, adir: f32, adif: f32, light: vec3<f32>) -> vec3<f32> {
  if (!SCATTER) { let rest = streamEscape(cloudDepth, keep, mu, adir, adif); return vec3<f32>(light.x * rest, (1.0 - light.x) * rest, light.x * reflection(cloudDepth, keep, mu)); }
  return vec3<f32>(light.x * visibleEscape(cloudDepth, keep, mu, adir, adif, light), (1.0 - light.x) * streamEscape(cloudDepth + NIR_RAY * light.y, keep, mu, adir, adif), light.x * visibleReflection(cloudDepth, keep, mu, light));
}
fn clearLight(beam: f32, mu: f32, ozoneHeating: f32, incident: f32, pi: f32, aerosol: f32) -> vec4<f32> {
  let visible = max(0.0, VIS_FRAC * beam - ozoneHeating);
  var taken = 0.0;
  if (mu > 0.0 && aerosol > 0.0) { taken = visible * (1.0 - exp(-AER_ABS * aerosol * 35.0 / sqrt(1224.0 * mu * mu + 1.0))); }
  let left = incident - taken;
  return vec4<f32>(taken, select(0.0, min(1.0, (visible - taken) / left), left > 0.0), pi, AER_SCAT * aerosol);
}
fn stream(cloudDepth: f32, keep: f32, mu: f32, adir: f32, adif: f32) -> vec4<f32> {
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
fn saturatedTemperature(energy: f32, pressure: f32, guess: f32) -> f32 {
  var t = guess;
  for (var n = 0; n < 4; n++) {
    let qs = qsat(t, pressure);
    t -= (CP * t + LHEAT * qs - energy) / (CP + LHEAT * LHEAT * qs / (RVAP * t * t));
  }
  return t;
}
fn relaxedSeries(x: f32) -> f32 { return x * (1.0 - 0.5 * x * (1.0 - x / 3.0 * (1.0 - 0.25 * x))); }
fn relaxedFraction(x: f32) -> f32 {
  if (x < 1e-2) { return relaxedSeries(x); }
  return 1.0 - exp(-x);
}
// both forms, no branch: the longwave's neighbouring columns straddle the cutoff often enough that branching costs more
fn relaxedFractionBoth(x: f32) -> f32 { return select(1.0 - exp(-x), relaxedSeries(x), x < 1e-2); }
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
fn thomasRate(n: i32, upper: ptr<function, array<f32, K>>, lower: ptr<function, array<f32, K>>, own: ptr<function, array<f32, K>>, rhs: ptr<function, array<f32, K>>, dt: f32) {
  var gain: array<f32, K>;
  var denominator = 1.0 + (*upper)[0] + (*lower)[0] + dt * (*own)[0];
  gain[0] = -(*lower)[0] / denominator;
  (*rhs)[0] = (*rhs)[0] / denominator;
  for (var j = 1; j < n; j++) {
    denominator = 1.0 + (*upper)[j] + (*lower)[j] + dt * (*own)[j] + (*upper)[j] * gain[j - 1];
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
struct MlmSun { incident: f32, mu: f32, adir: f32, adif: f32, path: f32, layer: f32, clear: f32, light: vec3<f32>, depth: f32, unit: f32 }
fn mlmFinite(x: f32) -> bool { return (bitcast<u32>(x) & 0x7f800000u) != 0x7f800000u; }
fn mlmFresh(x: f32) -> f32 {
  if (x < 1e-2) { return x * (1.0 - 0.5 * x * (1.0 - x / 3.0 * (1.0 - 0.25 * x))); }
  return 1.0 - exp(-x);
}
fn mlmRest(i: i32, resting: f32, dt: f32) {
  let h = PH[PH_MLMH + i];
  if (!(h > 0.0)) { return; }
  var settled = resting;
  if (MLM_HMEM > 0.0) { settled = h + (resting - h) * mlmFresh(dt / MLM_HMEM); }
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
  let sw = shortwave(select(sun.depth + sun.unit * water, CLOUD_SCAT * total, GRAY_SW), cloudKeep(total), sun.mu, sun.adir, sun.adif, sun.light);
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
fn mlmCapping(i: i32, floor: f32) -> i32 {
  for (var k = K - 2; k >= 1; k--) {
    if ((D[D_GEO + (k + 1) * C + i] + LV[L_GABS + k + 1]) / GRAV >= MLM_HMAX) { break; }
    if ((D[D_GEO + k * C + i] + LV[L_GABS + k]) / GRAV > floor && D[D_THV + k * C + i] - D[D_THV + (k + 1) * C + i] >= MLM_CEILINV) { return k; }
  }
  return -1;
}
fn mlmColumn(i: i32, pi: f32, mixedDepth: f32, sensible: f32, evaporation: f32, dt: f32, sun: MlmSun) -> MlmDeck {
  let none = MlmDeck(false, 0.0, 0.0, 0.0, 0.0);
  let bottom = (K - 1) * C + i;
  let depth = mixedDepth + (D[D_GEO + bottom] + LV[L_GABS + K - 1]) / GRAV;
  let regime = PH[PH_REGIME + i];
  let lifted = MLM_REST_REGIME && (regime == 1.0 || regime == 2.0);
  var capping = -1;
  if (MLM_PROGNOSTIC || lifted) { capping = mlmCapping(i, depth); }
  var ceiling = MLM_HMAX;
  if (MLM_PROGNOSTIC && capping >= 0) { ceiling = min(MLM_HMAX, (D[D_GEO + capping * C + i] + LV[L_GABS + capping]) / GRAV - 1.0); }
  var standDown = lifted;
  if (lifted && capping >= 0) { standDown = !(mlmInterface(i, capping + 1) <= MLM_CUCEIL); }
  let resting = select(depth, ceiling, MLM_REST_INVERSION && !(MLM_REST_REGIME && regime == 3.0) && !standDown && ceiling < MLM_HMAX);
  var h = depth;
  if (MLM_PROGNOSTIC) { h = resting; }
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
  if (MLM_FRACTION && k >= 1 && k < K - 1) {
    let upper = mlmInterface(i, k + 1);
    var part = k;
    if (h < upper) { part = k + 1; }
    if (part != capping) {
      let pdx = part * C + i;
      var lowerEdge = (D[D_GEO + bottom] + LV[L_GABS + K - 1] - CP * D[D_THV + bottom] * (D[D_EXL + bottom] - D[D_EXM + bottom])) / GRAV;
      if (part < K - 1) { lowerEdge = mlmInterface(i, part + 1); }
      let share = clamp((h - lowerEdge) / (mlmInterface(i, part) - lowerEdge), 0.0, 1.0) - select(0.0, 1.0, part == k + 1);
      let cloud = max(0.0, IN[S_QC + pdx]);
      heat += share * LV[L_DS + part] * (IN[S_TH + pdx] - LHEAT * cloud / (CP * D[D_EXM + pdx]));
      water += share * LV[L_DS + part] * (max(0.0, IN[S_Q + pdx]) + cloud);
      weight += share * LV[L_DS + part];
    }
  }
  if (k < 1) { mlmRest(i, resting, dt); return none; }
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
  var thetaAbove = IN[S_TH + above] - LHEAT * aboveCloud / (CP * D[D_EXM + above]);
  var qtAbove = max(0.0, IN[S_Q + above]) + aboveCloud;
  if (MLM_INTERPOLATE && k < K - 1) {
    let higher = above - C; let higherCloud = max(0.0, IN[S_QC + higher]);
    let below = D[D_GEO + above + C] + LV[L_GABS + k + 1];
    let share = clamp((GRAV * h - below) / (D[D_GEO + above] + LV[L_GABS + k] - below), 0.0, 1.0);
    thetaAbove += share * (IN[S_TH + higher] - LHEAT * higherCloud / (CP * D[D_EXM + higher]) - thetaAbove);
    qtAbove += share * (max(0.0, IN[S_Q + higher]) + higherCloud - qtAbove);
  }
  let start = MlmState(h, heat / weight, water / weight);
  var now = MlmOut(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  var passed = 0.0;
  if (sinking && !standDown) {
    now = mlmDiagnose(start, pi, sensible, evaporation, thetaAbove, qtAbove, sun);
    if (MLM_BL_GATE) { if (PH[PH_REGIME + i] == 3.0) { passed = 1.0; } } else if (now.jump >= MLM_MININV) { passed = 1.0; }
  }
  var gate = passed;
  if (MLM_GATEMEM > 0.0) { gate = PH[PH_MLMGATE + i] + (passed - PH[PH_MLMGATE + i]) * mlmFresh(dt / MLM_GATEMEM); }
  if (standDown) { gate = 0.0; }
  PH[PH_MLMGATE + i] = gate;
  if (!(gate > MLM_UNDECIDED || (gate == MLM_UNDECIDED && passed > 0.0)) || MLM_BYPASS) { mlmRest(i, resting, dt); return none; }
  if (!sinking) { now = mlmDiagnose(start, pi, sensible, evaporation, thetaAbove, qtAbove, sun); }
  var next = now; var top = h;
  if (dt > 0.0) {
    let deepened = h + dt * (now.entrainment + subsidence);
    top = deepened;
    if (MLM_PROGNOSTIC) { top = max(depth, min(ceiling, deepened)); }
    next = mlmDiagnose(MlmState(top, start.thetaL + dt * now.heat / deepened, start.qt + dt * now.water / deepened), pi, sensible, evaporation, thetaAbove, qtAbove, sun);
  }
  if (!(mlmFinite(next.lwp) && mlmFinite(next.cover) && mlmFinite(now.entrainment))) { mlmRest(i, resting, dt); return none; }
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
    PH[PH_SNOWALB + i] = ALB_FRESH;
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
    PH[PH_SNOWALB + i] = select(ALB_FRESH, agedSnow(PH[PH_SNOWALB + i], T, dt, ICE_SNOWFLOOR), snow > 0.0);
  }`;

/*
 * Snow of `amount` kg/m² falling on sea cell i, on its ice or into its water.
 */
export const snowOnSea = (amount) => `if (IN[S_ICE + i] > 0.0) {
      let conc = PH[PH_CONC + i]; let cover = select(conc, 1.0, conc <= 0.0);
      PH[PH_SNOW + i] += ${amount};
      PH[PH_SNOWALB + i] = refreshedSnow(PH[PH_SNOWALB + i], ${amount});
      IN[S_ICE + i] += (1.0 - cover) * (${amount}) / (RHOICE * cover);
    } else { IN[S_TS + i] -= LFUS * (${amount}) / PH[PH_CAP + i]; }`;

/*
 * The column physics runs as four kernels: physics (the surface's state,
 * the exchange and the deck), radiation (the clouds and the shortwave, on
 * full calls only when held), longwave (one thread per column, sweeping
 * the g-points a chunk at a time) and physicsSurface (the surface's energy
 * budget, the heating and the surface's own state). Each leaves what the
 * later ones read in the K1 register, free between the stages of a step:
 * per layer the ozone, the cloud emissivity, the overlap share and weight
 * and (unheld) the net flux, which the longwave turns into the net flux
 * after it; and per column the scalars stash() numbers.
 */
const LONGWAVE_STASH = `const LWS_CLOUD: i32 = 0; const LWS_SHARE: i32 = K * C; const LWS_OZONE: i32 = 2 * K * C; const LWS_FLUX: i32 = 3 * K * C; const LWS_ALPHA: i32 = 4 * K * C; const LWS_SCALAR: i32 = 5 * K * C;
fn stashed(j: i32, i: i32) -> f32 { return OUT[LWS_SCALAR + j * C + i]; }
fn stash(j: i32, i: i32, x: f32) { OUT[LWS_SCALAR + j * C + i] = x; }
`;
export const LONGWAVE_STASH_FLOATS = (K, C) => 5 * K * C + 26 * C;
/*
 * The correlated longwave's sweeps with the g-points taken a chunk at a
 * time, each chunk's fluxes held in scalars: the chunk's code is unrolled
 * here, so the g-points' fluxes stay in registers where arrays indexed by
 * g would be private memory the kernel streams through on every layer.
 * Each layer's paths and overlap weights are found again for every chunk,
 * and the layer's heat carries from chunk to chunk in g-point order, so
 * every sum adds the same terms in the same order as one pass over all
 * g-points would.
 */
const LONGWAVE_DOWN = 12, LONGWAVE_UP = 6;
function longwaveSweeps(held, ng, { layer, chainF, cloudE, alpha, downEnd, upEnd }) {
  const split = (size) => Array.from({ length: Math.ceil(ng / size) }, (_, n) => Array.from({ length: Math.min(size, ng - n * size) }, (_, j) => n * size + j));
  const downChunks = split(LONGWAVE_DOWN), upChunks = split(LONGWAVE_UP);
  const each = (G, line) => G.map(line).join('');
  const chainDown = `
        let f = ${chainF('k')}; let a = select(0.0, ${chainF('max(k - 1, 0)')}, k > 0); let alpha = ${alpha('k')};
        let cloudIn = select(0.0, ${cloudE} / f, f > 0.0); let clearIn = select(${cloudE}, 0.0, f > 0.0);
        let both = a + f - (alpha * max(a, f) + (1.0 - alpha) * (a + f - a * f));
        let p00 = 1.0 - a - f + both; let p01 = f - both; let p10 = a - both;
        let clearShare = select(0.0, 1.0 / (1.0 - a), a < 1.0); let cloudShare = select(0.0, 1.0 / a, a > 0.0);`;
  const chainUp = `
        let f = ${chainF('k')}; let b = select(0.0, ${chainF('min(k + 1, K - 1)')}, k < K - 1); let alpha = ${alpha('min(k + 1, K - 1)')};
        let cloudIn = select(0.0, ${cloudE} / f, f > 0.0); let clearIn = select(${cloudE}, 0.0, f > 0.0);
        let both = f + b - (alpha * max(f, b) + (1.0 - alpha) * (f + b - f * b));
        let p00 = 1.0 - f - b + both; let p01 = b - both; let p10 = f - both;
        let clearShare = select(0.0, 1.0 / (1.0 - b), b < 1.0); let cloudShare = select(0.0, 1.0 / b, b > 0.0);`;
  const carried = downChunks.length > 1 || upChunks.length > 1 ? `    var lwHeat: array<f32, K>;${held && upChunks.length > 1 ? ' var lwTaken: array<f32, K>;' : ''}\n` : '';
  const open = (n, name) => (n === 0 ? '0.0' : `lw${name}[k]`);
  const close = (n, last) => (n === upChunks.length - 1 ? last : `lwHeat[k] = heat;${held ? ' lwTaken[k] = taken;' : ''}`);
  const sweeps = (chain) => {
    let text = carried;
    downChunks.forEach((G, n) => {
      text += `    {
${each(G, (g) => (chain ? `      var down0_${g} = 0.0; var down1_${g} = 0.0;\n` : `      var down_${g} = 0.0;\n`))}      for (var k = 0; k < K; k++) {
        ${layer}${chain ? chainDown : ` let clear = 1.0 - ${cloudE};`}
        var heat = ${open(n, 'Heat')};
${each(G, (g) => (chain ? `        {
          let gas = relaxedFractionBoth(longwaveDepth(${g}, row)); let src = planckShare(${g}, t) * hot;
          let e0 = 1.0 - (1.0 - gas) * (1.0 - clearIn); let e1 = 1.0 - (1.0 - gas) * (1.0 - cloudIn);
          let from0 = down0_${g} * clearShare; let from1 = down1_${g} * cloudShare;
          let in0 = from0 * p00 + from1 * p10; let in1 = from0 * p01 + from1 * both;
          heat += e0 * (in0 - src * (1.0 - f)) + e1 * (in1 - src * f);
          down0_${g} = in0 * (1.0 - e0) + e0 * src * (1.0 - f); down1_${g} = in1 * (1.0 - e1) + e1 * src * f;
        }
` : `        {
          let e = 1.0 - (1.0 - relaxedFractionBoth(longwaveDepth(${g}, row))) * clear; let src = planckShare(${g}, t) * hot;
          heat += e * (down_${g} - src);
          down_${g} = down_${g} * (1.0 - e) + e * src;
        }
`))}        ${n === downChunks.length - 1 ? downEnd : 'lwHeat[k] = heat;'}
      }
${each(G, (g) => (chain ? `      back += down0_${g} + down1_${g};\n` : `      back += down_${g};\n`))}    }
`;
    });
    text += '    let surfaceScaled = (ts - 250.0) / 100.0;\n';
    upChunks.forEach((G, n) => {
      text += `    {
${each(G, (g) => (chain ? `      var up0_${g} = planckShare(${g}, surfaceScaled) * surfaceEmission; var up1_${g} = 0.0; var clear_${g} = up0_${g};\n` : `      var up_${g} = planckShare(${g}, surfaceScaled) * surfaceEmission; var clear_${g} = up_${g};\n`))}${held ? each(G, (g) => (chain ? `      var lift0_${g} = planckShare(${g}, surfaceScaled); var lift1_${g} = 0.0; var clearLift_${g} = lift0_${g};\n` : `      var lift_${g} = planckShare(${g}, surfaceScaled); var clearLift_${g} = lift_${g};\n`)) : ''}      for (var k = K - 1; k >= 0; k--) {
        ${layer}${chain ? chainUp : ` let clear = 1.0 - ${cloudE};`}
        var heat = ${open(n, 'Heat')};${held ? ` var taken = ${open(n, 'Taken')};` : ''}
${each(G, (g) => (chain ? `        {
          let gas = relaxedFractionBoth(longwaveDepth(${g}, row)); let src = planckShare(${g}, t) * hot;
          let e0 = 1.0 - (1.0 - gas) * (1.0 - clearIn); let e1 = 1.0 - (1.0 - gas) * (1.0 - cloudIn);
          let from0 = up0_${g} * clearShare; let from1 = up1_${g} * cloudShare;
          let in0 = from0 * p00 + from1 * p01; let in1 = from0 * p10 + from1 * both;
          heat += e0 * (in0 - src * (1.0 - f)) + e1 * (in1 - src * f);
          up0_${g} = in0 * (1.0 - e0) + e0 * src * (1.0 - f); up1_${g} = in1 * (1.0 - e1) + e1 * src * f;
          if (CLEAR_SKY) { clear_${g} = clear_${g} * (1.0 - gas) + gas * src; }
${held ? `          let rise0 = lift0_${g} * clearShare; let rise1 = lift1_${g} * cloudShare;
          let reach0 = rise0 * p00 + rise1 * p01; let reach1 = rise0 * p10 + rise1 * both;
          taken += e0 * reach0 + e1 * reach1; lift0_${g} = reach0 * (1.0 - e0); lift1_${g} = reach1 * (1.0 - e1); clearLift_${g} *= 1.0 - gas;
` : ''}        }
` : `        {
          let gas = relaxedFractionBoth(longwaveDepth(${g}, row)); let e = 1.0 - (1.0 - gas) * clear; let src = planckShare(${g}, t) * hot;
          heat += e * (up_${g} - src);
          up_${g} = up_${g} * (1.0 - e) + e * src;
          if (CLEAR_SKY) { clear_${g} = clear_${g} * (1.0 - gas) + gas * src; }
${held ? `          taken += e * lift_${g}; lift_${g} *= 1.0 - e; clearLift_${g} *= 1.0 - gas;\n` : ''}        }
`))}        ${close(n, `${upEnd}${held ? ' PH[PH_RADDF + k * C + i] = taken;' : ''}`)}
      }
${each(G, (g) => (chain ? `      outgoing += up0_${g} + up1_${g}; clearOutgoing += clear_${g};\n` : `      outgoing += up_${g}; clearOutgoing += clear_${g};\n`))}${held ? each(G, (g) => (chain ? `      escape += lift0_${g} + lift1_${g}; clearEscape += clearLift_${g};\n` : `      escape += lift_${g}; clearEscape += clearLift_${g};\n`)) : ''}    }
`;
    });
    return text;
  };
  return `  if (LW_CORRELATED && LW_CHAIN) {
    // each region's heat is its absorption less its emission: the difference of whole fluxes loses the thin top layers to f32 cancellation
${sweeps(true)}  } else if (LW_CORRELATED) {
${sweeps(false)}  }`;
}
// The kernels' held variants hold the radiation between full calls as radiation.module.js's applyHeld does.
const HELD_GATE = `  if (!(PH[PH_RADP] > 0.5 || !(PH[PH_RADEMIT + i] > 0.0))) { stash(0, i, 0.0); return; }
`;
const HELD_OPEN = `  var shone = 0.0; var weighted = 0.0;
  for (var n = 0; n < i32(PH[PH_RADP + 1]); n++) {
    let turn = f32(n) * PH[PH_RADP + 2]; let c = cos(turn); let s = sin(turn);
    let ahead = MF[F_XC + 3 * i] * (sun.x * c + sun.y * s) + MF[F_XC + 3 * i + 1] * (sun.y * c - sun.x * s) + MF[F_XC + 3 * i + 2] * sun.z;
    if (ahead > 0.0) { shone += ahead; weighted += ahead * ahead; }
  }
  let mu = select(0.0, weighted / shone, shone > 0.0); let beam = S0 * mu;
  PH[PH_RADMU + i] = mu;
  let waterDir = openWaterAlbedo(mu); let iceDir = surfaceAlbedo(ice, waterDir, snow0, skin, PH[PH_SNOWALB + i]);
  let adir = select(cover * iceDir + (1.0 - cover) * waterDir, landAlbedo, onLand);
  var ozoneTaken: array<f32, K>; var gasTaken: array<f32, K>; var solar = vec4<f32>(0.0, 0.0, 0.0, 0.0);
  if (SOLAR_CLIRAD && mu > 0.0) { solar = solarGases(i, pi, mu, &layerOzone, &ozoneTaken, &gasTaken); }
  let ozoneHeating = select(beam * OZONE_ABS, beam * ozoneTaken[K - 1], SOLAR_CLIRAD);
  let visibleTaken = ozoneHeating + beam * solar.x;
`;
const HELD_STORE = `  for (var k = 0; k < K; k++) { PH[PH_RADSW + k * C + i] = beforeBands[k]; }
  let atmosphereSolar = ozoneHeating + vaporHeating + aerosolHeating + cloudHeating + upwardHeating;
  PH[PH_RADABS + i] = absorbed; PH[PH_RADEMIT + i] = surfaceEmission;
  PH[PH_RADSWDN + i] = incident * sw.y; PH[PH_RADDIR + i] = incident * sw.z;
  PH[PH_RADATM + i] = atmosphereSolar; PH[PH_RADTOA + i] = absorbed + atmosphereSolar; PH[PH_RADREFL + i] = incident - absorbed - cloudHeating - upwardHeating;
  if (CLEAR_SKY) {
    let clearSw = shortwave(0.0, 1.0, mu, adir, adif, light);
    var clearUp = 0.0;
    if (UPWARD && mu > 0.0) { let clearEsc = escapes(0.0, 1.0, mu, adir, adif, light); clearUp = clearEsc.y * restLoss + clearEsc.x * aerosolLoss + select(0.0, (clearEsc.x * (1.0 - aerosolLoss) + clearEsc.z) * ozoneLoss, ozoneLoss > 0.0); }
    PH[PH_RADCLRSW + i] = incident * (clearSw.x + clearUp) + ozoneHeating + vaporHeating + aerosolHeating;
  }
`;
const HELD_GRAY = `    let outgoing = v.x + g.x + w.x; let back = v.y + g.y + w.y;
    var through = vec3<f32>(VAPOR_FRAC, GAS_FRAC, WINDOW);
    for (var k = K - 1; k >= 0; k--) {
      PH[PH_RADDF + k * C + i] = through.x * vaporE[k] + through.y * mixedE[k] + through.z * cloudE[k];
      through *= vec3<f32>(1.0 - vaporE[k], 1.0 - mixedE[k], 1.0 - cloudE[k]);
    }
    let escape = through.x + through.y + through.z;
    for (var k = 0; k < K; k++) { PH[PH_LWH + k * C + i] = netFlux[k] - beforeBands[k]; }
    PH[PH_RADBACK + i] = back; PH[PH_RADOLR + i] = outgoing; PH[PH_RADT + i] = escape;
    if (CLEAR_SKY) {
      var upVapor = VAPOR_FRAC * surfaceEmission; var upGas = GAS_FRAC * surfaceEmission; var throughVapor = VAPOR_FRAC; var throughGas = GAS_FRAC;
      for (var k = K - 1; k >= 0; k--) {
        let t = temperature[k]; let ev = clearE[k]; let eg = LV[L_GASE + k];
        upVapor = upVapor * (1.0 - ev) + VAPOR_FRAC * ev * STEFAN * t * t * t * t;
        upGas = upGas * (1.0 - eg) + GAS_FRAC * eg * STEFAN * t * t * t * t;
        throughVapor *= 1.0 - ev; throughGas *= 1.0 - eg;
      }
      PH[PH_RADCLROLR + i] = upVapor + upGas + WINDOW * surfaceEmission;
      PH[PH_RADTC + i] = throughVapor + throughGas + WINDOW;
    }
`;
const HELD_REST = `  let held = PH[PH_RADMU + i]; let sunScale = select(0.0, mu / held, held > 0.0); let lift = surfaceEmission - PH[PH_RADEMIT + i];
  let absorbed = sunScale * PH[PH_RADABS + i]; let back = PH[PH_RADBACK + i]; let outgoing = PH[PH_RADOLR + i] + PH[PH_RADT + i] * lift;
`;
const HELD_THETA = `  for (var k = 0; k < K; k++) {
    let idx = k * C + i;
    let mass = pi * LV[L_DS + k] / GRAV;
    let heat = sunScale * PH[PH_RADSW + idx] + PH[PH_LWH + idx] + PH[PH_RADDF + idx] * lift + select(0.0, sensible, k == K - 1);
    IN[S_TH + idx] += dt * heat / (CP * mass) / D[D_EXM + idx];
  }
`;
const HELD_TOA = `  let atmosphereSolar = sunScale * PH[PH_RADATM + i]; let absorbedSolar = sunScale * PH[PH_RADTOA + i]; let reflectedSolar = sunScale * PH[PH_RADREFL + i];
`;
const HELD_CLEAR = `  if (CLEAR_SKY) { PH[PH_ABSCLRSUM + i] += sunScale * PH[PH_RADCLRSW + i]; PH[PH_OLRCLRSUM + i] += PH[PH_RADCLROLR + i] + PH[PH_RADTC + i] * lift; }
`;
const PHYSICS_KERNEL = `${MIXED_LAYER_WGSL}${EXCHANGE_WGSL}${LONGWAVE_STASH}fn cumulusCloud(k: i32, i: i32) -> vec2<f32> {
  if (!CU_CLOUD || k < CU_K0) { return vec2<f32>(0.0, 0.0); }
  let slot = (k - CU_K0) * C + i;
  return vec2<f32>(PH[PH_CUCOVER + slot], PH[PH_CUCOVER + slot] * PH[PH_CUWATER + slot]);
}
@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= C) { return; }
  let pi = IN[S_PI + i]; let skin = IN[S_TS + i]; let ice = IN[S_ICE + i]; let bottom = (K - 1) * C + i;
  let wind = cellWind(i, K - 1);
  let ws = length(wind);
  D[D_WIND + i] = ws;
  let sun = vec3<f32>(P[2], P[3], P[4]);
  let mu = max(0.0, MF[F_XC + 3 * i] * sun.x + MF[F_XC + 3 * i + 1] * sun.y + MF[F_XC + 3 * i + 2] * sun.z);
  let beam = S0 * mu;
  let onLand = PH[PH_LAND + i] > 0.5; let onIceSheet = PH[PH_LAND + i] > 1.5; let continental = onLand && !onIceSheet;
  let soil0 = PH[PH_SOIL + i]; let snow0 = PH[PH_SNOW + i]; let veg0 = PH[PH_VEG + i]; let surf0 = PH[PH_SURF + i];
  let bucket = select(BUCKET, ROOTCAP, VEGETATED);
  let dryAlbedo = select(ALB_BARE, ALB_MINERAL - (ALB_MINERAL - ALB_HUMUS) * (1.0 - exp(-HUMUS_SCALE * max(0.0, PH[PH_SOILC + i]))), HUMIC);
  let soilAlbedo = select(dryAlbedo, dryAlbedo - (dryAlbedo - ALB_WETSOIL * (dryAlbedo / ALB_BARE)) * clamp((select(soil0 / ROOTCAP, surf0 / SURFCAP, DARK_SURFACE) - DARK_FROM) / DARK_SPAN, 0.0, 1.0), DARKENING);
  let trees0 = PH[PH_CANOPY + i];
  let coverAlbedo = select(ALB_VEG, select(ALB_GRASS, ALB_GRASS + (ALB_FORESTV - ALB_GRASS) * min(1.0, trees0 / veg0), veg0 > 0.0), GRASSY);
  let bareAlbedo = select(ALB_LAND, soilAlbedo + (coverAlbedo - soilAlbedo) * veg0, VEGETATED);
  let ownSnow = select(ALB_SNOW, PH[PH_SNOWALB + i], AGEING) - select(0.0, GRASS_SNOW * max(0.0, veg0 - trees0), GRASSY);
  let coveredSnow = select(ownSnow, ownSnow + (ALB_FOREST - ownSnow) * min(1.0, PH[PH_CANOPY + i] / CLOSED_CANOPY), MASKED);
  let landAlbedo = select(bareAlbedo + min(1.0, snow0 / FULLSNOW) * (coveredSnow - bareAlbedo), ALB_ICESHEET, onIceSheet);
  let conc0 = PH[PH_CONC + i];
  let cover = select(0.0, select(conc0, 1.0, conc0 <= 0.0), ice > 0.0);
  let waterDir = openWaterAlbedo(mu);
  let iceDif = surfaceAlbedo(ice, ALB_DIF_WATER, snow0, skin, PH[PH_SNOWALB + i]); let iceDir = surfaceAlbedo(ice, waterDir, snow0, skin, PH[PH_SNOWALB + i]);
  let adif = select(cover * iceDif + (1.0 - cover) * ALB_DIF_WATER, landAlbedo, onLand);
  let adir = select(cover * iceDir + (1.0 - cover) * waterDir, landAlbedo, onLand);
  let ts = select(skin, cover * skin + (1.0 - cover) * FREEZING, !onLand && ice > 0.0 && cover < 1.0);
  let warmth = clamp((ts - GROWCOLD) / (GROWWARM - GROWCOLD), 0.0, 1.0);
  let roots = min(1.0, soil0 / (WETT * bucket));
  let bareWet = (1.0 - veg0) * min(1.0, surf0 / SURFCAP);
  if (ROUGH) {
    let before = select(0.0, xWetness(PH[PH_HEATX + i] * xWind(i, ws), roots, bareWet, veg0, snow0, warmth), onLand);
    let exchange = surfaceExchange(i, pi, ts, ws, select(0.0, cover, !onLand), snow0, veg0, trees0, onLand, onIceSheet, before, PH[PH_DEPTH + i] - (D[D_GEO + bottom] + LV[L_GABS + K - 1]) / GRAV);
    PH[PH_DRAG + i] = exchange.x; PH[PH_HEATX + i] = exchange.y; PH[PH_REFX + i] = exchange.z;
  }
  let aero = PH[PH_HEATX + i] * xWind(i, ws);
  let canopyWet = veg0 * roots / (1.0 + RSTOM * aero / max(0.05, warmth));
  let landWet = select(roots, bareWet + canopyWet, VEGETATED);
  let wetness = select(1.0, select(landWet, 1.0, snow0 > ${TRACE_SNOW}), onLand);
  let bareShare = select(0.0, bareWet / max(1e-12, bareWet + canopyWet), VEGETATED && snow0 <= ${TRACE_SNOW});
  var gases = vec3<f32>(0.0, 0.0, 0.0); var layerOzone: array<f32, K>;
  let ozoneColumn = columnOzone(i, pi, &layerOzone);
  if (SOLAR_CLIRAD && mu > 0.0) { gases = solarGasTotals(i, pi, mu, &layerOzone); }
  let ozoneHeating = select(beam * OZONE_ABS, beam * gases.y, SOLAR_CLIRAD);
  let visibleTaken = ozoneHeating + beam * gases.x;
  let aerosol = select(SEA_AER, LAND_AER, onLand && !onIceSheet);
  let surfaceEmission = STEFAN * ts * ts * ts * ts;
  var deck = 0.0; var fraction = 0.0; var mlmCover = 0.0; var mlmWater = 0.0; var mlmEntrainment = 0.0; var mlmTop = 0.0;
  let mixedDepth = PH[PH_DEPTH + i] - (D[D_GEO + bottom] + LV[L_GABS + K - 1]) / GRAV;
  let mixTop = PH[PH_MIXTOP + i] - (D[D_GEO + bottom] + LV[L_GABS + K - 1]) / GRAV;
  let inversionShare = stratiformShare(i, bottom, pi);
  PH[PH_STRAT + i] = select(0.0, inversionShare, !onLand && 1.0 - cover > 0.0);
  let stratiform = select(0.0, inversionShare, PDF_COVER && BOUND_WIDTH);
  if (STRATUS && MLM_DECK && !onLand && 1.0 - cover > 0.0 && mixedDepth > 0.0) {
    let airT = IN[S_TH + bottom] * D[D_EXM + bottom];
    let rho = pi * LV[L_SM + K - 1] / (RGAS * airT);
    let exchange = rho * PH[PH_HEATX + i] * xWind(i, ws);
    let sensible = select(exchange * CP * (ts - airT), exchange * (CP * (ts - airT) - CP * D[D_THV + bottom] * (D[D_EXL + bottom] - D[D_EXM + bottom])), ROUGH);
    let evap = wetness * max(0.0, exchange * (qsat(ts, pi) - IN[S_Q + bottom]));
    var deckSun = MlmSun(0.0, mu, adir, adif, 0.0, 0.0, 0.0, vec3<f32>(0.0, 0.0, 0.0), 0.0, 0.0);
    if (STRATUS_SOLAR && CLOUD_SW > 0.0) {
      var path = 0.0; var depth = 0.0; var unit = 0.0; var layer = 0.0; var shadeBlocks = vec3<f32>(0.0, 1.0, 0.0); var shadeLayered = vec2<f32>(0.0, 0.0);
      for (var k = 0; k < K; k++) {
        let mass = pi * LV[L_DS + k] / GRAV;
        var water = max(0.0, IN[S_QC + k * C + i]) * mass;
        var f = layerCover(k * C + i, k, bottom, pi, water, mixedDepth, stratiform, mixTop);
        let cu = cumulusCloud(k, i);
        if (cu.y * mass > 0.0) { f = select(cu.x, max(f, cu.x), water > 0.0); water += cu.y * mass; }
        path += water;
        let optics = cloudOptics(IN[S_TH + k * C + i] * D[D_EXM + k * C + i], continental);
        depth += optics.x * water;
        if (k == STRATUS_K) { layer = water; unit = optics.x; }
        let shadeSeen = plumeSeen(select(0.0, f * smallRate(water / VISIBLE_PATH), PDF_COVER && water > 0.0), cu, mass);
        overlap(shadeSeen, k == K - 1, &shadeBlocks);
        if (EXP_OVERLAP) { overlapLayer(shadeSeen, overlapAlpha(i, k), &shadeLayered); }
      }
      var shade = overlapCover(shadeBlocks, shadeLayered);
      if (!(shade > 0.0)) { shade = 1.0; }
      var lit = beam - ozoneHeating;
      if (SOLAR_CLIRAD) { lit -= beam * gases.z; }
      else if (VAPOR_ABS > 0.0 && mu > 0.0) {
        let magnification = 35.0 / sqrt(1224.0 * mu * mu + 1.0);
        var vapor = 0.0;
        for (var k = 0; k < K; k++) { vapor += max(0.0, IN[S_Q + k * C + i]) * pi * LV[L_DS + k] / GRAV * sqrt(LV[L_SM + k]) * 0.1 * magnification; }
        lit -= lit * (VAPOR_ABS * 2.9 * vapor / (pow(1.0 + 141.5 * vapor, 0.635) + 5.925 * vapor));
      }
      var deckLight = vec3<f32>(0.0, 0.0, 0.0);
      if (SCATTER) { let clear = clearLight(beam, mu, visibleTaken, lit, pi, aerosol); lit -= clear.x; deckLight = clear.yzw; }
      var sky = shortwave(select(depth / shade, CLOUD_SCAT * path / shade, GRAY_SW), cloudKeep(path / shade), mu, adir, adif, deckLight);
      if (PDF_COVER && shade < 1.0) { sky = shade * sky + (1.0 - shade) * shortwave(0.0, 1.0, mu, adir, adif, deckLight); }
      deckSun = MlmSun(lit, mu, adir, adif, path, layer, select(0.0, lit * sky.w / path, path > 0.0) * layer, deckLight, depth, unit);
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
  for (var k = 0; k < K; k++) { OUT[LWS_OZONE + k * C + i] = layerOzone[k]; }
  stash(1, i, ws); stash(2, i, mu); stash(3, i, beam); stash(4, i, cover); stash(5, i, waterDir); stash(6, i, iceDif); stash(7, i, iceDir); stash(8, i, adif);
  stash(9, i, ts); stash(10, i, warmth); stash(11, i, wetness); stash(12, i, bareShare); stash(13, i, surfaceEmission);
  stash(22, i, landAlbedo); stash(23, i, stratiform); stash(24, i, ozoneColumn); stash(25, i, adir);
}`;
export const radiationKernel = (held) => `${LONGWAVE_STASH}fn cumulusCloud(k: i32, i: i32) -> vec2<f32> {
  if (!CU_CLOUD || k < CU_K0) { return vec2<f32>(0.0, 0.0); }
  let slot = (k - CU_K0) * C + i;
  return vec2<f32>(PH[PH_CUCOVER + slot], PH[PH_CUCOVER + slot] * PH[PH_CUWATER + slot]);
}
@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= C) { return; }
${held ? HELD_GATE : ''}  let pi = IN[S_PI + i]; let skin = IN[S_TS + i]; let ice = IN[S_ICE + i]; let bottom = (K - 1) * C + i;
  let sun = vec3<f32>(P[2], P[3], P[4]);
  let onLand = PH[PH_LAND + i] > 0.5; let onIceSheet = PH[PH_LAND + i] > 1.5; let continental = onLand && !onIceSheet;
  let snow0 = PH[PH_SNOW + i];
${held ? '' : `  let mu = stashed(2, i); let beam = stashed(3, i); let adir = stashed(25, i);
`}  let cover = stashed(4, i); let adif = stashed(8, i); let ts = stashed(9, i); let surfaceEmission = stashed(13, i);
  let landAlbedo = stashed(22, i); let stratiform = stashed(23, i); let ozoneColumn = stashed(24, i);
  var layerOzone: array<f32, K>;
  for (var k = 0; k < K; k++) { layerOzone[k] = OUT[LWS_OZONE + k * C + i]; }
  let aerosol = select(SEA_AER, LAND_AER, onLand && !onIceSheet);
  let tau0 = PH[PH_TAU + i]; let deck = PH[PH_DECK + i]; let fraction = PH[PH_DECKF + i];
  let mixedDepth = PH[PH_DEPTH + i] - (D[D_GEO + bottom] + LV[L_GABS + K - 1]) / GRAV;
  let mixTop = PH[PH_MIXTOP + i] - (D[D_GEO + bottom] + LV[L_GABS + K - 1]) / GRAV;
  var vaporE: array<f32, K>; var mixedE: array<f32, K>; var cloudE: array<f32, K>; var clearE: array<f32, K>; var temperature: array<f32, K>; var netFlux: array<f32, K>; var chainF: array<f32, K>;
  var cloudPath = 0.0; var cloudDepth = 0.0; var deckUnit = 0.0; var blocks = vec3<f32>(0.0, 1.0, 0.0); var layered = vec2<f32>(0.0, 0.0);
${held ? '' : `  var ozoneTaken: array<f32, K>; var gasTaken: array<f32, K>; var solar = vec4<f32>(0.0, 0.0, 0.0, 0.0);
  if (SOLAR_CLIRAD && mu > 0.0) { solar = solarGases(i, pi, mu, &layerOzone, &ozoneTaken, &gasTaken); }
  let ozoneHeating = select(beam * OZONE_ABS, beam * ozoneTaken[K - 1], SOLAR_CLIRAD);
  let visibleTaken = ozoneHeating + beam * solar.x;
`}${held ? HELD_OPEN : ''}  for (var k = 0; k < K; k++) {
    let idx = k * C + i;
    let mass = pi * LV[L_DS + k] / GRAV;
    var eps = 1.0 - exp(-tau0 * LV[L_SHAPE + k]);
    if (COUPLED) { eps = 1.0 - exp(-VCOUP * max(0.0, IN[S_Q + idx]) * mass); }
    if (CLEAR_SKY) { clearE[k] = eps; }
    var water = max(0.0, IN[S_QC + idx]) * mass;
    var f = layerCover(idx, k, bottom, pi, water, mixedDepth, stratiform, mixTop);
    let cu = cumulusCloud(k, i);
    if (cu.y * mass > 0.0) { f = select(cu.x, max(f, cu.x), water > 0.0); water += cu.y * mass; }
    cloudPath += water;
    let optics = cloudOptics(IN[S_TH + idx] * D[D_EXM + idx], continental);
    cloudDepth += optics.x * water;
    if (k == STRATUS_K) { deckUnit = optics.x; }
    let seen = plumeSeen(select(0.0, f * smallRate(water / VISIBLE_PATH), PDF_COVER && water > 0.0), cu, mass);
    overlap(seen, k == K - 1, &blocks);
    let alpha = overlapAlpha(i, k);
    if (EXP_OVERLAP) { overlapLayer(seen, alpha, &layered); }
    if (LW_CORRELATED && LW_CHAIN) { OUT[LWS_ALPHA + idx] = alpha; }
    cloudE[k] = select(0.0, f * (1.0 - exp(-optics.y * water / f)), water > 0.0);
    if (STRATUS && k == STRATUS_K && deck > 0.0) { cloudE[k] = fraction * (1.0 - exp(-optics.y * (water + deck))) + (1.0 - fraction) * cloudE[k]; }
    if (LW_CHAIN) { chainF[k] = select(0.0, f, water > 0.0 && !(STRATUS && k == STRATUS_K && deck > 0.0)); }
    let clear = 1.0 - cloudE[k];
    vaporE[k] = 1.0 - (1.0 - eps) * clear;
    mixedE[k] = 1.0 - (1.0 - LV[L_GASE + k]) * clear;
    temperature[k] = IN[S_TH + idx] * D[D_EXM + idx];
    netFlux[k] = select(ozoneHeating * LV[L_OZ + k], beam * (ozoneTaken[k] - select(0.0, ozoneTaken[max(k - 1, 0)], k > 0)), SOLAR_CLIRAD);
  }
  var incident = beam - ozoneHeating;
  var vaporHeating = 0.0; var path = 0.0;
  if (SOLAR_CLIRAD) {
    var taken = 0.0;
    for (var k = 0; k < K; k++) { netFlux[k] += beam * (gasTaken[k] - taken); taken = gasTaken[k]; }
    vaporHeating = beam * taken;
    incident -= vaporHeating;
  } else if (VAPOR_ABS > 0.0 && mu > 0.0) {
    let magnification = 35.0 / sqrt(1224.0 * mu * mu + 1.0);
    var taken = 0.0;
    for (var k = 0; k < K; k++) {
      path += max(0.0, IN[S_Q + k * C + i]) * pi * LV[L_DS + k] / GRAV * sqrt(LV[L_SM + k]) * 0.1 * magnification;
      let through = VAPOR_ABS * 2.9 * path / (pow(1.0 + 141.5 * path, 0.635) + 5.925 * path);
      netFlux[k] += incident * (through - taken);
      taken = through;
    }
    vaporHeating = incident * taken;
    incident -= vaporHeating;
  }
  var light = vec3<f32>(0.0, 0.0, 0.0); var aerosolHeating = 0.0;
  if (SCATTER || UPWARD) {
    let clear = clearLight(beam, mu, visibleTaken, incident, pi, aerosol);
    aerosolHeating = clear.x; light = clear.yzw;
    incident -= aerosolHeating;
    if (aerosolHeating > 0.0) { for (var k = 0; k < K; k++) { netFlux[k] += aerosolHeating * LV[L_AER + k]; } }
  }
  var columnCover = overlapCover(blocks, layered);
  if (!(columnCover > 0.0)) { columnCover = 1.0; }
  var sw = shortwave(select(cloudDepth / columnCover, CLOUD_SCAT * cloudPath / columnCover, GRAY_SW), cloudKeep(cloudPath / columnCover), mu, adir, adif, light);
  if (PDF_COVER && columnCover < 1.0) { sw = columnCover * sw + (1.0 - columnCover) * shortwave(0.0, 1.0, mu, adir, adif, light); }
  var clearShare = select(0.0, incident * sw.w / cloudPath, cloudPath > 0.0);
  var deckShare = 0.0;
  if (deck > 0.0) {
    let decked = shortwave(select(cloudDepth + deckUnit * deck, CLOUD_SCAT * (cloudPath + deck), GRAY_SW), cloudKeep(cloudPath + deck), mu, adir, adif, light);
    deckShare = fraction * incident * decked.w / (cloudPath + deck);
    clearShare *= 1.0 - fraction;
    sw = fraction * decked + (1.0 - fraction) * sw;
  }
  var cloudHeating = 0.0;
  if (CLOUD_SW > 0.0 && incident > 0.0 && sw.w > 0.0) {
    for (var k = 0; k < K; k++) {
      let mass = pi * LV[L_DS + k] / GRAV;
      let water = max(0.0, IN[S_QC + k * C + i]) * mass + cumulusCloud(k, i).y * mass;
      let share = clearShare * water + deckShare * (water + select(0.0, deck, k == STRATUS_K));
      netFlux[k] += share;
      cloudHeating += share;
    }
  }
  var upwardHeating = 0.0; var upwardAerosol = 0.0; var restLoss = 0.0; var visibleLoss = 0.0; var aerosolLoss = 0.0; var ozoneLoss = 0.0;
  if (UPWARD && mu > 0.0) {
    aerosolLoss = select(0.0, 1.0 - exp(-AER_ABS * aerosol * DIFFUSE_PATH), aerosol > 0.0);
    if (SOLAR_CLIRAD) { ozoneLoss = relaxedFraction(VIS_O3 * ozoneColumn * DIFFUSE_PATH); }
    visibleLoss = select(aerosolLoss, 1.0 - (1.0 - aerosolLoss) * (1.0 - ozoneLoss), ozoneLoss > 0.0);
    let sunlit = beam - ozoneHeating;
    let restAfter = sunlit - max(0.0, VIS_FRAC * beam - visibleTaken) - vaporHeating;
    var esc = escapes(select(cloudDepth / columnCover, CLOUD_SCAT * cloudPath / columnCover, GRAY_SW), cloudKeep(cloudPath / columnCover), mu, adir, adif, light);
    if (PDF_COVER && columnCover < 1.0) { esc = columnCover * esc + (1.0 - columnCover) * escapes(0.0, 1.0, mu, adir, adif, light); }
    if (deck > 0.0) { esc = fraction * escapes(select(cloudDepth + deckUnit * deck, CLOUD_SCAT * (cloudPath + deck), GRAY_SW), cloudKeep(cloudPath + deck), mu, adir, adif, light) + (1.0 - fraction) * esc; }
    if (SOLAR_CLIRAD && restAfter > 0.0) {
      var upGas: array<f32, K>;
      restLoss = beam * nearInfraredUpward(i, pi, solar.yzw, &upGas) / restAfter;
      let scale = select(1.0, 1.0 / restLoss, restLoss > 1.0);
      restLoss = min(restLoss, 1.0);
      let rest = incident * esc.y;
      for (var k = 0; k < K; k++) { netFlux[k] += rest * scale * beam * upGas[k] / restAfter; }
      upwardHeating = rest * restLoss;
    } else if (VAPOR_ABS > 0.0 && restAfter > 0.0) {
      let start = 2.9 * path / (pow(1.0 + 141.5 * path, 0.635) + 5.925 * path);
      var up = path;
      for (var k = 0; k < K; k++) { up += max(0.0, IN[S_Q + k * C + i]) * pi * LV[L_DS + k] / GRAV * sqrt(LV[L_SM + k]) * 0.1 * DIFFUSE_PATH; }
      restLoss = sunlit * VAPOR_ABS * (2.9 * up / (pow(1.0 + 141.5 * up, 0.635) + 5.925 * up) - start) / restAfter;
      let scale = select(1.0, 1.0 / restLoss, restLoss > 1.0);
      restLoss = min(restLoss, 1.0);
      let rest = incident * esc.y;
      var climbed = path; var before = start;
      for (var k = K - 1; k >= 0; k--) {
        climbed += max(0.0, IN[S_Q + k * C + i]) * pi * LV[L_DS + k] / GRAV * sqrt(LV[L_SM + k]) * 0.1 * DIFFUSE_PATH;
        let through = 2.9 * climbed / (pow(1.0 + 141.5 * climbed, 0.635) + 5.925 * climbed);
        netFlux[k] += rest * scale * sunlit * VAPOR_ABS * (through - before) / restAfter;
        before = through;
      }
      upwardHeating = rest * restLoss;
    }
    upwardAerosol = incident * esc.x * aerosolLoss;
    if (upwardAerosol > 0.0) { for (var k = 0; k < K; k++) { netFlux[k] += upwardAerosol * LV[L_AER + k]; } }
    upwardHeating += upwardAerosol;
    if (ozoneLoss > 0.0) {
      let upwardOzone = incident * (esc.x * (1.0 - aerosolLoss) + esc.z) * ozoneLoss;
      for (var k = 0; k < K; k++) { netFlux[k] += upwardOzone * layerOzone[k] / ozoneColumn; }
      upwardHeating += upwardOzone;
    }
  }
  let absorbed = incident * sw.x;
  var beforeBands: array<f32, K>;
  if (MOIST_BL${held ? ' || true' : ''}) { for (var k = 0; k < K; k++) { beforeBands[k] = netFlux[k]; } }
  if (LW_CORRELATED) {
    for (var k = 0; k < K; k++) { let idx = k * C + i; OUT[LWS_CLOUD + idx] = cloudE[k]; OUT[LWS_SHARE + idx] = chainF[k];${held ? '' : ' OUT[LWS_FLUX + idx] = netFlux[k];'} }
  } else {
    let v = band(VAPOR_FRAC, &vaporE, &temperature, &netFlux, surfaceEmission);
    let g = band(GAS_FRAC, &mixedE, &temperature, &netFlux, surfaceEmission);
    let w = band(WINDOW, &cloudE, &temperature, &netFlux, surfaceEmission);
${held ? HELD_GRAY : `    stash(20, i, v.y + g.y + w.y); stash(21, i, v.x + g.x + w.x);
    for (var k = 0; k < K; k++) { OUT[LWS_FLUX + k * C + i] = netFlux[k]; }
    if (MOIST_BL) { for (var k = 0; k < K; k++) { PH[PH_LWH + k * C + i] = netFlux[k] - beforeBands[k]; } }
`}  }
${held ? HELD_STORE : `  if (CLEAR_SKY) {
    let clearSw = shortwave(0.0, 1.0, mu, adir, adif, light);
    var clearUp = 0.0;
    if (UPWARD && mu > 0.0) { let clearEsc = escapes(0.0, 1.0, mu, adir, adif, light); clearUp = clearEsc.y * restLoss + clearEsc.x * aerosolLoss + select(0.0, (clearEsc.x * (1.0 - aerosolLoss) + clearEsc.z) * ozoneLoss, ozoneLoss > 0.0); }
    PH[PH_ABSCLRSUM + i] += incident * (clearSw.x + clearUp) + ozoneHeating + vaporHeating + aerosolHeating;
    if (!LW_CORRELATED) {
      var upVapor = VAPOR_FRAC * surfaceEmission; var upGas = GAS_FRAC * surfaceEmission;
      for (var k = K - 1; k >= 0; k--) {
        let t = temperature[k]; let ev = clearE[k]; let eg = LV[L_GASE + k];
        upVapor = upVapor * (1.0 - ev) + VAPOR_FRAC * ev * STEFAN * t * t * t * t;
        upGas = upGas * (1.0 - eg) + GAS_FRAC * eg * STEFAN * t * t * t * t;
      }
      PH[PH_OLRCLRSUM + i] += upVapor + upGas + WINDOW * surfaceEmission;
    }
  }
  let atmosphereSolar = ozoneHeating + vaporHeating + aerosolHeating + cloudHeating + upwardHeating; let absorbedSolar = absorbed + ozoneHeating + vaporHeating + aerosolHeating + cloudHeating + upwardHeating; let reflectedSolar = incident - absorbed - cloudHeating - upwardHeating;
  stash(14, i, absorbed); stash(15, i, incident * sw.y); stash(16, i, incident * sw.z); stash(17, i, atmosphereSolar); stash(18, i, absorbedSolar); stash(19, i, reflectedSolar);
`}  stash(0, i, 1.0);
}`;
export const longwaveKernel = (held, ng = DEFAULT_LONGWAVE_TABLE.points.length) => `${LONGWAVE_STASH}
@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= C || !LW_CORRELATED || !(stashed(0, i) > 0.0)) { return; }
  let pi = IN[S_PI + i]; let ts = stashed(9, i); let surfaceEmission = stashed(13, i);
  var back = 0.0; var outgoing = 0.0; var clearOutgoing = 0.0;${held ? ' var escape = 0.0; var clearEscape = 0.0;' : ''}
${longwaveSweeps(held, ng, {
    layer: 'let T = IN[S_TH + k * C + i] * D[D_EXM + k * C + i]; let row = longwavePaths(i, k, pi, T, OUT[LWS_OZONE + k * C + i]); let hot = STEFAN * T * T * T * T; let t = (T - 250.0) / 100.0;',
    chainF: (k) => `OUT[LWS_SHARE + ${k} * C + i]`, cloudE: 'OUT[LWS_CLOUD + k * C + i]', alpha: (k) => `OUT[LWS_ALPHA + ${k} * C + i]`,
    downEnd: held ? 'PH[PH_LWH + k * C + i] = PH[PH_RADSW + k * C + i] + heat;' : 'let b = OUT[LWS_FLUX + k * C + i]; if (MOIST_BL) { PH[PH_LWH + k * C + i] = b; } OUT[LWS_FLUX + k * C + i] = b + heat;',
    upEnd: held ? 'PH[PH_LWH + k * C + i] = PH[PH_LWH + k * C + i] + heat - PH[PH_RADSW + k * C + i];' : 'let n = OUT[LWS_FLUX + k * C + i] + heat; OUT[LWS_FLUX + k * C + i] = n; if (MOIST_BL) { PH[PH_LWH + k * C + i] = n - PH[PH_LWH + k * C + i]; }',
  })}
${held ? `  PH[PH_RADBACK + i] = back; PH[PH_RADOLR + i] = outgoing; PH[PH_RADT + i] = escape;
  if (CLEAR_SKY) { PH[PH_RADCLROLR + i] = clearOutgoing; PH[PH_RADTC + i] = clearEscape; }` : `  stash(20, i, back); stash(21, i, outgoing);
  if (CLEAR_SKY) { PH[PH_OLRCLRSUM + i] += clearOutgoing; }`}
}`;
export const physicsSurfaceKernel = (held) => `${EXCHANGE_WGSL}${LONGWAVE_STASH}
@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= C) { return; }
  let pi = IN[S_PI + i]; let skin = IN[S_TS + i]; let ice = IN[S_ICE + i]; let bottom = (K - 1) * C + i;
  let onLand = PH[PH_LAND + i] > 0.5; let onIceSheet = PH[PH_LAND + i] > 1.5;
  let soil0 = PH[PH_SOIL + i]; let snow0 = PH[PH_SNOW + i]; let veg0 = PH[PH_VEG + i]; let surf0 = PH[PH_SURF + i];
  let bucket = select(BUCKET, ROOTCAP, VEGETATED);
  let ws = stashed(1, i); let mu = stashed(2, i); let beam = stashed(3, i); let cover = stashed(4, i); let waterDir = stashed(5, i); let iceDif = stashed(6, i); let iceDir = stashed(7, i); let adif = stashed(8, i);
  let ts = stashed(9, i); let warmth = stashed(10, i); let wetness = stashed(11, i); let bareShare = stashed(12, i); let surfaceEmission = stashed(13, i);
${held ? HELD_REST : `  let absorbed = stashed(14, i); let swdn = stashed(15, i); let directDown = stashed(16, i); let atmosphereSolar = stashed(17, i); let absorbedSolar = stashed(18, i); let reflectedSolar = stashed(19, i);
  let back = stashed(20, i); let outgoing = stashed(21, i);
`}  let airT = IN[S_TH + bottom] * D[D_EXM + bottom];
  let rho = pi * LV[L_SM + K - 1] / (RGAS * airT);
  let exchange = rho * PH[PH_HEATX + i] * xWind(i, ws);
  let sensible = select(exchange * CP * (ts - airT), exchange * (CP * (ts - airT) - CP * D[D_THV + bottom] * (D[D_EXL + bottom] - D[D_EXM + bottom])), ROUGH);
  let evap = wetness * max(0.0, exchange * (qsat(ts, pi) - IN[S_Q + bottom]));
  if (ROUGH) { PH[PH_BUOY + i] = GRAV / IN[S_TH + bottom] * (sensible / (rho * CP * D[D_EXM + bottom]) + 0.61 * IN[S_TH + bottom] * evap / rho); }
  let airQs = qsat(airT, pi); let airSlope = airQs * 4302.645 / ((airT - 29.65) * (airT - 29.65)); let conductance = PH[PH_REFX + i] * xWind(i, ws);
  let potential = (airSlope * (absorbed - surfaceEmission + back) + rho * CP * conductance * (airQs - IN[S_Q + bottom])) / (LHEAT * airSlope + CP * (1.0 + REF_RESIST * conductance));
  let net = absorbed - surfaceEmission + back - sensible - LHEAT * evap;
  let dt = P[0];
${held ? HELD_THETA : `  for (var k = 0; k < K; k++) {
    let idx = k * C + i;
    let mass = pi * LV[L_DS + k] / GRAV;
    let flux = OUT[LWS_FLUX + idx];
    IN[S_TH + idx] += dt * select(flux, flux + sensible, k == K - 1) / (CP * mass) / D[D_EXM + idx];
  }
`}  IN[S_Q + bottom] += dt * evap * GRAV / (pi * LV[L_DS + K - 1]);
${held ? '  let swdn = sunScale * PH[PH_RADSWDN + i];\n' : ''}  PH[PH_SWDN + i] = swdn;
${held ? '  let directDown = sunScale * PH[PH_RADDIR + i];\n' : ''}  let contrast = directDown * (iceDir - waterDir) + (swdn - directDown) * (iceDif - ALB_DIF_WATER);
${held ? HELD_TOA : ''}  PH[PH_SFLUX + i] = net; PH[PH_ABS + i] = absorbedSolar; PH[PH_ATMSW + i] = atmosphereSolar; PH[PH_OLR + i] = outgoing; PH[PH_SH + i] = sensible; PH[PH_EVAP + i] = evap; PH[PH_INS + i] = beam; PH[PH_REFL + i] = reflectedSolar; PH[PH_ADIF + i] = adif;
  PH[PH_ABSSUM + i] += absorbedSolar; PH[PH_ATMSUM + i] += atmosphereSolar; PH[PH_OLRSUM + i] += outgoing; PH[PH_INSSUM + i] += beam; PH[PH_REFLSUM + i] += reflectedSolar; PH[PH_LWSFCSUM + i] += back - surfaceEmission;
${held ? HELD_CLEAR : ''}  let ocean = PH[PH_OFLUX + i]; let capacity = PH[PH_CAP + i];
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
      if (snow <= 0.0 || onIceSheet) { PH[PH_SNOWFREEV + i] = veg; }
      let air = select(IN[S_TH + bottom] * D[D_EXM + bottom], airT, ROUGH);
      let keep = P[5]; let hold = P[7];
      let seasonLength = PH[PH_SEASONL + i] + (select(0.0, 1.0, air >= SEASON_K) - PH[PH_SEASONL + i]) * keep;
      let seasonWarmth = PH[PH_SEASONW + i] + (max(0.0, air - SEASON_K) - PH[PH_SEASONW + i]) * keep;
      PH[PH_SEASONL + i] = seasonLength; PH[PH_SEASONW + i] = seasonWarmth;
      let demand = PH[PH_DEMAND + i] + (86400.0 * potential - PH[PH_DEMAND + i]) * P[6];
      PH[PH_DEMAND + i] = demand;
      let standing = PH[PH_CANOPY + i];
      var trees = max(veg, standing + (veg - standing) * (1.0 - exp(-dt / CANOPY_MEM)));
      if (TREELINE) {
        let moist = select(1.0, clamp((PH[PH_RAINMEAN + i] / max(demand, 1e-3) - ARID_LO) / ARID_SPAN, 0.0, 1.0), GATED);
        let factor = clamp((SEASON_C + seasonWarmth / max(seasonLength, SEASON_SHORTEST) - TREE_LO) / TREE_SPAN, 0.0, 1.0) * moist;
        let goal = select(factor * veg, min(standing, factor), snow > 0.0);
        trees = standing + (goal - standing) * smallRate(dt / select(TREE_DECLINE, TREE_GROW, goal > standing));
        if (hold > 0.5) { trees = select(select(factor * veg, 0.0, hold < 2.5), HELD_SHARE * veg, hold < 1.5); }
      }
      PH[PH_CANOPY + i] = select(trees, 0.0, onIceSheet);
      if (HUMIC) {
        let fill = clamp(soil / ROOTCAP, 0.0, 1.0); let airC = air - MELTING;
        let input = LITTER_IN * (LITTER_TREE * min(1.0, trees) + LITTER_GRASS * max(0.0, veg - trees)) * min(1.0, soil / (WETT * ROOTCAP)) * select(0.0, 1.0 / (1.0 + exp(MIAMI_A - MIAMI_B * airC)), airC >= SEASON_C);
        let moistDecay = select(select(select(1.0 - 0.8 * (fill - DECAY_OPT), 0.2 + 0.8 * (fill - DECAY_WILT) / (DECAY_OPT - DECAY_WILT), fill <= DECAY_OPT), 0.2, fill <= DECAY_WILT), 0.2, air < MELTING);
        let decayFactor = select(0.0, exp(LT_E * (LT_REF - 1.0 / max(air - LT_T0, 1e-3))), air > LT_T0) * moistDecay;
        let decay = decayFactor * DECAY_RATE;
        let litterFactor = min(1.0, soil / (WETT * ROOTCAP)) * select(0.0, 1.0 / (1.0 + exp(MIAMI_A - MIAMI_B * airC)), airC >= SEASON_C);
        PH[PH_LITTERM + i] = select(PH[PH_LITTERM + i] + (litterFactor - PH[PH_LITTERM + i]) * keep, 0.0, onIceSheet);
        PH[PH_DECAYM + i] = select(PH[PH_DECAYM + i] + (decayFactor - PH[PH_DECAYM + i]) * keep, 0.0, onIceSheet);
        let x = CARBON_ACC * decay * dt; let carbon0 = PH[PH_SOILC + i];
        let heldCarbon = select(select(HELD_GREEN_C, 0.0, hold < 2.5), HELD_NEUTRAL_C, hold < 1.5);
        PH[PH_SOILC + i] = select(select(max(0.0, carbon0 + (input - decay * carbon0) * CARBON_ACC * dt * select(1.0 - x * (0.5 - x / 6.0), (1.0 - exp(-x)) / x, x > 1e-2)), heldCarbon, hold > 0.5), 0.0, onIceSheet);
      }
      cap = ROOTCAP;
    }
    if (soil > cap) { PH[PH_RUNOFF + i] += soil - cap; soil = cap; }
    PH[PH_SOIL + i] = soil; PH[PH_SNOW + i] = snow; PH[PH_SURF + i] = surf;
    PH[PH_SNOWALB + i] = select(ALB_FRESH, agedSnow(PH[PH_SNOWALB + i], T, dt, ALB_OLDSNOW), snow > 0.0);
  } else ${SEA_SURFACE_WGSL}
  IN[S_TS + i] = T; IN[S_ICE + i] = h;
}`;

const ADJUST_FUNCTIONS = `fn upperInterface(i: i32, k: i32) -> f32 {
  let idx = k * C + i;
  return (D[D_GEO + idx] + LV[L_GABS + k] + CP * D[D_THV + idx] * (D[D_EXM + idx] - D[D_EXL + idx - C])) / GRAV;
}
fn longCloudShare(i: i32, k: i32, pi: f32, iced: f32) -> f32 {
  if (PH[PH_CUMF + i] > 0.0 && pi * LV[L_SM + k] >= PH[PH_CUTOP + i]) { return 0.0; }
  var share = 0.0;
  if (MOIST_BL) {
    if (!((D[D_GEO + k * C + i] + LV[L_GABS + k]) / GRAV < PH[PH_MIXTOP + i])) { share = PH[PH_STRAT + i]; }
    else if (PH[PH_REGIME + i] == 3.0) { share = 1.0; }
  }
  return max(share, iced);
}
fn saturateColumn(i: i32, pi: f32) {
  for (var k = 0; k < K; k++) {
    let idx = k * C + i;
    let ex = D[D_EXM + idx];
    let temperature = IN[S_TH + idx] * ex;
    let pressure = pi * LV[L_SM + k];
    var change = 0.0;
    if (UNIFORM && (BL_UNIFORM || (BL_CLOUDLAYER && f32(k) <= PH[PH_CLOUDK + i]) || !(MOIST_BL && (D[D_GEO + idx] + LV[L_GABS + k]) / GRAV < PH[PH_MIXTOP + i]))) {
      let water = IN[S_QC + idx];
      let liquidT = temperature - LHEAT * water / CP;
      let saturated = cloudSat(liquidT, pressure);
      let a = 1.0 / (1.0 + LHEAT * saturated.y / CP);
      let b = uniformWidth(saturated, pressure, pi);
      let total = IN[S_Q + idx] + water;
      let Q = a * (total - saturated.x);
      var held = max(0.0, Q);
      if (b > 0.0) { held = select(select(0.0, (Q + b) * (Q + b) / (4.0 * b), Q > -b), Q, Q >= b); }
      if (NUCLEATION && !(water > CLEAR_AIR) && liquidT < ICE_T) {
        let reference = min(qsat(liquidT, pressure), (2.583 - liquidT / 207.8) * saturated.x);
        let width = b / (a * saturated.x) * reference;
        let lowest = max(reference, total - width);
        held = select(0.0, max(0.0, a * (total + width - lowest) / (2.0 * width) * (0.5 * (lowest + total + width) - saturated.x)), total + width > reference);
      }
      change = held - water;
    } else {
      let saturated = cloudSat(temperature, pressure);
      change = (IN[S_Q + idx] - saturated.x) / (1.0 + LHEAT * saturated.y / CP);
    }
    if (change < 0.0) { change = max(change, -IN[S_QC + idx]); }
    if (change == 0.0) { continue; }
    IN[S_Q + idx] -= change; IN[S_QC + idx] += change; IN[S_TH + idx] += LHEAT * change / (CP * ex);
  }
}
// Never runs (dt > 0): two more call sites keep Metal from inlining saturateColumn, which in these kernels is slower and rounds differently.
fn saturateOutOfLine(i: i32, pi: f32, dt: f32) { if (dt < 0.0) { saturateColumn(i, pi); saturateColumn(i, pi + 1.0); } }
fn clearCumulus(i: i32) {
  for (var k = CU_K0; k < K; k++) { PH[PH_CUCOVER + (k - CU_K0) * C + i] = 0.0; PH[PH_CUWATER + (k - CU_K0) * C + i] = 0.0; }
  PH[PH_CUMF + i] = 0.0; PH[PH_CUTOP + i] = 0.0;
}
fn plumeState(energy: f32, water: f32, height: f32, pressure: f32, guess: f32) -> vec2<f32> {
  if (PL_MIXED) { return frozenState(energy, water, height, pressure, guess).xy; }
  let dry = (energy - GRAV * height) / CP;
  if (!(water > qsat(dry, pressure))) { return vec2<f32>(dry, 0.0); }
  let t = saturatedTemperature(energy - GRAV * height + LHEAT * water, pressure, max(dry, guess));
  return vec2<f32>(t, max(0.0, water - qsat(t, pressure)));
}
fn mixedSat(T: f32, p: f32) -> vec3<f32> {
  let alpha = clamp((T - ICE_T) / (LIQUID_T - ICE_T), 0.0, 1.0);
  var es = esat(T);
  if (alpha < 1.0) { es = alpha * es + (1.0 - alpha) * 611.21 * exp(22.587 * (T - 273.16) / (T + 0.7)); }
  let dry = p - (1.0 - EPSILON) * es;
  let qs = select(1.0, EPSILON * es / dry, dry > 0.0);
  return vec3<f32>(qs, qs * (LHEAT + (1.0 - alpha) * LFUSION) / (RVAP * T * T), alpha);
}
// moist.module.js's frozenPlumeState: temperature, condensate and its ice
fn frozenState(energy: f32, water: f32, height: f32, pressure: f32, guess: f32) -> vec3<f32> {
  let dry = (energy - GRAV * height) / CP;
  if (!(water > mixedSat(dry, pressure).x)) { return vec3<f32>(dry, 0.0, 0.0); }
  let goal = energy - GRAV * height + LHEAT * water; let span = 1.0 / (LIQUID_T - ICE_T);
  var t = max(dry, guess);
  for (var n = 0; n < 4; n++) {
    let m = mixedSat(t, pressure); let held = water - m.x; let frozen = 1.0 - m.z;
    t -= (CP * t + LHEAT * m.x - LFUSION * frozen * held - goal) / (CP + (LHEAT + LFUSION * frozen) * m.y + select(0.0, LFUSION * span * held, m.z > 0.0 && m.z < 1.0));
  }
  let m = mixedSat(t, pressure); let l = max(0.0, water - m.x);
  return vec3<f32>(t, l, (1.0 - m.z) * l);
}
// Dawn compiles with Metal's relaxed math, which fuses and reassociates unguarded products and sums with their consumers; the max, min and fma here fix each value's rounding wherever it is inlined
fn layerP(pi: f32, k: i32) -> f32 { return max(pi * LV[L_SM + k], 0.0); }
fn layerDp(pi: f32, k: i32) -> f32 { return max(pi * LV[L_DS + k], 0.0); }
fn plumeT(i: i32, k: i32) -> f32 { let idx = k * C + i; return max(IN[S_TH + idx] * D[D_EXM + idx], 0.0); }
fn plumeZ(i: i32, k: i32) -> f32 { return min((D[D_GEO + k * C + i] + LV[L_GABS + k]) / GRAV, 3.4e38); }
fn plumeEnv(i: i32, k: i32) -> vec2<f32> {
  let idx = k * C + i; let cloud = max(0.0, IN[S_QC + idx]);
  return vec2<f32>(fma(-LHEAT, cloud, fma(CP, IN[S_TH + idx] * D[D_EXM + idx], D[D_GEO + idx] + LV[L_GABS + k])), max(0.0, IN[S_Q + idx]) + cloud);
}
fn plumeVirtual(T: f32, air: f32, cloud: f32) -> f32 { return max(T * max(fma(PARCEL_VIRT, air, 1.0) - CU_LOADING * cloud, -3.4e38), -3.4e38); }
fn layerTop(i: i32, k: i32) -> f32 {
  let idx = k * C + i;
  return max(fma(max(CP * D[D_THV + idx], -3.4e38), max(D[D_EXM + idx] - D[D_EXL + idx - C], -3.4e38), max(D[D_GEO + idx] + LV[L_GABS + k], -3.4e38)) / GRAV, -3.4e38);
}
fn layerDepth(i: i32, k: i32) -> f32 { return max(layerTop(i, k) - layerTop(i, k + 1), -3.4e38); }
// the shallow cumulus mass flux of moist.module.js's cumulusColumn; returns its rain
fn cumulusColumn(i: i32, pi: f32, dt: f32) -> f32 {
  clearCumulus(i);
  cuSnow = 0.0;
  let open = select(1.0, clamp((DECK_CLOSED - PH[PH_MLMGATE + i]) / (DECK_CLOSED - DECK_OPEN), 0.0, 1.0), DECK_VETO) * select(1.0, 0.0, COUPLED_VETO && PH[PH_REGIME + i] == 3.0);
  let buoyancy = PH[PH_BUOY + i];
  if (!(open > 0.0) || !(buoyancy > 0.0)) { return 0.0; }
  let bottom = K - 1;
  let depth = PH[PH_DEPTH + i];
  var mass = 0.0; var energy = 0.0; var water = 0.0; var source = bottom;
  for (var k = bottom; k >= 0; k--) {
    let pk = layerP(pi, k);
    if (k < bottom && ((!(upperInterface(i, k + 1) < depth) && !(pk >= pi - CU_SOURCE)) || !(pk > SHALLOW_TOP))) { break; }
    let dpk = layerDp(pi, k); let env = plumeEnv(i, k);
    mass += dpk; energy += dpk * env.x; water += dpk * env.y; source = k;
  }
  var sourceS = energy / mass; var sourceQ = water / mass;
  if (CU_LOWEST) { let env = plumeEnv(i, bottom); sourceS = env.x; sourceQ = env.y; }
  if (!(sourceQ > 0.0)) { return 0.0; }
  let lcl = condensationLevel((sourceS - GRAV * plumeZ(i, bottom)) / CP, sourceQ, layerP(pi, bottom));
  if (!(lcl.y > SHALLOW_TOP)) { return 0.0; }
  var plumeS: array<f32, K + 1>; var plumeQ: array<f32, K + 1>; var liquid: array<f32, K>; var flux: array<f32, K + 1>; var fallout: array<f32, K>; var frozenOut: array<f32, K>;
  var s = sourceS; var w = sourceQ; var inhibition = 0.0; var cloudy = false; var top = -1; var guess = 0.0;
  plumeS[source] = s; plumeQ[source] = w;
  for (var k = source - 1; k >= 0; k--) {
    let pk = layerP(pi, k);
    if (!(pk > SHALLOW_TOP)) { if (cloudy) { top = k + 1; } break; }
    let below = upperInterface(i, k + 1); let above = upperInterface(i, k);
    let mixes = pi * LV[L_SL + k] <= lcl.y; let epsilon = select(0.0, CU_EPS, mixes);
    let zk = plumeZ(i, k); let env = plumeEnv(i, k);
    let half = exp(-epsilon * (zk - below));
    let midS = env.x + (s - env.x) * half; let midQ = env.y + (w - env.y) * half;
    let mid = plumeState(midS, midQ, zk, pk, guess);
    guess = mid.x; liquid[k] = mid.y;
    let idx = k * C + i; let air = max(0.0, IN[S_Q + idx]); let cloud = max(0.0, IN[S_QC + idx]);
    let work = RGAS * (mid.x * (1.0 + PARCEL_VIRT * (midQ - mid.y) - CU_LOADING * mid.y) - plumeT(i, k) * (1.0 + PARCEL_VIRT * air - CU_LOADING * cloud)) * layerDp(pi, k) / pk;
    if (mid.y > 0.0) { cloudy = true; }
    if (!cloudy) { if (work < 0.0) { inhibition -= work; } } else if (!(work > 0.0)) { top = k; break; }
    let full = exp(-epsilon * (above - below));
    s = env.x + (s - env.x) * full; w = env.y + (w - env.y) * full;
    if (CU_RAIN) {
      fallout[k] = 0.0; frozenOut[k] = 0.0;
      if (PL_MIXED) {
        let at = frozenState(s, w, above, pi * LV[L_SU + k], guess);
        let excess = at.y - CU_RAIN_Q;
        if (excess > 0.0) { w -= excess; s += LHEAT * excess; fallout[k] = excess; frozenOut[k] = excess * at.z / at.y; s += LFUSION * frozenOut[k]; }
      } else {
        let liquidAt = plumeState(s, w, above, pi * LV[L_SU + k], guess);
        let liquidExcess = liquidAt.y - CU_RAIN_Q;
        if (liquidExcess > 0.0) { w -= liquidExcess; s += LHEAT * liquidExcess; fallout[k] = liquidExcess; }
      }
    }
    plumeS[k] = s; plumeQ[k] = w;
    flux[k] = exp((epsilon - select(0.0, CU_DEL, mixes)) * (above - below));
  }
  if (top < 0) { return 0.0; }
  flux[top] = 0.0; plumeS[top] = 0.0; plumeQ[top] = 0.0;
  let lift = buoyancy * max(0.0, depth - plumeZ(i, bottom));
  let velocity = max(select(0.0, pow(lift, 1.0 / 3.0), lift > 0.0), CU_FRIC * PH[PH_USTAR + i]);
  if (!(velocity > 0.0)) { return 0.0; }
  var base = open * CU_C * lcl.y / (RGAS * lcl.x) * velocity * exp(-inhibition / (velocity * velocity));
  base = min(base, CU_LOSS * mass / (GRAV * dt));
  if (!(base > CU_FLOOR)) { return 0.0; }
  var sourceBelow = 0.0;
  for (var j = bottom; j > source; j--) { sourceBelow += layerDp(pi, j); flux[j] = base * sourceBelow / mass; }
  flux[source] = base;
  for (var k = source - 1; k > top; k--) { flux[k] = flux[k + 1] * flux[k]; }
  flux[top + 1] *= CU_OVER;
  var scale = 1.0;
  for (var k = top; k <= bottom; k++) {
    let courant = max(flux[k], flux[k + 1]) * GRAV * dt / layerDp(pi, k);
    if (courant * scale > 1.0) { scale = 1.0 / courant; }
  }
  for (var j = top + 1; j <= bottom; j++) { flux[j] *= scale; }
  var belowMass = 0.0; var belowS = 0.0; var belowQ = 0.0;
  for (var j = bottom; j > top; j--) {
    var upS = plumeS[j]; var upQ = plumeQ[j];
    if (j > source) {
      let dpj = layerDp(pi, j); let env = plumeEnv(i, j);
      belowMass += dpj; belowS += dpj * env.x; belowQ += dpj * env.y;
      upS = select(belowS / belowMass, sourceS, CU_LOWEST); upQ = select(belowQ / belowMass, sourceQ, CU_LOWEST);
    }
    let above = plumeEnv(i, j - 1);
    plumeS[j] = flux[j] * (upS - above.x); plumeQ[j] = flux[j] * (upQ - above.y);
  }
  for (var k = max(top, CU_K0); k < source; k++) {
    if (!(liquid[k] > 0.0)) { continue; }
    let slot = (k - CU_K0) * C + i;
    PH[PH_CUCOVER + slot] = min(1.0, 0.5 * (flux[k] + flux[k + 1]) * RGAS * plumeT(i, k) / (layerP(pi, k) * CU_WU));
    PH[PH_CUWATER + slot] = liquid[k];
  }
  var rain = 0.0;
  for (var k = top; k <= bottom; k++) {
    let idx = k * C + i; let dpk = layerDp(pi, k); let per = GRAV * dt / dpk;
    var dS = (plumeS[k + 1] - plumeS[k]) * per; var dQ = (plumeQ[k + 1] - plumeQ[k]) * per;
    if (CU_RAIN && k > top && k < source && fallout[k] > 0.0) {
      let fallen = flux[k] * fallout[k] * dt;
      rain += fallen; dQ -= fallen * GRAV / dpk; dS += LHEAT * fallen * GRAV / dpk;
      if (PL_MIXED && frozenOut[k] > 0.0) { let frozen = flux[k] * frozenOut[k] * dt; cuSnow += frozen; dS += LFUSION * frozen * GRAV / dpk; }
    }
    IN[S_TH + idx] += dS / (CP * D[D_EXM + idx]); IN[S_Q + idx] += dQ;
  }
  PH[PH_CUMF + i] = base * scale; PH[PH_CUTOP + i] = pi * LV[L_SU + top];
  return rain;
}
var<private> cuFall: array<f32, K>;
var<private> cuReserve: array<f32, K>;
var<private> cuDeep: bool;
var<private> cuBase: i32;
var<private> cuShallowRain: f32;
var<private> cuFrozen: array<f32, K>;
var<private> cuSnow: f32;
// moist.module.js's plumeColumn with the cumulus cloud's memory of rememberedPlume
fn plumeColumn(i: i32, pi: f32, dt: f32) -> f32 {
  if (!(CU_MEMORY > 0.0)) { return diagnosedPlume(i, pi, dt); }
  var cover: array<f32, K - CU_K0 + 1>; var path: array<f32, K - CU_K0 + 1>;
  for (var k = CU_K0; k < K; k++) { let slot = (k - CU_K0) * C + i; cover[k - CU_K0] = PH[PH_CUCOVER + slot]; path[k - CU_K0] = PH[PH_CUCOVER + slot] * PH[PH_CUWATER + slot]; }
  let produced = diagnosedPlume(i, pi, dt);
  let keep = exp(-dt / CU_MEMORY);
  for (var k = CU_K0; k < K; k++) {
    let slot = (k - CU_K0) * C + i; let now = PH[PH_CUCOVER + slot]; let made = now * PH[PH_CUWATER + slot];
    let kept = now + (cover[k - CU_K0] - now) * keep;
    if (!(kept >= CU_TRACE)) { PH[PH_CUCOVER + slot] = 0.0; PH[PH_CUWATER + slot] = 0.0; continue; }
    PH[PH_CUCOVER + slot] = kept; PH[PH_CUWATER + slot] = (made + (path[k - CU_K0] - made) * keep) / kept;
  }
  return produced;
}
// the convective mass flux of moist.module.js's plumeColumn; returns the rain it leaves falling, per layer in cuFall
fn diagnosedPlume(i: i32, pi: f32, dt: f32) -> f32 {
  cuDeep = false; cuBase = K; cuShallowRain = 0.0; cuSnow = 0.0;
  if (PL_MOMENTUM) {
    PH[PH_MOMS + i] = f32(K);
    for (var k = 0; k <= K; k++) { PH[PH_MOMU + k * C + i] = 0.0; PH[PH_MOMD + k * C + i] = 0.0; }
    for (var k = 0; k < K; k++) { PH[PH_MOMK + k * C + i] = 1.0; PH[PH_MOMKD + k * C + i] = 1.0; }
  }
  let open = select(1.0, clamp((DECK_CLOSED - PH[PH_MLMGATE + i]) / (DECK_CLOSED - DECK_OPEN), 0.0, 1.0), DECK_VETO) * select(1.0, 0.0, COUPLED_VETO && PH[PH_REGIME + i] == 3.0);
  if (!(open > 0.0)) { return cumulusColumn(i, pi, dt); }
  let bottom = K - 1;
  let depthTop = PH[PH_DEPTH + i];
  var mass = 0.0; var energy = 0.0; var water = 0.0; var source = bottom;
  for (var k = bottom; k >= 0; k--) {
    let pk = layerP(pi, k);
    if (k < bottom && (((PL_SURFACE50 || !(upperInterface(i, k + 1) < depthTop)) && !(pk >= pi - CU_SOURCE)) || !(pk > SHALLOW_TOP))) { break; }
    let dpk = layerDp(pi, k); let env = plumeEnv(i, k);
    mass += dpk; energy += dpk * env.x; water += dpk * env.y; source = k;
  }
  var sourceS = energy / mass; var sourceQ = water / mass;
  if (PL_LOWEST) { let env = plumeEnv(i, bottom); sourceS = env.x; sourceQ = env.y; }
  let surfaceT = plumeT(i, bottom); let surfaceP = layerP(pi, bottom); let surfaceZ = plumeZ(i, bottom);
  if (PL_SURFACE50) {
    let density = surfaceP / (RGAS * surfaceT);
    var velocity = 0.0;
    if (EX_LAYER) {
      let b = bottom * C + i; let height = CP * D[D_THV + b] * (D[D_EXL + b] - D[D_EXM + b]) / GRAV;
      let flux = (PH[PH_SH + i] / CP + EX_VIRT * surfaceT * PH[PH_EVAP + i]) / density;
      if (flux > 0.0) { velocity = EX_SCALE * pow(EX_USTAR * EX_USTAR * EX_USTAR + EX_STAB * GRAV * height * EX_KARMAN / surfaceT * flux, 1.0 / 3.0); }
    } else {
      let lift = PH[PH_BUOY + i] * max(0.0, depthTop - surfaceZ);
      velocity = max(select(0.0, pow(lift, 1.0 / 3.0), lift > 0.0), PH[PH_USTAR + i]);
    }
    if (velocity > 0.0) {
      var dT = min(EX_T, EX_COEF * PH[PH_SH + i] / (density * CP * velocity)); var dQ = min(EX_Q, EX_COEF * PH[PH_EVAP + i] / (density * velocity));
      if (EX_LAYER) { dT = max(0.0, dT); dQ = max(0.0, dQ); }
      sourceS += CP * dT; sourceQ += dQ;
    }
  }
  if (!(sourceQ > 0.0)) { return cumulusColumn(i, pi, dt); }
  let lcl = condensationLevel((sourceS - GRAV * surfaceZ) / CP, sourceQ, surfaceP);
  if (!(lcl.y > pi * LV[L_SL + 0])) { return cumulusColumn(i, pi, dt); }
  if (PL_PARCEL) {
    let surfaceSat = cloudSat(surfaceT, surfaceP).x;
    var ts = sourceS; var tq = sourceQ; var tw2 = PL_W0 * PL_W0; var tguess = 0.0; var tbase = 0.0; var typed = false; var decided = false;
    for (var k = source - 1; k > 0; k--) {
      let lower = upperInterface(i, k + 1); let upper = upperInterface(i, k); let depth = upper - lower;
      let Tk = plumeT(i, k); let pk = layerP(pi, k); let zk = plumeZ(i, k); let env = plumeEnv(i, k);
      let ratio = cloudSat(Tk, pk).x / surfaceSat;
      let epsilon = TP_EPS * IFS_EPS * ratio * ratio * ratio;
      let half = exp(-epsilon * (zk - lower));
      let midQ = env.y + (tq - env.y) * half;
      let mid = plumeState(env.x + (ts - env.x) * half, midQ, zk, pk, tguess);
      tguess = mid.x;
      if (!(tbase > 0.0) && mid.y > 0.0) { tbase = pi * LV[L_SL + k]; }
      let idx = k * C + i; let environment = Tk * (1.0 + PARCEL_VIRT * max(0.0, IN[S_Q + idx]) - CU_LOADING * max(0.0, IN[S_QC + idx]));
      let buoyancy = GRAV * (mid.x * (1.0 + PARCEL_VIRT * (midQ - mid.y) - CU_LOADING * mid.y) - environment) / environment;
      let mixing = IFS_DRAG * epsilon; let x = 2.0 * mixing * depth;
      let next = tw2 * exp(-x) + 2.0 * PL_ACC * buoyancy * depth * select(1.0, relaxedFraction(x) / x, x > 0.0);
      if (!(next > 0.0)) {
        if (tbase > 0.0) {
          let lift = -PL_ACC * buoyancy; let y = mixing * tw2 / lift;
          let reach = select(tw2 / (2.0 * lift), select(y * (1.0 - 0.5 * y), log(1.0 + y), y > 1e-3) / (2.0 * mixing), mixing > 0.0);
          typed = tbase - pi * LV[L_SL + k] * pow(LV[L_SU + k] / LV[L_SL + k], min(1.0, reach / depth)) > DEEP_DEPTH;
        }
        decided = true; break;
      }
      tw2 = next;
      if (tbase > 0.0 && tbase - pi * LV[L_SU + k] > DEEP_DEPTH) { typed = true; decided = true; break; }
      let full = exp(-epsilon * depth);
      ts = env.x + (ts - env.x) * full; tq = env.y + (tq - env.y) * full;
      if (PL_MIXED) { let at = frozenState(ts, tq, upper, pi * LV[L_SU + k], tguess); tq -= TP_REMOVE * at.y; ts += LHEAT * TP_REMOVE * at.y + LFUSION * TP_REMOVE * at.z; }
      else { let liquidAt = plumeState(ts, tq, upper, pi * LV[L_SU + k], tguess); tq -= TP_REMOVE * liquidAt.y; ts += LHEAT * TP_REMOVE * liquidAt.y; }
    }
    if (!decided) { typed = tbase > 0.0 && tbase - pi * LV[L_SL + 0] > DEEP_DEPTH; }
    if (!typed) { return cumulusColumn(i, pi, dt); }
  }
  var plumeS: array<f32, K>; var plumeQ: array<f32, K>; var liquid: array<f32, K>; var speed: array<f32, K + 1>; var flux: array<f32, K + 1>;
  var entrained: array<f32, K>; var thick: array<f32, K>; var counted: array<f32, K>;
  var sourceBelow = 0.0;
  for (var j = bottom; j > source; j--) { sourceBelow += layerDp(pi, j); flux[j] = sourceBelow / mass; }
  flux[source] = 1.0;
  var s = sourceS; var w = sourceQ; var w2 = 0.0; var below = 0.0; var inhibition = 0.0; var cloudy = false; var started = false; var top = -1; var guess = 0.0; var cape = 0.0; var base = -1; var neutral = -1; var neutralB = 0.0; var aboveB = 0.0; var pcape = 0.0; var baseSat = 0.0;
  var rainAbove = 0.0; var rainOne = 0.0;
  plumeS[source] = s; plumeQ[source] = w;
  for (var k = source - 1; k >= 0; k--) {
    let lower = upperInterface(i, k + 1);
    var upper = lower;
    if (k > 0) { upper = upperInterface(i, k); }
    let depth = upper - lower;
    let Tk = plumeT(i, k); let pk = layerP(pi, k); let dpk = layerDp(pi, k); let zk = plumeZ(i, k); let env = plumeEnv(i, k);
    let mixes = pi * LV[L_SL + k] <= lcl.y;
    if (mixes && !started) { started = true; base = k + 1; w2 = PL_W0 * PL_W0; speed[k + 1] = w2; if (PL_IFS) { baseSat = cloudSat(Tk, pk).x; } }
    var epsilon = 0.0; var mixing = 0.0; var rh = 0.0; var detrained = 0.0;
    if (PL_IFS) {
      let qsk = cloudSat(Tk, pk).x;
      rh = min(1.0, max(0.0, IN[S_Q + k * C + i]) / qsk);
      detrained = select(0.0, IFS_DEL * (IFS_DRH - rh), mixes);
      if (mixes && below > 0.0) { let ratio = qsk / baseSat; epsilon = IFS_EPS * (IFS_RH - rh) * ratio * ratio * ratio; }
      mixing = IFS_DRAG * select(detrained, epsilon, epsilon > 0.0);
    } else {
      if (mixes) { epsilon = max(PL_FLOOR, PL_EPS * max(0.0, below) / w2); }
      mixing = PL_DRAG * epsilon;
    }
    let half = exp(-epsilon * (zk - lower));
    let midS = env.x + (s - env.x) * half; let midQ = env.y + (w - env.y) * half;
    let mid = plumeState(midS, midQ, zk, pk, guess);
    guess = mid.x; liquid[k] = mid.y;
    let idx = k * C + i; let air = max(0.0, IN[S_Q + idx]); let cloud = max(0.0, IN[S_QC + idx]);
    let environment = Tk * (1.0 + PARCEL_VIRT * air - CU_LOADING * cloud); let rising = mid.x * (1.0 + PARCEL_VIRT * (midQ - mid.y) - CU_LOADING * mid.y);
    let wk = RGAS * (rising - environment) * dpk / pk; let buoyancy = GRAV * (rising - environment) / environment;
    if (mid.y > 0.0) { cloudy = true; }
    if (!cloudy && wk < 0.0) { inhibition -= wk; }
    var work = wk;
    if (!PL_IFS || PL_MOMENTUM) { entrained[k] = epsilon; thick[k] = depth; }
    if (PL_UNDILUTE) {
      let parcel = plumeState(sourceS, sourceQ, zk, pk, mid.x);
      work = RGAS * (parcel.x * (1.0 + PARCEL_VIRT * (sourceQ - parcel.y)) - environment) * dpk / pk;
    }
    if (mixes) {
      if (k == 0) { top = k; if (neutral == k + 1) { aboveB = buoyancy; } break; }
      let x = 2.0 * mixing * depth;
      w2 = w2 * exp(-x) + 2.0 * PL_ACC * buoyancy * depth * select(1.0, relaxedFraction(x) / x, x > 0.0);
      if (!(w2 > 0.0)) { top = k; if (neutral == k + 1) { aboveB = buoyancy; } break; }
    }
    speed[k] = select(0.0, w2, mixes);
    if (buoyancy > 0.0) { neutral = k; neutralB = buoyancy; aboveB = 0.0; } else if (neutral == k + 1) { aboveB = buoyancy; }
    if (cloudy && work > 0.0) { cape += work; if (PL_BUOYANT_F) { counted[k] = 1.0; } if (PL_BECHTOLD) { pcape += work * pk / (RGAS * environment); } }
    below = buoyancy;
    let full = exp(-epsilon * depth);
    s = env.x + (s - env.x) * full; w = env.y + (w - env.y) * full;
    var fallout = 0.0;
    if (PL_MIXED) { cuFrozen[k] = 0.0; }
    if (mixes) {
      var at = vec3<f32>(0.0, 0.0, 0.0);
      if (PL_MIXED) { at = frozenState(s, w, upper, pi * LV[L_SU + k], guess); } else { let liquidAt = plumeState(s, w, upper, pi * LV[L_SU + k], guess); at = vec3<f32>(liquidAt, 0.0); }
      var fallen = 0.0;
      if (PL_SUNDQVIST) {
        if (at.y > select(SQ_SEA, SQ_LAND, PH[PH_LAND + i] > 0.5)) {
          let alpha = clamp((at.x - ICE_T) / (LIQUID_T - ICE_T), 0.0, 1.0);
          var bergeron = 1.0;
          if (at.x < SQ_BF) { bergeron = 1.0 + 0.5 * sqrt(min(SQ_BF - at.x, SQ_BF - SQ_ICE)); }
          let ratio = at.y * bergeron / SQ_CRIT;
          let rate = SQ_C00 * (SQ_LIQ * alpha + 1.0 - alpha) * bergeron / (SQ_VSCALE * min(SQ_SPEED, max(PL_W0, sqrt(w2)))) * relaxedFraction(ratio * ratio);
          fallen = at.y * relaxedFraction(rate * depth);
        }
      } else {
        let excess = at.y - PL_RAIN_Q;
        if (excess > 0.0) { fallen = excess * relaxedFraction(PL_RAIN_RATE * depth); }
      }
      if (fallen > 0.0) {
        w -= fallen; s += LHEAT * fallen; fallout = fallen;
        if (PL_MIXED && at.z > 0.0) { let frozenOut = fallen * at.z / at.y; cuFrozen[k] = frozenOut; s += LFUSION * frozenOut; }
      }
    }
    cuFall[k] = fallout;
    plumeS[k] = s; plumeQ[k] = w;
    if (PL_IFS) {
      if (k >= base) { flux[k] = flux[k + 1]; }
      else if (buoyancy > 0.0) { flux[k] = flux[k + 1] * exp((epsilon - detrained) * depth); }
      else { flux[k] = flux[k + 1] * exp(-detrained * depth) * min(1.0, (IFS_DRH - rh) * sqrt(speed[k] / speed[k + 1])); }
      if (k > 1) { rainAbove += flux[k] * fallout; } else { rainOne = flux[k] * fallout; }
    }
  }
  if (top == 0) { top = 1; flux[1] = 0.0; } else { rainAbove += rainOne; }
  if (!cloudy || top < 1 || !(PL_PARCEL || select((LV[L_SU + top] * DEEP_REFERENCE < SHALLOW_TOP), (pi * (LV[L_SU + base] - LV[L_SU + top]) > DEEP_DEPTH), PL_BY_DEPTH))) { return cumulusColumn(i, pi, dt); }
  clearCumulus(i);
  let topHeight = upperInterface(i, top);
  var neutralHeight = upperInterface(i, top + 1);
  if (neutral > top && !(aboveB > 0.0)) { let zn = plumeZ(i, neutral); neutralHeight = min(neutralHeight, zn + (plumeZ(i, neutral - 1) - zn) * neutralB / (neutralB - aboveB)); }
  if (!PL_IFS) {
    var kn = source - 1;
    loop {
      if (!(kn > top) || upperInterface(i, kn) > neutralHeight) { break; }
      flux[kn] = flux[kn + 1] * exp((entrained[kn] - max(0.0, entrained[kn] - PL_GROWTH)) * thick[kn]);
      kn--;
    }
    let anchor = flux[kn + 1];
    for (; kn > top; kn--) { flux[kn] = anchor * (topHeight - upperInterface(i, kn)) / (topHeight - neutralHeight); }
    for (var k = top + 1; k < source; k++) { rainAbove += flux[k] * cuFall[k]; }
  }
  for (var k = 0; k < top; k++) { cuFall[k] = 0.0; if (PL_MIXED) { cuFrozen[k] = 0.0; } }
  var dflux: array<f32, K + 1>; var dS: array<f32, K + 1>; var dQ: array<f32, K + 1>; var devap: array<f32, K>;
  var start = -1; var share = 0.0; var lowest = 0.0;
  if (DD_SHARE > 0.0 && rainAbove > 0.0) {
    for (var k = top + 1; k < base; k++) {
      let env = plumeEnv(i, k); let h = env.x + LHEAT * env.y;
      if (start < 0 || h < lowest) { start = k; lowest = h; }
    }
  }
  if (start >= 0) {
    let env = plumeEnv(i, start);
    var hd = env.x + LHEAT * env.y; var fd = 1.0; var subcloud = 0.0;
    for (var k = base + 1; k <= bottom; k++) { subcloud += layerDp(pi, k); }
    let td = saturatedTemperature(hd - GRAV * upperInterface(i, start + 1), pi * LV[L_SL + start], plumeT(i, start));
    let qd = qsat(td, pi * LV[L_SL + start]);
    dflux[start + 1] = -1.0; dS[start + 1] = hd - LHEAT * qd; dQ[start + 1] = qd;
    devap[start] = qd - env.y;
    var wd = qd; var left = subcloud; var atBase = 1.0;
    for (var k = start + 1; k < bottom; k++) {
      var mixedQ = wd;
      if (k <= base) {
        let env = plumeEnv(i, k);
        let keep = exp(-DD_EPS * (upperInterface(i, k) - upperInterface(i, k + 1))); let ambient = env.x + LHEAT * env.y;
        fd /= keep;
        hd = ambient + (hd - ambient) * keep;
        mixedQ = env.y + (wd - env.y) * keep;
        atBase = fd;
      } else {
        left -= layerDp(pi, k);
        fd = atBase * left / subcloud;
      }
      let pk = pi * LV[L_SL + k];
      let tk = saturatedTemperature(hd - GRAV * upperInterface(i, k + 1), pk, plumeT(i, k));
      let qk = qsat(tk, pk);
      dflux[k + 1] = -fd; dS[k + 1] = hd - LHEAT * qk; dQ[k + 1] = qk;
      devap[k] = (qk - mixedQ) * fd;
      wd = qk;
    }
    share = DD_SHARE;
  }
  var flying = 0.0; var produced = 0.0; var taken = 0.0;
  for (var k = top; k <= bottom; k++) {
    var melted = 0.0; var frozenOut = 0.0;
    if (PL_MIXED) {
      if (flying > 0.0 && plumeT(i, k) >= LIQUID_T) { melted = flying; flying = 0.0; }
      frozenOut = cuFrozen[k];
      cuFrozen[k] = -melted;
      if (k > top && k < source && frozenOut > 0.0) { flying += flux[k] * frozenOut; cuFrozen[k] += flux[k] * frozenOut; }
    }
    if (start < 0 || k >= bottom) { continue; }
    if (k > top && k < source) { produced += select(flux[k] * cuFall[k], flux[k] * (cuFall[k] - frozenOut), PL_MIXED); }
    if (PL_MIXED && k >= top) { produced += melted; }
    taken += devap[k];
    if (taken > 0.0 && share * taken > produced) { share = produced / taken; }
  }
  if (start >= 0 && !(share > 0.0)) { share = 0.0; }
  var belowMass = 0.0; var belowS = 0.0; var belowQ = 0.0;
  var envJ = plumeEnv(i, bottom);
  for (var j = bottom; j > top; j--) {
    var upS = plumeS[j]; var upQ = plumeQ[j];
    if (j > source) {
      let dpj = layerDp(pi, j);
      belowMass += dpj; belowS += dpj * envJ.x; belowQ += dpj * envJ.y;
      upS = select(belowS / belowMass, sourceS, PL_LOWEST); upQ = select(belowQ / belowMass, sourceQ, PL_LOWEST);
    }
    let envUp = plumeEnv(i, j - 1);
    dS[j] = flux[j] * (upS - envUp.x) + share * dflux[j] * (dS[j] - envJ.x);
    dQ[j] = flux[j] * (upQ - envUp.y) + share * dflux[j] * (dQ[j] - envJ.y);
    envJ = envUp;
  }
  var consumption = 0.0; var consumptionP = 0.0;
  for (var k = top; k <= bottom; k++) {
    let dpk = layerDp(pi, k); let per = GRAV / dpk;
    var made = 0.0;
    if (k > top && k < source) { made = flux[k] * cuFall[k]; }
    let evaporated = share * devap[k];
    if (k < source && k > top && (!PL_BUOYANT_F || counted[k] > 0.0)) {
      var tS = (dS[k + 1] - dS[k] + LHEAT * (made - evaporated)) * per;
      if (PL_MIXED) { tS = (dS[k + 1] - dS[k] + LHEAT * (made - evaporated) + LFUSION * cuFrozen[k]) * per; } let tQ = (dQ[k + 1] - dQ[k] - made + evaporated) * per;
      let idx = k * C + i; let air = max(0.0, IN[S_Q + idx]); let cloud = max(0.0, IN[S_QC + idx]);
      let Tk = plumeT(i, k);
      let warming = tS / CP * (1.0 + PARCEL_VIRT * air - CU_LOADING * cloud) + PARCEL_VIRT * Tk * tQ;
      consumption += RGAS * warming * dpk / layerP(pi, k);
      if (PL_BECHTOLD) { consumptionP += warming * dpk / plumeVirtual(Tk, air, cloud); }
    }
  }
  var relaxed = 0.0;
  if (PL_BECHTOLD) {
    let baseHeight = upperInterface(i, base); let cloudDepth = topHeight - baseHeight;
    var weighted = 0.0; var thickness = 0.0;
    for (var k = top; k < base; k++) { let depth = layerDepth(i, k); weighted += depth * sqrt(max(0.0, 0.5 * (speed[k] + speed[k + 1]))); thickness += depth; }
    let turnover = cloudDepth / (weighted / thickness);
    let tau = clamp((1.0 + BT_SCALE * sqrt(MF[F_AREA + i])) * turnover, BT_SHORT, BT_LONG);
    let b = bottom * C + i; let ground = (D[D_GEO + b] + LV[L_GABS + bottom] - CP * D[D_THV + b] * (D[D_EXL + b] - D[D_EXM + b])) / GRAV;
    var forcing = 0.0; var layerMass = 0.0; var wind = 0.0;
    for (var k = max(base, K - SUB_K); k <= bottom; k++) {
      let idx = k * C + i; let saved = PH[PH_SUBTV + (k - K + SUB_K) * C + i]; let dpk = layerDp(pi, k);
      if (saved > 0.0) { forcing += (plumeT(i, k) * (1.0 + PARCEL_VIRT * max(0.0, IN[S_Q + idx]) - CU_LOADING * max(0.0, IN[S_QC + idx])) - saved) / dt * dpk; }
      wind += dpk * length(cellWind(i, k)); layerMass += dpk;
    }
    let boundaryTime = select((baseHeight - ground) / max(BT_WIND, wind / layerMass), turnover, PH[PH_LAND + i] > 0.5);
    let pcapeBoundary = boundaryTime / BT_TSTAR * select(forcing, max(0.0, forcing), BT_POSITIVE);
    if (consumptionP > 0.0) { relaxed = max(0.0, pcape - pcapeBoundary) / (tau * consumptionP); }
  } else if (consumption > 0.0 && cape > PL_CAPE) { relaxed = (cape - PL_CAPE) / (PL_TAU * consumption); }
  let gate = clamp(0.5 + (CIN_MAX - inhibition) / max(1.0, CIN_MAX), 0.0, 1.0);
  let buoyancyFlux = PH[PH_BUOY + i];
  var shallowBase = 0.0;
  if (buoyancyFlux > 0.0 && !PL_RELAXED) {
    let lift = buoyancyFlux * max(0.0, depthTop - surfaceZ);
    let velocity = max(select(0.0, pow(lift, 1.0 / 3.0), lift > 0.0), CU_FRIC * PH[PH_USTAR + i]);
    if (velocity > 0.0) { shallowBase = CU_C * lcl.y / (RGAS * lcl.x) * velocity * exp(-inhibition / (velocity * velocity)); }
  }
  var baseFlux = min(open * max(shallowBase, gate * relaxed), CU_LOSS * mass / (GRAV * dt));
  for (var k = top; k <= bottom; k++) {
    let courant = baseFlux * (max(flux[k], flux[k + 1]) + share * max(-dflux[k], -dflux[k + 1])) * GRAV * dt / layerDp(pi, k);
    if (courant > 1.0) { baseFlux /= courant; }
  }
  if (!(baseFlux > CU_FLOOR)) {
    if (PL_RELAXED) { return cumulusColumn(i, pi, dt); }
    return 0.0;
  }
  cuDeep = true; cuBase = base;
  var held: array<f32, K>;
  for (var k = max(top, CU_K0); k < source; k++) {
    if (!(liquid[k] > 0.0)) { continue; }
    let up = max(PL_W0, sqrt(max(0.0, 0.5 * (speed[k] + speed[k + 1]))));
    let cover = min(1.0, 0.5 * (flux[k] + flux[k + 1]) * baseFlux * RGAS * plumeT(i, k) / (layerP(pi, k) * up));
    if (PL_SEPARATE) { held[k] = cover; continue; }
    let slot = (k - CU_K0) * C + i;
    if (cover > PH[PH_CUCOVER + slot]) { PH[PH_CUCOVER + slot] = cover; PH[PH_CUWATER + slot] = liquid[k]; }
  }
  var fallen = 0.0;
  for (var k = top; k <= bottom; k++) {
    let per = GRAV / layerDp(pi, k);
    var made = 0.0;
    if (k > top && k < source) { made = flux[k] * cuFall[k]; }
    let evaporated = share * devap[k];
    let idx = k * C + i;
    if (PL_MIXED) { IN[S_TH + idx] += baseFlux * dt * (dS[k + 1] - dS[k] + LHEAT * (made - evaporated) + LFUSION * cuFrozen[k]) * per / (CP * D[D_EXM + idx]); cuFrozen[k] *= baseFlux * dt; }
    else { IN[S_TH + idx] += baseFlux * dt * (dS[k + 1] - dS[k] + LHEAT * (made - evaporated)) * per / (CP * D[D_EXM + idx]); }
    IN[S_Q + idx] += baseFlux * dt * (dQ[k + 1] - dQ[k] - made + evaporated) * per;
    cuFall[k] = (made - evaporated) * baseFlux * dt;
    fallen += cuFall[k];
  }
  for (var k = bottom - 1; k >= 0; k--) { cuReserve[k] = max(0.0, cuReserve[k + 1] - cuFall[k + 1]); }
  if (PL_MOMENTUM) {
    for (var k = 0; k <= K; k++) { PH[PH_MOMU + k * C + i] = baseFlux * flux[k]; PH[PH_MOMD + k * C + i] = baseFlux * share * dflux[k]; }
    for (var k = 0; k < K; k++) {
      var keep = 1.0;
      if (k > top && k < source) { keep = exp(-entrained[k] * thick[k]); }
      PH[PH_MOMK + k * C + i] = keep;
      var keepD = 1.0;
      if (k == start) { keepD = 0.0; } else if (start >= 0 && k > start && k <= base) { keepD = exp(-DD_EPS * (upperInterface(i, k) - upperInterface(i, k + 1))); }
      PH[PH_MOMKD + k * C + i] = keepD;
    }
    PH[PH_MOMS + i] = f32(source);
  }
  var shallowRain = 0.0;
  cuSnow = 0.0;
  if (PL_SEPARATE) {
    shallowRain = cumulusColumn(i, pi, dt);
    for (var k = max(top, CU_K0); k < source; k++) {
      if (!(liquid[k] > 0.0)) { continue; }
      let slot = (k - CU_K0) * C + i;
      if (held[k] > PH[PH_CUCOVER + slot]) { PH[PH_CUCOVER + slot] = held[k]; PH[PH_CUWATER + slot] = liquid[k]; }
    }
  }
  PH[PH_CUMF + i] += baseFlux; PH[PH_CUTOP + i] = pi * LV[L_SU + top];
  cuShallowRain = shallowRain;
  return fallen;
}
fn mixConserved(i: i32, pi: f32, dt: f32) {
  var upper: array<f32, K>; var lower: array<f32, K>; var level: array<f32, K>; var total: array<f32, K>;
  let n = K - KTOP;
  for (var j = 0; j < n; j++) {
    let k = KTOP + j; let idx = k * C + i;
    let mass = pi * LV[L_DS + k] / GRAV;
    upper[j] = select(0.0, dt * PH[PH_MIX + (k - 1) * C + i] / mass, j > 0);
    lower[j] = select(0.0, dt * PH[PH_MIX + k * C + i] / mass, k < K - 1);
    level[j] = IN[S_TH + idx] - LHEAT * IN[S_QC + idx] / (CP * D[D_EXM + idx]);
    total[j] = IN[S_Q + idx] + IN[S_QC + idx];
  }
  thomas(n, &upper, &lower, &level);
  thomas(n, &upper, &lower, &total);
  for (var j = 0; j < n; j++) {
    let k = KTOP + j;
    if ((j > 0 && PH[PH_MIX + (k - 1) * C + i] > 0.0) || (k < K - 1 && PH[PH_MIX + k * C + i] > 0.0)) {
      let idx = k * C + i;
      IN[S_TH + idx] = level[j]; IN[S_Q + idx] = total[j]; IN[S_QC + idx] = 0.0;
    }
  }
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
}`;

export const ADJUST_KERNELS = ['adjustMix', 'adjustPlume', 'adjustRain'];

export const PHYSICS_PASSES = ['physics', 'radiation', 'longwave', 'physicsSurface'];
export const PHYSICS_KERNELS = {
  physics: PHYSICS_KERNEL,
  radiation: radiationKernel(false),
  longwave: longwaveKernel(false),
  physicsSurface: physicsSurfaceKernel(false),
  pblDiagnose: `fn blHeight(k: i32, i: i32) -> f32 { return (D[D_GEO + k * C + i] + LV[L_GABS + k]) / GRAV; }
fn blInterface(k: i32, i: i32, zb: f32) -> f32 { return 0.5 * (blHeight(k, i) + blHeight(k + 1, i)) - zb; }
fn blDensity(k: i32, i: i32, pi: f32) -> f32 { let idx = k * C + i; return pi * LV[L_SM + k] / (RGAS * IN[S_TH + idx] * D[D_EXM + idx]); }
fn parcelVirtual(level: f32, total: f32, k: i32, i: i32, pi: f32) -> f32 {
  let ex = D[D_EXM + k * C + i]; let liquidT = level * ex; let qs = qsat(liquidT, pi * LV[L_SM + k]);
  if (!(total > qs)) { return level * (1.0 + VIRT * total); }
  let slope = qs * LHEAT / (RVAP * liquidT * liquidT); let liquid = (total - qs) / (1.0 + LHEAT * slope / CP);
  return (liquidT + LHEAT * liquid / CP) / ex * (1.0 + VIRT * (total - liquid) - liquid);
}
fn blSurfaceInterface(i: i32, zb: f32, depthAbove: f32) -> i32 {
  var kE = -1;
  for (var k = KTOP; k < K - 1; k++) { if (blInterface(k, i, zb) >= depthAbove) { kE = k; } }
  return kE;
}
fn blEntrain(i: i32, pi: f32, kE: i32, lowest: i32, h: f32, buoyant: f32, sheared: f32) -> f32 {
  if (kE < KTOP || !(h > 0.0)) { return 0.0; }
  var weight = 0.0; var sumV = 0.0; var sumL = 0.0; var sumQ = 0.0;
  for (var k = kE + 1; k <= lowest; k++) {
    let idx = k * C + i;
    weight += LV[L_DS + k]; sumV += LV[L_DS + k] * D[D_THV + idx];
    sumL += LV[L_DS + k] * (IN[S_TH + idx] - LHEAT * IN[S_QC + idx] / (CP * D[D_EXM + idx])); sumQ += LV[L_DS + k] * (IN[S_Q + idx] + IN[S_QC + idx]);
  }
  var above = kE * C + i;
  if (BL_JUMP2 && kE > KTOP && D[D_THV + above - C] > D[D_THV + above]) { above -= C; }
  let mean = sumV / weight; let jump = GRAV * (D[D_THV + above] - mean) / mean;
  var efficiency = BL_A;
  let top = (kE + 1) * C + i;
  if (IN[S_QC + top] > CT_THRESH && BL_A2 > 0.0) {
    let ex = D[D_EXM + top]; let T = IN[S_TH + top] * ex; let qs = qsat(T, pi * LV[L_SM + kE + 1]); let dqs = qs * LHEAT / (RVAP * T * T);
    let gamma = LHEAT / CP * dqs; let c = 1.0 + VIRT * IN[S_Q + top] - IN[S_QC + top] + (1.0 + VIRT) * T * dqs;
    let jumpL = IN[S_TH + above] - LHEAT * IN[S_QC + above] / (CP * D[D_EXM + above]) - sumL / weight; let jumpQ = IN[S_Q + above] + IN[S_QC + above] - sumQ / weight;
    let virtualJump = max(D[D_THV + above] - mean, BL_BMIN * mean / GRAV);
    let saturatedJump = c / (1.0 + gamma) * jumpL + (c * LHEAT / (CP * ex * (1.0 + gamma)) - IN[S_TH + top]) * jumpQ;
    let demand = dqs * ex * jumpL - jumpQ;
    var share = 1.0;
    if (demand > 0.0) { share = min(1.0, IN[S_QC + top] * (1.0 + gamma) / demand); }
    efficiency = min(BL_AMAX, BL_A * (1.0 + BL_A2 * max(0.0, share * (1.0 - saturatedJump / virtualJump))));
  }
  var velocity = min(BL_WEMAX, (efficiency * buoyant + BL_AS * sheared) / (h * max(jump, BL_BMIN)));
  if (BL_TAPER) { velocity *= clamp((DECK_CLOSED - PH[PH_MLMGATE + i]) / (DECK_CLOSED - DECK_OPEN), 0.0, 1.0) * (1.0 - PH[PH_STRAT + i]); }
  if (!(velocity > 0.0)) { return 0.0; }
  PH[PH_MIX + kE * C + i] += 0.5 * (blDensity(kE, i, pi) + blDensity(kE + 1, i, pi)) * velocity;
  return velocity;
}
fn blMoist(i: i32, pi: f32, richardsonDepth: f32, zb: f32, buoyancy: f32, friction: f32) {
  let bottom = K - 1;
  var surfaceDepth = richardsonDepth; var surfaceLevel = -1;
  if (buoyancy > 0.0 && CT_CUMULUS > 0.0) {
    let base = bottom * C + i;
    let mixed = pow(friction * friction * friction + 0.6 * buoyancy * max(0.0, surfaceDepth), 1.0 / 3.0);
    let level = IN[S_TH + base] - LHEAT * IN[S_QC + base] / (CP * D[D_EXM + base]) + CT_EXCESS * buoyancy * D[D_THV + base] / (GRAV * mixed);
    let total = IN[S_Q + base] + IN[S_QC + base];
    var k = bottom - 1; var condensation = -1;
    loop {
      if (k < KTOP) { break; }
      if (condensation < 0 && total > qsat(level * D[D_EXM + k * C + i], pi * LV[L_SM + k])) { condensation = k; }
      if (!(parcelVirtual(level, total, k, i, pi) + CT_TOLERANCE > D[D_THV + k * C + i])) { break; }
      k--;
    }
    if (condensation > k && k >= KTOP) {
      let parcelTop = blInterface(k, i, zb); let cloudBase = blInterface(condensation, i, zb);
      if (parcelTop - cloudBase <= CT_CUMULUS && parcelTop <= CT_HMAX && parcelTop > surfaceDepth) { surfaceDepth = parcelTop; surfaceLevel = k; }
    }
  }
  var top = -1; var cooling = 0.0;
  for (var k = bottom; k > KTOP; k--) {
    if (blInterface(k - 1, i, zb) > CT_HMAX) { break; }
    if (IN[S_QC + k * C + i] > CT_THRESH && !(IN[S_QC + (k - 1) * C + i] > CT_THRESH)) { top = k; break; }
  }
  var runBottom = -1;
  if (top >= 0) {
    for (var k = top; k <= bottom; k++) { if (!(IN[S_QC + k * C + i] > CT_THRESH)) { break; } cooling -= PH[PH_LWH + k * C + i]; runBottom = k; }
    if (!(cooling > 0.0)) { top = -1; cooling = 0.0; runBottom = -1; }
  }
  if (BL_CLOUDLAYER) { PH[PH_CLOUDK + i] = f32(runBottom); }
  var coupled = false; var lowest = bottom; var base0 = 0.0; var cloudTopZ = 0.0;
  if (top >= 0) {
    let idx = top * C + i;
    let level = IN[S_TH + idx] - LHEAT * IN[S_QC + idx] / (CP * D[D_EXM + idx]) - CT_PERT; let total = IN[S_Q + idx] + IN[S_QC + idx];
    var k = top + 1;
    loop {
      if (k > bottom) { break; }
      if (!(parcelVirtual(level, total, k, i, pi) < D[D_THV + k * C + i])) { break; }
      k++;
    }
    lowest = k - 1;
    if (k <= bottom) { base0 = blInterface(k - 1, i, zb); }
    coupled = k > bottom || select(base0 <= surfaceDepth, k - 1 >= surfaceLevel, surfaceLevel >= 0);
    if (coupled) { base0 = 0.0; }
    cloudTopZ = blInterface(top - 1, i, zb);
  }
  var h = select(surfaceDepth, cloudTopZ, coupled);
  PH[PH_DEPTH + i] = zb + h;
  if (MLM_PROGNOSTIC && PH[PH_MLMTOP + i] > 0.0) { h = max(h, PH[PH_MLMTOP + i] - zb); }
  PH[PH_REGIME + i] = select(select(0.0, 1.0, buoyancy > 0.0), select(2.0, 3.0, coupled), top >= 0);
  PH[PH_MIXTOP + i] = zb + max(h, cloudTopZ);
  PH[PH_CTCOOL + i] = cooling;
  let layerDepth = cloudTopZ - base0;
  var velocityCubed = 0.0;
  if (top >= 0 && layerDepth > 0.0) { velocityCubed = GRAV / D[D_THV + top * C + i] * cooling / (blDensity(top, i, pi) * CP) * layerDepth; }
  var velocity = 0.0;
  if (velocityCubed > 0.0) { velocity = pow(velocityCubed, 1.0 / 3.0); }
  PH[PH_VRAD + i] = velocity;
  var scale = friction;
  if (STABILITY && buoyancy > 0.0 && h > 0.0) { scale = friction * pow(1.0 - 15.0 * max(-2.0, -0.1 * h * KARMAN * buoyancy / (friction * friction * friction)), 0.25); }
  for (var k = KTOP; k < K - 1; k++) {
    let z = blInterface(k, i, zb);
    var diffusivity = 0.0;
    if (z < h) { diffusivity += KARMAN * scale * z * (1.0 - z / h) * (1.0 - z / h); }
    if (velocity > 0.0 && z > base0 && z < cloudTopZ) { let x = (z - base0) / layerDepth; diffusivity += CT_PROFILE * KARMAN * velocity * layerDepth * x * x * sqrt(1.0 - x); }
    if (diffusivity > 0.0) { PH[PH_MIX + k * C + i] = 0.5 * (blDensity(k, i, pi) + blDensity(k + 1, i, pi)) * diffusivity / (blHeight(k, i) - blHeight(k + 1, i)); }
  }
  if (!BL_ENTRAIN) { return; }
  var onset = 1.0;
  if (BL_ONSET > 0.0) { onset = min(1.0, buoyancy / BL_ONSET); }
  let sheared = select(0.0, onset * friction * friction * friction, buoyancy > 0.0);
  if (top >= 0) {
    let driven = coupled && buoyancy > 0.0;
    PH[PH_ENTRAIN + i] = blEntrain(i, pi, top - 1, lowest, select(layerDepth, cloudTopZ, coupled), velocityCubed + select(0.0, buoyancy * cloudTopZ, driven), select(0.0, sheared, driven));
    if (!coupled && buoyancy > 0.0 && surfaceDepth > 0.0) {
      let kE = select(blSurfaceInterface(i, zb, surfaceDepth), surfaceLevel, surfaceLevel >= 0);
      if (kE >= lowest) { _ = blEntrain(i, pi, kE, K - 1, surfaceDepth, buoyancy * surfaceDepth, sheared); }
    }
  } else if (buoyancy > 0.0 && h > 0.0) {
    PH[PH_ENTRAIN + i] = blEntrain(i, pi, select(blSurfaceInterface(i, zb, h), surfaceLevel, surfaceLevel >= 0 && h == surfaceDepth), K - 1, h, buoyancy * h, sheared);
  }
}
@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= C) { return; }
  diagnoseColumnMid(i);
  let pi = IN[S_PI + i]; let base = (K - 1) * C + i;
  var edges = edgesOf(i);
  let bottomWind = edgeWind(&edges, i, K - 1);
  let speed = length(bottomWind);
  let friction = sqrt(PH[PH_DRAG + i]) * xWind(i, speed);
  let zb = (D[D_GEO + base] + LV[L_GABS + K - 1]) / GRAV;
  var found = false; var riPrev = 0.0; var zPrev = zb; var depth = zb;
  for (var k = K - 2; k >= KTOP; k--) {
    if (found) { continue; }
    let idx = k * C + i;
    let z = (D[D_GEO + idx] + LV[L_GABS + k]) / GRAV;
    let dw = edgeWind(&edges, i, k) - bottomWind;
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
  PH[PH_ENTRAIN + i] = 0.0;
  let moisture = select(0.61 * IN[S_TH + base] * (qsat(IN[S_TS + i], pi) - IN[S_Q + base]), 0.0, PH[PH_LAND + i] > 0.5);
  if (IMPLICIT_DRAG) { PH[PH_SDRAG + i] = blDensity(K - 1, i, pi) * PH[PH_DRAG + i] * xWind(i, speed); }
  if (FORM_DRAG && PH[PH_OFLT + i] > 0.0) {
    let sflt = PH[PH_OFLT + i];
    for (var k = KTOP; k < K; k++) {
      var rate = 0.0;
      let z = (D[D_GEO + k * C + i] + LV[L_GABS + k]) / GRAV;
      if (sflt > 0.0 && z > 0.0) { rate = TOFD_SCALE * sflt * sflt * exp(-pow(z / TOFD_DECAY, 1.5)) * pow(z, -1.2) * length(edgeWind(&edges, i, k)); }
      PH[PH_TOFD + (k - KTOP) * C + i] = rate;
    }
  }
  let buoyancy = select(GRAV / IN[S_TH + base] * PH[PH_DRAG + i] * max(speed, GUST) * (IN[S_TS + i] * pow(LV[L_SM + K - 1], KAPPA) / D[D_EXM + base] - IN[S_TH + base] + moisture), PH[PH_BUOY + i], ROUGH);
  PH[PH_BUOY + i] = buoyancy; PH[PH_USTAR + i] = friction;
  if (MOIST_BL) { blMoist(i, pi, depth - zb, zb, buoyancy, friction); return; }
  if (h <= 0.0) { return; }
  var scale = friction;
  if (STABILITY && buoyancy > 0.0) { scale = friction * pow(1.0 - 15.0 * max(-2.0, -0.1 * h * KARMAN * buoyancy / (friction * friction * friction)), 0.25); }
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
  let taper = clamp((DECK_CLOSED - PH[PH_MLMGATE + i]) / (DECK_CLOSED - DECK_OPEN), 0.0, 1.0) * (1.0 - PH[PH_STRAT + i]);
  if (BL_ENTRAIN && entrainK >= KTOP && buoyancy > 0.0 && taper > 0.0) {
    let idx = entrainK * C + i; let below = idx + C;
    var weight = 0.0; var sum = 0.0;
    for (var k = entrainK + 1; k < K; k++) { weight += LV[L_DS + k]; sum += LV[L_DS + k] * D[D_THV + k * C + i]; }
    let mean = sum / weight; let jump = GRAV * (D[D_THV + idx] - mean) / mean;
    var onset = 1.0;
    if (BL_ONSET > 0.0) { onset = min(1.0, buoyancy / BL_ONSET); }
    let velocity = taper * min(BL_WEMAX, (BL_A * buoyancy + BL_AS * onset * friction * friction * friction / h) / max(jump, BL_BMIN));
    PH[PH_ENTRAIN + i] = velocity;
    PH[PH_MIX + idx] = 0.5 * (pi * LV[L_SM + entrainK] / (RGAS * IN[S_TH + idx] * D[D_EXM + idx]) + pi * LV[L_SM + entrainK + 1] / (RGAS * IN[S_TH + below] * D[D_EXM + below])) * velocity;
  }
}`,
  adjustMix: `${ADJUST_FUNCTIONS}
@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= C) { return; }
  let pi = IN[S_PI + i]; let dt = P[0];
  var mixes = false;
  for (var k = KTOP; k < K - 1; k++) { if (PH[PH_MIX + k * C + i] > 0.0) { mixes = true; } }
  if (mixes) {
    if (MOIST_BL) { mixConserved(i, pi, dt); }
    else { mixField(S_TH, i, pi, dt); mixField(S_Q, i, pi, dt); mixField(S_QC, i, pi, dt); }
  }
  diagnoseColumn(i);
  saturateColumn(i, pi);
  saturateOutOfLine(i, pi, dt);
}`,
  adjustPlume: `${ADJUST_FUNCTIONS}
@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= C) { return; }
  let pi = IN[S_PI + i]; let dt = P[0];
  let produced = plumeColumn(i, pi, dt);
  if (PH[PH_CUMF + i] > 0.0) { saturateColumn(i, pi); }
  PH[PH_CURAIN + i] = produced; PH[PH_CUDEEP + i] = select(0.0, 1.0, cuDeep); PH[PH_CUBASE + i] = f32(cuBase); PH[PH_CUSHALLOW + i] = cuShallowRain; PH[PH_CUSNOW + i] = cuSnow;
  if (cuDeep) {
    for (var k = 0; k < K; k++) {
      let idx = k * C + i;
      PH[PH_CUFALL + idx] = cuFall[k]; PH[PH_CURESERVE + idx] = cuReserve[k];
      if (PL_MIXED) { PH[PH_CUFROZEN + idx] = cuFrozen[k]; }
    }
  }
  saturateOutOfLine(i, pi, dt);
}`,
  adjustRain: `${ADJUST_FUNCTIONS}
@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= C) { return; }
  let pi = IN[S_PI + i]; let dt = P[0]; let bottom = K - 1;
  let produced = PH[PH_CURAIN + i]; let deep = PH[PH_CUDEEP + i] > 0.0; let plumeBase = i32(PH[PH_CUBASE + i]);
  // autoconversion, and the rain evaporating as it falls
  var rained = 0.0; var convective = 0.0; var descending = 0.0; var moved = false; var frozen = 0.0; var settling = 0.0; var fogged = 0.0;
  let floor = PH[PH_DEPTH + i];
  let iceConc = PH[PH_CONC + i]; let iced = select(0.0, select(iceConc, 1.0, iceConc <= 0.0), IN[S_ICE + i] > 0.0);
  for (var k = 0; k < K; k++) {
    let idx = k * C + i;
    let mass = pi * LV[L_DS + k] / GRAV;
    if (rained > 0.0 && (EVAP_IN_CLOUD || !(IN[S_QC + idx] > CLEAR_AIR)) && RAIN_EVAP > 0.0) {
      let ex = D[D_EXM + idx];
      let temperature = IN[S_TH + idx] * ex;
      let saturated = cloudSat(temperature, pi * LV[L_SM + k]);
      let qs = saturated.x;
      let slope = saturated.y;
      let deficit = max(0.0, (qs - IN[S_Q + idx]) / (1.0 + LHEAT * slope / CP)) * mass;
      let evaporated = min(rained, RAIN_EVAP * deficit);
      if (evaporated > 0.0) {
        rained = max(0.0, rained - evaporated);
        IN[S_Q + idx] += evaporated / mass;
        IN[S_TH + idx] -= LHEAT * evaporated / (mass * CP * ex);
      }
    }
    if (deep) {
      convective = max(0.0, convective + PH[PH_CUFALL + idx]);
      if (PL_MIXED) {
        let arriving = frozen + PH[PH_CUFROZEN + idx];
        frozen = min(convective, max(0.0, arriving));
        let fusion = max(0.0, arriving) - frozen + min(0.0, arriving);
        if (fusion != 0.0) { IN[S_TH + idx] -= LFUSION * fusion / (mass * CP * D[D_EXM + idx]); }
      }
      let spare = convective - PH[PH_CURESERVE + idx];
      if (spare > 0.0 && k > plumeBase && PL_EVAP > 0.0 && RAIN_EVAP > 0.0 && (EVAP_IN_CLOUD || !(IN[S_QC + idx] > CLEAR_AIR))) {
        let ex = D[D_EXM + idx];
        let temperature = IN[S_TH + idx] * ex;
        let qs = qsat(temperature, pi * LV[L_SM + k]);
        let slope = qs * LHEAT / (RVAP * temperature * temperature);
        let airborne = convective * relaxedFraction(PL_EVAP * max(0.0, 1.0 - IN[S_Q + idx] / qs) * RGAS * temperature * LV[L_DS + k] / (LV[L_SM + k] * GRAV));
        let evaporated = min(spare, min(airborne, RAIN_EVAP * max(0.0, (qs - IN[S_Q + idx]) / (1.0 + LHEAT * slope / CP)) * mass));
        if (evaporated > 0.0) {
          var sublimated = 0.0;
          if (frozen > 0.0) { sublimated = evaporated * frozen / convective; }
          convective -= evaporated;
          IN[S_Q + idx] += evaporated / mass;
          if (PL_MIXED) { frozen = max(0.0, frozen - sublimated); IN[S_TH + idx] -= (LHEAT * evaporated + LFUSION * sublimated) / (mass * CP * ex); }
          else { IN[S_TH + idx] -= LHEAT * evaporated / (mass * CP * ex); }
        }
      }
    }
    var liquid = IN[S_QC + idx];
    if (ICE_FALL) {
      if (descending > 0.0) {
        let melted = descending * clamp((IN[S_TH + idx] * D[D_EXM + idx] - ICE_T) / (LIQUID_T - ICE_T), 0.0, 1.0);
        IN[S_QC + idx] += (descending - melted) / mass;
        rained += melted;
        descending = 0.0;
      }
      liquid = IN[S_QC + idx];
      if (liquid > 0.0) {
        let temperature = IN[S_TH + idx] * D[D_EXM + idx]; let pressure = pi * LV[L_SM + k];
        let share = clamp((temperature - ICE_T) / (LIQUID_T - ICE_T), 0.0, 1.0);
        let iceWater = (1.0 - share) * liquid;
        let total = liquid;
        liquid = share * total;
        if (iceWater > 0.0) {
          let cover = max(CLEAR_AIR, uniformCover(total, uniformWidth(cloudSat(temperature, pressure), pressure, pi)));
          let speed = FALL_C * pow(pressure / (RGAS * temperature) * iceWater / cover, FALL_EXP);
          let courant = speed * dt * LV[L_SM + k] * GRAV / (RGAS * temperature * LV[L_DS + k]);
          let leaving = iceWater * courant / (1.0 + courant);
          IN[S_QC + idx] = total - leaving;
          descending = leaving * mass;
          moved = true;
        }
      }
    }
    if (FOG && settling > 0.0) {
      IN[S_QC + idx] += settling / mass;
      liquid += settling / mass;
      settling = 0.0;
    }
    let qc = IN[S_QC + idx];
    if (!(qc > 0.0)) { continue; }
    var floored = false;
    if (AUTO_NONE) { } else if (AUTO_BL) { floored = k > 0 && upperInterface(i, k) < floor; } else { floored = k >= K - 2; }
    if (floored) {
      if ((FOG || (FOG_DEP && k == K - 1)) && liquid > 0.0) {
        let temperature = IN[S_TH + idx] * D[D_EXM + idx]; let pressure = pi * LV[L_SM + k];
        let continental = PH[PH_LAND + i] > 0.5 && !(PH[PH_LAND + i] > 1.5);
        var courant = 0.0;
        if (FOG) {
          let speed = select(FOG_SETTLE_SEA, FOG_SETTLE_LAND, continental) * pow(pressure / (RGAS * temperature) * liquid, FOG_THIRDS);
          courant = speed * dt * LV[L_SM + k] * GRAV / (RGAS * temperature * LV[L_DS + k]);
        }
        if (FOG_DEP && k == K - 1) { courant += min(FOG_E * PH[PH_SDRAG + i], FOG_VMAX * pressure / (RGAS * temperature)) * dt * GRAV / (pi * LV[L_DS + k]); }
        let leaving = liquid * courant / (1.0 + courant);
        var water = qc - leaving;
        if (k == K - 1) { rained += leaving * mass; fogged += leaving * mass; }
        else { settling = leaving * mass; moved = true; }
        if (FOG) {
          let kept = liquid - leaving;
          if (kept > 0.0) {
            let converted = kept * relaxedFraction(select(FOG_KK_SEA, FOG_KK_LAND, continental) * pow(kept, 1.47) * dt);
            water -= converted;
            rained += mass * converted;
          }
        }
        IN[S_QC + idx] = water;
      }
      continue;
    }
    let excess = max(0.0, liquid - AUTO_T);
    var lifetime = select(CLOUD_LIFE, UPPER_LIFE, UPPER_SPLIT && pi * LV[L_SM + k] < SHALLOW_TOP);
    if (STRAT_SPLIT) { lifetime += longCloudShare(i, k, pi, iced) * (STRAT_LIFE - lifetime); }
    let converted = min(liquid, excess * (1.0 - exp(-AUTO_R * dt)) + liquid * (1.0 - exp(-dt / lifetime)));
    IN[S_QC + idx] = qc - converted;
    rained += mass * converted;
  }
  rained += descending;
  if (moved) { saturateColumn(i, pi); }
  let convected = select(produced, convective + PH[PH_CUSHALLOW + i], deep);
  var snowed = 0.0;
  if (PL_MIXED) { snowed = select(PH[PH_CUSNOW + i], min(frozen, convective) + PH[PH_CUSNOW + i], deep); }
  // filler
  for (var f = 0; f < 2; f++) {
    let off = select(S_Q, S_QC, f == 1);
    for (var k = 0; k < K - 1; k++) {
      let idx = k * C + i;
      if (IN[off + idx] < 0.0) { IN[off + idx + C] += IN[off + idx] * LV[L_DS + k] / LV[L_DS + k + 1]; IN[off + idx] = 0.0; }
    }
    if (IN[off + bottom * C + i] < 0.0) { IN[off + bottom * C + i] = 0.0; }
  }
  if (PL_BECHTOLD) {
    for (var k = K - SUB_K; k < K; k++) {
      let idx = k * C + i;
      PH[PH_SUBTV + (k - K + SUB_K) * C + i] = IN[S_TH + idx] * D[D_EXM + idx] * (1.0 + PARCEL_VIRT * max(0.0, IN[S_Q + idx]) - CU_LOADING * max(0.0, IN[S_QC + idx]));
    }
  }
  PH[PH_RAIN + i] += rained + convected; PH[PH_COND + i] += rained; PH[PH_CONV + i] += convected; PH[PH_STEPRAIN + i] = rained + convected; PH[PH_FOG + i] += fogged;
  let airT = IN[S_TH + bottom * C + i] * D[D_EXM + bottom * C + i];
  if (PH[PH_LAND + i] < 0.5 && airT < MELTING && rained + convected > 0.0) {
    ${snowOnSea('rained + convected')}
    IN[S_TH + bottom * C + i] += LFUS * (rained + convected - snowed) * GRAV / (CP * pi * LV[L_DS + K - 1] * D[D_EXM + bottom * C + i]);
  }
  if (PL_MIXED && snowed > 0.0 && !(airT < MELTING)) { IN[S_TH + bottom * C + i] -= LFUS * snowed * GRAV / (CP * pi * LV[L_DS + K - 1] * D[D_EXM + bottom * C + i]); }
  if (PH[PH_LAND + i] > 0.5) {
    if (VEGETATED) { PH[PH_RAINMEAN + i] += ((rained + convected) * 86400.0 / dt - PH[PH_RAINMEAN + i]) * P[6]; }
    if (airT < MELTING) {
      PH[PH_SNOW + i] += rained + convected;
      PH[PH_SNOWALB + i] = refreshedSnow(PH[PH_SNOWALB + i], rained + convected);
      IN[S_TH + bottom * C + i] += LFUS * (rained + convected - snowed) * GRAV / (CP * pi * LV[L_DS + K - 1] * D[D_EXM + bottom * C + i]);
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
  saturateOutOfLine(i, pi, dt);
}`,
  mixMomentum: `@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let e = i32(id.x); if (e >= E) { return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  var mixes = false;
  for (var k = KTOP; k < K - 1; k++) { if (PH[PH_MIX + k * C + a] + PH[PH_MIX + k * C + b] > 0.0) { mixes = true; } }
  if (FORM_DRAG) {
    PH[PH_FSTRESS + e] = 0.0;
    if (PH[PH_OFLT + a] + PH[PH_OFLT + b] > 0.0) { mixes = true; }
  }
  if (mixes || IMPLICIT_DRAG) { mixEdge(e, a, b); }
  if (PL_MOMENTUM) { transportEdge(e, a, b); }
}
fn transportEdge(e: i32, a: i32, b: i32) {
  let bottom = K - 1; let dt = P[0];
  var up: array<f32, K + 1>; var down: array<f32, K + 1>; var flux: array<f32, K + 1>; var before: array<f32, K>;
  var moving = false;
  for (var k = 1; k < K; k++) {
    up[k] = 0.5 * (PH[PH_MOMU + k * C + a] + PH[PH_MOMU + k * C + b]);
    down[k] = 0.5 * (PH[PH_MOMD + k * C + a] + PH[PH_MOMD + k * C + b]);
    if (up[k] > 0.0 || down[k] < 0.0) { moving = true; }
  }
  if (!moving) { return; }
  for (var k = 0; k < K; k++) { before[k] = IN[S_U + k * E + e]; }
  let source = min(i32(PH[PH_MOMS + a]), i32(PH[PH_MOMS + b])); let columnMass = 0.5 * (IN[S_PI + a] + IN[S_PI + b]);
  var rising = 0.0; var below = 0.0; var weight = 0.0;
  for (var j = bottom; j >= 1; j--) {
    if (j >= source) { below += LV[L_DS + j] * before[j]; weight += LV[L_DS + j]; rising = below / weight; }
    else { rising = before[j] + (rising - before[j]) * 0.5 * (PH[PH_MOMK + j * C + a] + PH[PH_MOMK + j * C + b]); }
    flux[j] = up[j] * (rising - before[j - 1]);
  }
  var sinking = 0.0;
  for (var j = 1; j < K; j++) {
    sinking = before[j - 1] + (sinking - before[j - 1]) * 0.5 * (PH[PH_MOMKD + (j - 1) * C + a] + PH[PH_MOMKD + (j - 1) * C + b]);
    flux[j] += down[j] * (sinking - before[j]);
  }
  var share: array<f32, K>;
  var loss = 0.0; var total = 0.0;
  for (var k = 0; k < K; k++) {
    let mass = columnMass * LV[L_DS + k] / GRAV; let now = before[k] + (flux[k + 1] - flux[k]) * dt / mass;
    IN[S_U + k * E + e] = now;
    loss += mass * (before[k] * before[k] - now * now);
    share[k] = mass * (now - before[k]) * (now - before[k]);
    total += share[k];
  }
  if (!(loss > 0.0) || !(total > 0.0)) { return; }
  for (var k = 0; k < K; k++) { D[D_DISS + k * E + e] += loss * share[k] / (total * columnMass * LV[L_DS + k] / GRAV); }
}
fn mixOwn(j: i32, a: i32, b: i32) -> f32 { return 0.5 * (PH[PH_TOFD + j * C + a] + PH[PH_TOFD + j * C + b]); }
/*
 * The share of the dissipated kinetic energy layer j takes, its own change
 * and the shear against its neighbours, from the mixed winds in rhs and
 * the winds still in the state; summed over the column it is the whole.
 */
fn mixShare(j: i32, n: i32, e: i32, a: i32, b: i32, columnMass: f32, dt: f32, surfaceDrag: f32, forming: bool, rhs: ptr<function, array<f32, K>>) -> f32 {
  let k = KTOP + j; let mass = columnMass * LV[L_DS + k] / GRAV;
  let change = (*rhs)[j] - IN[S_U + k * E + e];
  var share = mass * change * change;
  if (j > 0) { let shear = (*rhs)[j - 1] - (*rhs)[j]; share += dt * 0.5 * (PH[PH_MIX + (k - 1) * C + a] + PH[PH_MIX + (k - 1) * C + b]) * shear * shear; }
  if (j < n - 1) { let shear = (*rhs)[j] - (*rhs)[j + 1]; share += dt * 0.5 * (PH[PH_MIX + k * C + a] + PH[PH_MIX + k * C + b]) * shear * shear; }
  if (IMPLICIT_DRAG && j == n - 1) { share += dt * surfaceDrag * (*rhs)[j] * (*rhs)[j]; }
  if (forming) { share += dt * columnMass * LV[L_DS + k] / GRAV * mixOwn(j, a, b) * (*rhs)[j] * (*rhs)[j]; }
  return share;
}
/*
 * The implicit mixing of the edge's wind down the column: the tridiagonal
 * solve keeps only the winds and the back-sweep gains per layer, the
 * coefficients read again from the mixing rates where they are used, and
 * the state's winds stand until the dissipation has been shared out.
 */
fn mixEdge(e: i32, a: i32, b: i32) {
  let columnMass = 0.5 * (IN[S_PI + a] + IN[S_PI + b]); let dt = P[0];
  let n = K - KTOP;
  let surfaceDrag = 0.5 * (PH[PH_SDRAG + a] + PH[PH_SDRAG + b]);
  var forming = false;
  if (FORM_DRAG) { for (var j = 0; j < n; j++) { if (mixOwn(j, a, b) > 0.0) { forming = true; } } }
  var rhs: array<f32, K>; var gain: array<f32, K>;
  var mixAbove = 0.0;
  for (var j = 0; j < n; j++) {
    let k = KTOP + j;
    let mass = columnMass * LV[L_DS + k] / GRAV;
    let mixHere = PH[PH_MIX + k * C + a] + PH[PH_MIX + k * C + b];
    let upper = select(0.0, dt * 0.5 * mixAbove / mass, j > 0);
    let lower = select(select(0.0, dt * 0.5 * (PH[PH_SDRAG + a] + PH[PH_SDRAG + b]) / mass, IMPLICIT_DRAG), dt * 0.5 * mixHere / mass, k < K - 1);
    let u = IN[S_U + k * E + e];
    var denominator = 0.0;
    if (j == 0) {
      if (forming) { denominator = 1.0 + upper + lower + dt * mixOwn(j, a, b); } else { denominator = 1.0 + upper + lower; }
      rhs[j] = u / denominator;
    } else {
      if (forming) { denominator = 1.0 + upper + lower + dt * mixOwn(j, a, b) + upper * gain[j - 1]; } else { denominator = 1.0 + upper + lower + upper * gain[j - 1]; }
      rhs[j] = (u + upper * rhs[j - 1]) / denominator;
    }
    gain[j] = -lower / denominator;
    mixAbove = mixHere;
  }
  for (var j = n - 2; j >= 0; j--) { rhs[j] -= gain[j] * rhs[j + 1]; }
  if (IMPLICIT_DRAG) { PH[PH_STRESS + e] = surfaceDrag * rhs[n - 1]; PH[PH_STRESSOK] = 1.0; }
  if (forming) {
    var formed = 0.0;
    for (var j = 0; j < n; j++) { formed += columnMass * LV[L_DS + KTOP + j] / GRAV * mixOwn(j, a, b) * rhs[j]; }
    PH[PH_FSTRESS + e] = formed;
  }
  var loss = 0.0; var total = 0.0;
  for (var j = 0; j < n; j++) {
    let k = KTOP + j; let mass = columnMass * LV[L_DS + k] / GRAV; let before = IN[S_U + k * E + e];
    loss += mass * (before * before - rhs[j] * rhs[j]);
    total += mixShare(j, n, e, a, b, columnMass, dt, surfaceDrag, forming, &rhs);
  }
  if (total > 0.0) {
    for (var j = 0; j < n; j++) {
      let k = KTOP + j;
      D[D_DISS + k * E + e] += loss * mixShare(j, n, e, a, b, columnMass, dt, surfaceDrag, forming, &rhs) / (total * columnMass * LV[L_DS + k] / GRAV);
    }
  }
  for (var j = 0; j < n; j++) { IN[S_U + (KTOP + j) * E + e] = rhs[j]; }
}`,
};
