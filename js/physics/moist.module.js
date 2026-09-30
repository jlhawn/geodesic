import { R_DRY, VIRTUAL_FACTOR } from '../dynamics/sigmaCore.module.js';

export const LATENT_HEAT = 2.5e6;
export const EPSILON = 0.622;
export const R_VAPOR = R_DRY / EPSILON;
export const ACTIVITY_UNDECIDED = 0.5;
export const CLEAR_AIR = 1e-7;
const STABILITY_SIGMA = 0.7;
export const DECK_CLOSED = 0.6;
export const CUMULUS_FLOOR = 1e-6;

export function saturationVaporPressure(T) {
  return 611.2 * Math.exp(17.67 * (T - 273.15) / (T - 29.65));
}

export function saturationHumidity(T, p) {
  const es = saturationVaporPressure(T);
  const dry = p - (1 - EPSILON) * es;
  return dry > 0 ? EPSILON * es / dry : 1;
}

/*
 * Lifting condensation level of a parcel (T, q, p) by Bolton (1980):
 * the dew point from the vapour pressure, the LCL temperature from his
 * eq. 15, and the pressure along the dry adiabat. Returns null when the
 * parcel is already saturated or has no vapour.
 */
export function liftingCondensationLevel(T, q, p, kappa) {
  if (q <= 0) return null;
  const e = q * p / (EPSILON + (1 - EPSILON) * q);
  const y = Math.log(e / 611.2);
  const dewPoint = (273.15 * 17.67 - 29.65 * y) / (17.67 - y);
  if (dewPoint >= T) return { temperature: T, pressure: p };
  const temperature = 1 / (1 / (dewPoint - 56) + Math.log(T / dewPoint) / 800) + 56;
  return { temperature, pressure: p * Math.pow(temperature / T, 1 / kappa) };
}

/*
 * Moist physics for the sigma core, applied to the state after each
 * step in the order: saturation adjustment, convection, rain, filler.
 *
 * Saturation adjustment: supersaturated vapour condenses into cloud
 * water and cloud water evaporates into subsaturated air, latent heat to
 * the layer.
 *
 * Convection is the Betts–Miller relaxation of Frierson (2007) behind a
 * trigger. Its parcel is the lowest layer's air or, with `boundaryParcel`,
 * the mass-weighted mean θ and q of the boundary layer, the layers whose
 * lower interfaces lie below `boundaryDepth` (per cell, the boundary
 * layer's Richardson depth in the height of the core's geopotential),
 * or of the layers whose midpoints lie within `parcelDepth` of the
 * surface where that reaches higher. It rises dry to its lifting
 * condensation level, then saturated, its moist static energy diluted
 * toward the air's at the fractional rate `entrainmentRate` per metre,
 * and is buoyant where its virtual temperature exceeds the air's. Its
 * inhibition is the negative buoyant energy from the top of its source
 * layers to its level of free convection, the first buoyant layer above
 * the condensation level; its CAPE the positive energy above that; the
 * top the highest buoyant layer, stable layers in between ignored, the
 * search ending where the parcel is 10 K colder than the air. Tops below
 * `shallowTop` are shallow, the others deep.
 *
 * A column's pass is the product of two ramps from 0 to 1, one as the
 * CAPE rises from half `capeThreshold` to one and a half times it, the
 * other as the inhibition falls from one and a half `inhibitionThreshold`
 * to half of it, and of the deck's opening, 1 where the mixed-layer
 * deck's gate `deckGate` (radiation.mlmGate) is at most one half,
 * falling linearly to 0 at DECK_CLOSED. The per-cell `activity` relaxes
 * toward the pass over `activityMemory`, and the column fires while it
 * is above one half, or at one half on a pass above one half, so a
 * column fires while its CAPE has stood above the threshold and its
 * inhibition below it for a while; `activity` is saved with the state
 * (key `convectiveActivity`, one half in a state saved without it).
 * With `shallowScheme` 'massFlux' only deep tops fire; the shallow
 * cumulus mass flux below takes the shallow ones, and runs after the deep
 * branch in columns where that fires too (without `cumulusWithDeep`, not
 * there). With 'bettsMiller' a shallow top fires on
 * the deep trigger too and also vents without memory: at the vent, the same product
 * of ramps with `shallowCape` and `shallowInhibition` in place of the
 * deep thresholds, times the deck's opening, and 0 where the column's
 * estimated inversion strength (the radiation's EIS) exceeds
 * `shallowStability` (null: no bound; `shallowCape` null: no venting).
 * A firing column relaxes over relaxationTime, a venting one over
 * relaxationTime divided by its vent. Only the layers from cloud base,
 * the layer holding the condensation level, to the top relax; the
 * subcloud layers are the boundary layer's.
 *   Toward the parcel (deep tops, and shallow ones with `shallowReference`
 *   'parcel'): the parcel's temperature and `referenceHumidity` of its
 *   saturation, the temperature shifted so the enthalpy change equals
 *   the latent heat of the water removed, which rains, or, where the
 *   layers would have to moisten, both shifted so that neither heat nor
 *   water changes and nothing rains.
 *   Toward the mixing line (shallow tops with 'mixingLine'): the mixing
 *   line of Betts (1986) between the parcel at its condensation level
 *   and the air of the layer above the top, moist static energy and
 *   water mixed linearly in pressure, the water at most
 *   `shallowHumidity` of saturation (the rest of the mixture's energy in
 *   its temperature).
 *   Shallow tops rain as deep ones do with `shallowRain`; without it they
 *   are shifted so that neither heat nor water changes.
 *   Of any rain the fraction `detrainment` stays as cloud water spread
 *   through the anvil, the top `anvilDepth` of pressure below the top,
 *   and `downdraftEvaporation` of the rest may evaporate into the
 *   subcloud layers as it falls, the proxy of a downdraft, offered to
 *   each in proportion to its mass (`downdraftSpread` 'mass') or all of
 *   it to each in turn from cloud base down ('fall').
 *
 * The shallow cumulus mass flux, a bulk entraining plume after
 * Bretherton, McCaa and Grenier (2004), simplified. Its source is the
 * mass-weighted mean liquid-water static energy s_l = c_p T + g z − L q_c
 * and total water q_t of the layers whose lower interfaces lie below
 * the boundary-layer top, or whose midpoints lie within
 * `cumulusSourceDepth` of the surface, and below `shallowTop`. The plume
 * leaves the source's top interface with its air, rises without mixing
 * to the Bolton LCL of that air and above it entrains the layer's s_l
 * and q_t at `cumulusEntrainment` and detrains at `cumulusDetrainment`
 * (per metre), its mass flux growing by exp((ε − δ) Δz) through each
 * layer. At each layer's midpoint the plume's temperature and condensate
 * follow from its s_l and q_t, and its buoyancy is its virtual
 * temperature against the air's, condensate loading included in both.
 * Its inhibition is the negative buoyant energy from the source's top to
 * the first layer where it is saturated; from there on the plume ends in
 * the first layer where it is not buoyant, or in the last layer below
 * `shallowTop`, and detrains all it carries there, `cumulusOvershoot` of
 * the mass flux entering that layer and the rest in the layer below. The
 * base mass flux is ρ_LCL c w exp(−inhibition / w²), c the
 * `cumulusClosure`, w the larger of (B₀ h)^⅓ and `cumulusFriction` u*,
 * B₀ the boundary layer's surface buoyancy flux and h its depth above the
 * lowest layer, times the deck's opening, and zero where B₀ ≤ 0, below
 * CUMULUS_FLOOR (kg/m²/s), where
 * the LCL lies above `shallowTop` or where the plume does not saturate
 * below it; at most `cumulusBoundaryLoss` of the source
 * layers' mass leaves in a step. Inside the source the flux draws on
 * each layer in proportion to its mass. The layers change in flux form:
 * the flux M (X_u − X_above) of s_l and q_t through each interface, X_u
 * the plume's and X_above the layer above's, and each layer's X by
 * g Δt/Δp times the flux through its lower interface less that through
 * its upper one, the mass flux scaled down where needed so that
 * M g Δt/Δp ≤ 1 in every layer; the change in s_l goes to the
 * temperature, that in q_t to the vapour, and a column the plume moved
 * is saturation-adjusted again, so detrained condensate evaporates into
 * dry layers and stays as cloud water in saturated ones. Column enthalpy
 * and water are exact. With `cumulusRain` (kg/kg; null: none) the
 * plume's condensate above it rains at each interface, convective rain.
 * `cumulusCover` holds each plume layer's cumulus fraction, its mean
 * mass flux over ρ `cumulusUpdraft`, and `cumulusWater` the plume's
 * condensate at its midpoint, for the radiation; `cumulusBaseFlux` the
 * base mass flux (kg/m²/s) and `cumulusTop` the pressure of the plume's
 * top interface.
 *
 * Rain: Kessler autoconversion of cloud water above the threshold at
 * autoconversionRate, and of all cloud water over cloudLifetime, except
 * in the lowest two layers (`autoconversionFloor` 'lowest') or in the
 * layers wholly below the boundary-layer top ('boundaryLayer'; the
 * lowest two without a boundary layer). The rain falls through the layers below within the
 * step and evaporates into each cloud-free (at most CLEAR_AIR of cloud
 * water) subsaturated one up to `rainEvaporation` of what would saturate
 * it, latent cooling included,
 * as does the downdraft's share of the convective rain in the subcloud
 * layers. A filler removes negative humidity by borrowing from the layer
 * below.
 *
 * Precipitation accumulates per cell (kg/m²), and so do its two parts:
 * convective, the Betts–Miller rain less its detrained share and less
 * what its downdraft evaporates; large-scale, the autoconversion rain
 * less what evaporates on the way down, which includes the detrained
 * anvil water once it rains out. `readRain` turns the parts' sums into
 * their means over the interval they cover, in mm/d. The budget sums are
 * area-weighted masses (kg). `trace`, when its arrays (K·C) are set,
 * receives each layer's temperature change (K) from convection, its
 * downdraft included, and from the large-scale condensation,
 * autoconversion and rain evaporation.
 *
 * Defaults: relaxationTime 2 h, referenceHumidity 0.6, the lowest
 * layer's parcel (boundaryParcel false, parcelDepth 0), entrainmentRate
 * 5e-5 /m, capeThreshold 100 J/kg, inhibitionThreshold 50 J/kg,
 * activityMemory 2 h, shallowTop 700 hPa, shallowCape 10 J/kg,
 * shallowInhibition 15 J/kg, no shallowStability, shallowReference
 * 'parcel', shallowRain, shallowHumidity 0.8, detrainment 0.1,
 * anvilDepth 150 hPa, downdraftEvaporation 0.01 spread by mass,
 * autoconversionThreshold 2e-4, autoconversionRate 1e-3 /s,
 * cloudLifetime 3 h, autoconversionFloor 'lowest', rainEvaporation 1,
 * shallowScheme 'massFlux', cumulusClosure 0.04, cumulusEntrainment
 * 1.5e-3 /m, cumulusDetrainment 2e-3 /m, cumulusSourceDepth 50 hPa,
 * cumulusBoundaryLoss 0.1, cumulusFriction 1, cumulusOvershoot 1,
 * cumulusUpdraft 1 m/s, no cumulusRain, cumulusSource 'mean' (or
 * 'lowest': the plume leaves with the lowest layer's air), cumulusWithDeep.
 */
export const MOIST_DEFAULTS = {
  latentHeat: LATENT_HEAT, relaxationTime: 7200, referenceHumidity: 0.6, parcelDepth: 0, entrainmentRate: 5e-5,
  capeThreshold: 100, inhibitionThreshold: 50, activityMemory: 2 * 3600, shallowTop: 700e2, detrainment: 0.1, anvilDepth: 150e2,
  downdraftEvaporation: 0.01, autoconversionThreshold: 2e-4, autoconversionRate: 1e-3, cloudLifetime: 3 * 3600, rainEvaporation: 1, autoconversionFloor: 'lowest', shallowHumidity: 0.8,
  shallowCape: 10, shallowInhibition: 15, shallowStability: null, shallowReference: 'parcel', shallowRain: true,
  boundaryParcel: false, adjustFrom: 'cloudBase', deckVeto: true, evaporationInCloud: false, downdraftSpread: 'mass', virtualBuoyancy: true,
  shallowScheme: 'massFlux', cumulusClosure: 0.04, cumulusEntrainment: 1.5e-3, cumulusDetrainment: 2e-3, cumulusSourceDepth: 50e2, cumulusBoundaryLoss: 0.1,
  cumulusFriction: 1, cumulusOvershoot: 1, cumulusUpdraft: 1, cumulusRain: null, cumulusSource: 'mean', cumulusWithDeep: true,
};

export function createMoistPhysics(mesh, core, { boundaryDepth = null, deckGate = null, surfaceBuoyancy = null, frictionVelocity = null, buffers = null, ...options } = {}) {
  const {
    latentHeat, relaxationTime, referenceHumidity, parcelDepth, entrainmentRate, capeThreshold, inhibitionThreshold, activityMemory, shallowTop,
    detrainment, anvilDepth, downdraftEvaporation, autoconversionThreshold, autoconversionRate, cloudLifetime, rainEvaporation, autoconversionFloor, shallowHumidity,
    shallowCape, shallowInhibition, shallowStability, shallowReference, shallowRain, boundaryParcel, adjustFrom, deckVeto, evaporationInCloud, downdraftSpread, virtualBuoyancy,
    shallowScheme, cumulusClosure, cumulusEntrainment, cumulusDetrainment, cumulusSourceDepth, cumulusBoundaryLoss, cumulusFriction, cumulusOvershoot, cumulusUpdraft, cumulusRain, cumulusSource, cumulusWithDeep,
  } = { ...MOIST_DEFAULTS, ...options };
  if (cumulusSource !== 'mean' && cumulusSource !== 'lowest') throw new Error(`cumulusSource must be 'mean' or 'lowest', not ${cumulusSource}`);
  if (shallowScheme !== 'massFlux' && shallowScheme !== 'bettsMiller') throw new Error(`shallowScheme must be 'massFlux' or 'bettsMiller', not ${shallowScheme}`);
  const massFlux = shallowScheme === 'massFlux';
  if (shallowReference !== 'mixingLine' && shallowReference !== 'parcel') throw new Error(`shallowReference must be 'mixingLine' or 'parcel', not ${shallowReference}`);
  if (autoconversionFloor !== 'lowest' && autoconversionFloor !== 'boundaryLayer' && autoconversionFloor !== 'none') throw new Error(`autoconversionFloor must be 'lowest' or 'boundaryLayer', not ${autoconversionFloor}`);
  const { K, C, levels, dSigma, sigmaMid, cp, R, g, kappa, exnerLayer, exnerLower, geopotential } = core.diagnostics;
  const thetaV = core.arrays.thetaV;
  let stabilityLayer = 0;
  for (let k = 1; k < K; k++) if (Math.abs(sigmaMid[k] - STABILITY_SIGMA) < Math.abs(sigmaMid[stabilityLayer] - STABILITY_SIGMA)) stabilityLayer = k;
  const upperInterface = (i, k) => (geopotential[k * C + i] + cp * thetaV[k * C + i] * (exnerLayer[k * C + i] - exnerLower[(k - 1) * C + i])) / g;
  const shared = (name, n) => (buffers && buffers[name] ? buffers[name] : new SharedArrayBuffer(8 * n));
  const precipBuffer = shared('precipitation', C), rainBuffer = shared('rain', C), convectiveBuffer = shared('convectivePrecipitation', C);
  const largeScaleBuffer = shared('largeScalePrecipitation', C), activityBuffer = shared('convectiveActivity', C);
  const cumulusCoverBuffer = shared('cumulusCover', K * C), cumulusWaterBuffer = shared('cumulusWater', K * C), baseFluxBuffer = shared('cumulusBaseFlux', C), cumulusTopBuffer = shared('cumulusTop', C);
  const cumulusCover = new Float64Array(cumulusCoverBuffer), cumulusWater = new Float64Array(cumulusWaterBuffer), cumulusBaseFlux = new Float64Array(baseFluxBuffer), cumulusTop = new Float64Array(cumulusTopBuffer);
  const precipitation = new Float64Array(precipBuffer), rain = new Float64Array(rainBuffer);
  const convectivePrecipitation = new Float64Array(convectiveBuffer), largeScalePrecipitation = new Float64Array(largeScaleBuffer);
  const activity = new Float64Array(activityBuffer);
  if (!(buffers && buffers.convectiveActivity)) activity.fill(ACTIVITY_UNDECIDED);
  const convectiveRain = new Float64Array(C), largeScaleRain = new Float64Array(C);
  const T = new Float64Array(K), p = new Float64Array(K), dp = new Float64Array(K), z = new Float64Array(K), Tref = new Float64Array(K), qref = new Float64Array(K);
  const downdraftCooling = new Float64Array(K);
  const plumeS = new Float64Array(K), plumeQ = new Float64Array(K), plumeLiquid = new Float64Array(K), plumeGrowth = new Float64Array(K), plumeRain = new Float64Array(K);
  const envS = new Float64Array(K), envQ = new Float64Array(K), cumulusFlux = new Float64Array(K + 1), fluxS = new Float64Array(K + 1), fluxQ = new Float64Array(K + 1);
  const plume = { T: 0, liquid: 0 };
  const cumulus = { top: -1, source: K - 1, inhibition: 0, lclPressure: 0, velocity: 0, baseFlux: 0 };
  const parcel = { top: -1, base: -1, source: K - 1, cape: 0, inhibition: 0, lclPressure: 0, theta: 0, q: 0, energy: 0, fired: false };
  const falling = { base: K, downdraft: 0, evaporated: 0 };
  const budget = { condensation: 0, convection: 0, lost: 0 };
  const trace = { convection: null, largeScale: null };
  const marked = new Float64Array(K);
  const ramp = (x) => Math.min(1, Math.max(0, x));
  function mark(i, theta) { for (let k = 0; k < K; k++) marked[k] = theta[k * C + i]; }
  function charge(into, i, theta) {
    for (let k = 0; k < K; k++) { const idx = k * C + i; if (into) into[idx] += (theta[idx] - marked[k]) * exnerLayer[idx]; marked[k] = theta[idx]; }
  }

  /*
   * The temperature at which air at `pressure` and `humidity` of its
   * saturation holds `energy` as cp T + L h q_s(T, p), by four Newton
   * steps from `guess`.
   */
  function saturatedTemperature(energy, pressure, guess, humidity = 1) {
    let t = guess;
    for (let n = 0; n < 4; n++) {
      const qs = humidity * saturationHumidity(t, pressure);
      t -= (cp * t + latentHeat * qs - energy) / (cp + latentHeat * latentHeat * qs / (R_VAPOR * t * t));
    }
    return t;
  }

  /*
   * Saturation adjustment: supersaturated vapour condenses into cloud
   * water and cloud water evaporates into subsaturated air, each with
   * one implicit step, so afterwards a layer is either saturated or
   * cloud-free. Returns the condensate formed (kg/m², negative when
   * cloud evaporated); nothing rains here.
   */
  function condenseColumn(i, pi, theta, q, qc) {
    let formed = 0;
    for (let k = 0; k < K; k++) {
      const idx = k * C + i;
      const ex = exnerLayer[idx];
      const temperature = theta[idx] * ex;
      const pressure = pi[i] * sigmaMid[k];
      const qs = saturationHumidity(temperature, pressure);
      const slope = qs * latentHeat / (R_VAPOR * temperature * temperature);
      let change = (q[idx] - qs) / (1 + latentHeat * slope / cp);
      if (change < 0) change = Math.max(change, -qc[idx]);
      if (change === 0) continue;
      q[idx] -= change;
      qc[idx] += change;
      theta[idx] += latentHeat * change / (cp * ex);
      formed += pi[i] * dSigma[k] / g * change;
    }
    return formed;
  }

  /*
   * Autoconversion and the rain's fall (see the header). `downdraft` is
   * the convective rain (kg/m²) that may evaporate below cloud base
   * `base`; falling.evaporated receives what did. Returns the
   * autoconversion rain that reaches the ground (kg/m²).
   */
  function autoconvertColumn(i, pi, theta, q, qc, dt, downdraft = 0, base = K) {
    let rain = 0, left = downdraft;
    const floor = autoconversionFloor === 'boundaryLayer' && boundaryDepth ? boundaryDepth[i] : null;
    if (trace.convection) downdraftCooling.fill(0);
    let subcloud = 0;
    if (downdraftSpread === 'mass') for (let k = base + 1; k < K; k++) subcloud += dSigma[k];
    for (let k = 0; k < K; k++) {
      const idx = k * C + i;
      const offer = downdraftSpread === 'mass' && k > base ? Math.min(left, downdraft * dSigma[k] / subcloud) : left;
      const open = rain > 0 && (evaporationInCloud || !(qc[idx] > CLEAR_AIR)), draft = offer > 0 && k > base;
      if ((open || draft) && rainEvaporation > 0) {
        const ex = exnerLayer[idx], mass = pi[i] * dSigma[k] / g;
        const temperature = theta[idx] * ex;
        const qs = saturationHumidity(temperature, pi[i] * sigmaMid[k]);
        const slope = qs * latentHeat / (R_VAPOR * temperature * temperature);
        const deficit = Math.max(0, (qs - q[idx]) / (1 + latentHeat * slope / cp)) * mass;
        const available = (open ? rain : 0) + (draft ? offer : 0);
        const evaporated = Math.min(available, rainEvaporation * deficit);
        if (evaporated > 0) {
          let fromDraft = 0, fromRain = 0;
          if (evaporated >= available) { fromDraft = draft ? offer : 0; fromRain = open ? rain : 0; } else { fromDraft = draft ? evaporated * offer / available : 0; fromRain = evaporated - fromDraft; }
          left -= fromDraft;
          rain = Math.max(0, rain - fromRain);
          q[idx] += (fromRain + fromDraft) / mass;
          theta[idx] -= latentHeat * (fromRain + fromDraft) / (mass * cp * ex);
          if (trace.convection) downdraftCooling[k] = latentHeat * fromDraft / (mass * cp);
        }
      }
      if (!(qc[idx] > 0)) continue;
      if (autoconversionFloor !== 'none' && (floor === null ? k >= K - 2 : k > 0 && upperInterface(i, k) < floor)) continue;
      const excess = Math.max(0, qc[idx] - autoconversionThreshold);
      const converted = Math.min(qc[idx], excess * (1 - Math.exp(-autoconversionRate * dt)) + qc[idx] * (1 - Math.exp(-dt / cloudLifetime)));
      qc[idx] -= converted;
      rain += pi[i] * dSigma[k] / g * converted;
    }
    falling.evaporated = downdraft - left;
    return rain;
  }

  /*
   * The parcel of column i (see the header): fills T, p, dp, z and, from
   * the bottom up to where the search ends, the parcel's temperature in
   * Tref, and `parcel` with its source, cloud base, CAPE, inhibition and
   * top. Returns the top, or -1 when the parcel never saturates or is
   * never buoyant above its condensation level.
   */
  function diagnoseParcel(i, pi, theta, q) {
    const bottom = K - 1;
    parcel.top = -1; parcel.base = -1; parcel.cape = 0; parcel.inhibition = 0;
    for (let k = 0; k < K; k++) {
      const idx = k * C + i;
      T[k] = theta[idx] * exnerLayer[idx];
      p[k] = pi[i] * sigmaMid[k];
      dp[k] = pi[i] * dSigma[k];
      z[k] = geopotential[idx] / g;
    }
    const depth = boundaryDepth && boundaryParcel ? boundaryDepth[i] : -Infinity;
    let weight = 0, heat = 0, water = 0, source = bottom;
    for (let k = bottom; k >= 0; k--) {
      if (k < bottom && !(upperInterface(i, k + 1) < depth) && !(p[k] >= pi[i] - parcelDepth)) break;
      const idx = k * C + i;
      weight += dSigma[k]; heat += dSigma[k] * theta[idx]; water += dSigma[k] * Math.max(0, q[idx]);
      source = k;
    }
    const thetaP = heat / weight, qP = water / weight;
    parcel.source = source; parcel.theta = thetaP; parcel.q = qP;
    const T0 = thetaP * exnerLayer[bottom * C + i];
    const lcl = liftingCondensationLevel(T0, qP, p[bottom], kappa);
    if (!lcl || lcl.pressure < p[0]) return -1;
    let base = 0;
    for (let k = bottom; k >= 0; k--) if (pi[i] * levels[k] < lcl.pressure) { base = k; break; }
    parcel.base = base; parcel.lclPressure = lcl.pressure;
    let height = z[bottom] + cp * (T0 - lcl.temperature) / g;
    let energy = cp * lcl.temperature + g * height + latentHeat * qP;
    parcel.energy = energy;
    let temperature = lcl.temperature, free = false, top = -1, cape = 0, inhibition = 0;
    for (let k = bottom; k >= 0; k--) {
      const idx = k * C + i, air = Math.max(0, q[idx]);
      const saturated = p[k] < lcl.pressure;
      let vapour = qP;
      if (!saturated) Tref[k] = thetaP * exnerLayer[idx];
      else {
        const environment = cp * T[k] + g * z[k] + latentHeat * air;
        energy = environment + (energy - environment) * Math.exp(-entrainmentRate * (z[k] - height));
        height = z[k];
        temperature = saturatedTemperature(energy - g * z[k], p[k], temperature);
        Tref[k] = temperature;
        vapour = saturationHumidity(temperature, p[k]);
      }
      if (k >= source) continue;
      const virtual = virtualBuoyancy ? VIRTUAL_FACTOR : 0;
      const work = R * (Tref[k] * (1 + virtual * vapour) - T[k] * (1 + virtual * air)) * dp[k] / p[k];
      if (!free && saturated && work > 0) free = true;
      if (!free) { if (work < 0) inhibition -= work; }
      else if (work > 0) { cape += work; top = k; }
      if (saturated && T[k] - Tref[k] > 10) break;
    }
    parcel.top = top; parcel.cape = cape; parcel.inhibition = inhibition;
    return top;
  }

  /*
   * The estimated inversion strength of column i as the radiation's deck
   * finds it (radiation.module.js), from the T, p and z that
   * diagnoseParcel filled.
   */
  function inversionStrength(i, theta, q) {
    const lo = K - 1, hi = stabilityLayer, air = q[lo * C + i];
    const lcl = air > 0 ? liftingCondensationLevel(T[lo], air, p[lo], kappa) : null;
    const base = lcl ? Math.max(0, cp * (T[lo] - lcl.temperature) / g) : 0;
    const mean = 0.5 * (T[lo] + T[hi]), qs = saturationHumidity(mean, 85000);
    const lapse = g / cp * (1 + latentHeat * qs / (R * mean)) / (1 + latentHeat * latentHeat * qs / (cp * R_VAPOR * mean * mean));
    return theta[hi * C + i] - theta[lo * C + i] - (g / cp - lapse) * (z[hi] - z[lo] - base);
  }

  function referenceProfile(i, pi, theta, q) {
    const top = diagnoseParcel(i, pi, theta, q);
    if (top >= 0) for (let k = top; k <= parcel.base; k++) qref[k] = referenceHumidity * saturationHumidity(Tref[k], p[k]);
    return top;
  }

  /*
   * Convection in column i over dt (see the header), which also moves
   * its activity. Returns the rain it produces less the detrained share
   * (kg/m²); falling.base and falling.downdraft carry its cloud base and
   * the share its downdraft offers to the subcloud layers.
   */
  function convectColumn(i, pi, theta, q, dt, qc = null) {
    falling.base = K; falling.downdraft = 0; parcel.fired = false;
    const open = deckVeto && deckGate !== null ? ramp((DECK_CLOSED - deckGate[i]) / (DECK_CLOSED - ACTIVITY_UNDECIDED)) : 1;
    const top = open > 0 ? diagnoseParcel(i, pi, theta, q) : -1;
    const pass = top >= 0 ? open * ramp(0.5 + (parcel.cape - capeThreshold) / Math.max(1, capeThreshold)) * ramp(0.5 + (inhibitionThreshold - parcel.inhibition) / Math.max(1, inhibitionThreshold)) : 0;
    const now = activityMemory > 0 ? activity[i] - (pass - activity[i]) * Math.expm1(-dt / activityMemory) : pass;
    activity[i] = now;
    const shallow = top > 0 && p[top] > shallowTop;
    const firing = now > ACTIVITY_UNDECIDED || (now === ACTIVITY_UNDECIDED && pass > ACTIVITY_UNDECIDED);
    let vent = !massFlux && shallow && shallowCape !== null ? open * ramp(0.5 + (parcel.cape - shallowCape) / Math.max(1, shallowCape)) * ramp(0.5 + (shallowInhibition - parcel.inhibition) / Math.max(1, shallowInhibition)) : 0;
    if (vent > 0 && shallowStability !== null && inversionStrength(i, theta, q) > shallowStability) vent = 0;
    const acting = firing && !(massFlux && shallow);
    if (top < 0 || !(acting || vent > 0)) return 0;
    parcel.fired = acting;
    const base = adjustFrom === 'surface' ? K - 1 : parcel.base;
    const mixing = shallow && shallowReference === 'mixingLine', raining = !shallow || shallowRain;
    if (mixing) {
      const above = top - 1, aboveQ = Math.max(0, q[above * C + i]);
      const aboveEnergy = cp * T[above] + g * z[above] + latentHeat * aboveQ;
      const span = parcel.lclPressure - p[above];
      for (let k = top; k <= base; k++) {
        const chi = Math.min(1, Math.max(0, (parcel.lclPressure - p[k]) / span));
        const energy = parcel.energy + chi * (aboveEnergy - parcel.energy) - g * z[k];
        const water = parcel.q + chi * (aboveQ - parcel.q);
        let t = (energy - latentHeat * water) / cp;
        const most = shallowHumidity * saturationHumidity(t, p[k]);
        if (water > most) {
          t = saturatedTemperature(energy, p[k], t, shallowHumidity);
          qref[k] = shallowHumidity * saturationHumidity(t, p[k]);
        } else qref[k] = water;
        Tref[k] = t;
      }
    } else {
      for (let k = top; k <= base; k++) qref[k] = referenceHumidity * saturationHumidity(Tref[k], p[k]);
    }
    let heating = 0, drying = 0, depth = 0;
    for (let k = top; k <= base; k++) {
      heating += cp * (Tref[k] - T[k]) * dp[k];
      drying -= (qref[k] - q[k * C + i]) * dp[k];
      depth += dp[k];
    }
    if (!mixing && raining && heating <= 0) return 0;
    const rate = dt / relaxationTime * (acting ? 1 : vent);
    let rain = 0;
    if (raining && drying > 0) {
      const shift = (latentHeat * drying - heating) / (cp * depth);
      for (let k = top; k <= base; k++) Tref[k] += shift;
      rain = drying / g * rate;
    } else {
      const shiftQ = drying / depth, shiftT = -heating / (cp * depth);
      for (let k = top; k <= base; k++) { qref[k] += shiftQ; Tref[k] += shiftT; }
    }
    for (let k = top; k <= base; k++) {
      const idx = k * C + i;
      theta[idx] += (Tref[k] - T[k]) * rate / exnerLayer[idx];
      q[idx] += (qref[k] - q[idx]) * rate;
    }
    if (rain <= 0) return 0;
    let detrained = 0;
    if (qc !== null && detrainment > 0) {
      let anvilMass = 0, anvilBottom = top;
      for (let k = top; k <= base && (k === top || anvilMass < anvilDepth); k++) { anvilMass += dp[k]; anvilBottom = k; }
      detrained = detrainment * rain;
      for (let k = top; k <= anvilBottom; k++) qc[k * C + i] += detrained * g / anvilMass;
    }
    falling.base = parcel.base;
    falling.downdraft = downdraftEvaporation * (rain - detrained);
    return rain - detrained;
  }

  function clearCumulus(i) {
    cumulusBaseFlux[i] = 0; cumulusTop[i] = 0;
    for (let k = 0; k < K; k++) { cumulusCover[k * C + i] = 0; cumulusWater[k * C + i] = 0; }
  }

  /*
   * The temperature and condensate of plume air of liquid-water static
   * energy `energy` and total water `water` at `height` and `pressure`.
   */
  function plumeState(energy, water, height, pressure, guess) {
    const dry = (energy - g * height) / cp;
    if (!(water > saturationHumidity(dry, pressure))) { plume.T = dry; plume.liquid = 0; return; }
    const t = saturatedTemperature(energy - g * height + latentHeat * water, pressure, Math.max(dry, guess));
    plume.T = t;
    plume.liquid = Math.max(0, water - saturationHumidity(t, pressure));
  }

  /*
   * The shallow cumulus mass flux of column i over dt (see the header).
   * Returns its rain (kg/m²).
   */
  function cumulusColumn(i, pi, theta, q, qc, dt) {
    const bottom = K - 1;
    cumulus.top = -1; cumulus.inhibition = 0; cumulus.baseFlux = 0; cumulus.velocity = 0;
    clearCumulus(i);
    const open = deckVeto && deckGate !== null ? ramp((DECK_CLOSED - deckGate[i]) / (DECK_CLOSED - ACTIVITY_UNDECIDED)) : 1;
    const buoyancy = surfaceBuoyancy ? surfaceBuoyancy[i] : 0;
    if (!(open > 0) || !(buoyancy > 0) || !boundaryDepth) return 0;
    for (let k = 0; k < K; k++) {
      const idx = k * C + i, cloud = qc ? Math.max(0, qc[idx]) : 0;
      T[k] = theta[idx] * exnerLayer[idx];
      p[k] = pi[i] * sigmaMid[k];
      dp[k] = pi[i] * dSigma[k];
      z[k] = geopotential[idx] / g;
      envS[k] = cp * T[k] + g * z[k] - latentHeat * cloud;
      envQ[k] = Math.max(0, q[idx]) + cloud;
    }
    const depth = boundaryDepth[i];
    let mass = 0, energy = 0, water = 0, source = bottom;
    for (let k = bottom; k >= 0; k--) {
      if (k < bottom && ((!(upperInterface(i, k + 1) < depth) && !(p[k] >= pi[i] - cumulusSourceDepth)) || !(p[k] > shallowTop))) break;
      mass += dp[k]; energy += dp[k] * envS[k]; water += dp[k] * envQ[k];
      source = k;
    }
    const lowest = cumulusSource === 'lowest';
    const sourceS = lowest ? envS[bottom] : energy / mass, sourceQ = lowest ? envQ[bottom] : water / mass;
    const lcl = liftingCondensationLevel((sourceS - g * z[bottom]) / cp, sourceQ, p[bottom], kappa);
    if (!lcl || !(lcl.pressure > shallowTop)) return 0;
    cumulus.source = source; cumulus.lclPressure = lcl.pressure;
    const virtual = virtualBuoyancy ? VIRTUAL_FACTOR : 0, loading = virtualBuoyancy ? 1 : 0;
    let s = sourceS, w = sourceQ, inhibition = 0, cloudy = false, top = -1, guess = 0;
    plumeS[source] = s; plumeQ[source] = w;
    for (let k = source - 1; k >= 0; k--) {
      if (!(p[k] > shallowTop)) { if (cloudy) top = k + 1; break; }
      const below = upperInterface(i, k + 1), above = upperInterface(i, k);
      const mixes = pi[i] * levels[k + 1] <= lcl.pressure, epsilon = mixes ? cumulusEntrainment : 0;
      const half = Math.exp(-epsilon * (z[k] - below));
      const midS = envS[k] + (s - envS[k]) * half, midQ = envQ[k] + (w - envQ[k]) * half;
      plumeState(midS, midQ, z[k], p[k], guess);
      guess = plume.T;
      plumeLiquid[k] = plume.liquid;
      const air = Math.max(0, q[k * C + i]), cloud = qc ? Math.max(0, qc[k * C + i]) : 0;
      const work = R * (plume.T * (1 + virtual * (midQ - plume.liquid) - loading * plume.liquid) - T[k] * (1 + virtual * air - loading * cloud)) * dp[k] / p[k];
      if (plume.liquid > 0) cloudy = true;
      if (!cloudy) { if (work < 0) inhibition -= work; } else if (!(work > 0)) { top = k; break; }
      const full = Math.exp(-epsilon * (above - below));
      s = envS[k] + (s - envS[k]) * full; w = envQ[k] + (w - envQ[k]) * full;
      plumeRain[k] = 0;
      if (cumulusRain !== null) {
        plumeState(s, w, above, pi[i] * levels[k], guess);
        const excess = plume.liquid - cumulusRain;
        if (excess > 0) { w -= excess; s += latentHeat * excess; plumeRain[k] = excess; }
      }
      plumeS[k] = s; plumeQ[k] = w;
      plumeGrowth[k] = Math.exp((epsilon - (mixes ? cumulusDetrainment : 0)) * (above - below));
    }
    cumulus.inhibition = inhibition;
    if (top < 0) return 0;
    const convective = Math.cbrt(buoyancy * Math.max(0, depth - z[bottom]));
    const velocity = Math.max(convective, cumulusFriction * (frictionVelocity ? frictionVelocity[i] : 0));
    cumulus.velocity = velocity;
    if (!(velocity > 0)) return 0;
    let base = open * cumulusClosure * lcl.pressure / (R * lcl.temperature) * velocity * Math.exp(-inhibition / (velocity * velocity));
    base = Math.min(base, cumulusBoundaryLoss * mass / (g * dt));
    if (!(base > CUMULUS_FLOOR)) return 0;
    cumulusFlux.fill(0);
    let sourceBelow = 0;
    for (let j = bottom; j > source; j--) { sourceBelow += dp[j]; cumulusFlux[j] = base * sourceBelow / mass; }
    cumulusFlux[source] = base;
    for (let k = source - 1; k > top; k--) cumulusFlux[k] = cumulusFlux[k + 1] * plumeGrowth[k];
    cumulusFlux[top + 1] *= cumulusOvershoot;
    let scale = 1;
    for (let k = top; k <= bottom; k++) {
      const courant = Math.max(cumulusFlux[k], cumulusFlux[k + 1]) * g * dt / dp[k];
      if (courant * scale > 1) scale = 1 / courant;
    }
    for (let j = top + 1; j <= bottom; j++) cumulusFlux[j] *= scale;
    let belowMass = 0, belowS = 0, belowQ = 0;
    fluxS[K] = 0; fluxQ[K] = 0; fluxS[top] = 0; fluxQ[top] = 0;
    for (let j = bottom; j > top; j--) {
      let upS = plumeS[j], upQ = plumeQ[j];
      if (j > source) {
        belowMass += dp[j]; belowS += dp[j] * envS[j]; belowQ += dp[j] * envQ[j];
        upS = lowest ? sourceS : belowS / belowMass; upQ = lowest ? sourceQ : belowQ / belowMass;
      }
      fluxS[j] = cumulusFlux[j] * (upS - envS[j - 1]);
      fluxQ[j] = cumulusFlux[j] * (upQ - envQ[j - 1]);
    }
    let rain = 0;
    for (let k = top; k <= bottom; k++) {
      const idx = k * C + i, per = g * dt / dp[k];
      let dS = (fluxS[k + 1] - fluxS[k]) * per, dQ = (fluxQ[k + 1] - fluxQ[k]) * per;
      if (cumulusRain !== null && k > top && k < source && plumeRain[k] > 0) {
        const fallen = cumulusFlux[k] * plumeRain[k] * dt;
        rain += fallen; dQ -= fallen * g / dp[k]; dS += latentHeat * fallen * g / dp[k];
      }
      theta[idx] += dS / (cp * exnerLayer[idx]);
      q[idx] += dQ;
    }
    for (let k = top; k < source; k++) {
      if (!(plumeLiquid[k] > 0)) continue;
      const idx = k * C + i;
      cumulusCover[idx] = Math.min(1, 0.5 * (cumulusFlux[k] + cumulusFlux[k + 1]) * R * T[k] / (p[k] * cumulusUpdraft));
      cumulusWater[idx] = plumeLiquid[k];
    }
    cumulus.top = top; cumulus.baseFlux = base * scale;
    cumulusBaseFlux[i] = base * scale; cumulusTop[i] = pi[i] * levels[top];
    return rain;
  }

  function fillColumn(i, pi, q) {
    if (!q) return;
    for (let k = 0; k < K - 1; k++) {
      const idx = k * C + i;
      if (q[idx] < 0) {
        q[idx + C] += q[idx] * dSigma[k] / dSigma[k + 1];
        q[idx] = 0;
      }
    }
    const bottom = (K - 1) * C + i;
    if (q[bottom] < 0) {
      budget.lost -= mesh.areaCell[i] * pi[i] * dSigma[K - 1] / g * q[bottom];
      q[bottom] = 0;
    }
  }

  function adjust(state, iFrom, iTo, dt) {
    const [pi, theta, , , q, qc] = state;
    for (let i = iFrom; i < iTo; i++) {
      core.diagnoseColumn(i, pi, theta, q, qc);
      const traced = trace.convection || trace.largeScale;
      if (traced) mark(i, theta);
      condenseColumn(i, pi, theta, q, qc);
      if (traced) charge(trace.largeScale, i, theta);
      let produced = convectColumn(i, pi, theta, q, dt, qc);
      if (massFlux) {
        if (parcel.fired && !cumulusWithDeep) clearCumulus(i);
        else produced += cumulusColumn(i, pi, theta, q, qc, dt);
      }
      if (traced) charge(trace.convection, i, theta);
      if (cumulusBaseFlux[i] > 0) condenseColumn(i, pi, theta, q, qc);
      const rained = autoconvertColumn(i, pi, theta, q, qc, dt, falling.downdraft, falling.base);
      const convected = produced - falling.evaporated;
      if (traced) {
        charge(trace.largeScale, i, theta);
        if (trace.convection) {
          for (let k = 0; k < K; k++) {
            const idx = k * C + i;
            trace.convection[idx] -= downdraftCooling[k];
            if (trace.largeScale) trace.largeScale[idx] += downdraftCooling[k];
          }
        }
      }
      fillColumn(i, pi, q);
      fillColumn(i, pi, qc);
      precipitation[i] += rained + convected;
      convectivePrecipitation[i] += convected;
      largeScalePrecipitation[i] += rained;
      rain[i] = rained + convected;
      budget.condensation += mesh.areaCell[i] * rained;
      budget.convection += mesh.areaCell[i] * convected;
    }
  }

  function readRain(interval) {
    const scale = 86400 / interval;
    for (let i = 0; i < C; i++) {
      convectiveRain[i] = scale * convectivePrecipitation[i];
      largeScaleRain[i] = scale * largeScalePrecipitation[i];
    }
  }

  function columnWater(pi, q, i) {
    let water = 0;
    for (let k = 0; k < K; k++) water += pi[i] * dSigma[k] / g * q[k * C + i];
    return water;
  }

  return {
    adjust, condenseColumn, autoconvertColumn, convectColumn, cumulusColumn, fillColumn, referenceProfile, diagnoseParcel, inversionStrength, columnWater, readRain,
    precipitation, rain, convectivePrecipitation, largeScalePrecipitation, convectiveRain, largeScaleRain, activity, convectiveActivity: activity, budget, latentHeat, trace, parcel, falling,
    cumulus, cumulusCover, cumulusWater, cumulusBaseFlux, cumulusTop, cumulusFlux,
    shared: { precipitation: precipBuffer, rain: rainBuffer, convectivePrecipitation: convectiveBuffer, largeScalePrecipitation: largeScaleBuffer, convectiveActivity: activityBuffer, cumulusCover: cumulusCoverBuffer, cumulusWater: cumulusWaterBuffer, cumulusBaseFlux: baseFluxBuffer, cumulusTop: cumulusTopBuffer },
    reference: { T: Tref, q: qref },
  };
}
