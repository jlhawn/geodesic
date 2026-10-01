import { cellVector } from '../dynamics/operators.module.js';
import { VIRTUAL_FACTOR } from '../dynamics/sigmaCore.module.js';
import { saturationHumidity, DECK_OPEN, DECK_CLOSED, LATENT_HEAT, R_VAPOR } from './moist.module.js';
import { SEA_DRAG } from './surface.module.js';

/*
 * A diffusive planetary boundary layer in the manner of Troen and Mahrt
 * (1986), as used by the simple moist GCMs that this model follows.
 * `diagnose` finds, per cell, the boundary-layer top as the height where
 * the bulk Richardson number of the lowest layer's virtual potential
 * temperature and wind, with the convective floor 100 u*², first exceeds
 * richardsonCritical, and lays the K-profile κ u* z (1 − z/h)² over the
 * layer interfaces below it (u* from the bulk drag on the lowest layer's
 * wind, with the gustiness floor). Where the surface is warmer than the
 * lowest layer the profile's velocity scale is the unstable one of
 * Holtslag and Boville (1993), u* (1 − 15 ζ)^¼ with ζ = 0.1 h / L from
 * the bulk surface buoyancy flux (virtual, with the saturation humidity
 * of a sea surface; dry over `land`) and floored at −2, so a convective
 * marine boundary layer mixes momentum down to the surface; stable
 * columns keep the neutral profile. The interface coefficients ρK/Δz are
 * kept in a shared array so the cell units of the adjust phase can mix
 * θ, q and qc down each column and the edge units can mix the normal
 * velocity down each edge, both by implicit Euler on the same
 * tridiagonal system, which conserves each column's mass-weighted
 * total exactly. Nothing mixes above the boundary-layer top; the search
 * stops at searchTop in σ. The surface fluxes and drag remain explicit
 * sources on the lowest layer, which the diffusion then spreads upward.
 *
 * A stratocumulus deck mixes its layer from cloud top, which the bulk
 * Richardson number of the surface-driven search does not see. With
 * `deckTop`, per cell the inversion height of the mixed-layer deck in the
 * height coordinate of `depth` where the deck ran this step and 0 where
 * it did not (radiation.module.js's mlmTop), the K-profile of such a
 * cell is laid over max(depth, deckTop) instead of depth — the same
 * profile, shape and velocity scale, for the deeper layer — so the
 * deck's layer is mixed through to its inversion. `depth` itself stays
 * the Richardson depth: the deck starts from it and relaxes toward it.
 *
 * The K-profile vanishes at h, so its top entrains nothing; an explicit
 * entrainment flux closes it (`entrainment`). Where the surface buoyancy
 * flux B0 (the one above, positive upward) is positive,
 *   w_e = o (1 − s) min(cap, (A B0 + A_s r u*³ / h) / max(Δb, bMin)),
 * Δb = g Δθv / θv between the first layer whose base lies above h and
 * the mass mean of the layers below it: the buoyancy and friction-velocity
 * sources of Tennekes's (1973) inversion model with Driedonks's (1982)
 * constants A 0.2 and A_s 5, bMin 0.015 m/s² (a jump of about 0.5 K) and
 * the cap 0.05 m/s. r = min(1, B0 / `shearOnset`)
 * (5·10⁻⁵ m²/s³, about 2 W/m² of virtual heat flux) brings the shear
 * term in continuously as the surface turns unstable. The stratocumulus
 * regime belongs to the mixed-layer deck: o is the deck's opening, 1 at a
 * gate `deckGate` (radiation.mlmGate) of one half or less and 0 at
 * DECK_CLOSED, as convection takes it, and s the radiation's stratiform
 * share `stratiform` (radiation.stratiform), 0 below an estimated
 * inversion strength of 8 K and 1 above 12 K. It enters as the
 * coefficient ρ w_e of that layer's lower interface, so the exchange of
 * θ, q, qc and momentum across h goes through the same conservative
 * implicit solve; `entrainment` keeps each cell's w_e (m/s). The depth
 * stays diagnostic. All of the above is `turbulence` 'dry'.
 *
 * `turbulence` 'moist' (the default) is a moist closure after Lock et al.
 * (2000), simplified for a coarse column. It mixes the liquid-water
 * potential temperature θ_l = θ − L q_c/(c_p Π) and the total water
 * q_t = q + q_c (and momentum) with the eddy diffusivity; each layer the
 * solve touches leaves it with θ = θ_l, q = q_t and no cloud water, and
 * the saturation adjustment that follows returns the cloud, so a
 * well-mixed layer forms its stratus at its top. Untouched layers keep
 * their values exactly. Two profiles add:
 *  - the surface-driven K-profile above, over h_s: the Richardson depth,
 *    or where B0 > 0 the top of a surface parcel (the lowest layer's θ_l
 *    with Holtslag and Boville's excess cloudTop.excess (8.5) B0 θ_v/(g w_m),
 *    w_m³ = u*³ + 0.6 B0 h, and its q_t) rising with its condensate in
 *    equilibrium until its θ_v falls short of the layer's by more than
 *    cloudTop.tolerance (0.5 K, for a layer that straddles the
 *    inversion), where it has
 *    condensed and stops no more than cloudTop.cumulusDepth (400 m) above
 *    the base of the layer where it first saturated and below
 *    cloudTop.maximumHeight, and lies above the Richardson depth: a
 *    stratocumulus-capped layer, whose turbulence reaches its inversion
 *    (Lock et al.'s parcel test); a parcel that rises further is cumulus,
 *    left to the plume, and its column keeps the Richardson depth;
 *  - a cloud-top-driven profile where the lowest run of cloudy layers
 *    (q_c above cloudTop.threshold 10⁻⁶) tops out below
 *    cloudTop.maximumHeight (3 km) and cools: ΔF, the longwave cooling
 *    summed over the run's layers (`longwave`, the radiation's per-layer
 *    longwave heating, W/m²), gives V³ = (g/θ_v) ΔF/(ρ c_p) z_ml and
 *    K = cloudTop.profile κ V z_ml x² (1 − x)^½, x = (z − z_b)/z_ml, over
 *    z_b < z < h_c (Lock et al.'s 0.85), h_c the cloud top's upper
 *    interface and z_ml = h_c − z_b. z_b is where a parcel of the cloud
 *    top's θ_l less cloudTop.perturbation (0.2 K) and q_t, descending
 *    with its condensate in equilibrium, stops being negatively buoyant
 *    in θ_v (the lower interface of the last layer it sinks through; 0
 *    when it reaches the lowest layer).
 * Regimes (`regime`, REGIME): stable (no cloud-top layer, B0 ≤ 0);
 * surface-driven (no cloud-top layer, B0 > 0); decoupled (a cloud-top
 * layer with z_b above h_s: the subcloud layer mixes from the surface
 * and the plume carries air into the cloud layer); coupled (z_b at the
 * surface or at or below h_s: the surface profile reaches h_c as well).
 * `depth` is h_c in the coupled regime and h_s otherwise; `mixingTop` is
 * the top of all mixing, h_c wherever there is a cloud-top layer.
 * Entrainment at the inversion follows the closure:
 *   w_e = min(cap, (A (w_s³ + V³) + A_s r u*³) / (h max(Δb, bMin)))
 * across the interface above the mixed layer, Δb = g Δθ_v/θ_v between
 * the layer above, or the one above that where its θ_v is higher (with
 * `entrainment.jumpLayers` 2, the default: the layer above is the
 * inversion's own grid layer and holds part of the jump, as Lock et
 * al. take it), and the mass mean of the mixed layers below, w_s³ =
 * B0 h for a surface-driven top (coupled or clear; 0 for a decoupled
 * cloud top, where h = z_ml), and A Nicholls and Turton's efficiency as
 * the mixed-layer deck takes it: a_1 [1 + a_2 χ* (1 − Δθ_vs/Δθ_v)] with
 * a_1 `efficiency` (0.2), a_2 `evaporativeEnhancement` (25), at most
 * `maximumEfficiency` (1), χ* and Δθ_vs from the cloudy top layer's
 * state and the jumps in θ_l and q_t, and a_1 alone where the top layer
 * holds no cloud, which is the dry scheme's form above. A decoupled
 * column also entrains across its surface-driven top at that form. Nothing is
 * tapered by the deck unless `entrainment.taper`. The diagnosis keeps
 * per cell the cloud-top cooling ΔF (`cloudTopCooling`, W/m²), V
 * (`radiativeVelocity`) and the decoupling height z_b (`decoupling`, in
 * the coordinate of `depth`; 0 where coupled or without a cloud top).
 */
export const ENTRAINMENT_DEFAULTS = { efficiency: 0.2, shear: 5, cap: 0.05, jumpFloor: 0.015, shearOnset: 5e-5, evaporativeEnhancement: 25, maximumEfficiency: 1, taper: false, jumpLayers: 2 };
export const CLOUD_TOP_DEFAULTS = { threshold: 1e-6, maximumHeight: 3000, perturbation: 0.2, profile: 0.85, excess: 8.5, tolerance: 0.5, cumulusDepth: 400 };
export const REGIME = { STABLE: 0, SURFACE: 1, DECOUPLED: 2, COUPLED: 3 };

export function createBoundaryLayer(mesh, core, {
  dragCoefficient = SEA_DRAG, dragCoefficients = null, gustiness = 3, richardsonCritical = 0.5, vonKarman = 0.4, searchTop = 0.5, stability = true, land = null, deckTop = null, deckGate = null, stratiform = null, buffers = null,
  entrainment: entrainmentOptions = {}, turbulence = 'moist', cloudTop: cloudTopOptions = {}, longwave = null, latentHeat = LATENT_HEAT,
} = {}) {
  if (turbulence !== 'moist' && turbulence !== 'dry') throw new Error(`turbulence must be 'moist' or 'dry', not ${turbulence}`);
  const { efficiency, shear, cap, jumpFloor, shearOnset, evaporativeEnhancement, maximumEfficiency, taper: tapered, jumpLayers } = { ...ENTRAINMENT_DEFAULTS, ...entrainmentOptions };
  const { threshold: cloudThreshold, maximumHeight: cloudTopHeight, perturbation, profile, excess: excessCoefficient, tolerance, cumulusDepth } = { ...CLOUD_TOP_DEFAULTS, ...cloudTopOptions };
  const moistScheme = turbulence === 'moist';
  const { K, C, E, dSigma, sigmaMid, R, g, cp, kappa, exnerLayer, geopotential } = core.diagnostics;
  const thetaV = core.arrays.thetaV;
  const { cellsOnEdge } = mesh;
  const bottom = K - 1;
  let kTop = 0;
  while (kTop < bottom && sigmaMid[kTop] <= searchTop) kTop++;
  const n = K - kTop;
  const mixingBuffer = buffers && buffers.mixing ? buffers.mixing : new SharedArrayBuffer(8 * K * C);
  const depthBuffer = buffers && buffers.depth ? buffers.depth : new SharedArrayBuffer(8 * C);
  const buoyancyBuffer = buffers && buffers.buoyancyFlux ? buffers.buoyancyFlux : new SharedArrayBuffer(8 * C);
  const frictionBuffer = buffers && buffers.friction ? buffers.friction : new SharedArrayBuffer(8 * C);
  const entrainmentBuffer = buffers && buffers.entrainment ? buffers.entrainment : new SharedArrayBuffer(8 * C);
  const mixing = new Float64Array(mixingBuffer);
  const depth = new Float64Array(depthBuffer);
  const buoyancyFlux = new Float64Array(buoyancyBuffer), friction = new Float64Array(frictionBuffer);
  const entrainmentVelocity = new Float64Array(entrainmentBuffer);
  const extra = (name) => (buffers && buffers[name] ? buffers[name] : new SharedArrayBuffer(8 * C));
  const regimeBuffer = extra('regime'), mixingTopBuffer = extra('mixingTop'), coolingBuffer = extra('cloudTopCooling'), velocityBuffer = extra('radiativeVelocity'), decouplingBuffer = extra('decoupling');
  const regime = new Float64Array(regimeBuffer), mixingTop = new Float64Array(mixingTopBuffer), cloudTopCooling = new Float64Array(coolingBuffer), radiativeVelocity = new Float64Array(velocityBuffer), decoupling = new Float64Array(decouplingBuffer);
  const thetaL = new Float64Array(K), totalWater = new Float64Array(K);
  const entraining = efficiency > 0 || shear > 0;
  const vector = new Float64Array(3 * C), bottomVector = new Float64Array(3 * C);
  const speed = new Float64Array(C), riPrev = new Float64Array(C), zPrev = new Float64Array(C);
  const found = new Uint8Array(C);
  const upper = new Float64Array(K), lower = new Float64Array(K), gain = new Float64Array(K), rhs = new Float64Array(K), mass = new Float64Array(K);

  function diagnose(state, iFrom = 0, iTo = C) {
    const [pi, theta, u, surfaceT, q = null, qc = null] = state;
    for (let i = iFrom; i < iTo; i++) core.diagnoseColumn(i, pi, theta, q, qc);
    cellVector(mesh, u.subarray(bottom * E, K * E), bottomVector, iFrom, iTo);
    for (let i = iFrom; i < iTo; i++) {
      speed[i] = Math.hypot(bottomVector[3 * i], bottomVector[3 * i + 1], bottomVector[3 * i + 2]);
      friction[i] = Math.sqrt(dragCoefficients ? dragCoefficients[i] : dragCoefficient) * Math.max(speed[i], gustiness);
      found[i] = 0;
      riPrev[i] = 0;
      zPrev[i] = geopotential[bottom * C + i] / g;
      depth[i] = zPrev[i];
    }
    for (let k = bottom - 1; k >= kTop; k--) {
      cellVector(mesh, u.subarray(k * E, (k + 1) * E), vector, iFrom, iTo);
      for (let i = iFrom; i < iTo; i++) {
        if (found[i]) continue;
        const idx = k * C + i, base = bottom * C + i;
        const z = geopotential[idx] / g, zb = geopotential[base] / g;
        const du = vector[3 * i] - bottomVector[3 * i], dv = vector[3 * i + 1] - bottomVector[3 * i + 1], dw = vector[3 * i + 2] - bottomVector[3 * i + 2];
        const shear = du * du + dv * dv + dw * dw + 100 * friction[i] * friction[i];
        const ri = g * (thetaV[idx] - thetaV[base]) * (z - zb) / (thetaV[base] * shear);
        if (ri > richardsonCritical) {
          depth[i] = zPrev[i] + (z - zPrev[i]) * (richardsonCritical - riPrev[i]) / (ri - riPrev[i]);
          found[i] = 1;
        } else {
          riPrev[i] = ri;
          zPrev[i] = z;
          if (k === kTop) depth[i] = z;
        }
      }
    }
    if (moistScheme) {
      for (let i = iFrom; i < iTo; i++) moistColumn(i, pi, theta, surfaceT, q, qc);
      return;
    }
    for (let i = iFrom; i < iTo; i++) {
      const zb = geopotential[bottom * C + i] / g, h = (deckTop && deckTop[i] > 0 ? Math.max(depth[i], deckTop[i]) : depth[i]) - zb;
      for (let k = kTop; k < K; k++) mixing[k * C + i] = 0;
      entrainmentVelocity[i] = 0;
      const base = bottom * C + i;
      const moisture = q && !(land && land[i]) ? 0.61 * theta[base] * (saturationHumidity(surfaceT[i], pi[i]) - q[base]) : 0;
      const buoyancy = g / theta[base] * (dragCoefficients ? dragCoefficients[i] : dragCoefficient) * Math.max(speed[i], gustiness) * (surfaceT[i] * Math.pow(sigmaMid[bottom], kappa) / exnerLayer[base] - theta[base] + moisture);
      buoyancyFlux[i] = buoyancy;
      if (h <= 0) continue;
      let scale = friction[i];
      if (stability && buoyancy > 0) scale = friction[i] * Math.pow(1 - 15 * Math.max(-2, -0.1 * h * vonKarman * buoyancy / friction[i] ** 3), 0.25);
      let entrainK = -1;
      for (let k = kTop; k < bottom; k++) {
        const idx = k * C + i, below = idx + C;
        const zAbove = geopotential[idx] / g, zBelow = geopotential[below] / g;
        const z = 0.5 * (zAbove + zBelow) - zb;
        if (z >= h) { entrainK = k; continue; }
        const diffusivity = vonKarman * scale * z * (1 - z / h) ** 2;
        const rhoAbove = pi[i] * sigmaMid[k] / (R * theta[idx] * exnerLayer[idx]);
        const rhoBelow = pi[i] * sigmaMid[k + 1] / (R * theta[below] * exnerLayer[below]);
        mixing[idx] = 0.5 * (rhoAbove + rhoBelow) * diffusivity / (zAbove - zBelow);
      }
      const open = deckGate ? Math.min(1, Math.max(0, (DECK_CLOSED - deckGate[i]) / (DECK_CLOSED - DECK_OPEN))) : 1;
      const taper = open * (stratiform ? 1 - stratiform[i] : 1);
      if (entraining && entrainK >= kTop && buoyancy > 0 && taper > 0) {
        const idx = entrainK * C + i, below = idx + C;
        let weight = 0, sum = 0;
        for (let k = entrainK + 1; k < K; k++) { weight += dSigma[k]; sum += dSigma[k] * thetaV[k * C + i]; }
        const mean = sum / weight, jump = g * (thetaV[idx] - mean) / mean;
        const onset = shearOnset > 0 ? Math.min(1, buoyancy / shearOnset) : 1;
        const velocity = taper * Math.min(cap, (efficiency * buoyancy + shear * onset * friction[i] ** 3 / h) / Math.max(jump, jumpFloor));
        entrainmentVelocity[i] = velocity;
        mixing[idx] = 0.5 * (pi[i] * sigmaMid[entrainK] / (R * theta[idx] * exnerLayer[idx]) + pi[i] * sigmaMid[entrainK + 1] / (R * theta[below] * exnerLayer[below])) * velocity;
      }
    }
  }

  function surfaceBuoyancy(i, pi, theta, surfaceT, q) {
    const base = bottom * C + i;
    const moisture = q && !(land && land[i]) ? 0.61 * theta[base] * (saturationHumidity(surfaceT[i], pi[i]) - q[base]) : 0;
    return g / theta[base] * (dragCoefficients ? dragCoefficients[i] : dragCoefficient) * Math.max(speed[i], gustiness) * (surfaceT[i] * Math.pow(sigmaMid[bottom], kappa) / exnerLayer[base] - theta[base] + moisture);
  }

  function density(pi, k, i, theta) {
    const idx = k * C + i;
    return pi[i] * sigmaMid[k] / (R * theta[idx] * exnerLayer[idx]);
  }

  function parcelVirtual(level, total, k, i, pi) {
    const ex = exnerLayer[k * C + i], liquidT = level * ex, qs = saturationHumidity(liquidT, pi[i] * sigmaMid[k]);
    if (!(total > qs)) return level * (1 + VIRTUAL_FACTOR * total);
    const slope = qs * latentHeat / (R_VAPOR * liquidT * liquidT), liquid = (total - qs) / (1 + latentHeat * slope / cp);
    return (liquidT + latentHeat * liquid / cp) / ex * (1 + VIRTUAL_FACTOR * (total - liquid) - liquid);
  }

  /*
   * w_e across interface kE into the mixed layers kE + 1 … lowest (see
   * the header); `buoyant` is w_s³ + V³, `sheared` r u*³, h the depth the
   * velocity scales are spread over.
   */
  function entrain(i, pi, theta, q, qc, kE, lowest, h, buoyant, sheared) {
    if (kE < kTop || !(h > 0)) return 0;
    let weight = 0, sumV = 0, sumL = 0, sumQ = 0;
    for (let k = kE + 1; k <= lowest; k++) {
      const idx = k * C + i;
      weight += dSigma[k]; sumV += dSigma[k] * thetaV[idx];
      sumL += dSigma[k] * (theta[idx] - latentHeat * qc[idx] / (cp * exnerLayer[idx])); sumQ += dSigma[k] * (q[idx] + qc[idx]);
    }
    const above = jumpLayers > 1 && kE > kTop && thetaV[(kE - 1) * C + i] > thetaV[kE * C + i] ? (kE - 1) * C + i : kE * C + i;
    const mean = sumV / weight, jump = g * (thetaV[above] - mean) / mean;
    let efficiencyNow = efficiency;
    const top = (kE + 1) * C + i;
    if (qc[top] > cloudThreshold && evaporativeEnhancement > 0) {
      const ex = exnerLayer[top], T = theta[top] * ex, qs = saturationHumidity(T, pi[i] * sigmaMid[kE + 1]), dqs = qs * latentHeat / (R_VAPOR * T * T);
      const gamma = latentHeat / cp * dqs, c = 1 + VIRTUAL_FACTOR * q[top] - qc[top] + (1 + VIRTUAL_FACTOR) * T * dqs;
      const jumpL = theta[above] - latentHeat * qc[above] / (cp * exnerLayer[above]) - sumL / weight, jumpQ = q[above] + qc[above] - sumQ / weight;
      const virtualJump = Math.max(thetaV[above] - mean, jumpFloor * mean / g);
      const saturatedJump = c / (1 + gamma) * jumpL + (c * latentHeat / (cp * ex * (1 + gamma)) - theta[top]) * jumpQ;
      const demand = dqs * ex * jumpL - jumpQ;
      const share = demand > 0 ? Math.min(1, qc[top] * (1 + gamma) / demand) : 1;
      efficiencyNow = Math.min(maximumEfficiency, efficiency * (1 + evaporativeEnhancement * Math.max(0, share * (1 - saturatedJump / virtualJump))));
    }
    let velocity = Math.min(cap, (efficiencyNow * buoyant + shear * sheared) / (h * Math.max(jump, jumpFloor)));
    if (tapered) {
      const open = deckGate ? Math.min(1, Math.max(0, (DECK_CLOSED - deckGate[i]) / (DECK_CLOSED - DECK_OPEN))) : 1;
      velocity *= open * (stratiform ? 1 - stratiform[i] : 1);
    }
    if (!(velocity > 0)) return 0;
    mixing[kE * C + i] += 0.5 * (density(pi, kE, i, theta) + density(pi, kE + 1, i, theta)) * velocity;
    return velocity;
  }

  function moistColumn(i, pi, theta, surfaceT, q, qc) {
    const base = bottom * C + i, zb = geopotential[base] / g;
    const interfaceZ = (k) => 0.5 * (geopotential[k * C + i] / g + geopotential[(k + 1) * C + i] / g) - zb;
    for (let k = kTop; k < K; k++) mixing[k * C + i] = 0;
    entrainmentVelocity[i] = 0;
    const buoyancy = surfaceBuoyancy(i, pi, theta, surfaceT, q);
    buoyancyFlux[i] = buoyancy;
    let surfaceDepth = depth[i] - zb;
    if (buoyancy > 0 && q && qc && cumulusDepth > 0) {
      const mixed = Math.cbrt(friction[i] ** 3 + 0.6 * buoyancy * Math.max(0, surfaceDepth));
      const level = theta[base] - latentHeat * qc[base] / (cp * exnerLayer[base]) + excessCoefficient * buoyancy * thetaV[base] / (g * mixed), total = q[base] + qc[base];
      let k = bottom - 1, condensation = -1;
      for (; k >= kTop; k--) {
        const ex = exnerLayer[k * C + i];
        if (condensation < 0 && total > saturationHumidity(level * ex, pi[i] * sigmaMid[k])) condensation = k;
        if (!(parcelVirtual(level, total, k, i, pi) + tolerance > thetaV[k * C + i])) break;
      }
      if (condensation > k && k >= kTop) {
        const parcelTop = interfaceZ(k), cloudBase = interfaceZ(condensation);
        if (parcelTop - cloudBase <= cumulusDepth && parcelTop <= cloudTopHeight && parcelTop > surfaceDepth) {
          surfaceDepth = parcelTop;
          depth[i] = zb + parcelTop;
        }
      }
    }
    let top = -1, cooling = 0;
    if (longwave && q && qc) {
      for (let k = bottom; k > kTop; k--) {
        if (interfaceZ(k - 1) > cloudTopHeight) break;
        if (qc[k * C + i] > cloudThreshold && !(qc[(k - 1) * C + i] > cloudThreshold)) { top = k; break; }
      }
      if (top >= 0) {
        for (let k = top; k <= bottom && qc[k * C + i] > cloudThreshold; k++) cooling -= longwave[k * C + i];
        if (!(cooling > 0)) { top = -1; cooling = 0; }
      }
    }
    let coupled = false, lowest = bottom, base0 = 0, cloudTopZ = 0;
    if (top >= 0) {
      const idx = top * C + i, level = theta[idx] - latentHeat * qc[idx] / (cp * exnerLayer[idx]) - perturbation, total = q[idx] + qc[idx];
      let k = top + 1;
      while (k <= bottom && parcelVirtual(level, total, k, i, pi) < thetaV[k * C + i]) k++;
      lowest = k - 1;
      base0 = k > bottom ? 0 : interfaceZ(k - 1);
      coupled = k > bottom || base0 <= surfaceDepth;
      if (coupled) base0 = 0;
      cloudTopZ = interfaceZ(top - 1);
    }
    let h = coupled ? cloudTopZ : surfaceDepth;
    depth[i] = zb + h;
    if (deckTop && deckTop[i] > 0) h = Math.max(h, deckTop[i] - zb);
    regime[i] = top >= 0 ? (coupled ? REGIME.COUPLED : REGIME.DECOUPLED) : buoyancy > 0 ? REGIME.SURFACE : REGIME.STABLE;
    mixingTop[i] = zb + Math.max(h, cloudTopZ);
    cloudTopCooling[i] = cooling;
    decoupling[i] = top >= 0 && !coupled ? zb + base0 : 0;
    const layerDepth = cloudTopZ - base0;
    let velocityCubed = 0;
    if (top >= 0 && layerDepth > 0) velocityCubed = g / thetaV[top * C + i] * cooling / (density(pi, top, i, theta) * cp) * layerDepth;
    const velocity = Math.cbrt(velocityCubed);
    radiativeVelocity[i] = velocity;
    let scale = friction[i];
    if (stability && buoyancy > 0 && h > 0) scale = friction[i] * Math.pow(1 - 15 * Math.max(-2, -0.1 * h * vonKarman * buoyancy / friction[i] ** 3), 0.25);
    for (let k = kTop; k < bottom; k++) {
      const z = interfaceZ(k);
      let diffusivity = 0;
      if (z < h) diffusivity += vonKarman * scale * z * (1 - z / h) ** 2;
      if (velocity > 0 && z > base0 && z < cloudTopZ) { const x = (z - base0) / layerDepth; diffusivity += profile * vonKarman * velocity * layerDepth * x * x * Math.sqrt(1 - x); }
      if (diffusivity > 0) mixing[k * C + i] = 0.5 * (density(pi, k, i, theta) + density(pi, k + 1, i, theta)) * diffusivity / (geopotential[k * C + i] / g - geopotential[(k + 1) * C + i] / g);
    }
    if (!entraining) return;
    const onset = shearOnset > 0 ? Math.min(1, buoyancy / shearOnset) : 1;
    const sheared = buoyancy > 0 ? onset * friction[i] ** 3 : 0;
    const surfaceInterface = (depthAbove) => { let kE = -1; for (let k = kTop; k < bottom; k++) if (interfaceZ(k) >= depthAbove) kE = k; return kE; };
    if (top >= 0) {
      const driven = coupled && buoyancy > 0;
      entrainmentVelocity[i] = entrain(i, pi, theta, q, qc, top - 1, lowest, coupled ? cloudTopZ : layerDepth, velocityCubed + (driven ? buoyancy * cloudTopZ : 0), driven ? sheared : 0);
      if (!coupled && buoyancy > 0 && surfaceDepth > 0) {
        const kE = surfaceInterface(surfaceDepth);
        if (kE >= lowest) entrain(i, pi, theta, q, qc, kE, bottom, surfaceDepth, buoyancy * surfaceDepth, sheared);
      }
    } else if (buoyancy > 0 && h > 0) {
      entrainmentVelocity[i] = entrain(i, pi, theta, q, qc, surfaceInterface(h), bottom, h, buoyancy * h, sheared);
    }
  }

  function solve(field, offset, stride, coefficient, coefficientStride, dt, columnMass) {
    let active = false;
    for (let k = kTop; k < bottom; k++) if (coefficient[k * coefficientStride] > 0) { active = true; break; }
    if (!active) return false;
    for (let j = 0; j < n; j++) {
      const k = kTop + j;
      mass[j] = columnMass * dSigma[k] / g;
      upper[j] = j > 0 ? dt * coefficient[(k - 1) * coefficientStride] / mass[j] : 0;
      lower[j] = k < bottom ? dt * coefficient[k * coefficientStride] / mass[j] : 0;
      rhs[j] = field[offset + k * stride];
    }
    let denominator = 1 + upper[0] + lower[0];
    gain[0] = -lower[0] / denominator;
    rhs[0] /= denominator;
    for (let j = 1; j < n; j++) {
      denominator = 1 + upper[j] + lower[j] + upper[j] * gain[j - 1];
      gain[j] = -lower[j] / denominator;
      rhs[j] = (rhs[j] + upper[j] * rhs[j - 1]) / denominator;
    }
    field[offset + (kTop + n - 1) * stride] = rhs[n - 1];
    for (let j = n - 2; j >= 0; j--) {
      rhs[j] -= gain[j] * rhs[j + 1];
      field[offset + (kTop + j) * stride] = rhs[j];
    }
    return true;
  }

  function mixColumn(i, pi, theta, q, qc, dt) {
    const coefficient = mixing.subarray(i);
    if (moistScheme && q && qc) {
      for (let k = kTop; k < K; k++) {
        const idx = k * C + i;
        thetaL[k] = theta[idx] - latentHeat * qc[idx] / (cp * exnerLayer[idx]);
        totalWater[k] = q[idx] + qc[idx];
      }
      if (!solve(thetaL, 0, 1, coefficient, C, dt, pi[i])) return;
      solve(totalWater, 0, 1, coefficient, C, dt, pi[i]);
      for (let k = kTop; k < K; k++) {
        if (!((k > kTop && coefficient[(k - 1) * C] > 0) || (k < bottom && coefficient[k * C] > 0))) continue;
        const idx = k * C + i;
        theta[idx] = thetaL[k]; q[idx] = totalWater[k]; qc[idx] = 0;
      }
      return;
    }
    if (!solve(theta, i, C, coefficient, C, dt, pi[i])) return;
    if (q) solve(q, i, C, coefficient, C, dt, pi[i]);
    if (qc) solve(qc, i, C, coefficient, C, dt, pi[i]);
  }

  const edgeCoefficient = new Float64Array(K), before = new Float64Array(K), share = new Float64Array(K);
  /*
   * The implicit mixing removes kinetic energy at each interface's shear
   * and in each layer's own increment; `dissipation` receives each
   * layer's loss of u² in those proportions, scaled so that the column's
   * mass-weighted loss is exact.
   */
  function mixEdges(pi, u, eFrom, eTo, dt, dissipation = null) {
    for (let e = eFrom; e < eTo; e++) {
      const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1], columnMass = 0.5 * (pi[a] + pi[b]);
      for (let k = kTop; k < K; k++) { edgeCoefficient[k] = 0.5 * (mixing[k * C + a] + mixing[k * C + b]); before[k] = u[k * E + e]; }
      if (!solve(u, e, E, edgeCoefficient, 1, dt, columnMass) || !dissipation) continue;
      let loss = 0, total = 0;
      for (let k = kTop; k < K; k++) {
        const m = columnMass * dSigma[k] / g, now = u[k * E + e], change = now - before[k];
        loss += m * (before[k] * before[k] - now * now);
        share[k] = m * change * change;
      }
      for (let k = kTop; k < bottom; k++) {
        const shear = u[k * E + e] - u[(k + 1) * E + e], part = dt * edgeCoefficient[k] * shear * shear;
        share[k] += part; share[k + 1] += part;
      }
      for (let k = kTop; k < K; k++) total += share[k];
      if (total <= 0) continue;
      for (let k = kTop; k < K; k++) dissipation[k * E + e] += loss * share[k] / (total * columnMass * dSigma[k] / g);
    }
  }

  return {
    diagnose, mixColumn, mixEdges, mixing, depth, buoyancyFlux, friction, entrainment: entrainmentVelocity, regime, mixingTop, cloudTopCooling, radiativeVelocity, decoupling, kTop, turbulence,
    shared: { mixing: mixingBuffer, depth: depthBuffer, buoyancyFlux: buoyancyBuffer, friction: frictionBuffer, entrainment: entrainmentBuffer, regime: regimeBuffer, mixingTop: mixingTopBuffer, cloudTopCooling: coolingBuffer, radiativeVelocity: velocityBuffer, decoupling: decouplingBuffer },
  };
}
