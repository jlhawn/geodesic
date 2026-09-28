import { LATENT_HEAT, EPSILON, saturationHumidity, saturationVaporPressure } from './moist.module.js';
import { CP_DRY, R_DRY, GRAVITY } from '../dynamics/sigmaCore.module.js';

/*
 * Bulk mixed-layer model of a stratocumulus-topped boundary layer
 * (Lilly 1968; Bretherton & Wyant 1997, BW97; Stevens 2002). The state
 * is the inversion height h and the liquid-water potential temperature
 * θ_l = θ − L q_l/(c_p Π) and total water q_t (specific) of a layer well
 * mixed from the surface to h. The layer is Boussinesq with the
 * reference density ρ of the forcing (by default the mean density of the
 * layer, (p_s − p_h)/(g h)):
 *
 *   dh/dt   = w_e + w_s(h)
 *   dθ_l/dt = [w_e Δθ_l + SH/(ρ c_p) − ΔF/(ρ c_p) + L P/(ρ c_p Π_b)]/h
 *   dq_t/dt = [w_e Δq_t + E/ρ − P/ρ]/h
 *
 * with Δ the jump from the layer to the free troposphere just above h,
 * w_s the (negative) subsidence, SH and E the surface sensible heat and
 * evaporation, ΔF = F(h) − F(0) the net longwave flux divergence across
 * the layer and P the drizzle leaving it. `step` advances h, h θ_l and
 * h q_t by forward Euler in that flux form, so the column's heat and
 * water change by exactly the fluxes of the step.
 *
 * Profile: below cloud base θ_v = θ_l (1 + δ q_t) is uniform, so the
 * Exner function falls linearly, Π(z) = Π_s − g z/(c_p θ_v), and cloud
 * base z_b is where q_s(θ_l Π, p) = q_t, found by Newton's method in Π.
 * Above it the pressure follows the hydrostatic equation (Heun steps on
 * `cloudLevels` intervals) with the saturation-adjusted θ_v, and the
 * liquid water path is the trapezoidal integral of ρ q_l over those
 * levels; the radiation sees the path below and above each level.
 *
 * Entrainment, closure 'radiative' (default): the Lilly form
 * w_e = A ΔF/(ρ c_p Δθ_v), Δθ_v the virtual jump between the free
 * troposphere at h and the cloud-top air; with no radiative cooling it
 * does not entrain. Closure 'buoyancy' is BW97's
 * w_e Δθ_v = 2.5 A ⟨w'θ_v'⟩, the layer-mean buoyancy flux, whose profile
 * is linear in w_e (below). Both carry the evaporative enhancement
 * of Nicholls & Turton (1986), A = a_1 [1 + a_2 χ* (1 − Δθ_vs/Δθ_v)]:
 * χ* is the fraction of free-tropospheric air that evaporates the
 * cloud-top liquid of a mixture and Δθ_vs the jump a saturated parcel
 * would feel (the saturated coefficients below applied to Δθ_l and
 * Δq_t), so χ* (1 − Δθ_vs/Δθ_v) is the mean evaporative cooling over all
 * mixtures in units of the jump. a_1 = entrainmentEfficiency = 0.2, the
 * efficiency of a dry convective layer; a_2 = evaporativeEnhancement =
 * 25, Caldwell & Bretherton's (2009) fit to DYCOMS-II (Nicholls &
 * Turton had 60).
 *
 * Limiter: rules of this kind, driven by the net forcing, become
 * singular under strong buoyancy reversal as Δθ_v → 0 (Stevens 2002,
 * his Eq. 22) — the enhancement grows as 1/Δθ_v and w_e as 1/Δθ_v², or
 * without bound once entrainment in the buoyancy closure supplies as
 * much buoyancy as it costs. A is therefore at most `maximumEfficiency`
 * = 1, the larger of the two efficiencies Stevens (2002) runs his
 * minimal model w_e = A ΔF/Δb with, and w_e at most `maximumEntrainment`
 * = 20 mm/s, several times the 3.8–5.9 mm/s of the RF01 simulations
 * (Stevens et al. 2005). Δθ_v enters both denominators no smaller than
 * `minimumJump` = 0.1 K, a numerical floor.
 *
 * Buoyancy flux: the θ_l and q_t fluxes are linear in z between the
 * surface fluxes and the entrainment fluxes −w_e Δ at h, less the
 * radiative (and drizzle) flux divergence below each level, and
 * w'θ_v' = A w'θ_l' + B w'q_t' with A = 1 + δ q_t, B = δ θ_l below
 * cloud base and the saturated coefficients A = c/(1 + γ),
 * B = c L/(c_p Π (1 + γ)) − θ, c = 1 + δ q_s − q_l + (1 + δ) T dq_s/dT,
 * γ = (L/c_p) dq_s/dT in the cloud.
 *
 * Decoupling: the buoyancy integral ratio of BW97, BIR = −(negative
 * part of the subcloud buoyancy-flux integral)/(the rest of the
 * integral). The layer is coupled below BW97's threshold 0.15
 * (`decouplingOnset`) and cloud cover is 1 there; it falls linearly to
 * the trade-cumulus cover `decoupledCover` (0.3) at Turton & Nicholls'
 * (1987) threshold 0.4 (`decoupledRatio`). Cover is 0 without cloud.
 * A layer with no virtual jump (Δθ_v ≤ 0) is uncapped: it neither
 * entrains nor carries cover.
 *
 * Shortwave: with forcing.solar, the shortwave flux (W/m²) incident at
 * the cloud top, the cloud absorbs S = solar × min(0.15, 0.4 W) for a
 * liquid water path W in kg/m² — about 4 % per 100 g/m², at most 15 %
 * (Stephens 1978; `shortwaveAbsorption`, `maximumShortwaveAbsorption`) —
 * spread through the cloud in proportion to its water. The layer then
 * sees the net forcing ΔF − S in its heat budget, as the heating
 * S/(ρ c_p h) in dθ_l/dt, and in the buoyancy-flux profile that the
 * buoyancy closure and the BIR integrate, so sunlight slows entrainment
 * under that closure and decouples the layer by day. The radiative
 * closure keeps the longwave ΔF: the longwave cooling that drives
 * entrainment sits in the cloud's top few tens of metres, where little
 * of the sunlight is absorbed. With forcing.absorbedSolar, a function of
 * W, S = absorbedSolar(W) instead; the coupled deck of
 * radiation.module.js passes the column's own absorption this way.
 * Without either the model is the nocturnal one.
 *
 * Drizzle (off by default): the cloud-base rate of Comstock et al.
 * (2004), 0.37 (LWP/N)^1.75 mm/day with LWP in g/m² and the droplet
 * number N in cm⁻³, leaves the layer at the surface; in the flux
 * profiles it forms uniformly in the cloud.
 *
 * Forcing: surfacePressure (Pa); optionally density; sensibleHeat
 * (W/m²) and evaporation (kg/m²/s), or seaSurfaceTemperature and
 * transferVelocity (C_T V, m/s) for bulk fluxes from the layer's surface
 * air; thetaLAbove and qtAbove, numbers or functions of height;
 * divergence D (w_s = −D z) or subsidence(z); radiation(below, above),
 * the net upward longwave flux (W/m²) at a level with liquid water paths
 * (kg/m²) below and above it; and optionally solar or absorbedSolar
 * (above).
 *
 * Carried height: a caller that keeps h from one step to the next (the
 * coupled deck of radiation.module.js) holds it with
 * `bound(h, floor, ceiling)` to [floor, min(ceiling, maximumHeight)]
 * (3000 m) — floor, the depth of the surface-driven boundary layer the
 * layer must at least span, winning over both caps — and, while the
 * layer does not run, lets it fall back with `relax(h, floor, dt)`,
 * h ← floor + (h − floor) e^(−dt/heightMemory) (1 day; 0 resets it to
 * floor). A height of 0 is unset and stays so.
 *
 * `dycomsLongwave` is the idealised longwave of the DYCOMS-II RF01 case
 * (Stevens et al. 2005): the net upward flux F0 e^(−κ W_above) +
 * F1 e^(−κ W_below) of a level with liquid paths W below and above it.
 * Its free-tropospheric term vanishes at the inversion and is left out.
 */
export const DYCOMS_LONGWAVE = { F0: 70, F1: 22, kappa: 85 };

export function dycomsLongwave(options = {}) {
  const { F0, F1, kappa } = { ...DYCOMS_LONGWAVE, ...options };
  return (below, above) => F0 * Math.exp(-kappa * above) + F1 * Math.exp(-kappa * below);
}

export const MIXED_LAYER_DEFAULTS = {
  closure: 'radiative', entrainmentEfficiency: 0.2, evaporativeEnhancement: 25,
  maximumEfficiency: 1, maximumEntrainment: 0.02, minimumJump: 0.1,
  decouplingOnset: 0.15, decoupledRatio: 0.4, decoupledCover: 0.3,
  drizzle: false, dropletNumber: 100, cloudLevels: 20, shortwaveAbsorption: 0.4, maximumShortwaveAbsorption: 0.15,
  heightMemory: 86400, maximumHeight: 3000,
};

export function createMixedLayer({ cp = CP_DRY, R = R_DRY, g = GRAVITY, latentHeat = LATENT_HEAT, referencePressure = 1e5, ...options } = {}) {
  const {
    closure, entrainmentEfficiency, evaporativeEnhancement, maximumEfficiency, maximumEntrainment, minimumJump,
    decouplingOnset, decoupledRatio, decoupledCover, drizzle, dropletNumber, cloudLevels, shortwaveAbsorption, maximumShortwaveAbsorption,
    heightMemory, maximumHeight,
  } = { ...MIXED_LAYER_DEFAULTS, ...options };
  if (closure !== 'radiative' && closure !== 'buoyancy') throw new Error(`closure must be 'radiative' or 'buoyancy', not ${closure}`);
  const kappa = R / cp, delta = 1 / EPSILON - 1, Lc = latentHeat / cp;
  const n = cloudLevels;
  const zNode = new Float64Array(n + 1), piNode = new Float64Array(n + 1), rhoNode = new Float64Array(n + 1), qlNode = new Float64Array(n + 1);
  const aNode = new Float64Array(n + 1), bNode = new Float64Array(n + 1), pathNode = new Float64Array(n + 1);
  const air = { T: 0, ql: 0, qv: 0, qs: 0, dqs: 0 };

  const pressure = (pi) => referencePressure * Math.pow(pi, 1 / kappa);
  const value = (f, z) => (typeof f === 'function' ? f(z) : f);

  function slope(T, p) {
    const e = saturationVaporPressure(T), dry = p - (1 - EPSILON) * e;
    return EPSILON * p * e * 17.67 * 243.5 / ((T - 29.65) * (T - 29.65) * dry * dry);
  }

  function saturate(thetaL, qt, pi, out = air) {
    const p = pressure(pi), Tl = thetaL * pi;
    let qs = saturationHumidity(Tl, p);
    if (qt <= qs) { out.T = Tl; out.ql = 0; out.qv = qt; out.qs = qs; out.dqs = slope(Tl, p); return out; }
    let T = Tl + Lc * (qt - qs) / (1 + Lc * slope(Tl, p));
    for (let it = 0; it < 30; it++) {
      qs = saturationHumidity(T, p);
      const step = (T - Tl - Lc * (qt - qs)) / (1 + Lc * slope(T, p));
      T -= step;
      if (Math.abs(step) < 1e-11 * T) break;
    }
    qs = saturationHumidity(T, p);
    out.T = T; out.ql = Math.max(0, qt - qs); out.qv = qt - out.ql; out.qs = qs; out.dqs = slope(T, p);
    return out;
  }

  const virtualTheta = (a, pi) => a.T / pi * (1 + delta * a.qv - a.ql);

  function cloudBaseExner(thetaL, qt, piSurface) {
    let pi = piSurface;
    for (let it = 0; it < 50; it++) {
      const T = thetaL * pi, p = pressure(pi), e = saturationVaporPressure(T), dry = p - (1 - EPSILON) * e;
      const qs = EPSILON * e / dry;
      const dlnT = 17.67 * 243.5 / ((T - 29.65) * (T - 29.65)) * p / dry;
      const derivative = dlnT * thetaL - p / (kappa * pi) / dry;
      const step = (Math.log(qs) - Math.log(qt)) / derivative;
      pi -= step;
      if (Math.abs(step) < 1e-14) break;
    }
    return pi;
  }

  function diagnose(state, forcing) {
    const { h, thetaL, qt } = state;
    const ps = forcing.surfacePressure, piS = Math.pow(ps / referencePressure, kappa);
    const thetaVDry = thetaL * (1 + delta * qt);
    let zb, piB;
    if (qt >= saturationHumidity(thetaL * piS, ps)) { zb = 0; piB = piS; }
    else { piB = cloudBaseExner(thetaL, qt, piS); zb = cp * thetaVDry * (piS - piB) / g; }
    const cloudy = zb < h;
    let lwp = 0, piH, thetaVTop, qlTop = 0, topA = 0, topB = 0, topGamma = 0, topSlope = 0;
    if (cloudy) {
      const dz = (h - zb) / n;
      let pi = piB;
      saturate(thetaL, qt, pi);
      let tv = virtualTheta(air, pi);
      for (let j = 0; j <= n; j++) {
        if (j > 0) {
          const predicted = pi - g * dz / (cp * tv);
          saturate(thetaL, qt, predicted);
          const tvPredicted = virtualTheta(air, predicted);
          pi -= 0.5 * g * dz / cp * (1 / tv + 1 / tvPredicted);
          saturate(thetaL, qt, pi);
          tv = virtualTheta(air, pi);
        }
        const p = pressure(pi), gamma = Lc * air.dqs, theta = air.T / pi;
        const c = 1 + delta * air.qv - air.ql + (1 + delta) * air.T * air.dqs;
        zNode[j] = zb + j * dz; piNode[j] = pi; qlNode[j] = air.ql;
        rhoNode[j] = p / (R * air.T * (1 + delta * air.qv - air.ql));
        aNode[j] = c / (1 + gamma);
        bNode[j] = c * Lc / (pi * (1 + gamma)) - theta;
        pathNode[j] = j === 0 ? 0 : pathNode[j - 1] + 0.5 * dz * (rhoNode[j - 1] * qlNode[j - 1] + rhoNode[j] * qlNode[j]);
        if (j === n) { topGamma = gamma; topSlope = air.dqs; }
      }
      lwp = pathNode[n];
      piH = piNode[n]; thetaVTop = tv; qlTop = qlNode[n]; topA = aNode[n]; topB = bNode[n];
    } else {
      piH = piS - g * h / (cp * thetaVDry);
      thetaVTop = thetaVDry;
    }
    const pH = pressure(piH);
    const density = forcing.density ?? (ps - pH) / (g * h);

    const thetaAbove = value(forcing.thetaLAbove, h), qtAbove = value(forcing.qtAbove, h);
    const jumpTheta = thetaAbove - thetaL, jumpQ = qtAbove - qt;
    const above = saturate(thetaAbove, qtAbove, piH, {});
    const jumpVirtual = virtualTheta(above, piH) - thetaVTop, capped = jumpVirtual > 0, jump = Math.max(minimumJump, jumpVirtual);

    let sensible, evaporation;
    if (forcing.seaSurfaceTemperature !== undefined) {
      const transfer = forcing.transferVelocity;
      sensible = density * cp * transfer * (forcing.seaSurfaceTemperature - thetaL * piS);
      evaporation = density * transfer * (saturationHumidity(forcing.seaSurfaceTemperature, ps) - qt);
    } else {
      sensible = forcing.sensibleHeat;
      evaporation = forcing.evaporation;
    }
    const subsidence = forcing.subsidence ? forcing.subsidence(h) : -forcing.divergence * h;
    const radiation = forcing.radiation;
    const fluxSurface = radiation(0, lwp), fluxTop = radiation(lwp, 0);
    const divergence = fluxTop - fluxSurface;
    const absorbed = forcing.absorbedSolar ? forcing.absorbedSolar(lwp) : (forcing.solar ?? 0) * Math.min(maximumShortwaveAbsorption, shortwaveAbsorption * lwp);
    const netDivergence = divergence - absorbed;
    const rain = drizzle && lwp > 0 ? 0.37 * Math.pow(1000 * lwp / dropletNumber, 1.75) / 86400 : 0;
    const drizzleHeat = Lc / piB;

    let chi = 0, efficiency = entrainmentEfficiency;
    if (cloudy && qlTop > 0 && capped) {
      const saturatedJump = topA * jumpTheta + topB * jumpQ;
      const demand = topSlope * piH * jumpTheta - jumpQ;
      chi = demand > 0 ? Math.min(1, qlTop * (1 + topGamma) / demand) : 1;
      efficiency = Math.min(maximumEfficiency, entrainmentEfficiency * (1 + evaporativeEnhancement * Math.max(0, chi * (1 - saturatedJump / jump))));
    }

    const heat0 = sensible / (density * cp), water0 = evaporation / density;
    const heatRate = (heat0 - netDivergence / (density * cp) + drizzleHeat * rain / density) / h;
    const waterRate = (water0 - rain / density) / h;
    const aDry = 1 + delta * qt, bDry = delta * thetaL, top = Math.min(zb, h);
    const dryFlux0 = (z) => aDry * (heat0 - z * heatRate) + bDry * (water0 - z * waterRate);
    const dryFlux1 = (z) => -(z / h) * (aDry * jumpTheta + bDry * jumpQ);
    let I0 = 0.5 * top * (dryFlux0(0) + dryFlux0(top)), I1 = 0.5 * top * (dryFlux1(0) + dryFlux1(top));
    const cloudFlux0 = (j) => {
      const z = zNode[j], precipitation = rain * (h - z) / (h - zb);
      const sunBelow = absorbed > 0 ? absorbed * pathNode[j] / lwp : 0;
      const heat = heat0 - z * heatRate - (radiation(pathNode[j], lwp - pathNode[j]) - fluxSurface - sunBelow) / (density * cp) - drizzleHeat * (precipitation - rain) / density;
      const water = water0 - z * waterRate + (precipitation - rain) / density;
      return aNode[j] * heat + bNode[j] * water;
    };
    const cloudFlux1 = (j) => -(zNode[j] / h) * (aNode[j] * jumpTheta + bNode[j] * jumpQ);
    if (cloudy) {
      const dz = (h - zb) / n;
      for (let j = 0; j < n; j++) {
        I0 += 0.5 * dz * (cloudFlux0(j) + cloudFlux0(j + 1));
        I1 += 0.5 * dz * (cloudFlux1(j) + cloudFlux1(j + 1));
      }
    }

    let entrainment = 0;
    if (capped) {
      if (closure === 'radiative') entrainment = Math.min(maximumEntrainment, Math.max(0, efficiency * divergence / (density * cp * jump)));
      else {
        const denominator = h * jump - 2.5 * efficiency * I1;
        entrainment = denominator > 0 ? Math.min(maximumEntrainment, Math.max(0, 2.5 * efficiency * I0 / denominator)) : maximumEntrainment;
      }
    }

    const integral = I0 + entrainment * I1;
    const dryLow = dryFlux0(0) + entrainment * dryFlux1(0), dryHigh = dryFlux0(top) + entrainment * dryFlux1(top);
    let negative = 0;
    if (dryLow < 0 && dryHigh < 0) negative = 0.5 * top * (dryLow + dryHigh);
    else if (dryLow < 0) negative = 0.5 * dryLow * top * dryLow / (dryLow - dryHigh);
    else if (dryHigh < 0) negative = 0.5 * dryHigh * top * dryHigh / (dryHigh - dryLow);
    const rest = integral - negative;
    const ratio = negative < 0 ? (rest > 0 ? -negative / rest : Infinity) : 0;
    const decoupled = ratio <= decouplingOnset ? 1 : ratio >= decoupledRatio ? decoupledCover : 1 - (1 - decoupledCover) * (ratio - decouplingOnset) / (decoupledRatio - decouplingOnset);
    const cover = cloudy && lwp > 0 && capped ? decoupled : 0;

    const thetaSource = entrainment * thetaAbove + subsidence * thetaL + heat0 - netDivergence / (density * cp) + drizzleHeat * rain / density;
    const waterSource = entrainment * qtAbove + subsidence * qt + water0 - rain / density;
    return {
      h, thetaL, qt, cloudBase: Math.min(zb, h), liquidWaterPath: lwp, topLiquid: qlTop, cover, cloudy,
      density, surfacePressure: ps, topPressure: pH, cloudBaseExner: piB,
      entrainment, subsidence, efficiency, mixingFraction: chi,
      thetaLAbove: thetaAbove, qtAbove, thetaLJump: jumpTheta, qtJump: jumpQ, virtualJump: jumpVirtual,
      sensibleHeat: sensible, evaporation, drizzle: rain, drizzleHeating: drizzleHeat * rain,
      radiativeDivergence: divergence, absorbedShortwave: absorbed, buoyancyIntegral: integral, buoyancyIntegralRatio: ratio,
      convectiveVelocity: Math.cbrt(Math.max(0, 2.5 * g / thetaVDry * integral)),
      heatSource: thetaSource, waterSource,
    };
  }

  function step(state, forcing, dt, d = diagnose(state, forcing)) {
    const h = state.h + dt * (d.entrainment + d.subsidence);
    const heat = state.h * state.thetaL + dt * d.heatSource;
    const water = state.h * state.qt + dt * d.waterSource;
    return { h, thetaL: heat / h, qt: water / h };
  }

  function bound(h, floor, ceiling = maximumHeight) {
    return Math.max(floor, Math.min(maximumHeight, ceiling, h));
  }

  function relax(h, floor, dt) {
    if (!(h > 0)) return h;
    return heightMemory > 0 ? floor + (h - floor) * Math.exp(-dt / heightMemory) : floor;
  }

  return { diagnose, step, bound, relax, maximumHeight };
}
