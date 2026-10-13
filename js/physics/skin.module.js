import { JERLOV } from '../ocean/layered.module.js';

/*
 * The sea's skin over open water: a diurnal warm layer and a cool skin,
 * a profile inside the slab's heat that stores none of its own.
 *
 * The warm layer is Zeng and Beljaars (2005) with Takaya et al. (2010)'s
 * stable φ and Langmuir factor, as NEMO's diu_bulk codes it: over the top
 * d = min(depth, h₀) the profile is T_f + ΔT_w[1 − (z/d)^ν], and
 *   ΔT_w ← clamp((ΔT_w + Δt Q_w (ν+1)/(ν d ρc_p)) / (1 + Δt (ν+1) κ u*_w f_La/(d φ(ζ))), 0, cap)
 * with Q_w = Q_ns + f_d R_s, f_d the Jerlov type's light absorbed above d,
 * ζ = d κ g α Q_w/(ρc_p u*_w³), φ = 1 + (5ζ + 4ζ²)/(1 + 3ζ + ζ²/4) when
 * ζ ≥ 0 and (1 − 16ζ)^(−1/2) when not, and α COARE's
 * max(1e-5, 2.1e-5 (T − 273.15 + 3.2)^0.79).
 * The cool skin is Fairall et al. (1996a) as COARE 3.6 codes it, two
 * passes from δ = 1 mm: Q_c = −(Q_ns + f_s(δ) R_s), λ = 6 when α Q_c ≤ 0
 * and else 6 [1 + (B_c α Q_c/u*_a⁴)^(3/4)]^(−1/3), δ = min(0.01, λ ν_w/u*_w),
 * ΔT_c = clamp(Q_c δ/k_w, 0, 1).
 * The atmosphere sees T_slab + s ΔT_w − ΔT_c, s = 1 − ν d/((ν+1) h₀), so
 * that the profile's mean over h₀ is the slab's temperature.
 * GPU twin: SKIN_WGSL in gpu/physics.gpu.js.
 */
export const SKIN_DEFAULTS = { warmLayer: false, coolSkin: false, warmLayerDepth: 3, warmLayerProfile: 0.3, warmLayerLangmuir: 2.23, warmLayerCap: 5 };
export const SKIN = { karman: 0.4, gravity: 9.81, waterDensity: 1025, specificHeat: 3985, conductivity: 0.6, viscosity: 1e-6, coolCap: 1, firstDepth: 1e-3, maxDepth: 0.01, ustarFloor: 1e-6 };
export const SEA_RHO_CP = SKIN.waterDensity * SKIN.specificHeat;

export function skinOptions(surface = {}) {
  const o = Object.fromEntries(Object.keys(SKIN_DEFAULTS).map((key) => [key, surface[key] ?? SKIN_DEFAULTS[key]]));
  if (!(o.warmLayerDepth > 0 && o.warmLayerProfile > 0 && o.warmLayerLangmuir >= 1 && o.warmLayerCap >= 0)) throw new Error('the warm layer takes a depth and profile above 0, a Langmuir factor of 1 or more and a cap of 0 or more');
  return o;
}

export function skinRequiresRoughness(o, mode) {
  if ((o.warmLayer || o.coolSkin) && mode !== 'roughness') throw new Error('warmLayer and coolSkin need exchange \'roughness\'');
}

export const expansion = (T) => Math.max(1e-5, 2.1e-5 * Math.pow(Math.max(0.01, T - 273.15 + 3.2), 0.79));

export const absorbedAbove = ([R, z1, z2], d) => 1 - (R * Math.exp(-d / z1) + (1 - R) * Math.exp(-d / z2));

export function stability(zeta) {
  return zeta >= 0 ? 1 + (5 * zeta + 4 * zeta * zeta) / (1 + 3 * zeta + 0.25 * zeta * zeta) : 1 / Math.sqrt(1 - 16 * zeta);
}

export const waterFriction = (ustarAir, airDensity) => Math.max(SKIN.ustarFloor, ustarAir * Math.sqrt(airDensity / SKIN.waterDensity));

/*
 * One step of ΔT_w over a slab of capacity `capacity` (J/m²/K) at T_slab.
 */
export function warmLayerStep(warm, nonSolar, solar, ustarAir, airDensity, slabT, capacity, dt, o, light = JERLOV.IB) {
  const { karman, gravity } = SKIN, nu = o.warmLayerProfile;
  const d = Math.min(o.warmLayerDepth, capacity / SEA_RHO_CP);
  const flux = nonSolar + absorbedAbove(light, d) * solar, usw = waterFriction(ustarAir, airDensity);
  const zeta = d * karman * gravity * expansion(slabT) * flux / (SEA_RHO_CP * usw * usw * usw);
  const next = (warm + dt * flux * (nu + 1) / (nu * d * SEA_RHO_CP)) / (1 + dt * (nu + 1) * karman * usw * o.warmLayerLangmuir / (d * stability(zeta)));
  return Math.min(o.warmLayerCap, Math.max(0, next));
}

export function coolSkinDepression(nonSolar, solar, ustarAir, airDensity, slabT) {
  const { gravity, waterDensity, specificHeat, conductivity, viscosity, coolCap, firstDepth, maxDepth } = SKIN;
  const alpha = expansion(slabT), usw = waterFriction(ustarAir, airDensity), ua = Math.max(SKIN.ustarFloor, ustarAir);
  const bigC = 16 * gravity * specificHeat * (waterDensity * viscosity) ** 3 / (conductivity * conductivity * airDensity * airDensity);
  let delta = firstDepth, loss = 0;
  for (let pass = 0; pass < 2; pass++) {
    const share = 0.065 + 11 * delta - 6.6e-5 / delta * (1 - Math.exp(-delta / 8e-4));
    loss = -(nonSolar + share * solar);
    const lambda = alpha * loss > 0 ? 6 / Math.cbrt(1 + Math.pow(bigC * alpha * loss / (ua * ua * ua * ua), 0.75)) : 6;
    delta = Math.min(maxDepth, lambda * viscosity / usw);
  }
  return Math.min(coolCap, Math.max(0, loss * delta / conductivity));
}

export const surfaceShare = (capacity, o) => {
  const h0 = capacity / SEA_RHO_CP, d = Math.min(o.warmLayerDepth, h0);
  return 1 - o.warmLayerProfile * d / ((o.warmLayerProfile + 1) * h0);
};

/*
 * The per-cell skin: `warm` (ΔT_w) and `cool` (ΔT_c) in K, on shared
 * buffers; offset(i, capacity) is what the atmosphere adds to T_slab,
 * step(...) advances cell i after the slab's own update, clear(i) where
 * the cell has ice.
 */
export function createSkin(C, surface = {}, { light = JERLOV.IB, buffers = null } = {}) {
  const o = skinOptions(surface), on = o.warmLayer || o.coolSkin;
  const array = (name) => new Float64Array(buffers && buffers[name] ? buffers[name] : new SharedArrayBuffer(8 * C));
  const warm = array('warm'), cool = array('cool');
  return {
    on, options: o, warm, cool, light,
    offset(i, capacity) { return on ? (o.warmLayer ? surfaceShare(capacity, o) * warm[i] : 0) - cool[i] : 0; },
    step(i, nonSolar, solar, ustarAir, airDensity, slabT, capacity, dt) {
      if (o.warmLayer) warm[i] = warmLayerStep(warm[i], nonSolar, solar, ustarAir, airDensity, slabT, capacity, dt, o, light);
      if (o.coolSkin) cool[i] = coolSkinDepression(nonSolar, solar, ustarAir, airDensity, slabT);
    },
    clear(i) { warm[i] = 0; cool[i] = 0; },
    load(saved) {
      warm.fill(0); cool.fill(0);
      if (saved && saved.warmLayer && saved.warmLayer.length === C) warm.set(saved.warmLayer);
      if (saved && saved.coolSkin && saved.coolSkin.length === C) cool.set(saved.coolSkin);
    },
    serialize() { return on ? { warmLayer: Float64Array.from(warm), coolSkin: Float64Array.from(cool) } : {}; },
    shared: { warm: warm.buffer, cool: cool.buffer },
  };
}
