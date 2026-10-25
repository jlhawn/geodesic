import { LAYERED_TABLE } from './layeredTable.module.js';
import { NEAR_INFRARED_TABLE } from './nearInfraredTable.module.js';
import { OZONE_SHARES, OZONE_COEFFICIENTS, VISIBLE_VAPOR } from './shortwaveGases.module.js';

/*
 * The shortwave layer by layer (radiation option shortwave 'layered'): in each of
 * 32 near-infrared g-points (layeredTable.module.js: the k-distributions of
 * nearInfraredTable.module.js folded by scripts/layeredFit.mjs) and two visible
 * slices, every layer of the column has a depth that is its gases' (the g-point's
 * absorber in the near infrared; ozone and water vapor in the visible), its
 * Rayleigh scattering, its aerosol's (visible) and its cloud's, with the single
 * scattering albedo and asymmetry of the mixture. The layer's reflection and
 * transmission are the delta-Eddington two-stream's (Joseph, Wiscombe and Weinman
 * 1976, in Meador and Weaver's 1980 solution, its Eddington coefficient gamma2 held
 * at or above 0, which the Eddington closure takes below 0 where the albedo is below
 * 1/(4 - 3g)) for the beam at the cosine 1/m of the airmass m = 35 / sqrt(1224 mu^2 + 1)
 * that the gases' paths take, and for diffuse light. The layers are added from the
 * surface up (the albedo and source of everything below each interface) and the
 * fluxes found from the top down, so that a layer's heating is the divergence of the
 * net flux and the surface and the top take the rest: absorbed plus reflected is the
 * beam to rounding, slice by slice.
 *
 * The cloudy columns are one region, of cover c, holding each layer's cloud water
 * over c with its depth times each of the model's lognormal nodes x; the clear region
 * is the rest. Cloud optics: the model's visible extinction per kg of condensate
 * (cloudOptics: 1500 / r_e for liquid, 1000 (3.448e-3 + 2.431 / r_e) for ice) times
 * the liquid extinction of the slice's band over the first band's in Slingo (1989),
 * and the single scattering albedo and asymmetry of Slingo's liquid and Ebert and
 * Curry's (1992) ice in the band of their four that holds the slice's centre
 * wavelength.
 *
 * The visible slices hold the beam that is neither the near infrared's nor the
 * ultraviolet's: CLIRAD-SW's ozone bands 1-6 (ULTRAVIOLET_BANDS) take their share before any
 * cloud, in the column's own ozone as ozoneAbsorptivity finds it, and are left to the
 * clear column they are added to; what remains weighs the visible slices, whose ozone
 * (bands 7 and 8) and vapor coefficients are raised to take what the sequential scheme
 * takes per unit path and whose Rayleigh depths are VISIBLE_RAYLEIGH_SCALE times its
 * own, the mean over their band (0.279-0.607 um) being that much deeper than over
 * its 0.297-0.683 um (scripts/rayleighReference.mjs, OZONE=0.0208 VISIBLE=0.3844).
 * The clear region then follows the sequential scheme's clear column to 0.3 % of the
 * beam at zenith angles to 70 degrees, so that its difference from the cloudy region is
 * the cloud's.
 */
export const ULTRAVIOLET_BANDS = 6, VISIBLE_RAYLEIGH_SCALE = 1.329;
export const VISIBLE_OZONE_ABSORPTION = OZONE_SHARES[6] * OZONE_COEFFICIENTS[6] + OZONE_SHARES[7] * OZONE_COEFFICIENTS[7];
const expRatio = (x) => (x < 1e-9 ? 1 : -Math.expm1(-x) / x);
export const SLINGO_LIQUID = {
  extinction: [[2.817e-02, 1.305], [2.682e-02, 1.346], [2.264e-02, 1.454], [1.281e-02, 1.641]],
  absorption: [[-5.62e-08, 1.63e-07], [-6.94e-06, 2.35e-05], [4.64e-04, 1.24e-03], [0.201, 7.56e-03]],
  asymmetry: [[0.829, 2.482e-03], [0.794, 4.226e-03], [0.754, 6.560e-03], [0.826, 4.353e-03]],
};
export const EBERT_CURRY_ICE = {
  absorption: [[1.00e-05, 0], [1.10e-04, 1.405e-05], [1.861e-02, 8.328e-04], [0.46658, 2.05e-05]],
  asymmetry: [[0.7661, 5.851e-04], [0.7730, 5.665e-04], [0.794, 7.267e-04], [0.9595, 1.076e-04]],
};
export const NEAR_INFRARED_POINTS = LAYERED_TABLE.bands.reduce((s, b) => s + b.points.length, 0);
export const NEAR_INFRARED_SHARE = LAYERED_TABLE.bands.reduce((s, b) => s + b.share, 0);

/*
 * The layer of extinction depth tau, co-albedo c = 1 - omega and asymmetry g for a
 * beam of cosine mu0: out.reflect the share of the beam it reflects, out.scatter the
 * share it transmits scattered, out.direct the share it transmits unscattered (the
 * forward peak, delta-scaled, counts as direct), and out.reflectDiffuse and
 * out.transmitDiffuse the reflection and total transmission of diffuse light.
 * The solution is written in the delta-scaled depth tau', the diffuse attenuation k and
 * Q = tau' (1 - exp(-2 k tau')) / (2 k tau'), with the factor 1 - k mu0 divided out of
 * the beam's terms: it is regular where the layer does not absorb (k = 0) and where the
 * beam's cosine meets the attenuation (k mu0 = 1), and accurate in single precision,
 * where the exponential form of the same solution cancels terms of order 1 in a thin
 * layer that absorbs little and errs by 1e-4 of the beam.
 */
export function twoStreamLayer(tau, c, g, mu0, out) {
  const f = g * g, q = 1 - f + f * c, scaled = q * tau, cs = c / q, gs = g / (1 + g), w = 1 - cs;
  const g1 = (3 - 3 * gs + cs * (4 + 3 * gs)) / 4, g2 = Math.max(0, (3 - 3 * gs - cs * (4 - 3 * gs)) / 4), k = Math.sqrt(g2 > 0 ? 2 * cs * (g1 + g2) : g1 * g1);
  const u = k * scaled, E = Math.exp(-u), Q = scaled * expRatio(2 * u), den = 1 + (g1 - k) * Q;
  out.reflectDiffuse = g2 * Q / den;
  out.transmitDiffuse = E / den;
  if (w < 1e-6) {
    out.reflect = 0; out.scatter = 0; out.direct = Math.exp(-scaled / mu0);
    return;
  }
  const g3 = (2 - 3 * gs * mu0) / 4, g4 = 1 - g3, a1 = g1 * g4 + g2 * g3, a2 = g1 * g3 + g2 * g4, km = k * mu0;
  const v = scaled / mu0, t = Math.exp(-v), d = scaled * (1 - km) / mu0, ratio = expRatio(Math.abs(d)), F = d >= 0 ? E * E : E * t, G = d >= 0 ? E : t, pref = w / ((1 + km) * den);
  out.reflect = pref * ((a2 + k * g3) * Q - (a2 * mu0 - g3) * v * F * ratio);
  out.scatter = pref * ((g4 + a1 * mu0) * v * G * ratio - t * (a1 - k * g4) * Q);
  out.direct = t;
}

/*
 * Layers 0..K-1 from the top (the arrays reflect, scatter, direct, reflectDiffuse and
 * transmitDiffuse of twoStreamLayer) added over a surface of albedo surfaceAlbedo for the beam
 * and diffuseAlbedo for diffuse light, under a unit beam: adds `weight` times, to
 * out.heat[j], the share of the beam that layer j absorbs and, to out.surface,
 * out.reflected, out.down and out.direct, the share that the surface absorbs, the top
 * reflects and the surface receives in all and as beam. `albedo` and `source` (K + 1
 * long) hold the albedo and the beam-sourced upward flux over each interface.
 */
export function addLayers(layers, K, surfaceAlbedo, diffuseAlbedo, weight, out, albedo, source) {
  const { reflect, scatter, direct, reflectDiffuse, transmitDiffuse } = layers;
  albedo[K] = diffuseAlbedo; source[K] = surfaceAlbedo;
  for (let j = K - 1; j >= 0; j--) {
    const rest = 1 / (1 - reflectDiffuse[j] * albedo[j + 1]);
    albedo[j] = reflectDiffuse[j] + transmitDiffuse[j] * transmitDiffuse[j] * albedo[j + 1] * rest;
    source[j] = reflect[j] + transmitDiffuse[j] * (albedo[j + 1] * scatter[j] + source[j + 1] * direct[j]) * rest;
  }
  let beam = 1, down = 0, net = 1 - source[0];
  out.reflected += weight * source[0];
  for (let j = 0; j < K; j++) {
    const next = beam * direct[j], rest = 1 / (1 - reflectDiffuse[j] * albedo[j + 1]);
    const below = (scatter[j] * beam + transmitDiffuse[j] * down + reflectDiffuse[j] * source[j + 1] * next) * rest;
    const through = next + below - (albedo[j + 1] * below + source[j + 1] * next);
    out.heat[j] += weight * (net - through);
    net = through; beam = next; down = below;
  }
  out.surface += weight * net;
  out.down += weight * (beam + down);
  out.direct += weight * beam;
}

/*
 * The layered shortwave of K layers. solve(inputs) fills `clear` and `cloudy`, each
 * { heat[k], surface, reflected, down, direct }, shares of the beam summed over the
 * slices of the near infrared and the visible. inputs:
 *   mu, airmass: the sun's cosine and the beam's magnification
 *   air, mixing: the near infrared's layers and gases (nearInfraredBands')
 *   dSigma, pressureScale: the layers' share of the mass and the surface pressure over 1013.25 hPa
 *   ozone[k] (cm-atm), vaporPath[k] (g/cm2 scaled, 0.1 vaporScaling q mass), aerosol, aerosolFraction[k]
 *   surfaceAlbedo, diffuseAlbedo
 *   cloud: null or { cover, count, layers[], water[k] (kg/m2 of the cell), liquid[k], reIce[k], reLiquid, nodes }
 */
export function createLayered(K, { rayleighBands, nearInfraredScale = 1, visibleVaporStrength, visibleVaporCoefficient, aerosolAlbedo, aerosolAsymmetry }) {
  const field = () => new Float64Array(K);
  const makeLayers = () => ({ reflect: field(), scatter: field(), direct: field(), reflectDiffuse: field(), transmitDiffuse: field() });
  const clearLayers = makeLayers(), cloudLayers = makeLayers();
  const albedo = new Float64Array(K + 1), source = new Float64Array(K + 1);
  const depth = field(), aerosolDepth = field(), one = {};
  const clear = { heat: field(), surface: 0, reflected: 0, down: 0, direct: 0 }, cloudy = { heat: field(), surface: 0, reflected: 0, down: 0, direct: 0 };
  const base = [Array.from({ length: 4 }, field), Array.from({ length: 4 }, field)], self = field();
  const sourceBands = new Map(NEAR_INFRARED_TABLE.bands.map((b) => [b.band, b]));
  const reset = (out) => { out.heat.fill(0); out.surface = 0; out.reflected = 0; out.down = 0; out.direct = 0; };

  function setLayer(layers, k, tau, c, g, mu) {
    twoStreamLayer(tau, c, g, mu, one);
    layers.reflect[k] = one.reflect; layers.scatter[k] = one.scatter; layers.direct[k] = one.direct; layers.reflectDiffuse[k] = one.reflectDiffuse; layers.transmitDiffuse[k] = one.transmitDiffuse;
  }

  function clearLayer(k, gas, rayleigh, aerosol, mu) {
    const tau = gas + rayleigh + aerosol, scattered = rayleigh + aerosolAlbedo * aerosol;
    setLayer(clearLayers, k, tau, (gas + (1 - aerosolAlbedo) * aerosol) / tau, scattered > 0 ? aerosolAsymmetry * aerosolAlbedo * aerosol / scattered : 0, mu);
  }

  function cloudLayer(k, band, water, cloud, gas, rayleigh, aerosol, mu) {
    const liquidShare = cloud.liquid[k], reLiquid = cloud.reLiquid, reIce = cloud.reIce[k];
    const [a0, b0] = SLINGO_LIQUID.extinction[0], [a, b] = SLINGO_LIQUID.extinction[band];
    const liquidTau = liquidShare * water * 1500 / reLiquid * (a + b / reLiquid) / (a0 + b0 / reLiquid), iceTau = (1 - liquidShare) * water * 1000 * (3.448e-3 + 2.431 / reIce);
    const liquidCo = Math.max(0, SLINGO_LIQUID.absorption[band][0] + SLINGO_LIQUID.absorption[band][1] * reLiquid), iceCo = EBERT_CURRY_ICE.absorption[band][0] + EBERT_CURRY_ICE.absorption[band][1] * reIce;
    const liquidG = SLINGO_LIQUID.asymmetry[band][0] + SLINGO_LIQUID.asymmetry[band][1] * reLiquid, iceG = EBERT_CURRY_ICE.asymmetry[band][0] + EBERT_CURRY_ICE.asymmetry[band][1] * reIce;
    const cloudTau = liquidTau + iceTau, cloudScattered = liquidTau * (1 - liquidCo) + iceTau * (1 - iceCo);
    const tau = gas + rayleigh + aerosol + cloudTau, absorbed = gas + (1 - aerosolAlbedo) * aerosol + liquidTau * liquidCo + iceTau * iceCo;
    const scattered = rayleigh + aerosolAlbedo * aerosol + cloudScattered;
    setLayer(cloudLayers, k, tau, absorbed / tau, scattered > 0 ? (aerosolAsymmetry * aerosolAlbedo * aerosol + liquidTau * (1 - liquidCo) * liquidG + iceTau * (1 - iceCo) * iceG) / scattered : 0, mu);
  }

  function slice(inputs, gas, rayleigh, visible, band, weight) {
    const { surfaceAlbedo, diffuseAlbedo, dSigma, pressureScale, cloud } = inputs, mu = 1 / inputs.airmass;
    for (let k = 0; k < K; k++) clearLayer(k, gas[k], rayleigh * pressureScale * dSigma[k], visible ? aerosolDepth[k] : 0, mu);
    addLayers(clearLayers, K, surfaceAlbedo, diffuseAlbedo, weight, clear, albedo, source);
    if (!cloud || cloud.count === 0) return;
    for (const key in clearLayers) cloudLayers[key].set(clearLayers[key]);
    for (const [x, nodeWeight] of cloud.nodes) {
      for (let n = 0; n < cloud.count; n++) {
        const k = cloud.layers[n];
        cloudLayer(k, band, x * cloud.water[k] / cloud.cover, cloud, gas[k], rayleigh * pressureScale * dSigma[k], visible ? aerosolDepth[k] : 0, mu);
      }
      addLayers(cloudLayers, K, surfaceAlbedo, diffuseAlbedo, weight * nodeWeight, cloudy, albedo, source);
    }
  }

  function solve(inputs) {
    reset(clear); reset(cloudy);
    const { air, mixing } = inputs;
    for (let k = 0; k < K; k++) aerosolDepth[k] = inputs.aerosol * inputs.aerosolFraction[k];
    for (const band of LAYERED_TABLE.bands) {
      const data = sourceBands.get(band.band);
      self.fill(0);
      band.gases.forEach((name, i) => {
        const par = data.gases[name], isVapor = name === 'h2o';
        for (let k = 0; k < K; k++) {
          const u = isVapor ? air.vapour[k] : mixing[name] * air.dry[k], lnV = Math.log(Math.sqrt(air.p[k] * air.p[k] + par.pD * par.pD) / LAYERED_TABLE.pRef), dT = air.T[k] - 250;
          for (let j = 0; j < par.w.length; j++) base[i][j][k] = par.k[j] * u * Math.exp(par.n[j] * lnV + par.t[j] * dT);
          if (isVapor && par.self > 0) self[k] = par.self * u * air.e[k] / LAYERED_TABLE.pRef * Math.exp(par.selfT * Math.log(296 / air.T[k]));
        }
      });
      for (const point of band.points) {
        for (let k = 0; k < K; k++) {
          let d = self[k];
          for (let i = 0; i < band.gases.length; i++) for (let j = 0; j < point.mix[i].length; j++) d += point.mix[i][j] * base[i][j][k];
          depth[k] = point.scale * d;
        }
        slice(inputs, depth, band.rayleigh * nearInfraredScale, false, band.slingo, band.share * point.weight);
      }
    }
    let ozonePath = 0;
    for (let k = 0; k < K; k++) ozonePath += inputs.ozone[k];
    ozonePath *= inputs.airmass;
    let ultraviolet = 0;
    for (let b = 0; b < ULTRAVIOLET_BANDS; b++) ultraviolet -= OZONE_SHARES[b] * Math.expm1(-OZONE_COEFFICIENTS[b] * ozonePath);
    const visibleShare = 1 - NEAR_INFRARED_SHARE - ultraviolet, ozoneStrength = VISIBLE_OZONE_ABSORPTION / visibleShare, vaporStrength = VISIBLE_VAPOR.share * visibleVaporCoefficient * visibleVaporStrength / visibleShare;
    for (const [bandWeight, rayleigh] of rayleighBands) {
      for (let k = 0; k < K; k++) depth[k] = ozoneStrength * inputs.ozone[k] + vaporStrength * inputs.vaporPath[k];
      slice(inputs, depth, VISIBLE_RAYLEIGH_SCALE * rayleigh, true, 0, visibleShare * bandWeight);
    }
  }
  return { solve, clear, cloudy };
}
