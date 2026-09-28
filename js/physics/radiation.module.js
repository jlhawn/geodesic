import { LATENT_HEAT, EPSILON, saturationHumidity, liftingCondensationLevel } from './moist.module.js';
import { createMixedLayer, dycomsLongwave } from './mixedLayer.module.js';
export const STEFAN_BOLTZMANN = 5.670374419e-8;
export const SOLAR_CONSTANT = 1362;
export const AXIAL_TILT = 23.44 * Math.PI / 180;
export const DAY = 86400;
export const YEAR = 365 * DAY;

/*
 * Unit vector toward the sun at model time t. t = 0 is the spring equinox
 * with the sun over the +x meridian; the subsolar latitude follows
 * AXIAL_TILT·sin(2πt/YEAR) and the sun circles westward once per day.
 */
export function sunDirection(t, out = new Float64Array(3)) {
  const tilt = -AXIAL_TILT * Math.sin((t % YEAR) * 2 * Math.PI / YEAR);
  const x = Math.cos(tilt), z = -Math.sin(tilt);
  const spin = -((t % DAY) * 2 * Math.PI / DAY);
  out[0] = x * Math.cos(spin);
  out[1] = x * Math.sin(spin);
  out[2] = z;
  return out;
}

/*
 * Three-band gray longwave column with a slab-ocean surface and a bulk
 * sensible heat flux. A window band carrying the fraction `window` of
 * blackbody emission is transparent: the surface radiates it straight
 * to space. A vapour band whose optical depth follows Frierson et al.
 * (2006) — tau0(lat) = tauEquator + (tauPole − tauEquator) sin²lat,
 * distributed in the vertical as tau0 (f·σ + (1 − f)·σ⁴) — concentrates
 * near the surface like water vapour. A well-mixed-gas band carrying
 * `gasFraction` of the emission has the optical depth gasOpticalDepth
 * spread uniformly per unit mass, so thin high layers keep an
 * emissivity they can cool with, as CO₂'s 15 µm band lets the
 * stratosphere do. In each absorbing band every emission is either
 * absorbed on its way or leaves through the top or reaches the surface,
 * so the layer and surface energy fluxes sum exactly to absorbed solar
 * minus outgoing longwave. With vaporCoupling > 0 (m²/kg) and a
 * humidity field, the vapour band's optical depth is instead
 * vaporCoupling times each layer's water mass, so the greenhouse
 * follows the model's own humidity.
 *
 * Clouds: each layer's cloud water path gives it a gray emissivity
 * 1 − exp(−cloudAbsorption × path) that joins every longwave band —
 * including the window, which is transparent only where there is no
 * cloud. In the shortwave the column's cloud optical depth
 * cloudScattering × path reflects the beam by the two-stream
 * reflectance τ / (τ + 2μ). What reaches the surface is direct beam,
 * exp(−τ/μ) of it less the clear-sky `skylight` fraction, and diffuse
 * light, the rest; the surface reflects each with its own albedo, and
 * the multiple reflections between surface and cloud base (diffuse,
 * at the mean cosine DIFFUSE_MU) are summed. The two albedos are given
 * per cell (open water or sea ice); `albedo` is the default for both.
 *
 * Marine stratocumulus: over the part of a cell that is ice-free sea
 * (`openSea`, the per-cell fraction the caller passes; no deck without it)
 * a diagnostic deck covers the fraction f of the column. f is
 * 0.19 + 0.08 (EIS − 1) clamped to [0, 1] — 0.2 at the warm pool's EIS
 * of about 1 K, 0.67 at the south-east Pacific deck's 7 K, the 6–8 %
 * per K of Wood & Bretherton (2006) — times a ramp from 0 at a 5 °C
 * surface to 1 at 10 °C that keeps the deck off polar seas, times
 * `openSea`. EIS is their estimated inversion strength
 * LTS − Γ_θ (z_700 − z_LCL): the lower-tropospheric stability LTS is the
 * potential temperature of the layer nearest σ = 0.7 less the lowest
 * layer's, z_700 that layer's height above the lowest layer, z_LCL the
 * height of the lifting condensation level of the lowest layer's air
 * (Bolton's, as the moist physics finds it, reached along the dry
 * adiabat), and Γ_θ = g/c_p − Γ_m the potential-temperature gradient of
 * the moist adiabat at 850 hPa and the mean temperature of the two
 * layers, so EIS counts only the θ at σ = 0.7 beyond what a moist
 * adiabat from cloud base reaches. With stratusIndex 'ectei' the fit
 * takes instead the estimated cloud-top entrainment index of Kawai,
 * Koshiro & Webb (2017), ECTEI = EIS − 0.23 (L/c_p)(q_lowest − q_700)
 * with q_700 the humidity of the σ = 0.7 layer, which lowers the cover
 * where the air the deck entrains is dry. The deck fills the boundary
 * layer from z_LCL to the boundary-layer top, `mixedDepth` metres above
 * the lowest layer, and its water path is stratusScale × ½ Γ_l Δz² for
 * that thickness Δz, at most stratusWaterMax, with Γ_l the adiabatic
 * liquid-water lapse rate at cloud base. The water sits in the layer
 * nearest σ = stratusSigma and never enters qc. The deck and the clear
 * part of the column are two independent columns: every shortwave
 * quantity is the f-weighted mean of the column with and without the
 * deck's water, and the deck layer's cloud emissivity is the f-weighted
 * mean of its emissivity with and without it. `stratus: false` removes
 * the deck.
 *
 * With mixedLayerDeck the deck's cover and water path come instead from
 * the mixed-layer model (mixedLayer.module.js, options `mixedLayer`),
 * started afresh in each column and advanced one physics step: h is the
 * boundary-layer top above the surface, θ_l and q_t the dσ-weighted
 * means of the layers whose midpoints lie below it, the free troposphere
 * the first layer above it, the subsidence −πσ̇/(ρ g) of the last
 * dynamics stage interpolated to h, the surface fluxes this column's
 * bulk sensible heat and evaporation, and the longwave the DYCOMS-II
 * form (dycomsLongwave) driven by the mixed layer's own liquid water.
 * Where that first layer's θ_v exceeds the cloud top's by less than
 * `minimumJump` (1 K) the layer counts as uncapped and carries no deck:
 * under so weak a jump the Nicholls–Turton efficiency grows as 1/Δθ_v
 * and one step would carry h through the layer it entrains from. The
 * deck covers the mixed layer's cover times `openSea`, with its water
 * path (at most stratusWaterMax) in the same layer and the same
 * two-column blend; the EIS is still diagnosed. The mixed layer's cover,
 * water path and entrainment rate are kept per cell in mlmCover,
 * mlmWater and mlmEntrainment (0 where it does not run).
 *
 * Γ_l: a saturated parcel conserves q_s + q_l, so it condenses −dq_s/dz
 * per metre of ascent. With d ln q_s/dT = L/(R_v T²), d ln q_s/d ln p
 * = −1, dp/dz = −p g/(R T) and the moist lapse rate
 * Γ_m = (g/c_p)(1 + L q_s/(R T))/(1 + L² q_s/(c_p R_v T²)),
 * −dq_s/dz = q_s (L Γ_m/(R_v T²) − g/(R T)), and Γ_l is that times the
 * air density p/(R T): 2.44e-6 kg/m³ per m at 290 K and 950 hPa.
 *
 * Shortwave: the fraction `ozoneAbsorption` of the incoming beam is
 * absorbed aloft. The ozone column follows Lacis & Hansen (1974)
 * (centred at ozoneHeight with width ozoneWidth, heights from σ with the
 * scale height) and the absorbing part of the beam decays through it
 * with the optical depth ozoneOpacity, so the heating peaks above the
 * ozone maximum as it does at the stratopause. Water vapour absorbs the
 * beam below it by the Lacis & Hansen (1974) absorptivity of the water
 * path the beam has crossed — pressure-scaled by √σ and lengthened by
 * their magnification 35/√(1224μ² + 1) — times `vaporAbsorption`, each
 * layer taking what its own vapour adds to the path above it; what is
 * left goes on to the clouds and the surface. Dry air absorbs nothing.
 */
export function waterVaporAbsorptivity(path) {
  return 2.9 * path / (Math.pow(1 + 141.5 * path, 0.635) + 5.925 * path);
}

const DIFFUSE_MU = 0.6;
export const STABILITY_SIGMA = 0.7;

export function nearestLayer(sigmaMid, sigma) {
  let best = 0;
  for (let k = 1; k < sigmaMid.length; k++) if (Math.abs(sigmaMid[k] - sigma) < Math.abs(sigmaMid[best] - sigma)) best = k;
  return best;
}

export function stratusFraction(index, surfaceT) {
  return Math.min(1, Math.max(0, 0.19 + 0.08 * (index - 1))) * Math.min(1, Math.max(0, (surfaceT - 278.15) / 5));
}

function moistLapse(T, qs, cp, R, g, latentHeat) {
  return g / cp * (1 + latentHeat * qs / (R * T)) / (1 + latentHeat * latentHeat * qs / (cp * (R / EPSILON) * T * T));
}

export function inversionStrength(stability, lowerT, upperT, depth, cp, R, g, latentHeat = LATENT_HEAT) {
  const T = 0.5 * (lowerT + upperT);
  return stability - (g / cp - moistLapse(T, saturationHumidity(T, 85000), cp, R, g, latentHeat)) * depth;
}

export function entrainmentIndex(inversion, lowerQ, upperQ, cp, latentHeat = LATENT_HEAT) {
  return inversion - 0.23 * latentHeat / cp * (lowerQ - upperQ);
}

export function adiabaticWaterLapse(T, p, cp, R, g, latentHeat = LATENT_HEAT) {
  const qs = saturationHumidity(T, p), vaporR = R / EPSILON;
  const moist = moistLapse(T, qs, cp, R, g, latentHeat);
  return p / (R * T) * qs * (latentHeat * moist / (vaporR * T * T) - g / (R * T));
}

export function createRadiation(mesh, core, {
  solarConstant = SOLAR_CONSTANT, albedo = 0.07, cloudAbsorption = 130, cloudScattering = 55, stratus = true, stratusIndex = 'eis', stratusScale = 0.15, stratusWaterMax = 0.15, stratusSigma = 0.92,
  mixedLayerDeck = false, mixedLayer: mixedLayerOptions = {},
  window = 0.25, tauEquator = 5.3, tauPole = 1.325, linearFraction = 0.1, gasFraction = 0.2, gasOpticalDepth = 7,
  ozoneAbsorption = 0.03, ozoneHeight = 25e3, ozoneWidth = 5e3, ozoneOpacity = 4, scaleHeight = 7e3, vaporAbsorption = 1,
  exchangeCoefficient = 1.5e-3, exchangeCoefficients = null, gustiness = 3, latentHeat = LATENT_HEAT, vaporCoupling = 0.55, skylight = 0.15, buffers = null,
} = {}) {
  const { K, C, dSigma, sigmaMid, cp, R, g, kappa, exnerLayer, exnerLower, geopotential, piSigmaDot, p0 } = core.diagnostics;
  const { thetaV } = core.arrays;
  const levels = core.levels;
  const vaporFraction = 1 - window - gasFraction;
  const opticalDepth = (lat) => tauEquator + (tauPole - tauEquator) * Math.sin(lat) ** 2;
  const tauCell = Float64Array.from({ length: C }, (_, i) => opticalDepth(mesh.latCell[i]));
  const shape = Float64Array.from({ length: K }, (_, k) => linearFraction * (levels[k + 1] - levels[k]) + (1 - linearFraction) * (levels[k + 1] ** 4 - levels[k] ** 4));
  const ozoneAbove = (sigma) => (sigma <= 0 ? 0 : (1 + Math.exp(-ozoneHeight / ozoneWidth)) / (1 + Math.exp((-scaleHeight * Math.log(sigma) - ozoneHeight) / ozoneWidth)));
  const beamLeft = (sigma) => Math.exp(-ozoneOpacity * ozoneAbove(sigma));
  const ozoneFraction = Float64Array.from({ length: K }, (_, k) => (beamLeft(levels[k]) - beamLeft(levels[k + 1])) / (1 - Math.exp(-ozoneOpacity)));
  const emissivity = new Float64Array(K);
  const cloudEmissivity = new Float64Array(K);
  const vaporEmissivity = new Float64Array(K);
  const mixedEmissivity = new Float64Array(K);
  const surfaceFlux = new Float64Array(C), surfaceDirect = new Float64Array(C);
  const outgoingBuffer = buffers && buffers.outgoing ? buffers.outgoing : new SharedArrayBuffer(8 * C);
  const shortwaveBuffer = buffers && buffers.surfaceShortwave ? buffers.surfaceShortwave : new SharedArrayBuffer(8 * C);
  const outgoing = new Float64Array(outgoingBuffer), surfaceShortwave = new Float64Array(shortwaveBuffer);
  const evaporationBuffer = buffers && buffers.evaporation ? buffers.evaporation : new SharedArrayBuffer(8 * C);
  const evaporation = new Float64Array(evaporationBuffer);
  const stratusBuffer = buffers && buffers.stratus ? buffers.stratus : new SharedArrayBuffer(8 * C);
  const stratusPath = new Float64Array(stratusBuffer);
  const coverBuffer = buffers && buffers.stratusFraction ? buffers.stratusFraction : new SharedArrayBuffer(8 * C);
  const stratusCover = new Float64Array(coverBuffer);
  const indexBuffer = buffers && buffers.stabilityIndex ? buffers.stabilityIndex : new SharedArrayBuffer(8 * C);
  const stabilityIndex = new Float64Array(indexBuffer);
  const mlmCoverBuffer = buffers && buffers.mlmCover ? buffers.mlmCover : new SharedArrayBuffer(8 * C);
  const mlmWaterBuffer = buffers && buffers.mlmWater ? buffers.mlmWater : new SharedArrayBuffer(8 * C);
  const mlmEntrainmentBuffer = buffers && buffers.mlmEntrainment ? buffers.mlmEntrainment : new SharedArrayBuffer(8 * C);
  const mlmCover = new Float64Array(mlmCoverBuffer), mlmWater = new Float64Array(mlmWaterBuffer), mlmEntrainment = new Float64Array(mlmEntrainmentBuffer);
  const shadow = mixedLayerDeck ? createMixedLayer({ cp, R, g, latentHeat, referencePressure: p0, cloudLevels: 8, minimumJump: 1, ...mixedLayerOptions }) : null;
  const shadowLongwave = dycomsLongwave();
  if (stratusIndex !== 'eis' && stratusIndex !== 'ectei') throw new Error(`stratusIndex must be 'eis' or 'ectei', not ${stratusIndex}`);
  const entraining = stratusIndex === 'ectei';
  const stratusLayer = nearestLayer(sigmaMid, stratusSigma), stabilityLayer = nearestLayer(sigmaMid, STABILITY_SIGMA);
  const gasEmissivity = Float64Array.from({ length: K }, (_, k) => 1 - Math.exp(-gasOpticalDepth * (levels[k + 1] - levels[k])));
  const temperature = new Float64Array(K);
  const emitted = new Float64Array(K);
  const netFlux = new Float64Array(K);
  const sun = new Float64Array([1, 0, 0]);
  const budget = { absorbedSolar: 0, outgoingLongwave: 0, sensibleHeat: 0, evaporation: 0, surfaceFlux: 0, insolation: 0, reflectedSolar: 0, cloudReflectance: 0, stratus: 0, stratusFraction: 0, stabilityIndex: NaN, mlmCover: 0, mlmWater: 0, mlmEntrainment: 0 };
  const sky = { absorbed: 0, down: 0, direct: 0, reflectance: 0 }, decked = { absorbed: 0, down: 0, direct: 0, reflectance: 0 };

  function setTime(t) {
    sunDirection(t, sun);
  }

  function cosZenith(i) {
    return Math.max(0, mesh.xCell[3 * i] * sun[0] + mesh.xCell[3 * i + 1] * sun[1] + mesh.xCell[3 * i + 2] * sun[2]);
  }

  function insolation(i) {
    return solarConstant * cosZenith(i);
  }

  function band(fraction, eps, surfaceEmission) {
    for (let k = 0; k < K; k++) emitted[k] = fraction * eps[k] * STEFAN_BOLTZMANN * temperature[k] ** 4;
    let down = 0;
    for (let k = 0; k < K; k++) {
      netFlux[k] += eps[k] * down - 2 * emitted[k];
      down = down * (1 - eps[k]) + emitted[k];
    }
    let up = fraction * surfaceEmission;
    for (let k = K - 1; k >= 0; k--) {
      netFlux[k] += eps[k] * up;
      up = up * (1 - eps[k]) + emitted[k];
    }
    return [up, down];
  }

  function deckWater(lcl, base, mixedDepth) {
    const thickness = mixedDepth - base;
    if (thickness <= 0) return 0;
    return Math.min(stratusWaterMax, stratusScale * 0.5 * adiabaticWaterLapse(lcl.temperature, lcl.pressure, cp, R, g, latentHeat) * thickness * thickness);
  }

  function shadowDeck(i, pi, theta, q, qc, mixedDepth, sensible, evaporation, dt) {
    const bottom = (K - 1) * C + i;
    const surface = geopotential[bottom] - cp * thetaV[bottom] * (exnerLower[bottom] - exnerLayer[bottom]);
    const h = mixedDepth + (geopotential[bottom] - surface) / g;
    let weight = 0, heat = 0, water = 0, k = K - 1;
    for (; k >= 0 && geopotential[k * C + i] - surface < g * h; k--) {
      const idx = k * C + i, cloud = qc ? Math.max(0, qc[idx]) : 0;
      heat += dSigma[k] * (theta[idx] - latentHeat * cloud / (cp * exnerLayer[idx]));
      water += dSigma[k] * (Math.max(0, q[idx]) + cloud);
      weight += dSigma[k];
    }
    if (k < 1) return false;
    const above = k * C + i, aboveCloud = qc ? Math.max(0, qc[above]) : 0;
    const interfaceHeight = (m) => (geopotential[m * C + i] + cp * thetaV[m * C + i] * (exnerLayer[m * C + i] - exnerLower[(m - 1) * C + i]) - surface) / g;
    let lowerHeight = 0, lowerFlow = 0, m = K - 1;
    for (; m > k && interfaceHeight(m) < h; m--) { lowerHeight = interfaceHeight(m); lowerFlow = piSigmaDot[m * C + i]; }
    const upperHeight = interfaceHeight(m);
    const flow = lowerFlow + (piSigmaDot[m * C + i] - lowerFlow) * (h - lowerHeight) / (upperHeight - lowerHeight);
    const density = pi * sigmaMid[m] / (R * thetaV[m * C + i] * exnerLayer[m * C + i]);
    const subsidence = -flow / (density * g);
    const forcing = {
      surfacePressure: pi, sensibleHeat: sensible, evaporation, radiation: shadowLongwave, subsidence: () => subsidence,
      thetaLAbove: theta[above] - latentHeat * aboveCloud / (cp * exnerLayer[above]), qtAbove: Math.max(0, q[above]) + aboveCloud,
    };
    const start = { h, thetaL: heat / weight, qt: water / weight };
    const now = shadow.diagnose(start, forcing);
    const next = dt > 0 ? shadow.diagnose(shadow.step(start, forcing, dt, now), forcing) : now;
    if (!(Number.isFinite(next.liquidWaterPath) && Number.isFinite(next.cover) && Number.isFinite(now.entrainment))) return false;
    budget.mlmCover = next.cover;
    budget.mlmWater = next.liquidWaterPath;
    budget.mlmEntrainment = now.entrainment;
    return true;
  }

  function shortwave(out, cloudDepth, mu, surfaceAlbedo, diffuseAlbedo) {
    const reflectance = mu > 0 && cloudDepth > 0 ? cloudDepth / (cloudDepth + 2 * mu) : 0;
    const direct = (1 - skylight) * (cloudDepth > 0 && mu > 0 ? Math.exp(-cloudDepth / mu) : 1);
    const diffuse = 1 - reflectance - direct;
    const returned = cloudDepth > 0 ? cloudDepth / (cloudDepth + 2 * DIFFUSE_MU) : 0;
    const upward = surfaceAlbedo * direct + diffuseAlbedo * diffuse;
    const reflections = returned * upward / (1 - diffuseAlbedo * returned);
    out.absorbed = (1 - surfaceAlbedo) * direct + (1 - diffuseAlbedo) * (diffuse + reflections);
    out.down = direct + diffuse + reflections;
    out.direct = direct;
    out.reflectance = reflectance;
  }

  function column(i, pi, theta, surfaceT, windSpeed, tau0 = tauCell[i], beam = insolation(i), qAir = null, q = null, qc = null, surfaceAlbedo = albedo, diffuseAlbedo = surfaceAlbedo, wetness = 1, exchangeCoefficientAt = exchangeCoefficient, openSea = 0, mixedDepth = 0, dt = 0) {
    const ozoneHeating = beam * ozoneAbsorption;
    const surfaceEmission = STEFAN_BOLTZMANN * surfaceT * surfaceT * surfaceT * surfaceT;
    const coupled = vaporCoupling > 0 && q !== null;
    const bottom = K - 1;
    const airTemperature = theta[bottom * C + i] * exnerLayer[bottom * C + i];
    const airDensity = pi * sigmaMid[bottom] / (R * airTemperature);
    const exchange = airDensity * exchangeCoefficientAt * Math.max(windSpeed, gustiness);
    const sensible = exchange * cp * (surfaceT - airTemperature);
    const evaporation = qAir === null ? 0 : wetness * Math.max(0, exchange * (saturationHumidity(surfaceT, pi) - qAir));
    let fraction = 0, deck = 0, index = NaN;
    budget.mlmCover = 0; budget.mlmWater = 0; budget.mlmEntrainment = 0;
    if (stratus && openSea > 0 && qAir !== null && mixedDepth > 0) {
      const lower = bottom * C + i, upper = stabilityLayer * C + i, lowerT = theta[lower] * exnerLayer[lower];
      const lcl = liftingCondensationLevel(lowerT, qAir, pi * sigmaMid[bottom], kappa);
      if (lcl) {
        const base = Math.max(0, cp * (lowerT - lcl.temperature) / g);
        const inversion = inversionStrength(theta[upper] - theta[lower], lowerT, theta[upper] * exnerLayer[upper], (geopotential[upper] - geopotential[lower]) / g - base, cp, R, g, latentHeat);
        index = entraining ? entrainmentIndex(inversion, qAir, q ? q[upper] : qAir, cp, latentHeat) : inversion;
        if (!shadow) {
          fraction = stratusFraction(index, surfaceT) * openSea;
          if (fraction > 0) deck = deckWater(lcl, base, mixedDepth);
        }
      }
      if (shadow && q && shadowDeck(i, pi, theta, q, qc, mixedDepth, sensible, evaporation, dt)) {
        fraction = budget.mlmCover * openSea;
        if (fraction > 0) deck = Math.min(stratusWaterMax, budget.mlmWater);
      }
      if (deck <= 0) fraction = 0;
    }
    let cloudPath = 0;
    for (let k = 0; k < K; k++) {
      const mass = pi * dSigma[k] / g;
      emissivity[k] = 1 - Math.exp(coupled ? -vaporCoupling * Math.max(0, q[k * C + i]) * mass : -tau0 * shape[k]);
      const water = qc ? Math.max(0, qc[k * C + i]) * mass : 0;
      cloudPath += water;
      cloudEmissivity[k] = water > 0 ? 1 - Math.exp(-cloudAbsorption * water) : 0;
      if (deck > 0 && k === stratusLayer) cloudEmissivity[k] = fraction * (1 - Math.exp(-cloudAbsorption * (water + deck))) + (1 - fraction) * cloudEmissivity[k];
      const clear = 1 - cloudEmissivity[k];
      vaporEmissivity[k] = 1 - (1 - emissivity[k]) * clear;
      mixedEmissivity[k] = 1 - (1 - gasEmissivity[k]) * clear;
      temperature[k] = theta[k * C + i] * exnerLayer[k * C + i];
      netFlux[k] = ozoneHeating * ozoneFraction[k];
    }
    const mu = beam / solarConstant;
    let incident = beam - ozoneHeating, vaporHeating = 0;
    if (vaporAbsorption > 0 && q !== null && mu > 0) {
      const magnification = 35 / Math.sqrt(1224 * mu * mu + 1);
      let path = 0, taken = 0;
      for (let k = 0; k < K; k++) {
        path += Math.max(0, q[k * C + i]) * pi * dSigma[k] / g * Math.sqrt(sigmaMid[k]) * 0.1 * magnification;
        const through = vaporAbsorption * waterVaporAbsorptivity(path);
        netFlux[k] += incident * (through - taken);
        taken = through;
      }
      vaporHeating = incident * taken;
      incident -= vaporHeating;
    }
    shortwave(sky, cloudScattering * cloudPath, mu, surfaceAlbedo, diffuseAlbedo);
    if (deck > 0) {
      shortwave(decked, cloudScattering * (cloudPath + deck), mu, surfaceAlbedo, diffuseAlbedo);
      sky.absorbed = fraction * decked.absorbed + (1 - fraction) * sky.absorbed;
      sky.down = fraction * decked.down + (1 - fraction) * sky.down;
      sky.direct = fraction * decked.direct + (1 - fraction) * sky.direct;
      sky.reflectance = fraction * decked.reflectance + (1 - fraction) * sky.reflectance;
    }
    const absorbedSolar = incident * sky.absorbed;
    const [outVapor, backVapor] = band(vaporFraction, vaporEmissivity, surfaceEmission);
    const [outGas, backGas] = band(gasFraction, mixedEmissivity, surfaceEmission);
    const [outWindow, backWindow] = band(window, cloudEmissivity, surfaceEmission);
    const outgoing = outVapor + outGas + outWindow;
    const back = backVapor + backGas + backWindow;
    netFlux[bottom] += sensible;
    const net = absorbedSolar - surfaceEmission + back - sensible - latentHeat * evaporation;
    budget.absorbedSolar = absorbedSolar + ozoneHeating + vaporHeating;
    budget.outgoingLongwave = outgoing;
    budget.sensibleHeat = sensible;
    budget.evaporation = evaporation;
    budget.surfaceFlux = net;
    budget.insolation = beam;
    budget.reflectedSolar = incident - absorbedSolar;
    budget.surfaceShortwave = incident * sky.down;
    budget.surfaceDirect = incident * sky.direct;
    budget.cloudReflectance = sky.reflectance;
    budget.stratus = deck;
    budget.stratusFraction = fraction;
    budget.stabilityIndex = index;
    return net;
  }

  /*
   * Heating tendencies of the layers and the evaporation tendency of the
   * lowest layer for the cells in range; the net surface flux of each
   * cell is left in `surfaceFlux` for the surface model to apply, with
   * the sunlight reaching the surface in `surfaceShortwave`, of which
   * `surfaceDirect` is the direct beam, and the stratocumulus deck's
   * water path and cover in `stratus` and `stratusFraction`, with the
   * EIS or ECTEI its cover follows in `stabilityIndex` (NaN where the
   * deck is not diagnosed: over land or full ice, without humidity or a
   * boundary layer, or with `stratus: false`). `depth` is the height of
   * each cell's boundary-layer top; `dt` is the physics step the
   * mixed-layer deck advances by.
   */
  function apply(state, out, windSpeed, totals, iFrom = 0, iTo = C, surfaceAlbedo = null, diffuseAlbedo = null, wetness = null, openSea = null, depth = null, dt = 0) {
    const [pi, theta, , surfaceT] = state;
    const [, dTheta] = out;
    const q = state[4] ?? null, dQ = out[4] ?? null, qc = state[5] ?? null;
    const bottom = (K - 1) * C;
    if (totals) for (const name of ['absorbedSolar', 'outgoingLongwave', 'sensibleHeat', 'evaporation', 'insolation', 'reflectedSolar']) totals[name] = 0;
    for (let i = iFrom; i < iTo; i++) {
      surfaceFlux[i] = column(i, pi[i], theta, surfaceT[i], windSpeed[i], tauCell[i], insolation(i), q && dQ ? q[bottom + i] : null, q && dQ ? q : null, q && dQ ? qc : null, surfaceAlbedo ? surfaceAlbedo[i] : albedo, diffuseAlbedo ? diffuseAlbedo[i] : surfaceAlbedo ? surfaceAlbedo[i] : albedo, wetness ? wetness[i] : 1, exchangeCoefficients ? exchangeCoefficients[i] : exchangeCoefficient, openSea ? openSea[i] : 0, depth ? depth[i] - geopotential[bottom + i] / g : 0, dt);
      outgoing[i] = budget.outgoingLongwave;
      stratusPath[i] = budget.stratus;
      stratusCover[i] = budget.stratusFraction;
      stabilityIndex[i] = budget.stabilityIndex;
      mlmCover[i] = budget.mlmCover;
      mlmWater[i] = budget.mlmWater;
      mlmEntrainment[i] = budget.mlmEntrainment;
      evaporation[i] = budget.evaporation;
      surfaceShortwave[i] = budget.surfaceShortwave;
      surfaceDirect[i] = budget.surfaceDirect;
      for (let k = 0; k < K; k++) {
        const massPerArea = pi[i] * dSigma[k] / g;
        dTheta[k * C + i] += netFlux[k] / (cp * massPerArea) / exnerLayer[k * C + i];
      }
      if (q && dQ) dQ[bottom + i] += budget.evaporation * g / (pi[i] * dSigma[K - 1]);
      if (totals) {
        const a = mesh.areaCell[i];
        totals.absorbedSolar += a * budget.absorbedSolar;
        totals.outgoingLongwave += a * budget.outgoingLongwave;
        totals.sensibleHeat += a * budget.sensibleHeat;
        totals.evaporation += a * budget.evaporation;
        totals.insolation += a * budget.insolation;
        totals.reflectedSolar += a * budget.reflectedSolar;
      }
    }
  }

  return { setTime, sun, cosZenith, insolation, column, apply, layerFlux: netFlux, surfaceFlux, outgoing, surfaceShortwave, surfaceDirect, evaporation, stratus: stratusPath, stratusFraction: stratusCover, stabilityIndex, mlmCover, mlmWater, mlmEntrainment, stratusLayer, stabilityLayer, budget, emissivity, opticalDepth, ozoneFraction, shared: { outgoing: outgoingBuffer, surfaceShortwave: shortwaveBuffer, evaporation: evaporationBuffer, stratus: stratusBuffer, stratusFraction: coverBuffer, stabilityIndex: indexBuffer, mlmCover: mlmCoverBuffer, mlmWater: mlmWaterBuffer, mlmEntrainment: mlmEntrainmentBuffer } };
}
