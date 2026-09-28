import { LATENT_HEAT, saturationHumidity } from './moist.module.js';
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
 * the column carries a diagnostic deck of cloud water stratusWater × f
 * in the layer nearest σ = stratusSigma, radiating like the condensed
 * water beside it but never added to qc. f follows Klein & Hartmann
 * (1993), 0.057 LTS − 0.556 clamped to [0, 1], with the lower-
 * tropospheric stability LTS the potential temperature of the layer
 * nearest σ = 0.7 less the lowest layer's, times a ramp from 0 at a
 * 5 °C surface to 1 at 10 °C that keeps the deck off polar seas, times
 * `openSea`. `stratus: false` removes it.
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

export function stratusFraction(stability, surfaceT) {
  return Math.min(1, Math.max(0, 0.057 * stability - 0.556)) * Math.min(1, Math.max(0, (surfaceT - 278.15) / 5));
}

export function createRadiation(mesh, core, {
  solarConstant = SOLAR_CONSTANT, albedo = 0.07, cloudAbsorption = 130, cloudScattering = 55, stratus = true, stratusWater = 0.004, stratusSigma = 0.92,
  window = 0.25, tauEquator = 5.3, tauPole = 1.325, linearFraction = 0.1, gasFraction = 0.2, gasOpticalDepth = 7,
  ozoneAbsorption = 0.03, ozoneHeight = 25e3, ozoneWidth = 5e3, ozoneOpacity = 4, scaleHeight = 7e3, vaporAbsorption = 1,
  exchangeCoefficient = 1.5e-3, exchangeCoefficients = null, gustiness = 3, latentHeat = LATENT_HEAT, vaporCoupling = 0.55, skylight = 0.15, buffers = null,
} = {}) {
  const { K, C, dSigma, sigmaMid, cp, R, g, exnerLayer } = core.diagnostics;
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
  const stratusLayer = nearestLayer(sigmaMid, stratusSigma), stabilityLayer = nearestLayer(sigmaMid, STABILITY_SIGMA);
  const gasEmissivity = Float64Array.from({ length: K }, (_, k) => 1 - Math.exp(-gasOpticalDepth * (levels[k + 1] - levels[k])));
  const temperature = new Float64Array(K);
  const emitted = new Float64Array(K);
  const netFlux = new Float64Array(K);
  const sun = new Float64Array([1, 0, 0]);
  const budget = { absorbedSolar: 0, outgoingLongwave: 0, sensibleHeat: 0, evaporation: 0, surfaceFlux: 0, insolation: 0, reflectedSolar: 0, cloudReflectance: 0, stratus: 0 };

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

  function column(i, pi, theta, surfaceT, windSpeed, tau0 = tauCell[i], beam = insolation(i), qAir = null, q = null, qc = null, surfaceAlbedo = albedo, diffuseAlbedo = surfaceAlbedo, wetness = 1, exchangeCoefficientAt = exchangeCoefficient, openSea = 0) {
    const ozoneHeating = beam * ozoneAbsorption;
    const surfaceEmission = STEFAN_BOLTZMANN * surfaceT * surfaceT * surfaceT * surfaceT;
    const coupled = vaporCoupling > 0 && q !== null;
    const deck = stratus && openSea > 0 ? stratusWater * stratusFraction(theta[stabilityLayer * C + i] - theta[(K - 1) * C + i], surfaceT) * openSea : 0;
    let cloudPath = 0;
    for (let k = 0; k < K; k++) {
      const mass = pi * dSigma[k] / g;
      emissivity[k] = 1 - Math.exp(coupled ? -vaporCoupling * Math.max(0, q[k * C + i]) * mass : -tau0 * shape[k]);
      let water = qc ? Math.max(0, qc[k * C + i]) * mass : 0;
      if (deck > 0 && k === stratusLayer) water += deck;
      cloudPath += water;
      cloudEmissivity[k] = water > 0 ? 1 - Math.exp(-cloudAbsorption * water) : 0;
      const clear = 1 - cloudEmissivity[k];
      vaporEmissivity[k] = 1 - (1 - emissivity[k]) * clear;
      mixedEmissivity[k] = 1 - (1 - gasEmissivity[k]) * clear;
      temperature[k] = theta[k * C + i] * exnerLayer[k * C + i];
      netFlux[k] = ozoneHeating * ozoneFraction[k];
    }
    const mu = beam / solarConstant;
    const cloudDepth = cloudScattering * cloudPath;
    const reflectance = mu > 0 && cloudDepth > 0 ? cloudDepth / (cloudDepth + 2 * mu) : 0;
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
    const direct = (1 - skylight) * (cloudDepth > 0 && mu > 0 ? Math.exp(-cloudDepth / mu) : 1);
    const diffuse = 1 - reflectance - direct;
    const returned = cloudDepth > 0 ? cloudDepth / (cloudDepth + 2 * DIFFUSE_MU) : 0;
    const upward = surfaceAlbedo * direct + diffuseAlbedo * diffuse;
    const absorbedSolar = incident * ((1 - surfaceAlbedo) * direct + (1 - diffuseAlbedo) * (diffuse + returned * upward / (1 - diffuseAlbedo * returned)));
    const [outVapor, backVapor] = band(vaporFraction, vaporEmissivity, surfaceEmission);
    const [outGas, backGas] = band(gasFraction, mixedEmissivity, surfaceEmission);
    const [outWindow, backWindow] = band(window, cloudEmissivity, surfaceEmission);
    const outgoing = outVapor + outGas + outWindow;
    const back = backVapor + backGas + backWindow;
    const bottom = K - 1;
    const airTemperature = temperature[bottom];
    const airDensity = pi * sigmaMid[bottom] / (R * airTemperature);
    const exchange = airDensity * exchangeCoefficientAt * Math.max(windSpeed, gustiness);
    const sensible = exchange * cp * (surfaceT - airTemperature);
    const evaporation = qAir === null ? 0 : wetness * Math.max(0, exchange * (saturationHumidity(surfaceT, pi) - qAir));
    netFlux[bottom] += sensible;
    const net = absorbedSolar - surfaceEmission + back - sensible - latentHeat * evaporation;
    budget.absorbedSolar = absorbedSolar + ozoneHeating + vaporHeating;
    budget.outgoingLongwave = outgoing;
    budget.sensibleHeat = sensible;
    budget.evaporation = evaporation;
    budget.surfaceFlux = net;
    budget.insolation = beam;
    budget.reflectedSolar = incident - absorbedSolar;
    budget.surfaceShortwave = incident * (direct + diffuse + returned * upward / (1 - diffuseAlbedo * returned));
    budget.surfaceDirect = incident * direct;
    budget.cloudReflectance = reflectance;
    budget.stratus = deck;
    return net;
  }

  /*
   * Heating tendencies of the layers and the evaporation tendency of the
   * lowest layer for the cells in range; the net surface flux of each
   * cell is left in `surfaceFlux` for the surface model to apply, with
   * the sunlight reaching the surface in `surfaceShortwave`, of which
   * `surfaceDirect` is the direct beam, and the stratocumulus deck's
   * water path in `stratus`.
   */
  function apply(state, out, windSpeed, totals, iFrom = 0, iTo = C, surfaceAlbedo = null, diffuseAlbedo = null, wetness = null, openSea = null) {
    const [pi, theta, , surfaceT] = state;
    const [, dTheta] = out;
    const q = state[4] ?? null, dQ = out[4] ?? null, qc = state[5] ?? null;
    const bottom = (K - 1) * C;
    if (totals) for (const name of ['absorbedSolar', 'outgoingLongwave', 'sensibleHeat', 'evaporation', 'insolation', 'reflectedSolar']) totals[name] = 0;
    for (let i = iFrom; i < iTo; i++) {
      surfaceFlux[i] = column(i, pi[i], theta, surfaceT[i], windSpeed[i], tauCell[i], insolation(i), q && dQ ? q[bottom + i] : null, q && dQ ? q : null, q && dQ ? qc : null, surfaceAlbedo ? surfaceAlbedo[i] : albedo, diffuseAlbedo ? diffuseAlbedo[i] : surfaceAlbedo ? surfaceAlbedo[i] : albedo, wetness ? wetness[i] : 1, exchangeCoefficients ? exchangeCoefficients[i] : exchangeCoefficient, openSea ? openSea[i] : 0);
      outgoing[i] = budget.outgoingLongwave;
      stratusPath[i] = budget.stratus;
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

  return { setTime, sun, cosZenith, insolation, column, apply, layerFlux: netFlux, surfaceFlux, outgoing, surfaceShortwave, surfaceDirect, evaporation, stratus: stratusPath, stratusLayer, stabilityLayer, budget, emissivity, opticalDepth, ozoneFraction, shared: { outgoing: outgoingBuffer, surfaceShortwave: shortwaveBuffer, evaporation: evaporationBuffer, stratus: stratusBuffer } };
}
