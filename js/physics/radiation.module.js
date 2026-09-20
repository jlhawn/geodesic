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
 * Shortwave: the fraction `ozoneAbsorption` of the incoming beam is
 * absorbed aloft. The ozone column follows Lacis & Hansen (1974)
 * (centred at ozoneHeight with width ozoneWidth, heights from σ with the
 * scale height) and the absorbing part of the beam decays through it
 * with the optical depth ozoneOpacity, so the heating peaks above the
 * ozone maximum as it does at the stratopause.
 */
const DIFFUSE_MU = 0.6;

export function createRadiation(mesh, core, {
  solarConstant = SOLAR_CONSTANT, albedo = 0.07, cloudAbsorption = 130, cloudScattering = 120,
  window = 0.25, tauEquator = 5.3, tauPole = 1.325, linearFraction = 0.1, gasFraction = 0.2, gasOpticalDepth = 5,
  ozoneAbsorption = 0.03, ozoneHeight = 25e3, ozoneWidth = 5e3, ozoneOpacity = 4, scaleHeight = 7e3,
  exchangeCoefficient = 1.5e-3, gustiness = 3, latentHeat = LATENT_HEAT, vaporCoupling = 0.55, skylight = 0.15,
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
  const surfaceFlux = new Float64Array(C);
  const gasEmissivity = Float64Array.from({ length: K }, (_, k) => 1 - Math.exp(-gasOpticalDepth * (levels[k + 1] - levels[k])));
  const temperature = new Float64Array(K);
  const emitted = new Float64Array(K);
  const netFlux = new Float64Array(K);
  const sun = new Float64Array([1, 0, 0]);
  const budget = { absorbedSolar: 0, outgoingLongwave: 0, sensibleHeat: 0, evaporation: 0, surfaceFlux: 0, insolation: 0, reflectedSolar: 0, cloudReflectance: 0 };

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
    let carry = fraction * surfaceEmission;
    for (let k = K - 1; k >= 0; k--) {
      netFlux[k] += eps[k] * carry;
      carry *= 1 - eps[k];
    }
    let outgoing = carry, back = 0;
    for (let k = 0; k < K; k++) {
      netFlux[k] -= 2 * emitted[k];
      let down = emitted[k];
      for (let j = k + 1; j < K; j++) {
        netFlux[j] += eps[j] * down;
        down *= 1 - eps[j];
      }
      back += down;
      let up = emitted[k];
      for (let j = k - 1; j >= 0; j--) {
        netFlux[j] += eps[j] * up;
        up *= 1 - eps[j];
      }
      outgoing += up;
    }
    return [outgoing, back];
  }

  function column(i, pi, theta, surfaceT, windSpeed, tau0 = tauCell[i], beam = insolation(i), qAir = null, q = null, qc = null, surfaceAlbedo = albedo, diffuseAlbedo = surfaceAlbedo) {
    const ozoneHeating = beam * ozoneAbsorption;
    const surfaceEmission = STEFAN_BOLTZMANN * surfaceT * surfaceT * surfaceT * surfaceT;
    const coupled = vaporCoupling > 0 && q !== null;
    let cloudPath = 0;
    for (let k = 0; k < K; k++) {
      const mass = pi * dSigma[k] / g;
      emissivity[k] = 1 - Math.exp(coupled ? -vaporCoupling * Math.max(0, q[k * C + i]) * mass : -tau0 * shape[k]);
      const water = qc ? Math.max(0, qc[k * C + i]) * mass : 0;
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
    const incident = beam - ozoneHeating;
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
    const exchange = airDensity * exchangeCoefficient * Math.max(windSpeed, gustiness);
    const sensible = exchange * cp * (surfaceT - airTemperature);
    const evaporation = qAir === null ? 0 : Math.max(0, exchange * (saturationHumidity(surfaceT, pi) - qAir));
    netFlux[bottom] += sensible;
    const net = absorbedSolar - surfaceEmission + back - sensible - latentHeat * evaporation;
    budget.absorbedSolar = absorbedSolar + ozoneHeating;
    budget.outgoingLongwave = outgoing;
    budget.sensibleHeat = sensible;
    budget.evaporation = evaporation;
    budget.surfaceFlux = net;
    budget.insolation = beam;
    budget.reflectedSolar = incident - absorbedSolar;
    budget.cloudReflectance = reflectance;
    return net;
  }

  /*
   * Heating tendencies of the layers and the evaporation tendency of the
   * lowest layer for the cells in range; the net surface flux of each
   * cell is left in `surfaceFlux` for the surface model to apply.
   */
  function apply(state, out, windSpeed, totals, iFrom = 0, iTo = C, surfaceAlbedo = null, diffuseAlbedo = null) {
    const [pi, theta, , surfaceT] = state;
    const [, dTheta] = out;
    const q = state[4] ?? null, dQ = out[4] ?? null, qc = state[5] ?? null;
    const bottom = (K - 1) * C;
    if (totals) for (const name of ['absorbedSolar', 'outgoingLongwave', 'sensibleHeat', 'evaporation', 'insolation', 'reflectedSolar']) totals[name] = 0;
    for (let i = iFrom; i < iTo; i++) {
      surfaceFlux[i] = column(i, pi[i], theta, surfaceT[i], windSpeed[i], tauCell[i], insolation(i), q && dQ ? q[bottom + i] : null, q && dQ ? q : null, q && dQ ? qc : null, surfaceAlbedo ? surfaceAlbedo[i] : albedo, diffuseAlbedo ? diffuseAlbedo[i] : surfaceAlbedo ? surfaceAlbedo[i] : albedo);
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

  return { setTime, sun, cosZenith, insolation, column, apply, layerFlux: netFlux, surfaceFlux, budget, emissivity, opticalDepth, ozoneFraction };
}
