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
 * Shortwave: the fraction `ozoneAbsorption` of the incoming beam is
 * absorbed aloft. The ozone column follows Lacis & Hansen (1974)
 * (centred at ozoneHeight with width ozoneWidth, heights from σ with the
 * scale height) and the absorbing part of the beam decays through it
 * with the optical depth ozoneOpacity, so the heating peaks above the
 * ozone maximum as it does at the stratopause. The surface absorbs
 * (1 − albedo) of what remains.
 */
export function createRadiation(mesh, core, {
  solarConstant = SOLAR_CONSTANT, albedo = 0.3, surfaceHeatCapacity = 2.1e7,
  window = 0.25, tauEquator = 5.3, tauPole = 1.325, linearFraction = 0.1, gasFraction = 0.2, gasOpticalDepth = 5,
  ozoneAbsorption = 0.03, ozoneHeight = 25e3, ozoneWidth = 5e3, ozoneOpacity = 4, scaleHeight = 7e3,
  exchangeCoefficient = 1.5e-3, gustiness = 3, latentHeat = LATENT_HEAT, vaporCoupling = 2,
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
  const gasEmissivity = Float64Array.from({ length: K }, (_, k) => 1 - Math.exp(-gasOpticalDepth * (levels[k + 1] - levels[k])));
  const temperature = new Float64Array(K);
  const emitted = new Float64Array(K);
  const netFlux = new Float64Array(K);
  const sun = new Float64Array([1, 0, 0]);
  const budget = { absorbedSolar: 0, outgoingLongwave: 0, sensibleHeat: 0, evaporation: 0, surfaceFlux: 0 };

  function setTime(t) {
    sunDirection(t, sun);
  }

  function insolation(i) {
    const cosZenith = mesh.xCell[3 * i] * sun[0] + mesh.xCell[3 * i + 1] * sun[1] + mesh.xCell[3 * i + 2] * sun[2];
    return solarConstant * Math.max(0, cosZenith);
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

  function column(i, pi, theta, surfaceT, windSpeed, tau0 = tauCell[i], beam = insolation(i), qAir = null, q = null) {
    const ozoneHeating = beam * ozoneAbsorption;
    const absorbedSolar = (1 - albedo) * (beam - ozoneHeating);
    const surfaceEmission = STEFAN_BOLTZMANN * surfaceT * surfaceT * surfaceT * surfaceT;
    const coupled = vaporCoupling > 0 && q !== null;
    for (let k = 0; k < K; k++) {
      emissivity[k] = 1 - Math.exp(coupled ? -vaporCoupling * Math.max(0, q[k * C + i]) * pi * dSigma[k] / g : -tau0 * shape[k]);
      temperature[k] = theta[k * C + i] * exnerLayer[k * C + i];
      netFlux[k] = ozoneHeating * ozoneFraction[k];
    }
    const [outVapor, backVapor] = band(vaporFraction, emissivity, surfaceEmission);
    const [outGas, backGas] = band(gasFraction, gasEmissivity, surfaceEmission);
    const outgoing = outVapor + outGas + window * surfaceEmission;
    const back = backVapor + backGas;
    const bottom = K - 1;
    const airTemperature = temperature[bottom];
    const airDensity = pi * sigmaMid[bottom] / (R * airTemperature);
    const exchange = airDensity * exchangeCoefficient * Math.max(windSpeed, gustiness);
    const sensible = exchange * cp * (surfaceT - airTemperature);
    const evaporation = qAir === null ? 0 : Math.max(0, exchange * (saturationHumidity(surfaceT, pi) - qAir));
    netFlux[bottom] += sensible;
    const surfaceFlux = absorbedSolar - surfaceEmission + back - sensible - latentHeat * evaporation;
    budget.absorbedSolar = absorbedSolar + ozoneHeating;
    budget.outgoingLongwave = outgoing;
    budget.sensibleHeat = sensible;
    budget.evaporation = evaporation;
    budget.surfaceFlux = surfaceFlux;
    return surfaceFlux;
  }

  function apply(state, out, windSpeed, totals, iFrom = 0, iTo = C) {
    const [pi, theta, , surfaceT] = state;
    const [, dTheta, , dSurfaceT] = out;
    const q = state[4] ?? null, dQ = out[4] ?? null;
    const bottom = (K - 1) * C;
    if (totals) { totals.absorbedSolar = 0; totals.outgoingLongwave = 0; totals.sensibleHeat = 0; totals.evaporation = 0; }
    for (let i = iFrom; i < iTo; i++) {
      const surfaceFlux = column(i, pi[i], theta, surfaceT[i], windSpeed[i], tauCell[i], insolation(i), q && dQ ? q[bottom + i] : null, q && dQ ? q : null);
      for (let k = 0; k < K; k++) {
        const massPerArea = pi[i] * dSigma[k] / g;
        dTheta[k * C + i] += netFlux[k] / (cp * massPerArea) / exnerLayer[k * C + i];
      }
      if (q && dQ) dQ[bottom + i] += budget.evaporation * g / (pi[i] * dSigma[K - 1]);
      dSurfaceT[i] = surfaceFlux / surfaceHeatCapacity;
      if (totals) {
        const a = mesh.areaCell[i];
        totals.absorbedSolar += a * budget.absorbedSolar;
        totals.outgoingLongwave += a * budget.outgoingLongwave;
        totals.sensibleHeat += a * budget.sensibleHeat;
        totals.evaporation += a * budget.evaporation;
      }
    }
  }

  return { setTime, sun, insolation, column, apply, layerFlux: netFlux, budget, emissivity, opticalDepth, ozoneFraction };
}
