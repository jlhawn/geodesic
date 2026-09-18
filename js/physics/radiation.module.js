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
 * Gray longwave column with a slab-ocean surface and a bulk sensible heat
 * flux. Each layer's emissivity follows its mass so a standard atmosphere
 * has the prescribed top-of-atmosphere emissivity; every emission is
 * either absorbed on its way or leaves through the top or reaches the
 * surface, so the layer and surface energy fluxes sum exactly to absorbed
 * solar minus outgoing longwave.
 */
export function createRadiation(mesh, core, {
  solarConstant = SOLAR_CONSTANT, albedo = 0.3, surfaceHeatCapacity = 2.1e7, toaEmissivity = 0.78,
  emissivityReferencePressure = 101325, exchangeCoefficient = 1.5e-3, gustiness = 3,
} = {}) {
  const { K, C, dSigma, sigmaMid, cp, R, g, exnerLayer } = core.diagnostics;
  const decay = Math.log(1 - toaEmissivity) / emissivityReferencePressure;
  const emissivity = new Float64Array(K);
  const temperature = new Float64Array(K);
  const emitted = new Float64Array(K);
  const netFlux = new Float64Array(K);
  const sun = new Float64Array([1, 0, 0]);
  const budget = { absorbedSolar: 0, outgoingLongwave: 0, sensibleHeat: 0, surfaceFlux: 0 };

  function setTime(t) {
    sunDirection(t, sun);
  }

  function insolation(i) {
    const cosZenith = mesh.xCell[3 * i] * sun[0] + mesh.xCell[3 * i + 1] * sun[1] + mesh.xCell[3 * i + 2] * sun[2];
    return solarConstant * Math.max(0, cosZenith);
  }

  function column(i, pi, theta, surfaceT, windSpeed) {
    const absorbedSolar = (1 - albedo) * insolation(i);
    const surfaceEmission = STEFAN_BOLTZMANN * surfaceT * surfaceT * surfaceT * surfaceT;
    for (let k = 0; k < K; k++) {
      emissivity[k] = 1 - Math.exp(decay * pi * dSigma[k]);
      const ex = exnerLayer[k * C + i];
      temperature[k] = theta[k * C + i] * ex;
      emitted[k] = emissivity[k] * STEFAN_BOLTZMANN * temperature[k] ** 4;
    }
    let carry = surfaceEmission;
    for (let k = K - 1; k >= 0; k--) {
      netFlux[k] = emissivity[k] * carry;
      carry *= 1 - emissivity[k];
    }
    let outgoing = carry;
    let back = 0;
    for (let k = 0; k < K; k++) {
      netFlux[k] -= 2 * emitted[k];
      let down = emitted[k];
      for (let j = k + 1; j < K; j++) {
        netFlux[j] += emissivity[j] * down;
        down *= 1 - emissivity[j];
      }
      back += down;
      let up = emitted[k];
      for (let j = k - 1; j >= 0; j--) {
        netFlux[j] += emissivity[j] * up;
        up *= 1 - emissivity[j];
      }
      outgoing += up;
    }
    const bottom = K - 1;
    const airTemperature = temperature[bottom];
    const airDensity = pi * sigmaMid[bottom] / (R * airTemperature);
    const sensible = airDensity * cp * exchangeCoefficient * Math.max(windSpeed, gustiness) * (surfaceT - airTemperature);
    netFlux[bottom] += sensible;
    const surfaceFlux = absorbedSolar - surfaceEmission + back - sensible;
    budget.absorbedSolar = absorbedSolar;
    budget.outgoingLongwave = outgoing;
    budget.sensibleHeat = sensible;
    budget.surfaceFlux = surfaceFlux;
    return surfaceFlux;
  }

  function apply(state, out, windSpeed, totals) {
    const [pi, theta, , surfaceT] = state;
    const [, dTheta, , dSurfaceT] = out;
    if (totals) { totals.absorbedSolar = 0; totals.outgoingLongwave = 0; totals.sensibleHeat = 0; }
    for (let i = 0; i < C; i++) {
      const surfaceFlux = column(i, pi[i], theta, surfaceT[i], windSpeed[i]);
      for (let k = 0; k < K; k++) {
        const massPerArea = pi[i] * dSigma[k] / g;
        dTheta[k * C + i] += netFlux[k] / (cp * massPerArea) / exnerLayer[k * C + i];
      }
      dSurfaceT[i] = surfaceFlux / surfaceHeatCapacity;
      if (totals) {
        const a = mesh.areaCell[i];
        totals.absorbedSolar += a * budget.absorbedSolar;
        totals.outgoingLongwave += a * budget.outgoingLongwave;
        totals.sensibleHeat += a * budget.sensibleHeat;
      }
    }
  }

  return { setTime, sun, insolation, column, apply, layerFlux: netFlux, budget, emissivity };
}
