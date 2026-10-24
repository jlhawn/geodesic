/*
 * One step of the model's own moist adjustment (moist.adjust: large-scale
 * condensation, the test parcel, the entraining plume and its closure,
 * autoconversion and rain) on cell 0 of a 42-cell createModel mesh shrunk
 * so its mean cell is `spacing` across, on the bl36 layers, with Jordan's
 * (1958) mean West Indies hurricane-season sounding (hPa, °C, and g/kg
 * taken as specific humidity, as test/convection.test.mjs places it) over
 * that test's plume surface. `humidity` scales the vapor above 850 hPa and
 * `warming` heats the lowest kilometer at its relative humidity (at its
 * specific humidity with moistenWarming false), both capped at saturation.
 * Every call starts from a column with no past: the closure's saved
 * subcloud state and the cumulus cloud are cleared, so the Bechtold
 * closure sees no boundary-layer forcing and the cover is the step's own,
 * which the model's 30-minute cumulus memory relaxes toward. Large-scale
 * heating is L times the condensate it makes, so the rest of each layer's
 * vapor change is convection's. Layers run top first; heating is K/day,
 * moistening g/kg/day, mass fluxes kg/m²/s at the interfaces, rain mm/day.
 */
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createModel } from '../js/model.module.js';
import { sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { saturationHumidity, LATENT_HEAT } from '../js/physics/moist.module.js';
import { PUBLISHED_RADIATION, PUBLISHED_MOIST } from './published.module.js';

export const JORDAN = [[1015, 26.3, 18.9], [1000, 25.4, 18.4], [950, 22.2, 16.6], [900, 19.6, 14.4], [850, 17.3, 12.4], [800, 14.8, 10.5], [750, 12.0, 8.8],
  [700, 8.9, 7.4], [650, 5.6, 6.0], [600, 1.8, 4.7], [550, -2.3, 3.6], [500, -6.9, 2.6], [450, -11.9, 1.8], [400, -17.3, 1.2], [350, -24.4, 0.66],
  [300, -32.4, 0.35], [250, -42.1, 0.14], [200, -54.0, 0.04], [175, -60.2, 0.02], [150, -67.2, 0.01], [125, -73.0, 0.005], [100, -73.3, 0.003],
  [70, -69, 0.003], [50, -62, 0.003], [30, -55, 0.003], [10, -45, 0.003], [1, -20, 0.003]];
export const SURFACE_PRESSURE = 101500;
export const SPACING = 55.8e3;
export const SPACINGS = [55.8e3, 111.6e3];
export const STEP = 168.75;
export const HUMID_ABOVE = 850e2;
export const WARMED_DEPTH = 1000;
export const SURFACE = { depthAbove: 500, buoyancyFlux: 4e-4, friction: 0.25, deckGate: 0.3 };
export const LEVELS = sigmaInterfaces('bl36');

export function jordan(p) {
  if (p >= 100 * JORDAN[0][0]) return { T: 273.15 + JORDAN[0][1], q: 1e-3 * JORDAN[0][2] };
  const x = Math.log(p);
  for (let n = 0; n < JORDAN.length - 1; n++) {
    const [p0, t0, q0] = JORDAN[n], [p1, t1, q1] = JORDAN[n + 1];
    if ((p <= 100 * p0 && p >= 100 * p1) || n === JORDAN.length - 2) {
      const f = Math.max(0, (x - Math.log(100 * p0)) / (Math.log(100 * p1) - Math.log(100 * p0)));
      return { T: 273.15 + t0 + f * (t1 - t0), q: 1e-3 * (q0 + f * (q1 - q0)) };
    }
  }
  return null;
}

export function sounding({ humidity = 1, warming = 0, moistenWarming = true } = {}) {
  return (z, p) => {
    const air = jordan(p), T = air.T + (z <= WARMED_DEPTH ? warming : 0);
    const q = moistenWarming ? air.q * saturationHumidity(T, p) / saturationHumidity(air.T, p) : air.q;
    return { T, q: Math.min(q * (p < HUMID_ABOVE ? humidity : 1), saturationHumidity(T, p)) };
  };
}

const models = new Map();
let meanCell = 0;
function modelFor(spacing, ice) {
  const key = `${spacing}:${ice}`;
  if (models.has(key)) return models.get(key);
  const grid = new Grid(2);
  if (!meanCell) { const probe = buildMesh(grid); meanCell = Math.sqrt(probe.areaCell.reduce((a, b) => a + b, 0) / probe.nCells); }
  const model = createModel(buildMesh(grid, { radius: 6.371e6 * spacing / meanCell }), { ocean: false, levels: LEVELS, radiation: PUBLISHED_RADIATION, moist: { ...PUBLISHED_MOIST, cumulusMemory: 0, plumePhase: ice ? 'mixed' : 'liquid' } });
  const { K, C } = model.core.diagnostics;
  model.moist.trace.convection = new Float64Array(K * C);
  model.moist.trace.largeScale = new Float64Array(K * C);
  models.set(key, model);
  return model;
}

function place(model, profile) {
  const { K, C, sigmaMid, exnerLayer, geopotential, g } = model.core.diagnostics;
  const [pi, theta, , , q, qc] = model.state;
  pi[0] = SURFACE_PRESSURE;
  for (let k = 0; k < K; k++) { theta[k * C] = 300; q[k * C] = 0; qc[k * C] = 0; }
  for (let pass = 0; pass < 4; pass++) {
    model.core.diagnoseColumn(0, pi, theta, q, qc);
    for (let k = K - 1; k >= 0; k--) {
      const idx = k * C, air = profile(geopotential[idx] / g, pi[0] * sigmaMid[k]);
      theta[idx] = air.T / exnerLayer[idx];
      q[idx] = air.q;
    }
  }
  model.core.diagnoseColumn(0, pi, theta, q, qc);
}

function forget(model) {
  const { K, C } = model.core.diagnostics, moist = model.moist;
  for (let k = 0; k < moist.subcloudLayers; k++) moist.subcloudVirtual[k * C] = 0;
  for (let k = 0; k < K; k++) {
    moist.cumulusCover[k * C] = 0; moist.cumulusWater[k * C] = 0;
    moist.trace.convection[k * C] = 0; moist.trace.largeScale[k * C] = 0;
  }
  moist.precipitation[0] = 0; moist.convectivePrecipitation[0] = 0; moist.largeScalePrecipitation[0] = 0; moist.fogPrecipitation[0] = 0;
  moist.plumeBuoyancy.fill(NaN);
  moist.plumeInterfaceT.fill(NaN);
}

export function plumeColumn({ humidity = 1, warming = 0, moistenWarming = true, spacing = SPACING, ice = true, dt = STEP * spacing / SPACING } = {}) {
  const model = modelFor(spacing, ice), moist = model.moist;
  const { K, C, sigmaMid, levels, exnerLayer, geopotential, g, cp, dSigma } = model.core.diagnostics;
  const [pi, theta, , , q, qc] = model.state;
  place(model, sounding({ humidity, warming, moistenWarming }));
  model.boundaryLayer.depth[0] = geopotential[(K - 1) * C] / g + SURFACE.depthAbove;
  model.boundaryLayer.buoyancyFlux[0] = SURFACE.buoyancyFlux;
  model.boundaryLayer.friction[0] = SURFACE.friction;
  model.radiation.mlmGate[0] = SURFACE.deckGate;
  forget(model);
  const ps = pi[0], pressure = new Float64Array(K), height = new Float64Array(K), temperature = new Float64Array(K), vapor = new Float64Array(K), relative = new Float64Array(K), mass = new Float64Array(K);
  for (let k = 0; k < K; k++) {
    const idx = k * C;
    pressure[k] = ps * sigmaMid[k]; height[k] = geopotential[idx] / g; temperature[k] = theta[idx] * exnerLayer[idx];
    vapor[k] = q[idx]; relative[k] = q[idx] / saturationHumidity(temperature[k], pressure[k]); mass[k] = ps * dSigma[k] / g;
  }
  moist.adjust(model.state, 0, 1, dt);
  const perDay = 86400 / dt, deep = moist.deep, ran = moist.cumulusBaseFlux[0] > 0;
  const type = deep.deep ? 'deep' : ran ? 'shallow' : 'none';
  const heating = { convection: new Float64Array(K), largeScale: new Float64Array(K) }, moistening = { convection: new Float64Array(K), largeScale: new Float64Array(K) };
  const cloudWater = Float64Array.from({ length: K }, (_, k) => qc[k * C]), cover = new Float64Array(K), buoyancy = new Float64Array(K), massFlux = new Float64Array(K + 1), downdraft = new Float64Array(K + 1), plumeT = new Float64Array(K + 1).fill(NaN);
  for (let k = 0; k < K; k++) {
    const idx = k * C, largeScale = moist.trace.largeScale[idx], dq = q[idx] - vapor[k];
    heating.convection[k] = moist.trace.convection[idx] * perDay;
    heating.largeScale[k] = largeScale * perDay;
    moistening.largeScale[k] = -cp * largeScale / LATENT_HEAT * perDay * 1e3;
    moistening.convection[k] = (dq + cp * largeScale / LATENT_HEAT) * perDay * 1e3;
    cover[k] = moist.cumulusCover[idx];
    buoyancy[k] = deep.deep ? moist.plumeBuoyancy[k] : NaN;
  }
  for (let k = 0; k <= K; k++) {
    massFlux[k] = !ran ? 0 : deep.deep ? deep.baseFlux * moist.cumulusFlux[k] : moist.cumulusFlux[k];
    downdraft[k] = deep.deep ? deep.downdraft * moist.draftFlux[k] : 0;
    if (deep.deep && k < K) plumeT[k] = moist.plumeInterfaceT[k];
  }
  return {
    type, parcelDeep: deep.parcelDeep, cellSize: Math.sqrt(model.mesh.areaCell[0]), dt, surfacePressure: ps,
    pressure, height, temperature, vapor, relativeHumidity: relative, mass, interfaces: Float64Array.from(levels, (s) => ps * s),
    top: ran ? moist.cumulusTop[0] : null,
    base: deep.deep ? ps * levels[deep.base] : ran ? moist.cumulus.lclPressure : null,
    cape: deep.deep ? deep.cape : null, inhibition: deep.deep ? deep.inhibition : null, tau: deep.deep ? deep.tau : null,
    baseFlux: moist.cumulusBaseFlux[0], downdraftFlux: deep.deep ? deep.downdraft : 0,
    rain: moist.rain[0] * perDay,
    convectiveRain: moist.convectivePrecipitation[0] * perDay,
    largeScaleRain: moist.largeScalePrecipitation[0] * perDay,
    snow: ice ? moist.convectiveSnow[0] * perDay : 0,
    heating, moistening, cloudWater, cover, massFlux, downdraft, plume: { buoyancy, interfaceTemperature: plumeT },
  };
}
