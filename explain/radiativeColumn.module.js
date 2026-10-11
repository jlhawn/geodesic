/*
 * One column of the model's own radiation on the bl36 layers over a slab
 * ground: createRadiation(...).column with all-zero humidity (the shipped
 * correlated-k and CLIRAD gases, never the gray q = null path) on cell 0 of
 * a 42-cell mesh, θ stepped by dt·Q/(c_p Δm Π) and the ground by dt·F/C,
 * then the model's dry convective adjustment. Layers run top first; fluxes
 * are W/m², heating K/day at constant pressure; every column closes
 * absorbed sunlight − OLR = Σ layer terms + the ground's net flux.
 */
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createSigmaCore, sigmaInterfaces, P0 } from '../js/dynamics/sigmaCore.module.js';
import { createRadiation, SOLAR_CONSTANT, STEFAN_BOLTZMANN, GREENHOUSE_GASES, SUMMED, CLEAR_SUMMED, DAY, YEAR } from '../js/physics/radiation.module.js';
import { createSurface, SEA_DRAG } from '../js/physics/surface.module.js';
import { ozoneWeights, ozoneAbove } from '../js/physics/ozone.module.js';
import { openWaterAlbedo } from '../js/physics/ice.module.js';
import { sunTrack } from './insolation.module.js';

export const LEVELS = sigmaInterfaces('bl36');
export const K = LEVELS.length - 1;
export const GLOBAL_MU = 2 / 3;
export const GLOBAL_BEAM = SOLAR_CONSTANT / 4;
export const HEAT_CAPACITY = 2.1e7;
export const SEA_DIFFUSE_ALBEDO = 0.06;
export const STANDARD_URL = new URL('./data/mls_bl36.json', import.meta.url).href;

const CORE_BUFFERS = ['flux', 'dissipation', 'divFlux', 'piSigmaDot', 'exnerLower', 'exnerLayer', 'dExnerDpi', 'thetaLower', 'qLower', 'qcLower', 'thetaV', 'geopotential', 'piVertex'];
const CELL_BUFFERS = ['outgoing', 'surfaceShortwave', 'evaporation', 'sensibleHeat', 'surfaceBuoyancy', 'stratus', 'stratusFraction', 'stabilityIndex', 'mlmCover', 'mlmWater', 'mlmEntrainment', 'mlmSubsidence', 'mlmHeight', 'mlmGate', 'mlmTop', 'stratiform'];
const GASES = { carbonDioxide: true, ozone: true, methane: true, nitrousOxide: true };

let world = null;
function setting() {
  if (world) return world;
  const mesh = buildMesh(new Grid(2)), C = mesh.nCells, E = mesh.nEdges;
  const sizes = { flux: K * E, dissipation: K * E, piSigmaDot: (K + 1) * C, piVertex: mesh.nVertices };
  const core = createSigmaCore(mesh, { levels: LEVELS, buffers: Object.fromEntries(CORE_BUFFERS.map((name) => [name, new ArrayBuffer(8 * (sizes[name] ?? K * C))])) });
  const surface = createSurface(mesh, core, { buffers: { windSpeed: new ArrayBuffer(8 * C) } });
  world = { mesh, core, surface, C, pi: new Float64Array(C), theta: new Float64Array(K * C) };
  return world;
}

function radiationBuffers(C) {
  const cell = () => new ArrayBuffer(8 * C);
  return {
    ...Object.fromEntries(CELL_BUFFERS.map((name) => [name, cell()])),
    summed: Object.fromEntries([...SUMMED, ...CLEAR_SUMMED].map((name) => [name, cell()])),
    longwave: new ArrayBuffer(8 * K * C),
  };
}

export function layerPressures(surfacePressure = P0) {
  return Float64Array.from({ length: K }, (_, k) => 0.5 * (LEVELS[k] + LEVELS[k + 1]) * surfacePressure);
}

export function interfacePressures(surfacePressure = P0) {
  return Float64Array.from(LEVELS, (s) => s * surfacePressure);
}

function ozoneLayers(weights, surfacePressure) {
  const out = new Float64Array(K);
  let above = 0;
  for (let k = 0; k < K; k++) {
    const below = ozoneAbove(surfacePressure * LEVELS[k + 1], weights);
    out[k] = below - above;
    above = below;
  }
  return out;
}

export function ozoneAt(latitude, day, surfacePressure = P0) {
  return ozoneLayers(ozoneWeights(latitude * Math.PI / 180, ((day * DAY) % YEAR) / YEAR), surfacePressure);
}

let sphere = null;
export function meanOzone(latitude = null, surfacePressure = P0) {
  if (latitude !== null) return ozoneLayers(yearWeights([latitude * Math.PI / 180]), surfacePressure);
  sphere ??= yearWeights(Array.from({ length: 3600 }, (_, b) => Math.asin(-1 + (b + 0.5) / 1800)));
  return ozoneLayers(sphere, surfacePressure);
}

function yearWeights(latitudes, months = 12) {
  const weights = new Float64Array(5), sample = new Float64Array(5), n = latitudes.length * months;
  for (const latitude of latitudes) for (let m = 0; m < months; m++) {
    ozoneWeights(latitude, m / months, sample);
    for (let r = 0; r < 5; r++) weights[r] += sample[r] / n;
  }
  return weights;
}

export function meanSunlight(latitude, { from = 0, days = 365, perDay = 288 } = {}) {
  const c = Math.cos(latitude * Math.PI / 180), s = Math.sin(latitude * Math.PI / 180), track = { x: new Float64Array(perDay), z: new Float64Array(perDay) };
  let shone = 0, weighted = 0, lit = 0;
  for (let d = 0; d < days; d++) {
    sunTrack(from + d, perDay, track);
    for (let k = 0; k < perDay; k++) {
      const mu = c * track.x[k] + s * track.z[k];
      if (mu > 0) { shone += mu; weighted += mu * mu; lit++; }
    }
  }
  const n = days * perDay;
  return { insolation: SOLAR_CONSTANT * shone / n, mu: shone > 0 ? weighted / shone : 0, daylight: lit / n };
}

export function readColumn(data) {
  if (data.levels.length !== LEVELS.length || data.levels.some((s, k) => Math.abs(s - LEVELS[k]) > 1e-12)) throw new Error(`${data.name} is not on the bl36 levels`);
  const array = (values) => Float64Array.from(values);
  return { ...data, temperature: array(data.temperature), ozone: array(data.ozone), humidity: array(data.humidity), pressure: array(data.pressure), height: array(data.height) };
}

const loading = new Map();
export function loadColumn(url = STANDARD_URL) {
  if (!loading.has(url)) loading.set(url, fetch(url).then((response) => {
    if (!response.ok) throw new Error(`${url}: ${response.status}`);
    return response.json();
  }).then(readColumn));
  return loading.get(url);
}

export function makeColumn(options = {}) {
  const { mesh, core, surface, C, pi } = setting(), { dSigma, cp, g, exnerLayer, geopotential } = core.diagnostics;
  const theta = new Float64Array(K * C), q = new Float64Array(K * C), before = new Float64Array(K);
  const settings = { surfacePressure: P0, albedo: 0.07, diffuseAlbedo: null, mu: GLOBAL_MU, beam: GLOBAL_BEAM, carbonDioxide: 1, gases: GASES, ozoneProfile: null, humidity: null, wind: 5, heatCapacity: HEAT_CAPACITY };
  const T = new Float64Array(K), potential = new Float64Array(K), height = new Float64Array(K), mass = new Float64Array(K);
  const field = () => new Float64Array(K);
  const terms = { sunlight: field(), infrared: field(), ground: field(), convection: field() };
  const heating = { sunlight: field(), infrared: field(), ground: field(), convection: field(), total: field() };
  const budget = {};
  const column = { T, theta: potential, height, pressure: field(), interfacePressure: new Float64Array(K + 1), terms, heating, budget, Ts: 288, time: 0, settings, set, step, evaluate, restart, coldPoint };
  let radiation = null, ozone = null, sensible = true, start = null;

  function albedos() {
    const mu = settings.beam > 0 ? settings.mu : 1;
    if (settings.albedo === 'sea') return [openWaterAlbedo(mu), settings.diffuseAlbedo ?? SEA_DIFFUSE_ALBEDO];
    return [settings.albedo, settings.diffuseAlbedo ?? settings.albedo];
  }

  function build() {
    const { surfacePressure, mu, beam, carbonDioxide, gases: on } = settings;
    ozone = Float64Array.from(settings.ozoneProfile ?? meanOzone(null, surfacePressure), (x) => (on.ozone ? x : 0));
    radiation = createRadiation(mesh, core, {
      solarConstant: beam > 0 ? beam / mu : SOLAR_CONSTANT,
      carbonDioxide: on.carbonDioxide ? carbonDioxide * GREENHOUSE_GASES.carbonDioxide : 0,
      methane: on.methane ? GREENHOUSE_GASES.methane : 0,
      nitrousOxide: on.nitrousOxide ? GREENHOUSE_GASES.nitrousOxide : 0,
      ozoneProfile: ozone,
      buffers: radiationBuffers(C),
    });
    pi[0] = surfacePressure;
    column.pressure.set(layerPressures(surfacePressure));
    column.interfacePressure.set(interfacePressures(surfacePressure));
    for (let k = 0; k < K; k++) mass[k] = surfacePressure * dSigma[k] / g;
    q.fill(0);
    if (settings.humidity) for (let k = 0; k < K; k++) q[k * C] = settings.humidity[k];
  }

  function load(temperature) {
    pi[0] = settings.surfacePressure;
    core.diagnoseColumn(0, pi, theta);
    for (let k = 0; k < K; k++) theta[k * C] = (typeof temperature === 'number' ? temperature : temperature[k]) / exnerLayer[k * C];
  }

  function refresh() {
    pi[0] = settings.surfacePressure;
    core.diagnoseColumn(0, pi, theta);
    for (let k = 0; k < K; k++) {
      potential[k] = theta[k * C];
      T[k] = theta[k * C] * exnerLayer[k * C];
      height[k] = geopotential[k * C] / g;
    }
  }

  function radiate(part) {
    refresh();
    const [direct, diffuse] = albedos(), { beam } = settings, Ts = column.Ts;
    const net = radiation.column(0, settings.surfacePressure, theta, Ts, settings.wind, undefined, beam, null, q, null, direct, diffuse, 1, sensible ? SEA_DRAG : 0, 0, 0, 0, 0, part);
    const b = radiation.budget, flux = radiation.layerFlux, longwave = radiation.longwave;
    for (let k = 0; k < K; k++) {
      terms.infrared[k] = longwave[k * C];
      terms.ground[k] = k === K - 1 ? b.sensibleHeat : 0;
      terms.sunlight[k] = flux[k] - terms.infrared[k] - terms.ground[k];
    }
    const emission = STEFAN_BOLTZMANN * Ts ** 4;
    Object.assign(budget, {
      insolation: beam, reflected: b.reflectedSolar, albedo: beam > 0 ? b.reflectedSolar / beam : 0, absorbed: b.absorbedSolar,
      air: b.atmosphereSolar, ozone: b.ozoneSolar, oxygen: b.oxygenSolar, carbonDioxide: b.carbonDioxideSolar, vapor: b.vaporSolar, aerosol: b.aerosolSolar,
      upwardNearInfrared: b.atmosphereSolar - b.ozoneSolar - b.oxygenSolar - b.carbonDioxideSolar - b.vaporSolar - b.aerosolSolar - b.cloudSolar,
      ground: b.surfaceSolar, reachingGround: b.surfaceShortwave, directBeam: b.surfaceDirect,
      olr: b.outgoingLongwave, dlr: b.downwardLongwave, groundEmission: emission, greenhouse: emission - b.outgoingLongwave, escaping: part === 'radiation' ? b.escape : NaN,
      sensible: b.sensibleHeat, groundNet: net, imbalance: b.absorbedSolar - b.outgoingLongwave,
    });
    return net;
  }

  function toHeating() {
    for (let k = 0; k < K; k++) {
      const perDay = DAY / (cp * mass[k]);
      heating.sunlight[k] = terms.sunlight[k] * perDay;
      heating.infrared[k] = terms.infrared[k] * perDay;
      heating.ground[k] = terms.ground[k] * perDay;
      heating.convection[k] = terms.convection[k] * perDay;
      heating.total[k] = heating.sunlight[k] + heating.infrared[k] + heating.ground[k] + heating.convection[k];
    }
  }

  function evaluate(choice = {}) {
    if (choice.sensible !== undefined) sensible = choice.sensible;
    radiate('radiation');
    terms.convection.fill(0);
    toHeating();
    return column;
  }

  function step(dt, { convection = true, sensible: groundHeat = convection } = {}) {
    sensible = groundHeat;
    const net = radiate('all'), flux = radiation.layerFlux;
    for (let k = 0; k < K; k++) theta[k * C] += dt * flux[k] / (cp * mass[k] * exnerLayer[k * C]);
    column.Ts += dt * net / settings.heatCapacity;
    if (convection) {
      for (let k = 0; k < K; k++) before[k] = theta[k * C];
      surface.convectiveAdjustColumn(0, pi, theta);
      for (let k = 0; k < K; k++) terms.convection[k] = cp * mass[k] * exnerLayer[k * C] * (theta[k * C] - before[k]) / dt;
    } else terms.convection.fill(0);
    column.time += dt;
    refresh();
    toHeating();
    return column;
  }

  function set(changes = {}) {
    const { profile, temperature, groundTemperature, gases, ...rest } = changes;
    if (profile) {
      settings.surfacePressure = profile.surfacePressure;
      settings.ozoneProfile = profile.ozone;
    }
    Object.assign(settings, rest);
    settings.gases = { ...settings.gases, ...gases };
    build();
    const startT = temperature ?? profile?.temperature;
    if (startT !== undefined) load(startT);
    const startTs = groundTemperature ?? profile?.groundTemperature;
    if (startTs !== undefined) column.Ts = startTs;
    return evaluate();
  }

  function restart() {
    load(start.T);
    column.Ts = start.Ts;
    column.time = 0;
    return evaluate();
  }

  function coldPoint(above = 500) {
    let best = -1;
    for (let k = 0; k < K; k++) if (column.pressure[k] > above && (best < 0 || T[k] < T[best])) best = k;
    return { k: best, pressure: column.pressure[best], temperature: T[best], height: height[best] };
  }

  const { profile } = options;
  set({ ...options, temperature: options.temperature ?? profile?.temperature ?? 288, groundTemperature: options.groundTemperature ?? profile?.groundTemperature ?? 288 });
  start = { T: Float64Array.from(T), Ts: column.Ts };
  return column;
}
