// Precomputes the layered-model frames that explain/sphere.html,
// explain/primitive.html and explain/radiation.html play back:
//   node scripts/explainLayers.mjs jw06|heldsuarez|radiation|tracks|state [OUT]
//   node scripts/explainLayers.mjs wet [SPINUP_DAYS] [RECORD_DAYS] [OUT]
// writing explain/data/<case>_N16.bin, jw06state_N16.bin for state and
// water_N16.bin for wet (365 days of spin-up and 365 recorded by default).
import { writeFileSync, mkdirSync } from 'node:fs';
import { dirname } from 'node:path';
import { execFileSync } from 'node:child_process';
import { Grid } from '../js/grid.module.js';
import { buildMesh } from '../js/mesh.module.js';
import { createSigmaCore, sigmaInterfaces } from '../js/dynamics/sigmaCore.module.js';
import { createRK4Arrays } from '../js/dynamics/integrators.module.js';
import { cellVector } from '../js/dynamics/operators.module.js';
import { createRadiation, YEAR } from '../js/physics/radiation.module.js';
import { createSurface } from '../js/physics/surface.module.js';
import { saturationHumidity } from '../js/physics/moist.module.js';
import { FREEZING_POINT, MELTING_POINT } from '../js/physics/ice.module.js';
import { initializeState } from '../js/physics/init.module.js';
import { createParallelModel } from '../js/parallel.module.js';
import { withCadence } from '../js/cadence.module.js';
import { levelFields } from '../js/levels.module.js';
import { EARTH, edgeNormalVelocity, hyperdiffusion, cellVelocity } from '../explain/swCases.module.js';
import { PUBLISHED_RADIATION, PUBLISHED_MOIST } from '../explain/published.module.js';

const which = process.argv[2], N = 16, DAY = 86400, dt = 450 * 16 / N;
const named = { state: 'jw06state', wet: 'water' };
const out = (which === 'wet' ? process.argv[5] : process.argv[3]) ?? new URL(`../explain/data/${named[which] ?? which}_N${N}.bin`, import.meta.url).pathname;
const { a, omega, g } = EARTH, R = 287, cp = 1004.5, p0 = 1e5;
const mesh = buildMesh(new Grid(N), { radius: a, omega }), C = mesh.nCells, E = mesh.nEdges;
const EQUINOX_DAY_OF_YEAR = 78;
const MONTHS = [['January', 31], ['February', 28], ['March', 31], ['April', 30], ['May', 31], ['June', 30], ['July', 31], ['August', 31], ['September', 30], ['October', 31], ['November', 30], ['December', 31]];
const monthStart = MONTHS.map((_, m) => MONTHS.slice(0, m).reduce((sum, [, days]) => sum + days, 0));
const monthOf = (day) => { const date = (Math.floor(day) + EQUINOX_DAY_OF_YEAR) % 365; let m = 11; while (monthStart[m] > date) m--; return m; };
const monthDays = (m) => { const from = (monthStart[m] - EQUINOX_DAY_OF_YEAR + 365) % 365, to = from + MONTHS[m][1]; return to <= 365 ? [[from, to]] : [[from, 365], [0, to - 365]]; };

function jw06() {
  const eta0 = 0.252, u0 = 35, T0 = 288, gamma = 0.005, deltaT = 4.8e5, etaT = 0.2;
  const etaV = (eta) => (eta - eta0) * Math.PI / 2;
  const shape = (lat) => { const s = Math.sin(lat), c = Math.cos(lat); return { A: -2 * s ** 6 * (c * c + 1 / 3) + 10 / 63, B: 8 / 5 * c ** 3 * (s * s + 2 / 3) - Math.PI / 4 }; };
  const zonalWind = (eta, lat) => u0 * Math.cos(etaV(eta)) ** 1.5 * Math.sin(2 * lat) ** 2;
  const geopotential = (eta, lat) => {
    const exp = R * gamma / g;
    let mean = T0 * g / gamma * (1 - eta ** exp);
    if (eta < etaT) mean -= R * deltaT * ((Math.log(eta / etaT) + 137 / 60) * etaT ** 5 - 5 * etaT ** 4 * eta + 5 * etaT ** 3 * eta ** 2 - (10 / 3) * etaT ** 2 * eta ** 3 + (5 / 4) * etaT * eta ** 4 - eta ** 5 / 5);
    const cv = Math.cos(etaV(eta)) ** 1.5, { A, B } = shape(lat);
    return mean + u0 * cv * (A * u0 * cv + B * a * omega);
  };
  const perturbation = (lon, lat) => {
    const lonC = Math.PI / 9, latC = 2 * Math.PI / 9, radius = a / 10;
    const r = a * Math.acos(Math.min(1, Math.sin(latC) * Math.sin(lat) + Math.cos(latC) * Math.cos(lat) * Math.cos(lon - lonC)));
    return Math.exp(-((r / radius) ** 2));
  };
  const surface = Float64Array.from(mesh.latCell, (lat) => geopotential(1, lat));
  const core = createSigmaCore(mesh, { g, cp, R, p0, surfaceGeopotential: surface, nu4: hyperdiffusion(mesh, 3), nu4Theta: hyperdiffusion(mesh, 3) });
  const K = core.K, pi = new Float64Array(C).fill(p0), theta = new Float64Array(K * C);
  core.diagnose(pi, theta);
  const { exnerLayer, exnerLower } = core.arrays;
  for (let i = 0; i < C; i++) {
    let below = surface[i], thetaBelow = 0;
    for (let k = K - 1; k >= 0; k--) {
      const idx = k * C + i, target = geopotential(core.sigmaMid[k], mesh.latCell[i]);
      const fromBelow = k === K - 1 ? 0 : cp * thetaBelow * (exnerLayer[idx + C] - exnerLower[idx]);
      theta[idx] = (target - below - fromBelow) / (cp * (exnerLower[idx] - exnerLayer[idx]));
      below = target; thetaBelow = theta[idx];
    }
  }
  const u = new Float64Array(K * E);
  for (let k = 0; k < K; k++) u.set(edgeNormalVelocity(mesh, (lon, lat) => ({ zonal: zonalWind(core.sigmaMid[k], lat) + perturbation(lon, lat), meridional: 0 })), k * E);
  return { core, state: [pi, theta, u], tendency: core.tendency, surface };
}

function heldSuarez() {
  const core = createSigmaCore(mesh, { g, cp, R, p0, nu4: hyperdiffusion(mesh, 3), nu4Theta: hyperdiffusion(mesh, 3) });
  const { K, sigmaMid, kappa, exnerLayer } = core.diagnostics;
  const kf = 1 / DAY, ka = 1 / (40 * DAY), ks = 1 / (4 * DAY), sigmaB = 0.7, deltaTy = 60, deltaThetaZ = 10, tMin = 200, tMax = 315;
  const cosLat = Float64Array.from(mesh.latCell, Math.cos), sinLat = Float64Array.from(mesh.latCell, Math.sin);
  const nudging = (state, outs) => {
    const [pi, theta] = state, [, dTheta] = outs;
    for (let k = 0; k < K; k++) {
      const s = sigmaMid[k], weight = Math.max(0, (s - sigmaB) / (1 - sigmaB));
      for (let i = 0; i < C; i++) {
        const idx = k * C + i, p = pi[i] * s, c2 = cosLat[i] * cosLat[i];
        const tEq = Math.max(tMin, (tMax - deltaTy * sinLat[i] * sinLat[i] - deltaThetaZ * Math.log(p / p0) * c2) * Math.pow(p / p0, kappa));
        dTheta[idx] -= (ka + (ks - ka) * weight * c2 * c2) * (theta[idx] - tEq / exnerLayer[idx]);
      }
    }
  };
  const friction = (state, outs) => {
    const u = state[2], dU = outs[2];
    for (let k = 0; k < K; k++) {
      const weight = Math.max(0, (sigmaMid[k] - sigmaB) / (1 - sigmaB));
      if (weight > 0) for (let e = 0; e < E; e++) dU[k * E + e] -= kf * weight * u[k * E + e];
    }
  };
  const tendency = (state, outs) => { core.tendency(state, outs); nudging(state, outs); friction(state, outs); };
  const frictionOnly = (state, outs) => { core.tendency(state, outs); friction(state, outs); };
  const pi = new Float64Array(C).fill(p0), theta = new Float64Array(K * C);
  core.diagnose(pi, theta);
  let seed = 12345;
  const noise = () => { seed = (seed * 1664525 + 1013904223) >>> 0; return seed / 2 ** 32 - 0.5; };
  for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) {
    const p = pi[i] * sigmaMid[k], c2 = cosLat[i] ** 2;
    const tEq = Math.max(tMin, (tMax - deltaTy * sinLat[i] ** 2 - deltaThetaZ * Math.log(p / p0) * c2) * (p / p0) ** kappa);
    theta[k * C + i] = tEq / exnerLayer[k * C + i] + 0.1 * noise();
  }
  return { core, state: [pi, theta, new Float64Array(K * E)], tendency, frictionOnly, surface: null };
}

if (which === 'wet') {
  const SPINUP_DAYS = Number(process.argv[3] ?? 365), RECORD_DAYS = Number(process.argv[4] ?? 365);
  if (![SPINUP_DAYS, RECORD_DAYS].every((days) => Number.isInteger(days) && days >= 0)) { console.error('usage: node scripts/explainLayers.mjs wet [SPINUP_DAYS] [RECORD_DAYS] [OUT]'); process.exit(1); }
  const STEP = 450, WORKERS = 6, GRID = 'bl36', FRAME_DAYS = 5, BINS = 36, ICE = JSON.parse(process.env.WET_ICE ?? '{}');
  const weatherFrom = Math.max(0, Math.min(287, RECORD_DAYS - 30)), WEATHER_DAYS = [weatherFrom, Math.min(RECORD_DAYS, weatherFrom + 30)];
  const cadence = withCadence({ ocean: false }, STEP, N);
  const options = { ocean: false, levels: sigmaInterfaces(GRID), ...cadence, radiation: { ...cadence.radiation, clearSkyPass: true, orbit: PUBLISHED_RADIATION.orbit, shortwave: PUBLISHED_RADIATION.shortwave, aerosol: PUBLISHED_RADIATION.aerosol }, moist: { bergeronAutoconversion: PUBLISHED_MOIST.bergeronAutoconversion }, ice: ICE };
  const t0 = performance.now(), seconds = () => (performance.now() - t0) / 1000;
  const model = await createParallelModel(new Grid(N), options, WORKERS);
  initializeState(model).forEach((values, s) => model.state[s].set(values));
  const startTime = (Math.ceil(SPINUP_DAYS / 365) * 365 - SPINUP_DAYS) * DAY, recordStart = startTime + SPINUP_DAYS * DAY;
  model.time = startTime;
  model.restartPrecipitation();

  const { mesh: grid, core, radiation, moist, seaIce } = model, [pi, theta, u, surfaceT, q, qc, ice] = model.state;
  const cells = grid.nCells, edges = grid.nEdges, K = core.K, { sigmaMid, dSigma } = core.diagnostics, gravity = core.diagnostics.g, airCp = core.diagnostics.cp;
  const look = createSigmaCore(grid, { levels: core.levels }), { exnerLayer: lookExner } = look.diagnostics;
  const area = grid.areaCell.reduce((sum, x) => sum + x, 0), weight = grid.areaCell;
  const bin = Int32Array.from(grid.latCell, (lat) => Math.min(BINS - 1, Math.floor((lat + Math.PI / 2) / Math.PI * BINS)));
  const binArea = new Float64Array(BINS);
  for (let i = 0; i < cells; i++) binArea[bin[i]] += weight[i];
  const east = new Float64Array(3 * cells), north = new Float64Array(3 * cells);
  for (let i = 0; i < cells; i++) {
    const x = grid.xCell[3 * i], y = grid.xCell[3 * i + 1], z = grid.xCell[3 * i + 2], rho = Math.hypot(x, y) || 1e-12;
    east.set([-y / rho, x / rho, 0], 3 * i); north.set([-z * x / rho, -z * y / rho, rho], 3 * i);
  }
  const columnMass = (i, k) => pi[i] * dSigma[k] / gravity;
  const zonal = (value, digits = 3) => {
    const sum = new Float64Array(BINS);
    for (let i = 0; i < cells; i++) sum[bin[i]] += weight[i] * value(i);
    return Array.from(sum, (v, b) => (digits === null ? v / binArea[b] : +(v / binArea[b]).toFixed(digits)));
  };
  const global = (value) => { let sum = 0; for (let i = 0; i < cells; i++) sum += weight[i] * value(i); return sum / area; };

  const perDay = Math.round(DAY / STEP), steps = (SPINUP_DAYS + RECORD_DAYS) * perDay, recordStep = SPINUP_DAYS * perDay, sixHours = Math.round(6 * 3600 / STEP), end = RECORD_DAYS * DAY;
  const FIELDS = [
    { name: 'surfaceTemperature', unit: 'K', per: 50, offset: 260 }, { name: 'outgoingLongwave', unit: 'W/m²', per: 50, offset: 0 },
    { name: 'temperature850', unit: 'K', per: 50, offset: 260 }, { name: 'wind250', unit: 'm/s', per: 100, offset: 0 },
    { name: 'rain', unit: 'mm/d', per: 100, offset: 0 }, { name: 'precipitableWater', unit: 'kg/m²', per: 200, offset: 0 }, { name: 'cloudWater', unit: 'g/m²', per: 20, offset: 0 },
  ];
  const pack = (values, { per, offset }) => Int16Array.from(values, (v) => Math.max(-32768, Math.min(32767, Math.round((v - offset) * per))));
  const rainBank = new Float64Array(cells);
  const sequences = [{ name: 'year', every: FRAME_DAYS * DAY, from: 0, to: end }, { name: 'weather', every: 6 * 3600, from: WEATHER_DAYS[0] * DAY, to: WEATHER_DAYS[1] * DAY }]
    .map((s) => ({ ...s, stride: Math.round(s.every / STEP), times: [], frames: [], last: new Float64Array(cells), lastTime: -SPINUP_DAYS * DAY }));
  function wetFrame(s, time) {
    look.diagnose(pi, theta);
    const winds = new Array(K).fill(null), layerWind = (k) => (winds[k] ??= cellVector(grid, u.subarray(k * edges, (k + 1) * edges), new Float64Array(3 * cells)));
    const t850 = levelFields(look, pi, theta, layerWind, 850, q).temperature, w250 = levelFields(look, pi, theta, layerWind, 250, q).speed;
    const values = [
      surfaceT, radiation.outgoing, t850, w250,
      Float64Array.from({ length: cells }, (_, i) => (rainBank[i] + moist.precipitation[i] - s.last[i]) * DAY / Math.max(STEP, time - s.lastTime)),
      Float64Array.from({ length: cells }, (_, i) => moist.columnWater(pi, q, i)),
      Float64Array.from({ length: cells }, (_, i) => 1000 * model.cloudWater(i)),
    ];
    return values.map((v, f) => pack(v, FIELDS[f]));
  }

  const SECTION = ['T', 'u', 'v', 'q', 'RH', 'qc'];
  const sectionSums = () => ({ ...Object.fromEntries(SECTION.map((name) => [name, new Float64Array(BINS * K)])), surfacePressure: new Float64Array(BINS), samples: 0 });
  const sections = MONTHS.map(sectionSums), wind = new Float64Array(3 * cells);
  function sampleSection(sums) {
    look.diagnose(pi, theta);
    for (let i = 0; i < cells; i++) sums.surfacePressure[bin[i]] += weight[i] * pi[i];
    for (let k = 0; k < K; k++) {
      cellVector(grid, u.subarray(k * edges, (k + 1) * edges), wind);
      for (let i = 0; i < cells; i++) {
        const n = k * cells + i, at = bin[i] * K + k, w = weight[i], T = theta[n] * lookExner[n];
        const ue = wind[3 * i] * east[3 * i] + wind[3 * i + 1] * east[3 * i + 1] + wind[3 * i + 2] * east[3 * i + 2];
        const vn = wind[3 * i] * north[3 * i] + wind[3 * i + 1] * north[3 * i + 1] + wind[3 * i + 2] * north[3 * i + 2];
        sums.T[at] += w * T; sums.u[at] += w * ue; sums.v[at] += w * vn; sums.q[at] += w * q[n]; sums.qc[at] += w * qc[n];
        sums.RH[at] += w * q[n] / saturationHumidity(T, pi[i] * sigmaMid[k]);
      }
    }
    sums.samples++;
  }
  const SECTION_DIGITS = { T: [2, 'fixed'], u: [2, 'fixed'], v: [3, 'fixed'], q: [4, 'precision'], RH: [3, 'fixed'], qc: [3, 'precision'] };
  const SECTION_SCALE = { T: 1, u: 1, v: 1, q: 1000, RH: 1, qc: 1e6 };
  const finishSection = (sums) => ({
    samples: sums.samples,
    ...Object.fromEntries(SECTION.map((name) => {
      const [digits, kind] = SECTION_DIGITS[name];
      return [name, Array.from(sums[name], (v, at) => { const x = SECTION_SCALE[name] * v / (binArea[Math.floor(at / K)] * sums.samples); return kind === 'fixed' ? +x.toFixed(digits) : +x.toPrecision(digits); })];
    })),
    surfacePressure: Array.from(sums.surfacePressure, (v, b) => +(v / (binArea[b] * sums.samples) / 100).toFixed(2)),
  });

  const BUDGET = ['absorbedSolar', 'atmosphereSolar', 'outgoingLongwave', 'clearAbsorbedSolar', 'clearOutgoingLongwave', 'reflectedSolar', 'insolation', 'surfaceTemperature', 'precipitation', 'convectiveRain', 'largeScaleRain', 'evaporation', 'sensibleHeat', 'precipitableWater', 'cloudCondensate', 'cloudWater', 'ice'];
  const budgetSums = () => ({ ...Object.fromEntries(BUDGET.map((name) => [name, new Float64Array(BINS)])), days: 0 });
  const budgets = MONTHS.map(budgetSums);
  const today = Object.fromEntries(['evaporation', 'sensibleHeat', 'surfaceTemperature', 'precipitableWater', 'cloudCondensate', 'cloudWater', 'ice'].map((name) => [name, new Float64Array(cells)]));
  const cellDay = Object.fromEntries(BUDGET.map((name) => [name, new Float64Array(cells)]));
  const finishBudget = (sums) => {
    const mean = (name, digits) => Array.from(sums[name], (v, b) => +(v / (binArea[b] * sums.days)).toFixed(digits));
    return {
      days: sums.days,
      ...Object.fromEntries(['absorbedSolar', 'atmosphereSolar', 'outgoingLongwave', 'clearAbsorbedSolar', 'clearOutgoingLongwave', 'surfaceTemperature', 'sensibleHeat', 'precipitableWater', 'cloudCondensate', 'cloudWater'].map((name) => [name, mean(name, 2)])),
      ...Object.fromEntries(['precipitation', 'convectiveRain', 'largeScaleRain', 'evaporation', 'ice'].map((name) => [name, mean(name, 4)])),
      planetaryAlbedo: Array.from(sums.reflectedSolar, (v, b) => (sums.insolation[b] > 0 ? +(v / sums.insolation[b]).toFixed(4) : null)),
    };
  };
  const airEnergy = (i) => { let sum = 0; for (let k = 0; k < K; k++) sum += airCp * theta[k * cells + i] * lookExner[k * cells + i] * columnMass(i, k); return sum; };
  const storage = () => {
    look.diagnose(pi, theta);
    return {
      time: model.time - recordStart,
      surfaceEnergy: zonal((i) => seaIce.energy(surfaceT[i], ice[i], seaIce.snow[i], seaIce.cover(i, ice[i])), null),
      airEnergy: zonal(airEnergy, null),
      precipitableWater: zonal((i) => moist.columnWater(pi, q, i), null), cloudCondensate: zonal((i) => moist.columnWater(pi, qc, i), null),
      surfaceTemperature: zonal((i) => surfaceT[i], null), iceVolume: zonal((i) => seaIce.cover(i, ice[i]) * ice[i], null),
    };
  };

  const DAILY = ['surfaceTemperature', 'absorbedSolar', 'outgoingLongwave', 'clearAbsorbedSolar', 'clearOutgoingLongwave', 'planetaryAlbedo', 'precipitation', 'convectiveRain', 'largeScaleRain', 'evaporation', 'sensibleHeat', 'precipitableWater', 'cloudCondensate', 'cloudWater', 'ice', 'maxWind', 'waterResidual'];
  const daily = Object.fromEntries(DAILY.map((name) => [name, []]));
  const totalWater = () => global((i) => moist.columnWater(pi, q, i) + moist.columnWater(pi, qc, i));
  const water0 = totalWater();
  let evaporated = 0, rained = 0, storedStart = null, storedEnd = null, lastDay = performance.now();

  for (let n = 0; n <= steps; n++) {
    const time = (n - recordStep) * STEP;
    for (const s of sequences) {
      if ((n - recordStep) % s.stride !== 0) continue;
      if (time >= s.from && time <= s.to) { s.times.push(time); s.frames.push(wetFrame(s, time)); }
      for (let i = 0; i < cells; i++) s.last[i] = rainBank[i] + moist.precipitation[i];
      s.lastTime = time;
    }
    if (time >= 0 && time < end && (n - recordStep) % sixHours === 0) sampleSection(sections[monthOf(time / DAY)]);
    if (n === recordStep) storedStart = storage();
    if (n === steps) { storedEnd = storage(); break; }
    model.step(STEP);
    for (let i = 0; i < cells; i++) {
      today.evaporation[i] += radiation.evaporation[i] * STEP; today.sensibleHeat[i] += radiation.sensibleHeat[i] / perDay; today.surfaceTemperature[i] += surfaceT[i] / perDay;
      today.precipitableWater[i] += moist.columnWater(pi, q, i) / perDay; today.cloudCondensate[i] += 1000 * moist.columnWater(pi, qc, i) / perDay;
      today.cloudWater[i] += 1000 * model.cloudWater(i) / perDay;
      today.ice[i] += seaIce.cover(i, ice[i]) / perDay;
    }
    if ((n + 1) % perDay !== 0) continue;

    const meanSteps = model.radiationSteps, summed = radiation.summed;
    for (let i = 0; i < cells; i++) {
      for (const name of ['absorbedSolar', 'atmosphereSolar', 'outgoingLongwave', 'clearAbsorbedSolar', 'clearOutgoingLongwave']) cellDay[name][i] = summed[name][i] / meanSteps;
      cellDay.reflectedSolar[i] = summed.reflectedSolar[i]; cellDay.insolation[i] = summed.insolation[i];
      cellDay.precipitation[i] = moist.precipitation[i]; cellDay.convectiveRain[i] = moist.convectivePrecipitation[i]; cellDay.largeScaleRain[i] = moist.largeScalePrecipitation[i];
      for (const name of Object.keys(today)) cellDay[name][i] = today[name][i];
      rainBank[i] += moist.precipitation[i];
    }
    const d = model.diagnostics();
    for (const array of Object.values(today)) array.fill(0);
    evaporated += global((i) => cellDay.evaporation[i]); rained += global((i) => cellDay.precipitation[i]);
    const residual = totalWater() - water0 - evaporated + rained;
    const dayStart = (n + 1 - perDay - recordStep) * STEP;
    if (dayStart >= 0) {
      const sums = budgets[monthOf(dayStart / DAY)];
      for (let i = 0; i < cells; i++) for (const name of BUDGET) sums[name][bin[i]] += weight[i] * cellDay[name][i];
      sums.days++;
    }
    const means = Object.fromEntries(['surfaceTemperature', 'precipitation', 'convectiveRain', 'largeScaleRain', 'evaporation', 'sensibleHeat', 'precipitableWater', 'cloudCondensate', 'cloudWater', 'ice'].map((name) => [name, global((i) => cellDay[name][i])]));
    const values = {
      ...means, absorbedSolar: d.absorbedSolar, outgoingLongwave: d.outgoingLongwave, clearAbsorbedSolar: d.clearAbsorbedSolar, clearOutgoingLongwave: d.clearOutgoingLongwave, planetaryAlbedo: d.planetaryAlbedo,
      ice: 100 * means.ice, maxWind: d.maxWind, waterResidual: residual,
    };
    for (const name of DAILY) daily[name].push(Number.isFinite(values[name]) ? +values[name].toPrecision(6) : null);
    let low = Infinity, high = -Infinity, northIce = 0, southIce = 0, northArea = 0, southArea = 0;
    for (let i = 0; i < cells; i++) {
      low = Math.min(low, surfaceT[i]); high = Math.max(high, surfaceT[i]);
      const cover = seaIce.cover(i, ice[i]) * weight[i];
      if (grid.latCell[i] >= 0) { northIce += cover; northArea += weight[i]; } else { southIce += cover; southArea += weight[i]; }
    }
    let fastest = 0;
    for (let m = 1; m < K * edges; m++) if (Math.abs(u[m]) > Math.abs(u[fastest])) fastest = m;
    const fastLayer = Math.floor(fastest / edges), fastEdge = fastest % edges, fastCell = grid.cellsOnEdge[2 * fastEdge];
    const day = (n + 1) / perDay, pace = (performance.now() - lastDay) / 1000;
    lastDay = performance.now();
    console.error(`wet day ${day} of ${steps / perDay}${dayStart < 0 ? ' (spin-up)' : ''}: Ts ${values.surfaceTemperature.toFixed(2)} K (${low.toFixed(1)}–${high.toFixed(1)}), ASR ${d.absorbedSolar.toFixed(1)} OLR ${d.outgoingLongwave.toFixed(1)} clear ${d.clearAbsorbedSolar.toFixed(1)}/${d.clearOutgoingLongwave.toFixed(1)} W/m², albedo ${d.planetaryAlbedo.toFixed(3)}, P ${means.precipitation.toFixed(2)} (conv ${means.convectiveRain.toFixed(2)}, ls ${means.largeScaleRain.toFixed(2)}) E ${means.evaporation.toFixed(2)} mm/d, SH ${means.sensibleHeat.toFixed(1)} W/m², water ${means.precipitableWater.toFixed(1)} kg/m², cloud ${means.cloudCondensate.toFixed(0)} g/m² (${means.cloudWater.toFixed(0)} with cumulus and deck), ice ${(100 * d.iceFraction).toFixed(2)} % (north ${(100 * northIce / northArea).toFixed(1)}, south ${(100 * southIce / southArea).toFixed(1)}), max |u| ${d.maxWind.toFixed(0)} m/s at ${(sigmaMid[fastLayer] * pi[fastCell] / 100).toPrecision(3)} hPa ${(grid.latEdge[fastEdge] * 180 / Math.PI).toFixed(0)}°, lost ${residual.toExponential(2)} kg/m², ${pace.toFixed(1)} s/day, ${seconds().toFixed(0)} s`);
    if (![values.surfaceTemperature, d.absorbedSolar, d.outgoingLongwave, means.precipitableWater, d.maxWind, low, high].every(Number.isFinite) || d.maxWind > 400 || low < 100 || high > 400) { console.error('unstable: stopping'); await model.close(); process.exit(1); }
  }

  const sumMonths = (list, make) => {
    const total = make();
    for (const sums of list) for (const [name, value] of Object.entries(sums)) { if (typeof value === 'number') total[name] += value; else value.forEach((v, m) => { total[name][m] += v; }); }
    return total;
  };
  const annualSection = sumMonths(sections, sectionSums), annualBudget = sumMonths(budgets, budgetSums);
  const recorded = (list, finish) => MONTHS.map(([name], m) => ({ name, days: monthDays(m), ...finish(list[m]) }));
  const git = (...args) => { try { return execFileSync('git', args, { cwd: new URL('..', import.meta.url).pathname, encoding: 'utf8' }).trim(); } catch { return null; } };
  const probe = (h, snow, T) => +seaIce.albedo(h, null, snow, h > 0 ? 1 : 0, T).toFixed(3);
  const header = {
    version: 1, case: 'wet', N, K, dt: STEP, levels: Array.from(core.levels), grid: GRID, workers: WORKERS,
    fields: FIELDS.map(({ name, unit, per, offset }) => ({ name, unit, scale: 1 / per, offset })),
    sequences: sequences.map((s) => ({ name: s.name, every: s.every, times: s.times })),
    time: 'seconds from the March equinox (March 20) that begins the recorded year; the model clock reads recordStart + time and the sun is radiation.setTime(model.time)',
    rain: 'each frame\'s rain is the mean over the interval since the sequence\'s previous frame time (5 days for year, 6 hours for weather), or since the run began when that is shorter',
    spinupDays: SPINUP_DAYS, recordDays: RECORD_DAYS, startTime, recordStart,
    options: { ...options, levels: GRID }, cadence: { radiationEvery: options.radiation.radiationEvery, dragEvery: options.dragEvery },
    initial: 'initializeState(model): sea at 30 cos²φ − 2 °C (28 °C at the equator) less 4 K over 50–70°, 1.5 m of ice poleward of 72°N and 0.7 m poleward of 68°S, RH 0.7σ², winds geostrophic outside ±15°',
    ice: {
      options: ICE, slabHeatCapacity: seaIce.slabHeatCapacity, latentHeatFusion: seaIce.latentHeatFusion, freezingPoint: FREEZING_POINT, meltingPoint: MELTING_POINT,
      diffuseAlbedo: { water: probe(0, 0, 0), coldIce: probe(2, 0, 250), meltingIce: probe(2, 0, MELTING_POINT), freshSnowOnIce: probe(2, 100, 250) },
    },
    constants: { cp: airCp, R: core.diagnostics.R, g: gravity, latentHeat: moist.latentHeat, radius: Math.sqrt(area / (4 * Math.PI)) },
    commit: git('rev-parse', '--short', 'HEAD'), modified: (git('status', '--porcelain', '--', 'js', 'scripts') ?? '').split('\n').filter(Boolean),
    zonalMean: { bins: BINS, from: 0, to: end, sigma: Array.from(sigmaMid), weighting: 'area, 6-hourly samples; q in g/kg, qc in mg/kg, RH over liquid water, surfacePressure in hPa', ...finishSection(annualSection), months: recorded(sections, finishSection) },
    zonalBudget: { bins: BINS, from: 0, to: end, units: 'W/m², K, mm/d for precipitation and evaporation, kg/m² for precipitableWater, g/m² for clouds, area fraction for ice; cloudWater adds the cumulus and the deck to the resolved condensate', ...finishBudget(annualBudget), months: recorded(budgets, finishBudget) },
    storage: { units: 'J/m² for surfaceEnergy (slab and ice, relative to freezing water) and airEnergy (cp T over the column), kg/m² for water, K, m for iceVolume', start: storedStart, end: storedEnd },
    daily: { from: -SPINUP_DAYS * DAY, units: 'waterResidual is the global water gained beyond evaporation less precipitation since the start, kg/m²; ice in percent; maxWind is the largest edge-normal |u|', ...daily },
  };
  await model.close();
  const json = Buffer.from(JSON.stringify(header)), pad = (4 - ((8 + json.length) % 4)) % 4;
  const blocks = [];
  for (const s of sequences) for (const f of s.frames) for (const field of f) blocks.push(Buffer.from(field.buffer));
  const head = Buffer.alloc(8); head.write('EXF1', 0); head.writeUInt32LE(json.length + pad, 4);
  mkdirSync(dirname(out), { recursive: true });
  writeFileSync(out, Buffer.concat([head, json, Buffer.alloc(pad, 32), ...blocks]));
  console.error(`wrote ${out}: ${sequences.map((s) => `${s.name} ${s.frames.length} frames`).join(', ')}, ${(json.length / 1024).toFixed(0)} KB header, ${(seconds() / 60).toFixed(1)} min`);
  process.exit(0);
}

const run = which === 'jw06' || which === 'tracks' || which === 'state' ? jw06() : which === 'heldsuarez' || which === 'radiation' ? heldSuarez() : null;
if (!run) { console.error('usage: node scripts/explainLayers.mjs jw06|heldsuarez|radiation|tracks|state [OUT], or wet [SPINUP_DAYS] [RECORD_DAYS] [OUT]'); process.exit(1); }
const { core, state, tendency } = run, K = core.K, step = createRK4Arrays([C, K * C, K * E]);
const [pi, theta, u] = state, { exnerLayer } = core.arrays, sigmaMid = core.sigmaMid;

function columnAt(i, pressure, values) {
  for (let k = K - 1; k > 0; k--) {
    const pk = sigmaMid[k] * pi[i], pu = sigmaMid[k - 1] * pi[i];
    if (pressure <= pk && pressure >= pu) { const w = Math.log(pk / pressure) / Math.log(pk / pu); return values(k) * (1 - w) + values(k - 1) * w; }
  }
  return values(pressure > sigmaMid[K - 1] * pi[i] ? K - 1 : 0);
}
const velocity = [0, 0];
function frame() {
  core.diagnose(pi, theta);
  const ps = new Int16Array(C), t850 = new Int16Array(C), wind250 = new Int16Array(C);
  for (let i = 0; i < C; i++) {
    ps[i] = Math.round((pi[i] / 100 - 1000) * 20);
    t850[i] = Math.round((columnAt(i, 85000, (k) => theta[k * C + i] * exnerLayer[k * C + i]) - 260) * 50);
    wind250[i] = Math.round(columnAt(i, 25000, (k) => { cellVelocity(mesh, u, i, velocity, 0, k * E); return Math.hypot(velocity[0], velocity[1]); }) * 100);
  }
  return [ps, t850, wind250];
}

const BINS = 36;
const zonalSums = () => ({ u: new Float64Array(BINS * K), T: new Float64Array(BINS * K), count: new Float64Array(BINS * K), samples: 0 });
const means = zonalSums();
const binOf = (i) => Math.min(BINS - 1, Math.floor((mesh.latCell[i] + Math.PI / 2) / Math.PI * BINS));
function accumulate(target = means) {
  core.diagnose(pi, theta);
  for (let i = 0; i < C; i++) {
    const bin = binOf(i);
    for (let k = 0; k < K; k++) {
      cellVelocity(mesh, u, i, velocity, 0, k * E);
      target.u[bin * K + k] += velocity[0] * mesh.areaCell[i]; target.T[bin * K + k] += theta[k * C + i] * exnerLayer[k * C + i] * mesh.areaCell[i]; target.count[bin * K + k] += mesh.areaCell[i];
    }
  }
  target.samples++;
}

if (which === 'tracks') {
  const { maxEdges, nEdgesOnCell, cellsOnCell, verticesOnCell, cellsOnVertex, cellsOnEdge, dcEdge, dvEdge, areaCell, xCell } = mesh, dSigma = core.diagnostics.dSigma, levels = core.levels;
  const surface = run.surface;
  const velocities = new Float64Array(2 * K * C), sigmaDot = new Float64Array((K + 1) * C), v2 = [0, 0];
  const scratch = [new Float64Array(C), new Float64Array(K * C), new Float64Array(K * E)];
  const parcels = [];
  for (const sigma of [0.9, 0.7]) for (let latDeg = 32; latDeg <= 62; latDeg += 6) for (let lonDeg = 80; lonDeg <= 230; lonDeg += 10) parcels.push({ lat: latDeg * Math.PI / 180, lon: lonDeg * Math.PI / 180, sigma, cell: 0 });
  const energy = [], times = [], START = 6 * DAY, END = 11 * DAY;
  const nearest = (x, y, z, start) => { let best = start, bestDot = x * xCell[3 * start] + y * xCell[3 * start + 1] + z * xCell[3 * start + 2]; for (;;) { let next = best; for (let m = 0; m < nEdgesOnCell[best]; m++) { const j = cellsOnCell[maxEdges * best + m], d = x * xCell[3 * j] + y * xCell[3 * j + 1] + z * xCell[3 * j + 2]; if (d > bestDot) { bestDot = d; next = j; } } if (next === best) return best; best = next; } };
  const triple = (x, y, z, i, j) => { const ax = xCell[3 * i], ay = xCell[3 * i + 1], az = xCell[3 * i + 2], bx = xCell[3 * j], by = xCell[3 * j + 1], bz = xCell[3 * j + 2]; return x * (ay * bz - az * by) + y * (az * bx - ax * bz) + z * (ax * by - ay * bx); };
  function weights(p) {
    const x = Math.cos(p.lat) * Math.cos(p.lon), y = Math.cos(p.lat) * Math.sin(p.lon), z = Math.sin(p.lat);
    p.cell = nearest(x, y, z, p.cell);
    let best = null, bestMin = -Infinity;
    for (let m = 0; m < nEdgesOnCell[p.cell]; m++) {
      const v = verticesOnCell[maxEdges * p.cell + m], [i, j, k] = [cellsOnVertex[3 * v], cellsOnVertex[3 * v + 1], cellsOnVertex[3 * v + 2]];
      const w = [triple(x, y, z, j, k), triple(x, y, z, k, i), triple(x, y, z, i, j)], total = w[0] + w[1] + w[2];
      const smallest = Math.min(...w) / total;
      if (smallest > bestMin) { bestMin = smallest; best = [[i, w[0] / total], [j, w[1] / total], [k, w[2] / total]]; }
    }
    return best;
  }
  function sample(p) {
    const ws = weights(p), mids = core.sigmaMid;
    let k = 0; while (k < K - 2 && p.sigma > mids[k + 1]) k++;
    const f = Math.min(1, Math.max(0, (p.sigma - mids[k]) / (mids[k + 1] - mids[k])));
    let u = 0, v = 0, th = 0, ps = 0, sd = 0;
    let kk = 0; while (kk < K - 1 && p.sigma > levels[kk + 1]) kk++;
    const g2 = Math.min(1, Math.max(0, (p.sigma - levels[kk]) / (levels[kk + 1] - levels[kk])));
    for (const [i, w] of ws) {
      u += w * ((1 - f) * velocities[2 * (k * C + i)] + f * velocities[2 * ((k + 1) * C + i)]);
      v += w * ((1 - f) * velocities[2 * (k * C + i) + 1] + f * velocities[2 * ((k + 1) * C + i) + 1]);
      th += w * ((1 - f) * theta[k * C + i] + f * theta[(k + 1) * C + i]);
      ps += w * pi[i];
      sd += w * ((1 - g2) * sigmaDot[kk * C + i] + g2 * sigmaDot[(kk + 1) * C + i]);
    }
    return { u, v, theta: th, ps, sigmaDot: sd };
  }
  function measure(time) {
    let kinetic = 0, internal = 0, mass = 0, thetaMass = 0;
    core.diagnose(pi, theta);
    for (let k = 0; k < K; k++) for (let e = 0; e < E; e++) { const i = cellsOnEdge[2 * e], j = cellsOnEdge[2 * e + 1], uu = u[k * E + e] ** 2; kinetic += 0.25 * dcEdge[e] * dvEdge[e] * uu * dSigma[k] * (pi[i] + pi[j]) / g; }
    for (let i = 0; i < C; i++) { mass += areaCell[i] * pi[i] / g; internal += areaCell[i] * surface[i] * pi[i] / g; for (let k = 0; k < K; k++) { const layer = areaCell[i] * pi[i] * dSigma[k] / g; internal += layer * cp * theta[k * C + i] * exnerLayer[k * C + i]; thetaMass += layer * theta[k * C + i]; } }
    const area = 4 * Math.PI * a * a;
    energy.push({ day: +(time / DAY).toFixed(3), kinetic: kinetic / area, internal: internal / area, mass: mass / area, thetaMass: thetaMass / area });
  }
  const steps = Math.round(12 * DAY / dt), records = [], t0 = performance.now();
  for (let n = 0; n <= steps; n++) {
    const time = n * dt;
    if (n % Math.round(3 * 3600 / dt) === 0) measure(time);
    if (time >= START - 1 && time <= END + 1) {
      core.diagnose(pi, theta);
      for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) { cellVelocity(mesh, u, i, v2, 0, k * E); velocities[2 * (k * C + i)] = v2[0]; velocities[2 * (k * C + i) + 1] = v2[1]; }
      core.tendency(state, scratch);
      for (let k = 0; k <= K; k++) for (let i = 0; i < C; i++) sigmaDot[k * C + i] = core.arrays.piSigmaDot[k * C + i] / pi[i];
      if (Math.round((time - START) / dt) % Math.round(3600 / dt) === 0) {
        times.push(time);
        records.push(parcels.map((p) => { const s = sample(p); return [p.lat * 180 / Math.PI, p.lon * 180 / Math.PI, p.sigma * s.ps / 100, s.theta]; }));
      }
      if (time < END) for (const p of parcels) {
        const s = sample(p);
        p.lat += s.v / a * dt; p.lon += s.u / (a * Math.cos(p.lat)) * dt; p.sigma = Math.min(0.995, Math.max(0.05, p.sigma + s.sigmaDot * dt));
      }
    }
    if (n % Math.round(DAY / dt) === 0) console.error(`tracks day ${(time / DAY).toFixed(0)}: ${((performance.now() - t0) / 1000).toFixed(0)} s`);
    if (n < steps) step(tendency, state, dt);
  }
  const header = { version: 1, case: 'tracks', N, K, parcels: parcels.length, times, fields: ['lat', 'lon', 'pressure', 'theta'], energy };
  const json = Buffer.from(JSON.stringify(header)), pad = (4 - ((8 + json.length) % 4)) % 4, data = new Float32Array(records.length * parcels.length * 4);
  records.forEach((r, h) => r.forEach((v, q) => data.set(v, (h * parcels.length + q) * 4)));
  const head = Buffer.alloc(8); head.write('EXT1', 0); head.writeUInt32LE(json.length + pad, 4);
  mkdirSync(dirname(out), { recursive: true });
  writeFileSync(out, Buffer.concat([head, json, Buffer.alloc(pad, 32), Buffer.from(data.buffer)]));
  console.error(`wrote ${out}: ${parcels.length} parcels, ${records.length} hours, ${energy.length} energy samples`);
  process.exit(0);
}

if (which === 'state') {
  const DAYS = 9, steps = Math.round(DAYS * DAY / dt), t0 = performance.now();
  for (let n = 1; n <= steps; n++) {
    step(tendency, state, dt);
    if (n % Math.round(DAY / dt) === 0) console.error(`state day ${(n * dt / DAY).toFixed(0)}: ${((performance.now() - t0) / 1000).toFixed(0)} s`);
  }
  const arrays = [['pi', 'Pa', ['C'], 'i', pi], ['theta', 'K', ['K', 'C'], 'k * C + i', theta], ['u', 'm/s', ['K', 'E'], 'k * E + e', u], ['surfaceGeopotential', 'm²/s²', ['C'], 'i', run.surface]];
  let offset = 0;
  const layout = arrays.map(([name, unit, shape, index, values]) => { const entry = { name, unit, shape, index, offset, length: values.length }; offset += values.length; return entry; });
  const header = {
    version: 1, case: 'jw06', N, K, C, E, day: DAYS, time: steps * dt, dt, levels: Array.from(core.levels), R, cp, p0, g, a, omega, nu4: core.nu4, nu4Theta: core.nu4Theta,
    data: 'float32, little-endian, the arrays one after another; offset and length count floats from the start of the data, after the header; layer k = 0 is the top; cells and edges in buildMesh(new Grid(N)) order; u is the wind across each edge, along its normal nEdge', arrays: layout,
  };
  const json = Buffer.from(JSON.stringify(header)), pad = (4 - ((8 + json.length) % 4)) % 4, data = new Float32Array(offset);
  for (const [n, entry] of layout.entries()) data.set(arrays[n][4], entry.offset);
  const head = Buffer.alloc(8); head.write('EXS1', 0); head.writeUInt32LE(json.length + pad, 4);
  mkdirSync(dirname(out), { recursive: true });
  writeFileSync(out, Buffer.concat([head, json, Buffer.alloc(pad, 32), Buffer.from(data.buffer)]));
  let low = 0; for (let i = 0; i < C; i++) if (pi[i] < pi[low]) low = i;
  console.error(`wrote ${out}: day ${DAYS}, ${offset} floats, deepest low ${(pi[low] / 100).toFixed(1)} hPa at ${(mesh.latCell[low] * 180 / Math.PI).toFixed(1)}° ${(mesh.lonCell[low] * 180 / Math.PI).toFixed(1)}°, ${((performance.now() - t0) / 1000).toFixed(0)} s`);
  process.exit(0);
}

if (which === 'radiation') {
  const HELD_DAYS = 60, SPINUP_DAYS = 150, RECORD_DAYS = 365, FRAME_DAYS = 5, WEATHER_DAYS = [287, 317];
  const ALBEDO = 0.07, HEAT_CAPACITY = 2.1e7, EVERY = 2;
  const surface = createSurface(mesh, core), radiation = createRadiation(mesh, core, { albedo: ALBEDO, radiationEvery: EVERY, orbit: PUBLISHED_RADIATION.orbit, shortwave: PUBLISHED_RADIATION.shortwave, aerosol: PUBLISHED_RADIATION.aerosol });
  const Ts = new Float64Array(C), q0 = new Float64Array(K * C), qc0 = new Float64Array(K * C), dry = new Float64Array(C), albedo = new Float64Array(C).fill(ALBEDO);
  const columns = [pi, theta, u, Ts, q0, qc0], heating = [null, new Float64Array(K * C), null, null, new Float64Array(K * C)], totals = {};
  const area = mesh.areaCell.reduce((sum, x) => sum + x, 0), perDay = Math.round(DAY / dt), bottom = (K - 1) * C, t0 = performance.now();
  const seconds = () => ((performance.now() - t0) / 1000).toFixed(0);
  const strongest = () => {
    let speed = 0, at = 0;
    for (let k = 0; k < K; k++) for (let i = 0; i < C; i++) {
      cellVelocity(mesh, u, i, velocity, 0, k * E);
      const s = Math.hypot(velocity[0], velocity[1]);
      if (!Number.isFinite(s)) return { speed: NaN, where: 'not finite' };
      if (s > speed) { speed = s; at = k * C + i; }
    }
    const k = Math.floor(at / C), i = at % C;
    return { speed, where: `${(sigmaMid[k] * pi[i] / 100).toFixed(1)} hPa ${(mesh.latCell[i] * 180 / Math.PI).toFixed(0)}°` };
  };

  for (let n = 1; n <= HELD_DAYS * perDay; n++) {
    step(run.tendency, state, dt);
    if (n % (10 * perDay) === 0) console.error(`held–suarez day ${n / perDay}: max wind ${strongest().speed.toFixed(0)} m/s, ${seconds()} s`);
  }
  core.diagnose(pi, theta);
  for (let i = 0; i < C; i++) Ts[i] = theta[bottom + i] * exnerLayer[bottom + i];

  const start = -SPINUP_DAYS * DAY, end = RECORD_DAYS * DAY, steps = Math.round((end - start) / dt), recordStep = Math.round(-start / dt), sixHours = Math.round(6 * 3600 / dt);
  const sequences = [{ name: 'year', every: FRAME_DAYS * DAY, from: 0, to: end }, { name: 'weather', every: 6 * 3600, from: WEATHER_DAYS[0] * DAY, to: WEATHER_DAYS[1] * DAY }].map((s) => ({ ...s, times: [], frames: [] }));
  const due = (time, s) => time >= s.from && time <= s.to + 1 && Math.abs(time / s.every - Math.round(time / s.every)) < 1e-6;
  const radiationFrame = () => {
    const [, t850, wind250] = frame(), ts = new Int16Array(C), olr = new Int16Array(C);
    for (let i = 0; i < C; i++) { ts[i] = Math.round((Ts[i] - 260) * 50); olr[i] = Math.round(radiation.outgoing[i] * 50); }
    return [ts, olr, t850, wind250];
  };
  const monthly = MONTHS.map(() => zonalSums());
  const daily = { absorbedSolar: [], outgoingLongwave: [], surfaceTemperature: [] }, today = { absorbedSolar: 0, outgoingLongwave: 0, surfaceTemperature: 0 };
  const surfaceSum = new Float64Array(C);
  let recorded = 0;
  for (let n = 0; n <= steps; n++) {
    const time = start + n * dt;
    for (const s of sequences) if (due(time, s)) { s.times.push(time); s.frames.push(radiationFrame()); }
    if (time >= 0 && time < end && (n - recordStep) % sixHours === 0) accumulate(monthly[monthOf(time / DAY)]);
    if (n === recordStep) radiation.restartSums();
    if (n === steps) break;
    step(run.frictionOnly, state, dt);
    heating[1].fill(0); heating[4].fill(0);
    core.diagnose(pi, theta, q0, qc0);
    surface.lowestWindSpeed(u);
    radiation.setTime(YEAR + time);
    radiation.apply(columns, heating, surface.windSpeed, totals, 0, C, albedo, albedo, dry, null, null, dt);
    for (let m = 0; m < K * C; m++) theta[m] += dt * heating[1][m];
    for (let i = 0; i < C; i++) Ts[i] += dt * radiation.surfaceFlux[i] / HEAT_CAPACITY;
    surface.convectiveAdjustment(pi, theta);
    let meanTs = 0;
    for (let i = 0; i < C; i++) meanTs += mesh.areaCell[i] * Ts[i];
    today.absorbedSolar += totals.absorbedSolar / area; today.outgoingLongwave += totals.outgoingLongwave / area; today.surfaceTemperature += meanTs / area;
    if (time >= 0) { for (let i = 0; i < C; i++) surfaceSum[i] += Ts[i]; recorded++; }
    if ((n + 1) % perDay === 0) {
      for (const key of Object.keys(today)) { daily[key].push(+(today[key] / perDay).toFixed(2)); today[key] = 0; }
      let low = Infinity, high = -Infinity;
      for (const t of Ts) { low = Math.min(low, t); high = Math.max(high, t); }
      const wind = strongest(), day = (n + 1) / perDay;
      console.error(`radiation day ${day} of ${steps / perDay}${time < 0 ? ' (spin-up)' : ''}: Ts ${daily.surfaceTemperature.at(-1).toFixed(2)} K (${low.toFixed(1)}–${high.toFixed(1)}), ASR ${daily.absorbedSolar.at(-1).toFixed(2)}, OLR ${daily.outgoingLongwave.at(-1).toFixed(2)} W/m², max wind ${wind.speed.toFixed(0)} m/s at ${wind.where}, ${seconds()} s`);
      if (!(wind.speed < 400 && low > 100 && high < 400)) { console.error('unstable: stopping'); process.exit(1); }
    }
  }

  const annual = zonalSums();
  for (const sums of monthly) { for (const key of ['u', 'T', 'count']) sums[key].forEach((v, m) => { annual[key][m] += v; }); annual.samples += sums.samples; }
  const finish = (sums) => ({ samples: sums.samples, u: Array.from(sums.u, (v, m) => +(v / sums.count[m]).toFixed(2)), T: Array.from(sums.T, (v, m) => +(v / sums.count[m]).toFixed(2)) });
  const zonal = (value) => {
    const sum = new Float64Array(BINS), weight = new Float64Array(BINS);
    for (let i = 0; i < C; i++) { sum[binOf(i)] += value(i) * mesh.areaCell[i]; weight[binOf(i)] += mesh.areaCell[i]; }
    return Array.from(sum, (v, b) => +(v / weight[b]).toFixed(2));
  };
  const header = {
    version: 1, case: 'radiation', N, K, dt, levels: Array.from(core.levels),
    fields: [{ name: 'surfaceTemperature', unit: 'K', scale: 1 / 50, offset: 260 }, { name: 'outgoingLongwave', unit: 'W/m²', scale: 1 / 50, offset: 0 }, { name: 'temperature850', unit: 'K', scale: 1 / 50, offset: 260 }, { name: 'wind250', unit: 'm/s', scale: 1 / 100, offset: 0 }],
    sequences: sequences.map((s) => ({ name: s.name, times: s.times })),
    time: 'seconds from the March equinox (March 20) that begins the recorded year; the sun is radiation.setTime(YEAR + time)',
    ground: { albedo: ALBEDO, heatCapacity: HEAT_CAPACITY }, radiationEvery: EVERY, heldSuarezDays: HELD_DAYS, spinupDays: SPINUP_DAYS,
    zonalMean: { bins: BINS, from: 0, to: end, sigma: Array.from(sigmaMid), ...finish(annual), months: MONTHS.map(([name], m) => ({ name, days: monthDays(m), ...finish(monthly[m]) })) },
    zonalBudget: { bins: BINS, from: 0, to: end, steps: recorded, absorbedSolar: zonal((i) => radiation.summed.absorbedSolar[i] / recorded), outgoingLongwave: zonal((i) => radiation.summed.outgoingLongwave[i] / recorded), surfaceTemperature: zonal((i) => surfaceSum[i] / recorded) },
    daily: { from: start, ...daily },
  };
  const json = Buffer.from(JSON.stringify(header)), pad = (4 - ((8 + json.length) % 4)) % 4;
  const blocks = [];
  for (const s of sequences) for (const f of s.frames) for (const field of f) blocks.push(Buffer.from(field.buffer));
  const head = Buffer.alloc(8); head.write('EXF1', 0); head.writeUInt32LE(json.length + pad, 4);
  mkdirSync(dirname(out), { recursive: true });
  writeFileSync(out, Buffer.concat([head, json, Buffer.alloc(pad, 32), ...blocks]));
  console.error(`wrote ${out}: ${sequences.map((s) => `${s.name} ${s.frames.length} frames`).join(', ')}, ${(json.length / 1024).toFixed(0)} KB header, ${((performance.now() - t0) / 60000).toFixed(1)} min`);
  process.exit(0);
}

const plan = which === 'jw06'
  ? { days: 12, sequences: [{ name: 'wave', every: 3 * 3600, from: 0, to: 12 * DAY }] }
  : { days: 200, sequences: [{ name: 'spinup', every: DAY, from: 0, to: 120 * DAY }, { name: 'settled', every: 6 * 3600, from: 170 * DAY, to: 200 * DAY }], meanFrom: 100 * DAY };
const sequences = plan.sequences.map((s) => ({ ...s, times: [], frames: [] }));
const t0 = performance.now(), steps = Math.round(plan.days * DAY / dt), mass0 = core.mass(pi);
for (let n = 0; n <= steps; n++) {
  const time = n * dt;
  for (const s of sequences) if (time >= s.from && time <= s.to + 1 && Math.abs(time / s.every - Math.round(time / s.every)) < 1e-6) { s.times.push(time); s.frames.push(frame()); }
  if (plan.meanFrom !== undefined && time >= plan.meanFrom && Math.round(time / dt) % Math.round(6 * 3600 / dt) === 0) accumulate();
  if (n % Math.round(DAY / dt) === 0) {
    let min = Infinity, max = -Infinity; for (const p of pi) { min = Math.min(min, p); max = Math.max(max, p); }
    console.error(`${which} day ${(time / DAY).toFixed(0)}: surface pressure ${(min / 100).toFixed(1)}–${(max / 100).toFixed(1)} hPa, mass drift ${((core.mass(pi) - mass0) / mass0).toExponential(1)}, ${((performance.now() - t0) / 1000).toFixed(0)} s`);
  }
  if (n < steps) step(tendency, state, dt);
}

const header = {
  version: 1, case: which, N, K, dt, levels: Array.from(core.levels), fields: [{ name: 'surfacePressure', unit: 'hPa', scale: 1 / 20, offset: 1000 }, { name: 'temperature850', unit: 'K', scale: 1 / 50, offset: 260 }, { name: 'wind250', unit: 'm/s', scale: 1 / 100, offset: 0 }],
  sequences: sequences.map((s) => ({ name: s.name, times: s.times })),
};
if (plan.meanFrom !== undefined) {
  header.zonalMean = { bins: BINS, from: plan.meanFrom, to: plan.days * DAY, samples: means.samples, sigma: Array.from(sigmaMid), u: Array.from(means.u, (v, n) => +(v / means.count[n]).toFixed(3)), T: Array.from(means.T, (v, n) => +(v / means.count[n]).toFixed(3)) };
}
const json = Buffer.from(JSON.stringify(header)), pad = (4 - ((8 + json.length) % 4)) % 4;
const blocks = [];
for (const s of sequences) for (const f of s.frames) for (const field of f) blocks.push(Buffer.from(field.buffer));
const head = Buffer.alloc(8); head.write('EXF1', 0); head.writeUInt32LE(json.length + pad, 4);
mkdirSync(new URL('../explain/data/', import.meta.url), { recursive: true });
writeFileSync(out, Buffer.concat([head, json, Buffer.alloc(pad, 32), ...blocks]));
console.error(`wrote ${out}: ${sequences.map((s) => `${s.name} ${s.frames.length} frames`).join(', ')}, ${((performance.now() - t0) / 60000).toFixed(1)} min`);
