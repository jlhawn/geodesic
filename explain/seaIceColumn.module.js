/*
 * One cell of the model's own sea ice, createSeaIce on a cell of unit area
 * over a mixed layer `depth` meters deep, stepped hourly from open water a
 * little above freezing (or from `startIce` meters of ice) on model day
 * `startDay`. Each hour the cell takes the day's mean sunlight from the
 * model's sun at `latitude` (the samples of sunTrack that dailyMean
 * averages), of which `transmission` reaches it, less what the model's
 * albedos reflect: the ice's, the ageing snow's and the open water's
 * openWaterAlbedo averaged over the day's sun weighted by its light, which
 * is exact because every albedo of the cell is linear in the water's. A toy
 * atmosphere, not the model's, takes A + B(T − 273.15) from the cell, T its
 * area-weighted surface temperature as the model's fluxes see it, and the
 * winds bring `convergence`; the currents bring `oceanHeat` to the water,
 * under ice to the ice's base as the ocean's oceanFlux does. Snow falls at
 * `snowfall` kg/m² of water a year through the hours the surface is below
 * the melting point. The mixed layer's salinity changes only by brine, the
 * change in the ice's volume per unit area (cover times thickness) at the
 * ocean's ICE_SALINITY and ICE_DENSITY, as the ocean books it: no
 * evaporation, rain or melt water freshens it, and the layer never mixes
 * with the water below. The record is daily, states at the end of each day
 * and fluxes as day means; temperatures are kelvin.
 */
import { createSeaIce, openWaterAlbedo, FREEZING_POINT, MELTING_POINT } from '../js/physics/ice.module.js';
import { SOLAR_CONSTANT } from '../js/physics/radiation.module.js';
import { seawaterDensity } from '../js/ocean/seawater.module.js';
import { LAYER_DENSITIES } from '../js/ocean/layered.module.js';
import { sunTrack, DAYS } from './insolation.module.js';

export const DAY = 86400, STEP = 3600, RHO_CP = 1025 * 3985, ICE_SALINITY = 5, ICE_DENSITY = 917;
export const DEEP_CLASS = LAYER_DENSITIES[LAYER_DENSITIES.indexOf(1026.95)];
export const DARK_ICE = { iceAlbedo: 0.07, meltingIceAlbedo: 0.07, iceSnowAlbedo: 0.07 };
export const TOY = { A: 200, B: 2.1, transmission: 0.7 };
export const DEPTHS = { slab: 4.2e7 / RHO_CP, floor: 50, deep: 200 };
const unit = () => ({ nCells: 1, areaCell: [1] });

let darkest = 0;
function darkestSun() {
  if (darkest) return darkest;
  let best = 1;
  for (let m = 1; m >= 0.01; m -= 1e-4) if (openWaterAlbedo(m) < openWaterAlbedo(best)) best = m;
  return (darkest = best);
}

export function equivalentSun(albedo) {
  let low = 1e-3, high = darkestSun();
  if (albedo <= openWaterAlbedo(high)) return high;
  for (let n = 0; n < 60; n++) {
    const mid = 0.5 * (low + high);
    if (openWaterAlbedo(mid) > albedo) low = mid; else high = mid;
  }
  return 0.5 * (low + high);
}

const tables = new Map();
export function dailySunlight(latitude, perDay = 288) {
  const key = `${latitude}:${perDay}`;
  if (tables.has(key)) return tables.get(key);
  const c = Math.cos(latitude * Math.PI / 180), s = Math.sin(latitude * Math.PI / 180);
  const mean = new Float64Array(DAYS), water = new Float64Array(DAYS), sun = new Float64Array(DAYS);
  for (let d = 0; d < DAYS; d++) {
    const { x, z } = sunTrack(d, perDay);
    let lit = 0, reflected = 0;
    for (let k = 0; k < perDay; k++) {
      const mu = c * x[k] + s * z[k];
      if (mu > 0) { lit += mu; reflected += mu * openWaterAlbedo(mu); }
    }
    mean[d] = SOLAR_CONSTANT * lit / perDay;
    water[d] = lit > 0 ? reflected / lit : NaN;
    sun[d] = lit > 0 ? equivalentSun(water[d]) : NaN;
  }
  const table = { mean, water, sun };
  tables.set(key, table);
  return table;
}

export function brine(salinity, grown, depth) {
  return (salinity - ICE_SALINITY) * grown * ICE_DENSITY / 1000 / depth;
}

export function iceForDensity({ salinity = 34, depth = 50, density = DEEP_CLASS } = {}) {
  let low = 0, high = 100;
  for (let n = 0; n < 60; n++) {
    const mid = 0.5 * (low + high);
    if (seawaterDensity(FREEZING_POINT, salinity + brine(salinity, mid, depth)) < density) low = mid; else high = mid;
  }
  return 0.5 * (low + high);
}

export function iceGrowth({ thickness = 1, below = 20, days = 1, depth = 50, dt = STEP } = {}) {
  const sea = createSeaIce(unit(), { heatCapacity: Float64Array.of(RHO_CP * depth) });
  const T = new Float64Array(1), ice = Float64Array.of(thickness), still = new Float64Array(1);
  sea.load(ice);
  let conduction = null;
  for (let n = 0; n < Math.round(days * DAY / dt); n++) {
    const before = ice[0];
    T[0] = FREEZING_POINT - below;
    sea.update(T, ice, still, 0, dt);
    conduction ??= (ice[0] - before) * sea.latent / dt;
  }
  return { conduction, grown: ice[0] - thickness, perDay: (ice[0] - thickness) / days, thickness: ice[0] };
}

export function seaIceColumn({
  latitude = 75, depth = 50, ice: iceAlbedos = 'bright', convergence = 85, oceanHeat = 0, snowfall = 0, salinity = 34, years = 24,
  A = TOY.A, B = TOY.B, transmission = TOY.transmission, dt = STEP, start = FREEZING_POINT + 2, startDay = 0, startIce = 0, iceOptions = {},
} = {}) {
  const sea = createSeaIce(unit(), { heatCapacity: Float64Array.of(RHO_CP * depth), ...(iceAlbedos === 'dark' ? DARK_ICE : {}), ...iceOptions });
  const light = dailySunlight(latitude), perDay = Math.round(DAY / dt), days = years * DAYS;
  const T = Float64Array.of(startIce > 0 ? FREEZING_POINT : start), ice = Float64Array.of(startIce), flux = new Float64Array(1);
  sea.load(ice);
  const fields = ['thickness', 'cover', 'coverMax', 'volume', 'snow', 'surface', 'water', 'salinity', 'density', 'sunlight', 'absorbed', 'outgoing', 'net', 'albedo'];
  const record = { day: new Int16Array(days), ...Object.fromEntries(fields.map((name) => [name, new Float32Array(days)])) };
  const snowRate = snowfall / (DAYS * DAY);
  let S = salinity, volume = startIce;
  const waterT = () => (ice[0] > 0 ? FREEZING_POINT : T[0]);

  for (let d = 0; d < days; d++) {
    const day = (startDay + d) % DAYS, sun = transmission * light.mean[day], mu = Number.isNaN(light.sun[day]) ? null : light.sun[day];
    let absorbed = 0, outgoing = 0, net = 0, coverMax = 0;
    for (let n = 0; n < perDay; n++) {
      const h = ice[0], area = sea.cover(0, h), skin = T[0], snow = sea.snow[0], snowy = sea.snowAlbedo[0];
      const taken = sun * (1 - sea.albedo(h, mu, snow, area, skin, snowy));
      const surface = h > 0 && area < 1 ? area * skin + (1 - area) * FREEZING_POINT : skin;
      const lost = A + B * (surface - MELTING_POINT);
      flux[0] = taken - lost + convergence;
      sea.oceanFlux[0] = oceanHeat;
      sea.update(T, ice, flux, 0, dt, h > 0 && area < 1 ? sun * sea.albedoContrast(h, mu, snow, skin, snowy) : 0);
      if (snowRate > 0) sea.deposit(0, snowRate * dt, T[0], ice, T);
      const now = sea.cover(0, ice[0]) * ice[0];
      S += brine(S, now - volume, depth);
      volume = now;
      absorbed += taken / perDay; outgoing += lost / perDay; net += (flux[0] + oceanHeat) / perDay;
      coverMax = Math.max(coverMax, sea.cover(0, ice[0]));
    }
    const area = sea.cover(0, ice[0]);
    record.day[d] = day;
    record.thickness[d] = ice[0]; record.cover[d] = area; record.coverMax[d] = coverMax; record.volume[d] = volume; record.snow[d] = sea.snow[0];
    record.surface[d] = ice[0] > 0 && area < 1 ? area * T[0] + (1 - area) * FREEZING_POINT : T[0];
    record.water[d] = waterT(); record.salinity[d] = S; record.density[d] = seawaterDensity(waterT(), S);
    record.sunlight[d] = sun; record.absorbed[d] = absorbed; record.outgoing[d] = outgoing; record.net[d] = net;
    record.albedo[d] = sun > 0 ? 1 - absorbed / sun : NaN;
  }

  function year(y = years - 1) {
    const from = y * DAYS, to = from + DAYS;
    const span = (name) => {
      let low = Infinity, high = -Infinity;
      for (let d = from; d < to; d++) { low = Math.min(low, record[name][d]); high = Math.max(high, record[name][d]); }
      return [low, high];
    };
    let iceFree = 0, open = 0;
    for (let d = from; d < to; d++) { if (record.coverMax[d] === 0) iceFree++; if (record.cover[d] < 0.15) open++; }
    const [thinnest, thickest] = span('volume'), [freshest, saltiest] = span('salinity'), [waterLow, waterHigh] = span('water'), [surfaceLow, surfaceHigh] = span('surface'), [lightest, densest] = span('density');
    return {
      year: y, thinnest, thickest, thickestIce: span('thickness')[1], deepestSnow: span('snow')[1], iceFreeDays: iceFree, openDays: open,
      waterLow, waterHigh, waterRange: waterHigh - waterLow, surfaceLow, surfaceHigh,
      salinityLow: freshest, salinityHigh: saltiest, salinityGained: saltiest - freshest, salinityDrift: record.salinity[to - 1] - (from > 0 ? record.salinity[from - 1] : salinity),
      densityLow: lightest, densityHigh: densest, againstClass: densest - DEEP_CLASS,
    };
  }

  return { latitude, depth, years, days, record, sea, budget: sea.budget, year, sunlight: light };
}
