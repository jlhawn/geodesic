import { MELTING_POINT, SNOW_AGEING, agedSnowAlbedo, refreshedSnowAlbedo } from './ice.module.js';
import { SOLAR_CONSTANT, AXIAL_TILT } from './radiation.module.js';

/*
 * The land surface: a skin of small heat capacity, a bucket of soil
 * water in the manner of Manabe (1969) that limits evaporation as it
 * dries and spills what it cannot hold into runoff, and a snow cover in
 * water equivalent that precipitation builds when the lowest air is
 * below freezing and the surface energy melts. Evaporation draws on the
 * snow while there is any, else on the bucket; snow raises the albedo
 * toward snowAlbedo over fullSnow kg/m². update() applies the surface
 * flux and the evaporation of the physics phase; deposit() adds the
 * step's precipitation in the adjustment phase, when it is known. The
 * ocean cells of `snow` belong to the sea ice, which keeps its snow there.
 *
 * With `vegetation` each land cell also carries a vegetation cover v
 * between 0 (bare ground) and 1 (dense forest) that the climate grows:
 * snow-free, it relaxes toward a goal set by how full the bucket is (0
 * below dryWetness of the capacity, 1 above wetWetness), over
 * growthTime when rising and declineTime when falling, and under snow it
 * decays toward 0 over snowDeclineTime. The snow-free albedo runs from
 * the bare soil's to vegetatedAlbedo with v. With soilDarkening 'surface'
 * (the default) the bare soil darkens linearly with the surface layer's
 * fill, as Idso et al. (1975) found it darken with the water of the top
 * centimetres, from bareAlbedo at the first fill of darkeningWetness to
 * wetSoilAlbedo at the second (by default 0 and 1); 'rootZone' reads the
 * bucket's fill instead (by default 0.2 to 0.5), and false keeps
 * bareAlbedo. The bucket holds a fixed rootZoneCapacity whatever the
 * cover, because a soil keeps its water capacity when its plants die and
 * a browned region can therefore regreen when the rain returns. Without
 * it the bucket is bucketCapacity and the albedo `albedo` everywhere. A
 * cell of the geography's `iceSheet` grows no vegetation and keeps
 * iceSheetAlbedo whatever lies on it.
 *
 * With snowAgeing (the default) the snow's albedo is per cell, in
 * `snowAlbedo`: agedSnowAlbedo (ice.module.js) ages it toward
 * oldSnowAlbedo after each update, a snowfall refreshes it
 * (refreshedSnowAlbedo), and ground without snow holds freshSnowAlbedo,
 * so that the next snow starts fresh; without it the snow's albedo is
 * snowAlbedo. With snowMasking and vegetation, trees standing above the
 * snow darken it: the albedo of a full snow cover falls linearly from the
 * snow's own to forestSnowAlbedo as the tree cover `canopy` rises to
 * closedCanopy, and stays there above it, the shape of the MODIS
 * snow-covered albedo against tree cover (Moody et al. 2007).
 *
 * With `treeline` (the default) the canopy is a tree cover that the
 * growing season admits. Each land cell keeps two running means of the
 * lowest air's temperature over seasonMemory: `seasonLength`, the share
 * of the time it is at least seasonThreshold (°C), and `seasonWarmth`,
 * its mean excess over that threshold (K). Their ratio plus the threshold
 * is the growing season's mean temperature (Paulsen and Körner 2014: the
 * season of days at least 0.9 °C, at least 94 days long and at least
 * 6.4 °C on the mean at the treeline), taken over a season of at least
 * minimumSeason days, so a shorter one counts as cooler. The treeline
 * factor f rises from 0 to 1 as that mean rises across treelineWarmth.
 * Snow-free, the tree cover relaxes toward f v, over treeGrowthTime
 * rising and treeDeclineTime falling; under snow it holds, save that it
 * falls toward f over treeDeclineTime where it stands above it. Without
 * the treeline the canopy is the standing cover: v where that is higher,
 * otherwise relaxing toward it over canopyMemory.
 *
 * With treeMoisture (the default) water decides between forest and
 * grass where the season admits trees: each land cell keeps running
 * means over moistureMemory of its rain and snow (`rainMean`, mm/d, from
 * deposit) and of the FAO-56 reference evapotranspiration the radiation
 * gives it (`demandMean`, mm/d, from update), and f is multiplied by
 * aridityFactor of their ratio, the aridity index P/PET, ramping from 0
 * at forestAridity[0] to 1 at forestAridity[1]. With grassland (the
 * default) the cover not under trees is grass: the vegetated albedo runs
 * from grassAlbedo to forestAlbedo with the trees' share of v, and snow
 * on the grass takes its own albedo less grassSnowDarkening times the
 * grass's share of the cell before the trees mask it.
 *
 * With vegetation the soil has two stores. Rain fills a surface layer
 * of surfaceCapacity kg/m² first; what it cannot hold infiltrates the
 * root zone (the bucket), a share (soil/capacity)⁴ of it running off,
 * and the surface layer seeps into the root zone over percolationTime.
 * Bare ground evaporates from the surface layer alone, so a desert
 * dries within a day of rain; the cover transpires from the root zone
 * through stomata, at the aerodynamic rate times 1/(1 + r_s g_a) with
 * r_s = stomatalResistance divided by a warmth that rises from 0 at
 * growthColdest to 1 at growthWarmest, so cold or dry roots close
 * them. wetness(i, g_a, T) gives the radiation column that factor and
 * remembers how much of it was bare ground, which update() draws from
 * the surface layer. Growth toward the goal needs the same warmth;
 * decline does not.
 */
export const DARKENING_WETNESS = { surface: [0, 1], rootZone: [0.2, 0.5] };

/*
 * The start of the season means where a state carries none: a year of
 * the lowest air's temperature as a sine, its mean a + b Q̄ less the
 * lapse rate times the ground's height and its amplitude
 * min(k ΔQ, c₀ + c₁ ΔQ), with Q̄ and ΔQ the mean and the annual
 * harmonic of the daily insolation at the top of the atmosphere at the
 * cell's latitude (W/m²). The constants are rounded least-squares fits to the
 * model's own first year: the lowest air at the four quarter days of the
 * one-year runs nine64 and eight64, harmonically fitted per land cell
 * off the ice sheets (the amplitude's second branch over 40–85N).
 */
export const SEASON_ESTIMATE = { mean: [-31.8, 0.148], amplitude: [0.069, 19.7, -0.0205], lapseRate: 0.0065 };

export function insolationCycle(lat, samples = 73) {
  let mean = 0, cosine = 0, sine = 0;
  for (let d = 0; d < samples; d++) {
    const w = 2 * Math.PI * d / samples, declination = AXIAL_TILT * Math.sin(w);
    const c = -Math.tan(lat) * Math.tan(declination), hour = c >= 1 ? 0 : c <= -1 ? Math.PI : Math.acos(c);
    const q = SOLAR_CONSTANT / Math.PI * (hour * Math.sin(lat) * Math.sin(declination) + Math.cos(lat) * Math.cos(declination) * Math.sin(hour));
    mean += q / samples; cosine += 2 * q * Math.cos(w) / samples; sine += 2 * q * Math.sin(w) / samples;
  }
  return { mean, amplitude: Math.hypot(cosine, sine) };
}

/*
 * Of a year whose temperature (°C) is mean + amplitude cos(2πt): the
 * share at or above `threshold` and the mean excess over it (K).
 */
export function sineSeason(mean, amplitude, threshold) {
  const x = amplitude > 0 ? (threshold - mean) / amplitude : (mean >= threshold ? -1 : 1);
  if (x >= 1) return { length: 0, warmth: 0 };
  if (x <= -1) return { length: 1, warmth: mean - threshold };
  const half = Math.acos(x);
  return { length: half / Math.PI, warmth: ((mean - threshold) * half + amplitude * Math.sin(half)) / Math.PI };
}

export function seasonEstimate(lat, elevation = 0, threshold = 0.9) {
  const { mean, amplitude } = insolationCycle(lat);
  const [a, b] = SEASON_ESTIMATE.mean, [k, c0, c1] = SEASON_ESTIMATE.amplitude;
  return sineSeason(a + b * mean - SEASON_ESTIMATE.lapseRate * Math.max(0, elevation), Math.max(0, Math.min(k * amplitude, c0 + c1 * amplitude)), threshold);
}

export function treelineFactor(length, warmth, { seasonThreshold = 0.9, minimumSeason = 94, treelineWarmth = [6.4, 8.0] } = {}) {
  const [treeless, treed] = treelineWarmth;
  return Math.min(1, Math.max(0, (seasonThreshold + warmth / Math.max(length, minimumSeason / 365) - treeless) / (treed - treeless)));
}

export const FOREST_ARIDITY = [0.2, 1.0];

export function aridityFactor(rain, demand, [dry, wet] = FOREST_ARIDITY) {
  return Math.min(1, Math.max(0, (rain / Math.max(demand, 1e-3) - dry) / (wet - dry)));
}

/*
 * The start of the moisture means where a state carries none: the
 * demand a + b Q̄ (mm/d) of the annual mean insolation at the top of the
 * atmosphere (W/m²), a rounded least-squares fit per land cell off the
 * ice sheets to four ten-day means of the model's first year (nine64
 * and eight64 at days 91, 183, 274 and 365), and the rain that demand
 * times an aridity index c + d × the bucket's fill + e × the cover,
 * fitted to the last year's rain of the 21 regions in five64's log over
 * the estimated demand against the year-six state's fill and cover.
 */
export const MOISTURE_ESTIMATE = { demand: [-2.14, 0.0134], aridity: [0.01, 0.79, 0.63] };

export function moistureEstimate(lat, fill, cover = 0.5) {
  const [d0, d1] = MOISTURE_ESTIMATE.demand, [a0, a1, a2] = MOISTURE_ESTIMATE.aridity;
  const demand = Math.max(0, d0 + d1 * insolationCycle(lat).mean);
  return { rain: demand * Math.max(0, a0 + a1 * Math.min(1, Math.max(0, fill)) + a2 * Math.min(1, Math.max(0, cover))), demand };
}

export function createLandSurface(mesh, geography, {
  heatCapacity = 1e6, bucketCapacity = 150, wetnessThreshold = 0.75, albedo = 0.2, snowAlbedo = 0.55, fullSnow = 20,
  latentHeatFusion = 3.34e5, vegetation: vegetated = true, bareAlbedo = 0.30, vegetatedAlbedo = 0.13, rootZoneCapacity = 300, dryWetness = 0.1, wetWetness = 0.6, growthTime = 180 * 86400, declineTime = 365 * 86400,
  snowDeclineTime = 720 * 86400, iceSheetAlbedo = 0.8, surfaceCapacity = 15, percolationTime = 86400, stomatalResistance = 70,
  growthColdest = 278.15, growthWarmest = 288.15, soilDarkening = 'surface', wetSoilAlbedo = 0.15, darkeningWetness = null,
  snowAgeing = true, oldSnowAlbedo = 0.50, freshSnowAlbedo = SNOW_AGEING.freshSnowAlbedo, coldSnowAgeing = SNOW_AGEING.coldSnowAgeing, meltingSnowAgeing = SNOW_AGEING.meltingSnowAgeing,
  refreshSnowfall = SNOW_AGEING.refreshSnowfall, wetSnowRange = SNOW_AGEING.wetSnowRange, ageingActivation = SNOW_AGEING.ageingActivation,
  snowMasking = true, forestSnowAlbedo = 0.27, closedCanopy = 0.7, canopyMemory = 365 * 86400,
  treeline = true, seasonThreshold = 0.9, minimumSeason = 94, treelineWarmth = [6.4, 8.0], seasonMemory = 3 * 365 * 86400, treeGrowthTime = 10 * 365 * 86400, treeDeclineTime = 3 * 365 * 86400,
  treeMoisture = true, moistureMemory = 3 * 365 * 86400, forestAridity = FOREST_ARIDITY, grassland = true, forestAlbedo = 0.13, grassAlbedo = 0.20, grassSnowDarkening = 0.06, buffers = null,
} = {}) {
  const C = mesh.nCells;
  const shared = (name) => new Float64Array(buffers && buffers[name] ? buffers[name] : new SharedArrayBuffer(8 * C));
  const soil = shared('soil'), snow = shared('snow'), runoff = shared('runoff'), vegetation = shared('vegetation'), surface = shared('surface');
  const snowAlbedoField = shared('snowAlbedo'), canopy = shared('canopy'), seasonLength = shared('seasonLength'), seasonWarmth = shared('seasonWarmth');
  const rainMean = shared('rainMean'), demandMean = shared('demandMean');
  if (!(buffers && buffers.snowAlbedo)) snowAlbedoField.fill(freshSnowAlbedo);
  const ageing = { freshSnowAlbedo, coldSnowAgeing, meltingSnowAgeing, refreshSnowfall, wetSnowRange, ageingActivation };
  const masked = snowMasking && vegetated, treed = treeline && vegetated, gated = treed && treeMoisture, grassy = grassland && vegetated;
  const bareShare = new Float64Array(C);
  const warmth = (t) => Math.min(1, Math.max(0, (t - growthColdest) / (growthWarmest - growthColdest)));
  const { land, iceSheet = null } = geography;
  const onIceSheet = (i) => iceSheet !== null && iceSheet[i] > 0;
  const budget = { runoff: 0, melt: 0 };

  const capacity = () => (vegetated ? rootZoneCapacity : bucketCapacity);
  const darkening = soilDarkening === true ? 'surface' : soilDarkening;
  if (darkening !== false && !DARKENING_WETNESS[darkening]) throw new Error(`soilDarkening is 'surface', 'rootZone' or false, not ${soilDarkening}`);
  const [darkeningFrom, darkeningTo] = darkeningWetness ?? DARKENING_WETNESS[darkening || 'surface'];
  if (!(darkeningTo > darkeningFrom && darkeningFrom >= 0)) throw new Error(`darkeningWetness must rise from its first to its second fill, not ${darkeningWetness}`);
  const wetFill = (i) => (darkening === 'surface' ? surface[i] / surfaceCapacity : soil[i] / capacity(i));
  const soilAlbedo = (i) => (darkening ? bareAlbedo - (bareAlbedo - wetSoilAlbedo) * Math.min(1, Math.max(0, (wetFill(i) - darkeningFrom) / (darkeningTo - darkeningFrom))) : bareAlbedo);
  const [treelessWarmth, treedWarmth] = treelineWarmth;
  if (!(treedWarmth > treelessWarmth)) throw new Error(`treelineWarmth must rise from its first to its second temperature, not ${treelineWarmth}`);
  const seasonKelvin = MELTING_POINT + seasonThreshold;
  if (!(forestAridity[1] > forestAridity[0])) throw new Error(`forestAridity must rise from its first to its second index, not ${forestAridity}`);
  const moistureFactor = (i) => (gated ? aridityFactor(rainMean[i], demandMean[i], forestAridity) : 1);
  const treeFactor = (i) => treelineFactor(seasonLength[i], seasonWarmth[i], { seasonThreshold, minimumSeason, treelineWarmth }) * moistureFactor(i);
  const grass = (i) => Math.max(0, vegetation[i] - canopy[i]);
  const coverAlbedo = (i) => (grassy ? (vegetation[i] > 0 ? grassAlbedo + (forestAlbedo - grassAlbedo) * Math.min(1, canopy[i] / vegetation[i]) : grassAlbedo) : vegetatedAlbedo);
  const bareGround = (i) => { if (!vegetated) return albedo; const s = soilAlbedo(i); return s + (coverAlbedo(i) - s) * vegetation[i]; };

  function overflow(i) {
    const cap = capacity(i);
    if (soil[i] > cap) {
      const excess = soil[i] - cap;
      runoff[i] += excess;
      budget.runoff += mesh.areaCell[i] * excess;
      soil[i] = cap;
    }
  }

  function wetness(i, aero = 0, temperature = growthWarmest) {
    if (snow[i] > 0) { bareShare[i] = 0; return 1; }
    const roots = Math.min(1, soil[i] / (wetnessThreshold * capacity(i)));
    if (!vegetated) { bareShare[i] = 0; return roots; }
    const bare = (1 - vegetation[i]) * Math.min(1, surface[i] / surfaceCapacity);
    const canopy = vegetation[i] * roots / (1 + stomatalResistance * aero / Math.max(0.05, warmth(temperature)));
    const total = bare + canopy;
    bareShare[i] = total > 0 ? bare / total : 0;
    return total;
  }

  function snowCovered(i) {
    const own = (snowAgeing ? snowAlbedoField[i] : snowAlbedo) - (grassy ? grassSnowDarkening * grass(i) : 0);
    return masked ? own + (forestSnowAlbedo - own) * Math.min(1, canopy[i] / closedCanopy) : own;
  }

  function surfaceAlbedo(i) {
    if (onIceSheet(i)) return iceSheetAlbedo;
    const bare = bareGround(i);
    return bare + Math.min(1, snow[i] / fullSnow) * (snowCovered(i) - bare);
  }

  function moisture(i, potential, dt) {
    demandMean[i] += (86400 * potential - demandMean[i]) * (1 - Math.exp(-dt / moistureMemory));
  }

  function season(i, airTemperature, dt) {
    const keep = 1 - Math.exp(-dt / seasonMemory);
    seasonLength[i] += ((airTemperature >= seasonKelvin ? 1 : 0) - seasonLength[i]) * keep;
    seasonWarmth[i] += (Math.max(0, airTemperature - seasonKelvin) - seasonWarmth[i]) * keep;
  }

  function grow(i, dt, temperature) {
    if (onIceSheet(i)) { vegetation[i] = 0; canopy[i] = 0; return; }
    if (snow[i] > 0) vegetation[i] *= Math.exp(-dt / snowDeclineTime);
    else {
      const cap = capacity(i);
      const goal = Math.min(1, Math.max(0, (Math.min(soil[i], cap) / cap - dryWetness) / (wetWetness - dryWetness)));
      if (goal > vegetation[i]) vegetation[i] += (goal - vegetation[i]) * (1 - Math.exp(-dt * warmth(temperature) / growthTime));
      else vegetation[i] += (goal - vegetation[i]) * (1 - Math.exp(-dt / declineTime));
    }
    if (!treed) { canopy[i] = Math.max(vegetation[i], canopy[i] + (vegetation[i] - canopy[i]) * (1 - Math.exp(-dt / canopyMemory))); return; }
    const f = treeFactor(i), goal = snow[i] > 0 ? Math.min(canopy[i], f) : f * vegetation[i];
    canopy[i] += (goal - canopy[i]) * (1 - Math.exp(-dt / (goal > canopy[i] ? treeGrowthTime : treeDeclineTime)));
  }

  function update(i, surfaceT, flux, evaporation, dt, airTemperature = null, potential = null) {
    const temperature = surfaceT[i];
    surfaceT[i] += dt * flux[i] / heatCapacity;
    let left = evaporation * dt;
    const fromSnow = Math.min(snow[i], left);
    snow[i] -= fromSnow; left -= fromSnow;
    surfaceT[i] -= latentHeatFusion * fromSnow / heatCapacity;
    if (vegetated) {
      const fromSurface = Math.min(surface[i], left * bareShare[i]);
      surface[i] -= fromSurface; left -= fromSurface;
      const seep = surface[i] * (1 - Math.exp(-dt / percolationTime));
      surface[i] -= seep; soil[i] += seep;
    }
    soil[i] = Math.max(0, soil[i] - left);
    if (snow[i] > 0 && surfaceT[i] > MELTING_POINT) {
      const energy = (surfaceT[i] - MELTING_POINT) * heatCapacity;
      const melt = Math.min(snow[i], energy / latentHeatFusion);
      snow[i] -= melt;
      soil[i] += melt;
      budget.melt += mesh.areaCell[i] * melt;
      surfaceT[i] = MELTING_POINT + (energy - melt * latentHeatFusion) / heatCapacity;
    }
    if (vegetated && airTemperature !== null) season(i, airTemperature, dt);
    if (vegetated && potential !== null) moisture(i, potential, dt);
    if (vegetated) grow(i, dt, temperature);
    snowAlbedoField[i] = snow[i] > 0 ? agedSnowAlbedo(snowAlbedoField[i], surfaceT[i], dt, oldSnowAlbedo, ageing) : freshSnowAlbedo;
    overflow(i);
  }

  function deposit(i, rain, airTemperature, dt = null) {
    if (vegetated && dt) rainMean[i] += (86400 * rain / dt - rainMean[i]) * (1 - Math.exp(-dt / moistureMemory));
    if (airTemperature < MELTING_POINT) { snow[i] += rain; snowAlbedoField[i] = refreshedSnowAlbedo(snowAlbedoField[i], rain, ageing); return; }
    if (!vegetated) { soil[i] += rain; overflow(i); return; }
    surface[i] += rain;
    if (surface[i] > surfaceCapacity) {
      const infiltration = surface[i] - surfaceCapacity;
      surface[i] = surfaceCapacity;
      const shed = infiltration * Math.pow(Math.min(1, soil[i] / capacity(i)), 4);
      runoff[i] += shed;
      budget.runoff += mesh.areaCell[i] * shed;
      soil[i] += infiltration - shed;
      overflow(i);
    }
  }

  let estimated = null;
  function estimateSeason(i) {
    estimated ??= Array.from({ length: C }, (_, n) => (land[n] ? seasonEstimate(mesh.latCell[n], geography.elevation ? geography.elevation[n] : 0, seasonThreshold) : { length: 0, warmth: 0 }));
    seasonLength[i] = estimated[i].length; seasonWarmth[i] = estimated[i].warmth;
  }
  function estimateMoisture(i) {
    if (!land[i]) { rainMean[i] = 0; demandMean[i] = 0; return; }
    const e = moistureEstimate(mesh.latCell[i], soil[i] / capacity(i), vegetation[i]);
    rainMean[i] = e.rain; demandMean[i] = e.demand;
  }
  const startingTrees = (i) => (treed ? treeFactor(i) * vegetation[i] : vegetation[i]);

  /*
   * A fresh land surface is half vegetated with half-full buckets, so
   * that forests and deserts both have to emerge; its season means are
   * SEASON_ESTIMATE's and its tree cover the cover times their treeline
   * factor.
   */
  function initialize() {
    for (let i = 0; i < C; i++) {
      vegetation[i] = land[i] && vegetated && !onIceSheet(i) ? 0.5 : 0;
      soil[i] = land[i] ? 0.5 * capacity(i) : 0;
      estimateSeason(i);
      estimateMoisture(i);
      canopy[i] = startingTrees(i);
      snow[i] = 0; runoff[i] = 0; surface[i] = 0; snowAlbedoField[i] = freshSnowAlbedo;
    }
  }

  function water() {
    let total = 0;
    for (let i = 0; i < C; i++) if (land[i]) total += mesh.areaCell[i] * (soil[i] + surface[i] + snow[i] + runoff[i]);
    return total;
  }

  /*
   * A saved land state without vegetation starts as initialize() would
   * where it is free of snow (vegetated, bucket full) and bare under snow.
   * A sea cell keeps the saved snow only where `ice` has ice to hold it.
   * A state saved without a snow albedo starts its snow fresh
   * (freshSnowAlbedo). One saved without season means starts them at
   * SEASON_ESTIMATE's and its tree cover, as a fresh start does, at the
   * cover times their treeline factor; without the treeline a state
   * without a standing cover stands it at the cover. One saved without
   * the moisture means starts them at moistureEstimate's from its own
   * bucket's fill and, under the moisture gate, keeps its trees where
   * they stand below the cover times their potential and lowers them to
   * it elsewhere.
   */
  function load(saved, ice = null) {
    const seasoned = !!(saved.seasonLength && saved.seasonWarmth), moistened = !!(saved.rainMean && saved.demandMean);
    for (let i = 0; i < C; i++) {
      soil[i] = land[i] ? saved.soil[i] : 0;
      surface[i] = land[i] && saved.surface ? Math.min(surfaceCapacity, saved.surface[i]) : 0;
      snow[i] = land[i] || (ice && ice[i] > 0) ? saved.snow[i] : 0;
      snowAlbedoField[i] = snow[i] > 0 && saved.snowAlbedo ? saved.snowAlbedo[i] : freshSnowAlbedo;
      if (!land[i] || !vegetated || onIceSheet(i)) vegetation[i] = 0;
      else if (saved.vegetation) vegetation[i] = Math.min(1, Math.max(0, saved.vegetation[i]));
      else { vegetation[i] = snow[i] > 0 ? 0 : 1; if (snow[i] <= 0) soil[i] = capacity(i); }
      if (!land[i]) { seasonLength[i] = 0; seasonWarmth[i] = 0; }
      else if (seasoned) { seasonLength[i] = Math.min(1, Math.max(0, saved.seasonLength[i])); seasonWarmth[i] = Math.max(0, saved.seasonWarmth[i]); }
      else estimateSeason(i);
      if (land[i] && moistened) { rainMean[i] = Math.max(0, saved.rainMean[i]); demandMean[i] = Math.max(0, saved.demandMean[i]); }
      else estimateMoisture(i);
      if (!treed) canopy[i] = vegetation[i] > 0 && saved.canopy ? Math.min(1, Math.max(vegetation[i], saved.canopy[i])) : vegetation[i];
      else if (!land[i] || !vegetated || onIceSheet(i)) canopy[i] = 0;
      else if (!seasoned || !saved.canopy) canopy[i] = startingTrees(i);
      else canopy[i] = moistened || !gated ? Math.min(1, Math.max(0, saved.canopy[i])) : Math.min(1, Math.max(0, saved.canopy[i]), startingTrees(i));
    }
    runoff.fill(0);
  }

  return {
    soil, surface, snow, runoff, vegetation, snowAlbedo: snowAlbedoField, canopy, seasonLength, seasonWarmth, rainMean, demandMean, treeFactor, moistureFactor, land, budget, heatCapacity, latentHeatFusion, bucketCapacity, capacity, wetness, albedo: surfaceAlbedo, update, deposit, initialize, water, load,
    serialize() { return { soil: Float64Array.from(soil), snow: Float64Array.from(snow), snowAlbedo: Float64Array.from(snowAlbedoField), ...(vegetated ? { vegetation: Float64Array.from(vegetation), surface: Float64Array.from(surface), canopy: Float64Array.from(canopy), seasonLength: Float64Array.from(seasonLength), seasonWarmth: Float64Array.from(seasonWarmth), rainMean: Float64Array.from(rainMean), demandMean: Float64Array.from(demandMean) } : {}) }; },
    shared: { soil: soil.buffer, surface: surface.buffer, snow: snow.buffer, runoff: runoff.buffer, vegetation: vegetation.buffer, snowAlbedo: snowAlbedoField.buffer, canopy: canopy.buffer, seasonLength: seasonLength.buffer, seasonWarmth: seasonWarmth.buffer, rainMean: rainMean.buffer, demandMean: demandMean.buffer },
  };
}
