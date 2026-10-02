export const FREEZING_POINT = 271.35;
export const MELTING_POINT = 273.15;
export const MINIMUM_CONCENTRATION = 0.01;
export const MINIMUM_VOLUME = 1e-4;

/*
 * The slab ocean and its sea ice, a zero-layer thermodynamic model in
 * the manner of Semtner (1976). Open water is a mixed layer of heat
 * capacity slabHeatCapacity; when it cools to the freezing point the
 * deficit freezes into ice of that latent heat. Ice has a skin of small
 * heat capacity whose temperature answers the surface flux and the
 * conduction through the ice from its base at the freezing point; the
 * conducted heat freezes water onto the base, and a skin that would
 * pass the melting point melts the ice from the top instead. Ice that
 * melts away returns its leftover energy to the mixed layer. The surface
 * energy — mixed-layer heat over the freezing point, skin heat, minus
 * the ice's latent heat — changes by exactly the surface flux plus
 * `oceanFlux`, the heat above freezing that the ocean's mixed layer
 * hands to the base of the ice. surfaceT is the skin temperature of
 * the ice, or of the water where there is none. Open water reflects
 * the direct beam with the zenith-angle albedo of Briegleb et al.
 * (1986), 0.02 under a high sun and 0.3 near the horizon, and diffuse
 * light (`albedo` without a zenith cosine) with diffuseWaterAlbedo,
 * unless oceanAlbedo fixes both. With `heatCapacity` given, that per-cell
 * array (a dynamic ocean's upper layer) replaces slabHeatCapacity in the
 * cell update.
 *
 * Ice covers the fraction `concentration` of its cell (A), with the
 * thickness `ice` over that part, so the volume per cell area is A·h;
 * the open part, the leads, is water held at the freezing point. The
 * atmosphere sees the area-weighted albedo and, for its fluxes, the
 * area-weighted surface temperature. The cell's net flux is split
 * between the parts by the sunlight the leads absorb beyond the ice
 * (`contrast` in update) and by the heat the leads at the freezing
 * point lose beyond the ice, leadExchange per kelvin of the ice skin's
 * difference from freezing, so that the ice's flux per unit ice area
 * and the water's per unit water area average back to it. The ice part
 * evolves as above; the heat the leads gain or lose, with the ocean's
 * flux under them, melts ice or freezes new ice. Melting takes area as
 * Hibler (1979) does, half the relative loss of volume from the area;
 * ice frozen in the leads closes them as new ice leadClosing thick. Open
 * water that cools below freezing forms ice leadClosing thick over the
 * area its volume covers. Ice below MINIMUM_CONCENTRATION of its cell
 * or MINIMUM_VOLUME of volume melts away, so water cools a little below
 * freezing before its first ice forms. When the area shrinks the lost
 * area's snow melts, its latent heat drawn from the water at freezing,
 * which freezes the same mass onto the ice, and its skin heat goes to
 * the ice volume; new area joins at the freezing point without snow.
 * A cell whose ice carries no concentration is fully covered.
 *
 * Snow lies on the ice: deposit() adds precipitation that falls on an
 * iced cell from air below the melting point, in water equivalent per
 * unit ice area, to `snow` — the ocean cells of the array the land
 * surface keeps its own snow in, when its buffer is shared. Snow
 * brightens the surface toward iceSnowAlbedo over iceFullSnow kg/m²,
 * conducts in series with the ice (snowConductivity over its depth at
 * snowDensity), melts before the ice does, and goes into the water when
 * the ice is gone, its latent heat drawn from the mixed layer. Snow
 * that falls on open water melts at once, cooling the water by that
 * latent heat; what falls on the leads of an iced cell freezes the same
 * mass of water onto the ice. The ocean's freshwater counts
 * precipitation when it falls, snow or not. Snow heavier than the ice's
 * freeboard (waterDensity − iceDensity per metre of ice) floods and
 * freezes into snow-ice: the surplus mass leaves `snow` and joins the
 * ice at iceDensity, which conserves both mass and energy.
 *
 * The bare ice's albedo depends on its skin temperature as CCSM3's sea
 * ice does (Briegleb et al. 2004): iceAlbedo while cold, falling linearly
 * over the last iceMeltingRange kelvin below the melting point to
 * meltingIceAlbedo, which stands for melting ice with its ponds; ice
 * thinner than fullAlbedoThickness blends toward the water's albedo. With
 * snowAgeing the snow on the ice has its own albedo in `snowAlbedo` (the
 * ocean cells of the land's array when its buffer is shared), aged by
 * agedSnowAlbedo toward iceSnowFloor and refreshed by snowfall
 * (refreshedSnowAlbedo); without it the snow's albedo is iceSnowAlbedo.
 */
export function openWaterAlbedo(mu) {
  return 0.026 / (Math.pow(mu, 1.7) + 0.065) + 0.15 * (mu - 0.1) * (mu - 0.5) * (mu - 1);
}

/*
 * The snow-albedo ageing of Douville et al. (1995) as the ECMWF land
 * scheme carries it (Dutra et al. 2010, appendix eq. A7 and eq. 9): snow
 * whose surface lies within wetSnowRange kelvin of the melting point
 * relaxes toward `floor` at meltingSnowAgeing per day; colder snow loses
 * coldSnowAgeing per day, down to `floor`, scaled with ageingActivation
 * by exp(ageingActivation (1/T_melt − 1/T)), the temperature dependence of
 * grain growth in BATS (Dickinson et al. 1993), so that cold dry snow
 * keeps its brightness. A snowfall of `fall` kg/m² moves the albedo
 * min(1, fall / refreshSnowfall) of the way to freshSnowAlbedo.
 */
export const SNOW_AGEING = { freshSnowAlbedo: 0.85, coldSnowAgeing: 0.008, meltingSnowAgeing: 0.24, refreshSnowfall: 10, wetSnowRange: 2, ageingActivation: 5000 };
export function agedSnowAlbedo(albedo, temperature, dt, floor, { coldSnowAgeing, meltingSnowAgeing, wetSnowRange, ageingActivation } = SNOW_AGEING) {
  const days = dt / 86400;
  if (temperature >= MELTING_POINT - wetSnowRange) return floor + (albedo - floor) * Math.exp(-meltingSnowAgeing * days);
  const pace = ageingActivation > 0 ? Math.min(1, Math.exp(ageingActivation * (temperature - MELTING_POINT) / (MELTING_POINT * temperature))) : 1;
  return Math.max(floor, albedo - coldSnowAgeing * pace * days);
}
export function refreshedSnowAlbedo(albedo, fall, { freshSnowAlbedo, refreshSnowfall } = SNOW_AGEING) {
  return albedo + Math.min(1, fall / refreshSnowfall) * (freshSnowAlbedo - albedo);
}

export function createSeaIce(mesh, {
  slabHeatCapacity = 2.1e7, skinHeatCapacity = 2e5, conductivity = 2.0, minimumThickness = 0.1,
  iceDensity = 917, latentHeatFusion = 3.34e5, oceanAlbedo = null, diffuseWaterAlbedo = 0.06, iceAlbedo = 0.62, meltingIceAlbedo = 0.48, iceMeltingRange = 1, fullAlbedoThickness = 0.5,
  iceSnowAlbedo = 0.75, iceFullSnow = 20, snowConductivity = 0.31, snowDensity = 300, waterDensity = 1026, leadClosing = 0.3, leadExchange = 10,
  snowAgeing = true, iceSnowFloor = 0.70, freshSnowAlbedo = SNOW_AGEING.freshSnowAlbedo, coldSnowAgeing = SNOW_AGEING.coldSnowAgeing, meltingSnowAgeing = SNOW_AGEING.meltingSnowAgeing,
  refreshSnowfall = SNOW_AGEING.refreshSnowfall, wetSnowRange = SNOW_AGEING.wetSnowRange, ageingActivation = SNOW_AGEING.ageingActivation,
  heatCapacity = null, buffers = null,
} = {}) {
  const C = mesh.nCells;
  const latent = iceDensity * latentHeatFusion;
  const budget = { frozen: 0, melted: 0, snowfall: 0, snowMelted: 0, snowIce: 0, leadFrozen: 0, lateralMelted: 0 };
  const oceanFlux = new Float64Array(buffers && buffers.oceanFlux ? buffers.oceanFlux : new SharedArrayBuffer(8 * C));
  const snow = new Float64Array(buffers && buffers.snow ? buffers.snow : new SharedArrayBuffer(8 * C));
  const concentration = new Float64Array(buffers && buffers.concentration ? buffers.concentration : new SharedArrayBuffer(8 * C));
  const snowAlbedo = new Float64Array(buffers && buffers.snowAlbedo ? buffers.snowAlbedo : new SharedArrayBuffer(8 * C));
  if (!(buffers && buffers.snowAlbedo)) snowAlbedo.fill(freshSnowAlbedo);
  const ageing = { freshSnowAlbedo, coldSnowAgeing, meltingSnowAgeing, refreshSnowfall, wetSnowRange, ageingActivation };

  const waterAlbedo = (mu) => oceanAlbedo ?? (mu === null ? diffuseWaterAlbedo : openWaterAlbedo(mu));
  const bareIceAlbedo = (temperature) => iceAlbedo + (meltingIceAlbedo - iceAlbedo) * Math.min(1, Math.max(0, (temperature - MELTING_POINT + iceMeltingRange) / iceMeltingRange));
  function coverAlbedo(thickness, water, snowCover, temperature, snowy) {
    if (thickness <= 0) return water;
    const bare = water + (bareIceAlbedo(temperature) - water) * Math.min(1, thickness / fullAlbedoThickness);
    return bare + ((snowAgeing ? snowy : iceSnowAlbedo) - bare) * Math.min(1, snowCover / iceFullSnow);
  }

  /*
   * temperature is the ice's skin temperature (cold when not given) and
   * snowy the albedo of the snow on it (fresh when not given).
   */
  function albedo(thickness, mu = null, snowCover = 0, fraction = thickness > 0 ? 1 : 0, temperature = -Infinity, snowy = freshSnowAlbedo) {
    const water = waterAlbedo(mu);
    return fraction * coverAlbedo(thickness, water, snowCover, temperature, snowy) + (1 - fraction) * water;
  }

  function albedoContrast(thickness, mu = null, snowCover = 0, temperature = -Infinity, snowy = freshSnowAlbedo) {
    const water = waterAlbedo(mu);
    return coverAlbedo(thickness, water, snowCover, temperature, snowy) - water;
  }

  function cover(i, thickness) {
    return thickness > 0 ? (concentration[i] > 0 ? concentration[i] : 1) : 0;
  }

  function load(ice, saved = null) {
    for (let i = 0; i < C; i++) concentration[i] = ice[i] > 0 ? (saved && saved[i] > 0 ? Math.min(1, saved[i]) : 1) : 0;
  }

  function energy(surfaceT, ice, snowCover = 0, area = ice > 0 ? 1 : 0) {
    return ice > 0 ? area * (skinHeatCapacity * (surfaceT - FREEZING_POINT) - latent * ice - latentHeatFusion * snowCover) : slabHeatCapacity * (surfaceT - FREEZING_POINT);
  }

  function deposit(i, amount, airTemperature, ice, surfaceT) {
    if (airTemperature >= MELTING_POINT || amount <= 0) return false;
    budget.snowfall += mesh.areaCell[i] * amount;
    if (ice[i] > 0) {
      const area = cover(i, ice[i]), leads = (1 - area) * amount;
      snow[i] += amount;
      ice[i] += leads / (iceDensity * area);
      snowAlbedo[i] = refreshedSnowAlbedo(snowAlbedo[i], amount, ageing);
      budget.snowMelted += mesh.areaCell[i] * leads;
      budget.frozen += mesh.areaCell[i] * leads / iceDensity;
    } else { surfaceT[i] -= latentHeatFusion * amount / (heatCapacity ? heatCapacity[i] : slabHeatCapacity); budget.snowMelted += mesh.areaCell[i] * amount; }
    return true;
  }

  function update(surfaceT, ice, flux, i, dt, contrast = 0) {
    const ocean = oceanFlux[i];
    const capacity = heatCapacity ? heatCapacity[i] : slabHeatCapacity;
    const cellArea = mesh.areaCell[i];
    if (ice[i] <= 0) {
      concentration[i] = 0;
      if (snow[i] > 0) { surfaceT[i] -= latentHeatFusion * snow[i] / capacity; budget.snowMelted += cellArea * snow[i]; snow[i] = 0; }
      snowAlbedo[i] = freshSnowAlbedo;
      surfaceT[i] += dt * (flux[i] + ocean) / capacity;
      if (surfaceT[i] < FREEZING_POINT) {
        const volume = (FREEZING_POINT - surfaceT[i]) * capacity / latent, area = Math.min(1, volume / leadClosing);
        if (area >= MINIMUM_CONCENTRATION && volume >= MINIMUM_VOLUME) {
          ice[i] = volume / area;
          concentration[i] = area;
          surfaceT[i] = FREEZING_POINT;
          budget.frozen += cellArea * volume;
        }
      }
      return;
    }
    const h = ice[i], A = cover(i, h);
    let T = surfaceT[i], s = snow[i];
    const split = contrast - leadExchange * (FREEZING_POINT - T);
    const iceFlux = flux[i] - (1 - A) * split, waterFlux = flux[i] + A * split;
    const conduction = (FREEZING_POINT - T) / (Math.max(h, minimumThickness) / conductivity + s / (snowDensity * snowConductivity));
    T += dt * (iceFlux + conduction) / skinHeatCapacity;
    let thickness = h + dt * (conduction - ocean) / latent;
    if (T > MELTING_POINT) {
      let excess = (T - MELTING_POINT) * skinHeatCapacity;
      T = MELTING_POINT;
      const fromSnow = Math.min(s, excess / latentHeatFusion);
      s -= fromSnow;
      budget.snowMelted += cellArea * A * fromSnow;
      excess -= fromSnow * latentHeatFusion;
      thickness -= excess / latent;
    }
    const leadHeat = (1 - A) * (waterFlux + ocean) * dt;
    const grown = A * (thickness - h), leadMelt = Math.max(0, leadHeat / latent), leadIce = Math.max(0, -leadHeat / latent);
    const melted = Math.max(0, -grown) + leadMelt;
    budget.leadFrozen += cellArea * leadIce;
    budget.lateralMelted += cellArea * leadMelt;
    const energyLeft = A * skinHeatCapacity * (T - FREEZING_POINT) - latent * (A * thickness - leadHeat / latent) - latentHeatFusion * A * s;
    let volume = A * thickness - leadHeat / latent;
    const area = Math.min(1, A - melted / (2 * h) + (1 - A) * leadIce / leadClosing);
    let skin = T, snowOnIce = s;
    if (area < A) volume += (A - area) * (s / iceDensity - skinHeatCapacity * (T - FREEZING_POINT) / latent);
    else if (area > A) { skin = FREEZING_POINT + A * (T - FREEZING_POINT) / area; snowOnIce = A * s / area; }
    if (area < MINIMUM_CONCENTRATION || volume < MINIMUM_VOLUME) {
      budget.melted += cellArea * A * h;
      budget.snowMelted += cellArea * A * s;
      surfaceT[i] = FREEZING_POINT + energyLeft / capacity;
      ice[i] = 0;
      snow[i] = 0;
      snowAlbedo[i] = freshSnowAlbedo;
      concentration[i] = 0;
      return;
    }
    if (area < A) budget.snowMelted += cellArea * (A - area) * s;
    budget.frozen += cellArea * (volume - A * h);
    const newThickness = volume / area;
    const flooded = Math.max(0, snowOnIce - (waterDensity - iceDensity) * newThickness) * iceDensity / waterDensity;
    surfaceT[i] = skin;
    snow[i] = snowOnIce - flooded;
    snowAlbedo[i] = snow[i] > 0 ? agedSnowAlbedo(snowAlbedo[i], skin, dt, iceSnowFloor, ageing) : freshSnowAlbedo;
    budget.snowIce += cellArea * area * flooded;
    ice[i] = newThickness + flooded / iceDensity;
    concentration[i] = area;
  }

  return {
    albedo, albedoContrast, cover, load, energy, update, deposit, budget, slabHeatCapacity, latent, latentHeatFusion, leadClosing, leadExchange, oceanFlux, snow, snowAlbedo, concentration,
    shared: { oceanFlux: oceanFlux.buffer, snow: snow.buffer, snowAlbedo: snowAlbedo.buffer, concentration: concentration.buffer },
  };
}
