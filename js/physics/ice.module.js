export const FREEZING_POINT = 271.35;
export const MELTING_POINT = 273.15;

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
 * hands to the base of the ice. surfaceT is the skin
 * temperature the atmosphere sees in both states. Open water reflects
 * the direct beam with the zenith-angle albedo of Briegleb et al.
 * (1986), 0.02 under a high sun and 0.3 near the horizon, and diffuse
 * light (`albedo` without a zenith cosine) with diffuseWaterAlbedo,
 * unless oceanAlbedo fixes both. With `heatCapacity` given, that per-cell
 * array (a dynamic ocean's upper layer) replaces slabHeatCapacity in the
 * cell update.
 *
 * Snow lies on the ice: deposit() adds precipitation that falls on an
 * iced cell from air below the melting point, in water equivalent, to
 * `snow` — the ocean cells of the array the land surface keeps its own
 * snow in, when its buffer is shared. Snow brightens the surface toward
 * iceSnowAlbedo over iceFullSnow kg/m², conducts in series with the ice
 * (snowConductivity over its depth at snowDensity), melts before the ice
 * does, and goes into the water when the ice is gone, its latent heat
 * drawn from the mixed layer. Snow that falls on open water melts at
 * once, cooling the water by that latent heat. The ocean's freshwater
 * counts precipitation when it falls, snow or not. Snow heavier than
 * the ice's freeboard (waterDensity − iceDensity per metre of ice)
 * floods and freezes into snow-ice: the surplus mass leaves `snow` and
 * joins the ice at iceDensity, which conserves both mass and energy.
 */
export function openWaterAlbedo(mu) {
  return 0.026 / (Math.pow(mu, 1.7) + 0.065) + 0.15 * (mu - 0.1) * (mu - 0.5) * (mu - 1);
}

export function createSeaIce(mesh, {
  slabHeatCapacity = 2.1e7, skinHeatCapacity = 2e5, conductivity = 2.0, minimumThickness = 0.1,
  iceDensity = 917, latentHeatFusion = 3.34e5, oceanAlbedo = null, diffuseWaterAlbedo = 0.06, iceAlbedo = 0.5, fullAlbedoThickness = 0.5,
  iceSnowAlbedo = 0.75, iceFullSnow = 20, snowConductivity = 0.31, snowDensity = 300, waterDensity = 1026,
  heatCapacity = null, buffers = null,
} = {}) {
  const C = mesh.nCells;
  const latent = iceDensity * latentHeatFusion;
  const budget = { frozen: 0, melted: 0, snowfall: 0, snowMelted: 0, snowIce: 0 };
  const oceanFlux = new Float64Array(buffers && buffers.oceanFlux ? buffers.oceanFlux : new SharedArrayBuffer(8 * C));
  const snow = new Float64Array(buffers && buffers.snow ? buffers.snow : new SharedArrayBuffer(8 * C));

  function albedo(thickness, mu = null, snowCover = 0) {
    const water = oceanAlbedo ?? (mu === null ? diffuseWaterAlbedo : openWaterAlbedo(mu));
    if (thickness <= 0) return water;
    const bare = water + (iceAlbedo - water) * Math.min(1, thickness / fullAlbedoThickness);
    return bare + (iceSnowAlbedo - bare) * Math.min(1, snowCover / iceFullSnow);
  }

  function energy(surfaceT, ice, snowCover = 0) {
    return ice > 0 ? skinHeatCapacity * (surfaceT - FREEZING_POINT) - latent * ice - latentHeatFusion * snowCover : slabHeatCapacity * (surfaceT - FREEZING_POINT);
  }

  function deposit(i, amount, airTemperature, ice, surfaceT) {
    if (airTemperature >= MELTING_POINT || amount <= 0) return false;
    budget.snowfall += mesh.areaCell[i] * amount;
    if (ice[i] > 0) snow[i] += amount;
    else { surfaceT[i] -= latentHeatFusion * amount / (heatCapacity ? heatCapacity[i] : slabHeatCapacity); budget.snowMelted += mesh.areaCell[i] * amount; }
    return true;
  }

  function update(surfaceT, ice, flux, i, dt) {
    const ocean = oceanFlux[i];
    const capacity = heatCapacity ? heatCapacity[i] : slabHeatCapacity;
    if (ice[i] <= 0) {
      if (snow[i] > 0) { surfaceT[i] -= latentHeatFusion * snow[i] / capacity; budget.snowMelted += mesh.areaCell[i] * snow[i]; snow[i] = 0; }
      surfaceT[i] += dt * (flux[i] + ocean) / capacity;
      if (surfaceT[i] < FREEZING_POINT) {
        ice[i] = (FREEZING_POINT - surfaceT[i]) * capacity / latent;
        surfaceT[i] = FREEZING_POINT;
        budget.frozen += mesh.areaCell[i] * ice[i];
      }
      return;
    }
    const conduction = (FREEZING_POINT - surfaceT[i]) / (Math.max(ice[i], minimumThickness) / conductivity + snow[i] / (snowDensity * snowConductivity));
    surfaceT[i] += dt * (flux[i] + conduction) / skinHeatCapacity;
    let thickness = ice[i] + dt * (conduction - ocean) / latent;
    if (surfaceT[i] > MELTING_POINT) {
      let excess = (surfaceT[i] - MELTING_POINT) * skinHeatCapacity;
      surfaceT[i] = MELTING_POINT;
      const fromSnow = Math.min(snow[i], excess / latentHeatFusion);
      snow[i] -= fromSnow;
      budget.snowMelted += mesh.areaCell[i] * fromSnow;
      excess -= fromSnow * latentHeatFusion;
      thickness -= excess / latent;
    }
    if (thickness <= 0) {
      budget.melted += mesh.areaCell[i] * ice[i];
      budget.snowMelted += mesh.areaCell[i] * snow[i];
      surfaceT[i] = FREEZING_POINT + (-thickness * latent + skinHeatCapacity * (surfaceT[i] - FREEZING_POINT) - latentHeatFusion * snow[i]) / capacity;
      ice[i] = 0;
      snow[i] = 0;
    } else {
      budget.frozen += mesh.areaCell[i] * (thickness - ice[i]);
      const flooded = Math.max(0, snow[i] - (waterDensity - iceDensity) * thickness) * iceDensity / waterDensity;
      snow[i] -= flooded;
      budget.snowIce += mesh.areaCell[i] * flooded;
      ice[i] = thickness + flooded / iceDensity;
    }
  }

  return { albedo, energy, update, deposit, budget, slabHeatCapacity, latent, latentHeatFusion, oceanFlux, snow, shared: { oceanFlux: oceanFlux.buffer, snow: snow.buffer } };
}
