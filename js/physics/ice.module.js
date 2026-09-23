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
 */
export function openWaterAlbedo(mu) {
  return 0.026 / (Math.pow(mu, 1.7) + 0.065) + 0.15 * (mu - 0.1) * (mu - 0.5) * (mu - 1);
}

export function createSeaIce(mesh, {
  slabHeatCapacity = 2.1e7, skinHeatCapacity = 2e5, conductivity = 2.0, minimumThickness = 0.1,
  iceDensity = 917, latentHeatFusion = 3.34e5, oceanAlbedo = null, diffuseWaterAlbedo = 0.06, iceAlbedo = 0.5, fullAlbedoThickness = 0.5,
  heatCapacity = null, buffers = null,
} = {}) {
  const C = mesh.nCells;
  const latent = iceDensity * latentHeatFusion;
  const budget = { frozen: 0, melted: 0 };
  const oceanFlux = new Float64Array(buffers && buffers.oceanFlux ? buffers.oceanFlux : new SharedArrayBuffer(8 * C));

  function albedo(thickness, mu = null) {
    const water = oceanAlbedo ?? (mu === null ? diffuseWaterAlbedo : openWaterAlbedo(mu));
    return thickness > 0 ? water + (iceAlbedo - water) * Math.min(1, thickness / fullAlbedoThickness) : water;
  }

  function energy(surfaceT, ice) {
    return ice > 0 ? skinHeatCapacity * (surfaceT - FREEZING_POINT) - latent * ice : slabHeatCapacity * (surfaceT - FREEZING_POINT);
  }

  function update(surfaceT, ice, flux, i, dt) {
    const ocean = oceanFlux[i];
    const capacity = heatCapacity ? heatCapacity[i] : slabHeatCapacity;
    if (ice[i] <= 0) {
      surfaceT[i] += dt * (flux[i] + ocean) / capacity;
      if (surfaceT[i] < FREEZING_POINT) {
        ice[i] = (FREEZING_POINT - surfaceT[i]) * capacity / latent;
        surfaceT[i] = FREEZING_POINT;
        budget.frozen += mesh.areaCell[i] * ice[i];
      }
      return;
    }
    const conduction = conductivity * (FREEZING_POINT - surfaceT[i]) / Math.max(ice[i], minimumThickness);
    surfaceT[i] += dt * (flux[i] + conduction) / skinHeatCapacity;
    let thickness = ice[i] + dt * (conduction - ocean) / latent;
    if (surfaceT[i] > MELTING_POINT) {
      const excess = (surfaceT[i] - MELTING_POINT) * skinHeatCapacity;
      surfaceT[i] = MELTING_POINT;
      thickness -= excess / latent;
    }
    if (thickness <= 0) {
      budget.melted += mesh.areaCell[i] * ice[i];
      surfaceT[i] = FREEZING_POINT + (-thickness * latent + skinHeatCapacity * (surfaceT[i] - FREEZING_POINT)) / capacity;
      ice[i] = 0;
    } else {
      budget.frozen += mesh.areaCell[i] * (thickness - ice[i]);
      ice[i] = thickness;
    }
  }

  return { albedo, energy, update, budget, slabHeatCapacity, latent, oceanFlux, shared: { oceanFlux: oceanFlux.buffer } };
}
