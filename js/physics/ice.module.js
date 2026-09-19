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
 * the ice's latent heat — changes by exactly the surface flux.
 * surfaceT is the skin temperature the atmosphere sees in both states.
 */
export function createSeaIce(mesh, {
  slabHeatCapacity = 2.1e7, skinHeatCapacity = 2e5, conductivity = 2.0, minimumThickness = 0.1,
  iceDensity = 917, latentHeatFusion = 3.34e5, oceanAlbedo = 0.07, iceAlbedo = 0.6, fullAlbedoThickness = 0.5,
} = {}) {
  const latent = iceDensity * latentHeatFusion;
  const budget = { frozen: 0, melted: 0 };

  function albedo(thickness) {
    return thickness > 0 ? oceanAlbedo + (iceAlbedo - oceanAlbedo) * Math.min(1, thickness / fullAlbedoThickness) : oceanAlbedo;
  }

  function energy(surfaceT, ice) {
    return ice > 0 ? skinHeatCapacity * (surfaceT - FREEZING_POINT) - latent * ice : slabHeatCapacity * (surfaceT - FREEZING_POINT);
  }

  function update(surfaceT, ice, flux, i, dt) {
    if (ice[i] <= 0) {
      surfaceT[i] += dt * flux[i] / slabHeatCapacity;
      if (surfaceT[i] < FREEZING_POINT) {
        ice[i] = (FREEZING_POINT - surfaceT[i]) * slabHeatCapacity / latent;
        surfaceT[i] = FREEZING_POINT;
        budget.frozen += mesh.areaCell[i] * ice[i];
      }
      return;
    }
    const conduction = conductivity * (FREEZING_POINT - surfaceT[i]) / Math.max(ice[i], minimumThickness);
    surfaceT[i] += dt * (flux[i] + conduction) / skinHeatCapacity;
    let thickness = ice[i] + dt * conduction / latent;
    if (surfaceT[i] > MELTING_POINT) {
      const excess = (surfaceT[i] - MELTING_POINT) * skinHeatCapacity;
      surfaceT[i] = MELTING_POINT;
      thickness -= excess / latent;
    }
    if (thickness <= 0) {
      budget.melted += mesh.areaCell[i] * ice[i];
      surfaceT[i] = FREEZING_POINT + (-thickness * latent + skinHeatCapacity * (surfaceT[i] - FREEZING_POINT)) / slabHeatCapacity;
      ice[i] = 0;
    } else {
      budget.frozen += mesh.areaCell[i] * (thickness - ice[i]);
      ice[i] = thickness;
    }
  }

  return { albedo, energy, update, budget, slabHeatCapacity, latent };
}
