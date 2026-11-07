/*
 * The surface exchange of a run on a model whose saved state is loaded and
 * whose columns are diagnosed (core.phaseColumn): exchange.cell for every
 * cell as the physics phase calls it, so that exchange.drag, heat and wind,
 * and with them surface.stress, are the run's rather than the constant
 * SEA_DRAG and LAND_DRAG a new model starts with. A save carries no
 * boundary-layer depth, which sets the convective gust: the first pass
 * stands in for the step before, the boundary layer diagnoses its depth
 * from it, and the second pass uses that depth. A save carries no sea skin
 * either, so the skin offset is left out.
 */
import { FREEZING_POINT } from '../js/physics/ice.module.js';

export function diagnoseExchange(model, dt) {
  const { core, state, surface, exchange, seaIce, land, boundaryLayer, geography } = model;
  if (!exchange || exchange.fixed) return;
  const { K, C, g, geopotential } = core.diagnostics, bottom = (K - 1) * C;
  const [pi, theta, u, skin, q, qc, ice] = state, moist = !!model.moist, landMask = geography ? geography.land : null;
  surface.lowestWindSpeed(u);
  const pass = () => {
    for (let i = 0; i < C; i++) {
      const onLand = land && landMask[i], area = onLand ? 0 : seaIce.cover(i, ice[i]);
      const fluxT = !onLand && ice[i] > 0 && area < 1 ? area * skin[i] + (1 - area) * FREEZING_POINT : skin[i];
      const before = onLand ? land.wetness(i, exchange.heat[i] * exchange.wind[i], skin[i]) : 0;
      const depthAbove = boundaryLayer ? boundaryLayer.depth[i] - geopotential[bottom + i] / g : 0;
      exchange.cell(i, pi, theta, moist ? q : null, moist ? qc : null, fluxT, surface.windSpeed[i], area, land ? land.snow[i] : 0, land ? land.vegetation[i] : 0, land ? land.canopy[i] : 0, before, depthAbove);
    }
  };
  pass();
  if (boundaryLayer) boundaryLayer.diagnose(state, 0, C, dt);
  pass();
}
