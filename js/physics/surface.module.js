import { cellVector } from '../dynamics/operators.module.js';

/*
 * Boundary-layer momentum sinks and the dry convective adjustment. Bulk
 * aerodynamic drag acts on the lowest layer with a gustiness floor on the
 * wind speed; Rayleigh drag ramps from zero at pblTop to pblRate at the
 * ground. Both are applied to edge velocities with the rate averaged from
 * the two adjacent cells, so they only ever remove kinetic energy.
 */
export function createSurface(mesh, core, { dragCoefficient = 1.5e-3, gustiness = 3, pblTop = 0.7, pblRate = 1 / 86400 } = {}) {
  const { K, C, E, dSigma, sigmaMid, R, g, exnerLayer } = core.diagnostics;
  const { cellsOnEdge } = mesh;
  const bottom = K - 1;
  const vector = new Float64Array(3 * C);
  const windSpeed = new Float64Array(C);
  const dragRate = new Float64Array(C);

  function lowestWindSpeed(u) {
    cellVector(mesh, u.subarray(bottom * E, K * E), vector);
    for (let i = 0; i < C; i++) windSpeed[i] = Math.hypot(vector[3 * i], vector[3 * i + 1], vector[3 * i + 2]);
    return windSpeed;
  }

  function apply(state, out) {
    const [pi, theta, u] = state;
    const [, , dU] = out;
    for (let i = 0; i < C; i++) {
      const idx = bottom * C + i;
      const airTemperature = theta[idx] * exnerLayer[idx];
      const airDensity = pi[i] * sigmaMid[bottom] / (R * airTemperature);
      const massPerArea = pi[i] * dSigma[bottom] / g;
      dragRate[i] = dragCoefficient * airDensity * Math.max(windSpeed[i], gustiness) / massPerArea;
    }
    for (let e = 0; e < E; e++) {
      const rate = 0.5 * (dragRate[cellsOnEdge[2 * e]] + dragRate[cellsOnEdge[2 * e + 1]]);
      dU[bottom * E + e] -= rate * u[bottom * E + e];
    }
    for (let k = 0; k < K; k++) {
      if (sigmaMid[k] <= pblTop) continue;
      const rate = pblRate * (sigmaMid[k] - pblTop) / (1 - pblTop);
      for (let e = 0; e < E; e++) dU[k * E + e] -= rate * u[k * E + e];
    }
  }

  /*
   * Mixes every statically unstable layer pair of column i to a common θ,
   * weighted by Exner·Δσ so column enthalpy Σ cp θ Π Δσ is unchanged,
   * sweeping until the column is stable. Refreshes the column's Exner
   * ratios for the current π first.
   */
  function convectiveAdjustColumn(i, pi, theta) {
    core.diagnoseColumn(i, pi, theta);
    let mixes = 0, dirty = true, guard = 0;
    while (dirty && guard < K * K) {
      dirty = false;
      guard++;
      for (let k = K - 2; k >= 0; k--) {
        const above = k * C + i, below = above + C;
        if (theta[below] > theta[above] * (1 + 1e-9)) {
          const wAbove = exnerLayer[above] * dSigma[k];
          const wBelow = exnerLayer[below] * dSigma[k + 1];
          const mixed = (theta[above] * wAbove + theta[below] * wBelow) / (wAbove + wBelow);
          theta[above] = mixed;
          theta[below] = mixed;
          dirty = true;
          mixes++;
        }
      }
    }
    return mixes;
  }

  function convectiveAdjustment(pi, theta) {
    let mixes = 0;
    for (let i = 0; i < C; i++) mixes += convectiveAdjustColumn(i, pi, theta);
    return mixes;
  }

  return { lowestWindSpeed, apply, convectiveAdjustment, convectiveAdjustColumn, windSpeed };
}
