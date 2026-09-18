import { cellVector } from '../dynamics/operators.module.js';

/*
 * Boundary-layer momentum sinks and the dry convective adjustment. Bulk
 * aerodynamic drag acts on the lowest layer with a gustiness floor on the
 * wind speed; Rayleigh drag ramps from zero at pblTop to pblRate at the
 * ground. Both are applied to edge velocities with the rate averaged from
 * the two adjacent cells, so they only ever remove kinetic energy. An
 * optional Rayleigh drag above topSigma, ramping to 1/topDragDays at the
 * model top, absorbs what reaches the lid.
 */
export function createSurface(mesh, core, { dragCoefficient = 1.5e-3, gustiness = 3, pblTop = 0.7, pblRate = 1 / 86400, topSigma = 0.05, topDragDays = 0, buffers = null } = {}) {
  const { K, C, E, dSigma, sigmaMid, R, g, exnerLayer } = core.diagnostics;
  const { cellsOnEdge } = mesh;
  const bottom = K - 1;
  const vector = new Float64Array(3 * C);
  const windBuffer = buffers && buffers.windSpeed ? buffers.windSpeed : new SharedArrayBuffer(8 * C);
  const windSpeed = new Float64Array(windBuffer);
  const dragRate = new Float64Array(C);

  function lowestWindSpeed(u, iFrom = 0, iTo = C) {
    cellVector(mesh, u.subarray(bottom * E, K * E), vector, iFrom, iTo);
    for (let i = iFrom; i < iTo; i++) windSpeed[i] = Math.hypot(vector[3 * i], vector[3 * i + 1], vector[3 * i + 2]);
    return windSpeed;
  }

  function applyLayers(state, out, kFrom = 0, kTo = K) {
    const [pi, theta, u] = state;
    const [, , dU] = out;
    if (kFrom <= bottom && bottom < kTo) {
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
    }
    for (let k = kFrom; k < kTo; k++) {
      if (sigmaMid[k] <= pblTop) continue;
      const rate = pblRate * (sigmaMid[k] - pblTop) / (1 - pblTop);
      for (let e = 0; e < E; e++) dU[k * E + e] -= rate * u[k * E + e];
    }
    if (topDragDays > 0) {
      for (let k = kFrom; k < kTo; k++) {
        if (sigmaMid[k] >= topSigma) continue;
        const rate = (topSigma - sigmaMid[k]) / topSigma / (topDragDays * 86400);
        for (let e = 0; e < E; e++) dU[k * E + e] -= rate * u[k * E + e];
      }
    }
  }

  function apply(state, out) {
    applyLayers(state, out, 0, K);
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

  function convectiveAdjustment(pi, theta, iFrom = 0, iTo = C) {
    let mixes = 0;
    for (let i = iFrom; i < iTo; i++) mixes += convectiveAdjustColumn(i, pi, theta);
    return mixes;
  }

  return { lowestWindSpeed, apply, applyLayers, convectiveAdjustment, convectiveAdjustColumn, windSpeed, shared: { windSpeed: windBuffer } };
}
