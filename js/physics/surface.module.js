import { cellVector } from '../dynamics/operators.module.js';

/*
 * Surface drag and the dry convective adjustment. Bulk aerodynamic drag
 * acts on the lowest layer with a gustiness floor on the wind speed,
 * applied to edge velocities with the rate averaged from the two adjacent
 * cells, so it only ever removes kinetic energy. An optional Rayleigh drag
 * above topSigma, ramping to 1/topDragDays at the model top, absorbs what
 * reaches the lid. heatLayers returns the kinetic energy both drags
 * remove as heat in the cells whose edges lost it.
 */
export function createSurface(mesh, core, { dragCoefficient = 1.5e-3, dragCoefficients = null, gustiness = 3, topSigma = 0.05, topDragDays = 0, buffers = null } = {}) {
  const { K, C, E, dSigma, sigmaMid, R, g, cp, exnerLayer } = core.diagnostics;
  const { cellsOnEdge, nEdgesOnCell, edgesOnCell, maxEdges, dcEdge, dvEdge, areaCell } = mesh;
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

  function surfaceDragRates(pi, theta, rates) {
    for (let i = 0; i < C; i++) {
      const idx = bottom * C + i;
      const airTemperature = theta[idx] * exnerLayer[idx];
      const airDensity = pi[i] * sigmaMid[bottom] / (R * airTemperature);
      const massPerArea = pi[i] * dSigma[bottom] / g;
      rates[i] = (dragCoefficients ? dragCoefficients[i] : dragCoefficient) * airDensity * Math.max(windSpeed[i], gustiness) / massPerArea;
    }
  }
  const topDragRate = (k) => (topDragDays > 0 && sigmaMid[k] < topSigma ? (topSigma - sigmaMid[k]) / topSigma / (topDragDays * 86400) : 0);

  function applyLayers(state, out, kFrom = 0, kTo = K) {
    const [pi, theta, u] = state;
    const [, , dU] = out;
    if (kFrom <= bottom && bottom < kTo) {
      surfaceDragRates(pi, theta, dragRate);
      for (let e = 0; e < E; e++) {
        const rate = 0.5 * (dragRate[cellsOnEdge[2 * e]] + dragRate[cellsOnEdge[2 * e + 1]]);
        dU[bottom * E + e] -= rate * u[bottom * E + e];
      }
    }
    for (let k = kFrom; k < kTo; k++) {
      const rate = topDragRate(k);
      if (rate > 0) for (let e = 0; e < E; e++) dU[k * E + e] -= rate * u[k * E + e];
    }
  }

  const heatRate = new Float64Array(C);
  function heatLayers(state, out, kFrom = 0, kTo = K) {
    const [pi, theta, u] = state;
    const dTheta = out[1];
    for (let k = kFrom; k < kTo; k++) {
      const top = topDragRate(k);
      if (k !== bottom && top === 0) continue;
      if (k === bottom) surfaceDragRates(pi, theta, heatRate);
      for (let i = 0; i < C; i++) {
        let power = 0;
        for (let m = 0; m < nEdgesOnCell[i]; m++) {
          const e = edgesOnCell[maxEdges * i + m], ue = u[k * E + e];
          const rate = top + (k === bottom ? 0.5 * (heatRate[cellsOnEdge[2 * e]] + heatRate[cellsOnEdge[2 * e + 1]]) : 0);
          power += 0.5 * dcEdge[e] * dvEdge[e] * rate * ue * ue;
        }
        dTheta[k * C + i] += power / areaCell[i] / (cp * exnerLayer[k * C + i]);
      }
    }
  }

  const aeroFactor = new Float64Array(C);

  function stress(state, out = new Float64Array(E)) {
    const [pi, theta, u] = state;
    for (let i = 0; i < C; i++) {
      const idx = bottom * C + i;
      const airDensity = pi[i] * sigmaMid[bottom] / (R * theta[idx] * exnerLayer[idx]);
      aeroFactor[i] = (dragCoefficients ? dragCoefficients[i] : dragCoefficient) * airDensity * Math.max(windSpeed[i], gustiness);
    }
    for (let e = 0; e < E; e++) {
      out[e] = 0.5 * (aeroFactor[cellsOnEdge[2 * e]] + aeroFactor[cellsOnEdge[2 * e + 1]]) * u[bottom * E + e];
    }
    return out;
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
  /*
   * Dry convective adjustment in one pass up the column: each layer
   * starts a block, and a block cooler in θ than the block beneath it
   * (by more than one part in 10⁶) merges into it, θ averaged with the
   * Exner-and-mass weights that conserve enthalpy and q, qc with the
   * mass weights, until the stack is stable; merged blocks are written
   * back. This is the state repeated pairwise mixing converges to.
   */
  const blockTop = new Int32Array(K), blockHeat = new Float64Array(K), blockWeight = new Float64Array(K), blockMass = new Float64Array(K), blockQ = new Float64Array(K), blockQc = new Float64Array(K);
  function convectiveAdjustColumn(i, pi, theta, q = null, qc = null) {
    core.diagnoseColumn(i, pi, theta, q, qc);
    let n = 0, mixes = 0;
    for (let k = K - 1; k >= 0; k--) {
      const idx = k * C + i, w = exnerLayer[idx] * dSigma[k];
      blockTop[n] = k; blockHeat[n] = theta[idx] * w; blockWeight[n] = w; blockMass[n] = dSigma[k];
      blockQ[n] = q ? q[idx] * dSigma[k] : 0; blockQc[n] = qc ? qc[idx] * dSigma[k] : 0;
      n++;
      while (n > 1 && blockHeat[n - 2] / blockWeight[n - 2] > (blockHeat[n - 1] / blockWeight[n - 1]) * (1 + 1e-6)) {
        blockHeat[n - 2] += blockHeat[n - 1]; blockWeight[n - 2] += blockWeight[n - 1]; blockMass[n - 2] += blockMass[n - 1];
        blockQ[n - 2] += blockQ[n - 1]; blockQc[n - 2] += blockQc[n - 1]; blockTop[n - 2] = blockTop[n - 1];
        n--; mixes++;
      }
    }
    if (!mixes) return 0;
    let lowest = K - 1;
    for (let b = 0; b < n; b++) {
      const top = blockTop[b];
      if (top < lowest) {
        const mixed = blockHeat[b] / blockWeight[b], mixedQ = blockQ[b] / blockMass[b], mixedQc = blockQc[b] / blockMass[b];
        for (let k = top; k <= lowest; k++) {
          const idx = k * C + i;
          theta[idx] = mixed;
          if (q) q[idx] = mixedQ;
          if (qc) qc[idx] = mixedQc;
        }
      }
      lowest = top - 1;
    }
    return mixes;
  }

  function convectiveAdjustment(pi, theta, iFrom = 0, iTo = C, q = null, qc = null) {
    let mixes = 0;
    for (let i = iFrom; i < iTo; i++) mixes += convectiveAdjustColumn(i, pi, theta, q, qc);
    return mixes;
  }

  return { lowestWindSpeed, apply, applyLayers, heatLayers, stress, convectiveAdjustment, convectiveAdjustColumn, windSpeed, shared: { windSpeed: windBuffer } };
}
