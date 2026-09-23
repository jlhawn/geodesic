import { cellVector } from '../dynamics/operators.module.js';

/*
 * A diffusive planetary boundary layer in the manner of Troen and Mahrt
 * (1986), as used by the simple moist GCMs that this model follows.
 * `diagnose` finds, per cell, the boundary-layer top as the height where
 * the bulk Richardson number of the lowest layer's virtual potential
 * temperature and wind, with the convective floor 100 u*², first exceeds
 * richardsonCritical, and lays the K-profile κ u* z (1 − z/h)² over the
 * layer interfaces below it (u* from the bulk drag on the lowest layer's
 * wind, with the gustiness floor). The interface coefficients ρK/Δz are
 * kept in a shared array so the cell units of the adjust phase can mix
 * θ, q and qc down each column and the edge units can mix the normal
 * velocity down each edge, both by implicit Euler on the same
 * tridiagonal system, which conserves each column's mass-weighted
 * total exactly. Nothing mixes above the boundary-layer top; the search
 * stops at searchTop in σ. The surface fluxes and drag remain explicit
 * sources on the lowest layer, which the diffusion then spreads upward.
 */
export function createBoundaryLayer(mesh, core, {
  dragCoefficient = 1.5e-3, dragCoefficients = null, gustiness = 3, richardsonCritical = 0.5, vonKarman = 0.4, searchTop = 0.5, buffers = null,
} = {}) {
  const { K, C, E, dSigma, sigmaMid, R, g, exnerLayer, geopotential } = core.diagnostics;
  const thetaV = core.arrays.thetaV;
  const { cellsOnEdge } = mesh;
  const bottom = K - 1;
  let kTop = 0;
  while (kTop < bottom && sigmaMid[kTop] <= searchTop) kTop++;
  const n = K - kTop;
  const mixingBuffer = buffers && buffers.mixing ? buffers.mixing : new SharedArrayBuffer(8 * K * C);
  const depthBuffer = buffers && buffers.depth ? buffers.depth : new SharedArrayBuffer(8 * C);
  const mixing = new Float64Array(mixingBuffer);
  const depth = new Float64Array(depthBuffer);
  const vector = new Float64Array(3 * C), bottomVector = new Float64Array(3 * C);
  const speed = new Float64Array(C), friction = new Float64Array(C), riPrev = new Float64Array(C), zPrev = new Float64Array(C);
  const found = new Uint8Array(C);
  const upper = new Float64Array(K), lower = new Float64Array(K), gain = new Float64Array(K), rhs = new Float64Array(K), mass = new Float64Array(K);

  function diagnose(state, iFrom = 0, iTo = C) {
    const [pi, theta, u, , q = null, qc = null] = state;
    for (let i = iFrom; i < iTo; i++) core.diagnoseColumn(i, pi, theta, q, qc);
    cellVector(mesh, u.subarray(bottom * E, K * E), bottomVector, iFrom, iTo);
    for (let i = iFrom; i < iTo; i++) {
      speed[i] = Math.hypot(bottomVector[3 * i], bottomVector[3 * i + 1], bottomVector[3 * i + 2]);
      friction[i] = Math.sqrt(dragCoefficients ? dragCoefficients[i] : dragCoefficient) * Math.max(speed[i], gustiness);
      found[i] = 0;
      riPrev[i] = 0;
      zPrev[i] = geopotential[bottom * C + i] / g;
      depth[i] = zPrev[i];
    }
    for (let k = bottom - 1; k >= kTop; k--) {
      cellVector(mesh, u.subarray(k * E, (k + 1) * E), vector, iFrom, iTo);
      for (let i = iFrom; i < iTo; i++) {
        if (found[i]) continue;
        const idx = k * C + i, base = bottom * C + i;
        const z = geopotential[idx] / g, zb = geopotential[base] / g;
        const du = vector[3 * i] - bottomVector[3 * i], dv = vector[3 * i + 1] - bottomVector[3 * i + 1], dw = vector[3 * i + 2] - bottomVector[3 * i + 2];
        const shear = du * du + dv * dv + dw * dw + 100 * friction[i] * friction[i];
        const ri = g * (thetaV[idx] - thetaV[base]) * (z - zb) / (thetaV[base] * shear);
        if (ri > richardsonCritical) {
          depth[i] = zPrev[i] + (z - zPrev[i]) * (richardsonCritical - riPrev[i]) / (ri - riPrev[i]);
          found[i] = 1;
        } else {
          riPrev[i] = ri;
          zPrev[i] = z;
          if (k === kTop) depth[i] = z;
        }
      }
    }
    for (let i = iFrom; i < iTo; i++) {
      const zb = geopotential[bottom * C + i] / g, h = depth[i] - zb;
      for (let k = kTop; k < K; k++) mixing[k * C + i] = 0;
      if (h <= 0) continue;
      for (let k = kTop; k < bottom; k++) {
        const idx = k * C + i, below = idx + C;
        const zAbove = geopotential[idx] / g, zBelow = geopotential[below] / g;
        const z = 0.5 * (zAbove + zBelow) - zb;
        if (z >= h) continue;
        const diffusivity = vonKarman * friction[i] * z * (1 - z / h) ** 2;
        const rhoAbove = pi[i] * sigmaMid[k] / (R * theta[idx] * exnerLayer[idx]);
        const rhoBelow = pi[i] * sigmaMid[k + 1] / (R * theta[below] * exnerLayer[below]);
        mixing[idx] = 0.5 * (rhoAbove + rhoBelow) * diffusivity / (zAbove - zBelow);
      }
    }
  }

  function solve(field, offset, stride, coefficient, coefficientStride, dt, columnMass) {
    let active = false;
    for (let k = kTop; k < bottom; k++) if (coefficient[k * coefficientStride] > 0) { active = true; break; }
    if (!active) return false;
    for (let j = 0; j < n; j++) {
      const k = kTop + j;
      mass[j] = columnMass * dSigma[k] / g;
      upper[j] = j > 0 ? dt * coefficient[(k - 1) * coefficientStride] / mass[j] : 0;
      lower[j] = k < bottom ? dt * coefficient[k * coefficientStride] / mass[j] : 0;
      rhs[j] = field[offset + k * stride];
    }
    let denominator = 1 + upper[0] + lower[0];
    gain[0] = -lower[0] / denominator;
    rhs[0] /= denominator;
    for (let j = 1; j < n; j++) {
      denominator = 1 + upper[j] + lower[j] + upper[j] * gain[j - 1];
      gain[j] = -lower[j] / denominator;
      rhs[j] = (rhs[j] + upper[j] * rhs[j - 1]) / denominator;
    }
    field[offset + (kTop + n - 1) * stride] = rhs[n - 1];
    for (let j = n - 2; j >= 0; j--) {
      rhs[j] -= gain[j] * rhs[j + 1];
      field[offset + (kTop + j) * stride] = rhs[j];
    }
    return true;
  }

  function mixColumn(i, pi, theta, q, qc, dt) {
    const coefficient = mixing.subarray(i);
    if (!solve(theta, i, C, coefficient, C, dt, pi[i])) return;
    if (q) solve(q, i, C, coefficient, C, dt, pi[i]);
    if (qc) solve(qc, i, C, coefficient, C, dt, pi[i]);
  }

  const edgeCoefficient = new Float64Array(K), before = new Float64Array(K), share = new Float64Array(K);
  /*
   * The implicit mixing removes kinetic energy at each interface's shear
   * and in each layer's own increment; `dissipation` receives each
   * layer's loss of u² in those proportions, scaled so that the column's
   * mass-weighted loss is exact.
   */
  function mixEdges(pi, u, eFrom, eTo, dt, dissipation = null) {
    for (let e = eFrom; e < eTo; e++) {
      const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1], columnMass = 0.5 * (pi[a] + pi[b]);
      for (let k = kTop; k < K; k++) { edgeCoefficient[k] = 0.5 * (mixing[k * C + a] + mixing[k * C + b]); before[k] = u[k * E + e]; }
      if (!solve(u, e, E, edgeCoefficient, 1, dt, columnMass) || !dissipation) continue;
      let loss = 0, total = 0;
      for (let k = kTop; k < K; k++) {
        const m = columnMass * dSigma[k] / g, now = u[k * E + e], change = now - before[k];
        loss += m * (before[k] * before[k] - now * now);
        share[k] = m * change * change;
      }
      for (let k = kTop; k < bottom; k++) {
        const shear = u[k * E + e] - u[(k + 1) * E + e], part = dt * edgeCoefficient[k] * shear * shear;
        share[k] += part; share[k + 1] += part;
      }
      for (let k = kTop; k < K; k++) total += share[k];
      if (total <= 0) continue;
      for (let k = kTop; k < K; k++) dissipation[k * E + e] += loss * share[k] / (total * columnMass * dSigma[k] / g);
    }
  }

  return { diagnose, mixColumn, mixEdges, mixing, depth, kTop, shared: { mixing: mixingBuffer, depth: depthBuffer } };
}
