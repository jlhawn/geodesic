import { cellVector } from '../dynamics/operators.module.js';
import { saturationHumidity, DECK_CLOSED, ACTIVITY_UNDECIDED } from './moist.module.js';

/*
 * A diffusive planetary boundary layer in the manner of Troen and Mahrt
 * (1986), as used by the simple moist GCMs that this model follows.
 * `diagnose` finds, per cell, the boundary-layer top as the height where
 * the bulk Richardson number of the lowest layer's virtual potential
 * temperature and wind, with the convective floor 100 u*², first exceeds
 * richardsonCritical, and lays the K-profile κ u* z (1 − z/h)² over the
 * layer interfaces below it (u* from the bulk drag on the lowest layer's
 * wind, with the gustiness floor). Where the surface is warmer than the
 * lowest layer the profile's velocity scale is the unstable one of
 * Holtslag and Boville (1993), u* (1 − 15 ζ)^¼ with ζ = 0.1 h / L from
 * the bulk surface buoyancy flux (virtual, with the saturation humidity
 * of a sea surface; dry over `land`) and floored at −2, so a convective
 * marine boundary layer mixes momentum down to the surface; stable
 * columns keep the neutral profile. The interface coefficients ρK/Δz are
 * kept in a shared array so the cell units of the adjust phase can mix
 * θ, q and qc down each column and the edge units can mix the normal
 * velocity down each edge, both by implicit Euler on the same
 * tridiagonal system, which conserves each column's mass-weighted
 * total exactly. Nothing mixes above the boundary-layer top; the search
 * stops at searchTop in σ. The surface fluxes and drag remain explicit
 * sources on the lowest layer, which the diffusion then spreads upward.
 *
 * A stratocumulus deck mixes its layer from cloud top, which the bulk
 * Richardson number of the surface-driven search does not see. With
 * `deckTop`, per cell the inversion height of the mixed-layer deck in the
 * height coordinate of `depth` where the deck ran this step and 0 where
 * it did not (radiation.module.js's mlmTop), the K-profile of such a
 * cell is laid over max(depth, deckTop) instead of depth — the same
 * profile, shape and velocity scale, for the deeper layer — so the
 * deck's layer is mixed through to its inversion. `depth` itself stays
 * the Richardson depth: the deck starts from it and relaxes toward it.
 *
 * The K-profile vanishes at h, so its top entrains nothing; an explicit
 * entrainment flux closes it (`entrainment`). Where the surface buoyancy
 * flux B0 (the one above, positive upward) is positive,
 *   w_e = o (1 − s) min(cap, (A B0 + A_s r u*³ / h) / max(Δb, bMin)),
 * Δb = g Δθv / θv between the first layer whose base lies above h and
 * the mass mean of the layers below it: the buoyancy and friction-velocity
 * sources of Tennekes's (1973) inversion model with Driedonks's (1982)
 * constants A 0.2 and A_s 5, bMin 0.015 m/s² (a jump of about 0.5 K) and
 * the cap 0.05 m/s. r = min(1, B0 / `shearOnset`)
 * (5·10⁻⁵ m²/s³, about 2 W/m² of virtual heat flux) brings the shear
 * term in continuously as the surface turns unstable. The stratocumulus
 * regime belongs to the mixed-layer deck: o is the deck's opening, 1 at a
 * gate `deckGate` (radiation.mlmGate) of one half or less and 0 at
 * DECK_CLOSED, as convection takes it, and s the radiation's stratiform
 * share `stratiform` (radiation.stratiform), 0 below an estimated
 * inversion strength of 8 K and 1 above 12 K. It enters as the
 * coefficient ρ w_e of that layer's lower interface, so the exchange of
 * θ, q, qc and momentum across h goes through the same conservative
 * implicit solve; `entrainment` keeps each cell's w_e (m/s). The depth
 * stays diagnostic.
 */
export const ENTRAINMENT_DEFAULTS = { efficiency: 0.2, shear: 5, cap: 0.05, jumpFloor: 0.015, shearOnset: 5e-5 };

export function createBoundaryLayer(mesh, core, {
  dragCoefficient = 1.5e-3, dragCoefficients = null, gustiness = 3, richardsonCritical = 0.5, vonKarman = 0.4, searchTop = 0.5, stability = true, land = null, deckTop = null, deckGate = null, stratiform = null, buffers = null,
  entrainment: entrainmentOptions = {},
} = {}) {
  const { efficiency, shear, cap, jumpFloor, shearOnset } = { ...ENTRAINMENT_DEFAULTS, ...entrainmentOptions };
  const { K, C, E, dSigma, sigmaMid, R, g, kappa, exnerLayer, geopotential } = core.diagnostics;
  const thetaV = core.arrays.thetaV;
  const { cellsOnEdge } = mesh;
  const bottom = K - 1;
  let kTop = 0;
  while (kTop < bottom && sigmaMid[kTop] <= searchTop) kTop++;
  const n = K - kTop;
  const mixingBuffer = buffers && buffers.mixing ? buffers.mixing : new SharedArrayBuffer(8 * K * C);
  const depthBuffer = buffers && buffers.depth ? buffers.depth : new SharedArrayBuffer(8 * C);
  const buoyancyBuffer = buffers && buffers.buoyancyFlux ? buffers.buoyancyFlux : new SharedArrayBuffer(8 * C);
  const frictionBuffer = buffers && buffers.friction ? buffers.friction : new SharedArrayBuffer(8 * C);
  const entrainmentBuffer = buffers && buffers.entrainment ? buffers.entrainment : new SharedArrayBuffer(8 * C);
  const mixing = new Float64Array(mixingBuffer);
  const depth = new Float64Array(depthBuffer);
  const buoyancyFlux = new Float64Array(buoyancyBuffer), friction = new Float64Array(frictionBuffer);
  const entrainmentVelocity = new Float64Array(entrainmentBuffer);
  const entraining = efficiency > 0 || shear > 0;
  const vector = new Float64Array(3 * C), bottomVector = new Float64Array(3 * C);
  const speed = new Float64Array(C), riPrev = new Float64Array(C), zPrev = new Float64Array(C);
  const found = new Uint8Array(C);
  const upper = new Float64Array(K), lower = new Float64Array(K), gain = new Float64Array(K), rhs = new Float64Array(K), mass = new Float64Array(K);

  function diagnose(state, iFrom = 0, iTo = C) {
    const [pi, theta, u, surfaceT, q = null, qc = null] = state;
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
      const zb = geopotential[bottom * C + i] / g, h = (deckTop && deckTop[i] > 0 ? Math.max(depth[i], deckTop[i]) : depth[i]) - zb;
      for (let k = kTop; k < K; k++) mixing[k * C + i] = 0;
      entrainmentVelocity[i] = 0;
      const base = bottom * C + i;
      const moisture = q && !(land && land[i]) ? 0.61 * theta[base] * (saturationHumidity(surfaceT[i], pi[i]) - q[base]) : 0;
      const buoyancy = g / theta[base] * (dragCoefficients ? dragCoefficients[i] : dragCoefficient) * Math.max(speed[i], gustiness) * (surfaceT[i] * Math.pow(sigmaMid[bottom], kappa) / exnerLayer[base] - theta[base] + moisture);
      buoyancyFlux[i] = buoyancy;
      if (h <= 0) continue;
      let scale = friction[i];
      if (stability && buoyancy > 0) scale = friction[i] * Math.pow(1 - 15 * Math.max(-2, -0.1 * h * vonKarman * buoyancy / friction[i] ** 3), 0.25);
      let entrainK = -1;
      for (let k = kTop; k < bottom; k++) {
        const idx = k * C + i, below = idx + C;
        const zAbove = geopotential[idx] / g, zBelow = geopotential[below] / g;
        const z = 0.5 * (zAbove + zBelow) - zb;
        if (z >= h) { entrainK = k; continue; }
        const diffusivity = vonKarman * scale * z * (1 - z / h) ** 2;
        const rhoAbove = pi[i] * sigmaMid[k] / (R * theta[idx] * exnerLayer[idx]);
        const rhoBelow = pi[i] * sigmaMid[k + 1] / (R * theta[below] * exnerLayer[below]);
        mixing[idx] = 0.5 * (rhoAbove + rhoBelow) * diffusivity / (zAbove - zBelow);
      }
      const open = deckGate ? Math.min(1, Math.max(0, (DECK_CLOSED - deckGate[i]) / (DECK_CLOSED - ACTIVITY_UNDECIDED))) : 1;
      const taper = open * (stratiform ? 1 - stratiform[i] : 1);
      if (entraining && entrainK >= kTop && buoyancy > 0 && taper > 0) {
        const idx = entrainK * C + i, below = idx + C;
        let weight = 0, sum = 0;
        for (let k = entrainK + 1; k < K; k++) { weight += dSigma[k]; sum += dSigma[k] * thetaV[k * C + i]; }
        const mean = sum / weight, jump = g * (thetaV[idx] - mean) / mean;
        const onset = shearOnset > 0 ? Math.min(1, buoyancy / shearOnset) : 1;
        const velocity = taper * Math.min(cap, (efficiency * buoyancy + shear * onset * friction[i] ** 3 / h) / Math.max(jump, jumpFloor));
        entrainmentVelocity[i] = velocity;
        mixing[idx] = 0.5 * (pi[i] * sigmaMid[entrainK] / (R * theta[idx] * exnerLayer[idx]) + pi[i] * sigmaMid[entrainK + 1] / (R * theta[below] * exnerLayer[below])) * velocity;
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

  return { diagnose, mixColumn, mixEdges, mixing, depth, buoyancyFlux, friction, entrainment: entrainmentVelocity, kTop, shared: { mixing: mixingBuffer, depth: depthBuffer, buoyancyFlux: buoyancyBuffer, friction: frictionBuffer, entrainment: entrainmentBuffer } };
}
