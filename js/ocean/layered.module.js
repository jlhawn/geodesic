import { divergence, gradient, curl, kineticEnergy, laplacianVelocity, cellVector } from '../dynamics/operators.module.js';
import { FREEZING_POINT } from '../physics/ice.module.js';

/*
 * A layered ocean on the C-grid: a bulk mixed layer with its own
 * temperature and salinity over interior layers of fixed density, the
 * hybrid isopycnal design of MICOM, on the real bathymetry. Each layer
 * is a TRiSK shallow-water layer in the vector-invariant form, carrying
 * thickness, edge velocity, heat h·T and salt h·S. The pressure force
 * in an interior layer k is the gradient of a cell potential,
 *   Φ_k = g η + (g/ρ₀)[ρ_ml h₀ − ρ_k h₀ + Σ_{j<k} (ρ_j − ρ_k) h_j],
 * which is exact on the mesh whatever the bathymetry; the mixed layer's
 * force is g∇η plus the depth-mean of its own density gradient. A layer
 * that has outcropped, or lies below the bottom, keeps a token thickness
 * and follows the velocity of the layer above.
 *
 * The free surface η = Σh − D moves at √(gD), too fast for the ocean's
 * step, so the depth-integrated flow is sub-stepped: the baroclinic
 * layers take one RK4 step with η frozen at its starting value, the
 * barotropic transport and η take RK4 sub-steps within the same
 * interval forced by the depth integral of the slow tendencies, and the
 * layers are then rescaled to the sub-stepped η and shifted to the
 * averaged transport.
 *
 * The mixed layer exchanges mass with the interior after each step: it
 * swallows any interior layer lighter than itself (convection, within
 * `maximumMixedDepth`, and keeps at most that depth: the excess returns
 * to the interior each step, so the column keeps mixing while the
 * mixed layer itself stays a surface layer), entrains the layer below at the Kraus–Turner
 * wind-stirring rate, the stirring fading with depth over
 * `stirringDepth`, and detrains when it is deeper than the maximum,
 * when it is as dense as the water beneath it (convectively neutral),
 * or when surface warming makes it deeper than the Monin–Obukhov depth.
 * Detrained water goes to the interior layer whose density is nearest
 * its own, so water swallowed from a layer returns to that layer; a
 * mass-conserving split between the two bracketing layers instead
 * ratcheted a fraction of every swallow-and-return cycle into the
 * denser class. Tracers are carried by the flux with the donor cell's
 * value. Surface freshwater (evaporation
 * minus rain, runoff spread over the sea) and ice growth or melt act on
 * its salinity as virtual salt fluxes. As before, the mixed layer's
 * temperature is the sea surface temperature the atmosphere sees, its
 * heat capacity is published per cell, and the heat converged under ice
 * is handed to the ice base.
 */
export const LAYER_DENSITIES = [1024.0, 1025.5, 1026.5, 1027.2, 1027.7];
export const LAYER_BOTTOMS = [250, 600, 1200, 2500];
const EPS = 0.01, THIN = 5, PV_FLOOR = 20, SPEED_LIMIT = 5, DENSITY_TOLERANCE = 0.005;

/*
 * The model's bathymetry: the cell-mean ETOPO depth of every sea cell,
 * at least `minimumDepth`, and never shallower than `neighbourRatio`
 * times its deepest sea neighbour, so that no shelf break or trench
 * wall drops by more than that ratio across one edge. A cell 110 km
 * across cannot hold a real shelf, and a ten-to-one step in depth
 * between two cells makes the free surface swing wildly at the coast.
 */
export function bathymetryFrom(mesh, geography, { minimumDepth = 50, neighbourRatio = 0.5, flatDepth = 4000 } = {}) {
  const { nCells: C, maxEdges, nEdgesOnCell, cellsOnCell } = mesh;
  const D = new Float64Array(C);
  for (let i = 0; i < C; i++) D[i] = geography ? (geography.land[i] ? 0 : Math.max(minimumDepth, -geography.elevation[i])) : flatDepth;
  if (!geography) return D;
  for (let sweep = 0; sweep < 50; sweep++) {
    let changed = false;
    for (let i = 0; i < C; i++) {
      if (geography.land[i]) continue;
      let deepest = 0;
      for (let m = 0; m < nEdgesOnCell[i]; m++) { const j = cellsOnCell[maxEdges * i + m]; if (!geography.land[j]) deepest = Math.max(deepest, D[j]); }
      const floor = neighbourRatio * deepest;
      if (D[i] < floor) { D[i] = floor; changed = true; }
    }
    if (!changed) break;
  }
  return D;
}

export function createOcean(mesh, {
  densities = LAYER_DENSITIES, bottoms = LAYER_BOTTOMS, mixedDepth = 60, minimumDepth = 50, flatDepth = 4000, thermoclineTilt = 0.3,
  salinityProfile = (lat) => 34.5 + 1.5 * Math.exp(-(((Math.abs(lat) * 180 / Math.PI - 25) / 15) ** 2)),
  density = 1025, specificHeat = 3985, thermalExpansion = 2e-4, halineContraction = 7.6e-4, referenceT = 283.15, referenceS = 35, gravity = 9.81,
  minimumThickness = 20, shallowestMixedDepth = 50, maximumMixedDepth = 200, stirring = 0.8, stirringDepth = 100, detrainmentTime = 86400, iceSalinity = 5, iceDensity = 917,
  interfacialDrag = 2e-4, bottomDrag = 3e-3, closureHours = 12, diffusivity = 0.3, everySteps = 4,
  geography = null, bathymetry = null, buffers = null,
} = {}) {
  const {
    nCells: C, nEdges: E, nVertices: V, maxEdgesOnEdge, nEdgesOnEdge, edgesOnEdge, weightsOnEdge,
    cellsOnEdge, verticesOnEdge, areaCell, dcEdge, dvEdge, fVertex, fEdge, radius, latCell,
    maxEdges, nEdgesOnCell, edgesOnCell, cellsOnCell,
  } = mesh;
  const L = densities.length + 1, g = gravity, rho0 = density, rhoCp = density * specificHeat;
  const rho = [rho0, ...densities];
  const labelT = rho.map((r) => Math.max(FREEZING_POINT, referenceT - (r / rho0 - 1) / thermalExpansion));
  const diffusion = diffusivity * radius * radius / rhoCp;
  let spacing = 0, minSpacing = Infinity;
  for (let e = 0; e < E; e++) { spacing += dcEdge[e]; minSpacing = Math.min(minSpacing, dcEdge[e]); }
  spacing /= E;
  const nu4 = closureHours > 0 ? Math.pow(spacing / Math.PI, 4) / (closureHours * 3600) : 0;

  const edgeOcean = geography ? geography.edgeOcean : new Uint8Array(E).fill(1);
  const cellOcean = geography ? Uint8Array.from(geography.land, (l) => (l ? 0 : 1)) : new Uint8Array(C).fill(1);
  const D = bathymetry ? Float64Array.from(bathymetry, (d, i) => (cellOcean[i] ? d : 0)) : bathymetryFrom(mesh, geography, { minimumDepth, flatDepth });
  let deepest = 0;
  for (let i = 0; i < C; i++) deepest = Math.max(deepest, D[i]);
  const substepLimit = 0.35 * minSpacing / Math.sqrt(g * Math.max(deepest, 1));

  const h = new Float64Array(L * C), u = new Float64Array(L * E), Q = new Float64Array(L * C), W = new Float64Array(L * C);
  const state = [h, u, Q, W];
  const eta = new Float64Array(C);
  const T0 = new Float64Array(C), S0 = new Float64Array(C), rhoMl = new Float64Array(C), previousT0 = new Float64Array(C), surfaceIn = new Float64Array(C), previousIce = new Float64Array(C);
  const capacity = new Float64Array(buffers && buffers.capacity ? buffers.capacity : new SharedArrayBuffer(8 * C)).fill(rhoCp * mixedDepth);
  const stress = new Float64Array(E), fresh = new Float64Array(C);
  const iced = new Uint8Array(C);
  const hEdge = new Float64Array(L * E), flux = new Float64Array(E), tracerFlux = new Float64Array(E);
  const T = new Float64Array(C), S = new Float64Array(C), lapT = new Float64Array(C);
  const zeta = new Float64Array(V), qEdge = new Float64Array(E);
  const K = new Float64Array(C), phi = new Float64Array(C), gradPhi = new Float64Array(E), gradEta = new Float64Array(E), gradRho = new Float64Array(E);
  const lap = new Float64Array(E), lap2 = new Float64Array(E), divScratch = new Float64Array(C), curlScratch = new Float64Array(V);
  const slow = new Float64Array(E), U = new Float64Array(E), depthEdge = new Float64Array(E), etaB = new Float64Array(C), avgU = new Float64Array(E), avgEta = new Float64Array(C), divU = new Float64Array(C);
  const stages = [0, 1, 2, 3].map(() => state.map((a) => new Float64Array(a.length)));
  const trial = state.map((a) => new Float64Array(a.length));
  const tauCell = new Float64Array(3 * C);
  let counter = 0, limited = 0, relaxRate = 1 / 3600, initialised = false;

  const eos = (t, s) => rho0 * (1 - thermalExpansion * (t - referenceT) + halineContraction * (s - referenceS));
  const at = (k, i) => k * C + i;
  const ae = (k, e) => k * E + e;

  function maskedLaplacian(field, out) {
    for (let i = 0; i < C; i++) {
      let sum = 0;
      for (let m = 0; m < nEdgesOnCell[i]; m++) {
        const e = edgesOnCell[maxEdges * i + m];
        if (edgeOcean[e]) sum += dvEdge[e] * (field[cellsOnCell[maxEdges * i + m]] - field[i]) / dcEdge[e];
      }
      out[i] = sum / areaCell[i];
    }
  }
  /*
   * The mixed layer, present everywhere, takes the centred thickness at
   * an edge; an interior layer takes the smaller of the two, so it never
   * flows into a cell where it has no water, whether the layer has
   * outcropped there or the bottom lies above it. A column can only
   * flow through the water that exists on both sides of an edge, so at a
   * shelf break or a coast the thicknesses are scaled to the shallower
   * side's depth.
   */
  function edgeThicknesses(hIn) {
    for (let e = 0; e < E; e++) {
      const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
      hEdge[e] = 0.5 * (hIn[a] + hIn[b]);
      for (let k = 1; k < L; k++) hEdge[ae(k, e)] = Math.min(hIn[at(k, a)], hIn[at(k, b)]);
    }
    for (let e = 0; e < E; e++) {
      const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
      const sill = Math.max(EPS, Math.min(D[a], D[b]) + 0.5 * (eta[a] + eta[b]));
      let sum = 0;
      for (let k = 0; k < L; k++) sum += hEdge[ae(k, e)];
      if (sum > sill) { const f = sill / sum; for (let k = 0; k < L; k++) hEdge[ae(k, e)] *= f; }
    }
  }
  function coriolisTransport(Ue, out) {
    for (let e = 0; e < E; e++) {
      let sum = 0;
      for (let s = 0; s < nEdgesOnEdge[e]; s++) {
        const other = edgesOnEdge[maxEdgesOnEdge * e + s];
        sum += weightsOnEdge[maxEdgesOnEdge * e + s] * dvEdge[other] * Ue[other] * 0.5 * (fEdge[e] + fEdge[other]);
      }
      out[e] = sum / dcEdge[e];
    }
  }

  function surfaceDensity(hIn, QIn, WIn) {
    for (let i = 0; i < C; i++) {
      const h0 = Math.max(EPS, hIn[i]);
      rhoMl[i] = eos(QIn[i] / h0, WIn[i] / h0);
    }
  }

  /*
   * One tendency evaluation. Interior layer potentials telescope through
   * the layers above, so a token-thickness layer contributes nothing and
   * the layers it separates feel the density step across it.
   */
  function freeSurface() {
    for (let i = 0; i < C; i++) {
      let sum = 0;
      for (let k = 0; k < L; k++) sum += h[at(k, i)];
      eta[i] = cellOcean[i] ? sum - D[i] : 0;
    }
  }

  function tendency(input, out) {
    const [hIn, uIn, QIn, WIn] = input;
    const [dh, du, dQ, dW] = out;
    surfaceDensity(hIn, QIn, WIn);
    gradient(mesh, eta, gradEta);
    gradient(mesh, rhoMl, gradRho);
    edgeThicknesses(hIn);
    for (let k = 0; k < L; k++) {
      const oc = k * C, oe = k * E;
      for (let e = 0; e < E; e++) {
        let he = hEdge[oe + e];
        if (k === 0) he = Math.min(he, Math.max(0, hIn[cellsOnEdge[2 * e + (uIn[e] > 0 ? 0 : 1)]]));
        flux[e] = edgeOcean[e] ? he * uIn[oe + e] : 0;
      }
      divergence(mesh, flux, divScratch);
      for (let i = 0; i < C; i++) {
        dh[oc + i] = -divScratch[i];
        const hh = Math.max(EPS, hIn[oc + i]);
        T[i] = QIn[oc + i] / hh; S[i] = WIn[oc + i] / hh;
      }
      for (let e = 0; e < E; e++) tracerFlux[e] = flux[e] * T[cellsOnEdge[2 * e + (flux[e] > 0 ? 0 : 1)]];
      divergence(mesh, tracerFlux, divScratch);
      for (let i = 0; i < C; i++) dQ[oc + i] = -divScratch[i];
      for (let e = 0; e < E; e++) tracerFlux[e] = flux[e] * S[cellsOnEdge[2 * e + (flux[e] > 0 ? 0 : 1)]];
      divergence(mesh, tracerFlux, divScratch);
      for (let i = 0; i < C; i++) dW[oc + i] = -divScratch[i];
      if (k === 0 && diffusion > 0) {
        maskedLaplacian(T, lapT);
        for (let i = 0; i < C; i++) dQ[i] += diffusion * lapT[i];
        maskedLaplacian(S, lapT);
        for (let i = 0; i < C; i++) dW[i] += diffusion * lapT[i];
      }
      curl(mesh, uIn.subarray(oe, oe + E), zeta);
      for (let e = 0; e < E; e++) qEdge[e] = 0.5 * (zeta[verticesOnEdge[2 * e]] + fVertex[verticesOnEdge[2 * e]] + zeta[verticesOnEdge[2 * e + 1]] + fVertex[verticesOnEdge[2 * e + 1]]) / Math.max(hEdge[oe + e], PV_FLOOR);
      kineticEnergy(mesh, uIn.subarray(oe, oe + E), K);
      if (k === 0) {
        for (let i = 0; i < C; i++) phi[i] = K[i] + g * eta[i];
      } else {
        for (let i = 0; i < C; i++) {
          let p = (rhoMl[i] - rho[k]) * hIn[i];
          for (let j = 1; j < k; j++) p += (rho[j] - rho[k]) * hIn[at(j, i)];
          phi[i] = K[i] + g * eta[i] + g * p / rho0;
        }
      }
      gradient(mesh, phi, gradPhi);
      for (let e = 0; e < E; e++) {
        let sum = 0;
        for (let s = 0; s < nEdgesOnEdge[e]; s++) {
          const other = edgesOnEdge[maxEdgesOnEdge * e + s];
          sum += weightsOnEdge[maxEdgesOnEdge * e + s] * dvEdge[other] * flux[other] * 0.5 * (qEdge[e] + qEdge[other]);
        }
        du[oe + e] = sum / dcEdge[e] - gradPhi[e];
      }
      if (k === 0) for (let e = 0; e < E; e++) du[e] -= g / rho0 * 0.5 * hEdge[e] * gradRho[e];
      for (let e = 0; e < E; e++) {
        const he = Math.max(hEdge[oe + e], minimumThickness);
        let force = 0;
        if (k === 0) force += stress[e] / rho0;
        if (k > 0) force += interfacialDrag * (uIn[ae(k - 1, e)] - uIn[oe + e]);
        if (k < L - 1) force -= interfacialDrag * (uIn[oe + e] - uIn[ae(k + 1, e)]);
        let bottom = k === L - 1;
        if (!bottom) { bottom = true; for (let j = k + 1; j < L; j++) if (hEdge[ae(j, e)] >= THIN) { bottom = false; break; } }
        if (bottom) force -= bottomDrag * Math.abs(uIn[oe + e]) * uIn[oe + e];
        du[oe + e] += force / he;
      }
      if (nu4 > 0) {
        laplacianVelocity(mesh, uIn.subarray(oe, oe + E), lap, divScratch, curlScratch);
        laplacianVelocity(mesh, lap, lap2, divScratch, curlScratch);
        for (let e = 0; e < E; e++) du[oe + e] -= nu4 * lap2[e];
      }
      if (k > 0) for (let e = 0; e < E; e++) if (hEdge[oe + e] < THIN) du[oe + e] = (uIn[ae(k - 1, e)] - uIn[oe + e]) * relaxRate;
      for (let e = 0; e < E; e++) if (!edgeOcean[e]) du[oe + e] = 0;
    }
    for (let i = 0; i < C; i++) if (!cellOcean[i]) for (let k = 0; k < L; k++) { dh[at(k, i)] = 0; dQ[at(k, i)] = 0; dW[at(k, i)] = 0; }
  }

  function combine(target, base, stage, dt) {
    for (let a = 0; a < state.length; a++) {
      const t = target[a], b = base[a], s = stage[a];
      for (let i = 0; i < t.length; i++) t[i] = b[i] + dt * s[i];
    }
  }

  function barotropic(dt) {
    const [k1] = stages;
    const [, du] = k1;
    for (let e = 0; e < E; e++) {
      let sumH = 0, transport = 0, forcing = 0;
      for (let k = 0; k < L; k++) { const he = hEdge[ae(k, e)]; sumH += he; transport += he * u[ae(k, e)]; forcing += he * (du[ae(k, e)] + g * gradEta[e]); }
      U[e] = edgeOcean[e] ? transport : 0;
      slow[e] = edgeOcean[e] ? forcing : 0;
      depthEdge[e] = sumH;
    }
    coriolisTransport(U, lap);
    for (let e = 0; e < E; e++) slow[e] -= lap[e];
    const M = Math.max(1, Math.ceil(dt / substepLimit)), dtb = dt / M;
    etaB.set(eta); avgU.fill(0); avgEta.fill(0);
    for (let m = 0; m < M; m++) {
      rk4Barotropic(dtb);
      for (let i = 0; i < C; i++) avgEta[i] += etaB[i] / M;
      for (let e = 0; e < E; e++) avgU[e] += U[e] / M;
    }
  }

  const bStages = [0, 1, 2, 3].map(() => [new Float64Array(C), new Float64Array(E)]);
  const bTrial = [new Float64Array(C), new Float64Array(E)];
  function barotropicTendency(etaIn, uIn, dEta, dU) {
    gradient(mesh, etaIn, gradPhi);
    coriolisTransport(uIn, lap);
    for (let e = 0; e < E; e++) {
      if (!edgeOcean[e]) { dU[e] = 0; continue; }
      dU[e] = -g * depthEdge[e] * gradPhi[e] + lap[e] + slow[e];
    }
    divergence(mesh, uIn, divU);
    for (let i = 0; i < C; i++) dEta[i] = cellOcean[i] ? -divU[i] : 0;
  }
  function rk4Barotropic(dtb) {
    const [s1, s2, s3, s4] = bStages, [tEta, tU] = bTrial;
    barotropicTendency(etaB, U, s1[0], s1[1]);
    for (let i = 0; i < C; i++) tEta[i] = etaB[i] + 0.5 * dtb * s1[0][i];
    for (let e = 0; e < E; e++) tU[e] = U[e] + 0.5 * dtb * s1[1][e];
    barotropicTendency(tEta, tU, s2[0], s2[1]);
    for (let i = 0; i < C; i++) tEta[i] = etaB[i] + 0.5 * dtb * s2[0][i];
    for (let e = 0; e < E; e++) tU[e] = U[e] + 0.5 * dtb * s2[1][e];
    barotropicTendency(tEta, tU, s3[0], s3[1]);
    for (let i = 0; i < C; i++) tEta[i] = etaB[i] + dtb * s3[0][i];
    for (let e = 0; e < E; e++) tU[e] = U[e] + dtb * s3[1][e];
    barotropicTendency(tEta, tU, s4[0], s4[1]);
    for (let i = 0; i < C; i++) etaB[i] += dtb / 6 * (s1[0][i] + 2 * s2[0][i] + 2 * s3[0][i] + s4[0][i]);
    for (let e = 0; e < E; e++) U[e] += dtb / 6 * (s1[1][e] + 2 * s2[1][e] + 2 * s3[1][e] + s4[1][e]);
  }

  function step(dt) {
    const [k1, k2, k3, k4] = stages;
    freeSurface();
    tendency(state, k1);
    barotropic(dt);
    combine(trial, state, k1, dt / 2);
    tendency(trial, k2);
    combine(trial, state, k2, dt / 2);
    tendency(trial, k3);
    combine(trial, state, k3, dt);
    tendency(trial, k4);
    const w = dt / 6;
    for (let a = 0; a < state.length; a++) {
      const s = state[a], a1 = k1[a], a2 = k2[a], a3 = k3[a], a4 = k4[a];
      for (let i = 0; i < s.length; i++) s[i] += w * (a1[i] + 2 * a2[i] + 2 * a3[i] + a4[i]);
    }
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) continue;
      let sum = 0;
      for (let k = 0; k < L; k++) {
        const n = at(k, i);
        if (h[n] < EPS) { const t = h[n] > 1e-9 ? Q[n] / h[n] : labelT[k], s = h[n] > 1e-9 ? W[n] / h[n] : referenceS; h[n] = EPS; Q[n] = EPS * t; W[n] = EPS * s; }
        sum += h[n];
      }
      const scale = (D[i] + avgEta[i]) / sum;
      for (let k = 0; k < L; k++) { const n = at(k, i); h[n] *= scale; Q[n] *= scale; W[n] *= scale; }
      eta[i] = avgEta[i];
    }
    edgeThicknesses(h);
    for (let e = 0; e < E; e++) {
      if (!edgeOcean[e]) continue;
      let sumH = 0, transport = 0;
      for (let k = 0; k < L; k++) { sumH += hEdge[ae(k, e)]; transport += hEdge[ae(k, e)] * u[ae(k, e)]; }
      const shift = (avgU[e] - transport) / Math.max(sumH, EPS);
      for (let k = 0; k < L; k++) { const n = ae(k, e); u[n] += shift; if (Math.abs(u[n]) > SPEED_LIMIT) { u[n] = Math.sign(u[n]) * SPEED_LIMIT; limited++; } }
      for (let k = 1; k < L; k++) if (hEdge[ae(k, e)] < THIN) u[ae(k, e)] = u[ae(k - 1, e)];
    }
  }

  function move(i, from, to, amount) {
    const a = at(from, i), b = at(to, i);
    const f = amount / h[a];
    const dQ = Q[a] * f, dW = W[a] * f;
    h[a] -= amount; Q[a] -= dQ; W[a] -= dW;
    h[b] += amount; Q[b] += dQ; W[b] += dW;
  }

  function detrain(i, amount, rm) {
    if (amount <= 0) return;
    let k = 1;
    for (let j = 2; j < L; j++) if (Math.abs(rho[j] - rm) < Math.abs(rho[k] - rm)) k = j;
    move(i, 0, k, amount);
  }

  function mixedLayer(dt) {
    cellVector(mesh, stress, tauCell);
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) continue;
      const tau = Math.hypot(tauCell[3 * i], tauCell[3 * i + 1], tauCell[3 * i + 2]);
      const ustar3 = Math.pow(tau / rho0, 1.5), stir = stirring * Math.exp(-h[i] / stirringDepth);
      let rm = eos(Q[i] / h[i], W[i] / h[i]);
      for (let k = 1; k < L && h[i] < maximumMixedDepth; k++) {
        if (h[at(k, i)] <= EPS || rho[k] > rm) continue;
        move(i, k, 0, Math.min(h[at(k, i)] - EPS, maximumMixedDepth - h[i]));
        rm = eos(Q[i] / h[i], W[i] / h[i]);
      }
      let below = -1;
      for (let k = 1; k < L; k++) if (h[at(k, i)] > THIN) { below = k; break; }
      if (below > 0 && h[i] < maximumMixedDepth) {
        const db = Math.max(1e-4, g * (rho[below] - rm) / rho0);
        const entrain = Math.min(2 * stir * ustar3 / (h[i] * db) * dt, h[at(below, i)] - EPS);
        if (entrain > 0) { move(i, below, 0, entrain); rm = eos(Q[i] / h[i], W[i] / h[i]); }
      }
      below = -1;
      for (let k = 1; k < L; k++) if (h[at(k, i)] > THIN) { below = k; break; }
      let excess = Math.max(0, h[i] - maximumMixedDepth);
      if (below > 0 && rm >= rho[below] - DENSITY_TOLERANCE) excess = Math.max(excess, h[i] - shallowestMixedDepth);
      detrain(i, excess, rm);
      const buoyancy = g * thermalExpansion * (previousT0[i] - surfaceIn[i]) * h[i] / dt;
      if (buoyancy < -1e-9) {
        const monin = Math.max(shallowestMixedDepth, 2 * stir * ustar3 / -buoyancy);
        if (h[i] > monin) detrain(i, (h[i] - monin) * Math.min(1, dt / detrainmentTime), rm);
      }
      if (h[i] < minimumThickness) {
        for (let k = 1; k < L && h[i] < minimumThickness; k++) {
          const available = h[at(k, i)] - EPS;
          if (available > 0) move(i, k, 0, Math.min(available, minimumThickness - h[i]));
        }
      }
    }
  }

  function setStress(total, ice) {
    for (let e = 0; e < E; e++) stress[e] = !edgeOcean[e] || ice[cellsOnEdge[2 * e]] > 0 || ice[cellsOnEdge[2 * e + 1]] > 0 ? 0 : total[e];
  }
  function readSurface(surfaceT, ice) {
    for (let i = 0; i < C; i++) {
      iced[i] = ice[i] > 0 ? 1 : 0;
      T0[i] = iced[i] ? FREEZING_POINT : surfaceT[i];
      surfaceIn[i] = T0[i];
      Q[i] = h[i] * T0[i];
    }
  }
  function accumulate(evaporation, rain, dt, runoffTotal = 0) {
    let area = 0;
    for (let i = 0; i < C; i++) if (cellOcean[i]) area += areaCell[i];
    const spread = area > 0 ? runoffTotal / area : 0;
    for (let i = 0; i < C; i++) if (cellOcean[i]) fresh[i] += (evaporation ? evaporation[i] * dt : 0) - (rain ? rain[i] : 0) - spread;
  }
  function salt(dt, ice) {
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) continue;
      const s = W[i] / h[i];
      W[i] += s * fresh[i] / 1000;
      const grown = (ice[i] - previousIce[i]) * iceDensity / 1000;
      W[i] += (s - iceSalinity) * grown;
      W[i] = Math.max(0, W[i]);
      fresh[i] = 0;
      previousIce[i] = ice[i];
    }
  }
  function writeSurface(surfaceT, oceanFlux, dt) {
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) { oceanFlux[i] = 0; continue; }
      capacity[i] = rhoCp * Math.max(h[i], 1);
      S0[i] = W[i] / h[i];
      if (iced[i]) {
        oceanFlux[i] = rhoCp * (Q[i] - h[i] * FREEZING_POINT) / dt;
        Q[i] = h[i] * FREEZING_POINT;
        T0[i] = FREEZING_POINT;
      } else {
        T0[i] = Q[i] / h[i];
        surfaceT[i] = T0[i];
        oceanFlux[i] = 0;
      }
      previousT0[i] = T0[i];
    }
  }

  function advance(surfaceT, ice, oceanFlux, totalStress, dt) {
    if (++counter % everySteps !== 0) return false;
    if (!initialised) initialize(surfaceT, ice);
    const dtOcean = everySteps * dt;
    relaxRate = Math.min(1 / 3600, 1 / dtOcean);
    readSurface(surfaceT, ice);
    setStress(typeof totalStress === 'function' ? totalStress() : totalStress, ice);
    step(dtOcean);
    mixedLayer(dtOcean);
    salt(dtOcean, ice);
    writeSurface(surfaceT, oceanFlux, dtOcean);
    return true;
  }

  /*
   * The climatological start: a mixed layer over interior layers whose
   * bases sit at their subtropical depths, shallower toward the poles
   * and rising to the surface over one layer's density step as the
   * surface water approaches each layer's density, with the deepest
   * layer filling to the bottom.
   * The free surface starts at its steric height, so the pressure below
   * the thermocline is level and the ocean does not begin with a
   * barotropic shock.
   */
  function initialize(surfaceT, ice) {
    u.fill(0); eta.fill(0); fresh.fill(0);
    for (let i = 0; i < C; i++) {
      previousIce[i] = ice[i];
      iced[i] = ice[i] > 0 ? 1 : 0;
      const lat = latCell[i];
      const s0 = salinityProfile(lat);
      T0[i] = iced[i] ? FREEZING_POINT : surfaceT[i];
      S0[i] = s0;
      if (!cellOcean[i]) { for (let k = 0; k < L; k++) { h[at(k, i)] = 0; Q[at(k, i)] = 0; W[at(k, i)] = 0; } continue; }
      const rm = eos(T0[i], s0);
      const stretch = 1 - thermoclineTilt + 2 * thermoclineTilt * Math.cos(lat) ** 2;
      let cumulative = Math.min(mixedDepth, D[i]);
      h[i] = cumulative; Q[i] = cumulative * T0[i]; W[i] = cumulative * s0;
      for (let k = 1; k < L; k++) {
        let hk;
        if (k === L - 1) hk = Math.max(0, D[i] - cumulative);
        else {
          const taper = Math.max(0, Math.min(1, (rho[k] - rm) / (rho[k + 1] - rho[k])));
          hk = Math.max(0, Math.min(D[i], bottoms[k - 1] * stretch * taper) - cumulative);
        }
        hk = Math.max(EPS, hk);
        cumulative += hk;
        h[at(k, i)] = hk; Q[at(k, i)] = hk * labelT[k]; W[at(k, i)] = hk * referenceS;
      }
      const scale = D[i] / cumulative;
      for (let k = 0; k < L; k++) { h[at(k, i)] *= scale; Q[at(k, i)] *= scale; W[at(k, i)] *= scale; }
      previousT0[i] = T0[i];
      capacity[i] = rhoCp * Math.max(h[i], 1);
    }
    surfaceDensity(h, Q, W);
    stericSurface();
    counter = 0;
    initialised = true;
  }

  /*
   * The free surface that levels the pressure at `referenceDepth` in every
   * column deep enough to reach it, so the deep ocean starts without a
   * barotropic pressure gradient; shallower columns take the free surface
   * of the deep water around them, found by relaxation, as a shelf's sea
   * level follows the ocean beside it. The ocean-mean height is zero, and
   * the height goes into the deepest layer that holds water.
   */
  function stericSurface(referenceDepth = 3500) {
    const rhoRef = rho[L - 1];
    const deep = new Uint8Array(C);
    for (let i = 0; i < C; i++) {
      eta[i] = 0;
      if (!cellOcean[i] || D[i] < referenceDepth) continue;
      deep[i] = 1;
      let budget = referenceDepth, anomaly = 0;
      for (let k = 0; k < L && budget > 0; k++) {
        const part = Math.min(h[at(k, i)], budget);
        anomaly += ((k === 0 ? rhoMl[i] : rho[k]) - rhoRef) * part;
        budget -= part;
      }
      eta[i] = -anomaly / rho0;
    }
    const next = new Float64Array(C);
    for (let sweep = 0; sweep < 2000; sweep++) {
      let moved = 0;
      for (let i = 0; i < C; i++) {
        if (!cellOcean[i] || deep[i]) { next[i] = eta[i]; continue; }
        let sum = 0, weight = 0;
        for (let m = 0; m < nEdgesOnCell[i]; m++) {
          const j = cellsOnCell[maxEdges * i + m];
          if (!cellOcean[j]) continue;
          sum += eta[j]; weight++;
        }
        next[i] = weight ? sum / weight : eta[i];
        moved = Math.max(moved, Math.abs(next[i] - eta[i]));
      }
      eta.set(next);
      if (moved < 1e-7) break;
    }
    let area = 0, mean = 0;
    for (let i = 0; i < C; i++) if (cellOcean[i]) { area += areaCell[i]; mean += areaCell[i] * eta[i]; }
    mean = area > 0 ? mean / area : 0;
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) { eta[i] = 0; continue; }
      eta[i] -= mean;
      let deepest = 0;
      for (let k = L - 1; k >= 1; k--) if (h[at(k, i)] > THIN - eta[i]) { deepest = k; break; }
      const n = at(deepest, i), t = Q[n] / h[n], sal = W[n] / h[n];
      h[n] += eta[i]; Q[n] = h[n] * t; W[n] = h[n] * sal;
    }
  }

  function load(saved, surfaceT, ice) {
    if (!saved.h || saved.h.length !== L * C) {
      initialize(surfaceT, ice);
      if (saved.h1 && saved.u1) {
        for (let i = 0; i < C; i++) {
          if (!cellOcean[i]) continue;
          const wanted = Math.max(minimumThickness, Math.min(saved.h1[i], D[i] - EPS * (L - 1)));
          let below = -1;
          for (let k = 1; k < L; k++) if (h[at(k, i)] > THIN) { below = k; break; }
          if (below > 0) { const delta = Math.max(-(h[i] - minimumThickness), Math.min(wanted - h[i], h[at(below, i)] - EPS)); if (delta > 0) move(i, below, 0, delta); else if (delta < 0) move(i, 0, below, -delta); }
        }
        for (let e = 0; e < E; e++) u[e] = edgeOcean[e] ? saved.u1[e] : 0;
      }
      return;
    }
    h.set(saved.h); u.set(saved.u); eta.set(saved.eta);
    for (let n = 0; n < L * C; n++) { Q[n] = h[n] * saved.T[n]; W[n] = h[n] * saved.S[n]; }
    for (let k = 0; k < L; k++) for (let e = 0; e < E; e++) if (!edgeOcean[e]) u[ae(k, e)] = 0;
    for (let i = 0; i < C; i++) { previousIce[i] = ice[i]; iced[i] = ice[i] > 0 ? 1 : 0; T0[i] = Q[i] / Math.max(EPS, h[i]); S0[i] = W[i] / Math.max(EPS, h[i]); previousT0[i] = T0[i]; capacity[i] = rhoCp * Math.max(h[i], 1); }
    fresh.fill(0);
    counter = 0;
    initialised = true;
  }

  function serialize() {
    const Tall = new Array(L * C), Sall = new Array(L * C);
    for (let n = 0; n < L * C; n++) { const hh = Math.max(EPS, h[n]); Tall[n] = Q[n] / hh; Sall[n] = W[n] / hh; }
    return { h: Array.from(h), u: Array.from(u), T: Tall, S: Sall, eta: Array.from(eta) };
  }

  const thermoclineDepth = new Float64Array(C), sst = T0, sss = S0;
  function fields() {
    for (let i = 0; i < C; i++) thermoclineDepth[i] = cellOcean[i] ? h[i] + h[at(1, i)] + h[at(2, i)] : NaN;
    return { h1: h.subarray(0, C), T1: T0, S1: S0, u1: u.subarray(0, E), T2: thermoclineDepth, eta, thermoclineDepth };
  }

  function diagnostics() {
    let area = 0, depth = 0, heat = 0, thermo = 0, salinity = 0, ssh = 0, speed = 0, transport = 0, interiorT = 0, interiorH = 0;
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) continue;
      const a = areaCell[i];
      area += a; depth += a * h[i]; salinity += a * (W[i] / h[i]); ssh = Math.max(ssh, Math.abs(eta[i]));
      thermo += a * (h[i] + h[at(1, i)] + h[at(2, i)]);
      for (let k = 0; k < L; k++) heat += a * rhoCp * Q[at(k, i)];
      interiorT += Q[at(1, i)] * a; interiorH += h[at(1, i)] * a;
    }
    for (let e = 0; e < E; e++) {
      speed = Math.max(speed, Math.abs(u[e]));
      let t = 0;
      for (let k = 0; k < L; k++) t += hEdge[ae(k, e)] * u[ae(k, e)];
      transport = Math.max(transport, Math.abs(t) * dvEdge[e]);
    }
    return { oceanUpperDepth: depth / area, oceanHeat: heat / area, oceanThermoclineT: interiorH > 0 ? interiorT / interiorH : 0, oceanSpeed: speed, oceanThermoclineDepth: thermo / area, oceanSalinity: salinity / area, oceanSSH: ssh, oceanTransport: transport / 1e6, oceanLimited: limited };
  }

  initialize(new Float64Array(C).fill(288), new Float64Array(C));
  initialised = false;
  return { state, layers: L, h, u, Q, W, eta, D, T0, S0, capacity, stress, fresh, tendency, advance, accumulate, initialize, load, serialize, diagnostics, fields, setStress, readSurface, rhoCp, edgeOcean, cellOcean, densities: rho, shared: { capacity: capacity.buffer } };
}
