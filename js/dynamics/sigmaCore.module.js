import { divergence, gradient, curl, kineticEnergy, laplacianVelocity, laplacianScalar } from './operators.module.js';

export const R_DRY = 287.06;
export const VIRTUAL_FACTOR = 0.608;
export const CP_DRY = 1003.5;
export const P0 = 101325;
export const GRAVITY = 9.806;

const CAM_L26_HYAI = [
  0.00219406700000001, 0.00489520900000001, 0.009882418, 0.01805201,
  0.02983724, 0.0446233400000002, 0.0616058700000002, 0.0785124300000004,
  0.0773127100000002, 0.0759013100000003, 0.0742408600000002,
  0.0722874400000002, 0.0699893299999998, 0.06728574, 0.06410509,
  0.0603632200000002, 0.0559611100000001, 0.0507822500000001,
  0.0446896000000001, 0.0375219099999999, 0.0290894900000001, 0.02084739,
  0.01334443, 0.00708499000000001, 0.00252136, 0, 0,
];
const CAM_L26_HYBI = [
  0, 0, 0, 0, 0, 0, 0, 0, 0.01505309, 0.03276228, 0.05359622,
  0.0781062700000006, 0.1069411, 0.140863700000001, 0.180772, 0.227722,
  0.282956200000001, 0.347936400000002, 0.4243822, 0.514316800000003,
  0.620120200000002, 0.723535500000004, 0.817676800000001,
  0.896215300000001, 0.953476100000003, 0.9851122, 1,
];

/*
 * Interface sigma values, top (0) to ground (1): the CAM 26-level hybrid
 * grid that Jablonowski & Williamson 2006 ran on, read as sigma with
 * p_s = p0 (hyai + hybi from HOMME's cami-26.ascii), plus the cap above
 * CAM's 2.19 hPa lid that a sigma coordinate needs. 27 layers.
 */
export function sigmaInterfaces() {
  const levels = [0];
  for (let k = 0; k < CAM_L26_HYAI.length; k++) levels.push(CAM_L26_HYAI[k] + CAM_L26_HYBI[k]);
  levels[levels.length - 1] = 1;
  return Float64Array.from(levels);
}

/*
 * Hydrostatic primitive equations in sigma coordinates on the C-grid.
 * State: pi[C] surface pressure, theta[K*C] layer potential temperature,
 * u[K*E] layer normal velocity (layer k occupies [k*C, (k+1)*C) and
 * [k*E, (k+1)*E)), and optionally state[4] = q[K*C], specific humidity,
 * and state[5] = qc[K*C], cloud condensate, both carried with the same
 * flux-form transport as theta. Interface
 * quantities use K+1 slots per column, top first. Steps follow the
 * A-grid column: mass flux and divergence per layer, dpi/dt,
 * pi*sigma-dot telescoping to zero at the ground, Exner ratios from the
 * exact layer integral of sigma^kappa, interface theta interpolated in
 * Exner, geopotential integrated upward with each layer's virtual theta
 * over its own Exner span, then flux-form theta and q transport and the
 * vector-invariant momentum equation with the RTSK PV flux. The
 * pressure-gradient force uses virtual theta (with condensate loading);
 * with no q it equals theta.
 */
export function createSigmaCore(mesh, options = {}) {
  const {
    levels = sigmaInterfaces(), g = GRAVITY, cp = CP_DRY, R = R_DRY, p0 = P0,
    nu4 = 0, nu4Theta = 0, forcing = null, surfaceGeopotential = null, buffers = null,
  } = options;
  const {
    nCells: C, nEdges: E, nVertices: V, maxEdgesOnEdge, nEdgesOnEdge, edgesOnEdge, weightsOnEdge,
    cellsOnEdge, verticesOnEdge, cellsOnVertex, kiteAreasOnVertex, areaTriangle, dcEdge, dvEdge, fVertex,
  } = mesh;
  const K = levels.length - 1;
  const kappa = R / cp;
  let applyForcing = forcing;
  const sigmaUpper = levels.subarray(0, K);
  const sigmaLower = levels.subarray(1, K + 1);
  const dSigma = Float64Array.from(sigmaLower, (s, k) => s - sigmaUpper[k]);
  const sigmaMid = Float64Array.from(sigmaLower, (s, k) => 0.5 * (s + sigmaUpper[k]));

  const shared = {};
  const sharedArray = (name, n) => {
    const buffer = buffers && buffers[name] ? buffers[name] : new SharedArrayBuffer(8 * n);
    shared[name] = buffer;
    return new Float64Array(buffer);
  };
  const flux = sharedArray('flux', K * E);
  const divFlux = sharedArray('divFlux', K * C);
  const piSigmaDot = sharedArray('piSigmaDot', (K + 1) * C);
  const exnerLower = sharedArray('exnerLower', K * C);
  const exnerLayer = sharedArray('exnerLayer', K * C);
  const dExnerDpi = sharedArray('dExnerDpi', K * C);
  const thetaLower = sharedArray('thetaLower', K * C);
  const qLower = sharedArray('qLower', K * C);
  const qcLower = sharedArray('qcLower', K * C);
  const thetaV = sharedArray('thetaV', K * C);
  const geopotential = sharedArray('geopotential', K * C);
  const piVertex = sharedArray('piVertex', V);
  const piEdge = new Float64Array(E);
  const thetaFlux = new Float64Array(E);
  const divThetaFlux = new Float64Array(C);
  const zeta = new Float64Array(V);
  const qVertex = new Float64Array(V);
  const qEdge = new Float64Array(E);
  const kinetic = new Float64Array(C);
  const phi = new Float64Array(C);
  const gradPhi = new Float64Array(E);
  const gradPi = new Float64Array(E);
  const lap = new Float64Array(E);
  const lap2 = new Float64Array(E);
  const divScratch = new Float64Array(C);
  const curlScratch = new Float64Array(V);
  const lapTheta = new Float64Array(C);
  const lapTheta2 = new Float64Array(C);
  const qFlux = new Float64Array(E);
  const divQFlux = new Float64Array(C);
  const qcFlux = new Float64Array(E);
  const divQcFlux = new Float64Array(C);
  const massField = new Float64Array(C);

  function diagnoseColumn(i, pi, theta, q = null, qc = null) {
    for (let k = 0; k < K; k++) {
      const idx = k * C + i;
      const exLower = Math.pow(pi[i] * sigmaLower[k] / p0, kappa);
      const exUpper = k === 0 ? 0 : exnerLower[idx - C];
      const span = exLower * sigmaLower[k] - exUpper * sigmaUpper[k];
      exnerLower[idx] = exLower;
      exnerLayer[idx] = span / ((1 + kappa) * dSigma[k]);
      dExnerDpi[idx] = (kappa / (1 + kappa)) * span / (pi[i] * dSigma[k]);
    }
    for (let k = 0; k < K; k++) {
      const idx = k * C + i;
      thetaV[idx] = q ? theta[idx] * (1 + VIRTUAL_FACTOR * q[idx] - (qc ? qc[idx] : 0)) : theta[idx];
    }
    for (let k = 0; k < K - 1; k++) {
      const idx = k * C + i;
      const t = (exnerLower[idx] - exnerLayer[idx]) / (exnerLayer[idx + C] - exnerLayer[idx]);
      thetaLower[idx] = theta[idx] + t * (theta[idx + C] - theta[idx]);
      qLower[idx] = q ? q[idx] + t * (q[idx + C] - q[idx]) : 0;
      qcLower[idx] = qc ? qc[idx] + t * (qc[idx + C] - qc[idx]) : 0;
    }
    const bottom = (K - 1) * C + i;
    geopotential[bottom] = (surfaceGeopotential ? surfaceGeopotential[i] : 0) + cp * thetaV[bottom] * (exnerLower[bottom] - exnerLayer[bottom]);
    for (let k = K - 2; k >= 0; k--) {
      const idx = k * C + i;
      const below = idx + C;
      geopotential[idx] = geopotential[below]
        + cp * thetaV[below] * (exnerLayer[below] - exnerLower[idx])
        + cp * thetaV[idx] * (exnerLower[idx] - exnerLayer[idx]);
    }
  }

  function diagnose(pi, theta, q = null, qc = null) {
    for (let i = 0; i < C; i++) diagnoseColumn(i, pi, theta, q, qc);
  }

  function edgePi(pi) {
    for (let e = 0; e < E; e++) piEdge[e] = 0.5 * (pi[cellsOnEdge[2 * e]] + pi[cellsOnEdge[2 * e + 1]]);
  }

  function phaseFlux(state, kFrom, kTo) {
    const [pi, , u] = state;
    edgePi(pi);
    for (let k = kFrom; k < kTo; k++) {
      const fk = flux.subarray(k * E, (k + 1) * E);
      const uk = u.subarray(k * E, (k + 1) * E);
      for (let e = 0; e < E; e++) fk[e] = piEdge[e] * uk[e];
      divergence(mesh, fk, divFlux.subarray(k * C, (k + 1) * C));
    }
  }

  function phaseColumn(state, out, iFrom, iTo) {
    const [pi, theta] = state;
    const q = state[4] ?? null, qc = state[5] ?? null;
    const [dPi] = out;
    for (let i = iFrom; i < iTo; i++) {
      let sum = 0;
      for (let k = 0; k < K; k++) sum += divFlux[k * C + i] * dSigma[k];
      dPi[i] = -sum;
      let cumulative = 0;
      piSigmaDot[i] = 0;
      for (let k = 0; k < K; k++) {
        cumulative += divFlux[k * C + i] * dSigma[k];
        piSigmaDot[(k + 1) * C + i] = -cumulative - sigmaLower[k] * dPi[i];
      }
      piSigmaDot[K * C + i] = 0;
      diagnoseColumn(i, pi, theta, q, qc);
    }
  }

  function phaseVertex(state, vFrom, vTo) {
    const [pi] = state;
    for (let v = vFrom; v < vTo; v++) {
      let sum = 0;
      for (let m = 0; m < 3; m++) sum += kiteAreasOnVertex[3 * v + m] * pi[cellsOnVertex[3 * v + m]];
      piVertex[v] = sum / areaTriangle[v];
    }
  }

  /*
   * Flux-form transport of a layer scalar by the mass fluxes, with the
   * ∇⁴ closure. `conservative` applies the closure to the mass-weighted
   * field so the column's total is preserved exactly under it (used for
   * water); theta keeps the closure on the field itself.
   */
  function transportLayer(k, pi, field, fieldLower, edgeFlux, divField, dField, conservative = false) {
    const off = k * C;
    const fk = flux.subarray(k * E, (k + 1) * E);
    for (let e = 0; e < E; e++) {
      edgeFlux[e] = fk[e] * 0.5 * (field[off + cellsOnEdge[2 * e]] + field[off + cellsOnEdge[2 * e + 1]]);
    }
    divergence(mesh, edgeFlux, divField);
    for (let i = 0; i < C; i++) {
      const idx = off + i;
      const lowerFlow = piSigmaDot[(k + 1) * C + i];
      const upperFlow = piSigmaDot[k * C + i];
      const lower = k === K - 1 ? 0 : fieldLower[idx];
      const upper = k === 0 ? 0 : fieldLower[idx - C];
      const vertical = (lowerFlow * lower - upperFlow * upper - field[idx] * (lowerFlow - upperFlow)) / (pi[i] * dSigma[k]);
      dField[idx] = -(divField[i] - field[idx] * divFlux[idx]) / pi[i] - vertical;
    }
    if (nu4Theta > 0) {
      if (conservative) {
        for (let i = 0; i < C; i++) massField[i] = pi[i] * field[off + i];
        laplacianScalar(mesh, massField, lapTheta);
        laplacianScalar(mesh, lapTheta, lapTheta2);
        for (let i = 0; i < C; i++) dField[off + i] -= nu4Theta * lapTheta2[i] / pi[i];
      } else {
        laplacianScalar(mesh, field.subarray(off, off + C), lapTheta);
        laplacianScalar(mesh, lapTheta, lapTheta2);
        for (let i = 0; i < C; i++) dField[off + i] -= nu4Theta * lapTheta2[i];
      }
    }
  }

  function phaseLayer(state, out, kFrom, kTo) {
    const [pi, theta, u] = state;
    const [, dTheta, dU] = out;
    const q = state[4] ?? null, dQ = out[4] ?? null, qc = state[5] ?? null, dQc = out[5] ?? null;
    edgePi(pi);
    gradient(mesh, pi, gradPi);
    for (let k = kFrom; k < kTo; k++) {
      const off = k * C;
      const fk = flux.subarray(k * E, (k + 1) * E);
      transportLayer(k, pi, theta, thetaLower, thetaFlux, divThetaFlux, dTheta);
      if (q && dQ) transportLayer(k, pi, q, qLower, qFlux, divQFlux, dQ, true);
      if (qc && dQc) transportLayer(k, pi, qc, qcLower, qcFlux, divQcFlux, dQc, true);
      const uk = u.subarray(k * E, (k + 1) * E);
      const dUk = dU.subarray(k * E, (k + 1) * E);
      curl(mesh, uk, zeta);
      for (let v = 0; v < V; v++) qVertex[v] = (zeta[v] + fVertex[v]) / piVertex[v];
      for (let e = 0; e < E; e++) qEdge[e] = 0.5 * (qVertex[verticesOnEdge[2 * e]] + qVertex[verticesOnEdge[2 * e + 1]]);
      kineticEnergy(mesh, uk, kinetic);
      for (let i = 0; i < C; i++) phi[i] = geopotential[off + i] + kinetic[i];
      gradient(mesh, phi, gradPhi);
      for (let e = 0; e < E; e++) {
        const i = cellsOnEdge[2 * e], j = cellsOnEdge[2 * e + 1];
        let pv = 0;
        for (let s = 0; s < nEdgesOnEdge[e]; s++) {
          const other = edgesOnEdge[maxEdgesOnEdge * e + s];
          pv += weightsOnEdge[maxEdgesOnEdge * e + s] * dvEdge[other] * fk[other] * 0.5 * (qEdge[e] + qEdge[other]);
        }
        const pgfPi = cp * 0.5 * (thetaV[off + i] * dExnerDpi[off + i] + thetaV[off + j] * dExnerDpi[off + j]) * gradPi[e];
        const lowerFlow = 0.5 * (piSigmaDot[(k + 1) * C + i] + piSigmaDot[(k + 1) * C + j]);
        const upperFlow = 0.5 * (piSigmaDot[k * C + i] + piSigmaDot[k * C + j]);
        const lowerU = k === K - 1 ? 0 : 0.5 * (uk[e] + u[(k + 1) * E + e]);
        const upperU = k === 0 ? 0 : 0.5 * (uk[e] + u[(k - 1) * E + e]);
        const vertical = (lowerFlow * lowerU - upperFlow * upperU - uk[e] * (lowerFlow - upperFlow)) / (piEdge[e] * dSigma[k]);
        dUk[e] = pv / dcEdge[e] - gradPhi[e] - pgfPi - vertical;
      }
      if (nu4 > 0) {
        laplacianVelocity(mesh, uk, lap, divScratch, curlScratch);
        laplacianVelocity(mesh, lap, lap2, divScratch, curlScratch);
        for (let e = 0; e < E; e++) dUk[e] -= nu4 * lap2[e];
      }
    }
  }

  function tendency(state, out) {
    phaseFlux(state, 0, K);
    phaseColumn(state, out, 0, C);
    phaseVertex(state, 0, V);
    phaseLayer(state, out, 0, K);
    if (applyForcing) applyForcing(state, out, diagnostics);
  }

  function setForcing(fn) {
    applyForcing = fn;
  }

  const diagnostics = { K, C, E, V, levels, sigmaMid, sigmaLower, sigmaUpper, dSigma, kappa, cp, R, g, p0, exnerLayer, exnerLower, dExnerDpi, geopotential, piSigmaDot, diagnose, diagnoseColumn };

  function mass(pi) {
    let m = 0;
    for (let i = 0; i < C; i++) m += mesh.areaCell[i] * pi[i];
    return m / g;
  }

  return { K, levels, sigmaMid, tendency, phaseFlux, phaseColumn, phaseVertex, phaseLayer, diagnose, diagnoseColumn, diagnostics, mass, setForcing, shared, arrays: { exnerLayer, exnerLower, dExnerDpi, geopotential, piSigmaDot, thetaLower, qLower, qcLower, thetaV } };
}

/*
 * Held & Suarez 1994 forcing: Newtonian relaxation of temperature toward
 * the prescribed radiative-equilibrium profile and Rayleigh friction in
 * the boundary layer, as tendencies added to theta and u.
 */
export function createHeldSuarez(mesh, core, {
  kf = 1 / 86400, ka = 1 / (40 * 86400), ks = 1 / (4 * 86400), sigmaB = 0.7,
  deltaTy = 60, deltaThetaZ = 10, tMin = 200, tMax = 315, pRef = 1e5,
} = {}) {
  const { K, C, E, sigmaMid, kappa, cp, exnerLayer } = core.diagnostics;
  const cosLat = Float64Array.from(mesh.latCell, Math.cos);
  const sinLat = Float64Array.from(mesh.latCell, Math.sin);
  const cosLatEdge = Float64Array.from(mesh.latEdge, Math.cos);
  return function forcing(state, out) {
    const [pi, theta, u] = state;
    const [, dTheta, dU] = out;
    for (let k = 0; k < K; k++) {
      const s = sigmaMid[k];
      const weight = Math.max(0, (s - sigmaB) / (1 - sigmaB));
      for (let i = 0; i < C; i++) {
        const idx = k * C + i;
        const p = pi[i] * s;
        const c2 = cosLat[i] * cosLat[i];
        const tEq = Math.max(tMin, (tMax - deltaTy * sinLat[i] * sinLat[i] - deltaThetaZ * Math.log(p / pRef) * c2) * Math.pow(p / pRef, kappa));
        const kT = ka + (ks - ka) * weight * c2 * c2;
        dTheta[idx] -= kT * (theta[idx] - tEq / exnerLayer[idx]);
      }
      if (weight > 0) {
        for (let e = 0; e < E; e++) dU[k * E + e] -= kf * weight * u[k * E + e];
      }
    }
  };
}
