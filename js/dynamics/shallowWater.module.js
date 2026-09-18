import { divergence, gradient, curl, kineticEnergy, laplacianVelocity } from './operators.module.js';

/*
 * Nonlinear shallow water on the C-grid in the vector-invariant form of
 * Ringler et al. 2010:
 *   dh/dt = -div(h_e u)
 *   du_e/dt = Q⊥_e - grad(g(h + b) + K)_e - ν4 ∇⁴u_e
 * with q_v = (ζ_v + f_v) / h_v on the dual mesh, h_v the kite-weighted
 * cell thickness, and the PV flux Q⊥ reconstructed with the PV averaged
 * inside the TRiSK sum (their eq. 49), which is what makes it
 * energy-conserving.
 */
export function createShallowWater(mesh, { g = 9.80616, bottom = null, nu4 = 0 } = {}) {
  const {
    nCells, nEdges, nVertices, maxEdgesOnEdge, nEdgesOnEdge, edgesOnEdge, weightsOnEdge,
    cellsOnEdge, verticesOnEdge, cellsOnVertex, kiteAreasOnVertex, areaTriangle, areaCell,
    dcEdge, dvEdge, fVertex,
  } = mesh;

  const hEdge = new Float64Array(nEdges);
  const flux = new Float64Array(nEdges);
  const zeta = new Float64Array(nVertices);
  const qVertex = new Float64Array(nVertices);
  const qEdge = new Float64Array(nEdges);
  const K = new Float64Array(nCells);
  const phi = new Float64Array(nCells);
  const gradPhi = new Float64Array(nEdges);
  const lap = new Float64Array(nEdges);
  const lap2 = new Float64Array(nEdges);
  const divScratch = new Float64Array(nCells);
  const curlScratch = new Float64Array(nVertices);

  function thicknessOnVertices(h, out) {
    for (let v = 0; v < nVertices; v++) {
      let sum = 0;
      for (let m = 0; m < 3; m++) sum += kiteAreasOnVertex[3 * v + m] * h[cellsOnVertex[3 * v + m]];
      out[v] = sum / areaTriangle[v];
    }
    return out;
  }

  function tendency(h, u, dh, du) {
    for (let e = 0; e < nEdges; e++) {
      hEdge[e] = 0.5 * (h[cellsOnEdge[2 * e]] + h[cellsOnEdge[2 * e + 1]]);
      flux[e] = hEdge[e] * u[e];
    }
    divergence(mesh, flux, dh);
    for (let i = 0; i < nCells; i++) dh[i] = -dh[i];

    curl(mesh, u, zeta);
    thicknessOnVertices(h, qVertex);
    for (let v = 0; v < nVertices; v++) qVertex[v] = (zeta[v] + fVertex[v]) / qVertex[v];
    for (let e = 0; e < nEdges; e++) {
      qEdge[e] = 0.5 * (qVertex[verticesOnEdge[2 * e]] + qVertex[verticesOnEdge[2 * e + 1]]);
    }

    kineticEnergy(mesh, u, K);
    for (let i = 0; i < nCells; i++) phi[i] = g * (h[i] + (bottom ? bottom[i] : 0)) + K[i];
    gradient(mesh, phi, gradPhi);

    for (let e = 0; e < nEdges; e++) {
      let sum = 0;
      for (let s = 0; s < nEdgesOnEdge[e]; s++) {
        const other = edgesOnEdge[maxEdgesOnEdge * e + s];
        sum += weightsOnEdge[maxEdgesOnEdge * e + s] * dvEdge[other] * flux[other] * 0.5 * (qEdge[e] + qEdge[other]);
      }
      du[e] = sum / dcEdge[e] - gradPhi[e];
    }

    if (nu4 > 0) {
      laplacianVelocity(mesh, u, lap, divScratch, curlScratch);
      laplacianVelocity(mesh, lap, lap2, divScratch, curlScratch);
      for (let e = 0; e < nEdges; e++) du[e] -= nu4 * lap2[e];
    }
  }

  function diagnostics(h, u) {
    let mass = 0, potential = 0, kinetic = 0, enstrophy = 0;
    for (let i = 0; i < nCells; i++) {
      const b = bottom ? bottom[i] : 0;
      mass += areaCell[i] * h[i];
      potential += areaCell[i] * g * (0.5 * h[i] * h[i] + h[i] * b);
    }
    for (let e = 0; e < nEdges; e++) {
      const he = 0.5 * (h[cellsOnEdge[2 * e]] + h[cellsOnEdge[2 * e + 1]]);
      kinetic += 0.5 * dcEdge[e] * dvEdge[e] * he * u[e] * u[e];
    }
    curl(mesh, u, zeta);
    thicknessOnVertices(h, curlScratch);
    for (let v = 0; v < nVertices; v++) {
      const q = (zeta[v] + fVertex[v]) / curlScratch[v];
      enstrophy += areaTriangle[v] * 0.5 * curlScratch[v] * q * q;
    }
    return { mass, energy: potential + kinetic, kinetic, potential, enstrophy };
  }

  return { tendency, diagnostics, thicknessOnVertices };
}
