import { divergence, gradient, curl, kineticEnergy, laplacianVelocity, laplacianScalar } from '../dynamics/operators.module.js';
import { createRK4Arrays } from '../dynamics/integrators.module.js';
import { FREEZING_POINT } from '../physics/ice.module.js';

/*
 * A reduced-gravity two-layer ocean on the C-grid: an upper (mixed)
 * layer and a thermocline layer of prescribed density contrasts over a
 * motionless abyss. Each layer is a TRiSK shallow-water layer in the
 * vector-invariant form of shallowWater.module.js, driven by the
 * Montgomery potential of a resting abyss,
 *   M1 = g'12 h1 + g'23 (h1 + h2),   M2 = g'23 (h1 + h2),
 * so the fastest wave is internal (√(g'h), 2–3 m/s) and the ocean takes
 * everySteps atmosphere steps at once. Momentum: the surface stress on
 * the upper layer is the whole momentum the atmosphere loses to the
 * surface, aerodynamic drag and boundary-layer damping together (zero
 * under ice), interfacial and bottom drag as linear
 * stresses ρ r Δu, and a ∇⁴ closure of e-folding time closureHours at
 * the grid scale, long enough for the ocean's longer step. Heat: each layer carries h·T in flux
 * form with centred edge temperatures; the upper layer diffuses with
 * the energy-balance diffusivity of M11; a layer thinner than
 * minimumThickness entrains from below over entrainmentTime, the
 * thermocline layer from the abyss at abyssTemperature, which is the one
 * exchange the ocean's heat budget does not close. The upper layer's
 * temperature is the sea surface temperature: `advance` reads it from
 * surfaceT over open water (the freezing point under ice), steps the
 * ocean, writes it back, publishes the upper layer's heat capacity per
 * cell, and hands the heat converged under ice to oceanFlux for the ice
 * base over the atmosphere steps until the next ocean step.
 */
export function createOcean(mesh, {
  upperDepth = 50, lowerDepth = 350, reducedGravity = 0.02, abyssReducedGravity = 0.01, abyssTemperature = 275,
  minimumThickness = 10, entrainmentTime = 86400, density = 1025, specificHeat = 3985, interfacialDrag = 2e-4, bottomDrag = 2e-4,
  closureHours = 12, diffusivity = 0.45, everySteps = 4, buffers = null,
} = {}) {
  const {
    nCells: C, nEdges: E, nVertices: V, maxEdgesOnEdge, nEdgesOnEdge, edgesOnEdge, weightsOnEdge,
    cellsOnEdge, verticesOnEdge, cellsOnVertex, kiteAreasOnVertex, areaTriangle, areaCell, dcEdge, dvEdge, fVertex, radius,
  } = mesh;
  const rhoCp = density * specificHeat;
  const diffusion = diffusivity * radius * radius / rhoCp;
  let spacing = 0;
  for (let e = 0; e < E; e++) spacing += dcEdge[e];
  spacing /= E;
  const nu4 = closureHours > 0 ? Math.pow(spacing / Math.PI, 4) / (closureHours * 3600) : 0;
  const h1 = new Float64Array(C).fill(upperDepth), h2 = new Float64Array(C).fill(lowerDepth);
  const u1 = new Float64Array(E), u2 = new Float64Array(E);
  const H1 = new Float64Array(C), H2 = new Float64Array(C);
  const state = [h1, h2, u1, u2, H1, H2];
  const T1 = new Float64Array(C), T2 = new Float64Array(C);
  const capacity = new Float64Array(buffers && buffers.capacity ? buffers.capacity : new SharedArrayBuffer(8 * C)).fill(rhoCp * upperDepth);
  const stress = new Float64Array(E);
  const iced = new Uint8Array(C);
  const hEdge = new Float64Array(E), flux = new Float64Array(E), heatFlux = new Float64Array(E), uDiff = new Float64Array(E);
  const zeta = new Float64Array(V), qVertex = new Float64Array(V), qEdge = new Float64Array(E);
  const K = new Float64Array(C), phi = new Float64Array(C), M1 = new Float64Array(C), M2 = new Float64Array(C), T = new Float64Array(C), lapT = new Float64Array(C);
  const gradPhi = new Float64Array(E), lap = new Float64Array(E), lap2 = new Float64Array(E), divScratch = new Float64Array(C), curlScratch = new Float64Array(V);
  const rk4 = createRK4Arrays(state.map((a) => a.length));
  let counter = 0;

  function thicknessOnVertices(h, out) {
    for (let v = 0; v < V; v++) {
      let sum = 0;
      for (let m = 0; m < 3; m++) sum += kiteAreasOnVertex[3 * v + m] * h[cellsOnVertex[3 * v + m]];
      out[v] = sum / areaTriangle[v];
    }
  }

  function layer(h, u, H, M, dh, du, dH, upper) {
    for (let e = 0; e < E; e++) {
      hEdge[e] = 0.5 * (h[cellsOnEdge[2 * e]] + h[cellsOnEdge[2 * e + 1]]);
      flux[e] = hEdge[e] * u[e];
    }
    divergence(mesh, flux, dh);
    for (let i = 0; i < C; i++) { dh[i] = -dh[i]; T[i] = H[i] / h[i]; }
    for (let e = 0; e < E; e++) heatFlux[e] = flux[e] * 0.5 * (T[cellsOnEdge[2 * e]] + T[cellsOnEdge[2 * e + 1]]);
    divergence(mesh, heatFlux, dH);
    for (let i = 0; i < C; i++) dH[i] = -dH[i];
    if (upper && diffusion > 0) {
      laplacianScalar(mesh, T, lapT);
      for (let i = 0; i < C; i++) dH[i] += diffusion * lapT[i];
    }
    curl(mesh, u, zeta);
    thicknessOnVertices(h, qVertex);
    for (let v = 0; v < V; v++) qVertex[v] = (zeta[v] + fVertex[v]) / qVertex[v];
    for (let e = 0; e < E; e++) qEdge[e] = 0.5 * (qVertex[verticesOnEdge[2 * e]] + qVertex[verticesOnEdge[2 * e + 1]]);
    kineticEnergy(mesh, u, K);
    for (let i = 0; i < C; i++) phi[i] = M[i] + K[i];
    gradient(mesh, phi, gradPhi);
    for (let e = 0; e < E; e++) {
      let sum = 0;
      for (let s = 0; s < nEdgesOnEdge[e]; s++) {
        const other = edgesOnEdge[maxEdgesOnEdge * e + s];
        sum += weightsOnEdge[maxEdgesOnEdge * e + s] * dvEdge[other] * flux[other] * 0.5 * (qEdge[e] + qEdge[other]);
      }
      du[e] = sum / dcEdge[e] - gradPhi[e];
    }
    if (upper) for (let e = 0; e < E; e++) du[e] += (stress[e] / density - interfacialDrag * uDiff[e]) / hEdge[e];
    else for (let e = 0; e < E; e++) du[e] += (interfacialDrag * uDiff[e] - bottomDrag * u[e]) / hEdge[e];
    if (nu4 > 0) {
      laplacianVelocity(mesh, u, lap, divScratch, curlScratch);
      laplacianVelocity(mesh, lap, lap2, divScratch, curlScratch);
      for (let e = 0; e < E; e++) du[e] -= nu4 * lap2[e];
    }
  }

  function tendency(input, out) {
    const [a1, a2, v1, v2, Q1, Q2] = input;
    const [dh1, dh2, du1, du2, dH1, dH2] = out;
    for (let e = 0; e < E; e++) uDiff[e] = v1[e] - v2[e];
    for (let i = 0; i < C; i++) {
      M2[i] = abyssReducedGravity * (a1[i] + a2[i]);
      M1[i] = M2[i] + reducedGravity * a1[i];
    }
    layer(a1, v1, Q1, M1, dh1, du1, dH1, true);
    layer(a2, v2, Q2, M2, dh2, du2, dH2, false);
    for (let i = 0; i < C; i++) {
      if (a1[i] < minimumThickness) {
        const w = (minimumThickness - a1[i]) / entrainmentTime, t2 = Q2[i] / a2[i];
        dh1[i] += w; dh2[i] -= w; dH1[i] += w * t2; dH2[i] -= w * t2;
      }
      if (a2[i] < minimumThickness) {
        const w = (minimumThickness - a2[i]) / entrainmentTime;
        dh2[i] += w; dH2[i] += w * abyssTemperature;
      }
    }
  }

  function setStress(total, ice) {
    for (let e = 0; e < E; e++) stress[e] = ice[cellsOnEdge[2 * e]] > 0 || ice[cellsOnEdge[2 * e + 1]] > 0 ? 0 : total[e];
  }

  function readSurface(surfaceT, ice) {
    for (let i = 0; i < C; i++) {
      iced[i] = ice[i] > 0 ? 1 : 0;
      T1[i] = iced[i] ? FREEZING_POINT : surfaceT[i];
      H1[i] = h1[i] * T1[i];
    }
  }

  function advance(surfaceT, ice, oceanFlux, totalStress, dt) {
    if (++counter % everySteps !== 0) return false;
    const dtOcean = everySteps * dt;
    readSurface(surfaceT, ice);
    setStress(typeof totalStress === 'function' ? totalStress() : totalStress, ice);
    rk4(tendency, state, dtOcean);
    for (let i = 0; i < C; i++) {
      capacity[i] = rhoCp * h1[i];
      T2[i] = H2[i] / h2[i];
      if (iced[i]) {
        oceanFlux[i] = rhoCp * (H1[i] - h1[i] * FREEZING_POINT) / dtOcean;
        H1[i] = h1[i] * FREEZING_POINT;
        T1[i] = FREEZING_POINT;
      } else {
        T1[i] = H1[i] / h1[i];
        surfaceT[i] = T1[i];
        oceanFlux[i] = 0;
      }
    }
    return true;
  }

  function initialize(surfaceT, ice) {
    h1.fill(upperDepth); h2.fill(lowerDepth); u1.fill(0); u2.fill(0);
    readSurface(surfaceT, ice);
    for (let i = 0; i < C; i++) {
      T2[i] = Math.max(FREEZING_POINT, abyssTemperature + 0.5 * (T1[i] - abyssTemperature));
      H2[i] = h2[i] * T2[i];
      capacity[i] = rhoCp * h1[i];
    }
    counter = 0;
  }

  function load(saved, surfaceT, ice) {
    h1.set(saved.h1); h2.set(saved.h2); u1.set(saved.u1); u2.set(saved.u2); T2.set(saved.T2);
    readSurface(surfaceT, ice);
    for (let i = 0; i < C; i++) { H2[i] = h2[i] * T2[i]; capacity[i] = rhoCp * h1[i]; }
    counter = 0;
  }

  function serialize() {
    return { h1: Array.from(h1), h2: Array.from(h2), u1: Array.from(u1), u2: Array.from(u2), T2: Array.from(T2) };
  }

  function diagnostics() {
    let area = 0, depth = 0, heat = 0, thermocline = 0, speed = 0;
    for (let i = 0; i < C; i++) {
      const a = areaCell[i];
      area += a; depth += a * h1[i]; heat += a * rhoCp * (H1[i] + H2[i]); thermocline += a * T2[i];
    }
    for (let e = 0; e < E; e++) speed = Math.max(speed, Math.abs(u1[e]));
    return { oceanUpperDepth: depth / area, oceanHeat: heat / area, oceanThermoclineT: thermocline / area, oceanSpeed: speed };
  }

  return { state, h1, h2, u1, u2, H1, H2, T1, T2, capacity, stress, tendency, advance, initialize, load, serialize, diagnostics, setStress, readSurface, rhoCp, shared: { capacity: capacity.buffer } };
}
