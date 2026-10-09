/*
 * A row of hydrostatic σ-coordinate columns, periodic in x: the model's
 * dynamical core in one horizontal dimension. Prognostic surface pressure
 * π, layer mass-weighted potential temperature πθ and the wind on the
 * edges between columns; geopotential and Exner from the hydrostatic
 * integration, σ̇ from continuity, upwind transport, RK3 in time, then a
 * dry convective adjustment. Heating goes into one column.
 */
import { R, CP, KAPPA, G, P0, thetaAt, heightOf } from './physics.module.js';

export function createSlice({ M = 64, K = 12, dx = 1e5, heated = M >> 1, heatedHalf = 1, dt = 60 } = {}) {
  const dSigma = 1 / K, n = M * K;
  const sigmaK = new Float64Array(K + 1), sigma1K = new Float64Array(K + 1);
  for (let k = 0; k <= K; k++) { sigmaK[k] = (k / K) ** KAPPA; sigma1K[k] = (k / K) ** (1 + KAPPA); }
  const theta0 = Float64Array.from({ length: K }, (_, k) => thetaAt(heightOf((k + 0.5) * dSigma * P0)));
  const z0 = new Float64Array(K + 1);
  const sponge = Float64Array.from({ length: M }, (_, i) => { const d = Math.min(i, M - 1 - i); return (1 / 1800) * Math.max(0, 1 - d / 8) ** 2; });

  const pi = new Float64Array(M), Theta = new Float64Array(n), u = new Float64Array(n);
  const exnerLayer = new Float64Array(n), phiLayer = new Float64Array(n), phiInterface = new Float64Array(M * (K + 1));
  const piSigmaDot = new Float64Array(M * (K + 1)), theta = new Float64Array(n), F = new Float64Array(n), D = new Float64Array(n);
  const stages = [0, 1, 2].map(() => ({ pi: new Float64Array(M), Theta: new Float64Array(n), u: new Float64Array(n) }));
  const tend = [0, 1, 2].map(() => ({ pi: new Float64Array(M), Theta: new Float64Array(n), u: new Float64Array(n) }));
  const state = { time: 0 };

  function reset() {
    pi.fill(P0); u.fill(0);
    for (let k = 0; k < K; k++) for (let i = 0; i < M; i++) Theta[k * M + i] = P0 * theta0[k];
    state.time = 0;
    tendency({ pi, Theta, u }, tend[0], 0, 'column', false);
    for (let k = 0; k <= K; k++) z0[k] = phiInterface[k * M] / G;
  }

  function hydrostatic(pi, Theta) {
    for (let i = 0; i < M; i++) {
      const s = (pi[i] / P0) ** KAPPA;
      phiInterface[K * M + i] = 0;
      for (let k = K - 1; k >= 0; k--) {
        const idx = k * M + i, th = Theta[idx] / pi[i];
        const lower = s * sigmaK[k + 1], upper = s * sigmaK[k], layer = s * (sigma1K[k + 1] - sigma1K[k]) / ((1 + KAPPA) * dSigma);
        theta[idx] = th;
        exnerLayer[idx] = layer;
        phiLayer[idx] = phiInterface[(k + 1) * M + i] + CP * th * (lower - layer);
        phiInterface[k * M + i] = phiLayer[idx] + CP * th * (layer - upper);
      }
    }
  }

  function tendency(s, out, heating, mode, friction) {
    const { pi, Theta, u } = s;
    hydrostatic(pi, Theta);
    for (let k = 0; k < K; k++) for (let i = 0; i < M; i++) {
      const j = (i + 1) % M;
      F[k * M + i] = 0.5 * (pi[i] + pi[j]) * u[k * M + i];
    }
    for (let i = 0; i < M; i++) {
      const h = (i + M - 1) % M;
      let sum = 0;
      for (let k = 0; k < K; k++) { const d = (F[k * M + i] - F[k * M + h]) / dx; D[k * M + i] = d; sum += d * dSigma; }
      out.pi[i] = -sum;
      let acc = 0;
      piSigmaDot[i] = 0;
      for (let k = 0; k < K; k++) { acc -= D[k * M + i] * dSigma; piSigmaDot[(k + 1) * M + i] = acc - (k + 1) * dSigma * out.pi[i]; }
      piSigmaDot[K * M + i] = 0;
    }
    for (let k = 0; k < K; k++) for (let i = 0; i < M; i++) {
      const idx = k * M + i, h = (i + M - 1) % M, j = (i + 1) % M;
      const fluxRight = F[idx] * (F[idx] > 0 ? theta[idx] : theta[k * M + j]);
      const fluxLeft = F[k * M + h] * (F[k * M + h] > 0 ? theta[k * M + h] : theta[idx]);
      const below = piSigmaDot[(k + 1) * M + i], above = piSigmaDot[k * M + i];
      const fluxDown = k === K - 1 ? 0 : below * (below > 0 ? theta[idx] : theta[idx + M]);
      const fluxUp = k === 0 ? 0 : above * (above > 0 ? theta[idx - M] : theta[idx]);
      let q = -(theta[idx] - theta0[k]) / 28800;
      if (Math.abs(i - heated) <= heatedHalf) q += mode === 'ground' ? (k === K - 1 ? heating * K : 0) / exnerLayer[idx] : heating / exnerLayer[idx];
      out.Theta[idx] = -(fluxRight - fluxLeft) / dx - (fluxDown - fluxUp) / dSigma + pi[i] * q;
    }
    for (let k = 0; k < K; k++) for (let i = 0; i < M; i++) {
      const idx = k * M + i, h = (i + M - 1) % M, j = (i + 1) % M;
      const ue = u[idx], piE = 0.5 * (pi[i] + pi[j]);
      const pgf = (phiLayer[k * M + j] - phiLayer[idx]) / dx + CP * 0.5 * (theta[idx] + theta[k * M + j]) * (exnerLayer[k * M + j] - exnerLayer[idx]) / dx;
      const adv = ue > 0 ? ue * (ue - u[k * M + h]) / dx : ue * (u[k * M + j] - ue) / dx;
      const sdUp = k === 0 ? 0 : 0.5 * (piSigmaDot[k * M + i] + piSigmaDot[k * M + j]) / piE;
      const sdDown = k === K - 1 ? 0 : 0.5 * (piSigmaDot[(k + 1) * M + i] + piSigmaDot[(k + 1) * M + j]) / piE;
      const vadv = (Math.max(sdUp, 0) * (ue - (k === 0 ? ue : u[idx - M])) + Math.min(sdDown, 0) * ((k === K - 1 ? ue : u[idx + M]) - ue)) / dSigma;
      const rate = (k === K - 1 ? 1 / 21600 : 0) + (friction ? 1 / 5400 : 0) + sponge[i];
      out.u[idx] = -adv - vadv - pgf - rate * ue;
    }
  }

  function combine(target, base, slope, h) {
    for (let i = 0; i < M; i++) target.pi[i] = base.pi[i] + h * slope.pi[i];
    for (let i = 0; i < n; i++) { target.Theta[i] = base.Theta[i] + h * slope.Theta[i]; target.u[i] = base.u[i] + h * slope.u[i]; }
  }

  function adjust() {
    for (let i = 0; i < M; i++) {
      for (let pass = 0; pass < K; pass++) {
        let mixed = false;
        for (let k = K - 2; k >= 0; k--) {
          const upper = k * M + i, lower = upper + M;
          if (Theta[lower] > Theta[upper] + 1e-9) { const wl = sigma1K[k + 2] - sigma1K[k + 1], wu = sigma1K[k + 1] - sigma1K[k], mean = (wl * Theta[lower] + wu * Theta[upper]) / (wl + wu); Theta[lower] = mean; Theta[upper] = mean; mixed = true; }
        }
        if (!mixed) break;
      }
    }
  }

  function step({ heating = 0, mode = 'column', friction = false } = {}, h = dt) {
    const s0 = { pi, Theta, u };
    tendency(s0, tend[0], heating, mode, friction);
    combine(stages[0], s0, tend[0], h / 3);
    tendency(stages[0], tend[1], heating, mode, friction);
    combine(stages[1], s0, tend[1], h / 2);
    tendency(stages[1], tend[2], heating, mode, friction);
    combine(s0, s0, tend[2], h);
    adjust();
    state.time += h;
    hydrostatic(pi, Theta);
    tendency(s0, tend[0], heating, mode, friction);
  }

  const temperature = (k, i) => theta[k * M + i] * exnerLayer[k * M + i];
  function verticalVelocity(k, i) {
    if (k === 0 || k === K) return 0;
    const h = (i + M - 1) % M, sigma = k * dSigma, p = sigma * pi[i], t = 0.5 * (temperature(k - 1, i) + temperature(k, i));
    const advection = (l) => D[l * M + i] - pi[i] * (u[l * M + i] - u[l * M + h]) / dx;
    const omega = piSigmaDot[k * M + i] + sigma * (tend[0].pi[i] + 0.5 * (advection(k - 1) + advection(k)));
    return -omega * R * t / (p * G);
  }

  reset();
  return {
    M, K, dx, dt, heated, heatedHalf, dSigma, theta0, z0, state, pi, u, theta, exnerLayer, phiLayer, phiInterface, piSigmaDot,
    reset, step, temperature, verticalVelocity,
    interfaceHeight: (k, i) => phiInterface[k * M + i] / G,
    layerHeight: (k, i) => phiLayer[k * M + i] / G,
    initialTemperature: (k) => theta0[k] * (sigma1K[k + 1] - sigma1K[k]) / ((1 + KAPPA) * dSigma),
    mass: () => { let sum = 0; for (let i = 0; i < M; i++) sum += pi[i]; return sum; },
  };
}
