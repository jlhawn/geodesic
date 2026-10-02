import { cellVector } from '../dynamics/operators.module.js';

/*
 * The drag of the mountains the mesh does not resolve: Lott and Miller
 * (1997) as the IFS documents it (Cy47r3 Part IV, Chapter 4), from the
 * subgrid orography of each cell (geography.module.js's
 * subgridOrography: standard deviation μ, anisotropy γ, orientation θ,
 * slope σ).
 *
 * Per column, from the state after the physics: the incident flow U_L, its
 * density ρ_L and Brunt–Väisälä frequency N_L (N² = g Δθ/(θ̄ Δz) between
 * layer midpoints, IFS eq. 4.26) are the layer-mass means over μ < z < 2μ
 * above the model's ground (4.27). ψ is the angle from a wind to the
 * principal axis, B = 1 − 0.18γ − 0.04γ², C = 0.48γ + 0.3γ² (Phillips
 * 1984), D1 = B cos²ψ_L + C sin²ψ_L, D2 = (B − C) sin ψ_L cos ψ_L.
 *  - The blocking height Z_b is the highest level below 3μ where
 *    ∫ N/U_p dz from it to 3μ reaches criticalHeight (H_n,crit), U_p the
 *    wind along U_L, and a level where U_p ≤ 0 blocks everything below it
 *    (4.9, 4.32). Each layer whose midpoint z lies below Z_b feels
 *      ∂u/∂t = −C_d max(2 − 1/r, 0) (σ/2μ) √((Z_b − z)/(z + μ))
 *               (B cos²ψ + C sin²ψ) |U| u / 2
 *    with r = (cos²ψ + γ sin²ψ)/(γ cos²ψ + sin²ψ) (4.14, 4.40), solved
 *    implicitly with |U| at the start of the step (4.41).
 *  - The gravity waves launched above it carry the stress
 *      τ_0 = ρ_L (H_eff²/9) (σ/μ) G |U_L| √(D1² + D2²) N_L,
 *    H_eff = effectiveHeight × (3μ − Z_b) (4.37, which for Z_b = 0 and
 *    effectiveHeight 1 is LM97's eq. with H = 2μ),
 *    along the direction (D1, D2) in the frame of U_L (D2 towards the
 *    left of the flow), constant from the ground to Z_b. Up the column the
 *    stress stays until the wave Richardson number
 *      Ri_w = N² (1 − α)/(S + N α)²,  α = N δz / V,
 *    with V the wind in the plane of the stress (4.31), S = |∂V/∂z| and
 *    the displacement δz from ρ N V δz² ∝ τ (δz = H_eff at launch), falls
 *    below criticalRichardson, where δz is cut to the amplitude that
 *    holds it there and τ with it; V ≤ 0 is a critical level, where all
 *    of it goes. Where it breaks below Z_b + Δz, the depth with
 *    ∫ N/U_p dz = π/2 above Z_b and at least 4μ, the stress falls
 *    linearly in pressure over that depth (4.33, 4.35); what reaches the
 *    model top goes into the top layer. A layer takes
 *    ∂u/∂t = −g ∂τ/∂p along the stress, at most what stops its own wind
 *    along it within the step, the rest passing to the layer above.
 *
 * Constants (Lott and Miller 1997): blockingDrag C_d 1, waveDrag G 1,
 * criticalHeight 0.5, criticalRichardson 0.25. The IFS Cy47r3 documents
 * C_d 2 and H_eff doubled (from Cy32r2) with its own G, not taken here.
 *
 * `diagnose` lays the per-cell rates, `apply` steps each edge's normal
 * velocity as u ← (u + Δt a)/(1 + Δt β), with a and β the two cells'
 * means projected on the edge, counts the kinetic energy it removes into
 * each layer's dissipation and keeps in `stress` the momentum the column
 * lost per unit area and time, the stress on the ground along the edge's
 * normal (N/m²).
 */
export const OROGRAPHY_DEFAULTS = { blockingDrag: 1, waveDrag: 1, criticalHeight: 0.5, criticalRichardson: 0.25, effectiveHeight: 1 };

/*
 * One column, k = 0 at the top: z (above the ground), p, rho, theta, east
 * and north wind at the K layer midpoints, pTop and pBottom each layer's
 * interface pressures. Writes the blocking rate out.beta (s⁻¹), the wave
 * acceleration out.wave (m/s², along out.direction, the stress's unit
 * vector east and north), the blocking height, the launched stress and
 * the stress the column takes; with dt > 0 the wave's limit per step.
 */
export function orographicColumn(sub, column, options = OROGRAPHY_DEFAULTS, out = null, dt = 0) {
  const { deviation: mu, anisotropy: gamma, orientation: theta0, slope: sigma } = sub;
  const { z, p, rho, theta, east, north, pTop, pBottom, g = 9.80616 } = column;
  const K = z.length;
  const o = { ...OROGRAPHY_DEFAULTS, ...options };
  out ??= { beta: new Float64Array(K), wave: new Float64Array(K), direction: [0, 0], blocking: 0, launch: 0, deposited: 0 };
  out.beta.fill(0); out.wave.fill(0); out.blocking = 0; out.launch = 0; out.direction[0] = 0; out.direction[1] = 0;
  if (!(mu > 0) || !(sigma > 0)) return out;
  const bottom = K - 1;
  const zi = (k) => (k >= K ? 0 : k <= 0 ? 2 * z[0] - 0.5 * (z[0] + z[1]) : 0.5 * (z[k] + z[k - 1]));
  const n2 = (k) => g * (theta[k - 1] - theta[k]) / (0.5 * (theta[k - 1] + theta[k]) * (z[k - 1] - z[k]));
  let weight = 0, uL = 0, vL = 0, rhoL = 0, nL2 = 0;
  for (let k = bottom; k >= 0; k--) {
    const lo = zi(k + 1), hi = zi(k), overlap = Math.min(hi, 2 * mu) - Math.max(lo, mu);
    if (overlap > 0) {
      const w = rho[k] * overlap;
      const nk = k === 0 ? n2(1) : k === bottom ? n2(bottom) : 0.5 * (n2(k) + n2(k + 1));
      weight += w; uL += w * east[k]; vL += w * north[k]; rhoL += w * rho[k]; nL2 += w * nk;
    }
    if (lo >= 2 * mu) break;
  }
  uL /= weight; vL /= weight; rhoL /= weight; nL2 /= weight;
  const speedL = Math.hypot(uL, vL);
  if (!(speedL > 1e-3)) return out;
  const ax = uL / speedL, ay = vL / speedL;
  const B = 1 - 0.18 * gamma - 0.04 * gamma * gamma, Cc = 0.48 * gamma + 0.3 * gamma * gamma;
  const psiL = theta0 - Math.atan2(vL, uL), sL = Math.sin(psiL), cL = Math.cos(psiL);
  const D1 = B * cL * cL + Cc * sL * sL, D2 = (B - Cc) * sL * cL, D = Math.hypot(D1, D2);
  const tx = (D1 * ax - D2 * ay) / D, ty = (D1 * ay + D2 * ax) / D;
  out.direction[0] = tx; out.direction[1] = ty;
  const along = (k) => east[k] * ax + north[k] * ay;
  const plane = (k) => (along(k) * D1 + (north[k] * ax - east[k] * ay) * D2) / D;
  const top = 3 * mu;
  let blocking = 0, integral = 0;
  {
    let k = bottom;
    while (k > 0 && z[k - 1] < top) k--;
    let upper = top;
    for (; k <= bottom + 1; k++) {
      const lower = k <= bottom ? z[k] : 0;
      if (lower >= upper) continue;
      const kk = Math.max(1, Math.min(bottom, k));
      const nk = Math.sqrt(Math.max(0, n2(kk)));
      const up = k <= bottom ? (k > 0 ? 0.5 * (along(k) + along(k - 1)) : along(0)) : along(bottom);
      if (!(up > 0)) { blocking = upper; break; }
      const step = nk / up * (upper - lower);
      if (integral + step >= o.criticalHeight) { blocking = upper - (o.criticalHeight - integral) / (nk / up); break; }
      integral += step;
      upper = lower;
    }
  }
  out.blocking = blocking;
  for (let k = bottom; k >= 0 && z[k] < blocking; k--) {
    const psi = theta0 - Math.atan2(north[k], east[k]), s = Math.sin(psi), c = Math.cos(psi);
    const across = c * c + gamma * s * s, shape = across > 0 ? Math.max(2 - (gamma * c * c + s * s) / across, 0) : 0;
    const speed = Math.hypot(east[k], north[k]);
    out.beta[k] = o.blockingDrag * shape * sigma / (2 * mu) * Math.sqrt((blocking - z[k]) / (z[k] + mu)) * (B * c * c + Cc * s * s) * speed / 2;
  }
  const height = o.effectiveHeight * (top - blocking);
  if (!(nL2 > 0) || !(height > 0)) return out;
  const nL = Math.sqrt(nL2);
  const launch = rhoL * height * height / 9 * sigma / mu * o.waveDrag * speedL * D * nL;
  out.launch = launch;
  const vLaunch = speedL * D1 / D, flux0 = rhoL * nL * vLaunch;
  const tau = new Float64Array(K + 1);
  tau[K] = launch;
  let first = bottom;
  while (first > 0 && zi(first) <= blocking) { tau[first] = launch; first--; }
  let breakTop = blocking, phase = 0;
  for (let k = bottom; k >= 1; k--) {
    if (z[k - 1] <= blocking) continue;
    const lower = Math.max(z[k], blocking), nk = Math.sqrt(Math.max(0, n2(k))), up = 0.5 * (along(k) + along(k - 1));
    if (!(up > 0)) { breakTop = z[k - 1]; break; }
    const step = nk / up * (z[k - 1] - lower);
    if (phase + step >= Math.PI / 2) { breakTop = lower + (Math.PI / 2 - phase) / (nk / up); break; }
    phase += step; breakTop = z[k - 1];
  }
  breakTop = Math.max(breakTop, 4 * mu);
  let current = launch;
  for (let k = first; k >= 1; k--) {
    const N2 = n2(k), V = 0.5 * (plane(k) + plane(k - 1)), rhoI = 0.5 * (rho[k] + rho[k - 1]);
    if (!(V > 0) || !(current > 0)) { current = 0; tau[k] = 0; continue; }
    if (!(N2 > 0)) { current = 0; tau[k] = 0; continue; }
    const N = Math.sqrt(N2), S = Math.abs(plane(k - 1) - plane(k)) / (z[k - 1] - z[k]);
    const dz = height * Math.sqrt(current / launch * flux0 / (rhoI * N * V)), alpha = N * dz / V;
    const rc = o.criticalRichardson;
    if (N2 * (1 - alpha) < rc * (S + N * alpha) ** 2) {
      const critical = (-(2 * rc * S * N + N2) + Math.pow(N, 1.5) * Math.sqrt(N * (1 + 4 * rc) + 4 * rc * S)) / (2 * rc * N2);
      current = critical > 0 ? current * (critical / alpha) ** 2 : 0;
    }
    tau[k] = current;
  }
  tau[0] = 0;
  const pAt = (height) => {
    if (height <= 0) return pBottom[bottom];
    if (height <= z[bottom]) return Math.exp(Math.log(pBottom[bottom]) + height / z[bottom] * (Math.log(p[bottom]) - Math.log(pBottom[bottom])));
    let k = bottom; while (k > 0 && z[k - 1] < height) k--;
    if (k === 0) return p[0];
    return Math.exp(Math.log(p[k]) + (height - z[k]) / (z[k - 1] - z[k]) * (Math.log(p[k - 1]) - Math.log(p[k])));
  };
  let kb = first;
  while (kb > 0 && zi(kb) <= breakTop) kb--;
  if (kb < first && tau[kb] < launch) {
    const pBlock = pAt(blocking), pBreak = pTop[kb];
    for (let k = first; k > kb; k--) tau[k] = launch + (tau[kb] - launch) * Math.min(1, Math.max(0, (pTop[k] - pBlock) / (pBreak - pBlock)));
  }
  let carried = 0;
  for (let k = bottom; k >= 0; k--) {
    const mass = (pBottom[k] - pTop[k]) / g;
    let accel = (tau[k + 1] - tau[k]) / mass + carried / mass;
    const limit = dt > 0 ? Math.max(0, east[k] * tx + north[k] * ty) / dt : Infinity;
    carried = 0;
    if (accel > limit) { carried = (accel - limit) * mass; accel = limit; }
    out.wave[k] = -accel;
  }
  out.deposited = launch - carried;
  return out;
}

export function createOrographicDrag(mesh, core, sub, { buffers = null, ...options } = {}) {
  const o = { ...OROGRAPHY_DEFAULTS, ...options };
  const { K, C, E, dSigma, sigmaMid, levels, R, g, exnerLayer, geopotential } = core.diagnostics;
  const { cellsOnEdge, nEdge, lonCell, latCell } = mesh;
  const shared = (name, n) => (buffers && buffers[name] ? buffers[name] : new SharedArrayBuffer(8 * n));
  const betaBuffer = shared('beta', K * C), waveBuffer = shared('wave', K * C), directionBuffer = shared('direction', 3 * C), stressBuffer = shared('stress', E);
  const beta = new Float64Array(betaBuffer), wave = new Float64Array(waveBuffer), direction = new Float64Array(directionBuffer), stress = new Float64Array(stressBuffer);
  const vector = new Float64Array(3 * K * C), layer = new Float64Array(3 * C);
  const column = { z: new Float64Array(K), p: new Float64Array(K), rho: new Float64Array(K), theta: new Float64Array(K), east: new Float64Array(K), north: new Float64Array(K), pTop: new Float64Array(K), pBottom: new Float64Array(K), g };
  const result = { beta: new Float64Array(K), wave: new Float64Array(K), direction: [0, 0], blocking: 0, launch: 0, deposited: 0 };
  const blocking = new Float64Array(C), launch = new Float64Array(C);

  function diagnose(state, iFrom, iTo, dt) {
    const [pi, theta, u] = state;
    for (let k = 0; k < K; k++) { cellVector(mesh, u.subarray(k * E, (k + 1) * E), layer, iFrom, iTo); for (let i = 3 * iFrom; i < 3 * iTo; i++) vector[3 * k * C + i] = layer[i]; }
    for (let i = iFrom; i < iTo; i++) {
      for (let k = 0; k < K; k++) { beta[k * C + i] = 0; wave[k * C + i] = 0; }
      direction[3 * i] = 0; direction[3 * i + 1] = 0; direction[3 * i + 2] = 0;
      blocking[i] = 0; launch[i] = 0;
      if (!(sub.deviation[i] > 0)) continue;
      const lon = lonCell[i], lat = latCell[i];
      const ex = [-Math.sin(lon), Math.cos(lon), 0], ny = [-Math.sin(lat) * Math.cos(lon), -Math.sin(lat) * Math.sin(lon), Math.cos(lat)];
      const ground = geopotential[(K - 1) * C + i] - core.diagnostics.cp * core.arrays.thetaV[(K - 1) * C + i] * (core.diagnostics.exnerLower[(K - 1) * C + i] - exnerLayer[(K - 1) * C + i]);
      for (let k = 0; k < K; k++) {
        const idx = k * C + i, v = 3 * idx;
        column.z[k] = (geopotential[idx] - ground) / g;
        column.p[k] = pi[i] * sigmaMid[k];
        column.pTop[k] = pi[i] * levels[k]; column.pBottom[k] = pi[i] * levels[k + 1];
        column.theta[k] = theta[idx];
        column.rho[k] = column.p[k] / (R * theta[idx] * exnerLayer[idx]);
        column.east[k] = vector[v] * ex[0] + vector[v + 1] * ex[1];
        column.north[k] = vector[v] * ny[0] + vector[v + 1] * ny[1] + vector[v + 2] * ny[2];
      }
      orographicColumn({ deviation: sub.deviation[i], anisotropy: sub.anisotropy[i], orientation: sub.orientation[i], slope: sub.slope[i] }, column, o, result, dt);
      for (let k = 0; k < K; k++) { beta[k * C + i] = result.beta[k]; wave[k * C + i] = result.wave[k]; }
      const [tx, ty] = result.direction;
      for (let c = 0; c < 3; c++) direction[3 * i + c] = tx * ex[c] + ty * ny[c];
      blocking[i] = result.blocking; launch[i] = result.launch;
    }
  }

  function apply(pi, u, eFrom, eTo, dt, dissipation = null) {
    for (let e = eFrom; e < eTo; e++) {
      const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1];
      const na = direction[3 * a] * nEdge[3 * e] + direction[3 * a + 1] * nEdge[3 * e + 1] + direction[3 * a + 2] * nEdge[3 * e + 2];
      const nb = direction[3 * b] * nEdge[3 * e] + direction[3 * b + 1] * nEdge[3 * e + 1] + direction[3 * b + 2] * nEdge[3 * e + 2];
      const columnMass = 0.5 * (pi[a] + pi[b]);
      let lost = 0;
      for (let k = 0; k < K; k++) {
        const rate = 0.5 * (beta[k * C + a] + beta[k * C + b]), push = 0.5 * (wave[k * C + a] * na + wave[k * C + b] * nb);
        if (rate === 0 && push === 0) continue;
        const idx = k * E + e, before = u[idx], after = (before + dt * push) / (1 + dt * rate);
        u[idx] = after;
        lost += columnMass * dSigma[k] / g * (before - after);
        if (dissipation) dissipation[idx] += before * before - after * after;
      }
      stress[e] = lost / dt;
    }
  }

  return { diagnose, apply, beta, wave, direction, stress, blocking, launch, fields: sub, options: o, shared: { beta: betaBuffer, wave: waveBuffer, direction: directionBuffer, stress: stressBuffer } };
}
