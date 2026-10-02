import { OROGRAPHY_DEFAULTS } from '../physics/orography.module.js';
import { FORM_DRAG_DEFAULTS, formDragScale } from '../physics/formDrag.module.js';

/*
 * The subgrid orography's drag of physics/orography.module.js in WGSL, a
 * line-by-line port: `orography` one thread per column after the
 * boundary layer's diagnosis, `orographyApply` one per edge after the
 * momentum mixing. With the scheme off both kernels return at once.
 */
export function orographyConstants(o) {
  const on = o.orography !== false && o.orography !== undefined;
  const c = { ...OROGRAPHY_DEFAULTS, ...(on ? o.orography : {}) };
  const form = { ...FORM_DRAG_DEFAULTS, ...(o.formDrag || {}) };
  return `
const OROGRAPHY: bool = ${on}; const ORO_CD: f32 = ${c.blockingDrag}; const ORO_G: f32 = ${c.waveDrag}; const ORO_HN: f32 = ${c.criticalHeight}; const ORO_RIC: f32 = ${c.criticalRichardson}; const ORO_HEFF: f32 = ${c.effectiveHeight};
const FORM_DRAG: bool = ${!!o.formDrag}; const TOFD_SCALE: f32 = ${formDragScale(form)}; const TOFD_DECAY: f32 = ${form.decayHeight};
`;
}

export const OROGRAPHY_KERNELS = {
  orography: `fn oroN2(k: i32, z: ptr<function, array<f32, K>>, th: ptr<function, array<f32, K>>) -> f32 {
  return GRAV * ((*th)[k - 1] - (*th)[k]) / (0.5 * ((*th)[k - 1] + (*th)[k]) * ((*z)[k - 1] - (*z)[k]));
}
fn oroInterface(k: i32, z: ptr<function, array<f32, K>>) -> f32 {
  if (k >= K) { return 0.0; }
  if (k <= 0) { return 2.0 * (*z)[0] - 0.5 * ((*z)[0] + (*z)[1]); }
  return 0.5 * ((*z)[k] + (*z)[k - 1]);
}
@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= C) { return; }
  if (!OROGRAPHY) { return; }
  for (var k = 0; k < K; k++) { PH[PH_OBETA + k * C + i] = 0.0; PH[PH_OWAVE + k * C + i] = 0.0; }
  PH[PH_ODIR + 3 * i] = 0.0; PH[PH_ODIR + 3 * i + 1] = 0.0; PH[PH_ODIR + 3 * i + 2] = 0.0;
  PH[PH_OBLOCK + i] = 0.0; PH[PH_OLAUNCH + i] = 0.0;
  let mu = PH[PH_OSTD + i]; let sigma = PH[PH_OSLP + i];
  if (!(mu > 0.0) || !(sigma > 0.0)) { return; }
  let gamma = PH[PH_OANI + i]; let theta0 = PH[PH_OORI + i];
  let pi = IN[S_PI + i]; let dt = P[0]; let bottom = K - 1;
  let xc = vec3<f32>(MF[F_XC + 3 * i], MF[F_XC + 3 * i + 1], MF[F_XC + 3 * i + 2]);
  let ex = vec3<f32>(-xc.y, xc.x, 0.0) / length(xc.xy);
  let ny = cross(xc, ex);
  var z: array<f32, K>; var p: array<f32, K>; var rho: array<f32, K>; var th: array<f32, K>; var ue: array<f32, K>; var vn: array<f32, K>;
  for (var k = 0; k < K; k++) {
    let idx = k * C + i;
    z[k] = (D[D_GEO + idx] + LV[L_GABS + k]) / GRAV;
    p[k] = pi * LV[L_SM + k];
    th[k] = IN[S_TH + idx];
    rho[k] = p[k] / (RGAS * th[k] * D[D_EXM + idx]);
    let w = cellWind(i, k);
    ue[k] = dot(w, ex); vn[k] = dot(w, ny);
  }
  var weight = 0.0; var uL = 0.0; var vL = 0.0; var rhoL = 0.0; var nL2 = 0.0;
  for (var k = bottom; k >= 0; k--) {
    let lo = oroInterface(k + 1, &z); let hi = oroInterface(k, &z);
    let overlap = min(hi, 2.0 * mu) - max(lo, mu);
    if (overlap > 0.0) {
      let w = rho[k] * overlap;
      var nk = 0.0;
      if (k == 0) { nk = oroN2(1, &z, &th); } else if (k == bottom) { nk = oroN2(bottom, &z, &th); } else { nk = 0.5 * (oroN2(k, &z, &th) + oroN2(k + 1, &z, &th)); }
      weight += w; uL += w * ue[k]; vL += w * vn[k]; rhoL += w * rho[k]; nL2 += w * nk;
    }
    if (lo >= 2.0 * mu) { break; }
  }
  uL /= weight; vL /= weight; rhoL /= weight; nL2 /= weight;
  let speedL = sqrt(uL * uL + vL * vL);
  if (!(speedL > 1e-3)) { return; }
  let ax = uL / speedL; let ay = vL / speedL;
  let B = 1.0 - 0.18 * gamma - 0.04 * gamma * gamma; let Cc = 0.48 * gamma + 0.3 * gamma * gamma;
  let psiL = theta0 - atan2(vL, uL); let sL = sin(psiL); let cL = cos(psiL);
  let D1 = B * cL * cL + Cc * sL * sL; let D2 = (B - Cc) * sL * cL; let Dn = sqrt(D1 * D1 + D2 * D2);
  let tx = (D1 * ax - D2 * ay) / Dn; let ty = (D1 * ay + D2 * ax) / Dn;
  let dir = tx * ex + ty * ny;
  PH[PH_ODIR + 3 * i] = dir.x; PH[PH_ODIR + 3 * i + 1] = dir.y; PH[PH_ODIR + 3 * i + 2] = dir.z;
  var along: array<f32, K>; var plane: array<f32, K>;
  for (var k = 0; k < K; k++) { along[k] = ue[k] * ax + vn[k] * ay; plane[k] = (along[k] * D1 + (vn[k] * ax - ue[k] * ay) * D2) / Dn; }
  let top = 3.0 * mu;
  var blocking = 0.0; var integral = 0.0;
  var kStart = bottom;
  while (kStart > 0 && z[kStart - 1] < top) { kStart--; }
  var upper = top;
  for (var k = kStart; k <= bottom + 1; k++) {
    var lower = 0.0;
    if (k <= bottom) { lower = z[k]; }
    if (lower >= upper) { continue; }
    let kk = max(1, min(bottom, k));
    let nk = sqrt(max(0.0, oroN2(kk, &z, &th)));
    var up = along[bottom];
    if (k <= bottom) { if (k > 0) { up = 0.5 * (along[k] + along[k - 1]); } else { up = along[0]; } }
    if (!(up > 0.0)) { blocking = upper; break; }
    let step = nk / up * (upper - lower);
    if (integral + step >= ORO_HN) { blocking = upper - (ORO_HN - integral) / (nk / up); break; }
    integral += step;
    upper = lower;
  }
  PH[PH_OBLOCK + i] = blocking;
  for (var k = bottom; k >= 0; k--) {
    if (!(z[k] < blocking)) { break; }
    let psi = theta0 - atan2(vn[k], ue[k]); let s = sin(psi); let c = cos(psi);
    let across = c * c + gamma * s * s;
    var shape = 0.0;
    if (across > 0.0) { shape = max(2.0 - (gamma * c * c + s * s) / across, 0.0); }
    let speed = sqrt(ue[k] * ue[k] + vn[k] * vn[k]);
    PH[PH_OBETA + k * C + i] = ORO_CD * shape * sigma / (2.0 * mu) * sqrt((blocking - z[k]) / (z[k] + mu)) * (B * c * c + Cc * s * s) * speed / 2.0;
  }
  let height = ORO_HEFF * (top - blocking);
  if (!(nL2 > 0.0) || !(height > 0.0)) { return; }
  let nL = sqrt(nL2);
  let launch = rhoL * height * height / 9.0 * sigma / mu * ORO_G * speedL * Dn * nL;
  PH[PH_OLAUNCH + i] = launch;
  let flux0 = rhoL * nL * speedL * D1 / Dn;
  var tau: array<f32, K + 1>;
  for (var k = 0; k <= K; k++) { tau[k] = 0.0; }
  tau[K] = launch;
  var first = bottom;
  while (first > 0 && oroInterface(first, &z) <= blocking) { tau[first] = launch; first--; }
  var breakTop = blocking; var phase = 0.0;
  for (var k = bottom; k >= 1; k--) {
    if (z[k - 1] <= blocking) { continue; }
    let lower = max(z[k], blocking); let nk = sqrt(max(0.0, oroN2(k, &z, &th))); let up = 0.5 * (along[k] + along[k - 1]);
    if (!(up > 0.0)) { breakTop = z[k - 1]; break; }
    let step = nk / up * (z[k - 1] - lower);
    if (phase + step >= 1.5707963) { breakTop = lower + (1.5707963 - phase) / (nk / up); break; }
    phase += step; breakTop = z[k - 1];
  }
  breakTop = max(breakTop, 4.0 * mu);
  var current = launch;
  for (var k = first; k >= 1; k--) {
    let N2 = oroN2(k, &z, &th); let V = 0.5 * (plane[k] + plane[k - 1]); let rhoI = 0.5 * (rho[k] + rho[k - 1]);
    if (!(V > 0.0) || !(current > 0.0) || !(N2 > 0.0)) { current = 0.0; tau[k] = 0.0; continue; }
    let N = sqrt(N2); let S = abs(plane[k - 1] - plane[k]) / (z[k - 1] - z[k]);
    let dz = height * sqrt(current / launch * flux0 / (rhoI * N * V)); let alpha = N * dz / V;
    if (N2 * (1.0 - alpha) < ORO_RIC * (S + N * alpha) * (S + N * alpha)) {
      let critical = (-(2.0 * ORO_RIC * S * N + N2) + pow(N, 1.5) * sqrt(N * (1.0 + 4.0 * ORO_RIC) + 4.0 * ORO_RIC * S)) / (2.0 * ORO_RIC * N2);
      current = select(0.0, current * (critical / alpha) * (critical / alpha), critical > 0.0);
    }
    tau[k] = current;
  }
  tau[0] = 0.0;
  var kb = first;
  while (kb > 0 && oroInterface(kb, &z) <= breakTop) { kb--; }
  if (kb < first && tau[kb] < launch) {
    var pBlock = pi * LV[L_SL + bottom];
    if (blocking > 0.0) {
      if (blocking <= z[bottom]) { pBlock = exp(log(pBlock) + blocking / z[bottom] * (log(p[bottom]) - log(pBlock))); }
      else {
        var k = bottom; while (k > 0 && z[k - 1] < blocking) { k--; }
        if (k == 0) { pBlock = p[0]; } else { pBlock = exp(log(p[k]) + (blocking - z[k]) / (z[k - 1] - z[k]) * (log(p[k - 1]) - log(p[k]))); }
      }
    }
    let pBreak = pi * LV[L_SU + kb];
    for (var k = first; k > kb; k--) { tau[k] = launch + (tau[kb] - launch) * clamp((pi * LV[L_SU + k] - pBlock) / (pBreak - pBlock), 0.0, 1.0); }
  }
  var carried = 0.0;
  for (var k = bottom; k >= 0; k--) {
    let mass = pi * LV[L_DS + k] / GRAV;
    var accel = (tau[k + 1] - tau[k]) / mass + carried / mass;
    carried = 0.0;
    if (dt > 0.0) {
      let limit = max(0.0, ue[k] * tx + vn[k] * ty) / dt;
      if (accel > limit) { carried = (accel - limit) * mass; accel = limit; }
    }
    PH[PH_OWAVE + k * C + i] = -accel;
  }
}`,
  orographyApply: `@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let e = i32(id.x); if (e >= E) { return; }
  if (!OROGRAPHY) { return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1]; let dt = P[0];
  let n = vec3<f32>(MF[F_NEDGE + 3 * e], MF[F_NEDGE + 3 * e + 1], MF[F_NEDGE + 3 * e + 2]);
  let na = dot(vec3<f32>(PH[PH_ODIR + 3 * a], PH[PH_ODIR + 3 * a + 1], PH[PH_ODIR + 3 * a + 2]), n);
  let nb = dot(vec3<f32>(PH[PH_ODIR + 3 * b], PH[PH_ODIR + 3 * b + 1], PH[PH_ODIR + 3 * b + 2]), n);
  let columnMass = 0.5 * (IN[S_PI + a] + IN[S_PI + b]);
  var lost = 0.0;
  for (var k = 0; k < K; k++) {
    let rate = 0.5 * (PH[PH_OBETA + k * C + a] + PH[PH_OBETA + k * C + b]);
    let push = 0.5 * (PH[PH_OWAVE + k * C + a] * na + PH[PH_OWAVE + k * C + b] * nb);
    if (rate == 0.0 && push == 0.0) { continue; }
    let idx = k * E + e; let before = IN[S_U + idx]; let after = (before + dt * push) / (1.0 + dt * rate);
    IN[S_U + idx] = after;
    lost += columnMass * LV[L_DS + k] / GRAV * (before - after);
    D[D_DISS + idx] += before * before - after * after;
  }
  PH[PH_OSTRESS + e] = lost / dt;
}`,
};
