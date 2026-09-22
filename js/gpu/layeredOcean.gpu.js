import { emptyBuffer, readBuffer } from './device.module.js';
import { LAYER_DENSITIES, LAYER_BOTTOMS, bathymetryFrom } from '../ocean/layered.module.js';
import { FREEZING_POINT } from '../physics/ice.module.js';

/*
 * The layered ocean of ocean/layered.module.js on the GPU: a bulk mixed
 * layer over interior isopycnal layers, each a TRiSK shallow-water layer
 * carrying thickness, edge velocity, heat h·T and salt h·S. The RK4
 * baroclinic state (h, u, h·T, h·S for all L layers) lives in its own
 * ping-pong buffer set (S/T/K1-K4), exactly the shape of the two-layer
 * ocean.gpu.js generalized to L layers; everything else the step needs
 * (edge thicknesses, PV, the interior potential, the ∇⁴ closure scratch,
 * the barotropic sub-stepping state, the surface staging) lives in one
 * big scratch/static buffer (OD) bound alongside it. The barotropic mode
 * takes its own small RK4 with fixed-offset stage blocks inside OD,
 * looped on the host the M times a step needs; every other multi-stage
 * loop (the baroclinic RK4, the ∇⁴ closure's two Laplacian passes) mirrors
 * the dispatch-ordering-is-execution-ordering pattern the sigma core and
 * the two-layer ocean already rely on within one compute pass.
 *
 * Initialization, loading and serialization run once per model build and
 * are cheap relative to a step, so they stay in JavaScript at double
 * precision, exactly porting the CPU functions over plain arrays before
 * uploading; only the per-step dynamics (tendency, barotropic, mixed
 * layer, salt, surface write-back) are WGSL.
 */
const WORKGROUP = 64;
const EPS = 0.01, THIN = 5, PV_FLOOR = 20, SPEED_LIMIT = 5, DENSITY_TOLERANCE = 0.005;

export const OCEAN_DEFAULTS = {
  densities: LAYER_DENSITIES, bottoms: LAYER_BOTTOMS, mixedDepth: 60, minimumDepth: 50, flatDepth: 4000, thermoclineTilt: 0.3,
  density: 1025, specificHeat: 3985, thermalExpansion: 2e-4, halineContraction: 7.6e-4, referenceT: 283.15, referenceS: 35, gravity: 9.81,
  minimumThickness: 20, shallowestMixedDepth: 50, stirringDepth: 100, maximumMixedDepth: 200, stirring: 0.8, detrainmentTime: 86400, iceSalinity: 5, iceDensity: 917,
  interfacialDrag: 2e-4, bottomDrag: 3e-3, closureHours: 12, diffusivity: 0.3, everySteps: 4,
  dragCoefficient: 1.5e-3, gustiness: 3,
};
const defaultSalinityProfile = (lat) => 34.5 + 1.5 * Math.exp(-(((Math.abs(lat) * 180 / Math.PI - 25) / 15) ** 2));

function seq(names) { const out = {}; let off = 0; for (const [name, n] of names) { out[name] = off; off += n; } out.total = off; return out; }

function oceanKernels(o) {
  const { L, C, E, V, OS, OD, B } = o;
  const constLine = (name, value) => `const ${name}: f32 = ${Number(value).toExponential(10)};`;
  const rhoLine = `const RHO: array<f32, ${L}> = array<f32, ${L}>(${o.rho.map((v) => v.toFixed(6)).join(', ')});`;
  const labelLine = `const LABEL_T: array<f32, ${L}> = array<f32, ${L}>(${o.labelT.map((v) => v.toFixed(6)).join(', ')});`;
  const offsetLines = Object.entries(OD).filter(([k]) => k !== 'total').map(([k, v]) => `const O_${k}: i32 = ${v};`).join('\n');
  const bLines = Object.entries(B).filter(([k]) => k !== 'total').map(([k, v]) => `const B_${k}: i32 = ${v};`).join('\n');
  const head = `
const L: i32 = ${L};
const OH: i32 = ${OS.OH}; const OU: i32 = ${OS.OU}; const OQ: i32 = ${OS.OQ}; const OW: i32 = ${OS.OW}; const OSTOTAL: i32 = ${OS.total};
${offsetLines}
${bLines}
${rhoLine}
${labelLine}
${constLine('RHO0', o.density)} ${constLine('RHOCP', o.density * o.specificHeat)}
${constLine('THERMAL_EXP', o.thermalExpansion)} ${constLine('HALINE_CONTRACT', o.halineContraction)}
${constLine('REF_T', o.referenceT)} ${constLine('REF_S', o.referenceS)} ${constLine('OGRAV', o.gravity)}
${constLine('EPSO', EPS)} ${constLine('THINO', THIN)} ${constLine('PVFLOOR', PV_FLOOR)} ${constLine('SPEEDLIM', SPEED_LIMIT)} ${constLine('DENSTOL', DENSITY_TOLERANCE)}
${constLine('MINTHICK', o.minimumThickness)} ${constLine('SHALLOWMIXED', o.shallowestMixedDepth)} ${constLine('MAXMIXED', o.maximumMixedDepth)}
${constLine('STIRRING', o.stirring)} ${constLine('STIRDEPTH', o.stirringDepth)} ${constLine('DETRAINT', o.detrainmentTime)} ${constLine('ICESAL', o.iceSalinity)} ${constLine('ICEDENS', o.iceDensity)}
${constLine('RINT', o.interfacialDrag)} ${constLine('RBOT', o.bottomDrag)} ${constLine('NU4O', o.nu4)} ${constLine('DIFFUSION', o.diffusion)}
${constLine('FREEZE', FREEZING_POINT)} ${constLine('CDO', o.dragCoefficient)} ${constLine('GUSTO', o.gustiness)}
@group(0) @binding(0) var<storage, read_write> MI: array<i32>;
@group(0) @binding(1) var<storage, read_write> MF: array<f32>;
@group(0) @binding(2) var<storage, read_write> LV: array<f32>;
@group(0) @binding(3) var<storage, read_write> IN: array<f32>;
@group(0) @binding(4) var<storage, read_write> OUT: array<f32>;
@group(0) @binding(5) var<storage, read_write> OD: array<f32>;
@group(0) @binding(6) var<storage, read_write> P: array<f32>;
@group(0) @binding(7) var<storage, read_write> PH: array<f32>;
@group(0) @binding(8) var<storage, read_write> S: array<f32>;
@group(0) @binding(9) var<storage, read_write> D: array<f32>;
fn hOff(k: i32) -> i32 { return OH + k * C; }
fn uOff(k: i32) -> i32 { return OU + k * E; }
fn qOff(k: i32) -> i32 { return OQ + k * C; }
fn wOff(k: i32) -> i32 { return OW + k * C; }
fn eos(t: f32, s: f32) -> f32 { return RHO0 * (1.0 - THERMAL_EXP * (t - REF_T) + HALINE_CONTRACT * (s - REF_S)); }
fn detrain(i: i32, amount: f32, rm: f32) {
  if (amount <= 0.0) { return; }
  var k = 1;
  for (var j = 2; j < L; j++) { if (abs(RHO[j] - rm) < abs(RHO[k] - rm)) { k = j; } }
  moveLayer(i, 0, k, amount);
}
fn moveLayer(i: i32, srcK: i32, dstK: i32, amount: f32) {
  let ha = hOff(srcK) + i; let hb = hOff(dstK) + i;
  let qa = qOff(srcK) + i; let qb = qOff(dstK) + i;
  let wa = wOff(srcK) + i; let wb = wOff(dstK) + i;
  let f = amount / IN[ha];
  let dq = IN[qa] * f; let dw = IN[wa] * f;
  IN[ha] -= amount; IN[qa] -= dq; IN[wa] -= dw;
  IN[hb] += amount; IN[qb] += dq; IN[wb] += dw;
}
`;
  const idx = `(i32(id.x) + i32(id.y) * 4194240)`;
  const K = `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {\n`;
  return {
    head,
    oFreeSurface: `${K}  let i = ${idx}; if (i >= C) { return; }
  var sum = 0.0;
  for (var k = 0; k < L; k++) { sum += IN[hOff(k) + i]; }
  OD[O_ETA + i] = select(0.0, sum - OD[O_BATH + i], OD[O_CMASK + i] > 0.5);
}`,
    oGradEta: `${K}  let e = ${idx}; if (e >= E) { return; }
  OD[O_GRADETA + e] = (OD[O_ETA + MI[COE + 2 * e + 1]] - OD[O_ETA + MI[COE + 2 * e]]) / MF[F_DC + e];
}`,
    oSurfaceDensity: `${K}  let i = ${idx}; if (i >= C) { return; }
  let h0 = max(EPSO, IN[hOff(0) + i]);
  OD[O_RHOML + i] = eos(IN[qOff(0) + i] / h0, IN[wOff(0) + i] / h0);
}`,
    oGradRho: `${K}  let e = ${idx}; if (e >= E) { return; }
  OD[O_GRADRHO + e] = (OD[O_RHOML + MI[COE + 2 * e + 1]] - OD[O_RHOML + MI[COE + 2 * e]]) / MF[F_DC + e];
}`,
    oEdgeThicknessRaw: `${K}  let n = ${idx}; if (n >= L * E) { return; }
  let k = n / E; let e = n % E;
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  if (k == 0) { OD[O_HEDGE + e] = 0.5 * (IN[hOff(0) + a] + IN[hOff(0) + b]); }
  else { OD[O_HEDGE + k * E + e] = min(IN[hOff(k) + a], IN[hOff(k) + b]); }
}`,
    oEdgeThicknessSill: `${K}  let e = ${idx}; if (e >= E) { return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  let sill = max(EPSO, min(OD[O_BATH + a], OD[O_BATH + b]) + 0.5 * (OD[O_ETA + a] + OD[O_ETA + b]));
  var sum = 0.0;
  for (var k = 0; k < L; k++) { sum += OD[O_HEDGE + k * E + e]; }
  if (sum > sill) { let f = sill / sum; for (var k = 0; k < L; k++) { OD[O_HEDGE + k * E + e] *= f; } }
}`,
    oFlux: `${K}  let n = ${idx}; if (n >= L * E) { return; }
  let k = n / E; let e = n % E;
  var he = OD[O_HEDGE + n];
  if (k == 0) {
    let donor = MI[COE + 2 * e + select(1, 0, IN[uOff(0) + e] > 0.0)];
    he = min(he, max(0.0, IN[hOff(0) + donor]));
  }
  OD[O_FLUX + n] = select(0.0, he * IN[uOff(k) + e], OD[O_EMASK + e] > 0.5);
}`,
    oCellTendency: `${K}  let n = ${idx}; if (n >= L * C) { return; }
  let k = n / C; let i = n % C;
  let hh = max(EPSO, IN[hOff(k) + i]); let Ti = IN[qOff(k) + i] / hh; let Si = IN[wOff(k) + i] / hh;
  var divH = 0.0; var divQ = 0.0; var divW = 0.0; var lapQ = 0.0; var lapW = 0.0;
  for (var m = 0; m < MI[NEC + i]; m++) {
    let e = MI[EOC + MAXE * i + m]; let j = MI[COC + MAXE * i + m];
    let f = f32(MI[ESC + MAXE * i + m]) * OD[O_FLUX + k * E + e] * MF[F_DV + e];
    divH += f;
    let hhj = max(EPSO, IN[hOff(k) + j]); let Tj = IN[qOff(k) + j] / hhj; let Sj = IN[wOff(k) + j] / hhj;
    divQ += f * select(Tj, Ti, f > 0.0); divW += f * select(Sj, Si, f > 0.0);
    if (k == 0 && OD[O_EMASK + e] > 0.5) { lapQ += MF[F_DV + e] * (Tj - Ti) / MF[F_DC + e]; lapW += MF[F_DV + e] * (Sj - Si) / MF[F_DC + e]; }
  }
  let area = MF[F_AREA + i];
  var dh = -divH / area; var dQ = -divQ / area; var dW = -divW / area;
  if (k == 0 && DIFFUSION > 0.0) { dQ += DIFFUSION * lapQ / area; dW += DIFFUSION * lapW / area; }
  if (OD[O_CMASK + i] < 0.5) { dh = 0.0; dQ = 0.0; dW = 0.0; }
  OUT[hOff(k) + i] = dh; OUT[qOff(k) + i] = dQ; OUT[wOff(k) + i] = dW;
}`,
    oVertexVort: `${K}  let n = ${idx}; if (n >= L * V) { return; }
  let k = n / V; let v = n % V;
  var zeta = 0.0;
  for (var m = 0; m < 3; m++) { let e = MI[EOV + 3 * v + m]; zeta += f32(MI[ESV + 3 * v + m]) * IN[uOff(k) + e] * MF[F_DC + e]; }
  OD[O_AVORT + n] = zeta / MF[F_ATRI + v] + MF[F_FV + v];
}`,
    oEdgePV: `${K}  let n = ${idx}; if (n >= L * E) { return; }
  let k = n / E; let e = n % E;
  let v1 = MI[VOE + 2 * e]; let v2 = MI[VOE + 2 * e + 1];
  OD[O_QE + n] = 0.5 * (OD[O_AVORT + k * V + v1] + OD[O_AVORT + k * V + v2]) / max(OD[O_HEDGE + n], PVFLOOR);
}`,
    oKineticPhi: `${K}  let n = ${idx}; if (n >= L * C) { return; }
  let k = n / C; let i = n % C;
  var kinetic = 0.0;
  for (var m = 0; m < MI[NEC + i]; m++) { let e = MI[EOC + MAXE * i + m]; let u = IN[uOff(k) + e]; kinetic += 0.25 * MF[F_DC + e] * MF[F_DV + e] * u * u; }
  kinetic = kinetic / MF[F_AREA + i];
  var phi = kinetic + OGRAV * OD[O_ETA + i];
  if (k > 0) {
    var p = (OD[O_RHOML + i] - RHO[k]) * IN[hOff(0) + i];
    for (var j = 1; j < k; j++) { p += (RHO[j] - RHO[k]) * IN[hOff(j) + i]; }
    phi += OGRAV * p / RHO0;
  }
  OD[O_PHI + n] = phi;
}`,
    oDivCurl: `${K}  let n = ${idx};
  let fromLap = P[1] > 0.5;
  if (n < L * C) {
    let k = n / C; let i = n % C;
    var sum = 0.0;
    for (var m = 0; m < MI[NEC + i]; m++) { let e = MI[EOC + MAXE * i + m]; let u = select(IN[uOff(k) + e], OD[O_LAPA + k * E + e], fromLap); sum += f32(MI[ESC + MAXE * i + m]) * u * MF[F_DV + e]; }
    OD[O_DIVS + n] = sum / MF[F_AREA + i];
  }
  if (n < L * V) {
    let k = n / V; let v = n % V;
    var sum = 0.0;
    for (var m = 0; m < 3; m++) { let e = MI[EOV + 3 * v + m]; let u = select(IN[uOff(k) + e], OD[O_LAPA + k * E + e], fromLap); sum += f32(MI[ESV + 3 * v + m]) * u * MF[F_DC + e]; }
    OD[O_CURLS + n] = sum / MF[F_ATRI + v];
  }
}`,
    oLapVelocity: `${K}  let n = ${idx}; if (n >= L * E) { return; }
  let k = n / E; let e = n % E;
  let lap = (OD[O_DIVS + k * C + MI[COE + 2 * e + 1]] - OD[O_DIVS + k * C + MI[COE + 2 * e]]) / MF[F_DC + e]
    - (OD[O_CURLS + k * V + MI[VOE + 2 * e + 1]] - OD[O_CURLS + k * V + MI[VOE + 2 * e]]) / MF[F_DV + e];
  if (P[1] > 0.5) { OD[O_LAPB + n] = lap; } else { OD[O_LAPA + n] = lap; }
}`,
    oMomentum: `${K}  let n = ${idx}; if (n >= L * E) { return; }
  let k = n / E; let e = n % E;
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  let qHere = 0.5 * OD[O_QE + n];
  var pv = 0.0;
  for (var s = 0; s < MI[NEE + e]; s++) { let slot = MAXEE * e + s; let other = MI[EOE + slot]; pv += MF[F_PVW + slot] * OD[O_FLUX + k * E + other] * (qHere + 0.5 * OD[O_QE + k * E + other]); }
  let dc = MF[F_DC + e];
  let gradPhi = (OD[O_PHI + k * C + b] - OD[O_PHI + k * C + a]) / dc;
  var du = pv / dc - gradPhi;
  if (k == 0) { du -= OGRAV / RHO0 * 0.5 * OD[O_HEDGE + e] * OD[O_GRADRHO + e]; }
  let he = max(OD[O_HEDGE + n], MINTHICK);
  var force = 0.0;
  if (k == 0) { force += OD[O_STRESS + e] / RHO0; }
  if (k > 0) { force += RINT * (IN[uOff(k - 1) + e] - IN[uOff(k) + e]); }
  if (k < L - 1) { force -= RINT * (IN[uOff(k) + e] - IN[uOff(k + 1) + e]); }
  var bottom = k == L - 1;
  if (!bottom) { bottom = true; for (var j = k + 1; j < L; j++) { if (OD[O_HEDGE + j * E + e] >= THINO) { bottom = false; break; } } }
  if (bottom) { force -= RBOT * abs(IN[uOff(k) + e]) * IN[uOff(k) + e]; }
  du += force / he;
  if (NU4O > 0.0) { du -= NU4O * OD[O_LAPB + n]; }
  if (k > 0 && OD[O_HEDGE + n] < THINO) { du = (IN[uOff(k - 1) + e] - IN[uOff(k) + e]) * P[2]; }
  if (OD[O_EMASK + e] < 0.5) { du = 0.0; }
  OUT[uOff(k) + e] = du;
}`,
    oAdvance: `${K}  let n = ${idx}; if (n >= OSTOTAL) { return; }
  OUT[n] = IN[n] + P[0] * OD[n];
}`,
    oCombine: `${K}  let n = ${idx}; if (n >= OSTOTAL) { return; }
  IN[n] += P[0] * (OUT[n] + 2.0 * OD[n] + 2.0 * PH[n] + S[n]);
}`,
    oBarotropicSetup: `${K}  let n = ${idx};
  if (n < C) { OD[B_BCUR + n] = OD[O_ETA + n]; }
  if (n < E) {
    let e = n;
    var sumH = 0.0; var transport = 0.0; var forcing = 0.0;
    for (var k = 0; k < L; k++) { let he = OD[O_HEDGE + k * E + e]; sumH += he; transport += he * IN[uOff(k) + e]; forcing += he * (OUT[uOff(k) + e] + OGRAV * OD[O_GRADETA + e]); }
    let onOcean = OD[O_EMASK + e] > 0.5;
    OD[B_BCUR + C + e] = select(0.0, transport, onOcean);
    OD[O_SLOW + e] = select(0.0, forcing, onOcean);
    OD[O_DEPTHEDGE + e] = sumH;
  }
}`,
    oBarotropicCoriolis: `${K}  let e = ${idx}; if (e >= E) { return; }
  if (OD[O_EMASK + e] < 0.5) { return; }
  var sum = 0.0;
  for (var s = 0; s < MI[NEE + e]; s++) { let slot = MAXEE * e + s; let other = MI[EOE + slot]; sum += MF[F_PVW + slot] * OD[B_BCUR + C + other] * 0.5 * (OD[O_FEDGE + e] + OD[O_FEDGE + other]); }
  OD[O_SLOW + e] -= sum / MF[F_DC + e];
}`,
    /*
     * Four fixed-block variants (rather than one kernel parametrized at
     * runtime) so the whole M-substep barotropic loop can be recorded
     * into one pass with parameters set once, not once per dispatch:
     * within a single un-submitted command encoder, queue.writeBuffer
     * calls to the params buffer all land before any dispatch actually
     * runs on the device, so a dispatch sequence that needs a different
     * parameter per step must either resubmit between them or, as here,
     * bake the varying offsets into the shader at build time.
     */
    ...Object.fromEntries([1, 2, 3, 4].map((stage) => {
      const inBase = stage === 1 ? 'B_BCUR' : 'B_BTRIAL';
      const outBase = `B_BK${stage}`;
      return [`oBarotropicTendency${stage}`, `${K}  let n = ${idx};
  if (n < C) {
    let i = n;
    var sum = 0.0;
    for (var m = 0; m < MI[NEC + i]; m++) { let e = MI[EOC + MAXE * i + m]; sum += f32(MI[ESC + MAXE * i + m]) * OD[${inBase} + C + e] * MF[F_DV + e]; }
    OD[${outBase} + i] = select(0.0, -(sum / MF[F_AREA + i]), OD[O_CMASK + i] > 0.5);
  }
  if (n < E) {
    let e = n;
    if (OD[O_EMASK + e] < 0.5) { OD[${outBase} + C + e] = 0.0; } else {
      let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
      let gradEtaB = (OD[${inBase} + b] - OD[${inBase} + a]) / MF[F_DC + e];
      var cor = 0.0;
      for (var s = 0; s < MI[NEE + e]; s++) { let slot = MAXEE * e + s; let other = MI[EOE + slot]; cor += MF[F_PVW + slot] * OD[${inBase} + C + other] * 0.5 * (OD[O_FEDGE + e] + OD[O_FEDGE + other]); }
      OD[${outBase} + C + e] = -OGRAV * OD[O_DEPTHEDGE + e] * gradEtaB + cor / MF[F_DC + e] + OD[O_SLOW + e];
    }
  }
}`];
    })),
    ...Object.fromEntries([1, 2, 3].map((stage) => {
      const stageBase = `B_BK${stage}`;
      const factorSlot = stage === 3 ? 'P[1]' : 'P[0]';
      return [`oBarotropicCombine${stage}`, `${K}  let n = ${idx}; if (n >= C + E) { return; }
  OD[B_BTRIAL + n] = OD[B_BCUR + n] + ${factorSlot} * OD[${stageBase} + n];
}`];
    })),
    oBarotropicFinalCombine: `${K}  let n = ${idx}; if (n >= C + E) { return; }
  OD[B_BCUR + n] += P[2] * (OD[B_BK1 + n] + 2.0 * OD[B_BK2 + n] + 2.0 * OD[B_BK3 + n] + OD[B_BK4 + n]);
}`,
    oBarotropicAccumulate: `${K}  let n = ${idx}; if (n >= C + E) { return; }
  OD[B_BAVG + n] += OD[B_BCUR + n] * P[3];
}`,
    oRescale: `${K}  let i = ${idx}; if (i >= C) { return; }
  if (OD[O_CMASK + i] < 0.5) { return; }
  var sum = 0.0;
  for (var k = 0; k < L; k++) {
    var hv = IN[hOff(k) + i];
    if (hv < EPSO) {
      var t = LABEL_T[k]; var s = REF_S;
      if (hv > 1e-9) { t = IN[qOff(k) + i] / hv; s = IN[wOff(k) + i] / hv; }
      hv = EPSO;
      IN[hOff(k) + i] = EPSO; IN[qOff(k) + i] = EPSO * t; IN[wOff(k) + i] = EPSO * s;
    }
    sum += hv;
  }
  let scale = (OD[O_BATH + i] + OD[B_BAVG + i]) / sum;
  for (var k = 0; k < L; k++) { IN[hOff(k) + i] *= scale; IN[qOff(k) + i] *= scale; IN[wOff(k) + i] *= scale; }
  OD[O_ETA + i] = OD[B_BAVG + i];
}`,
    oVelocityShiftClamp: `${K}  let e = ${idx}; if (e >= E) { return; }
  if (OD[O_EMASK + e] < 0.5) { return; }
  var sumH = 0.0; var transport = 0.0;
  for (var k = 0; k < L; k++) { let he = OD[O_HEDGE + k * E + e]; sumH += he; transport += he * IN[uOff(k) + e]; }
  let shift = (OD[B_BAVG + C + e] - transport) / max(sumH, EPSO);
  for (var k = 0; k < L; k++) { let n = uOff(k) + e; IN[n] = clamp(IN[n] + shift, -SPEEDLIM, SPEEDLIM); }
  for (var k = 1; k < L; k++) { if (OD[O_HEDGE + k * E + e] < THINO) { IN[uOff(k) + e] = IN[uOff(k - 1) + e]; } }
}`,
    oReadSurface: `${K}  let i = ${idx}; if (i >= C) { return; }
  let iced = OD[O_SURFICE + i] > 0.0;
  OD[O_ICED + i] = select(0.0, 1.0, iced);
  let t0 = select(OD[O_SURFT + i], FREEZE, iced);
  OD[O_SURFACEIN + i] = t0;
  IN[qOff(0) + i] = IN[hOff(0) + i] * t0;
}`,
    oStressFromAtmosphere: `${K}  let e = ${idx}; if (e >= E) { return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  if (OD[O_EMASK + e] < 0.5 || OD[O_ICED + a] > 0.5 || OD[O_ICED + b] > 0.5) { OD[O_STRESS + e] = 0.0; return; }
  let bottom = (K - 1) * C;
  let rhoA = S[S_PI + a] * LV[L_SM + K - 1] / (RGAS * S[S_TH + bottom + a] * D[D_EXM + bottom + a]);
  let rhoB = S[S_PI + b] * LV[L_SM + K - 1] / (RGAS * S[S_TH + bottom + b] * D[D_EXM + bottom + b]);
  let fa = CDO * rhoA * max(D[D_WIND + a], GUSTO); let fb = CDO * rhoB * max(D[D_WIND + b], GUSTO);
  OD[O_STRESS + e] = 0.5 * (fa + fb) * S[S_U + (K - 1) * E + e];
}`,
    oMixedLayer: `${K}  let i = ${idx}; if (i >= C) { return; }
  if (OD[O_CMASK + i] < 0.5) { return; }
  var wv = vec3<f32>(0.0, 0.0, 0.0);
  for (var m = 0; m < MI[NEC + i]; m++) {
    let e = MI[EOC + MAXE * i + m];
    let s = 0.5 * MF[F_DC + e] * MF[F_DV + e] * OD[O_STRESS + e];
    wv += s * vec3<f32>(MF[F_NEDGE + 3 * e], MF[F_NEDGE + 3 * e + 1], MF[F_NEDGE + 3 * e + 2]);
  }
  wv = wv / MF[F_AREA + i];
  let tau = length(wv);
  let ustar3 = pow(tau / RHO0, 1.5); let stir = STIRRING * exp(-IN[hOff(0) + i] / STIRDEPTH);
  var rm = eos(IN[qOff(0) + i] / IN[hOff(0) + i], IN[wOff(0) + i] / IN[hOff(0) + i]);
  for (var k = 1; k < L; k++) {
    if (IN[hOff(0) + i] >= MAXMIXED) { break; }
    if (IN[hOff(k) + i] <= EPSO || RHO[k] > rm) { continue; }
    moveLayer(i, k, 0, min(IN[hOff(k) + i] - EPSO, MAXMIXED - IN[hOff(0) + i]));
    rm = eos(IN[qOff(0) + i] / IN[hOff(0) + i], IN[wOff(0) + i] / IN[hOff(0) + i]);
  }
  var below = -1;
  for (var k = 1; k < L; k++) { if (IN[hOff(k) + i] > THINO) { below = k; break; } }
  if (below > 0 && IN[hOff(0) + i] < MAXMIXED) {
    let db = max(1e-4, OGRAV * (RHO[below] - rm) / RHO0);
    let entrain = min(2.0 * stir * ustar3 / (IN[hOff(0) + i] * db) * P[6], IN[hOff(below) + i] - EPSO);
    if (entrain > 0.0) { moveLayer(i, below, 0, entrain); rm = eos(IN[qOff(0) + i] / IN[hOff(0) + i], IN[wOff(0) + i] / IN[hOff(0) + i]); }
  }
  below = -1;
  for (var k = 1; k < L; k++) { if (IN[hOff(k) + i] > THINO) { below = k; break; } }
  var excess = max(0.0, IN[hOff(0) + i] - MAXMIXED);
  if (below > 0 && rm >= RHO[below] - DENSTOL) { excess = max(excess, IN[hOff(0) + i] - SHALLOWMIXED); }
  detrain(i, excess, rm);
  let buoyancy = OGRAV * THERMAL_EXP * (OD[O_PREVT0 + i] - OD[O_SURFACEIN + i]) * IN[hOff(0) + i] / P[6];
  if (buoyancy < -1e-9) {
    let monin = max(SHALLOWMIXED, 2.0 * stir * ustar3 / -buoyancy);
    if (IN[hOff(0) + i] > monin) { detrain(i, (IN[hOff(0) + i] - monin) * min(1.0, P[6] / DETRAINT), rm); }
  }
  if (IN[hOff(0) + i] < MINTHICK) {
    for (var k = 1; k < L; k++) {
      if (IN[hOff(0) + i] >= MINTHICK) { break; }
      let available = IN[hOff(k) + i] - EPSO;
      if (available > 0.0) { moveLayer(i, k, 0, min(available, MINTHICK - IN[hOff(0) + i])); }
    }
  }
}`,
    oSalt: `${K}  let i = ${idx}; if (i >= C) { return; }
  if (OD[O_CMASK + i] < 0.5) { return; }
  let s = IN[wOff(0) + i] / IN[hOff(0) + i];
  var w = IN[wOff(0) + i] + s * OD[O_FRESH + i] / 1000.0;
  let grown = (OD[O_SURFICE + i] - OD[O_PREVICE + i]) * ICEDENS / 1000.0;
  w += (s - ICESAL) * grown;
  w = max(0.0, w);
  IN[wOff(0) + i] = w;
  OD[O_FRESH + i] = 0.0;
  OD[O_PREVICE + i] = OD[O_SURFICE + i];
}`,
    oWriteSurface: `${K}  let i = ${idx}; if (i >= C) { return; }
  if (OD[O_CMASK + i] < 0.5) { PH[PH_OFLUX + i] = 0.0; return; }
  PH[PH_CAP + i] = RHOCP * max(IN[hOff(0) + i], 1.0);
  OD[O_S0 + i] = IN[wOff(0) + i] / IN[hOff(0) + i];
  if (OD[O_ICED + i] > 0.5) {
    PH[PH_OFLUX + i] = RHOCP * (IN[qOff(0) + i] - IN[hOff(0) + i] * FREEZE) / P[6];
    IN[qOff(0) + i] = IN[hOff(0) + i] * FREEZE;
    OD[O_T0 + i] = FREEZE;
  } else {
    OD[O_T0 + i] = IN[qOff(0) + i] / IN[hOff(0) + i];
    S[S_TS + i] = OD[O_T0 + i];
    PH[PH_OFLUX + i] = 0.0;
  }
  OD[O_PREVT0 + i] = OD[O_T0 + i];
}`,
    oAccumulateFresh: `${K}  let i = ${idx}; if (i >= C) { return; }
  let rain = PH[PH_RAIN + i]; let seen = OD[O_RAINSEEN + i];
  let delta = select(rain, rain - seen, rain >= seen);
  if (OD[O_CMASK + i] > 0.5) { OD[O_FRESH + i] += PH[PH_EVAP + i] * P[6] - delta; }
  OD[O_RAINSEEN + i] = rain;
}`,
  };
}

export function createLayeredOcean(core, options = {}) {
  const o = { ...OCEAN_DEFAULTS, ...options };
  const { device, buffers, mesh, meshSpacing } = core;
  const C = core.C, E = core.E, V = core.V;
  const L = o.densities.length + 1;
  const rho = [o.density, ...o.densities];
  const labelT = rho.map((r) => Math.max(FREEZING_POINT, o.referenceT - (r / o.density - 1) / o.thermalExpansion));
  const nu4 = o.closureHours > 0 ? Math.pow(meshSpacing / Math.PI, 4) / (o.closureHours * 3600) : 0;
  const diffusion = o.diffusivity * mesh.radius * mesh.radius / (o.density * o.specificHeat);
  let minSpacing = Infinity;
  for (let e = 0; e < E; e++) minSpacing = Math.min(minSpacing, mesh.dcEdge[e]);
  const geography = o.geography || null;
  const edgeOcean = geography ? geography.edgeOcean : new Uint8Array(E).fill(1);
  const cellOcean = geography ? Uint8Array.from(geography.land, (l) => (l ? 0 : 1)) : new Uint8Array(C).fill(1);
  const D = new Float64Array(C);
  const smoothed = o.bathymetry ? null : bathymetryFrom(mesh, geography, { minimumDepth: o.minimumDepth, flatDepth: o.flatDepth });
  for (let i = 0; i < C; i++) D[i] = !cellOcean[i] ? 0 : o.bathymetry ? o.bathymetry[i] : smoothed[i];
  let deepest = 0;
  for (let i = 0; i < C; i++) deepest = Math.max(deepest, D[i]);
  const substepLimit = 0.35 * minSpacing / Math.sqrt(o.gravity * Math.max(deepest, 1));
  const salinityProfile = o.salinityProfile || defaultSalinityProfile;

  const OS = seq([['OH', L * C], ['OU', L * E], ['OQ', L * C], ['OW', L * C]]);
  const OD = seq([
    ['FLUX', L * E], ['HEDGE', L * E], ['AVORT', L * V], ['QE', L * E], ['PHI', L * C],
    ['LAPA', L * E], ['LAPB', L * E], ['DIVS', L * C], ['CURLS', L * V],
    ['RHOML', C], ['GRADRHO', E], ['GRADETA', E], ['SLOW', E], ['DEPTHEDGE', E],
    ['ETA', C], ['FRESH', C], ['PREVT0', C], ['PREVICE', C], ['ICED', C], ['STRESS', E],
    ['SURFT', C], ['SURFICE', C], ['SURFACEIN', C], ['T0', C], ['S0', C], ['RAINSEEN', C],
    ['EMASK', E], ['CMASK', C], ['BATH', C], ['FEDGE', E],
  ]);
  // The barotropic RK4's fixed-offset blocks index into the same OD storage
  // array (see oceanKernels' barotropic tendency/combine kernels), so their
  // offsets must continue on from OD's own sections rather than starting a
  // second, colliding index space at 0.
  const B = seq([['BCUR', C + E], ['BK1', C + E], ['BK2', C + E], ['BK3', C + E], ['BK4', C + E], ['BTRIAL', C + E], ['BAVG', C + E]]);
  const bTotal = B.total;
  for (const k of Object.keys(B)) if (k !== 'total') B[k] += OD.total;
  B.total = bTotal;
  const ODTOTAL = OD.total + bTotal;

  const kernels = oceanKernels({ ...o, L, C, E, V, OS, OD, B, rho, labelT, nu4, diffusion });
  const ob = {
    S: emptyBuffer(device, 4 * OS.total), T: emptyBuffer(device, 4 * OS.total),
    K1: emptyBuffer(device, 4 * OS.total), K2: emptyBuffer(device, 4 * OS.total), K3: emptyBuffer(device, 4 * OS.total), K4: emptyBuffer(device, 4 * OS.total),
    OD: emptyBuffer(device, 4 * ODTOTAL),
  };
  for (const [name, b] of Object.entries(ob)) b.label = 'layeredOcean' + name;
  const bindLayout = device.createBindGroupLayout({ entries: Array.from({ length: 10 }, (_, binding) => ({ binding, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } })) });
  const pipelineLayout = device.createPipelineLayout({ bindGroupLayouts: [bindLayout] });
  const head = core.preludeConstants + kernels.head;
  const pipelines = {};
  for (const [name, body] of Object.entries(kernels)) {
    if (name === 'head') continue;
    const module = device.createShaderModule({ code: head + body, label: name });
    pipelines[name] = device.createComputePipeline({ label: name, layout: pipelineLayout, compute: { module, entryPoint: 'main' } });
  }
  const groups = new Map();
  function group(IN, OUT, ODbuf = ob.OD, PH = buffers.PH, Sbuf = buffers.S) {
    const key = [IN.label, OUT.label, ODbuf.label, PH.label, Sbuf.label].join('|');
    let g = groups.get(key);
    if (!g) {
      g = device.createBindGroup({ layout: bindLayout, entries: [buffers.MI, buffers.MF, buffers.LV, IN, OUT, ODbuf, buffers.P, PH, Sbuf, buffers.D].map((buffer, binding) => ({ binding, resource: { buffer } })) });
      groups.set(key, g);
    }
    return g;
  }
  function dispatch(pass, name, bindGroup, count) {
    pass.setPipeline(pipelines[name]);
    pass.setBindGroup(0, bindGroup);
    const n = Math.ceil(count / WORKGROUP);
    pass.dispatchWorkgroups(Math.min(n, 65535), Math.ceil(n / 65535));
  }
  const params = new Float32Array(8);
  const setParams = (values) => { params.fill(0); params.set(values); device.queue.writeBuffer(buffers.P, 0, params); };

  function tendency(IN, OUT) {
    const g = group(IN, OUT);
    setParams([0, 0]);
    let encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
    dispatch(pass, 'oSurfaceDensity', g, C);
    dispatch(pass, 'oGradRho', g, E);
    dispatch(pass, 'oEdgeThicknessRaw', g, L * E);
    dispatch(pass, 'oEdgeThicknessSill', g, E);
    dispatch(pass, 'oFlux', g, L * E);
    dispatch(pass, 'oCellTendency', g, L * C);
    dispatch(pass, 'oVertexVort', g, L * V);
    dispatch(pass, 'oEdgePV', g, L * E);
    dispatch(pass, 'oKineticPhi', g, L * C);
    dispatch(pass, 'oDivCurl', g, Math.max(L * C, L * V));
    dispatch(pass, 'oLapVelocity', g, L * E);
    pass.end();
    device.queue.submit([encoder.finish()]);
    setParams([0, 1, relaxRate]);
    encoder = device.createCommandEncoder(); pass = encoder.beginComputePass();
    dispatch(pass, 'oDivCurl', g, Math.max(L * C, L * V));
    dispatch(pass, 'oLapVelocity', g, L * E);
    dispatch(pass, 'oMomentum', g, L * E);
    pass.end();
    device.queue.submit([encoder.finish()]);
  }
  function advanceState(next, stage, factor) {
    setParams([factor]);
    const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
    dispatch(pass, 'oAdvance', group(ob.S, next, stage), OS.total);
    pass.end();
    device.queue.submit([encoder.finish()]);
  }

  let relaxRate = 1 / 3600;
  function barotropicStep(dt, k1) {
    const g1 = group(ob.S, k1);
    {
      const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
      dispatch(pass, 'oBarotropicSetup', g1, Math.max(C, E));
      dispatch(pass, 'oBarotropicCoriolis', g1, E);
      pass.end();
      device.queue.submit([encoder.finish()]);
    }
    device.queue.writeBuffer(ob.OD, 4 * B.BAVG, new Float32Array(C + E));
    const M = Math.max(1, Math.ceil(dt / substepLimit)), dtb = dt / M;
    const gAny = group(ob.S, ob.K1);
    // P0=dtb/2 (stage1,2 trial factor), P1=dtb (stage3 trial factor), P2=dtb/6 (final combine), P3=1/M (average accumulation).
    setParams([dtb / 2, dtb, dtb / 6, 1 / M]);
    const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
    for (let m = 0; m < M; m++) {
      dispatch(pass, 'oBarotropicTendency1', gAny, Math.max(C, E));
      dispatch(pass, 'oBarotropicCombine1', gAny, C + E);
      dispatch(pass, 'oBarotropicTendency2', gAny, Math.max(C, E));
      dispatch(pass, 'oBarotropicCombine2', gAny, C + E);
      dispatch(pass, 'oBarotropicTendency3', gAny, Math.max(C, E));
      dispatch(pass, 'oBarotropicCombine3', gAny, C + E);
      dispatch(pass, 'oBarotropicTendency4', gAny, Math.max(C, E));
      dispatch(pass, 'oBarotropicFinalCombine', gAny, C + E);
      dispatch(pass, 'oBarotropicAccumulate', gAny, C + E);
    }
    pass.end();
    device.queue.submit([encoder.finish()]);
  }

  function step(dt) {
    relaxRate = Math.min(1 / 3600, 1 / dt);
    {
      const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
      dispatch(pass, 'oFreeSurface', group(ob.S, ob.T), C);
      dispatch(pass, 'oGradEta', group(ob.S, ob.T), E);
      pass.end();
      device.queue.submit([encoder.finish()]);
    }
    tendency(ob.S, ob.K1);
    barotropicStep(dt, ob.K1);
    advanceState(ob.T, ob.K1, dt / 2);
    tendency(ob.T, ob.K2); advanceState(ob.T, ob.K2, dt / 2);
    tendency(ob.T, ob.K3); advanceState(ob.T, ob.K3, dt);
    tendency(ob.T, ob.K4);
    setParams([dt / 6]);
    {
      const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
      dispatch(pass, 'oCombine', group(ob.S, ob.K1, ob.K2, ob.K3, ob.K4), OS.total);
      pass.end();
      device.queue.submit([encoder.finish()]);
    }
    {
      const g = group(ob.S, ob.T);
      const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
      dispatch(pass, 'oRescale', g, C);
      dispatch(pass, 'oEdgeThicknessRaw', g, L * E);
      dispatch(pass, 'oEdgeThicknessSill', g, E);
      dispatch(pass, 'oVelocityShiftClamp', g, E);
      pass.end();
      device.queue.submit([encoder.finish()]);
    }
  }

  function readSurface(surfaceT, ice) {
    device.queue.writeBuffer(ob.OD, 4 * OD.SURFT, Float32Array.from(surfaceT));
    device.queue.writeBuffer(ob.OD, 4 * OD.SURFICE, Float32Array.from(ice));
    const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
    dispatch(pass, 'oReadSurface', group(ob.S, ob.T), C);
    pass.end();
    device.queue.submit([encoder.finish()]);
  }
  function readSurfaceFromAtmosphere() {
    const encoder = device.createCommandEncoder();
    encoder.copyBufferToBuffer(buffers.S, 4 * core.layout.S.TS, ob.OD, 4 * OD.SURFT, 4 * C);
    encoder.copyBufferToBuffer(buffers.S, 4 * core.layout.S.ICE, ob.OD, 4 * OD.SURFICE, 4 * C);
    device.queue.submit([encoder.finish()]);
    const pass2 = device.createCommandEncoder(), p2 = pass2.beginComputePass();
    dispatch(p2, 'oReadSurface', group(ob.S, ob.T), C);
    p2.end();
    device.queue.submit([pass2.finish()]);
  }
  function stressFromAtmosphere() {
    const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
    dispatch(pass, 'oStressFromAtmosphere', group(ob.S, ob.T), E);
    pass.end();
    device.queue.submit([encoder.finish()]);
  }
  function setStress(total, ice) {
    const masked = new Float32Array(E);
    for (let e = 0; e < E; e++) masked[e] = !edgeOcean[e] || ice[mesh.cellsOnEdge[2 * e]] > 0 || ice[mesh.cellsOnEdge[2 * e + 1]] > 0 ? 0 : total[e];
    device.queue.writeBuffer(ob.OD, 4 * OD.STRESS, masked);
  }
  function mixedLayer(dt) {
    setParams([0, 0, relaxRate, 0, 0, 0, dt]);
    const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
    dispatch(pass, 'oMixedLayer', group(ob.S, ob.T), C);
    pass.end();
    device.queue.submit([encoder.finish()]);
  }
  function salt(dt) {
    setParams([0, 0, 0, 0, 0, 0, dt]);
    const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
    dispatch(pass, 'oSalt', group(ob.S, ob.T), C);
    pass.end();
    device.queue.submit([encoder.finish()]);
  }
  function writeSurface(dt) {
    setParams([0, 0, 0, 0, 0, 0, dt]);
    const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
    dispatch(pass, 'oWriteSurface', group(ob.S, ob.T), C);
    pass.end();
    device.queue.submit([encoder.finish()]);
  }
  function accumulateFreshwater(dt) {
    setParams([0, 0, 0, 0, 0, 0, dt]);
    const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
    dispatch(pass, 'oAccumulateFresh', group(ob.S, ob.T), C);
    pass.end();
    device.queue.submit([encoder.finish()]);
  }

  /*
   * A single testable/coupled ocean step at double-precision-matched
   * physics but single-precision arithmetic: readSurface, setStress,
   * step, mixedLayer, salt, writeSurface, unconditionally (the caller
   * decides when to call it; model.gpu.js gates on its own everySteps
   * counter so it only pays for the atmosphere-state read and the step
   * on the steps that need it).
   */
  async function advance(surfaceT, ice, totalStress, dt) {
    readSurface(surfaceT, ice);
    setStress(typeof totalStress === 'function' ? totalStress() : totalStress, ice);
    step(dt);
    mixedLayer(dt);
    salt(dt);
    writeSurface(dt);
    await device.queue.onSubmittedWorkDone();
    return true;
  }
  async function advanceCoupled(dt) {
    readSurfaceFromAtmosphere();
    stressFromAtmosphere();
    step(dt);
    mixedLayer(dt);
    salt(dt);
    writeSurface(dt);
    await device.queue.onSubmittedWorkDone();
  }

  /*
   * Initialization, loading, serialization and diagnostics run once per
   * model build or at the (infrequent) cadence a caller asks for, so
   * they port the CPU functions over plain double-precision arrays and
   * push the result to the GPU, rather than duplicating that logic in
   * WGSL.
   */
  const eos = (t, s) => o.density * (1 - o.thermalExpansion * (t - o.referenceT) + o.halineContraction * (s - o.referenceS));
  /*
   * The free surface that levels the pressure at `referenceDepth` in
   * every column deep enough to reach it, so the deep ocean starts
   * without a barotropic pressure gradient; shallower columns take the
   * free surface of the deep water around them, found by relaxation.
   * Ported over plain double-precision arrays exactly as the CPU module's
   * initialize()/stericSurface() do, since this runs once per model
   * build rather than every ocean step.
   */
  function stericSurface(h, Q, W, eta, referenceDepth = 3500) {
    const at = (k, i) => k * C + i;
    const rhoMl = new Float64Array(C);
    for (let i = 0; i < C; i++) { const h0 = Math.max(EPS, h[i]); rhoMl[i] = eos(Q[i] / h0, W[i] / h0); }
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
      eta[i] = -anomaly / o.density;
    }
    const next = new Float64Array(C);
    for (let sweep = 0; sweep < 2000; sweep++) {
      let moved = 0;
      for (let i = 0; i < C; i++) {
        if (!cellOcean[i] || deep[i]) { next[i] = eta[i]; continue; }
        let sum = 0, weight = 0;
        for (let m = 0; m < mesh.nEdgesOnCell[i]; m++) {
          const j = mesh.cellsOnCell[mesh.maxEdges * i + m];
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
    for (let i = 0; i < C; i++) if (cellOcean[i]) { area += mesh.areaCell[i]; mean += mesh.areaCell[i] * eta[i]; }
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
  function initializeArrays(surfaceT, ice) {
    const h = new Float64Array(L * C), u = new Float64Array(L * E), Q = new Float64Array(L * C), W = new Float64Array(L * C);
    const eta = new Float64Array(C), T0 = new Float64Array(C), previousT0 = new Float64Array(C);
    const at = (k, i) => k * C + i;
    for (let i = 0; i < C; i++) {
      const iced = ice[i] > 0;
      const lat = mesh.latCell[i];
      const s0 = salinityProfile(lat);
      T0[i] = iced ? FREEZING_POINT : surfaceT[i];
      if (!cellOcean[i]) continue;
      const rm = eos(T0[i], s0);
      const stretchReal = 1 - o.thermoclineTilt + 2 * o.thermoclineTilt * Math.cos(lat) ** 2;
      let cumulative = Math.min(o.mixedDepth, D[i]);
      h[i] = cumulative; Q[i] = cumulative * T0[i]; W[i] = cumulative * s0;
      for (let k = 1; k < L; k++) {
        let hk;
        if (k === L - 1) hk = Math.max(0, D[i] - cumulative);
        else {
          const taper = Math.max(0, Math.min(1, (rho[k] - rm) / (rho[k + 1] - rho[k])));
          hk = Math.max(0, Math.min(D[i], o.bottoms[k - 1] * stretchReal * taper) - cumulative);
        }
        hk = Math.max(EPS, hk);
        cumulative += hk;
        h[at(k, i)] = hk; Q[at(k, i)] = hk * labelT[k]; W[at(k, i)] = hk * o.referenceS;
      }
      const scale = D[i] / cumulative;
      for (let k = 0; k < L; k++) { h[at(k, i)] *= scale; Q[at(k, i)] *= scale; W[at(k, i)] *= scale; }
      previousT0[i] = T0[i];
    }
    stericSurface(h, Q, W, eta);
    return { h, u, Q, W, eta, T0, previousT0 };
  }
  function uploadArrays({ h, u, Q, W, eta }, surfaceT, ice) {
    const packed = new Float32Array(OS.total);
    packed.set(h, OS.OH); packed.set(u, OS.OU); packed.set(Q, OS.OQ); packed.set(W, OS.OW);
    for (let e = 0; e < E; e++) if (!edgeOcean[e]) packed[OS.OU + e] = 0;
    device.queue.writeBuffer(ob.S, 0, packed);
    device.queue.writeBuffer(ob.T, 0, packed);
    const zero = new Float32Array(OS.total);
    for (const b of [ob.K1, ob.K2, ob.K3, ob.K4]) device.queue.writeBuffer(b, 0, zero);
    const odZero = new Float32Array(ODTOTAL);
    device.queue.writeBuffer(ob.OD, 0, odZero);
    device.queue.writeBuffer(ob.OD, 4 * OD.ETA, Float32Array.from(eta));
    device.queue.writeBuffer(ob.OD, 4 * OD.EMASK, Float32Array.from(edgeOcean));
    device.queue.writeBuffer(ob.OD, 4 * OD.CMASK, Float32Array.from(cellOcean));
    device.queue.writeBuffer(ob.OD, 4 * OD.BATH, Float32Array.from(D));
    device.queue.writeBuffer(ob.OD, 4 * OD.FEDGE, Float32Array.from(mesh.fEdge));
    const previousIce = Float64Array.from(ice);
    device.queue.writeBuffer(ob.OD, 4 * OD.PREVICE, Float32Array.from(previousIce));
    const capacity = Float64Array.from({ length: C }, (_, i) => o.density * o.specificHeat * Math.max(h[i], 1));
    core.uploadPhysics({ capacity });
  }
  function initialize(surfaceT, ice) {
    const arrays = initializeArrays(surfaceT, ice);
    uploadArrays(arrays, surfaceT, ice);
    device.queue.writeBuffer(ob.OD, 4 * OD.PREVT0, Float32Array.from(arrays.previousT0));
  }
  function loadArraysFromOldFormat(surfaceT, ice, saved) {
    const arrays = initializeArrays(surfaceT, ice);
    if (saved.h1 && saved.u1) {
      const { h } = arrays;
      const at = (k, i) => k * C + i;
      for (let i = 0; i < C; i++) {
        if (!cellOcean[i]) continue;
        const wanted = Math.max(o.minimumThickness, Math.min(saved.h1[i], D[i] - EPS * (L - 1)));
        let below = -1;
        for (let k = 1; k < L; k++) if (h[at(k, i)] > THIN) { below = k; break; }
        if (below > 0) {
          const delta = Math.max(-(h[i] - o.minimumThickness), Math.min(wanted - h[i], h[at(below, i)] - EPS));
          if (delta !== 0) {
            const from = delta > 0 ? below : 0, to = delta > 0 ? 0 : below, amt = Math.abs(delta);
            const ff = amt / h[at(from, i)];
            const dQ = arrays.Q[at(from, i)] * ff, dW = arrays.W[at(from, i)] * ff;
            h[at(from, i)] -= amt; arrays.Q[at(from, i)] -= dQ; arrays.W[at(from, i)] -= dW;
            h[at(to, i)] += amt; arrays.Q[at(to, i)] += dQ; arrays.W[at(to, i)] += dW;
          }
        }
      }
      for (let e = 0; e < E; e++) arrays.u[e] = edgeOcean[e] ? saved.u1[e] : 0;
    }
    return arrays;
  }
  function upload(saved, surfaceT, ice) {
    if (!saved || !saved.h || saved.h.length !== L * C) {
      const arrays = loadArraysFromOldFormat(surfaceT, ice, saved || {});
      uploadArrays(arrays, surfaceT, ice);
      device.queue.writeBuffer(ob.OD, 4 * OD.PREVT0, Float32Array.from(arrays.previousT0));
      return;
    }
    const h = Float64Array.from(saved.h), u = Float64Array.from(saved.u), eta = Float64Array.from(saved.eta);
    const Q = new Float64Array(L * C), W = new Float64Array(L * C);
    for (let n = 0; n < L * C; n++) { Q[n] = h[n] * saved.T[n]; W[n] = h[n] * saved.S[n]; }
    for (let e = 0; e < E; e++) if (!edgeOcean[e]) for (let k = 0; k < L; k++) u[k * E + e] = 0;
    uploadArrays({ h, u, Q, W, eta }, surfaceT, ice);
    const previousT0 = new Float64Array(C);
    for (let i = 0; i < C; i++) previousT0[i] = Q[i] / Math.max(EPS, h[i]);
    device.queue.writeBuffer(ob.OD, 4 * OD.PREVT0, Float32Array.from(previousT0));
  }

  async function download() {
    const s = await readBuffer(device, ob.S, 4 * OS.total);
    const od = await readBuffer(device, ob.OD, 4 * OD.total);
    const h = Float64Array.from(s.subarray(OS.OH, OS.OH + L * C));
    const u = Float64Array.from(s.subarray(OS.OU, OS.OU + L * E));
    const Qraw = s.subarray(OS.OQ, OS.OQ + L * C), Wraw = s.subarray(OS.OW, OS.OW + L * C);
    const T = new Float64Array(L * C), S = new Float64Array(L * C);
    for (let n = 0; n < L * C; n++) { const hh = Math.max(EPS, h[n]); T[n] = Qraw[n] / hh; S[n] = Wraw[n] / hh; }
    const eta = Float64Array.from(od.subarray(OD.ETA, OD.ETA + C));
    const h1 = h.subarray(0, C), T1 = T.subarray(0, C), S1 = S.subarray(0, C), u1 = u.subarray(0, E);
    const thermoclineDepth = new Float64Array(C);
    for (let i = 0; i < C; i++) thermoclineDepth[i] = cellOcean[i] ? h[i] + h[C + i] + h[2 * C + i] : NaN;
    return { h, u, T, S, eta, h1, T1, S1, u1, T2: thermoclineDepth, thermoclineDepth, layers: L };
  }
  function serializeFrom(d) {
    return { h: Array.from(d.h), u: Array.from(d.u), T: Array.from(d.T), S: Array.from(d.S), eta: Array.from(d.eta) };
  }
  async function serialize() { return serializeFrom(await download()); }

  function diagnosticsFrom(d) {
    let area = 0, depth = 0, heat = 0, thermo = 0, salinity = 0, ssh = 0, speed = 0, transport = 0, interiorT = 0, interiorH = 0, limited = 0;
    const rhoCp = o.density * o.specificHeat;
    const at = (k, i) => k * C + i, ae = (k, e) => k * E + e;
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) continue;
      const a = mesh.areaCell[i];
      area += a; depth += a * d.h[i]; salinity += a * d.S[i]; ssh = Math.max(ssh, Math.abs(d.eta[i]));
      thermo += a * (d.h[i] + d.h[at(1, i)] + d.h[at(2, i)]);
      for (let k = 0; k < L; k++) heat += a * rhoCp * (d.h[at(k, i)] * d.T[at(k, i)]);
      interiorT += d.h[at(1, i)] * d.T[at(1, i)] * a; interiorH += d.h[at(1, i)] * a;
    }
    for (let e = 0; e < E; e++) {
      speed = Math.max(speed, Math.abs(d.u[e]));
      if (Math.abs(Math.abs(d.u[e]) - SPEED_LIMIT) < 1e-4) limited++;
      const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1];
      const sill = Math.max(EPS, Math.min(D[a], D[b]) + 0.5 * (d.eta[a] + d.eta[b]));
      let sum = 0, t = 0;
      const hEdgeK = new Float64Array(L);
      hEdgeK[0] = 0.5 * (d.h[a] + d.h[b]);
      for (let k = 1; k < L; k++) hEdgeK[k] = Math.min(d.h[at(k, a)], d.h[at(k, b)]);
      for (let k = 0; k < L; k++) sum += hEdgeK[k];
      const f = sum > sill ? sill / sum : 1;
      for (let k = 0; k < L; k++) t += hEdgeK[k] * f * d.u[ae(k, e)];
      transport = Math.max(transport, Math.abs(t) * mesh.dvEdge[e]);
    }
    return { oceanUpperDepth: depth / area, oceanHeat: heat / area, oceanThermoclineT: interiorH > 0 ? interiorT / interiorH : 0, oceanSpeed: speed, oceanThermoclineDepth: thermo / area, oceanSalinity: salinity / area, oceanSSH: ssh, oceanTransport: transport / 1e6, oceanLimited: limited };
  }
  async function diagnostics() { return diagnosticsFrom(await download()); }

  return {
    layers: L, everySteps: o.everySteps, options: o,
    initialize, upload, download, serialize, serializeFrom, diagnostics, diagnosticsFrom,
    advance, advanceCoupled, accumulateFreshwater,
    setStress, readSurface, readSurfaceFromAtmosphere, stressFromAtmosphere, mixedLayer, salt, writeSurface, step,
    buffers: ob, layout: { OS, OD, B },
  };
}
