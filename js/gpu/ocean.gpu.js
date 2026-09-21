import { storageBuffer, emptyBuffer, readBuffer } from './device.module.js';

/*
 * The two-layer reduced-gravity ocean of ocean/reducedGravity.module.js
 * on the GPU. It keeps its own state buffer (h₁, h₂, u₁, u₂, h₁T₁,
 * h₂T₂), RK4 stages and scratch, and binds the atmosphere's state,
 * diagnostics and physics arrays beside them so the coupling kernels
 * can read the surface temperature, ice and lowest-layer wind and write
 * back the sea surface temperature, the upper layer's heat capacity and
 * the heat handed to the ice base. One tendency evaluation is the edge
 * flux, vertex PV, edge PV, cell (thickness, heat, Montgomery potential
 * plus kinetic energy, entrainment) and momentum kernels, with the ∇⁴
 * closure's two Laplacian passes before the momentum kernel; the ocean
 * steps once every `everySteps` atmosphere steps with the same RK4.
 */
const WORKGROUP = 64;
const FREEZING = 271.35;

export const OCEAN_DEFAULTS = {
  upperDepth: 50, lowerDepth: 350, reducedGravity: 0.02, abyssReducedGravity: 0.01, abyssTemperature: 275,
  minimumThickness: 10, entrainmentTime: 3600, density: 1025, specificHeat: 3985, interfacialDrag: 2e-4, bottomDrag: 2e-4,
  closureHours: 12, diffusivity: 0.3, everySteps: 4, dragCoefficient: 1.5e-3, gustiness: 3,
};

function oceanKernels(o) {
  const head = `
const OH1: i32 = ${o.OS.OH1}; const OH2: i32 = ${o.OS.OH2}; const OU1: i32 = ${o.OS.OU1}; const OU2: i32 = ${o.OS.OU2}; const OQ1: i32 = ${o.OS.OQ1}; const OQ2: i32 = ${o.OS.OQ2}; const OL: i32 = ${o.OS.total};
const O_FLUX: i32 = ${o.OD.FLUX}; const O_QV: i32 = ${o.OD.QV}; const O_QE: i32 = ${o.OD.QE}; const O_PHI: i32 = ${o.OD.PHI}; const O_STRESS: i32 = ${o.OD.STRESS}; const O_ICED: i32 = ${o.OD.ICED}; const O_EMASK: i32 = ${o.OD.EMASK}; const O_CMASK: i32 = ${o.OD.CMASK};
const O_DIVS: i32 = ${o.OD.DIVS}; const O_CURLS: i32 = ${o.OD.CURLS}; const O_LAPA: i32 = ${o.OD.LAPA}; const O_LAPB: i32 = ${o.OD.LAPB}; const O_T2: i32 = ${o.OD.T2};
const G12: f32 = ${o.reducedGravity}; const G23: f32 = ${o.abyssReducedGravity}; const TABYSS: f32 = ${o.abyssTemperature}; const HMIN: f32 = ${o.minimumThickness}; const TENTRAIN: f32 = ${o.entrainmentTime};
const RHO: f32 = ${o.density}; const RHOCP: f32 = ${o.density * o.specificHeat}; const RINT: f32 = ${o.interfacialDrag}; const RBOT: f32 = ${o.bottomDrag}; const NU4O: f32 = ${o.nu4}; const DIFFUSION: f32 = ${o.diffusion};
const TF: f32 = ${FREEZING}; const CDO: f32 = ${o.dragCoefficient}; const GUSTO: f32 = ${o.gustiness};
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
fn hOff(l: i32) -> i32 { return select(OH2, OH1, l == 0); }
fn uOff(l: i32) -> i32 { return select(OU2, OU1, l == 0); }
fn qOff(l: i32) -> i32 { return select(OQ2, OQ1, l == 0); }
`;
  const idx = `(i32(id.x) + i32(id.y) * 4194240)`;
  return {
    head,
    oRead: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = ${idx}; if (i >= C) { return; }
  let iced = S[S_ICE + i] > 0.0;
  OD[O_ICED + i] = select(0.0, 1.0, iced);
  let t1 = select(S[S_TS + i], TF, iced);
  IN[OQ1 + i] = IN[OH1 + i] * t1;
}`,
    oStress: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let e = ${idx}; if (e >= E) { return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  if (OD[O_EMASK + e] < 0.5 || OD[O_ICED + a] > 0.5 || OD[O_ICED + b] > 0.5) { OD[O_STRESS + e] = 0.0; return; }
  let bottom = (K - 1) * C;
  let rhoA = S[S_PI + a] * LV[L_SM + K - 1] / (RGAS * S[S_TH + bottom + a] * D[D_EXM + bottom + a]);
  let rhoB = S[S_PI + b] * LV[L_SM + K - 1] / (RGAS * S[S_TH + bottom + b] * D[D_EXM + bottom + b]);
  let fa = CDO * rhoA * max(D[D_WIND + a], GUSTO); let fb = CDO * rhoB * max(D[D_WIND + b], GUSTO);
  OD[O_STRESS + e] = 0.5 * (fa + fb) * S[S_U + (K - 1) * E + e];
}`,
    oFlux: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let n = ${idx}; if (n >= 2 * E) { return; }
  let l = n / E; let e = n % E;
  let h = hOff(l);
  let hEdge = 0.5 * (IN[h + MI[COE + 2 * e]] + IN[h + MI[COE + 2 * e + 1]]);
  OD[O_FLUX + n] = select(0.0, hEdge * IN[uOff(l) + e], OD[O_EMASK + e] > 0.5);
}`,
    oVertex: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let n = ${idx}; if (n >= 2 * V) { return; }
  let l = n / V; let v = n % V;
  var zeta = 0.0; var hv = 0.0;
  for (var m = 0; m < 3; m++) {
    let e = MI[EOV + 3 * v + m];
    zeta += f32(MI[ESV + 3 * v + m]) * IN[uOff(l) + e] * MF[F_DC + e];
    hv += MF[F_KAV + 3 * v + m] * IN[hOff(l) + MI[COV + 3 * v + m]];
  }
  OD[O_QV + n] = (zeta / MF[F_ATRI + v] + MF[F_FV + v]) / (hv / MF[F_ATRI + v]);
}`,
    oEdgePv: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let n = ${idx}; if (n >= 2 * E) { return; }
  let l = n / E; let e = n % E;
  OD[O_QE + n] = 0.5 * (OD[O_QV + l * V + MI[VOE + 2 * e]] + OD[O_QV + l * V + MI[VOE + 2 * e + 1]]);
}`,
    oCell: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let n = ${idx}; if (n >= 2 * C) { return; }
  let l = n / C; let i = n % C;
  let h = IN[hOff(l) + i]; let H = IN[qOff(l) + i]; let T = H / h;
  var div = 0.0; var heat = 0.0; var lap = 0.0; var kinetic = 0.0;
  for (var m = 0; m < MI[NEC + i]; m++) {
    let e = MI[EOC + MAXE * i + m]; let j = MI[COC + MAXE * i + m];
    let f = f32(MI[ESC + MAXE * i + m]) * OD[O_FLUX + l * E + e] * MF[F_DV + e];
    let Tj = IN[qOff(l) + j] / IN[hOff(l) + j];
    div += f; heat += f * 0.5 * (T + Tj);
    if (OD[O_EMASK + e] > 0.5) { lap += MF[F_DV + e] * (Tj - T) / MF[F_DC + e]; }
    let u = IN[uOff(l) + e];
    kinetic += 0.25 * MF[F_DC + e] * MF[F_DV + e] * u * u;
  }
  let area = MF[F_AREA + i];
  var dh = -div / area; var dH = -heat / area;
  if (l == 0) { dH += DIFFUSION * lap / area; }
  let h1 = IN[OH1 + i]; let h2 = IN[OH2 + i]; let t2 = IN[OQ2 + i] / h2;
  if (h1 < HMIN) { let w = (HMIN - h1) / TENTRAIN; if (l == 0) { dh += w; dH += w * t2; } else { dh -= w; dH -= w * t2; } }
  if (h2 < HMIN && l == 1) { let w = (HMIN - h2) / TENTRAIN; dh += w; dH += w * TABYSS; }
  let m2 = G23 * (h1 + h2);
  OD[O_PHI + n] = select(m2, m2 + G12 * h1, l == 0) + kinetic / area;
  OUT[hOff(l) + i] = dh; OUT[qOff(l) + i] = dH;
}`,
    oDivCurl: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let n = ${idx};
  let fromLap = P[1] > 0.5;
  if (n < 2 * C) {
    let l = n / C; let i = n % C;
    var sum = 0.0;
    for (var m = 0; m < MI[NEC + i]; m++) {
      let e = MI[EOC + MAXE * i + m];
      let u = select(IN[uOff(l) + e], OD[O_LAPA + l * E + e], fromLap);
      sum += f32(MI[ESC + MAXE * i + m]) * u * MF[F_DV + e];
    }
    OD[O_DIVS + n] = sum / MF[F_AREA + i];
  }
  if (n < 2 * V) {
    let l = n / V; let v = n % V;
    var sum = 0.0;
    for (var m = 0; m < 3; m++) {
      let e = MI[EOV + 3 * v + m];
      let u = select(IN[uOff(l) + e], OD[O_LAPA + l * E + e], fromLap);
      sum += f32(MI[ESV + 3 * v + m]) * u * MF[F_DC + e];
    }
    OD[O_CURLS + n] = sum / MF[F_ATRI + v];
  }
}`,
    oLapVelocity: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let n = ${idx}; if (n >= 2 * E) { return; }
  let l = n / E; let e = n % E;
  let lap = (OD[O_DIVS + l * C + MI[COE + 2 * e + 1]] - OD[O_DIVS + l * C + MI[COE + 2 * e]]) / MF[F_DC + e]
    - (OD[O_CURLS + l * V + MI[VOE + 2 * e + 1]] - OD[O_CURLS + l * V + MI[VOE + 2 * e]]) / MF[F_DV + e];
  if (P[1] > 0.5) { OD[O_LAPB + n] = lap; } else { OD[O_LAPA + n] = lap; }
}`,
    oMomentum: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let n = ${idx}; if (n >= 2 * E) { return; }
  let l = n / E; let e = n % E;
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  let qHere = 0.5 * OD[O_QE + n];
  var pv = 0.0;
  for (var s = 0; s < MI[NEE + e]; s++) {
    let slot = MAXEE * e + s; let other = MI[EOE + slot];
    pv += MF[F_PVW + slot] * OD[O_FLUX + l * E + other] * (qHere + 0.5 * OD[O_QE + l * E + other]);
  }
  let dc = MF[F_DC + e];
  let gradPhi = (OD[O_PHI + l * C + b] - OD[O_PHI + l * C + a]) / dc;
  var du = pv / dc - gradPhi;
  let hEdge = 0.5 * (IN[hOff(l) + a] + IN[hOff(l) + b]);
  let uDiff = IN[OU1 + e] - IN[OU2 + e];
  if (l == 0) { du += (OD[O_STRESS + e] / RHO - RINT * uDiff) / hEdge; }
  else { du += (RINT * uDiff - RBOT * IN[OU2 + e]) / hEdge; }
  du -= NU4O * OD[O_LAPB + n];
  OUT[uOff(l) + e] = select(0.0, du, OD[O_EMASK + e] > 0.5);
}`,
    oAdvance: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let n = ${idx}; if (n >= OL) { return; }
  OUT[n] = IN[n] + P[0] * OD[n];
}`,
    oCombine: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let n = ${idx}; if (n >= OL) { return; }
  IN[n] += P[0] * (OUT[n] + 2.0 * OD[n] + 2.0 * PH[n] + S[n]);
}`,
    oWrite: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = ${idx}; if (i >= C) { return; }
  if (OD[O_CMASK + i] < 0.5) { PH[PH_OFLUX + i] = 0.0; return; }
  let h1 = IN[OH1 + i]; let h2 = IN[OH2 + i];
  PH[PH_CAP + i] = RHOCP * max(h1, 1.0);
  OD[O_T2 + i] = IN[OQ2 + i] / h2;
  if (OD[O_ICED + i] > 0.5) {
    PH[PH_OFLUX + i] = RHOCP * (IN[OQ1 + i] - h1 * TF) / P[0];
    IN[OQ1 + i] = h1 * TF;
  } else {
    S[S_TS + i] = IN[OQ1 + i] / h1;
    PH[PH_OFLUX + i] = 0.0;
  }
}`,
  };
}

export function createGpuOcean(core, options = {}) {
  const o = { ...OCEAN_DEFAULTS, ...options };
  const { device, buffers, layout: L, C, E, V, mesh, meshSpacing } = core;
  const nu4 = o.closureHours > 0 ? Math.pow(meshSpacing / Math.PI, 4) / (o.closureHours * 3600) : 0;
  const diffusion = o.diffusivity * mesh.radius * mesh.radius / (o.density * o.specificHeat);
  const seq = (names) => { const out = {}; let off = 0; for (const [name, n] of names) { out[name] = off; off += n; } out.total = off; return out; };
  const OS = seq([['OH1', C], ['OH2', C], ['OU1', E], ['OU2', E], ['OQ1', C], ['OQ2', C]]);
  const OD = seq([['FLUX', 2 * E], ['QV', 2 * V], ['QE', 2 * E], ['PHI', 2 * C], ['STRESS', E], ['ICED', C], ['DIVS', 2 * C], ['CURLS', 2 * V], ['LAPA', 2 * E], ['LAPB', 2 * E], ['T2', C], ['EMASK', E], ['CMASK', C]]);
  const kernels = oceanKernels({ ...o, OS, OD, nu4, diffusion });
  const ob = {
    S: emptyBuffer(device, 4 * OS.total), T: emptyBuffer(device, 4 * OS.total),
    K1: emptyBuffer(device, 4 * OS.total), K2: emptyBuffer(device, 4 * OS.total), K3: emptyBuffer(device, 4 * OS.total), K4: emptyBuffer(device, 4 * OS.total),
    D: emptyBuffer(device, 4 * Math.max(OD.total, OS.total)),
  };
  for (const [name, b] of Object.entries(ob)) b.label = 'ocean' + name;
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
  function group(IN, OUT, D = ob.D, PH = buffers.PH, S = buffers.S) {
    const key = [IN.label, OUT.label, D.label, PH.label, S.label].join('|');
    let g = groups.get(key);
    if (!g) {
      g = device.createBindGroup({ layout: bindLayout, entries: [buffers.MI, buffers.MF, buffers.LV, IN, OUT, D, buffers.P, PH, S, buffers.D].map((buffer, binding) => ({ binding, resource: { buffer } })) });
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
    dispatch(pass, 'oFlux', g, 2 * E);
    dispatch(pass, 'oVertex', g, 2 * V);
    dispatch(pass, 'oEdgePv', g, 2 * E);
    dispatch(pass, 'oCell', g, 2 * C);
    dispatch(pass, 'oDivCurl', g, 2 * Math.max(C, V));
    dispatch(pass, 'oLapVelocity', g, 2 * E);
    pass.end();
    device.queue.submit([encoder.finish()]);
    setParams([0, 1]);
    encoder = device.createCommandEncoder(); pass = encoder.beginComputePass();
    dispatch(pass, 'oDivCurl', g, 2 * Math.max(C, V));
    dispatch(pass, 'oLapVelocity', g, 2 * E);
    dispatch(pass, 'oMomentum', g, 2 * E);
    pass.end();
    device.queue.submit([encoder.finish()]);
  }
  function advance(next, stage, factor) {
    setParams([factor]);
    const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
    dispatch(pass, 'oAdvance', group(ob.S, next, stage), OS.total);
    pass.end();
    device.queue.submit([encoder.finish()]);
  }

  let counter = 0;
  async function step(dt) {
    if (++counter % o.everySteps !== 0) return false;
    const dtOcean = o.everySteps * dt;
    const g = group(ob.S, ob.K1);
    {
      const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
      dispatch(pass, 'oRead', g, C);
      dispatch(pass, 'oStress', g, E);
      pass.end();
      device.queue.submit([encoder.finish()]);
    }
    tendency(ob.S, ob.K1); advance(ob.T, ob.K1, dtOcean / 2);
    tendency(ob.T, ob.K2); advance(ob.T, ob.K2, dtOcean / 2);
    tendency(ob.T, ob.K3); advance(ob.T, ob.K3, dtOcean);
    tendency(ob.T, ob.K4);
    setParams([dtOcean / 6]);
    {
      const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
      dispatch(pass, 'oCombine', group(ob.S, ob.K1, ob.K2, ob.K3, ob.K4), OS.total);
      pass.end();
      device.queue.submit([encoder.finish()]);
    }
    setParams([dtOcean]);
    {
      const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
      dispatch(pass, 'oWrite', g, C);
      pass.end();
      device.queue.submit([encoder.finish()]);
    }
    return true;
  }

  function upload({ h1, h2, u1, u2, T2 }, surfaceT, ice) {
    const packed = new Float32Array(OS.total);
    packed.set(h1, OS.OH1); packed.set(h2, OS.OH2); packed.set(u1, OS.OU1); packed.set(u2, OS.OU2);
    if (o.geography) for (let e = 0; e < E; e++) if (!o.geography.edgeOcean[e]) { packed[OS.OU1 + e] = 0; packed[OS.OU2 + e] = 0; }
    for (let i = 0; i < C; i++) {
      const t1 = ice[i] > 0 ? FREEZING : surfaceT[i];
      packed[OS.OQ1 + i] = h1[i] * t1;
      packed[OS.OQ2 + i] = h2[i] * T2[i];
    }
    device.queue.writeBuffer(ob.S, 0, packed);
    device.queue.writeBuffer(ob.T, 0, packed);
    device.queue.writeBuffer(ob.D, 0, new Float32Array(Math.max(OD.total, OS.total)));
    device.queue.writeBuffer(ob.D, 4 * OD.EMASK, Float32Array.from(o.geography ? o.geography.edgeOcean : new Uint8Array(E).fill(1)));
    device.queue.writeBuffer(ob.D, 4 * OD.CMASK, Float32Array.from(o.geography ? o.geography.land : new Uint8Array(C), (l) => (o.geography ? 1 - l : 1)));
    core.uploadPhysics({ capacity: Float64Array.from(h1, (h) => o.density * o.specificHeat * h) });
    counter = 0;
  }
  function initialize(surfaceT, ice) {
    const h1 = new Float64Array(C).fill(o.upperDepth), h2 = new Float64Array(C).fill(o.lowerDepth);
    const T2 = Float64Array.from(surfaceT, (t, i) => Math.max(FREEZING, o.abyssTemperature + 0.5 * ((ice[i] > 0 ? FREEZING : t) - o.abyssTemperature)));
    upload({ h1, h2, u1: new Float64Array(E), u2: new Float64Array(E), T2 }, surfaceT, ice);
  }
  async function download() {
    const s = await readBuffer(device, ob.S, 4 * OS.total);
    const h1 = Float64Array.from(s.subarray(OS.OH1, OS.OH1 + C)), h2 = Float64Array.from(s.subarray(OS.OH2, OS.OH2 + C));
    const u1 = Float64Array.from(s.subarray(OS.OU1, OS.OU1 + E)), u2 = Float64Array.from(s.subarray(OS.OU2, OS.OU2 + E));
    const T1 = Float64Array.from({ length: C }, (_, i) => s[OS.OQ1 + i] / h1[i]), T2 = Float64Array.from({ length: C }, (_, i) => s[OS.OQ2 + i] / h2[i]);
    return { h1, h2, u1, u2, T1, T2 };
  }
  return { step, upload, initialize, download, layout: { OS, OD }, options: o, everySteps: o.everySteps };
}
