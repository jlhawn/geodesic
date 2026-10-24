/* FLake (http://www.flake.igb-berlin.de/source/src_flake_sfcflx.tar.gz) is distributed under the MIT license; see the notice in js/physics/lake.module.js. */
import { FLAKE, LAKE_FIELDS, SETTLE_SLACK } from '../physics/lake.module.js';
import { SEA_SURFACE_WGSL } from './physics.gpu.js';

export const LAKE_SLOTS = ['depth', 'coriolis', ...LAKE_FIELDS, 'passed', 'pending'];
const slot = (name) => LAKE_SLOTS.indexOf(name);
const SNOW_FUSION = 3.34e5;

const literal = (v) => (Number.isInteger(v) ? v.toFixed(1) : String(v));
const k = FLAKE;

export const LAKE_WGSL = `const LK_CBL1: f32 = ${literal(k.cbl1)}; const LK_CBL2: f32 = ${literal(k.cbl2)}; const LK_ZMN: f32 = ${literal(k.zmN)}; const LK_ZMS: f32 = ${literal(k.zmS)}; const LK_ZMI: f32 = ${literal(k.zmI)};
const LK_RELAXH: f32 = ${literal(k.relaxH)}; const LK_RELAXC: f32 = ${literal(k.relaxC)}; const LK_CTMIN: f32 = ${literal(k.shapeMin)}; const LK_CTMAX: f32 = ${literal(k.shapeMax)};
const LK_PHI1: f32 = ${literal(k.phi1)}; const LK_PHI2: f32 = ${literal(k.phi2)}; const LK_CTT1: f32 = ${literal(k.ctt1)}; const LK_CTT2: f32 = ${literal(k.ctt2)};
const LK_ISHAPE: f32 = ${literal(k.iceShape)}; const LK_ISLOPE: f32 = ${literal(k.iceSlope)}; const LK_IAST: f32 = ${literal(k.iceAst)}; const LK_IMR: f32 = ${literal(k.iceMR)}; const LK_IMAX: f32 = ${literal(k.iceMax)};
const LK_IMIN: f32 = ${literal(k.iceMin)}; const LK_HMIN: f32 = ${literal(k.mixedMin)}; const LK_HMAX: f32 = ${literal(k.mixedMax)}; const LK_USMIN: f32 = ${literal(k.ustarMin)}; const LK_SMALL: f32 = ${literal(k.small)};
const LK_G: f32 = ${literal(k.g)}; const LK_TR: f32 = ${literal(k.maximumDensityT)}; const LK_TF: f32 = ${literal(k.freezing)}; const LK_AT: f32 = ${literal(k.expansion)}; const LK_RHOW: f32 = ${literal(k.rhoWater)}; const LK_RHOI: f32 = ${literal(k.rhoIce)};
const LK_LF: f32 = ${literal(k.fusion)}; const LK_CW: f32 = ${literal(k.cWater)}; const LK_CI: f32 = ${literal(k.cIce)}; const LK_KW: f32 = ${literal(k.kWater)}; const LK_KI: f32 = ${literal(k.kIce)}; const LK_EXT: f32 = ${literal(k.extinction)};
const LK_SNOW_FUSION: f32 = ${literal(SNOW_FUSION)}; const LK_SLACK: f32 = ${literal(SETTLE_SLACK)};
fn lkAt(n: i32, i: i32) -> i32 { return PH_LAKE + n * C + i; }
fn lkBuoyancy(T: f32) -> f32 { return LK_G * LK_AT * (T - LK_TR); }
fn lkEnthalpy(D: f32, mean: f32, ice: f32, iceT: f32) -> f32 {
  let CI = LK_ISHAPE - LK_IMR * (1.0 + LK_IAST) * min(1.0, ice / LK_IMAX);
  return LK_RHOW * LK_CW * D * (mean - LK_TF) - LK_RHOI * ice * (LK_LF + LK_CI * CI * (LK_TF - iceT));
}
fn lakeSettle(i: i32, D: f32, goal: f32) {
  let rc = LK_RHOW * LK_CW;
  let left = goal - rc * D * (PH[lkAt(${slot('mean')}, i)] - LK_TF); let need = -left / LK_RHOI;
  if (need <= LK_LF * LK_IMIN + select(rc * D * LK_SLACK / LK_RHOI, 0.0, PH[lkAt(${slot('ice')}, i)] > 0.0)) {
    let shift = left / rc / D;
    PH[lkAt(${slot('mixed')}, i)] += shift; PH[lkAt(${slot('mean')}, i)] += shift; PH[lkAt(${slot('bottom')}, i)] += shift; PH[lkAt(${slot('ice')}, i)] = 0.0; PH[lkAt(${slot('iceT')}, i)] = LK_TF;
    return;
  }
  let ice = PH[lkAt(${slot('ice')}, i)];
  let deficit = select(0.0, LK_TF - PH[lkAt(${slot('iceT')}, i)], ice > 0.0); let slope = LK_IMR * (1.0 + LK_IAST) / LK_IMAX;
  let a = LK_LF + LK_CI * deficit * LK_ISHAPE; let b = LK_CI * deficit * slope; let root = a * a - 4.0 * b * need;
  var h = select(LK_IMAX + 1.0, 2.0 * need / (a + sqrt(root)), root > 0.0);
  if (h > LK_IMAX) {
    h = LK_IMAX;
    PH[lkAt(${slot('iceT')}, i)] = LK_TF - (need / h - LK_LF) / (LK_CI * (LK_ISHAPE - slope * h));
  } else if (ice == 0.0) { PH[lkAt(${slot('iceT')}, i)] = LK_TF; }
  PH[lkAt(${slot('ice')}, i)] = h;
}
fn lakeIceAlbedo(T: f32) -> f32 { let blue = exp(-95.6 * (LK_TF - T) / LK_TF); return 0.6 * (1.0 - blue) + 0.1 * blue; }
`;


const LAKE_STEP_WGSL = `fn lakeStep(i: i32, absorbed: f32, flux: f32, ustar: f32, dt: f32) -> vec2<f32> {
  let D = PH[lkAt(${slot('depth')}, i)]; let f = PH[lkAt(${slot('coriolis')}, i)];
  if (PH[lkAt(${slot('mean')}, i)] == 0.0) {
    let t = max(LK_TF, IN[S_TS + i]);
    PH[lkAt(${slot('mixed')}, i)] = t; PH[lkAt(${slot('mean')}, i)] = t; PH[lkAt(${slot('bottom')}, i)] = t; PH[lkAt(${slot('shape')}, i)] = LK_CTMIN;
    PH[lkAt(${slot('mixedDepth')}, i)] = D; PH[lkAt(${slot('ice')}, i)] = 0.0; PH[lkAt(${slot('iceT')}, i)] = LK_TF;
  }
  let Tml = PH[lkAt(${slot('mixed')}, i)]; let Tmn = PH[lkAt(${slot('mean')}, i)]; let Tb = PH[lkAt(${slot('bottom')}, i)]; let CT = PH[lkAt(${slot('shape')}, i)];
  let h = PH[lkAt(${slot('mixedDepth')}, i)]; let hi = PH[lkAt(${slot('ice')}, i)]; let Ti = PH[lkAt(${slot('iceT')}, i)];
  let goal = lkEnthalpy(D, Tmn, hi, Ti) + (absorbed + flux) * dt;
  let rc = LK_RHOW * LK_CW;
  let iced = hi >= LK_IMIN;
  let Iw = select(absorbed, 0.0, iced); let Iice = select(0.0, absorbed, iced);
  let eh = exp(-LK_EXT * h); let eD = exp(-LK_EXT * D);
  let Ih = select(Iw, Iw * eh, h >= LK_HMIN); let Ibot = Iw * eD;
  let I0h = select(Ih, Iw * (1.0 - eh) / LK_EXT / h, h >= LK_HMIN);
  let IhD = select(Ih, Iw * (eh - eD) / LK_EXT / (D - h), h <= D - LK_HMIN);
  var Qw = select(flux, 0.0, iced);
  let Qice = select(0.0, flux, iced);
  if (iced) { Qw = select(0.0, -LK_KW * (Tb - Tml) / D * max(LK_PHI1 * CT - LK_PHI2, 1.0), h <= LK_HMIN); }
  let Qstar = Qw + Iw + Ih - 2.0 * I0h;
  var hiN = hi; var TiN = Ti; var hN = h; var CTN = CT; var TbN = Tb; var TmlN = Tml; var create = false;
  if (!iced) {
    create = Tml <= LK_TF + LK_SMALL && Qw < 0.0;
    if (create) {
      hiN = hi - Qw / LK_RHOI / LK_LF * dt;
      TiN = LK_TF + hiN * Qw / LK_KI / LK_ISLOPE;
    }
  } else {
    var melting = false;
    if (Ti >= LK_TF - LK_SMALL) {
      let gain = Qice + Iice - Iw - Qw;
      if (gain >= 0.0) { melting = true; hiN = hi - gain / LK_LF / LK_RHOI * dt; TiN = LK_TF; }
    }
    if (!melting) {
      let ratio = hi / LK_IMAX; let CI = LK_ISHAPE - LK_IMR * (1.0 + LK_IAST) * ratio; let slope = LK_ISLOPE - ratio;
      var threshold = max(1.0, 2.0 * CI * LK_CI * (LK_TF - Ti) / LK_LF);
      threshold = sqrt(slope / CI * LK_KI / LK_RHOI / LK_CI * threshold * dt);
      threshold = min(0.9 * LK_IMAX, max(threshold, LK_IMIN));
      if (hi < threshold) {
        let top = Qice + Iice - Iw;
        hiN = hi - (top - Qw) / LK_LF / LK_RHOI * dt;
        TiN = LK_TF + hiN * top / LK_KI / slope;
      } else {
        let conducted = LK_KI * (LK_TF - Ti) / hi * slope;
        hiN = hi + (Qw + conducted) / LK_LF / LK_RHOI * dt;
        let star = (1.0 - CI) * LK_CI * (LK_TF - Ti) / LK_LF;
        let rate = (conducted * (1.0 - star) - star * Qw + (Qice + Iice - Iw)) / LK_RHOI / LK_CI / (CI * hi);
        TiN = Ti + rate * dt;
      }
    }
  }
  hiN = min(hiN, LK_IMAX);
  TiN = max(73.15, min(TiN, LK_TF));
  if (hiN < LK_IMIN) { hiN = 0.0; TiN = LK_TF; create = false; }
  if (create) { Qw = 0.0; }
  var TmnN = max(LK_TF, Tmn + (Qw + Iw) / rc / D * dt);
  if (hiN >= LK_IMIN) {
    TmnN = min(TmnN, LK_TR);
    TmlN = LK_TF;
    if (create) {
      if (h >= D - LK_HMIN) { hN = 0.0; CTN = LK_CTMIN; } else { hN = h; CTN = CT; }
      TbN = TmlN - (TmlN - TmnN) / CTN / (1.0 - hN / D);
    } else if (Tb < LK_TR) {
      TbN = TmlN - (TmlN - TmnN) / CTN / (1.0 - hN / D);
    } else {
      TbN = LK_TR;
      if (h >= LK_SMALL) { hN = max(0.0, D * (1.0 - (TmlN - TmnN) / (TmlN - TbN) / CT)); }
      else { CTN = min(LK_CTMAX, max((TmlN - TmnN) / (TmlN - TbN), LK_CTMIN)); }
    }
    TbN = min(TbN, LK_TR);
  } else {
    let B = lkBuoyancy(Tml) * Qstar / rc;
    let wstar = select(0.0, pow(-B * h, 1.0 / 3.0), B < 0.0 && h > 0.0);
    var equilibrium = -Qw / max(Iw, LK_SMALL);
    if (equilibrium > 0.0 && equilibrium < 1.0 && Tml > LK_TR) { equilibrium = min(D, (sqrt(6.0 * equilibrium) + 2.0 * equilibrium / (1.0 - equilibrium)) / LK_EXT); }
    else { equilibrium = 0.0; }
    let N = select(0.0, sqrt(max(0.0, lkBuoyancy(0.5 * (Tml + Tb)) * (Tml - Tb)) / (D - h)), h <= D - LK_HMIN);
    let mixing = max(max(wstar, ustar), LK_USMIN);
    var dCT = (LK_CTMAX - LK_CTMIN) / max(N * (D - h) * (D - h) / LK_RELAXC / (mixing * mixing), LK_SMALL);
    let CTT = LK_CTT1 * CT - LK_CTT2; let CQ = 2.0 * CTT / CT;
    var dh = 0.0;
    if (B < 0.0) {
      CTN = min(LK_CTMAX, max(CT + dCT * dt, LK_CTMIN));
      dCT = (CTN - CT) / dt;
      if (h <= D - LK_HMIN) {
        if (h <= LK_HMIN) { dh = LK_CBL1 / LK_CBL2 * max(wstar, LK_SMALL); }
        else {
          let RH = D / h; let Rr = RH - 1.0; let RT = CT / CTT; let Rs = (RT / 2.0 - 1.0) * Rr + 1.0;
          dh = -Qstar * (Rs * (1.0 + LK_CBL1) - 1.0) / rc;
          dh += (D - h) * (Tml - Tb) * LK_CTT2 / CTT * dCT;
          dh += (Ibot + (RT - 1.0) * Ih - RT * IhD + (RT - 2.0) * Rr * (Ih - I0h)) / rc;
          dh /= -LK_CBL2 * Rs * Qstar / rc / max(wstar, LK_SMALL) + CT * (Tml - Tb);
        }
        dh = max(dh, LK_SMALL);
        hN = max(LK_HMIN, min(h + dh * dt, D));
      } else { hN = D; }
    } else {
      let u = max(ustar, LK_USMIN);
      var scale = max((abs(f) / LK_ZMN + N / LK_ZMI) * u * u + B / LK_ZMS, LK_SMALL);
      scale = max(max(LK_HMIN, min(u * u * u / scale, LK_HMAX)), equilibrium);
      hN = max(LK_HMIN, min(scale - (scale - h) * exp(-LK_RELAXH * u / scale * dt), D));
      if (hN <= h) { dCT = -dCT; }
      CTN = min(LK_CTMAX, max(CT + dCT * dt, LK_CTMIN));
      dCT = (CTN - CT) / dt;
    }
    dh = (hN - h) / dt;
    if (hN <= D - LK_HMIN) {
      var dTb = 0.0;
      if (hN > h) {
        let RH = h / D; let Rr = 1.0 - RH; let R0 = 0.5 * CT * Rr + CTT * (2.0 * RH - 1.0); let Rs = (0.5 + CTT - CQ) / R0; let RT = (1.0 - CT * Rr) / R0;
        dTb = ((Qw + Iw - Ibot) / rc - CT * (Tml - Tb) * dh) * Rs / D;
        dTb += (IhD - (1.0 - CQ) * Ih - CQ * Ibot) * RT / (D - h) / rc;
        dTb += (1.0 - LK_CTT2 * RT) / CT * (Tml - Tb) * dCT;
      }
      TbN = max(LK_TF, Tb + dTb * dt);
      if ((TbN - LK_TR) * lkBuoyancy(TmnN) < 0.0) { TbN = LK_TR; }
      let w = CTN * (1.0 - hN / D);
      TmlN = max(LK_TF, (TmnN - TbN * w) / (1.0 - w));
    } else { hN = D; TmlN = TmnN; TbN = TmnN; CTN = LK_CTMIN; }
  }
  if ((TmlN - TbN) * lkBuoyancy(TmnN) < 0.0) { hN = D; TmlN = TmnN; TbN = TmnN; CTN = LK_CTMIN; }
  PH[lkAt(${slot('mixed')}, i)] = TmlN; PH[lkAt(${slot('mean')}, i)] = TmnN; PH[lkAt(${slot('bottom')}, i)] = TbN; PH[lkAt(${slot('shape')}, i)] = CTN;
  PH[lkAt(${slot('mixedDepth')}, i)] = hN; PH[lkAt(${slot('ice')}, i)] = hiN; PH[lkAt(${slot('iceT')}, i)] = TiN;
  lakeSettle(i, D, goal);
  let iceN = PH[lkAt(${slot('ice')}, i)];
  return vec2<f32>(select(PH[lkAt(${slot('mixed')}, i)], PH[lkAt(${slot('iceT')}, i)], iceN > 0.0), iceN);
}
`;

const ON_LAKE = 'PH[PH_LAND + i] < -0.5';

function replaceOnce(body, from, to, where) {
  const at = body.indexOf(from);
  if (at < 0 || body.indexOf(from, at + 1) >= 0) throw new Error(`the lake's ${where} anchor is not in its kernel exactly once`);
  return body.slice(0, at) + to + body.slice(at + from.length);
}

const iceAlbedo = (body, name, open, where) => replaceOnce(body, `let ${name} = ${open};`, `let ${name} = select(${open}, lakeIceAlbedo(skin), ${ON_LAKE} && ice > 0.0);`, where);
const ADIF = 'select(cover * iceDif + (1.0 - cover) * ALB_DIF_WATER, landAlbedo, onLand)';
const ADIR = 'select(cover * iceDir + (1.0 - cover) * waterDir, landAlbedo, onLand)';
const DEPOSIT = '  if (PH[PH_LAND + i] < 0.5 && airT < MELTING && rained + convected > 0.0) {';
const SKIN_STEP = '  if (!onLand) { skinStep(';

export function lakeKernels(bodies, { held = false, skin = false } = {}) {
  const physics = iceAlbedo(iceAlbedo(bodies.physics, 'adif', ADIF, 'diffuse albedo'), 'adir', ADIR, 'direct albedo');
  const radiation = held ? iceAlbedo(bodies.radiation, 'adir', ADIR, 'held direct albedo') : bodies.radiation;
  let surface = replaceOnce(bodies.physicsSurface, `} else ${SEA_SURFACE_WGSL}`, `} else if (${ON_LAKE}) {
    let lake = lakeStep(i, absorbed, net - absorbed, sqrt(rho / LK_RHOW * PH[PH_DRAG + i]) * xWind(i, ws), dt);
    T = lake.x; h = lake.y;
    PH[PH_CONC + i] = select(0.0, 1.0, h > 0.0);
    PH[lkAt(${slot('passed')}, i)] -= evap * dt; PH[lkAt(${slot('pending')}, i)] -= evap * dt;
${skin ? '    PH[PH_WARM + i] = 0.0; PH[PH_COOL + i] = 0.0;\n' : ''}  } else ${SEA_SURFACE_WGSL}`, 'surface');
  if (skin) surface = replaceOnce(surface, SKIN_STEP, `  if (!onLand && !(${ON_LAKE})) { skinStep(`, 'skin');
  const rain = replaceOnce(bodies.adjustRain, DEPOSIT, `  if (${ON_LAKE}) {
    PH[lkAt(${slot('passed')}, i)] += rained + convected; PH[lkAt(${slot('pending')}, i)] += rained + convected;
    if (airT < MELTING && rained + convected > 0.0) {
      if (PH[lkAt(${slot('mean')}, i)] != 0.0) { let depth = PH[lkAt(${slot('depth')}, i)]; lakeSettle(i, depth, lkEnthalpy(depth, PH[lkAt(${slot('mean')}, i)], PH[lkAt(${slot('ice')}, i)], PH[lkAt(${slot('iceT')}, i)]) - LK_SNOW_FUSION * (rained + convected)); }
      IN[S_TH + bottom * C + i] += LFUS * (rained + convected - snowed) * GRAV / (CP * pi * LV[L_DS + K - 1] * D[D_EXM + bottom * C + i]);
    }
  } else${DEPOSIT.slice(1)}`, 'deposit');
  return {
    physics: LAKE_WGSL + physics,
    radiation: held ? LAKE_WGSL + radiation : radiation,
    physicsSurface: LAKE_WGSL + LAKE_STEP_WGSL + surface,
    adjustRain: LAKE_WGSL + rain,
  };
}

export const LAKE_DRAIN_WGSL = `      if (PH[PH_LAND + j] < -0.5) { let lakeIndex = PH_LAKE + ${slot('pending')} * C + j; fresh -= PH[lakeIndex] * OD[O_DRAINW + p]; PH[lakeIndex] = 0.0; }
`;

export function lakeUpload(lake, C) {
  const out = new Float32Array(LAKE_SLOTS.length * C);
  LAKE_SLOTS.forEach((name, n) => { if (name !== 'pending') out.set(lake[name], n * C); });
  return out;
}

export function lakeDownload(lake, values, C) {
  LAKE_SLOTS.forEach((name, n) => { if (name !== 'pending' && name !== 'depth' && name !== 'coriolis') lake[name].set(values.subarray(n * C, (n + 1) * C)); });
}
