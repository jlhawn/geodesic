/*
 * The column physics of the model as WGSL, one thread per column (or per
 * edge for momentum mixing), sharing the core's bindings and layouts:
 * the three-band gray radiation with clouds and the zenith/diffuse
 * surface reflection, bulk surface fluxes, the zero-layer sea ice, the
 * boundary-layer diagnosis, and the adjustment phase — boundary-layer
 * mixing by the implicit tridiagonal solve, saturation adjustment,
 * Betts–Miller convection with anvil detrainment, autoconversion, the
 * filler and the dry convective adjustment. Each is a line-by-line port
 * of the JavaScript module it names; the physics reads the Exner ratios
 * the last RK4 stage left in the diagnostic buffer, as the CPU does, and
 * the adjustment re-diagnoses the column first.
 */
export function physicsConstants(o) {
  return `
const S0: f32 = ${o.solarConstant}; const STEFAN: f32 = 5.670374419e-8; const LHEAT: f32 = ${o.latentHeat}; const EPSILON: f32 = 0.622; const RVAP: f32 = ${o.R / 0.622};
const CLOUD_ABS: f32 = ${o.cloudAbsorption}; const CLOUD_SCAT: f32 = ${o.cloudScattering}; const WINDOW: f32 = ${o.window}; const GAS_FRAC: f32 = ${o.gasFraction};
const VAPOR_FRAC: f32 = ${1 - o.window - o.gasFraction}; const OZONE_ABS: f32 = ${o.ozoneAbsorption}; const CEX: f32 = ${o.exchangeCoefficient};
const VCOUP: f32 = ${o.vaporCoupling}; const COUPLED: bool = ${o.vaporCoupling > 0}; const SKYLIGHT: f32 = ${o.skylight}; const DIFFUSE_MU: f32 = 0.6;
const ALB_ICE: f32 = ${o.iceAlbedo}; const FULLALB: f32 = ${o.fullAlbedoThickness}; const ALB_DIF_WATER: f32 = ${o.diffuseWaterAlbedo};
const FREEZING: f32 = 271.35; const MELTING: f32 = 273.15; const SKINC: f32 = ${o.skinHeatCapacity}; const COND: f32 = ${o.conductivity}; const HMIN: f32 = ${o.minimumThickness}; const LATENT_ICE: f32 = ${o.iceDensity * o.latentHeatFusion};
const RELAX: f32 = ${o.relaxationTime}; const RH_REF: f32 = ${o.referenceHumidity}; const AUTO_T: f32 = ${o.autoconversionThreshold}; const AUTO_R: f32 = ${o.autoconversionRate}; const CLOUD_LIFE: f32 = ${o.cloudLifetime};
const DETRAIN: f32 = ${o.detrainment}; const ANVIL: f32 = ${o.anvilDepth};
const RIC: f32 = ${o.richardsonCritical}; const KARMAN: f32 = ${o.vonKarman}; const KTOP: i32 = ${o.kTop};
`;
}

export const PHYSICS_FUNCTIONS = `
fn esat(T: f32) -> f32 { return 611.2 * exp(17.67 * (T - 273.15) / (T - 29.65)); }
fn qsat(T: f32, p: f32) -> f32 { let es = esat(T); let dry = p - (1.0 - EPSILON) * es; return select(1.0, EPSILON * es / dry, dry > 0.0); }
fn openWaterAlbedo(mu: f32) -> f32 { return 0.026 / (pow(mu, 1.7) + 0.065) + 0.15 * (mu - 0.1) * (mu - 0.5) * (mu - 1.0); }
fn surfaceAlbedo(h: f32, water: f32) -> f32 { return select(water, water + (ALB_ICE - water) * min(1.0, h / FULLALB), h > 0.0); }
fn cellWind(i: i32, k: i32) -> vec3<f32> {
  var w = vec3<f32>(0.0, 0.0, 0.0);
  for (var m = 0; m < MI[NEC + i]; m++) {
    let e = MI[EOC + MAXE * i + m];
    let s = 0.5 * MF[F_DC + e] * MF[F_DV + e] * IN[S_U + k * E + e];
    w += s * vec3<f32>(MF[F_NEDGE + 3 * e], MF[F_NEDGE + 3 * e + 1], MF[F_NEDGE + 3 * e + 2]);
  }
  return w / MF[F_AREA + i];
}
fn band(fraction: f32, eps: ptr<function, array<f32, K>>, temperature: ptr<function, array<f32, K>>, netFlux: ptr<function, array<f32, K>>, surfaceEmission: f32) -> vec2<f32> {
  var emitted: array<f32, K>;
  for (var k = 0; k < K; k++) { let t = (*temperature)[k]; emitted[k] = fraction * (*eps)[k] * STEFAN * t * t * t * t; }
  var carry = fraction * surfaceEmission;
  for (var k = K - 1; k >= 0; k--) { (*netFlux)[k] += (*eps)[k] * carry; carry *= 1.0 - (*eps)[k]; }
  var outgoing = carry; var back = 0.0;
  for (var k = 0; k < K; k++) {
    (*netFlux)[k] -= 2.0 * emitted[k];
    var down = emitted[k];
    for (var j = k + 1; j < K; j++) { (*netFlux)[j] += (*eps)[j] * down; down *= 1.0 - (*eps)[j]; }
    back += down;
    var up = emitted[k];
    for (var j = k - 1; j >= 0; j--) { (*netFlux)[j] += (*eps)[j] * up; up *= 1.0 - (*eps)[j]; }
    outgoing += up;
  }
  return vec2<f32>(outgoing, back);
}
fn moistLapse(temperature: f32, pressure: f32) -> f32 {
  let qs = qsat(temperature, pressure);
  return (RGAS * temperature + LHEAT * qs) / (CP + LHEAT * LHEAT * EPSILON * qs / (RGAS * temperature * temperature));
}
fn thomas(n: i32, upper: ptr<function, array<f32, K>>, lower: ptr<function, array<f32, K>>, rhs: ptr<function, array<f32, K>>) {
  var gain: array<f32, K>;
  var denominator = 1.0 + (*upper)[0] + (*lower)[0];
  gain[0] = -(*lower)[0] / denominator;
  (*rhs)[0] = (*rhs)[0] / denominator;
  for (var j = 1; j < n; j++) {
    denominator = 1.0 + (*upper)[j] + (*lower)[j] + (*upper)[j] * gain[j - 1];
    gain[j] = -(*lower)[j] / denominator;
    (*rhs)[j] = ((*rhs)[j] + (*upper)[j] * (*rhs)[j - 1]) / denominator;
  }
  for (var j = n - 2; j >= 0; j--) { (*rhs)[j] -= gain[j] * (*rhs)[j + 1]; }
}
`;

export const PHYSICS_KERNELS = {
  physics: `@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= C) { return; }
  let pi = IN[S_PI + i]; let ts = IN[S_TS + i]; let ice = IN[S_ICE + i]; let bottom = (K - 1) * C + i;
  let wind = cellWind(i, K - 1);
  let ws = length(wind);
  D[D_WIND + i] = ws;
  let sun = vec3<f32>(P[2], P[3], P[4]);
  let mu = max(0.0, MF[F_XC + 3 * i] * sun.x + MF[F_XC + 3 * i + 1] * sun.y + MF[F_XC + 3 * i + 2] * sun.z);
  let beam = S0 * mu;
  let adif = surfaceAlbedo(ice, ALB_DIF_WATER);
  let adir = surfaceAlbedo(ice, openWaterAlbedo(mu));
  let ozoneHeating = beam * OZONE_ABS;
  let surfaceEmission = STEFAN * ts * ts * ts * ts;
  var vaporE: array<f32, K>; var mixedE: array<f32, K>; var cloudE: array<f32, K>; var temperature: array<f32, K>; var netFlux: array<f32, K>;
  var cloudPath = 0.0;
  let tau0 = PH[PH_TAU + i];
  for (var k = 0; k < K; k++) {
    let idx = k * C + i;
    let mass = pi * LV[L_DS + k] / GRAV;
    var eps = 1.0 - exp(-tau0 * LV[L_SHAPE + k]);
    if (COUPLED) { eps = 1.0 - exp(-VCOUP * max(0.0, IN[S_Q + idx]) * mass); }
    let water = max(0.0, IN[S_QC + idx]) * mass;
    cloudPath += water;
    cloudE[k] = select(0.0, 1.0 - exp(-CLOUD_ABS * water), water > 0.0);
    let clear = 1.0 - cloudE[k];
    vaporE[k] = 1.0 - (1.0 - eps) * clear;
    mixedE[k] = 1.0 - (1.0 - LV[L_GASE + k]) * clear;
    temperature[k] = IN[S_TH + idx] * D[D_EXM + idx];
    netFlux[k] = ozoneHeating * LV[L_OZ + k];
  }
  let cloudDepth = CLOUD_SCAT * cloudPath;
  let reflectance = select(0.0, cloudDepth / (cloudDepth + 2.0 * mu), mu > 0.0 && cloudDepth > 0.0);
  let incident = beam - ozoneHeating;
  let direct = (1.0 - SKYLIGHT) * select(1.0, exp(-cloudDepth / mu), cloudDepth > 0.0 && mu > 0.0);
  let diffuse = 1.0 - reflectance - direct;
  let returned = select(0.0, cloudDepth / (cloudDepth + 2.0 * DIFFUSE_MU), cloudDepth > 0.0);
  let upward = adir * direct + adif * diffuse;
  let absorbed = incident * ((1.0 - adir) * direct + (1.0 - adif) * (diffuse + returned * upward / (1.0 - adif * returned)));
  let v = band(VAPOR_FRAC, &vaporE, &temperature, &netFlux, surfaceEmission);
  let g = band(GAS_FRAC, &mixedE, &temperature, &netFlux, surfaceEmission);
  let w = band(WINDOW, &cloudE, &temperature, &netFlux, surfaceEmission);
  let outgoing = v.x + g.x + w.x; let back = v.y + g.y + w.y;
  let airT = temperature[K - 1];
  let rho = pi * LV[L_SM + K - 1] / (RGAS * airT);
  let exchange = rho * CEX * max(ws, GUST);
  let sensible = exchange * CP * (ts - airT);
  let evap = max(0.0, exchange * (qsat(ts, pi) - IN[S_Q + bottom]));
  netFlux[K - 1] += sensible;
  let net = absorbed - surfaceEmission + back - sensible - LHEAT * evap;
  let dt = P[0];
  for (var k = 0; k < K; k++) {
    let idx = k * C + i;
    let mass = pi * LV[L_DS + k] / GRAV;
    IN[S_TH + idx] += dt * netFlux[k] / (CP * mass) / D[D_EXM + idx];
  }
  IN[S_Q + bottom] += dt * evap * GRAV / (pi * LV[L_DS + K - 1]);
  PH[PH_SFLUX + i] = net; PH[PH_ABS + i] = absorbed + ozoneHeating; PH[PH_OLR + i] = outgoing; PH[PH_SH + i] = sensible; PH[PH_EVAP + i] = evap; PH[PH_INS + i] = beam; PH[PH_REFL + i] = incident - absorbed; PH[PH_ADIF + i] = adif;
  let ocean = PH[PH_OFLUX + i]; let capacity = PH[PH_CAP + i];
  var T = ts; var h = ice;
  if (h <= 0.0) {
    T += dt * (net + ocean) / capacity;
    if (T < FREEZING) { h = (FREEZING - T) * capacity / LATENT_ICE; T = FREEZING; }
  } else {
    let conduction = COND * (FREEZING - T) / max(h, HMIN);
    T += dt * (net + conduction) / SKINC;
    var thickness = h + dt * (conduction - ocean) / LATENT_ICE;
    if (T > MELTING) { let excess = (T - MELTING) * SKINC; T = MELTING; thickness -= excess / LATENT_ICE; }
    if (thickness <= 0.0) { T = FREEZING + (-thickness * LATENT_ICE + SKINC * (T - FREEZING)) / capacity; h = 0.0; } else { h = thickness; }
  }
  IN[S_TS + i] = T; IN[S_ICE + i] = h;
}`,
  pblDiagnose: `@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= C) { return; }
  diagnoseColumn(i);
  let pi = IN[S_PI + i]; let base = (K - 1) * C + i;
  let bottomWind = cellWind(i, K - 1);
  let speed = length(bottomWind);
  let friction = sqrt(CDRAG) * max(speed, GUST);
  let zb = (D[D_GEO + base] + LV[L_GABS + K - 1]) / GRAV;
  var found = false; var riPrev = 0.0; var zPrev = zb; var depth = zb;
  for (var k = K - 2; k >= KTOP; k--) {
    if (found) { continue; }
    let idx = k * C + i;
    let z = (D[D_GEO + idx] + LV[L_GABS + k]) / GRAV;
    let dw = cellWind(i, k) - bottomWind;
    let shear = dot(dw, dw) + 100.0 * friction * friction;
    let ri = GRAV * (D[D_THV + idx] - D[D_THV + base]) * (z - zb) / (D[D_THV + base] * shear);
    if (ri > RIC) { depth = zPrev + (z - zPrev) * (RIC - riPrev) / (ri - riPrev); found = true; }
    else { riPrev = ri; zPrev = z; if (k == KTOP) { depth = z; } }
  }
  PH[PH_DEPTH + i] = depth;
  let h = depth - zb;
  for (var k = KTOP; k < K; k++) { PH[PH_MIX + k * C + i] = 0.0; }
  if (h <= 0.0) { return; }
  for (var k = KTOP; k < K - 1; k++) {
    let idx = k * C + i; let below = idx + C;
    let zAbove = (D[D_GEO + idx] + LV[L_GABS + k]) / GRAV; let zBelow = (D[D_GEO + below] + LV[L_GABS + k + 1]) / GRAV;
    let z = 0.5 * (zAbove + zBelow) - zb;
    if (z >= h) { continue; }
    let diffusivity = KARMAN * friction * z * (1.0 - z / h) * (1.0 - z / h);
    let rhoAbove = pi * LV[L_SM + k] / (RGAS * IN[S_TH + idx] * D[D_EXM + idx]);
    let rhoBelow = pi * LV[L_SM + k + 1] / (RGAS * IN[S_TH + below] * D[D_EXM + below]);
    PH[PH_MIX + idx] = 0.5 * (rhoAbove + rhoBelow) * diffusivity / (zAbove - zBelow);
  }
}`,
  adjust: `fn mixField(fieldOff: i32, i: i32, pi: f32, dt: f32) {
  var upper: array<f32, K>; var lower: array<f32, K>; var rhs: array<f32, K>;
  let n = K - KTOP;
  for (var j = 0; j < n; j++) {
    let k = KTOP + j;
    let mass = pi * LV[L_DS + k] / GRAV;
    upper[j] = select(0.0, dt * PH[PH_MIX + (k - 1) * C + i] / mass, j > 0);
    lower[j] = select(0.0, dt * PH[PH_MIX + k * C + i] / mass, k < K - 1);
    rhs[j] = IN[fieldOff + k * C + i];
  }
  thomas(n, &upper, &lower, &rhs);
  for (var j = 0; j < n; j++) { IN[fieldOff + (KTOP + j) * C + i] = rhs[j]; }
}
@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= C) { return; }
  let pi = IN[S_PI + i]; let dt = P[0]; let bottom = K - 1;
  var mixes = false;
  for (var k = KTOP; k < K - 1; k++) { if (PH[PH_MIX + k * C + i] > 0.0) { mixes = true; } }
  if (mixes) { mixField(S_TH, i, pi, dt); mixField(S_Q, i, pi, dt); mixField(S_QC, i, pi, dt); }
  diagnoseColumn(i);
  // saturation adjustment
  for (var k = 0; k < K; k++) {
    let idx = k * C + i;
    let ex = D[D_EXM + idx];
    let temperature = IN[S_TH + idx] * ex;
    let pressure = pi * LV[L_SM + k];
    let qs = qsat(temperature, pressure);
    let slope = qs * LHEAT / (RVAP * temperature * temperature);
    var change = (IN[S_Q + idx] - qs) / (1.0 + LHEAT * slope / CP);
    if (change < 0.0) { change = max(change, -IN[S_QC + idx]); }
    if (change == 0.0) { continue; }
    IN[S_Q + idx] -= change; IN[S_QC + idx] += change; IN[S_TH + idx] += LHEAT * change / (CP * ex);
  }
  // Betts–Miller convection
  var T: array<f32, K>; var p: array<f32, K>; var dp: array<f32, K>; var Tref: array<f32, K>; var qref: array<f32, K>;
  for (var k = 0; k < K; k++) { let idx = k * C + i; T[k] = IN[S_TH + idx] * D[D_EXM + idx]; p[k] = pi * LV[L_SM + k]; dp[k] = pi * LV[L_DS + k]; }
  var convected = 0.0;
  let Tp = T[bottom]; let pb = p[bottom];
  let qp = min(IN[S_Q + bottom * C + i], qsat(Tp, pb));
  var top = -1;
  if (qp > 0.0) {
    let e = qp * pb / (EPSILON + (1.0 - EPSILON) * qp);
    let y = log(e / 611.2);
    let dewPoint = (273.15 * 17.67 - 29.65 * y) / (17.67 - y);
    var lclT = Tp; var lclP = pb;
    if (dewPoint < Tp) { lclT = 1.0 / (1.0 / (dewPoint - 56.0) + log(Tp / dewPoint) / 800.0) + 56.0; lclP = pb * pow(lclT / Tp, 1.0 / KAPPA); }
    if (lclP >= p[0]) {
      var temperature = lclT; var pressure = lclP;
      for (var k = bottom; k >= 0; k--) {
        if (p[k] >= lclP) { Tref[k] = Tp * pow(p[k] / pb, KAPPA); }
        else {
          let dlnp = (log(p[k]) - log(pressure)) / 2.0;
          for (var n = 0; n < 2; n++) {
            let k1 = moistLapse(temperature, pressure);
            let k2 = moistLapse(temperature + 0.5 * dlnp * k1, pressure * exp(0.5 * dlnp));
            temperature += dlnp * k2;
            pressure *= exp(dlnp);
          }
          Tref[k] = temperature;
          if (Tref[k] > T[k]) { top = k; } else if (T[k] - Tref[k] > 10.0) { break; }
        }
      }
    }
  }
  if (top >= 0 && top != bottom) {
    for (var k = top; k <= bottom; k++) { qref[k] = RH_REF * qsat(Tref[k], p[k]); }
    var heating = 0.0; var drying = 0.0; var depth = 0.0;
    for (var k = top; k <= bottom; k++) { heating += CP * (Tref[k] - T[k]) * dp[k]; drying -= (qref[k] - IN[S_Q + k * C + i]) * dp[k]; depth += dp[k]; }
    if (heating > 0.0) {
      let rate = dt / RELAX;
      var rain = 0.0;
      if (drying > 0.0) {
        let shift = (LHEAT * drying - heating) / (CP * depth);
        for (var k = top; k <= bottom; k++) { Tref[k] += shift; }
        rain = drying / GRAV * rate;
      } else {
        let shiftQ = drying / depth; let shiftT = -heating / (CP * depth);
        for (var k = top; k <= bottom; k++) { qref[k] += shiftQ; Tref[k] += shiftT; }
      }
      for (var k = top; k <= bottom; k++) {
        let idx = k * C + i;
        IN[S_TH + idx] += (Tref[k] - T[k]) * rate / D[D_EXM + idx];
        IN[S_Q + idx] += (qref[k] - IN[S_Q + idx]) * rate;
      }
      if (rain > 0.0 && DETRAIN > 0.0) {
        var anvilMass = 0.0; var anvilBottom = top;
        for (var k = top; k <= bottom; k++) { if (k == top || anvilMass < ANVIL) { anvilMass += dp[k]; anvilBottom = k; } else { break; } }
        let detrained = DETRAIN * rain;
        for (var k = top; k <= anvilBottom; k++) { IN[S_QC + k * C + i] += detrained * GRAV / anvilMass; }
        rain -= detrained;
      }
      convected = rain;
    }
  }
  // autoconversion
  var rained = 0.0;
  for (var k = 0; k < K; k++) {
    let idx = k * C + i;
    let qc = IN[S_QC + idx];
    if (qc <= 0.0) { continue; }
    let excess = max(0.0, qc - AUTO_T);
    let converted = min(qc, excess * (1.0 - exp(-AUTO_R * dt)) + qc * (1.0 - exp(-dt / CLOUD_LIFE)));
    IN[S_QC + idx] = qc - converted;
    rained += pi * LV[L_DS + k] / GRAV * converted;
  }
  // filler
  for (var f = 0; f < 2; f++) {
    let off = select(S_Q, S_QC, f == 1);
    for (var k = 0; k < K - 1; k++) {
      let idx = k * C + i;
      if (IN[off + idx] < 0.0) { IN[off + idx + C] += IN[off + idx] * LV[L_DS + k] / LV[L_DS + k + 1]; IN[off + idx] = 0.0; }
    }
    if (IN[off + bottom * C + i] < 0.0) { IN[off + bottom * C + i] = 0.0; }
  }
  PH[PH_RAIN + i] += rained + convected; PH[PH_COND + i] += rained; PH[PH_CONV + i] += convected;
  // dry convective adjustment
  var dirty = true; var guard = 0;
  while (dirty && guard < K * K) {
    dirty = false; guard++;
    for (var k = K - 2; k >= 0; k--) {
      let above = k * C + i; let below = above + C;
      if (IN[S_TH + below] > IN[S_TH + above] * (1.0 + 1e-9)) {
        let wAbove = D[D_EXM + above] * LV[L_DS + k]; let wBelow = D[D_EXM + below] * LV[L_DS + k + 1];
        let mixed = (IN[S_TH + above] * wAbove + IN[S_TH + below] * wBelow) / (wAbove + wBelow);
        IN[S_TH + above] = mixed; IN[S_TH + below] = mixed;
        let mq = (IN[S_Q + above] * LV[L_DS + k] + IN[S_Q + below] * LV[L_DS + k + 1]) / (LV[L_DS + k] + LV[L_DS + k + 1]);
        IN[S_Q + above] = mq; IN[S_Q + below] = mq;
        let mc = (IN[S_QC + above] * LV[L_DS + k] + IN[S_QC + below] * LV[L_DS + k + 1]) / (LV[L_DS + k] + LV[L_DS + k + 1]);
        IN[S_QC + above] = mc; IN[S_QC + below] = mc;
        dirty = true;
      }
    }
  }
}`,
  mixMomentum: `@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let e = i32(id.x); if (e >= E) { return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  var mixes = false;
  for (var k = KTOP; k < K - 1; k++) { if (PH[PH_MIX + k * C + a] + PH[PH_MIX + k * C + b] > 0.0) { mixes = true; } }
  if (!mixes) { return; }
  let columnMass = 0.5 * (IN[S_PI + a] + IN[S_PI + b]); let dt = P[0];
  var upper: array<f32, K>; var lower: array<f32, K>; var rhs: array<f32, K>;
  let n = K - KTOP;
  for (var j = 0; j < n; j++) {
    let k = KTOP + j;
    let mass = columnMass * LV[L_DS + k] / GRAV;
    upper[j] = select(0.0, dt * 0.5 * (PH[PH_MIX + (k - 1) * C + a] + PH[PH_MIX + (k - 1) * C + b]) / mass, j > 0);
    lower[j] = select(0.0, dt * 0.5 * (PH[PH_MIX + k * C + a] + PH[PH_MIX + k * C + b]) / mass, k < K - 1);
    rhs[j] = IN[S_U + k * E + e];
  }
  thomas(n, &upper, &lower, &rhs);
  for (var j = 0; j < n; j++) { IN[S_U + (KTOP + j) * E + e] = rhs[j]; }
}`,
};
