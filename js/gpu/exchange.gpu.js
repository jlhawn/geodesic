import { EXCHANGE_DEFAULTS, ROUGHNESS, ANDREAS, KARMAN, FAO_KARMAN } from '../physics/exchange.module.js';

/*
 * The surface layer of js/physics/exchange.module.js in WGSL, line by
 * line: surfaceExchange(i, ...) returns C_D, C_H and the FAO-56 reference
 * C_H of column i.
 */
export function exchangeConstants(o) {
  if (o.surfaceExchange !== 'roughness' && o.surfaceExchange !== 'fixed') throw new Error(`exchange must be 'roughness' or 'fixed', not ${o.surfaceExchange}`);
  const x = { ...EXCHANGE_DEFAULTS, ...o.exchangeOptions };
  const r = { ...ROUGHNESS, ...x.roughness };
  const [slope, offset, ceiling] = x.charnock;
  const gusty = o.surfaceExchange === 'roughness' && !!x.convectiveGust, [seaGust, landGust, gustFloor] = gusty ? x.convectiveGust : [0, 0, 0];
  return `
const ROUGH: bool = ${o.surfaceExchange === 'roughness'}; const IMPLICIT_DRAG: bool = ${o.surfaceExchange === 'roughness' && !!o.implicitDrag};
const XKARMAN: f32 = ${KARMAN}; const XGRAV: f32 = 9.81; const X_ITER: i32 = ${x.iterations}; const X_BLEND: f32 = ${x.blendingHeight}; const X_SNOWCOVER: f32 = ${x.snowCover};
const Z0_FOREST: f32 = ${r.forest[0]}; const ZH_FOREST: f32 = ${r.forest[1]}; const Z0_GRASS: f32 = ${r.grass[0]}; const ZH_GRASS: f32 = ${r.grass[1]}; const Z0_BARE: f32 = ${r.bare[0]}; const ZH_BARE: f32 = ${r.bare[1]}; const Z0_SNOW: f32 = ${r.snow};
const CH_SLOPE: f32 = ${slope}; const CH_OFFSET: f32 = ${offset}; const CH_CEIL: f32 = ${ceiling}; const SMOOTH: f32 = ${x.smoothFlow}; const REF_CROP: f32 = ${x.referenceCrop}; const REF_KARMAN: f32 = ${FAO_KARMAN};
const AND_SMOOTH: f32 = ${ANDREAS[0][0]}; const AND_TRANS: f32 = ${ANDREAS[1][0]};
const X_GUSTY: bool = ${gusty}; const X_GUST_ITER: i32 = ${x.gustIterations}; const X_GUST_SEA: f32 = ${seaGust}; const X_GUST_LAND: f32 = ${landGust}; const X_GUST_FLOOR: f32 = ${gustFloor}; const X_WET_SURFACE: bool = ${x.landHumidity === 'wetness'};
fn xWind(i: i32, wind: f32) -> f32 { return select(max(wind, GUST), PH[PH_XWIND + i], X_GUSTY); }
const AND0: vec3<f32> = vec3<f32>(${ANDREAS[0].slice(1).join(', ')}); const AND1: vec3<f32> = vec3<f32>(${ANDREAS[1].slice(1).join(', ')}); const AND2: vec3<f32> = vec3<f32>(${ANDREAS[2].slice(1).join(', ')});
`;
}

export const EXCHANGE_WGSL = `
fn xPsiM(zeta: f32) -> f32 {
  if (zeta < 0.0) { let x = pow(1.0 - 16.0 * zeta, 0.25); return 1.5707963267948966 - 2.0 * atan(x) + log((1.0 + x) * (1.0 + x) * (1.0 + x * x) / 8.0); }
  return -(2.0 / 3.0) * (zeta - 5.0 / 0.35) * exp(-0.35 * zeta) - zeta - (2.0 / 3.0) * 5.0 / 0.35;
}
fn xPsiH(zeta: f32) -> f32 {
  if (zeta < 0.0) { return 2.0 * log((1.0 + sqrt(1.0 - 16.0 * zeta)) / 2.0); }
  return -(2.0 / 3.0) * (zeta - 5.0 / 0.35) * exp(-0.35 * zeta) - pow(1.0 + 2.0 * zeta / 3.0, 1.5) - (2.0 / 3.0) * 5.0 / 0.35 + 1.0;
}
fn xSlopeM(zeta: f32) -> f32 { return -(2.0 / 3.0) * exp(-0.35 * zeta) * (6.0 - 0.35 * zeta) - 1.0; }
fn xSlopeH(zeta: f32) -> f32 { return -(2.0 / 3.0) * exp(-0.35 * zeta) * (6.0 - 0.35 * zeta) - sqrt(1.0 + 2.0 * zeta / 3.0); }
fn xProfiles(zeta: f32, z: f32, z0m: f32, z0h: f32) -> vec2<f32> {
  let top = zeta * (1.0 + z0m / z);
  return vec2<f32>(log((z + z0m) / z0m) - xPsiM(top) + xPsiM(zeta * z0m / z), log((z + z0m) / z0h) - xPsiH(top) + xPsiH(zeta * z0h / z));
}
fn xTransfer(ri: f32, z: f32, z0m: f32, z0h: f32) -> vec2<f32> {
  var f = xProfiles(0.0, z, z0m, z0h);
  var zeta = ri * f.x * f.x / f.y;
  if (ri > 0.0) {
    let goal = log(ri);
    var s = log(zeta);
    for (var n = 0; n < X_ITER; n++) {
      zeta = exp(s);
      f = xProfiles(zeta, z, z0m, z0h);
      let top = zeta * (1.0 + z0m / z); let lowM = zeta * z0m / z; let lowH = zeta * z0h / z;
      let dm = -xSlopeM(top) * top + xSlopeM(lowM) * lowM; let dh = -xSlopeH(top) * top + xSlopeH(lowH) * lowH;
      s -= (s - goal - log(f.x * f.x / f.y)) / (1.0 - 2.0 * dm / f.x + dh / f.y);
    }
    zeta = exp(s);
  } else {
    for (var n = 0; n < X_ITER; n++) { f = xProfiles(zeta, z, z0m, z0h); zeta = ri * f.x * f.x / f.y; }
  }
  f = xProfiles(zeta, z, z0m, z0h);
  return vec2<f32>(XKARMAN * XKARMAN / (f.x * f.x), XKARMAN * XKARMAN / (f.x * f.y));
}
fn xAndreas(z0: f32, ustar: f32, viscosity: f32) -> f32 {
  let reynolds = ustar * z0 / viscosity; let l = log(min(1000.0, reynolds));
  var b = AND2;
  if (reynolds <= AND_TRANS) { b = AND1; }
  if (reynolds <= AND_SMOOTH) { b = AND0; }
  return z0 * exp(b.x + b.y * l + b.z * l * l);
}
struct XBlend { drag: f32, heat: f32 }
fn xAdd(blend: ptr<function, XBlend>, share: f32, z0m: f32, z0h: f32) {
  if (!(share > 0.0)) { return; }
  let lm = log(X_BLEND / z0m);
  (*blend).drag += share * XKARMAN * XKARMAN / (lm * lm);
  (*blend).heat += share * XKARMAN * XKARMAN / (lm * log(X_BLEND / z0h));
}
fn xWetness(aero: f32, roots: f32, bareWet: f32, veg: f32, snow: f32, warmth: f32) -> f32 {
  let canopyWet = veg * roots / (1.0 + RSTOM * aero / max(0.05, warmth));
  return select(select(roots, bareWet + canopyWet, VEGETATED), 1.0, snow > 0.0);
}
fn surfaceExchange(i: i32, pi: f32, skin: f32, wind: f32, concentration: f32, snow: f32, cover: f32, trees: f32, onLand: bool, onIceSheet: bool, wetness: f32, depth: f32) -> vec3<f32> {
  let b = (K - 1) * C + i;
  let z = CP * D[D_THV + b] * (D[D_EXL + b] - D[D_EXM + b]) / GRAV;
  let scale = select(X_GUST_SEA, X_GUST_LAND, onLand); let mixed = max(depth, z);
  var speed = select(max(wind, GUST), sqrt(wind * wind + X_GUST_FLOOR * X_GUST_FLOOR), X_GUSTY);
  let celsius = IN[S_TH + b] * D[D_EXM + b] - 273.15;
  let viscosity = 1.326e-5 * (1.0 + 6.542e-3 * celsius + 8.301e-6 * celsius * celsius - 4.84e-9 * celsius * celsius * celsius);
  let exS = D[D_EXL + b];
  let airV = IN[S_TH + b] * (1.0 + VIRT * IN[S_Q + b] - IN[S_QC + b]);
  let airQ = IN[S_Q + b];
  let surfaceQ = select(qsat(skin, pi), select(airQ, airQ + wetness * max(0.0, qsat(skin, pi) - airQ), X_WET_SURFACE), onLand);
  let surfaceV = skin / exS * (1.0 + VIRT * surfaceQ);
  var ri = 0.0; var c = vec2<f32>(0.0, 0.0);
  var passes = 1;
  if (X_GUSTY) { passes = X_GUST_ITER + 1; }
  for (var sweep = 0; sweep < passes; sweep++) {
    if (sweep > 0) {
      let buoyancy = -ri * speed * speed * speed * c.y / z;
      var gust = X_GUST_FLOOR;
      if (buoyancy > 0.0) { gust = scale * pow(buoyancy * mixed, 1.0 / 3.0); }
      speed = sqrt(wind * wind + gust * gust);
    }
    var blend = XBlend(0.0, 0.0);
    let snowUstar = XKARMAN * speed / log((z + Z0_SNOW) / Z0_SNOW);
    let snowH = xAndreas(Z0_SNOW, snowUstar, viscosity);
    if (onLand && onIceSheet) { xAdd(&blend, 1.0, Z0_SNOW, snowH); }
    else if (onLand) {
      let t = select(0.0, min(1.0, trees), VEGETATED); let gr = select(1.0, max(0.0, cover - t), VEGETATED); let bare = max(0.0, 1.0 - t - gr);
      let covered = min(1.0, snow / X_SNOWCOVER);
      xAdd(&blend, t, Z0_FOREST, ZH_FOREST); xAdd(&blend, gr * (1.0 - covered), Z0_GRASS, ZH_GRASS); xAdd(&blend, bare * (1.0 - covered), Z0_BARE, ZH_BARE);
      xAdd(&blend, (gr + bare) * covered, Z0_SNOW, snowH);
    } else {
      var z0 = 1.0e-4; var ustar = 0.0;
      for (var n = 0; n < 4; n++) {
        ustar = XKARMAN * speed / log((z + z0) / z0);
        let wind10 = ustar / XKARMAN * log((10.0 + z0) / z0);
        let alpha = max(0.0, CH_SLOPE * min(wind10, CH_CEIL) + CH_OFFSET);
        z0 = alpha * ustar * ustar / XGRAV + SMOOTH * viscosity / ustar;
      }
      ustar = XKARMAN * speed / log((z + z0) / z0);
      xAdd(&blend, 1.0 - concentration, z0, min(1.6e-4, 5.8e-5 / pow(z0 * ustar / viscosity, 0.72)));
      if (concentration > 0.0) {
        let zi = max(1.0e-3, 0.93e-3 * (1.0 - concentration) + 6.05e-3 * exp(-17.0 * (concentration - 0.5) * (concentration - 0.5)));
        xAdd(&blend, concentration, zi, xAndreas(zi, XKARMAN * speed / log((z + zi) / zi), viscosity));
      }
    }
    let lm = XKARMAN / sqrt(blend.drag);
    let z0m = X_BLEND * exp(-lm); let z0h = X_BLEND * exp(-XKARMAN * XKARMAN / (blend.heat * lm));
    ri = GRAV * z * (airV - surfaceV) / (0.5 * (airV + surfaceV) * speed * speed);
    c = xTransfer(ri, z, z0m, z0h);
  }
  PH[PH_XWIND + i] = speed;
  let d = 2.0 / 3.0 * REF_CROP; let rm = 0.123 * REF_CROP; let rh = 0.1 * rm;
  return vec3<f32>(c.x, c.y, REF_KARMAN * REF_KARMAN / (log((z - d) / rm) * log((z - d) / rh)));
}
`;
