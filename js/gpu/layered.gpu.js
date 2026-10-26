import { LAYERED_TABLE } from '../physics/layeredTable.module.js';
import { NEAR_INFRARED_TABLE } from '../physics/nearInfraredTable.module.js';
import { NEAR_INFRARED_MOLAR, OXYGEN, OZONE_SHARES, OZONE_COEFFICIENTS, VISIBLE_VAPOR } from '../physics/shortwaveGases.module.js';
import { SLINGO_LIQUID, EBERT_CURRY_ICE, NEAR_INFRARED_SHARE, ULTRAVIOLET_BANDS, VISIBLE_RAYLEIGH_SCALE, VISIBLE_OZONE_ABSORPTION } from '../physics/layeredShortwave.module.js';
import { REFERENCE_PRESSURE, NEAR_INFRARED_RAYLEIGH } from '../physics/radiation.module.js';

/*
 * The layered shortwave (js/physics/layeredShortwave.module.js) as a pass between the
 * radiation kernel and the longwave kernel. The radiation kernel, given radiation option
 * shortwave 'layered', takes the clear sky for the sunlight of the cells where there is
 * cloud outside the deck and leaves, for each such cell, the slots of PH_LYR: whether the
 * pass applies, the sun's cosine, the beam, the surface's albedos, the cover of the
 * cloud and the relative variance of its depth (0 for one node), the deck's fraction
 * and the aerosol's depth. This pass solves each slice's
 * layers (a clear region and each node of the cloud), adds the cloud's difference from
 * the clear region, weighted by the cover outside the deck, to the stashed sunlight
 * (the layers' heating, the surface's absorption and downward fluxes, the top's
 * absorption and reflection) and leaves the rest as the radiation kernel stashed it.
 * The near infrared's g-points are the table's, each with its band's cloud optics.
 */
export const LAYERED_SLOTS = 9;

const x = (v) => { if (!Number.isFinite(v)) throw new Error(`the layered table holds ${v}`); return v.toPrecision(9); };

function slices(o) {
  const mixing = { co2: o.carbonDioxide * NEAR_INFRARED_MOLAR.co2, ch4: o.methane * NEAR_INFRARED_MOLAR.ch4, o2: OXYGEN.mixingRatio * NEAR_INFRARED_MOLAR.o2 };
  const rayleigh = o.rayleighDepth != null ? [[1, o.rayleighDepth]] : o.rayleighBands;
  const out = [];
  for (const band of LAYERED_TABLE.bands) {
    const source = NEAR_INFRARED_TABLE.bands.find((b) => b.band === band.band);
    const parameters = band.gases.map((name) => source.gases[name]);
    const h2o = parameters[0];
    band.points.forEach((point) => {
      const terms = [];
      band.gases.forEach((name, i) => {
        const par = parameters[i], air = name === 'h2o' ? 'airV[k]' : `${x(mixing[name])} * airD[k]`;
        point.mix[i].forEach((mix, j) => { if (mix > 0) terms.push(`d += ${x(mix * par.k[j])} * ${air} * exp(${x(par.n[j])} * lv${i} + ${x(par.t[j])} * airDT[k]);`); });
      });
      const logs = band.gases.map((name, i) => `let lv${i} = log(sqrt(airP2[k] + ${x(parameters[i].pD * parameters[i].pD)}) / ${x(NEAR_INFRARED_TABLE.pRef)});`).join(' ');
      const self = h2o.self > 0 ? `var d = ${x(h2o.self)} * airE[k] * exp(${x(h2o.selfT)} * airLn[k]);` : 'var d = 0.0;';
      out.push({
        text: `w = ${x(band.share * point.weight)}; sb = ${band.slingo}; ray = ${x(band.rayleigh * o.nearInfraredRayleigh / NEAR_INFRARED_RAYLEIGH)}; visible = false;
        for (var k = 0; k < K; k++) { ${logs} ${self} ${terms.join(' ')} gas[k] = ${x(point.scale)} * d; }`,
      });
    });
  }
  rayleigh.forEach(([weight, depth]) => out.push({ text: `w = visShare * ${x(weight)}; sb = 0; ray = ${x(depth * VISIBLE_RAYLEIGH_SCALE)}; visible = true;
        for (var k = 0; k < K; k++) { gas[k] = vis[k]; }` }));
  return out;
}

const chain = (list) => list.map((s, n) => `${n ? ' else ' : ''}if (s == ${n}) {\n        ${s.text}\n      }`).join('');

const COEFFICIENTS = (name, table, field, pick) => `if (sb == 0) { ${name} = ${x(table[field][0][pick])}; } else if (sb == 1) { ${name} = ${x(table[field][1][pick])}; } else if (sb == 2) { ${name} = ${x(table[field][2][pick])}; } else { ${name} = ${x(table[field][3][pick])}; }`;

export const layeredKernel = (o, held, stashText) => {
  const list = slices(o);
  const ultraviolet = OZONE_SHARES.slice(0, ULTRAVIOLET_BANDS).map((share, b) => `${x(share)} * (1.0 - exp(-${x(OZONE_COEFFICIENTS[b])} * ozonePath))`).join(' + ');
  const aerosolAlbedo = o.aerosolAlbedo, aerosolAsymmetry = o.aerosolAsymmetry;
  const outputs = held ? `    PH[PH_RADSW + k * C + i] += weightL * dh[k];` : `    OUT[LWS_FLUX + k * C + i] += weightL * dh[k];`;
  const scalars = held
    ? `  let dAbsKept = max(dAbs, -PH[PH_RADABS + i]);
  PH[PH_RADABS + i] += dAbsKept; PH[PH_RADSWDN + i] += max(weightL * down, -PH[PH_RADSWDN + i]); PH[PH_RADDIR + i] += max(weightL * dir, -PH[PH_RADDIR + i]);
  PH[PH_RADATM + i] += dAtm; PH[PH_RADTOA + i] = PH[PH_RADABS + i] + PH[PH_RADATM + i]; PH[PH_RADREFL + i] = beam - PH[PH_RADTOA + i];`
    : `  let dAbsKept = max(dAbs, -stashed(14, i)); let atmosphereNow = stashed(17, i) + dAtm;
  stash(14, i, stashed(14, i) + dAbsKept); stash(15, i, stashed(15, i) + max(weightL * down, -stashed(15, i))); stash(16, i, stashed(16, i) + max(weightL * dir, -stashed(16, i)));
  stash(17, i, atmosphereNow); stash(18, i, stashed(18, i) + dAbsKept + dAtm); stash(19, i, beam - stashed(18, i));`;
  return `${stashText}fn cumulusCloud(k: i32, i: i32) -> vec2<f32> {
  if (!CU_CLOUD || k < CU_K0) { return vec2<f32>(0.0, 0.0); }
  let slot = (k - CU_K0) * C + i;
  return vec2<f32>(PH[PH_CUCOVER + slot], PH[PH_CUCOVER + slot] * PH[PH_CUWATER + slot]);
}
struct Layer { rd: f32, td: f32, rb: f32, tb: f32, tt: f32 }
fn expRatio(x: f32) -> f32 {
  if (x < 0.25) { return 1.0 - x * 0.5 * (1.0 - x / 3.0 * (1.0 - x * 0.25 * (1.0 - x * 0.2 * (1.0 - x / 6.0 * (1.0 - x / 7.0))))); }
  return (1.0 - exp(-x)) / x;
}
fn twoStream(tau: f32, co: f32, g: f32, mu0: f32) -> Layer {
  let f = g * g; let q = 1.0 - f + f * co; let scaled = q * tau; let cs = co / q; let gs = g / (1.0 + g); let w = 1.0 - cs;
  let g1 = (3.0 - 3.0 * gs + cs * (4.0 + 3.0 * gs)) * 0.25; let g2 = max(0.0, (3.0 - 3.0 * gs - cs * (4.0 - 3.0 * gs)) * 0.25);
  let k = sqrt(select(g1 * g1, 2.0 * cs * (g1 + g2), g2 > 0.0));
  let u = k * scaled; let e = exp(-u); let big = scaled * expRatio(2.0 * u); let den = 1.0 + (g1 - k) * big;
  var result: Layer;
  result.rd = g2 * big / den; result.td = e / den;
  if (w < 1e-6) { result.rb = 0.0; result.tb = 0.0; result.tt = exp(-scaled / mu0); return result; }
  let g3 = (2.0 - 3.0 * gs * mu0) * 0.25; let g4 = 1.0 - g3; let a1 = g1 * g4 + g2 * g3; let a2 = g1 * g3 + g2 * g4; let km = k * mu0;
  let v = scaled / mu0; let t = exp(-v); let d = scaled * (1.0 - km) / mu0; let ratio = expRatio(abs(d));
  let fa = select(e * t, e * e, d >= 0.0); let fb = select(t, e, d >= 0.0); let pref = w / ((1.0 + km) * den);
  result.rb = pref * ((a2 + k * g3) * big - (a2 * mu0 - g3) * v * fa * ratio);
  result.tb = pref * ((g4 + a1 * mu0) * v * fb * ratio - t * (a1 - k * g4) * big);
  result.tt = t;
  return result;
}
@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= C) { return; }
  if (!(stashed(0, i) > 0.0) || !(PH[PH_LYR + i] > 0.5)) { return; }
  let mu = PH[PH_LYR + C + i]; let beam = PH[PH_LYR + 2 * C + i]; let adir = PH[PH_LYR + 3 * C + i]; let adif = PH[PH_LYR + 4 * C + i];
  let cover = PH[PH_LYR + 5 * C + i]; let variance = PH[PH_LYR + 6 * C + i]; let fraction = PH[PH_LYR + 7 * C + i]; let aerosol = PH[PH_LYR + 8 * C + i];
  let pi = IN[S_PI + i];
  let muEff = sqrt(1224.0 * mu * mu + 1.0) / 35.0; let pressureScale = pi / ${x(REFERENCE_PRESSURE)};
  let continental = PH[PH_LAND + i] > 0.5 && !(PH[PH_LAND + i] > 1.5);
  let liquidRadius = select(DROP_SEA, DROP_LAND, continental);
  var airV: array<f32, K>; var airD: array<f32, K>; var airP2: array<f32, K>; var airDT: array<f32, K>; var airLn: array<f32, K>; var airE: array<f32, K>;
  var vis: array<f32, K>; var wat: array<f32, K>; var gas: array<f32, K>;
  var rd: array<f32, K>; var td: array<f32, K>; var rb: array<f32, K>; var tb: array<f32, K>; var tt: array<f32, K>;
  var alb: array<f32, K + 1>; var src: array<f32, K + 1>; var dh: array<f32, K>;
  var ozoneSum = 0.0;
  for (var k = 0; k < K; k++) {
    let a = nirAir(i, k, pi);
    airV[k] = a.vapour; airD[k] = a.dry; airP2[k] = a.p2; airDT[k] = a.dT; airLn[k] = a.lnT; airE[k] = a.vapourE;
    let idx = k * C + i; let mass = pi * LV[L_DS + k] / GRAV; let q = max(0.0, IN[S_Q + idx]); let p = pi * LV[L_SM + k];
    vis[k] = q * mass * 0.1 * pow(p / SCALE_P, SCALE_N) * (1.0 + 0.00135 * (a.dT + 10.0)); ozoneSum += OUT[LWS_OZONE + idx];
    wat[k] = max(0.0, IN[S_QC + idx]) * mass + cumulusCloud(k, i).y * mass;
  }
  let ozonePath = ozoneSum / muEff; let ultraviolet = ${ultraviolet};
  let visShare = 1.0 - ${x(NEAR_INFRARED_SHARE)} - ultraviolet; let ozoneStrength = ${x(VISIBLE_OZONE_ABSORPTION)} / visShare; let vaporStrength = ${x(VISIBLE_VAPOR.share)} * VIS_H2O_K * VIS_STRENGTH / visShare;
  for (var k = 0; k < K; k++) { vis[k] = ozoneStrength * OUT[LWS_OZONE + k * C + i] + vaporStrength * vis[k]; }
  let nodes = select(1, 4, variance > 0.0);
  let sigma = sqrt(log(1.0 + max(variance, 0.0))); let shift = 0.5 * sigma * sigma;
  var surf = 0.0; var down = 0.0; var dir = 0.0;
  for (var s = 0; s < ${list.length}; s++) {
    var w = 0.0; var sb = 0; var ray = 0.0; var visible = false;
    ${chain(list)}
    var la = 0.0; var lb = 0.0; var lc = 0.0; var ld = 0.0; var le = 0.0; var lf = 0.0; var ic = 0.0; var id2 = 0.0; var ie = 0.0; var if2 = 0.0; var la0 = 0.0; var lb0 = 0.0;
    ${COEFFICIENTS('la', SLINGO_LIQUID, 'extinction', 0)}
    ${COEFFICIENTS('lb', SLINGO_LIQUID, 'extinction', 1)}
    ${COEFFICIENTS('lc', SLINGO_LIQUID, 'absorption', 0)}
    ${COEFFICIENTS('ld', SLINGO_LIQUID, 'absorption', 1)}
    ${COEFFICIENTS('le', SLINGO_LIQUID, 'asymmetry', 0)}
    ${COEFFICIENTS('lf', SLINGO_LIQUID, 'asymmetry', 1)}
    ${COEFFICIENTS('ic', EBERT_CURRY_ICE, 'absorption', 0)}
    ${COEFFICIENTS('id2', EBERT_CURRY_ICE, 'absorption', 1)}
    ${COEFFICIENTS('ie', EBERT_CURRY_ICE, 'asymmetry', 0)}
    ${COEFFICIENTS('if2', EBERT_CURRY_ICE, 'asymmetry', 1)}
    la0 = ${x(SLINGO_LIQUID.extinction[0][0])}; lb0 = ${x(SLINGO_LIQUID.extinction[0][1])};
    let extinctionRatio = (la + lb / liquidRadius) / (la0 + lb0 / liquidRadius);
    for (var region = 0; region <= nodes; region++) {
      var node = 1.0; var nodeWeight = 1.0;
      if (region > 0 && nodes == 4) {
        let n = region - 1;
        let z = select(select(-0.7419637843027259, 0.7419637843027259, n == 2), select(-2.3344142183389778, 2.3344142183389778, n == 3), n == 0 || n == 3);
        node = exp(sigma * z - shift); nodeWeight = select(0.45412414523193156, 0.0458758547680685, n == 0 || n == 3);
      }
      let polarity = select(1.0, -1.0, region == 0); let wt = select(w * nodeWeight, w, region == 0);
      for (var k = 0; k < K; k++) {
        let cw = select(0.0, wat[k] / cover * node, region > 0);
        if (region > 0 && !(cw > 0.0)) { continue; }
        let rayleigh = ray * pressureScale * LV[L_DS + k];
        let aer = select(0.0, aerosol * LV[L_AER + k], visible);
        var cloudTau = 0.0; var cloudAbsorbed = 0.0; var cloudScattered = 0.0; var cloudAsymmetry = 0.0;
        if (cw > 0.0) {
          let T = airDT[k] + 250.0;
          let liquid = clamp((T - ICE_T) / (LIQUID_T - ICE_T), 0.0, 1.0);
          let c = clamp(T - 273.15, ICE_COLDEST, ICE_WARMEST);
          let crystal = 0.5 * (326.3 + c * (12.42 + c * (0.197 + c * 0.0012)));
          let liquidTau = liquid * cw * 1500.0 / liquidRadius * extinctionRatio; let iceTau = (1.0 - liquid) * cw * 1000.0 * (3.448e-3 + 2.431 / crystal);
          let liquidCo = max(0.0, lc + ld * liquidRadius); let iceCo = ic + id2 * crystal;
          cloudTau = liquidTau + iceTau; cloudAbsorbed = liquidTau * liquidCo + iceTau * iceCo; cloudScattered = cloudTau - cloudAbsorbed;
          cloudAsymmetry = liquidTau * (1.0 - liquidCo) * (le + lf * liquidRadius) + iceTau * (1.0 - iceCo) * (ie + if2 * crystal);
        }
        let tau = gas[k] + rayleigh + aer + cloudTau;
        let absorbed = gas[k] + ${x(1 - aerosolAlbedo)} * aer + cloudAbsorbed; let scattered = rayleigh + ${x(aerosolAlbedo)} * aer + cloudScattered;
        let asymmetry = select(0.0, (${x(aerosolAsymmetry * aerosolAlbedo)} * aer + cloudAsymmetry) / scattered, scattered > 0.0);
        let layer = twoStream(tau, absorbed / tau, asymmetry, muEff);
        rd[k] = layer.rd; td[k] = layer.td; rb[k] = layer.rb; tb[k] = layer.tb; tt[k] = layer.tt;
      }
      var up = adif; var source = adir;
      alb[K] = up; src[K] = source;
      for (var j = 0; j < K; j++) {
        let k = K - 1 - j;
        let rest = 1.0 / (1.0 - rd[k] * up);
        let nextUp = rd[k] + td[k] * td[k] * up * rest;
        let nextSource = rb[k] + td[k] * (up * tb[k] + source * tt[k]) * rest;
        up = nextUp; source = nextSource; alb[k] = up; src[k] = source;
      }
      var direct = 1.0; var diffuse = 0.0; var net = 1.0 - source;
      for (var k = 0; k < K; k++) {
        let next = direct * tt[k]; let rest = 1.0 / (1.0 - rd[k] * alb[k + 1]);
        let below = (tb[k] * direct + td[k] * diffuse + rd[k] * src[k + 1] * next) * rest;
        let through = next + below - (alb[k + 1] * below + src[k + 1] * next);
        dh[k] += polarity * wt * (net - through);
        net = through; direct = next; diffuse = below;
      }
      surf += polarity * wt * net; down += polarity * wt * (direct + diffuse); dir += polarity * wt * direct;
    }
  }
  let weightL = beam * (1.0 - fraction) * cover;
  var dAtm = 0.0;
  for (var k = 0; k < K; k++) {
${outputs}
    dAtm += weightL * dh[k];
  }
  let dAbs = weightL * surf;
${scalars}
}
`;
};
