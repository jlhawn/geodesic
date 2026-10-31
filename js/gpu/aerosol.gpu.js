import { AEROSOL_SPECTRA } from '../physics/aerosolSpectra.module.js';
import { AEROSOL_FIELD, AEROSOL_FIELDS, WAVELENGTHS, MIXED_DEPTH, MIXED_CAP, HEIGHT_RANGE, ASYMMETRY_CEILING } from '../aerosol.module.js';

/*
 * Radiation option aerosol 'mac' on the GPU (js/aerosol.module.js has the arithmetic, which this
 * repeats in single precision): the climatology's eight fields of a cell are the slots
 * PH_AER + field * C + i of the physics buffer, written once a model day by the core, and these
 * functions make the optics of a cell's two modes from them. Without the option none of this
 * text is generated. The functions here belong to the constants of physicsConstants; the
 * profile of the modes over the layers is code for the kernels that need it (macProfiles),
 * which keep its two arrays, macF and macG, to themselves.
 */
const x = (v) => { if (!Number.isFinite(v)) throw new Error(`the aerosol table holds ${v}`); return v.toPrecision(9); };
const table = (name, values) => `const ${name}: array<f32, ${values.length}> = array<f32, ${values.length}>(${values.map(x).join(', ')});`;
const { fine, coarse } = AEROSOL_SPECTRA;

/*
 * Shader text with the light of a column a vec4, whose fourth number is the near infrared
 * aerosol's scattering depth (the physics functions have the three of the visible share, the
 * pressure ratio and the visible aerosol's); `near` is whether the text adds the near
 * infrared's Rayleigh depth, to which this one is added.
 */
export function aerosolLight(text, near = true) {
  const edits = [['light: vec3<f32>', 'light: vec4<f32>'], ...(near ? [['NIR_RAY * light.y', 'NIR_RAY * light.y + light.w']] : [])];
  return edits.reduce((out, [from, to]) => {
    if (!out.includes(from)) throw new Error(`the shader text has no ${from} for the aerosol's light`);
    return out.replaceAll(from, to);
  }, text);
}

export function aerosolConstants(o) {
  if (!(Number.isFinite(o.aerosolScale) && o.aerosolScale >= 0)) throw new Error(`aerosolScale must be a factor of the aerosol's depth, 0 or more, not ${o.aerosolScale}`);
  const T = (name, key) => table(name, [...fine[key], ...coarse[key]]);
  const stream = (m) => `{
    let tau = MAC_SCALE * macField(${m}, i); let x = macField(${2 + m}, i); let g = macField(${4 + m}, i);
    for (var b = 0; b < ${WAVELENGTHS}; b++) {
      let one = macBand(${m}, b, tau, x, g); let scattering = one.x * (1.0 - one.y);
      let wv = MAC_VIS[b]; let wn = MAC_NIR[b];
      if (wv > 0.0) { tV += wv * one.x; sV += wv * scattering; pV += wv * scattering * one.z; aV${m} += wv * one.x * one.y; dV${m} += wv * one.x; }
      if (wn > 0.0) { tN += wn * one.x; sN += wn * scattering; pN += wn * scattering * one.z; aN${m} += wn * one.x * one.y; dN${m} += wn * one.x; }
    }
  }`;
  return `
const MAC_SCALE: f32 = ${x(o.aerosolScale)}; const MAC_HMIN: f32 = ${x(HEIGHT_RANGE[0])}; const MAC_HMAX: f32 = ${x(HEIGHT_RANGE[1])}; const MAC_G_MAX: f32 = ${x(ASYMMETRY_CEILING)};
const MAC_MIX0: f32 = ${x(MIXED_DEPTH[0])}; const MAC_MIX1: f32 = ${x(MIXED_DEPTH[1])}; const MAC_MIX_CAP: f32 = ${x(MIXED_CAP)};
${T('MAC_R', 'r')}
${T('MAC_A', 'a')}
${T('MAC_C', 'c')}
${T('MAC_D', 'd')}
${T('MAC_E', 'e')}
${table('MAC_VIS', AEROSOL_SPECTRA.visible)}
${table('MAC_NIR', AEROSOL_SPECTRA.nearInfrared)}
fn macField(f: i32, i: i32) -> f32 { return PH[PH_AER + f * C + i]; }
fn macHeight(k: i32, i: i32) -> f32 { return (D[D_GEO + k * C + i] + LV[L_GABS + k]) / (1000.0 * GRAV); }
fn macBand(mode: i32, b: i32, tau: f32, x: f32, g: f32) -> vec3<f32> {
  let n = mode * ${WAVELENGTHS} + b;
  return vec3<f32>(MAC_R[n] * tau, clamp(MAC_A[n] + MAC_C[n] * x, 0.0, 1.0), clamp(MAC_D[n] + MAC_E[n] * g, 0.0, MAC_G_MAX));
}
struct MacStreams { tauV: f32, omegaV: f32, gV: f32, fineV: f32, tauN: f32, omegaN: f32, gN: f32, fineN: f32 }
fn macStreams(i: i32) -> MacStreams {
  var tV = 0.0; var sV = 0.0; var pV = 0.0; var tN = 0.0; var sN = 0.0; var pN = 0.0;
  var aV0 = 0.0; var aV1 = 0.0; var dV0 = 0.0; var dV1 = 0.0; var aN0 = 0.0; var aN1 = 0.0; var dN0 = 0.0; var dN1 = 0.0;
  ${stream(0)}
  ${stream(1)}
  let fineV = select(dV0, aV0, aV0 + aV1 > 0.0); let otherV = select(dV1, aV1, aV0 + aV1 > 0.0);
  let fineN = select(dN0, aN0, aN0 + aN1 > 0.0); let otherN = select(dN1, aN1, aN0 + aN1 > 0.0);
  var out: MacStreams;
  out.tauV = tV; out.omegaV = select(1.0, sV / tV, tV > 0.0); out.gV = select(0.0, pV / sV, sV > 0.0);
  out.tauN = tN; out.omegaN = select(1.0, sN / tN, tN > 0.0); out.gN = select(0.0, pN / sN, sN > 0.0);
  out.fineV = select(0.0, fineV / (fineV + otherV), fineV + otherV > 0.0);
  out.fineN = select(0.0, fineN / (fineN + otherN), fineN + otherN > 0.0);
  return out;
}
struct MacLight { taken: f32, takenN: f32, light: vec4<f32> }
fn macLight(beam: f32, mu: f32, ozoneHeating: f32, incident: f32, pi: f32, seq: MacStreams, nirOn: bool, gasAbove: f32) -> MacLight {
  let visible = max(0.0, VIS_FRAC * beam - ozoneHeating);
  let magnification = 35.0 / sqrt(1224.0 * mu * mu + 1.0);
  var taken = 0.0; var takenN = 0.0; var scatterN = 0.0;
  if (mu > 0.0 && seq.tauV > 0.0) { taken = visible * (1.0 - exp(-(1.0 - seq.omegaV) * seq.tauV * magnification)); }
  var left = incident - taken;
  if (nirOn && seq.tauN > 0.0) {
    takenN = max(0.0, beam * (1.0 - VIS_FRAC - gasAbove)) * (1.0 - exp(-(1.0 - seq.omegaN) * seq.tauN * magnification));
    left = left - takenN;
    scatterN = (1.0 - seq.gN) * seq.omegaN * seq.tauN;
  }
  var out: MacLight;
  out.taken = taken; out.takenN = takenN;
  out.light = vec4<f32>(select(0.0, min(1.0, (visible - taken) / left), left > 0.0), pi, (1.0 - seq.gV) * seq.omegaV * seq.tauV, scatterN);
  return out;
}
`;
}

/*
 * Statements of a kernel's main for the cell i: macF and macG, the fine and the coarse mode's
 * share of its column in each layer, filled where `guard` holds, counting up and read in place.
 */
export const macProfiles = (guard = 'true') => `var macF: array<f32, K>; var macG: array<f32, K>;
  if (${guard}) {
${[0, 1].map((m) => `  {
    let mean = clamp(macField(${AEROSOL_FIELD.fineHeight + m}, i), MAC_HMIN, MAC_HMAX); let mixed = min(MAC_MIX${m}, MAC_MIX_CAP * mean);
    let scale = 0.5 * ((mean - mixed) + sqrt(mean * mean + 2.0 * mean * mixed - mixed * mixed)); let total = mixed + scale;
    var upper = 0.0; var zHere = macHeight(0, i);
    for (var k = 0; k < K; k++) {
      var height = 0.0;
      if (k < K - 1) { let zNext = macHeight(k + 1, i); height = 0.5 * (zHere + zNext); zHere = zNext; }
      let lower = select(1.0 - height / total, scale / total * exp(-(height - mixed) / scale), height >= mixed);
      mac${m === 0 ? 'F' : 'G'}[k] = lower - upper; upper = lower;
    }
  }`).join('\n')}
  }
`;

/*
 * The near infrared gas the beam has lost above the aerosol, in `gasAbove` (0 unless `guard`
 * holds): the cumulative absorption of the near infrared gases in OUT at `nirAt`, at the
 * layers' mid-levels, weighted by the aerosol's share of the layers in the near infrared, of
 * which `fine` is the fine mode's.
 */
export const macGasAbove = (fine, nirAt, guard) => `var gasAbove = 0.0;
    if (${guard}) {
      var nirBefore = 0.0;
      for (var k = 0; k < K; k++) { let nirNow = OUT[${nirAt} + k * C + i]; gasAbove += (${fine} * macF[k] + (1.0 - ${fine}) * macG[k]) * 0.5 * (nirBefore + nirNow); nirBefore = nirNow; }
    }`;

export const AEROSOL_SLOTS = AEROSOL_FIELDS;

/*
 * The layered kernel's per-slice statements: mq, ma, ms and mp (0 the fine mode, 1 the coarse)
 * the slice's depth, absorption depth, scattering depth and asymmetry-weighted scattering depth
 * of each mode, of the wavelength index `index` alone or, for the visible slices, of the
 * wavelengths of AEROSOL_SPECTRA.slices. The cell's 550 nm numbers are mt, mx and mg.
 */
export const macCell = () => `let mt0 = MAC_SCALE * macField(0, i); let mx0 = macField(2, i); let mg0 = macField(4, i);
  let mt1 = MAC_SCALE * macField(1, i); let mx1 = macField(3, i); let mg1 = macField(5, i);
  var mq0 = 0.0; var ma0 = 0.0; var ms0 = 0.0; var mp0 = 0.0; var mq1 = 0.0; var ma1 = 0.0; var ms1 = 0.0; var mp1 = 0.0;
`;
const band = (m, index) => `{ let one = macBand(${m}, ${index}, mt${m}, mx${m}, mg${m}); mq${m} = one.x; ma${m} = one.x * one.y; ms${m} = one.x * (1.0 - one.y); mp${m} = one.x * (1.0 - one.y) * one.z; }`;
export const macSliceBand = (index) => ` ${band(0, index)} ${band(1, index)}`;
export const macSliceVisible = () => {
  const weights = AEROSOL_SPECTRA.slices.map((w, b) => [w, b]).filter(([w]) => w > 0);
  const mode = (m) => `{ ${weights.map(([w, b], n) => `let one${n} = macBand(${m}, ${b}, mt${m}, mx${m}, mg${m}); let t${n} = ${x(w)} * one${n}.x;`).join(' ')}
    mq${m} = ${weights.map((_, n) => `t${n}`).join(' + ')}; ma${m} = ${weights.map((_, n) => `t${n} * one${n}.y`).join(' + ')};
    ms${m} = ${weights.map((_, n) => `t${n} * (1.0 - one${n}.y)`).join(' + ')}; mp${m} = ${weights.map((_, n) => `t${n} * (1.0 - one${n}.y) * one${n}.z`).join(' + ')}; }`;
  return ` ${mode(0)} ${mode(1)}`;
};
