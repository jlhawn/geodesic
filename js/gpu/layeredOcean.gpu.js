import { emptyBuffer, readRanges, reductionKernel, finishReduction, reductionGroups, writeInPieces } from './device.module.js';
import { LAKE_DRAIN_WGSL } from './lake.gpu.js';
import { oceanCells, restingColumns, displacedWeight, beyondCells, SHEAR_DIFFUSIVITY, HYCOM_DIFFUSION_SPEED, TRACER_ADVECTIONS, EDGE_TRANSPORTS, LAYER_DENSITIES, LAYER_SALINITIES, THERMOCLINE_DENSITY, EPS, THIN, PV_FLOOR, SPEED_LIMIT, CLOSURE_SPACING, EDDY_BOTTOM_TAPER, EDDY_SLACK, closureCoefficient, inertialStepLimit, eddyDiffusivities, eddyDiffusionLimit, bathymetryFrom, edgeSills, runoffOutlets, abyssalCells, savedDensities, sameDensities, JERLOV, PENETRATION_DEPTH, FRICTION_PER_WIND, CALM_WIND } from '../ocean/layered.module.js';
import { SEAWATER, SEAWATER_WGSL, seawaterDensity, labelTemperature } from '../ocean/seawater.module.js';
import { HYBRID_LAYERS, H_MIN, MODE_DEPTH, GENERATORS, ROTATION_LIMITS, RH18_NEUTRAL, LANGMUIR_ICE, ICE_WAVE_DAMPING, labelBoundaries, columnScratch, spaceInterfaces, remapColumn, mixedLayerDepth, transmitted as transmittedShare } from '../ocean/column.module.js';
import { profileAt } from '../ocean/climatology.module.js';
import { FREEZING_POINT, MINIMUM_CONCENTRATION, SNOW_AGEING, seaIceRemoval } from '../physics/ice.module.js';
import { iceDriftOptions, iceDriftPasses, ICE_REMNANT } from '../ocean/iceDrift.module.js';
import { iceDriftHead, iceDriftKernels, driftStressKernel, driftSections } from './iceDrift.gpu.js';
import { basalMeltWgsl, baseFreezing, SHELF_MELT } from '../ocean/iceShelfMelt.module.js';

const waveFactor = (o) => (o.langmuirIce === 'tanh' ? 'select(1.0, 1.0 + (langmuirFactor(ustar, max(depth, 1.0)) - 1.0) * waves, ustar > 0.0)' : 'select(1.0, langmuirFactor(ustar, max(depth, 1.0)), ustar > 0.0 && waves)');

/*
 * The ocean of ocean/layered.module.js on the GPU: the ALE hybrid layers,
 * each a TRiSK shallow-water layer carrying thickness, edge velocity,
 * heat h·(T − HEAT_REFERENCE) and salt h·(S − SALT_REFERENCE), each
 * measured from its reference so that a thick layer's float32 store
 * keeps the eddy transport's small changes, and the sunlight it deposits
 * below the top layer (`penetration`), whose remainder under each store's
 * spacing waits in PH_PEN. The RK4 baroclinic state (h, u, heat, salt
 * for all L layers) lives in a ping-pong buffer set (S/T/K1-K4);
 * everything else the step needs (edge thicknesses, PV, the pressure
 * pairs, the ∇⁴ closure scratch, the barotropic sub-stepping state, the
 * regrid targets, the viscosities, the surface staging) lives in one big
 * scratch/static buffer (OD) bound alongside it. The barotropic mode
 * takes its own small RK4 with fixed-offset stage blocks inside OD,
 * looped on the host the M times a step needs; every other multi-stage
 * loop (the baroclinic RK4, the ∇⁴ closure's two Laplacian passes) mirrors
 * the dispatch-ordering-is-execution-ordering pattern the sigma core
 * relies on within one compute pass.
 *
 * Initialization, loading and serialization run once per model build and
 * are cheap relative to a step, so they stay in JavaScript at double
 * precision, exactly porting the CPU functions over plain arrays before
 * uploading; only the per-step work (tendency, barotropic, eddy
 * transport, regrid, column physics, salt, surface write-back) is WGSL.
 */
const WORKGROUP = 64;
// Float32 values, so the host's conversions use the numbers the kernels do.
export const HEAT_REFERENCE = Math.fround(273.15), SALT_REFERENCE = Math.fround(35);

export const OCEAN_DEFAULTS = {
  minimumDepth: 50, flatDepth: 4000, thermoclineTilt: 0.3,
  density: 1025, specificHeat: 3985, referenceS: 35, gravity: 9.81,
  vorticityCentring: 0.5, iceSalinity: 5, iceStressTransmission: 0.8, iceDensity: 917, latentHeatFusion: 3.34e5, langmuir: 2, langmuirIce: 'tanh',
  momentumDepth: 50, shearViscosity: 1e-2, backgroundViscosity: 1e-4, bottomDrag: 3e-3, topDrag: 2.5e-3, closureHours: 12, closureSpacing: CLOSURE_SPACING, everySteps: 4,
  eddyDiffusivity: 1000, eddyTaperDepth: 200, penetration: 'IB',
  gustiness: 3,
  generator: 'hybrid', regridTime: 0, minimumLayer: H_MIN, modeDepth: MODE_DEPTH, mstar: 1.25, ekmanFactor: 0.5, rotationLimit: 'rh18', nstar: 0.2, backgroundDiffusivity: 1e-5, convectiveDiffusivity: 0.1,
  tracerAdvection: 'third', shearDiffusivity: SHEAR_DIFFUSIVITY, shearDiffusivityCap: 4000, lateralDiffusionSpeed: HYCOM_DIFFUSION_SPEED,
  edgeTransport: 'kept', shelfMelt: true, tidalSpeed: SHELF_MELT.tidalSpeed, shelfHeatTransfer: SHELF_MELT.heatTransfer, shelfSaltTransfer: SHELF_MELT.saltTransfer, shelfMeltLayer: 20,
};
const defaultSalinityProfile = (lat) => 34 + 2 * Math.exp(-(((Math.abs(lat) * 180 / Math.PI - 25) / 20) ** 2));

export function oceanStateLength(C, E, options = {}) {
  const o = { ...OCEAN_DEFAULTS, ...options };
  return (o.layerTable ?? HYBRID_LAYERS).K * (3 * C + E);
}

function seq(names) { const out = {}; let off = 0; for (const [name, n] of names) { for (const alias of [].concat(name)) out[alias] = off; off += n; } out.total = off; return out; }

/*
 * The column kernels (docs/ocean-hybrid-coordinate.md; the CPU's
 * ocean/layered.module.js and the column library in
 * ocean/column.module.js). Every per-column array is a local of the kernel's
 * entry point of at most L + 1 floats, and the column routines (the
 * tridiagonal solve, the limited reconstruction, the remap, the boundary
 * layer's work) are written into each kernel by the generators below
 * rather than called as WGSL functions, since Tint hands a function its
 * arrays by pointer and Apple's compilers then lose writes to them (the
 * forms in CLAUDE.md). Every loop over an array counts up; the solve's
 * elimination carries the layer above's coefficient and gain in scalars
 * and its back substitution runs over j with k = n − 2 − j carrying the
 * value it just wrote, since Metal loses the read of an element written
 * in the previous iteration. The
 * remap gives each old layer's pieces to the new layers and its rounding
 * residual to the last of them, so a column's water, heat and salt are
 * kept to the stores' rounding. Heat and salt are the stores' anomalies
 * h·(T − TREF) and h·(S − SREF) throughout; only the equation of state
 * sees the absolute values.
 */
const declareArrays = (names, n = 'L') => names.map((name) => `var ${name}: array<f32, ${n}>;`).join(' ');
// solveLayers: the implicit diffusion of the first N layers of H with KAPPA at their interfaces, the one below layer N − 1 closed, from Q1 (and Q2) into D1 (and D2), an output may be its input: the elimination carries the layer above's coefficient and gain in scalars, the substitution the value it just wrote; GAIN holds the gains.
const solveInline = ({ h, kappa, n, dt, q1, d1, q2 = null, d2 = null, gain = 'SG' }, body = '') => `{
      var wPrev = 0.0; var gPrev = 0.0; var r1 = 0.0;${q2 ? ' var r2 = 0.0;' : ''}
      for (var k = 0; k < ${n}; k++) {
        let w = select(0.0, ${dt} * ${kappa}[k] / max(0.5 * (${h}[k] + ${h}[min(k + 1, L - 1)]), DZFLOOR), k < ${n} - 1);
        let m = ${h}[k] + wPrev + w + wPrev * gPrev; let g = -w / m;
        r1 = (${h}[k] * ${q1}[k] + wPrev * r1) / m; ${d1}[k] = r1;${q2 ? ` r2 = (${h}[k] * ${q2}[k] + wPrev * r2) / m; ${d2}[k] = r2;` : ''}
        ${gain}[k] = g; wPrev = w; gPrev = g;
      }
      var n1 = ${d1}[${n} - 1];${q2 ? ` var n2 = ${d2}[${n} - 1];` : ''}
      for (var j = 0; j < ${n} - 1; j++) { let k = ${n} - 2 - j; n1 = ${d1}[k] - ${gain}[k] * n1; ${d1}[k] = n1;${q2 ? ` n2 = ${d2}[k] - ${gain}[k] * n2; ${d2}[k] = n2;` : ''} }
${body}    }`;
/*
 * The regrid kernels sweep a column once, old layer by old layer, with
 * the layers an interface value needs held in scalars (the four around it
 * for the quartic, the three around a layer for the linear slope) and the
 * new layer being filled in one accumulator, so they keep no private
 * array. A sweep reads old layer j + 3 (j + 2 on an edge) after finishing
 * old layer j, so a new layer up to that index is written in place and a
 * later one into OUT, free after the dynamics, and copied back at the end.
 */
const ppmEdgeInline = ({ k, hm1, h0, h1, h2, tracers }) => `    ${tracers.map(({ q0, out }) => `${out} = ${q0};`).join(' ')}
    if (${k} < L) {
      if (${k} >= 2 && ${k} <= L - 2 && ${hm1} > HNEGLECT && ${h0} > HNEGLECT && ${h1} > HNEGLECT && ${h2} > HNEGLECT) {
        let x0 = -(${hm1} + ${h0}); let x1 = -${h0}; let x3 = ${h1}; let x4 = ${h1} + ${h2};
        let c1 = -(x0 * x3 * x4) / ((x1 - x0) * x1 * (x1 - x3) * (x1 - x4));
        let c2 = -(1.0 / x0 + 1.0 / x1 + 1.0 / x3 + 1.0 / x4);
        let c3 = -(x0 * x1 * x4) / ((x3 - x0) * (x3 - x1) * x3 * (x3 - x4));
        let c4 = -(x0 * x1 * x3) / ((x4 - x0) * (x4 - x1) * (x4 - x3) * x4);
${tracers.map(({ qm1, q0, q1, q2, out }) => `        { let p1 = ${hm1} * ${qm1}; let p2 = p1 + ${h0} * ${q0}; let p3 = p2 + ${h1} * ${q1}; let p4 = p3 + ${h2} * ${q2}; ${out} = c1 * p1 + c2 * p2 + c3 * p3 + c4 * p4; }`).join('\n')}
      } else {
        let d = max(${h0} + ${h1}, HNEGLECT);
        ${tracers.map(({ q0, q1, out }) => `${out} = (${h1} * ${q0} + ${h0} * ${q1}) / d;`).join(' ')}
      }
      ${tracers.map(({ q0, q1, out }) => `${out} = clamp(${out}, min(${q0}, ${q1}), max(${q0}, ${q1}));`).join(' ')}
    }
`;
const ppmLimitInline = ({ left, right, mean, out }) => `    l${out} = ${left}; r${out} = ${right};
    if ((r${out} - ${mean}) * (${mean} - l${out}) <= 0.0) { l${out} = ${mean}; r${out} = ${mean}; } else {
      let d = r${out} - l${out}; let excess = d * (${mean} - 0.5 * (l${out} + r${out}));
      if (excess > d * d / 6.0) { l${out} = 3.0 * ${mean} - 2.0 * r${out}; } else if (excess < -d * d / 6.0) { r${out} = 3.0 * ${mean} - 2.0 * l${out}; }
    }
    s${out} = r${out} - l${out}; c${out} = 6.0 * (${mean} - 0.5 * (l${out} + r${out}));
`;
const limiterVars = (out) => `var l${out} = 0.0; var r${out} = 0.0; var s${out} = 0.0; var c${out} = 0.0;`;
const pendingBit = `if (cur < 32) { pendLo |= 1u << u32(cur); } else { pendHi |= 1u << u32(cur - 32); }`;
const pendingTest = `let pending = select((pendHi >> u32(max(k - 32, 0))) & 1u, (pendLo >> u32(min(k, 31))) & 1u, k < 32);`;
// The solved tracer's implicit fluxes through the interfaces, taken from one store and given to the next, so the column's amount is kept.
const fluxInline = ({ h, kappa, q, amount, dt }) => `  for (var k = 0; k < L - 1; k++) {
    let moved = ${dt} * ${kappa}[k] * (${q}[k] - ${q}[k + 1]) / max(0.5 * (${h}[k] + ${h}[k + 1]), DZFLOOR);
    ${amount}[k] -= moved; ${amount}[k + 1] += moved;
  }
`;

function hybridHead(o) {
  const { L, table } = o;
  const arr = (name, values) => `const ${name}: array<f32, ${values.length}> = array<f32, ${values.length}>(${Array.from(values, (v) => Number(v).toExponential(10)).join(', ')});`;
  const constLine = (name, value) => `const ${name}: f32 = ${Number(value).toExponential(10)};`;
  const boundaries = labelBoundaries(o.rho);
  const floors = Array.from(table.minimum, (m) => Math.min(PV_FLOOR, m));
  return `
${arr('MINIMUM', table.minimum)}
const FIXEDTOP: i32 = ${table.fixed ?? 0};
${arr('TOKENEDGE', Array.from(table.minimum, (m) => 2 * Math.min(m, table.spacing ?? m)))}
${arr('DEPTHS', table.depths)}
${arr('PVFLOORK', floors)}
${arr('BOUNDARY', boundaries.subarray(1, L))}
const GEN: i32 = ${GENERATORS.indexOf(o.generator)};
const LP1: i32 = ${L + 1};
${constLine('HMIN', o.minimumLayer)} ${constLine('THINEDGE', 2 * o.minimumLayer)} ${constLine('MODEDEPTH', o.modeDepth)} ${constLine('MSTAR', o.mstar)} ${constLine('EKMANF', o.ekmanFactor)} ${constLine('RH18C2', RH18_NEUTRAL.c2)} ${constLine('RH18C3', RH18_NEUTRAL.c3)} ${constLine('ICEWAVES', ICE_WAVE_DAMPING)} ${constLine('NSTAR', o.nstar)}
${constLine('BACKKAPPA', o.backgroundDiffusivity)} ${constLine('CONVKAPPA', o.convectiveDiffusivity)} ${constLine('VONK', 0.4)} ${constLine('BLFLOOR', 0.1)} ${constLine('DZFLOOR', 1e-3)} ${constLine('HNEGLECT', 1e-6)} ${constLine('MLEXCESS', 0.03)}
${constLine('BLTOL', 1e-6)} const BLITER: i32 = 24;
fn eosA(t: f32, s: f32) -> f32 { return eos(t + TREF, s + SREF); }
fn modeWeight(k: i32, i: i32) -> f32 {
  if (GEN == 0) { return 1.0; }
  if (GEN == 1) { return 0.0; }
  return clamp((IN[hOff(k - 1) + i] - MINIMUM[k - 1]) / MODEDEPTH, 0.0, 1.0);
}
fn crossing(qL: f32, d: f32, c: f32, b: f32) -> f32 {
  let aim = b - qL;
  if (abs(c) < 1.0e-6 * max(abs(d), 1.0e-30)) { return select(0.0, clamp(aim / d, 0.0, 1.0), d > 0.0); }
  let p = d + c; let discriminant = p * p - 4.0 * c * aim;
  var xi = select(0.5, (p - sqrt(max(discriminant, 0.0))) / (2.0 * c), discriminant >= 0.0);
  if (!(xi >= 0.0 && xi <= 1.0)) {
    var low = 0.0; var high = 1.0;
    for (var it = 0; it < 40; it++) { xi = 0.5 * (low + high); if (xi * (d + c * (1.0 - xi)) < aim) { low = xi; } else { high = xi; } }
  }
  return xi;
}
fn blEnergy(depth: f32, ustar: f32, f: f32, dt: f32, convective: f32, waves: ${o.langmuirIce === 'tanh' ? 'f32' : 'bool'}) -> f32 {
${o.rotationLimit === 'rh18' ? `  let decay = select(0.0, exp(RH18C3 * depth * abs(f) / max(ustar, 1.0e-12)), ustar > 0.0);
  let factor = ${o.langmuir > 0 ? waveFactor(o) : '1.0'};
  return RHO0 * dt * MSTAR * factor * ustar * ustar * ustar * (1.0 + RH18C2) * decay / (1.0 + RH18C2 * decay) + NSTAR * convective;` : `  let rotation = select(1.0e30, 1.0 + (depth * f / (EKMANF * ustar)) * (depth * f / (EKMANF * ustar)), ustar > 0.0);
  let factor = ${o.langmuir > 0 ? waveFactor(o) : '1.0'};
  return RHO0 * dt * MSTAR * factor * ustar * ustar * ustar / rotation + NSTAR * convective;`}
}
`;
}

function hybridKernels(o, K, idx) {
  const { L, C, E } = o, kept = o.edgeTransport === 'kept';
  const edgeFlow = (thickness, sum) => `  { var zAbove = 0.0;
    for (var k = 0; k < L; k++) {
      let ha = ${thickness('a')}; let hb = ${thickness('b')}; let hk = select(min(ha, hb), 0.5 * (ha + hb), k == 0);
      let he = min(hk, max(0.0, sillTop - zAbove)); zAbove += hk; ${sum}
    } }
`;
  const sillTop = kept ? `  let sillTop = max(EPSO, OD[O_SILL + e] + 0.5 * (OD[O_ETA + a] + OD[O_ETA + b]));
` : '';
  const stateThickness = (c) => `IN[hOff(k) + ${c}]`, targetThickness = (c) => `OD[O_TARGET + (k + 1) * C + ${c}] - OD[O_TARGET + k * C + ${c}]`;
  const keepFlow = (thickness, before) => `  var keptAfter = 0.0; var keptDepth = 0.0;
${edgeFlow(thickness, 'keptAfter += he * IN[uOff(k) + e]; keptDepth += he;')}  if (keptDepth > EPSO) { let shift = (${before} - keptAfter) / keptDepth; for (var k = 0; k < L; k++) { IN[uOff(k) + e] += shift; } }
`;
  const cellVector = (field) => `vec3<f32>(0.0, 0.0, 0.0);
  for (var m = 0; m < MAXE; m++) {
    let e = MI[EOC + MAXE * i + m];
    let s = abs(f32(MI[ESC + MAXE * i + m])) * 0.5 * MF[F_DC + e] * MF[F_DV + e] * ${field};
    wv += s * vec3<f32>(MF[F_NEDGE + 3 * e], MF[F_NEDGE + 3 * e + 1], MF[F_NEDGE + 3 * e + 2]);
  }
  wv = wv / MF[F_AREA + i];`;
  // boundaryLayer's work: the elimination down the layers the depth couples, setting CK at each interface it reaches and stopping at the first it does not, then the substitution back up summing the potential energy the trial profile raises; NL the layers solved.
  const blWork = `      var z = 0.0; var wPrev = 0.0; var gPrev = 0.0; var rT = 0.0; var rS = 0.0;
      NL = 0;
      for (var k = 0; k < L; k++) {
        let hk = CH[k]; z += hk;
        let inside = k < L - 1 && z < depth;
        let kap = select(0.0, VONK * vstar * z * (1.0 - z / depth) * (1.0 - z / depth), inside);
        CK[k] = kap;
        let w = select(0.0, dt * kap / max(0.5 * (hk + CH[min(k + 1, L - 1)]), DZFLOOR), inside);
        let m = hk + wPrev + w + wPrev * gPrev; let g = -w / m;
        rT = (hk * CT[k] + wPrev * rT) / m; CTT[k] = rT;
        rS = (hk * CS[k] + wPrev * rS) / m; CTS[k] = rS;
        SG[k] = g; wPrev = w; gPrev = g;
        NL = k + 1;
        if (!inside) { break; }
      }
      var raised = 0.0; var nT = CTT[NL - 1]; var nS = CTS[NL - 1];
      for (var j = 0; j < NL - 1; j++) {
        let k = NL - 2 - j; let tBelow = nT; let sBelow = nS;
        nT = CTT[k] - SG[k] * nT; nS = CTS[k] - SG[k] * nS;
        let lift = 0.5 * (CH[k] + CH[k + 1]); let dz = max(lift, DZFLOOR);
        let heatDown = dt * CK[k] * (nT - tBelow) / dz; let saltDown = dt * CK[k] * (nS - sBelow) / dz;
        raised += OGRAV * RHO0 * lift * (CA[k] * heatDown - CB[k] * saltDown);
      }
      let supplied = blEnergy(depth, ustar, f, dt, convective, ${o.langmuirIce === 'tanh' ? 'waves' : '!capped'}); let worked = max(raised, 0.0);
      let excess = (worked - supplied) / (worked + supplied);
`;
  return {
    oKineticPhi: `${K}  let i = ${idx}; if (i >= C) { return; }
  var edges: array<i32, MAXE>; var signs: array<f32, MAXE>; var dcs: array<f32, MAXE>; var dvs: array<f32, MAXE>;
  for (var m = 0; m < MAXE; m++) { let e = MI[EOC + MAXE * i + m]; edges[m] = e; signs[m] = f32(MI[ESC + MAXE * i + m]); dcs[m] = MF[F_DC + e]; dvs[m] = MF[F_DV + e]; }
  let area = MF[F_AREA + i];
  for (var k = 0; k < L; k++) {
    var kinetic = 0.0;
    for (var m = 0; m < MAXE; m++) { let u = IN[uOff(k) + edges[m]]; kinetic += abs(signs[m]) * 0.25 * dcs[m] * dvs[m] * u * u; }
    OD[O_PHI + k * C + i] = kinetic / area;
  }
}`,
    /*
     * hybridTendencyLayers' momentum: the PV flux, the kinetic energy and
     * free-surface gradients, the pressure force at the depth the layer
     * occupies in both cells or the middle of the gap between them, the
     * stress on layer 0, the bottom drag on the part of the layer within
     * MOMDEPTH of the edge column's floor, the ∇⁴ closure, and the pin to
     * the layer above of a layer no thicker than THINEDGE, or no thicker
     * than TOKENEDGE while its two cells' spans, each measured down from its
     * column's top, share no depth; P[2] is relax.
     */
    oMomentum: `${K}  let e = ${idx}; if (e >= E) { return; }
  if (OD[O_EMASK + e] < 0.5) { for (var k = 0; k < L; k++) { OUT[uOff(k) + e] = 0.0; } return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  let dc = MF[F_DC + e];
  var others: array<i32, MAXEE>; var weights: array<f32, MAXEE>; var live: array<bool, MAXEE>;
  for (var s = 0; s < MAXEE; s++) { let slot = MAXEE * e + s; others[s] = MI[EOE + slot]; weights[s] = MF[F_PVW + slot]; live[s] = OD[O_EMASK + others[s]] > 0.5; }
  var floor = 0.0;
  for (var k = 0; k < L; k++) { floor += OD[O_HEDGE + k * E + e]; }
  let da = OD[O_DRAFT + a]; let db = OD[O_DRAFT + b]; var layerTop = 0.0; var uAbove = 0.0; var Ha = da; var Hb = db;
  for (var k = 0; k < L; k++) {
    let n = k * E + e;
    let u = IN[uOff(k) + e];
    let qHere = 0.5 * OD[O_QE + n];
    var pv = 0.0;
    for (var s = 0; s < MAXEE; s++) { let other = others[s]; let fo = select(0.0, OD[O_HEDGE + k * E + other] * IN[uOff(k) + other], live[s]); pv += weights[s] * fo * (qHere + 0.5 * OD[O_QE + k * E + other]); }
    let gradK = (OD[O_PHI + k * C + b] - OD[O_PHI + k * C + a]) / dc;
    let pa = O_PR + 2 * (k * C + a); let pb = O_PR + 2 * (k * C + b); let ha = IN[hOff(k) + a]; let hb = IN[hOff(k) + b];
    let top = max(Ha, Hb); let bottom = min(Ha + ha, Hb + hb);
    let zR = select(0.5 * (Ha + 0.5 * ha + Hb + 0.5 * hb), 0.5 * (top + bottom), bottom > top);
    let force = OD[pb] + OD[pb + 1] * (zR - Hb) - OD[pa] - OD[pa + 1] * (zR - Ha);
    var du = pv / dc - gradK - OGRAV * OD[O_GRADETA + e] - OGRAV / RHO0 * force / dc;
    let he = OD[O_HEDGE + n];
    if (k == 0) { du += OD[O_STRESS + e] / RHO0 / max(he, HMIN) - TOPDRAG * 0.5 * (OD[O_SHELF + a] + OD[O_SHELF + b]) * abs(u) * u / max(he, HMIN); }
    let overlap = min(layerTop + he, floor) - max(layerTop, floor - MOMDEPTH);
    if (overlap > 0.0) { du -= RBOT * abs(u) * u * (overlap / MOMDEPTH) / max(he, HMIN); }
    if (NU4O > 0.0) { du -= NU4O * OD[O_LAPB + n]; }
    let ra = Ha - da; let rb = Hb - db;
    if (k > 0 && (he <= THINEDGE || (he <= TOKENEDGE[k] && !(min(ra + ha, rb + hb) > max(ra, rb))))) { du = (uAbove - u) * P[2]; }
    OUT[uOff(k) + e] = du;
    uAbove = u; layerTop += he; Ha += ha; Hb += hb;
  }
}`,
    /*
     * oEdgeThickness, oFlux and oEdgePV in one sweep of the edge column:
     * the edge thicknesses clipped to the sill, the mass flux from the
     * donor cell, and the edge PV from the vertex vorticities oVertexVort
     * has left in AVORT.
     */
    oEdgeFluxPV: `${K}  let e = ${idx}; if (e >= E) { return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1]; let v1 = MI[VOE + 2 * e]; let v2 = MI[VOE + 2 * e + 1];
  let sill = max(EPSO, OD[O_SILL + e] + 0.5 * (OD[O_ETA + a] + OD[O_ETA + b]));
  let wet = OD[O_EMASK + e] > 0.5;
  var above = 0.0;
  for (var k = 0; k < L; k++) {
    let n = k * E + e; let ha = IN[hOff(k) + a]; let hb = IN[hOff(k) + b];
    let hk = select(min(ha, hb), 0.5 * (ha + hb), k == 0);
    let he = min(hk, max(0.0, sill - above)); above += hk;
    OD[O_HEDGE + n] = he;
    let ue = IN[uOff(k) + e];
    OD[O_FLUX + n] = select(0.0, min(he, max(0.0, select(hb, ha, ue > 0.0))) * ue, wet);
    let hq = max(he, CENTRING * 0.5 * (ha + hb));
    OD[O_QE + n] = 0.5 * (OD[O_AVORT + k * V + v1] + OD[O_AVORT + k * V + v2]) / max(hq, PVFLOORK[k]);
  }
}`,
    oRescale: `${K}  let i = ${idx}; if (i >= C) { return; }
  if (OD[O_CMASK + i] < 0.5) { return; }
  var sum = 0.0; var dried = 0.0;
  for (var k = 0; k < L; k++) {
    var hv = IN[hOff(k) + i];
    if (!(hv > 0.0)) { if (hv < 0.0) { dried = 1.0; } hv = 0.0; IN[hOff(k) + i] = 0.0; IN[qOff(k) + i] = 0.0; IN[wOff(k) + i] = 0.0; }
    sum += hv;
  }
  let excess = (OD[B_BAVG + i] - (sum - OD[O_BATH + i])) / sum;
  for (var k = 0; k < L; k++) { IN[hOff(k) + i] += IN[hOff(k) + i] * excess; IN[qOff(k) + i] += IN[qOff(k) + i] * excess; IN[wOff(k) + i] += IN[wOff(k) + i] * excess; }
  OD[O_ETA + i] = OD[B_BAVG + i];
  OD[O_DRIED + i] = dried;
}`,
    /*
     * regridTargets: the column's target interfaces from its generator
     * (targetInterfaces in ocean/column.module.js) into TARGET; P[0] is
     * the regrid's relax. The sorted density's reconstruction is produced
     * layer by layer as the boundaries' crossing search descends.
     */
    oRegridTargets: (() => {
      const load = (n, s) => `  { let n = ${n}; var hv = 0.0;
    if (n < L) { hv = IN[hOff(n) + i]; if (hv > 1e-9) { tLoad = IN[qOff(n) + i] / hv; sLoad = IN[wOff(n) + i] / hv; } rMax = max(rMax, eosA(tLoad, sLoad)); }
    h${s} = hv; r${s} = rMax; }
`;
      const produce = `  {
${ppmEdgeInline({ k: 'j + 1', hm1: 'hm1', h0: 'h0', h1: 'h1', h2: 'h2', tracers: [{ qm1: 'rm1', q0: 'r0', q1: 'r1', q2: 'r2', out: 'eR1' }] })}${ppmLimitInline({ left: 'eR', right: 'eR1', mean: 'r0', out: 'R' })}    let below = select(0.0, select(clamp((h0 - MINIMUM[j]) / MODEDEPTH, 0.0, 1.0), 1.0, GEN == 0), j + 1 < L);
    let aboveW = select(0.0, select(clamp((hm1 - MINIMUM[max(j - 1, 0)]) / MODEDEPTH, 0.0, 1.0), 1.0, GEN == 0), j > 0);
    let keep = 1.0 - min(aboveW, below);
    lR = r0 + keep * (lR - r0); sR *= keep; cR *= keep;
  }
`;
      return `${K}  let i = ${idx}; if (i >= C) { return; }
  if (OD[O_CMASK + i] < 0.5) { return; }
  var H = 0.0;
  for (var k = 0; k < L; k++) { H += IN[hOff(k) + i]; }
  let depth = OD[O_BATH + i]; let spacing = min(HMIN, H / f32(L));
  OD[O_TARGET + i] = 0.0; OD[O_TARGET + L * C + i] = H;
  var prev = 0.0;
  if (GEN == 1) {
    for (var k = 1; k < L; k++) {
      let fwd = max(min(DEPTHS[k], depth) * H / depth, prev + spacing);
      OD[O_TARGET + k * C + i] = min(fwd, H - f32(L - k) * spacing); prev = fwd;
    }
    return;
  }
  var tLoad = LABEL_T[0] - TREF; var sLoad = LABEL_S[0] - SREF; var rMax = -1.0e30;
  var hm1 = 0.0; var rm1 = 0.0; var h0 = 0.0; var r0 = 0.0; var h1 = 0.0; var r1 = 0.0; var h2 = 0.0; var r2 = 0.0;
${load('0', '0')}${load('1', '1')}${load('2', '2')}  var j = 0; var zj = 0.0; var eR = r0; var eR1 = 0.0; ${limiterVars('R')}
${produce}  var above = 0.0; var zk = 0.0;
  for (var k = 1; k < L; k++) {
    zk += IN[hOff(k - 1) + i];
    let bk = BOUNDARY[k - 1];
    while (j < L && lR + sR < bk) {
      j++;
      if (j < L) {
        zj += h0; eR = eR1;
        hm1 = h0; rm1 = r0; h0 = h1; r0 = r1; h1 = h2; r1 = r2;
${load('j + 2', '2')}${produce}      }
    }
    var d = H;
    if (j < L) { if (bk <= lR) { d = zj; } else { d = zj + h0 * crossing(lR, sR, cR, bk); } }
    let floor = select(0.0, above + MINIMUM[k - 1], GEN == 2);
    var wanted = max(d, floor);
    if (GEN == 2 && k <= FIXEDTOP) { wanted = floor; } else if (wanted > floor) { wanted = max(floor, zk + P[0] * (wanted - zk)); }
    above = wanted;
    let fwd = max(wanted, prev + spacing);
    OD[O_TARGET + k * C + i] = min(fwd, H - f32(L - k) * spacing); prev = fwd;
  }
}`;
    })(),
    // regridEdges: each edge's velocity remapped piecewise-linearly onto the mean of its cells' new thicknesses.
    oRegridEdges: (() => {
      const load = (n, s) => `  { let n = ${n}; var hv = 0.0; var uv = 0.0; if (n < L) { hv = 0.5 * (IN[hOff(n) + a] + IN[hOff(n) + b]); uv = IN[uOff(n) + e]; } h${s} = hv; u${s} = uv; }
`;
      const flush = (safe) => `      { let value = select(IN[uOff(cur) + e], acc / hn, hn > 0.0);
        if (${safe}) { IN[uOff(cur) + e] = value; } else { OUT[uOff(cur) + e] = value; ${pendingBit} }
        cur++; znLo = znHi; taL = taH; tbL = tbH; let kn = min(cur + 1, L); taH = OD[O_TARGET + kn * C + a]; tbH = OD[O_TARGET + kn * C + b];
        hn = 0.5 * (taH - taL + tbH - tbL); znHi = znLo + hn; acc = 0.0; }
`;
      return `${K}  let e = ${idx}; if (e >= E) { return; }
  if (OD[O_EMASK + e] < 0.5) { return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
${kept ? `${sillTop}  var keptBefore = 0.0;
${edgeFlow(stateThickness, 'keptBefore += he * IN[uOff(k) + e];')}` : ''}  var hm1 = 0.0; var um1 = 0.0; var h0 = 0.0; var u0 = 0.0; var h1 = 0.0; var u1 = 0.0;
${load('0', '0')}${load('1', '1')}  var z = 0.0; var cur = 0; var znLo = 0.0;
  var taL = 0.0; var tbL = 0.0; var taH = OD[O_TARGET + C + a]; var tbH = OD[O_TARGET + C + b];
  var hn = 0.5 * (taH - taL + tbH - tbL); var znHi = hn;
  var acc = 0.0; var pendLo = 0u; var pendHi = 0u;
  for (var j = 0; j < L; j++) {
    var s = 0.0;
    if (j > 0 && j < L - 1) {
      let below = u1 - u0; let above = u0 - um1;
      if (below * above > 0.0) {
        let central = h0 * (u1 - um1) / max(0.5 * (hm1 + 2.0 * h0 + h1), DZFLOOR);
        s = sign(central) * min(abs(central), min(2.0 * abs(below), 2.0 * abs(above)));
      }
    }
    let lU = u0 - 0.5 * s;
    let zLo = z; let zHi = z + h0;
    if (h0 > 0.0) {
      while (cur < L - 1 && znHi <= zLo) {
${flush('cur <= j + 1')}      }
      var given = 0.0;
      loop {
        let lo = max(zLo, znLo); let hi = min(zHi, znHi);
        if (hi > lo) {
          let pa = clamp((lo - zLo) / h0, 0.0, 1.0); let pb = clamp((hi - zLo) / h0, 0.0, 1.0);
          let p = h0 * (lU * (pb - pa) + s * 0.5 * (pb * pb - pa * pa));
          acc += p; given += p;
        }
        if (znHi >= zHi || cur >= L - 1) { break; }
${flush('cur <= j + 1')}      }
      acc += h0 * u0 - given;
    }
    z = zHi; hm1 = h0; um1 = u0; h0 = h1; u0 = u1;
${load('j + 2', '1')}  }
  while (cur < L) {
${flush('true')}  }
  if ((pendLo | pendHi) != 0u) {
    for (var k = 0; k < L; k++) { ${pendingTest} if (pending != 0u) { IN[uOff(k) + e] = OUT[uOff(k) + e]; } }
  }
${kept ? keepFlow(targetThickness, 'keptBefore') : ''}}`;
    })(),
    // regridCells: heat and salt remapped by the parabolas onto the targets, which set the thicknesses; the interface geometry is shared between the two tracers.
    oRegridCells: (() => {
      const load = (n, s) => `  { let n = ${n}; var hv = 0.0; var qv = 0.0; var wv = 0.0;
    if (n < L) { hv = IN[hOff(n) + i]; qv = IN[qOff(n) + i]; wv = IN[wOff(n) + i]; if (hv > 1e-9) { tLoad = qv / hv; sLoad = wv / hv; } }
    h${s} = hv; q${s} = qv; w${s} = wv; t${s} = tLoad; s${s} = sLoad; }
`;
      const flush = (safe) => `      { let hn = znHi - znLo;
        if (${safe}) { IN[hOff(cur) + i] = hn; IN[qOff(cur) + i] = accT; IN[wOff(cur) + i] = accS; }
        else { OUT[hOff(cur) + i] = hn; OUT[qOff(cur) + i] = accT; OUT[wOff(cur) + i] = accS; ${pendingBit} }
        cur++; znLo = znHi; znHi = OD[O_TARGET + min(cur + 1, L) * C + i]; accT = 0.0; accS = 0.0; }
`;
      return `${K}  let i = ${idx}; if (i >= C) { return; }
  if (OD[O_CMASK + i] < 0.5) { return; }
  var tLoad = LABEL_T[0] - TREF; var sLoad = LABEL_S[0] - SREF;
  var hm1 = 0.0; var tm1 = 0.0; var sm1 = 0.0;
  var h0 = 0.0; var t0 = 0.0; var s0 = 0.0; var q0 = 0.0; var w0 = 0.0;
  var h1 = 0.0; var t1 = 0.0; var s1 = 0.0; var q1 = 0.0; var w1 = 0.0;
  var h2 = 0.0; var t2 = 0.0; var s2 = 0.0; var q2 = 0.0; var w2 = 0.0;
${load('0', '0')}${load('1', '1')}${load('2', '2')}  var eT = t0; var eS = s0; var eT1 = 0.0; var eS1 = 0.0; var z = 0.0;
  ${limiterVars('T')} ${limiterVars('S')}
  var cur = 0; var znLo = 0.0; var znHi = OD[O_TARGET + C + i];
  var accT = 0.0; var accS = 0.0; var pendLo = 0u; var pendHi = 0u;
  for (var j = 0; j < L; j++) {
${ppmEdgeInline({ k: 'j + 1', hm1: 'hm1', h0: 'h0', h1: 'h1', h2: 'h2', tracers: [{ qm1: 'tm1', q0: 't0', q1: 't1', q2: 't2', out: 'eT1' }, { qm1: 'sm1', q0: 's0', q1: 's1', q2: 's2', out: 'eS1' }] })}${ppmLimitInline({ left: 'eT', right: 'eT1', mean: 't0', out: 'T' })}${ppmLimitInline({ left: 'eS', right: 'eS1', mean: 's0', out: 'S' })}    let zLo = z; let zHi = z + h0;
    if (h0 > 0.0) {
      while (cur < L - 1 && znHi <= zLo) {
${flush('cur <= j + 2')}      }
      var givenT = 0.0; var givenS = 0.0;
      loop {
        let lo = max(zLo, znLo); let hi = min(zHi, znHi);
        if (hi > lo) {
          let pa = clamp((lo - zLo) / h0, 0.0, 1.0); let pb = clamp((hi - zLo) / h0, 0.0, 1.0);
          let d1 = pb - pa; let d2 = 0.5 * (pb * pb - pa * pa); let d3 = (pb * pb * pb - pa * pa * pa) / 3.0;
          let pT = h0 * (lT * d1 + sT * d2 + cT * (d2 - d3)); let pS = h0 * (lS * d1 + sS * d2 + cS * (d2 - d3));
          accT += pT; givenT += pT; accS += pS; givenS += pS;
        }
        if (znHi >= zHi || cur >= L - 1) { break; }
${flush('cur <= j + 2')}      }
      accT += q0 - givenT; accS += w0 - givenS;
    }
    z = zHi; eT = eT1; eS = eS1;
    hm1 = h0; tm1 = t0; sm1 = s0;
    h0 = h1; t0 = t1; s0 = s1; q0 = q1; w0 = w1;
    h1 = h2; t1 = t2; s1 = s2; q1 = q2; w1 = w2;
${load('j + 3', '2')}  }
  while (cur < L) {
${flush('true')}  }
  if ((pendLo | pendHi) != 0u) {
    for (var k = 0; k < L; k++) { ${pendingTest} if (pending != 0u) { IN[hOff(k) + i] = OUT[hOff(k) + i]; IN[qOff(k) + i] = OUT[qOff(k) + i]; IN[wOff(k) + i] = OUT[wOff(k) + i]; } }
  }
}`;
    })(),
    /*
     * columnCells: mixColumn of ocean/column.module.js on each sea column
     * with the cell's stress and the shear of its reconstructed velocities:
     * convective adjustment, the Pacanowski–Philander interior and the
     * interfaces' expansion coefficients computed as the velocities are
     * gathered, the energetic boundary layer, one implicit solve of heat
     * and salt together whose substitution moves the stores' amounts and
     * writes them; the viscosity profile into VISC for the edges, the
     * mixed-layer and boundary-layer depths into MIXED and BLD; P[0] is
     * the step. The boundary layer's loop is the search of boundaryLayer
     * in ocean/column.module.js, phase for phase: the excess
     * (work − energy)/(work + energy) at the whole column (0), the floor
     * (1), the previous step's depth from BLD (2) and half (3) or twice (4)
     * it, the whole column again (5), the Anderson–Björck secant points
     * until the bracket is BLTOL wide in log depth or BLITER have been
     * tried (6), and the depth settled on (7), so that CK holds its
     * diffusivity over the NL layers that trial solved. The slabs of the
     * convective adjustment borrow CK (their first layer), CTT (thickness),
     * CTS (heat), CA (salt) and CB (density), which are filled afterward.
     */
    oColumnCells: `${K}  let i = ${idx}; if (i >= C) { return; }
  if (OD[O_CMASK + i] < 0.5) { OD[O_MIXED + i] = 0.0; OD[O_BLD + i] = 0.0; return; }
  ${declareArrays(['CH', 'CT', 'CS', 'CQ', 'CW', 'CA', 'CB', 'CK', 'CIK', 'CINU', 'CTT', 'CTS', 'SG'])}
  var edges: array<i32, MAXE>; var weights: array<f32, MAXE>; var normals: array<vec3<f32>, MAXE>;
  let dt = P[0];
  var tLoad = LABEL_T[0] - TREF; var sLoad = LABEL_S[0] - SREF; var H = 0.0;
  for (var k = 0; k < L; k++) {
    let hv = IN[hOff(k) + i]; let qv = IN[qOff(k) + i]; let wv = IN[wOff(k) + i];
    CH[k] = hv; CQ[k] = qv; CW[k] = wv; H += hv;
    if (hv > 1e-9) { tLoad = qv / hv; sLoad = wv / hv; }
    CT[k] = tLoad; CS[k] = sLoad;
  }
  let area = MF[F_AREA + i];
  for (var m = 0; m < MAXE; m++) {
    let e = MI[EOC + MAXE * i + m];
    edges[m] = e; weights[m] = abs(f32(MI[ESC + MAXE * i + m])) * 0.5 * MF[F_DC + e] * MF[F_DV + e];
    normals[m] = vec3<f32>(MF[F_NEDGE + 3 * e], MF[F_NEDGE + 3 * e + 1], MF[F_NEDGE + 3 * e + 2]);
  }
  var sv = vec3<f32>(0.0, 0.0, 0.0); var vt = vec3<f32>(0.0, 0.0, 0.0);
  for (var m = 0; m < MAXE; m++) { sv += weights[m] * OD[O_STRESS + edges[m]] * normals[m]; vt += weights[m] * IN[uOff(0) + edges[m]] * normals[m]; }
  sv = sv / area; vt = vt / area;
  let capped = OD[O_SHELF + i] > 0.5; let ustar = select(sqrt(length(sv) / RHO0), sqrt(TOPDRAG) * length(vt), capped); let f = abs(OD[O_FCELL + i]);${o.langmuirIce === 'tanh' ? ' let waves = select(1.0 - tanh(ICEWAVES * select(0.0, select(1.0, PH[PH_CONC + i], PH[PH_CONC + i] > 0.0), OD[O_ICED + i] > 0.5)), 0.0, capped);' : ''}
${o.shelfMelt ? `  var freezing = false;
  if (capped) {
${o.shelfMeltLayer > 0 ? `    var slab = 0.0; var sumT = 0.0; var sumS = 0.0;
    for (var k = 0; k < L; k++) { if (slab >= MELTLAYER) { break; } let part = min(CH[k], MELTLAYER - slab); sumT += part * CT[k]; sumS += part * CS[k]; slab += part; }
    let meanT = sumT / slab; let meanS = sumS / slab + SREF;
    let melt = basalMelt(meanT, meanS, OD[O_LIQUIDUS + i], sqrt(TOPDRAG * (dot(vt, vt) + TIDE2)));
    let m = melt.x * dt; let water = m / FRESHRHO;
    let meltHeat = -(OLFUS * m / RHOCP + water * (meanT - melt.y)); let meltSalt = -meanS * water;
    var given = 0.0;
    for (var k = 0; k < L; k++) {
      if (given >= slab) { break; }
      let part = min(CH[k], slab - given);
      if (part > 0.0) { CQ[k] += meltHeat * part / slab; CW[k] += meltSalt * part / slab; CT[k] = CQ[k] / CH[k]; CS[k] = CW[k] / CH[k]; }
      given += part;
    }
` : `    let melt = basalMelt(CT[0], CS[0] + SREF, OD[O_LIQUIDUS + i], sqrt(TOPDRAG * (dot(vt, vt) + TIDE2)));
    let m = melt.x * dt; let water = m / FRESHRHO;
    CQ[0] -= OLFUS * m / RHOCP + water * (CT[0] - melt.y); CW[0] -= (CS[0] + SREF) * water;
    CT[0] = CQ[0] / CH[0]; CS[0] = CW[0] / CH[0];
`}    OD[O_MELT + i] = m; OD[O_MELTSUM + i] += m; freezing = m < 0.0;
  }
` : ''}  var n = 0;
  for (var k = 0; k < L; k++) {
    CK[n] = f32(k); CTT[n] = CH[k]; CTS[n] = CQ[k]; CA[n] = CW[k]; CB[n] = eosA(CT[k], CS[k]);
    n++;
    loop {
      if (!(n >= 2 && CB[n - 2] > CB[n - 1])) { break; }
      CTT[n - 2] += CTT[n - 1]; CTS[n - 2] += CTS[n - 1]; CA[n - 2] += CA[n - 1];
      CB[n - 2] = select(CB[n - 2], eosA(CTS[n - 2] / CTT[n - 2], CA[n - 2] / CTT[n - 2]), CTT[n - 2] > 0.0);
      n--;
    }
  }
  var released = 0.0; var zTop = 0.0;
  for (var s = 0; s < n; s++) {
    let first = i32(CK[s]); let past = select(L, i32(CK[min(s + 1, L - 1)]), s + 1 < n);
    if (past - first < 2) { zTop += CH[first]; continue; }
    let t = CTS[s] / CTT[s]; let w = CA[s] / CTT[s]; let mixed = eosA(t, w);
    var drop = 0.0; var sumQ = 0.0; var sumW = 0.0;
    for (var k = first; k < past; k++) {
      let hk = CH[k]; drop += OGRAV * (zTop + 0.5 * hk) * hk * (mixed - eosA(CT[k], CS[k])); zTop += hk;
      CT[k] = t; CS[k] = w; CQ[k] = hk * t; CW[k] = hk * w; sumQ += CQ[k]; sumW += CW[k];
    }
    CQ[past - 1] += CTS[s] - sumQ; CW[past - 1] += CA[s] - sumW;
    if (first == 0) { released = drop; }
  }
  var previous = vec3<f32>(0.0, 0.0, 0.0); var rPrev = 0.0; var hPrev = 0.0; var aPrev = 0.0; var bPrev = 0.0; var mixing = false;
  for (var k = 0; k < L; k++) {
    let hk = CH[k]; let tk = CT[k] + TREF; let sk = CS[k] + SREF;
    let rk = eos(tk, sk); let ak = alphaT(tk, sk); let bk = betaS(tk, sk);
    var wv = vec3<f32>(0.0, 0.0, 0.0);
    for (var m = 0; m < MAXE; m++) { wv += weights[m] * IN[uOff(k) + edges[m]] * normals[m]; }
    wv = wv / area;
    if (k > 0) {
      let j = k - 1; let d = previous - wv;
      let dz = max(0.5 * (hPrev + hk), DZFLOOR); let buoyancy = OGRAV * (rk - rPrev) / RHO0;
      var kap = CONVKAPPA; var nu = CONVKAPPA;
      if (buoyancy >= 0.0) {
        let richardson = buoyancy * dz / (dot(d, d) + 1e-12);
        nu = SHEARNU / ((1.0 + 5.0 * richardson) * (1.0 + 5.0 * richardson)) + BACKNU;
        kap = nu / (1.0 + 5.0 * richardson) + BACKKAPPA;
      }
      CIK[j] = kap; CINU[j] = nu; CA[j] = 0.5 * (aPrev + ak); CB[j] = 0.5 * (bPrev + bk);
      if (kap > 0.0) { mixing = true; }
    }
    previous = wv; rPrev = rk; hPrev = hk; aPrev = ak; bPrev = bk;
  }
  CIK[L - 1] = 0.0; CINU[L - 1] = 0.0;
  let convective = select(max(0.0, released), 0.0, capped${o.shelfMelt ? ' && !freezing' : ''}); let vstar = pow(ustar * ustar * ustar + convective / (RHO0 * dt), 1.0 / 3.0);
  var depth = 0.0; var NL = 0;
  if (vstar > 0.0) {
    let lowest = min(BLFLOOR, H); let guess = OD[O_BLD + i]; let half = max(lowest, 0.5 * guess); let twice = min(H, 2.0 * guess);
    var phase = select(0, 2, guess > lowest && guess < H);
    var a = log(lowest); var b = log(H); var fa = 0.0; var fb = 0.0; var x = 0.0; var side = 0; var count = 0;
    depth = H;
    for (var it = 0; it < BLITER + 5; it++) {
      if (phase == 1) { depth = lowest; } else if (phase == 2) { depth = guess; } else if (phase == 3) { depth = half; } else if (phase == 4) { depth = twice; } else if (phase == 0 || phase == 5) { depth = H; }
      else if (phase == 6) {
        x = (a * fb - b * fa) / (fb - fa);
        if (!(x > a && x < b)) { x = 0.5 * (a + b); }
        if (!(x > a && x < b)) { phase = 7; } else { depth = exp(x); }
      }
      if (phase == 7) { depth = exp(a); }
${blWork}      if (phase == 0 || phase == 5) { if (excess <= 0.0) { break; } fb = excess; phase = select(6, 1, phase == 0); }
      else if (phase == 1) { if (excess >= 0.0) { depth = lowest; break; } fa = excess; phase = 6; }
      else if (phase == 2) { if (excess > 0.0) { b = log(guess); fb = excess; phase = 3; } else { a = log(guess); fa = excess; phase = 4; } }
      else if (phase == 3) { if (excess < 0.0) { a = log(half); fa = excess; phase = 6; } else { b = log(half); fb = excess; if (half <= lowest) { depth = lowest; break; } phase = 1; } }
      else if (phase == 4) { if (excess > 0.0) { b = log(twice); fb = excess; phase = 6; } else { a = log(twice); fa = excess; if (twice >= H) { depth = H; break; } phase = 5; } }
      else if (phase == 6) {
        if (excess > 0.0) { if (side > 0) { let m = 1.0 - excess / fb; fa *= select(0.5, m, m > 0.0); } b = x; fb = excess; side = 1; }
        else { if (side < 0) { let m = 1.0 - excess / fa; fb *= select(0.5, m, m > 0.0); } a = x; fa = excess; side = -1; }
        count++;
        if (b - a <= BLTOL || count >= BLITER) { phase = 7; }
      } else { break; }
    }
    if (depth <= BLFLOOR) { depth = 0.0; NL = 0; }
    if (NL > 1) { mixing = true; }
  }
  {
    var wPrev = 0.0; var gPrev = 0.0; var rT = 0.0; var rS = 0.0;
    for (var k = 0; k < L; k++) {
      let hk = CH[k]; let boundary = select(0.0, CK[k], k < NL); let kap = CIK[k] + boundary;
      CIK[k] = kap; OD[O_VISC + k * C + i] = CINU[k] + boundary;
      let w = select(0.0, dt * kap / max(0.5 * (hk + CH[min(k + 1, L - 1)]), DZFLOOR), k < L - 1);
      let m = hk + wPrev + w + wPrev * gPrev; let g = -w / m;
      rT = (hk * CT[k] + wPrev * rT) / m; rS = (hk * CS[k] + wPrev * rS) / m;
      CTT[k] = rT; CTS[k] = rS; SG[k] = g; wPrev = w; gPrev = g;
    }
    if (mixing) {
      var nT = CTT[L - 1]; var nS = CTS[L - 1]; var carriedQ = 0.0; var carriedW = 0.0;
      CT[L - 1] = nT; CS[L - 1] = nS;
      for (var j = 0; j < L - 1; j++) {
        let k = L - 2 - j; let tBelow = nT; let sBelow = nS;
        nT = CTT[k] - SG[k] * nT; nS = CTS[k] - SG[k] * nS; CT[k] = nT; CS[k] = nS;
        let dz = max(0.5 * (CH[k] + CH[k + 1]), DZFLOOR);
        let movedQ = dt * CIK[k] * (nT - tBelow) / dz; let movedW = dt * CIK[k] * (nS - sBelow) / dz;
        IN[qOff(k + 1) + i] = (CQ[k + 1] + movedQ) - carriedQ; IN[wOff(k + 1) + i] = (CW[k + 1] + movedW) - carriedW;
        carriedQ = movedQ; carriedW = movedW;
      }
      IN[qOff(0) + i] = CQ[0] - carriedQ; IN[wOff(0) + i] = CW[0] - carriedW;
    } else {
      for (var k = 0; k < L; k++) { IN[qOff(k) + i] = CQ[k]; IN[wOff(k) + i] = CW[k]; }
    }
  }
  var mixed = H; var hAbove = CH[0]; var zc = 0.5 * hAbove; var prior = eosA(CT[0], CS[0]); let limit = prior + MLEXCESS; var found = false;
  for (var k = 1; k < L; k++) {
    if (found) { continue; }
    let hk = CH[k]; let center = zc + 0.5 * (hAbove + hk); let rk = eosA(CT[k], CS[k]);
    if (rk > limit) { mixed = zc + (limit - prior) / (rk - prior) * (center - zc); found = true; }
    zc = center; prior = rk; hAbove = hk;
  }
  OD[O_MIXED + i] = mixed; OD[O_BLD + i] = depth;
}`,
    // columnEdges: the implicit solve of each edge's velocity with the mean of its cells' viscosity profiles, the thicknesses and viscosities read as the elimination descends.
    oColumnEdges: `${K}  let e = ${idx}; if (e >= E) { return; }
  if (OD[O_EMASK + e] < 0.5) { return; }
  ${declareArrays(['SQ', 'SG'])}
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1]; let dt = P[0];
${kept ? `${sillTop}  var keptBefore = 0.0;
${edgeFlow(stateThickness, 'keptBefore += he * IN[uOff(k) + e];')}` : ''}  var mixing = false; var wPrev = 0.0; var gPrev = 0.0; var r = 0.0;
  var hk = 0.5 * (IN[hOff(0) + a] + IN[hOff(0) + b]);
  for (var k = 0; k < L; k++) {
    let kn = min(k + 1, L - 1); let below = k < L - 1;
    let hNext = 0.5 * (IN[hOff(kn) + a] + IN[hOff(kn) + b]);
    let kap = 0.5 * (OD[O_VISC + k * C + a] + OD[O_VISC + k * C + b]);
    if (below && kap > 0.0) { mixing = true; }
    let w = select(0.0, dt * kap / max(0.5 * (hk + hNext), DZFLOOR), below);
    let m = hk + wPrev + w + wPrev * gPrev; let g = -w / m;
    r = (hk * IN[uOff(k) + e] + wPrev * r) / m;
    SQ[k] = r; SG[k] = g; wPrev = w; gPrev = g; hk = hNext;
  }
  if (mixing) {
    var next = SQ[L - 1]; IN[uOff(L - 1) + e] = next;
    for (var j = 0; j < L - 1; j++) { let k = L - 2 - j; next = SQ[k] - SG[k] * next; IN[uOff(k) + e] = next; }
${kept ? keepFlow(stateThickness, 'keptBefore') : ''}  }
}`,
  };
}

function oceanKernels(o) {
  const { L, C, E, V, OS, OD, B } = o, third = o.tracerAdvection === 'third', shearOn = o.shearDiffusivity > 0, lateralOn = o.lateralDiffusionSpeed > 0, shear = shearOn || lateralOn;
  if (L > 64) throw new Error(`oMomentum marks an edge's thick classes in 64 bits, and the ocean has ${L} classes`);
  const constLine = (name, value) => `const ${name}: f32 = ${Number(value).toExponential(10)};`;
  const jerlov = JERLOV[o.penetration ?? 'IB'];
  const rhoLine = `const RHO: array<f32, ${L}> = array<f32, ${L}>(${o.rho.map((v) => v.toFixed(6)).join(', ')});`;
  const labelLine = `const LABEL_T: array<f32, ${L}> = array<f32, ${L}>(${o.labelT.map((v) => v.toFixed(6)).join(', ')});\nconst LABEL_S: array<f32, ${L}> = array<f32, ${L}>(${o.labelS.map((v) => v.toFixed(6)).join(', ')});\nconst RHOA: array<f32, ${L}> = array<f32, ${L}>(${o.rho.map((v) => (v - SEAWATER.rho0).toFixed(6)).join(', ')});`;
  const langmuirLines = `fn erfcApprox(x: f32) -> f32 {
  let t = 1.0 / (1.0 + 0.3275911 * x);
  return t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429)))) * exp(-x * x);
}
fn langmuirFactor(ustar: f32, depth: f32) -> f32 {
  let u10 = ustar / ${FRICTION_PER_WIND.toExponential(10)};
  if (!(u10 > ${CALM_WIND.toExponential(10)})) { return 1.0; }
  let drift = 0.0162 * u10; let height = 0.0246 * u10 * u10; let mean = ${(1.296 * 0.877 * 9.81 / (2 * Math.PI * 1.075)).toExponential(10)} / u10;
  let k = 0.176 * drift / (${(0.125 * Math.PI * 0.667).toExponential(10)} * mean * height * height); let ks = 2.56 * k; let z = 0.2 * depth;
  let r1 = (0.151 / (k * z) - 0.84) * (1.0 - exp(-2.0 * k * z));
  let r2 = -(0.84 + 0.0591 / (k * z)) * sqrt(${(2 * Math.PI).toExponential(10)} * k * z) * erfcApprox(sqrt(2.0 * k * z));
  let r3 = (0.0632 / (ks * z) + 0.125) * (1.0 - exp(-2.0 * ks * z));
  let r4 = (0.125 + 0.0946 / (ks * z)) * sqrt(${(2 * Math.PI).toExponential(10)} * ks * z) * erfcApprox(sqrt(2.0 * ks * z));
  let inverse = min(100.0, max(0.0, drift * (0.715 + r1 + r2 + r3 + r4) / ustar));
  return min(LANGMUIR, sqrt(1.0 + inverse / 2.25 + inverse * inverse / 850.3056));
}`;
  const langmuirHead = o.langmuir > 0 ? `${constLine('LANGMUIR', o.langmuir)}\n${langmuirLines}\n` : '';
  const offsetLines = Object.entries(OD).filter(([k]) => k !== 'total').map(([k, v]) => `const O_${k}: i32 = ${v};`).join('\n');
  const bLines = Object.entries(B).filter(([k]) => k !== 'total').map(([k, v]) => `const B_${k}: i32 = ${v};`).join('\n');
  const head = `
const L: i32 = ${L}; const DRY: f32 = -1.0e30;
const OH: i32 = ${OS.OH}; const OU: i32 = ${OS.OU}; const OQ: i32 = ${OS.OQ}; const OW: i32 = ${OS.OW}; const OSTOTAL: i32 = ${OS.total};
${offsetLines}
${bLines}
${rhoLine}
${labelLine}
${constLine('RHO0', o.density)} ${constLine('RHOCP', o.density * o.specificHeat)} ${constLine('OLFUS', o.latentHeatFusion)}
${constLine('OGRAV', o.gravity)}
${constLine('EPSO', EPS)} ${constLine('THINO', THIN)} ${constLine('PVFLOOR', PV_FLOOR)} ${constLine('SPEEDLIM', SPEED_LIMIT)}
${langmuirHead}${constLine('MOMDEPTH', o.momentumDepth)} ${constLine('CENTRING', o.vorticityCentring)} ${constLine('ICESAL', o.iceSalinity)} ${constLine('TRANSMIT', o.iceStressTransmission)} ${constLine('ICEDENS', o.iceDensity)}
${constLine('SHEARNU', o.shearViscosity)} ${constLine('BACKNU', o.backgroundViscosity)} ${constLine('RBOT', o.bottomDrag)} ${constLine('TOPDRAG', o.topDrag)} ${constLine('NU4O', o.nu4)}
${constLine('TIDE2', o.tidalSpeed * o.tidalSpeed)} ${constLine('FRESHRHO', SHELF_MELT.freshDensity)} ${constLine('MELTLAYER', o.shelfMeltLayer)}
${basalMeltWgsl({ rho0: o.density, rhoCp: o.density * o.specificHeat, latentHeat: o.latentHeatFusion, heatTransfer: o.shelfHeatTransfer, saltTransfer: o.shelfSaltTransfer })}
${constLine('FREEZE', FREEZING_POINT)} ${constLine('TREF', HEAT_REFERENCE)} ${constLine('SREF', SALT_REFERENCE)} ${constLine('GUSTO', o.gustiness)}
${constLine('JER_R1', jerlov[0])} ${constLine('JER_R2', 1 - jerlov[0])} ${constLine('JER_Z1', jerlov[1])} ${constLine('JER_Z2', jerlov[2])} ${constLine('PENDEPTH', PENETRATION_DEPTH)} ${constLine('JER_FLOOR', jerlov[0] * Math.exp(-PENETRATION_DEPTH / jerlov[1]) + (1 - jerlov[0]) * Math.exp(-PENETRATION_DEPTH / jerlov[2]))}
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
${shear ? `${constLine('SHEARDIFF', shearOn ? o.shearDiffusivity : 0)} ${constLine('SHEARCAP', o.shearDiffusivityCap)} ${constLine('LATSPEED', lateralOn ? o.lateralDiffusionSpeed : 0)} const FIXEDL: i32 = ${shearOn ? o.table.fixed : 0};
` : ''}${third ? `fn limitedValue(u: f32, d: f32, uu: f32, ratio: f32) -> f32 {
  let central = d - u;
  if (central == 0.0) { return u; }
  let r = (u - uu) * ratio / central;
  return u + 0.5 * max(0.0, min(min(2.0 * r, (2.0 + r) / 3.0), 2.0)) * central;
}
` : ''}fn iceVolume(i: i32) -> f32 {
  let h = OD[O_SURFICE + i]; let a = OD[O_SURFCONC + i];
  return select(0.0, h * select(1.0, a, a > 0.0), h > 0.0);
}
${SEAWATER_WGSL}
fn grain(x: f32) -> f32 { return ldexp(1.0, max(frexp(x).exp - 24, -126)); }
fn settle(x: f32, change: f32) -> f32 {
  let g = grain(max(abs(x), abs(x + change))); let r = round(change / g) * g;
  return select(r, r - (x / g - round(x / g)) * g, abs(r) <= abs(x));
}
fn transmitted(z: f32) -> f32 { let d = min(z, PENDEPTH); return JER_R1 * exp(-d / JER_Z1) + JER_R2 * exp(-d / JER_Z2) - JER_FLOOR; }
${hybridHead(o)}${o.drift ? iceDriftHead(o.drift) : ''}`;
  const idx = `(i32(id.x) + i32(id.y) * 4194240)`;
  const K = `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {\n`;
  const eddyShared = `${constLine('EDDYSLACK', EDDY_SLACK)}
fn eddyAllowance(k: i32, i: i32) -> f32 {
  return (max(0.0, IN[hOff(k) + i] - EPSO) + EDDYSLACK) * MF[F_AREA + i] / (f32(MI[NEC + i]) * P[0]);
}
fn eddyShare(xa: f32, xb: f32, sa: f32, sb: f32, moved: f32) -> f32 {
  let ga = grain(max(abs(xa), abs(xa - moved / sa))); let gb = grain(max(abs(xb), abs(xb + moved / sb)));
  let q = max(ga, gb); let area = select(sb, sa, ga > gb || (ga == gb && moved > 0.0));
  return round(moved / area / q) * q * area;
}
`;
  const eddyApply = (fromK) => `${K}  let n = ${idx}; if (n >= L * C) { return; }
  let k = n / C; let i = n % C;
  if (${fromK ? 'k == 0 || ' : ''}OD[O_CMASK + i] < 0.5) { return; }
  var volume = 0.0; var heat = 0.0; var salt = 0.0;
  for (var m = 0; m < MI[NEC + i]; m++) {
    let e = MI[EOC + MAXE * i + m]; let s = f32(MI[ESC + MAXE * i + m]);
    volume += s * OD[O_EDDYV + k * E + e]; heat += s * OD[O_EDDYQ + k * E + e]; salt += s * OD[O_EDDYW + k * E + e];
  }
  let factor = P[0] / MF[F_AREA + i];
  let dq = -factor * heat; let dw = -factor * salt;
  let rh = settle(IN[hOff(k) + i], -factor * volume); let rq = settle(IN[qOff(k) + i], dq); let rw = settle(IN[wOff(k) + i], dw);
  IN[hOff(k) + i] += rh; IN[qOff(k) + i] += rq; IN[wOff(k) + i] += rw;
  OUT[hOff(k) + i] = rh; OUT[qOff(k) + i] = dq - rq; OUT[wOff(k) + i] = dw - rw;
}`;
  const eddyLeftover = (fromK) => `${K}  let i = ${idx}; if (i >= C) { return; }
  if (OD[O_CMASK + i] < 0.5) { return; }
  var kept = 0.0; var owedQ = 0.0; var owedW = 0.0; var holder = ${fromK - 1}; var finest = 0.0;
  for (var k = ${fromK}; k < L; k++) {
    kept += OUT[hOff(k) + i]; owedQ += OUT[qOff(k) + i]; owedW += OUT[wOff(k) + i];
    let hk = IN[hOff(k) + i];
    if (hk > THINO && (holder == ${fromK - 1} || grain(hk) < finest)) { holder = k; finest = grain(hk); }
  }
  if (holder == ${fromK - 1}) { return; }
  let hh = IN[hOff(holder) + i];
  let rh = settle(hh, -kept); let rq = settle(IN[qOff(holder) + i], rh * IN[qOff(holder) + i] / hh); let rw = settle(IN[wOff(holder) + i], rh * IN[wOff(holder) + i] / hh);
  IN[hOff(holder) + i] = hh + rh; IN[qOff(holder) + i] += rq; IN[wOff(holder) + i] += rw;
  owedQ -= rq; owedW -= rw;
  var heatHolder = holder; var saltHolder = holder; var finestQ = grain(IN[qOff(holder) + i]); var finestW = grain(IN[wOff(holder) + i]);
  for (var k = ${fromK}; k < L; k++) {
    if (IN[hOff(k) + i] <= THINO) { continue; }
    let gq = grain(IN[qOff(k) + i]); let gw = grain(IN[wOff(k) + i]);
    if (gq < finestQ) { heatHolder = k; finestQ = gq; }
    if (gw < finestW) { saltHolder = k; finestW = gw; }
  }
  IN[qOff(heatHolder) + i] += settle(IN[qOff(heatHolder) + i], owedQ);
  IN[wOff(saltHolder) + i] += settle(IN[wOff(saltHolder) + i], owedW);
}`;
  return {
    head,
    oGradEta: `${K}  let e = ${idx}; if (e >= E) { return; }
  OD[O_GRADETA + e] = (OD[O_ETA + MI[COE + 2 * e + 1]] - OD[O_ETA + MI[COE + 2 * e]]) / MF[F_DC + e];
}`,
    oEdgeThickness: `${K}  let e = ${idx}; if (e >= E) { return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  OD[O_HEDGE + e] = 0.5 * (IN[hOff(0) + a] + IN[hOff(0) + b]);
  for (var k = 1; k < L; k++) { OD[O_HEDGE + k * E + e] = min(IN[hOff(k) + a], IN[hOff(k) + b]); }
  let sill = max(EPSO, OD[O_SILL + e] + 0.5 * (OD[O_ETA + a] + OD[O_ETA + b]));
  var above = 0.0;
  for (var k = 0; k < L; k++) { let hk = OD[O_HEDGE + k * E + e]; OD[O_HEDGE + k * E + e] = min(hk, max(0.0, sill - above)); above += hk; }
}`,
    /*
     * oCellTendency also prepares the cells for oMomentum:
     * each cell-layer's running sum of ρ'h above it and its density
     * anomaly as the adjacent pair PR (the sum of h above is re-summed
     * where it is read).
     */
    oCellTendency: `${K}  let i = ${idx}; if (i >= C) { return; }
  var edges: array<i32, MAXE>; var nbrs: array<i32, MAXE>;
  for (var m = 0; m < MAXE; m++) { edges[m] = MI[EOC + MAXE * i + m]; nbrs[m] = MI[COC + MAXE * i + m]; }
  var weighted = OD[O_LOAD + i] + OD[O_DRAFT + i] * eosAnomaly(OD[O_LOADT + i] + TREF, OD[O_LOADS + i] + SREF);
  for (var k = 0; k < L; k++) {
    let hv = IN[hOff(k) + i]; let Ti = select(LABEL_T[k] - TREF, IN[qOff(k) + i] / hv, hv > 1e-6); let Si = select(LABEL_S[k] - SREF, IN[wOff(k) + i] / hv, hv > 1e-6);
    let r = select(RHOA[k], eosAnomaly(IN[qOff(k) + i] / hv + TREF, IN[wOff(k) + i] / hv + SREF), hv > 1e-6);
    let n = k * C + i; OD[O_PR + 2 * n] = weighted; OD[O_PR + 2 * n + 1] = r; weighted += r * hv;
    var divH = 0.0; var divQ = 0.0; var divW = 0.0; var difQ = 0.0; var difW = 0.0;
    for (var m = 0; m < MAXE; m++) {
      let e = edges[m]; let j = nbrs[m];
      let f = f32(MI[ESC + MAXE * i + m]) * OD[O_FLUX + k * E + e] * MF[F_DV + e];
      divH += f;
      let hvj = IN[hOff(k) + j]; let Tj = select(LABEL_T[k] - TREF, IN[qOff(k) + j] / hvj, hvj > 1e-6); let Sj = select(LABEL_S[k] - SREF, IN[wOff(k) + j] / hvj, hvj > 1e-6);
      ${shear ? `if (OD[O_EMASK + e] > 0.5) {
        let dc = MF[F_DC + e]; var kappa = LATSPEED * dc;
        if (k < FIXEDL) {
          let v1 = MI[VOE + 2 * e]; let v2 = MI[VOE + 2 * e + 1];
          let zeta = 0.5 * ((OD[O_AVORT + k * V + v1] - MF[F_FV + v1]) + (OD[O_AVORT + k * V + v2] - MF[F_FV + v2]));
          kappa += min(SHEARCAP, SHEARDIFF * abs(zeta) * dc * dc);
        }
        let wd = kappa * OD[O_HEDGE + k * E + e] / dc * MF[F_DV + e];
        difQ += wd * (Tj - Ti); difW += wd * (Sj - Si);
      }` : ''}
      var Te = select(Tj, Ti, f > 0.0); var Se = select(Sj, Si, f > 0.0);
      ${third ? `if (f != 0.0) {
        let side = select(1, 0, OD[O_FLUX + k * E + e] > 0.0); let uu = i32(OD[O_BEY + 2 * e + side]); let huu = IN[hOff(k) + uu];
        if (OD[O_CMASK + uu] > 0.5 && huu > 1e-6 && select(hv, hvj, f > 0.0) > 1e-6) {
          let ratio = MF[F_DC + e] / MF[F_DC + i32(OD[O_BED + 2 * e + side])];
          Te = limitedValue(Te, select(Ti, Tj, f > 0.0), IN[qOff(k) + uu] / huu, ratio); Se = limitedValue(Se, select(Si, Sj, f > 0.0), IN[wOff(k) + uu] / huu, ratio);
        }
      }` : ''}
      divQ += f * Te; divW += f * Se;
    }
    let area = MF[F_AREA + i];
    var dh = -divH / area; var dQ = -divQ / area; var dW = -divW / area;
    ${shear ? 'dQ += difQ / area; dW += difW / area;' : ''}
    if (OD[O_CMASK + i] < 0.5) { dh = 0.0; dQ = 0.0; dW = 0.0; }
    OUT[hOff(k) + i] = dh; OUT[qOff(k) + i] = dQ; OUT[wOff(k) + i] = dW;
  }
}`,
    oVertexVort: `${K}  let v = ${idx}; if (v >= V) { return; }
  var edges: array<i32, 3>; var signs: array<f32, 3>; var dcs: array<f32, 3>;
  for (var m = 0; m < 3; m++) { let e = MI[EOV + 3 * v + m]; edges[m] = e; signs[m] = f32(MI[ESV + 3 * v + m]); dcs[m] = MF[F_DC + e]; }
  let atri = MF[F_ATRI + v]; let fv = MF[F_FV + v];
  for (var k = 0; k < L; k++) {
    var zeta = 0.0;
    for (var m = 0; m < 3; m++) { zeta += signs[m] * IN[uOff(k) + edges[m]] * dcs[m]; }
    OD[O_AVORT + k * V + v] = zeta / atri + fv;
  }
}`,
    oDivCurl: `${K}  let n = ${idx};
  let fromLap = P[1] > 0.5;
  if (n < C) {
    let i = n;
    var edges: array<i32, MAXE>; var signs: array<f32, MAXE>; var dvs: array<f32, MAXE>;
    for (var m = 0; m < MAXE; m++) { let e = MI[EOC + MAXE * i + m]; edges[m] = e; signs[m] = f32(MI[ESC + MAXE * i + m]); dvs[m] = MF[F_DV + e]; }
    let area = MF[F_AREA + i];
    if (fromLap) { for (var k = 0; k < L; k++) { var sum = 0.0; for (var m = 0; m < MAXE; m++) { let e = edges[m]; sum += signs[m] * OD[O_LAPA + k * E + e] * dvs[m]; } OD[O_DIVS + k * C + i] = sum / area; } } else { for (var k = 0; k < L; k++) { var sum = 0.0; for (var m = 0; m < MAXE; m++) { let e = edges[m]; sum += signs[m] * IN[uOff(k) + e] * dvs[m]; } OD[O_DIVS + k * C + i] = sum / area; } }
  } else if (n < C + V) {
    let v = n - C;
    var edges: array<i32, 3>; var signs: array<f32, 3>; var dcs: array<f32, 3>;
    for (var m = 0; m < 3; m++) { let e = MI[EOV + 3 * v + m]; edges[m] = e; signs[m] = f32(MI[ESV + 3 * v + m]); dcs[m] = MF[F_DC + e]; }
    let area = MF[F_ATRI + v];
    if (fromLap) { for (var k = 0; k < L; k++) { var sum = 0.0; for (var m = 0; m < 3; m++) { let e = edges[m]; sum += signs[m] * OD[O_LAPA + k * E + e] * dcs[m]; } OD[O_CURLS + k * V + v] = sum / area; } } else { for (var k = 0; k < L; k++) { var sum = 0.0; for (var m = 0; m < 3; m++) { let e = edges[m]; sum += signs[m] * IN[uOff(k) + e] * dcs[m]; } OD[O_CURLS + k * V + v] = sum / area; } }
  }
}`,
    oLapVelocity: `${K}  let e = ${idx}; if (e >= E) { return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1]; let va = MI[VOE + 2 * e]; let vb = MI[VOE + 2 * e + 1];
  let dc = MF[F_DC + e]; let dv = MF[F_DV + e];
  let dst = select(O_LAPA, O_LAPB, P[1] > 0.5);
  for (var k = 0; k < L; k++) {
    let lap = (OD[O_DIVS + k * C + b] - OD[O_DIVS + k * C + a]) / dc - (OD[O_CURLS + k * V + vb] - OD[O_CURLS + k * V + va]) / dv;
    OD[dst + k * E + e] = lap;
  }
}`,
    /*
     * The eddy transport of ocean/layered.module.js, one edge a thread,
     * leaving each class's volume, heat and salt fluxes in the EDDYV, EDDYQ
     * and EDDYW scratch (see OD's layout); P[0] is the step and P[1] the
     * stability limit on κ. Each class moves whole multiples of the coarser
     * cell's spacing, measured in that cell, and the thinnest class thick
     * in both cells takes back what the rounding left of the edge's net
     * volume. oEddyApply rounds each class's change to its own spacing,
     * leaving what it did not record in the trial state OUT, free after the
     * combine, and oEddyLeftover gives each column's leftovers to its finest
     * stores, so the column keeps its thickness to its finest layer's
     * rounding and gains the heat and salt its edges carried.
     */
    oEddyFlux: `${constLine('EDDYTOP', o.eddyTaperDepth)} ${constLine('EDDYBOTTOM', EDDY_BOTTOM_TAPER)} ${eddyShared}${K}  let e = ${idx}; if (e >= E) { return; }
  for (var k = 0; k < L; k++) { OD[O_EDDYV + k * E + e] = 0.0; OD[O_EDDYQ + k * E + e] = 0.0; OD[O_EDDYW + k * E + e] = 0.0; }
  let kappa = min(OD[O_EDDYK + e], P[1]);
  if (OD[O_EMASK + e] < 0.5 || !(kappa > 0.0)) { return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  var overlap = 0.0;
  for (var k = 1; k < L; k++) { overlap += max(0.0, min(IN[hOff(k) + a], IN[hOff(k) + b]) - EPSO); }
  let sill = OD[O_SILL + e] + 0.5 * (OD[O_ETA + a] + OD[O_ETA + b]);
  let coefficient = kappa * MF[F_DV + e] / MF[F_DC + e];
  var za = IN[hOff(0) + a]; var zb = IN[hOff(0) + b];
  var aboveA = 0.0; var aboveB = 0.0; var overlapAbove = 0.0; var upper = 0.0;
  var flux: array<f32, ${L}>; var share: array<f32, ${L}>;
  for (var k = 1; k < L; k++) {
    let ha = IN[hOff(k) + a]; let hb = IN[hOff(k) + b]; let both = max(0.0, min(ha, hb) - EPSO);
    za += ha; zb += hb;
    aboveA += max(0.0, ha - EPSO); aboveB += max(0.0, hb - EPSO); overlapAbove += both;
    var lower = 0.0;
    if (k < L - 1) {
      let top = clamp(min(za, zb) / EDDYTOP, 0.0, 1.0);
      let bottom = clamp(min(sill - max(za, zb), overlap - overlapAbove) / EDDYBOTTOM, 0.0, 1.0);
      let ceiling = clamp(max(aboveA, aboveB) / THINO, 0.0, 1.0);
      lower = -coefficient * top * bottom * ceiling * (zb - za);
      lower *= min(modeWeight(k + 1, a), modeWeight(k + 1, b));
    }
    let weight = clamp((max(ha, hb) - EPSO) / THINO, 0.0, 1.0);
    flux[k] = weight * (lower - upper);
    share[k] = weight * both;
    upper = lower;
  }
  var residual = 0.0; var carriers = 0.0;
  for (var round = 0; ; round++) {
    residual = 0.0; carriers = 0.0;
    for (var k = 1; k < L; k++) { residual += flux[k]; carriers += share[k]; }
    if (round == 3 || !(carriers > THINO)) { break; }
    var held = false;
    for (var k = 1; k < L; k++) {
      let f = flux[k] - residual * share[k] / carriers;
      let allowed = eddyAllowance(k, select(b, a, f > 0.0));
      if (abs(f) > allowed) { flux[k] = select(-allowed, allowed, f > 0.0); share[k] = 0.0; held = true; }
    }
    if (!held) { break; }
  }
  if (!(carriers > THINO)) { return; }
  var scale = 1.0;
  for (var k = 1; k < L; k++) {
    let f = flux[k] - residual * share[k] / carriers;
    flux[k] = f;
    let allowed = eddyAllowance(k, select(b, a, f > 0.0));
    if (abs(f) > allowed) { scale = min(scale, allowed / abs(f)); }
  }
  let sa = MF[F_AREA + a]; let sb = MF[F_AREA + b];
  var net = 0.0; var taker = 0; var finest = 0.0;
  for (var k = 1; k < L; k++) {
    let ha = IN[hOff(k) + a]; let hb = IN[hOff(k) + b];
    if (flux[k] != 0.0) { flux[k] = eddyShare(ha, hb, sa, sb, scale * flux[k] * P[0]); net += flux[k]; }
    if (min(ha, hb) > THINO) { let g = grain(max(ha, hb)); if (taker == 0 || g < finest) { taker = k; finest = g; } }
  }
  if (taker > 0) { flux[taker] -= net; }
  for (var k = 1; k < L; k++) {
    let moved = flux[k];
    if (moved == 0.0) { continue; }
    let qa = IN[qOff(k) + a]; let qb = IN[qOff(k) + b]; let wa = IN[wOff(k) + a]; let wb = IN[wOff(k) + b];
    let toB = moved > 0.0; let carried = moved / IN[hOff(k) + select(b, a, toB)];
    OD[O_EDDYV + k * E + e] = moved / P[0];
    OD[O_EDDYQ + k * E + e] = eddyShare(qa, qb, sa, sb, carried * select(qb, qa, toB)) / P[0];
    OD[O_EDDYW + k * E + e] = eddyShare(wa, wb, sa, sb, carried * select(wb, wa, toB)) / P[0];
  }
}`,
    oEddyApply: eddyApply(1),
    oEddyLeftover: eddyLeftover(1),
    oAdvance: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>, @builtin(num_workgroups) groups: vec3<u32>) {
  let stride = i32(groups.x * groups.y) * ${WORKGROUP};
  for (var n = ${idx}; n < OSTOTAL; n += stride) { OUT[n] = IN[n] + P[0] * OD[n]; }
}`,
    oCombine: `@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>, @builtin(num_workgroups) groups: vec3<u32>) {
  let stride = i32(groups.x * groups.y) * ${WORKGROUP};
  for (var n = ${idx}; n < OSTOTAL; n += stride) { IN[n] += P[0] * (OUT[n] + 2.0 * OD[n] + 2.0 * PH[n] + S[n]); }
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
  for (var s = 0; s < MAXEE; s++) { let slot = MAXEE * e + s; let other = MI[EOE + slot]; sum += MF[F_PVW + slot] * OD[B_BCUR + C + other] * 0.5 * (OD[O_FEDGE + e] + OD[O_FEDGE + other]); }
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
     * bake the varying offsets into the shader at build time. Each stage
     * reads its trial state across neighbors but writes only its own
     * cell and edge, so the trial alternates between two blocks and no
     * stage writes the block it reads.
     */
    ...Object.fromEntries([1, 2, 3, 4].map((stage) => {
      const inBase = ['B_BCUR', 'B_BTRIALA', 'B_BTRIALB', 'B_BTRIALA'][stage - 1];
      const finish = (m) => (stage < 4 ? `OD[B_BK${stage} + ${m}] = t;
    OD[${stage === 2 ? 'B_BTRIALB' : 'B_BTRIALA'} + ${m}] = OD[B_BCUR + ${m}] + ${stage === 3 ? 'P[1]' : 'P[0]'} * t;` : `let next = OD[B_BCUR + ${m}] + P[2] * (OD[B_BK1 + ${m}] + 2.0 * OD[B_BK2 + ${m}] + 2.0 * OD[B_BK3 + ${m}] + t);
    OD[B_BCUR + ${m}] = next;
    OD[B_BAVG + ${m}] += next * P[3];`);
      return [`oBarotropicStage${stage}`, `${K}  let n = ${idx};
  if (n < C) {
    let i = n;
    var sum = 0.0;
    for (var m = 0; m < MAXE; m++) { let e = MI[EOC + MAXE * i + m]; sum += f32(MI[ESC + MAXE * i + m]) * OD[${inBase} + C + e] * MF[F_DV + e]; }
    let t = select(0.0, -(sum / MF[F_AREA + i]), OD[O_CMASK + i] > 0.5);
    ${finish('i')}
  }
  if (n < E) {
    let e = n;
    var t = 0.0;
    if (OD[O_EMASK + e] >= 0.5) {
      let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
      let gradEtaB = (OD[${inBase} + b] - OD[${inBase} + a]) / MF[F_DC + e];
      var cor = 0.0;
      for (var s = 0; s < MAXEE; s++) { let slot = MAXEE * e + s; let other = MI[EOE + slot]; cor += MF[F_PVW + slot] * OD[${inBase} + C + other] * 0.5 * (OD[O_FEDGE + e] + OD[O_FEDGE + other]); }
      t = -OGRAV * OD[O_DEPTHEDGE + e] * gradEtaB + cor / MF[F_DC + e] + OD[O_SLOW + e];
    }
    ${finish('C + e')}
  }
}`];
    })),
    oVelocityShiftClamp: `${K}  let e = ${idx}; if (e >= E) { return; }
  if (OD[O_EMASK + e] < 0.5) { return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  var Ha = IN[hOff(0) + a]; var Hb = IN[hOff(0) + b];
  for (var k = 1; k < L; k++) {
    let ha = IN[hOff(k) + a]; let hb = IN[hOff(k) + b]; let he = OD[O_HEDGE + k * E + e];
    if (he < THINEDGE || (he <= TOKENEDGE[k] && min(Ha + ha, Hb + hb) <= max(Ha, Hb))) { IN[uOff(k) + e] = IN[uOff(k - 1) + e]; }
    Ha += ha; Hb += hb;
  }
  var sumH = 0.0; var transport = 0.0;
  for (var k = 0; k < L; k++) { let he = OD[O_HEDGE + k * E + e]; sumH += he; transport += he * IN[uOff(k) + e]; }
  let shift = (OD[B_BAVG + C + e] - transport) / max(sumH, EPSO);
  for (var k = 0; k < L; k++) { let n = uOff(k) + e; IN[n] = clamp(IN[n] + shift, -SPEEDLIM, SPEEDLIM); }
}`,
    oReadSurface: `${K}  let i = ${idx}; if (i >= C) { return; }
  if (OD[O_SHELF + i] > 0.5) { OD[O_ICED + i] = 0.0; return; }
  let iced = OD[O_SURFICE + i] > 0.0;
  OD[O_ICED + i] = select(0.0, 1.0, iced);
${o.drift ? `  let arrived = OD[O_ARRIVED + i] > 0.5;
  let t0 = select(OD[O_SURFT + i], select(FREEZE, OD[O_WATERT + i], arrived), iced);
  // Open water, and ice that arrived over it, adds the surface's change since PREVT0: reading T0 = Q/h back whole gains heat under fast-math division.
  if (iced && !arrived) { IN[qOff(0) + i] = IN[hOff(0) + i] * (t0 - TREF); } else { IN[qOff(0) + i] += IN[hOff(0) + i] * (t0 - OD[O_PREVT0 + i]); }
` : `  let t0 = select(OD[O_SURFT + i], FREEZE, iced);
  // Open water adds the surface's change since PREVT0: reading T0 = Q/h back whole gains heat under fast-math division.
  if (iced) { IN[qOff(0) + i] = IN[hOff(0) + i] * (t0 - TREF); } else { IN[qOff(0) + i] += IN[hOff(0) + i] * (t0 - OD[O_PREVT0 + i]); }
`}  let light = PH[PH_PEN + i];
  if (light == 0.0) { return; }
  let taken = PH[PH_TRANS + i];
  var given = 0.0; var last = 0;
  if (taken > 0.0) {
    var z = IN[hOff(0) + i]; var top = transmitted(z); var carry = 0.0;
    for (var k = 1; k < L; k++) {
      let hk = IN[hOff(k) + i];
      z += hk;
      let bottom = transmitted(z);
      let want = carry + light * (top - bottom) / taken;
      if (hk > 0.0) {
        let dq = settle(IN[qOff(k) + i], want / RHOCP);
        IN[qOff(k) + i] += dq;
        given += RHOCP * dq; carry = want - RHOCP * dq; last = k;
      } else { carry = want; }
      top = bottom;
      if (z >= PENDEPTH && last > 0) { break; }
    }
  }
  let rest = settle(IN[qOff(last) + i], (light - given) / RHOCP);
  IN[qOff(last) + i] += rest;
  PH[PH_PEN + i] = light - given - RHOCP * rest;
}`,
    oStressFromAtmosphere: o.drift ? driftStressKernel(K, idx) : `${K}  let e = ${idx}; if (e >= E) { return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  if (OD[O_EMASK + e] < 0.5) { OD[O_STRESS + e] = 0.0; return; }
  let coverA = select(0.0, select(1.0, PH[PH_CONC + a], PH[PH_CONC + a] > 0.0), OD[O_ICED + a] > 0.5);
  let coverB = select(0.0, select(1.0, PH[PH_CONC + b], PH[PH_CONC + b] > 0.0), OD[O_ICED + b] > 0.5);
  let through = (1.0 - 0.5 * (coverA + coverB) * (1.0 - TRANSMIT)) * (1.0 - 0.5 * (OD[O_SHELF + a] + OD[O_SHELF + b]));
${o.implicitStress ? '  if (PH[PH_STRESSOK] > 0.5) { OD[O_STRESS + e] = through * PH[PH_STRESS + e]; return; }\n' : ''}  let bottom = (K - 1) * C;
  let rhoA = S[S_PI + a] * LV[L_SM + K - 1] / (RGAS * S[S_TH + bottom + a] * D[D_EXM + bottom + a]);
  let rhoB = S[S_PI + b] * LV[L_SM + K - 1] / (RGAS * S[S_TH + bottom + b] * D[D_EXM + bottom + b]);
  let fa = PH[PH_DRAG + a] * rhoA * max(D[D_WIND + a], GUSTO); let fb = PH[PH_DRAG + b] * rhoB * max(D[D_WIND + b], GUSTO);
  OD[O_STRESS + e] = through * 0.5 * (fa + fb) * S[S_U + (K - 1) * E + e];
}`,
    oShelfBudget: `var<workgroup> POOLED: array<vec4<f32>, ${WORKGROUP}>;
${K}  let t = i32(id.x);
  var sums = vec4<f32>(0.0, 0.0, 0.0, 0.0);
  for (var i = t; i < C; i += ${WORKGROUP}) {
    let area = MF[F_AREA + i]; let south = OD[O_FCELL + i] < 0.0;
    if (OD[O_SHELF + i] > 0.5) { let m = OD[O_MELT + i] * area; sums += select(vec4<f32>(0.0, 0.0, 0.0, m), vec4<f32>(0.0, 0.0, m, 0.0), south); }
    else if (OD[O_CMASK + i] > 0.5) { let d = OD[O_CALVED + i] * area; sums += select(vec4<f32>(0.0, d, 0.0, 0.0), vec4<f32>(d, 0.0, 0.0, 0.0), south); }
  }
  POOLED[t] = sums;
  workgroupBarrier();
  for (var half = ${WORKGROUP / 2}; half > 0; half = half / 2) {
    if (t < half) { POOLED[t] += POOLED[t + half]; }
    workgroupBarrier();
  }
  if (t == 0) {
    for (var p = 0; p < 2; p++) {
      let d = POOLED[0][p]; let m = POOLED[0][2 + p]; let moved = select(0.0, min(m, d), d > 0.0);
      OD[O_POOLS + p] = d; OD[O_POOLS + 2 + p] = m; OD[O_POOLS + 4 + p] = moved; OD[O_POOLS + 6 + p] = m - moved;
      OD[O_POOLS + 8 + p] = select(1.0, (d - moved) / d, d > 0.0); OD[O_POOLS + 10 + p] = select(0.0, (m - moved) / m, m != 0.0);
    }
  }
}`,
    oSalt: `${K}  let i = ${idx}; if (i >= C) { return; }
  if (OD[O_CMASK + i] < 0.5) { return; }
  let pool = select(1, 0, OD[O_FCELL + i] < 0.0);
  if (OD[O_SHELF + i] > 0.5) { OD[O_DRAW + i] += OD[O_POOLS + 10 + pool] * OD[O_MELT + i]; return; }
  let keep = OD[O_POOLS + 8 + pool]; let calved = OD[O_CALVED + i];
  let s = IN[wOff(0) + i] / IN[hOff(0) + i] + SREF;
  var w = IN[wOff(0) + i] + s * (OD[O_FRESH + i] + (1.0 - keep) * calved) / 1000.0;
  let volume = iceVolume(i);
  let grown = (volume - OD[O_PREVVOL + i]) * ICEDENS / 1000.0;
  w += (s - ICESAL) * grown;
  w = max(-SREF * IN[hOff(0) + i], w);
  IN[wOff(0) + i] = w;
  OD[O_FRESH + i] = 0.0;
  OD[O_PREVVOL + i] = volume;
  IN[qOff(0) + i] -= OLFUS * keep * calved / RHOCP;
  OD[O_CALVED + i] = 0.0;
}`,
    oWriteSurface: `${K}  let i = ${idx}; if (i >= C) { return; }
  if (OD[O_CMASK + i] < 0.5) { PH[PH_OFLUX + i] = 0.0; return; }
  if (OD[O_SHELF + i] > 0.5) {
    PH[PH_OFLUX + i] = 0.0;
    OD[O_S0 + i] = IN[wOff(0) + i] / IN[hOff(0) + i] + SREF;
    OD[O_T0 + i] = IN[qOff(0) + i] / IN[hOff(0) + i] + TREF;
    OD[O_PREVT0 + i] = OD[O_T0 + i];
    return;
  }
  PH[PH_CAP + i] = RHOCP * max(IN[hOff(0) + i], 1.0);
${o.penetration ? '  PH[PH_TRANS + i] = transmitted(IN[hOff(0) + i]);\n' : ''}  OD[O_S0 + i] = IN[wOff(0) + i] / IN[hOff(0) + i] + SREF;
  if (OD[O_ICED + i] > 0.5) {
    PH[PH_OFLUX + i] = RHOCP * (IN[qOff(0) + i] - IN[hOff(0) + i] * (FREEZE - TREF)) / P[6];
    IN[qOff(0) + i] = IN[hOff(0) + i] * (FREEZE - TREF);
    OD[O_T0 + i] = FREEZE;
  } else {
    OD[O_T0 + i] = IN[qOff(0) + i] / IN[hOff(0) + i] + TREF;
    S[S_TS + i] = OD[O_T0 + i];
    PH[PH_OFLUX + i] = 0.0;
  }
  OD[O_PREVT0 + i] = OD[O_T0 + i];
}`,
    oAccumulateFresh: `${K}  let i = ${idx}; if (i >= C) { return; }
  let rain = PH[PH_RAIN + i]; let seen = OD[O_RAINSEEN + i];
  let delta = select(rain, rain - seen, rain >= seen);
  if (OD[O_CMASK + i] > 0.5 && OD[O_SHELF + i] < 0.5) {
    var fresh = PH[PH_EVAP + i] * P[6] - delta;
    var calved = 0.0;
    for (var p = i32(OD[O_DRAINSTART + i]); p < i32(OD[O_DRAINSTART + i + 1]); p++) {
      let j = i32(OD[O_DRAINCELL + p]);
      let runoff = PH[PH_RUNOFF + j]; let ran = OD[O_RUNOFFSEEN + j];
      fresh -= select(runoff, runoff - ran, runoff >= ran) * OD[O_DRAINW + p];
      OD[O_RUNOFFSEEN + j] = runoff;
      let shed = PH[PH_DISCHARGE + j]; let taken = OD[O_DISCHARGESEEN + j];
      calved += select(shed, shed - taken, shed >= taken) * OD[O_DRAINW + p];
      OD[O_DISCHARGESEEN + j] = shed;
${o.lakes ? LAKE_DRAIN_WGSL : ''}    }
    OD[O_FRESH + i] += fresh;
    OD[O_CALVED + i] += calved;
  }
  OD[O_RAINSEEN + i] = rain;
}`,
    ...hybridKernels(o, K, idx),
    ...(o.drift ? iceDriftKernels({ K, idx, options: o.drift.options, third }) : {}),
  };
}

/*
 * The global sums behind the ocean diagnostics: [name, how the
 * workgroup partials combine, the cell's term]. Edge terms are taken by
 * the edge's first cell.
 */
const HYBRID_REDUCED = [
  ['area', 'sum', 'wet * a'], ['depth', 'sum', 'wet * a * OD[O_MIXED + i]'], ['salinity', 'sum', 'wet * a * (IN[wOff(0) + i] / max(EPSO, h0) + SREF)'],
  ['thermocline', 'sum', 'wet * a * thermo'], ['heat', 'sum', 'wet * a * RHOCP * heat'],
  ['interiorT', 'sum', 'wet * a * interiorQ'], ['interiorH', 'sum', 'wet * a * interiorH'],
  ['ssh', 'max', 'wet * abs(OD[O_ETA + i])'], ['speed', 'max', 'speed'], ['limited', 'sum', 'limited'], ['transport', 'max', 'transport'],
  ['dried', 'sum', 'wet * OD[O_DRIED + i]'],
];
// The depth at which the column's density, its layer means held from decreasing with depth, reaches THERMOCLINE_DENSITY.
const HYBRID_THERMOCLINE = `const THERMOA: f32 = ${(THERMOCLINE_DENSITY - SEAWATER.rho0).toExponential(10)};
fn thermoclineDepthOf(i: i32) -> f32 {
  var z = 0.0; var previous = -1.0e30; var center = 0.0;
  for (var k = 0; k < L; k++) {
    let hv = IN[hOff(k) + i];
    let r = select(RHOA[k], eosAnomaly(IN[qOff(k) + i] / hv + TREF, IN[wOff(k) + i] / hv + SREF), hv > 1e-6);
    let rk = max(previous, r); let c = z + 0.5 * hv;
    if (rk >= THERMOA) { return select(center + (THERMOA - previous) / max(rk - previous, 1e-12) * (c - center), 0.0, k == 0); }
    previous = rk; center = c; z += hv;
  }
  return z;
}
`;
const hybridReducedSetup = `    let a = MF[F_AREA + i];
    let wet = select(0.0, 1.0, OD[O_CMASK + i] > 0.5);
    let h0 = IN[hOff(0) + i];
    let thermo = select(0.0, thermoclineDepthOf(i), wet > 0.5);
    var heat = 0.0; var interiorQ = 0.0; var interiorH = 0.0;
    for (var k = 0; k < L; k++) {
      heat += IN[qOff(k) + i] + TREF * IN[hOff(k) + i];
      if (k > 0) { interiorQ += IN[qOff(k) + i] + TREF * IN[hOff(k) + i]; interiorH += IN[hOff(k) + i]; }
    }
    var speed = 0.0; var limited = 0.0; var transport = 0.0;
    for (var m = 0; m < MAXE; m++) {
      let e = MI[EOC + MAXE * i + m];
      if (MI[ESC + MAXE * i + m] == 0 || MI[COE + 2 * e] != i) { continue; }
      let u0 = IN[uOff(0) + e];
      speed = max(speed, abs(u0));
      var clamped = false;
      for (var k = 0; k < L; k++) { clamped = clamped || abs(abs(IN[uOff(k) + e]) - SPEEDLIM) < 1e-4; }
      limited += select(0.0, 1.0, clamped);
      let ca = MI[COE + 2 * e]; let cb = MI[COE + 2 * e + 1];
      let sill = max(EPSO, OD[O_SILL + e] + 0.5 * (OD[O_ETA + ca] + OD[O_ETA + cb]));
      var above = 0.0; var flow = 0.0;
      for (var k = 0; k < L; k++) { let hk = select(min(IN[hOff(k) + ca], IN[hOff(k) + cb]), 0.5 * (IN[hOff(0) + ca] + IN[hOff(0) + cb]), k == 0); flow += min(hk, max(0.0, sill - above)) * IN[uOff(k) + e]; above += hk; }
      transport = max(transport, abs(flow) * MF[F_DV + e]);
    }`;

export function createLayeredOcean(core, { registers = null, ...options } = {}) {
  const o = { ...OCEAN_DEFAULTS, ...options };
  if (o.vertical !== undefined && o.vertical !== 'hybrid') throw new Error(`the ocean's vertical coordinate is 'hybrid'; a saved slab ocean is converted on upload, not run (vertical ${o.vertical})`);
  if (!GENERATORS.includes(o.generator)) throw new Error(`the generator is one of ${GENERATORS.join(', ')}, not ${o.generator}`);
  const table = o.layerTable ?? HYBRID_LAYERS;
  if (!TRACER_ADVECTIONS.includes(o.tracerAdvection)) throw new Error(`tracerAdvection is one of ${TRACER_ADVECTIONS.join(', ')}, not ${o.tracerAdvection}`);
  if (!EDGE_TRANSPORTS.includes(o.edgeTransport)) throw new Error(`edgeTransport is one of ${EDGE_TRANSPORTS.join(', ')}, not ${o.edgeTransport}`);
  if (!ROTATION_LIMITS.includes(o.rotationLimit)) throw new Error(`rotationLimit is one of ${ROTATION_LIMITS.join(', ')}, not ${o.rotationLimit}`);
  if (!LANGMUIR_ICE.includes(o.langmuirIce)) throw new Error(`langmuirIce is one of ${LANGMUIR_ICE.join(', ')}, not ${o.langmuirIce}`);
  const given = !!o.iceDrift && typeof o.iceDrift === 'object' && o.iceDrift.airStress === 'given';
  const driftOptions = o.iceDrift ? iceDriftOptions(typeof o.iceDrift === 'object' ? Object.fromEntries(Object.entries(o.iceDrift).filter(([name]) => name !== 'airStress')) : o.iceDrift) : null;
  const third = o.tracerAdvection === 'third';
  if (third && !(o.shearDiffusivity > 0)) throw new Error('tracerAdvection third needs shearDiffusivity above zero');
  if (o.penetration !== null && !JERLOV[o.penetration]) throw new Error(`penetration is one of ${Object.keys(JERLOV).join(', ')} or null, not ${o.penetration}`);
  if (o.langmuir !== 0 && !(o.langmuir >= 1)) throw new Error(`langmuir caps a factor of at least 1, so it is 0 or at least 1, not ${o.langmuir}`);
  const { device, buffers, mesh, meshSpacing } = core;
  const C = core.C, E = core.E, V = core.V;
  const L = table.K;
  const rho = Array.from(table.densities);
  const labels = rho;
  const labelS = rho.map((r) => { const j = LAYER_DENSITIES.indexOf(r); return j >= 0 ? LAYER_SALINITIES[j] : o.referenceS; });
  const labelT = rho.map((r, k) => Math.max(FREEZING_POINT, labelTemperature(r, labelS[k])));
  const nu4 = closureCoefficient(meshSpacing, o.closureHours, o.closureSpacing);
  const eddyKappa = eddyDiffusivities(mesh, o.eddyDiffusivity), eddyLimit = o.eddyDiffusivity > 0 ? eddyDiffusionLimit(mesh) : 0;
  let minSpacing = Infinity;
  for (let e = 0; e < E; e++) minSpacing = Math.min(minSpacing, mesh.dcEdge[e]);
  const geography = o.geography || null, { water, cellOcean, shelf, cellOpen, draft } = oceanCells(C, geography);
  const edgeOcean = water ? water.edgeOcean : new Uint8Array(E).fill(1);
  const bed = new Float64Array(C);
  const smoothed = o.bathymetry ? null : bathymetryFrom(mesh, water, { minimumDepth: o.minimumDepth, flatDepth: o.flatDepth });
  for (let i = 0; i < C; i++) bed[i] = !cellOcean[i] ? 0 : o.bathymetry ? o.bathymetry[i] : smoothed[i];
  const { D, sillEdge } = restingColumns(mesh, bed, draft, edgeOcean, o.edgeSill === undefined ? (water ? water.edgeSill ?? null : null) : o.edgeSill, o.sillFloor ?? o.minimumDepth);
  const hasShelves = shelf.some((s) => s), shelfLoad = new Float64Array(C), shelfMelted = new Float64Array(C), shelfDraw = new Float64Array(C), meltOn = o.shelfMelt && hasShelves;
  let deepest = 0;
  for (let i = 0; i < C; i++) deepest = Math.max(deepest, D[i]);
  const substepLimit = 0.35 * minSpacing / Math.sqrt(o.gravity * Math.max(deepest, 1));
  const inertialLimit = inertialStepLimit(mesh);
  const salinityProfile = o.salinityProfile || defaultSalinityProfile;
  let seaCover = null;

  const OS = seq([['OH', L * C], ['OU', L * E], ['OQ', L * C], ['OW', L * C]]);
  const outlet = runoffOutlets(mesh, geography);
  const drainStart = new Float32Array(C + 1), drainCells = [], drainWeights = [];
  {
    const byOutlet = Array.from({ length: C }, () => []);
    for (let j = 0; j < C; j++) if (!cellOpen[j] && outlet[j] >= 0) byOutlet[outlet[j]].push(j);
    for (let i = 0; i < C; i++) {
      drainStart[i] = drainCells.length;
      for (const j of byOutlet[i]) { drainCells.push(j); drainWeights.push(mesh.areaCell[j] / mesh.areaCell[i]); }
    }
    drainStart[C] = drainCells.length;
  }
  /*
   * The L-sized scratch lives within one tendency (tendency below, in
   * dispatch order): FLUX from oEdgeFluxPV to oCellTendency, LAPA from
   * the first oLapVelocity to the second oDivCurl and LAPB from the
   * second oLapVelocity to oMomentum share a slot with the regrid's
   * TARGET and the column physics' VISC; AVORT from oVertexVort to
   * oEdgeFluxPV, CURLS from each oDivCurl to the oLapVelocity after it
   * and PHI from oKineticPhi, dispatched just before oMomentum, to
   * oMomentum share another. After the combine the eddy passes leave
   * their fluxes in EDDYV, EDDYQ and EDDYW over FLUX, HEDGE and QE, which
   * oEdgeThickness and every tendency recompute before reading them.
   */
  const OD = seq([
    [['FLUX', 'LAPA', 'LAPB', 'EDDYV', 'TARGET', 'VISC'], L * E],
    [['HEDGE', 'EDDYQ'], L * E], [['AVORT', 'CURLS', 'PHI'], L * V], [['QE', 'EDDYW'], L * E], ['DIVS', L * C],
    ['BUOY', C], ['GRADETA', E], ['SLOW', E], ['DEPTHEDGE', E],
    ['ETA', C], ['FRESH', C], ['PREVT0', C], ['PREVVOL', C], ['ICED', C], ['STRESS', E],
    ['SURFT', C], ['SURFICE', C], ['SURFCONC', C], ['T0', C], ['S0', C], ['RAINSEEN', C],
    ['EMASK', E], ['CMASK', C], ['BATH', C], ['FEDGE', E],
    ['RUNOFFSEEN', C], ['DISCHARGESEEN', C], ['CALVED', C], ['DRAINSTART', C + 1], ['DRAINCELL', Math.max(1, drainCells.length)], ['DRAINW', Math.max(1, drainCells.length)],
    ['EDDYK', E],
    ['PR', 2 * L * C], ['MIXED', C], ['BLD', C], ['DRIED', C], ['FCELL', C], ['SILL', E], ['SHELF', C], ['DRAFT', C], ['LOAD', C], ['LOADT', C], ['LOADS', C], ['LIQUIDUS', C], ['MELT', C], ['MELTSUM', C], ['DRAW', C], ['POOLS', 12],
    ...(third ? [['BEY', 2 * E], ['BED', 2 * E]] : []),
    ...(driftOptions ? driftSections(C, E, given) : []),
  ]);
  // The barotropic RK4's fixed-offset blocks index into the same OD storage
  // array, so their offsets must continue on from OD's own sections rather
  // than starting a second, colliding index space at 0.
  const B = seq([['BCUR', C + E], ['BK1', C + E], ['BK2', C + E], ['BK3', C + E], ['BTRIALA', C + E], ['BTRIALB', C + E], ['BAVG', C + E]]);
  const bTotal = B.total;
  for (const k of Object.keys(B)) if (k !== 'total') B[k] += OD.total;
  B.total = bTotal;
  const ODTOTAL = OD.total + bTotal;

  const implicitStress = !!core.physics && core.physics.surfaceExchange === 'roughness' && !!core.physics.implicitDrag;
  const iceThermo = { skinHeatCapacity: core.physics?.skinHeatCapacity ?? 2e5, iceDensity: core.physics?.iceDensity ?? 917, latentHeatFusion: core.physics?.latentHeatFusion ?? 3.34e5, freshSnowAlbedo: core.physics?.freshSnowAlbedo ?? SNOW_AGEING.freshSnowAlbedo };
  const driftPlan = driftOptions ? { options: driftOptions, waterDensity: o.density, thermo: iceThermo, minimumConcentration: MINIMUM_CONCENTRATION, removal: seaIceRemoval(!!core.physics?.keepSmallIce), remnant: ICE_REMNANT, freezing: FREEZING_POINT, implicitStress, given } : null;
  const kernels = oceanKernels({ ...o, table, L, C, E, V, OS, OD, B, rho, labelT, labelS, nu4, omega: mesh.omega, implicitStress, drift: driftPlan, lakes: core.layout.PH.LAKE !== undefined });
  const reduced = HYBRID_REDUCED;
  const OF = seq([['SST', C], ['SSS', C], ['H1', C], ['THD', C], ['ETA', C], ['CUR', 3 * C], ['CSPD', C], ['UPW', C], ['PART', reduced.length * reductionGroups(C)]]);
  kernels.oFrame = `${HYBRID_THERMOCLINE}@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= C) { return; }
  let h0 = max(EPSO, IN[hOff(0) + i]);
  OUT[${OF.SSS} + i] = IN[wOff(0) + i] / h0 + SREF;
  OUT[${OF.H1} + i] = OD[O_MIXED + i];
  OUT[${OF.THD} + i] = thermoclineDepthOf(i);
  OUT[${OF.ETA} + i] = OD[O_ETA + i];
  let depth = P[4]; let upwellingDepth = P[5];
  var layer = -1; var top = 0.0;
  for (var k = 0; k < L; k++) { let hk = IN[hOff(k) + i]; if (layer < 0 && depth < top + hk) { layer = k; } if (layer < 0) { top += hk; } }
  if (layer < 0) {
    OUT[${OF.SST} + i] = DRY; OUT[${OF.UPW} + i] = DRY;
    OUT[${OF.CUR} + 3 * i] = 0.0; OUT[${OF.CUR} + 3 * i + 1] = 0.0; OUT[${OF.CUR} + 3 * i + 2] = 0.0; OUT[${OF.CSPD} + i] = DRY;
    return;
  }
  OUT[${OF.SST} + i] = IN[qOff(layer) + i] / max(EPSO, IN[hOff(layer) + i]) + TREF;
  var w = vec3<f32>(0.0, 0.0, 0.0); var upwelling = 0.0;
  for (var m = 0; m < MAXE; m++) {
    let e = MI[EOC + MAXE * i + m]; let s = f32(MI[ESC + MAXE * i + m]);
    w += abs(s) * 0.5 * MF[F_DC + e] * MF[F_DV + e] * IN[uOff(layer) + e] * vec3<f32>(MF[F_NEDGE + 3 * e], MF[F_NEDGE + 3 * e + 1], MF[F_NEDGE + 3 * e + 2]);
    let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
    var above = 0.0; var transport = 0.0;
    for (var k = 0; k < L; k++) {
      if (above >= upwellingDepth) { break; }
      let he = 0.5 * (IN[hOff(k) + a] + IN[hOff(k) + b]);
      transport += IN[uOff(k) + e] * min(he, upwellingDepth - above);
      above += he;
    }
    upwelling += s * MF[F_DV + e] * transport;
  }
  w = w / MF[F_AREA + i];
  OUT[${OF.CUR} + 3 * i] = w.x; OUT[${OF.CUR} + 3 * i + 1] = w.y; OUT[${OF.CUR} + 3 * i + 2] = w.z;
  OUT[${OF.CSPD} + i] = length(w);
  OUT[${OF.UPW} + i] = upwelling / MF[F_AREA + i];
}`;
  kernels.oReduce = HYBRID_THERMOCLINE + reductionKernel(HYBRID_REDUCED, { count: C, base: OF.PART, setup: hybridReducedSetup });
  const ob = {}, own = (name, size) => { const b = emptyBuffer(device, size); b.label = 'layeredOcean' + name; return b; };
  ob.S = own('S', 4 * OS.total);
  for (const name of ['T', 'K1', 'K2', 'K3', 'K4']) ob[name] = registers && registers[name] && registers[name].size >= 4 * OS.total ? registers[name] : own(name, 4 * OS.total);
  ob.OD = own('OD', 4 * ODTOTAL); ob.OF = own('OF', 4 * OF.total);
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
  const setParams = (values) => { params.fill(0); params.set(values); core.writeParams(params); };
  const { compute } = core;

  function tendency(IN, OUT) {
    const g = group(IN, OUT);
    setParams([0, 0]);
    compute((pass) => {
      dispatch(pass, 'oVertexVort', g, V);
      dispatch(pass, 'oEdgeFluxPV', g, E);
      dispatch(pass, 'oCellTendency', g, C);
      dispatch(pass, 'oDivCurl', g, C + V);
      dispatch(pass, 'oLapVelocity', g, E);
    });
    setParams([0, 1, relaxRate, stepDt]);
    compute((pass) => {
      dispatch(pass, 'oDivCurl', g, C + V);
      dispatch(pass, 'oLapVelocity', g, E);
      dispatch(pass, 'oKineticPhi', g, C);
      dispatch(pass, 'oMomentum', g, E);
    });
  }
  function advanceState(next, stage, factor) {
    setParams([factor]);
    compute((pass) => dispatch(pass, 'oAdvance', group(ob.S, next, stage), Math.ceil(OS.total / 4)));
  }

  let relaxRate = 1 / 3600, stepDt = 3600;
  function barotropicStep(dt, k1) {
    const g1 = group(ob.S, k1);
    compute((pass) => {
      dispatch(pass, 'oBarotropicSetup', g1, Math.max(C, E));
      dispatch(pass, 'oBarotropicCoriolis', g1, E);
    });
    core.clearBuffer(ob.OD, 4 * B.BAVG, 4 * (C + E));
    const M = Math.max(1, Math.ceil(dt / substepLimit)), dtb = dt / M;
    const gAny = group(ob.S, ob.K1);
    // P0=dtb/2 (stage1,2 trial factor), P1=dtb (stage3 trial factor), P2=dtb/6 (final combine), P3=1/M (average accumulation).
    setParams([dtb / 2, dtb, dtb / 6, 1 / M]);
    compute((pass) => {
      for (let m = 0; m < M; m++) for (let stage = 1; stage <= 4; stage++) dispatch(pass, `oBarotropicStage${stage}`, gAny, Math.max(C, E));
    });
  }

  function dynamics(dt) {
    relaxRate = Math.min(1 / 3600, 1 / dt); stepDt = dt;
    compute((pass) => dispatch(pass, 'oGradEta', group(ob.S, ob.T), E));
    tendency(ob.S, ob.K1);
    barotropicStep(dt, ob.K1);
    advanceState(ob.T, ob.K1, dt / 2);
    tendency(ob.T, ob.K2); advanceState(ob.T, ob.K2, dt / 2);
    tendency(ob.T, ob.K3); advanceState(ob.T, ob.K3, dt);
    tendency(ob.T, ob.K4);
    setParams([dt / 6]);
    compute((pass) => dispatch(pass, 'oCombine', group(ob.S, ob.K1, ob.K2, ob.K3, ob.K4), Math.ceil(OS.total / 2)));
  }
  const rescalePass = (pass, g) => dispatch(pass, 'oRescale', g, C);
  const shiftPass = (pass, g) => { dispatch(pass, 'oEdgeThickness', g, E); dispatch(pass, 'oVelocityShiftClamp', g, E); };
  let edgesCurrent = false;
  function step(dt) {
    const parts = Math.max(1, Math.ceil(dt / inertialLimit)), part = dt / parts;
    for (let m = 0; m < parts; m++) {
      dynamics(part);
      const g = group(ob.S, ob.T);
      if (eddyLimit > 0) setParams([part, eddyLimit / part]);
      compute((pass) => {
        rescalePass(pass, g);
        if (eddyLimit > 0) eddyPasses(pass, g);
        shiftPass(pass, g);
      });
    }
    edgesCurrent = true;
  }
  /*
   * step's parts in its order, for the tests that follow heat and salt
   * through it: dynamics(dt), rescale(), eddy(dt), shift(), then the
   * regrid and the column physics.
   */
  const phases = {
    dynamics(dt) { dynamics(dt); edgesCurrent = false; },
    rescale() { compute((pass) => rescalePass(pass, group(ob.S, ob.T))); edgesCurrent = false; },
    eddy: eddyTransport,
    shift() { compute((pass) => shiftPass(pass, group(ob.S, ob.T))); edgesCurrent = true; },
    regrid, column: columnPhysics,
  };
  /*
   * The ALE step and column physics after the shift: the targets
   * of every column, the edges' velocities and the cells' heat and salt
   * remapped onto them (P[0] the relax), then each column's mixing and
   * each edge's implicit solve (P[0] the step).
   */
  function regrid(dt) {
    setParams([o.regridTime > 0 && Number.isFinite(dt) ? Math.min(1, dt / o.regridTime) : 1]);
    const g = group(ob.S, ob.T);
    compute((pass) => {
      dispatch(pass, 'oRegridTargets', g, C);
      dispatch(pass, 'oRegridEdges', g, E);
      dispatch(pass, 'oRegridCells', g, C);
    });
    edgesCurrent = false;
  }
  function columnPhysics(dt) {
    setParams([dt]);
    const g = group(ob.S, ob.T);
    compute((pass) => {
      dispatch(pass, 'oColumnCells', g, C);
      dispatch(pass, 'oColumnEdges', g, E);
    });
  }

  function eddyPasses(pass, g) {
    dispatch(pass, 'oEddyFlux', g, E);
    dispatch(pass, 'oEddyApply', g, L * C);
    dispatch(pass, 'oEddyLeftover', g, C);
  }
  function eddyTransport(dt) {
    if (!(eddyLimit > 0)) return;
    setParams([dt, eddyLimit / dt]);
    compute((pass) => eddyPasses(pass, group(ob.S, ob.T)));
    edgesCurrent = false;
  }
  function readSurface(surfaceT, ice, concentration = seaCover) {
    device.queue.writeBuffer(ob.OD, 4 * OD.SURFT, Float32Array.from(surfaceT));
    device.queue.writeBuffer(ob.OD, 4 * OD.SURFICE, Float32Array.from(ice));
    device.queue.writeBuffer(ob.OD, 4 * OD.SURFCONC, concentration ? Float32Array.from(concentration) : new Float32Array(C));
    compute((pass) => dispatch(pass, 'oReadSurface', group(ob.S, ob.T), C));
  }
  function readSurfaceFromAtmosphere() {
    core.encode((encoder) => {
      encoder.copyBufferToBuffer(buffers.S, 4 * core.layout.S.TS, ob.OD, 4 * OD.SURFT, 4 * C);
      encoder.copyBufferToBuffer(buffers.S, 4 * core.layout.S.ICE, ob.OD, 4 * OD.SURFICE, 4 * C);
      encoder.copyBufferToBuffer(buffers.PH, 4 * core.layout.PH.CONC, ob.OD, 4 * OD.SURFCONC, 4 * C);
    });
    compute((pass) => dispatch(pass, 'oReadSurface', group(ob.S, ob.T), C));
  }
  function stressFromAtmosphere() {
    compute((pass) => dispatch(pass, 'oStressFromAtmosphere', group(ob.S, ob.T), E));
  }
  const coverOf = (ice, concentration, i) => (ice[i] > 0 ? (concentration && concentration[i] > 0 ? concentration[i] : 1) : 0);
  function setStress(total, ice, concentration = null) {
    const masked = new Float32Array(E), cover = (i) => (shelf[i] ? 0 : coverOf(ice, concentration, i));
    for (let e = 0; e < E; e++) masked[e] = edgeOcean[e] ? total[e] * (1 - 0.5 * (cover(mesh.cellsOnEdge[2 * e]) + cover(mesh.cellsOnEdge[2 * e + 1])) * (1 - o.iceStressTransmission)) * (1 - 0.5 * (shelf[mesh.cellsOnEdge[2 * e]] + shelf[mesh.cellsOnEdge[2 * e + 1]])) : 0;
    device.queue.writeBuffer(ob.OD, 4 * OD.STRESS, masked);
  }
  function salt(dt) {
    setParams([0, 0, 0, 0, 0, 0, dt]);
    compute((pass) => { if (meltOn) dispatch(pass, 'oShelfBudget', group(ob.S, ob.T), WORKGROUP); dispatch(pass, 'oSalt', group(ob.S, ob.T), C); });
  }
  function writeSurface(dt) {
    setParams([0, 0, 0, 0, 0, 0, dt]);
    compute((pass) => dispatch(pass, 'oWriteSurface', group(ob.S, ob.T), C));
  }
  function accumulateFreshwater(dt) {
    setParams([0, 0, 0, 0, 0, 0, dt]);
    compute((pass) => dispatch(pass, 'oAccumulateFresh', group(ob.S, ob.T), C));
  }
  function forgetAccumulated(which) {
    device.queue.writeBuffer(ob.OD, 4 * (which === 'rain' ? OD.RAINSEEN : OD.RUNOFFSEEN), new Float32Array(C));
    if (which !== 'rain') device.queue.writeBuffer(ob.OD, 4 * OD.DISCHARGESEEN, new Float32Array(C));
  }

  /*
   * A single testable/coupled ocean step at double-precision-matched
   * physics but single-precision arithmetic: readSurface, setStress,
   * step, regrid, columnPhysics, salt, writeSurface, unconditionally (the caller
   * decides when to call it; model.gpu.js gates on its own everySteps
   * counter so it only pays for the atmosphere-state read and the step
   * on the steps that need it).
   */
  async function advance(surfaceT, ice, totalStress, dt, concentration = seaCover) {
    if (driftOptions) throw new Error('the ice drift runs in advanceCoupled, where the ice lives on the device');
    readSurface(surfaceT, ice, concentration);
    setStress(typeof totalStress === 'function' ? totalStress() : totalStress, ice, concentration);
    step(dt);
    regrid(dt); columnPhysics(dt);
    salt(dt);
    writeSurface(dt);
    await device.queue.onSubmittedWorkDone();
    return true;
  }
  let driftCalls = 0, driftSpacing = 0;
  for (let e = 0; e < E; e++) driftSpacing += mesh.dcEdge[e] / E;
  const driftMask = driftOptions ? Float32Array.from({ length: E }, (_, e) => {
    const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1];
    const side = driftOptions.hemisphere === 'both' || (driftOptions.hemisphere === 'south' ? mesh.latEdge[e] < 0 : mesh.latEdge[e] > 0);
    return edgeOcean[e] && cellOpen[a] && cellOpen[b] && side ? 1 : 0;
  }) : null;
  function driftStep(dt) {
    const d = driftOptions, g = group(ob.S, ob.T);
    if (++driftCalls % d.dynamicsEvery !== 0) { compute((pass) => dispatch(pass, 'iClearArrived', g, C)); return; }
    const interval = dt * d.dynamicsEvery, { alpha, beta, passes } = iceDriftPasses(d, interval, o.iceDensity, driftSpacing);
    const sub = Math.max(1, Math.ceil(d.speedLimit * interval / (minSpacing * d.courantLimit)));
    const moving = d.rheology !== 'frozen', cavitating = d.rheology === 'cavitating';
    setParams([interval, alpha, beta, interval / sub]);
    compute((pass) => {
      dispatch(pass, 'iPrepareCells', g, C);
      dispatch(pass, 'iPrepareEdges', g, E);
      if (moving) {
        for (let p = 0; p < passes; p++) {
          if (cavitating) dispatch(pass, p % 2 ? 'iStressB' : 'iStressA', g, C);
          dispatch(pass, p % 2 ? 'iMomentumBA' : 'iMomentumAB', g, E);
        }
        if (passes % 2) dispatch(pass, 'iCopy', g, E);
        dispatch(pass, 'iStressOut', g, E);
      }
      for (let s = 0; s < sub; s++) { dispatch(pass, 'iFlux', g, E); dispatch(pass, s === sub - 1 ? 'iUpdateLast' : 'iUpdateMid', g, C); }
    });
  }
  function setAirStress(values) {
    if (!given) throw new Error('setAirStress needs iceDrift: { airStress: \'given\' }');
    device.queue.writeBuffer(ob.OD, 4 * OD.TAUA, Float32Array.from(values));
  }
  async function advanceCoupled(dt) {
    if (driftOptions) driftStep(dt);
    readSurfaceFromAtmosphere();
    stressFromAtmosphere();
    step(dt);
    regrid(dt); columnPhysics(dt);
    salt(dt);
    writeSurface(dt);
  }

  /*
   * Initialization, loading, serialization and diagnostics run once per
   * model build or at the (infrequent) cadence a caller asks for, so
   * they port the CPU functions over plain double-precision arrays and
   * push the result to the GPU, rather than duplicating that logic in
   * WGSL.
   */
  const eos = seawaterDensity;
  /*
   * The free surface that levels the pressure at `referenceDepth` in
   * every column of the open ocean's abyss, so the deep ocean starts
   * without a barotropic pressure gradient; every other column takes the
   * free surface of the water around it, found by relaxation from the
   * abyss's mean. Ported over plain double-precision arrays exactly as
   * the CPU module's initialize()/stericSurface() do, since this runs
   * once per model build rather than every ocean step.
   */
  function stericSurface(h, Q, W, eta, referenceDepth = 3500) {
    const at = (k, i) => k * C + i;
    const densityAt = (k, i) => { const n = at(k, i); return h[n] > 1e-6 ? eos(Q[n] / h[n], W[n] / h[n]) : rho[k]; };
    const deep = abyssalCells(mesh, D, cellOcean, referenceDepth);
    let deepArea = 0, deepMean = 0;
    for (let i = 0; i < C; i++) {
      eta[i] = 0;
      if (!deep[i]) continue;
      let budget = referenceDepth, anomaly = 0;
      for (let k = 0; k < L && budget > 0; k++) {
        const part = Math.min(h[at(k, i)], budget);
        anomaly += (densityAt(k, i) - o.density) * part;
        budget -= part;
      }
      eta[i] = -anomaly / o.density;
      deepArea += mesh.areaCell[i]; deepMean += mesh.areaCell[i] * eta[i];
    }
    if (deepArea > 0) for (let i = 0; i < C; i++) if (deep[i]) eta[i] -= deepMean / deepArea;
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
  /*
   * The start, build of ocean/layered.module.js: every
   * sea column on the depth table clipped to its floor, the atlas's layer
   * means where the atlas covers the cell, the nearest covered column's
   * where it does not, under a shelf the nearest open column's at the same
   * depth (its bottom water below its floor, no colder than the in-situ
   * freezing point), the analytic profile without an atlas; the steric
   * free surface; then one regrid on the device once uploaded.
   */
  function hybridGrid(i, eta, zCol) {
    const H = D[i] + eta[i];
    for (let k = 1; k < L; k++) zCol[k] = Math.min(table.depths[k], H);
    spaceInterfaces(zCol, L, H, o.minimumLayer);
    return H;
  }
  function initializeArrays(surfaceT, ice, atlas = null) {
    const h = new Float64Array(L * C), u = new Float64Array(L * E), Q = new Float64Array(L * C), W = new Float64Array(L * C);
    const eta = new Float64Array(C), T0 = new Float64Array(C), previousT0 = new Float64Array(C);
    const at = (k, i) => k * C + i, zCol = new Float64Array(L + 1), filled = new Uint8Array(C), zs = [], ts = [], ss = [], pending = [];
    const integral = (v, a, b) => {
      let sum = 0;
      for (let j = 1; j < zs.length; j++) {
        const lo = Math.max(a, zs[j - 1]), hi = Math.min(b, zs[j]);
        if (hi <= lo) continue;
        sum += (v[j - 1] + (v[j] - v[j - 1]) / (zs[j] - zs[j - 1]) * (0.5 * (lo + hi) - zs[j - 1])) * (hi - lo);
      }
      return sum;
    };
    for (let i = 0; i < C; i++) {
      const iced = ice[i] > 0 || shelf[i] === 1, lat = mesh.latCell[i];
      T0[i] = iced ? FREEZING_POINT : surfaceT[i];
      if (!cellOcean[i]) continue;
      hybridGrid(i, eta, zCol);
      const lid = draft[i], bed = lid + D[i];
      let column = atlas && !shelf[i] ? atlas.columnAt(lat, mesh.lonCell[i]) : null;
      if (column) {
        const reached = column.depths.length;
        if (atlas.columnReaching && atlas.depths && reached < atlas.depths.length && atlas.depths[reached] <= bed) {
          const deep = atlas.columnReaching(lat, mesh.lonCell[i], bed);
          if (deep && deep.depths.length > reached) {
            const Td = Float64Array.from(deep.T), Sd = Float64Array.from(deep.S);
            Td.set(column.T); Sd.set(column.S);
            column = { depths: deep.depths, T: Td, S: Sd };
          }
        }
        filled[i] = 1;
        zs.length = ts.length = ss.length = 0;
        for (let j = 0; j < column.depths.length && column.depths[j] < bed; j++) { zs.push(column.depths[j]); ts.push(column.T[j]); ss.push(column.S[j]); }
        const [tb, sb] = profileAt(column, bed);
        zs.push(bed); ts.push(tb); ss.push(sb);
        for (let k = 0; k < L; k++) { const n = at(k, i), hk = zCol[k + 1] - zCol[k]; h[n] = hk; Q[n] = integral(ts, lid + zCol[k], lid + zCol[k + 1]); W[n] = integral(ss, lid + zCol[k], lid + zCol[k + 1]); }
        if (lid > 0) shelfLoad[i] = displacedWeight(lid, (z) => profileAt(column, z));
      } else if (atlas) {
        pending.push(i);
        for (let k = 0; k < L; k++) h[at(k, i)] = zCol[k + 1] - zCol[k];
      } else {
        const s0 = salinityProfile(lat), t0 = T0[i], scale = 400 * (1 - o.thermoclineTilt + 2 * o.thermoclineTilt * Math.cos(lat) ** 2);
        const tDeep = 275.15, sDeep = 34.7, decay = (a, b) => scale * (Math.exp(-a / scale) - Math.exp(-b / scale)) / (b - a);
        if (lid > 0) shelfLoad[i] = displacedWeight(lid, (z) => [tDeep + (t0 - tDeep) * Math.exp(-z / scale), sDeep + (s0 - sDeep) * Math.exp(-z / scale)]);
        for (let k = 0; k < L; k++) {
          const n = at(k, i), hk = zCol[k + 1] - zCol[k], f = decay(lid + zCol[k], lid + zCol[k + 1]);
          h[n] = hk; Q[n] = hk * (tDeep + (t0 - tDeep) * f); W[n] = hk * (sDeep + (s0 - sDeep) * f);
        }
      }
    }
    for (const i of pending) {
      let source = -1, ring = [i];
      const seen = new Set(ring);
      for (let r = 0; r < 40 && source < 0 && ring.length; r++) {
        const next = [];
        for (const j of ring) for (let m = 0; m < mesh.nEdgesOnCell[j]; m++) {
          const c = mesh.cellsOnCell[mesh.maxEdges * j + m];
          if (seen.has(c)) continue;
          seen.add(c);
          if (filled[c]) { source = c; break; }
          if (cellOcean[c]) next.push(c);
        }
        ring = next;
      }
      const sourceAt = (y) => {
        let top = draft[source], j = 0, last = at(0, source);
        while (j < L - 1 && top + h[at(j, source)] < y) { if (h[at(j, source)] > 1e-9) last = at(j, source); top += h[at(j, source)]; j++; }
        return shelf[i] && !(h[at(j, source)] > 1e-9) ? last : at(j, source);
      };
      if (draft[i] > 0) shelfLoad[i] = displacedWeight(draft[i], (y) => { const m = source >= 0 ? sourceAt(y) : -1; return m >= 0 && h[m] > 1e-9 ? [Q[m] / h[m], W[m] / h[m]] : [FREEZING_POINT, o.referenceS]; });
      let z = draft[i];
      for (let k = 0; k < L; k++) {
        const n = at(k, i), hk = h[n], center = z + 0.5 * hk;
        z += hk;
        let t = T0[i], sal = o.referenceS;
        if (source >= 0) {
          const m = sourceAt(center);
          if (h[m] > 1e-9) { t = Q[m] / h[m]; sal = W[m] / h[m]; }
        }
        if (shelf[i]) t = Math.max(t, SHELF_MELT.liquidusSlope * sal + baseFreezing(o.gravity * (shelfLoad[i] + o.density * (center - draft[i]))));
        Q[n] = hk * t; W[n] = hk * sal;
      }
    }
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) continue;
      T0[i] = shelf[i] ? Q[i] / h[i] : ice[i] > 0 ? FREEZING_POINT : Math.max(FREEZING_POINT, Q[i] / h[i]);
      if (!shelf[i]) Q[i] = h[i] * T0[i];
      previousT0[i] = T0[i];
    }
    stericSurface(h, Q, W, eta);
    return { h, u, Q, W, eta, T0, previousT0, filled: atlas ? filled : null };
  }
  function uploadArrays({ h, u, Q, Qref, W, Wref, eta }, surfaceT, ice, restart = null, concentration = seaCover, regridOnDevice = !restart) {
    edgesCurrent = false;
    const packed = new Float32Array(OS.total);
    packed.set(h, OS.OH); packed.set(u, OS.OU);
    packed.set(Qref ?? Float64Array.from(Q, (q, n) => q - HEAT_REFERENCE * h[n]), OS.OQ);
    packed.set(Wref ?? Float64Array.from(W, (w, n) => w - SALT_REFERENCE * h[n]), OS.OW);
    for (let e = 0; e < E; e++) if (!edgeOcean[e]) packed[OS.OU + e] = 0;
    writeInPieces(device, ob.S, 0, packed);
    const encoder = device.createCommandEncoder();
    encoder.copyBufferToBuffer(ob.S, 0, ob.T, 0, 4 * OS.total);
    for (const b of [ob.K1, ob.K2, ob.K3, ob.K4]) encoder.clearBuffer(b, 0, 4 * OS.total);
    encoder.clearBuffer(ob.OD, 0, 4 * ODTOTAL);
    device.queue.submit([encoder.finish()]);
    device.queue.writeBuffer(ob.OD, 4 * OD.ETA, Float32Array.from(eta));
    device.queue.writeBuffer(ob.OD, 4 * OD.EMASK, Float32Array.from(edgeOcean));
    device.queue.writeBuffer(ob.OD, 4 * OD.CMASK, Float32Array.from(cellOcean));
    device.queue.writeBuffer(ob.OD, 4 * OD.SHELF, Float32Array.from(shelf));
    device.queue.writeBuffer(ob.OD, 4 * OD.DRAFT, Float32Array.from(draft));
    const load = new Float32Array(C), loadT = new Float32Array(C), loadS = new Float32Array(C);
    for (let i = 0; i < C; i++) {
      if (!shelf[i]) continue;
      loadT[i] = Math.fround(packed[OS.OQ + i] / packed[OS.OH + i]); loadS[i] = Math.fround(packed[OS.OW + i] / packed[OS.OH + i]);
      load[i] = shelfLoad[i] - seawaterDensity(loadT[i] + HEAT_REFERENCE, loadS[i] + SALT_REFERENCE) * draft[i];
    }
    device.queue.writeBuffer(ob.OD, 4 * OD.LOAD, load); device.queue.writeBuffer(ob.OD, 4 * OD.LOADT, loadT); device.queue.writeBuffer(ob.OD, 4 * OD.LOADS, loadS);
    device.queue.writeBuffer(ob.OD, 4 * OD.LIQUIDUS, Float32Array.from(shelfLoad, (weight, i) => (shelf[i] ? baseFreezing(o.gravity * weight, 273.15 - HEAT_REFERENCE) : 0)));
    device.queue.writeBuffer(ob.OD, 4 * OD.MELTSUM, Float32Array.from(shelfMelted)); device.queue.writeBuffer(ob.OD, 4 * OD.DRAW, Float32Array.from(shelfDraw));
    device.queue.writeBuffer(ob.OD, 4 * OD.POOLS, Float32Array.from([0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0]));
    device.queue.writeBuffer(ob.OD, 4 * OD.BATH, Float32Array.from(D));
    device.queue.writeBuffer(ob.OD, 4 * OD.FEDGE, Float32Array.from(mesh.fEdge));
    if (third) { const beyond = beyondCells(mesh); device.queue.writeBuffer(ob.OD, 4 * OD.BEY, Float32Array.from(beyond.cell)); device.queue.writeBuffer(ob.OD, 4 * OD.BED, Float32Array.from(beyond.edge)); }
    device.queue.writeBuffer(ob.OD, 4 * OD.EDDYK, Float32Array.from(eddyKappa));
    device.queue.writeBuffer(ob.OD, 4 * OD.FCELL, Float32Array.from(mesh.fCell));
    device.queue.writeBuffer(ob.OD, 4 * OD.SILL, Float32Array.from(sillEdge));
    if (driftOptions) {
      device.queue.writeBuffer(ob.OD, 4 * OD.IDYN, driftMask);
      if (restart && restart.iceVelocity && restart.iceVelocity.length === E) device.queue.writeBuffer(ob.OD, 4 * OD.UICE, Float32Array.from(restart.iceVelocity));
      if (restart && restart.iceSigma && restart.iceSigma.length === C) device.queue.writeBuffer(ob.OD, 4 * OD.SIGMA, Float32Array.from(restart.iceSigma));
      driftCalls = 0;
    }
    device.queue.writeBuffer(ob.OD, 4 * OD.MIXED, restart && restart.mixedDepths && restart.mixedDepths.length === C ? Float32Array.from(restart.mixedDepths) : mixedDepthsOf(h, Q, W));
    if (restart && restart.boundaryDepths && restart.boundaryDepths.length === C) device.queue.writeBuffer(ob.OD, 4 * OD.BLD, Float32Array.from(restart.boundaryDepths));
    device.queue.writeBuffer(ob.OD, 4 * OD.DRAINSTART, drainStart);
    if (drainCells.length) { device.queue.writeBuffer(ob.OD, 4 * OD.DRAINCELL, Float32Array.from(drainCells)); device.queue.writeBuffer(ob.OD, 4 * OD.DRAINW, Float32Array.from(drainWeights)); }
    device.queue.writeBuffer(ob.OD, 4 * OD.PREVVOL, restart ? Float32Array.from(restart.previousIceVolume) : Float32Array.from(ice, (_, i) => coverOf(ice, concentration, i) * ice[i]));
    if (restart) device.queue.writeBuffer(ob.OD, 4 * OD.FRESH, Float32Array.from(restart.fresh));
    if (restart && restart.calved && restart.calved.length === C) device.queue.writeBuffer(ob.OD, 4 * OD.CALVED, Float32Array.from(restart.calved));
    if (restart) device.queue.writeBuffer(ob.OD, 4 * OD.BUOY, Float32Array.from(restart.buoyancy));
    const capacity = restart ? restart.capacity : Float64Array.from({ length: C }, (_, i) => o.density * o.specificHeat * Math.max(h[i], 1));
    const kept = (name) => (restart && restart[name] && restart[name].length === C ? restart[name] : null);
    core.uploadPhysics({ capacity, oceanFlux: restart ? restart.flux : null, transmission: o.penetration ? kept('transmission') ?? transmissionOf(h) : null, penetrated: kept('penetrated') });
    if (regridOnDevice) regrid(Infinity);
  }
  function mixedDepthsOf(h, Q, W) {
    const out = new Float32Array(C), hCol = new Float64Array(L), rhoCol = new Float64Array(L);
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) continue;
      for (let k = 0; k < L; k++) { const n = k * C + i, hk = Math.max(EPS, h[n]); hCol[k] = h[n]; rhoCol[k] = eos(Q[n] / hk, W[n] / hk); }
      out[i] = mixedLayerDepth(hCol, rhoCol, L);
    }
    return out;
  }
  function transmissionOf(h) {
    const out = new Float64Array(C);
    for (let i = 0; i < C; i++) if (cellOpen[i]) out[i] = transmittedShare(h[i], o.penetration);
    return out;
  }
  function initialize(surfaceT, ice, { climatology: atlas = o.climatology ?? null, concentration = seaCover } = {}) {
    const arrays = initializeArrays(surfaceT, ice, atlas);
    uploadArrays(arrays, surfaceT, ice, null, concentration);
    device.queue.writeBuffer(ob.OD, 4 * OD.PREVT0, Float32Array.from(arrays.previousT0));
    const { filled, T0 } = arrays;
    if (!filled) return null;
    let covered = 0, other = 0;
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) continue;
      if (filled[i]) covered++; else other++;
      if (filled[i] && !(ice[i] > 0) && !shelf[i]) surfaceT[i] = T0[i];
    }
    return { atlas: covered, analytic: other };
  }
  const RESTART = [['previousT0', C], ['previousIceVolume', C], ['fresh', C], ['buoyancy', C], ['flux', C], ['capacity', C]];
  function restartable(saved) {
    if (!saved || !saved.h || saved.h.length !== L * C || ![saved.Qref ?? saved.Q, saved.Wref ?? saved.W].every((store) => store && store.length === L * C) || !RESTART.every(([name, length]) => saved[name] && saved[name].length === length)) return false;
    for (let i = 0; i < C; i++) {
      let sum = 0;
      for (let k = 0; k < L; k++) sum += saved.h[k * C + i];
      if (cellOcean[i] ? !(Math.abs(sum - D[i] - saved.eta[i]) < 1) : sum !== 0) return false;
    }
    return true;
  }
  /*
   * A saved ocean on these labels is uploaded as it is, its columns
   * stretched to this bathymetry's depth; a saved slab ocean is converted as
   * convertSlab of ocean/layered.module.js does, each column remapped by the
   * parabolas onto the depth table at its floor and each edge's velocity
   * piecewise-linearly onto the mean of its cells' new thicknesses, and
   * regridded on the device; any other state starts the ocean afresh. A
   * state that carries the restart arrays (restartArrays) on these labels
   * goes up as it stands, stores and all.
   */
  function upload(saved, surfaceT, ice, concentration = seaCover) {
    const surfaceAt = (i) => (ice[i] > 0 ? FREEZING_POINT : surfaceT[i]);
    const at = (k, i) => k * C + i;
    const start = () => { const arrays = initializeArrays(surfaceT, ice, o.climatology ?? null); uploadArrays(arrays, surfaceT, ice, null, concentration); device.queue.writeBuffer(ob.OD, 4 * OD.PREVT0, Float32Array.from(arrays.previousT0)); };
    if (!saved || !saved.h) { start(); return; }
    if (hasShelves) for (let i = 0; i < C; i++) shelfLoad[i] = !shelf[i] ? 0 : saved.shelfLoad && saved.shelfLoad.length === C ? saved.shelfLoad[i] : draft[i] * seawaterDensity(saved.T[i], saved.S[i]);
    for (const [name, field] of [['shelfMelted', shelfMelted], ['shelfDraw', shelfDraw]]) field.set(saved[name] && saved[name].length === C ? saved[name] : new Float64Array(C));
    const classes = savedDensities(saved, C), own = classes && sameDensities(classes, labels) && saved.h.length === L * C;
    if (own && restartable(saved)) {
      uploadArrays(saved, surfaceT, ice, saved, concentration, false);
      device.queue.writeBuffer(ob.OD, 4 * OD.PREVT0, Float32Array.from(saved.previousT0));
      return;
    }
    const h = new Float64Array(L * C), u = new Float64Array(L * E), Q = new Float64Array(L * C), W = new Float64Array(L * C), eta = Float64Array.from(saved.eta);
    const zCol = new Float64Array(L + 1);
    if (own) {
      h.set(saved.h); u.set(saved.u);
      const Qs = saved.Qref ? Float64Array.from(saved.Qref, (q, n) => q + HEAT_REFERENCE * h[n]) : null, Ws = saved.Wref ? Float64Array.from(saved.Wref, (w, n) => w + SALT_REFERENCE * h[n]) : null;
      for (let n = 0; n < L * C; n++) { Q[n] = Qs ? Qs[n] : h[n] * saved.T[n]; W[n] = Ws ? Ws[n] : h[n] * saved.S[n]; }
    } else {
      const Ls = classes ? classes.length + 1 : Math.round(saved.h.length / C);
      if (saved.h.length !== Ls * C) throw new Error(`a saved ocean on ${Ls} layers holds ${Ls * C} thicknesses, not ${saved.h.length}`);
      const scratch = columnScratch(Math.max(Ls, L)), hs = new Float64Array(Ls), qs = new Float64Array(Ls), hNew = new Float64Array(L), amount = new Float64Array(L);
      for (let i = 0; i < C; i++) {
        if (!cellOcean[i]) continue;
        const H = hybridGrid(i, eta, zCol);
        let sum = 0;
        for (let j = 0; j < Ls; j++) { hs[j] = Math.max(0, saved.h[j * C + i]); sum += hs[j]; }
        if (!(sum > 0)) { for (let k = 0; k < L; k++) { const n = at(k, i), hk = zCol[k + 1] - zCol[k]; h[n] = hk; Q[n] = hk * surfaceAt(i); W[n] = hk * o.referenceS; } continue; }
        for (let j = 0; j < Ls; j++) hs[j] *= H / sum;
        for (let k = 0; k < L; k++) hNew[k] = zCol[k + 1] - zCol[k];
        for (let j = 0; j < Ls; j++) qs[j] = saved.T[j * C + i];
        remapColumn(hs, qs, Ls, hNew, amount, scratch, { layers: L });
        for (let k = 0; k < L; k++) Q[at(k, i)] = amount[k];
        for (let j = 0; j < Ls; j++) qs[j] = saved.S[j * C + i];
        remapColumn(hs, qs, Ls, hNew, amount, scratch, { layers: L });
        for (let k = 0; k < L; k++) { W[at(k, i)] = amount[k]; h[at(k, i)] = hNew[k]; }
      }
      const us = new Float64Array(Ls), hEdgeNew = new Float64Array(L);
      for (let e = 0; e < E; e++) {
        if (!edgeOcean[e]) continue;
        const a = mesh.cellsOnEdge[2 * e], b = mesh.cellsOnEdge[2 * e + 1];
        for (let j = 0; j < Ls; j++) { hs[j] = 0.5 * (Math.max(0, saved.h[j * C + a]) + Math.max(0, saved.h[j * C + b])); us[j] = saved.u[j * E + e]; }
        for (let k = 0; k < L; k++) hEdgeNew[k] = 0.5 * (h[at(k, a)] + h[at(k, b)]);
        remapColumn(hs, us, Ls, hEdgeNew, amount, scratch, { order: 2, layers: L });
        for (let k = 0; k < L; k++) u[k * E + e] = hEdgeNew[k] > 0 ? amount[k] / hEdgeNew[k] : 0;
      }
    }
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) { for (let k = 0; k < L; k++) { h[at(k, i)] = 0; Q[at(k, i)] = 0; W[at(k, i)] = 0; } eta[i] = 0; continue; }
      let sum = 0;
      for (let k = 0; k < L; k++) sum += h[at(k, i)];
      const H = D[i] + eta[i];
      if (!(sum > 0)) { hybridGrid(i, eta, zCol); for (let k = 0; k < L; k++) { const n = at(k, i), hk = zCol[k + 1] - zCol[k]; h[n] = hk; Q[n] = hk * surfaceAt(i); W[n] = hk * o.referenceS; } continue; }
      if (Math.abs(sum - H) > 1e-9) { const scale = H / sum; for (let k = 0; k < L; k++) { const n = at(k, i); h[n] *= scale; Q[n] *= scale; W[n] *= scale; } }
    }
    for (let e = 0; e < E; e++) if (!edgeOcean[e]) for (let k = 0; k < L; k++) u[k * E + e] = 0;
    const restart = own && restartable({ ...saved, h, eta }) ? saved : null;
    uploadArrays({ h, u, Q, W, eta }, surfaceT, ice, restart, concentration, !own);
    const previousT0 = new Float64Array(C);
    for (let i = 0; i < C; i++) previousT0[i] = cellOcean[i] ? Q[i] / h[i] : surfaceAt(i);
    device.queue.writeBuffer(ob.OD, 4 * OD.PREVT0, Float32Array.from(restart ? restart.previousT0 : previousT0));
  }
  async function readShelfMelt() {
    const [melted, shelfMelted, shelfDraw, pools] = await readRanges(device, ob.OD, [{ offset: OD.MELT, length: C }, { offset: OD.MELTSUM, length: C }, { offset: OD.DRAW, length: C }, { offset: OD.POOLS, length: 12 }]);
    return { melted, shelfMelted, shelfDraw, pools };
  }
  async function download() {
    const stateRead = readRanges(device, ob.S, [{ offset: 0, length: OS.total }]), etaRead = readRanges(device, ob.OD, [{ offset: OD.ETA, length: C }]);
    const [s] = await stateRead;
    const [etaRange] = await etaRead;
    const h = Float64Array.from(s.subarray(OS.OH, OS.OH + L * C));
    const u = Float64Array.from(s.subarray(OS.OU, OS.OU + L * E));
    const Qraw = s.subarray(OS.OQ, OS.OQ + L * C), Wraw = s.subarray(OS.OW, OS.OW + L * C);
    const T = new Float64Array(L * C), S = new Float64Array(L * C);
    for (let n = 0; n < L * C; n++) { const hh = Math.max(EPS, h[n]); T[n] = (Qraw[n] + HEAT_REFERENCE * h[n]) / hh; S[n] = (Wraw[n] + SALT_REFERENCE * h[n]) / hh; }
    const eta = Float64Array.from(etaRange);
    const shelfFields = hasShelves ? await readRanges(device, ob.OD, [{ offset: OD.MELTSUM, length: C }, { offset: OD.DRAW, length: C }]) : null;
    const h1 = Float64Array.from((await readRanges(device, ob.OD, [{ offset: OD.MIXED, length: C }]))[0]), T1 = T.subarray(0, C), S1 = S.subarray(0, C), u1 = u.subarray(0, E);
    const thermoclineDepth = new Float64Array(C);
    for (let i = 0; i < C; i++) thermoclineDepth[i] = cellOcean[i] ? thermoclineOf(h, T, S, i) : NaN;
    return { h, u, T, S, eta, h1, T1, S1, u1, thermoclineDepth, layers: L, ...(shelfFields ? { shelfMelted: Float64Array.from(shelfFields[0]), shelfDraw: Float64Array.from(shelfFields[1]) } : {}) };
  }
  function thermoclineOf(h, T, S, i) {
    let z = 0, previous = -Infinity, center = 0;
    for (let k = 0; k < L; k++) {
      const n = k * C + i, hv = h[n], rk = Math.max(previous, hv > 1e-6 ? eos(T[n], S[n]) : rho[k]), c = z + 0.5 * hv;
      if (rk >= THERMOCLINE_DENSITY) return k === 0 ? 0 : center + (THERMOCLINE_DENSITY - previous) / Math.max(rk - previous, 1e-12) * (c - center);
      previous = rk; center = c; z += hv;
    }
    return z;
  }
  function serializeFrom(d) {
    return { h: Array.from(d.h), u: Array.from(d.u), T: Array.from(d.T), S: Array.from(d.S), eta: Array.from(d.eta), densities: Array.from(labels), vertical: 'hybrid', generator: o.generator, ...(hasShelves ? { shelfLoad: Array.from(shelfLoad), ...(d.shelfMelted ? { shelfMelted: Array.from(d.shelfMelted), shelfDraw: Array.from(d.shelfDraw) } : {}) } : {}) };
  }
  async function serialize() { return serializeFrom(await download()); }
  /*
   * What upload needs beside serialize()'s fields to put the device back
   * exactly as it stands between two ocean steps: the heat and salt
   * contents as stored, Qref and Wref (upload also takes a Q of h·T and a
   * W of h·S), the surface's last temperature and the ice volume it
   * last saw, the freshwater not yet taken, the land's ice not yet
   * melted, `calved`, the heat flux and capacity
   * the surface update is using, the sunlight not yet deposited,
   * `penetrated`, the mixed-layer and boundary-layer depths the column
   * physics left, and with penetration the share of the sunlight that
   * passes the top layer,
   * `transmission`. A saved ocean that carries all of them at this mesh's
   * sizes, the last two and `calved` optional, its columns within a meter of this
   * bathymetry and dry on land, is uploaded as it is, without fitting its
   * columns.
   */
  async function restartArrays() {
    const PHL = core.layout.PH;
    const [[Qref, Wref], [previousT0, previousIceVolume, fresh, buoyancy, mixedDepths, boundaryDepths, calved], [flux, capacity, penetrated]] = await Promise.all([
      readRanges(device, ob.S, [{ offset: OS.OQ, length: L * C }, { offset: OS.OW, length: L * C }]),
      readRanges(device, ob.OD, [{ offset: OD.PREVT0, length: C }, { offset: OD.PREVVOL, length: C }, { offset: OD.FRESH, length: C }, { offset: OD.BUOY, length: C }, { offset: OD.MIXED, length: C }, { offset: OD.BLD, length: C }, { offset: OD.CALVED, length: C }]),
      readRanges(device, buffers.PH, [{ offset: PHL.OFLUX, length: C }, { offset: PHL.CAP, length: C }, { offset: PHL.PEN, length: C }]),
    ]);
    const transmission = o.penetration ? (await readRanges(device, buffers.PH, [{ offset: PHL.TRANS, length: C }]))[0] : null;
    const [iceVelocity, iceSigma] = driftOptions ? await readRanges(device, ob.OD, [{ offset: OD.UICE, length: E }, { offset: OD.SIGMA, length: C }]) : [null, null];
    return Object.fromEntries(Object.entries({ Qref, Wref, previousT0, previousIceVolume, fresh, calved, buoyancy, flux, capacity, penetrated, mixedDepths, boundaryDepths, ...(transmission ? { transmission } : {}), ...(iceVelocity ? { iceVelocity, iceSigma } : {}) }).map(([name, values]) => [name, Float32Array.from(values)]));
  }

  function diagnosticsFrom(sums) {
    const { area, depth, heat, interiorT, interiorH, speed, thermocline, salinity, ssh, transport, limited, dried } = sums;
    return { oceanUpperDepth: depth / area, oceanHeat: heat / area, oceanInteriorT: interiorH > 0 ? interiorT / interiorH : 0, oceanSpeed: speed, oceanThermoclineDepth: thermocline / area, oceanSalinity: salinity / area, oceanSSH: ssh, oceanTransport: transport / 1e6, oceanLimited: limited, oceanDried: dried };
  }

  /*
   * Queues the fields the page asked for, those that vary with depth at
   * `depth` metres (0 for the top layer), the vertical velocity at
   * `upwellingDepth` (`depth` by default), and, when asked, the
   * diagnostics reduction, and reads back only those; fields over land
   * or below the sea floor are NaN and current vectors zero there.
   */
  const FIELDS = { sst: ['SST', 1], sss: ['SSS', 1], layerDepth: ['H1', 1], thermocline: ['THD', 1], ssh: ['ETA', 1], currents: ['CUR', 3], current: ['CSPD', 1], upwelling: ['UPW', 1] };
  function frame({ fields = [], diagnostics: summarize = false, depth = 0, upwellingDepth = depth } = {}) {
    const wanted = fields.filter((name) => name in FIELDS);
    if (!wanted.length && !summarize) return Promise.resolve({ fields: {}, diagnostics: null });
    const g = group(ob.S, ob.OF);
    setParams([0, 0, 0, 0, depth, upwellingDepth]);
    const encoder = device.createCommandEncoder(), pass = encoder.beginComputePass();
    if (wanted.length) dispatch(pass, 'oFrame', g, C);
    if (summarize) dispatch(pass, 'oReduce', g, C);
    pass.end();
    device.queue.submit([encoder.finish()]);
    const ranges = wanted.map((name) => ({ name, offset: OF[FIELDS[name][0]], length: FIELDS[name][1] * C }));
    if (summarize) ranges.push({ name: null, offset: OF.PART, length: reduced.length * reductionGroups(C) });
    return readRanges(device, ob.OF, ranges).then((views) => {
      const out = { fields: {}, diagnostics: null };
      views.forEach((view, n) => {
        const name = ranges[n].name;
        if (!name) { out.diagnostics = diagnosticsFrom(finishReduction(reduced, view, C)); return; }
        const width = FIELDS[name][1];
        for (let i = 0; i < C; i++) if (!cellOcean[i] || (width === 1 && view[i] <= -1e30)) for (let c = 0; c < width; c++) view[width * i + c] = width === 1 ? NaN : 0;
        out.fields[name] = view;
      });
      return out;
    });
  }
  async function diagnostics() { return (await frame({ diagnostics: true })).diagnostics; }

  return {
    layers: L, everySteps: o.everySteps, options: o, D, sillEdge, cellOcean, vertical: 'hybrid', generator: o.generator, table,
    useSeaIce(concentration) { seaCover = concentration; },
    initialize, upload, download, serialize, serializeFrom, restartArrays, diagnostics, frame,
    advance, advanceCoupled, accumulateFreshwater, forgetAccumulated, readShelfMelt, driftStep, setAirStress, iceDrift: driftOptions,
    setStress, readSurface, readSurfaceFromAtmosphere, stressFromAtmosphere, salt, writeSurface, step, phases, eddyTransport,
    buffers: ob, layout: { OS, OD, B },
  };
}
