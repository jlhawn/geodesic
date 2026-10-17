import { emptyBuffer, readRanges, reductionKernel, finishReduction, reductionGroups, writeInPieces } from './device.module.js';
import { LAYER_DENSITIES, LAYER_SALINITIES, LAYER_BOTTOMS, THERMOCLINE_DENSITY, EPS, THIN, PV_FLOOR, SPEED_LIMIT, DENSITY_TOLERANCE, RESTORE_TOLERANCE, THINNED_MARGIN, CLOSURE_SPACING, CLOSURE_RIDGE, EDDY_BOTTOM_TAPER, EDDY_SLACK, closureCoefficient, eddyDiffusivities, eddyDiffusionLimit, bathymetryFrom, fitColumns, runoffOutlets, interiorWater, atlasColumns, abyssalCells, savedDensities, sameDensities, rebinOcean, JERLOV, PENETRATION_DEPTH, penetrationShares, FRICTION_PER_WIND, CALM_WIND, RETREAT_EKMAN, mixingOptionCheck, closeDay, DAILY } from '../ocean/layered.module.js';
import { SEAWATER, SEAWATER_WGSL, seawaterDensity, labelTemperature } from '../ocean/seawater.module.js';
import { HYBRID_LAYERS, H_MIN, MODE_DEPTH, GENERATORS, labelBoundaries, columnScratch, spaceInterfaces, remapColumn, mixedLayerDepth, transmitted as transmittedShare } from '../ocean/column.module.js';
import { profileAt } from '../ocean/climatology.module.js';
import { FREEZING_POINT } from '../physics/ice.module.js';

/*
 * The layered ocean of ocean/layered.module.js on the GPU: a bulk mixed
 * layer over interior isopycnal layers, each a TRiSK shallow-water layer
 * carrying thickness, edge velocity, heat h·(T − HEAT_REFERENCE) and salt
 * h·(S − SALT_REFERENCE), each measured from its reference so that a thick
 * layer's float32 store keeps the eddy transport's small changes, and the
 * sunlight it deposits below the mixed layer (`penetration`), whose
 * remainder under each store's spacing waits in PH_PEN. The RK4
 * baroclinic state (h, u, heat, salt for all L layers) lives in a
 * ping-pong buffer set (S/T/K1-K4); everything else the step needs
 * (edge thicknesses, PV, the interior potential, the ∇⁴ closure scratch,
 * the barotropic sub-stepping state, the surface staging) lives in one
 * big scratch/static buffer (OD) bound alongside it. The barotropic mode
 * takes its own small RK4 with fixed-offset stage blocks inside OD,
 * looped on the host the M times a step needs; every other multi-stage
 * loop (the baroclinic RK4, the ∇⁴ closure's two Laplacian passes) mirrors
 * the dispatch-ordering-is-execution-ordering pattern the sigma core
 * relies on within one compute pass.
 *
 * Initialization, loading and serialization run once per model build and
 * are cheap relative to a step, so they stay in JavaScript at double
 * precision, exactly porting the CPU functions over plain arrays before
 * uploading; only the per-step dynamics (tendency, barotropic, eddy
 * transport, Ekman block, mixed layer, salt, surface write-back) are WGSL.
 */
const WORKGROUP = 64;
// Float32 values, so the host's conversions use the numbers the kernels do.
export const HEAT_REFERENCE = Math.fround(273.15), SALT_REFERENCE = Math.fround(35);

export const OCEAN_DEFAULTS = {
  densities: LAYER_DENSITIES, salinities: LAYER_SALINITIES, bottoms: LAYER_BOTTOMS, mixedDepth: 60, minimumDepth: 50, flatDepth: 4000, thermoclineTilt: 0.3,
  density: 1025, specificHeat: 3985, referenceS: 35, gravity: 9.81,
  minimumThickness: 50, shallowestMixedDepth: 50, stirringDepth: 100, maximumMixedDepth: 600, convectiveRate: 100 / 86400, neutralSnap: false, convectiveErosion: true, buoyancyMemory: 86400, mixedNeighbourRatio: 0, vorticityCentring: 0.5, stirring: 0.8, detrainmentTime: 86400, restoreTime: 2 * 86400, iceSalinity: 5, iceStressTransmission: 0.8, iceDensity: 917,
  momentumDepth: 50, carryMomentum: true, refillRate: 0.00027778, refillFloor: 5, shearEntrainment: 0.65,
  ekmanDepth: 0.4, momentumFloor: 10, stressThickness: PV_FLOOR, retreatStirring: 2.5, langmuir: 2, retreatEkman: RETREAT_EKMAN, retreatDaily: false, retreatLight: 0,
  interfacialDrag: 2e-4, shearMixing: false, interiorShearMixing: false, blockDrag: false, shearViscosity: 1e-2, backgroundViscosity: 1e-4, bottomDrag: 3e-3, closureHours: 12, closureSpacing: CLOSURE_SPACING, closureFill: 0, closureTokens: 'interior', diffusivity: 0.01, mixedDiffusivity: 0, everySteps: 4,
  eddyDiffusivity: 1000, eddyTaperDepth: 200, ekmanBlock: false, upwellingSource: 0, penetration: 'IB',
  gustiness: 3,
  generator: 'hybrid', regridTime: 0, minimumLayer: H_MIN, modeDepth: MODE_DEPTH, mstar: 1.25, ekmanFactor: 0.5, nstar: 0.2, backgroundDiffusivity: 1e-5, convectiveDiffusivity: 0.1,
};
const HYBRID_REFUSED = ['ekmanBlock', 'intakeLedger', 'retreatDaily', 'shearMixing', 'interiorShearMixing', 'blockDrag'];
const SLAB_ONLY_KERNELS = ['oSurfaceDensity', 'oGradRho', 'oDeepestEdge', 'oClosureFill', 'oClosureRing', 'oClosureBack1', 'oClosureBack2', 'oMixedReach', 'oMixedLayer', 'oExchangeMomentum', 'oShearEdge', 'oWindDepth', 'oEkmanFlux', 'oEkmanApply', 'oEkmanLeftover'];
const defaultSalinityProfile = (lat) => 34 + 2 * Math.exp(-(((Math.abs(lat) * 180 / Math.PI - 25) / 20) ** 2));

export function oceanStateLength(C, E, options = {}) {
  const o = { ...OCEAN_DEFAULTS, ...options };
  return (o.vertical === 'hybrid' ? (o.layerTable ?? HYBRID_LAYERS).K : o.densities.length + 1) * (3 * C + E);
}

function seq(names) { const out = {}; let off = 0; for (const [name, n] of names) { for (const alias of [].concat(name)) out[alias] = off; off += n; } out.total = off; return out; }

/*
 * The hybrid ocean's kernels (docs/ocean-hybrid-coordinate.md; the CPU's
 * hybrid path in ocean/layered.module.js and the column library in
 * ocean/column.module.js), used in place of the slab's under
 * `vertical: 'hybrid'`. Every per-column array is a local of the kernel's
 * entry point of at most L + 1 floats, and the column routines (the
 * tridiagonal solve, the limited reconstruction, the remap, the boundary
 * layer's work) are written into each kernel by the generators below
 * rather than called as WGSL functions, since Tint hands a function its
 * arrays by pointer and Apple's compilers then lose writes to them (the
 * forms in CLAUDE.md). Every loop over an array counts up; the solve's
 * back substitution runs over j with k = L − 2 − j and carries the value
 * it just wrote in a scalar, since Metal loses the read of an element
 * written in the previous iteration. The
 * remap gives each old layer's pieces to the new layers and its rounding
 * residual to the last of them, so a column's water, heat and salt are
 * kept to the stores' rounding. Heat and salt are the stores' anomalies
 * h·(T − TREF) and h·(S − SREF) throughout; only the equation of state
 * sees the absolute values.
 */
const declareArrays = (names, n = 'L') => names.map((name) => `var ${name}: array<f32, ${n}>;`).join(' ');
// solveColumn: the implicit diffusion of Q over the H layers with KAPPA at their interfaces, when any is positive.
const solveInline = ({ h, kappa, q, dt }, body = '') => `{
    var mixing = false;
    for (var k = 0; k < L - 1; k++) { if (${kappa}[k] > 0.0) { mixing = true; } }
    if (mixing) {
      for (var k = 0; k < L; k++) {
        SUP[k] = select(0.0, ${dt} * ${kappa}[max(k - 1, 0)] / max(0.5 * (${h}[max(k - 1, 0)] + ${h}[k]), DZFLOOR), k > 0);
        SLO[k] = select(0.0, ${dt} * ${kappa}[k] / max(0.5 * (${h}[k] + ${h}[min(k + 1, L - 1)]), DZFLOOR), k < L - 1);
      }
      var m = ${h}[0] + SUP[0] + SLO[0];
      SG[0] = -SLO[0] / m;
      SR[0] = ${h}[0] * ${q}[0] / m;
      for (var k = 1; k < L; k++) {
        m = ${h}[k] + SUP[k] + SLO[k] + SUP[k] * SG[k - 1];
        SG[k] = -SLO[k] / m;
        SR[k] = (${h}[k] * ${q}[k] + SUP[k] * SR[k - 1]) / m;
      }
      var next = SR[L - 1];
      ${q}[L - 1] = next;
      for (var j = 0; j < L - 1; j++) { let k = L - 2 - j; next = SR[k] - SG[k] * next; ${q}[k] = next; }
${body}    }
  }`;
// The limited reconstruction of RQ on the RH layers, RZO their interfaces: order 2 the monotonized central slope, order 3 the parabolas with the interface values from the quartic through the running integrals.
const reconstructInline = (order) => (order === 2 ? `  for (var k = 0; k < L; k++) {
    var s = 0.0;
    if (k > 0 && k < L - 1) {
      let below = RQ[min(k + 1, L - 1)] - RQ[k]; let above = RQ[k] - RQ[max(k - 1, 0)];
      if (below * above > 0.0) {
        let central = RH[k] * (RQ[min(k + 1, L - 1)] - RQ[max(k - 1, 0)]) / max(0.5 * (RH[max(k - 1, 0)] + 2.0 * RH[k] + RH[min(k + 1, L - 1)]), DZFLOOR);
        s = sign(central) * min(abs(central), min(2.0 * abs(below), 2.0 * abs(above)));
      }
    }
    RL[k] = RQ[k] - 0.5 * s; RS[k] = s; RC[k] = 0.0;
  }
` : `  RE[0] = RQ[0]; RE[L] = RQ[L - 1];
  for (var k = 1; k < L; k++) {
    let a = k - 1; let low = min(RQ[a], RQ[k]); let high = max(RQ[a], RQ[k]);
    var value = 0.0;
    let wide = k >= 2 && k <= L - 2 && RH[max(k - 2, 0)] > HNEGLECT && RH[a] > HNEGLECT && RH[k] > HNEGLECT && RH[min(k + 1, L - 1)] > HNEGLECT;
    if (wide) {
      RK[0] = RZO[k - 2] - RZO[k]; RK[1] = RZO[k - 1] - RZO[k]; RK[2] = 0.0; RK[3] = RZO[k + 1] - RZO[k]; RK[4] = RZO[k + 2] - RZO[k];
      RP[0] = 0.0; RP[1] = RH[k - 2] * RQ[k - 2]; RP[2] = RP[1] + RH[k - 1] * RQ[k - 1]; RP[3] = RP[2] + RH[k] * RQ[k]; RP[4] = RP[3] + RH[k + 1] * RQ[k + 1];
      var sum = 0.0;
      for (var p = 0; p < 5; p++) {
        if (p == 2) {
          var s = 0.0;
          for (var j = 0; j < 5; j++) { if (j != 2) { s += 1.0 / (RK[2] - RK[j]); } }
          sum += RP[2] * s;
        } else {
          var numerator = 1.0; var denominator = 1.0;
          for (var j = 0; j < 5; j++) {
            if (j == p) { continue; }
            denominator *= RK[p] - RK[j];
            if (j != 2) { numerator *= RK[2] - RK[j]; }
          }
          sum += RP[p] * numerator / denominator;
        }
      }
      value = sum;
    } else {
      value = (RH[k] * RQ[a] + RH[a] * RQ[k]) / max(RH[a] + RH[k], HNEGLECT);
    }
    RE[k] = clamp(value, low, high);
  }
  for (var k = 0; k < L; k++) {
    var qL = RE[k]; var qR = RE[k + 1]; let mean = RQ[k];
    if ((qR - mean) * (mean - qL) <= 0.0) { qL = mean; qR = mean; } else {
      let d = qR - qL; let excess = d * (mean - 0.5 * (qL + qR));
      if (excess > d * d / 6.0) { qL = 3.0 * mean - 2.0 * qR; } else if (excess < -d * d / 6.0) { qR = 3.0 * mean - 2.0 * qL; }
    }
    RL[k] = qL; RS[k] = qR - qL; RC[k] = 6.0 * (mean - 0.5 * (qL + qR));
  }
`);
const interfacesInline = (h, z) => `  ${z}[0] = 0.0;
  for (var k = 0; k < L; k++) { ${z}[k + 1] = ${z}[k] + ${h}[k]; }
`;
// remapColumn: RQ on the RH layers, RO their amounts, integrated over the RHN layers into RA, old layer by old layer.
const remapInline = (order) => `${interfacesInline('RH', 'RZO')}  let H = RZO[L];
${interfacesInline('RHN', 'RZN')}  RZN[L] = H;
${reconstructInline(order)}  for (var i = 0; i < L; i++) { RA[i] = 0.0; }
  var inew = 0;
  for (var j = 0; j < L; j++) {
    let hj = RH[j];
    if (!(hj > 0.0)) { continue; }
    let zLo = RZO[j]; let zHi = RZO[j + 1];
    while (inew < L - 1 && RZN[inew + 1] <= zLo) { inew++; }
    var given = 0.0; var last = inew;
    for (var ii = inew; ii < L; ii++) {
      let lo = max(zLo, RZN[ii]); let hi = min(zHi, RZN[ii + 1]);
      if (hi > lo) {
        let a = clamp((lo - zLo) / hj, 0.0, 1.0); let b = clamp((hi - zLo) / hj, 0.0, 1.0);
        let d1 = b - a; let d2 = 0.5 * (b * b - a * a); let d3 = (b * b * b - a * a * a) / 3.0;
        let p = hj * (RL[j] * d1 + RS[j] * d2 + RC[j] * (d2 - d3));
        RA[ii] += p; given += p; last = ii;
      }
      if (RZN[ii + 1] >= zHi) { break; }
    }
    RA[last] += RO[j] - given;
  }
`;
const REMAP_ARRAYS = ['RH', 'RQ', 'RO', 'RHN', 'RA', 'RL', 'RS', 'RC'], REMAP_INTERFACES = ['RZO', 'RZN', 'RE'], SOLVE_ARRAYS = ['SUP', 'SLO', 'SG', 'SR'];
// cellColumn: the cell's layers into CH, CT, CS, the tracers as the stores' anomalies.
const loadColumnInline = `  var tLoad = LABEL_T[0] - TREF; var sLoad = LABEL_S[0] - SREF;
  for (var k = 0; k < L; k++) {
    let hv = IN[hOff(k) + i];
    CH[k] = hv;
    if (hv > 1e-9) { tLoad = IN[qOff(k) + i] / hv; sLoad = IN[wOff(k) + i] / hv; }
    CT[k] = tLoad; CS[k] = sLoad;
  }
`;
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
${arr('DEPTHS', table.depths)}
${arr('PVFLOORK', floors)}
${arr('BOUNDARY', boundaries.subarray(1, L))}
const GEN: i32 = ${GENERATORS.indexOf(o.generator)};
const LP1: i32 = ${L + 1};
${constLine('HMIN', o.minimumLayer)} ${constLine('THINEDGE', 2 * o.minimumLayer)} ${constLine('MODEDEPTH', o.modeDepth)} ${constLine('MSTAR', o.mstar)} ${constLine('EKMANF', o.ekmanFactor)} ${constLine('NSTAR', o.nstar)}
${constLine('BACKKAPPA', o.backgroundDiffusivity)} ${constLine('CONVKAPPA', o.convectiveDiffusivity)} ${constLine('VONK', 0.4)} ${constLine('BLFLOOR', 0.1)} ${constLine('DZFLOOR', 1e-3)} ${constLine('HNEGLECT', 1e-6)} ${constLine('MLEXCESS', 0.03)}
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
fn blEnergy(depth: f32, ustar: f32, f: f32, dt: f32, convective: f32) -> f32 {
  let rotation = select(1.0e30, 1.0 + (depth * f / (EKMANF * ustar)) * (depth * f / (EKMANF * ustar)), ustar > 0.0);
  let factor = ${o.langmuir > 0 ? 'select(1.0, langmuirFactor(ustar, max(depth, 1.0)), ustar > 0.0)' : '1.0'};
  return RHO0 * dt * MSTAR * factor * ustar * ustar * ustar / rotation + NSTAR * convective;
}
`;
}

function hybridKernels(o, K, idx) {
  const { L, C, E } = o;
  const cellVector = (field) => `vec3<f32>(0.0, 0.0, 0.0);
  for (var m = 0; m < MAXE; m++) {
    let e = MI[EOC + MAXE * i + m];
    let s = abs(f32(MI[ESC + MAXE * i + m])) * 0.5 * MF[F_DC + e] * MF[F_DV + e] * ${field};
    wv += s * vec3<f32>(MF[F_NEDGE + 3 * e], MF[F_NEDGE + 3 * e + 1], MF[F_NEDGE + 3 * e + 2]);
  }
  wv = wv / MF[F_AREA + i];`;
  // boundaryLayer's work: the potential energy the diffusivity to `depth` raises in a step, the diffusivity left in CK and the trial tracers in CTT, CTS.
  const blWork = `      for (var k = 0; k < L - 1; k++) { let z = CI[k + 1]; CK[k] = select(0.0, VONK * vstar * z * (1.0 - z / depth) * (1.0 - z / depth), z < depth); }
      CK[L - 1] = 0.0;
      for (var k = 0; k < L; k++) { CTT[k] = CT[k]; CTS[k] = CS[k]; }
      ${solveInline({ h: 'CH', kappa: 'CK', q: 'CTT', dt: 'dt' })}
      ${solveInline({ h: 'CH', kappa: 'CK', q: 'CTS', dt: 'dt' })}
      var raised = 0.0;
      for (var k = 0; k < L - 1; k++) {
        let dz = max(0.5 * (CH[k] + CH[k + 1]), DZFLOOR); let lift = CC[k + 1] - CC[k];
        let heatDown = dt * CK[k] * (CTT[k] - CTT[k + 1]) / dz; let saltDown = dt * CK[k] * (CTS[k] - CTS[k + 1]) / dz;
        raised += OGRAV * RHO0 * lift * (CA[k] * heatDown - CB[k] * saltDown);
      }
`;
  return {
    /*
     * prepareCells: each cell's layer densities as anomalies and the
     * running sums of ρ'h and h above each layer.
     */
    oLayerDensity: `${K}  let i = ${idx}; if (i >= C) { return; }
  var weighted = 0.0; var total = 0.0;
  for (var k = 0; k < L; k++) {
    let n = k * C + i; let hv = IN[hOff(k) + i];
    let r = select(RHOA[k], eosAnomaly(IN[qOff(k) + i] / hv + TREF, IN[wOff(k) + i] / hv + SREF), hv > 1e-6);
    OD[O_RHOL + n] = r; OD[O_PRESS + n] = weighted; OD[O_HABOVE + n] = total;
    weighted += r * hv; total += hv;
  }
}`,
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
     * MOMDEPTH of the edge column's floor, the ∇⁴ closure, and the pin of
     * a layer no thicker than THINEDGE to the one above; P[2] is relax.
     */
    oMomentum: `${K}  let e = ${idx}; if (e >= E) { return; }
  if (OD[O_EMASK + e] < 0.5) { for (var k = 0; k < L; k++) { OUT[uOff(k) + e] = 0.0; } return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  let dc = MF[F_DC + e];
  var others: array<i32, MAXEE>; var weights: array<f32, MAXEE>; var live: array<bool, MAXEE>;
  for (var s = 0; s < MAXEE; s++) { let slot = MAXEE * e + s; others[s] = MI[EOE + slot]; weights[s] = MF[F_PVW + slot]; live[s] = OD[O_EMASK + others[s]] > 0.5; }
  var floor = 0.0;
  for (var k = 0; k < L; k++) { floor += OD[O_HEDGE + k * E + e]; }
  var layerTop = 0.0; var uAbove = 0.0;
  for (var k = 0; k < L; k++) {
    let n = k * E + e;
    let u = IN[uOff(k) + e];
    let qHere = 0.5 * OD[O_QE + n];
    var pv = 0.0;
    for (var s = 0; s < MAXEE; s++) { let other = others[s]; let fo = select(0.0, OD[O_HEDGE + k * E + other] * IN[uOff(k) + other], live[s]); pv += weights[s] * fo * (qHere + 0.5 * OD[O_QE + k * E + other]); }
    let gradK = (OD[O_PHI + k * C + b] - OD[O_PHI + k * C + a]) / dc;
    let na = k * C + a; let nb = k * C + b;
    let Ha = OD[O_HABOVE + na]; let Hb = OD[O_HABOVE + nb]; let ha = IN[hOff(k) + a]; let hb = IN[hOff(k) + b];
    let top = max(Ha, Hb); let bottom = min(Ha + ha, Hb + hb);
    let zR = select(0.5 * (Ha + 0.5 * ha + Hb + 0.5 * hb), 0.5 * (top + bottom), bottom > top);
    let force = OD[O_PRESS + nb] + OD[O_RHOL + nb] * (zR - Hb) - OD[O_PRESS + na] - OD[O_RHOL + na] * (zR - Ha);
    var du = pv / dc - gradK - OGRAV * OD[O_GRADETA + e] - OGRAV / RHO0 * force / dc;
    let he = OD[O_HEDGE + n];
    if (k == 0) { du += OD[O_STRESS + e] / RHO0 / max(he, HMIN); }
    let overlap = min(layerTop + he, floor) - max(layerTop, floor - MOMDEPTH);
    if (overlap > 0.0) { du -= RBOT * abs(u) * u * (overlap / MOMDEPTH) / max(he, HMIN); }
    if (NU4O > 0.0) { du -= NU4O * OD[O_LAPB + n]; }
    if (k > 0 && he <= THINEDGE) { du = (uAbove - u) * P[2]; }
    OUT[uOff(k) + e] = du;
    uAbove = u; layerTop += he;
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
     * the regrid's relax.
     */
    oRegridTargets: `${K}  let i = ${idx}; if (i >= C) { return; }
  if (OD[O_CMASK + i] < 0.5) { return; }
  ${declareArrays(['CH', 'CT', 'CS', 'RH', 'RQ', 'RL', 'RS', 'RC', 'FL'])} ${declareArrays(['RZO', 'RZN', 'RE'], 'LP1')} ${declareArrays(['RK', 'RP'], 5)}
${loadColumnInline}  for (var k = 0; k < L; k++) { RH[k] = CH[k]; RQ[k] = eosA(CT[k], CS[k]); }
${interfacesInline('RH', 'RZO')}  let H = RZO[L];
  let depth = OD[O_BATH + i];
  RZN[0] = 0.0; RZN[L] = H;
  if (GEN == 1) {
    for (var k = 1; k < L; k++) { RZN[k] = min(DEPTHS[k], depth) * H / depth; }
  } else {
    for (var k = 0; k < L; k++) {
      let below = select(0.0, select(clamp((RH[k] - MINIMUM[k]) / MODEDEPTH, 0.0, 1.0), 1.0, GEN == 0), k + 1 < L);
      let above = select(0.0, select(clamp((RH[max(k - 1, 0)] - MINIMUM[max(k - 1, 0)]) / MODEDEPTH, 0.0, 1.0), 1.0, GEN == 0), k > 0);
      FL[k] = min(above, below);
    }
    for (var k = 1; k < L; k++) { RQ[k] = max(RQ[k - 1], RQ[k]); }
${reconstructInline(3)}    for (var k = 0; k < L; k++) { let keep = 1.0 - FL[k]; RL[k] = RQ[k] + keep * (RL[k] - RQ[k]); RS[k] *= keep; RC[k] *= keep; }
    var j = 0;
    for (var k = 1; k < L; k++) {
      let bk = BOUNDARY[k - 1];
      while (j < L && RL[min(j, L - 1)] + RS[min(j, L - 1)] < bk) { j++; }
      if (j >= L) { RZN[k] = H; continue; }
      if (bk <= RL[j]) { RZN[k] = RZO[j]; continue; }
      RZN[k] = RZO[j] + RH[j] * crossing(RL[j], RS[j], RC[j], bk);
    }
    var above = 0.0;
    for (var k = 1; k < L; k++) {
      let floor = select(0.0, above + MINIMUM[k - 1], GEN == 2);
      var wanted = max(RZN[k], floor);
      if (wanted > floor) { wanted = max(floor, RZO[k] + P[0] * (wanted - RZO[k])); }
      RZN[k] = wanted;
      above = wanted;
    }
  }
  let spacing = min(HMIN, H / f32(L));
  for (var k = 1; k < L; k++) { RZN[k] = max(RZN[k], RZN[k - 1] + spacing); }
  for (var j = 0; j < L - 1; j++) { let k = L - 1 - j; RZN[k] = min(RZN[k], H - f32(L - k) * spacing); }
  for (var k = 0; k <= L; k++) { OD[O_TARGET + k * C + i] = RZN[k]; }
}`,
    // regridEdges: each edge's velocity remapped piecewise-linearly onto the mean of its cells' new thicknesses.
    oRegridEdges: `${K}  let e = ${idx}; if (e >= E) { return; }
  if (OD[O_EMASK + e] < 0.5) { return; }
  ${declareArrays(REMAP_ARRAYS)} ${declareArrays(REMAP_INTERFACES, 'LP1')}
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  for (var k = 0; k < L; k++) {
    RH[k] = 0.5 * (IN[hOff(k) + a] + IN[hOff(k) + b]);
    RHN[k] = 0.5 * (OD[O_TARGET + (k + 1) * C + a] - OD[O_TARGET + k * C + a] + OD[O_TARGET + (k + 1) * C + b] - OD[O_TARGET + k * C + b]);
    RQ[k] = IN[uOff(k) + e]; RO[k] = RH[k] * RQ[k];
  }
${remapInline(2)}  for (var k = 0; k < L; k++) { IN[uOff(k) + e] = select(RQ[k], RA[k] / RHN[k], RHN[k] > 0.0); }
}`,
    // regridCells: heat and salt remapped by the parabolas onto the targets, which set the thicknesses.
    oRegridCells: `${K}  let i = ${idx}; if (i >= C) { return; }
  if (OD[O_CMASK + i] < 0.5) { return; }
  ${declareArrays(['CH', 'CT', 'CS', ...REMAP_ARRAYS])} ${declareArrays(REMAP_INTERFACES, 'LP1')} ${declareArrays(['RK', 'RP'], 5)}
${loadColumnInline}  for (var k = 0; k < L; k++) { RH[k] = CH[k]; RHN[k] = OD[O_TARGET + (k + 1) * C + i] - OD[O_TARGET + k * C + i]; RQ[k] = CT[k]; RO[k] = IN[qOff(k) + i]; }
  {
${remapInline(3)}  }
  for (var k = 0; k < L; k++) { IN[qOff(k) + i] = RA[k]; RQ[k] = CS[k]; RO[k] = IN[wOff(k) + i]; }
  {
${remapInline(3)}  }
  for (var k = 0; k < L; k++) { IN[wOff(k) + i] = RA[k]; IN[hOff(k) + i] = RHN[k]; }
}`,
    /*
     * columnCells: mixColumn of ocean/column.module.js on each sea column
     * with the cell's stress and the shear of its reconstructed velocities:
     * convective adjustment, the energetic boundary layer by bisection, the
     * Pacanowski–Philander interior, one implicit solve of heat and salt;
     * the viscosity profile into VISC for the edges, the mixed-layer and
     * boundary-layer depths into MIXED and BLD; P[0] is the step. The
     * bisection's loop evaluates the work at the whole column first (it
     * -1), then halves the interval 24 times, then once more at the depth
     * it settled on.
     */
    oColumnCells: `${K}  let i = ${idx}; if (i >= C) { return; }
  if (OD[O_CMASK + i] < 0.5) { OD[O_MIXED + i] = 0.0; OD[O_BLD + i] = 0.0; return; }
  ${declareArrays(['CH', 'CT', 'CS', 'CQ', 'CW', 'CSH', 'CK', 'CIK', 'CINU', 'CC', 'CA', 'CB', 'CTT', 'CTS', 'CR', 'CTH', 'CHE', 'CSA', 'CSR', ...SOLVE_ARRAYS])} ${declareArrays(['CI'], 'LP1')} var CST: array<i32, LP1>;
  let dt = P[0];
${loadColumnInline}  for (var k = 0; k < L; k++) { CQ[k] = IN[qOff(k) + i]; CW[k] = IN[wOff(k) + i]; }
  var previous = vec3<f32>(0.0, 0.0, 0.0);
  for (var k = 0; k < L; k++) {
    var wv = ${cellVector('IN[uOff(k) + e]')}
    if (k > 0) { let d = previous - wv; CSH[k - 1] = dot(d, d); }
    previous = wv;
  }
  var wv = ${cellVector('OD[O_STRESS + e]')}
  let ustar = sqrt(length(wv) / RHO0); let f = abs(OD[O_FCELL + i]);
  var z = 0.0;
  for (var k = 0; k < L; k++) { CC[k] = z + 0.5 * CH[k]; z += CH[k]; CR[k] = eosA(CT[k], CS[k]); }
  var n = 0;
  for (var k = 0; k < L; k++) {
    CST[n] = k; CTH[n] = CH[k]; CHE[n] = CQ[k]; CSA[n] = CW[k]; CSR[n] = CR[k];
    n++;
    loop {
      if (!(n >= 2 && CSR[n - 2] > CSR[n - 1])) { break; }
      CTH[n - 2] += CTH[n - 1]; CHE[n - 2] += CHE[n - 1]; CSA[n - 2] += CSA[n - 1];
      CSR[n - 2] = select(CSR[n - 2], eosA(CHE[n - 2] / CTH[n - 2], CSA[n - 2] / CTH[n - 2]), CTH[n - 2] > 0.0);
      n--;
    }
  }
  CST[n] = L;
  var released = 0.0;
  for (var s = 0; s < n; s++) {
    let first = CST[s]; let past = CST[s + 1];
    if (past - first < 2) { continue; }
    let t = CHE[s] / CTH[s]; let w = CSA[s] / CTH[s]; let mixed = eosA(t, w);
    var drop = 0.0; var sumQ = 0.0; var sumW = 0.0;
    for (var k = first; k < past; k++) { drop += OGRAV * CC[k] * CH[k] * (mixed - CR[k]); CT[k] = t; CS[k] = w; CQ[k] = CH[k] * t; CW[k] = CH[k] * w; sumQ += CQ[k]; sumW += CW[k]; }
    CQ[past - 1] += CHE[s] - sumQ; CW[past - 1] += CSA[s] - sumW;
    if (first == 0) { released = drop; }
  }
${interfacesInline('CH', 'CI')}  let H = CI[L];
  for (var k = 0; k < L - 1; k++) {
    CA[k] = 0.5 * (alphaT(CT[k] + TREF, CS[k] + SREF) + alphaT(CT[k + 1] + TREF, CS[k + 1] + SREF));
    CB[k] = 0.5 * (betaS(CT[k] + TREF, CS[k] + SREF) + betaS(CT[k + 1] + TREF, CS[k + 1] + SREF));
  }
  let convective = max(0.0, released); let vstar = pow(ustar * ustar * ustar + convective / (RHO0 * dt), 1.0 / 3.0);
  for (var k = 0; k < L; k++) { CK[k] = 0.0; }
  var depth = 0.0;
  if (vstar > 0.0) {
    depth = H;
    var low = log(min(BLFLOOR, H)); var high = log(H); var searching = false;
    for (var it = -1; it <= 24; it++) {
      if (it >= 0 && !searching) { break; }
      let middle = 0.5 * (low + high);
      if (it >= 0 && it < 24) { depth = exp(middle); } else if (it == 24) { depth = exp(low); }
${blWork}      if (it == -1) { searching = raised > blEnergy(depth, ustar, f, dt, convective); } else if (it < 24) { if (raised <= blEnergy(depth, ustar, f, dt, convective)) { low = middle; } else { high = middle; } }
    }
    if (searching && depth <= BLFLOOR) { depth = 0.0; for (var k = 0; k < L; k++) { CK[k] = 0.0; } }
  }
  for (var k = 0; k < L; k++) { CR[k] = eosA(CT[k], CS[k]); }
  for (var k = 0; k < L - 1; k++) {
    let dz = max(0.5 * (CH[k] + CH[k + 1]), DZFLOOR); let buoyancy = OGRAV * (CR[k + 1] - CR[k]) / RHO0;
    if (buoyancy < 0.0) { CIK[k] = CONVKAPPA; CINU[k] = CONVKAPPA; continue; }
    let richardson = buoyancy * dz / (CSH[k] + 1e-12);
    let viscosity = SHEARNU / ((1.0 + 5.0 * richardson) * (1.0 + 5.0 * richardson)) + BACKNU;
    CINU[k] = viscosity;
    CIK[k] = viscosity / (1.0 + 5.0 * richardson) + BACKKAPPA;
  }
  CIK[L - 1] = 0.0; CINU[L - 1] = 0.0;
  for (var k = 0; k < L; k++) { CIK[k] += CK[k]; CINU[k] += CK[k]; }
  ${solveInline({ h: 'CH', kappa: 'CIK', q: 'CT', dt: 'dt' }, fluxInline({ h: 'CH', kappa: 'CIK', q: 'CT', amount: 'CQ', dt: 'dt' }))}
  ${solveInline({ h: 'CH', kappa: 'CIK', q: 'CS', dt: 'dt' }, fluxInline({ h: 'CH', kappa: 'CIK', q: 'CS', amount: 'CW', dt: 'dt' }))}
  for (var k = 0; k < L; k++) {
    let m = k * C + i;
    IN[qOff(k) + i] = CQ[k]; IN[wOff(k) + i] = CW[k]; OD[O_VISC + m] = CINU[k];
    CR[k] = eosA(CT[k], CS[k]);
  }
  var mixed = H; var zc = 0.5 * CH[0]; var prior = CR[0]; let limit = CR[0] + MLEXCESS; var found = false;
  for (var k = 1; k < L; k++) {
    if (found) { continue; }
    let center = zc + 0.5 * (CH[k - 1] + CH[k]);
    if (CR[k] > limit) { mixed = zc + (limit - prior) / (CR[k] - prior) * (center - zc); found = true; }
    zc = center; prior = CR[k];
  }
  OD[O_MIXED + i] = mixed; OD[O_BLD + i] = depth;
}`,
    // columnEdges: the implicit solve of each edge's velocity with the mean of its cells' viscosity profiles.
    oColumnEdges: `${K}  let e = ${idx}; if (e >= E) { return; }
  if (OD[O_EMASK + e] < 0.5) { return; }
  ${declareArrays(['SH', 'SK', 'SQ', ...SOLVE_ARRAYS])}
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  for (var k = 0; k < L; k++) { SH[k] = 0.5 * (IN[hOff(k) + a] + IN[hOff(k) + b]); SK[k] = 0.5 * (OD[O_VISC + k * C + a] + OD[O_VISC + k * C + b]); SQ[k] = IN[uOff(k) + e]; }
  ${solveInline({ h: 'SH', kappa: 'SK', q: 'SQ', dt: 'P[0]' }, '      for (var k = 0; k < L; k++) { IN[uOff(k) + e] = SQ[k]; }\n')}
}`,
  };
}

function oceanKernels(o) {
  const { L, C, E, V, OS, OD, B } = o, hy = !!o.hybrid;
  if (L > 64) throw new Error(`oMomentum marks an edge's thick classes in 64 bits, and the ocean has ${L} classes`);
  const constLine = (name, value) => `const ${name}: f32 = ${Number(value).toExponential(10)};`;
  const jerlov = JERLOV[o.penetration ?? 'IB'];
  const rhoLine = `const RHO: array<f32, ${L}> = array<f32, ${L}>(${o.rho.map((v) => v.toFixed(6)).join(', ')});`;
  const spread = (sign) => o.rho.map((r, k) => r + sign * 0.5 * (sign < 0 ? (k > 1 ? r - o.rho[k - 1] : o.rho[k + 1] - r) : (k < L - 1 ? o.rho[k + 1] - r : r - o.rho[k - 1])));
  const classLine = `const LIGHTEST: array<f32, ${L}> = array<f32, ${L}>(${spread(-1).map((v) => v.toFixed(6)).join(', ')});\nconst DENSEST: array<f32, ${L}> = array<f32, ${L}>(${spread(1).map((v) => v.toFixed(6)).join(', ')});`;
  const labelLine = `const LABEL_T: array<f32, ${L}> = array<f32, ${L}>(${o.labelT.map((v) => v.toFixed(6)).join(', ')});\nconst LABEL_S: array<f32, ${L}> = array<f32, ${L}>(${o.labelS.map((v) => v.toFixed(6)).join(', ')});\nconst RHOA: array<f32, ${L}> = array<f32, ${L}>(${o.rho.map((v) => (v - SEAWATER.rho0).toFixed(6)).join(', ')});`;
  const { a0, b0, lambda1, lambda2, nu, t0, s0 } = SEAWATER;
  const shearDensity = `fn heldAnomaly(q: f32, w: f32) -> f32 {
  let ta = q + ${(HEAT_REFERENCE - t0).toExponential(10)}; let sa = w + ${(SALT_REFERENCE - s0).toExponential(10)};
  return -${a0} * (1.0 + ${0.5 * lambda1} * ta) * ta + ${b0} * (1.0 - ${0.5 * lambda2} * sa) * sa - ${nu} * ta * sa;
}`;
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
  const mixingLines = [
    o.ekmanDepth > 0 ? `${constLine('EKMAN', o.ekmanDepth)} ${constLine('MOMFLOOR', o.momentumFloor)}` : '',
    o.stressThickness !== PV_FLOOR ? constLine('STRESSTHICK', o.stressThickness) : '',
    o.retreatStirring > 0 ? constLine('RETREATSTIR', o.retreatStirring) : '',
    o.langmuir > 0 ? `${constLine('LANGMUIR', o.langmuir)}\n${langmuirLines}` : '',
    o.retreatEkman > 0 ? `${constLine('RETREATEKMAN', o.retreatEkman)} ${constLine('TWOOMEGA', 2 * o.omega)}` : '',
    o.upwellingSource > 0 ? constLine('UPSRC', o.upwellingSource) : '',
    o.retreatLight > 0 ? `${constLine('RETREATLIGHT', o.retreatLight)}\nfn deposited(z: f32) -> f32 { return select(0.0, JER_R1 / JER_Z1 * exp(-z / JER_Z1) + JER_R2 / JER_Z2 * exp(-z / JER_Z2), z < PENDEPTH); }` : '',
  ].filter((line) => line).map((line) => `${line}\n`).join('');
  const spanLine = `  let dz = max(THINO, 0.5 * (OD[O_HEDGE + up * E + e] + OD[O_HEDGE + down * E + e]));\n`;
  const viscosityLines = (dv) => `  let upper = select(RHO[up], 0.5 * (OD[O_RHOML + MI[COE + 2 * e]] + OD[O_RHOML + MI[COE + 2 * e + 1]]), up == 0);
  let buoyancy = max(0.0, OGRAV * (RHO[down] - upper) / RHO0);
  let du = IN[uOff(up) + e] - IN[uOff(down) + e]; let dv = ${dv};
  let ri = buoyancy * dz / (du * du + dv * dv + 1e-12);
  let nu = SHEARNU / ((1.0 + 5.0 * ri) * (1.0 + 5.0 * ri)) + BACKNU;
`;
  const rateFunction = o.blockDrag ? `fn depthTo(e: i32, first: i32, past: i32, start: f32) -> f32 {
  var z = start;
  for (var j = first; j < past; j++) { z += OD[O_HEDGE + j * E + e]; }
  return z;
}
fn interfaceRate(e: i32, up: i32, down: i32, dz: f32, held: f32, tangentUp: f32, tangentDown: f32) -> f32 {
${viscosityLines('tangentUp - tangentDown')}  let viscous = min(nu / dz, 0.5 * min(dragThickness(up, e), dragThickness(down, e)) / P[3]);
  return viscous + (RINT - viscous) * held;
}
` : `fn interfaceRate(e: i32, up: i32, down: i32) -> f32 {
  let interior = up > 0 && INTSHEAR > 0.5;
  if (SHEARMIX < 0.5 && !interior) { return RINT; }
${spanLine}${viscosityLines('tangentialU(up, e) - tangentialU(down, e)')}  return min(max(select(RINT, 0.0, interior), nu / dz), 0.5 * min(dragThickness(up, e), dragThickness(down, e)) / P[3]);
}
`;
  const blockTop = o.blockDrag ? `  let block = ${o.ekmanDepth > 0 ? 'OD[O_HM + e]' : 'MOMDEPTH'};\n  var zTop = 0.0; var rateAbove = 0.0; var tangent = 0.0; var tangentOf = -1;\n` : '';
  const dragLines = o.blockDrag ? `    let below = firstThickBelow(k, thickLo, thickHi);
    if (k == 0 || hk >= THINO) {
      if (k > 0) { drag += rateAbove * (IN[uOff(above) + e] - u); }
      if (below >= 0) {
        let span = max(THINO, 0.5 * (hk + OD[O_HEDGE + below * E + e]));
        let held = min(1.0, (block - depthTo(e, k + 1, below, zTop + hk)) / span);
        if (held >= 1.0) {
          rateAbove = RINT;
        } else {
          if (tangentOf != k) { tangent = tangentialU(k, e); }
          let tangentBelow = tangentialU(below, e);
          rateAbove = interfaceRate(e, k, below, span, max(held, 0.0), tangent, tangentBelow);
          tangent = tangentBelow; tangentOf = below;
        }
        drag -= rateAbove * (u - IN[uOff(below) + e]);
      }
    }
` : `    if (k > 0) { drag += interfaceRate(e, above, k) * (IN[uOff(above) + e] - u); }
    let below = firstThickBelow(k, thickLo, thickHi);
    if (below >= 0) { drag -= interfaceRate(e, k, below) * (u - IN[uOff(below) + e]); }
`;
  const topStep = o.blockDrag ? '    zTop += hk;\n' : '';
  const momentumReach = o.ekmanDepth > 0 ? 'hm' : 'MOMDEPTH', stressFloor = o.stressThickness !== PV_FLOOR ? 'STRESSTHICK' : 'PVFLOOR';
  const retreatLoss = o.retreatDaily ? 'dayDrive' : 'loss';
  const unboundedRetreat = o.retreatStirring > 0 || o.langmuir > 0
    ? `${o.retreatStirring > 0 ? 'RETREATSTIR' : '2.0 * stir'} * ${o.langmuir > 0 ? '(1.0 + (langmuirFactor(sqrt(tau / RHO0), began) - 1.0) * openSea)' : '1.0'} * ustar3 / -${retreatLoss}`
    : `2.0 * stir * ustar3 / -${retreatLoss}`;
  const dailySums = `  let lightBuoy = OGRAV * alphaT(OD[O_SURFACEIN + i], s0) / RHOCP;
  let sunCall = OD[O_SUNSTEP + i];
  var sunSum = OD[O_SUNSUM + i] + sunCall;
  var darkSum = OD[O_DARKSUM + i] + buoyancy * P[6] + lightBuoy * sunCall * (1.0 - PH[PH_TRANS + i]);
  if (P[7] != 0.0) {
    if (P[7] > 0.0) { OD[O_SUNDAY + i] = sunSum * P[7]; OD[O_DARKDAY + i] = darkSum * P[7]; }
    sunSum = 0.0; darkSum = 0.0;
  }
  OD[O_SUNSUM + i] = sunSum; OD[O_DARKSUM + i] = darkSum;
`;
  const dailyDrive = `  let dayDepth = IN[hOff(0) + i];
${o.penetration ? `  var dayReceiver = false;
  for (var k = 1; k < L; k++) { if (IN[hOff(k) + i] > THINO) { dayReceiver = true; break; } }
  let dayPassing = select(0.0, transmitted(dayDepth) - transmitted(OD[O_BATH + i] + OD[O_ETA + i]), dayReceiver);
` : '  let dayPassing = 0.0;\n'}  let daySun = lightBuoy * OD[O_SUNDAY + i];
  let dayDrive = OD[O_DARKDAY + i] - daySun * (1.0 - dayPassing);
  let dayGate = ${o.retreatLight > 0 ? 'select(dayDrive, dayDrive + RETREATLIGHT * daySun * dayDepth * deposited(dayDepth), dayPassing > 0.0)' : 'dayDrive'};
`;
  const upwellingFunctions = o.upwellingSource > 0 ? `fn upwellingBudget(i: i32) -> vec2<f32> {
  var coast = false; var divergence = 0.0; var wind = 0.0; var sea = 0.0;
  for (var m = 0; m < MAXE; m++) {
    let sign = MI[ESC + MAXE * i + m];
    if (sign == 0) { continue; }
    let e = MI[EOC + MAXE * i + m];
    if (OD[O_CMASK + MI[COC + MAXE * i + m]] < 0.5) { coast = true; }
    divergence += f32(sign) * OD[O_UEK + e] * MF[F_DV + e];
    if (OD[O_EMASK + e] > 0.5) { wind += ${o.ekmanDepth > 0 ? 'OD[O_HM + e]' : 'MOMDEPTH'}; sea += 1.0; }
  }
  if (!coast || sea < 0.5) { return vec2<f32>(0.0, 0.0); }
  return vec2<f32>(UPSRC * max(0.0, divergence / MF[F_AREA + i]) * P[6], wind / sea);
}
fn upwellingSource(i: i32, wind: f32) -> i32 {
  let depth = max(wind, IN[hOff(0) + i]);
  var top = IN[hOff(0) + i];
  for (var k = 1; k < L; k++) {
    let hk = IN[hOff(k) + i];
    if (hk > THINO && top + hk > depth) { return k; }
    top += hk;
  }
  return -1;
}
` : '';
  const drawFromSource = (want) => `var drawn = 0.0;
      if (upLeft > 0.0 && ${want} > 0.0) {
        let source = upwellingSource(i, upwelling.y);
        if (source > 0) {
          let part = min(min(${want}, upLeft), IN[hOff(source) + i] - EPSO);
          if (part > 0.0) { moveLayer(i, source, 0, part); upLeft -= part; drawn = part; }
        }
      }`;
  const retreatTerm = o.retreatEkman > 0 ? `min(${unboundedRetreat}, select(1.0e30, RETREATEKMAN * sqrt(tau / RHO0) / fCell, fCell > 0.0))` : unboundedRetreat;
  const offsetLines = Object.entries(OD).filter(([k]) => k !== 'total').map(([k, v]) => `const O_${k}: i32 = ${v};`).join('\n');
  const bLines = Object.entries(B).filter(([k]) => k !== 'total').map(([k, v]) => `const B_${k}: i32 = ${v};`).join('\n');
  const head = `
const L: i32 = ${L}; const DRY: f32 = -1.0e30;
const OH: i32 = ${OS.OH}; const OU: i32 = ${OS.OU}; const OQ: i32 = ${OS.OQ}; const OW: i32 = ${OS.OW}; const OSTOTAL: i32 = ${OS.total};
${offsetLines}
${bLines}
${rhoLine}
${classLine}
${labelLine}
${constLine('RHO0', o.density)} ${constLine('RHOCP', o.density * o.specificHeat)}
${constLine('OGRAV', o.gravity)}
${constLine('EPSO', EPS)} ${constLine('THINO', THIN)} ${constLine('PVFLOOR', PV_FLOOR)} ${constLine('SPEEDLIM', SPEED_LIMIT)} ${constLine('DENSTOL', DENSITY_TOLERANCE)} ${constLine('RESTTOL', RESTORE_TOLERANCE)} ${constLine('RESTORET', o.restoreTime)}
${o.shearEntrainment > 0 ? `${constLine('RICRIT', o.shearEntrainment)}\n${shearDensity}\n` : ''}${mixingLines}${constLine('MINTHICK', o.minimumThickness)} ${constLine('REFILLRATE', o.refillRate)} ${constLine('REFILLFLOOR', o.refillFloor)} ${constLine('MOMDEPTH', o.momentumDepth)} ${constLine('SHALLOWMIXED', o.shallowestMixedDepth)} ${constLine('MAXMIXED', o.maximumMixedDepth)} ${constLine('CONVRATE', o.convectiveRate)}
${constLine('NEUTRALSNAP', o.neutralSnap ? 1 : 0)} ${constLine('EROSION', o.convectiveErosion ? 1 : 0)} ${constLine('BUOYMEM', o.buoyancyMemory)} ${constLine('NBRRATIO', o.mixedNeighbourRatio)} ${constLine('CENTRING', o.vorticityCentring)}
${constLine('STIRRING', o.stirring)} ${constLine('STIRDEPTH', o.stirringDepth)} ${constLine('DETRAINT', o.detrainmentTime)} ${constLine('ICESAL', o.iceSalinity)} ${constLine('TRANSMIT', o.iceStressTransmission)} ${constLine('ICEDENS', o.iceDensity)}
${constLine('RINT', o.interfacialDrag)} ${constLine('SHEARMIX', o.shearMixing ? 1 : 0)} ${constLine('INTSHEAR', o.interiorShearMixing ? 1 : 0)} ${constLine('SHEARNU', o.shearViscosity)} ${constLine('BACKNU', o.backgroundViscosity)} ${constLine('RBOT', o.bottomDrag)} ${constLine('NU4O', o.nu4)} ${constLine('DIFFUSION', o.diffusion)} ${constLine('MIXDIFF', o.mixedDiffusivity)}
${constLine('FREEZE', FREEZING_POINT)} ${constLine('TREF', HEAT_REFERENCE)} ${constLine('SREF', SALT_REFERENCE)} ${constLine('GUSTO', o.gustiness)} ${constLine('CLOSURERIDGE', CLOSURE_RIDGE)} ${constLine('CLOSUREFILL', o.closureFill || 0)} ${constLine('CLOSURERINGS', o.closureRings ? 1 : 0)}
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
fn iceVolume(i: i32) -> f32 {
  let h = OD[O_SURFICE + i]; let a = OD[O_SURFCONC + i];
  return select(0.0, h * select(1.0, a, a > 0.0), h > 0.0);
}
${SEAWATER_WGSL}
fn detrain(i: i32, amount: f32, rm: f32) {
  if (amount <= 0.0) { return; }
  var k = 1;
  for (var j = 2; j < L; j++) { if (abs(RHO[j] - rm) < abs(RHO[k] - rm)) { k = j; } }
  moveLayer(i, 0, k, amount);
}
struct ClosureFit { saa: f32, sab: f32, sbb: f32, sau: f32, sbu: f32, found: bool }
fn closureAdd(fit: ClosureFit, o: i32, nv: vec3<f32>, tv: vec3<f32>, value: f32) -> ClosureFit {
  let no = vec3<f32>(MF[F_NEDGE + 3 * o], MF[F_NEDGE + 3 * o + 1], MF[F_NEDGE + 3 * o + 2]);
  let a = dot(no, nv); let b = dot(no, tv);
  return ClosureFit(fit.saa + a * a, fit.sab + a * b, fit.sbb + b * b, fit.sau + a * value, fit.sbu + b * value, true);
}
fn closureWeight(k: i32, e: i32, o: i32, usableBelow: f32) -> f32 {
  let nv = vec3<f32>(MF[F_NEDGE + 3 * e], MF[F_NEDGE + 3 * e + 1], MF[F_NEDGE + 3 * e + 2]);
  let tv = vec3<f32>(OD[O_TEDGE + 3 * e], OD[O_TEDGE + 3 * e + 1], OD[O_TEDGE + 3 * e + 2]);
  var saa = 0.0; var sab = 0.0; var sbb = 0.0;
  for (var s = 0; s < MI[NEE + e]; s++) {
    let x = MI[EOE + MAXEE * e + s]; let flag = OD[O_FLUX + k * E + x];
    if (flag < 0.5 || flag > usableBelow) { continue; }
    let nx = vec3<f32>(MF[F_NEDGE + 3 * x], MF[F_NEDGE + 3 * x + 1], MF[F_NEDGE + 3 * x + 2]);
    let a = dot(nx, nv); let b = dot(nx, tv);
    saa += a * a; sab += a * b; sbb += b * b;
  }
  let no = vec3<f32>(MF[F_NEDGE + 3 * o], MF[F_NEDGE + 3 * o + 1], MF[F_NEDGE + 3 * o + 2]);
  let p = saa + CLOSURERIDGE; let q = sbb + CLOSURERIDGE;
  return (q * dot(no, nv) - sab * dot(no, tv)) / (p * q - sab * sab);
}
fn closureSolve(fit: ClosureFit) -> f32 {
  let p = fit.saa + CLOSURERIDGE; let q = fit.sbb + CLOSURERIDGE;
  return (q * fit.sau - fit.sab * fit.sbu) / (p * q - fit.sab * fit.sab);
}
fn grain(x: f32) -> f32 { return ldexp(1.0, max(frexp(x).exp - 24, -126)); }
fn settle(x: f32, change: f32) -> f32 {
  let g = grain(max(abs(x), abs(x + change))); let r = round(change / g) * g;
  return select(r, r - (x / g - round(x / g)) * g, abs(r) <= abs(x));
}
fn transmitted(z: f32) -> f32 { let d = min(z, PENDEPTH); return JER_R1 * exp(-d / JER_Z1) + JER_R2 * exp(-d / JER_Z2) - JER_FLOOR; }
// Returns the amount moved and what the destination's new spacing cannot hold; Metal's fast math folds the (b + a) − b form to a.
fn transfer(src: i32, dst: i32, amount: f32) -> vec2<f32> {
  let a = IN[src]; let b = IN[dst];
  var q = max(max(grain(a), grain(b)), max(grain(a - amount), grain(b + amount)));
  var moved = round(amount / q) * q;
  let g = max(grain(b), grain(b + moved));
  if (g > q) { q = g; moved = round(amount / q) * q; }
  let kept = round(b / g) * g;
  IN[src] = a - moved; IN[dst] = kept + moved;
  return vec2<f32>(moved, (b / g - round(b / g)) * g);
}
var<private> LEFT: vec3<f32>;
fn moveLayer(i: i32, srcK: i32, dstK: i32, amount: f32) {
  let ha = hOff(srcK) + i; let hb = hOff(dstK) + i;
  let qa = qOff(srcK) + i; let qb = qOff(dstK) + i;
  let wa = wOff(srcK) + i; let wb = wOff(dstK) + i;
  let h = IN[ha];
  let mh = transfer(ha, hb, amount); let f = mh.x / h;
  let mq = transfer(qa, qb, IN[qa] * f); let mw = transfer(wa, wb, IN[wa] * f);
  LEFT += vec3<f32>(mh.y, mq.y, mw.y);
}
fn placeLeft(at: i32, change: f32) -> f32 {
  let x = IN[at]; let g = grain(x);
  if (change == 0.0 || grain(x + change) > g) { return change; }
  let r = round(change / g) * g;
  IN[at] = x + r;
  return change - r;
}
${upwellingFunctions}${hy ? hybridHead(o) : ''}`;
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
    oSurfaceDensity: `${K}  let i = ${idx}; if (i >= C) { return; }
  let h0 = max(EPSO, IN[hOff(0) + i]);
  OD[O_RHOML + i] = eos(IN[qOff(0) + i] / h0 + TREF, IN[wOff(0) + i] / h0 + SREF);
}`,
    oGradRho: `${K}  let e = ${idx}; if (e >= E) { return; }
  OD[O_GRADRHO + e] = (OD[O_RHOML + MI[COE + 2 * e + 1]] - OD[O_RHOML + MI[COE + 2 * e]]) / MF[F_DC + e];
}`,
    oEdgeThickness: `${K}  let e = ${idx}; if (e >= E) { return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  OD[O_HEDGE + e] = 0.5 * (IN[hOff(0) + a] + IN[hOff(0) + b]);
  for (var k = 1; k < L; k++) { OD[O_HEDGE + k * E + e] = min(IN[hOff(k) + a], IN[hOff(k) + b]); }
  let sill = max(EPSO, min(OD[O_BATH + a], OD[O_BATH + b]) + 0.5 * (OD[O_ETA + a] + OD[O_ETA + b]));
  var sum = 0.0;
  for (var k = 0; k < L; k++) { sum += OD[O_HEDGE + k * E + e]; }
  if (sum > sill) { let f = sill / sum; for (var k = 0; k < L; k++) { OD[O_HEDGE + k * E + e] *= f; } }
}`,
    oFlux: hy ? `${K}  let e = ${idx}; if (e >= E) { return; }
  let wet = OD[O_EMASK + e] > 0.5;
  for (var k = 0; k < L; k++) {
    let n = k * E + e; let ue = IN[uOff(k) + e];
    let donor = MI[COE + 2 * e + select(1, 0, ue > 0.0)];
    OD[O_FLUX + n] = select(0.0, min(OD[O_HEDGE + n], max(0.0, IN[hOff(k) + donor])) * ue, wet);
  }
}` : `${K}  let e = ${idx}; if (e >= E) { return; }
  let wet = OD[O_EMASK + e] > 0.5;
  for (var k = 0; k < L; k++) {
    let n = k * E + e;
    var he = OD[O_HEDGE + n];
    if (k == 0) {
      let donor = MI[COE + 2 * e + select(1, 0, IN[uOff(0) + e] > 0.0)];
      he = min(he, max(0.0, IN[hOff(0) + donor]));
    }
    OD[O_FLUX + n] = select(0.0, he * IN[uOff(k) + e], wet);
  }
}`,
    oCellTendency: `${K}  let i = ${idx}; if (i >= C) { return; }
  var edges: array<i32, MAXE>; var nbrs: array<i32, MAXE>;
  for (var m = 0; m < MAXE; m++) { edges[m] = MI[EOC + MAXE * i + m]; nbrs[m] = MI[COC + MAXE * i + m]; }
  for (var k = 0; k < L; k++) {
    let hv = IN[hOff(k) + i]; let Ti = select(LABEL_T[k] - TREF, IN[qOff(k) + i] / hv, hv > 1e-6); let Si = select(LABEL_S[k] - SREF, IN[wOff(k) + i] / hv, hv > 1e-6);
    var divH = 0.0; var divQ = 0.0; var divW = 0.0; var lapQ = 0.0; var lapW = 0.0;
    for (var m = 0; m < MAXE; m++) {
      let e = edges[m]; let j = nbrs[m];
      let f = f32(MI[ESC + MAXE * i + m]) * OD[O_FLUX + k * E + e] * MF[F_DV + e];
      divH += f;
      let hvj = IN[hOff(k) + j]; let Tj = select(LABEL_T[k] - TREF, IN[qOff(k) + j] / hvj, hvj > 1e-6); let Sj = select(LABEL_S[k] - SREF, IN[wOff(k) + j] / hvj, hvj > 1e-6);
      divQ += f * select(Tj, Ti, f > 0.0); divW += f * select(Sj, Si, f > 0.0);
      if (k == 0 && OD[O_EMASK + e] > 0.5) { let wd = select(1.0, 2.0 * hv * hvj / max(hv + hvj, 1e-6), MIXDIFF > 0.0); lapQ += wd * MF[F_DV + e] * (Tj - Ti) / MF[F_DC + e]; lapW += wd * MF[F_DV + e] * (Sj - Si) / MF[F_DC + e]; }
    }
    let area = MF[F_AREA + i];
    var dh = -divH / area; var dQ = -divQ / area; var dW = -divW / area;
    if (k == 0) { let coefficient = select(DIFFUSION, MIXDIFF, MIXDIFF > 0.0); if (coefficient > 0.0) { dQ += coefficient * lapQ / area; dW += coefficient * lapW / area; } }
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
    oEdgePV: `${K}  let e = ${idx}; if (e >= E) { return; }
  let v1 = MI[VOE + 2 * e]; let v2 = MI[VOE + 2 * e + 1]; let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  for (var k = 0; k < L; k++) {
    let n = k * E + e;
    var hq = OD[O_HEDGE + n];
    if (${hy ? 'true' : 'k > 0'}) { hq = max(hq, CENTRING * 0.5 * (IN[hOff(k) + a] + IN[hOff(k) + b])); }
    OD[O_QE + n] = 0.5 * (OD[O_AVORT + k * V + v1] + OD[O_AVORT + k * V + v2]) / max(hq, ${hy ? 'PVFLOORK[k]' : 'select(PVFLOOR, THINO, k == 0)'});
  }
}`,
    oKineticPhi: `${K}  let i = ${idx}; if (i >= C) { return; }
  var edges: array<i32, MAXE>; var signs: array<f32, MAXE>; var dcs: array<f32, MAXE>; var dvs: array<f32, MAXE>;
  for (var m = 0; m < MAXE; m++) { let e = MI[EOC + MAXE * i + m]; edges[m] = e; signs[m] = f32(MI[ESC + MAXE * i + m]); dcs[m] = MF[F_DC + e]; dvs[m] = MF[F_DV + e]; }
  let area = MF[F_AREA + i];
  for (var k = 0; k < L; k++) {
    var kinetic = 0.0;
    for (var m = 0; m < MAXE; m++) { let u = IN[uOff(k) + edges[m]]; kinetic += abs(signs[m]) * 0.25 * dcs[m] * dvs[m] * u * u; }
    kinetic = kinetic / area;
    var phi = kinetic + OGRAV * OD[O_ETA + i];
    if (k > 0) {
      var p = (OD[O_RHOML + i] - RHO[k]) * IN[hOff(0) + i];
      for (var j = 1; j < k; j++) { p += (RHO[j] - RHO[k]) * IN[hOff(j) + i]; }
      phi += OGRAV * p / RHO0;
    }
    OD[O_PHI + k * C + i] = phi;
  }
}`,
    /*
     * The ∇⁴ closure's input, closureVelocity of ocean/layered.module.js:
     * under closureTokens 'beside' in LAPA, which the closure's first
     * Laplacian then reads; under 'interior' the first ring in LAPB with
     * FLUX, free once the cell tendency has read it, marking the thick and
     * fitted edges (1 thick, 2 the first ring), and oClosureRing the second
     * into LAPA. DEEPEST is the deepest class holding more than THIN metres
     * in both cells of an edge. After the second Laplacian, closureAdjoint
     * as two gathers: oClosureBack1 carries the second ring's ∇⁴ to the
     * first in LAPA, free by then, and oClosureBack2 the first ring's to
     * the thick edges in LAPB.
     */
    oDeepestEdge: `${K}  let e = ${idx}; if (e >= E) { return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  var ka = 0; var kb = 0;
  for (var k = L - 1; k > 0; k--) { if (ka == 0 && IN[hOff(k) + a] > THINO) { ka = k; } if (kb == 0 && IN[hOff(k) + b] > THINO) { kb = k; } }
  OD[O_DEEPEST + e] = f32(min(ka, kb));
}`,
    oClosureFill: `${K}  let n = ${idx}; if (n >= L * E) { return; }
  let k = n / E; let e = n % E;
  let own = IN[uOff(k) + e];
  let dst = select(O_LAPA, O_LAPB, CLOSURERINGS > 0.5);
  let thick = OD[O_EMASK + e] > 0.5 && OD[O_HEDGE + n] >= THINO;
  OD[dst + n] = own;
  if (CLOSURERINGS > 0.5) { OD[O_FLUX + n] = select(0.0, 1.0, thick); }
  if (k == 0 || OD[O_EMASK + e] < 0.5 || thick) { return; }
  if (CLOSURERINGS > 0.5 && !(OD[O_DEEPEST + e] > f32(k) + 0.5)) { return; }
  let nv = vec3<f32>(MF[F_NEDGE + 3 * e], MF[F_NEDGE + 3 * e + 1], MF[F_NEDGE + 3 * e + 2]);
  let tv = vec3<f32>(OD[O_TEDGE + 3 * e], OD[O_TEDGE + 3 * e + 1], OD[O_TEDGE + 3 * e + 2]);
  var fit = ClosureFit(0.0, 0.0, 0.0, 0.0, 0.0, false);
  for (var s = 0; s < MI[NEE + e]; s++) {
    let o = MI[EOE + MAXEE * e + s];
    if (OD[O_EMASK + o] < 0.5 || OD[O_HEDGE + k * E + o] < THINO) { continue; }
    fit = closureAdd(fit, o, nv, tv, IN[uOff(k) + o]);
  }
  if (fit.found) {
    OD[dst + n] = CLOSUREFILL * closureSolve(fit) + (1.0 - CLOSUREFILL) * own;
    if (CLOSURERINGS > 0.5) { OD[O_FLUX + n] = 2.0; }
  }
}`,
    oClosureRing: `${K}  let n = ${idx}; if (n >= L * E) { return; }
  let k = n / E; let e = n % E;
  OD[O_LAPA + n] = OD[O_LAPB + n];
  if (k == 0 || OD[O_EMASK + e] < 0.5 || OD[O_FLUX + n] > 0.5 || !(OD[O_DEEPEST + e] > f32(k) + 0.5)) { return; }
  let nv = vec3<f32>(MF[F_NEDGE + 3 * e], MF[F_NEDGE + 3 * e + 1], MF[F_NEDGE + 3 * e + 2]);
  let tv = vec3<f32>(OD[O_TEDGE + 3 * e], OD[O_TEDGE + 3 * e + 1], OD[O_TEDGE + 3 * e + 2]);
  var fit = ClosureFit(0.0, 0.0, 0.0, 0.0, 0.0, false);
  for (var s = 0; s < MI[NEE + e]; s++) {
    let o = MI[EOE + MAXEE * e + s];
    if (OD[O_EMASK + o] < 0.5 || OD[O_FLUX + k * E + o] < 0.5) { continue; }
    fit = closureAdd(fit, o, nv, tv, OD[O_LAPB + k * E + o]);
  }
  if (fit.found) { OD[O_LAPA + n] = CLOSUREFILL * closureSolve(fit) + (1.0 - CLOSUREFILL) * IN[uOff(k) + e]; }
}`,
    oClosureBack1: `${K}  let n = ${idx}; if (n >= L * E) { return; }
  let k = n / E; let o = n % E;
  if (k == 0 || abs(OD[O_FLUX + n] - 2.0) > 0.5) { return; }
  var sum = 0.0;
  for (var s = 0; s < MI[NEE + o]; s++) {
    let e = MI[EOE + MAXEE * o + s];
    if (OD[O_EMASK + e] < 0.5 || OD[O_FLUX + k * E + e] > 0.5 || !(OD[O_DEEPEST + e] > f32(k) + 0.5)) { continue; }
    sum += closureWeight(k, e, o, 2.5) * OD[O_LAPB + k * E + e] * MF[F_DC + e] * MF[F_DV + e];
  }
  OD[O_LAPA + n] = CLOSUREFILL * sum / (MF[F_DC + o] * MF[F_DV + o]);
}`,
    oClosureBack2: `${K}  let n = ${idx}; if (n >= L * E) { return; }
  let k = n / E; let o = n % E;
  if (k == 0 || abs(OD[O_FLUX + n] - 1.0) > 0.5) { return; }
  var sum = 0.0;
  for (var s = 0; s < MI[NEE + o]; s++) {
    let e = MI[EOE + MAXEE * o + s];
    if (abs(OD[O_FLUX + k * E + e] - 2.0) > 0.5) { continue; }
    sum += closureWeight(k, e, o, 1.5) * (OD[O_LAPB + k * E + e] + OD[O_LAPA + k * E + e]) * MF[F_DC + e] * MF[F_DV + e];
  }
  OD[O_LAPB + n] += CLOSUREFILL * sum / (MF[F_DC + o] * MF[F_DV + o]);
}`,
    oDivCurl: `${K}  let n = ${idx};
  let fromLap = P[1] > 0.5 || CLOSUREFILL > 0.0;
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
     * interfaceRate is ocean/layered.module.js's; P[3] is the step.
     */
    oMomentum: `fn tangentialU(k: i32, e: i32) -> f32 {
  var sum = 0.0;
  for (var s = 0; s < MAXEE; s++) { let slot = MAXEE * e + s; sum += MF[F_PVW + slot] * IN[uOff(k) + MI[EOE + slot]]; }
  return sum / MF[F_DC + e];
}
// Folding the k == 0 branch into the general one changes the drag's float32 rounding under Metal's fast math.
fn dragThickness(k: i32, e: i32) -> f32 { return select(max(OD[O_HEDGE + k * E + e], THINO), max(OD[O_HEDGE + e], THINO), k == 0); }
${rateFunction}fn firstThickBelow(k: i32, lo: u32, hi: u32) -> i32 {
  let start = k + 1;
  if (start < 32) {
    let m = lo & (0xffffffffu << u32(start));
    if (m != 0u) { return i32(countTrailingZeros(m)); }
    if (hi != 0u) { return 32 + i32(countTrailingZeros(hi)); }
    return -1;
  }
  if (start >= 64) { return -1; }
  let m = hi & (0xffffffffu << u32(start - 32));
  if (m != 0u) { return 32 + i32(countTrailingZeros(m)); }
  return -1;
}
${K}  let e = ${idx}; if (e >= E) { return; }
  if (OD[O_EMASK + e] < 0.5) { for (var k = 0; k < L; k++) { OUT[uOff(k) + e] = 0.0; } return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  let dc = MF[F_DC + e];
  var others: array<i32, MAXEE>; var weights: array<f32, MAXEE>; var live: array<bool, MAXEE>;
  for (var s = 0; s < MAXEE; s++) { let slot = MAXEE * e + s; others[s] = MI[EOE + slot]; weights[s] = MF[F_PVW + slot]; live[s] = OD[O_EMASK + others[s]] > 0.5; }
  var thickLo = 0u; var thickHi = 0u;
  for (var j = 0; j < L; j++) { if (OD[O_HEDGE + j * E + e] >= THINO) { if (j < 32) { thickLo |= 1u << u32(j); } else { thickHi |= 1u << u32(j - 32); } } }
  var above = 0; var uAbove = 0.0;
${blockTop}  for (var k = 0; k < L; k++) {
    let n = k * E + e;
    let u = IN[uOff(k) + e];
    let qHere = 0.5 * OD[O_QE + n];
    var pv = 0.0;
    for (var s = 0; s < MAXEE; s++) { let other = others[s]; let fo = select(0.0, OD[O_HEDGE + k * E + other] * IN[uOff(k) + other], live[s]); pv += weights[s] * fo * (qHere + 0.5 * OD[O_QE + k * E + other]); }
    let gradPhi = (OD[O_PHI + k * C + b] - OD[O_PHI + k * C + a]) / dc;
    var du = pv / dc - gradPhi;
    if (k == 0) { du -= OGRAV / RHO0 * 0.5 * OD[O_HEDGE + e] * OD[O_GRADRHO + e]; }
    let hk = OD[O_HEDGE + n];
    let he = max(hk, MOMDEPTH);
    var force = 0.0; var drag = 0.0;
${o.ekmanDepth > 0 ? '    if (k == 0 && hk >= OD[O_HM + e]) { force += OD[O_STRESS + e] / RHO0 * he / hk; }' : '    if (k == 0 && hk >= MOMDEPTH) { force += OD[O_STRESS + e] / RHO0; }'}
${dragLines}    if (below < 0) { force -= RBOT * abs(u) * u; }
    du += force / he + drag / dragThickness(k, e);
    if (NU4O > 0.0) { du -= NU4O * OD[O_LAPB + n]; }
    if (k > 0 && hk < THINO) { du = (uAbove - u) * P[2]; }
    OUT[uOff(k) + e] = du;
    if (k >= 1 && hk >= THINO) { above = k; }
    uAbove = u;
${topStep}  }
${o.ekmanDepth > 0 ? '  let hm = OD[O_HM + e];\n' : ''}  if (OD[O_HEDGE + e] < ${momentumReach}) {
    var top = 0.0; var depth = 0.0; var last = 0; var fraction = 1.0;
    for (var j = 0; j < L; j++) {
      let hj = OD[O_HEDGE + j * E + e]; let part = min(hj, max(0.0, ${momentumReach} - top));
      if (j == 0 || hj >= ${stressFloor}) { depth += part; }
      top += hj; last = j; fraction = select(1.0, part / hj, hj > 0.0);
      if (top >= ${momentumReach}) { break; }
    }
    if (depth > 0.0) {
      let push = OD[O_STRESS + e] / RHO0 / depth;
      for (var k = 0; k <= last; k++) { if (k == 0 || OD[O_HEDGE + k * E + e] >= ${stressFloor}) { OUT[uOff(k) + e] += push * select(fraction, 1.0, k < last); } }
    }
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
  let sill = min(OD[O_BATH + a], OD[O_BATH + b]) + 0.5 * (OD[O_ETA + a] + OD[O_ETA + b]);
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
${hy ? '      lower *= min(modeWeight(k + 1, a), modeWeight(k + 1, b));\n' : ''}    }
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
    /*
     * The Ekman block of ocean/layered.module.js, one edge a thread, in the
     * eddy transport's scratch: each class's change in the edge's volume
     * flux, in whole multiples of the coarser cell's spacing and with the
     * thinnest class thick in both cells taking back what the rounding left
     * of the edge's net, and the heat and salt it carries from its upwind
     * cell. The stress shares and the block's parts wait in EDDYQ and EDDYW
     * until the changes are made; UEK keeps the edge's Ekman transport.
     * oEkmanApply and oEkmanLeftover are the eddy pass's, mixed layer
     * included.
     */
    ...(o.ekmanBlock ? {
      oEkmanFlux: `fn ekmanColumn(k: i32, a: i32, b: i32, squeeze: f32) -> f32 {
  return select(min(IN[hOff(k) + a], IN[hOff(k) + b]), 0.5 * (IN[hOff(0) + a] + IN[hOff(0) + b]), k == 0) * squeeze;
}
${eddyShared}${K}  let e = ${idx}; if (e >= E) { return; }
  OD[O_UEK + e] = 0.0;
  for (var k = 0; k < L; k++) { OD[O_EDDYV + k * E + e] = 0.0; OD[O_EDDYQ + k * E + e] = 0.0; OD[O_EDDYW + k * E + e] = 0.0; }
  if (OD[O_EMASK + e] < 0.5) { return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  let sill = max(EPSO, min(OD[O_BATH + a], OD[O_BATH + b]) + 0.5 * (OD[O_ETA + a] + OD[O_ETA + b]));
  var sum = 0.0;
  for (var k = 0; k < L; k++) { sum += ekmanColumn(k, a, b, 1.0); }
  let squeeze = select(1.0, sill / sum, sum > sill);
  let wind = ${o.ekmanDepth > 0 ? 'OD[O_HM + e]' : 'MOMDEPTH'};
  let reach = max(wind, ekmanColumn(0, a, b, squeeze));
  let rate = RINT / reach; let f = OD[O_FEDGE + e];
  let denominator = RHO0 * (rate * rate + f * f);
  if (!(denominator > 0.0)) { return; }
  var along = 0.0;
  for (var s = 0; s < MAXEE; s++) { let slot = MAXEE * e + s; along += MF[F_PVW + slot] * OD[O_STRESS + MI[EOE + slot]]; }
  let transport = (rate * OD[O_STRESS + e] + f * (along / MF[F_DC + e])) / denominator;
  OD[O_UEK + e] = transport;
  if (transport == 0.0) { return; }
  let donor = select(b, a, transport > 0.0);
  let blockDepth = max(wind, IN[hOff(0) + donor]);
  var top = 0.0; var depth = 0.0; var reaching = true; var z = 0.0; var total = 0.0;
  for (var k = 0; k < L; k++) {
    var share = 0.0;
    if (reaching) {
      let hk = ekmanColumn(k, a, b, squeeze);
      let part = min(hk, max(0.0, reach - top));
      if (k == 0 || hk >= ${stressFloor}) { share = part; depth += part; }
      top += hk; reaching = top < reach;
    }
    let hd = IN[hOff(k) + donor];
    let inBlock = min(select(max(0.0, hd - EPSO), hd, k == 0), max(0.0, blockDepth - z));
    z += hd; total += inBlock;
    OD[O_EDDYQ + k * E + e] = share; OD[O_EDDYW + k * E + e] = inBlock;
  }
  if (!(depth > 0.0) || !(total > 0.0)) {
    for (var k = 0; k < L; k++) { OD[O_EDDYQ + k * E + e] = 0.0; OD[O_EDDYW + k * E + e] = 0.0; }
    return;
  }
  let volumeFlux = transport * MF[F_DV + e];
  var scale = 1.0;
  for (var k = 0; k < L; k++) {
    let change = volumeFlux * (OD[O_EDDYW + k * E + e] / total - OD[O_EDDYQ + k * E + e] / depth);
    OD[O_EDDYV + k * E + e] = change;
    if (change != 0.0) {
      let allowed = eddyAllowance(k, select(b, a, change > 0.0));
      if (abs(change) > allowed) { scale = min(scale, allowed / abs(change)); }
    }
  }
  let sa = MF[F_AREA + a]; let sb = MF[F_AREA + b];
  var net = 0.0; var taker = -1; var finest = 0.0;
  for (var k = 0; k < L; k++) {
    let ha = IN[hOff(k) + a]; let hb = IN[hOff(k) + b];
    let change = scale * OD[O_EDDYV + k * E + e];
    var moved = 0.0;
    if (change != 0.0) { moved = eddyShare(ha, hb, sa, sb, change * P[0]); net += moved; }
    OD[O_EDDYV + k * E + e] = moved;
    if (min(ha, hb) > THINO) { let g = grain(max(ha, hb)); if (taker < 0 || g < finest) { taker = k; finest = g; } }
  }
  if (taker >= 0) { OD[O_EDDYV + taker * E + e] -= net; }
  for (var k = 0; k < L; k++) {
    let moved = OD[O_EDDYV + k * E + e];
    if (moved == 0.0) { OD[O_EDDYQ + k * E + e] = 0.0; OD[O_EDDYW + k * E + e] = 0.0; continue; }
    let qa = IN[qOff(k) + a]; let qb = IN[qOff(k) + b]; let wa = IN[wOff(k) + a]; let wb = IN[wOff(k) + b];
    let toB = moved > 0.0; let carried = moved / IN[hOff(k) + select(b, a, toB)];
    OD[O_EDDYV + k * E + e] = moved / P[0];
    OD[O_EDDYQ + k * E + e] = eddyShare(qa, qb, sa, sb, carried * select(qb, qa, toB)) / P[0];
    OD[O_EDDYW + k * E + e] = eddyShare(wa, wb, sa, sb, carried * select(wb, wa, toB)) / P[0];
  }
}`,
      oEkmanApply: eddyApply(0),
      oEkmanLeftover: eddyLeftover(0),
    } : {}),
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
    oRescale: `${K}  let i = ${idx}; if (i >= C) { return; }
  if (OD[O_CMASK + i] < 0.5) { return; }
  var sum = 0.0; var dh = 0.0; var dQ = 0.0; var dW = 0.0;
  for (var k = 1; k < L; k++) {
    var hv = IN[hOff(k) + i];
    if (hv < EPSO) {
      let held = hv > 1e-9;
      let t = LABEL_T[k]; let s = LABEL_S[k];
      let tHeld = select(t, clamp(IN[qOff(k) + i] / hv + TREF, t - 30.0, t + 30.0), held); let sHeld = select(s, clamp(IN[wOff(k) + i] / hv + SREF, s - 5.0, s + 5.0), held);
      dh += EPSO - hv; dQ += EPSO * (t - TREF) - hv * (tHeld - TREF); dW += EPSO * (s - SREF) - hv * (sHeld - SREF);
      hv = EPSO;
      IN[hOff(k) + i] = EPSO; IN[qOff(k) + i] = EPSO * (t - TREF); IN[wOff(k) + i] = EPSO * (s - SREF);
    }
    sum += hv;
  }
  IN[hOff(0) + i] -= dh; IN[qOff(0) + i] -= dQ; IN[wOff(0) + i] -= dW;
  sum += IN[hOff(0) + i];
  let excess = (OD[B_BAVG + i] - (sum - OD[O_BATH + i])) / sum;
  for (var k = 0; k < L; k++) { IN[hOff(k) + i] += IN[hOff(k) + i] * excess; IN[qOff(k) + i] += IN[qOff(k) + i] * excess; IN[wOff(k) + i] += IN[wOff(k) + i] * excess; }
  OD[O_ETA + i] = OD[B_BAVG + i];
}`,
    oVelocityShiftClamp: `${K}  let e = ${idx}; if (e >= E) { return; }
  if (OD[O_EMASK + e] < 0.5) { return; }
  for (var k = 1; k < L; k++) { if (OD[O_HEDGE + k * E + e] < ${hy ? 'THINEDGE' : 'THINO'}) { IN[uOff(k) + e] = IN[uOff(k - 1) + e]; } }
  var sumH = 0.0; var transport = 0.0;
  for (var k = 0; k < L; k++) { let he = OD[O_HEDGE + k * E + e]; sumH += he; transport += he * IN[uOff(k) + e]; }
  let shift = (OD[B_BAVG + C + e] - transport) / max(sumH, EPSO);
  for (var k = 0; k < L; k++) { let n = uOff(k) + e; IN[n] = clamp(IN[n] + shift, -SPEEDLIM, SPEEDLIM); }
}`,
    /*
     * The mixed layer's momentum through its exchanges, exchangeMomentum
     * of ocean/layered.module.js, one edge a thread: HEDGE holds the edge
     * thicknesses the exchanges started from and is left holding the new.
     */
    oExchangeMomentum: `${K}  let e = ${idx}; if (e >= E) { return; }
  if (OD[O_EMASK + e] < 0.5) { return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  let sill = max(EPSO, min(OD[O_BATH + a], OD[O_BATH + b]) + 0.5 * (OD[O_ETA + a] + OD[O_ETA + b]));
  var sum = 0.5 * (IN[hOff(0) + a] + IN[hOff(0) + b]); var carried = OD[O_HEDGE + e] * IN[uOff(0) + e];
  for (var k = 1; k < L; k++) { sum += min(IN[hOff(k) + a], IN[hOff(k) + b]); carried += OD[O_HEDGE + k * E + e] * IN[uOff(k) + e]; }
  let scale = select(1.0, sill / sum, sum > sill);
  let reach = max(OD[O_HEDGE + e], 0.5 * (IN[hOff(0) + a] + IN[hOff(0) + b]) * scale);
  var top = 0.0; var within = 0.0; var mixed = 0.0;
  for (var k = 0; k < L; k++) {
    if (top >= reach) { break; }
    let old = OD[O_HEDGE + k * E + e]; let part = min(old, reach - top);
    mixed += part * IN[uOff(k) + e]; within += part; top += old;
  }
  let uMixed = select(IN[uOff(0) + e], mixed / within, within > 0.0);
  top = 0.0;
  var total = 0.0; var held = 0.0; var above = 0.0;
  for (var k = 0; k < L; k++) {
    let n = uOff(k) + e;
    let hk = select(min(IN[hOff(k) + a], IN[hOff(k) + b]), 0.5 * (IN[hOff(0) + a] + IN[hOff(0) + b]), k == 0) * scale;
    let inside = max(0.0, min(hk, reach - top));
    var v = IN[n];
    if (k > 0 && hk < THINO) { v = above; } else if (hk > 0.0) { v = (inside * uMixed + (hk - inside) * v) / hk; }
    IN[n] = v; OD[O_HEDGE + k * E + e] = hk;
    total += hk; held += hk * v; top += hk; above = v;
  }
  let shift = (carried - held) / max(total, EPSO);
  for (var k = 0; k < L; k++) { let n = uOff(k) + e; IN[n] = clamp(IN[n] + shift, -SPEEDLIM, SPEEDLIM); }
}`,
    oReadSurface: `${K}  let i = ${idx}; if (i >= C) { return; }
  let iced = OD[O_SURFICE + i] > 0.0;
  OD[O_ICED + i] = select(0.0, 1.0, iced);
  let t0 = select(OD[O_SURFT + i], FREEZE, iced);
  OD[O_SURFACEIN + i] = t0;
  // Open water adds the surface's change since PREVT0: reading T0 = Q/h back whole gains heat under fast-math division.
  if (iced) { IN[qOff(0) + i] = IN[hOff(0) + i] * (t0 - TREF); } else { IN[qOff(0) + i] += IN[hOff(0) + i] * (t0 - OD[O_PREVT0 + i]); }
  let light = PH[PH_PEN + i];
${o.retreatDaily ? '  OD[O_SUNSTEP + i] = select(0.0, light / max(PH[PH_TRANS + i], 1.0e-6), PH[PH_TRANS + i] > 1.0e-6);\n' : ''}  if (light == 0.0) { return; }
  let taken = PH[PH_TRANS + i];
  var given = 0.0; var last = 0;
  if (taken > 0.0) {
    var z = IN[hOff(0) + i]; var top = transmitted(z); var carry = 0.0;
    for (var k = 1; k < L; k++) {
      let hk = IN[hOff(k) + i];
      z += hk;
      let bottom = transmitted(z);
      let want = carry + light * (top - bottom) / taken;
      if (hk > ${hy ? '0.0' : 'THINO'}) {
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
    oStressFromAtmosphere: `${K}  let e = ${idx}; if (e >= E) { return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  if (OD[O_EMASK + e] < 0.5) { OD[O_STRESS + e] = 0.0; return; }
  let coverA = select(0.0, select(1.0, PH[PH_CONC + a], PH[PH_CONC + a] > 0.0), OD[O_ICED + a] > 0.5);
  let coverB = select(0.0, select(1.0, PH[PH_CONC + b], PH[PH_CONC + b] > 0.0), OD[O_ICED + b] > 0.5);
  let through = 1.0 - 0.5 * (coverA + coverB) * (1.0 - TRANSMIT);
${o.implicitStress ? '  if (PH[PH_STRESSOK] > 0.5) { OD[O_STRESS + e] = through * PH[PH_STRESS + e]; return; }\n' : ''}  let bottom = (K - 1) * C;
  let rhoA = S[S_PI + a] * LV[L_SM + K - 1] / (RGAS * S[S_TH + bottom + a] * D[D_EXM + bottom + a]);
  let rhoB = S[S_PI + b] * LV[L_SM + K - 1] / (RGAS * S[S_TH + bottom + b] * D[D_EXM + bottom + b]);
  let fa = PH[PH_DRAG + a] * rhoA * max(D[D_WIND + a], GUSTO); let fb = PH[PH_DRAG + b] * rhoB * max(D[D_WIND + b], GUSTO);
  OD[O_STRESS + e] = through * 0.5 * (fa + fb) * S[S_U + (K - 1) * E + e];
}`,
    ...(o.ekmanDepth > 0 ? { oWindDepth: `${K}  let e = ${idx}; if (e >= E) { return; }
  var hm = MOMDEPTH;
  let f = abs(OD[O_FEDGE + e]); let normal = OD[O_STRESS + e];
  if (OD[O_EMASK + e] > 0.5 && EKMAN * sqrt(abs(normal) / RHO0) < MOMDEPTH * f) {
    var along = 0.0;
    for (var s = 0; s < MAXEE; s++) { let slot = MAXEE * e + s; along += MF[F_PVW + slot] * OD[O_STRESS + MI[EOE + slot]]; }
    along = along / MF[F_DC + e];
    hm = min(MOMDEPTH, max(MOMFLOOR, EKMAN * sqrt(sqrt(normal * normal + along * along) / RHO0) / f));
  }
  OD[O_HM + e] = hm;
}` } : {}),
    ...(o.shearEntrainment > 0 ? { oShearEdge: `${K}  let e = ${idx}; if (e >= E) { return; }
  var du = 0.0;
  if (OD[O_EMASK + e] > 0.5) {
    for (var k = 1; k < L; k++) { if (OD[O_HEDGE + k * E + e] >= THINO) { du = IN[uOff(0) + e] - IN[uOff(k) + e]; break; } }
  }
  OD[O_DU + e] = du;
}` } : {}),
    oMixedReach: `${K}  let i = ${idx}; if (i >= C) { return; }
  var sum = 0.0; var count = 0.0;
  for (var m = 0; m < MAXE; m++) {
    if (MI[ESC + MAXE * i + m] == 0) { continue; }
    let j = MI[COC + MAXE * i + m];
    if (OD[O_CMASK + j] > 0.5) { sum += IN[hOff(0) + j]; count += 1.0; }
  }
  OD[O_REACH + i] = select(1.0e30, NBRRATIO * sum / max(count, 1.0), NBRRATIO > 0.0 && count > 0.0);
}`,
    oMixedLayer: `${K}  let i = ${idx}; if (i >= C) { return; }
  if (OD[O_CMASK + i] < 0.5) { return; }
${o.upwellingSource > 0 ? '  let upwelling = upwellingBudget(i);\n  var upLeft = upwelling.x;\n' : ''}  var wv = vec3<f32>(0.0, 0.0, 0.0);
  for (var m = 0; m < MAXE; m++) {
    let e = MI[EOC + MAXE * i + m];
    let s = abs(f32(MI[ESC + MAXE * i + m])) * 0.5 * MF[F_DC + e] * MF[F_DV + e] * OD[O_STRESS + e];
    wv += s * vec3<f32>(MF[F_NEDGE + 3 * e], MF[F_NEDGE + 3 * e + 1], MF[F_NEDGE + 3 * e + 2]);
  }
  wv = wv / MF[F_AREA + i];
  let tau = length(wv);
  let ustar3 = pow(tau / RHO0, 1.5); let stir = STIRRING * exp(-IN[hOff(0) + i] / STIRDEPTH);
${o.langmuir > 0 ? '  let began = IN[hOff(0) + i];\n  let openSea = 1.0 - select(0.0, select(1.0, OD[O_SURFCONC + i], OD[O_SURFCONC + i] > 0.0), OD[O_SURFICE + i] > 0.0);\n' : ''}${o.retreatEkman > 0 ? '  let fCell = abs(TWOOMEGA * sin(MF[F_LAT + i]));\n' : ''}  let reach = OD[O_REACH + i];
  let deepest = min(MAXMIXED, reach);
  let s0 = IN[wOff(0) + i] / IN[hOff(0) + i] + SREF;
  let salted = s0 * OD[O_FRESH + i] / 1000.0 + (s0 - ICESAL) * (iceVolume(i) - OD[O_PREVVOL + i]) * ICEDENS / 1000.0;
  let buoyancy = OGRAV * (alphaT(OD[O_SURFACEIN + i], s0) * (OD[O_PREVT0 + i] - OD[O_SURFACEIN + i]) * PH[PH_CAP + i] / RHOCP + betaS(OD[O_SURFACEIN + i], s0) * salted) / P[6];
  let loss = OD[O_BUOY + i] + (buoyancy - OD[O_BUOY + i]) * select(1.0, min(1.0, P[6] / BUOYMEM), BUOYMEM > 0.0);
  OD[O_BUOY + i] = loss;
${o.retreatDaily ? dailySums : ''}  var rm = eos(IN[qOff(0) + i] / IN[hOff(0) + i] + TREF, IN[wOff(0) + i] / IN[hOff(0) + i] + SREF);
  var budget = CONVRATE * P[6];
  for (var k = 1; k < L; k++) {
    if (IN[hOff(0) + i] >= deepest || budget <= 0.0) { break; }
    let hk = IN[hOff(k) + i];
    if (hk <= EPSO) { continue; }
    var take = min(min(hk - EPSO, deepest - IN[hOff(0) + i]), budget);
    let eroding = EROSION > 0.5 && rm < DENSEST[k];
    if (eroding) {
      if (rm >= LIGHTEST[k] && loss > 0.0) { take = min(take, loss * P[6] * RHO0 * hk / (OGRAV * IN[hOff(0) + i] * (DENSEST[k] - rm))); } else { take = 0.0; }
    } else if (EROSION < 0.5 && RHO[k] > rm) { continue; }
    if (take > 0.0) {
      moveLayer(i, k, 0, take); budget -= take;
      rm = eos(IN[qOff(0) + i] / IN[hOff(0) + i] + TREF, IN[wOff(0) + i] / IN[hOff(0) + i] + SREF);
    }
    if (eroding) { break; }
  }
  var below = -1;
  for (var k = 1; k < L; k++) { if (IN[hOff(k) + i] > THINO) { below = k; break; } }
${o.shearEntrainment > 0 ? `  var hold = 0.0;
  if (below > 0) {
    var sv = vec3<f32>(0.0, 0.0, 0.0);
    for (var m = 0; m < MAXE; m++) {
      let e = MI[EOC + MAXE * i + m];
      let s = abs(f32(MI[ESC + MAXE * i + m])) * 0.5 * MF[F_DC + e] * MF[F_DV + e] * OD[O_DU + e];
      sv += s * vec3<f32>(MF[F_NEDGE + 3 * e], MF[F_NEDGE + 3 * e + 1], MF[F_NEDGE + 3 * e + 2]);
    }
    sv = sv / MF[F_AREA + i];
    let h0 = IN[hOff(0) + i];
    let db = max(OGRAV * (RHOA[below] - heldAnomaly(IN[qOff(0) + i] / h0, IN[wOff(0) + i] / h0)) / RHO0, OGRAV * DENSTOL / RHO0);
    let ri = db * h0 / (dot(sv, sv) + 1e-12);
    if (ri < RICRIT) {
      let reached = h0 * sqrt(RICRIT / ri);
${o.upwellingSource > 0 ? `      let want = min(reached - h0, deepest - h0);
      ${drawFromSource('want')}
      let take = min(want - drawn, IN[hOff(below) + i] - EPSO);
      if (take > 0.0) { moveLayer(i, below, 0, take); }
      if (drawn + take > 0.0) { rm = eos(IN[qOff(0) + i] / IN[hOff(0) + i] + TREF, IN[wOff(0) + i] / IN[hOff(0) + i] + SREF); }
` : `      let take = min(min(reached - h0, IN[hOff(below) + i] - EPSO), deepest - h0);
      if (take > 0.0) { moveLayer(i, below, 0, take); rm = eos(IN[qOff(0) + i] / IN[hOff(0) + i] + TREF, IN[wOff(0) + i] / IN[hOff(0) + i] + SREF); }
`}      hold = min(reached, IN[hOff(0) + i]);
      below = -1;
      for (var k = 1; k < L; k++) { if (IN[hOff(k) + i] > THINO) { below = k; break; } }
    }
  }
` : ''}  if (below > 0 && IN[hOff(0) + i] < deepest) {
    let db = max(1e-3, OGRAV * (RHO[below] - rm) / RHO0);
${o.upwellingSource > 0 ? `    let want = min(2.0 * stir * ustar3 / (IN[hOff(0) + i] * db) * P[6], reach - IN[hOff(0) + i]);
    ${drawFromSource('want')}
    let entrain = min(want - drawn, IN[hOff(below) + i] - EPSO);
    if (entrain > 0.0) { moveLayer(i, below, 0, entrain); }
    if (drawn + entrain > 0.0) { rm = eos(IN[qOff(0) + i] / IN[hOff(0) + i] + TREF, IN[wOff(0) + i] / IN[hOff(0) + i] + SREF); }
` : `    let entrain = min(min(2.0 * stir * ustar3 / (IN[hOff(0) + i] * db) * P[6], IN[hOff(below) + i] - EPSO), reach - IN[hOff(0) + i]);
    if (entrain > 0.0) { moveLayer(i, below, 0, entrain); rm = eos(IN[qOff(0) + i] / IN[hOff(0) + i] + TREF, IN[wOff(0) + i] / IN[hOff(0) + i] + SREF); }
`}  }
  below = -1;
  for (var k = 1; k < L; k++) { if (IN[hOff(k) + i] > THINO) { below = k; break; } }
  var excess = max(0.0, IN[hOff(0) + i] - MAXMIXED);
  if (NEUTRALSNAP > 0.5 && below > 0 && rm >= RHO[below] - DENSTOL) { excess = max(excess, IN[hOff(0) + i] - SHALLOWMIXED); }
  detrain(i, excess, rm);
  if (IN[hOff(0) + i] > reach) { detrain(i, (IN[hOff(0) + i] - max(reach, SHALLOWMIXED)) * min(1.0, P[6] / DETRAINT), rm); }
${o.retreatDaily ? dailyDrive : ''}  if (${o.retreatDaily ? 'dayDrive < -1e-9 && dayGate' : 'loss'} < -1e-9) {
    let monin = ${o.shearEntrainment > 0 ? `max(max(SHALLOWMIXED, ${retreatTerm}), hold)` : `max(SHALLOWMIXED, ${retreatTerm})`};
    if (IN[hOff(0) + i] > monin) { detrain(i, (IN[hOff(0) + i] - monin) * min(1.0, P[6] / DETRAINT), rm); }
  }
  if (IN[hOff(0) + i] < MINTHICK) {
    let refillTo = ${o.refillRate > 0 ? 'max(min(MINTHICK, IN[hOff(0) + i] + REFILLRATE * P[6]), min(MINTHICK, REFILLFLOOR))' : 'MINTHICK'};
${o.upwellingSource > 0 ? `    let wanted = refillTo - IN[hOff(0) + i];
    ${drawFromSource('wanted')}
` : ''}    for (var k = 1; k < L; k++) {
      if (IN[hOff(0) + i] >= refillTo) { break; }
      let available = IN[hOff(k) + i] - EPSO;
      if (available > 0.0) { moveLayer(i, k, 0, min(available, refillTo - IN[hOff(0) + i])); }
    }
  }
  let fraction = min(1.0, P[6] / RESTORET);
  for (var k = 1; k < L; k++) {
    let hk = IN[hOff(k) + i];
    if (hk <= THINO) { continue; }
    let r = eosAnomaly(IN[qOff(k) + i] / hk + TREF, IN[wOff(k) + i] / hk + SREF);
    let label = RHOA[k];
    if (abs(r - label) <= RESTTOL) { continue; }
    let stepK = select(-1, 1, r < label);
    var j = k + stepK;
    loop {
      if (j < 1 || j >= L) { break; }
      let hj = IN[hOff(j) + i];
      if (hj > THINO) {
        let rd = eosAnomaly(IN[qOff(j) + i] / hj + TREF, IN[wOff(j) + i] / hj + SREF);
        if (f32(stepK) * (rd - label) > RESTTOL) {
          let amount = min(fraction * hk * (label - r) / (rd - label), hj - EPSO);
          if (amount > 0.0) { moveLayer(i, j, k, amount); }
          break;
        }
      }
      j += stepK;
    }
  }
  let left = LEFT + vec3<f32>(OD[O_REM + i], OD[O_REM + C + i], OD[O_REM + 2 * C + i]);
  if (any(left != vec3<f32>(0.0))) {
    var kh = 0; var kq = 0; var kw = 0; var gh = 0.0; var gq = 0.0; var gw = 0.0;
    for (var k = 0; k < L; k++) {
      let hk = IN[hOff(k) + i];
      if (hk <= THINO) { continue; }
      let qk = IN[qOff(k) + i]; let wk = IN[wOff(k) + i];
      let a = grain(hk); let b = grain(qk); let c = grain(wk);
      if (grain(hk + left.x) <= a && (gh == 0.0 || a < gh)) { kh = k; gh = a; }
      if (grain(qk + left.y) <= b && (gq == 0.0 || b < gq)) { kq = k; gq = b; }
      if (grain(wk + left.z) <= c && (gw == 0.0 || c < gw)) { kw = k; gw = c; }
    }
    OD[O_REM + i] = select(left.x, placeLeft(hOff(kh) + i, left.x), gh > 0.0);
    OD[O_REM + C + i] = select(left.y, placeLeft(qOff(kq) + i, left.y), gq > 0.0);
    OD[O_REM + 2 * C + i] = select(left.z, placeLeft(wOff(kw) + i, left.z), gw > 0.0);
  }
}`,
    oSalt: `${K}  let i = ${idx}; if (i >= C) { return; }
  if (OD[O_CMASK + i] < 0.5) { return; }
  let s = IN[wOff(0) + i] / IN[hOff(0) + i] + SREF;
  var w = IN[wOff(0) + i] + s * OD[O_FRESH + i] / 1000.0;
  let volume = iceVolume(i);
  let grown = (volume - OD[O_PREVVOL + i]) * ICEDENS / 1000.0;
  w += (s - ICESAL) * grown;
  w = max(-SREF * IN[hOff(0) + i], w);
  IN[wOff(0) + i] = w;
  OD[O_FRESH + i] = 0.0;
  OD[O_PREVVOL + i] = volume;
}`,
    oWriteSurface: `${K}  let i = ${idx}; if (i >= C) { return; }
  if (OD[O_CMASK + i] < 0.5) { PH[PH_OFLUX + i] = 0.0; return; }
  PH[PH_CAP + i] = RHOCP * max(IN[hOff(0) + i], 1.0);
${o.penetration ? hy ? '  PH[PH_TRANS + i] = transmitted(IN[hOff(0) + i]);\n' : `  var receiver = false;
  for (var k = 1; k < L; k++) { if (IN[hOff(k) + i] > THINO) { receiver = true; break; } }
  PH[PH_TRANS + i] = select(0.0, transmitted(IN[hOff(0) + i]) - transmitted(OD[O_BATH + i] + OD[O_ETA + i]), receiver);
` : ''}  OD[O_S0 + i] = IN[wOff(0) + i] / IN[hOff(0) + i] + SREF;
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
  if (OD[O_CMASK + i] > 0.5) {
    var fresh = PH[PH_EVAP + i] * P[6] - delta;
    for (var p = i32(OD[O_DRAINSTART + i]); p < i32(OD[O_DRAINSTART + i + 1]); p++) {
      let j = i32(OD[O_DRAINCELL + p]);
      let runoff = PH[PH_RUNOFF + j]; let ran = OD[O_RUNOFFSEEN + j];
      fresh -= select(runoff, runoff - ran, runoff >= ran) * OD[O_DRAINW + p];
      OD[O_RUNOFFSEEN + j] = runoff;
    }
    OD[O_FRESH + i] += fresh;
  }
  OD[O_RAINSEEN + i] = rain;
}`,
    ...(hy ? hybridKernels(o, K, idx) : {}),
  };
}

/*
 * The global sums behind the ocean diagnostics: [name, how the
 * workgroup partials combine, the cell's term]. Edge terms are taken by
 * the edge's first cell.
 */
const OCEAN_REDUCED = [
  ['area', 'sum', 'wet * a'], ['depth', 'sum', 'wet * a * h0'], ['salinity', 'sum', 'wet * a * (IN[wOff(0) + i] / max(EPSO, h0) + SREF)'],
  ['thermocline', 'sum', 'wet * a * thermo'], ['heat', 'sum', 'wet * a * RHOCP * heat'],
  ['interiorT', 'sum', 'wet * a * interiorQ'], ['interiorH', 'sum', 'wet * a * interiorH'],
  ['ssh', 'max', 'wet * abs(OD[O_ETA + i])'], ['speed', 'max', 'speed'], ['limited', 'sum', 'limited'], ['transport', 'max', 'transport'],
  ['thinned', 'sum', `wet * a * select(0.0, 1.0, h0 < min(MINTHICK, h0 + interiorH - f32(L - 1) * EPSO) - ${THINNED_MARGIN})`],
];
const oceanReducedSetup = (thermoclineLayers) => `    let a = MF[F_AREA + i];
    let wet = select(0.0, 1.0, OD[O_CMASK + i] > 0.5);
    let h0 = IN[hOff(0) + i];
    var thermo = 0.0; var heat = 0.0; var interiorQ = 0.0; var interiorH = 0.0;
    for (var k = 0; k < L; k++) {
      heat += IN[qOff(k) + i] + TREF * IN[hOff(k) + i];
      if (k <= ${thermoclineLayers}) { thermo += IN[hOff(k) + i]; }
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
      let sill = max(EPSO, min(OD[O_BATH + ca], OD[O_BATH + cb]) + 0.5 * (OD[O_ETA + ca] + OD[O_ETA + cb]));
      var sum = 0.5 * (IN[hOff(0) + ca] + IN[hOff(0) + cb]); var flow = sum * u0;
      for (var k = 1; k < L; k++) { let hk = min(IN[hOff(k) + ca], IN[hOff(k) + cb]); sum += hk; flow += hk * IN[uOff(k) + e]; }
      transport = max(transport, abs(select(1.0, sill / sum, sum > sill) * flow) * MF[F_DV + e]);
    }`;

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
      let sill = max(EPSO, min(OD[O_BATH + ca], OD[O_BATH + cb]) + 0.5 * (OD[O_ETA + ca] + OD[O_ETA + cb]));
      var sum = 0.5 * (IN[hOff(0) + ca] + IN[hOff(0) + cb]); var flow = sum * u0;
      for (var k = 1; k < L; k++) { let hk = min(IN[hOff(k) + ca], IN[hOff(k) + cb]); sum += hk; flow += hk * IN[uOff(k) + e]; }
      transport = max(transport, abs(select(1.0, sill / sum, sum > sill) * flow) * MF[F_DV + e]);
    }`;

export function createLayeredOcean(core, { registers = null, ...options } = {}) {
  const o = { ...OCEAN_DEFAULTS, ...options };
  const hybrid = o.vertical === 'hybrid';
  if (o.vertical !== undefined && o.vertical !== 'slab' && !hybrid) throw new Error(`vertical is 'slab' or 'hybrid', not ${o.vertical}`);
  if (hybrid) {
    const refused = HYBRID_REFUSED.filter((name) => o[name]);
    if (o.upwellingSource > 0) refused.push('upwellingSource');
    if (o.closureFill > 0) refused.push('closureFill');
    if (refused.length) throw new Error(`${refused.join(', ')} belong to the slab ocean`);
    if (!GENERATORS.includes(o.generator)) throw new Error(`the generator is one of ${GENERATORS.join(', ')}, not ${o.generator}`);
  }
  const table = hybrid ? o.layerTable ?? HYBRID_LAYERS : null;
  if (o.closureTokens !== 'interior' && o.closureTokens !== 'beside') throw new Error(`closureTokens is 'interior' or 'beside', not ${o.closureTokens}`);
  if (o.penetration !== null && !JERLOV[o.penetration]) throw new Error(`penetration is one of ${Object.keys(JERLOV).join(', ')} or null, not ${o.penetration}`);
  if (o.shearEntrainment > 0 && !o.carryMomentum) throw new Error('shearEntrainment needs carryMomentum: without the entrained water\'s momentum mixed in, the base\'s Richardson number never rises and the mixed layer swallows a class a step');
  if (!(o.upwellingSource >= 0) || (o.upwellingSource > 0 && !o.ekmanBlock)) throw new Error(`upwellingSource is a share of the Ekman divergence, not negative, and needs ekmanBlock for the transport it divides; not ${o.upwellingSource} with ekmanBlock ${o.ekmanBlock}`);
  mixingOptionCheck(o);
  const day = { span: 0 };
  o.closureRings = o.closureTokens === 'interior' && o.closureFill > 0;
  const { device, buffers, mesh, meshSpacing } = core;
  const C = core.C, E = core.E, V = core.V;
  const L = hybrid ? table.K : o.densities.length + 1;
  const rho = hybrid ? Array.from(table.densities) : [o.density, ...o.densities];
  const labels = hybrid ? rho : o.densities;
  const labelS = hybrid ? rho.map((r) => { const j = LAYER_DENSITIES.indexOf(r); return j >= 0 ? LAYER_SALINITIES[j] : o.referenceS; }) : [o.referenceS, ...o.salinities];
  const labelT = rho.map((r, k) => Math.max(FREEZING_POINT, labelTemperature(r, labelS[k])));
  const thermoclineLayers = rho.filter((r, k) => k > 0 && r < THERMOCLINE_DENSITY).length;
  const nu4 = closureCoefficient(meshSpacing, o.closureHours, o.closureSpacing);
  const eddyKappa = eddyDiffusivities(mesh, o.eddyDiffusivity), eddyLimit = o.eddyDiffusivity > 0 ? eddyDiffusionLimit(mesh) : 0;
  const diffusion = hybrid ? 0 : o.diffusivity * mesh.radius * mesh.radius / (o.density * o.specificHeat);
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
  let seaCover = null;

  const OS = seq([['OH', L * C], ['OU', L * E], ['OQ', L * C], ['OW', L * C]]);
  const outlet = runoffOutlets(mesh, geography);
  const drainStart = new Float32Array(C + 1), drainCells = [], drainWeights = [];
  {
    const byOutlet = Array.from({ length: C }, () => []);
    for (let j = 0; j < C; j++) if (!cellOcean[j] && outlet[j] >= 0) byOutlet[outlet[j]].push(j);
    for (let i = 0; i < C; i++) {
      drainStart[i] = drainCells.length;
      for (const j of byOutlet[i]) { drainCells.push(j); drainWeights.push(mesh.areaCell[j] / mesh.areaCell[i]); }
    }
    drainStart[C] = drainCells.length;
  }
  /*
   * The L-sized scratch lives within one tendency (tendency below, in
   * dispatch order): FLUX from oFlux to oCellTendency, LAPA from
   * oClosureFill (or the first oLapVelocity without it) to the second
   * oDivCurl and LAPB from the second oLapVelocity to oMomentum share a
   * slot; AVORT from oVertexVort to oEdgePV, CURLS from each oDivCurl to
   * the oLapVelocity after it and PHI from oKineticPhi, dispatched just
   * before oMomentum, to oMomentum share another. With the closure's
   * rings FLUX marks the fitted edges and LAPA and LAPB hold the rings
   * until oClosureBack2, so those three stay apart. After the combine the
   * eddy and Ekman passes leave their fluxes in EDDYV, EDDYQ and EDDYW over FLUX,
   * HEDGE and QE (over FLUX, LAPA and LAPB with the rings), which
   * oEdgeThickness and every tendency recompute before reading them.
   * DU, the shear at the mixed layer's base from oShearEdge to
   * oMixedLayer, sits in FLUX's slot. UEK, the Ekman transport of each
   * edge, has a slot of its own.
   */
  const sheared = o.shearEntrainment > 0 ? ['DU'] : [];
  const OD = seq([
    ...(o.closureRings ? [[['FLUX', 'EDDYV', ...sheared], L * E], [['LAPA', 'EDDYQ'], L * E], [['LAPB', 'EDDYW'], L * E]] : [[['FLUX', 'LAPA', 'LAPB', 'EDDYV', ...sheared, ...(hybrid ? ['TARGET', 'VISC'] : [])], L * E]]),
    [o.closureRings ? 'HEDGE' : ['HEDGE', 'EDDYQ'], L * E], [['AVORT', 'CURLS', 'PHI'], L * V], [o.closureRings ? 'QE' : ['QE', 'EDDYW'], L * E], ['DIVS', L * C],
    ['RHOML', C], ['REACH', C], ['BUOY', C], ['GRADRHO', E], ['GRADETA', E], ['SLOW', E], ['DEPTHEDGE', E],
    ['ETA', C], ['FRESH', C], ['PREVT0', C], ['PREVVOL', C], ['ICED', C], ['STRESS', E],
    ['SURFT', C], ['SURFICE', C], ['SURFCONC', C], ['SURFACEIN', C], ['T0', C], ['S0', C], ['RAINSEEN', C],
    ['EMASK', E], ['CMASK', C], ['BATH', C], ['FEDGE', E],
    ['RUNOFFSEEN', C], ['DRAINSTART', C + 1], ['DRAINCELL', Math.max(1, drainCells.length)], ['DRAINW', Math.max(1, drainCells.length)],
    ['EDDYK', E], ['TEDGE', 3 * E], ['DEEPEST', E], ['REM', 3 * C],
    ...(o.ekmanDepth > 0 ? [['HM', E]] : []),
    ...(o.retreatDaily ? [['SUNSTEP', C], ['SUNSUM', C], ['DARKSUM', C], ['SUNDAY', C], ['DARKDAY', C]] : []),
    ...(o.ekmanBlock ? [['UEK', E]] : []),
    ...(hybrid ? [['RHOL', L * C], ['PRESS', L * C], ['HABOVE', L * C], ['MIXED', C], ['BLD', C], ['DRIED', C], ['FCELL', C]] : []),
  ]);
  // The barotropic RK4's fixed-offset blocks index into the same OD storage
  // array, so their offsets must continue on from OD's own sections rather
  // than starting a second, colliding index space at 0.
  const B = seq([['BCUR', C + E], ['BK1', C + E], ['BK2', C + E], ['BK3', C + E], ['BTRIALA', C + E], ['BTRIALB', C + E], ['BAVG', C + E]]);
  const bTotal = B.total;
  for (const k of Object.keys(B)) if (k !== 'total') B[k] += OD.total;
  B.total = bTotal;
  const ODTOTAL = OD.total + bTotal;

  const kernels = oceanKernels({ ...o, hybrid, table, L, C, E, V, OS, OD, B, rho, labelT, labelS, nu4, diffusion, omega: mesh.omega, implicitStress: !!core.physics && core.physics.surfaceExchange === 'roughness' && !!core.physics.implicitDrag });
  if (hybrid) for (const name of SLAB_ONLY_KERNELS) delete kernels[name];
  const reduced = hybrid ? HYBRID_REDUCED : OCEAN_REDUCED;
  const OF = seq([['SST', C], ['SSS', C], ['H1', C], ['THD', C], ['ETA', C], ['CUR', 3 * C], ['CSPD', C], ['UPW', C], ['PART', reduced.length * reductionGroups(C)]]);
  kernels.oFrame = `${hybrid ? HYBRID_THERMOCLINE : ''}@compute @workgroup_size(${WORKGROUP}) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
  let i = i32(id.x); if (i >= C) { return; }
  let h0 = max(EPSO, IN[hOff(0) + i]);
  OUT[${OF.SSS} + i] = IN[wOff(0) + i] / h0 + SREF;
${hybrid ? `  OUT[${OF.H1} + i] = OD[O_MIXED + i];
  OUT[${OF.THD} + i] = thermoclineDepthOf(i);` : `  OUT[${OF.H1} + i] = IN[hOff(0) + i];
  var thermo = 0.0;
  for (var k = 0; k <= ${thermoclineLayers}; k++) { thermo += IN[hOff(k) + i]; }
  OUT[${OF.THD} + i] = thermo;`}
  OUT[${OF.ETA} + i] = OD[O_ETA + i];
  let depth = P[4];
  var layer = -1; var top = 0.0;
  for (var k = 0; k < L; k++) { let hk = IN[hOff(k) + i]; if (${hybrid ? 'false' : 'k > 0 && hk <= THINO'}) { continue; } if (layer < 0 && depth < top + hk) { layer = k; } if (layer < 0) { top += hk; } }
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
      if (above >= depth) { break; }
      let he = 0.5 * (IN[hOff(k) + a] + IN[hOff(k) + b]);
      if (${hybrid ? 'false' : 'k > 0 && he <= THINO'}) { continue; }
      transport += IN[uOff(k) + e] * min(he, depth - above);
      above += he;
    }
    upwelling += s * MF[F_DV + e] * transport;
  }
  w = w / MF[F_AREA + i];
  OUT[${OF.CUR} + 3 * i] = w.x; OUT[${OF.CUR} + 3 * i + 1] = w.y; OUT[${OF.CUR} + 3 * i + 2] = w.z;
  OUT[${OF.CSPD} + i] = length(w);
  OUT[${OF.UPW} + i] = upwelling / MF[F_AREA + i];
}`;
  kernels.oReduce = hybrid ? HYBRID_THERMOCLINE + reductionKernel(HYBRID_REDUCED, { count: C, base: OF.PART, setup: hybridReducedSetup }) : reductionKernel(OCEAN_REDUCED, { count: C, base: OF.PART, setup: oceanReducedSetup(thermoclineLayers) });
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
      if (hybrid) dispatch(pass, 'oLayerDensity', g, C);
      else { dispatch(pass, 'oSurfaceDensity', g, C); dispatch(pass, 'oGradRho', g, E); }
      dispatch(pass, 'oEdgeThickness', g, E);
      dispatch(pass, 'oFlux', g, E);
      dispatch(pass, 'oCellTendency', g, C);
      dispatch(pass, 'oVertexVort', g, V);
      dispatch(pass, 'oEdgePV', g, E);
      if (o.closureRings) dispatch(pass, 'oDeepestEdge', g, E);
      if (o.closureFill > 0) dispatch(pass, 'oClosureFill', g, L * E);
      if (o.closureRings) dispatch(pass, 'oClosureRing', g, L * E);
      dispatch(pass, 'oDivCurl', g, C + V);
      dispatch(pass, 'oLapVelocity', g, E);
    });
    setParams([0, 1, relaxRate, stepDt]);
    compute((pass) => {
      dispatch(pass, 'oDivCurl', g, C + V);
      dispatch(pass, 'oLapVelocity', g, E);
      if (o.closureRings) {
        dispatch(pass, 'oClosureBack1', g, L * E);
        dispatch(pass, 'oClosureBack2', g, L * E);
      }
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
    compute((pass) => {
      if (o.ekmanDepth > 0 && !hybrid) dispatch(pass, 'oWindDepth', group(ob.S, ob.T), E);
      dispatch(pass, 'oGradEta', group(ob.S, ob.T), E);
    });
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
    dynamics(dt);
    const g = group(ob.S, ob.T);
    if (eddyLimit > 0 || o.ekmanBlock) setParams([dt, eddyLimit / dt]);
    compute((pass) => {
      rescalePass(pass, g);
      if (eddyLimit > 0) eddyPasses(pass, g);
      if (o.ekmanBlock) ekmanPasses(pass, g);
      shiftPass(pass, g);
    });
    edgesCurrent = true;
  }
  /*
   * step's parts in its order, for the tests that follow heat and salt
   * through it: dynamics(dt), rescale(), eddy(dt), ekman(dt), shift().
   */
  const phases = {
    dynamics(dt) { dynamics(dt); edgesCurrent = false; },
    rescale() { compute((pass) => rescalePass(pass, group(ob.S, ob.T))); edgesCurrent = false; },
    eddy: eddyTransport,
    ekman: ekmanTransport,
    shift() { compute((pass) => shiftPass(pass, group(ob.S, ob.T))); edgesCurrent = true; },
    ...(hybrid ? { regrid, column: columnPhysics } : {}),
  };
  /*
   * The hybrid's ALE step and column physics after the shift: the targets
   * of every column, the edges' velocities and the cells' heat and salt
   * remapped onto them (P[0] the relax), then each column's mixing and
   * each edge's implicit solve (P[0] the step).
   */
  function regrid(dt) {
    if (!hybrid) return;
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
    if (!hybrid) return;
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
  function ekmanPasses(pass, g) {
    dispatch(pass, 'oEkmanFlux', g, E);
    dispatch(pass, 'oEkmanApply', g, L * C);
    dispatch(pass, 'oEkmanLeftover', g, C);
  }
  function ekmanTransport(dt) {
    if (!o.ekmanBlock) return;
    setParams([dt, eddyLimit / dt]);
    compute((pass) => ekmanPasses(pass, group(ob.S, ob.T)));
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
    const masked = new Float32Array(E), cover = (i) => coverOf(ice, concentration, i);
    for (let e = 0; e < E; e++) masked[e] = edgeOcean[e] ? total[e] * (1 - 0.5 * (cover(mesh.cellsOnEdge[2 * e]) + cover(mesh.cellsOnEdge[2 * e + 1])) * (1 - o.iceStressTransmission)) : 0;
    device.queue.writeBuffer(ob.OD, 4 * OD.STRESS, masked);
  }
  function mixedLayer(dt, closing = 0) {
    setParams([0, 0, relaxRate, 0, 0, 0, dt, closing]);
    const g = group(ob.S, ob.T);
    compute((pass) => {
      if (o.carryMomentum && !edgesCurrent) dispatch(pass, 'oEdgeThickness', g, E);
      if (o.shearEntrainment > 0) dispatch(pass, 'oShearEdge', g, E);
      dispatch(pass, 'oMixedReach', g, C);
      dispatch(pass, 'oMixedLayer', g, C);
      if (o.carryMomentum) dispatch(pass, 'oExchangeMomentum', g, E);
    });
    edgesCurrent = o.carryMomentum;
  }
  function salt(dt) {
    setParams([0, 0, 0, 0, 0, 0, dt]);
    compute((pass) => dispatch(pass, 'oSalt', group(ob.S, ob.T), C));
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
  }

  /*
   * A single testable/coupled ocean step at double-precision-matched
   * physics but single-precision arithmetic: readSurface, setStress,
   * step, mixedLayer, salt, writeSurface, unconditionally (the caller
   * decides when to call it; model.gpu.js gates on its own everySteps
   * counter so it only pays for the atmosphere-state read and the step
   * on the steps that need it).
   */
  function dayClosing(time, dt) { return o.retreatDaily ? closeDay(day, time, dt) : 0; }
  async function advance(surfaceT, ice, totalStress, dt, concentration = seaCover, time = null) {
    const closing = dayClosing(time, dt);
    readSurface(surfaceT, ice, concentration);
    setStress(typeof totalStress === 'function' ? totalStress() : totalStress, ice, concentration);
    step(dt);
    if (hybrid) { regrid(dt); columnPhysics(dt); } else mixedLayer(dt, closing);
    salt(dt);
    writeSurface(dt);
    await device.queue.onSubmittedWorkDone();
    return true;
  }
  async function advanceCoupled(dt, time = null) {
    const closing = dayClosing(time, dt);
    readSurfaceFromAtmosphere();
    stressFromAtmosphere();
    step(dt);
    if (hybrid) { regrid(dt); columnPhysics(dt); } else mixedLayer(dt, closing);
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
    const rhoMl = new Float64Array(C);
    for (let i = 0; i < C; i++) { const h0 = Math.max(EPS, h[i]); rhoMl[i] = eos(Q[i] / h0, W[i] / h0); }
    const densityAt = hybrid ? (k, i) => { const n = at(k, i); return h[n] > 1e-6 ? eos(Q[n] / h[n], W[n] / h[n]) : rho[k]; } : (k, i) => (k === 0 ? rhoMl[i] : rho[k]);
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
   * The hybrid ocean's start, buildHybrid of ocean/layered.module.js: every
   * sea column on the depth table clipped to its floor, the atlas's layer
   * means where the atlas covers the cell, the nearest covered column's
   * where it does not, the analytic profile without an atlas; the steric
   * free surface; then one regrid on the device once uploaded.
   */
  function hybridGrid(i, eta, zCol) {
    const H = D[i] + eta[i];
    for (let k = 1; k < L; k++) zCol[k] = Math.min(table.depths[k], H);
    spaceInterfaces(zCol, L, H, o.minimumLayer);
    return H;
  }
  function initializeHybridArrays(surfaceT, ice, atlas = null) {
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
      const iced = ice[i] > 0, lat = mesh.latCell[i];
      T0[i] = iced ? FREEZING_POINT : surfaceT[i];
      if (!cellOcean[i]) continue;
      hybridGrid(i, eta, zCol);
      let column = atlas ? atlas.columnAt(lat, mesh.lonCell[i]) : null;
      if (column) {
        const reached = column.depths.length;
        if (atlas.columnReaching && atlas.depths && reached < atlas.depths.length && atlas.depths[reached] <= D[i]) {
          const deep = atlas.columnReaching(lat, mesh.lonCell[i], D[i]);
          if (deep && deep.depths.length > reached) {
            const Td = Float64Array.from(deep.T), Sd = Float64Array.from(deep.S);
            Td.set(column.T); Sd.set(column.S);
            column = { depths: deep.depths, T: Td, S: Sd };
          }
        }
        filled[i] = 1;
        zs.length = ts.length = ss.length = 0;
        for (let j = 0; j < column.depths.length && column.depths[j] < D[i]; j++) { zs.push(column.depths[j]); ts.push(column.T[j]); ss.push(column.S[j]); }
        const [tb, sb] = profileAt(column, D[i]);
        zs.push(D[i]); ts.push(tb); ss.push(sb);
        for (let k = 0; k < L; k++) { const n = at(k, i), hk = zCol[k + 1] - zCol[k]; h[n] = hk; Q[n] = integral(ts, zCol[k], zCol[k + 1]); W[n] = integral(ss, zCol[k], zCol[k + 1]); }
      } else if (atlas) {
        pending.push(i);
        for (let k = 0; k < L; k++) h[at(k, i)] = zCol[k + 1] - zCol[k];
      } else {
        const s0 = salinityProfile(lat), t0 = T0[i], scale = 400 * (1 - o.thermoclineTilt + 2 * o.thermoclineTilt * Math.cos(lat) ** 2);
        const tDeep = 275.15, sDeep = 34.7, decay = (a, b) => scale * (Math.exp(-a / scale) - Math.exp(-b / scale)) / (b - a);
        for (let k = 0; k < L; k++) {
          const n = at(k, i), hk = zCol[k + 1] - zCol[k], f = decay(zCol[k], zCol[k + 1]);
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
      let z = 0;
      for (let k = 0; k < L; k++) {
        const n = at(k, i), hk = h[n], center = z + 0.5 * hk;
        z += hk;
        let t = T0[i], sal = o.referenceS;
        if (source >= 0) {
          let top = 0, j = 0;
          while (j < L - 1 && top + h[at(j, source)] < center) { top += h[at(j, source)]; j++; }
          const m = at(j, source);
          if (h[m] > 1e-9) { t = Q[m] / h[m]; sal = W[m] / h[m]; }
        }
        Q[n] = hk * t; W[n] = hk * sal;
      }
    }
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) continue;
      T0[i] = ice[i] > 0 ? FREEZING_POINT : Math.max(FREEZING_POINT, Q[i] / h[i]);
      Q[i] = h[i] * T0[i];
      previousT0[i] = T0[i];
    }
    stericSurface(h, Q, W, eta);
    return { h, u, Q, W, eta, T0, previousT0, filled: atlas ? filled : null };
  }
  function initializeArrays(surfaceT, ice, atlas = null) {
    if (hybrid) return initializeHybridArrays(surfaceT, ice, atlas);
    const h = new Float64Array(L * C), u = new Float64Array(L * E), Q = new Float64Array(L * C), W = new Float64Array(L * C);
    const eta = new Float64Array(C), T0 = new Float64Array(C), previousT0 = new Float64Array(C);
    const at = (k, i) => k * C + i;
    const filled = atlas ? atlasColumns(mesh, atlas, { D, cellOcean, ice, rho, labelT, labelS, h, Q, W, T0, shallowestMixedDepth: o.shallowestMixedDepth, maximumMixedDepth: o.maximumMixedDepth }) : null;
    for (let i = 0; i < C; i++) {
      if (filled && filled[i]) { previousT0[i] = T0[i]; continue; }
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
        const [tk, sk] = hk > EPS ? interiorWater(rho[k], labelT[k], labelS[k], lat) : [labelT[k], labelS[k]];
        h[at(k, i)] = hk; Q[at(k, i)] = hk * tk; W[at(k, i)] = hk * sk;
      }
      const excess = cumulative - D[i];
      let giver = L - 1;
      while (giver > 0 && !(h[at(giver, i)] > excess + EPS)) giver--;
      const n = at(giver, i), f = (h[n] - excess) / h[n];
      h[n] -= excess; Q[n] *= f; W[n] *= f;
      previousT0[i] = T0[i];
    }
    stericSurface(h, Q, W, eta);
    return { h, u, Q, W, eta, T0, previousT0, filled };
  }
  function uploadArrays({ h, u, Q, Qref, W, Wref, eta }, surfaceT, ice, restart = null, concentration = seaCover, regridOnDevice = hybrid && !restart) {
    edgesCurrent = false;
    day.span = 0;
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
    device.queue.writeBuffer(ob.OD, 4 * OD.BATH, Float32Array.from(D));
    device.queue.writeBuffer(ob.OD, 4 * OD.FEDGE, Float32Array.from(mesh.fEdge));
    device.queue.writeBuffer(ob.OD, 4 * OD.EDDYK, Float32Array.from(eddyKappa));
    device.queue.writeBuffer(ob.OD, 4 * OD.TEDGE, Float32Array.from(mesh.tEdge));
    if (hybrid) device.queue.writeBuffer(ob.OD, 4 * OD.FCELL, Float32Array.from(mesh.fCell));
    if (hybrid) device.queue.writeBuffer(ob.OD, 4 * OD.MIXED, mixedDepthsOf(h, Q, W));
    device.queue.writeBuffer(ob.OD, 4 * OD.DRAINSTART, drainStart);
    if (drainCells.length) { device.queue.writeBuffer(ob.OD, 4 * OD.DRAINCELL, Float32Array.from(drainCells)); device.queue.writeBuffer(ob.OD, 4 * OD.DRAINW, Float32Array.from(drainWeights)); }
    device.queue.writeBuffer(ob.OD, 4 * OD.PREVVOL, restart ? Float32Array.from(restart.previousIceVolume) : Float32Array.from(ice, (_, i) => coverOf(ice, concentration, i) * ice[i]));
    if (restart) device.queue.writeBuffer(ob.OD, 4 * OD.FRESH, Float32Array.from(restart.fresh));
    if (restart) device.queue.writeBuffer(ob.OD, 4 * OD.BUOY, Float32Array.from(restart.buoyancy));
    if (restart && restart.remainder && restart.remainder.length === 3 * C) device.queue.writeBuffer(ob.OD, 4 * OD.REM, Float32Array.from(restart.remainder));
    if (o.retreatDaily && restart && DAILY.every((name) => restart[name] && restart[name].length === C)) {
      for (const [name, at] of [['sunSum', OD.SUNSUM], ['darkSum', OD.DARKSUM], ['sunDay', OD.SUNDAY], ['darkDay', OD.DARKDAY]]) device.queue.writeBuffer(ob.OD, 4 * at, Float32Array.from(restart[name]));
      day.span = restart.daySpan ? restart.daySpan[0] : 0;
    }
    const capacity = restart ? restart.capacity : Float64Array.from({ length: C }, (_, i) => o.density * o.specificHeat * Math.max(h[i], 1));
    const kept = (name) => (restart && restart[name] && restart[name].length === C ? restart[name] : null);
    core.uploadPhysics({ capacity, oceanFlux: restart ? restart.flux : null, transmission: o.penetration ? kept('transmission') ?? transmissionOf(h, eta) : null, penetrated: kept('penetrated') });
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
  function transmissionOf(h, eta) {
    const out = new Float64Array(C), column = new Float64Array(L);
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) continue;
      if (hybrid) { out[i] = transmittedShare(h[i], o.penetration); continue; }
      for (let k = 0; k < L; k++) column[k] = h[k * C + i];
      out[i] = penetrationShares(column, L, D[i] + eta[i], o.penetration);
    }
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
      if (filled[i] && !(ice[i] > 0)) surfaceT[i] = T0[i];
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
   * A saved hybrid ocean on these labels is uploaded as it is, its columns
   * stretched to this bathymetry's depth; a saved slab ocean is converted as
   * convertSlab of ocean/layered.module.js does, each column remapped by the
   * parabolas onto the depth table at its floor and each edge's velocity
   * piecewise-linearly onto the mean of its cells' new thicknesses, and
   * regridded on the device; any other state starts the ocean afresh.
   */
  function uploadHybrid(saved, surfaceT, ice, concentration) {
    const surfaceAt = (i) => (ice[i] > 0 ? FREEZING_POINT : surfaceT[i]);
    const at = (k, i) => k * C + i;
    const start = () => { const arrays = initializeArrays(surfaceT, ice, o.climatology ?? null); uploadArrays(arrays, surfaceT, ice, null, concentration); device.queue.writeBuffer(ob.OD, 4 * OD.PREVT0, Float32Array.from(arrays.previousT0)); };
    if (!saved || !saved.h) { start(); return; }
    const classes = savedDensities(saved, C), own = classes && sameDensities(classes, labels) && saved.h.length === L * C;
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
  function upload(saved, surfaceT, ice, concentration = seaCover) {
    if (hybrid) return uploadHybrid(saved, surfaceT, ice, concentration);
    if (saved && saved.vertical === 'hybrid') throw new Error("a hybrid ocean state loads into the ocean with vertical: 'hybrid'");
    const classes = savedDensities(saved, C);
    if (classes && !sameDensities(classes, o.densities)) saved = rebinOcean(saved, classes, o.densities, C, { labelT, labelS });
    if (restartable(saved)) {
      uploadArrays(saved, surfaceT, ice, saved);
      device.queue.writeBuffer(ob.OD, 4 * OD.PREVT0, Float32Array.from(saved.previousT0));
      return;
    }
    if (!saved || !saved.h || saved.h.length !== L * C) {
      const arrays = initializeArrays(surfaceT, ice, o.climatology ?? null);
      uploadArrays(arrays, surfaceT, ice, null, concentration);
      device.queue.writeBuffer(ob.OD, 4 * OD.PREVT0, Float32Array.from(arrays.previousT0));
      return;
    }
    const start = () => initializeArrays(surfaceT, ice, o.climatology ?? null);
    const h = Float64Array.from(saved.h), u = Float64Array.from(saved.u), eta = Float64Array.from(saved.eta);
    const Q = new Float64Array(L * C), W = new Float64Array(L * C);
    for (let n = 0; n < L * C; n++) { Q[n] = h[n] * saved.T[n]; W[n] = h[n] * saved.S[n]; }
    const { fromNeighbors } = fitColumns({ h, Q, W, eta }, start, { D, cellOcean, L, C, labelT, labelS, minimumThickness: o.minimumThickness, mesh });
    for (let e = 0; e < E; e++) if (!edgeOcean[e]) for (let k = 0; k < L; k++) u[k * E + e] = 0;
    uploadArrays({ h, u, Q, W, eta }, surfaceT, ice, null, concentration);
    const previousT0 = new Float64Array(C);
    for (let i = 0; i < C; i++) previousT0[i] = Q[i] / Math.max(EPS, h[i]);
    device.queue.writeBuffer(ob.OD, 4 * OD.PREVT0, Float32Array.from(previousT0));
    let surfaced = false;
    for (const i of fromNeighbors) if (!(ice[i] > 0)) { surfaceT[i] = previousT0[i]; surfaced = true; }
    return surfaced;
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
    const h1 = hybrid ? Float64Array.from((await readRanges(device, ob.OD, [{ offset: OD.MIXED, length: C }]))[0]) : h.subarray(0, C), T1 = T.subarray(0, C), S1 = S.subarray(0, C), u1 = u.subarray(0, E);
    const thermoclineDepth = new Float64Array(C);
    for (let i = 0; i < C; i++) {
      if (!cellOcean[i]) { thermoclineDepth[i] = NaN; continue; }
      if (hybrid) { thermoclineDepth[i] = hybridThermocline(h, T, S, i); continue; }
      let depth = 0;
      for (let k = 0; k <= thermoclineLayers; k++) depth += h[k * C + i];
      thermoclineDepth[i] = depth;
    }
    return { h, u, T, S, eta, h1, T1, S1, u1, thermoclineDepth, layers: L };
  }
  function hybridThermocline(h, T, S, i) {
    let z = 0, previous = -Infinity, center = 0;
    for (let k = 0; k < L; k++) {
      const n = k * C + i, hv = h[n], rk = Math.max(previous, hv > 1e-6 ? eos(T[n], S[n]) : rho[k]), c = z + 0.5 * hv;
      if (rk >= THERMOCLINE_DENSITY) return k === 0 ? 0 : center + (THERMOCLINE_DENSITY - previous) / Math.max(rk - previous, 1e-12) * (c - center);
      previous = rk; center = c; z += hv;
    }
    return z;
  }
  function serializeFrom(d) {
    return { h: Array.from(d.h), u: Array.from(d.u), T: Array.from(d.T), S: Array.from(d.S), eta: Array.from(d.eta), densities: Array.from(labels), ...(hybrid ? { vertical: 'hybrid', generator: o.generator } : {}) };
  }
  async function serialize() { return serializeFrom(await download()); }
  /*
   * What upload needs beside serialize()'s fields to put the device back
   * exactly as it stands between two ocean steps: the heat and salt
   * contents as stored, Qref and Wref (upload also takes a Q of h·T and a
   * W of h·S), the mixed layer's last temperature and the ice volume it
   * last saw, the freshwater not yet taken, its remembered surface
   * buoyancy loss, the heat flux and capacity the surface update is
   * using, the sunlight not yet deposited, `penetrated`, and with
   * penetration the share of the sunlight that passes the mixed layer,
   * `transmission`. A saved ocean that carries all of them at this mesh's
   * sizes, the last two optional, its columns within a meter of this
   * bathymetry and dry on land, is uploaded as it is, without fitting its
   * columns.
   */
  async function restartArrays() {
    const PHL = core.layout.PH;
    const [[Qref, Wref], [previousT0, previousIceVolume, fresh, buoyancy, remainder], [flux, capacity, penetrated]] = await Promise.all([
      readRanges(device, ob.S, [{ offset: OS.OQ, length: L * C }, { offset: OS.OW, length: L * C }]),
      readRanges(device, ob.OD, [{ offset: OD.PREVT0, length: C }, { offset: OD.PREVVOL, length: C }, { offset: OD.FRESH, length: C }, { offset: OD.BUOY, length: C }, { offset: OD.REM, length: 3 * C }]),
      readRanges(device, buffers.PH, [{ offset: PHL.OFLUX, length: C }, { offset: PHL.CAP, length: C }, { offset: PHL.PEN, length: C }]),
    ]);
    const transmission = o.penetration ? (await readRanges(device, buffers.PH, [{ offset: PHL.TRANS, length: C }]))[0] : null;
    const days = o.retreatDaily ? await readRanges(device, ob.OD, [OD.SUNSUM, OD.DARKSUM, OD.SUNDAY, OD.DARKDAY].map((offset) => ({ offset, length: C }))) : null;
    const daily = days ? { ...Object.fromEntries(DAILY.map((name, n) => [name, days[n]])), daySpan: [day.span] } : {};
    return Object.fromEntries(Object.entries({ Qref, Wref, previousT0, previousIceVolume, fresh, buoyancy, remainder, flux, capacity, penetrated, ...(transmission ? { transmission } : {}), ...daily }).map(([name, values]) => [name, Float32Array.from(values)]));
  }

  function diagnosticsFrom(sums) {
    const { area, depth, heat, interiorT, interiorH, speed, thermocline, salinity, ssh, transport, limited, thinned = 0, dried } = sums;
    return { oceanUpperDepth: depth / area, oceanHeat: heat / area, oceanInteriorT: interiorH > 0 ? interiorT / interiorH : 0, oceanSpeed: speed, oceanThermoclineDepth: thermocline / area, oceanSalinity: salinity / area, oceanSSH: ssh, oceanTransport: transport / 1e6, oceanLimited: limited, oceanThinned: thinned / area, ...(hybrid ? { oceanDried: dried } : {}) };
  }

  /*
   * Queues the fields the page asked for, those that vary with depth at
   * `depth` metres (0 for the mixed layer), and, when asked, the
   * diagnostics reduction, and reads back only those; fields over land
   * or below the sea floor are NaN and current vectors zero there.
   */
  const FIELDS = { sst: ['SST', 1], sss: ['SSS', 1], layerDepth: ['H1', 1], thermocline: ['THD', 1], ssh: ['ETA', 1], currents: ['CUR', 3], current: ['CSPD', 1], upwelling: ['UPW', 1] };
  function frame({ fields = [], diagnostics: summarize = false, depth = 0 } = {}) {
    const wanted = fields.filter((name) => name in FIELDS);
    if (!wanted.length && !summarize) return Promise.resolve({ fields: {}, diagnostics: null });
    const g = group(ob.S, ob.OF);
    setParams([0, 0, 0, 0, depth]);
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
    layers: L, everySteps: o.everySteps, options: o, D, cellOcean, vertical: hybrid ? 'hybrid' : 'slab', ...(hybrid ? { generator: o.generator, table } : {}),
    useSeaIce(concentration) { seaCover = concentration; },
    initialize, upload, download, serialize, serializeFrom, restartArrays, diagnostics, frame,
    advance, advanceCoupled, accumulateFreshwater, forgetAccumulated,
    setStress, readSurface, readSurfaceFromAtmosphere, stressFromAtmosphere, mixedLayer, dayClosing, salt, writeSurface, step, phases, eddyTransport, ekmanTransport,
    buffers: ob, layout: { OS, OD, B },
  };
}
