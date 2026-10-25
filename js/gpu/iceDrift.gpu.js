/*
 * The sea-ice drift of ocean/iceDrift.module.js on the GPU, kernel for
 * kernel and in the same order within an ocean call: iPrepareCells and
 * iPrepareEdges (the cell totals, the edge forcing and masks), then the
 * pseudo-time iteration, a cell kernel for σ1 and an edge kernel for the
 * velocity, ping-ponging between UICE and UALT within one compute pass,
 * iStressOut (the ice-water stress the ocean feels), and the transport,
 * iFlux on the edges and iUpdate on the cells, which carries the ice's
 * area, volume, snow, skin enthalpy and snow albedo, shifts the ocean's
 * stored ice volume, closes the cover, melts away what the transport left
 * under the sea ice's removal floors or out of ICE_REMNANT's bounds and flags
 * the cells the ice reached while open. The cells'
 * totals are gathered from the edge fluxes by each cell's own six slots, the
 * form that needs no atomics, where the CPU scatters them from the edges.
 *
 * The WGSL rules of CLAUDE.md hold throughout: loops count up over the
 * fixed MAXE and MAXEE slots with a zero weight or sign for the unused
 * ones, nothing is held in a thread-private array, and every running sum is
 * a scalar. The Courant subcycling is fixed by the host from the speed
 * limit, where the CPU measures the fastest edge.
 */
const lit = (x) => Number(x).toExponential(10);
const constant = (name, value) => `const ${name}: f32 = ${lit(value)};`;

export const DRIFT_EDGE_SECTIONS = ['UICE', 'UALT', 'ICESTRESS', 'VALID', 'IDYN', 'IACT', 'IUPD', 'IAE', 'IMU', 'IUW', 'IUWT', 'ITILT', 'ITAU', 'IUN', 'IFLX', 'IFA', 'IFV', 'IFS', 'IFE', 'IFSA'];
export const DRIFT_CELL_SECTIONS = ['SIGMA', 'WATERT', 'ARRIVED', 'IA', 'IV', 'ISN', 'IEN', 'ISA', 'IMASS', 'IPRESS', 'IA0', 'IV0', 'ITOUCH'];
export const driftSections = (C, E) => [...DRIFT_EDGE_SECTIONS.map((name) => [name, E]), ['TAUA', E], ...DRIFT_CELL_SECTIONS.map((name) => [name, C])];

/*
 * The thermodynamics' constants and the atmosphere stress the drift reads:
 * `thermo` the sea ice's skinHeatCapacity, latentHeatFusion, iceDensity and
 * freshSnowAlbedo as ice.module.js has them, `given` the stress as a
 * section of the ocean's scratch the host fills, else the surface layer's
 * as oStressFromAtmosphere computes it.
 */
export function iceDriftHead({ options: d, waterDensity, thermo, minimumConcentration, removal, remnant, freezing, implicitStress, given }) {
  const turn = d.waterTurn * Math.PI / 180;
  return `
${constant('ID_STRENGTH', d.strength)} ${constant('ID_DECAY', d.strengthDecay)} ${constant('ID_DRAG', waterDensity * d.waterDrag)} ${constant('ID_COS', Math.cos(turn))} ${constant('ID_SIN', Math.sin(turn))}
${constant('ID_DMIN', d.minimumDeformation)} ${constant('ID_LIMIT', d.speedLimit)} ${constant('ID_FLOOR', d.minimumRelativeSpeed)} ${constant('ID_DYNCONC', d.dynamicConcentration)} ${constant('ID_MAXCONC', d.maximumConcentration)}
${constant('ID_SKINC', thermo.skinHeatCapacity)} ${constant('ID_LATENT', thermo.iceDensity * thermo.latentHeatFusion)} ${constant('ID_LFUS', thermo.latentHeatFusion)} ${constant('ID_FRESH', thermo.freshSnowAlbedo)}
${constant('ID_MINCONC', minimumConcentration)} ${constant('ID_GONECONC', removal.concentration)} ${constant('ID_GONEVOL', removal.volume)} ${constant('ID_FREEZE', freezing)}
${constant('ID_MAXH', remnant.thickness)} ${constant('ID_MAXSNOW', remnant.snow)} ${constant('ID_EMIN', thermo.skinHeatCapacity * (remnant.skinLow - freezing))} ${constant('ID_EMAX', thermo.skinHeatCapacity * (remnant.skinHigh - freezing))}
fn openCell(i: i32) -> bool { return OD[O_CMASK + i] > 0.5 && OD[O_SHELF + i] < 0.5; }
fn iceCover(i: i32) -> f32 {
  let c = PH[PH_CONC + i];
  return select(0.0, select(1.0, c, c > 0.0), S[S_ICE + i] > 0.0 && openCell(i));
}
fn hiOf(x: i32) -> f32 { let a = OD[O_IA + x]; return select(0.0, OD[O_IV + x] / a, a > 0.0); }
fn siOf(x: i32) -> f32 { let a = OD[O_IA + x]; return select(0.0, OD[O_ISN + x] / a, a > 0.0); }
fn diOf(x: i32) -> f32 { let a = OD[O_IA + x]; return select(0.0, OD[O_IEN + x] / a, a > 0.0); }
fn aiOf(x: i32) -> f32 { let sn = OD[O_ISN + x]; return select(ID_FRESH, OD[O_ISA + x] / sn, sn > 0.0); }
fn rawStress(e: i32, a: i32, b: i32) -> f32 {
${given ? '  return OD[O_TAUA + e];' : `${implicitStress ? '  if (PH[PH_STRESSOK] > 0.5) { return PH[PH_STRESS + e]; }\n' : ''}  let bottom = (K - 1) * C;
  let rhoA = S[S_PI + a] * LV[L_SM + K - 1] / (RGAS * S[S_TH + bottom + a] * D[D_EXM + bottom + a]);
  let rhoB = S[S_PI + b] * LV[L_SM + K - 1] / (RGAS * S[S_TH + bottom + b] * D[D_EXM + bottom + b]);
  let fa = PH[PH_DRAG + a] * rhoA * max(D[D_WIND + a], GUSTO); let fb = PH[PH_DRAG + b] * rhoB * max(D[D_WIND + b], GUSTO);
  return 0.5 * (fa + fb) * S[S_U + (K - 1) * E + e];`}
}
`;
}

/*
 * oStressFromAtmosphere with the drift: (1 − Ā) τ_a + Ā τ_io on the edges
 * the solver covered this call, the transmission formula elsewhere, times
 * the cavity factor.
 */
export const driftStressKernel = (K, idx) => `${K}  let e = ${idx}; if (e >= E) { return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  if (OD[O_EMASK + e] < 0.5) { OD[O_STRESS + e] = 0.0; OD[O_TAUA + e] = 0.0; return; }
  let coverA = select(0.0, select(1.0, PH[PH_CONC + a], PH[PH_CONC + a] > 0.0), OD[O_ICED + a] > 0.5);
  let coverB = select(0.0, select(1.0, PH[PH_CONC + b], PH[PH_CONC + b] > 0.0), OD[O_ICED + b] > 0.5);
  let area = 0.5 * (coverA + coverB);
  let raw = rawStress(e, a, b);
  OD[O_TAUA + e] = raw;
  let open = 1.0 - 0.5 * (OD[O_SHELF + a] + OD[O_SHELF + b]);
  OD[O_STRESS + e] = select(raw * (1.0 - area * (1.0 - TRANSMIT)), (1.0 - area) * raw + area * OD[O_ICESTRESS + e], OD[O_VALID + e] > 0.5) * open;
}`;

export function iceDriftKernels({ K, idx, options: d, third }) {
  const moving = d.rheology !== 'frozen', cavitating = d.rheology === 'cavitating';
  const sinTurn = 'select(select(0.0, -ID_SIN, fe < 0.0), ID_SIN, fe > 0.0)';
  const tangential = (field) => `  var tang = 0.0;
  for (var s = 0; s < MAXEE; s++) { let slot = MAXEE * e + s; tang += MF[F_PVW + slot] * ${field(`MI[EOE + slot]`)}; }
  tang = tang / dc;`;
  const kernels = {
    iClearArrived: `${K}  let i = ${idx}; if (i >= C) { return; }
  OD[O_ARRIVED + i] = 0.0;
}`,
    iPrepareCells: `${K}  let i = ${idx}; if (i >= C) { return; }
  OD[O_ARRIVED + i] = 0.0; OD[O_ITOUCH + i] = 0.0;
  let a = iceCover(i);
  if (a > 0.0) {
    let v = a * S[S_ICE + i]; let sn = a * PH[PH_SNOW + i];
    OD[O_IA + i] = a; OD[O_IV + i] = v; OD[O_ISN + i] = sn; OD[O_IA0 + i] = a; OD[O_IV0 + i] = v;
    OD[O_IEN + i] = a * ID_SKINC * (S[S_TS + i] - ID_FREEZE); OD[O_ISA + i] = sn * PH[PH_SNOWALB + i];
    OD[O_IMASS + i] = ICEDENS * v + sn;
${cavitating ? '    OD[O_IPRESS + i] = ID_STRENGTH * v * exp(-ID_DECAY * (1.0 - a));\n' : ''}  } else {
    OD[O_IA + i] = 0.0; OD[O_IV + i] = 0.0; OD[O_ISN + i] = 0.0; OD[O_IEN + i] = 0.0; OD[O_ISA + i] = 0.0; OD[O_IMASS + i] = 0.0; OD[O_IA0 + i] = 0.0; OD[O_IV0 + i] = 0.0;
    OD[O_SIGMA + i] = 0.0;
  }
}`,
    iPrepareEdges: `${K}  let e = ${idx}; if (e >= E) { return; }
  OD[O_VALID + e] = 0.0; OD[O_IACT + e] = 0.0; OD[O_IUPD + e] = 0.0;
  if (OD[O_IDYN + e] < 0.5) { return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  let area = OD[O_IA + a] + OD[O_IA + b];
  if (!(area > 0.0)) {
${moving ? '    OD[O_UICE + e] = 0.0;\n' : ''}    return;
  }
  OD[O_IACT + e] = 1.0;
${moving ? `  let dc = MF[F_DC + e];
  let ae = 0.5 * area;
  OD[O_IAE + e] = ae; OD[O_IMU + e] = (OD[O_IMASS + a] + OD[O_IMASS + b]) / area;
${tangential((other) => `IN[uOff(0) + ${other}]`)}
  let uw = IN[uOff(0) + e];
  OD[O_IUW + e] = uw; OD[O_IUWT + e] = tang;
  OD[O_ITILT + e] = OGRAV * (OD[O_ETA + b] - OD[O_ETA + a]) / dc;
  OD[O_ITAU + e] = rawStress(e, a, b);
  if (ae < ID_DYNCONC) { OD[O_UICE + e] = uw; } else { OD[O_IUPD + e] = 1.0; OD[O_IUN + e] = OD[O_UICE + e]; }
` : ''}}`,
  };
  if (moving) {
    const momentum = (cur, next) => `${K}  let e = ${idx}; if (e >= E) { return; }
  let c = OD[O_${cur} + e];
  if (OD[O_IUPD + e] < 0.5) { OD[O_${next} + e] = c; return; }
  let dc = MF[F_DC + e]; let fe = OD[O_FEDGE + e];
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  var tang = 0.0; var cor = 0.0;
  for (var s = 0; s < MAXEE; s++) {
    let slot = MAXEE * e + s; let other = MI[EOE + slot]; let w = MF[F_PVW + slot] * OD[O_${cur} + other];
    tang += w; cor += w * (fe + OD[O_FEDGE + other]);
  }
  tang = tang / dc; cor = cor * 0.5 / dc;
  let uw = OD[O_IUW + e]; let dn = uw - c; let dtn = OD[O_IUWT + e] - tang;
  let dr = ID_DRAG * max(sqrt(dn * dn + dtn * dtn), ID_FLOOR);
  let sinT = ${sinTurn};
  let mu = OD[O_IMU + e]; let k = P[0] / mu; let drag = k * dr * ID_COS; let bet = max(P[2], drag);
  ${cavitating ? 'let force = (OD[O_SIGMA + b] - OD[O_SIGMA + a]) / (2.0 * dc * OD[O_IAE + e]);' : 'let force = 0.0;'}
  let rhs = OD[O_ITAU + e] + dr * (uw * ID_COS - dtn * sinT) + mu * (cor - OD[O_ITILT + e]) + force;
  let un = (bet * c + OD[O_IUN + e] + k * rhs) / (1.0 + bet + drag);
  OD[O_${next} + e] = clamp(un, -ID_LIMIT, ID_LIMIT);
}`;
    kernels.iMomentumAB = momentum('UICE', 'UALT');
    kernels.iMomentumBA = momentum('UALT', 'UICE');
    kernels.iCopy = `${K}  let e = ${idx}; if (e >= E) { return; }
  OD[O_UICE + e] = OD[O_UALT + e];
}`;
    kernels.iStressOut = `${K}  let e = ${idx}; if (e >= E) { return; }
  if (OD[O_IACT + e] < 0.5) { return; }
  let dc = MF[F_DC + e]; let fe = OD[O_FEDGE + e];
${tangential((other) => `OD[O_UICE + ${other}]`)}
  let dn = OD[O_IUW + e] - OD[O_UICE + e]; let dtn = OD[O_IUWT + e] - tang;
  let dr = ID_DRAG * max(sqrt(dn * dn + dtn * dtn), ID_FLOOR);
  OD[O_ICESTRESS + e] = -dr * (dn * ID_COS - dtn * ${sinTurn});
  OD[O_VALID + e] = 1.0;
}`;
  }
  if (cavitating) {
    const stress = (cur) => `${K}  let i = ${idx}; if (i >= C) { return; }
  if (OD[O_IA + i] <= 0.0) { return; }
  var div = 0.0;
  for (var m = 0; m < MAXE; m++) { let e = MI[EOC + MAXE * i + m]; div += f32(MI[ESC + MAXE * i + m]) * MF[F_DV + e] * OD[O_${cur} + e]; }
  div = div / MF[F_AREA + i];
  let want = select(0.0, 2.0 * OD[O_IPRESS + i] * div / max(-div, ID_DMIN), div < 0.0);
  OD[O_SIGMA + i] += (want - OD[O_SIGMA + i]) / P[1];
}`;
    kernels.iStressA = stress('UICE');
    kernels.iStressB = stress('UALT');
  }
  kernels.iFlux = `${K}  let e = ${idx}; if (e >= E) { return; }
  OD[O_IFLX + e] = 0.0; OD[O_IFA + e] = 0.0; OD[O_IFV + e] = 0.0; OD[O_IFS + e] = 0.0; OD[O_IFE + e] = 0.0; OD[O_IFSA + e] = 0.0;
  if (OD[O_IACT + e] < 0.5) { return; }
  let u = OD[O_UICE + e];
  if (u == 0.0) { return; }
  let a = MI[COE + 2 * e]; let b = MI[COE + 2 * e + 1];
  let forward = u > 0.0;
  let up = select(b, a, forward); let down = select(a, b, forward);
  let wet = OD[O_IA + up];
  if (!(wet > 0.0)) { return; }
  var cover = wet; var thick = hiOf(up); var snowy = siOf(up); var skin = diOf(up); var albedo = aiOf(up);
${third ? `  let beyond = i32(OD[O_BEY + 2 * e + select(1, 0, forward)]);
  if (openCell(beyond)) {
    let ratio = MF[F_DC + e] / MF[F_DC + i32(OD[O_BED + 2 * e + select(1, 0, forward)])];
    cover = max(0.0, limitedValue(wet, OD[O_IA + down], OD[O_IA + beyond], ratio));
    if (OD[O_IA + down] > 0.0 && OD[O_IA + beyond] > 0.0) {
      thick = limitedValue(thick, hiOf(down), hiOf(beyond), ratio); snowy = limitedValue(snowy, siOf(down), siOf(beyond), ratio);
      skin = limitedValue(skin, diOf(down), diOf(beyond), ratio); albedo = limitedValue(albedo, aiOf(down), aiOf(beyond), ratio);
    }
  }
` : ''}  let fA = u * cover;
  OD[O_IFLX + e] = 1.0; OD[O_IFA + e] = fA; OD[O_IFV + e] = fA * thick; OD[O_IFS + e] = fA * snowy; OD[O_IFE + e] = fA * skin; OD[O_IFSA + e] = fA * snowy * albedo;
}`;
  const update = (last) => `${K}  let i = ${idx}; if (i >= C) { return; }
  var moved = false;
  var dA = 0.0; var dV = 0.0; var dS = 0.0; var dE = 0.0; var dSa = 0.0;
  for (var m = 0; m < MAXE; m++) {
    let e = MI[EOC + MAXE * i + m]; let sg = f32(MI[ESC + MAXE * i + m]);
    if (sg != 0.0 && OD[O_IFLX + e] > 0.5) {
      moved = true;
      let scale = sg * P[3] * MF[F_DV + e];
      dA -= scale * OD[O_IFA + e]; dV -= scale * OD[O_IFV + e]; dS -= scale * OD[O_IFS + e]; dE -= scale * OD[O_IFE + e]; dSa -= scale * OD[O_IFSA + e];
    }
  }
  if (!(moved${last ? ' || OD[O_ITOUCH + i] > 0.5' : ''})) { return; }
  let omega = MF[F_AREA + i];
  let A1 = max(0.0, OD[O_IA + i] + dA / omega); let V1 = max(0.0, OD[O_IV + i] + dV / omega); let S1 = max(0.0, OD[O_ISN + i] + dS / omega);
  let E1 = OD[O_IEN + i] + dE / omega; let Sa1 = max(0.0, OD[O_ISA + i] + dSa / omega);
${last ? `  OD[O_PREVVOL + i] += V1 - OD[O_IV0 + i];
  let started = OD[O_IA0 + i];
  let full = min(A1, ID_MAXCONC);
  let tidy = full >= ID_MINCONC || (V1 <= ID_MAXH * full && S1 <= ID_MAXSNOW * full && E1 >= ID_EMIN * full && E1 <= ID_EMAX * full);
  if (full < ID_GONECONC || V1 < ID_GONEVOL || !tidy) {
    let base = select(S[S_TS + i], ID_FREEZE, started > 0.0);
    S[S_TS + i] = base + (E1 - ID_LATENT * V1 - ID_LFUS * S1) / PH[PH_CAP + i];
    S[S_ICE + i] = 0.0; PH[PH_SNOW + i] = 0.0; PH[PH_SNOWALB + i] = ID_FRESH; PH[PH_CONC + i] = 0.0;
    return;
  }
  if (started == 0.0) { OD[O_WATERT + i] = S[S_TS + i]; OD[O_ARRIVED + i] = 1.0; }
  S[S_ICE + i] = V1 / full; PH[PH_CONC + i] = full; PH[PH_SNOW + i] = S1 / full;
  S[S_TS + i] = ID_FREEZE + E1 / (full * ID_SKINC);
  PH[PH_SNOWALB + i] = select(ID_FRESH, min(1.0, Sa1 / S1), S1 > 0.0);
}` : `  OD[O_IA + i] = A1; OD[O_IV + i] = V1; OD[O_ISN + i] = S1; OD[O_IEN + i] = E1; OD[O_ISA + i] = Sa1; OD[O_ITOUCH + i] = 1.0;
}`}`;
  kernels.iUpdateMid = update(false);
  kernels.iUpdateLast = update(true);
  return kernels;
}
