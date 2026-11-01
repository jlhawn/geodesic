import { ANVIL_TRACE } from '../physics/moist.module.js';

/*
 * moist.anvilLifetime on the device. The anvil's cover and ice are the PH
 * sections ANVA and ANVQ, K0 × C each, for the layers above cumulusK0. Each
 * piece below ports the JavaScript module's and is patched into the kernels
 * by anchor, which leaves the kernels of a model without the option as they
 * are: the plume's source into adjustPlume, the decay and fall into
 * adjustRain, the radiation's reading of the anvil as cumulus cloud into
 * every kernel that reads cumulusCloud, and the page's cloud totals into
 * frameFields.
 */
export const anvilOn = (o) => o.anvilLifetime != null;

export function anvilConstants(o) {
  if (!anvilOn(o)) return '';
  if (!(o.anvilLifetime > 0 && o.anvilLifetime < Infinity)) throw new Error(`anvilLifetime must be null or a lifetime above 0 in seconds, not ${o.anvilLifetime}`);
  if (o.anvilFall != null && !(o.anvilFall > 0 && o.anvilFall < Infinity)) throw new Error(`anvilFall must be null or a speed above 0 in m/s, not ${o.anvilFall}`);
  if (o.anvilSpread != null && !(o.anvilSpread > 0 && o.anvilSpread < Infinity)) throw new Error(`anvilSpread must be a factor above 0, not ${o.anvilSpread}`);
  return `
const ANVIL_LIFE: f32 = ${o.anvilLifetime}; const ANVIL_TRACE: f32 = ${ANVIL_TRACE}; const ANVIL_FALL_C: f32 = ${o.anvilFall ?? o.iceFall ?? 0}; const ANVIL_SPREAD: f32 = ${o.anvilSpread ?? 1};`;
}

function replaceOnce(body, from, to, where) {
  const at = body.indexOf(from);
  if (at < 0 || body.indexOf(from, at + 1) >= 0) throw new Error(`the anvil's anchor is not in the ${where} kernel exactly once: ${from.slice(0, 60)}`);
  return body.slice(0, at) + to + body.slice(at + from.length);
}

const CLOUD_READ = `  if (!CU_CLOUD || k < CU_K0) { return vec2<f32>(0.0, 0.0); }
`;
const CLOUD_ANVIL = `  if (!CU_CLOUD) { return vec2<f32>(0.0, 0.0); }
  if (k < CU_K0) { return vec2<f32>(PH[PH_ANVA + k * C + i], PH[PH_ANVQ + k * C + i]); }
`;

const PLUME_ARRAYS = `var flux: array<f32, K + 1>;
  var entrained: array<f32, K>;`;
const PLUME_RAIN = `      if (fallen > 0.0) {
        w -= fallen; s += LHEAT * fallen; fallout = fallen;`;
const PLUME_DONE = `  for (var k = bottom - 1; k >= 0; k--) { cuReserve[k] = max(0.0, cuReserve[k + 1] - cuFall[k + 1]); }
  if (PL_MOMENTUM) {`;

const SOURCE = `  for (var k = top; k < CU_K0; k++) {
    let detrained = baseFlux * (flux[k + 1] - flux[k]); let supply = iceCarried[k + 1];
    if (!(detrained > 0.0 && supply > CLEAR_AIR)) { continue; }
    let idx = k * C + i; let rate = detrained * GRAV * dt / layerDp(pi, k);
    let added = min(rate * supply, max(0.0, IN[S_Q + idx]));
    if (!(added > 0.0)) { continue; }
    PH[PH_ANVA + idx] += (1.0 - PH[PH_ANVA + idx]) * min(1.0, rate * ANVIL_SPREAD);
    PH[PH_ANVQ + idx] += added;
    IN[S_Q + idx] -= added;
    IN[S_TH + idx] += LHEAT * added / (CP * D[D_EXM + idx]);
  }
`;

const RAIN_BEFORE = `  for (var k = 0; k < K; k++) {
    let idx = k * C + i;
    let mass = pi * LV[L_DS + k] / GRAV;
    if (rained > 0.0 && (EVAP_IN_CLOUD`;
const FALL = `  var anvilFall = 0.0;
  {
    let decay = exp(-dt / ANVIL_LIFE);
    var fall = 0.0; var fallCover = 0.0;
    for (var k = 0; k < CU_K0; k++) {
      let idx = k * C + i; let mass = pi * LV[L_DS + k] / GRAV; let temperature = IN[S_TH + idx] * D[D_EXM + idx];
      var ice = PH[PH_ANVQ + idx]; var cover = PH[PH_ANVA + idx];
      let arriving = fall / mass; let arrivingCover = fallCover;
      fall = 0.0; fallCover = 0.0;
      if (!(ice > 0.0 || arriving > 0.0)) { PH[PH_ANVQ + idx] = 0.0; PH[PH_ANVA + idx] = 0.0; continue; }
      moved = true;
      if (temperature >= LIQUID_T) {
        fall = (ice + arriving) * mass; fallCover = max(cover, arrivingCover);
        PH[PH_ANVQ + idx] = 0.0; PH[PH_ANVA + idx] = 0.0;
        continue;
      }
      var sublimated = ice * (1.0 - decay);
      ice += arriving - sublimated; cover = max(cover * decay, arrivingCover);
      if (ICE_FALL && ice > 0.0 && cover > 0.0) {
        let speed = ANVIL_FALL_C * pow(pi * LV[L_SM + k] / (RGAS * temperature) * ice / cover, FALL_EXP);
        let courant = speed * dt * LV[L_SM + k] * GRAV / (RGAS * temperature * LV[L_DS + k]); let leaving = courant / (1.0 + courant);
        fall = leaving * ice * mass; fallCover = cover;
        ice -= leaving * ice;
      }
      if (!(ice >= ANVIL_TRACE && cover >= CU_TRACE)) { sublimated += ice; ice = 0.0; cover = 0.0; }
      IN[S_Q + idx] += sublimated;
      IN[S_TH + idx] -= LHEAT * sublimated / (CP * D[D_EXM + idx]);
      PH[PH_ANVQ + idx] = ice; PH[PH_ANVA + idx] = cover;
    }
    anvilFall = fall;
  }
`;
const RAIN_MASS = `    let mass = pi * LV[L_DS + k] / GRAV;
    if (rained > 0.0 && (EVAP_IN_CLOUD`;

const FRAME_TOTAL = `    cloud += path; cumulus += path;
  }
  let deck = PH[PH_DECKF + i] * PH[PH_DECK + i];`;
const FRAME_BASE = `          if (j >= CU_K0) { let cu = (j - CU_K0) * C + i; pj += mass * PH[PH_CUCOVER + cu] * PH[PH_CUWATER + cu]; }`;

/*
 * `bodies` with the anvil in the adjustment kernels, in the kernels that
 * read the cumulus cloud for the radiation and in the page's cloud totals,
 * where the CPU model counts it with the plumes' cumulus.
 */
export function anvilKernels(bodies) {
  const out = { ...bodies };
  out.adjustPlume = replaceOnce(out.adjustPlume, PLUME_ARRAYS, `${PLUME_ARRAYS.replace('var entrained', 'var iceCarried: array<f32, K + 1>; var entrained')}`, 'adjustPlume');
  out.adjustPlume = replaceOnce(out.adjustPlume, PLUME_RAIN, `      iceCarried[k] = (at.y - fallen) * (1.0 - clamp((at.x - ICE_T) / (LIQUID_T - ICE_T), 0.0, 1.0));
${PLUME_RAIN}`, 'adjustPlume');
  out.adjustPlume = replaceOnce(out.adjustPlume, PLUME_DONE, `${SOURCE}${PLUME_DONE}`, 'adjustPlume');
  out.adjustRain = replaceOnce(out.adjustRain, RAIN_BEFORE, `${FALL}${RAIN_BEFORE.replace(RAIN_MASS, `    let mass = pi * LV[L_DS + k] / GRAV;
    if (k == CU_K0) { descending += anvilFall; }
    if (rained > 0.0 && (EVAP_IN_CLOUD`)}`, 'adjustRain');
  for (const name of ['physics', 'radiation', 'layered']) if (out[name]) out[name] = replaceOnce(out[name], CLOUD_READ, CLOUD_ANVIL, name);
  if (out.frameFields) {
    out.frameFields = replaceOnce(out.frameFields, FRAME_TOTAL, `    cloud += path; cumulus += path;
  }
  for (var j = 0; j < CU_K0; j++) {
    let path = pi * LV[L_DS + j] / GRAV * PH[PH_ANVQ + j * C + i];
    cloud += path; cumulus += path;
  }
  let deck = PH[PH_DECKF + i] * PH[PH_DECK + i];`, 'frameFields');
    out.frameFields = replaceOnce(out.frameFields, FRAME_BASE, `${FRAME_BASE}
          else { pj += mass * PH[PH_ANVQ + j * C + i]; }`, 'frameFields');
  }
  return out;
}
