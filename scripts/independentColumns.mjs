// The day-mean sunlight a saved state's cloud reflects under three readings
// of the same cloud, on the CPU engine after one step (ocean off):
//   M=128 FSD=0 node scripts/independentColumns.mjs <state.bin>
// Each column's cloud is read as the radiation reads it (the layer covers
// of radiation.layerCover, the visibility rule, the plumes' cumulus, the
// deck) and reflects the sun at 24 hour angles of the state's declination
// in one band through the radiation's own two-stream (stream(): τ'/(τ' +
// 2μ), the beam, the multiple reflections with the surface, the cloud's
// absorption), no gases or Rayleigh scattering. "pooled": every layer's
// water spread over the overlapped cover, as the shortwave does without
// cloudColumns. "independent": M columns per cell drawn from the same
// exponential-random overlap as a Markov chain on each layer's cloud
// indicator whose pair covers are those of overlapped(), each layer's
// water in its own cover, with FSD gamma-distributed water in each cloudy
// layer of each column. "moment-matched": the lognormal of the chain's
// first two moments at GH (4, 6 or 8) Gauss-Hermite values, as the
// shortwave does with cloudColumns. Printed as cloud reflection (all-sky
// less clear-sky, W/m2), globally, by latitude band and by the columns
// with cloud only below 680 hPa, only above, or both.
import { readFileSync } from 'node:fs';
import { topographyFromInt16 } from '../js/geography.module.js';
import { createModel, STATE_NAMES } from '../js/model.module.js';
import { decodeState, savedLevels } from '../js/stateFile.module.js';
import { savedDeckField, DECK_FIELDS } from '../js/physics/regrid.module.js';
import { CLOUD_OPTICS, cloudOptics, overlapped, sunDistance, SOLAR_CONSTANT, DECORRELATION_LENGTH, DECORRELATION_SLOPE, LOW_CLOUD_PRESSURE } from '../js/physics/radiation.module.js';
import { stateGrid } from '../js/stateGrid.module.js';

const FILE = process.argv[2], M = Number(process.env.M ?? 128), FSD = Number(process.env.FSD ?? 0);
const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const model = createModel(stateGrid(saved, topography).mesh, { topography, levels: savedLevels(saved), ocean: false, radiation: { clearSkyPass: true } });
const { mesh, core, state, radiation, boundaryLayer: bl, moist, seaIce, land, geography } = model;
const { K, C, sigmaMid, dSigma, g, exnerLayer, geopotential } = core.diagnostics;
STATE_NAMES.forEach((name, a) => state[a].set(saved[name]));
seaIce.load(state[6], saved.concentration ?? null);
for (const field of Object.keys(DECK_FIELDS)) radiation[field].set(savedDeckField(saved, field, model));
land.load({ soil: Float64Array.from(saved.land.soil), snow: Float64Array.from(saved.land.snow), ...(saved.land.vegetation ? { vegetation: Float64Array.from(saved.land.vegetation) } : {}), ...(saved.land.surface ? { surface: Float64Array.from(saved.land.surface) } : {}) });
model.time = saved.time;
const zs = Float64Array.from({ length: C }, (_, i) => (model.surfaceGeopotential ? model.surfaceGeopotential[i] / g : 0));
if (saved.boundaryDepth) bl.depth.set(Float64Array.from(saved.boundaryDepth, (z, i) => z + zs[i]));
if (saved.mixingTop) bl.mixingTop.set(Float64Array.from(saved.mixingTop, (z, i) => z + zs[i]));
if (saved.boundaryRegime) bl.regime.set(saved.boundaryRegime);
if (saved.boundaryBuoyancy) bl.buoyancyFlux.set(saved.boundaryBuoyancy);
model.step(1350 * 16 / saved.N);
const [pi, theta, , surfaceT, q, qc, ice] = state;
const landMask = geography.land, area = mesh.areaCell, deg = 180 / Math.PI;
const fluxT = new Float64Array(C), adir = new Float64Array(C), adif = new Float64Array(C), openSea = new Float64Array(C);
const fluxState = [pi, theta, state[2], fluxT, q, qc];
const out = [null, new Float64Array(K * C), null, null, new Float64Array(K * C)];
for (let i = 0; i < C; i++) {
  fluxT[i] = surfaceT[i];
  if (landMask[i]) { adir[i] = adif[i] = land.albedo(i); openSea[i] = 0; continue; }
  const h = ice[i], cover = seaIce.cover(i, h);
  adir[i] = adif[i] = seaIce.albedo(h, null, seaIce.snow[i], cover);
  openSea[i] = 1 - cover;
  if (h > 0 && cover < 1) fluxT[i] = cover * surfaceT[i] + (1 - cover) * 271.35;
}
radiation.setTime(model.time);
const sun = radiation.sun, dec = Math.asin(sun[2]);
const cuCover = moist.cumulusCover, cuWater = moist.cumulusWater;
const S0 = SOLAR_CONSTANT * sunDistance(model.time), HOURS = 24, ABS = 0.4, DIFFUSE_MU = 0.6;
function toa(tau, path, mu, a) {
  if (!(mu > 0)) return 0;
  const keep = Math.exp(-ABS * path);
  const R = tau > 0 ? tau / (tau + 2 * mu) : 0;
  let direct = tau > 0 ? Math.exp(-tau / mu) : 1, diffuse = 1 - R - direct, returned = tau > 0 ? tau / (tau + 2 * DIFFUSE_MU) : 0;
  direct *= keep; diffuse *= keep; returned *= keep;
  const upward = a * direct + a * diffuse;
  const absorbed = (1 - a) * (direct + diffuse + returned * upward / (1 - a * returned));
  const cloud = (1 - keep) * (1 + upward / (1 - a * returned));
  return 1 - absorbed - cloud;
}
let seed = 12345;
const rand = () => { seed = (seed * 1664525 + 1013904223) >>> 0; return seed / 4294967296; };
const gammaDraw = (nu) => { // Marsaglia-Tsang, mean 1
  const d = nu - 1 / 3, c = 1 / Math.sqrt(9 * d);
  for (;;) { let x, v; do { const u1 = rand(), u2 = rand(); x = Math.sqrt(-2 * Math.log(u1 + 1e-300)) * Math.cos(2 * Math.PI * u2); v = 1 + c * x; } while (v <= 0); v = v * v * v; const u = rand(); if (Math.log(u + 1e-300) < 0.5 * x * x + d - d * v + d * Math.log(v)) return d * v / nu; }
};
const GH8 = [[-2.930637420257244, 1.996040722113676e-4], [-1.981656756695843, 1.707798300741348e-2], [-1.157193712446780, 2.078023258148919e-1], [-0.381186990207322, 6.611470125582413e-1], [0.381186990207322, 6.611470125582413e-1], [1.157193712446780, 2.078023258148919e-1], [1.981656756695843, 1.707798300741348e-2], [2.930637420257244, 1.996040722113676e-4]];
const GH4 = [[-1.650680123885785, 0.08131283544724518], [-0.5246476232752903, 0.8049140900055128], [0.5246476232752903, 0.8049140900055128], [1.650680123885785, 0.08131283544724518]];
const GH6 = [[-2.350604973674492, 0.004530009905508846], [-1.335849074013697, 0.1570673203228566], [-0.4360774119276165, 0.7246295952243925], [0.4360774119276165, 0.7246295952243925], [1.335849074013697, 0.1570673203228566], [2.350604973674492, 0.004530009905508846]];
const GH = { 4: GH4, 6: GH6, 8: GH8 }[process.env.GH ?? 8];
const nu = FSD > 0 ? 1 / (FSD * FSD) : 0;
const optics = { liquid: 0, visible: 0, solar: 0, infrared: 0 };
const tauK = new Float64Array(K), pathK = new Float64Array(K), fK = new Float64Array(K), alphaK = new Float64Array(K), lowK = new Uint8Array(K);
const sums = { ln: 0, clear: 0, model: 0, ica: 0, area: 0, modelNoDeck: 0, icaNoDeck: 0 };
const bands = Array.from({ length: 6 }, () => ({ clear: 0, model: 0, ica: 0, area: 0 }));
const kinds = { lowOnly: { model: 0, ica: 0, area: 0 }, layered: { model: 0, ica: 0, area: 0 }, upperOnly: { model: 0, ica: 0, area: 0 } };
const sub = new Float64Array(M), subPath = new Float64Array(M);
const t0 = performance.now();
for (let i = 0; i < C; i++) {
  radiation.apply(fluxState, out, model.surface.windSpeed, null, i, i + 1, adir, adif, null, openSea, bl.depth, 0);
  const fraction = radiation.stratusFraction[i], deck = radiation.stratus[i];
  const lat = mesh.latCell[i], a = adif[i];
  const decorrelation = DECORRELATION_LENGTH - DECORRELATION_SLOPE * Math.abs(lat * deg);
  let tauSum = 0, pathSum = 0, cumulative = 0, previous = 0, deckTau = 0, hasLow = false, hasUpper = false;
  for (let k = 0; k < K; k++) {
    const idx = k * C + i, m = pi[i] * dSigma[k] / g;
    cloudOptics(theta[idx] * exnerLayer[idx], landMask[i] === 1 && !(geography.iceSheet && geography.iceSheet[i]), CLOUD_OPTICS, optics);
    const resolved = Math.max(0, qc[idx]) * m, cumulus = cuCover[idx] * cuWater[idx] * m, w = resolved + cumulus;
    const cover = radiation.layerCover[k];
    let seen = w > 0 ? cover * -Math.expm1(-w / 1e-3) : 0;
    if (cumulus > 0) seen = Math.max(seen, cuCover[idx] * -Math.expm1(-cumulus / (cuCover[idx] * 1e-3)));
    alphaK[k] = k > 0 ? Math.exp(-(geopotential[(k - 1) * C + i] - geopotential[idx]) / (g * decorrelation)) : 0;
    cumulative = overlapped(cumulative, previous, seen, alphaK[k]);
    previous = seen;
    fK[k] = seen; tauK[k] = optics.solar * w; pathK[k] = w;
    tauSum += tauK[k]; pathSum += w;
    if (seen > 0.01) { if (pi[i] * sigmaMid[k] > LOW_CLOUD_PRESSURE) hasLow = true; else hasUpper = true; }
    if (k === radiation.stratusLayer) deckTau = optics.solar * deck;
  }
  const cover = cumulative > 0 ? cumulative : 1;
  let A0 = 1, A1 = 0, B0 = 0, B1 = 0, D0 = 0, D1 = 0, P0 = 0, P1 = 0;
  for (let k = 0; k < K; k++) {
    const f = fK[k], fp = k > 0 ? fK[k - 1] : 0;
    let p11 = 0, p10 = f;
    if (fp > 0) { const pair = alphaK[k] * Math.max(fp, f) + (1 - alphaK[k]) * (fp + f - fp * f), joint = fp + f - pair; p11 = joint / fp; p10 = fp < 1 ? (f - joint) / (1 - fp) : 0; }
    const t = f > 0 ? tauK[k] / f : 0, w = f > 0 ? pathK[k] / f : 0;
    const a1 = A1 * p11 + A0 * p10, a0 = A1 * (1 - p11) + A0 * (1 - p10);
    const b1 = B1 * p11 + B0 * p10, b0 = B1 * (1 - p11) + B0 * (1 - p10);
    const d1 = D1 * p11 + D0 * p10, d0 = D1 * (1 - p11) + D0 * (1 - p10);
    const q1 = P1 * p11 + P0 * p10, q0 = P1 * (1 - p11) + P0 * (1 - p10);
    A1 = a1; A0 = a0; B0 = b0; D0 = d0; P0 = q0;
    B1 = b1 + t * a1; D1 = d1 + 2 * t * b1 + (1 + FSD * FSD) * t * t * a1; P1 = q1 + w * a1;
  }
  const meanT = cumulative > 0 ? (B0 + B1) / cumulative : 0, cv2 = cumulative > 0 && B0 + B1 > 0 ? Math.max(0, cumulative * (D0 + D1) / ((B0 + B1) ** 2) - 1) : 0, pathPerT = B0 + B1 > 0 ? (P0 + P1) / (B0 + B1) : 0;
  const sig = Math.sqrt(Math.log(1 + cv2));
  for (let s = 0; s < M; s++) {
    let t = 0, p = 0, prev = false;
    for (let k = 0; k < K; k++) {
      const f = fK[k];
      if (!(f > 0)) { prev = false; continue; }
      const fp = k > 0 ? fK[k - 1] : 0;
      let pCloudy;
      if (fp > 0) {
        const pair = alphaK[k] * Math.max(fp, f) + (1 - alphaK[k]) * (fp + f - fp * f), joint = fp + f - pair;
        pCloudy = prev ? joint / fp : (fp < 1 ? (f - joint) / (1 - fp) : 0);
      } else pCloudy = f;
      prev = rand() < pCloudy;
      if (prev) { const x = nu > 0 ? gammaDraw(nu) : 1; t += x * tauK[k] / f; p += x * pathK[k] / f; }
    }
    sub[s] = t; subPath[s] = p;
  }
  let clear = 0, mod = 0, ica = 0, modND = 0, icaND = 0, ln = 0;
  for (let h = 0; h < HOURS; h++) {
    const hour = 2 * Math.PI * (h + 0.5) / HOURS - Math.PI;
    const mu = Math.sin(lat) * Math.sin(dec) + Math.cos(lat) * Math.cos(dec) * Math.cos(hour);
    if (!(mu > 0)) continue;
    const beam = S0 * mu / HOURS;
    const c0 = toa(0, 0, mu, a);
    const twoCol = pathSum > 0 ? cover * toa(tauSum / cover, pathSum / cover, mu, a) + (1 - cover) * c0 : c0;
    const decked = toa(tauSum + deckTau, pathSum + deck, mu, a);
    let lnSum = 0, lnDeck = 0;
    for (let n = 0; n < GH.length; n++) { const T = meanT * Math.exp(Math.SQRT2 * sig * GH[n][0] - sig * sig / 2), wgt = GH[n][1] / Math.sqrt(Math.PI); lnSum += wgt * toa(T, T * pathPerT, mu, a); if (fraction > 0) lnDeck += wgt * toa(T + deckTau, T * pathPerT + deck, mu, a); }
    const lnCol = pathSum > 0 ? cumulative * lnSum + (1 - cumulative) * c0 : c0, lnDecked = pathSum > 0 ? cumulative * lnDeck + (1 - cumulative) * toa(deckTau, deck, mu, a) : toa(deckTau, deck, mu, a);
    ln += beam * (fraction > 0 ? (1 - fraction) * lnCol + fraction * lnDecked : lnCol);
    let icaSum = 0, icaDeck = 0;
    for (let s = 0; s < M; s++) { icaSum += toa(sub[s], subPath[s], mu, a); if (fraction > 0) icaDeck += toa(sub[s] + deckTau, subPath[s] + deck, mu, a); }
    icaSum /= M; icaDeck /= M;
    clear += beam * c0;
    modND += beam * twoCol; icaND += beam * icaSum;
    mod += beam * (fraction > 0 ? (1 - fraction) * twoCol + fraction * decked : twoCol);
    ica += beam * (fraction > 0 ? (1 - fraction) * icaSum + fraction * icaDeck : icaSum);
  }
  sums.ln += area[i] * ln; sums.clear += area[i] * clear; sums.model += area[i] * mod; sums.ica += area[i] * ica; sums.area += area[i];
  sums.modelNoDeck += area[i] * modND; sums.icaNoDeck += area[i] * icaND;
  const b = Math.min(5, Math.floor((lat * deg + 90) / 30));
  bands[b].clear += area[i] * clear; bands[b].model += area[i] * mod; bands[b].ica += area[i] * ica; bands[b].area += area[i];
  const kind = hasLow && hasUpper ? 'layered' : hasLow ? 'lowOnly' : hasUpper ? 'upperOnly' : null;
  if (kind) { kinds[kind].model += area[i] * (mod - clear); kinds[kind].ica += area[i] * (ica - clear); kinds[kind].area += area[i]; }
}
const A = sums.area, fx = (x) => x.toFixed(1);
console.log(`${FILE.split('/').pop()} M=${M} FSD=${FSD} declination ${(dec * deg).toFixed(1)}; ${((performance.now() - t0) / 1000).toFixed(0)} s`);
console.log(`global day-mean reflected (simple single-band): clear ${fx(sums.clear / A)}, model two-column ${fx(sums.model / A)}, independent columns ${fx(sums.ica / A)} W/m2`);
console.log(`moment-matched lognormal (${GH.length}-node Gauss-Hermite): cloud reflection ${fx((sums.ln - sums.clear) / A)}, ratio ${((sums.ln - sums.clear) / (sums.model - sums.clear)).toFixed(3)}`);
console.log(`cloud reflection (all-sky minus clear): model ${fx((sums.model - sums.clear) / A)}, independent ${fx((sums.ica - sums.clear) / A)}; ratio ${((sums.ica - sums.clear) / (sums.model - sums.clear)).toFixed(3)}; without deck ${fx((sums.modelNoDeck - sums.clear) / A)} vs ${fx((sums.icaNoDeck - sums.clear) / A)}`);
for (let b = 0; b < 6; b++) { const x = bands[b]; console.log(`  lat ${-90 + 30 * b}..${-60 + 30 * b}: cloud reflection model ${fx((x.model - x.clear) / x.area)} independent ${fx((x.ica - x.clear) / x.area)}`); }
for (const [name, x] of Object.entries(kinds)) console.log(`  ${name}: area ${(x.area / A).toFixed(3)}, cloud reflection model ${fx(x.model / x.area)} independent ${fx(x.ica / x.area)} (global contribution ${fx((x.ica - x.model) / A)})`);
