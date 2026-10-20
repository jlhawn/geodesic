// The low cloud of regime boxes layer by layer, as the radiation covers it,
// on the CPU engine after one step (ocean off):
//   RADIATION=… MOIST=… node scripts/tradeCover.mjs <state.bin>
// Per box (the trades of cloudRegimes.mjs, the SE Pacific, the open sea of
// 20–40S and 20–40N): the radiation's low cover and water, the mixing top,
// the plumes' share of columns and mean overshoot share, the stratiform
// ramp's share and the deck's fraction. Per layer below 3.5 km, area means
// over the box's ice-free sea: the share of columns whose layer lies below
// the mixing top, the share holding resolved cloud (q_c above 10⁻⁶), the
// grid-mean q_c, the radiated cover (radiation.layerCover, the plumes'
// cumulus included) as a grid mean; over the cloudy mixed layers the
// medians of the variance cover's spread σ_s at the defaults, the moist
// physics' uniform half-width b = a (1 − RH_c) q_s, the deficit s =
// a (q_t − q_sl) and Q₁ = s/σ_s, and the means of the variance cover
// Φ(Q₁) (f_var), the uniform distribution's sqrt(q_c/b) that made the
// water (f_unif), the Gaussian's of spread σ_s that holds the layer's q_c
// (f_held: σ_s G(Q) = q_c, G(Q) = Q Φ(Q) + φ(Q)) and the radiated cover
// (c_in); the radiated cover of the cloudy layers above the mixing top
// (c_out); σ_s with the length bounded by 0.53 w*/N in stable layers under
// a buoyant surface (σ_N, median) and the variance cover with it at the
// scales 5 and 3.2 (f_N, f_N32); the shares of columns whose plume ended
// in the layer or the layer above (top, below); and the cloudy mixed layers
// split by what mixes them (surface: below the surface-driven top or in a
// column without a cloud top; between: the decoupled column's cumulus
// layer between that top and the cloud-top layer's base; cloud-top;
// coupled), each with its share and mean radiated cover.
import { readFileSync } from 'node:fs';
import { topographyFromInt16 } from '../js/geography.module.js';
import { createModel, STATE_NAMES } from '../js/model.module.js';
import { decodeState, savedLevels } from '../js/stateFile.module.js';
import { savedDeckField, DECK_FIELDS } from '../js/physics/regrid.module.js';
import { erf } from '../js/physics/radiation.module.js';
import { cloudSaturation, criticalHumidityAt, LATENT_HEAT, MOIST_DEFAULTS } from '../js/physics/moist.module.js';
import { inLongitudes } from '../js/audit.module.js';
import { stateGrid } from '../js/stateGrid.module.js';

const FILE = process.argv[2];
const env = (name) => JSON.parse(process.env[name] ?? '{}');
const RADIATION = env('RADIATION'), MOIST = env('MOIST'), BOUNDARY_LAYER = env('BOUNDARY_LAYER'), SURFACE = env('SURFACE');
const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const model = createModel(stateGrid(saved, topography).mesh, { topography, levels: savedLevels(saved), ocean: false, radiation: RADIATION, moist: MOIST, boundaryLayer: BOUNDARY_LAYER, surface: SURFACE });
const { mesh, core, state, radiation, boundaryLayer: bl, moist, seaIce, land, geography } = model;
const { K, C, sigmaMid, g, exnerLayer, geopotential } = core.diagnostics;
const exnerLower = core.diagnostics.exnerLower ?? core.arrays.exnerLower, thetaV = core.arrays.thetaV;
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

const cp = core.diagnostics.cp ?? 1004, L = RADIATION.latentHeat ?? LATENT_HEAT, M = { ...MOIST_DEFAULTS, ...MOIST };
const scale = RADIATION.varianceScale ?? 5, floor = RADIATION.varianceFloor ?? 0.002, lambda = RADIATION.mixingLength ?? 300, stableLambda = RADIATION.stableMixingLength ?? 30;
const [pi, theta, , surfaceT, q, qc, ice] = state;
const landMask = geography.land, deg = 180 / Math.PI, area = mesh.areaCell;
const lat = Float64Array.from(mesh.latCell, (x) => x * deg), lon = Float64Array.from(mesh.lonCell, (x) => x * deg);
const fluxT = new Float64Array(C), adir = new Float64Array(C), adif = new Float64Array(C), openSea = new Float64Array(C);
const fluxState = [pi, theta, state[2], fluxT, q, qc], out = [null, new Float64Array(K * C), null, null, new Float64Array(K * C)];
for (let i = 0; i < C; i++) {
  fluxT[i] = surfaceT[i];
  if (landMask[i]) { adir[i] = adif[i] = land.albedo(i); continue; }
  const h = ice[i], cover = seaIce.cover(i, h), mu = radiation.cosZenith(i);
  adir[i] = seaIce.albedo(h, mu, seaIce.snow[i], cover); adif[i] = seaIce.albedo(h, null, seaIce.snow[i], cover);
  openSea[i] = 1 - cover;
  if (h > 0 && cover < 1) fluxT[i] = cover * surfaceT[i] + (1 - cover) * 271.35;
}
radiation.setTime(model.time);

const Phi = (x) => 0.5 * (1 + erf(x / Math.SQRT2)), phi = (x) => Math.exp(-0.5 * x * x) / Math.sqrt(2 * Math.PI);
function heldCover(ratio) {
  let lo = -9, hi = 9;
  for (let n = 0; n < 60; n++) { const m = 0.5 * (lo + hi); if (m * Phi(m) + phi(m) < ratio) lo = m; else hi = m; }
  return Phi(0.5 * (lo + hi));
}
const sat = { qs: 0, slope: 0 };
function spread(i, k, mixingDepth) {
  const bottom = (K - 1) * C + i, idx = k * C + i, ex = exnerLayer[idx], water = Math.max(0, qc[idx]), total = Math.max(0, q[idx]) + water;
  const level = theta[idx] - L * water / (cp * ex), p = pi[i] * sigmaMid[k];
  cloudSaturation(level * ex, p, M.iceSaturation, M.liquidTemperature, M.iceTemperature, sat);
  const qs = sat.qs, slope = sat.slope, a = 1 / (1 + L / cp * slope);
  let gq = 0, gl = 0, gv = 0, n = 0;
  for (const j of [k - 1, k + 1]) {
    if (j < 0 || j > K - 1 || (j < k && (geopotential[j * C + i] - geopotential[bottom]) / g >= mixingDepth)) continue;
    const jdx = j * C + i, jw = Math.max(0, qc[jdx]), dz = (geopotential[jdx] - geopotential[idx]) / g;
    gq += (Math.max(0, q[jdx]) + jw - total) / dz; gl += (theta[jdx] - L * jw / (cp * exnerLayer[jdx]) - level) / dz; gv += (thetaV[jdx] - thetaV[idx]) / dz; n++;
  }
  if (n > 1) { gq /= n; gl /= n; gv /= n; }
  const surface = geopotential[bottom] - cp * thetaV[bottom] * (exnerLower[bottom] - exnerLayer[bottom]), z = (geopotential[idx] - surface) / g;
  const length = 0.4 * z / (1 + 0.4 * z / (bl.buoyancyFlux[i] > 0 ? lambda : stableLambda));
  const b = a * (1 - criticalHumidityAt(p, pi[i], M.surfaceCriticalHumidity, M.topCriticalHumidity, M.criticalExponent)) * qs;
  const n2 = g * gv / thetaV[idx], wstar = Math.cbrt(Math.max(0, bl.buoyancyFlux[i]) * Math.max(0, bl.depth[i] - zs[i])), limited = n2 > 0 ? Math.min(length, 0.53 * wstar / Math.sqrt(n2)) : length;
  const gradient = a * Math.abs(gq - ex * slope * gl);
  return { s: a * (total - qs), sigma: Math.max(floor * qs, scale * length * gradient), sigmaN: Math.max(floor * qs, scale * limited * gradient), sigmaN32: Math.max(floor * qs, 3.2 * limited * gradient), b, z };
}
const median = (xs) => { if (!xs.length) return NaN; const s = Float64Array.from(xs).sort(); return s[s.length >> 1]; };
const levels = savedLevels(saved);
const BOXES = [
  ['N Pacific trades 15-25N 170-130W', [15, 25, -170, -130]],
  ['S Pacific trades 10-20S 160-120W', [-20, -10, -160, -120]],
  ['Atlantic trades 10-20N 50-25W', [10, 20, -50, -25]],
  ['SE Pacific 10-30S 110-80W', [-30, -10, -110, -80]],
  ['sea 20-40S', [-40, -20, -180, 180]],
  ['sea 20-40N', [20, 40, -180, 180]],
];
for (const [name, [south, north, west, east]] of BOXES) {
  const cells = [];
  for (let i = 0; i < C; i++) if (!landMask[i] && !(ice[i] > 0) && lat[i] >= south && lat[i] <= north && inLongitudes(lon[i], west, east)) cells.push(i);
  const rows = Array.from({ length: K }, () => ({ w: 0, z: 0, p: 0, mix: 0, cloudy: 0, qc: 0, cover: 0, top: 0, below: 0, sig: [], b: [], s: [], Q: [], fv: 0, fu: 0, fh: 0, fn: 0, fn32: 0, sigN: [], n: 0, cin: 0, cout: 0, nout: 0, part: [0, 0, 0, 0], partCover: [0, 0, 0, 0] }));
  let W = 0, strat = 0, deck = 0, low = 0, lowWater = 0, share = 0, plumes = 0, mixTop = 0;
  for (const i of cells) {
    radiation.apply(fluxState, out, model.surface.windSpeed, null, i, i + 1, adir, adif, null, openSea, bl.depth, 0);
    const w = area[i], bottom = (K - 1) * C + i, mixingDepth = bl.mixingTop[i] - geopotential[bottom] / g;
    W += w; strat += w * radiation.stratiform[i]; deck += w * radiation.stratusFraction[i]; low += w * radiation.lowCover[i]; lowWater += w * radiation.lowWater[i]; mixTop += w * (bl.mixingTop[i] - zs[i]);
    const topP = moist.cumulusTop[i];
    let kTop = -1;
    if (topP > 0) { let best = Infinity; for (let k = 0; k < K; k++) { const d = Math.abs(pi[i] * levels[k] - topP); if (d < best) { best = d; kTop = k; } } plumes += w; share += w * (moist.cumulusOvershootShare ? moist.cumulusOvershootShare[i] : 1); }
    for (let k = K - 1; k >= 0; k--) {
      const idx = k * C + i, r = rows[k], d = spread(i, k, mixingDepth);
      if (d.z > 3500) break;
      r.w += w; r.z += w * d.z; r.p += w * pi[i] * sigmaMid[k] / 100;
      const inside = (geopotential[idx] - geopotential[bottom]) / g < mixingDepth, water = Math.max(0, qc[idx]);
      if (inside) r.mix += w;
      if (k === kTop) r.top += w;
      if (k === kTop + 1) r.below += w;
      r.qc += w * water * 1e3; r.cover += w * (water > 0 || moist.cumulusCover[idx] > 0 ? radiation.layerCover[k] : 0);
      if (water > 1e-6) r.cloudy += w;
      if (water > 1e-6 && !inside) { r.cout += w * radiation.layerCover[k]; r.nout += w; }
      if (water > 1e-6 && inside) {
        r.sig.push(d.sigma * 1e3); r.b.push(d.b * 1e3); r.s.push(d.s * 1e3); r.Q.push(d.s / d.sigma);
        const hs = bl.depth[i] - zs[i], zb = bl.decoupling ? bl.decoupling[i] - zs[i] : 0, regime = bl.regime[i];
        const part = regime === 3 ? 3 : regime === 2 ? (d.z < hs ? 0 : d.z < zb ? 1 : 2) : 0;
        r.part[part] += w; r.partCover[part] += w * radiation.layerCover[k];
        r.fv += w * Phi(d.s / d.sigma); r.fn += w * Phi(d.s / d.sigmaN); r.fn32 += w * Phi(d.s / d.sigmaN32); r.sigN.push(d.sigmaN * 1e3); r.fu += w * Math.min(1, Math.sqrt(water / d.b)); r.fh += w * heldCover(water / d.sigma); r.cin += w * radiation.layerCover[k]; r.n += w;
      }
    }
  }
  console.log(`\n${name}: ${cells.length} cells, low cover ${(low / W).toFixed(3)}, low water ${(1e3 * lowWater / W).toFixed(1)} g/m², mixing top ${(mixTop / W).toFixed(0)} m, plumes on ${(plumes / W).toFixed(2)}, overshoot share ${(share / Math.max(plumes, 1e-30)).toFixed(2)}, stratiform share ${(strat / W).toFixed(2)}, deck ${(deck / W).toFixed(3)}`);
  console.log('   z m    hPa   mixed cloudy  qc g/kg cover |  σ_s    b     s    Q₁  (g/kg, medians) | f_var f_unif f_held  c_in | c_out | σ_N f_N f_N32 | top below | cloudy mixed by part: surface / between / cloud-top / coupled (share, cover)');
  for (let k = K - 1; k >= 0; k--) {
    const r = rows[k];
    if (!(r.w > 0)) continue;
    const f = (x, n = 3) => (Number.isFinite(x) ? x.toFixed(n) : '  -  ').padStart(6);
    console.log(`${(r.z / r.w).toFixed(0).padStart(6)} ${(r.p / r.w).toFixed(0).padStart(6)} ${f(r.mix / r.w, 2)} ${f(r.cloudy / r.w, 2)} ${f(r.qc / r.w, 4)} ${f(r.cover / r.w)} | ${f(median(r.sig), 2)} ${f(median(r.b), 2)} ${f(median(r.s), 2)} ${f(median(r.Q), 2)} | ${f(r.fv / r.n)} ${f(r.fu / r.n)} ${f(r.fh / r.n)} ${f(r.cin / r.n)} | ${f(r.cout / r.nout)} | ${f(median(r.sigN), 2)} ${f(r.fn / r.n)} ${f(r.fn32 / r.n)} | ${f(r.top / r.w, 2)} ${f(r.below / r.w, 2)} | ${r.part.map((x, j) => `${(x / Math.max(r.n, 1e-30)).toFixed(2)} ${(r.partCover[j] / Math.max(x, 1e-30)).toFixed(2)}`).join(' / ')}`);
  }
}
