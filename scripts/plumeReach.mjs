// The shallow plume's speed and reach in regime boxes, on the CPU engine
// after one step (ocean off), and the same ascent through BOMEX's initial
// profile (Siebesma et al. 2003, Table B1):
//   MOIST=… node scripts/plumeReach.mjs <state.bin>
// The ascent is cumulusColumn's (moist.module.js), checked column by
// column against its top layer and overshoot share. Per box, over the
// ice-free sea columns where the shallow plume runs: the number whose
// first cloudy layer is already not buoyant, and medians [quartiles] of
// the first cloudy layer's base, the number of buoyant cloudy layers, the
// LCL height, w* = (B₀ h)^⅓, z_b the base of the first non-buoyant cloudy
// layer, that layer's depth and the plume's virtual temperature deficit
// there, the interface of strongest θ_l stability between the LCL and
// 4 km (the inversion) and its gradient; then for each w equation dw²/dz =
// 2 a B − 2 b ε w² from w0 at the first cloudy layer, the speed at z_b,
// the overshoot share, and the implied top z_b + share × depth less the
// inversion's height; the same ascent with ε = 1/z; and the median speed
// by height; and the model's own plume: its implied top (the top layer's
// base + share × depth), that less the inversion, and its share. The check
// and the first set follow MOIST's cumulusEntrainmentHeight and
// cumulusAscent; the other sets and z_b assume the buoyant ascent. For
// BOMEX: the speed every 100 m and the top where w² reaches zero for each
// equation on 20 m layers, and the share and implied top on the state's
// levels.
import { readFileSync } from 'node:fs';
import { topographyFromInt16 } from '../js/geography.module.js';
import { createModel, STATE_NAMES } from '../js/model.module.js';
import { decodeState, savedLevels } from '../js/stateFile.module.js';
import { savedDeckField, DECK_FIELDS } from '../js/physics/regrid.module.js';
import { liftingCondensationLevel, MOIST_DEFAULTS } from '../js/physics/moist.module.js';
import { VIRTUAL_FACTOR } from '../js/dynamics/sigmaCore.module.js';
import { inLongitudes } from '../js/audit.module.js';
import { stateGrid } from '../js/stateGrid.module.js';

const FILE = process.argv[2];
const env = (name) => JSON.parse(process.env[name] ?? '{}');
const RADIATION = env('RADIATION'), MOIST = env('MOIST'), BOUNDARY_LAYER = env('BOUNDARY_LAYER'), SURFACE = env('SURFACE');
const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const model = createModel(stateGrid(saved, topography).mesh, { topography, levels: savedLevels(saved), ocean: false, radiation: RADIATION, moist: MOIST, boundaryLayer: BOUNDARY_LAYER, surface: SURFACE });
const { mesh, core, state, radiation, boundaryLayer: bl, moist, seaIce, land, geography } = model;
STATE_NAMES.forEach((name, a) => state[a].set(saved[name]));
seaIce.load(state[6], saved.concentration ?? null);
for (const field of Object.keys(DECK_FIELDS)) radiation[field].set(savedDeckField(saved, field, model));
land.load({ soil: Float64Array.from(saved.land.soil), snow: Float64Array.from(saved.land.snow), ...(saved.land.vegetation ? { vegetation: Float64Array.from(saved.land.vegetation) } : {}), ...(saved.land.surface ? { surface: Float64Array.from(saved.land.surface) } : {}) });
model.time = saved.time;
const { K, C, levels, dSigma, sigmaMid, cp, R, g, kappa, exnerLayer, geopotential } = core.diagnostics;
const exnerLower = core.diagnostics.exnerLower ?? core.arrays.exnerLower, thetaV = core.arrays.thetaV;
const zs = Float64Array.from({ length: C }, (_, i) => (model.surfaceGeopotential ? model.surfaceGeopotential[i] / g : 0));
if (saved.boundaryDepth) bl.depth.set(Float64Array.from(saved.boundaryDepth, (z, i) => z + zs[i]));
if (saved.mixingTop) bl.mixingTop.set(Float64Array.from(saved.mixingTop, (z, i) => z + zs[i]));
if (saved.boundaryRegime) bl.regime.set(saved.boundaryRegime);
if (saved.boundaryBuoyancy) bl.buoyancyFlux.set(saved.boundaryBuoyancy);
const dt = 1350 * 16 / saved.N;
model.step(dt);

const M = { ...MOIST_DEFAULTS, ...MOIST }, Lv = M.latentHeat;
const virtual = M.virtualBuoyancy ? VIRTUAL_FACTOR : 0, loading = M.virtualBuoyancy ? 1 : 0;
const entrainmentRate = M.cumulusEntrainmentHeight != null ? (z) => M.cumulusEntrainmentHeight / Math.max(z, 1) : () => M.cumulusEntrainment, carried = M.cumulusAscent === 'kinetic';
const [pi, theta, , , q, qc, ice] = state;

const kinetic = (w2, work, depth, mixing, a) => { const x = 2 * mixing * depth; return w2 * Math.exp(-x) + 2 * a * work * (x > 0 ? -Math.expm1(-x) / x : 1); };
function overshoot(w2, work, depth, mixing, a) {
  if (!(work < 0)) return 1;
  if (!(w2 > 0)) return 0;
  const lift = -a * work / depth, y = mixing * w2 / lift;
  const reach = mixing > 0 ? (y > 1e-3 ? Math.log1p(y) : y * (1 - 0.5 * y)) / (2 * mixing) : w2 / (2 * lift);
  return Math.min(1, reach / depth);
}

// col: arrays over layers k = 0 (top) … n − 1 (bottom): T, p, dp, z (midpoint heights), S, Q, air, cloud,
// upper (height of layer k's upper interface), lowerP (pressure of its lower interface); entrainment(z) per meter.
function ascend(col, source, sourceS, sourceQ, lcl, entrainment, through = false, carry = false) {
  const n = col.T.length, layers = [];
  let s = sourceS, w = sourceQ, cloudy = false, guess = 0, top = -1, capped = false, base = NaN;
  for (let k = source - 1; k >= 0; k--) {
    if (!(col.p[k] > M.shallowTop) && !through) { if (cloudy) { top = k + 1; capped = true; } break; }
    const below = col.upper[k + 1], above = col.upper[k];
    const mixes = col.lowerP[k] <= lcl.pressure, epsilon = mixes ? entrainment(col.z[k] - col.surface) : 0;
    const half = Math.exp(-epsilon * (col.z[k] - below));
    const midS = col.S[k] + (s - col.S[k]) * half, midQ = col.Q[k] + (w - col.Q[k]) * half;
    const plume = moist.plumeAir(midS, midQ, col.z[k], col.p[k], guess);
    guess = plume.T;
    const work = R * (plume.T * (1 + virtual * (midQ - plume.liquid) - loading * plume.liquid) - col.T[k] * (1 + virtual * col.air[k] - loading * col.cloud[k])) * col.dp[k] / col.p[k];
    if (plume.liquid > 0 && !cloudy) { cloudy = true; base = below; }
    if (cloudy) layers.push({ k, work, depth: above - below, epsilon, below, above, deficit: work / (above - below) * col.T[k] / g });
    if (cloudy && !(work > 0) && !through && !carry) { top = k; break; }
    const full = Math.exp(-epsilon * (above - below));
    s = col.S[k] + (s - col.S[k]) * full; w = col.Q[k] + (w - col.Q[k]) * full;
  }
  return { top, capped, base, layers };
}

// Speed at each layer top, and the reach: within the first non-buoyant layer (the model's) or, with
// `through`, the height where w² reaches zero.
function reach(layers, w0, a, b, through = false) {
  let w2 = w0 * w0, speeds = [];
  for (const L of layers) {
    if (!(L.work > 0)) {
      const share = overshoot(w2, L.work, L.depth, b * L.epsilon, a);
      if (!through || share < 1) return { atTop: Math.sqrt(w2), share, top: L.below + share * L.depth, speeds, k: L.k };
    }
    w2 = kinetic(w2, L.work, L.depth, b * L.epsilon, a);
    if (!(w2 > 0)) return { atTop: 0, share: 0, top: L.above, speeds };
    speeds.push([L.above, Math.sqrt(w2)]);
  }
  return { atTop: Math.sqrt(w2), share: 1, top: NaN, speeds };
}

const SETS = [
  ['model: w0 = plumeVelocity, a, b = plumeAcceleration, plumeDrag', null],
  ['Gregory 2001 (IFS): a 1/3, b 2, w0 = plumeVelocity', { a: 1 / 3, b: 2 }],
  ['Gregory 2001: a 1/3, b 2, w0 = w*', { a: 1 / 3, b: 2, star: true }],
  ['Simpson-Wiggert 1969: a 2/3, b 2, w0 = w*', { a: 2 / 3, b: 2, star: true }],
];
const coefficients = (set, star) => (set === null ? { w0: M.plumeVelocity, a: M.plumeAcceleration, b: M.plumeDrag } : { w0: set.star ? star : M.plumeVelocity, a: set.a, b: set.b });

function modelColumn(i) {
  const n = K, col = { T: new Float64Array(n), p: new Float64Array(n), dp: new Float64Array(n), z: new Float64Array(n), S: new Float64Array(n), Q: new Float64Array(n), air: new Float64Array(n), cloud: new Float64Array(n), upper: new Float64Array(n + 1), lowerP: new Float64Array(n), level: new Float64Array(n) };
  for (let k = 0; k < n; k++) {
    const idx = k * C + i, cloud = Math.max(0, qc[idx]);
    col.T[k] = theta[idx] * exnerLayer[idx]; col.p[k] = pi[i] * sigmaMid[k]; col.dp[k] = pi[i] * dSigma[k]; col.z[k] = geopotential[idx] / g;
    col.S[k] = cp * col.T[k] + g * col.z[k] - Lv * cloud; col.Q[k] = Math.max(0, q[idx]) + cloud; col.air[k] = Math.max(0, q[idx]); col.cloud[k] = cloud;
    col.upper[k] = k > 0 ? (geopotential[idx] + cp * thetaV[idx] * (exnerLayer[idx] - exnerLower[(k - 1) * C + i])) / g : Infinity;
    col.lowerP[k] = pi[i] * levels[k + 1];
    col.level[k] = theta[idx] - Lv * cloud / (cp * exnerLayer[idx]);
  }
  const bottom = n - 1;
  col.surface = (geopotential[bottom * C + i] - cp * thetaV[bottom * C + i] * (exnerLower[bottom * C + i] - exnerLayer[bottom * C + i])) / g;
  col.upper[n] = col.surface;
  return col;
}

function lclHeight(col, pressure) {
  const n = col.T.length;
  for (let k = n - 2; k >= 0; k--) if (col.p[k] <= pressure) { const f = Math.log(col.p[k + 1] / pressure) / Math.log(col.p[k + 1] / col.p[k]); return col.z[k + 1] + f * (col.z[k] - col.z[k + 1]) - col.surface; }
  return NaN;
}

const quart = (xs) => { const s = Float64Array.from(xs.filter(Number.isFinite)).sort(); if (!s.length) return [NaN, NaN, NaN]; return [s[Math.floor(0.25 * (s.length - 1))], s[Math.floor(0.5 * (s.length - 1))], s[Math.floor(0.75 * (s.length - 1))]]; };
const show = (xs, digits = 0) => { const [lo, mid, hi] = quart(xs); return `${mid.toFixed(digits)} [${lo.toFixed(digits)}, ${hi.toFixed(digits)}]`; };

const landMask = geography.land, deg = 180 / Math.PI;
const BOXES = [
  ['N Pacific trades 15-25N 170-130W', [15, 25, -170, -130]],
  ['S Pacific trades 10-20S 160-120W', [-20, -10, -160, -120]],
  ['Atlantic trades 10-20N 50-25W', [10, 20, -50, -25]],
  ['SE Pacific 10-30S 110-80W', [-30, -10, -110, -80]],
  ['sea 20-40S', [-40, -20, -180, 180]],
  ['sea 20-40N', [20, 40, -180, 180]],
];
let checked = 0, mismatched = 0, worst = 0;
for (const [name, [south, north, west, east]] of BOXES) {
  const rows = { n: 0, cells: 0, capped: 0, forced: 0, buoyant: [], base: [], lcl: [], star: [], zb: [], depth: [], deficit: [], inversion: [], gradient: [], wmax: [], topz: [], sets: SETS.map(() => ({ atTop: [], share: [], over: [] })), model: { top: [], over: [], share: [] }, alt: { n: 0, forced: 0, capped: 0, share: [], over: [], top: [] }, profile: new Map() };
  for (let i = 0; i < C; i++) {
    const lat = mesh.latCell[i] * deg, lon = mesh.lonCell[i] * deg;
    if (landMask[i] || ice[i] > 0 || lat < south || lat > north || !inLongitudes(lon, west, east)) continue;
    rows.cells++;
    const saveTheta = new Float64Array(K), saveQ = new Float64Array(K);
    for (let k = 0; k < K; k++) { saveTheta[k] = theta[k * C + i]; saveQ[k] = q[k * C + i]; }
    moist.cumulusColumn(i, pi, theta, q, qc, dt);
    const topModel = moist.cumulus.top, shareModel = moist.cumulus.overshoot, ran = moist.cumulus.baseFlux > 0;
    for (let k = 0; k < K; k++) { theta[k * C + i] = saveTheta[k]; q[k * C + i] = saveQ[k]; }
    if (!ran || topModel < 0) continue;
    const col = modelColumn(i), bottom = K - 1, depth = bl.depth[i];
    let mass = 0, energy = 0, water = 0, source = bottom;
    for (let k = bottom; k >= 0; k--) {
      if (k < bottom && ((!(col.upper[k + 1] < depth) && !(col.p[k] >= pi[i] - M.cumulusSourceDepth)) || !(col.p[k] > M.shallowTop))) break;
      mass += col.dp[k]; energy += col.dp[k] * col.S[k]; water += col.dp[k] * col.Q[k]; source = k;
    }
    const sourceS = M.cumulusSource === 'lowest' ? col.S[bottom] : energy / mass, sourceQ = M.cumulusSource === 'lowest' ? col.Q[bottom] : water / mass;
    const lcl = liftingCondensationLevel((sourceS - g * col.z[bottom]) / cp, sourceQ, col.p[bottom], kappa);
    const run = ascend(col, source, sourceS, sourceQ, lcl, entrainmentRate, false, carried);
    const mine = reach(run.layers, M.plumeVelocity, M.plumeAcceleration, M.plumeDrag, carried);
    checked++;
    const stopped = carried && mine.k !== undefined, myTop = stopped ? mine.k : run.top, myShare = carried ? (stopped ? mine.share : 1) : (run.capped ? 1 : mine.share);
    const miss = myTop !== topModel ? 1 : Math.abs(myShare - shareModel);
    if (miss > 1e-9) mismatched++;
    worst = Math.max(worst, miss);
    rows.n++;
    const modelTop = col.upper[topModel + 1] + shareModel * (col.upper[topModel] - col.upper[topModel + 1]) - col.surface;
    rows.model.top.push(modelTop); rows.model.share.push(shareModel);
    if (run.capped && !stopped) { rows.capped++; continue; }
    const star = Math.cbrt(Math.max(0, bl.buoyancyFlux[i]) * Math.max(0, depth - col.z[bottom]));
    const first = run.layers.findIndex((L) => !(L.work > 0)), buoyant = first < 0 ? run.layers.length : first;
    const zLcl = lclHeight(col, lcl.pressure), last = run.layers[first < 0 ? run.layers.length - 1 : first];
    if (!buoyant) rows.forced++; rows.buoyant.push(buoyant); rows.base.push(run.base - col.surface);
    rows.lcl.push(zLcl); rows.star.push(star); rows.zb.push(last.below - col.surface); rows.depth.push(last.depth); rows.deficit.push(last.deficit);
    let best = -Infinity, at = NaN;
    for (let k = K - 2; k >= 0; k--) {
      if (col.z[k + 1] - col.surface < zLcl) continue;
      if (col.z[k] - col.surface > 4000) break;
      const gradient = (col.level[k] - col.level[k + 1]) / (col.z[k] - col.z[k + 1]);
      if (gradient > best) { best = gradient; at = col.upper[k + 1] - col.surface; }
    }
    rows.inversion.push(at); rows.gradient.push(1e3 * best); rows.model.over.push(modelTop - at);
    SETS.forEach(([, set], n) => {
      const { w0, a, b } = coefficients(set, star), r = reach(run.layers, w0, a, b);
      rows.sets[n].atTop.push(r.atTop); rows.sets[n].share.push(r.share); rows.sets[n].over.push(r.top - col.surface - at);
      if (n === 0) {
        rows.topz.push(r.top - col.surface);
        rows.wmax.push(Math.max(w0, ...r.speeds.map(([, v]) => v)));
        for (const [z, v] of r.speeds) { const bin = Math.round((z - col.surface) / 250) * 250; if (!rows.profile.has(bin)) rows.profile.set(bin, []); rows.profile.get(bin).push(v); }
      }
    });
    const alt = ascend(col, source, sourceS, sourceQ, lcl, (z) => 1 / Math.max(z, 1));
    rows.alt.n++;
    if (alt.capped || !alt.layers.length) { rows.alt.capped++; continue; }
    const r = reach(alt.layers, M.plumeVelocity, M.plumeAcceleration, M.plumeDrag);
    if (!alt.layers.some((L) => L.work > 0)) rows.alt.forced++;
    rows.alt.share.push(r.share); rows.alt.top.push(r.top - col.surface); rows.alt.over.push(r.top - col.surface - at);
  }
  console.log(`\n${name}: ${rows.cells} sea cells, shallow plume on ${rows.n}, ${rows.capped} capped at shallowTop, ${rows.forced} (${(rows.forced / Math.max(1, rows.n)).toFixed(2)}) with no buoyant cloudy layer; first cloudy layer base ${show(rows.base)} m, buoyant cloudy layers ${show(rows.buoyant)}`);
  console.log(`  LCL ${show(rows.lcl)} m, w* ${show(rows.star, 2)} m/s, z_b (base of the first non-buoyant cloudy layer) ${show(rows.zb)} m, next layer depth ${show(rows.depth)} m, plume Tv deficit there ${show(rows.deficit.map((x) => -x), 2)} K`);
  console.log(`  strongest θ_l interface above the LCL ${show(rows.inversion)} m at ${show(rows.gradient, 1)} K/km; model's largest speed ${show(rows.wmax, 2)} m/s`);
  SETS.forEach(([label], n) => console.log(`  ${label}: w at z_b ${show(rows.sets[n].atTop, 2)} m/s, share ${show(rows.sets[n].share, 2)}, implied top less inversion ${show(rows.sets[n].over)} m${n === 0 ? `, implied top ${show(rows.topz)} m` : ''}`));
  console.log(`  the ascent with ε = 1/z (Siebesma et al. 2003), the model's w equation: ${rows.alt.forced} of ${rows.alt.n} with no buoyant cloudy layer, ${rows.alt.capped} capped at shallowTop; share ${show(rows.alt.share, 2)}, implied top ${show(rows.alt.top)} m, less inversion ${show(rows.alt.over)} m`);
  console.log(`  the model's plume: implied top ${show(rows.model.top)} m, less inversion ${show(rows.model.over)} m, share ${show(rows.model.share, 2)}`);
  console.log(`  model's speed by height (m: median m/s, columns): ${[...rows.profile.keys()].sort((x, y) => x - y).filter((z) => z <= 3000).map((z) => `${z}: ${quart(rows.profile.get(z))[1].toFixed(2)} (${rows.profile.get(z).length})`).join(', ')}`);
}
console.log(`\nascent check: ${checked} columns, ${mismatched} differ from cumulusColumn's top or share (largest ${worst.toExponential(1)})`);

// BOMEX, Siebesma et al. (2003) Table B1, surface 1015 hPa, no cloud water in the environment.
const PROFILE = { z: [0, 520, 1480, 2000, 3000, 4000], thetaL: [298.7, 298.7, 302.4, 308.2, 311.85, 315.5], qt: [17.0, 16.3, 10.7, 4.2, 3.0, 2.5] };
const interp = (xs, ys, x) => { for (let j = 1; j < xs.length; j++) if (x <= xs[j]) return ys[j - 1] + (ys[j] - ys[j - 1]) * (x - xs[j - 1]) / (xs[j] - xs[j - 1]); return ys[ys.length - 1]; };
const fine = 20, top = 4000, heights = [], pressures = [];
{
  let p = 101500;
  for (let z = 0; z <= top; z += fine) {
    heights.push(z); pressures.push(p);
    const zm = z + fine / 2, qt = 1e-3 * interp(PROFILE.z, PROFILE.qt, zm);
    let mid = p;
    for (let n = 0; n < 4; n++) { const tv = interp(PROFILE.z, PROFILE.thetaL, zm) * Math.pow(mid / 1e5, kappa) * (1 + VIRTUAL_FACTOR * qt); mid = p * Math.exp(-g * fine / (2 * R * tv)); }
    p = mid * mid / p;
  }
}
const heightAt = (pr) => { for (let j = 1; j < pressures.length; j++) if (pr >= pressures[j]) { const f = Math.log(pressures[j - 1] / pr) / Math.log(pressures[j - 1] / pressures[j]); return heights[j - 1] + f * fine; } return top; };
function bomexColumn(interfaces) {
  const n = interfaces.length - 1, col = { T: new Float64Array(n), p: new Float64Array(n), dp: new Float64Array(n), z: new Float64Array(n), S: new Float64Array(n), Q: new Float64Array(n), air: new Float64Array(n), cloud: new Float64Array(n), upper: new Float64Array(n + 1), lowerP: new Float64Array(n), surface: 0 };
  for (let k = 0; k < n; k++) {
    const pTop = interfaces[k], pBottom = interfaces[k + 1], pm = Math.sqrt(pTop * pBottom), zm = heightAt(pm), qt = 1e-3 * interp(PROFILE.z, PROFILE.qt, zm);
    col.p[k] = pm; col.dp[k] = pBottom - pTop; col.z[k] = zm; col.T[k] = interp(PROFILE.z, PROFILE.thetaL, zm) * Math.pow(pm / 1e5, kappa);
    col.S[k] = cp * col.T[k] + g * zm; col.Q[k] = qt; col.air[k] = qt; col.upper[k] = heightAt(pTop); col.lowerP[k] = pBottom;
  }
  col.upper[n] = 0;
  return col;
}
// `air`: the source's θ_l (K) and q_t (g/kg) in place of the mean of the layers below 500 m.
function bomexRun(col, entrainment, through, air = null) {
  const n = col.T.length, bottom = n - 1;
  let mass = 0, energy = 0, water = 0, source = bottom;
  for (let k = bottom; k >= 0 && col.upper[k + 1] < 500; k--) { mass += col.dp[k]; energy += col.dp[k] * col.S[k]; water += col.dp[k] * col.Q[k]; source = k; }
  const sourceS = air ? cp * air.thetaL * Math.pow(col.p[bottom] / 1e5, kappa) + g * col.z[bottom] : energy / mass, sourceQ = air ? 1e-3 * air.qt : water / mass;
  const lcl = liftingCondensationLevel((sourceS - g * col.z[bottom]) / cp, sourceQ, col.p[bottom], kappa);
  return { lcl: heightAt(lcl.pressure), run: ascend(col, source, sourceS, sourceQ, lcl, entrainment, through) };
}
const CORE_AIR = { thetaL: 298.8, qt: 17.25 }, perHeight = (z) => 1 / Math.max(z, 1);
const BOMEX_SETS = [
  ['model: w0 1, a 1/3, b 1, ε 2.5e-3', 1, 1 / 3, 1, () => M.cumulusEntrainment],
  ['model, w0 0.6 (the LES at cloud base)', 0.6, 1 / 3, 1, () => M.cumulusEntrainment],
  ['Gregory 2001: w0 0.6, a 1/3, b 2, ε 2.5e-3', 0.6, 1 / 3, 2, () => M.cumulusEntrainment],
  ['Simpson-Wiggert: w0 0.6, a 2/3, b 2, ε 2.5e-3', 0.6, 2 / 3, 2, () => M.cumulusEntrainment],
  ['Siebesma 2003 eq. 15: w0 0.6, a 1/3, b 2, ε 1/z', 0.6, 1 / 3, 2, perHeight],
  ['model, ε 1/z', 1, 1 / 3, 1, perHeight],
  ['model, source the LES core at cloud base (θ_l 298.8 K, q_t 17.25 g/kg)', 1, 1 / 3, 1, () => M.cumulusEntrainment, CORE_AIR],
  ['model, ε 1/z and the core source', 1, 1 / 3, 1, perHeight, CORE_AIR],
  ['Siebesma 2003 eq. 15 with the core source', 0.6, 1 / 3, 2, perHeight, CORE_AIR],
];
const fineInterfaces = []; for (let z = top; z >= 0; z -= fine) fineInterfaces.push(pressures[heights.indexOf(z)]);
const stateInterfaces = Array.from(levels, (s) => s * 101500);
console.log('\nBOMEX (Siebesma et al. 2003, Table B1): w every 100 m from 600 m, the top where w² = 0 (20 m layers), and on the state\'s levels the share and implied top');
for (const [label, w0, a, b, entrainment, air] of BOMEX_SETS) {
  const { lcl, run } = bomexRun(bomexColumn(fineInterfaces), entrainment, true, air), r = reach(run.layers, w0, a, b, true);
  const speeds = []; for (let z = 600; z <= 2000; z += 100) { const s = r.speeds.find(([h]) => h >= z); speeds.push(s && s[0] < r.top ? s[1].toFixed(2) : '0'); }
  const coarse = bomexRun(bomexColumn(stateInterfaces), entrainment, false, air), cr = reach(coarse.run.layers, w0, a, b), last = coarse.run.layers[coarse.run.layers.length - 1];
  console.log(`  ${label}: LCL ${lcl.toFixed(0)} m; w ${speeds.join(' ')}; top ${r.top.toFixed(0)} m | levels: first non-buoyant layer ${last.below.toFixed(0)}-${last.above.toFixed(0)} m, w at its base ${cr.atTop.toFixed(2)}, share ${cr.share.toFixed(2)}, implied top ${cr.top.toFixed(0)} m`);
}
process.exit(0);
