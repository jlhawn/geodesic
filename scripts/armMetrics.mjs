// Time-mean equatorial Pacific metrics of an arm over several saved states,
// on the CPU with the ocean off:
//   MOIST='{...}' node scripts/armMetrics.mjs <label> <state.bin>...
import { readFileSync } from 'node:fs';
const R = new URL('../', import.meta.url).pathname;
const [{ stateGeography }, { createModel, STATE_NAMES }, { decodeState, savedLevels }, { savedDeckField, DECK_FIELDS }, { stateGrid }, { saturationHumidity }] = await Promise.all([
  import(R + 'js/geography.module.js'), import(R + 'js/model.module.js'), import(R + 'js/stateFile.module.js'), import(R + 'js/physics/regrid.module.js'), import(R + 'js/stateGrid.module.js'), import(R + 'js/physics/moist.module.js')]);
const [LABEL, ...FILES] = process.argv.slice(2), LAT = Number(process.env.LAT ?? 3);
const env = (n) => JSON.parse(process.env[n] ?? '{}');
const first = await decodeState(new Uint8Array(readFileSync(FILES[0])));
const { topography, geography, land: landOptions } = stateGeography(first);
const model = createModel(stateGrid(first, topography).mesh, { topography, geography, land: landOptions, levels: savedLevels(first), ocean: false, radiation: env('RADIATION'), moist: env('MOIST'), boundaryLayer: env('BOUNDARY_LAYER'), surface: env('SURFACE') });
const { mesh, core, state, radiation, seaIce, land, surface } = model;
const { K, C, E, sigmaMid, R: Rd, g, exnerLayer, geopotential } = core.diagnostics;
const { thetaV } = core.arrays;
const deg = 180 / Math.PI, { latCell, lonCell, areaCell, dcEdge, dvEdge } = mesh, landMask = model.geography.land;
const BINS = [[140, 160], [160, 180], [180, 200], [200, 220], [220, 240], [240, 260], [260, 280]], NB = BINS.length, REF = 5;
const binName = ([a, b]) => `${a <= 180 ? a + 'E' : 360 - a + 'W'}-${b <= 180 ? b + 'E' : 360 - b + 'W'}`.replace('180E', '180').replace('180W', '180');
const binOf = new Int8Array(C).fill(-1);
for (let i = 0; i < C; i++) {
  if (landMask[i] || first.ice[i] > 0 || Math.abs(latCell[i] * deg) > LAT) continue;
  const lon = ((lonCell[i] * deg) % 360 + 360) % 360;
  binOf[i] = BINS.findIndex(([a, b]) => lon >= a && lon < b);
}
const cells = []; for (let i = 0; i < C; i++) if (binOf[i] >= 0) cells.push(i);
const binMean = (fn) => { const s = new Float64Array(NB), a = new Float64Array(NB); for (const i of cells) { const v = fn(i); if (!Number.isFinite(v)) continue; s[binOf[i]] += areaCell[i] * v; a[binOf[i]] += areaCell[i]; } return Array.from(s, (x, b) => x / a[b]); };
const vec = new Float64Array(3 * C);
const east = (i) => -Math.sin(lonCell[i]) * vec[3 * i] + Math.cos(lonCell[i]) * vec[3 * i + 1];
function eastAt(edgeField, out) {
  for (const i of cells) {
    let x = 0, y = 0, z = 0;
    for (let m = 0; m < mesh.nEdgesOnCell[i]; m++) {
      const e = mesh.edgesOnCell[mesh.maxEdges * i + m], w = 0.5 * dcEdge[e] * dvEdge[e] * edgeField[e];
      x += w * mesh.nEdge[3 * e]; y += w * mesh.nEdge[3 * e + 1]; z += w * mesh.nEdge[3 * e + 2];
    }
    vec[3 * i] = x / areaCell[i]; vec[3 * i + 1] = y / areaCell[i]; vec[3 * i + 2] = z / areaCell[i];
    out[i] = east(i);
  }
  return out;
}
const HEIGHTS = [500, 1000, 1500, 2000, 3000, 4000, 6000, 9000, 12000], BANDS = [0, 1000, 2000, 3000, 5000, 8000, 12000, 18000], LAYERS = [35, 30, 26, 25, 24, 23, 22, 21, 19];
const sums = {}, add = (name, vals) => { if (!sums[name]) sums[name] = new Float64Array(NB); vals.forEach((v, b) => sums[name][b] += v); };
const zs = Float64Array.from({ length: C }, (_, i) => (model.surfaceGeopotential ? model.surfaceGeopotential[i] / g : 0));
const zMid = new Float64Array(K * C), uE = new Float64Array(C), tauE = new Float64Array(C);
const atHeight = (arr, i, zq) => { for (let k = K - 1; k > 0; k--) { const za = zMid[(k - 1) * C + i] - zs[i], zb = zMid[k * C + i] - zs[i]; if (zq <= za) { const t = Math.max(0, Math.min(1, (zq - zb) / (za - zb))); return arr[k * C + i] + t * (arr[(k - 1) * C + i] - arr[k * C + i]); } } return arr[i]; };
let files = 0;
for (const file of FILES) {
  const saved = file === FILES[0] ? first : await decodeState(new Uint8Array(readFileSync(file)));
  STATE_NAMES.forEach((name, a) => state[a].set(saved[name]));
  seaIce.load(state[6], saved.concentration ?? null);
  for (const field of Object.keys(DECK_FIELDS)) radiation[field].set(savedDeckField(saved, field, model));
  land.load(saved.land, state[6]);
  model.time = saved.time;
  const [pi, theta, u, , q] = state;
  core.phaseFlux(state, 0, K); core.phaseColumn(state, [new Float64Array(C)], 0, C);
  surface.lowestWindSpeed(u); surface.applyLayers(state, [null, null, new Float64Array(K * E)], 0, K);
  const stressEdge = surface.stress(state, new Float64Array(E)), speed = surface.windSpeed;
  const geo = Float64Array.from(geopotential), exL = Float64Array.from(exnerLayer), thv = Float64Array.from(thetaV);
  for (let k = 0; k < K; k++) for (const i of cells) zMid[k * C + i] = geo[k * C + i] / g;
  eastAt(u.subarray((K - 1) * E, K * E), uE); eastAt(stressEdge, tauE);
  const tv = new Float64Array(K * C); for (let k = 0; k < K; k++) for (const i of cells) tv[k * C + i] = thv[k * C + i] * exL[k * C + i];
  const T = (k, i) => theta[k * C + i] * exL[k * C + i];
  add('sea-level pressure (hPa)', binMean((i) => pi[i] * Math.exp(g * zs[i] / (Rd * T(K - 1, i))) / 100));
  add('u lowest layer (m/s)', binMean((i) => uE[i]));
  add('|V| lowest layer (m/s)', binMean((i) => speed[i]));
  add('stress tau_x (mN/m2)', binMean((i) => 1000 * tauE[i]));
  for (const k of LAYERS) { add(`T k=${k} (C)`, binMean((i) => T(k, i) - 273.15)); add(`z k=${k} (m)`, binMean((i) => zMid[k * C + i] - zs[i])); }
  for (const k of [26, 24, 23]) add(`RH k=${k}`, binMean((i) => q[k * C + i] / saturationHumidity(T(k, i), pi[i] * sigmaMid[k])));
  if (saved.mlmHeight) add('deck height mlmHeight (m)', binMean((i) => saved.mlmHeight[i]));
  if (saved.mlmSubsidence) add('deck subsidence mlmSubsidence (mm/s)', binMean((i) => 1000 * saved.mlmSubsidence[i]));
  if (saved.mlmGate) add('deck gate mlmGate', binMean((i) => saved.mlmGate[i]));
  if (saved.cumulusCover) { const cc = saved.cumulusCover, off = K * C - cc.length; add('cumulus top, highest layer with cover > 0.01 (m)', binMean((i) => { for (let k = 0; k < K; k++) { const j = k * C + i - off; if (j >= 0 && cc[j] > 0.01) return zMid[k * C + i] - zs[i]; } return NaN; })); add('share of cells with cumulus cover > 0.01', binMean((i) => { for (let k = 0; k < K; k++) { const j = k * C + i - off; if (j >= 0 && cc[j] > 0.01) return 1; } return 0; })); add('cumulus top above 2 km, share', binMean((i) => { for (let k = 0; k < K; k++) { const j = k * C + i - off; if (j >= 0 && cc[j] > 0.01) return zMid[k * C + i] - zs[i] > 2000 ? 1 : 0; } return 0; })); }
  for (let k = K - 1; k >= 18; k--) { for (const b of [2, 5]) { const name = binName(BINS[b]); add(`profile k=${k} ${name} z (m)`, binMean((i) => zMid[k * C + i] - zs[i])); add(`profile k=${k} ${name} T (C)`, binMean((i) => T(k, i) - 273.15)); add(`profile k=${k} ${name} RH`, binMean((i) => q[k * C + i] / saturationHumidity(T(k, i), pi[i] * sigmaMid[k]))); add(`profile k=${k} ${name} q (g/kg)`, binMean((i) => 1000 * q[k * C + i])); } }
  for (const zq of HEIGHTS) { const v = binMean((i) => atHeight(tv, i, zq)); add(`Tv at ${zq} m less ${binName(BINS[REF])}'s (K)`, v.map((x) => x - v[REF])); }
  const bandOf = new Float64Array(C * (BANDS.length - 1));
  for (const i of cells) {
    const zk = [0], lr = [Math.log(pi[i] / (Rd * tv[(K - 1) * C + i]))];
    for (let k = K - 1; k >= 0; k--) { zk.push(zMid[k * C + i] - zs[i]); lr.push(Math.log(pi[i] * sigmaMid[k] / (Rd * tv[k * C + i]))); }
    let j = 0;
    for (let b = 0; b < BANDS.length - 1; b++) {
      let sum = 0;
      for (let z = BANDS[b] + 25; z < BANDS[b + 1]; z += 50) { while (j < zk.length - 2 && z > zk[j + 1]) j++; const t = (z - zk[j]) / (zk[j + 1] - zk[j]); sum += Math.exp(lr[j] + t * (lr[j + 1] - lr[j])) * 50; }
      bandOf[b * C + i] = g * sum / 100;
    }
  }
  for (let b = 0; b < BANDS.length - 1; b++) { const v = binMean((i) => bandOf[b * C + i]); add(`SLP weight ${BANDS[b] / 1000}-${BANDS[b + 1] / 1000} km less ${binName(BINS[REF])}'s (hPa)`, v.map((x) => x - v[REF])); }
  files++;
}
const f = (x, d = 2) => (Number.isFinite(x) ? x.toFixed(d) : '—');
console.log(`== ${LABEL}: mean of ${files} states (${FILES.map((p) => p.replace(/.*day/, '')).join(', ')}), sea cells within ${LAT} deg`);
console.log('| | ' + BINS.map(binName).join(' | ') + ' |\n|---|' + BINS.map(() => '---|').join(''));
for (const [name, s] of Object.entries(sums)) if (!name.startsWith('profile')) console.log(`| ${name} | ${Array.from(s, (x) => f(x / files)).join(' | ')} |`);
const m = (name) => Array.from(sums[name], (x) => x / files);
const slp = m('sea-level pressure (hPa)'), tau = m('stress tau_x (mN/m2)'), uu = m('u lowest layer (m/s)');
console.log(`summary ${LABEL}: SLP 180-160W less 120W-100W ${f(slp[2] - slp[5])} hPa; stress 160W-140W ${f(tau[3], 1)}, 140W-120W ${f(tau[4], 1)}, 120W-100W ${f(tau[5], 1)} mN/m2 (mean ${f((tau[3] + tau[4] + tau[5]) / 3, 1)}); u ${f(uu[3])}, ${f(uu[4])}, ${f(uu[5])} m/s (mean ${f((uu[3] + uu[4] + uu[5]) / 3)})`);

console.log('profile (bin means): layer, then for 180-160W and 120W-100W: z m, T C, RH, q g/kg');
for (let k = K - 1; k >= 18; k--) { const cols = []; for (const b of [2, 5]) { const name = binName(BINS[b]); const v = (what) => sums[`profile k=${k} ${name} ${what}`][b] / files; cols.push(`${f(v('z (m)'), 0).padStart(5)} ${f(v('T (C)'), 1).padStart(5)} ${f(v('RH'), 2)} ${f(v('q (g/kg)'), 1).padStart(4)}`); } console.log(`k=${String(k).padStart(2)}  ${cols.join('   |   ')}`); }
