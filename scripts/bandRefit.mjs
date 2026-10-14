// A saved state with the interior classes of the tropical band refit to WOA23.
//
//   node scripts/bandRefit.mjs SOURCE.bin OUT.bin [--inner 10] [--outer 15] [--json REPORT.json]
//
// Every sea column within `outer` degrees of the equator takes the weight
// w = 1 within `inner` degrees, falling linearly to 0 at `outer`. Its
// classes become the blend w·(WOA23's classes, atlasColumns' way) +
// (1 − w)·(the state's) in thickness, heat and salt. The mixed layer
// (h₀, Q₀, W₀), the sea level η and every velocity are kept, and the
// atlas classes are cut from the top to begin below h₀ and rescaled to
// the water the state's classes hold, so each column still sums to
// D + η with its tokens at EPS. Columns beyond `outer`, land, and
// columns the atlas does not cover are left bit for bit. Heat and salt
// are not conserved; the report gives the heat change by band and the
// thickness of the classes of 15–23 °C in the eastern equatorial box.
// The class list is the state's own, and the output loads like any state
// (a 44-class state rebins at load).
import { readFileSync, writeFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { atlasColumns, savedDensities, EPS, LAYER_DENSITIES, LAYER_SALINITIES } from '../js/ocean/layered.module.js';
import { labelTemperature } from '../js/ocean/seawater.module.js';
import { loadClimatology, CLIMATOLOGY_FILE } from '../js/ocean/climatology.module.js';

const DEG = 180 / Math.PI, RHO_CP = 1025 * 3985, KELVIN = 273.15;
export const BAND = { inner: 10, outer: 15 };
export const EUC_BOX = { south: -5, north: 5, west: -140, east: -90, coldest: 15, warmest: 23 };

export const bandWeight = (latDegrees, { inner = BAND.inner, outer = BAND.outer } = {}) => {
  const a = Math.abs(latDegrees);
  return a <= inner ? 1 : a >= outer ? 0 : (outer - a) / (outer - inner);
};

function classSalinities(classes) {
  return classes.map((rho) => {
    const j = LAYER_DENSITIES.findIndex((d) => Math.abs(d - rho) < 1e-6);
    if (j < 0) throw new Error(`class ${rho} is not one of LAYER_DENSITIES; pass its salinities`);
    return LAYER_SALINITIES[j];
  });
}

export function refitBand(ocean, mesh, atlas, { ice = null, inner = BAND.inner, outer = BAND.outer, salinities = null, referenceS = 35, density = 1025 } = {}) {
  const C = mesh.nCells, classes = savedDensities(ocean, C);
  if (!classes) throw new Error('the ocean carries no class list');
  const L = classes.length + 1;
  if (ocean.h.length !== L * C) throw new Error(`h holds ${ocean.h.length} values, not ${L} layers of ${C} cells`);
  const at = (k, i) => k * C + i, { h, T, S, eta } = ocean;
  const rho = [density, ...classes], labelS = [referenceS, ...(salinities ?? classSalinities(classes))];
  const labelT = rho.map((r, k) => labelTemperature(r, labelS[k]));
  const column = new Float64Array(C), weight = new Float64Array(C), D = new Float64Array(C), sea = new Uint8Array(C);
  for (let i = 0; i < C; i++) {
    let sum = 0;
    for (let k = 0; k < L; k++) sum += h[at(k, i)];
    column[i] = sum; sea[i] = sum > 1 ? 1 : 0; D[i] = sum - eta[i];
    weight[i] = sea[i] ? bandWeight(mesh.latCell[i] * DEG, { inner, outer }) : 0;
  }
  const candidate = Uint8Array.from(weight, (w) => (w > 0 ? 1 : 0));
  const hA = new Float64Array(L * C), QA = new Float64Array(L * C), WA = new Float64Array(L * C), T0A = new Float64Array(C);
  const filled = atlasColumns(mesh, atlas, { D, cellOcean: candidate, ice: ice ?? new Float64Array(C), rho, labelT, labelS, h: hA, Q: QA, W: WA, T0: T0A, shallowestMixedDepth: 1 });
  const stats = { candidates: 0, refit: 0, noAtlas: 0, thin: 0, cut: 0, gap: 0, tokenized: 0 };
  const dQ = new Float64Array(C), dW = new Float64Array(C), refit = new Uint8Array(C);
  const euc = classes.map((r, k) => ({ rho: r, label: labelT[k + 1] - KELVIN, area: 0, before: 0, after: 0 }));
  const inEuc = (i) => {
    const lat = mesh.latCell[i] * DEG, lon = ((mesh.lonCell[i] * DEG + 540) % 360) - 180;
    return lat >= EUC_BOX.south && lat <= EUC_BOX.north && lon >= EUC_BOX.west && lon <= EUC_BOX.east;
  };
  const boxCells = [];
  for (let i = 0; i < C; i++) if (sea[i] && inEuc(i)) boxCells.push(i);
  let boxArea = 0;
  for (const i of boxCells) {
    boxArea += mesh.areaCell[i];
    for (let k = 1; k < L; k++) euc[k - 1].before += mesh.areaCell[i] * h[at(k, i)];
  }
  const a = new Float64Array(L), qa = new Float64Array(L), wa = new Float64Array(L);
  for (let i = 0; i < C; i++) {
    if (!candidate[i]) continue;
    stats.candidates++;
    if (!filled[i]) { stats.noAtlas++; continue; }
    const h0 = h[at(0, i)], total = column[i] - h0;
    if (!(total > (L - 1) * EPS + 1e-6)) { stats.thin++; continue; }
    for (let k = 1; k < L; k++) { a[k] = hA[at(k, i)]; qa[k] = QA[at(k, i)]; wa[k] = WA[at(k, i)]; }
    let cut = h0 - hA[i];
    if (cut > 0) stats.cut++; else if (cut < 0) stats.gap++;
    for (let k = 1; k < L && cut > 0; k++) {
      const take = Math.min(Math.max(0, a[k] - EPS), cut);
      if (take <= 0) continue;
      const f = (a[k] - take) / a[k];
      a[k] -= take; qa[k] *= f; wa[k] *= f; cut -= take;
    }
    let scalable = true;
    for (let pass = 0; pass < L && scalable; pass++) {
      let tokens = 0, mass = 0;
      for (let k = 1; k < L; k++) { if (a[k] <= EPS + 1e-12) tokens++; else mass += a[k]; }
      const target = total - tokens * EPS, s = target / mass;
      if (!(target > 0) || !(mass > 0)) { scalable = false; break; }
      let violated = false;
      for (let k = 1; k < L; k++) if (a[k] > EPS + 1e-12 && a[k] * s < EPS) { a[k] = EPS; qa[k] = labelT[k] * EPS; wa[k] = labelS[k] * EPS; violated = true; stats.tokenized++; }
      if (violated) continue;
      for (let k = 1; k < L; k++) if (a[k] > EPS + 1e-12) { a[k] *= s; qa[k] *= s; wa[k] *= s; }
      break;
    }
    if (!scalable) { stats.thin++; continue; }
    const w = weight[i];
    let sum = 0, thickest = 1;
    const nh = new Float64Array(L), nq = new Float64Array(L), nw = new Float64Array(L);
    for (let k = 1; k < L; k++) {
      const n = at(k, i), oh = h[n], oq = oh * T[n], ow = oh * S[n];
      nh[k] = w * a[k] + (1 - w) * oh; nq[k] = w * qa[k] + (1 - w) * oq; nw[k] = w * wa[k] + (1 - w) * ow;
      sum += nh[k];
      if (nh[k] > nh[thickest]) thickest = k;
    }
    const residual = total - sum;
    if (residual !== 0) { const f = (nh[thickest] + residual) / nh[thickest]; nh[thickest] += residual; nq[thickest] *= f; nw[thickest] *= f; }
    let heat = 0, salt = 0;
    for (let k = 1; k < L; k++) {
      const n = at(k, i);
      heat += nq[k] - h[n] * T[n]; salt += nw[k] - h[n] * S[n];
      h[n] = nh[k]; T[n] = nq[k] / nh[k]; S[n] = nw[k] / nh[k];
    }
    dQ[i] = heat; dW[i] = salt; refit[i] = 1; stats.refit++;
  }
  for (const i of boxCells) for (let k = 1; k < L; k++) euc[k - 1].after += mesh.areaCell[i] * h[at(k, i)];
  const bands = [], edges = [-15, -10, -5, 0, 5, 10, 15];
  for (let b = 0; b + 1 < edges.length; b++) {
    let area = 0, heat = 0, columns = 0;
    for (let i = 0; i < C; i++) {
      const lat = mesh.latCell[i] * DEG;
      if (!refit[i]) continue;
      if (lat < edges[b] || lat >= edges[b + 1]) continue;
      area += mesh.areaCell[i]; heat += mesh.areaCell[i] * dQ[i] * RHO_CP; columns++;
    }
    bands.push({ south: edges[b], north: edges[b + 1], columns, meanJPerM2: area ? heat / area : 0, totalJ: heat });
  }
  let totalJ = 0, extreme = { value: 0, lat: null, lon: null };
  for (let i = 0; i < C; i++) {
    const q = dQ[i] * RHO_CP; totalJ += mesh.areaCell[i] * q;
    if (Math.abs(q) > Math.abs(extreme.value)) extreme = { value: q, lat: +(mesh.latCell[i] * DEG).toFixed(1), lon: +(((mesh.lonCell[i] * DEG + 540) % 360) - 180).toFixed(1) };
  }
  const warm = euc.filter((c) => c.label >= EUC_BOX.coldest && c.label <= EUC_BOX.warmest).map((c) => ({ rho: c.rho, labelC: +c.label.toFixed(1), beforeM: c.before / boxArea, afterM: c.after / boxArea }));
  return { stats, bands, totalJ, extreme, box: { ...EUC_BOX, cells: boxCells.length }, warmClasses: warm, dQ, dW, weight, refit };
}

export function summary(report) {
  const f = (x, n = 2) => (Number.isFinite(x) ? x.toFixed(n) : '–');
  const out = [`columns: ${report.stats.candidates} within the outer latitude, ${report.stats.refit} refit, ${report.stats.noAtlas} without WOA23, ${report.stats.thin} too thin; atlas cut at the top in ${report.stats.cut}, gap under a thinner mixed layer in ${report.stats.gap}, ${report.stats.tokenized} classes made tokens by the rescale`];
  out.push(`heat change ${f(report.totalJ / 1e21, 3)} ZJ in all; largest column ${f(report.extreme.value / 1e9, 2)} GJ/m² at ${report.extreme.lat}, ${report.extreme.lon}`);
  out.push(...report.bands.map((b) => `  ${String(b.south).padStart(3)} to ${String(b.north).padStart(3)}: ${b.columns} columns, mean ${f(b.meanJPerM2 / 1e9, 3)} GJ/m²`));
  out.push(`classes of ${report.box.coldest}-${report.box.warmest} °C, mean thickness over ${report.box.south}..${report.box.north}N ${report.box.west}..${report.box.east}E (m), before to after:`);
  out.push(...report.warmClasses.map((c) => `  ${c.rho} (${c.labelC} °C): ${f(c.beforeM)} to ${f(c.afterM)}`));
  return out.join('\n');
}

if (process.argv[1] && fileURLToPath(import.meta.url) === process.argv[1]) {
  const args = process.argv.slice(2), files = [], options = {};
  let jsonOut = null;
  for (let n = 0; n < args.length; n++) {
    if (args[n] === '--inner') options.inner = Number(args[++n]);
    else if (args[n] === '--outer') options.outer = Number(args[++n]);
    else if (args[n] === '--json') jsonOut = args[++n];
    else files.push(args[n]);
  }
  if (files.length !== 2 || !((options.inner ?? BAND.inner) < (options.outer ?? BAND.outer))) {
    console.error('usage: node scripts/bandRefit.mjs SOURCE.bin OUT.bin [--inner 10] [--outer 15] [--json REPORT.json]');
    process.exit(2);
  }
  const { decodeState, encodeState } = await import('../js/stateFile.module.js');
  const { stateGrid } = await import('../js/stateGrid.module.js');
  const { topographyFromInt16 } = await import('../js/geography.module.js');
  const state = await decodeState(new Uint8Array(readFileSync(files[0])));
  const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
  const mesh = stateGrid(state, topography).mesh;
  const atlas = await loadClimatology(fileURLToPath(new URL(`../${CLIMATOLOGY_FILE}`, import.meta.url)));
  const report = refitBand(state.ocean, mesh, atlas, { ice: state.ice, ...options });
  const f64 = [];
  for (const [k, v] of Object.entries(state)) {
    if (v instanceof Float64Array) f64.push(k);
    else if ((k === 'ocean' || k === 'land') && v) for (const [i, x] of Object.entries(v)) if (x instanceof Float64Array) f64.push(`${k}.${i}`);
  }
  const bytes = encodeState(state, { f64 });
  writeFileSync(files[1], bytes);
  const back = await decodeState(new Uint8Array(readFileSync(files[1])));
  const same = back.ocean.h.length === state.ocean.h.length && back.ocean.h.every((x, n) => Object.is(x, state.ocean.h[n])) && back.ocean.T.every((x, n) => Object.is(x, state.ocean.T[n]));
  console.log(`N=${state.N} day ${state.day}: ${files[1]} (${bytes.length} bytes, round trip ${same ? 'identical' : 'DIFFERS'})`);
  console.log(summary(report));
  if (jsonOut) {
    const { dQ, dW, weight, refit, ...plain } = report;
    writeFileSync(jsonOut, JSON.stringify({ source: files[0], out: files[1], day: state.day, N: state.N, inner: options.inner ?? BAND.inner, outer: options.outer ?? BAND.outer, ...plain }, null, 1));
  }
  if (!same) process.exit(1);
}
