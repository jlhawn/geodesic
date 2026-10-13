// One coupled CPU day from a saved state, judging the mixed layer's mixing:
//   [WORKERS=6] [DAYS=1] [GATE=2.3,3.6] [OCEAN_OPTS='{"shallowestMixedDepth":50,...}'] node scripts/mixingCheck.mjs [runs/eleven64_day1825.bin]
// OCEAN_OPTS (or OCEAN) sets the ocean's options. It reports the share of
// sea cells entraining by shear each ocean step and the volume a day, the
// momentum depth H_m by band, the depth the retreat asks for where the
// floor binds, the mixed layer's share of the stress where the stress
// depth S falls under 0.8 H, and the cells whose pending light exceeds its
// bound. The asked depth is recomputed at the first ocean step from the
// state's buoyancy memory and the coupled stress, under the options in
// force and under the retreat without the package, on the cells where the
// latter binds. It fails when the median asked depth there falls below the
// package-off value and, with the package on at a 50 m floor, when their
// ratio leaves GATE. The default 2.3,3.6 is the model's own form, the
// package-off rule's faded constant taken at the step's depth (the floor);
// taken at the asked depth, as the probe does, the band is 1.4,2.2. It
// refuses to start beside another `node scripts/` or `node --test` job.
import { readFileSync } from 'node:fs';
import { execSync } from 'node:child_process';
import { Grid } from '../js/grid.module.js';
import { createParallelModel } from '../js/parallel.module.js';
import { withCadence } from '../js/cadence.module.js';
import { decodeState, savedLevels } from '../js/stateFile.module.js';
import { topographyFromInt16, rebalanceSurfacePressure } from '../js/geography.module.js';
import { savedDeckField, DECK_FIELDS, savedMoistField, MOIST_FIELDS, savedRadiationField, RADIATION_FIELDS } from '../js/physics/regrid.module.js';
import { PV_FLOOR } from '../js/ocean/layered.module.js';
import { cellVector } from '../js/dynamics/operators.module.js';
import { stateGrid, gridMesh } from '../js/stateGrid.module.js';

const table = execSync('ps -axo pid=,ppid=,command=', { encoding: 'utf8' }).split('\n').filter((row) => row.trim()).map((row) => { const [pid, ppid, ...command] = row.trim().split(/\s+/); return { pid: Number(pid), ppid: Number(ppid), command: command.join(' ') }; });
const parent = new Map(table.map((row) => [row.pid, row.ppid])), mine = new Set();
for (let pid = process.pid; pid > 1 && !mine.has(pid); pid = parent.get(pid) ?? 0) mine.add(pid);
const others = table.filter((row) => !mine.has(row.pid) && /node (scripts\/|--test)/.test(row.command)).map((row) => `${row.pid} ${row.command.slice(0, 160)}`).join('\n');
if (others) {
  console.error(`refusing to start: another node job is running, and the GPU rule allows one at a time:\n${others}`);
  process.exit(2);
}

const FILE = process.argv[2] ?? new URL('../runs/eleven64_day1825.bin', import.meta.url).pathname;
const WORKERS = Number(process.env.WORKERS ?? 6), DAYS = Number(process.env.DAYS ?? 1);
const OCEAN = JSON.parse(process.env.OCEAN_OPTS ?? process.env.OCEAN ?? '{}');
const BINS = [0, 2, 5, 10, 20, 30, 40, 50], GATE = (process.env.GATE ?? '2.3,3.6').split(',').map(Number), U10_RATIO = 1.19e-3, LA_MIN = 0.1;
const t0 = performance.now(), seconds = () => ((performance.now() - t0) / 1000).toFixed(0);

const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
const N = saved.N, STEP = 1350 * 16 / N, found = process.env.RELAX ? { mesh: gridMesh(new Grid(N, { relax: Number(process.env.RELAX) })), relax: Number(process.env.RELAX) } : stateGrid(saved, topography);
const options = { topography, terrain: true, levels: savedLevels(saved), ...(Object.keys(OCEAN).length ? { ocean: OCEAN } : {}) };
Object.assign(options, withCadence(options, STEP, N));
const model = await createParallelModel(found.mesh, options, WORKERS);
const { mesh, ocean, state } = model, C = mesh.nCells, E = mesh.nEdges, L = ocean.layers;

['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'].forEach((name, a) => state[a].set(saved[name]));
model.seaIce.load(state[6], saved.concentration ?? null);
for (const field of Object.keys(DECK_FIELDS)) model.radiation[field].set(savedDeckField(saved, field, model));
for (const field of Object.keys(MOIST_FIELDS)) model.moist[field].set(savedMoistField(saved, field, model));
for (const field of Object.keys(RADIATION_FIELDS)) model.radiation[field].set(savedRadiationField(saved, field, model));
const restore = (target, values) => { if (target && values && values.length === target.length) target.set(values); };
const bl = model.boundaryLayer;
restore(bl.depth, saved.boundaryDepth); restore(bl.mixingTop, saved.mixingTop); restore(bl.regime, saved.boundaryRegime); restore(bl.buoyancyFlux, saved.boundaryBuoyancy);
restore(model.surface.windSpeed, saved.windSpeed); restore(model.radiation.evaporation, saved.evaporation);
if (model.exchange && !model.exchange.fixed) restore(model.exchange.heat, saved.exchangeHeat);
if (model.exchange) restore(model.exchange.wind, saved.exchangeWind);
const K0 = model.moist.cumulusK0 ?? 0, K = model.core.K;
for (const name of ['cumulusCover', 'cumulusWater']) {
  const target = model.moist[name], values = saved[name];
  if (!target || !values) continue;
  target.fill(0);
  if (values.length === target.length) target.set(values); else if (values.length === (K - K0) * C) target.set(values, K0 * C);
}
restore(model.moist.subcloudVirtual, saved.subcloudVirtual);
restore(model.seaIce.snow, saved.seaSnow); restore(model.seaIce.snowAlbedo, saved.seaSnowAlbedo);
model.time = saved.time;
ocean.load(saved.ocean, state[3], state[6]);
restore(model.seaIce.oceanFlux, saved.ocean.flux);
if (model.land && saved.land) model.land.load(saved.land, state[6]);
model.dragsDue();
model.restartPrecipitation();
model.core.diagnose(state[0], state[1], null, null);
rebalanceSurfacePressure(model.core, state[0], state[1], model.surfaceGeopotential, model.surfaceGeopotential);

const deg = 180 / Math.PI, rho0 = OCEAN.density ?? 1025, rhoCp = rho0 * (OCEAN.specificHeat ?? 3985);
const F = OCEAN.shallowestMixedDepth ?? 50, stirring = OCEAN.stirring ?? 0.8, stirringDepth = OCEAN.stirringDepth ?? 100;
const retreatStirring = OCEAN.retreatStirring ?? 0, langmuirCap = OCEAN.langmuir ?? 0, packaged = retreatStirring > 0 || langmuirCap > 0;
const momentumDepth = OCEAN.momentumDepth ?? 50, ekmanDepth = OCEAN.ekmanDepth ?? 0, stressThickness = OCEAN.stressThickness ?? PV_FLOOR;
const { latCell, lonCell, latEdge, areaCell, cellsOnEdge } = mesh;
const lat = (i) => latCell[i] * deg, lon = (i) => ((((lonCell[i] * deg) + 540) % 360) - 180);
const ice = state[6];
const column = new Float64Array(C);
const columns = () => { column.fill(0); for (let k = 0; k < L; k++) for (let i = 0; i < C; i++) column[i] += ocean.h[k * C + i]; };
columns();
const box = (s, n, w, e) => (i) => lat(i) >= s && lat(i) <= n && (w <= e ? lon(i) >= w && lon(i) <= e : lon(i) >= w || lon(i) <= e);
const REGIONS = {
  'central Pacific 5S–5N 170E–120W': box(-5, 5, 170, -120),
  'W Pacific 5S–5N 140–170E': box(-5, 5, 140, 170),
  'cold tongue 5S–5N 130–90W': box(-5, 5, -130, -90),
  'Arctic shelf (66N, 50–180E, D < 400 m)': (i) => lat(i) > 66 && lon(i) >= 50 && ocean.D[i] < 400,
  'other shelves (D < 200 m)': (i) => ocean.D[i] < 200 && !(lat(i) > 66 && lon(i) >= 50),
};
const regionCells = Object.fromEntries(Object.entries(REGIONS).map(([name, test]) => [name, Array.from({ length: C }, (_, i) => i).filter((i) => ocean.cellOcean[i] && test(i))]));
const seaCells = Array.from({ length: C }, (_, i) => i).filter((i) => ocean.cellOcean[i]);
const BANDS = [['60–20S', -60, -20], ['15S–15N', -15, 15], ['20–60N', 20, 60], ['60–80N', 60, 80]];

function erfc(x) {
  const z = Math.abs(x), t = 1 / (1 + 0.5 * z);
  const r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
  return x >= 0 ? r : 2 - r;
}
function stokesSurfaceLayer(U10, depth) {
  const us = 0.0162 * U10, hm0 = 0.0246 * U10 * U10, fp = 0.877 * 9.81 / (2 * Math.PI * 1.075 * U10), fm = 1.296 * fp;
  const transport = 0.125 * Math.PI * 0.667 * fm * hm0 * hm0, kphil = 0.176 * us / transport, kstar = 2.56 * kphil, z0 = 0.2 * Math.abs(depth);
  const r1 = (0.151 / (kphil * z0) - 0.84) * (1 - Math.exp(-2 * kphil * z0));
  const r2 = -(0.84 + 0.0591 / (kphil * z0)) * Math.sqrt(2 * Math.PI * kphil * z0) * erfc(Math.sqrt(2 * kphil * z0));
  const r3 = (0.0632 / (kstar * z0) + 0.125) * (1 - Math.exp(-2 * kstar * z0));
  const r4 = (0.125 + 0.0946 / (kstar * z0)) * Math.sqrt(2 * Math.PI * kstar * z0) * erfc(Math.sqrt(2 * kstar * z0));
  return us * (0.715 + r1 + r2 + r3 + r4);
}
function langmuirFactor(ustar, h0, cap) {
  if (!(cap > 0) || !(ustar > 0)) return 1;
  const la = Math.max(LA_MIN, Math.sqrt(ustar / stokesSurfaceLayer(ustar / U10_RATIO, h0)));
  return Math.min(cap, Math.sqrt(1 + la ** -2 / 1.5 ** 2 + la ** -4 / 5.4 ** 4));
}
const fadeOff = (h0) => 2 * stirring * Math.exp(-h0 / stirringDepth);
function askedDepth(tau, loss, h0, inForce) {
  const ustar3 = Math.pow(tau / rho0, 1.5);
  if (!inForce) return fadeOff(h0) * ustar3 / -loss;
  const constant = retreatStirring > 0 ? retreatStirring : fadeOff(h0);
  return constant * langmuirFactor(Math.sqrt(tau / rho0), h0, langmuirCap) * ustar3 / -loss;
}

const fieldOf = (value) => (value == null ? null : ArrayBuffer.isView(value) ? value : value.water ?? value.volume ?? null);
const ulp = (x) => { const a = Math.abs(Math.fround(x)); return a === 0 ? 2 ** -149 : 2 ** (Math.floor(Math.log2(a)) - 23); };
const tauCell = new Float64Array(3 * C);
let previousH0 = Float64Array.from(ocean.h.subarray(0, C)), previousTime = model.time, oceanSeconds = 0;
let shearedBefore = fieldOf(ocean.sheared) ? Float64Array.from(fieldOf(ocean.sheared)) : null;
const shear = { steps: 0, share: 0, maxShare: 0, volume: 0, regions: Object.fromEntries(Object.keys(REGIONS).map((name) => [name, 0])) };
const windDepths = BANDS.map(() => []);
const modelAsked = [];
const stress = { edges: 0, low: 0, tropics: 0, lowTropics: 0, shareLow: [], ratioLow: [], shareThin: [] };
const pen = { stepsOver: 0, maxOver: 0, causes: { iced: 0, entrained: 0, binade: 0, other: 0 }, worst: 0, worstAt: null };
let first = null, oceanSteps = 0;

function sampleAsked() {
  const loss = saved.ocean.buoyancy;
  if (!loss || loss.length !== C) return { missing: 'the state carries no buoyancy memory' };
  cellVector(mesh, ocean.stress, tauCell);
  const rows = [];
  for (const i of seaCells) {
    if (ice[i] > 0 || !(loss[i] < -1e-9) || column[i] < F + 10) continue;
    const h0 = ocean.h[i];
    if (h0 > F + 1) continue;
    const tau = Math.hypot(tauCell[3 * i], tauCell[3 * i + 1], tauCell[3 * i + 2]);
    if (!(tau > 0)) continue;
    const off = askedDepth(tau, loss[i], h0, false);
    if (!(off < F)) continue;
    rows.push({ i, off, on: askedDepth(tau, loss[i], h0, true), model: fieldOf(ocean.retreatDepth)?.[i], tau, loss: loss[i], h0 });
  }
  return { rows };
}

function edgeColumn(e, out) {
  const a = cellsOnEdge[2 * e], b = cellsOnEdge[2 * e + 1], h = ocean.h, eta = ocean.eta;
  out[0] = 0.5 * (h[a] + h[b]);
  let sum = out[0];
  for (let k = 1; k < L; k++) { out[k] = Math.min(h[k * C + a], h[k * C + b]); sum += out[k]; }
  const sill = Math.max(0.01, Math.min(ocean.D[a], ocean.D[b]) + 0.5 * (eta[a] + eta[b]));
  if (sum > sill) for (let k = 0; k < L; k++) out[k] *= sill / sum;
  return out;
}

function sampleOceanStep() {
  oceanSteps++;
  columns();
  const sheared = fieldOf(ocean.sheared);
  if (sheared) {
    let cells = 0;
    for (const i of seaCells) { const d = sheared[i] - shearedBefore[i]; if (d > 1e-9) { cells++; shear.volume += areaCell[i] * d; } }
    for (const [name, list] of Object.entries(regionCells)) if (list.length) shear.regions[name] += list.filter((i) => sheared[i] - shearedBefore[i] > 1e-9).length / list.length;
    shear.steps++; shear.share += cells / seaCells.length; shear.maxShare = Math.max(shear.maxShare, cells / seaCells.length);
    shearedBefore = Float64Array.from(sheared);
  }
  const windDepth = fieldOf(ocean.windDepth);
  if (windDepth) for (let e = 0; e < E; e++) {
    if (!ocean.edgeOcean[e]) continue;
    const l = latEdge[e] * deg;
    BANDS.forEach(([, lo, hi], b) => { if (l >= lo && l < hi) windDepths[b].push(windDepth[e]); });
  }
  const retreat = fieldOf(ocean.retreatDepth);
  if (retreat) for (const i of seaCells) {
    if (ice[i] > 0 || column[i] < F + 10 || ocean.h[i] > F + 1) continue;
    const x = retreat[i];
    if (Number.isFinite(x) && x > 0 && x < F) modelAsked.push(x);
  }
  const he = new Float64Array(L);
  for (let e = 0; e < E; e++) {
    if (!ocean.edgeOcean[e]) continue;
    edgeColumn(e, he);
    const reach = Math.max(ekmanDepth > 0 && windDepth ? windDepth[e] : momentumDepth, he[0]);
    let top = 0, depth = 0;
    for (let k = 0; k < L; k++) { const part = Math.min(he[k], Math.max(0, reach - top)); if (k === 0 || he[k] >= stressThickness) depth += part; top += he[k]; if (top >= reach) break; }
    const share = Math.min(he[0], reach) / depth, tropical = Math.abs(latEdge[e] * deg) <= 20;
    stress.edges++; if (tropical) stress.tropics++;
    if (he[0] < reach) stress.shareThin.push(share);
    if (depth < 0.8 * reach) { stress.low++; if (tropical) stress.lowTropics++; stress.shareLow.push(share); stress.ratioLow.push(depth / reach); }
  }
  let over = 0;
  for (const i of seaCells) {
    const p = ocean.penetrated ? ocean.penetrated[i] : 0;
    if (!p) continue;
    const trans = ocean.transmission ? ocean.transmission[i] : 0;
    let z = ocean.h[i], spacing = ulp(ocean.Q[i]);
    for (let k = 1; k < L && z < 300; k++) { if (ocean.h[k * C + i] > 5) spacing = Math.max(spacing, ulp(ocean.Q[k * C + i])); z += ocean.h[k * C + i]; }
    const light = oceanSeconds * 1100 * trans, bound = 0.5 * rhoCp * spacing + light;
    if (!(Math.abs(p) > bound)) continue;
    over++;
    const cause = ice[i] > 0 ? 'iced' : trans < 1e-3 ? 'entrained' : Math.abs(p) <= rhoCp * spacing + light ? 'binade' : 'other';
    pen.causes[cause]++;
    if (Math.abs(p) / bound > pen.worst) { pen.worst = Math.abs(p) / bound; pen.worstAt = `${lat(i).toFixed(1)}, ${lon(i).toFixed(1)} (${cause}, h0 ${ocean.h[i].toFixed(1)} m)`; }
  }
  pen.stepsOver += over; pen.maxOver = Math.max(pen.maxOver, over);
}

console.log(`${FILE}: N=${N} on ${found.relax} Lloyd passes; floor ${F} m; ${packaged ? `retreat constant ${retreatStirring || 'faded'}, Langmuir cap ${langmuirCap || 'none'}` : 'retreat without the package'}; ekmanDepth ${ekmanDepth}, stressThickness ${stressThickness}; OCEAN ${JSON.stringify(OCEAN)}; setup ${seconds()} s`);
const steps = Math.round(DAYS * 86400 / STEP);
for (let n = 0; n < steps; n++) {
  model.step(STEP);
  let moved = false;
  for (let i = 0; i < C && !moved; i++) if (ocean.h[i] !== previousH0[i]) moved = true;
  if (!moved) continue;
  oceanSeconds = model.time - previousTime; previousTime = model.time;
  previousH0 = Float64Array.from(ocean.h.subarray(0, C));
  if (!first) { columns(); first = sampleAsked(); }
  sampleOceanStep();
}
const days = steps * STEP / 86400;

const f = (x, n = 1) => (x == null || !Number.isFinite(x) ? 'n/a' : x.toFixed(n));
const sorted = (values) => Float64Array.from(values).sort();
const quantile = (s, q) => (s.length ? s[Math.min(s.length - 1, Math.floor(q * s.length))] : NaN);
const histogram = (values) => {
  const s = sorted(values), bins = BINS.slice(0, -1).map((lo, b) => `${lo}–${BINS[b + 1]} ${f(100 * values.filter((x) => x >= lo && x < BINS[b + 1]).length / values.length, 0)} %`);
  return `${bins.join(', ')}, ≥ ${BINS.at(-1)} ${f(100 * values.filter((x) => x >= BINS.at(-1)).length / values.length, 0)} %; median ${f(quantile(s, 0.5))} m, under 5 m ${f(100 * values.filter((x) => x < 5).length / values.length, 0)} %, under 10 m ${f(100 * values.filter((x) => x < 10).length / values.length, 0)} %`;
};
console.log(`${f(days, 2)} days, ${oceanSteps} ocean steps of ${f(oceanSeconds, 0)} s, ${seconds()} s`);
if (shear.steps) {
  const regions = Object.entries(shear.regions).map(([name, s]) => `${name} ${f(100 * s / shear.steps)} % of ${regionCells[name].length} cells`).join('; ');
  console.log(`shear entrainment: ${f(100 * shear.share / shear.steps, 2)} % of sea cells per ocean step (largest ${f(100 * shear.maxShare, 2)} %), ${f(shear.volume / 1e9 / days, 0)} km³ a day; ${regions}`);
} else console.log('shear entrainment: n/a (the ocean has no sheared diagnostic)');
if (windDepths.some((v) => v.length)) console.log(`H_m p5/p50/p95 (m): ${BANDS.map(([name], b) => { const s = sorted(windDepths[b]); return `${name} ${[0.05, 0.5, 0.95].map((q) => f(quantile(s, q))).join('/')}`; }).join('; ')}`);
else console.log('H_m p5/p50/p95: n/a (the ocean has no windDepth)');

const checks = [];
if (first.missing) console.log(`asked depth where the floor binds: n/a (${first.missing})`);
else if (!first.rows.length) console.log(`asked depth where the floor binds: no ice-free heated column on its ${F} m floor asks for less`);
else {
  const rows = first.rows, off = sorted(rows.map((r) => r.off)), on = sorted(rows.map((r) => r.on)), ratio = quantile(on, 0.5) / quantile(off, 0.5);
  const perCell = sorted(rows.map((r) => r.on / r.off));
  console.log(`asked depth at the first ocean step, ${rows.length} ice-free heated columns on the ${F} m floor where the retreat without the package asks for less (state's buoyancy memory, coupled stress):`);
  console.log(`  without the package: ${histogram(rows.map((r) => r.off))}`);
  console.log(`  options in force:    ${histogram(rows.map((r) => r.on))}`);
  console.log(`  median ${f(quantile(on, 0.5))} against ${f(quantile(off, 0.5))} m, ratio ${f(ratio, 2)}; per-column ratio p5/p50/p95 ${[0.05, 0.5, 0.95].map((q) => f(quantile(perCell, q), 2)).join('/')}; median stress ${f(quantile(sorted(rows.map((r) => r.tau)), 0.5), 4)} N/m², buoyancy memory ${quantile(sorted(rows.map((r) => r.loss)), 0.5).toExponential(2)} m²/s³, h0 ${f(quantile(sorted(rows.map((r) => r.h0)), 0.5))} m`);
  const both = rows.filter((r) => Number.isFinite(r.model) && r.model > 0);
  if (both.length) { const m = sorted(both.map((r) => r.model / r.on)); console.log(`  ocean.retreatDepth against the recomputed depth on ${both.length} of them: ratio p5/p50/p95 ${[0.05, 0.5, 0.95].map((q) => f(quantile(m, q), 3)).join('/')}`); }
  else console.log('  ocean.retreatDepth: n/a');
  checks.push([`the median asked depth where the floor binds is not below the package-off value (${f(quantile(on, 0.5))} against ${f(quantile(off, 0.5))} m)`, quantile(on, 0.5) >= quantile(off, 0.5) - 1e-9]);
  if (packaged && F === 50) checks.push([`gate: with the package the median asked depth rises by ${GATE[0]}–${GATE[1]} (ratio ${f(ratio, 2)})`, ratio >= GATE[0] && ratio <= GATE[1]]);
}
console.log(modelAsked.length ? `in-model asked depth over the day where under the floor (ocean.retreatDepth, ${modelAsked.length} column-steps): ${histogram(modelAsked)}` : 'in-model asked depth over the day: n/a (the ocean has no retreatDepth)');
const ss = sorted(stress.shareLow), sr = sorted(stress.ratioLow), st = sorted(stress.shareThin);
console.log(`stress depth (H = ${ekmanDepth > 0 && fieldOf(ocean.windDepth) ? 'H_m' : `${momentumDepth} m`}, classes of ${stressThickness} m or more): S < 0.8 H at ${f(100 * stress.low / stress.edges, 2)} % of sea edge-steps, ${f(100 * stress.lowTropics / stress.tropics, 2)} % at 20S–20N; there the mixed layer's share of the stress p10/p50/p90 ${[0.1, 0.5, 0.9].map((q) => f(quantile(ss, q), 2)).join('/')} (mean ${f(ss.reduce((a, b) => a + b, 0) / ss.length, 2)}), S/H median ${f(quantile(sr, 0.5), 2)}; under every layer thinner than H, mixed-layer share median ${f(quantile(st, 0.5), 2)}`);
console.log(`PEN over its bound: ${pen.stepsOver} cell-steps, at most ${pen.maxOver} in a step (iced ${pen.causes.iced}, entrained ${pen.causes.entrained}, binade ${pen.causes.binade}, other ${pen.causes.other}); worst ${f(pen.worst, 2)} of the bound${pen.worstAt ? ` at ${pen.worstAt}` : ''}`);
for (const [text, ok] of checks) console.log(`${ok ? 'PASS' : 'FAIL'} ${text}`);
await model.close();
process.exit(checks.every(([, ok]) => ok) ? 0 : 1);
