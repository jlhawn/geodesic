// How the page's cloud overlay blinks from one step to the next, and
// which part of the cloud and which process makes it. A saved state is
// loaded on the GPU as a continuing spin-up segment (OCEAN as for
// spinup.mjs, '{"everySteps":8}' unless given) and stepped STEPS times;
// after every step the frame's 'cloud' field (the column cloud water plus
// the deck's cover × liquid water path, as the page's Cloud cover overlay
// paints it) is read with the opacity the page gives it,
// a = 1 − exp(−g / 40) with g in g/m², and its parts per cell.
//
// The parts come from a copy of the adjust kernel with extra stores into
// the K1 buffer (free between a step's adjust and the next step's first
// Runge–Kutta stage); the stores read the state and write nothing the
// model reads. CHECK steps run twice with the plain kernel and once with
// the copy, and the frames are compared value by value (the copy's code
// may round a few values differently, by under 0.1 g/m²). Per cell and step:
// the column cloud water by layer group (below 800 hPa, 800–500, above
// 500) on entry to the adjust kernel (after the dynamics), after the
// boundary layer's mixing, after the saturation adjustment, after the
// plume convection, after the conversion to rain and the falling ice, and
// at the step's end; the conversion to rain and the ice falling out of
// each group; the deck's cover × LWP, cover, LWP, gate, carried height;
// the shallow cumulus' cover × water (which the radiation sees and the
// frame does not); the boundary layer's mixing top, depth and regime;
// the plume's deep and shallow firing; the largest relative humidity
// below 800 hPa; and the cell's per-step reflected and incident sunlight
// and outgoing longwave. A second pass from the same state repeats the
// steps with the profiles stored by layer of the CELLS worst-blinking
// cells and of the three cells whose blinks the deck most often drives.
//
// A blink is a cell whose opacity moves by more than JUMP (0.4) in one
// step and back by more than JUMP within the next three steps. Writes
// <out>/flicker_<tag>_N<N>.json (the numbers and the worst cells'
// histories) and <out>/flicker_<tag>_N<N>_plot.json, and draws the
// blink-frequency map and the strip of overlay frames with
// scripts/figures/cloudFlicker.py.
//   node scripts/cloudFlicker.mjs <state.bin> <outdir>
// Environment: STEPS (64), CELLS (10), JUMP (0.4), STRIDE (the page's
// steps per frame, max(2, round(384 / N))), CHECK (4; 0 skips), PYTHON.
import { writeFileSync, mkdirSync } from 'node:fs';
import { spawnSync } from 'node:child_process';
import { getDevice, readRanges } from '../js/gpu/device.module.js';
import { readState, tagOf, figureHeader, gpuModelFrom, dt as stepOf, DEG } from './figures/figureState.mjs';
import { DECK_OPEN, DECK_CLOSED } from '../js/physics/moist.module.js';

const [file, outDir] = process.argv.slice(2);
if (!file || !outDir) { console.error('usage: node scripts/cloudFlicker.mjs <state.bin> <outdir>'); process.exit(2); }
process.env.OCEAN ??= '{"everySteps":8}';
const STEPS = Number(process.env.STEPS ?? 64), CELLS = Number(process.env.CELLS ?? 10), JUMP = Number(process.env.JUMP ?? 0.4), CHECK = Number(process.env.CHECK ?? 4);
const OPACITY_SCALE = 40, WINDOW = 3, LOW = 80000, MID = 50000;
mkdirSync(outDir, { recursive: true });

const SLOT = {
  entry: 0, mixed: 3, saturated: 6, convected: 9, rained: 12, end: 15, conversion: 18, fall: 21,
  rain: 24, convectiveRain: 25, flags: 26, cumulus: 27, deck: 28, mixTop: 29, regime: 30, gate: 31, depth: 32, strat: 33,
  rhEntry: 34, rhEnd: 35, deckWater: 36, deckCover: 37, cumulusFlux: 38, deckHeight: 39, belowTop: 40, zeroed: 41, mixedLayers: 42, topMixed: 43, cumulusCover: 44, total: 45, entryBelowTop: 46, deckLayers: 47,
};
const NSLOT = 48, ROWS = 15;
const ROW = { qc: 0, conversion: 6, fall: 7, rhEntry: 8, rhEnd: 9, qtEntry: 10, mix: 11, belowTop: 12, temperature: 13, pressure: 14 };
const STAGES = ['transport', 'mixing', 'saturation', 'plume', 'rainout', 'filler'];
const GROUPS = ['low', 'mid', 'high'];

/*
 * The adjust kernel's text with the stores; `cells` are the cells whose
 * profiles are stored by layer (none: [-1]).
 */
function instrumented(code, C, cells) {
  const base = NSLOT * C, n = cells.length;
  const helpers = `
const FLK_BASE: i32 = ${base}; const FLK_ROWS: i32 = ${ROWS}; const FLK_N: i32 = ${n};
const FLK_CELLS = array<i32, ${n}>(${cells.join(', ')});
var<private> flkConv: vec3<f32>;
var<private> flkFall: vec3<f32>;
fn flkGroup(pi: f32, k: i32) -> i32 { let p = pi * LV[L_SM + k]; return select(select(2, 1, p > ${MID}.0), 0, p > ${LOW}.0); }
fn flkSlot(i: i32) -> i32 { for (var s = 0; s < FLK_N; s++) { if (FLK_CELLS[s] == i) { return s; } } return -1; }
fn flkLayer(i: i32, row: i32, k: i32, value: f32) { let s = flkSlot(i); if (s >= 0) { OUT[FLK_BASE + (s * FLK_ROWS + row) * K + k] = value; } }
fn flkRH(i: i32, pi: f32, k: i32) -> f32 { let idx = k * C + i; return IN[S_Q + idx] / cloudSat(IN[S_TH + idx] * D[D_EXM + idx], pi * LV[L_SM + k]).x; }
fn flkLowRH(i: i32, pi: f32) -> f32 { var r = 0.0; for (var k = 0; k < K; k++) { if (pi * LV[L_SM + k] > ${LOW}.0) { r = max(r, flkRH(i, pi, k)); } } return r; }
fn flkStage(i: i32, pi: f32, st: i32) {
  var s = vec3<f32>(0.0, 0.0, 0.0);
  let slot = flkSlot(i);
  for (var k = 0; k < K; k++) {
    let idx = k * C + i;
    s[flkGroup(pi, k)] += pi * LV[L_DS + k] / GRAV * IN[S_QC + idx];
    if (slot >= 0) { OUT[FLK_BASE + (slot * FLK_ROWS + st) * K + k] = IN[S_QC + idx]; }
  }
  OUT[(3 * st) * C + i] = s.x; OUT[(3 * st + 1) * C + i] = s.y; OUT[(3 * st + 2) * C + i] = s.z;
}
fn flkEntry(i: i32, pi: f32) {
  flkConv = vec3<f32>(0.0, 0.0, 0.0); flkFall = vec3<f32>(0.0, 0.0, 0.0);
  flkStage(i, pi, 0);
  var zeroed = 0.0; var count = 0.0; var topMixed = -1; var inside = 0.0;
  for (var k = 0; k < K; k++) {
    let idx = k * C + i;
    if ((D[D_GEO + idx] + LV[L_GABS + k]) / GRAV < PH[PH_MIXTOP + i]) { inside += pi * LV[L_DS + k] / GRAV * IN[S_QC + idx]; }
  }
  OUT[${SLOT.entryBelowTop} * C + i] = inside;
  for (var k = KTOP; k < K; k++) {
    let idx = k * C + i;
    let mixed = (k > KTOP && PH[PH_MIX + (k - 1) * C + i] > 0.0) || (k < K - 1 && PH[PH_MIX + k * C + i] > 0.0);
    if (mixed) { zeroed += pi * LV[L_DS + k] / GRAV * IN[S_QC + idx]; count += 1.0; if (topMixed < 0) { topMixed = k; } }
  }
  OUT[${SLOT.zeroed} * C + i] = zeroed; OUT[${SLOT.mixedLayers} * C + i] = count; OUT[${SLOT.topMixed} * C + i] = f32(topMixed);
  OUT[${SLOT.rhEntry} * C + i] = flkLowRH(i, pi);
  if (flkSlot(i) >= 0) {
    for (var k = 0; k < K; k++) {
      let idx = k * C + i;
      flkLayer(i, ${ROW.rhEntry}, k, flkRH(i, pi, k));
      flkLayer(i, ${ROW.qtEntry}, k, IN[S_Q + idx] + IN[S_QC + idx]);
      flkLayer(i, ${ROW.mix}, k, PH[PH_MIX + k * C + i]);
      flkLayer(i, ${ROW.belowTop}, k, select(0.0, 1.0, (D[D_GEO + idx] + LV[L_GABS + k]) / GRAV < PH[PH_MIXTOP + i]));
      flkLayer(i, ${ROW.pressure}, k, pi * LV[L_SM + k]);
      flkLayer(i, ${ROW.conversion}, k, 0.0); flkLayer(i, ${ROW.fall}, k, 0.0);
    }
  }
}
fn flkExit(i: i32, pi: f32, rained: f32, convected: f32, mixes: bool) {
  flkStage(i, pi, 5);
  let zb = (D[D_GEO + (K - 1) * C + i] + LV[L_GABS + K - 1]) / GRAV;
  OUT[${SLOT.conversion} * C + i] = flkConv.x; OUT[${SLOT.conversion + 1} * C + i] = flkConv.y; OUT[${SLOT.conversion + 2} * C + i] = flkConv.z;
  OUT[${SLOT.fall} * C + i] = flkFall.x; OUT[${SLOT.fall + 1} * C + i] = flkFall.y; OUT[${SLOT.fall + 2} * C + i] = flkFall.z;
  OUT[${SLOT.rain} * C + i] = rained; OUT[${SLOT.convectiveRain} * C + i] = convected;
  OUT[${SLOT.flags} * C + i] = select(0.0, 1.0, cuDeep) + select(0.0, 2.0, PH[PH_CUMF + i] > 0.0) + select(0.0, 4.0, mixes);
  var cumulus = 0.0; var cumulusCover = 0.0;
  for (var k = CU_K0; k < K; k++) {
    let slot = (k - CU_K0) * C + i;
    cumulus += pi * LV[L_DS + k] / GRAV * PH[PH_CUCOVER + slot] * PH[PH_CUWATER + slot];
    cumulusCover = max(cumulusCover, PH[PH_CUCOVER + slot]);
  }
  OUT[${SLOT.cumulus} * C + i] = cumulus; OUT[${SLOT.cumulusCover} * C + i] = cumulusCover;
  OUT[${SLOT.deck} * C + i] = PH[PH_DECKF + i] * PH[PH_DECK + i];
  OUT[${SLOT.deckWater} * C + i] = PH[PH_DECK + i]; OUT[${SLOT.deckCover} * C + i] = PH[PH_DECKF + i];
  OUT[${SLOT.mixTop} * C + i] = PH[PH_MIXTOP + i] - zb; OUT[${SLOT.depth} * C + i] = PH[PH_DEPTH + i] - zb;
  OUT[${SLOT.regime} * C + i] = PH[PH_REGIME + i]; OUT[${SLOT.gate} * C + i] = PH[PH_MLMGATE + i];
  OUT[${SLOT.strat} * C + i] = PH[PH_STRAT + i]; OUT[${SLOT.cumulusFlux} * C + i] = PH[PH_CUMF + i]; OUT[${SLOT.deckHeight} * C + i] = PH[PH_MLMH + i];
  OUT[${SLOT.rhEnd} * C + i] = flkLowRH(i, pi);
  var below = 0.0; var column = 0.0; var deckLayers = 0.0;
  for (var k = 0; k < K; k++) {
    let idx = k * C + i; let m = pi * LV[L_DS + k] / GRAV * IN[S_QC + idx];
    column += m;
    if (D[D_GEO + idx] + LV[L_GABS + k] < GRAV * PH[PH_MLMH + i]) { deckLayers += 1.0; }
    if ((D[D_GEO + idx] + LV[L_GABS + k]) / GRAV < PH[PH_MIXTOP + i]) { below += m; }
    flkLayer(i, ${ROW.rhEnd}, k, flkRH(i, pi, k));
    flkLayer(i, ${ROW.temperature}, k, IN[S_TH + idx] * D[D_EXM + idx]);
  }
  OUT[${SLOT.belowTop} * C + i] = below; OUT[${SLOT.deckLayers} * C + i] = deckLayers;
  OUT[${SLOT.total} * C + i] = column + PH[PH_DECKF + i] * PH[PH_DECK + i];
}
`;
  const swaps = [
    ['  if (mixes) {\n    if (MOIST_BL)', '  flkEntry(i, pi);\n  if (mixes) {\n    if (MOIST_BL)'],
    ['  diagnoseColumn(i);\n  saturateColumn(i, pi);\n  let produced = plumeColumn(i, pi, dt);\n  if (PH[PH_CUMF + i] > 0.0) { saturateColumn(i, pi); }\n',
      '  diagnoseColumn(i);\n  flkStage(i, pi, 1);\n  saturateColumn(i, pi);\n  flkStage(i, pi, 2);\n  let produced = plumeColumn(i, pi, dt);\n  if (PH[PH_CUMF + i] > 0.0) { saturateColumn(i, pi); }\n  flkStage(i, pi, 3);\n'],
    ['          descending = leaving * mass;\n', '          descending = leaving * mass;\n          flkFall[flkGroup(pi, k)] += leaving * mass; flkLayer(i, 7, k, leaving * mass);\n'],
    ['    IN[S_QC + idx] = qc - converted;\n', '    IN[S_QC + idx] = qc - converted;\n    flkConv[flkGroup(pi, k)] += mass * converted; flkLayer(i, 6, k, mass * converted);\n'],
    ['  if (moved) { saturateColumn(i, pi); }\n', '  if (moved) { saturateColumn(i, pi); }\n  flkStage(i, pi, 4);\n'],
  ];
  let out = code;
  for (const [from, to] of swaps) {
    if (out.split(from).length !== 2) throw new Error(`the adjust kernel no longer holds exactly one ${JSON.stringify(from.slice(0, 40))}`);
    out = out.replace(from, to);
  }
  const main = out.lastIndexOf('@compute');
  out = out.slice(0, main) + helpers + out.slice(main);
  const end = out.lastIndexOf('}');
  return `${out.slice(0, end)}  flkExit(i, pi, rained, convected, mixes);\n}${out.slice(end + 1)}`;
}

const saved = await readState(file);
const header = figureHeader(file, saved);
const { device } = await getDevice();
const plainModule = device.createShaderModule.bind(device);
const probe = { on: false, cells: [-1], C: 10 * saved.N * saved.N + 2 };
device.createShaderModule = (descriptor) => (probe.on && descriptor.label === 'adjust' ? plainModule({ ...descriptor, code: instrumented(descriptor.code, probe.C, probe.cells) }) : plainModule(descriptor));

const N = saved.N, dt = stepOf(N), STRIDE = Number(process.env.STRIDE ?? Math.max(2, Math.round(24 * 16 / N)));
const toc = (t0) => ((performance.now() - t0) / 1000).toFixed(1);

async function run(on, cells, steps, record) {
  probe.on = on; probe.cells = cells.length ? cells : [-1];
  const t0 = performance.now();
  const model = await gpuModelFrom(saved);
  const { gpu } = model, C = gpu.C, PH = gpu.layout.PH;
  if (on && gpu.layout.S.total < NSLOT * C + probe.cells.length * ROWS * gpu.K) throw new Error('K1 is too small for the stores');
  const loaded = toc(t0);
  for (let n = 1; n <= steps; n++) {
    await model.step(dt);
    const frame = model.beginFrame({ fields: ['cloud'] });
    const extra = on ? readRanges(device, gpu.buffers.K1, [{ offset: 0, length: NSLOT * C }, ...(cells.length ? [{ offset: NSLOT * C, length: cells.length * ROWS * gpu.K }] : [])]) : null;
    const flux = on ? readRanges(device, gpu.buffers.PH, ['REFL', 'INS', 'OLR'].map((name) => ({ offset: PH[name], length: C }))) : null;
    const [f, e, x] = await Promise.all([frame, extra, flux]);
    record(n, { cloud: Float32Array.from(f.fields.cloud), diag: e ? Float32Array.from(e[0]) : null, profiles: e && e[1] ? Float32Array.from(e[1]) : null, flux: x ? x.map((v) => Float32Array.from(v)) : null });
  }
  await model.settle();
  console.log(`${on ? 'probed' : 'plain'} run: ${steps} steps, loaded in ${loaded} s, done in ${toc(t0)} s`);
  return model;
}

let check = null;
if (CHECK > 0) {
  const plain = [], again = [], probed = [];
  for (const [on, list] of [[false, plain], [false, again], [true, probed]]) (await run(on, [], CHECK, (n, r) => list.push(r.cloud))).destroy();
  const compare = (x, y) => {
    let differing = 0, largest = 0;
    x.forEach((c, n) => { for (let i = 0; i < c.length; i++) if (c[i] !== y[n][i]) { differing++; largest = Math.max(largest, Math.abs(c[i] - y[n][i])); } });
    return { differingValues: differing, largestKgPerM2: largest };
  };
  check = { steps: CHECK, plainAgainstPlain: compare(plain, again), plainAgainstProbed: compare(plain, probed) };
  console.log(`probe check over ${CHECK} steps: ${JSON.stringify(check)}`);
}

const frames = [], diags = [], fluxes = [];
let model = await run(true, [], STEPS, (n, r) => { frames.push(r.cloud); diags.push(r.diag); fluxes.push(r.flux); });
const { mesh, geography } = model, C = mesh.nCells, K = model.gpu.K;
if (C !== probe.C) throw new Error(`the mesh has ${C} cells, not ${probe.C}`);
model.destroy();

const lat = Float64Array.from(mesh.latCell, (v) => v * DEG), lon = Float64Array.from(mesh.lonCell, (v) => { const d = v * DEG; return d > 180 ? d - 360 : d; });
const land = Uint8Array.from({ length: C }, (_, i) => (geography && geography.land[i] ? 1 : 0));
const inBox = (i, [s, n, w, e]) => lat[i] >= s && lat[i] <= n && (w <= e ? lon[i] >= w && lon[i] <= e : lon[i] >= w || lon[i] <= e);
const REGIONS = {
  global: () => true, sea: (i) => !land[i], land: (i) => !!land[i], tropics: (i) => Math.abs(lat[i]) <= 30,
  sePacific: (i) => inBox(i, [-30, -10, -110, -80]), peru: (i) => inBox(i, [-20, -10, -90, -80]), namibia: (i) => inBox(i, [-20, -10, 0, 10]),
  california: (i) => inBox(i, [20, 30, -130, -120]), nPacific: (i) => inBox(i, [35, 50, 160, -160]),
  northStormTrack: (i) => lat[i] >= 40 && lat[i] <= 60, southStormTrack: (i) => lat[i] >= -60 && lat[i] <= -40, itcz: (i) => inBox(i, [5, 12, 160, -100]),
};
const regionCells = Object.fromEntries(Object.entries(REGIONS).map(([name, test]) => [name, Array.from({ length: C }, (_, i) => i).filter(test)]));

const T = frames.length;
const opacity = frames.map((g) => Float32Array.from(g, (v) => 1 - Math.exp(-Math.max(0, 1000 * v) / OPACITY_SCALE)));
const at = (t, slot, i) => diags[t][slot * C + i];

/* Blinks in a sequence of opacity fields: onsets with their return transitions. */
function blinksOf(seq) {
  const n = seq.length, onsets = [], plain = new Float64Array(n - 1), counted = Math.max(0, n - 1 - WINDOW);
  for (let t = 0; t < n - 1; t++) {
    const a0 = seq[t], a1 = seq[t + 1];
    let jumps = 0;
    for (let i = 0; i < C; i++) {
      const d = a1[i] - a0[i];
      if (Math.abs(d) <= JUMP) continue;
      jumps++;
      if (t >= counted) continue;
      const sign = Math.sign(d);
      for (let j = 1; j <= WINDOW; j++) {
        if ((seq[t + 1 + j][i] - a1[i]) * sign < -JUMP) {
          let back = t + 1, largest = 0;
          for (let s = t + 1; s < t + 1 + j; s++) { const step = (seq[s + 1][i] - seq[s][i]) * -sign; if (step > largest) { largest = step; back = s; } }
          onsets.push({ t, i, sign, j, back });
          break;
        }
      }
    }
    plain[t] = jumps / C;
  }
  return { onsets, plain, counted };
}

const steps = blinksOf(opacity);
const pageFrames = opacity.filter((_, t) => t % STRIDE === STRIDE - 1);
const page = blinksOf(pageFrames);

const share = (count, cells, transitions) => (cells && transitions ? count / (cells * transitions) : null);
const round = (v, n = 4) => (v === null || !Number.isFinite(v) ? null : +v.toFixed(n));
const blinkCount = new Uint16Array(C), plainCount = new Uint16Array(C), magnitude = new Float64Array(C);
for (const b of steps.onsets) { blinkCount[b.i]++; magnitude[b.i] += Math.abs(opacity[b.t + 1][b.i] - opacity[b.t][b.i]); }
for (let t = 0; t < T - 1; t++) for (let i = 0; i < C; i++) if (Math.abs(opacity[t + 1][i] - opacity[t][i]) > JUMP) plainCount[i]++;
const pageCount = new Uint16Array(C);
for (const b of page.onsets) pageCount[b.i]++;

/* The parts' changes over the transition t → t + 1, and the step's stages for one resolved group. */
const PARTS = ['low', 'mid', 'high', 'deck'];
function partChanges(t, i) {
  const d = PARTS.map((_, p) => (p < 3 ? at(t + 1, SLOT.end + p, i) - at(t, SLOT.end + p, i) : at(t + 1, SLOT.deck, i) - at(t, SLOT.deck, i)));
  return { d, cumulus: at(t + 1, SLOT.cumulus, i) - at(t, SLOT.cumulus, i) };
}
function stagesOf(t, i, g) {
  const s = (st) => at(t + 1, 3 * st + g, i);
  return [s(0) - at(t, SLOT.end + g, i), s(1) - s(0), s(2) - s(1), s(3) - s(2), s(4) - s(3), s(5) - s(4)];
}

function attribute(onsets) {
  const tally = () => ({ count: 0, driver: Object.fromEntries(PARTS.map((p) => [p, 0])), cumulusLarger: 0 });
  const byRegion = Object.fromEntries(Object.keys(REGIONS).map((r) => [r, tally()]));
  const transitions = {
    appear: { count: 0, stage: Object.fromEntries(STAGES.map((s) => [s, 0])), mixPlusSaturation: 0, cloudAboveMixTop: 0, previousCloudInsideMixTop: 0, mixTops: [], previousMixTops: [] },
    vanish: { count: 0, stage: Object.fromEntries(STAGES.map((s) => [s, 0])), mixPlusSaturation: 0, conversionOverHalf: 0, mixingZeroedHalf: 0, mixingZeroedNotRestored: 0, iceFallOverHalf: 0, cloudInsideNewMixTop: 0, previousCloudAboveMixTop: 0, schemeSwitch: 0, mixTops: [], previousMixTops: [] },
  };
  const deck = { transitions: 0, gateCrossedUndecided: 0, gateBelowOpenEitherSide: 0, gateAboveClosedBothSides: 0, regimeChanged: 0, heightWithin5m: 0, carriedLayersChanged: 0, layersChanged: 0, coverUp: 0, waterUp: 0, gates: [] };
  const context = { transitions: 0, regimeChanged: 0, regimeAlternates: 0, topMixedChanged: 0, mixTopJumpOver200m: 0, deepPlumeToggled: 0, shallowCumulusToggled: 0, mixedToggled: 0, driverBelowMixTop: 0, regimePairs: {} };
  const radiation = { albedo: [], olr: [] };
  const level = { lowestTwo: 0 };
  const isIn = Object.fromEntries(Object.entries(regionCells).map(([r, cells]) => { const m = new Uint8Array(C); for (const i of cells) m[i] = 1; return [r, m]; }));
  function transition(t, i) {
    const { d, cumulus } = partChanges(t, i);
    let p = 0;
    for (let q = 1; q < 4; q++) if (Math.abs(d[q]) > Math.abs(d[p])) p = q;
    context.transitions++;
    const r0 = at(t, SLOT.regime, i), r1 = at(t + 1, SLOT.regime, i);
    if (r0 !== r1) context.regimeChanged++;
    const key = `${r0}->${r1}`;
    context.regimePairs[key] = (context.regimePairs[key] ?? 0) + 1;
    if (t >= 1 && at(t - 1, SLOT.regime, i) === r1 && r0 !== r1) context.regimeAlternates++;
    if (at(t, SLOT.topMixed, i) !== at(t + 1, SLOT.topMixed, i)) context.topMixedChanged++;
    if (Math.abs(at(t + 1, SLOT.mixTop, i) - at(t, SLOT.mixTop, i)) > 200) context.mixTopJumpOver200m++;
    const f0 = at(t, SLOT.flags, i), f1 = at(t + 1, SLOT.flags, i);
    if ((f0 & 1) !== (f1 & 1)) context.deepPlumeToggled++;
    if ((f0 & 2) !== (f1 & 2)) context.shallowCumulusToggled++;
    if ((f0 & 4) !== (f1 & 4)) context.mixedToggled++;
    const below = at(t + 1, SLOT.belowTop, i) - at(t, SLOT.belowTop, i);
    if (p < 3 && Math.abs(below) > 0.5 * Math.abs(d[p])) context.driverBelowMixTop++;
    if (t + 2 < T) {
      const [r1x, i1x, o1x] = fluxes[t + 1].map((v) => v[i]), [r2x, i2x, o2x] = fluxes[t + 2] ? fluxes[t + 2].map((v) => v[i]) : [NaN, NaN, NaN];
      if (i1x > 50 && i2x > 50) radiation.albedo.push(Math.abs(r2x / i2x - r1x / i1x));
      if (Number.isFinite(o2x)) radiation.olr.push(Math.abs(o2x - o1x));
    }
    if (p === 3) {
      deck.transitions++;
      const g0 = at(t, SLOT.gate, i), g1 = at(t + 1, SLOT.gate, i);
      if ((g0 > 0.5) !== (g1 > 0.5)) deck.gateCrossedUndecided++;
      if (Math.min(g0, g1) < DECK_OPEN) deck.gateBelowOpenEitherSide++;
      if (g0 > DECK_CLOSED && g1 > DECK_CLOSED) deck.gateAboveClosedBothSides++;
      if (r0 !== r1) deck.regimeChanged++;
      if (Math.abs(at(t + 1, SLOT.deckHeight, i) - at(t, SLOT.deckHeight, i)) < 5) deck.heightWithin5m++;
      if (t >= 1 && at(t, SLOT.deckLayers, i) !== at(t - 1, SLOT.deckLayers, i)) deck.carriedLayersChanged++;
      if (at(t + 1, SLOT.deckLayers, i) !== at(t, SLOT.deckLayers, i)) deck.layersChanged++;
      if (Math.sign(at(t + 1, SLOT.deckCover, i) - at(t, SLOT.deckCover, i)) === Math.sign(d[3])) deck.coverUp++;
      if (Math.sign(at(t + 1, SLOT.deckWater, i) - at(t, SLOT.deckWater, i)) === Math.sign(d[3])) deck.waterUp++;
      if (deck.gates.length < 2000) deck.gates.push([round(g0, 3), round(g1, 3)]);
      return { p, d, cumulus };
    }
    const st = stagesOf(t, i, p), sign = Math.sign(d[p]), kind = sign > 0 ? 'appear' : 'vanish', bucket = transitions[kind];
    let best = 0;
    for (let s = 1; s < STAGES.length; s++) if (st[s] * sign > st[best] * sign) best = s;
    bucket.count++; bucket.stage[STAGES[best]]++;
    const sum = (u, slot) => at(u, slot, i) + at(u, slot + 1, i) + at(u, slot + 2, i);
    const entry = sum(t + 1, SLOT.entry), entryInside = at(t + 1, SLOT.entryBelowTop, i);
    const before = sum(t, SLOT.end), beforeInside = at(t, SLOT.belowTop, i), after = sum(t + 1, SLOT.end), afterInside = at(t + 1, SLOT.belowTop, i);
    bucket.mixTops.push(at(t + 1, SLOT.mixTop, i)); bucket.previousMixTops.push(at(t, SLOT.mixTop, i));
    if (kind === 'appear') {
      if (after > 0 && after - afterInside > 0.5 * after) bucket.cloudAboveMixTop++;
      if (before > 0 && beforeInside > 0.5 * before) bucket.previousCloudInsideMixTop++;
    } else {
      const inside = entry > 0 && entryInside > 0.5 * entry, above = before > 0 && before - beforeInside > 0.5 * before;
      if (inside) bucket.cloudInsideNewMixTop++;
      if (above) bucket.previousCloudAboveMixTop++;
      if (inside && above) bucket.schemeSwitch++;
    }
    const others = [st[0], st[3], st[4], st[5]].map((v) => v * sign);
    if ((st[1] + st[2]) * sign > Math.max(...others)) bucket.mixPlusSaturation++;
    if (kind === 'vanish') {
      const before = at(t + 1, SLOT.convected + p, i), entry = at(t + 1, SLOT.entry + p, i);
      if (before > 0 && at(t + 1, SLOT.conversion + p, i) > 0.5 * before) bucket.conversionOverHalf++;
      if (before > 0 && at(t + 1, SLOT.fall + p, i) > 0.5 * before) bucket.iceFallOverHalf++;
      if (entry > 0 && st[1] < -0.5 * entry) { bucket.mixingZeroedHalf++; if (st[1] + st[2] < -0.5 * entry) bucket.mixingZeroedNotRestored++; }
    }
    return { p, d, cumulus };
  }
  for (const b of onsets) {
    const { p, d, cumulus } = transition(b.t, b.i);
    transition(b.back, b.i);
    for (const [r, m] of Object.entries(isIn)) if (m[b.i]) { const x = byRegion[r]; x.count++; x.driver[PARTS[p]]++; if (Math.abs(cumulus) > Math.abs(d[p])) x.cumulusLarger++; }
  }
  const median = (v) => { if (!v.length) return null; const s = [...v].sort((x, y) => x - y); return s[s.length >> 1]; };
  for (const bucket of Object.values(transitions)) {
    const quartiles = (v) => { const s = [...v].sort((x, y) => x - y); return s.length ? [0.25, 0.5, 0.75].map((q) => Math.round(s[Math.floor(q * (s.length - 1))])) : null; };
    bucket.mixTopQuartiles = quartiles(bucket.mixTops); bucket.previousMixTopQuartiles = quartiles(bucket.previousMixTops);
    delete bucket.mixTops; delete bucket.previousMixTops;
  }
  return { byRegion, transitions, deck: { ...deck, gates: deck.gates.slice(0, 40) }, context, radiation: { medianAlbedoChange: round(median(radiation.albedo), 3), medianOlrChange: round(median(radiation.olr), 1), albedoSamples: radiation.albedo.length } };
}
const attribution = attribute(steps.onsets);

/* Radiation's step-to-step change in cells that do not blink, for comparison. */
const quiet = (() => {
  const albedo = [], olr = [];
  for (let t = 1; t + 1 < T && albedo.length < 200000; t += 3) for (let i = 0; i < C; i += 7) {
    if (blinkCount[i]) continue;
    const [r1, i1, o1] = fluxes[t].map((v) => v[i]), [r2, i2, o2] = fluxes[t + 1].map((v) => v[i]);
    if (i1 > 50 && i2 > 50) albedo.push(Math.abs(r2 / i2 - r1 / i1));
    olr.push(Math.abs(o2 - o1));
  }
  const median = (v) => { const s = [...v].sort((x, y) => x - y); return s[s.length >> 1]; };
  return { medianAlbedoChange: round(median(albedo), 3), medianOlrChange: round(median(olr), 1) };
})();

/* Lag-1 autocorrelation of the opacity's anomalies about each cell's mean, and of its step changes, pooled. */
function lagOne(cells, seq = opacity) {
  let n1 = 0, d1 = 0, n2 = 0, d2 = 0;
  const n = seq.length;
  for (const i of cells) {
    let mean = 0;
    for (let t = 0; t < n; t++) mean += seq[t][i];
    mean /= n;
    for (let t = 0; t < n; t++) { const x = seq[t][i] - mean; d1 += x * x; if (t + 1 < n) n1 += x * (seq[t + 1][i] - mean); }
    let dm = 0;
    for (let t = 0; t + 1 < n; t++) dm += seq[t + 1][i] - seq[t][i];
    dm /= n - 1;
    for (let t = 0; t + 1 < n; t++) { const x = seq[t + 1][i] - seq[t][i] - dm; d2 += x * x; if (t + 2 < n) n2 += x * (seq[t + 2][i] - seq[t + 1][i] - dm); }
  }
  return { anomaly: round(n1 / d1, 3), change: round(n2 / d2, 3) };
}
const allCells = regionCells.global, blinking = allCells.filter((i) => blinkCount[i] > 0);
const autocorrelation = { all: lagOne(allCells), blinkingCells: lagOne(blinking), pageCadenceAll: lagOne(allCells, pageFrames), pageCadenceBlinking: lagOne(blinking, pageFrames) };

/* How many steps a blinking cell stays on (from an up-jump to the next down-jump) or off. */
function runLengths(cells) {
  const on = {}, off = {};
  let censoredOn = 0, censoredOff = 0;
  for (const i of cells) {
    let state = null, since = 0;
    for (let t = 0; t + 1 < T; t++) {
      const d = opacity[t + 1][i] - opacity[t][i];
      if (Math.abs(d) <= JUMP) continue;
      const now = d > 0 ? 'on' : 'off';
      if (state !== null && state !== now) { const tally = state === 'on' ? on : off; const len = t + 1 - since; tally[len] = (tally[len] ?? 0) + 1; }
      state = now; since = t + 1;
    }
    if (state === 'on') censoredOn++; else if (state === 'off') censoredOff++;
  }
  return { on, off, censoredOn, censoredOff };
}
const runs = runLengths(blinking);

/* The worst cells, at least two degrees apart. */
const order = Array.from({ length: C }, (_, i) => i).filter((i) => blinkCount[i] > 0).sort((x, y) => blinkCount[y] - blinkCount[x] || magnitude[y] - magnitude[x]);
const worst = [];
for (const i of order) {
  if (worst.length >= CELLS) break;
  if (worst.every((j) => Math.acos(Math.min(1, Math.sin(mesh.latCell[i]) * Math.sin(mesh.latCell[j]) + Math.cos(mesh.latCell[i]) * Math.cos(mesh.latCell[j]) * Math.cos(mesh.lonCell[i] - mesh.lonCell[j]))) * DEG > 2)) worst.push(i);
}

const deckDriven = new Uint16Array(C);
for (const b of steps.onsets) { const { d } = partChanges(b.t, b.i); if (Math.abs(d[3]) >= Math.max(Math.abs(d[0]), Math.abs(d[1]), Math.abs(d[2]))) deckDriven[b.i]++; }
const deckOrder = Array.from({ length: C }, (_, i) => i).filter((i) => deckDriven[i] > 0).sort((x, y) => deckDriven[y] - deckDriven[x] || magnitude[y] - magnitude[x]);
const deckWorst = [];
for (const i of deckOrder) {
  if (deckWorst.length >= 3) break;
  if ([...worst, ...deckWorst].every((j) => Math.acos(Math.min(1, Math.sin(mesh.latCell[i]) * Math.sin(mesh.latCell[j]) + Math.cos(mesh.latCell[i]) * Math.cos(mesh.latCell[j]) * Math.cos(mesh.lonCell[i] - mesh.lonCell[j]))) * DEG > 2)) deckWorst.push(i);
}
const traced = [...worst, ...deckWorst];
const profiles = [];
let replayDiffers = 0;
model = await run(true, traced, STEPS, (n, r) => {
  profiles.push(r.profiles);
  for (let i = 0; i < C; i++) if (r.cloud[i] !== frames[n - 1][i]) replayDiffers++;
});
model.destroy();

const history = traced.map((i, s) => {
  const series = (slot, n = 4) => diags.map((d) => round(d[slot * C + i], n));
  const layers = (row, t) => Array.from(profiles[t].subarray((s * ROWS + row) * K, (s * ROWS + row + 1) * K));
  const pressure = layers(ROW.pressure, 0);
  const shown = pressure.map((p, k) => k).filter((k) => pressure[k] > 55000);
  const qcAt = (t, st) => shown.map((k) => round(1000 * layers(ROW.qc + st, t)[k], 3));
  return {
    cell: i, lat: round(lat[i], 2), lon: round(lon[i], 2), land: !!land[i], blinks: blinkCount[i], deckDrivenBlinks: deckDriven[i], jumps: plainCount[i], chosenFor: s < worst.length ? 'blinks' : 'deck-driven blinks',
    regions: Object.keys(REGIONS).filter((r) => r !== 'global' && REGIONS[r](i)),
    layersHPa: shown.map((k) => Math.round(pressure[k] / 100)),
    opacity: opacity.map((a) => round(a[i], 3)), cloudGrams: frames.map((g) => round(1000 * g[i], 2)),
    lowGrams: diags.map((d) => round(1000 * d[(SLOT.end) * C + i], 2)), midGrams: diags.map((d) => round(1000 * d[(SLOT.end + 1) * C + i], 2)), highGrams: diags.map((d) => round(1000 * d[(SLOT.end + 2) * C + i], 2)),
    deckGrams: diags.map((d) => round(1000 * d[SLOT.deck * C + i], 2)), deckLayers: series(SLOT.deckLayers, 0), deckCover: series(SLOT.deckCover, 3), deckLwpGrams: diags.map((d) => round(1000 * d[SLOT.deckWater * C + i], 1)), gate: series(SLOT.gate, 3), deckHeight: series(SLOT.deckHeight, 0),
    cumulusGrams: diags.map((d) => round(1000 * d[SLOT.cumulus * C + i], 2)), cumulusCover: series(SLOT.cumulusCover, 3),
    regime: series(SLOT.regime, 0), mixTop: series(SLOT.mixTop, 0), entryInsideMixTopGrams: diags.map((d) => round(1000 * d[SLOT.entryBelowTop * C + i], 2)), endInsideMixTopGrams: diags.map((d) => round(1000 * d[SLOT.belowTop * C + i], 2)), depth: series(SLOT.depth, 0), topMixedLayer: series(SLOT.topMixed, 0), mixedLayers: series(SLOT.mixedLayers, 0),
    deepPlume: diags.map((d) => d[SLOT.flags * C + i] & 1), shallowCumulus: diags.map((d) => (d[SLOT.flags * C + i] & 2) >> 1),
    lowRhEntry: series(SLOT.rhEntry, 3), lowRhEnd: series(SLOT.rhEnd, 3),
    conversionGrams: diags.map((d) => [0, 1, 2].map((g) => round(1000 * d[(SLOT.conversion + g) * C + i], 2))),
    rainMm: diags.map((d) => round(d[SLOT.rain * C + i], 4)), convectiveMm: diags.map((d) => round(d[SLOT.convectiveRain * C + i], 4)),
    stagesGrams: diags.map((d, t) => (t === 0 ? null : [0, 1, 2].map((g) => [0, 1, 2, 3, 4, 5].map((st) => round(1000 * (d[(3 * st + g) * C + i] - (st === 0 ? diags[t - 1][(SLOT.end + g) * C + i] : d[(3 * (st - 1) + g) * C + i])), 2))))),
    qcEndMilligramsPerKg: profiles.map((_, t) => qcAt(t, 5)),
    qcEntryMilligramsPerKg: profiles.map((_, t) => qcAt(t, 0)),
    qcAfterMixingMilligramsPerKg: profiles.map((_, t) => qcAt(t, 1)),
    qcAfterSaturationMilligramsPerKg: profiles.map((_, t) => qcAt(t, 2)),
    rhEntry: profiles.map((_, t) => shown.map((k) => round(layers(ROW.rhEntry, t)[k], 3))),
    rhEnd: profiles.map((_, t) => shown.map((k) => round(layers(ROW.rhEnd, t)[k], 3))),
    belowMixTop: profiles.map((_, t) => shown.map((k) => layers(ROW.belowTop, t)[k])),
    mixing: profiles.map((_, t) => shown.map((k) => round(layers(ROW.mix, t)[k], 3))),
    conversionByLayerGrams: profiles.map((_, t) => shown.map((k) => round(1000 * layers(ROW.conversion, t)[k], 3))),
  };
});

let worstTotal = 0;
for (let t = 0; t < T; t++) for (let i = 0; i < C; i++) worstTotal = Math.max(worstTotal, Math.abs(diags[t][SLOT.total * C + i] - frames[t][i]));

const regions = Object.fromEntries(Object.entries(regionCells).map(([name, cells]) => {
  const blinks = attribution.byRegion[name];
  let plainJumps = 0, pageBlinks = 0;
  for (const i of cells) { plainJumps += plainCount[i]; pageBlinks += pageCount[i]; }
  return [name, {
    cells: cells.length, blinkShare: round(share(blinks.count, cells.length, steps.counted), 5), jumpShare: round(share(plainJumps, cells.length, T - 1), 5),
    pageBlinkShare: round(share(pageBlinks, cells.length, page.counted), 5),
    driver: Object.fromEntries(Object.entries(blinks.driver).map(([k, v]) => [k, blinks.count ? round(v / blinks.count, 3) : null])),
    cumulusChangeLarger: blinks.count ? round(blinks.cumulusLarger / blinks.count, 3) : null,
  }];
}));

const summary = {
  ...header, steps: STEPS, dt, stride: STRIDE, jump: JUMP, window: WINDOW, cells: C, check, replayDiffers, probeTotalMaxError: worstTotal,
  blinkOnsets: steps.onsets.length, countedTransitions: steps.counted, pageFrames: pageFrames.length, pageCounted: page.counted,
  perStepBlinkShare: steps.plain.length ? round(steps.onsets.length / (C * steps.counted), 5) : null,
  perStepJumpShare: Array.from(steps.plain, (v) => round(v, 5)),
  regions, attribution, quietRadiation: quiet, autocorrelation, runs,
  cellsBlinking: blinking.length, cellsBlinkingShare: round(blinking.length / C, 4),
  worst: history,
};
const name = `flicker_${tagOf(file)}_N${N}`;
writeFileSync(`${outDir}/${name}.json`, JSON.stringify(summary));

/* The plot's data: the blink frequency map and the strip over the densest blinking box. */
const box = (() => {
  let best = null;
  for (let s = -60; s <= 50; s += 2.5) for (let w = -180; w < 180; w += 2.5) {
    const e = w + 15 > 180 ? w + 15 - 360 : w + 15, b = [s, s + 10, w, e];
    let count = 0;
    for (const x of steps.onsets) if (inBox(x.i, b)) count++;
    if (!best || count > best.count) best = { box: b, count };
  }
  return best;
})();
const boxCells = Array.from({ length: C }, (_, i) => i).filter((i) => inBox(i, box.box));
let start = 0, most = -1;
for (let t = 0; t + 8 <= T; t++) {
  let count = 0;
  for (const x of steps.onsets) if (x.t >= t && x.t < t + 7 && inBox(x.i, box.box)) count++;
  if (count > most) { most = count; start = t; }
}
const plot = {
  ...header, steps: STEPS, dt, stride: STRIDE, jump: JUMP, transitions: steps.counted,
  lon: Array.from(lon, (v) => round(v, 3)), lat: Array.from(lat, (v) => round(v, 3)), land: Array.from(land), blinks: Array.from(blinkCount),
  strip: { box: box.box, start, cells: boxCells, lon: boxCells.map((i) => round(lon[i], 3)), lat: boxCells.map((i) => round(lat[i], 3)), land: boxCells.map((i) => land[i]),
    frames: Array.from({ length: 8 }, (_, n) => boxCells.map((i) => round(opacity[start + n][i], 3))),
    pageFrames: Array.from({ length: 8 }, (_, n) => (STRIDE * n + STRIDE - 1 < T ? boxCells.map((i) => round(opacity[STRIDE * n + STRIDE - 1][i], 3)) : null)).filter(Boolean),
    hours: Array.from({ length: 8 }, (_, n) => round((start + n + 1) * dt / 3600, 2)) },
  worst: worst.map((i) => [round(lon[i], 2), round(lat[i], 2)]),
};
writeFileSync(`${outDir}/${name}_plot.json`, JSON.stringify(plot));
const python = process.env.PYTHON ?? 'python3';
const drawn = spawnSync(python, [new URL('./figures/cloudFlicker.py', import.meta.url).pathname, `${outDir}/${name}_plot.json`, `${outDir}/${name}_map.png`, `${outDir}/${name}_strip.png`], { stdio: 'inherit' });
if (drawn.status !== 0) console.error('cloudFlicker.py failed');
console.log(JSON.stringify({ ...summary, worst: history.map((h) => ({ cell: h.cell, lat: h.lat, lon: h.lon, blinks: h.blinks, deckDriven: h.deckDrivenBlinks, regions: h.regions })) }, null, 1));
console.log(`wrote ${outDir}/${name}.json, ${name}_plot.json, ${name}_map.png, ${name}_strip.png`);
