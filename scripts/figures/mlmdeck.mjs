// The mixed-layer stratocumulus deck from a saved state, for mlmdeck.py:
// the state loaded on the GPU as a continuing segment and four steps
// taken; the physics kernel's inputs of the fourth are read before it
// runs and its deck after, and every ice-free sea column is re-diagnosed
// on the host by mixedLayer.module.js along the GPU's own path through
// mlmColumn in js/gpu/physics.gpu.js (the regime rest, the inversion
// ceiling, the ring-smoothed subsidence and its running mean, the gate
// memory, the carried height). The host gives the mixed layer no
// sunlight, so the liquid water path, cover and the heights after the
// step are compared on the night side alone, where the GPU's absorbed
// sunlight is zero too; the gates, the running means and the entrainment
// (the radiative closure's, driven by the longwave alone) everywhere.
// The cloud base the maps' thickness takes is the host's, by day without
// the step's solar heating. OCEAN, RADIATION and the rest as for
// spinup.mjs; the run uses OCEAN='{"everySteps":8}'.
//   node scripts/figures/mlmdeck.mjs <state.bin> <out.json>
import { writeFileSync } from 'node:fs';
import { readRanges } from '../../js/gpu/device.module.js';
import { createMixedLayer, dycomsLongwave, MIXED_LAYER_DEFAULTS } from '../../js/physics/mixedLayer.module.js';
import { DECK_CLOUD_LEVELS, UNDECIDED, sunDirection, ringMean } from '../../js/physics/radiation.module.js';
import { CP_DRY as cp, GRAVITY as g, P0 } from '../../js/dynamics/sigmaCore.module.js';
import { saturationHumidity, EPSILON } from '../../js/physics/moist.module.js';
import { readState, figureHeader, gpuModelFrom, dt as stepOf, DEG, rounded } from './figureState.mjs';

const [file, outFile] = process.argv.slice(2);
if (!file || !outFile) { console.error('usage: node scripts/figures/mlmdeck.mjs <state.bin> <out.json>'); process.exit(2); }
const saved = await readState(file);
const model = await gpuModelFrom(saved);
const { mesh, gpu } = model, C = mesh.nCells, dt = stepOf(saved.N), STEPS = 4;
const L = gpu.layout, K = gpu.K, KC = K * C, ph = gpu.physics, R = ph.R;
const read = async (buffer, names) => Object.fromEntries((await readRanges(gpu.device, gpu.buffers[buffer], names.map(([n, len]) => ({ offset: L[buffer][n], length: len })))).map((v, j) => [names[j][0], Float64Array.from(v)]));
let calls = 0, input = null;
const ocean = gpu.hooks.beforePhysics;
gpu.hooks.beforePhysics = async (step, n) => {
  if (ocean) await ocean(step, n);
  if (++calls !== STEPS) return;
  const [S, D, P, LV] = await Promise.all([
    read('S', [['PI', C], ['TH', KC], ['TS', C], ['Q', KC], ['QC', KC], ['ICE', C]]),
    read('D', [['GEO', KC], ['EXM', KC], ['EXL', KC], ['THV', KC], ['PSD', (K + 1) * C]]),
    read('PH', [['DEPTH', C], ['MLMSUB', C], ['MLMH', C], ['MLMGATE', C], ['CONC', C], ['SNOW', C], ['LAND', C], ['REGIME', C]]),
    read('LV', [['DS', K], ['SM', K], ['GABS', K]]),
  ]);
  input = { time: model.time, S, D, P, LV };
};
for (let n = 0; n < STEPS; n++) await model.step(dt);
await model.settle();
const out = await read('PH', [['MLMCOVER', C], ['MLMWATER', C], ['MLMENT', C], ['MLMSUB', C], ['DECK', C], ['DECKF', C], ['SH', C], ['EVAP', C], ['MLMH', C], ['MLMGATE', C], ['MLMTOP', C]]);
const { S, D, P, LV } = input;

const shadow = createMixedLayer({ cp, R, g, latentHeat: ph.latentHeat, referencePressure: P0, cloudLevels: DECK_CLOUD_LEVELS, ...ph.mixedLayer });
const longwave = dycomsLongwave(), Lh = ph.latentHeat, HMAX = shadow.maximumHeight, HMEM = { ...MIXED_LAYER_DEFAULTS, ...ph.mixedLayer }.heightMemory;
const prognostic = !!ph.prognosticHeight, restInversion = ph.deckRest !== 'depth', restRegime = ph.deckRest === 'regime' && ph.turbulence === 'moist';
const ceilingJump = ph.ceilingInversion ?? ph.minimumInversion, blGate = ph.deckRegime === 'boundaryLayer';
const sun = sunDirection(input.time);
const fresh = (x) => -Math.expm1(-x);
const zOf = (k, i) => (D.GEO[k * C + i] + LV.GABS[k]) / g;
const interfaceHeight = (i, m) => (D.GEO[m * C + i] + LV.GABS[m] + cp * D.THV[m * C + i] * (D.EXM[m * C + i] - D.EXL[(m - 1) * C + i])) / g;
function cappingOf(i, floor) {
  for (let k = K - 2; k >= 1; k--) {
    if (zOf(k + 1, i) >= HMAX) break;
    if (zOf(k, i) > floor && D.THV[k * C + i] - D.THV[(k + 1) * C + i] >= ceilingJump) return k;
  }
  return -1;
}
const GATES = ['active (G > 0.5 and cloudy)', 'land or full sea ice', 'boundary-layer top within the lowest layer', 'no layer above h', 'off: running-mean ascent faster than 1 mm/s', 'off: virtual jump below minimumInversion', 'runs: no cloud (base above h)', 'non-finite', 'off: G ≤ 0.5 although both tests pass now', 'off: lifted cumulus layer (regime rest)'];

function column(i) {
  const pi = S.PI[i], zBottom = zOf(K - 1, i), mixedDepth = P.DEPTH[i] - zBottom;
  if (!(mixedDepth > 0)) return { gate: 2 };
  const depth = mixedDepth + zBottom, regime = P.REGIME[i];
  const lifted = restRegime && (regime === 1 || regime === 2);
  const capping = prognostic || lifted ? cappingOf(i, depth) : -1;
  const ceiling = prognostic && capping >= 0 ? Math.min(HMAX, zOf(capping, i) - 1) : HMAX;
  const standDown = lifted && !(capping >= 0 && interfaceHeight(i, capping + 1) <= ph.cumulusCeiling);
  const resting = restInversion && !(restRegime && regime === 3) && !standDown && ceiling < HMAX ? ceiling : depth;
  const h = prognostic && P.MLMH[i] > 0 ? Math.max(depth, Math.min(ceiling, P.MLMH[i])) : prognostic ? resting : depth;
  const rested = P.MLMH[i] > 0 ? (HMEM > 0 ? P.MLMH[i] + (resting - P.MLMH[i]) * fresh(dt / HMEM) : resting) : P.MLMH[i];
  let weight = 0, heat = 0, water = 0, k = K - 1;
  for (; k >= 0 && zOf(k, i) < h; k--) {
    const idx = k * C + i, cloud = Math.max(0, S.QC[idx]);
    heat += LV.DS[k] * (S.TH[idx] - Lh * cloud / (cp * D.EXM[idx]));
    water += LV.DS[k] * (Math.max(0, S.Q[idx]) + cloud);
    weight += LV.DS[k];
  }
  if (k < 1) return { gate: 3, height: rested };
  let lowerHeight = 0, lower = K, m = K - 1;
  for (; m > k; m--) { const z = interfaceHeight(i, m); if (!(z < h)) break; lowerHeight = z; lower = m; }
  const flowAt = (n) => ringMean(mesh, D.PSD, n * C, i, ph.subsidenceSmoothing);
  const lowerFlow = lower < K ? flowAt(lower) : 0;
  const flow = lowerFlow + (flowAt(m) - lowerFlow) * (h - lowerHeight) / (interfaceHeight(i, m) - lowerHeight);
  const subsidence = -flow / (pi * LV.SM[m] / (R * D.THV[m * C + i] * D.EXM[m * C + i]) * g);
  const mean = P.MLMSUB[i] + (subsidence - P.MLMSUB[i]) * fresh(dt / ph.subsidenceMemory), sinking = !(mean > -ph.stratusSubsidence);
  const above = k * C + i, aboveCloud = Math.max(0, S.QC[above]);
  const forcing = {
    surfacePressure: pi, sensibleHeat: out.SH[i], evaporation: out.EVAP[i], radiation: longwave, subsidence: () => subsidence, absorbedSolar: null,
    thetaLAbove: S.TH[above] - Lh * aboveCloud / (cp * D.EXM[above]), qtAbove: Math.max(0, S.Q[above]) + aboveCloud,
  };
  const start = { h, thetaL: heat / weight, qt: water / weight };
  let now = null, pass = 0;
  if (sinking && !standDown) {
    now = shadow.diagnose(start, forcing);
    pass = (blGate ? regime === 3 : now.virtualJump >= ph.minimumInversion) ? 1 : 0;
  }
  const G = standDown ? 0 : ph.gateMemory > 0 ? P.MLMGATE[i] + (pass - P.MLMGATE[i]) * fresh(dt / ph.gateMemory) : pass;
  const piS = Math.pow(pi / P0, R / cp);
  const diag = { thetaL: start.thetaL, qt: start.qt, h, jump: now ? now.virtualJump : NaN, thetaLAbove: forcing.thetaLAbove, qtAbove: forcing.qtAbove, w: subsidence, ceiling, pass, regime };
  if (now) diag.base = cp * start.thetaL * (1 + (1 / EPSILON - 1) * start.qt) * (piS - now.cloudBaseExner) / g;
  if (!(G > UNDECIDED || (G === UNDECIDED && pass === 1)) || ph.deckBypass) return { mean, G, gate: standDown ? 9 : !sinking ? 4 : pass ? 8 : 5, diag, height: rested };
  now ??= shadow.diagnose(start, forcing);
  const stepped = shadow.step(start, forcing, dt, now), top = prognostic ? Math.max(depth, Math.min(ceiling, stepped.h)) : stepped.h;
  const next = shadow.diagnose({ h: top, thetaL: stepped.thetaL, qt: stepped.qt }, forcing);
  if (!(Number.isFinite(next.liquidWaterPath) && Number.isFinite(next.cover) && Number.isFinite(now.entrainment))) return { mean, G, gate: 7, diag, height: rested };
  return { mean, G, ok: true, diag, gate: next.cover > 0 || next.liquidWaterPath > 0 ? 0 : 6, lwp: next.liquidWaterPath, cover: next.cover, entrainment: now.entrainment, base: next.cloudBase, height: top };
}

const meanArea = mesh.areaCell.reduce((a, b) => a + b, 0) / C;
const r = rounded;
const rows = [], diagnosed = [], rel = [], relBig = [], coverDiff = [], entRel = [], entBig = [], subDiff = [], gateDiff = [], topDiff = [], onlyGpu = [], onlyHost = [];
let lwpMax = 0, both = 0, decisions = 0, nightSea = 0, gpuNight = 0;
for (let i = 0; i < C; i++) {
  const land = P.LAND[i] > 0.5, ice = S.ICE[i] > 0 ? (P.CONC[i] <= 0 ? 1 : P.CONC[i]) : 0;
  const mu = mesh.xCell[3 * i] * sun[0] + mesh.xCell[3 * i + 1] * sun[1] + mesh.xCell[3 * i + 2] * sun[2], night = !(mu > 0);
  const host = !land && 1 - ice > 0 ? column(i) : null;
  const gpuOn = out.MLMCOVER[i] > 0 || out.MLMWATER[i] > 0, hostOn = !!(host && host.ok && (host.cover > 0 || host.lwp > 0));
  if (host && host.mean !== undefined) subDiff.push(Math.abs(host.mean - out.MLMSUB[i]) / Math.max(1e-6, Math.abs(out.MLMSUB[i])));
  if (host && host.G !== undefined) { gateDiff.push(Math.abs(host.G - out.MLMGATE[i])); if ((host.G > UNDECIDED) !== (out.MLMGATE[i] > UNDECIDED)) decisions++; }
  if (host && host.ok && out.MLMENT[i] > 0) { entRel.push(Math.abs(host.entrainment - out.MLMENT[i]) / out.MLMENT[i]); if (out.MLMENT[i] >= 1e-4) entBig.push(entRel[entRel.length - 1]); }
  if (host && night) {
    nightSea++;
    if (gpuOn) gpuNight++;
    if (host.height !== undefined && (host.ok || out.MLMH[i] > 0)) topDiff.push(Math.abs(host.height - out.MLMH[i]));
    if (gpuOn && !hostOn) onlyGpu.push(i);
    if (hostOn && !gpuOn) onlyHost.push(i);
    if (gpuOn && hostOn) {
      both++;
      const d = Math.abs(host.lwp - out.MLMWATER[i]); lwpMax = Math.max(lwpMax, d);
      rel.push(d / Math.max(1e-9, out.MLMWATER[i])); coverDiff.push(Math.abs(host.cover - out.MLMCOVER[i]));
      if (out.MLMWATER[i] >= 1e-3) relBig.push(rel[rel.length - 1]);
    }
  }
  let lon = mesh.lonCell[i] * DEG; if (lon > 180) lon -= 360;
  const sea = !land, x = host && host.diag, bottom = (K - 1) * C + i, tLow = S.TH[bottom] * D.EXM[bottom];
  const rh = 100 * Math.max(0, S.Q[bottom]) / saturationHumidity(tLow, S.PI[i] * LV.SM[K - 1]);
  const base = host && host.ok ? host.base : null, top = gpuOn ? out.MLMH[i] : null;
  if (x && ice === 0 && [0, 4, 5, 6, 8, 9].includes(host.gate)) diagnosed.push({ lon, lat: mesh.latCell[i] * DEG, ...x, rh, sst: S.TS[i], sh: out.SH[i], evap: 86400 * out.EVAP[i] });
  rows.push([r(lon, 2), r(mesh.latCell[i] * DEG, 2), r(mesh.areaCell[i] / meanArea, 4), land ? 1 : 0, r(ice, 3), gpuOn ? 1 : 0, night ? 1 : 0,
    gpuOn && base !== null ? r(Math.max(0, top - base), 1) : null, gpuOn && base !== null ? r(base, 1) : null, gpuOn ? r(top, 1) : null, sea ? r(P.DEPTH[i], 1) : null,
    gpuOn ? r(1000 * out.MLMWATER[i], 2) : null, gpuOn ? r(out.MLMCOVER[i], 4) : null, gpuOn ? r(1000 * out.MLMENT[i], 3) : null, sea ? r(1000 * out.MLMSUB[i], 3) : null,
    sea ? r(out.DECKF[i], 4) : null, sea ? r(1000 * out.DECK[i], 2) : null, host && host.ok ? r(1000 * host.lwp, 2) : null, host ? host.gate : 1,
    x ? r(x.thetaL, 3) : null, x ? r(1000 * x.qt, 3) : null, x ? r(x.h, 1) : null, x ? r(x.base, 1) : null, x ? r(x.jump, 3) : null, x ? r(x.thetaLAbove, 3) : null, x ? r(1000 * x.qtAbove, 3) : null, x ? r(1000 * x.w, 4) : null,
    sea ? r(rh, 2) : null, sea ? r(tLow, 3) : null, sea ? r(S.TS[i], 3) : null, sea ? r(out.SH[i], 2) : null, sea ? r(86400 * out.EVAP[i], 3) : null,
    sea ? r(P.MLMH[i], 1) : null, sea ? r(out.MLMGATE[i], 4) : null, sea ? r(out.MLMH[i], 1) : null, x ? r(x.ceiling, 1) : null, x ? x.pass : null, sea ? P.REGIME[i] : null]);
}
const q = (a, p) => { if (!a.length) return NaN; const s = [...a].sort((x, y) => x - y); return s[Math.min(s.length - 1, Math.floor(p * s.length))]; };
const rms = (a) => Math.sqrt(a.reduce((s, x) => s + x * x, 0) / Math.max(1, a.length));
const e = (v) => v.toExponential(2);
const validation = [
  `day ${saved.day} N=${saved.N}: night side ${nightSea} ice-free sea columns, GPU deck on ${gpuNight}, host on ${onlyHost.length + both}, both ${both}; only GPU ${onlyGpu.length}, only host ${onlyHost.length}`,
  `night LWP host vs GPU over ${both}: max |diff| ${e(1000 * lwpMax)} g/m², relative max ${e(q(rel, 1))} median ${e(q(rel, 0.5))} rms ${e(rms(rel))} p99 ${e(q(rel, 0.99))}; over the ${relBig.length} with LWP ≥ 1 g/m² relative max ${e(q(relBig, 1))} median ${e(q(relBig, 0.5))} rms ${e(rms(relBig))}`,
  `night cover |diff| max ${e(q(coverDiff, 1))} median ${e(q(coverDiff, 0.5))}; carried height after the step |host − GPU| max ${e(q(topDiff, 1))} m median ${e(q(topDiff, 0.5))} m over ${topDiff.length}`,
  `all ice-free sea: entrainment relative rms ${e(rms(entRel))} max ${e(q(entRel, 1))} over ${entRel.length} (≥ 0.1 mm/s: rms ${e(rms(entBig))} max ${e(q(entBig, 1))} over ${entBig.length}); running-mean subsidence relative max ${e(q(subDiff, 1))} over ${subDiff.length}; gate memory |host − GPU| max ${e(q(gateDiff, 1))} over ${gateDiff.length}, columns whose pass differs (|ΔG| > dt/gateMemory/2) ${gateDiff.filter((x) => x > 0.5 * dt / ph.gateMemory).length}, decisions (G > ${UNDECIDED}) differing ${decisions}`,
];
for (const line of validation) console.log(line);
for (const [name, list] of [['only GPU', onlyGpu], ['only host', onlyHost]]) for (const i of list.slice(0, 5)) {
  console.log(`  ${name} cell ${i} (${(mesh.latCell[i] * DEG).toFixed(1)}, ${(mesh.lonCell[i] * DEG).toFixed(1)}): GPU cover ${out.MLMCOVER[i]} lwp ${out.MLMWATER[i]} sub ${out.MLMSUB[i]} gate ${out.MLMGATE[i]}; host ${JSON.stringify(column(i))}`);
}
const subsolar = [Math.atan2(sun[1], sun[0]) * DEG, Math.asin(sun[2]) * DEG];
writeFileSync(outFile, JSON.stringify({
  ...figureHeader(file, saved), time: input.time, dt, subsolar, validation,
  columns: ['lon', 'lat', 'area', 'land', 'ice', 'active', 'night', 'thickness', 'cloudBase', 'h', 'pblDepth', 'lwp', 'cover', 'entrainment', 'subsidence', 'deckf', 'deck', 'lwpHost', 'gate',
    'thetaL', 'qt', 'hNow', 'cloudBaseNow', 'virtualJump', 'thetaLAbove', 'qtAbove', 'subsidenceNow', 'rhLow', 'tLow', 'sst', 'sh', 'evap', 'mlmHeight', 'mlmGate', 'mlmHeightOut', 'ceiling', 'passNow', 'regime'],
  units: { thickness: 'm, GPU height after the step less the host cloud base', cloudBase: 'm, host', h: 'm, GPU height after the step', lwp: 'g/m², GPU', entrainment: 'mm/s, GPU', subsidence: 'mm/s, GPU running mean', deck: 'g/m², deck water the radiation uses', thetaL: 'K', qt: 'g/kg', hNow: 'm', cloudBaseNow: 'm, unclipped (may exceed h)', virtualJump: 'K', thetaLAbove: 'K', qtAbove: 'g/kg', subsidenceNow: 'mm/s at h, instantaneous, as used by step', rhLow: '%', tLow: 'K', sst: 'K', sh: 'W/m²', evap: 'mm/day (kg/m²/day)', mlmHeight: 'm, carried into the step (0 unset)', mlmGate: 'gate memory G after the step; the deck runs where G > 0.5', mlmHeightOut: 'm, carried height after the step', ceiling: 'm, inversion ceiling (capped at maximumHeight)', passNow: '1 where both gate tests pass at this step', regime: 'boundary-layer regime of the step before: 0 stable, 1 surface-driven, 2 decoupled, 3 coupled' },
  gates: GATES, rows,
}));
console.log(`wrote ${outFile} (${rows.length} cells; subsolar point ${subsolar[1].toFixed(1)}°, ${subsolar[0].toFixed(1)}°)`);
const TABLE = [['deck band', -180, -100, -15, 0], ['SE Pacific', -100, -80, -30, -10], ['Peru', -90, -75, -15, -5], ['Namibia', 5, 15, -30, -10], ['California', -135, -120, 20, 35], ['Canaries', -30, -15, 15, 35], ['warm pool', 140, 180, -15, 0]];
const QUANTITIES = [['h m', (c) => c.h, 0], ['cloud base m', (c) => c.base, 0], ['base − h m', (c) => c.base - c.h, 0], ['Δθv K', (c) => c.jump, 2], ['qt g/kg', (c) => 1000 * c.qt, 2], ['qt above g/kg', (c) => 1000 * c.qtAbove, 2],
  ['θl K', (c) => c.thetaL, 2], ['θl above K', (c) => c.thetaLAbove, 2], ['RH lowest %', (c) => c.rh, 1], ['SST K', (c) => c.sst, 2], ['SH W/m²', (c) => c.sh, 1], ['EVAP mm/d', (c) => c.evap, 2], ['w at h mm/s', (c) => 1000 * c.w, 3]];
const inBox = ([, x0, x1, y0, y1]) => (c) => c.lon >= x0 && c.lon <= x1 && c.lat >= y0 && c.lat <= y1;
const cells = TABLE.map((b) => diagnosed.filter(inBox(b)));
const finite = (list) => list.filter(Number.isFinite);
console.log(`\nmedians over ice-free sea cells with gate 0/4/5/6/8/9 (Δθv and cloud base where the gates were tested)\n${''.padEnd(16)}${TABLE.map((b) => b[0].padStart(12)).join('')}`);
for (const [name, f, n] of QUANTITIES) console.log(name.padEnd(16) + cells.map((list) => q(finite(list.map(f)), 0.5).toFixed(n).padStart(12)).join(''));
console.log('base > h'.padEnd(16) + cells.map((list) => (list.filter((c) => c.base > c.h).length / Math.max(1, list.length)).toFixed(2).padStart(12)).join(''));
console.log('cells'.padEnd(16) + cells.map((list) => String(list.length).padStart(12)).join(''));
process.exit(0);
