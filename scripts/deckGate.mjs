// The mixed-layer deck's gate over the stratocumulus boxes, column by column,
// on the single-thread CPU engine with the ocean off:
//   STEPS=4 RADIATION=… MOIST=… BOUNDARY_LAYER=… SURFACE=… node scripts/deckGate.mjs <state.bin>
// At the start of each step's physics the deck's decision is replicated from
// the inputs the radiation reads (shadowDeck of js/physics/radiation.module.js:
// the regime of the step before, the inversion ceiling, the carried height and
// its slab, the running-mean subsidence, the virtual jump at h, the gate's
// memory), checked against the model's own decision after the step, and each
// column-step is put in one class: runs; stood down because its regime is
// surface-driven or decoupled and it has no ceiling (no interface below the
// deck's maximum height across which θ_v rises by the ceiling's jump), or a
// ceiling above `cumulusCeiling`; off because the running-mean subsidence
// rises; off because the jump at h is below `minimumInversion`; off because
// the gate's memory has not turned. Beside the classes, each test's own
// failing share and the regimes' shares.
// The column at the same instant: the resolved inversion (the interface of
// largest θ_v jump between 300 m and 3 km above the surface, its height, its
// jump and the jump over the interfaces within 500 m above it), the largest
// one-interface jump under 2 km, the lowest run of cloudy layers (q_c above
// 10⁻⁶) below 3 km, its base, top and liquid water path, the boundary layer's
// depth, mixing top, decoupling height, entrainment velocity and cloud-top
// cooling, the shallow plume's share of columns, base flux and top, the sink
// at the inversion (π σ̇ ring-meaned as the deck reads it) and at 700 hPa,
// the saved running-mean sink, the slab's cloud base and the deck's own
// diagnosis at its h with the column's fluxes of the step before and no
// sunlight (cover, liquid water path, buoyancy integral ratio, entrainment),
// whether or not it runs; after the step the radiation's EIS and its share
// above MOIST's `stratocumulusInversion` (7 K without one), low cover and
// water, the deck's fraction and water, and the surface fluxes. Per box the
// area-weighted means (medians in brackets) and a profile of the box means
// by layer below 4 km.
import { readFileSync } from 'node:fs';
import { topographyFromInt16 } from '../js/geography.module.js';
import { createModel, STATE_NAMES } from '../js/model.module.js';
import { decodeState, savedLevels } from '../js/stateFile.module.js';
import { savedDeckField, DECK_FIELDS } from '../js/physics/regrid.module.js';
import { createMixedLayer, dycomsLongwave } from '../js/physics/mixedLayer.module.js';
import { DECK_CLOUD_LEVELS, ringMean } from '../js/physics/radiation.module.js';
import { LATENT_HEAT, saturationHumidity } from '../js/physics/moist.module.js';
import { inLongitudes } from '../js/audit.module.js';
import { REGIME } from '../js/physics/boundaryLayer.module.js';
import { stateGrid } from '../js/stateGrid.module.js';

const FILE = process.argv[2];
if (!FILE) { console.error('usage: node scripts/deckGate.mjs <state.bin>'); process.exit(1); }
const env = (name) => JSON.parse(process.env[name] ?? '{}');
const RADIATION = env('RADIATION'), MOIST = env('MOIST'), EIS_SWITCH = MOIST.stratocumulusInversion ?? 7, BOUNDARY_LAYER = env('BOUNDARY_LAYER'), SURFACE = env('SURFACE');
const STEPS = Number(process.env.STEPS ?? 4);
const saved = await decodeState(new Uint8Array(readFileSync(FILE)));
const topography = topographyFromInt16(readFileSync(new URL('../data/topography_0p25.bin', import.meta.url)).buffer);
const model = createModel(stateGrid(saved, topography).mesh, { topography, levels: savedLevels(saved), ocean: false, radiation: RADIATION, moist: MOIST, boundaryLayer: BOUNDARY_LAYER, surface: SURFACE });
const { mesh, core, state, phases, radiation, boundaryLayer: bl, moist, seaIce, land } = model;
const { K, C, sigmaMid, R, g, cp, p0, exnerLayer, exnerLower, geopotential, piSigmaDot } = core.diagnostics;
const { thetaV } = core.arrays;
STATE_NAMES.forEach((name, a) => state[a].set(saved[name]));
seaIce.load(state[6], saved.concentration ?? null);
for (const field of Object.keys(DECK_FIELDS)) radiation[field].set(savedDeckField(saved, field, model));
land.load(saved.land, state[6]);
for (const field of ['cumulusCover', 'cumulusWater', 'subcloudVirtual']) if (saved[field] && moist[field] && saved[field].length === moist[field].length) moist[field].set(saved[field]);
model.time = saved.time;
const zs = Float64Array.from({ length: C }, (_, i) => (model.surfaceGeopotential ? model.surfaceGeopotential[i] / g : 0));
if (saved.boundaryDepth) bl.depth.set(Float64Array.from(saved.boundaryDepth, (z, i) => z + zs[i]));
if (saved.mixingTop) bl.mixingTop.set(Float64Array.from(saved.mixingTop, (z, i) => z + zs[i]));
if (saved.boundaryRegime) bl.regime.set(saved.boundaryRegime);
if (saved.boundaryBuoyancy) bl.buoyancyFlux.set(saved.boundaryBuoyancy);
if (saved.evaporation) radiation.evaporation.set(saved.evaporation);
const [pi, theta, , surfaceT, q, qc, ice] = state;
const gates = radiation.deckGates, L = RADIATION.latentHeat ?? LATENT_HEAT;
const prognosticHeight = RADIATION.prognosticHeight ?? true, fractionalSlab = (RADIATION.deckSlab ?? 'fraction') === 'fraction', interpolatedReference = (RADIATION.deckReference ?? 'layer') === 'interpolate';
const shadow = createMixedLayer({ cp, R, g, latentHeat: L, referencePressure: p0, cloudLevels: DECK_CLOUD_LEVELS, ...(RADIATION.mixedLayer ?? {}) });
const longwave = dycomsLongwave();
const deg = 180 / Math.PI, area = mesh.areaCell, landMask = model.geography.land, dt = 1350 * 16 / saved.N, bottom = (K - 1) * C;
const lat = Float64Array.from(mesh.latCell, (x) => x * deg), lon = Float64Array.from(mesh.lonCell, (x) => x * deg);

const BOXES = [
  ['Peruvian 10-20S 90-80W (Klein and Hartmann 1993)', [-20, -10, -90, -80]],
  ['SE Pacific 10-30S 110-80W', [-30, -10, -110, -80]],
  ['Peru 5-20S 90-75W', [-20, -5, -90, -75]],
  ['20S 85-75W (VOCALS-REx)', [-22, -18, -85, -75]],
  ['Namibia 10-20S 0-10E', [-20, -10, 0, 10]],
  ['California 20-30N 130-120W', [20, 30, -130, -120]],
];
const masks = BOXES.map(([, [south, north, west, east]]) => Uint8Array.from({ length: C }, (_, i) => (!landMask[i] && !(ice[i] > 0) && lat[i] >= south && lat[i] <= north && inLongitudes(lon[i], west, east) ? 1 : 0)));
const traced = [];
for (let i = 0; i < C; i++) if (masks.some((m) => m[i])) traced.push(i);

const CLASSES = ['runs', 'stood down: no ceiling', 'stood down: ceiling above cumulusCeiling', 'off: subsidence', 'off: jump', 'off: gate memory', 'no slab'];
const FIELDS = ['inversion', 'inversionJump', 'inversionJump500', 'inversionGradient', 'jump2km', 'jump2kmHeight', 'ceilingHeight', 'ceilingJump', 'cloudBase', 'cloudTop', 'cloudPath', 'depth', 'mixingTop', 'decoupling', 'entrainment', 'cooling', 'plume', 'plumeFlux', 'plumeTop',
  'sinkInversion', 'omega700', 'sinkMean', 'h', 'slabBase', 'jumpAtH', 'deckCover', 'deckPath', 'deckRatio', 'deckEntrainment', 'eis', 'eisAbove', 'lowCover', 'lowWater', 'deckFraction', 'deckWater', 'sensible', 'latent', 'sst', 'airGap', 'surfaceRH', 'lts'];
const take = () => ({ weight: 0, classes: new Float64Array(CLASSES.length), regimes: new Float64Array(4), fails: { stood: 0, sink: 0, jump: 0 }, sums: Object.fromEntries(FIELDS.map((f) => [f, { s: 0, w: 0, xs: [] }])), runs: Object.fromEntries(['h', 'deckCover', 'deckPath', 'slabBase', 'deckRatio', 'deckEntrainment', 'jumpAtH'].map((f) => [f, { s: 0, w: 0, xs: [] }])), profile: Array.from({ length: K }, () => ({ w: 0, z: 0, theta: 0, thetaV: 0, q: 0, rh: 0, qc: 0, cloudy: 0 })) });
const boxes = BOXES.map(take);
const replica = { checked: 0, decisions: 0, gate: 0 };
const record = (acc, field, w, x) => { if (!Number.isFinite(x)) return; const a = acc[field]; a.s += w * x; a.w += w; a.xs.push(x); };

function interfaceHeightOf(i, surface) {
  return (m) => (geopotential[m * C + i] + cp * thetaV[m * C + i] * (exnerLayer[m * C + i] - exnerLower[(m - 1) * C + i]) - surface) / g;
}

function column(i) {
  const b = bottom + i, surface = geopotential[b] - cp * thetaV[b] * (exnerLower[b] - exnerLayer[b]), interfaceHeight = interfaceHeightOf(i, surface);
  const mixedDepth = bl.depth[i] - geopotential[b] / g, depth = mixedDepth + (geopotential[b] - surface) / g;
  const out = { regime: bl.regime[i], class: 6, pass: 0, gate: NaN, sinking: true, standDown: false, jumpAtH: NaN, values: {} };
  const v = out.values, height = (k) => (geopotential[k * C + i] - surface) / g;
  let best = -1, bestJump = -Infinity, jump2 = -Infinity, jump2Height = NaN;
  for (let k = K - 2; k >= 1; k--) {
    const zi = interfaceHeight(k + 1);
    if (zi > 3000) break;
    const jump = thetaV[k * C + i] - thetaV[(k + 1) * C + i];
    if (zi >= 300 && jump > bestJump) { bestJump = jump; best = k; }
    if (zi <= 2000 && zi > depth && jump > jump2) { jump2 = jump; jump2Height = zi; }
  }
  if (best >= 0) {
    const zi = interfaceHeight(best + 1);
    let top = best;
    while (top > 1 && interfaceHeight(top) <= zi + 500) top--;
    v.inversion = zi; v.inversionJump = bestJump; v.inversionJump500 = thetaV[top * C + i] - thetaV[(best + 1) * C + i];
    v.inversionGradient = 1000 * bestJump / (height(best) - height(best + 1));
  }
  v.jump2km = jump2; v.jump2kmHeight = jump2Height;
  let base = NaN, cloudTop = NaN, path = 0, inRun = false, ended = false;
  for (let k = K - 1; k >= 1 && height(k) < 3000; k--) {
    const idx = k * C + i, water = Math.max(0, qc[idx]);
    if (water > 1e-6 && !ended) { if (!inRun) { base = k < K - 1 ? interfaceHeight(k + 1) : 0; inRun = true; } cloudTop = interfaceHeight(k); }
    else if (inRun) ended = true;
    path += water * pi[i] * core.diagnostics.dSigma[k] / g;
  }
  v.cloudBase = base; v.cloudTop = cloudTop; v.cloudPath = 1e3 * path;
  v.depth = depth; v.mixingTop = bl.mixingTop[i] - geopotential[b] / g + (geopotential[b] - surface) / g; v.decoupling = bl.regime[i] === REGIME.DECOUPLED ? bl.decoupling[i] - zs[i] : NaN;
  v.entrainment = 1e3 * bl.entrainment[i]; v.cooling = bl.cloudTopCooling[i];
  v.sst = surfaceT[i] - 273.15; v.airGap = surfaceT[i] - theta[b] * exnerLayer[b]; v.surfaceRH = q[b] / saturationHumidity(theta[b] * exnerLayer[b], pi[i] * sigmaMid[K - 1]);
  let k700 = 0; for (let k = 1; k < K; k++) if (Math.abs(pi[i] * sigmaMid[k] - 7e4) < Math.abs(pi[i] * sigmaMid[k700] - 7e4)) k700 = k;
  v.lts = theta[k700 * C + i] - theta[b];
  v.omega700 = 0.5 * (piSigmaDot[k700 * C + i] + piSigmaDot[(k700 + 1) * C + i]);
  if (best >= 0) {
    const density = pi[i] * sigmaMid[best] / (R * thetaV[best * C + i] * exnerLayer[best * C + i]);
    v.sinkInversion = 1e3 * ringMean(mesh, piSigmaDot, (best + 1) * C, i, gates.subsidenceSmoothing) / (density * g);
  }
  v.sinkMean = -1e3 * radiation.mlmSubsidence[i];

  const regime = gates.deckRest === 'regime' ? bl.regime[i] : -1;
  const lifted = regime === REGIME.SURFACE || regime === REGIME.DECOUPLED;
  let capping = -1;
  if (prognosticHeight || lifted) for (let k = K - 2; k >= 1; k--) {
    if ((geopotential[(k + 1) * C + i] - surface) / g >= shadow.maximumHeight) break;
    if ((geopotential[k * C + i] - surface) / g > depth && thetaV[k * C + i] - thetaV[(k + 1) * C + i] >= gates.ceilingInversion) { capping = k; break; }
  }
  const ceiling = prognosticHeight && capping >= 0 ? (geopotential[capping * C + i] - surface) / g - 1 : Infinity;
  if (capping >= 0) { v.ceilingHeight = interfaceHeight(capping + 1); v.ceilingJump = thetaV[capping * C + i] - thetaV[(capping + 1) * C + i]; }
  const standDown = lifted && !(capping >= 0 && interfaceHeight(capping + 1) <= gates.cumulusCeiling);
  const resting = gates.deckRest !== 'depth' && regime !== REGIME.COUPLED && !standDown && ceiling < shadow.maximumHeight ? ceiling : depth;
  const h = prognosticHeight && radiation.mlmHeight[i] > 0 ? shadow.bound(radiation.mlmHeight[i], depth, ceiling) : prognosticHeight ? resting : depth;
  out.standDown = standDown; out.capping = capping;
  let weight = 0, heat = 0, water = 0, k = K - 1;
  for (; k >= 0 && geopotential[k * C + i] - surface < g * h; k--) {
    const idx = k * C + i, cloud = Math.max(0, qc[idx]);
    heat += core.diagnostics.dSigma[k] * (theta[idx] - L * cloud / (cp * exnerLayer[idx]));
    water += core.diagnostics.dSigma[k] * (Math.max(0, q[idx]) + cloud);
    weight += core.diagnostics.dSigma[k];
  }
  if (fractionalSlab && k >= 1 && k < K - 1) {
    const part = h < interfaceHeight(k + 1) ? k + 1 : k;
    if (part !== capping) {
      const pdx = part * C + i, lowerEdge = part < K - 1 ? interfaceHeight(part + 1) : 0, cloud = Math.max(0, qc[pdx]);
      const share = Math.min(1, Math.max(0, (h - lowerEdge) / (interfaceHeight(part) - lowerEdge))) - (part === k + 1 ? 1 : 0);
      heat += share * core.diagnostics.dSigma[part] * (theta[pdx] - L * cloud / (cp * exnerLayer[pdx]));
      water += share * core.diagnostics.dSigma[part] * (Math.max(0, q[pdx]) + cloud);
      weight += share * core.diagnostics.dSigma[part];
    }
  }
  if (k < 1) return out;
  const above = k * C + i, aboveCloud = Math.max(0, qc[above]);
  let lowerHeight = 0, lower = K, m = K - 1;
  for (; m > k && interfaceHeight(m) < h; m--) { lowerHeight = interfaceHeight(m); lower = m; }
  const upperHeight = interfaceHeight(m);
  const lowerFlow = lower < K ? ringMean(mesh, piSigmaDot, lower * C, i, gates.subsidenceSmoothing) : 0;
  const flow = lowerFlow + (ringMean(mesh, piSigmaDot, m * C, i, gates.subsidenceSmoothing) - lowerFlow) * (h - lowerHeight) / (upperHeight - lowerHeight);
  const density = pi[i] * sigmaMid[m] / (R * thetaV[m * C + i] * exnerLayer[m * C + i]);
  const subsidence = -flow / (density * g), keep = Math.exp(-dt / gates.subsidenceMemory);
  const mean = radiation.mlmSubsidence[i] * keep + subsidence * (1 - keep);
  const sinking = !(mean > -gates.stratusSubsidence);
  let thetaLAbove = theta[above] - L * aboveCloud / (cp * exnerLayer[above]), qtAbove = Math.max(0, q[above]) + aboveCloud;
  if (interpolatedReference && k < K - 1) {
    const higher = above - C, higherCloud = Math.max(0, qc[higher]), below = geopotential[above + C] - surface;
    const share = Math.min(1, Math.max(0, (g * h - below) / (geopotential[above] - surface - below)));
    thetaLAbove += share * (theta[higher] - L * higherCloud / (cp * exnerLayer[higher]) - thetaLAbove);
    qtAbove += share * (Math.max(0, q[higher]) + higherCloud - qtAbove);
  }
  const forcing = { surfacePressure: pi[i], sensibleHeat: radiation.sensibleHeat[i], evaporation: radiation.evaporation[i], radiation: longwave, subsidence: () => subsidence, thetaLAbove, qtAbove };
  const d = shadow.diagnose({ h, thetaL: heat / weight, qt: water / weight }, forcing);
  const regimeTest = gates.deckRegime === 'boundaryLayer' ? bl.regime[i] === REGIME.COUPLED : d.virtualJump >= gates.minimumInversion;
  const pass = sinking && !standDown && regimeTest ? 1 : 0;
  const gate = standDown ? 0 : gates.gateMemory > 0 ? radiation.mlmGate[i] - (pass - radiation.mlmGate[i]) * Math.expm1(-dt / gates.gateMemory) : pass;
  const runs = (gate > 0.5 || (gate === 0.5 && pass === 1)) && !gates.deckBypass;
  out.class = runs ? 0 : standDown ? (capping < 0 ? 1 : 2) : pass ? 5 : !sinking ? 3 : 4;
  out.pass = pass; out.gate = gate; out.sinking = sinking; out.jumpAtH = d.virtualJump; out.runs = runs; out.regimeTest = regimeTest;
  v.h = h; v.slabBase = d.cloudBase; v.jumpAtH = d.virtualJump;
  v.deckCover = d.cover; v.deckPath = 1e3 * d.liquidWaterPath; v.deckRatio = Math.min(2, d.buoyancyIntegralRatio); v.deckEntrainment = 1e3 * d.entrainment;
  const topP = moist.cumulusTop[i];
  v.plume = moist.cumulusBaseFlux[i] > 0 ? 1 : 0;
  if (v.plume) {
    v.plumeFlux = 1e3 * moist.cumulusBaseFlux[i];
    let kTop = 0, nearest = Infinity;
    for (let n = 0; n <= K; n++) { const dp = Math.abs(pi[i] * core.levels[n] - topP); if (dp < nearest) { nearest = dp; kTop = n; } }
    v.plumeTop = kTop >= 1 && kTop < K ? interfaceHeight(kTop) : NaN;
  }
  return out;
}

const decided = new Map();
const physicsPhase = phases.physics;
phases.physics = (...args) => {
  decided.clear();
  for (const i of traced) {
    const out = column(i);
    decided.set(i, out);
    const w = area[i];
    masks.forEach((mask, n) => {
      if (!mask[i]) return;
      const box = boxes[n];
      box.weight += w; box.classes[out.class] += w;
      if (out.regime >= 0 && out.regime <= 3) box.regimes[out.regime] += w;
      if (out.standDown) box.fails.stood += w;
      if (!out.sinking) box.fails.sink += w;
      if (out.regimeTest === false) box.fails.jump += w;
      for (const [field, x] of Object.entries(out.values)) if (box.sums[field]) record(box.sums, field, w, x);
      if (out.runs) for (const field of Object.keys(box.runs)) record(box.runs, field, w, out.values[field]);
      const surface = geopotential[bottom + i] - cp * thetaV[bottom + i] * (exnerLower[bottom + i] - exnerLayer[bottom + i]);
      for (let k = K - 1; k >= 0; k--) {
        const idx = k * C + i, z = (geopotential[idx] - surface) / g;
        if (z > 4000) break;
        const p = box.profile[k], T = theta[idx] * exnerLayer[idx];
        p.w += w; p.z += w * z; p.theta += w * theta[idx]; p.thetaV += w * thetaV[idx]; p.q += w * 1e3 * q[idx];
        p.rh += w * q[idx] / saturationHumidity(T, pi[i] * sigmaMid[k]); p.qc += w * 1e3 * Math.max(0, qc[idx]); p.cloudy += w * (qc[idx] > 1e-6 ? 1 : 0);
      }
    });
  }
  physicsPhase(...args);
  for (const i of traced) {
    const out = decided.get(i);
    if (out.class === 6) continue;
    replica.checked++;
    if (out.runs !== radiation.mlmTop[i] > 0) replica.decisions++;
    replica.gate = Math.max(replica.gate, Math.abs(out.gate - radiation.mlmGate[i]));
    const w = area[i];
    masks.forEach((mask, n) => {
      if (!mask[i]) return;
      const sums = boxes[n].sums;
      record(sums, 'eis', w, radiation.stabilityIndex[i]); record(sums, 'eisAbove', w, radiation.stabilityIndex[i] > EIS_SWITCH ? 1 : 0); record(sums, 'lowCover', w, radiation.lowCover[i]); record(sums, 'lowWater', w, 1e3 * radiation.lowWater[i]);
      record(sums, 'deckFraction', w, radiation.stratusFraction[i]); record(sums, 'deckWater', w, 1e3 * radiation.stratusFraction[i] * radiation.stratus[i]);
      record(sums, 'sensible', w, radiation.sensibleHeat[i]); record(sums, 'latent', w, L * radiation.evaporation[i]);
    });
  }
};

const t0 = performance.now();
for (let n = 0; n < STEPS; n++) model.step(dt);
const f = (x, n = 2) => (Number.isFinite(x) ? x.toFixed(n) : '-');
const median = (xs) => { if (!xs.length) return NaN; const s = Float64Array.from(xs).sort(); return s[s.length >> 1]; };
console.log(`${FILE.replace(/.*\//, '')}: day ${saved.day}, ${STEPS} CPU steps of ${dt} s in ${((performance.now() - t0) / 1000).toFixed(0)} s; replica: ${replica.decisions} of ${replica.checked} column-step decisions apart, gate within ${replica.gate.toExponential(1)}`);
console.log(`gates: ceiling jump ${gates.ceilingInversion} K, cumulusCeiling ${gates.cumulusCeiling} m, minimumInversion ${gates.minimumInversion} K, sink at least ${1e3 * gates.stratusSubsidence} mm/s, gate memory ${gates.gateMemory / 3600} h, maximum height ${shadow.maximumHeight} m`);
const LABELS = {
  inversion: 'resolved inversion, m', inversionJump: 'its θv jump, K', inversionJump500: 'θv jump over the interfaces within 500 m above it, K', inversionGradient: 'its dθv/dz, K/km', jump2km: 'largest one-interface θv jump above the depth below 2 km, K', jump2kmHeight: 'its height, m',
  ceilingHeight: "the deck's ceiling interface (4 K), m", ceilingJump: 'its jump, K', cloudBase: 'lowest cloudy run below 3 km: base, m', cloudTop: 'its top, m', cloudPath: 'resolved cloud water below 3 km, g/m²', depth: 'boundary-layer depth, m', mixingTop: 'mixing top, m', decoupling: 'decoupling height z_b (decoupled columns), m',
  entrainment: "the boundary layer's entrainment velocity, mm/s", cooling: 'cloud-top cooling ΔF, W/m²', plume: 'shallow plume on (share)', plumeFlux: 'its base flux, g/m²/s', plumeTop: 'its top, m', sinkInversion: 'sink at the resolved inversion (ring-meaned), mm/s', omega700: 'π σ̇ at 700 hPa (descent > 0), Pa/s', sinkMean: "the deck's running-mean sink, mm/s",
  h: "the deck's h at the start, m", slabBase: "the slab's cloud base, m", jumpAtH: 'virtual jump at h, K', deckCover: "the deck's diagnosed cover (night, own fluxes)", deckPath: 'its liquid water path, g/m²', deckRatio: 'its buoyancy integral ratio (at most 2)', deckEntrainment: 'its entrainment, mm/s',
  eis: 'EIS, K', eisAbove: `EIS above ${EIS_SWITCH} K (share)`, lowCover: "the radiation's low cover", lowWater: 'low water, g/m²', deckFraction: "the deck's fraction", deckWater: "the deck's water (fraction × path), g/m²", sensible: 'sensible heat, W/m²', latent: 'latent heat, W/m²', sst: 'SST, °C', airGap: 'SST less the lowest air, K', surfaceRH: 'lowest layer RH', lts: 'LTS (θ at σ 0.7 less the lowest layer), K',
};
BOXES.forEach(([name], n) => {
  const box = boxes[n], W = box.weight;
  if (!(W > 0)) return;
  console.log(`\n== ${name}: ${masks[n].reduce((a, b) => a + b, 0)} sea cells`);
  console.log(`  classes: ${CLASSES.map((c, j) => `${c} ${f(box.classes[j] / W, 3)}`).join('; ')}`);
  console.log(`  failing alone: stood down ${f(box.fails.stood / W, 3)}, subsidence ${f(box.fails.sink / W, 3)}, jump at h ${f(box.fails.jump / W, 3)}; regimes stable ${f(box.regimes[0] / W, 3)}, surface-driven ${f(box.regimes[1] / W, 3)}, decoupled ${f(box.regimes[2] / W, 3)}, coupled ${f(box.regimes[3] / W, 3)}`);
  for (const field of FIELDS) {
    const a = box.sums[field];
    if (!(a.w > 0)) continue;
    console.log(`  ${LABELS[field].padEnd(62)} ${f(a.s / a.w, field === 'omega700' ? 3 : 2).padStart(9)} [${f(median(a.xs), field === 'omega700' ? 3 : 2)}] (${f(a.w / W, 2)})`);
  }
  const ran = Object.entries(box.runs).filter(([, a]) => a.w > 0);
  if (ran.length) console.log(`  where it runs: ${ran.map(([field, a]) => `${field} ${f(a.s / a.w, 1)} [${f(median(a.xs), 1)}]`).join(', ')}`);
  console.log('     z m   θ K    θv K   q g/kg   RH   qc g/kg cloudy');
  for (let k = K - 1; k >= 0; k--) {
    const p = box.profile[k];
    if (!(p.w > 0)) continue;
    console.log(`  ${(p.z / p.w).toFixed(0).padStart(6)} ${(p.theta / p.w).toFixed(2).padStart(7)} ${(p.thetaV / p.w).toFixed(2).padStart(7)} ${(p.q / p.w).toFixed(2).padStart(7)} ${(p.rh / p.w).toFixed(3).padStart(6)} ${(p.qc / p.w).toFixed(4).padStart(8)} ${(p.cloudy / p.w).toFixed(2).padStart(6)}`);
  }
});
process.exit(0);
