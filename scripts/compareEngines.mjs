// The CPU engine against the GPU engine at a saved state, field by field.
//
//   node scripts/compareEngines.mjs STATE.bin [STEPS=1,2] [--out FILE] [--cpu-only]
//
// Builds the CPU (f64) and GPU (f32) models of one state, steps both and, after each count in STEPS,
// compares every state array, the ocean (h, u, T, S, eta and the heat and salt contents h·T, h·S), the
// land and every region of the GPU's PH buffer against the CPU array it stands for. The regions are
// matched by the fixed name map PH_MAP below, not by content, so a region that is all zero on one side
// and not on the other cannot pass as a best match of something else. A PH region that is in neither
// PH_MAP nor PH_INPUT is reported first as UNMAPPED, and the run exits 1 when any region is unmapped, is
// missing on the CPU, is all zero on one side only, is all zero on the GPU without being listed in
// MAY_BE_ZERO, holds a non-finite value, differs from its CPU array by more than TOL (env, relative rms,
// default 1e-3; a row in ABSOLUTE_FLOOR passes too while its largest gap is under the floor, and once the
// ocean has stepped the row in CARRIED (PEN) counts only each cell's gap past ρc_p times twice the largest
// f32 grain of the heat stores h·(T − HEAT_REFERENCE) above PENETRATION_DEPTH, under which the GPU's
// deposit of the light leaves a remainder that waits for the next ocean step), or, for the
// DECISION_ROWS below, differs by more than TOL times the row's largest CPU magnitude in a larger share of
// its entries than the row's bound, or has its gaps under that threshold exceed, as a root sum of squares,
// 3·TOL of the root sum of squares of the CPU's values over all its entries. The rows in ABSOLUTE_BOUNDS
// (T2M) are judged in their own units instead, by the rms gap and the share of cells past a gap, whatever
// TOL is. A row reads n/a, and cannot fail, where the model does not compute it: ROW_NEEDS gives T2M n/a
// without the roughness exchange, STRESS and SDRAG without the implicit drag, WARM without the warm
// layer, COOL without the cool skin and UO without the relative wind, and the clear-sky, form-drag and
// held-radiation rows read n/a without the clear-sky pass, the form drag or the held radiation. The run
// also exits 1 when the CPU model takes the raster sub-grid orography, which has no form drag and is not
// the configuration a saved state runs; env NO_SUBGRID builds both models without sub-grid orography and
// lifts this. The state arrays are only reported, with each engine's own rms change since the last
// report. The ocean's rows, OFLUX and UO change only at an ocean step, the first at step 8 at N=64 and
// N=128, so a run through steps 1 and 2 compares them at their loaded values, UO at zero from a save
// without a sea flow.
//   --cpu-only   builds only the CPU model (no GPU is imported), steps it, and checks that every
//                mapped CPU array exists and is not all zero: the dry run of this tool
//   --gpu-dump P builds only the GPU model, steps it and writes its states at 0 and at each count in STEPS to
//                P.step<n> (each array in its own type): the GPU half, a minute on the GPU at N=128
//   --against P  builds only the CPU model and compares it at each count in STEPS with P.step<n> from a
//                --gpu-dump run of the same state, options and STEPS: the CPU half, which needs no GPU
//   --out FILE   also writes the report to FILE;  env FLAGS sets Dawn toggles for the GPU adapter,
//                env OCEAN a JSON object of ocean options for both models, env BOUNDARY_LAYER one of
//                boundary-layer options, env RADIATION one of
//                radiation options, env MOIST one of moist options, env SURFACE one of surface options and
//                env LAND one of land options
//   The default runs the clear-sky pass, which the page does not: RADIATION='{"clearSkyPass":false}' compares the page's shader.
// Run it with the GPU free (CLAUDE.md): one GPU model at a time.
import { readFileSync, writeFileSync } from 'node:fs';
import { serialize, deserialize } from 'node:v8';
import { createRequire } from 'node:module';
import { fileURLToPath } from 'node:url';
import { stateTerrain } from '../js/stateFile.module.js';
import { refuseUnreadOptions } from './optionEnvironment.mjs';

refuseUnreadOptions(import.meta.url);
const REPO = fileURLToPath(new URL('..', import.meta.url)).replace(/\/$/, '');
const flagValue = (flag) => (process.argv.includes(flag) ? process.argv[process.argv.indexOf(flag) + 1] : null);
const args = process.argv.slice(2).filter((a) => !a.startsWith('--') && !['--out', '--gpu-dump', '--against'].some((flag) => a === flagValue(flag)));
const FILE = args[0], STEPS = (args[1] ?? '1,2').split(',').map(Number), CPU_ONLY = process.argv.includes('--cpu-only');
const GPU_DUMP = flagValue('--gpu-dump'), AGAINST = flagValue('--against'), NO_GPU = CPU_ONLY || !!AGAINST;
const OUT = process.argv.includes('--out') ? process.argv[process.argv.indexOf('--out') + 1] : null, TOL = Number(process.env.TOL ?? 1e-3);
if (!FILE) { console.error('usage: node scripts/compareEngines.mjs STATE.bin [STEPS=1,2] [--out FILE] [--cpu-only]'); process.exit(2); }

// PH region -> the CPU model's array ([tail]: the CPU array is K·C long and the region holds its last layers)
export const PH_MAP = {
  SFLUX: 'radiation.surfaceFlux', OFLUX: 'seaIce.oceanFlux', CAP: 'ocean.capacity', MIX: 'boundaryLayer.mixing', MIXH: 'boundaryLayer.scalarMixing', FSTRESS: 'boundaryLayer.formStress', TOFD: 'boundaryLayer.formRate[tail]', OFLT: 'boundaryLayer.formDrag.sigma', DEPTH: 'boundaryLayer.depth', MIXTOP: 'boundaryLayer.mixingTop',
  RAIN: 'moist.precipitation', STEPRAIN: 'moist.rain', CONV: 'moist.convectivePrecipitation', COND: 'moist.largeScalePrecipitation', FOG: 'moist.fogPrecipitation',
  SH: 'radiation.sensibleHeat', EVAP: 'radiation.evaporation', SWDN: 'radiation.surfaceShortwave', BUOY: 'boundaryLayer.buoyancyFlux', USTAR: 'boundaryLayer.friction',
  STRESS: 'boundaryLayer.surfaceStress', SDRAG: 'boundaryLayer.surfaceDrag', HEATX: 'exchange.heat', REFX: 'exchange.reference', DRAG: 'exchange.drag', XWIND: 'exchange.wind', T2M: 'exchange.temperature2m',
  UO: 'seaFlow', WARM: 'skin.warm', COOL: 'skin.cool',
  REGIME: 'boundaryLayer.regime', ENTRAIN: 'boundaryLayer.entrainment', VRAD: 'boundaryLayer.radiativeVelocity', CTCOOL: 'boundaryLayer.cloudTopCooling',
  DECK: 'radiation.stratus', DECKF: 'radiation.stratusFraction', STRAT: 'radiation.stratiform', MLMCOVER: 'radiation.mlmCover', MLMWATER: 'radiation.mlmWater', MLMENT: 'radiation.mlmEntrainment',
  MLMTOP: 'radiation.mlmTop', MLMH: 'radiation.mlmHeight', MLMSUB: 'radiation.mlmSubsidence', MLMGATE: 'radiation.mlmGate', LWH: 'radiation.longwave',
  REFLSUM: 'radiation.summed.reflectedSolar', ABSSUM: 'radiation.summed.absorbedSolar',
  OLR: 'radiation.outgoing', OLRSUM: 'radiation.summed.outgoingLongwave', ATMSUM: 'radiation.summed.atmosphereSolar',
  INSSUM: 'radiation.summed.insolation', ABSCLRSUM: 'radiation.summed.clearAbsorbedSolar', OLRCLRSUM: 'radiation.summed.clearOutgoingLongwave',
  ASRMEAN: 'radiation.meanAbsorbedSolar', OLRMEAN: 'radiation.meanOutgoingLongwave', ALBMEAN: 'radiation.meanPlanetaryAlbedo', SWCREMEAN: 'radiation.meanShortwaveCloudEffect', LWCREMEAN: 'radiation.meanLongwaveCloudEffect',
  CUMF: 'moist.cumulusBaseFlux', CUTOP: 'moist.cumulusTop', CUCOVER: 'moist.cumulusCover[tail]', CUWATER: 'moist.cumulusWater[tail]', SUBTV: 'moist.subcloudVirtual', CONVMEAN: 'moist.convectiveRain', CONDMEAN: 'moist.largeScaleRain',
  PEN: 'ocean.penetrated', TRANS: 'ocean.transmission', CONC: 'seaIce.concentration', SNOW: 'seaIce.snow', SNOWALB: 'seaIce.snowAlbedo',
  SURF: 'land.surface', SOIL: 'land.soil', RUNOFF: 'land.runoff', DISCHARGE: 'land.discharge', ICELOSS: 'land.iceLoss', MLMK: 'radiation.mlmLayer', MLMM: 'radiation.mlmShare', DECKK: 'radiation.deckLayer', DECKFROM: 'radiation.deckFrom', VEG: 'land.vegetation', CANOPY: 'land.canopy', SEASONL: 'land.seasonLength', SEASONW: 'land.seasonWarmth', DEMAND: 'land.demandMean', RAINMEAN: 'land.rainMean',
  LITTERM: 'land.litterMean', DECAYM: 'land.decayMean', SOILC: 'land.soilCarbon', SNOWFREEV: 'land.snowFreeCover', SOILT2: 'land.soilT2', SOILT3: 'land.soilT3', SOILT4: 'land.soilT4', DEEP: 'land.deepWater',
  GWE: 'gravityWaves.east', GWN: 'gravityWaves.north', OLAUNCH: 'orography.launch', OWAVE: 'orography.wave', OSTRESS: 'orography.stress', OBETA: 'orography.beta', OBLOCK: 'orography.blocking', ODIR: 'orography.direction',
  OSTD: 'orography.fields.deviation', OORI: 'orography.fields.orientation', OANI: 'orography.fields.anisotropy', OSLP: 'orography.fields.slope',
  RADT: 'radiation.held.escape', RADTC: 'radiation.held.clearEscape', RADDF: 'radiation.held.share', RADDIR: 'radiation.held.direct', RADREFL: 'radiation.held.reflected', RADABS: 'radiation.held.absorbed',
  RADSWDN: 'radiation.held.down', RADSW: 'radiation.held.shortwave', RADTOA: 'radiation.held.top', RADBACK: 'radiation.held.back', RADOLR: 'radiation.held.outgoing', RADATM: 'radiation.held.atmosphere',
  RADCLROLR: 'radiation.held.clearOutgoing', RADCLRSW: 'radiation.held.clearAbsorbed', RADEMIT: 'radiation.held.emission', RADMU: 'radiation.held.mu',
  AER: 'radiation.aerosolFields',
};
// regions that are inputs or hand-offs with no CPU array of their own (they are checked for non-finite values only)
export const PH_INPUT = {
  TAU: 'relaxation time, an input', LAKE: 'lake state slots, some compared by test/lakeGpu.test.mjs', ADIF: 'diffuse albedo, a local of model.module.js', LAND: 'surface flags 0/1/2, an input', CUFALL: 'plume-to-rain hand-off', CURESERVE: 'plume-to-rain hand-off', CUFROZEN: 'plume-to-rain hand-off',
  CUDEEP: 'plume-to-rain hand-off', CUBASE: 'plume-to-rain hand-off', CURAIN: 'cumulus part of CONV', LWSFCSUM: 'running sum for the means', CLOUDK: 'cloud layer index, written inside the physics kernel',
  CUSHALLOW: 'plume diagnostic', CUSNOW: 'plume diagnostic', MOMU: 'momentum-layer scratch', MOMK: 'momentum-layer scratch', MOMD: 'momentum-layer scratch', MOMKD: 'momentum-layer scratch', MOMS: 'momentum-layer scratch',
  ABS: 'per-step value whose CPU counterpart is the running sum ABSSUM compares', INS: 'per-step value whose CPU counterpart is the running sum INSSUM compares', REFL: 'per-step value whose CPU counterpart is the running sum REFLSUM compares', ATMSW: 'per-step value whose CPU counterpart is the running sum ATMSUM compares',
  FOGMEAN: 'running mean', RADP: 'held radiation scratch', LYR: 'layered shortwave hand-off from the radiation kernel to its pass',
};
// mapped regions that may legitimately be all zero in a configuration (env MAY_BE_ZERO_EXTRA adds names, comma separated); the land
// regions are exempt on a planet without land and the radiation.held ones where held radiation is off (no CPU array, no GPU region)
export const MAY_BE_ZERO = new Set(['UO', 'OFLUX', 'FSTRESS', 'TOFD', 'OFLT', 'OSTRESS', 'OLAUNCH', 'OWAVE', 'OBETA', 'OBLOCK', 'ODIR', 'GWE', 'GWN', 'RUNOFF', 'DISCHARGE', 'ICELOSS', 'MLMK', 'MLMM', 'SNOW', 'SNOWALB', 'FOG', 'FOGMEAN', ...(process.env.MAY_BE_ZERO_EXTRA ?? '').split(',').filter(Boolean)]);
// mapped regions set by a threshold decision in each column (mixing regime, cumulus trigger, rain, fog, launch): f32 and
// f64 part at a few columns near the threshold, so these are judged by the share of cells that differ by more than TOL of
// the row's largest CPU value and fail past their bound here, max(3 × the worst share measured on
// runs/sixteen/sixteen128_day0730.bin through 17 steps, 1e-3), the rows from COND on measured on seventeen128 day 1095 and
// combined128 day 1460, and the rms of the cells under that threshold stays within 3 TOL; env SHARE sets one bound for every row, env
// DECISION_EXTRA adds rows (comma separated) under it
export const DECISION_BOUNDS = { MIX: 0.016, DEPTH: 0.048, MIXTOP: 0.048, ENTRAIN: 0.21, REGIME: 0.03, VRAD: 0.15, CTCOOL: 0.14, RAIN: 0.001, STEPRAIN: 0.026, CONV: 0.02, FOG: 0.0033, CUMF: 0.22, CUTOP: 0.072, CUCOVER: 0.051, CUWATER: 0.036, GWE: 0.001, GWN: 0.001, OWAVE: 0.001, COND: 0.001, MIXH: 0.0078, STRESS: 0.0165, DRAG: 0.0066, HEATX: 0.0048, DECKF: 0.001, MLMCOVER: 0.001, MLMENT: 0.001, OBETA: 0.001, OBLOCK: 0.001, OLAUNCH: 0.00156, OSTRESS: 0.00126, CONC: 0.001, STRAT: 0.001, USTAR: 0.001 };
// f32 rounding of the soil store: eight spacings at 300 kg/m²
const ABSOLUTE_FLOOR = { RUNOFF: 2.5e-4, DISCHARGE: 2.5e-4 };
const CARRIED = new Set(['PEN']);
const f32Grain = (x) => 2 ** (x ? Math.max(Math.floor(Math.log2(Math.abs(x))) - 23, -126) : -24);
const carriedLight = (ocean, rhoCp) => {
  const L = ocean.h.length / C, bound = new Float64Array(C);
  for (let i = 0; i < C; i++) for (let k = 0, top = 0; k < L && top < PENETRATION_DEPTH; k++) {
    const h = ocean.h[k * C + i];
    if (h > 0) bound[i] = Math.max(bound[i], 2 * rhoCp * f32Grain(h * (ocean.T[k * C + i] - HEAT_REFERENCE)));
    top += h;
  }
  return bound;
};
const pastCarry = (g, c, bound) => {
  const n = Math.min(g.length, c.length, bound.length);
  let ee = 0, cc = 0, past = 0;
  for (let x = 0; x < n; x++) { const e = Math.abs(g[x] - c[x]) - bound[x]; if (e > 0) { ee += e * e; past++; } cc += c[x] * c[x]; }
  return { rms: cc > 0 ? Math.sqrt(ee / cc) : ee > 0 ? Infinity : 0, past };
};
let oceanStepped = false;
// judged in kelvin, by the rms gap and the share of cells parting past `over`: about ten times the parity of seventeen64 day 365 and combined128 day 1460 through 2 steps
export const ABSOLUTE_BOUNDS = { T2M: { rms: 0.01, over: 0.05, share: 5e-4 } };
const absoluteGap = (g, c, over) => {
  const n = Math.min(g.length, c.length);
  let dd = 0, past = 0;
  for (let x = 0; x < n; x++) { const d = g[x] - c[x]; dd += d * d; if (Math.abs(d) > over) past++; }
  return { rms: Math.sqrt(dd / n), past, share: past / n };
};
for (const name of (process.env.DECISION_EXTRA ?? '').split(',').filter(Boolean)) DECISION_BOUNDS[name] ??= 1e-3;
export const DECISION_ROWS = new Set(Object.keys(DECISION_BOUNDS));
const NO_CLEAR_SKY = JSON.parse(process.env.RADIATION ?? '{}').clearSkyPass === false, CLEAR_SKY_ROWS = new Set(['ABSCLRSUM', 'OLRCLRSUM', 'RADCLROLR', 'RADCLRSW', 'RADTC']);
const drags = (model) => !!model.boundaryLayer?.implicitDrag;
const ROW_NEEDS = { T2M: ['no roughness exchange', (model) => model.exchange?.mode === 'roughness'], STRESS: ['no implicit drag', drags], SDRAG: ['no implicit drag', drags], WARM: ['no warm layer', (model) => !!model.skin?.options.warmLayer], COOL: ['no cool skin', (model) => !!model.skin?.options.coolSkin], UO: ['no relative wind', (model) => !!model.seaFlow] };
const decisionBound = (name) => (process.env.SHARE ? Number(process.env.SHARE) : DECISION_BOUNDS[name]);

const typed = (v) => ArrayBuffer.isView(v) && !(v instanceof DataView);
export function cpuArray(cpu, path) {
  const tail = path.endsWith('[tail]'), parts = path.replace('[tail]', '').split('.');
  let v = cpu;
  for (const key of parts) { if (v == null) return null; try { v = v[key]; } catch { return null; } }
  return typed(v) ? { array: v, tail } : null;
}

if (!NO_GPU) {
  const require = createRequire(`${REPO}/package.json`);
  const { create, globals } = await import(require.resolve('webgpu'));
  Object.assign(globalThis, globals);
  Object.defineProperty(globalThis, 'navigator', { value: { gpu: create((process.env.FLAGS ?? '').split(' ').filter(Boolean)) }, configurable: true });
}
const { createModel } = await import(`${REPO}/js/model.module.js`);
const { createGpuModel } = NO_GPU ? {} : await import(`${REPO}/js/gpu/model.gpu.js`);
const { readBuffer } = NO_GPU ? {} : await import(`${REPO}/js/gpu/device.module.js`);
const { HEAT_REFERENCE, PENETRATION_DEPTH } = await import(`${REPO}/js/ocean/layered.module.js`);

async function build(file) {
  const fs = await import(`${REPO}/scripts/figures/figureState.mjs`);
  const { withCadence } = await import(`${REPO}/js/cadence.module.js`);
  const { savedLevels } = await import(`${REPO}/js/stateFile.module.js`);
  const { savedDeckField, DECK_FIELDS, savedMoistField, MOIST_FIELDS, savedRadiationField, RADIATION_FIELDS } = await import(`${REPO}/js/physics/regrid.module.js`);
  const { stateGrid } = await import(`${REPO}/js/stateGrid.module.js`);
  const { stateGeography } = await import(`${REPO}/js/geography.module.js`);
  const saved = await fs.readState(file), sg = stateGeography(saved), topography = sg.topography, dt = fs.dt(saved.N);
  const mesh = stateGrid(saved, topography).mesh;
  const subgrid = process.env.NO_SUBGRID ? { subgrid: false } : {};
  const options = { terrain: stateTerrain(saved), topography, geography: sg.geography, land: { ...sg.land, ...JSON.parse(process.env.LAND ?? '{}') }, levels: savedLevels(saved), ...(process.env.MOIST ? { moist: JSON.parse(process.env.MOIST) } : {}), ...(process.env.SURFACE ? { surface: JSON.parse(process.env.SURFACE) } : {}), ...withCadence({ ocean: JSON.parse(process.env.OCEAN ?? '{}'), radiation: { clearSkyPass: true, ...JSON.parse(process.env.RADIATION ?? '{}') } }, dt, saved.N) };
  const load = (model) => {
    const { state } = model, C = model.mesh.nCells, K = model.core.K;
    ['pi', 'theta', 'u', 'surfaceT', 'q', 'qc', 'ice'].forEach((name, a) => state[a].set(saved[name]));
    model.seaIce.load(state[6], saved.concentration ?? null);
    for (const field of Object.keys(DECK_FIELDS)) model.radiation[field].set(savedDeckField(saved, field, model));
    for (const field of Object.keys(MOIST_FIELDS)) model.moist[field].set(savedMoistField(saved, field, model));
    for (const field of Object.keys(RADIATION_FIELDS)) model.radiation[field].set(savedRadiationField(saved, field, model));
    if (saved.boundaryDepth) model.boundaryLayer.depth.set(saved.boundaryDepth);
    if (saved.mixingTop) model.boundaryLayer.mixingTop.set(saved.mixingTop);
    if (saved.boundaryRegime) model.boundaryLayer.regime.set(saved.boundaryRegime);
    if (saved.boundaryBuoyancy) model.boundaryLayer.buoyancyFlux.set(saved.boundaryBuoyancy);
    if (saved.windSpeed) model.surface.windSpeed.set(saved.windSpeed);
    if (model.seaFlow && saved.seaFlow?.length === model.seaFlow.length) model.seaFlow.set(saved.seaFlow);
    if (saved.evaporation) model.radiation.evaporation.set(saved.evaporation);
    if (saved.exchangeHeat && model.exchange && !model.exchange.fixed) model.exchange.heat.set(saved.exchangeHeat);
    if (saved.exchangeWind && model.exchange) model.exchange.wind.set(saved.exchangeWind);
    if (model.skin) model.skin.load(saved);
    for (const field of ['cumulusCover', 'cumulusWater']) {
      const v = saved[field], t = model.moist[field];
      if (!v) continue;
      if (v.length === t.length) t.set(v); else if (t.length === K * C) { t.fill(0); t.set(v, K * C - v.length); }
    }
    const sv = model.moist.subcloudVirtual;
    if (saved.subcloudVirtual && saved.subcloudVirtual.length === sv.length) sv.set(saved.subcloudVirtual); else sv.fill(0);
    if (!model.gpu && model.surfaceGeopotential) for (let i = 0; i < C; i++) { const zs = model.surfaceGeopotential[i] / model.core.diagnostics.g; if (saved.boundaryDepth) model.boundaryLayer.depth[i] += zs; if (saved.mixingTop) model.boundaryLayer.mixingTop[i] += zs; }
    model.time = saved.time;
    if (model.load) model.load();
    model.ocean.load(saved.ocean, state[3], state[6]);
    if (saved.land) model.land.load(saved.land, state[6], saved.time, state[3]);
    return model;
  };
  if (process.env.BOUNDARY_LAYER) options.boundaryLayer = JSON.parse(process.env.BOUNDARY_LAYER);
  const cpu = GPU_DUMP ? null : load(createModel(mesh, { ...options, ...subgrid }));
  const gpu = NO_GPU ? null : load(await createGpuModel(mesh, { ...options, ...subgrid }));
  return { cpu, gpu, dt, label: `${file.replace(/.*\//, '')} N=${saved.N}` };
}

const { cpu, gpu, dt, label } = await build(FILE);
const C = (cpu ?? gpu).mesh.nCells;
const lines = [`# ${label}: CPU against GPU, dt ${dt} s (${REPO.replace(/.*\//, '')})`];
let failures = 0;
const fail = (text) => { failures++; lines.push(`**${text}**`); };
if (cpu && cpu.orography?.fields?.raster && !process.env.NO_SUBGRID) fail('the model took the raster sub-grid orography, with no form drag, which is not the configuration a saved state runs');

function stats(g, c, base = null, n = Math.min(g.length, c.length), gOffset = 0, cOffset = 0, gBase = null) {
  let dd = 0, cc = 0, bb = 0, gg = 0, mc = 0, md = 0, gz = 0, cz = 0, nan = 0;
  for (let x = 0; x < n; x++) {
    const a = g[gOffset + x], b = c[cOffset + x];
    if (!Number.isFinite(a) || !Number.isFinite(b)) { nan++; continue; }
    const d = a - b;
    dd += d * d; cc += b * b; md = Math.max(md, Math.abs(d)); mc = Math.max(mc, Math.abs(b));
    if (a === 0) gz++; if (b === 0) cz++;
    if (base) { const s = b - base[cOffset + x]; bb += s * s; }
    if (gBase) { const s = a - gBase[gOffset + x]; gg += s * s; }
  }
  let over = 0, under = 0;
  for (let x = 0; x < n; x++) { const a = g[gOffset + x], b = c[cOffset + x]; if (!Number.isFinite(a) || !Number.isFinite(b)) continue; const d = a - b; if (Math.abs(d) > TOL * mc) over++; else under += d * d; }
  return { own: base && gBase ? [Math.sqrt(bb / n), Math.sqrt(gg / n)] : null, share: over / n, rest: cc > 0 ? Math.sqrt(under / cc) : 0, rms: cc > 0 ? Math.sqrt(dd / cc) : dd > 0 ? Infinity : 0, max: mc > 0 ? md / mc : md > 0 ? Infinity : 0, step: base ? (bb > 0 ? Math.sqrt(dd / bb) : dd > 0 ? Infinity : 0) : null, gz: gz / n, cz: cz / n, scale: mc, nan, n };
}
const fmt = (v) => (v === null ? '-' : v === Infinity ? 'inf' : v === 0 ? '0' : v.toExponential(1));
const row = (name, s) => `| ${name} | ${fmt(s.rms)} | ${fmt(s.max)} | ${s.step === null ? '-' : fmt(s.step)} | ${(100 * s.gz).toFixed(0)}/${(100 * s.cz).toFixed(0)} | ${s.own ? `${fmt(s.own[0])} / ${fmt(s.own[1])}` : '-'} | ${fmt(s.scale)}${s.nan ? ` (${s.nan} non-finite on one side)` : ''} |`;
const HEAD = '| field | rel rms | rel max | rms / change since last report | zero % gpu/cpu | rms change cpu / gpu | max abs cpu |\n|---|---|---|---|---|---|---|';

function nonFinite(region) {
  let count = 0, first = -1;
  for (let x = 0; x < region.length; x++) if (!Number.isFinite(region[x])) { if (first < 0) first = x; count++; }
  return count ? `NON-FINITE (${count} of ${region.length}, first at layer ${Math.floor(first / C)} cell ${first % C})` : null;
}

function mappedReport(phNames, phArray, step, carry = null) {
  lines.push(`\n### PH regions against the fixed name map\n\n| PH | CPU array | rel rms | rel max | zero % gpu/cpu | max abs cpu | verdict |\n|---|---|---|---|---|---|---|`);
  const decisions = [];
  const unmapped = phNames.filter((n) => !(n in PH_MAP) && !(n in PH_INPUT));
  if (unmapped.length) fail(`UNMAPPED PH regions (add them to PH_MAP or PH_INPUT): ${unmapped.join(', ')}`);
  for (const [name, region] of phNames.map((n) => [n, phArray(n)])) {
    const path = PH_MAP[name];
    if (NO_CLEAR_SKY && CLEAR_SKY_ROWS.has(name)) { lines.push(`| ${name} | ${path} | | | | | n/a (no clear-sky pass) |`); continue; }
    if (ROW_NEEDS[name] && !ROW_NEEDS[name][1](cpu)) { lines.push(`| ${name} | ${path} | | | | | n/a (${ROW_NEEDS[name][0]}) |`); continue; }
    if (!path) { const bad = region ? nonFinite(region) : null; lines.push(`| ${name} | (input: ${PH_INPUT[name] ?? 'unmapped'}) | | | | | ${bad ?? 'ok'} |`); if (bad) failures++; continue; }
    const found = cpuArray(cpu, path);
    if (!found && path.startsWith('boundaryLayer.form') && !cpu.boundaryLayer.formDrag) { lines.push(`| ${name} | ${path} | | | | | n/a (form drag is off on this grid) |`); continue; }
    if (!found && path.startsWith('radiation.held.') && !cpu.radiation?.held) { lines.push(`| ${name} | ${path} | | | | | n/a (held radiation is off at this N) |`); continue; }
    if (!found) { fail(`MISSING CPU array ${path} for PH ${name}`); lines.push(`| ${name} | ${path} | | | | | MISSING |`); continue; }
    const c = found.tail && region ? found.array.subarray(found.array.length - region.length) : found.array;
    const landless = path.startsWith('land.') && (!cpu.geography || !cpu.geography.land || !Array.from(cpu.geography.land).some(Boolean));
    let cNonzero = false; for (let x = 0; x < c.length; x++) if (c[x] !== 0) { cNonzero = true; break; }
    if (!region) { lines.push(`| ${name} | ${path} | | | | ${cNonzero ? '' : '0'} | ${cNonzero ? 'cpu nonzero' : 'CPU ALL ZERO'} |`); if (!cNonzero && !MAY_BE_ZERO.has(name) && !landless) failures++; continue; }
    const s = stats(region, c, null, Math.min(region.length, c.length));
    let gNonzero = false; for (let x = 0; x < region.length; x++) if (region[x] !== 0) { gNonzero = true; break; }
    let verdict = 'ok', note = '';
    if (s.nan) verdict = `NON-FINITE (${s.nan} of ${Math.min(region.length, c.length)})`; else if (!gNonzero && cNonzero) verdict = 'GPU ALL ZERO'; else if (gNonzero && !cNonzero) verdict = 'CPU ALL ZERO'; else if (!gNonzero && !cNonzero) verdict = MAY_BE_ZERO.has(name) ? 'zero on both (listed)' : 'ZERO ON BOTH (not exercised)';
    else if (ABSOLUTE_BOUNDS[name]) { const b = ABSOLUTE_BOUNDS[name], a = absoluteGap(region, c, b.over); note = ` (rms ${fmt(a.rms)}, ${a.past} cells past ${b.over})`; verdict = a.rms > b.rms ? `OFF (rms > ${b.rms})` : a.share > b.share ? `OFF (share past ${b.over} > ${b.share})` : 'ok'; }
    else if (carry && CARRIED.has(name)) { const p = pastCarry(region, c, carry); note = ` (rel rms ${fmt(p.rms)} past the GPU's carry, ${p.past} cells past it)`; verdict = p.rms > TOL ? `OFF (past the carry > ${TOL})` : 'ok'; }
    else if (DECISION_ROWS.has(name)) verdict = s.share > decisionBound(name) ? `OFF (share > ${decisionBound(name)})` : s.rest > 3 * TOL ? `OFF (rms under the threshold > ${3 * TOL})` : 'ok';
    else if (!(s.rms <= TOL) && !(s.max * s.scale <= ABSOLUTE_FLOOR[name])) verdict = `OFF (> ${TOL})`;
    if (verdict !== 'ok' && verdict !== 'zero on both (listed)' && !((MAY_BE_ZERO.has(name) || landless) && verdict.startsWith('ZERO'))) failures++;
    if (DECISION_ROWS.has(name) && gNonzero && cNonzero) decisions.push(`| ${name} | ${path} | ${fmt(s.share)} | ${decisionBound(name)} | ${fmt(s.rms)} | ${fmt(s.rest)} | ${fmt(s.max)} | ${fmt(s.scale)} | ${verdict} |`);
    else lines.push(`| ${name} | ${path} | ${fmt(s.rms)} | ${fmt(s.max)} | ${(100 * s.gz).toFixed(0)}/${(100 * s.cz).toFixed(0)} | ${fmt(s.scale)} | ${verdict}${note} |`);
  }
  lines.push(`\n### decision divergence: share of cells\n\n| PH | CPU array | share > TOL | bound | rel rms | rel rms under TOL | rel max | max abs cpu | verdict |\n|---|---|---|---|---|---|---|---|---|`, ...decisions);
}

const held = (v) => (typed(v) ? v.slice() : Array.isArray(v) && v.every((x) => typeof x === 'number') ? Float64Array.from(v) : v);
const heldAll = (o) => Object.fromEntries(Object.entries(o).map(([k, v]) => [k, held(v)]));
async function gpuSnapshot() {
  await gpu.sync();
  const layout = gpu.gpu.layout.PH, ph = await readBuffer(gpu.gpu.device, gpu.gpu.buffers.PH, gpu.gpu.buffers.PH.size, Float32Array);
  const entries = Object.entries(layout).filter(([k, v]) => k !== 'total' && typeof v === 'number').sort((x, y) => x[1] - y[1]);
  const regions = Object.fromEntries(entries.map(([name, start], e) => [name, ph.slice(start, e + 1 < entries.length ? entries[e + 1][1] : layout.total ?? ph.length)]).filter(([, r]) => r.length >= C));
  return {
    state: gpu.state.map((a) => Float32Array.from(a)), concentration: Float32Array.from(gpu.seaIce.concentration), regions,
    ocean: gpu.ocean ? heldAll(await gpu.ocean.serialize()) : null, land: gpu.land ? heldAll(await gpu.land.serialize()) : null,
  };
}
const gpuAt = async (n) => (AGAINST ? deserialize(readFileSync(`${AGAINST}.step${n}`)) : gpuSnapshot());
let before = cpu ? cpu.state.map((a) => Float64Array.from(a)) : null;
let oceanBefore = cpu && cpu.ocean ? cpu.ocean.serialize() : null;
let gBefore = CPU_ONLY ? null : await gpuAt(0);
if (GPU_DUMP) writeFileSync(`${GPU_DUMP}.step0`, serialize(gBefore));
let done = 0;
for (const target of STEPS) {
  while (done < target) { const t0 = Date.now(); if (cpu) cpu.step(dt); const t1 = Date.now(); if (gpu) { await gpu.step(dt); await gpu.settle(); } done++; console.error(`step ${done}${cpu ? `: cpu ${((t1 - t0) / 1000).toFixed(1)} s` : ''}`); }
  if (GPU_DUMP) { writeFileSync(`${GPU_DUMP}.step${done}`, serialize(await gpuSnapshot())); continue; }
  lines.push(`\n## after ${done} step${done > 1 ? 's' : ''}`);
  if (CPU_ONLY) { mappedReport(Object.keys(PH_MAP), () => null, done); continue; }
  const g = await gpuAt(done);
  lines.push(`\n### state\n\n${HEAD}`);
  ['pi', 'theta', 'u', 'Ts', 'q', 'qc', 'ice'].forEach((name, a) => lines.push(row(name, stats(g.state[a], cpu.state[a], before[a], undefined, 0, 0, gBefore.state[a]))));
  lines.push(row('seaIce.concentration', stats(g.concentration, cpu.seaIce.concentration)));
  if (cpu.ocean && g.ocean) {
    const co = cpu.ocean.serialize(), go = g.ocean, gb = gBefore.ocean;
    oceanStepped ||= co.T.some((v, n) => v !== oceanBefore.T[n]);
    lines.push(`\n### ocean\n\n${HEAD}`);
    for (const key of ['h', 'u', 'T', 'S', 'eta']) lines.push(row(key, stats(go[key], co[key], oceanBefore[key], undefined, 0, 0, gb[key])));
    const content = (o, key) => Float64Array.from(o[key], (v, n) => v * o.h[n]);
    for (const key of ['T', 'S']) lines.push(row(`h·${key}`, stats(content(go, key), content(co, key), content(oceanBefore, key), undefined, 0, 0, content(gb, key))));
  }
  if (cpu.land && g.land) {
    const cl = cpu.land.serialize(), gl = g.land;
    lines.push(`\n### land\n\n${HEAD}`);
    for (const key of Object.keys(cl)) if (gl[key] && cl[key] && cl[key].length === gl[key].length && cl[key].length) lines.push(row(key, stats(gl[key], cl[key])));
  }
  mappedReport(Object.keys(g.regions), (n) => {
    let r = g.regions[n];
    if ((n === 'DEPTH' || n === 'MIXTOP') && cpu.surfaceGeopotential) { r = Float32Array.from(r); for (let i = 0; i < C; i++) r[i] += cpu.surfaceGeopotential[i] / cpu.core.diagnostics.g; }
    return r;
  }, done, oceanStepped && g.ocean ? carriedLight(g.ocean, cpu.ocean.rhoCp) : null);
  before = cpu.state.map((x) => Float64Array.from(x));
  if (cpu.ocean) oceanBefore = cpu.ocean.serialize();
  gBefore = g;
}
const text = lines.join('\n');
if (OUT) writeFileSync(OUT, text + '\n');
console.log(text);
console.log(failures ? `\n${failures} failure${failures > 1 ? 's' : ''}` : '\nno failures');
process.exit(failures ? 1 : 0);
